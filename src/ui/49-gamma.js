/* ====== gamma-ray spectrum of a time slice (needs ENDF-6 decay data: discrete gamma lines of MF 8 / MT 457) ======
   Decays of every nuclide of the chain in the slice [t1, t2] come from the exact balance  D_i = N_i(t1) − N_i(t2) + R_i·(overlap with the irradiation) + Σ_p b_p→i D_p
   (parents first). Each decay emits the lines of the nuclide with their absolute intensities; 1 keV per channel; Gaussian resolution; detector efficiency. */
const GAM_DEF = {
  fwhm: 1.5,
  fmode: 'const',
  eff: 1,
  effMode: 'const',
  curve: '',
  curveDef: '',
  ylog: 'log',
  rate: false,
  xray: true,
  top: 8,
  emax: 0,
  tau: 1,
  tauD: 0,
  ftp: 1,
  puCurve: false,
  ftpMode: 'auto',
  geD: 60,
  geL: 60,
  geDist: 25,
  geGap: 5,
  geDead: 0.7,
  geWin: 1,
  geScale: 1
};
let GAM = null; // state of the open window
const GAM_OLD_CURVE = '100 0.05\n300 0.04\n662 0.025\n1332 0.013\n3000 0.006'; // former default text, saved in old setups: it is now replaced by the (empty) default = Ge geometry model
function gamGeCurveText(g) {
  return [50, 100, 150, 200, 300, 400, 500, 662, 800, 1000, 1332, 1500, 2000, 3000]
    .map(E => {
      const v = geEffC(E, g);
      return E + '  ' + +v.tot.toPrecision(4) + '  ' + +v.peak.toPrecision(4);
    })
    .join('\n');
} // energy (keV), total, peak from the Ge geometry model
function gamDefaultCurve(g) {
  return String(g.curveDef || '').trim() ? g.curveDef : gamGeCurveText(g);
} // the user's default, else the Ge geometry model
const GAM_BOOL = ['rate', 'puCurve', 'xray']; // settings stored as true/false
let GAM_XR = true; // X-ray lines (e.xr) are part of the γ lines of the spectra
const gamHas = e => !!e && !!(e.gam || e.gc || (GAM_XR && e.xr));
function gamSet() {
  const s = irrSet(),
    g = (s.gam = Object.assign({}, GAM_DEF, s.gam || {}));
  if (String(g.curve).trim() === GAM_OLD_CURVE) g.curve = '';
  GAM_XR = g.xray !== false;
  return g;
}

// ---- numerics ----
function gamErfc(x) {
  // complementary error function, relative error < 1.2e-7 (Numerical Recipes, erfc by Chebyshev fit)
  const z = Math.abs(x),
    t = 1 / (1 + 0.5 * z),
    r =
      t *
      Math.exp(
        -z * z -
          1.26551223 +
          t *
            (1.00002368 +
              t *
                (0.37409196 +
                  t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))
      );
  return x >= 0 ? r : 2 - r;
}
const gamPhi = x => 0.5 * gamErfc(-x / Math.SQRT2);
function gamFwhm(E, g) {
  return g.fmode === 'sqrt' ? g.fwhm * Math.sqrt(Math.max(E, 0) / 1000) : g.fwhm;
} // keV; 'sqrt': the value is the FWHM at 1 MeV, ∝ √E
function gamCurve(text) {
  // per line "energy_keV total peak" (3 values) or "energy_keV peak" (2 values) → sorted [E, peak, total|null]
  const pts = [];
  for (const ln of String(text || '').split(/\r?\n/)) {
    const a = ln
      .trim()
      .split(/[\s,;]+/)
      .map(Number);
    if (a.length >= 2 && a[0] > 0 && isFinite(a[0])) {
      if (a.length >= 3 && isFinite(a[2]) && a[2] >= 0 && a[1] >= 0 && isFinite(a[1])) pts.push([a[0], a[2], a[1]]);
      else if (a[1] >= 0 && isFinite(a[1])) pts.push([a[0], a[1], null]);
    }
  }
  return pts.sort((a, b) => a[0] - b[0]);
}
function gamInterp(p, k) {
  // log-log interpolation of column k of the points (constant outside the range); p must not be empty
  if (p.length === 1) return () => p[0][k];
  return E => {
    if (E <= p[0][0]) return p[0][k];
    if (E >= p[p.length - 1][0]) return p[p.length - 1][k];
    let lo = 0,
      hi = p.length - 1;
    while (hi - lo > 1) {
      const m = (lo + hi) >> 1;
      if (p[m][0] <= E) lo = m;
      else hi = m;
    }
    const x0 = p[lo][0],
      y0 = p[lo][k],
      x1 = p[hi][0],
      y1 = p[hi][k],
      f = (Math.log(E) - Math.log(x0)) / (Math.log(x1) - Math.log(x0));
    return y0 > 0 && y1 > 0 ? Math.exp(Math.log(y0) + f * (Math.log(y1) - Math.log(y0))) : y0 + ((E - x0) / (x1 - x0)) * (y1 - y0);
  };
}
function gamEffTotFn(g) {
  // total efficiency: Ge model; curve model: typed totals, else the Ge-model total (never below the peak efficiency)
  if (g.effMode === 'ge') return E => geEffC(E, g).tot;
  if (g.effMode !== 'curve') return gamEffFn(g);
  const p = gamCurve(g.curve);
  if (!p.length) return E => geEffC(E, g).tot;
  const pt = p.filter(q => q[2] !== null);
  if (pt.length) return gamInterp(pt, 2);
  const pk = gamInterp(p, 1);
  return E => Math.max(pk(E), geEffC(E, g).tot);
}
function gamEffFn(g) {
  // E (keV) → peak efficiency (0..1+): constant, or log-log interpolation of the curve (constant outside); empty curve = Ge geometry model
  if (g.effMode === 'ge') return E => geEffC(E, g).peak;
  if (g.effMode !== 'curve') {
    const v = +g.eff;
    return () => (v >= 0 ? v : 0);
  }
  const p = gamCurve(g.curve);
  if (!p.length) return E => geEffC(E, g).peak;
  return gamInterp(p, 1);
}
function gamAdd(arr, E, counts, fwhm) {
  // adds `counts` of a line at E (keV) to the 1 keV channels, Gaussian of the given FWHM
  const n = arr.length,
    sg = fwhm / 2.3548;
  if (!(sg > 0.02)) {
    const c = Math.floor(E);
    if (c >= 0 && c < n) arr[c] += counts;
    return;
  }
  const c0 = Math.max(0, Math.floor(E - 5 * sg)),
    c1 = Math.min(n - 1, Math.ceil(E + 5 * sg));
  let prev = gamPhi((c0 - E) / sg);
  for (let c = c0; c <= c1; c++) {
    const nx = gamPhi((c + 1 - E) / sg);
    arr[c] += counts * (nx - prev);
    prev = nx;
  }
}

// ---- decays in the slice of every nuclide of a solved layer ----
function gamNearest(sol, t) {
  // index of the output time closest to t
  let lo = 0,
    hi = sol.t.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (sol.t[m] <= t) lo = m;
    else hi = m;
  }
  return Math.abs(sol.t[lo] - t) <= Math.abs(sol.t[hi] - t) ? lo : hi;
}
function gamDecays(sol, prod, Tirr, t1, t2) {
  // → { D: Float64Array per node, t1e, t2e }
  const k1 = t1 > 0 ? gamNearest(sol, t1) : -1,
    k2 = gamNearest(sol, t2),
    t1e = k1 >= 0 ? sol.t[k1] : 0,
    t2e = sol.t[k2],
    n = sol.ids.length,
    D = new Float64Array(n);
  const ov = Math.max(0, Math.min(t2e, Tirr) - Math.max(t1e, 0));
  for (const i of sol.order) {
    let d = (k1 >= 0 ? sol.num[i][k1] : 0) - sol.num[i][k2] + (prod.get(sol.ids[i]) || 0) * ov;
    for (const [p, f] of sol.par[i]) d += f * D[p];
    D[i] = d > 0 ? d : 0;
  }
  return { D, t1e, t2e };
}

// ---- the spectrum ----
function gamAll(e) {
  // discrete lines + the continuous gamma spectrum cut into 2 keV pseudo-lines (E = bin centre, I = photons per decay in the bin)
  if (!e) return [];
  const xr = GAM_XR && e.xr ? e.xr : null;
  if (!e.gc) return xr ? (e.gam || []).concat(xr) : e.gam || [];
  if (!e._gcl) {
    const { E, y } = e.gc,
      w = 2,
      out = [];
    let k = 0;
    for (let x = Math.ceil(E[0] / w) * w - w / 2; x < E[E.length - 1]; x += w) {
      while (k < E.length - 2 && E[k + 1] < x) k++;
      const dE = E[k + 1] - E[k],
        f = dE > 0 ? Math.min(1, Math.max(0, (x - E[k]) / dE)) : 0,
        v = y[k] + f * (y[k + 1] - y[k]);
      if (x > 0 && v * w > 1e-12) out.push([x, v * w]);
    }
    e._gcl = out;
  }
  return (e.gam || []).concat(e._gcl, xr || []);
}
function gamBuild(sol, prod, Tirr, t1, t2, g) {
  const { D, t1e, t2e } = gamDecays(sol, prod, Tirr, t1, t2),
    live = Math.max(t2e - t1e, 1e-12),
    eff = gamEffFn(g),
    nodes = [];
  let Emax = 0;
  sol.ids.forEach((id, i) => {
    const e = sol.info[i];
    if (gamHas(e) && D[i] > 0) for (const l of gamAll(e)) if (l[0] > Emax) Emax = l[0];
  });
  const nch = Math.max(100, Math.min(12000, Math.ceil(g.emax > 0 ? +g.emax : Emax + (6 * gamFwhm(Emax, g)) / 2.3548 + 2))),
    tot = new Float64Array(nch),
    totE = new Float64Array(nch);
  let sumD = 0,
    emitted = 0,
    detected = 0;
  sol.ids.forEach((id, i) => {
    const e = sol.info[i];
    if (!(D[i] > 0)) return;
    sumD += D[i];
    if (!gamHas(e)) return;
    const arr = new Float64Array(nch);
    let em = 0,
      de = 0;
    for (const [E, I] of gamAll(e)) {
      const c = D[i] * I,
        ef = Math.max(0, eff(E));
      em += c;
      de += c * ef;
      gamAdd(arr, E, c * ef, gamFwhm(E, g));
      gamAdd(totE, E, c, gamFwhm(E, g));
    }
    for (let c = 0; c < nch; c++) tot[c] += arr[c];
    emitted += em;
    detected += de;
    const top = gamAll(e)
      .map(l => [l[0], l[1] * D[i]])
      .sort((a, b) => b[1] - a[1])
      .slice(0, 3);
    nodes.push({ id, i, D: D[i], gper: (e.gtot || 0) + (e.gcTot || 0), em, de, arr, top, T: e.T });
  });
  nodes.sort((a, b) => b.de - a.de);
  return { nodes, tot, totE, nch, t1e, t2e, live, sumD, emitted, detected };
}

// ---- pile-up (see the Instructions, section 11.7) ----
// n = rate of pulses in the detector, τ = pulse width (resolving time): another pulse within ±τ distorts the event.
//   μ = n·τ;  P(event piled up) = 1 − exp(−2μ)  (two-sided; one-sided 1 − exp(−μ));  pairs per second = n²·τ;
//   spectrum: S_pu(E) = exp(−2μ)·S(E) + (pairs with both pulses full-energy)·(p ⊛ p)(E), p = S/ΣS, sum of the two energies.
function gamDetPerDecay(e, eff) {
  let s = 0;
  if (gamHas(e)) for (const l of gamAll(e)) s += l[1] * Math.max(0, eff(l[0]));
  return s;
}
function gamPile(B, sol, g) {
  const eff = gamEffFn(g),
    ge = g.effMode === 'ge',
    et = gamEffTotFn(g),
    dpt = ge ? sol.ids.map((id, i) => gamDetPerDecay(sol.info[i], et)) : null;
  let f = Math.max(1, +g.ftp || 1);
  if (ge && g.ftpMode !== 'man') {
    let a = 0,
      b = 0;
    sol.ids.forEach((id, i) => {
      const r = sol.act[i][gamNearest(sol, (B.t1e + B.t2e) / 2)];
      a += r * dpt[i];
      b += r * gamDetPerDecay(sol.info[i], eff);
    });
    f = b > 0 ? Math.max(1, a / b) : 1;
  }
  const tau = Math.max(0, +g.tau || 0) * 1e-6,
    dpd = sol.ids.map((id, i) => gamDetPerDecay(sol.info[i], eff));
  const rateAt = t => {
    const k = gamNearest(sol, t);
    let r = 0;
    for (let i = 0; i < dpd.length; i++) if (dpd[i] > 0) r += sol.act[i][k] * dpd[i];
    return r * f;
  }; // pulses per second at the node nearest to t
  const nbar = (B.detected / B.live) * f,
    n0 = rateAt(B.t1e),
    n1 = rateAt(B.t2e),
    nmax = Math.max(n0, n1, nbar);
  const mu = nbar * tau,
    mu0 = n0 * tau,
    mu1 = n1 * tau,
    mumax = nmax * tau;
  const pairs = (tau * B.live * (n0 * n0 + n0 * n1 + n1 * n1)) / 3; // ∫ n²τ dt with n linear between the limits
  const P = {
    f,
    tau,
    nbar,
    n0,
    n1,
    mu,
    mu0,
    mu1,
    mumax,
    Pm: 1 - Math.exp(-2 * mu),
    P0: 1 - Math.exp(-2 * mu0),
    Pmax: 1 - Math.exp(-2 * mumax),
    Pone: 1 - Math.exp(-mu),
    pairs,
    spec: null,
    sums: null,
    note: ''
  };
  if (g.puCurve && B.detected > 0 && tau > 0) {
    const nch = B.nch,
      p = new Float64Array(nch),
      idx = [];
    let mx = 0;
    for (let c = 0; c < nch; c++) {
      p[c] = B.tot[c] / B.detected;
      if (p[c] > mx) mx = p[c];
    }
    for (let c = 0; c < nch; c++) if (p[c] > 1e-6 * mx) idx.push(c);
    if (idx.length * idx.length > 4e7) P.note = 'the sum spectrum was skipped (too many channels with counts)';
    else {
      const nch2 = Math.min(2 * nch, 24000),
        sums = new Float64Array(nch2),
        w = pairs / (f * f);
      for (let a = 0; a < idx.length; a++) {
        const ia = idx[a],
          pa = p[ia] * w;
        for (let b = 0; b < idx.length; b++) {
          const c = ia + idx[b] + 1;
          if (c >= nch2) break;
          sums[c] += pa * p[idx[b]];
        }
      }
      const k0 = Math.exp(-2 * mu),
        spec = new Float64Array(nch2);
      for (let c = 0; c < nch2; c++) spec[c] = (c < nch ? k0 * B.tot[c] : 0) + sums[c];
      P.sums = sums;
      P.spec = spec;
    }
  }
  return P;
}
