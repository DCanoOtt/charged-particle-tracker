/* ====== Neutron detector: efficiency and pulse-height Monte Carlo of a cylindrical organic scintillator (BC-501A / EJ-301, EJ-309, BC-400) ======
   Point source on the axis at the flight-path distance; free paths from the macroscopic cross sections of H and C; 1H(n,n) and 12C(n,n), the 12C inelastic, (n,p), (n,α)
   and break-up channels; light from the proton / α / carbon light-output curves (NRESP style: L_p = a1 E - a2 [1 - exp(-a3 E^a4)]); the light is smeared with the resolution of the scintillator (FWHM/L ∝ 1/√L); events above a threshold in MeVee are counted.
   Cross sections: ENDF/B-VIII.0 files for n+1H and n+12C (MF3 and the MF4 elastic angular distribution) when they are loaded, otherwise rough built-in tables (see ND_NOTE). */
const ND_NA = 6.02214076e23,
  ND_MN = 939.56542, // neutron mass (MeV)
  ND_A_C = 11.8969, // 12C / neutron mass ratio
  ND_A_H = 0.99862, // proton / neutron mass ratio
  ND_3A = 7.2748; // 12C -> 3 alpha threshold (excitation energy, MeV)
const ND_MAT = {
  BC501A: { name: 'BC-501A / EJ-301 / NE-213 (liquid)', rho: 0.874, hc: 1.212, curve: 'liq' },
  EJ309: { name: 'EJ-309 (liquid)', rho: 0.959, hc: 1.25, curve: 'liq' },
  STILBENE: { name: 'Stilbene (C14H12 crystal)', rho: 1.16, hc: 0.857, curve: 'liq' },
  BC400: { name: 'BC-400 (plastic)', rho: 1.032, hc: 1.103, curve: 'pla' }
};
// light-output parameters (MeVee, E in MeV): protons a1..a4 (Cecil et al. for NE-213; Madey / Verbinski for the plastic), alphas L = ka E^kb, carbon recoils L = ca E^cb
const ND_LP = {
  liq: { a1: 0.83, a2: 2.82, a3: 0.25, a4: 0.93, ka: 0.0382, kb: 1.6, ca: 0.0035, cb: 1.7 },
  pla: { a1: 0.95, a2: 8.0, a3: 0.1, a4: 0.9, ka: 0.0382, kb: 1.6, ca: 0.0035, cb: 1.7 }
};
const ND_NOTE =
  'Built-in tables: 1H(n,n) from effective-range theory (about 2 % up to 15 MeV), 12C values are smooth estimates without resonances (±20 %). Load the ENDF/B-VIII.0 files n-1H and n-12C for final numbers.';

function ndAtoms(rho, hc) {
  const nC = (rho * ND_NA) / (hc * 1.00794 + 12.0107);
  return { nC, nH: hc * nC };
}
const ndLp = (E, p) => Math.max(0, p.a1 * E - p.a2 * (1 - Math.exp(-p.a3 * Math.pow(E, p.a4))));
const ndLa = (E, p) => p.ka * Math.pow(Math.max(E, 0), p.kb);
const ndLc = (E, p) => (p.carbon === false ? 0 : p.ca * Math.pow(Math.max(E, 0), p.cb));
// proton energy that gives the light L (monotonic part of the curve)
function ndInvLp(L, p) {
  let lo = 0,
    hi = 200;
  for (let i = 0; i < 60; i++) {
    const m = 0.5 * (lo + hi);
    if (ndLp(m, p) < L) lo = m;
    else hi = m;
  }
  return 0.5 * (lo + hi);
}

// ---- built-in cross sections (barn, MeV) ----
// n-p scattering from the effective-range expansion: triplet a=5.424 fm, r=1.759 fm; singlet a=-23.74 fm, r=2.77 fm
function ndNpXs(E) {
  const k2 = E / 82.94,
    k = Math.sqrt(k2),
    cotT = (-1 / 5.424 + 0.5 * 1.759 * k2) / k,
    cotS = (-1 / -23.74 + 0.5 * 2.77 * k2) / k,
    sT = 1 / (1 + cotT * cotT),
    sS = 1 / (1 + cotS * cotS);
  return ((4 * Math.PI) / k2) * (0.75 * sT + 0.25 * sS) * 0.01;
}
// n-p total cross section above the end of the ACE table (20 MeV), barn: evaluated-data values (SAID / ENDF-type) at 20 ... 250 MeV, log-log interpolation
const ND_NPH = { E: [20, 30, 50, 80, 100, 150, 200, 250], y: [0.482, 0.3, 0.17, 0.095, 0.073, 0.05, 0.041, 0.036] };
function ndNpHigh(E) {
  const x = ND_NPH.E,
    y = ND_NPH.y;
  if (E <= x[0]) return y[0];
  if (E >= x[x.length - 1]) return y[y.length - 1];
  let i = 0;
  while (x[i + 1] < E) i++;
  const f = Math.log(E / x[i]) / Math.log(x[i + 1] / x[i]);
  return Math.exp(Math.log(y[i]) + f * Math.log(y[i + 1] / y[i]));
}
const ND_C = {
  E: [0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 60],
  el: [4.7, 4.7, 3.7, 2.6, 1.9, 1.75, 1.7, 1.6, 1.5, 1.3, 1.1, 0.95, 0.82, 0.78, 0.78, 0.77, 0.65, 0.55, 0.45, 0.35],
  ch: [
    {
      name: '(n,n′) 4.44 MeV',
      kind: 'inel',
      Q: -4.439,
      thr: 4.81,
      E: [4.81, 5, 6, 8, 10, 12, 14, 16, 20, 25, 30, 40, 60],
      y: [0, 0.03, 0.18, 0.28, 0.25, 0.23, 0.22, 0.18, 0.12, 0.08, 0.06, 0.04, 0.02]
    },
    {
      name: '(n,n′) 7.65 MeV → 3α',
      kind: 'brk',
      par: 'a',
      Q: -7.654,
      thr: 8.3,
      E: [8.3, 10, 12, 14, 16, 20, 25, 30, 40, 60],
      y: [0, 0.06, 0.1, 0.11, 0.1, 0.08, 0.06, 0.05, 0.03, 0.02]
    },
    {
      name: '(n,n′) 9.64 MeV → 3α',
      kind: 'brk',
      par: 'a',
      Q: -9.641,
      thr: 10.5,
      E: [10.5, 12, 14, 16, 20, 25, 30, 40, 60],
      y: [0, 0.05, 0.08, 0.08, 0.07, 0.05, 0.04, 0.03, 0.02]
    },
    {
      name: '(n,n′3α) continuum',
      kind: 'brk',
      par: 'a',
      Q: -11.5,
      thr: 12.5,
      E: [12.5, 14, 16, 18, 20, 25, 30, 40, 60],
      y: [0, 0.05, 0.12, 0.16, 0.18, 0.16, 0.12, 0.09, 0.06]
    },
    {
      name: '(n,α)',
      kind: 'chg',
      par: 'a',
      Q: -5.702,
      thr: 6.2,
      E: [6.2, 8, 10, 12, 14, 16, 20, 30, 40, 60],
      y: [0, 0.03, 0.05, 0.06, 0.065, 0.05, 0.03, 0.02, 0.012, 0.006]
    },
    { name: '(n,p)', kind: 'chg', par: 'p', Q: -12.587, thr: 13.7, E: [13.7, 15, 16, 20, 25, 30, 40, 60], y: [0, 0.01, 0.02, 0.03, 0.025, 0.02, 0.015, 0.01] }
  ]
};
const ndTable = (x, y) => ({ x: Float64Array.from(x), y: Float64Array.from(y), nbt: [x.length], int: [2] });
const ND_GRID = { l0: Math.log(0.005), dl: 0.004, n: Math.ceil(Math.log(250 / 0.005) / 0.004) + 1 };
const ndGridE = i => Math.exp(ND_GRID.l0 + i * ND_GRID.dl);
const ndPdfNp = (E, mu) => 1 + Math.max(0, 0.012 * (E - 10)) * mu * mu; // built-in: isotropic in the CM below 10 MeV, then forward-backward peaked
const ndPdfC = (E, mu) => (E < 1 ? 1 : Math.exp(((1.8 * E) / (E + 4)) * mu));

// angular sampler: CDF of the pdf f(E, mu) (CM cosine of the scattered neutron) on 120 energies x 65 cosines
function ndAngTable(fn) {
  const NE = 150,
    NM = 65,
    l0 = Math.log(0.005),
    dl = Math.log(250 / 0.005) / (NE - 1),
    cdf = new Float32Array(NE * NM);
  for (let i = 0; i < NE; i++) {
    const E = Math.exp(l0 + i * dl);
    let s = 0;
    const row = i * NM;
    cdf[row] = 0;
    let f0 = Math.max(fn(E, -1), 0);
    for (let j = 1; j < NM; j++) {
      const f1 = Math.max(fn(E, -1 + (2 * j) / (NM - 1)), 0);
      s += 0.5 * (f0 + f1);
      cdf[row + j] = s;
      f0 = f1;
    }
    if (!(s > 0)) for (let j = 0; j < NM; j++) cdf[row + j] = j / (NM - 1);
    else for (let j = 0; j < NM; j++) cdf[row + j] /= s;
  }
  return { cdf, NE, NM, l0, dl };
}
function ndSampleMu(t, E, r) {
  let i = Math.round((Math.log(Math.max(E, 0.005)) - t.l0) / t.dl);
  i = Math.min(Math.max(i, 0), t.NE - 1);
  const row = i * t.NM;
  let lo = 0,
    hi = t.NM - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (t.cdf[row + m] > r) hi = m;
    else lo = m;
  }
  const a = t.cdf[row + lo],
    b = t.cdf[row + hi],
    f = b > a ? (r - a) / (b - a) : 0.5;
  return -1 + (2 * (lo + f)) / (t.NM - 1);
}

// ---- reader of ENDF-6 incident-neutron files (MF3 cross sections, MF4 angular distribution of MT2) ----
function ndReadMF4(lines, a, b) {
  const cur = new EndfCursor(lines, a, b),
    h = cur.cont(),
    ltt = h[3] | 0,
    c2 = cur.cont(),
    E = [],
    L = [];
  const leg = () => {
      const ne = cur.skipTab2()[5] | 0;
      for (let e = 0; e < ne; e++) {
        const { c, v } = cur.list();
        E.push(c[1] * 1e-6);
        L.push({ lang: 0, a: v });
      }
    },
    tab = () => {
      const ne = cur.skipTab2()[5] | 0;
      for (let e = 0; e < ne; e++) {
        const t = cur.tab1();
        E.push(t.c[1] * 1e-6);
        L.push({ lang: 12, mu: t.x, f: t.y });
      }
    };
  if (ltt === 1) leg();
  else if (ltt === 2) tab();
  else if (ltt === 3) {
    leg();
    tab();
  }
  return { E: Float64Array.from(E), L, lct: c2[3] | 0, ltt };
}
const ND_ABS = new Set([16, 17, 102]),
  ND_CHG = { 103: 'p', 104: 'p', 105: 'p', 106: 'a', 107: 'a', 108: 'a', 111: 'p' };
function ndParseN(text) {
  const { lines, sec } = endfSections(text),
    s451 = sec.get('1/451');
  let ZA = 0,
    AWR = 0;
  if (s451) {
    const c = new EndfCursor(lines, s451[0], s451[1]).cont();
    ZA = Math.round(c[0]);
    AWR = c[1];
  } else {
    const c = endfFields(lines[1] || '');
    ZA = Math.round(c[0]);
    AWR = c[1];
  }
  const mts = {};
  let ang = null;
  for (const [key, [a, b]] of sec) {
    const [mf, mt] = key.split('/').map(Number);
    if (mf === 3 && (mt === 2 || mt === 4 || (mt >= 51 && mt <= 91) || mt === 22 || mt === 28 || ND_ABS.has(mt) || ND_CHG[mt])) {
      const cur = new EndfCursor(lines, a, b);
      cur.cont();
      const t = cur.tab1();
      if (!t.x.length) continue;
      mts[mt] = {
        mt,
        Q: (t.c[1] || t.c[0]) * 1e-6,
        QM: t.c[0] * 1e-6,
        xs: ndTable(
          Array.from(t.x, v => v * 1e-6),
          Array.from(t.y)
        )
      };
      mts[mt].xs.nbt = t.nbt;
      mts[mt].xs.int = t.int;
    } else if (mf === 4 && mt === 2) {
      try {
        ang = ndReadMF4(lines, a, b);
      } catch (e) {
        ang = null;
      }
    }
  }
  return { ZA, AWR, mts, ang, nMt: Object.keys(mts).length };
}

// ---- reader of ACE files (continuous-energy neutron, NJOY/MCNP format, ACE 1 and 2.0.x headers) ----
// returns the same structure as ndParseN, plus kT (MeV): MT2 from the ESZ block, further MTs from MTR/LSIG/SIG,
// the elastic angular distribution (CM cosine, LC<0 tabulated, LC>0 32 equiprobable bins) from the AND block
const ND_KT_ROOM = 2.5301e-8; // MeV, 293.6 K
function ndIsAce(text) {
  const L = text.slice(0, 600).split(/\r?\n/);
  return /^\s*\d+\.\d+\.\d+\s+\S+/.test(L[0] || '') || /^\s*\d+\.\d+[a-z]\s+[\d.]+\s+[\d.eE+-]+\s+\d+\/\d+\/\d+/.test(L[0] || '');
}
function ndParseAce(text) {
  const L = text.split(/\r?\n/),
    v2 = /^\s*\d+\.\d+\.\d+\s+\S+/.test(L[0] || '');
  let awr, kT;
  if (v2) {
    const t = L[1].trim().split(/\s+/);
    awr = +t[0];
    kT = +t[1];
  } else {
    const t = L[0].trim().split(/\s+/);
    awr = +t[1];
    kT = +t[2];
  }
  const h0 = v2 ? 8 : 6,
    nxs = L.slice(h0, h0 + 2)
      .join(' ')
      .trim()
      .split(/\s+/)
      .map(Number),
    jxs = L.slice(h0 + 2, h0 + 6)
      .join(' ')
      .trim()
      .split(/\s+/)
      .map(Number),
    xs = [];
  for (let i = h0 + 6; i < L.length && xs.length < nxs[0]; i++) {
    const t = L[i].trim();
    if (!t) continue;
    for (const w of t.split(/\s+/)) xs.push(parseFloat(w.replace(/([0-9.])([+-]\d)/, '$1e$2')));
  }
  if (xs.length < nxs[0] || !(nxs[1] > 0)) throw new Error('not a valid ACE neutron file');
  const X = i => xs[i - 1]; // 1-based XSS access
  const ZA = nxs[1],
    NES = nxs[2],
    NTR = nxs[3],
    NR = nxs[4],
    esz = jxs[0],
    mtr = jxs[2],
    lqr = jxs[3],
    lsig = jxs[5],
    sig = jxs[6],
    land = jxs[7],
    and = jxs[8];
  const Eg = Float64Array.from({ length: NES }, (_, i) => X(esz + i)),
    mts = {};
  const mk = (mt, Q, E, y) => {
    const t = ndTable(E, y);
    mts[mt] = { mt, Q, QM: Q, xs: t };
  };
  mk(
    2,
    0,
    Eg,
    Array.from({ length: NES }, (_, i) => X(esz + 3 * NES + i))
  );
  const keep = mt => mt === 4 || mt === 5 || mt === 22 || mt === 28 || (mt >= 51 && mt <= 91) || ND_ABS.has(mt) || ND_CHG[mt];
  for (let i = 0; i < NTR; i++) {
    const mt = X(mtr + i);
    if (!keep(mt)) continue;
    const p = sig + X(lsig + i) - 1,
      ie = X(p),
      ne = X(p + 1);
    if (!(ne > 0)) continue;
    mk(
      mt,
      X(lqr + i),
      Eg.subarray(ie - 1, ie - 1 + ne),
      Array.from({ length: ne }, (_, k) => X(p + 2 + k))
    );
  }
  // elastic angular distribution (LAND(1))
  let ang = null;
  const l1 = X(land);
  if (land > 0 && l1 > 0) {
    const p = and + l1 - 1,
      ne = X(p),
      E = [],
      Ls = [];
    let iso = true;
    for (let j = 0; j < ne; j++) {
      E.push(X(p + 1 + j));
      const lc = X(p + 1 + ne + j);
      if (lc === 0) {
        Ls.push({ lang: 12, mu: [-1, 1], f: [0.5, 0.5] });
        continue;
      }
      iso = false;
      const q = and + Math.abs(lc) - 1;
      if (lc < 0) {
        const np = X(q + 1); // LC < 0: tabulated (JJ, NP, cos, pdf, cdf); LC > 0: 32 equiprobable bins
        Ls.push({ lang: 12, mu: Array.from({ length: np }, (_, k) => X(q + 2 + k)), f: Array.from({ length: np }, (_, k) => X(q + 2 + np + k)) });
      } else {
        const b = Array.from({ length: 33 }, (_, k) => X(q + k)),
          mu = [b[0]],
          f = [];
        for (let k = 0; k < 32; k++) f.push(1 / (32 * Math.max(b[k + 1] - b[k], 1e-9)));
        // piecewise-constant pdf: sample it at the bin edges (average of the neighbouring bins)
        const fe = [f[0]];
        for (let k = 1; k < 32; k++) mu.push(b[k]), fe.push(0.5 * (f[k - 1] + f[k]));
        mu.push(b[32]);
        fe.push(f[31]);
        Ls.push({ lang: 12, mu, f: fe });
      }
    }
    ang = { E: Float64Array.from(E), L: Ls, lct: 2, ltt: 2, iso };
  }
  return { ZA, AWR: awr, kT, mts, ang, nMt: Object.keys(mts).length, ace: true };
}

// embedded evaluated data (51-nd-data.js) in the structure of ndParseN
const ndEmbC = {};
function ndEmb(k) {
  if (ndEmbC[k]) return ndEmbC[k];
  const e = ND_EMB[k],
    mts = {},
    NM = 33,
    mu = Array.from({ length: NM }, (_, j) => -1 + (2 * j) / (NM - 1));
  for (const m of Object.keys(e.mts)) {
    const a = e.mts[m],
      v = a.slice(2);
    mts[m] = { mt: +m, Q: a[0], QM: a[0], xs: ndTable(e.E.slice(a[1], a[1] + v.length), v) };
  }
  const L = e.aE.map((_, i) => ({ lang: 12, mu, f: e.aF.slice(i * NM, (i + 1) * NM) }));
  return (ndEmbC[k] = {
    ZA: e.ZA,
    AWR: e.AWR,
    kT: ND_KT_ROOM,
    mts,
    ang: { E: Float64Array.from(e.aE), L, lct: 2, ltt: 2 },
    nMt: Object.keys(mts).length,
    ace: true
  });
}
// reaction channels of 12C from the file: levels, break-up channels and charged-particle channels
function ndChannelsC(d) {
  const out = [],
    m = d.mts,
    lev = Object.keys(m)
      .map(Number)
      .filter(k => k >= 51 && k <= 91);
  const add = (k, name, kind, par, Q) => out.push({ name, kind, par, Q, xs: m[k].xs, mt: k });
  for (const k of lev) {
    const Q = m[k].Q;
    add(k, '(n,n′) MT ' + k + ' (' + (-Q).toFixed(2) + ' MeV)', -Q > ND_3A + 0.05 ? 'brk' : 'inel', 'a', Q);
  }
  if (!lev.length && m[4]) add(4, '(n,n′) MT 4', 'inel', 'a', -4.439);
  if (m[5]) add(5, '(n,x) MT 5 above 20 MeV (lumped nonelastic, Q ≈ −12 MeV)', 'brk', 'a', -12);
  if (m[22]) add(22, '(n,n′α) MT 22', 'brk', 'a', m[22].QM);
  if (m[28]) add(28, '(n,n′p) MT 28', 'brk', 'p', m[28].QM);
  for (const k of Object.keys(m).map(Number)) {
    if (ND_CHG[k]) add(k, '(n,x) MT ' + k, 'chg', ND_CHG[k], m[k].QM);
    else if (ND_ABS.has(k)) add(k, 'MT ' + k + ' (no light)', 'abs', '', m[k].QM);
  }
  return out;
}

// n+1H angular distribution: tabulated up to the last energy of the file, above it the last distribution times (1 + b mu^2), b = 0.025 (E - E_last) up to 3
// (forward and backward peaking of the n-p scattering in the CM at 50 - 200 MeV)
function ndAngHi(a, E, mu) {
  const El = a.E[a.E.length - 1];
  return E <= El ? d2Eval(a, E, mu) : d2Eval(a, El, mu) * (1 + Math.min(3, 0.025 * (E - El)) * mu * mu);
}
// ---- medium: macroscopic cross sections on a fine log grid ----
// data: { H: ndParseN() result | null, C: ndParseN() result | null }
function ndMedium(rho, hc, data) {
  const at = ndAtoms(rho, hc),
    nH = at.nH * 1e-24, // atoms per barn-cm: the cross sections are in barn
    nC = at.nC * 1e-24,
    G = ND_GRID,
    n = G.n,
    H = data && data.H,
    C = data && data.C,
    sH = new Float64Array(n),
    sEl = new Float64Array(n);
  const chans = C ? ndChannelsC(C) : ND_C.ch.map(c => ({ ...c, xs: ndTable(c.E, c.y) })),
    elTab = C && C.mts[2] ? C.mts[2].xs : ndTable(ND_C.E, ND_C.el),
    hTab = H && H.mts[2] ? H.mts[2].xs : null,
    cs = chans.map(() => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    const E = ndGridE(i);
    // the ACE tables end at 20 MeV: above, 1H follows an evaluated-data n-p cross section table (ND_NPH) normalised at the last point,
    // 12C channels keep their last value
    const xe = (t, e) => {
      const m = t.x.length;
      return e <= t.x[m - 1] ? tab1Eval(t, e) : t.y[m - 1];
    };
    sH[i] = hTab
      ? E <= hTab.x[hTab.x.length - 1]
        ? tab1Eval(hTab, E)
        : (hTab.y[hTab.x.length - 1] * ndNpHigh(E)) / ndNpHigh(hTab.x[hTab.x.length - 1])
      : ndNpXs(E);
    sEl[i] = xe(elTab, E);
    chans.forEach((c, k) => (cs[k][i] = Math.max(0, c.xs.x[c.xs.x.length - 1] < 19.9 ? tab1Eval(c.xs, E) : xe(c.xs, E))));
  }
  const tot = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = nH * sH[i] + nC * sEl[i];
    for (let k = 0; k < cs.length; k++) s += nC * cs[k][i];
    tot[i] = s;
  }
  const angH = ndAngTable(H && H.ang && H.ang.lct !== 1 ? (E, mu) => ndAngHi(H.ang, E, mu) : ndPdfNp),
    angC = ndAngTable(C && C.ang && C.ang.lct !== 1 ? (E, mu) => d2Eval(C.ang, E, mu) : ndPdfC);
  return { nH, nC, sH, sEl, cs, chans, tot, angH, angC, fromFile: { H: !!(H && hTab), C: !!(C && C.mts[2]), angH: !!(H && H.ang), angC: !!(C && C.ang) } };
}
function ndIdx(E) {
  const G = ND_GRID;
  let x = (Math.log(E) - G.l0) / G.dl;
  if (x < 0) x = 0;
  if (x > G.n - 1.001) x = G.n - 1.001;
  const i = Math.floor(x);
  return [i, x - i];
}
const ndInterp = (a, i, f) => a[i] + (a[i + 1] - a[i]) * f;

// mulberry32 random numbers
function ndRng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
// new neutron energy, lab cosine and the energy left in the target nucleus for a two-body collision with Q (target/neutron mass ratio a, CM cosine mu)
function ndKin(E, a, Q, mu) {
  const kap = 1 + (Q * (1 + a)) / (E * a);
  if (!(kap > 0)) return null;
  const sk = Math.sqrt(kap),
    den = 1 + kap * a * a + 2 * a * sk * mu,
    En = (E * den) / ((1 + a) * (1 + a));
  return { En, c: (1 + a * sk * mu) / Math.sqrt(Math.max(den, 1e-30)), Tres: E + Q - En };
}
function ndRot(ux, uy, uz, c, phi) {
  const s = Math.sqrt(Math.max(0, 1 - c * c)),
    sp = Math.sin(phi),
    cp = Math.cos(phi);
  if (Math.abs(uz) > 0.99999) return [s * cp, s * sp, uz > 0 ? c : -c];
  const d = Math.sqrt(1 - uz * uz);
  return [ux * c + (s * (ux * uz * cp - uy * sp)) / d, uy * c + (s * (uy * uz * cp + ux * sp)) / d, uz * c - s * cp * d];
}

// one history: a neutron of energy E enters at (x, y, 0) with direction (ux, uy, uz) into the cylinder (radius R, thickness T, cm); returns the total light (MeVee)
function ndHistory(med, lp, geo, E, x, y, ux, uy, uz, rng) {
  const R = geo.R,
    T = geo.T;
  let z = 0,
    L = 0,
    guard = 0;
  while (E > 0.01 && guard++ < 200) {
    const [i, f] = ndIdx(E),
      St = ndInterp(med.tot, i, f);
    if (!(St > 0)) break;
    const s = -Math.log(1 - rng()) / St;
    let dB = Infinity;
    if (uz > 0) dB = (T - z) / uz;
    else if (uz < 0) dB = -z / uz;
    const aa = ux * ux + uy * uy;
    if (aa > 1e-14) {
      const bb = x * ux + y * uy,
        cc = x * x + y * y - R * R,
        disc = bb * bb - aa * cc;
      if (disc >= 0) dB = Math.min(dB, (-bb + Math.sqrt(disc)) / aa);
    }
    if (s >= dB) break;
    x += ux * s;
    y += uy * s;
    z += uz * s;
    let u = rng() * St;
    const sH = med.nH * ndInterp(med.sH, i, f);
    let c, En;
    if (u < sH) {
      // 1H(n,n): the proton takes E - E'
      const k = ndKin(E, ND_A_H, 0, ndSampleMu(med.angH, E, rng()));
      L += ndLp(E - k.En, lp);
      En = k.En;
      c = k.c;
    } else {
      u -= sH;
      const el = med.nC * ndInterp(med.sEl, i, f);
      if (u < el) {
        const k = ndKin(E, ND_A_C, 0, ndSampleMu(med.angC, E, rng()));
        L += ndLc(E - k.En, lp);
        En = k.En;
        c = k.c;
      } else {
        u -= el;
        let ch = null;
        for (let k = 0; k < med.chans.length; k++) {
          const sk = med.nC * ndInterp(med.cs[k], i, f);
          if (u < sk) {
            ch = med.chans[k];
            break;
          }
          u -= sk;
        }
        if (!ch) ch = med.chans[med.chans.length - 1];
        if (!ch) break;
        if (ch.kind === 'inel' || ch.kind === 'brk') {
          const k = ndKin(E, ND_A_C, ch.Q, ndSampleMu(med.angIso, E, rng()));
          if (!k) break;
          const Tr = Math.max(k.Tres, 0);
          if (ch.kind === 'inel') L += ndLc(Tr, lp);
          else if (ch.par === 'p') L += ndLp(Tr * (11 / 12), lp);
          else {
            const Ech = Tr + Math.max(0, -ch.Q - ND_3A);
            L += 3 * ndLa(Ech / 3, lp);
          }
          En = k.En;
          c = k.c;
        } else {
          // charged-particle emission or absorption: the neutron disappears
          if (ch.kind === 'chg') {
            const Et = (E * ND_A_C) / (1 + ND_A_C) + ch.Q;
            if (Et > 0) L += ch.par === 'a' ? ndLa(Et * (9 / 13), lp) : ndLp(Et * (12 / 13), lp);
          }
          break;
        }
      }
    }
    const nd = ndRot(ux, uy, uz, c, 2 * Math.PI * rng());
    ux = nd[0];
    uy = nd[1];
    uz = nd[2];
    E = En;
  }
  return L;
}
// isotropic CM table shared by all inelastic channels
const ND_ISO = ndAngTable(() => 1);
// entry point on the front face and direction from the point source at distance Ld (cm); weight (Ld/r)^3 = solid-angle element relative to the axis
function ndEntry(geo, rng) {
  const r = geo.R * Math.sqrt(rng()),
    ph = 2 * Math.PI * rng(),
    x = r * Math.cos(ph),
    y = r * Math.sin(ph),
    d = Math.sqrt(x * x + y * y + geo.Ld * geo.Ld);
  return { x, y, ux: x / d, uy: y / d, uz: geo.Ld / d, w: Math.pow(geo.Ld / d, 3) };
}
// light-output resolution: FWHM/L = res / sqrt(L / 1 MeVee) (statistical, res = FWHM/L at 1 MeVee); Gaussian smearing of the total light of one history
function ndSmear(L, res, rng) {
  if (!(res > 0) || !(L > 0)) return L;
  const g = Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng()),
    sg = (res * Math.sqrt(L)) / 2.3548;
  return Math.max(0, L + g * sg);
}
/* ---- PPAC with a U-238 target: the neutron-induced fission is detected through its two back-to-back fragments ----
   Efficiency (events / neutrons on the target spot) = [1 - exp(-n sigma_f(E))] x eta, with n the U-238 atoms per cm2 of the layer and eta the probability that the fragments are detected: a Monte Carlo
   over the point of the fission in the spot (uniform disc of the beam-spot radius) and in depth (uniform in the areal thickness), isotropic back-to-back fragments, escape from the layer when the areal path
   (t - s)/|cos theta| (forward) or s/|cos theta| (backward) is shorter than the effective range, and impact on the active disc of the PPAC (distance d to the target plane, radius Ra). */
const ND_NU = (6.02214076e23 * 1e-6) / 238.0508; // U-238 atoms per cm2 and per microgram/cm2
function ppacSigma(E) {
  const x = ND_U238.E,
    y = ND_U238.sf,
    n = x.length;
  if (!(E > x[0])) return E > 0 ? y[0] * (E / x[0]) : 0;
  if (E >= x[n - 1]) return y[n - 1]; // constant above the end of the evaluation
  let lo = 0,
    hi = n - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (x[m] > E) hi = m;
    else lo = m;
  }
  return y[lo] + ((y[hi] - y[lo]) * (E - x[lo])) / (x[hi] - x[lo]);
}
const ppacProb = (S, E) => 1 - Math.exp(-S.n * 1e-24 * ppacSigma(E));
function ppacEta(P, N, seed) {
  const rng = ndRng(seed || 99),
    both = P.mode !== 'any',
    pe = P.eff;
  let hit = 0;
  for (let h = 0; h < N; h++) {
    const rr = P.rs * Math.sqrt(rng()),
      ph = 2 * Math.PI * rng(),
      x = rr * Math.cos(ph),
      y = rr * Math.sin(ph),
      s = rng() * P.tg,
      c = 2 * rng() - 1,
      sn = Math.sqrt(1 - c * c),
      az = 2 * Math.PI * rng(),
      ca = Math.cos(az),
      sa = Math.sin(az),
      lat = c !== 0 ? (P.dist * sn) / Math.abs(c) : Infinity,
      // fragment A flies along cos theta = c, fragment B opposite
      one = (cz, sg) => {
        const path = cz > 0 ? (P.tg - s) / cz : s / -cz;
        if (!(path < P.range)) return false;
        const px = x + sg * lat * ca,
          py = y + sg * lat * sa;
        return px * px + py * py <= P.rad * P.rad && rng() < pe;
      },
      a = one(c, 1),
      b = one(-c, -1);
    if (both ? a && b : a || b) hit++;
  }
  const eta = hit / N;
  return { eta, err: Math.sqrt((eta * (1 - eta)) / N) };
}
function ppacSetup(det) {
  const P = { rs: det.prs, tg: det.ptg, range: det.prg, dist: det.pdist, rad: det.prad, eff: Math.min(1, Math.max(0, det.peff / 100)), mode: det.pmode },
    r = ppacEta(P, Math.max(10000, det.pN | 0), 4711);
  return {
    ppac: true,
    P,
    n: det.ptg * ND_NU,
    eta: r.eta,
    etaErr: r.err,
    geo: { R: det.prs, T: 0, Ld: det.Lm * 100 },
    lp: Object.assign({}, ND_LP.liq),
    Lth: 0,
    res: 0
  };
}
function ndSetup(det, data) {
  if (det.kind === 'ppac') return ppacSetup(det);
  const med = ndMedium(det.rho, det.hc, data);
  med.angIso = ND_ISO;
  const lp = Object.assign({}, ND_LP[(ND_MAT[det.mat] || ND_MAT.BC501A).curve], det.lp || {}),
    geo = { R: det.D / 2, T: det.T, Ld: det.Lm * 100 };
  return { med, lp, geo, Lth: det.th / 1000, res: det.res > 0 ? det.res / 100 : 0 };
}
// efficiency at one energy: { eps, err, eps1 }, eps1 = single-collision estimate for normal incidence on the axis
// w > 0: the neutron energy is drawn log-uniformly in E exp(±w), i.e. the efficiency is the average over the energy bin of the curve point
// (resonances of 12C narrower than the grid would otherwise alias into jagged noise)
function ndEff(S, E, N, seed, w) {
  if (S.ppac) {
    // analytic: fission probability (averaged over the energy window of the curve point) x detection probability of the fragments
    const ns = w > 0 ? 41 : 1;
    let p = 0;
    for (let k = 0; k < ns; k++) p += ppacProb(S, ns > 1 ? E * Math.exp(((2 * (k + 0.5)) / ns - 1) * w) : E);
    p /= ns;
    return { eps: S.eta * p, err: S.etaErr * p, eps1: p };
  }
  const rng = ndRng(seed || 12345),
    E0 = E;
  let sw = 0,
    sw2 = 0,
    sd = 0;
  for (let h = 0; h < N; h++) {
    const e = ndEntry(S.geo, rng),
      Eh = w > 0 ? E0 * Math.exp((2 * rng() - 1) * w) : E0,
      L = ndSmear(ndHistory(S.med, S.lp, S.geo, Eh, e.x, e.y, e.ux, e.uy, e.uz, rng), S.res, rng);
    sw += e.w;
    sw2 += e.w * e.w;
    if (L > S.Lth) sd += e.w;
  }
  const eps = sw > 0 ? sd / sw : 0,
    neff = sw2 > 0 ? (sw * sw) / sw2 : N,
    Eth = ndInvLp(S.Lth, S.lp),
    ns = w > 0 ? 41 : 1;
  let e1 = 0;
  for (let k = 0; k < ns; k++) {
    const Ek = ns > 1 ? E0 * Math.exp(((2 * (k + 0.5)) / ns - 1) * w) : E0,
      [i, f] = ndIdx(Ek),
      St = ndInterp(S.med.tot, i, f),
      SH = S.med.nH * ndInterp(S.med.sH, i, f);
    e1 += St > 0 ? (1 - Math.exp(-St * S.geo.T)) * (SH / St) * Math.max(0, (Ek - Eth) / Ek) : 0;
  }
  e1 /= ns;
  return { eps, err: Math.sqrt(Math.max(eps * (1 - eps), 0) / Math.max(neff, 1)), eps1: e1 };
}
// time of flight (ns) over Lm metres of a neutron of kinetic energy E (MeV), and the energy of a neutron with that time of flight
function ndTof(E, Lm) {
  const g = 1 + E / ND_MN,
    b = Math.sqrt(1 - 1 / (g * g));
  return Lm / (b * 0.299792458);
}
function ndEfromT(t, Lm) {
  const b = Lm / (0.299792458 * t);
  return b >= 1 ? Infinity : ND_MN * (1 / Math.sqrt(1 - b * b) - 1);
}
// full-statistics Monte Carlo: energies drawn from the cumulative distribution cdf (over Ecent); returns the pulse-height histogram (nb bins up to Lmax) and, for the events above the
// threshold, the reconstructed energy (nE bins up to Eo) and time-of-flight (nT bins in [t0, t1]) histograms; the time of flight of the source energy over the flight path Lm
// is smeared with a Gaussian of sigma st (ns). Weights sum to sw for all histories.
// PPAC: every history carries the detection probability of its neutron energy as a weight (no event-by-event sampling of the rare fission), so that the reconstructed spectra are smooth
function ndPulsePpac(S, Ecent, cdf, N, nb, seed, rec) {
  const rng = ndRng(seed || 777),
    r = rec || { nE: 1, Eo: 1, nT: 1, t0: 0, t1: 1, Lm: S.geo.Ld / 100, st: 0 },
    recE = new Float64Array(r.nE),
    recT = new Float64Array(r.nT);
  let sd = 0,
    sd2 = 0,
    nd = 0;
  for (let h = 0; h < N; h++) {
    const u = rng();
    let lo = 0,
      hi = cdf.length - 1;
    while (hi - lo > 1) {
      const m = (lo + hi) >> 1;
      if (cdf[m] > u) hi = m;
      else lo = m;
    }
    const En = r.dE > 0 ? Ecent[lo] + (rng() - 0.5) * r.dE : Ecent[lo],
      w = S.eta * ppacProb(S, En);
    sd += w;
    sd2 += w * w;
    if (w > 0) nd++;
    if (rec) {
      const g = r.st > 0 ? Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng()) : 0,
        t = ndTof(En, r.Lm) + g * r.st,
        bt = Math.floor(((t - r.t0) / (r.t1 - r.t0)) * r.nT);
      if (bt >= 0 && bt < r.nT) recT[bt] += w;
      const Er = t > 0 ? ndEfromT(t, r.Lm) : Infinity,
        be = Math.floor((Er / r.Eo) * r.nE);
      if (be >= 0 && be < r.nE) recE[be] += w;
    }
  }
  return { all: new Float64Array(nb), recE, recT, sw: N, sw2: N, sd, sd2, nd, N, ppac: true, relEta: S.eta > 0 ? S.etaErr / S.eta : 0 };
}
function ndPulse(S, Ecent, cdf, N, nb, Lmax, seed, rec) {
  if (S.ppac) return ndPulsePpac(S, Ecent, cdf, N, nb, seed, rec);
  const rng = ndRng(seed || 777),
    all = new Float64Array(nb),
    r = rec || { nE: 1, Eo: 1, nT: 1, t0: 0, t1: 1, Lm: S.geo.Ld / 100, st: 0 },
    recE = new Float64Array(r.nE),
    recT = new Float64Array(r.nT);
  let sw = 0,
    sw2 = 0,
    sd = 0,
    sd2 = 0,
    nd = 0;
  for (let h = 0; h < N; h++) {
    const u = rng();
    let lo = 0,
      hi = cdf.length - 1;
    while (hi - lo > 1) {
      const m = (lo + hi) >> 1;
      if (cdf[m] > u) hi = m;
      else lo = m;
    }
    const e = ndEntry(S.geo, rng),
      En = r.dE > 0 ? Ecent[lo] + (rng() - 0.5) * r.dE : Ecent[lo], // energy uniform inside the bin of the histogram
      L = ndSmear(ndHistory(S.med, S.lp, S.geo, En, e.x, e.y, e.ux, e.uy, e.uz, rng), S.res, rng);
    sw += e.w;
    sw2 += e.w * e.w;
    if (L > S.Lth) {
      sd += e.w;
      sd2 += e.w * e.w;
      nd++;
    }
    if (L > 0) {
      const b = Math.floor((L / Lmax) * nb);
      if (b >= 0 && b < nb) all[b] += e.w;
    }
    if (L > S.Lth && rec) {
      const g = r.st > 0 ? Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng()) : 0,
        t = ndTof(En, r.Lm) + g * r.st,
        bt = Math.floor(((t - r.t0) / (r.t1 - r.t0)) * r.nT);
      if (bt >= 0 && bt < r.nT) recT[bt] += e.w;
      const Er = t > 0 ? ndEfromT(t, r.Lm) : Infinity,
        be = Math.floor((Er / r.Eo) * r.nE);
      if (be >= 0 && be < r.nE) recE[be] += e.w;
    }
  }
  return { all, recE, recT, sw, sw2, sd, sd2, nd, N };
}

// pile-up of the detected-neutron rate n (s⁻¹) for the pulse width tau (ns):
//   'par' (paralysable, window ±τ): recorded n·exp(−2nτ), loss 1 − exp(−2nτ)
//   'non' (non-paralysable, dead time τ): recorded n/(1 + nτ), loss nτ/(1 + nτ)
//   'off': no loss.   Returns {mu: nτ, loss, rec}
function ndPile(n, tauNs, mode) {
  const mu = Math.max(n, 0) * Math.max(+tauNs || 0, 0) * 1e-9;
  if (mode === 'off') return { mu, loss: 0, rec: Math.max(n, 0) };
  const loss = mode === 'non' ? mu / (1 + mu) : 1 - Math.exp(-2 * mu);
  return { mu, loss, rec: Math.max(n, 0) * (1 - loss) };
}
