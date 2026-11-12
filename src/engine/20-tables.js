// ===================== Tables and Monte Carlo transport =====================
function elasticSet(sp, mat, T, Rtot, Nh) {
  // per element: hard-scattering cross section, cut u_c, screening 2A, mass ratio; and the soft <theta^2> per cm
  const pc = Math.sqrt(T * (T + 2 * sp.M)), beta = pc / (T + sp.M), isE = sp.isE, Z1 = sp.Z1;
  const out = { sig: [], uc: [], A2: [], r: [], th2: 0, sigT: 0, dSh: 0, vs: 0 };
  let wsum = 0; const Cs = [];
  for (const e of mat.elems) {
    const Z2 = e.Z, z2 = isE ? Z2 * (Z2 + 1) : Z2 * Z2;
    let aScr, r = 0, pcm = pc;
    if (isE) aScr = 0.8853 * A0FM * Math.pow(Z2, -1 / 3);
    else { aScr = 0.8854 * A0FM / (Math.pow(Z1, 0.23) + Math.pow(Z2, 0.23)); r = sp.M / (e.A * AMU); pcm = pc / (1 + r); }
    const chi0 = HBARC / (pcm * aScr), aa = ALPHA * Z1 * Z2 / Math.max(beta, 1e-4);
    const A2 = 0.5 * chi0 * chi0 * (1.13 + 3.76 * aa * aa);
    const C = Math.pow(Z1 * E2 / (pcm * beta), 2) * z2;                  // fm2
    Cs.push({ C, A2, r }); wsum += e.n * z2;
  }
  mat.elems.forEach((e, k) => {
    const { C, A2, r } = Cs[k], z2 = sp.isE ? e.Z * (e.Z + 1) : e.Z * e.Z;
    const N = e.n, k0 = 2 * Math.PI * C * N * 1e-26;           // 1/cm per unit of 1/(u+A2)
    const Sigt = Nh > 0 && Rtot > 0 ? (N * z2 / wsum) * Nh / Rtot : 0;
    let uc = 2, sig = 0;
    if (Sigt > 0) {
      const X = Sigt / k0 + 1 / (2 + A2);
      uc = 1 / X - A2;
      if (uc < 1e-12) uc = 1e-12;
      if (uc >= 2) { uc = 2; sig = 0; } else sig = k0 * (1 / (uc + A2) - 1 / (2 + A2));
    }
    const th2 = 2 * k0 * (Math.log((uc + A2) / A2) + A2 / (uc + A2) - 1) / ((1 + r) * (1 + r));
    if (!isE && r > 0) {                                  // nuclear recoil energy: hard events (u > uc) deposit it explicitly, the soft part stays continuous
      const Tn = 4 * r / ((1 + r) * (1 + r)) * T, w1 = Math.min(uc, 2) + A2, w0 = A2, F = w => w - 2 * A2 * Math.log(w) - A2 * A2 / w;
      if (uc < 2) out.dSh += k0 * (Tn / 2) * (Math.log((2 + A2) / (uc + A2)) + A2 / (2 + A2) - A2 / (uc + A2));
      out.vs += k0 * Tn * Tn / 4 * (F(w1) - F(w0));
    }
    out.sig.push(sig); out.uc.push(uc); out.A2.push(A2); out.r.push(r); out.th2 += th2; out.sigT += sig;
  });
  return out;
}
function lsChi(sp, mat, T) {                       // Lindhard-Scharff straggling factor
  const v = betaOf(T, sp.M) / ALPHA, chi2 = v * v / Math.max(mat.Zbar, 1), chi = Math.sqrt(chi2);
  if (chi2 >= 3) return 1;
  return Math.max(0.5 * (1.36 * chi - 0.016 * chi * chi * chi), 0.03);
}
function buildTable(mat, sp, ph, Emin, Emax, ppd) {
  STOPSET = ph.stopSet || 'auto';
  ppd = ppd || 36;
  const n = Math.max(12, Math.ceil(Math.log10(Emax / Emin) * ppd) + 1), l0 = Math.log(Emin), dl = Math.log(Emax / Emin) / (n - 1);
  const nE = mat.elems.length, rho = mat.rho, Tc = ph.Tcut;
  const tb = { n, l0, dl, Emin, Emax, nE, E: new Float64Array(n), Se: new Float64Array(n), Sn: new Float64Array(n), Srad: new Float64Array(n), Sel: new Float64Array(n),
    Stot: new Float64Array(n), R: new Float64Array(n), vr: new Float64Array(n), sigD: new Float64Array(n), sigN: new Float64Array(n), th2: new Float64Array(n),
    sigEl: new Float64Array(n), fh: new Float64Array(n).fill(0), sigElE: [], ucE: [], A2E: [], rE: [] };
  for (let k = 0; k < nE; k++) { tb.sigElE.push(new Float64Array(n)); tb.ucE.push(new Float64Array(n)); tb.A2E.push(new Float64Array(n)); tb.rE.push(new Float64Array(n)); }
  // pass 1: stopping powers (per cm)
  for (let i = 0; i < n; i++) {
    const T = Emin * Math.exp(dl * i); tb.E[i] = T;
    let se, sn = 0, sr = 0, hard = { sig: 0, dS: 0 };
    if (sp.isE) {
      se = eCol(mat, T) * rho; sr = eRad(mat, T) * rho;
      if (ph.delta) hard = moellerHard(mat, T, Tc);
    } else {
      const s = stopHeavy(sp, mat, T); se = s.Se * rho; sn = ph.nucStop ? s.Sn * rho : 0;
      if (ph.delta) hard = heavyHard(mat, sp, T, Tc);
    }
    tb.Se[i] = se; tb.Sn[i] = sn; tb.Srad[i] = sr; tb.Stot[i] = se + sn + sr;
    tb.Sel[i] = Math.max(se - hard.dS, 0.02 * se);              // restricted (soft) collision stopping power
    tb.sigD[i] = hard.sig;
    // fluctuation variance of the soft part (MeV2/cm)
    const g = 1 + T / sp.M, b2 = 1 - 1 / (g * g);
    let v;
    if (sp.isE) { const Tce = Math.min(Tc, T / 2); v = KC / 2 * mat.ZA * rho / b2 * (ph.delta ? Tce : T / 2); }
    else {
      const z2 = zEff2(sp, mat, T), Tm = tmaxHeavy(T, sp.M), Tce = ph.delta ? Math.min(Tc, Tm) : Tm;
      v = KC / 2 * mat.ZA * rho * z2 / b2 * (Tce - b2 * Tce * Tce / (2 * Tm)) * lsChi(sp, mat, T);
    }
    tb.vr[i] = v;
    // nuclear reactions
    let sN = 0;
    if (!sp.isE && ph.reactions) {
      for (const e of mat.elems) {
        let sg = 0; for (const is of e.iso) sg += is.ab * sigmaIso(sp, e.Z, is.A, T);
        sN += e.n * sg * 1e-27 * ph.rxScale;
      }
    }
    tb.sigN[i] = sN;
  }
  // CSDA range
  const p0 = 0.5;
  tb.R[0] = tb.E[0] / ((1 + p0) * tb.Stot[0]);
  for (let i = 1; i < n; i++) {
    const f0 = tb.E[i - 1] / tb.Stot[i - 1], f1 = tb.E[i] / tb.Stot[i];
    tb.R[i] = tb.R[i - 1] + 0.5 * (f0 + f1) * dl;
  }
  // pass 2: elastic scattering (needs the range)
  for (let i = 0; i < n; i++) {
    const es = ph.scatter ? elasticSet(sp, mat, tb.E[i], tb.R[i], ph.Nh) : { sig: new Array(nE).fill(0), uc: new Array(nE).fill(2), A2: new Array(nE).fill(1), r: new Array(nE).fill(0), th2: 0, sigT: 0 };
    tb.th2[i] = es.th2; tb.sigEl[i] = es.sigT;
    if (!sp.isE && ph.nucStop && ph.scatter) {          // split nuclear stopping: hard recoils (explicit) + soft remainder (continuous, with its own small variance)
      const Sn = tb.Sn[i], fh = es.dSh > 0 ? Math.min(1, 0.95 * Sn / es.dSh) : 1; tb.fh[i] = fh; tb.Sn[i] = Sn - fh * es.dSh; tb.vr[i] += es.vs;
    }
    for (let k = 0; k < nE; k++) { tb.sigElE[k][i] = es.sig[k]; tb.ucE[k][i] = es.uc[k]; tb.A2E[k][i] = es.A2[k]; tb.rE[k][i] = es.r[k]; }
  }
  return tb;
}
