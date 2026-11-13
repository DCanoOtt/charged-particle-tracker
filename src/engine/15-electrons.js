// ---------- electrons ----------
function eCol(mat, T) {                          // collision stopping power: ICRU 90 table (air, water, graphite at the reference density) or the ICRU 37 formula below
  const f = eColF(mat, T);
  if (STOPSET === 'auto') {
    const i9 = id90(mat);
    if (i9 && Math.abs(mat.rho / RHO90[i9] - 1) < 0.15) { const c = tabCurve('E90', 0, i9); if (c) { if (T >= c.Emin && T <= c.Emax) return tabEval(c, T); const ed = T < c.Emin ? c.Emin : c.Emax; return f * tabEval(c, ed) / eColF(mat, ed); } }
  }
  return f;
}
function eColF(mat, T) {                         // ICRU 37 / Rohrlich-Carlson, Joy-Luo low-energy form
  const tau = T / ME, g = tau + 1, b2 = 1 - 1 / (g * g), bg = Math.sqrt(b2) * g;
  const J = mat.I * 1e-6 / (1 + 0.85 * mat.I * 1e-6 / T);
  const F = 1 - b2 + (tau * tau / 8 - (2 * tau + 1) * Math.LN2) / (g * g);
  let term = Math.log(tau * tau * (tau + 2) / (2 * (J / ME) * (J / ME))) + F - densityEffect(mat, bg);
  if (term < 0.05) term = 0.05;
  return KC / 2 * mat.ZA / b2 * term;
}
function eRad(mat, T) {                          // radiative stopping power, MeV cm2/g
  const ALPHA_re2 = ALPHA * (2.8179403262e-13) ** 2;
  const f = T / (T + 12);
  const phi = 16 / 3 + (mat.phiInf - 16 / 3) * f;
  return NA * ALPHA_re2 * mat.radK * (T + ME) * phi;
}
function moellerHard(mat, T, Tc) {               // {sig (1/cm), dS (MeV/cm)} for delta production T_delta > Tc, electrons
  const E = T; if (Tc >= E / 2) return { sig: 0, dS: 0 };
  const g = 1 + T / ME, b2 = 1 - 1 / (g * g), e0 = Tc / E;
  const gg = (g - 1) / g, h = (2 * g - 1) / (g * g);
  const sig = 1 / e0 - 1 / (1 - e0) + gg * gg * (0.5 - e0) - h * Math.log((1 - e0) / e0);
  const dS = Math.log(0.5 / e0) + (2 + Math.log(0.5)) - (1 / (1 - e0) + Math.log(1 - e0)) + gg * gg * (0.125 - e0 * e0 / 2) - h * (Math.LN2 + Math.log(1 - e0));
  const c = KC / 2 * mat.ZA / b2 * mat.rho;
  return { sig: c / E * sig, dS: c * dS };
}
function heavyHard(mat, sp, T, Tc) {               // delta rays from heavy particles
  const Tm = tmaxHeavy(T, sp.M); if (Tc >= Tm) return { sig: 0, dS: 0, Tm };
  const g = 1 + T / sp.M, b2 = 1 - 1 / (g * g), z2 = zEff2(sp, mat, T);
  const c = KC / 2 * mat.ZA * z2 / b2 * mat.rho;
  return { sig: c * (1 / Tc - 1 / Tm - b2 / Tm * Math.log(Tm / Tc)), dS: 2 * c * (0.5 * Math.log(Tm / Tc) - 0.5 * b2 * (1 - Tc / Tm)), Tm };
}
