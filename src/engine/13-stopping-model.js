// ---------- heavy charged particles: electronic stopping ----------
const SPT = {            // tuning constants of the semi-empirical low/intermediate-energy model
  etaMin: 0.13, shellPow: 1.0, shellLow: 1.0, regA: 1.0, lsFac: 1.0, blendN: 1.0, Lmin: 0.02, bk: 0, E1: 0.8, E2: 3.0, phK: 0.16, phE0: 100, phW: 1.5, gI:0.3, gLS:0.8, gN:1
};
function shellC(I, eta) {                          // Barkas/Bichsel shell correction (total C of the atom), I in eV
  const e2 = 1 / (eta * eta), e4 = e2 * e2, e6 = e4 * e2;
  return (0.422377 * e2 + 0.0304043 * e4 - 0.00038106 * e6) * 1e-6 * I * I + (3.858019 * e2 - 0.1667989 * e4 + 0.00157955 * e6) * 1e-9 * I * I * I;
}
function shellTerm(mat, eta) {                       // sum_i f_i C_i / Z_i
  let s = 0; const em = Math.max(eta, SPT.etaMin);
  const sc = eta >= SPT.etaMin ? 1 : SPT.shellLow * Math.pow(eta / SPT.etaMin, SPT.shellPow);
  for (const e of mat.elems) s += e.fel * shellC(e.I, em) / e.Z;
  return s * sc;
}

function blochL2(y) {                              // L2 = Re psi(1 + i y) - psi(1)
  if (y < 0.3) { const y2 = y * y; return -y2 * (1.2020569 - y2 * (1.0369278 - 1.0083493 * y2 + 1.0020084 * y2 * y2)); }
  let re = 1, im = y, acc = 0;                       // z = 1 + iy; shift upward 8 steps
  for (let k = 0; k < 8; k++) { const d = re * re + im * im; acc += re / d; re += 1; }
  const d = re * re + im * im, lnw = 0.5 * Math.log(d);
  // 1/w, 1/w^2, 1/w^4, 1/w^6 real parts
  const ir = re / d, ii = -im / d;
  const w2r = ir * ir - ii * ii, w2i = 2 * ir * ii;
  const w4r = w2r * w2r - w2i * w2i, w4i = 2 * w2r * w2i;
  const w6r = w4r * w2r - w4i * w2i;
  const psi = lnw - 0.5 * ir - w2r / 12 + w4r / 120 - w6r / 252 - acc;
  return -(psi + 0.5772156649);
}
function effCharge2(sp, T) {
  const bet = betaOf(T, sp.M), vr = bet / ALPHA;            // v/v0
  if (sp.Z1 === 1) return 1;
  if (sp.Z1 === 2) {
    const B = Math.log(Math.max(T / sp.A * 1000, 1e-3));
    const arg = 0.7446 + 0.1429 * B + 0.01562 * B * B - 0.00267 * B ** 3 + 1.325e-6 * B ** 8;
    const g2 = 1 - Math.exp(-Math.min(Math.max(arg, 0), 30));
    return 4 * g2;
  }
  const Z1 = sp.Z1, y = vr / Math.pow(Z1, 2 / 3);
  let q = 1 - Math.exp(0.803 * Math.pow(y, 0.3) - 1.3167 * Math.pow(y, 0.6) - 0.38157 * y - 0.008983 * y * y);
  if (!(q > 0)) q = 0; if (q > 1) q = 1;
  const vF = 1.3, Lam = 2 * Math.pow(1 - q, 2 / 3) / (Math.pow(Z1, 1 / 3) * (1 - (1 - q) / 7));
  const g = q + (1 - q) / (2 * vF * vF) * Math.log(1 + Math.pow(4 * Lam * vF / 1.919, 2));
  const z = Z1 * Math.min(g, 1);
  return z * z;
}
function betaOf(T, M) { const g = 1 + T / M; return Math.sqrt(Math.max(1 - 1 / (g * g), 0)); }
function tmaxHeavy(T, M) {
  const g = 1 + T / M, b2 = 1 - 1 / (g * g);
  return 2 * ME * b2 * g * g / (1 + 2 * g * ME / M + (ME / M) * (ME / M));
}
// Lindhard-Scharff-type low-velocity electronic stopping, per g/cm2 (MeV cm2/g)
function stopHeavyHigh(sp, mat, T) {               // Bethe-Bloch with corrections, MeV cm2/g
  const g = 1 + T / sp.M, b2 = 1 - 1 / (g * g), bg2 = b2 * g * g, bg = Math.sqrt(bg2);
  const Tm = tmaxHeavy(T, sp.M), I = mat.I * 1e-6;
  const z2 = effCharge2(sp, T), zr = Math.sqrt(z2);
  const X = 2 * ME * bg2 * Tm / (I * I);
  let L = 0.5 * Math.log(1 + SPT.regA * X) - b2 - densityEffect(mat, bg) / 2 - shellTerm(mat, bg);
  const y = zr * ALPHA / Math.max(Math.sqrt(b2), 1e-6);            // Bloch correction (exact: Re psi(1+iy) + gamma)
  L += blochL2(y) + SPT.bk * y * Math.exp(-y * y);
  if (L < SPT.Lmin) L = SPT.Lmin;
  return KC * mat.ZA * z2 / b2 * L;
}
function stopNuclear(sp, mat, T) {               // ZBL universal nuclear stopping, MeV cm2/g
  const Z1 = sp.Z1, M1 = sp.M / AMU, Ek = T * 1000;     // keV
  let s = 0;
  for (const e of mat.elems) {
    const Z2 = e.Z, M2 = e.A, zz = Math.pow(Z1, 0.23) + Math.pow(Z2, 0.23);
    const eps = 32.53 * M2 * Ek / (Z1 * Z2 * (M1 + M2) * zz);
    let sn;
    if (eps <= 30) sn = Math.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * Math.pow(eps, 0.21226) + 0.19593 * Math.sqrt(eps)));
    else sn = Math.log(eps) / (2 * eps);
    const Sat = 8.462e-15 * Z1 * Z2 * M1 * sn / ((M1 + M2) * zz);      // eV cm2 / atom
    s += e.w / e.A * Sat * 1e-6;
  }
  return s * NA;
}
// Andersen-Ziegler (1977) proton coefficients, Z2 = 1..14: [A1..A5], E in keV, S in eV/(1e15 atoms/cm2)
const AZ = { 1: [1.262, 1.44, 242.6, 12000, 0.1159], 2: [1.229, 1.397, 484.5, 5873, 0.05225], 3: [1.411, 1.6, 725.6, 3013, 0.04578], 4: [2.248, 2.59, 966, 153.8, 0.03475],
  5: [2.474, 2.815, 1206, 1060, 0.02855], 6: [2.631, 2.601, 1701, 1279, 0.01638], 7: [2.954, 3.35, 1683, 1900, 0.02513], 8: [2.652, 3.0, 1920, 2000, 0.0230],
  9: [2.085, 2.352, 2157, 2634, 0.01816], 10: [1.951, 2.199, 2393, 2699, 0.01568], 11: [2.542, 3.0, 2628, 1854, 0.01472], 12: [3.791, 4.46, 2862, 1009, 0.01397],
  13: [4.154, 4.739, 2766, 164.5, 0.02023], 14: [4.914, 5.598, 3193, 232.7, 0.01419] };
function pAtomAZ(Z, Ek) {                         // eV/(1e15 atoms/cm2)
  const a = AZ[Z];
  if (Ek < 10) return a[0] * Math.sqrt(Ek);
  const sl = a[1] * Math.pow(Ek, 0.45), sh = a[2] / Ek * Math.log(1 + a[3] / Ek + a[4] * Ek);
  return sl * sh / (sl + sh);
}
const PSPE = { M: 938.27208816 };
function pAtomGeneric(e, Ep) {                    // generic proton stopping for one element (Z > 14), eV/(1e15 atoms/cm2)  (Ep in MeV)
  const T = Ep, M = PSPE.M, g = 1 + T / M, b2 = 1 - 1 / (g * g), bg2 = b2 * g * g, I = SPT.gI * e.Iraw * 1e-6;
  const u = 2 * ME * bg2 / I;
  let L = 0.5 * Math.log(1 + u * u) - b2; if (L < 0.02) L = 0.02;
  const shAtom = KC * (e.Z / e.A) / b2 * L * (e.A / NA);                   // MeV cm2 per atom
  const Sh = shAtom * 1e21;                                                   // MeV cm2 -> eV per (1e15 atoms/cm2)
  const vr = Math.sqrt(b2) / ALPHA, Zs = Math.pow(1 + Math.pow(e.Z, 2 / 3), 1.5);
  const Sl = SPT.gLS * 191.5e-16 * e.Z / Zs * vr * 1e15;                      // eV cm2 -> eV per 1e15 at/cm2
  return Math.pow(Math.pow(Sl, -SPT.gN) + Math.pow(Sh, -SPT.gN), -1 / SPT.gN);
}
function phaseFactor(mat, Ek) {                   // chemical/phase effect: condensed compounds stop less than Bragg sums of gas-phase data
  if (mat.phase === 'gas' || mat.elems.length < 2) return 1;
  const l = Math.log(Ek / SPT.phE0) / SPT.phW; return 1 - SPT.phK * Math.exp(-l * l);
}
function stopPLow(mat, Ep) {                      // proton stopping power at proton energy Ep (MeV), low-energy empirical branch, MeV cm2/g
  const Ek = Ep * 1000; let s = 0;
  for (const e of mat.elems) {
    const sa = AZ[e.Z] ? pAtomAZ(e.Z, Ek) : pAtomGeneric(e, Ep);
    s += e.w / e.A * sa;
  }
  return s * 1e-21 * NA * phaseFactor(mat, Ek);       // eV/(1e15 at/cm2) -> MeV cm2/atom: 1e-6*1e-15
}
function stopHeavyOld(sp, mat, T) {              // semi-empirical model (Andersen-Ziegler / Bethe-Bloch blend); used only to interpolate between tabulated elements
  const g = 1 + T / sp.M, Epr = (g - 1) * PSPE.M;                    // proton energy at the same velocity
  const z2 = effCharge2(sp, T);
  let Se;
  if (Epr <= SPT.E1) Se = z2 * stopPLow(mat, Epr);
  else if (Epr >= SPT.E2) Se = stopHeavyHigh(sp, mat, T);
  else {
    const w = Math.log(Epr / SPT.E1) / Math.log(SPT.E2 / SPT.E1);
    Se = (1 - w) * z2 * stopPLow(mat, Epr) + w * stopHeavyHigh(sp, mat, T);
  }
  return { Se, Sn: stopNuclear(sp, mat, T) };
}
