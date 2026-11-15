/* global ELTAB */
// ===================== Physics: constants, materials, stopping powers =====================
export const ME = 0.51099895,
  AMU = 931.49410242,
  ALPHA = 1 / 137.035999084,
  NA = 6.02214076e23;
export const KC = 0.307075; // 4 pi NA re^2 me c^2, MeV cm2 / mol
export const A0FM = 52917.7210903; // Bohr radius, fm
export const HBARC = 197.3269804; // MeV fm
export const E2 = 1.43996454; // e^2, MeV fm
export const RGAS = 8.314462618;
export const LN10 = Math.LN10;
export const SYMZ = {};
for (const z in ELTAB) SYMZ[ELTAB[z].s] = +z;
export const SYMS = Object.keys(ELTAB).map(z => ELTAB[z].s);
// mean excitation energies (eV), ICRU 37 / NIST, index = Z
export const IVAL = [
  0, 19.2, 41.8, 40, 63.7, 76, 78, 82, 95, 115, 137, 149, 156, 166, 173, 173, 180, 174, 188, 190, 191, 216, 233, 245, 257, 272, 286, 297, 311, 322, 330, 334,
  350, 347, 348, 343, 352, 363, 366, 379, 393, 417, 424, 428, 441, 449, 470, 470, 469, 488, 488, 487, 485, 491, 482, 488, 491, 501, 523, 535, 546, 560, 574,
  580, 591, 614, 628, 650, 658, 674, 684, 694, 705, 718, 727, 736, 746, 757, 790, 790, 800, 810, 823, 823, 830, 825, 794, 827, 826, 841, 847, 878, 890
];
const EXTRA_MASS = { '1-3': 3.01604928 }; // tritium (not a natural isotope)

export function isoMass(Z, A) {
  const e = ELTAB[Z];
  if (e) for (const i of e.iso) if (i[0] === A) return i[2];
  if (Z === 1 && A === 3) return EXTRA_MASS['1-3'];
  return A * 1.0; // fallback: A u
}
export function elemMean(Z, isoList) {
  // mean atomic mass (g/mol) for an isotope list [[A, %]]
  const e = ELTAB[Z];
  if (!isoList) return e.w;
  let s = 0,
    m = 0;
  for (const [A, ab] of isoList) {
    s += ab;
    m += ab * isoMass(Z, A);
  }
  return s > 0 ? m / s : e.w;
}
export function isoListOf(Z, ov) {
  const e = ELTAB[Z];
  const l = ov && ov[e.s] ? ov[e.s] : e.iso.map(i => [i[0], i[1]]);
  let s = 0;
  for (const i of l) s += i[1];
  if (!(s > 0)) return e.iso.map(i => [i[0], i[1]]);
  return l.map(i => [i[0], (i[1] / s) * 100]);
}
