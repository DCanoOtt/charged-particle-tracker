// ---------- species ----------
function makeSpecies(s) {
  if (s.kind === 'e') return { kind: 'e', Z1: 1, A: 0, M: ME, name: 'e⁻', isE: true };
  const Z = s.Z, A = s.A;
  const m = isoMass(Z, A);
  const M = m * AMU - Z * ME;
  return { kind: 'ion', Z1: Z, A, M, name: (ELTAB[Z] ? ELTAB[Z].s : 'Z' + Z) + '-' + A, isE: false };
}
