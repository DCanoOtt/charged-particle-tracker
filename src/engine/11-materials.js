// ---------- materials ----------
// def: { name, phase:'solid'|'liquid'|'gas', mode:'atoms'|'mass', rows:[{kind:'el',sym,amt}|{kind:'mat',key,amt}], rho, rhoAuto:bool,
//        gas:{P(kPa),T(K)}, I (eV)|null, iso:{sym:[[A,%],..]} }
function buildMaterial(def, lib, depth) {
  depth = depth || 0; if (depth > 6) throw new Error('Material nesting too deep');
  const rows = (def.rows || []).filter(r => r && +r.amt > 0 && (r.kind === 'mat' ? r.key : r.sym));
  if (!rows.length) throw new Error('The material has no components.');
  const mode = def.mode === 'mass' ? 'mass' : 'atoms';
  const wZ = {}, parts = [];       // wZ: element mass fraction; parts: {w, ZA, lnI, rho}
  let rhoMix = 0, formulaMass = 0;
  if (mode === 'atoms') {
    let tot = 0; const m = {};
    for (const r of rows) {
      if (r.kind === 'mat') throw new Error('Materials can only be mixed in mass-fraction mode.');
      const Z = SYMZ[r.sym]; if (!Z) throw new Error('Unknown element: ' + r.sym);
      const A = elemMean(Z, def.iso && def.iso[r.sym] ? isoListOf(Z, def.iso) : null);
      m[Z] = (m[Z] || 0) + r.amt * A; tot += r.amt * A;
    }
    let fm = tot;
    { // amounts such as Ar 100 or C 100 / O 200 are ratios: reduce to the smallest integer formula
      const am = rows.map(r => +r.amt), gcd = (a, b) => b < 1e-9 ? a : gcd(b, a % b);
      if (am.every(v => Math.abs(v - Math.round(v)) < 1e-9)) { const g = am.reduce((a, v) => gcd(a, v)); if (g >= 3) fm = tot / g; }
    }
    formulaMass = fm; for (const Z in m) wZ[Z] = m[Z] / tot;
  } else {
    let tot = 0; for (const r of rows) tot += +r.amt;
    for (const r of rows) {
      const w = r.amt / tot;
      if (r.kind === 'mat') {
        const sd = lib && lib[r.key]; if (!sd) throw new Error('Unknown material: ' + r.key);
        const sm = buildMaterial(sd, lib, depth + 1);
        for (const el of sm.elems) wZ[el.Z] = (wZ[el.Z] || 0) + w * el.w;
        parts.push({ w, ZA: sm.ZA, lnI: Math.log(sm.I), rho: sm.rho });
      } else {
        const Z = SYMZ[r.sym]; if (!Z) throw new Error('Unknown element: ' + r.sym);
        wZ[Z] = (wZ[Z] || 0) + w;
        const A = elemMean(Z, def.iso && def.iso[r.sym] ? isoListOf(Z, def.iso) : null);
        parts.push({ w, ZA: Z / A, lnI: Math.log(IVAL[Z]), rho: ELDENS[ELTAB[Z].s] || NaN, Z });
      }
    }
    let inv = 0, ok = true; for (const p of parts) { if (!(p.rho > 0)) ok = false; inv += p.w / p.rho; }
    rhoMix = ok && inv > 0 ? 1 / inv : NaN;
  }
  // element list
  const Zs = Object.keys(wZ).map(Number).sort((a, b) => wZ[b] - wZ[a]);
  const elems = []; let sumN = 0;
  for (const Z of Zs) {
    const il = isoListOf(Z, def.iso), A = elemMean(Z, il);
    elems.push({ Z, sym: ELTAB[Z].s, w: wZ[Z], A, iso: il.map(i => ({ A: i[0], ab: i[1] / 100, mass: isoMass(Z, i[0]) })), Iraw: IVAL[Z] });
    sumN += wZ[Z] / A;
  }
  let ZA = 0; for (const e of elems) { e.x = (e.w / e.A) / sumN; ZA += e.w * e.Z / e.A; }
  const Abar = 1 / sumN;
  // Bragg additivity of ln I (by electron density)
  let num = 0, den = 0;
  if (mode === 'atoms') { for (const e of elems) { const q = e.w * e.Z / e.A; num += q * Math.log(e.Iraw); den += q; } }
  else { for (const p of parts) { const q = p.w * p.ZA; num += q * p.lnI; den += q; } }
  const Ibragg = Math.exp(num / den);
  const I = def.I > 0 ? +def.I : Ibragg;
  // density
  const phase = def.phase || 'solid';
  let rho = NaN, rhoSrc = 'user';
  const ru = +def.rho;
  if (phase === 'gas' && def.rhoAuto !== false) {
    const g = def.gas || {}, P = (g.P > 0 ? g.P : 101.325), T = (g.T > 0 ? g.T : 293.15);
    const M = def.molMass > 0 ? +def.molMass : (mode === 'atoms' ? formulaMass : Abar);
    rho = P * 1000 * M / (RGAS * T) / 1e6; rhoSrc = 'ideal gas';
  } else if (ru > 0) rho = ru;
  else if (mode === 'mass' && rhoMix > 0) { rho = rhoMix; rhoSrc = 'ideal mixing'; }
  else if (elems.length === 1 && ELDENS[elems[0].sym] > 0) { rho = ELDENS[elems[0].sym]; rhoSrc = 'element'; }
  if (!(rho > 0)) throw new Error('Enter the density of the material.');
  const Ifac = I / Ibragg;
  for (const e of elems) { e.I = e.Iraw * Ifac; e.n = rho * NA * e.w / e.A; }
  const nAt = rho * NA / Abar, nEl = rho * NA * ZA;
  // density effect (Sternheimer-Peierls parametrisation, Sternheimer-Berger-Seltzer 1984 defaults)
  const hwp = 28.8159 * Math.sqrt(rho * ZA);
  const Cb = 1 + 2 * Math.log(I / hwp);
  let x0, x1;
  if (phase === 'gas') {
    if (Cb < 10) x0 = 1.6; else if (Cb < 10.5) x0 = 1.7; else if (Cb < 11) x0 = 1.8; else if (Cb < 11.5) x0 = 1.9; else if (Cb < 12.25) x0 = 2.0; else if (Cb < 13.804) x0 = 2.0; else x0 = 0.326 * Cb - 2.5;
    x1 = Cb < 12.25 ? 4.0 : 5.0;
  } else if (I < 100) { x1 = 2.0; x0 = Cb <= 3.681 ? 0.2 : 0.326 * Cb - 1.0; }
  else { x1 = 3.0; x0 = Cb <= 5.215 ? 0.2 : 0.326 * Cb - 1.5; }
  const m = 3, a = (Cb - 2 * LN10 * x0) / Math.pow(x1 - x0, m);
  // radiation length (Tsai), g/cm2
  let invX0 = 0;
  for (const e of elems) {
    const Z = e.Z, aa = ALPHA * Z;
    const fZ = aa * aa * (1 / (1 + aa * aa) + 0.20206 - 0.0369 * aa * aa + 0.0083 * aa ** 4 - 0.002 * aa ** 6);
    const LR = Z === 1 ? 5.31 : Z === 2 ? 4.79 : Z === 3 ? 4.74 : Z === 4 ? 4.71 : Math.log(184.15 * Math.pow(Z, -1 / 3));
    const LP = Z === 1 ? 6.144 : Z === 2 ? 5.621 : Z === 3 ? 5.805 : Z === 4 ? 5.924 : Math.log(1194 * Math.pow(Z, -2 / 3));
    invX0 += 4 * ALPHA * 2.8179403262e-13 ** 2 * NA * (e.w / e.A) * (Z * Z * (LR - fZ) + Z * LP);
  }
  // shell-correction weights (electron fractions)
  for (const e of elems) e.fel = e.w * e.Z / e.A / ZA;
  // radiative stopping constants for electrons: sum n Z(Z+1.2), phi_inf
  let sRadNum = 0, phiNum = 0;
  for (const e of elems) {
    const q = (e.w / e.A) * e.Z * (e.Z + 1.2);
    sRadNum += q; phiNum += q * (4 * Math.log(183 * Math.pow(e.Z, -1 / 3)) + 2 / 9);
  }
  return {
    name: def.name || 'Material', phase, rho, rhoSrc, elems, Abar, ZA, nAt, nEl, I, Ibragg, hwp,
    den: { Cb, x0, x1, a, m }, X0: 1 / invX0, formulaMass, mode,
    radK: sRadNum, phiInf: phiNum / sRadNum, Zbar: ZA * Abar, vac: false
  };
}

function densityEffect(mat, bg) {
  const x = Math.log(bg) / LN10, p = mat.den;
  if (x < p.x0) return 0;
  let d = 2 * LN10 * x - p.Cb;
  if (x < p.x1) d += p.a * Math.pow(p.x1 - x, p.m);
  return d > 0 ? d : 0;
}
