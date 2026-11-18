// @ts-check
/* ====== X-rays from atomic vacancies (ENSDF decay data) ======
   A decay leaves vacancies in the K and L shells of the daughter atom: electron capture (E records: CK, CL, CM+) and internal conversion (G records: KC, LC, MC, or the total CC).
   Each vacancy emits X-rays with the yields of data/xray.js (xraylib: fluorescence yield × radiative rate); K vacancies create L vacancies (K→L X-rays and Auger electrons);
   L vacancies are redistributed by Coster–Kronig transitions before they emit. Shells outside K and L (M…) are not followed. */
const XR_EC_DEF = { K: 0.86, L: 0.12 }; // electron capture without CK/CL in the data: shares of the captures from the K and L shells (approximate)
const XR_AUGER_L = [0.5, 0.5, 0.7]; // L1, L2, L3 vacancies left by the Auger decay of one K vacancy (1.7 in total)
/** @param {number} ck @param {number} cl @returns {Vacancies} */
function xrEcSplit(ck, cl) {
  // capture probabilities → vacancies per capture {K, L1, L2, L3}
  const k = ck === ck ? ck : XR_EC_DEF.K,
    l = cl === cl ? cl : ck === ck ? Math.max(0, 1 - ck) * 0.85 : XR_EC_DEF.L;
  return { K: k, L1: l * 0.9, L2: l * 0.1, L3: 0 };
}
/** @param {number} Z @param {number} Eg @param {number} kc @param {number} lc @param {number} cc @param {number=} l1c @param {number=} l2c @param {number=} l3c @returns {Vacancies | null} */
function xrIcSplit(Z, Eg, kc, lc, cc, l1c, l2c, l3c) {
  // conversion coefficients → vacancies per conversion-electron-emitting transition {K, L1, L2, L3} (as coefficients, to be divided by 1+α)
  const d = typeof XRAY !== 'undefined' ? XRAY[Z] : null;
  if (!d || !(cc > 0)) return null;
  const EK = d.e[0],
    EL = d.e[1],
    o = { K: 0, L1: 0, L2: 0, L3: 0 };
  const ok = v => typeof v === 'number' && v === v;
  if (ok(kc) || ok(lc) || ok(l1c) || ok(l2c) || ok(l3c)) {
    if (ok(kc)) o.K = kc;
    if (ok(l1c) || ok(l2c) || ok(l3c)) {
      o.L1 = l1c || 0;
      o.L2 = l2c || 0;
      o.L3 = l3c || 0;
    } else if (ok(lc)) {
      o.L1 = lc * 0.5;
      o.L2 = lc * 0.3;
      o.L3 = lc * 0.2;
    }
    return o;
  }
  if (Eg > EK) {
    o.K = cc * 0.8;
    o.L1 = cc * 0.1;
    o.L2 = cc * 0.05;
    o.L3 = cc * 0.02;
  } else if (Eg > EL) {
    o.L1 = cc * 0.5;
    o.L2 = cc * 0.3;
    o.L3 = cc * 0.2;
  } // total CC only: typical shell shares (K ~ 80 % of α above the K edge)
  return o;
}
/** @param {number} Z @param {Partial<Vacancies>} v @returns {[number, number][]} */
function xrEmit(Z, v) {
  // primary vacancies per decay {K, L1, L2, L3} → X-ray lines [[E (keV), photons per decay]]
  const d = typeof XRAY !== 'undefined' ? XRAY[Z] : null;
  if (!d) return [];
  let nK = v.K || 0,
    n1 = v.L1 || 0,
    n2 = v.L2 || 0,
    n3 = v.L3 || 0;
  if (!(nK + n1 + n2 + n3 > 0)) return [];
  const au = Math.max(0, 1 - d.w[0]);
  n1 += nK * au * XR_AUGER_L[0];
  n2 += nK * (d.kl[0] + au * XR_AUGER_L[1]);
  n3 += nK * (d.kl[1] + au * XR_AUGER_L[2]);
  const m1 = n1,
    m2 = n2 + n1 * d.f[0],
    m3 = n3 + n1 * d.f[1] + m2 * d.f[2]; // Coster–Kronig: L1→L2, L1→L3, L2→L3
  /** @type {[number, number][]} */ const out = [];
  /** @type {[[number, number][], number][]} */ const parts = [
    [d.K, nK],
    [d.L1, m1],
    [d.L2, m2],
    [d.L3, m3]
  ];
  parts.forEach(([ls, n]) => {
    if (n > 0) for (const l of ls) if (l[1] * n > 1e-9) out.push([l[0], l[1] * n]);
  });
  return out.sort((a, b) => a[0] - b[0]);
}
/** @param {[number, number][]} lines @returns {[number, number][]} */
function xrMerge(lines) {
  // sums lines of the same energy (several datasets / transitions)
  const m = new Map();
  for (const l of lines) {
    const k = l[0].toFixed(3);
    const o = m.get(k);
    if (o) o[1] += l[1];
    else m.set(k, [l[0], l[1]]);
  }
  return Array.from(m.values()).sort((a, b) => a[0] - b[0]);
}
