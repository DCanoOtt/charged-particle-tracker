// ---------- tabulated electronic stopping powers: NIST PSTAR (p), NIST ASTAR (alpha), ICRU 73 (Li..Ar) ----------
const SP_P = { kind: 'ion', Z1: 1, A: 1, M: PSPE.M }, ION18 = { kind: 'ion', Z1: 18, A: 40, M: 39.9623831 * AMU - 18 * ME };
const TCACHE = {}, EMATS = {};
let STOPSET = 'auto';                            // 'auto': newest tables (ICRU 90 / ICRU 73 update) where available, 'nist': PSTAR/ASTAR/ICRU 73, 'semi': semi-empirical model only
const RHO90 = { 276: 1.0, 104: 1.2048e-3, 906: 1.7 };
function id90(mat) {                               // material with an ICRU 90 table (liquid water, dry air, graphite)
  if (mat._id90 !== undefined) return mat._id90;
  const t = tabIdOf(mat); let id = null;
  if (t === 276 && mat.phase !== 'gas') id = 276; else if (t === 104) id = 104; else if (mat.elems.length === 1 && mat.elems[0].Z === 6 && mat.phase !== 'gas') id = 906;
  return (mat._id90 = id);
}
function tabCurve(kind, Z1, id) {                  // log-log curve with Hermite slopes, or null
  const key = kind + (kind === 'h' || kind === 'HN' ? Z1 : '') + ':' + id; if (key in TCACHE) return TCACHE[key];
  const E = { p: STAB.pE, a: STAB.aE, h: STAB.hE, HN: STAB.hE, P90: STAB.pE90, A90: STAB.aE90, C90: STAB.cE90, E90: STAB.eE90 }[kind];
  const V = kind === 'p' ? STAB.p[id] : kind === 'a' ? STAB.a[id] : kind === 'h' ? (STAB.h[Z1] ? STAB.h[Z1][id] : null) : kind === 'HN' ? (STAB.hn[Z1] ? STAB.hn[Z1][id] : null) : kind === 'P90' ? STAB.p90[id] : kind === 'A90' ? STAB.a90[id] : kind === 'C90' ? STAB.c90[id] : STAB.e90[id];
  if (!V) return (TCACHE[key] = null);
  const n = E.length, x = E.map(Math.log), y = V.map(v => Math.log(Math.max(v, 1e-30))), m = new Array(n);
  for (let i = 0; i < n; i++) { const a = Math.max(i - 1, 0), b = Math.min(i + 1, n - 1); m[i] = (y[b] - y[a]) / (x[b] - x[a]); }
  return (TCACHE[key] = { n, x, y, m, Emin: E[0], Emax: E[n - 1] });
}
function tabEval(c, E) {                           // MeV cm2/g at energy E (proton MeV or MeV/u); -1 above the table; below it S ~ v
  if (E > c.Emax * 1.0000001) return -1;
  if (E <= c.Emin) return Math.exp(c.y[0]) * Math.sqrt(E / c.Emin);
  const lx = Math.log(E); let lo = 0, hi = c.n - 1; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (c.x[mid] <= lx) lo = mid; else hi = mid; }
  const h = c.x[hi] - c.x[lo], t = (lx - c.x[lo]) / h, t2 = t * t, t3 = t2 * t;
  return Math.exp((2 * t3 - 3 * t2 + 1) * c.y[lo] + (t3 - 2 * t2 + t) * h * c.m[lo] + (-2 * t3 + 3 * t2) * c.y[hi] + (t3 - t2) * h * c.m[hi]);
}
function elemMat(Z) { if (!EMATS[Z]) EMATS[Z] = buildMaterial({ mode: 'atoms', rows: [{ kind: 'el', sym: ELTAB[Z].s, amt: 1 }] }, null); return EMATS[Z]; }
function tabZlist(kind, Z1) {
  const key = 'Z' + kind + (kind === 'h' ? Z1 : ''); if (TCACHE[key]) return TCACHE[key];
  const src = kind === 'p' ? STAB.p : kind === 'a' ? STAB.a : (STAB.h[Z1] || {}); return (TCACHE[key] = Object.keys(src).map(Number).filter(i => i <= 98).sort((a, b) => a - b));
}
function elemTab(kind, spx, Tx, Z, E) {            // mass stopping power of the element Z, MeV cm2/g (-1 above the table)
  const c = tabCurve(kind, spx.Z1, Z); if (c) return tabEval(c, E);
  const zl = tabZlist(kind, spx.Z1); let a = 0, b = 0; for (const z of zl) { if (z < Z) a = z; if (z > Z && !b) b = z; }
  const old = z => stopHeavyOld(spx, elemMat(z), Tx).Se;
  const rat = z => { const s = tabEval(tabCurve(kind, spx.Z1, z), E); return s < 0 ? -1 : Math.log(s / old(z)); };
  let lr;
  if (a && b) { const ra = rat(a), rb = rat(b); if (ra === -1 && rb === -1) return -1; lr = ra + (rb - ra) * (Z - a) / (b - a); }
  else { const r = rat(a || b); if (r === -1) return -1; lr = r; }
  return Math.exp(lr) * old(Z);
}
function tabIdOf(mat) {                            // ICRU/NIST material id whose composition equals this material (or its single element)
  if (mat._tabId !== undefined) return mat._tabId;
  let id = null;
  if (mat.elems.length === 1) id = mat.elems[0].Z;
  else {
    const w = {}; mat.elems.forEach(e => w[e.Z] = e.w); let best = null;
    for (const k in STAB.comp) {
      const cm = {}; STAB.comp[k].forEach(([z, f]) => cm[z] = f); let ok = true;
      for (const z of new Set([...Object.keys(w), ...Object.keys(cm)])) if (Math.abs((w[z] || 0) - (cm[z] || 0)) > 0.0025) { ok = false; break; }
      if (ok) { const gas = +k === 277, isg = mat.phase === 'gas'; if (best === null || (gas === isg)) best = +k; }
    }
    id = best;
  }
  return (mat._tabId = id);
}
function matTab(kind, spx, Tx, mat, E) {           // Bragg sum / compound table, MeV cm2/g (-1 above the table)
  if (STOPSET === 'auto') {
    const i9 = id90(mat);
    if (i9) { const k9 = spx.Z1 === 1 && kind === 'p' ? 'P90' : spx.Z1 === 2 && kind === 'a' ? 'A90' : spx.Z1 === 6 && kind === 'h' ? 'C90' : null; const c9 = k9 ? tabCurve(k9, 0, i9) : null; if (c9) return tabEval(c9, E); }
    if (kind === 'h' && (tabIdOf(mat) === 276 || tabIdOf(mat) === 104) && STAB.hn[spx.Z1]) { const cn = tabCurve('HN', spx.Z1, tabIdOf(mat)); if (cn) return tabEval(cn, E); }
  }
  const id = tabIdOf(mat);
  if (id !== null) { const c = tabCurve(kind, spx.Z1, id); if (c) return tabEval(c, E); }
  let s = 0; for (const e of mat.elems) { const v = elemTab(kind, spx, Tx, e.Z, E); if (v < 0) return -1; s += e.w * v; }
  return s;
}
function seTab(sp, mat, T) {
  if (STOPSET === 'semi') return stopHeavyOld(sp, mat, T).Se;
  const g1 = T / sp.M, Eu = g1 * AMU, Ep = g1 * PSPE.M;
  let kind, E, spx = sp, Tx = T, mul = 1;
  if (sp.Z1 === 1) { kind = 'p'; E = Ep; }
  else if (sp.Z1 === 2) { kind = 'a'; E = Eu; }
  else if (sp.Z1 <= 18) { kind = 'h'; E = Eu; }
  else {                                           // Z1 > 18: Brandt-Kitagawa scaling from the tabulated Ar ion at the same reduced velocity v / Z1^(2/3)
    kind = 'p'; E = Ep; spx = SP_P; Tx = Ep;
    const Z1 = sp.Z1, EuA = Eu * Math.pow(18 / Z1, 4 / 3);
    if (EuA >= STAB.hE[STAB.hE.length - 1]) mul = Z1 * Z1;
    else { const ar = ION18, sI = matTab('h', ar, EuA / AMU * ar.M, mat, EuA), sPp = matTab('p', SP_P, EuA * PSPE.M / AMU, mat, EuA * PSPE.M / AMU); mul = (Z1 / 18) * (Z1 / 18) * (sI > 0 && sPp > 0 ? sI / sPp : 324); }
  }
  const s = matTab(kind, spx, Tx, mat, E);
  if (s >= 0) return mul * s;
  // above the table: Bethe-Bloch, matched to the table at its upper edge
  const key = STOPSET + kind + sp.Z1 + ':' + (id90(mat) || 0); mat._cf = mat._cf || {};
  if (mat._cf[key] === undefined) {
    const Eb = (kind === 'p' ? STAB.pE : kind === 'a' ? STAB.aE : STAB.hE).slice(-1)[0]; let Em = Eb, tab = -1;
    for (let k = 0; k < 80 && tab < 0; k++, Em *= 0.97) { const gm0 = Em / (kind === 'p' ? PSPE.M : AMU); tab = matTab(kind, spx, Tx === T ? gm0 * sp.M : Em, mat, Em); }
    Em /= 0.97; const gm = Em / (kind === 'p' ? PSPE.M : AMU), Tm = gm * sp.M, hi = stopHeavyHigh(sp, mat, Tm);
    mat._cf[key] = tab > 0 && hi > 0 ? tab * mul / hi : 1;
  }
  return mat._cf[key] * stopHeavyHigh(sp, mat, T);
}
function zEff2(sp, mat, T) {
  if (STOPSET === 'semi') return effCharge2(sp, T);                       // effective charge squared (ratio of the ion's to the proton's stopping power at equal velocity), for straggling and delta rays
  if (sp.Z1 === 1) return 1;
  const Ep = T / sp.M * PSPE.M, sPp = matTab('p', SP_P, Ep, mat, Ep);
  if (sPp <= 0) return effCharge2(sp, T);
  return seTab(sp, mat, T) / sPp;
}
function stopHeavy(sp, mat, T) { return { Se: seTab(sp, mat, T), Sn: stopNuclear(sp, mat, T) }; }
