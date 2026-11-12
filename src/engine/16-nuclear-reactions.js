// ---------- nuclear reaction model (simple; replace by databases later) ----------
function sigmaReaction(Zp, Ap, Zt, At, T) {      // mb
  const cp = Math.cbrt(Ap), ct = Math.cbrt(At);
  const Ecm = T * At / (Ap + At);
  const Vb = E2 * Zp * Zt / (1.2 * (cp + ct) + 2.0);
  let s;
  if (Ap === 1 && Zp === 1) {
    if (At === 1) return T > 290 ? 25 : 0;
    const base = 45 * Math.pow(At, 0.7) * (1 + 0.016 * Math.sin(5.3 - 2.63 * Math.log(At)));
    const Te = Math.max(T, 10);
    const f = Te < 2000 ? 1 - 0.62 * Math.exp(-Te / 200) * Math.sin(10.9 * Math.pow(Te, -0.28)) : 1;
    s = base * f;
    if (T < 10) s *= Math.max(0, 1 - Vb / Ecm) / Math.max(1 - Vb / (10 * At / (Ap + At)), 1e-3);
    else s *= 1;
    return Math.max(s, 0);
  }
  const b = Ap === 2 ? 1.0 : Ap === 3 ? 0.6 : 0.2;
  s = Math.PI * 1.36 * 1.36 * 10 * Math.pow(cp + ct - b, 2);        // fm2 -> mb (1 fm2 = 10 mb)
  const bar = 1 - Vb / Ecm; if (!(bar > 0)) return 0;
  return s * bar;
}

// ---------- evaluated cross sections (p, d, t, ³He, α) (TENDL, ENDF-6) used instead of the model when available ----------
// XSLIB['Z-A'] = { E: Float64Array (MeV), s: Float64Array (non-elastic cross section, barn) }
let XSLIB = null, XSLIB_SP = { Z: 2, A: 4 };    // library of evaluated cross sections and the projectile (Z, A) it belongs to
function sigmaIso(sp, Z, A, T) {            // mb
  if (XSLIB && sp.Z1 === XSLIB_SP.Z && sp.A === XSLIB_SP.A) {
    const t = XSLIB[Z + '-' + A];
    if (t) {
      const E = t.E, s = t.s, n = E.length;
      if (T <= E[0]) return 0; if (T >= E[n - 1]) return s[n - 1] * 1000;
      let lo = 0, hi = n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (E[m] > T) hi = m; else lo = m; }
      return (s[lo] + (s[hi] - s[lo]) * (T - E[lo]) / (E[hi] - E[lo])) * 1000;
    }
  }
  return sigmaReaction(sp.Z1, sp.A, Z, A, T);
}
