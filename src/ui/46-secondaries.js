/* ====== energy and angular distributions of the secondary particles (ENDF-6 MF 6) ======
   For every layer and secondary species the reaction rate in the energy bins of the transported ion flux, p_k = n σ(E_k) L_k y(E_k), is
   distributed with the MF 6 data of each product:  LAW 1 (continuum, Legendre or Kalbach-Mann), LAW 2 (two-body, discrete level)
   and LAW 6 (N-body phase space).  Data given in the centre-of-mass system are transformed to the laboratory by the velocity of the CM:
   with E_c = A_p·A_proj·E/(A_proj+A_target)²,  E_cm' = (√E − √E_c·…)  and f_lab(E,Ω) = f_cm(E',Ω')·√(E/E') (phase-space invariance).
   Results are per primary: dN/dE (all angles, per MeV) and d²N/(dE dΩ) at the lab angles in SEC_ANGLES (per MeV per sr). */
const SEC_ANGLES = [0, 5, 10, 30, 60, 120, 160];

const SEC = new Map(), SECST = new Map();   // SEC: data of the plots on screen; SECST: their view state (zoom, log axis, hidden curves)
function secLeg(l, mu, st) {                // Legendre polynomial P_l(mu) by recurrence
  let p0 = 1, p1 = mu; if (l === 0) return 1; for (let k = 1; k < l; k++) { const p2 = ((2 * k + 1) * mu * p1 - k * p0) / (k + 1); p0 = p1; p1 = p2; } return p1;
}
function secLoc(x, v) { const n = x.length; let lo = 0, hi = n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (x[m] > v) hi = m; else lo = m; } return lo; }
// MF6 LAW 1: f(E', mu) of the table of one incident energy (per MeV and per unit mu); mu === null returns the energy distribution f0(E')
function d1List(d, L, Ep, mu, aCtx) {
  const x = L.Ep, n = x.length; if (n < 2 || Ep < x[0] || Ep > x[n - 1]) return 0;
  const lo = secLoc(x, Ep), hi = Math.min(n - 1, lo + 1), t = d.lep === 1 || x[hi] === x[lo] ? 0 : (Ep - x[lo]) / (x[hi] - x[lo]), g = a => a[lo] + (a[hi] - a[lo]) * t;
  const f0 = g(L.f0); if (mu === null) return f0;
  if (d.lang === 1 && !L.p.length) return 0.5 * f0;                      // no angular coefficients: isotropic
  if (d.lang === 1) { let s = 0.5 * f0; for (let l = 1; l <= L.p.length; l++) s += (2 * l + 1) / 2 * g(L.p[l - 1]) * secLeg(l, mu); return Math.max(s, 0); }
  const r = L.p.length ? g(L.p[0]) : 0, a = L.p[1] ? g(L.p[1]) : (aCtx ? aCtx(Ep) : 0);   // NA = 1: a from the Kalbach systematics
  if (Math.abs(a) < 1e-6) return 0.5 * f0;
  return Math.max(0.5 * f0 * a / Math.sinh(a) * (Math.cosh(a * mu) + r * Math.sinh(a * mu)), 0);
}
function d1Eval(d, Ein, Ep, mu, ac) {
  const E = d.E, n = E.length; if (!n) return 0; if (Ein <= E[0]) return d1List(d, d.L[0], Ep, mu, ac); if (Ein >= E[n - 1]) return d1List(d, d.L[n - 1], Ep, mu, ac);
  const lo = secLoc(E, Ein), t = (Ein - E[lo]) / (E[lo + 1] - E[lo]); return (1 - t) * d1List(d, d.L[lo], Ep, mu, ac) + t * d1List(d, d.L[lo + 1], Ep, mu, ac);
}
// MF6 LAW 2: normalised angular pdf f(mu)
function d2List(L, mu) {
  if (L.lang === 0) { let s = 0.5; for (let l = 1; l <= L.a.length; l++) s += (2 * l + 1) / 2 * L.a[l - 1] * secLeg(l, mu); return Math.max(s, 0); }
  const x = L.mu, n = x.length; if (n < 2 || mu < x[0] || mu > x[n - 1]) return 0; const lo = secLoc(x, mu), hi = Math.min(n - 1, lo + 1); return x[hi] === x[lo] ? L.f[lo] : L.f[lo] + (L.f[hi] - L.f[lo]) * (mu - x[lo]) / (x[hi] - x[lo]);
}
function d2Eval(d, Ein, mu) {
  const E = d.E, n = E.length; if (!n) return 0.5; if (Ein <= E[0]) return d2List(d.L[0], mu); if (Ein >= E[n - 1]) return d2List(d.L[n - 1], mu);
  const lo = secLoc(E, Ein), t = (Ein - E[lo]) / (E[lo + 1] - E[lo]); return (1 - t) * d2List(d.L[lo], mu) + t * d2List(d.L[lo + 1], mu);
}
// Kalbach-Mann slope a(E') when the file gives only f0 and r (ENDF-6 NA = 1): Kalbach systematics for an incident alpha (M_a = 0),
// with the separation energies from a liquid-drop binding energy and S_b = S_a - Q_M
function wzB(Z, A) { const N = A - Z, pr = (Z % 2 === 0 && N % 2 === 0) ? 1 : (Z % 2 === 1 && N % 2 === 1) ? -1 : 0; return 15.75 * A - 17.8 * Math.pow(A, 2 / 3) - 0.711 * Z * (Z - 1) / Math.cbrt(A) - 23.7 * (A - 2 * Z) * (A - 2 * Z) / A + pr * 34 / Math.pow(A, 0.75); }
function kmSlope(tgt, ch, sp, E, AWR) {
  const P = PJ(), Sa = wzB(tgt.Z + P.Z, tgt.A + P.A) - wzB(tgt.Z, tgt.A) - P.Bp, Sb = Sa - (ch.QM || 0), ea = E * AWR / (AWR + P.awi) + Sa, mb = { n: 0.5, a: 2 }[sp] || 1;
  return Ep => { const eb = Ep + Sb; if (!(ea > 0 && eb > 0)) return 0; const X1 = Math.min(ea, 130) * eb / ea, X3 = Math.min(ea, 41) * eb / ea; return 0.04 * X1 + 1.8e-6 * X1 * X1 * X1 + 6.7e-7 * P.Ma * mb * Math.pow(X3, 4); };
}
// discrete photon lines of a LAW 1 product at incident energy Ein: { E: line energies (MeV), f: probabilities }
function d1Disc(d, Ein) {
  const E = d.E, n = E.length; if (!n) return null; let lo = 0, t = 0;
  if (Ein >= E[n - 1]) lo = n - 1; else if (Ein > E[0]) { lo = secLoc(E, Ein); t = (Ein - E[lo]) / (E[lo + 1] - E[lo]); }
  const A = d.L[lo].disc, B = t > 0 ? d.L[lo + 1].disc : null; if (!A) return null;
  if (B && B.E.length === A.E.length && A.E.every((e, i) => Math.abs(e - B.E[i]) <= 1e-9 * Math.max(e, 1e-9))) return { E: A.E, f: A.f.map((v, i) => v + (B.f[i] - v) * t) };
  return t > 0.5 && B ? B : A;
}
// photons: lines + continuum of LAW 1, taken as isotropic (no Doppler shift)
function secPhoton(d, E, W, tot, nbo, dEo) {
  const D = d1Disc(d, E); if (D) for (let i = 0; i < D.E.length; i++) { const j = Math.floor(D.E[i] / dEo); if (j >= 0 && j < nbo) tot[j] += W * D.f[i] / dEo; }
  for (let j = 0; j < nbo; j++) tot[j] += W * d1Eval(d, E, (j + 0.5) * dEo, null);
}
// N-body phase space (LAW 6): energy pdf in the CM, isotropic
function d6Make(d, Ap, q) {
  const Em = (d.apsx - Ap) / d.apsx * q; if (!(Em > 0)) return null; const ex = 1.5 * d.n - 4;
  const raw = e => e <= 0 || e >= Em ? 0 : Math.sqrt(e) * Math.pow(Em - e, ex); let I = 0; const N = 400; for (let i = 0; i < N; i++) I += raw((i + 0.5) * Em / N) * Em / N;
  return { Em, f0: e => I > 0 ? raw(e) / I : 0 };
}

// emit the products of one incident-energy block (weight W reactions, mean energy E) into the histograms tot[nbo] and ang[i][nbo]
function secEmit(pr, ch, tgt, E, W, tot, ang, nbo, dEo, cosT, sinT, kin) {
  const AWR = tgt.AWR;
  if (pr.sp === 'g') { if (pr.dist.law === 1) secPhoton(pr.dist, E, W, tot, nbo, dEo); return; }
  const d = pr.dist, Ap = pr.awp || 1, Ec = Ap * PJ().awi * E / Math.pow(PJ().awi + AWR, 2), w = Math.sqrt(Ec), Ecm = E * AWR / (AWR + PJ().awi), TWO_PI = 2 * Math.PI, nA = cosT.length;
  if (d.law === 2) {
    const q = Ecm + ch.QI; if (!(q > 0)) return; const AR = AWR + PJ().awi - Ap, E0 = AR / (AR + Ap) * q, u = Math.sqrt(E0); if (!(u > 0)) return;
    for (let i = 0; i < nA; i++) {
      const c = cosT[i], sn = sinT[i], disc = u * u - w * w * sn * sn; if (disc <= 1e-14) continue; const sq = Math.sqrt(disc);
      for (const sg of [1, -1]) {
        const s = w * c + sg * sq; if (!(s > 1e-9)) continue; const j = Math.floor(s * s / dEo); if (j >= nbo) continue;
        const mu = Math.max(-1, Math.min(1, (s * c - w) / u)); ang[i][j] += W * d2Eval(d, E, mu) / TWO_PI * s * s / (u * sq) / dEo;
      }
    }
    const N = 64; for (let m = 0; m < N; m++) { const mu = -1 + (m + 0.5) * 2 / N, j = Math.floor((E0 + Ec + 2 * u * w * mu) / dEo); if (j >= 0 && j < nbo) tot[j] += W * d2Eval(d, E, mu) * (2 / N) / dEo; }
    return;
  }
  let fn, f0fn;
  if (d.law === 1) { const ac = d.lang === 2 ? kmSlope(tgt, ch, pr.sp, E, AWR) : null; fn = (Ep, mu) => d1Eval(d, E, Ep, mu, ac); f0fn = Ep => d1Eval(d, E, Ep, null); }
  else if (d.law === 6) { const q = Ecm + ch.QM, ps = d6Make(d, Ap, q); if (!ps) return; f0fn = ps.f0; fn = (Ep) => ps.f0(Ep) / 2; }
  else return;
  const lab = d.law === 1 && d.frame === 'lab' && kin === 'table';   // 'twobody': E' is the CM energy and the lab spectrum follows from the kinematics
  for (let j = 0; j < nbo; j++) {
    const Eo = (j + 0.5) * dEo, s = Math.sqrt(Eo);
    if (lab) { tot[j] += W * f0fn(Eo); for (let i = 0; i < nA; i++) ang[i][j] += W * fn(Eo, cosT[i]) / TWO_PI; continue; }
    for (let i = 0; i < nA; i++) { const par = s * cosT[i] - w, perp = s * sinT[i], Ep = par * par + perp * perp; if (Ep > 0) ang[i][j] += W * fn(Ep, par / Math.sqrt(Ep)) / TWO_PI * Math.sqrt(Eo / Ep); }
    let sum = 0; const N = 48;
    for (let m = 0; m < N; m++) { const c = -1 + (m + 0.5) * 2 / N, par = s * c - w, perp = s * Math.sqrt(1 - c * c), Ep = par * par + perp * perp; if (Ep > 0) sum += fn(Ep, par / Math.sqrt(Ep)) * Math.sqrt(Eo / Ep); }
    tot[j] += W * sum * (2 / N);
  }
}
// short description of the selected reactions for titles
function secSelLabel(sel) {
  const mts = sel.filter(v => v !== 'sum').map(Number).sort((a, b) => a - b), parts = [];
  if (sel.includes('sum')) parts.push('all channels');
  if (mts.length && mts.length <= 2) mts.forEach(m => parts.push('MT ' + m + ' ' + mtName(m, PJ())));
  else if (mts.length) parts.push('MT ' + (mts.length <= 6 ? mts.join(', ') : mts.slice(0, 6).join(', ') + ' …'));
  return parts.join(' + ');
}
// channels of a row that belong to the user selection (the sum = all primary channels)
function secMts(row, sel) { const set = new Set(); if (sel.includes('sum')) for (const m of row.d.order) if (row.d.channels[m].primary) set.add(m); for (const v of sel) if (v !== 'sum' && row.d.channels[+v]) set.add(+v); return Array.from(set); }
function secProducts(F, li, sp, sel) {
  const out = []; for (const row of F.rows) { if (row.t.layer !== li) continue; for (const mt of secMts(row, sel)) { const ch = row.d.channels[mt]; for (const pr of (ch.prods || [])) if (sp === 'r' ? (!pr.sp && pr.zap > 2004) : pr.sp === sp) out.push({ row, ch, pr }); } }
  return out;
}
// description of the angular data actually used for a species (shown under each plot)
function secKinds(F, li, sp, sel) {
  const k = new Set(); for (const { pr } of secProducts(F, li, sp, sel)) {
    const d = pr.dist; if (!d) { k.add('no usable distribution'); continue; } const fr = d.frame === 'lab' ? 'lab' : 'CM';
    if (d.law === 1 && d.lang === 2) k.add('Kalbach-Mann, ' + fr + (d.L.some(L => L.p[1]) ? ' (a from the file)' : ' (a from the Kalbach systematics: the file gives only r)'));
    else if (d.law === 1) k.add(d.L.some(L => L.p.length) ? 'Legendre, ' + fr : 'isotropic, ' + fr + ' (no angular coefficients in the file)');
    else if (d.law === 2) k.add('two-body discrete level, CM, ' + (d.L.some(L => L.lang === 0 && !L.a.length) ? 'isotropic' : 'angular distribution from the file'));
    else if (d.law === 6) k.add('N-body phase space, isotropic in the CM');
  } return Array.from(k).join('; ');
}
// upper edge of the secondary energy axis (MeV)
function secBound(F, li, sp, sel, kin) {
  let Emax = 0; for (const row of F.rows) if (row.t.layer === li) for (let k = F.nLb - 1; k >= 0; k--) if (row.Lh[k] > 0) { Emax = Math.max(Emax, F.Ec[k]); break; }
  let b = 0;
  for (const { row, ch, pr } of secProducts(F, li, sp, sel)) {
    const d = pr.dist; if (!d) continue;
    if (sp === 'g') { if (d.law === 1) for (let i = 0; i < d.L.length; i++) { const L = d.L[i]; for (const e of L.disc.E) b = Math.max(b, e); if (L.Ep.length) b = Math.max(b, L.Ep[L.Ep.length - 1]); if (d.E[i] >= Emax) break; } continue; }
    const AWR = row.d.AWR, Ap = pr.awp || 1, w = Math.sqrt(Ap * PJ().awi * Emax / Math.pow(PJ().awi + AWR, 2)), Ecm = Emax * AWR / (AWR + PJ().awi);
    if (d.law === 2) { const q = Ecm + ch.QI; if (q > 0) { const AR = AWR + PJ().awi - Ap; b = Math.max(b, Math.pow(Math.sqrt(AR / (AR + Ap) * q) + w, 2)); } }
    else if (d.law === 6) { const ps = d6Make(d, Ap, Ecm + ch.QM); if (ps) b = Math.max(b, Math.pow(Math.sqrt(ps.Em) + w, 2)); }
    else if (d.law === 1) { let m = 0; for (let i = 0; i < d.L.length; i++) { const L = d.L[i]; if (L.Ep.length) m = Math.max(m, L.Ep[L.Ep.length - 1]); if (d.E[i] >= Emax) break; } b = Math.max(b, d.frame === 'lab' && kin === 'table' ? m : Math.pow(Math.sqrt(m) + w, 2)); }
  }
  return b * 1.03;
}
// angles: lab angles in degrees; returns { tot (per MeV), ang[i] (per MeV per sr), W (reactions weighted with the yield), nUnsup }
function secCompute(F, li, sp, sel, nbo, Eo, angles, kin) {
  angles = angles || SEC_ANGLES; kin = kin || 'twobody'; const dEo = Eo / nbo, tot = new Float64Array(nbo), ang = angles.map(() => new Float64Array(nbo)), G = 8;
  const cosT = angles.map(a => Math.cos(a * Math.PI / 180)), sinT = angles.map(a => Math.abs(Math.sin(a * Math.PI / 180)));
  let Wtot = 0, nUnsup = 0;
  for (const { row, ch, pr } of secProducts(F, li, sp, sel)) {
    if (!pr.dist) { nUnsup++; continue; } const n = row.t.n * 1e-24, AWR = row.d.AWR;
    for (let kb = 0; kb < F.nLb; kb += G) {
      let W = 0, WE = 0;
      for (let k = kb; k < Math.min(F.nLb, kb + G); k++) { const Lk = row.Lh[k]; if (!(Lk > 0)) continue; const E = F.Ec[k], sg = chSigma(ch, E); if (!(sg > 0)) continue; const y = tab1Eval(pr.ytab, E); if (!(y > 0)) continue; const p = n * sg * Lk * y; W += p; WE += p * E; }
      if (!(W > 0)) continue; Wtot += W; secEmit(pr, ch, row.d, WE / W, W, tot, ang, nbo, dEo, cosT, sinT, kin);
    }
  }
  return { tot, ang, W: Wtot, nUnsup };
}

// HTML (charts + legends) of the secondary distributions of the selected channels, one chart per layer and species
let SECGEN = 0;
// one species in one layer: computes the spectra and stores the plot in SEC (the expensive part)
function secJob(F, li, sp, nm, sel, kin, s) {
  const R = RES, f = s.eU === 'keV' ? 1e3 : 1;
    const nbo = sp === 'g' ? 300 : 120, key = [ALIB.ver || 0, li, sp, sel.join(','), kin].join('|');
    let c = R._sec.get(key);
    if (!c) {
      let Eo = secBound(F, li, sp, sel, kin), res = Eo > 0 ? secCompute(F, li, sp, sel, nbo, Eo, null, kin) : null;
      if (res) {                                           // shrink the axis to the range that is actually populated (the bound is conservative)
        const mx = Math.max(...res.tot); let last = -1; for (let j = 0; j < nbo; j++) if (res.tot[j] > 1e-3 * mx) last = j;
        if (last >= 0 && last < 0.75 * nbo) { Eo = (last + 2) * Eo / nbo; res = secCompute(F, li, sp, sel, nbo, Eo, null, kin); }
      }
      c = { Eo, res }; R._sec.set(key, c);
    }
    if (!c.res || !(c.res.W > 0) || !(Math.max(...c.res.tot) > 0)) return { ok: false, nUn: secProducts(F, li, sp, sel).some(q => !q.pr.dist) ? 1 : 0 };
    const edges = new Float64Array(nbo + 1); for (let i = 0; i <= nbo; i++) edges[i] = i * c.Eo / nbo * f;
    const series = sp === 'g' ? [{ name: 'any angle (isotropic, per sr)', type: 'hist', edges, y: Float64Array.from(c.res.tot, v => v / f / (4 * Math.PI)), color: scol(0) }] : SEC_ANGLES.map((a, i) => ({ name: 'θ = ' + a + '°', type: 'hist', edges, y: Float64Array.from(c.res.ang[i], v => v / f), color: scol(i) }));
    series.push({ name: 'all angles (integrated, ÷4π)', type: 'hist', edges, y: Float64Array.from(c.res.tot, v => v / f / (4 * Math.PI)), color: scol(SEC_ANGLES.length) });
    const pk = li + '|' + sp; SEC.set(pk, { series, hint: fmt(c.res.W, 4) + ' ' + nm + ' per primary from the selected channels · kinematics: ' + (sp === 'g' ? 'photons are emitted isotropically (no Doppler shift); discrete lines and continuum of MF 6' : kin === 'twobody' ? 'two-body (CM energy-angle distribution boosted to the lab)' : 'lab double-differential as filed') + ' · angular data: ' + secKinds(F, li, sp, sel) + ' · ' + (sp === 'g' ? 'per sr: d²N/(dE dΩ) = dN/dE / 4π' : 'θ curves: d²N/(dE dΩ) [per ' + s.eU + ' per sr]; last curve: the angle-integrated spectrum dN/dE divided by 4π, i.e. the angle-averaged value per sr, directly comparable with the θ curves (multiply by 4π for dN/dE)') + '.',
      o: { xlabel: ({ r: 'Recoil nucleus', g: 'Photon', n: 'Neutron', p: 'Proton', d: 'Deuteron', t: 'Triton', h: '³He', a: 'α particle' })[sp] + ' energy (' + s.eU + ')', ylabel: 'per primary, per ' + s.eU + ' and per sr', title: 'Layer ' + (li + 1) + ' · ' + R.names[li] + ' · ' + nm + ' · ' + secSelLabel(sel), H: 380, xlog: s.xlog === 'log' } });
    if (!SECST.has(pk)) SECST.set(pk, { ylog: s.ylog === 'log', hidden: {}, xmin: '', xmax: '', ymin: '', ymax: '' });
    return { ok: true, nUn: c.res.nUnsup };
}
// the plots are computed one after the other, with the page staying responsive; a new render cancels the previous run
function secondaryHtml(F, sel, layers, s) {
  const kin = s.kin === 'table' ? 'table' : 'twobody', jobs = []; RES._sec = RES._sec || new Map(); SEC.clear(); const gen = ++SECGEN;
  for (const li of layers) for (const [sp, nm] of SPECIES.concat([['r', 'recoil nuclei']])) if (secProducts(F, li, sp, sel).length) jobs.push({ li, sp, nm, pk: li + '|' + sp });
  const none = nU => '<p class="hint">No secondary-particle spectra to show: the selected channels produce nothing at this beam energy (e.g. below their threshold), or the library files have no energy-angle distribution (MF 6, LAW 1, 2 or 6) for their products' + (nU ? ' (' + nU + ' product(s) have unsupported distribution laws)' : '') + '.</p>';
  if (!jobs.length) return none(0);
  let k = 0, any = 0, nUn = 0;
  const step = () => {
    const blk = document.getElementById('secBlock'); if (gen !== SECGEN || !blk) return;
    const prog = document.getElementById('secProg');
    if (k >= jobs.length) {
      if (prog) prog.remove(); if (!any) { blk.innerHTML = none(nUn); return; }
      const u = document.getElementById('secUn'); if (u && nUn) u.textContent = ' Products with unsupported laws (' + nUn + ') are not shown.'; return;
    }
    const j = jobs[k++]; if (prog) prog.textContent = 'Computing the secondary spectra… ' + (k - 1) + ' of ' + jobs.length + ' done (now: layer ' + (j.li + 1) + ', ' + j.nm + ')';
    setTimeout(() => {
      if (gen !== SECGEN) return; const el = blk.querySelector('.secplot[data-sk="' + j.pk + '"]'), r = secJob(F, j.li, j.sp, j.nm, sel, kin, s); nUn += r.nUn;
      if (r.ok) { any++; if (el) secDraw(el); } else if (el) el.remove();
      const n = document.getElementById('secN'); if (n) n.textContent = any; step();
    }, 15);
  };
  setTimeout(step, 0);
  return '<div id="secBlock"><h3 class="sech">Secondary particles: energy distributions at fixed laboratory angles</h3><div class="btns" id="secAll"><span class="lbl">Save all spectra (<span id="secN">0</span> plots, as shown)</span>' + SEC_SAVE_BTNS + '<span class="hint" id="secMsg"></span></div><p class="hint" id="secProg">Computing the secondary spectra…</p>' +
    jobs.map(j => '<div class="secplot" data-sk="' + j.pk + '"><p class="hint">Layer ' + (j.li + 1) + ', ' + j.nm + ': waiting…</p></div>').join('') +
    '<p class="hint">Built from the MF 6 data of the library (LAW 1 Legendre / Kalbach-Mann, LAW 2 two-body, LAW 6 phase space). Data in the centre-of-mass frame are transformed to the laboratory frame. Photons (discrete lines and continuum of MF 6, e.g. for (' + PJ().sym + ',γ)) are taken as isotropic, without Doppler shift; for them the plot shows the spectrum per sr and the integrated one. Between the tabulated incident energies the distributions are interpolated linearly. Drag a box on a plot to zoom (double-click resets); click a legend entry to hide or show that angle.<span id="secUn"></span></p></div>';
}

// draw (or redraw) one secondary plot with its controls: zoom by dragging, legend toggles, linear / log y axis
function secDraw(el) {
  const key = el.dataset.sk, d = SEC.get(key), st = SECST.get(key); if (!d || !st) return;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }), r = renderChart(o), g = r.geom;
  const zoomed = st.xmin !== '' || st.ymin !== '';
  el.innerHTML = '<div class="btns"><label class="lbl" style="display:inline-flex;gap:6px;align-items:center">Y axis <select data-sy style="width:auto;min-width:140px;height:30px;padding-right:26px"><option value="lin"' + (st.ylog ? '' : ' selected') + '>Linear</option><option value="log"' + (st.ylog ? ' selected' : '') + '>Logarithmic</option></select></label>' +
    '<button class="btn sm" data-sreset' + (zoomed ? '' : ' disabled') + '>Reset zoom</button><span class="hint">Drag a box to zoom · double-click resets</span></div>' +
    '<div class="chart">' + r.svg + '</div><div class="legend">' + series.map(se => '<button data-sleg="' + esc(se.name) + '" class="' + (se.hidden ? 'off' : '') + '"><i style="background:' + se.color + '"></i>' + esc(se.name) + '</button>').join('') + '</div>' +
    '<div class="readout" data-sro>Move the pointer over the plot to read values.</div><p class="hint">' + esc(d.hint) + '</p>' +
    '<div class="btns"><span class="lbl">Save this plot</span>' + SEC_SAVE_BTNS + '<span class="hint" id="smsg_' + key.replace(/\W/g, '_') + '"></span></div>';
  const again = () => secDraw(el), svg = $('svg', el);
  $('[data-sy]', el).onchange = ev => { st.ylog = ev.target.value === 'log'; st.ymin = st.ymax = ''; again(); };
  $('[data-sreset]', el).onclick = () => { st.xmin = st.xmax = st.ymin = st.ymax = ''; again(); };
  $$('[data-sleg]', el).forEach(b => b.onclick = () => { const n = b.dataset.sleg; st.hidden[n] = !st.hidden[n]; again(); });
  const pt = ev => { const q = svg.getBoundingClientRect(); return [(ev.clientX - q.left) / q.width * g.W, (ev.clientY - q.top) / q.height * g.H]; }, clamp = ([px, py]) => [Math.min(Math.max(px, g.m.l), g.m.l + g.iw), Math.min(Math.max(py, g.m.t), g.m.t + g.ih)];
  let a0 = null, box = null; const ro = $('[data-sro]', el);
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch') return; const [px, py] = pt(ev); if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
    a0 = [px, py]; svg.setPointerCapture(ev.pointerId); box = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); box.setAttribute('fill', 'var(--accent)'); box.setAttribute('fill-opacity', '.15'); box.setAttribute('stroke', 'var(--accent)'); svg.appendChild(box);
  });
  svg.addEventListener('pointermove', ev => {
    const p = clamp(pt(ev));
    if (a0) { const dx = Math.abs(p[0] - a0[0]), dy = Math.abs(p[1] - a0[1]); box.setAttribute('x', dx < 6 ? g.m.l : Math.min(p[0], a0[0])); box.setAttribute('width', dx < 6 ? g.iw : dx); box.setAttribute('y', dy < 6 ? g.m.t : Math.min(p[1], a0[1])); box.setAttribute('height', dy < 6 ? g.ih : dy); }
    const [X, Y] = pxData(g, p[0], p[1]); ro.textContent = 'x = ' + fmt(X, 5) + '   |   y = ' + fmt(Y, 4);
  });
  svg.addEventListener('pointerup', ev => {
    if (!a0) return; const p = clamp(pt(ev)), s0 = a0, dx = Math.abs(p[0] - s0[0]), dy = Math.abs(p[1] - s0[1]); a0 = null; if (box) { box.remove(); box = null; } if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(g, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])), [xb, yb] = pxData(g, Math.max(p[0], s0[0]), Math.min(p[1], s0[1]));
    if (dx >= 6) { st.xmin = fmtAx(xa); st.xmax = fmtAx(xb); } if (dy >= 6) { st.ymin = fmtAx(ya); st.ymax = fmtAx(yb); } again();
  });
  svg.addEventListener('pointercancel', () => { a0 = null; if (box) { box.remove(); box = null; } });
  svg.addEventListener('dblclick', () => { st.xmin = st.xmax = st.ymin = st.ymax = ''; again(); });
}
function secHydrate() { $$('.secplot[data-sk]').forEach(secDraw); }
