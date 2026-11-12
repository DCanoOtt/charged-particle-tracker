/* ====== tab "Activation": residual nuclei produced by a constant beam current during an irradiation, their half-lives and activities ====== */
const ACT_E = 1.602176634e-19;
const ACT_TU = { s: 1, min: 60, h: 3600, d: 86400, y: 31556926 };
const ACT_IU = { pps: 0, pA: 1e-12, nA: 1e-9, uA: 1e-6, mA: 1e-3 };
let ACT_LAST = null;
function actNuc(za) { return (ELSYM[Math.floor(za / 1000)] || '?') + '-' + (za % 1000); }
function actHl(za) {                                                    // {T (s) | null, st: 'stable' | 'rad' | 'unk', mode}
  if (HL_STABLE.has(za)) return { T: null, st: 'stable', mode: '' };
  if (HL_T[za] > 0) return { T: HL_T[za], st: 'rad', mode: HL_MODE[za] || '' };
  return { T: null, st: HL_MODE[za] ? 'rad' : 'unk', mode: HL_MODE[za] || '' };
}
function fmtHl(t) {
  if (!(t > 0)) return '–'; const u = [[1e-6, 'ns', 1e-9], [1e-3, 'µs', 1e-6], [1, 'ms', 1e-3], [60, 's', 1], [3600, 'min', 60], [86400, 'h', 3600], [31556926 * 0.999, 'd', 86400], [Infinity, 'y', 31556926]];
  for (const [lim, nm, f] of u) if (t < lim) return fmt(t / f, 4) + ' ' + nm; return '';
}
function actTime(s) { return fmtHl(s) || '0'; }
const actEv = (R, lam, T) => lam > 0 ? R * (1 - Math.exp(-lam * T)) : 0;      // activity (Bq) at time T of constant production R (atoms/s)

const ACT_DEF = { I: 1, IU: 'nA', q: 0, tirr: 1, tirrU: 'h', tcool: 0, tcoolU: 'h', lay: 'all', show: 'rad', top: 8, ylog: 'log', xlog: 'lin', qty: 'act', src: 'int', fb: true, chain: true };
function irrSet() { return (S.irr = Object.assign({}, ACT_DEF, S.irr || {})); }
const actUnits = [['s', 's'], ['min', 'min'], ['h', 'hours'], ['d', 'days'], ['y', 'years']];
let ACT_CTX = null, ACT_TMR = null, ACT_GEN = 0;
const ACTP = new Map(), ACTST = new Map();          // plot data and view state (hidden curves, zoom) per layer

// the panel "Irradiation and activation" below the plots: inputs, summary, one plot per layer and the table
function actRender() {
  const card = $('#cAct'); if (!card) return; const P = PJ();
  if (!P) { card.hidden = true; return; } card.hidden = false;
  const s = irrSet(), key = P.letter + '|' + S.layers.length + '|' + S.layers.map(l => layerName(l)).join(',') + '|' + (DEC.busy ? 'b' : DEC.nuc) + '|' + DEC.prog;
  if (card.dataset.key !== key) {                                                // controls are built once, so typing is never interrupted
    card.dataset.key = key; const nI = (k, extra) => '<input type="number" data-ak="' + k + '" step="any" ' + (extra || '') + ' value="' + (s[k] === undefined ? '' : s[k]) + '">', sel = (k, o) => '<select data-ak="' + k + '">' + o.map(x => '<option value="' + x[0] + '"' + (String(s[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') + '</select>';
    const fld = (l, h) => '<label class="fld"><span>' + l + '</span>' + h + '</label>', chk = (k, t) => '<label class="chk"><input type="checkbox" data-ak="' + k + '"' + (s[k] ? ' checked' : '') + '><span>' + t + '</span></label>';
    card.innerHTML = '<h2>Irradiation and activation</h2><p class="hint">Residual nuclei produced by a constant beam of ' + P.plural + ' in all reactions of the ' + P.sym + ' library (flux and cross sections of the simulation above), their decay chains (grow-in and decay solved with the Bateman equations), half-lives, production rates, atoms and activities.</p>' +
      '<div class="row3">' + fld('Beam current (constant)', nI('I', 'min="0"')) + fld('Current unit', sel('IU', [['pps', 'particles / s'], ['pA', 'pA'], ['nA', 'nA'], ['uA', 'µA'], ['mA', 'mA']])) + fld('Charge state q of the ions (e) — blank = fully stripped (' + P.Z + ')', nI('q', 'min="1" step="1" placeholder="' + P.Z + '"')) + '</div>' +
      '<div class="row3">' + fld('Irradiation time', nI('tirr', 'min="0"')) + fld('Unit', sel('tirrU', actUnits)) + fld('Layer', sel('lay', layOpts())) + '</div>' +
      '<div class="row3">' + fld('Cooling time after the irradiation', nI('tcool', 'min="0"')) + fld('Unit', sel('tcoolU', actUnits)) + fld('Table', sel('show', [['rad', 'radioactive nuclei'], ['all', 'all nuclei, stable too']])) + '</div>' +
      '<div class="row3">' + fld('Nuclei in each plot', nI('top', 'min="1" max="30" step="1"')) + fld('Activity axis', sel('ylog', [['log', 'logarithmic'], ['lin', 'linear']])) + fld('Time axis', sel('xlog', [['lin', 'linear'], ['log', 'logarithmic']])) + '</div>' +
      '<div class="row3">' + fld('Quantity in the plots', sel('qty', [['act', 'Decays per second = activity (radioactive nuclei)'], ['gam', 'Gamma emission rate, γ/s (needs ENDF-6 decay data)'], ['num', 'Number of atoms (all isotopes produced, stable ones too)']])) + '</div>' +
      '<div class="row3">' + fld('Decay data', sel('src', [['int', 'Internal database (ICRP-107 / NUBASE2020 with isomers)'], ['file', 'Uploaded ENDF-6 decay data' + (DEC.nuc ? ' (' + DEC.nuc.toLocaleString('en') + ' nuclides)' : ' — none loaded')]])) + '<div class="fld"><span>Uploaded decay data (ENDF-6, MF 8 / MT 457)</span><div class="btns" style="margin:0"><button class="btn sm" id="decChoose">Choose folder…</button><button class="btn sm" id="decFiles">Choose files…</button></div></div>' + '<div class="fld">' + chk('chain', 'Follow the decay chains (grow-in of daughters)') + chk('fb', 'Uploaded data: internal values for nuclides missing in the files') + '</div></div>' +
      '<p class="hint" id="decMsg">' + esc(DEC.busy ? DEC.prog : (DEC.msg || 'The internal database is used. Evaluated decay sub-library files (e.g. ENDF/B-VIII.0 decay, JEFF-3.3/4.0 RDD, JENDL-5 DDF) can be loaded from a folder or as individual files.')) + '</p>' +
      '<input type="file" id="decInput" webkitdirectory directory multiple hidden><input type="file" id="decInputF" multiple hidden><div id="actRes"></div>';
  } else card.querySelectorAll('[data-ak]').forEach(el => { if (document.activeElement !== el) { if (el.type === 'checkbox') el.checked = !!s[el.dataset.ak]; else el.value = s[el.dataset.ak] === undefined ? '' : s[el.dataset.ak]; } });
  actUpdate();
}
function actUpdate() { const my = ++ACT_GEN; clearTimeout(ACT_TMR); actCompute(my).catch(e => { const b = $('#actRes'); if (b && my === ACT_GEN) b.innerHTML = '<p class="warn">Calculation failed: ' + esc(e && e.message || String(e)) + '</p>'; }); }
async function actCompute(my) {
  const box = $('#actRes'), P = PJ(); if (!box || !P) return; const R = RES, s = irrSet(), msg = m => { ACT_CTX = null; ACT_LAST = null; box.innerHTML = '<p class="hint">' + m + '</p>'; };
  if (!R) return msg('Run the simulation first.');
  if (!alibF().size) return msg('Load the cross-section library for ' + P.plural + ' first (left panel “Nuclear reaction data (ENDF-6)”: choose the TENDL2025 folder or ' + P.lib + ', files ' + P.letter + '-Xx000.tendl).');
  if (R.pj !== P.letter) return msg('The results on screen were computed with another projectile — press Run for the ' + P.name + ' beam.');
  if (ALIB.busy) return msg('Reading the cross-section files…');
  if (!ALIB.data.size || ![...ALIB.data.values()].some(Boolean)) return msg('No cross-section file was found for the isotopes of the target.');
  if (!R.Lh) return msg('This run has no flux tally — run the simulation again.');
  const I = +s.I, Tirr = (+s.tirr) * ACT_TU[s.tirrU || 'h'], Tc = Math.max(0, (+s.tcool) || 0) * ACT_TU[s.tcoolU || 'h'];
  if (!(I > 0) || !(Tirr > 0)) return msg('Enter a beam current and an irradiation time greater than zero.');
  if (s.src === 'file' && !DEC.nuc) return msg('Load ENDF-6 decay data (button “Choose folder…” or “Choose files…”) or select the internal database.');
  const q = +s.q > 0 ? +s.q : P.Z, rate = s.IU === 'pps' ? I : I * ACT_IU[s.IU] / (q * ACT_E), Ntot = rate * Tirr;       // primaries per second
  box.innerHTML = '<p class="hint">Solving the decay chains…</p>'; await new Promise(r => setTimeout(r, 20)); if (my !== ACT_GEN) return;
  const F = alibFold(R), layers = s.lay === 'all' || +s.lay >= S.layers.length ? S.layers.map((l, i) => i) : [+s.lay], src = s.src === 'file' ? 'file' : 'int', info = id => dkInfo(id, src, !!s.fb), rows = [], plots = [];
  let nFile = 0, nInt = 0;
  for (const li of layers) {
    const m = new Map();
    for (const row of F.rows) { if (row.t.layer !== li) continue; const r = alibTotal(F, row); for (const [id, y] of r.iso) { let k = id; if (k % 10 && !info(k)) k -= k % 10; m.set(k, (m.get(k) || 0) + y); } }       // an isomer without decay data is counted as ground state
    const prod = new Map(); for (const [id, y] of m) if (y > 0) prod.set(id, y * rate);
    if (!prod.size) continue;
    if (src === 'file' && s.chain) { await decPrefetch(Array.from(prod.keys()), !!s.fb, m2 => { if (my === ACT_GEN) box.innerHTML = '<p class="hint">' + esc(m2) + '</p>'; }); if (my !== ACT_GEN) return; }
    else if (src === 'file') { await decPrefetch(Array.from(prod.keys()).filter(() => false), !!s.fb); for (const id of prod.keys()) if (!DEC.map.has(id) && !DEC.miss.has(id) && DEC.idx.has(id)) await decReadId(id); }
    const sol = await dkSolve({ prod, info, Tirr, Tc, chain: !!s.chain }); if (my !== ACT_GEN) return;
    const last = sol.nI - 1, lastC = sol.nT - 1;
    sol.ids.forEach((id, i) => {
      const e = sol.info[i], lam = sol.lam[i], direct = prod.get(id) || 0, A0 = sol.act[i][last], A1 = Tc > 0 ? sol.act[i][lastC] : A0, N0 = sol.Nirr[i], N1 = sol.Ncool ? sol.Ncool[i] : N0;
      if (e && e.src === 'file') nFile++; else if (e) nInt++;
      let peak = 0; for (let k = 0; k < sol.nT; k++) if (sol.act[i][k] > peak) peak = sol.act[i][k];
      const gt = e && e.gtot || 0; rows.push({ li, id, y: m.get(id) || 0, Rt: direct, e, lam, direct: direct > 0, A0, A1, N0, N1, peak, gt });
    });
    plots.push({ li, sol, prod });
  }
  if (my !== ACT_GEN) return;
  // plots: the strongest nuclei of every layer (peak activity), time in the unit of the irradiation
  const uT = s.tirrU || 'h', fT = ACT_TU[uT], nTop = Math.max(1, Math.min(30, Math.round(+s.top || 8))), pkey = [P.letter, Tirr, Tc, s.src, s.qty].join('|'); ACTP.clear();
  for (const { li, sol, prod } of plots) {
    const hasGam = src === 'file' && sol.info.some(e => e && e.gam), gm = s.qty === 'gam' && hasGam, atoms = s.qty === 'num', gt = sol.info.map(e => e && e.gtot || 0);
    const val = atoms ? sol.num : gm ? sol.act.map((a, i) => gt[i] > 0 ? Float64Array.from(a, v => v * gt[i]) : new Float64Array(a.length)) : sol.act;
    const rk = sol.ids.map((id, i) => { let pk = 0; const v = val[i]; for (let k = 0; k < v.length; k++) if (v[k] > pk) pk = v[k]; return { id, i, pk }; }).filter(r => r.pk > 0 && (atoms || sol.lam[r.i] > 0)).sort((a, b) => b.pk - a.pk);
    const yMax = rk.length ? rk[0].pk : 0, yMin = yMax * 1e-9, series = [], X = sol.t.map(t => t / fT);
    const mk = (name, f, color, extra) => { const x = [], y = []; for (let k = 0; k < X.length; k++) { const v = f(k); if (v > 0 && (s.ylog !== 'log' || v >= yMin)) { x.push(X[k]); y.push(v); } } if (x.length) series.push(Object.assign({ name, x, y, color, w: 1.8 }, extra || {})); };
    rk.slice(0, nTop).forEach((r, k) => { const e = sol.info[r.i]; mk(dkNuc(r.id) + '  (' + (e && e.st === 'stable' ? 'stable' : 'T½ ' + fmtHl(e && e.T)) + ')', kk => val[r.i][kk], scol(k)); });
    if (rk.length > 1 && !atoms) mk(gm ? 'total γ emission rate' : 'sum of all radioactive nuclei', k => { let a = 0; for (const r of rk) a += val[r.i][k]; return a; }, 'var(--ink)', { dash: '6 3', w: 2 });
    ACTP.set(li, { series, fT, atoms, gm, hasGam, unit: atoms ? 'atoms' : gm ? 'γ/s' : 'Bq', Q: { prod, info, Tirr, Tc, chain: !!s.chain }, o: { series, xlabel: 'Time since the start of the irradiation (' + ({ s: 's', min: 'min', h: 'hours', d: 'days', y: 'years' })[uT] + ')', ylabel: atoms ? 'Number of atoms' : gm ? 'Gamma emission rate (γ/s)' : 'Activity (decays/s = Bq)', title: 'Layer ' + (li + 1) + ' · ' + R.names[li] + ' — ' + P.sym + ' beam', xlog: false, ylog: false }, hint: 'Irradiation 0 – ' + actTime(Tirr) + (Tc > 0 ? ', cooling ' + actTime(Tirr) + ' – ' + actTime(Tirr + Tc) : '') + '. Time grid (rebuilt for these times and for the half-lives of the chain): exact Bateman propagators at ' + sol.nI + ' times in the irradiation' + (Tc > 0 ? ' and ' + (sol.nT - sol.nI) + ' in the cooling' : '') + ' (6 per octave from T/' + Math.pow(2, DK_K).toExponential(1) + ' up to T, 48 uniform, and ½, 1, 2, 3, 5 half-lives of the leading nuclei); ' + sol.nn + ' nuclei in the chain network' + (sol.fast.some(Boolean) ? ', ' + sol.fast.filter(Boolean).length + ' very short-lived ones in secular equilibrium' : '') + '. Click a legend entry to show/hide it.' + (s.qty === 'gam' && !hasGam ? ' No γ-line data: the activity is shown (upload ENDF-6 decay data to plot γ/s).' : '') });
    if (!ACTST.has(li) || ACTST.get(li).key !== pkey) ACTST.set(li, { key: pkey, hidden: {}, xmin: '', xmax: '', ymin: '', ymax: '' });
    ACTST.get(li).ylog = s.ylog === 'log'; ACTST.get(li).xlog = s.xlog === 'log';
  }
  const cool = Tc > 0, order = rows.slice().sort((a, b) => (a.li - b.li) || (b.A0 - a.A0) || (b.N0 - a.N0));
  const sorted = order.filter(r => s.show === 'all' || (r.e && r.e.st === 'stable' ? false : true));
  const stTxt = r => !r.e ? 'not in table' : r.e.st === 'stable' ? 'stable' : r.e.T ? fmtHl(r.e.T) : 'T½ n/a', orig = r => r.direct ? (r.A0 > 0 || r.N0 > 0) && r.Rt > 0 ? 'direct' : 'direct' : 'ingrowth';
  const gamAvail = sorted.some(r => r.gt > 0);
  const tr = sorted.map(r => '<tr><td>' + (r.li + 1) + ' · ' + esc(R.names[r.li]) + '</td><td><b>' + dkNuc(r.id) + '</b></td><td class="n">' + stTxt(r) + '</td><td>' + esc(r.e && r.e.br ? dkLab(r.e.br) : '') + '</td><td>' + orig(r) + '</td><td class="n">' + (r.direct ? fmt(r.y, 4) : '–') + '</td><td class="n">' + (r.direct ? fmt(r.Rt, 4) : '–') + '</td><td class="n">' + fmt(r.N0, 4) + '</td><td class="n">' + (r.lam > 0 ? fmt(r.A0, 4) : '–') + '</td>' + (cool ? '<td class="n">' + (r.lam > 0 ? fmt(r.A1, 4) : '–') + '</td>' : '') + (gamAvail ? '<td class="n">' + (r.gt > 0 ? fmt(r.gt, 4) : '–') + '</td><td class="n">' + (r.gt > 0 && r.lam > 0 ? fmt(r.A0 * r.gt, 4) : '–') + '</td>' + (cool ? '<td class="n">' + (r.gt > 0 && r.lam > 0 ? fmt(r.A1 * r.gt, 4) : '–') + '</td>' : '') : '') + '<td class="n">' + (r.lam > 0 ? fmt(r.peak, 4) : '–') + '</td></tr>').join('');
  const sumA0 = rows.reduce((a, r) => a + (r.lam > 0 ? r.A0 : 0), 0), sumA1 = rows.reduce((a, r) => a + (r.lam > 0 ? r.A1 : 0), 0), nb = rows.filter(r => r.lam > 0).length, nUnk = rows.filter(r => r.e && r.e.st !== 'stable' && !(r.lam > 0)).length, nNone = rows.filter(r => !r.e).length, nGrow = rows.filter(r => !r.direct).length;
  ACT_LAST = { rows: sorted, rate, Tirr, Tc, base: 'cpt_activation_' + P.letter };
  const srcTxt = src === 'file' ? 'uploaded ENDF-6 decay data “' + esc(DEC.name) + '” (' + nFile + ' nuclides from the files' + (nInt ? ', ' + nInt + ' from the internal database' : '') + ')' : 'internal database (ICRP-107 decay data with isomers, NUBASE2020 half-lives; mendeleev 1.3.0 for the other nuclides)';
  box.innerHTML = '<p class="hint"><b>Irradiation:</b> ' + fmt(I, 4) + ' ' + (s.IU === 'pps' ? 'particles/s' : ({ pA: 'pA', nA: 'nA', uA: 'µA', mA: 'mA' })[s.IU]) + ' of ' + P.plural + (s.IU === 'pps' ? '' : ' (charge ' + q + 'e → ' + fmt(rate, 4) + ' particles/s)') + ' during ' + actTime(Tirr) + ' = ' + fmt(Ntot, 4) + ' particles' + (Tc > 0 ? ', then ' + actTime(Tc) + ' of cooling' : '') + '. <b>Decay data:</b> ' + srcTxt + '; decay chains ' + (s.chain ? 'followed (' + nGrow + ' nuclei only from ingrowth)' : 'not followed') + '. <b>Total activity of the radioactive nuclei:</b> ' + fmt(sumA0, 4) + ' Bq at the end of the irradiation' + (Tc > 0 ? ', ' + fmt(sumA1, 4) + ' Bq after the cooling' : '') + '.' + (sorted.some(r => r.gt > 0) ? ' <b>Total γ emission rate (discrete lines of the decay data):</b> ' + fmt(rows.reduce((a, r) => a + (r.lam > 0 ? r.A0 * r.gt : 0), 0), 4) + ' γ/s at the end of the irradiation' + (Tc > 0 ? ', ' + fmt(rows.reduce((a, r) => a + (r.lam > 0 ? r.A1 * r.gt : 0), 0), 4) + ' γ/s after the cooling' : '') + '.' : (src === 'file' ? ' (the uploaded decay data have no γ lines for these nuclei)' : ' (γ emission needs uploaded ENDF-6 decay data; the internal database has none)')) + '</p>' +
    (plots.length ? plots.map(p => '<div class="secplot actp" data-al="' + p.li + '"></div>').join('') : '<p class="hint">No residual nuclei are produced.</p>') + '<span class="hint" id="actMsg"></span>' +
    '<details class="cw" open><summary><b>Residual nuclei and decay products</b> — ' + nb + ' radioactive with known half-life' + (nUnk ? ', ' + nUnk + ' without half-life' : '') + (nNone ? ', ' + nNone + ' not in the decay data' : '') + ' · sorted by activity</summary><div class="btns"><button class="btn sm" data-actcsv>Save table (CSV)</button></div>' +
    '<div class="tscroll" style="max-height:560px;overflow:auto"><table class="t lc"><tr><th>Layer</th><th>Nuclide</th><th>Half-life</th><th>Decay modes (branching)</th><th>Origin</th><th class="n">Yield per primary</th><th class="n">Direct production (atoms/s)</th><th class="n">Atoms at end of irradiation</th><th class="n">Decays/s at end (Bq)</th>' + (cool ? '<th class="n">Decays/s after cooling (Bq)</th>' : '') + (gamAvail ? '<th class="n">γ per decay</th><th class="n">γ/s at end</th>' + (cool ? '<th class="n">γ/s after cooling</th>' : '') : '') + '<th class="n">Peak decays/s (Bq)</th></tr>' + (tr || '<tr><td colspan="14">No residual nuclei.</td></tr>') + '</table></div></details>' +
    '<details><summary>How this is computed</summary><p class="hint">Residual nuclei are those of the MF 6 residual-nucleus entries of <i>all</i> non-overlapping reaction channels of the library (like the “Everything produced” table), per primary, from the same flux and cross sections as the Reactions (ENDF) tab; production rate R = yield × particles per second, constant during the irradiation (the charge state only converts an electric current into particles per second). Each layer is solved on its own. With the decay constants λ and the branching fractions b of the selected decay data, N<sub>i</sub>′ = R<sub>i</sub> + Σ<sub>p</sub> b<sub>p→i</sub> λ<sub>p</sub> N<sub>p</sub> − λ<sub>i</sub> N<sub>i</sub> is solved for the whole chain network (all daughters of all products, down to stable nuclei) during the irradiation (R on) and the cooling (R off) with the exact matrix-exponential propagator; activities are A = λN. Isomers: the production of metastable states of the products is used when the cross-section file gives it (LIP of the MF 6 products or the state cross sections of MF 10); the isomers and the isomeric transitions (IT) of the decay data are followed in the chains. Not included: burn-up of the target, reactions of secondary particles, variations of the beam current, delayed-neutron emission as a particle source.</p></details>';
  for (const p of plots) { const el = $('.actp[data-al="' + p.li + '"]', box); if (el) actDraw(el); }
  ACT_CTX = null;
}
// one plot of a layer with its controls: zoom by dragging, legend toggles, save
function actDraw(el) {
  const li = +el.dataset.al, d = ACTP.get(li), st = ACTST.get(li); if (!d || !st) return;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, xlog: st.xlog, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }), r = renderChart(o), g = r.geom, zoomed = st.xmin !== '' || st.ymin !== '';
  el.innerHTML = '<div class="btns"><button class="btn sm" data-areset' + (zoomed ? '' : ' disabled') + '>Reset zoom</button>' + (d.hasGam ? '<button class="btn sm' + (st.slice ? ' on' : '') + '" data-aslice title="Drag over a time range of the plot: a window with the γ spectrum of the isotopes decaying in it opens">' + (st.slice ? 'Drag a time range on the plot…' : 'Select a time slice → γ spectrum') + '</button><button class="btn sm" data-agam>γ spectrum of the visible range</button>' : '<span class="hint">γ spectra of a time slice: upload ENDF-6 decay data</span>') + '<span class="hint">' + (st.slice ? 'Drag over the time range of interest' : 'Drag a box to zoom · double-click resets · click the legend to show/hide a nucleus') + '</span></div>' +
    '<div class="chart">' + r.svg + '</div><div class="legend">' + series.map(se => '<button data-aleg="' + esc(se.name) + '" class="' + (se.hidden ? 'off' : '') + '"><i style="background:' + se.color + '"></i>' + esc(se.name) + '</button>').join('') + '</div>' +
    '<div class="readout" data-aro>Move the pointer over the plot to read values.</div><p class="hint">' + esc(d.hint) + '</p>' +
    '<div class="btns"><span class="lbl">Save this plot</span><button class="btn sm" data-asave="png">PNG</button><button class="btn sm" data-asave="svg">SVG</button><button class="btn sm" data-asave="txt">ASCII (.txt)</button><button class="btn sm" data-asave="csv">CSV</button><button class="btn sm" data-asave="all-png">All layers (PNG, zip)</button><button class="btn sm" data-asave="all-csv">All layers (CSV)</button></div>';
  const again = () => actDraw(el), svg = $('svg', el);
  $('[data-areset]', el).onclick = () => { st.xmin = st.xmax = st.ymin = st.ymax = ''; again(); };
  const bs = $('[data-aslice]', el); if (bs) bs.onclick = () => { st.slice = !st.slice; again(); };
  const bg = $('[data-agam]', el); if (bg) bg.onclick = () => gamOpenAt(li, Math.max(0, g.x0 * d.fT), Math.min(g.x1 * d.fT, d.Q.Tirr + d.Q.Tc));
  $$('[data-aleg]', el).forEach(b => b.onclick = () => { const n = b.dataset.aleg; st.hidden[n] = !st.hidden[n]; again(); });
  const pt = ev => { const q = svg.getBoundingClientRect(); return [(ev.clientX - q.left) / q.width * g.W, (ev.clientY - q.top) / q.height * g.H]; }, clamp = ([px, py]) => [Math.min(Math.max(px, g.m.l), g.m.l + g.iw), Math.min(Math.max(py, g.m.t), g.m.t + g.ih)];
  let a0 = null, box = null; const ro = $('[data-aro]', el);
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch') return; const [px, py] = pt(ev); if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
    a0 = [px, py]; svg.setPointerCapture(ev.pointerId); box = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); box.setAttribute('fill', 'var(--accent)'); box.setAttribute('fill-opacity', '.15'); box.setAttribute('stroke', 'var(--accent)'); svg.appendChild(box);
  });
  svg.addEventListener('pointermove', ev => {
    const p = clamp(pt(ev));
    if (a0) { const dx = Math.abs(p[0] - a0[0]), dy = Math.abs(p[1] - a0[1]); box.setAttribute('x', dx < 6 ? g.m.l : Math.min(p[0], a0[0])); box.setAttribute('width', dx < 6 ? g.iw : dx); box.setAttribute('y', dy < 6 ? g.m.t : Math.min(p[1], a0[1])); box.setAttribute('height', dy < 6 ? g.ih : dy); }
    const [X, Y] = pxData(g, p[0], p[1]); ro.textContent = 't = ' + fmt(X, 5) + '   |   ' + fmt(Y, 4) + ' ' + d.unit;
  });
  svg.addEventListener('pointerup', ev => {
    if (!a0) return; const p = clamp(pt(ev)), s0 = a0, dx = Math.abs(p[0] - s0[0]), dy = Math.abs(p[1] - s0[1]); a0 = null; if (box) { box.remove(); box = null; } if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(g, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])), [xb, yb] = pxData(g, Math.max(p[0], s0[0]), Math.min(p[1], s0[1]));
    if (st.slice) { if (dx < 6) return; st.slice = false; again(); gamOpenAt(li, Math.max(0, xa * d.fT), Math.min(xb * d.fT, d.Q.Tirr + d.Q.Tc)); return; }
    if (dx >= 6) { st.xmin = fmtAx(xa); st.xmax = fmtAx(xb); } if (dy >= 6) { st.ymin = fmtAx(ya); st.ymax = fmtAx(yb); } again();
  });
  svg.addEventListener('pointercancel', () => { a0 = null; if (box) { box.remove(); box = null; } });
  svg.addEventListener('dblclick', () => { st.xmin = st.xmax = st.ymin = st.ymax = ''; again(); });
}
function actCtx(li) {                                                       // chart context of one activity plot as shown now
  const d = ACTP.get(li), st = ACTST.get(li); if (!d || !st) return null;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, xlog: st.xlog, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }), r = renderChart(o), div = document.createElement('div'); div.innerHTML = r.svg;
  return { o, geom: r.geom, svg: div.querySelector('svg'), base: 'cpt_activation_' + PJ().letter + '_layer' + (li + 1), key: 'act' + li };
}
document.addEventListener('input', ev => {
  const el = ev.target.closest && ev.target.closest('#cAct [data-ak]'); if (!el || el.type === 'checkbox') return; const k = el.dataset.ak, s = irrSet();
  s[k] = el.type === 'number' ? (el.value === '' ? '' : +el.value) : el.value; if (k === 'q' && !(s.q > 0)) s.q = 0; saveLocal();
  clearTimeout(ACT_TMR); ACT_TMR = setTimeout(actUpdate, 250);
});
document.addEventListener('change', ev => {
  const el = ev.target.closest && ev.target.closest('#cAct [data-ak]'); if (el) { if (el.type === 'checkbox') { irrSet()[el.dataset.ak] = el.checked; saveLocal(); } clearTimeout(ACT_TMR); ACT_TMR = setTimeout(actUpdate, 0); return; }
  if (ev.target.id === 'decInput' || ev.target.id === 'decInputF') { const fl = Array.from(ev.target.files); if (fl.length) dkLoadFiles(fl, fl[0].webkitRelativePath ? fl[0].webkitRelativePath.split('/')[0] : (fl.length > 1 ? fl.length + ' files' : fl[0].name)); ev.target.value = ''; }
});
document.addEventListener('click', async ev => {
  if (ev.target.id === 'decChoose') { const i = $('#decInput'); if (i) i.click(); return; }
  if (ev.target.id === 'decFiles') { const i = $('#decInputF'); if (i) i.click(); return; }
  const sv = ev.target.closest('[data-asave]');
  if (sv) {
    const el = sv.closest('.actp'), m = $('#actMsg'), kind = sv.dataset.asave; if (m) m.textContent = 'Preparing…';
    try {
      if (kind.startsWith('all-')) {
        const k = kind.slice(4), items = []; for (const li of ACTP.keys()) { const c = actCtx(li); if (c) items.push(await secBuild(c, k)); }
        if (k === 'csv') await saveFile('cpt_activation_' + PJ().letter + '_layers.csv', items.map(it => '# ===== ' + it.name + ' =====\n' + it.data).join('\n'), 'CSV', '#actMsg');
        else { const files = []; for (const it of items) files.push({ name: it.name, data: await secBytes(it.data) }); await saveFile('cpt_activation_' + PJ().letter + '_layers_' + k + '.zip', zipMulti(files), k, '#actMsg'); }
      } else { const c = actCtx(+el.dataset.al), it = await secBuild(c, kind); await saveFile(it.name, it.data, kind, '#actMsg'); }
    } catch (e) { if (m) m.textContent = 'Export failed: ' + (e && e.message || e); }
    return;
  }
  if (!ev.target.closest('[data-actcsv]') || !ACT_LAST) return; const A = ACT_LAST, cool = A.Tc > 0;
  let t = '# Charged Particle Tracker — activation; ' + fmt(A.rate, 6) + ' particles/s for ' + A.Tirr + ' s' + (cool ? ', cooling ' + A.Tc + ' s' : '') + '\nlayer,nuclide,half_life_s,decay_modes,origin,yield_per_primary,direct_production_atoms_per_s,atoms_end_irradiation,activity_end_Bq' + (cool ? ',activity_after_cooling_Bq' : '') + ',gamma_per_decay,gamma_per_s_end' + (cool ? ',gamma_per_s_after_cooling' : '') + ',peak_activity_Bq\n';
  for (const r of A.rows) t += [r.li + 1, dkNuc(r.id), !r.e ? '' : r.e.st === 'stable' ? 'stable' : (r.e.T || ''), '"' + (r.e && r.e.br ? dkLab(r.e.br) : '') + '"', r.direct ? 'direct' : 'ingrowth', r.y, r.Rt, r.N0, r.lam > 0 ? r.A0 : 0].concat(cool ? [r.lam > 0 ? r.A1 : 0] : [], [r.gt || 0, r.lam > 0 ? r.A0 * (r.gt || 0) : 0], cool ? [r.lam > 0 ? r.A1 * (r.gt || 0) : 0] : [], [r.peak]).join(',') + '\n';
  saveFile(A.base + '.csv', t, 'CSV', '#actMsg');
});
