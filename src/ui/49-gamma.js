/* ====== gamma-ray spectrum of a time slice (needs ENDF-6 decay data: discrete gamma lines of MF 8 / MT 457) ======
   Decays of every nuclide of the chain in the slice [t1, t2] come from the exact balance  D_i = N_i(t1) − N_i(t2) + R_i·(overlap with the irradiation) + Σ_p b_p→i D_p
   (parents first). Each decay emits the lines of the nuclide with their absolute intensities; 1 keV per channel; Gaussian resolution; detector efficiency. */
const GAM_DEF = { fwhm: 1.5, fmode: 'const', eff: 1, effMode: 'const', curve: '100 0.05\n300 0.04\n662 0.025\n1332 0.013\n3000 0.006', ylog: 'log', rate: false, top: 8, emax: 0 };
let GAM = null;                                              // state of the open window
function gamSet() { const s = irrSet(); return (s.gam = Object.assign({}, GAM_DEF, s.gam || {})); }

// ---- numerics ----
function gamErfc(x) {                                        // complementary error function, relative error < 1.2e-7 (Numerical Recipes, erfc by Chebyshev fit)
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z), r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
const gamPhi = x => 0.5 * gamErfc(-x / Math.SQRT2);
function gamFwhm(E, g) { return g.fmode === 'sqrt' ? g.fwhm * Math.sqrt(Math.max(E, 0) / 1000) : g.fwhm; }       // keV; 'sqrt': the value is the FWHM at 1 MeV, ∝ √E
function gamCurve(text) {                                     // "energy_keV efficiency" per line → sorted points
  const pts = []; for (const ln of String(text || '').split(/\r?\n/)) { const a = ln.trim().split(/[\s,;]+/).map(Number); if (a.length >= 2 && a[0] > 0 && a[1] >= 0 && isFinite(a[0]) && isFinite(a[1])) pts.push([a[0], a[1]]); }
  return pts.sort((a, b) => a[0] - b[0]);
}
function gamEffFn(g) {                                        // E (keV) → efficiency (0..1+): constant, or log-log interpolation of the curve (constant outside)
  if (g.effMode !== 'curve') { const v = +g.eff; return () => (v >= 0 ? v : 0); }
  const p = gamCurve(g.curve); if (!p.length) return () => 0; if (p.length === 1) return () => p[0][1];
  return E => {
    if (E <= p[0][0]) return p[0][1]; if (E >= p[p.length - 1][0]) return p[p.length - 1][1];
    let lo = 0, hi = p.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (p[m][0] <= E) lo = m; else hi = m; }
    const [x0, y0] = p[lo], [x1, y1] = p[hi], f = (Math.log(E) - Math.log(x0)) / (Math.log(x1) - Math.log(x0));
    return y0 > 0 && y1 > 0 ? Math.exp(Math.log(y0) + f * (Math.log(y1) - Math.log(y0))) : y0 + (E - x0) / (x1 - x0) * (y1 - y0);
  };
}
function gamAdd(arr, E, counts, fwhm) {                       // adds `counts` of a line at E (keV) to the 1 keV channels, Gaussian of the given FWHM
  const n = arr.length, sg = fwhm / 2.3548;
  if (!(sg > 0.02)) { const c = Math.floor(E); if (c >= 0 && c < n) arr[c] += counts; return; }
  const c0 = Math.max(0, Math.floor(E - 5 * sg)), c1 = Math.min(n - 1, Math.ceil(E + 5 * sg)); let prev = gamPhi((c0 - E) / sg);
  for (let c = c0; c <= c1; c++) { const nx = gamPhi((c + 1 - E) / sg); arr[c] += counts * (nx - prev); prev = nx; }
}

// ---- decays in the slice of every nuclide of a solved layer ----
function gamNearest(sol, t) {                                 // index of the output time closest to t
  let lo = 0, hi = sol.t.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (sol.t[m] <= t) lo = m; else hi = m; }
  return Math.abs(sol.t[lo] - t) <= Math.abs(sol.t[hi] - t) ? lo : hi;
}
function gamDecays(sol, prod, Tirr, t1, t2) {                 // → { D: Float64Array per node, t1e, t2e }
  const k1 = t1 > 0 ? gamNearest(sol, t1) : -1, k2 = gamNearest(sol, t2), t1e = k1 >= 0 ? sol.t[k1] : 0, t2e = sol.t[k2], n = sol.ids.length, D = new Float64Array(n);
  const ov = Math.max(0, Math.min(t2e, Tirr) - Math.max(t1e, 0));
  for (const i of sol.order) {
    let d = (k1 >= 0 ? sol.num[i][k1] : 0) - sol.num[i][k2] + (prod.get(sol.ids[i]) || 0) * ov;
    for (const [p, f] of sol.par[i]) d += f * D[p];
    D[i] = d > 0 ? d : 0;
  }
  return { D, t1e, t2e };
}

// ---- the spectrum ----
function gamBuild(sol, prod, Tirr, t1, t2, g) {
  const { D, t1e, t2e } = gamDecays(sol, prod, Tirr, t1, t2), live = Math.max(t2e - t1e, 1e-12), eff = gamEffFn(g), nodes = [];
  let Emax = 0; sol.ids.forEach((id, i) => { const e = sol.info[i]; if (e && e.gam && D[i] > 0) for (const l of e.gam) if (l[0] > Emax) Emax = l[0]; });
  const nch = Math.max(100, Math.min(12000, Math.ceil(g.emax > 0 ? +g.emax : Emax + 6 * gamFwhm(Emax, g) / 2.3548 + 2))), tot = new Float64Array(nch), totE = new Float64Array(nch);
  let sumD = 0, emitted = 0, detected = 0;
  sol.ids.forEach((id, i) => {
    const e = sol.info[i]; if (!(D[i] > 0)) return; sumD += D[i]; if (!e || !e.gam) return;
    const arr = new Float64Array(nch); let em = 0, de = 0;
    for (const [E, I] of e.gam) { const c = D[i] * I, ef = Math.max(0, eff(E)); em += c; de += c * ef; gamAdd(arr, E, c * ef, gamFwhm(E, g)); gamAdd(totE, E, c, gamFwhm(E, g)); }
    for (let c = 0; c < nch; c++) tot[c] += arr[c];
    emitted += em; detected += de;
    const top = e.gam.map(l => [l[0], l[1] * D[i]]).sort((a, b) => b[1] - a[1]).slice(0, 3);
    nodes.push({ id, i, D: D[i], gper: e.gtot, em, de, arr, top, T: e.T });
  });
  nodes.sort((a, b) => b.de - a.de);
  return { nodes, tot, totE, nch, t1e, t2e, live, sumD, emitted, detected };
}

// ---- window ----
function gamNum(v) { return v === 0 ? '0' : fmt(v, 4); }
function gamOpenAt(li, t1, t2) {                              // t in seconds since the start of the irradiation
  const d = ACTP.get(li); if (!d || !d.Q) return; const w = $('#gamWin'); w.hidden = false;
  w.innerHTML = '<div class="mbox" style="width:min(1100px,100%)"><h3>γ-ray spectrum of a time slice</h3><p class="hint">Solving the chains for the slice limits…</p></div>';
  dkSolve(Object.assign({}, d.Q, { extra: [t1, t2] })).then(sol => { GAM = { li, sol, Q: d.Q, t1, t2, hidden: {}, xmin: '', xmax: '', ymin: '', ymax: '' }; gamRender(); }).catch(e => { w.innerHTML = '<div class="mbox"><h3>γ-ray spectrum of a time slice</h3><p class="warn">' + esc(e && e.message || String(e)) + '</p><div class="btns"><button class="btn" data-gclose>Close</button></div></div>'; });
}
function gamRender() {
  const w = $('#gamWin'), G = GAM; if (!G) return; const g = gamSet(), s = irrSet(), fT = ACT_TU[s.tirrU || 'h'];
  const B = gamBuild(G.sol, G.Q.prod, G.Q.Tirr, G.t1, G.t2, g); G.B = B;
  const rate = !!g.rate, div = rate ? B.live : 1, top = Math.max(1, Math.min(30, Math.round(+g.top || 8)));
  const edges = new Float64Array(B.nch + 1); for (let c = 0; c <= B.nch; c++) edges[c] = c;
  const sc = a => rate ? Float64Array.from(a, v => v / div) : a;
  const series = [{ name: 'total (all isotopes)', type: 'hist', edges, y: sc(B.tot), color: 'var(--ink)', w: 1.6 }];
  B.nodes.slice(0, top).forEach((n, k) => series.push({ name: dkNuc(n.id), type: 'hist', edges, y: sc(n.arr), color: scol(k), w: 1.2 }));
  series.forEach(se => { se.hidden = !!G.hidden[se.name]; });
  const stats = ['slice ' + actTime(B.t1e) + ' – ' + actTime(B.t2e) + ' (live time ' + actTime(B.live) + '), layer ' + (G.li + 1),
    'γ emitted: ' + fmt(B.emitted / B.live, 4) + ' γ/s  (' + fmt(B.emitted, 4) + ' in the slice)',
    'detected: ' + fmt(B.detected / B.live, 4) + ' counts/s  (' + fmt(B.detected, 4) + ' counts)  ·  mean efficiency ' + (B.emitted > 0 ? fmt(100 * B.detected / B.emitted, 3) + ' %' : '–'),
    'decays in the slice: ' + fmt(B.sumD / B.live, 4) + ' Bq (mean)  ·  FWHM ' + (g.fmode === 'sqrt' ? fmt(g.fwhm, 3) + ' keV at 1 MeV, ∝ √E' : fmt(g.fwhm, 3) + ' keV')];
  const o = { series, xlabel: 'Energy (keV, 1 keV per channel)', ylabel: rate ? 'Counts per second per keV' : 'Counts per keV in the slice', title: 'γ spectrum · layer ' + (G.li + 1) + ' · detected', stats, xlog: false, ylog: g.ylog === 'log', xmin: G.xmin, xmax: G.xmax, ymin: G.ymin, ymax: G.ymax };
  const r = renderChart(o), gg = r.geom, zoomed = G.xmin !== '' || G.ymin !== '';
  G.o = o; G.geom = gg;
  const inp = (k, ex) => '<input type="number" data-gk="' + k + '" step="any" value="' + g[k] + '" ' + (ex || '') + '>', sel = (k, o2) => '<select data-gk="' + k + '">' + o2.map(x => '<option value="' + x[0] + '"' + (String(g[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') + '</select>';
  const fld = (l, h, wide) => '<label class="fld"' + (wide ? ' style="grid-column:1/-1"' : '') + '><span>' + l + '</span>' + h + '</label>';
  const rows = B.nodes.map(n => '<tr><td><b>' + dkNuc(n.id) + '</b></td><td class="n">' + (n.T ? fmtHl(n.T) : '–') + '</td><td class="n">' + gamNum(n.D) + '</td><td class="n">' + fmt(n.gper, 4) + '</td><td class="n">' + gamNum(n.em) + '</td><td class="n">' + gamNum(n.em / B.live) + '</td><td class="n">' + gamNum(n.de) + '</td><td class="n">' + gamNum(n.de / B.live) + '</td><td style="white-space:normal;text-align:left">' + n.top.map(l => fmt(l[0], 6) + ' keV (' + gamNum(l[1]) + ')').join(', ') + '</td></tr>').join('');
  const keep = w.querySelector('.mbox'), scrollY = keep ? keep.scrollTop : 0;
  w.innerHTML = '<div class="mbox" style="width:min(1100px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:10px"><h3>γ-ray spectrum of the time slice</h3><button class="btn sm" data-gclose>Close ✕</button></div>' +
    '<div class="mgrid">' + fld('Slice start (s)', '<input type="number" data-gt="t1" step="any" value="' + +(+B.t1e).toPrecision(8) + '" min="0">') + fld('Slice end (s)', '<input type="number" data-gt="t2" step="any" value="' + +(+B.t2e).toPrecision(8) + '" min="0">') +
    fld('Energy resolution (Gaussian), FWHM in keV', inp('fwhm', 'min="0"')) + fld('FWHM model', sel('fmode', [['const', 'constant'], ['sqrt', 'FWHM given at 1 MeV, ∝ √E']])) +
    fld('Detector efficiency', sel('effMode', [['const', 'constant'], ['curve', 'curve ε(E)']])) + (g.effMode === 'curve' ? '' : fld('Constant efficiency (0 – 1)', inp('eff', 'min="0" step="any"'))) +
    fld('Y axis', sel('ylog', [['log', 'logarithmic'], ['lin', 'linear']])) + fld('Y units', sel('rate', [['false', 'counts per keV'], ['true', 'counts per second per keV']])) + fld('Isotopes in the plot', inp('top', 'min="1" max="30" step="1"')) +
    (g.effMode === 'curve' ? '<label class="fld" style="grid-column:1/-1"><span>Efficiency curve: one pair “energy in keV, efficiency” per line (log-log interpolation, constant outside the range)</span><textarea data-gk="curve" rows="5" spellcheck="false" style="font-family:var(--font-mono);width:100%">' + esc(g.curve) + '</textarea></label>' : '') + '</div>' +
    '<div class="btns"><button class="btn sm" data-greset' + (zoomed ? '' : ' disabled') + '>Reset zoom</button><span class="hint">Drag a box to zoom · double-click resets · click a legend entry to show/hide an isotope</span></div>' +
    '<div class="chart">' + r.svg + '</div><div class="legend">' + series.map(se => '<button data-gleg="' + esc(se.name) + '" class="' + (se.hidden ? 'off' : '') + '"><i style="background:' + se.color + '"></i>' + esc(se.name) + '</button>').join('') + '</div>' +
    '<div class="readout" data-gro>Move the pointer over the plot to read values.</div>' +
    '<div class="btns"><span class="lbl">Save</span><button class="btn sm" data-gsave="png">PNG</button><button class="btn sm" data-gsave="svg">SVG</button><button class="btn sm" data-gsave="txt">ASCII (.txt)</button><button class="btn sm" data-gsave="csv">Plot CSV</button><button class="btn sm" data-gsave="spec">Spectrum CSV (all isotopes)</button><button class="btn sm" data-gsave="tab">Isotope table CSV</button><span class="hint" id="gamMsg"></span></div>' +
    '<details class="cw" open><summary><b>Isotopes contributing in the slice</b> — ' + B.nodes.length + ' with γ lines · sorted by detected counts</summary><div class="tscroll" style="max-height:340px;overflow:auto"><table class="t lc"><tr><th>Nuclide</th><th class="n">Half-life</th><th class="n">Decays</th><th class="n">γ per decay</th><th class="n">γ emitted</th><th class="n">γ emitted per s</th><th class="n">Detected counts</th><th class="n">Detected per s</th><th>Strongest lines (emitted γ in the slice)</th></tr>' + (rows || '<tr><td colspan="9">No nuclide with γ lines decays in this slice.</td></tr>') + '</table></div></details>' +
    '<p class="hint">Decays of every nuclide in the slice are exact (atom balance with the ingrowth from the parents); each decay emits the discrete γ lines of the ENDF-6 decay file (MF 8 / MT 457, absolute intensity per decay). Not included: X-rays, annihilation radiation, conversion electrons, the continuum, coincidence summing, dead time and the detector response other than the Gaussian peak and the efficiency.</p></div>';
  const nm = w.querySelector('.mbox'); if (nm) nm.scrollTop = scrollY;
  const svg = $('svg', w), ro = $('[data-gro]', w);
  const pt = ev => { const q = svg.getBoundingClientRect(); return [(ev.clientX - q.left) / q.width * gg.W, (ev.clientY - q.top) / q.height * gg.H]; }, clamp = ([px, py]) => [Math.min(Math.max(px, gg.m.l), gg.m.l + gg.iw), Math.min(Math.max(py, gg.m.t), gg.m.t + gg.ih)];
  let a0 = null, box = null;
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch') return; const [px, py] = pt(ev); if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) return;
    a0 = [px, py]; svg.setPointerCapture(ev.pointerId); box = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); box.setAttribute('fill', 'var(--accent)'); box.setAttribute('fill-opacity', '.15'); box.setAttribute('stroke', 'var(--accent)'); svg.appendChild(box);
  });
  svg.addEventListener('pointermove', ev => {
    const p = clamp(pt(ev));
    if (a0) { const dx = Math.abs(p[0] - a0[0]), dy = Math.abs(p[1] - a0[1]); box.setAttribute('x', dx < 6 ? gg.m.l : Math.min(p[0], a0[0])); box.setAttribute('width', dx < 6 ? gg.iw : dx); box.setAttribute('y', dy < 6 ? gg.m.t : Math.min(p[1], a0[1])); box.setAttribute('height', dy < 6 ? gg.ih : dy); }
    const [X, Y] = pxData(gg, p[0], p[1]); const c = Math.floor(X); ro.textContent = 'E = ' + fmt(X, 6) + ' keV   |   channel ' + c + ' (' + c + '–' + (c + 1) + ' keV): ' + series.filter(se => !se.hidden && c >= 0 && c < se.y.length).map(se => se.name + ' = ' + fmt(se.y[c], 4)).join('   |   ');
  });
  svg.addEventListener('pointerup', ev => {
    if (!a0) return; const p = clamp(pt(ev)), s0 = a0, dx = Math.abs(p[0] - s0[0]), dy = Math.abs(p[1] - s0[1]); a0 = null; if (box) { box.remove(); box = null; } if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(gg, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])), [xb, yb] = pxData(gg, Math.max(p[0], s0[0]), Math.min(p[1], s0[1]));
    if (dx >= 6) { G.xmin = fmtAx(xa); G.xmax = fmtAx(xb); } if (dy >= 6) { G.ymin = fmtAx(ya); G.ymax = fmtAx(yb); } gamRender();
  });
  svg.addEventListener('pointercancel', () => { a0 = null; if (box) { box.remove(); box = null; } });
  svg.addEventListener('dblclick', () => { G.xmin = G.xmax = G.ymin = G.ymax = ''; gamRender(); });
}
function gamCtx() {                                           // chart context of the spectrum as shown (for the exports)
  const G = GAM; if (!G || !G.o) return null; const div = document.createElement('div'); div.innerHTML = renderChart(G.o).svg;
  return { o: G.o, geom: G.geom, svg: div.querySelector('svg'), base: 'cpt_gamma_layer' + (G.li + 1) + '_' + Math.round(G.B.t1e) + 's_' + Math.round(G.B.t2e) + 's', key: 'gam' };
}
function gamCsvSpectrum() {
  const B = GAM.B, ns = B.nodes.filter(n => n.de > 0); let t = '# gamma spectrum, layer ' + (GAM.li + 1) + ', slice ' + B.t1e + ' - ' + B.t2e + ' s; detected counts per 1 keV channel in the slice (efficiency and resolution applied)\nchannel_low_keV,total_detected,total_emitted_with_resolution,' + ns.map(n => dkNuc(n.id)).join(',') + '\n';
  for (let c = 0; c < B.nch; c++) t += c + ',' + B.tot[c] + ',' + B.totE[c] + (ns.length ? ',' + ns.map(n => n.arr[c]).join(',') : '') + '\n';
  return t;
}
function gamCsvTable() {
  const B = GAM.B; let t = '# gamma emission of the isotopes in the slice ' + B.t1e + ' - ' + B.t2e + ' s (live time ' + B.live + ' s), layer ' + (GAM.li + 1) + '\nnuclide,half_life_s,decays,gamma_per_decay,gamma_emitted,gamma_emitted_per_s,detected_counts,detected_per_s,strongest_lines_keV_emitted\n';
  for (const n of B.nodes) t += [dkNuc(n.id), n.T || '', n.D, n.gper, n.em, n.em / B.live, n.de, n.de / B.live, '"' + n.top.map(l => l[0] + ' (' + l[1] + ')').join('; ') + '"'].join(',') + '\n';
  return t + '# total,,' + B.sumD + ',,' + B.emitted + ',' + B.emitted / B.live + ',' + B.detected + ',' + B.detected / B.live + '\n';
}
document.addEventListener('click', async ev => {
  if (ev.target.closest('[data-gclose]') || (ev.target.id === 'gamWin')) { $('#gamWin').hidden = true; return; }
  const b = ev.target.closest('[data-gsave]'); if (b && GAM) {
    const k = b.dataset.gsave, m = $('#gamMsg'); if (m) m.textContent = 'Preparing…';
    try {
      if (k === 'spec') await saveFile(gamCtx().base + '_spectrum.csv', gamCsvSpectrum(), 'CSV', '#gamMsg');
      else if (k === 'tab') await saveFile(gamCtx().base + '_isotopes.csv', gamCsvTable(), 'CSV', '#gamMsg');
      else { const it = await secBuild(gamCtx(), k); await saveFile(it.name, it.data, k, '#gamMsg'); }
    } catch (e) { if (m) m.textContent = 'Export failed: ' + (e && e.message || e); }
    return;
  }
  if (!GAM) return;
  if (ev.target.closest('[data-greset]')) { GAM.xmin = GAM.xmax = GAM.ymin = GAM.ymax = ''; gamRender(); return; }
  const lg = ev.target.closest('[data-gleg]'); if (lg) { const n = lg.dataset.gleg; GAM.hidden[n] = !GAM.hidden[n]; gamRender(); }
});
document.addEventListener('keydown', ev => { if (ev.key === 'Escape' && $('#gamWin') && !$('#gamWin').hidden) $('#gamWin').hidden = true; });
document.addEventListener('change', ev => {
  const e = ev.target.closest && ev.target.closest('#gamWin [data-gk]'); if (e && GAM) { const g = gamSet(), k = e.dataset.gk; g[k] = k === 'rate' ? e.value === 'true' : (e.type === 'number' ? +e.value : e.value); if (k === 'ylog') { GAM.ymin = GAM.ymax = ''; } saveLocal(); gamRender(); return; }
  const t = ev.target.closest && ev.target.closest('#gamWin [data-gt]'); if (t && GAM) {
    const a = +$('#gamWin [data-gt=t1]').value, b2 = +$('#gamWin [data-gt=t2]').value; if (isFinite(a) && isFinite(b2) && b2 > a && a >= 0) gamOpenAt(GAM.li, a, b2);
  }
});
document.addEventListener('input', ev => {                                    // the efficiency curve is applied while typing (after a pause)
  const e = ev.target.closest && ev.target.closest('#gamWin textarea[data-gk]'); if (!e || !GAM) return; gamSet().curve = e.value; clearTimeout(gamRender._t); gamRender._t = setTimeout(() => { saveLocal(); const pos = e.selectionStart; gamRender(); const n = $('#gamWin textarea[data-gk]'); if (n) { n.focus(); n.setSelectionRange(pos, pos); } }, 700);
});
