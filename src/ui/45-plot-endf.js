/* ====== tab "Reactions (ENDF)": reaction probability and secondary yields from the ion flux and the library cross sections ====== */
PLOTS.endf = {
  name: 'Reactions (ENDF)',
  defaults: { kin: 'twobody', chs: ['sum'], what: 'prob', lay: 'all', bins: 100, eU: 'MeV', ylog: 'lin', xlog: 'lin' },
  ctrl: s => {
    const chs = alibChannelList();
    const opts = [['sum', 'All reaction channels (sum)']].concat(chs.map(c => [String(c.mt), 'MT ' + c.mt + ' · ' + c.name]));
    const lay = [['all', 'All layers']].concat(S.layers.map((l, i) => [String(i), (i + 1) + ' · ' + layerName(l)]));
    return [{ k: 'chs', label: 'Reaction channels (tick one or more)', type: 'multi', opts }, { k: 'what', label: 'Plot', type: 'select', opts: [['prob', 'Reaction probability dP/dE'], ['xs', 'Cross section σ(E)']] },
      { k: 'lay', label: 'Layer', type: 'select', opts: lay }, { k: 'eU', label: 'Energy unit', type: 'select', opts: [['MeV', 'MeV'], ['keV', 'keV']] }, { k: 'kin', label: 'Secondary energies from', type: 'select', opts: [['twobody', 'two-body kinematics (CM distribution → lab)'], ['table', 'lab double-differential as filed']] }, { k: 'bins', label: 'Bins', type: 'number', min: 10, max: 1000 }];
  },
  build(s) {
    const R = RES, noData = m => ({ series: [], xlabel: '', ylabel: '', title: '', note: '<p class="hint">' + m + '</p>' });
    const P = PJ(); if (!P) return noData('This tab needs a proton, deuteron, triton, ³He or α particle as the projectile.');
    if (!alibF().size) return noData('Load the cross-section library for ' + P.plural + ' first (left panel “Nuclear reaction data (ENDF-6)”: choose the TENDL2025 folder or ' + P.lib + ', files ' + P.letter + '-Xx000.tendl).');
    if (R.pj !== P.letter) return noData('The results on screen were computed with another projectile — press Run for the ' + P.name + ' beam.');
    if (!ALIB.data.size || ![...ALIB.data.values()].some(Boolean)) return noData('No cross-section file was found for the isotopes of the target. ' + (ALIB.report && ALIB.report.miss.length ? 'Missing: ' + esc(ALIB.report.miss.join(', ')) + '.' : ''));
    if (!R.Lh) return noData('This run has no flux tally — run the simulation again.');
    const avail = new Set(['sum'].concat(alibChannelList().map(c => String(c.mt))));
    let sel = (Array.isArray(s.chs) ? s.chs : (s.ch ? [s.ch] : ['sum'])).map(String).filter(v => avail.has(v));
    if (!sel.length) return noData('Tick at least one reaction channel in the menu above.');
    const f = s.eU === 'keV' ? 1e3 : 1, F = alibFold(R), iso = F.rows, nch = sel.length, series = [];
    const chName = v => { if (v === 'sum') return 'all channels'; const c = iso.map(r => r.d.channels[+v]).find(Boolean); return 'MT ' + v + (c ? ' ' + c.name : ''); };
    const result = (row, v) => v === 'sum' ? alibTotal(F, row) : alibChannel(F, row, +v);
    const title = nch === 1 ? chName(sel[0]) : nch + ' channels';
    if (s.what === 'xs') {
      let k = 0;
      for (const v of sel) {
        const seen = new Set();
        for (const row of iso) {
          const key = row.t.Z + '-' + row.t.A; if (seen.has(key)) continue; seen.add(key);
          const ch = row.d.channels[v === 'sum' ? 3 : +v]; if (!ch || !(Math.max(...ch.s) > 0)) continue;     // no empty curves
          series.push({ name: row.t.sym + '-' + row.t.A + (nch > 1 || v === 'sum' ? ' ' + (v === 'sum' ? 'non-elastic' : 'MT ' + v) : ''), type: 'line', x: Array.from(ch.E, e => e * f), y: Array.from(ch.s), color: scol(k++), w: 1.5 });
        }
      }
      return { series, xlabel: P.sym + ' energy (' + s.eU + ')', ylabel: 'Cross section (b)', title: 'Cross section · ' + (nch === 1 && sel[0] === 'sum' ? 'non-elastic (MT 3)' : title), note: '<p class="hint">Evaluated cross sections read from the library files (TENDL, ENDF-6 MF 3).</p>' };
    }
    const nb0 = F.nLb, nb = Math.max(2, Math.min(nb0, Math.round(s.bins))), g = nb0 / nb, edges = new Float64Array(nb + 1); for (let i = 0; i <= nb; i++) edges[i] = i * F.w * g * f;
    const layers = s.lay === 'all' ? S.layers.map((l, i) => i) : [+s.lay], tab = [], skipped = []; let ci = 0;
    layers.forEach(li => {
      const resT = new Map(), rowsHtml = [], totRows = [];
      for (const v of sel) {
        const y = new Float64Array(nb); let Ptot = 0; const Y = {}; for (const [sp] of SPECIES) Y[sp] = 0;
        for (const row of iso) {
          if (row.t.layer !== li) continue; const r = result(row, v); if (!r) continue;
          for (let k = 0; k < nb0; k++) y[Math.min(nb - 1, Math.floor(k / g))] += r.dPdE[k] / (g * f);
          Ptot += r.P; for (const [sp] of SPECIES) Y[sp] += r.Y[sp]; for (const [za, q] of r.res) resT.set(za, (resT.get(za) || 0) + q);
          rowsHtml.push('<tr><td>' + (li + 1) + ' · ' + esc(R.names[li]) + '</td><td>' + esc(chName(v)) + '</td><td>' + row.t.sym + '-' + row.t.A + '</td><td>' + fmt(100 * row.t.ab, 4) + ' %</td><td class="n">' + fmt(r.P, 4) + '</td>' + SPECIES.map(([sp]) => '<td class="n">' + (r.Y[sp] ? fmt(r.Y[sp], 4) : '–') + '</td>').join('') + '</tr>');
        }
        totRows.push('<tr class="tot"><td colspan="4"><b>Layer ' + (li + 1) + ' · ' + esc(chName(v)) + ' (all isotopes)</b></td><td class="n"><b>' + fmt(Ptot, 4) + '</b></td>' + SPECIES.map(([sp]) => '<td class="n"><b>' + (Y[sp] ? fmt(Y[sp], 4) : '–') + '</b></td>').join('') + '</tr>');
        if (!(Ptot > 0)) skipped.push('layer ' + (li + 1) + ', ' + chName(v)); else series.push({ name: (nch > 1 ? chName(v) : '') + (nch > 1 && layers.length > 1 ? ' · ' : '') + (layers.length > 1 || nch === 1 ? (li + 1) + ' · ' + R.names[li] : ''), type: 'hist', edges, y, color: scol(ci++) });
      }
      const nonel = iso.filter(r => r.t.layer === li).reduce((a, r) => a + (alibChannel(F, r, 3) || { P: 0 }).P, 0), sumAll = iso.filter(r => r.t.layer === li).reduce((a, r) => a + alibTotal(F, r).P, 0), ana = R.lay.nReact[li] / R.N;
      tab.push({ li, rowsHtml, totRows, resT, nonel, sumAll, ana });
    });
    const at = endfAllTable(F, layers, R); let html = secondaryHtml(F, sel, layers, s) + '<div class="btns"><button class="btn sm" data-tcsv="all">Save all tables (CSV)</button><span class="hint" id="tm_all"></span></div>' + endfSelTable(F, layers, R, sel, chName) + at.html;
    const resAll = at.res;
    if (skipped.length) html += '<p class="hint">No reactions at this beam energy (zero probability, e.g. below the threshold) for: ' + esc(skipped.join('; ')) + ' — not plotted.</p>';
    if (nch > 1 && sel.includes('sum')) html += '<p class="warn">“All channels” overlaps with the individual channels you ticked — do not add them together.</p>';
    const chk = tab.map(t => 'layer ' + (t.li + 1) + ': Σ channels ' + fmt(t.sumAll, 4) + ', MT 3 ' + fmt(t.nonel, 4) + (ALIB.useXs && t.ana >= 0 ? ', transport ' + fmt(t.ana, 4) : '')).join(' · ');
    const off = at.chk.filter(c => c.sum > 0 && c.mt3 > 0 && Math.abs(c.sum / c.mt3 - 1) > 0.03), bad = off.filter(c => c.mt3 >= 1e-9), tiny = off.filter(c => c.mt3 < 1e-9);
    if (bad.length) html += '<p class="warn">The sum of the channels differs from MT 3 by more than 3 % (' + bad.map(c => 'layer ' + (c.li + 1) + ': ' + fmt(100 * (c.sum / c.mt3 - 1), 3) + ' %').join(', ') + '): some channels of this file may overlap or be missing, so yields could be over- or under-counted.</p>';
    if (tiny.length) html += '<p class="hint">' + tiny.map(c => 'Layer ' + (c.li + 1) + ': channel sum ' + fmt(100 * (c.sum / c.mt3 - 1), 3) + ' % from MT 3, but the reaction probability is only ' + fmt(c.mt3, 2) + ' per primary (below 10<sup>−9</sup>; cross sections far below the Coulomb barrier, where evaluations are numerically imprecise) — irrelevant for the yields.').join(' ') + '</p>';
    html += '<p class="hint">Reaction probability = Σ<sub>isotopes</sub> n<sub>i</sub> ∫ σ<sub>i</sub>(E) Φ(E) dE, with Φ(E) the path length of the transported ' + P.plural + ' per energy bin (the ion flux of the “Ion flux” tab times the layer thickness) and n<sub>i</sub> the atom density of the isotope. Yields are the same sum weighted with the multiplicity of each particle (MF 6 product yields, or the fixed multiplicities of the reaction). ' +
      'Photon yields are only counted where the file gives them (MF 6); capture photons are not. “All channels” adds the non-overlapping channels MT 5, 16, 17, 22, 28, 51-91, 102 and 600-849 (MT 3, 4 and 103-107 are sums of these). Check — ' + esc(chk) + '. Secondaries are tallied, not transported.</p>';
    return { series, xlabel: P.sym + ' energy (' + s.eU + ')', ylabel: 'dP/dE (reactions per primary per ' + s.eU + ')', title: 'Reaction probability · ' + title, note: html };
  },
};
// table of everything produced per primary, for all reaction channels (independent of the channels selected for the plots)
function endfAllTable(F, layers, R) {
  const keys = ['n', 'p', 'd', 't', 'h', 'a', 'g'], heads = ['reactions', 'neutrons', 'protons', 'deuterons', 'tritons', 'He-3', 'α', 'γ', 'rest'], res = new Map(); let rows = '';
  const cells = (P, Y, rest, b) => [P].concat(keys.map(k => Y[k]), [rest]).map(v => '<td class="n" data-v="' + (v || 0) + '">' + (b ? '<b>' : '') + (v ? fmt(v, 4) : '–') + (b ? '</b>' : '') + '</td>').join('');
  for (const li of layers) {
    let P = 0, rest = 0; const Y = {}; keys.forEach(k => Y[k] = 0);
    for (const row of F.rows) {
      if (row.t.layer !== li) continue; const r = alibTotal(F, row); let rr = 0; for (const [za, q] of r.res) { rr += q; res.set(za, (res.get(za) || 0) + q); }
      P += r.P; rest += rr; keys.forEach(k => Y[k] += r.Y[k] || 0);
      rows += '<tr><td>' + (li + 1) + ' · ' + esc(R.names[li]) + '</td><td>' + row.t.sym + '-' + row.t.A + '</td><td data-v="' + 100 * row.t.ab + '">' + fmt(100 * row.t.ab, 4) + ' %</td>' + cells(r.P, r.Y, rr, false) + '</tr>';
    }
    rows += '<tr class="tot"><td colspan="3"><b>Layer ' + (li + 1) + ' total</b></td>' + cells(P, Y, rest, true) + '</tr>';
  }
  const chk = layers.map(li => ({ li, sum: F.rows.filter(r => r.t.layer === li).reduce((a, r) => a + alibTotal(F, r).P, 0), mt3: F.rows.filter(r => r.t.layer === li).reduce((a, r) => a + (alibChannel(F, r, 3) || { P: 0 }).P, 0) }));
  // per-isotope diagnostics of the channel sum
  let dg = '';
  for (const row of F.rows) {
    if (!layers.includes(row.t.layer)) continue; const c3 = row.d.channels[3]; if (!c3) continue; const prim = row.d.order.filter(m => row.d.channels[m].primary);
    let P3 = 0, Ps = 0, worst = 0, wE = 0, w3 = 0, ws = 0; const n = row.t.n * 1e-24;
    for (let k = 0; k < F.nLb; k++) { const L = row.Lh[k]; if (!(L > 0)) continue; const E = F.Ec[k], s3 = chSigma(c3, E); let ss = 0; for (const m of prim) ss += chSigma(row.d.channels[m], E); P3 += n * s3 * L; Ps += n * ss * L; const dd = n * (s3 - ss) * L; if (Math.abs(dd) > Math.abs(worst)) { worst = dd; wE = E; w3 = s3; ws = ss; } }
    dg += '<tr><td>Layer ' + (row.t.layer + 1) + ' · ' + row.t.sym + '-' + row.t.A + '</td><td class="n" data-v="' + (P3 > 0 ? Ps / P3 : '') + '">' + (P3 > 0 ? fmt(Ps / P3, 4) : '–') + '</td><td class="n">' + (worst ? fmt(wE, 3) + ' MeV: σ(MT 3) = ' + fmt(w3, 3) + ' b, Σ channels = ' + fmt(ws, 3) + ' b' : '–') + '</td><td style="white-space:normal;text-align:left">' + esc(row.d.lumpNotes.join('; ') || '–') + '</td><td style="white-space:normal;text-align:left">MT ' + prim.join(', ') + '</td></tr>';
  }
  const diag = dg ? '<details class="cw" data-tk="check"><summary><b>Check of the channel sum, per isotope</b> — Σ channels against the non-elastic cross section MT 3</summary>' + tblBtns('check') + '<div class="tscroll"><table class="t"><tr><th>Layer · isotope</th><th class="n">Σ channels / MT 3 (flux-weighted)</th><th>Largest deficit or excess</th><th>Lumped sums</th><th>Channels added</th></tr>' + dg + '</table></div></details>' : '';
  return { res, chk, html: '<details class="cw" data-tk="everything"><summary><b>Everything produced, per primary (per incident ' + PJ().sym + '), all reaction channels</b> — independent of the channels ticked</summary>' + tblBtns('everything') + '<p class="hint">Values are per primary. “rest” = residual (recoil) nuclei.</p><div class="tscroll"><table class="t lc"><tr><th>Layer</th><th>Isotope</th><th>Abundance</th>' + heads.map(h => '<th class="n">' + h + '</th>').join('') + '</tr>' + rows + '</table></div></details>' + diag };
}
// table of the selected reaction channels only (same columns), one block per channel
function endfSelTable(F, layers, R, sel, chName) {
  const keys = ['n', 'p', 'd', 't', 'h', 'a', 'g'], heads = ['reactions', 'neutrons', 'protons', 'deuterons', 'tritons', 'He-3', 'α', 'γ', 'rest']; let rows = '';
  const cells = (P, Y, rest, b) => [P].concat(keys.map(k => Y[k]), [rest]).map(v => '<td class="n" data-v="' + (v || 0) + '">' + (b ? '<b>' : '') + (v ? fmt(v, 4) : '–') + (b ? '</b>' : '') + '</td>').join('');
  for (const li of layers) for (const v of sel) {
    let P = 0, rest = 0; const Y = {}; keys.forEach(k => Y[k] = 0);
    for (const row of F.rows) {
      if (row.t.layer !== li) continue; const r = v === 'sum' ? alibTotal(F, row) : alibChannel(F, row, +v); if (!r) continue; let rr = 0; for (const [, q] of r.res) rr += q;
      P += r.P; rest += rr; keys.forEach(k => Y[k] += r.Y[k] || 0); if (!(r.P > 0 || rr > 0 || keys.some(k => r.Y[k] > 0))) continue;
      rows += '<tr><td>' + (li + 1) + ' · ' + esc(R.names[li]) + '</td><td style="white-space:normal;text-align:left">' + esc(chName(v)) + '</td><td>' + row.t.sym + '-' + row.t.A + '</td><td data-v="' + 100 * row.t.ab + '">' + fmt(100 * row.t.ab, 4) + ' %</td>' + cells(r.P, r.Y, rr, false) + '</tr>';
    }
    if (P > 0 || rest > 0) rows += '<tr class="tot"><td colspan="4"><b>Layer ' + (li + 1) + ' · ' + esc(chName(v)) + ' (all isotopes)</b></td>' + cells(P, Y, rest, true) + '</tr>';
  }
  return '<details class="cw" data-tk="selected"><summary><b>Selected reaction channels, per primary (per incident ' + PJ().sym + ')</b></summary>' + tblBtns('selected') + '<p class="hint">Values are per primary, for the channels ticked in the menu only. “rest” = residual (recoil) nuclei.</p><div class="tscroll"><table class="t lc"><tr><th>Layer</th><th>Channel</th><th>Isotope</th><th>Abundance</th>' + heads.map(h => '<th class="n">' + h + '</th>').join('') + '</tr>' + rows + '</table></div></details>';
}
TABS.splice(TABS.indexOf('react') + 1, 0, 'endf');

// CSV of the tables of this tab: the numbers come from the data-v attribute (full precision), the text from the cells
function tblBtns(k) { return '<div class="btns"><button class="btn sm" data-tcsv="' + k + '">Save table (CSV)</button><span class="hint" id="tm_' + k + '"></span></div>'; }
function tableCsv(d) {
  const q = v => /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v, out = [], t = d.querySelector('table'); if (!t) return '';
  for (const tr of t.rows) { const r = []; for (const c of tr.cells) { r.push(q(c.dataset.v !== undefined ? c.dataset.v : c.textContent.replace(/\s+/g, ' ').trim())); for (let k = 1; k < (c.colSpan || 1); k++) r.push(''); } out.push(r.join(',')); }
  return out.join('\n') + '\n';
}
document.addEventListener('click', ev => {
  const b = ev.target.closest('[data-tcsv]'); if (!b) return; const k = b.dataset.tcsv, root = $('#plotNote') || document, base = plotBase() + '_table_';
  const one = d => '# ' + d.querySelector('summary').textContent.replace(/\s+/g, ' ').trim() + '\n' + tableCsv(d);
  if (k === 'all') { const ds = Array.from(root.querySelectorAll('details.cw[data-tk]')); if (!ds.length) return; saveFile(plotBase() + '_tables.csv', ds.map(one).join('\n'), 'CSV', '#tm_all'); }
  else { const d = root.querySelector('details.cw[data-tk="' + k + '"]'); if (d) saveFile(base + k + '.csv', one(d), 'CSV', '#tm_' + k); }
});
