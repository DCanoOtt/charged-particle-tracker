/* ====== card "Production cross sections": σ(E) for the production of chosen residual nuclei, from the library of the “Nuclear reaction data (ENDF-6)” panel ======
   Input: the product (e.g. Ar-37, 38Cl, Cl-38m) and the layer (all or one). Output: σ_prod(E) in mb versus the energy of the incident particle, one curve per entry.
   σ_prod = Σ over the isotopes of the selected layer(s), weighted by their atom fraction, of Σ over the primary channels of σ_channel(E) · yield of the residual · state fraction
   (the same channels, residual yields and isomer fractions as the production rates of the activation panel). */
const XS = { hidden: {}, yl: 'log', xl: 'lin', z: {}, cache: new Map(), msg: '' };
const xsItems = () => {
  const s = irrSet();
  return Array.isArray(s.xsItems) ? s.xsItems : (s.xsItems = []);
};
// “Ar-37”, “37Ar”, “Ar37”, “Cl-38m”, “Cl38m2”, “Cl-38g”: → { ZA, st (-1 = all states, 0 = ground, n = isomer n), label } or null
function xsParse(str) {
  const t = String(str || '')
    .trim()
    .replace(/\s+/g, '');
  let m = /^(\d{1,3})-?([A-Za-z]{1,2})(g|m\d?)?$/.exec(t),
    A,
    sym,
    sfx;
  if (m) {
    A = +m[1];
    sym = m[2];
    sfx = m[3];
  } else if ((m = /^([A-Za-z]{1,2})-?(\d{1,3})(g|m\d?)?$/.exec(t))) {
    sym = m[1];
    A = +m[2];
    sfx = m[3];
  } else return null;
  const Z = ELSYM.indexOf(sym[0].toUpperCase() + sym.slice(1).toLowerCase());
  if (!(Z > 0) || !(A >= Z) || A > 300) return null;
  const st = !sfx ? -1 : sfx === 'g' ? 0 : +sfx.slice(1) || 1;
  return { ZA: Z * 1000 + A, st, label: actNuc(Z * 1000 + A) + (st < 0 ? '' : st === 0 ? 'g' : 'm' + (st > 1 ? st : '')) };
}
function xsLayerLabel(lay) {
  return lay === 'all' || +lay >= S.layers.length ? 'all layers' : 'layer ' + (+lay + 1) + ' · ' + layerName(S.layers[+lay], +lay);
}
// rows (isotope × layer) of the selection with their weights (atom fractions of the selected layers)
function xsRows(lay) {
  const all = alibTargetIsotopes().filter(t => lay === 'all' || +lay >= S.layers.length || t.layer === +lay);
  let ntot = 0;
  for (const t of all) ntot += t.n;
  const rows = [],
    miss = [];
  for (const t of all) {
    const d = ALIB.data.get(alibKey(t.Z, t.A));
    if (d) rows.push({ t, d, w: ntot > 0 ? t.n / ntot : 0 });
    else if (t.n > 0) miss.push(t.sym + '-' + t.A);
  }
  return { rows, miss };
}
// energy grid (MeV) of one entry: N equidistant points, or the energy points of the library channels that make the product (within 0 … Emax)
function xsGrid(it, rows, Emax) {
  const N = xsN();
  if (!irrSet().xsLib) {
    const x = new Float64Array(N + 1);
    for (let k = 0; k <= N; k++) x[k] = (Emax * k) / N;
    return x;
  }
  const set = new Set([0, Emax]);
  for (const row of rows)
    for (const mt of row.d.order) {
      const ch = row.d.channels[mt];
      if (!ch.primary || !ch.res.some(r => r.ZA === it.ZA)) continue;
      for (const e of ch.E) if (e > 0 && e < Emax) set.add(e);
    }
  return Float64Array.from(Array.from(set).sort((p, q) => p - q));
}
// σ (mb) of one entry on the energy grid
function xsCurve(it, Emax) {
  const { rows, miss } = xsRows(it.lay),
    x = xsGrid(it, rows, Emax),
    y = new Float64Array(x.length);
  for (let k = 0; k < x.length; k++) {
    const E = x[k];
    let v = 0;
    for (const row of rows) {
      let sum = 0;
      for (const mt of row.d.order) {
        const ch = row.d.channels[mt];
        if (!ch.primary || !ch.res.length) continue;
        const s = chSigma(ch, E);
        if (!(s > 0)) continue;
        for (const r of ch.res) {
          if (r.ZA !== it.ZA) continue;
          const m = chResYield(r, E);
          if (!(m > 0)) continue;
          if (it.st < 0) sum += s * m;
          else for (const [st, f] of chResStates(ch, r, E, s)) if (st === it.st) sum += s * m * f;
        }
      }
      v += row.w * sum;
    }
    y[k] = v * 1000; // barn → mb
  }
  return { x, y, miss, nRows: rows.length };
}
function xsN() {
  const v = Math.round(+irrSet().xsN);
  return v >= 10 ? Math.min(v, 20000) : 300;
}
function xsEmax() {
  const v = +irrSet().xsEmax;
  return v > 0 ? v : Math.max(0.1, +beamE0().toPrecision(4) || 10);
}
function xsKey(it) {
  return [it.ZA, it.st, it.lay, xsEmax(), xsN(), !!irrSet().xsLib, ALIB.ver, ALIB.data.size, S.layers.map(l => layerName(l)).join(',')].join('|');
}
function xsData(it) {
  const k = xsKey(it);
  let c = XS.cache.get(k);
  if (!c) {
    c = xsCurve(it, xsEmax());
    XS.cache.set(k, c);
    if (XS.cache.size > 60) XS.cache.delete(XS.cache.keys().next().value);
  }
  return c;
}
// materials of the layers used by the curves: “layer 1 Argon gas (Ar); layer 2 Silicon (Si)”
function xsMaterials(items) {
  const use = new Set();
  items.forEach(it => {
    if (it.lay === 'all' || +it.lay >= S.layers.length) S.layers.forEach((l, i) => use.add(i));
    else use.add(+it.lay);
  });
  return Array.from(use)
    .sort((a, b) => a - b)
    .map(i => {
      const q = setupLayerSummary(S.layers[i], i);
      return (S.layers.length > 1 ? 'layer ' + (i + 1) + ' ' : '') + layerName(S.layers[i], i) + (q.symbols ? ' (' + q.symbols + ')' : '');
    })
    .join('; ');
}
const xsName = it => it.label + '  (' + xsLayerLabel(it.lay) + ')';
function xsOpts() {
  const items = xsItems(),
    ylog = XS.yl === 'log',
    xlog = XS.xl === 'log',
    z = XS.z,
    series = [];
  items.forEach((it, i) => {
    const c = xsData(it),
      x = [],
      y = [];
    for (let k = 0; k < c.x.length; k++) {
      if ((ylog && !(c.y[k] > 0)) || (xlog && !(c.x[k] > 0))) continue;
      x.push(c.x[k]);
      y.push(c.y[k]);
    }
    if (x.length) series.push({ name: xsName(it), type: 'line', x, y, color: scol(i), w: 1.8, hidden: !!XS.hidden[xsName(it)] });
  });
  const P = PJ(),
    mats = xsMaterials(items),
    base = 'Production cross sections · ' + (P ? P.lib : 'library'),
    ttl = base + ' · target: ' + mats;
  return {
    series,
    xlog,
    ylog,
    xmin: z.xmin || (xlog ? '' : '0'),
    xmax: z.xmax || '',
    ymin: z.ymin || '',
    ymax: z.ymax || '',
    xlabel: 'Energy of the incident ' + (P ? P.name : 'particle') + ' (MeV)',
    ylabel: 'Production cross section (mb)',
    title: ttl.length > 100 ? base : ttl,
    stats: ttl.length > 100 ? ['Target: ' + mats] : []
  };
}
function xsStatus() {
  const P = PJ();
  if (!alibF().size)
    return (
      'Load the cross-section library for ' +
      (P ? P.plural : 'the projectile') +
      ' first (left panel “Nuclear reaction data (ENDF-6)”: choose the TENDL2025 folder or ' +
      (P ? P.lib : 'the library') +
      ').'
    );
  if (ALIB.busy) return 'Reading the cross-section files…';
  if (!ALIB.data.size || ![...ALIB.data.values()].some(Boolean)) return 'No cross-section file was found for the isotopes of the target.';
  return '';
}
function xsRender() {
  const card = $('#cXs');
  if (!card) return;
  foldApply(card);
  const P = PJ();
  if (!P) {
    card.hidden = true;
    return;
  }
  card.hidden = false;
  const s = irrSet(),
    key = [P.letter, S.layers.length, S.layers.map(l => layerName(l)).join(',')].join('|');
  if (card.dataset.key !== key) {
    card.dataset.key = key;
    if (s.xsYl === 'log' || s.xsYl === 'lin') XS.yl = s.xsYl;
    if (s.xsXl === 'log' || s.xsXl === 'lin') XS.xl = s.xsXl;
    card.innerHTML =
      '<h2 data-fold tabindex="0" role="button" title="Click to collapse or expand">Production cross sections</h2><p class="hint">Cross section for the production of a residual nucleus as a function of the energy of the incident ' +
      P.name +
      ', computed from the library read in the panel “Nuclear reaction data (ENDF-6)” (sum of the primary channels, MF 3 cross section × MF 6 residual yield; isomeric states from the MF 6 / MF 10 data). For a layer with several target isotopes the value is per atom of the layer: the isotopes are weighted by their atom fraction. Give the product as <code>Ar-37</code>, <code>37Ar</code> (all states) or <code>Cl-38m</code>, <code>Cl-38g</code> (isomer or ground state), choose the layer and press “Add to the plot”; add as many curves as you wish.</p>' +
      '<div class="row3"><label class="fld"><span>Product</span><input type="text" id="xsProd" placeholder="e.g. Ar-37" autocomplete="off" spellcheck="false"></label>' +
      '<label class="fld"><span>Layer</span><select id="xsLay">' +
      layOpts()
        .map(o => '<option value="' + o[0] + '">' + esc(o[1]) + '</option>')
        .join('') +
      '</select></label><div class="fld"><span>&nbsp;</span><div class="btns" style="margin:0"><button class="btn sm pri" id="xsAdd">Add to the plot</button></div></div></div>' +
      '<div class="row3"><label class="fld"><span>Maximum energy (MeV) — blank = beam energy</span><input type="number" min="0" step="any" data-xk="xsEmax" value="' +
      (s.xsEmax === undefined ? '' : s.xsEmax) +
      '"></label><label class="fld"><span>Number of energy points (10 – 20000)</span><input type="number" min="10" max="20000" step="1" data-xk="xsN" value="' +
      xsN() +
      '"></label><label class="chk"><input type="checkbox" data-xk="xsLib"' +
      (s.xsLib ? ' checked' : '') +
      '><span>use the energy points of the library instead (exact, uneven spacing)</span></label></div>' +
      '<p class="hint" id="xsMsg"></p><div class="btns" id="xsChips" style="margin:4px 0"></div><div id="xsPlot"></div>';
  }
  xsDraw();
}
function xsChips() {
  const b = $('#xsChips');
  if (!b) return;
  b.innerHTML = xsItems()
    .map(
      (it, i) =>
        '<span class="chip" style="display:inline-flex;gap:6px;align-items:center;border:1px solid var(--line);border-radius:999px;padding:2px 4px 2px 10px"><i style="width:10px;height:10px;border-radius:50%;background:' +
        scol(i) +
        '"></i>' +
        esc(xsName(it)) +
        '<button class="btn sm" data-xsdel="' +
        i +
        '" title="Remove this curve">×</button></span>'
    )
    .join('');
}
function xsDraw() {
  const box = $('#xsPlot'),
    msg = $('#xsMsg');
  if (!box) return;
  xsChips();
  const st = xsStatus();
  if (msg) msg.textContent = st;
  if (st) {
    box.innerHTML = '';
    return;
  }
  const items = xsItems();
  if (!items.length) {
    box.innerHTML = '<p class="hint">No curve yet: enter a product above.</p>';
    return;
  }
  const notes = [];
  items.forEach(it => {
    const c = xsData(it);
    if (c.miss.length) notes.push(xsName(it) + ': no library file for ' + c.miss.join(', '));
    else if (!c.nRows) notes.push(xsName(it) + ': no target isotope in this layer');
    else if (!c.y.some(v => v > 0)) notes.push(xsName(it) + ': not produced in the channels of the library');
  });
  if (msg) msg.textContent = notes.join(' · ');
  const o = xsOpts();
  if (!o.series.length) {
    box.innerHTML = '<p class="hint">Nothing to plot with these axes (no positive values).</p>';
    return;
  }
  XS.o = o;
  const r = renderChart(o);
  XS.g = r.geom;
  const names = o.series.map(se => [se.name, se.color]),
    zoomed = gamZoomed(XS.z);
  box.innerHTML =
    gamRen(
      pcTop({
        y: ['data-gyl="xs"', XS.yl],
        x: ['data-xgxl="xs"', XS.xl],
        reset: ['data-gzreset="xs"', zoomed],
        hint: 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click a legend entry to show/hide a curve'
      }),
      'x'
    ) +
    '<div class="chart" id="xsChart">' +
    r.svg +
    '</div>' +
    gamRen(gamLegend(names, XS.hidden, 'xs'), 'x') +
    '<div class="readout" id="xsRo">Move the pointer over the plot to read values.</div>' +
    pcSave({ attr: 'data-xsave', formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'rootbin', 'root'], msg: '<span class="hint" id="xsSaveMsg"></span>' });
  const sv = $('#xsChart svg', box);
  gamZoom(sv, r.geom, 'z', XS, xsDraw);
  const ro = $('#xsRo', box),
    g = r.geom;
  sv.addEventListener('pointermove', ev => {
    const q = sv.getBoundingClientRect(),
      px = ((ev.clientX - q.left) / q.width) * g.W,
      py = ((ev.clientY - q.top) / q.height) * g.H;
    if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
    const E = pxData(g, px, py)[0],
      parts = [];
    o.series.forEach(se => {
      if (se.hidden) return;
      let k = 0;
      while (k < se.x.length - 1 && se.x[k + 1] < E) k++;
      const v =
        se.x.length > 1 && se.x[k + 1] > se.x[k]
          ? se.y[k] + ((se.y[k + 1] - se.y[k]) * (Math.min(Math.max(E, se.x[k]), se.x[k + 1]) - se.x[k])) / (se.x[k + 1] - se.x[k])
          : se.y[k];
      parts.push(se.name.replace(/\s+\(.*/, '') + ' ' + fmt(v, 4) + ' mb');
    });
    ro.textContent = 'E = ' + fmt(E, 5) + ' MeV   |   ' + parts.join('   |   ');
  });
}
async function xsSave(kind) {
  const m = $('#xsSaveMsg'),
    o = XS.o;
  if (!o) return;
  if (m) m.textContent = 'Preparing…';
  const r = renderChart(o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  const base = 'cpt_production_cross_sections_' + (PJ() ? PJ().letter : ''),
    ctx = { o, geom: r.geom, svg: div.querySelector('svg'), base, key: 'xs' };
  try {
    const it = await secBuild(ctx, kind);
    if (kind === 'rootbin') {
      const rn = it.name,
        rb = RootIO.file(rn, it.hs),
        inClaude = !!(window.claude && claude.use && (await claude.use('downloads')));
      if (inClaude) await saveFile(rn.replace(/\.root$/, '') + '_histograms.zip', zipStore(rn, rb), 'ROOT file', '#xsSaveMsg');
      else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', '#xsSaveMsg');
      if (m) m.textContent = 'Saved ' + rn + (inClaude ? ' (inside a zip)' : '') + ' with ' + it.hs.length + ' histograms.';
    } else await saveFile(it.name, it.data, kind, '#xsSaveMsg');
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  }
}
function xsAdd() {
  const inp = $('#xsProd'),
    lay = $('#xsLay'),
    msg = $('#xsMsg');
  if (!inp) return;
  const p = xsParse(inp.value);
  if (!p) {
    if (msg) msg.textContent = 'Product not understood: write the isotope as Ar-37, 37Ar or Cl-38m (symbol and mass number).';
    return;
  }
  const it = { ZA: p.ZA, st: p.st, label: p.label, lay: lay.value },
    items = xsItems();
  if (items.some(q => xsName(q) === xsName(it))) {
    if (msg) msg.textContent = xsName(it) + ' is already in the plot.';
    return;
  }
  items.push(it);
  saveLocal();
  inp.value = '';
  xsDraw();
}
document.addEventListener('click', ev => {
  const c = ev.target.closest && ev.target.closest('#cXs');
  if (!c) return;
  if (ev.target.closest('#xsAdd')) return xsAdd();
  const d = ev.target.closest('[data-xsdel]');
  if (d) {
    xsItems().splice(+d.dataset.xsdel, 1);
    saveLocal();
    return xsDraw();
  }
  if (ev.target.closest('[data-xgzr]')) {
    XS.z = {};
    return xsDraw();
  }
  const ga = ev.target.closest('[data-xgall]');
  if (ga) {
    xsItems().forEach(it => {
      XS.hidden[xsName(it)] = ga.dataset.gstate === 'off';
    });
    return xsDraw();
  }
  const lg = ev.target.closest('[data-xleg]');
  if (lg) {
    XS.hidden[lg.dataset.xleg] = !XS.hidden[lg.dataset.xleg];
    return xsDraw();
  }
  const sv = ev.target.closest('[data-xsave]');
  if (sv) xsSave(sv.dataset.xsave);
});
document.addEventListener('keydown', ev => {
  if (ev.key === 'Enter' && ev.target && ev.target.id === 'xsProd') xsAdd();
});
document.addEventListener('change', ev => {
  const t = ev.target;
  if (!t.closest || !t.closest('#cXs')) return;
  if (t.matches('[data-xgyl]')) {
    XS.yl = irrSet().xsYl = t.value;
    saveLocal();
    XS.z = Object.assign({}, XS.z, { ymin: '', ymax: '' });
    return xsDraw();
  }
  if (t.matches('[data-xgxl]')) {
    XS.xl = irrSet().xsXl = t.value;
    saveLocal();
    XS.z = Object.assign({}, XS.z, { xmin: '', xmax: '' });
    return xsDraw();
  }
  if (t.matches('[data-xk]')) {
    irrSet()[t.dataset.xk] = t.type === 'checkbox' ? t.checked : t.value === '' ? undefined : +t.value;
    XS.z = {};
    saveLocal();
    xsDraw();
  }
});
