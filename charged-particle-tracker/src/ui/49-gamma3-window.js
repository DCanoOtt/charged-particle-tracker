/* ====== γ spectrum: the window (continued from 49-gamma.js) ====== */
// ---- window ----
function gamNum(v) {
  return v === 0 ? '0' : fmt(v, 4);
}
function gamOpenAt(li, t1, t2) {
  // t in seconds since the start of the irradiation
  const d = ACTP.get(li);
  if (!d || !d.Q) return;
  const w = $('#gamWin');
  w.hidden = false;
  w.innerHTML =
    '<div class="mbox" style="width:min(1100px,100%)"><h3>γ-ray spectrum of a time slice</h3><p class="hint">Solving the chains for the slice limits…</p></div>';
  dkSolve(Object.assign({}, d.Q, { extra: [t1, t2] }))
    .then(sol => {
      GAM = { li, sol, Q: d.Q, t1, t2, hidden: {}, xmin: '', xmax: '', ymin: '', ymax: '' };
      gamRender();
    })
    .catch(e => {
      w.innerHTML =
        '<div class="mbox"><h3>γ-ray spectrum of a time slice</h3><p class="warn">' +
        esc((e && e.message) || String(e)) +
        '</p><div class="btns"><button class="btn" data-gclose>Done</button></div></div>';
    });
}
function gamRender() {
  const w = $('#gamWin'),
    G = GAM;
  if (!G) return;
  const g = gamSet(),
    s = irrSet(),
    fT = ACT_TU[s.tirrU || 'h'];
  const B = gamBuild(G.sol, G.Q.prod, G.Q.Tirr, G.t1, G.t2, g);
  G.B = B;
  const rate = !!g.rate,
    div = rate ? B.live : 1,
    top = Math.max(1, Math.min(30, Math.round(+g.top || 8)));
  const edges = new Float64Array(B.nch + 1);
  for (let c = 0; c <= B.nch; c++) edges[c] = c;
  const sc = a => (rate ? Float64Array.from(a, v => v / div) : a);
  const series = [{ name: 'total (all isotopes)', type: 'hist', edges, y: sc(B.tot), color: 'var(--ink)', w: 1.6 }];
  B.nodes.slice(0, top).forEach((n, k) => series.push({ name: dkNuc(n.id), type: 'hist', edges, y: sc(n.arr), color: scol(k), w: 1.2 }));
  const PU = gamPile(B, G.sol, g);
  B.pu = PU;
  if (PU.spec) {
    const e2 = new Float64Array(PU.spec.length + 1);
    for (let c = 0; c < e2.length; c++) e2[c] = c;
    series.push({ name: 'total with pile-up', type: 'hist', edges: e2, y: sc(PU.spec), color: 'var(--bad)', w: 1.4 });
    series.push({ name: 'pile-up sum pulses only', type: 'hist', edges: e2, y: sc(PU.sums), color: 'var(--warn)', w: 1.2 });
  }
  series.forEach(se => {
    se.hidden = !!G.hidden[se.name];
  });
  const stats = [
    'slice ' + actTime(B.t1e) + ' – ' + actTime(B.t2e) + ' (live time ' + actTime(B.live) + '), layer ' + (G.li + 1),
    'γ emitted: ' + fmt(B.emitted / B.live, 4) + ' γ/s  (' + fmt(B.emitted, 4) + ' in the slice)',
    'detected: ' +
      fmt(B.detected / B.live, 4) +
      ' counts/s  (' +
      fmt(B.detected, 4) +
      ' counts)  ·  mean efficiency ' +
      (B.emitted > 0 ? fmt((100 * B.detected) / B.emitted, 3) + ' %' : '–'),
    'pile-up (τ = ' +
      fmt(PU.tau * 1e6, 3) +
      ' µs' +
      (PU.f > 1 ? ', ' + fmt(PU.f, 3) + ' pulses per full-energy count' : '') +
      '): pulse rate ' +
      fmt(PU.n0, 3) +
      ' → ' +
      fmt(PU.n1, 3) +
      ' /s, mean ' +
      fmt(PU.nbar, 3) +
      ' /s',
    '   PILE-UP PROBABILITY = 1 − e^(−2μ) = ' +
      fmt(100 * PU.Pm, 3) +
      ' % (mean rate)  ·  ' +
      fmt(100 * PU.P0, 3) +
      ' % (slice start)  ·  ' +
      fmt(100 * PU.Pmax, 3) +
      ' % (highest rate)  ·  trailing only 1 − e^(−μ) = ' +
      fmt(100 * PU.Pone, 3) +
      ' %',
    '   μ = nτ = ' +
      fmt(PU.mu, 3) +
      ' (mean), ' +
      fmt(PU.mu0, 3) +
      ' (start)  ·  piled-up pairs in the slice: ' +
      fmt(PU.pairs, 3) +
      ' = ' +
      fmt((100 * PU.pairs) / Math.max(PU.nbar * B.live, 1e-300), 3) +
      ' % of the pulses' +
      (PU.mumax > 0.1 ? '  ⚠ μ > 0.1: first-order model, not reliable' : ''),
    'decays in the slice: ' +
      fmt(B.sumD / B.live, 4) +
      ' Bq (mean)  ·  FWHM ' +
      (g.fmode === 'sqrt' ? fmt(g.fwhm, 3) + ' keV at 1 MeV, ∝ √E' : fmt(g.fwhm, 3) + ' keV')
  ];
  const o = {
    series,
    xlabel: 'Energy (keV, 1 keV per channel)',
    ylabel: rate ? 'Counts per second per keV' : 'Counts per keV in the slice',
    title: 'γ spectrum · layer ' + (G.li + 1) + ' · detected',
    stats,
    xlog: false,
    ylog: g.ylog === 'log',
    xmin: G.xmin,
    xmax: G.xmax,
    ymin: G.ymin,
    ymax: G.ymax
  };
  const r = renderChart(o),
    gg = r.geom,
    zoomed = G.xmin !== '' || G.ymin !== '';
  G.o = o;
  G.geom = gg;
  G.names = { main: series.map(se => se.name) };
  const inp = (k, ex) => '<input type="number" data-gk="' + k + '" step="any" value="' + g[k] + '" ' + (ex || '') + '>',
    sel = (k, o2) =>
      '<select data-gk="' +
      k +
      '">' +
      o2.map(x => '<option value="' + x[0] + '"' + (String(g[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') +
      '</select>';
  const fld = (l, h, wide) => '<label class="fld"' + (wide ? ' style="grid-column:1/-1"' : '') + '><span>' + l + '</span>' + h + '</label>';
  const rows = B.nodes
    .map(
      n =>
        '<tr><td><b>' +
        dkNuc(n.id) +
        '</b></td><td class="n">' +
        (n.T ? fmtHl(n.T) : '–') +
        '</td><td class="n">' +
        gamNum(n.D) +
        '</td><td class="n">' +
        fmt(n.gper, 4) +
        '</td><td class="n">' +
        gamNum(n.em) +
        '</td><td class="n">' +
        gamNum(n.em / B.live) +
        '</td><td class="n">' +
        gamNum(n.de) +
        '</td><td class="n">' +
        gamNum(n.de / B.live) +
        '</td><td style="white-space:normal;text-align:left">' +
        n.top.map(l => fmt(l[0], 6) + ' keV (' + gamNum(l[1]) + ')').join(', ') +
        '</td></tr>'
    )
    .join('');
  let geBlock = '';
  if (g.effMode === 'ge') geBlock = gamGeBlock(G, B, g);
  const keep = w.querySelector('.mbox'),
    scrollY = keep ? keep.scrollTop : 0;
  w.innerHTML =
    '<div class="mbox" style="width:min(1100px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:10px"><h3>γ-ray spectrum of the time slice</h3><button class="btn sm" data-gclose>Done</button></div>' +
    '<div class="mgrid">' +
    fld('Slice start (s)', '<input type="number" data-gt="t1" step="any" value="' + +(+B.t1e).toPrecision(8) + '" min="0">') +
    fld('Slice end (s)', '<input type="number" data-gt="t2" step="any" value="' + +(+B.t2e).toPrecision(8) + '" min="0">') +
    fld(
      'Y units',
      sel('rate', [
        ['false', 'counts per keV'],
        ['true', 'counts per second per keV']
      ])
    ) +
    fld(
      'X-rays',
      sel('xray', [
        ['true', 'included as lines'],
        ['false', 'not included']
      ])
    ) +
    fld('Isotopes in the plot', inp('top', 'min="1" max="30" step="1"')) +
    '<div class="hint" style="grid-column:1/-1"><b>Ge detector and DAQ:</b> <span data-gcsum>' +
    esc(gamCfgSummary(g)) +
    '</span> <button class="btn sm" data-gecfg>Configure the Ge detector and DAQ…</button> <span>The calculations of this window use these parameters.</span></div>' +
    '</div>' +
    '<details class="cw" open><summary><b>γ-ray spectrum</b> — layer ' +
    (G.li + 1) +
    ' · detected counts · statistics in the plot</summary>' +
    pcTop({
      y: ['data-gk="ylog"', g.ylog],
      reset: ['data-greset', zoomed],
      hint: 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click a legend entry to show/hide an isotope'
    }) +
    '<div class="chart" id="gamMainChart">' +
    r.svg +
    '</div><div class="legend">' +
    series
      .map(
        se =>
          '<button data-gleg="' +
          esc(se.name) +
          '" class="' +
          (se.hidden ? 'off' : '') +
          '"><i style="background:' +
          se.color +
          '"></i>' +
          esc(se.name) +
          '</button>'
      )
      .join('') +
    '' +
    pcLegendTools('data-gall="main" data-gstate="off"', 'data-gall="main" data-gstate="on"') +
    '</div>' +
    '<div class="readout" data-gro>Move the pointer over the plot to read values.</div>' +
    '</details>' +
    geBlock +
    pcSave({
      attr: 'data-gsave',
      extra: [
        ['spec', 'Spectrum CSV (all isotopes)', ''],
        ['tab', 'Isotope table CSV', '']
      ],
      msg: '<span class="hint" id="gamMsg"></span>'
    }) +
    '<details class="cw" open><summary><b>Isotopes contributing in the slice</b> — ' +
    B.nodes.length +
    ' with γ lines · sorted by detected counts</summary><div class="tscroll" style="max-height:340px;overflow:auto"><table class="t lc"><tr><th>Nuclide</th><th class="n">Half-life</th><th class="n">Decays</th><th class="n">γ per decay</th><th class="n">γ emitted</th><th class="n">γ emitted per s</th><th class="n">Detected counts</th><th class="n">Detected per s</th><th>Strongest lines (emitted γ in the slice)</th></tr>' +
    (rows || '<tr><td colspan="9">No nuclide with γ lines decays in this slice.</td></tr>') +
    '</table></div></details>' +
    '<p class="hint">Decays of every nuclide in the slice are exact (atom balance with the ingrowth from the parents); each decay emits the discrete γ lines of the ENDF-6 decay file (MF 8 / MT 457, absolute intensity per decay). X-rays (ENSDF: computed from the electron-capture and conversion vacancies; ENDF-6: the tabulated lines) are shown as lines when “X-rays” is on. Not included: conversion electrons, the continuum, coincidence summing, dead time and the detector response other than the Gaussian peak and the efficiency.</p></div>';
  const nm = w.querySelector('.mbox');
  if (nm) nm.scrollTop = scrollY;
  const svg = $('svg', w),
    ro = $('[data-gro]', w);
  const pt = ev => {
      const q = svg.getBoundingClientRect();
      return [((ev.clientX - q.left) / q.width) * gg.W, ((ev.clientY - q.top) / q.height) * gg.H];
    },
    clamp = ([px, py]) => [Math.min(Math.max(px, gg.m.l), gg.m.l + gg.iw), Math.min(Math.max(py, gg.m.t), gg.m.t + gg.ih)];
  let a0 = null,
    box = null;
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch' || ev.shiftKey || G.pick) return;
    const [px, py] = pt(ev);
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) return;
    a0 = [px, py];
    svg.setPointerCapture(ev.pointerId);
    box = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    box.setAttribute('fill', 'var(--accent)');
    box.setAttribute('fill-opacity', '.15');
    box.setAttribute('stroke', 'var(--accent)');
    svg.appendChild(box);
  });
  svg.addEventListener('pointermove', ev => {
    const p = clamp(pt(ev));
    if (a0) {
      const dx = Math.abs(p[0] - a0[0]),
        dy = Math.abs(p[1] - a0[1]);
      box.setAttribute('x', dx < 6 ? gg.m.l : Math.min(p[0], a0[0]));
      box.setAttribute('width', dx < 6 ? gg.iw : dx);
      box.setAttribute('y', dy < 6 ? gg.m.t : Math.min(p[1], a0[1]));
      box.setAttribute('height', dy < 6 ? gg.ih : dy);
    }
    const [X, Y] = pxData(gg, p[0], p[1]);
    const c = Math.floor(X);
    ro.textContent =
      'E = ' +
      fmt(X, 6) +
      ' keV   |   channel ' +
      c +
      ' (' +
      c +
      '–' +
      (c + 1) +
      ' keV): ' +
      series
        .filter(se => !se.hidden && c >= 0 && c < se.y.length)
        .map(se => se.name + ' = ' + fmt(se.y[c], 4))
        .join('   |   ');
  });
  svg.addEventListener('pointerup', ev => {
    if (!a0) return;
    const p = clamp(pt(ev)),
      s0 = a0,
      dx = Math.abs(p[0] - s0[0]),
      dy = Math.abs(p[1] - s0[1]);
    a0 = null;
    if (box) {
      box.remove();
      box = null;
    }
    if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(gg, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])),
      [xb, yb] = pxData(gg, Math.max(p[0], s0[0]), Math.min(p[1], s0[1]));
    if (dx >= 6) {
      G.xmin = fmtAx(xa);
      G.xmax = fmtAx(xb);
    }
    if (dy >= 6) {
      G.ymin = fmtAx(ya);
      G.ymax = fmtAx(yb);
    }
    gamRender();
  });
  svg.addEventListener('pointercancel', () => {
    a0 = null;
    if (box) {
      box.remove();
      box = null;
    }
  });
  svg.addEventListener('dblclick', () => {
    G.xmin = G.xmax = G.ymin = G.ymax = '';
    gamRender();
  });
  if (G.geoms) {
    gamZoom($('#geEffChart svg', w), G.geoms.eff, 'ze');
    gamZoom($('#geSimChart svg', w), G.geoms.sim, 'zs');
  }
  specTip(svg, gg, () => G.o);
  if (G.geoms) specTip($('#geSimChart svg', w), G.geoms.sim, () => G.subO && G.subO.sim);
}
