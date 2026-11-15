/* ====== activation: the plots (continued from 48-activation.js) ====== */
// one plot of a layer with its controls: zoom by dragging, legend toggles, save
function actDraw(el) {
  const key = el.dataset.al,
    d = ACTP.get(isNaN(+key) ? key : +key),
    st = ACTST.get(isNaN(+key) ? key : +key);
  if (!d || !st) return;
  const li = d.li;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, xlog: st.xlog, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }),
    r = renderChart(o),
    g = r.geom,
    zoomed = st.xmin !== '' || st.ymin !== '';
  el.innerHTML =
    pcTop({
      y: ['data-ak="ylog"', irrSet().ylog],
      x: ['data-ak="xlog"', irrSet().xlog],
      reset: ['data-areset', zoomed],
      hint: st.slice
        ? 'Drag over the time range of interest'
        : 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click the legend to show/hide a nucleus',
      extra: d.hasGam
        ? '<button class="btn sm' +
          (st.slice ? ' on' : '') +
          '" data-aslice title="Drag over a time range of the plot: a window with the γ spectrum of the isotopes decaying in it opens">' +
          (st.slice ? 'Drag a time range on the plot…' : 'Select a time slice → γ spectrum') +
          '</button>'
        : '<span class="hint">γ spectra of a time slice: upload ENDF-6 decay data</span>'
    }) +
    '<div class="chart">' +
    r.svg +
    '</div><div class="legend">' +
    series
      .map(
        se =>
          '<button data-aleg="' +
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
    '<button class="btn sm" data-aall="off" title="hide all curves of this plot">Unselect all</button><button class="btn sm" data-aall="on" title="show all curves of this plot">Select all</button></div>' +
    '<div class="readout" data-aro>Move the pointer over the plot to read values.</div><p class="hint">' +
    esc(d.hint) +
    '</p>' +
    pcSave({
      attr: 'data-asave',
      extra: [
        ['all-png', 'All layers (PNG, zip)', ''],
        ['all-csv', 'All layers (CSV)', '']
      ]
    });
  const again = () => actDraw(el),
    svg = $('svg', el);
  $('[data-areset]', el).onclick = () => {
    st.xmin = st.xmax = st.ymin = st.ymax = '';
    again();
  };
  const bs = $('[data-aslice]', el);
  if (bs)
    bs.onclick = () => {
      st.slice = !st.slice;
      again();
    };
  $$('[data-aall]', el).forEach(
    b =>
      (b.onclick = () => {
        d.series.forEach(se => {
          st.hidden[se.name] = b.dataset.aall === 'off';
        });
        again();
      })
  );
  $$('[data-aleg]', el).forEach(
    b =>
      (b.onclick = () => {
        const n = b.dataset.aleg;
        st.hidden[n] = !st.hidden[n];
        again();
      })
  );
  const pt = ev => {
      const q = svg.getBoundingClientRect();
      return [((ev.clientX - q.left) / q.width) * g.W, ((ev.clientY - q.top) / q.height) * g.H];
    },
    clamp = ([px, py]) => [Math.min(Math.max(px, g.m.l), g.m.l + g.iw), Math.min(Math.max(py, g.m.t), g.m.t + g.ih)];
  let a0 = null,
    box = null;
  const ro = $('[data-aro]', el);
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch') return;
    const [px, py] = pt(ev);
    if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
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
      box.setAttribute('x', dx < 6 ? g.m.l : Math.min(p[0], a0[0]));
      box.setAttribute('width', dx < 6 ? g.iw : dx);
      box.setAttribute('y', dy < 6 ? g.m.t : Math.min(p[1], a0[1]));
      box.setAttribute('height', dy < 6 ? g.ih : dy);
    }
    const [X, Y] = pxData(g, p[0], p[1]);
    ro.textContent = 't = ' + fmt(X, 5) + '   |   ' + fmt(Y, 4) + ' ' + d.unit;
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
    const [xa, ya] = pxData(g, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])),
      [xb, yb] = pxData(g, Math.max(p[0], s0[0]), Math.min(p[1], s0[1]));
    if (st.slice) {
      if (dx < 6) return;
      st.slice = false;
      again();
      gamOpenAt(li, Math.max(0, xa * d.fT), Math.min(xb * d.fT, d.Q.Tirr + d.Q.Tc));
      return;
    }
    if (dx >= 6) {
      st.xmin = fmtAx(xa);
      st.xmax = fmtAx(xb);
    }
    if (dy >= 6) {
      st.ymin = fmtAx(ya);
      st.ymax = fmtAx(yb);
    }
    again();
  });
  svg.addEventListener('pointercancel', () => {
    a0 = null;
    if (box) {
      box.remove();
      box = null;
    }
  });
  svg.addEventListener('dblclick', () => {
    st.xmin = st.xmax = st.ymin = st.ymax = '';
    again();
  });
}
function actCtx(key) {
  // chart context of one activity plot as shown now (key: layer index, or 'g' + layer index for the γ emission plot)
  const d = ACTP.get(key),
    st = ACTST.get(key);
  if (!d || !st) return null;
  const li = d.li;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, xlog: st.xlog, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }),
    r = renderChart(o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  return {
    o,
    geom: r.geom,
    svg: div.querySelector('svg'),
    base: 'cpt_activation_' + PJ().letter + '_layer' + (li + 1) + (d.gam ? '_gamma' : ''),
    key: 'act' + key
  };
}
