/* ====== zoom ====== */
function fmtAx(v) { return String(+v.toPrecision(6)); }
function pxData(g, px, py) {
  const fx = (px - g.m.l) / g.iw, fy = 1 - (py - g.m.t) / g.ih;
  const X = g.xlog ? Math.pow(10, Math.log10(g.x0) + fx * (Math.log10(g.x1) - Math.log10(g.x0))) : g.x0 + fx * (g.x1 - g.x0);
  const Y = g.ylog ? Math.pow(10, Math.log10(g.y0) + fy * (Math.log10(g.y1) - Math.log10(g.y0))) : g.y0 + fy * (g.y1 - g.y0);
  return [X, Y];
}
function setRange(x0, x1, y0, y1) {
  const s = pset(S.tab);
  if (x0 != null) { s.xmin = fmtAx(x0); s.xmax = fmtAx(x1); }
  if (y0 != null) { s.ymin = fmtAx(y0); s.ymax = fmtAx(y1); }
  ['xmin', 'xmax', 'ymin', 'ymax'].forEach(k => { const i = $('[data-pk="' + k + '"]'); if (i) i.value = s[k]; });
  drawChart(); saveLocal();
}
function zoomBy(f) {   // f<1 zooms in, f>1 zooms out (both axes, about the centre)
  if (!_chart) return; const g = _chart.geom;
  const rg = (a, b, lg) => { if (lg) { const la = Math.log10(a), lb = Math.log10(b), c = (la + lb) / 2, h = (lb - la) / 2 * f; return [Math.pow(10, c - h), Math.pow(10, c + h)]; } const c = (a + b) / 2, h = (b - a) / 2 * f; return [c - h, c + h]; };
  const [a, b] = rg(g.x0, g.x1, g.xlog), [c, d] = rg(g.y0, g.y1, g.ylog); setRange(a, b, c, d);
}
function resetZoom() { const s = pset(S.tab); ['xmin', 'xmax', 'ymin', 'ymax'].forEach(k => { s[k] = ''; const i = $('[data-pk="' + k + '"]'); if (i) i.value = ''; }); drawChart(); saveLocal(); }
function attachZoom(svg) {
  let st = null, box = null;
  const svgPt = ev => { const r = svg.getBoundingClientRect(), g = _chart.geom; return [(ev.clientX - r.left) / r.width * g.W, (ev.clientY - r.top) / r.height * g.H]; };
  svg.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch' || !_chart) return;
    const g = _chart.geom, [px, py] = svgPt(ev); if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
    st = [px, py]; svg.setPointerCapture(ev.pointerId);
    box = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); box.setAttribute('fill', 'var(--accent)'); box.setAttribute('fill-opacity', '.15'); box.setAttribute('stroke', 'var(--accent)'); box.setAttribute('stroke-dasharray', '4 3'); svg.appendChild(box);
  });
  svg.addEventListener('pointermove', ev => {
    if (!st) return; const g = _chart.geom; let [px, py] = svgPt(ev);
    px = Math.min(Math.max(px, g.m.l), g.m.l + g.iw); py = Math.min(Math.max(py, g.m.t), g.m.t + g.ih);
    const dx = Math.abs(px - st[0]), dy = Math.abs(py - st[1]);
    const x = Math.min(px, st[0]), w = dx, y = dy < 6 ? g.m.t : Math.min(py, st[1]), h = dy < 6 ? g.ih : dy;
    const xx = dx < 6 ? g.m.l : x, ww = dx < 6 ? g.iw : w;
    box.setAttribute('x', xx); box.setAttribute('y', y); box.setAttribute('width', ww); box.setAttribute('height', h);
  });
  const end = ev => {
    if (!st) return; const g = _chart.geom; let [px, py] = svgPt(ev); px = Math.min(Math.max(px, g.m.l), g.m.l + g.iw); py = Math.min(Math.max(py, g.m.t), g.m.t + g.ih);
    const dx = Math.abs(px - st[0]), dy = Math.abs(py - st[1]), s0 = st; st = null; if (box) { box.remove(); box = null; }
    if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(g, Math.min(px, s0[0]), Math.max(py, s0[1])), [xb, yb] = pxData(g, Math.max(px, s0[0]), Math.min(py, s0[1]));
    setRange(dx >= 6 ? xa : null, dx >= 6 ? xb : null, dy >= 6 ? ya : null, dy >= 6 ? yb : null);
  };
  svg.addEventListener('pointerup', end); svg.addEventListener('pointercancel', () => { st = null; if (box) { box.remove(); box = null; } });
  svg.addEventListener('dblclick', resetZoom);
}
