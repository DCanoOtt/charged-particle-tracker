/* ====== SVG chart renderer ====== */
const SCOL = ['--fc-blue', '--fc-red', '--fc-green', '--fc-orange', '--fc-purple', '--fc-magenta', '--fc-brown', '--fc-black'];
const scol = i => 'var(' + SCOL[i % SCOL.length] + ')';
function niceStep(r, n) { const raw = r / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p; }
function pnum(s) { if (s === '' || s == null) return null; const v = parseFloat(String(s).replace(',', '.')); return isFinite(v) ? v : null; }
function tickLabel(v) { if (v === 0) return '0'; const a = Math.abs(v); if (a >= 1e5 || a < 1e-3) return v.toExponential(1).replace('.0e', 'e').replace('e+', 'e'); return String(+v.toPrecision(6)); }

/* o: {series:[{name,color,type:'line'|'hist'|'pts',x,y,edges,dash,hidden}], xlog,ylog,xmin,xmax,ymin,ymax,xlabel,ylabel,title,bands:[{x0,x1,label}],H} */
function renderChart(o) {
  const W = 900, H = o.H || 440, m = { l: 80, r: 22, t: o.title ? 38 : 16, b: 58 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const vis = o.series.filter(s => !s.hidden);
  let dx0 = Infinity, dx1 = -Infinity, dy0 = Infinity, dy1 = -Infinity, ypos = Infinity, xpos = Infinity;
  for (const s of vis) {
    const xs = s.type === 'hist' ? s.edges : s.x, ys = s.y;
    for (let i = 0; i < xs.length; i++) { const v = xs[i]; if (isFinite(v)) { if (v < dx0) dx0 = v; if (v > dx1) dx1 = v; if (v > 0 && v < xpos) xpos = v; } }
  }
  if (!isFinite(dx0)) { dx0 = 0; dx1 = 1; }
  let x0 = pnum(o.xmin), x1 = pnum(o.xmax), y0 = pnum(o.ymin), y1 = pnum(o.ymax);
  if (o.xlog) { if (x0 == null || x0 <= 0) x0 = isFinite(xpos) ? Math.pow(10, Math.floor(Math.log10(xpos))) : 1; if (x1 == null || x1 <= x0) x1 = Math.pow(10, Math.ceil(Math.log10(Math.max(dx1, x0 * 1.01)))); if (x1 <= x0) x1 = x0 * 10; }
  else { if (x0 == null) x0 = dx0; if (x1 == null) x1 = dx1; if (x1 <= x0) { x1 = x0 + (Math.abs(x0) || 1) * 0.5; } }
  // automatic y range: only the data inside the visible x window count (a zoomed plot must not be scaled by points outside it)
  for (const s of vis) {
    const hs = s.type === 'hist', xs = hs ? s.edges : s.x, ys = s.y;
    for (let i = 0; i < ys.length; i++) {
      const v = ys[i]; if (!isFinite(v)) continue;
      const a = hs ? xs[i] : xs[i], b = hs ? xs[i + 1] : xs[i];
      if (b < x0 || a > x1) continue;
      if (v < dy0) dy0 = v; if (v > dy1) dy1 = v; if (v > 0 && v < ypos) ypos = v;
    }
  }
  if (!isFinite(dy0)) { dy0 = 0; dy1 = 1; }
  if (o.ylog) {
    if (y0 == null || y0 <= 0) { y0 = isFinite(ypos) ? Math.pow(10, Math.floor(Math.log10(ypos))) : 1; if (y1 != null && y1 > 0 && y0 >= y1) y0 = Math.pow(10, Math.floor(Math.log10(y1)) - 1); }
    if (y1 == null || y1 <= y0) y1 = Math.pow(10, Math.ceil(Math.log10(Math.max(dy1, y0 * 1.01))));
    if (y1 <= y0) y1 = y0 * 10;
  } else {
    if (y0 == null) y0 = dy0 >= 0 ? 0 : dy0 - 0.05 * (dy1 - dy0);
    if (y1 == null) { const st = niceStep((dy1 - y0) || 1, 6); y1 = Math.ceil(dy1 / st * 1.0001) * st; if (y1 <= dy1) y1 += st; }
    if (y1 <= y0) y1 = y0 + 1;
  }
  const lx = v => o.xlog ? Math.log10(v) : v, ly = v => o.ylog ? Math.log10(v) : v;
  const X0 = lx(x0), X1 = lx(x1), Y0 = ly(y0), Y1 = ly(y1);
  const fx = v => m.l + (lx(v) - X0) / (X1 - X0) * iw, fy = v => m.t + ih - (ly(v) - Y0) / (Y1 - Y0) * ih;
  const ticks = (a, b, lg) => {
    const maj = [], min = [];
    if (lg) {
      const k0 = Math.floor(Math.log10(a)), k1 = Math.ceil(Math.log10(b)), span = k1 - k0, stp = Math.max(1, Math.ceil(span / 8));
      for (let k = k0; k <= k1; k++) { const v = Math.pow(10, k); if ((k - k0) % stp === 0 && v >= a * 0.999 && v <= b * 1.001) maj.push(v); if (span <= 8) for (let q = 2; q < 10; q++) { const w = v * q; if (w > a && w < b) min.push(w); } }
    } else { const st = niceStep(b - a, 7); for (let v = Math.ceil(a / st - 1e-9) * st; v <= b + st * 1e-9; v += st) maj.push(Math.abs(v) < st * 1e-9 ? 0 : v); }
    return { maj, min };
  };
  const tx = ticks(x0, x1, o.xlog), ty = ticks(y0, y1, o.ylog);
  const lg = (v, lgf) => { if (!lgf) return tickLabel(v); const k = Math.round(Math.log10(v)); return Math.abs(k) <= 2 ? tickLabel(v) : '10<tspan dy="-5" font-size="9">' + k + '</tspan>'; };
  let g = '';
  // bands
  (o.bands || []).forEach((b, i) => {
    const a = Math.max(b.x0, x0), c = Math.min(b.x1, x1); if (!(c > a)) return;
    const xa = fx(a), xc = fx(c);
    g += '<rect class="band' + (i % 2 ? ' alt' : '') + '" x="' + xa.toFixed(1) + '" y="' + m.t + '" width="' + (xc - xa).toFixed(1) + '" height="' + ih + '"/>';
    if (xc - xa > 40) g += '<text class="bandl" x="' + ((xa + xc) / 2).toFixed(1) + '" y="' + (m.t + 13) + '" text-anchor="middle">' + esc(b.label) + '</text>';
  });
  for (const v of tx.min) g += '<line class="grid minor" x1="' + fx(v).toFixed(1) + '" x2="' + fx(v).toFixed(1) + '" y1="' + m.t + '" y2="' + (m.t + ih) + '"/>';
  for (const v of ty.min) g += '<line class="grid minor" x1="' + m.l + '" x2="' + (m.l + iw) + '" y1="' + fy(v).toFixed(1) + '" y2="' + fy(v).toFixed(1) + '"/>';
  for (const v of tx.maj) g += '<line class="grid" x1="' + fx(v).toFixed(1) + '" x2="' + fx(v).toFixed(1) + '" y1="' + m.t + '" y2="' + (m.t + ih) + '"/>';
  for (const v of ty.maj) g += '<line class="grid" x1="' + m.l + '" x2="' + (m.l + iw) + '" y1="' + fy(v).toFixed(1) + '" y2="' + fy(v).toFixed(1) + '"/>';
  // series (clipped)
  const cid = 'c' + Math.random().toString(36).slice(2, 8);
  let sg = '';
  o.series.forEach((s, si) => {
    if (s.hidden) return;
    const col = s.color || scol(si);
    if (s.type === 'hist') {
      const e = s.edges, n = s.y.length; let d = '', started = false;
      const base = fy(o.ylog ? y0 : Math.max(y0, 0)).toFixed(1);
      for (let i = 0; i < n; i++) {
        const xa = fx(e[i]).toFixed(1), xb = fx(e[i + 1]).toFixed(1), yy = fy(o.ylog ? Math.max(s.y[i] > 0 ? s.y[i] : y0, y0) : Math.max(s.y[i], y0)).toFixed(1);
        if (!started) { d += 'M' + xa + ',' + base; started = true; }
        d += 'L' + xa + ',' + yy + 'L' + xb + ',' + yy;
      }
      if (started) d += 'L' + fx(e[n]).toFixed(1) + ',' + base;
      sg += '<path d="' + d + '" fill="' + col + '" fill-opacity="' + (s.fill == null ? 0.18 : s.fill) + '" stroke="' + col + '" stroke-width="1.6" stroke-linejoin="round"/>';
    } else if (s.type === 'pts') {
      let d = ''; for (let i = 0; i < s.x.length; i++) { const a = s.x[i], b = s.y[i]; if (!isFinite(a) || !isFinite(b) || (o.ylog && b <= 0) || (o.xlog && a <= 0)) continue; d += '<circle cx="' + fx(a).toFixed(1) + '" cy="' + fy(b).toFixed(1) + '" r="2.4"/>'; }
      sg += '<g fill="' + col + '">' + d + '</g>';
    } else {
      let d = '', pen = false;
      for (let i = 0; i < s.x.length; i++) {
        const a = s.x[i], b = s.y[i];
        if (!isFinite(a) || !isFinite(b) || (o.ylog && b <= 0) || (o.xlog && a <= 0)) { pen = false; continue; }
        d += (pen ? 'L' : 'M') + fx(a).toFixed(1) + ',' + fy(b).toFixed(1); pen = true;
      }
      sg += '<path d="' + d + '" fill="none" stroke="' + col + '" stroke-width="' + (s.w || 1.9) + '"' + (s.dash ? ' stroke-dasharray="' + s.dash + '"' : '') + (s.op ? ' stroke-opacity="' + s.op + '"' : '') + ' stroke-linejoin="round" stroke-linecap="round"/>';
    }
  });
  g += '<clipPath id="' + cid + '"><rect x="' + m.l + '" y="' + m.t + '" width="' + iw + '" height="' + ih + '"/></clipPath><g clip-path="url(#' + cid + ')">' + sg + '</g>';
  (o.vlines || []).forEach(v => { if (v.x >= x0 && v.x <= x1) g += '<line x1="' + fx(v.x).toFixed(1) + '" x2="' + fx(v.x).toFixed(1) + '" y1="' + m.t + '" y2="' + (m.t + ih) + '" stroke="var(--warn)" stroke-width="1.4" stroke-dasharray="6 4"/>' + (v.label ? '<text class="bandl" x="' + (fx(v.x) + 4).toFixed(1) + '" y="' + (m.t + ih - 6) + '" style="fill:var(--warn)">' + esc(v.label) + '</text>' : ''); });
  // axes
  let ax = '<rect x="' + m.l + '" y="' + m.t + '" width="' + iw + '" height="' + ih + '" fill="none" stroke="var(--ink2)"/>';
  for (const v of tx.maj) ax += '<line x1="' + fx(v).toFixed(1) + '" x2="' + fx(v).toFixed(1) + '" y1="' + (m.t + ih) + '" y2="' + (m.t + ih + 5) + '"/><text x="' + fx(v).toFixed(1) + '" y="' + (m.t + ih + 19) + '" text-anchor="middle">' + lg(v, o.xlog) + '</text>';
  for (const v of ty.maj) ax += '<line x1="' + (m.l - 5) + '" x2="' + m.l + '" y1="' + fy(v).toFixed(1) + '" y2="' + fy(v).toFixed(1) + '"/><text x="' + (m.l - 9) + '" y="' + (fy(v) + 4).toFixed(1) + '" text-anchor="end">' + lg(v, o.ylog) + '</text>';
  g += '<g class="ax">' + ax + '</g>';
  g += '<text class="lab" x="' + (m.l + iw / 2) + '" y="' + (H - 12) + '" text-anchor="middle">' + esc(o.xlabel || '') + '</text>';
  g += '<text class="lab" transform="translate(18,' + (m.t + ih / 2) + ') rotate(-90)" text-anchor="middle">' + esc(o.ylabel || '') + '</text>';
  if (o.stats && o.stats.length) {                                        // statistics box inside the plot (top left)
    const wst = Math.min(iw - 20, 14 + 6.4 * Math.max(...o.stats.map(t => t.length)));
    g += '<rect x="' + (m.l + 6) + '" y="' + (m.t + 5) + '" width="' + wst + '" height="' + (o.stats.length * 15 + 8) + '" fill="var(--bg)" fill-opacity=".82" stroke="var(--line)"/>' + o.stats.map((t, i) => '<text x="' + (m.l + 12) + '" y="' + (m.t + 21 + i * 15) + '" font-size="11.5" fill="var(--ink)">' + esc(t) + '</text>').join('');
  }
  if (o.title) g += '<text class="ttl" x="' + (m.l + iw / 2) + '" y="22" text-anchor="middle">' + esc(o.title) + '</text>';
  g += '<line id="xh" y1="' + m.t + '" y2="' + (m.t + ih) + '" stroke="var(--ink2)" stroke-dasharray="3 3" style="display:none"/>';
  return { svg: '<svg viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(o.title || o.ylabel || 'plot') + '">' + g + '</svg>', geom: { W, H, m, iw, ih, x0, x1, y0, y1, xlog: o.xlog, ylog: o.ylog } };
}
function chartReadout(o, geom, svgEl, ev) {
  const r = svgEl.getBoundingClientRect(), px = (ev.clientX - r.left) / r.width * geom.W;
  if (px < geom.m.l || px > geom.m.l + geom.iw) return null;
  const f = (px - geom.m.l) / geom.iw, X = geom.xlog ? Math.pow(10, Math.log10(geom.x0) + f * (Math.log10(geom.x1) - Math.log10(geom.x0))) : geom.x0 + f * (geom.x1 - geom.x0);
  let out = [];
  o.series.forEach((s, si) => {
    if (s.hidden || s.noread) return;
    let val = null;
    if (s.type === 'hist') { const e = s.edges; if (X >= e[0] && X <= e[e.length - 1]) { let lo = 0, hi = e.length - 2; while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (e[mid] <= X) lo = mid; else hi = mid - 1; } val = s.y[lo]; } }
    else if (s.type !== 'pts') {
      const xs = s.x; let best = -1, bd = Infinity; const n = xs.length; if (n < 4000) { for (let i = 0; i < n; i++) { const d = Math.abs(xs[i] - X); if (d < bd && isFinite(s.y[i])) { bd = d; best = i; } } if (best >= 0) val = s.y[best]; }
    }
    if (val != null && isFinite(val)) out.push(s.name + ' = ' + fmt(val, 4));
  });
  return { X, text: 'x = ' + fmt(X, 5) + (out.length ? '   |   ' + out.join('   |   ') : ''), px: geom.xlog ? px : px };
}
function histo(vals, nb, lo, hi) {
  const c = new Float64Array(nb), w = (hi - lo) / nb; if (!(w > 0)) return c;
  for (let i = 0; i < vals.length; i++) { const v = vals[i]; if (v < lo || v > hi || !isFinite(v)) continue; let k = Math.floor((v - lo) / w); if (k >= nb) k = nb - 1; c[k]++; }
  return c;
}
function edgesOf(lo, hi, nb, log) { const e = new Float64Array(nb + 1); for (let i = 0; i <= nb; i++) e[i] = log ? lo * Math.pow(hi / lo, i / nb) : lo + (hi - lo) * i / nb; return e; }
function describe(v) {
  const n = v.length; if (!n) return null; let s = 0, mn = Infinity, mx = -Infinity;
  for (let i = 0; i < n; i++) { s += v[i]; if (v[i] < mn) mn = v[i]; if (v[i] > mx) mx = v[i]; }
  const mean = s / n; let q = 0; for (let i = 0; i < n; i++) q += (v[i] - mean) ** 2; const sd = n > 1 ? Math.sqrt(q / (n - 1)) : 0;
  const srt = Float64Array.from(v).sort(); const med = n % 2 ? srt[(n - 1) / 2] : 0.5 * (srt[n / 2 - 1] + srt[n / 2]);
  return { n, mean, sd, min: mn, max: mx, med, fwhm: 2.3548 * sd, sem: sd / Math.sqrt(n) };
}
