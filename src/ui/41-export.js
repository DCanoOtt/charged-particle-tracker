/* ====== export ====== */
const LIGHT = { '--bg': '#ffffff', '--ink': '#111c26', '--ink2': '#52626f', '--line': '#d3dae1', '--accent': '#0a6f82', '--warn': '#b4540a', '--fc-red': '#d6281e', '--fc-black': '#1b1f24', '--fc-orange': '#e8741f', '--fc-blue': '#2340d6', '--fc-magenta': '#c62fd0', '--fc-green': '#2f9e2f', '--fc-purple': '#8a35cf', '--fc-brown': '#8a5a2b' };
const EXP_STYLE = 'text{font-family:Arial,Helvetica,sans-serif}.ax text{font-size:11.5px;fill:#52626f}.ax line,.ax path{stroke:#52626f;fill:none}.grid{stroke:#d3dae1;stroke-width:1}.grid.minor{stroke-opacity:.45;stroke-dasharray:2 3}.lab{font-size:12.5px;fill:#111c26}.ttl{font-size:14px;font-weight:600;fill:#111c26}.band{fill:#0a6f82;opacity:.08}.band.alt{opacity:.16}.bandl{font-size:11px;fill:#52626f}';
function exportSvgString() {
  const { o, geom } = _chart, src = _chart.svg || $('#chartBox svg'); const c = src.cloneNode(true);
  const xh = c.querySelector('#xh'); if (xh) xh.remove(); c.querySelectorAll('rect[stroke-dasharray="4 3"]').forEach(e => e.remove());
  c.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); c.setAttribute('width', geom.W); c.setAttribute('height', geom.H);
  // legend inside the plot, top right
  const items = o.series.filter(s => s.name && !s.hidden); let lg = '', y = geom.m.t + 8, wmax = 0;
  items.forEach(s => { wmax = Math.max(wmax, s.name.length); });
  if (items.length) {
    const bw = Math.min(260, 40 + wmax * 6.6), bx = geom.m.l + geom.iw - bw - 8;
    lg += '<rect x="' + bx + '" y="' + (y - 4) + '" width="' + bw + '" height="' + (items.length * 17 + 6) + '" fill="#fff" fill-opacity=".85" stroke="#d3dae1"/>';
    items.forEach((s, i) => { const yy = y + 8 + i * 17, col = s.color || 'var(--accent)'; lg += '<line x1="' + (bx + 8) + '" x2="' + (bx + 28) + '" y1="' + yy + '" y2="' + yy + '" stroke="' + col + '" stroke-width="2.4"' + (s.dash ? ' stroke-dasharray="' + s.dash + '"' : '') + '/><text x="' + (bx + 34) + '" y="' + (yy + 4) + '" style="font-size:11.5px;fill:#111c26">' + esc(s.name) + '</text>'; });
  }
  c.insertAdjacentHTML('afterbegin', '<style>' + EXP_STYLE + '</style><rect width="' + geom.W + '" height="' + geom.H + '" fill="#fff"/>');
  c.insertAdjacentHTML('beforeend', lg);
  let str = new XMLSerializer().serializeToString(c);
  str = str.replace(/var\((--[a-z0-9-]+)\)/g, (m, k) => LIGHT[k] || '#000');
  return str;
}
function svgToCanvas(widthPx) {
  return new Promise((res, rej) => {
    const str = exportSvgString(), g = _chart.geom, img = new Image();
    img.onload = () => { const cv = document.createElement('canvas'); cv.width = widthPx; cv.height = Math.round(widthPx * g.H / g.W); const cx = cv.getContext('2d'); cx.fillStyle = '#fff'; cx.fillRect(0, 0, cv.width, cv.height); cx.drawImage(img, 0, 0, cv.width, cv.height); res(cv); };
    img.onerror = () => rej(new Error('render failed')); img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(str);
  });
}
const canvasBlob = (cv, type, q) => new Promise((res, rej) => cv.toBlob(b => b ? res(b) : rej(new Error('encode failed')), type, q));
async function makePdf(cv, g) {
  const jpg = new Uint8Array(await (await canvasBlob(cv, 'image/jpeg', 0.95)).arrayBuffer());
  const PW = 720, PH = +(720 * g.H / g.W).toFixed(2), enc = new TextEncoder(), chunks = [], offs = []; let pos = 0;
  const push = u => { chunks.push(u); pos += u.length; }, ps = s => push(enc.encode(s));
  ps('%PDF-1.4\n%\xE2\xE3\xCF\xD3\n');
  const obj = (n, body) => { offs[n] = pos; ps(n + ' 0 obj\n' + body + '\nendobj\n'); };
  obj(1, '<< /Type /Catalog /Pages 2 0 R >>'); obj(2, '<< /Type /Pages /Kids [3 0 R] /Count 1 >>');
  obj(3, '<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ' + PW + ' ' + PH + '] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >>');
  const cs = 'q ' + PW + ' 0 0 ' + PH + ' 0 0 cm /Im0 Do Q'; obj(4, '<< /Length ' + cs.length + ' >>\nstream\n' + cs + '\nendstream');
  offs[5] = pos; ps('5 0 obj\n<< /Type /XObject /Subtype /Image /Width ' + cv.width + ' /Height ' + cv.height + ' /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ' + jpg.length + ' >>\nstream\n'); push(jpg); ps('\nendstream\nendobj\n');
  const xr = pos; ps('xref\n0 6\n0000000000 65535 f \n'); for (let i = 1; i <= 5; i++) ps(String(offs[i]).padStart(10, '0') + ' 00000 n \n');
  ps('trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n' + xr + '\n%%EOF\n');
  return new Blob(chunks, { type: 'application/pdf' });
}
function plotBase() { if (_chart && _chart.base) return _chart.base; const R = RES; return ('cpt_' + S.tab + '_' + (R ? R.beamName.replace(/α/g, 'alpha').replace(/³He/g, 'He3') : 'plot') + '_' + (R ? fmt(R.E0nom, 4) + 'MeV' : '')).replace(/[^A-Za-z0-9._-]+/g, '_'); }
function asciiOf() {
  const { o } = _chart, R = RES; let t = '# Charged Particle Tracker — ' + (PLOTS[S.tab].name) + '\n# title : ' + (o.title || '') + '\n# x     : ' + o.xlabel + '\n# y     : ' + o.ylabel + '\n';
  if (R) t += '# beam  : ' + R.beamName + ' ' + fmtE(R.E0nom) + ', ' + R.N + ' primaries, layers: ' + R.names.map((n, i) => (i + 1) + '=' + n + ' (' + fmtLen(R.edges[i + 1] - R.edges[i]) + ')').join('; ') + '\n';
  o.series.forEach(se => {
    if (!se.name && !se.x && !se.edges) return; t += '\n# series: ' + (se.name || 'unnamed') + '\n';
    if (se.type === 'hist') { t += '# bin_low  bin_high  bin_centre  value\n'; for (let i = 0; i < se.y.length; i++) t += se.edges[i] + '  ' + se.edges[i + 1] + '  ' + (se.edges[i] + se.edges[i + 1]) / 2 + '  ' + se.y[i] + '\n'; }
    else { t += '# x  y\n'; for (let i = 0; i < se.x.length; i++) if (isFinite(se.x[i]) && isFinite(se.y[i])) t += se.x[i] + '  ' + se.y[i] + '\n'; }
  });
  return t;
}
function csvOf() {
  let txt = ''; _chart.o.series.forEach(se => { if (!se.name && !se.x && !se.edges) return; txt += '# ' + (se.name || 'series') + '\n'; if (se.type === 'hist') { txt += 'bin_low,bin_high,value\n'; for (let i = 0; i < se.y.length; i++) txt += se.edges[i] + ',' + se.edges[i + 1] + ',' + se.y[i] + '\n'; } else { txt += 'x,y\n'; for (let i = 0; i < se.x.length; i++) if (isFinite(se.x[i]) && isFinite(se.y[i])) txt += se.x[i] + ',' + se.y[i] + '\n'; } });
  return txt;
}
function rootMacroOf() {
  const { o } = _chart, base = plotBase(), q = s => '"' + String(s).replace(/\\/g, '\\\\').replace(/"/g, '\\"') + '"';
  const ttl = (o.title || '') + ';' + (o.xlabel || '') + ';' + (o.ylabel || '');
  let c = '// ROOT macro generated by Charged Particle Tracker — ' + PLOTS[S.tab].name + '\n// Save this file as ' + base + '.C and run:  root -l -b -q ' + base + '.C   → writes ' + base + '.root\n// (histograms: TH1D, curves: TGraph; every object is also drawn on a TCanvas when run without -b)\n#include "TFile.h"\n#include "TH1D.h"\n#include "TGraph.h"\n#include "TCanvas.h"\n#include "TLegend.h"\n\nvoid ' + base.replace(/[^A-Za-z0-9_]/g, '_') + '() {\n  TFile *f = TFile::Open(' + q(base + '.root') + ', "RECREATE");\n  TCanvas *cv = new TCanvas("cv", ' + q(o.title || 'plot') + ', 900, 600);\n  bool first = true; int col = 1;\n';
  if (o.xlog) c += '  cv->SetLogx();\n'; if (o.ylog) c += '  cv->SetLogy();\n';
  let n = 0;
  o.series.forEach(se => {
    if (se.hidden) return; const id = 'o' + (n++), nm = (se.name || 'series' + n).replace(/[^A-Za-z0-9_]+/g, '_');
    if (se.type === 'hist') {
      c += '  { Double_t e[] = {' + Array.from(se.edges).map(v => +(+v).toPrecision(9)).join(',') + '};\n    TH1D *' + id + ' = new TH1D(' + q(nm) + ', ' + q(ttl) + ', ' + se.y.length + ', e);\n    Double_t v[] = {' + Array.from(se.y).map(v => +(+v).toPrecision(9)).join(',') + '};\n    for (int i = 0; i < ' + se.y.length + '; i++) ' + id + '->SetBinContent(i + 1, v[i]);\n    ' + id + '->SetStats(0); ' + id + '->SetLineColor(col++); ' + id + '->Write(); ' + id + '->Draw(first ? "HIST" : "HIST SAME"); first = false; }\n';
    } else if (se.type !== 'pts' || true) {
      const xs = [], ys = []; for (let i = 0; i < se.x.length; i++) if (isFinite(se.x[i]) && isFinite(se.y[i])) { xs.push(+(+se.x[i]).toPrecision(9)); ys.push(+(+se.y[i]).toPrecision(9)); }
      if (!xs.length) return;
      c += '  { Double_t x[] = {' + xs.join(',') + '}, y[] = {' + ys.join(',') + '};\n    TGraph *' + id + ' = new TGraph(' + xs.length + ', x, y); ' + id + '->SetName(' + q(nm) + '); ' + id + '->SetTitle(' + q(ttl) + ');\n    ' + id + '->SetLineColor(col++); ' + id + '->Write(); ' + id + '->Draw(first ? "AL" : "L SAME"); first = false; }\n';
    }
  });
  c += '  f->Close();\n}\n'; return c;
}
async function saveFile(name, data, label, msgSel) {
  const m = $(msgSel || '#csvMsg'); if (!m) return;
  try {
    const d = window.claude && claude.use ? await claude.use('downloads') : null;
    if (!d) { const blob = data instanceof Blob ? data : new Blob([data]), u = URL.createObjectURL(blob), l = document.createElement('a'); l.href = u; l.download = name; document.body.appendChild(l); l.click(); l.remove(); setTimeout(() => URL.revokeObjectURL(u), 4000); m.textContent = 'Download started: ' + name + ' (if nothing happens, use “Copy CSV”).'; return; }
    const r = await d.save({ filename: name, data }); m.textContent = r.status === 'saved' ? 'Saved ' + name + '.' : 'Sent.';
  } catch (e) { m.textContent = e && e.code === 'declined' ? 'Save cancelled.' : 'Could not save (' + ((e && (e.code || e.message)) || e) + ').'; }
}
async function doExport(kind) {
  if (!_chart) return; const m = $('#csvMsg'); m.textContent = 'Preparing…'; const base = plotBase();
  try {
    const px = +($('#expRes') ? $('#expRes').value : 1800) || 1800;
    if (kind === 'png') { const cv = await svgToCanvas(px); await saveFile(base + '.png', await canvasBlob(cv, 'image/png'), 'PNG'); }
    else if (kind === 'pdf') { const cv = await svgToCanvas(Math.max(px, 2400)); await saveFile(base + '.pdf', await makePdf(cv, _chart.geom), 'PDF'); }
    else if (kind === 'svg') await saveFile(base + '.svg', exportSvgString(), 'SVG');
    else if (kind === 'txt') await saveFile(base + '.txt', asciiOf(), 'ASCII');
    else if (kind === 'csv') await saveFile(base + '.csv', csvOf(), 'CSV');
    else if (kind === 'rootbin') {
      const hs = rootHistograms(); if (!hs.length) { m.textContent = 'This plot has no series that can be stored as histograms.'; return; }
      const rn = base + '.root', rb = RootIO.file(rn, hs), inClaude = !!(window.claude && claude.use && await claude.use('downloads'));
      // the viewer only accepts standard media types, so there the binary file is delivered inside a zip
      if (inClaude) await saveFile(base + '_histograms.zip', zipStore(rn, rb), 'ROOT file', '#csvMsg'); else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file');
      m.textContent = 'Saved ' + (inClaude ? base + '_histograms.zip (containing ' + rn + ')' : rn) + ' with ' + hs.length + ' histograms.';
    }
    else if (kind === 'root') await saveFile(base + '.C.txt', rootMacroOf(), 'ROOT macro');
  } catch (e) { m.textContent = 'Export failed: ' + e.message; }
}
document.addEventListener('click', ev => {
  const t = ev.target.closest('[data-exp]'); if (t) { doExport(t.dataset.exp); return; }
  const z = ev.target.closest('[data-zoom]'); if (z) { const a = z.dataset.zoom; if (a === 'in') zoomBy(0.5); else if (a === 'out') zoomBy(2); else resetZoom(); }
});
