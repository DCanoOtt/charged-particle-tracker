/* ====== controls shared by all plots ======
   One look and one order everywhere: [Y axis] [X axis] [Reset zoom] hint  →  plot  →  legend [Unselect all] [Select all]  →  readout  →  Save plot: PNG PDF SVG ASCII CSV ROOT file ROOT macro [PNG width].
   Each plot family keeps its own data-* hooks (passed in as attribute strings); only the markup, the labels and the list of export formats are common. */
const PC_FORMATS = [
  ['png', 'PNG', 'Raster image (width set by “PNG width”)'],
  ['pdf', 'PDF', 'Single-page PDF'],
  ['svg', 'SVG', 'Vector graphics'],
  ['txt', 'ASCII (.txt)', 'Plain-text columns'],
  ['csv', 'CSV', 'Comma-separated values'],
  ['rootbin', 'ROOT file', 'Binary ROOT file with one histogram per curve'],
  ['root', 'ROOT macro', 'C++ macro that creates the .root file when run with ROOT']
];
let EXPPX = 1800; // PNG width of all plots (px); every “PNG width” selector shows and sets it
const expPx = () => EXPPX;
function pcSel(label, attrs, cur) {
  // “Y axis [Linear|Logarithmic]”; cur: 'log' | 'lin' (or true for log)
  const lg = cur === true || cur === 'log';
  return (
    '<label class="pc-sel"><span>' +
    label +
    '</span><select ' +
    attrs +
    '><option value="lin"' +
    (lg ? '' : ' selected') +
    '>Linear</option><option value="log"' +
    (lg ? ' selected' : '') +
    '>Logarithmic</option></select></label>'
  );
}
function pcReset(attr, zoomed) {
  return '<button class="btn sm" ' + attr + (zoomed ? '' : ' disabled') + ' title="Reset the zoom of this plot (or double-click it)">Reset zoom</button>';
}
function pcTop(o) {
  // o: { y: [attrs, cur], x: [attrs, cur], reset: [attr, zoomed], extra: html, hint: text | null }
  return (
    '<div class="pcbar">' +
    (o.y ? pcSel('Y axis', o.y[0], o.y[1]) : '') +
    (o.x ? pcSel('X axis', o.x[0], o.x[1]) : '') +
    (o.reset ? pcReset(o.reset[0], o.reset[1]) : '') +
    (o.extra || '') +
    (o.hint === null ? '' : '<span class="hint">' + (o.hint || 'Drag a box to zoom (a thin band zooms one axis) · double-click resets') + '</span>') +
    '</div>'
  );
}
function pcLegendTools(offAttr, onAttr) {
  // the two buttons that follow the legend entries
  return (
    '<button class="btn sm" ' +
    offAttr +
    ' title="Hide all curves of this plot">Unselect all</button><button class="btn sm" ' +
    onAttr +
    ' title="Show all curves of this plot">Select all</button>'
  );
}
function pcSave(o) {
  // o: { attr: 'data-asave', title, formats: [kinds] | null, extra: [[kind, label, title]], msg: html span | '', selId: id of the PNG-width select }
  const fm = PC_FORMATS.filter(f => !o.formats || o.formats.includes(f[0]));
  return (
    '<div class="pcsave"><span class="lbl">' +
    (o.title || 'Save plot') +
    '</span>' +
    fm.map(f => '<button class="btn sm" ' + o.attr + '="' + f[0] + '" title="' + f[2] + '">' + f[1] + '</button>').join('') +
    (o.extra || [])
      .map(f => '<button class="btn sm" ' + o.attr + '="' + f[0] + '"' + (f[2] ? ' title="' + f[2] + '"' : '') + '>' + f[1] + '</button>')
      .join('') +
    '<label class="pc-sel"><span>PNG width</span><select data-expres' +
    (o.selId ? ' id="' + o.selId + '"' : '') +
    '>' +
    [1200, 1800, 3000].map(w => '<option value="' + w + '"' + (w === EXPPX ? ' selected' : '') + '>' + w + ' px</option>').join('') +
    '</select></label>' +
    (o.msg || '') +
    '</div>'
  );
}
document.addEventListener('change', ev => {
  // the PNG-width selectors stay in step
  const t = ev.target.closest && ev.target.closest('[data-expres]');
  if (!t) return;
  EXPPX = +t.value || 1800;
  document.querySelectorAll('[data-expres]').forEach(s => {
    s.value = String(EXPPX);
  });
});
