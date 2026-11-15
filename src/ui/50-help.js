/* The text of the instructions is in src/help/*.html (one file per chapter, joined by build.py into HELP_HTML); this file has the code that shows it. */
document.getElementById('help').innerHTML = HELP_HTML;
(function () {
  try {
    const g = Object.assign({}, GE_DEF),
      c = geCurves(g, 60),
      el = document.getElementById('geHelpPlot');
    if (el)
      el.innerHTML = renderChart({
        series: [
          { name: 'ε_peak', type: 'line', x: c.x, y: c.pk, color: 'var(--accent)', w: 2 },
          { name: 'ε_tot', type: 'line', x: c.x, y: c.t, color: 'var(--bad)', w: 2, dash: '6 4' }
        ],
        xlog: true,
        ylog: true,
        xmin: '20',
        xmax: '3000',
        ymin: '1e-4',
        ymax: '0.2',
        xlabel: 'Energy (keV)',
        ylabel: 'Efficiency (per emitted γ)',
        title: 'Ge detector efficiency (default geometry)'
      }).svg;
  } catch (e) {
    /* help figure is optional */
  }
})();
document.addEventListener('click', ev => {
  const g = ev.target.closest('[data-goto]');
  if (!g) return;
  ev.preventDefault();
  const t = document.getElementById(g.dataset.goto);
  if (t) t.scrollIntoView({ block: 'start' });
});
document.addEventListener('click', ev => {
  if (ev.target.id === 'helpBtn') {
    $('#help').hidden = false;
    $('#help .mbox').scrollTop = 0;
    $('#help').scrollTop = 0;
  } else if (ev.target.closest('[data-help="close"]') || (ev.target.id === 'help' && _mdownId === 'help')) $('#help').hidden = true;
});
document.addEventListener('keydown', ev => {
  if (ev.key === 'Escape' && !$('#help').hidden) $('#help').hidden = true;
});
