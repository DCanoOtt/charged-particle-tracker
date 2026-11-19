/* ====== γ spectrum: Ge detector and DAQ configuration window (continued from 49-gamma.js) ====== */
// ---- window "Configure the Ge detector and DAQ" (parameters used by the γ-spectrum windows) ----
function gamCfgSummary(g) {
  const eff =
    g.effMode === 'ge'
      ? 'Ge detector Ø' + g.geD + ' × ' + g.geL + ' mm, ' + +g.geDist + ' mm from the endcap'
      : g.effMode === 'curve'
        ? 'efficiency curve'
        : 'constant efficiency ' + g.eff;
  return (
    'FWHM ' +
    g.fwhm +
    ' keV' +
    (g.fmode === 'sqrt' ? ' at 1 MeV (∝ √E)' : '') +
    ' · ' +
    eff +
    ' · pulse width τ = ' +
    g.tau +
    ' μs' +
    (g.puCurve ? ' · pile-up spectrum on' : '')
  );
}
function geSketch(g) {
  // side view of the point source and the detector (schematic, not to scale); the labels show the current values
  const T = (x, y, t, a, c) =>
      '<text x="' + x + '" y="' + y + '" text-anchor="' + (a || 'middle') + '" font-size="11.5" fill="' + (c || 'var(--ink)') + '">' + t + '</text>',
    dim = (x0, x1, y, lab, val) =>
      '<line x1="' +
      x0 +
      '" y1="' +
      y +
      '" x2="' +
      x1 +
      '" y2="' +
      y +
      '" stroke="var(--accent)" stroke-width="1.4" marker-start="url(#gsa)" marker-end="url(#gsa)"/>' +
      T((x0 + x1) / 2, y + 15, lab) +
      T((x0 + x1) / 2, y + 28, val, 'middle', 'var(--ink2)'),
    vdim = (x, y0, y1, lab, val) =>
      '<line x1="' +
      x +
      '" y1="' +
      y0 +
      '" x2="' +
      x +
      '" y2="' +
      y1 +
      '" stroke="var(--accent)" stroke-width="1.4" marker-start="url(#gsa)" marker-end="url(#gsa)"/>' +
      T(x + 8, (y0 + y1) / 2 - 4, lab, 'start') +
      T(x + 8, (y0 + y1) / 2 + 11, val, 'start', 'var(--ink2)'),
    ext = (x, y0, y1) => '<line x1="' + x + '" y1="' + y0 + '" x2="' + x + '" y2="' + y1 + '" stroke="var(--ink2)" stroke-width=".8" stroke-dasharray="2 3"/>';
  return (
    '<svg viewBox="0 0 660 340" role="img" aria-label="Geometry of the point source and the Ge detector" style="width:100%;max-width:660px;display:block;margin:4px auto 0">' +
    '<defs><marker id="gsa" viewBox="0 0 8 8" refX="4" refY="4" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L8 4 L0 8 z" fill="var(--accent)"/></marker></defs>' +
    '<rect x="190" y="52" width="320" height="196" rx="6" fill="none" stroke="var(--ink2)" stroke-width="1.2" stroke-dasharray="5 3"/>' + // endcap housing
    '<rect x="186" y="52" width="5" height="196" fill="var(--ink2)" opacity=".55"/>' + // Al window (front of the endcap)
    '<rect x="230" y="70" width="240" height="160" fill="var(--accent-soft, rgba(60,140,160,.18))" stroke="var(--accent)" stroke-width="1.6"/>' + // Ge crystal
    '<rect x="230" y="70" width="9" height="160" fill="var(--warn, #c77d1a)" opacity=".75"/><rect x="230" y="70" width="240" height="6" fill="var(--warn, #c77d1a)" opacity=".75"/><rect x="230" y="224" width="240" height="6" fill="var(--warn, #c77d1a)" opacity=".75"/>' + // dead layer (front and side)
    '<line x1="20" y1="150" x2="500" y2="150" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="8 4 2 4"/>' + // axis
    '<line x1="40" y1="150" x2="230" y2="70" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="3 3"/>' + // limiting photon
    '<circle cx="40" cy="150" r="5" fill="var(--bad, #c0392b)"/>' +
    T(40, 135, 'point source') +
    T(40, 170, 'on the axis') +
    T(349, 126, 'Ge crystal (active volume)') +
    '<path d="M300 76 L300 40 L330 40" fill="none" stroke="var(--ink2)" stroke-width=".9"/>' +
    T(335, 44, 'dead layer t<tspan baseline-shift="sub" font-size="9">d</tspan> = ' + +g.geDead + ' mm (front and side, Ge)', 'start', 'var(--ink)') +
    '<path d="M188 252 L188 268 L150 268" fill="none" stroke="var(--ink2)" stroke-width=".9"/>' +
    T(146, 272, 'Al window w = ' + +g.geWin + ' mm', 'end') +
    ext(40, 160, 262) +
    ext(190, 248, 262) +
    ext(230, 232, 262) +
    ext(470, 232, 262) +
    dim(40, 190, 280, 'source – endcap distance', +g.geDist + ' mm') +
    dim(190, 230, 280, 'gap', +g.geGap + ' mm') +
    dim(230, 470, 280, 'crystal length L', +g.geL + ' mm') +
    vdim(500, 70, 230, 'crystal diameter Ø D', +g.geD + ' mm') +
    '<line x1="470" y1="70" x2="506" y2="70" stroke="var(--ink2)" stroke-width=".8" stroke-dasharray="2 3"/><line x1="470" y1="230" x2="506" y2="230" stroke="var(--ink2)" stroke-width=".8" stroke-dasharray="2 3"/>' +
    T(350, 334, 'Side view, schematic (not to scale). The source is at z = source–endcap distance + gap from the crystal front.', 'middle', 'var(--ink2)') +
    '</svg>'
  );
}
function gamCfgRender() {
  const w = $('#geCfgWin'),
    g = gamSet();
  if (g.effMode === 'curve' && !String(g.curve).trim()) {
    g.curve = gamDefaultCurve(g);
    saveLocal();
  }
  const keep = w.querySelector('.mbox'),
    scrollY = keep ? keep.scrollTop : 0;
  const inp = (k, ex) => '<input type="number" data-gk="' + k + '" step="any" value="' + g[k] + '" ' + (ex || '') + '>',
    sel = (k, o2) =>
      '<select data-gk="' +
      k +
      '">' +
      o2.map(x => '<option value="' + x[0] + '"' + (String(g[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') +
      '</select>';
  const fld = (l, h) => '<label class="fld"><span>' + l + '</span>' + h + '</label>';
  w.innerHTML =
    '<div class="mbox" style="width:min(820px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:10px"><h3>Configure the Ge detector and DAQ</h3><button class="btn sm" data-gcfgclose>Done</button></div>' +
    '<p class="hint">These parameters are used by the γ-ray spectra of the time slices (window opened from the activation plots). Changes are applied at once and saved with the setup.</p>' +
    '<h4>Energy resolution</h4><div class="mgrid">' +
    fld('Energy resolution (Gaussian), FWHM in keV', inp('fwhm', 'min="0"')) +
    fld(
      'FWHM model',
      sel('fmode', [
        ['const', 'constant'],
        ['sqrt', 'FWHM given at 1 MeV, ∝ √E']
      ])
    ) +
    '</div>' +
    '<h4>Detector efficiency</h4><div class="mgrid">' +
    fld(
      'Detector efficiency model',
      sel('effMode', [
        ['const', 'constant'],
        ['curve', 'curve ε(E)'],
        ['ge', 'Ge detector (geometry model)']
      ])
    ) +
    (g.effMode === 'const' ? fld('Constant efficiency (0 – 1)', inp('eff', 'min="0" step="any"')) : '') +
    (g.effMode === 'ge'
      ? fld('Crystal diameter (mm)', inp('geD', 'min="1"')) +
        fld('Crystal length (mm)', inp('geL', 'min="1"')) +
        fld('Source – endcap distance (mm)', inp('geDist', 'min="0"')) +
        fld('Endcap – crystal gap (mm)', inp('geGap', 'min="0"')) +
        fld('Dead layer, Ge (mm)', inp('geDead', 'min="0"')) +
        fld('Window, Al (mm)', inp('geWin', 'min="0"')) +
        fld('Efficiency scale factor', inp('geScale', 'min="0"')) +
        '<div style="grid-column:1/-1">' +
        geSketch(g) +
        '</div><div class="hint" style="grid-column:1/-1" id="geCfgInfo">' +
        geInfo(g) +
        '</div>'
      : '') +
    (g.effMode === 'curve'
      ? '<label class="fld" style="grid-column:1/-1"><span>Efficiency curve: one line per energy, “energy (keV)  total efficiency  peak efficiency” (log-log interpolation, constant outside the range). A line with only two values is “energy  peak efficiency”; a missing total is taken from the Ge geometry model. Leave it empty to use the Ge geometry model (its parameters are the ones set under “Ge detector (geometry model)”).</span><textarea data-gk="curve" rows="6" spellcheck="false" placeholder="100 0.06 0.05&#10;662 0.03 0.025&#10;1332 0.016 0.013" style="font-family:var(--font-mono);width:100%">' +
        esc(g.curve) +
        '</textarea></label><div class="btns" style="grid-column:1/-1"><button class="btn sm" data-gcurve="default">Load the default curve</button><button class="btn sm" data-gcurve="savedef">Save this curve as my default</button><button class="btn sm" data-gcurve="ge">Fill from the Ge geometry model</button><button class="btn sm" data-gcurve="cleardef"' +
        (String(g.curveDef || '').trim() ? '' : ' disabled') +
        '>Forget my default</button><span class="hint" id="gcurveMsg">Default: ' +
        (String(g.curveDef || '').trim() ? 'your saved curve' : 'the values of the Ge geometry model (energy, total, peak), listed when the box is empty') +
        '.</span></div>'
      : '') +
    '</div>' +
    '<div id="geCfgEff"></div>' +
    '<h4>DAQ: pulse width and pile-up</h4><div class="mgrid">' +
    fld('Pulse width τ (μs), pile-up window ±τ', inp('tau', 'min="0"')) +
    fld('DAQ dead time per pulse (μs), non-paralysable — 0 = none', inp('tauD', 'min="0"')) +
    (g.effMode === 'ge'
      ? fld(
          'Pulses per full-energy count',
          sel('ftpMode', [
            ['auto', 'from the model: ε_tot / ε_peak'],
            ['man', 'manual']
          ])
        ) + (g.ftpMode === 'man' ? fld('Pulses per full-energy count (≥ 1)', inp('ftp', 'min="1"')) : '')
      : fld('Pulses per full-energy count (total/peak ≥ 1)', inp('ftp', 'min="1"'))) +
    fld(
      'Pile-up spectrum',
      sel('puCurve', [
        ['false', 'only probabilities'],
        ['true', 'add the pile-up spectrum']
      ])
    ) +
    '</div>' +
    (g.effMode === 'ge'
      ? '<p class="hint">With the Ge detector model the γ-spectrum window also shows the efficiency curves and a simulated Ge spectrum (peaks, Compton continuum, escape peaks). Models and formulas: Instructions, sections 11.8 and 11.9.</p>'
      : '') +
    '<div class="btns"><button class="btn" data-gcfgclose>Done</button></div></div>';
  const nm = w.querySelector('.mbox');
  if (nm) nm.scrollTop = scrollY;
  gamCfgEff();
}
function gamCfgChanged(structure) {
  // after a parameter changed: refresh the other windows that use it
  const g = gamSet();
  saveLocal();
  if (structure) gamCfgRender();
  else {
    const i = $('#geCfgInfo');
    if (i && g.effMode === 'ge') {
      i.innerHTML = geInfo(g);
      const sk = i.previousElementSibling;
      if (sk) sk.innerHTML = geSketch(g);
    }
    gamCfgEff();
  }
  document.querySelectorAll('[data-gcsum]').forEach(e => {
    e.textContent = gamCfgSummary(g);
  });
  const sm = $('#geCfgSum');
  if (sm) sm.textContent = gamCfgSummary(g);
  if (GAM && $('#gamWin') && !$('#gamWin').hidden) gamRender();
  if (typeof measOnCfg === 'function') measOnCfg();
}
