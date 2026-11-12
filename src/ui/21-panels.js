/* ====== panels ====== */
function renderBeam() {
  const b = S.beam, p = curPart(), isE = p.kind === 'e';
  const eu = isE ? ['eV', 'keV', 'MeV', 'GeV'] : ['keV', 'MeV', 'GeV', 'keV/u', 'MeV/u'];
  if (!eu.includes(b.Eu)) b.Eu = 'MeV';
  const d = dirNorm(), th = Math.acos(Math.min(1, Math.max(-1, d[2]))) * 180 / Math.PI;
  $('#gBeam').innerHTML = '<h2>Beam</h2>' +
    '<label class="fld"><span>Particle</span><select data-b="beam.pid" data-re="1">' + PARTICLES.map(q => '<option value="' + q.id + '"' + (q.id === b.pid ? ' selected' : '') + '>' + esc(q.name) + '</option>').join('') + '</select></label>' +
    (p.custom ? '<div class="row2"><label class="fld"><span>Atomic number Z</span><input type="number" min="1" max="92" step="1" data-b="beam.Z" value="' + b.Z + '"></label><label class="fld"><span>Mass number A</span><input type="number" min="1" max="260" step="1" data-b="beam.A" value="' + b.A + '"></label></div>' : '') +
    '<div class="row2"><label class="fld"><span>Kinetic energy</span><input type="number" min="0" step="any" data-b="beam.E" value="' + b.E + '"></label><label class="fld"><span>Unit</span><select data-b="beam.Eu" data-re="1">' + eu.map(u => '<option' + (u === b.Eu ? ' selected' : '') + '>' + u + '</option>').join('') + '</select></label></div>' +
    '<p class="hint" id="beamInfo"></p>' +
    '<div class="row3"><label class="fld"><span>Energy spread FWHM (%)</span><input type="number" min="0" step="any" data-b="beam.fwhm" value="' + b.fwhm + '"></label><label class="fld"><span>Spot radius (mm)</span><input type="number" min="0" step="any" data-b="beam.spot" value="' + b.spot + '"></label><label class="fld"><span>Divergence FWHM (mrad)</span><input type="number" min="0" step="any" data-b="beam.div" value="' + b.div + '"></label></div>' +
    '<div class="lbl">Direction vector (dx, dy, dz) — default 0, 0, 1 along +z</div>' +
    '<div class="row3"><label class="fld"><span>dx</span><input type="number" step="any" data-b="beam.dir.0" value="' + b.dir[0] + '"></label><label class="fld"><span>dy</span><input type="number" step="any" data-b="beam.dir.1" value="' + b.dir[1] + '"></label><label class="fld"><span>dz</span><input type="number" step="any" data-b="beam.dir.2" value="' + b.dir[2] + '"></label></div>' +
    '<div class="row2"><label class="fld"><span>or tilt θ from z (°, in x–z plane)</span><input type="number" step="any" id="tiltIn" value="' + (+th.toFixed(3)) + '"></label><div class="hint" id="dirInfo"></div></div>';
  updBeamInfo();
}
function updBeamInfo() {
  const el_ = $('#beamInfo'); if (!el_) return;
  const sp = beamSpec(), E0 = beamE0(); let t = 'Total kinetic energy: ' + fmtE(E0);
  try { const s = PT.makeSpecies(sp); const g = 1 + E0 / s.M; t += ' · β = ' + fmt(Math.sqrt(1 - 1 / (g * g)), 4) + (sp.kind === 'e' ? '' : ' · ' + fmt(E0 / sp.A, 4) + ' MeV/u'); } catch (e) { }
  el_.textContent = t;
  const d = dirNorm(), di = $('#dirInfo'); if (di) di.textContent = 'Normalised: (' + d.map(v => fmt(v, 3)).join(', ') + ')' + (d[2] <= 1e-6 ? ' — dz must be > 0' : '');
}
function renderPhys() {
  const p = S.phys;
  const ck = (k, t) => '<label class="chk"><input type="checkbox" data-b="phys.' + k + '"' + (p[k] ? ' checked' : '') + '> ' + t + '</label>';
  $('#gPhys').innerHTML = '<h2>Stopping powers and models</h2>' +
    '<label class="fld"><span>Stopping-power data</span><select data-b="phys.stopSet">' + [['auto', 'Newest available (ICRU 90 → ICRU 73 update → NIST/ICRU 73)'], ['nist', 'NIST PSTAR / ASTAR + ICRU 73'], ['semi', 'Semi-empirical model (no tables)']].map(o => '<option value="' + o[0] + '"' + ((p.stopSet || 'auto') === o[0] ? ' selected' : '') + '>' + o[1] + '</option>').join('') + '</select></label>' +
    '<p class="hint">ICRU 90 (2016) covers p, ⁴He, ¹²C and e⁻ in liquid water, dry air and graphite; ICRU 73 update: Li–Ar ions in water and air. Elsewhere the newest table is NIST PSTAR/ASTAR and ICRU 73. Compare the sets in “Stopping power &amp; range”.</p>' +
    ck('delta', 'Knock-on electron (δ-ray) production') +
    '<div class="row2"><label class="fld"><span>δ-ray threshold T<sub>cut</sub> (keV)</span><input type="number" min="0.1" step="any" data-b="phys.Tcut" value="' + p.Tcut + '"></label>' + '<label class="fld"><span>Step size (fraction of range)</span><input type="number" min="0.005" max="0.5" step="any" data-b="phys.stepFrac" value="' + p.stepFrac + '"></label></div>' +
    ck('deltaTrack', 'Transport the δ-electrons (otherwise deposit locally)') + ck('scatter', 'Multiple + single Coulomb scattering') + ck('straggle', 'Energy-loss straggling') + ck('nucStop', 'Nuclear (elastic recoil) stopping') + ck('reactions', 'Simple nuclear reactions (removal model)') +
    '<div class="row3"><label class="fld"><span>Reaction σ scale</span><input type="number" min="0" step="any" data-b="phys.rxScale" value="' + p.rxScale + '"></label><label class="fld"><span>Ion cut-off (keV)</span><input type="number" min="0.01" step="any" data-b="phys.EcutI" value="' + p.EcutI + '"></label><label class="fld"><span>e⁻ cut-off (keV)</span><input type="number" min="0.01" step="any" data-b="phys.EcutE" value="' + p.EcutE + '"></label></div>' +
    '<p class="hint">Reaction channels use a geometric + Coulomb-barrier cross section now; tabulated databases will replace it later.</p>';
}
function renderRun() {
  const s = S.sim;
  $('#gRun').innerHTML = '<h2>Simulation</h2>' +
    '<div class="row2"><label class="fld"><span>Primary particles</span><input type="number" min="1" max="5000000" step="1" data-b="sim.N" value="' + s.N + '"></label><label class="fld"><span>Random seed</span><input type="number" step="1" data-b="sim.seed" value="' + s.seed + '"></label></div>' +
    '<details><summary>More options</summary><div class="row3" style="margin-top:10px"><label class="fld"><span>Tracks drawn</span><input type="number" min="0" max="500" step="1" data-b="sim.dispTracks" value="' + s.dispTracks + '"></label><label class="fld"><span>Depth bins</span><input type="number" min="20" max="4000" step="1" data-b="sim.depthBins" value="' + s.depthBins + '"></label><label class="fld"><span>Depth max (mm)</span><input type="text" placeholder="auto" data-b="sim.depthMax" value="' + esc(s.depthMax) + '"></label></div><div class="row2" style="margin-top:8px"><label class="fld"><span>Stored track points (millions)</span><input type="number" min="0.1" max="12" step="any" data-b="sim.trackCapM" value="' + s.trackCapM + '"></label><p class="hint" style="align-self:end">Radial and lateral plots use these stored primary tracks; ~16 bytes per point.</p></div></details>' +
    '<div class="btns"><button class="btn pri" id="runBtn">Run simulation</button><button class="btn" id="cancelBtn" disabled>Cancel</button></div>' +
    '<div class="prog"><i id="progBar"></i></div><p class="hint" id="runStat">Ready.</p><p class="err" id="runErr" hidden></p>';
}
function renderSetup() {
  $('#gSetup').innerHTML = '<h2>Setup</h2><div class="btns" style="margin:8px 0"><button class="btn sm" id="setupSave">Save setup to file</button><button class="btn sm" id="setupOpen">Load setup from file…</button><input type="file" id="setupFile" accept=".json,application/json,text/plain" hidden></div><p class="hint" id="setupMsg2"></p><details><summary>Copy / paste the setup as JSON text</summary><div style="display:flex;flex-direction:column;gap:8px;margin-top:10px"><textarea id="setupTxt" spellcheck="false" placeholder="Paste a saved setup here, then press Load"></textarea><div class="btns"><button class="btn sm" id="setupCopy">Copy current</button><button class="btn sm" id="setupLoad">Load pasted</button><button class="btn sm" id="setupReset">Reset all</button></div><p class="hint" id="setupMsg"></p></div></details><p class="hint">Settings and your saved materials are remembered in this browser.</p>';
}
function setupJSON() { return JSON.stringify({ app: 'charged-particle-tracker', v: 1, layers: S.layers.map(l => ({ t: l.t, tu: l.tu, def: l.def })), beam: S.beam, phys: S.phys, sim: S.sim, irr: S.irr, userMaterials: USERLIB }, null, 1); }

/* generic data-b binding */
function setPath(path, val) { const k = path.split('.'); let o = S; for (let i = 0; i < k.length - 1; i++) o = o[k[i]]; o[k[k.length - 1]] = val; }
document.addEventListener('input', onBind); document.addEventListener('change', onBind);
function onBind(ev) {
  const t = ev.target; if (t.id === 'tiltIn') { if (ev.type === 'change') { const a = (parseFloat(t.value) || 0) * Math.PI / 180; S.beam.dir = [+Math.sin(a).toFixed(6), S.beam.dir[1] ? S.beam.dir[1] : 0, +Math.cos(a).toFixed(6)]; if (Math.cos(a) <= 0) S.beam.dir[2] = 0.0001; renderBeam(); drawStack(); saveLocal(); } return; }
  const path = t.dataset && t.dataset.b; if (!path) return;
  if (path.startsWith('plot.')) return;
  let v; if (t.type === 'checkbox') v = t.checked; else if (t.type === 'number') { v = parseFloat(t.value); if (!isFinite(v)) return; } else v = t.value;
  setPath(path, v);
  if (path.startsWith('beam.')) renderAlpha();
  if (path === 'phys.stopSet') {   // results were computed with other stopping-power data: clear them
    if (RES && ev.type === 'change') { RES = null; renderAll(); const rs = $('#runStat'); if (rs) rs.textContent = 'Stopping-power data changed: previous results cleared, press Run.'; }
    else if (S.tab === 'stop' && typeof drawChart === 'function') { try { drawChart(); } catch (e) { } }
  }
  if (path === 'beam.pid') { const q = curPart(); if (!q.custom && q.Z) { S.beam.Z = q.Z; S.beam.A = q.A; } renderAlpha(); if (alibF().size) alibAfterIndex(); }
  if (t.dataset.re && ev.type === 'change') { renderBeam(); }
  updBeamInfo(); drawStack(); saveLocal();
}
