/* ====== material modal ====== */
const PU = { kPa: 1, bar: 100, atm: 101.325, mbar: 0.1, Pa: 0.001, Torr: 0.133322368 };   // unit → kPa
let M = null;      // working copy: {idx, t, tu, def|null}
function openModal(i) {
  const l = S.layers[i];
  M = { idx: i, t: l.t, tu: l.tu, def: l.def ? clone(l.def) : null, open: {}, msg: '' };
  if (!M.def) M.def = presetDef('Silicon'), M.vac = true;
  renderModal(); $('#modal').hidden = false;
}
// the window closes on a click on the dark backdrop only when the mouse was pressed AND released there:
// selecting text in a field and letting go outside the window must not close it
let _mdownOnBackdrop = false;
let _mdownId = '';                                                         // a click closes a window only if the mouse button also went down on the backdrop (selecting text by dragging out of a box must not close it)
document.addEventListener('mousedown', ev => { _mdownId = ev.target.id || ''; _mdownOnBackdrop = ev.target.id === 'modal'; }, true);
function closeModal() { $('#modal').hidden = true; M = null; }
function normRows(d) { d.rows = d.rows || []; return d; }
const ELNAMES = 'Hydrogen Helium Lithium Beryllium Boron Carbon Nitrogen Oxygen Fluorine Neon Sodium Magnesium Aluminium Silicon Phosphorus Sulfur Chlorine Argon Potassium Calcium Scandium Titanium Vanadium Chromium Manganese Iron Cobalt Nickel Copper Zinc Gallium Germanium Arsenic Selenium Bromine Krypton Rubidium Strontium Yttrium Zirconium Niobium Molybdenum Technetium Ruthenium Rhodium Palladium Silver Cadmium Indium Tin Antimony Tellurium Iodine Xenon Caesium Barium Lanthanum Cerium Praseodymium Neodymium Promethium Samarium Europium Gadolinium Terbium Dysprosium Holmium Erbium Thulium Ytterbium Lutetium Hafnium Tantalum Tungsten Rhenium Osmium Iridium Platinum Gold Mercury Thallium Lead Bismuth Polonium Astatine Radon Francium Radium Actinium Thorium Protactinium Uranium Neptunium Plutonium Americium Curium Berkelium Californium Einsteinium Fermium Mendelevium Nobelium Lawrencium'.split(' ');   // English element names, index = Z-1
// a material made of one single element takes the element's name (until the user types another one)
// force = the user just changed the composition: then a preset name or another element's name is replaced too
function autoName(d, force) {
  const one = d.rows.length === 1 && d.rows[0].kind !== 'mat' ? PT.SYMZ[d.rows[0].sym] : 0, nm = one ? ELNAMES[one - 1] : '';
  if (!M) return;
  if (force || M.force) { M.force = false; if (nm && d.name !== nm && (!d.name || d.name === 'Material' || d.name === M.auto || ELNAMES.includes(d.name) || PRESETS.some(p => p.name === d.name))) { d.name = nm; M.auto = nm; return; } }
  if (nm && (!d.name || d.name === 'Material' || d.name === M.auto)) { d.name = nm; M.auto = nm; }
  else if (!nm && M.auto && d.name === M.auto) { d.name = 'Material'; M.auto = ''; }
}
function syncNameFields(d) { const a = $('#modal [data-m="name"]'), b = $('#saveName'); if (a && a !== document.activeElement) a.value = d.name || ''; if (b && b !== document.activeElement) b.value = d.name || ''; }
let _mTarget = null;
function renderModal() {
  const d = normRows(M.def); autoName(d);
  const isGas = d.phase === 'gas', hasMat = d.rows.some(r => r.kind === 'mat');
  const syms = PT.SYMS;
  let prev = null, perr = '';
  try { prev = matOf(d); } catch (e) { perr = e.message; }
  const lib = libMap();
  const groups = ['Solids', 'Liquids', 'Gases'];
  const opts = '<option value="">Load a preset…</option>' + groups.map(g => '<optgroup label="' + g + '">' + PRESETS.filter(p => p.g === g).map(p => '<option value="p:' + esc(p.name) + '">' + esc(p.name) + '</option>').join('') + '</optgroup>').join('') +
    (USERLIB.length ? '<optgroup label="My materials">' + USERLIB.map(p => '<option value="p:' + esc(p.name) + '">' + esc(p.name) + '</option>').join('') + '</optgroup>' : '') + '<option value="vac">Vacuum (empty gap)</option>';
  const rowsH = d.rows.map((r, i) => {
    const kind = r.kind === 'mat' ? 'mat' : 'el';
    const comp = kind === 'el' ?
      '<input type="text" list="symlist" data-m="sym" data-i="' + i + '" value="' + esc(r.sym || '') + '" placeholder="Element, e.g. Fe" autocomplete="off" style="text-transform:none">' :
      '<select data-m="key" data-i="' + i + '"><option value="">Choose material…</option>' + Object.keys(lib).map(k => '<option' + (k === r.key ? ' selected' : '') + ' value="' + esc(k) + '">' + esc(k) + '</option>').join('') + '</select>';
    const sym = r.sym, Z = PT.SYMZ[sym];
    let iso = '';
    if (kind === 'el' && Z && M.open[i]) {
      const nat = PT.ELTAB[Z].iso, cur = (d.iso && d.iso[sym]) ? d.iso[sym] : nat.map(a => [a[0], a[1]]);
      const sum = cur.reduce((s, a) => s + a[1], 0);
      iso = '<div class="isot"><div class="hint" style="grid-column:1/-1">Isotopic abundances of ' + esc(sym) + ' (atom %, normalised on use; natural by default)' + (d.iso && d.iso[sym] ? ' <b style="color:var(--warn)">— modified</b>' : '') + ' · sum = ' + fmt(sum, 6) + '</div>' +
        nat.map(a => { const c = cur.find(q => q[0] === a[0]); return '<div class="ir"><span>' + a[0] + esc(sym) + '</span><input type="number" min="0" step="any" data-m="iso" data-i="' + i + '" data-a="' + a[0] + '" value="' + (c ? c[1] : 0) + '"><button class="btn sm" data-ma="only" data-i="' + i + '" data-a="' + a[0] + '">100%</button></div>'; }).join('') +
        '<div style="grid-column:1/-1"><button class="btn sm" data-ma="natural" data-i="' + i + '">Reset to natural</button></div></div>';
    }
    return '<div class="crow"><select data-m="kind" data-i="' + i + '"><option value="el"' + (kind === 'el' ? ' selected' : '') + '>Element</option><option value="mat"' + (kind === 'mat' ? ' selected' : '') + (d.mode === 'atoms' ? ' disabled' : '') + '>Material</option></select>' + comp +
      '<input type="number" min="0" step="any" data-m="amt" data-i="' + i + '" value="' + r.amt + '" title="' + (d.mode === 'atoms' ? 'Number of atoms per formula unit' : 'Mass fraction or mass (normalised)') + '">' +
      (kind === 'el' ? '<button class="btn sm" data-ma="iso" data-i="' + i + '">' + (M.open[i] ? 'Isotopes ▾' : 'Isotopes ▸') + '</button>' : '<span></span>') +
      '<button class="btn sm danger" data-ma="del" data-i="' + i + '" aria-label="Remove component">✕</button></div>' + iso;
  }).join('');
  let pv = '';
  if (perr) pv = '<p class="err">' + esc(perr) + '</p>';
  else if (prev) {
    pv = '<div class="kpis" style="grid-template-columns:repeat(auto-fit,minmax(140px,1fr))"><div class="kpi"><small>Density</small><b>' + fmt(prev.rho, 4) + '</b><span>g/cm³ · ' + esc(prev.rhoSrc) + '</span></div>' +
      '<div class="kpi"><small>⟨Z/A⟩</small><b>' + fmt(prev.ZA, 5) + '</b><span>electrons: ' + fmt(prev.nEl, 3) + ' /cm³</span></div>' +
      '<div class="kpi"><small>Mean excitation I</small><b>' + fmt(prev.I, 4) + ' eV</b><span>Bragg: ' + fmt(prev.Ibragg, 4) + ' eV</span></div>' +
      '<div class="kpi"><small>Radiation length</small><b>' + fmt(prev.X0, 4) + '</b><span>g/cm² = ' + fmt(prev.X0 / prev.rho, 4) + ' cm</span></div></div>' +
      '<div class="tscroll"><table class="t"><tr><th>Element</th><th>Mass fraction</th><th>Atom fraction</th><th>Mean A</th><th>Isotopes (abundance %)</th></tr>' +
      prev.elems.map(e => '<tr><td>' + esc(e.sym) + ' (Z=' + e.Z + ')</td><td>' + fmt(e.w, 5) + '</td><td>' + fmt(e.x, 5) + '</td><td>' + fmt(e.A, 6) + '</td><td style="white-space:normal;text-align:left">' + e.iso.filter(a => a.ab > 0).map(a => a.A + ': ' + fmt(a.ab * 100, 4)).join(' · ') + '</td></tr>').join('') + '</table></div>';
  }
  const vac = M.vac;
  (_mTarget || $('#modal')).innerHTML = '<div class="mbox" role="dialog" aria-modal="true" aria-label="Layer editor">' +
    '<h3>Layer ' + (M.idx + 1) + '</h3>' +
    '<div class="mgrid"><label class="fld"><span>Thickness</span><input type="number" min="0" step="any" data-m="t" value="' + M.t + '"></label>' +
    '<label class="fld"><span>Unit</span><select data-m="tu">' + LENU.concat(['mg/cm²', 'g/cm²']).map(u => '<option' + (u === M.tu ? ' selected' : '') + '>' + u + '</option>').join('') + '</select></label>' +
    '<label class="fld"><span>Material preset</span><select data-m="preset">' + opts + '</select></label></div>' +
    (vac ? '<p class="hint">This layer is a <b>vacuum gap</b> (no interaction). Choose a preset or press “Define a material” to give it a material.</p><div class="btns"><button class="btn" data-ma="unvac">Define a material</button></div>' :
      '<div class="mgrid"><label class="fld wide"><span>Material name</span><input type="text" data-m="name" value="' + esc(d.name || '') + '"></label></div>' +
      '<div class="btns"><span class="lbl">Form</span><div class="seg" role="group" aria-label="Phase">' + ['solid', 'liquid', 'gas'].map(ph => '<button data-ma="phase" data-v="' + ph + '" class="' + (d.phase === ph ? 'sel' : '') + '">' + ph[0].toUpperCase() + ph.slice(1) + '</button>').join('') + '</div>' +
      '<span class="lbl" style="margin-left:10px">Composition by</span><div class="seg" role="group" aria-label="Composition mode">' + [['atoms', 'Atoms (formula)'], ['mass', 'Mass fractions']].map(q => '<button data-ma="mode" data-v="' + q[0] + '" class="' + (d.mode === q[0] ? 'sel' : '') + '">' + q[1] + '</button>').join('') + '</div></div>' +
      '<div class="mgrid"><label class="fld wide"><span>Quick entry from chemical formula, e.g. Al2O3, Ca(OH)2, C8H8</span><div style="display:flex;gap:8px"><input type="text" id="formulaIn" placeholder="H2O" style="text-transform:none"><button class="btn" data-ma="formula">Set composition</button></div></label></div>' +
      '<datalist id="symlist">' + syms.map(s => '<option value="' + s + '">').join('') + '</datalist>' +
      '<div style="display:flex;flex-direction:column;gap:8px"><div class="lbl">Components (natural isotopic composition unless edited)</div>' + rowsH + '<div class="btns"><button class="btn sm" data-ma="addel">+ Element</button><button class="btn sm" data-ma="addmat"' + (d.mode === 'atoms' ? ' disabled title="Switch to mass fractions to mix materials"' : '') + '>+ Material (composite)</button></div>' + (hasMat && d.mode === 'atoms' ? '<p class="err">Composites need mass-fraction mode.</p>' : '') + '</div>' +
      '<div class="mgrid">' +
      (isGas ? '<label class="chk wide"><input type="checkbox" data-m="rhoAuto"' + (d.rhoAuto !== false ? ' checked' : '') + '> Density from the ideal-gas law (P, T)</label>' +
        '<label class="fld"><span>Pressure</span><div style="display:flex;gap:6px"><input type="number" min="0" step="any" data-m="P" value="' + (+(((d.gas && d.gas.P) || 101.325) / (PU[(d.gas && d.gas.Pu) || 'kPa'] || 1)).toPrecision(7)) + '"><select data-m="Pu" style="width:90px">' + Object.keys(PU).map(u => '<option' + (u === ((d.gas && d.gas.Pu) || 'kPa') ? ' selected' : '') + '>' + u + '</option>').join('') + '</select></div></label>' +
        '<label class="fld"><span>Temperature (K)</span><input type="number" min="0" step="any" data-m="T" value="' + ((d.gas && d.gas.T) || 293.15) + '"></label>' +
        '<label class="fld"><span>Molar mass (g/mol)</span><input type="number" min="0" step="any" data-m="molMass" value="' + (d.molMass || '') + '" placeholder="' + (d.mode === 'atoms' && prev ? fmt(prev.formulaMass, 5) : 'mean atomic mass') + '"></label>' : '') +
      '<label class="fld"><span>Density (g/cm³)' + (isGas ? ' — used when ideal-gas is off' : '') + '</span><input type="number" min="0" step="any" data-m="rho" value="' + (d.rho || '') + '" placeholder="' + (prev ? 'auto ' + fmt(prev.rho, 4) : 'required') + '"></label>' +
      '<label class="fld"><span>Mean excitation energy I (eV), optional</span><input type="number" min="0" step="any" data-m="I" value="' + (d.I || '') + '" placeholder="' + (prev ? fmt(prev.Ibragg, 4) + ' (Bragg)' : '') + '"></label></div>' +
      '<div class="preview" id="mpreview">' + pv + '</div>' +
      '<div class="mgrid"><label class="fld"><span>Save as “My material” under name</span><div style="display:flex;gap:8px"><input type="text" id="saveName" value="' + esc(d.name || '') + '"><button class="btn" data-ma="save">Save</button></div></label></div>' +
      (USERLIB.length ? '<div class="chips">' + USERLIB.map((p, i) => '<span class="chip">' + esc(p.name) + '<button data-ma="dellib" data-i="' + i + '" aria-label="Delete saved material ' + esc(p.name) + '">✕</button></span>').join('') + '</div>' : '')) +
    (M.msg ? '<p class="' + (M.msgBad ? 'err' : 'hint') + '">' + esc(M.msg) + '</p>' : '') +
    '<div class="btns" style="justify-content:flex-end"><button class="btn" data-ma="cancel">Cancel</button><button class="btn pri" data-ma="apply"' + (!vac && perr ? ' disabled' : '') + '>Apply</button></div></div>';
}
let _modalTimer = null;
function modalInput(t) {
  const d = M.def, k = t.dataset.m, i = +t.dataset.i; const num = parseFloat(t.value);
  let rerender = false;
  switch (k) {
    case 't': M.t = isFinite(num) ? num : 0; break;
    case 'tu': M.tu = t.value; break;
    case 'name': { const sv = $('#saveName'); if (sv && sv.value === (d.name || '')) sv.value = t.value; d.name = t.value; M.auto = ''; break; }
    case 'sym': d.rows[i].sym = t.value.trim().replace(/^./, c => c.toUpperCase()); break;
    case 'key': d.rows[i].key = t.value; rerender = true; break;
    case 'amt': d.rows[i].amt = isFinite(num) ? num : 0; break;
    case 'kind': { const r = d.rows[i]; if (t.value === 'mat') { r.kind = 'mat'; r.key = ''; delete r.sym; } else { r.kind = 'el'; r.sym = r.sym || 'H'; delete r.key; } rerender = true; break; }
    case 'iso': { const sym = d.rows[i].sym, Z = PT.SYMZ[sym]; d.iso = d.iso || {}; const nat = PT.ELTAB[Z].iso; const cur = d.iso[sym] || nat.map(a => [a[0], a[1]]); const a = +t.dataset.a; const c = cur.find(q => q[0] === a); c[1] = isFinite(num) && num > 0 ? num : 0; d.iso[sym] = cur; break; }
    case 'rhoAuto': d.rhoAuto = t.checked; break;
    case 'P': d.gas = d.gas || { P: 101.325, T: 293.15 }; d.gas.P = isFinite(num) ? num * (PU[d.gas.Pu || 'kPa'] || 1) : 0; d.rhoAuto = true; break;
    case 'T': d.gas = d.gas || { P: 101.325, T: 293.15 }; d.gas.T = isFinite(num) ? num : 0; d.rhoAuto = true; break;
    case 'Pu': d.gas = d.gas || { P: 101.325, T: 293.15 }; d.gas.Pu = t.value; d.rhoAuto = true; break;
    case 'molMass': d.molMass = isFinite(num) && num > 0 ? num : 0; break;
    case 'rho': d.rho = isFinite(num) && num > 0 ? num : 0; break;
    case 'I': d.I = isFinite(num) && num > 0 ? num : 0; break;
  }
  return rerender;
}
function updatePreviewOnly() {                    // refresh the derived parts (preview, 'auto' placeholders, messages) WITHOUT touching the inputs being typed in
  if (!M || M.vac || !$('#mpreview')) return;
  const tmp = document.createElement('div'); _mTarget = tmp; try { renderModal(); } finally { _mTarget = null; }
  const pv = $('#mpreview', tmp); if (pv) $('#mpreview').innerHTML = pv.innerHTML;
  const oi = $$('#modal input, #modal button[data-ma="apply"]'), ni = $$('input, button[data-ma="apply"]', tmp);
  if (oi.length === ni.length) oi.forEach((o, k) => { const n = ni[k]; if (n.placeholder !== o.placeholder) o.placeholder = n.placeholder; if (o.type === 'checkbox' && o.checked !== n.checked) o.checked = n.checked; if (o.dataset.ma === 'apply') o.disabled = n.disabled; });
}
document.addEventListener('input', ev => {
  const t = ev.target; if (!M || !t.dataset || !t.dataset.m || t.dataset.m === 'preset') return;
  if (t.tagName === 'SELECT') return;
  const re = modalInput(t); if (t.dataset.m === 'sym') { autoName(M.def, true); syncNameFields(M.def); } clearTimeout(_modalTimer); _modalTimer = setTimeout(re ? renderModal : updatePreviewOnly, 250);
});
document.addEventListener('change', ev => {
  const t = ev.target; if (!M || !t.dataset || !t.dataset.m) return;
  if (t.dataset.m === 'preset') {
    const v = t.value; if (!v) return;
    if (v === 'vac') { M.vac = true; M.msg = ''; }
    else { const nm = v.slice(2); M.def = presetDef(nm); M.vac = false; M.msg = ''; }
    renderModal(); return;
  }
  if (t.tagName === 'SELECT' || t.type === 'checkbox') { if (['kind', 'key'].includes(t.dataset.m)) M.force = true; modalInput(t); clearTimeout(_modalTimer); renderModal(); }
});
document.addEventListener('click', ev => {
  const b = ev.target.closest('[data-ma]'); if (!b || !M) { if (ev.target.id === 'modal' && _mdownOnBackdrop) closeModal(); return; }
  const a = b.dataset.ma, d = M.def, i = +b.dataset.i;
  if (a === 'cancel') { closeModal(); return; }
  if (a === 'apply') {
    const l = S.layers[M.idx];
    if (M.vac) l.def = null; else { try { matOf(d); } catch (e) { M.msg = e.message; M.msgBad = true; renderModal(); return; } l.def = clone(d); if (!l.def.name) l.def.name = 'Material'; }
    l.t = M.t > 0 ? M.t : 1; l.tu = M.tu; closeModal(); drawStack(); saveLocal(); return;
  }
  if (a === 'unvac') { M.vac = false; }
  else if (a === 'phase') { d.phase = b.dataset.v; if (d.phase === 'gas' && !d.gas) d.gas = { P: 101.325, T: 293.15 }; }
  else if (a === 'mode') {
    if (b.dataset.v === 'atoms' && d.rows.some(r => r.kind === 'mat')) { M.msg = 'Remove composite components before switching to atom counts.'; M.msgBad = true; renderModal(); return; }
    if (b.dataset.v !== d.mode) {
      if (b.dataset.v === 'mass') { try { const m = matOf(d); d.rows = d.rows.filter(r => r.sym).map(r => { const e = m.elems.find(q => q.sym === r.sym); return { kind: 'el', sym: r.sym, amt: e ? +(e.w * 100).toPrecision(5) : r.amt }; }); } catch (e) { } }
      else { try { const m = matOf(d); d.rows = m.elems.map(e => ({ kind: 'el', sym: e.sym, amt: +(e.x * 100).toPrecision(5) })); } catch (e) { } }
      d.mode = b.dataset.v;
    }
  }
  else if (a === 'formula') {
    try { const f = parseFormula($('#formulaIn').value); d.rows = Object.keys(f).map(s => ({ kind: 'el', sym: s, amt: f[s] })); d.mode = 'atoms'; d.molMass = 0; autoName(d, true); if (!d.name || d.name === 'Material') d.name = $('#formulaIn').value.trim(); M.msg = ''; } catch (e) { M.msg = 'Formula: ' + e.message; M.msgBad = true; }
  }
  else if (a === 'addel') { M.force = true; d.rows.push({ kind: 'el', sym: 'H', amt: 1 }); }
  else if (a === 'addmat') d.rows.push({ kind: 'mat', key: '', amt: 1 });
  else if (a === 'del') { M.force = true; d.rows.splice(i, 1); M.open = {}; }
  else if (a === 'iso') M.open[i] = !M.open[i];
  else if (a === 'only') { const sym = d.rows[i].sym, Z = PT.SYMZ[sym]; d.iso = d.iso || {}; d.iso[sym] = PT.ELTAB[Z].iso.map(q => [q[0], q[0] === +b.dataset.a ? 100 : 0]); }
  else if (a === 'natural') { const sym = d.rows[i].sym; if (d.iso) { delete d.iso[sym]; if (!Object.keys(d.iso).length) delete d.iso; } }
  else if (a === 'save') {
    const nm = ($('#saveName').value || '').trim(); if (!nm) { M.msg = 'Give the material a name first.'; M.msgBad = true; renderModal(); return; }
    try { matOf(d); } catch (e) { M.msg = e.message; M.msgBad = true; renderModal(); return; }
    const c = clone(d); c.name = nm; delete c.g; const k = USERLIB.findIndex(p => p.name === nm); if (k >= 0) USERLIB[k] = c; else USERLIB.push(c);
    d.name = nm; M.msg = 'Saved “' + nm + '” to My materials (kept in this browser).'; M.msgBad = false; saveLocal();
  }
  else if (a === 'dellib') { USERLIB.splice(i, 1); saveLocal(); }
  else return;
  if (a !== 'save' && a !== 'formula') { /* keep message */ }
  renderModal();
});
