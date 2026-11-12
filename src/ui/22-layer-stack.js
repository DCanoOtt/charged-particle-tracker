/* ====== layer stack ====== */
function drawStack() {
  const L = S.layers; if (S.sel >= L.length) S.sel = L.length - 1;
  const H = 210, y0 = 40, bh = 120; let x = 86; let cum = 0;
  const infos = L.map((l, i) => {
    const cm = safeLayerCm(l); const w = isFinite(cm) && cm > 0 ? Math.max(116, Math.min(152, 116 + 6 * Math.log10(Math.max(cm / 1e-6, 1)))) : 116;
    let rho = null; try { if (l.def) rho = matOf(l.def).rho; } catch (e) { }
    return { l, cm, w, rho };
  });
  const d = dirNorm(), ang = Math.atan2(d[0], d[2]);
  let svg = '<defs><pattern id="pLiq" width="14" height="8" patternUnits="userSpaceOnUse"><path d="M0 4 Q3.5 0 7 4 T14 4" fill="none" stroke="var(--fc-blue)" stroke-width="1" opacity=".55"/></pattern><pattern id="pGas" width="12" height="12" patternUnits="userSpaceOnUse"><circle cx="3" cy="3" r="1" fill="var(--ink2)" opacity=".5"/><circle cx="9" cy="9" r="1" fill="var(--ink2)" opacity=".5"/></pattern><pattern id="pVac" width="8" height="8" patternUnits="userSpaceOnUse"><path d="M0 8L8 0" stroke="var(--line)" stroke-width="1"/></pattern></defs>';
  // beam arrow
  const ay = y0 + bh / 2;
  svg += '<line x1="12" y1="' + (ay - 28 * Math.sin(ang)) + '" x2="' + (12 + 56 * Math.cos(ang)) + '" y2="' + (ay + 28 * Math.sin(ang)) + '" stroke="var(--fc-red)" stroke-width="3"/>';
  const tipx = 12 + 56 * Math.cos(ang), tipy = ay + 28 * Math.sin(ang), ca = Math.cos(ang), sa = Math.sin(ang);
  const hp = [[tipx, tipy], [tipx - 12 * ca + 5 * sa, tipy - 12 * sa - 5 * ca], [tipx - 12 * ca - 5 * sa, tipy - 12 * sa + 5 * ca]];
  svg += '<polygon points="' + hp.map(p => p[0].toFixed(1) + ',' + p[1].toFixed(1)).join(' ') + '" fill="var(--fc-red)"/>';
  svg += '<text x="12" y="' + (y0 - 12) + '" class="sm" style="font-weight:600;fill:var(--fc-red)">' + esc(beamName()) + ' · ' + esc(fmtE(beamE0())) + '</text>';
  svg += '<text x="12" y="' + (y0 + bh + 32) + '" class="sm">z = 0</text>';
  infos.forEach((o, i) => {
    const l = o.l, def = l.def; let fill = 'url(#pVac)';
    if (def) { if (def.phase === 'gas') fill = 'url(#pGas)'; else if (def.phase === 'liquid') fill = 'url(#pLiq)'; }
    let base = 'var(--field)';
    if (def && def.phase === 'solid') { const r = o.rho || 2; base = 'color-mix(in srgb, var(--accent) ' + Math.round(14 + 52 * Math.min(Math.log10(1 + r) / 1.35, 1)) + '%, var(--field))'; }
    else if (def && def.phase === 'liquid') base = 'color-mix(in srgb, var(--fc-blue) 16%, var(--field))';
    const nm = layerName(l, i); const lab = nm.length > Math.floor((o.w - 8) / 6.8) ? nm.slice(0, Math.floor((o.w - 8) / 6.8) - 1) + '…' : nm;
    let dep = ''; if (RES && RES.L === L.length && RES.layerE) { const e = RES.layerE[i]; dep = fmtE(e / RES.N); }
    svg += '<clipPath id="lc' + i + '"><rect x="' + x + '" y="' + y0 + '" width="' + o.w + '" height="' + bh + '"/></clipPath><g class="lay' + (i === S.sel ? ' sel' : '') + '" data-i="' + i + '" tabindex="0" role="button" aria-label="Edit layer ' + (i + 1) + ': ' + esc(nm) + '"><title>Click to edit layer ' + (i + 1) + '</title>' +
      '<rect class="body" x="' + x + '" y="' + y0 + '" width="' + o.w + '" height="' + bh + '" fill="' + base + '" rx="2"/><rect x="' + x + '" y="' + y0 + '" width="' + o.w + '" height="' + bh + '" fill="' + fill + '" rx="2" pointer-events="none"/>' +
      '<text clip-path="url(#lc' + i + ')" x="' + (x + o.w / 2) + '" y="' + (y0 + 34) + '" text-anchor="middle" style="font-weight:600;font-size:12px">' + esc(lab) + '</text>' +
      '<text clip-path="url(#lc' + i + ')" x="' + (x + o.w / 2) + '" y="' + (y0 + 52) + '" text-anchor="middle" class="sm">' + esc(fmt(l.t, 4) + ' ' + l.tu) + '</text>' +
      '<text clip-path="url(#lc' + i + ')" x="' + (x + o.w / 2) + '" y="' + (y0 + 68) + '" text-anchor="middle" class="sm">' + (def ? esc(def.phase) : 'gap') + '</text>' +
      '<text clip-path="url(#lc' + i + ')" x="' + (x + o.w / 2) + '" y="' + (y0 + 83) + '" text-anchor="middle" class="sm">' + (o.rho ? esc(fmt(o.rho, 3) + ' g/cm³' + (o.P ? ' · ' + o.P : '')) : '') + '</text>' +
      '<text x="' + (x + o.w / 2) + '" y="' + (y0 + 108) + '" text-anchor="middle" class="sm" style="font-weight:600">' + (i + 1) + '</text></g>';
    if (dep) svg += '<text x="' + (x + o.w / 2) + '" y="' + (y0 + bh + 16) + '" text-anchor="middle" class="sm" style="fill:var(--accent)">E<tspan baseline-shift="sub" font-size="8">dep</tspan> ' + esc(dep) + '</text>';
    cum += isFinite(o.cm) ? o.cm : 0;
    svg += '<text x="' + (x + o.w) + '" y="' + (y0 + bh + 32) + '" text-anchor="middle" class="sm">' + esc(fmtLen(cum)) + '</text>';
    x += o.w;
  });
  svg += '<g class="addb" id="addLayer" tabindex="0" role="button" aria-label="Add layer"><rect x="' + (x + 14) + '" y="' + y0 + '" width="84" height="' + bh + '" rx="4"/><text x="' + (x + 56) + '" y="' + (y0 + bh / 2 - 2) + '" text-anchor="middle" style="font-size:26px">+</text><text x="' + (x + 57) + '" y="' + (y0 + bh / 2 + 18) + '" text-anchor="middle">Add layer</text></g>';
  const W = x + 112;
  const sel = S.layers[S.sel];
  $('#cStack').innerHTML = '<h2>Target layers <span class="hint" style="text-transform:none;letter-spacing:0">beam enters at z = 0 and travels towards +z · E<sub>dep</sub> = energy deposited per primary · click a layer to define its material</span></h2>' +
    '<div class="stackwrap"><svg class="stack" viewBox="0 0 ' + W + ' ' + H + '" style="width:' + W + 'px;max-width:none;height:auto" role="group" aria-label="Layer stack">' + svg + '</svg></div>' +
    '<div class="laybar"><span>Selected: <b>' + (S.sel + 1) + ' · ' + esc(layerName(sel, S.sel)) + '</b></span>' +
    '<button class="btn sm" data-lay="edit">Edit material…</button><button class="btn sm" data-lay="left"' + (S.sel === 0 ? ' disabled' : '') + '>◀ Move</button><button class="btn sm" data-lay="right"' + (S.sel === L.length - 1 ? ' disabled' : '') + '>Move ▶</button><button class="btn sm" data-lay="dup">Duplicate</button><button class="btn sm danger" data-lay="del"' + (L.length === 1 ? ' disabled' : '') + '>Remove</button></div>';
}
function addLayer() {
  S.layers.push({ id: uid(), t: 100, tu: 'µm', def: presetDef('Aluminium') });
  S.sel = S.layers.length - 1; drawStack(); saveLocal(); openModal(S.sel);
}
document.addEventListener('click', ev => {
  const t = ev.target;
  const lay = t.closest('.lay'); if (lay && $('#cStack').contains(lay)) { S.sel = +lay.dataset.i; drawStack(); openModal(S.sel); return; }
  if (t.closest('#addLayer')) { addLayer(); return; }
  const b = t.closest('[data-lay]'); if (b) {
    const a = b.dataset.lay, L = S.layers;
    if (a === 'edit') openModal(S.sel);
    else if (a === 'left' && S.sel > 0) { [L[S.sel - 1], L[S.sel]] = [L[S.sel], L[S.sel - 1]]; S.sel--; }
    else if (a === 'right' && S.sel < L.length - 1) { [L[S.sel + 1], L[S.sel]] = [L[S.sel], L[S.sel + 1]]; S.sel++; }
    else if (a === 'dup') { const c = clone(L[S.sel]); c.id = uid(); L.splice(S.sel + 1, 0, c); S.sel++; }
    else if (a === 'del' && L.length > 1) { L.splice(S.sel, 1); S.sel = Math.max(0, S.sel - 1); }
    drawStack(); saveLocal();
  }
});
document.addEventListener('keydown', ev => {
  if ((ev.key === 'Enter' || ev.key === ' ') && ev.target.closest) {
    const lay = ev.target.closest('.lay'); if (lay) { ev.preventDefault(); S.sel = +lay.dataset.i; drawStack(); openModal(S.sel); }
    else if (ev.target.closest('#addLayer')) { ev.preventDefault(); addLayer(); }
  }
  if (ev.key === 'Escape' && !$('#modal').hidden) closeModal();
});
