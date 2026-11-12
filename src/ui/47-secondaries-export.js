/* ====== saving the secondary-particle spectra: every plot on its own, or all of them at once ====== */
const SEC_NAMES = { n: 'neutrons', p: 'protons', d: 'deuterons', t: 'tritons', h: 'He3', a: 'alphas', g: 'gammas', r: 'recoils' };
function zipMulti(files) {                                              // several files, stored without compression
  const enc = new TextEncoder(), t = new Date(), parts = [], cds = [], dosT = (t.getHours() << 11) | (t.getMinutes() << 5) | (t.getSeconds() >> 1), dosD = ((t.getFullYear() - 1980) << 9) | ((t.getMonth() + 1) << 5) | t.getDate(); let off = 0, cdLen = 0;
  for (const f of files) {
    const nb = enc.encode(f.name), data = f.data, crc = crc32(data);
    const lh = new DataView(new ArrayBuffer(30)); lh.setUint32(0, 0x04034b50, true); lh.setUint16(4, 20, true); lh.setUint16(6, 0x0800, true); lh.setUint16(8, 0, true); lh.setUint16(10, dosT, true); lh.setUint16(12, dosD, true); lh.setUint32(14, crc, true); lh.setUint32(18, data.length, true); lh.setUint32(22, data.length, true); lh.setUint16(26, nb.length, true); lh.setUint16(28, 0, true);
    const ch = new DataView(new ArrayBuffer(46)); ch.setUint32(0, 0x02014b50, true); ch.setUint16(4, 20, true); ch.setUint16(6, 20, true); ch.setUint16(8, 0x0800, true); ch.setUint16(10, 0, true); ch.setUint16(12, dosT, true); ch.setUint16(14, dosD, true); ch.setUint32(16, crc, true); ch.setUint32(20, data.length, true); ch.setUint32(24, data.length, true); ch.setUint16(28, nb.length, true); ch.setUint32(42, off, true);
    parts.push(lh.buffer, nb, data); cds.push(ch.buffer, nb); off += 30 + nb.length + data.length; cdLen += 46 + nb.length;
  }
  const ed = new DataView(new ArrayBuffer(22)); ed.setUint32(0, 0x06054b50, true); ed.setUint16(8, files.length, true); ed.setUint16(10, files.length, true); ed.setUint32(12, cdLen, true); ed.setUint32(16, off, true);
  return new Blob(parts.concat(cds, [ed.buffer]), { type: 'application/zip' });
}
// chart context of one secondary plot as it is shown now (hidden angles, log axis, zoom)
function secCtx(key) {
  const d = SEC.get(key), st = SECST.get(key); if (!d || !st) return null;
  const series = d.series.map(se => Object.assign({}, se, { hidden: !!st.hidden[se.name] }));
  const o = Object.assign({}, d.o, { series, ylog: st.ylog, xmin: st.xmin, xmax: st.xmax, ymin: st.ymin, ymax: st.ymax }), r = renderChart(o), div = document.createElement('div'); div.innerHTML = r.svg;
  const [li, sp] = key.split('|');
  return { o, geom: r.geom, svg: div.querySelector('svg'), base: plotBase() + '_layer' + (+li + 1) + '_' + (SEC_NAMES[sp] || sp), key };
}
async function secBuild(ctx, kind) {                                    // → { name, data } for one plot
  const old = _chart; _chart = ctx;
  try {
    const px = +($('#expRes') ? $('#expRes').value : 1800) || 1800, b = ctx.base;
    if (kind === 'png') return { name: b + '.png', data: await canvasBlob(await svgToCanvas(px), 'image/png') };
    if (kind === 'svg') return { name: b + '.svg', data: exportSvgString() };
    if (kind === 'txt') return { name: b + '.txt', data: asciiOf() };
    if (kind === 'csv') return { name: b + '.csv', data: csvOf() };
    if (kind === 'rootbin') { const hs = rootHistograms(); return { name: b + '.root', hs }; }
  } finally { _chart = old; }
}
const secBytes = async d => d instanceof Blob ? new Uint8Array(await d.arrayBuffer()) : new TextEncoder().encode(d);
async function secSave(kind, keys, msgSel, all) {
  const m = $(msgSel); if (!m) return; m.textContent = 'Preparing…';
  try {
    const items = []; for (const k of keys) { const c = secCtx(k); if (c) items.push(await secBuild(c, kind)); }
    if (!items.length) { m.textContent = 'Nothing to save.'; return; }
    const base = plotBase() + '_secondaries', inClaude = !!(window.claude && claude.use && await claude.use('downloads'));
    if (kind === 'rootbin') {
      const hs = []; items.forEach((it, i) => it.hs.forEach(h => hs.push(all ? Object.assign({}, h, { name: 'p' + (i + 1) + '_' + h.name }) : h)));
      if (!hs.length) { m.textContent = 'No histograms to store.'; return; }
      const rn = (all ? base : items[0].name.replace(/\.root$/, '')) + '.root', rb = RootIO.file(rn, hs);
      if (inClaude) await saveFile(rn.replace(/\.root$/, '') + '_histograms.zip', zipStore(rn, rb), 'ROOT file', msgSel); else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', msgSel);
      m.textContent = 'Saved ' + rn + (inClaude ? ' (inside a zip)' : '') + ' with ' + hs.length + ' histograms.'; return;
    }
    if (!all) { const it = items[0]; await saveFile(it.name, it.data, kind, msgSel); return; }
    if (kind === 'txt' || kind === 'csv') {
      const txt = items.map(it => '# ===== ' + it.name + ' =====\n' + it.data).join('\n'); await saveFile(base + '.' + kind, txt, kind, msgSel);
      m.textContent = 'Saved ' + base + '.' + kind + ' with ' + items.length + ' spectra.'; return;
    }
    const files = []; for (const it of items) files.push({ name: it.name, data: await secBytes(it.data) });
    await saveFile(base + '_' + kind + '.zip', zipMulti(files), kind, msgSel); m.textContent = 'Saved ' + base + '_' + kind + '.zip with ' + files.length + ' files.';
  } catch (e) { m.textContent = 'Export failed: ' + (e && e.message || e); }
}
const SEC_SAVE_BTNS = '<button class="btn sm" data-ssave="png">PNG</button><button class="btn sm" data-ssave="svg">SVG</button><button class="btn sm" data-ssave="txt">ASCII (.txt)</button><button class="btn sm" data-ssave="csv">CSV</button><button class="btn sm" data-ssave="rootbin" title="Binary ROOT file with one TH1D per curve">ROOT file</button>';
document.addEventListener('click', ev => {
  const b = ev.target.closest('[data-ssave]'); if (!b) return; const all = !!b.closest('#secAll'), kind = b.dataset.ssave;
  if (all) secSave(kind, Array.from(SEC.keys()), '#secMsg', true);
  else { const el = b.closest('.secplot[data-sk]'); if (el) secSave(kind, [el.dataset.sk], '#smsg_' + el.dataset.sk.replace(/\W/g, '_'), false); }
});
