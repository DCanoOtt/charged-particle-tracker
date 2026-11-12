/* ====== run ====== */
let RES = null, RUN = null, RUNSIG = '';
const WORKER_TAIL = "\nlet __sim=null,__cancel=false;\nfunction __loop(){ if(__cancel){postMessage({type:'result',res:__sim.result(),partial:true});return;} let done; try{done=__sim.runChunk(60);}catch(e){postMessage({type:'error',msg:String(e&&e.message||e)});return;} postMessage({type:'progress',done:__sim.done(),total:__sim.total}); if(done){postMessage({type:'result',res:__sim.result()});} else setTimeout(__loop,0);}\nonmessage=function(e){const m=e.data; if(m.cmd==='cancel'){__cancel=true;return;} try{const layers=m.layers.map(l=>({t:l.t,mat:l.def?PT.buildMaterial(l.def,m.lib):null})); const cfg=m.cfg; cfg.layers=layers; cfg.beam.sp=m.spec; __sim=PT.simulate(cfg); __cancel=false; __loop();}catch(err){postMessage({type:'error',msg:String(err&&err.message||err)});}};\n";
function buildPayload() {
  const b = S.beam, p = S.phys, s = S.sim, spec = beamSpec();
  if (spec.kind === 'ion' && !(spec.Z >= 1 && spec.Z <= 92 && spec.A >= spec.Z)) throw new Error('Ion needs 1 ≤ Z ≤ 92 and A ≥ Z.');
  const E0 = beamE0(); if (!(E0 > 0)) throw new Error('Beam energy must be positive.');
  if (spec.kind === 'e' && E0 > 100) throw new Error('Electron energy is limited to 100 MeV.');
  const d = dirNorm(); if (!(d[2] > 1e-6)) throw new Error('The beam direction needs dz > 0 (it must travel towards +z).');
  const layers = S.layers.map((l, i) => { const t = layerCm(l); if (!(t > 0)) throw new Error('Layer ' + (i + 1) + ' needs a positive thickness.'); return { t, def: l.def }; });
  const N = Math.max(1, Math.min(5000000, Math.round(s.N)));
  const memMB = N * (22 + 8 * layers.length) / 1e6; if (memMB > 1500) throw new Error('This many primaries needs about ' + Math.round(memMB) + ' MB of per-primary tallies with ' + layers.length + ' layers; reduce the number of primaries (or run several seeds).');
  const lim = Math.round(s.N) > 5000000 ? 'Number of primaries limited to 5,000,000. ' : '';
  if (lim) runErr(lim);
  const dm = pnum(s.depthMax);
  const cfg = {
    beam: { E0, fwhm: b.fwhm / 100, spot: b.spot / 10, div: b.div / 1000, dir: d }, N, seed: Math.round(s.seed) || 1, dispTracks: Math.round(s.dispTracks), trackCap: Math.round(Math.min(Math.max(s.trackCapM || 3, 0.1), 12) * 1e6), depthBins: Math.round(s.depthBins), depthMax: dm > 0 ? dm / 10 : 0,
    phys: { Tcut: p.Tcut / 1000, delta: p.delta, deltaTrack: p.deltaTrack, reactions: p.reactions, rxScale: p.rxScale, scatter: p.scatter, straggle: p.straggle, stepFrac: p.stepFrac, EcutI: p.EcutI / 1000, EcutE: p.EcutE / 1000, nucStop: p.nucStop, stopSet: p.stopSet || 'auto' },
    xslib: alphaXsLib(), xslibSp: PJ() ? { Z: PJ().Z, A: PJ().A } : null,
  };
  return { cmd: 'run', layers, cfg, spec, lib: libMap() };
}
function setRunUI(running) { $('#runBtn').disabled = running; $('#cancelBtn').disabled = !running; }
function runErr(msg) { const e = $('#runErr'); e.hidden = !msg; e.textContent = msg || ''; }
async function startRun() {
  if (alibIsAlpha() && alibF().size) { $('#runStat').textContent = 'Reading cross sections…'; try { await alibEnsure(); } catch (e) { } }
  runErr(''); let pl;
  try { pl = buildPayload(); RUNSIG = resSig(); } catch (e) { runErr(e.message); return; }
  const t0 = performance.now(); setRunUI(true); $('#runStat').textContent = 'Building stopping tables…'; $('#progBar').style.width = '0%';
  const onMsg = m => {
    if (m.type === 'progress') { $('#progBar').style.width = (100 * m.done / m.total).toFixed(1) + '%'; $('#runStat').textContent = m.done.toLocaleString('en') + ' / ' + m.total.toLocaleString('en') + ' primaries'; }
    else if (m.type === 'result') { finishRun(m.res, pl, (performance.now() - t0) / 1000, m.partial); }
    else if (m.type === 'error') { runErr(m.msg); setRunUI(false); $('#runStat').textContent = 'Failed.'; RUN = null; }
    else if (m.type === 'cancelled') { setRunUI(false); $('#runStat').textContent = 'Cancelled.'; RUN = null; }
  };
  let usedWorker = false;
  try {
    const blob = new Blob([ENGINE_SRC + WORKER_TAIL], { type: 'text/javascript' }), url = URL.createObjectURL(blob), w = new Worker(url);
    w.onmessage = e => { onMsg(e.data); if (e.data.type !== 'progress') { w.terminate(); URL.revokeObjectURL(url); } };
    w.onerror = e => { onMsg({ type: 'error', msg: 'Worker error: ' + (e.message || 'unknown') }); try { w.terminate(); } catch (x) { } };
    w.postMessage(pl); RUN = { cancel: () => { w.postMessage({ cmd: 'cancel' }); } }; usedWorker = true;
  } catch (e) { usedWorker = false; }
  if (!usedWorker) {   // fallback: run on the main thread in slices
    let sim, cancelled = false;
    try { const layers = pl.layers.map(l => ({ t: l.t, mat: l.def ? PT.buildMaterial(l.def, pl.lib) : null })); pl.cfg.layers = layers; pl.cfg.beam.sp = pl.spec; sim = PT.simulate(pl.cfg); } catch (e) { onMsg({ type: 'error', msg: e.message }); return; }
    RUN = { cancel: () => { cancelled = true; } };
    const loop = () => { if (cancelled) { onMsg({ type: 'result', res: sim.result(), partial: true }); return; } let done; try { done = sim.runChunk(40); } catch (e) { onMsg({ type: 'error', msg: e.message }); return; } onMsg({ type: 'progress', done: sim.done(), total: sim.total }); if (done) onMsg({ type: 'result', res: sim.result() }); else setTimeout(loop, 0); };
    setTimeout(loop, 10);
  }
}
function finishRun(res, pl, secs, partial) {
  RUN = null; setRunUI(false);
  const L = res.L; res.layerE = new Float64Array(L); for (let i = 0; i < L; i++) res.layerE[i] = res.lay.depEl[i] + res.lay.depNu[i] + res.lay.depSe[i];
  res.names = S.layers.map(layerName); res.units = S.layers.map(l => l.tu); res.secs = secs;
  res.beamName = beamName(); res.E0nom = pl.cfg.beam.E0; res.rho = S.layers.map(l => { try { return l.def ? matOf(l.def).rho : 0; } catch (e) { return 0; } });
  RES = res; RES.sig = RUNSIG || resSig(); RES.pj = PJ() ? PJ().letter : ''; $('#runStat').textContent = (partial ? 'Stopped early — results from the first ' : 'Done: ') + res.N.toLocaleString('en') + ' primaries in ' + secs.toFixed(2) + ' s.'; $('#progBar').style.width = partial ? (100 * res.N / pl.cfg.N).toFixed(1) + '%' : '100%';
  drawStack(); renderSummary(); renderPlots(true);
}
function loadSetupText(txt, msg) {
  try { const o = JSON.parse(txt); if (!o.layers || !o.beam) throw new Error('Not a setup file.'); if (Array.isArray(o.userMaterials)) o.userMaterials.forEach(m => { const k = USERLIB.findIndex(p => p.name === m.name); if (k >= 0) USERLIB[k] = m; else USERLIB.push(m); }); applyState(o); RES = null; renderAll(); $(msg).textContent = 'Setup loaded.'; saveLocal(); } catch (e) { $(msg).textContent = 'Could not load: ' + e.message; }
}
document.addEventListener('change', ev => {
  if (ev.target.id !== 'setupFile') return; const f = ev.target.files[0]; if (!f) return;
  const r = new FileReader(); r.onload = () => { loadSetupText(String(r.result), '#setupMsg2'); ev.target.value = ''; }; r.onerror = () => { $('#setupMsg2').textContent = 'Could not read the file.'; }; r.readAsText(f);
});
document.addEventListener('click', ev => {
  const t = ev.target;
  if (t.id === 'runBtn') startRun();
  else if (t.id === 'cancelBtn' && RUN) RUN.cancel();
  else if (t.id === 'setupCopy') { const txt = setupJSON(); $('#setupTxt').value = txt; $('#setupTxt').select(); (navigator.clipboard && navigator.clipboard.writeText ? navigator.clipboard.writeText(txt) : Promise.reject()).then(() => $('#setupMsg').textContent = 'Copied to clipboard.', () => $('#setupMsg').textContent = 'Select the text above and copy it (Ctrl/Cmd+C).'); }
  else if (t.id === 'setupLoad') loadSetupText($('#setupTxt').value, '#setupMsg');
  else if (t.id === 'setupSave') { const nm = ((S.layers.map(l => l.def && l.def.name).filter(Boolean).join('_') || 'setup').replace(/[^\w.-]+/g, '_').slice(0, 50)) + '.setup.json'; saveFile(nm, setupJSON(), 'setup', '#setupMsg2'); }
  else if (t.id === 'setupOpen') $('#setupFile').click();
  else if (t.id === 'setupReset') { S = defaultState(); RES = null; renderAll(); saveLocal(); }
});
