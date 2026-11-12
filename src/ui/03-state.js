/* ====== state ====== */
const PARTICLES = [
  { id: 'e', name: 'Electron (e⁻)', kind: 'e' }, { id: 'p', name: 'Proton', Z: 1, A: 1 }, { id: 'd', name: 'Deuteron', Z: 1, A: 2 },
  { id: 't', name: 'Triton', Z: 1, A: 3 }, { id: 'he3', name: 'Helion (³He)', Z: 2, A: 3 }, { id: 'a', name: 'α particle (⁴He)', Z: 2, A: 4 },
  { id: 'li7', name: '⁷Li', Z: 3, A: 7 }, { id: 'c12', name: '¹²C', Z: 6, A: 12 }, { id: 'n14', name: '¹⁴N', Z: 7, A: 14 }, { id: 'o16', name: '¹⁶O', Z: 8, A: 16 },
  { id: 'ne20', name: '²⁰Ne', Z: 10, A: 20 }, { id: 'ar40', name: '⁴⁰Ar', Z: 18, A: 40 }, { id: 'fe56', name: '⁵⁶Fe', Z: 26, A: 56 },
  { id: 'kr84', name: '⁸⁴Kr', Z: 36, A: 84 }, { id: 'xe132', name: '¹³²Xe', Z: 54, A: 132 }, { id: 'au197', name: '¹⁹⁷Au', Z: 79, A: 197 }, { id: 'u238', name: '²³⁸U', Z: 92, A: 238 },
  { id: 'ion', name: 'Other ion (set Z, A)…', Z: 2, A: 4, custom: true },
];
function defaultState() {
  return {
    layers: [{ id: uid(), t: 100, tu: 'µm', def: presetDef('Silicon') }],
    beam: { pid: 'a', Z: 2, A: 4, E: 5.486, Eu: 'MeV', fwhm: 0, spot: 0, div: 0, dir: [0, 0, 1] },
    phys: { delta: true, Tcut: 20, deltaTrack: true, reactions: true, rxScale: 1, scatter: true, straggle: true, nucStop: true, stopSet: 'auto', stepFrac: 0.1, EcutI: 1, EcutE: 1 },
    sim: { N: 2000, seed: 12345, dispTracks: 40, depthBins: 400, depthMax: '', trackCapM: 3 },
    plots: {}, irr: {}, tab: 'depth', sel: 0,
  };
}
let S = defaultState();
// Signature of everything that changes the transport result: when it differs from the one of the run on screen, the results are out of date and are cleared
// (otherwise the ENDF-6 reaction folding and the activation would combine the old tracks with the new layers, beam or physics).
function resSig() { return JSON.stringify([S.layers.map(l => [l.t, l.tu, l.def]), S.beam, S.phys, [S.sim.seed, S.sim.depthBins, S.sim.depthMax, S.sim.trackCapM]]); }
function resCheck() {
  if (typeof RES === 'undefined' || !RES || typeof RUN !== 'undefined' && RUN) return; if (RES.sig === resSig()) return;
  RES = null; setTimeout(() => { try { renderAll(); const rs = $('#runStat'); if (rs) rs.textContent = 'The setup changed: the previous results were cleared, press Run.'; } catch (e) { } }, 0);
}
function saveLocal() { try { resCheck(); } catch (e) { } try { if (typeof alibSync === 'function') alibSync(); } catch (e) { }
  try { localStorage.setItem('chargedtracker.v1', JSON.stringify({ S: { layers: S.layers, beam: S.beam, phys: S.phys, sim: S.sim, plots: S.plots, irr: S.irr, tab: S.tab }, USERLIB })); } catch (e) { } }
function loadLocal() {
  try {
    const j = JSON.parse(localStorage.getItem('chargedtracker.v1') || 'null'); if (!j) return;
    if (Array.isArray(j.USERLIB)) USERLIB = j.USERLIB;
    if (j.S && j.S.layers && j.S.layers.length) { applyState(j.S); }
  } catch (e) { }
}
function applyState(o) {
  const d = defaultState();
  S.layers = (o.layers || d.layers).map(l => ({ id: uid(), t: +l.t > 0 ? +l.t : 1, tu: LEN[l.tu] || l.tu === 'mg/cm²' || l.tu === 'g/cm²' ? l.tu : 'µm', def: l.def ? l.def : null }));
  S.beam = Object.assign(d.beam, o.beam || {}); S.phys = Object.assign(d.phys, o.phys || {}); S.sim = Object.assign(d.sim, o.sim || {});
  S.plots = o.plots || {}; S.irr = o.irr || {}; S.tab = o.tab === 'espec' ? 'flux' : (o.tab || 'depth'); S.sel = 0;
}
