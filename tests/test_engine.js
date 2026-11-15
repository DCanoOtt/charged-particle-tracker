// node tests/test_engine.js [--update] — transport engine: determinism, reference stopping powers and ranges, straggling, energy balance, golden output
// The engine is loaded from the build (python3 build.py first). The golden file tests/golden/engine.json holds results of fixed-seed runs; --update rewrites it.
const fs = require('fs'), vm = require('vm'), path = require('path'), assert = require('assert');
// default: the engine as built by build.py (modules joined into one script). With --esm (and esbuild installed: ESBUILD=/path/to/esbuild) the ES modules of src/engine are bundled by esbuild instead — proof that the imports/exports are complete and give the same results.
if (process.argv.includes('--esm')) {
  const { execFileSync } = require('child_process'), R = p => path.join(__dirname, '..', p), esb = process.env.ESBUILD || 'esbuild';
  const bundle = execFileSync(esb, [R('src/engine/30-exports.js'), '--bundle', '--format=iife', '--global-name=PTM', '--log-level=warning'], { encoding: 'utf8' });
  vm.runInThisContext(['data/elements.js', 'data/densities.js', 'data/stopping-tables.js'].map(f => fs.readFileSync(R(f), 'utf8')).join('') + bundle + ';globalThis.PT=PTM.PT;'); console.log('(engine bundled by esbuild from the ES modules)');
} else vm.runInThisContext(fs.readFileSync(path.join(__dirname, '../dist/engine.js'), 'utf8') + ';globalThis.PT=PT;');
const el = (sym, amt) => ({ kind: 'el', sym, amt });
const MAT = {
  Si: () => PT.buildMaterial({ mode: 'atoms', rows: [el('Si', 1)], phase: 'solid' }, {}),
  Au: () => PT.buildMaterial({ mode: 'atoms', rows: [el('Au', 1)], phase: 'solid' }, {}),
  water: () => PT.buildMaterial({ mode: 'atoms', rows: [el('H', 2), el('O', 1)], phase: 'liquid', rho: 1, I: 78 }, {}),
};
const ALPHA = { kind: 'ion', Z: 2, A: 4 }, PROTON = { kind: 'ion', Z: 1, A: 1 };
function run(sp, E0, mat, tCm, N, seed, phys) {
  const sim = PT.simulate({ beam: { sp, E0, fwhm: 0, spot: 0, div: 0, dir: [0, 0, 1] }, N, seed, depthBins: 100, phys: phys || {}, layers: [{ t: tCm, mat }] });
  while (!sim.runChunk(500)) { }
  return sim.result();
}
const stats = (r) => { let m = 0, s = 0; for (let k = 0; k < r.N; k++) m += r.endZ[k]; m /= r.N; for (let k = 0; k < r.N; k++) s += (r.endZ[k] - m) ** 2; return { mean: m, sigma: Math.sqrt(s / r.N) }; };
const near = (a, b, tol, msg) => assert.ok(Math.abs(a / b - 1) < tol, msg + ': ' + a + ' vs ' + b + ' (tolerance ' + tol * 100 + ' %)');
const at = (t, E) => { const i = Math.max(1, t.E.findIndex(e => e >= E)), f = Math.log(E / t.E[i - 1]) / Math.log(t.E[i] / t.E[i - 1]), ll = a => Math.exp(Math.log(a[i - 1]) + f * (Math.log(a[i]) - Math.log(a[i - 1]))); return { Se: ll(t.Se), R: ll(t.R), E }; };

// 1. reference stopping powers and CSDA ranges (NIST PSTAR / ASTAR values, electronic stopping, MeV cm²/g and g/cm²)
{ const w = MAT.water(), tp = PT.stopTable(w, PROTON, 0.5, 100, 80, {}), a = at(tp, 10); near(a.Se, 45.67, 0.03, 'proton 10 MeV in water, electronic stopping power'); near(a.R, 0.1230, 0.03, 'proton 10 MeV in water, CSDA range');
  const b = at(PT.stopTable(w, PROTON, 0.5, 300, 80, {}), 100); near(b.Se, 7.289, 0.03, 'proton 100 MeV in water'); near(b.R, 7.718, 0.03, 'proton 100 MeV in water, CSDA range');
  const ta = PT.stopTable(w, ALPHA, 0.5, 20, 80, {}), c = at(ta, 5); near(c.Se, 889, 0.05, 'alpha 5 MeV in water, electronic stopping power'); }
// 2. determinism: the same seed gives the same run, another seed another one
{ const si = MAT.Si(), r1 = run(ALPHA, 5, si, 0.01, 300, 11), r2 = run(ALPHA, 5, si, 0.01, 300, 11), r3 = run(ALPHA, 5, si, 0.01, 300, 12);
  assert.deepStrictEqual(Array.from(r1.endZ), Array.from(r2.endZ), 'same seed, same result'); assert.notDeepStrictEqual(Array.from(r1.endZ), Array.from(r3.endZ), 'another seed, another result'); }
// 3. mean projected range of 5 MeV alphas in Si ≈ CSDA range; straggling of 1–2.5 % of the range
{ const si = MAT.Si(), t = PT.stopTable(si, ALPHA, 0.1, 10, 60, {}), R = at(t, 5).R / si.rho, r = run(ALPHA, 5, si, 0.01, 3000, 5), s = stats(r);
  near(s.mean, R, 0.02, 'mean end depth ≈ CSDA range'); assert.ok(s.sigma / s.mean > 0.008 && s.sigma / s.mean < 0.025, 'range straggling ' + s.sigma / s.mean);
  const r0 = run(ALPHA, 5, si, 0.01, 1000, 5, { straggle: false, scatter: false, delta: false }), s0 = stats(r0); assert.ok(s0.sigma < 0.2 * s.sigma, 'without straggling and scattering the spread collapses: ' + s0.sigma + ' vs ' + s.sigma); }
// 4. energy balance: deposited + escaping = incoming; stopped particles all end inside the layer
{ const si = MAT.Si(), r = run(ALPHA, 5, si, 0.01, 1000, 3); let dep = 0; for (let k = 0; k < r.N; k++) dep += r.dEdep[k];
  near(dep, 5 * r.N, 1e-4, 'all the energy is deposited in a layer thicker than the range'); assert.strictEqual(r.esc.nStopP, r.N, 'all primaries stop'); assert.strictEqual(r.esc.nTransP + r.esc.nBackP, 0);
  const thin = run(ALPHA, 5, si, 0.001, 500, 3); assert.ok(thin.esc.nTransP >= 485 && thin.esc.nStopP <= 3, 'a 10 µm layer is crossed by 5 MeV alphas: ' + JSON.stringify(thin.esc)); }
// 5. energy-loss ordering: a heavier target stops more per cm than a lighter one; range of protons scales ~ E^1.8
{ const rSi = at(PT.stopTable(MAT.Si(), PROTON, 0.5, 50, 60, {}), 5).R, rAu = at(PT.stopTable(MAT.Au(), PROTON, 0.5, 50, 60, {}), 5).R; assert.ok(rAu / MAT.Au().rho < rSi / MAT.Si().rho, 'protons stop in a shorter distance in gold than in silicon');
  const w = MAT.water(), t = PT.stopTable(w, PROTON, 0.5, 300, 80, {}), p = Math.log(at(t, 100).R / at(t, 10).R) / Math.log(10); assert.ok(p > 1.7 && p < 1.9, 'proton range exponent ' + p); }
// 6. golden output: fixed-seed runs (Node/V8 floating-point is deterministic; tolerance 1e-6 relative)
{ const gold = path.join(__dirname, 'golden/engine.json'), cases = [['alpha 5 MeV Si', ALPHA, 5, 'Si', 0.01, 400, 101], ['proton 20 MeV water', PROTON, 20, 'water', 0.6, 300, 102], ['alpha 8 MeV Au thin', ALPHA, 8, 'Au', 0.0003, 400, 103]], out = {};
  for (const [name, sp, E, m, t, N, seed] of cases) { const r = run(sp, E, MAT[m](), t, N, seed), s = stats(r); let dep = 0, path_ = 0; for (let k = 0; k < r.N; k++) { dep += r.dEdep[k]; path_ += r.pathL[k]; } out[name] = { meanZ: s.mean, sigmaZ: s.sigma, meanPath: path_ / r.N, depTotal: dep, nTrans: r.esc.nTransP, nBack: r.esc.nBackP, nStop: r.esc.nStopP, nDelta: r.lay.nDelta[0] || 0 }; }
  if (process.argv.includes('--update') || !fs.existsSync(gold)) { fs.writeFileSync(gold, JSON.stringify(out, null, 1) + '\n'); console.log('golden file written:', gold); }
  else { const ref = JSON.parse(fs.readFileSync(gold, 'utf8')); for (const k in ref) for (const q in ref[k]) { const a = out[k][q], b = ref[k][q]; assert.ok(Math.abs(a - b) <= 1e-6 * Math.max(Math.abs(b), 1e-12) + 1e-12, 'golden ' + k + ' / ' + q + ': ' + a + ' vs ' + b + '  (intended change? run with --update)'); } } }
console.log('engine: all checks passed');
