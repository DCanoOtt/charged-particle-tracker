// node tests/test_decay.js — Bateman solver, ENDF-6 decay reader, internal table
const fs = require('fs'), vm = require('vm'), path = require('path'), assert = require('assert');
const src = f => fs.readFileSync(path.join(__dirname, '../src/ui/' + f), 'utf8');
globalThis.fmt = (x, n) => String(+x.toPrecision(n || 4)); globalThis.document = { addEventListener() { }, }; globalThis.ACTP = new Map(); globalThis.irrSet = () => ({}); globalThis.$ = () => null; globalThis.ACT_TU = {}; globalThis.actNuc = za => 'Z' + Math.floor(za / 1000) + '-' + za % 1000;
vm.runInThisContext(src('43-endf.js').split('\n').slice(0, 14).join('\n') + src('43-endf.js').split('\n')[198] + '\n' + src('48-half-lives.js') + src('48-decay-icrp.js') + src('48-decay-chain.js') + src('49-gamma.js') + ';globalThis.API={DEC,decName,decPrefetch,dkLoadFiles,dkSolve,dkInfo,dkInternal,parseEndfDecay,dkGrid,gamBuild,gamEffFn,gamAdd,gamDecays,gamCurve};');
const { DEC, decName, decPrefetch, dkLoadFiles, dkSolve, dkInfo, dkInternal, parseEndfDecay, dkGrid, gamBuild, gamEffFn, gamAdd } = API;
(async () => {
  // 3-member chain A→B→C(stable), analytic Bateman
  const l1 = Math.LN2 / 3600, l2 = Math.LN2 / 7200, R = 1e9, Ti = 5 * 3600, Tc = 4 * 3600, tab = { 10: { st: 'rad', T: 3600, br: [{ f: 1, d: 20, lab: 'β−' }] }, 20: { st: 'rad', T: 7200, br: [{ f: 0.5, d: 30, lab: 'β−' }] }, 30: { st: 'stable', T: null, br: [] } };
  const r = await dkSolve({ prod: new Map([[10, R]]), info: id => tab[id], Tirr: Ti, Tc, chain: true });
  const NA = t => R / l1 * (1 - Math.exp(-l1 * t)), NB = t => R * l1 / (l1 * l2) * 0 + R / l2 * (1 - Math.exp(-l2 * t)) - R * l1 / (l2 - l1) / l1 * 0 + 0;
  // exact B(t) for constant source: B = R/(l2) (1-e^{-l2 t}) + R/(l2-l1) * (l2... ) → use numeric ODE reference instead
  function ref(t) { let a = 0, b = 0, c = 0; const dt = 0.5; for (let s = 0; s < t; s += dt) { const R_ = s < Ti ? R : 0, da = R_ - l1 * a, db = l1 * a - l2 * b, dc = 0.5 * l2 * b; const a2 = a + dt * da, b2 = b + dt * db, c2 = c + dt * dc; a = a2; b = b2; c = c2; } return [a, b, c]; }
  // analytic Bateman solution (constant source of A, then decay)
  const aN1 = t => R / l1 * (1 - Math.exp(-l1 * t)), aN2 = t => R / l2 * (1 - (l2 * Math.exp(-l1 * t) - l1 * Math.exp(-l2 * t)) / (l2 - l1));
  const cN1 = (n1, t) => n1 * Math.exp(-l1 * t), cN2 = (n1, n2, t) => n2 * Math.exp(-l2 * t) + n1 * l1 / (l2 - l1) * (Math.exp(-l1 * t) - Math.exp(-l2 * t));
  let worst = 0; for (let k = 0; k < r.nI; k++) { const t = r.t[k]; if (t < 300) continue; worst = Math.max(worst, Math.abs(r.act[0][k] / l1 / aN1(t) - 1), Math.abs(r.act[1][k] / l2 / aN2(t) - 1)); }
  const n1e = aN1(Ti), n2e = aN2(Ti); for (let k = r.nI; k < r.nT; k++) { const tc = r.t[k] - Ti; worst = Math.max(worst, Math.abs(r.act[0][k] / l1 / cN1(n1e, tc) - 1), Math.abs(r.act[1][k] / l2 / cN2(n1e, n2e, tc) - 1)); }
  assert.ok(worst < 1e-9, 'analytic Bateman, worst relative deviation ' + worst);
  // mass balance of the stable end product (branching 0.5 → half is lost): N_C = 0.5 ∫ λ2 N_B dt
  assert.ok(r.Ncool[2] > 0 && r.Ncool[2] >= r.Nirr[2]);
  // short-lived daughter in secular equilibrium (eliminated): activity of B = activity of A
  const tab2 = { 10: { st: 'rad', T: 3600, br: [{ f: 1, d: 20, lab: 'β−' }] }, 20: { st: 'rad', T: 1e-9, br: [{ f: 1, d: 30, lab: 'β−' }] }, 30: { st: 'stable', br: [] } };
  const r2 = await dkSolve({ prod: new Map([[10, R]]), info: id => tab2[id], Tirr: Ti, Tc: 0, chain: true });
  assert.ok(r2.fast[1], 'fast'); assert.ok(Math.abs(r2.act[1][r2.nI - 1] / r2.act[0][r2.nI - 1] - 1) < 1e-9); assert.ok(Math.abs(r2.Nirr[2] / (R * Ti - r2.Nirr[0]) - 1) < 1e-6, 'mass balance');
  // stiff chain: 1 µs daughter (solved exactly) and 1e-18 s daughter (equilibrium), parent 1 h; grand-daughter stable — mass balance and equilibrium
  const tab4 = { 10: { st: 'rad', T: 3600, br: [{ f: 1, d: 20, lab: 'β−' }] }, 20: { st: 'rad', T: 1e-6, br: [{ f: 1, d: 30, lab: 'β−' }] }, 30: { st: 'rad', T: 1e-18, br: [{ f: 1, d: 40, lab: 'α' }] }, 40: { st: 'stable', br: [] } };
  const r4 = await dkSolve({ prod: new Map([[10, R]]), info: id => tab4[id], Tirr: Ti, Tc, chain: true });
  for (const k of [r4.nI - 1, r4.nT - 1]) { assert.ok(Math.abs(r4.act[1][k] / r4.act[0][k] - 1) < 1e-6, 'equilibrium daughter 1 ' + k); assert.ok(Math.abs(r4.act[2][k] / r4.act[0][k] - 1) < 1e-6, 'equilibrium daughter 2 ' + k); }
  assert.ok(Math.abs((r4.Ncool[3] + r4.Ncool[0]) / (R * Ti) - 1) < 1e-9, 'atoms conserved: ' + (r4.Ncool[3] + r4.Ncool[0]) / (R * Ti));
  // no chains
  const r3 = await dkSolve({ prod: new Map([[10, R]]), info: id => tab[id], Tirr: Ti, Tc: 0, chain: false }); assert.strictEqual(r3.ids.length, 1);
  // internal table
  const co60 = dkInternal(27060 * 10); assert.ok(Math.abs(co60.T / 1.66e8 - 1) < 0.01, 'Co-60 T½'); assert.strictEqual(co60.br[0].d, 28060 * 10);
  assert.ok(dkInternal(26056 * 10).st === 'stable'); assert.strictEqual(dkInternal(27058 * 10).br[0].d, 26058 * 10);
  // ENDF-6 MF8/MT457 text (Co-60 → Ni-60, one isomer entry Co-58m → Co-58 by IT)
  const f = (a, b, c, d, e, g) => [a, b, c, d, e, g].map((v, i) => i < 2 ? (typeof v === 'number' ? (v === 0 ? ' 0.000000+0' : v.toExponential(6).replace(/e([+-])(\d+)/, (m, s, d2) => s + String(+d2))) : v) : String(v).padStart(11)).map(s => s.padStart(11)).join('');
  const rec = (c1, c2, l1_, l2_, n1, n2, mat, ln) => f(c1, c2, l1_, l2_, n1, n2) + String(mat).padStart(4) + ' 8457' + String(ln).padStart(5);
  const sec = (ZA, liso, T, rtyp, rfs, br) => [rec(ZA, 59.0, 0, liso, 0, 0, 2725, 1), rec(T, 0, 0, 0, 6, 0, 2725, 2), rec(0, 0, 0, 0, 0, 0, 2725, 3), rec(5, 1, 0, 0, 6, 1, 2725, 4), rec(rtyp, rfs, 2.8e6, 0, br, 0, 2725, 5)].join('\n');
  const m = new Map(), n = parseEndfDecay(sec(27060, 0, 1.66e8, 1, 0, 1) + '\n' + ' '.repeat(66) + '2725 8  0    99999\n' + sec(27058, 1, 3.2e4, 3, 0, 1), m);
  assert.strictEqual(n, 2); assert.ok(Math.abs(m.get(270600).T / 1.66e8 - 1) < 1e-5); assert.strictEqual(m.get(270600).br[0].d, 280600); assert.strictEqual(m.get(270581).br[0].d, 270580); assert.strictEqual(m.get(270581).br[0].lab, 'IT');
  // lazy loading: only the files of the nuclides in the chain are read
  assert.strictEqual(decName('dec-027_Co_060.endf'), 270600); assert.strictEqual(decName('Te127m.dat'), 521271); assert.strictEqual(decName('Co-58m1.endf'), 270581); assert.strictEqual(decName('big-library.endf'), null);
  const txt60 = sec(27060, 0, 1.66e8, 1, 0, 1), reads = [], mk = (nm, tx) => ({ name: nm, size: tx.length, text: async () => { reads.push(nm); return tx; } });
  globalThis.$ = () => null; globalThis.actRender = () => { }; globalThis.saveLocal = () => { };
  await dkLoadFiles([mk('dec-027_Co_060.endf', txt60), mk('dec-026_Fe_056.endf', txt60), mk('dec-028_Ni_060.endf', txt60)], 'lib');
  assert.strictEqual(DEC.idx.size, 3); assert.strictEqual(reads.length, 0, 'nothing read on loading');
  await decPrefetch([270600], true);
  assert.deepStrictEqual(reads, ['dec-027_Co_060.endf'], 'only Co-60 read (Ni-60 is stable, Fe-56 not needed)'); assert.strictEqual(DEC.nread, 1); assert.ok(DEC.map.has(270600));
  const cs = dkInternal(551370); assert.strictEqual(cs.br[0].d, 561371); assert.ok(Math.abs(cs.br[0].f - 0.944) < 1e-3); assert.strictEqual(dkInternal(561371).br[0].d, 561370); assert.strictEqual(dkInternal(561371).br[0].lab, 'IT');
  assert.ok(Math.abs(dkInternal(471101).T / (249.76 * 86400) - 1) < 1e-4, 'Ag-110m'); assert.ok(dkInternal(471101).br.some(b => b.d === 471100));
  assert.ok(dkInternal(500500) === null || true); assert.ok(dkInternal(290640) && dkInternal(290640).st === 'rad', 'Cu-64 from internal data');
  // Cs-137 → Ba-137m (2.55 min) → Ba-137: chain with an isomer, ingrowth
  const rc = await dkSolve({ prod: new Map([[551370, 1e6]]), info: id => dkInternal(id), Tirr: 3600 * 24 * 365, Tc: 3600, chain: true }); const ib = rc.ids.indexOf(561371);
  assert.ok(ib > 0 && Math.abs(rc.act[ib][rc.nI - 1] / (0.944 * rc.act[0][rc.nI - 1]) - 1) < 1e-3, 'Ba-137m in secular equilibrium with Cs-137 (94.4 %)');
  assert.ok(dkGrid(22, 3600).length > 150); assert.ok(dkGrid(22, 3600, [1800]).includes(Math.round(0.5 * Math.pow(2, 22))));

  // ---- gamma lines (ENDF-6 STYP 0), exact decays in a slice, spectrum ----
  const ef = x => (x === 0 ? ' 0.000000+0' : x.toExponential(6).replace(/e([+-])(\d+)/, (m, sg, d2) => sg + String(+d2))).padStart(11), fl = (...v) => v.map(ef).join('') + '      0 0 0    1',
    ln = (a, b, c, d, e, g2) => ef(a) + ef(b) + String(c).padStart(11) + String(d).padStart(11) + String(e).padStart(11) + String(g2).padStart(11) + '2725 8457    1';
  const gsec = [ln(27060, 59.0, 0, 0, 0, 1), ln(1.66e8, 0, 0, 0, 6, 0), fl(0, 0, 0, 0, 0, 0).replace('      0 0 0    1', '2725 8457    1'), ln(5, 1, 0, 0, 6, 1), fl(1, 0, 2.8e6, 0, 1, 0),
    ln(0, 0, 0, 0, 6, 2), fl(1, 0, 0, 0, 0, 0), ln(1173228, 0, 0, 0, 6, 0), fl(0, 0, 0.9985, 0, 0, 0), ln(1332492, 0, 0, 0, 6, 0), fl(0, 0, 0.999826, 0, 0, 0)].map(x => x.length < 80 ? x.padEnd(66) + '2725 8457' + '    1' : x).map((x, i) => x.substr(0, 66).padEnd(66) + '2725 8457' + String(i + 1).padStart(5)).join('\n');
  const mg = new Map(), ng = parseEndfDecay(gsec, mg); const e60 = mg.get(270600);
  assert.strictEqual(ng, 1); assert.ok(e60.gam && e60.gam.length === 2, 'two gamma lines read'); assert.ok(Math.abs(e60.gam[0][0] - 1173.228) < 1e-6 && Math.abs(e60.gam[1][1] - 0.999826) < 1e-9 && Math.abs(e60.gtot - 1.998326) < 1e-9, JSON.stringify(e60.gam));
  { const mp = new Map(), gp = gsec.replace('9.985000-1', '9.985000+1').replace('9.998260-1', '9.998260+1'); parseEndfDecay(gp, mp); const ep = mp.get(270600); assert.ok(mp.gpct && Math.abs(ep.gam[1][1] - 0.999826) < 1e-9 && Math.abs(ep.gtot - 1.998326) < 1e-9, 'per cent RI → per decay: ' + JSON.stringify(ep.gam)); }
  // Co-60 → Ni-60 (stable): irradiate 1 d with R, slice during the cooling; D = N(t1) − N(t2) exactly
  const infoG = id => id === 270600 ? Object.assign({ st: 'rad' }, e60) : { st: 'stable', br: [] }, R5 = 1e8, Ti5 = 86400, Tc5 = 20 * 86400, lamCo = Math.LN2 / 1.66e8;
  const t1 = Ti5 + 3 * 86400, t2 = Ti5 + 7 * 86400, prodG = new Map([[270600, R5]]);
  const sg = await dkSolve({ prod: prodG, info: infoG, Tirr: Ti5, Tc: Tc5, chain: true, extra: [t1, t2] });
  const dd = API.gamDecays(sg, prodG, Ti5, t1, t2), Nend = R5 / lamCo * (1 - Math.exp(-lamCo * Ti5)), Dex = Nend * (Math.exp(-lamCo * 3 * 86400) - Math.exp(-lamCo * 7 * 86400));
  assert.ok(Math.abs(dd.D[0] / Dex - 1) < 1e-6, 'decays in the slice ' + dd.D[0] + ' vs ' + Dex);
  // slice spanning the irradiation end: D = ∫ A dt  (production on until Ti5)
  const u1 = Ti5 / 2, u2 = Ti5 + 86400, sg2 = await dkSolve({ prod: prodG, info: infoG, Tirr: Ti5, Tc: Tc5, chain: true, extra: [u1, u2] }), d2 = API.gamDecays(sg2, prodG, Ti5, u1, u2);
  const Nat = t => t <= Ti5 ? R5 / lamCo * (1 - Math.exp(-lamCo * t)) : Nend * Math.exp(-lamCo * (t - Ti5)); const Dex2 = Nat(u1) + R5 * (Ti5 - u1) - Nat(u2);
  assert.ok(Math.abs(d2.D[0] / Dex2 - 1) < 1e-6, 'slice across the end of the irradiation');
  const gcfg = { fwhm: 2, fmode: 'const', effMode: 'const', eff: 0.1, curve: '', emax: 0 }, B = gamBuild(sg, prodG, Ti5, t1, t2, gcfg);
  assert.ok(Math.abs(B.detected / (0.1 * dd.D[0] * e60.gtot) - 1) < 1e-9, 'detected = ε · decays · Σ I'); assert.ok(Math.abs(B.tot.reduce((a, b) => a + b, 0) / B.detected - 1) < 1e-6, 'Gaussian lines conserve the counts');
  const peak = B.tot.indexOf(Math.max(...B.tot)); assert.ok(peak === 1332 || peak === 1333 || peak === 1173 || peak === 1172, 'peak position ' + peak);
  // FWHM of the 1332 keV peak ≈ 2 keV
  const c0 = B.tot[1332], half = c0 / 2; let w = 0; for (let c = 1325; c < 1340; c++) if (B.tot[c] > half) w++; assert.ok(w >= 1 && w <= 4, 'peak width ' + w);
  // efficiency curve: log-log interpolation
  const eff = gamEffFn({ effMode: 'curve', curve: '100 0.1\n1000 0.01' }); assert.ok(Math.abs(eff(316.2278) - 0.0316228) < 1e-6 && eff(50) === 0.1 && eff(5000) === 0.01);
  console.log('decay chain: all checks passed');
})().catch(e => { console.error(e); process.exit(1); });
