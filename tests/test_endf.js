// node tests/test_endf.js  — parser test on the synthetic ENDF-6 file (python3 tests/make_endf.py writes it)
const fs = require('fs'), vm = require('vm'), path = require('path'), assert = require('assert');
vm.runInThisContext(fs.readFileSync(path.join(__dirname, '../src/ui/43-endf.js'), 'utf8') + ';globalThis.API={parseEndfAlpha,chYield,chSigma,endfNum,PROJ};');
const { parseEndfAlpha, chYield, chSigma, endfNum, PROJ } = API;
assert.strictEqual(endfNum(' 8.921400+4'), 89214); assert.ok(Math.abs(endfNum('-1.234567-5') + 1.234567e-5) < 1e-12); assert.strictEqual(endfNum('        283'), 283);
const r = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/a-Ar040.tendl'), 'utf8'));
assert.strictEqual(r.Z, 18); assert.strictEqual(r.A, 40); assert.deepStrictEqual(r.order, [2, 3, 5, 16, 102, 103, 600, 649]);
assert.ok(Math.abs(chSigma(r.channels[3], 10) - 0.7886) < 1e-3, 'MT3 value');
const c16 = r.channels[16]; assert.strictEqual(chYield(c16, 'n', 10), 2); assert.ok(Math.abs(chYield(c16, 'g', 10) - 1.6667) < 1e-3, 'photon yield interpolated from MF6');
const c5 = r.channels[5]; assert.strictEqual(chYield(c5, 'n', 5), 0.5); assert.strictEqual(chYield(c5, 'p', 5), 0.3); assert.strictEqual(chYield(c5, 'a', 5), 0.2);
assert.strictEqual(chYield(r.channels[103], 'p', 5), 1); assert.strictEqual(r.channels[102].res[0].ZA, 20044);
assert.strictEqual(c5.prods.length, 3); assert.strictEqual(c5.prods[0].dist.law, 1); assert.strictEqual(c5.prods[0].dist.lang, 2); assert.strictEqual(c5.prods[0].dist.frame, 'lab');
assert.strictEqual(c16.prods[0].dist.lang, 1); assert.strictEqual(c16.prods[0].dist.frame, 'cm'); assert.strictEqual(r.channels[600].prods[0].dist.law, 2);
const rp = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/a-Ar040.tendl'), 'utf8'), PROJ.p);
assert.strictEqual(rp.channels[102].res[0].ZA, 19041, 'MF6-less fallback: compound nucleus follows the projectile (p + Ar-40 -> K-41, not Ca-44)');
const rd = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/a-Ar040.tendl'), 'utf8'), PROJ.d);
assert.strictEqual(rd.channels[102].res[0].ZA, 19042, 'd + Ar-40 -> K-42');
// JENDL-5 style file: MF3 MT2/MT5 only, MF6/MT5 with LAW 7 (laboratory angle-energy) products and bare-yield residuals, incl. the unbound H-4 (ZAP 1004)
const rj = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/d_003-Li-006.dat'), 'utf8'), PROJ.d);
assert.strictEqual(rj.Z, 3); assert.strictEqual(rj.A, 6); assert.deepStrictEqual(rj.order, [2, 5]);
const j5 = rj.channels[5]; assert.strictEqual(j5.name, '(d,anything)'); assert.ok(j5.primary && !rj.channels[2].primary);
assert.strictEqual(j5.prods.length, 7, 'all products read although the first ones have LAW 7 (the neutron is listed three times)');
assert.ok(Math.abs(chYield(j5, 'n', 5) - 0.5) < 1e-9 && Math.abs(chYield(j5, 'a', 5) - 0.4) < 1e-9 && Math.abs(chYield(j5, 'p', 5) - 0.3) < 1e-9, 'yields of the LAW 7 products');
assert.deepStrictEqual(j5.res.map(r => r.ZA), [3007, 1004], 'residual nuclei, H-4 included');
assert.ok(Math.abs(chSigma(j5, 5) - 0.98) < 1e-9);
// ENDF/B-VIII.0 d + Li-6 (file d-003_Li_006.endf): MF3/MF6 for MT 2, 50 (d,n0), 600 (d,p0), 800 (d,α0), no MT 3/5; evaluation ends at 5 MeV
const re = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/d-003_Li_006.endf'), 'utf8'), PROJ.d);
assert.deepStrictEqual(re.order, [2, 50, 600, 800]); assert.strictEqual(re.channels[50].name, '(d,n) to level 0'); assert.ok(re.channels[50].primary && re.channels[600].primary && re.channels[800].primary && !re.channels[2].primary);
assert.strictEqual(chYield(re.channels[50], 'n', 3), 1); assert.strictEqual(chYield(re.channels[800], 'a', 3), 2, 'two α per (d,α) reaction'); assert.strictEqual(re.channels[50].res[0].ZA, 4007, 'Be-7 residual of (d,n)');
assert.ok(Math.abs(chSigma(re.channels[600], 1) - 0.0675) < 1e-3 && chSigma(re.channels[600], 6) === 0, 'cross sections are zero above the 5 MeV of the evaluation');
// LAW 7 data are read (not skipped) and evaluated in the secondaries module: normalisation ∫∫ f dE' dμ = 1, isotropic here
vm.runInThisContext(fs.readFileSync(path.join(__dirname, '../src/ui/46-secondaries.js'), 'utf8') + ';globalThis.API7={d7Tot,d7Eval};');
const dn = j5.prods[0].dist; assert.ok(dn && dn.law === 7 && dn.frame === 'lab' && dn.L.length === 2 && dn.L[0].mu.length === 2, 'LAW 7 distribution of the neutron');
let I7 = 0, I7a = 0; const E7 = 1, top = 1.05, N7 = 4000;
for (let i = 0; i < N7; i++) { const e = (i + 0.5) * top / N7; I7 += API7.d7Tot(dn, E7, e) * top / N7; I7a += API7.d7Eval(dn, E7, 0.3, e) * top / N7 }
assert.ok(Math.abs(I7 - 1) < 2e-3, 'energy distribution integrates to 1: ' + I7); assert.ok(Math.abs(I7a - 0.5) < 2e-3, 'f(E\', mu) integrates to 1/2 per unit mu: ' + I7a);
// JENDL-5 alpha file: neutron production MT201 (MF3 σ, MF6 yield 1, LAW 7) is merged into MT5 as the multiplicity σ201/σ5 (MT5 itself has no neutrons)
const ra = parseEndfAlpha(fs.readFileSync(path.join(__dirname, 'data/a_013-Al-027.dat'), 'utf8'), PROJ.a);
assert.deepStrictEqual(ra.order, [2, 5, 201]); const a5 = ra.channels[5];
assert.ok(ra.channels[201].merged && !ra.channels[201].primary && ra.channels[201].name === '(α,xn) neutron production');
assert.ok(Math.abs(chYield(a5, 'n', 10) - 0.27 / 0.9) < 1e-9 && chYield(a5, 'n', 3) === 0, 'neutron multiplicity = σ201 / σ5 (zero at the threshold of MT201)');
assert.ok(a5.prods.some(p => p.sp === 'n' && p.dist && p.dist.law === 7), 'the neutron distribution travels with the merged product');
assert.ok(Math.abs(chSigma(a5, 10) * chYield(a5, 'n', 10) - 0.27) < 1e-9, 'σ5 × multiplicity = the neutron production cross section');
console.log('ENDF parser: all checks passed');
