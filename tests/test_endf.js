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
console.log('ENDF parser: all checks passed');
