// node tests/test_data.js — audits of the reference data: elements, densities, half-lives, ICRP decay table (like test_stab.js for the stopping tables)
const fs = require('fs'), vm = require('vm'), path = require('path'), assert = require('assert');
const rd = f => fs.readFileSync(path.join(__dirname, '..', f), 'utf8');
vm.runInThisContext(rd('data/elements.js') + rd('data/densities.js') + (rd('src/engine/10-constants.js').match(/^(?:export )?const (?:SYMS|IVAL) =[\s\S]*?;\n/gm) || []).map(l => l.replace(/^export /, '')).join('\n') + rd('src/ui/48-half-lives.js') + rd('src/ui/48-decay-icrp.js') + ';globalThis.D={ELTAB,ELDENS,IVAL,HL_T,HL_MODE,HL_STABLE,DK_ICRP,DK_ICRP_STABLE};');
const { ELTAB, ELDENS, IVAL, HL_T, HL_MODE, HL_STABLE, DK_ICRP, DK_ICRP_STABLE } = D;
const near = (a, b, tol, msg) => assert.ok(Math.abs(a / b - 1) < tol, msg + ': ' + a + ' vs ' + b);
const zs = Object.keys(ELTAB).map(Number), mean = z => ELTAB[z].iso.reduce((a, i) => a + i[1] * i[2], 0) / ELTAB[z].iso.reduce((a, i) => a + i[1], 0);
// ---- elements
assert.strictEqual(new Set(zs.map(z => ELTAB[z].s)).size, zs.length, 'unique symbols');
for (const z of zs) {
  const e = ELTAB[z], A = e.iso.map(i => i[0]); assert.ok(e.iso.length > 0, e.s + ': isotopes');
  assert.ok(A.every((a, k) => k === 0 || a > A[k - 1]), e.s + ': mass numbers increase'); assert.ok(A.every(a => a >= z), e.s + ': A ≥ Z');
  assert.ok(Math.abs(e.iso.reduce((a, i) => a + i[1], 0) - 100) < 0.1, e.s + ': abundances add up to 100 % — ' + e.iso.reduce((a, i) => a + i[1], 0));
  for (const [a, ab, m] of e.iso) { assert.ok(ab >= 0 && ab <= 100, e.s + a + ': abundance'); assert.ok(Math.abs(m - a) < Math.max(0.3, 0.0012 * a), e.s + '-' + a + ': atomic mass ' + m + ' is not close to A'); }
  assert.ok(IVAL[z] > 0 && IVAL[z] < 1000, e.s + ': mean excitation energy');
  if (z >= 3) { const r = IVAL[z] / z; assert.ok(r > 8 && r < 17, e.s + ': I/Z = ' + r + ' eV is outside 8–17 eV'); if (IVAL[z - 1] && z > 3) assert.ok(IVAL[z] / IVAL[z - 1] > 0.8 && IVAL[z] / IVAL[z - 1] < 1.7, e.s + ': I jumps from the previous element'); }
}
const inv = new Set(['Ar-K', 'Co-Ni', 'Te-I', 'Th-Pa', 'Pa-U']); for (let k = 1; k < zs.length; k++) { const a = zs[k - 1], b = zs[k]; if (b === a + 1 && mean(b) < mean(a)) assert.ok(inv.has(ELTAB[a].s + '-' + ELTAB[b].s), 'mean atomic mass decreases at ' + ELTAB[a].s + '-' + ELTAB[b].s); }
for (const [s, m] of [['Fe', 55.845], ['Au', 196.967], ['Pb', 207.2], ['U', 238.029], ['Si', 28.0855], ['O', 15.999], ['Cu', 63.546]]) near(mean(zs.find(z => ELTAB[z].s === s)), m, 0.0006, 'mean atomic mass of ' + s);
// ---- densities (g/cm³; gases at 20 °C, 1 atm)
for (const z of zs) { const s = ELTAB[z].s; assert.ok(s in ELDENS, s + ': density missing'); assert.ok(ELDENS[s] > 5e-5 && ELDENS[s] < 23, s + ': density ' + ELDENS[s]); }
assert.ok(Object.keys(ELDENS).every(s => zs.some(z => ELTAB[z].s === s)), 'every density belongs to an element');
for (const [s, d] of [['Al', 2.699], ['Fe', 7.874], ['Cu', 8.96], ['Au', 19.32], ['Pb', 11.34], ['W', 19.25], ['Si', 2.329], ['Os', 22.59], ['Ge', 5.323], ['Ag', 10.49]]) near(ELDENS[s], d, 0.01, 'density of ' + s);
for (const s of ['H', 'He', 'N', 'O', 'F', 'Ne', 'Cl', 'Ar', 'Kr', 'Xe']) assert.ok(ELDENS[s] < 0.01, s + ' is a gas');
// ---- half-lives of the internal table (ids ZA = 1000 Z + A; seconds)
const ids = Object.keys(HL_T).map(Number); assert.ok(ids.length > 2500);
for (const id of ids) { const Z = Math.floor(id / 1000), A = id % 1000; assert.ok(HL_T[id] > 0 && A >= Z && Z >= 1, 'half-life entry ' + id); assert.ok(!HL_STABLE.has(id), id + ' is both stable and has a half-life'); }
{ const no = Object.keys(HL_MODE).filter(id => !(id in HL_T)); assert.ok(no.length <= 30, 'decay modes without half-life (unbound nuclei without a measured T½): ' + no.length); }
const Y = 3.15576e7;
for (const [id, T, tol] of [[27060, 5.2711 * Y, 0.002], [55137, 30.04 * Y, 0.01], [11022, 2.6029 * Y, 0.002], [19040, 1.248e9 * Y, 0.01], [92238, 4.468e9 * Y, 0.01], [88226, 1600 * Y, 0.01], [53131, 8.0252 * 86400, 0.001], [26055, 2.744 * Y, 0.005], [38090, 28.8 * Y, 0.01], [95241, 432.6 * Y, 0.005], [94239, 24110 * Y, 0.005], [43097, 4.21e6 * Y, 0.35]]) near(HL_T[id], T, tol, 'half-life of ' + id);
// ---- ICRP table: structure and consistency of the branches
const modes = new Set(); let nbr = 0;
for (const k in DK_ICRP) {
  const id = +k, za = Math.floor(id / 10), Z = Math.floor(za / 1000), A = za % 1000, [T, br] = DK_ICRP[k]; assert.ok(T > 0, k + ': half-life'); let sum = 0;
  for (const b of br.split(';')) { if (!b) continue; const [d, f, m] = b.split(':'), fr = +f, dza = Math.floor(+d / 10), dz = Math.floor(dza / 1000), da = dza % 1000; modes.add(m); sum += fr; nbr++; assert.ok(fr > 0 && fr <= 1.0000001, k + ': fraction ' + f);
    if (+d === 0) { assert.strictEqual(m, 'SF', k + ': daughter 0 is spontaneous fission'); continue; }
    const want = { 'B-': [Z + 1, A], 'B+': [Z - 1, A], EC: [Z - 1, A], A: [Z - 2, A - 4], IT: [Z, A] }[m]; assert.ok(want, k + ': unknown mode ' + m); assert.deepStrictEqual([dz, da], want, k + ' ' + m + ' → ' + d + ': daughter Z, A'); }
  assert.ok(sum > 0.96 && sum < 1.02, k + ': branches add up to ' + sum);          // ICRP omits some minor (cluster) branches: 5 nuclides add up to 0.97–0.99
  if (k % 10 === 0) assert.ok(!DK_ICRP_STABLE.has(id), k + ' is radioactive and stable');
}
assert.deepStrictEqual(Array.from(modes).sort(), ['A', 'B-', 'EC', 'IT', 'SF']);
// Bi-209 and similar are 'stable' in ICRP but have T½ > 10⁹ y in the internal table: allowed
for (const s of DK_ICRP_STABLE) { const za = Math.floor(s / 10); assert.ok(!(s in DK_ICRP), 'stable but decaying: ' + s); assert.ok(s % 10 !== 0 || !(za in HL_T) || HL_T[za] > 3e16, 'stable in ICRP but with a short half-life in the internal table: ' + za); }
// ---- internal table against the ICRP table (ground states): the two sets come from different evaluations; differences up to a factor 1.5 are known (evaluation updates);
//      larger ones are listed here — a new one means that a half-life of an isomer was taken for the ground state (as happened for Pr-134, Pm-154, Ta-178, Es-256)
const KNOWN = new Set([23050, 26060, 43097, 55135, 62146, 65146, 65154, 85215, 86216]);          // V-50, Fe-60, Tc-97, Cs-135, Sm-146, Tb-146, Tb-154, At-215, Rn-216: different evaluations
const odd = []; for (const k in DK_ICRP) { if (k % 10) continue; const za = Math.floor(k / 10), t = HL_T[za], t2 = DK_ICRP[k][0]; if (t && isFinite(t2) && (t / t2 > 1.5 || t / t2 < 1 / 1.5) && !KNOWN.has(za)) odd.push(za); }
assert.deepStrictEqual(odd, [], 'half-lives of the internal table that differ from the ICRP ground state by more than 1.5×: ' + odd.join(', '));
console.log('reference data: all checks passed');
