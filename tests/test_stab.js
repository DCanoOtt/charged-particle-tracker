// Audit of the tabulated stopping powers (STAB): catches the faults inherited from the libdedx ICRU 73 data
// (shifted rows, leaked header values, rows copied from another ion/target) so that they cannot come back.
const fs = require('fs'), path = require('path'), assert = require('assert');
const S = new Function(fs.readFileSync(process.env.STAB_FILE || path.join(__dirname, '..', 'data', 'stopping-tables.js'), 'utf8') + ';return STAB;')();
const L = Math.log, med = a => { a = [...a].sort((x, y) => x - y); return a[a.length >> 1]; };
let fails = 0; const bad = m => { fails++; console.log('FAIL', m); };
// 1. every array has the length of its energy grid, finite positive values, and no outlier in log-log (leaked first/last values)
function smooth(name, E, arr) {
  for (const k in arr) {
    const v = arr[k];
    if (!Array.isArray(v) || v.length !== E.length) { bad(`${name}/${k}: length ${v && v.length} != ${E.length}`); continue; }
    if (!v.every(x => isFinite(x) && x > 0)) { bad(`${name}/${k}: non-positive value`); continue; }
    for (let i = 0; i < v.length; i++) {
      const a = i === 0 ? 1 : i === v.length - 1 ? i - 2 : i - 1, b = i === 0 ? 2 : i === v.length - 1 ? i - 1 : i + 1;
      const t = (L(E[i]) - L(E[a])) / (L(E[b]) - L(E[a])), d = Math.abs(L(v[i]) - (L(v[a]) + t * (L(v[b]) - L(v[a]))));
      if (d > 0.12) { bad(`${name}/${k}: value ${v[i]} at E=${E[i]} deviates ${d.toFixed(2)} from its neighbours`); break; }
    }
  }
}
smooth('p', S.pE, S.p); smooth('a', S.aE, S.a); smooth('p90', S.pE90, S.p90); smooth('a90', S.aE90, S.a90); smooth('c90', S.cE90, S.c90); smooth('e90', S.eE90, S.e90);
for (const Z in S.h) smooth('h' + Z, S.hE, S.h[Z]); for (const Z in S.hn) smooth('hn' + Z, S.hE, S.hn[Z]);
// 2. no ion/target table is a copy of another ion's table for the same target
const ids = [...new Set(Object.values(S.h).flatMap(o => Object.keys(o)))];
for (const id of ids) { const seen = {}; for (const Z in S.h) { const v = S.h[Z][id]; if (!v) continue; const k = JSON.stringify(v); if (seen[k]) bad(`h: ion ${Z} / target ${id} duplicates ion ${seen[k]}`); else seen[k] = Z; } }
// 3. no row shifted by one energy bin: compare each row, relative to the silicon row of the same ion, with the shape found for Z1 = 5..8
const ref = {}; for (const id of ids) { ref[id] = []; for (let i = 0; i < S.hE.length; i++) { const a = []; for (const Z of [5, 6, 7, 8]) { const v = S.h[Z][id], w = S.h[Z][14]; if (v && w) a.push(L(v[i] / w[i])); } ref[id][i] = med(a); } }
for (const id of ids) if (id !== '14') for (let Z = 3; Z <= 18; Z++) {
  const v = S.h[Z] && S.h[Z][id], w = S.h[Z] && S.h[Z][14]; if (!v || !w) continue;
  let e0 = 0, e2 = 0; for (let i = 2; i <= 50; i++) { e0 += Math.abs(L(v[i] / w[i]) - ref[id][i]); e2 += Math.abs(L(v[i + 1] / w[i]) - ref[id][i]); }
  if (e2 < 0.5 * e0) bad(`h: ion ${Z} / target ${id} looks shifted by one energy bin`);
}
// 4. Bethe-scaled S/Z1^2 at 40-1000 MeV/u varies smoothly with Z1 (catches rows that belong to another ion)
for (let Z = 4; Z <= 17; Z++) for (const id of ids) {
  const v = S.h[Z] && S.h[Z][id], a = S.h[Z - 1] && S.h[Z - 1][id], b = S.h[Z + 1] && S.h[Z + 1][id]; if (!v || !a || !b) continue;
  let e = 0; for (let i = 40; i <= 52; i++) e += Math.abs(L((v[i] / (Z * Z)) / Math.sqrt((a[i] / ((Z - 1) ** 2)) * (b[i] / ((Z + 1) ** 2)))));
  if (e / 13 > 0.03) bad(`h: ion ${Z} / target ${id} inconsistent with its neighbours in Z1 (${(100 * e / 13).toFixed(1)} %)`);
}
assert.strictEqual(fails, 0, fails + ' stopping-table checks failed');
console.log('stopping tables: all checks passed');
