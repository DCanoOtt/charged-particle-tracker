// Neutron detector: cross sections, kinematics, ENDF reader, Monte Carlo efficiency checks
const fs = require('fs'), vm = require('vm'), path = require('path'), assert = require('assert');
const ui = f => fs.readFileSync(path.join(__dirname, '../src/ui', f), 'utf8');
vm.runInThisContext([ui('43-endf.js'), ui('46-secondaries.js'), ui('51-nd-data.js'), ui('52-nd-engine.js')].join('\n') +
  ';globalThis.A={ndTof,ndEfromT,ndNpXs,ndSetup,ndEff,ndKin,ND_MAT,ND_LP,ndInvLp,ndLp,ndParseN,ndChannelsC,ndMedium,ndSampleMu,ndRng,ndPulse,ndAtoms,ndEmb,tab1Eval,d2Eval,ndSmear};');
const A = globalThis.A;
const near = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, msg + ': ' + a + ' vs ' + b);
// n-p cross section (ENDF/B-VIII.0: 4.26 b at 1 MeV, 0.689 b at 14 MeV)
near(A.ndNpXs(1), 4.26, 0.1, 'sigma_np 1 MeV'); near(A.ndNpXs(14.1), 0.689, 0.03, 'sigma_np 14 MeV');
// atom densities of the scintillators
const bc = A.ND_MAT.BC501A, at = A.ndAtoms(bc.rho, bc.hc);
near(at.nH / 4.82e22, 1, 0.01, 'BC-501A hydrogen density'); near(at.nC / 3.98e22, 1, 0.01, 'BC-501A carbon density');
// kinematics: elastic scattering on a proton, energy conserved, lab cosine
const k = A.ndKin(2, 0.99862, 0, 0.3); near(k.En, 2 * (1 + 0.3 * 0.99862 * 0 + 0.3) / 2 * 1, 0.01, 'E\''); near(k.En + k.Tres, 2, 1e-9, 'energy balance (elastic: T_res = E - E\')');
const ki = A.ndKin(10, 11.8969, -4.439, 0); near(ki.En + ki.Tres + 4.439, 10, 1e-9, 'inelastic energy balance'); assert.strictEqual(A.ndKin(3, 11.8969, -4.439, 0), null, 'below threshold');
// light curves: Cecil NE-213 proton light at 1 MeV about 0.2 MeVee, 100 keVee threshold at 0.63 MeV protons
near(A.ndLp(1, A.ND_LP.liq), 0.206, 0.01, 'L_p(1 MeV)'); near(A.ndLp(A.ndInvLp(0.1, A.ND_LP.liq), A.ND_LP.liq), 0.1, 1e-6, 'inverse light curve');
// thin detector, parallel beam: Monte Carlo = single-collision analytic estimate below the carbon inelastic threshold
const mk = (D, T, Ld, th, data) => A.ndSetup({ mat: 'BC501A', D, T, th, rho: bc.rho, hc: bc.hc, lp: A.ND_LP.liq, Lm: Ld / 100 }, data || { H: null, C: null });
let S = mk(5, 0.05, 1e6, 100);
for (const E of [2, 5]) { const r = A.ndEff(S, E, 600000, 7); near(r.eps / r.eps1, 1, 0.04, 'thin cell, E=' + E); }
// thick cell: efficiency decreases with energy, multiple scattering raises it above the single-collision value
S = mk(12.7, 5.08, 2500, 100);
const e1 = A.ndEff(S, 2.5, 30000, 3), e2 = A.ndEff(S, 14, 30000, 3);
assert.ok(e1.eps > 0.3 && e1.eps < 0.5 && e2.eps > 0.12 && e2.eps < 0.26 && e1.eps > e2.eps && e1.eps > e1.eps1, 'thick 5"x2" cell: ' + e1.eps + ' ' + e2.eps);
assert.strictEqual(A.ndEff(S, 0.4, 5000, 3).eps, 0, 'below the threshold energy nothing is detected');
// higher threshold lowers the efficiency
near(A.ndEff(mk(12.7, 5.08, 2500, 400), 2.5, 30000, 3).eps < e1.eps ? 1 : 0, 1, 0, 'threshold');
// ENDF/B-VIII.0 layout files: n+1H (MF4 Legendre) and n+12C (MF3 channels, MF4 tabulated)
const H = A.ndParseN(fs.readFileSync(path.join(__dirname, 'data/n_001-H-001.endf'), 'utf8')), C = A.ndParseN(fs.readFileSync(path.join(__dirname, 'data/n_006-C-012.endf'), 'utf8'));
assert.strictEqual(H.ZA, 1001); assert.strictEqual(C.ZA, 6012); assert.ok(H.mts[2] && H.ang && H.ang.ltt === 1 && H.ang.lct === 2, 'H-1 elastic + angular distribution');
const ch = A.ndChannelsC(C), kinds = Object.fromEntries(ch.map(c => [c.mt, c.kind + (c.par ? ':' + c.par : '')]));
assert.deepStrictEqual(kinds, { 51: 'inel:a', 52: 'brk:a', 91: 'brk:a', 102: 'abs', 103: 'chg:p', 107: 'chg:a' }, JSON.stringify(kinds));
assert.ok(C.ang && C.ang.ltt === 2 && C.ang.L.length === 2, 'C-12 tabulated angular distribution');
near(H.mts[2].xs.y[H.mts[2].xs.y.length >> 1] > 0 ? 1 : 0, 1, 0, 'cross-section table');
const med = A.ndMedium(bc.rho, bc.hc, { H, C });
assert.ok(med.fromFile.H && med.fromFile.C && med.fromFile.angH && med.fromFile.angC, 'medium built from the files');
// forward-peaked tabulated C-12 distribution at 40 MeV: mean cosine > 0; isotropic at 1 MeV
const rng = A.ndRng(5); let m1 = 0, m40 = 0; for (let i = 0; i < 20000; i++) { m1 += A.ndSampleMu(med.angC, 1, rng()); m40 += A.ndSampleMu(med.angC, 40, rng()); }
near(m1 / 20000, 0, 0.02, 'isotropic at 1 MeV'); assert.ok(m40 / 20000 > 0.3, 'forward peaked at 40 MeV: ' + m40 / 20000);
// H-1 P1 coefficient a1 = 0.2 at 40 MeV: <mu> = a1 = 0.2
let mh = 0; for (let i = 0; i < 20000; i++) mh += A.ndSampleMu(med.angH, 40, rng()); near(mh / 20000, 0.2, 0.02, 'Legendre P1 sampled');
// efficiency with the file data is similar to the built-in values
const Sf = mk(12.7, 5.08, 2500, 100, { H, C }), ef = A.ndEff(Sf, 2.5, 30000, 3); near(ef.eps / e1.eps, 1, 0.1, 'file data vs built-in at 2.5 MeV');
// pulse-height spectrum: all weights end up in the histogram, none above the maximum light
const Ec = Float64Array.from({ length: 50 }, (_, i) => 1 + i * 0.1), cdf = Float64Array.from(Ec, (_, i) => (i + 1) / 50), P = A.ndPulse(S, Ec, cdf, 20000, 100, 8, 9);
assert.ok(P.sw > 19000 && P.all.reduce((a, b) => a + b, 0) / P.sw > 0.2, 'pulse-height spectrum filled');
// time of flight: 25 m, 14 MeV = 488 ns; round trip; the reconstructed energy of a detected event without time smearing is its source energy
near(A.ndTof(14, 25), 488.4, 0.2, 'TOF 14 MeV'); near(A.ndEfromT(A.ndTof(3.3, 25), 25), 3.3, 1e-9, 'E(t(E))');
const rec = { nE: 100, Eo: 10, nT: 50, t0: 500, t1: 2000, Lm: 25, st: 0 }, Ec2 = Float64Array.of(2.05), P2 = A.ndPulse(S, Ec2, Float64Array.of(1), 20000, 50, 3, 4, rec);
const be = P2.recE.findIndex(v => v > 0); assert.ok(be === 20 && P2.recE.every((v, i) => i === 20 || v === 0), 'rec bin ' + be);
assert.ok(P2.recE[20] > 0 && Math.abs(P2.recE[20] / A.ndEff(S, 2.05, 20000, 4).eps / P2.sw - 1) < 0.1, 'reconstructed counts = efficiency x neutrons');
const rec2 = { ...rec, Lm: 1, st: 3 / 2.3548 }, P3 = A.ndPulse(S, Ec2, Float64Array.of(1), 20000, 50, 3, 4, rec2);
assert.ok(P3.recE.filter(v => v > 0).length > 3, 'time smearing broadens the reconstructed energy');
console.log('neutron detector: all checks passed');

// embedded ENDF/B-VIII.0 ACE data at 293.6 K (1001.800nc, 6012.800nc)
{
  const H = A.ndEmb('H'), C = A.ndEmb('C');
  near(A.tab1Eval(H.mts[2].xs, 1), 4.2495, 0.01, 'embedded sigma_el(n,1H) 1 MeV'); near(A.tab1Eval(H.mts[2].xs, 14.1), 0.683, 0.005, 'embedded sigma_el(n,1H) 14.1 MeV');
  near(A.tab1Eval(C.mts[2].xs, 14.1), 0.827, 0.01, 'embedded sigma_el(n,12C) 14.1 MeV'); assert.ok(Object.keys(C.mts).length >= 15, 'C-12 channels');
  assert.ok(A.d2Eval(C.ang, 14.1, 1) > 3 * A.d2Eval(C.ang, 14.1, -1), 'C-12 elastic forward peaked at 14 MeV'); near(A.d2Eval(H.ang, 5, 0.9) / A.d2Eval(H.ang, 5, -0.9), 1, 0.1, 'H-1 almost isotropic in CM');
  const med = A.ndMedium(0.874, 1.212, { H, C }); assert.ok(med.fromFile.H && med.fromFile.C && med.fromFile.angH && med.fromFile.angC, 'medium uses the embedded files');
  console.log('embedded ACE data ok');
}

// light-output resolution: FWHM/L = 20 % at 1 MeVee, scaling as 1/sqrt(L); 0 = ideal
{
  const rng = A.ndRng(5);
  for (const [L, rel] of [[1, 0.2], [4, 0.1], [0.25, 0.4]]) {
    let m = 0, m2 = 0; const n = 40000;
    for (let i = 0; i < n; i++) { const v = A.ndSmear(L, 0.2, rng); m += v; m2 += v * v }
    m /= n; const sd = Math.sqrt(m2 / n - m * m);
    near(m / L, 1, 0.01, 'mean light unchanged at ' + L + ' MeVee'); near(sd / L * 2.3548, rel, 0.01 * rel + 0.004, 'FWHM/L at ' + L + ' MeVee');
  }
  assert.strictEqual(A.ndSmear(2, 0, rng), 2, 'resolution 0 leaves the light unchanged');
  assert.strictEqual(A.ndSetup({ mat: 'BC501A', D: 10, T: 10, th: 100, res: 20, rho: 0.874, hc: 1.212, lp: {}, Lm: 25 }, null).res, 0.2);
  console.log('light resolution ok');
}
