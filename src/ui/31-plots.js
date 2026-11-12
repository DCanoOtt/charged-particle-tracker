/* ====== plot machinery ====== */
let _chart = null;   // {opts, geom}
const lenSel = [['auto', 'auto'], ['nm', 'nm'], ['µm', 'µm'], ['mm', 'mm'], ['cm', 'cm']];
const normSel = [['prim', 'fraction of primaries per bin'], ['counts', 'counts per bin'], ['dens', 'probability density']];
function radialDepthCm(s) {
  if (s.depthCm > 0) return s.depthCm;
  const R = RES, zs = []; for (let i = 0; i < R.N; i++) if (R.status[i] === 0) zs.push(R.endZ[i]);
  return zs.length ? 0.5 * zs.reduce((a, b) => a + b, 0) / zs.length : 0.5 * R.edges[R.L];
}
function pset(id) { const sp = PLOTS[id]; const s = S.plots[id] = Object.assign({}, sp.defaults || {}, S.plots[id] || {}); return s; }
const layOpts = () => [['all', 'All layers']].concat(S.layers.map((l, i) => [String(i), (i + 1) + ' · ' + layerName(l, i)]));
function unitOf(s, key, cm) { return !s[key] || s[key] === 'auto' ? autoLen(cm) : s[key]; }
function binRange(s, dmin, dmax) {
  let lo = pnum(s.xmin), hi = pnum(s.xmax);
  if (lo == null) lo = dmin; if (hi == null) hi = dmax;
  if (!(hi > lo)) { const c = Math.abs(lo) || 1; lo -= 0.5 * c * 0.1; hi += 0.5 * c * 0.1; if (hi <= lo) { lo = lo - 1; hi = hi + 1; } }
  return [lo, hi];
}
function normC(counts, edges, nPrim, mode, areaFn) {
  const out = new Float64Array(counts.length);
  for (let i = 0; i < counts.length; i++) {
    const w = edges[i + 1] - edges[i];
    out[i] = mode === 'counts' ? counts[i] : mode === 'dens' ? counts[i] / (nPrim * (areaFn ? areaFn(i) : w)) : counts[i] / nPrim;
  }
  return out;
}
const ylabNorm = (mode, unit, what) => mode === 'counts' ? 'Counts per bin' : mode === 'dens' ? 'Probability density (1/' + unit + ')' : (what || 'Fraction of primaries') + ' per bin';
function statsTable(rows, unit) {
  return '<details class="cw"><summary><b>Statistics of the distribution</b></summary><div class="tscroll"><table class="t"><tr><th>Quantity</th><th>N</th><th>Mean</th><th>σ</th><th>FWHM (2.355σ)</th><th>Median</th><th>Min</th><th>Max</th></tr>' + rows.map(r => r.d ? '<tr><td>' + esc(r.name) + '</td><td>' + r.d.n.toLocaleString('en') + '</td><td>' + fmt(r.d.mean, 5) + ' ± ' + fmt(r.d.sem, 2) + '</td><td>' + fmt(r.d.sd, 4) + '</td><td>' + fmt(r.d.fwhm, 4) + '</td><td>' + fmt(r.d.med, 5) + '</td><td>' + fmt(r.d.min, 5) + '</td><td>' + fmt(r.d.max, 5) + '</td></tr>' : '<tr><td>' + esc(r.name) + '</td><td colspan="7">no entries</td></tr>').join('') + '</table></div><p class="hint">All values in ' + esc(unit) + '. Mean ± statistical error of the mean.</p></details>';
}
function crossings(R, d) {
  // First forward crossing of the plane z = d by every stored primary track. Vertices are stored as float32, so the plane is
  // compared in float32 too (a layer boundary and the point just beyond it can round to the same value).
  const pv = R.pv, off = R.pvOff, ax = R.d0[0] / R.d0[2], ay = R.d0[1] / R.d0[2], X = [], Y = [], E = [], dd = Math.fround(d); let tracked = 0;
  for (let k = 0; k < R.N; k++) {
    const a = off[k], b = off[k + 1]; if (b - a < 1) continue; tracked++;
    const sx = pv[a * 4], sy = pv[a * 4 + 1];
    if (dd <= 0) { X.push(0); Y.push(0); E.push(pv[a * 4 + 3]); continue; }
    for (let j = a + 1; j < b; j++) {
      const z0 = pv[(j - 1) * 4 + 2], z1 = pv[j * 4 + 2];
      if (z0 < dd && z1 >= dd) {
        const f = z1 > z0 ? (dd - z0) / (z1 - z0) : 1, q0 = (j - 1) * 4, q1 = j * 4;
        X.push(pv[q0] + f * (pv[q1] - pv[q0]) - (sx + ax * d)); Y.push(pv[q0 + 1] + f * (pv[q1 + 1] - pv[q0 + 1]) - (sy + ay * d)); E.push(pv[q0 + 3] + f * (pv[q1 + 3] - pv[q0 + 3])); break;
      }
    }
  }
  return { X, Y, E, tracked };
}
const gauss = (x, m, s) => Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI));
const NEED = '<p class="hint">Run a simulation to see this plot.</p>';

const PLOTS = {
  stop: {
    name: 'Stopping power & range', nores: true,
    defaults: { lay: 'all', quant: 'stop', sU: 'mass', rU: 'g/cm²', eax: 'MeV', emin: '', emax: '', cTot: true, cEl: true, cNu: true, cRad: false, cmp: false, ylog: 'log', xlog: 'log' },
    ctrl: s => [
      { k: 'quant', label: 'Quantity', type: 'select', opts: [['stop', 'Stopping power'], ['range', 'CSDA range']] }, { k: 'lay', label: 'Layer', type: 'select', opts: layOpts() },
      { k: 'sU', label: 'Stopping unit', type: 'select', opts: [['mass', 'MeV cm²/g'], ['kum', 'keV/µm'], ['mcm', 'MeV/cm'], ['srim', 'eV/(10¹⁵ atoms/cm²)']] },
      { k: 'rU', label: 'Range unit', type: 'select', opts: [['g/cm²', 'g/cm²'], ['mg/cm²', 'mg/cm²'], ['µm', 'µm'], ['mm', 'mm'], ['cm', 'cm']] },
      { k: 'eax', label: 'Energy axis', type: 'select', opts: [['MeV', 'total energy (MeV)'], ['MeV/u', 'per nucleon (MeV/u)']] },
      { k: 'emin', label: 'E min (MeV)', type: 'text' }, { k: 'emax', label: 'E max (MeV)', type: 'text' },
      { k: 'cTot', label: 'total', type: 'check' }, { k: 'cEl', label: 'electronic', type: 'check' }, { k: 'cNu', label: 'nuclear', type: 'check' }, { k: 'cRad', label: 'radiative (e⁻)', type: 'check' }, { k: 'cmp', label: 'compare the stopping-power data sets (electronic, first layer)', type: 'check' }],
    build(s) {
      const spec = beamSpec(), isE = spec.kind === 'e', A = isE ? 1 : spec.A, E0 = beamE0();
      const emin = pnum(s.emin) > 0 ? pnum(s.emin) : (isE ? 0.001 : 0.001 * A), emax = pnum(s.emax) > emin ? pnum(s.emax) : Math.max(E0 * 3, emin * 100);
      const idx = s.lay === 'all' ? S.layers.map((l, i) => i) : [+s.lay]; const series = []; let note = ''; let ci = 0;
      const eaxis = s.eax === 'MeV/u' && !isE;
      for (const i of idx) {
        const l = S.layers[i]; if (!l || !l.def) { continue; }
        let m, tb; try { m = matOf(l.def); tb = PT.stopTable(m, spec, emin, emax, 160, { stopSet: S.phys.stopSet || 'auto' }); } catch (e) { note += '<p class="err">' + esc(layerName(l, i)) + ': ' + esc(e.message) + '</p>'; continue; }
        const x = tb.E.map(v => eaxis ? v / A : v), rho = m.rho;
        const fS = s.sU === 'mass' ? 1 : s.sU === 'kum' ? rho * 0.1 : s.sU === 'mcm' ? rho : m.Abar * 1.66054e-3;
        const nm = idx.length > 1 ? (i + 1) + ' · ' + layerName(l, i) : '';
        if (s.quant === 'range') {
          const fR = s.rU === 'g/cm²' ? 1 : s.rU === 'mg/cm²' ? 1e3 : s.rU === 'µm' ? 1e4 / rho : s.rU === 'mm' ? 10 / rho : 1 / rho;
          series.push({ name: nm || 'CSDA range', x, y: tb.R.map(v => v * fR), color: scol(ci++) });
        } else if (idx.length > 1) series.push({ name: nm, x, y: tb.Stot.map(v => v * fS), color: scol(ci++) });
        else {
          if (s.cTot) series.push({ name: 'Total', x, y: tb.Stot.map(v => v * fS), color: scol(0), w: 2.4 });
          if (s.cEl) series.push({ name: 'Electronic', x, y: tb.Se.map(v => v * fS), color: scol(1), dash: '6 3' });
          if (s.cNu && !isE) series.push({ name: 'Nuclear', x, y: tb.Sn.map(v => v * fS), color: scol(2), dash: '2 3' });
          if (s.cRad && isE) series.push({ name: 'Radiative', x, y: tb.Srad.map(v => v * fS), color: scol(3), dash: '8 3 2 3' });
        }
      }
      if (s.cmp && s.quant !== 'range') {            // compare the data sets: electronic stopping of the first selected layer
        const l0 = S.layers[idx[0]]; series.length = 0; ci = 0;
        if (l0 && l0.def) try {
          const m = matOf(l0.def), fS0 = s.sU === 'mass' ? 1 : s.sU === 'kum' ? m.rho * 0.1 : s.sU === 'mcm' ? m.rho : m.Abar * 1.66054e-3;
          [['auto', 'Newest available'], ['nist', 'NIST / ICRU 73'], ['semi', 'Semi-empirical']].forEach(([k, nmS], j) => { const tb = PT.stopTable(m, spec, emin, emax, 160, { stopSet: k }); series.push({ name: nmS, x: tb.E.map(v => eaxis ? v / A : v), y: tb.Se.map(v => v * fS0), color: scol(j), dash: j === 1 ? '6 3' : j === 2 ? '2 3' : '', w: j === 0 ? 2.4 : 1.6 }); });
          note += '<p class="hint">Electronic stopping power in ' + esc(layerName(l0, idx[0])) + ' for the three data sets. Curves coincide where only one table exists for the material.</p>';
        } catch (e) { note += '<p class="err">' + esc(e.message) + '</p>'; }
      }
      const sUl = { mass: 'MeV cm²/g', kum: 'keV/µm', mcm: 'MeV/cm', srim: 'eV/(10¹⁵ atoms/cm²)' }[s.sU];
      if (!series.length && !note) note = '<p class="hint">Choose a layer that contains a material.</p>';
      return { series, xlabel: eaxis ? 'Energy (MeV/u)' : 'Kinetic energy (MeV)', ylabel: s.quant === 'range' ? 'CSDA range (' + s.rU + ')' : 'Stopping power (' + sUl + ')', title: beamName() + ' in ' + (idx.length > 1 ? 'all layers' : (S.layers[idx[0]] ? layerName(S.layers[idx[0]], idx[0]) : '')), note, vlines: [{ x: eaxis ? E0 / A : E0, label: 'beam' }] };
    },
  },
  depth: {
    name: 'Depth profile',
    defaults: { lenU: 'auto', yU: 'kum', style: 'line', pTot: true, pPrim: true, pEl: false, pNu: true, pSe: false, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'lenU', label: 'Depth unit', type: 'select', opts: lenSel }, { k: 'yU', label: 'Energy loss unit', type: 'select', opts: [['kum', 'keV/µm  (= eV/nm)'], ['mmm', 'MeV/mm'], ['mcm', 'MeV/cm']] }, { k: 'style', label: 'Draw as', type: 'select', opts: [['line', 'line'], ['step', 'histogram']] },
    { k: 'pTot', label: 'total deposited', type: 'check' }, { k: 'pPrim', label: 'lost by primary', type: 'check' }, { k: 'pEl', label: 'ion electronic', type: 'check' }, { k: 'pNu', label: 'ion nuclear', type: 'check' }, { k: 'pSe', label: 'δ-electron deposit', type: 'check' }],
    build(s) {
      const R = RES, u = unitOf(s, 'lenU', R.zmaxP), f = LEN[u], n = R.nb, fy = { kum: 0.1, mmm: 0.1, mcm: 1 }[s.yU] / (R.N * R.dzb);
      const x = new Float64Array(n), edges = new Float64Array(n + 1); for (let i = 0; i < n; i++) { x[i] = (i + 0.5) * R.dzb / f; edges[i] = i * R.dzb / f; } edges[n] = n * R.dzb / f;
      const mk = (name, arr, ci, dash) => { const y = new Float64Array(n); for (let i = 0; i < n; i++) y[i] = arr[i] * fy; return s.style === 'step' ? { name, type: 'hist', edges, y, color: scol(ci), fill: 0.08 } : { name, x, y, color: scol(ci), dash }; };
      const tot = new Float64Array(n); for (let i = 0; i < n; i++) tot[i] = R.depth.el[i] + R.depth.nu[i] + R.depth.se[i];
      const series = [];
      if (s.pTot) series.push(mk('Total deposited', tot, 0)); if (s.pPrim) series.push(mk('Lost by primary', R.depth.loss, 1, '6 3'));
      if (s.pEl) series.push(mk('Ion electronic', R.depth.el, 2, '2 3')); if (s.pNu) series.push(mk('Ion nuclear', R.depth.nu, 3, '8 3 2 3')); if (s.pSe) series.push(mk('δ-electron deposit', R.depth.se, 4));
      const bands = []; for (let i = 0; i < R.L; i++) bands.push({ x0: R.edges[i] / f, x1: R.edges[i + 1] / f, label: (i + 1) + ' · ' + R.names[i] });
      return { series, bands, xlabel: 'Depth z (' + u + ')', ylabel: 'Energy loss dE/dz per primary (' + { kum: 'keV/µm', mmm: 'MeV/mm', mcm: 'MeV/cm' }[s.yU] + ')', title: R.beamName + ' ' + fmtE(R.E0nom) };
    },
  },
  range: {
    name: 'Stopping positions',
    defaults: { lenU: 'auto', zOn: true, xOn: true, yOn: true, pOn: false, pop: 'stopped', norm: 'prim', bins: 80, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'lenU', label: 'Length unit', type: 'select', opts: lenSel }, { k: 'pop', label: 'Particles', type: 'select', opts: [['stopped', 'stopped in target'], ['all', 'all primaries (final position)']] }, { k: 'norm', label: 'Y axis', type: 'select', opts: normSel }, { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 1000 },
    { k: 'zOn', label: 'z (depth)', type: 'check' }, { k: 'xOn', label: 'x', type: 'check' }, { k: 'yOn', label: 'y', type: 'check' }, { k: 'pOn', label: 'path length', type: 'check' }],
    build(s) {
      const R = RES, u = unitOf(s, 'lenU', R.zmaxP), f = LEN[u], sel = []; for (let i = 0; i < R.N; i++) if (s.pop === 'all' || R.status[i] === 0) sel.push(i);
      const defs = [['zOn', 'z (depth)', R.endZ, 0], ['xOn', 'x', R.endX, 1], ['yOn', 'y', R.endY, 2], ['pOn', 'path length', R.pathL, 3]].filter(a => s[a[0]]);
      const data = defs.map(d => ({ name: d[1], ci: d[3], v: Float64Array.from(sel, i => d[2][i] / f) }));
      let dmin = Infinity, dmax = -Infinity; data.forEach(d => d.v.forEach(v => { if (v < dmin) dmin = v; if (v > dmax) dmax = v; }));
      if (!isFinite(dmin)) { dmin = 0; dmax = 1; } const pad = (dmax - dmin) * 0.03 || 1e-3; const [lo, hi] = binRange(s, dmin - pad, dmax + pad), nb = Math.max(2, Math.round(s.bins)), edges = edgesOf(lo, hi, nb);
      const series = data.map(d => ({ name: d.name, type: 'hist', edges, y: normC(histo(d.v, nb, lo, hi), edges, Math.max(sel.length, 1), s.norm), color: scol(d.ci) }));
      return { series, xlabel: 'Position / path length (' + u + ')', ylabel: ylabNorm(s.norm, u), title: 'Stopping-position distributions · ' + R.beamName + ' ' + fmtE(R.E0nom),
        note: statsTable(data.map(d => ({ name: d.name, d: describe(d.v) })), u) + '<p class="hint">x and y are absolute lateral coordinates (beam starts at the origin; tilted beams drift in x).' + (sel.length === 0 ? ' No particle stopped inside the target.' : '') + '</p>' };
    },
  },
  layerloss: {
    name: 'Energy loss per layer',
    defaults: { lay: '0', what: 'prim', pop: 'all', eU: 'auto', norm: 'prim', bins: 80, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'lay', label: 'Layer', type: 'select', opts: layOpts().slice(1) }, { k: 'what', label: 'Quantity', type: 'select', opts: [['prim', 'energy lost by the primary'], ['dep', 'energy deposited (incl. δ-rays)']] }, { k: 'pop', label: 'Particles', type: 'select', opts: [['all', 'all primaries'], ['hit', 'only those that lose energy in the layer']] },
    { k: 'eU', label: 'Energy unit', type: 'select', opts: [['auto', 'auto'], ['eV', 'eV'], ['keV', 'keV'], ['MeV', 'MeV']] }, { k: 'norm', label: 'Y axis', type: 'select', opts: normSel }, { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 1000 }],
    build(s) {
      const R = RES, l = Math.min(+s.lay, R.L - 1), arr = s.what === 'prim' ? R.dEprim : R.dEdep; let v = []; for (let k = 0; k < R.N; k++) { const e = arr[k * R.L + l]; if (s.pop === 'all' || e > 0) v.push(e); }
      const d0 = describe(v); let u = s.eU; if (u === 'auto') { const m = d0 ? Math.max(d0.mean, 1e-12) : 1; u = m < 1e-3 ? 'eV' : m < 1 ? 'keV' : 'MeV'; }
      const f = EUNITS[u]; const val = Float64Array.from(v, e => e / f); const dd = describe(val);
      let dmin = dd ? dd.min : 0, dmax = dd ? dd.max : 1; const pad = (dmax - dmin) * 0.03 || 1e-3; const [lo, hi] = binRange(s, Math.max(0, dmin - pad), dmax + pad), nb = Math.max(2, Math.round(s.bins)), edges = edgesOf(lo, hi, nb);
      const counts = histo(val, nb, lo, hi); let pk = 0; for (let i = 1; i < nb; i++) if (counts[i] > counts[pk]) pk = i;
      const nm = (s.what === 'prim' ? 'Energy lost by primary' : 'Energy deposited') + ' in layer ' + (l + 1);
      return { series: [{ name: nm, type: 'hist', edges, y: normC(counts, edges, Math.max(val.length, 1), s.norm), color: scol(0) }], vlines: dd ? [{ x: dd.mean, label: 'mean' }] : [], xlabel: nm.replace(/ in layer.*/, '') + ' (' + u + ')', ylabel: ylabNorm(s.norm, u), title: 'Layer ' + (l + 1) + ' · ' + R.names[l] + ' (' + fmtLen(R.edges[l + 1] - R.edges[l]) + ')',
        note: statsTable([{ name: nm, d: dd }], u) + (dd ? '<p class="hint">Most probable value (peak bin): ' + fmt((edges[pk] + edges[pk + 1]) / 2, 5) + ' ' + u + '. Thin layers give asymmetric (Landau/Vavilov-like) shapes with a high-energy tail from δ-rays; thick layers approach a Gaussian.</p>' : '') };
    },
  },
  radial: {
    name: 'Radial & lateral',
    defaults: { mode: 'r', lenU: 'auto', depth: '', bins: 60, norm: 'prim', fit: true, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'mode', label: 'Show', type: 'select', opts: [['r', 'radial distribution r at depth'], ['rdens', 'radial density (per area)'], ['x', 'x distribution at depth'], ['y', 'y distribution at depth'], ['E', 'energy spectrum at depth'], ['rms', 'lateral spread vs depth'], ['reach', 'fraction reaching depth']] },
    { k: 'lenU', label: 'Length unit', type: 'select', opts: lenSel }, { k: 'depth', label: 'Depth z (in length unit)', type: 'number', min: 0, step: 'any' }, { k: 'depthR', label: 'Depth slider', type: 'range' },
    { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 1000 }, { k: 'norm', label: 'Y axis', type: 'select', opts: normSel }, { k: 'fit', label: 'overlay Gaussian / Rayleigh', type: 'check' }],
    build(s) {
      const R = RES, u = unitOf(s, 'lenU', R.zmaxP), f = LEN[u]; let dcm = radialDepthCm(s);
      const info = c => c.tracked < R.N ? '<p class="warn">Only ' + c.tracked.toLocaleString('en') + ' of ' + R.N.toLocaleString('en') + ' primary tracks were stored (storage cap); statistics use these.</p>' : '';
      if (s.mode === 'rms' || s.mode === 'reach') {
        const nd = 100, xs = [], sx = [], sy = [], sr = [], mr = [], fr = [];
        for (let k = 0; k <= nd; k++) {
          const d = R.zmaxP * (k + 0.0001) / nd, c = crossings(R, d); xs.push(d / f); fr.push(c.tracked ? c.X.length / c.tracked : NaN);
          const n = c.X.length; let a = 0, b = 0, r2 = 0, rs = 0, mx = 0, my = 0; for (let i = 0; i < n; i++) { mx += c.X[i]; my += c.Y[i]; } mx /= n || 1; my /= n || 1;
          for (let i = 0; i < n; i++) { a += (c.X[i] - mx) ** 2; b += (c.Y[i] - my) ** 2; r2 += c.X[i] ** 2 + c.Y[i] ** 2; rs += Math.hypot(c.X[i], c.Y[i]); }
          sx.push(n > 1 ? Math.sqrt(a / (n - 1)) / f : NaN); sy.push(n > 1 ? Math.sqrt(b / (n - 1)) / f : NaN); sr.push(n ? Math.sqrt(r2 / n) / f : NaN); mr.push(n ? rs / n / f : NaN);
        }
        const bands = []; for (let i = 0; i < R.L; i++) bands.push({ x0: R.edges[i] / f, x1: R.edges[i + 1] / f, label: (i + 1) + '' });
        if (s.mode === 'reach') return { series: [{ name: 'Fraction of primaries reaching depth', x: xs, y: fr, color: scol(0), w: 2.2 }], bands, xlabel: 'Depth z (' + u + ')', ylabel: 'Fraction of primaries crossing z', title: 'Penetration · ' + R.beamName + ' ' + fmtE(R.E0nom) };
        return { series: [{ name: 'σ_x', x: xs, y: sx, color: scol(0) }, { name: 'σ_y', x: xs, y: sy, color: scol(1), dash: '6 3' }, { name: 'r rms', x: xs, y: sr, color: scol(2), w: 2.4 }, { name: '⟨r⟩', x: xs, y: mr, color: scol(3), dash: '2 3' }], bands, xlabel: 'Depth z (' + u + ')', ylabel: 'Lateral spread (' + u + ')', title: 'Lateral straggling vs depth · ' + R.beamName + ' ' + fmtE(R.E0nom), note: '<p class="hint">Measured relative to the undeflected beam axis; only primaries that reach each depth are counted.</p>' };
      }
      const c = crossings(R, dcm), n = c.X.length, nb = Math.max(2, Math.round(s.bins)); let note = info(c) + '<p class="hint">' + n.toLocaleString('en') + ' of ' + c.tracked.toLocaleString('en') + ' tracked primaries cross z = ' + fmt(dcm / f, 4) + ' ' + u + ' (' + fmt(100 * n / Math.max(c.tracked, 1), 3) + ' %).</p>';
      if (!n) return { series: [], xlabel: '', ylabel: '', note: note + '<p class="err">No primary reaches this depth. Choose a smaller depth.</p>' };
      let series = [], xlabel, ylabel, stats;
      if (s.mode === 'E') {
        const v = Float64Array.from(c.E); const dd = describe(v); const un = dd.mean < 1e-3 ? 'eV' : dd.mean < 1 ? 'keV' : 'MeV', fe = EUNITS[un]; const w = Float64Array.from(v, e => e / fe), d2 = describe(w);
        const pad = (d2.max - d2.min) * 0.03 || 1e-3; const [lo, hi] = binRange(s, Math.max(0, d2.min - pad), d2.max + pad), edges = edgesOf(lo, hi, nb);
        series = [{ name: 'Kinetic energy at depth', type: 'hist', edges, y: normC(histo(w, nb, lo, hi), edges, n, s.norm), color: scol(0) }]; xlabel = 'Kinetic energy at z (' + un + ')'; ylabel = ylabNorm(s.norm, un, 'Fraction of crossing primaries'); stats = statsTable([{ name: 'Energy', d: d2 }], un);
      } else if (s.mode === 'x' || s.mode === 'y') {
        const v = Float64Array.from(s.mode === 'x' ? c.X : c.Y, a => a / f), d2 = describe(v); const pad = (d2.max - d2.min) * 0.03 || 1e-3; const [lo, hi] = binRange(s, d2.min - pad, d2.max + pad), edges = edgesOf(lo, hi, nb);
        const cnt = histo(v, nb, lo, hi); series = [{ name: s.mode + ' at depth', type: 'hist', edges, y: normC(cnt, edges, n, s.norm), color: scol(0) }];
        if (s.fit && d2.sd > 0) { const g = new Float64Array(nb); for (let i = 0; i < nb; i++) { const xc = (edges[i] + edges[i + 1]) / 2; g[i] = n * (edges[i + 1] - edges[i]) * gauss(xc, d2.mean, d2.sd); } series.push({ name: 'Gaussian (same mean, σ)', x: Array.from(edges.slice(0, nb), (e, i) => (e + edges[i + 1]) / 2), y: normC(g, edges, n, s.norm), color: scol(1), dash: '6 3', w: 2 }); }
        xlabel = s.mode + ' at z (' + u + ')'; ylabel = ylabNorm(s.norm, u, 'Fraction of crossing primaries'); stats = statsTable([{ name: s.mode, d: d2 }], u);
      } else {
        const rr = Float64Array.from(c.X, (a, i) => Math.hypot(a, c.Y[i]) / f), d2 = describe(rr); const rms = Math.sqrt(rr.reduce((a, b) => a + b * b, 0) / n); const [lo, hi] = binRange(s, 0, d2.max * 1.03 || 1e-3), edges = edgesOf(Math.max(lo, 0), hi, nb);
        const cnt = histo(rr, nb, edges[0], hi); const area = i => Math.PI * (edges[i + 1] ** 2 - edges[i] ** 2);
        const rd = s.mode === 'rdens'; const y = rd ? Float64Array.from(cnt, (v, i) => s.norm === 'counts' ? v : v / (n * area(i))) : normC(cnt, edges, n, s.norm);
        series = [{ name: rd ? 'Radial density' : 'Radial distribution', type: 'hist', edges, y, color: scol(0) }];
        if (s.fit && rms > 0) { const sg = rms / Math.SQRT2, g = new Float64Array(nb); for (let i = 0; i < nb; i++) { const a = edges[i], b = edges[i + 1]; g[i] = n * (Math.exp(-a * a / (2 * sg * sg)) - Math.exp(-b * b / (2 * sg * sg))); }
          series.push({ name: 'Rayleigh (σ = r_rms/√2)', x: Array.from(edges.slice(0, nb), (e, i) => (e + edges[i + 1]) / 2), y: rd ? Float64Array.from(g, (v, i) => s.norm === 'counts' ? v : v / (n * area(i))) : normC(g, edges, n, s.norm), color: scol(1), dash: '6 3', w: 2 }); }
        xlabel = 'Radius r from the beam axis at z (' + u + ')'; ylabel = rd ? (s.norm === 'counts' ? 'Counts per bin' : 'Fraction of crossing primaries per ' + u + '²') : ylabNorm(s.norm, u, 'Fraction of crossing primaries'); stats = statsTable([{ name: 'r', d: d2 }], u) + '<p class="hint">r rms = ' + fmt(rms, 5) + ' ' + u + '</p>';
      }
      return { series, vlines: [], xlabel, ylabel, title: 'At z = ' + fmt(dcm / f, 4) + ' ' + u + ' · ' + R.beamName + ' ' + fmtE(R.E0nom), note: note + stats };
    },
  },
  flux: {
    name: 'Ion flux',
    defaults: { what: 'sel', lay: 'all', addExit: true, eU: 'auto', norm: 'phi', bins: 80, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'what', label: 'Show', type: 'select', opts: [['sel', 'inside the selected layer'], ['ov', 'inside every layer (overlaid)'], ['exit', 'at the exit of all layers only']] },
      { k: 'lay', label: 'Layer (for “inside the selected layer”)', type: 'select', opts: layOpts() }, { k: 'addExit', label: 'also show the exit spectrum (transmitted ions)', type: 'check' },
      { k: 'eU', label: 'Energy unit', type: 'select', opts: [['auto', 'auto'], ['eV', 'eV'], ['keV', 'keV'], ['MeV', 'MeV'], ['GeV', 'GeV']] },
      { k: 'norm', label: 'Y axis', type: 'select', opts: [['phi', 'differential flux Φ(E) per incident particle'], ['frac', 'fraction per bin'], ['dens', 'probability density']] }, { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 1000 }],
    build(s) {
      const R = RES, L = R.L, pv = R.pv, off = R.pvOff, ed = R.edges; let tracked = 0, emax = 0;
      for (let k = 0; k < R.N; k++) if (off[k + 1] - off[k] >= 1) tracked++;
      const nv = off[R.N]; for (let j = 0; j < nv; j++) { const e = pv[j * 4 + 3]; if (e > emax) emax = e; }
      for (let i = 0; i < R.N; i++) if (R.status[i] === 1 && R.endE[i] > emax) emax = R.endE[i];
      if (!(emax > 0)) return { series: [], xlabel: '', ylabel: '', note: '<p class="err">No stored tracks.</p>' };
      let u = s.eU; if (u === 'auto') u = emax < 1e-3 ? 'eV' : emax < 1 ? 'keV' : emax < 1e3 ? 'MeV' : 'GeV'; const f = EUNITS[u];
      const [lo, hi] = binRange(s, 0, emax / f * 1.02), nb = Math.max(2, Math.round(s.bins)), edges = edgesOf(lo, hi, nb), bw = (hi - lo) / nb;
      const cnt = Array.from({ length: L }, () => new Float64Array(nb)), W = new Float64Array(L), SX = new Float64Array(L), SX2 = new Float64Array(L);
      for (let k = 0; k < R.N; k++) { const a = off[k], b = off[k + 1]; if (b - a < 2) continue;
        for (let j = a + 1; j < b; j++) { const q0 = (j - 1) * 4, q1 = j * 4, zm = (pv[q0 + 2] + pv[q1 + 2]) / 2; let li = -1; for (let i = 0; i < L; i++) if (zm >= ed[i] && zm < ed[i + 1]) { li = i; break; } if (li < 0) continue;
          const len = Math.hypot(pv[q1] - pv[q0], pv[q1 + 1] - pv[q0 + 1], pv[q1 + 2] - pv[q0 + 2]); if (!(len > 0)) continue; const ea = Math.min(pv[q0 + 3], pv[q1 + 3]) / f, eb = Math.max(pv[q0 + 3], pv[q1 + 3]) / f;
          // energy falls linearly along the step: spread its path length uniformly over [ea, eb] (exact bin overlap, no sub-sampling)
          W[li] += len; SX[li] += len * (ea + eb) / 2; SX2[li] += len * (ea * ea + ea * eb + eb * eb) / 3;
          const span = eb - ea;
          if (span < 1e-12 * (eb + 1e-30)) { if (ea >= lo && ea <= hi) { let kk = Math.floor((ea - lo) / bw); if (kk >= nb) kk = nb - 1; cnt[li][kk] += len; } }
          else { let k0 = Math.max(0, Math.floor((ea - lo) / bw)), k1 = Math.min(nb - 1, Math.floor((eb - lo) / bw)); for (let kk = k0; kk <= k1; kk++) { const a = Math.max(ea, lo + kk * bw), b = Math.min(eb, lo + (kk + 1) * bw); if (b > a) cnt[li][kk] += len * (b - a) / span; } } } }
      const series = [], rows = []; const nT = Math.max(tracked, 1);
      const addVol = (name, ids, col) => { const c = new Float64Array(nb); let w = 0, sx = 0, sx2 = 0, th = 0; ids.forEach(i => { for (let b = 0; b < nb; b++) c[b] += cnt[i][b]; w += W[i]; sx += SX[i]; sx2 += SX2[i]; th += ed[i + 1] - ed[i]; }); if (!(w > 0)) return;
        const y = Float64Array.from(c, v => s.norm === 'frac' ? v / w : s.norm === 'dens' ? v / (w * bw) : v / (nT * th * bw)); series.push({ name, type: 'hist', edges, y, color: col });
        const m = sx / w, sd = Math.sqrt(Math.max(sx2 / w - m * m, 0)); rows.push([name, fmt(w / (nT * th), 4), fmt(m, 5), fmt(sd, 4)]); };
      const addExit = col => { const c = new Float64Array(nb); let sx = 0, sx2 = 0, n = 0; for (let i = 0; i < R.N; i++) if (R.status[i] === 1) { const e = R.endE[i] / f; n++; sx += e; sx2 += e * e; if (e >= lo && e <= hi) { let kk = Math.floor((e - lo) / bw); if (kk >= nb) kk = nb - 1; c[kk]++; } } if (!n) return;
        const y = Float64Array.from(c, v => s.norm === 'frac' ? v / n : s.norm === 'dens' ? v / (n * bw) : v / (R.N * bw)); series.push({ name: 'Exit (transmitted)', type: 'hist', edges, y, color: scol(L + 1) });
        const m = sx / n, sd = Math.sqrt(Math.max(sx2 / n - m * m, 0)); rows.push(['Exit (transmitted)', fmt(n / R.N, 4), fmt(m, 5), fmt(sd, 4)]); };
      if (s.what === 'ov') { for (let i = 0; i < L; i++) addVol((i + 1) + ' · ' + R.names[i], [i], scol(i)); if (L > 1) addVol('All layers', Array.from({ length: L }, (_, i) => i), scol(L)); }
      else if (s.what === 'sel') { if (s.lay === 'all') addVol('All layers', Array.from({ length: L }, (_, i) => i), scol(L)); else { const i = Math.min(+s.lay, L - 1); addVol((i + 1) + ' · ' + R.names[i], [i], scol(i)); } }
      if (s.what === 'exit' || s.addExit) addExit(scol(L + 1));
      if (!series.length) return { series: [], xlabel: '', ylabel: '', note: '<p class="err">No ions in this selection (no stored track crosses the layer / nothing is transmitted).</p>' };
      const ylabel = s.norm === 'frac' ? 'Fraction per bin' : s.norm === 'dens' ? 'Probability density (1/' + u + ')' : 'Ion flux Φ(E)/Φ₀ per ' + u + ' (relative to the incident flux)';
      const table = '<details class="cw"><summary><b>Integrals and moments of the spectra</b></summary><div class="tscroll"><table class="t"><tr><th>Spectrum</th><th>∫Φ dE</th><th>Mean E (' + u + ')</th><th>σ (' + u + ')</th></tr>' + rows.map(r => '<tr><td>' + esc(r[0]) + '</td><td>' + r[1] + '</td><td>' + r[2] + '</td><td>' + r[3] + '</td></tr>').join('') + '</table></div></details>';
      return { series, xlabel: 'Ion kinetic energy (' + u + ')', ylabel, title: 'Energy distribution of the ion flux · ' + R.beamName + ' ' + fmtE(R.E0nom),
        note: (tracked < R.N ? '<p class="warn">Only ' + tracked.toLocaleString('en') + ' of ' + R.N.toLocaleString('en') + ' primary tracks were stored; the in-layer spectra use these (the exit spectrum uses all primaries).</p>' : '') +
          '<p class="hint">Inside a layer Φ(E) is the track-length (volume) flux: total path length of the ions in the layer divided by its thickness, per incident particle, i.e. the ion fluence per unit incident fluence (laterally infinite slab). ∫Φ dE is the mean path length per unit thickness (1 for an undisturbed beam at normal incidence). At the exit it is the number of transmitted ions per incident ion (plane-crossing, without 1/cosθ weighting). Only primary ions; δ-electrons are not included.</p>' + table };
    },
  },
  delta: {
    name: 'Knock-on electrons',
    defaults: { lay: 'all', norm: 'dT', bins: 60, ylog: 'log', xlog: 'log' },
    ctrl: s => [{ k: 'lay', label: 'Layer', type: 'select', opts: layOpts() }, { k: 'norm', label: 'Y axis', type: 'select', opts: [['dT', 'dN/dT per primary (1/keV)'], ['dlnT', 'dN/dlnT per primary'], ['counts', 'counts per bin']] }, { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 500 }],
    build(s) {
      const R = RES, T = []; for (let i = 0; i < R.delE.length; i++) if (s.lay === 'all' || R.delL[i] === +s.lay) T.push(R.delE[i] * 1000);
      if (!T.length) return { series: [], xlabel: '', ylabel: '', note: '<p class="hint">No δ-electrons were produced. Lower the δ-ray threshold, use a heavier/faster projectile or a thicker layer.</p>' };
      let mn = Infinity, mx = 0; T.forEach(v => { if (v < mn) mn = v; if (v > mx) mx = v; }); mn = pnum(s.xmin) > 0 ? pnum(s.xmin) : mn * 0.999; mx = pnum(s.xmax) > mn ? pnum(s.xmax) : mx * 1.001;
      const nb = Math.max(2, Math.round(s.bins)), edges = edgesOf(mn, mx, nb, true), cnt = new Float64Array(nb), lr = Math.log(mx / mn);
      T.forEach(v => { if (v >= mn && v <= mx) { let k = Math.floor(Math.log(v / mn) / lr * nb); if (k >= nb) k = nb - 1; cnt[k]++; } });
      const y = Float64Array.from(cnt, (v, i) => s.norm === 'counts' ? v : s.norm === 'dT' ? v / R.N / (edges[i + 1] - edges[i]) : v / R.N / Math.log(edges[i + 1] / edges[i]));
      const xc = Array.from(cnt, (_, i) => Math.sqrt(edges[i] * edges[i + 1]));
      return { series: [{ name: 'δ-electrons', x: xc, y, type: 'pts', color: scol(0) }, { name: '', x: xc, y, color: scol(0), noread: false }], xlabel: 'δ-electron kinetic energy T (keV)', ylabel: s.norm === 'dT' ? 'dN/dT per primary (1/keV)' : s.norm === 'dlnT' ? 'dN/dlnT per primary' : 'Counts per bin', title: 'Knock-on electron spectrum · ' + R.beamName + ' ' + fmtE(R.E0nom),
        note: '<p class="hint">' + T.length.toLocaleString('en') + ' δ-electrons above the ' + fmtE(R.phys.Tcut) + ' threshold (' + fmt(T.length / R.N, 4) + ' per primary). Below the threshold the energy is part of the continuous (restricted) loss. The spectrum falls like 1/T² up to T<sub>max</sub>.</p>' };
    },
  },
  react: {
    name: 'Nuclear reactions',
    defaults: { lenU: 'auto', bins: 60, norm: 'prim', ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'lenU', label: 'Depth unit', type: 'select', opts: lenSel }, { k: 'norm', label: 'Y axis', type: 'select', opts: normSel }, { k: 'bins', label: 'Bins', type: 'number', min: 5, max: 500 }],
    build(s) {
      const R = RES, u = unitOf(s, 'lenU', R.zmaxP), f = LEN[u], z = Float64Array.from(R.rx.z, v => v / f);
      if (!z.length) return { series: [], xlabel: '', ylabel: '', note: '<p class="hint">' + (R.phys.reactions ? 'No nuclear reaction occurred in this run (cross sections are small for light, fast particles in thin targets — try more primaries or a higher σ scale).' : 'Nuclear reactions are switched off.') + '</p>' };
      const d = describe(z), pad = (d.max - d.min) * 0.03 || 1e-3, [lo, hi] = binRange(s, Math.max(0, d.min - pad), d.max + pad), nb = Math.max(2, Math.round(s.bins)), edges = edgesOf(lo, hi, nb);
      const rows = R.names.map((nm, i) => '<tr><td>' + (i + 1) + ' · ' + esc(nm) + '</td><td>' + fmt(R.lay.nReact[i], 4) + '</td><td>' + fmt(100 * R.lay.nReact[i] / R.N, 3) + ' %</td></tr>').join('');
      return { series: [{ name: 'Reaction depth', type: 'hist', edges, y: normC(histo(z, nb, lo, hi), edges, R.N, s.norm), color: scol(2) }], xlabel: 'Depth of the nuclear interaction (' + u + ')', ylabel: ylabNorm(s.norm, u), title: 'Nuclear interaction depth · ' + R.beamName + ' ' + fmtE(R.E0nom),
        note: '<details class="cw"><summary><b>Nuclear reactions per layer</b></summary><div class="tscroll"><table class="t"><tr><th>Layer</th><th>Reactions</th><th>per primary</th></tr>' + rows + '</table></div></details><p class="hint">Simple model: the projectile is removed at the reaction point (geometric + Coulomb-barrier cross section). Secondary-particle production will be added together with the cross-section databases.</p>' };
    },
  },
  tracks: {
    name: 'Tracks',
    defaults: { proj: 'x', lenU: 'auto', latU: 'auto', deltas: true, equal: false, ylog: 'lin', xlog: 'lin' },
    ctrl: s => [{ k: 'proj', label: 'Projection', type: 'select', opts: [['x', 'x – z'], ['y', 'y – z'], ['r', 'r – z']] }, { k: 'lenU', label: 'z unit', type: 'select', opts: lenSel }, { k: 'latU', label: 'Lateral unit', type: 'select', opts: lenSel }, { k: 'deltas', label: 'show δ-electron tracks', type: 'check' }, { k: 'equal', label: 'equal axis scales', type: 'check' }],
    build(s) {
      const R = RES, D = R.disp, u = unitOf(s, 'lenU', R.zmaxP), f = LEN[u], total = D.v.length / 3; const nT = D.off.length;
      let lat = 0; for (let i = 0; i < D.v.length; i += 3) lat = Math.max(lat, Math.abs(D.v[i]), Math.abs(D.v[i + 1])); const ul = unitOf(s, 'latU', lat || R.zmaxP / 10), fl = LEN[ul];
      const px = [], py = [], dx = [], dy = [];
      for (let t = 0; t < nT; t++) {
        const a = D.off[t], b = t + 1 < nT ? D.off[t + 1] : total, isE = D.type[t] === 1; if (isE && !s.deltas) continue;
        const X = isE ? dx : px, Y = isE ? dy : py;
        for (let j = a; j < b; j++) { X.push(D.v[j * 3 + 2] / f); const xx = D.v[j * 3], yy = D.v[j * 3 + 1]; Y.push((s.proj === 'x' ? xx : s.proj === 'y' ? yy : Math.hypot(xx, yy)) / fl); }
        X.push(NaN); Y.push(NaN);
      }
      const bands = []; for (let i = 0; i < R.L; i++) bands.push({ x0: R.edges[i] / f, x1: R.edges[i + 1] / f, label: (i + 1) + ' · ' + R.names[i] });
      const o = { series: [{ name: 'δ-electrons', x: dx, y: dy, color: 'var(--ink2)', w: 0.8, op: 0.55, noread: true, hidden: !s.deltas || !dx.length }, { name: 'Primary tracks', x: px, y: py, color: scol(1), w: 1.3, op: 0.85, noread: true }], bands, xlabel: 'Depth z (' + u + ')', ylabel: (s.proj === 'r' ? 'Radius r' : s.proj) + ' (' + ul + ')', title: 'First ' + (D.pid.length ? new Set(D.pid).size : 0) + ' primaries · ' + R.beamName + ' ' + fmtE(R.E0nom) };
      if (s.equal && u === ul) { const xr = pnum(s.xmax) != null ? pnum(s.xmax) - (pnum(s.xmin) || 0) : R.zmaxP / f; const h = xr * 0.484 / 2; if (s.proj === 'r') { o.ymin = 0; o.ymax = 2 * h; } else { o.ymin = -h; o.ymax = h; } }
      o.note = '<p class="hint">Showing ' + new Set(D.pid).size + ' stored primaries (change “Tracks drawn” under Simulation → More options). Lateral and depth scales differ unless “equal axis scales” is ticked.</p>';
      return o;
    },
  },
};
const TABS = ['stop', 'depth', 'range', 'layerloss', 'radial', 'flux', 'delta', 'react', 'tracks'];
const hidden = {};   // legend visibility per plot tab: hidden[tab][name]=true

function multiLabel(c, cur) { const n = c.opts.filter(o => cur.has(String(o[0]))).length; if (!n) return 'none selected'; if (n === 1) return esc(c.opts.find(o => cur.has(String(o[0])))[1]); return n + ' of ' + c.opts.length + ' selected'; }
function multiSet(k, vals) {
  const s = pset(S.tab), c = PLOTS[S.tab].ctrl(s).find(q => q.k === k); s[k] = vals; ['xmin', 'xmax', 'ymin', 'ymax'].forEach(q => s[q] = '');
  const sm = $('#msel-' + k + ' summary'); if (sm && c) sm.innerHTML = multiLabel(c, new Set(vals.map(String)));
  drawChart(); saveLocal();
}
document.addEventListener('change', ev => {
  const t = ev.target; if (!t.dataset || !t.dataset.pm) return; const k = t.dataset.pm;
  multiSet(k, $$('[data-pm="' + k + '"]').filter(x => x.checked).map(x => x.value));
});
document.addEventListener('click', ev => {
  const b = ev.target.closest && ev.target.closest('[data-mall],[data-mnone]'); if (!b) return; const k = b.dataset.mall || b.dataset.mnone, all = b.dataset.mall != null;
  $$('[data-pm="' + k + '"]').forEach(x => x.checked = all); multiSet(k, all ? $$('[data-pm="' + k + '"]').map(x => x.value) : []);
});
function renderPlots(full) {
  try { if (typeof actRender === 'function') actRender(); } catch (e) { console.error(e); }
  const tab = TABS.includes(S.tab) ? S.tab : 'depth'; S.tab = tab;
  const sp = PLOTS[tab], s = pset(tab);
  if (tab === 'layerloss' && +s.lay >= S.layers.length) s.lay = '0';
  if (tab === 'radial' && RES) { const u0 = unitOf(s, 'lenU', RES.zmaxP); s.depth = +(radialDepthCm(s) / LEN[u0]).toPrecision(4); }
  const ctl = sp.ctrl(s);
  const field = c => {
    if (c.type === 'check') return '<label class="chk"><input type="checkbox" data-pk="' + c.k + '"' + (s[c.k] ? ' checked' : '') + '> ' + c.label + '</label>';
    let inp;
    if (c.type === 'multi') {
      const cur = new Set((Array.isArray(s[c.k]) ? s[c.k] : []).map(String));
      inp = '<details class="msel" id="msel-' + c.k + '"><summary>' + multiLabel(c, cur) + '</summary><div class="mselbox"><div class="btns"><button type="button" class="btn sm" data-mall="' + c.k + '">All</button><button type="button" class="btn sm" data-mnone="' + c.k + '">Clear</button></div>' +
        c.opts.map(o => '<label class="chk"><input type="checkbox" data-pm="' + c.k + '" value="' + esc(o[0]) + '"' + (cur.has(String(o[0])) ? ' checked' : '') + '> ' + esc(o[1]) + '</label>').join('') + '</div></details>';
      return '<label class="fld wide"><span>' + c.label + '</span></label>' + inp;
    }
    if (c.type === 'select') inp = '<select data-pk="' + c.k + '">' + c.opts.map(o => '<option value="' + esc(o[0]) + '"' + (String(o[0]) === String(s[c.k]) ? ' selected' : '') + '>' + esc(o[1]) + '</option>').join('') + '</select>';
    else if (c.type === 'range') { const u = unitOf(s, 'lenU', RES ? RES.zmaxP : 0.1), mx = RES ? RES.zmaxP / LEN[u] : 1; inp = '<input type="range" min="0" max="' + mx + '" step="' + mx / 400 + '" value="' + (pnum(s.depth) || 0) + '" data-pk="depthR" style="width:100%">'; }
    else inp = '<input type="' + (c.type === 'number' ? 'number' : 'text') + '" data-pk="' + c.k + '" value="' + esc(s[c.k] == null ? '' : s[c.k]) + '"' + (c.min != null ? ' min="' + c.min + '"' : '') + (c.max != null ? ' max="' + c.max + '"' : '') + (c.step ? ' step="' + c.step + '"' : '') + (c.type === 'text' ? ' placeholder="auto"' : '') + '>';
    return '<label class="fld"><span>' + c.label + '</span>' + inp + '</label>';
  };
  const ax = [{ k: 'xmin', label: 'X min', type: 'text' }, { k: 'xmax', label: 'X max', type: 'text' }, { k: 'ymin', label: 'Y min', type: 'text' }, { k: 'ymax', label: 'Y max', type: 'text' },
  { k: 'xlog', label: 'X scale', type: 'select', opts: [['lin', 'Linear'], ['log', 'Logarithmic']] }, { k: 'ylog', label: 'Y scale', type: 'select', opts: [['lin', 'Linear'], ['log', 'Logarithmic']] }];
  const lab = [{ k: 'title', label: 'Title (blank = automatic)', type: 'text' }, { k: 'xlabel', label: 'X label', type: 'text' }, { k: 'ylabel', label: 'Y label', type: 'text' }];
  if (full !== false) {
    $('#cPlots').innerHTML = '<h2>Interactions with the target layers</h2><div class="ptabs" role="tablist">' + TABS.map(t => '<button role="tab" data-tab="' + t + '" class="' + (t === tab ? 'sel' : '') + '">' + PLOTS[t].name + '</button>').join('') + '</div>' +
      '<div class="pctl">' + ctl.map(field).join('') + '</div>' +
      '<details open><summary>Axes, ranges and labels</summary><div class="pctl" style="margin-top:10px">' + ax.map(field).join('') + lab.map((c, i) => field(c).replace('<label class="fld">', '<label class="fld wide">') + (i === 0 ? '<div style="grid-column:1/-1"><button class="btn sm" data-preset="reset">Reset axes &amp; labels</button></div>' : '')).join('') + '</div></details>' +
      '<div class="chart" id="chartBox"></div><div class="legend" id="legend"></div><div class="readout" id="readout">Move the pointer over the plot to read values.</div><div id="plotNote"></div>' +
      '<div class="btns"><span class="lbl">Save plot</span><button class="btn sm" data-exp="png">PNG</button><button class="btn sm" data-exp="pdf">PDF</button><button class="btn sm" data-exp="svg">SVG</button><button class="btn sm" data-exp="txt">ASCII (.txt)</button><button class="btn sm" data-exp="csv">CSV</button><button class="btn sm" data-exp="rootbin" title="Binary ROOT file with one TH1D per series">ROOT file</button><button class="btn sm" data-exp="root" title="C++ macro that creates the .root file when run with ROOT">ROOT macro</button><label class="lbl" style="display:inline-flex;gap:6px;align-items:center">PNG width <select id="expRes" style="width:auto;min-width:96px;height:30px;padding-right:22px"><option value="1200">1200 px</option><option value="1800" selected>1800 px</option><option value="3000">3000 px</option></select></label><button class="btn sm" id="csvBtn">Copy CSV</button></div><p class="hint" id="csvMsg"></p>';
  }
  drawChart();
}
function drawChart() {
  const tab = S.tab, sp = PLOTS[tab], s = pset(tab), box = $('#chartBox'); if (!box) return;
  if (!sp.nores && !RES) { box.innerHTML = NEED; $('#legend').innerHTML = ''; $('#plotNote').innerHTML = ''; _chart = null; return; }
  let o; try { o = sp.build(s); } catch (e) { box.innerHTML = '<p class="err">' + esc(e.message) + '</p>'; console.error(e); return; }
  const hid = hidden[tab] = hidden[tab] || {};
  o.series.forEach(se => { if (hid[se.name]) se.hidden = true; });
  o.xlog = s.xlog === 'log'; o.ylog = s.ylog === 'log'; o.xmin = s.xmin; o.xmax = s.xmax; o.ymin = o.ymin != null && !s.ymin ? o.ymin : s.ymin; o.ymax = o.ymax != null && !s.ymax ? o.ymax : s.ymax;
  if (s.xlabel) o.xlabel = s.xlabel; if (s.ylabel) o.ylabel = s.ylabel; if (s.title !== undefined && s.title !== '') o.title = s.title;
  if (!o.series.length) { box.innerHTML = ''; $('#legend').innerHTML = ''; $('#plotNote').innerHTML = o.note || ''; if (typeof secHydrate === 'function') secHydrate(); _chart = null; return; }
  const r = renderChart(o); box.innerHTML = r.svg; _chart = { o, geom: r.geom };
  $('#legend').innerHTML = o.series.filter(se => se.name).map(se => '<button data-leg="' + esc(se.name) + '" class="' + (hid[se.name] ? 'off' : '') + '"><i style="background:' + (se.color || 'var(--accent)') + '"></i>' + esc(se.name) + '</button>').join('');
  const man = ['xmin', 'xmax', 'ymin', 'ymax'].filter(k => pnum(s[k]) != null);
  let ban = ''; if (man.length) ban = '<p class="warn">Manual axis range active (' + man.map(k => k.replace('min', ' min').replace('max', ' max').toUpperCase() + ' = ' + s[k]).join(', ') + ') — data outside it is not shown. <button class="btn sm" data-zoom="reset">Reset range</button></p>';
  if (man.length && o.series.filter(se => se.type === 'hist' && !se.hidden).length && o.series.filter(se => se.type === 'hist' && !se.hidden).every(se => Array.from(se.y).every(v => !(v > 0)))) ban += '<p class="err">All bins inside this range are empty — the data lie outside the manual axis range.</p>';
  $('#plotNote').innerHTML = ban + (o.note || '');
  if (typeof secHydrate === 'function') secHydrate();
  const svg = $('svg', box), xh = $('#xh', svg); attachZoom(svg);
  svg.addEventListener('pointermove', ev => { const rd = chartReadout(_chart.o, _chart.geom, svg, ev); const ro = $('#readout'); if (!rd) { xh.style.display = 'none'; return; } xh.style.display = ''; const rect = svg.getBoundingClientRect(), px = (ev.clientX - rect.left) / rect.width * _chart.geom.W; xh.setAttribute('x1', px); xh.setAttribute('x2', px); ro.textContent = rd.text; });
  svg.addEventListener('pointerleave', () => { xh.style.display = 'none'; });
}
document.addEventListener('input', ev => {
  const t = ev.target, k = t.dataset && t.dataset.pk; if (!k) return;
  const s = pset(S.tab); if (t.type === 'checkbox') return;
  if (k === 'depthR') { s.depth = +parseFloat(t.value).toPrecision(4); const di = $('[data-pk="depth"]'); if (di) di.value = s.depth; }
  else s[k] = t.type === 'number' ? (t.value === '' ? '' : parseFloat(t.value)) : t.value;
  if ((k === 'depthR' || k === 'depth') && RES) { const v = parseFloat(s.depth); s.depthCm = v >= 0 ? v * LEN[unitOf(s, 'lenU', RES.zmaxP)] : 0; }
  drawChart(); saveLocal();
});
document.addEventListener('change', ev => {
  const t = ev.target, k = t.dataset && t.dataset.pk; if (!k) return;
  const s = pset(S.tab); if (t.tagName === 'SELECT' && !['xlog', 'ylog', 'norm', 'style'].includes(k)) ['xmin', 'xmax', 'ymin', 'ymax'].forEach(q => s[q] = '');
  if (t.type === 'checkbox') s[k] = t.checked; else if (t.tagName === 'SELECT') { s[k] = t.value;}
  else return;
  renderPlots(true); saveLocal();
});
document.addEventListener('click', ev => {
  const t = ev.target;
  const tb = t.closest('[data-tab]'); if (tb) { S.tab = tb.dataset.tab; renderPlots(true); saveLocal(); return; }
  const lg = t.closest('[data-leg]'); if (lg) { const h = hidden[S.tab] = hidden[S.tab] || {}; h[lg.dataset.leg] = !h[lg.dataset.leg]; drawChart(); return; }
  if (t.closest('[data-preset="reset"]')) { const s = pset(S.tab); ['xmin', 'xmax', 'ymin', 'ymax', 'xlabel', 'ylabel', 'title'].forEach(k => s[k] = ''); renderPlots(true); saveLocal(); return; }
  if (t.id === 'csvBtn' && _chart) {
    const txt = csvOf(); const m = $('#csvMsg'); (navigator.clipboard && navigator.clipboard.writeText ? navigator.clipboard.writeText(txt) : Promise.reject()).then(() => m.textContent = 'Copied (' + txt.length.toLocaleString('en') + ' characters).', () => { m.textContent = 'Clipboard blocked — select the text below.'; const ta = document.createElement('textarea'); ta.value = txt; ta.rows = 8; ta.style.marginTop = '8px'; $('#plotNote').appendChild(ta); ta.select(); });
  }
});
