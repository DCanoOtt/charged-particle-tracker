/* ====== Neutron detector: window "Configure neutron detector" (geometry sketch, parameters, efficiency, detected distribution, pulse-height spectrum) ====== */
const ND_DEF = {
  mat: 'BC501A',
  D: 10,
  T: 10,
  th: 100,
  N: 20000,
  Emax: 200,
  Ew: 2,
  Nph: 300000,
  rho: 0.874,
  hc: 1.212,
  lpo: {},
  dist: 'E',
  ia: 0,
  eyl: 'lin',
  exl: 'log',
  dyl: 'lin',
  dxl: 'lin',
  hyl: 'log',
  hxl: 'lin',
  ryl: 'lin',
  rxl: 'lin'
};
const ND = {
  data: { H: ndEmb('H'), C: ndEmb('C'), nH: 'ENDF/B-VIII.0 1001.800nc', nC: 'ENDF/B-VIII.0 6012.800nc' },
  eff: null,
  pl: { e: { hidden: {}, z: {} }, d: { hidden: {}, z: {} }, h: { hidden: {}, z: {} }, r: { hidden: {}, z: {} } },
  raw: {},
  gen: 0,
  busy: false,
  ph: null
};
function ndSet() {
  const s = neuSet(),
    o = s.det && typeof s.det === 'object' ? s.det : (s.det = {});
  for (const k of Object.keys(ND_DEF)) if (o[k] === undefined) o[k] = ND_DEF[k];
  if (o.Emax === 40) o.Emax = 200; // former default
  if (!ND_MAT[o.mat]) {
    o.mat = 'BC501A';
    o.rho = ND_MAT.BC501A.rho;
    o.hc = ND_MAT.BC501A.hc;
  }
  return o;
}
const ndLpAll = d => Object.assign({}, ND_LP[(ND_MAT[d.mat] || ND_MAT.BC501A).curve], d.lpo || {});
function ndDet() {
  const d = ndSet(),
    s = neuSet();
  return { mat: d.mat, D: +d.D, T: +d.T, th: +d.th, rho: +d.rho, hc: +d.hc, lp: ndLpAll(d), Lm: +s.L > 0 ? +s.L : 25 };
}

// side view of the point source and the cylindrical detector, labelled with the current values
function ndSketch() {
  const d = ndSet(),
    s = neuSet(),
    T = (x, y, t, a, c) =>
      '<text x="' + x + '" y="' + y + '" text-anchor="' + (a || 'middle') + '" font-size="11.5" fill="' + (c || 'var(--ink)') + '">' + t + '</text>',
    dim = (x0, x1, y, lab, val) =>
      '<line x1="' +
      x0 +
      '" y1="' +
      y +
      '" x2="' +
      x1 +
      '" y2="' +
      y +
      '" stroke="var(--accent)" stroke-width="1.4" marker-start="url(#nda)" marker-end="url(#nda)"/>' +
      T((x0 + x1) / 2, y + 15, lab) +
      T((x0 + x1) / 2, y + 28, val, 'middle', 'var(--ink2)'),
    ext = (x, y0, y1) => '<line x1="' + x + '" y1="' + y0 + '" x2="' + x + '" y2="' + y1 + '" stroke="var(--ink2)" stroke-width=".8" stroke-dasharray="2 3"/>';
  return (
    '<svg viewBox="0 0 660 300" role="img" aria-label="Geometry of the point source and the cylindrical neutron detector" style="width:100%;max-width:660px;display:block;margin:4px auto 0">' +
    '<defs><marker id="nda" viewBox="0 0 8 8" refX="4" refY="4" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L8 4 L0 8 z" fill="var(--accent)"/></marker></defs>' +
    '<rect x="400" y="60" width="130" height="150" fill="var(--accent-soft, rgba(60,140,160,.18))" stroke="var(--accent)" stroke-width="1.6"/>' +
    '<ellipse cx="530" cy="135" rx="22" ry="75" fill="var(--accent-soft, rgba(60,140,160,.18))" stroke="var(--accent)" stroke-width="1.6"/>' +
    '<ellipse cx="400" cy="135" rx="22" ry="75" fill="none" stroke="var(--accent)" stroke-width="1" stroke-dasharray="3 3"/>' +
    '<line x1="20" y1="135" x2="560" y2="135" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="8 4 2 4"/>' +
    '<line x1="40" y1="135" x2="400" y2="60" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="3 3"/><line x1="40" y1="135" x2="400" y2="210" stroke="var(--ink2)" stroke-width="1" stroke-dasharray="3 3"/>' +
    '<circle cx="40" cy="135" r="5" fill="var(--bad, #c0392b)"/>' +
    T(40, 120, 'point source') +
    T(12, 158, 'neutrons from the target', 'start', 'var(--ink2)') +
    T(465, 40, 'scintillator (cylinder)') +
    T(465, 100, '<tspan font-weight="700">' + esc((ND_MAT[d.mat] || {}).name || d.mat) + '</tspan>', 'middle', 'var(--ink)') +
    '<line x1="560" y1="60" x2="560" y2="210" stroke="var(--accent)" stroke-width="1.4" marker-start="url(#nda)" marker-end="url(#nda)"/>' +
    T(568, 130, 'diameter D', 'start') +
    T(568, 146, +d.D + ' cm', 'start', 'var(--ink2)') +
    ext(400, 215, 245) +
    ext(530, 215, 245) +
    ext(40, 140, 245) +
    dim(40, 400, 250, 'flight path L (source – front face)', +s.L + ' m') +
    dim(400, 530, 215, 'thickness t', +d.T + ' cm') +
    T(465, 292, 'the neutrons enter through the front face and may leave through the side wall', 'middle', 'var(--ink2)') +
    '</svg>'
  );
}
function ndInputs() {
  const d = ndSet(),
    s = neuSet(),
    lp = ndLpAll(d),
    nI = (a, k, v, ex) => '<input type="number" ' + a + '="' + k + '" step="any" ' + (ex || '') + ' value="' + v + '">',
    fld = (l, h) => '<label class="fld"><span>' + l + '</span>' + h + '</label>';
  return (
    '<div class="row3">' +
    fld(
      'Scintillator',
      '<select data-dk="mat">' +
        Object.keys(ND_MAT)
          .map(k => '<option value="' + k + '"' + (d.mat === k ? ' selected' : '') + '>' + esc(ND_MAT[k].name) + '</option>')
          .join('') +
        '</select>'
    ) +
    fld('Density (g/cm³)', nI('data-dk', 'rho', d.rho, 'min="0.1"')) +
    fld('H / C atom ratio', nI('data-dk', 'hc', d.hc, 'min="0.1"')) +
    '</div><div class="row3">' +
    fld('Diameter D (cm)', nI('data-dk', 'D', d.D, 'min="0.1"')) +
    fld('Thickness t (cm)', nI('data-dk', 'T', d.T, 'min="0.01"')) +
    fld('Threshold (keVee)', nI('data-dk', 'th', d.th, 'min="0"')) +
    '</div><div class="row3">' +
    fld('Flight path L (m) — source to front face', nI('data-dkn', 'L', s.L, 'min="0.01"')) +
    fld('Time resolution, FWHM (ns)', nI('data-dkn', 'dt', s.dt, 'min="0"')) +
    fld('Histories per energy', nI('data-dk', 'N', d.N, 'min="1000" step="1000"')) +
    '</div><div class="row3">' +
    fld('Largest neutron energy of the efficiency curve (MeV)', nI('data-dk', 'Emax', d.Emax, 'min="1" max="250"')) +
    fld(
      'Monte Carlo histories per angle (pulse-height spectra, reconstructed spectrum, counting statistics)',
      nI('data-dk', 'Nph', d.Nph, 'min="10000" max="5000000" step="10000"')
    ) +
    fld('Energy averaging of the efficiency curve (± %; 0 = point values)', nI('data-dk', 'Ew', d.Ew, 'min="0" max="50"')) +
    '</div>' +
    '<details class="cw"><summary>Light-output curves (NRESP type, editable)</summary><p class="hint">Protons: L<sub>p</sub>(E) = a1·E − a2·[1 − exp(−a3·E<sup>a4</sup>)] (MeVee; Cecil et al. for the NE-213 type liquids, Madey / Verbinski for the plastic). Alphas: L<sub>α</sub> = kα·E<sup>kb</sup>. Carbon recoils: L<sub>C</sub> = c·E<sup>cb</sup>. Electrons are linear (1 MeVee = 1 MeV of electron energy).</p><div class="row3">' +
    [
      ['a1', 'a1'],
      ['a2', 'a2'],
      ['a3', 'a3'],
      ['a4', 'a4'],
      ['ka', 'α: kα'],
      ['kb', 'α: kb'],
      ['ca', 'C: c'],
      ['cb', 'C: cb']
    ]
      .map(([k, l]) => fld(l, nI('data-dlp', k, +(+lp[k]).toPrecision(5))))
      .join('') +
    '<label class="chk"><input type="checkbox" data-dlpc' +
    (lp.carbon === false ? '' : ' checked') +
    '><span>light from carbon recoils</span></label></div></details>'
  );
}
function ndDataStatus() {
  const D = ND.data;
  return (
    'Cross sections and angular distributions (embedded): <b>' +
    esc(D.nH) +
    '</b> (n+¹H, elastic) and <b>' +
    esc(D.nC) +
    '</b> (n+¹²C, ' +
    ndChannelsC(D.C).length +
    ' reaction channels, lumped above 20 MeV), ACE files of the ENDF/B-VIII.0 library (lib80x) at room temperature, 293.6 K.'
  );
}
function ndBuild() {
  const w = $('#ndWin');
  w.innerHTML =
    '<div class="mbox" style="width:min(1100px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:12px;position:sticky;top:-18px;z-index:5;background:var(--panel);padding:4px 0"><h3>Configure neutron detector</h3><button class="btn pri" data-ndclose id="ndDoneTop">Done</button></div>' +
    '<p class="hint">Cylindrical organic scintillator on the beam axis; the neutrons of the “Neutron counting” card come from a point source at the flight-path distance. The efficiency is a Monte Carlo of the neutron transport in the cell (multiple scattering included): ¹H(n,n), ¹²C(n,n), the ¹²C inelastic, (n,α), (n,p) and 3α break-up channels of the library (lumped above 20 MeV); the light of every recoil particle is added and an event counts when the total light exceeds the threshold.</p>' +
    '<div id="ndSketch">' +
    ndSketch() +
    '</div>' +
    '<div id="ndIn">' +
    ndInputs() +
    '</div>' +
    '<p class="hint" id="ndDataLbl"></p>' +
    '<div class="btns"><button class="btn pri" id="ndRun">Calculate efficiency</button><button class="btn" data-ndclose>Done</button><span class="hint" id="ndMsg"></span></div>' +
    '<div id="ndRes"></div></div>';
  $('#ndDataLbl').innerHTML = ndDataStatus();
  ndShow();
}
function ndOpen() {
  $('#ndWin').hidden = false;
  ndBuild();
}

// ---- running ----
const ndKey = () =>
  JSON.stringify([ndDet(), neuSet().dt, ndSet().N, ndSet().Emax, ndSet().Ew, ND.data.nH, ND.data.nC, (ND.data.H || {}).nMt, (ND.data.C || {}).nMt]);
// log grid with a spacing of 1.5 % (resolves the resonance structure of 12C), 80 points at least
const ndEnergies = Emax => {
  const em = Math.max(Emax, 0.2),
    n = Math.min(900, Math.max(80, Math.ceil(Math.log(em / 0.1) / 0.015) + 1));
  return Array.from({ length: n }, (_, i) => 0.1 * Math.pow(em / 0.1, i / (n - 1)));
};
async function ndCalc() {
  const my = ++ND.gen,
    msg = $('#ndMsg'),
    d = ndSet();
  ND.busy = true;
  if (msg) msg.textContent = 'Preparing the cross sections…';
  await new Promise(r => setTimeout(r, 20));
  const det = ndDet(),
    S = ndSetup(det, ND.data),
    E = ndEnergies(+d.Emax || 200),
    ew = d.Ew > 0 ? Math.log(1 + d.Ew / 100) : 0,
    eps = [],
    err = [],
    e1 = [];
  for (let i = 0; i < E.length; i++) {
    const r = ndEff(S, E[i], Math.max(1000, +d.N || 20000), 1000 + i, ew);
    eps.push(r.eps);
    err.push(r.err);
    e1.push(r.eps1);
    if (i % 25 === 24) {
      if (msg) msg.textContent = 'Monte Carlo: ' + (i + 1) + ' / ' + E.length + ' energies…';
      await new Promise(r2 => setTimeout(r2, 0));
      if (my !== ND.gen) return;
    }
  }
  ND.eff = { E, eps, err, e1, key: ndKey(), S };
  ND.busy = false;
  if (msg) msg.textContent = '';
  ndShow();
  if (NEU.res) neuShow();
}
function ndErrAt(E) {
  const f = ND.eff;
  if (!f) return 0;
  const x = f.E;
  if (E <= x[0]) return f.err[0];
  if (E >= x[x.length - 1]) return f.err[x.length - 1];
  let lo = 0,
    hi = x.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (x[m] > E) hi = m;
    else lo = m;
  }
  return f.err[lo] + ((f.err[hi] - f.err[lo]) * (E - x[lo])) / (x[hi] - x[lo]);
}
function ndEffAt(E) {
  const f = ND.eff;
  if (!f) return 0;
  const x = f.E;
  if (E <= x[0]) return f.eps[0];
  if (E >= x[x.length - 1]) return f.eps[x.length - 1];
  let lo = 0,
    hi = x.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (x[m] > E) hi = m;
    else lo = m;
  }
  return f.eps[lo] + ((f.eps[hi] - f.eps[lo]) * (E - x[lo])) / (x[hi] - x[lo]);
}

// ---- results ----
// Monte Carlo of the detector for the neutron spectrum Y (per energy bin of the card): pulse-height, reconstructed energy and time-of-flight histograms
function ndPhPrep(S, Y, r) {
  const dE = r.Eo / r.nbo;
  let mx = 0;
  Y.forEach(v => (mx = Math.max(mx, v)));
  let eLast = 0,
    eFirst = 0;
  Y.forEach((v, j) => {
    if (v > 1e-4 * mx) eLast = j;
  });
  while (eFirst < r.nbo - 1 && !(Y[eFirst] > 1e-4 * mx)) eFirst++;
  return { eFirst, eLast, Lmax: Math.max(0.2, ndLp((eLast + 1) * dE, S.lp) * 1.2) };
}
function ndPhRun(S, Y, r, Lp, st, N, seed, LmaxIn) {
  const dE = r.Eo / r.nbo,
    cdf = new Float64Array(r.nbo),
    Ec = Float64Array.from({ length: r.nbo }, (_, j) => (j + 0.5) * dE);
  let acc = 0;
  for (let j = 0; j < r.nbo; j++) {
    acc += Y[j];
    cdf[j] = acc;
  }
  for (let j = 0; j < r.nbo; j++) cdf[j] /= acc || 1;
  const { eFirst, eLast, Lmax: L0 } = ndPhPrep(S, Y, r),
    Lmax = LmaxIn || L0,
    nb = 200,
    rec = { nE: r.nbo, Eo: r.Eo, nT: 240, Lm: Lp, st, dE, t0: ndTof((eLast + 1) * dE, Lp) - 4 * st, t1: ndTof(Math.max(0.2, eFirst * dE), Lp) + 4 * st },
    pr = ndPulse(S, Ec, cdf, N, nb, Lmax, seed, rec);
  return { Lmax, nb, all: pr.all, sw: pr.sw, recE: pr.recE, recT: pr.recT, rec, sw2: pr.sw2, sd: pr.sd, sd2: pr.sd2, nd: pr.nd, N: pr.N, Lth: S.Lth };
}
const ndRatio = P => {
  const p = P.sw > 0 ? P.sd / P.sw : 0,
    v = P.sw > 0 ? Math.max(0, P.sd2 * (1 - 2 * p) + p * p * P.sw2) / (P.sw * P.sw) : 0;
  return { p, sp: Math.sqrt(v) };
};
const ND_YU = { rate: ['counts s⁻¹', 's⁻¹'], tot: ['counts', ''], perq: ['counts µC⁻¹', 'µC⁻¹'] };
function ndShow() {
  const box = $('#ndRes');
  if (!box) return;
  const f = ND.eff,
    d = ndSet();
  if (!f) {
    box.innerHTML = '<p class="hint">Set the parameters and press <b>Calculate efficiency</b>.</p>';
    return;
  }
  const stale = f.key !== ndKey(),
    s = neuSet(),
    r = NEU.res,
    qk = NEU_QTY[s.qty] ? s.qty : 'rate',
    q = NEU_QTY[qk];
  // efficiency curve
  const nz = f.E.map((e, i) => i).filter(i => f.eps[i] > 0);
  ND.raw.e = {
    series: [
      { name: 'Monte Carlo (all collisions)', type: 'line', x: nz.map(i => f.E[i]), y: nz.map(i => f.eps[i]), color: scol(0), w: 2 },
      { name: 'single ¹H(n,n) collision (analytic)', type: 'line', x: nz.map(i => f.E[i]), y: nz.map(i => f.e1[i]), color: scol(3), w: 1.5 }
    ],
    xlabel: 'Neutron energy (MeV)',
    ylabel: 'Efficiency (events above threshold / neutrons on the face)',
    title: 'Detector efficiency · ' + ND_MAT[d.mat].name + ' Ø' + d.D + ' × ' + d.T + ' cm · threshold ' + d.th + ' keVee',
    stats: [
      'L = ' +
        neuSet().L +
        ' m (point source on the axis) · statistical error of the Monte Carlo: ≤ ' +
        fmt(Math.max(...f.err), 2) +
        ' (1σ) · ' +
        (d.Ew > 0
          ? 'averaged over ± ' + d.Ew + ' % in neutron energy at each point (smooths the resonance structure of ¹²C)'
          : 'point values, no energy averaging')
    ]
  };
  let html = stale ? '<p class="warn">The parameters changed since the last calculation — press <b>Calculate efficiency</b> to update the plots.</p>' : '';
  html += '<h3 class="sub">Efficiency</h3><div id="ndP_e"></div>';
  if (!r || !r.cs || !(r.sc.okI || s.qty === 'perq')) {
    html +=
      '<p class="hint">The detected-neutron distribution and the pulse-height spectrum need a neutron spectrum: run the simulation and enter a beam current in the “Neutron counting” card (the angles come from that card).</p>';
    box.innerHTML = html;
    ndPlot('e');
    return;
  }
  const ia = Math.min(Math.max(+d.ia || 0, 0), r.angs.length - 1),
    Rr = f.S.geo.R / 100,
    Lm = f.S.geo.Ld / 100,
    om = 2 * Math.PI * (1 - Lm / Math.sqrt(Lm * Lm + Rr * Rr)),
    sc = r.sc[q.f],
    dE = r.Eo / r.nbo,
    Y = r.cs.ang[ia],
    Lp = +s.L > 0 ? +s.L : 25,
    st = Math.max(0, +s.dt || 0) / 2.3548,
    fold = st > 0,
    det = Float64Array.from(Y, (v, j) => {
      let e = 0;
      for (let k = 0; k < 8; k++) e += ndEffAt((j + (k + 0.5) / 8) * dE); // efficiency averaged over the bin
      return (v * sc * om * e) / 8;
    }),
    detF = fold ? neuFold(det, r.nbo, r.Eo, Lp, st) : det,
    tof = d.dist === 't',
    edges = new Float64Array(r.nbo + 1);
  for (let i = 0; i <= r.nbo; i++) edges[i] = i * dE;
  let last = 0,
    mx = 0;
  detF.forEach(v => (mx = Math.max(mx, v)));
  detF.forEach((v, j) => {
    if (v > 1e-4 * mx) last = j;
  });
  const un = ND_YU[qk];
  const dx = [],
    dy = [];
  if (tof)
    for (let j = r.nbo - 1; j >= 0; j--) {
      const k = neuKin((j + 0.5) * dE, Lp);
      dx.push(k.t);
      dy.push(detF[j] * k.dEdt);
    }
  const angTxt = 'θ = ' + fmt(r.angs[ia], 5) + '°, solid angle of the front face ' + fmt(om, 3) + ' sr';
  ND.raw.d = {
    series: [
      tof
        ? { name: 'detected neutrons', type: 'line', x: dx, y: dy, color: scol(0), w: 1.8 }
        : { name: 'detected neutrons', type: 'hist', edges, y: detF, color: scol(0) }
    ],
    xlabel: tof ? 'Time of flight over ' + fmt(Lp, 4) + ' m (ns)' : 'Neutron energy (MeV)',
    ylabel: 'Detected neutrons (' + un[0] + (tof ? ' ns⁻¹' : ' MeV⁻¹') + ')',
    title: 'Detected neutron ' + (tof ? 'time-of-flight' : 'energy') + ' distribution',
    stats: [angTxt + (fold ? ' · time resolution ' + fmt(+s.dt, 3) + ' ns FWHM folded in' : ' · ideal resolution')],
    tof,
    xmaxAuto: tof ? '' : String(fmt((last + 2) * dE, 4))
  };
  // pulse-height spectrum: neutron energies drawn from the spectrum at this angle
  let acc = 0;
  for (let j = 0; j < r.nbo; j++) acc += Y[j];
  const Ytot = acc * dE,
    ntot = Ytot * sc * om;
  const key = [f.key, ia, r.cs.W, r.Eo, s.L, s.dt, d.Nph].join('|');
  if (!ND.ph || ND.ph.key !== key) ND.ph = { key, ...ndPhRun(f.S, Y, r, Lp, st, Math.max(10000, +d.Nph || 300000), 4242) };
  const P = ND.ph,
    dL = P.Lmax / P.nb,
    ledges = new Float64Array(P.nb + 1);
  for (let i = 0; i <= P.nb; i++) ledges[i] = i * dL;
  const phy = Float64Array.from(P.all, v => (P.sw > 0 ? (v / P.sw) * ntot : 0) / dL),
    phThr = Float64Array.from(phy, (v, i) => ((i + 0.5) * dL > f.S.Lth ? v : 0));
  let nAbove = 0;
  phThr.forEach(v => (nAbove += v * dL));
  // statistics of the detected neutrons: ratio estimator p = Σw_det / Σw of the Monte Carlo (weights = solid-angle weights of the entry points)
  const pe = P.sw > 0 ? P.sd / P.sw : 0,
    vp = P.sw > 0 ? Math.max(0, P.sd2 * (1 - 2 * pe) + pe * pe * P.sw2) / (P.sw * P.sw) : 0,
    sp = Math.sqrt(vp),
    nDet = ntot * pe,
    sDet = ntot * sp,
    expInt = (() => {
      let t = 0;
      for (let j = 0; j < r.nbo; j++) t += det[j] * dE;
      return t;
    })();
  const rc = P.rec,
    rtof = d.dist === 't',
    rE = new Float64Array(rc.nE + 1),
    rT = new Float64Array(rc.nT + 1);
  for (let i = 0; i <= rc.nE; i++) rE[i] = (i * rc.Eo) / rc.nE;
  for (let i = 0; i <= rc.nT; i++) rT[i] = rc.t0 + (i * (rc.t1 - rc.t0)) / rc.nT;
  const dT = (rc.t1 - rc.t0) / rc.nT,
    rcY = rtof
      ? Float64Array.from(P.recT, v => (P.sw > 0 ? (v / P.sw) * ntot : 0) / dT)
      : Float64Array.from(P.recE, v => (P.sw > 0 ? (v / P.sw) * ntot : 0) / dE);
  ND.raw.r = {
    series: [{ name: 'reconstructed from the Monte Carlo events', type: 'hist', edges: rtof ? rT : rE, y: rcY, color: scol(1) }],
    xlabel: rtof ? 'Measured time of flight over ' + fmt(Lp, 4) + ' m (ns)' : 'Reconstructed neutron energy (MeV)',
    ylabel: 'Detected neutrons (' + un[0] + (rtof ? ' ns⁻¹' : ' MeV⁻¹') + ')',
    title: 'Reconstructed ' + (rtof ? 'time-of-flight' : 'energy') + ' spectrum (full-statistics Monte Carlo)',
    stats: [
      angTxt +
        ' · ' +
        P.N +
        ' neutron histories · events above ' +
        d.th +
        ' keVee · time of flight over ' +
        fmt(Lp, 4) +
        ' m smeared with ' +
        fmt(+s.dt, 3) +
        ' ns FWHM'
    ],
    tof: rtof,
    xmaxAuto: rtof ? '' : String(fmt((last + 2) * dE, 4))
  };
  ND.raw.h = {
    series: [
      { name: 'all events (light > 0)', type: 'hist', edges: ledges, y: phy, color: scol(0) },
      { name: 'above the threshold (' + d.th + ' keVee)', type: 'hist', edges: ledges, y: phThr, color: scol(3) }
    ],
    xlabel: 'Light output (MeVee)',
    ylabel: 'Counts (' + un[0] + ' MeVee⁻¹)',
    title: 'Pulse-height spectrum · ' + ND_MAT[d.mat].name,
    stats: [angTxt + ' · ' + fmt(nAbove, 4) + ' ' + un[0].replace(' s⁻¹', ' per s') + ' above the threshold (neutrons reaching the face: ' + fmt(ntot, 4) + ')']
  };
  const angSel =
      '<label class="pc-sel"><span>Angle</span><select data-dk="ia">' +
      r.angs.map((a, i) => '<option value="' + i + '"' + (i === ia ? ' selected' : '') + '>' + fmt(a, 6) + '°</option>').join('') +
      '</select></label>',
    distSel =
      '<label class="pc-sel"><span>Distribution</span><select data-dk="dist"><option value="E"' +
      (d.dist === 'E' ? ' selected' : '') +
      '>Neutron energy</option><option value="t"' +
      (d.dist === 't' ? ' selected' : '') +
      '>Time of flight</option></select></label>';
  const isTot = qk === 'tot',
    uq = un[0];
  html +=
    '<h3 class="sub">Statistics of the detected neutrons (θ = ' +
    fmt(r.angs[ia], 5) +
    '°)</h3><div class="tscroll"><table class="t"><tbody>' +
    '<tr><td>Neutrons reaching the front face</td><td>' +
    fmt(ntot, 5) +
    ' ' +
    esc(uq.replace('counts', 'n')) +
    '</td></tr>' +
    '<tr><td><b>Detected neutrons (above ' +
    d.th +
    ' keVee)</b></td><td><b>' +
    fmt(nDet, 5) +
    ' ± ' +
    fmt(sDet, 2) +
    ' ' +
    esc(uq) +
    '</b> (' +
    fmt((100 * sDet) / (nDet || 1), 2) +
    ' %, Monte Carlo, 1σ)</td></tr>' +
    '<tr><td>Mean efficiency for this spectrum</td><td>' +
    fmt(pe, 4) +
    ' ± ' +
    fmt(sp, 2) +
    '</td></tr>' +
    '<tr><td>Monte Carlo</td><td>' +
    P.nd +
    ' detected events out of ' +
    P.N +
    ' histories (N_eff = ' +
    fmt((P.sw * P.sw) / (P.sw2 || 1), 5) +
    ')</td></tr>' +
    '<tr><td>Same number from the efficiency curve (expected distribution)</td><td>' +
    fmt(expInt, 5) +
    ' ' +
    esc(uq) +
    '</td></tr>' +
    (isTot
      ? '<tr><td>Counting (Poisson) uncertainty of a real measurement</td><td>± ' +
        fmt(Math.sqrt(Math.max(nDet, 0)), 3) +
        ' counts (' +
        fmt(100 / Math.sqrt(Math.max(nDet, 1)), 2) +
        ' %)</td></tr>'
      : '') +
    '</tbody></table></div><p class="hint">The ± is the statistical error of the efficiency Monte Carlo (' +
    P.N +
    ' histories; it falls as 1/√N). The neutron spectrum of the card comes from a deterministic calculation and has no statistical error; its systematic uncertainty (cross sections, stopping power) and the efficiency model uncertainty are not included.</p>';
  html +=
    '<h3 class="sub">Detected neutrons: expected distribution</h3><div class="pcbar">' +
    angSel +
    distSel +
    '<span class="hint">Quantity and angles follow the “Neutron counting” card (' +
    esc(q.t) +
    '); counts = neutron spectrum × solid angle × efficiency(E), time resolution folded in. The next plot is the same distribution reconstructed event by event.</span></div><div id="ndP_d"></div>' +
    '<h3 class="sub">Reconstructed energy / time-of-flight spectrum</h3><div id="ndP_r"></div>' +
    '<h3 class="sub">Pulse-height spectrum in the detector</h3><div id="ndP_h"></div>';
  box.innerHTML = html;
  ndPlot('e');
  ndPlot('d');
  ndPlot('r');
  ndPlot('h');
}
function ndOpts(k) {
  const d = ND.raw[k],
    st = ND.pl[k],
    D = ndSet(),
    ylog = D[k + 'yl'] === 'log',
    xlog = D[k + 'xl'] === 'log',
    z = st.z;
  const series = d.series
    .map(se => {
      let o = { ...se, hidden: !!st.hidden[se.name] };
      if (se.type === 'line' && (ylog || xlog)) {
        const x = [],
          y = [];
        for (let i = 0; i < se.x.length; i++) {
          if ((ylog && !(se.y[i] > 0)) || (xlog && !(se.x[i] > 0))) continue;
          x.push(se.x[i]);
          y.push(se.y[i]);
        }
        o = { ...o, x, y };
      }
      return o;
    })
    .filter(se => se.type === 'hist' || se.x.length);
  return {
    series,
    xlog,
    ylog,
    xmin: z.xmin || (!xlog && (k === 'e' || k === 'h' || !d.tof) ? '0' : ''),
    xmax: z.xmax || (k === 'd' || k === 'r' ? d.xmaxAuto || '' : ''),
    ymin: z.ymin || (ylog && k === 'h' ? '' : ''),
    ymax: z.ymax || '',
    xlabel: d.xlabel,
    ylabel: d.ylabel,
    title: d.title,
    stats: d.stats
  };
}
function ndPlot(k) {
  const box = $('#ndP_' + k);
  if (!box || !ND.raw[k]) return;
  const st = ND.pl[k],
    D = ndSet(),
    o = ndOpts(k);
  if (!o.series.length) {
    box.innerHTML = '<p class="hint">Nothing to plot with these axes (no positive values).</p>';
    return;
  }
  st.o = o;
  const r = renderChart(o),
    names = o.series.map(se => [se.name, se.color]);
  box.innerHTML =
    gamRen(
      pcTop({
        y: ['data-gyl="x"', D[k + 'yl']],
        x: ['data-' + k + 'xl', D[k + 'xl']],
        reset: ['data-gzreset="x"', gamZoomed(st.z)],
        hint: 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click a legend entry to show/hide a curve'
      }),
      k
    ) +
    '<div class="chart" id="ndC_' +
    k +
    '">' +
    r.svg +
    '</div>' +
    gamRen(gamLegend(names, st.hidden, 'x'), k) +
    '<div class="readout" id="ndRo_' +
    k +
    '">Move the pointer over the plot to read values.</div>' +
    pcSave({
      attr: 'data-' + k + 'save',
      formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'rootbin', 'root'],
      msg: '<span class="hint" id="ndSv_' + k + '"></span>'
    });
  const sv = $('#ndC_' + k + ' svg', box),
    g = r.geom,
    ro = $('#ndRo_' + k, box);
  gamZoom(sv, g, 'z', st, () => ndPlot(k));
  sv.addEventListener('pointermove', ev => {
    const q = sv.getBoundingClientRect(),
      px = ((ev.clientX - q.left) / q.width) * g.W,
      py = ((ev.clientY - q.top) / q.height) * g.H;
    if (px < g.m.l || px > g.m.l + g.iw || py < g.m.t || py > g.m.t + g.ih) return;
    const X = pxData(g, px, py)[0],
      parts = [];
    o.series.forEach(se => {
      if (se.hidden) return;
      let v = null;
      if (se.type === 'hist') {
        const e = se.edges;
        if (X >= e[0] && X <= e[e.length - 1]) {
          let j = 0;
          while (j < e.length - 2 && e[j + 1] <= X) j++;
          v = se.y[j];
        }
      } else if (se.x.length > 1 && X >= se.x[0] && X <= se.x[se.x.length - 1]) {
        let j = 0;
        while (j < se.x.length - 2 && se.x[j + 1] < X) j++;
        v = se.y[j] + ((se.y[j + 1] - se.y[j]) * (X - se.x[j])) / (se.x[j + 1] - se.x[j] || 1);
      }
      if (v != null) parts.push(se.name.replace(/\s+\(.*/, '') + ' ' + fmt(v, 4));
    });
    ro.textContent =
      (k === 'h' ? 'L = ' : (k === 'd' || k === 'r') && D.dist === 't' ? 't = ' : 'E = ') +
      fmt(X, 5) +
      (k === 'h' ? ' MeVee' : (k === 'd' || k === 'r') && D.dist === 't' ? ' ns' : ' MeV') +
      '   |   ' +
      parts.join('   |   ');
  });
}
async function ndSave(k, kind) {
  const m = $('#ndSv_' + k),
    st = ND.pl[k];
  if (!st.o) return;
  if (m) m.textContent = 'Preparing…';
  const r = renderChart(st.o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  const base = 'cpt_neutron_detector_' + { e: 'efficiency', d: 'detected', r: 'reconstructed', h: 'pulse_height' }[k],
    ctx = { o: st.o, geom: r.geom, svg: div.querySelector('svg'), base, key: 'nd' + k };
  try {
    const it = await secBuild(ctx, kind);
    if (kind === 'rootbin') {
      const rn = it.name,
        rb = RootIO.file(rn, it.hs),
        inClaude = !!(window.claude && claude.use && (await claude.use('downloads')));
      if (inClaude) await saveFile(rn.replace(/\.root$/, '') + '_histograms.zip', zipStore(rn, rb), 'ROOT file', '#ndSv_' + k);
      else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', '#ndSv_' + k);
      if (m) m.textContent = 'Saved ' + rn;
    } else await saveFile(it.name, it.data, kind, '#ndSv_' + k);
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  }
}

// ---- events ----
function ndRefreshIn() {
  const s = $('#ndSketch');
  if (s) s.innerHTML = ndSketch();
}
document.addEventListener('click', ev => {
  if (ev.target.closest('[data-ndopen]')) return ndOpen();
  const w = ev.target.closest && ev.target.closest('#ndWin');
  if (ev.target.closest('[data-ndclose]') || (ev.target.id === 'ndWin' && typeof _mdownId !== 'undefined' && _mdownId === 'ndWin')) {
    $('#ndWin').hidden = true;
    if (typeof neuUpdate === 'function') neuUpdate();
    return;
  }
  if (!w) return;
  if (ev.target.closest('#ndRun')) return ndCalc();
  for (const k of ['e', 'd', 'r', 'h']) {
    if (ev.target.closest('[data-' + k + 'gzr]')) {
      ND.pl[k].z = {};
      return ndPlot(k);
    }
    const ga = ev.target.closest('[data-' + k + 'gall]');
    if (ga) {
      (ND.raw[k] ? ND.raw[k].series : []).forEach(se => (ND.pl[k].hidden[se.name] = ga.dataset.gstate === 'off'));
      return ndPlot(k);
    }
    const lg = ev.target.closest('[data-' + k + 'leg]');
    if (lg) {
      const n = lg.dataset[k + 'leg'];
      ND.pl[k].hidden[n] = !ND.pl[k].hidden[n];
      return ndPlot(k);
    }
    const sv = ev.target.closest('[data-' + k + 'save]');
    if (sv) return ndSave(k, sv.dataset[k + 'save']);
  }
});
document.addEventListener('input', ev => {
  const el = ev.target.closest && ev.target.closest('#ndWin');
  if (!el) return;
  const t = ev.target;
  if (t.matches('[data-dkn]')) {
    const s = neuSet();
    s[t.dataset.dkn] = t.value === '' ? '' : +t.value;
    saveLocal();
    ndRefreshIn();
  } else if (t.matches('[data-dk]') && t.type === 'number') {
    ndSet()[t.dataset.dk] = t.value === '' ? '' : +t.value;
    saveLocal();
    ndRefreshIn();
  } else if (t.matches('[data-dlp]')) {
    const d = ndSet();
    d.lpo = Object.assign({}, d.lpo, { [t.dataset.dlp]: +t.value });
    saveLocal();
  }
});
document.addEventListener('change', ev => {
  const t = ev.target;
  if (!t.closest || !t.closest('#ndWin')) return;
  if (t.matches('[data-dlpc]')) {
    const d = ndSet();
    d.lpo = Object.assign({}, d.lpo, { carbon: t.checked });
    saveLocal();
    return;
  }
  if (t.matches('[data-dk]') && t.tagName === 'SELECT') {
    const d = ndSet(),
      k = t.dataset.dk;
    d[k] = k === 'ia' ? +t.value : t.value;
    if (k === 'mat') {
      d.rho = ND_MAT[t.value].rho;
      d.hc = ND_MAT[t.value].hc;
      d.lpo = {};
      $('#ndIn').innerHTML = ndInputs();
    }
    saveLocal();
    ndRefreshIn();
    return ndShow();
  }
  for (const k of ['e', 'd', 'r', 'h']) {
    if (t.hasAttribute('data-' + k + 'gyl')) {
      ndSet()[k + 'yl'] = t.value;
      saveLocal();
      ND.pl[k].z = Object.assign({}, ND.pl[k].z, { ymin: '', ymax: '' });
      return ndPlot(k);
    }
    if (t.hasAttribute('data-' + k + 'xl')) {
      ndSet()[k + 'xl'] = t.value;
      saveLocal();
      ND.pl[k].z = Object.assign({}, ND.pl[k].z, { xmin: '', xmax: '' });
      return ndPlot(k);
    }
  }
});
