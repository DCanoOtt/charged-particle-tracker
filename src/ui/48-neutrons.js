/* ====== Neutron counting: neutron rate, total number and angular distribution of the neutrons emitted in the reactions of the library ======
   Same beam / irradiation inputs as "Irradiation and activation". The neutrons come from the MF 6 data of the library (secCompute of the secondary-particle module):
   LAW 7 tables (JENDL-5, laboratory) are interpolated in cos θ, in the incident energy and in E' with the interpolation laws of the file; Legendre / Kalbach-Mann / two-body
   data are evaluated at any angle. */
const NEU_DEF = {
  I: 1,
  IU: 'nA',
  q: 0,
  tirr: 1,
  tirrU: 'h',
  lay: 'all',
  L: 25,
  dt: 3,
  qty: 'rate',
  emin: 0,
  emu: 'MeV',
  nyl: 'log',
  nxl: 'lin',
  vyl: 'lin',
  vxl: 'lin'
};
// neutron energy threshold of the emitted neutrons: value in the chosen unit → MeV
const NEU_EU = { eV: 1e-6, keV: 1e-3, MeV: 1 };
const neuEmin = s => (+s.emin > 0 ? +s.emin * (NEU_EU[s.emu] || 1) : 0);
// the spectra of a neuCompute result without the neutrons below Emin (MeV): a bin is kept in proportion to the part of it above Emin; W is scaled with the kept fraction
function neuMask(c, Emin) {
  if (!(Emin > 0)) return c;
  const dE = c.Eo / c.nbo,
    w = Float64Array.from({ length: c.nbo }, (_, j) => Math.min(1, Math.max(0, ((j + 1) * dE - Emin) / dE))),
    cut = a => Float64Array.from(a, (v, j) => v * w[j]),
    sum = a => a.reduce((t, v) => t + v, 0);
  return { ...c, ang: c.ang.map(cut), tot: cut(c.tot), W: c.W * (sum(c.tot) > 0 ? sum(cut(c.tot)) / sum(c.tot) : 0) };
}
const NEU_MC2 = 939.56542,
  NEU_C = 0.299792458; // neutron rest energy (MeV), speed of light (m/ns)
// time of flight (ns) and |dE/dt| (MeV/ns) of a neutron of kinetic energy E (MeV) over L metres; sigma_E of a time resolution st (ns, 1 sigma)
function neuKin(E, L) {
  const g = 1 + E / NEU_MC2,
    b = Math.sqrt(1 - 1 / (g * g)),
    t = L / (b * NEU_C);
  return { t, dEdt: (NEU_MC2 * g * g * g * b * b) / t };
}
const neuSigE = (E, L, st) => neuKin(E, L).dEdt * st;
// Gaussian smearing of a spectrum y (per MeV, nbo bins up to Eo) with the energy-dependent width of the time resolution; the integral is kept
function neuFold(y, nbo, Eo, L, st) {
  const dE = Eo / nbo,
    out = new Float64Array(nbo);
  for (let i = 0; i < nbo; i++) {
    if (!(y[i] > 0)) continue;
    const sg = neuSigE((i + 0.5) * dE, L, st);
    if (!(sg > dE / 2)) {
      out[i] += y[i];
      continue;
    }
    const j0 = Math.max(0, Math.floor(i - (5 * sg) / dE)),
      j1 = Math.min(nbo - 1, Math.ceil(i + (5 * sg) / dE)),
      w = [];
    let sum = 0;
    for (let j = j0; j <= j1; j++) {
      const v = Math.exp((-0.5 * ((j - i) * dE) ** 2) / (sg * sg));
      w.push(v);
      sum += v;
    }
    for (let j = j0; j <= j1; j++) out[j] += (y[i] * w[j - j0]) / sum;
  }
  return out;
}
const NEU_ANG0 = [0, 30, 60, 90, 120, 150, 180]; // default angles (replaced by the nearest ENDF angles when the file has a grid)
const NEU_INT = { 1: 'histogram', 2: 'linear-linear', 3: 'linear-log', 4: 'log-linear', 5: 'log-log' };
const NEU_QTY = {
  rate: { f: 'rate', u: 'n s⁻¹ sr⁻¹', t: 'Neutron rate', d: 'n/s/sr' },
  tot: { f: 'tot', u: 'n sr⁻¹', t: 'Total number of neutrons produced', d: 'n/sr' },
  perq: { f: 'perq', u: 'n sr⁻¹ μC⁻¹', t: 'Neutrons per unit of charge', d: 'n/sr/μC' }
};
const NEU = {
  pl: { n: { hidden: {}, z: {} }, v: { hidden: {}, z: {} }, w: { hidden: {}, z: {} }, p: { hidden: {}, z: {} } },
  data: {},
  res: null,
  msg: '',
  gen: 0,
  tmr: null
};
function neuSet() {
  // filled in place: callers keep a reference while calling other functions that also call neuSet()
  const o = S.neu && typeof S.neu === 'object' ? S.neu : (S.neu = {});
  for (const k of Object.keys(NEU_DEF)) if (o[k] === undefined) o[k] = NEU_DEF[k];
  return o;
}
const neuDeg = mu => (Math.acos(Math.max(-1, Math.min(1, mu))) * 180) / Math.PI;
const neuR6 = v => Math.round(v * 1e6) / 1e6;

// neutron products with LAW 7 in the loaded library files: the cos θ grid of the file, its interpolation laws
function neuLaw7() {
  const mu = new Set(),
    ie = new Set(),
    im = new Set(),
    iw = new Set();
  let other = false;
  for (const d of ALIB.data.values()) {
    if (!d) continue;
    for (const m of d.order || []) {
      const ch = d.channels[m];
      if (!ch || !ch.primary) continue;
      for (const pr of ch.prods || [])
        if (pr.sp === 'n' && pr.dist) {
          if (pr.dist.law !== 7) {
            other = true;
            continue;
          }
          (pr.dist.int || []).forEach(v => ie.add(v));
          if (!(pr.dist.int || []).length) ie.add(2);
          for (const L of pr.dist.L) {
            for (const v of L.mu) mu.add(neuR6(v));
            (L.int || []).forEach(v => im.add(v));
            if (!(L.int || []).length) im.add(2);
            for (const T of L.T) (T.int || []).forEach(v => iw.add(v));
          }
        }
    }
  }
  return { nodes: Array.from(mu).sort((a, b) => a - b), ie: Array.from(ie), im: Array.from(im), iw: Array.from(iw), other };
}
const neuLaws = a => a.map(v => NEU_INT[v] + ' (INT ' + v + ')').join(', ') || '–';
// ENDF angles in degrees (decreasing cos θ → increasing θ)
function neuEndfAngles(L7) {
  return L7.nodes
    .map(m => neuR6(neuDeg(m)))
    .sort((a, b) => a - b)
    .filter((a, i, v) => !i || a !== v[i - 1]);
}
function neuAngles(L7) {
  const s = neuSet();
  if (Array.isArray(s.ang)) return s.ang.slice();
  const en = neuEndfAngles(L7);
  if (!en.length) return NEU_ANG0.slice();
  const out = [];
  for (const a of NEU_ANG0) {
    let b = en[0];
    for (const v of en) if (Math.abs(v - a) < Math.abs(b - a)) b = v;
    if (!out.includes(b)) out.push(b);
  }
  return out.sort((a, b) => a - b);
}
const neuIsNode = (L7, deg) => L7.nodes.some(m => Math.abs(m - Math.cos((deg * Math.PI) / 180)) < 2e-6);
// cos θ grid of the angular-distribution plot: the nodes of the file (thinned to 61), otherwise 41 uniform points
function neuGrid(L7) {
  let mu = L7.nodes.slice();
  if (mu.length > 61) mu = Array.from({ length: 61 }, (_, i) => mu[Math.round((i * (mu.length - 1)) / 60)]);
  if (mu.length < 2) mu = Array.from({ length: 41 }, (_, i) => -1 + i / 20);
  return mu;
}

function neuStatus() {
  const P = PJ(),
    R = RES;
  if (!alibF().size)
    return 'Load the cross-section library for ' + P.plural + ' first (left panel “Nuclear reaction data (ENDF-6)”: the TENDL folder, or ' + P.lib + ' files).';
  if (ALIB.busy) return 'Reading the cross-section files…';
  if (!ALIB.data.size || ![...ALIB.data.values()].some(Boolean)) return 'No cross-section file was found for the isotopes of the target.';
  if (!R) return 'Run the simulation first.';
  if (R.pj !== P.letter) return 'The results on screen were computed with another projectile — press Run for the ' + P.name + ' beam.';
  if (R.stale) return 'The setup (layers, beam or physics) has changed since this run — press Run to update the neutron results.';
  if (!R.Lh) return 'This run has no flux tally — run the simulation again.';
  return '';
}

// per-ion neutron yields (per sr and per MeV) of the layers at the given lab angles (degrees): { ang[i][j], tot[j], W, nUn, Eo, nbo }
function neuCompute(F, layers, angs, Eo, nbo) {
  const R = RES;
  R._neu = R._neu || new Map();
  const key = [ALIB.ver || 0, layers.join('.'), Eo, angs.join(',')].join('|');
  let c = R._neu.get(key);
  if (c) return c;
  c = { ang: angs.map(() => new Float64Array(nbo)), tot: new Float64Array(nbo), W: 0, nUn: 0, Eo, nbo };
  for (const li of layers) {
    const r = secCompute(F, li, 'n', ['sum'], nbo, Eo, angs, 'twobody');
    for (let j = 0; j < nbo; j++) {
      c.tot[j] += r.tot[j];
      for (let i = 0; i < angs.length; i++) c.ang[i][j] += r.ang[i][j];
    }
    c.W += r.W;
    c.nUn += r.nUnsup;
  }
  R._neu.set(key, c);
  if (R._neu.size > 40) R._neu.delete(R._neu.keys().next().value);
  return c;
}

// scale factors of the three quantities: per ion → n/s (rate), n over the irradiation, n per μC
function neuScale(s, P) {
  const I = +s.I,
    q = +s.q > 0 ? +s.q : P.Z,
    Tirr = +s.tirr * ACT_TU[s.tirrU || 'h'],
    rate = s.IU === 'pps' ? I : (I * ACT_IU[s.IU]) / (q * ACT_E);
  return { rate, tot: rate * Tirr, perq: 1e-6 / (q * ACT_E), okI: I > 0, okT: Tirr > 0, q };
}

function neuRender() {
  const card = $('#cNeu');
  if (!card) return;
  foldApply(card);
  const P = PJ();
  if (!P) {
    card.hidden = true;
    return;
  }
  card.hidden = false;
  const s = neuSet(),
    key = [P.letter, S.layers.length, S.layers.map(l => layerName(l)).join(',')].join('|');
  if (card.dataset.key !== key) {
    card.dataset.key = key;
    const nI = (k, extra) => '<input type="number" data-nk="' + k + '" step="any" ' + (extra || '') + ' value="' + (s[k] === undefined ? '' : s[k]) + '">',
      sel = (k, o) =>
        '<select data-nk="' +
        k +
        '">' +
        o.map(x => '<option value="' + x[0] + '"' + (String(s[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') +
        '</select>',
      fld = (l, h) => '<label class="fld"><span>' + l + '</span>' + h + '</label>';
    card.innerHTML =
      '<h2 data-fold tabindex="0" role="button" title="Click to collapse or expand">Neutron counting</h2><p class="hint">Neutrons emitted in the reactions of a constant beam of ' +
      P.plural +
      ' (flux and cross sections of the simulation above, neutron yields and angular distributions of the ' +
      P.sym +
      ' library): rate, total number and angular distribution at any emission angle in the laboratory frame.</p>' +
      '<div class="row3">' +
      fld('Beam current (constant)', nI('I', 'min="0"')) +
      fld(
        'Current unit',
        sel('IU', [
          ['pps', 'particles / s'],
          ['pA', 'pA'],
          ['nA', 'nA'],
          ['uA', 'μA'],
          ['mA', 'mA']
        ])
      ) +
      fld('Charge state q of the ions (e) — blank = fully stripped (' + P.Z + ')', nI('q', 'min="1" step="1" placeholder="' + P.Z + '"')) +
      '</div>' +
      '<div class="row3">' +
      fld('Irradiation time', nI('tirr', 'min="0"')) +
      fld('Unit', sel('tirrU', actUnits)) +
      fld('Layer', sel('lay', layOpts())) +
      '</div>' +
      '<div class="row3">' +
      '<div class="fld"><span>Flight path, time resolution, detector efficiency</span><div class="btns" style="margin:0"><button class="btn pri" data-ndopen style="font-weight:700;padding:8px 14px">Configure n-detector &amp; DAQ</button></div></div>' +
      fld(
        'Quantity in the plots',
        sel('qty', [
          ['rate', 'Neutron rate, n/s/sr'],
          ['tot', 'Total number of neutrons produced, n/sr'],
          ['perq', 'Neutrons per unit of charge, n/sr/μC (= n/s/sr per μA)']
        ])
      ) +
      '</div>' +
      '<div class="row3">' +
      fld('Neutron energy threshold (emitted neutrons; 0 = none)', nI('emin', 'min="0"')) +
      fld(
        'Energy unit',
        sel('emu', [
          ['eV', 'eV'],
          ['keV', 'keV'],
          ['MeV', 'MeV']
        ])
      ) +
      '</div>' +
      '<div class="row3" id="neuAngRow"></div>' +
      '<p class="hint" id="neuMsg"></p><div class="btns" id="neuChips" style="margin:4px 0"></div><div id="neuRes"></div>';
  } else
    card.querySelectorAll('[data-nk]').forEach(el => {
      if (document.activeElement !== el) el.value = s[el.dataset.nk] === undefined ? '' : s[el.dataset.nk];
    });
  neuUpdate();
}
const neuInt = (a, f, dE) => {
  let t = 0;
  for (let j = 0; j < a.length; j++) t += a[j] * f * dE;
  return t;
};
function neuUpdate() {
  const my = ++NEU.gen;
  clearTimeout(NEU.tmr);
  neuCalc(my).catch(e => {
    const b = $('#neuRes');
    if (b && my === NEU.gen) b.innerHTML = '<p class="warn">Calculation failed: ' + esc((e && e.message) || String(e)) + '</p>';
  });
}
// angle controls (ENDF angles menu, custom angle) and the chips of the selected angles
function neuAngUi(L7, angs) {
  const row = $('#neuAngRow'),
    ch = $('#neuChips');
  if (row && row.dataset.k !== L7.nodes.length + '|' + (L7.nodes[0] || '') + '|' + (L7.nodes[L7.nodes.length - 1] || '')) {
    row.dataset.k = L7.nodes.length + '|' + (L7.nodes[0] || '') + '|' + (L7.nodes[L7.nodes.length - 1] || '');
    const en = neuEndfAngles(L7);
    row.innerHTML =
      (en.length
        ? '<label class="fld"><span>ENDF angle (' +
          en.length +
          ' in the file)</span><select id="neuEndf">' +
          en.map(a => '<option value="' + a + '">' + fmt(a, 6) + '°   (cos θ = ' + fmt(Math.cos((a * Math.PI) / 180), 5) + ')</option>').join('') +
          '</select></label><div class="fld"><span>&nbsp;</span><div class="btns" style="margin:0"><button class="btn sm pri" id="neuAddE">Add this angle</button><button class="btn sm" id="neuAddAll">Add all ENDF angles</button></div></div>'
        : '<div class="fld"><span>ENDF angles</span><span class="hint" style="margin:0">The library has no tabulated (LAW 7) neutron angles: the angular distributions (Legendre, Kalbach-Mann, two-body kinematics) are evaluated at any angle.</span></div>') +
      '<label class="fld"><span>Any angle θ (degrees, 0 – 180)' +
      (en.length ? ' — interpolated with the laws of the file' : '') +
      '</span><input type="number" id="neuCust" min="0" max="180" step="any" placeholder="e.g. 45"></label>' +
      '<div class="fld"><span>&nbsp;</span><div class="btns" style="margin:0"><button class="btn sm pri" id="neuAddC">Add this angle</button><button class="btn sm" id="neuReset">Default angles</button></div></div>';
  }
  if (ch)
    ch.innerHTML = angs
      .map(
        (a, i) =>
          '<span class="chip" style="display:inline-flex;gap:6px;align-items:center;border:1px solid var(--line);border-radius:999px;padding:2px 4px 2px 10px"><i style="width:10px;height:10px;border-radius:50%;background:' +
          scol(i) +
          '"></i>θ = ' +
          fmt(a, 6) +
          '°' +
          (L7.nodes.length ? (neuIsNode(L7, a) ? ' · ENDF' : ' · interp.') : '') +
          '<button class="btn sm" data-neudel="' +
          i +
          '" title="Remove this angle">×</button></span>'
      )
      .join('');
}

async function neuCalc(my) {
  const box = $('#neuRes'),
    msg = $('#neuMsg'),
    P = PJ();
  if (!box || !P) return;
  const s = neuSet(),
    L7 = neuLaw7(),
    angs = neuAngles(L7);
  neuAngUi(L7, angs);
  const st = neuStatus();
  if (st) {
    NEU.res = null;
    if (msg) msg.textContent = st;
    box.innerHTML = '';
    return;
  }
  const sc = neuScale(s, P);
  if (s.qty !== 'perq' && !sc.okI) {
    if (msg) msg.textContent = 'Enter a beam current greater than zero.';
    box.innerHTML = '';
    return;
  }
  if (s.qty === 'tot' && !sc.okT) {
    if (msg) msg.textContent = 'Enter an irradiation time greater than zero.';
    box.innerHTML = '';
    return;
  }
  if (msg) msg.textContent = '';
  box.innerHTML = '<p class="hint">Computing the neutron yields…</p>';
  await new Promise(r => setTimeout(r, 20));
  if (my !== NEU.gen) return;
  const R = RES,
    F = alibFold(R),
    layers = s.lay === 'all' || +s.lay >= S.layers.length ? S.layers.map((l, i) => i) : [+s.lay],
    nbo = 160;
  let Eo = 0;
  for (const li of layers) Eo = Math.max(Eo, secBound(F, li, 'n', ['sum'], 'twobody'));
  if (!(Eo > 0)) {
    NEU.res = null;
    box.innerHTML = '<p class="hint">No neutron is emitted in the reactions of this library for the target and energies of the run.</p>';
    return;
  }
  const grid = neuGrid(L7),
    gdeg = grid.map(m => neuR6(neuDeg(m))),
    Emin = neuEmin(s),
    csAll = neuCompute(F, layers, angs, Eo, nbo), // without the energy threshold (the detector window does not use it)
    cs = neuMask(csAll, Emin),
    cg = neuMask(neuCompute(F, layers, gdeg, Eo, nbo), Emin);
  if (my !== NEU.gen) return;
  if (!(cs.W > 0) && Emin > 0 && neuCompute(F, layers, angs, Eo, nbo).W > 0) {
    NEU.res = null;
    box.innerHTML =
      '<p class="hint">No neutron is emitted above the energy threshold of ' + fmt(Emin, 4) + ' MeV (the neutrons reach at most ' + fmt(Eo, 4) + ' MeV).</p>';
    return;
  }
  if (!(cs.W > 0)) {
    NEU.res = null;
    box.innerHTML =
      '<p class="hint">No neutron is produced in this target by the ' +
      P.name +
      ' beam (no neutron channel with a distribution in the library files' +
      (cs.nUn ? '; ' + cs.nUn + ' channel(s) without angular data' : '') +
      ').</p>';
    return;
  }
  const dE = Eo / nbo,
    ysr = a => a.reduce((t, v) => t + v, 0) * dE, // per ion and sr
    meanE = a => {
      let n = 0,
        d = 0;
      for (let j = 0; j < nbo; j++) {
        n += a[j] * (j + 0.5) * dE;
        d += a[j];
      }
      return d > 0 ? n / d : 0;
    },
    tot4 = cs.tot.reduce((t, v) => t + v, 0) * dE;
  // 4π integral of the plotted distribution (trapezoids in cos θ) as a consistency check against the energy-integrated total
  const Yg = cg.ang.map(ysr);
  let i4 = 0;
  for (let i = 1; i < grid.length; i++) i4 += 0.5 * (Yg[i] + Yg[i - 1]) * (grid[i] - grid[i - 1]);
  i4 *= 2 * Math.PI;
  NEU.res = { angs, L7, sc, cs, csAll, cg, grid, gdeg, Yg, tot4, i4, Eo, nbo, layers, ysr, meanE, key: s.qty };
  neuShow();
}
// table + plots from NEU.res (cheap: redone when the quantity or the current changes)
function neuShow() {
  const box = $('#neuRes'),
    r = NEU.res,
    P = PJ(),
    s = neuSet();
  if (!box || !r || !P) return;
  const { angs, L7, sc, cs, tot4, i4 } = r,
    rows = angs
      .map((a, i) => {
        const y = r.ysr(cs.ang[i]);
        return (
          '<tr><td>' +
          fmt(a, 6) +
          '</td><td>' +
          fmt(Math.cos((a * Math.PI) / 180), 5) +
          '</td><td>' +
          (L7.nodes.length ? (neuIsNode(L7, a) ? 'ENDF node' : 'interpolated') : 'analytic') +
          '</td><td>' +
          (sc.okI ? fmt(y * sc.rate, 4) : '–') +
          '</td><td>' +
          (sc.okI && sc.okT ? fmt(y * sc.tot, 4) : '–') +
          '</td><td>' +
          fmt(y * sc.perq, 4) +
          '</td><td>' +
          fmt(y, 4) +
          '</td><td>' +
          fmt(r.meanE(cs.ang[i]), 4) +
          '</td></tr>'
        );
      })
      .join('');
  neuPrep();
  box.innerHTML =
    '<p class="hint">' +
    (neuEmin(neuSet()) > 0
      ? '<b>Only the neutrons above the energy threshold of ' +
        fmt(neuEmin(neuSet()), 4) +
        ' MeV are included</b> in the spectra, yields, mean energies and detector results. '
      : '') +
    '<b>Per steradian at the angles chosen above.</b> Rate = n/s/sr; total = n/sr over the whole irradiation (rate × irradiation time); per charge = n/sr/μC (identical to n/s/sr per μA of beam). The “4π” row integrates over all directions.</p>' +
    '<div class="tscroll"><table class="t"><tr><th>θ lab (°)</th><th>cos θ</th><th>Angular data</th><th>n s⁻¹ sr⁻¹</th><th>n sr⁻¹</th><th>n sr⁻¹ μC⁻¹</th><th>n sr⁻¹ per ion</th><th>⟨E<sub>n</sub>⟩ (MeV)</th></tr>' +
    rows +
    '<tr><td><b>4π, all angles</b></td><td></td><td></td><td>' +
    (sc.okI ? fmt(tot4 * sc.rate, 4) + ' n/s' : '–') +
    '</td><td>' +
    (sc.okI && sc.okT ? fmt(tot4 * sc.tot, 4) + ' n' : '–') +
    '</td><td>' +
    fmt(tot4 * sc.perq, 4) +
    ' n/μC</td><td>' +
    fmt(tot4, 4) +
    ' n/ion</td><td>' +
    fmt(r.meanE(cs.tot), 4) +
    '</td></tr></table></div>' +
    tblBar('neu_sr', 'Neutron_yield_per_steradian') +
    '<p class="hint">Angle-integrated check: ∫ Y(θ) dω of the plotted distribution = ' +
    fmt(i4, 4) +
    ' n/ion, total from the energy spectrum = ' +
    fmt(tot4, 4) +
    ' n/ion (' +
    (tot4 > 0 ? fmt((100 * (i4 - tot4)) / tot4, 2) : '0') +
    ' %).' +
    (cs.nUn ? ' <b>' + cs.nUn + ' neutron-producing channel(s) have no angular/energy distribution in the file and are not counted.</b>' : '') +
    '</p>' +
    '<h3 class="sub">Neutron spectrum at the selected angles</h3><div id="neuP_n"></div>' +
    '<h3 class="sub">Angular distribution of the emitted neutrons</h3><div id="neuP_v"></div>' +
    neuDetHtml() +
    '<p class="hint" id="neuInfo"></p>';
  const info = $('#neuInfo');
  if (info) {
    const kinds = secKinds(RES ? alibFold(RES) : { rows: [] }, r.layers[0], 'n', ['sum']);
    info.textContent =
      'Angular data used: ' +
      (kinds || '–') +
      '.' +
      (L7.nodes.length
        ? ' LAW 7 tables of the file: ' +
          L7.nodes.length +
          ' cos θ values; interpolation laws read from the file — incident energy: ' +
          neuLaws(L7.ie) +
          '; cos θ: ' +
          neuLaws(L7.im) +
          '; emitted energy: ' +
          neuLaws(L7.iw) +
          '. An angle between two ENDF angles is interpolated with the cos θ law.'
        : '') +
      (L7.other && L7.nodes.length ? ' Other neutron channels of the library use Legendre / Kalbach-Mann / two-body data.' : '');
  }
  neuPlot('n');
  neuPlot('v');
  neuPlot('w');
  if (NEU.det) {
    neuDetRefresh();
    neuPhCalc();
  }
}
// raw series of the two plots for the current quantity
function neuPrep() {
  const r = NEU.res,
    s = neuSet(),
    P = PJ(),
    q = NEU_QTY[s.qty] || NEU_QTY.rate,
    f = r.sc[q.f],
    beam = RES ? RES.beamName + ' ' + fmt(RES.E0nom, 4) + ' MeV' : P.name,
    mats = (() => {
      const out = r.layers.map(i => (S.layers.length > 1 ? 'layer ' + (i + 1) + ' ' : '') + layerName(S.layers[i], i));
      return out.join('; ');
    })(),
    edges = new Float64Array(r.nbo + 1);
  for (let i = 0; i <= r.nbo; i++) edges[i] = (i * r.Eo) / r.nbo;
  let last = 0,
    mx = 0;
  r.cs.ang.forEach(a => a.forEach(v => (mx = Math.max(mx, v))));
  r.cs.ang.forEach(a => {
    for (let j = 0; j < r.nbo; j++) if (a[j] > 1e-4 * mx) last = Math.max(last, j);
  });
  const n = {
    series: r.cs.ang.map((a, i) => ({ name: 'θ = ' + fmt(r.angs[i], 5) + '°', type: 'hist', edges, y: Float64Array.from(a, v => v * f), color: scol(i) })),
    xlabel: 'Neutron energy (MeV)',
    ylabel: q.t + ' (' + q.u + ' MeV⁻¹)',
    title: 'Neutron spectra emitted from the target · ' + beam,
    stats: [
      'Target: ' + mats + ' · ideal spectra at the target (no flight path, no detector, no resolution)',
      'Bins: ' +
        r.nbo +
        ' of ' +
        fmt(r.Eo / r.nbo, 4) +
        ' MeV (0 – ' +
        fmt(r.Eo, 4) +
        ' MeV) · integral of the spectrum: ' +
        r.cs.ang.map((a, i) => 'θ = ' + fmt(r.angs[i], 5) + '°: ' + fmt(neuInt(a, f, r.Eo / r.nbo), 4)).join(' · ') +
        ' ' +
        q.u
    ],
    xmaxAuto: String(fmt(((last + 2) * r.Eo) / r.nbo, 4))
  };
  const xs = r.grid,
    ys = r.Yg.map(v => v * f);
  const m = {
    series: [
      { name: 'integrated over all neutron energies', type: 'line', x: Array.from(xs), y: ys, color: scol(0), w: 1.8 },
      {
        name: 'selected angles',
        type: 'pts',
        x: r.angs.map(a => Math.cos((a * Math.PI) / 180)),
        y: r.angs.map((a, i) => r.ysr(r.cs.ang[i]) * f),
        color: scol(3),
        r: 4
      }
    ],
    xlabel: 'cos θ (laboratory)',
    ylabel: q.t + ' (' + q.u + ')',
    title: 'Angular distribution of the neutrons · ' + beam,
    stats: ['Target: ' + mats + ' · forward = cos θ = +1']
  };
  NEU.data = { n, v: m };
  NEU.det = null;
  if (typeof ND !== 'undefined' && ND.eff && ND.eff.S) neuPrepDet(r, s, f, q, beam, edges, last);
}
// spectra of the neutrons detected by the configured detector: ideal spectrum x solid angle of the front face x intrinsic efficiency(E)
// (events above the threshold / neutrons entering the front face), time resolution folded in; the integral is the number of detected neutrons
function neuPrepDet(r, s, f, q, beam, edges, last) {
  const g = ND.eff.S.geo,
    Rr = g.R / 100,
    Lm = g.Ld / 100,
    om = ndOmega(g.rb > 0 ? g.rb / 100 : Rr, Lm),
    dE = r.Eo / r.nbo,
    Lp = +s.L > 0 ? +s.L : 25,
    st = Math.max(0, +s.dt || 0) / 2.3548,
    tof = ndSet().dist === 't',
    un = ND_YU[NEU_QTY[s.qty] ? s.qty : 'rate'],
    rows = [],
    series = [];
  r.cs.ang.forEach((a, i) => {
    let var2 = 0;
    const det = Float64Array.from(a, (v, j) => {
      let e = 0,
        er = 0;
      for (let k = 0; k < 8; k++) {
        e += ndEffAt((j + (k + 0.5) / 8) * dE);
        er += ndErrAt((j + (k + 0.5) / 8) * dE);
      }
      const y = (v * f * om * e) / 8;
      var2 += Math.pow((v * f * om * dE * er) / 8, 2);
      return y;
    });
    const detF = st > 0 ? neuFold(det, r.nbo, r.Eo, Lp, st) : det;
    let tot = 0;
    for (let j = 0; j < r.nbo; j++) tot += det[j] * dE;
    rows.push({ a: r.angs[i], tot, err: Math.sqrt(var2) });
    if (tof) {
      const x = [],
        y = [];
      for (let j = r.nbo - 1; j >= 0; j--) {
        const k = neuKin((j + 0.5) * dE, Lp);
        x.push(k.t);
        y.push(detF[j] * k.dEdt);
      }
      series.push({ name: 'θ = ' + fmt(r.angs[i], 5) + '°', type: 'line', x, y, color: scol(i), w: 1.8 });
    } else series.push({ name: 'θ = ' + fmt(r.angs[i], 5) + '°', type: 'hist', edges, y: detF, color: scol(i) });
  });
  NEU.det = { rows, om, un, tof };
  NEU.data.w = {
    series,
    tof,
    xlabel: tof ? 'Time of flight over ' + fmt(Lp, 4) + ' m (ns)' : 'Neutron energy (MeV)',
    ylabel: 'Detected neutrons (' + un[0] + (tof ? ' ns⁻¹' : ' MeV⁻¹') + ')',
    title: 'Neutrons detected · ' + ndLabel() + ' · ' + beam,
    stats: [
      'Bins: ' +
        r.nbo +
        ' of ' +
        fmt(dE, 4) +
        ' MeV' +
        (tof ? ' (converted to time of flight)' : '') +
        ' · integral of the spectrum (detected neutrons): ' +
        rows.map(x => 'θ = ' + fmt(x.a, 5) + '°: ' + fmt(x.tot, 4)).join(' · ') +
        ' ' +
        un[0],
      'spectrum × solid angle of the front face (' +
        fmt(om, 3) +
        ' sr) × intrinsic efficiency · time resolution ' +
        fmt(+s.dt, 3) +
        ' ns folded in · the integral is the number of detected neutrons'
    ],
    xmaxAuto: tof ? '' : String(fmt(((last + 2) * r.Eo) / r.nbo, 4))
  };
}
// counting statistics (all angles) and pulse-height spectra from the Monte Carlo of the detector at every angle
function neuStatHtml() {
  const r = NEU.res,
    D = NEU.det,
    sc = r.sc,
    ph = NEU.phc && NEU.phc.list,
    T = sc.okI && sc.okT ? sc.tot / sc.rate : 0,
    f2 = (v, e) => (v === null ? '…' : fmt(v, e)),
    fq = NEU_QTY[neuSet().qty] ? r.sc[NEU_QTY[neuSet().qty].f] : sc.rate;
  const rows = r.angs
    .map((a, i) => {
      const ysr = r.ysr(r.cs.ang[i]),
        nc = ysr * sc.tot * D.om,
        nr = ysr * sc.rate * D.om,
        R = ph ? ndRatio(ph[i]) : null,
        eps = R ? R.p : D.rows[i].tot / (ysr * fq * D.om || 1),
        N = nc * eps,
        Rt = nr * eps,
        mc = R && R.p > 0 ? R.sp / R.p : null,
        c = sc.okI && sc.okT,
        sN = Math.sqrt(Math.max(N, 0)),
        pu = sc.okI ? ndPile(Rt, ndSet().tau, ndSet().pu) : null;
      return (
        '<tr><td>' +
        fmt(a, 6) +
        '</td><td>' +
        fmt(ysr, 5) +
        '</td><td>' +
        (c ? fmt(nc, 5) : '–') +
        '</td><td>' +
        (c ? fmt(N, 5) : '–') +
        '</td><td>' +
        (c ? '± ' + fmt(sN, 3) : '–') +
        '</td><td>' +
        (c ? '± ' + f2(mc === null ? null : N * mc, 3) : '–') +
        '</td><td>' +
        (sc.okI ? fmt(Rt, 5) : '–') +
        '</td><td>' +
        (c ? '± ' + fmt(sN / T, 3) : '–') +
        '</td><td>' +
        (sc.okI ? '± ' + f2(mc === null ? null : Rt * mc, 3) : '–') +
        '</td><td' +
        (pu && pu.loss > 0.1 ? ' style="color:var(--bad);font-weight:700"' : pu && pu.loss > 0.01 ? ' style="color:var(--warn)"' : '') +
        '>' +
        (pu ? (ndSet().pu === 'off' ? 'off' : fmt(100 * pu.loss, 3) + ' % (nτ = ' + fmt(pu.mu, 2) + ')') : '–') +
        '</td><td>' +
        (pu ? fmt(pu.rec, 5) : '–') +
        '</td><td>' +
        fmt((eps * D.om) / (4 * Math.PI), 4) +
        (R ? ' ± ' + fmt((R.sp * D.om) / (4 * Math.PI), 2) : '') +
        '</td></tr>'
      );
    })
    .join('');
  return (
    '<p class="hint">' +
    'All numbers are for the detector placed at the flight-path distance L = ' +
    fmt(+neuSet().L > 0 ? +neuSet().L : 25, 4) +
    ' m from the neutron source (' +
    (ND.eff.S.ppac
      ? 'U-238 deposit on the axis, source in the direction θ, neutron beam spot Ø' + fmt(ND.eff.S.geo.R * 2, 4) + ' cm on the target'
      : 'front face on the axis, source in the direction θ, ' +
        (ND.eff.S.geo.rb > 0 && ND.eff.S.geo.rb < ND.eff.S.geo.R - 1e-9 ? 'neutron beam Ø' + fmt(ND.eff.S.geo.rb * 2, 4) + ' cm on the face, ' : '') +
        'cylinder ' +
        fmt(ND.eff.S.geo.R * 2, 4) +
        ' × ' +
        fmt(ND.eff.S.geo.T, 4) +
        ' cm') +
    '): the neutrons reaching the detector are the yield per sr at θ × the solid angle of the ' +
    (ND.eff.S.ppac ? 'neutron beam spot' : ND.eff.S.geo.rb > 0 && ND.eff.S.geo.rb < ND.eff.S.geo.R - 1e-9 ? 'neutron beam on the front face' : 'front face') +
    ', ω = ' +
    fmt(D.om, 4) +
    ' sr. ' +
    (sc.okI ? 'Ion rate = ' + fmt(sc.rate, 5) + ' s⁻¹ (beam current / (charge state × e)); count rate = ion rate × Y(θ) × ω × intrinsic efficiency. ' : '') +
    (T
      ? 'Irradiation time T = ' + fmt(T, 5) + ' s. '
      : 'Enter the beam current and the irradiation time to get the counts and the statistical uncertainty of the rate. ') +
    'Statistical = Poisson √N of the detected counts (√N / T for the rate): the uncertainty a real measurement would have. Monte Carlo = uncertainty of the efficiency calculation (' +
    (ph ? ph[0].N : '…') +
    ' histories per angle), reduced by running more histories. Pile-up level = fraction of the detected events lost at the rate n above (pulse width τ = ' +
    fmt(+ndSet().tau, 4) +
    ' ns, ' +
    (ndSet().pu === 'non' ? 'non-paralysable: n/(1 + nτ)' : ndSet().pu === 'off' ? 'model off' : 'paralysable, ±τ window: n·exp(−2nτ)') +
    '; set in the DAQ section of the detector window; red above 10 %, orange above 1 %); the count rate after pile-up is n × (1 − level). Total efficiency ε<sub>tot</sub> = detected counts / neutrons emitted into 4π at the same yield = intrinsic efficiency × ω / 4π, with ω = ' +
    fmt(D.om, 4) +
    ' sr the solid angle of the front face (mean intrinsic efficiencies: ' +
    r.angs.map((a, i) => 'θ = ' + fmt(a, 5) + '°: ' + fmt(ph ? ndRatio(ph[i]).p : 0, 4)).join(', ') +
    ').</p>' +
    '<div class="tscroll"><table class="t"><tr><th>θ lab (°)</th><th>Yield Y(θ) (n sr⁻¹ per ion)</th><th>Neutrons reaching the detector</th><th>Detected counts N</th><th>± statistical</th><th>± Monte Carlo</th><th>Count rate (s⁻¹)</th><th>± statistical</th><th>± Monte Carlo</th><th>Pile-up level</th><th>Count rate after pile-up (s⁻¹)</th><th>Total efficiency ε<sub>tot</sub></th></tr>' +
    rows +
    '</table></div>' +
    tblBar('neu_stat', 'Neutron_counting_statistics')
  );
}
function neuDetRefresh() {
  const box = $('#neuStat');
  if (!box || !NEU.res || !NEU.det) return;
  box.innerHTML = neuStatHtml();
  const ph = NEU.phc && NEU.phc.list,
    r = NEU.res;
  if (!ph || ph[0].ppac) return; // the PPAC has no pulse height
  const D = NEU.det,
    s = neuSet(),
    qk = NEU_QTY[s.qty] ? s.qty : 'rate',
    f = r.sc[NEU_QTY[qk].f],
    series = [],
    ints = [];
  ph.forEach((P, i) => {
    const dL = P.Lmax / P.nb,
      ntot = r.ysr(r.cs.ang[i]) * f * D.om,
      ledges = new Float64Array(P.nb + 1);
    for (let j = 0; j <= P.nb; j++) ledges[j] = j * dL;
    const y = Float64Array.from(P.all, v => (P.sw > 0 ? (v / P.sw) * ntot : 0) / dL);
    let tot = 0,
      thr = 0;
    y.forEach((v, j) => {
      tot += v * dL;
      if ((j + 0.5) * dL > P.Lth) thr += v * dL;
    });
    ints.push('θ = ' + fmt(r.angs[i], 5) + '°: ' + fmt(thr, 4) + ' above the threshold (' + fmt(tot, 4) + ' all events)');
    series.push({ name: 'θ = ' + fmt(r.angs[i], 5) + '°', type: 'hist', edges: ledges, y, color: scol(i) });
  });
  NEU.data.p = {
    series,
    xlabel: 'Light output (MeVee)',
    ylabel: 'Detected neutrons (' + D.un[0] + ' MeVee⁻¹)',
    title: 'Pulse-height spectra · ' + ndLabel(),
    stats: [
      'Bins: ' +
        ph[0].nb +
        ' of ' +
        fmt(ph[0].Lmax / ph[0].nb, 4) +
        ' MeVee (0 – ' +
        fmt(ph[0].Lmax, 4) +
        ' MeVee) · integral of the spectrum: ' +
        ints.join(' · '),
      'same vertical quantity as the detected-neutron plot · threshold ' + ndSet().th + ' keVee · ' + ph[0].N + ' histories per angle'
    ]
  };
  neuPlot('p');
}
let neuPhGen = 0;
async function neuPhCalc() {
  const r = NEU.res;
  if (!r || !NEU.det || !ND.eff || !ND.eff.S) return;
  const s = neuSet(),
    key = [ND.eff.key, r.Eo, r.nbo, r.angs.join(), r.cs.W, s.L, s.dt, r.sc.rate, ndSet().Nph].join('|');
  if (NEU.phc && NEU.phc.key === key) return neuDetRefresh();
  NEU.phc = null;
  const my = ++neuPhGen,
    Lp = +s.L > 0 ? +s.L : 25,
    st = Math.max(0, +s.dt || 0) / 2.3548,
    S = ND.eff.S,
    list = [];
  await new Promise(res => setTimeout(res, 30));
  for (let i = 0; i < r.angs.length; i++) {
    if (my !== neuPhGen) return;
    list.push(ndPhRun(S, r.cs.ang[i], r, Lp, st, Math.max(10000, +ndSet().Nph || 300000), 777 + i));
    await new Promise(res => setTimeout(res, 0));
  }
  if (my !== neuPhGen) return;
  NEU.phc = { key, list };
  neuDetRefresh();
}
function neuDetHtml() {
  const D = NEU.det,
    d = typeof ndSet === 'function' ? ndSet() : null;
  let h = '<h3 class="sub">Neutrons detected by the configured detector</h3>';
  if (!D)
    return (
      h +
      '<p class="hint">Open <i>Configure n-detector &amp; DAQ</i> and press <i>Calculate efficiency</i>: this plot then shows the spectra of each angle corrected by the intrinsic efficiency and the solid angle of the detector, so that their integral is the total number of neutrons detected.</p>'
    );
  h +=
    '<div class="pcbar"><label class="pc-sel"><span>Distribution</span><select data-neudist><option value="E"' +
    (d.dist === 'E' ? ' selected' : '') +
    '>Neutron energy</option><option value="t"' +
    (d.dist === 't' ? ' selected' : '') +
    '>Time of flight</option></select></label>' +
    '<span class="hint">Intrinsic efficiency = ' +
    (ND.eff && ND.eff.S && ND.eff.S.ppac
      ? 'detected fission events / neutrons on the target spot'
      : 'events above the threshold / neutrons entering the front face') +
    ' (the solid angle is applied separately: ' +
    fmt(D.om, 3) +
    ' sr).</span></div><div id="neuP_w"></div>' +
    (ND.eff && ND.eff.S && ND.eff.S.ppac
      ? ''
      : '<h3 class="sub">Pulse-height spectra of the detected neutrons</h3><div id="neuP_p"><p class="hint">Monte Carlo of the pulse height at every angle running…</p></div>') +
    '<h3 class="sub">Counting statistics</h3><div id="neuStat"></div>';
  return h;
}
function neuOpts(k) {
  const d = NEU.data[k],
    st = NEU.pl[k],
    s = neuSet(),
    ylog = s[k + 'yl'] === 'log',
    xlog = s[k + 'xl'] === 'log',
    z = st.z;
  const series = d.series
    .map(se => {
      let o = { ...se, hidden: !!st.hidden[se.name] };
      if (se.type !== 'hist' && (ylog || xlog)) {
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
    xmin: z.xmin || ((k === 'n' || k === 'w' || k === 'p') && !xlog && !d.tof ? '0' : ''),
    xmax: z.xmax || (k === 'n' || k === 'w' ? d.xmaxAuto || '' : ''),
    ymin: z.ymin || '',
    ymax: z.ymax || '',
    xlabel: d.xlabel,
    ylabel: d.ylabel,
    title: d.title,
    stats: d.stats
  };
}
function neuPlot(k) {
  const box = $('#neuP_' + k);
  if (!box || !NEU.data[k]) return;
  const st = NEU.pl[k],
    s = neuSet(),
    o = neuOpts(k);
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
        y: ['data-gyl="x"', s[k + 'yl']],
        x: ['data-' + k + 'xl', s[k + 'xl']],
        reset: ['data-gzreset="x"', gamZoomed(st.z)],
        hint: 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click a legend entry to show/hide a curve'
      }),
      k
    ) +
    '<div class="chart" id="neuC_' +
    k +
    '">' +
    r.svg +
    '</div>' +
    gamRen(gamLegend(names, st.hidden, 'x'), k) +
    '<div class="readout" id="neuRo_' +
    k +
    '">Move the pointer over the plot to read values.</div>' +
    pcSave({
      attr: 'data-' + k + 'save',
      formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'rootbin', 'root'],
      msg: '<span class="hint" id="neuSv_' + k + '"></span>'
    });
  const sv = $('#neuC_' + k + ' svg', box),
    g = r.geom,
    ro = $('#neuRo_' + k, box);
  gamZoom(sv, g, 'z', st, () => neuPlot(k));
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
      } else if (se.type === 'line' && se.x.length > 1) {
        let j = 0;
        while (j < se.x.length - 2 && se.x[j + 1] < X) j++;
        if (X >= se.x[0] && X <= se.x[se.x.length - 1]) v = se.y[j] + ((se.y[j + 1] - se.y[j]) * (X - se.x[j])) / (se.x[j + 1] - se.x[j] || 1);
      }
      if (v != null) parts.push(se.name + ' ' + fmt(v, 4));
    });
    const dd = NEU.data[k];
    ro.textContent =
      (k === 'v' ? 'cos θ = ' : k === 'p' ? 'L = ' : dd && dd.tof ? 't = ' : 'E = ') +
      fmt(X, 5) +
      (k === 'v' ? '' : k === 'p' ? ' MeVee' : dd && dd.tof ? ' ns' : ' MeV') +
      '   |   ' +
      parts.join('   |   ');
  });
}
async function neuSave(k, kind) {
  const m = $('#neuSv_' + k),
    st = NEU.pl[k];
  if (!st.o) return;
  if (m) m.textContent = 'Preparing…';
  const r = renderChart(st.o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  const base = 'cpt_neutrons_' + (k === 'n' ? 'spectrum_' : k === 'w' ? 'detected_' : k === 'p' ? 'pulseheight_' : 'angular_') + (PJ() ? PJ().letter : ''),
    ctx = { o: st.o, geom: r.geom, svg: div.querySelector('svg'), base, key: 'neu' + k };
  try {
    const it = await secBuild(ctx, kind);
    if (kind === 'rootbin') {
      const rn = it.name,
        rb = RootIO.file(rn, it.hs),
        inClaude = !!(window.claude && claude.use && (await claude.use('downloads')));
      if (inClaude) await saveFile(rn.replace(/\.root$/, '') + '_histograms.zip', zipStore(rn, rb), 'ROOT file', '#neuSv_' + k);
      else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', '#neuSv_' + k);
      if (m) m.textContent = 'Saved ' + rn + (inClaude ? ' (inside a zip)' : '') + '.';
    } else await saveFile(it.name, it.data, kind, '#neuSv_' + k);
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  }
}
function neuAddAngle(deg) {
  const msg = $('#neuMsg');
  if (!(deg >= 0 && deg <= 180)) {
    if (msg) msg.textContent = 'The angle must be between 0° and 180°.';
    return false;
  }
  const s = neuSet(),
    angs = neuAngles(neuLaw7()),
    a = neuR6(deg);
  if (angs.some(v => Math.abs(v - a) < 1e-6)) {
    if (msg) msg.textContent = 'θ = ' + fmt(a, 6) + '° is already in the plot.';
    return false;
  }
  angs.push(a);
  s.ang = angs.sort((x, y) => x - y);
  saveLocal();
  return true;
}
document.addEventListener('click', ev => {
  const c = ev.target.closest && ev.target.closest('#cNeu');
  if (!c) return;
  const t = ev.target;
  let redo = false;
  if (t.closest('#neuAddE')) redo = neuAddAngle(+$('#neuEndf').value);
  else if (t.closest('#neuAddAll')) {
    const s = neuSet(),
      en = neuEndfAngles(neuLaw7());
    s.ang = Array.from(new Set(neuAngles(neuLaw7()).concat(en))).sort((a, b) => a - b);
    saveLocal();
    redo = true;
  } else if (t.closest('#neuAddC')) {
    const v = $('#neuCust').value;
    if (v === '') return;
    if (neuAddAngle(+v)) $('#neuCust').value = '';
    redo = true;
  } else if (t.closest('#neuReset')) {
    delete neuSet().ang;
    saveLocal();
    redo = true;
  } else if (t.closest('[data-neudel]')) {
    const s = neuSet(),
      a = neuAngles(neuLaw7());
    a.splice(+t.closest('[data-neudel]').dataset.neudel, 1);
    s.ang = a;
    saveLocal();
    redo = true;
  }
  if (redo) return neuUpdate();
  for (const k of ['n', 'v', 'w', 'p']) {
    if (t.closest('[data-' + k + 'gzr]')) {
      NEU.pl[k].z = {};
      return neuPlot(k);
    }
    const ga = t.closest('[data-' + k + 'gall]');
    if (ga) {
      (NEU.data[k] ? NEU.data[k].series : []).forEach(se => (NEU.pl[k].hidden[se.name] = ga.dataset.gstate === 'off'));
      return neuPlot(k);
    }
    const lg = t.closest('[data-' + k + 'leg]');
    if (lg) {
      const n = lg.dataset[k + 'leg'];
      NEU.pl[k].hidden[n] = !NEU.pl[k].hidden[n];
      return neuPlot(k);
    }
    const sv = t.closest('[data-' + k + 'save]');
    if (sv) return neuSave(k, sv.dataset[k + 'save']);
  }
});
document.addEventListener('keydown', ev => {
  if (ev.key === 'Enter' && ev.target && ev.target.id === 'neuCust') $('#neuAddC').click();
});
document.addEventListener('input', ev => {
  const el = ev.target.closest && ev.target.closest('#cNeu [data-nk]');
  if (!el || el.tagName === 'SELECT') return;
  const k = el.dataset.nk,
    s = neuSet();
  s[k] = el.value === '' ? '' : +el.value;
  if (k === 'q' && !(s.q > 0)) s.q = 0;
  saveLocal();
  clearTimeout(NEU.tmr);
  NEU.tmr = setTimeout(neuUpdate, 250);
});
document.addEventListener('change', ev => {
  const t = ev.target;
  if (!t.closest || !t.closest('#cNeu')) return;
  if (t.matches('[data-nk]') && t.tagName === 'SELECT') {
    const k = t.dataset.nk,
      s = neuSet();
    s[k] = t.value;
    saveLocal();
    if (k === 'qty') {
      NEU.pl.n.z = {};
      NEU.pl.v.z = {};
      return neuUpdate();
    }
    return neuUpdate();
  }
  if (t.matches('[data-neudist]')) {
    ndSet().dist = t.value;
    saveLocal();
    return neuUpdate();
  }
  for (const k of ['n', 'v', 'w', 'p']) {
    if (t.hasAttribute('data-' + k + 'gyl')) {
      neuSet()[k + 'yl'] = t.value;
      saveLocal();
      NEU.pl[k].z = Object.assign({}, NEU.pl[k].z, { ymin: '', ymax: '' });
      return neuPlot(k);
    }
    if (t.hasAttribute('data-' + k + 'xl')) {
      neuSet()[k + 'xl'] = t.value;
      saveLocal();
      NEU.pl[k].z = Object.assign({}, NEU.pl[k].z, { xmin: '', xmax: '' });
      return neuPlot(k);
    }
  }
});
