/* ====== tab "Activation": residual nuclei produced by a constant beam current during an irradiation, their half-lives and activities ====== */
const ACT_E = 1.602176634e-19;
const ACT_TU = { s: 1, min: 60, h: 3600, d: 86400, y: 31556926 };
const ACT_IU = { pps: 0, pA: 1e-12, nA: 1e-9, uA: 1e-6, mA: 1e-3 };
let ACT_LAST = null;
function actNuc(za) {
  return (ELSYM[Math.floor(za / 1000)] || '?') + '-' + (za % 1000);
}
function actHl(za) {
  // {T (s) | null, st: 'stable' | 'rad' | 'unk', mode}
  if (HL_STABLE.has(za)) return { T: null, st: 'stable', mode: '' };
  if (HL_T[za] > 0) return { T: HL_T[za], st: 'rad', mode: HL_MODE[za] || '' };
  return { T: null, st: HL_MODE[za] ? 'rad' : 'unk', mode: HL_MODE[za] || '' };
}
function fmtHl(t) {
  if (!(t > 0)) return '–';
  const u = [
    [1e-6, 'ns', 1e-9],
    [1e-3, 'µs', 1e-6],
    [1, 'ms', 1e-3],
    [60, 's', 1],
    [3600, 'min', 60],
    [86400, 'h', 3600],
    [31556926 * 0.999, 'd', 86400],
    [Infinity, 'y', 31556926]
  ];
  for (const [lim, nm, f] of u) if (t < lim) return fmt(t / f, 4) + ' ' + nm;
  return '';
}
function actTime(s) {
  return fmtHl(s) || '0';
}
const actEv = (R, lam, T) => (lam > 0 ? R * (1 - Math.exp(-lam * T)) : 0); // activity (Bq) at time T of constant production R (atoms/s)

const ACT_DEF = {
  I: 1,
  IU: 'nA',
  q: 0,
  tirr: 1,
  tirrU: 'h',
  tcool: 0,
  tcoolU: 'h',
  lay: 'all',
  show: 'rad',
  top: 8,
  ylog: 'log',
  xlog: 'lin',
  qty: 'act',
  bplus: 0,
  src: 'int',
  fb: true,
  chain: true
};
function irrSet() {
  return (S.irr = Object.assign({}, ACT_DEF, S.irr || {}));
}
const actUnits = [
  ['s', 's'],
  ['min', 'min'],
  ['h', 'hours'],
  ['d', 'days'],
  ['y', 'years']
];
const DEC_MSG0 =
  'The internal database is used. Evaluated decay sub-library files (e.g. ENDF/B-VIII.0 decay, JEFF-3.3/4.0 RDD, JENDL-5 DDF) can be loaded from a folder or as individual files.';
let ACT_CTX = null,
  ACT_TMR = null,
  ACT_GEN = 0;
const ACTP = new Map(),
  ACTST = new Map(); // plot data and view state (hidden curves, zoom) per layer

// the panel "Irradiation and activation" below the plots: inputs, summary, one plot per layer and the table
function actRender() {
  const card = $('#cAct');
  if (!card) return;
  const P = PJ();
  foldApply(card);
  if (!P) {
    card.hidden = true;
    return;
  }
  card.hidden = false;
  const s = irrSet(),
    key =
      P.letter +
      '|' +
      S.layers.length +
      '|' +
      S.layers.map(l => layerName(l)).join(',') +
      '|' +
      (DEC.busy ? 'b' : DEC.nuc) +
      '|' +
      DEC.prog +
      '|' +
      (ENS.busy ? 'b' : ENS.nuc) +
      '|' +
      ENS.prog;
  if (card.dataset.key !== key) {
    // controls are built once, so typing is never interrupted
    card.dataset.key = key;
    const nI = (k, extra) => '<input type="number" data-ak="' + k + '" step="any" ' + (extra || '') + ' value="' + (s[k] === undefined ? '' : s[k]) + '">',
      sel = (k, o) =>
        '<select data-ak="' +
        k +
        '">' +
        o.map(x => '<option value="' + x[0] + '"' + (String(s[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') +
        '</select>';
    const fld = (l, h) => '<label class="fld"><span>' + l + '</span>' + h + '</label>',
      chk = (k, t) => '<label class="chk"><input type="checkbox" data-ak="' + k + '"' + (s[k] ? ' checked' : '') + '><span>' + t + '</span></label>';
    card.innerHTML =
      '<h2 data-fold tabindex="0" role="button" title="Click to collapse or expand">Irradiation and activation</h2><p class="hint">Residual nuclei produced by a constant beam of ' +
      P.plural +
      ' in all reactions of the ' +
      P.sym +
      ' library (flux and cross sections of the simulation above), their decay chains (grow-in and decay solved with the Bateman equations), half-lives, production rates, atoms and activities.</p>' +
      '<div class="row3">' +
      fld('Beam current (constant)', nI('I', 'min="0"')) +
      fld(
        'Current unit',
        sel('IU', [
          ['pps', 'particles / s'],
          ['pA', 'pA'],
          ['nA', 'nA'],
          ['uA', 'µA'],
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
      fld('Cooling time after the irradiation', nI('tcool', 'min="0"')) +
      fld('Unit', sel('tcoolU', actUnits)) +
      fld(
        'Table',
        sel('show', [
          ['rad', 'radioactive nuclei'],
          ['all', 'all nuclei, stable too']
        ])
      ) +
      '</div>' +
      '<div class="row3">' +
      fld('Nuclei in each plot', nI('top', 'min="1" max="30" step="1"')) +
      '</div>' +
      '<div class="row3">' +
      fld(
        'Quantity in the plots',
        sel('qty', [
          ['act', 'Decays per second = activity (radioactive nuclei)'],
          ['gam', 'Gamma emission rate, γ/s (needs ENDF-6 decay data)'],
          ['bet', 'Beta emission rate, β/s (β⁻ and β⁺ particles)'],
          ['xr', 'X-ray emission rate, X/s (needs ENDF-6 decay data with X-ray records)'],
          ['num', 'Number of atoms (all isotopes produced, stable ones too)']
        ])
      ) +
      '</div>' +
      '<div class="row3">' +
      fld(
        'Decay data',
        sel('src', [
          ['int', 'Internal database (ICRP-107 / NUBASE2020 with isomers)'],
          ['file', 'Uploaded ENDF-6 decay data' + (DEC.nuc ? ' (' + DEC.nuc.toLocaleString('en') + ' nuclides)' : ' — none loaded')],
          ['ens', 'Uploaded ENSDF mass-chain files' + (ENS.nuc ? ' (' + ENS.idx.size + ' mass chains)' : ' — none loaded')]
        ])
      ) +
      '<div class="fld"><span id="decLbl">Uploaded decay data</span><div class="btns" style="margin:0"><button class="btn sm pri" id="decChoose">Choose folder…</button><button class="btn sm" id="decFiles">Choose files…</button></div></div>' +
      '<div class="fld">' +
      chk('chain', 'Follow the decay chains (grow-in of daughters)') +
      chk('fb', 'Uploaded data: internal values for nuclides missing in the files') +
      '</div></div>' +
      '<div class="row3">' +
      fld(
        'β⁺ share of the combined EC/β⁺ decays (%) — used by the β emission rate where the decay data do not give it (the ENSDF files do, per nuclide)',
        nI('bplus', 'min="0" max="100" step="1"')
      ) +
      '</div>' +
      '<p class="hint" id="decMsg">' +
      esc(ENS.busy ? ENS.prog : DEC.busy ? DEC.prog : s.src === 'ens' && ENS.msg ? ENS.msg : DEC.msg || DEC_MSG0) +
      '</p>' +
      '<input type="file" id="decInput" webkitdirectory directory multiple hidden><input type="file" id="decInputF" multiple hidden><div id="actRes"></div>';
  } else
    card.querySelectorAll('[data-ak]').forEach(el => {
      if (document.activeElement !== el) {
        if (el.type === 'checkbox') el.checked = !!s[el.dataset.ak];
        else el.value = s[el.dataset.ak] === undefined ? '' : s[el.dataset.ak];
      }
    });
  actUpdate();
}
