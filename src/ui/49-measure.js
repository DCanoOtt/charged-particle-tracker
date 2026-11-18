/* ====== simulated decay measurement: irradiation + cooling counted with an HPGe detector ====== */
const BKG_PEAKS =
  '72.8 0.015 Pb Kα2\n75.0 0.025 Pb Kα1\n84.5 0.009 Pb Kβ1\n87.3 0.003 Pb Kβ2\n238.6 0.020 212Pb\n295.2 0.006 214Pb\n351.9 0.012 214Pb\n511.0 0.030 annihilation\n583.2 0.008 208Tl\n609.3 0.012 214Bi\n911.2 0.008 228Ac\n969.0 0.005 228Ac\n1120.3 0.003 214Bi\n1460.8 0.018 40K\n1764.5 0.003 214Bi\n2614.5 0.007 208Tl';
const BKG_EC = 25,
  BKG_ETH = 1.5,
  BKG_E1 = 3500,
  BKG_R0 = 5,
  BKG_SOFT = 150,
  BKG_HARD = 1500,
  BKG_FSOFT = 0.7,
  MEAS_MAXN = 1000;
const MEAS_DEF = {
  tU: 'h',
  tirr: 1,
  tcool: 4,
  tcU: 'h',
  ppu: 6,
  ppuC: 6,
  xU: 'auto',
  roiK: 1.5,
  bmode: 'cont',
  brate: BKG_R0,
  bL: 500,
  I: '',
  IU: '',
  q: '',
  bEc: BKG_EC,
  bpk: BKG_PEAKS,
  noise: true,
  seed: 1,
  newSeed: true,
  pile: true,
  infer: true,
  minI: 1,
  qty: 'cnt',
  xbar: true,
  layF: 'all',
  emax: 0,
  dI: 0,
  dE: 0,
  markSel: true,
  showTot: true,
  showNet: true,
  err: true,
  th: true,
  thPile: true,
  thBkg: true,
  ylog: 'lin',
  byl: 'log',
  bxl: 'lin',
  syl: 'log',
  sel: {}
};
let MEAS = { res: null, hidden: {}, z: {}, sol0: null, list: [], busy: false };
const MEAS_VER = 2;
function measSet() {
  // persistent settings object (mutated in place, so references stay valid)
  const st = irrSet(),
    t0 = +st.tirr > 0 ? +st.tirr : 1;
  if (st.meas && st.meas.ver !== MEAS_VER) st.meas = null; // settings saved by an older version of the window are discarded once (they can hold values that no longer mean the same)
  if (!st.meas) st.meas = { ver: MEAS_VER, tU: st.tirrU || 'h', tirr: t0, tcool: +st.tcool > 0 ? +st.tcool : 4 * t0 };
  const m = st.meas;
  if (m.tcU === undefined) m.tcU = m.tU || 'h';
  if (m.ppuC === undefined) m.ppuC = m.ppu !== undefined ? m.ppu : MEAS_DEF.ppu;
  for (const k in MEAS_DEF) if (m[k] === undefined) m[k] = k === 'sel' ? {} : MEAS_DEF[k];
  if (!m.sel || typeof m.sel !== 'object') m.sel = {};
  return m;
}
const MEAS_UN = { s: 's', min: 'min', h: 'h', d: 'd', y: 'y' };
function measAx() {
  // time unit of the plot axis: chosen, or the unit of the longer phase
  const m = measSet(),
    R = MEAS.res;
  let k = m.xU;
  if (!k || k === 'auto' || !ACT_TU[k]) k = R ? (R.Tc > R.Tirr ? R.kC : R.kI) : m.tU || 'h';
  return { k, f: ACT_TU[k] || 3600, u: MEAS_UN[k] || k };
}
const measPhi = z => 0.5 * gamErfc(-z / Math.SQRT2);
function measRng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function measPois(mu, rnd) {
  if (!(mu > 0)) return 0;
  if (mu < 40) {
    const L = Math.exp(-mu);
    let k = 0,
      p = 1;
    do {
      k++;
      p *= rnd();
    } while (p > L);
    return k - 1;
  }
  const z = Math.sqrt(-2 * Math.log(rnd() || 1e-12)) * Math.cos(2 * Math.PI * rnd());
  return Math.max(0, Math.round(mu + Math.sqrt(mu) * z));
}
function measSigma(E, g) {
  return Math.max(gamFwhm(E, g), 1e-6) / 2.354820045;
}
function measAB(g) {
  const s1 = measSigma(1000, g);
  return g.fmode === 'sqrt' ? { a: 0, b: (s1 * s1) / 1000 } : { a: s1 * s1, b: 0 };
} // σ²(E) = a + bE
function measParsePeaks(t) {
  const out = [];
  for (const ln of String(t || '').split(/\r?\n/)) {
    const w = ln.trim().split(/\s+/);
    const E = +w[0],
      R = +w[1];
    if (E > 0 && R >= 0 && isFinite(E) && isFinite(R)) out.push({ E, R, name: w.slice(2).join(' ') });
  }
  return out;
}
// Relative low-energy efficiency of the background: transmission of the end cap, dead layer and
// contacts (photoelectric absorption ∝ E⁻³) times a soft electronic threshold; → 0 for E → 0.
function measBkgEff(E, Ec) {
  if (!(E > 0)) return 0;
  return Math.exp(-Math.pow(Ec / E, 3)) * (1 - Math.exp(-E / BKG_ETH));
}
function measBkg(m, g) {
  // background model → { spec(E) cps/keV, roi(lo, hi) cps, total cps, ... }, defined from 0 to BKG_E1
  const Rt = Math.max(0, +m.brate || 0),
    Ec = Math.max(0.1, +m.bEc || BKG_EC),
    none = { spec: () => 0, roi: () => 0, total: 0, comps: [], Ec };
  if (m.bmode === 'none' || Rt <= 0) return none;
  const eps = E => measBkgEff(E, Ec),
    H = 0.25,
    N = Math.ceil(BKG_E1 / H);
  // cumulative integral of eps(E)·exp(−E/L) over 0…BKG_E1 (trapezoid, 0.25 keV)
  const cumExp = L => {
    const c = new Float64Array(N + 1);
    let prev = 0;
    for (let k = 1; k <= N; k++) {
      const E = k * H,
        f = eps(E) * Math.exp(-E / L);
      c[k] = c[k - 1] + 0.5 * H * (prev + f);
      prev = f;
    }
    return c;
  };
  const at = (c, E) => {
    const x = Math.min(Math.max(E, 0), BKG_E1) / H,
      k = Math.min(Math.floor(x), N - 1);
    return c[k] + (c[k + 1] - c[k]) * (x - k);
  };
  const expInt = (c, lo, hi) => (hi > lo ? at(c, hi) - at(c, lo) : 0);
  if (m.bmode === 'exp') {
    const L = Math.max(1, +m.bL || 500),
      c = cumExp(L),
      A = Rt / c[N];
    return {
      spec: E => (E >= 0 && E <= BKG_E1 ? A * eps(E) * Math.exp(-E / L) : 0),
      roi: (lo, hi) => A * expInt(c, lo, hi),
      total: Rt,
      comps: [{ A, L }],
      peaks: [],
      Ec
    };
  }
  const sc = Rt / BKG_R0,
    pk = measParsePeaks(m.bpk).map(p => ({ E: p.E, I: p.R * sc * eps(p.E), name: p.name })), // detected rate = intrinsic × efficiency
    inR = pk.filter(p => p.E <= BKG_E1),
    Rpk = inR.reduce((a, p) => a + p.I, 0),
    Rc = Math.max(0, Rt - Rpk),
    outR = pk.filter(p => p.E > BKG_E1).reduce((a, p) => a + p.I, 0),
    cumS = cumExp(BKG_SOFT),
    cumH = cumExp(BKG_HARD);
  const comps = [
    { L: BKG_SOFT, A: (BKG_FSOFT * Rc) / cumS[N], c: cumS },
    { L: BKG_HARD, A: ((1 - BKG_FSOFT) * Rc) / cumH[N], c: cumH }
  ];
  return {
    spec: E => {
      let v = 0;
      if (E >= 0 && E <= BKG_E1) for (const c of comps) v += c.A * eps(E) * Math.exp(-E / c.L);
      for (const p of pk) {
        const sg = measSigma(p.E, g);
        v += (p.I / (sg * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * ((E - p.E) / sg) ** 2);
      }
      return v;
    },
    roi: (lo, hi) => {
      let v = 0;
      for (const c of comps) v += c.A * expInt(c.c, lo, hi);
      for (const p of pk) {
        const sg = measSigma(p.E, g);
        v += p.I * (measPhi((hi - p.E) / sg) - measPhi((lo - p.E) / sg));
      }
      return v;
    },
    total: Rt + outR,
    comps,
    peaks: pk,
    Ec
  };
}
function measNm(c) {
  return (MEAS.keys && MEAS.keys.length > 1 ? 'L' + (c.li + 1) + ' · ' : '') + dkNuc(c.id);
}
function measLayName(li) {
  return 'Layer ' + (li + 1) + (RES && RES.names && RES.names[li] ? ' · ' + RES.names[li] : '');
}
function measSelScope() {
  const m = measSet();
  return m.layF === 'all' || m.layF === undefined || m.layF === '' ? null : +m.layF;
}
function measLayers() {
  const sc = measSelScope();
  return (MEAS.keys || []).filter(k => (sc === null || k === sc) && ACTP.get(k) && ACTP.get(k).Q);
} // layers seen by the detector: all, or the chosen one
function measSelCombos(m) {
  // [{ i, id, E, I, key }] from the selection and the loaded list
  const out = [];
  const sc = measSelScope();
  for (const n of MEAS.list)
    for (const l of n.lines) {
      if (sc !== null && n.li !== sc) continue;
      const key = n.li + '|' + n.id + '|' + l[0];
      if (m.sel[key]) out.push({ li: n.li, i: n.i, id: n.id, e: n.e, E: l[0], I: l[1], key });
    }
  return out;
}
function measSelChanged() {
  // a result that no longer matches the selected lines is removed (never show points of lines that are not selected)
  const R = MEAS.res;
  if (!R || MEAS.busy) return;
  const cur = new Set(measSelCombos(measSet()).map(c => c.li + '|' + c.id + '|' + c.E)),
    old = R.C.map(c => c.li + '|' + c.id + '|' + c.E);
  if (cur.size === old.length && old.every(k => cur.has(k))) return;
  MEAS.res = null;
  const o = $('#measOut');
  if (o) o.innerHTML = '';
  const msg = $('#measMsg');
  if (msg)
    msg.textContent = cur.size
      ? 'The selection of γ lines has changed: press “Simulate Experiment” to update the plots.'
      : 'No γ line selected: select at least one isotope and γ line (button “Select isotopes and γ lines”).';
}
function measBeam(m) {
  // beam of the measurement: empty fields = the values of the activation panel; → { I, IU, q, rate (primaries/s), f (factor on the production rates) }
  const a = irrSet(),
    I = m.I !== '' && m.I != null && +m.I >= 0 ? +m.I : +a.I,
    IU = m.IU || a.IU || 'pps',
    q = +m.q > 0 ? +m.q : +a.q > 0 ? +a.q : null,
    Z = PJ().Z,
    rate = IU === 'pps' ? I : (I * ACT_IU[IU]) / ((q || Z) * ACT_E);
  return { I, IU, q, rate };
}
function measProd(li, m) {
  // production rates of the layer scaled to the beam current of the measurement window (everything is linear in the beam)
  const Q = ACTP.get(li).Q,
    b = measBeam(m);
  if (!(Q.rate > 0) || !(b.rate >= 0) || Math.abs(b.rate / Q.rate - 1) < 1e-12) return Q.prod;
  const f = b.rate / Q.rate;
  return new Map(Array.from(Q.prod, ([k, v]) => [k, v * f]));
}
function measEffInfo(g) {
  const eff = gamEffFn(g),
    et = gamEffTotFn(g),
    ge = g.effMode === 'ge',
    manual = !ge || g.ftpMode === 'man',
    f = Math.max(1, +g.ftp || 1);
  return { eff, ppd: e => (manual ? gamDetPerDecay(e, eff) * f : gamDetPerDecay(e, et)) }; // pulses per decay: Ge model → total efficiency; otherwise peak efficiency × pulses per full-energy count
}

async function measRun() {
  if (MEAS.busy) {
    MEAS.again = true;
    return;
  }
  const m = measSet(),
    g = gamSet(),
    msg = $('#measMsg'),
    cmb0 = measSelCombos(m),
    keys = measLayers(); // the detector counts the chosen layer, or all layers
  if (!keys.length) {
    if (msg)
      msg.textContent = cmb0.length
        ? 'Run the irradiation and activation calculation first.'
        : 'Select at least one isotope and γ line (button “Select isotopes and γ lines”).';
    return;
  }
  const combos = measSelCombos(m);
  if (!combos.length) {
    if (msg) msg.textContent = 'Select at least one isotope and γ line (button “Select isotopes and γ lines”).';
    return;
  }
  const kI = ACT_TU[m.tU] ? m.tU : 'h',
    kC = ACT_TU[m.tcU] ? m.tcU : kI,
    fI0 = ACT_TU[kI],
    fC0 = ACT_TU[kC],
    Tirr = Math.max(0, +m.tirr) * fI0,
    Tc = Math.max(0, +m.tcool) * fC0;
  if (!(Tirr > 0)) {
    if (msg) msg.textContent = 'The irradiation time must be positive.';
    return;
  }
  const T = Tirr + Tc,
    pI = Math.max(0.01, +m.ppu || 1),
    pC = Math.max(0.01, +m.ppuC || 1); // each phase has its own time unit and its own number of points per unit
  let NI = Math.max(1, Math.ceil((Tirr / fI0) * pI - 1e-9)),
    NC = Tc > 0 ? Math.max(1, Math.ceil((Tc / fC0) * pC - 1e-9)) : 0,
    capped = false;
  if (NI + NC > MEAS_MAXN) {
    const sc = MEAS_MAXN / (NI + NC);
    NI = Math.max(1, Math.floor(NI * sc));
    NC = Tc > 0 ? Math.max(1, Math.floor(NC * sc)) : 0;
    capped = true;
  }
  const N = NI + NC,
    dtI = Tirr / NI,
    dtC = NC ? Tc / NC : 0,
    edges = [];
  for (let k = 0; k <= NI; k++) edges.push(k === NI ? Tirr : k * dtI);
  for (let k = 1; k <= NC; k++) edges.push(k === NC ? T : Tirr + k * dtC);
  MEAS.busy = true;
  if (msg) msg.textContent = 'Solving the Bateman equations for ' + N + ' time bins…';
  await new Promise(r => setTimeout(r, 20));
  try {
    const sols = {};
    for (const li of keys) sols[li] = await dkSolve(Object.assign({}, ACTP.get(li).Q, { prod: measProd(li, m), Tirr, Tc, extra: edges.slice(1, N) }));
    const prodM = Object.fromEntries(keys.map(li => [li, measProd(li, m)]));
    const seed = m.newSeed ? (Math.random() * 4e9) >>> 0 : +m.seed >>> 0,
      rnd = measRng(seed),
      dI = Math.max(0, +m.dI || 0) / 100,
      dE = Math.max(0, +m.dE || 0) / 100,
      K = Math.max(0.1, +m.roiK || 1.5);
    if (m.newSeed) {
      const s = irrSet().meas;
      s.seed = seed;
      const inp = $('#measWin [data-mk=seed]');
      if (inp) inp.value = seed;
      m.seed = seed;
    }
    const { eff, ppd } = measEffInfo(g),
      bk = measBkg(m, g),
      tau = Math.max(0, +g.tau || 0) * 1e-6,
      tauD = Math.max(0, +g.tauD || 0) * 1e-6,
      LF = new Float64Array(N),
      pps = {};
    for (const li of keys) {
      const so = sols[li];
      pps[li] = Float64Array.from(so.ids, (id, i) => ppd(so.info[i]));
    }
    const lines = [];
    for (const li of keys)
      sols[li].ids.forEach((id, i) => {
        const e = sols[li].info[i];
        if (gamHas(e)) for (const l of gamAll(e)) lines.push({ li, i, E: l[0], I: l[1], ep: Math.max(0, eff(l[0])) });
      });
    const C = combos.map((c, ci) => {
      const E = c.E,
        fw = gamFwhm(E, g),
        hw = K * fw,
        lo = E - hw,
        hi = E + hw,
        sg = measSigma(E, g),
        ep = Math.max(0, eff(E)),
        cov = measPhi(hw / sg) - measPhi(-hw / sg),
        idx = sols[c.li].ids.indexOf(c.id);
      const near = lines
        .filter(l => !(l.li === c.li && l.i === idx && l.E === E) && l.ep > 0 && Math.abs(l.E - E) < hw + 6 * measSigma(l.E, g))
        .map(l => ({ li: l.li, i: l.i, w: l.I * l.ep * (measPhi((hi - l.E) / measSigma(l.E, g)) - measPhi((lo - l.E) / measSigma(l.E, g))) }))
        .filter(l => l.w > 0);
      const z = () => new Float64Array(N);
      return Object.assign({}, c, {
        idx,
        fw,
        lo,
        hi,
        ep,
        cov,
        near,
        bkgRate: bk.roi(lo, hi),
        t: z(),
        t1: z(),
        t2: z(),
        tot: z(),
        bkg: z(),
        bkgExp: z(),
        net: z(),
        sig: z(),
        pu: z(),
        thTot: z(),
        thNet: z(),
        sTrue: z(),
        iOther: z(),
        sigStat: z(),
        sigSys: z(),
        rate: z(),
        dmax: 0
      });
    });
    for (let k = 0; k < N; k++) {
      const Ds = {};
      let t1e = 0,
        t2e = 0,
        pulses = 0;
      keys.forEach((li, q) => {
        const r = gamDecays(sols[li], prodM[li], Tirr, edges[k], edges[k + 1]);
        Ds[li] = r.D;
        if (q === 0) {
          t1e = r.t1e;
          t2e = r.t2e;
        }
        const pq = pps[li];
        for (let i = 0; i < pq.length; i++) pulses += r.D[i] * pq[i];
      });
      const live = Math.max(t2e - t1e, 1e-12);
      const R = pulses / live + bk.total,
        mu = R * tau,
        puLoss = 1 - Math.exp(-2 * mu),
        fD = tauD > 0 ? 1 / (1 + R * tauD) : 1, // non-paralysable DAQ dead time: fraction of the pulses recorded
        lossF = (m.pile ? 1 - puLoss : 1) * fD, // every recorded count — signal, interference and background — shares the pile-up and dead-time losses
        lossFth = m.thPile ? lossF : 1;
      LF[k] = lossF;
      for (const c of C) {
        let io = 0;
        if (m.infer) for (const l of c.near) io += Ds[l.li][l.i] * l.w;
        const S0 = Ds[c.li][c.idx] * c.I * c.ep * c.cov,
          S = S0,
          Bk = c.bkgRate * live * lossF,
          mean = (S + io) * lossF + Bk;
        const Tn = m.noise ? measPois(mean, rnd) : mean,
          Bn = m.noise ? measPois(Bk, rnd) : Bk;
        c.t1[k] = t1e;
        c.t2[k] = t2e;
        c.t[k] = (t1e + t2e) / 2;
        c.tot[k] = Tn;
        c.bkg[k] = Bn;
        c.bkgExp[k] = Bk;
        c.net[k] = Tn - Bn;
        c.sigStat[k] = Math.sqrt(Math.max(Tn + Bn, 0));
        c.sigSys[k] = Math.sqrt(dI * dI + dE * dE) * Math.abs(Tn - Bn);
        c.sig[k] = Math.hypot(c.sigStat[k], c.sigSys[k]);
        c.pu[k] = puLoss * 100;
        c.rate[k] = R;
        c.thNet[k] = (S0 + io) * lossFth;
        c.thTot[k] = c.thNet[k] + (m.thBkg ? c.bkgRate * live * lossFth : 0);
        c.sTrue[k] = S;
        if (Ds[c.li][c.idx] > c.dmax) c.dmax = Ds[c.li][c.idx];
        c.iOther[k] = io;
      }
    }
    MEAS.res = {
      keys,
      sols,
      prods: prodM,
      edges,
      live: Float64Array.from(C[0].t1, (v, k) => C[0].t2[k] - v),
      puLoss: Float64Array.from(C[0].pu, v => v / 100),
      lossF: LF,
      C,
      N,
      NI,
      NC,
      dtI,
      dtC,
      kI,
      kC,
      Tirr,
      Tc,
      T,
      seed,
      capped,
      dI,
      dE,
      bk,
      g: Object.assign({}, g),
      m: Object.assign({}, m)
    };
    if (msg)
      msg.textContent =
        N +
        ' bins (' +
        NI +
        ' of ' +
        fmtHl(dtI) +
        ' during the irradiation' +
        (NC ? ', ' + NC + ' of ' + fmtHl(dtC) + ' during the cooling' : '') +
        ') · random seed ' +
        seed +
        (capped ? ' · limited to ' + MEAS_MAXN + ' points' : '') +
        '.';
    measSpecCompute();
    measDraw();
  } catch (e) {
    if (msg) msg.textContent = 'Simulation failed: ' + ((e && e.message) || e);
  }
  MEAS.busy = false;
  if (MEAS.again) {
    MEAS.again = false;
    measRun();
  }
}

// ---- plot ----
function measOpts() {
  const R = MEAS.res,
    m = measSet(),
    ax = measAx(),
    f = ax.f,
    H = MEAS.hidden,
    ser = [];
  R.C.forEach((c, ci) => {
    const nm = measNm(c) + ' ' + fmt(c.E, 6) + ' keV',
      col = scol(ci),
      x = Array.from(c.t, v => v / f),
      rt = m.qty === 'rate',
      q = (a, k) => (rt ? a[k] / Math.max(c.t2[k] - c.t1[k], 1e-12) : a[k]),
      Y = a => Array.from(a, (v, k) => q(a, k)),
      grp = nm;
    if (m.showTot)
      ser.push({
        name: nm + ' · total in ROI',
        type: 'pts',
        shape: 'c',
        x,
        y: Y(c.tot),
        err: m.err ? Array.from(c.tot, (v, k) => Math.hypot(Math.sqrt(v), c.sigSys[k]) / (rt ? Math.max(c.t2[k] - c.t1[k], 1e-12) : 1)) : null,
        color: col,
        ci,
        kind: 'tot',
        grp,
        xlo: m.xbar ? Array.from(c.t1, v => v / f) : null,
        xhi: m.xbar ? Array.from(c.t2, v => v / f) : null
      });
    if (m.showNet)
      ser.push({
        name: nm + ' · net (background subtracted)',
        type: 'pts',
        shape: 's',
        x,
        y: Y(c.net),
        err: m.err ? Array.from(c.sig, (v, k) => v / (rt ? Math.max(c.t2[k] - c.t1[k], 1e-12) : 1)) : null,
        color: col,
        ci,
        kind: 'net',
        r: 2.2,
        grp,
        xlo: m.xbar ? Array.from(c.t1, v => v / f) : null,
        xhi: m.xbar ? Array.from(c.t2, v => v / f) : null
      });
    if (m.th && m.showTot)
      ser.push({
        name: nm + ' · expected total' + (m.thBkg ? '' : ' (no background)'),
        type: 'line',
        x,
        y: Y(c.thTot),
        color: col,
        dash: '6 4',
        w: 1.4,
        op: 0.9,
        ci,
        kind: 'ttot',
        noread: true,
        grp
      });
    if (m.th && m.showNet)
      ser.push({ name: nm + ' · expected net', type: 'line', x, y: Y(c.thNet), color: col, w: 1.4, op: 0.9, ci, kind: 'tnet', noread: true, grp });
  });
  ser.forEach(s => {
    s.hidden = !!H[s.grp];
  });
  const z = MEAS.z,
    uT = { s: 'seconds', min: 'minutes', h: 'hours', d: 'days', y: 'years' }[ax.k] || ax.u;
  return {
    series: ser,
    xlog: false,
    ylog: m.ylog === 'log',
    xmin: z.xmin || '',
    xmax: z.xmax || '',
    ymin: z.ymin || '',
    ymax: z.ymax || '',
    xlabel: 'Time since the start of the irradiation (' + uT + ')',
    ylabel:
      m.qty === 'rate'
        ? 'Count rate (counts/s)'
        : 'Counts per time bin (' + fmtHl(R.dtI) + (R.NC ? ' during the irradiation, ' + fmtHl(R.dtC) + ' during the cooling' : '') + ')',
    title: 'Simulated HPGe measurement · irradiation ' + fmtHl(R.Tirr) + (R.Tc > 0 ? ' + cooling ' + fmtHl(R.Tc) : ''),
    bands: MEAS.slice ? [{ x0: MEAS.slice.t1 / f, x1: MEAS.slice.t2 / f, label: 'γ spectrum' }] : [],
    vlines: R.Tc > 0 ? [{ x: R.Tirr / f, label: 'end of irradiation' }] : []
  };
}
function measDraw() {
  const box = $('#measOut'),
    R = MEAS.res;
  if (!box || !R) return;
  const m = measSet(),
    o = measOpts(),
    r = renderChart(o);
  MEAS.o = o;
  MEAS.geom = r.geom;
  const names = Array.from(new Set(o.series.map(s => s.grp))).map(n => [n, (o.series.find(s => s.grp === n) || {}).color]);
  MEAS.names = names.map(x => x[0]);
  const leg = gamRen(gamLegend(names, MEAS.hidden, 'meas'), 'm');
  const rows = R.C.map(c => {
    let sn = 0,
      mx = 0;
    for (let k = 0; k < R.N; k++) {
      sn += c.net[k];
      if (c.pu[k] > mx) mx = c.pu[k];
    }
    return (
      '<tr><td>' +
      esc(measNm(c)) +
      '</td><td>' +
      (c.e && c.e.T ? fmtHl(c.e.T) : '—') +
      '</td><td>' +
      fmt(c.E, 6) +
      '</td><td>' +
      fmt(c.I * 100, 4) +
      '</td><td>' +
      fmt(c.ep, 3) +
      '</td><td>' +
      fmt(c.fw, 3) +
      '</td><td>' +
      fmt(c.lo, 6) +
      ' – ' +
      fmt(c.hi, 6) +
      '</td><td>' +
      fmt(c.cov * 100, 4) +
      '</td><td>' +
      fmt(c.bkgRate, 3) +
      '</td><td>' +
      fmt(sn, 5) +
      '</td><td>' +
      fmt(mx, 3) +
      '</td></tr>'
    );
  }).join('');
  box.innerHTML =
    '<div class="btns"><span class="lbl">Plot</span><select data-mk="qty"><option value="cnt"' +
    (m.qty === 'rate' ? '' : ' selected') +
    '>counts per bin</option><option value="rate"' +
    (m.qty === 'rate' ? ' selected' : '') +
    '>count rates (counts/s)</option></select><span class="lbl">Show</span>' +
    [
      ['showTot', 'total counts in the ROI'],
      ['showNet', 'net counts in the peak'],
      ['err', 'vertical error bars (±1σ)'],
      ['xbar', 'horizontal bars (time-bin width)'],
      ['th', 'theoretical lines'],
      ['thPile', 'pile-up in the theory'],
      ['thBkg', 'background in the expected total']
    ]
      .map(x => '<label class="chk"><input type="checkbox" data-mk="' + x[0] + '"' + (m[x[0]] ? ' checked' : '') + '><span>' + x[1] + '</span></label>')
      .join('') +
    '</div>' +
    gamRen(
      pcTop({
        y: ['data-gyl="meas"', m.ylog],
        reset: ['data-gzreset="meas"', gamZoomed(MEAS.z)],
        extra:
          '<button class="btn sm' +
          (MEAS.pick ? ' on' : '') +
          '" data-mpick title="Then drag over a time range of the plot: the γ-ray spectrum in the detector for that range appears below">' +
          (MEAS.pick ? 'Drag a time range on the plot…' : 'Select a time range → γ spectrum') +
          '</button>',
        hint: 'Drag a box to zoom (a thin band zooms one axis) · Shift+drag pans · double-click resets · hover a point for its numbers'
      }),
      'm'
    ) +
    '<div class="chart" id="measChart" style="position:relative">' +
    r.svg +
    '</div>' +
    leg +
    pcSave({ attr: 'data-msave', formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'rootbin', 'root'], msg: '<span class="hint" id="measSaveMsg"></span>' }) +
    '<div class="tblwrap"><table class="tbl"><thead><tr><th>Isotope</th><th>T½</th><th>Eγ (keV)</th><th>Iγ (%)</th><th>ε peak</th><th>FWHM (keV)</th><th>ROI = Eγ ± ' +
    fmt(+m.roiK || 1.5, 3) +
    '·FWHM (keV)</th><th>peak in ROI (%)</th><th>background in ROI (cps)</th><th>Σ net counts</th><th>max pile-up (%)</th></tr></thead><tbody>' +
    rows +
    '</tbody></table></div>' +
    measBinNote(R) +
    measSatNote(R) +
    measZeroNote(R) +
    '<div id="measSpec"></div>';
  measAttach();
  measSpecDraw();
}
function measZeroNote(R) {
  // explains why the expected signal of a line is zero in every bin
  const out = [];
  for (const c of R.C) {
    let smax = 0;
    for (let k = 0; k < R.N; k++) if (c.sTrue[k] > smax) smax = c.sTrue[k];
    if (smax > 0) continue;
    const why = [];
    if (c.idx < 0) why.push('the nuclide is not in the solved decay chain (index −1)');
    else if (!(c.dmax > 0)) why.push('no decays of this nuclide in any time bin (decays per bin = ' + c.dmax + ')');
    if (!(c.I > 0)) why.push('γ intensity I = ' + c.I);
    if (!(c.ep > 0)) why.push('peak efficiency ε = ' + c.ep + ' at ' + fmt(c.E, 6) + ' keV (Ge/DAQ configuration)');
    if (!(c.cov > 0)) why.push('ROI covers ' + c.cov + ' of the peak (resolution σ = ' + fmt(measSigma(c.E, R.g), 4) + ' keV)');
    out.push('<li>' + esc(measNm(c)) + ' ' + fmt(c.E, 6) + ' keV: ' + (why.length ? why.join('; ') : 'unknown reason') + '</li>');
  }
  // zero or non-finite counts of the run (any line)
  let st = 0,
    sb = 0,
    bad = 0,
    lmin = Infinity,
    lmax = 0;
  for (const c of R.C)
    for (let k = 0; k < R.N; k++) {
      st += c.tot[k];
      sb += c.bkg[k];
      if (!isFinite(c.tot[k]) || !isFinite(c.rate[k]) || !isFinite(c.bkgExp[k])) bad++;
    }
  for (let k = 0; k < R.N; k++) {
    lmin = Math.min(lmin, R.live[k]);
    lmax = Math.max(lmax, R.live[k]);
  }
  const zero = !(st > 0) || bad > 0 || !(lmin > 1e-6);
  const diag = zero
    ? '<div class="warn"><b>Diagnostic (all counts are zero or not finite):</b> Σ total counts = ' +
      st +
      ', Σ background counts = ' +
      sb +
      ', non-finite values = ' +
      bad +
      ', live time per bin: min ' +
      fmt(lmin, 4) +
      ' s, max ' +
      fmt(lmax, 4) +
      ' s, background rate in the ROI of the first line = ' +
      fmt(R.C[0].bkgRate, 4) +
      ' cps, background total = ' +
      fmt(R.bk.total, 4) +
      ' cps, efficiency ε(first line) = ' +
      fmt(R.C[0].ep, 4) +
      ', FWHM = ' +
      fmt(R.C[0].fw, 4) +
      ' keV.</div>'
    : '';
  return diag + (out.length ? '<div class="warn"><b>Expected signal is exactly zero for:</b><ul>' + out.join('') + '</ul></div>' : '');
}
function measBinNote(R) {
  // the bins of the irradiation and of the cooling may have different widths: counts per bin then jump at the end of the irradiation although the count rate is continuous
  if (!R.NC || measSet().qty === 'rate' || Math.abs(R.dtC / R.dtI - 1) < 0.01) return '';
  const r = R.dtC / R.dtI;
  return (
    '<p class="warn">The time bins are ' +
    esc(fmtHl(R.dtI)) +
    ' wide during the irradiation and ' +
    esc(fmtHl(R.dtC)) +
    ' during the cooling, so the <b>counts per bin change by a factor ' +
    fmt(r, 3) +
    ' at the end of the irradiation</b> although the count rate is continuous. Use the same bin width in both phases (points per time unit) or plot the count rate: <button class="btn sm" data-mqtyrate>Plot counts/s</button></p>'
  );
}
function measSatNote(R) {
  // warning when the detector is saturated by the pile-up (the paralysable dead-time factor exp(−2Rτ) tends to zero)
  let rmax = 0,
    pmax = 0,
    nsat = 0;
  const c0 = R.C[0],
    k0 = R.NI < R.N ? R.NI : 0; // only the cooling bins count for the warning (the beam is on during the irradiation); without a cooling phase all bins are listed but no warning is given
  for (let k = k0; k < R.N; k++) {
    if (c0.rate[k] > rmax) rmax = c0.rate[k];
    if (c0.pu[k] > pmax) pmax = c0.pu[k];
    if (R.NI < R.N && c0.pu[k] > 99) nsat++;
  }
  const b = measBeam(R.m),
    tau = Math.max(0, +R.g.tau || 0),
    tauD = Math.max(0, +R.g.tauD || 0);
  return (
    '<p class="' +
    (nsat ? 'warn' : 'hint') +
    '">Beam ' +
    fmt(b.I, 4) +
    ' ' +
    esc(b.IU) +
    (b.q ? ', q = ' + b.q : '') +
    ' (' +
    fmt(b.rate, 4) +
    ' primaries/s) · maximum total pulse rate in the detector ' +
    (k0 ? 'during the cooling ' : '') +
    fmt(rmax, 4) +
    ' counts/s with τ = ' +
    fmt(tau, 3) +
    ' µs → pile-up loss up to ' +
    fmt(pmax, 4) +
    ' %' +
    (tauD > 0
      ? ' · DAQ dead time ' +
        fmt(tauD, 3) +
        ' µs (non-paralysable) → up to ' +
        fmt((100 * rmax * tauD * 1e-6) / (1 + rmax * tauD * 1e-6), 4) +
        ' % of the pulses lost'
      : '') +
    (nsat
      ? ' · <b>' +
        nsat +
        ' of ' +
        (R.N - R.NI) +
        ' cooling bins are saturated (loss above 99 %)</b>: the peak counts there vanish because the detector is blind. Reduce the beam current or the efficiency (distance) or the pulse width, or start counting later.'
      : '') +
    '</p>'
  );
}
function measCtx() {
  const R = MEAS.res;
  if (!R) return null;
  const o = measOpts(),
    r = renderChart(o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  return { o, geom: r.geom, svg: div.querySelector('svg'), base: 'cpt_measurement', key: 'meas' };
}
function measTipEl() {
  let t = document.getElementById('measTipWin');
  if (!t) {
    t = document.createElement('div');
    t.id = 'measTipWin';
    t.className = 'mtip';
    t.hidden = true;
    document.body.appendChild(t);
  }
  return t;
}
function measAttach() {
  const w = $('#measWin'),
    sv = $('#measChart svg', w),
    box = $('#measChart', w),
    tip = measTipEl(),
    gg = MEAS.geom,
    R = MEAS.res;
  if (!sv || !gg) return;
  gamZoom(sv, gg, 'z', MEAS, measDraw);
  const pt = ev => {
    const q = sv.getBoundingClientRect();
    return [((ev.clientX - q.left) / q.width) * gg.W, ((ev.clientY - q.top) / q.height) * gg.H];
  };
  const cx = v => gg.m.l + ((v - gg.x0) / (gg.x1 - gg.x0)) * gg.iw,
    cy = v =>
      gg.m.t + gg.ih - ((gg.ylog ? Math.log10(v) - Math.log10(gg.y0) : v - gg.y0) / (gg.ylog ? Math.log10(gg.y1) - Math.log10(gg.y0) : gg.y1 - gg.y0)) * gg.ih;
  let hl = null;
  const clear = () => {
    tip.hidden = true;
    if (hl) {
      hl.remove();
      hl = null;
    }
  };
  sv.addEventListener('pointermove', ev => {
    if (ev.buttons) {
      clear();
      return;
    }
    const [px, py] = pt(ev);
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) {
      clear();
      return;
    }
    let best = null,
      bd = 15 * 15;
    for (const s of MEAS.o.series) {
      if (s.hidden || s.type !== 'pts') continue;
      for (let i = 0; i < s.x.length; i++) {
        if (!isFinite(s.y[i]) || (gg.ylog && s.y[i] <= 0)) continue;
        const dx = cx(s.x[i]) - px,
          dy = cy(s.y[i]) - py,
          dd = dx * dx + dy * dy;
        if (dd < bd) {
          bd = dd;
          best = { s, i };
        }
      }
    }
    if (!best) {
      clear();
      return;
    }
    const { s, i } = best,
      c = R.C[s.ci],
      ax = measAx(),
      f = ax.f,
      uT = ax.u,
      net = c.net[i],
      sg = c.sig[i],
      lv = Math.max(c.t2[i] - c.t1[i], 1e-12);
    tip.innerHTML =
      '<b>' +
      esc(measNm(c)) +
      ' · ' +
      fmt(c.E, 6) +
      ' keV</b> (Iγ ' +
      fmt(c.I * 100, 4) +
      ' %)<br>' +
      (c.t2[i] <= R.Tirr * (1 + 1e-9) ? 'irradiation' : 'cooling') +
      ' · t = ' +
      fmt(c.t[i] / f, 5) +
      ' ' +
      uT +
      ' · bin ' +
      fmt(c.t1[i] / f, 5) +
      ' – ' +
      fmt(c.t2[i] / f, 5) +
      ' ' +
      uT +
      '<br>ROI ' +
      fmt(c.lo, 6) +
      ' – ' +
      fmt(c.hi, 6) +
      ' keV<br>Total counts in ROI: <b>' +
      fmt(c.tot[i], 6) +
      '</b><br>Background counts: <b>' +
      fmt(c.bkg[i], 6) +
      '</b> (expected ' +
      fmt(c.bkgExp[i], 4) +
      ')<br>Net counts: <b>' +
      fmt(net, 6) +
      ' ± ' +
      fmt(sg, 4) +
      '</b> (' +
      (Math.abs(net) > 0 ? fmt((sg / Math.abs(net)) * 100, 3) + ' %' : '—') +
      ')' +
      (c.sigSys[i] > 0 ? '<br>&nbsp;&nbsp;statistical ± ' + fmt(c.sigStat[i], 4) + ', beam current + efficiency ± ' + fmt(c.sigSys[i], 4) : '') +
      '<br>Rates (live ' +
      fmt(c.t2[i] - c.t1[i], 5) +
      ' s): total <b>' +
      fmt(c.tot[i] / lv, 5) +
      '</b>, background <b>' +
      fmt(c.bkg[i] / lv, 5) +
      '</b>, net (peak) <b>' +
      fmt(net / lv, 5) +
      ' ± ' +
      fmt(sg / lv, 4) +
      '</b> counts/s<br>Pile-up level: <b>' +
      fmt(c.pu[i], 3) +
      ' %</b> (total detector rate ' +
      fmt(c.rate[i], 5) +
      ' counts/s from all γ lines of ' +
      (R.keys.length > 1 ? 'layers ' + R.keys.map(k => k + 1).join(' + ') : 'layer ' + (R.keys[0] + 1)) +
      ' and the background)';
    tip.hidden = false;
    const tw = tip.offsetWidth,
      th = tip.offsetHeight,
      vw = document.documentElement.clientWidth,
      vh = document.documentElement.clientHeight,
      dx = 14,
      dy = 14; // lower-left corner of the window at the cursor plus a small shift (the cursor circle stays visible)
    let L = ev.clientX + dx,
      T = ev.clientY - dy - th;
    if (L + tw > vw - 4) L = ev.clientX - dx - tw;
    if (T < 4) T = ev.clientY + dy + 6;
    tip.style.left = Math.max(4, Math.min(L, vw - tw - 4)) + 'px';
    tip.style.top = Math.max(4, Math.min(T, vh - th - 4)) + 'px';
    if (hl) hl.remove();
    hl = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    hl.setAttribute('cx', cx(s.x[i]));
    hl.setAttribute('cy', cy(s.y[i]));
    hl.setAttribute('r', 6);
    hl.setAttribute('fill', 'none');
    hl.setAttribute('stroke', 'var(--ink)');
    hl.setAttribute('stroke-width', 1.6);
    sv.appendChild(hl);
  });
  sv.addEventListener('pointerleave', clear);
  // select a time range for the γ spectrum
  let pk = null;
  sv.addEventListener('pointerdown', ev => {
    if (!MEAS.pick || ev.button !== 0) return;
    const [px, py] = pt(ev);
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) return;
    const rc = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rc.setAttribute('fill', 'var(--accent)');
    rc.setAttribute('fill-opacity', '.18');
    rc.setAttribute('stroke', 'var(--accent)');
    rc.setAttribute('y', gg.m.t);
    rc.setAttribute('height', gg.ih);
    rc.setAttribute('x', px);
    rc.setAttribute('width', 0);
    sv.appendChild(rc);
    pk = { x0: px, rc };
    sv.setPointerCapture(ev.pointerId);
    clear();
  });
  sv.addEventListener('pointermove', ev => {
    if (!pk) return;
    const [px] = pt(ev),
      a = Math.min(Math.max(px, gg.m.l), gg.m.l + gg.iw);
    pk.rc.setAttribute('x', Math.min(a, pk.x0));
    pk.rc.setAttribute('width', Math.abs(a - pk.x0));
  });
  sv.addEventListener('pointerup', ev => {
    if (!pk) return;
    const [px] = pt(ev),
      a = Math.min(Math.max(px, gg.m.l), gg.m.l + gg.iw),
      x0 = pk.x0;
    pk.rc.remove();
    pk = null;
    if (Math.abs(a - x0) < 6) return;
    const X0 = pxData(gg, Math.min(a, x0), gg.m.t)[0],
      X1 = pxData(gg, Math.max(a, x0), gg.m.t)[0];
    MEAS.pick = false;
    measSetSlice(X0 * measAx().f, X1 * measAx().f);
  });
  // shift + drag pans
  let pan = null;
  sv.addEventListener('pointerdown', ev => {
    if (!ev.shiftKey || ev.button !== 0) return;
    const [px, py] = pt(ev);
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) return;
    pan = { px, py, r: [gg.x0, gg.x1, gg.y0, gg.y1] };
    sv.setPointerCapture(ev.pointerId);
    clear();
  });
  sv.addEventListener('pointermove', ev => {
    if (!pan) return;
    const [px, py] = pt(ev),
      fx = (px - pan.px) / gg.iw,
      fy = (py - pan.py) / gg.ih,
      [a, b, c, d] = pan.r,
      z = MEAS.z;
    z.xmin = fmtAx(a - fx * (b - a));
    z.xmax = fmtAx(b - fx * (b - a));
    if (gg.ylog) {
      const l0 = Math.log10(c),
        l1 = Math.log10(d);
      z.ymin = fmtAx(Math.pow(10, l0 + fy * (l1 - l0)));
      z.ymax = fmtAx(Math.pow(10, l1 + fy * (l1 - l0)));
    } else {
      z.ymin = fmtAx(c + fy * (d - c));
      z.ymax = fmtAx(d + fy * (d - c));
    }
    if (!pan.raf)
      pan.raf = requestAnimationFrame(() => {
        const p = pan;
        if (p) p.raf = 0;
        measDraw();
      });
  });
  sv.addEventListener('pointerup', () => {
    pan = null;
  });
}

// ---- window ----
function measBkgOpts() {
  // options of the background-preview plot (also used for the exports)
  const m = measSet(),
    g = gamSet(),
    bk = measBkg(m, g);
  if (!(bk.total > 0)) return null;
  const z = MEAS.bz || (MEAS.bz = {}),
    xl = m.bxl === 'log',
    H = MEAS.bh || (MEAS.bh = {}),
    x = [],
    y = [];
  for (let k = 0; k <= 700; k++) {
    const E = (3500 * k) / 700;
    x.push(E);
    y.push(bk.spec(E));
  }
  return {
    series: [{ name: 'background', type: 'line', x, y, color: 'var(--accent)', w: 1.6, hidden: !!H.background }],
    xlog: xl,
    ylog: (m.byl || 'log') === 'log',
    xmin: z.xmin || (xl ? '1' : '0'),
    xmax: z.xmax || '3500',
    ymin: z.ymin || '',
    ymax: z.ymax || '',
    xlabel: 'Energy (keV)',
    ylabel: 'Background rate (counts / s / keV)',
    title: 'Background model · total ' + fmt(bk.total, 3) + ' cps',
    H: 300
  };
}
function measBkgDraw() {
  const box = $('#measBkgChart');
  if (!box) return;
  const o = measBkgOpts();
  if (!o) {
    box.innerHTML = '<p class="hint">No background.</p>';
    return;
  }
  const r = renderChart(o);
  MEAS.bo = o;
  const top = gamRen(
      pcTop({
        y: ['data-gyl="bk"', measSet().byl || 'log'],
        x: ['data-bgxl="bk"', measSet().bxl || 'lin'],
        reset: ['data-gzreset="bk"', gamZoomed(MEAS.bz)]
      }),
      'b'
    ),
    leg = gamRen(gamLegend([['background', 'var(--accent)']], MEAS.bh || {}, 'bk'), 'b');
  box.innerHTML =
    top +
    '<div class="chart" id="measBkgPlot">' +
    r.svg +
    '</div>' +
    leg +
    pcSave({ attr: 'data-mbsave', formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'rootbin', 'root'], msg: '<span class="hint" id="measBkgMsg"></span>' });
  gamZoom($('#measBkgPlot svg', box), r.geom, 'bz', MEAS, measBkgDraw);
  specTip($('#measBkgPlot svg', box), r.geom, () => MEAS.bo);
}
async function measBkgSave(kind) {
  const m = $('#measBkgMsg'),
    o = measBkgOpts();
  if (!o) return;
  if (m) m.textContent = 'Preparing…';
  const r = renderChart(o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  const ctx = { o, geom: r.geom, svg: div.querySelector('svg'), base: 'cpt_background_spectrum', key: 'meas-bkg' },
    old = _chart;
  _chart = ctx;
  try {
    const px = expPx(),
      base = ctx.base,
      sel = '#measBkgMsg';
    if (kind === 'png') await saveFile(base + '.png', await canvasBlob(await svgToCanvas(px), 'image/png'), 'PNG', sel);
    else if (kind === 'pdf') await saveFile(base + '.pdf', await makePdf(await svgToCanvas(Math.max(px, 2400)), ctx.geom), 'PDF', sel);
    else if (kind === 'svg') await saveFile(base + '.svg', exportSvgString(), 'SVG', sel);
    else if (kind === 'csv') await saveFile(base + '.csv', csvOf(), 'CSV', sel);
    else if (kind === 'txt') await saveFile(base + '.txt', asciiOf(), 'ASCII', sel);
    else if (kind === 'root') await saveFile(base + '.C.txt', rootMacroOf(), 'ROOT macro', sel);
    else if (kind === 'rootbin') {
      const hs = rootHistograms(),
        rn = base + '.root',
        rb = RootIO.file(rn, hs),
        inClaude = !!(window.claude && claude.use && (await claude.use('downloads')));
      if (inClaude) await saveFile(base + '_histograms.zip', zipStore(rn, rb), 'ROOT file', sel);
      else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', sel);
      if (m) m.textContent = 'Saved ' + (inClaude ? base + '_histograms.zip (containing ' + rn + ')' : rn) + ' with ' + hs.length + ' histograms.';
    }
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  } finally {
    _chart = old;
  }
}
function measInfoHtml() {
  const m = measSet(),
    g = gamSet(),
    ab = measAB(g),
    K = +m.roiK || 1.5,
    E = 661.657,
    fw = gamFwhm(E, g);
  return (
    '<b>Resolution:</b> FWHM(E) = 2.355·σ(E), σ²(E) = a + b·E with a = ' +
    fmt(ab.a, 4) +
    ' keV², b = ' +
    fmt(ab.b, 4) +
    ' keV (from the Ge/DAQ configuration). <b>ROI</b> of a γ line = E<sub>γ</sub> ± ' +
    fmt(K, 3) +
    '·FWHM(E<sub>γ</sub>); e.g. at 661.7 keV: FWHM = ' +
    fmt(fw, 3) +
    ' keV → ROI = ' +
    fmt(E - K * fw, 6) +
    ' – ' +
    fmt(E + K * fw, 6) +
    ' keV (' +
    fmt(2 * K * fw, 3) +
    ' keV wide). A worse resolution widens the ROI, so it collects more background and the signal-to-background ratio drops.'
  );
}
function measBuild() {
  const w = $('#measWin'),
    m = measSet(),
    g = gamSet(),
    keys = MEAS.keys || [],
    R0 = RES;
  const nI = (k, ex) => '<input type="number" data-mk="' + k + '" step="any" value="' + m[k] + '" ' + (ex || '') + '>',
    sel = (k, o) =>
      '<select data-mk="' +
      k +
      '">' +
      o.map(x => '<option value="' + x[0] + '"' + (String(m[k]) === String(x[0]) ? ' selected' : '') + '>' + x[1] + '</option>').join('') +
      '</select>',
    fld = (l, h, c) => '<label class="fld' + (c ? ' ' + c : '') + '"><span>' + l + '</span>' + h + '</label>',
    chk = (k, t) => '<label class="chk"><input type="checkbox" data-mk="' + k + '"' + (m[k] ? ' checked' : '') + '><span>' + t + '</span></label>';
  w.innerHTML =
    '<div class="mbox" style="width:min(1200px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:10px"><h3>Simulate a decay measurement (HPGe)</h3><button class="btn sm" data-mclose>Done</button></div>' +
    '<p class="hint">Counts of selected γ lines in an HPGe detector, in time bins during the irradiation and the cooling. Source: the production rates and decay chains of every layer of the activation calculation (beam current of that panel), solved with the Bateman equations at every bin edge. The detector counts the chosen layer or all layers together: the pile-up and the γ spectrum depend on the total counting rate of all γ lines of those layers.</p>' +
    '<div class="mgrid">' +
    fld(
      'Beam current (empty = as in the activation panel)',
      '<input type="number" data-mk="I" step="any" min="0" value="' + (m.I === '' ? '' : m.I) + '" placeholder="' + esc(String(irrSet().I)) + '">'
    ) +
    fld(
      'Current unit',
      '<select data-mk="IU">' +
        [
          ['', 'as in the activation panel (' + (irrSet().IU || 'pps') + ')'],
          ['pps', 'particles / s'],
          ['pA', 'pA'],
          ['nA', 'nA'],
          ['uA', 'µA'],
          ['mA', 'mA']
        ]
          .map(x => '<option value="' + x[0] + '"' + (String(m.IU) === x[0] ? ' selected' : '') + '>' + x[1] + '</option>')
          .join('') +
        '</select>'
    ) +
    fld(
      'Charge state q of the ions (e), empty = activation panel',
      '<input type="number" data-mk="q" step="1" min="1" value="' +
        (m.q === '' ? '' : m.q) +
        '" placeholder="' +
        esc(String(irrSet().q || 'fully stripped')) +
        '">'
    ) +
    fld('Irradiation time', nI('tirr', 'min="0"')) +
    fld('Irradiation time unit', sel('tU', actUnits)) +
    fld('Points per irradiation time unit', nI('ppu', 'min="0.01"')) +
    fld('Cooling time', nI('tcool', 'min="0"')) +
    fld('Cooling time unit', sel('tcU', actUnits)) +
    fld('Points per cooling time unit', nI('ppuC', 'min="0.01"')) +
    fld('Time unit of the plot axis', sel('xU', [['auto', 'automatic (longer phase)']].concat(actUnits))) +
    fld('Beam-current uncertainty (%, relative)', nI('dI', 'min="0" step="any"')) +
    fld('Ge efficiency uncertainty (%, relative)', nI('dE', 'min="0" step="any"')) +
    '</div>' +
    '<div class="btns"><button class="btn pri" data-msel>Select isotopes and γ lines</button><span class="hint" id="measSelSum"></span><button class="btn pri" data-gecfg>Configure the Ge detector and DAQ</button><button class="btn" data-mreset title="Deletes the saved settings of this window (time, beam, background, statistics) and restores the defaults">Reset this window to the defaults</button><span class="hint" data-gcsum>' +
    esc(gamCfgSummary(g)) +
    '</span></div>' +
    '<h4>Background</h4><div class="mgrid">' +
    fld(
      'Background model',
      sel('bmode', [
        ['none', 'A · none'],
        ['exp', 'B · smooth exponential'],
        ['cont', 'C · continuum + γ peaks (Pb X-rays, 40K, 208Tl, 511 keV, …)']
      ])
    ) +
    fld('Total background rate, 0–' + BKG_E1 + ' keV (counts/s)', nI('brate', 'min="0"')) +
    fld('Low-energy efficiency cut-off E<sub>c</sub> (keV)', nI('bEc', 'min="0.1" step="1"')) +
    fld('Mode B: 1/λ (keV)', nI('bL', 'min="1"')) +
    '</div>' +
    '<details><summary>Mode C: peak list (energy in keV, area in counts/s at ' +
    BKG_R0 +
    ' cps total, label) and the model</summary><textarea data-mk="bpk" rows="8" spellcheck="false" style="font-family:var(--font-mono);width:100%">' +
    esc(m.bpk) +
    '</textarea>' +
    '<p class="hint">B(E) = ε(E)·[A₁e<sup>−E/L₁</sup> + A₂e<sup>−E/L₂</sup>] + Σ ε(E<sub>i</sub>)·I<sub>i</sub>/(σ<sub>i</sub>√2π)·exp(−(E−E<sub>i</sub>)²/2σ<sub>i</sub>²), defined from 0 to ' +
    BKG_E1 +
    ' keV, with L₁ = ' +
    BKG_SOFT +
    ' keV (soft, ambient Compton), L₂ = ' +
    BKG_HARD +
    ' keV (hard, cosmic rays / bremsstrahlung), ' +
    Math.round(BKG_FSOFT * 100) +
    ' % of the continuum in the soft term. The empirical shape is multiplied by the relative detection efficiency ε(E) = exp(−(E<sub>c</sub>/E)³)·(1−exp(−E/' +
    BKG_ETH +
    ' keV)): absorption in the end cap, dead layer and contacts (photoelectric, ∝ E<sup>−3</sup>; 50 % transmission at ≈ 1.1·E<sub>c</sub>) and a soft electronic threshold, so that B(E) → 0 for E → 0. The continuum area is the total rate (over 0–' +
    BKG_E1 +
    ' keV) minus the peaks; all areas scale with the total rate. Peaks above ' +
    BKG_E1 +
    ' keV are added on top of the total.</p></details>' +
    '<details id="measBkgPrev"><summary>Preview of the background spectrum</summary><div id="measBkgChart"></div></details>' +
    '<h4>Counting statistics</h4><div class="btns">' +
    chk('noise', 'Poisson fluctuations') +
    chk('pile', 'apply pile-up losses to the counts') +
    chk('infer', 'count neighbouring γ lines of other nuclides in the ROI') +
    '</div>' +
    '<details><summary>Advanced: ROI width and random seed</summary><div class="btns" style="margin-top:8px">' +
    '<label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>ROI half-width (× FWHM)</span><input type="number" data-mk="roiK" step="0.1" min="0.1" value="' +
    m.roiK +
    '" style="width:100px"></label>' +
    '<label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>Random seed</span><input type="number" data-mk="seed" step="1" value="' +
    m.seed +
    '" style="width:120px"></label>' +
    chk('newSeed', 'new seed at every “Simulate Experiment”') +
    '</details>' +
    '<p class="hint" id="measInfo">' +
    measInfoHtml() +
    '</p>' +
    '<div class="btns"><button class="btn pri" id="measRun" data-mrun style="font-size:15px;padding:9px 20px">Simulate Experiment</button><span class="hint" id="measMsg"></span></div>' +
    '<div id="measOut"></div></div>';
  measSelSummary();
  measBkgDraw();
  if (MEAS.res) measDraw();
}
function measLayerChanged() {
  // the layer chosen for the measurement changed
  const m = measSet(),
    s = irrSet().meas;
  if (!measSelCombos(m).length) {
    const eff = gamEffFn(gamSet()),
      sc = measSelScope();
    MEAS.list
      .filter(n => sc === null || n.li === sc)
      .sort((a, b) => b.pk * (b.e.gtot || 0) - a.pk * (a.e.gtot || 0))
      .slice(0, 3)
      .forEach(n => {
        const l = n.lines.find(l => l[0] >= 50 && eff(l[0]) > 0) || n.lines[0];
        m.sel[n.li + '|' + n.id + '|' + l[0]] = true;
      });
    s.sel = m.sel;
  }
  saveLocal();
  const lm = $('#measWin [data-mk=layF]');
  if (lm) lm.value = String(m.layF);
  if (!$('#measSelWin').hidden) measSelRender();
  measSelSummary();
  if (MEAS.res) measLater(50);
}
function measSelSummary() {
  const m = measSet(),
    n = measSelCombos(m).length,
    e = $('#measSelSum');
  if (e)
    e.textContent =
      n +
      ' γ line' +
      (n === 1 ? '' : 's') +
      ' selected (' +
      (measSelScope() === null ? 'all layers' : measLayName(measSelScope()) + ' only') +
      ')' +
      (MEAS.list.length ? '' : ' — no nuclide with γ data (load ENDF-6 decay data in the activation panel)');
}
async function measPrep() {
  // solve the chains of every layer once to list the nuclides with γ lines
  const m = measSet(),
    fI = ACT_TU[m.tU] || 3600,
    fC = ACT_TU[m.tcU] || fI,
    g = gamSet(),
    eff = gamEffFn(g),
    list = [];
  for (const li of MEAS.keys || []) {
    const d = ACTP.get(li);
    if (!d || !d.Q) continue;
    const sol = await dkSolve(Object.assign({}, d.Q, { Tirr: Math.max(1e-3, +m.tirr * fI), Tc: Math.max(0, +m.tcool * fC), extra: [] }));
    list.push(
      ...sol.ids
        .map((id, i) => {
          const e = sol.info[i];
          let pk = 0;
          for (const v of sol.act[i]) if (v > pk) pk = v;
          return { li, i, id, e, pk, lines: e && e.gam ? e.gam.slice().sort((a, b) => b[1] - a[1]) : [] };
        })
        .filter(n => n.e && n.lines.length && n.pk > 0)
        .sort((a, b) => b.pk * (b.e.gtot || 0) - a.pk * (a.e.gtot || 0))
        .slice(0, 40)
    );
  }
  MEAS.list = list;
  const valid = new Set();
  list.forEach(n => n.lines.forEach(l => valid.add(n.li + '|' + n.id + '|' + l[0])));
  const s = irrSet().meas;
  for (const k of Object.keys(m.sel)) if (!valid.has(k)) delete m.sel[k];
  if (!measSelCombos(m).length)
    list
      .filter(n => {
        const sc = measSelScope();
        return sc === null || n.li === sc;
      })
      .sort((a, b) => b.pk * (b.e.gtot || 0) - a.pk * (a.e.gtot || 0))
      .slice(0, 3)
      .forEach(n => {
        const l = n.lines.find(l => l[0] >= 50 && eff(l[0]) > 0) || n.lines[0];
        m.sel[n.li + '|' + n.id + '|' + l[0]] = true;
      });
  s.sel = m.sel;
}
async function measOpen() {
  const keys = Array.from(ACTP.keys()).filter(k => typeof k === 'number');
  if (!keys.length) return;
  MEAS.keys = keys;
  measSet();
  $('#measWin').hidden = false;
  $('#measWin').innerHTML = '<div class="mbox"><h3>Simulate a decay measurement (HPGe)</h3><p class="hint">Solving the decay chains…</p></div>';
  try {
    await measPrep();
  } catch (e) {
    $('#measWin').innerHTML =
      '<div class="mbox"><h3>Simulate a decay measurement (HPGe)</h3><p class="warn">' +
      esc((e && e.message) || String(e)) +
      '</p><div class="btns"><button class="btn" data-mclose>Done</button></div></div>';
    return;
  }
  saveLocal();
  measBuild();
}
function measSelRender() {
  const w = $('#measSelWin'),
    m = measSet(),
    g = gamSet(),
    eff = gamEffFn(g),
    minI = Math.max(0, +m.minI || 0) / 100,
    sc = measSelScope(),
    keys = (MEAS.keys || []).filter(k => sc === null || k === sc);
  const rows = keys
    .map(li => {
      const body = MEAS.list
        .filter(n => n.li === li)
        .map(n => {
          const ls = n.lines.filter(l => l[1] >= minI).slice(0, 14);
          if (!ls.length) return '';
          return (
            '<tr class="nuc"><td colspan="4"><b>' +
            esc(dkNuc(n.id)) +
            '</b> · T½ ' +
            fmtHl(n.e.T) +
            ' · peak activity ' +
            fmt(n.pk, 3) +
            ' Bq · Σ Iγ = ' +
            fmt((n.e.gtot || 0) * 100, 4) +
            ' %</td></tr>' +
            ls
              .map(
                l =>
                  '<tr><td><input type="checkbox" data-msk="' +
                  n.li +
                  '|' +
                  n.id +
                  '|' +
                  l[0] +
                  '"' +
                  (m.sel[n.li + '|' + n.id + '|' + l[0]] ? ' checked' : '') +
                  '></td><td>' +
                  fmt(l[0], 6) +
                  '</td><td>' +
                  fmt(l[1] * 100, 4) +
                  '</td><td>' +
                  fmt(Math.max(0, eff(l[0])), 3) +
                  '</td></tr>'
              )
              .join('')
          );
        })
        .join('');
      return (
        '<tr class="laysec"><td colspan="4" style="background:var(--accent);color:#fff;font-weight:700">' +
        esc(measLayName(li)) +
        '</td></tr>' +
        (body || '<tr><td colspan="4" class="hint">No nuclide with discrete γ lines in this layer.</td></tr>')
      );
    })
    .join('');
  w.innerHTML =
    '<div class="mbox" style="width:min(720px,100%)"><div style="display:flex;justify-content:space-between;align-items:center;gap:10px"><h3>Select isotopes and γ lines</h3><button class="btn sm" data-mselclose>Done</button></div>' +
    '<p class="hint">The nuclides of the decay chains, classified by layer (strongest 40 per layer, by activity × γ yield), and their discrete γ lines from the ENDF-6 decay data. Tick one or several lines; each one gets its own ROI and its own curves. Choose all layers or only one.</p>' +
    '<div class="btns"><label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>Layers shown</span><select data-mlayf><option value="all"' +
    (sc === null ? ' selected' : '') +
    '>All layers</option>' +
    (MEAS.keys || []).map(k => '<option value="' + k + '"' + (sc === k ? ' selected' : '') + '>' + esc(measLayName(k)) + '</option>').join('') +
    '</select></label>' +
    '<label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>Show lines with Iγ ≥ (%)</span><input type="number" data-mminI step="any" min="0" value="' +
    m.minI +
    '" style="width:80px"></label><button class="btn sm" data-mselall="strong">Strongest line of each isotope</button><button class="btn sm" data-mselall="all">Select all</button><button class="btn sm" data-mselall="none">Clear selection</button><span class="hint" id="measSelCnt"></span></div>' +
    (MEAS.list.length
      ? '<div class="tblwrap" style="max-height:60vh;overflow:auto"><table class="tbl"><thead><tr><th></th><th>Eγ (keV)</th><th>Iγ (%)</th><th>ε peak</th></tr></thead><tbody>' +
        rows +
        '</tbody></table></div>'
      : '<p class="warn">No nuclide with discrete γ lines. Load ENDF-6 decay data (e.g. ENDF/B-VIII.0 decay sublibrary) in the activation panel and select “Uploaded ENDF-6 decay data”.</p>') +
    '</div>';
  const c = $('#measSelCnt', w);
  if (c)
    c.textContent =
      measSelCombos(m).length + ' selected' + (sc === null ? ' (all layers)' : ' (' + measLayName(sc) + ' only; selections in other layers are ignored)');
}
let _measTmr = null;
function measLater(ms) {
  clearTimeout(_measTmr);
  _measTmr = setTimeout(() => {
    if (MEAS.res) measRun();
  }, ms || 350);
}
function measOnCfg() {
  // the Ge/DAQ configuration changed
  const w = $('#measWin');
  if (!w || w.hidden) return;
  const g = gamSet();
  const i = $('#measInfo');
  if (i) i.innerHTML = measInfoHtml();
  measBkgDraw();
  if (MEAS.res) measLater(500);
}
document.addEventListener('click', async ev => {
  if (ev.target.closest('[data-mreset]')) {
    irrSet();
    delete S.irr.meas;
    MEAS.res = null;
    MEAS.z = {};
    saveLocal();
    measOpen();
    return;
  }
  if (ev.target.closest('[data-measopen]')) {
    measOpen();
    return;
  }
  const mw = $('#measWin'),
    sw = $('#measSelWin');
  if (ev.target.closest('[data-mclose]') || (ev.target.id === 'measWin' && _mdownId === 'measWin')) {
    mw.hidden = true;
    measTipEl().hidden = true;
    return;
  }
  if (ev.target.closest('[data-mselclose]') || (ev.target.id === 'measSelWin' && _mdownId === 'measSelWin')) {
    sw.hidden = true;
    measSelSummary();
    if (MEAS.res) measLater(50);
    return;
  }
  if (ev.target.closest('[data-msel]')) {
    sw.hidden = false;
    measSelRender();
    return;
  }
  const sa = ev.target.closest('[data-mselall]');
  if (sa) {
    const m = measSet(),
      s = irrSet().meas,
      sc = measSelScope(),
      eff = gamEffFn(gamSet()),
      inSc = n => sc === null || n.li === sc;
    MEAS.list.filter(inSc).forEach(n =>
      n.lines.forEach(l => {
        delete m.sel[n.li + '|' + n.id + '|' + l[0]];
      })
    );
    if (sa.dataset.mselall === 'strong')
      MEAS.list.filter(inSc).forEach(n => {
        const l = n.lines.find(l => l[0] >= 50 && eff(l[0]) > 0) || n.lines[0];
        m.sel[n.li + '|' + n.id + '|' + l[0]] = true;
      });
    if (sa.dataset.mselall === 'all') {
      const minI = Math.max(0, +m.minI || 0) / 100;
      MEAS.list.filter(inSc).forEach(n =>
        n.lines
          .filter(l => l[1] >= minI)
          .slice(0, 14)
          .forEach(l => {
            m.sel[n.li + '|' + n.id + '|' + l[0]] = true;
          })
      );
    }
    s.sel = m.sel;
    saveLocal();
    measSelChanged();
    measSelRender();
    return;
  }
  if (ev.target.closest('[data-mrun]')) {
    measRun();
    return;
  }
  if (!mw || mw.hidden) return;
  if (ev.target.closest('[data-mgzr]')) {
    MEAS.z = {};
    measDraw();
    return;
  }
  const ga = ev.target.closest('[data-mgall]');
  if (ga) {
    (MEAS.names || []).forEach(n => {
      MEAS.hidden[n] = ga.dataset.gstate === 'off';
    });
    measDraw();
    return;
  }
  const lg = ev.target.closest('[data-mleg]');
  if (lg) {
    const n = lg.dataset.mleg;
    MEAS.hidden[n] = !MEAS.hidden[n];
    measDraw();
    return;
  }
  if (ev.target.closest('[data-mpick]')) {
    MEAS.pick = !MEAS.pick;
    measDraw();
    return;
  }
  const sl = ev.target.closest('[data-mslice]');
  if (sl && MEAS.res) {
    if (sl.dataset.mslice === 'whole') measSetSlice(0, MEAS.res.T);
    else {
      MEAS.slice = null;
      MEAS.spec = null;
      measDraw();
    }
    return;
  }
  if (ev.target.closest('[data-mqtyrate]')) {
    const q = $('#measWin [data-mk=qty]');
    if (q) {
      q.value = 'rate';
      q.dispatchEvent(new Event('change', { bubbles: true }));
    }
    return;
  }
  if (ev.target.closest('[data-bgzr]')) {
    MEAS.bz = {};
    measBkgDraw();
    return;
  }
  const bga = ev.target.closest('[data-bgall]');
  if (bga) {
    MEAS.bh = { background: bga.dataset.gstate === 'off' };
    measBkgDraw();
    return;
  }
  const blg = ev.target.closest('[data-bleg]');
  if (blg) {
    const H = MEAS.bh || (MEAS.bh = {});
    H[blg.dataset.bleg] = !H[blg.dataset.bleg];
    measBkgDraw();
    return;
  }
  const bsv = ev.target.closest('[data-mbsave]');
  if (bsv) {
    measBkgSave(bsv.dataset.mbsave);
    return;
  }
  if (ev.target.closest('[data-sgzr]')) {
    MEAS.zs = {};
    measSpecDraw();
    return;
  }
  const sga = ev.target.closest('[data-sgall]');
  if (sga) {
    (MEAS.spnames || []).forEach(n => {
      (MEAS.sh || (MEAS.sh = {}))[n] = sga.dataset.gstate === 'off';
    });
    measSpecDraw();
    return;
  }
  const slg = ev.target.closest('[data-sleg]');
  if (slg) {
    const n = slg.dataset.sleg,
      H = MEAS.sh || (MEAS.sh = {});
    H[n] = !H[n];
    measSpecDraw();
    return;
  }
  const sps = ev.target.closest('[data-mspsave]');
  if (sps) {
    measSpecSave(sps.dataset.mspsave);
    return;
  }
  const sv = ev.target.closest('[data-msave]');
  if (sv) measSave(sv.dataset.msave);
});
document.addEventListener('change', ev => {
  const sk = ev.target.closest && ev.target.closest('#measSelWin [data-msk]');
  if (sk) {
    const m = measSet();
    if (sk.checked) m.sel[sk.dataset.msk] = true;
    else delete m.sel[sk.dataset.msk];
    irrSet().meas.sel = m.sel;
    saveLocal();
    measSelChanged();
    const c = $('#measSelCnt');
    if (c)
      c.textContent =
        measSelCombos(m).length +
        ' selected' +
        (measSelScope() === null ? ' (all layers)' : ' (' + measLayName(measSelScope()) + ' only; selections in other layers are ignored)');
    return;
  }
  const lf = ev.target.closest && ev.target.closest('#measSelWin [data-mlayf]');
  if (lf) {
    irrSet().meas.layF = lf.value;
    saveLocal();
    measLayerChanged();
    return;
  }
  const mi = ev.target.closest && ev.target.closest('#measSelWin [data-mminI]');
  if (mi) {
    irrSet().meas.minI = +mi.value;
    saveLocal();
    measSelRender();
    return;
  }
  const by = ev.target.closest && ev.target.closest('#measWin [data-bgyl]');
  if (by) {
    measSet().byl = by.value;
    saveLocal();
    MEAS.bz = Object.assign({}, MEAS.bz, { ymin: '', ymax: '' });
    measBkgDraw();
    return;
  }
  const bx = ev.target.closest && ev.target.closest('#measWin [data-bgxl]');
  if (bx) {
    measSet().bxl = bx.value;
    saveLocal();
    MEAS.bz = Object.assign({}, MEAS.bz, { xmin: '', xmax: '' });
    measBkgDraw();
    return;
  }
  const sy = ev.target.closest && ev.target.closest('#measWin [data-sgyl]');
  if (sy) {
    measSet().syl = sy.value;
    saveLocal();
    MEAS.zs = Object.assign({}, MEAS.zs, { ymin: '', ymax: '' });
    measSpecDraw();
    return;
  }
  const st = ev.target.closest && ev.target.closest('#measWin [data-mst]');
  if (st && MEAS.res) {
    const f = measAx().f,
      a = +$('#measWin [data-mst="1"]').value,
      b = +$('#measWin [data-mst="2"]').value;
    if (isFinite(a) && isFinite(b) && b > a) measSetSlice(a * f, b * f);
    return;
  }
  const yl = ev.target.closest && ev.target.closest('#measWin [data-mgyl]');
  if (yl) {
    irrSet().meas.ylog = yl.value;
    MEAS.z = Object.assign({}, MEAS.z, { ymin: '', ymax: '' });
    saveLocal();
    measDraw();
    return;
  }
  const el = ev.target.closest && ev.target.closest('#measWin [data-mk]');
  if (!el) return;
  const k = el.dataset.mk,
    s = irrSet().meas || (measSet(), irrSet().meas);
  s[k] = el.type === 'checkbox' ? el.checked : el.type === 'number' ? +el.value : el.value;
  saveLocal();
  const plotOnly = ['showTot', 'showNet', 'err', 'xbar', 'th', 'qty', 'xU'].includes(k);
  if (k === 'qty' || k === 'xU') MEAS.z = {};
  if (plotOnly) {
    measDraw();
    return;
  }
  if (k === 'layF') {
    measLayerChanged();
    return;
  }
  if (k === 'markSel') {
    measSpecDraw();
    return;
  }
  if (k === 'emax') {
    measSpecCompute();
    measDraw();
    return;
  }
  if (k === 'tirr' || k === 'tcool' || k === 'tU' || k === 'tcU' || k === 'I' || k === 'IU' || k === 'q') {
    MEAS.z = {};
    measPrep().then(() => {
      measSelSummary();
      measBkgDraw();
      if (MEAS.res) measLater(50);
    });
    return;
  }
  if (k === 'bmode' || k === 'brate' || k === 'bL' || k === 'bEc' || k === 'bpk') {
    measBkgDraw();
  }
  if (k === 'roiK') {
    const i = $('#measInfo');
    if (i) i.innerHTML = measInfoHtml();
  }
  if (k === 'newSeed') return;
  if (MEAS.res) measLater(300);
});

// ---- γ-ray spectrum in the detector for a time range of the plot ----
function measSnap(R, t1, t2) {
  // limits in seconds → nearest edges of the time bins
  const e = R.edges,
    nr = t => {
      let b = 0;
      for (let k = 0; k < e.length; k++) if (Math.abs(e[k] - t) < Math.abs(e[b] - t)) b = k;
      return b;
    };
  let i1 = nr(Math.max(0, Math.min(t1, t2))),
    i2 = nr(Math.min(R.T, Math.max(t1, t2)));
  if (i2 <= i1) {
    if (i2 < e.length - 1) i2 = i1 + 1;
    else i1 = i2 - 1;
  }
  return { t1: e[i1], t2: e[i2], i1, i2 };
}
function measSetSlice(t1, t2) {
  const R = MEAS.res;
  if (!R) return;
  MEAS.slice = measSnap(R, t1, t2);
  measSpecCompute();
  measDraw();
}
function measSpecCompute() {
  const R = MEAS.res;
  MEAS.spec = null;
  if (!R || !MEAS.slice) return;
  const sl = (MEAS.slice = measSnap(R, MEAS.slice.t1, MEAS.slice.t2));
  const m = measSet(),
    g = R.g,
    ov = { peak: gamEffFn(g), tot: gamEffTotFn(g) };
  let emax = Math.round(+m.emax || 0);
  if (!(emax > 0)) {
    let mx = 0;
    for (const li of R.keys) {
      const so = R.sols[li];
      so.info.forEach(e => {
        if (gamHas(e)) for (const l of gamAll(e)) if (l[0] > mx) mx = l[0];
      });
    }
    emax = Math.max(3500, Math.ceil(mx * 1.05));
  }
  emax = Math.min(Math.max(emax, 100), 10000);
  const nch = emax,
    peaks = new Float64Array(nch),
    cont = new Float64Array(nch),
    dbl = new Float64Array(nch);
  let live = 0;
  for (const li of R.keys) {
    const { D, t1e, t2e } = gamDecays(R.sols[li], R.prods[li], R.Tirr, sl.t1, sl.t2);
    live = Math.max(live, t2e - t1e);
    const S = geSimSpectrum(R.sols[li], D, g, nch, ov);
    for (let c = 0; c < nch; c++) {
      peaks[c] += S.peaks[c];
      cont[c] += S.cont[c];
      dbl[c] += S.dbl[c];
    }
  }
  let wa = 0,
    wb = 0,
    ta = 0,
    tb = 0;
  for (let k = sl.i1; k < sl.i2; k++) {
    const lf = R.lossF ? R.lossF[k] : R.m.pile ? 1 - R.puLoss[k] : 1,
      w = R.live[k] * R.C[0].rate[k];
    wa += w * lf; // signal: weighted by the pulse rate of the bin
    wb += w;
    ta += R.live[k] * lf; // background: constant rate, weighted by time
    tb += R.live[k];
  }
  const pf = wb > 0 ? wa / wb : 1,
    pfB = tb > 0 ? ta / tb : 1,
    sc = pf,
    bk = R.bk,
    mu = [],
    tot = new Float64Array(nch),
    bkg = new Float64Array(nch),
    bkgE = new Float64Array(nch),
    expd = new Float64Array(nch),
    rnd = measRng((R.seed ^ (sl.i1 * 73856093) ^ (sl.i2 * 19349663)) >>> 0);
  for (let c = 0; c < nch; c++) {
    peaks[c] *= sc;
    cont[c] *= sc;
    dbl[c] *= sc;
    bkgE[c] = bk.roi(c, c + 1) * live * pfB; // the background counts share the pile-up and dead-time losses
    bkg[c] = R.m.noise ? measPois(bkgE[c], rnd) : bkgE[c];
    expd[c] = peaks[c] + cont[c] + dbl[c] + bkgE[c];
    tot[c] = (R.m.noise ? measPois(peaks[c] + cont[c] + dbl[c], rnd) : peaks[c] + cont[c] + dbl[c]) + bkg[c];
  } // background (any model) and signal fluctuate independently; the total is their sum
  MEAS.spec = { nch, peaks, cont, dbl, bkg, bkgE, expd, tot, live, pf, pfB, t1: sl.t1, t2: sl.t2 };
}
function measMarks(R) {
  // marks of the selected γ lines; lines closer than 0.6 keV (e.g. the 511 keV of many isotopes) share one mark and one label
  const u = Array.from(new Map(R.C.map(c => [c.li + '|' + c.id + '|' + c.E, c])).values()).sort((a, b) => a.E - b.E),
    grp = [];
  for (const c of u) {
    const q = grp[grp.length - 1];
    if (q && c.E - q.E1 < 0.6) {
      q.cs.push(c);
      q.E1 = c.E;
    } else grp.push({ cs: [c], E0: c.E, E1: c.E });
  }
  return grp.map(q => {
    const nm = Array.from(new Set(q.cs.map(measNm))),
      E = q.cs.reduce((a, c) => a + c.E, 0) / q.cs.length;
    return {
      x: E,
      label: (nm.length > 2 ? nm.slice(0, 2).join(', ') + ' +' + (nm.length - 2) : nm.join(', ')) + ' ' + fmt(E, 6),
      color: scol(R.C.indexOf(q.cs[0])),
      dash: '2 3',
      w: 1.3,
      op: 0.85,
      rot: true
    };
  });
}
function measSpecOpts() {
  const Sp = MEAS.spec,
    R = MEAS.res,
    H = MEAS.sh || (MEAS.sh = {}),
    edges = new Float64Array(Sp.nch + 1);
  for (let i = 0; i <= Sp.nch; i++) edges[i] = i;
  const z = MEAS.zs || (MEAS.zs = {}),
    ser = [
      { name: 'total spectrum (with fluctuations)', type: 'hist', edges, y: Sp.tot, color: 'var(--ink)', w: 1.4 },
      { name: 'full-energy peaks', type: 'hist', edges, y: Sp.peaks, color: 'var(--accent)', w: 1.1 },
      { name: 'single Compton + escape peaks', type: 'hist', edges, y: Sp.cont, color: 'var(--warn)', w: 1.1 },
      { name: 'double Compton (2 scatterings)', type: 'hist', edges, y: Sp.dbl, color: 'var(--bad)', w: 1.1 }
    ];
  if (R.m.bmode !== 'none') ser.push({ name: 'background (with fluctuations)', type: 'hist', edges, y: Sp.bkg, color: 'var(--ink2)', w: 1.1 });
  if (!MEAS.shInit) {
    MEAS.shInit = true;
    ['full-energy peaks', 'single Compton + escape peaks', 'double Compton (2 scatterings)'].forEach(n => {
      H[n] = true;
    });
  } // components off by default
  ser.forEach(se => {
    se.hidden = !!H[se.name];
  });
  return {
    series: ser,
    xlog: false,
    ylog: (measSet().syl || 'log') === 'log',
    xmin: z.xmin || '0',
    xmax: z.xmax || '',
    ymin: z.ymin || '',
    ymax: z.ymax || '',
    vlines: measSet().markSel ? measMarks(R) : [],
    xlabel: 'Energy (keV, 1 keV per channel)',
    ylabel: 'Counts per keV in the time range',
    title:
      'γ-ray spectrum in the detector · ' +
      (R.keys.length > 1 || (MEAS.keys || []).length > 1
        ? R.keys.length === (MEAS.keys || []).length
          ? 'all layers'
          : 'layer ' + (R.keys[0] + 1)
        : 'layer 1') +
      ' · ' +
      fmtHl(Sp.t1) +
      ' – ' +
      fmtHl(Sp.t2)
  };
}
function measSpecDraw() {
  const box = $('#measSpec');
  if (!box) return;
  const R = MEAS.res,
    m = measSet();
  if (!R) {
    box.innerHTML = '';
    return;
  }
  const ax = measAx(),
    uT = ax.u,
    Sp = MEAS.spec,
    sl = MEAS.slice;
  const head =
    '<div class="btns"><label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>Time range from</span><input type="number" step="any" data-mst="1" value="' +
    (sl ? +(sl.t1 / ax.f).toPrecision(6) : '') +
    '" style="width:110px"><span>to</span><input type="number" step="any" data-mst="2" value="' +
    (sl ? +(sl.t2 / ax.f).toPrecision(6) : '') +
    '" style="width:110px"><span>' +
    uT +
    '</span></label><button class="btn sm" data-mslice="whole">Whole run</button><button class="btn sm" data-mslice="clear">Remove</button><label class="chk"><input type="checkbox" data-mk="markSel"' +
    (m.markSel ? ' checked' : '') +
    '><span>mark the selected γ lines</span></label><label class="fld" style="flex-direction:row;align-items:center;gap:6px"><span>Upper energy (keV, 0 = automatic)</span><input type="number" step="1" min="0" data-mk="emax" value="' +
    (+m.emax || 0) +
    '" style="width:90px"></label></div>';
  if (!sl || !Sp) {
    box.innerHTML =
      '<details class="cw" open><summary><b>γ-ray spectrum in the detector</b></summary><p class="hint">Press “Select a time range → γ spectrum” and drag over the plot, or type the limits below. The spectrum contains <b>all</b> γ lines of the nuclides of the chosen layer (all layers if so chosen), not only the selected lines.</p>' +
      head +
      '</details>';
    return;
  }
  const o = measSpecOpts(),
    r = renderChart(o),
    names = o.series.map(x => [x.name, x.color]);
  MEAS.spo = o;
  MEAS.spg = r.geom;
  MEAS.spnames = names.map(x => x[0]);
  const leg = gamRen(gamLegend(names, MEAS.sh, 'sp'), 's'),
    stop = gamRen(gamTop('sp', measSet().syl || 'log', gamZoomed(MEAS.zs)), 's');
  let tp = 0,
    tb = 0;
  for (let c = 0; c < Sp.nch; c++) {
    tp += Sp.tot[c];
    tb += Sp.bkg[c];
  }
  const info =
    'Time range ' +
    fmt(sl.t1 / ax.f, 6) +
    ' – ' +
    fmt(sl.t2 / ax.f, 6) +
    ' ' +
    uT +
    ' (' +
    fmtHl(Sp.live) +
    ' live) · <b>total counts in the spectrum ' +
    fmt(tp, 7) +
    '</b> (± ' +
    fmt(Math.sqrt(tp), 3) +
    ' statistical; background ' +
    fmt(tb, 5) +
    ') · <b>average count rate ' +
    fmt(tp / Math.max(Sp.live, 1e-12), 6) +
    ' counts/s</b> (background ' +
    fmt(tb / Math.max(Sp.live, 1e-12), 4) +
    ' counts/s) · mean pile-up + dead-time loss ' +
    fmt((1 - Sp.pf) * 100, 3) +
    ' %' +
    ' · resolution FWHM(1 MeV) ' +
    fmt(gamFwhm(1000, R.g), 3) +
    ' keV · efficiency model ' +
    { ge: 'Ge geometry', curve: 'curve ε(E)', const: 'constant' }[R.g.effMode] +
    ' · background model ' +
    { none: 'A (none)', exp: 'B (exponential)', cont: 'C (continuum + peaks)' }[R.m.bmode] +
    (R.m.noise ? ' · Poisson fluctuations' : ' · no fluctuations') +
    '.';
  box.innerHTML =
    '<details class="cw" open><summary><b>γ-ray spectrum in the detector</b> — all γ lines of ' +
    (R.keys.length > 1 ? 'layers ' + R.keys.map(k => k + 1).join(' + ') : 'layer ' + (R.keys[0] + 1)) +
    ': peaks, single and double Compton, escape peaks, background · drag a box to zoom, double-click resets</summary>' +
    head +
    '<p class="hint">' +
    info +
    '</p>' +
    stop +
    '<div class="chart" id="measSpecChart">' +
    r.svg +
    '</div>' +
    leg +
    pcSave({ attr: 'data-mspsave', title: 'Save spectrum', formats: ['png', 'pdf', 'svg', 'csv', 'root'], msg: '<span class="hint" id="measSpMsg"></span>' }) +
    '</details>';
  gamZoom($('#measSpecChart svg', box), r.geom, 'zs', MEAS, measSpecDraw);
  specTip($('#measSpecChart svg', box), r.geom, () => MEAS.spo);
}
async function measSpecSave(kind) {
  const m = $('#measSpMsg'),
    Sp = MEAS.spec;
  if (!Sp) return;
  let tpS = 0;
  for (let c = 0; c < Sp.nch; c++) tpS += Sp.tot[c];
  if (m) m.textContent = 'Preparing…';
  const o = measSpecOpts(),
    r = renderChart(o),
    div = document.createElement('div');
  div.innerHTML = r.svg;
  const ctx = { o, geom: r.geom, svg: div.querySelector('svg'), base: 'cpt_measurement_spectrum', key: 'meas-sp' },
    old = _chart;
  _chart = ctx;
  try {
    const px = expPx(),
      base = ctx.base;
    if (kind === 'png') await saveFile(base + '.png', await canvasBlob(await svgToCanvas(px), 'image/png'), 'PNG', '#measSpMsg');
    else if (kind === 'pdf') await saveFile(base + '.pdf', await makePdf(await svgToCanvas(Math.max(px, 2400)), ctx.geom), 'PDF', '#measSpMsg');
    else if (kind === 'svg') await saveFile(base + '.svg', exportSvgString(), 'SVG', '#measSpMsg');
    else if (kind === 'root') await saveFile(base + '.C.txt', rootMacroOf(), 'ROOT macro', '#measSpMsg');
    else {
      let t =
        '# simulated gamma-ray spectrum in the detector, time range ' +
        Sp.t1 +
        ' - ' +
        Sp.t2 +
        ' s (live ' +
        Sp.live +
        ' s), total counts ' +
        tpS +
        ', average rate ' +
        tpS / Sp.live +
        ' counts/s\nchannel_low_keV,total_with_fluctuations,expected_total,full_energy_peaks,single_compton_escape,double_compton,background_with_fluctuations,expected_background\n';
      for (let c = 0; c < Sp.nch; c++) t += [c, Sp.tot[c], Sp.expd[c], Sp.peaks[c], Sp.cont[c], Sp.dbl[c], Sp.bkg[c], Sp.bkgE[c]].join(',') + '\n';
      await saveFile(base + '.csv', t, 'CSV', '#measSpMsg');
    }
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  } finally {
    _chart = old;
  }
}

// ---- export ----
function measCsv(sep) {
  const R = MEAS.res,
    ax = measAx(),
    f = ax.f,
    uT = ax.u;
  let t =
    '# Charged Particle Tracker — simulated HPGe measurement\n# irradiation ' +
    R.Tirr +
    ' s, cooling ' +
    R.Tc +
    ' s, bins ' +
    R.dtI +
    ' s (irradiation) and ' +
    R.dtC +
    ' s (cooling), ROI = E ± ' +
    R.m.roiK +
    ' × FWHM, background model ' +
    R.m.bmode +
    ' (' +
    R.m.brate +
    ' cps), random seed ' +
    R.seed +
    ', beam-current uncertainty ' +
    R.m.dI +
    ' %, efficiency uncertainty ' +
    R.m.dE +
    ' %\n';
  const head = [
    'isotope',
    'layer',
    'E_keV',
    'I_gamma',
    't_start_s',
    't_end_s',
    't_centre_s',
    't_centre_' + uT,
    'total_counts_ROI',
    'background_counts',
    'net_counts',
    'sigma_net_total',
    'sigma_net_percent',
    'sigma_stat',
    'sigma_syst',
    'pileup_percent',
    'detector_rate_cps',
    'total_rate_cps',
    'background_rate_cps',
    'net_rate_cps',
    'sigma_net_rate_cps',
    'expected_total',
    'expected_net',
    'ROI_low_keV',
    'ROI_high_keV'
  ];
  t += (sep === ',' ? '' : '# ') + head.join(sep) + '\n';
  for (const c of R.C)
    for (let k = 0; k < R.N; k++)
      t +=
        [
          dkNuc(c.id),
          c.li + 1,
          c.E,
          c.I,
          c.t1[k],
          c.t2[k],
          c.t[k],
          c.t[k] / f,
          c.tot[k],
          c.bkg[k],
          c.net[k],
          c.sig[k],
          Math.abs(c.net[k]) > 0 ? (c.sig[k] / Math.abs(c.net[k])) * 100 : '',
          c.sigStat[k],
          c.sigSys[k],
          c.pu[k],
          c.rate[k],
          c.tot[k] / (c.t2[k] - c.t1[k]),
          c.bkg[k] / (c.t2[k] - c.t1[k]),
          c.net[k] / (c.t2[k] - c.t1[k]),
          c.sig[k] / (c.t2[k] - c.t1[k]),
          c.thTot[k],
          c.thNet[k],
          c.lo,
          c.hi
        ].join(sep) + '\n';
  return t;
}
function measRootHists() {
  const R = MEAS.res,
    ax = measAx(),
    f = ax.f,
    uT = ax.u,
    edges = new Float64Array(R.N + 1);
  for (let k = 0; k < R.N; k++) edges[k] = R.C[0].t1[k] / f;
  edges[R.N] = R.C[0].t2[R.N - 1] / f;
  const out = [],
    A = asciiName;
  R.C.forEach(c => {
    const base = 'L' + (c.li + 1) + '_' + A(dkNuc(c.id)).replace(/[^A-Za-z0-9]+/g, '_') + '_' + String(Math.round(c.E * 10) / 10).replace('.', 'p') + 'keV',
      tt = 'layer ' + (c.li + 1) + ' ' + A(dkNuc(c.id)) + ' ' + c.E + ' keV, ROI ' + fmt(c.lo, 6) + '-' + fmt(c.hi, 6) + ' keV',
      xt = 'Time (' + uT + ')';
    out.push({
      name: 'h_' + base + '_total',
      title: tt + ' | total counts in ROI',
      xt,
      yt: 'counts per bin',
      edges,
      vals: Float64Array.from(c.tot),
      err: Float64Array.from(c.tot, v => Math.sqrt(v))
    });
    out.push({
      name: 'h_' + base + '_background',
      title: tt + ' | background counts in ROI',
      xt,
      yt: 'counts per bin',
      edges,
      vals: Float64Array.from(c.bkg),
      err: Float64Array.from(c.bkg, v => Math.sqrt(v))
    });
    out.push({
      name: 'h_' + base + '_net',
      title: tt + ' | net counts (total - background)',
      xt,
      yt: 'counts per bin',
      edges,
      vals: Float64Array.from(c.net),
      err: Float64Array.from(c.sig)
    });
    out.push({
      name: 'h_' + base + '_expected_net',
      title: tt + ' | expected net counts (no noise)',
      xt,
      yt: 'counts per bin',
      edges,
      vals: Float64Array.from(c.thNet)
    });
    const lvv = Float64Array.from(c.t1, (v, k) => Math.max(c.t2[k] - v, 1e-12)),
      rt = a => Float64Array.from(a, (v, k) => v / lvv[k]);
    out.push({ name: 'h_' + base + '_total_rate', title: tt + ' | total count rate in ROI', xt, yt: 'counts/s', edges, vals: rt(c.tot) });
    out.push({ name: 'h_' + base + '_background_rate', title: tt + ' | background count rate in ROI', xt, yt: 'counts/s', edges, vals: rt(c.bkg) });
    out.push({ name: 'h_' + base + '_net_rate', title: tt + ' | net (peak) count rate', xt, yt: 'counts/s', edges, vals: rt(c.net), err: rt(c.sig) });
    out.push({ name: 'h_' + base + '_pileup', title: tt + ' | pile-up level', xt, yt: 'percent', edges, vals: Float64Array.from(c.pu) });
  });
  return out;
}
async function measSave(kind) {
  const m = $('#measSaveMsg'),
    R = MEAS.res;
  if (!R) return;
  if (m) m.textContent = 'Preparing…';
  const ctx = measCtx(),
    base = ctx.base,
    old = _chart;
  _chart = ctx;
  try {
    const px = expPx();
    if (kind === 'png') await saveFile(base + '.png', await canvasBlob(await svgToCanvas(px), 'image/png'), 'PNG', '#measSaveMsg');
    else if (kind === 'pdf') await saveFile(base + '.pdf', await makePdf(await svgToCanvas(Math.max(px, 2400)), ctx.geom), 'PDF', '#measSaveMsg');
    else if (kind === 'svg') await saveFile(base + '.svg', exportSvgString(), 'SVG', '#measSaveMsg');
    else if (kind === 'csv') await saveFile(base + '.csv', measCsv(','), 'CSV', '#measSaveMsg');
    else if (kind === 'txt') await saveFile(base + '.txt', measCsv('\t'), 'ASCII', '#measSaveMsg');
    else if (kind === 'root') await saveFile(base + '.C.txt', rootMacroOf(), 'ROOT macro', '#measSaveMsg');
    else if (kind === 'rootbin') {
      const hs = measRootHists(),
        rn = base + '.root',
        rb = RootIO.file(rn, hs),
        inClaude = !!(window.claude && claude.use && (await claude.use('downloads')));
      if (inClaude) await saveFile(base + '_histograms.zip', zipStore(rn, rb), 'ROOT file', '#measSaveMsg');
      else await saveFile(rn, new Blob([rb], { type: 'application/octet-stream' }), 'ROOT file', '#measSaveMsg');
      if (m) m.textContent = 'Saved ' + (inClaude ? base + '_histograms.zip (containing ' + rn + ')' : rn) + ' with ' + hs.length + ' histograms.';
    }
  } catch (e) {
    if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
  } finally {
    _chart = old;
  }
}
