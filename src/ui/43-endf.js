/* ====== ENDF-6 reader for charged-particle (alpha) sublibraries, e.g. TENDL ======
   Reads MF3 (cross sections, TAB1), the product yields of MF6 (TAB1 of each product) and the MF6 energy-angle
   distributions LAW 1 (Legendre / Kalbach-Mann), LAW 2 (two-body) and LAW 6 (N-body); other laws are skipped by record structure.
   Units as in ENDF: energies in eV on file -> MeV here; cross sections in barn. */

// -- fixed-format helpers (11-character fields, 6 per line; MAT 67-70, MF 71-72, MT 73-75) --
function endfNum(s) {
  s = s.trim();
  if (!s) return 0;
  let k = -1;
  for (let i = s.length - 1; i > 0; i--) {
    const c = s[i];
    if ((c === '+' || c === '-') && s[i - 1] !== 'e' && s[i - 1] !== 'E') {
      k = i;
      break;
    }
  }
  if (k > 0) {
    const m = s.slice(0, k),
      e = s.slice(k);
    const v = parseFloat(m + 'e' + e);
    if (isFinite(v)) return v;
  }
  const v = parseFloat(s);
  return isFinite(v) ? v : 0;
}
function endfFields(line) {
  const o = [];
  for (let i = 0; i < 6; i++) o.push(endfNum(line.substr(i * 11, 11)));
  return o;
}

// -- TAB1 tabulated function with interpolation regions (INT: 1 histogram, 2 lin-lin, 3 lin-log, 4 log-lin, 5 log-log) --
function tab1Eval(t, x) {
  const X = t.x,
    Y = t.y,
    n = X.length;
  if (!n) return 0;
  if (x < X[0]) return 0;
  if (x > X[n - 1]) return 0;
  if (x === X[n - 1]) return Y[n - 1];
  let lo = 0,
    hi = n - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (X[m] > x) hi = m;
    else lo = m;
  }
  let law = 2;
  for (let r = 0; r < t.nbt.length; r++) {
    if (hi + 1 <= t.nbt[r]) {
      law = t.int[r];
      break;
    }
    law = t.int[r];
  }
  const x0 = X[lo],
    x1 = X[hi],
    y0 = Y[lo],
    y1 = Y[hi];
  if (x1 === x0) return y1;
  const f = (x - x0) / (x1 - x0);
  switch (law) {
    case 1:
      return y0;
    case 3:
      return y0 + (y1 - y0) * (Math.log(x / x0) / Math.log(x1 / x0));
    case 4:
      return y0 > 0 && y1 > 0 ? Math.exp(Math.log(y0) + (Math.log(y1) - Math.log(y0)) * f) : y0 + (y1 - y0) * f;
    case 5:
      return y0 > 0 && y1 > 0 && x0 > 0 ? Math.exp(Math.log(y0) + (Math.log(y1) - Math.log(y0)) * (Math.log(x / x0) / Math.log(x1 / x0))) : y0 + (y1 - y0) * f;
    default:
      return y0 + (y1 - y0) * f;
  }
}

// Reader over the lines of one MF/MT section
class EndfCursor {
  constructor(lines, a, b) {
    this.L = lines;
    this.i = a;
    this.end = b;
  }
  more() {
    return this.i < this.end;
  }
  cont() {
    if (this.i >= this.end) throw new Error('ENDF: unexpected end of section');
    return endfFields(this.L[this.i++]);
  }
  tab1() {
    // CONT + interpolation pairs + data pairs
    const c = this.cont(),
      NR = c[4] | 0,
      NP = c[5] | 0,
      nbt = [],
      intr = [];
    let vals = [],
      need = NR * 2;
    while (vals.length < need) {
      const f = this.cont();
      for (const v of f) vals.push(v);
    }
    for (let r = 0; r < NR; r++) {
      nbt.push(vals[2 * r] | 0);
      intr.push(vals[2 * r + 1] | 0);
    }
    vals = [];
    need = NP * 2;
    while (vals.length < need) {
      const f = this.cont();
      for (const v of f) vals.push(v);
    }
    const x = new Float64Array(NP),
      y = new Float64Array(NP);
    for (let p = 0; p < NP; p++) {
      x[p] = vals[2 * p];
      y[p] = vals[2 * p + 1];
    }
    return { c, nbt, int: intr, x, y };
  }
  list() {
    const c = this.cont(),
      n = c[4] | 0,
      v = [];
    for (let k = 0, m = Math.ceil(n / 6); k < m; k++) for (const x of this.cont()) v.push(x);
    v.length = n;
    return { c, v: Float64Array.from(v) };
  }
  skipTab1() {
    // TAB1 record without storing the data (CONT + interpolation lines + data lines)
    const c = this.cont(),
      n = Math.ceil((2 * (c[4] | 0)) / 6) + Math.ceil((2 * (c[5] | 0)) / 6);
    if (this.i + n > this.end) throw new Error('ENDF: unexpected end of section');
    this.i += n;
    return c;
  }
  skipTab2() {
    const c = this.cont(),
      NR = c[4] | 0;
    for (let k = 0, n = Math.ceil((NR * 2) / 6); k < n; k++) this.cont();
    return c;
  }
  skipList() {
    const c = this.cont(),
      NPL = c[4] | 0;
    for (let k = 0, n = Math.ceil(NPL / 6); k < n; k++) this.cont();
    return c;
  }
}

// -- MF6 energy-angle distributions (energies -> MeV, f in 1/MeV) --
// LAW 1 (continuum): LANG 1 Legendre b_l, LANG 2 Kalbach-Mann (f0, r, a); other LANG values are kept as unsupported (dist = null)
function readLaw1(cur) {
  const t2 = cur.skipTab2(),
    lang = t2[2] | 0,
    lep = t2[3] | 0,
    NE = t2[5] | 0,
    E = new Float64Array(NE),
    L = [];
  for (let e = 0; e < NE; e++) {
    const { c, v } = cur.list(),
      nd = c[2] | 0,
      na = c[3] | 0,
      nep = c[5] | 0,
      per = na + 2;
    E[e] = c[1] * 1e-6;
    const n = Math.max(0, nep - nd),
      Ep = new Float64Array(n),
      f0 = new Float64Array(n),
      p = [],
      disc = { E: new Float64Array(nd), f: new Float64Array(nd) }; // disc: discrete lines (photons): energy, probability
    for (let i = 0; i < nd; i++) {
      disc.E[i] = v[i * per] * 1e-6;
      disc.f[i] = v[i * per + 1];
    }
    for (let q = 0; q < na; q++) p.push(new Float64Array(n));
    for (let i = 0; i < n; i++) {
      const o = (i + nd) * per;
      Ep[i] = v[o] * 1e-6;
      f0[i] = v[o + 1] * 1e6;
      for (let q = 0; q < na; q++) p[q][i] = lang === 1 ? v[o + 2 + q] * 1e6 : v[o + 2 + q];
    }
    L.push({ Ep, f0, p, disc });
  }
  return lang === 1 || lang === 2 ? { law: 1, lang, lep, E, L } : null;
}
// LAW 7 (laboratory angle-energy, JENDL-5): for every incident energy a set of μ values, each with a TAB1 f(E') (E' in MeV, f per MeV and per unit μ; ∫∫ f dE' dμ = 1)
function readLaw7(cur) {
  const t2 = cur.skipTab2(),
    NE = t2[5] | 0,
    E = new Float64Array(NE),
    L = [];
  for (let e = 0; e < NE; e++) {
    const t3 = cur.skipTab2(),
      NM = t3[5] | 0,
      mu = new Float64Array(NM),
      T = [];
    E[e] = t3[1] * 1e-6;
    for (let m = 0; m < NM; m++) {
      const t = cur.tab1(),
        x = new Float64Array(t.x.length),
        y = new Float64Array(t.y.length);
      for (let i = 0; i < x.length; i++) {
        x[i] = t.x[i] * 1e-6;
        y[i] = t.y[i] * 1e6;
      }
      mu[m] = t.c[1];
      T.push({ x, y, nbt: t.nbt, int: t.int });
    }
    L.push({ mu, T });
  }
  return { law: 7, E, L, frame: 'lab' };
}
// LAW 2 (two-body discrete): angular distribution per incident energy; LANG 0 Legendre a_1..a_NL, LANG 12/14 tabulated (mu, f)
function readLaw2(cur) {
  const t2 = cur.skipTab2(),
    NE = t2[5] | 0,
    E = new Float64Array(NE),
    L = [];
  for (let e = 0; e < NE; e++) {
    const { c, v } = cur.list(),
      lang = c[2] | 0;
    E[e] = c[1] * 1e-6;
    if (lang === 0) L.push({ lang: 0, a: v });
    else {
      const n = v.length >> 1,
        mu = new Float64Array(n),
        f = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        mu[i] = v[2 * i];
        f[i] = v[2 * i + 1];
      }
      L.push({ lang, mu, f });
    }
  }
  return { law: 2, E, L };
}

// split a file into MF/MT sections: Map "mf/mt" -> [firstLine, endLine) of the section body (without SEND)
function endfSections(text) {
  const lines = text.split('\n'),
    sec = new Map();
  let curKey = null,
    a = 0;
  for (let i = 0; i < lines.length; i++) {
    const ln = lines[i];
    if (ln.length < 75) continue;
    const mf = +ln.substr(70, 2),
      mt = +ln.substr(72, 3);
    const key = mf + '/' + mt;
    if (mt === 0 || mf === 0) {
      if (curKey) {
        sec.set(curKey, [a, i]);
        curKey = null;
      }
      continue;
    }
    if (key !== curKey) {
      if (curKey) sec.set(curKey, [a, i]);
      curKey = key;
      a = i;
    }
  }
  if (curKey) sec.set(curKey, [a, lines.length]);
  return { lines, sec };
}

// -- particles, reaction names and channel definitions --
const ZAP_SPECIES = { 0: 'g', 1: 'n', 1001: 'p', 1002: 'd', 1003: 't', 2003: 'h', 2004: 'a' };
// projectiles with a TENDL library (file prefix a-, p-, d-, t-, h-): symbol, atomic mass (u), binding energy of the projectile (MeV), Kalbach M_a, first MT of its own-species level series (its ground-state channel is elastic scattering)
const PROJ = {};
[
  ['p', 'p', 'proton', 'protons', 1, 1, 1.00782503, 0, 1, 600],
  ['d', 'd', 'deuteron', 'deuterons', 1, 2, 2.01410178, 2.224573, 1, 650],
  ['t', 't', 'triton', 'tritons', 1, 3, 3.01604928, 8.481798, 1, 700],
  ['h', '³He', 'helion (³He)', '³He ions', 2, 3, 3.01602932, 7.718058, 1, 750],
  ['a', 'α', 'α particle', 'α particles', 2, 4, 4.00260325, 28.295674, 0, 800]
].forEach(([letter, sym, name, plural, Z, A, u, Bp, Ma, own0]) => {
  PROJ[letter] = { letter, sym, name, plural, Z, A, awi: u / 1.00866492, Bp, Ma, own0, lib: 'TENDL-' + letter };
});
function pjOf(Z, A) {
  for (const k in PROJ) if (PROJ[k].Z === Z && PROJ[k].A === A) return PROJ[k];
  return null;
}
const SPECIES = [
  ['n', 'neutrons'],
  ['p', 'protons'],
  ['d', 'deuterons'],
  ['t', 'tritons'],
  ['h', '³He'],
  ['a', 'α particles'],
  ['g', 'photons']
];
// particle-production reactions (MF 3 = production cross section): MT 201 (z,xn), 202 (z,xγ), 203 (z,xp), 204 (z,xd), 205 (z,xt), 206 (z,x³He), 207 (z,xα)
const PROD_MT = { 201: 'n', 202: 'g', 203: 'p', 204: 'd', 205: 't', 206: 'h', 207: 'a' };
const SPEC_ZA = { n: 1, p: 1001, d: 1002, t: 1003, h: 2003, a: 2004 };
// MT -> emitted particles for the channels that have no MF6 product list; null = lumped sum / not a reaction channel
function mtEmitted(mt) {
  const m = { 16: ['n', 'n'], 17: ['n', 'n', 'n'], 22: ['n', 'a'], 28: ['n', 'p'], 102: [], 103: ['p'], 104: ['d'], 105: ['t'], 106: ['h'], 107: ['a'] };
  if (m[mt]) return m[mt];
  if (mt >= 50 && mt <= 91) return ['n']; // MT 50 = (x,n0), ground state, as in the ENDF/B charged-particle sublibraries
  if (mt >= 600 && mt <= 649) return ['p'];
  if (mt >= 650 && mt <= 699) return ['d'];
  if (mt >= 700 && mt <= 749) return ['t'];
  if (mt >= 750 && mt <= 799) return ['h'];
  if (mt >= 800 && mt <= 849) return ['a'];
  return null;
}
// non-overlapping set of channels whose sum is the non-elastic cross section (MT3, MT4 and MT103-107 are lumped sums of these)
function mtIsPrimaryChannel(mt) {
  return mt === 5 || mt === 16 || mt === 17 || mt === 22 || mt === 28 || mt === 102 || (mt >= 50 && mt <= 91) || (mt >= 600 && mt <= 849);
}
function mtName(mt, pj) {
  const P = (pj || PROJ.a).sym,
    own = (pj || PROJ.a).letter,
    OUT = { n: 'n', p: 'p', d: 'd', t: 't', h: '³He', a: 'α' },
    x = k => '(' + P + ',' + OUT[k] + (k === own ? '′' : '') + ')';
  const lump = {
    2: 'elastic (nuclear + interference)',
    3: 'non-elastic (all reactions)',
    4: x('n') + ' total, sum of MT 50-91',
    5: '(' + P + ',anything)',
    16: '(' + P + ',2n)',
    17: '(' + P + ',3n)',
    22: '(' + P + ',nα)',
    28: '(' + P + ',np)',
    102: '(' + P + ',γ)',
    201: '(' + P + ',xn) neutron production',
    202: '(' + P + ',xγ) photon production',
    203: '(' + P + ',xp) proton production',
    204: '(' + P + ',xd) deuteron production',
    205: '(' + P + ',xt) triton production',
    206: '(' + P + ',x³He) ³He production',
    207: '(' + P + ',xα) α production',
    103: x('p') + ' total',
    104: x('d') + ' total',
    105: x('t') + ' total',
    106: x('h') + ' total',
    107: x('a') + ' total'
  };
  if (lump[mt]) return lump[mt];
  const lev = (m0, m1, k) => (mt >= m0 && mt <= m1 ? (mt === m1 ? x(k) + ' continuum' : x(k) + ' to level ' + (mt - m0)) : null);
  return lev(50, 91, 'n') || lev(600, 649, 'p') || lev(650, 699, 'd') || lev(700, 749, 't') || lev(750, 799, 'h') || lev(800, 849, 'a') || 'MT ' + mt;
}

// -- parse a whole alpha file; returns { ZA, Z, A, title, channels: {mt: ch}, order: [mt...] } --
// ch = { mt, name, QM, QI, E: Float64Array (MeV), s: Float64Array (barn), sT (tab1 raw),
//        y: { n,p,d,t,h,a,g : tab1 or constant number }, res: [{ZA, y: tab1|number}], mf6: bool }
function parseEndfAlpha(text, pj) {
  pj = pj || PROJ.a;
  const { lines, sec } = endfSections(text);
  const s451 = sec.get('1/451');
  let ZA = 0,
    AWR = 0,
    title = '';
  if (s451) {
    const c = new EndfCursor(lines, s451[0], s451[1]).cont();
    ZA = Math.round(c[0]);
    AWR = c[1];
    title = lines[s451[0]].substr(0, 66).trim();
  } else {
    const c = endfFields(lines[1] || '');
    ZA = Math.round(c[0]);
    AWR = c[1];
    title = (lines[0] || '').trim();
  }
  const Z = Math.floor(ZA / 1000),
    A = ZA % 1000,
    channels = {},
    order = [];
  const toMeV = t => {
    const x = new Float64Array(t.x.length);
    for (let i = 0; i < x.length; i++) x[i] = t.x[i] * 1e-6;
    return { x, y: t.y, nbt: t.nbt, int: t.int };
  };
  for (const [key, [a, b]] of sec) {
    const [mf, mt] = key.split('/').map(Number);
    if (mf !== 3) continue;
    const cur = new EndfCursor(lines, a, b);
    cur.cont(); // HEAD
    const t = cur.tab1();
    if (!t.x.length) continue;
    const tm = toMeV(t);
    let mx = 0;
    for (let i = 0; i < t.y.length; i++) if (Math.abs(t.y[i]) > mx) mx = Math.abs(t.y[i]);
    if (mt !== 2 && mt !== 3 && !(mx > 0)) continue;
    channels[mt] = {
      mt,
      name: mtName(mt, pj),
      QM: t.c[0] * 1e-6,
      QI: t.c[1] * 1e-6,
      tab: tm,
      E: tm.x,
      s: t.y,
      y: {},
      res: [],
      mf6: false,
      primary: mtIsPrimaryChannel(mt),
      mx
    };
    order.push(mt);
  }
  order.sort((p, q) => p - q);
  // non-overlapping (primary) channels: every reaction MT below 118 and 600-849 except the totals (MT 1, 2, 3, 27, 101).
  // MT 4 and 103-107 are sums of the level channels 51-91 / 600-849: when those are present and add up to the sum (3 %), the levels are used and the sum is not;
  // when they do not add up (incomplete set of levels), the sum is used instead of the levels, so that nothing is lost or counted twice.
  const integ = ch => {
    let q = 0;
    for (let i = 0; i + 1 < ch.E.length; i++) q += 0.5 * (ch.s[i] + ch.s[i + 1]) * (ch.E[i + 1] - ch.E[i]);
    return q;
  };
  for (const mt of order) channels[mt].primary = (mt <= 117 || (mt >= 600 && mt <= 849)) && ![1, 2, 3, 27, 101].includes(mt) && mt !== pj.own0; // the ground-state channel of the projectile's own species is elastic scattering
  const lumpNotes = [];
  for (const [L, m0, m1] of [
    [4, 50, 91],
    [103, 600, 649],
    [104, 650, 699],
    [105, 700, 749],
    [106, 750, 799],
    [107, 800, 849]
  ]) {
    if (!channels[L]) continue;
    let comps = 0,
      any = false;
    for (let m = m0; m <= m1; m++)
      if (channels[m] && m !== pj.own0) {
        any = true;
        comps += integ(channels[m]);
      }
    if (!any) continue;
    const lump = integ(channels[L]),
      ok = !(lump > 0) || Math.abs(comps / lump - 1) <= 0.03;
    channels[L].primary = !ok;
    for (let m = m0; m <= m1; m++) if (channels[m]) channels[m].primary = ok;
    lumpNotes.push(
      ok
        ? 'MT ' + L + ' replaced by its levels MT ' + m0 + '-' + m1
        : 'MT ' + L + ' used instead of its levels MT ' + m0 + '-' + m1 + ' (they add up to ' + ((100 * comps) / lump).toFixed(0) + ' % of it)'
    );
  }
  // MF6: product yields
  for (const [key, [a, b]] of sec) {
    const [mf, mt] = key.split('/').map(Number);
    if (mf !== 6 || !channels[mt]) continue;
    const ch = channels[mt],
      cur = new EndfCursor(lines, a, b);
    try {
      const h = cur.cont(),
        NK = h[4] | 0,
        lct = h[3] | 0;
      let ok = true;
      ch.prods = [];
      for (let k = 0; k < NK && ok; k++) {
        const y = cur.tab1(),
          zap = Math.round(y.c[0]),
          law = y.c[3] | 0,
          yt = toMeV(y);
        const ytab = { x: yt.x, y: y.y, nbt: y.nbt, int: y.int };
        const sp = ZAP_SPECIES[zap];
        if (sp)
          ch.y[sp] = ch.y[sp] ? { sum: [...(ch.y[sp].sum || [ch.y[sp]]), ytab] } : ytab; // a product listed several times: the yields add up
        else if (zap > 1003) ch.res.push({ ZA: zap, lip: y.c[2] | 0, y: ytab }); // residual nuclei (ZAP = 1000 Z + A; includes the unbound H-4, H-5 of the JENDL files)
        let dist = null;
        if (law === 1) dist = readLaw1(cur);
        else if (law === 2) dist = readLaw2(cur);
        else if (law === 5) {
          const t2 = cur.skipTab2(),
            NE = t2[5] | 0;
          for (let e = 0; e < NE; e++) cur.skipList();
        } else if (law === 6) {
          const c6 = cur.cont();
          dist = { law: 6, apsx: c6[0], n: c6[5] | 0 };
        } else if (law === 7) dist = readLaw7(cur);
        else if (law === 0 || law === 3 || law === 4) {
          /* no data follow */
        } else ok = false; // unknown LAW: stop reading this section
        if (dist && law !== 7) dist.frame = lct === 1 ? 'lab' : lct === 3 && y.c[1] > 4.1 ? 'lab' : 'cm'; // LCT 3: CM for light products (A ≤ 4), laboratory for the heavy recoils
        ch.prods.push({ sp, zap, awp: y.c[1], ytab, dist });
      }
      ch.mf6 = true;
    } catch (e) {
      /* truncated section: keep what was read */ ch.mf6 = ch.mf6 || Object.keys(ch.y).length > 0;
    }
  }
  // particle-production channels MT 201-207 ((z,xn), (z,xγ), (z,xp), … as in the JENDL-5 α file: MF 3 = production cross section, MF 6 = distribution with yield 1) when no
  // reaction channel carries that particle: merged into the total-reaction channel MT 5 as a multiplicity σ_prod/σ_MT5, so that every calculation (yields, spectra, activation) sees them
  const host = channels[5],
    merged = [];
  if (host)
    for (const mt of order) {
      const sp = PROD_MT[mt],
        cp = channels[mt];
      if (!sp || !cp || !cp.prods || order.some(m => m !== mt && channels[m].primary && channels[m].y[sp] != null)) continue;
      for (const pr of cp.prods) {
        if (pr.sp !== sp) continue;
        const x = cp.E,
          y = new Float64Array(x.length);
        for (let i = 0; i < x.length; i++) {
          const s5 = chSigma(host, x[i]);
          y[i] = s5 > 0 ? (cp.s[i] * tab1Eval(pr.ytab, x[i])) / s5 : 0;
        }
        const ytab = { x, y, nbt: [x.length], int: [2] };
        (host.prods = host.prods || []).push({ ...pr, ytab });
        host.y[sp] = host.y[sp] ? { sum: [...(host.y[sp].sum || [host.y[sp]]), ytab] } : ytab;
        host.mf6 = true;
      }
      cp.merged = true;
      merged.push(mt);
    }
  if (merged.length) lumpNotes.push('MT ' + merged.join(', ') + ' (particle production) merged into MT 5 as multiplicities');
  // MF10: cross sections for the production of individual states (isomers) of the residual nuclei
  for (const [key, [a, b]] of sec) {
    const [mf, mt] = key.split('/').map(Number);
    if (mf !== 10 || !channels[mt]) continue;
    try {
      const cur = new EndfCursor(lines, a, b),
        h = cur.cont(),
        NS = h[4] | 0,
        iso = [];
      for (let k = 0; k < NS; k++) {
        const t = cur.tab1();
        iso.push({ ZA: Math.round(t.c[2]), lfs: t.c[3] | 0, tab: toMeV(t) });
      }
      channels[mt].iso = iso;
    } catch (e) {
      /* truncated */
    }
  }
  // channels without MF6 products: yields follow from the reaction definition
  for (const mt of order) {
    const ch = channels[mt];
    if (ch.mf6 && Object.keys(ch.y).length) continue;
    const em = mtEmitted(mt);
    if (!em) continue;
    for (const sp of em) ch.y[sp] = (typeof ch.y[sp] === 'number' ? ch.y[sp] : 0) + 1;
    let za = ZA + pj.Z * 1000 + pj.A; // compound nucleus = target + projectile
    for (const sp of em) za -= SPEC_ZA[sp];
    if (mt !== 4 && !(mt >= 103 && mt <= 107)) ch.res.push({ ZA: za, y: 1 });
  }
  return { ZA, Z, A, AWR, title, channels, order, lumpNotes, pj: pj.letter };
}

// multiplicity of species `sp` of channel `ch` at energy E (MeV)
function chYield(ch, sp, E) {
  const y = ch.y[sp];
  if (y == null) return 0;
  if (typeof y === 'number') return y;
  if (y.sum) {
    let s = 0;
    for (const t of y.sum) s += typeof t === 'number' ? t : tab1Eval(t, E);
    return s;
  }
  return tab1Eval(y, E);
}
function chResYield(r, E) {
  return typeof r.y === 'number' ? r.y : tab1Eval(r.y, E);
}
function chSigma(ch, E) {
  return tab1Eval(ch.tab, E);
}
// fractions of the residual r of a channel going to each isomeric state at energy E (s = channel cross section, barn): [[state, fraction]].
// LIP of the MF 6 product, else the MF 10 state cross sections of the channel (fraction of the channel), else ground state; shared by the
// activation panel (alibChannel) and the production-cross-section card, so both give identical numbers.
function chResStates(ch, r, E, s) {
  if (r.lip > 0) return [[Math.min(9, r.lip), 1]];
  const out = [];
  let f = 0;
  if (ch.iso)
    for (const q of ch.iso) {
      if (q.ZA !== r.ZA || !(q.lfs > 0)) continue;
      const fq = Math.min(1, Math.max(0, tab1Eval(q.tab, E) / s));
      f += fq;
      out.push([Math.min(9, q.lfs), fq]);
    }
  out.push([0, 1 - Math.min(f, 1)]);
  return out;
}

// ENDF-style file name of a library file for an isotope:  a-Ar040.tendl
const ELSYM = [
  'n',
  'H',
  'He',
  'Li',
  'Be',
  'B',
  'C',
  'N',
  'O',
  'F',
  'Ne',
  'Na',
  'Mg',
  'Al',
  'Si',
  'P',
  'S',
  'Cl',
  'Ar',
  'K',
  'Ca',
  'Sc',
  'Ti',
  'V',
  'Cr',
  'Mn',
  'Fe',
  'Co',
  'Ni',
  'Cu',
  'Zn',
  'Ga',
  'Ge',
  'As',
  'Se',
  'Br',
  'Kr',
  'Rb',
  'Sr',
  'Y',
  'Zr',
  'Nb',
  'Mo',
  'Tc',
  'Ru',
  'Rh',
  'Pd',
  'Ag',
  'Cd',
  'In',
  'Sn',
  'Sb',
  'Te',
  'I',
  'Xe',
  'Cs',
  'Ba',
  'La',
  'Ce',
  'Pr',
  'Nd',
  'Pm',
  'Sm',
  'Eu',
  'Gd',
  'Tb',
  'Dy',
  'Ho',
  'Er',
  'Tm',
  'Yb',
  'Lu',
  'Hf',
  'Ta',
  'W',
  'Re',
  'Os',
  'Ir',
  'Pt',
  'Au',
  'Hg',
  'Tl',
  'Pb',
  'Bi',
  'Po',
  'At',
  'Rn',
  'Fr',
  'Ra',
  'Ac',
  'Th',
  'Pa',
  'U',
  'Np',
  'Pu',
  'Am',
  'Cm',
  'Bk',
  'Cf',
  'Es',
  'Fm'
];
