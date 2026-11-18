// @ts-check
/** @typedef {{ nuc: string, id: string, rec: string[], ct: Object<number, string> }} EnsSet   one decay dataset: header and the records (continuation records in ct, by index of the record they belong to) */
/** @typedef {{ E: { v: number, unk: boolean }, dE: number, feed: number, feedP?: number, g: any[], row: string, vac?: Vacancies, sc?: number }} EnsLevel */
/* ====== ENSDF decay data (alternative decay library) ======
   ENSDF mass-chain files (ensdf.NNN, 80-column records, one file per mass number A) are read as a third decay-data source, next to the internal database and the ENDF-6 decay files.
   A dataset "60NI  60CO B- DECAY (1925.28 D)" holds the decay of ONE parent state (P record) into the levels of the daughter (L records), with the feeding (B, E, A records) and the gamma rays (G records).
   Intensities are relative: photons per 100 decays of the parent = RI · NR · BR (N record: NR gamma normalisation, BR branching of the parent into this mode, NB feeding normalisation; the PN record, when filled, gives the products).
   From the datasets the following is built for each parent state (id = 10·ZA + isomer): half-life, branches into the states of the daughters and the discrete gamma lines per decay.
   Isomers: the levels flagged M (columns 78–79) of the Adopted Levels (and of the decay datasets) are numbered by energy (M1, M2, …). The beta/EC/alpha feeding of an isomeric daughter level (plus the gamma feeding from above, cascade analysis)
   goes to the isomer; the gamma rays emitted by the isomer itself are NOT part of the parent (they belong to the decay of the isomer, as in the ENDF/B decay files).
   Annihilation photons: 2 × β+ intensity (E records), as a 511 keV line. X-rays and conversion electrons are not tabulated in ENSDF files. */
const ENS = {
  idx: new Map(),
  done: new Set(),
  con: new Map(),
  iso: new Map(),
  stab: new Set(),
  map: new Map(),
  miss: new Set(),
  nread: 0,
  readA: [],
  nuc: 0,
  name: '',
  busy: false,
  prog: '',
  msg: '',
  rest: 0
};
const ENS_TMIN = 1e-3; // isomers (levels flagged M) with T½ ≥ 1 ms are numbered m1, m2, … by energy; shorter-lived ones are part of the cascade
const ENS_TU = { Y: 3.15576e7, D: 86400, H: 3600, M: 60, S: 1, MS: 1e-3, US: 1e-6, NS: 1e-9, PS: 1e-12, FS: 1e-15, AS: 1e-18 };
function ensNum(s) {
  s = String(s == null ? '' : s).trim();
  if (!s || !/^[+-]?(\d+\.?\d*|\.\d+)(E[+-]?\d+)?$/i.test(s)) return NaN;
  return parseFloat(s);
}
function ensEn(s) {
  // level energy string → {v, unk}; "271.21", "0+Y", "202.1+X", "X"
  s = String(s || '').trim();
  const m = /^([+-]?(?:\d+\.?\d*|\.\d+)(?:E[+-]?\d+)?)?\s*(\+?\s*[A-Z])?$/i.exec(s);
  if (!m) return { v: NaN, unk: true };
  return { v: m[1] ? parseFloat(m[1]) : 0, unk: !!m[2] };
}
function ensT(s) {
  // half-life text "1925.28 D" → seconds | null (stable or unknown)
  s = String(s || '')
    .trim()
    .toUpperCase();
  if (!s || s.indexOf('STABLE') === 0) return null;
  const m = /^([0-9.]+(?:E[+-]?\d+)?)\s*([A-Z]+)/.exec(s);
  if (!m) return null;
  const v = parseFloat(m[1]);
  if (!(v > 0)) return null;
  if (ENS_TU[m[2]]) return v * ENS_TU[m[2]];
  const w = { EV: 1, KEV: 1e3, MEV: 1e6 }[m[2]];
  if (w) return (6.582119569e-16 * Math.LN2) / (v * w); // level width Γ → T½ = ħ ln2 / Γ
  return null;
}
function ensMode(tok) {
  // ENSDF decay-mode token → {code (DK_LAB key | token), lab}
  const t = String(tok).toUpperCase(),
    M = { 'B-': 'B-', EC: 'EC', 'B+': 'EC', 'EC+B+': 'EC', A: 'A', IT: 'IT', P: 'p', N: 'n', '2B-': '2B-', '2EC': '2B+', '2N': '2n', '2P': '2p', SF: 'SF' };
  if (M[t]) return { code: M[t], lab: M[t] === 'IT' ? 'IT' : DK_LAB[M[t]] || t };
  const m = /^(B-|B\+|EC)((?:\d?[NPADT])+)$/.exec(t); // β-delayed particles: B-N, B-2N, B+P, B+3P, ECP, EC2P, ECA, B-A …
  if (m) {
    const pl = { N: 'n', P: 'p', A: 'α', D: 'd', T: 't' },
      parts = m[2].match(/\d?[NPADT]/g),
      lab = (m[1] === 'B-' ? 'β−' : 'β+') + parts.map(q => (q.length > 1 ? q[0] : '') + pl[q.slice(-1)]).join(''),
      code = (m[1] === 'B-' ? 'B-' : 'B+') + m[2].toLowerCase();
    return { code, lab };
  }
  return { code: t, lab: DK_LAB[t] || t };
}
function ensGZ(nuc) {
  // " 60CO" → ZA | 0
  const m = /^\s*(\d{1,3})\s*([A-Za-z]{1,2})\s*$/.exec(nuc);
  if (!m) return 0;
  const Z = ELSYM.indexOf(m[2][0].toUpperCase() + m[2].slice(1).toLowerCase());
  return Z > 0 ? Z * 1000 + parseInt(m[1], 10) : 0;
}
function ensIsoFind(za, E, tolFix) {
  // index (1…) of the isomer of ZA with energy E (keV), or 0 if none matches
  const a = ENS.iso.get(za);
  if (!a) return 0;
  let best = 0,
    bd = Infinity;
  a.forEach((x, k) => {
    const d = Math.abs(x.v - E.v);
    if (d < bd && x.unk === E.unk) {
      bd = d;
      best = k + 1;
    }
  });
  return bd <= Math.max(tolFix || 0.6, 0.0006 * E.v) ? best : 0;
}
function ensIsoAdd(za, E) {
  // adds an isomer energy (keeps the list sorted)
  if (!(E.v >= 0)) return;
  let a = ENS.iso.get(za);
  if (!a) ENS.iso.set(za, (a = []));
  if (a.some(x => x.unk === E.unk && Math.abs(x.v - E.v) <= Math.max(0.6, 0.0006 * E.v))) return;
  a.push(E);
  a.sort((x, y) => x.v - y.v || (x.unk ? 1 : 0) - (y.unk ? 1 : 0));
}
function ensParseText(text) {
  // reads the datasets of a text and stores isomers, stable nuclides and decay contributions; returns the number of decay datasets
  const L = text.split('\n');
  /** @type {EnsSet[]} */ const sets = [];
  /** @type {EnsSet | null | false} */ let cur = null;
  let nds = 0;
  for (let i = 0; i < L.length; i++) {
    const s = L[i].replace(/\r$/, '');
    if (!s.trim()) {
      cur = null;
      continue;
    }
    if (cur && s.length > 9 && s[5] !== ' ' && s[6] === ' ' && cur.rec.length) {
      const q = cur.rec.length - 1;
      cur.ct[q] = (cur.ct[q] || '') + ' $' + s.substr(9);
      continue;
    } // continuation record: its fields (CK, CL, KC, LC …) belong to the record before
    if (cur === null) {
      if (s.length > 12 && s.substr(5, 4) === '    ' && /^\s*\d{1,3}[A-Za-z]{1,2}\s*$/.test(s.substr(0, 5)) && s.substr(9, 30).trim()) {
        cur = { nuc: s.substr(0, 5), id: s.substr(9, 30).trim(), rec: [], ct: {} };
        sets.push(cur);
      } else cur = false;
      continue;
    }
    if (cur === false) continue;
    if (s[5] !== ' ' || (s[6] !== ' ' && !(s[6] === 'P' && s[7] === 'N'))) continue; // continuation, comment and delayed-particle records are not needed
    cur.rec.push(s);
  }
  // pass 1: isomers (levels flagged M) and stable ground states
  for (const d of sets) {
    const za = ensGZ(d.nuc);
    if (!za) continue;
    for (const s of d.rec) {
      if (s[7] !== 'L') continue;
      const E = ensEn(s.substr(9, 10));
      if (!(E.v >= 0)) continue;
      const ms = s.substr(77, 2).trim();
      if (E.v === 0 && !E.unk && /^STABLE/i.test(s.substr(39, 10).trim())) ENS.stab.add(za);
      if (/^M\d?$/.test(ms) && ensT(s.substr(39, 10)) >= ENS_TMIN) ensIsoAdd(za, E);
    }
  }
  // pass 2: decay datasets
  for (const d of sets) {
    const m = /^(\d+)\s*([A-Z]{1,2})\s+(\S+)\s+DECAY\b/i.exec(d.id);
    if (!m) continue;
    const dza = ensGZ(d.nuc),
      mode = ensMode(m[3]);
    if (!dza) continue;
    let c;
    try {
      c = ensDataset(d, dza, mode);
    } catch (e) {
      c = null;
      if (globalThis.ENS_DBG) console.log('ERR', d.id, e.message);
    }
    if (c) {
      nds++;
      let a = ENS.con.get(c.za);
      if (!a) ENS.con.set(c.za, (a = []));
      a.push(c);
    }
  }
  return nds;
}
function ensDataset(d, dza, mode) {
  const f = (s, a, b) => s.substring(a - 1, b).trim();
  let P = null,
    N = null,
    PN = null;
  /** @type {EnsLevel[]} */ const lv = [];
  for (let ri = 0; ri < d.rec.length; ri++) {
    const s = d.rec[ri],
      ct = d.ct[ri] || '',
      t = s[6] === 'P' ? 'PN' : s[7],
      cf = k => {
        const m = new RegExp('(?:^|[ $])' + k + '\\s*=\\s*([0-9.eE+-]+)').exec(ct);
        return m ? ensNum(m[1]) : NaN;
      };
    if (t === 'PN') {
      if (!PN) PN = { nrbr: ensNum(f(s, 10, 19)), ntbr: ensNum(f(s, 22, 29)), nbbr: ensNum(f(s, 42, 49)) };
    } else if (t === 'P' && !P) P = { za: ensGZ(s.substr(0, 5)), E: ensEn(f(s, 10, 19)), T: ensT(f(s, 40, 49)) };
    else if (t === 'N' && !N) N = { nr: ensNum(f(s, 10, 19)), nt: ensNum(f(s, 22, 29)), br: ensNum(f(s, 32, 39)), nb: ensNum(f(s, 42, 49)) };
    else if (t === 'L') lv.push({ E: ensEn(f(s, 10, 19)), dE: ensNum(f(s, 20, 21)), feed: 0, g: [], row: s });
    else if (lv.length) {
      const l = lv[lv.length - 1];
      if (t === 'B') l.feed += ensNum(f(s, 22, 29)) || 0;
      else if (t === 'E') {
        const ib = ensNum(f(s, 22, 29)) || 0,
          ie = ensNum(f(s, 32, 39)) || 0,
          ti = ensNum(f(s, 65, 74));
        l.feedP = (l.feedP || 0) + ib;
        l.feed += ib + ie > 0 ? ib + ie : ti > 0 ? ti : 0;
        if (ie > 0) {
          const v = xrEcSplit(cf('CK'), cf('CL'));
          l.vac = l.vac || { K: 0, L1: 0, L2: 0, L3: 0 };
          for (const q in v) l.vac[q] += ie * v[q];
        }
      } // electron capture: vacancies per 100 decays of the branch
      else if (t === 'A') l.feed += ensNum(f(s, 22, 29)) || 0;
      else if (t === 'G') {
        const E = ensNum(f(s, 10, 19));
        if (!(E > 0)) continue;
        const dri = f(s, 30, 31).toUpperCase(),
          ri = ensNum(f(s, 22, 29)),
          cc = ensNum(f(s, 56, 62)),
          ti = ensNum(f(s, 65, 74));
        l.g.push({
          E,
          kc: cf('KC'),
          lc: cf('LC'),
          l1c: cf('L1C'),
          l2c: cf('L2C'),
          l3c: cf('L3C'),
          ri: dri === 'LT' || dri === 'LE' ? NaN : ri,
          cc: cc > 0 ? cc : 0,
          hcc: cc > 0,
          e0: /E0/i.test(f(s, 32, 41)),
          ti: ti > 0 ? ti : NaN
        });
      }
    }
  }
  if (!P || !P.za || !lv.length) return null;
  const def = (x, y) => (x === x ? x : y),
    br = def(N && N.br, 1),
    nr = def(N && N.nr, 1),
    nt = def(N && N.nt, 1),
    nb = def(N && N.nb, 1);
  const kR = PN && PN.nrbr > 0 ? PN.nrbr : nr * br,
    kT = PN && PN.ntbr > 0 ? PN.ntbr : nt * br,
    kB = PN && PN.nbbr > 0 ? PN.nbbr : nb * br; // photons: RI·kR; transitions: TI·kT; feeding: IB·kB (per 100 parent decays)
  // no NR and no PN: the RI are relative (strongest line = 100): no absolute photon intensity can be derived, the lines of such a dataset are not used
  let maxRI = 0;
  for (const l of lv) for (const g of l.g) if (g.ri > maxRI) maxRI = g.ri;
  const noPh = !(N && N.nr === N.nr) && !(PN && PN.nrbr > 0) && maxRI >= 99.9;
  const n = lv.length;
  lv.sort((a, b) => a.E.v - b.E.v);
  // levels of the cascade: a transition goes to the level closest to E_i − Eγ
  const tr = lv.map((l, i) =>
    l.g.map(g => {
      let j = -1,
        bd = Infinity;
      const Ef = l.E.v - g.E;
      for (let k = 0; k < i; k++) {
        const dd = Math.abs(lv[k].E.v - Ef);
        if (dd < bd) {
          bd = dd;
          j = k;
        }
      }
      const T = g.ti === g.ti ? g.ti * kT : g.ri === g.ri ? g.ri * (1 + g.cc) * kR : 0;
      return { g, to: j >= 0 && bd <= Math.max(2, 0.004 * g.E) ? j : -1, T };
    })
  );
  let kp = -1; // IT decay: the parent level is the level of the daughter dataset nearest to the energy of the P record
  if (P.za === dza) {
    let bd = Infinity;
    lv.forEach((l, k) => {
      const dd = Math.abs(l.E.v - P.E.v);
      if (dd < bd) {
        bd = dd;
        kp = k;
      }
    });
    if (bd > Math.max(3, 0.01 * P.E.v)) return null;
  }
  const isParentLevel = l => kp >= 0 && lv[kp] === l;
  const isIso = lv.map(l => ensIsoFind(dza, l.E) > 0);
  const src = lv.map(l => l.feed * kB); // direct feeding (per 100 decays)
  if (kp >= 0) src[kp] += 100 * br; // IT: the parent level is the source
  const pop = new Float64Array(n),
    popM = new Float64Array(n); // pop: all flow; popM: flow that does not pass through an isomer (the isomer's own decay is part of the isomer)
  for (let k = n - 1; k >= 0; k--) {
    pop[k] += src[k];
    popM[k] += src[k];
    const stop = isIso[k] && !isParentLevel(lv[k]); // an isomeric daughter level collects its feeding and does not de-excite here
    let Tout = 0;
    for (const t of tr[k]) Tout += t.T;
    for (const t of tr[k])
      if (t.to >= 0 && t.T > 0) {
        const w = Tout > 0 ? t.T / Tout : 0;
        pop[t.to] += pop[k] * w;
        if (!stop) popM[t.to] += popM[k] * w;
      }
    lv[k].sc = stop ? 0 : pop[k] > 0 ? popM[k] / pop[k] : 1;
  }
  const gam = [],
    Zd = Math.floor(dza / 1000),
    vac = { K: 0, L1: 0, L2: 0, L3: 0 };
  let ann = 0,
    fIso = [],
    fsum = 0;
  lv.forEach((l, k) => {
    if (l.vac) for (const q in vac) vac[q] += (l.vac[q] * kB) / 100;
    for (const t of tr[k])
      if (t.T > 0 && l.sc > 0 && t.g.cc > 0) {
        const s = xrIcSplit(Zd, t.g.E, t.g.kc, t.g.lc, t.g.cc, t.g.l1c, t.g.l2c, t.g.l3c);
        if (s) for (const q in vac) vac[q] += (((t.T * l.sc) / 100) * s[q]) / (1 + t.g.cc);
      } // internal conversion: conversion electrons of the shell = transition intensity × α_shell / (1 + α)
    for (const t of tr[k]) {
      const I = noPh || t.g.e0 ? NaN : t.g.ri === t.g.ri ? t.g.ri * kR : t.g.ti === t.g.ti && t.g.hcc ? (t.g.ti * kT) / (1 + t.g.cc) : NaN;
      if (I === I && I > 0 && l.sc > 0) gam.push([t.g.E, (I * l.sc) / 100]);
    } // RI missing: photons from the total transition intensity TI/(1+α)
    if (l.feedP) ann += (2 * l.feedP * kB) / 100;
    if (isIso[k] && !isParentLevel(l) && popM[k] > 0) {
      fIso.push({ E: l.E, f: popM[k] / 100 });
      fsum += popM[k] / 100;
    }
  });
  const f0 = Math.max(0, br - fsum);
  let ib = 0,
    ft = 0;
  for (const l of lv) {
    ib += l.feedP || 0;
    ft += l.feed;
  } // β+ share of the EC/β+ feeding: Σ IB (β+) / Σ (IB + IE) over the E records
  const pp = mode.code === 'EC' && ft > 0 ? Math.min(1, ib / ft) : null;
  const xr = xrEmit(Zd, vac);
  return { za: P.za, pE: P.E, T: P.T, lab: mode.lab, code: mode.code, br, bg: !!(N && N.br === N.br), dza, f0, fIso, gam, xr, ann, pp };
}
function ensEntry(id) {
  // builds the entry of one parent state from the stored contributions | null
  const za = Math.floor(id / 10),
    iso = id % 10,
    has = ENS.iso.has(za);
  const cons = (ENS.con.get(za) || []).filter(c => (c.pE.v === 0 && !c.pE.unk ? 0 : ensIsoFind(za, c.pE, 1.5) || (!has && c.T >= ENS_TMIN ? 1 : -1)) === iso); // parent states shorter-lived than 1 ms (not isomers here) are left out
  if (!cons.length) return null;
  const e = { st: 'rad', T: null, br: [] };
  const brs = new Map(),
    gam = [];
  /** @type {[number, number][]} */ const xr = [];
  let ann = 0;
  // datasets without BR (the branching is not stated): they share what the datasets with a stated BR leave over (usually 1 − Σ BR)
  const st = c => c.bg && !(c.br === 1 && cons.length > 1); // BR = 1.0 written in one of several datasets of a parent is the default, not a statement
  const blank = cons.filter(c => !st(c)),
    given = cons.reduce((a, c) => a + (st(c) ? c.br : 0), 0),
    share = blank.length && given > 0 ? Math.max(0, 1 - given) / blank.length : 1;
  for (let c of cons) {
    const r = st(c) || share === 1 ? 1 : share / (c.br || 1);
    if (r !== 1)
      c = Object.assign({}, c, {
        f0: c.f0 * r,
        fIso: c.fIso.map(x => ({ E: x.E, f: x.f * r })),
        gam: c.gam.map(g => [g[0], g[1] * r]),
        xr: c.xr.map(g => [g[0], g[1] * r]),
        ann: c.ann * r
      });
    if (c.T > 0 && !(e.T > 0)) e.T = c.T;
    const add = (d, f, lab) => {
      if (!(f > 1e-12)) return;
      const k = d + '|' + lab;
      let b = brs.get(k);
      if (!b) brs.set(k, (b = { f: 0, d, lab, pf: 0, fk: 0 }));
      b.f += f;
      if (c.pp != null) {
        b.pf += f * c.pp;
        b.fk += f;
      } // fk: part of the branch whose β+ share is known
    };
    add(c.dza * 10, c.f0, c.lab);
    for (const x of c.fIso) {
      const im = ensIsoFind(c.dza, x.E, 1.5);
      add(c.dza * 10 + (im > 0 ? im : 0), x.f, c.lab);
    }
    for (const g of c.gam) gam.push([g[0], g[1]]);
    for (const g of c.xr) xr.push([g[0], g[1]]);
    ann += c.ann;
  }
  e.br = Array.from(brs.values()).map(b => {
    const o = { f: b.f, d: b.d, lab: b.lab };
    if (b.fk > 0) o.pp = b.pf / b.fk;
    return o;
  });
  let fs = e.br.reduce((a, b) => a + b.f, 0);
  if (fs > 1.0005) for (const b of e.br) b.f /= fs; // branching ratios not given in several datasets: renormalised
  if (xr.length) {
    e.xr = xrMerge(xr);
    e.xtot = e.xr.reduce((a, g) => a + g[1], 0);
  }
  if (ann > 0) gam.push([510.999, ann, 1]);
  gam.sort((a, b) => a[0] - b[0]);
  if (gam.length) {
    e.gam = gam;
    e.gtot = gam.reduce((a, g) => a + g[1], 0);
    if (ann > 0) e.ann = ann;
  }
  return e;
}
// ---- loading ----
function ensRanges(a) {
  // [22, 23, 24, 60] → "22–24, 60"
  a = Array.from(new Set(a)).sort((x, y) => x - y);
  const out = [];
  for (let i = 0; i < a.length; ) {
    let j = i;
    while (j + 1 < a.length && a[j + 1] === a[j] + 1) j++;
    out.push(j > i ? a[i] + '–' + a[j] : '' + a[i]);
    i = j + 1;
  }
  return out.join(', ');
}
function ensName(fname) {
  let m = /ensdf[._-]?(\d{1,3})(?:\D|$)/i.exec(fname);
  if (!m) m = /\.(\d{3})$/.exec(fname);
  return m ? parseInt(m[1], 10) : 0;
}
async function ensLoadFiles(files, name) {
  try {
    await ensLoadFiles0(files, name);
  } catch (e) {
    ENS.busy = false;
    ENS.prog = '';
    ENS.msg = 'Reading the ENSDF files failed: ' + ((e && e.message) || e);
    actRender();
  }
}
async function ensLoadFiles0(files, name) {
  decResetAll();
  ENS.name = name;
  ENS.busy = true;
  const rest = [];
  for (const f of files) {
    if (/\.(png|jpe?g|gif|pdf|zip|gz|html?|tendl|endf)$/i.test(f.name)) continue;
    const A = ensName(f.name);
    if (A) {
      if (!ENS.idx.has(A)) ENS.idx.set(A, f);
    } else if (f.size < 4e8) rest.push(f);
  }
  let k = 0;
  for (const f of rest) {
    // files without a mass number: read completely (checked for ENSDF datasets)
    k++;
    ENS.prog = 'Reading ENSDF… ' + k + ' of ' + rest.length;
    actRender();
    await new Promise(r => setTimeout(r, 0));
    try {
      const t = await f.text();
      if (/DECAY|ADOPTED LEVELS/.test(t) && ensParseText(t)) ENS.rest++;
    } catch (e) {}
  }
  ENS.nuc = ENS.idx.size + (ENS.rest ? ENS.con.size : 0);
  ENS.busy = false;
  ENS.prog = '';
  ENS.msg = ENS.nuc
    ? (ENS.idx.size ? ENS.idx.size + ' ENSDF mass-chain files indexed in “' + name + '”' : '') +
      (ENS.rest ? (ENS.idx.size ? ' and ' : '') + ENS.con.size + ' parent nuclides read from ' + ENS.rest + ' file(s) without a mass number' : '') +
      '. A mass chain is read when its nuclides appear in the chains.'
    : 'No ENSDF files (ensdf.NNN) found in “' + name + '”.';
  if (ENS.nuc) irrSet().src = 'ens';
  saveLocal();
  actRender();
  if (ENS.nuc) actUpdate();
}
async function ensReadId(id) {
  // reads the mass chains A−5 … A that contain decays of this parent and builds its entry
  const A = Math.floor(id / 10) % 1000;
  for (let a = A; a >= Math.max(1, A - 5); a--) {
    if (ENS.done.has(a)) continue;
    ENS.done.add(a);
    const f = ENS.idx.get(a);
    if (!f) continue;
    try {
      ensParseText(await f.text());
      ENS.nread++;
      ENS.readA.push(a);
    } catch (e) {}
  }
  const e = ensEntry(id);
  if (e) ENS.map.set(id, e);
  else ENS.miss.add(id);
}
async function ensPrefetch(ids, fb, prog) {
  // as decPrefetch for the ENSDF library
  const seen = new Set(),
    q = ids.slice();
  let t0 = Date.now();
  while (q.length) {
    const id = q.shift();
    if (seen.has(id)) continue;
    seen.add(id);
    if (!ENS.map.has(id) && !ENS.miss.has(id)) {
      const ie = dkInternal(id),
        A = Math.floor(id / 10) % 1000;
      if (ie && ie.st === 'stable') {
        /* nothing to read */
      } else if (ENS.rest || Array.from(ENS.idx.keys()).some(a => a <= A && a >= A - 5)) {
        if (prog) prog('Reading the ENSDF data of ' + dkNuc(id) + ' (' + ENS.nread + ' mass chains read so far)…');
        await ensReadId(id);
        if (Date.now() - t0 > 40) {
          await new Promise(r => setTimeout(r, 0));
          t0 = Date.now();
        }
      } else ENS.miss.add(id);
    }
    const e = dkInfo(id, 'ens', fb);
    if (e && e.br && e.T > 0) for (const b of e.br) if (b.d != null) q.push(b.d);
  }
  ENS.msg =
    ENS.idx.size.toLocaleString('en') +
    ' ENSDF mass-chain files indexed in “' +
    ENS.name +
    '”; only ' +
    ENS.nread +
    ' needed for these chains ' +
    (ENS.nread === 1 ? 'was' : 'were') +
    ' read (A = ' +
    ensRanges(ENS.readA) +
    ')' +
    (ENS.miss.size && fb ? '; ' + ENS.miss.size + ' nuclides without ENSDF decay data use the internal values' : '') +
    '.';
  const m = $('#decMsg');
  if (m) m.textContent = ENS.msg;
}
