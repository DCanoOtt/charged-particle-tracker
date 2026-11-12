/* ====== decay data and decay chains (Bateman equations) ======
   Nuclide id = 10·ZA + isomer number (ground state: ZA·10), ZA = 1000 Z + A.
   Decay-data sources: (1) internal data = ICRP Publication 107 decay data with isomers (Python package radioactivedecay 0.6.1, NUBASE2020 half-lives; 1512 nuclides) and, for the nuclides it does not contain, the mendeleev 1.3.0 table (half-lives, decay modes and intensities; ground states);
   (2) ENDF-6 decay sub-library files (MF 8 / MT 457) loaded by the user: half-life, decay modes RTYP with daughters and branching ratios BR, isomers. */
const DK_LN2 = Math.LN2;
const DK_SHIFT = { 'B-': [1, 0], 'EC': [-1, 0], 'A': [-2, -4], 'p': [-1, -1], 'n': [0, -1], '2n': [0, -2], '3n': [0, -3], '4n': [0, -4], '2p': [-2, -2], '3p': [-3, -3], 'B-n': [1, -1], 'B+p': [-2, -1], 'B-A': [-1, -4], 'B+A': [-3, -4], 'B-d': [0, -2], 'B-t': [0, -3], 'B-p': [0, -1], 'B+pA': [-4, -5], '2B-': [2, 0], '2B+': [-2, 0] };
const DK_LAB = { 'B-': 'β−', 'EC': 'EC/β+', 'A': 'α', 'p': 'p', 'n': 'n', '2n': '2n', '3n': '3n', '4n': '4n', '2p': '2p', '3p': '3p', 'B-n': 'β−n', 'B+p': 'β+p', 'B-A': 'β−α', 'B+A': 'β+α', 'B-d': 'β−d', 'B-t': 'β−t', 'B-p': 'β−p', 'B+pA': 'β+pα', '2B-': '2β−', '2B+': '2β+', 'SF': 'SF', 'B-SF': 'β−SF', 'B+SF': 'β+SF' };
const DK_RT = { 1: [1, 0, 'β−'], 2: [-1, 0, 'EC/β+'], 3: [0, 0, 'IT'], 4: [-2, -4, 'α'], 5: [0, -1, 'n'], 6: [0, 0, 'SF'], 7: [-1, -1, 'p'] };     // ENDF-6 RTYP digits: ΔZ, ΔA, label
const DEC = { map: new Map(), idx: new Map(), miss: new Set(), nread: 0, files: 0, nuc: 0, name: '', busy: false, prog: '', msg: '' };       // uploaded ENDF-6 decay data: parsed nuclides (map), files by nuclide (idx, read when needed)
const DK_INT = new Map();

function dkNuc(id) { const za = Math.floor(id / 10), m = id % 10; return actNuc(za) + (m ? 'm' + (m > 1 ? m : '') : ''); }
function dkLab(br) { return br.map(b => (b.f >= 0.9995 ? '' : fmt(100 * b.f, 3) + ' % ') + b.lab).join(', '); }

function dkInternal(id) {                                                 // {st:'stable'|'rad', T (s) | null, br:[{f, d (id | null), lab}]} | null
  if (DK_INT.has(id)) return DK_INT.get(id); let r = null;
  if (typeof DK_ICRP !== 'undefined' && DK_ICRP[id]) {                    // ICRP-107 / NUBASE2020 data, with isomers
    const e = DK_ICRP[id], br = e[1].split(';').map(q => { const a = q.split(':'); return { f: +a[1], d: +a[0] > 0 ? +a[0] : null, lab: DK_LAB[a[2]] || a[2] }; });
    r = { st: 'rad', T: e[0], br, db: 'icrp' }; DK_INT.set(id, r); return r;
  }
  if (typeof DK_ICRP_STABLE !== 'undefined' && DK_ICRP_STABLE.has(id)) { r = { st: 'stable', T: null, br: [], db: 'icrp' }; DK_INT.set(id, r); return r; }
  if (id % 10) { DK_INT.set(id, null); return null; } const za = id / 10;
  const Z = Math.floor(za / 1000), A = za % 1000;
  if (typeof HL_STABLE !== 'undefined' && HL_STABLE.has(za)) r = { st: 'stable', T: null, br: [], db: 'mendeleev' };
  else if (HL_T[za] > 0 || HL_BR[za] || HL_MODE[za]) {
    const br = [], s = HL_BR[za];
    if (s) for (const part of s.split(';')) {
      const i = part.indexOf(':'), code = part.slice(0, i), f = +part.slice(i + 1) / 100, sh = DK_SHIFT[code], lab = DK_LAB[code] || code;
      if (!sh) { br.push({ f, d: null, lab }); continue; } const z = Z + sh[0], a = A + sh[1];
      br.push({ f, d: z >= 1 && a >= z ? (z * 1000 + a) * 10 : null, lab });
    }
    r = { st: 'rad', T: HL_T[za] > 0 ? HL_T[za] : null, br, db: 'mendeleev' };
  }
  DK_INT.set(id, r); return r;
}
function dkInfo(id, src, fallback) {                                       // decay data of one nuclide from the selected source
  if (src === 'file') { const e = DEC.map.get(id); if (e) return Object.assign({ src: 'file' }, e); if (!fallback) return null; }
  const e = dkInternal(id); return e ? Object.assign({ src: 'int' }, e) : null;
}

// ---- ENDF-6 decay data: MF 8 / MT 457 sections of one or several materials in a text ----
function dkRtyp(r, zaParent, rfs) {                                       // → {d (id | null), lab}
  const ip = Math.floor(r + 1e-9), fr = Math.round((r - ip) * 1000), steps = [ip];
  if (fr > 0) { const s = String(fr).padStart(3, '0').replace(/0+$/, ''); for (const ch of s) steps.push(+ch); }
  if (steps.includes(6) || ip === 10 || !steps.every(k => DK_RT[k])) return { d: null, lab: steps.includes(6) ? 'SF' : '?' };
  let dz = 0, da = 0; const lab = []; for (const k of steps) { dz += DK_RT[k][0]; da += DK_RT[k][1]; lab.push(DK_RT[k][2]); }
  const z = Math.floor(zaParent / 1000) + dz, a = zaParent % 1000 + da; if (z < 1 || a < z) return { d: null, lab: lab.join('+') };
  return { d: (z * 1000 + a) * 10 + Math.max(0, Math.min(9, Math.round(rfs))), lab: lab.join('+') };
}
function gamUnits(out) {                                                   // applies the per cent → per decay conversion to the entries not yet converted
  if (!out.gpct) return;
  for (const e of out.values()) if (e.graw) { e.graw = false; for (const g of e.gam) g[1] /= 100; e.gtot /= 100; }
}
function parseEndfDecay(text, out) {                                       // adds the nuclides found to the Map out; returns their number
  if (text.indexOf(' 8457') < 0) return 0; const L = text.split('\n'); let n = 0;
  for (let i = 0; i < L.length; i++) {
    if (L[i].length < 75 || L[i].substr(70, 5) !== ' 8457') continue; let j = i; while (j < L.length && L[j].length >= 75 && L[j].substr(70, 5) === ' 8457') j++;
    try {
      const h = endfFields(L[i]), ZA = Math.round(h[0]), liso = Math.round(h[3]), nst = Math.round(h[4]); let p = i + 1;
      const r1 = endfFields(L[p]), T = r1[0], npl = Math.round(r1[4]); p += 1 + Math.ceil(npl / 6);
      const entry = { st: nst === 1 ? 'stable' : 'rad', T: nst !== 1 && T > 0 ? T : null, br: [] };
      if (nst !== 1 && p < j) {
        const r2 = endfFields(L[p]), nv = Math.round(r2[4]), v = []; for (let k = 1; k <= Math.ceil(nv / 6) && p + k < j; k++) for (const x of endfFields(L[p + k])) v.push(x);
        for (let k = 0; k + 5 < nv; k += 6) { const f = v[k + 4]; if (!(f > 0)) continue; const q = dkRtyp(v[k], ZA, v[k + 1]); entry.br.push({ f, d: q.d, lab: q.lab }); }
        // discrete gamma lines (first spectrum record, STYP = 0): energies in eV, absolute intensity per decay = FD · RI
        const pn = p + 1 + Math.ceil(nv / 6);
        if (Math.round(h[5]) > 0 && pn + 1 < j) {
          const sp = endfFields(L[pn]);
          if (Math.round(sp[1]) === 0) {
            const ner = Math.round(sp[5]), fd = endfFields(L[pn + 1])[0] || 1, gam = []; let p2 = pn + 2;
            for (let k = 0; k < ner && p2 + 1 < j; k++) { const hd = endfFields(L[p2]), nt = Math.round(hd[4]), ri = endfFields(L[p2 + 1])[2]; if (hd[0] > 0 && ri > 0) gam.push([hd[0] / 1000, ri * fd]); p2 += 1 + Math.ceil(nt / 6); }
            if (gam.length) { gam.sort((a, b) => a[0] - b[0]); entry.gam = gam; entry.gtot = gam.reduce((a, g) => a + g[1], 0); }
          }
        }
      }
      out.set(ZA * 10 + Math.max(0, Math.min(9, liso)), entry); n++;
      // RI is given in “relative units”: the ENDF/B, JEFF… decay sublibraries use per cent (per 100 decays), other files use per decay.
      // A line above 1.5 can only be per cent (a per-decay intensity of one line cannot exceed ~1): once seen, all lines of the library are divided by 100.
      if (entry.gam) { entry.graw = true; if (entry.gam.some(g => g[1] > 1.5)) out.gpct = true; }
    } catch (e) { /* damaged section: skipped */ }
    i = j - 1;
  }
  gamUnits(out);
  return n;
}
// The files are only indexed by name when chosen; a file is read when its nuclide is needed (the products and, step by step, their daughters).
// File names like dec-027_Co_060.endf, dec-052_Te_127m1.endf, Co-60.endf, Co60m.dat give the nuclide; files without such a name (one big file) are read at once.
function decName(fname) {
  const stem = fname.replace(/\.[A-Za-z0-9]+$/, ''), m = /([A-Za-z]{1,2})[-_ ]?(\d{1,3})(?:[-_ ]?(m\d*|n))?$/i.exec(stem); if (!m) return null;
  const Z = ELSYM.indexOf(m[1][0].toUpperCase() + m[1].slice(1).toLowerCase()), A = +m[2]; if (Z < 1 || A < Z) return null;
  const iso = m[3] ? (/^m\d+$/i.test(m[3]) ? +m[3].slice(1) : /^n$/i.test(m[3]) ? 2 : 1) : 0;
  return (Z * 1000 + A) * 10 + Math.min(9, iso);
}
async function dkLoadFiles(files, name) {                                  // files: File[] (a folder, several files or one big file)
  DEC.busy = true; DEC.map = new Map(); DEC.idx = new Map(); DEC.miss = new Set(); DEC.nread = 0; DEC.files = 0; DEC.nuc = 0; DEC.name = name; DEC.msg = '';
  const rest = []; for (const f of files) { if (/\.(png|jpe?g|gif|pdf|zip|gz|html?|tendl)$/i.test(f.name)) continue; const id = decName(f.name); if (id) { if (!DEC.idx.has(id)) DEC.idx.set(id, f); } else if (f.size < 2e8) rest.push(f); }
  let k = 0, t0 = Date.now();
  for (const f of rest) {                                                  // files without a nuclide name: read completely
    k++; if (Date.now() - t0 > 80) { DEC.prog = 'Reading decay data… ' + k + ' of ' + rest.length; actRender(); await new Promise(r => setTimeout(r, 0)); t0 = Date.now(); }
    try { if (parseEndfDecay(await f.text(), DEC.map)) DEC.files++; } catch (e) { }
  }
  DEC.nuc = DEC.idx.size + DEC.map.size; DEC.busy = false; DEC.prog = '';
  DEC.msg = DEC.nuc ? DEC.idx.size.toLocaleString('en') + ' decay files indexed by name in “' + name + '”' + (DEC.map.size ? ' and ' + DEC.map.size.toLocaleString('en') + ' nuclides read from files without a nuclide name' : '') + ' (ENDF-6, MF 8 / MT 457). A file is read only when its nuclide appears in the chains.' : 'No ENDF-6 decay data (MF 8, MT 457) found in “' + name + '”.';
  if (DEC.nuc) irrSet().src = 'file'; saveLocal(); actRender();
}
async function decReadId(id) {
  const f = DEC.idx.get(id); if (!f) { DEC.miss.add(id); return; }
  try { parseEndfDecay(await f.text(), DEC.map); DEC.nread++; if (!DEC.map.has(id)) DEC.miss.add(id); } catch (e) { DEC.miss.add(id); }
}
async function decPrefetch(ids, fb, prog) {                                // reads the files of the products and, step by step, of their daughters (not of stable nuclides)
  const seen = new Set(), q = ids.slice(); let t0 = Date.now();
  while (q.length) {
    const id = q.shift(); if (seen.has(id)) continue; seen.add(id);
    if (!DEC.map.has(id) && !DEC.miss.has(id)) {
      const ie = dkInternal(id);
      if (ie && ie.st === 'stable') { /* nothing to read */ }
      else if (DEC.idx.has(id)) { if (prog) prog('Reading the decay data of ' + dkNuc(id) + ' (' + DEC.nread + ' files read so far)…'); await decReadId(id); if (Date.now() - t0 > 40) { await new Promise(r => setTimeout(r, 0)); t0 = Date.now(); } }
      else DEC.miss.add(id);
    }
    const e = dkInfo(id, 'file', fb); if (e && e.br && e.T > 0) for (const b of e.br) if (b.d != null) q.push(b.d);
  }
  gamUnits(DEC.map);
  DEC.msg = DEC.idx.size.toLocaleString('en') + ' decay files indexed in “' + DEC.name + '”; only the ' + DEC.nread + ' needed for these chains ' + (DEC.nread === 1 ? 'was' : 'were') + ' read' + (DEC.miss.size && fb ? '; ' + DEC.miss.size + ' nuclides without a file use the internal values' : '') + '.';
  const m = $('#decMsg'); if (m) m.textContent = DEC.msg;
}

// ---- Bateman solution: constant production during the irradiation, then free decay ----
// The system dN/dt = P + M N is linear and has a constant source, which is written as a state x0 = 1: x' = A x with A lower triangular
// (the nuclides are ordered parents first). Its exact propagator exp(A h) is computed once for the smallest step h = T / 2^K of each phase and squared
// K times (exp(A·2^j h)); every output time n·h is reached by applying the factors of the binary digits of the increment. Nuclides with λh ≥ 1e5 are
// in secular equilibrium with their parents at every output time: their activity equals their feeding, which is eliminated analytically
// (this keeps the matrix well conditioned from T½ = 1e-22 s to 1e20 y). Only extreme decay constants (λh ≥ 1e5) are eliminated; everything else is solved exactly.
const DK_K = 22;
function dkGrid(K, T, marks) {                                              // integer multiples n of the smallest step h = T/2^K: geometric (6 per octave, resolves short half-lives) + uniform (48) + marks (times in s, e.g. multiples of the half-lives)
  const N = Math.pow(2, K), s = new Set([N]); if (marks) for (const t of marks) { const n = Math.round(t / T * N); if (n >= 1 && n <= N) s.add(n); }
  for (let o = 0; o < K; o++) for (let j = 0; j < 6; j++) s.add(Math.round(Math.pow(2, o + j / 6)));
  for (let m = 1; m <= 48; m++) s.add(Math.round(N * m / 48));
  return Array.from(s).sort((a, b) => a - b);
}
const dkIx = i => (i * (i + 1)) / 2;
function dkMul(A, B, n) {                                                   // packed lower-triangular product
  const C = new Float64Array(A.length);
  for (let i = 0; i < n; i++) { const ri = dkIx(i); for (let k = 0; k <= i; k++) { const a = A[ri + k]; if (a === 0) continue; const rk = dkIx(k); for (let j = 0; j <= k; j++) C[ri + j] += a * B[rk + j]; } }
  for (let q = 0; q < C.length; q++) if (C[q] < 1e-280 && C[q] > -1e-280) C[q] = 0;      // no denormals
  return C;
}
function dkExpChain(A, n, h, K) {                                           // [exp(A h), exp(A 2h), …, exp(A 2^K h)]
  // A is lower triangular, so the diagonal of every exp(A t) is exactly exp(A_ii t): it is imposed after every squaring, which removes the
  // amplification of the rounding error of (1 − λh) by the repeated squaring; the off-diagonal elements are sums of non-negative terms
  // (compartment model) and keep a relative error of a few ulp per squaring.
  let nrm = 0; const col = new Float64Array(n), dg = new Float64Array(n); for (let i = 0; i < n; i++) { dg[i] = A[dkIx(i) + i]; for (let j = 0; j <= i; j++) col[j] += Math.abs(A[dkIx(i) + j]); }
  for (let j = 0; j < n; j++) nrm = Math.max(nrm, col[j]); nrm *= h; let s = 0; while (nrm / Math.pow(2, s) > 0.125) s++;
  const setDiag = (T, t) => { for (let i = 0; i < n; i++) T[dkIx(i) + i] = Math.exp(dg[i] * t); };
  const B = new Float64Array(A.length); for (let q = 0; q < A.length; q++) B[q] = A[q] * h / Math.pow(2, s);
  let T = new Float64Array(A.length); for (let i = 0; i < n; i++) T[dkIx(i) + i] = 1;
  for (let k = 13; k >= 1; k--) { const M = dkMul(B, T, n); for (let q = 0; q < M.length; q++) M[q] /= k; for (let i = 0; i < n; i++) M[dkIx(i) + i] += 1; T = M; }       // Taylor series of the scaled matrix (norm ≤ 1/8, truncation < 1e-20)
  let e = h / Math.pow(2, s); setDiag(T, e);
  for (let q = 0; q < s; q++) { T = dkMul(T, T, n); e *= 2; setDiag(T, e); }
  const out = [T]; for (let j = 1; j <= K; j++) { T = dkMul(T, T, n); e *= 2; setDiag(T, e); out.push(T); } return out;
}
function dkApply(P, x, n) {                                                 // y = P x
  const y = new Float64Array(n); for (let i = 0; i < n; i++) { const ri = dkIx(i); let a = 0; for (let k = 0; k <= i; k++) a += P[ri + k] * x[k]; y[i] = a; } return y;
}
async function dkSolve(Q) {
  /* Q: { prod: Map(id → atoms/s), info(id) → {st,T,br}|null, Tirr, Tc, chain (bool), maxNodes }
     → { ids, lam, fast, info, t (s since the start of the irradiation), nI, A: Float64Array[] activity (Bq) of every node at the times t, Nirr, Ncool, direct: Map } */
  const ids = [], idx = new Map(), inf = [], lam = [], queue = []; let Rmax = 0; for (const v of Q.prod.values()) Rmax = Math.max(Rmax, v);
  const add = id => { if (idx.has(id)) return idx.get(id); const i = ids.length; idx.set(id, i); ids.push(id); queue.push(i); return i; };
  for (const id of Q.prod.keys()) add(id);
  const edges = [];                                                          // [parent, child, fraction]
  for (let q = 0; q < queue.length; q++) {
    const i = queue[q], e = Q.info(ids[i]); inf[i] = e; lam[i] = e && e.st !== 'stable' && e.T > 0 ? DK_LN2 / e.T : 0;
    if (lam[i] > 0 && Q.chain && e.br) for (const b of e.br) if (b.d != null && b.d !== ids[i] && b.f > 0 && (idx.has(b.d) || ids.length < (Q.maxNodes || 1500))) edges.push([i, add(b.d), b.f]);   // branches without daughter (SF, unknown) are lost
  }
  const n0 = ids.length, indeg = new Int32Array(n0), kids = Array.from({ length: n0 }, () => []);
  for (const [p, c, f] of edges) { kids[p].push([c, f]); indeg[c]++; }
  const order = [], st = []; for (let i = 0; i < n0; i++) if (!indeg[i]) st.push(i);
  while (st.length) { const i = st.pop(); order.push(i); for (const [c] of kids[i]) if (--indeg[c] === 0) st.push(c); }
  if (order.length < n0) { const seen = new Set(order); for (let i = 0; i < n0; i++) if (!seen.has(i)) order.push(i); }       // cycles (damaged data): back edges are ignored below
  const pos = new Int32Array(n0); order.forEach((i, k) => { pos[i] = k; });
  const par = Array.from({ length: n0 }, () => []); for (const [p, c, f] of edges) if (pos[p] < pos[c]) par[c].push([p, f]);
  const K = Q.K || DK_K, hI = Q.Tirr / Math.pow(2, K), hC = Q.Tc > 0 ? Q.Tc / Math.pow(2, K) : Infinity, hm = Math.min(hI, hC);
  const fast = lam.map(l => l * hm >= 1e5), sl = new Int32Array(n0).fill(-1); let ns = 0;
  for (const i of order) if (!fast[i]) sl[i] = ++ns; const nn = ns + 1;
  const A = new Float64Array(dkIx(nn + 1)), fvec = new Array(n0), scale = Rmax > 0 ? Rmax : 1;
  for (const i of order) {
    const v = new Float64Array(nn); if (Q.prod.has(ids[i])) v[0] += Q.prod.get(ids[i]) / scale;
    for (const [p, f] of par[i]) { if (!fast[p]) v[sl[p]] += f * lam[p]; else { const w = fvec[p]; for (let k = 0; k < nn; k++) v[k] += f * w[k]; } }
    if (fast[i]) fvec[i] = v; else { const r = dkIx(sl[i]); for (let k = 0; k < nn; k++) A[r + k] = v[k]; A[r + sl[i]] -= lam[i]; }
  }
  await new Promise(r => setTimeout(r, 0));
  // each phase gets its own grid, built from its own length (irradiation T, cooling Tc) and from the half-lives of the strongest nuclei of the chain: ½, 1, 2, 3, 5 half-lives
  const lead = ids.map((id, i) => i).filter(i => lam[i] > 0).sort((a, b) => (Q.prod.get(ids[b]) || 0) - (Q.prod.get(ids[a]) || 0) || lam[a] - lam[b]).slice(0, 30), marks = [];
  for (const i of lead) for (const m of [0.5, 1, 2, 3, 5]) marks.push(m * DK_LN2 / lam[i]);
  const mI = marks.slice(), mC = marks.slice(); for (const te of Q.extra || []) { if (te > 0 && te <= Q.Tirr) mI.push(te); else if (te > Q.Tirr) mC.push(te - Q.Tirr); }      // extra output times (slice limits of the γ spectra)
  const gI = dkGrid(K, Q.Tirr, mI), gC = Q.Tc > 0 ? dkGrid(K, Q.Tc, mC) : [], grid = gI, Npow = Math.pow(2, K), t = [], nI = gI.length, nC = gC.length, X = [], bitsRun = (chain, d, x) => { for (let j = 0; d > 0; j++, d = Math.floor(d / 2)) if (d % 2) x = dkApply(chain[j], x, nn); return x; };
  let chain = dkExpChain(A, nn, hI, K), x = new Float64Array(nn), prev = 0; x[0] = 1;
  for (const g of grid) { x = bitsRun(chain, g - prev, x); prev = g; t.push(g * hI); X.push(x); }
  chain = null;
  if (nC) {
    await new Promise(r => setTimeout(r, 0)); const c2 = dkExpChain(A, nn, hC, K); let y = Float64Array.from(X[X.length - 1]); y[0] = 0; prev = 0;
    for (const g of gC) { y = bitsRun(c2, g - prev, y); prev = g; t.push(Q.Tirr + g * hC); X.push(y); }
  }
  const nT = t.length, act = ids.map(() => new Float64Array(nT)), num = ids.map(() => new Float64Array(nT));      // activity (Bq) and number of atoms of every nuclide at every output time
  for (let i = 0; i < n0; i++) {
    const a = act[i], nu = num[i];
    if (lam[i] <= 0) { const s = sl[i]; for (let q = 0; q < nT; q++) nu[q] = X[q][s] * scale; continue; }
    if (fast[i]) { const w = fvec[i]; for (let q = 0; q < nT; q++) { let s = 0; const xq = X[q]; for (let k = 0; k < nn; k++) s += w[k] * xq[k]; a[q] = s * scale; } }
    else { const s = sl[i]; for (let q = 0; q < nT; q++) a[q] = lam[i] * X[q][s] * scale; }
    for (let q = 0; q < nT; q++) nu[q] = fast[i] ? a[q] / lam[i] : X[q][sl[i]] * scale;
  }
  const atoms = q => ids.map((id, i) => lam[i] > 0 ? act[i][q] / lam[i] : (fast[i] ? 0 : X[q][sl[i]] * scale));
  return { ids, lam, fast, info: inf, t, nI, nT, act, num, order, par, Nirr: atoms(nI - 1), Ncool: nC ? atoms(nT - 1) : null, nn: ns };
}
