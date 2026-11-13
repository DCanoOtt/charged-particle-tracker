/* ====== charged-particle cross-section libraries (TENDL, ENDF-6): p, d, t, ³He, α — folder, per-isotope loading, panel ======
   One folder is indexed by file name: TENDL2025 (with the sub-folders TENDL-a, TENDL-p, TENDL-d, TENDL-t, TENDL-h) or any single TENDL-x folder.
   The files with the prefix of the projectile selected (a-, p-, d-, t-, h-) are used; only the isotopes of the target are read. */
const ALIB = { ver: 0, reads: [], name: '', idx: new Map(), counts: {}, nFiles: 0, data: new Map(), status: '', busy: false, useXs: true, handle: null, need: 'none' };
const ALIB_RE = /^([apdth])-([A-Z][a-z]?)(\d{3})([mn]?)\.(?:tendl|endf|dat|txt)$/i;

function PJ() { const sp = beamSpec(); return sp.kind === 'ion' ? pjOf(sp.Z, sp.A) : null; }       // library projectile of the current beam (null: no library)
function alibIsAlpha() { return !!PJ(); }                                                          // (name kept) true when the beam has a library
function alibF() { const p = PJ(); return (p && ALIB.idx.get(p.letter)) || new Map(); }             // file index of the current projectile: 'Z-A' → entry
function alibKey(Z, A) { const p = PJ(); return (p ? p.letter : '?') + ':' + Z + '-' + A; }

// small IndexedDB store to remember the folder (Chromium)
const ALIB_DB = {
  open() { return new Promise((ok, no) => { const r = indexedDB.open('charged-particle-tracker', 1); r.onupgradeneeded = () => r.result.createObjectStore('h'); r.onsuccess = () => ok(r.result); r.onerror = () => no(r.error); }); },
  async get() { try { const db = await this.open(); return await new Promise(ok => { const r = db.transaction('h').objectStore('h').get('alpha'); r.onsuccess = () => ok(r.result); r.onerror = () => ok(null); }); } catch (e) { return null; } },
  async set(v) { try { const db = await this.open(); await new Promise(ok => { const t = db.transaction('h', 'readwrite'); t.objectStore('h').put(v, 'alpha'); t.oncomplete = ok; t.onerror = ok; }); } catch (e) { } },
};

function alibIndexEntries(entries, name) {       // entries: [{name, path, get: async () => File}]
  const idx = new Map(), counts = {}; let n = 0;
  for (const en of entries) {
    const m = ALIB_RE.exec(en.name); if (!m) continue; const Zi = ELSYM.indexOf(m[2][0].toUpperCase() + m[2].slice(1).toLowerCase()); if (Zi < 1) continue;
    const L = m[1].toLowerCase(); n++; counts[L] = (counts[L] || 0) + 1; if (m[4]) continue;     // isomer (m, n) files are not used for natural targets
    if (!idx.has(L)) idx.set(L, new Map()); const mp = idx.get(L), k = Zi + '-' + (+m[3]), old = mp.get(k);
    if (!old || (!(old.path || '').includes('TENDL-' + L) && (en.path || '').includes('TENDL-' + L))) mp.set(k, en);
  }
  ALIB.report = null; ALIB_SIG = ''; ALIB.idx = idx; ALIB.counts = counts; ALIB.nFiles = n; ALIB.name = name; ALIB.data.clear(); alibStatus();
}
function alibStatus() {
  const p = PJ(), n = ALIB.nFiles; if (!n && !ALIB.name) { ALIB.status = ''; return; }
  if (!n) { ALIB.status = 'No files named like a-Ar040.tendl, p-Ar040.tendl … were found in “' + ALIB.name + '”.'; return; }
  const per = ['a', 'p', 'd', 't', 'h'].filter(L => ALIB.counts[L]).map(L => PROJ[L].sym + ' ' + ALIB.counts[L].toLocaleString('en')).join(', ');
  if (p && !ALIB.counts[p.letter]) ALIB.status = 'The folder “' + ALIB.name + '” has no ' + p.name + ' files (' + p.letter + '-Ar040.tendl …, folder ' + p.lib + '). Found: ' + per + '. Choose the TENDL2025 folder or ' + p.lib + '.';
  else ALIB.status = n.toLocaleString('en') + ' files indexed in “' + ALIB.name + '” (names only; ' + per + ').' + (ALIB.report && p && ALIB.report.pj === p.letter ? ' ' + ALIB.report.ok.length + ' read: only the isotopes of the target.' : '');
}
async function alibDirEntries(h, depth, pre) {
  depth = depth || 0; pre = pre || ''; const out = [];
  for await (const [nm, eh] of h.entries()) {
    if (eh.kind === 'file') out.push({ name: nm, path: pre + nm, get: () => eh.getFile() });
    else if (eh.kind === 'directory' && depth < 3) { try { for (const e of await alibDirEntries(eh, depth + 1, pre + nm + '/')) out.push(e); } catch (e) { } }
  }
  return out;
}

// inside the Claude viewer the page runs in a cross-origin frame, where showDirectoryPicker is refused: the directory input is used there
const ALIB_IN_FRAME = (() => { try { return window.self !== window.top; } catch (e) { return true; } })();
async function alibChoose() {
  if (window.showDirectoryPicker && !ALIB_IN_FRAME) {
    try {
      let h; try { h = await window.showDirectoryPicker({ mode: 'read' }); } catch (e) { if (e && e.name === 'AbortError') return; throw e; }
      ALIB.status = 'Reading the folder…'; renderAlpha(); ALIB.handle = h;
      alibIndexEntries(await alibDirEntries(h), h.name); ALIB_DB.set({ handle: h, name: h.name });
      renderAlpha(); alibAfterIndex(); return;
    } catch (e) { /* picker refused (frame, policy): fall back to the folder input below */ }
  }
  $('#alibInput').click();
}
function alibFromInput(files) {
  const en = Array.from(files).map(f => ({ name: f.name, path: f.webkitRelativePath || f.name, get: async () => f })); const rel = files[0] && files[0].webkitRelativePath ? files[0].webkitRelativePath.split('/')[0] : 'selected folder';
  alibIndexEntries(en, rel); renderAlpha(); alibAfterIndex();
}
async function alibRestore(user) {            // reopen the remembered folder (permission may need a click)
  const rec = await ALIB_DB.get(); if (!rec || !rec.handle) return false;
  try {
    let pm = await rec.handle.queryPermission({ mode: 'read' });
    if (pm !== 'granted') { if (!user) { ALIB.need = rec.name; return false; } pm = await rec.handle.requestPermission({ mode: 'read' }); }
    if (pm !== 'granted') return false;
    ALIB.handle = rec.handle; alibIndexEntries(await alibDirEntries(rec.handle), rec.handle.name); ALIB.need = 'none'; return true;
  } catch (e) { return false; }
}

// isotopes of the current target: [{Z, A, ab, n (atoms/cm3), layer}]
function alibTargetIsotopes() {
  const out = [];
  S.layers.forEach((l, li) => { if (!l.def) return; let m; try { m = matOf(l.def); } catch (e) { return; } for (const e of m.elems) for (const is of e.iso) if (is.ab > 0) out.push({ Z: e.Z, A: is.A, ab: is.ab, n: e.n * is.ab, layer: li, sym: e.sym }); });
  return out;
}
async function alibLoadIso(Z, A) {
  const key = alibKey(Z, A); if (ALIB.data.has(key)) return ALIB.data.get(key);
  const en = alibF().get(Z + '-' + A); if (!en) { ALIB.data.set(key, null); return null; }
  const f = await en.get(); ALIB.reads.push(en.name); const txt = await f.text(), d = parseEndfAlpha(txt, PJ()); d.file = en.name; ALIB.data.set(key, d); return d;
}
async function alibEnsure() {                  // load every isotope the target needs; resolves a short report
  if (!alibF().size) { alibStatus(); return null; }
  ALIB.busy = true; ALIB.ver++; ALIB_SIG = alibSig(); const iso = alibTargetIsotopes(), seen = new Set(), miss = [], ok = [];
  const need = new Set(iso.map(t => alibKey(t.Z, t.A))); for (const k of Array.from(ALIB.data.keys())) if (!need.has(k)) ALIB.data.delete(k);
  renderAlpha();
  const todo = new Set(iso.map(t => alibKey(t.Z, t.A))).size; let done = 0;
  for (const t of iso) {
    const k = alibKey(t.Z, t.A); if (seen.has(k)) continue; seen.add(k); done++;
    if (!ALIB.data.has(k)) { ALIB.prog = 'Reading ' + t.sym + '-' + t.A + ' (' + done + ' of ' + todo + ')…'; renderAlpha(); await new Promise(r => setTimeout(r, 30)); }
    try { const d = await alibLoadIso(t.Z, t.A); (d ? ok : miss).push(t.sym + '-' + t.A); } catch (e) { miss.push(t.sym + '-' + t.A + ' (' + (e.message || 'read error') + ')'); }
  }
  ALIB.busy = false; ALIB.report = { ok, miss, pj: PJ().letter }; alibStatus(); renderAlpha(); setTimeout(alibSync, 0); return ALIB.report;
}
// the files follow the materials: when the isotopes of the target (or the projectile) change, the files of the new isotopes are read automatically
let ALIB_SIG = '', ALIB_TMR = null;
function alibSig() { const p = PJ(); return (p ? p.letter : '') + ':' + Array.from(new Set(alibTargetIsotopes().map(t => t.Z + '-' + t.A))).sort().join(','); }
function alibSync() {
  if (!alibF().size || ALIB.busy || alibSig() === ALIB_SIG) return;
  clearTimeout(ALIB_TMR); ALIB_TMR = setTimeout(() => { if (!ALIB.busy && alibSig() !== ALIB_SIG) alibAfterIndex(); }, 400);
}
function alibAfterIndex() { alibEnsure().then(() => { if (typeof drawChart === 'function') { if (S.tab === 'endf') renderPlots(true); else actRender(); } }); }

// cross sections handed to the transport (non-elastic, barn) for every loaded isotope of the current projectile
function alphaXsLib() {
  const p = PJ(); if (!ALIB.useXs || !p) return null;
  const o = {}; let n = 0;
  for (const [k, d] of ALIB.data) { if (!d || !d.channels[3] || k.charAt(0) !== p.letter) continue; const c = d.channels[3]; o[k.slice(2)] = { E: Float64Array.from(c.E), s: Float64Array.from(c.s) }; n++; }
  return n ? o : null;
}

function renderAlpha() {
  const g = $('#gAlpha'); if (!g) return;
  const p = PJ(); if (!p) { g.hidden = true; return; } g.hidden = false; alibStatus();
  const rep = ALIB.report && ALIB.report.pj === p.letter ? ALIB.report : null, has = alibF().size > 0;
  let info = '<p class="hint">' + esc(ALIB.busy ? (ALIB.prog || 'Reading the cross sections of the target isotopes…') : (ALIB.status || (ALIB.need !== 'none' ? 'Last folder: “' + ALIB.need + '” — click Reconnect.' : 'No library loaded.'))) + '</p>';
  if (rep && has) info += '<p class="hint">Read from the library: ' + (rep.ok.length ? esc(rep.ok.join(', ')) : 'none') + (rep.miss.length ? '<br><span class="warn">No file for: ' + esc(rep.miss.join(', ')) + ' (model cross section used for those)</span>' : '') + '</p>';
  g.innerHTML = '<h2>Nuclear reaction data (ENDF-6)</h2><p class="hint">Evaluated cross sections (TENDL, ENDF-6) for reactions induced by protons, deuterons, tritons, ³He and α particles: folders <code>TENDL-p</code>, <code>TENDL-d</code>, <code>TENDL-t</code>, <code>TENDL-h</code>, <code>TENDL-a</code>, files <code>p-Ar040.tendl</code>, <code>d-…</code>, <code>t-…</code>, <code>h-…</code>, <code>a-…</code>. Choose the TENDL2025 folder (it holds all five) or a single TENDL-x folder; the library of the projectile selected now (' + p.name + ': <code>' + p.lib + '</code>, files <code>' + p.letter + '-Ar040.tendl</code> …) is used and only the isotopes of the target are read. With the cross sections the reaction probability and the secondary-particle yields are obtained by folding the ion flux Φ(E) with them (tab “Reactions (ENDF)”).</p>' +
    '<div class="btns"><button class="btn sm" id="alibChoose">Choose library folder…</button>' + (ALIB.need !== 'none' && !ALIB.nFiles ? '<button class="btn sm" id="alibRestore">Reconnect “' + esc(ALIB.need) + '”</button>' : '') + (has ? '<button class="btn sm" id="alibReload">Re-read target isotopes</button>' : '') + '<input type="file" id="alibInput" webkitdirectory directory multiple hidden></div>' + info +
    '<label class="chk" style="margin-top:6px"><input type="checkbox" id="alibUse"' + (ALIB.useXs ? ' checked' : '') + '><span>Use the non-elastic cross section (MT 3) for the nuclear removal of the ' + p.plural + ' in the transport</span></label>';
}
document.addEventListener('click', ev => {
  const id = ev.target.id;
  if (id === 'alibChoose') alibChoose();
  else if (id === 'alibRestore') alibRestore(true).then(ok => { renderAlpha(); if (ok) alibAfterIndex(); });
  else if (id === 'alibReload') alibAfterIndex();
});
document.addEventListener('change', ev => {
  if (ev.target.id === 'alibInput') alibFromInput(ev.target.files);
  else if (ev.target.id === 'alibUse') { ALIB.useXs = ev.target.checked; }
});

/* ====== folding of the ion flux with the cross sections ======
   R.Lh[layer*nLb + k] = path length summed over all primaries (cm) of the primary ions in energy bin k (width R.LbW MeV).
   Reactions per primary in a layer = Σ_iso n_iso · 1e-24 · Σ_k σ(E_k) · L_k ;  yield of a species = same with the multiplicity. */
function alibFold(R) {
  if (!R || !R.Lh) return null;
  const iso = alibTargetIsotopes(), nLb = R.nLb, w = R.LbW, L = R.L, Ec = new Float64Array(nLb); for (let k = 0; k < nLb; k++) Ec[k] = (k + 0.5) * w;
  const rows = [];                                    // one entry per layer × isotope
  for (const t of iso) {
    const d = ALIB.data.get(alibKey(t.Z, t.A)); if (!d || t.layer >= L) continue;
    const Lh = Float64Array.from(R.Lh.subarray(t.layer * nLb, (t.layer + 1) * nLb), v => v / R.N);   // per primary
    rows.push({ t, d, Lh });
  }
  return { rows, Ec, w, nLb };
}
// per-channel result for one isotope row: { P, Y:{sp:count}, res: Map(ZA->count), dPdE: Float64Array }
function alibChannel(F, row, mt) {
  const ch = row.d.channels[mt]; if (!ch) return null;
  const { Ec, nLb, w } = F, n = row.t.n * 1e-24, Y = {}, res = new Map(), iso = new Map(), dP = new Float64Array(nLb); let P = 0;
  for (const [sp] of SPECIES) Y[sp] = 0;
  for (let k = 0; k < nLb; k++) {
    const Lk = row.Lh[k]; if (!(Lk > 0)) continue; const E = Ec[k], s = chSigma(ch, E); if (!(s > 0)) continue;
    const p = n * s * Lk; P += p; dP[k] = p / w;
    for (const [sp] of SPECIES) { const m = chYield(ch, sp, E); if (m) Y[sp] += p * m; }
    for (const r of ch.res) {
      const m = chResYield(r, E); if (!m) continue; res.set(r.ZA, (res.get(r.ZA) || 0) + p * m);
      // isomeric state of the residual nucleus: LIP of the MF 6 product, else the MF 10 state cross sections of this channel (fraction of the channel), else ground state
      const add = (st, v) => { const id = r.ZA * 10 + st; iso.set(id, (iso.get(id) || 0) + v); };
      if (r.lip > 0) { add(Math.min(9, r.lip), p * m); continue; }
      const is = ch.iso ? ch.iso.filter(q => q.ZA === r.ZA && q.lfs > 0) : []; let f = 0;
      for (const q of is) { const fq = Math.min(1, Math.max(0, tab1Eval(q.tab, E) / s)); f += fq; add(Math.min(9, q.lfs), p * m * fq); }
      if (f > 1) f = 1; add(0, p * m * (1 - f));
    }
  }
  return { P, Y, res, iso, dPdE: dP };
}
// all primary channels summed (the non-overlapping set whose sum is the non-elastic cross section)
function alibTotal(F, row) {
  const tot = { P: 0, Y: {}, res: new Map(), iso: new Map(), dPdE: new Float64Array(F.nLb) }; for (const [sp] of SPECIES) tot.Y[sp] = 0;
  for (const mt of row.d.order) { const ch = row.d.channels[mt]; if (!ch.primary) continue; const r = alibChannel(F, row, mt); if (!r) continue; tot.P += r.P; for (const [sp] of SPECIES) tot.Y[sp] += r.Y[sp]; for (const [za, v] of r.res) tot.res.set(za, (tot.res.get(za) || 0) + v); for (const [id, v] of r.iso) tot.iso.set(id, (tot.iso.get(id) || 0) + v); for (let k = 0; k < F.nLb; k++) tot.dPdE[k] += r.dPdE[k]; }
  return tot;
}
// channels available in at least one isotope of the target: [{mt, name}] (lumped sums included, primary channels flagged)
function alibChannelList() {
  const seen = new Map();
  for (const t of alibTargetIsotopes()) { const d = ALIB.data.get(alibKey(t.Z, t.A)); if (!d) continue; for (const mt of d.order) if (mt !== 2 && !seen.has(mt)) seen.set(mt, d.channels[mt]); }
  return Array.from(seen.keys()).sort((a, b) => a - b).map(mt => ({ mt, name: seen.get(mt).name }));
}
