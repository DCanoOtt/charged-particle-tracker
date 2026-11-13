// ---- random numbers ----
function makeRng(seed) {
  let a = seed >>> 0 || 123456789, b = 0x9E3779B9, c = 0x243F6A88, d = 0xB7E15162;
  function next() { a |= 0; b |= 0; c |= 0; d |= 0; const t = (((a + b) | 0) + d) | 0; d = (d + 1) | 0; a = b ^ (b >>> 9); b = (c + (c << 3)) | 0; c = (c << 21) | (c >>> 11); c = (c + t) | 0; return (t >>> 0) / 4294967296; }
  for (let i = 0; i < 20; i++) next();
  return next;
}
class VBuf {
  constructor(nc, cap) { this.nc = nc; this.cap = cap || 1 << 16; this.a = new Float32Array(this.cap * nc); this.n = 0; }
  push4(x, y, z, e) { if (this.n >= this.cap) { this.cap *= 2; const b = new Float32Array(this.cap * this.nc); b.set(this.a); this.a = b; } const o = this.n * this.nc; this.a[o] = x; this.a[o + 1] = y; this.a[o + 2] = z; if (this.nc > 3) this.a[o + 3] = e; this.n++; }
  done() { return this.a.slice(0, this.n * this.nc); }
}
class FBuf {
  constructor(cap) { this.cap = cap || 1 << 14; this.a = new Float32Array(this.cap); this.n = 0; }
  push(v) { if (this.n >= this.cap) { this.cap *= 2; const b = new Float32Array(this.cap); b.set(this.a); this.a = b; } this.a[this.n++] = v; }
  done() { return this.a.slice(0, this.n); }
}

function estimateDepth(layers, tabs, E0) {
  let E = E0, z = 0;
  for (let l = 0; l < layers.length; l++) {
    const t = layers[l].t, tb = tabs[l];
    if (!tb) { z += t; continue; }
    const R = interpR(tb, E);
    if (R <= t) return z + R;
    E = invR(tb, R - t); z += t;
  }
  return z;
}
function interpR(tb, E) { if (E <= tb.Emin) return tb.R[0] * E / tb.Emin; let x = (Math.log(E) - tb.l0) / tb.dl; if (x > tb.n - 1.0001) x = tb.n - 1.0001; const i = x | 0, f = x - i; return tb.R[i] + (tb.R[i + 1] - tb.R[i]) * f; }
function invR(tb, R) { let lo = 0, hi = tb.n - 1; if (R <= tb.R[0]) return tb.E[0] * R / tb.R[0]; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (tb.R[m] > R) hi = m; else lo = m; } const f = (R - tb.R[lo]) / (tb.R[hi] - tb.R[lo]); return tb.E[lo] * Math.pow(tb.E[hi] / tb.E[lo], f); }

function simulate(cfg, progress) {
  const ph = Object.assign({ Tcut: 0.02, delta: true, deltaTrack: true, reactions: true, rxScale: 1, scatter: true, straggle: true, stepFrac: 0.1, Nh: 8, EcutE: 0.001, EcutI: 0.001, nucStop: true }, cfg.phys || {});
  XSLIB = cfg.xslib || null; XSLIB_SP = cfg.xslibSp || { Z: 2, A: 4 };
  const rnd = makeRng(cfg.seed || 12345);
  let gspare = null;
  const gauss = () => { if (gspare !== null) { const g = gspare; gspare = null; return g; } let u, v, s; do { u = 2 * rnd() - 1; v = 2 * rnd() - 1; s = u * u + v * v; } while (s >= 1 || s === 0); const m = Math.sqrt(-2 * Math.log(s) / s); gspare = v * m; return u * m; };
  const layers = cfg.layers, L = layers.length; let N = cfg.N;
  const sp0 = makeSpecies(cfg.beam.sp), spE = makeSpecies({ kind: 'e' });
  const prim = sp0;
  const dir = cfg.beam.dir, dn = Math.hypot(dir[0], dir[1], dir[2]); const d0 = [dir[0] / dn, dir[1] / dn, dir[2] / dn];
  if (!(d0[2] > 1e-6)) throw new Error('The beam must travel towards +z (dz > 0).');
  const edges = [0]; for (const l of layers) edges.push(edges[edges.length - 1] + l.t);
  const zTot = edges[L];
  const E0 = cfg.beam.E0, sigE = (cfg.beam.fwhm || 0) / 2.3548 * E0;
  const EmaxP = (E0 + 6 * sigE) * 1.02 + 1e-6;
  const EcutP = prim.isE ? ph.EcutE : ph.EcutI;
  const Emin = Math.max(Math.min(EcutP, ph.EcutE) * 0.4, 1e-5);
  const needE = prim.isE || (ph.delta && ph.deltaTrack);
  const EmaxE = prim.isE ? EmaxP : Math.min(tmaxHeavy(EmaxP, prim.M), EmaxP) * 1.02 + 1e-6;
  const tabs = [[], []];
  layers.forEach((ly, l) => {
    if (ly.mat) {
      tabs[0][l] = buildTable(ly.mat, prim, ph, prim.isE ? Math.min(Emin, EcutP * 0.4) : Math.min(Emin, EcutP * 0.4), EmaxP);
      tabs[1][l] = needE ? (prim.isE ? tabs[0][l] : buildTable(ly.mat, spE, ph, Emin, Math.max(EmaxE, Emin * 4))) : null;
    } else { tabs[0][l] = null; tabs[1][l] = null; }
  });
  // depth binning
  const nb = cfg.depthBins || 400;
  const est = estimateDepth(layers.map((l, i) => ({ t: l.t })), tabs[0], E0);
  let zmaxP = cfg.depthMax > 0 ? cfg.depthMax : Math.min(zTot, est * 1.15 + 1e-7);
  if (!(zmaxP > 0)) zmaxP = zTot;
  const dzb = zmaxP / nb;
  const dLoss = new Float64Array(nb), dEl = new Float64Array(nb), dNu = new Float64Array(nb), dSe = new Float64Array(nb);
  function spread(arr, za, zb, e) {
    if (e <= 0) return;
    if (za > zb) { const t = za; za = zb; zb = t; }
    if (za >= zmaxP) return;
    const dz = zb - za;
    if (dz < dzb * 1e-3) { const k = Math.floor(za / dzb); if (k >= 0 && k < nb) arr[k] += e; return; }
    const k0 = Math.max(0, Math.floor(za / dzb)), k1 = Math.min(nb - 1, Math.floor(zb / dzb));
    for (let k = k0; k <= k1; k++) { const lo = Math.max(za, k * dzb), hi = Math.min(zb, (k + 1) * dzb); if (hi > lo) arr[k] += e * (hi - lo) / dz; }
  }
  // tallies
  const endX = new Float32Array(N), endY = new Float32Array(N), endZ = new Float32Array(N), pathL = new Float32Array(N), endE = new Float32Array(N);
  const status = new Int8Array(N), endLay = new Int8Array(N).fill(-1);
  const dEprim = new Float32Array(N * L), dEdep = new Float32Array(N * L);
  const lay = { depEl: new Float64Array(L), depNu: new Float64Array(L), depSe: new Float64Array(L), lossPrim: new Float64Array(L), nDelta: new Float64Array(L), nReact: new Float64Array(L) };
  const esc = { transP: 0, backP: 0, transS: 0, backS: 0, brems: 0, react: 0, in: 0, nTransP: 0, nBackP: 0, nStopP: 0, nReactP: 0 };
  const delE = new FBuf(1 << 15), delL = new FBuf(1 << 15);
  const rxZ = new FBuf(1 << 12), rxL = new FBuf(1 << 12), rxE = new FBuf(1 << 12);
  const pv = new VBuf(4, 1 << 16), pvOff = new Int32Array(N + 1);
  const nDisp = Math.min(cfg.dispTracks != null ? cfg.dispTracks : 60, N), dv = new VBuf(3, 1 << 14), dOff = [], dType = [], dPid = [];
  const trackCap = cfg.trackCap || 200000;     // max stored primary vertices
  let pvStored = 0;
  const eps = 1e-10;

  function layerOf(z) { if (z < 0) return -1; for (let i = 0; i < L; i++) if (z < edges[i + 1]) return i; return L; }
  function rotate(P, cosT, sinT, phi) {
    const cp = Math.cos(phi), sp = Math.sin(phi), ux = P.ux, uy = P.uy, uz = P.uz;
    if (Math.abs(uz) > 0.999999) { P.ux = sinT * cp; P.uy = sinT * sp; P.uz = uz > 0 ? cosT : -cosT; return; }
    const den = Math.sqrt(1 - uz * uz);
    P.ux = ux * cosT + sinT * (ux * uz * cp - uy * sp) / den;
    P.uy = uy * cosT + sinT * (uy * uz * cp + ux * sp) / den;
    P.uz = uz * cosT - sinT * cp * den;
    const nn = Math.hypot(P.ux, P.uy, P.uz); P.ux /= nn; P.uy /= nn; P.uz /= nn;
  }
  let ii = 0, ff = 0;
  const loc = (tb, E) => { let x = (Math.log(E) - tb.l0) / tb.dl; if (x < 0) x = 0; const mx = tb.n - 1.0001; if (x > mx) x = mx; ii = x | 0; ff = x - ii; };
  const gv = a => a[ii] + (a[ii + 1] - a[ii]) * ff;

  function dep(kind, P, li, e, za, zb) {
    if (!(e > 0)) return;
    if (kind === 0) { lay.depEl[li] += e; spread(dEl, za, zb, e); }
    else if (kind === 1) { lay.depNu[li] += e; spread(dNu, za, zb, e); }
    else { lay.depSe[li] += e; spread(dSe, za, zb, e); }
    dEdep[P.pid * L + li] += e;
  }
  // path length of the primary ions per layer and energy bin (cm per primary): the ion flux, to be folded with any cross section afterwards
  const nLb = 1000, LbW = EmaxP / nLb, Lh = new Float64Array(L * nLb);
  function tallyL(li, e0, e1, len) {
    let ea = Math.min(e0, e1), eb = Math.max(e0, e1); const o = li * nLb, span = eb - ea;
    if (ea >= EmaxP) ea = eb = EmaxP * (1 - 1e-12);
    if (span < 1e-12 * eb) { Lh[o + Math.min(nLb - 1, (ea / LbW) | 0)] += len; return; }
    const k0 = Math.min(nLb - 1, (ea / LbW) | 0), k1 = Math.min(nLb - 1, (eb / LbW) | 0);
    for (let k = k0; k <= k1; k++) { const a = Math.max(ea, k * LbW), b = Math.min(eb, (k + 1) * LbW); if (b > a) Lh[o + k] += len * (b - a) / span; }
  }
  function ploss(P, li, e, za, zb) {
    if (!(e > 0)) return;
    lay.lossPrim[li] += e; dEprim[P.pid * L + li] += e; spread(dLoss, za, zb, e);
  }
  const stack = [];
  function finish(P, kind) {         // kind: 'stop','trans','back','react'
    if (P.isPrim) {
      const k = P.pid;
      endX[k] = P.x; endY[k] = P.y; endZ[k] = P.z; endE[k] = P.E; pathL[k] = P.path;
      status[k] = kind === 'stop' ? 0 : kind === 'trans' ? 1 : kind === 'back' ? 2 : 3;
      endLay[k] = P.lay;
      if (kind === 'trans') { esc.transP += P.E; esc.nTransP++; } else if (kind === 'back') { esc.backP += P.E; esc.nBackP++; }
      else if (kind === 'react') { esc.react += P.E; esc.nReactP++; } else esc.nStopP++;
    } else {
      if (kind === 'trans') esc.transS += P.E; else if (kind === 'back') esc.backS += P.E;
    }
  }
  function track(P) {
    const spi = P.isE ? 1 : 0, Ecut = P.isE ? ph.EcutE : ph.EcutI, sp = P.isE ? spE : prim;
    const storeP = P.isPrim && pvStored < trackCap, disp = P.pid < nDisp;
    if (storeP) { pv.push4(P.x, P.y, P.z, P.E); pvStored++; }
    if (disp) { dOff.push(dv.n); dType.push(P.isE ? 1 : 0); dPid.push(P.pid); dv.push4(P.x, P.y, P.z); }
    let guard = 0;
    while (true) {
      if (++guard > 2e6) { finish(P, 'stop'); break; }
      const li = layerOf(P.z);
      if (li < 0) { finish(P, 'back'); break; }
      if (li >= L) { finish(P, 'trans'); break; }
      P.lay = li;
      const zlo = edges[li], zhi = edges[li + 1];
      const db = P.uz > 0 ? (zhi - P.z) / P.uz : P.uz < 0 ? (zlo - P.z) / P.uz : 1e30;
      const tb = tabs[spi][li];
      if (!tb) {           // vacuum
        if (db >= 1e29) { finish(P, 'trans'); break; }
        const zn = P.uz > 0 ? zhi + eps : zlo - eps, s = db;
        P.x += P.ux * s; P.y += P.uy * s; P.z = zn; P.path += s;
        if (storeP) { pv.push4(P.x, P.y, P.z, P.E); pvStored++; } if (disp) dv.push4(P.x, P.y, P.z);
        continue;
      }
      if (P.E <= Ecut) {   // absorbed
        const za = P.z; dep(P.isPrim ? 0 : 2, P, li, P.E, za, za); if (P.isPrim) ploss(P, li, P.E, za, za); P.E = 0;
        finish(P, 'stop'); break;
      }
      loc(tb, P.E);
      const Rr = gv(tb.R);
      if (P.isPrim && P.lastLay !== li) { P.lastLay = li; P.fsr = 0.02 + 0.98 * rnd(); }   // random first step in each layer: step grids of different primaries do not coincide (no comb in the flux spectrum)
      let smax = ph.stepFrac * Rr; if (P.fsr) { smax *= P.fsr; P.fsr = 0; } const smin = 1e-8; if (smax < smin) smax = smin;
      const sigD = gv(tb.sigD), sigEh = gv(tb.sigEl), sigNn = gv(tb.sigN), sigT = sigD + sigEh + sigNn;
      let s = Math.min(smax, db), ev = false;
      if (sigT > 0) { const sh = -Math.log(rnd() + 1e-300) / sigT; if (sh < s) { s = sh; ev = true; } }
      // energy loss
      let Em = P.E; const S1 = gv(tb.Sel) + gv(tb.Sn) + gv(tb.Srad);
      Em = P.E - 0.5 * Math.min(S1 * s, 0.5 * P.E);
      loc(tb, Em);
      const sel = gv(tb.Sel), sn = gv(tb.Sn), sr = gv(tb.Srad), vr = gv(tb.vr), th2 = gv(tb.th2);
      let dEel = sel * s, dEn = sn * s, dEr = sr * s;
      if (ph.straggle && vr > 0) { dEel += Math.sqrt(vr * s) * gauss(); if (dEel < 0) dEel = 0; }
      let dEt = dEel + dEn + dEr;
      if (dEt >= P.E) { const k = P.E / dEt; dEel *= k; dEn *= k; dEr *= k; dEt = P.E; }
      const za = P.z;
      P.x += P.ux * s; P.y += P.uy * s; P.z += P.uz * s; P.path += s;
      const zb = P.z;
      if (P.isPrim) { dep(0, P, li, dEel, za, zb); dep(1, P, li, dEn, za, zb); ploss(P, li, dEt, za, zb); if (!P.isE) tallyL(li, P.E, P.E - dEt, s); }
      else dep(2, P, li, dEel, za, zb);
      esc.brems += dEr;
      P.E -= dEt;
      if (ph.scatter && th2 > 0) {
        const t2 = th2 * s; if (t2 > 1e-14) { const th = Math.sqrt(-t2 * Math.log(rnd() + 1e-300)); rotate(P, Math.cos(th), Math.sin(th), 6.283185307 * rnd()); }
      }
      let ended = false;
      if (ev) {
        loc(tb, Math.max(P.E, tb.Emin));
        const r = rnd() * (gv(tb.sigD) + gv(tb.sigEl) + gv(tb.sigN)), sD = gv(tb.sigD), sE = gv(tb.sigEl);
        if (r < sD) {         // delta ray
          let T;
          if (P.isE) {
            const E = P.E, Tc = ph.Tcut;
            if (Tc < E / 2) {
              const g = 1 + E / ME, gg = (g - 1) / g, h = (2 * g - 1) / (g * g), e0 = Tc / E;
              let eps_;
              for (let k = 0; k < 200; k++) {
                eps_ = 1 / (1 / e0 - rnd() * (1 / e0 - 2));
                const f = 1 / (eps_ * eps_) + 1 / ((1 - eps_) * (1 - eps_)) + gg * gg - h / (eps_ * (1 - eps_));
                if (rnd() * 2.25 < f * eps_ * eps_) break;
              }
              T = eps_ * E;
            } else T = 0;
          } else {
            const Tm = tmaxHeavy(P.E, prim.M), Tc = ph.Tcut, b2 = Math.pow(betaOf(P.E, prim.M), 2);
            T = 0;
            if (Tc < Tm) for (let k = 0; k < 200; k++) { T = 1 / (1 / Tc - rnd() * (1 / Tc - 1 / Tm)); if (rnd() < 1 - b2 * T / Tm) break; }
          }
          if (T > 0 && T < P.E) {
            const Et = P.E + sp.M, p0 = Math.sqrt(P.E * (P.E + 2 * sp.M)), pd = Math.sqrt(T * (T + 2 * ME));
            let cd = T * (Et + ME) / (pd * p0); if (cd > 1) cd = 1;
            const sd = Math.sqrt(1 - cd * cd), phi = 6.283185307 * rnd();
            const Ep = P.E - T;
            // secondary
            if (P.isPrim) { lay.nDelta[li]++; delE.push(T); delL.push(li); ploss(P, li, T, P.z, P.z); }
            if (ph.deltaTrack && T > ph.EcutE) {
              const S2 = { isE: true, isPrim: false, pid: P.pid, x: P.x, y: P.y, z: P.z, ux: P.ux, uy: P.uy, uz: P.uz, E: T, path: 0, lay: li };
              rotate(S2, cd, sd, phi); stack.push(S2);
            } else dep(2, P, li, T, P.z, P.z);
            if (P.isE) {      // recoil of the primary electron (momentum conservation)
              const p1 = Math.sqrt(Ep * (Ep + 2 * ME)); let c1 = (p0 - pd * cd) / p1; if (c1 > 1) c1 = 1; if (c1 < -1) c1 = -1;
              rotate(P, c1, Math.sqrt(1 - c1 * c1), phi + Math.PI);
            }
            P.E = Ep;
          }
        } else if (r < sD + sE) {   // hard elastic scattering
          let acc = 0, ki = 0; const rr = (r - sD);
          for (let k = 0; k < tb.nE; k++) { acc += gv(tb.sigElE[k]); if (rr < acc) { ki = k; break; } ki = k; }
          const uc = tb.ucE[ki][ii], A2 = tb.A2E[ki][ii], rm = tb.rE[ki][ii];
          let u;
          for (let k = 0; k < 100; k++) {
            const a1 = 1 / (uc + A2), a2 = 1 / (2 + A2);
            u = 1 / (a1 - rnd() * (a1 - a2)) - A2;
            if (!P.isE) break;
            const b2 = Math.pow(betaOf(P.E, ME), 2);
            if (rnd() < 1 - b2 * u / 2) break;
          }
          const cm = 1 - u, sm = Math.sqrt(Math.max(1 - cm * cm, 0));
          let cl, sl;
          if (rm > 0) { const den = Math.sqrt(1 + rm * rm + 2 * rm * cm); cl = (cm + rm) / den; sl = sm / den; } else { cl = cm; sl = sm; }
          if (!P.isE && rm > 0) { const fh = gv(tb.fh); if (fh > 0) { const Tr = Math.min(fh * 4 * rm / ((1 + rm) * (1 + rm)) * P.E * u / 2, 0.999 * P.E); if (Tr > 0) { dep(1, P, li, Tr, P.z, P.z); ploss(P, li, Tr, P.z, P.z); P.E -= Tr; } } }
          rotate(P, cl, sl, 6.283185307 * rnd());
        } else {                      // nuclear reaction: the projectile is removed (secondaries not transported in this version)
          rxZ.push(P.z); rxL.push(li); rxE.push(P.E); lay.nReact[li]++;
          finish(P, 'react'); ended = true;
        }
      }
      if (ended) { if (storeP) { pv.push4(P.x, P.y, P.z, P.E); pvStored++; } if (disp) dv.push4(P.x, P.y, P.z); break; }
      if (!ev && s >= db * (1 - 1e-12)) { P.z = P.uz > 0 ? zhi + eps : zlo - eps; }
      if (storeP) { pv.push4(P.x, P.y, P.z, P.E); pvStored++; }
      if (disp) dv.push4(P.x, P.y, P.z);
    }
  }
  // main loop, run in chunks so that a worker can stay responsive
  const spot = cfg.beam.spot || 0, div = (cfg.beam.div || 0) / 2.3548;
  let kNext = 0;
  function runChunk(ms) {
    const t0 = Date.now();
    while (kNext < N && Date.now() - t0 < ms) {
      const k = kNext++;
      pvOff[k] = pv.n;
      let E = E0 + sigE * gauss(); if (!(E > 0)) E = E0;
      esc.in += E;
      let x = 0, y = 0;
      if (spot > 0) { const r = spot * Math.sqrt(rnd()), a = 6.283185307 * rnd(); x = r * Math.cos(a); y = r * Math.sin(a); }
      const P = { isE: prim.isE, isPrim: true, pid: k, x, y, z: 0, ux: d0[0], uy: d0[1], uz: d0[2], E, path: 0, lay: 0 };
      if (div > 0) { const th = Math.sqrt(-2 * div * div * Math.log(rnd() + 1e-300)); rotate(P, Math.cos(th), Math.sin(th), 6.283185307 * rnd()); if (P.uz <= 0) { P.ux = d0[0]; P.uy = d0[1]; P.uz = d0[2]; } }
      stack.length = 0; stack.push(P);
      while (stack.length) { const Q = stack.pop(); track(Q); }
    }
    return kNext >= N;
  }
  function result() {
    const Ndone = kNext;
    pvOff[Ndone] = pv.n;
    return {
      N: Ndone, L, edges, zmaxP, nb, dzb, E0, prim: sp0, d0,
      depth: { loss: dLoss, el: dEl, nu: dNu, se: dSe },
      endX: endX.slice(0, Ndone), endY: endY.slice(0, Ndone), endZ: endZ.slice(0, Ndone), endE: endE.slice(0, Ndone), pathL: pathL.slice(0, Ndone), status: status.slice(0, Ndone), endLay: endLay.slice(0, Ndone),
      dEprim: dEprim.slice(0, Ndone * L), dEdep: dEdep.slice(0, Ndone * L), lay, esc,
      delE: delE.done(), delL: delL.done(), rx: { z: rxZ.done(), l: rxL.done(), E: rxE.done() },
      Lh, nLb, LbW,
      pv: pv.done(), pvOff: pvOff.slice(0, Ndone + 1), pvStored,
      disp: { v: dv.done(), off: Int32Array.from(dOff), type: Int8Array.from(dType), pid: Int32Array.from(dPid) },
      phys: ph
    };
  }
  return { runChunk, result, total: N, done: () => kNext };
}

function stopTable(mat, spec, Emin, Emax, n, ph) {
  const sp = makeSpecies(spec), rho = mat.rho;
  const out = { E: [], Se: [], Sn: [], Srad: [], Stot: [], R: [] };
  const pp = Object.assign({ Tcut: 1e9, delta: false, nucStop: true, reactions: false, scatter: false, Nh: 0 }, ph || {});
  const tb = buildTable(mat, sp, pp, Emin, Emax, Math.max(8, n / Math.log10(Emax / Emin)));
  for (let i = 0; i < tb.n; i++) {
    out.E.push(tb.E[i]); out.Se.push(tb.Se[i] / rho); out.Sn.push(tb.Sn[i] / rho); out.Srad.push(tb.Srad[i] / rho); out.Stot.push(tb.Stot[i] / rho); out.R.push(tb.R[i] * rho);
  }
  return out;
}
