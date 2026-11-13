/* ====== helpers ====== */
function matOf(def) { return def ? PT.buildMaterial(def, libMap()) : null; }
function layerCm(l) {
  if (!l.def) { if (LEN[l.tu]) return l.t * LEN[l.tu]; throw new Error('Vacuum gaps need a length unit.'); }
  if (LEN[l.tu]) return l.t * LEN[l.tu];
  const m = matOf(l.def); return (l.tu === 'mg/cm²' ? l.t / 1000 : l.t) / m.rho;
}
function safeLayerCm(l) { try { return layerCm(l); } catch (e) { return NaN; } }
function layerName(l, i) { return l.def ? (l.def.name || 'Layer ' + (i + 1)) : 'Vacuum'; }
function curPart() { return PARTICLES.find(p => p.id === S.beam.pid) || PARTICLES[5]; }
function beamSpec() { const p = curPart(); if (p.kind === 'e') return { kind: 'e' }; if (p.custom) return { kind: 'ion', Z: Math.round(S.beam.Z), A: Math.round(S.beam.A) }; return { kind: 'ion', Z: p.Z, A: p.A }; }
function beamAZ() { const sp = beamSpec(); return sp.kind === 'e' ? { Z: -1, A: 0 } : { Z: sp.Z, A: sp.A }; }
function beamName() { const p = curPart(); if (p.custom) return 'Z=' + Math.round(S.beam.Z) + ', A=' + Math.round(S.beam.A) + ' ion'; return p.name.replace(/ \(.*/, ''); }
function beamE0() { const b = S.beam, perU = /\/u$/.test(b.Eu); const k = EUNITS[b.Eu] || 1; return b.E * k * (perU ? Math.max(beamAZ().A, 1) : 1); }
function dirNorm() { const d = S.beam.dir.map(Number), n = Math.hypot(d[0], d[1], d[2]); return n > 0 ? d.map(v => v / n) : [0, 0, 1]; }
