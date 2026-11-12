'use strict';
/* ====== engine ====== */
const ENGINE_SRC = document.getElementById('engine').textContent;
const PT = new Function(ENGINE_SRC + '\n;return PT;')();
const $ = (s, r) => (r || document).querySelector(s), $$ = (s, r) => [...(r || document).querySelectorAll(s)];
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
let _uid = 1; const uid = () => 'i' + (_uid++);
const clone = o => JSON.parse(JSON.stringify(o));
function fmt(v, p) {
  p = p || 4; if (v === null || v === undefined || !isFinite(v)) return '–'; if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-3) return v.toExponential(p - 1).replace(/\.?0+e/, 'e').replace('e+', 'e');
  return String(+v.toPrecision(p));
}
const LEN = { nm: 1e-7, 'µm': 1e-4, mm: 0.1, cm: 1, m: 100 };
const LENU = Object.keys(LEN);
function autoLen(cm) { return cm < 1e-5 ? 'nm' : cm < 0.1 ? 'µm' : cm < 10 ? 'mm' : 'cm'; }
function fmtLen(cm) { const u = autoLen(cm); return fmt(cm / LEN[u], 4) + ' ' + u; }
const EUNITS = { eV: 1e-6, keV: 1e-3, MeV: 1, GeV: 1e3, 'keV/u': 1e-3, 'MeV/u': 1 };
function fmtE(mev) { const a = Math.abs(mev); if (a === 0) return '0 MeV'; if (a < 1e-3) return fmt(mev * 1e6) + ' eV'; if (a < 1) return fmt(mev * 1e3) + ' keV'; if (a < 1e3) return fmt(mev) + ' MeV'; return fmt(mev / 1e3) + ' GeV'; }
