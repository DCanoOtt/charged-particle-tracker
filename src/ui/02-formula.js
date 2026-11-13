/* ====== formula parser ====== */
function parseFormula(str) {
  let i = 0; const s = str.replace(/\s+/g, '');
  function num() { let j = i; while (j < s.length && /[0-9.]/.test(s[j])) j++; const v = j > i ? parseFloat(s.slice(i, j)) : 1; i = j; return v; }
  function group() {
    const m = {};
    while (i < s.length && s[i] !== ')') {
      if (s[i] === '(') { i++; const sub = group(); if (s[i] !== ')') throw new Error('Missing )'); i++; const k = num(); for (const e in sub) m[e] = (m[e] || 0) + sub[e] * k; }
      else if (/[A-Z]/.test(s[i])) { let sym = s[i++]; if (i < s.length && /[a-z]/.test(s[i])) sym += s[i++]; if (!PT.SYMZ[sym]) throw new Error('Unknown element ' + sym); m[sym] = (m[sym] || 0) + num(); }
      else throw new Error('Unexpected "' + s[i] + '"');
    }
    return m;
  }
  const m = group(); if (i < s.length) throw new Error('Unbalanced )'); if (!Object.keys(m).length) throw new Error('Empty formula');
  return m;
}
