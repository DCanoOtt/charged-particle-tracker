/* ====== status bar: libraries and data loaded, run state ======
   A slim bar at the top of the page. Each chip has a dot: green = loaded / ready, orange = needed but missing or working, grey = not in use.
   It reads the state of the other modules (stopping-power set, reaction library, decay data, last run) and redraws only when the text changes. */
function sbChips() {
  const out = [],
    add = (state, label, text, tip) => out.push({ state, label, text, tip: tip || '' });
  // stopping powers
  const ss = (S.phys && S.phys.stopSet) || 'auto';
  add(
    'ok',
    'Stopping powers',
    ss === 'semi' ? 'semi-empirical model' : ss === 'nist' ? 'PSTAR / ASTAR / ICRU 73' : 'ICRU 90 / ICRU 73 / NIST',
    'Built-in tables; the choice is in “Stopping powers and models”.'
  );
  // charged-particle reaction library of the current projectile
  const p = PJ();
  if (!p) add('off', 'Reaction library', 'no library for this beam', 'ENDF-6 libraries exist for p, d, t, ³He and α beams.');
  else if (ALIB.busy) add('warn', 'Reaction library', 'indexing “' + ALIB.name + '”…');
  else if (ALIB.counts[p.letter])
    add(
      ALIB.useXs ? 'ok' : 'off',
      'Reaction library',
      (ALIB.name || 'library') + ' · ' + p.name + ' ' + ALIB.counts[p.letter].toLocaleString('en') + ' files' + (ALIB.useXs ? '' : ' (not used)'),
      'Folder “' + ALIB.name + '”; ' + ALIB.data.size + ' isotope(s) read so far.'
    );
  else if (ALIB.nFiles) add('warn', 'Reaction library', 'no ' + p.name + ' files in “' + ALIB.name + '”', ALIB.status);
  else add('warn', 'Reaction library', 'not loaded (' + p.name + ')', 'Choose the TENDL / ENDF-6 folder in “Nuclear reaction data (ENDF-6)”.');
  // decay data
  let src = 'int';
  try {
    src = irrSet().src || 'int';
  } catch (e) {}
  if (src === 'file')
    add(
      DEC.busy ? 'warn' : DEC.nuc ? 'ok' : 'warn',
      'Decay data',
      DEC.busy
        ? 'reading ENDF-6 decay files…'
        : DEC.nuc
          ? 'ENDF-6 decay · ' + (DEC.name || 'folder') + ' · ' + DEC.nuc.toLocaleString('en') + ' nuclides'
          : 'ENDF-6 decay files not loaded',
      DEC.msg
    );
  else if (src === 'ens')
    add(
      ENS.busy ? 'warn' : ENS.nuc ? 'ok' : 'warn',
      'Decay data',
      ENS.busy
        ? 'reading ENSDF files…'
        : ENS.nuc
          ? 'ENSDF · ' + (ENS.name || 'folder') + ' · ' + ENS.nuc.toLocaleString('en') + ' nuclides'
          : 'ENSDF files not loaded',
      ENS.msg
    );
  else add('ok', 'Decay data', 'internal (ICRP-107 / NUBASE2020)', 'Built-in database; uploaded ENDF-6 or ENSDF files can be chosen in “Decay data”.');
  // run
  const stale = !!(document.getElementById('staleBar') && !document.getElementById('staleBar').hidden);
  if (RUN) add('warn', 'Simulation', 'running…');
  else if (stale) add('warn', 'Simulation', 'setup changed: press Run');
  else if (RES) add('ok', 'Simulation', 'up to date');
  else add('off', 'Simulation', 'no run yet');
  return out;
}
let SB_LAST = '';
function sbRender() {
  const el = document.getElementById('sBar');
  if (!el) return;
  let chips;
  try {
    chips = sbChips();
  } catch (e) {
    return;
  }
  const html =
    chips
      .map(
        c =>
          '<span class="sbc ' + c.state + '"' + (c.tip ? ' title="' + esc(c.tip) + '"' : '') + '><i></i><b>' + esc(c.label) + '</b> ' + esc(c.text) + '</span>'
      )
      .join('') +
    '<span class="sbv">v' +
    esc(APP_VERSION) +
    '</span>';
  if (html === SB_LAST) return;
  SB_LAST = html;
  el.innerHTML = html;
}
sbRender();
setInterval(sbRender, 600);
