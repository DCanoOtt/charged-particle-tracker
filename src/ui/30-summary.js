/* ====== summary ====== */
function renderSummary() {
  const box = $('#cSum'); if (!RES) { box.hidden = true; return; }
  const R = RES, N = R.N, e = R.esc, L = R.L; box.hidden = false;
  const cnt = [0, 0, 0, 0]; for (let i = 0; i < N; i++) cnt[R.status[i]]++;
  const zs = []; for (let i = 0; i < N; i++) if (R.status[i] === 0) zs.push(R.endZ[i]);
  const zd = describe(zs), pct = k => fmt(100 * cnt[k] / N, 3) + ' %';
  let nD = 0; for (let i = 0; i < L; i++) nD += R.lay.nDelta[i];
  const sumA = a => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
  const dEl = sumA(R.lay.depEl) / N, dNu = sumA(R.lay.depNu) / N, dSe = sumA(R.lay.depSe) / N;
  const out = [['Electronic energy loss of the projectile (deposited locally)', dEl], ['Nuclear (elastic recoil) energy loss', dNu], ['Deposited by knock-on electrons', dSe],
  ['Escaping primaries — transmitted', e.transP / N], ['Escaping primaries — backscattered', e.backP / N], ['Escaping knock-on electrons — forward / backward', (e.transS + e.backS) / N], ['Bremsstrahlung (photons, not tracked)', e.brems / N], ['Projectile removed by a nuclear reaction', e.react / N]];
  let tot = 0; out.forEach(o => tot += o[1]); const inp = e.in / N;
  let kp = '<div class="kpis">' +
    '<div class="kpi"><small>Stopped in target</small><b>' + pct(0) + '</b><span>' + cnt[0].toLocaleString('en') + ' of ' + N.toLocaleString('en') + '</span></div>' +
    '<div class="kpi"><small>Transmitted</small><b>' + pct(1) + '</b><span>' + (cnt[1] ? 'mean E = ' + fmtE(e.transP / Math.max(e.nTransP, 1)) : '–') + '</span></div>' +
    '<div class="kpi"><small>Backscattered</small><b>' + pct(2) + '</b><span>' + (cnt[2] ? 'mean E = ' + fmtE(e.backP / Math.max(e.nBackP, 1)) : '–') + '</span></div>' +
    '<div class="kpi"><small>Nuclear reaction</small><b>' + pct(3) + '</b><span>projectile removed</span></div>' +
    '<div class="kpi"><small>Projected range ⟨z⟩</small><b>' + (zd ? fmtLen(zd.mean) : '–') + '</b><span>' + (zd ? 'σ = ' + fmtLen(zd.sd) + ' (' + fmt(100 * zd.sd / zd.mean, 3) + ' %)' : 'no particle stopped') + '</span></div>' +
    '<div class="kpi"><small>Mean total energy deposit</small><b>' + fmtE(sumA(R.layerE) / N) + '</b><span>per primary, incl. δ-rays</span></div>' +
    '<div class="kpi"><small>δ-rays produced</small><b>' + fmt(nD / N, 4) + '</b><span>per primary above ' + fmtE(R.phys.Tcut) + '</span></div>' +
    '<div class="kpi"><small>Run time</small><b>' + fmt(R.secs, 3) + ' s</b><span>' + R.beamName + ', ' + fmtE(R.E0nom) + '</span></div></div>';
  let lt = '<div class="tscroll"><table class="t"><tr><th>Layer</th><th>Thickness</th><th>E dep total</th><th>ion electronic</th><th>ion nuclear</th><th>δ-electrons</th><th>primary loss</th><th>δ-rays</th><th>reactions</th></tr>';
  for (let i = 0; i < L; i++) lt += '<tr><td>' + (i + 1) + ' · ' + esc(R.names[i]) + '</td><td>' + esc(fmtLen(R.edges[i + 1] - R.edges[i])) + '</td><td>' + fmtE(R.layerE[i] / N) + '</td><td>' + fmtE(R.lay.depEl[i] / N) + '</td><td>' + fmtE(R.lay.depNu[i] / N) + '</td><td>' + fmtE(R.lay.depSe[i] / N) + '</td><td>' + fmtE(R.lay.lossPrim[i] / N) + '</td><td>' + fmt(R.lay.nDelta[i] / N, 3) + '</td><td>' + fmt(R.lay.nReact[i], 4) + '</td></tr>';
  lt += '</table></div><p class="hint">Energies are per primary particle. “Primary loss” is the energy the projectile itself lost (including what it gave to δ-rays); “E dep total” is what stays in the layer.</p>';
  let bt = '<div class="tscroll"><table class="t"><tr><th>Energy balance (per primary)</th><th>MeV</th><th>% of input</th></tr><tr><td><b>Beam kinetic energy</b></td><td>' + fmt(inp, 6) + '</td><td>100</td></tr>' +
    out.map(o => '<tr><td>' + o[0] + '</td><td>' + fmt(o[1], 5) + '</td><td>' + fmt(100 * o[1] / inp, 4) + '</td></tr>').join('') +
    '<tr><td><b>Sum of outputs</b> (residual ' + fmt(inp - tot, 2) + ' MeV)</td><td>' + fmt(tot, 6) + '</td><td>' + fmt(100 * tot / inp, 5) + '</td></tr></table></div>';
  box.innerHTML = '<h2>Results summary</h2>' + kp + '<details open><summary>Layers</summary>' + lt + '</details><details><summary>Energy balance</summary>' + bt + '</details>';
}
