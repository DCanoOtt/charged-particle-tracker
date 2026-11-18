/* ====== activation: unit conversion, decay data sync and the calculation (continued from 48-activation.js) ====== */
function decUiSync() {
  // the load buttons follow the format chosen in “Decay data”
  const s = irrSet(),
    l = $('#decLbl'),
    b1 = $('#decChoose'),
    b2 = $('#decFiles');
  if (!l) return;
  l.textContent =
    s.src === 'ens'
      ? 'Load ENSDF mass-chain files (ensdf.NNN)'
      : s.src === 'file'
        ? 'Load ENDF-6 decay files (MF 8 / MT 457: ENDF/B, JEFF, JENDL)'
        : 'Load decay data: first choose the format in “Decay data”';
  for (const b of [b1, b2]) if (b) b.disabled = s.src !== 'ens' && s.src !== 'file';
}
function actUpdate() {
  decUiSync();
  const my = ++ACT_GEN;
  clearTimeout(ACT_TMR);
  actCompute(my).catch(e => {
    const b = $('#actRes');
    if (b && my === ACT_GEN) b.innerHTML = '<p class="warn">Calculation failed: ' + esc((e && e.message) || String(e)) + '</p>';
  });
}
async function actCompute(my) {
  const box = $('#actRes'),
    P = PJ();
  if (!box || !P) return;
  const R = RES,
    s = irrSet(),
    msg = m => {
      ACT_CTX = null;
      ACT_LAST = null;
      box.innerHTML = '<p class="hint">' + m + '</p>';
    };
  if (!R) return msg('Run the simulation first.');
  if (!alibF().size)
    return msg(
      'Load the cross-section library for ' +
        P.plural +
        ' first (left panel “Nuclear reaction data (ENDF-6)”: choose the TENDL2025 folder or ' +
        P.lib +
        ', files ' +
        P.letter +
        '-Xx000.tendl).'
    );
  if (R.pj !== P.letter) return msg('The results on screen were computed with another projectile — press Run for the ' + P.name + ' beam.');
  if (R.stale) return msg('The setup (layers, beam or physics) has changed since this run — press Run to update the activation.');
  if (ALIB.busy) return msg('Reading the cross-section files…');
  if (!ALIB.data.size || ![...ALIB.data.values()].some(Boolean)) return msg('No cross-section file was found for the isotopes of the target.');
  if (!R.Lh) return msg('This run has no flux tally — run the simulation again.');
  const I = +s.I,
    Tirr = +s.tirr * ACT_TU[s.tirrU || 'h'],
    Tc = Math.max(0, +s.tcool || 0) * ACT_TU[s.tcoolU || 'h'];
  if (!(I > 0) || !(Tirr > 0)) return msg('Enter a beam current and an irradiation time greater than zero.');
  if (s.src === 'file' && !DEC.nuc) return msg('Load ENDF-6 decay data (button “Choose folder…” or “Choose files…”) or select the internal database.');
  if (s.src === 'ens' && !ENS.nuc)
    return msg('Load the ENSDF mass-chain files (ensdf.NNN; buttons under “Uploaded ENSDF”) or select another decay data source.');
  const q = +s.q > 0 ? +s.q : P.Z,
    rate = s.IU === 'pps' ? I : (I * ACT_IU[s.IU]) / (q * ACT_E),
    Ntot = rate * Tirr; // primaries per second
  box.innerHTML = '<p class="hint">Solving the decay chains…</p>';
  await new Promise(r => setTimeout(r, 20));
  if (my !== ACT_GEN) return;
  const F = alibFold(R),
    layers = s.lay === 'all' || +s.lay >= S.layers.length ? S.layers.map((l, i) => i) : [+s.lay],
    src = s.src === 'file' ? 'file' : s.src === 'ens' ? 'ens' : 'int',
    info = id => dkInfo(id, src, !!s.fb),
    rows = [],
    plots = [];
  let nFile = 0,
    nInt = 0;
  for (const li of layers) {
    const m0 = new Map(),
      m = new Map();
    for (const row of F.rows) {
      if (row.t.layer !== li) continue;
      const r = alibTotal(F, row);
      for (const [id, y] of r.iso) m0.set(id, (m0.get(id) || 0) + y);
    }
    if (src !== 'int')
      for (const id of m0.keys())
        if (id % 10 && !info(id)) {
          // libraries read file by file: the isomer's data must be read before deciding that it has none
          if (src === 'file' && DEC.idx.has(id) && !DEC.map.has(id) && !DEC.miss.has(id)) await decReadId(id);
          else if (src === 'ens' && !ENS.map.has(id) && !ENS.miss.has(id) && (ENS.rest || ENS.idx.size)) await ensReadId(id);
        }
    if (src === 'file') gamUnits(DEC.map);
    for (const [id, y] of m0) {
      let k = id;
      if (k % 10 && !info(k)) k -= k % 10;
      m.set(k, (m.get(k) || 0) + y);
    } // an isomer without decay data is counted as ground state
    const prod = new Map();
    for (const [id, y] of m) if (y > 0) prod.set(id, y * rate);
    if (!prod.size) continue;
    if (src === 'file' && s.chain) {
      await decPrefetch(Array.from(prod.keys()), !!s.fb, m2 => {
        if (my === ACT_GEN) box.innerHTML = '<p class="hint">' + esc(m2) + '</p>';
      });
      if (my !== ACT_GEN) return;
    } else if (src === 'ens') {
      await ensPrefetch(Array.from(prod.keys()), !!s.fb, m2 => {
        if (my === ACT_GEN) box.innerHTML = '<p class="hint">' + esc(m2) + '</p>';
      });
      if (my !== ACT_GEN) return;
    } else if (src === 'file') {
      await decPrefetch(
        Array.from(prod.keys()).filter(() => false),
        !!s.fb
      );
      for (const id of prod.keys()) if (!DEC.map.has(id) && !DEC.miss.has(id) && DEC.idx.has(id)) await decReadId(id);
    }
    const sol = await dkSolve({ prod, info, Tirr, Tc, chain: !!s.chain });
    if (my !== ACT_GEN) return;
    const last = sol.nI - 1,
      lastC = sol.nT - 1;
    sol.ids.forEach((id, i) => {
      const e = sol.info[i],
        lam = sol.lam[i],
        direct = prod.get(id) || 0,
        A0 = sol.act[i][last],
        A1 = Tc > 0 ? sol.act[i][lastC] : A0,
        N0 = sol.Nirr[i],
        N1 = sol.Ncool ? sol.Ncool[i] : N0;
      if (e && e.src !== 'int') nFile++;
      else if (e) nInt++;
      let peak = 0;
      for (let k = 0; k < sol.nT; k++) if (sol.act[i][k] > peak) peak = sol.act[i][k];
      const gt = e ? (e.gtot || 0) + (e.gcTot || 0) : 0,
        xt = (e && e.xtot) || 0;
      rows.push({ li, id, y: m.get(id) || 0, Rt: direct, e, lam, direct: direct > 0, A0, A1, N0, N1, peak, gt, xt });
    });
    plots.push({ li, sol, prod });
  }
  if (my !== ACT_GEN) return;
  // plots: the strongest nuclei of every layer (peak activity), time in the unit of the irradiation
  const uT = s.tirrU || 'h',
    fT = ACT_TU[uT],
    nTop = Math.max(1, Math.min(30, Math.round(+s.top || 8))),
    pkey = [P.letter, Tirr, Tc, s.src, s.qty, s.bplus].join('|');
  ACTP.clear();
  for (const { li, sol, prod } of plots) {
    const hasGam = src !== 'int' && sol.info.some(e => e && (e.gam || e.gc)),
      gm = s.qty === 'gam' && hasGam,
      bm = s.qty === 'bet',
      hasXr = src !== 'int' && sol.info.some(e => e && e.xr),
      xm = s.qty === 'xr' && hasXr,
      xt = sol.info.map(e => (e && e.xtot) || 0),
      atoms = s.qty === 'num',
      gt = sol.info.map(e => (e ? (e.gtot || 0) + (e.gcTot || 0) : 0)),
      pp = Math.max(0, Math.min(1, (+s.bplus || 0) / 100)),
      bt = sol.info.map(e => dkBeta(e, pp));
    const val = atoms
      ? sol.num
      : xm
        ? sol.act.map((a, i) => (xt[i] > 0 ? Float64Array.from(a, v => v * xt[i]) : new Float64Array(a.length)))
        : bm
          ? sol.act.map((a, i) => (bt[i] > 0 ? Float64Array.from(a, v => v * bt[i]) : new Float64Array(a.length)))
          : gm
            ? sol.act.map((a, i) => (gt[i] > 0 ? Float64Array.from(a, v => v * gt[i]) : new Float64Array(a.length)))
            : sol.act;
    const mkSeries = (val, atoms, gm, bm, xm) => {
      const rk = sol.ids
        .map((id, i) => {
          let pk = 0;
          const v = val[i];
          for (let k = 0; k < v.length; k++) if (v[k] > pk) pk = v[k];
          return { id, i, pk };
        })
        .filter(r => r.pk > 0 && (atoms || sol.lam[r.i] > 0))
        .sort((a, b) => b.pk - a.pk);
      const yMax = rk.length ? rk[0].pk : 0,
        yMin = yMax * 1e-9,
        series = [],
        X = sol.t.map(t => t / fT);
      const mk = (name, f, color, extra) => {
        const x = [],
          y = [];
        for (let k = 0; k < X.length; k++) {
          const v = f(k);
          if (v > 0 && (s.ylog !== 'log' || v >= yMin)) {
            x.push(X[k]);
            y.push(v);
          }
        }
        if (x.length) series.push(Object.assign({ name, x, y, color, w: 1.8 }, extra || {}));
      };
      rk.slice(0, nTop).forEach((r, k) => {
        const e = sol.info[r.i];
        mk(dkNuc(r.id) + '  (' + (e && e.st === 'stable' ? 'stable' : 'T½ ' + fmtHl(e && e.T)) + ')', kk => val[r.i][kk], scol(k));
      });
      if ((rk.length > 1 || gm || bm || xm) && !atoms)
        mk(
          gm ? 'total γ emission rate' : xm ? 'total X-ray emission rate' : bm ? 'total β emission rate' : 'sum of all radioactive nuclei',
          k => {
            let a = 0;
            for (const r of rk) a += val[r.i][k];
            return a;
          },
          'var(--ink)',
          { dash: '6 3', w: 2 }
        );
      return series;
    };
    const series = mkSeries(val, atoms, gm, bm, xm);
    ACTP.set(li, {
      li,
      series,
      fT,
      atoms,
      gm,
      hasGam,
      unit: atoms ? 'atoms' : gm ? 'γ/s' : xm ? 'X/s' : bm ? 'β/s' : 'Bq',
      Q: { prod, info, Tirr, Tc, chain: !!s.chain, rate },
      o: {
        series,
        xlabel: 'Time since the start of the irradiation (' + { s: 's', min: 'min', h: 'hours', d: 'days', y: 'years' }[uT] + ')',
        ylabel: atoms
          ? 'Number of atoms'
          : gm
            ? 'Gamma emission rate (γ/s)'
            : xm
              ? 'X-ray emission rate (X-rays/s)'
              : bm
                ? 'Beta emission rate (β/s)'
                : 'Activity (decays/s = Bq)',
        title:
          'Layer ' +
          (li + 1) +
          ' · ' +
          R.names[li] +
          ' — ' +
          P.sym +
          ' beam' +
          (s.qty === 'xr' && !xm
            ? ' (no X-ray records in the decay data: activity shown)'
            : s.qty === 'gam' && !gm
              ? ' (no γ lines in the decay data: activity shown)'
              : ''),
        xlog: false,
        ylog: false
      },
      hint:
        'Irradiation 0 – ' +
        actTime(Tirr) +
        (Tc > 0 ? ', cooling ' + actTime(Tirr) + ' – ' + actTime(Tirr + Tc) : '') +
        '. Time grid (rebuilt for these times and for the half-lives of the chain): exact Bateman propagators at ' +
        sol.nI +
        ' times in the irradiation' +
        (Tc > 0 ? ' and ' + (sol.nT - sol.nI) + ' in the cooling' : '') +
        ' (6 per octave from T/' +
        Math.pow(2, DK_K).toExponential(1) +
        ' up to T, 48 uniform, and ½, 1, 2, 3, 5 half-lives of the leading nuclei); ' +
        sol.nn +
        ' nuclei in the chain network' +
        (sol.fast.some(Boolean) ? ', ' + sol.fast.filter(Boolean).length + ' very short-lived ones in secular equilibrium' : '') +
        '. Click a legend entry to show/hide it.' +
        (s.qty === 'gam' && !hasGam ? ' No γ-line data: the activity is shown (upload ENDF-6 decay data to plot γ/s).' : '')
    });
    if (hasGam && !gm) {
      // second plot: γ emission rate of the decays of this layer (discrete lines of the decay data), per nuclide and in total
      const gval = sol.act.map((a, i) => (gt[i] > 0 ? Float64Array.from(a, v => v * gt[i]) : new Float64Array(a.length))),
        d0 = ACTP.get(li),
        gs = mkSeries(gval, false, true);
      if (gs.length)
        ACTP.set(
          'g' + li,
          Object.assign({}, d0, {
            gam: true,
            series: gs,
            gm: true,
            unit: 'γ/s',
            o: Object.assign({}, d0.o, {
              series: gs,
              ylabel: 'Gamma emission rate (γ/s)',
              title: 'Layer ' + (li + 1) + ' · ' + R.names[li] + ' — γ emission of the decays, ' + P.sym + ' beam'
            })
          })
        );
    }
    for (const k of [li, 'g' + li]) {
      if (!ACTP.has(k)) continue;
      if (!ACTST.has(k) || ACTST.get(k).key !== pkey) ACTST.set(k, { key: pkey, hidden: {}, xmin: '', xmax: '', ymin: '', ymax: '' });
      ACTST.get(k).ylog = s.ylog === 'log';
      ACTST.get(k).xlog = s.xlog === 'log';
    }
  }
  const cool = Tc > 0,
    order = rows.slice().sort((a, b) => a.li - b.li || b.A0 - a.A0 || b.N0 - a.N0);
  const sorted = order.filter(r => s.show === 'all' || (r.e && r.e.st === 'stable' ? false : true));
  const stTxt = r => (!r.e ? 'not in table' : r.e.st === 'stable' ? 'stable' : r.e.T ? fmtHl(r.e.T) : 'T½ n/a'),
    orig = r => (r.direct ? ((r.A0 > 0 || r.N0 > 0) && r.Rt > 0 ? 'direct' : 'direct') : 'ingrowth');
  const gamAvail = sorted.some(r => r.gt > 0);
  const tr = sorted
    .map(
      r =>
        '<tr><td>' +
        (r.li + 1) +
        ' · ' +
        esc(R.names[r.li]) +
        '</td><td><b>' +
        dkNuc(r.id) +
        '</b></td><td class="n">' +
        stTxt(r) +
        '</td><td>' +
        esc(r.e && r.e.br ? dkLab(r.e.br) : '') +
        '</td><td>' +
        orig(r) +
        '</td><td class="n">' +
        (r.direct ? fmt(r.y, 4) : '–') +
        '</td><td class="n">' +
        (r.direct ? fmt(r.Rt, 4) : '–') +
        '</td><td class="n">' +
        fmt(r.N0, 4) +
        '</td><td class="n">' +
        (r.lam > 0 ? fmt(r.A0, 4) : '–') +
        '</td>' +
        (cool ? '<td class="n">' + (r.lam > 0 ? fmt(r.A1, 4) : '–') + '</td>' : '') +
        (gamAvail
          ? '<td class="n">' +
            (r.gt > 0 ? fmt(r.gt, 4) : '–') +
            '</td><td class="n">' +
            (r.gt > 0 && r.lam > 0 ? fmt(r.A0 * r.gt, 4) : '–') +
            '</td>' +
            (cool ? '<td class="n">' + (r.gt > 0 && r.lam > 0 ? fmt(r.A1 * r.gt, 4) : '–') + '</td>' : '')
          : '') +
        '<td class="n">' +
        (r.lam > 0 ? fmt(r.peak, 4) : '–') +
        '</td></tr>'
    )
    .join('');
  const sumA0 = rows.reduce((a, r) => a + (r.lam > 0 ? r.A0 : 0), 0),
    sumA1 = rows.reduce((a, r) => a + (r.lam > 0 ? r.A1 : 0), 0),
    nb = rows.filter(r => r.lam > 0).length,
    nUnk = rows.filter(r => r.e && r.e.st !== 'stable' && !(r.lam > 0)).length,
    nNone = rows.filter(r => !r.e).length,
    nGrow = rows.filter(r => !r.direct).length;
  ACT_LAST = { rows: sorted, rate, Tirr, Tc, base: 'cpt_activation_' + P.letter };
  const srcTxt =
    src === 'ens'
      ? 'uploaded ENSDF data “' + esc(ENS.name) + '” (' + nFile + ' nuclides from ENSDF' + (nInt ? ', ' + nInt + ' from the internal database' : '') + ')'
      : src === 'file'
        ? 'uploaded ENDF-6 decay data “' +
          esc(DEC.name) +
          '” (' +
          nFile +
          ' nuclides from the files' +
          (nInt ? ', ' + nInt + ' from the internal database' : '') +
          ')'
        : 'internal database (ICRP-107 decay data with isomers, NUBASE2020 half-lives; mendeleev 1.3.0 for the other nuclides)';
  const pp = Math.max(0, Math.min(1, (+s.bplus || 0) / 100)),
    cool2 = Tc > 0,
    actInfo = rs => {
      // totals of a set of rows (one layer, or all of them)
      const sA0 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A0 : 0), 0),
        sA1 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A1 : 0), 0),
        gA0 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A0 * r.gt : 0), 0),
        gA1 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A1 * r.gt : 0), 0),
        xA0 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A0 * r.xt : 0), 0),
        xA1 = rs.reduce((a, r) => a + (r.lam > 0 ? r.A1 * r.xt : 0), 0);
      const bal = dkModeBalance(rs, 'A0'),
        balTot = Object.values(bal).reduce((x, y) => x + y, 0),
        betaTot = rs.reduce((x, r) => x + (r.lam > 0 ? r.A0 * dkBeta(r.e, pp) : 0), 0);
      return (
        '<b>Total activity of the radioactive nuclei:</b> ' +
        fmt(sA0, 4) +
        ' Bq at the end of the irradiation' +
        (cool2 ? ', ' + fmt(sA1, 4) + ' Bq after the cooling' : '') +
        '.' +
        (rs.some(r => r.gt > 0)
          ? ' <b>Total γ emission rate (discrete lines, continuous γ spectrum and 511 keV annihilation photons of the decay data):</b> ' +
            fmt(gA0, 4) +
            ' γ/s at the end of the irradiation' +
            (cool2 ? ', ' + fmt(gA1, 4) + ' γ/s after the cooling' : '') +
            '.'
          : src === 'file'
            ? ' (the uploaded decay data have no γ emission for these nuclei)'
            : ' (γ emission needs uploaded ENDF-6 decay data; the internal database has none)') +
        (rs.some(r => r.xt > 0)
          ? ' <b>Total X-ray emission rate:</b> ' +
            fmt(xA0, 4) +
            ' X/s at the end of the irradiation' +
            (cool2 ? ', ' + fmt(xA1, 4) + ' X/s after the cooling' : '') +
            '.'
          : '') +
        (balTot > 0
          ? ' <b>Decays by mode at the end of the irradiation:</b> ' +
            Object.entries(bal)
              .sort((x, y) => y[1] - x[1])
              .map(([k, v]) => k + ' ' + fmt((v / balTot) * 100, 3) + ' %')
              .join(', ') +
            '. <b>Total β emission rate:</b> ' +
            fmt(betaTot, 4) +
            ' β/s (β⁺ share of the EC/β⁺ decays: ' +
            fmt(pp * 100, 3) +
            ' %).'
          : '')
      );
    },
    layIds = [...new Set(rows.map(r => r.li))].sort((a, b) => a - b),
    infoTxt =
      '<p class="hint"><b>' +
      (layIds.length > 1 ? 'All layers together' : 'Layer ' + ((layIds[0] || 0) + 1) + ' · ' + esc(R.names[layIds[0] || 0] || '')) +
      ':</b> ' +
      actInfo(rows) +
      '</p>' +
      (layIds.length > 1
        ? layIds
            .map(l => '<p class="hint"><b>Layer ' + (l + 1) + ' · ' + esc(R.names[l] || '') + ':</b> ' + actInfo(rows.filter(r => r.li === l)) + '</p>')
            .join('')
        : '') +
      '<p class="hint">Decays/s, γ/s, X/s and β/s are different counts and do not add up: a decay emits no β (EC, IT, α), one β or several γ, and nuclides without listed γ lines add decays but no γ/s.</p>';
  box.innerHTML =
    '<p class="hint"><b>Irradiation:</b> ' +
    fmt(I, 4) +
    ' ' +
    (s.IU === 'pps' ? 'particles/s' : { pA: 'pA', nA: 'nA', uA: 'µA', mA: 'mA' }[s.IU]) +
    ' of ' +
    P.plural +
    (s.IU === 'pps' ? '' : ' (charge ' + q + 'e → ' + fmt(rate, 4) + ' particles/s)') +
    ' during ' +
    actTime(Tirr) +
    ' = ' +
    fmt(Ntot, 4) +
    ' particles' +
    (Tc > 0 ? ', then ' + actTime(Tc) + ' of cooling' : '') +
    '. <b>Decay data:</b> ' +
    srcTxt +
    '; decay chains ' +
    (s.chain ? 'followed (' + nGrow + ' nuclei only from ingrowth)' : 'not followed') +
    '.</p>' +
    infoTxt +
    '<div class="btns" style="margin-top:16px"><button class="btn pri" data-gecfg>Configure the Ge detector and DAQ</button> <button class="btn pri" data-measopen title="Simulate a counting experiment with an HPGe detector: counts of selected γ lines vs time during the irradiation and the cooling">Simulate decay measurement</button><span class="hint" id="geCfgSum">' +
    esc(gamCfgSummary(gamSet())) +
    '</span></div>' +
    (plots.length
      ? plots
          .map(
            p =>
              '<div class="secplot actp" data-al="' +
              p.li +
              '"></div>' +
              (ACTP.has('g' + p.li) ? '<div class="secplot actp" data-al="g' + p.li + '"></div>' : '')
          )
          .join('')
      : '<p class="hint">No residual nuclei are produced.</p>') +
    '<span class="hint" id="actMsg"></span>' +
    '<details class="cw" open><summary><b>Residual nuclei and decay products</b> — ' +
    nb +
    ' radioactive with known half-life' +
    (nUnk ? ', ' + nUnk + ' without half-life' : '') +
    (nNone ? ', ' + nNone + ' not in the decay data' : '') +
    ' · sorted by activity</summary><div class="btns"><button class="btn sm" data-actcsv>Save table (CSV)</button></div>' +
    '<div class="tscroll" style="max-height:560px;overflow:auto"><table class="t lc"><tr><th>Layer</th><th>Nuclide</th><th>Half-life</th><th>Decay modes (branching)</th><th>Origin</th><th class="n">Yield per primary</th><th class="n">Direct production (atoms/s)</th><th class="n">Atoms at end of irradiation</th><th class="n">Decays/s at end (Bq)</th>' +
    (cool ? '<th class="n">Decays/s after cooling (Bq)</th>' : '') +
    (gamAvail ? '<th class="n">γ per decay</th><th class="n">γ/s at end</th>' + (cool ? '<th class="n">γ/s after cooling</th>' : '') : '') +
    '<th class="n">Peak decays/s (Bq)</th></tr>' +
    (tr || '<tr><td colspan="14">No residual nuclei.</td></tr>') +
    '</table></div></details>' +
    '<details><summary>How this is computed</summary><p class="hint">Residual nuclei are those of the MF 6 residual-nucleus entries of <i>all</i> non-overlapping reaction channels of the library (like the “Everything produced” table), per primary, from the same flux and cross sections as the Reactions (ENDF) tab; production rate R = yield × particles per second, constant during the irradiation (the charge state only converts an electric current into particles per second). Each layer is solved on its own. With the decay constants λ and the branching fractions b of the selected decay data, N<sub>i</sub>′ = R<sub>i</sub> + Σ<sub>p</sub> b<sub>p→i</sub> λ<sub>p</sub> N<sub>p</sub> − λ<sub>i</sub> N<sub>i</sub> is solved for the whole chain network (all daughters of all products, down to stable nuclei) during the irradiation (R on) and the cooling (R off) with the exact matrix-exponential propagator; activities are A = λN. Isomers: the production of metastable states of the products is used when the cross-section file gives it (LIP of the MF 6 products or the state cross sections of MF 10); the isomers and the isomeric transitions (IT) of the decay data are followed in the chains. Not included: burn-up of the target, reactions of secondary particles, variations of the beam current, delayed-neutron emission as a particle source.</p></details>';
  for (const k of ACTP.keys()) {
    const el = $('.actp[data-al="' + k + '"]', box);
    if (el) actDraw(el);
  }
  ACT_CTX = null;
}
