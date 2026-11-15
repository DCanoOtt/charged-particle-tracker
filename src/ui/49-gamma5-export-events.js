/* ====== γ spectrum: CSV export and event handlers (continued from 49-gamma.js) ====== */
function gamCtx() {
  // chart context of the spectrum as shown (for the exports)
  const G = GAM;
  if (!G || !G.o) return null;
  const div = document.createElement('div');
  div.innerHTML = renderChart(G.o).svg;
  return {
    o: G.o,
    geom: G.geom,
    svg: div.querySelector('svg'),
    base: 'cpt_gamma_layer' + (G.li + 1) + '_' + Math.round(G.B.t1e) + 's_' + Math.round(G.B.t2e) + 's',
    key: 'gam'
  };
}
function gamCsvSpectrum() {
  const B = GAM.B,
    ns = B.nodes.filter(n => n.de > 0);
  let t =
    '# gamma spectrum, layer ' +
    (GAM.li + 1) +
    ', slice ' +
    B.t1e +
    ' - ' +
    B.t2e +
    ' s; detected counts per 1 keV channel in the slice (efficiency and resolution applied)\nchannel_low_keV,total_detected,total_emitted_with_resolution,' +
    (B.pu && B.pu.spec ? 'total_with_pileup,pileup_sum_pulses,' : '') +
    ns.map(n => dkNuc(n.id)).join(',') +
    '\n';
  for (let c = 0, nr = B.pu && B.pu.spec ? B.pu.spec.length : B.nch; c < nr; c++)
    t +=
      c +
      ',' +
      (B.tot[c] || 0) +
      ',' +
      (B.totE[c] || 0) +
      (B.pu && B.pu.spec ? ',' + B.pu.spec[c] + ',' + B.pu.sums[c] : '') +
      (ns.length ? ',' + ns.map(n => n.arr[c] || 0).join(',') : '') +
      '\n';
  return t;
}
function gamCsvTable() {
  const B = GAM.B;
  let t =
    '# gamma emission of the isotopes in the slice ' +
    B.t1e +
    ' - ' +
    B.t2e +
    ' s (live time ' +
    B.live +
    ' s), layer ' +
    (GAM.li + 1) +
    '\nnuclide,half_life_s,decays,gamma_per_decay,gamma_emitted,gamma_emitted_per_s,detected_counts,detected_per_s,strongest_lines_keV_emitted\n';
  for (const n of B.nodes)
    t +=
      [dkNuc(n.id), n.T || '', n.D, n.gper, n.em, n.em / B.live, n.de, n.de / B.live, '"' + n.top.map(l => l[0] + ' (' + l[1] + ')').join('; ') + '"'].join(
        ','
      ) + '\n';
  return t + '# total,,' + B.sumD + ',,' + B.emitted + ',' + B.emitted / B.live + ',' + B.detected + ',' + B.detected / B.live + '\n';
}
document.addEventListener('click', async ev => {
  if (ev.target.closest('[data-gecfg]')) {
    $('#geCfgWin').hidden = false;
    gamCfgRender();
    return;
  }
  if (ev.target.closest('[data-gcfgclose]') || (ev.target.id === 'geCfgWin' && _mdownId === 'geCfgWin')) {
    $('#geCfgWin').hidden = true;
    return;
  }
  if (ev.target.closest('[data-gclose]') || (ev.target.id === 'gamWin' && _mdownId === 'gamWin')) {
    $('#gamWin').hidden = true;
    return;
  }
  const bs = ev.target.closest('[data-gsave-eff],[data-gsave-sim]');
  if (bs && GAM && GAM.subO) {
    const key = bs.hasAttribute('data-gsave-eff') ? 'eff' : 'sim',
      k = bs.getAttribute('data-gsave-' + key),
      sel = key === 'eff' ? '#gamMsgEff' : '#gamMsgSim',
      m = $(sel);
    if (m) m.textContent = 'Preparing…';
    try {
      const div = document.createElement('div');
      div.innerHTML = renderChart(GAM.subO[key]).svg;
      const ctx = {
          o: GAM.subO[key],
          geom: GAM.subG[key],
          svg: div.querySelector('svg'),
          base: 'cpt_gamma_layer' + (GAM.li + 1) + '_' + (key === 'eff' ? 'efficiency' : 'simulated_ge'),
          key: 'gam-' + key
        },
        it = await secBuild(ctx, k);
      await saveFile(it.name, it.data, k, sel);
    } catch (e) {
      if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
    }
    return;
  }
  const b = ev.target.closest('[data-gsave]');
  if (b && GAM) {
    const k = b.dataset.gsave,
      m = $('#gamMsg');
    if (m) m.textContent = 'Preparing…';
    try {
      if (k === 'spec') await saveFile(gamCtx().base + '_spectrum.csv', gamCsvSpectrum(), 'CSV', '#gamMsg');
      else if (k === 'tab') await saveFile(gamCtx().base + '_isotopes.csv', gamCsvTable(), 'CSV', '#gamMsg');
      else {
        const it = await secBuild(gamCtx(), k);
        await saveFile(it.name, it.data, k, '#gamMsg');
      }
    } catch (e) {
      if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
    }
    return;
  }
  {
    const w = ev.target.closest && ev.target.closest('#geCfgWin');
    if (w) {
      const gc = ev.target.closest('[data-gcurve]');
      if (gc) {
        const g = gamSet(),
          k = gc.dataset.gcurve;
        if (k === 'default') g.curve = gamDefaultCurve(g);
        else if (k === 'ge') g.curve = gamGeCurveText(g);
        else if (k === 'savedef') {
          g.curveDef = String(g.curve || '').trim() ? g.curve : '';
          saveLocal();
          gamCfgRender();
          const m = $('#gcurveMsg');
          if (m) m.textContent = g.curveDef ? 'Saved: this curve is now the default.' : 'The box is empty: the Ge-model values stay the default.';
          return;
        } else if (k === 'cleardef') {
          g.curveDef = '';
          saveLocal();
          gamCfgRender();
          return;
        }
        gamCfgChanged(true);
        return;
      }
      const zr = ev.target.closest('[data-cgzr]');
      if (zr) {
        GEC.z = {};
        gamCfgEff();
        return;
      }
      const ga = ev.target.closest('[data-cgall]');
      if (ga) {
        (GEC.names || []).forEach(n => {
          GEC.hidden[n] = ga.dataset.gstate === 'off';
        });
        gamCfgEff();
        return;
      }
      const lg = ev.target.closest('[data-cleg]');
      if (lg) {
        const n = lg.dataset.cleg;
        GEC.hidden[n] = !GEC.hidden[n];
        gamCfgEff();
        return;
      }
    }
  }
  if (!GAM) return;
  if (ev.target.closest('[data-greset]')) {
    GAM.xmin = GAM.xmax = GAM.ymin = GAM.ymax = '';
    gamRender();
    return;
  }
  const zr = ev.target.closest('[data-gzreset]');
  if (zr) {
    GAM[zr.dataset.gzreset === 'eff' ? 'ze' : 'zs'] = {};
    gamRender();
    return;
  }
  const ga = ev.target.closest('[data-gall]');
  if (ga) {
    const nm = (GAM.names && GAM.names[ga.dataset.gall]) || [];
    nm.forEach(n => {
      GAM.hidden[n] = ga.dataset.gstate === 'off';
    });
    gamRender();
    return;
  }
  const lg = ev.target.closest('[data-gleg]');
  if (lg) {
    const n = lg.dataset.gleg;
    GAM.hidden[n] = !GAM.hidden[n];
    gamRender();
  }
});
document.addEventListener('keydown', ev => {
  if (ev.key !== 'Escape') return;
  if ($('#geCfgWin') && !$('#geCfgWin').hidden) $('#geCfgWin').hidden = true;
  else if ($('#measSelWin') && !$('#measSelWin').hidden) $('#measSelWin').hidden = true;
  else if ($('#measWin') && !$('#measWin').hidden) $('#measWin').hidden = true;
  else if ($('#gamWin') && !$('#gamWin').hidden) $('#gamWin').hidden = true;
});
document.addEventListener('change', ev => {
  const cy = ev.target.closest && ev.target.closest('#geCfgWin [data-cgyl]');
  if (cy) {
    GEC.yl = cy.value;
    GEC.z = Object.assign({}, GEC.z, { ymin: '', ymax: '' });
    gamCfgEff();
    return;
  }
  const ce = ev.target.closest && ev.target.closest('#geCfgWin [data-gk]');
  if (ce) {
    const g = gamSet(),
      k = ce.dataset.gk;
    g[k] = GAM_BOOL.includes(k) ? ce.value === 'true' : ce.type === 'number' ? +ce.value : ce.value;
    gamCfgChanged(k === 'effMode' || k === 'ftpMode');
    return;
  }
  const yl = ev.target.closest && ev.target.closest('#gamWin [data-gyl]');
  if (yl && GAM) {
    GAM.yl = GAM.yl || {};
    GAM.yl[yl.dataset.gyl] = yl.value;
    GAM[yl.dataset.gyl === 'eff' ? 'ze' : 'zs'] = Object.assign({}, GAM[yl.dataset.gyl === 'eff' ? 'ze' : 'zs'], { ymin: '', ymax: '' });
    gamRender();
    return;
  }
  const e = ev.target.closest && ev.target.closest('#gamWin [data-gk]');
  if (e && GAM) {
    const g = gamSet(),
      k = e.dataset.gk;
    g[k] = GAM_BOOL.includes(k) ? e.value === 'true' : e.type === 'number' ? +e.value : e.value;
    if (k === 'ylog') {
      GAM.ymin = GAM.ymax = '';
    }
    saveLocal();
    gamRender();
    return;
  }
  const t = ev.target.closest && ev.target.closest('#gamWin [data-gt]');
  if (t && GAM) {
    const a = +$('#gamWin [data-gt=t1]').value,
      b2 = +$('#gamWin [data-gt=t2]').value;
    if (isFinite(a) && isFinite(b2) && b2 > a && a >= 0) gamOpenAt(GAM.li, a, b2);
  }
});
document.addEventListener('input', ev => {
  // efficiency curve typed in the configuration window
  const ce = ev.target.closest && ev.target.closest('#geCfgWin textarea[data-gk]');
  if (ce) {
    gamSet().curve = ce.value;
    clearTimeout(gamCfgChanged._t);
    gamCfgChanged._t = setTimeout(() => gamCfgChanged(false), 700);
    return;
  }
});
document.addEventListener('input', ev => {
  // the efficiency curve is applied while typing (after a pause)
  const e = ev.target.closest && ev.target.closest('#gamWin textarea[data-gk]');
  if (!e || !GAM) return;
  gamSet().curve = e.value;
  clearTimeout(gamRender._t);
  gamRender._t = setTimeout(() => {
    saveLocal();
    const pos = e.selectionStart;
    gamRender();
    const n = $('#gamWin textarea[data-gk]');
    if (n) {
      n.focus();
      n.setSelectionRange(pos, pos);
    }
  }, 700);
});
