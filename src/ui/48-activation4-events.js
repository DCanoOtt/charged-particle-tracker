/* ====== activation: events and collapsible cards (continued from 48-activation.js) ====== */
document.addEventListener('input', ev => {
  const el = ev.target.closest && ev.target.closest('#cAct [data-ak]');
  if (!el || el.type === 'checkbox') return;
  const k = el.dataset.ak,
    s = irrSet();
  s[k] = el.type === 'number' ? (el.value === '' ? '' : +el.value) : el.value;
  if (k === 'q' && !(s.q > 0)) s.q = 0;
  saveLocal();
  clearTimeout(ACT_TMR);
  ACT_TMR = setTimeout(actUpdate, 250);
});
document.addEventListener('change', ev => {
  const el = ev.target.closest && ev.target.closest('#cAct [data-ak]');
  if (el) {
    if (el.type === 'checkbox') {
      irrSet()[el.dataset.ak] = el.checked;
      saveLocal();
    }
    clearTimeout(ACT_TMR);
    ACT_TMR = setTimeout(actUpdate, 0);
    return;
  }
  if (ev.target.id === 'decInput' || ev.target.id === 'decInputF') {
    const fl = Array.from(ev.target.files),
      ld = irrSet().src === 'ens' ? ensLoadFiles : dkLoadFiles;
    if (fl.length) ld(fl, fl[0].webkitRelativePath ? fl[0].webkitRelativePath.split('/')[0] : fl.length > 1 ? fl.length + ' files' : fl[0].name);
    ev.target.value = '';
  }
});
document.addEventListener('click', async ev => {
  if (ev.target.id === 'decChoose') {
    const i = $('#decInput');
    if (i) i.click();
    return;
  }
  if (ev.target.id === 'decFiles') {
    const i = $('#decInputF');
    if (i) i.click();
    return;
  }
  const sv = ev.target.closest('[data-asave]');
  if (sv) {
    const el = sv.closest('.actp'),
      m = $('#actMsg'),
      kind = sv.dataset.asave;
    if (m) m.textContent = 'Preparing…';
    try {
      if (kind.startsWith('all-')) {
        const k = kind.slice(4),
          items = [];
        for (const k of ACTP.keys()) {
          const c = actCtx(k);
          if (c) items.push(await secBuild(c, k));
        }
        if (k === 'csv')
          await saveFile(
            'cpt_activation_' + PJ().letter + '_layers.csv',
            items.map(it => '# ===== ' + it.name + ' =====\n' + it.data).join('\n'),
            'CSV',
            '#actMsg'
          );
        else {
          const files = [];
          for (const it of items) files.push({ name: it.name, data: await secBytes(it.data) });
          await saveFile('cpt_activation_' + PJ().letter + '_layers_' + k + '.zip', zipMulti(files), k, '#actMsg');
        }
      } else {
        const c = actCtx(isNaN(+el.dataset.al) ? el.dataset.al : +el.dataset.al),
          it = await secBuild(c, kind);
        await saveFile(it.name, it.data, kind, '#actMsg');
      }
    } catch (e) {
      if (m) m.textContent = 'Export failed: ' + ((e && e.message) || e);
    }
    return;
  }
  if (!ev.target.closest('[data-actcsv]') || !ACT_LAST) return;
  const A = ACT_LAST,
    cool = A.Tc > 0;
  let t =
    '# Charged Particle Tracker — activation; ' +
    fmt(A.rate, 6) +
    ' particles/s for ' +
    A.Tirr +
    ' s' +
    (cool ? ', cooling ' + A.Tc + ' s' : '') +
    '\nlayer,nuclide,half_life_s,decay_modes,origin,yield_per_primary,direct_production_atoms_per_s,atoms_end_irradiation,activity_end_Bq' +
    (cool ? ',activity_after_cooling_Bq' : '') +
    ',gamma_per_decay,gamma_per_s_end' +
    (cool ? ',gamma_per_s_after_cooling' : '') +
    ',peak_activity_Bq\n';
  for (const r of A.rows)
    t +=
      [
        r.li + 1,
        dkNuc(r.id),
        !r.e ? '' : r.e.st === 'stable' ? 'stable' : r.e.T || '',
        '"' + (r.e && r.e.br ? dkLab(r.e.br) : '') + '"',
        r.direct ? 'direct' : 'ingrowth',
        r.y,
        r.Rt,
        r.N0,
        r.lam > 0 ? r.A0 : 0
      ]
        .concat(cool ? [r.lam > 0 ? r.A1 : 0] : [], [r.gt || 0, r.lam > 0 ? r.A0 * (r.gt || 0) : 0], cool ? [r.lam > 0 ? r.A1 * (r.gt || 0) : 0] : [], [r.peak])
        .join(',') + '\n';
  saveFile(A.base + '.csv', t, 'CSV', '#actMsg');
});

/* collapsible cards (Interactions with the target layers, Irradiation and activation): click the heading; the state is remembered */
const FOLD = (() => {
  try {
    return JSON.parse(localStorage.getItem('cpt_fold') || '{}') || {};
  } catch (e) {
    return {};
  }
})();
function foldApply(card) {
  if (!card) return;
  const c = !!FOLD[card.id];
  card.classList.toggle('collapsed', c);
  const h = card.querySelector(':scope > h2[data-fold]');
  if (h) h.setAttribute('aria-expanded', String(!c));
}
function foldToggle(h) {
  const c = h.parentElement;
  if (!c || !c.id) return;
  FOLD[c.id] = !FOLD[c.id];
  foldApply(c);
  try {
    localStorage.setItem('cpt_fold', JSON.stringify(FOLD));
  } catch (e) {
    /* storage may be unavailable */
  }
}
document.addEventListener('click', ev => {
  const h = ev.target.closest && ev.target.closest('h2[data-fold]');
  if (h) foldToggle(h);
});
document.addEventListener('keydown', ev => {
  if (ev.key !== 'Enter' && ev.key !== ' ') return;
  const h = ev.target.closest && ev.target.closest('h2[data-fold]');
  if (h && ev.target === h) {
    ev.preventDefault();
    foldToggle(h);
  }
});
