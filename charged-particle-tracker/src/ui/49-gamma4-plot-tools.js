/* ====== γ spectrum: zoom, tooltip, legends, efficiency and simulated-spectrum plots (continued from 49-gamma.js) ====== */
function specTipEl() {
  let t = document.getElementById('specTipWin');
  if (!t) {
    t = document.createElement('div');
    t.id = 'specTipWin';
    t.className = 'mtip';
    t.hidden = true;
    document.body.appendChild(t);
  }
  return t;
}
function specTip(sv, gg, getO) {
  // window next to the cursor with the energy and the counts of the channel under the pointer (all visible histogram curves of the plot); getO() → chart options as drawn
  const tip = specTipEl();
  tip.hidden = true;
  if (!sv || !gg) return;
  const hide = () => {
    tip.hidden = true;
  };
  sv.addEventListener('pointerleave', hide);
  sv.addEventListener('pointerdown', hide);
  sv.addEventListener('pointermove', ev => {
    if (ev.buttons) {
      hide();
      return;
    }
    const q = sv.getBoundingClientRect(),
      px = ((ev.clientX - q.left) / q.width) * gg.W,
      py = ((ev.clientY - q.top) / q.height) * gg.H;
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) {
      hide();
      return;
    }
    const o = getO();
    if (!o) {
      hide();
      return;
    }
    const X = pxData(gg, px, py)[0],
      rows = [];
    let lo = null,
      hi = null;
    for (const se of o.series) {
      if (se.hidden || se.type !== 'hist' || !se.edges) continue;
      const e = se.edges;
      let a = 0,
        b = e.length - 1;
      if (!(X >= e[0] && X < e[b])) continue;
      while (b - a > 1) {
        const m = (a + b) >> 1;
        if (e[m] <= X) a = m;
        else b = m;
      }
      lo = e[a];
      hi = e[a + 1];
      rows.push(
        '<tr><td style="padding-right:10px"><i style="display:inline-block;width:9px;height:9px;border-radius:2px;background:' +
          se.color +
          '"></i> ' +
          esc(se.name) +
          '</td><td style="text-align:right"><b>' +
          fmt(se.y[a], 5) +
          '</b></td></tr>'
      );
    }
    if (!rows.length) {
      hide();
      return;
    }
    tip.innerHTML =
      '<div><b>E = ' +
      fmt(X, 6) +
      ' keV</b> · channel ' +
      fmt(lo, 6) +
      '–' +
      fmt(hi, 6) +
      ' keV</div><table style="border-collapse:collapse;margin-top:3px">' +
      rows.slice(0, 10).join('') +
      '</table><div class="hint" style="margin-top:2px">' +
      esc(o.ylabel || 'counts') +
      '</div>';
    tip.hidden = false;
    const r = tip.getBoundingClientRect();
    let L = ev.clientX + 14,
      T = ev.clientY - r.height - 10;
    if (L + r.width > innerWidth - 6) L = ev.clientX - r.width - 14;
    if (T < 6) T = ev.clientY + 16;
    tip.style.left = Math.max(6, L) + 'px';
    tip.style.top = Math.max(6, T) + 'px';
  });
}
function gamZoom(sv, gg, key, st, redraw) {
  // mouse zoom of a plot below the spectrum: state in GAM[key] = {xmin, xmax, ymin, ymax}
  if (!sv || !gg) return;
  const G = st || GAM,
    rd = redraw || gamRender;
  const pt = ev => {
      const q = sv.getBoundingClientRect();
      return [((ev.clientX - q.left) / q.width) * gg.W, ((ev.clientY - q.top) / q.height) * gg.H];
    },
    clamp = ([px, py]) => [Math.min(Math.max(px, gg.m.l), gg.m.l + gg.iw), Math.min(Math.max(py, gg.m.t), gg.m.t + gg.ih)];
  let a0 = null,
    box = null;
  sv.addEventListener('pointerdown', ev => {
    if (ev.button !== 0 || ev.pointerType === 'touch' || ev.shiftKey || G.pick) return;
    const [px, py] = pt(ev);
    if (px < gg.m.l || px > gg.m.l + gg.iw || py < gg.m.t || py > gg.m.t + gg.ih) return;
    a0 = [px, py];
    sv.setPointerCapture(ev.pointerId);
    box = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    box.setAttribute('fill', 'var(--accent)');
    box.setAttribute('fill-opacity', '.15');
    box.setAttribute('stroke', 'var(--accent)');
    sv.appendChild(box);
  });
  sv.addEventListener('pointermove', ev => {
    if (!a0) return;
    const p = clamp(pt(ev)),
      dx = Math.abs(p[0] - a0[0]),
      dy = Math.abs(p[1] - a0[1]);
    box.setAttribute('x', dx < 6 ? gg.m.l : Math.min(p[0], a0[0]));
    box.setAttribute('width', dx < 6 ? gg.iw : dx);
    box.setAttribute('y', dy < 6 ? gg.m.t : Math.min(p[1], a0[1]));
    box.setAttribute('height', dy < 6 ? gg.ih : dy);
  });
  sv.addEventListener('pointerup', ev => {
    if (!a0) return;
    const p = clamp(pt(ev)),
      s0 = a0,
      dx = Math.abs(p[0] - s0[0]),
      dy = Math.abs(p[1] - s0[1]);
    a0 = null;
    if (box) {
      box.remove();
      box = null;
    }
    if (dx < 6 && dy < 6) return;
    const [xa, ya] = pxData(gg, Math.min(p[0], s0[0]), Math.max(p[1], s0[1])),
      [xb, yb] = pxData(gg, Math.max(p[0], s0[0]), Math.min(p[1], s0[1])),
      z = G[key] || (G[key] = {});
    if (dx >= 6) {
      z.xmin = fmtAx(xa);
      z.xmax = fmtAx(xb);
    }
    if (dy >= 6) {
      z.ymin = fmtAx(ya);
      z.ymax = fmtAx(yb);
    }
    rd();
  });
  sv.addEventListener('pointercancel', () => {
    a0 = null;
    if (box) {
      box.remove();
      box = null;
    }
  });
  sv.addEventListener('dblclick', () => {
    G[key] = {};
    rd();
  });
}
function geInfo(g) {
  const f = E => {
    const v = geEffC(E, g);
    return E + ' keV: ε<sub>peak</sub> ' + fmt(v.peak, 3) + ', ε<sub>tot</sub> ' + fmt(v.tot, 3);
  };
  return f(122) + ' · ' + f(662) + ' · ' + f(1332.5) + ' · relative efficiency (1332.5 keV, 25 cm, vs 3″×3″ NaI) = ' + fmt(geRelEff(g), 3) + ' %';
}
function gamZoomed(z) {
  return !!z && Object.values(z).some(v => v !== '' && v !== undefined);
}
function gamTop(key, ylSel, zoomed, hint) {
  return pcTop({
    y: ['data-gyl="' + key + '"', ylSel],
    reset: ['data-gzreset="' + key + '"', zoomed],
    hint: hint || 'Drag a box to zoom (a thin band zooms one axis) · double-click resets · click a legend entry to show/hide a curve'
  });
}
function gamLegend(names, hidden, key) {
  // legend buttons of one plot + Unselect all / Select all
  return (
    '<div class="legend">' +
    names
      .map(
        n =>
          '<button data-gleg="' + esc(n[0]) + '" class="' + (hidden[n[0]] ? 'off' : '') + '"><i style="background:' + n[1] + '"></i>' + esc(n[0]) + '</button>'
      )
      .join('') +
    pcLegendTools('data-gall="' + key + '" data-gstate="off"', 'data-gall="' + key + '" data-gstate="on"') +
    '</div>'
  );
}
const gamRen = (s, p) =>
  s
    .replace(/data-gleg/g, 'data-' + p + 'leg')
    .replace(/data-gall/g, 'data-' + p + 'gall')
    .replace(/data-gzreset/g, 'data-' + p + 'gzr')
    .replace(/data-gyl/g, 'data-' + p + 'gyl');
function gamEffPlot(g, Emax, ylE, ze, H) {
  // options of the efficiency plot (Ge geometry model: peak and total; curve model: the typed curve)
  const ge = g.effMode === 'ge',
    fp = gamEffFn(g),
    ft = gamEffTotFn(g),
    x = [],
    yp = [],
    yt = [],
    r = Math.max(2, Emax / 20);
  for (let k = 0; k < 140; k++) {
    const E = 20 * Math.pow(r, k / 139);
    x.push(E);
    yp.push(fp(E));
    yt.push(ft(E));
  }
  const c = { x },
    series = [
      { name: 'ε_peak (full energy)', type: 'line', x, y: yp, color: 'var(--accent)', w: 2 },
      { name: 'ε_tot (total)', type: 'line', x, y: yt, color: 'var(--bad)', w: 2, dash: '6 4' }
    ];
  const co = {
    series,
    xlog: true,
    ylog: ylE === 'log',
    xmin: ze.xmin || '20',
    xmax: ze.xmax || String(Emax),
    xlabel: 'Energy (keV)',
    ylabel: 'Efficiency (per emitted γ)',
    title: ge
      ? 'Ge detector efficiency · Ø' + g.geD + ' × ' + g.geL + ' mm, ' + +g.geDist + ' mm from the endcap'
      : 'Detector efficiency curve (typed values; Ge geometry model where empty)'
  };
  co.series.forEach(se => {
    se.hidden = !!H[se.name];
  });
  {
    const x0 = +co.xmin,
      x1 = +co.xmax,
      vis = co.series.filter(se => !se.hidden),
      mx = Math.max(0, ...vis.flatMap(se => se.y.filter((v, i) => c.x[i] >= x0 * 0.999 && c.x[i] <= x1 * 1.001)));
    if (mx > 0) {
      // Y range from the visible curves: the top is the next 1-2-5 value above the maximum, the log axis spans 4 decades (the tiny values at low energy are cut)
      const e = Math.pow(10, Math.floor(Math.log10(mx))),
        m = mx / e,
        top = (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * e;
      if (ylE === 'log') {
        co.ymax = String(+top.toPrecision(2));
        co.ymin = String(+Math.pow(10, Math.floor(Math.log10(top)) - 4).toPrecision(2));
      } else {
        co.ymin = '0';
        co.ymax = String(+top.toPrecision(2));
      }
    }
  }
  if (ze.ymin !== undefined && ze.ymin !== '') co.ymin = ze.ymin;
  if (ze.ymax !== undefined && ze.ymax !== '') co.ymax = ze.ymax;
  return co;
}
const GEC = { hidden: {}, yl: 'log', z: {} }; // state of the efficiency plot in the configuration window
function gamCfgEff() {
  // (re)draw the efficiency plot of the configuration window
  const box = $('#geCfgEff');
  if (!box) return;
  const g = gamSet();
  if (g.effMode !== 'ge' && g.effMode !== 'curve') {
    box.innerHTML = '';
    return;
  }
  const co = gamEffPlot(g, 3000, GEC.yl, GEC.z, GEC.hidden),
    r = renderChart(co),
    zm = gamZoomed(GEC.z);
  const leg = gamRen(gamTop('eff', GEC.yl, zm, null), 'c'),
    leg2 = gamRen(
      gamLegend(
        co.series.map(x => [x.name, x.color]),
        GEC.hidden,
        'eff'
      ),
      'c'
    );
  GEC.names = co.series.map(x => x.name);
  box.innerHTML =
    '<h4>Efficiency curve' +
    's — peak (full-energy) and total' +
    '</h4><p class="hint">Same plot as in the γ-ray spectrum of a time slice · drag a box to zoom, double-click resets</p>' +
    leg +
    '<div class="chart" id="geCfgEffChart">' +
    r.svg +
    '</div>' +
    leg2;
  gamZoom($('#geCfgEffChart svg', box), r.geom, 'z', GEC, gamCfgEff);
}
function gamGeBlock(G, B, g) {
  // efficiency curves and simulated Ge spectrum, shown below the γ plot
  G.yl = G.yl || { eff: 'log', sim: '' };
  const ylE = G.yl.eff || 'log',
    ylS = G.yl.sim || g.ylog,
    H = G.hidden;
  const ze = G.ze || (G.ze = {}),
    zs = G.zs || (G.zs = {}),
    Emax = Math.max(100, Math.min(10000, B.nch)),
    co = gamEffPlot(g, Emax, ylE, ze, H);
  const D = gamDecays(G.sol, G.Q.prod, G.Q.Tirr, G.t1, G.t2).D,
    S = geSimSpectrum(G.sol, D, g, B.nch),
    rate = !!g.rate,
    dv = rate ? B.live : 1,
    sc = a => Float64Array.from(a, v => v / dv),
    edges = new Float64Array(B.nch + 1);
  for (let i = 0; i <= B.nch; i++) edges[i] = i;
  const tot = new Float64Array(B.nch);
  for (let i = 0; i < B.nch; i++) tot[i] = S.peaks[i] + S.cont[i] + S.dbl[i];
  G.sim = { S, tot };
  const so = {
    series: [
      { name: 'simulated Ge spectrum', type: 'hist', edges, y: sc(tot), color: 'var(--ink)', w: 1.4 },
      { name: 'full-energy peaks', type: 'hist', edges, y: sc(S.peaks), color: 'var(--accent)', w: 1.1 },
      { name: 'single Compton + escape peaks', type: 'hist', edges, y: sc(S.cont), color: 'var(--warn)', w: 1.1 },
      { name: 'double Compton (2 scatterings)', type: 'hist', edges, y: sc(S.dbl), color: 'var(--bad)', w: 1.1 }
    ],
    xlog: false,
    ylog: ylS === 'log',
    xlabel: 'Energy (keV, 1 keV per channel)',
    ylabel: rate ? 'Counts per second per keV' : 'Counts per keV in the slice',
    title: 'Simulated Ge spectrum · layer ' + (G.li + 1),
    xmin: zs.xmin || G.xmin,
    xmax: zs.xmax || G.xmax,
    ymin: zs.ymin || '',
    ymax: zs.ymax || ''
  };
  so.series.forEach(se => {
    se.hidden = !!H[se.name];
  });
  G.names.eff = co.series.map(x => x.name);
  G.names.sim = so.series.map(x => x.name);
  const re = renderChart(co),
    rs = renderChart(so);
  G.geoms = { eff: re.geom, sim: rs.geom };
  G.subO = { eff: co, sim: so };
  G.subG = { eff: re.geom, sim: rs.geom };
  return (
    '<details class="cw" open><summary><b>Ge efficiency curves</b> — peak (full-energy) and total; they coincide at low energy, where every interaction ends in full absorption · drag a box to zoom, double-click resets</summary>' +
    gamTop('eff', ylE, gamZoomed(ze)) +
    '<div class="chart" id="geEffChart">' +
    re.svg +
    '</div>' +
    gamLegend(
      co.series.map(x => [x.name, x.color]),
      H,
      'eff'
    ) +
    pcSave({ attr: 'data-gsave-eff', formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'root'], msg: '<span class="hint" id="gamMsgEff"></span>' }) +
    '</details>' +
    '<details class="cw" open><summary><b>Simulated Ge spectrum</b> — peaks, Compton continuum and escape peaks (follows the x range of the plot above until zoomed here · drag a box to zoom, double-click resets)</summary>' +
    gamTop('sim', ylS, gamZoomed(zs)) +
    '<div class="chart" id="geSimChart">' +
    rs.svg +
    '</div>' +
    gamLegend(
      so.series.map(x => [x.name, x.color]),
      H,
      'sim'
    ) +
    pcSave({ attr: 'data-gsave-sim', formats: ['png', 'pdf', 'svg', 'txt', 'csv', 'root'], msg: '<span class="hint" id="gamMsgSim"></span>' }) +
    '</details>'
  );
}
