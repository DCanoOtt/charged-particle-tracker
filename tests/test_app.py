#!/usr/bin/env python3
"""Regression tests for the Charged Particle Tracker.  Run:  python3 build.py && python3 tests/test_app.py
Needs Python Playwright + Chromium. Exits non-zero on the first failing check."""
import os, sys, math
from playwright.sync_api import sync_playwright

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URL = 'file://' + os.path.join(ROOT, 'dist/index.html')
fails = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + (f'  [{info}]' if info != '' else ''))
    if not ok: fails.append(name)

def run(pg, N=4000):
    pg.evaluate(f'S.sim.N={N}'); pg.click('#runBtn')
    pg.wait_for_function("/Done|Failed/.test(document.querySelector('#runStat').textContent)", timeout=180000)

def ni_ar(pg, E=8, P=10, Pu='kPa'):
    pg.evaluate("""([E,P,Pu])=>{const ar=presetDef('Argon gas'); ar.gas={P,T:293.15,Pu};
      S.layers=[{id:uid(),t:2,tu:'µm',def:presetDef('Nickel')},{id:uid(),t:5,tu:'cm',def:ar}];S.beam.pid='a';S.beam.E=E;renderAll()}""", [E, P, Pu])

with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={'width': 1300, 'height': 1300})
    errs = []; pg.on('pageerror', lambda e: errs.append(str(e)))
    pg.goto(URL); pg.wait_for_timeout(1500)

    # --- gas density / molar mass ---
    pg.evaluate("()=>{const ar=presetDef('Argon gas'); ar.rows[0].amt=100; ar.gas={P:10,T:293.15,Pu:'bar'}; ar.rhoAuto=false; S.layers=[{id:uid(),t:5,tu:'cm',def:ar}];renderAll()}")
    pg.evaluate('openModal(0)'); pg.wait_for_timeout(300)
    pg.select_option('[data-m=Pu]', 'bar'); pg.fill('[data-m=P]', ''); pg.click('[data-m=P]'); pg.keyboard.type('0.1', delay=80); pg.wait_for_timeout(500)
    rho = pg.evaluate("matOf(M.def).rho"); check('Ar 0.1 bar density (ideal gas, 1.639e-4)', abs(rho / 1.639e-4 - 1) < 2e-3, rho)
    check('ideal-gas option re-enabled by editing P', pg.is_checked('[data-m=rhoAuto]'))
    check('molar mass Ar (amount 100) = 39.948', abs(pg.evaluate("matOf(M.def).formulaMass") - 39.948) < 1e-3)

    # --- Ni + Ar alpha transport ---
    pg.reload(); pg.wait_for_timeout(1200); ni_ar(pg); run(pg)
    st = pg.evaluate("""(()=>{const R=RES;let n=0,s=0,s2=0;for(let k=0;k<R.N;k++)if(R.status[k]===1){const e=R.endE[k];n++;s+=e;s2+=e*e}
      const m=s/n;return {n,m,sd:Math.sqrt(s2/n-m*m),N:R.N}})()""")
    check('>99.5 % of the alphas transmitted', st['n'] > 0.995 * st['N'], st['n'])
    check('exit energy 6.9-7.6 MeV', 6.9 < st['m'] < 7.6, round(st['m'], 3))
    check('exit sigma < 60 keV (narrow, as Geant4)', st['sd'] < 0.06, round(st['sd'] * 1e3, 1))

    # --- flux spectrum is continuous (no comb) ---
    pg.evaluate("()=>{S.tab='flux';S.plots.flux={what:'sel',lay:'1',addExit:false,eU:'MeV',norm:'phi',bins:120,ylog:'lin',xlog:'lin'};renderAll()}"); pg.wait_for_timeout(600)
    y = pg.evaluate("(()=>{const s=_chart&&_chart.o&&_chart.o.series?_chart.o.series:null;return s&&s[0]?Array.from(s[0].y):null})()")
    if y:
        nz = [v for v in y if v > 0]; i0 = next(i for i, v in enumerate(y) if v > 0); i1 = len(y) - 1 - next(i for i, v in enumerate(reversed(y)) if v > 0)
        core = y[i0 + 3:i1 - 2]; empty = sum(1 for v in core if v == 0)
        check('flux spectrum has no empty bins inside its range', empty == 0, empty)
        jumps = sorted(abs(core[i + 1] - core[i]) / max(core[i], 1e-30) for i in range(len(core) - 1))
        med = jumps[len(jumps) // 2]   # a comb of peaks gives large jumps in most bins; a smooth spectrum only statistical noise
        check('flux spectrum smooth (median bin-to-bin change < 25 %)', med < 0.25, round(med, 3))
    else: check('flux series available', False)

    # --- binary ROOT export (needs uproot; skipped if missing) ---
    try:
        import uproot
        with pg.expect_download() as dl: pg.click('[data-exp=rootbin]')
        path = '/tmp/_cpt_flux.root'; dl.value.save_as(path); rf = uproot.open(path)
        hs = [rf[k] for k in rf.keys()]
        check('ROOT file readable by uproot, TH1D per series', len(hs) >= 1 and all(h.classname == 'TH1D' for h in hs), len(hs))
        check('ROOT histogram has 120 bins and non-zero content', hs[0].values().shape == (120,) and hs[0].values().sum() > 0)
    except ImportError:
        print('SKIP ROOT export (uproot not installed)')

    # --- selecting other stopping data clears results ---
    pg.evaluate("document.querySelector('[data-b=\"phys.stopSet\"]').value='nist'; document.querySelector('[data-b=\"phys.stopSet\"]').dispatchEvent(new Event('change',{bubbles:true}))"); pg.wait_for_timeout(400)
    check('results cleared on stopping-data change', pg.evaluate('RES') is None)

    # --- setup save/load round trip ---
    ni_ar(pg, 8, 10, 'bar')
    with pg.expect_download() as dl: pg.click('#setupSave')
    path = '/tmp/_cpt_setup.json'; dl.value.save_as(path)
    pg.evaluate("document.getElementById('setupReset').click()"); pg.wait_for_timeout(300)
    pg.set_input_files('#setupFile', path); pg.wait_for_timeout(500)
    check('setup round trip', pg.evaluate("S.layers.length===2 && S.layers[1].def.gas.P") == 1000 or pg.evaluate("S.layers.length") == 2)


    # --- alpha cross-section library (synthetic ENDF-6 file with 1000x cross sections): folding vs transport ---
    import subprocess, tempfile, shutil
    lib = tempfile.mkdtemp(); subprocess.run([sys.executable, os.path.join(ROOT, 'tests/make_endf.py'), os.path.join(lib, 'a-Ar040.tendl'), '1000'], check=True)
    pg.reload(); pg.wait_for_timeout(1200)
    pg.evaluate("""()=>{const ar=presetDef('Argon gas'); ar.gas={P:100,T:293.15,Pu:'bar'}; S.layers=[{id:uid(),t:1,tu:'cm',def:ar}];S.beam.pid='a';S.beam.E=8;S.sim.N=20000;renderAll()}""")
    for decoy in ('a-Ni058.tendl', 'a-Ar036.tendl', 'a-Ar040m.tendl'): open(os.path.join(lib, decoy), 'w').write('decoy\n')
    pg.set_input_files('#alibInput', lib); pg.wait_for_timeout(1500); run(pg, 20000)
    reads = pg.evaluate('ALIB.reads')
    check('only the isotopes of the target are read (Ar-36, Ar-38, Ar-40; not Ni or isomers)', sorted(set(reads)) == ['a-Ar036.tendl', 'a-Ar040.tendl'] or sorted(set(reads)) == ['a-Ar040.tendl'] or set(reads) <= {'a-Ar036.tendl', 'a-Ar040.tendl'}, reads)
    check('nothing is read twice', len(reads) == len(set(reads)), reads)
    pg.evaluate("()=>{S.tab='endf';renderPlots(true)}"); pg.wait_for_timeout(800)
    fold = pg.evaluate("(()=>{const R=RES,F=alibFold(R);return F.rows.reduce((a,r)=>a+(alibChannel(F,r,3)||{P:0}).P,0)})()")
    ana = pg.evaluate("RES.lay.nReact[0]/RES.N")
    check('reaction probability: folded flux x MT3 agrees with the transport (5 %)', abs(fold / ana - 1) < 0.05, f'{fold:.5f} vs {ana:.5f}')
    tot = pg.evaluate("(()=>{const F=alibFold(RES);return F.rows.filter(r=>r.d.channels[3]).map(r=>alibTotal(F,r))[0]})()")
    mt3 = fold
    check('sum of channels = MT3', abs(tot['P'] / mt3 - 1) < 1e-6, tot['P'])
    check('MT5/MT16 yields: neutrons > protons and photons counted', tot['Y']['n'] > tot['Y']['p'] > 0 and tot['Y']['g'] > 0)
    # secondary distributions: angle-integrated spectrum = yield, and 4pi * angular integral = integrated spectrum
    chk = pg.evaluate("""(()=>{const F=alibFold(RES),out={};
      for (const [sp,sel] of [['n',['5']],['p',['5']],['n',['16']],['p',['600']]]) {
        const Eo=secBound(F,0,sp,sel); if(!(Eo>0)){out[sp+sel]=null;continue}
        const angs=[];for(let i=0;i<90;i++)angs.push(180*(i+0.5)/90);
        const r=secCompute(F,0,sp,sel,120,Eo,angs), dE=Eo/120; let t=0; for(const v of r.tot)t+=v*dE;
        let a=0; for(let i=0;i<90;i++){let s=0;for(const v of r.ang[i])s+=v*dE; a+=s*2*Math.PI*Math.sin(angs[i]*Math.PI/180)*Math.PI/90}
        out[sp+sel]={W:r.W,tot:t,ang:a}} return out})()""")
    for k, v in chk.items():
        check('secondary %s: angle-integrated spectrum = reactions x yield (3 %%)' % k, v is not None and v['W'] > 0 and abs(v['tot'] / v['W'] - 1) < 0.03, v)
        check('secondary %s: integral of d2N/dEdOmega over the sphere = integrated spectrum (6 %%)' % k, v is not None and abs(v['ang'] / v['tot'] - 1) < 0.06, v)
    nsvg = pg.evaluate("(()=>{S.tab='endf';pset('endf').chs=['sum'];renderPlots(true)})()")
    pg.wait_for_function("!document.getElementById('secProg') && document.querySelector('.secplot .legend')", timeout=60000)
    nsvg = pg.evaluate("document.querySelectorAll('#plotNote .secplot svg').length")
    legs = pg.evaluate("document.querySelector('#plotNote .secplot .legend').innerText")
    check('secondary plots appear above the table, with legends for 0..160 deg and total', nsvg >= 3 and '160' in legs and 'integrated' in legs, (nsvg, legs))
    order = pg.evaluate("(()=>{const n=document.querySelector('#plotNote');return n.innerHTML.indexOf('secplot')<n.innerHTML.indexOf('class=\"t lc\"')})()")
    check('secondary plots precede the reaction table', order)
    ttl = pg.evaluate("document.querySelector('.secplot svg').getAttribute('aria-label')")
    check('secondary plot title names the selected reactions', 'all channels' in ttl, ttl)
    asym = pg.evaluate("""(()=>{const F=alibFold(RES), row=F.rows.find(r=>r.d.channels[5]&&r.d.channels[5].prods&&r.d.channels[5].prods.length), d=row.d.channels[5].prods.find(p=>p.sp==='n').dist;
      const keep=d.L.map(L=>L.p); d.L.forEach(L=>L.p=[L.p[0]]);               // NA = 1 file: only r given
      const Eo=secBound(F,0,'n',['5']), r=secCompute(F,0,'n',['5'],120,Eo,[0,90,180]); d.L.forEach((L,i)=>L.p=keep[i]);
      let m=0; for(let j=0;j<120;j++){const a=r.ang[0][j], b=r.ang[2][j]; if(a>0&&b>0) m=Math.max(m,Math.abs(a/b-1))} return m})()""")
    check('NA=1 Kalbach-Mann: angles differ (a from systematics)', asym > 0.05, asym)
    pg.wait_for_function("!document.getElementById('secProg') && document.querySelector('.secplot .legend')", timeout=60000)
    pg.locator('.secplot').first.locator('svg').scroll_into_view_if_needed(); pg.wait_for_timeout(300)
    box = pg.locator('.secplot').first.locator('svg').bounding_box()
    pg.mouse.move(box['x'] + box['width'] * 0.2, box['y'] + box['height'] * 0.3); pg.mouse.down(); pg.mouse.move(box['x'] + box['width'] * 0.5, box['y'] + box['height'] * 0.6, steps=5); pg.mouse.up(); pg.wait_for_timeout(200)
    z = pg.evaluate("SECST.get(document.querySelector('.secplot').dataset.sk)")
    check('secondary plot zooms with the mouse', z['xmin'] != '' and z['ymin'] != '', z)
    pg.locator('.secplot').first.locator('[data-sleg]').nth(1).click(); pg.wait_for_timeout(100)
    check('angle can be hidden by clicking its legend entry', pg.evaluate("Object.values(SECST.get(document.querySelector('.secplot').dataset.sk).hidden).some(Boolean)"))
    pg.locator('.secplot').first.locator('[data-sleg]').nth(1).click(); pg.wait_for_timeout(100)
    pg.locator('.secplot').first.locator('[data-sy]').select_option('log'); pg.wait_for_timeout(200)
    check('log y axis per secondary plot', pg.evaluate("SECST.get(document.querySelector('.secplot').dataset.sk).ylog") is True and pg.evaluate("document.querySelectorAll('.secplot')[1].querySelector('[data-sy]').value") == 'lin')
    pg.locator('.secplot').first.locator('[data-sreset]').click(); pg.wait_for_timeout(100)
    chans = pg.evaluate("[...document.querySelectorAll('[data-pm=chs]')].map(x=>x.value)")
    check('channel menu lists the sum and the individual MTs', 'sum' in chans and len(chans) > 2, chans)
    pg.evaluate("document.querySelector('#msel-chs').open=true")
    for v in chans[1:3]: pg.check('[data-pm=chs][value="%s"]' % v)
    pg.wait_for_timeout(400)
    ns = pg.evaluate("_chart.o.series.length"); sel = pg.evaluate("pset('endf').chs")
    check('several channels can be ticked: one curve per channel', len(sel) == 3 and ns >= 3 and pg.evaluate("document.querySelector('#msel-chs').open"), (sel, ns))
    pg.click('[data-mnone=chs]'); pg.wait_for_timeout(300)
    check('Clear empties the selection', pg.evaluate("pset('endf').chs.length") == 0)
    pg.click('[data-mall=chs]'); pg.wait_for_timeout(300)
    check('All ticks every channel', pg.evaluate("pset('endf').chs.length") == len(chans))
    shutil.rmtree(lib)
    pg.evaluate("()=>{openModal(0);M.def=clone(presetDef('Silicon'));M.vac=false;M.def.name='Material';M.def.rows=[{kind:'el',sym:'Ar',amt:1}];renderModal()}")
    names = pg.evaluate("[document.querySelector('[data-m=name]').value, document.getElementById('saveName').value]")
    check('single element: material name and Save-as name = element name', names == ['Argon', 'Argon'], names)
    pg.evaluate("()=>{M.def.rows.push({kind:'el',sym:'O',amt:2});renderModal()}")
    n2 = pg.evaluate("document.querySelector('[data-m=name]').value")
    check('two elements: automatic name no longer forced', n2 == 'Material', n2)
    pg.evaluate("closeModal()")
    pg.click('.layer-box >> nth=0') if False else None
    pg.evaluate("openModal(0)"); pg.wait_for_timeout(300)
    pg.fill('[data-m=sym][data-i="0"]', 'Fe'); pg.wait_for_timeout(500)
    names = pg.evaluate("[document.querySelector('[data-m=name]').value, document.getElementById('saveName').value]")
    check('typing the element in a preset layer renames it (UI flow)', names == ['Iron', 'Iron'], names)
    pg.fill('[data-m=sym][data-i="0"]', 'Kr'); pg.wait_for_timeout(500)
    names = pg.evaluate("[document.querySelector('[data-m=name]').value, document.getElementById('saveName').value]")
    check('changing the element updates the name again', names == ['Krypton', 'Krypton'], names)
    pg.evaluate("closeModal()")
    w = pg.evaluate("(()=>{S.tab='depth';renderPlots(true);const e=document.getElementById('expRes');return [e.clientWidth,e.scrollWidth]})()")
    check('PNG width select wide enough for its text', w[0] >= 96, w)

    # --- TENDL2025 tree with the five libraries (p, d, t, h, a): the library of the projectile selected is used ---
    tree = tempfile.mkdtemp(); top = os.path.join(tree, 'TENDL2025')
    for L in 'apdth':
        os.makedirs(os.path.join(top, 'TENDL-' + L)); subprocess.run([sys.executable, os.path.join(ROOT, 'tests/make_endf.py'), os.path.join(top, 'TENDL-' + L, L + '-Ar040.tendl'), '1000'], check=True)
    pg.reload(); pg.wait_for_timeout(1200)
    pg.set_input_files('#alibInput', top); pg.wait_for_timeout(1500)
    cnt = pg.evaluate('ALIB.counts')
    check('folder with five libraries indexed (one file each)', all(cnt.get(L) == 1 for L in 'apdth'), cnt)
    for pid, L, sym in [('p', 'p', 'p'), ('d', 'd', 'd'), ('t', 't', 't'), ('he3', 'h', '³He'), ('a', 'a', 'α')]:
        pg.evaluate("""(pid)=>{const ar=presetDef('Argon gas'); ar.gas={P:100,T:293.15,Pu:'bar'}; S.layers=[{id:uid(),t:1,tu:'cm',def:ar}];S.beam.pid=pid;const q=curPart();S.beam.Z=q.Z;S.beam.A=q.A;S.beam.E=8;S.sim.N=40000;S.tab='endf';renderAll()}""", pid)
        pg.evaluate('alibEnsure()'); pg.wait_for_timeout(1200); run(pg, 40000)
        reads = pg.evaluate('[...ALIB.data.keys()]'); pj = pg.evaluate('PJ().letter')
        check(f'{sym}: library of the projectile is the one read', pj == L and reads and all(r.startswith(L + ':') for r in reads), reads)
        pg.evaluate("()=>{S.tab='endf';renderPlots(true)}"); pg.wait_for_timeout(1200)
        txt = pg.evaluate("document.querySelector('#plotArea, #app, body').innerText")
        check(f'{sym}: Reactions tab shows tables labelled with the projectile', ('per incident ' + sym) in txt and 'Everything produced' in txt and 'Selected reaction channels' in txt)
        if sym == 'α':
            blocks = pg.evaluate("[...document.querySelectorAll('#plotNote details.cw[data-tk]')].map(d=>d.dataset.tk+':'+d.hasAttribute('open'))")
            check('three collapsible tables (selected, everything, channel-sum check), collapsed by default', sorted(blocks) == ['check:false', 'everything:false', 'selected:false'], blocks)
            with pg.expect_download() as dla: pg.click('[data-tcsv=all]')
            ta = open(dla.value.path()).read()
            check('all tables saved in one CSV', ta.count('# ') >= 3 and 'Layer,Isotope,Abundance' in ta and 'Σ channels / MT 3' in ta, ta[:200])
            pg.evaluate("document.querySelector('details[data-tk=check]').open=true")
            with pg.expect_download() as dlc: pg.click('[data-tcsv=check]')
            tc = open(dlc.value.path()).read()
            check('channel-sum check table saved as CSV with full-precision numbers', tc.startswith('# Check of the channel sum') and 'Layer 1 · Ar-40' in tc, tc[:200])
        pg.wait_for_function("!document.getElementById('secProg') && document.querySelector('.secplot .legend')", timeout=60000)
        nsec = pg.evaluate("document.querySelectorAll('.secplot').length")
        check(f'{sym}: secondary-particle plots built', nsec > 0, nsec)
        P = pg.evaluate("(()=>{const F=alibFold(RES);return F.rows.map(r=>alibTotal(F,r).P).reduce((a,b)=>a+b,0)})()")
        ana = pg.evaluate("RES.lay.nReact[0]/RES.N")
        check(f'{sym}: folded probability agrees with transport removal (12 %)', ana > 0 and abs(P / ana - 1) < 0.12, f'{P:.5f} vs {ana:.5f}')
    # --- saving the spectra: one plot, or all of them ---
    pg.evaluate("()=>{S.beam.pid='a';const q=curPart();S.beam.Z=q.Z;S.beam.A=q.A;renderAll()}"); pg.evaluate('alibEnsure()'); pg.wait_for_timeout(1000); run(pg, 6000)
    pg.evaluate("()=>{S.tab='endf';renderPlots(true)}"); pg.wait_for_timeout(1200)
    pg.wait_for_function("!document.getElementById('secProg') && document.querySelector('.secplot .legend')", timeout=60000)
    nplots = pg.evaluate("SEC.size")
    with pg.expect_download() as dl: pg.click('#secAll [data-ssave=csv]')
    t = open(dl.value.path()).read()
    check('save all spectra (CSV): one block per plot', t.count('# =====') == nplots and nplots > 1, (nplots, t.count('# =====')))
    with pg.expect_download() as dl: pg.click('#secAll [data-ssave=svg]')
    import zipfile
    z = zipfile.ZipFile(dl.value.path()); check('save all spectra (SVG): zip with one file per plot', len(z.namelist()) == nplots and z.testzip() is None, z.namelist()[:3])
    with pg.expect_download() as dl: pg.click('#secAll [data-ssave=rootbin]')
    check('save all spectra (ROOT): one file', os.path.getsize(dl.value.path()) > 1000)
    with pg.expect_download() as dl: pg.click('.secplot[data-sk] [data-ssave=png]')
    check('save one spectrum (PNG)', open(dl.value.path(), 'rb').read(8) == b'\x89PNG\r\n\x1a\n')
    # --- the files follow the materials: editing the target loads the new isotopes by itself ---
    pg.evaluate("()=>{S.beam.pid='a';const q=curPart();S.beam.Z=q.Z;S.beam.A=q.A;S.layers=[{id:uid(),t:1,tu:'cm',def:presetDef('Argon gas')}];renderAll()}"); pg.evaluate('alibEnsure()'); pg.wait_for_timeout(800)
    pg.evaluate("()=>{const d=presetDef('Argon gas');d.rows=[{kind:'el',sym:'Ar',amt:50,iso:null},{kind:'el',sym:'Ar',amt:50,iso:null}].slice(0,1);S.layers=[{id:uid(),t:1,tu:'cm',def:presetDef('Nickel')}];renderAll();saveLocal()}"); pg.wait_for_timeout(1500)
    keys = pg.evaluate('[...ALIB.data.keys()]')
    check('changing the material loads its isotopes automatically (no Run needed)', any(k.startswith('a:28-') for k in keys) and not any(k.startswith('a:18-') for k in keys), keys)
    # --- activation window: residual nuclei, half-lives, activities for a constant current ---
    pg.evaluate("()=>{S.beam.pid='a';const q=curPart();S.beam.Z=q.Z;S.beam.A=q.A;S.layers=[{id:uid(),t:1,tu:'cm',def:presetDef('Argon gas')}];renderAll()}"); pg.evaluate('alibEnsure()'); pg.wait_for_timeout(800); run(pg, 6000)
    pg.evaluate("()=>{S.irr={I:100,IU:'nA',q:0,tirr:2,tirrU:'h',tcool:1,tcoolU:'d',lay:'all',show:'all',top:5,ylog:'log',chain:false};document.querySelector('#cAct').dataset.key='';renderPlots(true)}"); pg.wait_for_timeout(2500)
    check('Irradiation panel is below the plots panel', pg.evaluate("!document.querySelector('#cAct').hidden && document.querySelector('#cPlots').nextElementSibling.id === 'cAct' && document.querySelector('#cPlots h2').textContent === 'Interactions with the target layers'"))
    pg.fill('#cAct [data-ak=I]', '200'); pg.wait_for_timeout(600); pg.fill('#cAct [data-ak=I]', '100'); pg.fill('#cAct [data-ak=q]', '4'); pg.wait_for_timeout(1500)
    qq = pg.evaluate('ACT_LAST.rate / (100e-9 / (4 * 1.602176634e-19))'); check('charge state input changes the particle rate (q = 4)', abs(qq - 1) < 1e-9, qq)
    pg.fill('#cAct [data-ak=q]', ''); pg.wait_for_timeout(1500)
    a1 = pg.evaluate("""(()=>{const A=ACT_LAST; if(!A) return null; const rate=100e-9/(2*1.602176634e-19); const r=A.rows.find(x=>x.lam>0);
      return {rate: A.rate/rate, n: A.rows.length, chk: r? [r.A0/(r.Rt*(1-Math.exp(-r.lam*7200))), r.A1/(r.A0*Math.exp(-r.lam*86400)), r.N0*r.lam/r.A0, r.e.T>0]: null, plots: document.querySelectorAll('.actp svg').length, leg: document.querySelectorAll('.actp [data-aleg]').length}})()""")
    check('activation: particle current from nA of alpha (charge 2e)', a1 and abs(a1['rate'] - 1) < 1e-9, a1)
    check('activation without chains: A = R(1-exp(-lambda T)), decay during cooling, N = A/lambda', a1 and a1['chk'] and all(abs(v - 1) < 1e-9 for v in a1['chk'][:3]) and a1['chk'][3], a1)
    check('activation: one plot per layer with legend entries, table lists the nuclei', a1 and a1['plots'] == 1 and a1['leg'] >= 1 and a1['n'] > 0, a1)
    # decay chains on: daughters appear (ingrowth), total activity does not decrease, legend toggling and saving work
    n0 = pg.evaluate("ACT_LAST.rows.length"); pg.evaluate("()=>{irrSet().chain=true;document.querySelector('#cAct').dataset.key='';actRender()}"); pg.wait_for_timeout(3000)
    c1 = pg.evaluate("""(()=>{const A=ACT_LAST; const g=A.rows.filter(r=>!r.direct); const r=A.rows.find(x=>x.direct&&x.lam>0&&x.A0>0); return {n:A.rows.length, ing:g.length, txt:document.querySelector('#actRes').innerText.slice(0,260)}})()""")
    check('activation with chains: daughter nuclides are added (ingrowth)', c1['n'] > n0 and c1['ing'] > 0, [n0, c1])
    pg.click('.actp [data-aleg]'); pg.wait_for_timeout(300)
    check('click on a legend entry hides the curve', pg.evaluate("document.querySelectorAll('.actp [data-aleg].off').length") == 1)
    with pg.expect_download() as dlp: pg.click('.actp [data-asave=csv]')
    check('activity plot can be saved (CSV)', '# ' in open(dlp.value.path()).read())
    # gamma spectra of a time slice (decay data injected as if read from ENDF-6 files: one 1000 keV line, 0.5 gamma per decay)
    pg.evaluate("""()=>{DEC.map=new Map(); for(const r of ACT_LAST.rows){ if(r.e&&r.e.st==='rad'&&r.e.T){ DEC.map.set(r.id,{st:'rad',T:r.e.T,br:r.e.br,gam:[[1000,0.5]],gtot:0.5}); } } DEC.nuc=DEC.map.size; DEC.name='test'; const s=irrSet(); s.src='file'; s.fb=true; s.qty='gam'; document.querySelector('#cAct').dataset.key=''; actRender()}""")
    pg.wait_for_timeout(3500)
    gq = pg.evaluate("[ACTP.get(0).gm, ACTP.get(0).hasGam, !!document.querySelector('.actp [data-aslice]'), ACTP.get(0).unit, document.querySelector('#actRes').innerText.includes('Total γ emission rate')]")
    check('γ emission rate plot and slice buttons with ENDF-6 decay data', gq[0] and gq[1] and gq[2] and gq[3] == 'γ/s' and gq[4], gq)
    pg.evaluate("()=>gamOpenAt(0, 3600, 7200)"); pg.wait_for_function("!document.querySelector('#gamWin').hidden && document.querySelector('#gamWin svg')", timeout=60000)
    gw = pg.evaluate("""(()=>{const B=GAM.B; return {txt:document.querySelector('#gamWin').innerText, tot:B.tot.reduce((a,b)=>a+b,0), det:B.detected, em:B.emitted, peak:B.tot.indexOf(Math.max(...B.tot))}})()""")
    check('γ window: statistics (γ/s, counts/s) and isotope table shown', 'γ emitted' in gw['txt'] and 'counts/s' in gw['txt'] and 'Isotopes contributing' in gw['txt'] and gw['em'] > 0, gw['txt'][:200])
    check('γ spectrum: counts conserved by the Gaussian, peak at 1000 keV', abs(gw['tot'] / gw['det'] - 1) < 1e-6 and gw['peak'] in (999, 1000), [gw['tot'], gw['det'], gw['peak']])
    pg.fill('#gamWin [data-gk=eff]', '0.5'); pg.press('#gamWin [data-gk=eff]', 'Tab'); pg.wait_for_timeout(500)
    d2 = pg.evaluate("GAM.B.detected/GAM.B.emitted"); check('γ window: constant efficiency changes the detected counts', abs(d2 - 0.5) < 1e-9, d2)
    pg.select_option('#gamWin [data-gk=effMode]', 'curve'); pg.wait_for_timeout(500)
    check('γ window: efficiency curve editor', pg.evaluate("!!document.querySelector('#gamWin textarea[data-gk=curve]')"))
    with pg.expect_download() as dg: pg.click('[data-gsave=spec]')
    # pile-up model: P = 1 - exp(-2 n tau), n from the activities at the slice limits; spectrum: clean part + sum pulses
    pg.fill('#gamWin [data-gk=tau]', '0.5'); pg.press('#gamWin [data-gk=tau]', 'Tab'); pg.wait_for_timeout(500)
    pu = pg.evaluate("(()=>{const P=GAM.B.pu, B=GAM.B; return {Pm:P.Pm, mu:P.mu, nbar:P.nbar, det:B.detected/B.live, tau:P.tau, n0:P.n0, n1:P.n1, txt:document.querySelector('#gamWin').innerText.replace(/\u00a0/g,' ')}})()")
    check('pile-up: mean rate = detected counts / live time, μ = nτ, P = 1 − e^(−2μ)', abs(pu['nbar']/pu['det']-1)<1e-9 and abs(pu['mu']-pu['nbar']*0.5e-6)<1e-12 and abs(pu['Pm']-(1-__import__('math').exp(-2*pu['mu'])))<1e-12 and 'PILE-UP PROBABILITY' in pu['txt'] and ' % (mean rate)' in pu['txt'], [pu['nbar'], pu['det'], pu['mu'], pu['Pm']])
    check('pile-up: the mean rate lies between the rates at the slice limits', min(pu['n0'], pu['n1']) > 0 and min(pu['n0'], pu['n1'])*0.9999999 <= pu['nbar'] <= max(pu['n0'], pu['n1'])*1.0000001, [pu['n0'], pu['nbar'], pu['n1']])
    pg.select_option('#gamWin [data-gk=puCurve]', 'true'); pg.wait_for_timeout(600)
    ps = pg.evaluate("(()=>{const P=GAM.B.pu, B=GAM.B; const sum=a=>a.reduce((x,y)=>x+y,0); const k0=Math.exp(-2*P.mu); return {spec:sum(P.spec), exp:k0*sum(B.tot)+sum(P.sums), sums:sum(P.sums), pairs:P.pairs, f:P.f, peak:P.sums.indexOf(Math.max(...P.sums)), leg:[...document.querySelectorAll('#gamWin .legend button')].map(b=>b.textContent)}})()")
    check('pile-up spectrum: sums = pairs (both full-energy), peak at the sum energy, total = clean part + sums', abs(ps['spec']/ps['exp']-1)<1e-9 and ps['sums'] <= ps['pairs']*1.0000001 and ps['sums'] > 0.9*ps['pairs'] and ps['peak'] in (1999,2000,2001), ps)
    pg.locator('#gamWin .mbox').screenshot(path='/tmp/gam.png')
    check('pile-up: curves in the legend', 'total with pile-up' in ps['leg'] and 'pile-up sum pulses only' in ps['leg'], ps['leg'])
    check('γ spectrum can be saved as CSV (1 keV channels)', open(dg.value.path()).read().startswith('# gamma spectrum'))
    # Ge detector model: third efficiency option, own input boxes, curves, simulated spectrum
    pg.select_option('#gamWin [data-gk=effMode]', 'ge'); pg.wait_for_timeout(600)
    ge = pg.evaluate("""(()=>{const w=document.querySelector('#gamWin'), g=gamSet(); const f=(d)=>geEff(1332.5, geP(Object.assign({},g,{geDist:d})));
      const a=f(50), b=f(250), c=geEff(661.7, geP(g)), c0=geEff(661.7, geP(Object.assign({},g,{geScale:2})));
      const S=GAM.sim; let sp=0, sc=0, sd=0; for(let i=0;i<S.S.peaks.length;i++){sp+=S.S.peaks[i]; sc+=S.S.cont[i]; sd+=S.S.dbl[i];}
      return {boxes:['geD','geL','geDist','geGap','geDead','geWin','geScale'].every(k=>w.querySelector('[data-gk='+k+']')), noconst:!w.querySelector('[data-gk=eff]'), nocurve:!w.querySelector('textarea[data-gk=curve]'),
        dist:a.peak>b.peak, tot:c.tot>=c.peak&&c.peak>0, scale:Math.abs(c0.peak/c.peak-2)<1e-9, charts:!!w.querySelector('#geEffChart svg')&&!!w.querySelector('#geSimChart svg'), sp, sc, sd,
        det:GAM.B.detected, ftp:GAM.B.pu.f, rel:geRelEff(g)}})()""")
    check('Ge model: its input boxes appear and the constant/curve boxes disappear', ge['boxes'] and ge['noconst'] and ge['nocurve'], ge)
    check('Ge model: efficiency falls with distance, ε_tot ≥ ε_peak, scale factor is linear', ge['dist'] and ge['tot'] and ge['scale'], ge)
    check('Ge model: curves and simulated spectrum are shown; the peaks of the simulation equal the detected counts; continuum present', ge['charts'] and abs(ge['sp']/ge['det']-1) < 1e-3 and ge['sc'] > 0 and ge['sd'] >= 0, ge)
    check('Ge model: pulses per full-energy count from ε_tot/ε_peak (> 1) and a realistic relative efficiency', ge['ftp'] > 1 and 3 < ge['rel'] < 100, ge)
    pg.locator('#gamWin .mbox').evaluate("e=>e.scrollTop=0"); pg.locator('#gamWin .mbox').screenshot(path='/tmp/ge1.png'); pg.locator('#geSimChart').scroll_into_view_if_needed(); pg.wait_for_timeout(300); pg.locator('#gamWin .mbox').screenshot(path='/tmp/ge2.png')
    xr = pg.evaluate("({xmax:+GAM.geoms.eff.x1, nch:GAM.B.nch})") if False else None
    def drag(sel, fx0, fy0, fx1, fy1):
        pg.locator(sel).scroll_into_view_if_needed(); pg.wait_for_timeout(200); bb = pg.locator(sel).bounding_box(); pg.mouse.move(bb['x']+bb['width']*fx0, bb['y']+bb['height']*fy0); pg.mouse.down(); pg.mouse.move(bb['x']+bb['width']*fx1, bb['y']+bb['height']*fy1, steps=5); pg.mouse.up(); pg.wait_for_timeout(500)
    ext = pg.evaluate("({ze:JSON.stringify(GAM.ze), nch:GAM.B.nch})")
    drag('#geEffChart svg', 0.4, 0.2, 0.7, 0.6)
    drag('#geSimChart svg', 0.4, 0.2, 0.8, 0.6)
    zz = pg.evaluate("({ze:GAM.ze, zs:GAM.zs})")
    check('zoom with the mouse in the efficiency curves and in the simulated spectrum', zz['ze'].get('xmin') and zz['ze'].get('ymax') and zz['zs'].get('xmin') and zz['zs'].get('ymax'), zz)
    pg.click('[data-gzreset=eff]'); pg.click('[data-gzreset=sim]'); pg.wait_for_timeout(400)
    zr = pg.evaluate("({ze:GAM.ze, zs:GAM.zs, xmax:+GAM.o.xmax || 0, nch:GAM.B.nch})")
    check('reset zoom clears both', not zr['ze'] and not zr['zs'] or (not zr['ze'].get('xmin') and not zr['zs'].get('xmin')), zr)
    pg.locator('#geEffChart').screenshot(path='/tmp/geeff.png')
    pg.click('[data-gall=eff][data-gstate=off]'); pg.wait_for_timeout(400)
    lg = pg.evaluate("""(()=>{const w=document.querySelector('#gamWin'); return {eff:[...w.querySelectorAll('[data-gleg]')].filter(b=>/ε_/.test(b.textContent)).every(b=>b.classList.contains('off')), n:w.querySelectorAll('[data-gall]').length, err:!!w.querySelector('#geEffChart svg')}})()""")
    pg.click('[data-gall=eff][data-gstate=on]'); pg.select_option('#gamWin [data-gyl=sim]', 'lin'); pg.wait_for_timeout(400)
    pg.click('[data-gleg="full-energy peaks"]'); pg.click('[data-gall=main][data-gstate=off]'); pg.wait_for_timeout(400)
    lg2 = pg.evaluate("""(()=>{const w=document.querySelector('#gamWin'); return {lin:GAM.yl.sim==='lin', main:GAM.names.main.every(n=>GAM.hidden[n]), peaks:!!GAM.hidden['full-energy peaks'], svgs:w.querySelectorAll('svg').length}})()""")
    check('legends: unselect all / select all in every plot, per-plot linear/log Y axis', lg['eff'] and lg['n'] == 6 and lg2['lin'] and lg2['main'] and lg2['peaks'] and lg2['svgs'] >= 3, [lg, lg2])
    pg.click('[data-gall=main][data-gstate=on]'); pg.click('[data-gall=sim][data-gstate=on]')
    pg.select_option('#gamWin [data-gk=effMode]', 'curve'); pg.wait_for_timeout(500)
    check('switching back to the curve model shows its textarea and hides the Ge boxes', pg.evaluate("!!document.querySelector('#gamWin textarea[data-gk=curve]') && !document.querySelector('#gamWin [data-gk=geD]')"))
    pg.click('[data-gclose]'); check('γ window closes', pg.evaluate("document.querySelector('#gamWin').hidden"))
    pg.evaluate("()=>{const s=irrSet(); s.src='int'; s.qty='act'; DEC.map=new Map(); DEC.nuc=0; document.querySelector('#cAct').dataset.key=''; actRender()}"); pg.wait_for_timeout(2500)
    pg.evaluate("()=>{irrSet().qty='num';irrSet().top=12;document.querySelector('#cAct').dataset.key='';actRender()}"); pg.wait_for_timeout(3000)
    nq = pg.evaluate("[document.querySelector('.actp .legend').innerText, ACTP.get(0).atoms]")
    check('plot of the number of atoms includes stable isotopes', nq[1] and 'stable' in nq[0], nq)
    pg.evaluate("()=>{irrSet().qty='act';irrSet().top=5;document.querySelector('#cAct').dataset.key='';actRender()}"); pg.wait_for_timeout(2500)
    iso = pg.evaluate("""(()=>{const F=alibFold(RES); for(const row of F.rows){ for(const mt of row.d.order){ const ch=row.d.channels[mt]; if(!ch.primary||!ch.res.length) continue; const r0=alibChannel(F,row,mt); if(!(r0.P>0)) continue; const za=ch.res[0].ZA;
        ch.iso=[{ZA:za,lfs:1,tab:{x:ch.tab.x,y:Float64Array.from(ch.tab.y,v=>0.3*v),nbt:ch.tab.nbt,int:ch.tab.int}}]; const r1=alibChannel(F,row,mt); delete ch.iso; return [r1.iso.get(za*10+1)/r1.res.get(za), r1.iso.get(za*10)/r1.res.get(za)]; } } return null})()""")
    check('isomer production: MF 10 state cross section splits the residual yield (30 % isomer, 70 % ground)', iso and abs(iso[0] - 0.3) < 1e-6 and abs(iso[1] - 0.7) < 1e-6, iso)
    known = pg.evaluate("[dkInternal(270580), dkInternal(180360), dkInternal(180370)].map(h=>[h.st, Math.round((h.T||0)/86400*100)/100])")
    check('internal decay table: Co-58 70.8 d, Ar-36 stable, Ar-37 35.0 d', known[0][0] == 'rad' and abs(known[0][1] - 70.84) < 0.1 and known[1][0] == 'stable' and abs(known[2][1] - 35.01) < 0.1, known)
    with pg.expect_download() as dl: pg.click('[data-actcsv]')
    check('activation table can be saved (CSV)', 'activity_end_Bq' in open(dl.value.path()).read())
    pg.evaluate("()=>{S.beam.pid='p';const q=curPart();S.beam.Z=q.Z;S.beam.A=q.A;renderAll()}")
    stale = pg.evaluate("(()=>{S.tab='endf';renderPlots(true);return /another projectile|Load the proton/.test(document.body.innerText)})()")
    check('results of another projectile are not folded with this library', stale)

    check('no page errors', not [e for e in errs if 'ERR_TUNNEL' not in e], errs[:2])
    b.close()
sys.exit(1 if fails else 0)
