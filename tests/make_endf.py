#!/usr/bin/env python3
"""Writes small synthetic ENDF-6 alpha files in the TENDL layout (for parser tests; not real evaluated data)."""
import math, sys
def f(x):
    if x == 0: return ' 0.000000+0'
    s = '%.6E' % x; m, e = s.split('E'); e = int(e)
    return '%s%s%d' % (m, '+' if e >= 0 else '-', abs(e)) if len(m) + 2 <= 11 else None
def fld(v):
    if isinstance(v, int): return '%11d' % v
    return f(v).rjust(11)
class W:
    def __init__(s, mat): s.mat = mat; s.out = []; s.ns = 0
    def line(s, vals, mf, mt):
        s.ns += 1
        body = ''.join(fld(v) for v in vals) + ' ' * (11 * (6 - len(vals)))
        s.out.append('%66s%4d%2d%3d%5d' % (body, s.mat, mf, mt, s.ns % 100000))
    def end(s, mf, mt):            # mt != 0: SEND of that section (MT=0); mt == 0: FEND (MF=0)
        s.ns += 1; s.out.append('%66s%4d%2d%3d%5d' % ('', s.mat, mf if mt else 0, 0, s.ns))
    def tab1(s, c1, c2, l1, l2, xs, ys, mf, mt, inter=2):
        s.line([c1, c2, l1, l2, 1, len(xs)], mf, mt); s.line([len(xs), inter], mf, mt)
        pts = []
        for x, y in zip(xs, ys): pts += [x, y]
        for i in range(0, len(pts), 6): s.line(pts[i:i + 6], mf, mt)
    def lst(s, c1, c2, l1, l2, data, n2, mf, mt):
        s.line([c1, c2, l1, l2, len(data), n2], mf, mt)
        for i in range(0, len(data), 6): s.line(data[i:i + 6], mf, mt)
    def tab2(s, ne, mf, mt):
        s.line([0.0, 0.0, 0, 0, 1, ne], mf, mt); s.line([ne, 2], mf, mt)

def sigma(E, thr, peak, w=6.0):          # eV-free toy excitation function (barn), E in MeV
    return 0.0 if E < thr else peak * (1 - math.exp(-(E - thr) / 2.0)) * math.exp(-(E - thr) / w)

def make(Z=18, A=40, path=None, scale=1.0):
    ZA = Z * 1000 + A; mat = 8892; w = W(mat)
    Es = [0.25e6 + k * 0.25e6 for k in range(0, 64)]        # eV, 0.25 .. 16 MeV
    E_MeV = [e * 1e-6 for e in Es]
    nonel = [0.8 * (1 - math.exp(-(e - 1.5) / 2)) if e > 1.5 else 0.0 for e in E_MeV]
    mt102 = [0.02 * math.exp(-e / 4) if e > 0.5 else 0.0 for e in E_MeV]
    mt16 = [sigma(e, 6.0, 0.25) for e in E_MeV]
    mt103 = [sigma(e, 2.0, 0.30) for e in E_MeV]               # (a,p) total
    mt600 = [sigma(e, 2.0, 0.18) for e in E_MeV]
    mt649 = [sigma(e, 2.5, 0.12) for e in E_MeV]
    nonel = [v * scale for v in nonel]; mt102 = [v * scale for v in mt102]; mt16 = [v * scale for v in mt16]; mt103 = [v * scale for v in mt103]; mt600 = [v * scale for v in mt600]; mt649 = [v * scale for v in mt649]
    rest = [max(nonel[i] - mt102[i] - mt16[i] - mt600[i] - mt649[i], 0.0) for i in range(len(Es))]
    elastic = [1.5 if e > 0 else 0 for e in E_MeV]
    # MF1/MT451
    w.line([float(ZA), 39.62, -1, 0, 17, 1], 1, 451); w.line([0.0, 0.0, 0, 0, 0, 6], 1, 451); w.line([3.97, 2e8, 1, 0, 20040, 2023], 1, 451); w.line([0.0, 0.0, 0, 0, 283, 58], 1, 451)
    w.out.append('%-66s%4d%2d%3d%5d' % (' %d-Ar-%d IAEA       SYNTHETIC TEST' % (Z, A), mat, 1, 451, 5)); w.end(1, 451); w.end(1, 0)
    for mt, ys in [(2, elastic), (3, nonel), (5, rest), (16, mt16), (102, mt102), (103, mt103), (600, mt600), (649, mt649)]:
        w.line([float(ZA), 39.62, 0, 0, 0, 0], 3, mt); w.tab1(0.0, {16: -7.0e6, 600: -1.0e6, 649: -2.0e6}.get(mt, 0.0), 0, 0, Es, ys, 3, mt); w.end(3, mt)
    w.end(3, 0)
    # MF6 with real structures: MT5 lab frame (LCT 1) LAW 1 Kalbach-Mann; MT16 CM frame (LCT 2) LAW 1 Legendre; MT600/649 CM LAW 2 two-body
    def km_lists(w, mt, T=1.0e6):
        w.line([0.0, 0.0, 2, 2, 1, 2], 6, mt); w.line([2, 2], 6, mt)
        for E in (1e6, 2e7):
            eps = [E * 0.9 * i / 40 for i in range(41)]; fs = [math.exp(-e / T) for e in eps]
            nrm = sum(0.5 * (fs[i] + fs[i + 1]) * (eps[i + 1] - eps[i]) for i in range(40)); data = []
            for e, fv in zip(eps, fs): data += [e, fv / nrm, 0.3, 2.0]
            w.lst(0.0, E, 0, 2, data, 41, 6, mt)
    def leg_lists(w, mt, T=1.5e6):
        w.line([0.0, 0.0, 1, 2, 1, 2], 6, mt); w.line([2, 2], 6, mt)
        for E in (1e6, 2e7):
            eps = [E * 0.5 * i / 40 for i in range(41)]; fs = [math.exp(-e / T) for e in eps]
            nrm = sum(0.5 * (fs[i] + fs[i + 1]) * (eps[i + 1] - eps[i]) for i in range(40)); data = []
            for e, fv in zip(eps, fs): data += [e, fv / nrm, 0.3 * fv / nrm]
            w.lst(0.0, E, 0, 1, data, 41, 6, mt)
    w.line([float(ZA), 39.62, 0, 1, 3, 0], 6, 5)
    for zap, awp, y in [(1, 1.0, 0.5), (1001, 0.99917, 0.3), (2004, 3.9726, 0.2)]:
        w.tab1(float(zap), awp, 0, 1, [1e-5, 2e7], [y, y], 6, 5); km_lists(w, 5)
    w.end(6, 5)
    w.line([float(ZA), 39.62, 0, 2, 3, 0], 6, 16)
    w.tab1(1.0, 1.0, 0, 1, [1e-5, 2e7], [2.0, 2.0], 6, 16); leg_lists(w, 16)
    w.tab1(0.0, 0.0, 0, 2, [1e-5, 8e6, 2e7], [0.0, 1.5, 2.5], 6, 16); w.tab2(1, 6, 16); w.lst(0.0, 1e6, 0, 0, [1e6, 0.0, 0.0, 0.0, 0.0, 0.0], 1, 6, 16)
    w.tab1(float(ZA + 2004 - 2), 37.0, 0, 4, [1e-5, 2e7], [1.0, 1.0], 6, 16)
    w.end(6, 16)
    for mt in (600, 649):
        w.line([float(ZA), 39.62, 0, 2, 2, 0], 6, mt)
        w.tab1(1001.0, 0.99917, 0, 2, [1e-5, 2e7], [1.0, 1.0], 6, mt); w.tab2(2, 6, mt)
        for E in (1e6, 2e7): w.lst(0.0, E, 0, 0, [0.2, 0.05], 2, 6, mt)
        w.tab1(float(ZA + 2004 - 1001), 37.7, 0, 4, [1e-5, 2e7], [1.0, 1.0], 6, mt); w.end(6, mt)
    w.end(6, 0)
    w.out.append('%66s%4d%2d%3d%5d' % ('', 0, 0, 0, 0))
    txt = '\n'.join(w.out) + '\n'
    if path: open(path, 'w').write(txt)
    return txt
def make_jendl(path=None):
    """JENDL-5 style deuteron file for Li-6 (file name d_003-Li-006.dat): MF3 MT2 and MT5 only; MF6/MT5 holds a yield TAB1 for every product, laboratory
    angle-energy distributions (LAW 7) for n, p, a, bare yields (LAW 0) for residual nuclei including the unbound H-4 (ZAP 1004)."""
    ZA = 3006; mat = 325; w = W(mat)
    Es = [0.2e6, 1e6, 5e6, 2e7, 2e8]; E_MeV = [e * 1e-6 for e in Es]
    nonel = [0.3, 0.66, 0.98, 0.8, 0.5]; elastic = [-0.005, 0.2, 0.1, 0.05, 0.02]
    w.out.append('%-66s%4d%2d%3d%5d' % ('   JENDL-5  Li-  6 (synthetic test file)', 0, 0, 0, 0))
    w.line([float(ZA), 5.96345, -1, 0, 6, 1], 1, 451); w.line([0.0, 0.0, 0, 0, 0, 6], 1, 451); w.line([1.9968, 2e8, 0, 0, 10020, 5], 1, 451)
    w.line([0.0, 0.0, 0, 0, 39, 5], 1, 451); w.out.append('%-66s%4d%2d%3d%5d' % ('  3-Li-  6 JAEA       SYNTHETIC', mat, 1, 451, 5)); w.end(1, 451); w.end(1, 0)
    for mt, ys in [(2, elastic), (5, nonel)]:
        w.line([float(ZA), 5.96345, 0, 0, 0, 0], 3, mt); w.tab1(0.0, 0.0, 0, 0, Es, ys, 3, mt); w.end(3, mt)
    w.end(3, 0)
    w.line([float(ZA), 5.96345, 0, 1, 7, 0], 6, 5)
    for zap, awp, ys, law in [(1, 1.0, [0.25, 0.25, 0.25, 0.45, 0.65], 7), (1, 1.0, [0.15, 0.15, 0.15, 0.27, 0.39], 7), (1, 1.0, [0.1, 0.1, 0.1, 0.18, 0.26], 7), (1001, 0.99917, [0.2, 0.2, 0.3, 0.5, 0.9], 7), (2004, 3.9726, [0.1, 0.1, 0.4, 0.1, 0.1], 7),
                              (3007, 6.956, [0.4, 0.4, 0.2, 0.1, 0.0], 0), (1004, 3.992, [0.0, 0.0, 0.1, 0.1, 0.1], 0)]:
        w.tab1(float(zap), awp, 0, law, Es, ys, 6, 5)
        if law == 7:
            w.tab2(2, 6, 5)
            for E in (1e6, 2e8):
                w.line([0.0, E, 0, 0, 1, 2], 6, 5); w.line([2, 2], 6, 5)
                for mu in (-1.0, 1.0): w.tab1(0.0, mu, 0, 0, [0.0, E * 0.5, E], [0.0, 1.0 / E, 0.0], 6, 5)   # triangular E', isotropic; integral over E' and mu = 1
    w.end(6, 5); w.end(6, 0)
    w.out.append('%66s%4d%2d%3d%5d' % ('', -1, 0, 0, 0))
    txt = '\n'.join(w.out) + '\n'
    if path: open(path, 'w').write(txt)
    return txt
def make_jendl_alpha(path=None):
    """JENDL-5 style alpha file for Al-27 (a_013-Al-027.dat): MT5 total reaction without neutrons in MF6 (p, α as LAW 7, a residual with LAW 0) and the neutron production
    MT201 (MF3 cross section, MF6 neutron with yield 1 and LAW 7)."""
    ZA = 13027; mat = 1325; w = W(mat)
    Es = [0.5e6, 3.0e6, 6.0e6, 10.0e6, 15.0e6]
    s5 = [0.01, 0.2, 0.5, 0.9, 1.1]; s201 = [0.0, 0.0, 0.05, 0.27, 0.34]; sel = [-0.01, -0.05, -0.1, -0.15, -0.06]
    w.out.append('%-66s%4d%2d%3d%5d' % ('   JENDL-5  Al- 27 (synthetic test file)', 0, 0, 0, 0))
    w.line([float(ZA), 26.74975, -1, 0, 6, 1], 1, 451); w.line([0.0, 0.0, 0, 0, 0, 6], 1, 451); w.line([3.968219, 1.5e7, 0, 0, 20040, 5], 1, 451)
    w.line([0.0, 0.0, 0, 0, 3, 5], 1, 451); w.out.append('%-66s%4d%2d%3d%5d' % ('  13-Al- 27 JAEA       SYNTHETIC', mat, 1, 451, 5)); w.end(1, 451); w.end(1, 0)
    for mt, ys, q in [(2, sel, 0.0), (5, s5, 0.0), (201, s201, -2.6363)]:
        w.line([float(ZA), 26.74975, 0, 0, 0, 0], 3, mt); w.tab1(q * 1e6, q * 1e6, 0, 0, Es, ys, 3, mt); w.end(3, mt)
    w.end(3, 0)
    def law7(w, mt):
        w.tab2(2, 6, mt)
        for E in (1e6, 1.5e7):
            w.line([0.0, E, 0, 0, 1, 2], 6, mt); w.line([2, 2], 6, mt)
            for mu in (-1.0, 1.0): w.tab1(0.0, mu, 0, 0, [0.0, E * 0.5, E], [0.0, 1.0 / E, 0.0], 6, mt)
    w.line([float(ZA), 26.74975, 0, 1, 3, 0], 6, 5)
    for zap, awp, ys, law in [(1001, 0.99917, [0.8, 0.7, 0.6, 0.5, 0.5], 7), (2004, 3.9726, [0.0, 0.0, 0.1, 0.2, 0.3], 7), (13027, 26.75, [0.2, 0.3, 0.3, 0.3, 0.2], 0)]:
        w.tab1(float(zap), awp, 0, law, Es, ys, 6, 5)
        if law == 7: law7(w, 5)
    w.end(6, 5)
    w.line([float(ZA), 26.74975, 0, 1, 1, 0], 6, 201); w.tab1(1.0, 1.0, 0, 7, Es, [1.0] * 5, 6, 201); law7(w, 201); w.end(6, 201); w.end(6, 0)
    w.out.append('%66s%4d%2d%3d%5d' % ('', -1, 0, 0, 0))
    txt = '\n'.join(w.out) + '\n'
    if path: open(path, 'w').write(txt)
    return txt
def make_neutron(pathH=None, pathC=None):
    """Synthetic incident-neutron files in ENDF/B-VIII.0 layout: n+1H (MF3 MT2, MF4 MT2 Legendre) and n+12C (MF3 MT2, 4, 51, 52, 91, 102, 103, 107; MF4 MT2 tabulated)."""
    def head(w, ZA, AWR, name):
        w.out.append('%-66s%4d%2d%3d%5d' % (' ' + name, 0, 0, 0, 0))
        w.line([float(ZA), AWR, -1, 0, 0, 0], 1, 451); w.line([0.0, 0.0, 0, 0, 0, 0], 1, 451); w.end(1, 451); w.end(1, 0)
    Es = [1e3 * 1.15 ** k for k in range(0, 90)]; Es = [e for e in Es if e <= 4.0e7]; Em = [e * 1e-6 for e in Es]
    w = W(125); head(w, 1001, 0.99917, 'n+H1 synthetic test')
    sig = []
    for e in Em:
        k2 = e / 82.94; k = math.sqrt(k2)
        cT = (-1 / 5.424 + 0.5 * 1.759 * k2) / k; cS = (1 / 23.74 + 0.5 * 2.77 * k2) / k
        sig.append(4 * math.pi / k2 * (0.75 / (1 + cT * cT) + 0.25 / (1 + cS * cS)) * 0.01)
    w.line([1001.0, 0.99917, 0, 0, 0, 0], 3, 2); w.tab1(0.0, 0.0, 0, 0, Es, sig, 3, 2); w.end(3, 2); w.end(3, 0)
    w.line([1001.0, 0.99917, 0, 1, 0, 0], 4, 2); w.line([0.0, 0.99917, 0, 2, 0, 0], 4, 2); w.tab2(2, 4, 2)
    w.lst(0.0, 1e6, 0, 0, [0.0, 0.0], 0, 4, 2); w.lst(0.0, 4e7, 0, 0, [0.2, 0.0], 0, 4, 2)
    w.end(4, 2); w.end(4, 0); w.out.append('%66s%4d%2d%3d%5d' % ('', -1, 0, 0, 0))
    if pathH: open(pathH, 'w').write('\n'.join(w.out) + '\n')
    w = W(600); head(w, 6012, 11.8969, 'n+C12 synthetic test')
    def sg(thr, peak): return [0.0 if e < thr else peak * (1 - math.exp(-(e - thr) / 3.0)) * math.exp(-(e - thr) / 15.0) for e in Em]
    el = [1.5 + 3.0 * math.exp(-e) for e in Em]; i51 = sg(4.81, 0.3); i52 = sg(8.3, 0.12); i91 = sg(12.5, 0.2); c102 = [0.001 * math.exp(-e) for e in Em]; c103 = sg(13.7, 0.04); c107 = sg(6.2, 0.07)
    i4 = [a + b + c for a, b, c in zip(i51, i52, i91)]
    for mt, ys, q in [(2, el, 0.0), (4, i4, 0.0), (51, i51, -4.439e6), (52, i52, -7.654e6), (91, i91, -11.5e6), (102, c102, 4.9e6), (103, c103, -12.587e6), (107, c107, -5.702e6)]:
        w.line([6012.0, 11.8969, 0, 0, 0, 0], 3, mt); w.tab1(q, q, 0, 0, Es, ys, 3, mt); w.end(3, mt)
    w.end(3, 0)
    w.line([6012.0, 11.8969, 0, 2, 0, 0], 4, 2); w.line([0.0, 11.8969, 0, 2, 0, 0], 4, 2); w.tab2(2, 4, 2)
    for E in (1e6, 4e7): w.tab1(0.0, E, 0, 0, [-1.0, 0.0, 1.0], [0.1, 0.4, 1.5] if E > 1e6 else [0.5, 0.5, 0.5], 4, 2)
    w.end(4, 2); w.end(4, 0); w.out.append('%66s%4d%2d%3d%5d' % ('', -1, 0, 0, 0))
    if pathC: open(pathC, 'w').write('\n'.join(w.out) + '\n')
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'neutron': make_neutron('tests/data/n_001-H-001.endf', 'tests/data/n_006-C-012.endf'); sys.exit(0)
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'jendl': make_jendl(sys.argv[2] if len(sys.argv) > 2 else 'tests/data/d_003-Li-006.dat'); make_jendl_alpha(sys.argv[3] if len(sys.argv) > 3 else 'tests/data/a_013-Al-027.dat'); sys.exit(0)
if __name__ == '__main__': make(path=sys.argv[1] if len(sys.argv) > 1 else 'tests/data/a-Ar040.tendl', scale=float(sys.argv[2]) if len(sys.argv) > 2 else 1.0)
