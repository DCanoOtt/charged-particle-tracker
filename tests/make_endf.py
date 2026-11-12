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
if __name__ == '__main__': make(path=sys.argv[1] if len(sys.argv) > 1 else 'tests/data/a-Ar040.tendl', scale=float(sys.argv[2]) if len(sys.argv) > 2 else 1.0)
