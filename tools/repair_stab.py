#!/usr/bin/env python3
"""Repair of known faults in the ICRU 73 heavy-ion tables (STAB.h) inherited from the libdedx embedded data.

Faults found by audit (tests/test_stab.js re-checks all of them on every run):
 1. Rows Z1 = 11..18 of targets 18 (Ar), 215 (photographic emulsion) and 277 (water vapour) are shifted by one
    energy bin: element 0 is a leaked value (0, Z1*1000, or the first value of the Si row) and element i holds the
    stopping power of hE[i-1]; the 1000 MeV/u value is lost.  Fix: drop element 0, shift left, and rebuild the last
    point as S(900) * S10(1000)/S10(900) with S10 the (correct) Z1 = 10 row of the same target
    (validated on healthy targets: error < 0.03 %).
 2. Rows Z1 = 16 of targets 215 and 277 contain the Z1 = 8 table (shifted, header 8000).  Rebuilt from the Z1 = 15 and 17 rows
    of the same target and the Z1 = 16 row of Si (target 14) (sister-ratio interpolation in Z1, validated on healthy
    targets: mean error 0.05-0.3 %).
 3. Duplicated rows from another Z1: (16, Ti=22) = (14, Ti); (9, target 209) = (17, 209); (3, NaI=252) = (9, 252).
    Rebuilt with the same sister-ratio method (the Z1 = 3 row is extrapolated from Z1 = 4, 5).
Reconstructed values are rounded to 4 significant digits like the rest of the table.
"""
import json, math, re, sys
PATH = sys.argv[1] if len(sys.argv) > 1 else 'data/stopping-tables.js'
txt = open(PATH).read()
m = re.search(r'^const STAB = (.*);\s*$', txt, re.M)
S = json.loads(m.group(1)); H = S['h']; L = math.log
def r4(x): return float('%.4g' % x)
def get(Z, i): return H[str(Z)][str(i)]
done = []
# 1. shifted rows
SHIFT = [(Z, i) for Z in range(11, 19) for i in (18, 215, 277) if not (Z == 16 and i in (215, 277))]
for Z, i in SHIFT:
    v = get(Z, i)
    if v[0] in (0, Z * 1000) or v[0] == get(Z, 14)[0]:               # leaked first element (not yet repaired)
        v10 = get(10, i)
        new = v[1:] + [r4(v[52] * v10[52] / v10[51])]
        H[str(Z)][str(i)] = new; done.append('shift %d/%d' % (Z, i))
# 2/3. rebuilt rows: (Z, id, neighbour Zs)
def sister(Z, i, za, zb=None, zc=None):
    ra = [L(a / b) for a, b in zip(get(za, i), get(za, 14))]
    if zb is None: r = ra
    elif zc is None:
        rb = [L(a / b) for a, b in zip(get(zb, i), get(zb, 14))]; r = [(x + y) / 2 for x, y in zip(ra, rb)]
    else:                                                              # linear extrapolation from two rows
        rb = [L(a / b) for a, b in zip(get(zb, i), get(zb, 14))]; r = [2 * x - y for x, y in zip(ra, rb)]
    return [r4(math.exp(x) * s) for x, s in zip(r, get(Z, 14))]
REB = [(16, 215, 15, 17, None), (16, 277, 15, 17, None), (16, 22, 15, 17, None), (9, 209, 8, 10, None), (3, 252, 4, 5, 'x')]
def bad(Z, i):                                                         # still has the diagnosed fault?
    v = get(Z, i)
    if (Z, i) in ((16, 215), (16, 277)): return v[0] == 8000 or get(8, i)[:5] == v[1:6]
    if (Z, i) == (16, 22): return v == get(14, 22)
    if (Z, i) == (9, 209): return v == get(17, 209)
    if (Z, i) == (3, 252): return v == get(9, 252)
for Z, i, a, b, ex in REB:
    if bad(Z, i):
        H[str(Z)][str(i)] = sister(Z, i, a, b, ex); done.append('rebuild %d/%d' % (Z, i))
print(len(done), 'rows changed:', ', '.join(done))
if done:
    new = json.dumps(S, separators=(',', ':'))
    note = '// Repairs applied by tools/repair_stab.py (faults inherited from libdedx\'s ICRU 73 data): see that file for the list.\n'
    out = txt[:m.start(1)] + new + txt[m.end(1):]
    if 'tools/repair_stab.py' not in out: out = out.replace('const STAB = ', note + 'const STAB = ', 1)
    open(PATH, 'w').write(out)
