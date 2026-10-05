# Convention-free test: build the hemisphere frame banked toward (or away from) the turn centre,
# with the centre taken from the trajectory's own curvature, and see which reproduces PHIOBSAC.
import numpy as np, sys
sys.path.insert(0, '/Users/evg5332/.claude/jobs/a090ec0c/tmp/rollcheck')
from fit_roll import load_inp, load_his, A, wrap

for arg in sys.argv[1:]:
    case, name = arg.split(':')
    I = load_inp(A + case + "/SE/" + name + ".inp")
    M = load_his(A + case + "/SE/" + name + ".his")
    T = sorted(I); P = np.array([[I[t]['XAC_M'], I[t]['YAC_M']] for t in T])
    idx = {t: i for i, t in enumerate(T)}
    res = {'in': [], 'out': []}; summary = {}
    mm = M[sorted(M)[0]]; o = mm['obs']; O = np.array([o[0], o[1], o[2] + o[3]])
    for r in mm['rows']:
        t = round(r[3], 3); st = I.get(t)
        if st is None or abs(st['ROLL']) < 1e-6: continue
        i = idx[t]
        if i == 0 or i == len(T) - 1: continue
        a = P[i+1] - 2 * P[i] + P[i-1]            # centripetal acceleration direction (horizontal)
        if np.hypot(*a) < 1e-6: continue
        c = np.array([a[0], a[1], 0.]) / np.hypot(*a)
        v = P[i+1] - P[i-1]
        turn = 'R' if (v[0] * a[1] - v[1] * a[0]) < 0 else 'L'   # cross(v, a) < 0 in x-east/y-north => clockwise
        h = np.radians(st['HEADING']); f = np.array([np.sin(h), np.cos(h), 0.])
        mag = np.radians(abs(st['ROLL']))
        d = O - np.array(r[:3]); d /= np.linalg.norm(d)
        for k, sg in (('in', +1), ('out', -1)):
            up = np.cos(mag) * np.array([0, 0, 1.]) + sg * np.sin(mag) * c
            dn = -up; s = np.cross(dn, f)
            ph = np.degrees(np.arctan2(d @ s, d @ dn))
            res[k].append(wrap(ph - r[7]) if (abs(r[7]) < 89.99 and 20 < r[6] < 160) else np.nan)
        key = (turn, st['ROLL'])
        summary[key] = summary.get(key, 0) + 1
    for k in res:
        e = np.array(res[k]); e = e[~np.isnan(e)]
        print(f"{case} {name}: bank {k:3s} of turn: n={len(e)} rms dphi={np.sqrt((e**2).mean()):.4f} max={np.abs(e).max():.3f}")
    print('   (turn sense from curvature, ROLL) -> samples:', summary)
