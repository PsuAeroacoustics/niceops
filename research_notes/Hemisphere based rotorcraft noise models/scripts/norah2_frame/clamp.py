# Does SingleEvt clamp PHIOBSAC at +-90 for directions above the (rolled) hemisphere's equator?
import numpy as np, sys
sys.path.insert(0, '/Users/evg5332/.claude/jobs/a090ec0c/tmp/rollcheck')
from fit_roll import load_inp, load_his, angles, A, wrap
for arg in sys.argv[1:]:
    case, name = arg.split(':')
    I = load_inp(A + case + "/SE/" + name + ".inp"); M = load_his(A + case + "/SE/" + name + ".his")
    n_eq = 0; n_above = 0; agree = 0; worst_unclamped = 0
    for m, mm in M.items():
        o = mm['obs']; O = np.array([o[0], o[1], o[2] + o[3]])
        for r in mm['rows']:
            st = I.get(round(r[3], 3))
            if st is None: continue
            d = O - np.array(r[:3]); d /= np.linalg.norm(d)
            th, ph = angles(d, np.radians(st['HEADING']), np.radians(st['ROLL']), +1)
            if abs(ph) > 90:
                n_above += 1
                if abs(abs(r[7]) - 90) < 0.006 and np.sign(r[7]) == np.sign(ph) and abs(th - r[6]) < 0.01: agree += 1
            else:
                worst_unclamped = max(worst_unclamped, abs(wrap(ph - r[7])))
            if abs(abs(r[7]) - 90) < 0.006: n_eq += 1
    print(f"{case} {name}: rows with |phi_geom|>90: {n_above}, of which his reports phi=+-90.00 with theta unchanged: {agree}; "
          f"rows his |phi|=90: {n_eq}; max |dphi| on the rest: {worst_unclamped:.3f}")
