import numpy as np, sys
sys.path.insert(0, '/Users/evg5332/.claude/jobs/a090ec0c/tmp/rollcheck')
from fit_roll import *
case, name = sys.argv[1].split(':')
I = load_inp(A + case + "/SE/" + name + ".inp")
M = load_his(A + case + "/SE/" + name + ".his")
mm = M[1]; o = mm['obs']; O = np.array([o[0], o[1], o[2] + o[3]])
print('observer', O)
prev = None
for r in mm['rows']:
    st = I.get(round(r[3], 3))
    if st is None: continue
    d = O - np.array(r[:3]); d /= np.linalg.norm(d)
    th, ph = angles(d, np.radians(st['HEADING']), np.radians(st['ROLL']), +1)
    e = wrap(ph - r[7])
    if abs(e) > 0.05 or (st['ROLL'] != 0 and r[3] % 20 == 0):
        print(f"t={r[3]:7.2f} HEAD={st['HEADING']:6.1f} ROLL={st['ROLL']:5.2f} Vang={st['Vang']:5.2f} mode={st['iMode']:.0f} "
              f"theta_his={r[6]:6.2f} phi_his={r[7]:7.2f} phi_fit={ph:7.2f} err={e:6.3f}")
