# Independent re-derivation: which sense of ROLL reproduces SingleEvt's THETAOBS/PHIOBSAC?
import numpy as np, sys
A = "/Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/"

def load_inp(p):
    hdr = None; D = {}
    for l in open(p):
        t = l.split()
        if t and t[0] == 'TEM': hdr = t; continue
        if hdr and len(t) == len(hdr):
            try: v = [float(x) for x in t]
            except ValueError: continue
            D[round(v[0], 3)] = dict(zip(hdr, v))
    return D

def load_his(p):
    mics = {}; m = None
    for l in open(p):
        t = l.split()
        if len(t) == 2 and t[0] == 'Micro#':
            m = int(t[1]); mics[m] = {'obs': None, 'rows': []}; continue
        if m is None: continue
        try: v = [float(x) for x in t]
        except ValueError: continue
        if len(v) == 13 and mics[m]['obs'] is None: mics[m]['obs'] = v
        elif len(v) == 12: mics[m]['rows'].append(v)
    return mics

def frame(h, roll, sgn):
    # x east, y north, z up.  h = rotation angle clockwise from north.
    f = np.array([np.sin(h), np.cos(h), 0.])
    s = np.array([np.cos(h), -np.sin(h), 0.])   # starboard (right of nose)
    dn = np.array([0, 0, -1.])
    r = sgn * roll                               # sgn=+1: +ROLL puts starboard wing down
    s2 = np.cos(r) * s + np.sin(r) * dn
    dn2 = np.cos(r) * dn - np.sin(r) * s
    return f, s2, dn2

def angles(d, h, roll, sgn):
    f, s, dn = frame(h, roll, sgn)
    th = np.degrees(np.arccos(np.clip(d @ f, -1, 1)))
    ph = np.degrees(np.arctan2(d @ s, d @ dn))
    return th, ph

def wrap(a): return (a + 180) % 360 - 180

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        case, name = arg.split(':')
        I = load_inp(A + case + "/SE/" + name + ".inp")
        M = load_his(A + case + "/SE/" + name + ".his")
        stats = {+1: [], -1: [], 0: []}
        for m, mm in M.items():
            o = mm['obs']; O = np.array([o[0], o[1], o[2] + o[3]])
            for r in mm['rows']:
                st = I.get(round(r[3], 3))
                if st is None or abs(st["ROLL"]) > 1e-6: continue
                d = O - np.array(r[:3]); d /= np.linalg.norm(d)
                for sg in (+1, -1, 0):
                    th, ph = angles(d, np.radians(st['HEADING']), np.radians(st['ROLL']), sg)
                    stats[sg].append([th - r[6], wrap(ph - r[7])])
        for sg, e in stats.items():
            e = np.array(e)
            lab = {1: '+ROLL = starboard down', -1: '+ROLL = port down     ', 0: 'roll ignored          '}[sg]
            print(f"{case} {name} {lab}: n={len(e)} mics={len(M)} rms dth={np.sqrt((e[:,0]**2).mean()):.4f} "
                  f"rms dphi={np.sqrt((e[:,1]**2).mean()):.4f} max|dphi|={np.abs(e[:,1]).max():.3f}")
