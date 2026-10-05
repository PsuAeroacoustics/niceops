# Size of a 2*roll phi error: A-weighted level change per degree of phi between filled neighbouring bins.
import numpy as np, sys
D = "/Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/"
def aw(f):
    f2 = f * f
    ra = 12194**2 * f2**2 / ((f2 + 20.6**2) * np.sqrt((f2 + 107.7**2) * (f2 + 737.9**2)) * (f2 + 12194**2))
    return 20 * np.log10(ra) + 2.0
for name in sys.argv[1:]:
    L = open(D + name, errors='replace').read().splitlines()
    i = next(k for k, l in enumerate(L) if l.startswith('NFREQ'))
    freqs = np.array([float(x) for x in L[i + 1].split()])
    grid = {}; phi = None
    for l in L[i + 2:]:
        if l.strip().startswith('PHIOBSAC'):
            phi = float(l.split('=')[1]); continue
        t = l.split()
        if len(t) == len(freqs) + 1 and phi is not None:
            v = np.array([float(x) for x in t[1:]])
            if np.all(v > -998):
                grid[(phi, float(t[0]))] = 10 * np.log10(np.sum(10 ** ((v + aw(freqs)) / 10)))
    g = []
    for (p, th), la in grid.items():
        q = grid.get((p + 10, th))
        if q is not None: g.append(abs(q - la) / 10)
    g = np.array(g)
    print(f"{name}: filled bins {len(grid)}, |dLA/dphi| deg^-1: median {np.median(g):.3f}, 90th pct {np.percentile(g, 90):.3f}, max {g.max():.3f}"
          f" -> 5.6 deg error: median {5.6*np.median(g):.2f} dB, 90th pct {5.6*np.percentile(g, 90):.2f} dB")
