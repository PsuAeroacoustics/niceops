import numpy as np
from blend import Conditions
from direction_stats import random_states

z = np.load('B407_M2.npz')
L = z['grid']; C = Conditions(z['cond'])
rng = np.random.default_rng(1)
q = random_states(C, 20000, rng)
w, idx, x = C.barycentric(q)
recon = np.einsum('nk,nkd->nd', w, C.norm[idx])
e = np.abs(recon - x).max(1)
print('states with recon err > 1e-9:', (e > 1e-9).sum(), ' max', e.max())
s = C.tri.find_simplex(x)
vol = []
for k in np.nonzero(e > 1e-9)[0][:5]:
    v = C.verts[C.tri.simplices[s[k]]]
    vol.append(abs(np.linalg.det(v[1:] - v[0])) / 6)
    print(' state', q[k], 'w', np.round(w[k], 4), 'vol', vol[-1])
# distribution of simplex volumes overall
V = C.verts[C.tri.simplices]
vols = np.abs(np.linalg.det(V[:, 1:] - V[:, :1])) / 6
print('zero-volume simplices (<1e-14):', (vols < 1e-14).sum(), 'of', len(vols))
print('states landing in zero-volume simplices:', np.sum(vols[s] < 1e-14))
# low-level spheres
Lk = np.moveaxis(L, -1, 0)
mn = Lk.reshape(len(Lk), -1).min(1)
print('conditions with min level < 50 dB:', np.sum(mn < 50), C.raw[mn < 50][:10], mn[mn < 50][:10])
