import numpy as np
for f in ['S-76D_M3', 'B407_M2']:
    z = np.load(f + '.npz'); g = z['grid']; c = z['cond']
    fc = np.round(c[:, :2], 6)
    u, inv = np.unique(fc, axis=0, return_inverse=True)
    inv = inv.ravel()
    print(f, 'flight conditions', len(u), 'hover-like mu=0:', np.sum(u[:, 0] == 0), u[u[:, 0] == 0][:5])
    for k in [0, 5, 50, 100]:
        idx = np.nonzero(inv == k)[0]
        o = np.argsort(c[idx, 2]); idx = idx[o]
        ref = idx[np.argmin(np.abs(c[idx, 2] - c[idx, 2][1]))]
        print(' fc', u[k], 'CT', np.round(c[idx, 2], 5))
        # mean / spread over directions of level relative to CT index 1
        d = g[:, :, idx] - g[:, :, idx[1]][:, :, None]
        print('   mean dL vs CT[1]', np.round(d.reshape(-1, len(idx)).mean(0), 2),
              ' std', np.round(d.reshape(-1, len(idx)).std(0), 3))
        print('   20log10(CT/CT1)', np.round(20 * np.log10(np.maximum(c[idx, 2], 1e-12) / c[idx[1], 2]), 2))
