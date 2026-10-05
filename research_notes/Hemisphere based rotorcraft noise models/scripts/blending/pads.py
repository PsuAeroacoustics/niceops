import numpy as np
for db in ['S-76D_M3', 'AW139_M1', 'B205_M1', 'B407_M2']:
    z = np.load(db + '.npz'); c = z['cond']; L = np.moveaxis(z['grid'], -1, 0).reshape(len(c), -1)
    ctr = np.unique(c[:, 2]); ctr = ctr[ctr > 0].min()
    ref = np.isclose(c[:, 2], ctr, atol=1e-12, rtol=0)
    meas = ref & (c[:, 0] > 0) & (c[:, 1] > -23.9) & (c[:, 1] < 34.9)
    pad = ref & ~meas
    copies = [np.any(np.all(np.abs(L[meas] - L[i]) < 1e-9, 1)) for i in np.nonzero(pad)[0]]
    hov = ref & (c[:, 0] == 0)
    hcopy = [np.any(np.all(np.abs(L[meas] - L[i]) < 1e-9, 1)) for i in np.nonzero(hov)[0]]
    print(db, 'FCs', ref.sum(), 'measured', meas.sum(), 'padding', pad.sum(), 'of which exact copies of a measured sphere', sum(copies),
          '; mu=0 entries', hov.sum(), 'copies', sum(hcopy),
          '; measured gamma range %.1f..%.1f, mu %.3f..%.3f' % (c[meas, 1].min(), c[meas, 1].max(), c[meas, 0].min(), c[meas, 0].max()))
