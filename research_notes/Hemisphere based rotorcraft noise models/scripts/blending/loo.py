"""Experiment 3: leave-one-out.  Hold out one measured flight condition (all of
its CT copies, and any padding copies of the same sphere), rebuild the
triangulation from the rest, predict the held-out sphere at its own
(mu, gamma, CT_ref) by each blend, and compare with the stored sphere over the
footprint directions.  Errors include the held-out sphere's own measurement
scatter, so they are an upper bound on interpolation error, and the hull
shrinks where a held-out point was on it."""
import json
import numpy as np
from blend import Conditions, idw_weights

out = {}
for db in ['S-76D_M3', 'AW139_M1', 'B205_M1', 'B407_M2']:
    z = np.load(db + '.npz'); L = z['grid']; cond = z['cond']; phis = z['phis']; thetas = z['thetas']
    P, Th = np.meshgrid(phis, thetas, indexing='ij')
    fmask = ((np.abs(P) <= 60) & (Th >= 20) & (Th <= 160)).ravel()
    Lk = np.moveaxis(L, -1, 0).reshape(len(cond), -1)[:, fmask]
    C0 = Conditions(cond)
    ref = np.isclose(cond[:, 2], C0.ct_ref, rtol=0, atol=1e-12)
    meas = ref & (cond[:, 0] > 0) & (cond[:, 1] > -23.9) & (cond[:, 1] < 34.9)
    err = {k: [] for k in ['dB', 'energy', 'NORAH_IDW_energy', 'nearest3D', 'nearestFC']}
    inside = 0
    for h in np.nonzero(meas)[0]:
        same = np.all(np.abs(Lk - Lk[h]) < 1e-9, axis=1) | np.all(np.abs(Lk - Lk[h] - (Lk[:, :1] - Lk[h, :1])) < 1e-9, axis=1)
        # same sphere up to the CT offset: drop it and all its copies
        drop = same | np.all(np.isclose(cond[:, :2], cond[h, :2]), axis=1)
        keep = ~drop
        # the reduced database must keep the original normalisation box to be comparable;
        # Conditions normalises by its own min/max, which is what NICEOPS would do.
        C = Conditions(cond[keep])
        q = cond[h:h + 1]
        x = C.normalize(q)
        if C.tri.find_simplex(x)[0] < 0 or np.any((q - C.lo) / (C.hi - C.lo) < 0) or np.any((q - C.lo) / (C.hi - C.lo) > 1):
            continue  # held-out point outside the reduced hull: extrapolation, skip
        inside += 1
        w, idx, x = C.barycentric(q)
        Lr = Lk[keep]
        truth = Lk[h]
        a = w[0] @ Lr[idx[0]]
        b = 10 * np.log10(w[0] @ 10 ** (0.1 * Lr[idx[0]]))
        wi = idw_weights(C, x, w, idx)
        d = 10 * np.log10(wi[0] @ 10 ** (0.1 * Lr[idx[0]]))
        c = Lr[C.nearest(q)[0]]
        nfc, corr = C.nearest_flight_condition(q)
        c2 = Lr[nfc[0]] + corr[0]
        for k, v in zip(err, [a, b, d, c, c2]):
            err[k].append(v - truth)
    res = {}
    for k, v in err.items():
        e = np.array(v)
        per = np.sqrt(np.mean(e ** 2, 1))
        res[k] = dict(n_held_out=int(len(e)), bias=float(e.mean()), rms=float(np.sqrt(np.mean(e ** 2))),
                      median_abs=float(np.median(np.abs(e))), p95_abs=float(np.percentile(np.abs(e), 95)),
                      max_abs=float(np.abs(e).max()), median_sphere_rms=float(np.median(per)))
    # paired comparison: in how many held-out spheres is energy better than dB (rms)?
    ea = np.sqrt(np.mean(np.array(err['dB']) ** 2, 1)); eb = np.sqrt(np.mean(np.array(err['energy']) ** 2, 1))
    en = np.sqrt(np.mean(np.array(err['nearest3D']) ** 2, 1))
    res['energy_better_than_dB_fraction'] = float(np.mean(eb < ea))
    res['dB_better_than_nearest_fraction'] = float(np.mean(ea < en))
    res['measured_flight_conditions'] = int(meas.sum())
    out[db] = res
    print('==', db, 'measured FCs', meas.sum(), 'held out inside reduced hull', inside)
    for k, v in res.items():
        print('  ', k, v if not isinstance(v, dict) else
              ' '.join('%s=%.3f' % (kk, vv) if isinstance(vv, float) else '%s=%d' % (kk, vv) for kk, vv in v.items()))
json.dump(out, open('loo.json', 'w'), indent=1)
