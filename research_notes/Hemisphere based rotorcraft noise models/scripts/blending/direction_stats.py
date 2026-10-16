"""Experiment 1: per (state, sphere direction) differences between condition
blends, at the sphere grid nodes (no sphere-grid interpolation involved).

  a  barycentric dB blend (NICEOPS)
  b  barycentric energy blend 10 log10 sum w 10^(L/10)
  c  nearest stored condition, normalized (mu, gamma, CT) Euclidean
  c2 nearest (mu, gamma) flight condition + exact CT law of the database
  t  dB blend on a second valid Delaunay triangulation (Qhull QJ), vs a
"""
import glob
import json
import sys
import numpy as np
from blend import Conditions, idw_weights

STRIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 5
SETS = [('S-76D_M3', 'S-76D'), ('AW139_M1', 'AW139'), ('B205_M1', 'B205'), ('B407_M2', None)]
ST = '/Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend/states/'


def load_states(vehicle):
    qs = []
    for f in sorted(glob.glob(ST + vehicle + '/state_*.csv')):
        a = np.genfromtxt(f, delimiter=',', names=True)
        qs.append(np.c_[a['mu'], a['fpa_e'], a['CT']][::STRIDE])
    return np.vstack(qs)


def random_states(C, n, rng):
    """Uniform in the box of the measured (non-pad) flight conditions, LF 0.9-2."""
    c = C.raw
    meas = (c[:, 0] > 0) & (c[:, 1] > -23.9) & (c[:, 1] < 34.9)
    mlo = c[meas, :2].min(0); mhi = c[meas, :2].max(0)
    q = np.empty((n, 3))
    q[:, 0] = rng.uniform(mlo[0], mhi[0], n)
    q[:, 1] = rng.uniform(mlo[1], mhi[1], n)
    q[:, 2] = C.ct_ref * rng.uniform(0.9, 2.0, n)
    return q


def pct(x):
    x = x[np.isfinite(x)]
    return dict(n=int(x.size), mean=float(x.mean()), median=float(np.median(x)),
                p95=float(np.percentile(x, 95)), p99=float(np.percentile(x, 99)),
                max=float(x.max()), min=float(x.min()),
                absp95=float(np.percentile(np.abs(x), 95)), absmax=float(np.abs(x).max()))


def run(dbname, vehicle, rng):
    z = np.load(dbname + '.npz')
    L = z['grid']; phis = z['phis']; thetas = z['thetas']
    C = Conditions(z['cond'])
    CJ = Conditions(z['cond'], qhull_options='QJ Qbb Qc Q12')
    print(dbname, 'scipy simplices', len(C.tri.simplices), 'ext vertices', len(C.verts),
          ' QJ simplices', len(CJ.tri.simplices), flush=True)
    if vehicle:
        q = load_states(vehicle); source = 'trajectories (%s, every %d-th state)' % (vehicle, STRIDE)
    else:
        q = random_states(C, 20000, rng); source = 'random sample'
    w, idx, x = C.barycentric(q)
    wj, idxj, _ = CJ.barycentric(q)
    near = C.nearest(q)
    nfc, corr = C.nearest_flight_condition(q)
    # quality checks on the barycentric weights
    recon = np.einsum('nk,nkd->nd', w, C.norm[idx])
    err = np.abs(recon - x).max()
    single = np.mean(w[:, 0] == 1.0)
    print('  NaN weights:', np.isnan(w).any(1).sum(), np.isnan(wj).any(1).sum(), ' min weight', np.nanmin(w), np.nanmin(wj), flush=True)
    # direction masks
    P, Th = np.meshgrid(phis, thetas, indexing='ij')
    masks = {'all stored directions': np.ones_like(P, bool),
             'lower hemisphere |phi|<=90 (all lookup_angles can reach)': np.abs(P) <= 90,
             'footprint |phi|<=60, 20<=theta<=160': (np.abs(P) <= 60) & (Th >= 20) & (Th <= 160)}
    # level arrays (n_state, n_phi, n_theta)
    Lk = np.moveaxis(L, -1, 0)  # (cond, phi, theta)
    a = np.einsum('nk,nkpt->npt', w, Lk[idx])
    b = 10 * np.log10(np.einsum('nk,nkpt->npt', w, 10 ** (0.1 * Lk[idx])))
    c = Lk[near]
    c2 = Lk[nfc] + corr[:, None, None]
    t = np.einsum('nk,nkpt->npt', wj, Lk[idxj])
    wi = idw_weights(C, x, w, idx)
    dd = 10 * np.log10(np.einsum('nk,nkpt->npt', wi, 10 ** (0.1 * Lk[idx])))
    # spread of the four blended levels (max - min over conditions with weight)
    Lsel = Lk[idx]
    wpos = (w > 1e-12)[:, :, None, None]
    spread = np.where(wpos, Lsel, -np.inf).max(1) - np.where(wpos, Lsel, np.inf).min(1)
    res = dict(db=dbname, source=source, n_states=int(len(q)), weight_reconstruction_err=float(err),
               single_vertex_fraction=float(single), scipy_simplices=int(len(C.tri.simplices)),
               qhull_qj_simplices=int(len(CJ.tri.simplices)),
               states_outside_hull=float(np.mean(C.tri.find_simplex(x) < 0)),
               query_clamped_fraction=float(np.mean(np.any((q - C.lo) / (C.hi - C.lo) < 0, 1) | np.any((q - C.lo) / (C.hi - C.lo) > 1, 1))),
               stats={})
    for mname, m in masks.items():
        res['stats'][mname] = {
            'energy - dB (b-a)': pct((b - a)[:, m].ravel()),
            'nearest3D - dB (c-a)': pct((c - a)[:, m].ravel()),
            'nearestFC+CTlaw - dB (c2-a)': pct((c2 - a)[:, m].ravel()),
            'other triangulation - dB (t-a)': pct((t - a)[:, m].ravel()),
            'NORAH IDW energy - dB (d-a)': pct((dd - a)[:, m].ravel()),
            'spread of blended levels': pct(spread[:, m].ravel()),
        }
    # per-state maxima over footprint directions: how bad is the worst direction
    m = masks['footprint |phi|<=60, 20<=theta<=160']
    res['per_state_max_over_footprint'] = {
        'b-a': pct((b - a)[:, m].max(1)),
        '|c-a|': pct(np.abs(c - a)[:, m].max(1)),
        '|c2-a|': pct(np.abs(c2 - a)[:, m].max(1)),
        '|t-a|': pct(np.abs(t - a)[:, m].max(1)),
        '|d-a|': pct(np.abs(dd - a)[:, m].max(1)),
    }
    # The energy-dB gap against the spread of the blended levels (Jensen):
    sp = spread[:, m].ravel(); g = (b - a)[:, m].ravel()
    bins = [0, 1, 2, 3, 5, 8, 12, 20, 100]
    res['gap_vs_spread'] = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        k = (sp >= lo) & (sp < hi)
        if k.sum():
            res['gap_vs_spread'].append(dict(spread=f'{lo}-{hi}', frac=float(k.mean()),
                                             median=float(np.median(g[k])),
                                             p95=float(np.percentile(g[k], 95)), max=float(g[k].max())))
    return res


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    out = [run(db, v, rng) for db, v in SETS]
    json.dump(out, open('direction_stats.json', 'w'), indent=1)
    for r in out:
        print('\n==', r['db'], r['source'], 'states', r['n_states'], 'single-vertex', round(r['single_vertex_fraction'], 4),
              'recon err', r['weight_reconstruction_err'], 'clamped', round(r['query_clamped_fraction'], 4))
        for mname, s in r['stats'].items():
            print(' --', mname)
            for k, v in s.items():
                print('   %-34s med %7.3f mean %7.3f p95 %7.3f p99 %7.3f max %7.3f min %7.3f |p95| %6.3f' % (
                    k, v['median'], v['mean'], v['p95'], v['p99'], v['max'], v['min'], v['absp95']))
        print(' -- per-state max over footprint')
        for k, v in r['per_state_max_over_footprint'].items():
            print('   %-8s med %7.3f p95 %7.3f max %7.3f' % (k, v['median'], v['p95'], v['max']))
        print(' -- energy-dB gap vs spread')
        for g in r['gap_vs_spread']:
            print('   ', g)
