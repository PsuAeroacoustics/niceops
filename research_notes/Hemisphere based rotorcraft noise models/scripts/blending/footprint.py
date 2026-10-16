"""Experiment 2: SEL footprint of a trajectory under each condition blend.

Re-implements the broadband default path of NICEOPS (src/sphere_database.cpp
Database::eval, src/util/sel.hpp) in numpy: sphere frame, lookup_angles,
condition blend at the sphere grid nodes, bilinear on (phi, theta), spherical
spreading, EAA * (r - r_s)/1000, 10 dB-down window, reception-time weight.
Variant a (dB blend) is checked against the NICEOPS output on the same grid.
Only the dBA blend changes between variants; EAA is always dB-weight blended.

usage: footprint.py state.csv niceops_lae.csv out.npz
"""
import sys
import numpy as np
import netCDF4
from blend import Conditions, idw_weights

DBF = '/Users/evg5332/Code/D/NICEOPS/.claude/worktrees/hemisphere-review/databases/S-76D_M3.nod'
FT = 0.3048
C_SOUND = 343.2 * np.sqrt(300.71 / 293.15)  # Atmosphere::sound_speed, S-76D.json temperature
T = 0.02

state_file, lae_file, out = sys.argv[1:4]
dbname = sys.argv[4] if len(sys.argv) > 4 else 'S-76D_M3'
z = np.load(dbname + '.npz')
L = z['grid']; E = z['eaa']; phis = z['phis']; thetas = z['thetas']
ds = netCDF4.Dataset(DBF.replace('S-76D_M3', dbname))
g0 = ds.groups['sphere0']
rs = float(g0['radius'][0]) * float(g0['rotor_scale'][0])
ds.close()
C = Conditions(z['cond'])

st = np.genfromtxt(state_file, delimiter=',', names=True)
q = np.c_[st['mu'], st['fpa_e'], st['CT']]
w, idx, _ = C.barycentric(q)
near = C.nearest(q)
nfc, corr = C.nearest_flight_condition(q)

Lk = np.moveaxis(L, -1, 0); Ek = np.moveaxis(E, -1, 0)
S = {}
S['a'] = np.einsum('nk,nkpt->npt', w, Lk[idx])
S['b'] = 10 * np.log10(np.einsum('nk,nkpt->npt', w, 10 ** (0.1 * Lk[idx])))
S['c'] = Lk[near]
S['c2'] = Lk[nfc] + corr[:, None, None]
wi = idw_weights(C, C.normalize(q), w, idx)
S['d'] = 10 * np.log10(np.einsum('nk,nkpt->npt', wi, 10 ** (0.1 * Lk[idx])))
SE = np.einsum('nk,nkpt->npt', w, Ek[idx])
for k in S:
    S[k] = S[k].reshape(len(q), -1)
SE = SE.reshape(len(q), -1)
nth = len(thetas)

with open(lae_file) as f:
    f.readline()
    xs = np.array(f.readline().split(','), float) * FT
    ys = np.array(f.readline().split(','), float) * FT
    ref = np.array([list(map(float, l.split(','))) for l in f if l.strip()])
X, Y = np.meshgrid(xs, ys)  # rows = y, cols = x
obs = np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]

pos = np.c_[st['x'], st['y'], st['z']]
vel = np.c_[st['vx'], st['vy'], st['vz']]
nrm = np.c_[st['nx'], st['ny'], st['nz']]
nrm /= np.linalg.norm(nrm, axis=1)[:, None]
v = vel / np.linalg.norm(vel, axis=1)[:, None]
nose = v - np.sum(v * nrm, 1)[:, None] * nrm
nose /= np.linalg.norm(nose, axis=1)[:, None]
stb = np.cross(nose, nrm)

dphi = phis[1] - phis[0]; dth = thetas[1] - thetas[0]
nst, nobs = len(q), len(obs)
lev = {k: np.empty((nst, nobs)) for k in S}
rng = np.empty((nst, nobs))
for i in range(nst):
    d = obs - pos[i]
    r = np.sqrt(np.sum(d * d, 1))
    along = d @ nose[i]; up = d @ nrm[i]; across = d @ stb[i]
    th = np.degrees(np.arccos(np.clip(along / r, -1, 1)))
    with np.errstate(divide='ignore', invalid='ignore'):
        ph = np.degrees(np.arctan(across / np.abs(up)))
    ph = np.where((up == 0) & (across == 0), 90.0, ph)
    fp = (ph - phis[0]) / dphi
    pl = np.clip(np.floor(fp).astype(int), 0, len(phis) - 2); a_p = fp - pl
    ft = (th - thetas[0]) / dth
    tl = np.clip(np.floor(ft).astype(int), 0, nth - 2); a_t = ft - tl
    c00 = pl * nth + tl; c01 = c00 + 1; c10 = c00 + nth; c11 = c10 + 1
    w00 = (1 - a_p) * (1 - a_t); w01 = (1 - a_p) * a_t; w10 = a_p * (1 - a_t); w11 = a_p * a_t

    def bil(s):
        return s[c00] * w00 + s[c01] * w01 + s[c10] * w10 + s[c11] * w11
    eaa = bil(SE[i])
    loss = 20 * np.log10(r / rs) + eaa * (r - rs) * 1e-3
    for k in S:
        lev[k][i] = bil(S[k][i]) - loss
    rng[i] = r

# reception-time weight (sel.hpp reception_factor)
drdt = np.empty_like(rng)
drdt[1:-1] = (rng[2:] - rng[:-2]) / (2 * T)
drdt[0] = (rng[1] - rng[0]) / T
drdt[-1] = (rng[-1] - rng[-2]) / T
wt = np.maximum(1 + drdt / C_SOUND, 1e-3)


def sel(l):
    peak = l.max(0)
    above = l >= (peak - 10)[None, :]
    first = np.argmax(above, 0)
    last = nst - 1 - np.argmax(above[::-1], 0)
    ii = np.arange(nst)[:, None]
    m = (ii >= first[None, :]) & (ii <= last[None, :])
    return 10 * np.log10(T * np.sum(np.where(m, 10 ** (0.1 * l) * wt, 0.0), 0)), first, last


res = {}
for k in S:
    res[k], f0, l0 = sel(lev[k])
    res[k] = res[k].reshape(X.shape)
    res[k + '_window'] = (l0 - f0).reshape(X.shape)
dref = res['a'] - ref
print('validation a vs NICEOPS: max |diff| %.3e dB, mean %.3e' % (np.abs(dref).max(), dref.mean()))
# fraction of observers whose window changed
for k in ['b', 'c', 'c2', 'd']:
    dd = res[k] - res['a']
    print('%-3s - a SEL: mean %+.3f median %+.3f p5 %+.3f p95 %+.3f min %+.3f max %+.3f  |d|>0.5: %.1f%%  |d|>1: %.1f%%' % (
        k, dd.mean(), np.median(dd), np.percentile(dd, 5), np.percentile(dd, 95), dd.min(), dd.max(),
        100 * np.mean(np.abs(dd) > 0.5), 100 * np.mean(np.abs(dd) > 1)))
    # restricted to observers within 15 dB of the footprint max (the contoured region)
    m = res['a'] > res['a'].max() - 20
    print('      within 20 dB of peak (%d obs): mean %+.3f p95 %+.3f max %+.3f min %+.3f' % (
        m.sum(), dd[m].mean(), np.percentile(dd[m], 95), dd[m].max(), dd[m].min()))
# instantaneous-level differences at the peak-window samples
for k in ['b', 'c', 'c2', 'd']:
    dl = (lev[k] - lev['a'])
    print('%-3s - a instantaneous (all state-observer pairs): median %+.3f p95 %+.3f max %+.3f min %+.3f' % (
        k, np.median(dl), np.percentile(dl, 95), dl.max(), dl.min()))
np.savez_compressed(out, ref=ref, xs=xs, ys=ys, **res)
