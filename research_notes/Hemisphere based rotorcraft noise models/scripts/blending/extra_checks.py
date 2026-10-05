"""(1) Energy blend under two valid triangulations (Qhull Qt vs QJ); the dB blend
is invariant (prism structure), the energy blend need not be.
(2) The CT axis alone: gap between dB, energy and the exact 20 log10 CT law.
(3) Load factor distribution of the trajectory states."""
import numpy as np
from blend import Conditions
from direction_stats import load_states

for db, veh in [('S-76D_M3', 'S-76D'), ('AW139_M1', 'AW139'), ('B205_M1', 'B205')]:
    z = np.load(db + '.npz'); L = z['grid']
    C = Conditions(z['cond']); CJ = Conditions(z['cond'], qhull_options='QJ Qbb Qc Q12')
    q = load_states(veh)
    w, idx, _ = C.barycentric(q); wj, idxj, _ = CJ.barycentric(q)
    Lk = np.moveaxis(L, -1, 0).reshape(L.shape[2], -1)
    b = 10 * np.log10(np.einsum('nk,nkd->nd', w, 10 ** (0.1 * Lk[idx])))
    bj = 10 * np.log10(np.einsum('nk,nkd->nd', wj, 10 ** (0.1 * Lk[idxj])))
    same = np.mean(np.all(np.sort(idx, 1) == np.sort(idxj, 1), 1))
    d = np.abs(bj - b)
    print(db, 'states', len(q), 'same simplex vertices in both triangulations: %.3f' % same,
          ' energy blend Qt vs QJ: p95 %.4f max %.4f dB' % (np.percentile(d, 95), d.max()))
    # CT-only: level = g(CT) = 20 log10(max(CT, CTref)/CTref), weights as used
    ct = C.raw[idx, 2]
    g = np.where(ct > C.ct_ref, 20 * np.log10(np.maximum(ct, 1e-30) / C.ct_ref), 0.0)
    gdb = np.sum(w * g, 1)
    gen = 10 * np.log10(np.sum(w * 10 ** (0.1 * g), 1))
    ctq = np.clip(q[:, 2], C.lo[2], C.hi[2])
    gex = np.where(ctq > C.ct_ref, 20 * np.log10(ctq / C.ct_ref), 0.0)
    lf = q[:, 2] / C.ct_ref
    print('   LF (CT/CT_ref) percentiles 1/50/99: %.3f %.3f %.3f  max %.3f' % (*np.percentile(lf, [1, 50, 99]), lf.max()))
    print('   CT axis only: energy-dB max %.4f, dB-exact min %.4f, energy-exact max %.4f' % (
        (gen - gdb).max(), (gdb - gex).min(), (gen - gex).max()))

# analytic: two levels L1 < L2 = L1 + D blended at weight u on L2
print('\nTwo-sphere blend, worst-case energy - dB gap over the weight, by level difference D:')
u = np.linspace(0, 1, 100001)
for D in [1, 2, 3, 5, 10, 15, 20, 30]:
    gap = 10 * np.log10((1 - u) + u * 10 ** (D / 10)) - u * D
    print('  D=%4.0f dB  max gap %.3f dB at u=%.3f; gap at u=0.5 %.3f' % (D, gap.max(), u[gap.argmax()],
                                                                       gap[50000]))
# CT law interval LF 3 -> 5 (the widest stored CT step), midpoint
for lo, hi in [(1.0, 1.1), (1.2, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 5.0)]:
    uu = np.linspace(0, 1, 10001); lf = lo + uu * (hi - lo)
    ex = 20 * np.log10(lf)
    dbl = (1 - uu) * 20 * np.log10(lo) + uu * 20 * np.log10(hi)
    en = 10 * np.log10((1 - uu) * lo ** 2 + uu * hi ** 2)
    print('  LF %.1f-%.1f: dB blend - exact min %.3f, energy - exact max %.3f' % (lo, hi, (dbl - ex).min(), (en - ex).max()))
