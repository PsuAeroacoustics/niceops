import glob
import sys
import numpy as np

veh = sys.argv[1]
files = sorted(glob.glob('/Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend/fp/%s/fp_*.npz' % veh))
pool = {k: [] for k in ['b', 'c', 'c2', 'd']}
poolm = {k: [] for k in ['b', 'c', 'c2', 'd']}
per = {k: [] for k in ['b', 'c', 'c2', 'd']}
val = []
win = []
for f in files:
    z = np.load(f)
    a = z['a']
    val.append(np.abs(a - z['ref']).max())
    m = a > a.max() - 20
    win.append(np.mean(z['b_window'] != z['a_window']))
    for k in pool:
        d = z[k] - a
        pool[k].append(d.ravel()); poolm[k].append(d[m].ravel())
        per[k].append((d.mean(), d[m].mean(), np.abs(d[m]).max()))
print(veh, 'tracks', len(files), 'max validation |a - NICE-OPS| %.2e dB' % max(val))
print('fraction of observers whose 10 dB window changed under energy blend: median %.3f max %.3f' % (np.median(win), max(win)))
for k, name in [('b', 'energy - dB'), ('c', 'nearest 3D - dB'), ('c2', 'nearest FC + CT law - dB'), ('d', 'NORAH IDW energy - dB')]:
    for lab, P in [('all observers', pool[k]), ('within 20 dB of track peak', poolm[k])]:
        d = np.concatenate(P)
        print('%-26s %-28s n=%7d mean %+.3f med %+.3f p5 %+.3f p95 %+.3f min %+.3f max %+.3f |d|p95 %.3f |d|>1dB %.1f%%' % (
            name, lab, d.size, d.mean(), np.median(d), np.percentile(d, 5), np.percentile(d, 95), d.min(), d.max(),
            np.percentile(np.abs(d), 95), 100 * np.mean(np.abs(d) > 1)))
    p = np.array(per[k])
    print('   per-track mean (all obs): range %+.3f to %+.3f; per-track mean (<=20 dB): %+.3f to %+.3f; per-track max |d| (<=20 dB): median %.2f max %.2f' % (
        p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max(), np.median(p[:, 2]), p[:, 2].max()))
