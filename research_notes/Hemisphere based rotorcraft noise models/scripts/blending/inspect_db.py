import time
import numpy as np
import nod

DB = '/Users/evg5332/Code/D/NICEOPS/.claude/worktrees/hemisphere-review/databases/'
for f in ['S-76D_M3', 'B407_M2', 'AW139_M1', 'B205_M1']:
    t = time.time()
    d = nod.load(DB + f + '.nod')
    np.savez_compressed(f + '.npz', phis=d['phis'], thetas=d['thetas'], grid=d['grid'], eaa=d['eaa'], cond=d['cond'])
    g = d['grid']; c = d['cond']
    print(f, 'v', d['version'], 'flf', d['fixed_load_factor'], 'rep', d['repaired'], 'shape', g.shape,
          'nan', np.isnan(g).sum(), 'floor', (g <= -999).sum(), 'time %.1f' % (time.time() - t))
    print(' phis', d['phis'][[0, 1, -1]], len(d['phis']), ' thetas', d['thetas'][[0, 1, -1]], len(d['thetas']))
    for k, nm in enumerate(['mu', 'gamma', 'CT']):
        u = np.unique(np.round(c[:, k], 6))
        print('  ', nm, len(u), u[:25])
    print('  level range', np.nanmin(g), np.nanmax(g))
