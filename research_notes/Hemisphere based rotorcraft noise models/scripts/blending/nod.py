"""Read a NICE-OPS .nod sphere database into a (phi, theta, condition) grid,
mirroring Database::Database in src/sphere_database.cpp (read-only)."""
import numpy as np
import netCDF4

FLOOR = -1000.0  # stands in for no_energy_floor_db; any level <= FLOOR+1 is "empty"


def approx_equal(a, b, rel=1e-2, abs_=1e-5):
    # D approxEqual: relative against the SECOND argument
    return np.abs(a - b) <= np.maximum(rel * np.abs(b), abs_)


def repair_legacy(phis):
    n = len(phis)
    if np.any(phis < -90 - 1e-9) or np.any(phis > 90 + 1e-9):
        return phis, False
    if n < 2 or n % 2:
        return phis, False
    h = n // 2
    if np.any(np.abs(phis[h:] + phis[:h]) > 1e-9):
        return phis, False
    p = phis.copy()
    c = 180.0 + p[h:]
    c = np.fmod(c + 180.0, 360.0)
    c[c < 0] += 360.0
    p[h:] = c - 180.0
    return p, True


def load(path):
    ds = netCDF4.Dataset(path)
    version = int(ds.variables['database_version'][0]) if 'database_version' in ds.variables else 0
    flf = int(ds.variables['fixed_load_factor'][0]) if 'fixed_load_factor' in ds.variables else 0
    names = sorted(ds.groups.keys(), key=lambda s: int(s.replace('sphere', '')))
    mu, gam, ct = [], [], []
    grid = None
    repaired = 0
    for i, g in enumerate(names):
        G = ds.groups[g]
        mu.append(float(G['advance_ratio'][0]))
        gam.append(float(G['flight_path_angle'][0]))
        ct.append(float(G['thrust_coefficient'][0]))
        phi = np.asarray(G['phi'][:], float)
        th = np.asarray(G['theta'][:], float)
        if version < 1:
            phi, r = repair_legacy(phi)
            repaired += r
        dba = np.ma.filled(G['dBA'][:], np.nan).astype(float)
        eaa = np.ma.filled(G['EAA'][:], np.nan).astype(float)
        if i == 0:
            phis = np.unique(phi)
            thetas = np.unique(th)
            grid = np.full((len(phis), len(thetas), len(names)), np.nan)
            eaa_grid = np.full_like(grid, np.nan)
        for c in range(len(phi)):
            m = np.nonzero(approx_equal(phis, phi[c]))[0]
            n = np.nonzero(approx_equal(thetas, th[c]))[0]
            if len(m) == 0 or len(n) == 0:
                continue
            v = dba[c] if np.isfinite(dba[c]) else FLOOR
            for mm in m:
                grid[mm, n, i] = v  # later channel wins a shared cell
                eaa_grid[mm, n, i] = eaa[c] if np.isfinite(eaa[c]) else 0.0
    ds.close()
    cond = np.array([mu, gam, ct]).T
    return dict(phis=phis, thetas=thetas, grid=grid, eaa=eaa_grid, cond=cond, version=version,
                fixed_load_factor=flf, repaired=repaired, names=names)
