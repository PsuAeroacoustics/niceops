#!/usr/bin/env python3
"""The exact field of a line source behind a rigid block on rigid ground, by boundary elements.

The reference the point-source model's double diffraction is held to (tests/test_point_paths.cpp,
FIDELITY "Urban receivers"): a block `width` by `height` standing on the ground, two-dimensional,
a line source before it and receivers behind it.  The exterior Neumann problem is solved with
constant elements collocated at their midpoints, in the Kirchhoff-Helmholtz (double layer) form

    (1/2) p(x) - int p(y) dG(x, y)/dn_y ds = p_inc(x),   G = (i/4) H0(k |x - y|)

with the half-space Green's function (the source's image in the ground), so only the block's three
faces above the ground are meshed; the time convention is exp(-i w t), as NICEOPS's.  CHIEF points
inside the block (Schenck's method; a least-squares solve) remove the interior resonances, the
first of which a 2.44 x 5.18 m block has near 78 Hz.  Elements are a wavelength / `--per_wavelength`
long, at most 0.25 m, halved three times toward each corner.

Writes `w,f,sx,sz,x,z,re,im` per receiver: p, the total field, for the source (i/4) H0(k r).  The
level re free field is 20 lg |p / ((i/4) H0(k |x - s|))|.  Run twice, at 30 and 60 elements a
wavelength, and take a third of the difference as the error: constant elements converge as h^2.
The test's values are within 0.05 dB on that estimate.

    tools/thick_barrier_bem.py --check
    tools/thick_barrier_bem.py --width 2.44 --height 5.18 --per_wavelength 60 --out bem60.csv

With --ridge the block has a gable roof, its ridge over the middle (tests/test_pitched_roofs.cpp):

    tools/thick_barrier_bem.py --width 8 --height 6 --ridge 9 --sources='-10,2;-25,4.5' --behind 2,5,10 \
        --heights 1.2,4,7 --per_wavelength 60 --out gable60.csv

With --roof any roof: its corners between the eaves, from the near one to the far one (a mansard
there, rays diffracted at four edges and more):

    tools/thick_barrier_bem.py --width 12 --height 6 --roof '1.5,7.5;5.25,9;6.75,9;10.5,7.5' \
        --sources='-15,0.5;-30,2' --behind 2,6,15 --heights 0.5,1.5,4 --per_wavelength 60 --out mansard60.csv

With --blocks, more flat blocks beyond the first, and receivers --behind the last of them (a street, the
rays crossing it diffracted between two buildings' edges):

    tools/thick_barrier_bem.py --width 10 --height 10 --blocks '20,30,10' --sources='-100,12;-50,15' \
        --behind 2,5,10 --heights 1.2 --frequencies 125,250,500 --per_wavelength 30 --out street30.csv

A dip, one block whose middle stands lower than its ends (touching blocks 10/8/10 m), is a --roof;
--chief_top 0.7 keeps the CHIEF points below the dip's floor (tests/test_profile_pe.cpp):

    tools/thick_barrier_bem.py --width 30 --height 10 --roof '10,10;10,8;20,8;20,10' --chief_top 0.7 \
        --sources='-100,5;-100,12;-50,15' --behind 2,10 --heights 9.0,9.5,9.8,9.9,9.95,9.99,10.01,10.05,10.1,10.2,10.5,11.0 \
        --frequencies 250 --per_wavelength 45 --out dip45.csv

Needs numpy and scipy.  --check compares against the exact series for a rigid circular cylinder
(0.005 dB at 800 Hz with 800 elements).
"""

import argparse
import itertools
import os
from multiprocessing import Pool

import numpy as np
from scipy.linalg import lstsq as scipy_lstsq
from scipy.special import h1vp, hankel1, jv, jvp

C = 343.2
GAUSS = {n: np.polynomial.legendre.leggauss(n) for n in (3, 8, 32)}


def mesh(vertices, h, closed, grade):
    """Elements along a polyline traversed counterclockwise in (x, z), so that the right-hand
    normal of each points out of the body: (start, end, normal) arrays."""
    starts, ends = [], []
    count = len(vertices)
    for s in range(count if closed else count - 1):
        p, q = np.asarray(vertices[s], float), np.asarray(vertices[(s + 1) % count], float)
        t = np.linspace(0.0, 1.0, max(int(np.ceil(np.linalg.norm(q - p) / h)), 2) + 1)
        if grade:
            t = np.sort(np.concatenate([t, [t[1] / 2**g for g in range(1, grade + 1)],
                                        [1.0 - (1.0 - t[-2]) / 2**g for g in range(1, grade + 1)]]))
        points = p[None, :] + t[:, None] * (q - p)[None, :]
        starts.append(points[:-1])
        ends.append(points[1:])
    a, b = np.concatenate(starts), np.concatenate(ends)
    tangent = (b - a) / np.linalg.norm(b - a, axis=1)[:, None]
    return a, b, np.stack([tangent[:, 1], -tangent[:, 0]], axis=1)


def double_layer(x, a, b, normal, k, ground):
    """M[i, j] = integral over element j of dG(x_i, y)/dn_y, by Gauss-Legendre: 32 points within
    two element lengths, 8 within eight, 3 beyond."""
    x = np.atleast_2d(x)
    length = np.linalg.norm(b - a, axis=1)
    distance = np.linalg.norm(x[:, None, :] - (0.5 * (a + b))[None, :, :], axis=2)
    out = np.zeros((len(x), len(a)), complex)
    near = length[None, :]
    for order, mask in ((32, distance < 2 * near), (8, (distance >= 2 * near) & (distance < 8 * near)),
                        (3, distance >= 8 * near)):
        i, j = np.nonzero(mask)
        if len(i) == 0:
            continue
        total = np.zeros(len(i), complex)
        for u, weight in zip(*GAUSS[order]):
            y = a[j] + 0.5 * (u + 1.0) * (b[j] - a[j])
            for image in ([False, True] if ground else [False]):
                xi = x[i] * ([1.0, -1.0] if image else [1.0, 1.0])
                r = y - xi
                R = np.linalg.norm(r, axis=1)
                with np.errstate(invalid='ignore', divide='ignore'):
                    value = -0.25j * k * hankel1(1, k * R) * np.einsum('ij,ij->i', r, normal[j]) / R
                value[R < 1e-12] = 0.0
                total += 0.5 * weight * length[j] * value
        out[i, j] = total
    return out


def incident(x, source, k, ground):
    x = np.atleast_2d(x)
    p = 0.25j * hankel1(0, k * np.linalg.norm(x - source[None, :], axis=1))
    if ground:
        p += 0.25j * hankel1(0, k * np.linalg.norm(x - source * [1.0, -1.0], axis=1))
    return p


def solve(a, b, normal, k, source, ground, chief):
    """The boundary's pressure, and the field at given points from it."""
    middle = 0.5 * (a + b)
    A = 0.5 * np.eye(len(a)) - double_layer(middle, a, b, normal, k, ground)
    rhs = incident(middle, source, k, ground)
    A = np.vstack([A, double_layer(chief, a, b, normal, k, ground)])
    rhs = np.concatenate([rhs, -incident(chief, source, k, ground)])
    try:
        p = np.linalg.lstsq(A, rhs, rcond=None)[0]
    except np.linalg.LinAlgError:
        # LAPACK's SVD driver gives up on some large fine meshes; a
        # complete orthogonal factorization does not.
        p = scipy_lstsq(A, rhs, lapack_driver='gelsy')[0]
    return lambda x: incident(x, source, k, ground) + double_layer(x, a, b, normal, k, ground) @ p


def block_case(args):
    width, height, f, source, receivers, per_wavelength, roof, blocks, chief_top = args
    k = 2 * np.pi * f / C
    h = min(C / f / per_wavelength, 0.25)
    rng = np.random.default_rng(1)
    meshes, chief = [], []
    # The roof's corners between the eaves (none for a flat roof), far to near; then any more flat
    # blocks (x0, x1, height), each its own profile and CHIEF points.
    for x0, x1, top, corners in [(0.0, width, height, roof), *[(x0, x1, top, []) for x0, x1, top in blocks]]:
        profile = [(x1, 0.0), (x1, top), *corners[::-1], (x0, top), (x0, 0.0)]
        meshes.append(mesh(profile, h, closed=False, grade=3))
        chief.append(np.stack([rng.uniform(x0 + 0.1 * (x1 - x0), x1 - 0.1 * (x1 - x0), 12),
                               rng.uniform(0.1 * top, chief_top * top, 12)], axis=1))
    a, b, normal = (np.concatenate(parts) for parts in zip(*meshes))
    field = solve(a, b, normal, k, np.asarray(source, float), True, np.concatenate(chief))
    return [(width, f, source[0], source[1], x, z, p.real, p.imag) for (x, z), p in zip(receivers, field(receivers))]


def check():
    """Against the series for a rigid cylinder of radius 1 m, line source 5 m from its axis."""
    for f in (200.0, 800.0):
        k = 2 * np.pi * f / C
        theta = np.linspace(0, 2 * np.pi, 400, endpoint=False)
        a, b, normal = mesh(np.stack([np.cos(theta), np.sin(theta)], axis=1), 10.0, closed=True, grade=0)
        field = solve(a, b, normal, k, np.array([5.0, 0.0]), False,
                      np.array([[0.1, 0.2], [-0.3, 0.15], [0.05, -0.4], [0.33, -0.1]]))
        for phi in (0.0, 1.0, 2.5, np.pi):
            exact = sum(0.25j * (1 if m == 0 else 2) * hankel1(m, 5 * k) * np.cos(m * phi)
                        * (jv(m, 3 * k) - jvp(m, k) / h1vp(m, k) * hankel1(m, 3 * k)) for m in range(int(5 * k) + 40))
            p = field(np.array([[3 * np.cos(phi), 3 * np.sin(phi)]]))[0]
            print(f'{f:.0f} Hz, {np.degrees(phi):5.1f} deg: {20 * np.log10(abs(p / exact)):+.4f} dB, '
                  f'{np.angle(p / exact):+.4f} rad')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--check', action='store_true', help='check against the rigid cylinder and stop')
    parser.add_argument('--width', type=float, default=2.44)
    parser.add_argument('--height', type=float, default=5.18)
    parser.add_argument('--sources', default='-10,2;-25,4.5;-15,6.5', help='x,z;x,z;... (m, the block at 0 < x < width)')
    parser.add_argument('--ridge', type=float, help='a gable roof, its ridge this high over the middle (m)')
    parser.add_argument('--roof', help="the roof's corners between the eaves, near to far: x,z;x,z;... (m)")
    parser.add_argument('--blocks', help='more flat blocks beyond the first: x0,x1,height;... (m)')
    parser.add_argument('--behind', default='2,5,10,25', help='receivers this far past the last block (m)')
    parser.add_argument('--heights', default='1.2,4.8,5.1,5.26,5.6,6.5,8', help='receiver heights (m)')
    parser.add_argument('--frequencies', default='125,250,500,1000')
    parser.add_argument('--per_wavelength', type=float, default=60.0)
    parser.add_argument('--chief_top', type=float, default=0.9,
                        help='CHIEF points no higher than this fraction of each profile\'s height (0.7 keeps them '
                             'out of a roof\'s dip, which is air)')
    parser.add_argument('--out', default='bem.csv')
    parser.add_argument('-j', type=int, default=os.cpu_count() // 2)
    args = parser.parse_args()
    if args.check:
        check()
        return
    sources = [tuple(float(v) for v in s.split(',')) for s in args.sources.split(';')]
    blocks = [tuple(float(v) for v in b.split(',')) for b in args.blocks.split(';')] if args.blocks else []
    far = max([args.width] + [x1 for _, x1, _ in blocks])
    receivers = np.array([(far + float(d), float(z)) for d in args.behind.split(',') for z in args.heights.split(',')])
    roof = [] if args.ridge is None else [(0.5 * args.width, args.ridge)]
    if args.roof:
        roof = [tuple(float(v) for v in p.split(',')) for p in args.roof.split(';')]
    jobs = [(args.width, args.height, float(f), s, receivers, args.per_wavelength, roof, blocks, args.chief_top)
            for f, s in itertools.product(args.frequencies.split(','), sources)]
    with Pool(args.j) as pool:
        rows = sorted(row for rows in pool.imap_unordered(block_case, jobs) for row in rows)
    with open(args.out, 'w') as out:
        out.write('w,f,sx,sz,x,z,re,im\n')
        for row in rows:
            out.write(','.join(repr(float(v)) for v in row) + '\n')


if __name__ == '__main__':
    main()
