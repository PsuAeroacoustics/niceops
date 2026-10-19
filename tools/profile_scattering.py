#!/usr/bin/env python3
"""A facade relief's scattering coefficient per third-octave band, by boundary elements.

The relief's profile, extruded along its other dimension (vertical ribs, a sawtooth in plan), scatters in the plane
of the profile.  Its correlation scattering coefficient (Mommertz and Vorlaender 1995; ISO 17497-1's s by another
route) is

    s(psi) = 1 - |sum p1 p0*|^2 / (sum |p1|^2 sum |p0|^2),

p1 the far field the sample scatters toward each direction of the half-plane in front of it, p0 that of a flat
reference of the same size, for a plane wave from psi.  Both samples are rigid boxes of the same outline in free
field, the relief on the front of one; their edges scatter alike and leave s alone.  The fields come from
thick_barrier_bem.py's constant elements (a wavelength / --per_wavelength long, graded toward corners), CHIEF
points inside the box, and a line source 1e5 m away for the plane wave.  s is averaged over psi from -85 to 85
degrees, each weighted by cos(psi) (the 2-D diffuse field's share), and over frequencies in each band.

    tools/profile_scattering.py --profile sawtooth --length 3.72 --bands 1600,2000,2500 --out s.csv

Profiles (lengths in metres):
    sawtooth: EDGE's scatter_1, teeth 0.12 deep rising over 0.24 and falling over 0.12, a 0.12 half-tooth at the start
    ribs:     EDGE's scatter_2, 0.04 wide and 0.024 deep, 0.04 apart
    corrugated: a 20 ft container's side wall, trapezoids 0.036 deep on a 0.278 pitch (flats 0.072 out and 0.068
              in, slopes 0.069 each way), as one owner describes his (engineersedge.com post 8153; ISO 1496-1 does
              not fix it, and doors and ends differ)
"""

import argparse
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thick_barrier_bem import C, double_layer, incident, mesh  # noqa: E402


def profile(kind, length):
    """The relief as (along, out) vertices from (0, 0) to (length, 0), outward +z."""
    if kind == 'flat':
        return [(0.0, 0.0), (length, 0.0)]
    if kind == 'sawtooth':
        v = [(0.0, 0.0), (0.0, 0.12), (0.12, 0.0)]
        s = 0.12
        while s + 0.36 <= length + 1e-9:
            v += [(s + 0.24, 0.12), (s + 0.36, 0.0)]
            s += 0.36
        if s < length - 1e-9:
            v.append((length, 0.0))
        return v
    if kind == 'ribs':
        v = [(0.0, 0.0)]
        s = 0.0
        while s + 0.08 <= length + 1e-9:
            v += [(s, 0.024), (s + 0.04, 0.024), (s + 0.04, 0.0), (s + 0.08, 0.0)]
            s += 0.08
        if s < length - 1e-9:
            v.append((length, 0.0))
        return v
    if kind == 'corrugated':
        v = [(0.0, 0.0)]
        s = 0.0
        while s + 0.278 <= length + 1e-9:
            v += [(s + 0.034, 0.0), (s + 0.103, 0.036), (s + 0.175, 0.036), (s + 0.244, 0.0), (s + 0.278, 0.0)]
            s += 0.278
        if s < length - 1e-9:
            v.append((length, 0.0))
        return v
    raise ValueError(kind)


def body(kind, length, backing):
    """The sample as a closed polygon, counterclockwise in (x, z): the relief on top of a box `backing` deep."""
    top = profile(kind, length)
    # Counterclockwise: along the bottom, up the far end, back along the relief, down the near end.
    return [(0.0, -backing), (length, -backing)] + top[::-1]


def coefficient(args):
    kind, length, f, per_wavelength, psis, angles = args
    k = 2 * np.pi * f / C
    fields = []
    for shape in (kind, 'flat'):
        vertices = body(shape, length, 0.2)
        # Drop repeated vertices.
        clean = [vertices[0]]
        for v in vertices[1:]:
            if np.hypot(v[0] - clean[-1][0], v[1] - clean[-1][1]) > 1e-12:
                clean.append(v)
        if np.hypot(clean[0][0] - clean[-1][0], clean[0][1] - clean[-1][1]) < 1e-12:
            clean.pop()
        a, b, normal = mesh(clean, C / f / per_wavelength, closed=True, grade=2)
        middle = 0.5 * (a + b)
        rng = np.random.default_rng(2)
        chief = np.stack([rng.uniform(0.02 * length, 0.98 * length, 40), rng.uniform(-0.19, -0.01, 40)], axis=1)
        A = np.vstack([0.5 * np.eye(len(a)) - double_layer(middle, a, b, normal, k, False),
                       double_layer(chief, a, b, normal, k, False)])
        far = 1e5
        sources = [np.array([0.5 * length + far * np.sin(psi), far * np.cos(psi)]) for psi in psis]
        rhs = np.stack([np.concatenate([incident(middle, s, k, False), -incident(chief, s, k, False)]) for s in sources],
                       axis=1)
        p = np.linalg.lstsq(A, rhs, rcond=None)[0]
        points = np.stack([0.5 * length + far * np.sin(angles), far * np.cos(angles)], axis=1)
        fields.append(double_layer(points, a, b, normal, k, False) @ p)
    p1, p0 = fields
    s = 1.0 - np.abs(np.sum(p1 * np.conj(p0), axis=0)) ** 2 / (np.sum(np.abs(p1) ** 2, axis=0) * np.sum(np.abs(p0) ** 2, axis=0))
    return f, s


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--profile', choices=['sawtooth', 'ribs', 'corrugated'], required=True)
    parser.add_argument('--length', type=float, required=True, help='the sample, m')
    parser.add_argument('--bands', default='1600,2000,2500,3150,4000,5000,6300,8000')
    parser.add_argument('--per_band', type=int, default=5, help='frequencies averaged in each band')
    parser.add_argument('--per_wavelength', type=float, default=10.0)
    parser.add_argument('--out', default='scattering.csv')
    parser.add_argument('-j', type=int, default=max(1, os.cpu_count() // 2))
    args = parser.parse_args()
    psis = np.radians(np.arange(-85.0, 86.0, 5.0))
    angles = np.radians(np.arange(-89.5, 90.0, 0.5))
    weight = np.cos(psis) / np.sum(np.cos(psis))
    bands = [float(b) for b in args.bands.split(',')]
    half = 2 ** (1 / 6)
    jobs = []
    for band in bands:
        lo, hi = band / half, band * half
        for m in range(args.per_band):
            jobs.append((args.profile, args.length, lo * (hi / lo) ** ((m + 0.5) / args.per_band), args.per_wavelength,
                         psis, angles))
    with Pool(args.j) as pool:
        results = dict(pool.imap_unordered(coefficient, jobs))
    with open(args.out, 'w') as out:
        out.write('band_hz,s,s_normal\n')
        for band in bands:
            lo, hi = band / half, band * half
            in_band = [s for f, s in results.items() if lo <= f <= hi]
            s_random = float(np.mean([np.dot(weight, s) for s in in_band]))
            s_normal = float(np.mean([s[len(psis) // 2] for s in in_band]))
            out.write(f'{band:g},{s_random:.4f},{s_normal:.4f}\n')
            print(f'{band:6g} Hz: s {s_random:.3f} (normal incidence {s_normal:.3f})', flush=True)


if __name__ == '__main__':
    main()
