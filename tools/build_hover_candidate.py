#!/usr/bin/env python3
"""Build an experimental B407 hover and validate absolute levels on held-out cards.

Requires the local niceops-validation harness and its Python environment. Reads its
existing de-Dopplerized source and frame measurements; writes only to --out.
The low ring is HIGE and is not evidence for the OGE horizon. No uncertainty is
invented for the unmeasured 3--17.5 degree gap or elevations above 25 degrees.
"""

import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


HIGH = {f'H{k}': [283420 + k, 283430 + k] for k in range(7, 11)}
LOW = {f'H{k}': [283420 + k, 283430 + k] for k in range(3, 7)}


def design(azimuth, order):
    angle = np.radians(azimuth)
    columns = [np.ones_like(angle)]
    for m in range(1, order + 1):
        columns.extend((np.cos(m * angle), np.sin(m * angle)))
    return np.column_stack(columns)


def fit(rows, order):
    shape = design(rows['azimuth'], order)
    cards = np.unique(rows['card'])
    # Sum-to-zero nuisance offsets isolate the shape while retaining its mean.
    offsets = np.column_stack([(rows['card'] == card).astype(float) for card in cards])
    matrix = np.column_stack((shape, offsets[:, :-1] - offsets[:, -1:]))
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        raise ValueError('Insufficient angular coverage for this Fourier order')
    weight = np.sqrt(rows['count'])
    return np.linalg.lstsq(matrix * weight[:, None], rows['correction'] * weight, rcond=None)[0][:shape.shape[1]]


def energy_error(predicted, measured):
    return float(10 * np.log10(np.sum(10 ** (0.1 * predicted)) / np.sum(10 ** (0.1 * measured))))


def load_rows(directory, cards):
    parts = []
    for card, runs in cards.items():
        for run in runs:
            d = np.atleast_1d(np.genfromtxt(directory / f'{run}_frames.csv', delimiter=',', names=True))
            keep = (d['snr_db'] >= 10) & np.isfinite(d['predicted']) & np.isfinite(d['measured'])
            d = d[keep]
            p, t = np.radians(d['phi']), np.radians(d['theta'])
            azimuth = np.degrees(np.arctan2(np.sin(t) * np.sin(p), np.cos(t)))
            bins = np.minimum(((azimuth + 180) // 15).astype(int), 23)
            for band in np.unique(d['band']):
                if not 50 <= band <= 5000:
                    continue
                for cell in np.unique(bins):
                    mask = (d['band'] == band) & (bins == cell)
                    if mask.sum() >= 20:
                        parts.append((card, run, band, -172.5 + 15 * cell,
                                      -energy_error(d['predicted'][mask], d['measured'][mask]), int(mask.sum())))
    return np.array(parts, dtype=[('card', 'U3'), ('run', int), ('band', float),
                                 ('azimuth', float), ('correction', float), ('count', int)])


def validate(rows, order, ring, writer):
    # Order is specified before validation, never selected on these test folds.
    for split in ('heading', 'repeat'):
        held_out = np.unique(rows['card']) if split == 'heading' else [0, 1]
        for held in held_out:
            test = rows['card'] == held if split == 'heading' else (rows['run'] // 10) % 2 == held
            for band in np.unique(rows['band']):
                training = rows[(~test) & (rows['band'] == band)]
                testing = rows[test & (rows['band'] == band)]
                if not len(training) or not len(testing):
                    raise ValueError(f'Missing train/test data for {ring}, {split}, {held}, {band}')
                predicted = design(testing['azimuth'], order) @ fit(training, order)
                # Absolute error retains the held-out card's unknown wind offset.
                residual = predicted - testing['correction']
                weight = testing['count']
                bias = np.average(residual, weights=weight)
                rms = np.sqrt(np.average(residual ** 2, weights=weight))
                shape_rms = np.sqrt(np.average((residual - bias) ** 2, weights=weight))
                writer.writerow([ring, split, held, band, int(weight.sum()), bias, rms, shape_rms,
                                 np.sqrt(np.average(testing['correction'] ** 2, weights=weight))])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def high_only_correction(high_correction):
    def correct(phi, theta, frequency):
        p, t = np.radians(phi), np.radians(theta)
        elevation = np.degrees(np.arcsin(np.clip(np.sin(t) * np.cos(p), -1.0, 1.0)))
        outside = np.maximum(17.5 - elevation, elevation - 25.0)
        weight = 0.5 * (1.0 + np.cos(np.pi * np.clip(outside / 5.0, 0.0, 1.0)))
        return high_correction(phi, theta, frequency) * weight[:, :, None]
    correct.description = 'measured B407 high hovers only: full at 17.5--25 degrees, 5 degree cosine taper'
    return correct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--order', type=int, default=3, choices=range(0, 7))
    parser.add_argument('--score-only', action='store_true')
    parser.add_argument('--extent', choices=['high_only', 'two_ring'], default='high_only',
                        help='high_only preserves the low-ring synthesis; two_ring also applies the experimental HIGE proxy')
    args = parser.parse_args()
    harness, out = args.harness.resolve(), args.out.resolve()
    if out.exists():
        raise SystemExit('Use a fresh --out directory to preserve prior candidates')
    frames = harness / 'propagation/band_frames_hover_station_plate_dd_base'
    source = harness / 'spheres/out/hover_source/dedoppler/Be407/Be407103.nc'
    spheres = harness / 'spheres/out/min2_stratified_shepfloor5/Be407'
    recipe = json.loads((spheres / 'build.json').read_text())
    rows = {'high': load_rows(frames, HIGH), 'low': load_rows(frames, LOW)}
    out.mkdir(parents=True)
    with (out / 'held_out.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['ring', 'split', 'held_out', 'band_hz', 'frames', 'bias_db',
                         'absolute_rms_db', 'offset_removed_rms_db', 'baseline_rms_db'])
        for ring, data in rows.items():
            validate(data, args.order, ring, writer)
            bands = np.unique(data['band'])
            coefficients = np.array([fit(data[data['band'] == band], args.order) for band in bands])
            np.savez(out / f'{ring}.npz', bands=bands, coefficients=coefficients, order=args.order,
                     ring_elevation=21.5 if ring == 'high' else 1.8)
    inputs = [source, spheres / 'build.json', spheres / 'vehicle.cfg']
    inputs += sorted(spheres.glob('*.nc'))
    inputs += [frames / f'{run}_frames.csv' for cards in (HIGH, LOW) for runs in cards.values() for run in runs]
    provenance = dict(status='experimental', order=args.order, extent=args.extent, harness=str(harness),
                      recipe=recipe, inputs={str(path): sha256(path) for path in inputs},
                      limitations=['HIGE low ring is used as a proxy for the OGE horizon',
                                   '3--17.5 degrees unmeasured', 'above 25 degrees extrapolated',
                                   '63 Hz processing disagreement unresolved',
                                   'held-out scores use cached forward predictions and binned corrections'])
    for name, repo in [('harness', harness), ('panam', harness / 'external/panam')]:
        provenance[name + '_revision'] = subprocess.run(
            ['git', '-C', str(repo), 'describe', '--always', '--dirty'],
            check=True, capture_output=True, text=True).stdout.strip()
    if not args.score_only:
        sys.path.insert(0, str(harness))
        importlib.import_module('config')
        fa = importlib.import_module('flight_acoustics')
        sys.path.insert(0, str(harness / 'spheres'))
        hover = importlib.import_module('hover_sphere')
        from netCDF4 import Dataset
        # De-Dopplerization can empty a band at the recorded spectrum's edge.
        # Keep the previous synthesis at these cells rather than losing coverage.
        hover_source = out / 'hover_source.nc'
        shutil.copy2(source, hover_source)
        with Dataset(hover_source, 'r+') as target, Dataset(spheres / 'Be407103.nc') as control:
            emitted = np.asarray(target['AMPLITUDE'][:], dtype=float)
            original = np.asarray(control['AMPLITUDE'][:], dtype=float)
            for axis in ('PHI', 'THETA', 'FREQUENCY'):
                if not np.array_equal(target[axis][:], control[axis][:]):
                    raise ValueError('Hover source and control grids differ: ' + axis)
            missing = ~np.isfinite(fa.mask_missing_levels(emitted))
            supported = np.isfinite(fa.mask_missing_levels(original))
            fallback = missing & supported
            emitted[fallback] = original[fallback]
            target['AMPLITUDE'][:] = emitted
            target.coverage_fallback = 'original flight synthesis where emitted-frequency source has no data'
            provenance['source_coverage_fallback_cells'] = int(fallback.sum())
        provenance['hover_source_sha256'] = sha256(hover_source)
        directory = out / 'Be407'
        directory.mkdir()
        for path in sorted(spheres.glob('*.nc')):
            (directory / path.name).symlink_to(path)
        shutil.copy2(spheres / 'vehicle.cfg', directory / 'vehicle.cfg')
        database = directory / 'Be407_fixed_load_factor.nod'
        correction = (hover.two_ring_correction_function(str(out / 'high.npz'), str(out / 'low.npz'))
                      if args.extent == 'two_ring' else
                      high_only_correction(hover.fourier_correction_function(str(out / 'high.npz'))))
        fa.build_empirical_database(str(directory), str(database), load_factors=None, store_spectrum=True,
                                    extended_flight_path_angles=recipe['extended_flight_path_angles'],
                                    hover_source=str(hover_source),
                                    hover_correction=correction)
        provenance['database_sha256'] = sha256(database)
    (out / 'manifest.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'Candidate and held-out scores: {out}')


if __name__ == '__main__':
    main()
