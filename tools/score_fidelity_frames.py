#!/usr/bin/env python3
"""Compare two matching measurement-frame sets without selecting on model errors.

Outputs per-run energy-mean error by band group, elevation and range, plus pooled
scores. Common measurement-side SNR gates are applied once. Missing predictions
are counted and excluded from both models' score; new candidate gaps exit 2.
"""

import argparse
import csv
from pathlib import Path

import numpy as np


BANDS = [('50-160', 50, 160), ('200-630', 200, 630), ('315-1250', 315, 1250),
         ('800-2000', 800, 2000), ('2500-5000', 2500, 5000)]


def error(predicted, measured):
    # A shared scale prevents overflow without changing the energy ratio.
    shift = max(float(np.max(predicted)), float(np.max(measured)))
    return float(10 * np.log10(np.sum(10 ** ((predicted - shift) / 10))
                               / np.sum(10 ** ((measured - shift) / 10))))


def compare(baseline, candidate, hover=False):
    if baseline.shape != candidate.shape:
        raise ValueError('Frame counts differ')
    for name in ('mic', 'band', 'time', 'measured', 'snr_db', 'elevation', 'range'):
        if not np.array_equal(baseline[name], candidate[name], equal_nan=True):
            raise ValueError(f'Frame alignment or measurement differs: {name}')
    gate = (baseline['snr_db'] >= 10) & np.isfinite(baseline['measured'])
    if not hover:
        gate &= baseline['above_floor_db'] >= 3
    valid = (np.isfinite(baseline['predicted']) & np.isfinite(candidate['predicted'])
             & (baseline['predicted'] > -900) & (candidate['predicted'] > -900))
    return gate & valid, int(np.sum(gate & ~valid))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--hover', action='store_true', help='Do not apply a flyover floor gate to hovers')
    args = parser.parse_args()
    candidate_files = sorted(args.candidate.glob('*_frames.csv'))
    if not candidate_files:
        raise SystemExit('No candidate frames')
    rows = []
    gaps = 0
    new_gaps = 0
    coverage = []
    pooled = {}
    for path in candidate_files:
        a = np.atleast_1d(np.genfromtxt(args.baseline / path.name, delimiter=',', names=True))
        b = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        gate, missing = compare(a, b, args.hover)
        gaps += missing
        run = path.name[:6]
        measured = (a['snr_db'] >= 10) & np.isfinite(a['measured'])
        if not args.hover:
            measured &= a['above_floor_db'] >= 3
        baseline_valid = np.isfinite(a['predicted']) & (a['predicted'] > -900)
        candidate_valid = np.isfinite(b['predicted']) & (b['predicted'] > -900)
        introduced = int(np.sum(measured & baseline_valid & ~candidate_valid))
        new_gaps += introduced
        coverage.append([run, int(measured.sum()), int(np.sum(measured & ~baseline_valid)),
                         int(np.sum(measured & ~candidate_valid)), introduced])
        # Geometry bins are independent of the model's residual or shadow decision.
        selections = [('all', np.ones(a.size, dtype=bool))]
        for lo, hi in [(0, 2), (2, 5), (5, 10), (10, 17.5), (17.5, 25), (25, 90.1)]:
            selections.append((f'elevation_{lo:g}_{hi:g}', (a['elevation'] >= lo) & (a['elevation'] < hi)))
        for lo, hi in [(0, 2000), (2000, 4000), (4000, 8000), (8000, 20000)]:
            selections.append((f'range_{lo}_{hi}', (a['range'] >= lo) & (a['range'] < hi)))
        for name, geometry in selections:
            for label, lo, hi in BANDS:
                keep = gate & geometry & (a['band'] >= lo) & (a['band'] <= hi)
                if not np.any(keep):
                    continue
                baseline_error = error(a['predicted'][keep], a['measured'][keep])
                candidate_error = error(b['predicted'][keep], b['measured'][keep])
                rows.append([run, name, label, int(keep.sum()), baseline_error, candidate_error])
                pooled.setdefault((name, label), []).append((a['predicted'][keep], b['predicted'][keep], a['measured'][keep]))
    for (name, label), pieces in pooled.items():
        a, b, measured = [np.concatenate([piece[k] for piece in pieces]) for k in range(3)]
        rows.append(['pooled', name, label, measured.size, error(a, measured), error(b, measured)])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['run', 'geometry', 'band_group_hz', 'frames', 'baseline_error_db', 'candidate_error_db'])
        writer.writerows(rows)
    with args.output.with_name(args.output.stem + '_coverage.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['run', 'measurement_frames', 'baseline_missing', 'candidate_missing', 'new_candidate_gaps'])
        writer.writerows(coverage)
    print(f'{len(candidate_files)} runs, {gaps} missing predictions ({new_gaps} new candidate gaps), {len(rows)} scores: {args.output}')
    for row in rows:
        if row[0] == 'pooled' and row[1] in ('all', 'elevation_17.5_25'):
            print(row)
    return 2 if new_gaps else 0


if __name__ == '__main__':
    raise SystemExit(main())
