#!/usr/bin/env python3
"""Grid source samples under their emission-time flight conditions, in energy.

Writes experimental AAM spheres and coverage counts, not a production database.
Condition bins keep speed within --speed_step/2 and flight-path angle within
--angle_step/2 of their centers. Each sample enters exactly one condition.
Unsupported angular directions stay missing; no range-loss fit is introduced.
"""

import argparse
import csv
import importlib
import json
from pathlib import Path
import sys

import numpy as np


def condition_bins(speed, angle, speed_step, angle_step):
    if not (np.isfinite(speed_step) and speed_step > 0 and np.isfinite(angle_step) and angle_step > 0):
        raise ValueError('Condition steps must be finite and positive')
    if not np.all(np.isfinite(speed)) or not np.all(np.isfinite(angle)):
        raise ValueError('Non-finite flight conditions')
    return np.column_stack((np.floor(np.asarray(speed) / speed_step + 0.5),
                            np.floor(np.asarray(angle) / angle_step + 0.5))).astype(int)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--samples', type=Path, nargs='+', required=True, help='Directories of emission-labelled NPZ files')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--speed_step', type=float, default=5.0, help='knots (default 5)')
    parser.add_argument('--angle_step', type=float, default=2.0, help='degrees (default 2)')
    parser.add_argument('--min_samples', type=int, default=200)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit('Use a fresh --out directory')
    condition_bins([], [], args.speed_step, args.angle_step)
    if args.min_samples < 1:
        raise SystemExit('--min_samples must be positive')
    sys.path.insert(0, str(args.harness.resolve()))
    importlib.import_module('config')
    fa = importlib.import_module('flight_acoustics')
    parts, paths = [], []
    for directory in args.samples:
        for path in sorted(directory.glob('*.npz')):
            with np.load(path) as data:
                parts.append({key: data[key] for key in data.files})
                paths.append(path.resolve())
    if not parts:
        raise SystemExit('No sample files')
    for data in parts:
        if 'microphone_id' not in data:
            raise ValueError('Export samples again to retain physical microphone identifiers')
        if str(data['speed_reference']) != 'ground':
            raise ValueError('Only ground-referenced source samples are supported')
        if not np.array_equal(data['band_centers_hz'], parts[0]['band_centers_hz']):
            raise ValueError('Sample band grids differ')
        if float(data.get('reference_radius_ft', 100.0)) != float(parts[0].get('reference_radius_ft', 100.0)):
            raise ValueError('Sample reference radii differ')
    join = lambda key: np.concatenate([data[key] for data in parts])
    speed, angle = join('ground_speed_knots'), join('flight_path_angle_deg')
    elevation, azimuth = join('elevation_deg'), join('azimuth_deg')
    resolution = np.concatenate([data.get('resolution_deg', np.zeros(len(data['mic']))) for data in parts])
    # Count physical microphones, not repeated runs or channel positions.
    mic = join('microphone_id')
    levels = np.concatenate([data['bands_db'] for data in parts], axis=1)
    bands = parts[0]['band_centers_hz']
    bins = condition_bins(speed, angle, args.speed_step, args.angle_step)
    el_axis, az_axis = np.arange(0., 92., 2.), np.arange(0., 362., 2.)
    az_grid, el_grid = np.meshgrid(az_axis, el_axis)
    args.out.mkdir(parents=True)
    used = 0
    with (args.out / 'coverage.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['speed_center_knots', 'angle_center_deg', 'samples', 'status', 'missing_nodes',
                         'sphere_speed_knots', 'sphere_angle_deg'])
        for speed_bin, angle_bin in np.unique(bins, axis=0):
            keep = (bins[:, 0] == speed_bin) & (bins[:, 1] == angle_bin)
            center_speed, center_angle = speed_bin * args.speed_step, angle_bin * args.angle_step
            count = int(keep.sum())
            if count < args.min_samples or np.unique(mic[keep]).size < 3:
                writer.writerow([center_speed, center_angle, count, 'insufficient_samples', '', '', ''])
                continue
            weights, _, gap = fa.adaptive_idw_weights(
                el_grid.ravel(), az_grid.ravel(), elevation[keep], azimuth[keep], mic[keep], resolution[keep],
                aspect=5.0, shepard_floor=1.0, relax_min_mics=False)
            power = np.asarray(weights @ (10 ** (0.1 * levels[:, keep])).T)
            power[gap, :] = np.nan
            with np.errstate(divide='ignore', invalid='ignore'):
                grid_levels = 10 * np.log10(power.T.reshape((len(bands),) + el_grid.shape))
            grid_levels[:, :, -1] = grid_levels[:, :, 0]
            hemisphere = dict(azi_grid_deg=az_axis, elv_grid_deg=el_axis,
                              third_octave=dict(band_centers_hz=bands, bands_db=grid_levels))
            mean_speed, mean_angle = float(speed[keep].mean()), float(angle[keep].mean())
            name = f'state_speed{center_speed:g}_angle{center_angle:g}.nc'
            fa.write_aam_hemisphere_netcdf(str(args.out / name), hemisphere, mode='third_octave',
                                          radius_ft=float(parts[0].get('reference_radius_ft', 100.0)),
                                          speed_knots=mean_speed, flight_path_angle_deg=mean_angle,
                                          title='Experimental emission-condition source; coverage.csv describes gaps')
            writer.writerow([center_speed, center_angle, count, 'written', int(gap.sum()), mean_speed, mean_angle])
            used += count
            print(f'{name}: {count} samples, {int(gap.sum())} missing angular nodes', flush=True)
    (args.out / 'manifest.json').write_text(json.dumps(dict(
        status='experimental_source_spheres', sample_files=[str(path) for path in paths],
        total_samples=len(speed), gridded_samples=used, speed_step_knots=args.speed_step,
        angle_step_degrees=args.angle_step,
        limitations=['Condition bins narrow angular coverage; missing nodes are retained',
                     'Held-out flight validation required before database replacement']), indent=2) + '\n')


if __name__ == '__main__':
    main()
