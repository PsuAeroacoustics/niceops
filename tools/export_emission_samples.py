#!/usr/bin/env python3
"""Export depropagated samples with their flight condition at emission time.

Uses the existing niceops-validation recipe and preserves its 2,000 ft cap.
Adds emission times in memory to panam's scattered output; external files are
read only. Writes state-labelled NPZ samples and a guarded upstream patch to
--out. The accompanying CSV quantifies condition drift hidden by run means.
These samples are training data, not a validated replacement database.
"""

import argparse
import csv
import difflib
import importlib
import inspect
import json
from pathlib import Path
import sys

import numpy as np


def add_emission_time(source):
    changes = {
        '    fheight_list = []\n': '    fheight_list = []\n    femission_time_list = []\n',
        '        frange_list.append(r_v)\n':
            '        frange_list.append(r_v)\n        femission_time_list.append(np.asarray(track_time)[tidx][valid])\n',
        "            'range': np.concatenate(frange_list),\n":
            "            'range': np.concatenate(frange_list),\n            'emission_time': np.concatenate(femission_time_list),\n",
    }
    for old, new in changes.items():
        if source.count(old) != 1:
            raise ValueError('panam changed: cannot safely add emission times at ' + repr(old))
        source = source.replace(old, new)
    return source


def label_samples(times, track, tip_speed_mps):
    if not np.isfinite(tip_speed_mps) or tip_speed_mps <= 0:
        raise ValueError('Rotor tip speed must be finite and positive')
    times = np.asarray(times)
    if not np.all(np.isfinite(times)) or np.any(times < track['time'][0]) or np.any(times > track['time'][-1]):
        raise ValueError('Sample emission times lie outside the tracked record')
    speed = np.interp(times, track['time'], track['ground_speed_knots'])
    angle = np.interp(times, track['time'], track['fpa_deg'])
    return dict(emission_time_s=times, ground_speed_knots=speed, flight_path_angle_deg=angle,
                advance_ratio=speed * (1852.0 / 3600.0) / tip_speed_mps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--runs', nargs='+', required=True)
    args = parser.parse_args()
    harness, out = args.harness.resolve(), args.out.resolve()
    if out.exists():
        raise SystemExit('Use a fresh --out directory')
    sys.path.insert(0, str(harness))
    config = importlib.import_module('config')
    fa = importlib.import_module('flight_acoustics')
    na = importlib.import_module('noise_abatement_2017')
    sys.path.insert(0, str(harness / 'spheres'))
    rebuild = importlib.import_module('rebuild')
    build = harness / 'spheres/out/min2_stratified_shepfloor5/Be407'
    recipe = json.loads((build / 'build.json').read_text())
    test = na.NoiseAbatementTest('B407')
    # Only runs that were already used as training data may enter this export.
    with (build / 'manifest.csv').open(newline='') as handle:
        training = {row['run']: row for row in csv.DictReader(handle) if not row['error']}
    if any(run not in training or training[run]['condition'].startswith('A') for run in args.runs):
        raise SystemExit('Runs must be non-approach training runs in the original manifest')
    out.mkdir(parents=True)
    original_source = Path(fa.__file__).read_text()
    patched_source = add_emission_time(original_source)
    (out / 'panam_emission_time.patch').write_text(''.join(difflib.unified_diff(
        original_source.splitlines(True), patched_source.splitlines(True),
        fromfile='a/flight_acoustics.py', tofile='b/flight_acoustics.py')))
    # Compile only this function; keep the existing module and dependency state.
    function = add_emission_time(inspect.getsource(fa.depropagate_hemisphere))
    exec(compile(function, fa.__file__, 'exec'), fa.__dict__)
    depropagate = fa.depropagate_hemisphere
    captured = []
    channel_ids = []
    load_channels = na.load_run_channels

    def capture_channels(*positional, **keywords):
        result = load_channels(*positional, **keywords)
        channel_ids[:] = [str(mic) for mic in result[3]]
        return result

    na.load_run_channels = capture_channels

    def capture(*positional, **keywords):
        keywords['return_scattered'] = True
        result = depropagate(*positional, **keywords)
        captured.append(result)
        return result

    fa.depropagate_hemisphere = capture
    settings = dict(recipe['settings'])
    for key in ('steady_only', 'prefetch'):
        settings.pop(key, None)
    import configparser
    vehicle = configparser.ConfigParser()
    vehicle.read(build / 'vehicle.cfg')
    tip_speed = vehicle.getfloat('Main Rotor', 'tip speed')
    records = []
    with (out / 'condition_drift.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['run', 'samples', 'mean_speed_knots', 'sample_speed_min', 'sample_speed_max',
                         'mean_fpa_deg', 'sample_fpa_min', 'sample_fpa_max', 'sphere_control_max_db'])
        for run in args.runs:
            captured.clear()
            channel_ids.clear()
            rays = rebuild.ray_models('B407', str(out), recipe['rays'])(run) if recipe['rays'] != 'straight' else None
            reference = Path(rebuild.SPHERE_DATABASE) / 'Be407' / f'Be407{run[3:]}.nc'
            sphere = out / f'Be407{run[3:]}.nc'
            record = na.build_sphere(test, run, str(sphere), ray_model=rays,
                                     reference_sphere=str(reference) if reference.exists() else None, **settings)
            if len(captured) != 1:
                raise ValueError('Expected one depropagated sample set per run')
            hemisphere = captured[0]
            sc = hemisphere['scattered']
            labels = label_samples(sc['emission_time'], na.load_track(test.track_path(run)), tip_speed)
            powers = sc['third_octave']['bands_db']
            if powers.shape[1] != labels['emission_time_s'].size:
                raise ValueError('Sample labels and band levels lost alignment')
            np.savez(out / f'{run}.npz', azimuth_deg=sc['azi_deg'], elevation_deg=sc['elv_deg'],
                     band_centers_hz=sc['third_octave']['band_centers_hz'], bands_db=powers,
                     mic=sc['mic'], microphone_id=np.asarray(channel_ids)[sc['mic']],
                     range_ft=sc['range'], source_height_ft=sc['source_height'],
                     run=run, speed_reference='ground', tip_speed_mps=tip_speed,
                     reference_radius_ft=hemisphere['metadata']['r_ref'],
                     resolution_deg=sc.get('resolution_deg', np.zeros(len(sc['mic']))), **labels)
            old = np.asarray(fa.load_nc_sphere(str(build / sphere.name))[0], dtype=float)
            new = np.asarray(fa.load_nc_sphere(str(sphere))[0], dtype=float)
            if old.shape != new.shape or not np.array_equal(np.isfinite(old), np.isfinite(new)):
                raise ValueError(f'{run}: sphere control layout or finite mask differs from the original recipe')
            mask = np.isfinite(old)
            difference = float(np.max(np.abs(old[mask] - new[mask])))
            if difference > 1e-6:
                raise ValueError(f'{run}: recipe control changed by {difference} dB')
            speed, angle = labels['ground_speed_knots'], labels['flight_path_angle_deg']
            writer.writerow([run, len(speed), record['speed_knots'], speed.min(), speed.max(),
                             record['flight_path_angle_deg'], angle.min(), angle.max(), difference])
            records.append(record)
            print(f'{run}: {len(speed)} samples labelled; control max difference {difference:g} dB', flush=True)
    (out / 'manifest.json').write_text(json.dumps(dict(
        status='experimental_training_samples', recipe=recipe, records=records,
        panam=config.describe(config.PANAM), niceops=config.describe(config.NICEOPS_SRC),
        sample_condition='ground speed, mu and flight-path angle at emission time',
        limitations=['Source spheres remain labelled by run means until these samples are regridded',
                     'No new source condition or held-out accuracy claim is made by this export']), indent=2) + '\n')


if __name__ == '__main__':
    main()
