#!/usr/bin/env python3
"""How much a street microphone's SEL depends on the facades' scattering coefficient.

Synthetic LoD1 districts -- 7 by 7 blocks 60 m square, streets 20 m wide, roofs at H/W 0.5, 1
and 2 or of mixed heights, aligned with the track or turned 45 degrees -- are placed under and
beside the demo approach (tracks/S-76D/197421AC.csv, the B407 database, frame bearing 210).
Five microphones stand in each: on the centrelines of a street along the grid's x and one along
its y, mid-block; 2 m from a facade; a flush plate 1 m from a facade; and at an intersection.
Every wall absorbs 0.05 and scatters s of what it reflects, s in {0, 0.1, 0.2, 0.3, 0.5}; the
floor is rigid.  Each table is built with niceops_urban_table --buildings, and the run takes
plates for the flush microphones (as the urban CTests do).

Writes every SEL to <out>/results.csv and prints the summary FIDELITY quotes.  Run from the
repository root:

    tools/urban_scattering_study.py --build build --out /tmp/scattering_study -j 8
"""

import argparse
import csv
import math
import statistics
import subprocess
from pathlib import Path

FT = 0.3048
BLOCK, STREET, BLOCKS = 60.0, 20.0, 7
PITCH = BLOCK + STREET
SCATTERING = [0.0, 0.1, 0.2, 0.3, 0.5]
DISTRICTS = {'H/W 0.5': 10.0, 'H/W 1': 20.0, 'H/W 2': 40.0, 'mixed': None}
TURNS = [0.0, 45.0]
# District centres, ft in the run's frame: under the approach and 1000 ft aside, where the
# aircraft is near 400 ft and near 150 ft.
PLACES = [(-3000.0, 0.0), (-3000.0, 1000.0), (-1000.0, 0.0), (-1000.0, 1000.0)]
# Microphones in the grid's own frame (m from the district centre) and height (m).
MICROPHONES = [
    ('street along x', (0.0, PITCH / 2), 1.5),
    ('street along y', (PITCH / 2, 0.0), 1.5),
    ('2 m from a wall', (0.0, BLOCK / 2 + 2.0), 1.5),
    ('plate at a wall', (BLOCK / 2 + 1.0, 0.0), 0.0),
    ('intersection', (PITCH / 2, PITCH / 2), 1.5),
]
RUN = ['-d', 'databases/B407_fixed_load_factor.nod', '-v', 'databases/B407.json', '--humidity', '20',
       '-t', 'tracks/S-76D/197421AC.csv', '--frame_bearing', '210', '-L', 'warn',
       '--flow_resistance', '225', '--ground_effect', '--plate_table', 'plates/gr1425_db225.nc']


def place(centre_ft, turn, gx, gy):
    """A grid point (m) as local x, y (ft)."""
    c, s = math.cos(math.radians(turn)), math.sin(math.radians(turn))
    return (centre_ft[0] + (gx * c - gy * s) / FT, centre_ft[1] + (gx * s + gy * c) / FT)


def write_district(directory, height, turn, centre):
    buildings = directory / 'buildings.csv'
    with buildings.open('w', newline='') as f:
        out = csv.writer(f)
        out.writerow(['building', 'x_ft', 'y_ft', 'top_ft', 'absorption'])
        half = BLOCKS // 2
        n = 0
        for i in range(-half, half + 1):
            for j in range(-half, half + 1):
                n += 1
                top = height if height is not None else 10.0 + 30.0 * ((7 * i + 13 * j) % 10) / 9.0
                for gx, gy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                    x, y = place(centre, turn, i * PITCH + gx * BLOCK / 2, j * PITCH + gy * BLOCK / 2)
                    out.writerow([n, f'{x:.6f}', f'{y:.6f}', f'{top / FT:.6f}', 0.05])
    observers = directory / 'observers.csv'
    with observers.open('w', newline='') as f:
        out = csv.writer(f)
        out.writerow(['id', 'x_ft', 'y_ft', 'z_ft', 'height_ft'])
        for k, (_, (gx, gy), h) in enumerate(MICROPHONES):
            x, y = place(centre, turn, gx, gy)
            out.writerow([k + 1, f'{x:.6f}', f'{y:.6f}', 0, f'{h / FT:.6f}'])
    return buildings, observers


def levels(path):
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r]
    return [float(r[2]) for r in rows[1:]]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--build', default='build', help='The build directory (default build).')
    parser.add_argument('--out', required=True, help='Where to write the districts, tables and results.')
    parser.add_argument('-j', '--jobs', type=int, default=8, help='Threads for each table (default 8).')
    args = parser.parse_args()
    build, out = Path(args.build), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for name, height in DISTRICTS.items():
        for turn in TURNS:
            for centre in PLACES:
                tag = f"{name.replace('/', '').replace(' ', '')}_{int(turn)}_{int(centre[0])}_{int(centre[1])}"
                directory = out / tag
                directory.mkdir(exist_ok=True)
                buildings, observers = write_district(directory, height, turn, centre)
                open_file = directory / 'open.csv'
                subprocess.run([str(build / 'niceops'), *RUN, '-b', str(observers), '-o', str(open_file)],
                               check=True, capture_output=True)
                sel = {'open': levels(open_file)}
                for s in SCATTERING:
                    table = directory / f'table_{s}.nc'
                    level_file = directory / f'sel_{s}.csv'
                    subprocess.run([str(build / 'niceops_urban_table'), '--buildings', str(buildings), '--receivers',
                                    str(observers), '--database', 'databases/B407_fixed_load_factor.nod',
                                    '--frame_bearing', '210', '--facade_scattering', str(s), '-j', str(args.jobs),
                                    '-o', str(table)], check=True, capture_output=True)
                    subprocess.run([str(build / 'niceops'), *RUN, '-b', str(observers), '--urban_table', str(table),
                                    '-o', str(level_file)], check=True, capture_output=True)
                    sel[s] = levels(level_file)
                for k, (microphone, _, _) in enumerate(MICROPHONES):
                    rows.append({'district': name, 'turn': turn, 'x_ft': centre[0], 'y_ft': centre[1],
                                 'microphone': microphone, 'open': sel['open'][k],
                                 **{f's={s}': sel[s][k] for s in SCATTERING}})
                print(f'{tag} done', flush=True)

    with (out / 'results.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    def summary(key, label):
        print(f'\n{label}: SEL change from specular walls (s = 0), median [min, max] over 8 placements, dB')
        print('| ' + key + ' | street effect (s = 0) | ' + ' | '.join(f's = {s}' for s in SCATTERING[1:]) + ' |')
        print('| --- |' + ' --- |' * len(SCATTERING))
        for value in dict.fromkeys(r[key] for r in rows):
            group = [r for r in rows if r[key] == value]
            cells = []
            effect = [r['s=0.0'] - r['open'] for r in group]
            cells.append(f'{statistics.median(effect):+.2f} [{min(effect):+.2f}, {max(effect):+.2f}]')
            for s in SCATTERING[1:]:
                d = [r[f's={s}'] - r['s=0.0'] for r in group]
                cells.append(f'{statistics.median(d):+.2f} [{min(d):+.2f}, {max(d):+.2f}]')
            print(f'| {value} | ' + ' | '.join(cells) + ' |')

    summary('microphone', 'By microphone')
    summary('district', 'By district')
    all_d = {s: [r[f's={s}'] - r['s=0.0'] for r in rows] for s in SCATTERING[1:]}
    print('\nEvery microphone: ' + '; '.join(
        f's = {s}: median {statistics.median(d):+.2f}, max {max(d):+.2f}, min {min(d):+.2f} dB' for s, d in all_d.items()))


if __name__ == '__main__':
    main()
