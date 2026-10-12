"""Writes fit_run1.csv ... fit_run6.csv for the fit_z0_recovers_known CTest: five runs whose gates
(46, 125, 197 and 348 ft) follow the neutral log law exactly with z0 = 0.1 ft, each with its own
friction velocity and direction and a neutral temperature column (theta constant), and a sixth in
light wind that the screen leaves out.  The rows at 0 and 500 ft carry only temperature; their wind
cells are the fill value 9999.  python3 tests/data/atmosphere_builder/fit_runs.py"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FT = 0.3048
Z0 = 0.1 * FT
for run, (friction, direction) in enumerate([(0.30, 240.0), (0.42, 255.0), (0.25, 262.0),
                                             (0.55, 271.0), (0.36, 290.0), (0.04, 250.0)], 1):
    with open(os.path.join(HERE, f'fit_run{run}.csv'), 'w') as handle:
        handle.write('z_ft,wind_speed_mps,wind_direction,T_K\n')
        for ft in (0, 46, 125, 197, 348, 500):
            z = ft * FT
            t = 290.0 - 0.0098 * z
            if ft in (0, 500):
                handle.write(f'{ft},9999,9999,{t!r}\n')
            else:
                speed = friction / 0.4 * math.log(z / Z0)
                handle.write(f'{ft},{speed!r},{direction!r},{t!r}\n')
