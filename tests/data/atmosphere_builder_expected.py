"""Oracle for the build_atmosphere_is_harness_profile and build_atmosphere_units_agree CTests.

It runs the validation harness's own profile step -- profile_validation.Weather.profile in
~/Code/Python/niceops-validation, the code that wrote every 2017 run's --atmosphere -- on a
synthetic run laid out as the 2017 data are (seven ZephIR gates at 46-984 ft in knots, one of them
a fill value; the tethered balloon's string every 10 ft to 500 ft in deg F; a station humidity),
and writes:

  tests/data/atmosphere_builder/lidar_{ft,m}.csv    the gates, as niceops_build_atmosphere reads them
  tests/data/atmosphere_builder/string_{ft,m}.csv   the string, in deg C
  tests/data/atmosphere_builder_expected_{ft,m}.csv the harness's profile, 17 digits, SI values

The harness rounds the file it writes (%.3f); this writes the rows Weather.profile returns, unrounded,
with T in K and the wind in m/s.  Regenerate with the harness's interpreter (it needs numpy):

  ~/Code/Python/niceops-validation/.venv/bin/python tests/data/atmosphere_builder_expected.py \\
      ~/Code/Python/niceops-validation

never from niceops_build_atmosphere.  The tool's arguments that reproduce it are in CMakeLists.txt.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FT = 0.3048
KNOT = 1852.0 / 3600.0

harness = os.path.abspath(os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                                             else '~/Code/Python/niceops-validation'))
sys.path.insert(0, os.path.join(harness, 'propagation'))
import profile_validation as pv  # noqa: E402  (the harness's module)

# One run, as the 2017 files give it.  The 696 ft gate is a fill value.
GATES = {46: (8.4, 262.0), 125: (10.9, 268.0), 197: (12.1, 271.0), 348: (13.8, 276.0),
         495: (9999.0, 9999.0), 696: (15.2, 284.0), 984: (16.0, 290.0)}
STRING_F = {h: 61.3 - 0.012 * h + (1.4 if h < 60 else 0.0) * math.exp(-h / 25.0)
            for h in range(0, 501, 10)}
HUMIDITY = 62.5
assert sorted(GATES) == pv.LIDAR_FT

weather = pv.Weather.__new__(pv.Weather)
weather.lidar = {'run': {**{f'{h}_Vh': repr(s) for h, (s, d) in GATES.items()},
                         **{f'{h}_dir': repr(d) for h, (s, d) in GATES.items()}}}
weather.temperature = {'run': {str(h): repr(t) for h, t in STRING_F.items()}}
rows = weather.profile('run', HUMIDITY)


def write(path, header, lines):
    with open(path, 'w') as handle:
        handle.write(header + '\n')
        for line in lines:
            handle.write(','.join(repr(float(x)) for x in line) + '\n')


for unit, scale in (('ft', 1.0), ('m', FT)):
    write(os.path.join(HERE, 'atmosphere_builder', f'lidar_{unit}.csv'),
          f'z_{unit},wind_speed_kt,wind_direction',
          [(h * scale, s, d) for h, (s, d) in GATES.items()])
    write(os.path.join(HERE, 'atmosphere_builder', f'string_{unit}.csv'), f'z_{unit},T_C',
          [(h * scale, (t - 32.0) / 1.8) for h, t in STRING_F.items()])
    write(os.path.join(HERE, f'atmosphere_builder_expected_{unit}.csv'),
          f'z_{unit},T_K,RH,u_mps,v_mps',
          [(z * scale, t + 273.15, rh, u * KNOT, v * KNOT) for z, t, rh, u, v in rows])
