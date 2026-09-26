# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **OP**eration Planning **S**ystem — a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

NICE-OPS was developed by Rob Rau and Eric Greenwood. It began life in D; this is the C++17
implementation, which reproduces the D original's results on the demo case to within a few
units in the last place.

## Building

Requirements: a C++17 compiler and netCDF. Plotting additionally needs Python 3 with
development headers; matplotlib itself is installed by the build, not by you. On Ubuntu:

```sh
sudo apt-get install libnetcdf-dev python3-dev
```

Plotting can be dropped entirely, which removes the Python requirement with it:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_PLOTTING=OFF
```

CLI11, nlohmann_json, Catch2 and pybind11 are fetched at configure time, so the first configure
needs network access.

Build out of source, always with explicit `-S`/`-B`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

On macOS with Command Line Tools but no full Xcode, add
`-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`.

Link-time optimization is on wherever the toolchain supports it, since the acoustic loop and
the database it calls live in different translation units; `-DNICEOPS_LTO=OFF` turns it off if
a build needs to be faster to produce than to run.

Run the tests with `ctest --test-dir build -C Release`.

## Running

With no arguments, NICE-OPS runs a demo case: the trajectory `tracks/S-76D/197421AC.csv`
against the database `databases/S-76D_M3.nod` and vehicle `databases/S-76D.json`, writing
`integrated/197421_LAE.csv`.

```sh
./build/niceops                                     # the demo case
./build/niceops -t my_track.csv -d my_db.nod -v my_vehicle.json -o out.csv
./build/niceops --help                              # the full option list
```

Paths are relative to the working directory, so run from the repository root when using the
shipped data.

### Inputs

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV. Column 0 is UTC seconds, column 1 elapsed time; positions and velocities are in **feet**. |
| Trajectory (`-t`), WGS84 | A CSV whose header names `time`, `lat`, `lon` and `alt` columns — see [Geodetic coordinates](#geodetic-coordinates). Any column order, and optional `vE`/`vN`/`vUp`. |
| Database (`-d`) | netCDF `.nod`: one group per operating condition (advance ratio μ, flight path angle γ, thrust coefficient C<sub>T</sub>), each holding scattered (φ, θ) channels of dBA level and EAA excess-attenuation factor. |
| Vehicle (`-v`) | JSON with `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg) and `Atmosphere.density`. |
| Observers (`-b`) | CSV of scattered observer positions, in feet. Without it, a 120×120 grid spanning −1500…1475 m is used. |

Input and output files use **feet**; the computation is metric internally.

### Outputs

- `-o` — sound exposure level (SEL) per observer. Gridded output writes the x and y axes on the
  first two lines and then one row of levels per y.
- `-D` — day-night level, summed over every trajectory in the run (see below).
- `-i` — per-observer level time histories.
- `-m` — appends average, maximum and above-threshold counts for the case.
- `-S` — the derived trajectory state (t, position, velocity, load factor, C<sub>T</sub>, μ, flight path angles).
- `--powerfile` — per-state shaft power, if the database carries any.
- `--latlon` — writes a `_latlon` companion beside every level file above, as latitude, longitude, level.
- `--export_sphere` — a whole interpolated noise sphere, written as a one-group `.nod`.

## Units

`--units` selects the length unit for the whole interface: `ft` (the default) or `m`. It covers

- every length on the command line — `--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
  `--receiver_height`, `--path_smoothing`;
- every file carrying a length — trajectory input, observer lists, level output, time histories;
- the axis labels on the figures.

The computation is metric throughout. The conversion happens once, where the command line is
parsed, and nothing downstream of that seam sees anything else.

**`--units` says what your files already are; it does not convert them.** The shipped tracks and
observer sets are in feet, so running them with `--units m` reads feet as metres and flies the
aircraft three times too far. Nothing can detect that, so the flag is a declaration you have to
get right. Where a file names its own unit in a column header — `alt_ft`, `alt_m` — the name is
checked against `--units` and a disagreement fails the file rather than being silently believed;
a bare `alt` claims nothing and is always accepted.

Two things stay put whatever you select:

- **Defaults are physical quantities, not numbers in the current unit.** The observer grid
  defaults to ±1500 m in steps of 25 m whichever unit is in force — it reads as ±4921.26 ft or
  ±1500 m, the same ground either way. If defaults scaled with the flag, `--units m` would
  quietly shrink the footprint to a third of its area.
- **The `-S` state file is metric.** It is a dump of the internal state, and
  `golden/state_197421.csv` — the D-derived oracle — encodes it that way.

Non-length quantities keep the units of the standards they come from, because converting them
would make them harder to check against the source: temperature in K, pressure and flow
resistance in kPa (ISO 9613-1 is written in kPa), frequency in Hz, angles in degrees, time in
seconds, levels in dB.

> Before this flag existed, `--grid`, `--gridstep`, `--receiver_height` and `--path_smoothing`
> were metres while `--trim` and `--ground_elevation` were feet and every file was feet. Scripts
> written against that need those four converting, or `--units m` adding.

## Command line arguments
## Command line arguments

### Inputs and outputs

| Option | Meaning |
| --- | --- |
| `-d`, `--database` | netCDF noise database. Once, or once per trajectory. |
| `-t`, `--trajectoryfile` | CSV trajectory. Repeat for several flight paths. |
| `-v`, `--vehicle` | JSON vehicle definition. Once, or once per trajectory. |
| `-o`, `--output` | SEL output file. Once per trajectory; with several tracks and no `-o`, each output is named after its track. |
| `-b`, `--observer` | Scattered observer file; switches off the default grid. |
| `-S`, `--state` | Save the derived trajectory state. |
| `-i`, `--history` | Prefix for per-observer time history files. |
| `-m`, `--metricfile` | Append overall noise metrics to this file. |
| `--origin` | WGS84 reference for the local frame: latitude longitude, degrees. Defaults to the first position of the first geodetic track. |
| `--frame_bearing` | Compass bearing of the local +x axis. Default: x east, y north. The shipped AC tracks use 210. |
| `--ground_elevation` | Elevation of the observer plane above mean sea level, in `--units`. Subtracted from a WGS84 track's altitude column. |
| `--path_smoothing` | How far, RMS per axis in `--units`, the path fitted through a WGS84 track may sit from the given positions. Default 10 ft. Zero interpolates exactly, which on measured positions makes the load factor track coordinate rounding. |
| `--start_utc` | Seconds since UTC midnight at *t* = 0 of a WGS84 track, for the day-night window. |
| `--latlon` | Also write each level file as latitude, longitude, level. |
| `-l`, `--logfile` | Log file (default `log/niceops.log`). |
| `-L`, `--loglevel` | `trace`, `info`, `warn`, `error` or `fatal` (default `info`). |

### Physics and metrics

| Option | Meaning |
| --- | --- |
| `-a`, `--ambient` | Ambient level in dB, mixed into every sample. Off unless given. |
| `-q`, `--qsam` / `--no-qsam` | Use the flight path angle relative to the tip-path-plane normal rather than the vertical. On by default. |
| `-s`, `--step` | Resampling step for trajectory adaptation. **Not** the SEL integration step, which comes from the trajectory's own sampling. |
| `--units` | Length unit for the command line and for every file carrying a length: `ft` (default) or `m`. See **Units**. |
| `--trim` | Crop each trajectory from this x, in `--units`, to its lowest/slowest point. Off by default; it changes which part of the flight is integrated. |
| `--grid`, `--gridstep` | Observer grid bounds and spacing, in `--units`. Default ±1500 m (±4921.26 ft), step 25 m (82.021 ft) — a fixed extent either way. |
| `-D`, `--dnlfile` | Write the day-night level, summed over every trajectory. |
| `--utc_offset` | Hours from UTC to local time, for the DNL night window. |
| `--powerfile` | Write per-state shaft power, when the database has it. |
| `--spectral` | Propagate the database's per-band spectra to the ground with ISO 9613-1 absorption, instead of one broadband level with the database's EAA factor. Needs a database carrying spectra. |
| `--band_limits` | Frequency range integrated into the metric, Hz: `low high`. Default `20 20000`. |
| `--humidity`, `--temperature`, `--pressure` | Atmosphere for `--spectral`, overriding the vehicle file. Percent, K, kPa. |
| `--ground_effect` | Add excess ground attenuation (Chien–Soroka over a Delany–Bazley impedance) to the spectral path. Needs `--spectral`. |
| `--flow_resistance` | Ground flow resistance, kPa·s/m². 30 snow, 225 grass (default), 650 soil, 3000 dirt, 10000 concrete, 50000 asphalt. |
| `--receiver_height`, `--turbulence` | Receiver height (in `--units`, default 1.2 m) and turbulence (rad·s/√m, 0–16e-4) for `--ground_effect`. |
| `--impedance_model` | Ground impedance: `delany_bazley` (default), `miki`, `variable_porosity` (with `--porosity_rate`, 1/m) or `delany_bazley_layer` (with `--layer_depth`, in `--units`). `--flow_resistance` is the (effective) flow resistivity. |
| `--turbulence_model`, `--gamma_t` | `chessell` (default, `--turbulence`) or `harmonoise`, whose coherence takes `--gamma_t` (~1e-6 moderate, ~1e-5 strong). |
| `--roughness` | RMS roughness height, in `--units`: the coherent reflection falls by exp(−2(kσ cos θ)²). Default 0. |
| `--plate_table` | Observers are flush microphones on a rigid ground plate (GRAS 67AX in a GR1425 plate). Takes a netCDF table from panam; needs `--ground_effect` over the ground the table was built for. See below. |
| `--exact_ground` | Evaluate the ground model (or the plate's reflection coefficients) at every observer and state instead of interpolating per-trajectory tables of it: 3-4x slower, for checking the tables. |
| `--plate_mic_azimuth` | Which way the microphone is offset on its plate: `outboard` (default, away from the track) or degrees in the local frame from +x toward +y. |

### Performance

| Option | Meaning |
| --- | --- |
| `-j`, `--threads` | Worker threads for the observer loop. 0 (default) uses one per hardware thread; 1 runs serially. |
| `--db_cutoff` | Skip evaluations that provably fall below this level. See the caveat below. |

### Spheres and diagnostics

| Option | Meaning |
| --- | --- |
| `--export_sphere` | Interpolate a whole noise sphere and write it as a one-group `.nod` database. |
| `--export_aam` | The same sphere in AAM/RNM format. Needs a database with spectra. |
| `--sphere_condition` | The operating condition for it: advance ratio, flight path angle, C<sub>T</sub>. |
| `--self_test` | Validate every sphere in the database against its own stored levels. |
| `-p`, `--plotting` | Plot the footprint with the flight track over it; `--cmin`/`--cmax`/`--cstep` set the contour levels. |
| `-f`, `--figfile`, `-F`, `--statefigfile` | Save the footprint and trajectory-state figures to files. |

## Notes on particular features

### Several trajectories, and the day-night level

`-t` may be repeated. `-d`, `-v`, `-o` and `-S` then take either one value, shared by every
trajectory, or one per trajectory.

```sh
./build/niceops -t morning.csv -t evening.csv -t night.csv -D integrated/dnl.csv
```

DNL is a property of a set of operations rather than of one flight: each trajectory's SEL-style
integral gains a 10 dB penalty on samples flown at night, and the results are summed as energy.
Night is 22:00–07:00 local; `--utc_offset` shifts the trajectory's UTC seconds into local time.
With one daytime operation the DNL equals the SEL exactly.

Plotting is skipped for multi-trajectory runs, since one figure file cannot hold several paths.

A trajectory that cannot be read is logged and skipped rather than ending the run, so one bad
file in a batch costs that file alone. The exit status is nonzero if anything was skipped.

### `--db_cutoff`

On a large observer grid most of the work goes into observers that are never loud enough to
matter. `--db_cutoff L` skips any evaluation that provably falls below `L`, bounding the level
from the loudest sphere the state blends and spherical spreading alone.

The threshold sits 10 dB below `L`, which makes the result **exact wherever it matters**: the SEL
integral covers the window from the peak down to peak − 10, so for any observer that reaches `L`,
no sample in its window is ever culled. On the demo case, `--db_cutoff 85` skips 96% of the
evaluations and every cell at or above 85 dB is unchanged to the last bit.

Observers that never reach the cutoff come out as `-inf` — they contributed no energy, and those
cells will not plot. That is the trade the flag asks for; it is off by default.

### Spectra

A database may carry per-band spectra alongside the overall A-weighted level: a
`frequency` dimension and an `amplitude(phi, theta, frequency)` array per sphere. When it
does, `--export_aam` can write an interpolated sphere back out in the AAM/RNM hemisphere
format the source spheres arrive in. Without spectra the flag refuses, since
`AMPLITUDE(PHI, THETA, FREQUENCY)` is the whole payload of that format.

None of the databases in `databases/` carries spectra. Build one with panam's
`build_empirical_database(..., include_spectra=True)`, which now carries the per-band
amplitudes through from the source spheres. Band-limit it — `infreqs=(10.0, 3500.0)` —
because the top few bands of a measured sphere are usually noise-floor artifacts.

The AAM/RNM source spheres live in the group's OneDrive under `Data/Sphere_Database`, one
directory of `.nc` files per aircraft (AS350, Be206L3, Be407, Be430, EC130, MD-902, Mi-8, R-44,
R-66, SH-60B), each with its own `vehicle.cfg` and `vehicle.json`. `Be407_spectral.nod` sits
beside them — 26 bands, 10–3150 Hz, built from `Be407/` as above. Databases are not tracked in
this repository; they are far too large and are derived products.

#### A note on levels changing

**Results from before September 2026 are not comparable to results after.** `eval` used to
multiply the database's `EAA` by a constant 4.0, under a `//TODO Add ground effects` inherited
from the D original. That is not a ground effect, and it overstated atmospheric absorption
fourfold — 168 dB at 6 km where the physical value is 14. panam *writes* that `EAA` as
`SPLA - SPLAa` over a 1000 m reference and consumes it as `EAA / distance * (slant_range -
radius)`; the tool now does the same.

Removing it raised the demo footprint by **9.93 dB mean, up to 22.7 dB**, most of it far from
the track where the old factor was removing the most. `golden/197421_LAE.csv` was regenerated;
`golden/state_197421.csv` is trajectory state and unaffected.

It was measured against `--ground_effect`, the model it was standing in for, before being
removed: soft ground costs a roughly range-independent 2–4 dB relative to a rigid surface,
where the factor removed 7 dB at 250 m rising to 38 dB at 2 km. A real ground effect is bounded
by the surface impedance; a multiplier on a distance-proportional absorption is not.

#### Propagating them to the ground

`--spectral` carries the spectrum itself to the observer instead of a single number. Each
band is spread by the same `20 log10(r/r_sphere)` and absorbed at its own ISO 9613-1
coefficient over `r - r_sphere`, and only then are the bands A-weighted and summed into a
level. The two paths answer different questions and will not agree:

| | default | `--spectral` |
| --- | --- | --- |
| Attenuation | database's `EAA`, one dB/km factor | ISO 9613-1 α(f), per band |
| Fitted for | the spectrum the source happened to have | the atmosphere you state |
| Atmosphere | none — `EAA` is baked in | temperature, pressure, humidity |
| Needs | any database | one carrying spectra |

The absorption model is ported from panam's `panam_acoustics/iso_9613_1_1993.py` and pinned
against it in `tests/test_atmosphere.cpp`; the A-weighting is panam's `dBAw`, pinned in
`tests/test_weighting.cpp`.

**Band-limit deliberately.** Summing whatever bands a database happens to carry quietly
makes the answer a property of that database's bandwidth, so two databases are not
comparable. `--band_limits low high` states the range; bands outside it contribute exactly
nothing, and the run logs how many of the database's bands survived. A limit that excludes
every band is an error rather than silence.

The atmosphere comes from the vehicle JSON's `Atmosphere` section:

```json
"Atmosphere": {
    "density": "1.225",
    "temperature": "300.71",
    "pressure": "105.74",        // optional; else derived as rho*R*T
    "relative humidity": "20"    // optional; else 20%, panam's default
}
```

Only `density` is required, as before. Pressure defaults to the one implied by that file's
own density and temperature rather than to a standard day — `S-76D.json` pairs 1.225 kg/m³
with 300.71 K, which is not an ISA point, so assuming 101.325 kPa would contradict its own
numbers. `--humidity`, `--temperature` and `--pressure` override all of it.

#### Ground effect

`--ground_effect` adds the interference between the direct ray and the one reflected off the
ground — Chien–Soroka over a Delany–Bazley impedance, with Chessell's boundary-loss correction
and band averaging, ported from panam's `ega()` and pinned against it in
`tests/test_ground_effect.cpp`.

**It is a gain, not a loss**, despite the name: over rigid ground with the receiver on the
surface it tends to +6 dB, the pressure doubling of a perfect reflector. It goes sharply
negative only in the interference dip, where the reflected path arrives out of phase. On the
Be407 demo case:

| | mean | max |
| --- | --- | --- |
| spectral, no ground | 65.10 | 96.38 |
| + concrete (10000) | 68.45 (**+3.35**) | 99.00 |
| + grass (225) | 64.73 (**−0.36**) | 97.89 |

Hard ground is 3.7 dB louder than soft, and grass roughly cancels its own reflection gain
against the dip. `--flow_resistance` spans five orders of magnitude and the answer is sensitive
to it, so name the surface rather than guessing.

`--turbulence` decorrelates the two paths and fills the nulls in; leaving it at zero gives
unrealistically deep interference dips.

The impedance, turbulence and roughness models follow Attenborough & Taherzadeh, *Prediction
of outdoor ground effect*, Applied Acoustics 242 (2026) 111067, as panam's `ground_plane`
module does. They are pinned against `ground_plane.pole_level` in
`tests/test_ground_effect.cpp`. Delany–Bazley is poor for grassland. On the 2017 Amedee pole
microphones, variable porosity with σe ≈ 200 fitted best, with HARMONOISE γT ≈ 3e-5 and
10–15 mm roughness (panam `docs/ground_plane_corrections.md`).

#### Ground-plane microphones

A flush ground-plane microphone sits on a plate: rigid close in, soft ground beyond. The point
receiver can model one surface or the other, not both. `--plate_table` models the plate
itself, using panam's axisymmetric boundary-element solution for a GRAS GR1425 plate (400 mm,
8 mm thick, lying on the ground) with the microphone at 3/4 R per SAE ARP 4055. The solution is
stored as two complex transfer functions of the direct and image waves:

    p / p_free = P_d + Q P_r

P_d and P_r are tabulated over frequency (five sub-frequencies per band), elevation and
azimuth. Q is the ground's spherical-wave reflection coefficient, which NICE-OPS computes for
each geometry. The table depends on the ground, so it is built for one ground and NICE-OPS
refuses any other. Two tables are shipped in `plates/`:

| table | ground |
| --- | --- |
| `gr1425_db225.nc` | `--flow_resistance 225` (Delany–Bazley grass, the default) |
| `gr1425_vp200.nc` | `--impedance_model variable_porosity --flow_resistance 200` (the fitted Amedee ground) |

Build another with panam's `axisymmetric_bem.table(bands, …)` and `write_netcdf`. A table at
the databases' 31 bands takes about 20 s on 14 cores. Turbulence and roughness do not apply to
the plate, because its reflection is millimetres from its direct path.

Q is taken as quadratic in log frequency across each band, through the band's edges and
centre. This is within 0.07 dB of evaluating Q at every sub-frequency wherever the level is
above −15 dB. The exact per-sub-frequency form (`PlateReceiver::level_exact`) matches panam's
`board_level` to 1e-6 dB in `tests/test_plate_receiver.cpp`. Its cost on the full demo grid
(S-76D track, B407 fixed-load-factor database, 94.3 M evaluations, 14 threads):

| | acoustic loop |
| --- | --- |
| spectral, no ground | 1.6 s |
| `--ground_effect` | 2.1 s (7.1 s with `--exact_ground`) |
| `--ground_effect --plate_table` | 5.4 s (13.3 s with `--exact_ground`) |

By default both receivers interpolate per-trajectory tables instead of evaluating the model at
every observer and state. The tables are built exactly at their nodes in well under 0.1 s.
- **Point receiver:** the whole ground correction, over source elevation and log range. Level
  histories agree with the exact model within 0.001 dB.
- **Plate receiver:** its reflection coefficients. SEL agrees within 0.0005 dB rms.

The interpolated values never exceed the model's bounds, so the peak search and `--db_cutoff`
stay exact. `--exact_ground` evaluates the model directly. Either way the Faddeeva function uses
fast paths: an asymptotic series for |z| ≥ 6 and a Taylor table below, which agree with
Weideman's series within 1e-11.

The peak search and `--db_cutoff` bound the plate with `PlateReceiver::gain_bounds`, so both
keep working. The `plate_windowed_is_exact` ctest checks that skipped states change nothing.

**It is expensive.** Acoustic loop on the Be407 case (94.3 M evaluations, 26 bands, 14 threads,
default band limits): broadband 0.34 s, `--spectral` 1.44 s, `--spectral --ground_effect`
10.9 s. So carrying the spectrum costs about 4x, and the ground effect about 32x — a complex
reflection coefficient and a Faddeeva evaluation per band per observer per state, none of
which can be hoisted because it all depends on the geometry. The Faddeeva evaluation is about
half of the ground path on its own.

Those three are from before the September 2026 optimization round and all of them are now
lower. Measured against the same build on the S-76D demo track (94.3 M evaluations, 14
threads), over a synthetic single-condition 24-band database built the way
`write_flat_spectral_database` in `tests/test_spectral_propagation.cpp` builds one:

| | before | after |
| --- | --- | --- |
| broadband | 0.32 s | 0.16 s |
| `--spectral` | 1.36 s | 1.16 s |
| `--spectral --ground_effect` | 8.95 s | 7.81 s |

A single condition makes the operating-condition blend degenerate, so those two spectral
figures are optimistic against a real multi-condition database like the Be407 one in the row
above; the ground-effect term, which is most of the last row, does not depend on the condition
count. All three outputs agree with the previous build to 1.4e-12 dB.

`--db_cutoff` works on the spectral path through per-direction bounds on the band-limited,
weighted level (`Database::spectral_direction_bounds`). With `--ground_effect` or
`--plate_table`, each band is raised by the most the ground or the plate can add.

### How much gets evaluated

The sound exposure integral only uses the window from each observer's peak down to 10 dB below
it, which on a shipped case is a little over half the track. The observer loop therefore finds
the peak and walks out to the 10-dB-down points rather than marching every emission time, and
the run logs how much that skipped:

```
Evaluated 66041828 of 94334400 (70.0%); the rest is provably outside every observer's 10 dB window.
```

It is exact, not an approximation — the output is byte-identical to evaluating every state.
The path is skipped, and every state evaluated, with `--spectral`, `-a`, `-i`, `-D` or
`--db_cutoff`; see CLAUDE.md for why each one cannot use it.

### Threads

The observer loop runs on every core by default. Each observer is independent, so the output is
deterministic: runs at `-j 1`, `-j 2` and `-j 32` are byte-for-byte identical.

Observers are dealt out to the workers in chunks as they finish rather than split evenly up
front, because on a machine with both performance and efficiency cores — every Apple silicon
part, and most recent x86 laptops — an even split makes the slowest core the critical path.
On an M4 Pro the even split ran the demo case no faster on all 14 threads than on its 10
performance cores; pulling chunks recovers about 10%.

### Geodetic coordinates

A trajectory can be given in WGS84 coordinates — latitude, longitude and altitude — instead of
the 27-column AC state format. That is what a mission planner has: a route, not a flight-test
record.

```
time,lat,lon,alt
0.0,39.82914323,-74.40996141,1134.60
0.5,39.82903000,-74.41007000,1133.80
...
```

Both formats are read by `-t`, and which one a file is is decided by its width: 27 columns is
the positional AC format, anything else is read by the names in its header. Column order does
not matter and unrecognised headers are an error rather than a guess.

| Column | Meaning |
| --- | --- |
| `time` | Seconds from the start of the track. `utcsec` (seconds since UTC midnight) serves instead, and sets the clock for `-D`. |
| `lat`, `lon` | Degrees. |
| `alt` | **Feet above mean sea level.** `--ground_elevation` gives the elevation of the ground the observers stand on, which is subtracted to get height above them. |
| `vE`, `vN`, `vUp` | Optional, ft/s. Used when all three are present; otherwise velocity is differentiated from the path. |

Names are matched case-insensitively, ignoring punctuation and a parenthesised unit, so
`Latitude`, `LAT` and `lat (deg)` are one column. `alt_ft` is recognised and `alt_m` is
deliberately *not*: a file naming metres fails to load rather than being read as feet.

**Where the frame sits.** `--origin <lat> <lon>` fixes the local x/y origin; without it, the
first position of the first geodetic track is used. `--frame_bearing` turns the axes so +x runs
along a compass bearing — the default leaves them east/north, which is what a map wants. One
frame serves the whole run, because every trajectory in it shares the observers.

**What is derived, and how well.** A WGS84 track gives position only, so the velocity and the
load factor behind C<sub>T</sub> come from the first and second derivatives of a smoothing
spline through the path. `--path_smoothing` says how far, in feet RMS **per axis**, that fit
may sit from the given positions. It is an acceptable deviation, not an abstract weight: the
solver finds the smoothest curve whose RMS residual meets the bound, so asking for 10 ft means
the fit will use 10 ft of slack, not merely stay under it.

The default is 10 ft. Zero interpolates exactly, and exact interpolation is the wrong thing for
both kinds of input: measured positions make the second derivative track coordinate rounding,
and sparse waypoints are a statement of intent rather than a constraint to be met to the inch —
flying precisely through them manufactures aggressive manoeuvring, which is itself a noise
source.

Feeding the demo track's own lat/lon/alt columns back in, against the footprint computed from
its full inertial state:

| Input | Median | 99th pct | Max | Cells over 1 dB (of 14400) |
| --- | --- | --- | --- | --- |
| 50 Hz positions, default `--path_smoothing 10` | 0.085 dB | 1.00 | 2.88 | 145 |
| 50 Hz positions, `--path_smoothing 0.0656` (= 0.02 m) | 0.028 dB | 0.64 | 2.44 | 37 |
| 50 Hz positions, `--path_smoothing 0` (exact) | 11.6 dB | 15.51 | 16.32 | 14395 |
| 50 Hz positions and velocity columns, default | 0.014 dB | 0.42 | 2.25 | 18 |
| Waypoints 10 s apart, resampled to 0.1 s, default | 0.096 dB | 0.95 | 2.10 | 107 |

The projection itself is not the lossy part: it reproduces that file's own x/y to under a metre
over 3.2 km, and altitude passes through exactly. What costs is differentiating twice.

The default is not the best row in the table — a dense 50 Hz track does better at 0.0656 ft,
which is tuned to *that* file's coordinate resolution. It is the best *default*, because it is
the only setting that is good on every input shape at once: 0.085 dB on dense positions,
0.096 dB on sparse waypoints, 0.014 dB when velocities are given. Exact interpolation is
uniformly catastrophic on measured input, and a value tuned for one sampling rate is wrong for
another. Tune it down if you know your coordinate resolution; leave it alone otherwise.

Sparse waypoints are resampled onto the `-s` grid first — never coarsened — and their load
factor is an artefact of the interpolant, which is logged.

> The `--path_smoothing 0` row was, until recently, reported as 0.196 dB. That measured a bug:
> `SmoothingSpline` short-circuited at λ = 0 and returned zero curvature, so the load factor
> came out as exactly 1 g on every state and the rounding amplification never happened. With
> the interpolating spline carrying its own curvature, exact interpolation costs what the
> surrounding text always said it would.

**Note on velocity columns.** `--path_smoothing` is a distance and applies to the fitted
position only. A file carrying `ve`/`vn`/`vu` has those interpolated exactly; they are
measurements of the quantity wanted, not something to be re-derived.

**Output.** `--latlon` writes a companion beside every level file (`out.csv` →
`out_latlon.csv`) holding latitude, longitude and level, one row per observer, for dropping a
footprint onto a map. It needs a frame: either a geodetic track, which carries one, or
`--origin`. A grid regular in the local frame is not regular in latitude and longitude, so this
file is always scattered, never gridded.

**The transform.** `src/geodetic.hpp` is a transverse Mercator about the reference point, with
altitude carried through as z — a map projection with a height attached rather than a true ENU
frame. Over the few kilometres a footprint spans the difference is far below the modelling
error. `from_east_north_up` applies the frame's rotation without its translation, for
velocities; `to_geodetic` inverts the projection by Newton iteration, and is what `--latlon`
uses.

The shipped AC tracks are the case where the axes are not east/north: their x/y run along a
compass bearing of **210°**, the approach heading, not along the wandering heading in the
file's own column 8. Their ground plane sits 54.2651 ft below their altitude column. So
```sh
./build/niceops -t track.csv --origin 39.8104276 -74.4263086 --frame_bearing 210 \
                --ground_elevation 54.2651 --path_smoothing 0.02
```
is the demo case restated in WGS84, and is what the table above measures.

### Plotting

Plotting embeds Python through pybind11. The first configure creates a virtual environment at
`build/pyenv` and installs matplotlib into it, and the embedded interpreter puts that directory
at the front of `sys.path` — so a run depends on what the build made rather than on whatever
happens to be installed for the system interpreter, and needs no `PYTHONPATH`. The install
happens once per build directory and needs network access; afterwards configures are offline.

Two options control this:

| Option | Meaning |
| --- | --- |
| `-DNICEOPS_PLOTTING=OFF` | Drop plotting, and Python, pybind11 and matplotlib with it. The acoustics need none of them. `-p`, `-f` and `-F` are then refused with a message rather than silently doing nothing. |
| `-DNICEOPS_BUNDLE_PYTHON=OFF` | Keep plotting but skip the virtual environment, using the interpreter's own packages instead. |

Headless is handled: the facade selects the `Agg` backend itself.

The reference figures under `figures/` date from June 2020 and were produced with three settings
that differ from the defaults. All three are now reachable:

```sh
./build/niceops -t tracks/S-76D/197401AC.csv -v databases/S-76D.json -d databases/S-76D_M3.nod \
        --trim -6000 --grid -2000 1100 -1500 1600 -o out.csv -f out.eps
```

Matched that way the two runs land on the same 124 × 124 grid and agree cell for cell: the largest
difference anywhere is **1.05 dB**, 99% of cells are within 0.5 dB, and the peak levels differ by
0.07 dB. Contour levels now default to 65–100, the values the references were drawn with; the old
45–90 put its top level below the peak of the demo case, leaving the loudest patch unfilled.

## Licence

NICE-OPS is under the MIT License; see `LICENSE`. Every third-party component it
builds against is permissive — BSD, MIT or the Boost Software License — and none carries a
reciprocal or source-disclosure obligation. `THIRD_PARTY_NOTICES.md` reproduces their notices,
which a binary distribution needs to carry with it.

## Repository layout

- `src/` — the implementation. Every `.cpp` lives here.
- `tests/` — Catch2 unit tests and the `golden_compare` tool.
- `golden/` — reference outputs captured from the D implementation; the regression tests compare
  against these at 1e-9 relative tolerance.
- `databases/`, `tracks/` — shipped noise databases and flight tracks.

`CLAUDE.md` carries the developer-facing detail: the port's history, the shape of the hot path,
and the pitfalls worth knowing before changing the numerics.
