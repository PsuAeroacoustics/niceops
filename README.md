# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **OP**eration Planning **S**ystem — a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

NICE-OPS was developed by Rob Rau and Eric Greenwood. It began life in D; this is the C++17
implementation, which reproduces the D original's results on the demo case to within a few
units in the last place.

## Building

Requirements: a C++17 compiler, netCDF, CGAL, Eigen3, Boost, and Python 3 with development
headers plus matplotlib (pybind11 embeds the interpreter for plotting). On Ubuntu:

```sh
sudo apt-get install libcgal-dev libnetcdf-dev libeigen3-dev python3-dev python3-matplotlib
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
- `--export_sphere` — a whole interpolated noise sphere, written as a one-group `.nod`.

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
| `-l`, `--logfile` | Log file (default `log/niceops.log`). |
| `-L`, `--loglevel` | `trace`, `info`, `warn`, `error` or `fatal` (default `info`). |

### Physics and metrics

| Option | Meaning |
| --- | --- |
| `-a`, `--ambient` | Ambient level in dB, mixed into every sample. Off unless given. |
| `-q`, `--qsam` / `--no-qsam` | Use the flight path angle relative to the tip-path-plane normal rather than the vertical. On by default. |
| `-s`, `--step` | Resampling step for trajectory adaptation. **Not** the SEL integration step, which comes from the trajectory's own sampling. |
| `-D`, `--dnlfile` | Write the day-night level, summed over every trajectory. |
| `--utc_offset` | Hours from UTC to local time, for the DNL night window. |
| `--powerfile` | Write per-state shaft power, when the database has it. |

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
| `-p`, `--plotting` | Plot the footprint; `--cmin`/`--cmax`/`--cstep` set the contour levels. |
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

### Threads

The observer loop runs on every core by default. Each observer is independent, so the output is
deterministic: runs at `-j 1`, `-j 2` and `-j 32` are byte-for-byte identical.

## Repository layout

- `src/` — the implementation. Every `.cpp` lives here.
- `tests/` — Catch2 unit tests and the `golden_compare` tool.
- `golden/` — reference outputs captured from the D implementation; the regression tests compare
  against these at 1e-9 relative tolerance.
- `databases/`, `tracks/` — shipped noise databases and flight tracks.

`CLAUDE.md` carries the developer-facing detail: the port's history, the shape of the hot path,
and the pitfalls worth knowing before changing the numerics.
