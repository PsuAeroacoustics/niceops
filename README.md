# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **O**perations **P**lanning **S**ystem — a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

NICE-OPS was developed by Eric Greenwood and Rob Rau. It began life in D; this is the C++17
implementation, which reproduces the D original's results on the demo case to within a few
units in the last place.

For each state along the trajectory, NICE-OPS blends the database spheres that bracket the
aircraft's operating condition (advance ratio, flight path angle and thrust coefficient), and
reads off the level in the direction of each observer. It then carries that level to the
ground. A broadband database uses its own attenuation factor. A database with spectra is
propagated band by band, with ISO 9613-1 absorption, ground reflection and, optionally,
refraction through a stratified, windy atmosphere. Each observer's level history is integrated
into a sound exposure level (SEL), and a set of operations into a day-night level (DNL).

This file covers building NICE-OPS and running a first case. **[The manual](docs/MANUAL.md)**
covers every option, the theory behind each model, and how each model was checked.

## Building

Requirements: a C++17 compiler and netCDF. Plotting additionally needs Python 3 with
development headers and the `venv` module; matplotlib itself is installed by the build, not by
you. On Ubuntu:

```sh
sudo apt-get install libnetcdf-dev python3-dev python3-venv
```

Build out of source, always with explicit `-S`/`-B`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

CLI11, nlohmann_json, Catch2 and pybind11 are fetched at configure time, so the first configure
needs network access.

| CMake option | Effect |
| --- | --- |
| `-DNICEOPS_PLOTTING=OFF` | Drop plotting, and with it the Python requirement. The acoustics do not need it. |
| `-DNICEOPS_BUNDLE_PYTHON=OFF` | Keep plotting, but use the system interpreter's packages instead of the build's own environment. |
| `-DNICEOPS_ASSERTS=ON` | Keep `assert()` live in an optimized build, as CI does. The default Release build compiles the asserts out. |
| `-DNICEOPS_LTO=OFF` | Turn off link-time optimization, which is otherwise on wherever the toolchain supports it. |

On macOS with Command Line Tools but no full Xcode, add
`-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`.

Run the tests with `ctest --test-dir build -C Release`.

## Quick start

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

The S-76D demo database is broadband. A database that carries spectra runs the spectral path,
with ground effect, by default. That path needs a relative humidity, from the vehicle file,
`--humidity` or `--atmosphere`:

```sh
./build/niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json \
                --humidity 70 -o integrated/197421_B407.csv
```

### Inputs

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV, positions and velocities in feet. A WGS84 track (`time`, `lat`, `lon`, `alt` columns) or a list of static holds also works; see the manual's [Geodetic coordinates](docs/MANUAL.md#geodetic-coordinates) and [Input files](docs/MANUAL.md#input-files). |
| Database (`-d`) | netCDF `.nod`: one group per operating condition (advance ratio μ, flight path angle γ, thrust coefficient C<sub>T</sub>), each holding a sphere of levels. |
| Vehicle (`-v`) | JSON with `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg) and `Atmosphere.density`. |
| Observers (`-b`) | CSV of scattered observer positions: `x,y,z`, `id,x,y,z`, or `id,x,y,z,height` with each microphone's height above the ground (see [Receivers](docs/MANUAL.md#receivers)). The first line is a header. Without `-b`, a 120×120 grid spanning −1500…1475 m is used. |

### Outputs

- `-o`: the SEL per observer. Gridded output writes the x and y axes on the first two lines, then
  one row of levels per y.
- `-D`, `--lden_file`, `--cnel_file`: the day-night, day-evening-night and community noise
  equivalent levels, summed over every trajectory in the run.
- `--lamax_file`, `--nax_file`, `--ta_file`, `--csel_file`: LAFmax, NAx, time above a
  threshold and the C-weighted SEL.
- `-i`, `--band_history`: per-observer level and band time histories.
- `-S`: the derived trajectory state (position, velocity, load factor, C<sub>T</sub>, μ, flight
  path angles).
- `--latlon`: a latitude, longitude, level companion to every level file.

The manual's [Output files](docs/MANUAL.md#output-files) has the full list.

### Units

Input and output files use **feet** by default, and the computation is metric internally.
`--units m` switches every length on the command line and in every file to meters. It
declares what your files already are and converts nothing, so the shipped tracks, which are in
feet, must be run without it. The `-S` state file is always metric. See
[Units](docs/MANUAL.md#units).

### Options to know first

| Option | Meaning |
| --- | --- |
| `-t`, `-d`, `-v` | Trajectory, database and vehicle. `-t` may be repeated, with one or one-per-track `-d` and `-v`. |
| `-o`, `-b` | SEL output file, and a scattered observer file instead of the default grid. |
| `--grid`, `--gridstep` | Observer grid bounds and spacing, in `--units`. |
| `--frame_bearing` | Compass bearing of the local +x axis. The shipped AC tracks use 210; it matters for wind, WGS84 output and propagation tables. |
| `--spectral` / `--no-spectral` | Per-band propagation, on by default when the database carries spectra. |
| `--ground_effect` / `--no-ground_effect` | Ground reflection on the spectral path, on by default there. |
| `--atmosphere`, `--refraction` | A stratified atmosphere profile, and curved rays through it. |
| `-j` | Worker threads. 0 (the default) uses every core; output is identical at any count. |
| `-p`, `-f` | Plot the footprint, or save it to a file. |

[Command-line options](docs/MANUAL.md#command-line-options) in the manual lists every option.

## Shipped data

| Directory | Contents |
| --- | --- |
| `databases/` | Noise databases and vehicle files. `S-76D_M3.nod` is the demo database. The AS350, B407, Be206L3, EC130, RO-44 and RO-66 `_fixed_load_factor` databases carry spectra, built from the 2017 Noise Abatement flight test. |
| `tracks/` | Flight tracks in the AC state format, for the S-76D, AW139 and B205. |
| `arrays/` | `coyle_approach.csv`, a scattered microphone array. |
| `atmospheres/` | Profiles for `--atmosphere`: the ISA (`isa.csv`), a windy log law (`example_windy.csv`), a morning inversion (`example_inversion.csv`) and B407's own uniform atmosphere. |
| `plates/` | Ground-plate microphone tables for `--plate_table`. |

## Documentation

- [`docs/MANUAL.md`](docs/MANUAL.md): how a prediction is made, the full option reference, the
  propagation models, their checks against exact and full-wave solutions, and performance.
- [`CLAUDE.md`](CLAUDE.md): developer notes on the port's history, the shape of the hot path,
  and the pitfalls to know before changing the numerics.
- [`docs/FIDELITY.md`](docs/FIDELITY.md): fidelity status and validation, with current
  measurement evidence, experimental source improvements and propagation-table error checks.

**Results from before September 2026 are not comparable to results after.** See
[Changes that affect results](docs/MANUAL.md#changes-that-affect-results).

## Repository layout

- `src/`: the implementation. Every `.cpp` lives here.
- `tests/`: Catch2 unit tests and the `golden_compare` tool.
- `golden/`: reference outputs. The regression tests compare against these at 1e-9 relative
  tolerance.
- `tools/`: `pe_compare.cpp`, the refraction models against the parabolic equation;
  `benchmark.cpp`, the Attenborough et al. benchmark cases; `propagation_table.cpp`, which
  builds `--propagation_table` files; with plotting scripts for the first two.

## License

NICE-OPS is under the MIT License; see `LICENSE`. Every third-party component it
builds against is permissive — BSD, MIT or the Boost Software License — and none carries a
reciprocal or source-disclosure obligation. `THIRD_PARTY_NOTICES.md` reproduces their notices,
which a binary distribution needs to carry with it.
