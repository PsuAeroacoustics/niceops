# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **O**perations **P**lanning **S**ystem is a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

For each state along the trajectory, NICE-OPS blends the database spheres that bracket the
aircraft's operating condition (advance ratio, flight path angle and thrust coefficient). It
reads the level off the blend in the direction of each observer, then carries that level to the
ground:

- A broadband database uses its own excess-attenuation factor.
- A database with spectra is propagated band by band, with ISO 9613-1 absorption, ground
  reflection and, optionally, refraction through a stratified, windy atmosphere.

Each observer's level history is integrated into a sound exposure level (SEL), and a set of
operations into a day-night level (DNL). Lden, CNEL, LAFmax, NAx, time above, C-weighted SEL and
EPNL are available on request.

NICE-OPS was developed by Eric Greenwood and Rob Rau. It was written in D and is now
implemented in C++17. [docs/HISTORY.md](docs/HISTORY.md) covers the port and the changes that
moved results.

This file covers building NICE-OPS and running a first case. **[The manual](docs/MANUAL.md)**
covers every option, the reasoning behind each model, and how each model was checked.
**[The theory](docs/THEORY.md)** gives the equations and the references they come from.

## Building

The requirements are a C++17 compiler and netCDF. Plotting also needs Python 3 with
development headers and the `venv` module. The build installs matplotlib itself. On Ubuntu:

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
| `-DNICEOPS_ASSERTS=ON` | Keep `assert()` live in an optimized build. The default Release build compiles the asserts out. |
| `-DNICEOPS_LTO=OFF` | Turn off link-time optimization, which is otherwise on wherever the toolchain supports it. |

On macOS with Command Line Tools but no full Xcode, add
`-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`.

### Tests

Test with the asserts live, as CI does. The suite runs from the source root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON
cmake --build build -j8
ctest --test-dir build --output-on-failure
```

Most of the roughly 260 tests run `niceops` and compare its output against `golden/` or against
another run. [docs/internals/testing.md](docs/internals/testing.md) covers the goldens and
what each kind of test checks.

## Quick start

With no arguments, NICE-OPS runs a demo case. It runs the trajectory `tracks/S-76D/197421AC.csv`
against the database `databases/S-76D_M3.nod` and the vehicle `databases/S-76D.json`, and
writes `integrated/197421_LAE.csv`.

```sh
./build/niceops                                     # the demo case
./build/niceops -t my_track.csv -d my_db.nod -v my_vehicle.json -o out.csv
./build/niceops --help                              # the full option list
```

Paths are relative to the working directory, so run from the repository root when using the
shipped data. Output directories are created as needed.

The S-76D demo database is broadband. A database that carries spectra runs the spectral path,
with ground effect, by default. That path needs a relative humidity, which comes from the
vehicle file, `--humidity` or `--atmosphere`:

```sh
./build/niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json \
                --humidity 70 -o integrated/197421_B407.csv
```

Each run prints a one-line `Physics:` summary saying which path, ground model and atmosphere it
used.

### Inputs

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV, with positions and velocities in feet. A WGS84 track (`time`, `lat`, `lon`, `alt` columns) or a list of static holds also works. See the manual's [Geodetic coordinates](docs/MANUAL.md#geodetic-coordinates) and [Input files](docs/MANUAL.md#input-files). |
| Database (`-d`) | netCDF `.nod`, built by panam. It holds one group per operating condition (advance ratio μ, flight path angle γ, thrust coefficient C<sub>T</sub>), and each group holds a sphere of levels. |
| Vehicle (`-v`) | JSON with `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg) and `Atmosphere.density`. A vehicle file that disagrees with the database by more than 2 percent is refused. |
| Observers (`-b`) | CSV of scattered observer positions: `x,y,z`, `id,x,y,z`, or `id,x,y,z,height`, where `height` is each microphone's height above the ground (see [Receivers](docs/MANUAL.md#receivers)). The first line is a header. Without `-b`, NICE-OPS uses a 120×120 grid spanning −1500…1475 m. |
| Terrain (`--terrain`) | netCDF elevation grid (e.g. USGS 3DEP). See [Terrain](docs/MANUAL.md#terrain). |

### Outputs

- `-o`: the SEL per observer. Gridded output writes the x and y axes on the first two lines,
  then one row of levels per y.
- `-D`, `--lden_file`, `--cnel_file`: the day-night, day-evening-night and community noise
  equivalent levels, summed over every trajectory in the run.
- `--lamax_file`, `--nax_file`, `--ta_file`, `--csel_file`: LAFmax, NAx, time above a
  threshold, and the C-weighted SEL.
- `-i`, `--band_history`: per-observer level and band time histories.
- `--epnl_file`, `--band_frames`: per-observer EPNL (14 CFR 36 A36.4), and, as netCDF, the band
  levels in the 0.5 s reception-time frames it is computed from.
- `-S`: the derived trajectory state (position, velocity, load factor, C<sub>T</sub>, μ, flight
  path angles).
- `--latlon`: a latitude, longitude and level companion to every level file.
- `-p`, `-f`, `-F`: show the footprint plot, save it to a file, or save the state plots.

The manual's [Output files](docs/MANUAL.md#output-files) has the full list.

### Units

Input and output files use **feet** by default. The computation is metric internally.
`--units m` switches every length on the command line and in every file to meters. It declares
what your files already are and converts nothing, so the shipped tracks, which are in feet,
must be run without it. The `-S` state file is always metric. Speeds name their own unit
(`12kt`, `6.2mps`). See [Units](docs/MANUAL.md#units).

### Options to know first

| Option | Meaning |
| --- | --- |
| `-t`, `-d`, `-v` | Trajectory, database and vehicle. `-t` may be repeated, with either one `-d` and one `-v`, or one of each per track. |
| `-o`, `-b` | The SEL output file, and a scattered observer file to use instead of the default grid. |
| `--grid`, `--gridstep` | Observer grid bounds and spacing, in `--units`. |
| `--frame_bearing` | Compass bearing of the local +x axis. The default is x east, y north. The shipped AC tracks need 210. It matters for wind, WGS84 output, `--nose heading` and propagation tables. |
| `--spectral` / `--no-spectral` | Per-band propagation. On by default when the database carries spectra. |
| `--ground_effect` / `--no-ground_effect` | Ground reflection. On by default wherever the spectral path runs. |
| `--humidity`, `--temperature`, `--pressure` | The atmosphere for the spectral path, overriding the vehicle file. |
| `--atmosphere`, `--refraction` | A stratified atmosphere profile, and curved rays through it. |
| `--wind` | A profile built from one surface wind (`12kt,270`), used instead of an `--atmosphere` file. |
| `--terrain` | A terrain grid. Observers stand on it, and the ground model takes each path's mean ground plane. |
| `-j` | Worker threads. 0 (the default) uses every core. The output is identical at any count. |

[Command-line options](docs/MANUAL.md#command-line-options) in the manual lists every option.

## Shipped data

### Databases and vehicle files (`databases/`)

| Files | Aircraft | Contents |
| --- | --- | --- |
| `*_fixed_load_factor.nod` | AS350, B407, Be206L3, EC130, RO-44, RO-66 | **Spectral.** Built by panam from the 2017 Noise Abatement flight test. These are the current, validated databases. |
| `*_ambient_gated.nod` | the same six | Broadband versions of the same builds. |
| `*_fixed_load_factor.nod` | AW139, B205, S-76D | Broadband, converted from the `_M*` files below. |
| `*_M1`, `_M2`, `_M3`, `_M4.nod` | AW139, B205, B407, S-76D | Broadband, discrete-load-factor databases from the 2020–21 joint NASA/FAA flight tests. `S-76D_M3.nod` is the demo database and the one the goldens use. |
| `*.json` | AW139, B205, B407, S-76D | Vehicle files. |

Only four vehicle files ship. The vehicle files for the AS350, Be206L3, EC130, RO-44 and RO-66
are kept with the validation harness (`approach/vehicles/`). The manual's
[Noise databases](docs/MANUAL.md#noise-databases) section covers how the databases are built and
how well each one validates.

### Other data

| Directory | Contents |
| --- | --- |
| `tracks/` | Flight tracks in the AC state format: 30 S-76D, 41 AW139 and 13 B205. |
| `arrays/` | `coyle_approach.csv`, a scattered microphone array. |
| `atmospheres/` | Profiles for `--atmosphere`: the ISA (`isa.csv`), a windy log law (`example_windy.csv`, plus `example_windy_calm_floor.csv`, with a 0.1 ft calm floor, for a test), a morning inversion (`example_inversion.csv`) and B407's uniform build atmosphere (`B407_uniform.csv`). |
| `plates/` | Ground-plate microphone tables for `--plate_table`. |

## Documentation

- [`docs/MANUAL.md`](docs/MANUAL.md): how a prediction is made, the full option reference, the
  propagation models, their checks against exact and full-wave solutions, and performance.
- [`docs/THEORY.md`](docs/THEORY.md): the models as equations, from the source spheres through
  propagation to the metrics, with references.
- [`docs/FIDELITY.md`](docs/FIDELITY.md): fidelity status and validation against measurement,
  experimental source improvements, and propagation-table error checks.
- [`docs/HISTORY.md`](docs/HISTORY.md): the D original, the C++ port, the departures from D,
  and where to find the changes that moved results. **Results from before September 2026 are
  not comparable with results from after.**
- [`CLAUDE.md`](CLAUDE.md) and [`docs/internals/`](docs/internals/): developer notes. They
  cover the build, the shape of the hot path, the measurements behind each module, and the
  pitfalls to know before changing the numerics.
- `reports/` and `research_notes/`: literature reviews (grazing incidence, hemisphere models,
  urban propagation), the propagation roadmap, and the 2026-10-02 code review.

## Repository layout

- `src/`: the implementation. Every `.cpp` lives here.
- `tests/`: Catch2 unit tests and the `golden_compare` tool.
- `golden/`: reference outputs. The regression tests compare against these at a relative
  tolerance of 1e-9.
- `tools/`: developer programs. The C++ ones are built alongside `niceops`:
  - `niceops_pe_compare` and `niceops_pe_field` check the refraction models against the
    parabolic equation, and `niceops_ray_geometry` traces the same rays outside `niceops`.
    `plot_pe_compare.py` plots the comparison.
  - `niceops_benchmark` runs the Attenborough et al. benchmark cases, and `plot_benchmark.py`
    plots them.
  - `niceops_propagation_table` builds `--propagation_table` files.
  - `pnl_reference.py` generates the PNL/EPNL test references from panam.
  - `score_fidelity_frames.py`, `build_hover_candidate.py`, `export_emission_samples.py` and
    `regrid_emission_samples.py` are experimental source-model studies, run against the
    validation harness (see [FIDELITY.md](docs/FIDELITY.md)).
- `*batch.sh`: batch scripts from 2020. Each one runs every shipped track of one aircraft.
- `state/`, `integrated/`, `log/`, `figures/`: run outputs. These directories are not tracked,
  and each run creates what it needs.

## License

NICE-OPS is under the MIT License; see `LICENSE`. Every third-party component it builds against
is permissive (BSD, MIT or the Boost Software License), and none carries a reciprocal or
source-disclosure obligation. `THIRD_PARTY_NOTICES.md` reproduces their notices, which a binary
distribution must carry with it.
