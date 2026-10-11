# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **O**perations **P**lanning **S**ystem is a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

1. Blend the database spheres that bracket each state's operating condition (advance ratio,
   flight path angle, thrust coefficient), and read the level toward each observer.
2. Carry it to the ground: a broadband database uses its own excess-attenuation factor; a
   database with spectra is propagated band by band, with ISO 9613-1 absorption, ground
   reflection and, optionally, refraction through a stratified, windy atmosphere.
3. Integrate each observer's level history into a sound exposure level (SEL), and a set of
   operations into a day-night level (DNL). Lden, CNEL, LAmax, NAx, time above, C-weighted SEL
   and EPNL are available on request.

NICE-OPS was developed by Eric Greenwood and Rob Rau, written first in D and now implemented in
C++17; [docs/HISTORY.md](docs/HISTORY.md) covers the port and every change that moved results.

## Building

The requirements are a C++17 compiler and netCDF. Plotting also needs Python 3 with
development headers and the `venv` module. The build installs matplotlib itself. On Ubuntu:

```sh
sudo apt-get install libnetcdf-dev python3-dev python3-venv
```

Build out of source, always with explicit `-S`/`-B`. CLI11, nlohmann_json, Catch2 and pybind11
are fetched at configure time, so the first configure needs network access.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

| CMake option | Effect |
| --- | --- |
| `-DNICEOPS_PLOTTING=OFF` | Drop plotting, and with it the Python requirement. The acoustics do not need it. |
| `-DNICEOPS_BUNDLE_PYTHON=OFF` | Keep plotting, but use the system interpreter's packages instead of the build's own environment. |
| `-DNICEOPS_ASSERTS=ON` | Keep `assert()` live in an optimized build. The default Release build compiles the asserts out. |
| `-DNICEOPS_LTO=OFF` | Turn off link-time optimization, which is otherwise on wherever the toolchain supports it. |

On macOS with Command Line Tools but no full Xcode, add
`-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`. Run the
tests, from a build with `-DNICEOPS_ASSERTS=ON` as CI does: `ctest --test-dir build`.

## Quick start

With no arguments, NICE-OPS runs a demo case: the trajectory `tracks/S-76D/197421AC.csv`, the
database `databases/S-76D_M3.nod` and the vehicle `databases/S-76D.json`, writing
`integrated/197421_LAE.csv`. That database is broadband. A database that carries spectra runs
the spectral path, with ground effect, by default; it needs a relative humidity, from the
vehicle file, `--humidity` or `--atmosphere`.

```sh
./build/niceops                                     # the demo case
./build/niceops -t my_track.csv -d my_db.nod -v my_vehicle.json -o out.csv
./build/niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json \
                --humidity 70 -o integrated/197421_B407.csv    # spectral
./build/niceops --help                              # the full option list
```

`-t`, `-d` and `-v` give the trajectory, database and vehicle; `-o` is the SEL output file, and
`-D` adds the day-night level over every trajectory, or every [operations file](docs/MANUAL.md#an-operations-file)
row. Paths are relative to the working directory: run from the repository root to use the shipped
data. Each run prints a one-line `Physics:` summary of the path, ground model and atmosphere. The
manual lists every [option](docs/MANUAL.md#command-line-options), [input file](docs/MANUAL.md#input-files)
and [output file](docs/MANUAL.md#output-files).

### Units

Files and command-line lengths are in **feet** by default; the computation is metric inside.
`--units m` declares that your files are already in meters and converts nothing, so the shipped
tracks, which are in feet, must be run without it. The `-S` state file is always metric, and
speeds name their own unit (`12kt`, `6.2mps`). See [Units](docs/MANUAL.md#units).

## Shipped data

`tracks/` holds AC-format flight tracks (30 S-76D, 41 AW139, 13 B205). `databases/` holds the
sphere databases and four vehicle files; see [Shipped databases](docs/MANUAL.md#shipped-databases).
`atmospheres/` holds `--atmosphere` profiles, `arrays/` a scattered microphone array, and
`plates/` ground-plate tables for `--plate_table`.

## Documentation

- [`docs/MANUAL.md`](docs/MANUAL.md): using NICE-OPS: every option, input, output and model.
- [`docs/THEORY.md`](docs/THEORY.md): the models as equations, with references.
- [`docs/FIDELITY.md`](docs/FIDELITY.md): how accurate each model is, and the evidence.
- [`docs/HISTORY.md`](docs/HISTORY.md): the D original, the port, and the
  [changes that affect results](docs/HISTORY.md#changes-that-affect-results).
- [`docs/internals/README.md`](docs/internals/README.md): developer documentation.
- [`CLAUDE.md`](CLAUDE.md): rules for AI coding agents working in this repository.
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md): the notices of the fetched dependencies.

## License

NICE-OPS is under the MIT License; see [`LICENSE`](LICENSE). Every third-party component it
builds against is permissive (BSD, MIT or the Boost Software License), and none carries a
reciprocal or source-disclosure obligation. A binary distribution must carry their notices,
which `THIRD_PARTY_NOTICES.md` reproduces.
