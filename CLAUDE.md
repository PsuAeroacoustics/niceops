# CLAUDE.md

This file gives Claude Code (claude.ai/code) guidance for working with code in this repository.
It holds what applies everywhere. The measurements, the history, and the per-module "things
to preserve" are in **`docs/internals/`** (see the map below). **Read the matching file before
you change a module**: most of what looks over-engineered there was measured, and many of the
simplifications that look obvious were tried and failed.

## Project Overview

NICE-OPS (Noise Informed Community Environment Operation Planning System) is a helicopter
noise footprint tool. It predicts ground-level dBA noise exposure over an observer grid (or
scattered observer set) from a flight trajectory and a rotor-noise "sphere" database (`.nod`,
built by panam).

The implementation is C++17 and lives in `src/`. It began as a port of a D original, which is
why several comments reference D semantics. Those comments are load-bearing: the port reproduces
the D reference bit-for-bit on the demo case, and the goldens come from it. The D sources were
removed once the port was verified complete. They remain in history at `c74bab3` (on
`origin/no-stk`), which is what to check out if you need the oracle again.

- **Every `.cpp` lives in `src/` and is listed in CMakeLists.txt.** There are no C++ sources at
  the repo root. A new source has to go into the right target: a file that compiles but is never
  linked shows up as an undefined-symbol error far from its cause.
- `README.md` is short: a description, building, and a first run. **`docs/MANUAL.md`** is the
  user manual: every option, the theory behind each model, and how each was checked.
  User-facing detail goes there, not in the README. Developer detail goes in `docs/internals/`.
  Project history (the D original, the port, changes that moved results) goes in
  `docs/HISTORY.md`.

## Build

The requirements are a C++17 compiler and netCDF (shared library). Plotting also needs Python 3
**with development headers**; the build installs matplotlib into `build/pyenv`. On Ubuntu, as
CI does:

```sh
sudo apt-get install libnetcdf-dev python3-dev python3-venv
```

`python3-venv` is easy to miss. Without it the build's environment step only warns and falls
back to the system interpreter's packages, so CI checks the environment explicitly. CLI11,
nlohmann_json, Catch2 and pybind11 come in through `FetchContent`, so the first configure needs
network access.

Build out-of-source, **always with explicit `-S`/`-B`**:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON
cmake --build build -j8
```

Options:

- `-DNICEOPS_ASSERTS=ON` makes asserts **live when testing** by stripping `-DNDEBUG`, and
  `tests/test_build.cpp` fails if NDEBUG survives. A plain Release build is a shipped build
  with asserts off. Always test with asserts on. Every assert is side-effect free.
- `-DNICEOPS_PLOTTING=OFF` drops Python entirely (CI builds this too: it is the only build that
  compiles `src/plotting_disabled.cpp`). `-DNICEOPS_BUNDLE_PYTHON=OFF` keeps plotting but uses
  the interpreter's own packages.
- `-DNICEOPS_LTO=OFF` turns off link-time optimization, which is on by default and worth 18% of
  the acoustic loop.
- With Command Line Tools only (no Xcode), pass `-DCMAKE_CXX_COMPILER=/usr/bin/clang++
  -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`.

**CI builds on GCC/libstdc++, and a macOS clang build won't catch what GCC finds.** That
includes missing transitive headers (`<algorithm>`) and static-library link cycles, which ld64
resolves and GNU ld doesn't. Before you push anything structural, run the build and tests under
Homebrew GCC:

```sh
cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON \
      -DCMAKE_CXX_COMPILER=$(brew --prefix)/bin/g++-16 -DNICEOPS_PLOTTING=OFF
cmake --build build-gcc -j8 && ctest --test-dir build-gcc --output-on-failure
```

- **Keep the C++ sources in one static library.** A separate target that calls back into
  `niceops_core` makes a cycle, which breaks Linux only.
- **`state/`, `integrated/`, `log/` and `figures/` are not tracked.** Every writer must create
  its own parent directory (`ensure_parent_directory`, `src/util/paths.hpp`). A working tree is
  therefore not a faithful test of a checkout. To verify one, use a clone (`git clone
  --no-hardlinks . /tmp/fresh`, then build and run ctest there).
- **Don't assert bit-identity against a libm call on a literal argument** in a test. GCC
  constant-folds the literal through MPFR, and the two results differ by an ULP.

Details and history: [docs/internals/build.md](docs/internals/build.md).

## Tests

```sh
ctest --test-dir build --output-on-failure      # what CI runs: the whole suite, no -E
ctest --test-dir build -R unit_test -V          # one CTest test
build/niceops_unit_tests "[sphere_database]"    # Catch2 by tag (much faster than CTest)
build/niceops_unit_tests "lookup_angles"        # by test name
```

There are about 260 CTests (`ctest -N`). Most of them are `niceops` runs compared against
another run, and all of them run with the **source root** as their working directory. The main
ones:

- `regression_default` runs the demo case, and its output feeds `regression_golden` and
  `state_golden`. These two use `golden_compare` (elementwise, 1e-9 relative, header lines
  skipped).
- `self_test` checks every sphere against its own stored dBA and EAA at 1e-9.
- `unit_test` is the Catch2 binary. Catch2 rejects `&&` and chained comparisons inside `REQUIRE`
  unless the expression is wrapped in a second pair of parentheses: `REQUIRE((a && b));`.
- `*_windowed_is_exact` compares the windowed peak search against marching every state. Those
  runs must agree.
- Tests of refusals and warnings use `PASS_REGULAR_EXPRESSION` on the message text.

**The goldens:**

- `golden/state_197421.csv` is D-derived trajectory state.
- `golden/197421_LAE.csv` comes from the C++ build. It was regenerated deliberately twice: once
  when the 4x EAA multiplier was removed, and once when spreading was capped at the sphere
  surface.
- `golden/metrics_*`, `stratified_*_golden` and `terrain_hills_*_golden` also come from the C++
  build.

Regenerate a golden only for a deliberate numerics change, say so in the commit, and use the
`regenerate-goldens` skill when the D oracle is involved.

**Several tests exist only to kill a specific mutant**, because the rest of the suite passed
with that mutant in place. Keep them when you refactor. The list, and the mutant behind each
test, is in [docs/internals/testing.md](docs/internals/testing.md). If a mutation check restores
a file with `mv`, `touch` the file afterward or make won't rebuild it.

**No CTest exercises plotting** (`-p`/`-f`/`-F`). CI runs a separate figure step after the
golden comparisons.

## Running the Tool

With no arguments, the tool runs the demo case:

- trajectory `tracks/S-76D/197421AC.csv`
- database `databases/S-76D_M3.nod`
- vehicle file `databases/S-76D.json`
- output `integrated/197421_LAE.csv`

Every flag is in `docs/MANUAL.md`, or in the CLI11 setup in `src/niceops.cpp`. Trajectory CSVs
are the 27-column NASA/FAA AC state format, or WGS84 tracks read by their header.

**Units.** `--units ft|m` (default ft) sets the length unit for the whole interface: every
length on the command line, every file that carries a length, and the figure axes. Computation
is metric internally, and the conversion is one seam in `run()`. Three rules matter:

- **Defaults are physical quantities.** They are written in meters and skip the conversion.
- **`--units` declares what the files already are and converts nothing.**
- **`write_state` stays metric**, because the state golden encodes internal values.

Speeds name their own unit (`12kt`, `6.2mps`), and a bare number is refused.

**Input checking.** A bad input must fail or warn with its cause named. It must never give a
plausible footprint with exit status 0. Errors and warnings go to stderr as well as to
`log/niceops.log`. A skipped trajectory exits 1. Readers refuse what they used to swallow:
non-increasing time, non-finite values, unit mismatches, text after data, and observer files
that yield nothing. Follow the same rules when you add a reader. The full list is in
[docs/internals/inputs.md](docs/internals/inputs.md).

## Architecture

The pipeline in `run()`/`process_trajectory` (`src/niceops.cpp`):

1. Load the `Database` and the `Trajectory`.
2. Build a `Database::SourceState` once per trajectory state (frame, blend).
3. For each observer and state, evaluate the level at the microphone, either broadband dBA
   with EAA or per band with ISO 9613-1 absorption.
4. Integrate in reception time over the outermost 10-dB-down window (`src/util/sel.hpp`).
5. Write the CSVs and optional figures.

A default run makes about 94 million evaluations, so anything left in the inner loop costs
about 100 million times over.

### Invariants that cut across modules

- **The default path is byte-identical to the goldens**, and new options must leave it that
  way. A degenerate case of a new model, such as calm air, uniform air or level terrain, should
  reproduce the existing model, at 1e-9 or bit-for-bit. A CTest should pin that.
- **The windowed peak search is exact.** It skips only states that a bound *proves* are below
  the threshold. Any new gain (ground, table, focusing, plate) must enter the bound, or
  `*_windowed_is_exact` fails. Options that need every state (history outputs, time weighting,
  day metrics, EPNL, `-a`, and others) march instead. See
  [performance.md](docs/internals/performance.md).
- **`--spectral` is the default whenever the database carries spectra.** The goldens use the
  broadband path, because `S-76D_M3.nod` has no spectra. `eval` and `eval_spectrum` share their
  cell and stencil, so any disagreement between them is physics, not interpolation.
- **The ground model returns a gain, not a loss.**
- **Compass bearings** are turned into the local frame by `--frame_bearing` through one
  rotation, `CompassToLocal`. The shipped AC tracks need 210.
- **panam is the reference for ported physics**: ISO 9613-1, `ega()`, EPNL, and A-weighting are
  each pinned to panam at 1e-9 dB or similar. A change there is a change to panam first. A
  wrong result may also be bad database data, not bad code.

### Module map

| Area | Sources | Notes |
| --- | --- | --- |
| Hot loop, SEL window, windowed peak search, load time, rejected optimizations | `niceops.cpp`, `sphere_database.cpp`, `util/sel.hpp`, `util/parallel.hpp` | [performance.md](docs/internals/performance.md) |
| Sphere database, Delaunay and exact predicates, `--blend`, legacy repair, broadband vs spectral | `sphere_database`, `predicates`, `delaunay`, `interpolator`, `condition_selection` | [database.md](docs/internals/database.md) |
| ISO 9613-1, the `RayPath` seam, `--atmosphere`, `--wind` | `atmosphere`, `propagation.hpp`, `atmosphere_profile` | [atmosphere.md](docs/internals/atmosphere.md) |
| `--refraction linear` / `stratified`, shadow zones, `ShadowTransition`, ducts, `--tunneling` | `refraction`, `stratified`, `airy` | [refraction.md](docs/internals/refraction.md) |
| Parabolic equation (split-step Padé, Crank–Nicolson), benchmark exact solutions, `--propagation_table` | `pe`, `pade_propagator`, `exact_solutions`, `propagation_table`, `tools/` | [pe.md](docs/internals/pe.md) |
| Receivers (plates, poles), ground effect, Faddeeva, weighting, `--terrain` | `ground_effect`, `faddeeva`, `weighting`, `terrain`, `plate_receiver.hpp` | [ground.md](docs/internals/ground.md) |
| EPNL, band frames, LAFmax/NAx/TA, DNL/Lden/CNEL, C-SEL | `pnl`, `band_frames`, `util/time_weighting.hpp` | [metrics.md](docs/internals/metrics.md) |
| AC and geodetic tracks, smoothing spline, hover nose, `--nose`, `--speed_reference` | `trajectory`, `geodetic`, `smoothing_spline`, `util/csv.hpp` | [trajectory.md](docs/internals/trajectory.md) |
| `--units`, input checking | `run()`, readers | [inputs.md](docs/internals/inputs.md) |
| Goldens, `self_test`, mutation-guarded tests | `golden/`, `tests/` | [testing.md](docs/internals/testing.md) |
| D replacements, port quirks kept on purpose | | [porting.md](docs/internals/porting.md) |
| Plotting | `plotting` (pybind11-embedded matplotlib, created lazily) | (this file) |

The roadmap for propagation (refraction, then urban) is in the Claude Doc "NICE-OPS:
Atmospheric and Urban Propagation Options". `reports/` holds a snapshot of it.

## Code Conventions

- **Port fidelity outranks correctness, with a bar for departing from it.** Leave the `NOTE:`
  comments that reproduce D quirks alone unless you are also regenerating the goldens, and say
  so if you do. Departing from D takes *a measured case that the D behavior is wrong*, not a
  reading that it looks wrong. The 4x EAA multiplier and the ambient mix met that bar.
- **Porting pitfalls already hit (don't reintroduce them):**
  - D's `iota` excludes its end (`iota_exclusive`: 120 grid values ending at 1475).
  - Don't sort `phis_scattered`/`thetas_scattered` in place.
  - D doubles default to NaN, which is why there is no ambient unless `-a` is given.
  - `approx_equal` mirrors `approxEqual` (relative to the *second* argument, 1e-2/1e-5).
    Tightening it drops channels.

  When a D construct looks obvious, check it against Phobos with a throwaway `ldc2` program.
  See [porting.md](docs/internals/porting.md).
- `using namespace niceops;` at the top of `src/*.cpp` is intentional (`// NOLINT`).
- Logging goes through `src/util/logger.hpp` (`log/niceops.log`). `fatal()` is `[[noreturn]]`
  and aborts, so RAII cleanup doesn't run on fatal paths.
- **Licensing:** the project is MIT. `THIRD_PARTY_NOTICES.md` holds verbatim notices for the
  fetched dependencies; regenerate it from `build/_deps` when one changes. The `Agg` backend in
  `plotting.cpp` is load-bearing: the Qt backends would pull in GPL PyQt.
