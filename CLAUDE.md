# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NICE-OPS (Noise Informed Community Environment Operation Planning System) is a helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer grid (or scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

The implementation is C++17, under `src/`. It began as a port of a D original, which is why several comments reference D semantics — those are load-bearing: the port reproduces the D reference bit-for-bit on the demo case, and the goldens come from it.

- `src/` — the implementation. **Every `.cpp` lives here and is listed in CMakeLists.txt**; there are no C++ sources at the repo root.
- `README.md` — a short description, building and a first run. **`docs/MANUAL.md`** is the user manual: every option, the theory behind each model, and how each was checked. Put user-facing detail in the manual, not the README.

The D sources (`niceops.d`, `sphere_database.d`, `trajectory.d`, `splines.d`, `netcdf.d`), along with `cmake-d/` and the vendored `external/matplotlibd/`, were removed once the port was verified complete: every D function, all 19 CLI options and all 10 `unittest` blocks have C++ counterparts. **They remain in history at `c74bab3`** (on `origin/no-stk`), which is what to check out if you ever need the oracle again — see the golden-file section.

When you add a C++ source, put it in `src/` and add it to the right target — a file that compiles but is never linked in shows up as an undefined-symbol error far from its cause.

## Build

Requirements: a C++17 compiler and netCDF (shared library). Plotting also needs Python 3 **with development headers**; matplotlib is installed by the build into `build/pyenv` rather than being a prerequisite. CGAL is gone (and Boost, GMP and MPFR with it, which were only ever CGAL's transitive dependencies), as is Eigen -- `src/vec3.hpp` covers what was actually used of it. On Ubuntu, as CI does:

```sh
sudo apt-get install libnetcdf-dev python3-dev python3-venv
```

`python3-venv` is easy to miss: Debian and Ubuntu package `python3 -m venv` and its
ensurepip separately from `python3` itself, and without it the build's environment
step only *warns* before falling back to the system interpreter's packages. That
fallback is silent enough that CI checks the environment explicitly after building.

With `-DNICEOPS_PLOTTING=OFF` the only prerequisite left is netCDF.

`-DNICEOPS_PLOTTING=OFF` drops plotting and the whole Python dependency with it; CI builds this configuration too, since nothing else compiles `src/plotting_disabled.cpp`. `-DNICEOPS_BUNDLE_PYTHON=OFF` keeps plotting but uses the interpreter's own packages instead of the build's environment.

CLI11, nlohmann_json, Catch2, and pybind11 are pulled at configure time via `FetchContent`, so the first configure needs network access.

Link-time optimization is on by default (`-DNICEOPS_LTO=OFF` turns it off), where
`check_ipo_supported` says the toolchain has it. It is worth 18% of the acoustic loop,
because that loop and the database live in different translation units; it leaves the
golden comparison passing unchanged, and GCC 16 and Apple Clang produce identical output
with it.

Build out-of-source, and **always with explicit `-S`/`-B`**:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON
cmake --build build -j8
```

> **Use `-S`/`-B`, never `cmake <srcdir>` from a build directory.** The repo root still holds a stale `CMakeCache.txt` from the D era. When the directory you hand `cmake` contains a cache, CMake treats that argument as the *build* tree — so `cd build && cmake ..` silently configures **into the source root**, reusing the 2023 cache. That cache carries paths from an old Homebrew prefix (`/usr/local/...`) and produces baffling, unrelated errors such as `Imported target "CGAL" includes non-existent path "/usr/local/include"`. Deleting the root `CMakeCache.txt`/`CMakeFiles/` disarms the trap permanently.

Notes:
- The repo root contains stale in-source build artifacts from the D era (CMakeCache.txt, Makefile, `libniceops*.a`, binaries) that are gitignored — do not build in-source.
- **Asserts: live when testing, compiled out in shipped binaries.** `-DNICEOPS_ASSERTS=ON` strips CMake's `-DNDEBUG` from the optimized configurations and defines `NICEOPS_ASSERTS_REQUIRED`; `tests/test_build.cpp` fails if such a build still has NDEBUG. CI and the test recipe below pass the option. A plain `cmake -DCMAKE_BUILD_TYPE=Release` build is a shipped build, so its asserts are off (about 4% of the demo's CPU). The D original ran with asserts live, and `sphere_database.cpp`'s `assert(isFinite(...))` guards are what the suite exercises, so test with the option on. Every assert is side-effect free, so a failing assert is the only place the two builds can differ.
- If Xcode proper isn't installed (Command Line Tools only), pass `-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`. No special flags are needed; the vendored SPLINTER library that used to require them is gone.
- CGAL 5.0+ is header-only and exports `CGAL::CGAL`. Do **not** reintroduce `include(${CGAL_USE_FILE})` — `UseCGAL.cmake` was removed in CGAL 6, leaving `CGAL_USE_FILE` empty so the `include()` errors out.
- **CI builds on GCC/libstdc++; a macOS clang build will not catch what that finds.**
  Two classes of bug are invisible here and fatal there: headers libc++ drags in
  transitively but libstdc++ does not (`<algorithm>` is the usual one), and static-library
  link cycles, which ld64 resolves and GNU ld does not. Before pushing anything structural,
  build the tree with Homebrew GCC rather than waiting for a CI round-trip:

  ```sh
  cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON \
        -DCMAKE_CXX_COMPILER=$(brew --prefix)/bin/g++-16 -DNICEOPS_PLOTTING=OFF
  cmake --build build-gcc -j8 && ctest --test-dir build-gcc --output-on-failure
  ```

  The goldens pass under GCC too, so a compiler change does not move the 1e-9 comparison
  despite the different FMA contraction. **Do not assert bit-identity against a libm call
  on a literal argument** in a test: the literal gets constant-folded at build time (GCC
  folds via MPFR, correctly rounded) while the function under test reaches libm at run
  time, and the two differ by an ULP. `acosc` was written that way and only GCC noticed.
- **`arrays/coyle_approach.csv` is the only scattered-observer set in the repo**, and it is a
  working `-b` input again: every one of its rows was being dropped, because its fields are
  padded *before* the comma and `read_numeric_csv` counted trailing whitespace as junk. The run
  then evaluated zero observers and exited 0. An observer file that yields nothing is now an
  error. All 84 shipped AC tracks pad on the *left*, which `std::stod` always skipped, so the
  goldens never saw this.
- **`SmoothingSpline` at lambda = 0 is the natural cubic *interpolating* spline**, not a
  polyline. `PenalizedSystem::fit` used to short-circuit there and leave the curvature zero;
  the lambda = 0 limit of its own system is `R gamma = Q^T y`, which has to be solved. See the
  `--path_smoothing` note below for what that cost.
- **`state/`, `integrated/`, `log/` and `figures/` are not tracked**, so a fresh checkout
  has none of them and every writer must create its own parent directory
  (`ensure_parent_directory` in `src/util/paths.hpp`). `Trajectory::write_state` once
  opened its ofstream directly and so aborted on a clean CI checkout while passing on any
  working tree that had run the tool before. A working tree is therefore *not* a faithful
  test of a checkout — verify against a clone:

  ```sh
  git clone --no-hardlinks . /tmp/fresh && cd /tmp/fresh
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON && cmake --build build -j8
  ctest --test-dir build --output-on-failure
  ```

  Note what such a failure looks like: `fatal()` logs its reason and calls `abort()`, and
  CTest reports `Subprocess aborted`. The log line goes to `log/niceops.log` and, since
  warnings and errors are echoed to stderr as well (`Logger::write_line`), to the terminal --
  it used to reach the log alone, which is why CI still dumps the log on failure.
- **Keep the C++ sources in one static library.** `interpolator.cpp` used to be its own
  `interpolator` target, a D-era seam that walled off CGAL's heavy headers. Once CGAL went,
  that split was a *cycle*: `interpolator.cpp` calls into `Delaunay`, which lives in
  `niceops_core`, while `niceops_core` linked `interpolator`. GNU ld walks archives in a
  single pass and fails (`undefined reference to niceops::Delaunay::barycentric`); ld64
  resolves it either way, **so this breaks Linux only and a macOS build will not show it**.
  If you add a target, check that its dependencies point one direction.

### Where the database load time goes

A default run is now ~0.93 s, and ~0.61 s of that is `Database::Database`. Most of it is not
ours: `S-76D_M3.nod` has **one netCDF group per operating condition**, 1312 of them, and
netCDF-4 reads every group's metadata eagerly at `nc_open` -- it opens all ~8000 HDF5 datasets
in the file before returning. Measured split:

| | share of the load |
| --- | --- |
| `nc_open` (eager per-group metadata) | 41% |
| our variable reads (metadata lookup, not I/O) | 25% |
| `nc_close` (tearing that metadata down) | 12% |
| Delaunay build over the conditions | 14% |
| grid fill and the rest | 8% |

(The table is the load before the three changes below, when it was ~0.75 s; the first three rows
are unchanged and the last two shrank.) The first three -- 78% of that -- are netCDF's own cost
for this file layout and are not reachable through its API; the five tiny scalars per group cost more than the four big arrays, because
what they pay for is the lookup, not the bytes. What *was* ours, and is done: the Delaunay build runs on its own thread while `nc_close` tears
the file down (the build touches no netCDF state, so it is safe even though HDF5 here reports
thread-safety off; ~85 ms, `NcGuard::close`); each channel's grid cells are derived once and reused
while the angle labels are unchanged, which is every sphere but the legacy hover ones (the fill
and the spectral layout, 47 -> 10 ms); and the optional `ShaftPower` probe no longer logs a
line per group (1312 of a run's 1349 log lines). Output is byte-identical (golden compare).
The real fix is in the database builder
(panam): the same data in one group with `(condition, channel)` arrays would load in a
fraction of the time. Do not go looking for it in `sphere_database.cpp`.

### Golden reference files

`golden/state_197421.csv` is captured from the **D reference implementation**; the C++ port
matched it to **1e-9** (actual max relative difference ~1e-14, a few ULPs of LDC-vs-clang
round-off). It is trajectory state, which no propagation change touches, so it still stands as
the D-derived oracle.

`golden/197421_LAE.csv` **no longer comes from D.** Removing the spurious 4x EAA multiplier
(see below) deliberately departs from the D original, so the D binary can no longer produce
this file and regenerating it from D would only restore the bug. It was regenerated from the
C++ build at the commit that removed the multiplier. The port's agreement with D was verified
before that change and that verification stands on its own; D is retired as an oracle for the
propagation path.

`golden/197421_LAE.csv` was regenerated a second time when spherical spreading was capped at the
sphere's surface (`Database::eval` and `eval_spectrum` take `max(range, R)`, and absorption starts at
the surface): 27 of the 14400 cells change, by at most 2.0 dB in SEL, and they are exactly the
observers within about 45 m of the track, whose closest approach to a source state is inside the 30.48 m
sphere. Before that, the spreading law was extrapolated inward (1/r) and put those observers up to 4.4
dB above the sphere-surface level in the instantaneous level, and diverged at r = 0. It, too, was made
from the C++ build, the one that carries the change; `golden/state_197421.csv` did not move.

`regression_golden` and `state_golden` both pass.

To regenerate them, follow the `regenerate-goldens` skill
(`.claude/skills/regenerate-goldens/SKILL.md`) — it builds the D original in a scratch
worktree and covers the three traps (wrong `niceops.d` revision, D's 6-significant-figure
default output, hardcoded `/usr/local` paths) that otherwise make the result unusable.

## Tests

From the build directory:

```sh
ctest --output-on-failure -C Release   # what CI runs (no -E; the whole suite passes)
ctest -R unit_test -V              # one CTest test
```

The unit tests are a Catch2 binary, so you can run a single case or tag directly — far faster than going through CTest:

```sh
./niceops_unit_tests "[sphere_database]"      # by tag
./niceops_unit_tests "lookup_angles"          # by test name
./niceops_unit_tests --list-tests
```

CTest defines about ninety tests (`ctest -N`; most are `niceops` runs compared against another run), all with the **source root** as their working directory (the demo case's relative paths resolve there). The main ones:

- `regression_default` — runs `niceops -S state/ci_state.csv`; one run feeds both golden comparisons.
- `self_test` — runs `niceops --self_test`, validating every sphere in the database against its own stored levels (1,844,672 samples). **Passes.** It was `broken_self_test`, excluded from CI, while legacy hover directivity made it fail on 7,788 samples; that is fixed (see below), so it is now a real test and CI runs the full suite with no `-E` exclusion. It checks the stored dBA **and EAA** at 1e-9 (eval() returns the stored double exactly at a sphere's own condition and a grid node -- 0 difference over 22 million samples; the D original's 1% let a channel filed one grid cell off pass 60-90% of samples), and **`--self_test` now exits when it is done** instead of going on to run the demo case, which rewrote the file `regression_golden` reads under any test running beside it.

**Mutation-guarded behavior.** A review mutated the code one line at a time and wrote tests for what survived; these are the ones to keep when refactoring (each fails on the mutant it was written for): `an exported sphere reloads as the sphere it was, cell for cell` (swapped dBA/EAA, dropped load-factor correction), the AAM readback, `a negative excess-attenuation factor turns the bounds off` (the windowed run's exactness), `unevenly spaced knots give a continuous, natural spline` (every other spline test used `linspace`, where a transposed h[j]/h[j+1] is invisible), `the fit projects the wind onto the path's own direction` and `a short fit span is widened upward` and `a flush microphone sees the direct ray arrive at its own angle` (a flipped northward wind sign, widening downward, sin/cos swapped in the plate branch), `the plate's azimuth is the ray's bearing turned by the microphone's offset` and `outboard puts the microphone away from the track` (`plate_frame_azimuth` and `outboard_azimuth` live in `plate_receiver.hpp`, one definition for both callers), `add_bounce_loss is the plane-wave reflection loss`, the three `stratified_*_golden` CTests (windowed-versus-marched compares two runs of one binary, so dropping the extra arrivals from the sum moved the duct footprint 7.8 dB with every stratified test green; regenerate those goldens with the build that changed the arrivals on purpose), `scattering_windowed_is_exact` (no CTest ran `--scattering`), `propagation_table_bump_*` (an eight-bearing table with +25 dB beyond 1500 m and a cos(bearing) tilt, written by the hidden Catch case `[.bump_table]` from the ground-model table, run with `--frame_bearing 210`: a bound that left the table's gain out differed from the marched run by 0.57 dB in 21 of 23 cells and passed every CTest; the table-with-absorption twin must differ from it) and `CompassToLocal turns east/north ...` (the table's bearing is its `east_north` inverse, swapped atan2 arguments passed everything with only one-bearing tables), `every simplex has an empty circumsphere` (through the public API: whatever simplex a query lands in, no input point is inside its circumsphere) and `in_sphere agrees with an exact integer oracle` (128-bit Bareiss on lattice points, dimensions 1-5, ties included; the oracle needs no sign convention because the simplex's centroid is always inside). **A mutation check that restores a file with `mv` leaves it older than the objects built from the mutant, and make then considers it up to date**: `touch` it, or `cp` it back.

- `unit_test` — Catch2 cases in `tests/`, ported from the D `unittest` blocks. All cases pass (tens of thousands of assertions). Note Catch2 rejects `&&` and chained comparisons inside `REQUIRE` unless the whole expression is wrapped in a second pair of parens: `REQUIRE((a && b));`.
- `regression_golden`, `state_golden` — `golden_compare` diffs the run's CSV output against a golden file elementwise at 1e-9 relative tolerance, skipping non-numeric header lines.

**No CTest case passes `-p`/`-f`/`-F`, so none of them executes the plotting path.**
That is how `plotting.cpp` went years missing a `<pybind11/stl.h>` include: it
compiled clean and threw `cast_error` on every call. CI covers the gap with a
separate step that runs the demo for real figure output, after the golden
comparisons -- it rewrites `integrated/197421_LAE.csv`, which `regression_golden`
reads. Note the figure names come out as `state_CT_.png`: `niceops.cpp` builds them
as `base + "_CT_" + ext`, so the separator lands before the extension.

### Legacy databases (`database_version` < 1)

Databases carry a root-level `database_version`. Absent or 0 means the file predates the versioned builder, and its **hover spheres have broken directivity**: the upper surface was completed with `-phi` instead of `180 - phi`. Since the source spheres span only phi in [-90, 90] (phi = 0 straight down), negating is a left/right flip *within* the lower half that maps the range onto itself — so the mirrored channels wear lower-hemisphere labels and the upper surface is missing. Forward-flight spheres in those files were completed correctly; only hover is affected (in `S-76D_M3.nod`, exactly the 24 spheres at mu = 0).

Only the labels are wrong, so `Database::Database` repairs them on load via `repair_legacy_hover_directivity()`: for the mirrored half, `phi_correct = wrap(180 + phi_stored)`. It refuses to touch anything that does not match the broken shape (must be confined to [-90, 90], even channel count, halves actually mirrored) and logs how many spheres it fixed. `self_test` applies the same repair, or it would compare a grid built from repaired labels against queries at the raw ones.

Two consequences worth knowing:

- Channel angles are read **per sphere**, not once from sphere 0. Reusing sphere 0's array happened to give the right answer for `S-76D_M3.nod` — the repaired hover labels are byte-identical to sphere 0's — but only because every forward sphere there shares one layout. It would silently misfile data in any database with heterogeneous layouts.
- The repair leaves the demo footprint bit-identical, for the same reason, so `golden/` did not need regenerating.

Regenerating the database with the current builder (panam, `DATABASE_FORMAT_VERSION = 1`) makes the fix-up unnecessary; the loader then trusts the file as written.

## Running the Tool

Default demo case (no arguments): trajectory `tracks/S-76D/197421AC.csv`, database `databases/S-76D_M3.nod`, vehicle config `databases/S-76D.json`, output `integrated/197421_LAE.csv`. For the full option set see `docs/MANUAL.md`, which documents every flag, or the CLI11 setup in `src/niceops.cpp`.

**`--units` picks the length unit for the whole interface** -- `ft` (default) or `m` -- covering
every length on the command line (`--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
`--receiver_height`, `--path_smoothing`, `--roughness_length`, `--surface_layer_depth`, the
anemometer height in `--wind`), every file carrying a length, and the figure axes. Speeds are not
lengths: `--wind`'s speed names its unit (`12kt`, `6.2mps`, `14mph`), as an `--atmosphere`
file's wind columns must, and a bare number is refused rather than read in `--units`.
Internal computation is metric (`ft2m` and `LengthUnit` in `src/util/constants.hpp`); the
conversion is one seam in `run()`.  Four of those options used to be metric while every file was
feet, which is the inconsistency this removes.

Three things are load-bearing here:

- **A default is a physical quantity, not a number in the current unit.**  Defaults are written
  in meters and skip the conversion, which is why the seam tests CLI11's `count()` rather than
  scaling everything.  Scale the defaults too and `--units m` shrinks the observer grid to a
  third of its extent.
- **`--units` declares what the files already are; it converts nothing.**  The shipped tracks are
  feet, so `--units m` on the demo case reads feet as meters.  Nothing can detect that.  Where a
  column header names its own unit (`alt_ft`, `alt_m`) it is checked against `--units` and a
  disagreement fails the file; a bare `alt` is always accepted.
- **`write_state` stays metric** whatever `--units` says: `golden/state_197421.csv` is the
  D-derived oracle and encodes the internal values.

`fps2g` is chosen as a whole constant rather than composed as `unit.to_meters * mps2g`, because
the product of two constants is not the same double as two multiplies and the state golden is
pinned at 1e-9.  Verified end to end: a metric twin of the demo track under `--units m`
reproduces the golden footprint to 1.5e-12 dB. Trajectory CSVs are the 27-column NASA/FAA AC state format. Vehicle JSON keys: `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg, converted to Newtons), `Atmosphere.density`.

### Input checking and reporting

A bad input used to give a plausible-looking wrong footprint with exit status 0 and nothing on
screen; these now fail, or warn, with the cause named. Know them before adding a reader:

- **Errors and warnings reach stderr** as well as `log/niceops.log` (`Logger::write_line`); `-L`
  still decides what is written at all. A skipped trajectory is `Skipping <file>: <why>` on the
  terminal and exit status 1.
- **An exception in an observer-loop worker is carried to the caller** by `parallel_chunks` (it
  terminated the process at any thread count above one), so `input_error` skips one trajectory
  as `run()`'s handler intends.
- **Readers refuse what they used to swallow.** An AC track's time must increase and every
  column the states are built from must be finite; a geodetic track's time, altitude and
  velocity must be finite; a unit declared in a geodetic header (`Alt (ft)`, `Climb Rate (fpm)`,
  `Time (ms)`) is checked against `--units` (`check_declared_units`); a text line *after* the
  first row of data is an error (`CsvReport`; text above it is a header); observer files keep
  one width; `--atmosphere` rejects non-finite cells. `--trim` past the crop is diagnosed, and
  `-D` on a geodetic track with no clock needs `--start_utc` (it defaulted to UTC midnight, i.e.
  the night penalty, silently).
- **Atmosphere values are range-checked** (`temperature_problem` and friends in
  `atmosphere.hpp`, the bounds `--atmosphere` already used) from the command line and the vehicle
  file. The spectral path's log says where each value came from, and refuses a spectral run whose
  humidity is not given from any source. It is not assumed to be 20 percent:
  the CTest lines that used the old default pin `--humidity 20`. The vehicle file's tip speed,
  weight, rotor radius and C_T are compared with the database's, and a mismatch over 2 percent is
  refused unless `--accept_vehicle_mismatch` is given.
- **A database must vary where the triangulation needs it** (advance ratio, flight path angle,
  and CT unless `fixed_load_factor`) and may not hold two spheres at one condition; both were
  silent NaN or hang before. Test fixtures with one thrust coefficient carry the flag.
- **Several trajectories in one run** name `-i` / `--band_history` files for their track and
  derive `-S` names as `-o` does, and two trajectories asked to write one file are refused.
- **CTests that test these pass on the message text** (`PASS_REGULAR_EXPRESSION`), which
  replaces the exit-status check, so a usage error or a missing file cannot satisfy them;
  `golden_compare` itself is tested on NaN and infinity.

## Architecture

Pipeline in `main()` (`src/niceops.cpp`): load `Database` + `Trajectory` → for each observer and each trajectory state, compute instantaneous dBA via `Database::level_at_observer` → integrate the time history into a sound exposure level (SEL) using the 10-dB-down-from-peak window (`peak_range`/`sound_exposure_level` in `src/util/sel.hpp`) → write CSV output and optional matplotlib figures.

### Why the hot path is shaped the way it is

A default run evaluates 14,400 observers × 6,551 states = 94.3M times, so anything left in `eval` costs ~100 M×. Do not "simplify" these back:

- **`Database::SourceState` is built once per trajectory state, in `main`, and reused across every observer.** It carries the source position, the `SphereFrame` (unit TPP normal, nose, starboard) and the `Blend` (four bracketing spheres + barycentric weights). The `Blend` is a Delaunay `locate()` walk and was ~60% of the whole run when it sat inside `eval`. Anything that depends only on (μ, γ, CT) or on attitude belongs here, not in `eval`.
- **`lookup_angles` is split**: `sphere_frame(velocity, normal)` is per state, `lookup_angles(frame, direction)` is per observer. The four-argument overload still exists for the unit tests and composes the two.
- **`lookup_angles` takes three dot products and no projection.** `{nose, normal, starboard}` is orthonormal, so the propagation direction's components in the plane normal to the nose are just `normal.dot(d)` and `starboard.dot(d)`; forming the projected vector and normalizing it (a square root and three divides) scaled both by the same positive length, which the angle does not depend on. `phi` is then `atan(across / |up|)` below the plane, and `copysign(pi, across) - atan(across / up)` above it, which is `atan2(across, -up)` without the quadrant handling. Above the plane is read from the database's own upper half, when every condition group stores one (`upper_hemisphere_measured()`); a database without it folds the direction onto the lower half, as before. The 0/0 case at the poles is caught explicitly. The seam at 180 degrees is a real cell: `GridAxis::set_period(360)` makes the cell from the last phi knot to the first, one period on, and `locate()` interpolates across it.
- **The level grid is one flat `vector<Sample>`** (`{dBA, EAA}` pairs), indexed `(i_phi * n_theta + i_theta) * n_conditions + condition`, not `vector<vector<vector<double>>>`. Both metrics are read together, so they share a cache line. A condition-major layout was tried and is *worse* — it also breaks bit-exactness via different FMA contraction.
- **`GridAxis` replaces a binary search** with arithmetic when the axis is uniform, then corrects the estimate against the neighboring knots so the answer provably equals `upper_bound`'s. The searches alone were about a third of the loop. `tests/test_sphere_database.cpp` sweeps both paths against each other; keep that test if you touch it.
- **`DirectionBounds::bin` lives in the header.** The observer loop calls it ~94 million times from `niceops.cpp` while it is defined in `sphere_database.cpp`, so out of line it inlined only when `NICEOPS_LTO` was on and the toolchain supported it. It is fifteen instructions and should not depend on a build option; its tangent table is an `inline const` variable rather than a function-local static, which would test a guard on every call. The bin geometry — `bin`, the tangents, `narrowest_theta_bin_degrees`, `azimuth_bin_degrees` — sits on `DirectionBounds` together, because the table build has to pad by the same binning the hot path queries with, and three independent spellings of it is how that silently stops covering.
- **`GridAxis::locate` gives the cell and the interpolation weight together**, and is what `eval` and `eval_spectrum` call. Going through `bracket()` meant reading four knots back out of `values_`, four subtractions and a divide by the cell area, all to recover a number the uniform spacing already implies; `locate` returns the bracketing pair plus the fraction of the way across the cell, so the four corner weights are products of `(1-f)` and `f` and already sum to one. Worth 18% of the loop. It skips `bracket`'s correction loops, which is why it needs the tighter `exactly_uniform_` test in `reset()` rather than the tenth-of-a-step one: round-off can only put the estimate on the wrong side of a *knot*, where the two cells interpolate to the same value anyway. The wrap-around cell at the ends is not uniform and falls back to the knots. `tests/test_sphere_database.cpp` sweeps `locate` against `periodic_bilinear_interp` over uniform, irregular and wrapping cells — the two derivations have to agree, and that test is what says so.
- **Both `eval` and `eval_spectrum` build their stencil with `cell_stencil`.** It turns the two `GridAxis::Cell`s into the four corner indices and the four bilinear weights. Sharing `locate` alone got them the same *cell* but left each deriving its own weights, which is what keeps the claim below — that a broadband/spectral disagreement is the physics and never the interpolation — true rather than merely intended. `periodic_bilinear_interp` is no longer on either path; it is the reference `locate` is swept against in `tests/test_sphere_database.cpp`, and the cell/apply split it used to have (`BilinearCell`, `bilinear_apply`) existed only for the old `eval` and is gone.
- **The interpolation carries `{dBA, EAA}` as a two-lane vector** (`Pair` in `sphere_database.cpp`). Every step applies the same weight to both metrics, and a `Sample` is exactly a two-double vector in memory, so blending the four conditions at four corners and then the corners themselves is half as many instructions. Clang would not do it on its own from the scalar form or from a two-element array; the GCC/Clang `vector_size` spelling is explicit, with a plain-struct fallback that compiles to what the code said before.
- **Divides in the loop are multiplications by a stored reciprocal**: `inv_sphere_radius_`, `1e-3` for the EAA scaling, and one reciprocal for normalizing the propagation direction instead of three divides. A double divide is about four times a multiply. These move the result by an ULP, which is the only thing in this round that does.
- **The observer loop's chunks are dealt out on demand, not divided up front** (`src/util/parallel.hpp`). An equal split gives every worker the same slice, so the wall time is the slowest core's -- and the cores are not interchangeable. On an M4 Pro (10 performance + 4 efficiency) the old equal split ran no faster on all 14 threads than on 10, the efficiency cores' slices being the critical path; pulling chunks lets the fast cores take more of them and gets the 14-thread run back under the 10-thread one. Chunk size is `count / (8 * workers)`, which is fine-grained enough to balance and coarse enough that the per-chunk scratch allocation stays noise. The output still does not depend on the split -- every observer writes only its own slot -- so runs at any thread count remain byte-for-byte identical.
- **`trace()` is force-inlined and `Logger::write` is `noinline`** (`src/util/logger.hpp`). `eval` and `level_at_observer` make ~7 trace calls each; as ordinary calls that was a sixth of the loop spent entering functions that immediately returned. The level lives in `detail::min_level`, outside the `Logger` singleton, so testing it does not touch a function-local static's guard.
- **The database grid fill walks channels and places them**, rather than scanning all 1,406 channels for each of the 1,332 cells (2.5 billion `approx_equal` probes, most of the load time). Two semantics are load-bearing and preserved: a later channel still wins a shared cell, and `dba_max` covers matched channels only.

- **Link-time optimization is on** (`NICEOPS_LTO`, see CMakeLists.txt). The loop lives in `niceops.cpp` and `Database::level_at_observer` in `sphere_database.cpp`, so without it every one of ~100 million evaluations was an out-of-line call across a translation unit. Worth 18% of the loop on top of the source changes — and only 5% before them, because what it buys is hoisting and inlining that the divides and knot lookups were hiding.

Current cost on the demo case: database load ~0.75 s, acoustic loop 1.55 s of CPU (0.16 s wall on 14 threads). It was 2.2 s and 11.7 s originally, and 3.2 s of CPU before this work.

Where the 19.4 ns of an *evaluation* goes, from an ablation (`full`, `levels only`, `angles only`, `geometry only`, `eval at fixed angles`): 33% the two inverse trig calls, 35% `eval` itself (gather, blend, `log10`, arithmetic), 8% the vector geometry, 7% the SEL integral, and ~17% the difference between a real `eval` and one at synthetic angles -- which is *not* cache misses, see the gather item below.

Several further optimizations were measured and **rejected**; the numbers are here so they are not re-derived:

- **Hand-rolled `acos`, `atan` and `log10`.** They are a third of the loop, which makes them look like the obvious target, but this platform's libm does `acos` in 2.3 ns, `atan` in 2.0 and `log10` in 1.9 — around 8 cycles each, table-driven. A polynomial accurate enough to hold the golden (~1e-12 rad, so the level moves by less than 1e-10 dB) needs a square root and a degree-10 minimax and is *slower*; one cheap enough to win (degree 4, ~1e-5 rad, ~1e-3 dB on the footprint) saves only about 1.5 ns of the 3.55 an in-situ call costs. Bad trade in both directions.
- **`-fno-math-errno`.** The hope was that marking the libm calls pure would let the compiler overlap them. No measurable change: it will not software-pipeline across an opaque call, and the three calls in an evaluation are already independent of each other and so already overlap.
- **Storing the grid as `float`.** Would halve the 28 MB gather. It costs ~1e-5 dB, which is
  within engineering accuracy but ends the golden comparison as a bit-level regression test.
  Not taken without being asked -- and see the next item, which says it would not help anyway.
- **Anything aimed at the gather's locality.** The four blended conditions of a cell span 845 of
  the 1312 indices on average and touch 3.1 cache lines per corner, never 1, so an evaluation
  pulls ~800 bytes of lines to use 256 -- which looks like an obvious 3x to reclaim. It is not:
  the whole trajectory uses only **52 distinct blends**, so the four needed condition slices
  could be compacted into 4.4 MB of contiguous 64-byte rows (nothing like the rejected slab's
  +131 MB). Before building it, the ceiling was measured the cheap way -- run the same loop
  against a one-condition database, whose 21 KB grid makes the gather perfectly local -- and it
  came out at 1.82 s against the full database's 1.80 s. **The loop is compute-bound, not
  miss-bound**, which is what the slab note already said; this is an independent confirmation
  from the other direction. Measure that ceiling before touching the layout again.
- **Per-cell polynomial inverses for `acos`/`atan`.** The angles are never wanted for their own
  sake -- only a grid cell and a fraction across it -- so a polynomial per 5-degree cell mapping
  cos(theta) straight to the fraction would skip the inverse trig entirely. It cannot work:
  acos has a square-root singularity at the poles, so no polynomial in cos(theta) represents it.
  Degree 7 still leaves 0.31 degrees of error in the two polar cells, and 5e-5 degrees in the
  best interior ones; changing variables to sqrt(1 - a) to fix it puts a square root back.
- **Culling states the SEL window will discard, by bounding them one at a time.** 44% of
  (observer, state) pairs fall outside the window, but a bound that proves a *particular* pair
  is outside it has to be tighter than 10 dB, and `level <= sphere max - spreading` is 11.4 dB
  loose at the median -- it proves only 7.6%. That is the wrong question; asking it the other
  way round is what works, and is now the default path. See "The SEL window, the oversampling,
  and the windowed peak search" below.

### The SEL window, the oversampling, and the windowed peak search

**The window is the outermost 10 dB-down interval.** `ten_db_down_interval()` in
`src/util/sel.hpp` runs from the first to the last sample at or above peak - 10 dB, dips
below the threshold between them included. That is the duration of the aircraft noise metrics
(the EPNL duration of 14 CFR 36 A36.4.5, and what panam's `sound_exposure_level` uses, so a
prediction and a measurement are integrated alike). It replaced the D original's walk out from
the peak, which stopped at the *first* sample below the threshold on either side and so cut an
event off at its first dip: an approach heard from upstream lost its second rise, and on the
2017 approach data the two rules differ by up to 0.7 dB at a single microphone. The change was
made deliberately (commit `9d6e986`, which regenerated `golden/197421_LAE.csv`), so the old
first-crossing walk, `expand_window()` and the "left alone under the port-fidelity rule" note
are gone from the code; do not reintroduce them from an older copy of this file. The integral
is taken over *reception* time (`reception_factor`, 1 + dr/dt / c), with the step from the
trajectory itself (see the time-step note below).

**The SEL integral is oversampled, and the discrete window is why that cannot be used.** The
sound exposure integral is a Riemann sum over a 50 Hz trajectory, oversampled by roughly
thirty times. Decimating the demo case's 6551 states moves the footprint's **mean** by 0.008 dB
at stride 2 and 0.088 dB even at stride 64. What blocks a 2-8x cut in the whole acoustic loop
is the window: decimation steps over a transient dip, the window's end jumps, and the level
jumps with it. The error does not shrink with stride the way a quadrature error would --

| stride | max cell error, first-crossing (the retired rule) | outermost-crossing (current) |
| --- | --- | --- |
| 2 | 1.93 dB (104 cells > 0.1 dB) | 1.13 dB (86) |
| 8 | 2.34 dB (619) | 1.13 dB (271) |
| 64 | 2.48 dB (3150) | 1.54 dB (1346) |

-- which is the signature of a discrete window flip, not of the integral. So decimation is
**not** safe as a default under either definition: the mean is negligible but individual cells
move by 1-2 dB, and the cells are the deliverable. It stays available deliberately, through
`adapt()` and `-s` (both currently commented out in `main`).

**Skipping states the window cannot use: the windowed peak search.**
`windowed_sound_exposure_level` / `windowed_exposure` in `src/niceops.cpp` do not evaluate every
emission time. A coarse sweep (stride nt/128) seeds a best level; one pass then evaluates every
state `DirectionBounds` cannot prove is below that best level, handing the same geometry on to
the evaluation so a state that survives the test pays nothing for it; and
`ten_db_down_interval` finds the window's ends by scanning in from both sides, skipping a state
the bound puts below the threshold (`certainly_below`) and evaluating everything between the
ends. It is **exact** -- the output is byte-identical to evaluating every state, on every track,
grid, observer set, receiver and thread count tried -- because the integral and the window are
the shared `ten_db_down_interval()` and `sel_over_window()`, templated on a level accessor (a
filled array for `peak_range`, a memoizing lambda here), not a second copy.

It runs on the **spectral** path only: `--spectral` (the default whenever the database carries
spectra), `--ground_effect`, plates, the propagation table (whose largest value per band enters
the bound), and `--refraction` (the
bound is taken along each ray's launch direction at its ray tube's range; for a sum of
stratified arrivals, the loudest arrival's direction over the harmonic sum of their squared
ranges, `bound_geometry`). It is skipped, and the loop marches every state, for `-a` (the level
is transformed after eval), `-i` and `--band_history` (the whole history is written),
`-D` (the day-night window is found on a penalized history, which a window found on the raw one
need not cover), `--db_cutoff` (which culls on its own terms), and any database with a negative
excess-attenuation factor, where dropping that term is not a loss.

What it saves depends on how dear an evaluation is. Under the outermost-crossing rule it
evaluates most of what it is shown: on the broadband demo case it evaluated 93.9% of the pairs
(the window itself is 59%, and the bound proves only about 6% of the pairs outside it) and ran
about 4.6% *slower* than marching every state (measured twice, output identical), so **broadband
runs no longer build the table or use the search**; `Database::direction_bounds`, the broadband
builder, remains, tested, for a broadband search that pays. A spectral evaluation costs 150 ns
and a `--ground_effect` one a microsecond, and there the search pays. The main pass keeps each
skipped state's numerator and squared range (`skipped_numerator`, `skipped_squared_range`) so
that `certainly_below` does not recompute the geometry -- a full ray solve under refraction --
for the 30 to 40% of states it is asked about; the answers are the same two numbers compared the
same way, so the output is unchanged.

Three things had to hold, and each was measured before any of it was written:

- **The bound needs less than 10 dB of slack**, or the peak search reaches outside the window
  and there is nothing left to save. The whole-sphere maximum is 11.4 dB, just over. Binning it
  by direction -- `DirectionBounds`, bins uniform in cos^2(theta) crossed with azimuth sectors
  -- brings the median to 8.6 dB.
- **The test must cost nothing on the states it does not reject.** It forms the squared range
  and the three dot products with the sphere frame, which is exactly what an evaluation needs,
  so `Database::level_from_geometry` takes those and a surviving state pays nothing for having
  been tested. Letting the test compute its own geometry and then recomputing it inside
  `level_at_observer` was measured: it turns a 1.48x into a 1.08x.
- **The table has to stay in cache.** One table per *distinct blend*, not per state -- a shipped
  track resolves thousands of states onto a few dozen tetrahedra -- so it is tens of kilobytes.
  Built per state it was 3.3 MB, re-streamed for every observer, and that alone cost more than
  the scheme saved once the loop was threaded.

The bound maxes over the four conditions of a tetrahedron rather than over the particular blend
of them, so it is looser than a per-state table would be; buying that back means a table per
state, which is the cache problem above. Bin counts were swept under the earlier window rule
(8x8 to 32x16) and the choice hardly mattered; 16x8 is kept.

The whole thing rests on one property, that the bound is never exceeded, and it is tested
directly (`tests/test_sphere_database.cpp`, `tests/test_spectral_propagation.cpp`, tag
`[bounds]`) over sweeps of source states and observers, **including sweeps aimed within five
degrees of the nose and tail axes**. A cell is filed under the bin of its center direction, and
at theta = 0 or 180 that direction has no azimuth, so a pole cell is filed under every azimuth
bin of its theta bin (`Database::home_bins`); before that, a query a degree off the axis was
bounded up to 0.8 dB too low and the ordinary sweep, which samples hardly any observers there,
never saw it.

Two suggested optimizations were deliberately **not** taken:

- Storing reciprocals of the (max − min) normalizing ranges: a reciprocal multiply is not the same as a divide, so it would move the goldens, and the normalization is per state now anyway.
- Hoisting the output-directory check out of `write_time_history_file`: one `exists()` next to creating and writing a ~130 KB file is noise, and it would mean bypassing `open_output`'s guarantee.

**The per-state slab was built and measured; do not rebuild it.** Collapsing the database to one blended 36×37 slab per trajectory state (with the observer loop blocked so `state` can be the outer loop) works and is bit-exact, but it buys only 6.8% of the acoustic loop for +131 MB of resident memory, and the condition-major grid layout it requires makes the ordinary blend path — `self_test`, and any observer set below the break-even block size — about 6% *slower*. Block size barely moves it, and a single 768 MB block is the worst configuration of all, because the level buffer's own bandwidth then dominates. The work is preserved on the `slab-experiment` branch with the full numbers.

The reason the intuition fails is worth keeping: observers sweep a 120×120 ground grid, so consecutive evaluations land on neighboring (phi, θ) cells and the gather is already cache-resident. It was never miss-bound. Measuring the ceiling first — point `eval` at a fixed L1-resident slab and time it — costs ten minutes and would have shown this before any restructuring.

**The SEL integral takes its time step from the trajectory, not from `-s`.** `sound_exposure_level` weights every sample by a time step, and `-s` is only the resampling step for `adapt()` — the two coincide *only* when `adapt()` actually resamples onto that grid, and it is commented out in `main`. Every shipped track is sampled at 0.02 s, so using the 0.5 s default as the weight overstated every level by exactly `10*log10(0.5/0.02) = 13.98 dB`. `main` now derives the step via `uniform_sample_interval()` and logs it, warning if the sampling is not uniform. If you re-enable `adapt()`, it resamples to `-s` and the derived step follows automatically.

- `src/sphere_database.{hpp,cpp}` — `Database` class. Loads a netCDF noise database (`.nod`): one group per operating condition (μ advance ratio, γ flight path angle, CT thrust coefficient), each holding scattered (phi, θ) channels of dBA level and EAA excess-attenuation factor. `eval()` computes ground level: geometry → sphere angles (`lookup_angles`), barycentric interpolation between operating conditions, periodic bilinear interpolation on the sphere grid, then spherical spreading + excess attenuation, `EAA * (range - r_sphere) / 1000` — EAA being
the total A-weighted attenuation over the 1000 m reference distance the database was built at,
which is exactly how panam writes and consumes it. Ground effect: see `--ground_effect`, on by default wherever the spectral path runs. Only `same_grid == true` databases are supported (the interpolation path is unimplemented). Free functions (`bracket`, `periodic_bilinear_interp`, `lookup_angles`, the geometry helpers) are exported in the header so the unit tests can reach them.
  **The loader's guards**: scalars and arrays are read as `double` whatever the file stores (`nc_get_var` copies the external type's bytes unconverted, so a float32 `.nod` loaded as a radius of 0 and NaN footprints with exit 0), and a "scalar" reads element 0 only; a sphere layout that leaves grid cells without a channel (they would read as 0 dBA) or drops channels that match no angle of sphere 0's axes is warned about once per distinct layout, first three (none of the 20 shipped databases trips it; `B407_M3`'s 934 zero channels are written zeros, not unreached cells); `sphere_frame` with no velocity, or velocity along the thrust axis, takes +x's projection for its nose instead of a zero frame that looked every observer up at (phi 0, theta 90); and the sphere writers check `nc_close`'s status, which is where HDF5's buffering reports a full disk. **An exported sphere carries the root the loader interprets** (`Database::Sphere::speed_reference`, `tip_speed`, `rotor_radius`, `vehicle_weight_newtons`, written at the root where the source had them), so it reloads as air- or ground-referenced as its source was; `--export_sphere`/`--export_aam` refuse an output that is the database being read (they write with `NC_CLOBBER`), and `interpolate_sphere` warns when the requested condition is outside the envelope (`blend_at` clamps it, and the file is labeled with what was asked for). `an exported sphere reloads as the sphere it was, cell for cell` and the AAM readback test pin the payload; swapped dBA/EAA or a dropped load-factor correction used to pass the whole suite.
- `src/predicates.{hpp,cpp}` — exact orientation and in-sphere tests via Shewchuk expansions, with a filtered double fast path. The recursion, `exact_determinant_recurse`, takes its matrix as an **array of pointers** and builds its minors the same way: the entries are expansions, which are heap-allocated vectors, and copying `(size-1)^2` of them at every level of a factorial recursion put `malloc` and `free` at the top of the database-load profile. Callers use the `std::vector<Expansion>` overload of `exact_determinant`, which adapts to the pointer form once; it asserts `size <= MAX_PREDICATE_SIZE` so the fixed stack buffers the recursion uses guard themselves rather than relying on each entry point's `fatal()`. The exact path is taken often here — mass degeneracy is the point — so it is worth keeping allocation-free. **The exactness is not optional**: the operating conditions are a prism (164 flight conditions x 8 thrust coefficients), so whole families of points are *exactly* cospherical and a rounded in-sphere test answers the same question two different ways, which leaves Bowyer-Watson's cavity non-star-shaped. Assumes normal (non-subnormal) coordinates; conditions are normalized into [0, 1] first.
- `src/delaunay.{hpp,cpp}` — Delaunay triangulation in arbitrary dimension, Bowyer-Watson on those predicates. Dimension is a constructor argument, not a template parameter, so a fourth operating parameter costs nothing here. A simplex whose circumsphere passes exactly through the new point is left standing: a point *on* a circumsphere does not violate the Delaunay condition, and leaving it is what keeps the choice consistent under mass degeneracy.
- `src/interpolator.{hpp,cpp}` — barycentric weights over the operating-condition space, on top of `delaunay.hpp`. Replaced CGAL's `Delaunay_triangulation_3`, and Boost/GMP/MPFR went with it. Databases with fewer than `dim + 1` conditions (a single exported sphere) skip the triangulation and take the nearest condition.
- `src/atmosphere.{hpp,cpp}` — ISO 9613-1:1993 atmospheric absorption, ported term for term
  from panam's `panam_acoustics/iso_9613_1_1993.py` and pinned against values generated by
  running that module (`tests/test_atmosphere.cpp`). **Pressure is in kPa**, matching the
  standard and the Python, and is the one non-SI unit in the codebase. `absorption_coefficients()`
  is the form the hot path uses: it hoists the two relaxation frequencies out of the per-band
  loop and is bit-exact against the per-frequency call.
- `src/propagation.hpp` — the propagation seam. `straight_ray()` returns a `RayPath`: the
  direction the sphere is looked up by, the range that spreading and absorption take, and the
  ground model's (source height, ground distance). `spectral_level_at_observer` takes all of
  these from it and nothing from raw geometry. It went in bit-identical on every path. A
  refracting model returns the same record per eigenray. The header lists what such a model
  must *also* change and that deliberately does not live there: reception time
  (`reception_factor` differentiates straight-line range), the bounds (they subtract
  straight-line spreading; focusing adds gain), and the broadband path, which shares its
  geometry with the bound test and does not use the record.
- `src/atmosphere_profile.{hpp,cpp}` — `--atmosphere`: a stratified T/RH/p profile and
  `LayeredAbsorption`, per-band absorption integrated through it. Along a straight path a band
  loses the coefficient's mean over the heights spanned times the path length. The mean is
  (F_b(z_s) - F_b(z_o)) / (z_s - z_o), F being a 1 m table of alpha's running integral. The
  things to preserve:
  - **The hot loop pays nothing on a grid, and a subtract and a multiply per band otherwise.**
    F at each state's height is looked up once per trajectory. When every observer shares a
    height, as on any grid, each state's mean is formed once too (`state_absorption_mean`), and
    the loop only reads it. Without that the profile cost ~5% of the acoustic loop. Scattered
    observers at mixed heights take the per-pair path; `atmosphere_uniform_scattered_is_homogeneous`
    runs `arrays/coyle_approach.csv` through it.
  - **`Database` knows nothing about profiles.** The caller hands `spectral_level_at_observer`
    per-band coefficients exactly as before. Without `--atmosphere` the run-wide array goes
    through untouched, so the default output is byte-identical.
  - **F is quadratic within a node.** Alpha is linear between nodes and F is its exact
    integral, so the difference quotient tends smoothly to alpha as the heights meet. Within
    1 mm it switches to alpha at the midpoint.
  - **One level means uniform at every height.** With several levels, a trajectory above the
    top is refused, not extrapolated.
  - **The surface level supplies the sound speed** for the reception time and the ground
    impedance.
  - **Spectral path only.** `EAA` is a fitted broadband attenuation with nothing per band to
    integrate.
  - **`atmosphere_uniform_is_homogeneous` compares at 1e-9, not byte-for-byte.** A mean of a
    uniform integral is the coefficient to a few ULP. Do not "fix" it by comparing bytes.
  - **Wind columns are ignored** and named in the log.
  - **`--wind` builds the profile instead** (`synthetic_profile`, `SurfaceWind`): a neutral log
    law through one surface observation, held above `--surface_layer_depth`, a dry-adiabatic
    lapse from the run's surface air, hydrostatic pressure, water vapor held. On the held-out
    windy day 284 (stratified + plate) it scores 1.60-1.62 dB rms against the measured LIDAR and
    balloon profile's 1.69 and straight rays' 2.35; calm day 286 1.25 against 1.26. The things
    to preserve:
    - **Never emit a calm level above z = 0.** `--refraction stratified` takes one (rows at 0
      and 0.1 ft, both still) as a thin layer that put levels 22-28 dB low. The log law is zero
      at z0, so the first level above the ground is ten roughness lengths up (or the anemometer,
      if lower), and `synthetic_profile` throws `logic_error` if a level there comes out calm.
      Do not "refine near the ground" by adding levels at or below z0. A calm observation gives
      a profile with *no wind columns*, not columns of zeros; `refer_to_database_speed` treats a
      `--wind` of zero as calm stated, so `--speed_reference air` is not refused for it.
    - **Built per trajectory, in `process_trajectory`, not in `run()`.** Its surface air is that
      trajectory's vehicle file's (under `--temperature` and friends, which `--wind`, unlike
      `--atmosphere`, accepts) and its top has to clear that track. Inside `process_trajectory`
      the profile is the local `atmosphere`, never `options.atmosphere`, which is empty under
      `--wind`: anything new that reads the profile there has to take the local.
    - **The top is max(surface layer, 1.05 x highest path), and no lower.** Above both, the wind
      is held and the air lapses, so c + w.h falls with height in every direction and no ray
      returns from air the profile leaves out. A stable class's inversion runs to the top of the
      surface layer, which is why the top never stops at the trajectory.
    - **The anemometer, the top of the surface layer and the top are exact levels.** The
      profile passes through the observation, and the kink at the layer's top is a level rather
      than smeared into whichever layer straddles it.
    - **Discretization is checked through the real run** (`wind_is_discretization_independent`:
      16 -> 64 levels a decade, plate microphones over the Coyle array, 0.016 dB of a 0.1 dB
      allowance). It is on plates and not on a point-receiver grid for a reason: grid cells
      1500 m out downwind move +-0.2 dB with *any* change of layering and do not settle from 8
      levels a decade to 100. That is the stratified arrival search, whose samples follow the
      layers, not the profile. Do not tighten the tolerance onto that floor, or chase it here.
    - **`--write_atmosphere` replays to 1e-9** (`wind_is_its_written_profile`): 17 digits, the
      wind as u/v components rather than speed and direction (no trigonometric round trip), and
      the pressure written, so the replay does not re-integrate it. A replay compares the run with
      itself, so it is blind to the plumbing in `run()` and `wind_profile`: `wind_profile_is_its_law`
      holds a written profile to an independent closed-form oracle
      (`tests/data/wind_profile_expected.{py,csv}`, 1e-12) and fails on a length left in feet, a
      dropped `--temperature`/`--pressure` override, or `--profile_levels` not passed on. Change the
      generator on purpose and regenerate the CSV from the script, not from the binary.
    - **The bottom layer, 0 to 10 z0, is linear and `--profile_levels` does not refine it.** Filling
      it with the log law moves day 284 by 0.15 dB at most (0.014 rms). Do not "fix" that by
      starting the levels nearer z0: those levels carry almost no wind, the thin layer above.
    - **The stability classes are unvalidated, and warn.** TA Luft's exponents and RG 1.23's
      gradients, as SELENE has them. Class II scored 2.22 dB on day 284, IV 1.57; one day picked
      after the fact chooses nothing. Do not make one a default without measurements of the class.

  The roadmap this starts (refraction, then urban) is in the Claude Doc "NICE-OPS:
  Atmospheric and Urban Propagation Options"; `reports/` holds its 2026-09-26 snapshot, and
  `reports/Urban canyon propagation methods.md` works out the urban step.
- `src/refraction.{hpp,cpp}` — `--refraction linear` (roadmap step 4). Per path, a
  least-squares line through the effective sound speed c(z) + w(z).h over the heights the path
  spans (`EffectiveSoundSpeed`: running integrals of c, w_x, w_y and z times each on a 1 m
  table, so a fit costs a few flops per pair), then closed-form circular-arc eigenrays. The
  things to preserve:
  - **Numerically stable forms, not textbook ones.** The circles are centered L = c0/|a| away,
    kilometers to 1e11 m, so every difference of squares of h = L + sigma z is formed from
    kappa(z) = z(2L + sigma z). Reflected legs are kappa/(sum of d), travel time is one atanh of
    the addition formula, and arc length is 2 rho asin(chord/2 rho). Written the obvious way
    they cancel to nothing as the gradient shrinks. Gradients under 1e-9 /s take the straight
    path.
  - **Spreading is r / cos(launch), exactly.** In a linear profile the tube area
    r |dr/dtheta| sin|arrival| / cos(launch) reduces to (r rho / h_s)^2.
    `tests/test_refraction.cpp` checks it against finite differences of traced rays; keep that.
  - **The reflected ray is a 1-D solve** on launch angle, along which the range is monotonic
    (checked across both gradient signs; for a receiver *above* the source under downward
    refraction it is monotonic on each side of a window |theta| < acos(h_s/h_r) whose launches
    never climb to the receiver, and the solve takes the side the range at the window's edge
    says, or reports none: a root of the clamped continuation inside the window is not a ray,
    and one past the ray's apex is a second, descending crossing, which this is not). It starts from the direct ray to the receiver's image
    below the ground and takes secant steps from there, which usually close in one or two. That
    warm start took the point-receiver refraction run from 4.7 to 3.6 s. A bracketed Illinois
    solve is the fallback. Under downward refraction at long range the ray is launched
    *upward*, over its top: a separate branch of the range formula.
  - **Ground models get an equivalent straight geometry.** A point receiver gets
    h' = L2 sin(grazing) - h_r and d' = L2 cos(grazing) from the reflected ray; a plate gets
    the direct ray's arrival angle and length. The ground models themselves are untouched.
  - **The direct ray ends at the microphone**, the same point the straight path uses (see
    Receivers below), never clamped. Scattered flush mics on sloping ground sit below the
    plane, and clamping them once made plate runs disagree with the straight model. Uniform
    air reproduces the straight run to 1e-9 for plates, point receivers and bare observers
    (`refraction_uniform_*_is_straight`). The reflected ray goes to the microphone's height
    above the plane z = 0, where the ground models put it.
  - **The same bounds cover refraction.** The peak search and `--db_cutoff` take a refracted
    state's geometry as its launch direction scaled by its ray tube's range, so the existing
    DirectionBounds tables bound it unchanged. Absorption and the shadow's loss are losses, and
    the ground's gain is bounded independently of geometry. `direct_ray` gives that geometry
    and the travel time; `finish_ground`, with its reflected-ray solve, runs only for states
    that survive the bound. `refraction_*_windowed_is_exact` checks the windowed run against a
    march. On this database and track the bounds are loose even for straight rays (92–98% of
    states get evaluated), so the windowing saves little here. The overhead of refraction is
    the rays themselves: +0.7 s bare, +1.5 s with a point receiver, +0.7 s with a plate, on the
    demo grid.
  - **Every ray is checked against RK4 traces** of the ray equations: position to 1e-6 m,
    time and arc length to 1e-9, reflected grazing angle to 1e-8 rad.
  - **Wind directions are compass bearings** turned by `--frame_bearing` (default +x east).
    The shipped AC tracks need 210.
  - Shadow zones: the grazing ray's level, less the least attenuated creeping-wave mode,
    capped at `--shadow_floor`. `ShadowAttenuation` builds it: over rigid ground (no
    `--ground_effect`) it is `shadow_attenuation_rates` and `shadow_loss`'s old arithmetic,
    to the bit; over the ground's admittance it is one table per band in s = (pi f / |a|)^(1/3),
    of the first zero of Ai' + q Ai followed along q = -i w^2 beta s. Three things to keep:
    - **It is not monotonic in the ground.** The first-order shift -q/a'_1 lowers the rate
      for arg beta in (-90, -30) degrees, so ordinary grounds attenuate up to 12% *less*
      than rigid at small |q| before up to 2.3x more. Do not "fix" a soft ground coming out
      quieter than rigid in the shallow shadow.
    - **Bands whose Re q > 0 (arg beta < -60) take the rigid rate.** There the first zero
      turns into the surface wave and the follower loses it (measured: it runs to |t| ~ 100,
      near q^2, and throws). Only strong porosity rates and thin layers get there, at low
      frequency; the run logs the bands.
    - `tests/test_shadow_attenuation.cpp` holds the rate to the exact residue series' slope
      deep in the shadow, 0.02-0.3% on the benchmark grounds. Keep it.
  - **`ShadowTransition`: the level near a shadow boundary** (`--refraction linear`, point
    receivers only). Ray theory fails within a couple of diffraction lengths F = 2 k0 l^2 of
    the boundary; its ground term at the grazing geometry was 10-25 dB low over grass and
    ~5 dB high over rigid ground there. Per band and microphone height, a table V(xi, s) of
    the residue series relative to the ray's own spreading, xi = (r - r_b)/F. Things to keep:
    - **Built at a reference source height, 20 l** (0.1 s^2 when lower, keeping it under
      R/10). V collapses over source height (within ~1 dB, 5 l-100 l, from xi = -1.5) but
      the series itself fails for high sources at high frequency, so do not build it at
      the actual height.
    - **Each model at its own boundary**: the table places the exact field at the parabolic
      boundary of k^2-linear and the ray at the circle one; NICE-OPS uses the ray's.
    - **Ray model below xi = -2, blend to -1, V beyond**, floored at V(0) - floor and capped
      at 20 log10 2. `max_power_ratio` feeds the peak search's ground bound, or the windowed
      run stops being exact.
    - **A band whose boundary is nearer the source than two diffraction lengths keeps the
      ray model to the boundary**, blends to V over the next diffraction length, and takes V
      beyond (`cramped` in `apply`). There the source is a few l up at most, under the 5 l the
      collapse was checked to, and the table has no lit zone to hand over from. A 30-50 ft
      source under a weak upwind gradient has it at every band to 500 Hz: V put the whole
      path 2-5 dB under the PE, which matched straight rays. The blend leaves the creeping-wave
      loss out: it is valid only deep in the shadow, and near the boundary it put a 30 m
      source in the linear lapse 6 dB under the PE at 250 Hz. Ray model into the shadow too was
      worse deep in it (day-286 PE matrix, 50-160 Hz: 5.7 dB rms against 5.1), and so was the
      ray model with only V's fall past the boundary. Cost: that lapse case is 0.1 dB rms worse
      at 200-630 Hz (its old lit side was 1.3 dB low, its old shadow right); day 284's PE cases
      do not move.
    - `Database::spectral_level_along`'s `ground_weight` is how the ground term is faded out:
      a power of the power ratio, per band.
    - **Through layers, each band gets its own boundary**, r_b = r_s + w (r_l - r_s),
      w = x / (x + boundary_resolution), x the band's diffraction length at the fit's
      curvature over r_l - r_s (`band_boundary`). r_s is the layered boundary, r_l the linear
      fit's, and both ride on `DirectEigenray`. Scaled by r_s alone it fell 10-20 dB too
      fast at 125-250 Hz against the PE, because low bands do not see a log law's thin fast
      layer. The constant 0.363 was fitted on 12 log laws and held out on example_windy
      (6.8 dB edge rms, against 18.9 at r_s). Pooled over 16 log-law cases the edge error
      goes 30.7 -> 6.1 dB, and the lapse control stays at 0.4.
    - **Flush microphones take their own constant, 0.75** (`boundary_resolution_flush`,
      chosen by the table's receiver height or plate, never the observer's terrain z). Fitted
      per case against the PE, a flush microphone's best boundary is ~0.9 diffraction lengths
      at the fit's curvature past r_s on the measured day-284 profiles and ~1.2 on flush log
      laws, where 0.363 put it at 1.7-2.7. 1.2 m keeps 0.363 (0.387 refitted).
    - **Where the fit casts no shadow (r_l infinite), each band takes a multiple of r_s**:
      1.60 (f/1 kHz)^-0.08 flush, 1.90 (f/1 kHz)^-0.08 at 1.2 m (`fallback_multiple_*`).
      The law has nothing to interpolate toward there, and falling back to r_s itself was
      the cliff again. A strong inversion above the source does it: the fit over the path
      refracts downward while the shear under it still casts a shadow.
    - **The law cannot be made much better by fitting.** The per-case optimum is 2.7 dB
      pooled; every smooth law tried (constant multiples, power laws in f, offsets in
      diffraction lengths, tangent-extended profiles, regressions on f, z_s, r_l/r_s) sits
      at 6-9 over all cases. The optimum is sharp where the wall is steep. Do not spend another day on it
      without turbulence in the reference.
    - **Past 0.9 r_s the lit ray is the linear fit's** (`handover_fraction`, stratified.hpp).
      `layered_direct_ray(..., true)` hands it back (`from_linear_fit`, finished by
      `finish_ground`) whenever r > 0.9 r_s, *not only where the layered solve fails*, and
      keeps it past the fit's own boundary too, in its shadow. Short of r_s the layered ray
      focuses into a caustic (+10 dB) that the PE does not show; 0.9 was the best hand-over
      on both calibration sets (0.6 dB rms on day 284 over [0.5, 1] r_s, 1.4 at r_s). Going
      back to the layered shadow ray past r_l stepped 2-10 dB, where the fit's limiting ray
      is continuous. On a log law the dipping family's range is not monotonic
      and reaches tens of meters past its end, and those near-grazing rays have no
      reflected partner, so the ground term took the grazing geometry (-32 dB). It is
      formed in `form_ray`, so the peak search bounds the ray whose level it takes. Only
      for observers the transition serves; plates and no-ground runs are byte-identical.
    - **The calibration was confounded once, so beware.** An experiment that sets a
      boundary by hand must also set `linear_boundary` equal to it, or `apply()`'s own
      interpolation runs on top. That run falsely suggested w = 0.
    - Against the PE on a linear lapse: shadow-edge error 5-20 dB -> 0.1-1 dB.
      `tests/test_shadow_transition.cpp` holds it to the exact series, and pins the law
      against stored PE levels (example_windy, 100 m source, 500 Hz and 1 kHz).
    - **Plates** get their own table: the field at the plate's top (its thickness as the
      receiver height) plus K, the plate over bare ground, per band and azimuth, taken from
      `PlateReceiver::level` at 2 degrees. K barely depends on elevation near grazing
      (within 0.1 dB from 0.5 to 10 degrees), which is what makes one value hold through
      the shadow. The caller passes the pair's plate-frame azimuth, the same as
      `Database::spectral_level_along` forms it. `max_power_ratio` includes the loudest
      K. On the 2017 day 284 this takes `stratified`'s 5-10 degree rms from 4.5 to 3.6 dB.
    - The tables take the run's `-j` (`threads`). The harness runs a dozen `-j 1`
      processes at once, and a table building on every core in each of them took 7 s.
  - **`layered_direct_ray` solves for a microphone below the plane as if on it.** The
    layers hold their lowest speed below z = 0 and the direct-ray solve knows no ground,
    so it sent a ray grazing along under the plane: a "lit" arrival of -0.003 degrees
    whose tube was 46 times its range, 33 dB low. Ground boards 3-10 cm under the survey
    datum hit it on every run of the 2017 data. The ray still ends at the microphone,
    and uniform air never reaches the layered solve, so the straight-ray equivalence
    stands.
  - **The layered shadow's tube can be 0/0** for a receiver a few millimeters up, where the
    dipping family is a sliver that thin. It then falls back to r / cos(launch). Both
    fixes have tests in `test_shadow_transition.cpp`.
- `src/stratified.{hpp,cpp}` — `--refraction stratified` (roadmap step 6a): the effective
  sound speed piecewise linear between the profile's own levels, circular arcs within each
  layer, eigenrays solved per pair. The things to preserve:
  - **A ground-grazing boundary needs only that no ray from aloft reaches the receiver**:
    the air above the source must not outrun everything between receiver and source
    (`fastest_above(z_s) > fastest(z_r, z_s)`, the same rule families() uses for the
    upward arms). Testing against the source's own speed left the day-284 mornings -- an
    inversion over an upwind shear, 330.5 m/s aloft against 339 at the ground -- with no
    boundary, so no transition and no hand-over, and the shadow a 10-35 dB cliff the PE
    puts 30-50% further out. It moved day 284's upwind 5-10 degree bin most of all; see the
    manual (`docs/MANUAL.md`) for the table.
  - **And where the air aloft does outrun it, the boundary stands short of the skip
    distance** (`boundary_under_aloft`). Rays turned back from aloft come down only past it,
    kilometers out for an inversion a few hundred feet up; short of it the shadow is the same
    ground-grazing one. `layered_direct_ray` takes the boundary back where a ray from aloft
    or along a duct does arrive in the shadow, which leaves the duct's promotion as it was.
    Without it a 30 ft source hovering upwind on day 286 got no boundary 600 ft out: no
    transition, no hand-over, the layered shadow ray's tube focused 6.5 dB by the log law's
    lowest meters, and the whole ground term -- 12 dB over free field at 63 Hz where the PE
    gives 4. On day 286's measured frames (63-630 Hz) the states under 100 ft and 20 kt went
    from +8.9 to +4.4 dB against straight rays' +3.8, so what is left there is the source; no
    other height, speed, range or elevation cell moved by more than 0.2. On the day-286 PE matrix
    past 1000 ft it took 200-630 Hz from 9.5 to 7.5 dB rms and 800-2000 Hz from 11.7 to 9.2.
  - **Parametrize by the Snell invariant p, in four families.** Launched up or down, arriving
    on the first crossing of the receiver height or on the second, after dipping below it
    (`Arm`, `families`). The range is monotonic in p within a family. The dipping family is
    what the linear model's unique circle includes implicitly. Without it the shadow boundary
    comes out short by sqrt(2 L z_r), which was 170 m on a linear test profile.
  - **Layer crossings are closed forms written stably**
    (dx = p(c0+c1)dz/(q0+q1), plus log1p and atan2 forms for time and arc). They were checked
    against a 200,000-point quadrature to 3e-15. At a turning point the sine is pinned to 0:
    computing it as sqrt(1 - p^2 c^2) leaves ~1e-8 of rounding, which is 1e-7 in range.
  - **The RK4 trace in `tests/test_stratified.cpp` is fussy, for reasons worth knowing.** A
    finite-difference gradient rounds to 1e-8 /s, which is 1e-4 m over a km of ray. A step that
    straddles a level integrates across the gradient's jump and drops to first order. Even a
    step that only *lands* on a level takes its last RK4 stage from the next layer unless the
    step is held to the gradient of the layer it starts in. With all three fixed, the trace
    matches the layered model to 1e-10 m.
  - **The linear model stays alongside.** It gives the warm start (usually within a fraction
    of a degree), the fallback for a receiver at or above the source, the shadow's attenuation
    rate (its fitted gradient), and the refracted path's `profile`.
  - **Every arrival, summed in energy** (`layered_arrivals`, `Arrival` in refraction.hpp).
    The primary is the first family's root. The other families are not monotonic in general,
    so they are sampled in launch angle (32 + 8 per hop the path can hold, at most 512) and every
    crossing of r is solved. In a surface
    duct (air above the source faster than all of it below) the bounced rays come too:
    X_m(p) = X_0(p) + m L(p), from the same samples. Each earlier bounce costs the plane-wave
    |R|^2 (`add_bounce_loss`), and past every direct ray's reach the strongest is promoted to
    primary. Day 286's inversion was 8–24 dB under the PE past 4,000 ft without them, and is
    now mostly within 3 dB. Focusing is capped per ray at `max_focusing_gain_db` (+10 dB); distinct arrivals
    add in energy and can exceed it in total. All solved arrivals are retained. Things to keep:
    - **`solve` refuses a bracket that does not land** (residual over 1e-6 r). A family whose
      end ray never arrives stands in as an infinity, so it brackets falsely. In a duct that
      reported the duct's top as a lit ray, 55 dB low.
    - **`Jumps`: the range is discontinuous** wherever 1/p passes a local maximum of c that
      is the fastest air yet met above z_s, or below z_r. There the turning point leaps.
      `solve` splits its bracket there. Before that, a sign change at a jump alone ran the
      full 200 iterations to be refused, under a convergence test (`hi - lo <= 1e-17 hi`)
      that a double cannot meet. That was a third of the solves in a duct and 90% of their
      cost; the test is now that the estimate no longer falls inside the bracket.
    - **The samples are shared across families** (`Row`: down, up, dip, ground pieces at one
      p). Per family they repeated the walks 3–5 times in a duct.
    - **The sample count scales with the hops, and 32 is not generous.** A ray that bounces m
      times travels X_0(p) + m L(p), so the same interval holds m times the wiggles: aliasing
      in the hop count, not near-zero intervals (refining those helped 30%). Against a
      2048-sample reference (12 bearings x 8 source heights 5-300 m x 2 receiver heights x 120
      ranges to 10 km), pairs more than 1 dB out: inversion 11.3% at a fixed 32, 2.6% at 128,
      0.47% at 512, **2.9% at `32 + 8*m_max` (cap 512)**; windy 1.05%, 0.27%, 0.07%, **0.01%**.
      The grid is bunched toward the shallow end (t = t_a + (t_b - t_a) u^2), free, and worth
      most of the windy gain; the inversion's remainder is what is left at that cost (`+16/hop`
      gets 1.8% for another 20% of the arrival time, a cosine bunching 2.6% at +8). It costs
      4-8% of a stratified run. Against the PE (4 bands, 12 geometries, beyond 1000 ft) the
      inversion's agreement moved 0.1-0.2 dB rms for the better on average (100 ft source:
      4.1 -> 3.6, 3.8 -> 3.3, 3.4 -> 2.7, 3.6 -> 2.9; 30 ft: 5.1 -> 5.2, 3.5 -> 3.6, 2.5 -> 2.6,
      2.7 -> 2.8), the windy's not at all. `layered_arrivals`' `samples` argument exists for
      `the arrival search has converged in its sample count`. A base of 16 is far worse. Do not
      trade samples for speed.
    - **The primary's bracket is split around it, not solved whole.** Whole, it returned the
      primary or one of the caustic pair beside it, at random.
    - **Past the hand-over to the linear fit only bounced rays, and rays turned back from
      aloft, are searched** (`bounced_only`). The fit stands for every downward direct ray;
      rays from aloft exist there only under a `boundary_under_aloft`. Upwind under a lapse
      that skips the search entirely.
    - **The peak search bounds the sum** by the loudest arrival's direction over the harmonic
      sum of their squared ranges (`bound_geometry` in niceops.cpp).
      `stratified_duct_*_windowed_is_exact` holds it.
    - **Not fixed: an elevated duct whose returning rays turn just above the ground** (day
      286, run 286288, 100 ft source toward 270). It is not warm ground air: the air is
      fastest at z = 0 because the log-law headwind assumed below the 46 ft LIDAR gate takes
      0.2 m/s off in the lowest foot (the temperature rises there). Rays turned back from the
      ~330 ft inversion turn up again 1–20 ft up and reach neither a flush nor a 1.2 m
      microphone: 11–21 dB under the PE at 8,000–13,000 ft. The PE tunnels through: holding
      the wind below 46 ft flat moves it under 1 dB at 63–250 Hz and ~6 dB at 1 kHz, as a WKB
      factor exp(-2 int kappa dz) does for rays turning at 1–9 ft. The skip zone short of the
      duct's return is a second gap: with the shear removed, `stratified` sits 30–60 dB under
      the PE there and finds the return 1–2 kft late. On measured day-286 frames where
      `stratified` is >= 10 dB under straight rays (1.2% of frame-bands; 0.4–2 deg, 4–9 kft,
      sources under 300 ft) it is 10–20 dB too quiet and straight rays are within 4 dB at
      50–630 Hz; none of those frames is in an SEL window. Harness
      `propagation/micro_shadow.py`, `micro_shadow_frames.py`.
    - **Three opt-in treatments of it, compared** (`docs/MANUAL.md`, "Micro-shadows and skip zones";
      harness `micro_shadow_fixes.py`, `micro_shadow_frames.py`): `--tunneling` and
      `--skip_fill` (`layered_wave_arrivals`, stratified.cpp), and a PE `--propagation_table`
      beside the ray model. Both Airy options together take the deficit frames from -10.5 /
      -15.5 to -2.1 / -5.0 dB at 50-160 / 200-630 Hz and barely touch 800-2000 (-20.6 ->
      -18.9). The PE table puts the low sources (under 150 ft, 3-6 kft) within +0.3 to +2.5.
      The SEL moves under 0.02 dB under every option.
    - **Wave arrivals carry their loss per band (`WaveCorrection`, `add_wave_loss`)** on a
      tube taken at `wave_reference_frequency` (10 kHz), so the loss is never negative and
      the peak search's bound on the tube still bounds them. They stay out of the
      boundary's removal and the promotion to primary in `layered_direct_ray`, which only
      ray arrivals take. A tunneled one takes no ground term, plate or not.
    - **Near a tunneled family's caustic the coalescing pair is replaced by the caustic's
      band-averaged Airy field.** Left to the +10 dB focusing cap it put the return 5-12 dB
      over the PE at 250 Hz-2 kHz.
    - Cost: about +55% on a stratified run (demo grid, B407, `--ground_effect`: 23.6 -> 37 s
      on `example_windy.csv`, 57.5 -> 86 s on `example_inversion.csv`). It was +107% and +76%
      before the jumps, the shared samples and the hand-over skip.
  - **The log-law finding.** With zero wind at the ground, the upwind limiting ray is fixed by
    c(z_s)/c(0) alone. On `example_windy.csv` a 100 m source's upwind shadow starts at 520 m,
    where the linear fit puts it past 1.2 km; a test pins that. Ray theory probably overstates
    that shadow at low frequencies, so GFPE (step 7) should adjudicate before either model's
    shadow levels are trusted.
  - Costs: 12 s on the demo grid with a point receiver (linear 3.3 s), 35 ms over a
    96-microphone array.
- `src/pe.{hpp,cpp}`, `tools/pe_compare.cpp` — roadmap step 7a: a wide-angle parabolic
  equation as the full-wave reference, and `niceops_pe_compare`, which runs one path through
  straight, linear, stratified and PE. Load-bearing:
  - **A finite-difference PE, not the GFPE the roadmap named.** The impedance condition
    dpsi/dz = -i k beta psi goes straight into the grid (psi_0 eliminated by a one-sided
    difference), so no surface-wave pole has to be handled by hand. The convention is
    exp(-i omega t), as in `ground_effect.cpp`, whose admittance it takes. When speed became
    the question the answer was again not GFPE, which would have brought that pole back and
    splits refraction from diffraction (a step limit exactly where the gradient is strong),
    but a longer step for the same grid: the split-step Padé march below.
  - **The default march is split-step Padé** (`PeScheme::split_step_pade`,
    `src/pade_propagator.{hpp,cpp}`): psi(r + dr) = [c_0 + sum_j c_j (Q - p_j)^-1] psi(r), a
    rational stand-in for exp(i k_a dr (sqrt(1 + Q) - 1)). One tridiagonal solve per term,
    about ten terms for an eight-wavelength step, where Crank–Nicolson takes eighty
    steps. 39x per march at the defaults; a five-band `pe_compare` case's PE 25.4 s ->
    0.47 s, the benchmark tool 2.5 min -> 10 s. Things to keep:
    - **The coefficients are fitted at run time, not tabulated and not Collins' Taylor
      Padé.** The step has to divide `output_step`, so k_a dr is different for every
      frequency. r is the type (n, n) interpolant of the exact propagator at Chebyshev nodes
      on the interval of x the run needs: elevation to `max_angle` (50 degrees), shifted by
      the profile's range of k(z), plus room for the surface wave. A Taylor Padé is exact
      at x = 0 and degrades outward, which is the wrong shape for that.
    - **Two of the nodes sit among the evanescent modes (x = -1.5 and -6), and they are
      not optional.** Interpolating a function of modulus one gives an all-pass r, |r| = 1
      along the whole axis, whose phase is right only on the interval. Steeper modes then
      travel upward at a twentieth of their true speed, linger, and come back as several dB
      of noise (3-5 dB rms on the rigid still-air test, which is how it was found). Pinned
      out there r falls away below the interval and they die in a few steps. Two is the
      number: three or four overshoot, |r| > 1 by 1e-4 to more than 1 between the interval
      and the first node. This is Collins' stability constraint.
    - **Nothing in the construction guarantees stability, so it is measured.** Every fit is
      checked on dense samples: error on the interval, |r| <= 1 over the whole spectrum
      including the absorbing layer's imaginary parts, every pole below the axis. The fewest
      terms that pass are used; if none up to 16 does, the step is shortened. Do not remove
      the check to save the milliseconds.
    - **Results are never interpolated in range.** The step is `output_step` over a whole
      number, so the march lands on every output. Over eight wavelengths psi turns through
      radians and linear interpolation is meaningless.
    - **A path steeper than `max_angle` from the source is NaN**, because those modes are
      absorbed and the number would be hundreds of dB down. Crank–Nicolson reported a
      plausible-looking level there with the wrong phase. `max_angle` stops at 60 degrees
      (the interval is clamped at x = -0.75, short of the branch point at -1).
    - **Wide-angle accuracy took three things, and the propagator was the least of
      them.** With `max_angle` alone raised from 35 to 50 the error against the exact
      two-ray field at 30-50 degrees did not move (0.6-1.7 dB rms at 1 kHz). Doubling the
      grid with the Gaussian starter barely moved it either. What it took:
      - *The self-starter* (`PeStarter::self`, Collins 1992): the first step is the
        propagator times (1 + Q)^(-1/4), applied to a delta function shared between two
        grid points, psi = sqrt(2 pi i / k_a) (1 + Q)^(-1/4) exp(...) delta. It is fitted
        and checked exactly as the propagator is (`PadeRequest::starter`), except that
        its modulus is bounded by twice its own value at the interval's low end rather
        than by one: it rises toward the branch point and is applied once.
      - *Fourth order in height* (`vertical_order` 4, Numerov). (Q - p) y = psi becomes
        [a delta^2 + (1 + delta^2/12) D] y = (1 + delta^2/12) psi; in y' = c y, with
        c = a + D/12 the column's off-diagonal, the off-diagonals are one again and the
        one-number-a-row factorization still serves. It costs a weight per term per row
        and a third more time. Second order at 20 points a wavelength is a radian of
        phase out by 35 degrees at 1 kHz from a 100 m source.
      - Together, and with the grid's edges at fourth order (next): 0.07 dB rms or
        better at every elevation to 60 degrees at ten points a wavelength, 0.01 to 50 at
        twenty.
    - **The ground condition and the interpolation are fourth order too**
      (`boundary_order` 4), and that is what makes ten points a wavelength the default.
      With the operator alone at fourth order the grid still converged at second. On the
      windy profile against 80 points: second order throughout 0.02-0.04 dB rms at 20
      points, fourth-order operator with second-order edges 0.01-0.02, fourth throughout
      under 0.0005, and 0.001-0.002 at 10 (0.004-0.01 at 6). How:
      - *The ground row comes from the equation, not from a wider stencil*, which would
        leave the tridiagonal. With F = y'' = k_a^2 (psi - D y) and y' = lambda y at z = 0,
        y_1 - (1 + dz lambda) y_0 = dz^2 (7 F_0 + 6 F_1 - F_2) / 24 + O(dz^5). It is a row
        in y_0, y_1, y_2 that is linear in the pole, like Numerov's rows, so every solve is
        still a resolvent of one operator -- on the grid that includes psi_0.
      - *psi_0 is therefore carried, as `psi[m]`*, not derived from psi_1 and psi_2. Each
        step rebuilds it from the two rows above and the old field, a few operations per
        term; the screens turn it with the rest; the Gaussian starter sets it.
      - *Row 1's right side differs per term* (it loses c_0 H R, R = 7 psi_0 + 6 psi_1 -
        psi_2), so the forward sweep takes row 0 a term at a time before the paired loop.
      - *Cubic interpolation at both ends of the path.* The receivers read four nodes,
        and the self-starter's delta function is spread over four with the same weights:
        by reciprocity a source placed to second order is a receiver read to second order.
    - **`points_per_wavelength` 0 means 10 where all of that holds and 20 where anything
      is second order** (Crank–Nicolson, `vertical_order` 2, `boundary_order` 2). Do not
      run a second-order configuration at 10: it is 0.06-0.19 dB rms and 0.6 at worst.
    - **The starter's first step is taken in parts if its fit is refused.** Its fit is
      harder than the propagator's (the amplitude factor rises toward the branch point)
      and fails at long steps on the coarse grid. It used to fall back to the Gaussian
      silently, which showed as a 0.007 dB offset that depended on the range step; now it
      takes the starter over a half or a quarter of the step and that length's own
      propagator over the rest. `a long first step is still self-started` pins it.
    - **The Gaussian starter is biased for any source where k differs from k at the
      ground**, because it is normalized with k_a: 0.11-0.15 dB high at every range for a
      100 m source in the benchmark's upward column (against the residue series), 0.10 dB
      on `example_windy.csv`. Every PE level from before the self-starter carries that,
      so do not take a 0.1 dB offset against a stored run for a regression. The
      self-starter is within 0.04 dB. A source within a grid point of the ground, or a
      starter fit that is refused, falls back to the Gaussian; `PeResult::self_started`
      says which ran.
    - **Compare a turbulent run with a still one marched the same way.** Turbulence caps
      the step at the screen spacing, and the field at a point whose reflected ray is
      steeper than `max_angle` depends on the step. `pe_turbulent_mean_field_decays`
      takes its still field with `gamma_t = 1e-30` for that reason; with 0 the two
      differed by 8% before any turbulence.
    - **The result does not depend on the step**: 2 to 32 wavelengths moves it by under
      0.001 dB, 60 dB into a shadow, past the first few steps. Within them (the modes
      beyond `max_angle` still leaving) it is 0.1-0.4 dB. `tests/test_pade_propagator.cpp`
      holds that, agreement with Crank–Nicolson at 40 steps a wavelength (on its own
      second-order grid and Gaussian starter), the exact field to 50 degrees, the elevated
      source's strength, and fourth-order convergence in height (16x from doubling the
      grid, against 4x with `boundary_order` 2).
    - **The kernel pairs terms into two-lane vectors by hand** (`Lanes`, as `Pair` in
      `sphere_database.cpp`), with the term count a template parameter. From scalar source
      the compiler ran one term at a time: 11.7 ns a row against 5.9 (second order, eight
      terms). Each solve stores one
      complex number per row (e_r = a / pivot; the off-diagonals are all a, so the forward
      and back coefficients coincide), which halves what a step streams through the cache.
    - **Not twisted.** Eight interleaved terms are throughput-bound (16 cycles of
      arithmetic a row against 16 of latency), so the two-ended factorization below would
      buy nothing here.
    - **Turbulence caps the step at the screen spacing.**
    - **The 10 Hz benchmark legs start at 50 m, not 20, and case 2's 100 Hz leg at 25.**
      The first is inside a wavelength and a half of the source: the converged march is
      0.29 dB over the exact field at 25 m, and Crank–Nicolson's own error happened to
      bring it inside 0.2. The second runs at ten points a wavelength and closes on the
      series as 1/r, 0.100 dB under at 20 m against a 0.1 margin. No other tolerance
      moved.
  - **`PeScheme::crank_nicolson` is the march this replaced, kept as an independent
    check**: Padé(1,1) at ten steps a wavelength, good to about 30 degrees. It
    differs from the split-step by 0.03-0.12 dB rms on `example_windy.csv` to 1.5 km, and
    against the exact still-air field the error is its own (0.7 dB rms at 25-30 degrees,
    2.1 at 30-35, where the split-step has 0.1 and 0.6). Its left-hand tridiagonal matrix
    is factored once per run. It is latency-bound, and shaped for that (2.7 ns per point
    per step, from 9.0). The Thomas sweep is a recurrence; each row waited on two dependent
    complex multiplies. Three things, all needed:
    - The factorization is twisted (burn at both ends): rows below the middle eliminate
      downward, rows above upward. Each sweep then carries two independent chains.
    - The pivot is folded into the stored coefficients, so each row is one complex
      multiply-add.
    - That multiply-add is written as two explicit `std::fma` per part (`minus_mul`).
      Left to itself the compiler contracts inside the multiply and leaves a subtract at
      the end, which puts 11 cycles of latency on the chain where 8 will do; that alone
      was 5.1 -> 2.7 ns.

    Output matches the plain Thomas sweep to ~1e-9 dB (the same arithmetic reordered), and
    PE runs stored before the split-step reproduce to their printed digits under
    `--scheme cn`. `solve_pe_bands` pools every (band, frequency) march, highest frequency
    first, because cost goes as f^2 and a band at a time left the 2 kHz band on a 2-thread
    share. A five-band `pe_compare` case went from 194 s to 14 s that way. Summing is in
    frequency order afterwards, so any thread count gives the same bits, under either
    scheme.
  - **Crank–Nicolson's grid is 20 points per wavelength in height and 10 in range.**
    (The split-step's is 10 in height; above.) Against a
    30 x 30 reference, 10 in height was ~0.1 dB rms near the ground but 0.47 dB rms (2.2 at
    worst) for a 1.2 m receiver at 1 kHz; 20 is 0.09 dB, for twice the cost. Rerunning the
    stored calibration cases at 20 moved the PE by 0.16-0.39 dB rms (1.1 at worst) and did
    not move the shadow-edge calibration fitted against it, so the fits stand. The
    turbulent mean-field test runs at 10 in height with 1000 realizations; at 20 its
    tolerance was too tight for its sample size. Halving Crank–Nicolson's range steps costs
    most cases only 0.02-0.03 dB, but fails `test_benchmark` by 0.2-0.4 dB at long range. Do
    not coarsen its range step to buy speed; that is what the split-step is for.
  - **The band average takes the midpoints of `samples` equal slices of the third-octave band**
    (`solve_pe_bands`), one sample being the center. It once put the first and last sample on the band
    edges with a full 1/N weight each, which read 0.16 to 0.47 dB rms off the true band mean downwind
    (bias +0.1 to +0.35 dB) and 0.4 to 1.5 upwind, and the stored PE tables and the shadow-edge
    constants below were fitted against that. They stand within their tolerances (every PE-pinned test
    passes), but anything fitted to PE band levels should be re-scored against fresh PE runs
    (`--no_pe` re-scores the ray models against stored ones without re-solving).
  - **Normalization:** Salomons' second-order Gaussian starter times sqrt(i k_a) makes
    p = psi exp(i k_a r)/sqrt(r) a unit point source. `tests/test_pe.cpp` holds it to the
    exact two-ray field (0.3 dB, rigid) and to Weyl–van der Pol through NICE-OPS's own
    `spherical_reflection_coefficient` (0.5 dB, grass), which checks the ground condition's
    sign too.
  - **Results land on exact multiples of `output_step`.** Crank–Nicolson's range step
    follows the wavelength, so it interpolates between steps; without that each frequency
    had its own grid, and band averaging read past the end of one. The split-step lands on
    them (above).
  - **The shadow is checked against an exact solution** (`tests/test_pe_shadow.cpp`). For
    k^2 = k0^2 (1 + 2z/R) over rigid ground the field is a residue series of Airy creeping
    waves, and the PE matches it within 0.02 dB down to -47 dB at 250 Hz to 1 kHz. The complex
    Airy function it needs is `src/airy.{hpp,cpp}`: series for |z| <= 5, asymptotic beyond,
    with the connection formula near the negative axis. Five is where both are ~1e-7; near
    the positive real axis the series cancels exp(4/3 z^1.5) of its digits away. The
    reference values came from Bessel-integral quadrature, not memory, which had Ai(-10)
    wrong. **`airy_asymptotic` and `hankel1_0`'s large-argument branch are in real arithmetic**,
    for the reason `faddeeva` is (the library's complex divide, pow, exp and abs are calls with
    NaN-recovery branches that cannot fire here): `z^(3/2)` and `z^(1/4)` are `z sqrt(z)` and
    `sqrt(sqrt(z))` on the principal branch, the divides are by a reciprocal formed once, and the
    series stop on squared magnitudes. 2-3x per call (H0 209 -> 67 ns, Ai 151 -> 73), and a
    refraction run with a point receiver, whose ShadowTransition table is most of it, went from
    4.9 to 1.6 s. They agree with what they replaced to 1e-16 (H0; 2e-13 where exp(-Im z)
    underflows) and 4e-13 (Ai, where the new one is the closer to a long-double evaluation).
  - **The domain margin is at least 250 m, and the absorbing layer 100 wavelengths.** In a
    shadow the near-ground field is a weak remainder while energy leaks upward past it.
    With 60 wavelengths of margin (25 m at 1 kHz), reflection off the layer set a floor: the
    PE came out 3.8 dB high at -44 dB against the residue series.
  - **Turbulence is phase screens** (`PeSetup::gamma_t`, `outer_scale`, `screen_spacing`,
    `seed`; `pe_compare --gamma_t`). von Karman in the effective index with C^2 = gamma_t / 4
    -- HARMONOISE's gamma_t, the form Ostashev's and HARMONOISE's coherences both reduce to
    -- uniform in height, one screen per outer scale of range, drawn by FFT (`ScreenSource`).
    The screen's vertical spectrum is the 3-D spectrum integrated across the plane, so the
    vertical coherence is 3-D's. `solve_pe_band` makes each frequency its own realization.
    `tests/test_pe_turbulence.cpp` checks the screens' variance (0.06%), the mean field's
    exp(-sigma^2/2) decay (within sampling noise), and ShadowScattering against the PE.
    Still air is untouched: gamma_t 0 skips it all.
  - **`ShadowScattering` (`--scattering`) is fitted to it**, a + 10 lg gamma_t + b lg f +
    c lg R + d lg z_s dB re free field, added in energy to the transition's field where it
    runs, with the floor off. Single scattering is linear in gamma_t and the PE agrees (x3 ->
    4.82 dB, theory 4.77). Two things to know before trusting a strength:
    - **The pair fit's gamma_t = 3e-5 is far too strong for scattering.** It fills the
      day-284 shadows to -3 dB re free field; the measured SEL wants <= ~3e-6
      (HARMONOISE's "moderate"), and does as well with the floor. It is off by default for
      that reason.
    - The 1.2 m constants come from one source height; the height term is the flush fit's.
  - Without `--gamma_t` the deep shadow, 40–90 dB down, is not physical. Absorption is optional
    (`PeSetup::absorption`, Np/m, the imaginary part of k so it acts along the actual path)
    and off in `pe_compare`, because the ray-model comparison leaves it out of all four.
  - **Findings** (`docs/MANUAL.md` has the tables). Downwind, `stratified` is within 0.4–2.1 dB rms of
    the PE and the others are 1–8 dB off. Upwind, the PE's shadow onset lies between
    `stratified` (too early) and `linear` (too late), moving toward `stratified` with
    frequency; the per-band boundary law is that treatment.
  - **`--no_pe`** gives the ray models' columns alone, on the PE's range grid, for
    re-scoring a model change against stored PE runs without re-solving them.
- `src/exact_solutions.{hpp,cpp}`, `tools/benchmark.cpp`, `tests/test_benchmark.cpp` — step 7b:
  the Attenborough et al. benchmark (JASA 97, 173 (1995); the PDF is not in the repo) and its
  exact solutions: direct plus image ray (case 1), the downward normal-mode series (case 2)
  and the upward residue series (case 3), all over impedance ground with absorption as
  Im k0. The PE reproduces the paper's figs. 11–14 in every panel. `docs/MANUAL.md` has the
  tables. Load-bearing:
  - **The mode sums are derived in `exact_solutions.cpp`, not copied, because the paper's
    eq. (11) is wrong.** Its upward ground condition, Ai'(b) + i k0 l beta w Ai(b) = 0, should
    read Ai'(b) - i k0 l beta w² Ai(b) = 0, a factor exp(i pi/3) in the admittance term. The
    printed form misses the PE by 1.5–12 dB at 100 Hz, and the derived one agrees to 0.02 dB.
    Eq. (10) is right up to an overall sign. Over rigid ground the two forms coincide, so
    `test_pe_shadow` could never have caught it. Do not "fix" the code back to the paper.
  - **Case 2's exact profile, c0/sqrt(1 - 2z/R), turns every ray below z = R/2.** Near-vertical
    ones fall back within R = 3.4 km like mortar rounds. The PE cannot carry them (30–45° and
    an absorbing top), and neither can the benchmark's linear profile. Uncut, they put 5 dB
    between the series and the PE by 1 km at 100 Hz, and features at R and 2R that the paper's
    own normal-mode figure shows and its FFP/PE figure lacks. `downward_pressure`'s
    `max_elevation` drops modes steeper than the cut. The test uses 75°: paths reaching 1 km
    leave below 8.4° or above 81.6°. At 10 Hz only ~60 modes lie below 75°, too few for a
    clean cut within 6 wavelengths, so the test uses the full series there. At 1 kHz past
    500 m the PE's own grid is the limit, 0.8 dB at 10 points a wavelength and 0.3 at 20.
  - **`impedance_airy_zeros` is built for |q| <= 6 and throws past it.** It follows each zero
    of Ai' from q = 0, from a first-order seed once |a'| > 4|q|² + 10. It checks the count
    with the argument principle on a box, which caught a zero lost to the continuation in
    development. The benchmark reaches |q| = 4.9. At |q| = 8.5 the continuation loses its
    way, and at 10.8 the box's edge overflows a double. The low zero at 1 kHz is the ground's
    surface wave, near q², and the continuation reaches it. `first_impedance_airy_zero`
    follows the first zero alone, which stays well apart from its neighbors, to any |q|
    (checked to 1e5 against a_1 - 1/q - a_1/(3q^3)): what the shadow loss needs.
    **Past |q| = 6 with Re q <= 0** (every upward problem over ordinary ground) the zeros
    are seeded instead from zeta = (n - 3/4) pi - atan(q / sqrt(x)), b = -x, and polished by
    Newton; checked to step by pi in zeta. That removed the limit for `upward_pressure`.
    Re q > 0 past 6 (downward, and surface-wave grounds) still throws.
  - **A mode sum returns NaN where it cannot be trusted:** terms canceling by more than nine
    orders (the residue series near the source, whose terms grow before decaying), or a tail
    not yet dead at the 12000-mode cap. Callers skip NaN rather than tighten the cap.
  - **`hankel1_0` switches from the power series to the asymptotic expansion at
    |z| = 12 - 3 Im z / |z|.** The series loses exp(|z| + Im z) of its digits, and the
    expansion's smallest term is exp(-2|z|). A switch on |z| alone put H0(3+6i) off by 1e-6.
    Its test values are Mathematica's, as are the root-finder's.
  - **The benchmark's inputs are its own:** the tabulated impedances (Z = 38.79+38.41i,
    12.81+11.62i, 5.96+2.46i), not a NICE-OPS impedance model, and the absorption it quotes
    (0.18 and 5.4 dB/km, ANSI S1.26-1978; ISO 9613-1 gives 0.22 and 4.98). Its profiles stop
    at 1000 m, the height of its contour plots. Only case 2 at 10 km cares, and the paper
    notes the same dependence. A full run is ~10 s (2.5 min under Crank–Nicolson).
- `src/propagation_table.{hpp,cpp}`, `tools/propagation_table.cpp` -- `--propagation_table`: a
  table of each band's level re free field at the straight-line distance, by (bearing, source
  height, receiver height, range, band), for one stratified atmosphere over uniform ground.
  Built by `niceops_propagation_table` from the PE or from the straight-ray ground model; the
  format (netCDF, `propagation_table_version` 1) is the specification, so other methods can
  write it. Things to keep:
  - **It enters through the existing seam.** A covered pair is `spectral_level_along` on the
    straight ray with no ground model and `extra_loss` = minus the table: sphere by the
    straight-line direction, spreading over the straight-line range, absorption as the run
    applies it (zeros if the table `includes_absorption`). Nothing else in the evaluation knows
    the table exists.
  - **Uncovered pairs fall back, per pair, to the run's own model**: `--ground_effect` if on,
    free field if not. A pair is uncovered if any corner it would interpolate is NaN (the PE
    past its angle), or it is outside an axis. Zero-weight corners are not read, so a point
    on the grid beside a NaN is still covered. Microphone heights outside the table are
    refused up front, not dropped to the fallback.
  - **The peak search's bound takes max(ground model's gain, table's largest value) per band**
    (`path_gains` in the bound setup), since a pair may take either.
    `propagation_table_windowed_is_exact` holds it; the PE table was checked the same way by
    hand.
  - **`lookup` finds an even axis's cell by arithmetic and corrects it against the knots**
    (`AxisLocator`), so it is the cell a binary search gives, bit for bit; the fmod on the
    bearing is skipped inside [0, 360). -20% CPU on a table run, output byte-identical.
    `lookup finds the same cell on even and uneven axes` sweeps knots, midpoints and a hair
    either side of them against a plain search. What is still on the table: the source's cell
    is constant over a state's 14,400 observers and could be located once per state (the
    review measured 87 -> 13 ns per call that way); it needs the call site to carry it.
  - **`BandMap` is per run, not in the table**, because the table is loaded once
    (`RunOptions::propagation_table`, shared) and each trajectory can name its own database.
  - **Bearings are compass, toward.** The run turns the local separation back to east/north by
    inverting `CompassToLocal` at `--frame_bearing` (90, +x east, when unset).
  - **Microphone height comes from `--receiver_height` under a table**, as under
    `--ground_effect`; before, a run without either put microphones on the ground.
  - **The plumbing check is the ground-model table.** It reproduces `--ground_effect
    --exact_ground` to the table's interpolation, 0.19 dB at worst at 10 ft spacing, 0.02 at
    5. The mean error falls fourfold per halving, which is what says it is interpolation.
  - **Plates read the table at its lowest receiver height plus K**, the plate over bare
    ground there per band and azimuth at 2 degrees (`table_plate_offset` in niceops.cpp), as
    ShadowTransition places a plate. **With `--refraction` the table takes the pairs it
    covers and the ray model the rest** (`table_level` first in `refracted_level`); that
    run evaluates every state, since its bounds were not built for both geometries.
  - **Not done:** ground elevation (flat ground at z = 0), and a sphere lookup along a
    refracted launch direction (the table cannot say what it is).
- **Receivers** (`process_trajectory`). Every level is evaluated at the microphone: the
  observer's point raised by its height above the ground. A flush microphone (height 0) with
  `--plate_table` is on a plate; any other, with `--ground_effect`, is a point receiver at its
  height, whose ground model and `GroundTable` are built once per distinct height. An observer
  file's fifth column (`id,x,y,z,height`) gives heights per observer, so one run mixes plates
  and poles, as test arrays do. Without it everyone takes the run's receiver. The peak
  search's ground bound is the larger of the plate's and the point receiver's, which does not
  depend on height. `mixed_receivers_are_their_subsets` checks a mixed run against runs over
  each kind (with `--exact_ground`, because a table's extent follows its run's observers).
  The point receiver's direct path ended on the ground until this went in, while its ground
  model put the microphone 1.2 m up. Point-receiver SEL moved by up to ±0.8 dB on the demo
  grid, almost all from the low states; plates and the goldens did not move.
- `src/weighting.{hpp,cpp}` — A-weighting (panam's `dBAw`, the IEC 61672 pole-zero form rather
  than a table, so it evaluates at whatever band centers a database carries) and band-limited
  energy integration. Note the curve is **not** identically 0 dB at 1 kHz in this form but
  -1.49e-6 dB; a test asserting exact 0 there will fail. `a_weighting_gains()` folds the
  weighting and the band limit into one per-band multiply, with excluded bands exactly zero.
  `weighted_level()` returns -infinity for silence, not 0 dB, so an excluded spectrum cannot
  survive an energy sum as a quiet-but-real level.
- `src/faddeeva.{hpp,cpp}` — the Faddeeva function w(z) = exp(-z²)erfc(-iz), Weideman's 32-term
  rational approximation. The series is evaluated as **four interleaved Horner chains over
  t⁴** rather than one over t: 31 dependent complex multiply-adds is a latency chain, not a
  throughput problem, and it was running at about a quarter of the machine's floating-point
  rate. Splitting by the power of t modulo four gives chains a quarter as deep that issue
  together — 32% off `faddeeva`, 14% off the whole ground path. The coefficient table is
  left exactly as Weideman's construction produced it and indexed by stride, so it stays
  checkable against the source. It exists solely for the ground model's boundary-loss term, which panam
  writes as `exp(-w)(1 - erf(-i sqrt w))` using SciPy's complex `erf`; with z = sqrt(w) that *is*
  w(z), and evaluating it as written overflows. Max relative error 3e-13 over the domain the
  ground model reached while the boundary loss was cut at |w| = 500 (Re z ∈ [0.09, 22.3],
  Im z ∈ [-1.2, 15.6]), measured against SciPy's `wofz`; past it the tests hold it to `wofz` at
  1e-11 out to |z| = 700. The lower-half-plane reflection `w(z) = 2exp(-z²) - w(-z)` is safe
  because a passive ground keeps z within 45 degrees of the real axis wherever it is below it,
  so exp(-z²) is at most 1 in size.
- `src/ground_effect.{hpp,cpp}` — excess ground attenuation, ported from panam's `ega()` and pinned
  against it to 1e-9 dB. **It returns a gain, not a loss** — +6 dB over rigid ground with the
  receiver on the surface. `ground_impedance()` is hoisted out of the observer loop and the
  allocation-free overload is what the loop calls; returning a vector there allocated 94 million
  times and cost 20% of the run.
- **The spectral and ground paths were profiled, not guessed at.** Two rounds of intuition
  (cache layout, then `pow`) bought 14% between them; `sample` on a running process then showed
  the real costs and the same effort bought 2.5x. What it found, and what not to undo:
  - `std::sqrt` on a `std::complex` costs an **atan2**, because libc++ evaluates it in polar
    form. `ground_effect.cpp` called it three times per band; two were provably redundant
    (`sqrt(w)` *is* the numerical-distance parameter wherever its real part is positive, which is
    every passive ground but a mass-like one; see below) and the
    third is an algebraic root written out in the band loop. `std::abs`/`std::arg` on a complex cost a
    hypot and an atan2 and are gone the same way -- `|Q|cos(p + arg Q)` is
    `Re(Q)cos p - Im(Q)sin p`, and `|w| <= 500` is `norm(parameter) <= 500`.
  - `faddeeva` is written in real arithmetic. Every `std::complex` multiply carries a branch to
    recover NaN/infinity cases that cannot arise here, and the Horner loop does 31 of them.
    Its two trailing divides are folded into one: `2p/d^2 + c/d == (2p + c d)/d^2`.
  - `spectra_` is condition-innermost, matching `grid`. The netCDF file stores it the other way
    and the loader transposes once.
  - The whole band loop in `excess_ground_attenuation` is real arithmetic for the same reason,
    and its band-invariant branches (boundary loss on/off, the sinc, turbulence) are hoisted so
    the body stays straight-line.
  - Bands outside the metric's limits are not computed at all. `first_band`/`last_band` bound a
    contiguous run, which is always what the limits select because band centers are sorted.
  - The Weideman series behind `faddeeva` runs as four interleaved Horner chains rather than
    one; see `src/faddeeva.cpp`. It was 53% of the ground path on its own and is latency-
    bound, not throughput-bound.
  All of it is exact: the footprint moved by 2.8e-14 dB end to end, and the EGA tests still pin
  to panam at 1e-9 dB. Acoustic loop on the Be407 case, default band limits: broadband 0.34 s,
  `--spectral` 1.44 s, `--spectral --ground_effect` 10.9 s (was 41.4); the Horner split takes
  the last of those to about 9.4 s.

  The `*_fixed_load_factor.nod` databases carry spectra, so the path runs end to end on shipped
  data; it was also exercisable from a bare checkout when none did, and still is: `write_flat_spectral_database` in `tests/test_spectral_propagation.cpp`
  interpolates a sphere out of the shipped `S-76D_M3.nod`, hangs synthetic band amplitudes on
  it and writes it back with `write_sphere`, and the result is a one-group `.nod` that
  `--spectral` accepts. Its single condition makes the blend degenerate, so timings off it run
  optimistic, but the ground term -- which is most of the cost -- does not care. On the demo
  track that way: `--spectral` 1.36 s -> 1.16 s, `--spectral --ground_effect` 8.95 s -> 7.81 s.
  For the ground term alone, driving `excess_ground_attenuation` over a spread of source
  heights and ground distances with 26 synthetic third-octave bands reproduces the per-call
  cost the observer loop sees (~840 ns, was ~975), and `tests/test_ground_effect.cpp` -- pinned
  to panam at 1e-9 dB -- is what says it is still right.

  **The asymptotic boundary loss was tried and rejected.** For large `|w|` the term tends to
  `-(1/2w + 3/4w^2 + ...)`, which would skip `faddeeva` on ~20% of bands. It does not converge
  to the exact value -- the error plateaus at 1.6e-5 at `|w| >= 50` and 1.9e-6 at `|w| >= 100`
  no matter how many terms are taken -- so it cannot hold the 1e-9 dB pinning against panam at
  any cutoff in range. Roughly 9% of the ground path, for the loss of the property that proves
  the port. Not worth it; do not re-derive it.

  **There is no cutoff any more.** panam cut the boundary-loss term to zero past |w| = 500, and
  this followed, until panam e978703 evaluated it by `wofz` at any |w|; NICE-OPS now does the
  same, pinned to panam b758e3a (`ega` 1e-9 dB, `pole_level` and `board_level` 1e-6). Q moves by up
  to 2e-3 just past the old cut; the demo footprint with the plate receiver by at most 0.028 dB
  (0.08 dB with a 4 ft point receiver, where the two paths nearly cancel), means ~1e-4. Three
  things to keep:
  - **The gain bound past the edge is analytic** (`reflection_coefficient_bounds`). Past
    |w| = 500, |Q| <= |R_p| + |1 - R_p| (1.003/(2|w|) + |P|), P = 2i sqrt(pi) z exp(-z²) the
    surface wave, whose size peaks at |z|² = 1/(2 cos 2 theta), theta the parameter's direction,
    which does not depend on range. For soils that peak is far inside the edge. For a nearly
    lossless stiff surface (a thin layer's band, theta a hair under -45 degrees) the wave
    carries on: |Q| = 151 at |w| = 1250 at grazing, where the cut gave -1. The sweep stops at the
    edge and the formula takes over; `a surface wave carried past |w| = 500 ...` holds it to a
    sweep out to kappa = 1e10.
  - **A mass-like admittance (Im beta > 0) keeps NICE-OPS's branch, not panam's.** Its
    parameter is in the second quadrant; panam squares it and takes the principal root, the
    negative, whose w() grows as exp(|Re w|) and overflows by |w| ~ 1000. Such a surface has no
    surface wave, and the parameter itself gives none. A thin layer gives this in some bands.
  - The near-grazing three-term series that stood in past the cut (a refracted ray arriving
    level at a flush receiver, R_p -> -1) is gone: the exact term is what it approximated.
- **The 4x EAA multiplier is gone, and the goldens moved with it.** `eval` used to compute
  `4.0 * EAA * (range - r)/1000` under a `//TODO Add ground effects`, inherited verbatim from
  the D original (`c74bab3`). panam *writes* that EAA as `SPLA - SPLAa` over a 1000 m reference
  and consumes it as `EAA / distance * (slant_range - radius)`, with no such factor, so the
  multiplier overstated atmospheric absorption fourfold — 168 dB at 6 km where the physical
  value is 14. It was not a ground effect either: measured against the model it stood in for,
  soft ground costs a roughly range-independent 2–4 dB, where the factor removed 7 dB at 250 m
  rising to 38 dB at 2 km. Removing it raised the demo footprint by **9.93 dB mean, up to
  22.7 dB**. `golden/197421_LAE.csv` was regenerated; `golden/state_197421.csv` was not, being
  trajectory state.
- **Two propagation paths, and they do not agree.** The default carries one broadband dBA and
  attenuates it with the database's `EAA` factor; `--spectral` carries the per-band spectrum
  and attenuates each band at its own ISO 9613-1 coefficient. `Database::eval_spectrum` mirrors
  `eval`'s bracketing and bilinear cell exactly, so a difference between them is the physics and
  never the interpolation. `--spectral` is the default whenever the database carries spectra (and
  falls back to broadband when it does not), so the goldens -- made from `S-76D_M3.nod`, which has
  none -- are the broadband path. On B407, spectral minus `--no-spectral` is median -1.2 dB
  (-3.2 to +2.3) on the demo track: the humidity and temperature of the run against the ones the
  database's EAA was built at, and EAA being a linear-in-range rate anchored at 1 km where the
  true loss is concave. `spectra_` is condition-major (right for `interpolate_sphere`, wrong for this loop), so the
  spectral path pays a long stride per band rather than making the sphere export pay.
  `--db_cutoff` works on it through `spectral_direction_bounds`, the band-limited level's own bound.
- `src/geodetic.{hpp,cpp}` — the WGS84 transform: the American polyconic projection about a
  reference point (it is not a transverse Mercator, whatever older copies of this file say; it
  reproduces the shipped tracks' own x/y, and agrees with a Krueger transverse Mercator to 3 mm
  at 10 km), with altitude carried through as z. `LocalFrame(reference)` is east/north;
  `track_aligned(reference, bearing)` turns +x onto a compass bearing;
  `from_east_north_up` applies that rotation without the translation, for velocities; and
  `to_geodetic` inverts the projection by Newton iteration, which is what `--latlon` output
  uses. Ported from the D `worldtransforms.d`/`wgs84tocartesian.d` pair, which had no notion of
  a rotated frame and whose inverse walked fixed 1e-5 degree increments.
- `src/smoothing_spline.{hpp,cpp}` — natural cubic smoothing spline (Reinsch's algorithm; interior second derivatives solve a symmetric pentadiagonal system in O(n)). Two ways to pick how hard to smooth: `with_smoothing(lambda)` takes the penalty weight directly, `within_tolerance(tol)` solves for the lambda whose weighted squared residual just meets `tol`, as MATLAB's `spaps` does -- and returns the weighted least-squares line in closed form when `tol` reaches that line's residual (the search it replaced climbed lambda to 1e200 and, from ~30,000 samples, lost positive definiteness on the way and refused a straight, level leg; the system's condition number goes as n^4, so past ~40,000 points a `tol` just under the line's residual is unreachable and the last solvable lambda is used, with a warning). `within_tolerance_extended` does the same over a track padded along a local quadratic at each end, which is what a dense geodetic track uses (see "Velocity and load factor are derived" above). `trajectory.cpp` uses it directly. Replaced the vendored SPLINTER library, which brought a bundled Eigen 3.2.8 that predates C++17.

Both `extern "C"` seams above are leftovers of the old D FFI boundary, not a current requirement.
### Geodetic trajectory input

A trajectory may be given in WGS84 (latitude, longitude, altitude) instead of the 27-column AC
state format; `docs/MANUAL.md` documents the format and the options. Things worth knowing before
touching it:

- **The format is decided by width, not by the header.** 27 columns is the positional AC
  format; anything else is read by the names in its header (`read_labeled_csv` in
  `src/util/csv.hpp`). That order is load-bearing: the AC format's *own* header names `lat`,
  `lon` and `alt` columns, so a reader dispatching on the header alone takes a state file for a
  WGS84 track — and its x/y origin, which the file does not record, is about 2.5 km back along
  the track. `first_geodetic_position` applies the same width test, and a test pins it.
- **The AC path is byte-identical.** `read_labeled_csv` keeps the first record as the header
  instead of discarding it, but drops exactly the rows `read_numeric_csv(path, true)` dropped,
  so the goldens are untouched. Verified.
- **`alt_m` is deliberately not a recognized altitude spelling.** Altitude is read in feet;
  matching a meter spelling would read meters as feet and fly the aircraft three times too
  high, where leaving it unmatched fails the file with a message naming the units.
- **Altitude is MSL and z is height above the observers.** `--ground_elevation` is the
  difference. In `tracks/S-76D/197421AC.csv` it is a constant 54.2651 ft, which is how that
  file's own alt and z columns relate.
- **The shipped AC tracks' local frame runs +x along a compass bearing of 210°** — the approach
  heading, not the instantaneous heading in column 8, which wanders over 25°. With that bearing
  the projection reproduces that file's own x/y to under a meter over 3.2 km, which is the
  check that the port of the transform is the one its producer used.
- **One frame per run, settled before anything loads.** Every trajectory in a run shares the
  observers, so a position in one has to mean the same place as in another. `run()` resolves it
  from `--origin` or from the first geodetic track and hands it to each `Trajectory`.
- **Velocity and load factor are derived, and that is the whole accuracy story.** The
  projection is not the lossy part: altitude passes through exactly and x/y to under a meter.
  Differentiating a spline twice for the load factor is. **`--path_smoothing` defaults to 10 ft**
  and is an acceptable deviation in feet RMS per axis, not an abstract weight -- the solver
  spends the whole budget rather than merely staying inside it. Interpolating exactly
  (`--path_smoothing 0`) amplifies the coordinates' 1e-8 degree rounding by 1/step² -- 2500x on
  the 0.02 s demo track -- and moves the footprint by **8.0 dB at the median**; the 10 ft
  default takes that to 0.073 dB. Zero is wrong for both kinds of input, which is why it is no
  longer the default: measured positions turn into rounding noise, and sparse waypoints are a
  statement of intent, so flying exactly through them manufactures maneuvering that is itself a
  noise source. 10 ft is not the optimum for any one input shape (a dense 50 Hz track does
  better at 0.0656 ft, tuned to its own resolution) but it is the only value that is good on
  all of them at once.

  Those numbers were re-measured after `SmoothingSpline` was fixed to give the lambda = 0 case
  its curvature. Before that, `--path_smoothing 0` returned a spline with zero second
  derivative everywhere, so the load factor was exactly 1 g on every state, CT never moved, and
  the amplification above never happened -- the old "17% and 3 dB" measured that bug. The manual
  carries the re-measured table for all four input shapes; do not re-derive it.

  **A dense track's ends are freed from the natural condition** (`fit_path`,
  `SmoothingSpline::within_tolerance_extended`). A natural spline's second derivative is zero at
  the first and last knot and still a third low four seconds in, so the load factor, and with it
  CT, is wrong exactly where a track that starts or stops mid-manoeuvre is loudest. For positions
  at 0.5 s spacing or finer with a nonzero `--path_smoothing`, the fit is made over the track
  extended by 10 s each way along a local quadratic through the last 3 s, so the natural
  condition lands on the extension, and only the part over the data is kept. Dense demo
  footprint against the AC state: median 0.109 -> 0.073 dB, max 1.84 -> 1.42, cells over 1 dB
  84 -> 25. Sparse waypoints and `--path_smoothing 0` are unchanged, bit for bit, on purpose: a
  quadratic through a few waypoints is a guess, and the exact interpolant has no bias to
  remove. The AC path does no smoothing, so the goldens did not move.
- **Resampling never coarsens.** The sound exposure integral is a Riemann sum over whatever
  step the states end up with, so a 50 Hz file keeps its own sample times and only sparse or
  unevenly spaced input is refitted onto the `-s` grid.

- `src/trajectory.{hpp,cpp}` — `Trajectory` class. Parses the 27-column CSV into `State` records (metric), computing non-dimensional flight parameters: CT = LF·W/(ρ·V_tip²·πR²), μ = v/V_tip, and flight path angles `fpa` (vertical) vs `fpa_e` (relative to the TPP normal). The `qsam` CLI flag selects `fpa_e` (default on). `adapt()` resamples onto a uniform dt grid, smoothing via `src/smoothing_spline.hpp` and taking velocity and load factor from the spline's own derivatives (currently commented out in `main`; `watts_trim()` is reachable through `--trim`). It takes either a smoothing parameter or a `tolerance` bounding how far the smoothed path may sit from the measured one.
- **At hover the sphere's nose is the heading** (`nose_direction`, `heading_offset` in
  `src/trajectory.cpp`). Below 3 kt (`heading_nose_speed`) the velocity is drift, a few tenths
  of a knot in any direction, and the nose went wherever it pointed. Between 3 and 8 kt
  (`velocity_nose_speed`) the nose turns from the heading to the velocity on a smoothstep weight,
  the horizontal angle the short way round, so there is no step; the hard switch at 5 kt it
  replaced stepped the level 2.2-4.4 dB for a nose rotated 30-90 degrees. The AC format's column 8
  (`State::heading`, compass, NaN for the geodetic formats) is turned into the local frame by
  an offset inferred from the track's states over 20 kt -- the circular mean of the velocity's
  angle plus the heading, which folds in the INS heading's few degrees off the ground track --
  or, for a track with none (a hover), by `--frame_bearing`. Neither: the drift, with a warning.
  No shipped track goes under 10 kt and no day-286 approach under 13, so neither the goldens
  nor the approaches move; the 2017 hover runs (H1-H12) do.
- **`--speed_reference`** (roadmap step 3, the NICE-OPS half). It takes the flight condition
  against the air, not the ground. The things to preserve:
  - **Scope.** `Trajectory::relative_to_air` sets each state's wind (`State::wx`, `wy`) and
    recomputes μ, `fpa` and `fpa_e` from `air_velocity()`. `process_trajectory` hands the same
    air velocity to the sphere frame, so the nose follows it too. LF, CT and the load-factor
    normal stay inertial.
  - **Calm is bit-identical.** `vx - 0.0` is `vx`, and `speed_reference_air_calm_is_ground`
    checks it end to end.
  - **The database decides by default.** The root text attribute `speed_reference`, `ground`
    or `air` (absent is ground), is read into `Database::speed_reference()`. It is settled per
    trajectory in `refer_to_database_speed`, because each trajectory can name its own
    database. Air with no wind in `--atmosphere` is refused, not taken as calm.
  - **It runs after `--trim`,** whose "slowest point" is over the ground.
  - **One rotation.** `CompassToLocal` (atmosphere_profile.hpp) turns east/north into the
    local frame for this and for both ray models. Before it, refraction.cpp and
    stratified.cpp each spelled the rotation out, with the same expression.
  - **panam's side is not built yet.** It has to use the same definition, the manual's
    "Speed reference" section, and write the attribute. Until it does, every database is
    ground-referenced, and an air run against one mixes the two.
- `src/plotting.{hpp,cpp}` — `plt::` facade replacing the vendored matplotlib-d binding. Embeds Python via pybind11 and calls matplotlib directly; the interpreter is created lazily on first use, so non-plotting runs never pay the import cost.

### D standard-library replacements

The port swapped each D dependency for a small local equivalent; when a behavior looks odd, check the D-side original named in the header comment.

| D original | C++ replacement |
| --- | --- |
| `std.getopt` | CLI11 |
| `std.json` | nlohmann_json |
| `std.csv` `csvReader!(double)` | `src/util/csv.hpp` (`read_numeric_csv`, plus `read_labeled_csv` for the header-named formats, which have no D original) |
| `std.experimental.logger` | `src/util/logger.hpp` (`log`/`info`/`trace`/`fatal`) |
| `std.numeric.findRoot` | `src/util/rootfind.hpp` (bisection) |
| `std.math.approxEqual` | `approx_equal()` in `sphere_database.hpp` |
| `unittest` blocks | Catch2 cases in `tests/` |
| `scope(exit) nc_close` | `NcGuard` RAII in `sphere_database.cpp` |
| matplotlib-d | pybind11-embedded matplotlib |
| SPLINTER (vendored B-spline) | `src/smoothing_spline.{hpp,cpp}` |
| `netcdf.d` | the netCDF C API directly |

## Licensing

MIT (`LICENSE`), cherry-picked from master, where it had been added after this branch forked.
`THIRD_PARTY_NOTICES.md` carries the verbatim notices for netCDF, CLI11, nlohmann/json,
pybind11 and Catch2 -- all BSD/MIT/BSL, none copyleft. Regenerate it from the fetched sources
under `build/_deps` if a dependency is added, changed or dropped; the plotting packages in
`build/pyenv` are invoked rather than linked and are covered by a summary paragraph rather than
reproduced.

The `Agg` backend choice in `plotting.cpp` is load-bearing for this: matplotlib's Qt backends
would pull in PyQt, which is GPL or commercial.

## Code Conventions

### Porting pitfalls already hit (do not reintroduce)

The port was validated against the D reference and these four bugs were found and fixed. Each was a plausible-looking mistranslation:

- **`iota` is EXCLUSIVE of its end value.** D's `iota(-1500.0, 1500.0, 25.0)` yields **120** values ending at **1475**, not 121 ending at 1500 (verified by compiling D). This governs the observer grid and the contour levels. Comments in the port previously asserted the opposite.
- **Don't sort `phis_scattered`/`thetas_scattered` in place.** D builds the unique axis values from a *copy* (`.dup.sort().uniq().array()`). Sorting in place collapses the per-channel arrays that the grid fill pairs index-wise with `dBAs_scattered`, leaving most cells zero. Symptoms: `dba_max` too low, database load suspiciously fast (~0.3 s instead of ~2.9 s), `--self_test` returning exactly 0, levels tens of dB low.
- **D default-initializes `double` to NaN, not 0.** `ambient_level` relies on this: the `isNaN` guard means D adds no ambient unless `-a` is given. Initializing it to `0.0` adds a spurious 0 dB floor. (The mix itself is no longer D's: `-a` adds the ambient in
  *energy* (`level_with_ambient`, `sel.hpp`), where D summed pressure amplitudes, which added two
  uncorrelated levels coherently -- +6.02 dB over two equal levels, up to 3 dB high -- and was never
  tested. A departure from D by the port-fidelity bar: measured against the physics, not a reading.)
- **`approx_equal` mirrors D's `approxEqual`**: relative term against the *second* argument, defaults 1e-2 relative / 1e-5 absolute. It is used for database channel matching, so tightening it silently drops channels.

General rule: when a D construct looks obvious, check it against Phobos semantics rather than intuition — and prefer verifying with a throwaway `ldc2` program, which is how the `iota` question was settled.


- **Port fidelity outranks correctness, with one deliberate exception.** Several `NOTE:`
  comments mark places that faithfully reproduce D quirks; leave them alone unless you are also
  regenerating the golden outputs, and say so explicitly. The exception is the 4x EAA multiplier,
  removed above: it was a physics error rather than a quirk, it was measured against the model it
  stood in for before being removed, and `golden/197421_LAE.csv` was regenerated from the C++.
  That is the bar for departing from D — a measured case that the D behavior is wrong, not a
  reading that it looks wrong.
  - Contour defaults are 65/100, not c74bab3's 45/90: the older pair puts its top level below the demo case's 93.6 dB peak, leaving it unfilled, and the shipped reference figures use 65/100. Figures only -- the goldens are CSV.
  - `watts_trim` and the observer grid extent are now CLI options (`--trim`, `--grid`), both defaulting to the port's behavior so the goldens are untouched. With `--trim -6000 --grid -2000 1100 -1500 1600` the demo case reproduces the June 2020 reference figures to within 1.05 dB on every cell.
  - D's `iota` is EXCLUSIVE of the stop value; observer grids and contour levels reconstruct that off-by-one deliberately (`iota_exclusive` in `src/niceops.cpp`).
  - `write_levels` emits one extra trailing newline because the D original did.
  - `level_along_bearing` uses `heading[0]` twice where the second component wants `heading[1]`; that is the D original's latent bug, kept.
  - `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does not run on fatal paths.
- `using namespace niceops;` at the top of the `src/*.cpp` files is intentional (marked `// NOLINT`) so ported code reads like the D original.
- Logging goes through `src/util/logger.hpp`; `main` points it at `log/niceops.log`. The logger's `print` overloads handle vectors and arrays inline, as D's did.
