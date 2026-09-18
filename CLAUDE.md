# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NICE-OPS (Noise Informed Community Environment Operation Planning System) is a helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer grid (or scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

The implementation is C++17, under `src/`. It began as a port of a D original, which is why several comments reference D semantics — those are load-bearing: the port reproduces the D reference bit-for-bit on the demo case, and the goldens come from it.

- `src/` — the implementation. **Every `.cpp` lives here and is listed in CMakeLists.txt**; there are no C++ sources at the repo root.
- `tests/` — Catch2 unit tests plus the `golden_compare` tool.

The D sources (`niceops.d`, `sphere_database.d`, `trajectory.d`, `splines.d`, `netcdf.d`), along with `cmake-d/` and the vendored `external/matplotlibd/`, were removed once the port was verified complete: every D function, all 19 CLI options and all 10 `unittest` blocks have C++ counterparts. **They remain in history at `c74bab3`** (on `origin/no-stk`), which is what to check out if you ever need the oracle again — see the golden-file section.

When you add a C++ source, put it in `src/` and add it to the right target — a file that compiles but is never linked in shows up as an undefined-symbol error far from its cause.

## Build

Requirements: a C++17 compiler, netCDF (shared library), CGAL, Eigen3, Boost, and Python 3 **with development headers** plus matplotlib (pybind11 embeds the interpreter for plotting). On Ubuntu, as CI does:

```sh
sudo apt-get install libcgal-dev libnetcdf-dev libeigen3-dev python3-dev python3-matplotlib
```

CLI11, nlohmann_json, Catch2, and pybind11 are pulled at configure time via `FetchContent`, so the first configure needs network access.

Build out-of-source, and **always with explicit `-S`/`-B`**:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

> **Use `-S`/`-B`, never `cmake <srcdir>` from a build directory.** The repo root still holds a stale `CMakeCache.txt` from the D era. When the directory you hand `cmake` contains a cache, CMake treats that argument as the *build* tree — so `cd build && cmake ..` silently configures **into the source root**, reusing the 2023 cache. That cache carries paths from an old Homebrew prefix (`/usr/local/...`) and produces baffling, unrelated errors such as `Imported target "CGAL" includes non-existent path "/usr/local/include"`. Deleting the root `CMakeCache.txt`/`CMakeFiles/` disarms the trap permanently.

Notes:
- The repo root contains stale in-source build artifacts from the D era (CMakeCache.txt, Makefile, `libniceops*.a`, binaries) that are gitignored — do not build in-source.
- **Never add `-DNDEBUG` to Release.** The D original ran with asserts live, and `sphere_database.cpp` relies on `assert(isFinite(...))` guards; matching that is part of the port contract. This is called out in CMakeLists.txt.
- Vendored splinter bundles Eigen **3.2.8** (2015), which does not compile under C++17 (`std::binder1st`/`binder2nd` were removed). Point it at the system Eigen instead: `-DEIGEN_DIRECTORY=/path/to/eigen3` (e.g. `/opt/homebrew/include/eigen3`). This also avoids linking two different Eigen versions into one binary, since `niceops_core` uses `Eigen3::Eigen`.
- Splinter's own test target (`splinter-test`) still won't build under C++17 — its bundled Catch uses `std::random_shuffle`. It isn't needed, so build just what you want: `cmake --build build --target niceops niceops_unit_tests golden_compare`.
- If Xcode proper isn't installed (Command Line Tools only), pass `-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`. No `-DCMAKE_POLICY_VERSION_MINIMUM` is needed: vendored splinter's `cmake_minimum_required` was bumped from 2.8 to `3.5...3.27`, which CMake 4 accepts.
- CGAL 5.0+ is header-only and exports `CGAL::CGAL`. Do **not** reintroduce `include(${CGAL_USE_FILE})` — `UseCGAL.cmake` was removed in CGAL 6, leaving `CGAL_USE_FILE` empty so the `include()` errors out.

### Golden reference files

`golden/197421_LAE.csv` and `golden/state_197421.csv` are captured from the **D reference implementation** and the C++ port matches them to **1e-9** (actual max relative difference ~1e-14, i.e. a few ULPs of LDC-vs-clang round-off). `regression_golden` and `state_golden` both pass.

To regenerate them, build the D original in a scratch worktree — do NOT do this in the working tree:

```sh
# c74bab3 is the last commit that still contains the D sources and the D-capable
# CMakeLists; 9ccbf0a is the D-era build (c74bab3's CMakeLists is C++-only).
git worktree add --detach /tmp/dref 9ccbf0a
git show c74bab3:niceops.d > /tmp/dref/niceops.d   # see below: 9ccbf0a's copy is NOT the oracle
cmake -S /tmp/dref -B /tmp/drefbuild -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_D_COMPILER=$(which ldc2) -DNCDF=/opt/homebrew/lib/libnetcdf.dylib \
      -DCMAKE_CXX_FLAGS=-I/opt/homebrew/include/eigen3 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build /tmp/drefbuild --target niceops
```

Three things will bite you:

- **Take `niceops.d` from `c74bab3`, not from the D-era commit.** The port was derived from a later revision of it: contour defaults 45/90 (the D-era copy has 65/100), and `watts_trim()`/`adapt()` both commented out. The older copy runs trimming and spline adaptation, computes a different trajectory, and yields a golden the port can never match.
- **D prints only 6 significant figures** by default, so a stock D build agrees with the port to just 5e-6. For a 1e-9 golden, patch the writers to `%.17g` (`write_levels` in `niceops.d`, `write_state` in `trajectory.d`, both recovered from `c74bab3`) — formatting only, no numerics.
- The D CMakeLists hardcodes `/usr/local` paths for netCDF and Eigen; override as above on Apple Silicon. D also crashes if `state/` does not already exist (it does not create the directory, unlike the port).

## Tests

From the build directory:

```sh
ctest -V -C Release -E broken      # what CI runs
ctest -R unit_test -V              # one CTest test
```

The unit tests are a Catch2 binary, so you can run a single case or tag directly — far faster than going through CTest:

```sh
./niceops_unit_tests "[sphere_database]"      # by tag
./niceops_unit_tests "lookup_angles"          # by test name
./niceops_unit_tests --list-tests
```

CTest defines five tests, all with the **source root** as their working directory (the demo case's relative paths resolve there):

- `regression_default` — runs `niceops -S state/ci_state.csv`; one run feeds both golden comparisons.
- `self_test` — runs `niceops --self_test`, validating every sphere in the database against its own stored levels (1,844,672 samples). **Passes.** It was `broken_self_test`, excluded from CI, while legacy hover directivity made it fail on 7,788 samples; that is fixed (see below), so it is now a real test and CI runs the full suite with no `-E` exclusion.

- `unit_test` — Catch2 cases in `tests/`, ported from the D `unittest` blocks. All 7 cases pass (~26k assertions). Note Catch2 rejects `&&` and chained comparisons inside `REQUIRE` unless the whole expression is wrapped in a second pair of parens: `REQUIRE((a && b));`.
- `regression_golden`, `state_golden` — `golden_compare` diffs the run's CSV output against a golden file elementwise at 1e-9 relative tolerance, skipping non-numeric header lines. See the caveat above: the golden files are missing.

### Legacy databases (`database_version` < 1)

Databases carry a root-level `database_version`. Absent or 0 means the file predates the versioned builder, and its **hover spheres have broken directivity**: the upper surface was completed with `-phi` instead of `180 - phi`. Since the source spheres span only phi in [-90, 90] (phi = 0 straight down), negating is a left/right flip *within* the lower half that maps the range onto itself — so the mirrored channels wear lower-hemisphere labels and the upper surface is missing. Forward-flight spheres in those files were completed correctly; only hover is affected (in `S-76D_M3.nod`, exactly the 24 spheres at mu = 0).

Only the labels are wrong, so `Database::Database` repairs them on load via `repair_legacy_hover_directivity()`: for the mirrored half, `phi_correct = wrap(180 + phi_stored)`. It refuses to touch anything that does not match the broken shape (must be confined to [-90, 90], even channel count, halves actually mirrored) and logs how many spheres it fixed. `self_test` applies the same repair, or it would compare a grid built from repaired labels against queries at the raw ones.

Two consequences worth knowing:

- Channel angles are read **per sphere**, not once from sphere 0. Reusing sphere 0's array happened to give the right answer for `S-76D_M3.nod` — the repaired hover labels are byte-identical to sphere 0's — but only because every forward sphere there shares one layout. It would silently misfile data in any database with heterogeneous layouts.
- The repair leaves the demo footprint bit-identical, for the same reason, so `golden/` did not need regenerating.

Regenerating the database with the current builder (panam, `DATABASE_FORMAT_VERSION = 1`) makes the fix-up unnecessary; the loader then trusts the file as written.

## Running the Tool

Default demo case (no arguments): trajectory `tracks/S-76D/197421AC.csv`, database `databases/S-76D_M3.nod`, vehicle config `databases/S-76D.json`, output `integrated/197421_LAE.csv`. Key options (see the CLI11 setup in `src/niceops.cpp`): `-d` database, `-t` trajectory CSV, `-v` vehicle JSON, `-o` output, `-b` observer file (switches to scattered-observer mode), `-a` ambient level, `-i` time-history prefix, `-m` metrics file, `-S` state file, `-p`/`-f` plotting, `-q`/`--no-qsam`, `--self_test`.

Input/output files use **feet**; internal computation is metric (`ft2m = 0.3048`). Trajectory CSVs are the 27-column NASA/FAA AC state format. Vehicle JSON keys: `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg, converted to Newtons), `Atmosphere.density`.

## Architecture

Pipeline in `main()` (`src/niceops.cpp`): load `Database` + `Trajectory` → for each observer and each trajectory state, compute instantaneous dBA via `Database::level_at_observer` → integrate the time history into a sound exposure level (SEL) using the 10-dB-down-from-peak window (`peak_range`/`sound_exposure_level` in `src/util/sel.hpp`) → write CSV output and optional matplotlib figures.

The `niceops_core` static library (`sphere_database.cpp`, `trajectory.cpp`, `plotting.cpp` + the `interpolator` library) holds everything; `niceops.cpp` is only the CLI driver, and the tests link `niceops_core` directly.

- `src/sphere_database.{hpp,cpp}` — `Database` class. Loads a netCDF noise database (`.nod`): one group per operating condition (μ advance ratio, γ flight path angle, CT thrust coefficient), each holding scattered (phi, θ) channels of dBA level and EAA excess-attenuation factor. `eval()` computes ground level: geometry → sphere angles (`lookup_angles`), barycentric interpolation between operating conditions, periodic bilinear interpolation on the sphere grid, then spherical spreading + excess attenuation (EAA × range) + a fixed 4 dB ground effect. Only `same_grid == true` databases are supported (the interpolation path is unimplemented). Free functions (`bracket`, `periodic_bilinear_interp`, `lookup_angles`, the geometry helpers) are exported in the header so the unit tests can reach them.
- `src/interpolator.cpp` — C++/CGAL 3D Delaunay triangulation over the operating-condition space; exposes `set_interpolator3`/`barycentric`/`bilinear_interp` via `extern "C"` (declared in `sphere_database.hpp`). Built as its own `interpolator` library so CGAL's heavy headers compile once.
- `src/splinter.cpp` — `extern "C"` wrapper over the vendored SPLINTER B-spline library, declared by `src/splines.hpp`, which wraps it back up in an RAII `Spline` class. `splinter-static` supplies SPLINTER itself but *not* these `create_spline`/`eval_spline`/… entry points, so this file must stay in `niceops_core`'s sources or `trajectory.cpp` fails to link.

Both `extern "C"` seams above are leftovers of the old D FFI boundary, not a current requirement.
- `src/trajectory.{hpp,cpp}` — `Trajectory` class. Parses the 27-column CSV into `State` records (metric), computing non-dimensional flight parameters: CT = LF·W/(ρ·V_tip²·πR²), μ = v/V_tip, and flight path angles `fpa` (vertical) vs `fpa_e` (relative to the TPP normal). The `qsam` CLI flag selects `fpa_e` (default on). `adapt()` spline-smooths the trajectory via `src/splines.hpp` over vendored splinter (currently commented out in `main`), as is `watts_trim()`. Note `adapt()` is still *compiled* even though nothing calls it, which is why the splinter wrapper is a hard link-time dependency.
- `src/plotting.{hpp,cpp}` — `plt::` facade replacing the vendored matplotlib-d binding. Embeds Python via pybind11 and calls matplotlib directly; the interpreter is created lazily on first use, so non-plotting runs never pay the import cost.

### D standard-library replacements

The port swapped each D dependency for a small local equivalent; when a behavior looks odd, check the D-side original named in the header comment.

| D original | C++ replacement |
| --- | --- |
| `std.getopt` | CLI11 |
| `std.json` | nlohmann_json |
| `std.csv` `csvReader!(double)` | `src/util/csv.hpp` |
| `std.experimental.logger` | `src/util/logger.hpp` (`log`/`info`/`trace`/`fatal`) |
| `std.numeric.findRoot` | `src/util/rootfind.hpp` (bisection) |
| `std.math.approxEqual` | `approx_equal()` in `sphere_database.hpp` |
| `unittest` blocks | Catch2 cases in `tests/` |
| `scope(exit) nc_close` | `NcGuard` RAII in `sphere_database.cpp` |
| matplotlib-d | pybind11-embedded matplotlib |
| `netcdf.d` | the netCDF C API directly |

## Code Conventions

### Porting pitfalls already hit (do not reintroduce)

The port was validated against the D reference and these four bugs were found and fixed. Each was a plausible-looking mistranslation:

- **`iota` is EXCLUSIVE of its end value.** D's `iota(-1500.0, 1500.0, 25.0)` yields **120** values ending at **1475**, not 121 ending at 1500 (verified by compiling D). This governs the observer grid and the contour levels. Comments in the port previously asserted the opposite.
- **Don't sort `phis_scattered`/`thetas_scattered` in place.** D builds the unique axis values from a *copy* (`.dup.sort().uniq().array()`). Sorting in place collapses the per-channel arrays that the grid fill pairs index-wise with `dBAs_scattered`, leaving most cells zero. Symptoms: `dba_max` too low, database load suspiciously fast (~0.3 s instead of ~2.9 s), `--self_test` returning exactly 0, levels tens of dB low.
- **D default-initializes `double` to NaN, not 0.** `ambient_level` relies on this: the `isNaN` guard means D adds no ambient unless `-a` is given. Initializing it to `0.0` adds a spurious 0 dB floor.
- **`approx_equal` mirrors D's `approxEqual`**: relative term against the *second* argument, defaults 1e-2 relative / 1e-5 absolute. It is used for database channel matching, so tightening it silently drops channels.

General rule: when a D construct looks obvious, check it against Phobos semantics rather than intuition — and prefer verifying with a throwaway `ldc2` program, which is how the `iota` question was settled.


- **Port fidelity outranks correctness.** Several `NOTE:` comments mark places that faithfully reproduce D quirks; leave them alone unless you are also regenerating the golden outputs, and say so explicitly:
  - `src/niceops.cpp:208` — `ambient_level` is zero-initialized as D did, so the `!isNaN` guard is always true and every default run mixes in 0 dB "ambient" noise.
  - `src/trajectory.hpp:76` — `utcstart`'s isNaN guard is dead code for the same reason.
  - D's `iota` is inclusive of the stop value; observer grids and contour levels reconstruct that off-by-one deliberately (`src/niceops.cpp:423`).
  - `write_levels` emits one extra trailing newline because the D original did.
  - `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does not run on fatal paths.
- `using namespace niceops;` at the top of the `src/*.cpp` files is intentional (marked `// NOLINT`) so ported code reads like the D original.
- Logging goes through `src/util/logger.hpp`; `main` points it at `log/niceops.log`. The logger's `print` overloads handle vectors and arrays inline, as D's did.
