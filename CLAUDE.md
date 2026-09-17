# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NICE-OPS (Noise Informed Community Environment Operation Planning System) is a helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer grid (or scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

**The project is mid-port from D to C++17.** The layout encodes that split, and it is worth keeping clean:

- `src/` — the live C++ implementation, and the only thing CMake builds (`project(... LANGUAGES CXX)`). **Every `.cpp` in the port lives here and is listed in CMakeLists.txt**; there are no C++ sources at the repo root.
- `tests/` — Catch2 unit tests plus the `golden_compare` tool.
- Repo root — the D originals (`niceops.d`, `sphere_database.d`, `trajectory.d`, `splines.d`, `netcdf.d`), kept as the reference implementation. Read them to resolve intent, but don't edit them expecting a rebuild. `cmake-d/` and `external/matplotlibd/` are vestigial for the same reason.

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
git worktree add --detach /tmp/dref HEAD
cp niceops.d /tmp/dref/niceops.d          # see below: HEAD's copy is NOT the right oracle
cmake -S /tmp/dref -B /tmp/drefbuild -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_D_COMPILER=$(which ldc2) -DNCDF=/opt/homebrew/lib/libnetcdf.dylib \
      -DCMAKE_CXX_FLAGS=-I/opt/homebrew/include/eigen3 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build /tmp/drefbuild --target niceops
```

Three things will bite you:

- **Use the working-tree `niceops.d`, not `HEAD`'s.** The port was derived from the *uncommitted* version: contour defaults 45/90 (HEAD has 65/100), and `watts_trim()`/`adapt()` both commented out. HEAD's copy runs trimming and spline adaptation, computes a different trajectory, and yields a golden the port can never match.
- **D prints only 6 significant figures** by default, so a stock D build agrees with the port to just 5e-6. For a 1e-9 golden, patch the writers to `%.17g` (`write_levels` in `niceops.d`, `write_state` in `trajectory.d`) — formatting only, no numerics.
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
- `broken_self_test` — runs `niceops --self_test`; **known to fail, exactly as the D original does**. Both implementations fail on **7,788 of 1,844,672 samples (0.42%)** — the counts match exactly, so this is inherited behavior, not a port defect. The D-era TODO at `src/niceops.cpp:276` describes one such corner case. CI excludes it with `ctest -V -C Release -E broken`; don't rename it without updating `.github/workflows/ci.yml`.

  Historical note: this test briefly failed on 99.8% of samples with `eval` returning exactly 0. That was the in-place `phis_scattered`/`thetas_scattered` sort bug (see Code Conventions), now fixed. If the failure count jumps far above 7,788 again, suspect the database grid fill first.

  **Root cause of the remaining 7,788 (diagnosed, not yet fixed).** Every failure is at **μ = 0** — hover. All 24 hover groups in `S-76D_M3.nod` fail; no non-hover group fails. μ = 0 is the minimum of the advance-ratio range and is isolated by a large gap (next μ is 0.10617), so after normalization those points sit exactly on the convex-hull boundary of the operating-condition cloud that `interpolator.cpp` triangulates. `barycentric()` has an explicit branch for queries whose located cell is infinite (falling back to `nearest_vertex`), and a hull-boundary query is precisely where CGAL's `locate()` is ambiguous between that branch and the finite-cell volume-ratio path. Errors reach ±5.8 dB, median −0.85 dB, which is consistent with weights blending neighbouring conditions instead of returning the vertex's own value.

  Two concrete defects to address when fixing this (both inherited — `interpolator.cpp` was shared with the D build, which is why D fails identically):
  - `barycentric()` does not special-case a query that coincides with a triangulation vertex. Self-test queries always sit exactly on one, and should return weight 1 for it.
  - `total_volume` in `barycentric()` is declared `float` (`auto total_volume = 0.0f;`) inside an otherwise double-precision computation, capping every barycentric weight at ~7 significant digits. This is independent of the hover problem and worth fixing regardless.

  Fixing either changes numerical output and invalidates `golden/`, so it is deliberately staged after the faithful-port work.

- `unit_test` — Catch2 cases in `tests/`, ported from the D `unittest` blocks. All 7 cases pass (~26k assertions). Note Catch2 rejects `&&` and chained comparisons inside `REQUIRE` unless the whole expression is wrapped in a second pair of parens: `REQUIRE((a && b));`.
- `regression_golden`, `state_golden` — `golden_compare` diffs the run's CSV output against a golden file elementwise at 1e-9 relative tolerance, skipping non-numeric header lines. See the caveat above: the golden files are missing.

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
