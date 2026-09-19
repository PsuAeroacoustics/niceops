# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NICE-OPS (Noise Informed Community Environment Operation Planning System) is a helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer grid (or scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

The implementation is C++17, under `src/`. It began as a port of a D original, which is why several comments reference D semantics — those are load-bearing: the port reproduces the D reference bit-for-bit on the demo case, and the goldens come from it.

- `src/` — the implementation. **Every `.cpp` lives here and is listed in CMakeLists.txt**; there are no C++ sources at the repo root.

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

Build out-of-source, and **always with explicit `-S`/`-B`**:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

> **Use `-S`/`-B`, never `cmake <srcdir>` from a build directory.** The repo root still holds a stale `CMakeCache.txt` from the D era. When the directory you hand `cmake` contains a cache, CMake treats that argument as the *build* tree — so `cd build && cmake ..` silently configures **into the source root**, reusing the 2023 cache. That cache carries paths from an old Homebrew prefix (`/usr/local/...`) and produces baffling, unrelated errors such as `Imported target "CGAL" includes non-existent path "/usr/local/include"`. Deleting the root `CMakeCache.txt`/`CMakeFiles/` disarms the trap permanently.

Notes:
- The repo root contains stale in-source build artifacts from the D era (CMakeCache.txt, Makefile, `libniceops*.a`, binaries) that are gitignored — do not build in-source.
- **Never add `-DNDEBUG` to Release.** The D original ran with asserts live, and `sphere_database.cpp` relies on `assert(isFinite(...))` guards; matching that is part of the port contract. This is called out in CMakeLists.txt.
- If Xcode proper isn't installed (Command Line Tools only), pass `-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`. No special flags are needed; the vendored SPLINTER library that used to require them is gone.
- CGAL 5.0+ is header-only and exports `CGAL::CGAL`. Do **not** reintroduce `include(${CGAL_USE_FILE})` — `UseCGAL.cmake` was removed in CGAL 6, leaving `CGAL_USE_FILE` empty so the `include()` errors out.
- **CI builds on GCC/libstdc++; a macOS clang build will not catch what that finds.**
  Two classes of bug are invisible here and fatal there: headers libc++ drags in
  transitively but libstdc++ does not (`<algorithm>` is the usual one), and static-library
  link cycles, which ld64 resolves and GNU ld does not. Before pushing anything structural,
  build the tree with Homebrew GCC rather than waiting for a CI round-trip:

  ```sh
  cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER=$(brew --prefix)/bin/g++-16 -DNICEOPS_PLOTTING=OFF
  cmake --build build-gcc -j8 && ctest --test-dir build-gcc --output-on-failure
  ```

  The goldens pass under GCC too, so a compiler change does not move the 1e-9 comparison
  despite the different FMA contraction. **Do not assert bit-identity against a libm call
  on a literal argument** in a test: the literal gets constant-folded at build time (GCC
  folds via MPFR, correctly rounded) while the function under test reaches libm at run
  time, and the two differ by an ULP. `acosc` was written that way and only GCC noticed.
- **`state/`, `integrated/`, `log/` and `figures/` are not tracked**, so a fresh checkout
  has none of them and every writer must create its own parent directory
  (`ensure_parent_directory` in `src/util/paths.hpp`). `Trajectory::write_state` once
  opened its ofstream directly and so aborted on a clean CI checkout while passing on any
  working tree that had run the tool before. A working tree is therefore *not* a faithful
  test of a checkout — verify against a clone:

  ```sh
  git clone --no-hardlinks . /tmp/fresh && cd /tmp/fresh
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j8
  ctest --test-dir build --output-on-failure
  ```

  Note what such a failure looks like: `fatal()` logs its reason to `log/niceops.log` and
  calls `abort()`, so **nothing reaches stdout or stderr**. CTest reports only
  `Subprocess aborted` with an empty output block, and the reason is in the log alone.
  CI dumps the log on failure for this reason.
- **Keep the C++ sources in one static library.** `interpolator.cpp` used to be its own
  `interpolator` target, a D-era seam that walled off CGAL's heavy headers. Once CGAL went,
  that split was a *cycle*: `interpolator.cpp` calls into `Delaunay`, which lives in
  `niceops_core`, while `niceops_core` linked `interpolator`. GNU ld walks archives in a
  single pass and fails (`undefined reference to niceops::Delaunay::barycentric`); ld64
  resolves it either way, **so this breaks Linux only and a macOS build will not show it**.
  If you add a target, check that its dependencies point one direction.

### Golden reference files

`golden/197421_LAE.csv` and `golden/state_197421.csv` are captured from the **D reference implementation** and the C++ port matches them to **1e-9** (actual max relative difference ~1e-14, i.e. a few ULPs of LDC-vs-clang round-off). `regression_golden` and `state_golden` both pass.

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

CTest defines five tests, all with the **source root** as their working directory (the demo case's relative paths resolve there):

- `regression_default` — runs `niceops -S state/ci_state.csv`; one run feeds both golden comparisons.
- `self_test` — runs `niceops --self_test`, validating every sphere in the database against its own stored levels (1,844,672 samples). **Passes.** It was `broken_self_test`, excluded from CI, while legacy hover directivity made it fail on 7,788 samples; that is fixed (see below), so it is now a real test and CI runs the full suite with no `-E` exclusion.

- `unit_test` — Catch2 cases in `tests/`, ported from the D `unittest` blocks. All 7 cases pass (~26k assertions). Note Catch2 rejects `&&` and chained comparisons inside `REQUIRE` unless the whole expression is wrapped in a second pair of parens: `REQUIRE((a && b));`.
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

Default demo case (no arguments): trajectory `tracks/S-76D/197421AC.csv`, database `databases/S-76D_M3.nod`, vehicle config `databases/S-76D.json`, output `integrated/197421_LAE.csv`. For the full option set see `README.md`, which documents every flag, or the CLI11 setup in `src/niceops.cpp`.

Input/output files use **feet**; internal computation is metric (`ft2m = 0.3048`). Trajectory CSVs are the 27-column NASA/FAA AC state format. Vehicle JSON keys: `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg, converted to Newtons), `Atmosphere.density`.

## Architecture

Pipeline in `main()` (`src/niceops.cpp`): load `Database` + `Trajectory` → for each observer and each trajectory state, compute instantaneous dBA via `Database::level_at_observer` → integrate the time history into a sound exposure level (SEL) using the 10-dB-down-from-peak window (`peak_range`/`sound_exposure_level` in `src/util/sel.hpp`) → write CSV output and optional matplotlib figures.

### Why the hot path is shaped the way it is

A default run evaluates 14,400 observers × 6,551 states = 94.3M times, so anything left in `eval` costs ~100 M×. Do not "simplify" these back:

- **`Database::SourceState` is built once per trajectory state, in `main`, and reused across every observer.** It carries the source position, the `SphereFrame` (unit TPP normal, nose, starboard) and the `Blend` (four bracketing spheres + barycentric weights). The `Blend` is a CGAL Delaunay `locate()` walk and was ~60% of the whole run when it sat inside `eval`. Anything that depends only on (μ, γ, CT) or on attitude belongs here, not in `eval`.
- **`lookup_angles` is split**: `sphere_frame(velocity, normal)` is per state, `lookup_angles(frame, direction)` is per observer. The four-argument overload still exists for the unit tests and composes the two.
- **The level grid is one flat `vector<Sample>`** (`{dBA, EAA}` pairs), indexed `(i_phi * n_theta + i_theta) * n_conditions + condition`, not `vector<vector<vector<double>>>`. Both metrics are read together, so they share a cache line. A condition-major layout was tried and is *worse* — it also breaks bit-exactness via different FMA contraction.
- **`GridAxis` replaces a binary search** with arithmetic when the axis is uniform, then corrects the estimate against the neighbouring knots so the answer provably equals `upper_bound`'s. The searches alone were about a third of the loop. `tests/test_sphere_database.cpp` sweeps both paths against each other; keep that test if you touch it.
- **`trace()` is force-inlined and `Logger::write` is `noinline`** (`src/util/logger.hpp`). `eval` and `level_at_observer` make ~7 trace calls each; as ordinary calls that was a sixth of the loop spent entering functions that immediately returned. The level lives in `detail::min_level`, outside the `Logger` singleton, so testing it does not touch a function-local static's guard.
- **The database grid fill walks channels and places them**, rather than scanning all 1,406 channels for each of the 1,332 cells (2.5 billion `approx_equal` probes, most of the load time). Two semantics are load-bearing and preserved: a later channel still wins a shared cell, and `dba_max` covers matched channels only.

Current cost on the demo case: database load ~0.65 s, acoustic loop ~3.2 s. It was 2.2 s and 11.7 s. What is left is roughly 40% grid gather and arithmetic in `eval`, 20% `acos` (two per evaluation, irreducible without changing results), 10% `log10`.

Two suggested optimizations were deliberately **not** taken:

- Storing reciprocals of the (max − min) normalizing ranges: a reciprocal multiply is not the same as a divide, so it would move the goldens, and the normalization is per state now anyway.
- Hoisting the output-directory check out of `write_time_history_file`: one `exists()` next to creating and writing a ~130 KB file is noise, and it would mean bypassing `open_output`'s guarantee.

**The per-state slab was built and measured; do not rebuild it.** Collapsing the database to one blended 36×37 slab per trajectory state (with the observer loop blocked so `state` can be the outer loop) works and is bit-exact, but it buys only 6.8% of the acoustic loop for +131 MB of resident memory, and the condition-major grid layout it requires makes the ordinary blend path — `self_test`, and any observer set below the break-even block size — about 6% *slower*. Block size barely moves it, and a single 768 MB block is the worst configuration of all, because the level buffer's own bandwidth then dominates. The work is preserved on the `slab-experiment` branch with the full numbers.

The reason the intuition fails is worth keeping: observers sweep a 120×120 ground grid, so consecutive evaluations land on neighbouring (phi, θ) cells and the gather is already cache-resident. It was never miss-bound. Measuring the ceiling first — point `eval` at a fixed L1-resident slab and time it — costs ten minutes and would have shown this before any restructuring.

**The SEL integral takes its time step from the trajectory, not from `-s`.** `sound_exposure_level` weights every sample by a time step, and `-s` is only the resampling step for `adapt()` — the two coincide *only* when `adapt()` actually resamples onto that grid, and it is commented out in `main`. Every shipped track is sampled at 0.02 s, so using the 0.5 s default as the weight overstated every level by exactly `10*log10(0.5/0.02) = 13.98 dB`. `main` now derives the step via `uniform_sample_interval()` and logs it, warning if the sampling is not uniform. If you re-enable `adapt()`, it resamples to `-s` and the derived step follows automatically.

- `src/sphere_database.{hpp,cpp}` — `Database` class. Loads a netCDF noise database (`.nod`): one group per operating condition (μ advance ratio, γ flight path angle, CT thrust coefficient), each holding scattered (phi, θ) channels of dBA level and EAA excess-attenuation factor. `eval()` computes ground level: geometry → sphere angles (`lookup_angles`), barycentric interpolation between operating conditions, periodic bilinear interpolation on the sphere grid, then spherical spreading + excess attenuation (EAA × range) + a fixed 4 dB ground effect. Only `same_grid == true` databases are supported (the interpolation path is unimplemented). Free functions (`bracket`, `periodic_bilinear_interp`, `lookup_angles`, the geometry helpers) are exported in the header so the unit tests can reach them.
- `src/predicates.{hpp,cpp}` — exact orientation and in-sphere tests via Shewchuk expansions, with a filtered double fast path. **The exactness is not optional**: the operating conditions are a prism (164 flight conditions x 8 thrust coefficients), so whole families of points are *exactly* cospherical and a rounded in-sphere test answers the same question two different ways, which leaves Bowyer-Watson's cavity non-star-shaped. Assumes normal (non-subnormal) coordinates; conditions are normalized into [0, 1] first.
- `src/delaunay.{hpp,cpp}` — Delaunay triangulation in arbitrary dimension, Bowyer-Watson on those predicates. Dimension is a constructor argument, not a template parameter, so a fourth operating parameter costs nothing here. A simplex whose circumsphere passes exactly through the new point is left standing: a point *on* a circumsphere does not violate the Delaunay condition, and leaving it is what keeps the choice consistent under mass degeneracy.
- `src/interpolator.{hpp,cpp}` — barycentric weights over the operating-condition space, on top of `delaunay.hpp`. Replaced CGAL's `Delaunay_triangulation_3`, and Boost/GMP/MPFR went with it. Databases with fewer than `dim + 1` conditions (a single exported sphere) skip the triangulation and take the nearest condition.
- `src/atmosphere.{hpp,cpp}` — ISO 9613-1:1993 atmospheric absorption, ported term for term
  from panam's `panam_acoustics/iso_9613_1_1993.py` and pinned against values generated by
  running that module (`tests/test_atmosphere.cpp`). **Pressure is in kPa**, matching the
  standard and the Python, and is the one non-SI unit in the codebase. `absorption_coefficients()`
  is the form the hot path uses: it hoists the two relaxation frequencies out of the per-band
  loop and is bit-exact against the per-frequency call.
- `src/weighting.{hpp,cpp}` — A-weighting (panam's `dBAw`, the IEC 61672 pole-zero form rather
  than a table, so it evaluates at whatever band centres a database carries) and band-limited
  energy integration. Note the curve is **not** identically 0 dB at 1 kHz in this form but
  -1.49e-6 dB; a test asserting exact 0 there will fail. `a_weighting_gains()` folds the
  weighting and the band limit into one per-band multiply, with excluded bands exactly zero.
  `weighted_level()` returns -infinity for silence, not 0 dB, so an excluded spectrum cannot
  survive an energy sum as a quiet-but-real level.
- **Two propagation paths, and they do not agree.** The default carries one broadband dBA and
  attenuates it with the database's `EAA` factor; `--spectral` carries the per-band spectrum
  and attenuates each band at its own ISO 9613-1 coefficient. `Database::eval_spectrum` mirrors
  `eval`'s bracketing and bilinear cell exactly, so a difference between them is the physics and
  never the interpolation. `--spectral` is opt-in: the goldens are the broadband path, and
  `spectra_` is condition-major (right for `interpolate_sphere`, wrong for this loop), so the
  spectral path pays a long stride per band rather than making the sphere export pay.
  It is rejected together with `--db_cutoff`, whose range bound is derived from the broadband
  level and does not bound the band-limited one.
- `src/smoothing_spline.{hpp,cpp}` — natural cubic smoothing spline (Reinsch's algorithm; interior second derivatives solve a symmetric pentadiagonal system in O(n)). Two ways to pick how hard to smooth: `with_smoothing(lambda)` takes the penalty weight directly, `within_tolerance(tol)` solves for the lambda whose weighted squared residual just meets `tol`, as MATLAB's `spaps` does. `src/splines.hpp` wraps it as `Spline` for `trajectory.cpp`. Replaced the vendored SPLINTER library, which brought a bundled Eigen 3.2.8 that predates C++17.

Both `extern "C"` seams above are leftovers of the old D FFI boundary, not a current requirement.
- `src/trajectory.{hpp,cpp}` — `Trajectory` class. Parses the 27-column CSV into `State` records (metric), computing non-dimensional flight parameters: CT = LF·W/(ρ·V_tip²·πR²), μ = v/V_tip, and flight path angles `fpa` (vertical) vs `fpa_e` (relative to the TPP normal). The `qsam` CLI flag selects `fpa_e` (default on). `adapt()` resamples onto a uniform dt grid, smoothing via `src/splines.hpp` (currently commented out in `main`, as is `watts_trim()`). It takes either a smoothing parameter or a `tolerance` bounding how far the smoothed path may sit from the measured one.
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
- **D default-initializes `double` to NaN, not 0.** `ambient_level` relies on this: the `isNaN` guard means D adds no ambient unless `-a` is given. Initializing it to `0.0` adds a spurious 0 dB floor.
- **`approx_equal` mirrors D's `approxEqual`**: relative term against the *second* argument, defaults 1e-2 relative / 1e-5 absolute. It is used for database channel matching, so tightening it silently drops channels.

General rule: when a D construct looks obvious, check it against Phobos semantics rather than intuition — and prefer verifying with a throwaway `ldc2` program, which is how the `iota` question was settled.


- **Port fidelity outranks correctness.** Several `NOTE:` comments mark places that faithfully reproduce D quirks; leave them alone unless you are also regenerating the golden outputs, and say so explicitly:
  - Contour defaults are 65/100, not c74bab3's 45/90: the older pair puts its top level below the demo case's 93.6 dB peak, leaving it unfilled, and the shipped reference figures use 65/100. Figures only -- the goldens are CSV.
  - `watts_trim` and the observer grid extent are now CLI options (`--trim`, `--grid`), both defaulting to the port's behaviour so the goldens are untouched. With `--trim -6000 --grid -2000 1100 -1500 1600` the demo case reproduces the June 2020 reference figures to within 1.05 dB on every cell.
  - D's `iota` is EXCLUSIVE of the stop value; observer grids and contour levels reconstruct that off-by-one deliberately (`iota_exclusive` in `src/niceops.cpp`).
  - `write_levels` emits one extra trailing newline because the D original did.
  - `level_along_bearing` uses `heading[0]` twice where the second component wants `heading[1]`; that is the D original's latent bug, kept.
  - `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does not run on fatal paths.
- `using namespace niceops;` at the top of the `src/*.cpp` files is intentional (marked `// NOLINT`) so ported code reads like the D original.
- Logging goes through `src/util/logger.hpp`; `main` points it at `log/niceops.log`. The logger's `print` overloads handle vectors and arrays inline, as D's did.
