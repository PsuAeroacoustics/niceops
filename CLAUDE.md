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

Link-time optimization is on by default (`-DNICEOPS_LTO=OFF` turns it off), where
`check_ipo_supported` says the toolchain has it. It is worth 18% of the acoustic loop,
because that loop and the database live in different translation units; it leaves the
golden comparison passing unchanged, and GCC 16 and Apple Clang produce identical output
with it.

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

### Where the database load time goes

A default run is now ~1.0 s, and ~0.75 s of that is `Database::Database`. Most of it is not
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

The first three -- 78% -- are netCDF's own cost for this file layout and are not reachable
through its API; the five tiny scalars per group cost more than the four big arrays, because
what they pay for is the lookup, not the bytes. The real fix is in the database builder
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

**`--units` picks the length unit for the whole interface** -- `ft` (default) or `m` -- covering
every length on the command line (`--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
`--receiver_height`, `--path_smoothing`), every file carrying a length, and the figure axes.
Internal computation is metric (`ft2m` and `LengthUnit` in `src/util/constants.hpp`); the
conversion is one seam in `run()`.  Four of those options used to be metric while every file was
feet, which is the inconsistency this removes.

Three things are load-bearing here:

- **A default is a physical quantity, not a number in the current unit.**  Defaults are written
  in metres and skip the conversion, which is why the seam tests CLI11's `count()` rather than
  scaling everything.  Scale the defaults too and `--units m` shrinks the observer grid to a
  third of its extent.
- **`--units` declares what the files already are; it converts nothing.**  The shipped tracks are
  feet, so `--units m` on the demo case reads feet as metres.  Nothing can detect that.  Where a
  column header names its own unit (`alt_ft`, `alt_m`) it is checked against `--units` and a
  disagreement fails the file; a bare `alt` is always accepted.
- **`write_state` stays metric** whatever `--units` says: `golden/state_197421.csv` is the
  D-derived oracle and encodes the internal values.

`fps2g` is chosen as a whole constant rather than composed as `unit.to_metres * mps2g`, because
the product of two constants is not the same double as two multiplies and the state golden is
pinned at 1e-9.  Verified end to end: a metric twin of the demo track under `--units m`
reproduces the golden footprint to 1.5e-12 dB. Trajectory CSVs are the 27-column NASA/FAA AC state format. Vehicle JSON keys: `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg, converted to Newtons), `Atmosphere.density`.

## Architecture

Pipeline in `main()` (`src/niceops.cpp`): load `Database` + `Trajectory` → for each observer and each trajectory state, compute instantaneous dBA via `Database::level_at_observer` → integrate the time history into a sound exposure level (SEL) using the 10-dB-down-from-peak window (`peak_range`/`sound_exposure_level` in `src/util/sel.hpp`) → write CSV output and optional matplotlib figures.

### Why the hot path is shaped the way it is

A default run evaluates 14,400 observers × 6,551 states = 94.3M times, so anything left in `eval` costs ~100 M×. Do not "simplify" these back:

- **`Database::SourceState` is built once per trajectory state, in `main`, and reused across every observer.** It carries the source position, the `SphereFrame` (unit TPP normal, nose, starboard) and the `Blend` (four bracketing spheres + barycentric weights). The `Blend` is a CGAL Delaunay `locate()` walk and was ~60% of the whole run when it sat inside `eval`. Anything that depends only on (μ, γ, CT) or on attitude belongs here, not in `eval`.
- **`lookup_angles` is split**: `sphere_frame(velocity, normal)` is per state, `lookup_angles(frame, direction)` is per observer. The four-argument overload still exists for the unit tests and composes the two.
- **`lookup_angles` takes three dot products and no projection.** `{nose, normal, starboard}` is orthonormal, so the propagation direction's components in the plane normal to the nose are just `normal.dot(d)` and `starboard.dot(d)`; forming the projected vector and normalising it (a square root and three divides) scaled both by the same positive length, which the angle does not depend on. `phi` is then `atan(across / |up|)` — `atan`, not `atan2`, which costs twice as much for quadrant handling a non-negative second argument does not need, and gives the right half of pi from the infinity the division produces when `up` is zero. The 0/0 case at the poles is caught explicitly.
- **The level grid is one flat `vector<Sample>`** (`{dBA, EAA}` pairs), indexed `(i_phi * n_theta + i_theta) * n_conditions + condition`, not `vector<vector<vector<double>>>`. Both metrics are read together, so they share a cache line. A condition-major layout was tried and is *worse* — it also breaks bit-exactness via different FMA contraction.
- **`GridAxis` replaces a binary search** with arithmetic when the axis is uniform, then corrects the estimate against the neighbouring knots so the answer provably equals `upper_bound`'s. The searches alone were about a third of the loop. `tests/test_sphere_database.cpp` sweeps both paths against each other; keep that test if you touch it.
- **`DirectionBounds::bin` lives in the header.** The observer loop calls it ~94 million times from `niceops.cpp` while it is defined in `sphere_database.cpp`, so out of line it inlined only when `NICEOPS_LTO` was on and the toolchain supported it. It is fifteen instructions and should not depend on a build option; its tangent table is an `inline const` variable rather than a function-local static, which would test a guard on every call. The bin geometry — `bin`, the tangents, `narrowest_theta_bin_degrees`, `azimuth_bin_degrees` — sits on `DirectionBounds` together, because the table build has to pad by the same binning the hot path queries with, and three independent spellings of it is how that silently stops covering.
- **`GridAxis::locate` gives the cell and the interpolation weight together**, and is what `eval` and `eval_spectrum` call. Going through `bracket()` meant reading four knots back out of `values_`, four subtractions and a divide by the cell area, all to recover a number the uniform spacing already implies; `locate` returns the bracketing pair plus the fraction of the way across the cell, so the four corner weights are products of `(1-f)` and `f` and already sum to one. Worth 18% of the loop. It skips `bracket`'s correction loops, which is why it needs the tighter `exactly_uniform_` test in `reset()` rather than the tenth-of-a-step one: round-off can only put the estimate on the wrong side of a *knot*, where the two cells interpolate to the same value anyway. The wrap-around cell at the ends is not uniform and falls back to the knots. `tests/test_sphere_database.cpp` sweeps `locate` against `periodic_bilinear_interp` over uniform, irregular and wrapping cells — the two derivations have to agree, and that test is what says so.
- **Both `eval` and `eval_spectrum` build their stencil with `cell_stencil`.** It turns the two `GridAxis::Cell`s into the four corner indices and the four bilinear weights. Sharing `locate` alone got them the same *cell* but left each deriving its own weights, which is what keeps the claim below — that a broadband/spectral disagreement is the physics and never the interpolation — true rather than merely intended. `periodic_bilinear_interp` is no longer on either path; it is the reference `locate` is swept against in `tests/test_sphere_database.cpp`, and the cell/apply split it used to have (`BilinearCell`, `bilinear_apply`) existed only for the old `eval` and is gone.
- **The interpolation carries `{dBA, EAA}` as a two-lane vector** (`Pair` in `sphere_database.cpp`). Every step applies the same weight to both metrics, and a `Sample` is exactly a two-double vector in memory, so blending the four conditions at four corners and then the corners themselves is half as many instructions. Clang would not do it on its own from the scalar form or from a two-element array; the GCC/Clang `vector_size` spelling is explicit, with a plain-struct fallback that compiles to what the code said before.
- **Divides in the loop are multiplications by a stored reciprocal**: `inv_sphere_radius_`, `1e-3` for the EAA scaling, and one reciprocal for normalising the propagation direction instead of three divides. A double divide is about four times a multiply. These move the result by an ULP, which is the only thing in this round that does.
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
  way round is what works, and is now the default path. See "Walking out from the peak" below.

### The SEL integral is oversampled, and the window is why that cannot be used

The largest speedup still on the table is not a code change: the sound exposure integral is a
Riemann sum over a 50 Hz trajectory, and it is oversampled by roughly thirty times. Decimating
the demo case's 6551 states moves the footprint's **mean** by 0.008 dB at stride 2 and 0.088 dB
even at stride 64. A 2-8x cut in the whole acoustic loop is sitting there.

What blocks it is `peak_range`, which walks out from the peak and stops at the **first** sample
below peak - 10. One transient directivity dip therefore truncates the window, and the result
is a knife-edge: decimation steps over the dip, the window doubles, and the level jumps. The
tail does not shrink with stride the way a quadrature error would --

| stride | max cell error, first-crossing | outermost-crossing |
| --- | --- | --- |
| 2 | 1.93 dB (104 cells > 0.1 dB) | 1.13 dB (86) |
| 8 | 2.34 dB (619) | 1.13 dB (271) |
| 64 | 2.48 dB (3150) | 1.54 dB (1346) |

-- which is the signature of a discrete window flip, not of the integral. Widening the window
to the outermost crossing only halves it, and is itself a metric change: at full resolution it
moves 2758 of 14400 observers by more than 0.1 dB and one by 2.48 dB.

So decimation is **not** safe as a default under either definition: the mean is negligible but
individual cells move by 1-2 dB, and the cells are the deliverable. It stays available
deliberately, through `adapt()` and `-s` (both currently commented out in `main`).

The oversampling is spent a different way instead, and exactly: by not evaluating the states
the window was never going to use. That is the next section.

Worth knowing on its own, independent of performance: the first-crossing rule means ~0.7% of
observers' SEL depends on whether a single sample dips below the threshold. That is a faithful
port of the D original and is left alone under the port-fidelity rule, but it is a fragility,
not a feature -- and it is the reason the oversampling cannot be spent.

### Walking out from the peak instead of marching every emission time

`windowed_sound_exposure_level` in `src/niceops.cpp` is the default path and the reason the
loop is 1.2x faster than evaluating everything. It is **exact** -- byte-identical output on
every track, grid, observer set and thread count tried, not an approximation.

It does not carry its own copy of the window walk or the integral. `expand_window()` and
`sel_over_window()` in `src/util/sel.hpp` are templated on a level accessor, so `peak_range`
reads them out of a filled array and this path reads them through a memoizing lambda that
computes a state the first time it is asked for. That is deliberate: the backward walk's
`== threshold` tie rule decides ~0.7% of observers' SEL, and it now exists once rather than
in two places that have to be edited together. The accessor is a template, not a
`std::function`, so it inlines away and the comparisons are unchanged.

The idea is to stop asking "can I prove this state is outside the window?", which needs a bound
tighter than the 10 dB window and is hopeless (above), and ask "where is the peak?" instead.
`peak_range` walks out from the peak and stops at the **first** sample below peak - 10, so the
window is contiguous *by construction*: once the peak is known the walk itself decides where
the window ends, and nothing beyond it needs a bound at all. A bound is then needed for exactly
one thing -- being sure the walk started on the true peak -- and that is a far weaker demand,
because it only has to beat the best level found so far rather than that level minus ten.

Three things had to hold, and each was measured before any of it was written:

- **The bound needs less than 10 dB of slack**, or the peak search reaches outside the window
  and there is nothing left to save. The whole-sphere maximum is 11.4 dB, just over. Binning it
  by direction -- `DirectionBounds`, bins uniform in cos^2(theta) crossed with azimuth sectors
  -- brings the median to 8.6 dB and the search collapses inside the window.
- **The test must cost nothing on the states it does not reject.** It forms the squared range
  and the three dot products with the sphere frame, which is exactly what an evaluation needs,
  so `Database::level_from_geometry` takes those and a surviving state pays nothing for having
  been tested. Letting the test compute its own geometry and then recomputing it inside
  `level_at_observer` was measured: it turns a 1.48x into a 1.08x.
- **The table has to stay in cache.** One table per *distinct blend*, not per state -- a shipped
  track resolves thousands of states onto a few dozen tetrahedra -- so it is tens of kilobytes.
  Built per state it was 3.3 MB, re-streamed for every observer, and that alone cost more than
  the scheme saved once the loop was threaded.

What it gives up: the bound maxes over the four conditions of a tetrahedron rather than over
the particular blend of them, so it is looser than the per-state tables the prototype used --
70% of states evaluated against the prototype's 64%, where the window itself is the 56% floor.
Buying that back means a table per state, which is the cache problem above. Bin counts were
swept (8x8, 16x8, 16x16, 32x16): finer bins do cut the evaluated fraction, to 66%, but the
extra binning work cancels the saving almost exactly (1.585, 1.546, 1.557, 1.567 s), so the
choice among them hardly matters and 16x8 is kept only for being marginally ahead on both.

`tests/test_sphere_database.cpp` checks the bound directly -- that it is never exceeded over a
sweep of source states and observers -- because the whole thing rests on that one property.

The path is skipped and the loop marches every state whenever the bound does not apply:
`--spectral` (the bound is on the broadband dBA, which does not bound the band-limited weighted
level -- the same reason `--db_cutoff` is refused there), `-a` (the level is transformed after
eval), `-i` (the whole history is written out), `-D` (the day-night window is found on a
penalized history, which a window found on the raw one need not cover), `--db_cutoff` (already
culling on its own terms), and any database with a negative excess-attenuation factor, where
dropping that term is not a loss.

**Extending it to `--ground_effect` is the obvious next move and is not done.** An evaluation
costs ~1160 ns there against the broadband path's 19, so the test is free and the saving would
be close to the full 30% rather than being half eaten by overhead. It needs a bound on the
*A-weighted band-limited* level -- the largest weighted band sum per direction, which a
spectral database can supply -- and a bound on the ground term, which is a **gain** of up to
about +6 dB and so has to be bounded rather than dropped.

Two suggested optimizations were deliberately **not** taken:

- Storing reciprocals of the (max − min) normalizing ranges: a reciprocal multiply is not the same as a divide, so it would move the goldens, and the normalization is per state now anyway.
- Hoisting the output-directory check out of `write_time_history_file`: one `exists()` next to creating and writing a ~130 KB file is noise, and it would mean bypassing `open_output`'s guarantee.

**The per-state slab was built and measured; do not rebuild it.** Collapsing the database to one blended 36×37 slab per trajectory state (with the observer loop blocked so `state` can be the outer loop) works and is bit-exact, but it buys only 6.8% of the acoustic loop for +131 MB of resident memory, and the condition-major grid layout it requires makes the ordinary blend path — `self_test`, and any observer set below the break-even block size — about 6% *slower*. Block size barely moves it, and a single 768 MB block is the worst configuration of all, because the level buffer's own bandwidth then dominates. The work is preserved on the `slab-experiment` branch with the full numbers.

The reason the intuition fails is worth keeping: observers sweep a 120×120 ground grid, so consecutive evaluations land on neighbouring (phi, θ) cells and the gather is already cache-resident. It was never miss-bound. Measuring the ceiling first — point `eval` at a fixed L1-resident slab and time it — costs ten minutes and would have shown this before any restructuring.

**The SEL integral takes its time step from the trajectory, not from `-s`.** `sound_exposure_level` weights every sample by a time step, and `-s` is only the resampling step for `adapt()` — the two coincide *only* when `adapt()` actually resamples onto that grid, and it is commented out in `main`. Every shipped track is sampled at 0.02 s, so using the 0.5 s default as the weight overstated every level by exactly `10*log10(0.5/0.02) = 13.98 dB`. `main` now derives the step via `uniform_sample_interval()` and logs it, warning if the sampling is not uniform. If you re-enable `adapt()`, it resamples to `-s` and the derived step follows automatically.

- `src/sphere_database.{hpp,cpp}` — `Database` class. Loads a netCDF noise database (`.nod`): one group per operating condition (μ advance ratio, γ flight path angle, CT thrust coefficient), each holding scattered (phi, θ) channels of dBA level and EAA excess-attenuation factor. `eval()` computes ground level: geometry → sphere angles (`lookup_angles`), barycentric interpolation between operating conditions, periodic bilinear interpolation on the sphere grid, then spherical spreading + excess attenuation, `EAA * (range - r_sphere) / 1000` — EAA being
the total A-weighted attenuation over the 1000 m reference distance the database was built at,
which is exactly how panam writes and consumes it. No ground effect: see `--ground_effect`. Only `same_grid == true` databases are supported (the interpolation path is unimplemented). Free functions (`bracket`, `periodic_bilinear_interp`, `lookup_angles`, the geometry helpers) are exported in the header so the unit tests can reach them.
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

  The roadmap this starts (refraction, then urban) is in the Claude Doc "NICE-OPS:
  Atmospheric and Urban Propagation Options".
- `src/refraction.{hpp,cpp}` — `--refraction linear` (roadmap step 4). Per path, a
  least-squares line through the effective sound speed c(z) + w(z).h over the heights the path
  spans (`EffectiveSoundSpeed`: running integrals of c, w_x, w_y and z times each on a 1 m
  table, so a fit costs a few flops per pair), then closed-form circular-arc eigenrays. The
  things to preserve:
  - **Numerically stable forms, not textbook ones.** The circles are centred L = c0/|a| away,
    kilometres to 1e11 m, so every difference of squares of h = L + sigma z is formed from
    kappa(z) = z(2L + sigma z). Reflected legs are kappa/(sum of d), travel time is one atanh of
    the addition formula, and arc length is 2 rho asin(chord/2 rho). Written the obvious way
    they cancel to nothing as the gradient shrinks. Gradients under 1e-9 /s take the straight
    path.
  - **Spreading is r / cos(launch), exactly.** In a linear profile the tube area
    r |dr/dtheta| sin|arrival| / cos(launch) reduces to (r rho / h_s)^2.
    `tests/test_refraction.cpp` checks it against finite differences of traced rays; keep that.
  - **The reflected ray is a 1-D solve** on launch angle, along which the range is monotonic
    (checked across both gradient signs). It starts from the direct ray to the receiver's image
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
  - Shadow zones: the grazing ray's level, less the first creeping-wave mode over rigid ground
    (`shadow_attenuation_rates`), capped at `--shadow_floor`.
- `src/stratified.{hpp,cpp}` — `--refraction stratified` (roadmap step 6a): the effective
  sound speed piecewise linear between the profile's own levels, circular arcs within each
  layer, eigenrays solved per pair. The things to preserve:
  - **Parametrise by the Snell invariant p, in four families.** Launched up or down, arriving
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
  - **One arrival per pair.** Non-monotonic profiles give caustics and several arrivals; the
    solve returns the first family's root. Focusing is capped at `max_focusing_gain_db` (+10).
  - **The log-law finding.** With zero wind at the ground, the upwind limiting ray is fixed by
    c(z_s)/c(0) alone. On `example_windy.csv` a 100 m source's upwind shadow starts at 520 m,
    where the linear fit puts it past 1.2 km; a test pins that. Ray theory probably overstates
    that shadow at low frequencies, so GFPE (step 7) should adjudicate before either model's
    shadow levels are trusted.
  - Costs: 12 s on the demo grid with a point receiver (linear 3.3 s), 35 ms over a
    96-microphone array.
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
  than a table, so it evaluates at whatever band centres a database carries) and band-limited
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
  ground model reaches (Re z ∈ [0.09, 22.3], Im z ∈ [-1.2, 15.6]), measured against SciPy's `wofz`.
  The lower-half-plane reflection `w(z) = 2exp(-z²) - w(-z)` is safe **only** in that domain.
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
    (`sqrt(w)` *is* the numerical-distance parameter, verified over the whole domain) and the
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
    contiguous run, which is always what the limits select because band centres are sorted.
  - The Weideman series behind `faddeeva` runs as four interleaved Horner chains rather than
    one; see `src/faddeeva.cpp`. It was 53% of the ground path on its own and is latency-
    bound, not throughput-bound.
  All of it is exact: the footprint moved by 2.8e-14 dB end to end, and the EGA tests still pin
  to panam at 1e-9 dB. Acoustic loop on the Be407 case, default band limits: broadband 0.34 s,
  `--spectral` 1.44 s, `--spectral --ground_effect` 10.9 s (was 41.4); the Horner split takes
  the last of those to about 9.4 s.

  No database in `databases/` carries spectra, but the path is still exercisable end to end
  from a bare checkout: `write_flat_spectral_database` in `tests/test_spectral_propagation.cpp`
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
  never the interpolation. `--spectral` is opt-in: the goldens are the broadband path, and
  `spectra_` is condition-major (right for `interpolate_sphere`, wrong for this loop), so the
  spectral path pays a long stride per band rather than making the sphere export pay.
  It is rejected together with `--db_cutoff`, whose range bound is derived from the broadband
  level and does not bound the band-limited one.
- `src/geodetic.{hpp,cpp}` — the WGS84 transform: a transverse Mercator about a reference point
  with altitude carried through as z. `LocalFrame(reference)` is east/north;
  `track_aligned(reference, bearing)` turns +x onto a compass bearing;
  `from_east_north_up` applies that rotation without the translation, for velocities; and
  `to_geodetic` inverts the projection by Newton iteration, which is what `--latlon` output
  uses. Ported from the D `worldtransforms.d`/`wgs84tocartesian.d` pair, which had no notion of
  a rotated frame and whose inverse walked fixed 1e-5 degree increments.
- `src/smoothing_spline.{hpp,cpp}` — natural cubic smoothing spline (Reinsch's algorithm; interior second derivatives solve a symmetric pentadiagonal system in O(n)). Two ways to pick how hard to smooth: `with_smoothing(lambda)` takes the penalty weight directly, `within_tolerance(tol)` solves for the lambda whose weighted squared residual just meets `tol`, as MATLAB's `spaps` does. `trajectory.cpp` uses it directly. Replaced the vendored SPLINTER library, which brought a bundled Eigen 3.2.8 that predates C++17.

Both `extern "C"` seams above are leftovers of the old D FFI boundary, not a current requirement.
### Geodetic trajectory input

A trajectory may be given in WGS84 (latitude, longitude, altitude) instead of the 27-column AC
state format; `README.md` documents the format and the options. Things worth knowing before
touching it:

- **The format is decided by width, not by the header.** 27 columns is the positional AC
  format; anything else is read by the names in its header (`read_labelled_csv` in
  `src/util/csv.hpp`). That order is load-bearing: the AC format's *own* header names `lat`,
  `lon` and `alt` columns, so a reader dispatching on the header alone takes a state file for a
  WGS84 track — and its x/y origin, which the file does not record, is about 2.5 km back along
  the track. `first_geodetic_position` applies the same width test, and a test pins it.
- **The AC path is byte-identical.** `read_labelled_csv` keeps the first record as the header
  instead of discarding it, but drops exactly the rows `read_numeric_csv(path, true)` dropped,
  so the goldens are untouched. Verified.
- **`alt_m` is deliberately not a recognised altitude spelling.** Altitude is read in feet;
  matching a metre spelling would read metres as feet and fly the aircraft three times too
  high, where leaving it unmatched fails the file with a message naming the units.
- **Altitude is MSL and z is height above the observers.** `--ground_elevation` is the
  difference. In `tracks/S-76D/197421AC.csv` it is a constant 54.2651 ft, which is how that
  file's own alt and z columns relate.
- **The shipped AC tracks' local frame runs +x along a compass bearing of 210°** — the approach
  heading, not the instantaneous heading in column 8, which wanders over 25°. With that bearing
  the projection reproduces that file's own x/y to under a metre over 3.2 km, which is the
  check that the port of the transform is the one its producer used.
- **One frame per run, settled before anything loads.** Every trajectory in a run shares the
  observers, so a position in one has to mean the same place as in another. `run()` resolves it
  from `--origin` or from the first geodetic track and hands it to each `Trajectory`.
- **Velocity and load factor are derived, and that is the whole accuracy story.** The
  projection is not the lossy part: altitude passes through exactly and x/y to under a metre.
  Differentiating a spline twice for the load factor is. **`--path_smoothing` defaults to 10 ft**
  and is an acceptable deviation in feet RMS per axis, not an abstract weight -- the solver
  spends the whole budget rather than merely staying inside it. Interpolating exactly
  (`--path_smoothing 0`) amplifies the coordinates' 1e-8 degree rounding by 1/step² -- 2500x on
  the 0.02 s demo track -- and moves the footprint by **11.6 dB at the median**; the 10 ft
  default takes that to 0.085 dB. Zero is wrong for both kinds of input, which is why it is no
  longer the default: measured positions turn into rounding noise, and sparse waypoints are a
  statement of intent, so flying exactly through them manufactures manoeuvring that is itself a
  noise source. 10 ft is not the optimum for any one input shape (a dense 50 Hz track does
  better at 0.0656 ft, tuned to its own resolution) but it is the only value that is good on
  all of them at once.

  Those numbers were re-measured after `SmoothingSpline` was fixed to give the lambda = 0 case
  its curvature. Before that, `--path_smoothing 0` returned a spline with zero second
  derivative everywhere, so the load factor was exactly 1 g on every state, CT never moved, and
  the amplification above never happened -- the old "17% and 3 dB" measured that bug. README
  carries the re-measured table for all four input shapes; do not re-derive it.
  README carries the measured table for all four input shapes; do not re-derive it.
- **Resampling never coarsens.** The sound exposure integral is a Riemann sum over whatever
  step the states end up with, so a 50 Hz file keeps its own sample times and only sparse or
  unevenly spaced input is refitted onto the `-s` grid.

- `src/trajectory.{hpp,cpp}` — `Trajectory` class. Parses the 27-column CSV into `State` records (metric), computing non-dimensional flight parameters: CT = LF·W/(ρ·V_tip²·πR²), μ = v/V_tip, and flight path angles `fpa` (vertical) vs `fpa_e` (relative to the TPP normal). The `qsam` CLI flag selects `fpa_e` (default on). `adapt()` resamples onto a uniform dt grid, smoothing via `src/smoothing_spline.hpp` and taking velocity and load factor from the spline's own derivatives (currently commented out in `main`; `watts_trim()` is reachable through `--trim`). It takes either a smoothing parameter or a `tolerance` bounding how far the smoothed path may sit from the measured one.
- `src/plotting.{hpp,cpp}` — `plt::` facade replacing the vendored matplotlib-d binding. Embeds Python via pybind11 and calls matplotlib directly; the interpreter is created lazily on first use, so non-plotting runs never pay the import cost.

### D standard-library replacements

The port swapped each D dependency for a small local equivalent; when a behavior looks odd, check the D-side original named in the header comment.

| D original | C++ replacement |
| --- | --- |
| `std.getopt` | CLI11 |
| `std.json` | nlohmann_json |
| `std.csv` `csvReader!(double)` | `src/util/csv.hpp` (`read_numeric_csv`, plus `read_labelled_csv` for the header-named formats, which have no D original) |
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


- **Port fidelity outranks correctness, with one deliberate exception.** Several `NOTE:`
  comments mark places that faithfully reproduce D quirks; leave them alone unless you are also
  regenerating the golden outputs, and say so explicitly. The exception is the 4x EAA multiplier,
  removed above: it was a physics error rather than a quirk, it was measured against the model it
  stood in for before being removed, and `golden/197421_LAE.csv` was regenerated from the C++.
  That is the bar for departing from D — a measured case that the D behaviour is wrong, not a
  reading that it looks wrong.
  - Contour defaults are 65/100, not c74bab3's 45/90: the older pair puts its top level below the demo case's 93.6 dB peak, leaving it unfilled, and the shipped reference figures use 65/100. Figures only -- the goldens are CSV.
  - `watts_trim` and the observer grid extent are now CLI options (`--trim`, `--grid`), both defaulting to the port's behaviour so the goldens are untouched. With `--trim -6000 --grid -2000 1100 -1500 1600` the demo case reproduces the June 2020 reference figures to within 1.05 dB on every cell.
  - D's `iota` is EXCLUSIVE of the stop value; observer grids and contour levels reconstruct that off-by-one deliberately (`iota_exclusive` in `src/niceops.cpp`).
  - `write_levels` emits one extra trailing newline because the D original did.
  - `level_along_bearing` uses `heading[0]` twice where the second component wants `heading[1]`; that is the D original's latent bug, kept.
  - `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does not run on fatal paths.
- `using namespace niceops;` at the top of the `src/*.cpp` files is intentional (marked `// NOLINT`) so ported code reads like the D original.
- Logging goes through `src/util/logger.hpp`; `main` points it at `log/niceops.log`. The logger's `print` overloads handle vectors and arrays inline, as D's did.
