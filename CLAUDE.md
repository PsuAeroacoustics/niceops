# CLAUDE.md

Rules and invariants for Claude Code (claude.ai/code) in this repository. NICEOPS predicts
ground-level dBA noise exposure over observers from a helicopter trajectory and a rotor-noise
"sphere" database (`.nod`, built by panam). It is C++17 in `src/`, ported from a D original
whose semantics the comments reference; those comments are load-bearing. The D sources remain at
`c74bab3` (on `origin/no-stk`) for use as an oracle (see the `regenerate-goldens` skill and
`docs/HISTORY.md`).

**Read the matching `docs/internals/` file before you change a module.** Most of what looks
over-engineered there was measured, and many obvious-looking simplifications were tried and
failed. The module map, pipeline and repository layout are in
[docs/internals/README.md](docs/internals/README.md#module-map).

## Where documentation goes

- `README.md`: description, build, first run, documentation index. Keep it at 100 lines or fewer.
- `docs/MANUAL.md`: how to use each option (what, when, limits, cost, refusals). No equations
  beyond one-liners, no measurements, no test or code names, no history.
- `docs/THEORY.md`: the equations as run, with references. A change to a model's equations
  updates it in the same commit.
- `docs/FIDELITY.md`: every measured or validation result, once. Other docs link to it.
- `docs/HISTORY.md`: anything that moved results. "Used to" and "before" wording goes only here.
- `docs/internals/`: things to preserve, rationale, rejected designs, measured costs, guard tests.
- Never put dates, commit hashes, harness paths or session narrative in MANUAL or THEORY.

## Build

Commands and options are in `README.md`; rationale is in
[docs/internals/build.md](docs/internals/build.md). Build out-of-source, **always with explicit
`-S`/`-B`**.

- **Always test with asserts on** (`-DNICEOPS_ASSERTS=ON`, which strips `-DNDEBUG`;
  `tests/test_build.cpp` fails if NDEBUG survives). Every assert is side-effect free; keep it so.
- **CI builds on GCC/libstdc++, and a macOS clang build won't catch what GCC finds**: missing
  transitive headers (`<algorithm>`) and static-library link cycles, which ld64 resolves and GNU
  ld doesn't. Before you push anything structural, build and test under Homebrew GCC:

  ```sh
  cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON \
        -DCMAKE_CXX_COMPILER=$(brew --prefix)/bin/g++-16 -DNICEOPS_PLOTTING=OFF
  cmake --build build-gcc -j8 && ctest --test-dir build-gcc --output-on-failure
  ```

- **Every `.cpp` lives in `src/` (none at the repo root) and is listed in CMakeLists.txt**, in
  the right target. A file that compiles but is never linked shows up as an undefined-symbol
  error far from its cause.
- Plotting needs `python3-venv`. Without it the build's environment step only warns and falls
  back to the system interpreter's packages, so CI checks the environment explicitly.
- **Keep the C++ sources in one static library.** A separate target that calls back into
  `niceops_core` makes a cycle, which breaks Linux only. `-DNICEOPS_PLOTTING=OFF` is the only
  build that compiles `src/plotting_disabled.cpp`.
- **`state/`, `integrated/`, `log/` and `figures/` are not tracked.** Every writer must create
  its own parent directory (`ensure_parent_directory`, `src/util/paths.hpp`). A working tree is
  therefore not a faithful test of a checkout: verify in a clone (`git clone --no-hardlinks .
  /tmp/fresh`, then build and run ctest there).
- **Don't assert bit-identity against a libm call on a literal argument** in a test. GCC
  constant-folds the literal through MPFR, and the two results differ by an ULP.

## Tests

- Run the whole suite as CI does: `ctest --test-dir build --output-on-failure`, with no `-E`.
  CTests run with the **source root** as working directory. For speed, run Catch2 directly
  (`build/niceops_unit_tests "[sphere_database]"`, by tag or name). Never write a test count
  into a doc; `ctest -N` gives it. `regression_default` runs the demo case and feeds
  `regression_golden` and `state_golden`, which use `golden_compare` (elementwise, 1e-9
  relative, header lines skipped).
- Catch2 rejects `&&` and chained comparisons inside `REQUIRE` unless the expression is wrapped
  in a second pair of parentheses: `REQUIRE((a && b));`.
- **Goldens:** regenerate one only for a deliberate numerics change, and say so in the commit.
  `golden/state_197421.csv` is D-derived (use the `regenerate-goldens` skill); every other
  golden comes from the C++ build. The table is in
  [docs/internals/testing.md](docs/internals/testing.md#golden-files).
- **Several tests exist only to kill a specific mutant**, because the rest of the suite passed
  with that mutant in place. Keep them when you refactor; the list is in
  [testing.md](docs/internals/testing.md#mutation-guarded-tests). If a mutation check restores
  a file with `mv`, `touch` the file afterward or make won't rebuild it.
- Tests of refusals and warnings use `PASS_REGULAR_EXPRESSION` on the message text.
- **No CTest exercises plotting** (`-p`/`-f`/`-F`); CI runs a separate figure step.

## Invariants

- **The default path is byte-identical to the goldens**, and new options must leave it that
  way. A degenerate case of a new model, such as calm air, uniform air or level terrain, should
  reproduce the existing model, at 1e-9 or bit-for-bit. A CTest should pin that.
- **The windowed peak search is exact.** It skips only states that a bound *proves* are below
  the threshold. Any new gain (ground, table, focusing, plate) must enter the bound, or
  `*_windowed_is_exact` fails. Options that need every state (history outputs, time weighting,
  day metrics, EPNL, `-a`, and others) march instead. See
  [performance.md](docs/internals/performance.md).
- **The hot loop:** a default run makes about 94 million evaluations, so anything left in the
  inner loop costs about 100 million times over.
- **`--spectral` is the default whenever the database carries spectra.** The goldens use the
  broadband path, because `S-76D_M3.nod` has no spectra (`regression_default` also passes
  `--no-spectral --no-ground_effect` explicitly). `eval` and `eval_spectrum` share their
  cell and stencil, so any disagreement between them is physics, not interpolation.
- **The ground model returns a gain, not a loss.**
- **Compass bearings** are turned into the local frame by `--frame_bearing` through one
  rotation, `CompassToLocal`. The shipped AC tracks need 210.
- **panam is the reference for ported physics**: ISO 9613-1, `ega()`, EPNL, and A-weighting are
  each pinned to panam at 1e-9 dB or similar. A change there is a change to panam first. A
  wrong result may also be bad database data, not bad code.
- **Units.** Computation is metric internally, and the conversion is one seam in `run()`.
  Defaults are physical quantities written in meters and skip the conversion. `--units`
  declares what the files already are and converts nothing. `write_state` stays metric, because
  the state golden encodes internal values. Speeds name their own unit (`12kt`, `6.2mps`), and
  a bare number is refused.
- **Input checking.** A bad input must fail or warn with its cause named, never give a plausible
  footprint with exit status 0. Errors and warnings go to stderr as well as `log/niceops.log`; a
  skipped trajectory exits 1. Readers refuse non-increasing time, non-finite values, unit
  mismatches, text after data, and observer files that yield nothing. A new reader follows the
  same rules: [docs/internals/inputs.md](docs/internals/inputs.md).

## Code conventions

- **Port fidelity outranks correctness, with a bar for departing from it.** Leave the `NOTE:`
  comments that reproduce D quirks alone unless you are also regenerating the goldens, and say
  so if you do. Departing from D takes *a measured case that the D behavior is wrong*, not a
  reading that it looks wrong. The departures that met that bar are listed in
  [HISTORY](docs/HISTORY.md#departures-from-d).
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
- Flags are defined in the CLI11 setup in `src/niceops.cpp`; document a new one in the MANUAL.
- **Licensing:** the project is MIT. `THIRD_PARTY_NOTICES.md` holds verbatim notices for the
  fetched dependencies; regenerate it from `build/_deps` when one changes. The `Agg` backend in
  `plotting.cpp` is load-bearing: the Qt backends would pull in GPL PyQt.
