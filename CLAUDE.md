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
- **CI never builds the Metal backend** (it runs on Linux, where `NICEOPS_GPU=AUTO` is `OFF`).
  Build and test a GPU change on a Mac; the GCC check above covers the disabled backend.
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
- **GPU paths are on by default wherever there is a device, and the CPU is their reference.**
  `niceops --buildings`, `niceops_point_field`, `niceops_urban_table --wave_below` and
  `niceops_fdtd_field` take the GPU unless `--no-gpu`, and fall back to the CPU where the GPU does
  not take the case (`--edge_source`, `--fade`, `--ray_coherence measured`, `--phase_interval`, an
  FDTD grid the device cannot hold); `--gpu` asks for it and is refused where it cannot run. They work in float (phases in
  float-float), are tested against the CPU at a stated tolerance, never bit-for-bit, and never feed
  a golden: **a test pinned to CPU results passes `--no-gpu`**, or it fails on a Mac and passes on
  Linux CI. Kernels are written once, in the dialect of `src/gpu/kernels`, so another backend takes
  the same source: [docs/internals/gpu.md](docs/internals/gpu.md).
- **`--edge_source` is the point-source model's exact diffraction, and is held to exact
  solutions**, not to the UTD rays it replaces: Sommerfeld's half-plane, the wedge series, the
  long block's boundary elements. Its second order is localized by partner corrections that
  subtract the first-order model's own ray, and a vertical edge on the floor is integrated through
  it; both are what make it consistent with the ray set, and removing either moves real cases by
  decibels. See [ground.md](docs/internals/ground.md#the-edge-source-integral-srcedge_source---edge_source).
- **`niceops_fdtd_field` is a reference, never a prediction path.** It is held to exact
  solutions (two rays over a floor, the wedge, an impedance wall) and the GPU march to the CPU's.
  Judge a run's length band by band (`truncation_db`), not by the tail's raw energy, and give
  parallel walls their absorption: a rigid canyon's reflections die away as 1/n. Each step
  marches only its light cones' box, which must grow at the scheme's fastest group speed (1.23 c
  at fourth order), not at c; the GPU kernels are bound by memory, so cut cells, not arithmetic. The
  transfer divides by the march's own time difference, -i (2/dt) sin(omega dt/2), not -i omega, and
  the air's delay weighting goes on Q as on p (the transfer at a complex frequency); either undone
  fails the two-ray test at its tightened tolerance. See
  [ground.md](docs/internals/ground.md#the-full-wave-reference-srcfdtd-niceops_fdtd_field).
- **The point-source model's roofs are the floor at their own height**: a ray reflects in one
  horizontal surface at most, a roof's bounce is tested to round-off (`Buildings::on_roof`, not
  `inside()`) so it meets its edge's reflection boundary, and no ray off a roof diffracts at that
  roof's own edges, whose coefficients carry it already. Each ray's directivity is read along
  `point_ray_launch`, minus its length's gradient in the source position.
- **Buildings that touch stand against each other as one solid.** A side shared exactly drops the
  lower building's facade, roof edge and corners there; a side shared in part is split at the
  neighbour's vertices first (`split_at_junctions`), or rays diffract at a junction's corner and
  run along the seam between the two. Two level roofs across a shared side are one roof there
  (`ceded_`), or a bounce on the join counts twice.
- **`niceops --buildings` with one far building and a floor absorbing everything is the run without
  a ground model** (`buildings_far_is_open_without_ground`); keep each ray's sphere weight on the
  sphere's own spreading and absorption conventions. Its microphones are the run's observer points
  raised by their height column, as in any run without a ground model: `--buildings` turns the
  ground model off, and a database without spectra is refused, not run to `-inf`. `niceops.cpp` cannot include `buildings.hpp`
  (its `Ring` clashes with `contour.hpp`'s), so the district stays behind `UrbanDistrict`.
- **`--buildings` adds rays in energy beyond a 1° cone by default** (`src/source_coherence`;
  THEORY 6.5). Do not assume the source's phase: weighted by A-weighted power the flight data
  (harness `coherence/`) cannot tell any cone up to 2.5° from the energy sum. The cone is there
  because rays leaving together are one emission: without it the field jumps at shadow boundaries
  and a flush microphone loses its doubling. A database's record (`main_rotor_*`, `coherence_*`,
  written by panam's `directional_coherence`) switches to the measured model; the main rotor's
  delay comes from that record only, never the vehicle file. The delay rides on the rays as an
  effective length, the wrap a one-turn shift, and holds only to 4.5 blade-passage frequencies.
  `--ray_coherence coherent` keeps the old coherent sum to the bit. The record's spread fit and
  tonal fractions set `--phase_interval` under the cone as well (`cone_coherence(angle, record)`):
  the cone is the mean, the record measures the coherence the interval is about. Do not use `std::polar` with
  the sinc as its magnitude: it can be negative, and libc++ returns NaN.
- **`--buildings` samples states, not the source** (`--state_tolerance`, 0.02 dB; ground.md "In a
  run"). Each band's transfer -- its level less the direct ray's in free field at that state -- is
  evaluated at some states and interpolated; every state keeps its own sphere. The transfer is NaN
  where unknown (a band buried 120 dB or at band_level's -300 dB floor) and -infinity where no ray
  arrives; do not merge the two, and do not drop the outright evaluation of a state whose loud band
  has no known transfer (the B407 database's top bands jump by 1000 dB between states).
  `--state_tolerance 0` is the every-state path, byte for byte.
- **Shared faces are sampled** (`PointPathOptions::sample_shared_faces`): the cascade times the
  integral's ratio to it, taken two thirds of an octave apart while the transition parameter is
  under 30. An end's angle off the face's plane counts on either side: with an end in the plane
  the cascade has the wrong sign. `--exact_shared_faces` takes the integral at every wavenumber.
  The GPU samples by the same rule (`shared_face_samples`), its `face_ratios` kernel applying R.
- **Scattering walls in `niceops_point_field` reuse the urban tables' radiosity**, lit by the
  source's own specular paths to each patch; with s = 0 nothing is built and the run is byte for
  byte the specular one. The roofs scatter too (`DiffuseOptions::roofs`), with their building's
  first wall's coefficients as their reflection takes them, so a building then scatters as a whole
  where any wall is in reach; a roof patched while its reflection kept its whole rho, or the
  reverse, would make or lose energy.
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
