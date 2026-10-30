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
  build that links `src/plotting_disabled.cpp`; CI builds with plotting only, and compiles the stub
  there for `plotting_stub_is_complete`, which holds it to every `plt::` function `plotting.cpp`
  defines.
- **`state/`, `integrated/`, `log/` and `figures/` are not tracked.** Every writer must create
  its own parent directory (`ensure_parent_directory`, `src/util/paths.hpp`). A working tree is
  therefore not a faithful test of a checkout: verify in a clone (`git clone --no-hardlinks .
  /tmp/fresh`, then build and run ctest there).
- **Don't assert bit-identity against a libm call on a literal argument** in a test. GCC
  constant-folds the literal through MPFR, and the two results differ by an ULP.

## Tests

- Run the whole suite before a merge: `ctest --test-dir build --output-on-failure -j8`, with no
  `-E` or label filter. CI runs only `-LE extended`, a regression subset: `tests/extended.cmake`
  lists what it leaves out, and Catch2 cases tagged `[extended]` go with them. It also runs only
  the tests a change can reach (`.github/ci_select.py`): a doc-only push runs nothing, one tool's
  source runs that tool's tests. Tests run in parallel, so a comparison must name the fixture of
  every run it reads.
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
  broadband path, because `S-76D_load_factors.nod` has no spectra (`regression_default` also passes
  `--no-spectral --no-ground_effect` explicitly). `eval` and `eval_spectrum` share their
  cell and stencil, so any disagreement between them is physics, not interpolation.
- **GPU paths are on by default wherever there is a device, and the CPU is their reference.**
  `niceops --buildings`, `niceops_point_field`, `niceops_urban_table --wave_below` and
  `niceops_fdtd_field` take the GPU unless `--no-gpu`, and fall back to the CPU where the GPU does
  not take the case (`--edge_source`, `--fade`, `--ray_coherence measured`, `--phase_interval`, an
  FDTD grid the device cannot hold); `--gpu` asks for it and is refused where it cannot run. They
  work in float (phases in float-float), are tested against the CPU at a stated tolerance, never
  bit-for-bit, and never feed a golden: **a test pinned to CPU results passes `--no-gpu`**, or it
  fails on a Mac and passes on Linux CI. Kernels are written once, in the dialect of
  `src/gpu/kernels`, so another backend takes the same source:
  [docs/internals/gpu.md](docs/internals/gpu.md).
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

## Urban model invariants

Propagation among buildings: the plane-wave tables (`--urban_table`, `niceops_urban_table`), the
point-source model (`niceops --buildings`, `niceops_point_field`) and its full-wave reference
(`niceops_fdtd_field`). Read [urban.md](docs/internals/urban.md) before changing any of it; the GPU's share is in
[gpu.md](docs/internals/gpu.md).

Geometry:

- **Buildings that touch stand against each other as one solid.** A side shared exactly drops the
  lower building's facade, roof edge and corners there; a side shared in part is split at the
  neighbor's vertices first (`split_at_junctions`), or rays diffract at a junction's corner and
  run along the seam between the two. Two level roofs across a shared side are one roof there
  (`ceded_`), or a bounce on the join counts twice.
- **Pitched roofs are the least of planes** (`--roofs`, `Footprint::roof`;
  [urban.md](docs/internals/urban.md#pitched-roofs-footprintroof---roofs-srcbuilding_inputhpp)).
  - Convex by construction, so a valley or a step is several buildings sharing sides. A ridge needs
    its two planes least either side.
  - A wall's top is `facade_top`, never `Footprint::top`. The input's top caps the roof (a level
    plane added within 1 mm).
  - The pruning proof cuts the prism by the roof planes: the full prism hides real rays
    (`pruning finds every ray among pitched houses`).
  - A slope reflects as a mirror in its plane (`on_plane`; floor at either end, facade chains on
    both sides). Once-diffracted rays take a slope on either side of the edge: those are the
    partners that keep the field continuous where another edge cuts a slope's reflection off
    (8 dB steps without them). Twice-diffracted rays take a level roof or a slope between one bare
    end and its edge (`PathEnd::bounces`, the beam through the roof's box; one `twice` lambda over
    `TwiceEnd`s).
  - Slopes scatter (tilted patches in `add_roof_patches`), or their rays' rho (1 - s) would lose
    s; a level deck reflects and scatters (`on_roof`).
  - A shared side is judged along its profile: cut where either roof breaks or the two cross, each
    wall piece straight-topped and straight-bottomed (`Facade::bottom_start`/`bottom_end`,
    `facade_bottom`; the prefilter samples between them), with an edge wherever one roof is above
    the other.
  - The FDTD's slope is a staircase (5 dB loud at 6.3 kHz with EDGE's 3 mm cells), so pitched
    scenes are scored against the measurements and the BEM, not the FDTD.
- **The point-source model's roofs are the floor at their own height**: a ray reflects in one
  horizontal surface at most; a roof's bounce is tested to round-off (`Buildings::on_roof`, not
  `inside()`) so that it meets its edge's reflection boundary; and no ray off a roof diffracts at
  that roof's own edges, whose coefficients carry it already. Each ray's directivity is read along
  `point_ray_launch`, minus its length's gradient in the source position.
- **An end may lie on a wall or a level roof** (`Buildings::support`, within `contact`;
  [urban.md](docs/internals/urban.md#ends-on-a-surface-and-each-microphones-floor-buildingssupport-pathendsurface)).
  The floor's partners are traced; a wall's or a roof's are a factor `1 + R` in
  `reflected_amplitude`, halved where the end wedge has that surface as a face (judged
  geometrically, `edge_face`, never by an angular slack). A plate on a roof takes `P_d + R P_r`.
  An end on a slope takes the roof's coefficient; ends on edges are refused, as are plates on walls
  and slopes. In a run each microphone's floor is at its observer's z unless it stands on z = 0 or
  a surface of the district (`floor_height`); each floor is evaluated with the source lowered by
  it. The plane-wave tables take a microphone on a wall, a roof or a slope the same way
  (`BuildingReceiver::surface_factor`): the table then holds its first-order reflections (no run
  facades), a departure into the surface's own building is no ray (else a source behind it came
  through at full strength), and an end on a wedge's face is at that face's angle. The tables
  carry the roofs by default; they refuse a microphone on a building's edge, and, under
  `--no_roofs`, one on or over a roof.

Rays and the floor:

- **`niceops --buildings` with one far building and a floor absorbing everything is the run without
  a ground model** (`buildings_far_is_open_without_ground`); keep each ray's sphere weight on the
  sphere's own spreading and absorption conventions. Its microphones are the run's observer points
  raised by their height column, as in any run without a ground model: by default `--buildings`
  turns the ground model off. A database without spectra is refused, not run to `-inf`.
- **`--buildings --ground_effect` gives the floor the ground's Q, bounce by bounce**
  (`FloorModel`, `FloorBounce`; THEORY 6.5). Every place that builds a ray with a floor records
  each bounce's leg and grazing sine; a new ray kind that does not is silently rigid, and
  `each floor bounce is its leg's, unfolded` is the check. Held at the band center, Q gives the
  open-ground model exactly; the 0.06-0.14 dB `buildings_far_ground` leaves is the floor ray's
  own launch and Q across a band, not an error. Plates take P_d on the direct ray and Q P_r on
  its floor partner. The GPU's `reflected` is complex per problem (`BandProblem::slot`).
  `niceops.cpp` cannot include `buildings.hpp` (its `Ring` clashes with `contour.hpp`'s), so the
  district stays behind `UrbanDistrict`.
- **Rays reflect between two edges and cross three** (`max_between_reflections`,
  `triple_diffraction`). A reflection between edges unfolds about the facade or the floor, the
  second edge mirrored with every vector (`mirrored_edge`); under `--edge_source` it takes the
  integral unfolded about the facade (a floor bounce keeps UTD), which is what carries the vertex
  where the first edge runs into the wall. Three edges are roof edges only by default: every edge
  multiplied the demo run's cost by eight for nothing. A middle edge meeting an outer one at a
  corner is left out (a vertex, not a cascade).
- **Four roof edges and more are a chain** (`max_chain`, default 8): the corners of the taut
  string in the vertical plane through the ends, the points the convex least of the path's length.
  Only chains of more than three (fewer are the double and triple searches'); terms past the
  third go in `more_wedges`, read through `wedge_term`. On the GPU a ray's kind is its wedge count
  and the shared-face kind is `NK_KIND_FACE` (15), so a chain there has at most 14 edges. Over a
  mansard's nearly flat, face-sharing wedges the cascade stays 3-4 dB rms off the BEM, with three
  edges as with four (a gable's is 0.17).

Diffraction:

- **Shared faces are sampled** (`PointPathOptions::sample_shared_faces`): the cascade times the
  integral's ratio to it, taken two thirds of an octave apart while the transition parameter is
  under 30. An end's angle off the face's plane counts on either side: with an end in the plane
  the cascade has the wrong sign. `--exact_shared_faces` takes the integral at every wavenumber.
  The GPU samples by the same rule (`shared_face_samples`), its `face_ratios` kernel applying R.
- **`--edge_source` is the point-source model's exact diffraction, and is held to exact
  solutions**, not to the UTD rays it replaces: Sommerfeld's half-plane, the wedge series, the
  long block's boundary elements. Its second order is localized by partner corrections that
  subtract the first-order model's own ray, and a vertical edge on the floor is integrated through
  it; both are what make it consistent with the ray set, and removing either moves real cases by
  decibels. See [urban.md](docs/internals/urban.md#the-edge-source-integral-srcedge_source---edge_source); its cost
  is under "Speed" below.
- **`--tilted_step_split` is opt-in** (`PointPathOptions::tilted_step_split`, `TiltedStepSplit`):
  its shares were fitted on one step and gained less on a held-out one, so keep it off by default.
  A split ray carries its foot field's carry (`PointRay::jump_weight`, `summed_length`, which the
  phases, the sort and the air's absorption take); the foot fields do not split again, and there is
  no split under `RaySum::energy`.

The source's coherence and the states:

- **`--buildings` adds rays in energy beyond a 1° cone by default** (`src/source_coherence`;
  THEORY 6.5). Do not assume the source's phase: weighted by A-weighted power the flight data
  (harness `coherence/`) cannot tell any cone up to 2.5° from the energy sum. The cone is there
  because rays leaving together are one emission: without it the field jumps at shadow boundaries
  and a flush microphone loses its doubling.
  - A database's record (`main_rotor_*`, `coherence_*`, written by panam's
    `directional_coherence`) switches to the measured model; only the record selects it, never the
    vehicle file. The main rotor's delay, where the record gives one, turns at the database's rpm,
    else at the vehicle file's tip speed over radius. The delay rides on the rays as an effective
    length, the wrap as a one-turn shift, and holds only to 4.5 blade-passage frequencies.
  - The record's spread fit and tonal fractions set `--phase_interval` under the cone as well
    (`cone_coherence(angle, record)`): the cone is the mean, and the record measures the coherence
    the interval is about.
  - `--ray_coherence coherent` keeps the old coherent sum to the bit.
  - Do not use `std::polar` with the sinc as its magnitude: it can be negative, and libc++ returns
    NaN.
- **`--buildings` samples states, not the source** (`--state_tolerance`, 0.02 dB;
  [urban.md](docs/internals/urban.md#in-a-run---buildings-srcurban_point)). Each band's transfer (its level less the
  direct ray's in free field at that state) is evaluated at some states and interpolated; every
  state keeps its own sphere. The transfer is NaN where unknown (a band buried 120 dB, or at
  `band_level`'s -300 dB floor) and -infinity where no ray arrives: do not merge the two, and do
  not drop the outright evaluation of a state whose loud band has no known transfer (the B407
  database's top bands jump by 1000 dB between states). `--state_tolerance 0` is the every-state
  path, byte for byte.

Scattering:

- **Scattering walls in `niceops_point_field` and `niceops --buildings` (`--facade_scattering`)
  reuse the urban tables' radiosity**, lit by the source's own specular paths to each patch; with
  s = 0 nothing is built and the run is byte for byte the specular one. In a run each lighting path
  reads the sphere along its own launch (`DiffuseFacades::SourcePower`), which is what holds the
  run to the tool for an omnidirectional source at 1e-9 dB. The roofs scatter too
  (`DiffuseOptions::roofs`), with their building's first wall's coefficients as their reflection
  takes them, so a building then scatters as a whole where any wall is in reach; a roof patched
  while its reflection kept its whole rho, or the reverse, would make or lose energy.
- **A facade relief is a grating** (`--relief`, `--facade_relief`; `facade_relief`,
  `grating_paths`; THEORY 6.5).
  - Specular rays that diffract nowhere take its order 0 at the incidence recorded in
    `PointRay::hits`; rays that diffract reflect whole (the other way, DigiCity's mic 7 fell
    2.4 dB). Grating rays' floor variants add with their phases.
  - The default `hybrid` hands what a relief scatters onto other walls to the diffuse field
    (`first_scattering_elsewhere`), so its grating rays reach the microphone straight or off the
    floor only (`max_microphone_reflections` 0; more counts that energy twice).
  - The source tree goes `grating_reflections + 1` deep whatever `max_order`.
  - Shares are sub-band means; beyond the diffuse reach a relief reflects whole; a relief
    building's roof does not scatter.
- **`--scattering` with `--buildings` fills the buildings' shadows** (`turbulence_scattering`,
  Forssén et al. 2014 without its canyon term). gamma_T is taken as velocity turbulence of
  C_v^2 = (3/22) c0^2 gamma_T; the paper's temperature coefficients are not used (they contradict
  its own cross section). A pair the buildings do not shield is untouched to the bit, and the total
  is held to open ground's. Off by default: unvalidated, and the strength is the weak input.

The row-of-roofs parabolic equation:

- **A row of flat roofs with a street or a dip, crossed near grazing, is a parabolic equation**
  (`profile_pe`, `apply_profile`). Three or more roof edges in each other's transitions defeat
  every cascade (8 dB at 500 Hz over one street; up to 22 dB across a dip's roof plane).
  - `apply_profile`'s gate is a street between blocks or a dip (a touching roof lower than both
    its neighbors). Touching steps up or down keep their cascades: taking any touching row fails
    the stored step tables.
  - The row's rays give way to four, the parts of its 2-D field (source or floor image, to
    microphone or image), each at its own length so that a band averages their interference (one
    ray carrying the sum was 3.3 dB high at 2 kHz). They share one solve per wavenumber
    (`ProfileSolutions`).
  - It is weighed in by elevation, the line's clearance of the blocks' ends and the ends' angle
    off the vertical over the outer corners; the field is reciprocal and smooth through all of it
    (tests in `test_profile_pe.cpp`).
  - Split-step Padé in fixed step lengths, poles factored from the top; the domain grows as
    sqrt(8 lambda L), or the layer reflects grazing waves (3 dB at 8 kHz). Numerov, and sampling
    the parts in frequency, were tried and rejected ([urban.md](docs/internals/urban.md)).

The full-wave reference:

- **`niceops_fdtd_field` is a reference, never a prediction path.** It is held to exact
  solutions (two rays over a floor, the wedge, an impedance wall) and the GPU march to the CPU's.
  - Judge a run's length band by band (`truncation_db`), not by the tail's raw energy, and give
    parallel walls their absorption: a rigid canyon's reflections die away as 1/n.
  - Each step marches only its light cones' box, which must grow at the scheme's fastest group
    speed (1.23 c at fourth order), not at c. The GPU kernels are bound by memory, so cut cells,
    not arithmetic.
  - The transfer divides by the march's own time difference, -i (2/dt) sin(omega dt/2), not
    -i omega, and the air's delay weighting goes on Q as on p (the transfer at a complex
    frequency); either undone fails the two-ray test at its tightened tolerance.
  - See [urban.md](docs/internals/urban.md#the-full-wave-reference-srcfdtd-niceops_fdtd_field).

Speed:

- **The edge source's cost is its level curves', not its wavenumbers'** (`--edge_source`).
  Keep the tolerances 1e-3/1e-4 (`edge_source.cpp`), the cascade times `end_ratios` out of both
  transitions, and `PointPathOptions::edge_source_threads`, which point_field sets. The choice is
  per wavenumber (a band must not depend on the run's other bands); the ratios let the other
  edge's point slide (A - C^2/B), and a corner (both edges' ends near) or a partner's boundary
  takes the integral. Giving shared faces their uniform form was measured and rejected (0.74 dB off
  the BEM). It runs in a trajectory (`--buildings --edge_source`, CPU only) but is not the
  default, for its cost (about 0.8 s of CPU a pair on the demo district).

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
- `using namespace niceops;` at the top of a `src/*.cpp` is intentional (`// NOLINT`); the other
  files wrap their definitions in `namespace niceops { ... }`.
- Logging goes through `src/util/logger.hpp` (`log/niceops.log`). `fatal()` is `[[noreturn]]`
  and aborts, so RAII cleanup doesn't run on fatal paths.
- Flags are defined in the CLI11 setup in `src/niceops.cpp`; document a new one in the MANUAL.
- **Licensing:** the project is MIT. `THIRD_PARTY_NOTICES.md` holds verbatim notices for the
  fetched dependencies; regenerate it from `build/_deps` when one changes. The `Agg` backend in
  `plotting.cpp` is load-bearing: the Qt backends would pull in GPL PyQt.
