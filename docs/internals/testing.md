# Tests, goldens and mutation-guarded tests

## Running the tests

Test with the asserts live, as CI does (`-DNICEOPS_ASSERTS=ON`; see [build.md](build.md#asserts)),
and run the whole suite with no `-E` exclusion:

```sh
ctest --test-dir build --output-on-failure -j8
```

CI runs it in parallel (`-j$(nproc)`); the fixtures below are what make that safe, so a new
comparison must name the fixture of every run it reads.

`ctest -N` lists the tests; this file gives no count, because it goes stale. Every CTest runs
with the **source root** as its working directory. Most of them run `niceops` and compare its
output with a file in `golden/` or with another run, through `golden_compare`:

- **Tolerance.** The default is 1e-9 relative, scaled by max(1, |value|). An optional third
  argument is an absolute tolerance, for runs that are meant to agree only that closely; a value
  then agrees within the larger of the two. A fourth and fifth argument compare the candidate
  less an offset, after a number of leading numeric lines compared as they are (a grid file's x
  and y axes): a level that should sit a known number of dB from another. Header lines (anything
  non-numeric) are skipped.
- **Non-finite values.** A NaN agrees only with a NaN, and an infinity only with the same
  infinity (a culled cell is −inf in both runs). `golden_compare` is tested itself
  (`golden_compare_rejects_nan`, `golden_compare_rejects_inf`,
  `golden_compare_accepts_matching_inf`), because an earlier form let a footprint of NaN pass
  every comparison.
- **"Differs" tests** (a run that must move the footprint) pass on `golden_compare`'s
  `mismatch at value` message through `PASS_REGULAR_EXPRESSION`, not on its exit status, so a
  usage error or a missing file cannot satisfy them. Refusals and warnings are tested the same
  way, on their message text.
- **Runs matched on their log** go through `tests/run_logged.cmake` (`add_logged_run`). CTest
  ignores the exit status of a test with `PASS_REGULAR_EXPRESSION`, so `run_logged.cmake` prints
  `run_logged: the command exited with` on a nonzero exit, and every such test carries that as a
  `FAIL_REGULAR_EXPRESSION`. A test that sets its own `FAIL_REGULAR_EXPRESSION` must keep it in
  the list (`LOGGED_RUN_FAILED`).
- **Ordering.** A comparison reads files its producing run writes, so each pair is tied by a CTest
  fixture (`FIXTURES_SETUP`/`FIXTURES_REQUIRED`); `ctest -j` and `ctest -R` then run the
  producer first.

The Catch2 binary runs under CTest in four shards (`unit_test_0` to `unit_test_3`, Catch2's
`--shard-count`), so that `ctest -j` runs them side by side; each is its own process, so no case
may rely on another's leftovers. Run it directly by tag for speed
(`build/niceops_unit_tests "[sphere_database]"`). Hidden Catch cases (tags starting with `.`,
such as `[.bump_table]` and `[.terrain_files]`) write input files for CTests and run only when
a CTest names them.

`tests/test_fidelity_tools.py` (`python -B tests/test_fidelity_tools.py`, needs NumPy) checks the
experimental source-workflow scripts in `tools/`: emission-clock alignment, guarded source
patching, separation of heading offsets, absolute held-out scoring, correction extent,
condition-bin bounds, frame alignment and common prediction coverage. It is not a CTest.

## Golden files

| File | Made by | Compared by (tolerance) | Producing run | Regenerate |
| --- | --- | --- | --- | --- |
| `state_197421.csv` | D reference | `state_golden` (1e-9 rel.) | `regression_default` (`-S state/ci_state.csv`) | `regenerate-goldens` skill |
| `197421_LAE.csv` | C++ build | `regression_golden` (1e-9 rel.); also `blend_barycentric_is_golden`, `ambient_negligible_matches`, `mask_suspect_cells_clean_is_golden`, `config_demo_is_golden`, `operations_count_one_is_golden` | `regression_default` | copy `integrated/197421_LAE.csv` |
| `metrics_197421_lafmax.csv` | C++ build | `lafmax_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_lafmax.csv` |
| `metrics_197421_na64.csv` | C++ build | `nax_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_na64.csv` |
| `metrics_197421_ta70.csv` | C++ build | `ta_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_ta70.csv` |
| `stratified_windy_point.csv` | C++ build | `stratified_windy_point_golden` (1e-4 dB abs.) | `stratified_windy_point` | copy the run's output |
| `stratified_duct_point.csv` | C++ build | `stratified_duct_point_golden` (1e-4 dB abs.) | `stratified_duct_point` | copy the run's output |
| `stratified_duct_plate.csv` | C++ build | `stratified_duct_plate_golden` (1e-4 dB abs.) | `stratified_duct_plate` | copy the run's output |
| `stratified_calm_floor_plate_200hz.csv` | C++ build | `stratified_calm_floor_plate_200hz_golden` (1e-4 dB abs.) | `stratified_calm_floor_plate_200hz` | copy the run's output |
| `stratified_{windy_point,duct_point,duct_plate,calm_floor_plate_200hz}_exact.csv` | C++ build | `*_exact_golden` (1e-4 dB abs.) | `*_exact` (`--ray_stride 1`) | copy the run's output |
| `terrain_hills_point.csv` | C++ build | `terrain_hills_point_golden` (1e-6 dB abs.) | `terrain_hills_point` | copy the run's output |
| `terrain_hills_plate.csv` | C++ build | `terrain_hills_plate_golden` (1e-6 dB abs.) | `terrain_hills_plate` | copy the run's output |
| `terrain_steep_table_plate.csv` | C++ build | `terrain_steep_table_plate_golden` (1e-5 dB abs.) | `terrain_steep_table_plate` | copy the run's output |
| `terrain_ridge_screened_point.csv` | C++ build | `terrain_ridge_screened_point_golden` (1e-6 dB abs.) | `terrain_ridge_screened_point` | copy the run's output |
| `terrain_ridge_screened_plate.csv` | C++ build | `terrain_ridge_screened_plate_golden` (1e-6 dB abs.) | `terrain_ridge_screened_plate` | copy the run's output |
| `terrain_ridge_kirchhoff_point.csv` | C++ build | `terrain_ridge_kirchhoff_point_golden` (1e-6 dB abs.) | `terrain_ridge_kirchhoff_point` | copy the run's output |
| `terrain_ridge_kirchhoff_plate.csv` | C++ build | `terrain_ridge_kirchhoff_plate_golden` (1e-6 dB abs.) | `terrain_ridge_kirchhoff_plate` | copy the run's output |
| `terrain_ridge_kirchhoff_noground.csv` | C++ build | `terrain_ridge_kirchhoff_noground_golden` (1e-6 dB abs.) | `terrain_ridge_kirchhoff_noground` | copy the run's output |
| `terrain_ridge_kirchhoff_table.csv` | C++ build | `terrain_ridge_kirchhoff_table_golden` (1e-5 dB abs.) | `terrain_ridge_kirchhoff_table` | copy the run's output |
| `urban_canyon.csv` | C++ build | `urban_table_canyon_golden` (1e-6 dB abs.) | `urban_table_canyon` | copy the run's output |
| `urban_plate.csv` | C++ build | `urban_table_plate_golden` (1e-6 dB abs.) | `urban_table_plate` | copy the run's output |
| `urban_wave.csv` | C++ build | `urban_table_wave_golden` (1e-5 dB abs.) | `urban_table_wave` | copy the run's output |
| `urban_buildings.csv` | C++ build | `urban_buildings_golden` (1e-6 dB abs.) | `urban_buildings` | copy the run's output |
| `urban_buildings_plate.csv` | C++ build | `urban_buildings_plate_golden` (1e-6 dB abs.) | `urban_buildings_plate` | copy the run's output |
| `urban_buildings_diffuse.csv` | C++ build | `urban_buildings_diffuse_golden` (1e-6 dB abs.) | `urban_buildings_diffuse` | copy the run's output |
| `buildings_scattering.csv` | C++ build | `buildings_scattering_golden` (1e-6 dB abs.) | `buildings_scattering_run` (`--buildings --facade_scattering 0.3 --no-gpu`) | copy `integrated/buildings_scattering.csv` |

All files are in `golden/`; each producing run writes `integrated/<name>.csv` unless the table
says otherwise. To regenerate a C++ golden, run its comparison with the build that carries the
deliberate change (`ctest --test-dir build -R '^<comparing test>$'`, named in the table, pulls
the producing run in as a fixture), check that the change is meant, copy the output over the
golden, and say so in the commit. Only `state_197421.csv` goes through the `regenerate-goldens`
skill (`.claude/skills/regenerate-goldens/SKILL.md`), which builds the D original in a scratch
worktree and covers the three traps that otherwise make the result unusable: the wrong
`niceops.d` revision, D's 6-significant-figure default output, and hardcoded `/usr/local`
paths.

Notes on each:

- **`state_197421.csv`** is captured from the D reference implementation. The C++ port matched it
  to 1e-9 (actual maximum relative difference ~1e-14, a few ULPs of LDC against clang round-off).
  It is trajectory state, which no propagation change touches, so it still stands as the
  D-derived oracle. It is written by `write_state`, which stays metric under `--units ft`: the
  file encodes internal values.
- **`197421_LAE.csv`** no longer comes from D, and regenerating it from D would only restore a
  bug. It was regenerated from the C++ build twice, each time at the commit that made a
  deliberate departure (see [HISTORY](../HISTORY.md#departures-from-d)):
  - when the 4x EAA multiplier was removed. The port's agreement with D was verified before that
    change and stands on its own; D is retired as an oracle for the propagation path.
  - when spherical spreading was capped at the sphere's surface (`Database::eval` and
    `eval_spectrum` take `max(range, R)`, and absorption starts at the surface). Only the
    observers whose closest approach is inside the sphere moved; `state_197421.csv` did not.

  `blend_barycentric_is_golden` (naming the default `--blend`) and `ambient_negligible_matches`
  (`-a -200`, which exercises the marched, transformed-level path) must reproduce it, and
  `blend_nearest_moves_footprint` and `blend_aam_moves_footprint` must differ from it.
- **`metrics_197421_{lafmax,na64,ta70}.csv`** (LAFmax, NA64 and TA70 on the demo case at
  `--gridstep 1640`) were made by the C++ build that added `--lamax_weighting`. They were checked
  against an independent brute-force fast detector built from the `-S` state and the `-i`
  histories, not against D, which had no time weighting (the agreement is in
  [FIDELITY](../FIDELITY.md#regression)). They fail if the weighting hears emission time, or if
  NAx counts the instantaneous peak (64 dBA sits between the two at a few cells). They share
  `regression_golden`'s database and track, so a
  change that moves it moves them; see [metrics.md](metrics.md).
- **The four `stratified_*` goldens** exist because the windowed-versus-marched tests compare two
  runs of one binary, so a change to the arrivals' summation, bounce loss, ground geometry or
  completeness passes them (see [Mutation-guarded tests](#mutation-guarded-tests)). The
  tolerance is absolute because the root solves are not bit-stable across compilers.
  `stratified_calm_floor_plate_200hz` ends `--band_limits` at the 200 Hz band, so the layered
  solve's resolution frequency is not the 10 kHz default.
- **The four stratified goldens were regenerated when `--ray_stride 4` became the default**; the
  files from before, every state solved, are kept as `*_exact.csv` and pinned by runs with
  `--ray_stride 1`, which reproduce them to the bit. `stratified_stride_{windy,duct}_is_close`
  holds the default within 0.005 and 0.02 dB of them.
- **`terrain_hills_{point,plate}.csv`** were made by the build that added `--terrain`. They are
  pinned for the same reason: passing the plate's receiver as a point receiver at its call site
  (the plate then measures the source from the line, not from its microphone) moved cells by
  0.01 dB with every other terrain test green.
- **`terrain_steep_table_plate.csv`** pins six plate microphones under the approach on the
  5.7 degree plane, read through the ground-model table. A plate whose source's foot lies behind
  it on the path's mean line (d' < 0) takes K at the turned separation's azimuth; taking it from
  the true separation passed every other test, because such pairs are nearly overhead and too few
  to move the grids, and moves these plates by 2e-4 to 7e-4 dB. The tolerance is looser than the
  other terrain goldens' because the table is built by the test run.
- **`terrain_ridge_screened_{point,plate}.csv`** were made by the build that added
  `--terrain_screening`, behind the ridge fixture, and pin its call site: handing the screened
  pair the whole path's ground as well (ground counted twice) fails both, and giving a plate's
  microphone side the point model in place of the plate fails the plate's. The windowed/marched
  pairs compare a run with itself and pass with either.
- **`terrain_ridge_kirchhoff_{point,plate,noground,table}.csv`** were made by the builds that
  added `--terrain_screening kirchhoff` (point, plate) and its band-loss table (the other two;
  the table moved the first two by under 10⁻⁸ dB), behind the ridge fixture, and pin the
  knife edge's four call sites: the ground model over the plane, the plate, the free field with no
  ground model, and the table (looser, as the table is built by the test run). Dropping the loss
  from the free field's or the table's call passed every other test.
- **`urban_{canyon,plate}.csv`** were made by the build that added `--urban_table`, and
  regenerated when the ideal canyon learned the floor's rays (the low-pass microphone moved from
  -3.7 to -2.1 dB against open ground; flush microphones, their own floor images, did not move),
  and again, with `urban_wave.csv`, when UTD roof-edge diffraction replaced the CNOSSOS-EU
  placeholder (the low-pass microphone to +0.4 dB, the placeholder having left shadows 2-4 dB too
  deep).
  Seven microphones in ideal canyons on the demo track, poles and plates. They pin the call site:
  an image path whose separation is not mirrored back before the plate reads its azimuth fails
  the plate's, and a reflection evaluated with the ground model whatever its floor partner fails
  both, with every windowed/marched pair green.
- **`urban_wave.csv`** was made by the build that added `--wave_below`: the same microphones with
  the flush canyons' bands below 80 Hz from the wave solver (10 points per wavelength), blended
  through 125 Hz. Its tolerance is looser, the FDTD's round-off being compiler-dependent.
- **`urban_buildings{,_plate}.csv`** were made by the build that added `--buildings`: the same
  seven microphones among eleven LoD1 buildings (`tests/data/urban_buildings.csv`), as poles and
  as plates. They pin the builder's call site: passing the tool's facades as no `run_facades`
  (the run's first-order reflections then counted twice) fails both, with every windowed/marched
  pair green.
- **`urban_buildings{,_plate,_diffuse}.csv`** were regenerated when the table's directions moved
  1e-6 rad clockwise of its azimuths, off the walls' knife edges: the demo microphones moved by at
  most 0.006 dB. `urban_buildings_geojson_is_the_csv` runs the same district read from GeoJSON
  against the CSV's run, to 1e-7.
- **`urban_buildings_diffuse.csv`** was made by the build that added `--facade_scattering`: the
  same district with every wall scattering 0.2. It pins the run's reflections at rho (1 - s): written
  at rho, they count the scattered share twice and fail it, with the windowed/marched pair green.
- **`urban_buildings{,_plate,_diffuse}.csv`** moved again when UTD's edges took their ends (up to
  0.06 dB), and when the plane-wave tables took the roofs, which keeps a corner's end wave the
  foot prefilter had dropped (0.001 dB; 0.005 for the diffuse one, with the roofs scattering).
  The moves are in [HISTORY](../HISTORY.md#changes-that-affect-results).
- **`buildings_scattering.csv`** was made by the build that added `--facade_scattering` under
  `--buildings`: the test district's seven microphones in a trajectory run, every wall scattering
  0.3, on the CPU (`--no-gpu`). Like every run over the test district, it takes the demo track at
  every 4th state (`tests/data/urban_track_every4.csv`, a quarter of the cost; the golden moved
  0.023 dB when it left the full track). It pins the run's own diffuse field, which no table carries, and
  the point model's rays under it: it was regenerated with each deliberate change to those rays
  (reflections at a twice-diffracted ray's ends and between its edges, three roof edges, a floor
  bounce between two edges, UTD's edge ends, a shared face's aperture at its limit, a flat face's
  twice-diffracted ray through its plane). Beside it, `buildings_scattering_states_sampled` holds
  the sampled states within 0.01 dB of every state, and
  `buildings_scattering_out_of_reach_is_specular` holds a district out of every microphone's
  `--diffuse_radius` to the specular run.

## `self_test`

`self_test` runs `niceops --self_test`, which validates every sphere in the demo database
against its own stored levels (1,844,672 samples), dBA **and EAA**, at 1e-9. That is exact by
construction: `eval()` returns the stored double at a sphere's own condition and a grid node.
The D original's 1% tolerance let a channel filed one grid cell off pass 60-90% of samples.
Legacy hover directivity once made it fail; the repair is in
[database.md](database.md#legacy-databases-database_version--1).

`--self_test` exits when it is done. It must not go on to run the demo case, which would rewrite
the file `regression_golden` reads under any test running beside it.

## Mutation-guarded tests

A review mutated the code one line at a time and wrote tests for what survived. These are the
ones to keep when refactoring: each fails on the mutant it was written for, and the rest of the
suite passed with that mutant in place. Catch2 cases are named by their full `TEST_CASE` title;
CTests are defined in `CMakeLists.txt`.

| Test | Where | Mutant it kills |
| --- | --- | --- |
| `an exported sphere reloads as the sphere it was, cell for cell` | `test_sphere_database.cpp` | dBA and EAA swapped in the export; the load-factor correction dropped |
| `an AAM export lays out its axes and amplitudes as the format says` | `test_sphere_database.cpp` | the AAM export's readback (written with the test above) |
| `a dBA cell far from the rest of its sphere is reported, and dropped when asked` | `test_sphere_database.cpp` | `--mask_suspect_cells` leaving `sphere_max_` at the masked spike (its `cull_range` check) |
| `an island standing in a hole is its own polygon, and the hole stays with its outer ring` | `test_contour.cpp` | placing a hole by its vertex mean |
| `a negative excess-attenuation factor turns the bounds off` | `test_sphere_database.cpp` | bounds kept under a negative factor, which breaks the windowed run's exactness |
| `unevenly spaced knots give a continuous, natural spline` | `test_smoothing_spline.cpp` | `h[j]` and `h[j+1]` transposed (every other spline test used `linspace`, where it is invisible) |
| `the fit projects the wind onto the path's own direction, both components` | `test_refraction.cpp` | the northward wind's sign flipped |
| `a short fit span is widened upward, past a kink, not downward` | `test_refraction.cpp` | the span widened downward |
| `a flush microphone sees the direct ray arrive at its own angle` | `test_refraction.cpp` | sin and cos swapped in the plate branch |
| `the plate's azimuth is the ray's bearing turned by the microphone's offset`; `outboard puts the microphone away from the track` | `test_plate_receiver.cpp` | mutations of `plate_frame_azimuth` and `outboard_azimuth`, which live in `plate_receiver.hpp`: one definition for both callers |
| `add_bounce_loss is the plane-wave reflection loss, once per bounce` | `test_ground_effect.cpp` | the bounce loss's form or count |
| `stratified_windy_point_golden`, `stratified_duct_point_golden`, `stratified_duct_plate_golden` | CTest | the extra arrivals dropped from the energy sum, and the primary's bounce loss: moved the duct footprint 7.8 dB with every stratified test green. Regenerate them only with the build that changed the arrivals on purpose |
| `scattering_windowed_is_exact` | CTest | anything under `--scattering`, which no CTest ran before it |
| `propagation_table_bump_*` | CTest; the hidden Catch case `[.bump_table]` in `test_propagation_table.cpp` writes the table | a bound that left the table's gain out (0.57 dB off the marched run in 21 of 23 cells, passing every CTest). The table is eight bearings with +25 dB beyond 1500 m and a cos(bearing) tilt, made from the ground-model table and run with `--frame_bearing 210`; its with-absorption twin must differ from it (`propagation_table_absorption_is_not_added_twice`) |
| `CompassToLocal turns east/north into the frame whose +x runs along a compass bearing` | `test_atmosphere_profile.cpp` | swapped `atan2` arguments in the table's bearing (the `east_north` inverse), which passed everything while tables had one bearing |
| `every simplex has an empty circumsphere, in two to four dimensions` | `test_delaunay.cpp` | a triangulation defect seen through the public API: whatever simplex a query lands in, no input point is inside its circumsphere |
| `in_sphere agrees with an exact integer oracle, in one to five dimensions` | `test_predicates.cpp` | a wrong in-sphere predicate. The oracle is 128-bit Bareiss on lattice points, ties included; it needs no sign convention, because the simplex's centroid is always inside |
| `a calm floor thinner than a wavelength leaves the shadow boundary where it was`; `the resolution limit leaves every resolved turning alone`; `a thick flat floor is followed as far as the highest band resolves it`; `a path shorter than the turning region keeps its rays`; `a run's resolution frequency reaches every path's layers`; with the CTests `stratified_calm_floor_plate_is_floorless` and `stratified_calm_floor_plate_200hz_golden` | `test_stratified.cpp`; CTest | `resolved_pass` removed; its cap at k h² or back at 2 k h²; its frequency fixed at 10 kHz; the turning region stopped at the path's ends; the horizontal launch capped; the band dropped in `path_layers`, the `LayeredSoundSpeed` constructor or `process_trajectory` |
| `a source near or on the ground is where it is, as a receiver there is` | `test_pe.cpp` | `place_ground_source` adjusting the series' quadratic term after the source vectors are filled from it, as the code once did: its steep margin of 0.07 dB fails at 0.08 |
| the generated table in `tests/test_pnl.cpp` (values from panam through `tools/pnl_reference.py`) | `test_pnl.cpp` | the band-24 extrapolation on a slope; step 3(b) at a zero slope ("zero or negative"); the 500 Hz range; the band-sharing five-record mean truncated at the record's end; the duration constant's −13 taken at every step or none, rather than at 0.5 s only; either half of the tone masking |
| `the duration's limits are the records closest to PNLTM - 10 (panam)`; `the band-sharing adjustment never lowers EPNL (panam)` | `test_pnl.cpp` | panam F026's and F014's cases |
| `epnl_frames_are_consistent` | CTest | the reception factor dropped from the frames' weights; the frames' band offset; the run not passing `--epnl_masking` on |
| `epnl_frames_follow_reception_time` | CTest | states filed at t − r/c, or at their emission time: the energy check is blind to *when* |
| `terrain_hills_point_golden`, `terrain_hills_plate_golden` | CTest | the plate's receiver passed as a point receiver at its call site |
| `terrain_hills_table_point_is_the_ground_model`, `terrain_hills_table_plate_is_the_plate_model` | CTest | a table read over terrain taking the source's height above z = 0 for its height above the pair's mean plane: the flat-table and windowed tests compare a run with itself or with level ground, and passed with it |
| `terrain_steep_table_plate_golden` | CTest | a table-read plate taking K's azimuth from the true separation rather than the turned one (`TableRead::plate_separation`) |
| `terrain_hills_table_misses_are_free_field` | CTest | every pair over the plane counted toward the raised-height warning, including pairs the table misses that fall back to free field with no ground model; `terrain_hills_table_warns_unscreened` only shows that the pairs the table reads are counted |
| `database_cache_loads_each_once`, `database_spellings_run`, `database_spellings_load_once` | CTest | a one-entry database cache (the front database reloaded in `-d A -d B -d A -d B`); a database erased before its last trajectory; databases compared by the path as typed. Output cannot see any of them, so the check reads the run's own log for a second "Opening netCDF file". The front database erased after its last trajectory, which leaves the provenance a dangling reference, trips an assert after the loop instead, so this and every other run aborts with asserts on |
| `operations_count_three_equals_repeated_t` | CTest | the day-metric weight folded into the exponent (`pow(10, 0.1 L + lg w)`), which a count of 1 and the 1e-9 count-10 checks both pass; NAx adding 1 rather than the weight; the weight dropped from one coverage energy |
| `operations_count_three_provenance` | CTest | the weight dropped from the provenance's per-state energy |
| `operations_offset_mirror_left`, `operations_offset_mirror_right` | CTest; the hidden Catch case `[.offset_tracks]` in `test_trajectory.cpp` writes the tracks | the offset left unconverted at the `--units` seam (500 ft read as 500 m), or not handed to the job; the normal's sign. The `[lateral_offset]` cases build `Trajectory` directly and cannot see the seam |
| `an AC track's fold test reads the path's curvature, not its recorded acceleration` | `test_trajectory.cpp` | an AC track's κ taken from its recorded acceleration (`n * LF * g`), whose single-sample spikes refuse offsets well inside the turn, or from an exact (unsmoothed) fit of x and y |
| `ground_classes_tabled_matches_exact` | CTest | point receivers keyed by height alone, so a concrete microphone reads the grass table at its height: the windowed, subsets and single-class tests all passed with it |
| `ground_classes_impedance_is_the_observers`, `ground_classes_terrain_impedance_is_the_observers` | CTest | every microphone given the first class's admittance with its own class's model, on flat ground and along a terrain path. The second needs its first concrete microphone on a pole: a flush one over the hills reads its class's table, which was built with the right admittance |
| `the screening gate is the brute-force scan of the blended profile`; `a profile exactly on the line of sight does not block`; `a blended profile touching the line of sight does not block` | `test_terrain.cpp` | `GroundPlaneCache::may_block` reading the horizon at floor(D/step) rather than at the last interior sample; ray b skipped; `>=` in the quick test; `>=` in `blocked()`'s own scan, which is reached only after the quick test passes, so the last case puts one ray above the line of sight and the blend exactly on it |
| `terrain_ridge_kirchhoff_noground_golden`, `terrain_ridge_kirchhoff_table_golden` | CTest | the knife edge's loss dropped from the free-field fallback with no ground model, or from the table's level: each windowed/marched pair compares a run with itself, and `terrain_ridge_kirchhoff_is_not_cnossos` reads the plates |
| `terrain_ridge_screened_point_golden`, `terrain_ridge_screened_plate_golden` | CTest | the screened pair's call given the whole path's ground, plate and table besides its side terms; the plate's microphone side evaluated with the point model |
| `a screened pair never exceeds the bound's screened term` | `test_screening.cpp` | F's clamp to 1 dropped; `screen_gain_bound` without G_S, or without the 1/3 (it must come within 3% of the bound, which rigid ground near grazing does). The call site's screening term removed from the bound passes every `terrain_ridge_screened_*_windowed_is_exact`, although the term raises the bound in most bands, so the run logs how many bands it raises and `terrain_ridge_screened_bound` requires at least one. Which classes the term's source side takes is guarded by reading |
| `the point receiver's ground bound over several classes is their elementwise maximum` | `test_ground_classes.cpp` | the bound over the first class only, inside `point_ground_gain_bounds`. No CTest kills it: each class's bound is within about 3 dB of the others' in the bands A-weighting keeps, and `ground_classes{,_terrain}_windowed_is_exact` passed with it for grass, concrete, snow, water and passive thin layers. The call site passing every class (`point_classes`) is therefore guarded only by reading |
| `urban_table_bump_windowed_is_exact` | CTest | the urban table's largest gain left out of `path_gains`. The canyon tables' own gains, at most +3 dB, sit inside the bound's slack, so the canyon and plate pairs pass with it; the bump table adds +25 dB at 3 to 20 degrees of elevation. The first-order images left out of `urban_bound_geometry` pass every CTest, including 150 random tall canyons near the end of the approach, with and without the ground: the table's own gain and the slack cover reflections leaving 10 to 40 degrees from the direct path on this database. It is kept for directional spheres and guarded by reading |
| `urban_table_plate_golden` | CTest | an image path's separation not mirrored back before the plate reads its azimuth |
| `urban_table_canyon_golden`, `urban_table_plate_golden` | CTest | a first-order reflection evaluated with the ground model whatever its floor partner (`first_order_floor_partner` skipped) |
| `urban_buildings_golden`, `urban_buildings_plate_golden` | CTest | the builder's table keeping the first-order reflections the run evaluates (no `run_facades`), or their floor partners; the second also fails `[building_paths]` |
| `the beams find every diffracted ray the unpruned search does` | `test_building_paths.cpp` | `clip_to_beam` with the window's far side reversed, or clipping 1 mm short of its boundaries |
| `the field is continuous across a corner's shadow and reflection boundaries` | `test_point_paths.cpp` | `term_angles` taking every pole's limit a hair above it (9.3 dB on the corner's shadow boundary, at the receivers exactly on it) |
| `away from transitions the twice-diffracted cascade is the edge-source integral` | `test_edge_source.cpp` | the cascade's spreading for edges that are not parallel taken as the parallel one (1.45 dB on the perpendicular pair) |
| `with the edge source a corner's adjacent edges give a finite field` | `test_point_paths.cpp` | `where_along` taking a point on the edge's own line (the corner two edges share) for a wave's boundary, with no least length for a partner's part: a stretch of no length whose spectrum is not a number |
| `the level-curve spectrum is the nested integral` | `test_edge_source.cpp` | the level curves' breaks without the boundary points, where a directivity term's pole meets its ridge (the perpendicular pair 1-11% off) |
| `a point of an edge proved hidden is hidden` | `test_edge_visibility.cpp` | the buildings grown by 5 cm rather than shrunk by `visible_stretch`'s margin (209,070 points seen yet proved hidden). A margin of zero passes: no sampled path grazes a corner closely enough |
| `pruning the diffracted paths loses no ray` | `test_point_paths.cpp` | each visible stretch pulled in by 0.5 m at both ends (466 rays of 628 on one case); the buildings grown by 5 cm |
| `urban_buildings_diffuse_golden` | CTest | the run's first-order reflections written at rho rather than rho (1 - s) |
| `rays that leave together deviate together` | `test_source_coherence.cpp` | the interval's groups never linked (every pair alone); the rotor's delay with its sense reversed on the rays; the one-turn shift of a pair whose azimuth difference wraps left out |
| `rays leaving apart lose their coherence as frequency times angle grows` | `test_source_coherence.cpp` | every pair's mean coherency 1 (the coherent sum) under `--ray_coherence measured` |
| `within the cone rays add with their phases and beyond it in energy` | `test_source_coherence.cpp` | the pair sum not told it is the cone (every pair coherent) |
| `a roof under a reflected leg is found whatever stands far off` | `test_point_paths.cpp` | a reflected leg's roofs looked for in the plan box of its unfolded line rather than of its points and facades, which found or missed a real bounce by the plan grid's cell size, so by buildings a kilometer off |
| `a bounce on two level roofs' join is one ray` | `test_point_paths.cpp` | the later building not ceding a level shared side (`ceded_`): both roofs took the bounce, 1.4 dB loud |
| `the ray set is reciprocal at every edge order` | `test_point_paths.cpp` | the source's image tree at the specular order rather than the edges' (`tree_order`), which only `--max_edge_reflections 0` shows |
| `the fade adds partners and nothing else` | `test_point_paths.cpp` | `leg_edges` skipping its own buildings without testing them for blocking (a leg back through an L's other wing); `push_faded` truncating to the three nearest screens before keeping the ones across |
| `the edge source is the same with the district turned` | `test_edge_source.cpp` | `onto_air` left out of `angle_about` or `edge_frame_ends`, so a corner's sample lands at 0 or 2π by rounding (up to 1 dB, rotation by rotation) |
| `a partner's part reaches a shared vertex` | `test_edge_source.cpp` | the part's end not snapped to the edge's end, 1e-20 short of the vertex: the partner's whole field added, 10-16 dB |
| `a wall against part of another's stands against it there`; `a courtyard closed at T-junctions is reached over its roofs alone` | `test_buildings.cpp`; `test_point_paths.cpp` | `split_at_junctions` removed (rays along the seam into the courtyard); a side split where a neighbor only touches it at a vertex |
| `the tabulated transition function is the Faddeeva form` | `test_wedge.cpp` | `utd_transition(NaN)` indexing the table, which reads out of bounds on x86 |
| `over a rigid floor the FDTD's field is the two rays'` | `test_fdtd.cpp` | the transfer divided by -i omega (fourth order 0.0023 dB low at 1250 Hz, the tolerance 0.0015) or by the continuous -i omega + a c; Q left unweighted by the air (0.013 dB, the tolerance 0.003) |
| `the GPU's cone and energy sums are the CPU's`; `the GPU's cone sum is not below zero` | `test_gpu.cpp` | a cone set packed as coherent; `band_pair_sum` without its clamp at zero (three rays whose cone sum is 3 - 2 - 2) |
| `buildings_measured_interval_brackets` | CTest (`tests/check_interval.cmake`) | an interval collapsed onto the level, or its ends swapped |
| `the phase interval's ends are the gamma distribution's quantiles` | `test_source_coherence.cpp` | the gamma cdf's series taken past a shape of 1e4, where its 1000 terms stop short of the mean |

**A mutation check that restores a file with `mv` leaves it older than the objects built from the
mutant, and make then considers it up to date.** `touch` it, or `cp` it back.

## Other checks, by feature

These are what the user-facing features stand on. The modules' own internals files say what to
preserve; this is what fails if you don't.

### The load-factor law below 1 g

Unit tests in `test_sphere_database.cpp` (`two-sided load-factor correction follows the law below
the reference, to a floor`; `direction bounds hold for two-sided corrections and the other
policies`): two-sided follows the law on both sides of 1 g, holds at the floor through zero and
negative load factors, reaches an exported sphere, and leaves one-sided unchanged. CTests:
`load_factor_two_sided_moves_footprint` on the demo track, whose load factor crosses 1 g;
`blend_aam_windowed_is_exact`, run two-sided, for the windowed peak search; and the refusals
`min_load_factor_needs_two_sided` and `two_sided_needs_fixed_load_factor`, on their messages.

### EPNL and band frames

- **Against panam**, at 1e-6 dB: `tests/test_pnl.cpp` embeds values generated by running panam's
  own functions (`tools/pnl_reference.py`) on 29 spectra (every noy region; tones at every
  special band, 80 Hz, 100 Hz, 500 Hz, 5 and 6.3 kHz; a band-24 tone on a slope; a slope turning
  to exactly zero; F above 20 and between 1.5 and 3; silence; −∞ bands; random spectra) and 7
  histories, each with the band-sharing check on and off (it fires in three, one of them at the
  record's end), plus panam's own test cases one for one.
- **Masking**, against panam's `masked=True` (the same script): the 29 spectra (unchanged
  wherever every band is audible); six prediction-like ones (absorbed top bands, an inaudible
  10 kHz band rising out of the plunge, a −1064 dB band, `-inf` bands, an audible band before a
  plunge, and masked bands rising along their thresholds); and a pass whose EPNL masking lowers
  by 1.8 dB.
- **The frames**: unit tests (tag `[band_frames]`) hold their energy to the history's exposure
  (1e-12) through a reception factor below and above one, split a state straddling a boundary in
  proportion, and file negative times into negative frames. A one-off independent rebuild of the
  frames is in [FIDELITY](../FIDELITY.md#lamax-ta-and-epnl-against-the-2017-measurements).
- **End to end** (CTests `epnl_*`): `check_band_frames` recomputes every EPNL and PNLTM from the
  netCDF frames and compares them with the CSV, checks that each PNLTM record's tone band is
  audible, checks the files against each other and for plausibility, and holds the frames'
  A-weighted energy to a `--full_history` run's SEL within 1e-4 dB (float storage); dropping the
  reception factor, or filing states whole under their center's frame, fails that. The energy
  check cannot see
  where in time the energy went, so the frames are also rebuilt from a `--band_history` run of
  the same case (each state spread over its own reception-time column's step) and must agree
  frame by frame. The same case with `--no-epnl_masking` must read some EPNL more than 0.1 dB
  higher (`epnl_masking_lowers_tones`). Asking for either file must leave the SEL output
  unchanged (`epnl_leaves_footprint`), and the refusals are tested by their messages.

### Day-night metrics

The unit test `a day of steady noise gives each metric's period weighting`
(`test_time_weighting.cpp`) integrates a steady level over a whole day and recovers
10 lg((15 + 9·10)/24) for DNL, 10 lg((12 + 4·10^0.5 + 8·10)/24) for Lden and
10 lg((12 + 3·3 + 9·10)/24) for CNEL, to 1e-9 dB. CTests run the demo track by day, at 22:30 and
at 23:30 local, and check which metrics coincide there (`lden_by_day_is_dnl`,
`cnel_by_day_is_dnl`, `cnel_at_22h30_is_dnl`, `lden_at_23h30_is_dnl`) and which must not
(`lden_at_22h30_is_not_dnl`, `dnl_full_history_is_not_windowed`).

### The operations file

The reader's unit tests (tag `[operations]`, `test_operations.cpp`) cover the column matching,
comments, defaults, path resolution, the wind columns, the period counts, the clock (00:00 and
23:59:59 accepted, 24:00 and `7pm` refused) and every refusal by its message;
`an operations file's local_time replaces the clock a track carries` (`test_trajectory.cpp`)
covers the clock on AC and WGS84 tracks. The CTests `operations_*` run the fixtures in
`tests/data/operations/`, which name their files relative to that folder:

- a count of 1 is `regression_default`'s output byte for byte and the golden at 1e-9
  (`operations_count_one_is_bit_identical`, `operations_count_one_is_golden`), and every run
  total of a `-t` run byte for byte (`operations_count_one_totals_are_t`);
- a count of 10 is +10 dB on DNL, Lden and CNEL at 1e-9 and leaves the SEL file and LAmax
  alone; a count of 3 is three repeated `-t` byte for byte on the day metrics, NAx, TA, LAmax and
  the coverage, and at 1e-9 on the provenance;
- `local_time` at 23:30 is +10 dB of DNL on 12:00, 12:00 is the track's own daytime clock, and a
  static hold at 12:00 is the hold dated by `--start_utc 43200`;
- a fleet of rows A, B, A on two databases is `-t A -t B -t A` byte for byte, loads each
  database once (read from the run's log, as `database_cache_loads_each_once` does), and names
  each row's SEL file after it;
- day, evening and night counts of 1, 0, 0 are the demo track flown by day byte for byte on DNL,
  Lden and CNEL (`operations_period_day_is_clocked_day`), 0, 0, 1 is it flown at 23:30 at 1e-9
  (`operations_period_night_matches_clock_*`), 2, 1, 1 is 1, 0, 0 plus 10 lg(2 + 10^0.5 + 10)
  on Lden, 10 lg 13 on DNL and 10 lg 15 on CNEL (`operations_period_lden_recombines`,
  `operations_period_dnl_counts_evening_as_day`, `operations_period_cnel_recombines`), the log
  says the track's clock is not read, and a track with no clock runs under `-D`
  (`operations_period_clockless_runs`);
- an offset of 0 is `regression_default`'s output byte for byte, and a WGS84 track's `-t` run
  byte for byte (`operations_offset_zero_is_bit_identical`,
  `operations_offset_zero_geodetic_is_bit_identical`); a straight eastbound track offset 500 ft
  left or right is the same track written 500 ft north or south, at 1e-9, on a grid asymmetric in
  y (`operations_offset_mirror_*`); a subtrack through the spectral peak search is its marched run
  (`operations_windowed_is_exact_with_offset`) and differs from the centerline;
- the manifest records the operations file;
- each refusal that needs a run (`-t` with the file, two `-d`, a row without a database, a missing
  track, shared names, a bare wind speed, an empty air cell, a mismatched offset unit, `count`
  or `local_time` beside period counts, `--laeq_file` with period counts, an offset on a static
  hold, an offset that folds the path, an offset on a track never faster than 3 kt, provenance
  over two databases, a clockless row under `-D`), by its message.

The reader's side of the offset is in the `[lateral_offset]` cases (`test_trajectory.cpp`): a
straight track moves by exactly d with v and n unchanged, a turn has radius R − d and that
radius's load factor, zero is bit for bit, slow samples hold their direction, given velocities
are dropped, an AC track keeps `recorded` and its kinematics, and the fold, hover and static
refusals. Curved tracks with low and high recorded speeds must have the same fold radius.
Stops followed by a different travel direction refuse nonzero offsets; continuing in the
same direction retains the held offset. Both cases cover exact interpolation and default
smoothing.

### C-weighted SEL

`csel_at_1khz_is_sel` limits a run to the 1 kHz band, where both weighting curves are 0 dB (to
1.5e-6), and requires the A-weighted SEL to 1e-5 dB; that exercises the window, the reception
weights and the per-state bands. `csel_is_not_sel` requires a full-band run to differ. The refusals
(`csel_needs_spectra`, `csel_refuses_no_spectral`, `csel_refuses_db_cutoff`,
`csel_refuses_ambient`) are tested on their messages.

### The sphere's nose

Unit tests in `test_trajectory.cpp` pin the crab's sign against a synthetic crosswind whose crab
is known (and against `CompassToLocal`, which turns the wind, so the two conventions are one),
the heading-along-track nose against the velocity's sphere frame through a tilted tip-path
plane, the climb's fade at hover, and the frame's precedence over the fit. CTests:

- `nose_hover_takes_the_frame`: two tracks ending in the same hover after approaches that differ
  only in their crab;
- `nose_heading_on_track_is_velocity`: the demo track with its heading replaced by its own ground
  track, where heading mode reproduces the default to 1e-9;
- `nose_heading_crabbed_differs`: the heading turned 15 degrees moves the footprint;
- `nose_database_heading_is_default` and `nose_database_velocity_differs`: the database
  attribute selecting the mode by itself;
- `nose_heading_windowed_is_exact`: the windowed peak search under heading noses;
- `nose_ground_velocity_is_not_air`: `ground_velocity` differs from the air nose in a wind;
- each refusal, by its message.

### The wind profile

- **The law, the air and the conventions**, in `tests/test_atmosphere_profile.cpp`: the profile
  passes through the observation and is zero at z₀; the temperature, pressure and molar
  concentration match their closed forms; a wind from 270° blows toward the east and equals what
  `--atmosphere` makes of the same wind written as speed and direction; the speed's units; no
  calm level above the ground for any class, roughness or anemometer height; and the written
  profile reads back to the bit (heights to an ULP of the unit conversion).
- **Through the command line**, `wind_profile_is_its_law`: a run in feet that sets every shaping
  option and the three surface overrides writes its profile, which is compared with the same
  profile computed independently from the closed forms (`tests/data/wind_profile_expected.py`),
  1e-12 absolute at every level. That catches a length left in feet, an override not applied, or
  `--profile_levels` not passed on, none of which the replay can see.
- **The replay**, `wind_is_its_written_profile`: a `--wind` run must equal the run fed its
  written profile through `--atmosphere`.
- **Discretization**, `wind_is_discretization_independent`: refining the levels from 16 a decade
  to 64 must not move a plate microphone of the Coyle array under stratified refraction by
  0.1 dB.
- **Calm air**, `wind_calm_air_is_ground`; and the refusals (`wind_speed_needs_unit`,
  `wind_with_atmosphere_refused`, `wind_options_need_wind`, and others), on their messages.
- **The profile builder** (`niceops_build_atmosphere`). `build_atmosphere_surface_is_wind`: one
  surface observation with `--wind`'s choices writes the file `wind_profile_written` writes, byte
  for byte. `build_atmosphere_is_harness_profile`: the 2017 layout (seven LIDAR gates, one a fill
  value, a string every 10 ft, a station humidity, the harness's level list) equals the
  validation harness's own `Weather.profile`, which `tests/data/atmosphere_builder_expected.py`
  runs on the same synthetic run. Regenerate that CSV from the script with the harness checked
  out, never from the tool. `build_atmosphere_units_agree`: the same build in meters.
  `build_atmosphere_stratified_windowed_is_exact`: a built profile runs under stratified
  refraction, and the peak search over it equals a march. `fit_z0_recovers_known`: exact log-law
  runs (`tests/data/atmosphere_builder/fit_runs.py`) give z0 back, with the light-wind run left
  out. `build_atmosphere_z0_from_is_harness_profile`: the harness build taking z0 from that fit's
  JSON is the harness profile, so a meter-valued `z0_m` scaled by `--units`, or a renamed key,
  fails. The refusals pass on their messages. The Catch2 cases are listed in
  [atmosphere.md](atmosphere.md#wind-profile-tests).
- **The run's warnings on a profile**: `atmosphere_warns_no_ground_level`,
  `atmosphere_warns_calm_floor`, and `atmosphere_still_does_not_warn` (a profile calm throughout,
  with `FAIL_REGULAR_EXPRESSION`); `build_atmosphere_manifest` runs a built profile without either
  warning, and `build_atmosphere_manifest_records_provenance` finds its sidecar in the manifest.

## The plotting gap

**No CTest case passes `-p`/`-f`/`-F`, so none of them executes the plotting path.** That is how
`plotting.cpp` went years missing a `<pybind11/stl.h>` include: it compiled clean and threw
`cast_error` on every call. CI covers the gap with a separate step, after the tests, that runs
the demo with `-f` and `-F` and checks that `footprint.png` and `state_CT_.png` were written. It
runs after the golden comparisons because it rewrites `integrated/197421_LAE.csv`, which
`regression_golden` reads. The state figures' names are explained in
[build.md](build.md#plotting).
