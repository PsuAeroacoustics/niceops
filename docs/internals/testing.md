# Tests, goldens and mutation-guarded tests

## Running the tests

Test with the asserts live, as CI does (`-DNICEOPS_ASSERTS=ON`; see [build.md](build.md#asserts)),
and run the whole suite with no `-E` exclusion:

```sh
ctest --test-dir build --output-on-failure
```

`ctest -N` lists the tests; this file gives no count, because it goes stale. Every CTest runs
with the **source root** as its working directory. Most of them run `niceops` and compare its
output with a file in `golden/` or with another run, through `golden_compare`:

- **Tolerance.** The default is 1e-9 relative, scaled by max(1, |value|). An optional third
  argument is an absolute tolerance, for runs that are meant to agree only that closely; a value
  then agrees within the larger of the two. Header lines (anything non-numeric) are skipped.
- **Non-finite values.** A NaN agrees only with a NaN, and an infinity only with the same
  infinity (a culled cell is −inf in both runs). `golden_compare` is tested itself
  (`golden_compare_rejects_nan`, `golden_compare_rejects_inf`,
  `golden_compare_accepts_matching_inf`), because an earlier form let a footprint of NaN pass
  every comparison.
- **"Differs" tests** (a run that must move the footprint) pass on `golden_compare`'s
  `mismatch at value` message through `PASS_REGULAR_EXPRESSION`, not on its exit status, so a
  usage error or a missing file cannot satisfy them. Refusals and warnings are tested the same
  way, on their message text.
- **Ordering.** A comparison reads files its producing run writes, so each pair is tied by a CTest
  fixture (`FIXTURES_SETUP`/`FIXTURES_REQUIRED`); `ctest -j` and `ctest -R` then run the
  producer first.

The Catch2 binary is `unit_test` under CTest; run it directly by tag for speed
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
| `197421_LAE.csv` | C++ build | `regression_golden` (1e-9 rel.); also `blend_barycentric_is_golden`, `ambient_negligible_matches` | `regression_default` | copy `integrated/197421_LAE.csv` |
| `metrics_197421_lafmax.csv` | C++ build | `lafmax_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_lafmax.csv` |
| `metrics_197421_na64.csv` | C++ build | `nax_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_na64.csv` |
| `metrics_197421_ta70.csv` | C++ build | `ta_golden` (1e-9 rel.) | `metrics_broadband` | copy `integrated/metrics_bb_ta70.csv` |
| `stratified_windy_point.csv` | C++ build | `stratified_windy_point_golden` (1e-4 dB abs.) | `stratified_windy_point` | copy the run's output |
| `stratified_duct_point.csv` | C++ build | `stratified_duct_point_golden` (1e-4 dB abs.) | `stratified_duct_point` | copy the run's output |
| `stratified_duct_plate.csv` | C++ build | `stratified_duct_plate_golden` (1e-4 dB abs.) | `stratified_duct_plate` | copy the run's output |
| `stratified_calm_floor_plate_200hz.csv` | C++ build | `stratified_calm_floor_plate_200hz_golden` (1e-4 dB abs.) | `stratified_calm_floor_plate_200hz` | copy the run's output |
| `terrain_hills_point.csv` | C++ build | `terrain_hills_point_golden` (1e-6 dB abs.) | `terrain_hills_point` | copy the run's output |
| `terrain_hills_plate.csv` | C++ build | `terrain_hills_plate_golden` (1e-6 dB abs.) | `terrain_hills_plate` | copy the run's output |

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
- **`terrain_hills_{point,plate}.csv`** were made by the build that added `--terrain`. They are
  pinned for the same reason: passing the plate's receiver as a point receiver at its call site
  (the plate then measures the source from the line, not from its microphone) moved cells by
  0.01 dB with every other terrain test green.

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
| `database_cache_loads_each_once`, `database_spellings_run`, `database_spellings_load_once` | CTest | a one-entry database cache (the front database reloaded in `-d A -d B -d A -d B`); a database erased before its last trajectory; databases compared by the path as typed. Output cannot see any of them, so the check reads the run's own log for a second "Opening netCDF file". The front database erased after its last trajectory, which leaves the provenance a dangling reference, trips an assert after the loop instead, so this and every other run aborts with asserts on |

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

## The plotting gap

**No CTest case passes `-p`/`-f`/`-F`, so none of them executes the plotting path.** That is how
`plotting.cpp` went years missing a `<pybind11/stl.h>` include: it compiled clean and threw
`cast_error` on every call. CI covers the gap with a separate step, after the tests, that runs
the demo with `-f` and `-F` and checks that `footprint.png` and `state_CT_.png` were written. It
runs after the golden comparisons because it rewrites `integrated/197421_LAE.csv`, which
`regression_golden` reads. The state figures' names are explained in
[build.md](build.md#plotting).
