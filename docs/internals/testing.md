# Goldens and mutation-guarded tests

## Golden reference files

`golden/state_197421.csv` is captured from the **D reference implementation**; the C++ port
matched it to **1e-9** (actual max relative difference ~1e-14, a few ULPs of LDC-vs-clang
round-off). It is trajectory state, which no propagation change touches, so it still stands as
the D-derived oracle.

`golden/197421_LAE.csv` **no longer comes from D.** Removing the spurious 4x EAA multiplier
(see [database.md](database.md)) deliberately departs from the D original, so the D binary can no longer produce
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

`golden/metrics_197421_{lafmax,na64,ta70}.csv` (LAFmax, NA64 and TA70 on the demo case at
`--gridstep 1640`) were made from the C++ build that added `--lamax_weighting`, and checked against
an independent brute-force evaluation, not against D, which had no time weighting. They share
`regression_golden`'s database and track, so a change that moves it moves them; see "Output
metrics" in [metrics.md](metrics.md).

To regenerate them, follow the `regenerate-goldens` skill
(`.claude/skills/regenerate-goldens/SKILL.md`) — it builds the D original in a scratch
worktree and covers the three traps (wrong `niceops.d` revision, D's 6-significant-figure
default output, hardcoded `/usr/local` paths) that otherwise make the result unusable.

## `self_test`

- `self_test` — runs `niceops --self_test`, validating every sphere in the database against its own stored levels (1,844,672 samples). **Passes.** It was `broken_self_test`, excluded from CI, while legacy hover directivity made it fail on 7,788 samples; that is fixed (see below), so it is now a real test and CI runs the full suite with no `-E` exclusion. It checks the stored dBA **and EAA** at 1e-9 (eval() returns the stored double exactly at a sphere's own condition and a grid node -- 0 difference over 22 million samples; the D original's 1% let a channel filed one grid cell off pass 60-90% of samples), and **`--self_test` now exits when it is done** instead of going on to run the demo case, which rewrote the file `regression_golden` reads under any test running beside it.

## Mutation-guarded behavior

**Mutation-guarded behavior.** A review mutated the code one line at a time and wrote tests for what survived; these are the ones to keep when refactoring (each fails on the mutant it was written for): `an exported sphere reloads as the sphere it was, cell for cell` (swapped dBA/EAA, dropped load-factor correction), the AAM readback, `a negative excess-attenuation factor turns the bounds off` (the windowed run's exactness), `unevenly spaced knots give a continuous, natural spline` (every other spline test used `linspace`, where a transposed h[j]/h[j+1] is invisible), `the fit projects the wind onto the path's own direction` and `a short fit span is widened upward` and `a flush microphone sees the direct ray arrive at its own angle` (a flipped northward wind sign, widening downward, sin/cos swapped in the plate branch), `the plate's azimuth is the ray's bearing turned by the microphone's offset` and `outboard puts the microphone away from the track` (`plate_frame_azimuth` and `outboard_azimuth` live in `plate_receiver.hpp`, one definition for both callers), `add_bounce_loss is the plane-wave reflection loss`, the three `stratified_*_golden` CTests (windowed-versus-marched compares two runs of one binary, so dropping the extra arrivals from the sum moved the duct footprint 7.8 dB with every stratified test green; regenerate those goldens with the build that changed the arrivals on purpose), `scattering_windowed_is_exact` (no CTest ran `--scattering`), `propagation_table_bump_*` (an eight-bearing table with +25 dB beyond 1500 m and a cos(bearing) tilt, written by the hidden Catch case `[.bump_table]` from the ground-model table, run with `--frame_bearing 210`: a bound that left the table's gain out differed from the marched run by 0.57 dB in 21 of 23 cells and passed every CTest; the table-with-absorption twin must differ from it) and `CompassToLocal turns east/north ...` (the table's bearing is its `east_north` inverse, swapped atan2 arguments passed everything with only one-bearing tables), `every simplex has an empty circumsphere` (through the public API: whatever simplex a query lands in, no input point is inside its circumsphere) and `in_sphere agrees with an exact integer oracle` (128-bit Bareiss on lattice points, dimensions 1-5, ties included; the oracle needs no sign convention because the simplex's centroid is always inside), and `a calm floor thinner than a wavelength leaves the shadow boundary where it was`, `the resolution limit leaves every resolved turning alone`, `a thick flat floor is followed as far as the highest band resolves it`, `a path shorter than the turning region keeps its rays` and `a run's resolution frequency reaches every path's layers`, with the CTests `stratified_calm_floor_plate_is_floorless` and `stratified_calm_floor_plate_200hz_golden` (`resolved_pass` removed, its cap at k h^2 or back at 2 k h^2, its frequency fixed at 10 kHz, the turning region stopped at the path's ends, the horizontal launch capped, the band dropped in `path_layers`, the `LayeredSoundSpeed` constructor or `process_trajectory`). For EPNL: the generated table in `tests/test_pnl.cpp` (the band-24 extrapolation on a slope, step 3(b) at a zero slope, the band-sharing mean truncated at the record's end, the unrounded duration constant, either half of the tone masking), `epnl_frames_are_consistent` (the reception factor in the frames' weights, the frames' band offset, the run not passing `--epnl_masking` on), and `epnl_frames_follow_reception_time` (states filed at t - r/c, or at their emission time: the energy check is blind to *when*). **A mutation check that restores a file with `mv` leaves it older than the objects built from the mutant, and make then considers it up to date**: `touch` it, or `cp` it back.

## The plotting gap

**No CTest case passes `-p`/`-f`/`-F`, so none of them executes the plotting path.**
That is how `plotting.cpp` went years missing a `<pybind11/stl.h>` include: it
compiled clean and threw `cast_error` on every call. CI covers the gap with a
separate step that runs the demo for real figure output, after the golden
comparisons -- it rewrites `integrated/197421_LAE.csv`, which `regression_golden`
reads. Note the figure names come out as `state_CT_.png`: `niceops.cpp` builds them
as `base + "_CT_" + ext`, so the separator lands before the extension.
