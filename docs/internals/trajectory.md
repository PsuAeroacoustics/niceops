# Trajectories: geodetic input, smoothing, nose, speed reference

The formats and options are documented for users in the manual's
[Trajectories](../MANUAL.md#trajectories), [The sphere's nose](../MANUAL.md#the-spheres-nose) and
[Speed reference](../MANUAL.md#speed-reference); the equations are in [THEORY](../THEORY.md)
§2.4 and §3. This file is what to preserve in the code.

## Geodetic input

A trajectory may be given in WGS84 (latitude, longitude, altitude) instead of the 27-column AC
state format. Things worth knowing before touching it:

- **The format is decided by width, not by the header.** 27 columns is the positional AC
  format; anything else is read by the names in its header (`read_labeled_csv` in
  `src/util/csv.hpp`). That order is load-bearing: the AC format's *own* header names `lat`,
  `lon` and `alt` columns, so a reader dispatching on the header alone takes a state file for a
  WGS84 track, and its x/y origin, which the file does not record, is about 2.5 km back along
  the track. `first_geodetic_position` applies the same width test, and a test pins it.
- **The AC path is byte-identical.** `read_labeled_csv` keeps the first record as the header
  instead of discarding it, but drops exactly the rows `read_numeric_csv(path, true)` dropped,
  so the goldens are untouched. Verified.
- **An altitude spelling that names a unit is matched only under that unit**
  (`find_geodetic_columns`): `alt_ft` under `--units ft`, `alt_m` under `--units m`. Matching
  the other would read one unit as the other and fly the aircraft three times too high or too
  low; leaving it unmatched fails the file with a message naming the units. The bare spellings
  (`alt`, `altitude`, `height`) claim nothing and are always accepted.
- **Altitude is MSL and z is height above the observers.** `--ground_elevation` is the
  difference. In `tracks/S-76D/197421AC.csv` it is a constant 54.2651 ft, which is how that
  file's own alt and z columns relate.
- **The shipped AC tracks' local frame runs +x along a compass bearing of 210°**: the approach
  heading, not the instantaneous heading in column 9 (`data[8]`), which wanders over 25°. With
  that bearing the projection reproduces that file's own x/y
  ([FIDELITY](../FIDELITY.md#trajectory-smoothing)), which is the check that the port of the
  transform is the one its producer used.
- **One frame per run, settled before anything loads.** Every trajectory in a run shares the
  observers, so a position in one has to mean the same place as in another. `run()` resolves it
  from `--origin` or from the first geodetic track and hands it to each `Trajectory`.

### Velocity and load factor are derived

That is the whole accuracy story. The projection is not the lossy part; differentiating a
spline twice for the load factor is ([FIDELITY](../FIDELITY.md#trajectory-smoothing)).

- **`--path_smoothing` is an rms deviation per axis, in `--units`, not an abstract weight**
  (default 10 ft). `fit_path` scales it by the sample count, because `within_tolerance` bounds
  the *summed* squared residual, so it means the same thing at any sample rate. The solver
  spends the whole budget rather than merely staying inside it.
- **Zero is wrong for both kinds of input**, which is why it is not the default. Interpolating
  exactly amplifies the coordinates' 1e-8 degree rounding by 1/step² (2500x on the 0.02 s demo
  track), so measured positions turn into rounding noise. Sparse waypoints are a statement of
  intent, so flying exactly through them manufactures maneuvering that is itself a noise source.
  10 ft is not the optimum for any one input shape, but it is the only value that is good on all
  of them at once. The measured table for the four input shapes is in
  [FIDELITY](../FIDELITY.md#trajectory-smoothing); do not re-derive it here. The run logs a
  note for `--path_smoothing 0` on a dense track, and for a track sparser than 1 s, whose
  load factor is an artifact of the interpolant.
- **Velocity columns are measurements, not something to re-derive.** When a file gives `vE`,
  `vN` and `vUp`, they are turned into the frame by `from_east_north_up`, interpolated exactly
  (a velocity tolerance of 0), and the load factor comes from their first derivative.
  `--path_smoothing` is a distance and applies to the fitted position only.
- **The `--path_smoothing 0` measurement depends on `SmoothingSpline` at λ = 0** being the
  interpolating spline with its own curvature (see [Supporting modules](#supporting-modules)).

### A dense track's ends are freed from the natural condition

(`fit_path`, `SmoothingSpline::within_tolerance_extended`.) A natural spline's second
derivative is zero at the first and last knot and still a third low four seconds in, so the load
factor, and with it C_T, is wrong exactly where a track that starts or stops mid-maneuver is
loudest. On a constant-acceleration test path the freed fit stays close where the natural one
does not ([FIDELITY](../FIDELITY.md#trajectory-smoothing)).

- For positions at 0.5 s spacing or finer with a nonzero `--path_smoothing`, the fit is made
  over the track extended by 10 s each way along a local quadratic through the last 3 s, so the
  natural condition lands on the extension. Only the part over the data is kept.
- Sparse waypoints and `--path_smoothing 0` are unchanged, bit for bit, on purpose: a quadratic
  through a few waypoints is a guess, and the exact interpolant has no bias to remove.
- The AC path does no smoothing, so the goldens did not move. The effect on the dense demo
  footprint is in [FIDELITY](../FIDELITY.md#trajectory-smoothing).

### Resampling never coarsens

The sound exposure integral is a Riemann sum over whatever step the states end up with. A
geodetic track is resampled at min(`-s`, the file's mean step), so a 50 Hz file keeps its own
sample times (a uniform file at or below `-s` is used as sampled), and only sparse or unevenly
spaced input is refitted onto the `-s` grid. `-s` has no other effect: it does not touch an AC
track.

## `niceops_procedure`

`src/procedure.{hpp,cpp}`, `tools/procedure.cpp`. Things to keep:

- **It writes velocity and heading columns**, so niceops interpolates the flown velocity and
  differentiates it for the load factor, rather than differentiating positions twice. A steady
  30° turn reads back as 1/cos 30° to 1e-6 ("niceops reads a flown turn's load factor as
  1 / cos(bank)").
- **Transitions are rate-limited** (pull, roll rate, acceleration), never stepped: a step in the
  path angle is an infinite load factor once differentiated.
- **A turn rolls out within its own segment**, measured in the commanded direction, starting the
  roll-out when the course still to go equals what the roll-out from the bank reached turns,
  (g / (V p)) (−ln cos φ). Leaving the roll-out to the next segment let a reversal count the
  previous turn's roll-out as its own (an S-turn ended 44° off), skipped a turn shorter than the
  full-bank roll-out, and left a final turn banked and short. The tests fly an S-turn, a 20° turn
  under the full-bank roll-out, a final turn and a left turn through north to within 1°.
- **At rest the path angle and bank are zero**, so a descending deceleration to rest can hover.
- **Time counts internal steps** (`t = n h`, h dividing `step_s` exactly), so the written states are
  evenly spaced and niceops keeps them rather than resampling.
- **Units in keys, a bare number refused, an unknown key refused** (a typo flew the default);
  a segment below `ground_alt`, or a climb whose altitude is already passed, is refused. niceops
  itself warns of states below its ground plane without terrain (`track_below_ground_warns`), after
  `--trim` and beyond the smoothing a WGS84 fit is allowed, since a landing dips that far.

## Trajectory class

`src/trajectory.{hpp,cpp}`: the `Trajectory` class. It parses the 27-column CSV into `State`
records (metric) and computes the non-dimensional flight parameters:

- C_T = LF·W/(ρ·V_tip²·πR²) and μ = v/V_tip;
- the flight path angles `fpa` (to the horizontal, asin(v_z/|v|)) and `fpa_e` (to the tip-path
  plane, asin(v·n̂/|v|), with n̂ the load-factor direction). `--qsam` (on by default) selects
  `fpa_e`.

`adapt()` resamples onto a uniform dt grid, smoothing via `src/smoothing_spline.hpp` and taking
velocity and load factor from the spline's own derivatives. It takes either a smoothing
parameter or a `tolerance` bounding how far the smoothed path may sit from the measured one. Its
call is commented out in `process_trajectory`, so it is dead code; `watts_trim()` is reachable
through `--trim`.

**The AC vertical-speed sign is decided from the track** (`correct_vertical_speed_sign`). The AC
format does not agree with itself: the 2019 tracks store vz positive up, as z is, while the 2017
Noise Abatement tracks store it positive down (except the EC130 files of day 298). Read as up,
every 2017 descent becomes a climb of the same angle and selects climb spheres. The reader
takes the uncentered cosine between vz and the central-difference dz/dt; when it is below −0.5,
and both have more than 0.15 m/s rms, vz is negated with a warning and the path angles are
recomputed. A track with too little vertical motion is left alone. Keep it uncentered: a steady
descent has a constant dz/dt with no variance to correlate, and must still decide. The
`[trajectory][vz]` tests pin the positive-down, positive-up, level and shipped-track cases.

## Hover nose and `--nose`

**At hover the sphere's nose is the heading** (`nose_direction`, `heading_offset` in
`src/trajectory.cpp`).

- Below 3 kt (`heading_nose_speed`) the velocity is drift, a few tenths of a knot in any
  direction, so the nose takes the heading. Between 3 and 8 kt (`velocity_nose_speed`) the nose
  turns from the heading to the velocity on a smoothstep weight, the horizontal angle the short
  way round, so there is no step. A hard switch at 5 kt stepped the level by 2.2–4.4 dB for a
  nose rotated 30–90 degrees.
- The AC format's column 9 (`State::heading`, compass; a WGS84 track's `heading` column, NaN
  without one) is turned
  into the local frame by `--frame_bearing` (or east under a settled geodetic frame) whenever
  the run knows one. Only otherwise is it turned by an offset inferred from the track's states
  over 20 kt: the circular mean of the velocity's angle plus the heading, which folds in the INS
  heading's few degrees off the ground track and the mean crab (see below for why the frame
  wins). With neither, the nose follows the drift, with a warning.
- **A WGS84 heading is unwrapped as it is read and interpolated linearly**, never by spline (a
  heading through a sharp turn would overshoot), and wrapped back into [0, 360). Unwrapping as it
  is read is what sends 350 → 10 through north. Its round trip from the AC file's column 9 is
  exact ("A WGS84 track's heading column is the nose's compass heading").
- No shipped track goes under 12 kt (AW139 `192217AC.csv` reaches 12.5 kt, column 21, `VGk`) and
  no 2017 B407 flight track other than the day-283 hovers under 13, so neither the goldens nor
  the approaches depend on this; the 2017 hover runs (H1-H12) do.

**`--nose velocity|ground_velocity|heading`** (`NoseSource`, `nose_direction`; settled per
trajectory in `settle_nose`, `niceops.cpp`) says where the sphere frame's nose comes from. The
default follows the database's root `azimuth_reference` (`track` or `heading`, absent = track,
read beside `speed_reference` by `read_root_text`). The things to preserve:

- **Heading mode's offset comes from the frame and nowhere else** (`--frame_bearing`, or 90
  under a settled geodetic frame: `known_frame_bearing`). Never from `heading_offset`'s fit:
  that is the circular mean of heading against track, so it folds the track's mean crab into the
  offset, and heading mode read through it puts the nose back on the mean track. The option
  would then do nothing on a steady crabbed run, which is exactly the case it was measured on
  ([FIDELITY](../FIDELITY.md#the-spheres-nose)). No frame is refused.
- **The frame also beats the fit for the hover nose** (`heading_offset` returns it first).
  Otherwise an approach ending in hover has its hover nose turned by the approach's mean crab.
  Only fast-then-hover tracks are affected; none is in the repository or the 2017 B407 set.
  `nose_hover_takes_the_frame` runs two synthetic ones that differ only in the approach's crab,
  and fails under the fit-first precedence.
- **Heading mode keeps the velocity's climb** (vz over horizontal speed, faded below 8 kt by
  `nose_velocity_weight`, the same smoothstep the hover blend uses), so a heading along the track
  is parallel to the velocity, and `nose_heading_on_track_is_velocity` holds at 1e-9 through
  turns and descents. A horizontal-only heading would differ wherever the tip-path plane tilts
  sideways.
- **Only the nose moves.** μ, `fpa`, `fpa_e` and C_T come from the velocity as before. The
  experiment that motivated the option turned each track's velocity onto its heading instead,
  which also moved `fpa_e` by up to 1.4 degrees in turns (the load-factor normal tilts
  sideways); that is why the option and the emulation differ slightly
  ([FIDELITY](../FIDELITY.md#the-spheres-nose)).
- **Crab sign:** compass headings turn clockwise from north and the local frame counterclockwise
  from +x, so the local angle is bearing − heading, as `CompassToLocal` turns the wind. The test
  `--nose heading points the nose along the heading, crab and all` holds both to one synthetic
  crosswind.
- **Refusals.** A track with no heading at all is refused under an explicit `--nose heading`,
  and falls back to the velocity with a warning when only the database implied it. A single
  missing heading is refused. `--export_sphere` writes `azimuth_reference` only when the source
  recorded one.

## `--speed_reference`

It takes the flight condition against the air, not the ground: the source-side half of the
wind work in the [propagation roadmap](README.md#propagation-roadmap). The things to preserve:

- **Scope.** `Trajectory::relative_to_air` sets each state's wind (`State::wx`, `wy`) and
  recomputes μ, `fpa` and `fpa_e` from `air_velocity()`. `process_trajectory` hands the same air
  velocity to the sphere frame, so the nose follows it too. LF, C_T and the load-factor normal
  stay inertial.
- **Calm is bit-identical.** `vx - 0.0` is `vx`, and `speed_reference_air_calm_is_ground` checks
  it end to end.
- **The database decides by default.** The root text attribute `speed_reference`, `ground` or
  `air` (absent is ground), is read into `Database::speed_reference()`. It is settled per
  trajectory in `refer_to_database_speed`, because each trajectory can name its own database.
  Air with no wind profile (no wind columns in `--atmosphere`, and no `--wind`) is refused, not
  taken as calm. A `--wind` of zero is the exception: calm was stated, so the ground velocity is
  the air velocity.
- **It runs after `--trim`,** whose "slowest point" is over the ground.
- **One rotation.** `CompassToLocal` (`atmosphere_profile.hpp`) turns east/north into the local
  frame for this and for both ray models. Do not spell the rotation out again in a caller.
- **panam writes the attribute** (eeb040b and since, `build_empirical_database(...,
  speed_reference=...)`, default `ground`), but its 2017 recipes build against the ground, so
  every shipped database is ground-referenced, and an air run against one mixes the two.

## Supporting modules

- **`src/geodetic.{hpp,cpp}`: the WGS84 transform.** It is the American polyconic projection
  about a reference point, not a transverse Mercator: it reproduces the shipped tracks' own x/y,
  and agrees with a Krueger transverse Mercator to 3 mm at 10 km. Altitude is carried through as
  z, so this is a map projection with a height attached rather than a true ENU frame; over the
  few kilometers a footprint spans, the difference is far below the modeling error.
  - `LocalFrame(reference)` is east/north; `track_aligned(reference, bearing)` turns +x onto a
    compass bearing.
  - `from_east_north_up` applies that rotation without the translation, for velocities.
  - `to_geodetic` inverts the projection by Newton iteration, which is what `--latlon` output
    uses.
  - Ported from the D `worldtransforms.d`/`wgs84tocartesian.d` pair, which had no notion of a
    rotated frame and whose inverse walked fixed 1e-5 degree increments.
- **`src/smoothing_spline.{hpp,cpp}`: natural cubic smoothing spline** (Reinsch's algorithm;
  the interior second derivatives solve a symmetric pentadiagonal system in O(n)). It replaced
  the vendored SPLINTER library, which brought a bundled Eigen 3.2.8 that predates C++17.
  `trajectory.cpp` uses it directly. Two ways to pick how hard to smooth:
  - `with_smoothing(lambda)` takes the penalty weight directly.
  - `within_tolerance(tol)` solves for the lambda whose weighted squared residual just meets
    `tol`, as MATLAB's `spaps` does, and returns the weighted least-squares line in closed form
    when `tol` reaches that line's residual. A search that climbs lambda toward that line loses
    positive definiteness on the way from about 30,000 samples, and would refuse a straight,
    level leg. The system's condition number goes as n⁴, so past about 40,000 points a `tol`
    just under the line's residual is unreachable, and the last solvable lambda is used, with a
    warning.
  - `within_tolerance_extended` does the same over a track padded along a local quadratic at
    each end, which is what a dense geodetic track uses (see
    [A dense track's ends](#a-dense-tracks-ends-are-freed-from-the-natural-condition)).
- **`SmoothingSpline` at lambda = 0 is the natural cubic *interpolating* spline**, not a
  polyline. The lambda = 0 limit of `PenalizedSystem::fit`'s own system is `R gamma = Q^T y`,
  which has to be solved; short-circuiting it leaves the curvature zero, the load factor exactly
  1 g on every state, and the `--path_smoothing 0` measurement wrong. The earlier short circuit
  is in [HISTORY](../HISTORY.md#changes-that-affect-results).
