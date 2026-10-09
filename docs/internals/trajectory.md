# Trajectories: geodetic input, smoothing, nose, speed reference

## Geodetic trajectory input

A trajectory may be given in WGS84 (latitude, longitude, altitude) instead of the 27-column AC
state format; `docs/MANUAL.md` documents the format and the options. Things worth knowing before
touching it:

- **The format is decided by width, not by the header.** 27 columns is the positional AC
  format; anything else is read by the names in its header (`read_labeled_csv` in
  `src/util/csv.hpp`). That order is load-bearing: the AC format's *own* header names `lat`,
  `lon` and `alt` columns, so a reader dispatching on the header alone takes a state file for a
  WGS84 track — and its x/y origin, which the file does not record, is about 2.5 km back along
  the track. `first_geodetic_position` applies the same width test, and a test pins it.
- **The AC path is byte-identical.** `read_labeled_csv` keeps the first record as the header
  instead of discarding it, but drops exactly the rows `read_numeric_csv(path, true)` dropped,
  so the goldens are untouched. Verified.
- **`alt_m` is deliberately not a recognized altitude spelling.** Altitude is read in feet;
  matching a meter spelling would read meters as feet and fly the aircraft three times too
  high, where leaving it unmatched fails the file with a message naming the units.
- **Altitude is MSL and z is height above the observers.** `--ground_elevation` is the
  difference. In `tracks/S-76D/197421AC.csv` it is a constant 54.2651 ft, which is how that
  file's own alt and z columns relate.
- **The shipped AC tracks' local frame runs +x along a compass bearing of 210°** — the approach
  heading, not the instantaneous heading in column 9 (`data[8]`), which wanders over 25°. With that bearing
  the projection reproduces that file's own x/y to under a meter over 3.2 km, which is the
  check that the port of the transform is the one its producer used.
- **One frame per run, settled before anything loads.** Every trajectory in a run shares the
  observers, so a position in one has to mean the same place as in another. `run()` resolves it
  from `--origin` or from the first geodetic track and hands it to each `Trajectory`.
- **Velocity and load factor are derived, and that is the whole accuracy story.** The
  projection is not the lossy part: altitude passes through exactly and x/y to under a meter.
  Differentiating a spline twice for the load factor is. **`--path_smoothing` defaults to 10 ft**
  and is an acceptable deviation in feet RMS per axis, not an abstract weight -- the solver
  spends the whole budget rather than merely staying inside it. Interpolating exactly
  (`--path_smoothing 0`) amplifies the coordinates' 1e-8 degree rounding by 1/step² -- 2500x on
  the 0.02 s demo track -- and moves the footprint by **8.0 dB at the median**; the 10 ft
  default takes that to 0.073 dB. Zero is wrong for both kinds of input, which is why it is no
  longer the default: measured positions turn into rounding noise, and sparse waypoints are a
  statement of intent, so flying exactly through them manufactures maneuvering that is itself a
  noise source. 10 ft is not the optimum for any one input shape (a dense 50 Hz track does
  better at 0.0656 ft, tuned to its own resolution) but it is the only value that is good on
  all of them at once.

  Those numbers were re-measured after `SmoothingSpline` was fixed to give the lambda = 0 case
  its curvature. Before that, `--path_smoothing 0` returned a spline with zero second
  derivative everywhere, so the load factor was exactly 1 g on every state, CT never moved, and
  the amplification above never happened -- the old "17% and 3 dB" measured that bug. The manual
  carries the re-measured table for all four input shapes; do not re-derive it.

  **A dense track's ends are freed from the natural condition** (`fit_path`,
  `SmoothingSpline::within_tolerance_extended`). A natural spline's second derivative is zero at
  the first and last knot and still a third low four seconds in, so the load factor, and with it
  CT, is wrong exactly where a track that starts or stops mid-maneuver is loudest. For positions
  at 0.5 s spacing or finer with a nonzero `--path_smoothing`, the fit is made over the track
  extended by 10 s each way along a local quadratic through the last 3 s, so the natural
  condition lands on the extension, and only the part over the data is kept. Dense demo
  footprint against the AC state: median 0.109 -> 0.073 dB, max 1.84 -> 1.42, cells over 1 dB
  84 -> 25. Sparse waypoints and `--path_smoothing 0` are unchanged, bit for bit, on purpose: a
  quadratic through a few waypoints is a guess, and the exact interpolant has no bias to
  remove. The AC path does no smoothing, so the goldens did not move.
- **Resampling never coarsens.** The sound exposure integral is a Riemann sum over whatever
  step the states end up with, so a 50 Hz file keeps its own sample times and only sparse or
  unevenly spaced input is refitted onto the `-s` grid.

- `src/trajectory.{hpp,cpp}` — `Trajectory` class. Parses the 27-column CSV into `State` records (metric), computing non-dimensional flight parameters: CT = LF·W/(ρ·V_tip²·πR²), μ = v/V_tip, and flight path angles `fpa` (to the horizontal, asin(v_z/|v|)) vs `fpa_e` (to the tip-path plane, asin(v·n̂/|v|) with n̂ the load-factor direction). The `qsam` CLI flag selects `fpa_e` (default on). `adapt()` resamples onto a uniform dt grid, smoothing via `src/smoothing_spline.hpp` and taking velocity and load factor from the spline's own derivatives (currently commented out in `main`; `watts_trim()` is reachable through `--trim`). It takes either a smoothing parameter or a `tolerance` bounding how far the smoothed path may sit from the measured one.
- **At hover the sphere's nose is the heading** (`nose_direction`, `heading_offset` in
  `src/trajectory.cpp`). Below 3 kt (`heading_nose_speed`) the velocity is drift, a few tenths
  of a knot in any direction, and the nose went wherever it pointed. Between 3 and 8 kt
  (`velocity_nose_speed`) the nose turns from the heading to the velocity on a smoothstep weight,
  the horizontal angle the short way round, so there is no step; the hard switch at 5 kt it
  replaced stepped the level 2.2-4.4 dB for a nose rotated 30-90 degrees. The AC format's column 9
  (`State::heading`, compass, NaN for the geodetic formats) is turned into the local frame by
  `--frame_bearing` (or east under a settled geodetic frame) whenever the run knows one, and
  only otherwise by an offset inferred from the track's states over 20 kt -- the circular mean
  of the velocity's angle plus the heading, which folds in the INS heading's few degrees off the
  ground track and the mean crab (see `--nose` below for why the frame now wins). Neither: the
  drift, with a warning.
  No shipped track goes under 10 kt and no day-286 approach under 13, so neither the goldens
  nor the approaches move; the 2017 hover runs (H1-H12) do.
- **`--nose velocity|ground_velocity|heading`** (`NoseSource`, `nose_direction`; settled per
  trajectory in `settle_nose`, niceops.cpp). Where the sphere frame's nose comes from; the
  default follows the database's root `azimuth_reference` (`track`/`heading`, absent = track,
  read beside `speed_reference` by `read_root_text`). The things to preserve:
  - **Heading mode's offset comes from the frame and nowhere else** (`--frame_bearing`, or 90
    under a settled geodetic frame: `known_frame_bearing`). Never from `heading_offset`'s fit:
    that is the circular mean of heading against track, so it folds the track's mean crab into
    the offset, and heading mode read through it puts the nose back on the mean track -- the
    option would then do nothing on a steady crabbed run, which is exactly the case it was
    measured on (day 284, 1.69 -> 1.54 dB SEL rms, `stratified_plate`). No frame is refused.
  - **The frame also beats the fit for the hover nose** (`heading_offset` returns it first). The
    fit used to win whenever a track had 50 states over 20 kt, so an approach ending in hover had
    its hover nose turned by the approach's mean crab. Only fast-then-hover tracks move; none in
    the repo or the 2017 B407 set. `nose_hover_takes_the_frame` runs two synthetic ones that
    differ only in the approach's crab (the old precedence fails it).
  - **Heading mode keeps the velocity's climb** (vz over horizontal speed, faded below 8 kt by
    `nose_velocity_weight`, the same smoothstep the hover blend uses), so a heading along the
    track is parallel to the velocity and `nose_heading_on_track_is_velocity` holds at 1e-9
    through turns and descents. A horizontal-only heading would differ wherever the tip-path
    plane tilts sideways.
  - **Only the nose moves.** mu, `fpa`, `fpa_e` and CT come from the velocity as before. The
    experiment that motivated it turned each track's velocity onto its heading instead, which
    also moved `fpa_e` by up to 1.4 degrees in turns (the load-factor normal tilts sideways);
    that is why the option and the emulation differ by up to 0.5 dB at a microphone, 0.01 in rms.
  - **Crab sign:** compass headings turn clockwise from north and the local frame counterclockwise
    from +x, so the local angle is bearing - heading, as `CompassToLocal` turns the wind. The test
    `--nose heading points the nose along the heading, crab and all` holds both to one synthetic
    crosswind.
  - A track with no heading at all is refused under an explicit `--nose heading` and falls back
    to the velocity with a warning when only the database implied it; a single missing heading is
    refused. `--export_sphere` writes `azimuth_reference` only when the source recorded one.
- **`--speed_reference`** (roadmap step 3, the NICE-OPS half). It takes the flight condition
  against the air, not the ground. The things to preserve:
  - **Scope.** `Trajectory::relative_to_air` sets each state's wind (`State::wx`, `wy`) and
    recomputes μ, `fpa` and `fpa_e` from `air_velocity()`. `process_trajectory` hands the same
    air velocity to the sphere frame, so the nose follows it too. LF, CT and the load-factor
    normal stay inertial.
  - **Calm is bit-identical.** `vx - 0.0` is `vx`, and `speed_reference_air_calm_is_ground`
    checks it end to end.
  - **The database decides by default.** The root text attribute `speed_reference`, `ground`
    or `air` (absent is ground), is read into `Database::speed_reference()`. It is settled per
    trajectory in `refer_to_database_speed`, because each trajectory can name its own
    database. Air with no wind in `--atmosphere` is refused, not taken as calm.
  - **It runs after `--trim`,** whose "slowest point" is over the ground.
  - **One rotation.** `CompassToLocal` (atmosphere_profile.hpp) turns east/north into the
    local frame for this and for both ray models. Before it, refraction.cpp and
    stratified.cpp each spelled the rotation out, with the same expression.
  - **panam writes the attribute** (eeb040b and since, `build_empirical_database(...,
    speed_reference=...)`, default `ground`), but its 2017 recipes still build against the
    ground, so every database so far is ground-referenced and an air run against one mixes the
    two.

## Supporting modules

- `src/geodetic.{hpp,cpp}` — the WGS84 transform: the American polyconic projection about a
  reference point (it is not a transverse Mercator, whatever older copies of this file say; it
  reproduces the shipped tracks' own x/y, and agrees with a Krueger transverse Mercator to 3 mm
  at 10 km), with altitude carried through as z. `LocalFrame(reference)` is east/north;
  `track_aligned(reference, bearing)` turns +x onto a compass bearing;
  `from_east_north_up` applies that rotation without the translation, for velocities; and
  `to_geodetic` inverts the projection by Newton iteration, which is what `--latlon` output
  uses. Ported from the D `worldtransforms.d`/`wgs84tocartesian.d` pair, which had no notion of
  a rotated frame and whose inverse walked fixed 1e-5 degree increments.
- `src/smoothing_spline.{hpp,cpp}` — natural cubic smoothing spline (Reinsch's algorithm; interior second derivatives solve a symmetric pentadiagonal system in O(n)). Two ways to pick how hard to smooth: `with_smoothing(lambda)` takes the penalty weight directly, `within_tolerance(tol)` solves for the lambda whose weighted squared residual just meets `tol`, as MATLAB's `spaps` does -- and returns the weighted least-squares line in closed form when `tol` reaches that line's residual (the search it replaced climbed lambda to 1e200 and, from ~30,000 samples, lost positive definiteness on the way and refused a straight, level leg; the system's condition number goes as n^4, so past ~40,000 points a `tol` just under the line's residual is unreachable and the last solvable lambda is used, with a warning). `within_tolerance_extended` does the same over a track padded along a local quadratic at each end, which is what a dense geodetic track uses (see "Velocity and load factor are derived" above). `trajectory.cpp` uses it directly. Replaced the vendored SPLINTER library, which brought a bundled Eigen 3.2.8 that predates C++17.

- **`SmoothingSpline` at lambda = 0 is the natural cubic *interpolating* spline**, not a
  polyline. `PenalizedSystem::fit` used to short-circuit there and leave the curvature zero;
  the lambda = 0 limit of its own system is `R gamma = Q^T y`, which has to be solved. See the
  `--path_smoothing` note above for what that cost.
