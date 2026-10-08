# Refraction: `--refraction linear` and `stratified`

- `src/refraction.{hpp,cpp}` — `--refraction linear` (roadmap step 4). Per path, a
  least-squares line through the effective sound speed c(z) + w(z).h over the heights the path
  spans (`EffectiveSoundSpeed`: running integrals of c, w_x, w_y and z times each on a 1 m
  table, so a fit costs a few flops per pair), then closed-form circular-arc eigenrays. The
  things to preserve:
  - **Numerically stable forms, not textbook ones.** The circles are centered L = c0/|a| away,
    kilometers to 1e11 m, so every difference of squares of h = L + sigma z is formed from
    kappa(z) = z(2L + sigma z). Reflected legs are kappa/(sum of d), travel time is one atanh of
    the addition formula, and arc length is 2 rho asin(chord/2 rho). Written the obvious way
    they cancel to nothing as the gradient shrinks. Gradients under 1e-9 /s take the straight
    path.
  - **Spreading is r / cos(launch), exactly.** In a linear profile the tube area
    r |dr/dtheta| sin|arrival| / cos(launch) reduces to (r rho / h_s)^2.
    `tests/test_refraction.cpp` checks it against finite differences of traced rays; keep that.
  - **The reflected ray is a 1-D solve** on launch angle, along which the range is monotonic
    (checked across both gradient signs; for a receiver *above* the source under downward
    refraction it is monotonic on each side of a window |theta| < acos(h_s/h_r) whose launches
    never climb to the receiver, and the solve takes the side the range at the window's edge
    says, or reports none: a root of the clamped continuation inside the window is not a ray,
    and one past the ray's apex is a second, descending crossing, which this is not). It starts from the direct ray to the receiver's image
    below the ground and takes secant steps from there, which usually close in one or two. That
    warm start took the point-receiver refraction run from 4.7 to 3.6 s. A bracketed Illinois
    solve is the fallback. Under downward refraction at long range the ray is launched
    *upward*, over its top: a separate branch of the range formula.
  - **Ground models get an equivalent straight geometry.** A point receiver gets
    h' = L2 sin(grazing) - h_r and d' = L2 cos(grazing) from the reflected ray; a plate gets
    the direct ray's arrival angle and length. The ground models themselves are untouched.
  - **The direct ray ends at the microphone**, the same point the straight path uses (see
    Receivers in [ground.md](ground.md)), never clamped. Scattered flush mics on sloping ground sit below the
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
  - Shadow zones: the grazing ray's level, less the least attenuated creeping-wave mode,
    capped at `--shadow_floor`. `ShadowAttenuation` builds it: over rigid ground (no
    `--ground_effect`) it is `shadow_attenuation_rates` and `shadow_loss`'s old arithmetic,
    to the bit; over the ground's admittance it is one table per band in s = (pi f / |a|)^(1/3),
    of the first zero of Ai' + q Ai followed along q = -i w^2 beta s. Three things to keep:
    - **It is not monotonic in the ground.** The first-order shift -q/a'_1 lowers the rate
      for arg beta in (-90, -30) degrees, so ordinary grounds attenuate up to 12% *less*
      than rigid at small |q| before up to 2.3x more. Do not "fix" a soft ground coming out
      quieter than rigid in the shallow shadow.
    - **Bands whose Re q > 0 (arg beta < -60) take the rigid rate.** There the first zero
      turns into the surface wave and the follower loses it (measured: it runs to |t| ~ 100,
      near q^2, and throws). Only strong porosity rates and thin layers get there, at low
      frequency; the run logs the bands.
    - `tests/test_shadow_attenuation.cpp` holds the rate to the exact residue series' slope
      deep in the shadow, 0.02-0.3% on the benchmark grounds. Keep it.
  - **`ShadowTransition`: the level near a shadow boundary** (`--refraction linear` and
    `stratified`, with `--ground_effect`; point receivers and plates, which get their own
    table, below). Ray theory fails within a couple of diffraction lengths F = 2 k0 l^2 of
    the boundary; its ground term at the grazing geometry was 10-25 dB low over grass and
    ~5 dB high over rigid ground there. Per band and microphone height, a table V(xi, s) of
    the residue series relative to the ray's own spreading, xi = (r - r_b)/F. Things to keep:
    - **Built at a reference source height, 20 l** (0.1 s^2 when lower, keeping it under
      R/10). V collapses over source height (within ~1 dB, 5 l-100 l, from xi = -1.5) but
      the series itself fails for high sources at high frequency, so do not build it at
      the actual height.
    - **Each model at its own boundary**: the table places the exact field at the parabolic
      boundary of k^2-linear and the ray at the circle one; NICE-OPS uses the ray's.
    - **Ray model below xi = -2, blend to -1, V beyond**, floored at V(0) - floor and capped
      at 20 log10 2. `max_power_ratio` feeds the peak search's ground bound, or the windowed
      run stops being exact.
    - **A band whose boundary is nearer the source than two diffraction lengths keeps the
      ray model to the boundary**, blends to V over the next diffraction length, and takes V
      beyond (`cramped` in `apply`). There the source is a few l up at most, under the 5 l the
      collapse was checked to, and the table has no lit zone to hand over from. A 30-50 ft
      source under a weak upwind gradient has it at every band to 500 Hz: V put the whole
      path 2-5 dB under the PE, which matched straight rays. The blend leaves the creeping-wave
      loss out: it is valid only deep in the shadow, and near the boundary it put a 30 m
      source in the linear lapse 6 dB under the PE at 250 Hz. Ray model into the shadow too was
      worse deep in it (day-286 PE matrix, 50-160 Hz: 5.7 dB rms against 5.1), and so was the
      ray model with only V's fall past the boundary. Cost: that lapse case is 0.1 dB rms worse
      at 200-630 Hz (its old lit side was 1.3 dB low, its old shadow right); day 284's PE cases
      do not move.
    - `Database::spectral_level_along`'s `ground_weight` is how the ground term is faded out:
      a power of the power ratio, per band.
    - **Through layers, each band gets its own boundary**, r_b = r_s + w (r_l - r_s),
      w = x / (x + boundary_resolution), x the band's diffraction length at the fit's
      curvature over r_l - r_s (`band_boundary`). r_s is the layered boundary, r_l the linear
      fit's, and both ride on `DirectEigenray`. Scaled by r_s alone it fell 10-20 dB too
      fast at 125-250 Hz against the PE, because low bands do not see a log law's thin fast
      layer. The constant 0.363 was fitted on 12 log laws and held out on example_windy
      (6.8 dB edge rms, against 18.9 at r_s). Pooled over 16 log-law cases the edge error
      goes 30.7 -> 6.1 dB, and the lapse control stays at 0.4.
    - **Flush microphones take their own constant, 0.75** (`boundary_resolution_flush`,
      chosen by the table's receiver height or plate, never the observer's terrain z). Fitted
      per case against the PE, a flush microphone's best boundary is ~0.9 diffraction lengths
      at the fit's curvature past r_s on the measured day-284 profiles and ~1.2 on flush log
      laws, where 0.363 put it at 1.7-2.7. 1.2 m keeps 0.363 (0.387 refitted).
    - **Where the fit casts no shadow (r_l infinite), each band takes a multiple of r_s**:
      1.60 (f/1 kHz)^-0.08 flush, 1.90 (f/1 kHz)^-0.08 at 1.2 m (`fallback_multiple_*`).
      The law has nothing to interpolate toward there, and falling back to r_s itself was
      the cliff again. A strong inversion above the source does it: the fit over the path
      refracts downward while the shear under it still casts a shadow.
    - **The law cannot be made much better by fitting.** The per-case optimum is 2.7 dB
      pooled; every smooth law tried (constant multiples, power laws in f, offsets in
      diffraction lengths, tangent-extended profiles, regressions on f, z_s, r_l/r_s) sits
      at 6-9 over all cases. The optimum is sharp where the wall is steep. Do not spend another day on it
      without turbulence in the reference.
    - **Past 0.9 r_s the lit ray is the linear fit's** (`handover_fraction`, stratified.hpp).
      `layered_direct_ray(..., true)` hands it back (`from_linear_fit`, finished by
      `finish_ground`) whenever r > 0.9 r_s, *not only where the layered solve fails*, and
      keeps it past the fit's own boundary too, in its shadow. Short of r_s the layered ray
      focuses into a caustic (+10 dB) that the PE does not show; 0.9 was the best hand-over
      on both calibration sets (0.6 dB rms on day 284 over [0.5, 1] r_s, 1.4 at r_s). Going
      back to the layered shadow ray past r_l stepped 2-10 dB, where the fit's limiting ray
      is continuous. On a log law the dipping family's range is not monotonic
      and reaches tens of meters past its end, and those near-grazing rays have no
      reflected partner, so the ground term took the grazing geometry (-32 dB). It is
      formed in `form_ray`, so the peak search bounds the ray whose level it takes. Only
      for observers the transition serves; plates and no-ground runs are byte-identical.
    - **The calibration was confounded once, so beware.** An experiment that sets a
      boundary by hand must also set `linear_boundary` equal to it, or `apply()`'s own
      interpolation runs on top. That run falsely suggested w = 0.
    - Against the PE on a linear lapse: shadow-edge error 5-20 dB -> 0.1-1 dB.
      `tests/test_shadow_transition.cpp` holds it to the exact series, and pins the law
      against stored PE levels (example_windy, 100 m source, 500 Hz and 1 kHz).
    - **Plates** get their own table: the field at the plate's top (its thickness as the
      receiver height) plus K, the plate over bare ground, per band and azimuth, taken from
      `PlateReceiver::level` at 2 degrees. K barely depends on elevation near grazing
      (within 0.1 dB from 0.5 to 10 degrees), which is what makes one value hold through
      the shadow. The caller passes the pair's plate-frame azimuth, the same as
      `Database::spectral_level_along` forms it. `max_power_ratio` includes the loudest
      K. On the 2017 day 284 this takes `stratified`'s 5-10 degree rms from 4.5 to 3.6 dB.
    - The tables take the run's `-j` (`threads`). The harness runs a dozen `-j 1`
      processes at once, and a table building on every core in each of them took 7 s.
  - **`layered_direct_ray` solves for a microphone below the plane as if on it.** The
    layers hold their lowest speed below z = 0 and the direct-ray solve knows no ground,
    so it sent a ray grazing along under the plane: a "lit" arrival of -0.003 degrees
    whose tube was 46 times its range, 33 dB low. Ground boards 3-10 cm under the survey
    datum hit it on every run of the 2017 data. The ray still ends at the microphone,
    and uniform air never reaches the layered solve, so the straight-ray equivalence
    stands.
  - **The layered shadow's tube can be 0/0** for a receiver a few millimeters up, where the
    dipping family is a sliver that thin. It then falls back to r / cos(launch). Both
    fixes have tests in `test_shadow_transition.cpp`.
- `src/stratified.{hpp,cpp}` — `--refraction stratified` (roadmap step 6a): the effective
  sound speed piecewise linear between the profile's own levels, circular arcs within each
  layer, eigenrays solved per pair. The things to preserve:
  - **A ground-grazing boundary needs only that no ray from aloft reaches the receiver**:
    the air above the source must not outrun everything between receiver and source
    (`fastest_above(z_s) > fastest(z_r, z_s)`, the same rule families() uses for the
    upward arms). Testing against the source's own speed left the day-284 mornings -- an
    inversion over an upwind shear, 330.5 m/s aloft against 339 at the ground -- with no
    boundary, so no transition and no hand-over, and the shadow a 10-35 dB cliff the PE
    puts 30-50% further out. It moved day 284's upwind 5-10 degree bin most of all; see the
    manual (`docs/MANUAL.md`) for the table.
  - **And where the air aloft does outrun it, the boundary stands short of the skip
    distance** (`boundary_under_aloft`). Rays turned back from aloft come down only past it,
    kilometers out for an inversion a few hundred feet up; short of it the shadow is the same
    ground-grazing one. `layered_direct_ray` takes the boundary back where a ray from aloft
    or along a duct does arrive in the shadow, which leaves the duct's promotion as it was.
    Without it a 30 ft source hovering upwind on day 286 got no boundary 600 ft out: no
    transition, no hand-over, the layered shadow ray's tube focused 6.5 dB by the log law's
    lowest meters, and the whole ground term -- 12 dB over free field at 63 Hz where the PE
    gives 4. On day 286's measured frames (63-630 Hz) the states under 100 ft and 20 kt went
    from +8.9 to +4.4 dB against straight rays' +3.8, so what is left there is the source; no
    other height, speed, range or elevation cell moved by more than 0.2. On the day-286 PE matrix
    past 1000 ft it took 200-630 Hz from 9.5 to 7.5 dB rms and 800-2000 Hz from 11.7 to 9.2.
  - **Parametrize by the Snell invariant p, in four families.** Launched up or down, arriving
    on the first crossing of the receiver height or on the second, after dipping below it
    (`Arm`, `families`). The range is monotonic in p within a family. The dipping family is
    what the linear model's unique circle includes implicitly. Without it the shadow boundary
    comes out short by sqrt(2 L z_r), which was 170 m on a linear test profile.
  - **Layer crossings are closed forms written stably**
    (dx = p(c0+c1)dz/(q0+q1), plus log1p and atan2 forms for time and arc). They were checked
    against a 200,000-point quadrature to 3e-15. At a turning point the sine is pinned to 0:
    computing it as sqrt(1 - p^2 c^2) leaves ~1e-8 of rounding, which is 1e-7 in range.
  - **Flat air at the fastest height is followed only as far as the highest band resolves it**
    (`resolved_pass`, called by `families()` on the `pass` end every arrival family shares). The
    invariant 1/c_rs is horizontal where the air is fastest; with a gradient there the range of
    the rays approaching it converges, in flat air it does not. A log law calm below z0 and
    written to three decimals leaves an exactly flat floor at the ground, and on day 284 a 0.1 ft
    one put the upwind boundary 20 km out (the 1e-12 inset over a 3 cm floor, nothing physical):
    no hand-over, no transition, and flush microphones 22-28 dB low, SEL rms 1.63 -> 4.18. The
    turning region is the heights within one radian of vertical phase (k int sin dz, k at
    `Layers::resolution_frequency`, which a run sets to its highest band) of the fastest air,
    and the limiting ray's travel across it is capped at 1.5 k h^2: a wave confined to h leaves
    it within ~k h^2. Things to keep:
    - **The 1.5 is set by what it must leave alone and what it may add.** Ray theory travels
      (4/3) k h^2 across a linear gradient's region (sqrt(2 R h) against k^2 h^3 = 9R/8), less
      on exact arcs, so a cap above that never binds there and `pass` comes back bit for bit; at
      k h^2 it binds on every linear profile (`the resolution limit leaves every resolved turning
      alone` fails). What a floor can still add is the cap's slack over (4/3) k h^2, and it grows
      as k^(-1/3): at 2 k h^2, the implementer's first choice, the 0.1 ft floor moved the demo
      plates 0.26 dB, and 0.67 dB with the top band at 200 Hz; at 1.5, 0.12 and 0.19. Day 284's
      measured-profile runs move under 0.0007 dB either way. Do not raise it for margin.
    - **The highest band, not a typical one**: the geometry is one for all bands and the
      largest k removes only what no band resolves. ShadowTransition moves each band's
      boundary toward the fit's as before. A consequence worth knowing: `--band_limits` moves
      a floored profile's boundary.
    - **Do not replace it with merging thin layers.** Any thickness threshold leaves a cliff
      just above it, since a flat floor of any thickness diverges; the cap is continuous in the
      floor's thickness (0.1 ft: 186 -> 187 m; 1 m: 544 m at 10 kHz, 707 km uncapped), and the
      region, a set of heights, does not depend on how many levels describe the same air.
    - **The shortcut is exact, not a heuristic**: a single fastest height whose radian is gathered
      within the segment on each side is a linear turning, and q's trapezoid floors the phase
      (q is concave on a layer). Only past that does it integrate (closed form, `vertical_phase`)
      and solve for the edges (safeguarded Newton). Cost: nothing measurable on
      `example_windy.csv`, ~3% on `example_inversion.csv`.
    - **The region belongs to the air, not to the path** (`beyond`). Where the path ends before
      the radian, the walk carries on up through the profile and down to the ground, as far as
      the limiting ray could go (not into air faster than 1/pass), and only the travel inside the
      path is held to the cap. The first version stopped at the path's ends, which proves
      nothing: a source within a radian of the fastest height had a region as short as the path,
      and the cap bound on smooth gradients -- a 1.5 m source over a 1.2 m microphone under a
      plain lapse had its boundary brought in from 789 m to 34 m at 10 kHz, and at a 200 Hz top
      band it bound for sources to 10 m in any light wind. With the next item, 274 of 728 cases
      of a sweep of log laws differed from plain ray theory; with both fixed, none by more than
      1e-7 dB (interior maxima under a light tailwind). The ground is a hard stop.
    - **Not at the horizontal launch** (`launch_junction` in `families()`). Where c_rs is the
      source's own speed and the air above is faster, pass is where the downward family meets
      the families launched up to turn aloft, not an end of the rays. Capping it opened a gap
      between them inside the lit zone, receivers with no ray under any downward refraction.
      That, not the flat 180-190 ft of `example_inversion.csv` (a source there has faster air
      above it), is what moved the duct goldens by 9.2e-4 dB in the first version; with the
      junction left alone they are back to the bit. A flat stretch at the source's own height
      under such a profile is left to plain ray theory, as before.
    - **The goldens did not move**, but for five cells of `stratified_windy_point` by up to
      4e-6 dB: near crosswind the lowest layer of `example_windy.csv` is within 1e-4 m/s of
      flat, and its rays ran 1-40 km before the cap. The paper's measured-profile runs move under
      0.0006 dB (day 284) and 0.0002 dB (day 286).
    - Point receivers above a floor were never affected (the dipping family turns at the floor's
      top, in its gradient). A profile whose lowest level is above the ground has a floor of
      that height, constant below it; the cap follows a thick one for kilometers.
    - `stratified_calm_floor_plate_is_floorless` holds the floor to 0.3 dB of the profile
      without it (0.12 measured; 39 dB unfixed). `stratified_calm_floor_plate_200hz_golden` is
      the only run at a top band other than 10 kHz, so it is what fails if `process_trajectory`
      stops handing `LayeredSoundSpeed` its band (`a run's resolution frequency reaches every
      path's layers` covers `path_layers`).
  - **The RK4 trace in `tests/test_stratified.cpp` is fussy, for reasons worth knowing.** A
    finite-difference gradient rounds to 1e-8 /s, which is 1e-4 m over a km of ray. A step that
    straddles a level integrates across the gradient's jump and drops to first order. Even a
    step that only *lands* on a level takes its last RK4 stage from the next layer unless the
    step is held to the gradient of the layer it starts in. With all three fixed, the trace
    matches the layered model to 1e-10 m.
  - **The linear model stays alongside.** It gives the warm start (usually within a fraction
    of a degree), the fallback for a receiver at or above the source, the shadow's attenuation
    rate (its fitted gradient), and the refracted path's `profile`.
  - **Every arrival, summed in energy** (`layered_arrivals`, `Arrival` in refraction.hpp).
    The primary is the first family's root. The other families are not monotonic in general,
    so they are sampled in launch angle (32 + 8 per hop the path can hold, at most 512) and every
    crossing of r is solved. In a surface
    duct (air above the source faster than all of it below) the bounced rays come too:
    X_m(p) = X_0(p) + m L(p), from the same samples. Each earlier bounce costs the plane-wave
    |R|^2 (`add_bounce_loss`), and past every direct ray's reach the strongest is promoted to
    primary. Day 286's inversion was 8–24 dB under the PE past 4,000 ft without them, and is
    now mostly within 3 dB. Focusing is capped per ray at `max_focusing_gain_db` (+10 dB); distinct arrivals
    add in energy and can exceed it in total. All solved arrivals are retained. Things to keep:
    - **`solve` refuses a bracket that does not land** (residual over 1e-6 r). A family whose
      end ray never arrives stands in as an infinity, so it brackets falsely. In a duct that
      reported the duct's top as a lit ray, 55 dB low.
    - **`Jumps`: the range is discontinuous** wherever 1/p passes a local maximum of c that
      is the fastest air yet met above z_s, or below z_r. There the turning point leaps.
      `solve` splits its bracket there. Before that, a sign change at a jump alone ran the
      full 200 iterations to be refused, under a convergence test (`hi - lo <= 1e-17 hi`)
      that a double cannot meet. That was a third of the solves in a duct and 90% of their
      cost; the test is now that the estimate no longer falls inside the bracket.
    - **The samples are shared across families** (`Row`: down, up, dip, ground pieces at one
      p). Per family they repeated the walks 3–5 times in a duct.
    - **The sample count scales with the hops, and 32 is not generous.** A ray that bounces m
      times travels X_0(p) + m L(p), so the same interval holds m times the wiggles: aliasing
      in the hop count, not near-zero intervals (refining those helped 30%). Against a
      2048-sample reference (12 bearings x 8 source heights 5-300 m x 2 receiver heights x 120
      ranges to 10 km), pairs more than 1 dB out: inversion 11.3% at a fixed 32, 2.6% at 128,
      0.47% at 512, **2.9% at `32 + 8*m_max` (cap 512)**; windy 1.05%, 0.27%, 0.07%, **0.01%**.
      The grid is bunched toward the shallow end (t = t_a + (t_b - t_a) u^2), free, and worth
      most of the windy gain; the inversion's remainder is what is left at that cost (`+16/hop`
      gets 1.8% for another 20% of the arrival time, a cosine bunching 2.6% at +8). It costs
      4-8% of a stratified run. Against the PE (4 bands, 12 geometries, beyond 1000 ft) the
      inversion's agreement moved 0.1-0.2 dB rms for the better on average (100 ft source:
      4.1 -> 3.6, 3.8 -> 3.3, 3.4 -> 2.7, 3.6 -> 2.9; 30 ft: 5.1 -> 5.2, 3.5 -> 3.6, 2.5 -> 2.6,
      2.7 -> 2.8), the windy's not at all. `layered_arrivals`' `samples` argument exists for
      `the arrival search has converged in its sample count`. A base of 16 is far worse. Do not
      trade samples for speed.
    - **The primary's bracket is split around it, not solved whole.** Whole, it returned the
      primary or one of the caustic pair beside it, at random.
    - **Past the hand-over to the linear fit only bounced rays, and rays turned back from
      aloft, are searched** (`bounced_only`). The fit stands for every downward direct ray;
      rays from aloft exist there only under a `boundary_under_aloft`. Upwind under a lapse
      that skips the search entirely.
    - **The peak search bounds the sum** by the loudest arrival's direction over the harmonic
      sum of their squared ranges (`bound_geometry` in niceops.cpp).
      `stratified_duct_*_windowed_is_exact` holds it.
    - **Not fixed: an elevated duct whose returning rays turn just above the ground** (day
      286, run 286288, 100 ft source toward 270). It is not warm ground air: the air is
      fastest at z = 0 because the log-law headwind assumed below the 46 ft LIDAR gate takes
      0.2 m/s off in the lowest foot (the temperature rises there). Rays turned back from the
      ~330 ft inversion turn up again 1–20 ft up and reach neither a flush nor a 1.2 m
      microphone: 11–21 dB under the PE at 8,000–13,000 ft. The PE tunnels through: holding
      the wind below 46 ft flat moves it under 1 dB at 63–250 Hz and ~6 dB at 1 kHz, as a WKB
      factor exp(-2 int kappa dz) does for rays turning at 1–9 ft. The skip zone short of the
      duct's return is a second gap: with the shear removed, `stratified` sits 30–60 dB under
      the PE there and finds the return 1–2 kft late. On measured day-286 frames where
      `stratified` is >= 10 dB under straight rays (1.2% of frame-bands; 0.4–2 deg, 4–9 kft,
      sources under 300 ft) it is 10–20 dB too quiet and straight rays are within 4 dB at
      50–630 Hz; none of those frames is in an SEL window. Harness
      `propagation/micro_shadow.py`, `micro_shadow_frames.py`.
    - **Three opt-in treatments of it, compared** (`docs/MANUAL.md`, "Micro-shadows and skip zones";
      harness `micro_shadow_fixes.py`, `micro_shadow_frames.py`): `--tunneling` and
      `--skip_fill` (`layered_wave_arrivals`, stratified.cpp), and a PE `--propagation_table`
      beside the ray model. Both Airy options together take the deficit frames from -10.5 /
      -15.5 to -2.1 / -5.0 dB at 50-160 / 200-630 Hz and barely touch 800-2000 (-20.6 ->
      -18.9). The PE table puts the low sources (under 150 ft, 3-6 kft) within +0.3 to +2.5.
      The SEL moves under 0.02 dB under every option.
    - **Wave arrivals carry their loss per band (`WaveCorrection`, `add_wave_loss`)** on a
      tube taken at `wave_reference_frequency` (10 kHz), so the loss is never negative and
      the peak search's bound on the tube still bounds them. They stay out of the
      boundary's removal and the promotion to primary in `layered_direct_ray`, which only
      ray arrivals take. A tunneled one takes no ground term, plate or not.
    - **Near a tunneled family's caustic the coalescing pair is replaced by the caustic's
      band-averaged Airy field.** Left to the +10 dB focusing cap it put the return 5-12 dB
      over the PE at 250 Hz-2 kHz.
    - Cost: about +55% on a stratified run (demo grid, B407, `--ground_effect`: 23.6 -> 37 s
      on `example_windy.csv`, 57.5 -> 86 s on `example_inversion.csv`). It was +107% and +76%
      before the jumps, the shared samples and the hand-over skip.
  - **The log-law finding.** With zero wind at the ground, the upwind limiting ray is fixed by
    c(z_s)/c(0) alone. On `example_windy.csv` a 100 m source's upwind shadow starts at 516 m
    for a receiver on the ground (535 m at 1.2 m), where the linear fit puts it at 1.40 km
    (1.62 km); a test pins the first. Ray theory overstates that shadow at low frequencies:
    the PE (step 7) confirmed it, and each band's boundary now lies between the two
    (`ShadowTransition`, above).
  - Costs: 12 s on the demo grid with a point receiver (linear 3.3 s), 35 ms over a
    96-microphone array.
