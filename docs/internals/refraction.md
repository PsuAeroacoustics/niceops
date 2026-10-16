# Refraction: `--refraction linear` and `stratified`

The equations are in [THEORY §8](../THEORY.md); how far each model agrees with the parabolic
equation and with flight data is in [FIDELITY](../FIDELITY.md#refraction-models-against-the-pe).
Both models are steps of the [propagation roadmap](README.md#propagation-roadmap): the linear
fit first, then rays through the profile's own layers. This file is what to preserve.

## Linear model

`src/refraction.{hpp,cpp}`. Per path, a least-squares line through the effective sound speed
c(z) + w(z).h over the heights the path spans (`EffectiveSoundSpeed`: running integrals of c,
w_x, w_y and z times each on a 1 m table, so a fit costs a few flops per pair), then
closed-form circular-arc eigenrays.

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
  (checked across both gradient signs). For a receiver *above* the source under downward
  refraction the range is monotonic on each side of a window |theta| < acos(h_s/h_r) whose
  launches never climb to the receiver. The solve takes the side the range at the window's
  edge says, or reports none: a root of the clamped continuation inside the window is not a
  ray, and one past the ray's apex is a second, descending crossing, which this is not. The
  solve starts from the direct ray to the receiver's image below the ground and takes secant
  steps, which usually close in one or two; that warm start took the point-receiver refraction
  run from 4.7 to 3.6 s. A bracketed Illinois solve is the fallback. Under downward refraction
  at long range the ray is launched *upward*, over its top: a separate branch of the range
  formula.
- **Ground models get an equivalent straight geometry.** A point receiver gets
  h' = L2 sin(grazing) - h_r and d' = L2 cos(grazing) from the reflected ray; a plate gets the
  direct ray's arrival angle and length. The ground models themselves are untouched.
- **The direct ray ends at the microphone**, the same point the straight path uses (see
  Receivers in [ground.md](ground.md)), never clamped. Scattered flush microphones on sloping
  ground sit below the plane, and clamping them once made plate runs disagree with the
  straight model. Uniform air reproduces the straight run to 1e-9 for plates, point receivers
  and bare observers (`refraction_uniform_*_is_straight`). The reflected ray goes to the
  microphone's height above the plane z = 0, where the ground models put it.
- **The same bounds cover refraction.** The peak search and `--db_cutoff` take a refracted
  state's geometry as its launch direction scaled by its ray tube's range, so the existing
  `DirectionBounds` tables bound it unchanged. Absorption and the shadow's loss are losses, and
  the ground's gain is bounded independently of geometry. `direct_ray` gives that geometry and
  the travel time; `finish_ground`, with its reflected-ray solve, runs only for states that
  survive the bound. `refraction_*_windowed_is_exact` checks the windowed run against a march.
  On this database and track the bounds are loose even for straight rays (92–98% of states get
  evaluated), so windowing saves little here.
- **Every ray is checked against RK4 traces** of the ray equations: position to 1e-6 m, time
  and arc length to 1e-9, reflected grazing angle to 1e-8 rad.
- **Wind directions are compass bearings** turned by `--frame_bearing` (default +x east). The
  shipped AC tracks need 210.

## Shadow attenuation

In a shadow zone the level is the grazing ray's, less the least attenuated creeping-wave mode,
capped at `--shadow_floor`. `ShadowAttenuation` builds the rate. Over rigid ground (no
`--ground_effect`) it is `shadow_attenuation_rates` and `shadow_loss`'s original arithmetic, to
the bit. Over the ground's admittance it is one table per band in s = (pi f / |a|)^(1/3), of
the first zero of Ai' + q Ai followed along q = -i w^2 beta s.

- **It is not monotonic in the ground.** The first-order shift -q/a'_1 lowers the rate for
  arg beta in (-90, -30) degrees, so ordinary grounds attenuate up to 12% *less* than rigid at
  small |q| before up to 2.3x more. Do not "fix" a soft ground coming out quieter than rigid in
  the shallow shadow.
- **Bands whose Re q > 0 (arg beta < -60) take the rigid rate.** There the first zero turns
  into the surface wave and the follower loses it (measured: it runs to |t| ~ 100, near q^2,
  and throws). Only strong porosity rates and thin layers get there, at low frequency; the run
  logs the bands.
- `tests/test_shadow_attenuation.cpp` holds the rate to the exact residue series' slope deep in
  the shadow, 0.02–0.3% on the benchmark grounds. Keep it.

## ShadowTransition: the level near a shadow boundary

`--refraction linear` and `stratified`, with a ground model; point receivers, and plates, which
get their own table (below). Ray theory fails within a couple of diffraction lengths
F = 2 k0 l^2 of the boundary: its ground term at the grazing geometry is too low over grass and
too high over rigid ground there. Per band and microphone height, a table
V(xi, s) of the residue series relative to the ray's own spreading, xi = (r - r_b)/F. The law
is in [THEORY §8.1](../THEORY.md); the calibration against the PE is in
[FIDELITY](../FIDELITY.md#shadow-transition-calibration).

- **Built at a reference source height, 20 l** (0.1 s^2 when lower, keeping it under R/10).
  V collapses over source height (from 5 l to 100 l; how closely is in FIDELITY), but the series
  itself fails for high sources at high frequency, so do not build it at the actual height.
- **Each model at its own boundary**: the table places the exact field at the parabolic
  boundary of k^2-linear and the ray at the circle one; NICEOPS uses the ray's.
- **Ray model below xi = -2, blend to -1, V beyond**, floored at V(0) - floor and capped at
  20 log10 2. `max_power_ratio` feeds the peak search's ground bound, or the windowed run stops
  being exact.
- **V's deep lit side is rounding noise at large s.** Between xi = -2 and about -1.6, at s of a
  few hundred and up, the series' terms cancel by up to nine orders at k0 r of 10^6-10^7, so a
  mode zero's last bit moves V there by up to 0.4 dB: the same build under GCC and under clang
  differs by as much at the same 460 nodes (of 52,164), and by under 2e-6 dB everywhere from
  xi = -1 on. Those nodes enter only the blend, at a weight of 0.4 or less. The trusted test
  (`modal_sum`) bounds the terms' sizes, not their phases, so it passes them. A fix would move
  results; none has been made.
- **A cramped band keeps the ray model to the boundary.** When a band's boundary is nearer the
  source than two diffraction lengths, it keeps the ray model to the boundary, blends to V over
  the next diffraction length, and takes V beyond (`cramped` in `apply`). There the source is a
  few l up at most, under the 5 l the collapse was checked to, and the table has no lit zone to
  hand over from. A 30–50 ft source under a weak upwind gradient is cramped at every band to
  500 Hz, and V alone put the whole path under the PE, which straight rays matched. The
  blend leaves the creeping-wave loss out: it is valid only deep in the shadow. The
  alternatives (the ray model into the shadow too, or the ray model with only V's fall past the
  boundary) were both worse deep in it.
- `Database::spectral_level_along`'s `ground_weight` is how the ground term is faded out: a
  power of the power ratio, per band.
- **Through layers, each band gets its own boundary**, r_b = r_s + w (r_l - r_s),
  w = x / (x + `boundary_resolution`), x the band's diffraction length at the fit's curvature
  over r_l - r_s (`band_boundary`). r_s is the layered boundary, r_l the linear fit's, and both
  ride on `DirectEigenray`. Scaled by r_s alone the shadow fell much too fast at 125–250 Hz
  against the PE, because low bands do not see a log law's thin fast layer. The
  constant was fitted on twelve log laws and held out on `example_windy.csv`.
- **Flush microphones take their own constant** (`boundary_resolution_flush`), chosen by the
  table's receiver height or plate, never the observer's terrain z. The shared constant put a
  flush microphone's boundary about twice too far out.
- **Where the fit casts no shadow (r_l infinite), each band takes a multiple of r_s**,
  m (f/1 kHz)^-0.08, with m from `fallback_multiple_flush` or `fallback_multiple_raised`. The
  law has nothing to interpolate toward there, and falling back to r_s itself was the cliff
  again. A strong inversion above the source does it: the fit over the path refracts downward
  while the shear under it still casts a shadow.
- **The law cannot be made much better by fitting.** The best boundary per case beats every
  smooth law tried (constant multiples, power laws in f, offsets in diffraction lengths,
  tangent-extended profiles, regressions on f, z_s and r_l/r_s) by a wide margin, and the
  smooth laws all land together; the numbers are in FIDELITY. The optimum is sharp where the
  wall is steep. Do not spend another day on it without turbulence in the reference.
- **Past 0.9 r_s the lit ray is the linear fit's** (`handover_fraction`, stratified.hpp).
  `layered_direct_ray(..., true)` hands it back (`from_linear_fit`, finished by
  `finish_ground`) whenever r > 0.9 r_s, *not only where the layered solve fails*, and keeps it
  past the fit's own boundary too, in its shadow. Short of r_s the layered ray focuses into a
  caustic (+10 dB) that the PE does not show; 0.9 was the best hand-over on both calibration
  sets. Going back to the layered shadow ray past r_l stepped 2–10 dB, where the fit's
  limiting ray is continuous. On a log law the dipping family's range is not monotonic and
  reaches tens of meters past its end, and those near-grazing rays have no reflected partner,
  so the ground term took the grazing geometry (-32 dB). The hand-over is formed in `form_ray`,
  so the peak search bounds the ray whose level it takes. It applies only to observers the
  transition serves; plates and no-ground runs are byte-identical.
- **The calibration was confounded once, so beware.** An experiment that sets a boundary by
  hand must also set `linear_boundary` equal to it, or `apply()`'s own interpolation runs on
  top. That run falsely suggested w = 0.
- `tests/test_shadow_transition.cpp` holds the table to the exact series and pins the law
  against stored PE levels (`example_windy.csv`, 100 m source, 500 Hz and 1 kHz).
- **Plates** get their own table: the field at the plate's top (its thickness as the receiver
  height) plus K, the plate over bare ground, per band and azimuth, taken from
  `PlateReceiver::level` at 2 degrees. K barely depends on elevation near grazing (FIDELITY),
  which is what makes one value hold through the shadow. The
  caller passes the pair's plate-frame azimuth, the same as `Database::spectral_level_along`
  forms it. `max_power_ratio` includes the loudest K.
- The tables take the run's `-j` (`threads`). Validation runs a dozen `-j 1` processes at
  once, and a table building on every core in each of them took 7 s.
- **V is built when first wanted, not at setup** (`ShadowTransition::tables`, a
  `std::once_flag`): the first ray that comes near a boundary, the peak search's bound, or a
  test. A run that marches every state (any history output, so every validation run) and whose
  rays never near a boundary builds none of it. The values are the same computation whenever
  it runs; the tables live in a `shared_ptr` so that the vector of transitions can move them.
  `apply`'s clamp reads the largest power ratio, so V is built whole, never by the band.

## Stratified rays

`src/stratified.{hpp,cpp}`, `--refraction stratified`: the effective sound speed piecewise
linear between the profile's own levels, circular arcs within each layer, eigenrays solved per
pair.

- **Parametrize by the Snell invariant p, in four families.** Launched up or down, arriving on
  the first crossing of the receiver height or on the second, after dipping below it (`Arm`,
  `families`). The range is monotonic in p within a family. The dipping family is what the
  linear model's unique circle includes implicitly. Without it the shadow boundary comes out
  short by sqrt(2 L z_r), which was 170 m on a linear test profile.
- **Layer crossings are closed forms written stably** (dx = p(c0+c1)dz/(q0+q1), plus log1p and
  atan2 forms for time and arc). They were checked against a 200,000-point quadrature to
  3e-15. At a turning point the sine is pinned to 0: computing it as sqrt(1 - p^2 c^2) leaves
  about 1e-8 of rounding, which is 1e-7 in range.
- **A ground-grazing boundary needs only that no ray from aloft reaches the receiver**: the air
  above the source must not outrun everything between receiver and source
  (`fastest_above(z_s) > fastest(z_r, z_s)`, the same rule `families()` uses for the upward
  arms). Testing against the source's own speed left the day-284 mornings (an inversion over an
  upwind shear, 330.5 m/s aloft against 339 at the ground) with no boundary, so no transition
  and no hand-over, and the shadow a 10–35 dB cliff the PE puts 30–50% further out.
- **And where the air aloft does outrun it, the boundary stands short of the skip distance**
  (`boundary_under_aloft`). Rays turned back from aloft come down only past it, kilometers out
  for an inversion a few hundred feet up; short of it the shadow is the same ground-grazing
  one. `layered_direct_ray` takes the boundary back where a ray from aloft or along a duct does
  arrive in the shadow, which leaves the duct's promotion as it was. Without it a 30 ft source
  hovering upwind on day 286 got no boundary 600 ft out: the layered shadow ray's tube focused
  in the log law's lowest meters and took the whole ground term on top, far over the PE at
  63 Hz. The effect on the measured frames is in
  [FIDELITY](../FIDELITY.md#refraction-against-flight-data).
- **`layered_direct_ray` solves for a microphone below the plane as if on it.** The layers
  hold their lowest speed below z = 0 and the direct-ray solve knows no ground, so it sent a
  ray grazing along under the plane: a "lit" arrival of -0.003 degrees whose tube was 46 times
  its range, 33 dB low. Ground boards 3–10 cm under the survey datum hit it on every run of the
  2017 data. The ray still ends at the microphone, and uniform air never reaches the layered
  solve, so the straight-ray equivalence stands.
- **The layered shadow's tube can be 0/0** for a receiver a few millimeters up, where the
  dipping family is a sliver that thin. It then falls back to r / cos(launch). Both this and
  the previous fix have tests in `test_shadow_transition.cpp`.
- **The RK4 trace in `tests/test_stratified.cpp` is fussy, for reasons worth knowing.** A
  finite-difference gradient rounds to 1e-8 /s, which is 1e-4 m over a km of ray. A step that
  straddles a level integrates across the gradient's jump and drops to first order. Even a step
  that only *lands* on a level takes its last RK4 stage from the next layer unless the step is
  held to the gradient of the layer it starts in. With all three fixed, the trace matches the
  layered model to 1e-10 m.
- **The linear model stays alongside.** It gives the warm start (usually within a fraction of
  a degree), the fallback for a receiver at or above the source, the shadow's attenuation rate
  (its fitted gradient), and the refracted path's `profile`.
- **The log-law finding.** With zero wind at the ground, the upwind limiting ray is fixed by
  c(z_s)/c(0) alone, so the layered boundary on a log law lies far short of the linear fit's; a
  test pins one such boundary on `example_windy.csv`. Ray theory overstates that shadow at low
  frequencies, as the PE confirmed, and each band's boundary lies between the two
  (`ShadowTransition`, above).

## Flat air at the fastest height: `resolved_pass`

Flat air at the fastest height is followed only as far as the highest band resolves it
(`resolved_pass`, called by `families()` on the `pass` end every arrival family shares). The
invariant 1/c_rs is horizontal where the air is fastest; with a gradient there the range of the
rays approaching it converges, and in flat air it does not. A log law calm below z0 and written
to three decimals leaves an exactly flat floor at the ground. On day 284 a 0.1 ft one put the
upwind boundary 20 km out (the 1e-12 inset over a 3 cm floor, nothing physical): no hand-over,
no transition, and flush microphones 22–28 dB low. The turning region is the heights within one
radian of vertical phase (k int sin dz, k at `Layers::resolution_frequency`, which a run sets to
its highest band) of the fastest air, and the limiting ray's travel across it is capped at
1.5 k h^2: a wave confined to h leaves it within about k h^2 ([THEORY §8.2](../THEORY.md)).

- **The 1.5 is set by what it must leave alone and what it may add.** Ray theory travels
  (4/3) k h^2 across a linear gradient's region (sqrt(2 R h) against k^2 h^3 = 9R/8), less on
  exact arcs, so a cap above that never binds there and `pass` comes back bit for bit. At
  k h^2 it binds on every linear profile (`the resolution limit leaves every resolved turning
  alone` fails). What a floor can still add is the cap's slack over (4/3) k h^2, and it grows as
  k^(-1/3). 2 k h^2 was the first choice; at 1.5 a calm floor moves the demo plates less
  ([FIDELITY](../FIDELITY.md#refraction-models-against-the-pe), "A calm floor"). Do not raise it
  for margin.
- **The highest band, not a typical one**: the geometry is one for all bands, and the largest k
  removes only what no band resolves. `ShadowTransition` moves each band's boundary toward the
  fit's as before. A consequence worth knowing: `--band_limits` moves a floored profile's
  boundary.
- **Do not replace it with merging thin layers.** Any thickness threshold leaves a cliff just
  above it, since a flat floor of any thickness diverges. The cap is continuous in the floor's
  thickness (0.1 ft: 186 -> 187 m; 1 m: 544 m at 10 kHz, 707 km uncapped), and the region, a
  set of heights, does not depend on how many levels describe the same air.
- **The shortcut is exact, not a heuristic**: a single fastest height whose radian is gathered
  within the segment on each side is a linear turning, and q's trapezoid floors the phase (q is
  concave on a layer). Only past that does it integrate (closed form, `vertical_phase`) and
  solve for the edges (safeguarded Newton). Cost: nothing measurable on `example_windy.csv`,
  about 3% on `example_inversion.csv`.
- **The region belongs to the air, not to the path** (`beyond`). Where the path ends before the
  radian, the walk carries on up through the profile and down to the ground, as far as the
  limiting ray could go (not into air faster than 1/pass), and only the travel inside the path
  is held to the cap. Stopping at the path's ends proves nothing: a source within a radian of
  the fastest height had a region as short as the path, and the cap bound on smooth gradients.
  Under a plain lapse a 1.5 m source over a 1.2 m microphone had its boundary brought in from
  789 m to 34 m at 10 kHz. The ground is a hard stop.
- **Not at the horizontal launch** (`launch_junction` in `families()`). Where c_rs is the
  source's own speed and the air above is faster, pass is where the downward family meets the
  families launched up to turn aloft, not an end of the rays. Capping it opened a gap between
  them inside the lit zone: receivers with no ray under any downward refraction. A flat stretch
  at the source's own height under such a profile is left to plain ray theory.
- With both of the last two right, a sweep of 728 log-law cases differs from plain ray theory
  by no more than 1e-7 dB (interior maxima under a light tailwind). The goldens did not move,
  but for five cells of `stratified_windy_point` by up to 4e-6 dB: near crosswind the lowest
  layer of `example_windy.csv` is within 1e-4 m/s of flat, and its rays ran 1–40 km before the
  cap.
- Point receivers above a floor were never affected (the dipping family turns at the floor's
  top, in its gradient). A profile whose lowest level is above the ground has a floor of that
  height, constant below it, and the cap follows a thick one for kilometers.
- `stratified_calm_floor_plate_is_floorless` holds a 0.1 ft calm floor
  (`atmospheres/example_windy_calm_floor.csv`) to 0.3 dB of the profile without it (39 dB
  unfixed). `stratified_calm_floor_plate_200hz_golden` is the only run at a top band other than
  10 kHz, so it is what fails if `process_trajectory` stops handing `LayeredSoundSpeed` its band
  (`a run's resolution frequency reaches every path's layers` covers `path_layers`).

## Interpolated states (`--ray_stride`)

`LayeredSolution` (stratified.hpp) is what a pair's layered solves give: the direct eigenray
after the boundary's removal under aloft, the hand-over, the ray and wave arrivals, and the
reflected ray's grazing and length for a point receiver's ground. `layered_direct_ray` and
`layered_finish_ground` take one (`memo`): what it holds stands in place of solving, and what it
lacks the solves fill in. `StridedLayeredRays` keeps an observer's anchors, one per thread;
`process_trajectory` hands it each observer (`begin_observer`) and asks it for every state's ray
and ground finish. What to keep:

- **Interpolate the solves, not the levels.** The sphere's direction, the ground, the absorption
  and the shadow transition are formed at each state from its own geometry and the interpolated
  ray. Interpolating levels would miss the dips between anchors that move the 10 dB window, as
  decimating the trajectory does ([performance.md](performance.md#the-sel-window-the-oversampling-and-the-windowed-peak-search)).
- **Anchors are fixed states** (every K-th and the last), each solved exactly as without the
  option, so a state's ray does not depend on which state the windowed search asked for first:
  windowed equals marched byte for byte. `stratified_windowed_is_exact` and its twins run at the
  default stride; `a strided ray does not depend on the order the states are asked for` asks for
  the states shuffled, some twice, after another observer's anchors, and finishes them in
  another order.
- **Tubes, lengths and times as fractions of the chord.** Taken as they are, 14% of states
  between anchors 8 apart were refused on the 2% tube limit by nothing but the tube's growth with
  range; as fractions of the straight chord the growth is exact and only the focusing is
  interpolated.
- **The hand-over is decided afresh** from the interpolated boundary and the state's own linear
  boundary; where it differs from the anchors' the state is solved. A state whose ray was solved
  is finished exactly too (`RefractedRay::interpolated`).
- **A ray from aloft, a duct's promotion, a caustic coming or going** change the arrivals' count,
  bounces or kinds, or move angles and tubes past the limits, so those states are solved.
- **Cost:** see [performance.md](performance.md#interpolated-rays). What is left at stride 4 is
  the rest of the evaluation (the shadow transition, the band levels, the database load and the
  histories), so larger strides gain little; in a duct with wave arrivals the arrivals' kinds
  change often enough that most states are solved anyway.
- **4 is the default**; `--ray_stride 1` solves every state and reproduces the results from
  before it to the bit (the `*_exact` goldens).
- Tests: `two states' layered solves interpolate only where they agree`
  (`tests/test_stratified.cpp`), `stratified_stride_{windy,duct}_is_close` (against every state
  solved, at 0.005 and 0.02 dB), `stratified_stride_is_not_exact`, the `*_exact_golden` pins,
  `stratified_windowed_is_exact` (now at the default stride), `ray_stride_needs_stratified`.

## Arrivals and ducts

Every arrival is summed in energy (`layered_arrivals`, `Arrival` in refraction.hpp). The primary
is the first family's root. The other families are not monotonic in general, so they are sampled
in launch angle (32 + 8 per hop the path can hold, at most 512) and every crossing of r is
solved. In a surface duct (air above the source faster than all of it below) the bounced rays
come too: X_m(p) = X_0(p) + m L(p), from the same samples. Each earlier bounce costs the
plane-wave |R|^2 (`add_bounce_loss`), and past every direct ray's reach the strongest is promoted
to primary. Focusing is capped per ray at `max_focusing_gain_db` (+10 dB); distinct arrivals add
in energy and can exceed it in total. All solved arrivals are retained (see
[HISTORY](../HISTORY.md#changes-that-affect-results) for the earlier 24-arrival limit).

- **`solve` refuses a bracket that does not land** (residual over 1e-6 r). A family whose end
  ray never arrives stands in as an infinity, so it brackets falsely. In a duct that reported
  the duct's top as a lit ray, 55 dB low.
- **`Jumps`: the range is discontinuous** wherever 1/p passes a local maximum of c that is the
  fastest air yet met above z_s, or below z_r. There the turning point leaps, and `solve` splits
  its bracket. Without the split, a sign change at a jump alone ran the full 200 iterations to
  be refused, under a convergence test (`hi - lo <= 1e-17 hi`) that a double cannot meet. That
  was a third of the solves in a duct and 90% of their cost; the test is that the estimate no
  longer falls inside the bracket.
- **The samples are shared across families** (`Row`: down, up, dip, ground pieces at one p).
  Per family they repeated the walks 3–5 times in a duct.
- **The sample count scales with the hops, and 32 is not generous.** A ray that bounces m times
  travels X_0(p) + m L(p), so the same interval holds m times the wiggles: aliasing in the hop
  count, not near-zero intervals (refining those helped 30%). The comparison of sample counts
  against a 2048-sample reference is in [FIDELITY](../FIDELITY.md#refraction-models-against-the-pe)
  ("Arrivals missed at caustics"); 32 + 8 per hop, capped at 512, is the one used. The grid is
  bunched toward the shallow end (t = t_a + (t_b - t_a) u^2), free, and worth most
  of the windy gain. The inversion's remainder is what is left at that cost: 16 per hop gets
  1.8% for another 20% of the arrival time, and a cosine bunching 2.6% at 8. The scaling costs
  4–8% of a stratified run, and moved the PE agreement on the inversion slightly for the better
  and on the windy profile not at all. `layered_arrivals`' `samples` argument exists for
  `the arrival search has converged in its sample count`. A base of 16 is far worse. Do not
  trade samples for speed.
- **The primary's bracket is split around it, not solved whole.** Whole, it returned the
  primary or one of the caustic pair beside it, at random.
- **Past the hand-over to the linear fit only bounced rays, and rays turned back from aloft,
  are searched** (`bounced_only`). The fit stands for every downward direct ray; rays from
  aloft exist there only under a `boundary_under_aloft`. Upwind under a lapse that skips the
  search entirely.
- **The peak search bounds the sum** by the loudest arrival's direction over the harmonic sum
  of their squared ranges (`bound_geometry` in niceops.cpp).
  `stratified_duct_*_windowed_is_exact` holds it.

## Wave arrivals: `--tunneling` and `--skip_fill`

Ray theory leaves two gaps that these opt-in Airy terms address (the equations are in
[THEORY §8.2](../THEORY.md); the comparison with the PE and with flight data is in
[FIDELITY](../FIDELITY.md#micro-shadows-and-skip-zones)).

- **Not fixed by default: an elevated duct whose returning rays turn just above the ground**
  (day 286, run 286288, 100 ft source toward 270). It is not warm ground air: the air is
  fastest at z = 0 because of the log-law headwind assumed below the 46 ft LIDAR gate. Rays
  turned back from the inversion near 330 ft turn up again just above the ground and reach
  neither a flush nor a 1.2 m microphone. The PE tunnels through (the measurements are in
  [FIDELITY](../FIDELITY.md#micro-shadows-and-skip-zones)), as a WKB factor
  exp(-2 int kappa dz) does for rays turning at 1–9 ft. The skip zone short of the duct's
  return is a second gap. Neither moves the SEL. It is listed under
  [Known limitations](../FIDELITY.md#known-limitations).
- **Three opt-in treatments were built and compared**: `--tunneling` and `--skip_fill`
  (`layered_wave_arrivals`, stratified.cpp), and a PE `--propagation_table` beside the ray
  model for the pairs it does not cover.
- **Wave arrivals carry their loss per band (`WaveCorrection`, `add_wave_loss`)** on a tube
  taken at `wave_reference_frequency` (10 kHz), so the loss is never negative and the peak
  search's bound on the tube still bounds them. They stay out of the boundary's removal and the
  promotion to primary in `layered_direct_ray`, which only ray arrivals take. A tunneled one
  takes no ground term, plate or not.
- **Near a tunneled family's caustic the coalescing pair is replaced by the caustic's
  band-averaged Airy field.** Left to the +10 dB focusing cap it put the return well over the PE
  ([FIDELITY](../FIDELITY.md#micro-shadows-and-skip-zones)).
- **The band's mean of Ai^2 is a table in the center's argument** (`BandAiryTable`,
  stratified.cpp), since both corrections' arguments go as omega^(2/3): its logarithm and slope
  every 1/128 from -60 to 64, cubic Hermite, zero past 72 (every sub-band then under the
  1e-300 floor in every band), direct elsewhere. Seven complex Airy evaluations per band per
  arrival were half of a duct run with both options (222 s against 93 s on an eighth of the
  demo grid); the table is within 1e-10 dB of them
  (`a wave arrival's loss is its band's seven Airy evaluations`).
- Cost on a day-286 run: `--tunneling` about +20%, `--skip_fill` 2–4x, a PE table about the
  ray model's own. In `example_inversion.csv`'s duct both together are about 8x the plain
  stratified run, nearly all of it the arrival search's layer walks.

## Resolved shadow (`--shadow_model resolved`)

An opt-in change to how `--refraction stratified` places and fills a ground-grazing shadow,
written for low sources under light wind (the R44 landings at Eglin), where the layered model's
shadow was 5–11 dB deeper than the PE's at 1000–1500 ft. The default, `legacy`, is the
transition above, bit for bit (`stratified_resolved_is_not_legacy` checks that the option does
something; the goldens and every other CTest check that the default did not move). The
equations are in [THEORY §8.3](../THEORY.md); the scores against the PE and flight data, and why
it stays opt-in, are in [FIDELITY](../FIDELITY.md#resolved-shadow-against-the-pe).

Three causes of the deep shadow were found, and the option answers the first and third:

- **A layer the wave does not resolve sets the boundary.** Ray theory fixes a ground-grazing
  boundary by c(z_s)/c(0) alone. On run 229355 the log law extrapolated below the LIDAR's 46 ft
  gate takes 0.25 m/s off in the lowest foot, which sets r_s at 446 ft where the fit's r_l is
  1306 ft and the PE's onset follows r_l. Below 4 kHz that layer gathers less than a radian of
  vertical phase.
- **The lit side's ground term collapses.** The equivalent straight geometry
  h' = L sin(grazing) - h_r goes to 0 approaching the boundary, so the ground model sits in its
  grazing null, 8–20 dB under the true field. A cramped band keeps that ray model up to the
  boundary.
- **The decay is tied to the boundary**, a = c0 root^2 / r_b^2. With the fit casting no shadow,
  the fallback multiple of r_s gave 0.07–0.11 /s, far too fast.

What to preserve:

- **Each band sees only the profile it resolves.** For band f it finds h_f where the vertical
  phase reaches `--shadow_phase` (1 rad, `default_resolved_phase`, fixed before any scoring; it
  is the radian `resolved_pass` uses), and below h_f replaces the profile by its tangent
  (`tangent_below`). r_s(f) comes from that profile, and r_l(f) is the path's fit corrected by
  the least-squares line through the tangent's change (`resolved_shadow_bands`). Both go through
  the same `band_boundary` law. On a linear profile the resolved profile is the profile itself,
  so nothing moves.
- **The decay keeps the boundary's curvature.** Taking it from the tangent's gradient over the
  creeping-wave layer was built and scored, and is worse on every upwind class of the synthetic
  set. It is kept only as `niceops_pe_compare --tangent_decay` (`set_tangent_decay`). The decay
  rate, s and the diffraction length stay c0 root^2 / r_b(f)^2, each band at its own boundary.
- **A cramped band takes the exact field** (`--cramped_ground series`, the default). Where
  r_b/F < 2 (sqrt(y_s) + sqrt(y_r) < 2 in Airy lengths), it takes the residue series at the
  source's own height, with the curvature whose boundary is the band's. It is tabulated per band
  over (xi, s, sqrt(y_s)), with sqrt(y_s) in 0–2 every 0.4 and xi every 0.2
  (`ShadowTransition::cramped_level`); the modes and H0(kappa r) are shared across heights
  (`upward_pressure_heights`). Near the source, where the series is not trusted (NaN), the band
  falls back to the ray model. `straight` (the ground model along the straight path through the
  lit side and the blend) and `ray` (as legacy) remain as `--cramped_ground` values. The choice
  among the three was made on the synthetic set alone, before the measured set was scored.
- **The per-band boundaries depend only on bearing, source height and receiver**, so a run
  tabulates them once per trajectory (`ResolvedShadowTable`: every `bearing_step` = 2 degrees,
  source-height nodes even in sqrt(z_s), interpolated in 1/r_b so a band with no boundary is
  0). Every band whose h_f lies in one layer shares one profile, so the build costs one boundary
  solve per layer touched, about 0.01 s on `example_windy.csv`. The cramped series table is
  built a (band, s) column at a time, when `cramped_level` first reads it (a `std::once_flag`
  per column): a column is one series call, about 14 ms. Built whole it was 8.5 s at `-j 1` for
  an AS350 approach at Amedee (28 bands, one plate; 58 s before the mode sums were trimmed,
  [performance.md](performance.md#the-residue-series-tables)), and that run never reads it:
  none of its bands is ever cramped. The cramped values do not feed the bound (V's largest
  power ratio does), which is what lets them wait. A pair pays one table lookup per band.
- **The bound stays exact**: the series and V are held to the transition's largest power ratio,
  and the straight ground term is within the ground's own bound.
  `stratified_resolved_{point,plate,straight_point,duct_point}_windowed_is_exact` check it, and
  the `[resolved][bounds]` Catch2 case sweeps it. `tests/test_resolved_shadow.cpp` pins a linear
  -0.02 /s lapse case (20 ft source, flush) against the PE.
- **One approximation**: the table indexes receivers by their height above the ground (the
  run's point heights and the plate), so scattered observers whose surveyed z is not 0 are
  looked up at that height.
- `--shadow_phase` and `--cramped_ground` without `--shadow_model resolved` are refused.

## Cost

- **Linear**: on the demo grid (B407, `example_windy.csv`, `--frame_bearing 210`), the acoustic
  loop goes from 1.6 to 2.3 s bare, from 2.1 to 3.6 s with `--ground_effect`, and from 5.4 to
  6.1 s with `--plate_table`.
- **Stratified, primary ray only**: when the layered solve went in, the demo grid's acoustic
  loop with a point receiver took 12.0 s (3.3 s `linear`, 2.1 s without refraction); a
  96-microphone scattered array such as `arrays/coyle_approach.csv` takes 35 ms.
- **The search for every arrival** adds about half again. A whole demo-grid run (B407,
  `--ground_effect`, `--frame_bearing 210`, 14 threads) goes from 23.6 to 37 s on
  `example_windy.csv` and from 57.5 to 86 s on `example_inversion.csv`. It was +107% and +76%
  before the jump splits, the shared samples and the hand-over skip.
