# Receivers, ground effect, terrain and weighting

The user-facing description is in the manual's [Ground effect](../MANUAL.md#ground-effect),
[Receivers](../MANUAL.md#receivers),
[Ground-plane microphones](../MANUAL.md#ground-plane-microphones) and
[Terrain](../MANUAL.md#terrain). The equations are in [THEORY.md](../THEORY.md) §5 and §6. The
tables' interpolation accuracy is in
[FIDELITY.md](../FIDELITY.md#ground-and-plate-tables), and the measured effect of ground and
terrain under [Ground effect](../FIDELITY.md#ground-effect) and
[Terrain at Amedee](../FIDELITY.md#terrain-at-amedee). How the ground path was profiled, and the
optimizations measured and rejected there, are in [performance.md](performance.md).

## Receivers

Every level is evaluated at the microphone, the observer's point raised by its height above the
ground (`process_trajectory`). Which receiver a microphone gets is a user rule, given in the
manual. What to keep:

- **One ground model and one `GroundTable` per distinct (receiver height, ground class)**, built
  once per trajectory, not per observer. An observer file's fifth column lets one run mix plates
  and poles, and its sixth puts each microphone over a `--ground_classes` class. Keys are kept in
  the order first met, so with one class the tables, logs and the refraction setup's
  `point_heights` are what they were with heights alone.
- **The peak search's ground bound is the larger of the plate's and the point receiver's**, which
  does not depend on height, so a mixed run needs no per-observer bound. The point receiver's is
  `point_ground_gain_bounds` over every class point receivers stand on, the elementwise maximum
  of each class's `ground_gain_bounds`.
- **`mixed_receivers_are_their_subsets`** checks a mixed run (`tests/data/mixed_receivers.csv`)
  against runs over each kind. It runs with `--exact_ground`, because a table's extent follows its
  own run's observers, so tables built for different observer sets differ at the last digits.
- **Ground classes** (`src/ground_classes.{hpp,cpp}`, read as in
  [inputs.md](inputs.md#the-ground-classes-file)). Each class is the run's `GroundModel` with
  the file's overrides, and each has its admittance (`class_impedance`, one `ground_impedance()`
  call per class per trajectory). The hot loop picks the observer's admittance once per observer
  (`impedance_here`) and passes it where the run's was passed, on flat ground and over terrain, so
  an evaluation costs nothing more. What to keep:
  - **One class is the run without classes, byte for byte**: `class_models` is
    `{options.ground}`, `class_impedance[0]` is the one call made without classes, and
    `ground_class_single_is_scalar_{point,plate,terrain}` (and the explicit-columns twin) compare
    the files. A class file's numbers in those tests are whole, so the reader's and CLI11's
    parsing cannot differ by an ULP.
  - **`ground_impedance_bands` is the first class's**, and the refraction models and a table's
    plate offset take it. That is right only because `--refraction` and `--propagation_table`
    are refused with classes (an assert in the refraction setup says so). Per-class shadow
    transitions and bounce losses belong with ground along the path, not under the microphone.
  - **Flush microphones are on one class**, refused otherwise, and the plate is bound over it
    (`ground_model` points at it). The check is `PlateReceiver::ground_mismatch`, so `bound()`'s
    message on a run without classes is unchanged.
  - **Rejected:** a `DirectionBounds` table per class (more cache traffic for a bound that moves
    little, since each class's bound is (1 + max|Q|)² and |Q| nears 1 at grazing for every
    surface; revisit only if the measured evaluated fraction moves: half concrete moved it from
    98.48% to 98.57%, [performance.md](performance.md#terrain-ground-classes-and-screening));
    scaling the bound per observer (impossible, since the ground gain is folded into the
    A-weighted direction table); substituting a lone class for `options.ground` (it would enable
    refraction over a ground the user did not name); ids as row indexes (they move when rows are
    reordered, and do not match land-cover codes).
  - **A later blend over surfaces** (Fresnel-zone weighting along each path) stays under the same
    bound only if it blends power ratios convexly; otherwise the bound must be re-derived.
  - Tests: `ground_classes_mixed_are_their_subsets` (`--exact_ground`, as
    `mixed_receivers_are_their_subsets`), `ground_classes_tabled_matches_exact`,
    `ground_classes{,_terrain}_impedance_is_the_observers`,
    `ground_classes{,_terrain}_windowed_is_exact`, and the refusals.
- The direct path ends at the microphone, not on the ground. Moving it there changed
  point-receiver results; see [HISTORY](../HISTORY.md#changes-that-affect-results).

## Urban receivers (`--urban_table`)

`src/urban.{hpp,cpp}` holds the facades, the first-order image test, the ideal canyon under
plane-wave incidence, and `UrbanTable`; `tools/urban_table.cpp` builds tables; `urban_level` and
`urban_bound_geometry` in `process_trajectory` are the run's side. The design and its staged plan
are the Claude Doc "NICEOPS: Urban Propagation Plan"; the literature behind it is in
`reports/Urban canyon propagation methods.md` and `reports/Urban propagation survey 2026-10-11.md`.
What to keep:

- **A table of 0 dB with no facades is the run without one, byte for byte.** The direct path goes
  through `spectral_level_along` with the table's loss as `extra_loss`, and subtracting 0.0 is
  exact; with no valid reflection `urban_level` returns that level as it came, without an
  energy round trip; the bound multiplies by `max(1, ...)`. `urban_table{,_plate}_open_is_open_ground`
  compare the files with `compare_files`.
- **First-order reflections are paths, not table entries.** Each is
  `straight_ray(source, mirrored microphone)`: the sphere is read along the line leaving the
  source, so it gets its own directivity. Putting them in the table, keyed on the arrival
  direction alone, roughly doubles the error at 30 to 120 m (FIDELITY, Urban receivers).
- **An image's ground model is evaluated exactly** (no `GroundTable`): a table's extent follows
  the run's observers, and a mirrored microphone lies outside it. The vertical mirror keeps
  heights, so the ground geometry of the unfolded path is exact.
- **An image's separation is mirrored back** before the plate reads its azimuth: the sound
  crosses the plate travelling the mirrored way. `urban_table_plate_golden` fails without it.
- **The bound** multiplies `path_gains` by the table's largest gain, and the direct path and the
  images share one bound geometry, built by the eigenray rule (largest numerator over the
  harmonic sum of squared ranges). The image term has no CTest that fails without it (see
  [testing.md](testing.md#mutation-guarded-tests)); keep it, since a strongly directional sphere
  can make a reflection louder than the direct path's bin bound.
- **The canyon table's specular part comes from four traces, not a lattice sum.** A plane wave
  reaches a point by at most one ray per cross-canyon direction and per way of arriving, from
  above or off the floor, so `canyon_ratio` traces four rays back from the microphone and its
  floor image; `specular_paths` (brute force, with the walls carried below the floor for the
  floor's rays) is the test reference only.
- **The floor is the table's, not the ground model's, inside a canyon.** The first version left
  the floor reflection to the ground model, as if it shared the direct ray's fate. The wave
  solver showed otherwise: from the floor image, 2 z_r lower, a ray clears the roofs after a
  different number of reflections, and at low elevations it is often the only one open. At 15
  and 30 degrees over a 20 m canyon the table read 20-25 dB of shielding where the wave solution
  has 2-4 dB. The table is therefore relative to open ground's two rays, 1 + rho_f.
- **A first-order reflection's floor partner is checked at run time**
  (`first_order_floor_partner`, against `urban_walls_below`): the reflection is evaluated with
  the ground model only when its partner is open too, and in free field otherwise, and the
  table holds the partner exactly when the run does not. Leaving the check out overcounts by
  3 dB in the directions where the two differ; both pinned goldens fail without it.
- **Receivers are matched by position** (1 mm), because the observer reader keeps no ids. The
  table records its frame bearing, and the run refuses another: azimuths are compass bearings
  and the geometry is local.
- **Stage-1 refusals**: `--terrain`, `--refraction`, `--propagation_table`, `--db_cutoff`.
  Lifting refraction means looking the table up along each eigenray's arrival direction;
  `--db_cutoff` needs `urban_bound_geometry` in `is_culled`.
- Tests: `[urban]` in `tests/test_urban.cpp` (two-family count and first-order images against
  `specular_paths`, the Keller-cone path difference, the table's lookup, bound and file), the
  hidden `[.urban_separability]` study, and the `urban_table_*` CTests (open-sky identity,
  windowed against marched for canyons, plates and the +25 dB bump table, the pinned goldens and
  the refusals).

## Weighting

`src/weighting.{hpp,cpp}` holds A-weighting and band-limited energy integration.

- **The A curve is panam's `dBAw`**, the IEC 61672 pole-zero form rather than a table, so it
  evaluates at whatever band centers a database carries. It is pinned in
  `tests/test_weighting.cpp`. In this form the curve is **not** identically 0 dB at 1 kHz but
  −1.49e-6 dB; a test asserting exact 0 there will fail.
- **`a_weighting_gains()` folds the weighting and the band limit into one per-band multiply**, with
  excluded bands exactly zero.
- **`weighted_level()` returns −∞ for silence, not 0 dB**, so an excluded spectrum cannot survive an
  energy sum as a quiet-but-real level.
- `c_weighting()` / `c_weighting_gains()` are the C curve for `--csel_file`, from the same
  constants; see [metrics.md](metrics.md#time-weighting-day-metrics-c-sel).

## Faddeeva function

`src/faddeeva.{hpp,cpp}` computes w(z) = exp(−z²) erfc(−iz). It exists solely for the ground
model's boundary-loss term, which panam writes as `exp(-w)(1 - erf(-i sqrt w))` using SciPy's
complex `erf`. With z = √w that *is* w(z), and evaluating it as written overflows.

- **Three methods, chosen by the argument** (`faddeeva()`):
  - |z| ≥ 6 with Im z ≥ 0, or x² − y² ≥ 36: a 12-term asymptotic series, 1e-11 relative.
  - |z| < 6 with Re z ≥ 0 and −0.25 ≤ Im z ≤ 6: a sixth-order Taylor series about the nearest node
    of a 0.05 grid (15,246 nodes, filled once at start-up from the Weideman series), 3.6e-12
    relative, losing digits further below the real axis.
  - Everything else: Weideman's 32-term rational approximation (`faddeeva_reference` is that series
    everywhere).

  The domains were sampled from real runs (480,000 arguments from point and plate receivers on the
  demo case) and each method measured against SciPy's `wofz` on them.
- **The Weideman series runs as four interleaved Horner chains over t⁴** rather than one over t.
  31 dependent complex multiply-adds is a latency chain, not a throughput problem, and it ran at
  about a quarter of the machine's floating-point rate. Splitting by the power of t modulo four
  gives chains a quarter as deep that issue together (the gain is in
  [performance.md](performance.md#the-spectral-and-ground-paths)). The coefficient table is left
  exactly as Weideman's construction produced it and indexed by
  stride, so it stays checkable against the source.
- **It is written in real arithmetic.** Every `std::complex` multiply carries a branch to recover
  NaN/infinity cases that cannot arise here. The two trailing divides are folded into one:
  `2p/d^2 + c/d == (2p + c d)/d^2`.
- **Accuracy.** Max relative error 3e-13 over the domain the ground model reached while the
  boundary loss was cut at |w| = 500 (Re z ∈ [0.09, 22.3], Im z ∈ [−1.2, 15.6]), against `wofz`;
  past it the tests hold it to `wofz` at 1e-11 out to |z| = 700.
- **The lower-half-plane reflection `w(z) = 2exp(-z²) - w(-z)` is safe** because a passive ground
  keeps z within 45 degrees of the real axis wherever it is below it, so exp(−z²) is at most 1 in
  size.

## Ground effect

`src/ground_effect.{hpp,cpp}` is the excess ground attenuation, ported from panam's `ega()` and
pinned against it to 1e-9 dB in `tests/test_ground_effect.cpp`. The impedance, turbulence and
roughness models follow Attenborough & Taherzadeh (2026), as panam's `ground_plane` module does,
and are pinned against `ground_plane.pole_level` there (1e-6 dB). A change to any of them is a
change to panam first.

- **It returns a gain, not a loss**: +6 dB over rigid ground with the receiver on the surface
  (and no turbulence).
- **`ground_impedance()` is hoisted out of the observer loop**, and the allocation-free overload is
  what the loop calls. Returning a vector there allocated 94 million times and cost 20% of the run.
- **The band loop in `excess_ground_attenuation` is real arithmetic**, with its band-invariant
  branches (boundary loss on/off, the sinc, turbulence) hoisted so the body stays straight-line.
  Bands outside the metric's limits are not computed at all: `first_band`/`last_band` bound a
  contiguous run, which is always what the limits select because band centers are sorted. The
  details, and what not to undo, are in [performance.md](performance.md).
- **The boundary-loss term has no cutoff.** It is evaluated by the Faddeeva function at any |w|,
  as panam's `ega()` does, and pinned to panam b758e3a (`ega` 1e-9 dB, `pole_level` and
  `board_level` 1e-6). The cut at |w| = 500 it replaced is in
  [HISTORY](../HISTORY.md#changes-that-affect-results). Two things to keep:
  - **The gain bound past the edge is analytic** (`reflection_coefficient_bounds`). Past
    |w| = 500, |Q| ≤ |R_p| + |1 − R_p| (1.003/(2|w|) + |P|), with P = 2i√π z exp(−z²) the surface
    wave, whose size peaks at |z|² = 1/(2 cos 2θ), θ the parameter's direction, which does not
    depend on range. For soils that peak is far inside the edge. For a nearly lossless stiff
    surface (a thin layer's band, θ a hair under −45°) the wave carries on: |Q| = 151 at
    |w| = 1250 at grazing. The sweep stops at the edge and the formula takes over; `a surface wave
    carried past |w| = 500 ...` holds it to a sweep out to κ = 1e10.
  - **A mass-like admittance (Im β > 0) keeps NICEOPS's branch, not panam's.** Its parameter is
    in the second quadrant. panam squares it and takes the principal root, the negative, whose w()
    grows as exp(|Re w|) and overflows by |w| ~ 1000. Such a surface has no surface wave, and the
    parameter itself gives none. A thin layer gives this in some bands.

### The ground table (`GroundTable`)

By default the point receiver interpolates a per-trajectory table instead of evaluating the model
at every observer and state; `--exact_ground` bypasses it.

- **Nodes in s = √(h_s/r) and ln r**, r = hypot(h_s, d): 641 angle nodes and 40 per decade by
  default. s resolves grazing incidence, where the model changes fastest. Values are exact at the
  nodes and interpolated bilinearly, as a power ratio, in between. Ranges outside [r_min, r_max]
  are clamped to the table's edge.
- **The interpolant is a convex combination of exact values**, so `ground_gain_bounds()` bounds it
  as it bounds the model, and the peak search and `--db_cutoff` stay exact.
- **Why it exists:** the model costs about 20 ns per band and the observer loop reaches it 94
  million times on the demo case, where it was 5.5 of 7.1 s. A lookup is four reads per band. The
  tables are built in well under 0.1 s.
- **The table holds one receiver height.** A path with its own receiver height (terrain, below)
  evaluates the model exactly.

## Plate receiver (`--plate_table`)

`src/plate_receiver.{hpp,cpp}` is a flush ground-plane microphone on a GR1425 plate, from panam's
axisymmetric boundary-element tables (`axisymmetric_bem.write_netcdf()`): p/p_free = P_d + Q P_r,
with P_d and P_r tabulated over band, sub-frequency, elevation and azimuth, and Q the ground's
spherical-wave coefficient at the geometry. The file comment in `plate_receiver.hpp` has the netCDF
layout.

- **The table records its ground, and `load()`/`bound()` refuse a mismatch** with the run's ground
  model and band centers. P_d and P_r depend on the ground.
- **Q excludes roughness and turbulence** (`spherical_reflection_coefficient`): the plate's image
  path is millimeters longer than the direct one, so neither decorrelates it.
- **P_d and P_r are interpolated bilinearly** (real and imaginary parts) in elevation and azimuth,
  held below the lowest and above the highest elevation, with azimuth periodic.
- **Three evaluations exist.** `level_exact()` takes Q at every sub-frequency, as panam's
  `board_level` does, and matches it to 1e-6 dB in `tests/test_plate_receiver.cpp`. `level()` and
  `power_ratio()`, which the run calls, take Q quadratic in log frequency through each band's
  edges and center (2n + 1 evaluations for n bands, not 5n), about 2.5x cheaper for the whole
  run. Its error against `level_exact` is in
  [FIDELITY](../FIDELITY.md#ground-and-plate-tables).
- **`tabulate()` tables those Q values** over the image path's grazing cosine (as √cos) and log
  length, 161 angles and 20 per decade, coarser than `GroundTable` because Q is smooth and the
  plate's interference is in P_d and P_r. `level_exact()` never uses it.
- **Orientation.** The tables put the microphone on +y of the plate. `plate_frame_azimuth` turns
  the propagation direction by the microphone's offset direction, and it is the one function both
  the straight path (`Database::spectral_level_along`) and the refracted path use.
  `--plate_mic_azimuth` defaults to `outboard` (`outboard_azimuth`: away from the track's nearest
  state, or to the track's right directly under it).
- **`PlateReceiver::gain_bounds`** bounds the plate for the peak search and `--db_cutoff`;
  `plate_windowed_is_exact` checks that skipped states change nothing.

The acoustic loop's cost with each receiver, tabled and exact, is in
[performance.md](performance.md#the-spectral-and-ground-paths).

## Terrain (`--terrain`)

`src/terrain.{hpp,cpp}`, and `terrain_ray` in `src/propagation.hpp`. A bare-earth grid (netCDF,
heights in meters; the validation harness's `terrain/dem_to_netcdf.py` makes one from USGS 3DEP).
Grid observers stand on it, and each pair's ground geometry is taken over the least-squares mean
ground line along its ground projection, observer to the source's footprint (CNOSSOS-EU's mean
plane, NORAH2 and ECAC Doc 32's practice over terrain). Spreading, absorption and the sphere lookup
keep the true 3-D path. The geometry is in [THEORY.md](../THEORY.md) §6.3.

### Things to preserve

- **Level terrain is the flat run, byte for byte** (`terrain_flat_*_is_flat`: point, plate and
  scattered). That rests on three spellings:
  - `Terrain::elevation` and `TerrainRaster::height` interpolate from a corner (`a + f (b - a)`),
    because the product form's four weights do not sum to exactly 1 and a constant grid came out
    1.9e-13 ft off;
  - `plane_geometry` divides by √(1 + b²) rather than multiplying by its reciprocal, so a level
    line divides by exactly 1;
  - the test grid is written as `-((-1250 * 0.3048) - (54.2651 * 0.3048))`, the height
    `--terrain_offset -1250` and `--ground_elevation 54.2651` take back to zero in the program's
    own arithmetic.
- **The ground's z is H + `--terrain_offset` − `--ground_elevation`**, the same
  `--ground_elevation` that is a WGS84 track's z = 0. NICEOPS converts no datums. Leaving out
  either is a plausible wrong answer, so the run checks what it can (the list is in the manual's
  [Terrain](../MANUAL.md#terrain)). The AC track's own lat/lon check uses `Trajectory::recorded`,
  kept as read, before any trim. Every refusal has a CTest on its message.
- **A point receiver's height above each line is its own** (`RayPath::receiver_height`, NaN for
  every other path), so `spectral_level_impl` evaluates the ground model exactly when it differs
  from the model's, *before* the `GroundTable` branch, which holds one height. Do not read the
  table at the pole's height instead: a height error of δ moves the path difference by 2δ sin
  (grazing), a radian at 1 kHz for 0.3 m at 30°. A plate lies along the line under its microphone
  and keeps its table (source height h_s − h_r). `a path's own receiver height is the ground
  model's, past its table` holds the branch.
- **Heights below the line are raised to zero**, CNOSSOS-EU's rule: a receiver in a dip, or a
  source behind rising ground, which is screened only with `--terrain_screening` (below). The
  run counts them and warns.
- **A negative distance along the line is taken by its size, not clamped to zero.** A nearly
  overhead source over ground falling toward it has its foot behind the microphone's. The image
  construction depends on d², so |d| keeps the model's direct path equal to the real range (zero
  shortened it), and `terrain_ray` turns `separation` round so a plate is crossed from the side the
  arrival actually comes from. `a source whose foot lies behind the microphone's keeps the real
  range` holds both.
- **`terrain_hills_point_golden` and `terrain_hills_plate_golden` pin the hills runs** at 1e-6
  (`golden/terrain_hills_*.csv`, from the C++ build). The windowed/marched pair compares one binary
  with itself, so a call-site mutation that handed the plate's microphone to `terrain_ray` as a
  point receiver moved cells by 0.01 dB and passed every other terrain test. Regenerate them only
  with a deliberate change, and say so.
- **The bounds need nothing new** without screening (with it, see
  [Screening](#screening---terrain_screening)). `ground_gain_bounds` and
  `PlateReceiver::gain_bounds` hold for every geometry and receiver height, and the bound test's geometry is the true 3-D path, as the
  level's spreading is. `terrain_hills_point_windowed_is_exact` and
  `terrain_hills_plate_windowed_is_exact` hold it. The tables' reach (`reach()` in
  `process_trajectory`) extends below z = 0 under terrain; node placement does not depend on the
  reach, so the flat case is unmoved.
- **Refused with `--refraction`.** The ray models reflect off and shadow onto z = 0; a tilted
  plane per path handed to them would be another model. That refusal is the only one; the one
  with `--propagation_table` had no CTest, and `terrain_table_refuses_refraction` holds that
  `--refraction` still wins when both are given.
- **A propagation table reads the same path.** `straight_level` finds each pair's `terrain_ray`
  path once, when a ground model or a table will read it, and hands it to `table_level` and then,
  for a pair the table misses, to the ground model. The plane counts as taken, and a raised
  height as raised, only for a pair one of them read: a pair the table misses with no ground
  model is free field. Level terrain stays byte-identical through the same three spellings
  (`terrain_flat_table_*_is_flat`).
  - **The bearing is the true separation's; the plate's azimuth is the turned one's.** Read at
    |d′| on a fixed bearing, an anisotropic table is continuous as d′ passes through zero; the
    turned separation's bearing jumps 180° there. The plate is still crossed from the side the
    arrival comes from. `a table read over terrain is continuous through d' = 0` holds both.
  - **A point receiver outside the table's rows falls back and is counted, never clamped**, for
    the reason in the bullet above: a height error moves the path difference. The run warns with
    the heights seen, and up front for a table with one receiver row.
  - **The bound needs nothing new.** The table term is its largest node (times the largest K
    node for plates), which no interpolation exceeds, whichever cell a pair reads, and the
    pairs that fall back are covered by the ground bounds above
    (`terrain_hills_table_*_windowed_is_exact`, the bump table among them).
- **Observers' ground is read from the grid itself in `run()`**, once, not from the raster, which
  is a resampling of it.

### Screening (`--terrain_screening`)

`src/screening.{hpp,cpp}` holds the CNOSSOS-EU kernel and the knife edge,
`GroundPlaneCache::may_block`, `blocked` and `fit_range` their terrain, and `screened_level` in
`process_trajectory` the CNOSSOS call (the knife edge's is below). The equations are
in [THEORY §6.4](../THEORY.md#64-terrain-screening). What to keep:

- **An unblocked pair takes the untouched call.** `straight_level` asks `blocked()` first and
  otherwise goes on exactly as without screening, so level ground, a uniform slope and the hills
  (where nothing blocks) are the unscreened runs byte for byte
  (`terrain_flat_*_screened_is_flat`, `terrain_tilted_*_is_unscreened`,
  `terrain_hills_unblocked_is_unscreened`, each with a log check that nothing was screened). A
  pair with no ground model and no table is never given a `terrain_ray` path: it stays
  `spectral_level_at_observer`'s.
- **The onset is strict.** A sample exactly on the line of sight does not block. CNOSSOS's own
  rule screens down to δ = −λ/20, but over level ground a sample near a 1.2 m microphone has
  δ ≈ h²/2s, about 0.14 m at 5 m, under λ/20 at 63 Hz, so it would screen level terrain. (The
  Directive also asks δ > λ/4 − δ′, δ′ by the ground images, which keeps the low bands out; the
  conformance case TC06, a path clearing a plateau's edge by 0.67 m, is diffracted at 500 Hz
  and 1 kHz in the TR and unscreened here, `ISO/TR 17534-4 TC06` in `tests/test_screening.cpp`.)
- **A minimum protrusion** (`--screening_min_height`, 0.5 m): a sample blocks only when it is
  that far above the line of sight. At Sweetwater the 1 m lidar grid's texture 3–5 m from a
  flush board, 0.03–0.14 m above a path at 1–2° elevation, screened it, and the screened level came
  out up to 1.1 dB above the unscreened one: with δ near zero the 4.77 dB floor and the two side
  ground terms beat the whole path's plate term. The threshold leaves the quick test (a necessary
  condition) as it was; it moved the ridge goldens by up to 0.9 dB (points) and 1.6 dB (plates)
  in about 90 of 441 cells, all at the shadow's edge, where it was regenerated.
- **One spelling of the slope** (`sight_slope`, (z − h_m)/s) for the horizon, the line of sight
  and the scan. The quick test reads each ray's horizon at the last sample strictly inside the
  path (ceil(D/Δ) − 1, never the footprint's), and with no turn reads ray a only, as `line_from`
  does.
- **The profile is one per pair**, blended in bearing with `line_from`'s own turn fraction, and
  the microphone is at its height above the raster's ground at its foot, the base of the running
  sums, so a flush microphone is never under its own profile. An observer file's z kept without
  `--observers_on_terrain` is where the microphone is evaluated, so the horizon is taken from it
  (`set_observer_at`), held no lower than the ground at the foot; taking the height above the
  raster instead screened such a microphone from where it is not (`terrain_ridge_screened_file_z`:
  3000 ft up, nothing is screened). A source under the profile at its
  footprint (states down to 10 m under the terrain are let through) is raised onto it.
- **No ground for the whole path on a screened pair** (`ground = nullptr`, `table = nullptr`):
  the ground enters only through the two side terms, so it counts once. A propagation table's
  level is not read either, since its ground is the whole path's; that is why the combination
  needs a ground model (refused otherwise) and the log counts the pairs that left the table.
- **F is clamped to 1** for the bound, and the run counts how often it binds. The bound's term is
  `screen_gain_bound(G_S, G_R)` = max(1, G_S) max(1, G_R)/3, with G_S over every class in use,
  a plate's included (the source side is always the point model), and G_R the point and plate
  bounds. Without a ground model it is 1/3 and changes nothing. The 3 is `cnossos_floor`,
  written once. The run logs how many bands the term raises above the ground term: the
  windowed-is-exact tests pass without it, so that line is what `terrain_ridge_screened_bound`
  checks.
- **Side lines come from the running sums** (`fit_range`, re-based to k0), so a screened pair
  fits nothing sample by sample; raw heights and the horizon cost about 16 bytes per sample per
  ray per thread, only with the option.
- **Rejected:** Deygout on the resampled profile (every sample of a rounded crest is a hull
  vertex, so a path just grazing it took about 13.6 dB, depending on `--terrain_step`); blending
  the two rays' losses in dB (a pair partly screened by a profile not its own, and a duplicate of
  the four-way ground selection on half-blocked pairs); a persistent per-ray hull (about 20 bytes
  per sample, and only useful with the dB blend); the sphere looked up toward the first edge
  (breaks the direction bounds' bins); a fade-in below δ = 0 (breaks the level-terrain identity).
- Tests: `terrain_ridge_screened_{point,plate,noground,classes,table}_windowed_is_exact` behind
  the ridge fixture (`integrated/terrain_ridge.nc`, 25 m high, 250 m to the track's side), each
  with a log check that pairs were screened; the point and plate runs pinned at 1e-6
  (`terrain_ridge_screened_{point,plate}_golden`); the refusals; and the Catch2 `[screening]` cases,
  including the bound sweep and the half-plane comparison.

#### The knife edge (`--terrain_screening kirchhoff`)

`knife_edge`, `knife_edge_loss` and the band-loss table in `src/screening.cpp`. In
`straight_level` a blocked pair computes its loss (`knife`) and goes on down the unscreened call,
which takes it off whichever level that reaches. What to keep:

- **The same blocked pairs as CNOSSOS.** It reuses `blocked()`, `--screening_min_height` and
  `rubber_band`, so where nothing is blocked the runs are the unscreened ones byte for byte
  (`terrain_flat_*_kirchhoff_is_flat`, `terrain_tilted_plate_kirchhoff_is_unscreened`,
  `terrain_hills_kirchhoff_is_unscreened`, each with a log check that nothing was screened).
- **The loss goes on every unscreened level**: the ground model over the mean plane
  (`spectral_level_along`'s extra loss), the table (`table_level`'s `added_loss`), and with
  neither the free field along the straight ray, which replaces `spectral_level_at_observer` for
  a blocked pair. Each call site has its golden; the table's and the free field's were added
  because dropping the loss there passed every windowed/marched pair.
- **Keep the whole path's ground, not CNOSSOS's two side terms.** Behind the Sweetwater hill
  `cnossos` came 29 to 48 dB below the unscreened level per blocked state from 500 Hz to 2 kHz
  (median over 11 runs), past its 25 dB cap on Δdif, where the measured shielding was 4–7 dB: the
  rest is its side terms in place of the path's ground.
- **No bound term.** |F(v)|² ≤ 1/4 for v ≥ 0 and falls monotonically, so the band loss is at
  least 10 lg 4 and a blocked pair is below its unscreened level, which the bounds already cover.
  `screen_gain_bound` is skipped, and behind the ridge the run evaluates the unscreened run's
  fraction (98.44%, against 99.26% for `cnossos`). A lit-side extension (v < 0, where |F|² reaches
  1.37) would need a term.
- **The edge** is the higher of the two horizon rays where they cross the aperture plane, which
  stands at the sample deepest above the line of sight (the first, on a tie). δ is summed as
  H²/(√(d² + H²) + d), which keeps its digits for H ≪ d. `the equivalent edge is the lowest point
  both ends see` checks it against a bisection on a peak, a hill and 200 random profiles.
- **The band loss is a table in ln x**, x = δ f_c/c, which is all it depends on: built on first
  use (a function-local static, so thread-safe), every 1/32 from 10⁻⁸ to 10⁶ (1033 nodes, about
  2 ms), with each node's exact slope for cubic Hermite. It is within 10⁻⁹ dB of the quadrature
  for x from 10⁻⁶ to 10³ and 6·10⁻⁸ dB at its top end; beyond its ends the quadrature is
  evaluated directly.
- **Rejected:** the eight-point quadrature per band at run time. Eight Fresnel integrals per band,
  each a series or a continued fraction, came to 50–64 µs per screened evaluation for δ of 0.1 to
  1 m, and the full ridge grid took 74 s against 8 s unscreened; the table takes 0.12 µs for 31
  bands.
- Tests: `terrain_ridge_kirchhoff_{point,plate,noground,table}_windowed_is_exact`, each with a log
  check that pairs were screened (the table's that they kept the table); all four runs pinned
  (`terrain_ridge_kirchhoff_*_golden`, the table's at 1e-5); `terrain_ridge_kirchhoff_is_not_cnossos`;
  and the Catch2 `[kirchhoff]` cases: the Fresnel integrals against SciPy, |F|² at grazing and in
  the shadow, the half-plane comparison, the band loss's floor and order, the table against the
  quadrature and the quadrature against Simpson, and the edge against bisection.
- **Cost** behind the ridge on the full demo grid (point receivers on variable-porosity ground of
  200 kPa·s/m², 8 threads, one run each on a shared M4 Pro at load 5 to 10, so read the ratios):
  59.1 s of CPU unscreened, 65.7 s `kirchhoff`, 77.5 s `cnossos`; peak RSS 239, 250 and 249 MiB.
  The knife edge's 6.6 s over 8.92 million screened evaluations is about 0.7 µs each, mostly the
  blocked test, the rubber band and the edge.

### How the plane is found (`GroundPlaneCache`)

A direct fit per pair is about 300 samples for a 3 km path, a microsecond, 94 million times on the
demo grid. Instead:

1. **Per trajectory**, the grid is resampled onto the local frame every `--terrain_step` (default:
   the grid's own spacing, at least 5 m) as a `TerrainRaster`: an inverse projection per node, not
   per sample. Its margin is two steps plus the diagonal times the bearing step, which is how far a
   bracketing ray reaches.
2. **Per observer**, profiles run along 1440 fixed bearings (0.25°, 13 m apart at 3 km), each
   sampled every step, kept as running sums of z and k·z. The least-squares line over any whole
   number of steps then costs a few flops (`line_from_sums`). A bearing is sampled lazily, only when
   a pair first needs it and only as far as that pair needs. Heights are summed relative to the
   observer's, or a 1000 m site cancels in `sum_kz - n/2 sum_z`.
3. **Per pair**, the line is interpolated linearly in bearing and in length between the four fits
   around it, which keeps it continuous. A path shorter than one step takes the first step's fit
   along its bearing, not a one-sample "fit" with zero slope, which was wrong for every nearly
   overhead state on a slope.

Along a bearing, at a whole number of steps, this is the direct least-squares fit
(`tests/test_terrain.cpp`, 1e-9 in slope). On a uniform plane it is exact up to the cosine's
curvature between bearings: a 0.05 slope errs by 1e-7, or 0.4 mm at 3 km
(2e-7 in slope at worst along a ray at a whole step).
On rolling synthetic terrain with ±15 m hills it stays within 0.15 m of the direct fit at both ends
of every path tried. The mutation that read one bearing twice fails two cases in
`tests/test_terrain.cpp`.

Cost on the full demo grid (S-76D track, B407 database, 14,400 observers × 6,551 states, 14
threads, `--humidity 20`), CPU time:

| | level ground (no `--terrain`) | level `--terrain` | hilly `--terrain` |
| --- | --- | --- | --- |
| `--ground_effect`, 1.2 m point receivers | 32 s | 35 s | 127 s |
| `--plate_table` | 71 s | 75 s | 75 s |

The planes cost about 3 s of CPU, roughly 35 ns a pair. Point receivers over real terrain cost 4x,
because each pair's receiver height differs from the one height the `GroundTable` was built for,
so the model is evaluated exactly. On level terrain the height is the pole's and the table is read.
Plates keep their table either way. For the 2017 arrays (about 40 microphones) none of this is
noticeable.

With propagation tables and screening, on the same grid at 8 threads on an M4 Pro (asserts on,
variable-porosity ground of 200 kPa·s/m², best of three; the machine was shared, load 10 to 31,
so read the ratios rather than the seconds). The table is the five-row ground-model test table
(`tables/ground_model_terrain_test.nc`: rows 0, 1.2 m, 8, 16 and 25 ft, source heights and ranges
every 20 ft):

| run | CPU | acoustic loop | evaluated | peak RSS |
| --- | ---: | ---: | ---: | ---: |
| level ground (no `--terrain`), 1.2 m point receivers | 23.2 s | 2.90 s | 98.48% | 231 MiB |
| hills, point receivers | 93.4 s | 11.8 s | 98.90% | 240 MiB |
| hills, point receivers, `--propagation_table` | 34.2 s | 4.32 s | 98.89% | 311 MiB |
| hills, `--plate_table` | 61.3 s | 7.72 s | 98.91% | 275 MiB |
| hills, `--plate_table --propagation_table` | 35.4 s | 4.46 s | 98.92% | 343 MiB |
| hills, point receivers, `--terrain_screening cnossos` | 98.4 s | 12.4 s | 99.63% | 249 MiB |
| ridge, point receivers | 61.8 s | 7.84 s | 98.44% | 240 MiB |
| ridge, point receivers, `--terrain_screening cnossos` | 80.0 s | 10.2 s | 99.26% | 250 MiB |

- **A table over terrain is cheaper than the ground model it replaces**: the point receivers'
  exact evaluation (4.0x the level run) becomes a lookup at their height above the plane (1.5x).
  The table covered 99.999% of the evaluations; the rest took the ground model. The extra 70 MiB
  is the table itself. Its accuracy at this row spacing is in
  [FIDELITY](../FIDELITY.md#propagation-tables).
- **Screening with nothing to screen** (the hills, byte-identical to the unscreened run) cost
  5% in CPU: the quick test per pair, and the bound's screened term, up to 4.1 dB above the
  ground term in 28 of 31 bands, which evaluated 0.7 points more of the pairs.
- **Behind the ridge** 9.27 million evaluations were screened (9.9% of those made, at 6057 of
  14,400 observers; 8.95 million at 6048 with the minimum protrusion added since), 3.51 million of them (38%) reached the 25 dB cap in some band, and F was
  held at 1 on neither side. The 18 s of extra CPU, the pairs the looser bound adds included, is
  about 2 µs per screened evaluation, as at one thread on the 441-observer grid (0.55 s for
  289,159).
- **The screening cache** holds about 16 bytes more per sample per ray. At these fixtures' 42.8 m
  step a path of 4.4 km is about 100 samples, about 2.3 MB per thread with every ray sampled:
  the peak footprint rose 3.4 MiB at one thread and 7.5 MiB at eight (441 observers), and 7.1
  MiB at eight on the full grid. At a 5 m step and 3 km, 600 samples, it is about 14 MB per
  thread.

### Terrain tests

- *Level terrain reproduces the run without it, byte for byte*: point receivers, plates, and
  scattered observers stood on it (`terrain_flat_*_is_flat`).
- *A uniformly tilted plane gives the closed-form h_s, h_r and d*, nearly overhead included, for
  point receivers and plates, and the rest of the path record stays the straight line's. The
  plane's own geometry is a rotation, hand-checked for rising and falling ground
  (`tests/test_terrain.cpp`).
- *The hills runs are pinned* (`terrain_hills_point_golden`, `terrain_hills_plate_golden`), and
  the windowed runs over hills equal the march (`terrain_hills_*_windowed_is_exact`).
- *The terrain reaches the ground model, not just the observers.* Microphones stood on the hills
  through an observer file, without `--terrain`, give different levels, and the same pair with no
  ground model agrees (`terrain_reaches_the_ground_model`,
  `terrain_without_ground_places_observers`). A level that reads the `GroundTable` at the wrong
  height is caught in `tests/test_spectral_propagation.cpp`.
- *A propagation table over terrain* is the flat read on level terrain, byte for byte
  (`terrain_flat_table_*_is_flat`), keeps the windowed search exact over the hills
  (`terrain_hills_table_*_windowed_is_exact`), counts its pairs toward the raised-height warning
  with no ground model (`terrain_hills_table_warns_unscreened`), warns of receiver rows
  (`terrain_table_warns_one_row`, `terrain_table_warns_receiver_rows`), and warns of a steep
  plane for a refracting table but not for a ground-model one
  (`terrain_steep_refracting_table_warns`, `terrain_steep_ground_model_table_is_quiet`, on the
  5.7° plane `integrated/terrain_steep.nc`). The geometry is checked in `tests/test_terrain.cpp`
  (`a table read ...`). A ground-model table over the hills is the terrain ground model to its
  interpolation (`terrain_hills_table_point_is_the_ground_model`, from 20 to 250 Hz over rows
  every half foot, `tables/ground_model_terrain_low_test.nc`; and
  `terrain_hills_table_plate_is_the_plate_model`); the measured error is in
  [FIDELITY](../FIDELITY.md#propagation-tables).
- *Every refusal has a CTest* that passes on its message.

What terrain changes at Amedee, and the slope sensitivity, are in
[FIDELITY](../FIDELITY.md#terrain-at-amedee).
