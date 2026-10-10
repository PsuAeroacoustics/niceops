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

- **One ground model and one `GroundTable` per distinct receiver height**, built once per
  trajectory, not per observer. An observer file's fifth column lets one run mix plates and poles.
- **The peak search's ground bound is the larger of the plate's and the point receiver's**, which
  does not depend on height, so a mixed run needs no per-observer bound.
- **`mixed_receivers_are_their_subsets`** checks a mixed run (`tests/data/mixed_receivers.csv`)
  against runs over each kind. It runs with `--exact_ground`, because a table's extent follows its
  own run's observers, so tables built for different observer sets differ at the last digits.
- The direct path ends at the microphone, not on the ground. Moving it there changed
  point-receiver results; see [HISTORY](../HISTORY.md#changes-that-affect-results).

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
  - **A mass-like admittance (Im β > 0) keeps NICE-OPS's branch, not panam's.** Its parameter is
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
  `--ground_elevation` that is a WGS84 track's z = 0. NICE-OPS converts no datums. Leaving out
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
  source behind rising ground, which is screening and not modeled. The log counts them.
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
- **The bounds need nothing new.** `ground_gain_bounds` and `PlateReceiver::gain_bounds` hold for
  every geometry and receiver height, and the bound test's geometry is the true 3-D path, as the
  level's spreading is. `terrain_hills_point_windowed_is_exact` and
  `terrain_hills_plate_windowed_is_exact` hold it. The tables' reach (`reach()` in
  `process_trajectory`) extends below z = 0 under terrain; node placement does not depend on the
  reach, so the flat case is unmoved.
- **Refused with `--refraction` and `--propagation_table`.** The ray models reflect off and shadow
  onto z = 0, and the table's heights are over flat ground. A tilted plane per path handed to them
  would be another model.
- **Observers' ground is read from the grid itself in `run()`**, once, not from the raster, which
  is a resampling of it.

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
- *Every refusal has a CTest* that passes on its message.

What terrain changes at Amedee, and the slope sensitivity, are in
[FIDELITY](../FIDELITY.md#terrain-at-amedee).
