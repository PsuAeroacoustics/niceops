# Receivers, ground effect, terrain and weighting

- **Receivers** (`process_trajectory`). Every level is evaluated at the microphone: the
  observer's point raised by its height above the ground. A flush microphone (height 0) with
  `--plate_table` is on a plate; any other, with `--ground_effect`, is a point receiver at its
  height, whose ground model and `GroundTable` are built once per distinct height. An observer
  file's fifth column (`id,x,y,z,height`) gives heights per observer, so one run mixes plates
  and poles, as test arrays do. Without it everyone takes the run's receiver. The peak
  search's ground bound is the larger of the plate's and the point receiver's, which does not
  depend on height. `mixed_receivers_are_their_subsets` checks a mixed run against runs over
  each kind (with `--exact_ground`, because a table's extent follows its run's observers).
  The point receiver's direct path ended on the ground until this went in, while its ground
  model put the microphone 1.2 m up. Point-receiver SEL moved by up to ±0.8 dB on the demo
  grid, almost all from the low states; plates and the goldens did not move.

- `src/weighting.{hpp,cpp}` — A-weighting (panam's `dBAw`, the IEC 61672 pole-zero form rather
  than a table, so it evaluates at whatever band centers a database carries) and band-limited
  energy integration. Note the curve is **not** identically 0 dB at 1 kHz in this form but
  -1.49e-6 dB; a test asserting exact 0 there will fail. `a_weighting_gains()` folds the
  weighting and the band limit into one per-band multiply, with excluded bands exactly zero.
  `weighted_level()` returns -infinity for silence, not 0 dB, so an excluded spectrum cannot
  survive an energy sum as a quiet-but-real level. `c_weighting()` / `c_weighting_gains()` are
  the C curve for `--csel_file`, from the same constants (see "Output metrics" in [metrics.md](metrics.md)).
- `src/faddeeva.{hpp,cpp}` — the Faddeeva function w(z) = exp(-z²)erfc(-iz), Weideman's 32-term
  rational approximation. The series is evaluated as **four interleaved Horner chains over
  t⁴** rather than one over t: 31 dependent complex multiply-adds is a latency chain, not a
  throughput problem, and it was running at about a quarter of the machine's floating-point
  rate. Splitting by the power of t modulo four gives chains a quarter as deep that issue
  together — 32% off `faddeeva`, 14% off the whole ground path. The coefficient table is
  left exactly as Weideman's construction produced it and indexed by stride, so it stays
  checkable against the source. It exists solely for the ground model's boundary-loss term, which panam
  writes as `exp(-w)(1 - erf(-i sqrt w))` using SciPy's complex `erf`; with z = sqrt(w) that *is*
  w(z), and evaluating it as written overflows. Max relative error 3e-13 over the domain the
  ground model reached while the boundary loss was cut at |w| = 500 (Re z ∈ [0.09, 22.3],
  Im z ∈ [-1.2, 15.6]), measured against SciPy's `wofz`; past it the tests hold it to `wofz` at
  1e-11 out to |z| = 700. The lower-half-plane reflection `w(z) = 2exp(-z²) - w(-z)` is safe
  because a passive ground keeps z within 45 degrees of the real axis wherever it is below it,
  so exp(-z²) is at most 1 in size.
- `src/ground_effect.{hpp,cpp}` — excess ground attenuation, ported from panam's `ega()` and pinned
  against it to 1e-9 dB. **It returns a gain, not a loss** — +6 dB over rigid ground with the
  receiver on the surface. `ground_impedance()` is hoisted out of the observer loop and the
  allocation-free overload is what the loop calls; returning a vector there allocated 94 million
  times and cost 20% of the run.
- **The spectral and ground paths were profiled, not guessed at.** Two rounds of intuition
  (cache layout, then `pow`) bought 14% between them; `sample` on a running process then showed
  the real costs and the same effort bought 2.5x. What it found, and what not to undo:
  - `std::sqrt` on a `std::complex` costs an **atan2**, because libc++ evaluates it in polar
    form. `ground_effect.cpp` called it three times per band; two were provably redundant
    (`sqrt(w)` *is* the numerical-distance parameter wherever its real part is positive, which is
    every passive ground but a mass-like one; see below) and the
    third is an algebraic root written out in the band loop. `std::abs`/`std::arg` on a complex cost a
    hypot and an atan2 and are gone the same way -- `|Q|cos(p + arg Q)` is
    `Re(Q)cos p - Im(Q)sin p`, and `|w| <= 500` is `norm(parameter) <= 500`.
  - `faddeeva` is written in real arithmetic. Every `std::complex` multiply carries a branch to
    recover NaN/infinity cases that cannot arise here, and the Horner loop does 31 of them.
    Its two trailing divides are folded into one: `2p/d^2 + c/d == (2p + c d)/d^2`.
  - `spectra_` is condition-innermost, matching `grid`. The netCDF file stores it the other way
    and the loader transposes once.
  - The whole band loop in `excess_ground_attenuation` is real arithmetic for the same reason,
    and its band-invariant branches (boundary loss on/off, the sinc, turbulence) are hoisted so
    the body stays straight-line.
  - Bands outside the metric's limits are not computed at all. `first_band`/`last_band` bound a
    contiguous run, which is always what the limits select because band centers are sorted.
  - The Weideman series behind `faddeeva` runs as four interleaved Horner chains rather than
    one; see `src/faddeeva.cpp`. It was 53% of the ground path on its own and is latency-
    bound, not throughput-bound.
  All of it is exact: the footprint moved by 2.8e-14 dB end to end, and the EGA tests still pin
  to panam at 1e-9 dB. Acoustic loop on the Be407 case, default band limits: broadband 0.34 s,
  `--spectral` 1.44 s, `--spectral --ground_effect` 10.9 s (was 41.4); the Horner split takes
  the last of those to about 9.4 s.

  The `*_fixed_load_factor.nod` databases carry spectra, so the path runs end to end on shipped
  data; it was also exercisable from a bare checkout when none did, and still is: `write_flat_spectral_database` in `tests/test_spectral_propagation.cpp`
  interpolates a sphere out of the shipped `S-76D_M3.nod`, hangs synthetic band amplitudes on
  it and writes it back with `write_sphere`, and the result is a one-group `.nod` that
  `--spectral` accepts. Its single condition makes the blend degenerate, so timings off it run
  optimistic, but the ground term -- which is most of the cost -- does not care. On the demo
  track that way: `--spectral` 1.36 s -> 1.16 s, `--spectral --ground_effect` 8.95 s -> 7.81 s.
  For the ground term alone, driving `excess_ground_attenuation` over a spread of source
  heights and ground distances with 26 synthetic third-octave bands reproduces the per-call
  cost the observer loop sees (~840 ns, was ~975), and `tests/test_ground_effect.cpp` -- pinned
  to panam at 1e-9 dB -- is what says it is still right.

  **The asymptotic boundary loss was tried and rejected.** For large `|w|` the term tends to
  `-(1/2w + 3/4w^2 + ...)`, which would skip `faddeeva` on ~20% of bands. It does not converge
  to the exact value -- the error plateaus at 1.6e-5 at `|w| >= 50` and 1.9e-6 at `|w| >= 100`
  no matter how many terms are taken -- so it cannot hold the 1e-9 dB pinning against panam at
  any cutoff in range. Roughly 9% of the ground path, for the loss of the property that proves
  the port. Not worth it; do not re-derive it.

  **There is no cutoff any more.** panam cut the boundary-loss term to zero past |w| = 500, and
  this followed, until panam e978703 evaluated it by `wofz` at any |w|; NICE-OPS now does the
  same, pinned to panam b758e3a (`ega` 1e-9 dB, `pole_level` and `board_level` 1e-6). Q moves by up
  to 2e-3 just past the old cut; the demo footprint with the plate receiver by at most 0.028 dB
  (0.08 dB with a 4 ft point receiver, where the two paths nearly cancel), means ~1e-4. Three
  things to keep:
  - **The gain bound past the edge is analytic** (`reflection_coefficient_bounds`). Past
    |w| = 500, |Q| <= |R_p| + |1 - R_p| (1.003/(2|w|) + |P|), P = 2i sqrt(pi) z exp(-z²) the
    surface wave, whose size peaks at |z|² = 1/(2 cos 2 theta), theta the parameter's direction,
    which does not depend on range. For soils that peak is far inside the edge. For a nearly
    lossless stiff surface (a thin layer's band, theta a hair under -45 degrees) the wave
    carries on: |Q| = 151 at |w| = 1250 at grazing, where the cut gave -1. The sweep stops at the
    edge and the formula takes over; `a surface wave carried past |w| = 500 ...` holds it to a
    sweep out to kappa = 1e10.
  - **A mass-like admittance (Im beta > 0) keeps NICE-OPS's branch, not panam's.** Its
    parameter is in the second quadrant; panam squares it and takes the principal root, the
    negative, whose w() grows as exp(|Re w|) and overflows by |w| ~ 1000. Such a surface has no
    surface wave, and the parameter itself gives none. A thin layer gives this in some bands.
  - The near-grazing three-term series that stood in past the cut (a refracted ray arriving
    level at a flush receiver, R_p -> -1) is gone: the exact term is what it approximated.

- `src/terrain.{hpp,cpp}`, `terrain_ray` in `src/propagation.hpp` -- `--terrain`: a bare-earth
  grid (the harness's `terrain/dem_to_netcdf.py` netCDF; USGS 3DEP), compiled but called by
  nothing until this. Grid observers stand on it, and each pair's ground geometry is
  taken over the least-squares mean ground line along its ground projection, observer to the
  source's footprint (CNOSSOS-EU's mean plane, NORAH2 / ECAC Doc 32's practice over terrain).
  Spreading, absorption and the sphere lookup keep the true 3-D path. The things to preserve:
  - **Level terrain is the flat run, byte for byte** (`terrain_flat_*_is_flat`, point, plate
    and scattered). That rests on three spellings: `Terrain::elevation` and
    `TerrainRaster::height` interpolate from a corner (`a + f (b - a)`), because the product
    form's four weights do not sum to exactly 1 and a constant grid came out 1.9e-13 ft off;
    `plane_geometry` divides by sqrt(1 + b^2) rather than multiplying by its reciprocal, so a
    level line divides by exactly 1; and the test grid is written as
    `-((-1250 * 0.3048) - (54.2651 * 0.3048))`, the height `--terrain_offset -1250` and
    `--ground_elevation 54.2651` take back to zero in the program's own arithmetic.
  - **The ground's z is H + `--terrain_offset` - `--ground_elevation`**, the same
    `--ground_elevation` that is a WGS84 track's z = 0. NICE-OPS converts no datums. At Amedee
    3DEP is NAVD88 and the 2017 tracks and microphones are ellipsoidal: offset -23.6 m, and
    `--ground_elevation` 3926.7 ft (the AC tracks' alt - z). Leaving out either is a plausible
    wrong answer, so the run checks what it can: an AC track's own lat/lon against the frame
    (refused past 30 m rms, warned past 3; `Trajectory::recorded`, kept as read, before any
    trim), its alt - z against `--ground_elevation` (1 m), any state more than 10 m under the
    terrain (refused, a datum error; less is warned, 3DEP's own error), and scattered observers'
    z against the terrain (warned past 1 m; `--observers_on_terrain` takes the terrain's). Every
    refusal has a CTest on its message.
  - **A point receiver's height above each line is its own** (`RayPath::receiver_height`,
    NaN for every other path), so `spectral_level_impl` evaluates the ground model exactly
    when it differs from the model's -- *before* the `GroundTable` branch, which holds one
    height. That costs 4x on the demo grid over real terrain (127 against 32 s CPU). Do not
    read the table at the pole's height instead: a height error of delta moves the path
    difference by 2 delta sin(grazing), a radian at 1 kHz for 0.3 m at 30 degrees. A plate lies
    along the line under its microphone and keeps its table (source height h_s - h_r).
    `a path's own receiver height is the ground model's, past its table` holds the branch.
  - **Heights below the line are raised to zero**, CNOSSOS-EU's rule: a receiver in a dip, or
    a source behind rising ground -- which is screening, and not modeled. The log counts them.
  - **A negative distance along the line is taken by its size, not clamped to zero.** A nearly
    overhead source over ground falling toward it has its foot behind the microphone's; the
    image construction depends on d^2, so |d| keeps the model's direct path equal to the real
    range (zero shortened it), and `terrain_ray` turns `separation` round so a plate is crossed
    from the side the arrival actually comes from. `a source whose foot lies behind the
    microphone's keeps the real range` holds both.
  - **`terrain_hills_point_golden` and `terrain_hills_plate_golden` pin the hills runs** at 1e-6
    (`golden/terrain_hills_*.csv`, made by the build that added them). The windowed/marched
    pair compares one binary with itself, so a call-site mutation that handed the plate's
    microphone to `terrain_ray` as a point receiver moved cells by 0.01 dB and passed every other
    terrain test. Regenerate them only with a deliberate change, and say so.
  - **`GroundPlaneCache`, not a fit per pair.** A direct fit is ~300 samples for a 3 km path,
    a microsecond, 94 million times. The cache keeps running sums of z and k z along 1440 fixed
    bearings per observer, sampled lazily and only as far as a pair needs, so any whole-step fit
    is a few flops (`line_from_sums`), and a pair takes the four fits around it interpolated in
    bearing and length (continuous; exact along a ray at a whole step, to 2e-7 in slope on a
    plane). A path shorter than one step takes the first step's fit along its bearing, not a
    one-sample "fit" with zero slope, which was wrong for every nearly overhead state on a
    slope. Heights are summed relative to the observer's, or a 1000 m site cancels in
    `sum_kz - n/2 sum_z`. The terrain is resampled once per trajectory onto the local frame
    (`TerrainRaster`: an inverse projection per node, not per sample), with a margin of two
    steps plus the diagonal times the bearing step, which is how far a bracketing ray reaches.
    About 35 ns a pair; the mutation that read one bearing twice fails two cases in
    `tests/test_terrain.cpp`.
  - **The bounds need nothing new**: `ground_gain_bounds` and `PlateReceiver::gain_bounds`
    hold for every geometry and receiver height, and the bound test's geometry is the true
    3-D path, as the level's spreading is. `terrain_hills_*_windowed_is_exact` holds it. The
    tables' reach (`reach()` in `process_trajectory`) extends below z = 0 under terrain; node
    placement does not depend on the reach, so the flat case is unmoved.
  - **Refused with `--refraction` and `--propagation_table`**: the ray models reflect off and
    shadow onto z = 0, and the table's heights are over flat ground. A tilted plane per path
    handed to them would be another model.
  - **Observers' ground is read from the grid itself in `run()`**, once, not from the raster,
    which is a resampling of it.
  - At Amedee, on the paper's plate configuration (two day-286 approaches, two day-284 steady
    runs), terrain moves SEL by at most 0.026 dB (rms 0.007); 1.2 m poles by up to 0.19 dB,
    because h_r above each plane carries the grid's decimeter-scale error at the microphone's
    foot. The evaluation's estimate for a 1-3 degree slope is ~1 dB for plates at 5-10 degrees
    of peak elevation (manual, "Terrain").
