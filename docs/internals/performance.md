# Performance: why the hot paths are shaped as they are

Measured results, including the optimizations that were rejected. Read this before you optimize anything, or "simplify" code that looks over-engineered.

## Where the database load time goes

A default run is now ~0.93 s, and ~0.61 s of that is `Database::Database`. Most of it is not
ours: `S-76D_M3.nod` has **one netCDF group per operating condition**, 1312 of them, and
netCDF-4 reads every group's metadata eagerly at `nc_open` -- it opens all ~8000 HDF5 datasets
in the file before returning. Measured split:

| | share of the load |
| --- | --- |
| `nc_open` (eager per-group metadata) | 41% |
| our variable reads (metadata lookup, not I/O) | 25% |
| `nc_close` (tearing that metadata down) | 12% |
| Delaunay build over the conditions | 14% |
| grid fill and the rest | 8% |

(The table is the load before the three changes below, when it was ~0.75 s; the first three rows
are unchanged and the last two shrank.) The first three -- 78% of that -- are netCDF's own cost
for this file layout and are not reachable through its API; the five tiny scalars per group cost more than the four big arrays, because
what they pay for is the lookup, not the bytes. What *was* ours, and is done: the Delaunay build runs on its own thread while `nc_close` tears
the file down (the build touches no netCDF state, so it is safe even though HDF5 here reports
thread-safety off; ~85 ms, `NcGuard::close`); each channel's grid cells are derived once and reused
while the angle labels are unchanged, which is every sphere but the legacy hover ones (the fill
and the spectral layout, 47 -> 10 ms); and the optional `ShaftPower` probe no longer logs a
line per group (1312 of a run's 1349 log lines). Output is byte-identical (golden compare).
The real fix is in the database builder
(panam): the same data in one group with `(condition, channel)` arrays would load in a
fraction of the time. Do not go looking for it in `sphere_database.cpp`.

## Why the hot path is shaped the way it is

A default run evaluates 14,400 observers × 6,551 states = 94.3M times, so anything left in `eval` costs ~100 M×. Do not "simplify" these back:

- **`Database::SourceState` is built once per trajectory state, in `main`, and reused across every observer.** It carries the source position, the `SphereFrame` (unit TPP normal, nose, starboard) and the `Blend` (four bracketing spheres + barycentric weights, or under `--blend nearest|aam` the selection `condition_selection.cpp` makes, padded to four slots). The barycentric `Blend` is a Delaunay `locate()` walk and was ~60% of the whole run when it sat inside `eval`. Anything that depends only on (μ, γ, CT) or on attitude belongs here, not in `eval`.
- **`lookup_angles` is split**: `sphere_frame(velocity, normal)` is per state, `lookup_angles(frame, direction)` is per observer. The four-argument overload still exists for the unit tests and composes the two.
- **`lookup_angles` takes three dot products and no projection.** `{nose, normal, starboard}` is orthonormal, so the propagation direction's components in the plane normal to the nose are just `normal.dot(d)` and `starboard.dot(d)`; forming the projected vector and normalizing it (a square root and three divides) scaled both by the same positive length, which the angle does not depend on. `phi` is then `atan(across / |up|)` below the plane, and `copysign(pi, across) - atan(across / up)` above it, which is `atan2(across, -up)` without the quadrant handling. Above the plane is read from the database's own upper half, when every condition group stores one (`upper_hemisphere_measured()`); a database without it folds the direction onto the lower half, as before. The 0/0 case at the poles is caught explicitly. The seam at 180 degrees is a real cell: `GridAxis::set_period(360)` makes the cell from the last phi knot to the first, one period on, and `locate()` interpolates across it.
- **The level grid is one flat `vector<Sample>`** (`{dBA, EAA}` pairs), indexed `(i_phi * n_theta + i_theta) * n_conditions + condition`, not `vector<vector<vector<double>>>`. Both metrics are read together, so they share a cache line. A condition-major layout was tried and is *worse* — it also breaks bit-exactness via different FMA contraction.
- **`GridAxis` replaces a binary search** with arithmetic when the axis is uniform, then corrects the estimate against the neighboring knots so the answer provably equals `upper_bound`'s. The searches alone were about a third of the loop. `tests/test_sphere_database.cpp` sweeps both paths against each other; keep that test if you touch it.
- **`DirectionBounds::bin` lives in the header.** The observer loop calls it ~94 million times from `niceops.cpp` while it is defined in `sphere_database.cpp`, so out of line it inlined only when `NICEOPS_LTO` was on and the toolchain supported it. It is fifteen instructions and should not depend on a build option; its tangent table is an `inline const` variable rather than a function-local static, which would test a guard on every call. The bin geometry — `bin`, the tangents, `narrowest_theta_bin_degrees`, `azimuth_bin_degrees` — sits on `DirectionBounds` together, because the table build has to pad by the same binning the hot path queries with, and three independent spellings of it is how that silently stops covering.
- **`GridAxis::locate` gives the cell and the interpolation weight together**, and is what `eval` and `eval_spectrum` call. Going through `bracket()` meant reading four knots back out of `values_`, four subtractions and a divide by the cell area, all to recover a number the uniform spacing already implies; `locate` returns the bracketing pair plus the fraction of the way across the cell, so the four corner weights are products of `(1-f)` and `f` and already sum to one. Worth 18% of the loop. It skips `bracket`'s correction loops, which is why it needs the tighter `exactly_uniform_` test in `reset()` rather than the tenth-of-a-step one: round-off can only put the estimate on the wrong side of a *knot*, where the two cells interpolate to the same value anyway. The wrap-around cell at the ends is not uniform and falls back to the knots. `tests/test_sphere_database.cpp` sweeps `locate` against `periodic_bilinear_interp` over uniform, irregular and wrapping cells — the two derivations have to agree, and that test is what says so.
- **Both `eval` and `eval_spectrum` build their stencil with `cell_stencil`.** It turns the two `GridAxis::Cell`s into the four corner indices and the four bilinear weights. Sharing `locate` alone got them the same *cell* but left each deriving its own weights, which is what keeps the claim below — that a broadband/spectral disagreement is the physics and never the interpolation — true rather than merely intended. `periodic_bilinear_interp` is no longer on either path; it is the reference `locate` is swept against in `tests/test_sphere_database.cpp`, and the cell/apply split it used to have (`BilinearCell`, `bilinear_apply`) existed only for the old `eval` and is gone.
- **The interpolation carries `{dBA, EAA}` as a two-lane vector** (`Pair` in `sphere_database.cpp`). Every step applies the same weight to both metrics, and a `Sample` is exactly a two-double vector in memory, so blending the four conditions at four corners and then the corners themselves is half as many instructions. Clang would not do it on its own from the scalar form or from a two-element array; the GCC/Clang `vector_size` spelling is explicit, with a plain-struct fallback that compiles to what the code said before.
- **Divides in the loop are multiplications by a stored reciprocal**: `inv_sphere_radius_`, `1e-3` for the EAA scaling, and one reciprocal for normalizing the propagation direction instead of three divides. A double divide is about four times a multiply. These move the result by an ULP, which is the only thing in this round that does.
- **The observer loop's chunks are dealt out on demand, not divided up front** (`src/util/parallel.hpp`). An equal split gives every worker the same slice, so the wall time is the slowest core's -- and the cores are not interchangeable. On an M4 Pro (10 performance + 4 efficiency) the old equal split ran no faster on all 14 threads than on 10, the efficiency cores' slices being the critical path; pulling chunks lets the fast cores take more of them and gets the 14-thread run back under the 10-thread one. Chunk size is `count / (8 * workers)`, which is fine-grained enough to balance and coarse enough that the per-chunk scratch allocation stays noise. The output still does not depend on the split -- every observer writes only its own slot -- so runs at any thread count remain byte-for-byte identical.
- **`trace()` is force-inlined and `Logger::write` is `noinline`** (`src/util/logger.hpp`). `eval` and `level_at_observer` make ~7 trace calls each; as ordinary calls that was a sixth of the loop spent entering functions that immediately returned. The level lives in `detail::min_level`, outside the `Logger` singleton, so testing it does not touch a function-local static's guard.
- **The database grid fill walks channels and places them**, rather than scanning all 1,406 channels for each of the 1,332 cells (2.5 billion `approx_equal` probes, most of the load time). Two semantics are load-bearing and preserved: a later channel still wins a shared cell, and `dba_max` covers matched channels only.

- **Link-time optimization is on** (`NICEOPS_LTO`, see CMakeLists.txt). The loop lives in `niceops.cpp` and `Database::level_at_observer` in `sphere_database.cpp`, so without it every one of ~100 million evaluations was an out-of-line call across a translation unit. Worth 18% of the loop on top of the source changes — and only 5% before them, because what it buys is hoisting and inlining that the divides and knot lookups were hiding.

Current cost on the demo case: database load ~0.75 s, acoustic loop 1.55 s of CPU (0.16 s wall on 14 threads). It was 2.2 s and 11.7 s originally, and 3.2 s of CPU before this work.

Where the 19.4 ns of an *evaluation* goes, from an ablation (`full`, `levels only`, `angles only`, `geometry only`, `eval at fixed angles`): 33% the two inverse trig calls, 35% `eval` itself (gather, blend, `log10`, arithmetic), 8% the vector geometry, 7% the SEL integral, and ~17% the difference between a real `eval` and one at synthetic angles -- which is *not* cache misses, see the gather item below.

Several further optimizations were measured and **rejected**; the numbers are here so they are not re-derived:

- **Hand-rolled `acos`, `atan` and `log10`.** They are a third of the loop, which makes them look like the obvious target, but this platform's libm does `acos` in 2.3 ns, `atan` in 2.0 and `log10` in 1.9 — around 8 cycles each, table-driven. A polynomial accurate enough to hold the golden (~1e-12 rad, so the level moves by less than 1e-10 dB) needs a square root and a degree-10 minimax and is *slower*; one cheap enough to win (degree 4, ~1e-5 rad, ~1e-3 dB on the footprint) saves only about 1.5 ns of the 3.55 an in-situ call costs. Bad trade in both directions.
- **`-fno-math-errno`.** The hope was that marking the libm calls pure would let the compiler overlap them. No measurable change: it will not software-pipeline across an opaque call, and the three calls in an evaluation are already independent of each other and so already overlap.
- **Storing the grid as `float`.** Would halve the 28 MB gather. It costs ~1e-5 dB, which is
  within engineering accuracy but ends the golden comparison as a bit-level regression test.
  Not taken without being asked -- and see the next item, which says it would not help anyway.
- **Anything aimed at the gather's locality.** The four blended conditions of a cell span 845 of
  the 1312 indices on average and touch 3.1 cache lines per corner, never 1, so an evaluation
  pulls ~800 bytes of lines to use 256 -- which looks like an obvious 3x to reclaim. It is not:
  the whole trajectory uses only **52 distinct blends**, so the four needed condition slices
  could be compacted into 4.4 MB of contiguous 64-byte rows (nothing like the rejected slab's
  +131 MB). Before building it, the ceiling was measured the cheap way -- run the same loop
  against a one-condition database, whose 21 KB grid makes the gather perfectly local -- and it
  came out at 1.82 s against the full database's 1.80 s. **The loop is compute-bound, not
  miss-bound**, which is what the slab note already said; this is an independent confirmation
  from the other direction. Measure that ceiling before touching the layout again.
- **Per-cell polynomial inverses for `acos`/`atan`.** The angles are never wanted for their own
  sake -- only a grid cell and a fraction across it -- so a polynomial per 5-degree cell mapping
  cos(theta) straight to the fraction would skip the inverse trig entirely. It cannot work:
  acos has a square-root singularity at the poles, so no polynomial in cos(theta) represents it.
  Degree 7 still leaves 0.31 degrees of error in the two polar cells, and 5e-5 degrees in the
  best interior ones; changing variables to sqrt(1 - a) to fix it puts a square root back.
- **Culling states the SEL window will discard, by bounding them one at a time.** 44% of
  (observer, state) pairs fall outside the window, but a bound that proves a *particular* pair
  is outside it has to be tighter than 10 dB, and `level <= sphere max - spreading` is 11.4 dB
  loose at the median -- it proves only 7.6%. That is the wrong question; asking it the other
  way round is what works, and is now the default path. See "The SEL window, the oversampling,
  and the windowed peak search" below.

## The SEL window, the oversampling, and the windowed peak search

**The window is the outermost 10 dB-down interval.** `ten_db_down_interval()` in
`src/util/sel.hpp` runs from the first to the last sample at or above peak - 10 dB, dips
below the threshold between them included. That is the duration of the aircraft noise metrics
(the EPNL duration of 14 CFR 36 A36.4.5, and what panam's `sound_exposure_level` uses, so a
prediction and a measurement are integrated alike). It replaced the D original's walk out from
the peak, which stopped at the *first* sample below the threshold on either side and so cut an
event off at its first dip: an approach heard from upstream lost its second rise, and on the
2017 approach data the two rules differ by up to 0.7 dB at a single microphone. The change was
made deliberately (commit `9d6e986`, which regenerated `golden/197421_LAE.csv`), so the old
first-crossing walk, `expand_window()` and the "left alone under the port-fidelity rule" note
are gone from the code; do not reintroduce them from an older copy of this file. The integral
is taken over *reception* time (`reception_factor`, 1 + dr/dt / c), with the step from the
trajectory itself (see the time-step note below).

**The SEL integral is oversampled, and the discrete window is why that cannot be used.** The
sound exposure integral is a Riemann sum over a 50 Hz trajectory, oversampled by roughly
thirty times. Decimating the demo case's 6551 states moves the footprint's **mean** by 0.008 dB
at stride 2 and 0.088 dB even at stride 64. What blocks a 2-8x cut in the whole acoustic loop
is the window: decimation steps over a transient dip, the window's end jumps, and the level
jumps with it. The error does not shrink with stride the way a quadrature error would --

| stride | max cell error, first-crossing (the retired rule) | outermost-crossing (current) |
| --- | --- | --- |
| 2 | 1.93 dB (104 cells > 0.1 dB) | 1.13 dB (86) |
| 8 | 2.34 dB (619) | 1.13 dB (271) |
| 64 | 2.48 dB (3150) | 1.54 dB (1346) |

-- which is the signature of a discrete window flip, not of the integral. So decimation is
**not** safe as a default under either definition: the mean is negligible but individual cells
move by 1-2 dB, and the cells are the deliverable. It stays available deliberately, through
`adapt()` and `-s` (both currently commented out in `main`).

**Skipping states the window cannot use: the windowed peak search.**
`windowed_sound_exposure_level` / `windowed_exposure` in `src/niceops.cpp` do not evaluate every
emission time. A coarse sweep (stride nt/128) seeds a best level; one pass then evaluates every
state `DirectionBounds` cannot prove is below that best level, handing the same geometry on to
the evaluation so a state that survives the test pays nothing for it; and
`ten_db_down_interval` finds the window's ends by scanning in from both sides, skipping a state
the bound puts below the threshold (`certainly_below`) and evaluating everything between the
ends. It is **exact** -- the output is byte-identical to evaluating every state, on every track,
grid, observer set, receiver and thread count tried -- because the integral and the window are
the shared `ten_db_down_interval()` and `sel_over_window()`, templated on a level accessor (a
filled array for `peak_range`, a memoizing lambda here), not a second copy.

It runs on the **spectral** path only: `--spectral` (the default whenever the database carries
spectra), `--ground_effect`, plates, the propagation table (whose largest value per band enters
the bound), and `--refraction` (the
bound is taken along each ray's launch direction at its ray tube's range; for a sum of
stratified arrivals, the loudest arrival's direction over the harmonic sum of their squared
ranges, `bound_geometry`). It is skipped, and the loop marches every state, for `-a` (the level
is transformed after eval), `-i` and `--band_history` (the whole history is written),
`--epnl_file` and `--band_frames` (every state's bands go into the frames, and the EPNL's window
is found on PNLT, which the A-weighted bound says nothing about),
`-D`, `--lden_file` and `--cnel_file` (the window is found on a penalized history, which a
window found on the raw one need not cover), `--lamax_file`, `--nax_file` and `--ta_file` (the
time weighting runs over the whole history), `--csel_file` (the C-weighted window is not the
A-weighted one, and every state's bands are needed to find it), `--full_history`, `--db_cutoff`
(which culls on its own terms), and any database with a negative excess-attenuation factor,
where dropping that term is not a loss.

What it saves depends on how dear an evaluation is. Under the outermost-crossing rule it
evaluates most of what it is shown: on the broadband demo case it evaluated 93.9% of the pairs
(the window itself is 59%, and the bound proves only about 6% of the pairs outside it) and ran
about 4.6% *slower* than marching every state (measured twice, output identical), so **broadband
runs no longer build the table or use the search**; `Database::direction_bounds`, the broadband
builder, remains, tested, for a broadband search that pays. A spectral evaluation costs 150 ns
and a `--ground_effect` one a microsecond, and there the search pays. The main pass keeps each
skipped state's numerator and squared range (`skipped_numerator`, `skipped_squared_range`) so
that `certainly_below` does not recompute the geometry -- a full ray solve under refraction --
for the 30 to 40% of states it is asked about; the answers are the same two numbers compared the
same way, so the output is unchanged.

Three things had to hold, and each was measured before any of it was written:

- **The bound needs less than 10 dB of slack**, or the peak search reaches outside the window
  and there is nothing left to save. The whole-sphere maximum is 11.4 dB, just over. Binning it
  by direction -- `DirectionBounds`, bins uniform in cos^2(theta) crossed with azimuth sectors
  -- brings the median to 8.6 dB.
- **The test must cost nothing on the states it does not reject.** It forms the squared range
  and the three dot products with the sphere frame, which is exactly what an evaluation needs,
  so `Database::level_from_geometry` takes those and a surviving state pays nothing for having
  been tested. Letting the test compute its own geometry and then recomputing it inside
  `level_at_observer` was measured: it turns a 1.48x into a 1.08x.
- **The table has to stay in cache.** One table per *distinct blend*, not per state -- a shipped
  track resolves thousands of states onto a few dozen tetrahedra -- so it is tens of kilobytes.
  Built per state it was 3.3 MB, re-streamed for every observer, and that alone cost more than
  the scheme saved once the loop was threaded.

The bound maxes over the four conditions of a tetrahedron rather than over the particular blend
of them, so it is looser than a per-state table would be; buying that back means a table per
state, which is the cache problem above. Bin counts were swept under the earlier window rule
(8x8 to 32x16) and the choice hardly mattered; 16x8 is kept.

The whole thing rests on one property, that the bound is never exceeded, and it is tested
directly (`tests/test_sphere_database.cpp`, `tests/test_spectral_propagation.cpp`, tag
`[bounds]`) over sweeps of source states and observers, **including sweeps aimed within five
degrees of the nose and tail axes**. A cell is filed under the bin of its center direction, and
at theta = 0 or 180 that direction has no azimuth, so a pole cell is filed under every azimuth
bin of its theta bin (`Database::home_bins`); before that, a query a degree off the axis was
bounded up to 0.8 dB too low and the ordinary sweep, which samples hardly any observers there,
never saw it.

Two suggested optimizations were deliberately **not** taken:

- Storing reciprocals of the (max − min) normalizing ranges: a reciprocal multiply is not the same as a divide, so it would move the goldens, and the normalization is per state now anyway.
- Hoisting the output-directory check out of `write_time_history_file`: one `exists()` next to creating and writing a ~130 KB file is noise, and it would mean bypassing `open_output`'s guarantee.

**The per-state slab was built and measured; do not rebuild it.** Collapsing the database to one blended 36×37 slab per trajectory state (with the observer loop blocked so `state` can be the outer loop) works and is bit-exact, but it buys only 6.8% of the acoustic loop for +131 MB of resident memory, and the condition-major grid layout it requires makes the ordinary blend path — `self_test`, and any observer set below the break-even block size — about 6% *slower*. Block size barely moves it, and a single 768 MB block is the worst configuration of all, because the level buffer's own bandwidth then dominates. The work is preserved on the `slab-experiment` branch with the full numbers.

The reason the intuition fails is worth keeping: observers sweep a 120×120 ground grid, so consecutive evaluations land on neighboring (phi, θ) cells and the gather is already cache-resident. It was never miss-bound. Measuring the ceiling first — point `eval` at a fixed L1-resident slab and time it — costs ten minutes and would have shown this before any restructuring.

**The SEL integral takes its time step from the trajectory, not from `-s`.** `sound_exposure_level` weights every sample by a time step, and `-s` is only the resampling step for `adapt()` — the two coincide *only* when `adapt()` actually resamples onto that grid, and it is commented out in `main`. Every shipped track is sampled at 0.02 s, so using the 0.5 s default as the weight overstated every level by exactly `10*log10(0.5/0.02) = 13.98 dB`. `main` now derives the step via `uniform_sample_interval()` and logs it, warning if the sampling is not uniform. If you re-enable `adapt()`, it resamples to `-s` and the derived step follows automatically.

## Absorption setup and multi-band table lookup (2026-10-07)

Two further optimizations retain the formulas and per-band summation order:

- `AbsorptionTerms` computes the classical loss, two exponential numerators, and
  temperature power once per spectrum, alongside the existing relaxation frequencies.
  This reduces repeated transcendental work when building layered absorption profiles.
  Its evaluation method stays out of line so GCC uses the same compiled expression for
  scalar and batch calls; otherwise inlining contracts them differently and breaks their
  exact-equality test.
- `PropagationTable::lookup` accumulates four bands together, sharing each corner's
  address and weight. The scalar tail handles remaining bands. Each sum still starts
  at zero and visits corners in the original order. No approximation, reduced precision,
  changed table resolution, or change to which states are evaluated is involved.

A local Apple Clang `-O3 -std=c++17` microbenchmark, without LTO, measured the following
medians over seven alternating before/after runs:

| Work | Before | After | Time reduction |
| --- | ---: | ---: | ---: |
| 2 million 24-band table lookups | 0.454 s | 0.216 s | 52% |
| 200,000 24-band absorption spectra | 0.0881 s | 0.0225 s | 74% |

The synthetic table had four bearings, three source heights, three receiver heights,
four ranges, and sinusoidal cell values; queries varied bearing and range through the
interior. Absorption used temperatures from 280 to 299.9 K and bands from 50 Hz in
third-octave steps. These are routine timings, not whole-run speedups; database I/O and
other propagation work remain unchanged.

A separate before/after comparison covered 10,000 table queries (240,000 band values)
and 64 atmospheres (1,536 coefficients). Clang matched bit for bit throughout. GCC 16
matched table values bit for bit; absorption differed by at most 5.6e-17 dB/m due to
compiler rounding. Tests retain exact scalar/batch agreement and check the original
absorption formula across temperature, pressure, and humidity ranges. Table tests sweep
even/uneven grids, reordered columns, partial band ranges, and NaN gaps.

## Atmospheric profile setup and queries (2026-10-07, second pass)

`AtmosphereProfile::with_hydrostatic_pressure` previously integrated from the bottom
again for every output level, making N levels cost O(N²) segment evaluations. A running
integral now makes it O(N), with the same addition order and the same ground reference.
The original endpoint interpolation is preserved separately from the complete-segment
prefix so that even its rounding stays the same.

`AtmosphereProfile::at` now locates the common height cell once for temperature, humidity,
and pressure instead of three times. `wind_at` likewise shares one search between its
two components. Neither interpolation nor endpoint holding changes.

Seven alternating Apple Clang `-O3 -std=c++17` microbenchmarks on a synthetic 2,000-level
profile (uneven temperatures, heights 2.5 m apart) gave these medians:

| Work | Before | After |
| --- | ---: | ---: |
| 20 complete hydrostatic pressure reconstructions | 95.69 ms | 0.390 ms |
| 2 million paired atmosphere/wind queries | 142.34 ms | 78.32 ms |

The approximately 245x pressure speedup applies to a large profile; small profiles have
much less work to save. Field queries took 45% less time. These are setup/query timings,
not end-to-end acoustic speedups. Exact-equality tests against the original pressure
walk and interpolation formulas pass under both GCC 16 and Apple Clang, including
singleton profiles, uneven layers, nearly isothermal layers, and surface-reference offsets.

A source/receiver-height cache for propagation tables was rejected in this pass. Five
alternating GCC runs on B407/S-76D trajectory 197421AC with `tables/ground_model_test.nc`,
`--gridstep 328 -j 1 --humidity 20`, gave 2.953 versus 2.952 s median acoustic time
(6,208,943 evaluations); output was byte-identical, but the additional cache and public
API did not earn their complexity.
