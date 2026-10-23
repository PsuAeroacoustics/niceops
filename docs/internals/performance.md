# Performance: why the hot paths are shaped as they are

Measured results, including the optimizations that were rejected. Read this before you optimize
anything, or "simplify" code that looks over-engineered. Costs of the refraction models are in
[refraction.md](refraction.md#cost), and of the PE in [pe.md](pe.md).

## Where the database load time goes

A default run is about 0.93 s, and about 0.61 s of that is `Database::Database`. Most of it is
not ours: `S-76D_M3.nod` has **one netCDF group per operating condition**, 1312 of them, and
netCDF-4 reads every group's metadata eagerly at `nc_open`: it opens all of the file's ~8000
HDF5 datasets before returning. The measured split, taken when the load was about 0.75 s:

| | share of the load |
| --- | --- |
| `nc_open` (eager per-group metadata) | 41% |
| our variable reads (metadata lookup, not I/O) | 25% |
| `nc_close` (tearing that metadata down) | 12% |
| Delaunay build over the conditions | 14% |
| grid fill and the rest | 8% |

The first three rows, 78% of that, are netCDF's own cost for this file layout and are not
reachable through its API. The five tiny scalars per group cost more than the four big arrays,
because what they pay for is the lookup, not the bytes. Three changes took the last two rows
down to today's figure, with byte-identical output (golden compare):

- The Delaunay build runs on its own thread while `nc_close` tears the file down (`NcGuard::close`;
  about 85 ms). The build touches no netCDF state, so it is safe even though HDF5 here reports
  thread-safety off.
- Each channel's grid cells are derived once and reused while the angle labels are unchanged,
  which is every sphere but the legacy hover ones (the fill and the spectral layout,
  47 -> 10 ms).
- The optional `ShaftPower` probe logs once, not once per group (that was 1312 of a run's 1349
  log lines).

The real fix is in the database builder (panam): the same data in one group with
`(condition, channel)` arrays would load in a fraction of the time. Do not go looking for it in
`sphere_database.cpp`.

## Several databases in one run

`run()` holds its databases in a `std::map<std::string, Database>` keyed by `database_key`, the
weakly canonical path, so two spellings of one file are one database (and one database to the
`--provenance_file` and `--screen_database` refusals). The policy:

- **Each distinct database is loaded once, on its first trajectory, and erased after its last.**
  `last_use` is computed from the expanded `-d` list before the loop. The erase sits outside the
  trajectory's `try`, so a trajectory that fails still releases its database.
- **The front database is loaded before the loop and pinned.** `--export_sphere`,
  `--export_aam` and `--self_test` read it and stop. After the loop `--provenance_file` and
  `--provenance_spheres` read it, which is sound only because the provenance refuses several
  databases (asserted there).
- **An entry is made only once its load succeeds** (`try_emplace`, which inserts nothing when the
  constructor throws). An entry for a database that failed to open would let a later trajectory
  naming the same file skip the load. The manifest's attributes for a database are taken at
  that moment, while it is open; the front one's are read when the manifest is written.
- **Trajectories are never reordered by database.** The run totals are summed in `-t` order, and
  reordering would change them in the last bit.

What it costs, measured with `/usr/bin/time -l` on a coarse grid (`--gridstep 2000`): a
database held after its load is tens of megabytes (B407_fixed_load_factor about 33 MB), while
the load itself peaks far higher and returns most of it: about 540 MB for `S-76D_M3.nod` and
180 MB for `B407_fixed_load_factor.nod`. So holding a database until its last use is cheap
against loading it twice, and the peak is set by the largest single load, not by the number
held. Holding one database and reloading on every change would load A twice in
`-d A -d B -d A` (0.6 s each for the S-76D). There is no least-recently-used cap: with the
retained sizes above, a run would need hundreds of distinct databases before the held ones
rivaled one load's peak.

## The hot path

A default run evaluates 14,400 observers × 6,551 states = 94.3 M times, so anything left in
`eval` costs about 100 M times over. Do not "simplify" these back:

- **`Database::SourceState` is built once per trajectory state, in `process_trajectory`, and
  reused across every observer.** It carries the source position, the `SphereFrame` (unit TPP
  normal, nose, starboard) and the `Blend` (four bracketing spheres and barycentric weights, or
  under `--blend nearest|aam` the selection `condition_selection.cpp` makes, padded to four
  slots). The barycentric `Blend` is a Delaunay `locate()` walk and was about 60% of the whole
  run when it sat inside `eval`. Anything that depends only on (μ, γ, C_T) or on attitude
  belongs here, not in `eval`.
- **`lookup_angles` is split**: `sphere_frame(velocity, normal)` is per state,
  `lookup_angles(frame, direction)` is per observer. The four-argument overload still exists for
  the unit tests and composes the two.
- **`lookup_angles` takes three dot products and no projection.** `{nose, normal, starboard}`
  is orthonormal, so the propagation direction's components in the plane normal to the nose are
  just `normal.dot(d)` and `starboard.dot(d)`. Forming the projected vector and normalizing it
  (a square root and three divides) scaled both by the same positive length, which the angle
  does not depend on. `phi` is then `atan(across / |up|)` below the plane, and
  `copysign(pi, across) - atan(across / up)` above it, which is `atan2(across, -up)` without the
  quadrant handling. Above the plane is read from the database's own upper half when every
  condition group stores one (`upper_hemisphere_measured()`); a database without it folds the
  direction onto the lower half. The 0/0 case at the poles is caught explicitly. The seam at
  180 degrees is a real cell: `GridAxis::set_period(360)` makes the cell from the last phi knot
  to the first, one period on, and `locate()` interpolates across it.
- **The level grid is one flat `vector<Sample>`** (`{dBA, EAA}` pairs), indexed
  `(i_phi * n_theta + i_theta) * n_conditions + condition`, not
  `vector<vector<vector<double>>>`. Both metrics are read together, so they share a cache line.
  A condition-major layout was tried and is *worse*; it also breaks bit-exactness through
  different FMA contraction.
- **`GridAxis` replaces a binary search** with arithmetic when the axis is uniform, then
  corrects the estimate against the neighboring knots so the answer provably equals
  `upper_bound`'s. The searches alone were about a third of the loop.
  `tests/test_sphere_database.cpp` sweeps both paths against each other; keep that test if you
  touch it.
- **`GridAxis::locate` gives the cell and the interpolation weight together**, and is what
  `eval` and `eval_spectrum` call. Going through `bracket()` meant reading four knots back out
  of `values_`, four subtractions and a divide by the cell area, all to recover a number the
  uniform spacing already implies. `locate` returns the bracketing pair plus the fraction of
  the way across the cell, so the four corner weights are products of `(1-f)` and `f` and
  already sum to one. Worth 18% of the loop. It skips `bracket`'s correction loops, which is why
  it needs the tighter `exactly_uniform_` test in `reset()` rather than the tenth-of-a-step one:
  round-off can only put the estimate on the wrong side of a *knot*, where the two cells
  interpolate to the same value anyway. The wrap-around cell at the ends is not uniform and
  falls back to the knots. `tests/test_sphere_database.cpp` sweeps `locate` against
  `periodic_bilinear_interp` over uniform, irregular and wrapping cells; the two derivations
  have to agree, and that test is what says so.
- **Both `eval` and `eval_spectrum` build their stencil with `cell_stencil`.** It turns the two
  `GridAxis::Cell`s into the four corner indices and the four bilinear weights. Sharing
  `locate` alone got them the same *cell* but left each deriving its own weights. One stencil is
  what keeps the claim that a broadband/spectral disagreement is the physics, and never the
  interpolation, true rather than merely intended. `periodic_bilinear_interp` is on neither
  path; it is the reference `locate` is swept against.
- **`DirectionBounds::bin` lives in the header.** The observer loop calls it about 94 million
  times from `niceops.cpp`, and defined out of line in `sphere_database.cpp` it inlined only
  when `NICEOPS_LTO` was on and the toolchain supported it. It is fifteen instructions and
  should not depend on a build option. Its tangent table is an `inline const` variable rather
  than a function-local static, which would test a guard on every call. The bin geometry
  (`bin`, the tangents, `narrowest_theta_bin_degrees`, `azimuth_bin_degrees`) sits on
  `DirectionBounds` together, because the table build has to pad by the same binning the hot
  path queries with, and three independent spellings of it is how that silently stops covering.
- **The interpolation carries `{dBA, EAA}` as a two-lane vector** (`Pair` in
  `sphere_database.cpp`). Every step applies the same weight to both metrics, and a `Sample` is
  exactly a two-double vector in memory, so blending the four conditions at four corners and
  then the corners themselves is half as many instructions. Clang would not do it on its own
  from the scalar form or from a two-element array; the GCC/Clang `vector_size` spelling is
  explicit, with a plain-struct fallback that compiles to the scalar code.
- **Divides in the loop are multiplications by a stored reciprocal**: `inv_sphere_radius_`,
  `1e-3` for the EAA scaling, and one reciprocal for normalizing the propagation direction
  instead of three divides. A double divide is about four times a multiply. These moved the
  result by an ULP, the only change in that round that did.
- **`trace()` is force-inlined and `Logger::write` is `noinline`** (`src/util/logger.hpp`).
  `eval` and `level_at_observer` make about 7 trace calls each; as ordinary calls that was a
  sixth of the loop spent entering functions that immediately returned. The level lives in
  `detail::min_level`, outside the `Logger` singleton, so testing it does not touch a
  function-local static's guard.
- **The database grid fill walks channels and places them**, rather than scanning all 1,406
  channels for each of the 1,332 cells (2.5 billion `approx_equal` probes, most of the load time
  when it was written that way). Two semantics are load-bearing and preserved: a later channel
  still wins a shared cell, and `dba_max` covers matched channels only.
- **Link-time optimization is on** (`NICEOPS_LTO`, see CMakeLists.txt). The loop lives in
  `niceops.cpp` and `Database::level_at_observer` in `sphere_database.cpp`, so without it every
  one of about 100 million evaluations was an out-of-line call across a translation unit. Worth
  18% of the loop on top of the source changes, and only 5% before them, because what it buys is
  hoisting and inlining that the divides and knot lookups were hiding.

Cost on the demo case: database load about 0.61 s, broadband acoustic loop 1.55 s of CPU
(0.16 s wall on 14 threads). The first C++ build took 2.2 s to load and 11.7 s in the loop, and
the loop was 3.2 s of CPU before the work above.

Where the 19.4 ns of a broadband *evaluation* goes, from an ablation (`full`, `levels only`,
`angles only`, `geometry only`, `eval at fixed angles`): 33% the two inverse trig calls, 35%
`eval` itself (gather, blend, `log10`, arithmetic), 8% the vector geometry, 7% the SEL integral,
and about 17% the difference between a real `eval` and one at synthetic angles, which is *not*
cache misses (see the gather item under
[Measured and rejected optimizations](#measured-and-rejected-optimizations)).

## Threads

**The observer loop's chunks are dealt out on demand, not divided up front**
(`src/util/parallel.hpp`). An equal split gives every worker the same slice, so the wall time is
the slowest core's, and the cores are not interchangeable: on a machine with performance and
efficiency cores (every Apple silicon part, most recent x86 laptops) the efficiency cores'
slices become the critical path. On an M4 Pro (10 performance + 4 efficiency) the equal split
ran no faster on all 14 threads than on 10; pulling chunks lets the fast cores take more of them
and recovers about 10%. Chunk size is `count / (8 * workers)` (`chunks_per_worker`), fine-grained
enough to balance and coarse enough that the per-chunk scratch allocation stays noise. The output
does not depend on the split, since every observer writes only its own slot, so runs at any
thread count are byte-for-byte identical.

`niceops` takes its count from `-j` alone (0 means one per hardware thread). The
`NICEOPS_THREADS` environment variable is read only by the PE tools (`tool_thread_count`), which
refuse a value that is not a positive integer.

## The SEL window, the oversampling, and the windowed peak search

**The window is the outermost 10 dB-down interval.** `ten_db_down_interval()` in
`src/util/sel.hpp` runs from the first to the last sample at or above peak - 10 dB, dips below
the threshold between them included. That is the duration of the aircraft noise metrics (the
EPNL duration of 14 CFR 36 A36.4.5, and what panam's `sound_exposure_level` uses, so a
prediction and a measurement are integrated alike). It replaced the D original's walk out from
the peak, which stopped at the *first* sample below the threshold on either side (see
[HISTORY](../HISTORY.md#changes-that-affect-results)). The first-crossing walk,
`expand_window()` and its port-fidelity note are gone from the code; do not reintroduce them from
an older copy of this file. The integral is taken over *reception* time (`reception_factor`,
1 + dr/dt / c), with the step from the trajectory itself (below).

**The SEL integral is oversampled, and the discrete window is why that cannot be used.** The
sound exposure integral is a Riemann sum over a 50 Hz trajectory, oversampled by roughly thirty
times. Decimating the demo case's 6551 states moves the footprint's **mean** by 0.008 dB at
stride 2 and 0.088 dB even at stride 64. What blocks a 2–8x cut in the whole acoustic loop is the
window: decimation steps over a transient dip, the window's end jumps, and the level jumps with
it. The error does not shrink with stride the way a quadrature error would:

| stride | max cell error, first-crossing (the retired rule) | outermost-crossing (current) |
| --- | --- | --- |
| 2 | 1.93 dB (104 cells > 0.1 dB) | 1.13 dB (86) |
| 8 | 2.34 dB (619) | 1.13 dB (271) |
| 64 | 2.48 dB (3150) | 1.54 dB (1346) |

That is the signature of a discrete window flip, not of the integral. So decimation is **not**
safe as a default under either definition: the mean is negligible, but individual cells move by
1–2 dB, and the cells are the deliverable. It stays available deliberately through `adapt()`,
which is commented out in `process_trajectory`.

**The SEL integral takes its time step from the trajectory, not from `-s`.**
`sound_exposure_level` weights every sample by a time step. `-s` is the resampling step for a
WGS84 track that needs resampling (and for `adapt()`); it never coarsens a track, and an AC
track is used as sampled. Every shipped track is sampled at 0.02 s, so using the 0.5 s default as
the weight overstated every level by exactly `10*log10(0.5/0.02) = 13.98 dB`.
`process_trajectory` derives the step with `uniform_sample_interval()` and logs it, warning if
the sampling is not uniform. If you re-enable `adapt()`, it resamples to `-s` and the derived
step follows automatically.

**Skipping states the window cannot use: the windowed peak search.**
`windowed_sound_exposure_level` / `windowed_exposure` in `src/niceops.cpp` do not evaluate every
emission time. A coarse sweep (stride nt/128) seeds a best level. One pass then evaluates every
state `DirectionBounds` cannot prove is below that best level, handing the same geometry on to
the evaluation so a state that survives the test pays nothing for it. `ten_db_down_interval`
finds the window's ends by scanning in from both sides, skipping a state the bound puts below the
threshold (`certainly_below`) and evaluating everything between the ends. It is **exact**: the
output is byte-identical to evaluating every state, on every track, grid, observer set, receiver
and thread count tried, because the integral and the window are the shared
`ten_db_down_interval()` and `sel_over_window()`, templated on a level accessor (a filled array
for `peak_range`, a memoizing lambda here), not a second copy.

It runs on the **spectral** path only: `--spectral` (the default whenever the database carries
spectra), `--ground_effect`, plates, the propagation table (whose largest value per band enters
the bound), and `--refraction` (the bound is taken along each ray's launch direction at its ray
tube's range; for a sum of stratified arrivals, the loudest arrival's direction over the harmonic
sum of their squared ranges, `bound_geometry`). It is skipped, and the loop marches every state,
for:

- `-a` (the level is transformed after eval);
- `-i` and `--band_history` (the whole history is written);
- `--epnl_file` and `--band_frames` (every state's bands go into the frames, and the EPNL's
  window is found on PNLT, which the A-weighted bound says nothing about);
- `-D`, `--lden_file` and `--cnel_file` (the window is found on a penalized history, which a
  window found on the raw one need not cover);
- `--lamax_file`, `--nax_file` and `--ta_file` (the time weighting runs over the whole history);
- `--csel_file` (the C-weighted window is not the A-weighted one, and every state's bands are
  needed to find it);
- `--full_history`;
- `--db_cutoff` (which culls on its own terms);
- any database with a negative excess-attenuation factor, where dropping that term is not a
  loss.

What it saves depends on how dear an evaluation is. Under the outermost-crossing rule it
evaluates most of what it is shown: on the broadband demo case it evaluated 93.9% of the pairs
(the window itself is 59%, and the bound proves only about 6% of the pairs outside it) and ran
about 4.6% *slower* than marching every state (measured twice, output identical). So **broadband
runs do not build the table or use the search**. `Database::direction_bounds`, the broadband
builder, remains, tested, for a broadband search that pays. A spectral evaluation costs about
150 ns and a ground-effect one up to a microsecond, and there the search pays. The main pass
keeps each skipped state's numerator and squared range (`skipped_numerator`,
`skipped_squared_range`) so that `certainly_below` does not recompute the geometry (a full ray
solve under refraction) for the 30 to 40% of states it is asked about; the answers are the same
two numbers compared the same way, so the output is unchanged.

Three things had to hold, and each was measured before any of it was written:

- **The bound needs less than 10 dB of slack**, or the peak search reaches outside the window
  and there is nothing left to save. The whole-sphere maximum is 11.4 dB, just over. Binning it
  by direction (`DirectionBounds`, bins uniform in cos^2(theta) crossed with azimuth sectors)
  brings the median to 8.6 dB.
- **The test must cost nothing on the states it does not reject.** It forms the squared range
  and the three dot products with the sphere frame, which is exactly what an evaluation needs,
  so `Database::level_from_geometry` takes those and a surviving state pays nothing for having
  been tested. Letting the test compute its own geometry and then recomputing it inside
  `level_at_observer` was measured: it turns a 1.48x into a 1.08x.
- **The table has to stay in cache.** One table per *distinct blend*, not per state: a shipped
  track resolves thousands of states onto a few dozen tetrahedra, so it is tens of kilobytes.
  Built per state it was 3.3 MB, re-streamed for every observer, and that alone cost more than
  the scheme saved once the loop was threaded.

The bound maxes over the four conditions of a tetrahedron rather than over the particular blend
of them, so it is looser than a per-state table would be; buying that back means a table per
state, which is the cache problem above. Bin counts were swept under the first-crossing window
rule (8x8 to 32x16) and the choice hardly mattered; 16x8 is kept.

The whole thing rests on one property, that the bound is never exceeded, and it is tested
directly (`tests/test_sphere_database.cpp`, `tests/test_spectral_propagation.cpp`, tag
`[bounds]`) over sweeps of source states and observers, **including sweeps aimed within five
degrees of the nose and tail axes**. A cell is filed under the bin of its center direction, and
at theta = 0 or 180 that direction has no azimuth, so a pole cell is filed under every azimuth bin
of its theta bin (`Database::home_bins`). Without that, a query a degree off the axis was bounded
up to 0.8 dB too low, and the ordinary sweep, which samples hardly any observers there, never saw
it.

**`--db_cutoff` on the spectral path** bounds the band-limited, weighted level per direction
(`Database::spectral_direction_bounds`). With `--ground_effect` or `--plate_table`, each band is
raised by the most the ground or the plate can add.

## The spectral and ground paths

Current cost of the acoustic loop on the full demo grid (S-76D track, B407 fixed-load-factor
database, 94.3 M evaluations, 14 threads; the broadband demo case, on `S-76D_M3.nod`, is the
0.16 s wall above):

| | acoustic loop |
| --- | --- |
| spectral, no ground | 1.6 s |
| `--ground_effect` | 2.1 s (7.1 s with `--exact_ground`) |
| `--ground_effect --plate_table` | 5.4 s (13.3 s with `--exact_ground`) |

**By default both receivers interpolate per-trajectory tables** instead of evaluating the model
at every observer and state: the point receiver its whole ground correction over source
elevation and log range (`GroundTable`), the plate receiver its reflection coefficients. The
tables are built exactly at their nodes in well under 0.1 s. Their accuracy is in
[FIDELITY](../FIDELITY.md#ground-and-plate-tables). The interpolated values never exceed the
model's bounds, so the peak search and `--db_cutoff` stay exact; the plate is bounded by
`PlateReceiver::gain_bounds`, and `plate_windowed_is_exact` checks that skipped states change
nothing. `--exact_ground` evaluates the model directly.

**The exact path was profiled, not guessed at.** Two rounds of intuition (cache layout, then
`pow`) bought 14% between them; `sample` on a running process then showed the real costs, and
the same effort took the ground path (B407 database, default band limits) from 41.4 s to 10.9 s,
and the Faddeeva Horner split to about 9.4 s. At that point carrying the spectrum cost about 4x
the broadband loop and the ground effect about 32x: a complex reflection coefficient and a
Faddeeva evaluation per band per observer per state, none of which can be hoisted because it all
depends on the geometry. What it found, and what not to undo:

- `std::sqrt` on a `std::complex` costs an **atan2**, because libc++ evaluates it in polar form.
  `ground_effect.cpp` called it three times per band. Two were provably redundant (`sqrt(w)`
  *is* the numerical-distance parameter wherever its real part is positive, which is every
  passive ground but a mass-like one; see [ground.md](ground.md)), and the third is an algebraic
  root written out in the band loop. `std::abs`/`std::arg` on a complex cost a hypot and an
  atan2 and are gone the same way: `|Q|cos(p + arg Q)` is `Re(Q)cos p - Im(Q)sin p`, and
  `|w| <= 500` is `norm(parameter) <= 500`.
- `faddeeva` is written in real arithmetic. Every `std::complex` multiply carries a branch to
  recover NaN/infinity cases that cannot arise here, and the Horner loop does 31 of them. Its two
  trailing divides are folded into one: `2p/d^2 + c/d == (2p + c d)/d^2`.
- **The Weideman series runs as four interleaved Horner chains over t⁴** rather than one over t
  (`src/faddeeva.cpp`). It was 53% of the ground path on its own, and 31 dependent complex
  multiply-adds are a latency chain, not a throughput problem: it ran at about a quarter of the
  machine's floating-point rate. Chains a quarter as deep issue together: 32% off `faddeeva`,
  14% off the whole ground path. The coefficient table is left exactly as Weideman's
  construction produced it and indexed by stride, so it stays checkable against the source.
- `spectra_` is condition-innermost, matching `grid`. The netCDF file stores it the other way
  and the loader transposes once.
- The whole band loop in `excess_ground_attenuation` is real arithmetic for the same reason, and
  its band-invariant branches (boundary loss on/off, the sinc, turbulence) are hoisted so the
  body stays straight-line.
- Bands outside the metric's limits are not computed at all. `first_band`/`last_band` bound a
  contiguous run, which is always what the limits select because band centers are sorted.

All of it is exact: the footprint moved by 2.8e-14 dB end to end, and the EGA tests still pin to
panam at 1e-9 dB.

**A synthetic spectral database for timing.** `write_flat_spectral_database` in
`tests/test_spectral_propagation.cpp` interpolates a sphere out of the shipped `S-76D_M3.nod`,
hangs synthetic band amplitudes on it and writes it back with `write_sphere`, giving a one-group
`.nod` that `--spectral` accepts, so the path runs from a bare checkout. Its single condition
makes the blend degenerate, so timings off it run optimistic against a real multi-condition
database, but the ground term, which is most of the cost, does not care. Over it, the September
2026 round took broadband 0.32 -> 0.16 s, `--spectral` 1.36 -> 1.16 s and
`--spectral --ground_effect` 8.95 -> 7.81 s, all three agreeing with the previous build to
1.4e-12 dB. For the ground term alone, driving `excess_ground_attenuation` over a spread of
source heights and ground distances with 26 synthetic third-octave bands reproduces the per-call
cost the observer loop sees (about 840 ns, from about 975), and `tests/test_ground_effect.cpp`,
pinned to panam at 1e-9 dB, is what says it is still right.

**The asymptotic boundary loss was tried and rejected.** For large `|w|` the term tends to
`-(1/2w + 3/4w^2 + ...)`, which would skip `faddeeva` on about 20% of bands. It does not converge
to the exact value: the error plateaus at 1.6e-5 at `|w| >= 50` and 1.9e-6 at `|w| >= 100` no
matter how many terms are taken, so it cannot hold the 1e-9 dB pinning against panam at any
cutoff in range. Roughly 9% of the ground path, for the loss of the property that proves the
port. Not worth it; do not re-derive it.

### Terrain, ground classes and screening

Ground classes (the maximum over the classes present) and screening (`screen_gain_bound`) each
add a term to the peak search's bound, so each stays exact and evaluates a little more; a table
read over terrain adds none, since its term is the table's largest value whatever the geometry.
`--terrain_screening kirchhoff` adds none: its loss is at least 10 lg 4, so a blocked pair lies
below its unscreened level, and behind the ridge it evaluates what the unscreened run does.
None adds a march condition. The terrain runs (tables over terrain, screening) are timed
in [ground.md](ground.md#how-the-plane-is-found-groundplanecache). What the bound costs there,
in evaluated fraction on the full demo grid: a table over the hills nothing (98.89% against
98.90%), screening 0.7 points on the hills, where nothing is screened, and 0.8 behind the ridge.

**Ground classes**, on 14,400 1.2 m pole microphones at the demo grid's points (an observer file,
since grid observers have no class), variable-porosity grass of 200 kPa·s/m², and Delany–Bazley
concrete of 10,000 at every other microphone. 8 threads on an M4 Pro, asserts on, best of three,
on a shared machine (load 10 to 31):

| run | CPU | acoustic loop | evaluated | peak RSS |
| --- | ---: | ---: | ---: | ---: |
| all grass | 24.9 s | 3.10 s | 98.48% | 231 MiB |
| half concrete | 24.9 s | 3.10 s | 98.57% | 242 MiB |
| all grass, `--exact_ground` | 96.2 s | 12.0 s | 98.48% | 221 MiB |
| half concrete, `--exact_ground` | 96.2 s | 12.0 s | 98.57% | 221 MiB |

- **The per-observer admittance costs nothing measurable**: the loop time does not move.
- **The bound over both classes** evaluates 0.09 points more, which is why one
  `DirectionBounds` table serves every class ([ground.md](ground.md#receivers)).
- **Each further table** is 641 angles × 147 ranges × 28 bands of floats, 10.6 MB, built in
  0.07 s; peak RSS rose by 10.3 MiB. The run warns above 16 tables.
- Grass microphones give the same bits in both runs, tabled or exact, as the windowed search's
  exactness requires. The tables' accuracy is in
  [FIDELITY](../FIDELITY.md#ground-classes).

## Measured and rejected optimizations

### Rejected

The numbers are here so they are not re-derived.

- **Hand-rolled `acos`, `atan` and `log10`.** They are a third of the loop, which makes them
  look like the obvious target, but this platform's libm does `acos` in 2.3 ns, `atan` in 2.0
  and `log10` in 1.9: around 8 cycles each, table-driven. A polynomial accurate enough to hold
  the golden (about 1e-12 rad, so the level moves by less than 1e-10 dB) needs a square root and
  a degree-10 minimax and is *slower*. One cheap enough to win (degree 4, about 1e-5 rad, about
  1e-3 dB on the footprint) saves only about 1.5 ns of the 3.55 an in-situ call costs. Bad trade
  in both directions.
- **Per-cell polynomial inverses for `acos`/`atan`.** The angles are never wanted for their own
  sake, only a grid cell and a fraction across it, so a polynomial per 5-degree cell mapping
  cos(theta) straight to the fraction would skip the inverse trig entirely. It cannot work: acos
  has a square-root singularity at the poles, so no polynomial in cos(theta) represents it.
  Degree 7 still leaves 0.31 degrees of error in the two polar cells, and 5e-5 degrees in the
  best interior ones; changing variables to sqrt(1 - a) to fix it puts a square root back.
- **`-fno-math-errno`.** The hope was that marking the libm calls pure would let the compiler
  overlap them. No measurable change: it will not software-pipeline across an opaque call, and
  the three calls in an evaluation are already independent of each other and so already overlap.
- **Storing the grid as `float`.** It would halve the 28 MB gather. It costs about 1e-5 dB,
  which is within engineering accuracy but ends the golden comparison as a bit-level regression
  test. Not taken without being asked, and see the next item, which says it would not help
  anyway.
- **Anything aimed at the gather's locality.** The four blended conditions of a cell span 845 of
  the 1312 indices on average and touch 3.1 cache lines per corner, never 1, so an evaluation
  pulls about 800 bytes of lines to use 256, which looks like an obvious 3x to reclaim. It is
  not. The whole trajectory uses only **52 distinct blends**, so the four needed condition slices
  could be compacted into 4.4 MB of contiguous 64-byte rows. Before building it, the ceiling was
  measured the cheap way, running the same loop against a one-condition database whose 21 KB
  grid makes the gather perfectly local: 1.82 s against the full database's 1.80 s. **The loop
  is compute-bound, not miss-bound.** Measure that ceiling before touching the layout again.
- **The per-state slab; do not rebuild it.** Collapsing the database to one blended 36×37 slab
  per trajectory state (with the observer loop blocked so `state` can be the outer loop) works
  and is bit-exact, but it buys only 6.8% of the acoustic loop for +131 MB of resident memory.
  The condition-major grid layout it requires makes the ordinary blend path (`self_test`, and
  any observer set below the break-even block size) about 6% *slower*. Block size barely moves
  it, and a single 768 MB block is the worst configuration of all, because the level buffer's own
  bandwidth then dominates. The result is recorded in commit c720d23. The intuition fails
  because observers sweep a 120×120 ground grid, so consecutive evaluations land on neighboring
  (phi, θ) cells and the gather is already cache-resident. Pointing `eval` at a fixed
  L1-resident slab and timing it costs ten minutes and would have shown this before any
  restructuring.
- **Culling states the SEL window will discard, by bounding them one at a time.** 44% of
  (observer, state) pairs fall outside the window, but a bound that proves a *particular* pair is
  outside it has to be tighter than 10 dB, and `level <= sphere max - spreading` is 11.4 dB loose
  at the median: it proves only 7.6%. Asking the question the other way round is what works:
  the windowed peak search above.
- **Storing reciprocals of the (max − min) normalizing ranges.** A reciprocal multiply is not
  the same as a divide, so it would move the goldens, and the normalization is per state anyway.
- **Hoisting the output-directory check out of `write_time_history_file`.** One `exists()` next
  to creating and writing a ~130 KB file is noise, and it would mean bypassing `open_output`'s
  guarantee.
- **A source/receiver-height cache for propagation tables.** Caching source and receiver height
  cells in the caller, together with skipping bearings for isotropic tables, was tried after the
  four-band lookup below. Five alternating GCC runs on B407/S-76D trajectory 197421AC with
  `tables/ground_model_test.nc`, `--gridstep 328 -j 1 --humidity 20`, gave 2.953 versus 2.952 s
  median acoustic time (6,208,943 evaluations). Output was byte-identical, but the extra cache
  and public API did not earn their complexity, and were removed.
- **Sharing PE setup across turbulence realizations.** A realization repeats the vertical
  operator, the Padé fit, the factorizations and the receiver stencils. Timed as a 15 m march
  against a 1.2 km one (a receiver march at gamma_T 3e-5, L0 10 m), setup is 23% of a march at
  50 Hz but 4–5% above 1 kHz, where a table's time is spent: 6.0% over the 50 Hz–4 kHz
  third-octave bands. That share also includes the starter, one step and the per-realization
  screen source, which cannot be shared, so it is a ceiling. Not worth separating march closures
  that carry mutable scratch.
- **Reusing the straight-path range in the full-history loop.** The reception-time range and
  `level_at_observer` each form `observer - position` and its norm. Replacing the loop's range
  with a constant, which saves more than sharing it could, saved 0.05 s of 3.14 (1.6%) on the
  broadband demo and nothing measurable on a spectral `--full_history` run: LTO already merges
  most of it.

### The residue-series tables

`--refraction` with a ground model builds V, the residue series across each band's shadow
boundary (`ShadowTransition`), and `--shadow_model resolved` builds the cramped bands' table
(`set_cramped_ground`), both from `upward_pressure` and `upward_pressure_heights` in
`exact_solutions.cpp`. In the validation harness's wide study the cramped table was 75% of all
NICEOPS CPU (884 resolved runs at 74 s each, against 0.3–4 s for every other option), and V
about two thirds of every other refracting run. Profiled with `sample`, then taken, in order:

| change | cramped table, 24 bands, `-j 1` |
| --- | ---: |
| before | 47.9 s |
| each range summed for its own source height only (the table reads one of six) | 32.6 s |
| a range's terms skipped once they have stayed under 1e-17 of its largest sixteen in a row, the last sixteen still taken for the tail test | 11.4 s |
| `airy` near the negative real axis from one series: the connection formula's two expansions share it, since zeta(w^2 z) = -zeta(w z) | 8.6 s |
| real arithmetic for the mode weights' divides and `upper`'s square root, magnitudes as sqrt(re^2 + im^2) (std::abs outside 1e-290..1e290), Newton's divide; Ai(b) reused for a source on the ground | 6.5 s |

V for a flush plate (28 bands, `-j 1`, 329 m/s) went from 2.38 to 0.95 s and for a 1.2 m pole
from 1.30 to 0.62 s with the same changes to `modal_sum`. Whole harness runs: an AS350 approach
with `--shadow_model resolved` from 59 to 11 s at `-j 1` (setup 1.2 s at `-j 8`, from 8.1), and
the same with `--refraction stratified` alone from 4.5 to 2.8 s.

What it moved, measured: both runs' SEL, LAmax and every time history are byte-identical. In the
stored (float) cramped table 5 of 112,606 values moved, by at most 7.6e-6 dB, and the NaN pattern
is unchanged. V agrees within 2e-6 dB from xi = -1 on; its deep lit side moved up to 0.4 dB at
the nodes that are rounding noise anyway ([refraction.md](refraction.md#shadowtransition-the-level-near-a-shadow-boundary)),
as far as clang and GCC builds of the old code differ there. The skip rests on the terms' size
being unimodal in the mode number at a fixed range (|weight| grows as e^{c y n^{1/3}}, H0
decays as e^{-r Im kappa}, Im kappa growing as n^{2/3}), the same assumption that sets the
mode count. Guarded by `near the negative real axis Ai is the connection formula's two
expansions` and `the residue series over several source heights is each height's own`
(`tests/test_pe_shadow.cpp`), and by every resolved and stratified CTest.

Not taken:

- **A lower mode cap.** 566 of 576 cramped calls reach the 12,000-mode cap, driven by a
  near-source range that ends NaN. But 648 trusted ranges never go quiet within it and trusted
  ranges converge as late as mode 11,968, so their values and trust depend on the full count.
- **The zeros from their neighbors.** Seeding Newton from 2 t_{n-1} - t_{n-2} would skip
  `seeded_zero`'s fixed-point iteration, about 13% of what is left, but moves every zero's last
  bits, which reach V's lit side as above.
- **The PE.** A propagation table's build is 96% `pade_step`, already interleaved and fused
  (pe.md); the factorization's complex divides are the next 4%. The profile PE's, 8.6% of a
  terraced run, are written out (ground.md, the row of roofs).

### The propagation models, measured as a suite

Twenty-two cases, every propagation model at both of its shapes of use: the validation
harness's (one AS350 approach at Amedee, 46 flush plates, `-j 1`, histories written) and a
footprint's (the demo track on a 41 x 41 grid, an eighth of the demo grid, `-j 8`), plus one PE
table build. Each case was run with the build before and after, alternately, on the same machine
at the same load, and every output compared. The harness's own record of what its CPU went on
set the order: in the wide study `--shadow_model resolved` was 75% of all NICEOPS CPU, and
every refracting option's setup most of the rest.

| case | before, CPU | after | × | outputs |
| --- | ---: | ---: | ---: | --- |
| harness, `--shadow_model resolved` | 11.0 s | 1.38 s | 7.96 | identical |
| harness, `--refraction linear` | 1.42 s | 0.42 s | 3.35 | identical |
| harness, `--refraction stratified` | 2.54 s | 1.37 s | 1.86 | identical |
| harness, plate / soft ground, no refraction | 0.47 / 0.41 s | 0.41 / 0.35 s | 1.15 | identical |
| grid, duct with `--tunneling --skip_fill` | 1779 s | 586 s | 3.04 | within 3.5e-12 dB |
| grid, `--shadow_model resolved` | 41.4 s | 31.7 s | 1.31 | identical |
| grid, stratified in the duct | 94.2 s | 75.5 s | 1.25 | identical |
| grid, `--scattering` | 10.9 s | 9.9 s | 1.10 | identical |
| grid, stratified windy (point, plate) | 35.1, 28.7 s | 35.3, 27.6 s | 1.0 | identical |
| grid, ground, exact ground, plate, linear, table, terrain, screening, broadband | | | 1.0 | identical |
| PE table (8 bearings, 17 bands, 5 frequencies) | 116 s | 116 s | 1.0 | identical |

What did it, each bit-identical but the first:

- **The wave arrivals' band mean of Ai^2 is a table** in the band center's argument
  ([refraction.md](refraction.md#wave-arrivals---tunneling-and---skip_fill)): seven complex Airy
  evaluations per band per arrival had been half of a duct run with both options.
- **The shadow-boundary tables are built when first wanted**
  ([refraction.md](refraction.md#shadowtransition-the-level-near-a-shadow-boundary)): V whole,
  the first time a ray nears a boundary or the bound asks; the cramped series a (band, s)
  column at a time. A validation run marches every state, so builds no bound, and on the
  approach measured no ray came near a boundary: it built neither, and 8.5 s of the resolved
  run's 11 and 0.9 s of every other refracting run's went. Where they are built, the values are
  the same. `apply` fetches the tables once a pair: fetched per band (twice, through `level`
  and `plate_offset`) the `call_once` check cost a linear run 10%. And it places each band's s
  among the table's nodes once (`s_place`), for V and its floor alike, with the per-band helpers
  (`level_at`, `finish_level`, `plate_offset_in`) inline: a second log10 a band and the calls
  out of line cost a stratified grid run 5–10% of its CPU.
- **The layer walk** (`Path::between`, stratified.cpp) forms each level's sine once, as the top
  of one layer and the bottom of the next, and looks up its ends' level and speed once a pair
  (`Path::fixed`: the source, receiver and ground heights it walks between for every p, looked up
  when the path is formed) rather than by a binary search and an interpolation at both ends of
  every walk. 10–20% in a duct, where walks are long; nothing measurable on the windy profile's
  17 levels. A turning height moves with p, so it is looked up when asked, not kept.
- **`--scattering`'s level is summed from parts formed once**: the run's and band's in
  `set_scattering`, the pair's range and height before the band loop, in the original
  expression's order so its multiply-adds fuse the same way; and 10 lg of V's largest power
  ratio once per band. `ShadowScattering::level` is the same parts summed, so the two cannot
  drift apart.
- **Time and band histories are formatted by `std::to_chars`** (`append_number`), the
  standard's definition of the %.17g (and %.10g) the streams wrote, with every NaN written
  "nan" as the streams do. Writing 46 histories of 4780 rows was half of a plate run's loop.

Tried and not taken:

- **A warm start for the stratified eigenray solve from the previous state's root**: 2–4% on
  the windy grid, with outputs moved up to 1.6e-6 dB, and the result would depend on the order
  in which the windowed search evaluates states, so windowed and marched runs would stop
  agreeing to the bit. The linear model's guess is already close.
- **Secant-guided bracketing and Anderson–Björck in place of Illinois** in `solve`: 1.5%, and
  outputs moved up to 1e-4 dB: some pairs are sensitive to which root within the tolerance
  comes back.
- **A tighter bound for the windowed search.** Every grid case evaluates 97–99% of its pairs;
  the outermost 10 dB-down window covers most of the demo track for most observers, and even
  the broadband sphere bound alone (no ground term) evaluates 94%.
- **The PE table build** is 96% `pade_step`, already interleaved and fused (pe.md): cost goes
  as f^2 through the grid and the step, both fixed by accuracy.
- **The stratified solve itself** is about 11 range evaluations a solve (2 at the family's
  ends, 3.5 bracketing the guess, 5.3 iterating) and 30 a pair, each a walk whose inner loop is
  a square root, a divide and a dozen other instructions a layer. Batching two p values a walk
  across a vector's lanes is what is left, for the ends and the tube's two points at most.

### Interpolated rays

`--ray_stride` (4 by default; [refraction.md](refraction.md#interpolated-states---ray_stride)) on
the same suite, CPU against every state solved, same build, alternately:

| case | stride 1 | stride 4 | stride 8 |
| --- | ---: | ---: | ---: |
| harness, `--refraction stratified` (46 plates, `-j 1`) | 1.37 s | 0.77 s | 0.77 s |
| harness, `--shadow_model resolved` | 1.38 s | | 0.77 s |
| grid, `example_windy.csv`, point receivers | 35 s | 17.6 s | 16.1 s |
| grid, windy, flush plates | 27.6 s | | 18.3 s |
| grid, windy, `--shadow_model resolved` | 31.7 s | | 15.8 s |
| grid, `example_inversion.csv` duct | 75.5 s | 29.9 s | 27.7 s |
| grid, duct, `--tunneling --skip_fill` | 586 s | | 542 s |

Past stride 4 the solves are a small part and the rest of the evaluation is what is left.

### Taken, each exact

Each keeps the formulas and the per-band summation order; none involves an approximation,
reduced precision, a change of table resolution or a change to which states are evaluated.

- **`AbsorptionTerms`** computes the classical loss, two exponential numerators and the
  temperature power once per spectrum, alongside the relaxation frequencies, which cuts repeated
  transcendental work when building layered absorption profiles. Its evaluation method stays out
  of line so GCC uses the same compiled expression for scalar and batch calls; inlined, they
  contract differently and break their exact-equality test.
- **`PropagationTable::lookup` accumulates four bands together**, sharing each corner's address
  and weight. A scalar tail handles the remaining bands. Each band's sum still starts at zero and
  visits corners in the original order, so nothing is reassociated. Zero-weight corners remain
  unread, and excluded output bands untouched.
- **`AtmosphereProfile::with_hydrostatic_pressure` keeps a running integral**, so N levels cost
  O(N) segment evaluations instead of O(N²), with the same addition order and ground reference.
  The endpoint interpolation is kept separate from the complete-segment prefix so that even its
  rounding stays the same.
- **`AtmosphereProfile::at` locates the height cell once** for temperature, humidity and
  pressure, and `wind_at` shares one search between its two components. Neither interpolation
  nor endpoint holding changes.
- **LAmax (or NAx) with TA in one pass.** `weighted_maximum` and `time_above` each converted
  every sample to a mean square and ran the detector (or built the running energy and walked the
  window pieces) over the same history. `weighted_maximum_and_time_above`
  (`util/time_weighting.hpp`) does both in one pass with the same arithmetic, and is the only
  implementation: `weighted_maximum` and `time_above` call it (a threshold of +infinity is never
  exceeded, so LAmax alone costs what it did, within 0.5%), which keeps each detector and
  crossing formula in one place.

Measured medians (Apple Clang `-O3 -std=c++17` microbenchmarks without LTO, seven alternating
before/after runs, except the last row: clang Release with LTO, `-j 1`, five alternating runs):

| Work | Before | After |
| --- | ---: | ---: |
| 2 million 24-band table lookups (synthetic 4 × 3 × 3 × 4 table) | 0.454 s | 0.216 s |
| 200,000 24-band absorption spectra (280–299.9 K, from 50 Hz) | 0.0881 s | 0.0225 s |
| 20 hydrostatic pressure reconstructions, 2,000-level profile | 95.69 ms | 0.390 ms |
| 2 million paired atmosphere/wind queries, same profile | 142.34 ms | 78.32 ms |
| demo case, `--lamax_file --ta_file`, fast / `leq0.5` | 3.86 / 5.23 s | 3.54 / 4.28 s |

These are routine timings, not whole-run speedups; small atmosphere profiles have much less to
save. Exactness was checked directly: over 10,000 table queries (240,000 band values) and 64
atmospheres (1,536 coefficients), Clang matched bit for bit throughout, and GCC 16 matched the
table values bit for bit while absorption differed by at most 5.6e-17 dB/m (compiler rounding).
The tests keep exact scalar/batch agreement, check the absorption formula across temperature,
pressure and humidity, and check the pressure walk and interpolation against the original
formulas (singleton profiles, uneven and nearly isothermal layers, surface-reference offsets)
under both compilers. Table tests sweep even and uneven grids, reordered columns, partial band
ranges and NaN gaps. Every LAmax and TA output is byte-identical under all four weightings.
