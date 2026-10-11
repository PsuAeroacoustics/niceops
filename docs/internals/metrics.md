# Output metrics: EPNL, band frames, time weighting, day metrics, C-SEL

The user-facing description is in the manual's
[EPNL and reception-time frames](../MANUAL.md#epnl-and-reception-time-frames),
[Maximum level, time above and C-weighted SEL][max-level] and
[Several trajectories, and the day-night level][day-night].
The equations are in [THEORY.md](../THEORY.md) §10. Comparisons with the 2017 measurements are in
[FIDELITY.md](../FIDELITY.md#lamax-ta-and-epnl-against-the-2017-measurements).

## EPNL and band frames

`src/pnl.{hpp,cpp}` and `src/band_frames.{hpp,cpp}` implement `--epnl_file` and `--band_frames`:
EPNL (14 CFR 36 A36.4) from one-third-octave levels in 0.5 s reception-time frames, and the
frames themselves as an (observer, frame, band) netCDF file.

### Things to preserve

- **`pnl.cpp` is panam's `flight_acoustics.py` (eeb040b with F005, F014 and F026), term for
  term**: `noys`, `perceived_noise_level`, `tone_correction` and `effective_perceived_noise_level`.
  That includes panam's reading of the regulation where readings differ:
  - 500 Hz is in Table A36-2's *middle* range.
  - The duration's limits start at the first and last records with PNLT ≥ PNLTM − 10 (the longest
    duration when there are several peaks), and each moves one record outward when that record is
    strictly closer to PNLTM − 10 (A36.4.5.5). Dips inside are included.
  - D takes A36.4.5.4's −13 for 0.5 s records (not −13.0103) and 10 lg(Δt/10 s) for any other step.
  - The A36.4.4.2 band-sharing check is applied, with the five-record mean truncated at the
    history's ends. Its Δ_B is added to EPNL as a term of its own; PNLTM, the limits and D stay the
    unadjusted history's.
  - A −∞ band enters the tone correction at SPL(d) with no tone of its own. A NaN band propagates,
    and is refused in an EPNL history.

  panam's comments say the MATLAB `EPNLcalc.m` differs on the first three on purpose; do not "fix"
  the port toward it. A change here is a change to panam first.
- **The reference values in `tests/test_pnl.cpp` are generated, not typed**:
  `PANAM=... python tools/pnl_reference.py` prints the block between the `clang-format off/on`
  lines. Regenerate rather than edit. Two spectra exist only for mutants that every other case let
  through: a 10 kHz tone *on a slope* (on a flat spectrum s(23) = 0, so dropping it from
  SPL'(24) = SPL(23) + s(23) is invisible), and a rise that levels off to a slope of exactly 0 (the
  only place step 3(b)'s "zero or negative" differs from "negative"; a plateau after a sharp rise
  is encircled by 3(a) anyway). The inputs are 0.1 dB literals and the histories
  `base + envelope (+ add)`, so both sides parse and add the same doubles.
- **`--epnl_masking` (default on) is not optional for a prediction** (panam's `masked=True`, F005).
  A band below its noy threshold SPL(d) enters the tone correction's steps 1–7 at SPL(d) and
  carries no tone itself (`ToneMasking::below_noy_floor`; the library's default stays `none`,
  panam's). The literal procedure reads absorption's plunge in the top bands, and an inaudible band
  rising out of it, as tones of up to 20/3 dB, and EPNL reads high (how far, on the demo track:
  [FIDELITY](../FIDELITY.md#lamax-ta-and-epnl-against-the-2017-measurements)). A measurement
  never sees this (its noise floor), and certification excludes tones in masked bands (FAA AC
  36-4 App. 2).
  - **Both halves are needed.** Exclusion alone leaves an audible band before a plunge with a
    crashed background (20/3 dB). Reconstruction alone lets masked 8 and 10 kHz bands, rising 17
    and 21 dB along SPL(d), read as a tone. The reference table has a spectrum for each, and
    dropping either half fails it.
  - **Masking can raise an EPNL a little** (it did at one observer of the demo grid: the
    band-sharing check and the duration pick other records), so `epnl_masking_lowers_tones`
    asserts only that some EPNL falls, never that none rises.
- **A frame is an energy integral, not a mean of the states in it.** Each state's band energy is
  spread over [t_r − w/2, t_r + w/2], t_r = t + r/c and w = T × `reception_factor` (the SEL
  integral's own weight), and split between frames by overlap. That makes the frames sum to the
  whole-history exposure exactly. Frames sit at multiples of 0.5 s on the trajectory's clock, so a
  grid's observers share one axis.
- **A missing PNL band is refused, never filled.** The databases' bands are matched to the nominal
  centers within 5 percent (`pnl_band_indices`); a database without one, or `--band_limits`
  leaving one out, stops the trajectory with the bands named. Filling (the harness's 5 kHz minus
  5/10/15 dB for measured frames) belongs to the comparison, not the tool.
- **Without `--band_history`, `record_path` returns the band buffer and skips the geometry**
  (`lookup_angles`, an asin): it runs per state per observer. The accumulation uses
  exp(L ln10/10), as `weighted_level` does.
- **`--band_frames` holds every observer's frames, as float, until the write.** Each observer's
  record is one chunk, shuffled and deflated at level 1: +2.2 s and 28 percent smaller on the demo
  grid. `nc_close`'s status is checked, as the sphere writers check it.
- **`EpnlRow` keeps scalars, not `EpnlResult`'s per-record histories**, which on a grid were four
  vectors per observer kept only to read one element.

### Cost

Both options force every state to be evaluated (the windowed search bounds the A-weighted level,
which says nothing about where PNLT peaks), and hold every band of every state of an observer
while its frames are formed. On the demo grid with B407 (spectral, ground effect, 14 threads):

| | run time |
| --- | --- |
| neither | 4.2 s |
| `--epnl_file` | 5.6 s |
| `--epnl_file --band_frames` | 10.1 s, of which 2.2 s is the compression |

`--band_frames` holds every observer's frames until the end: about 450 MB of data, 0.9 GB
resident, and a 318 MB file (14,400 observers, 273 frames, 28 bands).

### How it was checked

The tests that guard EPNL and the frames (against panam's own functions at 1e-6 dB, the masking
cases, the frames' energy and timing, and the end-to-end `epnl_*` CTests) are described once, in
[testing.md](testing.md#epnl-and-band-frames); the mutants they were written for are in its
[table](testing.md#mutation-guarded-tests). The one-off independent rebuild of the frames is in
[FIDELITY](../FIDELITY.md#lamax-ta-and-epnl-against-the-2017-measurements).

## Source provenance

`--provenance_file` (`write_provenance` and the observer loop's provenance block in
`niceops.cpp`). Things to keep:

- **Each state is weighted as the SEL integrated it** (`sel_weight`: the 10 dB-down window, or
  the whole history, times the step and the reception factor), as is `--coverage_file`. Summing
  every sample alike counted states the SEL leaves out, and across operations of different steps
  weighted each by its sample count. `exposure_db` reproduces the SEL to 5e-9 dB.
- **Energy is attributed by nominal blend weight.** The blend is linear in dB, so a sphere's weight
  is not its share of the energy exactly, and a sphere empty toward the observer, which the
  evaluation drops, is still credited. It shares out all of the energy (the writer asserts the
  sums agree to 1e-9).
- **Per-state flags are taken where the run decides them**: the axis clamps from the blend, the
  load-factor holds where the run applies `--max_load_factor` and `--min_load_factor`.
- **Per-sphere sums are dense per thread and sparse per observer.** A dense observer-by-sphere
  table is 150 MB on the demo grid with a 1312-sphere database; the sparse list holds the few
  dozen spheres an observer hears. Only touched entries are reset.
- **It marches** (bounding off), since every state's energy is needed; the SEL is the windowed
  run's (`provenance_leaves_sel_alone`).
- **One database per run**: a sphere index means nothing across two.

## Contours

`src/contour.{hpp,cpp}` (`--contour_file`). Things to keep:

- **The grid is padded by one node on every side, each pad at its edge neighbor's coordinate**
  and below every level. The crossing toward a pad then sits on the edge node, so a region
  reaching the edge closes along it and every ring closes; the zero-width pad cells only repeat
  points, which are dropped.
- **Segments carry the high region on their left**, so rings around high ground run
  counterclockwise and holes clockwise, as RFC 7946 wants, and an outer ring is told from a hole
  by its signed area. A rotation (`--frame_bearing`) keeps orientation; only a reflection would not.
- **A saddle is decided by the mean of its four corners.** A non-finite node (culled) is below
  every level, and a crossing against one sits on the finite node.
- **A hole goes to the smallest outer ring holding a point strictly inside it**: a hair to the
  right of the midpoint of its longest edge (holes run clockwise). Rings of one level never
  cross, so that ring is its own. Its vertex mean was tried and failed: it lies outside a U-shaped
  hole, which was dropped (area 33% high), and inside an island standing in the hole, which took
  the hole and came out with negative area.
- **Rings are made simple** (`simple_rings`): repeats and out-and-back spikes removed, and a ring
  meeting itself at a vertex split there. They come from nodes exactly at a level, whose crossings
  land on the node; areas are unaffected, but GIS validity checks reject such rings.
- **Areas are planar, in the local frame**, before the projection to longitude and latitude.
- The unit tests pin the geometry against analytic fields (a cone's disc, an annulus with its
  hole, two peaks, the grid edge, culled nodes, both saddles), an island in a hole, a U-shaped
  hole, exact ties, and 300 random integer fields (simple rings, positive polygons, holes inside
  their outer rings); on the demo footprint the areas agree with a count of grid cells above each
  level to 1–2%.

## Time weighting, day metrics, C-SEL

**The period levels** (`DayMetric::lday` … `laeq` in `util/sel.hpp`) reuse the penalized-exposure
path: a sample outside the period takes a penalty of −∞, so it carries no energy and the 10 dB-down
window is found among the samples inside, and the energy is averaged over the period's length
(`averaging_seconds`). An operation wholly outside gives −∞, not NaN (`peak_range` of an all −∞
history takes every sample, whose sum is 0). The CTests pin them against each other on the
geodetic demo track: LAeq 0–24 equals Lden by day, Lday equals LAeq 7–19, Lnight equals LAeq 23–7,
and a night flight leaves Lday −∞; `period_*_over_*_hours` pin each period's normalization
against Lden (`golden_compare`'s offset), which the pairs above share and so cannot see; "Lden is
its period levels recombined" pins the arithmetic.


`--lamax_weighting`, `--nax_file`, `--ta_file`, `--lden_file`, `--cnel_file` and `--csel_file`
(`src/util/time_weighting.hpp`; the day schedules are in `src/util/sel.hpp`). What to keep:

- **The time weighting runs in reception time, on the SEL integral's own weights.**
  `reception_durations` is T times `reception_factor`, the very weights `sel_over_window` applies,
  so a maximum, a time above and an exposure from one history hear it on one clock. The samples
  are a piecewise-constant mean square over those durations, and the weightings are evaluated on
  it **exactly**: the exponential's maximum is at a sample's end (it relaxes monotonically within
  one), and the running 0.5 s Leq is linear between the instants where *either* edge of its window
  crosses a sample boundary (`for_each_window_piece` merges both sets). The detector starts settled
  on the first sample's level.
  - **Do not resample onto a grid or block it first**, as the validation harness does to mirror
    the measured 0.1 s frames. That blocking is what the harness's LAFmax differs by, not the
    filter.
  - Dropping the trailing-edge breakpoints passes every closed form but one, `the running Leq
    turns where the window's trailing edge crosses a sample`. Using one step for every sample
    instead of each duration fails `the weighting follows each sample's own duration`.
- **The unit tests hold the exponentials to their closed forms**
  (`tests/test_time_weighting.cpp`): the step response, with the time above half power, at three
  sampling steps; IEC 61672-1's tone-burst response 10 lg(1 − e^(−T_b/τ)), which its Table 4 lists
  as −1.0 dB (F, 200 ms), −18.0 (F, 2 ms), −27.0 (F, 0.25 ms) and −7.4 (S, 200 ms); a burst cut
  into unequal samples reading as the same burst; a source closing at Mach 0.5, whose 0.2 s burst
  the detector hears as 0.1 s; and the running Leq's ramps, including a crossing that falls where
  the window's trailing edge, not a sample boundary, turns.
- **The call site has its own guard**: `lafmax_golden`, `nax_golden` and `ta_golden` (demo case,
  broadband S-76D, `--gridstep 1640`), checked against an independent brute-force detector
  ([FIDELITY](../FIDELITY.md#regression)). Handing `reception_durations` no ranges (emission
  time) fails the first and last; NAx counting the instantaneous peak fails `nax_golden`.
  Nothing else in the suite caught either. Their provenance is in
  [testing.md](testing.md#golden-files); they move with anything that moves
  `regression_golden`, so regenerate them from the same build, and say so.
- **The τ/2 refusal is on the integration step (the mean), not the longest interval.** The 2017 AC
  tracks are 0.02 s with dropouts up to 0.14 s. Testing the longest interval refused every one of
  them, while the durations the detector actually sees are T times a reception factor and never
  contain the dropout. A coarser track is refused with the `-s` advice rather than given an
  exponential that settles within every sample.
- **`instantaneous` returns `*max_element(levels)`, not a round trip through the mean square**, so
  it reproduces the unweighted `--lamax_file` to the bit (checked against `6bbd6bb` on four 2017
  runs). `fast` is the default; the change is recorded in
  [HISTORY](../HISTORY.md#changes-that-affect-results).
- **The day metrics add a penalty per sample, 0.0 outside its periods**
  (`penalized_exposure_level`), so `-D` stayed byte-identical to the night-flag form it replaced
  (`day_night_level`, kept as a wrapper for the tests). The window is found after the penalty.
  `--full_history` integrates the whole history for DNL, Lden and CNEL too. The CTests
  `lden_at_22h30_is_not_dnl` and friends pin the period boundaries by which metrics coincide at
  22:30 and 23:30 local. The clock refusal names the option asked for (`dnl_needs_a_clock` matches
  "-D needs the time of day").
- **C-SEL reuses `--band_history`'s per-band recording** (`record_bands`, `record_path`). Every path
  that writes `bands_out` for the band history (straight, table, refracted, each extra arrival)
  feeds it, so a new propagation branch that forgets `bands_out` breaks both. Bands outside
  [first_band, last_band) are never written, which is safe only because the C gains are exactly
  zero there (same `--band_limits`). It is refused with `--db_cutoff` (a state culled on the
  A-weighted bound has stale bands and may carry C-weighted energy) and with `-a`. The tests
  that check the wiring are in [testing.md](testing.md#c-weighted-sel).
- **`c_weighting` is the A curve's second factor, constant and all** (2.242881e16 =
  12194.22⁴ 10^(0.0619/10), C1000 being −0.062 dB to IEC's three figures); its agreement with
  IEC 61672-1 is in [FIDELITY](../FIDELITY.md#regression). The test that A − C equals the inner
  pole pairs to 1e-12 dB keeps the two from drifting apart.

[max-level]: ../MANUAL.md#maximum-level-time-above-and-c-weighted-sel
[day-night]: ../MANUAL.md#several-trajectories-and-the-day-night-level
