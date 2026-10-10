# Output metrics: EPNL, band frames, time weighting, day metrics, C-SEL

- `src/pnl.{hpp,cpp}`, `src/band_frames.{hpp,cpp}` — `--epnl_file` and `--band_frames`: EPNL
  (14 CFR 36 A36.4) from one-third-octave levels in 0.5 s reception-time frames, and the frames
  themselves as an (observer, frame, band) netCDF file. The things to preserve:
  - **`pnl.cpp` is panam's `flight_acoustics.py` (eeb040b with F005, F014 and F026), term for
    term**, including panam's reading of the regulation where readings differ: 500 Hz is in Table
    A36-2's *middle* range; the duration's limits are the records closest to PNLTM - 10
    (A36.4.5.5), dips inside included; D takes A36.4.5.4's -13 for 0.5 s records and 10 lg(dt/10 s)
    for any other step; the A36.4.4.2 band-sharing check is applied, with the five-record mean
    truncated at the record's ends, its Delta_B added to EPNL as a term of its own (PNLTM, the
    limits and D stay the unadjusted history's); and a -inf band enters the tone correction at
    SPL(d) with no tone of its own, a NaN band propagates (and is refused in an EPNL history). panam's comments say the MATLAB
    `EPNLcalc.m` differs on the first three on purpose; do not "fix" the port toward it. A change
    here is a change to panam first.
  - **The reference values in `tests/test_pnl.cpp` are generated, not typed**:
    `PANAM=... python tools/pnl_reference.py` prints the block between the `clang-format off/on`
    lines. Regenerate rather than edit. Two spectra exist only for mutants that every other case
    let through: a 10 kHz tone *on a slope* (on a flat spectrum s(23) = 0, so dropping it from
    SPL'(24) = SPL(23) + s(23) is invisible) and a rise that levels off to a slope of exactly 0
    (the only place step 3(b)'s "zero or negative" differs from "negative"; a plateau after a
    sharp rise is encircled by 3(a) anyway). The inputs are 0.1 dB literals and the histories
    `base + envelope (+ add)`, so both sides parse and add the same doubles.
  - **`--epnl_masking` (default on) is not optional for a prediction** (panam's `masked=True`). A band
    below its noy threshold SPL(d) enters the tone correction's steps 1-7 at SPL(d) and carries no
    tone itself (`ToneMasking::below_noy_floor`; the library's default stays `none`, panam's). The
    literal procedure reads absorption's plunge in the top bands (10 kHz at a median -194 dB at
    PNLTM), and an inaudible band rising out of it, as tones of up to 20/3 dB: on the demo track over a
    100 m grid 615 of 961 observers took their PNLTM tone from an inaudible band and EPNL read a
    median 1.0 dB, up to 5.4, high. A measurement never sees this (its noise floor), and
    certification excludes tones in masked bands (FAA AC 36-4 App. 2). **Both halves are needed**:
    exclusion alone leaves an audible band before a plunge with a crashed background (20/3 dB),
    and reconstruction alone lets masked 8 and 10 kHz bands, rising 17 and 21 dB along SPL(d),
    read as a tone; the reference table has a spectrum for each, and dropping either half fails
    it. The references are panam's `masked=True` (F005, which adopted this definition). Masking can raise an EPNL a little (0.2 dB once on that grid: the
    band-sharing check and the duration pick other records), so `epnl_masking_lowers_tones`
    asserts only that some EPNL falls, never that none rises. On the 2017 measured spectra it moves
    EPNL 0.01 dB rms; on the predictions it is what made the 24-band EPNL agree with the 5 kHz-
    filled one (it read 0.04-0.19 dB above, all of it inaudible tones).
  - **A frame is an energy integral, not a mean of the states in it.** Each state's band energy
    is spread over [t_r - w/2, t_r + w/2], t_r = t + r/c and w = T x `reception_factor` -- the
    SEL integral's own weight -- and split between frames by overlap. That is what makes the
    frames sum to the whole-history exposure exactly, and `check_band_frames` holds the
    frames' A-weighted energy to a `--full_history` run's SEL at 1e-4 dB (float storage): drop
    the reception factor, or file states whole under their center's frame, and it fails. It is
    blind to *when* the energy arrives, so `check_band_frames --history` rebuilds every frame from
    a `--band_history` run of the same case (whose reception times are formed separately, where the
    history is written) and compares them frame by frame, the partial first and last frames
    excepted. Frames sit at multiples of 0.5 s on the trajectory's clock so a grid's observers share
    one axis.
  - **A missing PNL band is refused, never filled.** The databases' bands are matched to the
    nominal centers within 5 percent (`pnl_band_indices`); a database without one, or
    `--band_limits` leaving one out, stops the trajectory with the bands named. Filling (the
    harness's 5 kHz - 5/10/15 dB for measured frames) belongs to the comparison, not the tool.
  - **Without `--band_history`, `record_path` returns the band buffer and skips the geometry**
    (`lookup_angles`, an asin): it runs per state per observer. The accumulation uses
    exp(L ln10/10), as `weighted_level` does.
  - **`--band_frames` holds every observer's frames, as float, until the write** (about 450 MB on
    the demo grid, 0.9 GB resident). Each observer's record is one chunk, shuffled and deflated at
    level 1: +2.2 s and 28 percent smaller there. `nc_close`'s status is checked, as the sphere
    writers check it.
  - **`EpnlRow` keeps scalars, not `EpnlResult`'s per-record histories**, which on a grid were
    four vectors per observer kept only to read one element.
  - Against the 2017 data (four detail runs, `docs/MANUAL.md`): predicted minus measured EPNL,
    masked on both sides and both limited to 5 kHz, where the measured frames stop, -0.33 / -0.09
    dB median, 1.22 / 1.26 rms on the calm approaches (`station_plate` / `stratified_plate`), -0.22
    / -0.88 and 2.70 / 1.90 on the windy steady runs (unmasked: -0.25 / +0.02, 1.23 / 1.33, -0.22
    / -0.80, 2.85 / 1.87). PNLTM reads 0.75-1.4 dB low, about half of it a tone correction
    0.4-0.5 dB smaller than measured.

## Output metrics: time weighting, the day schedules, C-SEL

`--lamax_weighting`, `--ta_file`, `--lden_file`, `--cnel_file` and `--csel_file`
(`src/util/time_weighting.hpp`, the schedules in `src/util/sel.hpp`; the manual's "Maximum level,
time above and C-weighted SEL" has the definitions and the measured comparison). What to keep:

- **The time weighting runs in reception time, on the SEL integral's own weights.**
  `reception_durations` is T times `reception_factor`, the very weights `sel_over_window`
  applies, so a maximum, a time above and an exposure from one history hear it on one clock. The
  samples are a piecewise-constant mean square over those durations and the weightings are
  evaluated on it **exactly**: the exponential's maximum is at a sample's end (it relaxes
  monotonically within one), and the running 0.5 s Leq is linear between the instants where
  *either* edge of its window crosses a sample boundary (`for_each_window_piece` merges both
  sets). Do not resample onto a grid or block it first, as the harness does to mirror the
  measured 0.1 s frames: that is what the harness's LAFmax differs by (+0.03 to +0.06 dB median),
  not the filter. Dropping the trailing-edge breakpoints passes every closed form but one, `the
  running Leq turns where the window's trailing edge crosses a sample`; using one step for every
  sample instead of each duration fails `the weighting follows each sample's own duration`. The
  call site has its own guard: `lafmax_golden`, `nax_golden` and `ta_golden` (demo case, broadband
  S-76D, `--gridstep 1640`), checked independently with a brute-force reception-time detector to
  1e-4 dB. Handing `reception_durations` no ranges (emission time) fails the first and last; NAx
  counting the instantaneous peak fails `nax_golden` (64 dBA sits between them at three cells).
  Nothing else in the suite caught either. They move with anything that moves `regression_golden`;
  regenerate them from the same build, and say so.
- **The tau/2 refusal is on the integration step (the mean), not the longest interval.** The
  2017 AC tracks are 0.02 s with dropouts up to 0.14 s; testing the longest interval refused
  every one of them, while the durations the detector actually sees are T times a reception
  factor and never contain the dropout. A coarser track is refused with the -s advice rather
  than given an exponential that settles within every sample.
- **`instantaneous` returns `*max_element(levels)`, not a round trip through the mean square**,
  so it reproduces the pre-option `--lamax_file` to the bit (checked against `6bbd6bb` on four
  2017 runs). `fast` is the default since this change; the manual's "Changes that affect results"
  says so.
- **The day metrics add a penalty per sample, 0.0 outside its periods**
  (`penalized_exposure_level`), so `-D` stayed byte-identical to the night-flag form it replaced
  (`day_night_level`, kept as a wrapper for the tests). The window is found after the penalty;
  `--full_history` now integrates the whole history for DNL, Lden and CNEL too, where DNL used
  to ignore it. The CTests `lden_at_22h30_is_not_dnl` and friends pin the period boundaries by
  which metrics coincide at 22:30 and 23:30 local; the clock refusal names the option asked for
  (`dnl_needs_a_clock` matches "-D needs the time of day").
- **C-SEL reuses `--band_history`'s per-band recording** (`record_bands`, `record_path`): every
  path that writes `bands_out` for the band history -- straight, table, refracted, each extra
  arrival -- feeds it, so a new propagation branch that forgets `bands_out` breaks both. Bands
  outside [first_band, last_band) are never written, which is safe only because the C gains are
  exactly zero there (same `--band_limits`). It is refused with `--db_cutoff` (a state culled on
  the A-weighted bound has stale bands and may carry C-weighted energy) and with `-a`.
  `csel_at_1khz_is_sel` checks the wiring (both curves are 0 dB at 1 kHz, to 1.5e-6), and
  `csel_is_not_sel` that the C gains are what is applied.
- **`c_weighting` is the A curve's second factor, constant and all** (2.242881e16 =
  12194.22^4 10^(0.0619/10), C1000 being -0.062 dB to IEC's three figures). The test that A - C
  equals the inner pole pairs to 1e-12 dB keeps the two from drifting apart.
