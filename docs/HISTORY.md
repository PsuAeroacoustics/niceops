# NICE-OPS history

How NICE-OPS got to its current form, and why old results differ from new ones. For what the
tool does now, see the [README](../README.md) and the [manual](MANUAL.md); for how accurate it
is, see [FIDELITY](FIDELITY.md).

## Origins (2020–2021)

Eric Greenwood and Rob Rau developed NICE-OPS in D, starting in January 2020. It used a
CGAL-based barycentric interpolator, the SPLINTER B-spline library and matplotlib-d. The S-76D
case, with `S-76D_M3.nod`, became the default in March 2020. The AW139, B205 and B407 databases
and tracks from the joint NASA/FAA flight tests were added in January 2021, together with a
second S-76D database. These are the `_M` files in `databases/`: `AW139_M1`, `B205_M1`,
`B407_M2`, `S-76D_M3` and `S-76D_M4`, plus `B407_M3` from the same family, which was added to
the repository later. All are broadband, with discrete load factors.

The June 2020 reference figures came from that code; see
[Reference figures (June 2020)](#reference-figures-june-2020).

## The C++17 port (September 2026)

In September 2026 NICE-OPS was ported to C++17 (`666046b`) and checked against the D reference.
On the demo case the port matched D to 1e-9 relative, with an actual maximum difference of about
1e-14, a few ULPs of LDC-vs-clang floating point. Each D dependency was replaced with a small
local equivalent. [`internals/porting.md`](internals/porting.md) lists the replacements, the four
mistranslations found during validation, and the D quirks kept on purpose. The vendored SPLINTER
library, which needed special compiler flags on a Command Line Tools-only Mac, went with the
port.

The D sources were removed once the port was verified complete. They remain in history at
`c74bab3` on `origin/no-stk`. Check that commit out when the D oracle is needed again, for
example through the `regenerate-goldens` skill.

Two D settings became options during the port, both defaulting to the port's behavior so the
goldens did not move: `watts_trim` became `--trim`, and the observer grid's extent became
`--grid`. The contour levels changed default from D's 45–90 dB to 65–100 dB: the older pair put
its top level below the demo case's 93.6 dB peak, leaving the loudest patch unfilled, and the
reference figures were drawn at 65–100. That affects figures only; the goldens are CSV.

## Departures from D

Port fidelity is the rule: a D behavior is changed only when a measured case shows it is wrong,
not when a reading suggests it looks wrong. Three changes have met that bar.
`golden/state_197421.csv` is trajectory state, which none of them touches, so it is still
D-derived; `golden/197421_LAE.csv` now comes from the C++ build. D can no longer produce it, and
regenerating it from D would only restore the departures.
[`internals/testing.md`](internals/testing.md#golden-files) gives the provenance of each golden.

### The 4x EAA multiplier

**Results from before September 2026 are not comparable to results after.** `eval` used to
multiply the database's `EAA` by a constant 4.0, computing `4.0 * EAA * (range - r)/1000` under a
`//TODO Add ground effects` inherited verbatim from the D original (`c74bab3`). That is not a
ground effect, and it overstated atmospheric absorption fourfold: 168 dB at 6 km, where the
physical value is 14. panam *writes* that `EAA` as `SPLA - SPLAa` over a 1000 m reference and
consumes it as `EAA / distance * (slant_range - radius)`, with no such factor. NICE-OPS does the
same.

Before it was removed, it was measured against `--ground_effect`, the model it was standing in
for. Soft ground costs a roughly range-independent 2–4 dB relative to a rigid surface, whereas
the factor removed 7 dB at 250 m, rising to 38 dB at 2 km. A real ground effect is bounded by the
surface impedance; a multiplier on a distance-proportional absorption is not.

Removing it raised the demo footprint by **9.93 dB mean, up to 22.7 dB**, most of it far from the
track, where the old factor removed the most. `golden/197421_LAE.csv` was regenerated from the
C++ build at the commit that removed the multiplier. The port's agreement with D was verified
before that change, and that verification stands on its own.

### The ambient mix

D added the `-a` ambient level as pressure amplitude, which adds two uncorrelated levels
coherently: +6.02 dB over two equal levels, and up to 3 dB high in general. It had never been
tested. `-a` now adds the ambient in energy, so two equal levels give +3 dB.

### Spreading capped at the sphere surface

D extrapolated the spherical-spreading law inward (1/r) inside the database sphere, which put
observers close to the track up to 4.4 dB above the sphere-surface level in the instantaneous
level, and diverged at r = 0. Spreading is now capped at the sphere's surface (`eval` and
`eval_spectrum` take max(range, R)), and absorption starts at the surface. `golden/197421_LAE.csv`
was regenerated a second time, from the build that carries the change: 27 of its 14,400 cells
moved, by at most 2.0 dB in SEL. They are exactly the observers within about 45 m of the track,
whose closest approach to a source state is inside the 30.48 m sphere.

## Changes that affect results

Besides the departures from D above, these corrections have moved NICE-OPS's results. Compare
old and new results only with them in mind.

### Metrics

- **`-D` is normalized to a 24-hour day.** It used to write the raw sum, which read
  10 lg 86 400 = 49.37 dB high
  ([Several trajectories, and the day-night level](MANUAL.md#several-trajectories-and-the-day-night-level)).
- **`--start_utc` is required** with `-D`, `--lden_file` or `--cnel_file` for a WGS84 track
  without a `utcsec` column. It used to default to midnight, which is night, so the whole
  footprint took the 10 dB penalty.
- **`--full_history` governs `-D` as well as SEL.** `-D` used to leave each operation on its own
  10 dB-down window whatever `--full_history` said. The default is unchanged; how far the two
  windows differ is in [FIDELITY](FIDELITY.md#day-night-window-choice).
- **`--lamax_file` writes LAFmax**, the IEC 61672 fast-weighted maximum in reception time, by
  default. Before `--lamax_weighting` existed it wrote the instantaneous maximum of the 0.02 s
  history, which is not a meter quantity and reads high against a meter's LAFmax (by how much:
  [FIDELITY](FIDELITY.md#lamax-ta-and-epnl-against-the-2017-measurements)).
  `--lamax_weighting instantaneous` reproduces the old file to the bit, which was checked against
  the build before the change on four 2017 runs.
  `--nax_file` counts against the same weighted maximum
  ([Maximum level, time above and C-weighted SEL](MANUAL.md#maximum-level-time-above-and-c-weighted-sel)).
- **Coarse tracks are refused for time-weighted outputs.** A trajectory sampled more coarsely
  than half the weighting's time constant is refused, so WGS84 waypoints at the default
  `-s 0.5`, which were accepted before `--lamax_weighting` existed, are refused under `fast`.
  They need `-s 0.0625` or finer for `--lamax_file` or `--nax_file`, or
  `--lamax_weighting instantaneous` for the old number.
- **`-m` has a fifth column**, the area at or above 65 dB (October 2026). A script that read
  the line as four fields must take the fifth or ignore it; the first four are unchanged.
- **`--epnl_file` follows panam's corrected EPNL** (October 2026). The duration's limits are the
  records closest to PNLTM − 10; D, the duration constant, takes the regulation's −13 for the
  0.5 s records; and the band-sharing adjustment is added to EPNL rather than to a PNLTM that
  canceled it. Every EPNL rises, by 0.0103 dB from the constant alone. Over the demo grid at
  100 m with the B407 database (961 observers), the mean rise was 0.031 dB (0.040 with
  `--no-epnl_masking`) and the largest 0.41 dB (0.62). Without masking, a −∞ band no longer
  reads as a tone in its neighbor
  ([EPNL and reception-time frames](MANUAL.md#epnl-and-reception-time-frames)).

### Receivers and ground

- **The point receiver's microphone height sets the evaluation point.** Before, the point
  receiver's direct path ended on the ground, even though its ground model put the microphone
  1.2 m up. Correcting that moved point-receiver SEL on the demo grid by up to about ±0.8 dB
  (median 0). Almost all of it comes from the low states at the end of the approach, where
  1.2 m is a sizeable fraction of the source's height. Plate and free-field results, and the
  goldens, did not move ([Receivers](MANUAL.md#receivers)).
- **The ground model's boundary-loss term has no cutoff.** panam cut that term to zero past
  |w| = 500, and NICE-OPS followed, until panam `e978703` evaluated it by `wofz` at any |w|.
  NICE-OPS does the same, pinned to panam `b758e3a`. The reflection coefficient Q moves by up to
  2e-3 just past the old cut. The demo footprint with the plate receiver moved by at most
  0.028 dB, and by 0.08 dB with a 4 ft point receiver, where the two paths nearly cancel; the
  means moved by about 1e-4 dB. The near-grazing three-term series that stood in past the cut
  was removed, since the exact term is what it approximated.

- **`--coverage_file` shares out the energy the SEL integrated** (October 2026): its 10 dB-down
  window, each sample weighted by the step and its share of reception time. It summed every
  sample of the history alike, which counted states the SEL leaves out and, across operations,
  weighted each by its sample count. On the B407 demo grid the median observer's unmeasured
  share moved from 1.2% to 0.94%, the 95th percentile from 27% to 31%.

### Refraction and the parabolic equation

- **Stratified propagation keeps every arrival it solves.** It used to keep only the 24 arrivals
  with the smallest spreading ranges, which discarded energy without accounting for source
  directivity or ground losses. Deep ducts can therefore take longer to evaluate and predict
  more energy
  ([Refraction through the profile's layers](MANUAL.md#refraction-through-the-profiles-layers)).
- **`--refraction stratified` resolves a calm floor** only as finely as the run's highest band
  can. A profile with calm, isothermal air at the ground (a log law written below its roughness
  length) used to put the shadow boundary kilometers out and flush microphones 22–28 dB low.
  Profiles without such a floor moved by under 0.001 dB.
- **A PE source near the ground is placed as the adjoint of the receivers' read.** The
  fourth-order ground row differs from the rows above it, so a source within three grid points
  of the ground was misplaced, and a source on the ground (a flush microphone in a reciprocal
  march) had no self-starter at all. In propagation tables built before this, the 5 ft row's
  50–100 Hz bands were up to 0.39 dB off past the first range (0.12 dB mean); nothing else moved
  by more than 0.0002 dB.
- **Plate microphones under `--propagation_table` take the plate's offset at each pair's
  elevation** (October 2026). It was taken at 2° for every pair, which against the plate run the
  table should reproduce put footprints read through a flush receiver row up to 0.34 dB low
  ([FIDELITY](FIDELITY.md#propagation-tables)). Plate runs with a table that had no row at 0
  were refused, naming 0 ft as an uncovered height; they now run, with a warning.
- **The PE self-starter's ground source applies the second-order correction its derivation
  calls for** (October 2026). It had been computed and then discarded. It only affects sources
  within three grid spacings of the ground, mostly a PE table's flush-microphone rows in a
  reciprocal march. Against the exact field at 17–45°, those rows moved from 0.03–0.08 dB off to
  0.01–0.05 dB; at grazing nothing moved by more than 0.001 dB
  ([Propagation loss tables](MANUAL.md#propagation-loss-tables)).

### Trajectories and the source

- **The hover nose is turned by the frame** whenever the run knows one. With a frame known, the
  heading below 8 kt in `velocity` and `ground_velocity` modes used to be turned by a fit to the
  track, which absorbed the track's mean crab. On a crabbed approach ending in hover, that fit
  turned the hover nose by the approach's mean crab (how far on the day-286 approaches:
  [FIDELITY](FIDELITY.md#the-spheres-nose)). The fit is kept only for a run with no frame. This
  changes only tracks that both fly above 20 kt for 50 states and slow below 8 kt. No shipped
  track goes under 12 kt (the slowest, AW139 `192217AC.csv`, reaches 12.5 kt), no 2017 B407
  flight track other than the day-283 hovers goes under 13 kt, and the hovers have no fast
  states, so no golden or validation result moved ([The sphere's nose](MANUAL.md#the-spheres-nose)).
- **`--export_aam` writes SPEED with the exact knot** (October 2026). It converted m/s to knots
  with 1.94384, which wrote every sphere's speed 2.3e-6 (relative) low; no footprint moved.
- **A static operation with `--origin` and no `--frame_bearing`** reads its headings as compass
  bearings in the east/north frame. It used to orient the sphere by its zero velocity.
- **The `azimuth_reference` attribute.** panam's builder writes it on every database: `track`
  by default, `heading` for spheres built with its heading-frame option
  (`build_sphere(azimuth_reference='heading')`, `--azimuth-reference`), and `track` for a NORAH2
  import; it refuses to mix the two in one database. Every shipped database predates it and
  records nothing, which reads as `track`, so nothing shipped moved.
- **`--path_smoothing 0` carries the spline's curvature.** `SmoothingSpline` used to
  short-circuit at λ = 0 and return zero curvature, so the load factor came out as exactly 1 g
  on every state. With the interpolating spline carrying its own curvature, exact interpolation
  amplifies coordinate rounding into the load factor, as expected. The trajectory-smoothing
  study had reported that row as 0.196 dB, a figure that measured the bug
  ([FIDELITY](FIDELITY.md#trajectory-smoothing)).

### Inputs

- **Four length options were meters before `--units` existed**: `--grid`, `--gridstep`,
  `--receiver_height` and `--path_smoothing`, while `--trim` and `--ground_elevation` were feet
  and every file was feet. Scripts written against that need those four converted, or
  `--units m` added ([Units](MANUAL.md#units)).
- **`arrays/coyle_approach.csv` works as a `-b` input again.** Its fields are padded *before*
  the comma, and the CSV reader counted trailing whitespace as junk, so every row was dropped.
  The run then evaluated zero observers and exited 0. An observer file that yields nothing is
  now an error. All 84 shipped AC tracks pad on the left, which was always skipped, so the
  goldens never saw this.
- **Bad inputs are refused rather than run.** A bad input used to give a plausible but wrong
  footprint with exit status 0. Readers refuse non-increasing time, non-finite values, units
  that contradict `--units`, and text after the data. A spectral run whose relative humidity is
  given by no source is refused rather than run at an assumed 20 percent. A vehicle file whose
  tip speed, weight, rotor radius or C_T differs from the database's by more than 2 percent is
  refused unless `--accept_vehicle_mismatch` is given
  ([Errors and refusals](MANUAL.md#errors-and-refusals)).

## The 2017 Noise Abatement databases

The AS350, B407, Be206L3, EC130, RO-44 and RO-66 databases were built by panam from the 2017
Noise Abatement flight test. The current recipe is in the manual's
[Shipped databases](MANUAL.md#shipped-databases), and their scores are in
[FIDELITY](FIDELITY.md#shipped-databases-approach-scores-per-aircraft). The AW139, B205 and S-76D
`_fixed_load_factor` files were converted from those aircraft's discrete-load-factor `_M`
databases. They agree with their sources to 1e-6 dB at load factors from 0.5 to 3.0.

**The September 2026 files.** The first AS350, Be206L3, EC130, RO-44 and RO-66 builds had several
faults, each found and fixed in a later rebuild:

- **Unwritten last group.** A database written before panam `c53588f` (2026-09-26) may have its
  last group unwritten. Every entry is then the netCDF fill value (9.97e36), which the loader
  reads as a finite 1e37 dB level. The September AS350, EC130, RO-44 and RO-66 files had this,
  and spectral runs that reached the fastest +35-degree condition gave an infinite SEL.
- **Near-silent cells.** All four of those files carried finite levels far below the aircraft
  (as low as −70 dBA beside an 84 dBA neighbor), mostly in the abeam column near the horizon,
  which the loader blends as real. The September AS350 had 441 near-silent cells and 1,305 empty
  ones (the rebuild has 0 and 10); RO-66 had 62 and 3,962 empty ones below the horizon (the
  rebuild has 0 and 108).
- **EC130 depropagated at a standard day** (20 °C, 20% RH, 101.3 kPa), where the stations read
  −1 to 19 °C, 37–94% and 88–89.5 kPa. The dataset spells the folder `EC130B4_Ground Stations`,
  and panam looked for `EC130B4_Ground_Stations`. The September EC130 spheres also lacked
  day 296, which supplies 62 of the rebuilt file's spheres.
- **Missing days.** RO-44's September spheres lacked day 227, which holds 32 of its steady runs.

**The October 2026 rebuilds.** AS350, Be206L3, EC130, RO-44 and RO-66 were rebuilt on
2026-10-04 (the AS350 on 2026-10-05) on the B407's flight recipe and its de-Dopplerized hover
source, without a hover correction. Against each aircraft's approaches (NICE-OPS minus the
measured 10 dB-down SEL, station atmosphere, median / rms in dB), the rebuild moved the scores
as follows:

| Database | Spheres, September → rebuilt | Approaches | September | Rebuilt |
| --- | --- | --- | --- | --- |
| AS350 | — → 96 | day 291, 1328 mic-runs | +0.76 / 1.76 | +0.70 / 1.49 |
| Be206L3 | 110 → 110 | day 279, 906 mic-runs | +0.28 / 1.48 | +0.24 / 1.50 |
| EC130 | 38 → 112 | day 298, 1025 mic-runs | +0.19 / 3.17 | +0.70 / 2.70 |
| RO-44 | 50 → 83 | day 229, 585 mic-runs | +0.31 / 1.55 | +0.04 / 1.16 |
| RO-66 | 87 → 88 | day 236, 876 mic-runs | +0.28 / 1.06 | +0.87 / 1.56 |

The RO-44 and RO-66 rows are after the frame-bearing rebuild below.

- **AS350's gain is at low elevation.** At 5–10 degrees its over-prediction fell from
  +2.73 / 3.46 to +1.59 / 2.29. Above 20 degrees nothing moved, and the bias overhead remains
  in both files. 96 of its 107 steady runs built: six had no steady segment of 8 s, three had
  no valid samples, one had a window longer than its recording, and the track check refused
  291260, which passes 1.2 ft from a board.
- **RO-66 reads worse, because the September file was right for the wrong reason.** Its
  near-silent cells had hidden the low-elevation over-prediction every aircraft shows. The
  rebuild raised RO-66's predictions by 0.2 dB under the track and 1.3 dB at 1,200–2,000 ft off
  it, where the aircraft is heard abeam, while the cells both files measured moved by +0.05 dB
  at the median.

**The frame-bearing rebuild (2026-10-08).** RO-44 and RO-66 were rebuilt again for a
frame-bearing error in the validation harness's sphere build. Its stratified rays took the track
frame's +x bearing as 270° for every aircraft, but the Eglin layout these two flew runs along
140°, so the LIDAR winds were rotated by 130° in 55 of RO-44's 86 runs and 82 of RO-66's 98.
Nothing else changed: the recipe, panam `eeb040b` and NICE-OPS `e13577e`'s ray tool are the
2026-10-04 build's, and every straight-ray sphere came out bit-identical. Before the fix RO-44
scored +0.12 / 1.27 and RO-66 +0.86 / 1.57. RO-44's gain is below 20 degrees (0–5 degrees, rms
2.36 → 1.94; 5–10 degrees, 1.41 → 1.23); RO-66 is unchanged. NICE-OPS predictions were never
affected, since they fit each track's own bearing.

## Reference figures (June 2020)

The June 2020 reference figures were produced by the D code with three settings that differ from
the defaults: a trimmed trajectory, a shifted observer grid, and contour levels of 65–100 dB.
They were kept in a local, untracked `figures/` directory and are not in the repository. All
three settings are reachable from the C++ build; the contour levels are the default:

```sh
./build/niceops -t tracks/S-76D/197401AC.csv -v databases/S-76D.json -d databases/S-76D_M3.nod \
        --units m --trim -1828.8 --grid -2000 1100 -1500 1600 -o out.csv -f out.eps
```

`--units m` is needed because the grid bounds were meters when the comparison was made; `--trim`
was feet then (−6000 ft). Matched that way, the two runs land on the same 124 × 124 grid at the
default 25 m step and agree cell for cell: the largest difference anywhere is **1.05 dB**, 99% of
cells are within 0.5 dB, and the peak levels differ by 0.07 dB.
