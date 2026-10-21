# NICEOPS fidelity

How far each part of NICEOPS can be trusted, and the evidence for it. The equations are in
[THEORY.md](THEORY.md), how to run each option is in the [manual](MANUAL.md), and changes that
moved results are in [HISTORY.md](HISTORY.md).

## What is measured here

This page keeps three kinds of evidence apart:

- **Regression agreement.** A run reproduces a stored reference. Passing a golden comparison
  establishes stable behavior. It does not establish physical accuracy.
- **Numerical accuracy against a reference solver.** A model is set against an exact solution,
  or against a full-wave parabolic equation (PE). A PE reference is not a measurement.
- **Accuracy against measurements.** Predictions are scored against the 2017 Noise Abatement
  flight test ([Data sets](#data-sets)), and one urban case against a published measurement
  ([NLR Drone DigiCity](#a-loudspeaker-among-buildings-nlr-drone-digicity)).

The historical grazing-incidence reports in `reports/` record successive hypotheses and
rescorings; read each in the light of its latest addendum.

Each measured result appears once, here, with the aircraft, database, data set and metric it
was measured on. Unless a row says otherwise, levels are in dB, errors are predicted minus
measured (or model minus reference), and LAmax is fast-weighted.

## Regression

- **Trajectory state.** `golden/state_197421.csv` comes from the D original. The C++ port
  matched it to 1e-9 relative (about 1e-14 in fact, a few ULPs of compiler round-off). No
  propagation change touches it.
- **The demo footprint.** `golden/197421_LAE.csv` comes from the C++ build, because it departs
  from D on purpose (see [HISTORY](HISTORY.md#departures-from-d)). The port agreed with D on the
  demo case before those departures.
- **The other goldens** (the LAmax, NA and TA metrics, the stratified runs and the terrain
  runs) come from the C++ build. The metric goldens were checked against a brute-force
  detector run on a 5e-5 s reception-time grid built from the run's state and level histories.
  It agrees with them to 1e-4 dB; the same detector in emission time is up to 0.018 dB off. On
  that case NA64 counts three cells fewer than the instantaneous peak would.
- **`self_test`** checks every sphere of the demo database against its own stored dBA and EAA
  at 1e-9.
- **The windowed peak search is exact.** For each model, a `*_windowed_is_exact` test checks
  that skipping states changes nothing against marching every state.
- **Physics ported from panam is pinned to it:** ISO 9613-1 absorption, the ground effect,
  the A-weighting and EPNL. EPNL agrees with panam to 1e-6 dB (the cases are listed in
  [internals/testing.md](internals/testing.md#epnl-and-band-frames)). The plate receiver's exact
  form matches panam's
  board level to 1e-6 dB. The C weighting matches IEC 61672-1 Table 3 within its 0.1 dB
  rounding from 10 Hz to 20 kHz.

The goldens, how each is made and the tests that guard them are listed in
[internals/testing.md](internals/testing.md#golden-files).

## Numerical verification

### Parabolic equation against exact solutions

The PE (split-step Padé, [THEORY §9.1](THEORY.md#91-the-wide-angle-parabolic-equation)) is
checked against exact solutions:

| Case | Agreement |
| --- | --- |
| Two rays over rigid ground | within 0.3 dB |
| Weyl–van der Pol field over grass | within 0.5 dB |
| Residue series in an upward-refraction shadow, rigid ground | within 0.02 dB at 250, 500 and 1000 Hz, down to 47 dB below free field |
| Benchmark case 1: direct and image rays over the benchmark's ground | within 0.2 dB at 10 Hz from 50 m (0.29 dB at 25 m, inside a wavelength of the source), and 0.08 dB at 100 Hz and 1 kHz from 20 m, to 1 km |
| Benchmark case 2: normal-mode series, downward refraction | within 0.06 dB at 100 Hz to 1 km and at 1 kHz to 500 m, 0.2 dB at 10 Hz |
| Benchmark case 3: residue series, upward refraction, impedance ground | within 0.03 dB at 100 Hz and 1 kHz, to 60 dB below free field |

The shadow rows matter most for judging the ray models. Their profiles have k² linear in
height, and the field in the shadow is then exactly a sum of creeping waves. Doubling the grid
resolution moves no shadow onset, and no level past 100 m, by more than 0.7 dB.

**Range step.** Between 2 and 32 wavelengths per step the result moves by less than 0.001 dB,
60 dB into a shadow.

**Starter and grid order.** Against the exact two-ray field over soft ground at 1 kHz, from a
100 m source (rms error):

| Elevation | Gaussian starter, second order, 20 points per wavelength | self-starter, fourth order, 10 points | self-starter, fourth order, 20 points |
| --- | --- | --- | --- |
| 0–20° | 0.01–0.02 dB | 0.005 | 0.005 |
| 20–30° | 0.14 | 0.006 | 0.003 |
| 30–35° | 0.59 | 0.01 | 0.002 |
| 35–40° | 0.61 | 0.02 | 0.004 |
| 40–50° | 1.75 | 0.01–0.07 | 0.001–0.007 |

Grid convergence on `atmospheres/example_windy.csv` out to 1.5 km, 250 Hz – 2 kHz, against an
80-point run (rms, worst in parentheses):

| | 20 points | 10 points | 6 points |
| --- | --- | --- | --- |
| second order throughout | 0.02–0.04 dB (0.14) | 0.06–0.19 (0.64) | |
| derivative fourth order, ground and interpolation second | 0.01–0.02 (0.06) | 0.05–0.09 (0.29) | |
| fourth order throughout | under 0.0005 (0.001) | 0.001–0.002 (0.011) | 0.004–0.01 (0.08) |

The Gaussian starter, normalized with the wavenumber at the ground, puts a source aloft in a
refracting or windy column at the wrong strength: 0.11–0.15 dB high against the exact residue
series for a 100 m source where the sound speed is 3% lower than at the ground, and 0.10 dB on
`example_windy.csv`. The self-starter is within 0.04 dB.

**The Crank–Nicolson check.** `--scheme cn` (Crank–Nicolson on the Padé(1,1) equation) is kept
as an independent march. Through `example_windy.csv` from a 100 m source to 1.5 km,
125 Hz – 2 kHz, it differs from the split-step by 0.03–0.12 dB rms with the same starter and
grid order. Against the exact field in still air it is Crank–Nicolson that is off: 0.7 dB rms at
25–30° and 2.1 dB at 30–35°, at 1 kHz. Its own vertical grid was checked against a 30 × 30
points-per-wavelength Crank–Nicolson reference: 10 points per wavelength in height is within
about 0.1 dB rms for a receiver near the ground, but 0.47 dB rms (2.2 dB at worst) for one
1.2 m up at 1 kHz; 20 points brings that to 0.09 dB. Rerunning the stored shadow-calibration
cases at 20 points moved the PE by 0.16–0.39 dB rms (1.1 dB at worst) and left the calibration's
edge error against it unchanged.

**Limit.** In still air the PE's deep shadow, 40–80 dB down, is not physical: with no
turbulence nothing scatters energy into it ([Known limitations](#known-limitations)).

### The Attenborough et al. benchmark

The four cases of Attenborough et al. (JASA 97, 173 (1995)), run by `niceops_benchmark` to
200 m and to 10 km, set beside the paper's FFP and PE curves:

- All four cases agree with the paper in every panel.
- Where there is an exact solution, the PE is within 0.3 dB of it past 20 m on every 200 m
  panel.
- On the 10 km panels it is within 0.3 dB for case 1. For case 3 it is within 0.1 dB, or
  0.5 dB at 10 Hz, where the linear profile and the exact solution's profile part company
  within a few Airy lengths of the ground.
- Case 2 at 10 km depends on where the profile is cut off, as the paper says. The PE with the
  linear profile stopped at 1000 m matches the paper's curves in form, not null for null.
- Case 4 (a duct) has no exact solution. The PE reproduces the paper's curve feature for
  feature, down to the 100 Hz humps at 2.5–3.2 km and 4.3–4.7 km.

The benchmark's printed upward-refraction residue series (its eq. 11) has the ground condition
wrong, and case 2's exact solution carries steep returning rays no PE can. Both are described in
[THEORY §9.2](THEORY.md#92-analytic-reference-solutions).

### Refraction models against the PE

`niceops_pe_compare` runs one path through straight rays, `--refraction linear`,
`--refraction stratified` and the PE. Absorption and turbulence are left out of all four.

**Shadow onset.** Through `example_windy.csv` with a 1.2 m receiver, upwind, the range at which
each model falls 10 dB below free field and stays there:

| source | band | straight | linear | stratified | PE |
| --- | --- | --- | --- | --- | --- |
| 100 m | 125 Hz | never | 1.9 km | 1.2 km | 1.3 km |
| | 250 Hz | never | 1.4 km | 0.53 km | 1.07 km |
| | 500 Hz | never | 1.3 km | 0.53 km | 0.91 km |
| | 1 kHz | never | 1.4 km | 0.53 km | 0.84 km |
| 30 m | 125 Hz | never | 0.91 km | 0.72 km | 0.65 km |
| | 250 Hz | never | 0.46 km | 0.38 km | 0.50 km |
| | 500 Hz | 1.2 km | 0.42 km | 0.19 km | 0.38 km |
| | 1 kHz | never | 0.50 km | 0.23 km | 0.34 km |

**Downwind**, the rms difference from the PE beyond 400 m, 100 m source:

| band | straight | linear | stratified |
| --- | --- | --- | --- |
| 125 Hz | 2.8 dB | 4.0 dB | 0.9 dB |
| 250 Hz | 7.1 dB | 7.6 dB | 0.4 dB |
| 500 Hz | 3.1 dB | 0.9 dB | 0.9 dB |
| 1 kHz | 8.1 dB | 7.2 dB | 2.1 dB |

Downwind and in the lit region, `stratified` is right and the others are not. Upwind, the PE's
shadow begins between the two ray models' predictions, moving toward `stratified` as the
frequency rises: `stratified` starts the shadow too early and `linear` too late, by similar
amounts. The per-band boundary law of the shadow transition is the treatment
([Shadow-transition calibration](#shadow-transition-calibration)).

**The layered boundary.** A log-law wind is zero at the ground, so upwind the layered
model's limiting ray is fixed by the sound speeds at the source and the ground alone. For a
100 m source in `example_windy.csv`, the layered upwind shadow begins 535 m out for a
microphone 1.2 m up, where the linear fit puts it past 1.6 km (516 m and 1.40 km for a
microphone on the ground).

**Surface ducts.** Without the bounced rays of a duct, `stratified` was 8–24 dB under the PE
past 4,000 ft on day 286's profiles for 100–330 ft sources, and up to 54 dB under it at single
ranges. With them, the energy mean over 50 Hz – 2 kHz against the PE, per range bin, for a flush
receiver and a 100 ft source (the value without bounced rays in parentheses):

| day 286 profile | 1,000–4,000 ft | 4,000–8,000 ft | 8,000–13,000 ft |
| --- | --- | --- | --- |
| run 286264 | +0.4 | −1.1 (−12.2) | −0.5 (−24.3) |
| run 286270 | 0.0 | +0.8 (−7.2) | −2.5 (−8.1) |
| run 286276 | −0.7 | −1.5 (−10.7) | −1.1 (−9.2) |
| run 286282 | +4.4 | +0.2 | −3.1 (−17.1) |

Downwind on `example_windy.csv` and a log law, at 100 ft, the far third of a 10,000 ft path
comes within 1 dB of the PE (4.0–6.5 dB under it without bounced rays). Upwind, the linear lapse
and the day-284 profiles do not change.

**Arrivals missed at caustics.** Against a 2048-sample reference (12 bearings × 8 source heights
from 5 to 300 m × 2 receiver heights × 120 ranges from 50 m to 10 km), the share of pairs more
than 1 dB out:

| launch samples | `example_inversion.csv` | `example_windy.csv` |
| --- | --- | --- |
| fixed 32 | 11.3% | 1.05% |
| fixed 128 | 2.6% | 0.27% |
| fixed 512 | 0.47% | 0.07% |
| 32 + 8 per hop, capped at 512 (used) | 2.9% | 0.01% |

Scaling the sample count with the hops costs 4–8% of a stratified run.

**Geometry checks.** On a linear profile `stratified` reproduces the linear model's closed
forms, the shadow boundary included. On layered profiles its rays match numerically traced ones
to a micrometer. Uniform air gives the straight-ray result to 1e-9 for every receiver.

**A calm floor.** A profile with calm, isothermal air in its lowest centimeters leaves a flat
layer at the fastest height, along which plain ray theory pushes the boundary out without limit.
For a 30 m source upwind on a 7 kt log law at 10 kHz, the boundary:

| calm floor | none | 1 mm | 0.1 ft | 1 ft | 1 m |
| --- | --- | --- | --- | --- | --- |
| plain ray theory | 186 m | 893 m | 21.8 km | 216 km | 707 km |
| NICEOPS (travel capped at 1.5 k h²) | 186 m | 186 m | 187 m | 233 m | 544 m |

At 100 Hz the same source's boundary is 215 m (198 m under a 1 mm floor). The cap's margin over
ray theory's (4/3) k h² is what a floor can still add: with a cap of 2 k h², a 0.1 ft floor
moved the demo grid's flush microphones by up to 0.26 dB, and by 0.67 dB with the top band at
200 Hz; at 1.5 k h² it is 0.12 and 0.19 dB. The cap's phase region is integrated past the path's
ends. Cut off at the ends instead, a 1.5 m source over a 1.2 m microphone under a plain lapse had
its boundary brought in from 789 m to 34 m at 10 kHz.

Checks of the cap:

- On the demo grid, `example_windy.csv` with a 0.1 ft calm floor against the same profile
  without it: flush microphones differ by 0.12 dB at most (up to 39 dB with plain ray theory),
  and 1.2 m microphones by 0.05 dB.
- Over a sweep of log laws without a floor (headwinds of −1 to 3.6 m/s against the lapse,
  sources from 1.5 to 100 m, top bands from 200 Hz to 10 kHz) the rays are those of plain ray
  theory, except a few interior speed maxima under a light tailwind, which move by under
  1e-7 dB.
- Five cells of the `stratified_windy_point` golden move by up to 4e-6 dB, where the lowest
  layer of `example_windy.csv` is within 1e-4 m/s of flat near crosswind.

The measured-data checks are in [Refraction against flight data](#refraction-against-flight-data).

**Linear refraction on the demo track.** `example_windy.csv` is illustrative, not a measurement
(a neutral log-law westerly, 5.2 m/s at 10 m, over the ISA lapse). On the S-76D demo track with
the B407 database and `--frame_bearing 210`, `--refraction linear` against `--refraction none`
on the same profile:

| | mean | upwind half | downwind half | worst cell |
| --- | --- | --- | --- | --- |
| no ground effect | −0.77 dB | −1.51 dB | −0.04 dB | −7.85 dB |
| `--ground_effect` | −0.56 dB | −0.80 dB | −0.32 dB | −3.81 dB |

Most of the upwind loss comes from the low states near the end of the approach. From 10 m up,
the shadow begins about a kilometer upwind.

**The creeping-wave rate.** The shadow's decay rate (the least attenuated creeping-wave mode
over the ground's impedance) is held to the slope of the exact residue series deep in the
shadow to 0.3%.

### Shadow-transition calibration

The shadow transition replaces ray theory within a few diffraction lengths of a shadow
boundary ([THEORY §8.1](THEORY.md#81-linear-sound-speed-profile-circular-rays)). The evidence
behind its design and constants:

- **Ray theory near the boundary.** Its ground term at the grazing geometry was 10–25 dB low
  over grass and about 5 dB high over rigid ground there.
- **The table's reference source height.** V collapses over source height to within about
  1 dB, for sources from 5ℓ to 100ℓ (ℓ the Airy length), from ξ = −1.5 on.
- **On a linear lapse**, where the ray geometry is right, the boundary field takes the
  shadow-edge error against the PE from 5–20 dB to 0.1–1.0 dB (−0.018 /s, sources 30 and 100 m,
  125 Hz – 2 kHz).
- **On log laws** it needs the per-band boundary too. Over the stretch where the PE falls from
  −3 to −25 dB, the error against it (rms):

  | log-law winds, 16 cases, 1.2 m receiver | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz | pooled |
  | --- | --- | --- | --- | --- | --- | --- |
  | layered boundary, no transition | 7.1 | 22.6 | 34.6 | 37.4 | 39.3 | 30.7 |
  | boundary field and per-band boundary | 6.0 | 5.2 | 5.7 | 4.7 | 8.4 | 6.1 |

  Scaled by the layered boundary alone, the field falls 10–20 dB too fast at 125–250 Hz and the
  pooled error is 19.5 dB. The lapse control stays at 0.4–0.5 dB.
- **The boundary law's constant, 0.363**, was fitted on twelve log laws (u* 0.2–0.6 m/s,
  z₀ 0.01–0.1 m, sources 30 and 100 m, 125 Hz – 2 kHz) and checked on `example_windy.csv`, held
  out: 6.8 dB rms over the shadow edge against the PE, against 18.9 dB at the layered boundary.
  On the demo grid with a point receiver, `stratified` then averages 3.3 dB below `linear`
  (5th percentile −13.2 dB), against 6.1 dB (−19) with the layered boundary.
- **For a receiver on the ground** (3 mm up), over 16 log-law cases, the shadow-edge error is
  7.1 dB rms, against 10.0 for `linear`. By band, 125 Hz is the worst (10.2 dB) and 500 Hz the
  best (4.4). The constant that suits both heights together is 0.39, only 0.1 dB better.
- **Flush microphones take 0.75.** Fitted per case against the PE, a flush microphone's best
  boundary lies about 0.9 diffraction lengths (at the fit's curvature) past the layered one on
  the day-284 profiles, and 1.2 on flush log laws, where 0.363 put it at 1.7–2.7. 0.75 lies
  between the flush log laws' own optimum, 0.47, and day 284's, 0.94. A 1.2 m microphone keeps
  0.363 (0.387 when refitted with the day-284 profiles added).
- **Where the fit casts no shadow**, the fallback multiples 1.60 (flush) and 1.90 (1.2 m) times
  (f / 1 kHz)^−0.08 are the per-band optima pooled over the calibration set.
- **The refitted law against the PE** (shadow-edge rms, dB):

  | | before the flush refit | current law |
  | --- | --- | --- |
  | day-284 profiles, 18 cases (sources 50, 120, 250 m), flush | 12.9 | 8.1 |
  | log laws, 16 cases, 1.2 m | 6.2 | 6.0 |
  | log laws, 16 cases, flush | 7.3 | 7.6 |

  The best boundary chosen per case reaches 2.7 dB pooled. Every smooth law tried sits at about
  6–9 dB: the optimum is sharp where the shadow wall is steep, at high frequency.
- **The hand-over to the linear fit at 0.9 of the layered boundary.** The layered ray focuses
  short of its boundary into a spike of up to +10 dB the PE does not show. Taking the linear
  fit's ray from 0.9 on gives 0.6 dB rms against the PE over the last half of the way to the
  boundary on the day-284 profiles, and 1.7 dB on log laws (1.4 and 3.3 dB handing over at the
  boundary itself). Returning to the layered shadow ray past the fit's boundary stepped by 2 dB
  on a log law and 10 dB on day 284.
- **Cramped bands.** For a 30–50 ft source under a weak upwind gradient, every band to 500 Hz
  has its boundary within two diffraction lengths of the source. The table alone put the whole
  path 2–5 dB under the PE, which matched straight rays. Without the creeping-wave loss in the
  blend, a 30 m source in the linear lapse sat 6 dB under the PE at 250 Hz. Keeping the ray
  model into the shadow was worse deep in it (day-286 PE matrix, 50–160 Hz: 5.7 dB rms against
  5.1). The cramped rule costs the lapse case 0.1 dB rms at 200–630 Hz; day 284's PE cases do
  not move.
- **The day-286 PE matrix** past 1,000 ft, after the near-source boundary and cramped-band
  rules: `stratified` against the PE went from 9.5 to 7.5 dB rms at 200–630 Hz, from 11.7 to
  9.2 at 800–2,000 Hz, and from 5.4 to 5.1 at 50–160 Hz.
- **Plates.** The plate over bare ground, K, varies by under 0.1 dB from 0.5° to 10° of
  elevation near grazing. It is about 0 dB up to 1 kHz, +0.9 dB at 2 kHz, +2.6 at 4 kHz and +7
  at 8 kHz, and varies with azimuth by ±1 to ±2.5 dB.

**Caveat.** The PE band levels the constants were fitted against used an older band quadrature
that read 0.16–0.47 dB rms off the true band mean downwind (bias +0.1 to +0.35 dB) and 0.4–1.5
upwind. The constants stand within their tolerances, but anything fitted to PE band levels
should be rescored against fresh PE runs.

### Resolved shadow against the PE

`--shadow_model resolved` ([THEORY §8.3](THEORY.md#83-the-resolved-shadow)) was scored against
the PE before any flight data. The score is the rms of model minus PE over ranges of 300 ft and
beyond where the PE is within 20 dB of straight rays, per case and band, pooled. Bands are
125–4000 Hz, ranges to 4000 ft, ground variable porosity 200, with flush (3 mm) and 1.2 m
point receivers.

- **The synthetic set** (308 cases, generated and frozen before any variant was scored): linear
  lapses of −0.01 to −0.06 /s; neutral log laws through a 46 ft value of 2–12 kt, extrapolated
  below 46 ft; the same with +6 kt of shear from 46 to 200 ft; upwind, downwind and crosswind;
  sources of 5–200 ft.
- **The measured set** (194 cases, held out): R44 runs 229354, 229355, 229363, 229369 and
  229372 toward 140° (sources 5–80 ft); all 16 day-284 profiles toward their 46 ft upwind
  bearing (20–200 ft); and the calibration cases (`example_windy.csv` and three synthetic log
  laws, 30 and 100 m).

| class (case-bands) | legacy | resolved (series) | resolved (straight) |
| --- | --- | --- | --- |
| linear lapse (336) | 3.38 / −0.06 | 1.39 / +0.30 | 2.97 / +1.23 |
| log law, upwind (322) | 4.77 / −1.75 | 3.70 / +0.74 | 4.13 / +1.03 |
| log law + shear, upwind (167) | 6.73 / −2.71 | 4.25 / −0.20 | 4.82 / +0.16 |
| log law, crosswind (= adiabatic lapse) (336) | 4.01 / −0.89 | 0.73 / +0.04 | 2.66 / +1.18 |
| downwind (504) | unchanged | unchanged | unchanged |
| **synthetic, pooled (1833)** | **5.30 / −0.74** | **4.24 / +0.25** | 4.70 / +0.82 |
| R44 day 229 (300) | 6.58 / −2.81 | 5.19 / −0.47 | 5.70 / +0.28 |
| B407 day 284 (748) | 5.67 / −2.42 | 5.16 / +0.03 | 5.18 / +0.10 |
| calibration cases (96) | 4.28 / −1.93 | 3.02 / −0.18 | 3.02 / −0.17 |
| **measured, pooled (1144)** | **5.82 / −2.48** | **5.03 / −0.12** | 5.18 / +0.12 |

(rms / mean of model minus PE, dB.) Resampling whole profiles, the change in pooled rms is
−1.06 dB [−1.71, −0.61] on the synthetic set and −0.79 [−1.16, −0.41] on the measured one.

- **R44 run 229355**, 20 ft source, flush receiver, 1 kHz, at 1000, 1250 and 1500 ft: the PE
  gives −13.6, −17.3 and −20.1 dB; legacy −18.6, −28.2 and −28.9; resolved −14.2, −18.0 and
  −21.7. There the layered boundary is at 446 ft and the fit's at 1306 ft, and the PE's onset
  follows the fit's.
- **A linear −0.02 /s lapse**, 20 ft source, flush: resolved is within 0.3 dB of the PE from
  800 to 2000 ft at 500 Hz and 1 kHz; legacy is up to 17 dB low.
- **Decay from the resolved profile's tangent** (kept as `niceops_pe_compare --tangent_decay`)
  is worse on every upwind class: log law upwind 3.93 → 5.35 dB rms, log law with shear
  4.50 → 5.64, with the ray's cramped ground term.
- **Legacy with no fit shadow.** On day-284 run 200, 200 ft source, the legacy fallback
  boundary gave −27 dB at 2500 ft and 1 kHz, where the PE gives −5.6.
- **Where it is worse.** 58 of the 194 measured cases and 10 of the 308 synthetic ones worsen by
  more than 0.5 dB rms. Almost all are strong winds with a 1.2 m receiver: day 284's 20 ft
  sources go from 1.6–2.1 to 3.7–5.3 dB, and the 12 kt log laws worsen by up to 2.5 dB. There
  the resolved boundary lands 1–2 diffraction lengths too far out.
- **Beyond 4000 ft** (two day-284 profiles, 300 and 500 ft sources, to 8000 ft, flush): where
  the fit casts no shadow (run 200) the PE stays near straight rays and resolved is much closer
  (rms 19.5 → 13.8 and 14.0 → 7.9). In a true upwind shadow (run 205) resolved is 1.6 dB loud
  on average (rms 2.0 → 3.9 and 2.7 → 4.1).
- **Deep in shadows** (PE more than 20 dB under straight rays) neither model is right: the floor
  holds both 25–35 dB above a still-air PE.
- **`--shadow_phase` sensitivity:** 0.5 gives 5.05 dB on the measured set (day 284: 5.13); 2
  gives 5.28.

The flight-data scores are in [Refraction against flight data](#refraction-against-flight-data).

### Micro-shadows and skip zones

Under day 286's last profile (run 286288) the air toward 270° is fastest at the ground, because
the log-law headwind assumed below the LIDAR's 46 ft gate takes 0.2 m/s off in the lowest foot,
and an inversion near 330 ft makes an elevated duct above it. A 100 ft source's returning rays
turn up again 1–20 ft above the ground and reach neither a flush nor a 1.2 m microphone. At
8,000–13,000 ft those microphones stay on the shadow's floor, 11–21 dB under the PE. The PE
tunnels through the layer: flattening it moves the PE by under 1 dB to 250 Hz and about 6 dB at
1 kHz. Short of the duct's return, about 12,000 ft out, lies a skip zone no ray reaches, which
the PE fills from the caustic 1–2 kft early. Near a tunneled family's caustic, leaving the pair
of rays to the +10 dB focusing cap put the return 5–12 dB over the PE at 250 Hz – 2 kHz.

Three opt-in treatments against the PE (100 ft source, flush, 2,000–13,000 ft; `stratified`
minus PE, mean / rms dB):

| Day-286 PE matrix (five profiles) | 50–160 Hz | 200–630 Hz | 800–2,000 Hz |
| --- | --- | --- | --- |
| neither | −5.4 / 9.8 | −6.1 / 11.1 | −3.6 / 9.9 |
| `--tunneling` | −3.4 / 7.3 | −3.8 / 8.5 | −1.4 / 7.3 |
| `--skip_fill` | −4.5 / 8.9 | −5.6 / 10.5 | −2.7 / 9.2 |
| both | −1.7 / 5.1 | −2.3 / 6.5 | +0.3 / 6.7 |

Against measured day-286 frames where `stratified` predicts at least 10 dB under straight rays
(1.4% of frame-bands to 2 kHz; sources under 300 ft, 0.4–2° elevation, 4–9 kft), energy-mean
predicted minus measured, dB:

| | straight | `stratified` | `--tunneling` | `--skip_fill` | both | PE table (5 runs) |
| --- | --- | --- | --- | --- | --- | --- |
| 25 runs, 50–160 Hz | +3.8 | −10.5 | −7.3 | −3.6 | −2.1 | |
| 25 runs, 200–630 Hz | −0.1 | −15.5 | −10.7 | −7.2 | −5.0 | |
| 25 runs, 800–2,000 Hz | −7.6 | −20.6 | −20.6 | −20.7 | −18.9 | |
| 5 runs, sources under 150 ft at 3–6 kft, 50–630 Hz | +1.0 to +1.6 | −11.5 to −13.3 | | | −2.5 to −4.5 | +0.3 to +2.5 |

None moves the SEL: the exposure within 10 dB of each microphone's measured peak is unchanged to
0.02 dB under every option. What they leave:

- **The tails above 800 Hz.** The tunneled field at the ground, and the skip zone's edge short
  of the direct rays' reach, are still 10–20 dB low there.
- **Run 286276's far frames** (150–330 ft sources, 6–20 kft): every model but plain
  `stratified` is too loud, the PE table by 13.5 dB and straight rays by 8.8, so something the
  PE leaves out matters there.
- **Tunneling and skip-zone fill** improve some duct-shadow frames but can overpredict other far
  frames, which is why they stay optional.
- **The PE table** in this duct, read across 45° bearings and eight source heights, is within
  about 1 dB of a direct PE run off its nodes at 63–250 Hz, but 5 dB rms overall, almost all at
  1–2 kHz, where the return's caustic moves between nodes.

### Propagation tables

- **The plumbing check.** A table built from the straight-ray ground model reproduces the
  demo footprint under `--ground_effect --exact_ground` to −0.002 dB mean and 0.19 dB at worst
  with 10 ft source-height spacing, and to 0.02 dB at worst with 5 ft spacing. The ground effect
  itself moves cells by up to 4.3 dB. Over the ground model the interpolation error shrinks
  fourfold as the spacing halves.
- **A PE table through `atmospheres/isa.csv`** (21 bands from 20 Hz to 2 kHz, 23 source heights
  every 50 ft, ranges to 15,000 ft every 10 ft) covered 99.7% of the demo case's evaluations; the
  rest were steeper than the PE's 50° and took the ground model. Against the ray model's
  footprint it reads −0.32 dB mean (−0.9 to 0 dB, 5th to 95th percentile; −2.8 dB at worst):
  the standard lapse's mild upward refraction thins the field at long range.
- **Reciprocity.** Tables rebuilt both ways from the 2017 table study's recipe (9 source heights
  5–700 ft, flush and 1.2 m receivers, 8 bearings, 50 Hz – 4 kHz, ranges to 4000 ft every 50 ft).
  `--march receiver` minus `--march source`, by elevation, over cells above −40 dB (largest
  |Δ|):

  | run | below 5° | 5–10° | 10–20° | 20–35° | 35–50° |
  | --- | --- | --- | --- | --- | --- |
  | 229355 (R44, stable morning) | 0.013 dB | 0.001 | 0.003 | 0.025 | 0.11 |
  | 229369 (R44, unstable) | 0.008 | 0.001 | 0.003 | 0.023 | 0.10 |
  | 284200 (B407, day 284) | 0.013 | 0.003 | 0.005 | 0.026 | 0.11 |

  The rms over all covered cells is 0.003–0.010 dB. The 35–50° differences are the 1.2 m
  microphone at 63–80 Hz, where it sits on the second grid node, and the flush one's steep-angle
  residual (0.035 dB); the source-marched value is the closer to the exact field there. At the
  first range (a single self-starter step) the two differ by up to 0.29 dB at 50 Hz for the 5 ft
  source, and the receiver march is the closer to the exact field (6.1 dB, source march 6.5,
  exact 5.9). Still-air deep shadows below −40 dB re free field (2–10% of the cells, 2–4 kHz
  past 1 km) differ by up to 0.55 dB, because the absorbing layer follows each march's domain.
- **Turbulent tables** agree between the two marches only in the ensemble mean, to sampling
  noise.
- **Plate microphones through the ground-model table** against the plate run under
  `--ground_effect --exact_ground` (`gr1425_db225`, demo grid at 492 ft), by the table's lowest
  receiver row:

  | lowest row | mean | worst |
  | --- | --- | --- |
  | 0 (flush) | −0.001 dB | −0.012 dB |
  | 1.2 m | −0.003 | −0.13 / +0.11 |
  | 5 m | +0.008 | −0.46 / +0.60 |
  | 30 m | −0.024 | −1.28 / +0.55 |

  The raised rows' error does not shrink with a twelve times finer elevation grid for K. With K
  taken at 2° for every pair, the flush row was −0.014 dB mean (−0.34 to +0.08 dB).

- **Over terrain**, a ground-model table read in each pair's mean-plane frame against the terrain
  ground model itself (the synthetic hills grid, demo track, 492 ft grid, variable-porosity ground
  of 200 kPa·s/m²):
  - **Plates**, through the table's row at 0 plus K: 0.014 dB at worst (0.0016 dB rms) with
    source heights and ranges every 20 ft, and 0.0046 dB (0.0004 rms) every 10 ft.
  - **Point receivers** read between rows at their own height above the plane, 0 to 20.5 ft
    there. From 20 to 250 Hz, with rows every 2, 1 and 0.5 ft, the worst cell is 2.26, 0.83 and
    0.17 dB off (rms 0.31, 0.11 and 0.04 dB), so the error is the rows' interpolation. Over the
    full band, with rows every 4, 2 and 1 ft, the worst cell is 1.67, 1.92 and 1.96 dB off (rms
    0.37, 0.28 and 0.18 dB): the worst does not fall at these spacings, because the ground's
    interference fringes are λ/(2 sin ε) apart in receiver height, finer than the rows in the
    mid bands at steep elevations ε.
  - **On the full demo grid** (14,400 cells), through the five-row table the cost measurements
    use (rows 0, 1.2 m, 8, 16 and 25 ft, source heights and ranges every 20 ft): point receivers
    0.39 dB rms and 2.27 dB at worst, plates 0.003 dB rms and 0.17 dB at worst, with 8 cells
    over 0.05 dB.

  This checks the geometry and the plumbing.
- **Over measured terrain**, at Sweetwater in the altitude-variation test ([Data sets](#data-sets)):
  eight runs chosen to span its weather (three upwind, two with shadowing profiles, a calm
  inversion, two downwind), one PE table per run from its own LIDAR and balloon profile (44.7 Hz
  to 2.2 kHz, straight rays above), the 1 m 3DEP grid, flush boards on variable-porosity ground
  of 200 kPa·s/m². SEL, predicted minus measured:

  | | flat datum, layered rays | flat datum, PE table | terrain, PE table | terrain, straight rays |
  | --- | --- | --- | --- | --- |
  | mic 1 − mic 21 on the 3.4° slope (measured −2.54 dB), error | −2.09 | −0.70 | −0.70 | −0.28 |
  | mics 1–21, mean / rms | +1.40 / 1.66 | +1.49 / 1.77 | +1.55 / 1.85 | +1.55 / 1.84 |
  | mics 27–29, 300 ft below the datum: mean (MAE), 24 mic-runs | −8.76 (12.30) | −2.43 (5.89) | −0.32 (5.26) | +4.11 (4.49) |
  | the same, three upwind runs | −19.49 | −0.23 | +2.04 | +4.65 |
  | the same, the two shadowing profiles | −12.66 | −14.62 | −10.77 | +2.58 |

  Read over terrain, the PE table puts the off-datum microphones' source 280–310 ft above
  their paths' planes rather than 250 ft above the datum, a shallower shadow: 1–4 dB louder
  there, the best mean of the four. It cannot remove a shadow the table itself predicts and the
  measurements do not show (the two shadowing profiles, still 8–14 dB low). On the slope it
  changes nothing beyond what terrain does to straight rays: a flush board reads the table at
  the source's height above its plane, about the same as above the datum. No path's plane tilts
  past 3.7°, so the refracting-table warning stays silent, as it should. A ground-model table
  over the same terrain is the ground model to 0.04 dB. Eight runs; the run-to-run spread is
  large.
- **Turbulence and roughness built into a ground-model table** reproduce `--ground_effect
  --exact_ground` with the same settings on the plumbing check's grid: Chessell's 4e-4 rad·s/√m
  with a 0.3 ft roughness to 0.0001 dB mean and 0.021 dB at worst, HARMONOISE's γ<sub>T</sub>
  3e-5 to −0.001 dB mean and 0.077 dB at worst. Against the still ground model they move cells
  by −1.5 to +5.9 dB and −0.2 to +3.1 dB, filling the interference nulls.

A table checked with `--validation_csv` is compared with direct solves at cell midpoints. That
comparison includes both interpolation error and differences in resolution and band quadrature;
it does not bound the error everywhere, and its reference is not a measurement.

### Ground and plate tables

- **Point receiver.** The tabulated ground gain (641 angles, 40 points per decade of range)
  gives level histories within 0.001 dB of the exact model. On the demo grid 99.97% of SELs
  agree to 0.01 dB; the rest are points where the 10 dB-down window's first-to-last rule flips
  on a secondary hump at the threshold.
- **Plate receiver.** Interpolating its reflection coefficients from tables gives SEL within
  0.0005 dB rms of the exact evaluation.
- **The plate's band average.** Taking Q quadratic in log frequency through each band's edges and
  center is within 0.07 dB of evaluating it at every sub-frequency wherever the level is above
  −15 dB re free field, and up to about 1.5 dB off in deep near-grazing nulls 30–70 dB down.
- **Faddeeva fast paths.** The asymptotic series and the Taylor table agree with Weideman's
  series within 1e-11.

### Ground classes

`--ground_classes` on the demo track, over 14,400 1.2 m pole microphones at the demo grid's
points, every other one on Delany–Bazley concrete (10,000 kPa·s/m²) and the rest on
variable-porosity grass (200 kPa·s/m²), against the same microphones all on grass:

- **The grass microphones are unchanged to the bit**, tabled and exact, so a microphone's level
  depends only on its own class.
- **The tables against `--exact_ground`**, mixed classes: 0.009 dB at worst (0.003 dB on the
  concrete microphones), the same as with grass alone.
- **What concrete under a pole does**, SEL on concrete minus SEL on grass at the same
  microphone, exact, by horizontal distance from the ground track:

  | distance | under 500 ft | 500–1000 | 1000–2000 | 2000–3000 | 3000–5000 | over 5000 |
  | --- | --- | --- | --- | --- | --- | --- |
  | median | +0.74 dB | +1.28 | +2.41 | +4.11 | +5.13 | +6.33 |
  | range | +0.49 to +1.93 | +0.54 to +3.83 | +0.76 to +6.87 | +0.84 to +8.49 | +1.13 to +7.95 | +1.74 to +7.89 |

  Each microphone is evaluated as if the whole ground were its class, so this is the effect of
  a concrete plain. A pad of a few meters under a pole has less, most of all far from the track,
  where the reflection lands tens of meters out on whatever surrounds it.

**Against measurements**, two data sets with microphones on different surfaces:

- **The R66 static ring at Eglin** (2017, day 231, 18 hovers): boards 1 and 13 stand on cells the
  2021 ESA WorldCover land cover calls built-up along their whole path to the hover point, 3–11
  on grass. In ground effect (heard 0–2° below the horizon), relative to the median of boards
  4–11, boards 1 and 13 measure +9 to +11 dB louder than the grass plate predicts at 500 Hz to
  1 kHz (+6 to +8 at 1.6 kHz), at every heading. Taken as point receivers at the board's top over
  a hard class:

  | class for boards 1 and 13 | 500 Hz | 1 kHz | 1.6 kHz | 2.5 kHz | SEL, predicted − measured |
  | --- | --- | --- | --- | --- | --- |
  | measured excess, board 1 / 13 | +10.4 / +9.1 | +9.1 / +9.9 | +7.6 / +7.9 | +6.1 / +3.5 | |
  | grass plate (the model without classes) | 0 | 0 | 0 | 0 | −3.4 / −1.6 |
  | Delany–Bazley 20,000 kPa·s/m² (concrete) | +15.9 | +19.0 | +17.8 | +16.7 | +6.0 / +7.9 |
  | Delany–Bazley 2000 | +13.9 | +12.6 | +7.4 | +5.2 | +3.3 / +5.2 |
  | variable porosity 1000, α 100 | +11.9 | +9.9 | +6.2 | +4.7 | +2.0 / +3.9 |

  A hard class is right in kind and a class of about 1000–2000 kPa·s/m² fits 500 Hz–2.5 kHz to
  about 1.5 dB, but concrete overpredicts by 8–10 dB: "built-up" there is not a concrete plain.
  No uniform hard class fits 125–250 Hz (measured −5 to +1 dB, every hard class +4 to +7), and
  board 12, grass for its first 20 m and built-up beyond, measures soft: at grazing the ground
  next to the board and the mixed path decide, which a uniform class cannot represent. Out of
  ground effect (31–80° down) every class is within 2 dB in the bands and 1 dB in SEL, and the
  measurements cannot tell them apart. A ±5 ft error in the hover height moves the concrete
  prediction at 1 kHz by ±5 dB and the best-fitting class from about 1000 to 2000.
- **The T-38A runway test at Wallops** (NASA TP-1747): pairs of microphones on a concrete runway's
  centerline and on grass beside it, at 0, 1.2 and 9.14 m, five reference runs transcribed from
  its appendix C (100 Hz–2 kHz). Each pair's concrete-minus-grass difference against NICEOPS's,
  source spectrum cancelled, variable porosity 200 for grass and Delany–Bazley 20,000 for
  concrete, rms over bands: 1.2 m pairs near the track 3.3 dB (5.0 predicting no difference),
  far 7.2 dB (10.4), flush pairs 5.1 dB, 9.14 m pairs 5.9 dB. Willshire's own fits (62.5 and
  750) score worse in every group. About 6 dB of the scatter is not ground effect: the 9.14 m
  pairs, which should differ by under 1.5 dB, differ by −12 to +7 dB from run to run, so these
  data show the classes help at 1.2 m but cannot place the uniform class's limit.

### Terrain screening

`--terrain_screening`'s two diffraction terms, CNOSSOS-EU's 10 lg(3 + 40δ/λ) and the knife
edge's −10 lg |F(v)|², against Sommerfeld's exact solution for a rigid half-plane lit by a plane
wave, at kr = 2000 from the edge, at the same path difference (`Delta_dif against the exact
half-plane` and `the knife edge against the exact half-plane`, `tests/test_screening.cpp`):

| Angle into the shadow | 0° | 0.5° | 1° | 2° | 5° | 10° | 60° |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Exact − CNOSSOS, dB | +1.17 | +2.39 | +2.67 | +2.52 | +2.27 | +1.96 | −0.67 |
| Exact − knife edge, dB | −0.08 | −0.11 | −0.15 | −0.23 | −0.50 | −0.94 | −3.62 |

CNOSSOS-EU reads low by about 1.2 dB at grazing and by up to 2.7 dB just inside the shadow, so a
screened level errs high there. The knife edge is within 0.25 dB of the exact loss to 2° into the
shadow and deeper beyond, by 0.9 dB at 10°: the Fresnel-Kirchhoff edge is black, where the rigid
half-plane's faces also send sound into the shadow. Its band loss, read from its table, is within
10⁻⁸ dB of the quadrature it tabulates for δf/c from 10⁻⁶ to 10³ and within 10⁻⁷ dB to the table's
ends (`the knife edge's band loss table against direct quadrature`).

Against the CNOSSOS-EU conformance cases of ISO/TR 17534-4:2020 (`[iso17534]`,
`tests/test_screening.cpp`), `cnossos` on the three a vertical-plane terrain model can take, at
the TR's c = 340 m/s and homogeneous conditions:

- **TC23**, an earth berm 5 m high between a source 1 m and a receiver 4 m up, 70 m apart: the
  rubber band's two edges, every diffraction term (Δdif(S,R) uncapped to 29.05 dB at 8 kHz,
  Δdif(S′,R), Δdif(S,R′)), both side planes' heights and lengths, and, with the TR's own
  A_ground injected for NICEOPS's ground model, both Δground terms and A_dif (the 25 dB cap
  included), all within 0.01–0.03 dB or 0.01 m of the TR.
- **TC06**, a path clearing a plateau's edge by 0.67 m: unscreened here, where the TR diffracts
  its 500 Hz and 1 kHz bands (A_dif 4.31 and −0.83 dB); that is NICEOPS's strict onset, and the
  Directive's term at those bands is reproduced to 0.01 dB.
- **TC05 and TC20**, unobstructed paths over a 10 m rise: the whole path's mean plane (the
  geometry a table read over terrain takes) is the polyline's exact least-squares plane to
  1 mm. TC05's tabulated values agree to 0.01 m; TC20's, as NoiseModelling transcribes them,
  differ from the exact plane by up to 0.04 m.

The TR's sources are about 1 m up, so these test the method's arithmetic, not an elevated
source.

Behind the synthetic ridge (25 m high, 250 m to the track's side) on the full demo grid, with
1.2 m pole microphones on variable-porosity ground of 200 kPa·s/m², `cnossos` lowered 3044 of the
14,400 cells' SEL, 470 by more than 1 dB, 21 by more than 6 and two by more than 10 (17.6 dB at
most), and raised none. 9.6% of the evaluations were screened, 39% of those reached the 25 dB cap
in some band, and the ground weight F never needed clamping. `kirchhoff` lowered the same 3044
cells, 495 by more than 1 dB, 28 by more than 6 and two by more than 10 (14.1 dB at most), within
−0.7 and +3.5 dB of `cnossos` (a median 0.07 dB deeper), and raised none. Over the hills fixture
the full grid screened nothing and was the unscreened run byte for byte.

**At Sweetwater**, behind an isolated hill, in the AS350 altitude-variation test
([Data sets](#data-sets), with recordings 28 and 29 where they were made): all 166 runs (4088
mic-runs after the study's exclusions), SD1 database, flush boards on variable-porosity ground of
200 kPa·s/m², the 1 m 3DEP grid, the default minimum protrusion. Only position 28 is screened
behind the hill (19,439 evaluations in 119 runs, from the record's start to x ≈ −3216 ft), and
there the blocked states drop by a median 12.6 dBA with `cnossos` and 10.8 dBA with `kirchhoff`;
27 and 29 graze but are never screened, and the on-datum array's scores are unchanged. SEL and
LAmax, predicted minus measured, mean [95%] / rms, and the change from the unscreened run:

| position 28 | terrain, unscreened | `cnossos` | `kirchhoff` |
| --- | --- | --- | --- |
| SEL | +4.98 [4.35, 5.60] / 6.32 | +4.48 [3.85, 5.10] / 5.92 | +4.50 [3.84, 5.13] / 5.94 |
| SEL change | | −0.50 [−0.57, −0.44], rms −0.40 | −0.48 [−0.54, −0.41], rms −0.38 |
| LAmax | +0.13 [−0.59, 0.84] / 4.31 | −0.17 [−0.86, 0.51] / 4.25 | −0.17 [−0.87, 0.49] / 4.25 |
| LAmax change | | −0.31, rms −0.06 | −0.31, rms −0.06 |

The blocked segment comes early in the pass, at long range, so screening moves SEL by half a
decibel; the over-prediction left is shared by every off-datum microphone (27 and 29 read +2.6 and
+3.9 dB). Inside the blocked segment, against the rest of the window, the measured level is
−4.3 dB [−5.1, −3.3], the unscreened prediction +2.4, `cnossos` −9.7 and `kirchhoff` −7.7 (27
and 29 measure −0.2 and +0.2): the hill does shield position 28, and the screening overstates it
by about 5.4 dBA with `cnossos` and 3.4 dBA with `kirchhoff`. By band, the same contrast measured
against the unscreened prediction, median over the 119 runs:

| band | 63 Hz | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz |
| --- | --- | --- | --- | --- | --- | --- |
| position 28, measured | −3.6 | −8.2 | −6.0 | −7.1 | −4.8 | −4.3 |
| positions 27 / 29, measured (not screened) | −0.9 / −2.5 | −1.6 / −3.4 | −0.1 / −1.0 | −1.7 / −2.1 | −2.0 / −2.1 | −1.4 / −3.0 |
| position 28, `cnossos` | −4.8 | −4.9 | −10.2 | −25.5 | −37.9 | −39.4 |
| position 28, `kirchhoff` | −7.9 | −8.6 | −9.5 | −10.6 | −11.9 | −13.1 |

Net of the unscreened positions, the measured shielding is about 2–6 dB and nearly flat in
frequency. `cnossos` is right to within a few decibels up to about 200 Hz and too deep above it,
by more than 30 dB at 1–2 kHz. `kirchhoff` deepens by only about 1 dB per octave above 250 Hz;
it is within 4.5 dB of position 28's measured row up to 500 Hz and 7–9 dB too deep at 1–2 kHz. A
vertical-plane method leaves out both sound reaching the microphone around the hill's flanks and
turbulent scattering into the shadow; which of them fills the shadow here is not settled. The
same runs found the 1 m grid's texture 3–5 m from a flush board screening paths at 1–2°
elevation, by 0.03–0.14 m: `--screening_min_height` (0.5 m) now leaves those unscreened.

### Urban receivers

The ideal canyon's claim that a plane wave reaches a point by at most four specular paths -- from
above and off the floor, each with the wave's horizontal direction or its mirror -- and which
ones, holds against a brute-force image-source search over 3,000 random canyons (widths 8 to
38 m, walls 3 to 43 m, sources 10^8 m away; the floor's paths are those to the microphone's floor
image through the walls' images below it; `tests/test_urban.cpp`). The paths agree in number and
order every time; three or four arrive in over 100 cases, the direct ray is blocked in over 100,
and only the floor's rays arrive in over 50. The first-order reflections the run evaluates agree
with the same search in all 4,000 cases tried.

What evaluating the first-order reflections explicitly buys, for an omnidirectional source over a
20 m canyon with 15 and 18 m walls, a microphone 6 m from one wall at 1.5 m, walls reflecting 0.8
and a rigid floor (specular paths only, 1 kHz; the hidden Catch2 case `[.urban_separability]`).
The error is in the total energy at the microphone against every specular path of the point
source, from above and off the floor, as the absolute level difference over 4,000 random
directions above 5 degrees:

| Range | Explicit reflections + table: median / p95 / max, dB | Everything in the table: median / p95 / max, dB |
| --- | --- | --- |
| 30 m | 0.36 / 4.2 / 43 | 1.82 / 9.1 / 46 |
| 60 m | 0.21 / 3.9 / 40 | 0.71 / 7.1 / 40 |
| 120 m | 0.10 / 4.2 / 37 | 0.25 / 5.0 / 37 |
| 300 m | 0.04 / 2.2 / 10 | 0.07 / 2.7 / 10 |
| 1 km | 0.01 / 1.4 / 4.0 | 0.01 / 1.4 / 4.0 |
| 3 km | 0.00 / 0.21 / 3.7 | 0.00 / 0.28 / 3.7 |
| 10 km | 0.00 / 0.04 / 3.6 | 0.00 / 0.04 / 3.6 |

The largest errors are directions near a wall's top, where the plane wave and the point source
disagree on whether a path is open; with the diffraction left out, as here, that is a step of
tens of dB. A directional source would widen the gap between the two columns, since each
reflection reads the sphere in its own direction; that is not measured yet.

That table is the specular rays alone (`--no_diffraction`), so that the comparison is with a
lattice of specular paths.

On the demo track (B407 database), seven microphones in ideal canyons
(`tests/data/urban_canyons.csv`) move by −1.0 to +2.2 dB in SEL against open ground with the
roof edges diffracting: +2.2 dB where the helicopter passes about 550 ft up along the street,
reflections adding to an open direct path, and −1.0 dB at the most shielded. Nothing here is
validated against
measurements in streets; none of the measured data sets ([Data sets](#data-sets)) has buildings
near the microphones.

#### The point-source model in a run

`--buildings` against the plane-wave table (`--urban_table`, built by `niceops_urban_table` for the same district and microphones): the demo
track over the eleven-building test district (`tests/data/urban_buildings.csv`), B407 database,
seven microphones 4 ft up; the table run on ground of 10^6 kPa s/m^2, which the two models' open
references agree on to 0.16 dB. The two agree in SEL to 0.2 dB at every microphone (behind a
building −6.88 dB against −6.92, re open ground), but not state by state, within 10 dB of each
microphone's peak:

| Microphone | SEL, point − table | LAmax | Per state: rms, p95, worst |
| --- | --- | --- | --- |
| In a street | −0.03 | −2.30 | 0.94, 2.39, 3.26 |
| At an L-shaped building's corner | +0.04 | −0.02 | 0.15, 0.34, 0.63 |
| At an intersection | −0.18 | +0.19 | 1.67, 4.01, 5.32 |
| Beside a tall wall | +0.02 | −0.24 | 0.29, 0.43, 2.35 |
| In the open | +0.16 | +0.09 | 0.19, 0.32, 0.46 |
| In an alley | +0.15 | −0.17 | 4.21, 10.9, 11.7 |
| Behind a building | +0.09 | +0.45 | 0.34, 0.53, 2.71 |

The large per-state differences come at 60–150 m slant range, and over a moving pass they cancel
in the exposure; over a hover or a slow pass by the same buildings they would not. On this track
roofs move none of the seven SELs, the microphones standing below every roof.

**What the run's sampling costs** (THEORY §6.5), on the same demo track and district, against the
same run unsampled:

- **The states sampled** (`--state_tolerance 0.02`, the default) against every state evaluated:
  12,453 of the 45,857 pairs evaluated. SEL within 0.003 dB and LAmax within 0.007 dB at every
  microphone. State by state within each microphone's 10-dB window, the 99th percentile is 0.13 dB
  at worst; one state in the alley is 2.1 dB off, a feature narrower than the stretch's middle
  sees. With one building 1000 km off and the floor absorbing, the sampled run is the run without a
  ground model to 1e-12 dB.
- **The shared faces sampled** against their integral at every wavenumber: no state's level moves
  by more than 0.002 dB, nor any SEL by 0.0001. Ray by ray, on long blocks 1 to 10 m thick with
  sources low, grazing, in the roof's plane and just above it and microphones through its plane,
  from 20 Hz to 10 kHz, within 0.07 dB and 0.012 rad; the samples an octave apart left 0.13 dB.
  Where the transition parameter is 30 or more the cascade alone is within 0.04 dB and 0.02 rad of
  the integral; under 0.1 it is up to 16 dB off, and with an end in the roof's plane its sign is
  wrong. `niceops_point_field`'s band levels over the test district stay within 0.1 dB of the
  integral's.
- **Cost**, 14 threads on an M4 Pro (the development build, asserts on): 3.2 s with `--no-gpu`,
  1.2 s with the rays formed and summed on the GPU, against 0.3 s for the run without buildings.

**A ground model's floor** (`--buildings --ground_effect`, each bounce taking Q). Not validated
against measurement among buildings; what is checked is that it is the open-ground model where it
should be. Over open ground, with Q held at the band's centre, the ray sum is the open-ground
model's band average to 1e-5 dB, for grass and concrete and sources 3, 20 and 120 m up; taking Q at
each sub-band's middle instead moves a band by up to 0.12 dB (a source 3 m up over hard ground,
where Q changes most across a band). With one building 1000 km off, on the demo approach with every
phase known, the SELs are the open-ground run's (`--exact_ground`):

| Microphones | SEL, buildings − open ground |
| --- | --- |
| 4 ft up, grass (`--flow_resistance 225`) | −0.003 to +0.060 dB |
| 4 ft up, concrete (`--ground_classes`) | +0.012 to +0.141 dB |
| Flush on GR1425 plates (`--plate_table`) | −0.003 to 0.000 dB |

The rest is the floor's ray reading the sphere along its own launch, where the open-ground model
reads it along the direct ray for both. On the test district a grass floor takes 0.8 to 4.0 dB off
the rigid floor's SELs, the flush microphones most, for 5% more time; the GPU's SELs are the CPU's
to 1e-4 dB.

**Turbulence in the buildings' shadows** (`--scattering` with `--buildings`). Not validated against
measurement: EDGE's room is still air, and DigiCity's weather is not known. What is checked is the
implementation. The fit reproduces Forssén et al.'s own Figure 2 (their scattering integral at
1 kHz, C_v² = 1, d_S = 40 m) at its corners to 2 dB, less the 6 dB of the ground it doubles at both
ends; the screen's height and distances are the hand-worked ones for a box under a sloping line;
and in a run, behind a 10 m wall 15 m off with the source 200 m away and 60 m up, γ_T = 10⁻⁵
(HARMONOISE's strong, C_v² = 0.16) adds 0.07, 0.26 and 0.52 dB at 125, 500 and 2000 Hz, the total
under open ground's and a microphone in sight untouched to the bit. On the demo approach the same
strength moves no SEL by 0.003 dB, every microphone seeing the helicopter at some point, for under
1% of the run's time.

**Facades that scatter** (`--facade_scattering`). For a source the same in every direction, the
run's band levels are `niceops_point_field --scattering`'s to 1e-9 dB, rays and diffuse field
together, with one building in the diffuse field's reach and one beyond it; on the GPU, to 1e-6 dB
(`tests/test_urban_point.cpp`). A source whose power falls as cos² off a wall's normal lights the
wall's field as Lambert's law says to 0.5%, which checks that each lighting path reads the source
along its own launch. On the demo track and district, s = 0.3 on every wall and roof:

| Microphone | SEL, s = 0.3 − specular |
| --- | --- |
| In a street | +1.03 dB |
| At an L-shaped building's corner | +0.08 |
| At an intersection | +0.54 |
| Beside a tall wall | +0.07 |
| In the open | +0.00 |
| In an alley | −0.22 |
| Behind a building | +0.00 |

The states sampled agree with every state evaluated to 0.003 dB, as without scattering, and patches of
20 ft agree with the default 4 m to 0.03 dB. Cost: 1.7 s on the GPU against 1.5 s specular, and
3.8 s with `--no-gpu` against 3.4 s, each microphone's field built once for the run.

### Roof-edge diffraction

The UTD wedge kernel (`src/wedge`, `tests/test_wedge.cpp`), a unit plane wave on a rigid wedge,
geometrical optics plus the diffracted wave, against exact solutions, over every direction and
four incidences:

| Wedge | Reference | k r = 10 | k r = 40 | k r = 200 |
| --- | --- | --- | --- | --- |
| Half-plane (n = 2) | Sommerfeld's solution | 4e-7 dB | 2e-6 dB | 6e-6 dB |
| Right angle (n = 3/2), a roof edge | Eigenfunction series, 40 digits | 0.016 dB | 0.005 dB | |

For a half-plane under plane-wave incidence UTD's Fresnel form is the exact solution; the
right-angle series was summed in Wolfram Language and reproduces a flat plane's incident and
reflected waves exactly at n = 1, where the kernel diffracts nothing.

In the canyon (`tests/test_urban.cpp`): swept across the direct ray's shadow boundary in 0.05
degree steps at 1 kHz, the coherent field never steps by more than 0.1 dB, where the specular rays
alone jump by an infinite amount. A metre clear of every boundary, at 2 and 4 kHz, it is the
specular rays' energy to 0.23 dB on average over 106 direction-bands, 1.5 dB at worst: the
interference a third-octave band leaves between equal rays a few metres apart. On a boundary it
departs by more, rightly -- at 70 degrees over a 20 m canyon a floor-reflected ray strikes a wall
2 cm below its top, and the diffracted ray halves it.

### Canyon wave solver

The solver, its normalization and the far-field transform together, over flat rigid ground: a
recording line 3 m up transformed as the canyon's opening is, against the image solution
4 cos^2(k z sin phi) for a source 1.5 m up, band by band from 40 to 160 Hz
(`the open-ground far field from a recording line is the image solution`,
`tests/test_canyon_wave.cpp`). From 45 to 135 degrees the bands agree within 0.1 dB at 12 points
per wavelength and 0.03 dB at 20; at 15 to 30 degrees the line's own truncation errs by 0.2 to
0.8 dB whatever its length, which the canyon, whose transform runs over the opening alone between
rigid roofs, does not share. Moving the grid's edges (roof and sky margins of 4 and 10 m) moves a
canyon's bands by under 0.1 dB.

Convergence, a 20 m by 20 m canyon, the microphone 6 m from wall A at 1.5 m, walls absorbing 0.05
(the hidden case `[.canyon_wave_study]`): 15 against 20 points per wavelength agree within 0.3 dB
in every band and direction tried, 10 within 1.4 dB (at 500 Hz, 15 degrees). Each run went the
full 4 s and ended 28-29 dB down.

Against the ideal canyon's geometric table, with the floor's rays, the roof edges' diffraction,
the same walls and the first-order reflection the run adds:

| Direction from +y | 15° | 30° | 45° | 60° | 75° | 90° | 120° (shadow) | 150° |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 500 Hz, wave, dB | −2.39 | −3.65 | 0.13 | 2.32 | 0.06 | 0.04 | −10.86 | 0.25 |
| 500 Hz, geometric, dB | −2.11 | −3.53 | 0.16 | 2.35 | −0.18 | −0.08 | −10.81 | 0.73 |
| 250 Hz, wave, dB | 5.75 | −4.32 | 0.93 | 2.72 | 0.45 | 0.01 | −9.44 | 1.45 |
| 250 Hz, geometric, dB | 6.10 | −4.25 | 1.19 | 2.91 | 0.47 | −0.03 | −9.21 | 2.10 |

At 250 and 500 Hz (kW = 92 and 183) the two agree within 0.65 dB in every direction, the shadow
included. At 63 and 125 Hz they differ by up to 7.6 dB, the canyon's resonances, which is the
band the wave tables are for. What each piece of the geometric table is worth: without the
floor's rays it reads −24.6 and −20.4 dB at 15 and 30 degrees, and with the CNOSSOS-EU
Delta_dif over the blocking edge in place of UTD, −14.46 dB in the shadow.

**Rotor tones against band tables** (the stage-2 gate; the hidden case `[.canyon_tone_study]`).
Canyons 20 m wide with H/W = 0.5, 1 and 2, microphones at the centerline and at W/4, 1.5 m up,
walls absorbing 0.05, solved to 600 Hz at 12 points per wavelength. Main-rotor (27.5 Hz) and
tail-rotor (79.3 Hz) harmonics to 460 Hz, B407-like, each Doppler-swept (M = 0.175) over 18
level passes (60, 150 and 300 m up; along, oblique to and across the street; overhead and 150 m
aside). The error is in each harmonic's SEL change against open ground, the band table read as
the run reads it (bilinear on the builder's 72 by 33 grid) against the coherent canyon field:

| H/W | Microphone | Per harmonic and pass: median / p90 / max, dB | A-weighted, equal harmonics: median / max, dB |
| --- | --- | --- | --- |
| 0.5 | centerline | 0.17 / 0.68 / 3.9 | 0.05 / 0.15 |
| 0.5 | W/4 | 0.24 / 0.85 / 2.3 | 0.08 / 0.35 |
| 1 | centerline | 0.19 / 0.77 / 3.3 | 0.04 / 0.49 |
| 1 | W/4 | 0.26 / 1.03 / 2.0 | 0.06 / 0.43 |
| 2 | centerline | 0.19 / 0.71 / 2.6 | 0.07 / 0.80 |
| 2 | W/4 | 0.36 / 1.12 / 3.3 | 0.12 / 0.47 |

Band energy is therefore adequate for SEL: under 1 dB for the A-weighted level everywhere,
under about 1 dB for nine harmonics in ten. A single low harmonic can be off by 2-4 dB, which
matters for narrowband or tonal metrics and not for SEL. Without the Doppler sweep the p90
changes by −0.1 to +0.5 dB. Most of the error is the table's band averaging: with the coherent
open ground in place of the run's band-averaged one, the p90 moves by at most 0.2 dB.

On the demo track the wave band moves the seven test canyons' SEL changes against open ground by
at most 0.08 dB with `--wave_below 400`, the wave solution and the diffracting geometric table
agreeing across the crossover, and by under 0.03 dB with `--wave_below 80`, where A-weighting
leaves little.

**On the GPU** (`niceops_urban_table --gpu`, single precision, an M4 Pro). Against the march in
double: the band ratios of a 12 m canyon with walls of absorption 0.05 and 0.2, 25 to 200 Hz,
within 7.6e-6 dB, and open ground's gain within 5.7e-6 dB, each run stopping at the same step;
a 20 m canyon's table to 707 Hz within 3.4e-5 dB, and fourteen 20 m canyons' to 354 Hz within
1.6e-5 dB; the seven test canyons' SEL within 0.001 dB. Cost: the 20 m canyon to 707 Hz at 15
points per wavelength is 1.08 million cells and 66,764 steps, 10.8 s of GPU time, against about
67 s on one core; the fourteen canyons to 354 Hz take 22 s with `--gpu` against 34 s on 14
threads.

### Buildings

The general builder (`src/building_paths`, `tests/test_building_paths.cpp`) is checked against the
ideal canyon and against the wedge, and its pruning against itself:

- **The ideal canyon, from two buildings 40 km long.** With diffraction off, the specular table
  equals the canyon's to 1e-12 (relative) for three canyons over 48 azimuths and 6 elevations;
  with it on, to 1e-9 dB in three bands over 54 directions, through 40 reflections either way.
  `building_gains` equals `ideal_canyon_gains` through the compass mapping to 1e-9 dB.
- **Oblique incidence.** At a roof edge and at a vertical corner, each a right-angled wedge, the
  field from the edge's own rays equals the 2-D wedge's at k sin(beta0) to 1e-9 (relative), for
  three microphones and 28 directions each.
- **The beams lose nothing.** In a district of ten buildings, one of them skewed, the rays found
  with and without beam pruning, to three diffracted reflections, are the same rays (64
  directions, three microphones).

Convergence and cost on a synthetic district: 100 blocks 60 m square with roofs from 10 to 40 m,
streets 20 m wide, every facade reflecting 0.9 and the floor rigid, the microphone 1.5 m up in mid
street; 72 by 33 directions, three bands, one thread. The last column is the 95th percentile of
the change from the most thorough table, in the cells within 10 dB of open ground:

| Specular / diffracted reflections | Images | Time | Energy mean | p95 change |
| --- | --- | --- | --- | --- |
| 8 / 2 | 63 | 0.16 s | −0.67 dB | 13.9 dB |
| 64 / 2 | 63 | 0.14 s | −0.46 dB | 3.4 dB |
| 64 / 4 | 513 | 0.36 s | −0.36 dB | 1.4 dB |
| 64 / 6 (default) | 2763 | 1.01 s | −0.34 dB | 0.52 dB |
| 64 / 7 | 5969 | 1.80 s | −0.34 dB | 0.26 dB |
| 64 / 9 | 25597 | 6.03 s | −0.34 dB | (reference) |

Grazing directions need the specular bounces: a ray at 1.5 degrees climbs half a metre per
crossing of a 20 m street. Walls this reflective are the slow case; the energy mean has settled
by four diffracted reflections. Beam pruning made a 27-band table with three diffracted
reflections 25 times faster, 44 s to 1.8 s, with the same rays, and the plan-view grid behind
the visibility tests (identical to trying every building and facade on 4000 random segments and
rays) a further 2.8 to 4 times.

**How much SEL depends on the scattering coefficient** (`tools/urban_scattering_study.py`, 31
minutes on ten threads). Grids of 60 m blocks with 20 m streets -- roofs at H/W 0.5, 1 and 2, and
of mixed heights -- aligned with the demo approach and turned 45°, placed under it and 1000 ft
aside where it is near 400 ft and near 150 ft up; in each, five microphones: on the centrelines of
the two streets mid-block, 2 m from a wall, a flush plate 1 m from a wall, and at an
intersection. Walls absorb 0.05; the floor is rigid. The SEL change from specular walls, median
and range over the eight placements of each:

| Microphone | Street against open ground (s = 0) | s = 0.1 | s = 0.2 | s = 0.3 | s = 0.5 |
| --- | --- | --- | --- | --- | --- |
| Street along x | +0.43 [−1.08, +2.98] | +0.14 [−0.67, +0.57] | +0.22 [−1.50, +0.94] | +0.26 [−2.22, +1.19] | +0.31 [−3.34, +1.45] |
| Street along y | +0.30 [−0.61, +2.69] | +0.13 [−0.73, +0.67] | +0.24 [−1.65, +0.99] | +0.30 [−2.44, +1.19] | +0.32 [−3.74, +1.41] |
| 2 m from a wall | +0.63 [−3.76, +3.26] | +0.12 [−0.39, +0.88] | +0.17 [−0.68, +1.38] | +0.20 [−1.34, +1.78] | +0.22 [−2.59, +2.43] |
| Plate at a wall | −0.54 [−1.78, +2.41] | +0.16 [−0.15, +0.87] | +0.26 [−0.85, +1.12] | +0.42 [−1.66, +1.46] | +0.51 [−3.10, +1.98] |
| Intersection | +0.51 [−1.91, +2.20] | +0.17 [−0.47, +0.77] | +0.23 [−0.92, +1.29] | +0.34 [−1.41, +1.69] | +0.34 [−2.34, +2.22] |

The median microphone moves by +0.14 to +0.30 dB between s = 0.1 and 0.5, but single microphones
move as much as the street does: at s = 0.2, 35% of the 160 by more than 0.5 dB and 6% by more than
1 dB. Scattering pulls the street's effect toward nothing. Where the street raised the SEL -- the
deep streets (H/W 2) aside the approach, filled by specular bounces from its low, distant part --
it lowers it, by up to 3.7 dB at s = 0.5; where the street lowered it, the median rise is +1.0 dB
at s = 0.5, the shadow filling. So s is worth knowing for a particular site, and matters most in
deep streets and off the track; averaged over a district it shifts levels by a few tenths.

**Districts of many orientations.** 81 blocks 54 m square on an 80 m grid, each turned at random
within ±2° or ±10°, the microphone in a street; the cost of a direction's rays:

| Blocks turned | Tolerance | Orientations in the district | Candidate departures | Time per direction |
| --- | --- | --- | --- | --- |
| ±2° | 1e-5° | 162 | 4100 (the cap) | 14.5 ms |
| ±2° | 1° (default) | 6 | 62 | 0.8 ms |
| ±10° | 1° | 28 | 3463 | 12.8 ms |
| ±10° | 3° | 12 | 734 | 2.0 ms |

Exact orientations put every district at the cap, where the breadth-first candidates are cut off
arbitrarily; a degree's tolerance returns a near-grid district to a grid's cost.

**GeoJSON.** The demo district written as WGS84 polygons through the demo frame's inverse
projection (`tests/data/urban_buildings.geojson`) and read back with the run's origin gives every
table cell to 2.4e-5 dB of the CSV's, and the run to 1e-7 (`urban_buildings_geojson_is_the_csv`).
Evaluated exactly on its azimuths, the same pair of tables disagrees in 452 of 16632 cells, by up
to 23 dB, every one at an azimuth parallel or diagonal to the walls: there a wave grazes a wall or
passes through a corner exactly, and the projection's 1e-8 m of round-off decides on which side.
Taking each direction 1e-6 rad off its azimuth settles them, and moves the demo microphones' SEL
by at most 0.006 dB.

**Scattering walls** (`src/diffuse_facades`, `tests/test_diffuse_facades.cpp`):

- **A lit wall's diffuse field** at a microphone 200 m in front of a face 4 km by 2 km, in 20 m
  patches, is rho s cos(incidence) times the face's solid angle over pi to 0.2%.
- **The reverse solve is the forward one** to 1e-8, in three bands over twelve directions, for a
  street whose walls exchange energy through the air and off the floor.
- **Every specular reflection scatters its share.** Between two walls 30 m tall and 20 m apart,
  under a wave at 20 degrees, each of 600 patches' irradiance is cos(incidence) times the product
  of rho (1 - s) over the reflections it takes to climb out, to 1e-12, and nothing where the
  last reflection is on the wrong wall.
- **Form factors are reciprocal**, A_p F_pq = A_q F_qp, to 1e-12 where neither is subdivided.

What s = 0.2 on every wall does to the district above (2 kHz, the mean level over azimuth by
elevation):

| Elevation | 0° | 4.5° | 9° | 13.5° | 22.5° | 40° | 70° |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Specular walls | −38.9 dB | −4.1 | −4.5 | −4.6 | −5.6 | −6.4 | −5.8 |
| s = 0.2 | −18.8 dB | −4.8 | −5.0 | −3.8 | −3.0 | −2.2 | −2.1 |

At grazing the diffuse field fills the street's shadow by 20 dB; at low elevations the specular
bounces up the street keep 0.8 of their energy where they kept 1 and the diffuse field nearly
makes it up; higher up the lit facades scatter into the azimuths the street shades. Over every
direction and band the table's energy mean rises from −0.34 to +0.86 dB. Lighting the walls with
the plane wave and its floor reflection alone, the cheaper model, reads 2-3 dB low from 4.5 to
13.5 degrees: the specular bounces lose their scattered share and nothing receives it. On the demo track s = 0.2 moves the seven test microphones by −0.01 to +0.84 dB
against specular walls: most in the street (−2.09 to −1.25 dB against open ground) and the
intersection (−1.50 to −0.97), nothing for the microphone behind a building, which sees only that
building's shaded face. A receiver's diffuse field adds 1 s of setup and about 7 s of lookups on
one thread for the 100-block district.


On the demo track, the seven test microphones among eleven buildings
(`tests/data/urban_buildings.csv`), SEL against open ground: −2.09 dB in a street, +0.01 at an
L-shaped building's corner, −1.50 at an intersection, +1.39 beside a tall wall (+1.35 as a
plate), 0.00 in the open, −2.13 in an alley (−2.06 as a plate) and −7.25 behind a building. No
three-dimensional wave solution or measurement has checked the builder yet: that is stage 4 of
the plan and the Drone DigiCity data.

### A point source among buildings

The point-source model (`src/point_paths`, `src/double_wedge`, `tests/test_point_paths.cpp`) is
checked against the plane-wave builder and against an exact solution.

- **Far away it is the plane wave.** A source 1000 km off in three directions, one building, the
  microphone beside it, with the plane-wave builder's ray set (no reflections or floor between
  source and edge): scaled by the range, every band equals the table's ratio times open ground's
  energy to 2e-3 dB.
- **It is reciprocal**, source and microphone exchanged, to 1e-9 dB in a district of four
  buildings, through every kind of ray including the shared-face form.
- **The shared-face integral is converged.** Filon's rule, as the model takes it, is within
  0.016 dB (0.0015 dB rms) and 1.5e-3 rad of the midpoint rule at a 192nd of a wavelength over
  some 4000 shared faces of DigiCity and a 49-block district; the midpoint rule at a 24th of a
  wavelength is within 0.027 dB (0.014-0.018 dB rms).

**Double diffraction against boundary elements** (`tools/thick_barrier_bem.py`). A rigid block on
rigid ground, long enough to be two-dimensional, the source and microphone in one plane across it:
there every ray is a 2-D ray times 1/sqrt(its length), so the field in that plane is a
line source's, which boundary elements solve exactly. The solver reproduces the series for a rigid
cylinder to 0.005 dB at 800 Hz; at 60 elements a wavelength its values here move by under 0.2 dB
from 30 (a third of that is the error estimate, and points above 0.2 dB are left out). Blocks 1,
2.44 and 5 m thick and 5.18 m high; receivers 2 to 25 m behind, from 1.2 m up to 8 m, below and
above the roof's plane; 125 Hz to 1 kHz. Error re the exact level, range and rms, dB:

| Block | Source | Microphone | Cascade | Shared-face form |
| --- | --- | --- | --- | --- |
| 1 m | 18° below the roof (10 m out, 2 m up) | below the roof | −6.9 to +0.5 (2.3) | −0.7 to +0.3 (0.28) |
| 1 m | | above | −0.2 to +3.7 (1.4) | −0.4 to +1.7 (0.40) |
| 1 m | 1.6° below (25 m out, 4.5 m up) | below | −16.0 to +0.0 (6.0) | −2.9 to +0.6 (0.92) |
| 1 m | | above | −0.8 to +5.0 (2.3) | −0.3 to +1.8 (0.59) |
| 1 m | 5° above (15 m out, 6.5 m up) | below | −0.2 to +4.0 (1.8) | −0.4 to +1.3 (0.50) |
| 2.44 m | 18° below | below | −3.6 to +0.5 (1.1) | −0.2 to +0.4 (0.15) |
| 2.44 m | | above | −0.4 to +2.2 (0.8) | −0.4 to +0.8 (0.18) |
| 2.44 m | 1.6° below | below | −10.9 to +0.3 (3.9) | −0.9 to +0.4 (0.27) |
| 2.44 m | | above | −0.5 to +4.4 (1.8) | −0.2 to +0.8 (0.23) |
| 2.44 m | 5° above | below | −0.2 to +2.5 (1.1) | −0.2 to +0.5 (0.14) |
| 5 m | 18° below | below | −1.9 to +0.4 (0.5) | −0.3 to +0.2 (0.12) |
| 5 m | | above | −0.2 to +1.3 (0.4) | −0.1 to +0.6 (0.13) |
| 5 m | 1.6° below | below | −7.2 to +0.3 (2.5) | −0.4 to +0.2 (0.14) |
| 5 m | | above | −0.3 to +3.6 (1.3) | −0.1 to +0.6 (0.14) |
| 5 m | 5° above | below | −0.2 to +1.4 (0.6) | −0.2 to +0.5 (0.12) |

The cascade is low wherever the second edge sits in the first's transition region, by 11 dB for a
container stack under a grazing source, and the error does not fall away with frequency: at 1 kHz
the grazing case is still 6.7 dB low. With both ends above the roof's plane the cascade stands and
is off by −10 to +4 dB at 125-500 Hz, mostly the reflection on the roof between the two edges,
which a twice-diffracted ray does not take. Thinner
blocks, from the first study (source 18° below, receivers below the roof), mean error of the
shared-face form against the cascade's:

| Block | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz |
| --- | --- | --- | --- | --- | --- |
| 0.05 m | −5.8 (−8.4) | −4.4 (−7.2) | −3.8 (−5.6) | −2.0 (−4.8) | −1.1 (−2.9) |
| 0.2 m | −3.2 (−4.5) | −1.9 (−4.2) | −0.9 (−2.8) | −0.5 (−2.5) | +0.1 (−0.8) |
| 0.5 m | −1.2 (−2.6) | −0.9 (−2.9) | +0.0 (−1.4) | −0.2 (−1.4) | −0.1 (−0.2) |

Within a wavelength or so of each other the edges are in each other's near field, where neither
form holds; a block that thin is a screen, and LoD1 buildings are not.

**Across the roof's plane** (a block 2.44 m thick and 8 m high, no floor, the microphone 1 mm above
and below the plane 1 to 50 m behind), the level steps by at most 0.65 dB for a source 15 m out and
2 m up, 1.9 dB at 125 Hz for one 40 m out and 6 m up (under 0.25 dB from 1 kHz), and 0.43 dB for
a block 10 m thick. Against the boundary elements either side of the plane the shared-face form is
within ±0.6 dB, the exact field moving by 0.1 dB; the cascade steps by 5 to 15 dB. A source above
the roof steps by up to 4.9 dB, where the cascade takes over above the plane.

**Across a step** (THEORY §6.5; a long block 10 m across and 10 m high against one 10 m across
and 7 m high, rigid ground, boundary elements at 45 elements a wavelength, the cases not converged
to 0.15 dB against 30 left out; 125-500 Hz). The sound from above crosses the taller block's far
edge and then the lower one's, which share no face; the lower roof is 16.7° below the line between
them. Error re the exact level, dB:

| Microphones | Sources | Cases | Uniform form: rms, mean, worst | Cascade |
| --- | --- | --- | --- | --- |
| Below the lower roof, 2-10 m behind, 1.2-7 m up | 10.5-15 m up, 30-100 m off | 229 | 0.71, +0.18, +2.42 | 0.96, −0.58, −3.65 |
| The same | 1.5 m up, 20 m off | 30 | 0.88, −0.67, −2.20 | the same |
| Above it, 8-11 m up | 10.5-15 m up | 72 | 2.55, +0.36, +10.76 | 2.53, +0.40, +10.16 |

The first row includes the same step with the taller block 9 m high. From the low source the
pair is not tilted -- the sound crosses the taller block's roof as a shared face first -- and nothing
moves.
Above the lower roof both are as far off, nearly all of it at the lower roof's reflection
boundary for the first edge's wave (8.5 m up at 5 m behind, 10 m at 10), where the ray set steps
by up to 6 dB with either form: the first edge's ray reflected off the lower roof switching on
is not matched by the pair's partner. Across the line through the two edges the uniform form's
level steps by 0.45 to 2.1 dB from 1 mm below to 1 mm above, the cascade's by 0.02 to 1.5. The
worst is a source 0.26° over the taller block's roof: the end correction is 2.1 there, the
first wedge's coefficient taking that roof's reflection as if the roof went on, which the
three edges over it (its near edge too) would put right and two do not. Where the microphone's
image in the step's own wall lies behind the first edge, or the source sees the lower roof's plane
while the microphone takes that roof's reflection, the cascade stands; THEORY has why.

**EDGE's stepped terraces** (model scale; boundary elements at 45 elements a wavelength, converged
to 0.1 dB; the drone and the floor source, microphones on the floor 0.46 and 1.25 m behind,
1.6, 3.15 and 6.3 kHz, 12 cases each). Worst and rms, dB:

| Terrace | Uniform form | Edge source | Cascade |
| --- | --- | --- | --- |
| UE01, half gable: ridge, a 0.24 m step down to a flat roof | −0.42, 0.23 | −0.42, 0.23 | −1.63, 0.66 |
| UE02, pent: ridge, a 0.12 m step down to a deck, its edge | −2.75, 1.33 | −2.75, 1.33 | −4.99, 2.68 |
| UE03, gable: a shared face | −0.30, 0.13 | +0.20, 0.07 | the same |

UE02's drone row at 1.25 m goes from 4.0-5.0 dB low to 0.15-0.53; its rows at 0.46 m are left,
2.75 dB at worst: there the microphone is 6 mm inside the far eave's shadow from the deck's edge,
three edges whose pairs are each in transition. The edge source's own integral was 1-3 dB low across
these steps, and leaves them to the uniform form (THEORY §6.5).

**Across a street at grazing** (two blocks 10 m across and 10 m high, 10 m apart, a source 1.1°
above the roofs 100 m off, microphones 2-10 m behind the second at 1.2 m; 125-500 Hz). The model is
1.2-1.6 dB low at 125 and 250 Hz and 7.9-8.0 dB high at 500 Hz, the same at every microphone, so
the error is in what reaches the second block's far edge; 5.7° above the roofs it is within 1.3 dB,
and with the second block 2 m lower, 1.2-3.2 dB high. The four roof edges in a row are each in the
others' transition regions, and the cascades of two, three and more edges do not converge: leaving
out the thrice-diffracted rays moves the mean error from +1.1 to −3.8 dB. No uniform form for
more than two edges is in the model. Nor is a street a step: taken as one (the near block's far
edge and the lower block's far one, the street read as roof), the second block 2 m lower read
0.5 dB higher still, rms 2.0 dB against the cascade's 1.9, and the cascade stands there.

**On the GPU** (the default on a Mac, an M4 Pro), which works in single precision:

- **The sums.** Against double-precision sums of the same amplitudes, over ray sets of 200 to 1100
  rays with lengths from 120 m to 2 km and bands from 50 Hz to 10 kHz, the GPU's pair sums agree
  with the CPU's to 1.1e-5 dB, and its Gauss–Legendre integrals with every pair summed to
  3.4e-6 dB.
- **The amplitudes.** The UTD transition function agrees with the CPU's to 6e-7 and the wedge
  coefficient to 3e-5 (relative, wedges of n = 1 to 2 including angles within a hair of a
  boundary), and a Filon panel with its exact integral to 4e-5 of the panel. The shared face's
  aperture is worked out once for all of a run's wavenumbers rather than for each, its panels
  ending at each one's taper and extent: on the Filon test's faces it is within 0.006 dB of the
  fine midpoint rule from 50 Hz to 10 kHz, where the CPU's own is within 0.017 dB. Band by band,
  amplitudes the GPU forms give levels within 1.6e-4 dB of the same batch given the CPU's (five
  pairs among containers and blocks, 96 rays over shared faces, 50 Hz to 10 kHz), and with the
  shared faces sampled on both, within 2.4e-4 dB.
- **The partial sums.** The cone's and the energy sum's band levels, over the same five pairs and
  nine bands from 50 Hz to 10 kHz, with a cone of 1° and of 6°: within 3e-4 dB of
  `partially_coherent_energies`.
- **In a run.** `niceops --buildings` on the demo approach: SELs within 1e-6 dB of the CPU's, with
  the states sampled and with every state evaluated. `niceops_point_field` over every fifth of its
  source positions, the shared faces sampled on both: band levels within 0.04 dB.
- **The levels** differ from the CPU's by what the CPU's cutoff leaves out, the pairs whose sinc
  has fallen below 1/30, which the integral keeps:

| Case | Band levels | Median | 99th percentile | Worst |
| --- | --- | --- | --- | --- |
| 49-block district, a drone 30 m up, 20 microphones, 3 bands | 60 | 4e-5 dB | 0.012 dB | 0.012 dB |
| DigiCity, the loudspeaker at hover spot 5, 100 grid points, 3 bands | 300 | 5e-5 dB | 0.015 dB | 0.030 dB |
| 49-block district, 60 points of a drone's track 40 m up, 185 microphones, 24 bands | 266,400 | 2e-5 dB | 0.022 dB | 0.090 dB |

**Continuity.** Where every boundary's partner is traced the field is continuous: west of a tall
box, lines of microphones 2.5 cm apart across its corner's shadow and reflection boundaries step
by 0.32 dB at most between neighbours (125 Hz to 2 kHz). Across a city it is not. A sweep of
DigiCity (80 random lines of microphones 2.5 cm apart, 1.2 m up, the source at hover spot 5 and at
four points of flights 5 and 8.5 m up, 79,460 positions, 160 Hz to 5 kHz summed flat, walls
absorbing 0.05, 70% RH) finds 434 steps of more than 1 dB between neighbours, 0.55% of positions,
1.7 dB at the median and 9.5 dB at most. By the ray that appears or vanishes there, and why:

| Class | Steps |
| --- | --- |
| A diffracted ray's leg meets another edge, or its reflection leaves its wall: the ray that would take over needs a second edge with those reflections, or a third | 265 |
| A specular ray's first or last leg meets an edge: the diffracted ray that would take over needs 3 or 4 reflections on one side, beyond `--max_edge_reflections` 2 | 90 |
| A diffraction point leaves the end of its edge | 23 |
| A twice-diffracted ray's point reaches the floor at a wall's foot: it would carry on with the floor between its edges | 13 |
| The same rays, summed with a sign flipped | 15 |
| None of these found (the partner traced: 7) | 28 |

With `--max_edge_reflections 4` the second class goes but the diffracted rays it adds have their
own: 490 steps, 10.1 dB at most, for 20% more time.

That sweep predates three roof edges, a reflection between two edges (a facade or the floor) and two
reflections at each end of a twice-diffracted ray. On 20 of its lines (3771 positions, the source at
hover spot 5 only) those took the steps of more than 1 dB from 31 to 23, the median step from 1.8 to
1.5 dB and the largest from 9.5 to 4.9; the floor between two edges alone removed the two at a
wall's foot, for 8% of the demo `--buildings` run.

The cutoff itself, against every pair summed, moves the district's band levels by 0.002 dB at the
median and 0.037 dB at worst (40 microphones, the drone 30 m up, 24 bands). The ray dumps are
identical either way. The table's runs leave the air's absorption out; with it (70% RH, the
district's first row), the GPU's levels are within 0.012 dB of the CPU's, as without, and the
GPU's amplitudes within 2.4e-4 dB of the CPU's over the amplitude test's five pairs.

Cost: in the table builder's test district, a point source 10 m up in a street, the seven test
microphones take 0.44 s on one thread over B407's 31 bands. Among 49 blocks with roofs 10 to
40 m, a drone 30 m up costs 0.042 s per microphone on one thread over 24 bands. 60 points of a
drone's track 40 m up, over 185 microphones, take 30 s on 14 threads with `--sources` (34 s as a
run per source), and 3.9 s on the GPU. These were measured with every shared face's integral at
every wavenumber (`--exact_shared_faces`).

**An edge's ends in UTD** (THEORY §6.5, "An edge's ends"; `tests/test_point_paths.cpp`,
`test_wedge.cpp`, `test_buildings.cpp`). Against the edge-source integral, which takes each edge's
ends exactly: a box 10 by 20 m and 6 m tall, the source 3 m up 15 m before it, the microphone
swept up past its corner's top (15, −5, 1-14 m) and across past its end (15, −8 to 8, 9 m) in
0.25 m steps, rigid. Away from the planes of the corner's faces through the vertex (0.5 m or
more), rms and worst per band (dB):

| Band | Up, with ends | Up, unending | Across, with ends | Across, unending |
| --- | --- | --- | --- | --- |
| 125 Hz | 0.34, 0.72 | 3.08, 7.13 | 0.24, 0.59 | 2.05, 7.10 |
| 250 Hz | 0.28, 0.60 | 0.66, 1.39 | 0.22, 0.69 | 1.57, 4.01 |
| 500 Hz | 0.20, 0.86 | 0.88, 2.58 | 0.08, 0.30 | 1.05, 2.42 |
| 1 kHz | 0.05, 0.20 | 0.78, 2.39 | 0.06, 0.28 | 0.76, 2.39 |
| 2 kHz | 0.03, 0.08 | 0.70, 2.98 | 0.03, 0.11 | 0.73, 2.98 |

Within 0.5 m of those planes -- the double transition, where an end's transition meets a face's
reflection boundary and the per-term weights stand in for the vertex's own coefficient -- the ends
are no better than none: 1.21 dB rms and 3.68 at worst (125 Hz, on the roof's plane) going up,
against 1.03 and 3.13 unending; 0.82 and 1.96 going across, against 1.01 and 2.66. Albani et
al.'s vertex coefficient (IEEE TAP 57, 2009), re-derived for rigid faces, is what would put that
right.

Continuity: as the corner's Keller point leaves its top, and past the end as the microphone
crosses the plane of a face into the wedge's solid side, a micrometre's move changes the field by
1e-5 dB, where with the edges unending the first is a 3.8 dB step. Where a third edge cuts a
twice-diffracted ray the handover steps 0.063 dB at worst (125 Hz), against 0.046 with no ends
anywhere. Image theory (a low block's roof edges running into a wall filling a half-space) holds to
0.02 dB by UTD and by the edge source, the edges going on in their images at the wall. A microphone flush on a
wall facing the source is within 0.107 dB of one a hair off it (0.007 with the edges unending): the
backscatter off the far corner along the face, which only the one off the face has, is stronger
with the corners' ends.

On DigiCity's 3-D render (a loudspeaker flown along the central street 8.5 m up, 177 positions,
11,250 cells of floor, walls and roofs 0.5 m across, 160 Hz to 5 kHz summed), 0.17% of cell-frames
sit 3 dB or more off the median of their neighbours on the same surface, and 0.010% dip or peak
3 dB against both the frames either side; 0.17% of roof cells ever do so in the flight. With the
edges unending (and the model of the day before step 1a, which moved microphone 7 by under 0.1 dB)
those were 0.33%, 0.025% and 14.4%: the one-frame dips where two roof edges meet at a corner.

**The edge-source integral** (`--edge_source`; `src/edge_source`, `tests/test_edge_source.cpp`).
First order, against exact solutions:

- **The plane-wave limit** (a point source receding along the incidence direction, extrapolated from
  1e8 and 2e8 m): Sommerfeld's rigid half-plane within 8e-7 of the field (k r from 0.5 to 200, four
  incidences, every 7.5°, boundaries included); the right-angle wedge's eigenfunction series
  (n = 3/2, k r 10 and 40) within 1.5e-6 dB.
- **Against UTD** over 300 random wedges, a third within a few degrees of a boundary: 0.16, 0.021
  and 0.0005 dB at k = 5, 50 and 500 m⁻¹. Reciprocal to 1e-13. Across a boundary and as the apex crosses an edge's end, steps of
  2e-9 rad or m move the field by 8e-8 and 4e-8, its own slope.
- **The quadrature**: the steepest-descent tails within 1.4e-8 of direct integration along the edge;
  the 72-wavenumber spectrum (Filon's rule in the path) within 4.4e-5, 0.0004 dB, of the field at
  each.

Second order:

- **The level curves are the nested integral** to 4e-6 (two screens, a block's roof with the plane
  factor 1/2, and a screen's top with a perpendicular edge whose boundaries cross it; k = 1 to 30).
- **Over the long block** of the table above (2.44 m, all 80 cases), through the whole model, the
  field is within 0.28 dB of the boundary elements, 0.10 dB rms, where the shared-face form is
  within 0.94 dB, 0.21 rms. Without the plane factor of 1/2 over the roof the integral is the
  cascade without UTD's grazing half, to 0.01-0.1 dB at k = 20.
- **The cascade is the integral away from transitions**: two screens with the second edge at 0, 45
  and 90° to the first, k = 50 and 100, within 0.007 dB and 0.005 rad.
- **The partners' corrections are what make it localized.** On the perpendicular pair, in a window
  of 40, 80 and 160 radians about the least path, the second order alone gives −68.94, −69.84 and
  −70.44 dB (k = 20) as the boundary point enters; with the corrections, −68.944, −68.945 and
  −68.948, the corrected cascade −68.963.

**What it changes in DigiCity** (11 microphone positions across the city, the loudspeaker at hover
spot 5, 250 Hz and 2 kHz; 1241 once- and 4825 twice-diffracted rays). Once-diffracted, UTD against the
integral along the edge's infinite line: 0.001 dB at the median, 0.01 at the 90th percentile,
transitions included. Against the integral over the edge's own length, the ends counting (a
vertical edge on the floor taken through it): 0.25 dB at the median and 1.6 at the 90th percentile
at 250 Hz, 0.10 and 0.69 at 2 kHz; by the excess of the path through the nearer end, in radians:

| End's excess | Rays | Median | 90th percentile |
| --- | --- | --- | --- |
| under 1 | 179 | 3.0 dB | 5.5 dB |
| 1-3 | 191 | 1.0 | 1.7 |
| 3-10 | 547 | 0.29 | 0.99 |
| 10-30 | 642 | 0.14 | 0.58 |
| 30-100 | 592 | 0.10 | 0.33 |
| over 100 | 331 | 0.03 | 0.15 |

Twice-diffracted, over parallel edges, the cascade against the integral over infinite lines: 0.01
to 0.03 dB at the median and at most 0.06 at the 90th percentile wherever the smaller transition
parameter is 0.3 or more; below 0.3, 0.025 at the median and 3.3 dB at the 90th percentile, 15% of
rays over 1 dB. Over the edges' own lengths (in a window of 40 radians) the median ray moves by
0.89 dB at 250 Hz and 0.06 dB at 2 kHz. DigiCity's four microphones, the loudspeaker at hover
spot 5, 160 Hz to 6.3 kHz:

| Microphone | Up to 500 Hz | From 2 kHz | Most |
| --- | --- | --- | --- |
| 2 | −0.01 to +0.03 | 0.00 | +0.03 dB at 315 Hz |
| 7 | −0.50 to +0.64 | −0.16 to +0.28 | +0.64 dB at 250 Hz |
| 8 | −0.32 to +0.10 | 0.00 to +0.03 | −0.32 dB at 400 Hz |
| 10 | −0.45 to 0.00 | −0.01 to +0.01 | −0.45 dB at 200 Hz |

**Against measured scale models** (RWTH's EDGE database: 12 LoD1 scenes at 1/5 to 1/20 scale in a
hemi-anechoic room, rigid floor, receivers on the floor mostly in the buildings' shadow; each ray
given its loudspeaker's measured directivity at its emission angle and the air's absorption,
summed coherently over the band; pairs at least 15 dB over noise whose repeats agree to 0.5 dB,
1.6-20 kHz at model scale; one gain per loudspeaker fitted, the archive's being undocumented),
the two are close: 4.01 and 4.04 dB rms over 1124 band levels, 4.25 and 4.29 over the 845 in
shadow. At model scale the buildings are many wavelengths long, so their ends matter little, and
what is left is not diffraction's form. Read with the archive's UE06 and UE07 measurements swapped,
as the full-wave runs suggest (below; not confirmed), the two models are 3.74 and 3.72 dB rms over
all, 3.91 and 3.89 in shadow. The floor-source rows in UE06 and UE07 are 2.88 and 2.66 dB rms by
UTD and 2.12 and 2.53 by the edge source, against the FDTD's 2.30 and 2.32.

**Continuity.** On 20 of the sweep's lines (3771 positions, the loudspeaker at hover spot 5 only),
`--edge_source` left 30 steps of more than 1 dB where the model without it left 31 (9.1 and
9.5 dB at most), before either took a third edge or a reflection between two: the steps are
partners the orders leave out -- a leg meeting a third edge, a reflection past
`--max_edge_reflections` -- which neither model traces.

**Cost.** DigiCity's four microphones, 17 bands, on one thread: 0.22 s without `--edge_source`,
58 s with it when first written. Since (14 threads, against the integral as first written):

| Case | UTD | Edge source, before | Now | Now against before, A-weighted |
| --- | --- | --- | --- | --- |
| DigiCity, four microphones, 160 Hz-6.3 kHz | 1.3 s | 145 s (357 s CPU) | 16.9 s (155 s CPU) | within 0.0001 dB (0.001 in a band) |
| Demo district, 38 states × 7 microphones, 10 Hz-10 kHz | 0.27 s | 497 s (5177 s CPU) | 20 s (223 s CPU) | within 0.012 dB (266 pairs; 0.06 in a band) |
| EDGE, the 13 scene-source runs against the FDTD | | minutes | 3.6 s | each run's rms within 0.04 dB |

The integral's tolerances a hundred times looser, the cascade times the edges' end ratios out of
both transitions, and a pair's integrals over threads (THEORY §6.5; ground.md). The ratios once
held the other edge's point fixed and were chosen at the run's least wavenumber: 10.5 s on DigiCity
then, but 0.09 dB off in a band there and 0.75 on the demo district, up to 3 dB with one end near
on coupled edges, and a band moved with the run's band limits. Each wavenumber now decides alone,
the other edge's point slides, and a ray with both edges' ends near takes the integral, which
DigiCity's small containers pay for. What is left is the
twice-diffracted rays near a transition at low frequency and over shared faces, about 0.8 s of CPU a
pair on the demo district where UTD takes 4 ms: a trajectory's few thousand pairs in a run
(`niceops --buildings --edge_source`) take minutes, not seconds, so UTD stays the default.

**Where it departs from UTD most, unchecked.** On the demo district 6 of the 266 pairs move by more
than 1 dB A-weighted, all at one microphone 3500 ft from the source deep behind the tallest
building, where the levels are 40 dB below the district's: up to +20 dB. Two things carry it there,
once-diffracted rays kept past a roof edge's end (the corner's wave, which UTD has no ray for) and
a twice-diffracted ray over a roof edge and the vertical edge it meets, 29 dB above its cascade near
their shared vertex. Neither has been checked against a full-wave field: the case is a kilometre
across. The long block's 40 cases in its test: 4.9 s, against 0.01 s for the uniform form's 80.

Taken over [0, L], a vertical edge on the floor has a false end there, a metre or so from every
street-level ray, which moved microphone 7 by −1.5 dB at 160 Hz. Leaving the integral only to rays
within 45 dB of the band's energy moved microphone 7 by up to 0.48 dB: hundreds of rays each under
that carry half a decibel between them.

**The plane-wave table under a hover.** DigiCity (six container stacks 5.2-7.8 m tall,
`urban/digicity` in the harness), its 22 microphones at 1.2 m, a hover at each of its seven hover
spots at 15, 30 and 60 m, the B407's hover sphere as a directional source (there is no drone
database), rigid floor: the plane-wave table (`--urban_table`, on ground of 10^6 kPa s/m^2) differs
from the point model, its rays summed coherently, by 0.42 dB rms at 60 m, 0.89 at 30 m and 2.34 at
15 m. The table's distant-source approximation recovers only as the source rises well clear of the
buildings.

**Roofs as reflectors** (`niceops_point_field`, `--no_roofs` to compare). A microphone 1.5 m above
a 200 m roof, its edges 90 m off, the source 30 m above the roof: the direct ray and the roof's are
two rays over a rigid floor to within 0.05 dB in every band, and without roofs the level is 1.5 dB
or more lower. Where a roof's edge hides its reflection the field is continuous to 0.01 dB from
1 µm either side of the boundary. On EDGE's scenes (model − FDTD in [the full-wave reference](#the-full-wave-reference-among-buildings)), roofs move only the
drone over rows of unequal houses, where the low roof reflects up to its taller neighbours' edges:
UE07's drone from 2.54 to 2.29 dB rms with UTD and from 0.95 to 0.82 with the edge source, its
worst band from 2.9 to 2.1 dB; UE06's from 1.98 to 2.00 and from 1.43 to 1.47. The floor
sources, and the drone over the even rows, do not move: a ray off a roof reaches only what stands
above that roof, or an edge higher still. UE06's floor source, 5 dB low where the FDTD's
strongest arrival crosses a roof between two edges, needs a reflection between edges, which the
model does not have.

### Coherence between emission directions

How `--buildings` adds rays that leave the source apart (THEORY §6.5) was fitted to flight tests in
the harness's `coherence/` (its README has the method and every table): the 2017 Bell 407 level
passes 283104-05 (60 kt) and 283108-09 (100 kt), 47-48 ground boards within 1500 ft, its high
hovers H7-H10 (runs 283427-30 and 283437-40, 13 boards on a 500 ft semicircle), and PANAM's
hexacopter flyover (13 boards). Each board's pressure is de-Dopplerized onto emission time, and per
pair of boards, window and third-octave band the band coherency is the cross power over the root of
the two auto powers, segment lengths following the band so the geometry is taken while it holds.
A model is scored by the mean of |ρ − model|², which sets the mean square error of two rays' band
energy at an unknown relative path phase:

| | coherent (ρ = 1) | energy (ρ = 0) | the rotor's delay alone | the model | each band and separation's own mean |
| --- | --- | --- | --- | --- | --- |
| B407 passes | 1.156 | 0.294 | 0.439 | 0.278 | 0.276 |
| B407 hovers | 1.130 | 0.247 | 0.423 | 0.227 | 0.212 |
| Hexacopter, its rotors' delays half each way for the delay | 1.039 | 0.214 | 0.358 | 0.210 | 0.191 |

The model here is the B407 fit, applied to the hexacopter without its delay. Its gain is where a
ray sum feels it: B407 passes at 25-63 Hz and 10-20°, 1.26 coherent, 0.86 energy, 0.60 the model.
What the fit rests on:

- **The decorrelation is the source's.** At a fixed angle the passes' coherence does not fall with
  range from 100 to 500 m, though the boards' spacing grows fivefold: at 500 Hz and 5-10°, 0.37,
  0.47, 0.45, 0.37 and 0.41 from 0-100 m out to 300-500 m.
- **The rotor's delay holds to about 4.5 blade-passage frequencies.** Scaled by κ, the best κ is 1
  in hover and 0.75-1 in the passes in the main rotor's bands to 125 Hz, and 0 from 160 Hz, where it
  turns in-phase pairs out of phase (passes, 250 Hz, 1-5°: +0.64 without it, −0.16 with it). Applied
  at every frequency it scores worse than none (0.439 against 0.284 on the passes); where it stops,
  between 2.5 and 7 blade-passage frequencies, moves the score by 0.0005. In hover, 70 % of 1219
  pair-windows' fitted delays lie within 3 ms of δψ/Ω (chance: 17 %).
- **The decay.** g = exp(−(fθ/θ₀)^p): θ₀ = 850 Hz deg, p = 1.61 with the delay (passes 725/1.84,
  hovers 940/1.75, the joint fit weighting the two equally); 475/2.33 without it, fitted to all
  three, which scores 0.2836 on the passes, 0.2457 on the hovers and 0.2132 on the hexacopter where
  the with-delay pair would score 0.2866, 0.2588 (worse than energy) and 0.2104.
- **Not the hexacopter's layout.** Its rotors as mutually incoherent point sources on a ring of
  radius R, each coherent between directions, predict near-full coherence at 50-250 Hz and 10-40°
  for R = 0.5-1 m, where the measured mean is 0 to +0.5; the best R runs to the bound (20 m) and
  scores 0.221, worse than energy.
- **The spread** σ² = γ² − g², γ² = T + (1 − T) exp(−fθ/2120 Hz deg), T each band's coherence
  beyond 40° (the B407's: 0.89 at 25 Hz, 0.74 at 50, 0.39 at 80, 0.09 at 125, 0.14 at 250, 0.01 at
  500), reproduces the measured spread to 0.09-0.12 rms over bins of band and separation on all
  three cases. A uniformly random phase (σ² = 1) misses it by 0.76-0.82, and T(1 − g²) alone by 0.22
  on the passes: at 500 Hz within 10° the coherence is 0.4-0.6 where the mean is 0.2-0.3.
- **Limits.** Only the passes resolve separations below 10° (to under 2.5°); one helicopter and one
  multirotor; the hexacopter's T stays 0.23-0.40 from 125 to 315 Hz where the default law gives
  0.05-0.20. The measured coherence carries the air over the paths, so it bounds the source's own
  from below.

PANAM's `directional_coherence.fit_coherence`, given the passes and the hovers as two cases,
reproduces these parameters to the digits quoted. Written into a copy of the B407 database, its
record gives the DigiCity hovers below within 0.011 dB of the measured model with the B407's rotor
and the fits quoted.

**Why the default is energy beyond a cone.** Scored as above, but each pair-window weighted by its
A-weighted band power (for the B407 that power lies between 160 Hz and 2 kHz), over pairs under
20°: the energy sum scores 0.379 on the passes and 0.081 on the hovers, a cone of any half-angle
from 0.5° to 2.5° the same (0.378-0.379; 0.081), a 5° cone 0.421, and the measured decay 0.370
(0.081). Nothing simple adds A-weighted levels measurably better than the energy sum, which
assumes nothing of the source's phase. No fixed angle suits every band: unweighted, under 10°, a 5°
cone serves 25-63 Hz (0.61 against energy's 0.93) and costs at 800-2000 Hz (0.43 against 0.14).
The cone is not fitted. It keeps together what must be: rays that leave the source together are
one emission, below any separation the arrays resolve, and the point model's continuity at shadow
boundaries and a flush microphone's doubling rest on adding them with their phases.

**What it does to levels.** The DigiCity hovers above (154 microphone-hovers per height). The
default against the coherent, energy and measured sums (the last without a rotor, as a database
that records none takes it), and its 90 % interval:

| Hover height | cone − coherent, rms (max) dB | cone − energy, rms (max) | cone − measured, rms (max) | interval width, median (p90) |
| --- | --- | --- | --- | --- |
| 15 m | 0.77 (2.46) | 1.24 (7.08) | 0.49 (3.60) | 0.74 (1.34) |
| 30 m | 0.44 (2.10) | 0.73 (3.89) | 0.29 (1.61) | 0.62 (1.23) |
| 60 m | 0.22 (1.43) | 0.56 (3.11) | 0.15 (0.68) | 0.10 (0.66) |

The energy sum's 7 dB is the pairs the cone keeps together. The measured model with the B407's
rotor, as a database recording it would give, against the coherent sum: 0.69, 0.42 and 0.23 dB rms
at 15, 30 and 60 m (2.11 at worst), its interval 0.84, 0.58 and 0.35 dB wide (median). For
A-weighted levels the rotor's delay is a refinement: leaving it out moves the hovers by 0.07-0.21 dB
rms (0.64 at worst), mostly through the other decay, and reversing its sense by 0.05 at most; it
acts only below 124 Hz.

On the demo approach over `tests/data/urban_buildings.csv` the default's pass SEL is within 0.08 dB
of the coherent sum at all seven microphones, its 90 % interval 0.07-1.43 dB wide, where the energy
sum loses the floor's doubling at the flush microphones (up to 3 dB). The cone costs what the
coherent sum does (0.200 s against 0.202 s for the rays from three states to 22 microphones on 14
threads, at 15 m), the interval 2.7 times that; the measured sum 1.4 times.

### The full-wave reference among buildings

`niceops_fdtd_field` (`src/fdtd`, `tests/test_fdtd.cpp`) against exact solutions, at 1 cm cells
(25 cells a wavelength at the top band's upper edge, 500 Hz to 1.25 kHz), band levels, worst
error, second order and fourth:

| Case | Second order | Fourth order |
| --- | --- | --- |
| Two rays over a rigid floor: four receivers, near, far, on the floor and across the plane of symmetry | 0.080 dB | 0.0005 dB |
| The same with the plane of symmetry gridded as one, against the whole space | within 0.02 | within 0.02 |
| The air applied by delay (0.4 Np/m), against each ray absorbed along its length | 0.054 | 0.0018 |
| A right-angled wedge with no floor, against the edge-source integral over the whole line plus geometrical optics: three receivers deep in the shadow, one lit, one where the wall reflects | 0.077 | 0.031 |
| A concave right-angled corner, a floor and a wall through the whole domain, against the source and its three images | 0.035 | 0.0002 |
| A wall absorbing 0.36 at normal incidence, against its plane-wave reflection (ζ = 9), on its normal and off it | 0.094 | 0.073 |
| The same wall's absorption alone, from its reflected pulse gated from the direct one (1-3.15 kHz, 5 mm cells), at α = 0.05, 0.2 and 0.36 | 0.0497-0.0499, 0.1987-0.1996, 0.3579-0.3594 | 0.0495-0.0499, 0.1979-0.1994, 0.3567-0.3591 |
| The GPU's march against the CPU's (absorbing walls, every PML slab) | 1e-4 of the peak pressure; band levels within 1e-3 dB | the same |
| Marching only the light cones' boxes, against every cell (a block and a wall along a 1.8 m strip) | 3e-7 of each receiver's peak | 5e-7 |

Second order's bias is the source's and the receivers' trilinear spreading over their cells, low
by an amount that grows with frequency; fourth order's tricubic stencils remove it. The absorbing
wall's impedance condition is second order in either. Between EDGE's canyon walls (UE23, α = 0.05) the fourth
order at 6 mm agrees with the second at 4 mm to 0.03 dB mean, 0.28 rms, and the second at 6 mm
with it to 0.24 rms.

**Over many wavelengths** the scheme's phase error decides. Two rays over a floor out to 2.3 m,
2-5 kHz (13 to 33 wavelengths), the cell set by the 5 kHz band's upper edge; band levels' worst
error and rms over four receivers:

| Cells a wavelength | Second order | Fourth order | Fourth order, warped |
| --- | --- | --- | --- |
| 6 | 9.5 dB, 2.57 rms | 1.55, 0.43 | 1.14, 0.31 |
| 8 | 4.7, 1.30 | 0.58, 0.17 | 0.35, 0.10 |
| 12 | 1.74, 0.52 | 0.18, 0.05 | 0.07, 0.02 |

Warping second order makes it slightly worse (0.57 dB rms at 12): its error depends on direction.
A face off the grid is moved to the nearest cell face: the same corner with both planes 4 mm off,
at 8 mm cells, is out by up to 0.56 dB at 500-800 Hz.

**On a scene.** EDGE's UE00 (one row, four receivers on the floor behind it, 1.6-6.3 kHz, 28 band
levels), on the M4 Pro's GPU. For the drone, against the Richardson extrapolation of second order's
3 and 4 mm runs; for the floor source, against fourth order at 4 mm:

| Run | Drone: worst, mean, rms (dB) | Floor source | Time, floor and drone |
| --- | --- | --- | --- |
| Second order, 3 mm | 0.21, −0.04, 0.07 | 0.11, −0.02, 0.05 | 505 s, 1453 s |
| Second order, 4 mm | 0.38, −0.07, 0.12 | | 467 s (drone) |
| Fourth order, 4 mm | 0.08, +0.00, 0.02 | | 221 s, 681 s |
| Fourth order, 6 mm | 0.08, +0.02, 0.03 | 0.19, +0.05, 0.07 | 40 s, 116 s |
| The same, marching only the light cones' boxes | 0.08, +0.02, 0.03 | 0.19, +0.05, 0.07 | 28.5 s, 81 s |

Second order converges from below, worst at 6.3 kHz behind the row; fourth order at 6 mm, 8 cells
a wavelength at the top band's upper edge, is the nearer of the two to the limit at an eighteenth
of the time. Culling moves no band by more than 5e-6 dB. The tables that follow are second-order runs
at 3 or 4 mm; at 3 mm that is within 0.21 dB of the limit.

**The point-source model against it on EDGE's scale models** (omnidirectional source, no air,
the model's default orders; second order, 1.6-8 kHz at 3 mm cells, 1.6-6.3 kHz at 4 mm; four
receivers each, on the floor). Model − FDTD over every band level, mean and rms (range), dB:

| Scene, source | Run | UTD | Edge source | What it shows |
| --- | --- | --- | --- | --- |
| UE00 (one row 0.49 m high), floor | rigid, 3 mm | −0.25, 0.36 (−1.0 to 0.0) | −0.20, 0.34 | Over one long row the model is the full-wave field |
| UE00, drone | | +0.02, 0.06 (−0.1 to +0.1) | +0.01, 0.05 | |
| UE23 (deep canyon), floor | walls 0.05, 4 mm, 90 ms; line of sight | −0.25, 0.38 (−0.8 to +0.4) | −0.25, 0.38 | With 12 reflections −0.3 to +0.6: a canyon's field is right with orders enough |
| UE06 (five houses, the middle tall), floor | rigid, 3 mm | −0.02, 1.33 (−2.8 to +2.6) | −0.67, 1.00 (−1.6 to +2.3) | The low roofs' notches reflect between their edges: R2-R4 within 1.3 dB rms; behind the tall house −1.1 dB by UTD, −0.1 with the edge source (below) |
| UE06, drone | | +0.48, 2.07 (−3.9 to +5.1) | −0.10, 1.21 | |
| UE07 (the middle low), floor | | −0.21, 2.91 (−8.6 to +4.2) | −0.48, 2.55 (−7.5 to +4.9) | |
| UE07, drone | | −0.19, 2.27 (−9.1 to +3.3) | −0.02, 0.67 (−1.5 to +1.7) | A 1 m roof edge between taller houses: its ends count; the low roof reflects |
| UE13 (two rows, the first higher), floor | walls 0.05, 4 mm, 90 ms | −0.82, 2.73 (−7.6 to +3.1) | −0.79, 2.72 | Behind the second row, over three roof edges: +1.96, 2.20 |
| | | between the rows: −1.74, 2.97 (−7.4 to +0.8) | −1.68, 2.94 | The reflections between the rows after the roof, two each side; the floor's null at 5 kHz |
| UE12 (two rows of one height), floor | walls 0.05, 4 mm, 90 ms | +0.46, 2.30 (−3.4 to +4.7) | +0.47, 2.31 | Behind the second row 2.2-4.7 dB high: two roofs of one height taken for one |
| UE14 (one house 0.96 by 1 m), floor | rigid, 3 mm | +0.09, 0.84 (−2.1 to +1.7) | +0.06, 0.74 (−2.0 to +1.7) | A house a few wavelengths long: its edges' ends count, which UTD takes |
| UE14, drone | | +0.10, 0.45 (−1.0 to +1.2) | −0.01, 0.22 (−0.8 to +0.4) | From above too, the house's ends count |
| UE12, drone | | +0.22, 0.40 (−0.3 to +1.0) | +0.23, 0.40 | From above, one edge to each receiver: the model is the field |
| UE13, drone | | +0.28, 0.36 (−0.2 to +0.9) | +0.27, 0.36 | |

Behind UE12's second row the model's strongest ray diffracts at the first row's front edge and the
second row's back edge, over a "shared face": the two roofs lie in one plane, and the model takes
them for one roof, the street between them included, so the two edges at the street lose nothing.
The model tests only that its leg grazes the second edge's face, not that the face runs on between
the edges. Rows of one height across a street are common; of the twelve shared-face rays without
reflections at DigiCity's four microphones (the loudspeaker at hover spot 5), none crosses a gap.

Behind UE13's second row the model was 8-12 dB low at 1.6-3.15 kHz until it took three edges: the path
crosses the first row's two roof edges and the second's far one, its near one lit below the line.
With the third edge it is 1.2-2.2 dB high there. Between UE13's rows the model is within 1.7 dB at
1.6-3.15 kHz and 5.9-7.4 dB short at 5 kHz.
What reaches those receivers crosses the first roof and reflects between the two rows' facades.
With one reflection each side, as the twice-diffracted rays took until `--max_twice_reflections`,
it was 0.4-3.1 dB short at 1.6-3.15 kHz; three reflections gain 0.3 dB more. At 5 kHz the
rays all leave the floor source at nearly one angle, and the null of its floor reflection,
c / (4 h sin θ) ≈ 5.1 kHz at 29 mm, stays deep where the full-wave field's many paths fill it. Rigid, the same canyon had not settled in 36 ms:
1.6-2.1 dB still moved over the run's last tenth.

**Two of EDGE's scenes may be swapped.** UE06 (a row of five 1 m houses, the middle one tall) and
UE07 (the middle one low): the archive's drawings, photographs and geometry agree with each other
scene by scene, but each scene's measurements fit the other's geometry. Against UE00 (the uniform
low row, same source and receivers), the FDTD's change for UE06's geometry is the measured change
for UE07, and the reverse, within about 1 dB at all four receivers (floor source, 1.6-3.15 kHz,
8 mm cells); behind the tall middle house, for instance, +10.4, +11.4, +10.6 and +8.1 dB from
the FDTD against +10.0, +9.9, +10.1 and +7.4 measured in "UE07", where "UE06" measured −0.6 to
+6.9. This is an inference from the full-wave runs, not a confirmation: the archive's authors have
not said so, the models measured may not be the ones documented, and a different fault in the
experiment could give the same pattern. The comparisons here read the two swapped unless marked,
and the table below gives both readings.

**EDGE against the FDTD** (each ray reads the loudspeaker's directivity along the direction it
leaves the source, each FDTD run given the ray model's correction for the loudspeaker's
directivity and the air, the model's coherent ray sum with them less without; one calibration per
loudspeaker and band from the line-of-sight pairs, as above). Measured − predicted:

| Scene, source | Band levels | FDTD | UTD | Edge source |
| --- | --- | --- | --- | --- |
| UE00, floor | 32, shadow | +0.05 dB mean, 1.78 rms | +0.30, 1.96 | +0.26, 1.96 |
| UE00, drone | 32, shadow | −0.07, 0.62 | −0.09, 0.62 | −0.10, 0.62 |
| UE23, floor | 28, line of sight | +0.39, 0.86 | +0.39, 0.84 | +0.39, 0.84 |
| UE06 (read swapped), floor | 32, shadow | −1.26, 2.30 | −1.23, 2.88 | −0.39, 2.12 |
| UE06 (read swapped), drone | 31, shadow | −0.05, 0.75 | −0.51, 2.35 | +0.17, 1.35 |
| UE07 (read swapped), floor | 32, shadow | −0.98, 2.32 | −0.77, 2.66 | −0.56, 2.53 |
| UE07 (read swapped), drone | 32, shadow | −0.22, 1.18 | −0.03, 2.30 | −0.19, 1.11 |
| UE06 (as labelled), floor | 32, shadow | −2.01, 5.04 | −1.99, 5.31 | −1.15, 4.61 |
| UE06 (as labelled), drone | 32, shadow | +2.45, 3.86 | +1.97, 3.63 | +2.63, 4.00 |
| UE07 (as labelled), floor | 32, shadow | −0.23, 4.25 | −0.02, 4.04 | +0.19, 4.18 |
| UE07 (as labelled), drone | 31, shadow | −2.72, 4.09 | −2.50, 4.04 | −2.72, 3.97 |
| UE13, floor (walls 0.05, 90 ms) | 26, shadow | +2.04, 4.92 | +2.82, 5.91 | +2.79, 5.91 |
| UE12, floor (walls 0.05, 90 ms) | 27, shadow | +1.23, 3.82 | +0.66, 4.09 | +0.63, 4.08 |
| UE12, drone | 21, shadow | +0.36, 1.13 | −0.04, 1.13 | −0.04, 1.12 |
| UE14, floor | 32, shadow | −0.05, 1.34 | −0.14, 1.60 | −0.11, 1.48 |
| UE14, drone | 24, shadow | +0.16, 0.74 | +0.02, 0.97 | +0.15, 0.73 |
| UE13, drone | 28, shadow | −0.20, 1.57 | −0.60, 1.64 | −0.59, 1.65 |

Behind the second row of both UE12 and UE13 the measurement is louder than the full-wave field of
the LoD1 rows, by 3-5 dB at 1.6-2.5 kHz and more above (+4.4 dB mean in UE12, +8.5 in UE13); that
is not explained.

The drone's runs are the cleanest test, the source far from the floor and from anything that
might be under it: over UE00, UE06 and UE07 (the last two read swapped) the full-wave field is
the measurement to 0.7, 0.7 and 1.1 dB rms, band by band. Where one edge brings the drone's sound
to each receiver -- UE00, UE12, UE13 -- the ray model is as good; over UE06 and UE07, where the
paths cross several houses' edges and ends, UTD is 2.3 and 2.3 and the edge source 1.4 and 1.1. As labelled, the full-wave field
misses UE06 and UE07 by 3.9 and 4.1 dB rms and the ray model by 3.6 and 4.0, so the full-wave
field favours the swap. The ray model's floor-source rows favour it too (UE06 2.9 swapped against
5.3 as labelled, UE07 2.7 against 4.0 dB rms).

In UE06 the FDTD's strongest arrival behind the tall house is a reflection off the low neighbour's
roof and the next tall house's side wall, a corner whose double reflection misses its wall by 1.4 cm
and so arrives at about half its strength. It meets that roof between two edges. Until the model
took a reflection between two edges (`--max_between_reflections`) it was 5.7 dB low at the
receiver farther behind (R2) and 2.3 dB low at the nearer (R1); with it, R2-R4 are within 1.6 dB
rms. At R1 the model's ray reflects 3 cm inside the wall's end, where the FDTD has it at about
half: its partner would diffract at the tall house's corner, which meets the low roof's front edge
at a vertex. The front edge's end there is that vertex's wave, which UTD takes (THEORY §6.5, "An
edge's ends"): R1 is −1.1 dB mean (1.7 rms) by UTD. With `--edge_source` the ray takes the integral
unfolded about the wall, over the front edge as far as it runs: R1 is −0.1 dB mean (1.1 rms) and
the scene 1.0 dB rms. Over EDGE's 13 runs, 396 band levels, UTD is 1.61 dB rms against the FDTD
and the edge source 1.32.

What is left after the FDTD is mostly the floor source itself. Its own floor reflection makes a
null where the path over the first roof leaves at θ, near c / (4 h sin θ) (h = 29 mm): 4.9 kHz in
UE13, 6.6-7.0 kHz in UE06 and UE12, 10 kHz in UE00. In the nine floor-source pairs whose null
falls within 1.6-8 kHz, the measured level at the null's band is 7.8 dB above the pair's mean
residual (median), and in the low bands it keeps the reflection: leaving the source's floor
reflection out lowers 1.6-2.5 kHz by 4-5 dB, which the measurements do not do. Something under
the loudspeaker absorbs at high frequency and not at low; the archive's photographs show each
floor loudspeaker on a black pad.

### Scattering walls (point source)

`niceops_point_field --scattering` (`DiffuseFacades::point_energy`). A source and microphone
together 1 m in front of a 10 m square wall, patches of 0.2 m: within 0.2% of Lambert's law,
ρs d²/π ∫ dA/r⁶ over the wall. The same 2 m above the middle of a 20 m square roof: within 0.1% at
patches of 0.4 m, 1.6% low at 1 m and 10% low at 2 m, where a patch's parts are summed only within
three of its sizes.

**The roofs' share.** A district of sixteen blocks 20 m square and 9 to 27 m high on 12 m streets,
s = 0.5 on every wall and roof, four source positions 60 to 150 m up, 63 Hz to 4 kHz; the roofs'
scattered energy returned, against the roofs reflecting ρ(1 − s) and returning nothing:

| Microphones | Total, mean | Total, worst | Worst band |
| --- | --- | --- | --- |
| 2 m above a roof (16 pairs) | +1.16 dB | +1.94 dB | +7.41 dB |
| At a window, 0.5 m off a wall 4-20 m up (32) | +0.05 | +0.21 | +0.33 |
| In the street, 1.2 m up (16) | +0.02 | +0.04 | +0.12 |

A rooftop microphone sees the roof it stands on and its neighbours'; a street microphone sees no
roof, and what reaches it from one goes by a wall above it. Against walls and roofs specular
(s = 0) the rooftop microphones move by −0.15 to +1.06 dB with the roofs' share returned.

**The reliefs' coefficients** (`tools/profile_scattering.py`, boundary elements at 10 a
wavelength, four frequencies a band, incidence every 5° from −85 to 85°):

| Band (Hz) | 1600 | 2000 | 2500 | 3150 | 4000 | 5000 | 6300 | 8000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EDGE scatter_1, sawtooth 0.12 m deep on a 0.36 m pitch | 0.91 | 0.94 | 0.95 | 0.96 | 0.97 | 0.98 | 0.98 | 0.99 |
| EDGE scatter_2, ribs 24 mm deep on an 80 mm pitch | 0.05 | 0.06 | 0.30 | 0.62 | 0.76 | 0.53 | 0.19 | 0.15 |

A 20 ft container's side wall (trapezoids 36 mm deep on a 278 mm pitch, as one owner describes his;
the standard does not fix it): s under 0.01 to 500 Hz, then 0.10 at 630 Hz, 0.26 at 800, 0.35 at
1 kHz, 0.45 at 1.25, 0.55 at 1.6, 0.69 at 2, 0.83 at 2.5, 0.91 at 3.15, 0.88 at 4, 0.80 at 5 and
0.82 at 6.3 kHz.

**EDGE's deep canyon with relief** (UE24, UE25 and UE26 against UE23: walls 1.92 m high and 2 m
apart, sources and receivers on the floor in the canyon, all in line of sight; S2, 1.6-8 kHz; the
change from UE23, which needs no calibration). The model's walls take the relief's s(f), the rest
of the canyon its UE23 geometry, eight reflections:

| Scene | Pairs | Measured | Model | Specular rays alone (the diffuse energy left out) | Error, model | Error, no change |
| --- | --- | --- | --- | --- | --- | --- |
| UE24 (sawtooth on one wall) | 64 | −1.03 dB mean, 1.80 rms | −0.25, 0.82 | −3.46, 3.58 | +0.78 mean, 1.46 rms | +1.03, 1.80 |
| UE25 (ribs on the other) | 64 | −1.11, 1.92 | −0.04, 0.32 | −1.13, 1.47 | +1.07, 1.84 | +1.11, 1.92 |
| UE26 (both) | 64 | −1.62, 2.50 | −0.13, 0.83 | −4.05, 4.20 | +1.49, 2.23 | +1.62, 2.50 |

In a canyon Lambert's law returns almost all that the reliefs take from the specular rays, where
the measured levels fall by 1-1.6 dB: the model lies between keeping the energy and losing it, and
nearer keeping it. The full-wave field says why. With the sawtooth as rigid geometry (UE24, its
teeth facing one way along the canyon, which the archive does not give, walls absorbing 0.05, 4 mm cells, 1.6-6.3
kHz, the floor source at the canyon's middle; 28 band levels), the FDTD's change from UE23 is
−1.03 dB mean, 1.52 rms, against the measured −0.52 and 1.35: it is within 1.02 dB rms of the
measured change band by band, where Lambert's model is within 1.19 (+0.57 mean). The large
changes are a single receiver's in a single band -- −3.1 measured at 1.6 kHz against −3.0 from the
FDTD and −0.2 from Lambert's law; −2.5 at 2 kHz against −2.9 and −0.7 -- the relief redirecting
the reflections rather than diffusing them, which a coefficient cannot hold. The ribs (UE25, the
same way): the FDTD's change is −0.64 dB mean, 1.07 rms, against the measured −0.63 and 1.45,
within 1.30 dB rms band by band; Lambert's model, +0.04 mean, is within 1.46, which is no better
than assuming no change at all (1.45). The largest measured changes there, +2.8 and +3.3 dB at
1.6 kHz where the ribs scatter almost nothing, are in neither. The second source, half a metre
off the canyon's middle, says the same (28 band levels each): the sawtooth's measured change is
−1.77 dB mean, the FDTD's −1.67, within 1.45 dB rms, and Lambert's model −0.50, within 1.83; the
ribs' −1.68, −1.21 within 1.85, and −0.10 within 2.31. Both reliefs together (UE26, the first
source): −0.87 measured, −1.16 from the FDTD within 1.60, +0.16 from Lambert's within 1.89.

So the coefficient does what it is for where the diffuse energy is a path the rays lack -- a
microphone in a building's shadow, as DigiCity's below -- and not where it only redistributes
what the specular rays already carry, as in a street with the source and the microphone in it.

**EDGE's park with relief** (UE28 against UE27; scatter_2 on the two facades that face each other
across the park, read from the archive's photographs, each a 1 cm slab carrying s(f)): measured
changes are small, 0.86 dB rms in line of sight and +0.40 dB mean in shadow; the model's error is
0.75 dB rms and −0.38 dB mean. The largest measured change, the drone to receiver 3, −2 to −3 dB at
2.5-5 kHz, the model has at −0.5 to −1.4 dB, peaking in the same bands as the ribs' s.

**DigiCity's microphone 7** (the loudspeaker at hover spot 5; the source calibrated band by band
on microphones 2, 8 and 10, as the validation harness's DigiCity comparison does; walls absorbing
0.05, 70% RH, patches of 1 m):

| Walls | Microphone 7, predicted − measured, 160 Hz-5 kHz | Calibration microphones' spread |
| --- | --- | --- |
| Specular | −4.1 dB mean, 5.2 rms (−10.4 to +1.3) | 2.51 dB |
| s = 0.1 in every band | −3.4, 4.1 | 2.36 |
| s = 0.3 | −3.0, 3.4 | 2.34 |
| s = 0.5 | −3.1, 3.4 | 2.44 |
| The corrugation's s(f) | −2.3, 2.9 (−6.4 to +1.3) | 2.15 |

With the corrugation's coefficients microphone 7 rises by 1.0-7.7 dB from 800 Hz to 5 kHz (1.25 kHz
−7.8 → −2.9 dB, 5 kHz −10.4 → −2.7). Below 630 Hz nothing moves, so the shortfall at 200-250 Hz
(−3.4 dB) stays.

**Reliefs as gratings** (`--relief`, `--facade_relief`; THEORY §6.5). The same scenes with each
relief's grating orders (`tools/profile_scattering.py --orders`: the sawtooth 3.72 m long, the ribs
4 m, the container wall 6 m, four frequencies a band, incidence every 5°, the far field every
0.1°), walls absorbing 0.05 throughout; the validation harness's `urban/edge/relief.py` and
`urban/digicity/accuracy.py`. Each cell is the mean change from the scene without relief and, in
brackets, the error against the measured change, mean / rms (for DigiCity, microphone 7's
predicted − measured):

| Scene | Measured | FDTD | Random-incidence s(f) | Orders | Lambert at incidence | Hybrid (default) |
| --- | --- | --- | --- | --- | --- | --- |
| UE24, sawtooth (64 pairs) | −1.03 | −1.35 (−0.20 / 1.25) | −0.22 (+0.81 / 1.48) | −1.50 (−0.47 / 1.26) | −0.25 (+0.78 / 1.46) | −0.74 (+0.29 / 1.25) |
| UE25, ribs (64) | −1.11 | −0.92 (+0.24 / 1.60) | −0.02 (+1.09 / 1.85) | −0.50 (+0.61 / 1.59) | −0.05 (+1.05 / 1.80) | −0.23 (+0.88 / 1.70) |
| UE26, both (64) | −1.62 | −1.16 (−0.29 / 1.60) | −0.10 (+1.52 / 2.25) | −2.01 (−0.40 / 1.66) | −0.16 (+1.46 / 2.16) | −0.75 (+0.87 / 1.80) |
| UE28, in sight (72) | +0.23 | | +0.03 (−0.20 / 0.77) | +0.01 (−0.22 / 0.70) | +0.02 (−0.21 / 0.76) | +0.01 (−0.22 / 0.70) |
| UE28, shadow (22) | +0.40 | | +0.03 (−0.36 / 0.60) | +0.06 (−0.34 / 0.58) | +0.08 (−0.32 / 0.56) | +0.06 (−0.34 / 0.58) |
| DigiCity microphone 7 | | | −2.33 / 2.93 | −3.26 / 4.15 | −1.76 / 2.46 | −1.91 / 2.64 |

In the canyon the orders carry the relief's redirection, and their error is the FDTD's: what the
sawtooth sends out along the canyon or back is lost to the microphones, where Lambert's law returns
most of it. The sawtooth's teeth may face either way along the canyon (the archive does not say);
turned round, its mean change moves by 0.06 dB. Behind DigiCity's containers the order is reversed:
microphone 7 is in a shadow, and what the walls scatter onto other walls reaches it, which the
diffuse field carries and the orders alone do not. The hybrid, each wall's first scattering as
orders straight to the microphone and the rest diffusely, is better than the random-incidence s(f)
in every case, and in the canyon level with the FDTD (UE24 1.25 against its 1.25). It was 1.39,
1.76 and 1.98 while its grating rays could still reflect in another wall on the way to the
microphone: the diffuse field re-emits whatever a relief scatters onto other walls, so those
reached the microphone twice. The orders' shares are now their means over each band's three
sub-bands, as the specular share is taken, and a floor bounce's coefficient goes with its phase:
0.0-0.07 dB in the orders' column. Two choices were measured on the way: grating rays and their
floor partners add with their phases (in energy, the canyon's sawtooth came out −2.29 dB, the
partners near the floor being one length), and rays that diffract reflect in a relief whole
(taking its specular share, DigiCity's microphone 7 fell to −5.9 dB mean with the orders; EDGE did
not move). Grating rays with up to six other reflections change nothing in the canyon, every second
one being on the relief, whose specular share is a few per cent.

**A relief as geometry.** The sawtooth is 12 cm deep on a 36 cm pitch, so the canyon's microphones,
0.4-1.5 m from it, see its facets rather than a grating (their Fresnel zones span about a period).
Modelled as its own facets and ridges (the FDTD's solid, by UTD, eight reflections), UE24's mean
change is −0.88 dB, error +0.16 / 1.28, as good as the FDTD and the orders. The ribs that way are 24
million rays a pair, and sub-wavelength.

### Pitched roofs

`--roofs` (LoD2; THEORY §6.5). The geometry is checked directly (`tests/test_pitched_roofs.cpp`): a
gable house's walls, eaves, rakes and ridge and each wedge's angle; blocking against sampling the
segment at 4000 points over 4000 random segments, above and below the floor; shed, hip and mansard
roofs; the roofs file's units, cap and refusals; and pruning, which finds the unpruned search's rays
ray for ray among gable, hip, shed and mansard houses.

**Against exact solutions.** A house long enough to be two-dimensional, the source and the
microphone in one plane across it, every ray times the root of its length as for the long block
(boundary elements by `tools/thick_barrier_bem.py`, rigid house and floor, 90 elements a wavelength,
converged to 0.15 dB; the cases that were not are left out), 125 Hz to 1 kHz, sources 2-4.5 m up
and 10-30 m off, microphones 2-15 m behind and 0.5-7 m up:

| House | Cases | UTD: worst, rms (dB) | Edge source | Flat roofs instead, rms |
| --- | --- | --- | --- | --- |
| Gable 8 m across, eaves 6 m, ridge 9 m | 71 | 0.74, 0.17 | 0.53, 0.12 | 6.2 at the eaves, 8.9 at the mean height, 11.0 at the ridge |
| EDGE UE03's terrace at model scale (0.96 m across, eaves 0.492 m, ridge 0.732 m), the floor source, 1.6, 3.15 and 6.3 kHz, 45 elements a wavelength | 12 | 0.34, 0.15 | 0.47, 0.16 | 3.0 at the eaves (worst 6.5) |
| Mansard 12 m across, eaves 6 m, slopes 1 and 0.4, deck 9 m | 68 | 13.8, 3.9 | | |

The gable's path over it is eave, ridge and eave, or the ridge alone where the source sees the near
slope. The mansard's, from a source 0.5 m up to a microphone 2 m behind, bends at four edges -- the
near deck edge, the far one, the far break and the far eave -- and without a chain of them there is
no ray at all: chains (`max_chain`) bring those nine cases within 1.5 dB below 4 m (at 4 m the
microphone sits on the far lower slope's plane, at its eave's boundary). The mansard's other cases,
three edges, are 3.5 dB rms: its wedges are nearly flat ($n \approx 1.1$) and each pair shares the
face the ray runs along, where the cascade is weakest. The uniform form that fixes a flat roof's
pair (above) is for two edges; a cascade of transitions (Tzaras and Saunders's heuristic, or the
edge source's higher orders) would be the next step.

**Why not the FDTD.** A pitched roof in the FDTD is a staircase of cells. With the terrace's gable
at 6 mm cells, the FDTD behind it was 2-13 dB louder than the rays, rising with frequency, and
louder than the bare block at the eaves, while the same FDTD matched the rays within 0.5 dB with
the roof flat at the eaves or at the ridge. Built as 40 stacked slabs it gives the voxels' result to
0.3 dB, so the voxelization is consistent. The boundary elements, given the voxels' own staircase in
the profile's plane, put it there:

| Microphone, 6.3 kHz | Smooth slope (rays within 0.17) | Staircase of 6 mm | 3 mm | 1.5 mm |
| --- | --- | --- | --- | --- |
| 0.46 m behind, on the floor | −27.0 dB re free field | −16.8 | −21.9 | −24.5 |
| 1.25 m behind | −17.0 | −9.8 | −13.2 | −15.0 |

The error halves with the step: first order. A staircase's steps diffract a grazing wave down the
slope that a smooth slope does not, so the FDTD at 3 mm, the finest the EDGE runs used, is still
about 5 dB loud in the top bands; 1.5 mm does not fit the GPU. The FDTD is no reference for a sloped
roof at these cells, and the pitched scenes are scored against the measurements alone.

**Against EDGE's measurements.** The harness (`urban/edge/`) builds the pitched scenes from their
STEP files: the roof and sloping faces grouped into convex roofs, each a building, so a step (UE01,
UE02) or a valley (UE05) parts them into neighbours sharing a side; UE04's mansard is four planes
under its deck. UE19 stands on a carpet and is left out, as the carpeted LoD1 scenes are. Every
receiver is in a house's shadow, so each loudspeaker takes the per-band calibration of the LoD1
scenes' line-of-sight pairs. Measured − predicted, mean and rms (dB), 1.6-20 kHz at model scale:

| Scene, source | Band levels | Roof planes, UTD | Edge source | Flat at the mean height | Flat at the ridge |
| --- | --- | --- | --- | --- | --- |
| UE01 half gable, floor | 32 | +0.34, 2.43 | +0.44, 2.35 | +7.17, 8.43 | +9.04, 10.17 |
| UE01, drone | 48 | +1.38, 3.28 | +1.56, 3.31 | −4.14, 4.52 | +0.15, 3.08 |
| UE02 pent, floor | 32 | −0.33, 2.88 | −0.17, 2.83 | +5.03, 6.41 | +9.72, 10.67 |
| UE02, drone | 48 | +1.90, 3.30 | +2.12, 3.40 | −2.25, 3.00 | −0.19, 2.02 |
| UE03 gable, floor | 32 | −0.21, 2.00 | −0.27, 2.06 | +6.44, 7.58 | +9.50, 10.54 |
| UE03, drone | 48 | +0.70, 2.49 | +0.64, 2.46 | −0.67, 2.70 | +1.62, 3.92 |
| UE04 mansard, floor | 31 | −2.26, 4.47 | −2.26, 4.47 | +5.01, 6.20 | +7.10, 8.14 |
| UE04, drone | 48 | +1.27, 2.40 | +1.00, 2.37 | −0.01, 2.65 | +1.48, 3.52 |
| UE05 shed, floor | 31 | −0.52, 3.16 | −0.47, 3.17 | +3.07, 4.64 | +5.07, 6.41 |
| UE05, drone | 48 | +0.25, 1.85 | +0.20, 1.85 | −1.43, 2.21 | +0.18, 2.00 |
| UE15 detached gable, floor | 30 | +0.28, 1.83 | +0.24, 1.73 | +0.81, 2.71 | +1.24, 3.44 |
| UE15, drone | 48 | +0.07, 2.24 | +0.11, 2.18 | −0.40, 3.20 | +0.96, 4.15 |
| UE17 detached gable (1/5), floor | 31 | +0.58, 1.98 | +0.58, 1.98 | +2.25, 4.19 | +2.33, 4.11 |
| UE17, drone | 47 | +1.02, 4.69 | +1.01, 4.69 | +0.53, 3.85 | +1.99, 4.94 |
| All | 554 | +0.45, 2.94 | +0.47, 2.94 | +0.97, 4.57 | +3.04, 5.83 |

From the floor, where the roof's shape sets the path, the planes are 1.7-4.5 dB rms where a flat
block is 2.7-10.7; from the drone, which sees over the roofs, they and the flat block at the mean
height are alike. UE04's floor source was 7.8 dB low, 14.3 rms, before chains: its microphones
close behind are reached over four edges. UE17's drone source is the one row a flat block does
better, by 0.8 dB.

**What the drone rows leave, receiver by receiver.** The part of the drone rows' error that is the
model's is the pair over a step, UE01's and UE02's: by the cascade the boundary elements put it
0.2-1.6 and 4.0-5.0 dB low, against the measurements' 4.3-4.9 dB high at UE02's R2, and the
uniform form within 0.53 dB at R2's distance ("EDGE's stepped terraces", above). What remains is
receiver by receiver, and may be the receivers':

- **R3 and R4 above 8 kHz.** R4 reads up to 10-13 dB high at 20 kHz and R3 up to 5 dB low in every
  pitched scene, UE03 included, where the model matches the exact 2-D field. EDGE's equipment table
  gives R1 and R2 as B&K 4136 pressure-field capsules and R3 and R4 as G.R.A.S. 40BF free-field
  ones, and the calibration is per loudspeaker and band, not per receiver: the trends may be the
  40BF's response off its axis, or the receivers' positions. The data do not say which.
- **R1, 1-2 dB high in most pitched scenes.** Against the FDTD on UE15 (second order, 3 and 2 mm,
  extrapolated to zero cell for its stepped slopes, which are first order in the cell), the model is
  within 0.4 dB at R3 and R4 and ±2.3 dB at R1 and R2, the FDTD reproducing the measurements' band
  shape at R1 (+2.3 dB at 2.5 kHz, measured +3.6). At R1 the FDTD is 0.45 dB above the model on
  average over 1.6-8 kHz where the measurements are 2.25 dB above: about 1.8 dB lies outside the
  full-wave field.
- **UE17's R1** (the detached house at 1/5 scale, 1 m behind it on its axis): measured minus
  predicted +7 dB, rising to +12 at 20 kHz, its level falling with frequency like one diffraction
  where the model's paths diffract twice. The FDTD at these cells does not fit a scene that size;
  it is open.

**Cost.** Chains cost nothing measurable: a suburb of 64 gable houses (10 by 15 m, eaves 6 m, ridge
9 m, 25 m apart) to 162 microphones in its gardens took 166.1 s with them and 165.9 s without from
a source 1.5 m up, which made 26 chain rays and moved no band level by 0.1 dB. The pitched houses
cost 2.2 times the same houses flat (1408 against 637 s of CPU from 60 m up), most of it their
eleven edges against eight; before the pruning took their roof planes it was 3.2 times (2051 s),
with the same levels to the bit.

**Slopes reflect and scatter.** A slope is a mirror in its plane, with the floor at either end and
on either side of a once-diffracted ray, and scatters as the walls do (THEORY §6.5). EDGE's
pitched rows above do not move with it: every microphone there is low in a house's shadow, where no
slope reflects, and the drone is below the plane of the slope that faces them. In the suburb above,
from 60 m up, 2 of the 162 garden microphones take a slope's specular reflection -- a neighbour's
slope reflecting into the gap between two houses, the reflection point 0.6 m from its ridge -- and
move by up to 5.2 dB in a band; with the once-diffracted rays that reflect in a slope, 197 of the
810 band levels move by more than 0.1 dB, from −3.7 to +1.8 dB, the largest at 125 Hz. Those rays
are a slope's reflection reaching an edge -- a corner, an eave -- whose ridge's diffracted wave was
there without the reflection it is half of. 115 s on the GPU, 110 without either. The boundaries:

| Boundary | Step with the slope's rays | Without them |
| --- | --- | --- |
| The slope's own edge (eave, ridge) cuts off its reflection | under 0.02 dB | the reflection's level |
| Another building's edge cuts off the wave to the slope, or from it, 0.1 mm either side | under 0.05 dB | 8.1 dB |

(`a slope's reflection ends where its edges' boundaries are`, `an edge that cuts off a slope's
reflection leaves the field continuous`.) A slope's diffuse field is Lambert's law over the tilted
surface to 0.5% (`a point source lights a slope's diffuse field as Lambert's law says`).

**Walls with a slope.** Composite-image checks for a wall perpendicular to a slope and for
oblique mirrors hold the specular length to round-off, including the reciprocal path and the
same paths with diffraction disabled. Specular and once-diffracted floor partners have real
legs summing to their ray length, reciprocal lengths, and the correct floor-contact launch.
A screen cutting off a wall–slope reflection leaves less than 0.05 dB between microphones
0.1 mm either side of its boundary, in both directions, at 125–2000 Hz; this also holds with
only one reflection allowed on each side of the edge.

An on/off check of the added slope–wall paths on EDGE's UE01–05, UE15, UE17 and UE19, all
provided source/receiver pairs at 1.6, 3.15, 6.3, 12.5 and 20 kHz, leaves all 360 band levels
unchanged. This is an omnidirectional, rigid, no-air check of path impact, not a new fit to the
loudspeaker measurements.

**Twice-diffracted rays off a roof.** Where a second edge cuts off a once-diffracted ray's roof
reflection -- a tower's corner, then a sharp corner across the leg to a lower roof, flat or
sloped -- the ray diffracted at both corners and reflected in the same roof carries the field
across the boundary: the field steps by 1.9% (flat; the rays without the roof switch there too)
and 1.4% (slope) of the cut-off ray's amplitude between microphones 0.1 mm either side, at
125-2000 Hz, both ways round, and by 100% without these rays (`a roof's reflection after two
edges carries on where a second edge cuts it off`). The percent is UTD's own at a second edge's
shadow boundary: the same handover without a roof leaves the same. In dB it depends on the rest
of the field; at a 125 Hz null there it is 0.3 dB. On the gable suburb (162 street
microphones, source at 60 m) they move 236 of 810 band levels by more than 0.1 dB, at most
2.4 dB (125 Hz, at a microphone in a near-null), mean +0.07 dB, with 7.5% more rays and no
measurable cost. Over DigiCity's flat containers (20 drone positions, 5632 points, 17 bands) they
move 1.6% of the band levels by more than 0.1 dB, at most 3.1 dB, mean +0.0004 dB, for 6% more
time (122 to 130 s).

### Microphones on roofs and walls

A source or a microphone lying on a wall or a flat roof (THEORY §6.5, "Ends on a surface"), against
the same a hair (10 µm) off the surface, where the model traces the surface's partners itself.
A rigid building 20 m square and 10 m high, `niceops_point_field`, band levels re the source at 1 m:

| End | Before, on the surface | Now, on it | A hair off it |
| --- | --- | --- | --- |
| Microphone on the roof, source 50 m overhead (250 Hz) | −33.98 (free field) | −27.87 | −27.87 |
| Microphone on the roof 5 m behind its edge, source low beyond (250 Hz) | −300 (no ray) | −32.28 | −32.27 |
| Microphone on a wall, source facing it 40 m off (250 Hz) | −26.70 | −20.71 | −20.71 |
| Source on the roof, microphone 50 m overhead (250 Hz) | −33.98 | −27.87 | −27.87 |

Over eight geometries, overhead, behind and obliquely over the roof's edge, facing a wall, round
its corner, behind the building and oblique to it, six bands from 125 Hz to 4 kHz: within
0.062 dB, the largest at 125 Hz in the building's shadow, where a hair off the face carries an
edge's backscatter along it that the end on the face leaves out. An absorbing roof (ρ = 0.6)
overhead gives the free field times (1 + √ρ)² to 0.05 dB, as above it. The GPU's levels for a
microphone on a roof with a plate are the CPU's (`[gpu]`). On a gable's slope, five sources
overhead, in front, below the eave, beyond the ridge and oblique, the end on the slope and the
source on it, within 0.155 dB of a hair above it, the largest at 125 Hz beyond the ridge from the
same backscatter (the ridge's wave down the slope, off the eave and back). In the plane-wave
tables a receiver on a wall is one a hair off it, the wall's rays traced, to 0.001 dB over rigid
walls in every direction, either side of grazing included.

**Each microphone's floor.** On the demo district and a 4 s stretch of the demo track, observers
raised 5 ft (on no roof or wall) give the levels the same observers at z = 0 give under the track
lowered 5 ft, to 1e-10 dB (`buildings_floor_at_observer_z`); with half of them raised, each half
is its own floor's run to the digit.

### Trajectory smoothing

The demo track's own latitude, longitude and altitude columns, fed back in as a WGS84 track and
compared with the footprint from its full inertial state (S-76D demo case, 14,400 cells):

| Input | Median | 99th pct | Max | Cells over 1 dB |
| --- | --- | --- | --- | --- |
| 50 Hz positions, default `--path_smoothing 10` (ft) | 0.073 dB | 0.72 | 1.42 | 25 |
| 50 Hz positions, `--path_smoothing 0.0656` (0.02 m) | 0.031 dB | 0.36 | 1.44 | 6 |
| 50 Hz positions, `--path_smoothing 0` (exact) | 8.0 dB | 9.57 | 10.58 | 14400 |
| 50 Hz positions and velocity columns, default | 0.015 dB | 0.36 | 1.42 | 4 |
| Waypoints 10 s apart, resampled to 0.1 s, default | 0.121 dB | 0.86 | 2.10 | 76 |

The projection itself is not the lossy part: it reproduces the file's own x and y to under a
meter over 3.2 km, and altitude passes through exactly. What costs is differentiating twice.
The default is the only setting good on every input shape at once; 0.0656 ft does better on the
dense track only because it matches that file's coordinate resolution. Exact interpolation is
uniformly poor on measured input.

Fitting a dense track's ends freely (over a 10 s quadratic extension) takes the default row from
0.109 dB median and 84 cells over 1 dB to the values shown. On a constant-acceleration test path
(1.5 m/s²) a natural spline's acceleration is zero at the end, 0.50 m/s² four seconds in and
1.14 eight seconds in; the freed fit holds it within 11% throughout.

## Measured validation

### Data sets

The measured comparisons use the 2017 Noise Abatement flight test (data report
NASA/TM-2019-220264) unless they say otherwise. The others are the 2014–15 altitude-variation
test (NASA/TM-2016-219354: an AS350 SD1 over inverted ground boards at Amedee, Salton Sea and
Sweetwater, where microphones 27–29 stand behind a hill and 1–21 on a 3.4° slope) and the
1978 T-38A runway test at Wallops (NASA TP-1747, microphones on concrete and grass). At
Sweetwater the recordings labeled 28 and 29 were made at each other's surveyed positions, in every
run: their rotor tones de-Doppler only at the exchanged positions, and stationary sources heard with
no aircraft flying give delays no source could produce with them where labeled. Every Sweetwater
result here uses each recording where it was made.

- **Aircraft and sites.** The B407 flew at Amedee, where the ground is flat to within about a
  meter over the test range. The R44 (database RO-44) and R66 (RO-66) flew at Eglin, along a
  140° layout. The AS350, Be206L3 and EC130 complete the set.
- **Microphones.** Flush ground plates (GR1425 boards, whose top is 8 mm up; the boards sit
  3–10 cm below the survey datum) and 1.2 m pole microphones.
- **Weather.** Ground stations (anemometers SWS 1, 2 and 4), LIDAR winds (first gate at 46 ft)
  and balloon temperatures. The station files carry no unit and no anemometer height: knots
  (as the data report gives them) and 2 m are a reading, not a measurement.
- **B407 calm day, day 286.** 25 approaches (114,720 trajectory states), scored with the shipped
  database: 1111 microphone-runs. Its mornings have inversions:
  `atmospheres/example_inversion.csv` is run 286276's profile without its 1–3 kt of wind, 3.5 K
  of inversion capped near 350 ft. The "day-286 PE matrix" is the PE run on five of its profiles.
- **B407 windy day, day 284.** 16 steady runs flown west, scored with the held-out database
  `min2_no284`, built without them: 729 microphone-runs.
- **The four detail runs**: calm-day approaches 286270 and 286277 (shipped database) and
  windy-day steady runs 284205 and 284211 (`min2_no284`), with band histories at every
  microphone.
- **Other aircraft**: each aircraft's approaches, which no database was built from (AS350 day
  291, Be206L3 day 279, EC130 day 298, RO-44 day 229, RO-66 day 236), and the R44's held-out
  turns.

Configurations, as the validation study names them: `station_plate` is straight rays through
the station atmosphere to flush plates; `stratified_plate` is `--refraction stratified` through
the measured LIDAR and balloon profile to flush plates.

Scoring: SEL is the 10 dB-down exposure. Unless stated, scores are over the study's usable
microphone-runs (SNR of at least 10 dB, records that miss the flyover left out), and 95%
intervals come from a bootstrap that resamples whole runs.

### Shipped databases: approach scores per aircraft

NICEOPS minus the measured 10 dB-down SEL over each aircraft's approaches, station atmosphere,
with the shipped database:

| Database | Spheres | Approaches | Median / rms |
| --- | --- | --- | --- |
| AS350 | 96 | day 291, 1328 mic-runs | +0.70 / 1.49 dB |
| Be206L3 | 110 | day 279, 906 mic-runs | +0.24 / 1.50 |
| EC130 | 112 | day 298, 1025 mic-runs | +0.70 / 2.70 |
| RO-44 | 83 | day 229, 585 mic-runs | +0.04 / 1.16 |
| RO-66 | 88 | day 236, 876 mic-runs | +0.87 / 1.56 |

- **Every aircraft over-predicts at low elevation.** The AS350 is +1.59 / 2.29 dB at 5–10° of
  elevation, and carries a +0.6 dB bias overhead. The RO-44 is 1.94 dB rms at 0–5° and 1.23 dB
  at 5–10°.
- **EC130.** The day-298 approaches have no ambient run, so their measured SEL carries the
  background and cannot be filtered on SNR. Read its rms, not its median.
- **AS350 and RO-44** were scored on approaches flown on a day that also gave spheres (days 291
  and 229).
- **The B407** is scored on day 286 in
  [Source options](#source-options-blend-and-load-factor-law). Its hover's two measured rings
  score 2.8 dB held-out rms (18–25° below the horizon) and 5 dB (0.5–3°, in
  ground effect). Nothing measures 3–17° or the region under the aircraft.
- **Coverage.** On the demo grid, the B407 database's levels come 0.94% from unmeasured cells at
  the median observer and 31% at the 95th percentile. The builder's SNR gate drops the
  5–10 kHz bands in about half the measured cells, 0.4% of the lower hemisphere's energy.

How these scores compare with earlier builds of the same databases is in
[HISTORY](HISTORY.md#the-2017-noise-abatement-databases).

### Database screen

The load-time screen ([MANUAL](MANUAL.md#screening-a-database)) on every shipped database, and on
the NASA R44 spheres built for the cross-code comparison. Margins: a dBA cell more than 40 dB from
the median of its sphere, a band more than 45 dB above the median of its cell's spectrum.

| Database | dBA cells flagged (spheres) | Bands flagged (spheres) | Worst entry, against the median |
| --- | --- | --- | --- |
| `B407_M2` | 240 (24 of 1248), all low | — | 27.1 dBA against 89.3 |
| `B407_M3` | 3928 (72 of 1296), all low | — | 0 dBA (a written zero) against 89.6 |
| `S-76D_M4` | 1376 (32 of 1312), all high | — | 1463.5 dBA against 108.3 |
| `S-76D_fixed_load_factor` | 172 (4 of 164), all high | — | 1449.5 dBA against 94.3 (μ 0.208, γ −7.8°) |
| NASA R44 spheres | 2397 (86 of 122), all high | 5380 (102 of 122) | 170.2 dB at 10 kHz against 70.1 |
| The other 17 shipped files | 0 | 0 | |

- **The margins clear every sound database.** On the 17, no stored cell sits more than 29.1 dB
  above its sphere's median (B205) or 21.4 dB below it (B407), and no band more than 35.9 dB above
  its cell's median (RO-44). Most S-76D cells flagged are whole azimuth columns 40–60 dB high:
  sphere 33 of `S-76D_fixed_load_factor` reads 135–139 dBA along φ 100° against a median of 91.
- **Two rules that were tried do not separate.** Against the median of its eight grid neighbors, a
  sound cell beside a strip of written zeros reads 44 dB off (224 sound 88 dBA cells in
  `B407_M3`), and a strip two cells wide is its own median. Against the same band in neighboring
  cells, gating leaves sound bands in the AS350 file up to 169 dB below and 116 dB above.
- **The S-76D defects move the demo footprint.** On `S-76D_fixed_load_factor`, masking moves 958
  of the 14,400 nodes by more than 0.1 dB, from −1.8 to +7.5 dB (`S-76D_M4`: 957, the same
  range). The rises are where a spike set the observer's peak: at (3363, −1394) ft its kept peak
  is 69.5 dB, which narrows the 10 dB-down window to 105 samples (SEL 65.8 dB); masked, the peak
  is 62.0 dB, the window 2275 samples and the SEL 73.3 dB.
- **The fill value** (a group never written) is in none of the shipped files.

### Source options: blend and load-factor law

**Condition selection.** On the day-286 B407 approaches (all 114,720 states), AAM's bracketing
search interpolated in both axes at 9.4% of states, took three spheres at 35.8%, two at 34.3%
and one at 20.5%. The barycentric blend's spheres span a median 13.5° of path angle there, and
76% of states take a sphere further away than AAM's gates allow (more than 5° or 20 kt).

`--blend nearest` and `aam` paired against the default barycentric blend (rms over every scored
microphone):

| Case | Metric | barycentric rms | nearest Δrms | aam Δrms |
| --- | --- | --- | --- | --- |
| Day 286 approaches, `station_plate`, shipped B407 (25 runs, 1111 mics) | SEL | 1.26 dB | +0.47 [+0.37, +0.57] | +0.24 [+0.02, +0.48] |
| | LAmax | 1.94 dB | +0.45 [+0.32, +0.60] | +0.34 [+0.18, +0.50] |
| Day 284 steady runs, `stratified_plate`, `min2_no284` (16 runs, 729 mics) | SEL | 1.69 dB | +0.14 [+0.03, +0.23] | +0.42 [+0.17, +0.66] |
| | LAmax | 2.99 dB | −0.38 [−0.61, −0.19] | −0.15 [−0.43, +0.13] |

The barycentric blend is the best of the three on SEL on both days; nearest is worse on every
one of the 25 approaches. The exception is LAmax on day 284, where the barycentric prediction
runs 1.6 dB low at the median: a blend in dB of spheres whose directivity peaks lie in
different directions flattens the peak, and a single sphere keeps it. The SEL integrates over
the directions and does not show that. On average neither alternative moves the default's SEL
much (nearest +0.6 dB on day 286, aam +0.03), but single microphones move by up to 3 dB
(nearest) and 6 dB (aam), mostly below 20° of elevation.

**The load-factor law below 1 g** (`--load_factor_correction two_sided`):

- Normalizing every sphere to the reference C<sub>T</sub> and predicting each day-286 approach
  at its own C<sub>T</sub> with the full law took 0.047 dB off the SEL rms, interval
  [0.012, 0.084]. This was an experiment on the database's spheres, not the option.
- Run directly, the option takes 0.028 dB off the SEL rms on day 286, interval
  [0.007, 0.048], better on 72% of the runs, and 0.018 dB on day 284, interval
  [−0.017, +0.052]. LAmax is 0.05 dB worse on day 286 [+0.03, +0.08] and 0.02 dB worse on
  day 284 [+0.00, +0.04]. The level moves by 0.2 dB rms, at most 1.1 dB.
- On repeated steady runs of one flight condition (169 pairs of B407 spheres whose
  C<sub>T</sub> differs by weight and air density, 0.55 dB sd in 20 lg C<sub>T</sub>), the
  level difference follows the C<sub>T</sub> difference with a slope of 0.34 ± 0.17 per dB of
  the law, fitted together with speed: about a third of the law.

The evidence is modest and mixed, so the default stays one-sided.

### The sphere's nose

On day 284 (16 held-out runs, 729 microphones, `min2_no284`), every run flew west with the wind
6–10 kt from its right and crabbed 4–19° into it, 12° on average, most at the slowest speeds.
`--nose heading` against the same database, which is labeled along the track:

| Configuration | SEL rms, `velocity` | SEL rms, `heading` | Change [95% CI] | 5–10° elevation |
| --- | --- | --- | --- | --- |
| `stratified_plate`, day 284 | 1.69 dB | 1.54 dB | −0.15 [−0.28, −0.02] | 3.02 → 2.26 |
| `station_plate`, day 284 | 2.35 dB | 2.19 dB | −0.17 [−0.28, −0.03] | 5.76 → 5.22 |
| `station_plate`, day 286 (calm) | 1.26 dB | 1.24 dB | −0.02 [−0.05, +0.01] | 1.83 → 1.69 |

- LAmax rms moves by −0.14 dB on day 284 (`stratified_plate`) and −0.01 dB on day 286.
- **Placebo runs**, turning every run's nose by a fixed 13° rather than by its own crab, show
  that each run's own crab is what helps. On day 284 the fixed turn toward the crab gave −0.03
  [−0.22, +0.16] dB and the opposite turn +0.47; on the calm day it cost +0.29.
- Emulating the option by turning each track's velocity onto its heading agrees with it to
  0.01 dB in rms and up to 0.5 dB at a microphone, because turning the velocity also moved the
  path angle relative to the tip-path plane by up to 1.4° in a turn.
- **The hover nose's frame.** Taking the hover nose's offset from the track fit instead of the
  frame would have turned it by the approach's mean crab: 2–7° (median 4.6°) over the day-286
  approaches, had they ended in hover.

### Speed reference

On the 2017 approaches, taking the velocity through the air (`--speed_reference air`) instead of
over the ground moves the advance ratio by 1–3% and the nose by about 1° on calm days, and by
5–10% and 8–15° on a 15 kt day.

On day 284 (`profile_plate`, `min2_no284`), the along-track wind's share of an air reference,
which moves the advance ratio, took the SEL rms from 2.36 to 2.51 dB. The cross-track share,
which mostly turns the nose, took it to 2.25. The 2017 databases were built ground-referenced.

### Wind profile and stability

B407, plate microphones, `--refraction stratified` unless stated; SEL rms. For `--wind`, each
run's wind is the mean of the station anemometers over ±60 s of the run's start, its magnetic
direction turned by the 13.6° declination, taken as 2 m high. Station winds were 4–9 kt from
264–293° on day 284 and under 1.1 kt on day 286.

| Atmosphere | day 284, windy (16 steady runs, `min2_no284`, 729 mic-runs) | day 286, calm (25 approaches, 1111 mic-runs) |
| --- | --- | --- |
| straight rays, station air | 2.35 dB | 1.26 dB |
| measured profile (LIDAR and balloon) | 1.69 dB | 1.35 dB |
| `--wind`, station speed read as knots | **1.62 dB** | **1.25 dB** |
| `--wind`, station speed read as mph | **1.60 dB** | **1.25 dB** |
| `--wind`, knots, z₀ = 0.09 m | 1.71 dB | |
| `--wind --stability II` | 2.22 dB | |
| `--wind --stability IV` | 1.57 dB | |

- On the windy day the surface wind alone recovers all of what the measured profile gains over
  straight rays, and a little more: −0.07 dB against the measured profile, [−0.19, +0.05]. The
  gain is at the low-elevation microphones, 5.76 → 2.65 dB rms at 5–10°. On the calm day it
  does no harm.
- **The profile's discretization.** Refining the levels from 16 to 64 a decade moves the day-284
  results by 0.006 dB rms, 0.105 dB at most; the remainder is the stratified model's arrival
  search, whose samples follow the layers. The bottom layer, from the ground to ten roughness
  lengths, is linear in the built profile. Filling it with the log law (levels at 1.5, 2, 3, 5
  and 7 z₀) moves day 284 by 0.014 dB rms and 0.15 dB at most; dropping every level below 1 m
  moves it by 0.046 dB rms, 0.30 dB at most.
- **Stability classes are unvalidated.** Class II comes nearly back to straight rays: its
  inversion and shear deepen the upwind shadow, −3.3 dB median at 5–10°. An earlier trial of
  class II's exponent with a +3 K/100 m inversion, driven by the LIDAR's 46 ft wind, scored
  1.90 dB where the log law through that wind scored 1.72. One day, scored after the fact,
  cannot choose a class.
- **No veer with height.** An Ekman veer, turning the wind's direction with height, added
  nothing on the 2017 windy day, so `--wind` holds the direction at every height.

**`niceops_build_atmosphere` on the measured weather.** Built from the same measured columns
(LIDAR, balloon string, station air) with the matching settings, it reproduces the validation
study's own profiles for all 819 runs with measured tables (2017 and altitude-variation): 6e-14 K,
4e-14 % and 3e-14 kt at worst, pressure to 4e-13. Writing a profile to 3 decimals, as the study's
files do, moves 12 of 1796 mic-run SELs by more than 0.1 dB under `--refraction stratified`
(1.15 dB at most, at 5.6°): the stratified arrival search is that sensitive to perturbations of
1e-3.

`fit_z0` on the 46–197 ft gates, near-neutral runs only (|Ri| ≤ 0.05, 3 m/s at the top gate,
under 20° of veer):

- **Eglin**: z₀ = 0.34 m [0.10, 0.91] from 11 runs, rough to very rough by Wieringa's classes
  where forest stands within the gates' footprint; |Ri| ≤ 0.25 gives 0.42 [0.19, 0.80] from 72.
  Gates above 125 ft reach the search's bound, likely the forest's internal boundary layer. Carried
  down the log law, it predicts the station cups (2 m, near-neutral runs) to +0.03 m/s (rms 0.37)
  where the study's 0.1 ft gives +0.54 (0.59); the cups alone imply 0.47 m. With it, the R44's
  `stratified_plate` SEL rms at 0–5° falls 4.60 → 4.13 dB [−0.63, −0.31], and over all its
  mic-runs 1.77 → 1.66 (straight rays: 1.27).
- **Amedee, Salton Sea, Sweetwater**: no run meets the screen (the bulk Richardson number's
  median is +2.0 at Amedee), so the fit is refused. Relaxed screens give values the stability
  decides, not the surface: at Eglin the apparent z₀ runs from the search's lower bound (Ri < −1)
  to 2.3 m (weakly stable).

### Interpolated rays (`--ray_stride`)

`--ray_stride` against every state solved (`--ray_stride 1`), the same build. Stride 4 is the
default; the stratified results elsewhere in this file were made with every state solved, and
the default moves them by about what the tables below show. On the AS350 approach 291261 at
Amedee (46 flush plates, 4778 states, the day's measured profile, a history per microphone):

| stride | SEL, largest change | LAmax | one history sample, largest | samples over 0.1 dB |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0.0004 dB | 0.0007 dB | 0.41 dB | 4 |
| 4 | 0.0005 dB | 0.0016 dB | 1.15 dB | 21 |
| 8 | 0.0014 dB | 0.015 dB | 1.73 dB | |

of 219,788 samples. The largest sample changes fall where the solved history itself steps by
about 2 dB between consecutive states, near the end of the approach, which the interpolation
smooths. SEL on 1681 cells of the demo grid (41 x 41), largest change and 99th percentile:

| profile | stride 4 | stride 8 |
| --- | --- | --- |
| `example_windy.csv`, point receivers | 0.012 / 0.0005 dB | 0.004 / 0.0006 dB |
| `example_windy.csv`, flush plates | | 0.004 / 0.0004 dB |
| `example_windy.csv`, `--shadow_model resolved` | | 0.003 / 0.0006 dB |
| `example_inversion.csv` (a duct), point receivers | 0.057 / 0.003 dB | 0.27 / 0.002 dB |
| the duct, `--tunneling --skip_fill` | | 0.004 / 0.002 dB |

A windowed run equals its marched run byte for byte at any stride (`stratified_windowed_is_exact`).

### Refraction against flight data

**`--refraction stratified` on day 284** (`min2_no284`, plate, 729 mic-runs, SNR ≥ 10 dB,
records that miss the flyover left out), median / rms in dB, before and after the boundary law
was refitted for flush microphones
([Shadow-transition calibration](#shadow-transition-calibration)):

| | before the flush refit | current law |
| --- | --- | --- |
| upwind, peak elevation 5–10° (65) | +0.45 / 3.72 | +0.28 / 3.06 |
| upwind, 10–20° (80) | +0.99 / 2.24 | +0.71 / 1.87 |
| all | +0.09 / 1.85 | +0.09 / 1.70 |

The calm day 286 barely moves (1.46 → 1.45 dB rms overall): its winds cast few shadows.

- **The plate's own transition table**, with the offset K, takes `stratified`'s day-284 rms at
  5–10° from 4.5 to 3.6 dB.
- **Microphones just below the plane.** Before the direct-ray solve placed such microphones on
  the ground plane, a ray grazing under it read them 33 dB low on day 284.
- **The hovering end of the day-286 approaches** (sources under 100 ft and 20 kt).
  `stratified` was 6.5 dB over straight rays there on the measured frames, 63–250 Hz, at
  microphones 150–450 ft away, while the PE on those paths agreed with straight rays. In the
  shadow, the layered ray focused by the log law's lowest meters plus the whole ground term gave
  12 dB over free field at 63 Hz where the PE gives 4. With the near-source boundary and
  cramped-band rules, 63–630 Hz energy-mean error:

  | | without the rules | with them | straight rays |
  | --- | --- | --- | --- |
  | sources under 100 ft and 20 kt | +8.9 dB | +4.4 dB | +3.8 dB |
  | within 2,000 ft, 2–5° | +4.7 dB | +3.1 dB | +2.4 dB |
  | every other height, speed, range or elevation cell | | within 0.2 dB | |

  The 4 dB left at the hover is in straight rays too, so it is the source.
- **A calm floor on day 284.** A station-wind log law with z₀ = 0.1 ft, written with levels at
  0, 0.01, 0.03, 0.1, 0.3 and 1 ft (the first four calm), under plain ray theory predicted
  22–28 dB low at 2–5 low-elevation upwind flush microphones a run. Run 284202's microphone 3, at
  a 12° peak elevation, measured 78.3 dB and was predicted near 50 (78.6 without the floor). It
  took day 284's SEL rms from 1.63 to 4.18 dB. With the cap of 1.5 k h², the z₀ = 0.1 ft log law
  scores 1.63 dB, the same as z₀ = 0.003 ft, which has no floor; a z₀ = 1 ft log law, calm to
  1 ft, goes from 3.95 to 1.90. Run 284202's worst microphone is within 0.07 dB of its floorless
  level, with `--tunneling` and `--skip_fill` too, and every set of levels below 1 ft that was
  tried lands within 0.15 dB of the floorless run. The measured-profile runs move by at most
  0.0006 dB on day 284 and 0.0002 dB on day 286.
- **Day-286 frames in an elevated duct.** Where `stratified` predicts 10 dB or more under
  straight rays, it is 10–20 dB too quiet and straight rays are within 4 dB at 50–630 Hz. Those
  frames are far from each microphone's peak and do not move the SEL
  ([Micro-shadows and skip zones](#micro-shadows-and-skip-zones)).

**`--shadow_model resolved`** (plate receivers, measured profiles, paired with
`stratified_plate`):

| set | mic-runs | `stratified_plate` rms | resolved rms | Δrms [95% CI] |
| --- | --- | --- | --- | --- |
| R44 approaches, all | 585 | 1.84 | 1.53 | −0.31 [−0.42, −0.19] |
| R44 approaches, 0–5° | 45 | 4.74 (mean −3.89) | 3.55 (mean −2.88) | −1.19 [−1.55, −0.73] |
| R44 approaches, 5–10° | 118 | 2.09 | 1.70 | −0.38 [−0.63, −0.12] |
| B407 day 284, all | 729 | 1.69 | 1.78 | +0.09 [+0.02, +0.16] |
| B407 day 284, < 10° | 73 | 2.98 (mean +1.00) | 3.42 (mean +2.32) | +0.44 [+0.03, +0.75] |
| B407 day 286, all | 1111 | 1.35 | 1.36 | +0.00 [−0.00, +0.01] |
| B407 day 286, < 10° | 62 | 1.91 | 1.94 | +0.03 [−0.06, +0.11] |
| R66 approaches, all | 874 | 1.37 | 1.41 | +0.04 [+0.02, +0.05] |
| R66 approaches, < 10° | 99 | 1.99 | 2.13 | +0.15 [+0.10, +0.20] |
| R44 turns (held out), all | 4992 | 1.55 | 1.56 | +0.00 [+0.00, +0.00] |
| R44 turns, < 10° | 11 | 0.94 | 0.90 | −0.04 [−0.07, +0.00] |

- On the R44 it matches the per-run PE propagation tables (1.47 dB overall and 3.51 at 0–5°)
  with no table.
- It fails the windy day's do-no-harm gate (day 284 below 10°: rms ≤ 3.15, mean +1.00 ± 0.3,
  Δ upper bound ≤ +0.15; all ≤ 1.75). That is the strong-wind over-reach seen against the PE,
  and it is why the option stays opt-in.
- With per-microphone offsets removed, the R44 rms goes 1.13 → 1.02 (0–5°: 1.82 → 1.33). The
  approaches' offset of microphones 1, 2 and 4 as a group is not propagation and does not move.
- `--cramped_ground straight` scores the same on day 284 (3.44 below 10°). `--shadow_phase` 0.5,
  as a sensitivity only, gives 3.24 there (mean +2.00) and 1.58 on the R44 approaches (3.79 at
  0–5°). The boundary placement, not the cramped field, is what moves it.

### Turbulent scattering strength

At γT = 3e-5 (what the 2017 co-located pole and board pairs fit, in panam), the PE fills the
day-284 shadows to −3 dB re free field at a board, median, only 5–10 dB below the lit level. The
scattered-level fit behind `--scattering` reproduces the PE to 1.5 dB rms for flush receivers
(12 profiles, sources 50 and 120 m); the 1.2 m fit uses 6 profiles, sources 50 m.

The measurements want much less. Against the day-284 SEL, the upwind 5–10° bin (median / rms):

| shadow fill | median / rms |
| --- | --- |
| `--shadow_floor` (no scattering) | +0.28 / 3.06 dB |
| γT = 1e-6 | +0.49 / 3.04 |
| γT = 3e-6 | +0.95 / 3.02 |
| γT = 1e-5 | +1.55 / 3.16 |
| γT = 3e-5 | +2.96 / 4.00 |

Up to about 3e-6 is consistent with them, HARMONOISE's "moderate" range. The pair fit's 3e-5 is
ten times that and three times HARMONOISE's maximum, and probably carries more than turbulence
(roughness, the pole's height), or turbulence at the pole's height, where it is strongest. The
floor does as well on these data, so `--scattering` stays off by default.

### Ground effect

On the S-76D demo track with the B407 spectral database, SEL over the demo grid, dB:

| | mean | max |
| --- | --- | --- |
| spectral, no ground | 65.10 | 96.38 |
| + concrete (σ = 10000 kPa·s/m²) | 68.45 (**+3.35**) | 99.00 |
| + grass (σ = 225 kPa·s/m²) | 64.73 (**−0.36**) | 97.89 |

Hard ground is 3.7 dB louder than soft, and grass roughly cancels its own reflection gain against
the interference dip. The flow resistance spans five orders of magnitude, and the answer is
sensitive to it.

Delany–Bazley is poor for grassland. On the 2017 Amedee pole microphones, panam's ground-plane
study found variable porosity with σe ≈ 200 fitted best, with HARMONOISE γT ≈ 3e-5 and 10–15 mm
roughness.

### Terrain at Amedee

Along every track-to-microphone path within 3 km at Amedee, a mean plane departs from the ground
by at most 0.3 m. With the 3DEP grid, in the validation study's configuration (station air,
flush plates, the shipped B407 database for day 286 and `min2_no284` for day 284), with against
without `--terrain`:

| runs | microphone-runs | mean ΔSEL | largest \|ΔSEL\| |
| --- | --- | --- | --- |
| 286264, 286265 (approaches) | 92 | −0.003 dB | 0.026 dB |
| 284153, 284197 (steady) | 94 | −0.006 dB | 0.023 dB |
| all four | 186 | −0.004 dB (rms 0.007) | 0.026 dB |

Terrain does not explain the residuals at Amedee. The 1.2 m pole microphones move more: on
286264 and 284153, up to 0.11 dB with `--observers_on_terrain` and 0.19 dB with the surveyed z.
Their height above each plane follows the decimeter-scale relief and the grid's error at their
feet. For a sloped site, see [Sensitivity studies](#sensitivity-studies).

### LAmax, TA and EPNL against the 2017 measurements

**LAmax.** On the four detail runs (`station_plate` on 286270 and 286277, `stratified_plate` on
284205 and 284211, 188 microphones):

| In-tool LAFmax minus | Median, per run | 5th to 95th percentile | Largest difference |
| --- | --- | --- | --- |
| an independent evaluation of the same definition from `--band_history`'s reception times | — | — | 2e-6 dB (the file's printed digits) |
| the validation study's procedure without its 0.1 s blocks (emission time, 0.02 s) | −0.003 to −0.005 dB | −0.043 to +0.015 dB | 0.075 dB |
| the validation study's LAFmax (energy-averaged into 0.1 s blocks, then F) | +0.03 to +0.06 dB | −0.008 to +0.32 dB | 0.55 dB |

The last row is the 0.1 s blocking, which the study applies to put the prediction on the measured
frames' footing. Blocking the same reception-time signal into 0.1 s blocks moves it by +0.05 to
+0.07 dB at the median, while the clock alone (reception against emission time, both blocked)
moves it by −0.002 to +0.015. Within 0.05 dB of the blocked LAFmax are 44–74% of microphones per
run; within 0.08 dB of the unblocked procedure, all of them. Over every microphone-run of the
study's cases, the in-tool LAFmax is +0.033 to +0.046 dB from the blocked value at the median per
case (5th to 95th percentile −0.018 to +0.203 dB, 0.71 dB at worst, 53–62% within 0.05 dB).
Scored against the measurements the two differ by 0.02 dB rms or less. For scoring against
measured 0.1 s frames, keep the blocking; the in-tool value is the meter's reading of the
continuous prediction.

**Time weightings compared.** On the four detail runs, the instantaneous maximum of the 0.02 s
history is +0.72 to +0.93 dB above LAFmax at the median (95th percentiles +1.1 to +3.8 dB).
Against the study's LAFmax on the calm day it sits +0.95 dB high at the median (95th percentile
+2.29 dB). LAFmax is +0.52 to +0.92 dB above LASmax and +0.09 to +0.16 dB above the running 0.5 s
Leq, at the median.

**Time above and C-SEL.** On the detail runs, with instantaneous levels, TA65 is 11–36 s and TA75
4–10 s at the median microphone. C-SEL is 8.8–9.5 dB above the A-weighted SEL at the median
microphone, and on run 286270 40% of the C-weighted energy (median microphone) lies below 50 Hz,
so the default 20 Hz lower band limit leaves out a database's 10–16 Hz bands.

**EPNL.** On the four detail runs, over the microphones whose SEL passes the study's gates,
predicted EPNL against EPNL of the measured 0.5 s frames, both masked. The measured frames stop at
5 kHz, so for a like-for-like comparison the predicted frames' 6.3, 8 and 10 kHz bands were
replaced by the 5 kHz level less 5, 10 and 15 dB (a roll-off that triggers no tone correction).
Predicted minus measured, dB, with the full 24-band `--epnl_file` value beside it:

| day, model | mic-runs | EPNL median | EPNL rms | PNLTM median | PNLTM rms | `--epnl_file` EPNL median, rms |
| --- | --- | --- | --- | --- | --- | --- |
| 286 calm, `station_plate` | 88 | −0.33 | 1.22 | −0.89 | 2.00 | −0.23, 1.21 |
| 286 calm, `stratified_plate` | 88 | −0.09 | 1.26 | −0.75 | 1.96 | +0.01, 1.28 |
| 284 windy, `station_plate` | 87 | −0.22 | 2.70 | −0.85 | 3.07 | −0.11, 2.70 |
| 284 windy, `stratified_plate` | 87 | −0.88 | 1.90 | −1.36 | 2.91 | −0.78, 1.91 |

- **Without masking** (panam's procedure on both sides) the EPNL median / rms read −0.25 / 1.23,
  +0.02 / 1.33, −0.22 / 2.85 and −0.80 / 1.87. That reproduces an independent computation from
  per-state band histories framed by the states' reception times (1.24 and 1.88 dB rms for the
  calm `station_plate` and windy `stratified_plate` rows), so the two framings agree.
- **Masking** leaves the measured EPNL where it was (0.01 dB rms). Unmasked, the full 24-band
  predicted EPNL read 0.04–0.19 dB (median; 0.24–0.41 rms) above the 5 kHz-limited one; masked
  it reads the same (0.00–0.01 median, 0.11–0.21 rms), so the excess was inaudible tones.
- **About half the PNLTM bias is the tone correction:** the measured spectra carry a median
  1.2–1.3 dB of it at PNLTM, against 0.8–0.9 dB predicted, consistent with the blend of spheres
  and grid cells smoothing the tones. EPNL, whose duration correction integrates every record,
  is biased less.
- **Masking on predicted spectra.** On the S-76D demo track over a 100 m grid (961 observers,
  B407), the 10 kHz band's median level at PNLTM is −194 dB, and a near-silent database band
  reaches −1000 dB. Unmasked, 615 of the 961 observers took their PNLTM tone correction from a
  band below its noy threshold (10 kHz rises over 8 kHz for 81 of them), and EPNL read a median
  1.0 dB (up to 5.4 dB) high. Masking moved the PNLTM tone correction's median from 2.25 to
  1.20 dB and the count of 3 dB or more from 427 to 13, lowered EPNL by a median 1.0 dB, and
  raised it at one observer, by 0.2 dB.
- **The frames.** Rebuilt independently from a `--band_history` run's files, the frames agree
  with `--band_frames` to 6e-5 dB, and panam's EPNL on them with `--epnl_file` to 5e-8 dB.

### Day-night window choice

Integrating each operation over its whole history instead of its own 10 dB-down window: on the
four detail runs the whole history is +0.29 to +0.56 dB above the window at the median per
operation, and up to +1.4 dB at a single microphone. For a day-night level meant to sit beside an
AAM or NORAH2 contour, use `--full_history`; to match a measured SEL's convention, leave it off.

### A loudspeaker among buildings: NLR Drone DigiCity

Knepper, Admiraal and Lammers (ERF 2025 paper 518) measured a B&K OmniPower 4296 on a 2 m tripod
in NLR's mock city of 24 shipping containers on a concrete slab, at microphones 1.2 m up. They
compared it with NLR's CNOSSOS-EU ray model at four microphones. Their model agrees at the three the
source reaches directly and by reflection (2, 8 and 10), and is 10–15 dB low at microphone 7,
behind the central bar of containers, which the source reaches only by diffraction. The data here
are the paper's: the layout read off its plan figure (heights from its photographs) and the spectra
digitized from its plots. The comparison is therefore not blind, and the geometry is a
reconstruction.

The source's level at 1 m is calibrated per band on microphones 2, 8 and 10 with
`niceops_point_field` (walls absorbing 0.05, the floor rigid, the default ray orders, chosen before
any comparison). The three agree on it to 2.5 dB on average. Microphone 7, 160 Hz to 5 kHz,
predicted minus measured:

| Model | Mean | Rms | Range |
| --- | --- | --- | --- |
| `niceops_point_field` | −4.1 dB | 5.2 dB | −10.4 to +1.3 dB |
| NLR's CNOSSOS-EU model (as published) | −12.9 dB | 13.3 dB | −19.9 to −7.3 dB |

Below 1 kHz the prediction is within 3.4 dB in every band. It is 7.8–10.4 dB short
at 1.25–1.6 kHz and 5 kHz. Most of what reaches microphone 7 is diffracted once at a vertical
corner after one or two reflections off a stack's wall: the strongest ray, a third of the energy
with its floor partner, reflects off that wall and diffracts at the south-east stack's corner, and
the next at the bar's own end. Without that wall's stack microphone 7 drops by 8–22 dB. Double diffraction over the bar or around its end is 8–15 dB weaker,
and the cascade gives the same result here as the shared-face form. Walls rigid or absorbing 0.1, a
floor absorbing 0.1, the source 0.2 m higher, or calibration on fewer microphones each move the
mean by 0.4 dB at most. Allowing only one reflection either side of an edge loses 1.9 dB, and three
gain nothing. Over 60 reconstructions with the figure's scale within 3% and every container group
and microphone within 0.28 m, the mean ranges from −5.2 to −2.9 dB (10–90%).

The shortfall at 1–2 kHz is larger than that spread. The corrugated walls scatter: with the
corrugation's own s(f), computed from its profile and not fitted, the error comes to −2.3 dB mean,
2.9 rms, and the 1.25 kHz and 5 kHz bands to within 3 dB
([Scattering walls (point source)](#scattering-walls-point-source)). A run takes the same field with
`--facade_scattering_spectrum`. The source is not omnidirectional above about 500 Hz and is
calibrated in other directions, which has not been measured here.

## Sensitivity studies

**Terrain slope.** A 1–3° slope along the path changes the grazing angle by the slope. Through
the plate table and the B407 spectrum, that is about 1 dB at 5° peak elevation and 0.3–0.7 dB at
10° for flush plates, and tenths of a dB at 1.2 m. These come from a standalone model of the
plate and pole at the changed elevation, not from a NICEOPS run (a positive slope raises the
elevation over the plane):

| peak elevation | plate, slope −3° / −1° / +1° / +3° | 1.2 m pole, 500 m, the same slopes |
| --- | --- | --- |
| 5° | −5.0 / −1.05 / +0.77 / +1.83 dB | −0.52 / +0.16 / 0.00 / +0.23 dB |
| 7.5° | −2.1 / −0.53 / +0.42 / +1.06 | −0.11 / −0.12 / +0.15 / +0.56 |
| 10° | −1.2 / −0.32 / +0.27 / +0.69 | −0.49 / −0.20 / +0.23 / +0.50 |
| 20° | −0.29 / −0.09 / +0.08 / +0.21 | −0.43 / −0.04 / −0.03 / +0.04 |

**A stratified against a homogeneous atmosphere.** On the S-76D demo track with the B407
database, `atmospheres/isa.csv` moves 90% of the footprint by −0.10 to −0.005 dB (mean
−0.04 dB, extremes −0.57 and +1.11 dB) against a homogeneous run at ISA sea-level conditions.

## Validation required before replacing a source database

Keep complete runs out of training. Score the held-out measurements by source height, speed,
emission angle, range, frequency and microphone installation, alongside per-band SEL and total
exposure. Report coverage, per-run variation and clustered uncertainty: pooled energy means can
hide opposing errors. Repeat SNR-gate sensitivity checks. A propagation improvement on frames
selected by the old model's own shadows must also pass the other geometry regimes. A
range-dependent source fit would confound descent state drift with path loss.

## Known limitations

- **Elevated ducts whose returning rays turn just above the ground.** Ray theory misses both the
  tunneled field and the skip zone ([Micro-shadows and skip zones](#micro-shadows-and-skip-zones)).
  The opt-in Airy corrections recover part of it below 800 Hz, and none of it above.
- **Arrivals at caustics.** Two close arrivals within one pair of launch samples can still be
  missed, mostly under a strong inversion ("Arrivals missed at caustics" in
  [Refraction models against the PE](#refraction-models-against-the-pe)).
- **A profile whose lowest level stands above the ground**, such as one starting at a 46 ft
  LIDAR gate, keeps a constant floor that thick (14 m), which the calm-floor cap follows for
  kilometers. Give a profile a level at the ground.
- **The shadow boundary law** is several dB rms off the PE at the shadow edge, and no smooth law
  found does better ([Shadow-transition calibration](#shadow-transition-calibration)).
- **`--shadow_model resolved`** over-reaches in strong winds with a 1.2 m receiver, and fails the
  windy day's do-no-harm gate.
- **Stability classes** (`--stability`) are unvalidated, and the station wind's unit and height
  for `--wind` are a reading of the data report, not a measurement.
- **The PE's still-air deep shadow** (40–80 dB down) is not physical, and deep in shadows neither
  ray model is right: the floor holds them 25–35 dB above it.
- **Turbulent scattering** fitted to pole and board pairs is ten times what the SEL supports.
- **EPNL** is compared only to 5 kHz, where the measured frames stop, with no pseudo-tone
  exclusion; predicted spectra carry less tone correction than measured.
- **Terrain** is a plane fitted over the whole path, screened only with `--terrain_screening`.
  Behind the one measured hill `cnossos` is right to a few decibels below about 200 Hz and too
  deep above it, by more than 30 dB at 1–2 kHz; `kirchhoff` is within 4.5 dB to 500 Hz and 7–9 dB
  too deep at 1–2 kHz ([Terrain screening](#terrain-screening)).
- **Ground classes** treat each microphone's ground as a uniform half-space. At grazing a
  built-up patch measured like ground of about 1000–2000 kPa·s/m², not concrete, and the ground
  next to the microphone and along the path decides ([Ground classes](#ground-classes)).
- **The nose.** No shipped database records a heading-filed azimuth, and a WGS84 track's
  heading has to be given as a column.
- **A point source among buildings** diffracts at three edges at most but over roofs, where chains
  take up to eight; with both ends above a shared face, double diffraction is the cascade, off by
  up to 10 dB near the plane; blocks thinner than a wavelength or so are 2-6 dB low below
  500 Hz; and by UTD an edge's end, where its transition meets a face's plane through the vertex,
  is up to 3.7 dB off the edge-source integral, no better than with the edge unending
  ([A point source among buildings](#a-point-source-among-buildings)). A pitched roof's
  slopes reflect with walls, and a level roof or a slope at one bare end of a twice-diffracted ray,
  but not with walls on that ray, nor on rays diffracted three times or more; over a roof of nearly flat wedges whose pairs share the faces
  the ray runs along (a mansard's) the cascade is 3-4 dB rms off the exact field; the FDTD's
  staircase makes it no reference for a sloped roof below 1.5 mm cells at EDGE's scale
  ([Pitched roofs](#pitched-roofs)).
- **Source coverage.** Every aircraft over-predicts at low elevation; the B407 hover between 3°
  and 17° below the horizon, and under the aircraft, is an interpolation.
