# NICE-OPS fidelity

How far each part of NICE-OPS can be trusted, and the evidence for it. The equations are in
[THEORY.md](THEORY.md), how to run each option is in the [manual](MANUAL.md), and changes that
moved results are in [HISTORY.md](HISTORY.md).

## What is measured here

This page keeps three kinds of evidence apart:

- **Regression agreement.** A run reproduces a stored reference. Passing a golden comparison
  establishes stable behavior. It does not establish physical accuracy.
- **Numerical accuracy against a reference solver.** A model is set against an exact solution,
  or against a full-wave parabolic equation (PE). A PE reference is not a measurement.
- **Accuracy against measurements.** Predictions are scored against the 2017 Noise Abatement
  flight test ([Data sets](#data-sets)).

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
| NICE-OPS (travel capped at 1.5 k h²) | 186 m | 186 m | 187 m | 233 m | 544 m |

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

All measured comparisons use the 2017 Noise Abatement flight test (data report
NASA/TM-2019-220264).

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

NICE-OPS minus the measured 10 dB-down SEL over each aircraft's approaches, station atmosphere,
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
- **Coverage.** On the demo grid, the B407 database's levels come 1.2% from unmeasured cells at
  the median observer and 27% at the 95th percentile. The builder's SNR gate drops the
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

## Sensitivity studies

**Terrain slope.** A 1–3° slope along the path changes the grazing angle by the slope. Through
the plate table and the B407 spectrum, that is about 1 dB at 5° peak elevation and 0.3–0.7 dB at
10° for flush plates, and tenths of a dB at 1.2 m. These come from a standalone model of the
plate and pole at the changed elevation, not from a NICE-OPS run (a positive slope raises the
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
- **Terrain** is a plane fitted over the whole path, with no screening.
- **The nose.** No shipped database records a heading-filed azimuth, and WGS84 tracks carry no
  heading.
- **Source coverage.** Every aircraft over-predicts at low elevation; the B407 hover between 3°
  and 17° below the horizon, and under the aircraft, is an interpolation.
