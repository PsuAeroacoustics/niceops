# NICE-OPS manual

This manual describes what NICE-OPS computes, the models behind each option, and how each
model has been checked. To build the code and run a first case, start with the
[README](../README.md). Developer notes, including the shape of the hot path and the pitfalls
to know before changing the numerics, are in [`CLAUDE.md`](../CLAUDE.md).

## Contents

- [How a prediction is made](#how-a-prediction-is-made)
  - [The operating condition](#the-operating-condition)
  - [Blending the database's spheres](#blending-the-databases-spheres)
  - [Other ways to choose the spheres: `--blend`](#other-ways-to-choose-the-spheres---blend)
  - [The load-factor law below 1 g](#the-load-factor-law-below-1-g)
  - [Reading the sphere](#reading-the-sphere)
  - [Propagation to the microphone](#propagation-to-the-microphone)
  - [Exposure](#exposure)
  - [EPNL and reception-time frames](#epnl-and-reception-time-frames)
  - [Where the code is](#where-the-code-is)
- [Running NICE-OPS](#running-nice-ops)
  - [Command-line options](#command-line-options)
  - [Units](#units)
  - [Output files](#output-files)
  - [Several trajectories, and the day-night level](#several-trajectories-and-the-day-night-level)
  - [Maximum level, time above and C-weighted SEL](#maximum-level-time-above-and-c-weighted-sel)
- [Inputs](#inputs)
  - [Input files](#input-files)
  - [Geodetic coordinates](#geodetic-coordinates)
- [Noise databases](#noise-databases)
  - [Spectra](#spectra)
  - [Where the 2017 databases come from](#where-the-2017-databases-come-from)
  - [Speed reference](#speed-reference)
  - [The sphere's nose](#the-spheres-nose)
- [Propagation](#propagation)
  - [Broadband and spectral paths](#broadband-and-spectral-paths)
  - [Stratified atmosphere](#stratified-atmosphere)
  - [A profile from the surface wind](#a-profile-from-the-surface-wind)
  - [Ground effect](#ground-effect)
  - [Receivers](#receivers)
  - [Ground-plane microphones](#ground-plane-microphones)
  - [Terrain](#terrain)
  - [Refraction](#refraction)
  - [Refraction through the profile's layers](#refraction-through-the-profiles-layers)
  - [Micro-shadows and skip zones](#micro-shadows-and-skip-zones)
  - [Propagation loss tables](#propagation-loss-tables)
- [Reference solutions](#reference-solutions)
  - [Checking the refraction models against a parabolic equation](#checking-the-refraction-models-against-a-parabolic-equation)
  - [Turbulence, and the energy it scatters into shadows](#turbulence-and-the-energy-it-scatters-into-shadows)
  - [The Attenborough et al. benchmark](#the-attenborough-et-al-benchmark)
- [Performance](#performance)
  - [How much gets evaluated](#how-much-gets-evaluated)
  - [`--db_cutoff`](#--db_cutoff)
  - [Threads](#threads)
  - [Cost of the spectral and ground paths](#cost-of-the-spectral-and-ground-paths)
- [Plotting](#plotting)
- [Changes that affect results](#changes-that-affect-results)
  - [The 4x EAA multiplier](#the-4x-eaa-multiplier)

## How a prediction is made

NICE-OPS is a source–path–receiver model. The source is a database of noise spheres around
the aircraft, one per operating condition. The path carries the level from the aircraft to
each microphone. The receiver integrates the level history it sees into an exposure metric.
The steps below run for every trajectory state and every observer: on the demo case, 6,551
states by 14,400 observers.

### The operating condition

Each trajectory state is reduced to three non-dimensional parameters (`src/trajectory.cpp`):

| Parameter | Definition |
| --- | --- |
| Advance ratio μ | \|v\| / V<sub>tip</sub> |
| Flight path angle γ | asin(v · n̂ / \|v\|), the velocity's angle to the tip-path plane (the default, `--qsam`); or asin(v<sub>z</sub> / \|v\|), its angle to the horizontal (`--no-qsam`) |
| Thrust coefficient C<sub>T</sub> | n W / (ρ V<sub>tip</sub>² π R²) |

The load factor **n** is the aircraft's acceleration less gravity, in g; n is its magnitude
and n̂ its direction, which stands in for the rotor's thrust axis, the tip-path-plane normal.
The tip speed V<sub>tip</sub>, rotor radius R, weight W and air density ρ come from the
vehicle file. The velocity v is taken over the ground unless `--speed_reference air` takes it
through the air ([Speed reference](#speed-reference)). An AC state track carries its own
accelerations. For a WGS84 track they are differentiated from a smoothing spline
([Geodetic coordinates](#geodetic-coordinates)).

### Blending the database's spheres

The database's operating conditions are normalized into the unit cube and triangulated: a
Delaunay triangulation on exact predicates (`src/delaunay.cpp`, `src/predicates.cpp`). The
state's condition falls in one tetrahedron, and the spheres at its four corners are blended
with barycentric weights, in decibels (`src/interpolator.cpp`). `--blend` chooses another
way to pick the spheres (below). A `fixed_load_factor`
database stores 1 g spheres only. It is triangulated over (μ, γ) alone, and a load factor
above 1 g adds 20 lg(C<sub>T</sub> / C<sub>T,1g</sub>) dB, up to `--max_load_factor`
(`--load_factor_correction two_sided` takes the law below 1 g too). The
blend depends only on the operating condition, so it is formed once per state and shared by
every observer.

A condition outside the range of a database axis is clamped to that axis's end before
blending. Each run logs how many states were clamped on each axis, and how many took the
load-factor correction, and warns if any were clamped. The counts are about the axes'
ranges only: a condition inside every range can still lie outside the convex hull of the
database's conditions, and neither says which directions on the sphere were measured.

### Other ways to choose the spheres: `--blend`

The blend above is the default, `--blend barycentric`. Two others exist for comparing
NICE-OPS with other codes and with regulatory practice. Neither predicts better (see the
scores below); they answer the question of how much of a difference between two models is
the choice of sphere.

**`--blend nearest`** takes the one stored condition nearest the state, whole, at every
state. Nearest means in the space the triangulation is built in: the normalized (μ, γ, C<sub>T</sub>),
or (μ, γ) on a `fixed_load_factor` database, with the state clamped into it first. In raw
units the path angle, tens of degrees wide, would decide every distance against an advance
ratio a few tenths wide. A tie goes to the database's earlier sphere. ECAC Doc 32 and
models built on it assign one hemisphere to each constant segment of a procedure. A
trajectory tool has no segments, and taking the nearest sphere at every state is the
closest it can come; on a track that holds one condition it is the same thing.

**`--blend aam`** is the rotorcraft sphere search of AAM v3.1 (*Advanced Acoustic Model
Technical Reference and User's Guide*, Volpe, May 2021, §2.5.1.1, pp. 31–34):

1. The spheres are sorted into rows of one flight path angle, each angle rounded to the
   nearest 0.5° first (AAM has done so since RNM 7, p. 31; halves round away from zero, as
   Fortran's `NINT`). A row's angle is its rounded one; the state's is not rounded.
2. The rows bracketing the state's path angle are interpolated between, if they are no
   further apart than 3° on approach (the state descending) or 5° on departure (climbing or
   level), Table 2-7. A wider pair is not interpolated: the nearer row is taken alone. So is
   the end row for a state past either end. Path angle may disqualify spheres (Table 2-8).
3. Within each row, the speeds bracketing the state's are interpolated between if they are
   no more than 20 kt apart. Otherwise the nearest speed in the row is taken, however far:
   speed "may NOT disqualify" a sphere. A state exactly halfway takes the faster sphere,
   "generally a conservative result".
4. The two rows' speed blends are then blended in path angle: up to four spheres, linearly
   in dB, speed first and then path angle (Figure 2-20).

A pair exactly as wide as a limit is still interpolated (the limits are maxima, "dFPA deg
max"). The document leaves three things open, and NICE-OPS settles them this way:

- **Ties in path angle.** Table 2-8 says approaches "(FPA>0.0)" take the higher angle and
  departures "(FPA<0.0)" the lower, but the same document has positive angles climbing
  (Table 2-5, Table 4-37). NICE-OPS follows the labels: an approach takes the higher (less
  steep) angle, a departure the lower, toward level flight in both.
- **The nearest-sphere fallback.** "For vehicles with sparsely populated acoustic spheres
  flight conditions, the prior RNM methodology of selecting the closest matching flight
  condition available sphere without interpolation is utilized" (p. 31). With the rules
  above that is what a state past every gate gets: the nearest row in path angle, then the
  nearest speed in it. It is not the sphere nearest in a combined distance: path angle comes
  first.
- **Two spheres in one row at one speed** (their angles round alike). AAM stops on two
  spheres at one profile; NICE-OPS takes the one whose own angle is nearer the state.

Speeds are μ V<sub>tip</sub>, and the 20 kt gate is converted to advance ratio with the
database's tip speed, the one its spheres' μ were formed with, or the vehicle file's for a
database that records none. The path angle is the one the run uses, γ relative to the
tip-path plane under `--qsam` (the default; AAM's Q-SAM options do the same) and to the
horizontal under `--no-qsam` (AAM without them). AAM has no thrust axis. On a
`fixed_load_factor` database the load-factor law carries C<sub>T</sub> exactly as it does
under the barycentric blend. On a database whose load factors are spheres, each state takes
the thrust layer nearest its own C<sub>T</sub> (a tie takes the higher) and searches it
alone; the run logs how many layers there are, and warns if any holds a single sphere (a
database whose C<sub>T</sub> is scattered sphere by sphere rather than layered would leave
the search nothing to bracket). AAM's FRAME option, with load factor as a
third searched axis, its "duration adjustment" and its nacelle-angle axis are not
reproduced. The run logs how many states blended four, three, two and one sphere.

Two consequences of the rules are worth knowing. The gate changes at level flight, so a
state between rows 2° either side of level interpolates while it climbs and takes the nearer
row alone once it descends: AAM's output steps there, as it steps wherever a pair crosses a
gate. And the shipped databases pad their envelope with copies of the edge spheres at γ = −24°
and 35° (in `B407_fixed_load_factor.nod` all 60 are copies), so that the triangulation
covers steep states. Under `aam` those are ordinary rows, far outside every gate from the
measured ones: a state more than halfway from the outermost measured row to a padding row
takes the padding row alone, whose spheres repeat edge spheres.

The 2017 databases are a scattered cloud of conditions, not a grid, which is the case the
bracketing search handles worst. On the day-286 B407 approaches (all 114,720 states of the
25 runs) AAM interpolated in both axes at 9.4% of states, took three spheres at 35.8%, two
at 34.3% and one at 20.5%. The barycentric blend does not gate at all: its spheres span a
median 13.5° of path angle on those approaches, and 76% of states take a sphere further
than AAM's gates allow (more than 5° or 20 kt away).

How they score against the measurements, paired against the default on the paper's
held-out cases (run-level bootstrap 95% intervals; SEL and LAFmax error rms over every
scored microphone):

| Case | Metric | barycentric rms | nearest Δrms | aam Δrms |
| --- | --- | --- | --- | --- |
| Day 286 approaches, `station_plate`, B407 shipped (25 runs, 1111 mics) | SEL | 1.26 dB | +0.47 [+0.37, +0.57] | +0.24 [+0.02, +0.48] |
| | LAFmax | 1.94 dB | +0.45 [+0.32, +0.60] | +0.34 [+0.18, +0.50] |
| Day 284 steady runs, `stratified_plate`, held-out `min2_no284` (16 runs, 729 mics) | SEL | 1.69 dB | +0.14 [+0.03, +0.23] | +0.42 [+0.17, +0.66] |
| | LAFmax | 2.99 dB | −0.38 [−0.61, −0.19] | −0.15 [−0.43, +0.13] |

The barycentric blend is the best of the three on SEL on both days, which is why it stays
the default. Nearest is worse than it on every one of the 25 approaches. The one exception
is LAFmax on day 284, where the barycentric prediction runs 1.6 dB low at the median and
nearest takes 0.4 dB rms off: a blend in dB of spheres whose directivity peaks lie in
different directions flattens the peak, and a single sphere keeps it. The SEL, which
integrates over the directions, does not show that. Neither alternative changes the
default's SEL by much on average (nearest +0.6 dB on day 286, aam +0.03), but individual
microphones move by up to 3 dB (nearest) and 6 dB (aam), most of it below 20° of elevation,
where the event is long and many conditions are passed through. The runs are the
validation harness's fidelity-paper cases (`paper/fidelity_study.py`) with `--blend` added;
the barycentric runs reproduce the paper's pinned results (6bbd6bb) to the printed digit.

Checked by: unit tests of each rule on small synthetic condition sets
(`tests/test_condition_selection.cpp`: bracketing, both path-angle gates, a pair exactly as
wide as a gate, ties in path angle and speed, the fallback to the nearest row and then the
nearest speed, rounding, and duplicates within a row); database-level tests that nearest
searches the normalized space, that aam converts its gate with the tip speed, refuses a run
with none and picks the nearest thrust layer, and that the direction bounds hold under every
policy (`tests/test_sphere_database.cpp`, tag `[blend]`); and CTests that naming the
default reproduces the golden footprint, that nearest and aam each move it and differ from
each other, and that the windowed peak search stays exact against a march under each.

### The load-factor law below 1 g

A `fixed_load_factor` database stores only the 1 g condition, and the load factor enters
through the analytic law 20 lg(C<sub>T</sub> / C<sub>T,1g</sub>) dB. By default the law is
applied above 1 g only (`--load_factor_correction one_sided`), as panam applies it: a state
below 1 g, which an approach's flare and pushovers produce, reads the 1 g sphere unscaled.
`--load_factor_correction two_sided` applies the same law below 1 g too, down to
`--min_load_factor` (default 0.5 g), below which the correction is held: the law is a
logarithm, and near zero g (or below it, where C<sub>T</sub> is negative) it would take the
level tens of dB down or fail. The run warns with the number of states held at either end.
On a database whose load factors are spheres there is no law to apply, and two-sided is
refused.

The evidence for the law below 1 g is modest, and mixed:

- Normalizing every sphere to the reference C<sub>T</sub> and predicting each day-286 approach
  at its own C<sub>T</sub> with the full law took 0.047 dB off the SEL rms, interval
  [0.012, 0.084] (an experiment on the database's spheres, not this option).
- Run directly, the option takes 0.028 dB off the SEL rms on day 286, interval
  [0.007, 0.048], better on 72% of the runs, and 0.018 dB on day 284, interval
  [−0.017, +0.052]. LAFmax is 0.05 dB worse on day 286 [+0.03, +0.08] and 0.02 dB worse on
  day 284 [+0.00, +0.04]. The level moves by 0.2 dB rms, at most 1.1 dB.
- On repeated steady runs of one flight condition (169 pairs of the B407 spheres, whose
  C<sub>T</sub> differs by weight and air density: 0.55 dB sd in 20 lg C<sub>T</sub>), the
  level difference follows the C<sub>T</sub> difference with a slope of 0.34 ± 0.17 per dB of
  the law, fitted together with speed: about a third of the law. The law's full slope is
  not confirmed run by run.

So it is offered, and the default stays one-sided. Checked by: unit tests that two-sided
follows the law on both sides of 1 g, holds at the floor through zero and negative load
factors, reaches an exported sphere, and that one-sided is unchanged
(`tests/test_sphere_database.cpp`); a CTest that it moves the footprint on the demo track,
whose load factor crosses 1 g; the windowed-peak-search exactness CTest, run two-sided; and
the refusals on the message.

### Reading the sphere

Each sphere is fixed to the rotor. Its axes are the tip-path-plane normal, the nose and
starboard. The nose is the velocity projected into the tip-path plane. At hover it is the
heading instead: below 3 kt the heading alone, blending into the velocity by 8 kt.
`--nose heading` takes the heading at every speed, which in a crosswind is a crab angle off the
track ([The sphere's nose](#the-spheres-nose)). The
direction from the source to the observer gives two angles (`lookup_angles` in
`src/sphere_database.cpp`):

- θ, from the nose: 0° straight ahead, 180° straight behind.
- φ, around the nose axis: 0° straight down along the tip-path-plane normal, +90° to
  starboard and −90° to port in the tip-path plane, ±180° straight up.

The level in that direction is interpolated bilinearly between the sphere's grid nodes,
periodically in φ. Unless every sphere in the database stores its upper half, a direction
above the tip-path plane reads its mirror image below.

### Propagation to the microphone

The sphere gives the level at its surface, radius R<sub>s</sub>. It is carried over the range
r to the microphone, with r held at no less than R<sub>s</sub>. There are two paths
([Broadband and spectral paths](#broadband-and-spectral-paths)):

- **Broadband**, for a database without spectra or under `--no-spectral`:
  L = L<sub>sphere</sub> − 20 lg(r / R<sub>s</sub>) − EAA (r − R<sub>s</sub>) / 1000 m. EAA
  is the database's A-weighted attenuation over the 1000 m reference distance it was built at.
- **Spectral**, the default for a database with spectra. Each band is spread the same way and
  absorbed at its own ISO 9613-1 coefficient. The ground effect, the receiver, refraction and
  propagation tables act on each band. The bands are then A-weighted, limited to
  `--band_limits` and summed in energy. The ground is the plane z = 0 unless `--terrain` gives
  each path its own mean ground plane ([Terrain](#terrain)).

### Exposure

The level history is sampled at the trajectory's emission times, but a microphone
accumulates exposure in reception time. There each sample lasts T (1 + ṙ / c), compressed
while the source closes and stretched while it recedes. The sound exposure level is

    SEL = 10 lg( T Σ (1 + ṙ_i / c) 10^(L_i / 10) )

with T the trajectory's own sample interval (not `-s`). The sum runs over the 10 dB-down
window: from the first to the last sample within 10 dB of the peak, dips between them
included, which is the duration the aircraft noise metrics use (14 CFR 36 A36.4.5).
`--full_history` sums over the whole history instead. Both are in `src/util/sel.hpp`.

Across a run, `-D`, `--lden_file` and `--cnel_file` sum every trajectory's exposure into the
day-night, day-evening-night and community noise equivalent levels, each with its own penalties
on samples heard in the evening or at night
([Several trajectories, and the day-night level](#several-trajectories-and-the-day-night-level)).
`--lamax_file`, `--nax_file` and `--ta_file` give each observer's time-weighted maximum level
(LAFmax by default), the number of operations that reach a threshold, and the time spent above
one; `--csel_file` gives the C-weighted SEL
([Maximum level, time above and C-weighted SEL](#maximum-level-time-above-and-c-weighted-sel)).
`--epnl_file` gives each observer's effective perceived noise level, from one-third-octave
levels in 0.5 s reception-time frames, and `--band_frames` writes those frames themselves
(next section).

### EPNL and reception-time frames

`--epnl_file FILE` writes each observer's EPNL, the certification metric of 14 CFR 36
Appendix A, section A36.4 (identical to ICAO Annex 16 Vol. I, Appendix 2), and `--band_frames
FILE` the one-third-octave band levels it is computed from, for every observer at once.
Both need the spectral path. EPNL is a per-operation metric, so with several trajectories each
gets its own file: one name given for several tracks becomes `<track>_EPNL.csv` and
`<track>_band_frames.nc` beside it, as `-S` does.

**The records.** The regulation computes EPNL from one-third-octave spectra averaged linearly
over 0.5 s records (A36.3). NICE-OPS evaluates one spectrum per trajectory state, in emission
time, every 0.02 s on the shipped tracks. The state emitted at *t* is heard at *t* + *r*/*c*,
and its emission step T arrives spread over T (1 + ṙ / c) of reception time — the same weight the
SEL integral gives it ([Exposure](#exposure)). Each state's band energy is spread evenly over
that stretch, split between the frames it overlaps in proportion, and summed; a frame's level is
the energy over the frame's 0.5 s:

    L_frame = 10 lg( Σ_i 10^(L_i/10) × overlap_i / 0.5 s )

So the frames are linear-average records in reception time, and together they hold exactly the
whole history's exposure in every band: A-weighted and summed, they give the `--full_history`
SEL. Frame *k* covers [0.5 *k*, 0.5 (*k* + 1)) s on the trajectory's own clock, so every
observer's frames share one time axis. The first and last frames an observer hears are partly
silent and are averaged over the whole 0.5 s regardless, as a record of them would be. *r* is
the straight-line range, or, under `--refraction`, *c* times the ray's travel time. A state
culled by `--db_cutoff` carries no energy, as in the SEL. `-a` is refused: the ambient is an
A-weighted level with no spectrum to add to the bands.

**The bands.** EPNL is defined over the 24 bands from 50 Hz to 10 kHz. Each is matched to the
database band whose center is within 5 percent of the nominal one (nominal and exact centers
differ by at most 1.6 percent; neighboring bands are 26 percent apart). A database lacking any of
the 24, or `--band_limits` leaving any out, is refused with the bands named: a missing band is
not a quiet one, and filling it would move both the noy sum and the tone correction's
background. Bands outside 50 Hz–10 kHz are not part of EPNL and are ignored by it (they are
still in `--band_frames`). The shipped spectral databases (`*_fixed_load_factor.nod`) carry
10 Hz–10 kHz, so all 24 are there. The levels are the database's spectra at the microphone, with
spreading, absorption, the ground or plate and any shadow loss applied; no frequency shift is
applied (the spheres keep the measured, Doppler-shifted spectra). The microphone is whatever
receiver the run models: a flush microphone on a plate (`--plate_table`) carries the plate's
pressure doubling, as the 2017 measurements do, where 14 CFR 36 prescribes a microphone 1.2 m up
(A36.3); the default point receiver of `--ground_effect`, 1.2 m up, is that one.

**The metric**, per record *k* and band *i*, with SPL(*i*, *k*) the frame levels:

1. Noys *n*(*i*, *k*) from Table A36-3's formulas (A36.4.7.3): four regions bounded by SPL(a),
   SPL(b), SPL(e) and SPL(d), and 0 noy below SPL(d).
2. PNL(*k*) = 40 + (10 / lg 2) lg *N*, *N* = 0.85 *n*<sub>max</sub> + 0.15 Σ *n* (A36.4.2).
3. The tone correction C(*k*), the ten steps of A36.4.3.1: slopes from band 4 (80 Hz has none);
   encircle a slope that changes by more than 5 dB, and with it SPL(*i*) when the slope is
   positive and rising, SPL(*i* − 1) when it has turned from positive to zero or negative; replace
   an encircled level by its neighbors' mean (band 24 by SPL(23) + s(23)); recompute the slopes
   with s′(3) = s′(4) and an imaginary 25th band; average three adjacent slopes; build the
   background SPL″ from SPL″(3) = SPL(3); and take the largest Table A36-2 factor of F = SPL − SPL″
   over F ≥ 1.5 dB: F/3 − 1/2, F/6, 10/3 below 500 Hz and above 5 kHz; 2F/3 − 1, F/3, 20/3 from
   500 Hz to 5 kHz inclusive. Bands 1–3 never carry a tone correction. With `--epnl_masking`
   (the default; below), a band below its noy threshold enters these steps at the threshold and
   carries no tone itself.
4. PNLT(*k*) = PNL(*k*) + C(*k*), and PNLTM its maximum (A36.4.3, A36.4.4).
5. The band-sharing check (A36.4.4.2): if C at PNLTM is below the mean of C over the five records
   centered there (fewer at either end of the record), PNLTM is PNL there plus that mean.
6. The duration correction D = 10 lg[(1/T) Σ 10^(PNLT(*k*)/10) Δt] − PNLTM, T = 10 s,
   Δt = 0.5 s, summed from the first to the last record at or above PNLTM − 10, dips between them
   included (A36.4.5). The regulation's −13 for 10 lg(0.5/10) is that number rounded (0.01 dB);
   NICE-OPS uses it unrounded.
7. EPNL = PNLTM + D (A36.4.1).

This is panam's `flight_acoustics.py` implementation (`noys`, `perceived_noise_level`,
`tone_correction`, `effective_perceived_noise_level`), ported term for term
(`src/pnl.cpp`), including its reading of the regulation where readings differ: 500 Hz is in the
middle range of Table A36-2, the duration is the contiguous interval, and the band-sharing
check is applied.

**Masked bands (`--epnl_masking`, on by default).** A measured band never falls below the
analyzer's noise floor, and certification practice excludes tones in bands masked by it (FAA
AC 36-4, Appendix 2, item 11: the tone-correction description includes "exclusion of tones in
masked bands"). A predicted band falls as far as absorption takes it. On the demo track over a
100 m grid (961 observers, B407), the 10 kHz band's median level at PNLTM is −194 dB, and a
near-silent database band reaches −1000 dB. Taken literally, step 3 reads that plunge — or an
inaudible band rising out of it, as 10 kHz does over 8 kHz for 81 of those observers — as a tone
of up to 20/3 dB: 615 of the 961 took their PNLTM tone correction from a band below its noy
threshold, and EPNL read a median 1.0 dB (up to 5.4 dB) high. So a band whose level is below SPL(d) of Table A36-3 — where
it stops contributing noys, and which no measurement would resolve — is the prediction's masked
band: it enters steps 1–7 at SPL(d), and carries no tone correction itself. Both halves are
needed: without the first, an audible band just before a plunge has its background crashed and
reads as a 20/3 dB tone; without the second, masked bands taken at their thresholds, which rise
17 and 21 dB at 8 and 10 kHz, can read as one. PNL is unchanged (a masked band is 0 noy either
way), and so is any record whose bands are all at or above their thresholds. On the grid above,
masking moved the PNLTM tone correction's median from 2.25 to 1.20 dB and the count of 3 dB or more
from 427 to 13, lowered EPNL by a median 1.0 dB, and raised it at one observer, by 0.2 dB (the
band-sharing check and the duration can pick other records). On the 2017 measured spectra it
moves EPNL by 0.01 dB rms. `--no-epnl_masking` gives the procedure as written, which is panam's.

**Not included.** Everything A36 does to a *measured* spectrum before step 1: the
instrumentation and background-noise corrections (A36.3), the adjustment to reference
conditions (A36.9; a prediction is already at its own atmosphere and flight path), and the
identification of pseudo-tones from ground reflections, which the regulation lets an applicant
exclude from the tone correction — NICE-OPS takes every tone the predicted spectrum has in its
audible bands, ground-effect interference included (the plate and point-receiver models put real
ripple into the spectrum, as a microphone does). There is no 10 dB-down window check of the
record's own completeness beyond the `clipped` flag below.

**`--epnl_file`** writes one row per observer, in observer order: `observer` (the index `-i`
uses), the microphone's `x`, `y`, `z` in `--units`, `EPNL` (EPNdB), `PNLTM` (TPNdB),
`duration_correction` (D, dB), `tone_correction` (PNLTM − PNL at its record: the record's own C,
or the five-record mean when the band-sharing check raised it), `tone_band` (Hz, the band C came
from; `nan` when the record has none), `bandshare_adjusted` (1 when step 5 raised PNLTM),
`PNLTM_time`, `duration_first`, `duration_last` (the centers of the PNLTM record and of the
duration's first and last records, s, on the trajectory's clock, in reception time) and `clipped`
(1 when the duration reaches the first or last record the observer heard, so the event may not be
whole). An observer with no finite PNLT — every frame below the noy floor — is written as
`-inf,-inf,nan,...`.

**`--band_frames`** writes a netCDF-4 file: `level(observer, frame, band)`, float, dB re 20 µPa,
each band's unweighted linear-average level over the frame; `time(frame)`, the frame centers (s,
on the trajectory's clock); `band_center(band)` (Hz, the bands inside `--band_limits`); the
microphones' `x`, `y`, `z` and `microphone_height` in `--units`; and, for a grid, `grid_x` and
`grid_y` (observer = *iy* × size(`grid_x`) + *ix*). Levels are NaN outside the frames an observer
received and −inf in a band that received nothing. Global attributes record the trajectory, the
database, `frame_length_seconds` (0.5) and `band_frames_version` (1). Each observer's record is
one chunk, shuffled and deflated. It is the file for whole-grid post-processing — EPNL variants,
C-weighted exposure, comparison with measured frames — that `--band_history`'s file per observer
and row per state cannot serve on a grid: on the demo grid (14,400 observers, 273 frames, 28
bands) it is 318 MB, where `--band_history` would write 14,400 files of 6,551 rows.

**Cost.** Both force every state to be evaluated (the windowed search bounds the A-weighted
level, which says nothing about where PNLT peaks), and hold every band of every state of an
observer while its frames are formed. On the demo grid with B407 (spectral, ground effect,
14 threads) the run takes 4.2 s without them, 5.6 s with `--epnl_file`, and 10.1 s with
`--band_frames` too, of which 2.2 s is the compression; `--band_frames` holds every observer's
frames until the end, 0.9 GB resident there.

**How it was checked.**

- *Against panam*, at 1e-6 dB: `tests/test_pnl.cpp` embeds values generated by running panam's
  own functions (`tools/pnl_reference.py`) on 27 spectra — every noy region, tones at every
  special band (80 Hz, 100 Hz, 500 Hz, 5 and 6.3 kHz, a band-24 tone on a slope, a slope turning
  to exactly zero), F above 20 and between 1.5 and 3, silence, random spectra — and 7 histories,
  each with the band-sharing check on and off (it fires in three, one of them at the record's
  end), plus panam's own test cases one for one. Mutating the tone correction one line at a time
  (the 500 Hz range, the band-24 extrapolation, step 3(b)'s "zero or negative", the truncated
  five-record mean, the −13 rounding) fails it.
- *Masking*, which panam does not have, against panam's procedure run on the reconstructed
  spectrum with the masked bands' factors zeroed (the same script, which first checks its
  per-band form of panam's tone correction against panam's own on every spectrum): the 27 spectra
  (unchanged wherever every band is audible), six prediction-like ones — absorbed top bands, an
  inaudible 10 kHz band rising out of the plunge, a −1064 dB band, `-inf` bands, an audible band
  before a plunge and masked bands rising along their thresholds — and a pass whose EPNL masking
  lowers by 1.8 dB. Dropping either half of the masking fails it.
- *The frames*: unit tests hold their energy to the history's exposure (1e-12) through a
  reception factor below and above one, a state straddling a boundary to its proportional split,
  and negative times to negative frames. Rebuilt independently in Python from a `--band_history`
  run's CSV files, the frames agreed with `--band_frames` to 6e-5 dB (the CSV's ten digits and the
  file's float), and panam's EPNL on them with `--epnl_file` to 5e-8 dB.
- *End to end* (CTest `epnl_*`): `check_band_frames` recomputes every EPNL and PNLTM from the
  netCDF frames and compares them with the CSV, checks that each PNLTM record's tone band is
  audible, checks the files against each other and for plausibility, and holds the frames'
  A-weighted energy to a `--full_history` run's SEL within 1e-4 dB. The energy check cannot see
  where in time the energy went, so the frames are also rebuilt from a `--band_history` run of the
  same case (each state spread over its own reception-time column's step) and must agree frame by
  frame: filing states at their emission time, or at t − r/c, fails it. The same case with
  `--no-epnl_masking` must read some EPNL more than 0.1 dB higher. The SEL output is unchanged by
  asking for either file, and the refusals are tested by their messages.

**Against the 2017 measurements.** On the four detail sample runs (calm-day approaches 286270
and 286277 with the shipped B407 database; windy-day steady runs 284205 and 284211 with the
held-out `min2_no284` database), over the microphones the paper's SEL gate passes, NICE-OPS's
EPNL was set against EPNL of the measured 0.5 s frames, both with masking. The measured frames stop
at 5 kHz, so for a like-for-like comparison the predicted frames' 6.3, 8 and 10 kHz bands were
replaced, as the measured side's must be, by the 5 kHz level less 5, 10 and 15 dB (a roll-off
that triggers no tone correction); predicted minus measured, dB, with the full 24-band
`--epnl_file` value against the same measured EPNL beside it:

| day, model | mic-runs | EPNL median | EPNL rms | PNLTM median | PNLTM rms | `--epnl_file` EPNL median, rms |
| --- | --- | --- | --- | --- | --- | --- |
| 286 calm, `station_plate` | 88 | −0.33 | 1.22 | −0.89 | 2.00 | −0.23, 1.21 |
| 286 calm, `stratified_plate` | 88 | −0.09 | 1.26 | −0.75 | 1.96 | +0.01, 1.28 |
| 284 windy, `station_plate` | 87 | −0.22 | 2.70 | −0.85 | 3.07 | −0.11, 2.70 |
| 284 windy, `stratified_plate` | 87 | −0.88 | 1.90 | −1.36 | 2.91 | −0.78, 1.91 |

Without masking (panam's procedure on both sides) the first four columns read −0.25 / 1.23,
+0.02 / 1.33, −0.22 / 2.85 and −0.80 / 1.87 for EPNL, reproducing an earlier computation from the
harness's per-state band histories, framed there by the states' reception times (1.24 and 1.88 dB
rms for the calm `station_plate` and windy `stratified_plate` rows), so the two framings agree.
Masking leaves the measured EPNL where it was (0.01 dB rms); what it changes is the prediction's
top bands: unmasked, the full 24-band EPNL read 0.04–0.19 dB (median, 0.24–0.41 rms) above the
5 kHz-limited one, and masked it reads the same (0.00–0.01 median, 0.11–0.21 rms), so the excess
was inaudible tones. About half the PNLTM bias is the tone correction: the measured spectra carry
a median 1.2–1.3 dB of it at PNLTM against 0.8–0.9 dB predicted, consistent with the blend of
spheres and grid cells smoothing the tones. EPNL, whose duration correction integrates every
record, is biased less.

### Where the code is

| Step | Source |
| --- | --- |
| Command line, the observer loop, metrics | `src/niceops.cpp`, `src/util/sel.hpp` |
| EPNL and reception-time frames | `src/pnl.cpp`, `src/band_frames.hpp`, `src/band_frames.cpp` |
| Trajectories and their parameters | `src/trajectory.cpp`, `src/geodetic.cpp`, `src/smoothing_spline.cpp` |
| Databases, blending, the sphere lookup | `src/sphere_database.cpp`, `src/interpolator.cpp`, `src/delaunay.cpp`, `src/condition_selection.cpp` (`--blend nearest`/`aam`) |
| Absorption and weighting | `src/atmosphere.cpp`, `src/atmosphere_profile.cpp`, `src/weighting.cpp` |
| Ground and plate microphones | `src/ground_effect.cpp`, `src/faddeeva.cpp`, `src/plate_receiver.cpp` |
| Refraction | `src/propagation.hpp`, `src/refraction.cpp`, `src/stratified.cpp` |
| Parabolic equation and exact solutions | `src/pe.cpp`, `src/pade_propagator.cpp`, `src/exact_solutions.cpp`, `src/airy.cpp` |
| Propagation loss tables | `src/propagation_table.cpp` |
| Terrain | `src/terrain.cpp`, `terrain_ray` in `src/propagation.hpp` |

## Running NICE-OPS

### Command-line options

`niceops --help` prints the same list.

#### Inputs and outputs

| Option | Meaning |
| --- | --- |
| `-d`, `--database` | netCDF noise database. Once, or once per trajectory. |
| `-t`, `--trajectoryfile` | CSV trajectory. Repeat for several flight paths. |
| `-v`, `--vehicle` | JSON vehicle definition. Once, or once per trajectory. |
| `-o`, `--output` | SEL output file. Once per trajectory (default `integrated/197421_LAE.csv`); with several tracks and one name or none, each output is named after its track. Two trajectories asked to write one file are refused. |
| `-b`, `--observer` | Scattered observer file; switches off the default grid. |
| `-S`, `--state` | Save the derived trajectory state. With several tracks and one name, each is saved as `<track>_state.csv` beside it. |
| `-i`, `--history` | Prefix for per-observer time history files, `<prefix>_i<N>.csv`. With several tracks each file is also named for its track, `<prefix>_<track>_i<N>.csv`, so they do not overwrite one another. |
| `--band_history` | Prefix for per-observer band time history files (see Outputs), named as `-i` names its files. Needs the spectral path. |
| `--epnl_file` | Write each observer's EPNL, PNLTM, duration and tone corrections (14 CFR 36 A36.4) to this CSV, from one-third-octave levels in 0.5 s reception-time frames. Once, or once per trajectory; one name for several tracks becomes `<track>_EPNL.csv`. Needs the spectral path and the 24 bands 50 Hz–10 kHz. See [EPNL and reception-time frames](#epnl-and-reception-time-frames). |
| `--epnl_masking`, `--no-epnl_masking` | `--epnl_file`'s tone correction: a band below its noy threshold (Table A36-3's SPL(d)) enters it at the threshold and carries no tone itself, as tones in masked bands are excluded from a measured spectrum's. Default on; `--no-epnl_masking` is the procedure as written. See [Masked bands](#epnl-and-reception-time-frames). |
| `--band_frames` | Write each observer's band levels in 0.5 s reception-time frames to this netCDF file, (observer, frame, band). Once, or once per trajectory (`<track>_band_frames.nc`). Needs the spectral path. |
| `-m`, `--metricfile` | Append overall noise metrics to this file: the case name, the energy mean of the per-observer SELs (10 lg of the mean of 10^(L/10), not a mean of decibels), their maximum, and the share above 65 dB. |
| `--origin` | WGS84 reference for the local frame: latitude longitude, degrees. Defaults to the first position of the first geodetic track. |
| `--frame_bearing` | Compass bearing of the local +x axis. Default: x east, y north. The shipped AC tracks use 210. It also turns the AC format's compass heading into the local frame, for the hover nose and `--nose heading`. |
| `--ground_elevation` | Elevation of the observer plane above mean sea level, in `--units`. Subtracted from a WGS84 track's altitude column. With `--terrain` it is the height of the local z = 0 on the trajectory's datum, and an AC track's own alt − z must match it. |
| `--path_smoothing` | How far, RMS per axis in `--units`, the path fitted through a WGS84 track may sit from the given positions. Default 10 ft. Zero interpolates exactly, which on measured positions makes the load factor track coordinate rounding. |
| `--start_utc` | Seconds since UTC midnight at *t* = 0 of a WGS84 track, for the day-night window. Required with `-D`, `--lden_file` or `--cnel_file` when the track has no `utcsec` column (it used to default to midnight, which is night, so the whole footprint took the 10 dB penalty). |
| `--latlon` | Also write each level file as latitude, longitude, level. |
| `--terrain` | A bare-earth elevation grid (netCDF, e.g. USGS 3DEP). Grid observers stand on it, and the ground model and plate take each path's mean ground plane instead of the plane z = 0. Needs the run's frame (`--origin`, `--frame_bearing`). Not with `--refraction` or `--propagation_table`. See [Terrain](#terrain). |
| `--terrain_offset` | Added to every `--terrain` height, in `--units`, to put the grid on the trajectory's vertical datum. Default 0. At Amedee −23.6 m (−77.43 ft): 3DEP is NAVD88, the 2017 tracks ellipsoidal. |
| `--terrain_step` | Spacing, in `--units`, of the terrain samples each path's mean ground plane is fitted to. Default: the grid's own spacing, at least 5 m. |
| `--observers_on_terrain` | With `--terrain` and `-b`: take each observer's ground height from the terrain instead of the file's z. |
| `-l`, `--logfile` | Log file (default `log/niceops.log`). |
| `-L`, `--loglevel` | `trace`, `info`, `warn`, `error` or `fatal` (default `info`). |

#### Physics and metrics

| Option | Meaning |
| --- | --- |
| `-a`, `--ambient` | Ambient level in dB, added to every sample's level in energy (two equal levels give +3 dB). Off unless given. It used to be added as pressure amplitude, as the D original did, which read up to 3 dB high. |
| `-q`, `--qsam` / `--no-qsam` | Use the flight path angle relative to the tip-path-plane normal rather than the vertical. On by default. |
| `-s`, `--step` | Resampling step for trajectory adaptation. **Not** the SEL integration step, which comes from the trajectory's own sampling. |
| `--units` | Length unit for the command line and for every file carrying a length: `ft` (default) or `m`. See **Units**. |
| `--trim` | Crop each trajectory from this x, in `--units`, to its lowest/slowest point. Off by default; it changes which part of the flight is integrated. |
| `--grid`, `--gridstep` | Observer grid bounds and spacing, in `--units`. Default ±1500 m (±4921.26 ft), step 25 m (82.021 ft) — a fixed extent either way. |
| `-D`, `--dnlfile` | Write the day-night level (DNL), summed over every trajectory: the energy average over a 24-hour day, with the 22:00–07:00 night +10 dB. A daytime operation reads its SEL less 49.37 dB. |
| `--lden_file` | Write the day-evening-night level (Lden, EU Directive 2002/49/EC): the same 24-hour average with 19:00–23:00 +5 dB and 23:00–07:00 +10 dB. |
| `--cnel_file` | Write the community noise equivalent level (CNEL, California Title 21): 19:00–22:00 +10 lg 3 = 4.77 dB and 22:00–07:00 +10 dB. |
| `--lamax_file` | Write each observer's LAmax: the highest A-weighted level over every operation, time-weighted per `--lamax_weighting` (LAFmax by default), dB. Forces the whole history. **Before this option existed it wrote the instantaneous peak**, 0.7–0.9 dB above LAFmax at the median; `--lamax_weighting instantaneous` still gives that, to the bit. |
| `--lamax_weighting` | The time weighting `--lamax_file`, `--nax_file` and `--ta_file` read the level through, in reception time: `fast` (IEC 61672-1 F, τ = 0.125 s; default), `slow` (S, τ = 1 s), `leq0.5` (the largest running 0.5 s Leq, ECAC Doc 32's LASmax surrogate) or `instantaneous`. A trajectory sampled more coarsely than half the time constant is refused. Refused without one of the three outputs. |
| `--nax_file`, `--nax_threshold` | Write each observer's NAx: the number of operations whose own LAmax (per `--lamax_weighting`) reaches the threshold, dBA (required with `--nax_file`: there is no default, so the study's threshold is always stated). The count adds across trajectories. EPNL is `--epnl_file`. |
| `--ta_file`, `--ta_threshold` | Write each observer's time above: seconds of reception time the level, time-weighted per `--lamax_weighting`, exceeds the threshold, dBA (required, no default), summed over every trajectory. |
| `--csel_file` | Write the C-weighted sound exposure level (LCE) per observer, by SEL's window rule. One name, or one per trajectory; one name with several trajectories writes `<track>_LCE.csv` beside it. Spectral path only; refused with `--no-spectral`, a database without spectra, `--db_cutoff` or `-a`. |
| `--coverage_file` | Write, per observer, the energy-weighted share of its levels that came from cells the database did not measure. Needs a database that records `coverage`, as the 2017 aircraft databases do. With spectra, each cell counts by the share of its A-weighted band energy at the sphere that was measured, so a band the builder's SNR gate dropped counts by what it carries (on B407, 5–10 kHz in about half the measured cells, and 0.4% of the lower hemisphere's energy); channel coverage is all or nothing. On the demo grid B407 reads 1.2% unmeasured at the median, 27% at the 95th percentile. Straight rays only: refused with `--refraction` or `-a`. |
| `--full_history` | Integrate each exposure over the whole history rather than the event's 10 dB-down window: SEL, C-SEL, and each operation's share of DNL, Lden and CNEL. Slower, since no window is found. |
| `--utc_offset` | Hours from UTC to local time, for the DNL, Lden and CNEL periods. |
| `--powerfile` | Write per-state shaft power, when the database has it. |
| `--spectral`, `--no-spectral` | On by default when the database carries spectra: propagate its per-band spectra to the ground with ISO 9613-1 absorption. `--no-spectral` takes one broadband level with the database's EAA factor instead (databases without spectra always do). |
| `--band_limits` | Frequency range integrated into the metric, Hz: `low high`. Default `20 20000`. |
| `--humidity`, `--temperature`, `--pressure` | Atmosphere for `--spectral`, overriding the vehicle file. Percent, K, kPa. The spectral path refuses to run without a humidity (no `--humidity`, no vehicle-file value, no `--atmosphere`): absorption depends strongly on it, so there is no silent 20 %. |
| `--atmosphere` | A stratified atmosphere: a CSV of temperature, humidity, and optionally pressure and wind, against height. Absorption is integrated through it. Spectral path only. See below. |
| `--wind` | Build the stratified atmosphere from one surface wind instead of reading `--atmosphere`: `SPEED,DIRECTION[,HEIGHT]`. The speed names its unit (`12kt`, `6.2mps`, `14mph`); the direction is the compass bearing, degrees, it blows **from**; the height is the anemometer's, in `--units` (default 10 m). Neutral log law, dry-adiabatic lapse from the run's surface air. Spectral path only; not with `--atmosphere`. See below. |
| `--roughness_length` | `--wind`: aerodynamic roughness length z₀, in `--units`. Default 0.03 m (0.0984 ft), open flat terrain. Not `--roughness`, the acoustic ground's RMS height. |
| `--surface_layer_depth` | `--wind`: height, in `--units`, to which the wind follows its law; above it the wind is held. Default 200 m (656.2 ft). |
| `--stability` | `--wind`: `neutral` (default, the log law), or a TA Luft class `I`, `II`, `III/1`, `III/2`, `IV`, `V` (Pasquill `F` to `A`): a power-law wind and a temperature gradient in the surface layer. Unvalidated; a run that takes one warns. |
| `--profile_levels` | `--wind`: levels per factor of ten in height (default 16). Results do not depend on it; raising it checks that they do not. |
| `--write_atmosphere` | `--wind`: write each trajectory's built profile in the `--atmosphere` format, so the run can be repeated and inspected. With several tracks and one name, `<track>_atmosphere.csv` beside it. |
| `--refraction` | `none` (default); `linear`, curved rays through the `--atmosphere` (or `--wind`) profile's effective sound speed fitted to a straight line per path; or `stratified`, through the profile's own layers. Both give shadow zones under upward refraction. See below. |
| `--shadow_floor` | Most a refractive shadow zone attenuates, dB. Default 20. |
| `--shadow_model` | `--refraction stratified` with `--ground_effect`: how each band's shadow boundary is set. `legacy` (default): interpolated between the layered rays' boundary and the path's linear fit's. `resolved`: from the profile each band resolves (its tangent below the height where the band gathers `--shadow_phase` of vertical phase), with a cramped band's field from `--cramped_ground`. Opt-in; see "Low sources near grazing" below. |
| `--shadow_phase` | `--shadow_model resolved`: the vertical phase, radians, a band must gather to resolve the near-ground profile. Default 1, fixed a priori (the radian the layered rays' flat-floor limit takes); 0.5 and 2 are reported as sensitivities only. |
| `--cramped_ground` | `--shadow_model resolved`: a band whose boundary is within two diffraction lengths of the source takes the residue series at the source's own height (`series`, default); or, on its lit side, the ground model along the straight path (`straight`) or on the ray's equivalent geometry (`ray`). |
| `--tunneling` | `--refraction stratified`: also take the rays turned back from aloft that turn again just above the microphone, which the wave reaches through the layer below them (an Airy tail per band). Off by default; see "Micro-shadows and skip zones" below. |
| `--skip_fill` | `--refraction stratified`: fill the skip zone short of where rays from aloft first land with the caustic's Airy tail per band. Off by default; see below. |
| `--scattering` | Fill refraction shadows with the energy turbulence scatters into them, at HARMONOISE strength γT, instead of `--shadow_floor`. Default 0 (off). Needs `--refraction` and `--ground_effect`. |
| `--nose` | Where the sphere's nose points: `velocity`, the velocity `--speed_reference` selects; `ground_velocity`, the ground velocity whatever `--speed_reference` says; or `heading`, the track's compass heading at every speed, so a crabbing aircraft's sphere turns with its airframe. `heading` needs a frame (`--frame_bearing` for an AC track) and a heading column. The first two turn to the heading below 8 kt. Default: what the database's `azimuth_reference` records, `velocity` if it records nothing. See below. |
| `--speed_reference` | Which velocity sets the flight condition (μ, the path angle and, under `--nose velocity`, the sphere's nose): `ground`, or `air`, the ground velocity less the wind at the aircraft (from `--atmosphere` or `--wind`). Default: what the database records it was built with, ground if it records nothing. See below. |
| `--ground_effect`, `--no-ground_effect` | Add excess ground attenuation (Chien–Soroka over a Delany–Bazley impedance) to the spectral path. On by default wherever the spectral path runs on a database that carries spectra; `--no-ground_effect` turns it off. Asking for it with `--no-spectral` is an error. |
| `--max_load_factor` | Largest load factor, in g, that a `fixed_load_factor` database's load-factor correction is taken to (default 2). Above it the correction is held, and the run says how many states it held. |
| `--load_factor_correction` | `one_sided` (default): a `fixed_load_factor` database's law, 20 lg(C<sub>T</sub>/C<sub>T,1g</sub>), applies above 1 g only, as panam applies it. `two_sided`: below 1 g too, down to `--min_load_factor`. Refused on a database whose load factors are spheres. See [The load-factor law below 1 g](#the-load-factor-law-below-1-g). |
| `--min_load_factor` | Under `--load_factor_correction two_sided`, the smallest load factor, in g, the law is taken to (default 0.5, at most 1). Below it the correction is held, and the run says how many states it held. Refused without `two_sided`. |
| `--blend` | How a state's condition picks the spheres it is evaluated from: `barycentric` (default), `nearest` (the nearest stored condition at every state) or `aam` (AAM v3.1's bracketing search). See [Other ways to choose the spheres](#other-ways-to-choose-the-spheres---blend). |
| `--accept_vehicle_mismatch` | A vehicle file that disagrees with the database by more than 2 percent in tip speed, rotor radius, weight or thrust coefficient at 1 g is refused. Between 1 and 2 percent it is a warning; this flag turns the refusal into a warning. |
| `--flow_resistance` | Ground flow resistance, kPa·s/m². 30 snow, 225 grass (default), 650 soil, 3000 dirt, 10000 concrete, 50000 asphalt. |
| `--receiver_height`, `--turbulence` | Point receiver's microphone height (in `--units`, default 1.2 m; an observer file's fifth column overrides it per observer) and turbulence (rad·s/√m, 0–16e-4) for `--ground_effect`; the height also places the microphone for `--propagation_table`. |
| `--impedance_model` | Ground impedance: `delany_bazley` (default), `miki`, `variable_porosity` (with `--porosity_rate`, 1/m) or `delany_bazley_layer` (with `--layer_depth`, in `--units`). `--flow_resistance` is the (effective) flow resistivity. |
| `--turbulence_model`, `--gamma_t` | `chessell` (default, `--turbulence`) or `harmonoise`, whose coherence takes `--gamma_t` (~1e-6 moderate, ~1e-5 strong). |
| `--roughness` | RMS roughness height, in `--units`: the coherent reflection falls by exp(−2(kσ cos θ)²). Default 0. |
| `--plate_table` | Observers are flush microphones on a rigid ground plate (GRAS 67AX in a GR1425 plate). Takes a netCDF table from panam; needs `--ground_effect` over the ground the table was built for. See below. |
| `--propagation_table` | A propagation loss table (netCDF, from `niceops_propagation_table`): each band's level relative to free field by bearing, source height, receiver height and range, for one stratified atmosphere over uniform ground. Replaces the ground and ray models wherever it covers a pair; pairs it does not cover take the run's own model, `--refraction` included. A plate microphone reads it at its lowest receiver height plus the plate's offset over bare ground. Spectral path. See below. |
| `--exact_ground` | Evaluate the ground model (or the plate's reflection coefficients) at every observer and state instead of interpolating per-trajectory tables of it: 3-4x slower, for checking the tables. |
| `--plate_mic_azimuth` | Which way the microphone is offset on its plate: `outboard` (default, away from the track) or degrees in the local frame from +x toward +y. |

#### Threads and culling

| Option | Meaning |
| --- | --- |
| `-j`, `--threads` | Worker threads for the observer loop. 0 (default) uses one per hardware thread; 1 runs serially. |
| `--db_cutoff` | Skip evaluations that provably fall below this level. See the caveat below. |

#### Spheres and diagnostics

| Option | Meaning |
| --- | --- |
| `--export_sphere` | Interpolate a whole noise sphere and write it as a one-group `.nod` database. |
| `--export_aam` | The same sphere in AAM/RNM format. Needs a database with spectra. |
| `--sphere_condition` | The operating condition for it: advance ratio, flight path angle, C<sub>T</sub>. `--blend` and `--load_factor_correction` apply, so the sphere is what a state at that condition would be evaluated from. |
| `--self_test` | Validate every sphere in the database against its own stored levels. |
| `-p`, `--plotting` | Plot the footprint with the flight track over it; `--cmin`/`--cmax`/`--cstep` set the contour levels. |
| `-f`, `--figfile`, `-F`, `--statefigfile` | Save the footprint and trajectory-state figures to files. |

### Units

`--units` selects the length unit for the whole interface: `ft` (the default) or `m`. It covers

- every length on the command line — `--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
  `--receiver_height`, `--path_smoothing`, `--terrain_offset`, `--terrain_step`,
  `--roughness_length`, `--surface_layer_depth`, and the anemometer height in `--wind` (whose
  speed names its own unit, as the wind columns of an `--atmosphere` file do);
- every file carrying a length — trajectory input, observer lists, level output, time histories;
- the axis labels on the figures.

The computation is metric throughout. The conversion happens once, where the command line is
parsed, and nothing downstream of that seam sees anything else.

**`--units` says what your files already are; it does not convert them.** The shipped tracks and
observer sets are in feet, so running them with `--units m` reads feet as meters and flies the
aircraft three times too far. Nothing can detect that, so the flag is a declaration you have to
get right. Where a file names its own unit in a column header — `alt_ft`, `alt_m` — the name is
checked against `--units` and a disagreement fails the file rather than being silently believed;
a bare `alt` claims nothing and is always accepted.

Two things stay put whatever you select:

- **Defaults are physical quantities, not numbers in the current unit.** The observer grid
  defaults to ±1500 m in steps of 25 m whichever unit is in force — it reads as ±4921.26 ft or
  ±1500 m, the same ground either way. If defaults scaled with the flag, `--units m` would
  quietly shrink the footprint to a third of its area.
- **The `-S` state file is metric.** It is a dump of the internal state, and
  `golden/state_197421.csv` — the D-derived oracle — encodes it that way.

Non-length quantities keep the units of the standards they come from, because converting them
would make them harder to check against the source: temperature in K, pressure and flow
resistance in kPa (ISO 9613-1 is written in kPa), frequency in Hz, angles in degrees, time in
seconds, levels in dB.

> Before this flag existed, `--grid`, `--gridstep`, `--receiver_height` and `--path_smoothing`
> were meters while `--trim` and `--ground_elevation` were feet and every file was feet. Scripts
> written against that need those four converting, or `--units m` adding.

### Output files

- `-o` — sound exposure level (SEL) per observer. Gridded output writes the x and y axes on the
  first two lines and then one row of levels per y.
- `-D`, `--lden_file`, `--cnel_file` — day-night, day-evening-night and community noise
  equivalent levels, summed over every trajectory in the run (see below).
- `--lamax_file`, `--nax_file`, `--ta_file` — the time-weighted maximum level, the NAx count and
  the time above a threshold (seconds), over every trajectory in the run.
- `--csel_file` — the C-weighted SEL per observer, one file per trajectory.
- `-i` — per-observer level time histories.
- `--band_history` — per-observer band time histories, for comparing with measured spectra
  frame by frame. One file per observer, `PREFIX_i<N>.csv`, with a row per trajectory state:
  `emission_time`, `reception_time` (s), `range` (the ray tube's under `--refraction`),
  `elevation` (the emission direction below horizontal, deg), the sphere's `phi` and `theta`,
  `level_dBA` (the A-weighted, band-limited level `-i` writes), then `band_<Hz>`: each band's
  unweighted level at the microphone, dB re 20 µPa, with spreading, absorption, the ground or
  plate and any shadow loss applied. A-weighted and summed in energy the bands give
  `level_dBA`. Bands outside `--band_limits` and culled states are left empty; `-a`'s ambient
  goes into `level_dBA` only. Needs the spectral path.
- `--epnl_file` — per-observer EPNL, PNLTM, duration and tone corrections, and the times of the
  PNLTM record and the duration, one row per observer; `--band_frames` — every observer's band
  levels in 0.5 s reception-time frames, as netCDF. Both are described under
  [EPNL and reception-time frames](#epnl-and-reception-time-frames).
- `-m` — appends the energy-mean SEL, maximum and above-threshold counts for the case.
- `-S` — the derived trajectory state (t, position, velocity, load factor, C<sub>T</sub>, μ, flight path angles).
- `--powerfile` — per-state shaft power, if the database carries any.
- `--latlon` — writes a `_latlon` companion beside every level file above, as latitude, longitude, level.
- `--export_sphere` — a whole interpolated noise sphere, written as a one-group `.nod`.

### Several trajectories, and the day-night level

`-t` may be repeated. `-d`, `-v`, `-o`, `-S`, `--epnl_file` and `--band_frames` then take
either one value, shared by every trajectory (the output names are then derived per track), or
one per trajectory.

```sh
./build/niceops -t morning.csv -t evening.csv -t night.csv -D integrated/dnl.csv
```

DNL is a property of a set of operations rather than of one flight: each trajectory's SEL-style
integral gains a 10 dB penalty on samples flown at night, and the results are summed as energy.
Night is 22:00–07:00 local; `--utc_offset` shifts the trajectory's UTC seconds into local time.
The sum is normalized to a 24-hour day, as the DNL is defined: a single daytime operation
gives the SEL less 10 lg 86 400 = 49.37 dB. (Before the normalization, `-D` wrote the raw sum,
which read 49.37 dB high.)

`--lden_file` and `--cnel_file` are the same sum under other schedules. Each penalty is added to
the level of every sample heard in its period, local time:

| Metric | Option | Evening | Night | Definition |
| --- | --- | --- | --- | --- |
| DNL | `-D` | — | 22:00–07:00, +10 dB | 14 CFR Part 150 |
| Lden | `--lden_file` | 19:00–23:00, +5 dB | 23:00–07:00, +10 dB | EU Directive 2002/49/EC, Annex I |
| CNEL | `--cnel_file` | 19:00–22:00, +10 lg 3 = 4.77 dB | 22:00–07:00, +10 dB | California Code of Regulations, Title 21, §5001 |

The directive writes Lden = 10 lg[(12·10^(L<sub>day</sub>/10) + 4·10^((L<sub>evening</sub>+5)/10)
+ 8·10^((L<sub>night</sub>+10)/10)) / 24], with each L the equivalent level over its own period;
that is the same number as the penalized energies summed and divided by 86 400 s, which is how
it is computed. Title 21 weights each evening operation three times and each night operation ten
times, which is the same thing again. The periods are the directive's defaults; a member state's
shifted evening is not offered. A period's opening hour belongs to it: 07:00 is day, 19:00
evening. A unit test integrates a steady level over a whole day and recovers each metric's period
weighting, 10 lg((15 + 9·10)/24) for DNL, 10 lg((12 + 4·10^0.5 + 8·10)/24) for Lden and
10 lg((12 + 3·3 + 9·10)/24) for CNEL, to 1e-9 dB; CTests run the demo track by day, at 22:30 and
at 23:30 local and check which metrics coincide there.

**Each operation's window.** By default an operation contributes its own penalized 10 dB-down
window: the window is found on the penalized history, as SEL's is on the raw one, so a night
operation is integrated over the same samples a day one would be. `--full_history` integrates
each operation's whole history instead, as AAM and NORAH2 do. The window leaves out what lies
more than 10 dB down: on the four 2017 sample runs of the validation harness (calm-day
approaches and windy-day steady flight, band histories at every microphone) the whole history is
+0.29 to +0.56 dB above the window at the median per operation, up to +1.4 dB at a single
microphone. For a DNL meant to sit beside an AAM or NORAH2 contour, use `--full_history`; for one
meant to match a measured SEL's convention, leave it off. Before this release `-D` used the
window whatever `--full_history` said.

The clock is the AC format's `utcsec` column, or `--start_utc` for a WGS84 track, and a sample
takes the penalty of the period it is emitted in (the propagation delay, seconds, is not added);
a track without one is refused under any of the three, since classing it as day (or, from a
zero default, night) would be a made-up answer.

Plotting is skipped for multi-trajectory runs, since one figure file cannot hold several paths.

A trajectory that cannot be read is logged and skipped rather than ending the run, so one bad
file in a batch costs that file alone. The exit status is nonzero if anything was skipped.

### Maximum level, time above and C-weighted SEL

**Time-weighted maximum (`--lamax_file`, `--lamax_weighting`).** A sound level meter does not
report the instantaneous mean square; it reports it through a time weighting. IEC 61672-1:2013
defines F (fast) and S (slow) as first-order exponential averages of the squared,
frequency-weighted pressure,

    y(t) = (1/τ) ∫ p²(s) exp(−(t − s)/τ) ds,   τ = 0.125 s (F) or 1 s (S),

and LAFmax is the largest value of 10 lg y. `--lamax_weighting` picks the weighting:

| Value | Maximum | Coarsest trajectory step |
| --- | --- | --- |
| `fast` (default) | LAFmax, IEC 61672-1 F | 0.0625 s |
| `slow` | LASmax, IEC 61672-1 S | 0.5 s |
| `leq0.5` | the largest running 0.5 s Leq: ECAC Doc 32's "LASmax", which its footnote calls a surrogate for S, and NORAH2's | 0.25 s |
| `instantaneous` | the largest sample of the history: no time weighting | any |

The weighting runs in **reception time**. The history is a level per emission time, and each
sample lasts T (1 + ṙ / c) at the microphone, the same weights the SEL integral uses
([Exposure](#exposure)). The samples are taken as a piecewise-constant mean square over those
durations, back to back, and the weighting is evaluated on that signal exactly: across a sample
the exponential relaxes toward it as exp(−s/τ), so its maximum falls at a sample's end, and the
running average is linear between the instants where either edge of its window crosses a sample
boundary. Nothing is resampled or blocked, so the answer depends on neither a grid nor a phase.
The detector starts settled on the first sample's level. `--nax_file` counts the operations whose
own weighted maximum reaches `--nax_threshold`, and the maximum over several operations is the
largest of theirs.

A weighting is refused on a trajectory whose sample interval exceeds half its time constant
(the table's last column): a fast detector fed a 0.5 s history settles within every sample, and
its "maximum" is the largest sample again under another name. The check is on the mean interval,
the step every sample is weighted by; the 2017 AC tracks are 0.02 s with the odd dropout up to
0.14 s, and pass. A sparse WGS84 track can be refitted onto a finer `-s`; a 27-column AC track is
used as sampled, so a coarse one needs a finer file. WGS84 waypoints at the default `-s 0.5`
are refused under `fast`, which they were not before this option existed.

**What changed.** Before `--lamax_weighting`, `--lamax_file` wrote the instantaneous maximum of
the 0.02 s history, which is not a meter quantity. Against the validation harness's LAFmax on the
calm day it sat +0.95 dB high at the median (95th percentile +2.29 dB); on the four sample runs
below, the instantaneous maximum is +0.72 to +0.93 dB above the in-tool LAFmax at the median (95th
percentiles +1.1 to +3.8 dB). `--lamax_weighting instantaneous` reproduces the old file to the
bit (checked against the build before the change on all four runs). On those runs LAFmax is
+0.52 to +0.92 dB above LASmax and +0.09 to +0.16 dB above the running 0.5 s Leq, at the median.

**How it was checked.** Unit tests hold the exponentials to their closed forms: the step
response, with the time above half power, at three sampling steps; IEC 61672-1's tone-burst
response 10 lg(1 − e^(−T<sub>b</sub>/τ)), which its Table 4 lists as −1.0 dB (F, 200 ms),
−18.0 (F, 2 ms), −27.0 (F, 0.25 ms) and −7.4 (S, 200 ms); a burst cut into unequal samples
reading as the same burst; a source closing at Mach 0.5, whose 0.2 s burst the detector hears as
0.1 s; and the running Leq's ramps, including a crossing that falls where the window's trailing
edge, not a sample boundary, turns. CTests pin the call site on the demo case's broadband database
(`lafmax_golden`, `nax_golden`, `ta_golden`): the goldens were checked against a brute-force
detector run on a 5e-5 s reception-time grid built from the run's `-S` state and `-i` histories,
which agrees with them to 1e-4 dB, where the same detector in emission time is up to 0.018 dB
off; NA64 there counts three cells fewer than the instantaneous peak would. On the 2017 records
(`station_plate` on calm-day approaches 286270 and 286277, `stratified_plate` on windy-day runs
284205 and 284211, 188 microphones; the harness's own commands with this build):

| In-tool LAFmax minus | Median, per run | 5th to 95th percentile | Largest difference |
| --- | --- | --- | --- |
| an independent evaluation of the same definition from `--band_history`'s reception times | — | — | 2e-6 dB (the file's printed digits) |
| the harness's procedure without its 0.1 s blocks (emission time, 0.02 s) | −0.003 to −0.005 dB | −0.043 to +0.015 dB | 0.075 dB |
| the harness's LAFmax (energy-averaged into 0.1 s blocks, then F) | +0.03 to +0.06 dB | −0.008 to +0.32 dB | 0.55 dB |

The last row is the harness's 0.1 s blocking, which it applies to put the prediction on the
measured frames' footing: blocking the same reception-time signal into 0.1 s blocks moves it by
+0.05 to +0.07 at the median, while the clock alone (reception against emission time, both
blocked) moves it by −0.002 to +0.015. Within 0.05 dB of the harness's LAFmax are 44–74 % of
microphones per run; within 0.08 dB of its unblocked procedure, all of them. Scoring a
prediction against the measured 0.1 s frames should keep the harness's blocking; the in-tool
value is the meter's reading of the continuous prediction.

The harness repeated the comparison over every microphone-run of the paper's cases on this
build (the in-tool LAFmax against its own 0.1 s-blocked LAFmax, same microphone-runs): median +0.033 to +0.046 dB per case, 5th to 95th percentile −0.018 to
+0.203 dB, 0.71 dB at worst, 53–62 % within 0.05 dB. Scored against the measurements the two
differ by 0.02 dB rms or less, so the harness keeps its own blocked value for scoring.

**Time above (`--ta_file`, `--ta_threshold`).** Seconds of reception time during which the
level, through the same weighting, exceeds the threshold (strictly), over each operation's own
history, summed over operations. Within a sample the weighted level is followed exactly, so a
crossing inside a sample is credited at the instant it happens; under `instantaneous` each sample
counts whole. On the 2017 sample runs, with instantaneous levels, TA65 is 11–36 s and TA75
4–10 s at the median microphone.
There is no default threshold, as with NAx.

**C-weighted SEL (`--csel_file`).** The same exposure integral over the C-weighted level, LCE.
The C curve is IEC 61672-1:2013 Annex E's, in the pole-zero form the A curve already uses
(`src/weighting.cpp`): the A curve's outer pole pairs alone, normalized at 1 kHz, and in fact
exactly the second factor of panam's A-weighting, constant and all. It matches the standard's
Table 3 within its 0.1 dB rounding from 10 Hz to 20 kHz, and the difference between the A and C
curves is the A curve's two inner pole pairs to 1e-12 dB. Each state's C-weighted level is formed
from the per-band levels at the microphone that its A-weighted level is summed from (spreading,
absorption, ground or plate, refraction, every arrival), over the same `--band_limits`. **The
window rule is SEL's**: the 10 dB-down window of the C-weighted history itself, which is not the
A-weighted window, or the whole history under `--full_history`; the reception-time weights are
SEL's too. It needs band levels, so it is refused on the broadband path (`--no-spectral`, or a
database without spectra); with `--db_cutoff`, whose bound is proven on the A-weighted level and
could skip a state with C-weighted energy; and with `-a`, whose ambient is a single A-weighted
level.

C-SEL is dominated by low frequencies, so **the band limits matter**: on the 2017 sample runs it
is 8.8 to 9.5 dB above the A-weighted SEL at the median microphone, and on run 286270 40 % of
the C-weighted energy (median microphone) lies below 50 Hz, so the default 20 Hz lower limit leaves out a database's 10–16 Hz
bands. CTests check that a run limited to the 1 kHz band gives the A-weighted SEL to 1e-5 dB
(where both curves are 0 dB), which exercises the window, the reception weights and the per-state
bands, and that a full-band run does not.

## Inputs

### Input files

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV. Column 0 is UTC seconds, column 1 elapsed time; positions and velocities are in **feet**. |
| Trajectory (`-t`), WGS84 | A CSV whose header names `time`, `lat`, `lon` and `alt` columns — see [Geodetic coordinates](#geodetic-coordinates). Any column order, and optional `vE`/`vN`/`vUp`. |
| Trajectory (`-t`), static | A CSV whose header names `duration`, `x`, `y`, `z` and `heading`: each row is a hold at one point for `duration` seconds, nose on the compass `heading`. Every hold shares one position, in `--units`; a move is a trajectory. Sampled at 50 Hz, at 1 g with zero velocity. |
| Database (`-d`) | netCDF `.nod`: one group per operating condition (advance ratio μ, flight path angle γ, thrust coefficient C<sub>T</sub>), each holding scattered (φ, θ) channels of dBA level and EAA excess-attenuation factor. |
| Vehicle (`-v`) | JSON with `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg) and `Atmosphere.density`. |
| Observers (`-b`) | CSV of scattered observer positions, in feet: `x,y,z`, `id,x,y,z`, or `id,x,y,z,height` with each microphone's height above the ground (see [Receivers](#receivers)). The first line is always taken as a header. Without it, a 120×120 grid spanning −1500…1475 m is used. |
| Terrain (`--terrain`) | netCDF: 1-D `latitude` and `longitude` (degrees, strictly increasing) and `elevation(latitude, longitude)` in **meters**, whatever `--units` says; an optional global `vertical_datum` text attribute. Missing nodes are NaN, `_FillValue` or `missing_value`. See [Terrain](#terrain). |

Input and output files use **feet**; the computation is metric internally.

### Geodetic coordinates

A trajectory can be given in WGS84 coordinates — latitude, longitude and altitude — instead of
the 27-column AC state format. That is what a mission planner has: a route, not a flight-test
record.

```
time,lat,lon,alt
0.0,39.82914323,-74.40996141,1134.60
0.5,39.82903000,-74.41007000,1133.80
...
```

Both formats are read by `-t`, and which one a file is is decided by its width: 27 columns is
the positional AC format, anything else is read by the names in its header. Column order does
not matter and unrecognized headers are an error rather than a guess.

| Column | Meaning |
| --- | --- |
| `time` | Seconds from the start of the track. `utcsec` (seconds since UTC midnight) serves instead, and sets the clock for `-D`. |
| `lat`, `lon` | Degrees. |
| `alt` | **Feet above mean sea level.** `--ground_elevation` gives the elevation of the ground the observers stand on, which is subtracted to get height above them. |
| `vE`, `vN`, `vUp` | Optional, ft/s. Used when all three are present; otherwise velocity is differentiated from the path. |

Names are matched case-insensitively, ignoring punctuation and a parenthesized unit, so
`Latitude`, `LAT` and `lat (deg)` are one column. `alt_ft` is recognized and `alt_m` is
deliberately *not*: a file naming meters fails to load rather than being read as feet.

**Where the frame sits.** `--origin <lat> <lon>` fixes the local x/y origin; without it, the
first position of the first geodetic track is used. `--frame_bearing` turns the axes so +x runs
along a compass bearing — the default leaves them east/north, which is what a map wants. One
frame serves the whole run, because every trajectory in it shares the observers.

**What is derived, and how well.** A WGS84 track gives position only, so the velocity and the
load factor behind C<sub>T</sub> come from the first and second derivatives of a smoothing
spline through the path. `--path_smoothing` says how far, in feet RMS **per axis**, that fit
may sit from the given positions. It is an acceptable deviation, not an abstract weight: the
solver finds the smoothest curve whose RMS residual meets the bound, so asking for 10 ft means
the fit will use 10 ft of slack, not merely stay under it.

The default is 10 ft. Zero interpolates exactly, and exact interpolation is the wrong thing for
both kinds of input: measured positions make the second derivative track coordinate rounding,
and sparse waypoints are a statement of intent rather than a constraint to be met to the inch —
flying precisely through them manufactures aggressive maneuvering, which is itself a noise
source.

Feeding the demo track's own lat/lon/alt columns back in, against the footprint computed from
its full inertial state:

| Input | Median | 99th pct | Max | Cells over 1 dB (of 14400) |
| --- | --- | --- | --- | --- |
| 50 Hz positions, default `--path_smoothing 10` | 0.073 dB | 0.72 | 1.42 | 25 |
| 50 Hz positions, `--path_smoothing 0.0656` (= 0.02 m) | 0.031 dB | 0.36 | 1.44 | 6 |
| 50 Hz positions, `--path_smoothing 0` (exact) | 8.0 dB | 9.57 | 10.58 | 14400 |
| 50 Hz positions and velocity columns, default | 0.015 dB | 0.36 | 1.42 | 4 |
| Waypoints 10 s apart, resampled to 0.1 s, default | 0.121 dB | 0.86 | 2.10 | 76 |

The projection itself is not the lossy part: it reproduces that file's own x/y to under a meter
over 3.2 km, and altitude passes through exactly. What costs is differentiating twice.

The default is not the best row in the table — a dense 50 Hz track does better at 0.0656 ft,
which is tuned to *that* file's coordinate resolution. It is the best *default*, because it is
the only setting that is good on every input shape at once: 0.073 dB on dense positions,
0.121 dB on sparse waypoints, 0.015 dB when velocities are given. Exact interpolation is
uniformly catastrophic on measured input, and a value tuned for one sampling rate is wrong for
another. Tune it down if you know your coordinate resolution; leave it alone otherwise.

**The ends of a dense track are fitted freely.** A natural smoothing spline is pinned to zero
second derivative at the first and last knot, which is a statement that the aircraft is not
accelerating there. A track that starts or stops mid-maneuver is therefore fitted with a
velocity and a load factor that bend toward a straight line over the last few seconds: on a
constant-acceleration test path (1.5 m/s²) the natural fit's acceleration is zero at the end,
0.50 four seconds in and 1.14 eight seconds in, where the freed fit holds it within 11%
throughout. For positions sampled at
0.5 s or finer, with a nonzero `--path_smoothing`, the fit is therefore made over the track
extended by 10 s at each end along a local quadratic (fitted to the last 3 s), so that the
natural condition lands on the extension and not on the data. The extension is not used
afterwards. Sparser tracks, and `--path_smoothing 0`, keep the natural ends; there the
interpolant is the statement of intent, and a quadratic fitted to a few waypoints is not. On
the table above this takes the default row from 0.109 dB median and 84 cells over 1 dB to the
values shown, and leaves the exact and waypoint rows as they were. The shipped AC tracks carry
their own state and never go through this.

Sparse waypoints are resampled onto the `-s` grid first — never coarsened — and their load
factor is an artifact of the interpolant, which is logged.

> The `--path_smoothing 0` row was, until recently, reported as 0.196 dB. That measured a bug:
> `SmoothingSpline` short-circuited at λ = 0 and returned zero curvature, so the load factor
> came out as exactly 1 g on every state and the rounding amplification never happened. With
> the interpolating spline carrying its own curvature, exact interpolation costs what the
> surrounding text always said it would.

**Note on velocity columns.** `--path_smoothing` is a distance and applies to the fitted
position only. A file carrying `ve`/`vn`/`vu` has those interpolated exactly; they are
measurements of the quantity wanted, not something to be re-derived.

**Output.** `--latlon` writes a companion beside every level file (`out.csv` →
`out_latlon.csv`) holding latitude, longitude and level, one row per observer, for dropping a
footprint onto a map. It needs a frame: either a geodetic track, which carries one, or
`--origin`. A grid regular in the local frame is not regular in latitude and longitude, so this
file is always scattered, never gridded.

**The transform.** `src/geodetic.hpp` is the American polyconic projection about the reference
point (it agrees with a transverse Mercator to 3 mm at 10 km, and reproduces the shipped
tracks' own x/y), with altitude carried through as z — a map projection with a height attached rather than a true ENU
frame. Over the few kilometers a footprint spans the difference is far below the modeling
error. `from_east_north_up` applies the frame's rotation without its translation, for
velocities; `to_geodetic` inverts the projection by Newton iteration, and is what `--latlon`
uses.

The shipped AC tracks are the case where the axes are not east/north: their x/y run along a
compass bearing of **210°**, the approach heading, not along the wandering heading in the
file's own column 9. Their ground plane sits 54.2651 ft below their altitude column. So
```sh
./build/niceops -t track.csv --origin 39.8104276 -74.4263086 --frame_bearing 210 \
                --ground_elevation 54.2651
```
is the demo case restated in WGS84, and is what the first row of the table above measures
(`--path_smoothing 0.0656` gives the second).

## Noise databases

A database (`.nod`, netCDF) holds one group per operating condition: advance ratio μ, flight
path angle γ and thrust coefficient C<sub>T</sub>. Each group stores one sphere around the
rotor: scattered (φ, θ) channels of the overall A-weighted level and of the EAA
excess-attenuation factor, and optionally per-band spectra. Databases are built by panam, from
measured or computed hemispheres. [How a prediction is made](#how-a-prediction-is-made)
describes how NICE-OPS reads them.

### Spectra

A database may carry per-band spectra alongside the overall A-weighted level: a
`frequency` dimension and an `amplitude(phi, theta, frequency)` array per sphere. When it
does, `--export_aam` can write an interpolated sphere back out in the AAM/RNM hemisphere
format the source spheres arrive in. Without spectra the flag refuses, since
`AMPLITUDE(PHI, THETA, FREQUENCY)` is the whole payload of that format.

Six of the databases in `databases/` carry spectra: the AS350, B407, Be206L3, EC130, RO-44 and
RO-66 `_fixed_load_factor` files. The S-76D demo database does not, so the demo runs the
broadband path. To build one with spectra, use panam's
`build_empirical_database(..., include_spectra=True)`, which now carries the per-band
amplitudes through from the source spheres. Band-limit it — `infreqs=(10.0, 3500.0)` —
because the top few bands of a measured sphere are usually noise-floor artifacts.

The AAM/RNM source spheres live in the group's OneDrive under `Data/Sphere_Database`, one
directory of `.nc` files per aircraft (AS350, Be206L3, Be407, Be430, EC130, MD-902, Mi-8, R-44,
R-66, SH-60B), each with its own `vehicle.cfg` and `vehicle.json`. `Be407_spectral.nod` sits
beside them — 26 bands, 10–3150 Hz, built from `Be407/` as above. The built `.nod` files in
`databases/` are tracked; the source spheres are far too large and stay in OneDrive.

### Where the 2017 databases come from

The AS350, B407, Be206L3, EC130, RO-44 and RO-66 databases are built from the 2017 Noise
Abatement flight test by panam (`noise_abatement_2017.build_all`, then
`build_empirical_database`). The validation harness keeps the recipes (`spheres/rebuild.py`,
`spheres/hover_sphere.py`, `spheres/release_database.py`, `spheres/rebuild_database.py`). Each
file records its build at the root: `speed_reference`, the build atmosphere, `hover_source` when
the hover was synthesized from a sphere outside the build, and `DOPPLER_SHIFT_REMOVED` and
`coverage` per group.

- **B407** is the reference aircraft and the one validated most. Its spheres are depropagated
  along curved rays through each run's own atmosphere, with a 2-degree emission floor, an
  adaptive 2-degree grid and tone-aware third-octave bands. A rotor harmonic near a band edge is
  filed whole in its own band instead of being split by whole FFT bins.
- **The B407 hover** is synthesized from the slowest level pass with the Doppler shift taken
  out, then corrected by the measured hovers. Two rings were measured: 18–25 degrees below the
  horizon (out of ground effect) and 0.5–3 degrees (in ground effect). Each ring is a Fourier
  series in azimuth from the heading, with a sum-to-zero offset per heading for the wind. The
  two are blended linearly in elevation from 3 to 17.5 degrees. Below the high ring each
  azimuthal order fades toward the nadir as cos^m.
  - Nothing measures 3–17 degrees or the region under the aircraft, so treat the hover there
    as an interpolation.
  - The low ring is the noisier of the two (5 dB held-out rms against 2.8).
- **Be206L3, EC130, RO-44 and RO-66** (rebuilt 2026-10-04) take the B407's flight recipe and
  its de-Dopplerized hover source, but no hover correction: only the B407's hovers have been
  measured. Their hover is the slowest level pass averaged fore to aft (Be206L3 run 277103,
  EC130 298295, RO-44 229395, RO-66 234520). No approach card is in any of them. Against each
  aircraft's approaches, which are the validation set (NICE-OPS minus the measured 10 dB-down
  SEL, station atmosphere, median / rms in dB):

  | Database | Spheres, September → now | Approaches | September | Now |
  | --- | --- | --- | --- | --- |
  | Be206L3 | 110 → 110 | day 279, 906 mic-runs | +0.28 / 1.48 | +0.24 / 1.50 |
  | EC130 | 38 → 112 | day 298, 1025 mic-runs | +0.19 / 3.17 | +0.70 / 2.70 |
  | RO-44 | 50 → 83 | day 229, 585 mic-runs | +0.31 / 1.55 | +0.12 / 1.27 |
  | RO-66 | 87 → 88 | day 236, 876 mic-runs | +0.28 / 1.06 | +0.86 / 1.57 |

  - **RO-66 reads worse, and the September file was right for the wrong reason.** All four
    September files carried near-silent cells: finite levels far below the aircraft (as low as
    −70 dBA beside an 84 dBA neighbor), mostly in the abeam column near the horizon, which the
    loader blends as real. The rebuild raises RO-66's predictions by 0.2 dB under the track and
    1.3 dB at 1,200–2,000 ft off it, where the aircraft is heard abeam, while the cells both
    files measured moved by +0.05 dB at the median. The holes had hidden the low-elevation
    over-prediction every aircraft shows. The rebuilt files have none (RO-66 had 62 such cells
    and 3,962 empty ones below the horizon, and now has 0 and 108).
  - **EC130's September spheres were depropagated at a standard day** (20 °C, 20 % RH,
    101.3 kPa), where the stations read −1 to 19 °C, 37–94 % and 88–89.5 kPa. The dataset spells
    the folder `EC130B4_Ground Stations`, and panam looks for `EC130B4_Ground_Stations`. They
    also lacked day 296, 62 of the spheres now. The day-298 approaches have no ambient run, so
    their measured SEL carries the background and has no SNR to filter on: read EC130's rms,
    not its median.
  - RO-44's September spheres lacked day 227, which holds 32 of its steady runs.
- **AS350** keeps the September recipe and spheres: straight rays, a 10-degree floor and a
  25-degree grid.

A database written before panam c53588f (2026-09-26) may have its last group unwritten. Every
entry is then the netCDF fill value (9.97e36), which the loader reads as a finite 1e37 dB
level. The September AS350, EC130, RO-44 and RO-66 files had this, and spectral runs that
reached the fastest +35-degree condition gave an infinite SEL.

### Speed reference

A database's flight conditions come from some velocity: μ = |v| / V_tip, the path angle
asin(v_z / |v|), and, unless `--nose` says otherwise, the sphere's azimuth reference, which is v
projected into the tip-path plane. Every database so far took v over the ground. `--speed_reference air` takes it
through the air instead: the ground velocity less the profile's wind at each state's height,
turned into the local frame by `--frame_bearing`. The load factor, and so CT, stays inertial,
because the rotor's thrust balances the aircraft's acceleration over the ground whatever
the air does.

- **The database says which it is,** in a root text attribute `speed_reference`, `ground`
  or `air`, and NICE-OPS follows it by default. A file without one is ground. panam writes
  the attribute (`build_empirical_database(..., speed_reference=...)`, default `ground`; in
  panam eeb040b and since), and the 2017 databases it builds record `ground`.
  Giving `--speed_reference` against the database is logged as an override.
- **Air needs a wind.** With no wind columns in `--atmosphere`, and no `--wind`, an
  air-referenced run is refused rather than taken as calm. `--speed_reference ground` says so
  explicitly. A `--wind` of zero is calm stated rather than assumed, and the airspeed is then
  the ground speed.
- **It is after `--trim`,** which cuts an approach at its slowest point over the ground.
- **Calm air changes nothing**, to the bit, which a CTest case checks.

What it moves, on the 2017 Noise Abatement approaches: 1–3% in μ and about 1° of nose on
calm days. On a 15 kt day, 5–10% in μ and 8–15° of nose. Matching a ground-built database with
an air-referenced trajectory mixes the two, so change both together.

### The sphere's nose

The sphere's azimuth is measured from its nose, so where the nose points decides which side of
the sphere every observer reads. `--nose` chooses it:

- **`velocity`** (the default): the velocity `--speed_reference` selects, projected into the
  tip-path plane. At hover the velocity is drift, so below 3 kt the nose is the heading, turning
  into the velocity by 8 kt on a smoothstep.
- **`ground_velocity`**: the ground velocity, whatever `--speed_reference` says, with the same
  turn to the heading at hover. Under `--speed_reference ground` it is `velocity`. It lets μ be
  taken through the air while the nose stays on the ground track: on the windy day of the 2017
  test (day 284, `profile_plate`, the held-out database) the along-track wind's share of
  `--speed_reference air`, which moves μ, took the SEL rms from 2.36 to 2.51 dB, and the
  cross-track share, which mostly turns the nose, took it to 2.25.
- **`heading`**: the state's compass heading, at every speed. In a crosswind a helicopter flies
  with its nose turned into the wind by the crab angle, about asin(crosswind / airspeed), while
  its track runs straight. The rotor's and the airframe's directivity turn with the airframe, not
  with the track.

**How the heading enters.** The local frame's angle (counterclockwise from +x) of compass heading
h is *b* − h, *b* being the frame's bearing, since compass bearings turn clockwise from north.
*b* is `--frame_bearing`, or 90 (east) for a run whose geodetic frame was settled without one.
The nose is that horizontal direction carrying the velocity's climb, v<sub>z</sub> over the
horizontal speed, so a heading along the track gives exactly the nose `velocity` does, climbs and
turns included; below 8 kt the climb fades out on the hover smoothstep. The sphere frame projects
it into the tip-path plane as it does the velocity. μ, the path angles and C<sub>T</sub> do not
change: they are the flight condition, and come from the velocity as before.

**The frame has to be known, and it is never inferred.** When a run has no frame, the hover nose
takes its offset from the track itself: the circular mean of heading against track over the
states above 20 kt. That fit absorbs the track's mean crab, so heading mode read through it would
point the nose back along the mean track and undo itself. So:

- **No frame** (an AC track without `--frame_bearing`, and no `--origin`) is refused, with the
  reason. The shipped AC tracks need `--frame_bearing 210`; each 2017 track needs its own (the
  validation harness derives it per run, 270 for every day-284 steady run).
- **No heading column** (a WGS84 track) is refused under `--nose heading`. When only the database
  implies heading mode (below), the nose falls back to the velocity, as if the aircraft flew
  with no crab, with a warning: there is nothing the run could be given instead.
- **A state without a heading** in a track that has the column is refused, naming its time.
- The log gives the mean crab (heading minus ground track, positive with the nose right of the
  track) and its rms over the states above 8 kt, and warns if the mean is past 45 degrees, which
  is a frame turned the wrong way rather than a crosswind.

**The hover nose takes the frame too.** With a frame known, the heading below 8 kt in `velocity`
and `ground_velocity` modes is turned by it, rather than by the fit; the fit is kept only for a
run with no frame. On a crabbed approach ending in hover the fit turned the hover nose by the
approach's mean crab: 2-7 degrees (median 4.6) over the day-286 approaches, had they ended in
hover. This changes only tracks that both fly
above 20 kt for 50 states and slow below 8 kt. No shipped track and no 2017 B407 flight track
goes under 13 kt, and the 2017 hovers have no fast states (they took the frame already), so no
golden or validation result moved. A static operation with `--origin` and no
`--frame_bearing` now reads its headings as compass bearings in the east/north frame, where it
oriented the sphere by its zero velocity before.

**The database can say which it was built along**, in a root text attribute `azimuth_reference`:
`track` (the velocity, which is how every database so far was built) or `heading`. Absent is
`track`. NICE-OPS defaults the nose to it: `heading` for `heading`, `velocity` otherwise. An
explicit `--nose` wins, with a warning when the file records the other; a sphere read along a
nose it was not built along is turned by the crab. `--export_sphere` carries the attribute on
when the source records it. A value other than `track` or `heading` is refused.

On the panam side, the builder's integration branch of 2026-10-05 (not yet published) writes the
attribute on every database: `track` by default, `heading` for spheres built with its
heading-frame option (`build_sphere(azimuth_reference='heading')`, `--azimuth-reference`), and
`track` for a NORAH2 import; it refuses to mix the two in one database. Every shipped database
predates it and records nothing, which reads as `track`. Once a database records `heading`, an
AC-track run needs `--frame_bearing` (or an explicit `--nose velocity`), and a WGS84 track, which
has no heading column, falls back to the velocity with a warning.

**What it does, measured.** On the 2017 test's windy day (day 284 steady runs, 16 held-out runs,
729 microphones, the `min2_no284` database built without them) every run flew west with the wind
6-10 kt from its right, and crabbed 4 to 19 degrees into it, 12 on average, the most at the
slowest speeds. With
`--nose heading` against the same database, which is labeled along the track, and the paper's
inputs otherwise:

| Configuration | SEL rms, `velocity` | SEL rms, `heading` | Change (95% run bootstrap) | 5–10° elevation |
| --- | --- | --- | --- | --- |
| `stratified_plate`, day 284 | 1.69 dB | 1.54 dB | −0.15 [−0.28, −0.02] | 3.02 → 2.26 |
| `station_plate`, day 284 | 2.35 dB | 2.19 dB | −0.17 [−0.28, −0.03] | 5.76 → 5.22 |
| `station_plate`, day 286 (calm) | 1.26 dB | 1.24 dB | −0.02 [−0.05, +0.01] | 1.83 → 1.69 |

LAF<sub>max</sub> rms moves by −0.14 dB on day 284 (`stratified_plate`) and −0.01 on day 286.
Placebo runs, which turned every run's nose by a fixed 13 degrees rather than by its own crab,
show it is each run's own crab that helps: on day 284 the fixed turn in the crab's direction
gave −0.03 [−0.22, +0.16] dB, the opposite turn +0.47, and on the calm day it cost +0.29. The option was first emulated by turning
each track's velocity onto its heading; the option reproduces that to 0.01 dB in rms, and up to
0.5 dB at a microphone, because turning the velocity also moved the path angle relative to the
tip-path plane by up to 1.4 degrees in a turn, where the option leaves the flight condition
alone. The runs are the validation harness's paper cases (`paper/fidelity_study.py`), each
given `--nose heading` and its own `--frame_bearing`, scored as `paper/results.py` scores them,
with the bootstrap resampling whole runs.

**How it is checked.** Unit tests pin the crab's sign against a synthetic crosswind whose crab
is known (and against `CompassToLocal`, which turns the wind, so the two conventions are one),
the heading-along-track nose against the velocity's sphere frame through a tilted tip-path plane,
the climb's fade at hover, and the frame's precedence over the fit (and a CTest, two tracks
ending in the same hover after approaches that differ only in their crab). CTests run the demo track
with its heading replaced by its own ground track (heading mode reproduces the default to 1e-9),
and turned 15 degrees (the footprint moves); the database attribute selecting heading mode by
itself; each refusal by its message; the windowed peak search staying exact under heading noses;
and `ground_velocity` differing from the air nose in a wind.

## Propagation

### Broadband and spectral paths

The spectral path, which is the default whenever the database carries spectra (`--no-spectral`
turns it off), carries the spectrum itself to the observer instead of a single number. Each
band is spread by the same `20 log10(r/r_sphere)` and absorbed at its own ISO 9613-1
coefficient over `r - r_sphere`, and only then are the bands A-weighted and summed into a
level. The two paths answer different questions and will not agree:

| | `--no-spectral` | spectral (the default with spectra) |
| --- | --- | --- |
| Attenuation | database's `EAA`, one dB/km factor | ISO 9613-1 α(f), per band |
| Fitted for | the spectrum the source happened to have | the atmosphere you state |
| Atmosphere | none — `EAA` is baked in | temperature, pressure, humidity |
| Needs | any database | one carrying spectra |

The absorption model is ported from panam's `panam_acoustics/iso_9613_1_1993.py` and pinned
against it in `tests/test_atmosphere.cpp`; the A-weighting is panam's `dBAw`, pinned in
`tests/test_weighting.cpp`.

**Band-limit deliberately.** Summing whatever bands a database happens to carry quietly
makes the answer a property of that database's bandwidth, so two databases are not
comparable. `--band_limits low high` states the range; bands outside it contribute exactly
nothing, and the run logs how many of the database's bands survived. A limit that excludes
every band is an error rather than silence.

The atmosphere comes from the vehicle JSON's `Atmosphere` section:

```json
"Atmosphere": {
    "density": "1.225",
    "temperature": "300.71",
    "pressure": "105.74",        // optional; else derived as rho*R*T
    "relative humidity": "20"    // optional; a spectral run is refused without one
}
```

Only `density` is required, as before. Pressure defaults to the one implied by that file's
own density and temperature rather than to a standard day — `S-76D.json` pairs 1.225 kg/m³
with 300.71 K, which is not an ISA point, so assuming 101.325 kPa would contradict its own
numbers. `--humidity`, `--temperature` and `--pressure` override all of it.

### Stratified atmosphere

`--atmosphere file.csv` replaces that single atmosphere with a profile, and integrates each
band's absorption through it. Along a straight path, what a band loses is the ISO 9613-1
coefficient's mean over the heights the path spans, times the path length. That is the
fine-layered form of what certification does with 30 m layers (14 CFR 36 Appendix A), and a
unit test checks the two agree.

On the demo track with the B407 database, `atmospheres/isa.csv` moves 90% of the footprint by
-0.10 to -0.005 dB (mean -0.04 dB, extremes -0.57 and +1.11 dB), against a homogeneous run at
ISA sea-level conditions. On a grid it costs nothing measurable, because every observer stands
at one height and each state's mean is formed once. Scattered observers at different heights
pay a subtraction and a multiply per band per evaluation.

Columns are found by name, in any order, one row per level; rows may run top down:

| Quantity | Column names | Unit |
| --- | --- | --- |
| Height above the observers | `z`, `agl` | `--units`; `z_ft` / `z_m` must agree with it |
| or altitude above sea level | `alt`, `altitude`, `height` | `--units`, less `--ground_elevation` |
| Temperature | `T_K`, `T_C` | named; a bare `T` is told from the values |
| Relative humidity | `RH`, `humidity` | percent |
| Pressure (optional) | `p_kPa`, `p_hPa`, `p_Pa` | named; a bare `p` is told from the values |
| Wind (optional) | `wind_speed_mps` or `wind_speed_kt`, with `wind_direction` | direction is the compass bearing it blows **from**, degrees |
| or wind components | `u_mps` and `v_mps`, or `u_kt` and `v_kt` | toward east and toward north |

```csv
z_ft,T_K,RH,p_kPa
0,288.150,20,101.325
1000,286.169,20,97.717
2000,284.188,20,94.213
```

Without a pressure column, pressure is integrated hydrostatically up from the surface
pressure, which is the vehicle file's or `--pressure`'s. A single row is a uniform
atmosphere. `atmospheres/isa.csv` is the ISA to 10,000 ft at 20% RH, and
`atmospheres/B407_uniform.csv` holds `B407.json`'s own atmosphere. With the latter a run
reproduces the homogeneous one, which a CTest case checks.

The rules that keep it from guessing:

- **The profile has to reach the aircraft.** A trajectory climbing above the top level is
  refused, not extrapolated. Below the lowest level, that level is held down to the ground.
- **Surface properties stay at the surface.** The profile's level at the ground supplies the
  sound speed for the reception-time integral and for the ground impedance.
- **`--temperature` and `--humidity` cannot be combined with it.** Nor can `--pressure`
  unless the file gives no pressure, in which case it is the surface pressure.
- **Only the spectral path uses it.** The broadband path's `EAA` is the database's own
  fitted attenuation, and there is nothing per band to integrate.
- **Wind columns must name their unit.** Nothing in the values tells m/s from knots. Wind
  matters only to `--refraction` and to `--speed_reference air`; absorption does not depend
  on it. The log names every column it ignored.

### A profile from the surface wind

A measured profile — LIDAR, a balloon, a met tower — is rare. A surface wind is not: an
airport's METAR, a weather station beside the microphones. `--wind SPEED,DIRECTION[,HEIGHT]`
builds the profile from that one observation, in memory, in place of an `--atmosphere` file.
SELENE does the same from a surface wind and a stability class. It is the measured route to
refraction on a windy day without a sounding.

```sh
niceops ... --humidity 25 --wind 12kt,270,6.56 --refraction stratified --frame_bearing 210
```

**The conventions are the `--atmosphere` file's.**

- The speed names its unit, as a wind column has to: `12kt` (or `kts`, `knots`), `6.2mps`
  (`m/s`), `14mph`. A bare number is refused, since nothing in it tells m/s from knots, and it is
  not taken from `--units`.
- The direction is the compass bearing, in degrees from true north, that the wind blows
  **from** — the meteorological convention, the `wind_direction` column's. A METAR's 27012KT is
  `12kt,270`. Like every wind it is turned into the local frame by `--frame_bearing`, so the
  shipped AC tracks need `--frame_bearing 210`.
- The anemometer height is a length in `--units`, like every other length. Omitted, it is
  10 m, the WMO standard exposure for surface wind.

**The wind** follows the logarithmic law of the neutral surface layer,

u(z) = u(z_ref) ln(z / z₀) / ln(z_ref / z₀),

through the observation itself, up to `--surface_layer_depth` (default 200 m), and is held at
its value there above that. The neutral surface layer is about a tenth of a boundary layer one to
two kilometers deep; above it the wind turns and grows slowly toward the geostrophic wind, which
a single surface observation says nothing about. The direction is the same at every height: an
Ekman veer added nothing on the 2017 windy day. z₀ is `--roughness_length`, default 0.03 m: "open
flat terrain; grass, few isolated obstacles" in the Davenport–Wieringa classification (Wieringa
1992; WMO Guide to Instruments and Methods of Observation, WMO-No. 8, Annex 5.B), the class of an
airfield or of short desert grass. Mud flats and snow are 0.005 m, low crops 0.1 m, farmland
with hedges 0.25 m, suburbs 1 m.

**The air.** The surface temperature, humidity and pressure are the run's own: the vehicle file's,
or `--temperature`, `--humidity` and `--pressure`, which may be combined with `--wind` (they set
the surface air) where they may not with `--atmosphere`. As without a profile, a run whose
humidity nobody gave is refused. Above the ground:

- the temperature lapses at the dry-adiabatic rate, −9.8 K/km, which is what makes the layer
  neutral;
- the pressure is integrated hydrostatically through that temperature, as for a file without a
  pressure column;
- the water vapor is held: its molar concentration, which is what ISO 9613-1's absorption takes,
  keeps its surface value, as it does in a well-mixed layer. The relative humidity therefore rises
  as the air cools, about half a percent of itself per 10 m, and is held at 100 percent if the air
  saturates (the run warns).

**The levels.** The profile is linear between its levels, so they decide how faithfully it
follows the law. The ground, z = 0, has no wind: the law's limit. The next level is ten
roughness lengths up (0.3 m by default; the anemometer, if that is lower), and from there they
are spaced geometrically, `--profile_levels` (default 16) per factor of ten, to a top 5% above the
highest state or microphone, and never below the top of the surface layer. The anemometer, the
top of the surface layer and the top are levels themselves, so the profile passes through the
observation exactly. Above the top nothing is lost: the wind is held and the air lapses, so the
effective sound speed falls with height in every direction and no ray turns back from there.

**No level above the ground is ever calm.** That is deliberate. `--refraction stratified` reads a
calm level just above the ground (rows at 0 and 0.1 ft, both still) as a thin layer that puts the
levels 22–28 dB low. The law is zero at z₀, so a profile sampled there would hold such a level;
starting at ten roughness lengths keeps every level above the ground moving. A calm observation
(`0kt`) gives a profile with no wind at all, the temperature lapse alone, and a database referred
to the air then takes the ground speed as the airspeed.

**`--write_atmosphere FILE`** writes the profile each trajectory used, in the `--atmosphere`
format (z in `--units`, `T_K`, `RH`, `p_kPa`, `u_mps`, `v_mps`, to 17 significant digits). Fed back
through `--atmosphere`, it reproduces the run, which a CTest case holds to 1e-9. With several
trajectories each gets its own, `<track>_atmosphere.csv`, since each has its own surface air and
top. The one profile that does not replay is completely dry air (`--humidity 0`): `--atmosphere`
refuses a file whose humidity never exceeds 1 percent, taking it for a fraction, and the run warns
when it writes one.

**How it was checked.**

- *The law, the air and the conventions*, in `tests/test_atmosphere_profile.cpp`: the profile
  passes through the observation and is zero at z₀; the temperature, pressure and molar
  concentration against their closed forms; a wind from 270° blows toward the east and equals
  what `--atmosphere` makes of the same wind written as speed and direction; the speed's units;
  no calm level above the ground for any class, roughness or anemometer height; and the written
  profile reads back to the bit (heights to an ULP of the unit conversion).
- *Through the command line*, a run in feet that sets every shaping option and the three surface
  overrides writes its profile, and a CTest case compares it with the same profile computed
  independently from the closed forms (`tests/data/wind_profile_expected.py`): 1e-12 absolute at
  every level. That is what catches a length left in feet, an override not applied, or
  `--profile_levels` not passed on, none of which the replay can see.
- *Discretization.* Refining the levels from 16 to 64 a decade moves the plate microphones of the
  Coyle array under stratified refraction by at most 0.016 dB (0.053 from 8), which a CTest case
  holds to 0.1 dB. On the measured day-284 runs below it moves 729 microphone-runs by 0.006 dB
  rms, 0.105 at most. That last tenth is not the profile: on a point-receiver grid, cells 1500 m
  out downwind move by ±0.2 dB with any change of layering, without settling from 8 levels a
  decade to 100. It is the stratified model's arrival search, whose samples follow the layers.
  `--profile_levels` does not refine the bottom layer, from the ground to ten roughness lengths,
  where the profile is linear and the law logarithmic. Filling it with the law (levels at 1.5, 2,
  3, 5 and 7 z₀, fed through `--atmosphere`) moves the day-284 results by 0.014 dB rms and 0.15 dB
  at most; dropping every level below 1 m moves them 0.046 dB rms, 0.30 at most. Levels closer to
  z₀ would carry almost no wind, which is the thin layer the stratified model cannot take.
- *Against the 2017 measurements* (B407, plate microphones, `--refraction stratified`; SEL rms
  over the paper's usable microphone-runs). Each run's wind is the mean of the ground stations'
  anemometers (SWS 1, 2 and 4) over ±60 s of the run's start, its magnetic direction turned by the
  13.6° declination, taken as 2 m high. The data report (NASA/TM-2019-220264, §8.6.3–8.6.4)
  gives the station speeds in knots, as it does the LIDAR's, which the validation harness has
  since confirmed independently for the LIDAR (its `docs/lidar_units.md`: balloon, station cups
  and the aircraft's own turns); the files themselves carry no unit and the anemometer height is
  not recorded, so knots and 2 m are the reading, not a measurement. mph is given for the
  sensitivity.

  | | day 284, windy (16 steady runs, held-out database) | day 286, calm (25 approaches) |
  | --- | --- | --- |
  | straight rays, station air | 2.35 dB | 1.26 dB |
  | measured profile (LIDAR and balloon) | 1.69 dB | 1.35 dB |
  | `--wind`, station speed read as knots | **1.62 dB** | **1.25 dB** |
  | `--wind`, station speed read as mph | **1.60 dB** | **1.25 dB** |
  | `--wind`, knots, z₀ = 0.09 m | 1.71 dB | |

  The station winds were 4–9 kt from 264–293° on day 284 and under 1.1 kt on day 286. On the windy
  day the surface wind alone recovers all of what the measured profile buys over straight rays,
  and a little more; the gain is in the low-elevation microphones, 5.76 → 2.65 dB rms at 5–10°.
  On the calm day it does no harm. (The 284 sample is 729 microphone-runs, the 286 one 1111.)

**Stability classes (unvalidated).** `--stability` takes the Klug–Manier dispersion classes of TA
Luft (1986), or the Pasquill letters that correspond to them, as SELENE does. Each replaces the
log law with a power law through the observation, u ~ z^m, and the neutral lapse in the surface
layer with a temperature gradient: the middle of the class's range in US NRC Regulatory Guide
1.23's ΔT classification. Above the surface layer every class lapses dry-adiabatically.

| `--stability` | Pasquill | character | m (TA Luft) | dT/dz in the surface layer |
| --- | --- | --- | --- | --- |
| `I` | `F` | very stable | 0.42 | +2.75 K/100 m |
| `II` | `E` | stable | 0.35 | +0.5 K/100 m |
| `III/1` | `D` | neutral to stable | 0.22 | −1.0 K/100 m |
| `III/2` | `C` | neutral to unstable | 0.20 | −1.6 K/100 m |
| `IV` | `B` | unstable | 0.15 | −1.8 K/100 m |
| `V` | `A` | very unstable | 0.09 | −2.0 K/100 m |

None of them has been shown to help, and a run that takes one warns. On day 284, class II scored
2.22 dB, nearly back to straight rays (its inversion and shear deepen the upwind shadow: −3.3 dB
median at 5–10°), and IV 1.57; an earlier trial of II's exponent with a +3 K/100 m inversion,
driven by the LIDAR's 46 ft wind, scored 1.90 where the log law through that wind scored 1.72. One day, scored after the fact, cannot choose a class, and nothing on the day says
which it was. Use `neutral` unless the class is known from measurements.

### Ground effect

`--ground_effect` adds the interference between the direct ray and the one reflected off the
ground — Chien–Soroka over a Delany–Bazley impedance, with Chessell's boundary-loss correction
and band averaging, ported from panam's `ega()` and pinned against it in
`tests/test_ground_effect.cpp`.

**It is a gain, not a loss**, despite the name: over rigid ground with the receiver on the
surface it tends to +6 dB, the pressure doubling of a perfect reflector. It goes sharply
negative only in the interference dip, where the reflected path arrives out of phase. On the
Be407 demo case:

| | mean | max |
| --- | --- | --- |
| spectral, no ground | 65.10 | 96.38 |
| + concrete (10000) | 68.45 (**+3.35**) | 99.00 |
| + grass (225) | 64.73 (**−0.36**) | 97.89 |

Hard ground is 3.7 dB louder than soft, and grass roughly cancels its own reflection gain
against the dip. `--flow_resistance` spans five orders of magnitude and the answer is sensitive
to it, so name the surface rather than guessing.

`--turbulence` decorrelates the two paths and fills the nulls in; leaving it at zero gives
unrealistically deep interference dips.

The impedance, turbulence and roughness models follow Attenborough & Taherzadeh, *Prediction
of outdoor ground effect*, Applied Acoustics 242 (2026) 111067, as panam's `ground_plane`
module does. They are pinned against `ground_plane.pole_level` in
`tests/test_ground_effect.cpp`. Delany–Bazley is poor for grassland. On the 2017 Amedee pole
microphones, variable porosity with σe ≈ 200 fitted best, with HARMONOISE γT ≈ 3e-5 and
10–15 mm roughness (panam `docs/ground_plane_corrections.md`).

### Receivers

Each observer's level is evaluated at its microphone: the observer's point raised by the
microphone's height above the ground. What the microphone is decides the ground model:

| Microphone | Ground model | How to ask for it |
| --- | --- | --- |
| flush with the ground, height 0 | the plate (below) | `--plate_table` |
| on a pole, height h | the point receiver at h | `--ground_effect`, with h from `--receiver_height` or the observer file |
| anything, without `--ground_effect` | none: the free field at the microphone | |

With a four-column observer file, every observer takes the run's receiver: all on plates
with `--plate_table`, otherwise all at `--receiver_height`. A fifth column gives each
microphone's own height, in `--units`, so one run can hold an array of flush microphones with
a few on poles, as the test campaigns do. With `--plate_table`, height 0 is a plate and
anything else is a point receiver at that height:

```csv
id,x_ft,y_ft,z_ft,height_ft
2,-788.3,-1085.1,-1.1,0
10,-405.9,1940.8,-10.1,3.937
```

The log reports how many observers took each receiver. `tests/data/mixed_receivers.csv` is a
worked example, and a CTest case checks that it gives the same levels as separate runs over
each kind of microphone.

Before the microphone height set the evaluation point, the point receiver's direct path ended
on the ground even though its ground model put the microphone 1.2 m up. Correcting that moved
point-receiver SEL on the demo grid by up to about ±0.8 dB (median 0). Almost all of it comes
from the low states at the end of the approach, where 1.2 m is a sizeable fraction of the
source's height. Plate and free-field results did not change.

### Ground-plane microphones

A flush ground-plane microphone sits on a plate: rigid close in, soft ground beyond. The point
receiver can model one surface or the other, not both. `--plate_table` models the plate
itself, using panam's axisymmetric boundary-element solution for a GRAS GR1425 plate (400 mm,
8 mm thick, lying on the ground) with the microphone at 3/4 R per SAE ARP 4055. The solution is
stored as two complex transfer functions of the direct and image waves:

    p / p_free = P_d + Q P_r

P_d and P_r are tabulated over frequency (five sub-frequencies per band), elevation and
azimuth. Q is the ground's spherical-wave reflection coefficient, which NICE-OPS computes for
each geometry. The table depends on the ground, so it is built for one ground and NICE-OPS
refuses any other. Two tables are shipped in `plates/`:

| table | ground |
| --- | --- |
| `gr1425_db225.nc` | `--flow_resistance 225` (Delany–Bazley grass, the default) |
| `gr1425_vp200.nc` | `--impedance_model variable_porosity --flow_resistance 200` (the fitted Amedee ground) |

Build another with panam's `axisymmetric_bem.table(bands, …)` and `write_netcdf`. A table at
the databases' 31 bands takes about 20 s on 14 cores. Turbulence and roughness do not apply to
the plate, because its reflection is millimeters from its direct path.

Q is taken as quadratic in log frequency across each band, through the band's edges and
center. This is within 0.07 dB of evaluating Q at every sub-frequency wherever the level is
above −15 dB. The exact per-sub-frequency form (`PlateReceiver::level_exact`) matches panam's
`board_level` to 1e-6 dB in `tests/test_plate_receiver.cpp`. Its cost on the full demo grid
(S-76D track, B407 fixed-load-factor database, 94.3 M evaluations, 14 threads):

| | acoustic loop |
| --- | --- |
| spectral, no ground | 1.6 s |
| `--ground_effect` | 2.1 s (7.1 s with `--exact_ground`) |
| `--ground_effect --plate_table` | 5.4 s (13.3 s with `--exact_ground`) |

By default both receivers interpolate per-trajectory tables instead of evaluating the model at
every observer and state. The tables are built exactly at their nodes in well under 0.1 s.
- **Point receiver:** the whole ground correction, over source elevation and log range. Level
  histories agree with the exact model within 0.001 dB.
- **Plate receiver:** its reflection coefficients. SEL agrees within 0.0005 dB rms.

The interpolated values never exceed the model's bounds, so the peak search and `--db_cutoff`
stay exact. `--exact_ground` evaluates the model directly. Either way the Faddeeva function uses
fast paths: an asymptotic series for |z| ≥ 6 and a Taylor table below, which agree with
Weideman's series within 1e-11.

The peak search and `--db_cutoff` bound the plate with `PlateReceiver::gain_bounds`, so both
keep working. The `plate_windowed_is_exact` ctest checks that skipped states change nothing.

### Terrain

Without `--terrain` the ground is the plane z = 0 everywhere: the ground model and the plate take
the source's height above it and the horizontal distance, whatever height an observer stands
at. `--terrain FILE` gives the ground a shape. Two things change:

- **Grid observers stand on the terrain**: each one's z is the ground's there, and its
  microphone is `--receiver_height` (or the plate) above that.
- **Each source–receiver pair's ground geometry is taken over its mean ground plane**: the
  least-squares line through the terrain heights along the path's ground projection, from the
  observer to the source's footprint. This is the "mean plane" of CNOSSOS-EU (Directive (EU)
  2015/996, Annex II, 2.5.6), the construction NORAH2 and ECAC Doc 32 use over terrain.

Spreading, absorption and the sphere lookup still take the real three-dimensional path from the
source to the microphone. Only the ground term changes.

**The geometry.** In the vertical plane through the path, let the mean ground line be
z = a + b s, with s the horizontal distance from the observer. The source is at distance D and
height z<sub>s</sub>, and the microphone is over the observer at z<sub>m</sub>. Measured
perpendicular to the line:

    h_s = (z_s - (a + b D)) / sqrt(1 + b^2)
    h_r = (z_m - a) / sqrt(1 + b^2)
    d   = (D + b (z_s - z_m)) / sqrt(1 + b^2)

This is a rotation of the plane, so the source–microphone distance is unchanged. Ground that
rises toward the source lowers h<sub>s</sub> and lengthens d, which lowers the grazing angle.
Ground that falls away does the opposite. For a nearly overhead source over ground falling toward
it, d comes out negative: the source's foot on the line lies behind the microphone's. The ground
model takes its size, which keeps its direct path equal to the real range, and a plate is crossed
from the other side. What each receiver takes:

| Receiver | Ground model's geometry |
| --- | --- |
| point receiver (pole) | source height h<sub>s</sub>, receiver height h<sub>r</sub>, distance d. The microphone's height above each path's plane is its own, so the model is evaluated exactly for every pair (see cost below). |
| flush plate | the plate lies along the plane under its microphone: source height h<sub>s</sub> − h<sub>r</sub>, distance d. |
| no ground model | nothing changes but where the observers stand. |

A height that comes out below the plane is set to zero, as CNOSSOS-EU does with a negative
equivalent height. This happens for a receiver in a dip, and for a source behind rising ground.
The second case is screening, which is not modeled (below). The log counts how many evaluations
were affected.

On level ground at z = 0 every one of these is the plain geometry, bit for bit. A level terrain
file therefore reproduces the run without it, which is the first thing the tests check.

**Datums.** The local frame's z = 0 is `--ground_elevation` on the trajectory's own vertical
datum, as for a WGS84 track (see [Geodetic coordinates](#geodetic-coordinates)). The ground's z
is

    z_ground(x, y) = H(latitude, longitude) + --terrain_offset - --ground_elevation

where H is the grid's height, interpolated bilinearly, and the position is found by inverting
the run's frame projection. So the run needs a frame: `--origin`, and `--frame_bearing` for
rotated axes (or a WGS84 track, which brings its own frame). NICE-OPS does not convert between
datums. `--terrain_offset` is how you state the difference, and the log repeats it.

*The 2017 Amedee data* are the case to know. The 3DEP grid is NAVD88. The flight altitudes and
the surveyed microphone heights are ellipsoidal (WGS84): the microphone list's `Hgt(m)` sits
23.6 m below 3DEP at every microphone (standard deviation 0.03 m). The difference is the local
geoid height. The AC tracks' `alt` less `z` puts their z = 0 at 3926.7 ft ellipsoidal, which is
the reference point's `ref_elips_ft`. A run there is:

```sh
niceops -t 286264AC.csv -b 286264_observers.csv ... \
    --origin <ref_lat> <ref_lon> --frame_bearing <bearing of the track's +x> \
    --ground_elevation 3926.74 --terrain amedee_3dep.nc --terrain_offset -77.43
```

The validation harness's `terrain/dem_to_netcdf.py` crops a 3DEP GeoTIFF into this format, and
`terrain/README.md` there documents the datum finding.

**What the run checks**, since a datum or frame mistake gives a plausible-looking wrong answer:

- The grid must cover every observer and every state's footprint, with a margin, and have no
  missing node there. Otherwise the run is refused, and the message names the extent it needs.
  Heights must be in meters. A grid that declares another unit is refused.
- An AC track records each state's latitude, longitude and altitude beside its x/y/z. The run's
  frame must put the latitudes and longitudes within 30 m rms of the x/y, or it is refused,
  because the terrain would be read under the wrong place. Above 3 m it warns. The demo track
  agrees to 0.85 ft, and the 2017 tracks to 1–5 ft. Where the track's alt − z is one number
  (within 1 m), `--ground_elevation` must match it to 1 m.
- A state more than 10 m below the terrain under it is refused as a datum mismatch: leaving out
  `--ground_elevation` at Amedee puts the ground 1197 m above every state. A state up to 10 m
  below is warned about, since that is within the terrain's own error (3DEP's is 1–2 m) or a
  landing.
- Scattered observers (`-b`) keep their file's z. Where it disagrees with the terrain by more
  than 1 m, the run warns and gives the median and worst disagreement, since a uniform offset is
  a datum mistake. `--observers_on_terrain` stands them on the terrain instead. At Amedee the
  surveyed microphones sit 0.06–0.14 ft from 3DEP in the median and 0.84 ft at worst.

**What is not modeled.**

- Terrain does not screen or diffract. A ridge between source and receiver is a source below the
  mean plane, and its height is set to zero.
- Each path's ground is one plane. Reflections off slopes away from that plane, and the
  terrain's roughness about it, are left out.
- **`--refraction` is refused with `--terrain`.** The ray models reflect off, and cast their
  shadows onto, the flat ground z = 0. Handing them a tilted plane per path would be a different
  model, not this one. `--propagation_table` is refused for the same reason: its source and
  receiver heights are over flat ground.
- The broadband path (`--no-spectral`) has no ground model, so there `--terrain` only sets where
  the observers stand.

For an elevated source, the specular reflection point is close to the receiver, at
d·h<sub>r</sub>/(h<sub>s</sub> + h<sub>r</sub>). A plane fitted over the whole path also weighs
terrain that the reflection never sees. That is the standard practice, and on smooth ground the
two agree. A receiver's h<sub>r</sub> above the plane follows the ground at its foot relative to
the plane, so it carries the grid's own decimeter-scale error. That error matters more for a
point receiver than for a plate (see the Amedee figures below).

**How the plane is found, and what it costs.** Fitting each pair's profile directly would mean
about 300 terrain samples for a 3 km path, 94 million times on the demo grid. Instead:

1. Per trajectory, the grid is resampled onto the local frame every `--terrain_step` (default:
   the grid's own spacing, at least 5 m). Each later read is then bilinear in local meters, with
   no inverse projection.
2. Per observer, profiles run along 1440 fixed bearings (every 0.25°, 13 m apart at 3 km), each
   sampled every step. They are kept as running sums of z and k·z, from which the least-squares
   line over any whole number of steps costs a few flops. A bearing is sampled only when a pair
   first needs it, and only as far as that pair needs.
3. A path between two bearings, or ending between two samples, takes the line interpolated
   linearly in bearing and in length between the four fits around it. That keeps it continuous.
   A path shorter than one step takes the fit over the first step along its bearing.

Along a bearing, at a whole number of steps, this is the direct least-squares fit
(`tests/test_terrain.cpp`, to 1e-9 in slope). On a uniform plane it is exact up to the cosine's
curvature between bearings: a 0.05 slope errs by 1e-7, or 0.4 mm at 3 km. On rolling synthetic
terrain with ±15 m hills it stays within 0.15 m of the direct fit at both ends of every path
tried.

Cost on the full demo grid (S-76D track, B407 database, 14,400 observers × 6,551 states,
14 threads, `--humidity 20`):

| | CPU, level ground (no `--terrain`) | level `--terrain` | hilly `--terrain` |
| --- | --- | --- | --- |
| `--ground_effect`, 1.2 m point receivers | 32 s | 35 s | 127 s |
| `--plate_table` | 71 s | 75 s | 75 s |

The planes cost about 3 s of CPU, roughly 35 ns a pair. Point receivers over real terrain cost 4x,
because each pair's receiver height differs from the one height the `GroundTable` was built for,
so the model is evaluated exactly (as with `--exact_ground`). On level terrain the height is the
pole's and the table is read. Plates keep their table either way. For the 2017 arrays (about
40 microphones) none of this is noticeable.

**How it was checked.**

- *Level terrain reproduces the run without it, byte for byte*: point receivers, plates, and
  scattered observers stood on it (`terrain_flat_*_is_flat`).
- *A uniformly tilted plane gives the closed-form h<sub>s</sub>, h<sub>r</sub> and d*, nearly
  overhead included, for point receivers and plates. The rest of the path record stays the
  straight line's. The plane's own geometry is a rotation, and is hand-checked for rising and
  falling ground (`tests/test_terrain.cpp`).
- *The hills runs are pinned* (`terrain_hills_point_golden`, `terrain_hills_plate_golden`, at
  1e-6), because the windowed and marched runs come from one binary: treating the plate as a point
  receiver at its call site moved cells by 0.01 dB and passed every other terrain test.
- *The peak search stays exact.* The ground's gain bounds (`ground_gain_bounds`,
  `PlateReceiver::gain_bounds`) hold for any geometry and any receiver height, so the windowed run
  over hills equals the march (`terrain_hills_point_windowed_is_exact`,
  `terrain_hills_plate_windowed_is_exact`).
- *The terrain reaches the ground model, not just the observers.* Microphones stood on the hills
  through an observer file, without `--terrain`, give different levels, and the same pair with
  no ground model agrees (`terrain_reaches_the_ground_model`,
  `terrain_without_ground_places_observers`). A level that reads the `GroundTable` at the wrong
  height is caught in `tests/test_spectral_propagation.cpp`.
- *Every refusal above has a CTest* that passes on its message.

**What it changes.** At Amedee the ground is flat to within about a meter over the test range.
Along every track-to-microphone path within 3 km, a mean plane departs from the ground by at
most 0.3 m (harness `terrain/README.md`). With the 3DEP grid, on the paper's configuration
(station air, flush plates, the shipped B407 database for day 286 and the held-out one for
day 284), with against without `--terrain`:

| runs | microphone-runs | mean ΔSEL | largest \|ΔSEL\| |
| --- | --- | --- | --- |
| 286264, 286265 (approaches) | 92 | −0.003 dB | 0.026 dB |
| 284153, 284197 (steady) | 94 | −0.006 dB | 0.023 dB |
| all four | 186 | −0.004 dB (rms 0.007) | 0.026 dB |

So terrain does not explain the paper's residuals at Amedee. 1.2 m pole microphones move more:
on 286264 and 284153, up to 0.11 dB with `--observers_on_terrain` and 0.19 dB with the surveyed z. Their height above
each plane follows the decimeter-scale relief and the grid's error at their feet.

At a sloped site it matters. A 1–3° slope along the path changes the grazing angle by the slope.
Through the plate table and the B407 spectrum, that is about 1 dB at 5° peak elevation and 0.3–0.7
dB at 10°, for flush plates. At 1.2 m it is tenths of a dB:

| peak elevation | plate, slope −3° / −1° / +1° / +3° | 1.2 m pole, 500 m, the same slopes |
| --- | --- | --- |
| 5° | −5.0 / −1.05 / +0.77 / +1.83 dB | −0.52 / +0.16 / 0.00 / +0.23 dB |
| 7.5° | −2.1 / −0.53 / +0.42 / +1.06 | −0.11 / −0.12 / +0.15 / +0.56 |
| 10° | −1.2 / −0.32 / +0.27 / +0.69 | −0.49 / −0.20 / +0.23 / +0.50 |
| 20° | −0.29 / −0.09 / +0.08 / +0.21 | −0.43 / −0.04 / −0.03 / +0.04 |

(A positive slope raises the elevation over the plane. These come from the evaluation's
standalone model of the plate and pole at the changed elevation, not from a NICE-OPS run.)

### Refraction

`--refraction linear` bends the rays through the profile. For each path it fits a straight
line to the effective sound speed, c(z) + w(z)·h, where w is the wind and h the horizontal
direction of propagation. The fit is least squares over the heights the path spans, widened
to at least 10 m. This is the equivalent linear profile of Harmonoise and Nord2000. In a
linear profile every ray is an arc of a circle, so the eigenrays have closed forms, and the
path supplies:

- **The sphere's direction:** the direct ray's launch angle. At kilometer ranges it is a few
  degrees from the straight line, which is a large fraction of the database's 5° grid.
- **Spreading:** from the ray tube's area, which in a linear profile is exactly
  r / cos(launch angle).
- **Absorption:** along the arc.
- **Ground geometry:** a point receiver's ground model gets the reflected ray's grazing angle
  and length; a plate microphone gets the direct ray's arrival angle. Both ground models are
  otherwise unchanged.
- **Reception time:** from the ray's travel time.
- **Shadow zones:** under upward refraction, receivers beyond the ray that grazes the ground
  get the grazing ray's level, decaying with distance into the shadow. The decay is the least
  attenuated creeping-wave mode. With `--ground_effect` it runs over the ground's own
  impedance; without, over rigid ground, about 24 dB/km at 1 kHz for a gradient of 0.02 /s.
  Real ground is not always the deeper shadow: to first order it attenuates up to 12% less
  than rigid at low frequency and weak gradients, and up to 2.3 times as much once the ground
  counts as soft on the wave's scale. The benchmark's 1 kHz ground, at 0.1 /s, is 146 dB/km
  against rigid's 70. The rate is held to the slope of the exact residue series deep in the
  shadow to 0.3%. It is capped at `--shadow_floor`, where energy scattered in by turbulence
  takes over; `--scattering` puts that energy in directly instead (below). A ground reactive enough to carry a surface wave into the shadow (a porosity
  rate of tens per meter, a thin hard-backed layer, low bands) takes the rigid rate in those
  bands, and the log names them.
- **Near the shadow boundary**, with a point receiver or a flush plate microphone, the level
  comes from the exact field instead. Within a couple of diffraction lengths of the boundary
  (F = 2k₀ℓ², about 100–400 m) ray theory fails on both sides: its ground term, taken at
  the grazing ray's geometry, sinks toward the grazing ground wave, 10–25 dB below the exact
  field over grass, and over rigid ground it doubles where the exact field does not.
  NICE-OPS tabulates the residue series of the path's linear profile per band, relative to
  the ray's own spreading, as a function of distance past the boundary in diffraction
  lengths. That function barely depends on the source's height, which is Fock's result, so
  one table per band and microphone height serves every path. The ray model stands until
  2F before the boundary, blends into the table by 1F before it, and the table (floored 20 dB
  below its boundary level) takes over from there. It follows the exact field to 0.1–1 dB
  for sources from 30 to 300 m. With `--refraction stratified` each band gets its own
  boundary first, as the next section describes.

`atmospheres/example_windy.csv` is an illustrative profile, not a measurement: a neutral
log-law westerly (5.2 m/s at 10 m, 8.5 m/s at 1000 ft) over the ISA lapse. On the demo track
with the B407 database and `--frame_bearing 210`, compared with the same profile under
`--refraction none`:

| | mean | upwind half | downwind half | worst cell | acoustic loop |
| --- | --- | --- | --- | --- | --- |
| no ground effect | -0.77 dB | -1.51 dB | -0.04 dB | -7.85 dB | 2.3 s (1.6 s without) |
| `--ground_effect` | -0.56 dB | -0.80 dB | -0.32 dB | -3.81 dB | 3.6 s (2.1 s without) |
| `--plate_table` | | | | | 6.1 s (5.4 s without) |

Most of the upwind loss comes from the low states near the end of the approach. From 10 m up,
the shadow begins about a kilometer upwind.

What to know before relying on it:

- **Wind directions are compass bearings**, turned into the local frame by `--frame_bearing`,
  whose default is +x east. The shipped AC tracks run +x along 210°, so without
  `--frame_bearing 210` their wind is read 120° off. The log says so whenever wind is present
  and the bearing was not given.
- **One gradient per path.** A log-law profile's steep lowest meters are averaged into it.
  Equivalent-linear models have been found up to 10 dB off a parabolic-equation solution near
  the ground at 500–2000 Hz (van Maercke & Salomons 2009). `--refraction stratified`
  ([below](#refraction-through-the-profiles-layers)) removes this approximation.
- **The effective-sound-speed approximation** overstates the wind by 1/cos(elevation). That is
  worst on steep paths, where refraction matters least.
- **One direct and one reflected ray.** The extra ground-bouncing arrivals of strong downward
  refraction at long range are not included.
- **The source does not feel the wind.** Its advance ratio still comes from the track's ground
  speed.
- **The ray ends at the microphone,** as the straight path does (see [Receivers](#receivers)),
  so uniform air gives the straight result to 1e-9 for every receiver, scattered microphones
  on sloping ground included.
- **The peak search and `--db_cutoff` still work.** They bound a refracted state with the
  same direction tables as a straight one, taken along the ray's launch direction at its ray
  tube's range. A CTest case checks that the windowed run equals a full march.

### Refraction through the profile's layers

`--refraction stratified` drops the single straight-line fit. The effective sound speed is
taken as linear *between the profile's own levels*, and held constant below the lowest and
above the highest. Inside each layer a ray is still a circular arc, so crossing a layer has a
closed form, and a ray is a sum over its layers.

Eigenrays are solved per source–observer pair, over four families of rays:

- launched down, arriving from above;
- launched down, dipping below the microphone and turning back up to it;
- launched up, turning aloft, then arriving from above;
- launched up, turning aloft, then dipping.

The primary arrival is the first family's root, found from the linear model's answer. A
receiver at or above the source is left to the linear model.

The range need not be monotonic along a family, so the pair can have more arrivals than that.
A source at a local maximum of the sound speed is one way. Its shallowest rays skim the slow
layer and travel far, its steepest turn high and travel far, and the arrivals lie between. The
other families are sampled in launch angle, 32 times plus 8 for each bounce the path can hold
(at most 512), and every crossing of the range is solved for.
Each arrival is a level of its own, through its own launch direction on the sphere and its own
ray tube, and they are summed in energy.

**Every arrival is kept.** Stratified propagation retains all solved ray and wave arrivals. It
previously kept only the 24 arrivals with the smallest spreading ranges, which discarded energy
without accounting for source directivity or ground losses. The 10 dB focusing cap applies per
ray; distinct arrivals add in energy and can exceed it in total. Deep ducts can therefore take
longer to evaluate and predict more energy. Launch-angle sampling and the bounce search remain
approximate.

**Surface ducts.** When the air above the source outruns all of the air below it, a ray that
reaches the ground turns aloft and comes back down. A morning inversion does this. A ray that
bounces m times before its last descent travels X_0(p) + m L(p), where L is one hop from the
ground to its turning height and back. Those bounced rays are found the same way, from the same
samples. Each earlier bounce costs the plane-wave reflection loss, and the last is the ground
model's, as for the primary. Beyond the reach of every direct ray, about 2.4 km for a 30 m
source under day 286's morning inversion, the bounced rays are all there is. The strongest
becomes the primary. `atmospheres/example_inversion.csv` is that profile, from run 286276,
without its 1–3 kt of wind: 3.5 K of inversion capped near 350 ft.

Without them, `stratified` was 8–24 dB under the PE past 4,000 ft on day 286's profiles for
100–330 ft sources, and up to 54 dB under it at single ranges. The table compares the energy
mean over 50 Hz – 2 kHz against the PE, per range bin, for a flush receiver:

| day 286, source 100 ft | 1,000–4,000 ft | 4,000–8,000 ft | 8,000–13,000 ft |
| --- | --- | --- | --- |
| run 286264 | +0.4 | −1.1 (was −12.2) | −0.5 (was −24.3) |
| run 286270 | 0.0 | +0.8 (was −7.2) | −2.5 (was −8.1) |
| run 286276 | −0.7 | −1.5 (was −10.7) | −1.1 (was −9.2) |
| run 286282 | +4.4 | +0.2 | −3.1 (was −17.1) |

Downwind on `example_windy.csv` and a log law, at 100 ft, the far third of a 10,000 ft path goes
from 4.0–6.5 dB under the PE to within 1 dB of it. Upwind, the linear lapse and the day-284
profiles do not move.

Two limits remain. Within one pair of samples, two close arrivals at a caustic can still be
missed. Against a 2048-sample reference (12 bearings, sources 5–300 m, 50 m to 10 km), 2.9% of
pairs under `example_inversion.csv` are still more than 1 dB out, and 0.01% under
`example_windy.csv`; a fixed 32 samples left 11.3% and 1.05%. Scaling the count with the hops
cost 4–8% of a stratified run.
Day 286's last profile (run 286288) is an elevated duct. Toward 270 the air is fastest at the
ground, because the log-law headwind assumed below the LIDAR's first gate (46 ft) takes
0.2 m/s off in the lowest foot, and a ray turned back from the inversion near 330 ft turns up
again 1–20 ft above the ground. For a 100 ft source it reaches neither a flush nor a 1.2 m
microphone at 8,000–13,000 ft, which stay on the shadow's floor 11–21 dB under the PE, while
the wave tunnels through the layer (flattening it moves the PE by under 1 dB to 250 Hz and
about 6 dB at 1 kHz). Short of the duct's return lies a skip zone that no ray reaches, which
the PE fills from the caustic 1–2 kft early. Measured day-286 frames agree with the PE: where
`stratified` predicts 10 dB or more under straight rays, it is 10–20 dB too quiet. Those
frames are far from each microphone's peak and do not move the SEL.

It is checked two ways. On a linear profile it reproduces the linear model's closed forms,
the shadow boundary included. On layered profiles its rays match numerically traced ones to a
micrometer. Everything else follows `--refraction linear`: the ray tube, the ground geometry,
reception time, the shadow's decay and floor, and the peak search.

What the layers change is near the ground. A log-law wind is zero at the ground, so for
propagation upwind the fastest air is at the surface. Snell's law then fixes the limiting
ray's launch angle from the speeds at the source and at the ground alone, whatever lies
between. For a source at 100 m in `atmospheres/example_windy.csv`, the layered model's
upwind shadow begins 535 m out, where a fitted line puts it past 1.6 km.

**The wave does not take that boundary at face value.** Rays follow the thin near-ground
shear layer, but the wave resolves it only once its diffraction length is short against the
distance between the two boundaries. At 125 Hz the parabolic equation below shows no shadow
until about a kilometer, and even at 2 kHz it has not reached the layered boundary. So with a
point receiver each band gets an effective boundary between the two,
r_b = r_s + w (r_l − r_s), where w = x / (x + 0.363) and x is the band's diffraction length
at the linear fit's curvature over r_l − r_s. The boundary field is then scaled by the
curvature that makes r_b. Past the layered boundary the lit level is the linear fit's ray:
- The layered model's near-grazing rays there have no reflected partner.
- On a log law its dipping family reaches tens of meters past its own end.

The law's one constant was fitted on twelve log laws (u* 0.2–0.6 m/s, z0 0.01–0.1 m,
sources 30 and 100 m, 125 Hz – 2 kHz) and checked on `example_windy.csv`, which was held out:
6.8 dB rms over the shadow edge against the PE there, against 18.9 dB with the layered
boundary as it stands. On the demo grid with a point receiver, `stratified` then averaged
3.1 dB below `linear` (5th percentile −12.4 dB), where it was 6.1 (−19); with the fixes below,
3.3 (−13.2). Runs without
`--ground_effect` keep the layered boundary and the creeping-wave loss.

**Flush plate microphones** get the same treatment:
- **The field:** the boundary field for a receiver at the plate's top (8 mm for the GR1425).
- **The plate over bare ground, K:** added on top. Near grazing K barely depends on elevation
  (within 0.1 dB from 0.5° to 10°). It is about 0 dB up to 1 kHz, +0.9 dB at 2 kHz, +2.6 at
  4 kHz and +7 at 8 kHz, and it varies with azimuth by ±1 to ±2.5 dB, so it is tabulated per
  band and azimuth at 2° and applied at each pair's own azimuth.
- **At ground level** the boundary law holds up. Against the PE for a receiver 3 mm up, over
  16 log-law cases, the shadow-edge error is 7.1 dB rms, against 10.0 for `linear`. The
  constant that suits both heights together is 0.39, only 0.1 dB better than 0.363.

Two faults in the layered model surfaced with flush microphones, both fixed:
- **Microphones just below the plane.** Ground boards sit a few centimeters below the survey
  datum. The layers hold their lowest speed below z = 0 and the direct-ray solve knew no
  ground, so it sent a ray grazing along under the plane. On the 2017 test's windy day that
  was 33 dB low. Such a microphone is now solved for on the plane.
- **Receivers a few millimeters up.** The ray tube past the boundary came out 0/0. It is now
  spread as a linear profile's limiting ray.

**The measured profiles of the 2017 test's windy day** (LIDAR winds and balloon temperatures,
morning of day 284) exposed three more, all in how `stratified` meets its shadow. They were
found by running the PE on the geometries of the microphones the model got most wrong.
- **An inversion aloft hid the boundary.** Those mornings have an upwind shear near the
  ground under an inversion, so the air a hundred meters up is a little faster than at the
  source. The layered model took any faster air above the source to mean rays from aloft
  reach the receiver, and set no boundary. With none, there was no transition and no
  hand-over to the linear fit: the shadow came in as a 10–35 dB cliff at the layered
  boundary, which the PE puts 30–50% further out. A ray turned down aloft turns up again
  wherever the air is faster still, here in the shear, meters above a board. Rays from aloft
  reach the receiver only if the air aloft outruns everything between the two, which is the
  test now.
- **The layered ray focuses short of its boundary.** The steep lowest layer that sets the
  boundary also bends the rays that just clear it into a caustic, a spike of up to +10 dB the
  PE does not show. From 0.9 of the layered boundary on, the level is the linear fit's ray.
  Against the PE over the last half of the way to the boundary, that is 0.6 dB rms on the
  day-284 profiles (1.4 dB handing over at the boundary itself) and 1.7 on log laws (3.3).
- **Past the fit's own boundary the level stepped.** The model went back to the layered
  shadow ray, whose spreading differs, by 2 dB on a log law and 10 on day 284. The fit's
  limiting ray now carries on, with its spreading continuous.

**The day-286 approaches' hovering end** (sources 15–50 ft, under 20 kt) showed two more.
`stratified` was 6.5 dB over straight rays there on the measured frames, 63–250 Hz, at
microphones 150–450 ft away; the PE on those paths agreed with straight rays.
- **Faster air aloft still hid the boundary near the source.** The test above is right about
  whether rays from aloft reach the receiver at all, but they come down only past their skip
  distance, kilometers out. Short of it the boundary now stands, and is taken back where a ray
  from aloft or along a duct does arrive. Without it the shadow took the layered shadow ray,
  focused 6.5 dB by the log law's lowest meters, plus the whole ground term: 12 dB over free
  field at 63 Hz where the PE gives 4.
- **A low source leaves the transition no room.** When a band's boundary is nearer the
  source than two diffraction lengths, the table has no lit zone to hand over from, and the
  source is under the 5 l the table was checked to. Such a band keeps the ray model to its
  boundary and blends to the table over the next diffraction length.

Together, on the measured frames:

| 63–630 Hz energy-mean error | before | after | straight rays |
| --- | --- | --- | --- |
| sources under 100 ft and 20 kt | +8.9 dB | +4.4 dB | +3.8 dB |
| within 2,000 ft, 2–5° | +4.7 dB | +3.1 dB | +2.4 dB |
| every other height, speed, range or elevation cell | | within 0.2 dB of before | |

The 4 dB left at the hover is in straight rays too, so it is the source. On the day-286 PE
matrix past 1,000 ft, `stratified` against the PE went from 9.5 to 7.5 dB rms at 200–630 Hz
and 11.7 to 9.2 at 800–2,000 Hz (5.4 to 5.1 at 50–160).

**The boundary law, refitted for flush microphones.** With those fixed, the law put flush
microphones' boundaries about twice too far out. Fitted per case against the PE, the best
boundary is about 0.9 of a diffraction length at the fit's curvature past the layered one on
day 284 and 1.2 on flush log laws, where the law gave 1.7–2.7. A microphone 1.2 m up keeps
0.363 (refitted with the day-284 profiles added: 0.387). A flush one takes 0.75, between the
flush log laws' own optimum of 0.47 and day 284's 0.94. Where the linear fit casts no shadow
at all, as over a strong inversion, the law has nothing to interpolate toward. There each
band takes 1.60 (f / 1 kHz)^−0.08 of the layered boundary for a flush microphone, and 1.90
(f / 1 kHz)^−0.08 at 1.2 m, the per-band optima pooled over the calibration set.

Shadow-edge rms against the PE (from its last level above −3 dB to −25 dB), dB:

| | before | now |
| --- | --- | --- |
| day-284 profiles, 18 cases (sources 50, 120, 250 m), flush | 12.9 | 8.1 |
| log laws, 16 cases, 1.2 m | 6.2 | 6.0 |
| log laws, 16 cases, flush | 7.3 | 7.6 |

"Before" for day 284 already has the first fix; without it the cliff made those shadows
worse still. The best any boundary per case can do is 2.7 dB pooled, and no smooth law found
gets below about 6: the optimum is sharp where the shadow wall is steep, at high frequency.

Against the measured day-284 SEL, held-out spheres, `stratified` on the plate (729 mic-runs,
SNR ≥ 10 dB, records that miss the flyover left out), median / rms in dB:

| | before | now |
| --- | --- | --- |
| upwind, peak elevation 5–10° (65) | +0.45 / 3.72 | +0.28 / 3.06 |
| upwind, 10–20° (80) | +0.99 / 2.24 | +0.71 / 1.87 |
| all | +0.09 / 1.85 | +0.09 / 1.70 |

The calm day 286 barely moves (overall 1.46 → 1.45 dB rms): its winds cast few shadows.

**A calm floor at the ground, and flat air at the fastest height.** A log law written with no
wind below its roughness length leaves the lowest few centimeters calm, and a profile written to
three decimals leaves them isothermal too: a floor exactly as fast as the ground. On the 2017
test's windy day, a station-wind log law with z0 = 0.1 ft (levels at 0, 0.01, 0.03, 0.1, 0.3 and
1 ft, the first four calm) predicted 22–28 dB low at 2–5 low-elevation upwind flush microphones a
run: run 284202's microphone 3, at a 12° peak elevation, measured 78.3 dB and was predicted near
50, against 78.6 without the floor. It took day 284's SEL rms from 1.63 to 4.18 dB. A floor only
0.1 ft thick did it (levels 0, 0.1, 1 ft: −28.2 dB), where a level at 0.3 ft carrying some wind
moved the footprint by −0.15 dB. `--tunneling` and `--skip_fill` did not help.

*Cause.* A ray whose Snell invariant is 1/c, c the fastest air between microphone and source, is
horizontal where the air is that fast. With any gradient there it turns, and the range of the
rays approaching it converges. In flat air it runs along it without limit. The downward family's
last ray, 1e-12 short of that invariant, crossed a 3 cm floor at 1.4 microradians, and the
shadow boundary went 20 km out — a figure set by the inset alone. With the boundary gone there
was no hand-over to the linear fit and no shadow transition, and the pairs short of 20 km took
near-grazing rays whose tubes were tens of times their range. A microphone 1.2 m up was never
affected: its limiting ray dips and turns at the floor's top, in the gradient above it.

*Fix.* No wave can be steered along a stretch thinner than its own vertical wavelength at that
angle: one confined to a height h diffracts out of it within about k h². So the rays see the
air at the fastest height only as finely as the highest band the run evaluates resolves it. The
limiting ray's turning region is taken as the heights within one radian of vertical phase of the
fastest air, k ∫ sin θ dz = 1, and its horizontal travel across that region, of height h, is
capped at 1.5 k h². In a linear gradient ray theory travels (4/3) k h² there (√(2Rh) against an
Airy height h with k²h³ = 9R/8), 8/9 of the cap, so the cap never binds on a resolved turning,
and such profiles give the same rays bit for bit. Flat or nearly flat air gathers no phase and
is capped. Where it binds, the downward rays stop at the invariant whose travel across the
region is the cap, whose ray still crosses the region at a small but finite angle. The boundary
is where that ray lands, and past it the shadow, the hand-over and the transition run as they do
without the floor.

The cap's margin over (4/3) k h² is what a floor can still add, so it is kept small. With a cap
of 2 k h² the 0.1 ft floor moved the demo grid's flush microphones by up to 0.26 dB, and by
0.67 dB when the run's top band was 200 Hz (`--band_limits 20 224`), a floor a fiftieth of that
band's wavelength; at 1.5 k h² it is 0.12 and 0.19 dB. The margin scales as k^(−1/3) (h goes as
k^(−2/3)), so it is a run with a low top band that feels it most.

The region belongs to the air, not to the path. Where the path ends before the radian is
gathered — a source within a radian of the fastest height — the phase is integrated on past
the source through the profile, and past the microphone down to the ground, for as far as the
limiting ray could go, and only the ray's travel inside the path is held to the cap. Cut off at
the path's ends instead, a smooth gradient bound: under a plain lapse a 1.5 m source over a
1.2 m microphone had its boundary brought in from 789 m to 34 m at 10 kHz, and at a 200 Hz top
band sources below about 10 m had theirs moved in any light wind.

Nor is the cap applied where the fastest air is the source's own and there is faster air above
it, as under any downward refraction. There the downward family's last ray is the horizontal
launch, where the rays launched up to turn aloft take over: it is no end of the rays, and
capping it opened a gap between the two families inside the lit zone that no ray reached.

The highest band is used because the ray geometry is one for all bands, and the highest band's
k is the largest, so the cap removes only what no band of the run resolves. The shadow
transition then moves each band's boundary toward the linear fit's, as it does for any profile.
With `--band_limits` the top band changes, and a floored profile's boundary can move slightly
with it, by the margin above. The region is a set of heights, found by integrating the phase, so
the same air written with more levels gives the same boundary.

For a 30 m source upwind on a 7 kt log law at 10 kHz, the boundary moves as follows:

| calm floor | none | 1 mm | 0.1 ft | 1 ft | 1 m |
| --- | --- | --- | --- | --- | --- |
| plain ray theory | 186 m | 893 m | 21.8 km | 216 km | 707 km |
| now | 186 m | 186 m | 187 m | 233 m | 544 m |

A thick floor is followed as far as the highest band resolves it: a meter of calm air is ten
wavelengths thick at 10 kHz and a third of one at 100 Hz, where the same source's boundary is
215 m (198 m under a 1 mm floor, which only the margin moves). The limit cannot reach a profile
whose lowest level stands above the ground. The speed is held constant below the lowest level,
so a profile starting at a 46 ft LIDAR gate has a 14 m floor, which the cap follows for
kilometers. Give a profile a level at the ground.

Checks:
- **Day 284.** The z0 = 0.1 ft log law's SEL rms goes from 4.18 back to 1.63 dB (729 mic-runs,
  held-out spheres), the same as z0 = 0.003 ft, which has no floor. A z0 = 1 ft log law, calm
  to 1 ft, goes from 3.95 to 1.90.
- **Run 284202.** With the floor, the worst microphone is now 0.07 dB from its floorless level
  (it was −28.16), with `--tunneling` and `--skip_fill` too. Every set of levels below 1 ft that
  was tried, the full z0 = 0.1 ft set included, lands within 0.15 dB of the floorless run, as a
  single level at 0.3 ft carrying wind always did.
- **Unchanged runs.** The paper's measured-profile `stratified` runs move by at most 0.0006 dB
  on day 284 and 0.0002 dB on day 286, and their scores are unchanged. Over a sweep of log laws
  without a floor (headwinds of −1 to 3.6 m/s against the lapse, sources from 1.5 to 100 m, top
  bands from 200 Hz to 10 kHz) the rays are those of plain ray theory but for a few interior
  speed maxima under a light tailwind, which move by under 1e-7 dB.
- **CTests.** `stratified_calm_floor_plate_is_floorless` runs `example_windy.csv` with a 0.1 ft
  calm floor (`atmospheres/example_windy_calm_floor.csv`) against the profile without it on the
  demo grid. The flush microphones differ by 0.12 dB at most, against up to 39 dB before. The
  1.2 m microphones differ by 0.05 dB. Both are held to 0.3 dB.
  `stratified_calm_floor_plate_200hz_golden` pins the same floored run with the top band at
  200 Hz, which is what says the run's own top band reaches the solve.
- **The goldens did not move.** `example_inversion.csv` is exactly isothermal and calm from 180
  to 190 ft, but a source there has faster air above it, where the cap does not apply. Five cells
  of `stratified_windy_point` move by up to 4e-6 dB: near crosswind the lowest layer of
  `example_windy.csv` is within 1e-4 m/s of flat, and the rays along it are capped as they
  should be.
- **Cost.** None measurable on `example_windy.csv`, about 3% on `example_inversion.csv`.

Costs, with a point receiver: the demo grid's acoustic loop took 12.0 s when `stratified` went
in (3.3 s linear, 2.1 s without refraction). A 96-microphone scattered array like
`arrays/coyle_approach.csv` takes 35 ms. Searching for every arrival costs about half as much
again. A whole demo-grid run (B407, `--ground_effect`, `--frame_bearing 210`, 14 threads) goes
from 23.6 to 37 s on `example_windy.csv` and from 57.5 to 86 s on `example_inversion.csv`. The
search is skipped past the hand-over to the linear fit, unless a duct gives bounced rays to
look for. Focusing at a caustic is capped at +10 dB.

### Low sources near grazing: the resolved shadow (`--shadow_model resolved`)

`--shadow_model resolved` is an opt-in change to how `--refraction stratified` places and fills a
ground-grazing shadow. The default, `legacy`, is the transition described above, unchanged bit for
bit (the goldens and every existing CTest; the paper's flight runs reproduce to 0 dB). It was
written for low sources under light wind, the R44's landings at Eglin, where the layered ray
model's shadow was 5-11 dB deeper than a parabolic equation's at 1000-1500 ft. Three things
caused that (2026-10-05 code audit):

- **M2. A layer the wave does not resolve sets the boundary.** Ray theory fixes a ground-grazing
  boundary by c(z_s)/c(0) alone. On run 229355 the log law extrapolated below the lidar's 46 ft
  gate takes 0.25 m/s off in the lowest foot. That sets the layered boundary r_s at 446 ft,
  where the fit's r_l is 1306 ft and the PE's onset follows r_l. Below 4 kHz the layer gathers
  less than a radian of vertical phase.
- **M1. The lit side's ground term collapses.** The equivalent straight geometry
  h' = L sin(grazing) - h_r goes to 0 approaching the boundary, so the ground model sits in its
  grazing null, 8-20 dB under the true field. A band whose boundary is within two diffraction
  lengths of the source (cramped) keeps that ray model up to the boundary.
- **M3. The decay is tied to the boundary.** a = c0 root^2 / r_b^2. With the fit casting no
  shadow, the fallback r_b = 1.6 r_s gave 0.07-0.11 /s (day 284, run 200, 200 ft source:
  -27 dB at 2500 ft and 1 kHz, where the PE gives -5.6).

What the option does:

1. **Each band sees only the profile it resolves.** For band f it finds h_f, where
   k_f int_0^h sqrt(2 max(c(0) - c(z), 0) / c(0)) dz = `--shadow_phase` (1 rad, fixed before
   any scoring; it is the radian the flat-floor limit `resolved_pass` uses). Below h_f the
   effective sound speed is replaced by its tangent there, the slope of the layer containing
   h_f. That band's layered boundary r_s(f), and its linear-fit boundary r_l(f), come from that
   profile. r_l(f) is the path's fit corrected by the least-squares line through the tangent's
   change. Both go through the same `band_boundary` law as before. On a linear profile
   h_f = (3 phase / (2 k sqrt(2|a|/c0)))^(2/3), which is 1.31 creeping-wave layers at one radian.
   There the resolved profile is the profile itself, so nothing moves.
2. **The decay keeps the boundary's curvature.** The task asked for the decay to be taken from
   the resolved profile's gradient over the creeping-wave layer (the tangent below h_f). That
   was built and scored, and is kept in `niceops_pe_compare --tangent_decay`. It is worse on
   every upwind class of the synthetic set (log-law upwind 3.93 -> 5.35 dB rms against the PE,
   log law plus shear aloft 4.50 -> 5.64, with the ray's cramped ground term). The decay rate,
   s and the diffraction length therefore stay c0 root^2 / r_b(f)^2, each band at its own
   boundary.
3. **A cramped band takes the exact field.** For k^2 linear in height the residue series is the
   field. Where a band is cramped (r_b/F < 2, which means sqrt(y_s) + sqrt(y_r) < 2 in Airy
   lengths), it takes the series at the source's own height, with the curvature whose boundary
   is the band's. It is tabulated per band over (xi, s, sqrt(y_s)), with sqrt(y_s) in 0-2 every
   0.4 and xi every 0.2 (`ShadowTransition::cramped_level`). The modes and H0(kappa r) are
   shared across heights (`upward_pressure_heights`). Near the source, where the series is not
   trusted (NaN), the band falls back to the ray model. The choice among `series`, `straight`
   (the ground model along the straight path through the lit side and the blend) and `ray`
   (as legacy) was made on the synthetic set alone (below). `straight` and `ray` remain as
   `--cramped_ground` values.

The per-band boundaries depend only on bearing, source height and receiver. A run tabulates
them once per trajectory (`ResolvedShadowTable`: 2 degrees of bearing by nodes even in
sqrt(z_s), interpolated in 1/r_b). Every band whose h_f lies in one layer shares one profile,
so it costs one boundary solve per layer touched; the table builds in about 0.01 s on
`example_windy.csv`. The cramped series table costs about 2-3 s of CPU per band per receiver
height, so 24 bands take 5-8 s on 14 threads and about a minute at `-j 1`. A pair pays one
table lookup per band. The bound stays exact: the series and V are held to the transition's
largest power ratio, and the straight ground term is within the ground's own bound.
`stratified_resolved_{point,plate,straight_point,duct_point}_windowed_is_exact` check it, and
the `[resolved][bounds]` Catch2 case sweeps it. One approximation: the table indexes receivers
by their height above the ground (the run's point heights and the plate), so scattered
observers whose surveyed z is not 0 are looked up at that height.

**Validation against the PE.** `niceops_pe_compare --shadow_model resolved` writes a
`resolved_<f>` column. The score is code-audit's: the rms of (model - PE) over ranges of
300 ft and beyond where the PE is within 20 dB of straight rays, per case and band, pooled.
Bands 125-4000 Hz, ranges to 4000 ft, variable porosity 200, flush (3 mm) and 1.2 m point
receivers. The synthetic set was generated and frozen (sha256 manifest) before any variant was
scored: linear lapses of -0.01 to -0.06 /s; neutral log laws through a 46 ft value of 2-12 kt,
extrapolated as the harness does; the same plus +6 kt of shear from 46 to 200 ft; upwind,
downwind and crosswind; sources of 5-200 ft. 308 cases. The measured set is held out:
R44 229354/355/363/369/372 toward 140 (sources 5-80 ft); all 16 day-284 profiles toward their
46 ft upwind bearing (20-200 ft); and the calibration cases (`example_windy.csv` and three
synthetic log laws, 30 and 100 m). 194 cases.

| class (case-bands) | legacy | resolved (series) | resolved (straight) |
| --- | --- | --- | --- |
| linear lapse (336) | 3.38 / -0.06 | 1.39 / +0.30 | 2.97 / +1.23 |
| log law, upwind (322) | 4.77 / -1.75 | 3.70 / +0.74 | 4.13 / +1.03 |
| log law + shear, upwind (167) | 6.73 / -2.71 | 4.25 / -0.20 | 4.82 / +0.16 |
| log law, crosswind (= adiabatic lapse) (336) | 4.01 / -0.89 | 0.73 / +0.04 | 2.66 / +1.18 |
| downwind (504) | unchanged | unchanged | unchanged |
| **synthetic, pooled (1833)** | **5.30 / -0.74** | **4.24 / +0.25** | 4.70 / +0.82 |
| R44 day 229 (300) | 6.58 / -2.81 | 5.19 / -0.47 | 5.70 / +0.28 |
| B407 day 284 (748) | 5.67 / -2.42 | 5.16 / +0.03 | 5.18 / +0.10 |
| calibration cases (96) | 4.28 / -1.93 | 3.02 / -0.18 | 3.02 / -0.17 |
| **measured, pooled (1144)** | **5.82 / -2.48** | **5.03 / -0.12** | 5.18 / +0.12 |

(rms / mean of model - PE, dB.) Resampling whole profiles, the change in pooled rms is -1.06 dB
[-1.71, -0.61] on the synthetic set and -0.79 [-1.16, -0.41] on the measured one. On R44
229355, 20 ft source, flush receiver, at 1000/1250/1500 ft and 1 kHz the PE gives
-13.6/-17.3/-20.1 dB; legacy gives -18.6/-28.2/-28.9 and resolved -14.2/-18.0/-21.7. On a
linear -0.02 /s lapse, 20 ft source, flush, resolved is within 0.3 dB of the PE from 800 to
2000 ft at 500 Hz and 1 kHz. Legacy is up to 17 dB low there. That case is pinned in
`tests/test_resolved_shadow.cpp`.

**Where it is worse.** 58 of the 194 measured cases and 10 of the 308 synthetic ones worsen by
more than 0.5 dB rms. Almost all are strong winds with a 1.2 m receiver: day 284's 20 ft
sources go 1.6-2.1 -> 3.7-5.3 dB, and the 12 kt log laws worsen by up to 2.5 dB. On strong
log laws the PE resolves the lowest layers better than one radian implies. The resolved
boundary, which then still goes through `band_boundary`'s interpolation toward r_l, lands too
far out, by 1-2 diffraction lengths. Above 4000 ft, on two day-284 profiles with 300 and 500 ft
sources (8000 ft, flush), the effect splits by profile. Where the fit casts no shadow (run 200)
the PE stays near straight rays and resolved is much closer (rms 19.5 -> 13.8 and
14.0 -> 7.9). In a true upwind shadow (run 205) resolved is 1.6 dB loud on average
(rms 2.0 -> 3.9 and 2.7 -> 4.1). Deep in shadows (PE more than 20 dB under straight) neither
model is right: the floor holds both 25-35 dB above a still-air PE.
`--shadow_phase` 0.5 (sensitivity): measured 5.05 dB, day 284 5.13. `--shadow_phase` 2: 5.28.

**Flight data** (harness `fidelity_study` gates and common samples, plate receivers, profile
atmospheres; paired with `stratified_plate`; 95% intervals from resampling runs):

| set | mic-runs | stratified_plate rms | resolved rms | Δ rms [95% CI] |
| --- | --- | --- | --- | --- |
| R44 approaches, all | 585 | 1.84 | 1.53 | -0.31 [-0.42, -0.19] |
| R44 approaches, 0-5° | 45 | 4.74 (mean -3.89) | 3.55 (mean -2.88) | -1.19 [-1.55, -0.73] |
| R44 approaches, 5-10° | 118 | 2.09 | 1.70 | -0.38 [-0.63, -0.12] |
| B407 day 284, all | 729 | 1.69 | 1.78 | +0.09 [+0.02, +0.16] |
| B407 day 284, < 10° | 73 | 2.98 (mean +1.00) | 3.42 (mean +2.32) | +0.44 [+0.03, +0.75] |
| B407 day 286, all | 1111 | 1.35 | 1.36 | +0.00 [-0.00, +0.01] |
| B407 day 286, < 10° | 62 | 1.91 | 1.94 | +0.03 [-0.06, +0.11] |
| R66 approaches, all | 874 | 1.37 | 1.41 | +0.04 [+0.02, +0.05] |
| R66 approaches, < 10° | 99 | 1.99 | 2.13 | +0.15 [+0.10, +0.20] |
| R44 turns (held out), all | 4992 | 1.55 | 1.56 | +0.00 [+0.00, +0.00] |
| R44 turns, < 10° | 11 | 0.94 | 0.90 | -0.04 [-0.07, +0.00] |

On the R44 this matches the per-run PE tables, 1.47 overall and 3.51 at 0-5°, with no table. It
fails the windy day's do-no-harm gate (day 284 < 10°: rms <= 3.15, mean +1.00 ± 0.3, Δ upper
bound <= +0.15; all <= 1.75). That is the strong-wind over-reach seen against the PE, and it is
why the option stays opt-in. Removing per-microphone offsets, the R44 rms goes 1.13 -> 1.02
(0-5° 1.82 -> 1.33). The approaches' mic 1/2/4 group offset is not propagation and does not move.
`--cramped_ground straight` scores the same on day 284 (3.44 at < 10°). `--shadow_phase` 0.5, as a
sensitivity only, gives 3.24 there (mean +2.00) and 1.58 on the R44 approaches (3.79 at 0-5°). The boundary placement,
not the cramped field, is what moves it.

### Micro-shadows and skip zones

Under day 286's last profile (run 286288) the air toward 270 is fastest at the ground, because
the log-law headwind assumed below the lidar's 46 ft gate takes 0.2 m/s off in the lowest
foot, and the inversion near 330 ft makes an elevated duct above it. A 100 ft source's sound
turned back by the inversion turns up again 1-20 ft above the ground, so ray theory delivers it
to neither a flush nor a 1.2 m microphone; and short of where it first comes down, about
12,000 ft out, no ray arrives at all. The parabolic equation hardly notices the near-ground
layer (flattening it moves the PE by under 1 dB to 250 Hz) and fills the skip zone from the
caustic. Three treatments were built and compared (harness `propagation/micro_shadow*.py`):

- `--tunneling` adds the rays that turn within reach of the microphone, each with the Airy
  tail below its turning point, Ai^2((3/2 omega I)^(2/3)), I the integral of
  sqrt(p^2 - 1/c^2) from the microphone to the turning height. Near the family's caustic the
  pair of rays is replaced by the caustic's band-averaged Airy field.
- `--skip_fill` adds, for each family from aloft whose range has an interior extremum with the
  receiver on its dark side, the caustic's Airy tail, I = C Ai^2(omega^(2/3) (X_c - r)
  (2 / X'')^(1/3)), C matched to the two rays' energy on the lit side.
- `--propagation_table` built from the PE for each run's profile, beside `--refraction
  stratified` for the pairs it does not cover.

| 100 ft source, flush, 2,000-13,000 ft; `stratified` - PE, mean / rms dB | 50-160 Hz | 200-630 Hz | 800-2,000 Hz |
| --- | --- | --- | --- |
| Day-286 PE matrix (five profiles), neither | -5.4 / 9.8 | -6.1 / 11.1 | -3.6 / 9.9 |
| `--tunneling` | -3.4 / 7.3 | -3.8 / 8.5 | -1.4 / 7.3 |
| `--skip_fill` | -4.5 / 8.9 | -5.6 / 10.5 | -2.7 / 9.2 |
| both | -1.7 / 5.1 | -2.3 / 6.5 | +0.3 / 6.7 |

Against measured day-286 frames where `stratified` predicts at least 10 dB under straight rays
(1.4% of frame-bands to 2 kHz; sources under 300 ft, 0.4-2 degrees, 4-9 kft), energy-mean
predicted - measured:

| | straight | `stratified` | `--tunneling` | `--skip_fill` | both | PE table (5 runs) |
| --- | --- | --- | --- | --- | --- | --- |
| 25 runs, 50-160 Hz | +3.8 | -10.5 | -7.3 | -3.6 | -2.1 | |
| 25 runs, 200-630 Hz | -0.1 | -15.5 | -10.7 | -7.2 | -5.0 | |
| 25 runs, 800-2,000 Hz | -7.6 | -20.6 | -20.6 | -20.7 | -18.9 | |
| 5 runs, sources under 150 ft at 3-6 kft, 50-630 Hz | +1.0 to +1.6 | -11.5 to -13.3 | | | -2.5 to -4.5 | +0.3 to +2.5 |

None moves the SEL: the exposure within 10 dB of each microphone's measured peak is unchanged
to 0.02 dB under every option. All three are off by default. What they leave:

- The tails above 800 Hz. The tunneled field at the ground and the skip zone's edge short of
  the direct rays' reach (a duct's direct family ends without a caustic, and the PE diffracts
  past it) are still 10-20 dB low there.
- On run 286276's far frames (150-330 ft sources, 6-20 kft) every model but plain `stratified`
  is too loud, the PE table by 13.5 dB and straight rays by 8.8: something the PE leaves out.
- The table is read across 45 degree bearings and eight source heights; off those nodes it is
  within about 1 dB of a direct PE run at 63-250 Hz but 5 dB rms overall in this duct, almost
  all at 1-2 kHz, where the return's caustic moves between nodes. Each run's table took
  2.5 minutes to build.
- Cost on a day-286 run: `--tunneling` about +20%, `--skip_fill` 2-4x, a table about the
  ray model's own.

### Propagation loss tables

`--propagation_table` takes what the path does to each band from a table computed elsewhere,
so that a full-wave method -- the parabolic equation of [Reference solutions](#reference-solutions), or something else later -- can
carry a footprint. The table holds each band's level relative to a point source's free field
at the straight-line distance: the ground, refraction, shadows and diffraction, everything
beyond spherical spreading. The run then takes, per band,

    L = L_sphere(straight-line direction) - spreading(straight-line range) - absorption + table

which is the composition `--ground_effect` uses, with the table where the ground model and
the ray model were. Atmospheric absorption stays the run's (homogeneous, or layered through
`--atmosphere`) unless the table says it includes it.

It assumes one stratified atmosphere over flat, uniform ground, for now. The level then
depends on where the observer is only through the compass bearing and horizontal range from
the source and the two heights above the ground, which are the table's axes:

| Axis | |
| --- | --- |
| `bearing` | compass bearing the sound travels toward, degrees, uniform from 0 over 360; one value means the same in every direction |
| `source_height`, `receiver_height` | above the ground, m, ascending |
| `range` | horizontal, m, ascending |
| `band` | band centers, Hz; the run's bands within `--band_limits` must all be there |

Values are interpolated linearly in dB along each axis. A NaN means the method has no value
there. A pair the table does not cover, by a NaN among the values it would interpolate or by
falling outside an axis, takes the run's own propagation instead: `--refraction` if it is
on, `--ground_effect` if that is, free field if not. A flush microphone on a plate
(`--plate_table`) reads the table at its lowest receiver height and adds K, the plate's level
over bare ground there, per band and azimuth at 2 degrees, as the shadow transition places a
plate. With `--refraction` too, the run evaluates every state rather than windowing from the
peak, since the bounds were not built for the two geometries together. The log says how many did. Receivers are the exception: a microphone
height outside the table's is refused, since that table was built for another run. Bearings
are compass bearings, so a run in a rotated frame needs `--frame_bearing`, as with wind.

The file is netCDF, SI whatever `--units` says: a `relative_level(bearing, source_height,
receiver_height, range, band)` array, one coordinate variable per axis, and global attributes
`propagation_table_version` (1), `method`, `includes_absorption` and `provenance_*`.
Anything that can fill that array can be read; `src/propagation_table.hpp` is the
specification.

`niceops_propagation_table` writes one, by either of two methods:

```sh
# The parabolic equation through a profile, 36 bearings where it has wind
niceops_propagation_table --method pe --atmosphere atmospheres/example_windy.csv \
    --database databases/B407_fixed_load_factor.nod --band_limits 20 2000 \
    --source_heights 40 90 140 190 240 --receiver_heights 3.937 \
    --max_range 15000 --step 10 -o tables/windy.nc

# The straight-ray ground model --ground_effect uses: the same footprint, through a table
niceops_propagation_table --method ground_model --database databases/B407_fixed_load_factor.nod \
    --sound_speed 335.636 --flow_resistance 225 --source_heights 40 50 60 ... 1130 \
    --receiver_heights 3.937 --max_range 15000 --step 10 -o tables/ground_model.nc
```

The second is the check on the plumbing. On the demo case its footprint matches
`--ground_effect --exact_ground` to -0.002 dB mean and 0.19 dB at worst with 10 ft spacing,
and to 0.02 dB at worst with 5 ft spacing. The ground effect itself moves cells by up to
4.3 dB. CTest runs that comparison (`propagation_table_is_the_ground_model`), and checks that
the peak search stays exact with a table that adds gain (`propagation_table_windowed_is_exact`):
the search's bound takes each band's largest table value.

With `--method pe` through `atmospheres/isa.csv`, 21 bands from 20 Hz to 2 kHz, 23 source
heights every 50 ft and ranges to 15,000 ft every 10 ft, the table takes 103 s on 14 cores. The
footprint's acoustic loop goes from 0.06 to 0.08 s. The table covered 99.7% of evaluations; the
rest were steeper than the PE's 50 degrees and took the ground model. Against the ray model's
footprint it reads -0.32 dB mean (-0.9 to 0 dB, 5th to 95th percentile; -2.8 dB at worst), the
standard lapse rate's mild upward refraction thinning the field at long range.

**Threads.** `-j`/`--threads N` sets how many worker threads the table's marches (and the
ground model's pairs, and `--validation_csv`'s reference solves) run on. The default, 0, takes
the environment variable `NICEOPS_THREADS` when it is set and otherwise one thread per hardware
thread, as before; `--threads` wins over the variable. A `NICEOPS_THREADS` that is not a
positive whole number is refused rather than ignored. On a shared machine, `export
NICEOPS_THREADS=6` keeps every table build to six cores without editing scripts. The table is
the same bits at any count, turbulent tables included (each march's seeds are fixed before the
work is dealt out); CTest checks that (`pe_tools_thread_count_is_invisible`).
`niceops_pe_compare` takes the same option and variable for its bands.

Interpolating between source heights in dB is where the table's error is. Over the ground model it
shrinks fourfold as the spacing halves, but the PE's band-averaged field is only as smooth as
its band. Build tables a few wavelengths apart at the highest band, or densely at the lowest
heights, where the interference changes fastest.

**Checking a table between its nodes.** With `--validation_csv FILE`,
`niceops_propagation_table` solves its method again at the middle of every cell (source
height, receiver height and range, and the middle of each bearing sector) and compares the
table's interpolation there with it. The file has one row per point and band: the table's
level, the reference, the error and whether the table covered the point. Each band's bias, rms
and largest error are printed. The run exits with status 2, keeping the table, if a covered
point is off by more than `--validation_max_error` dB (default 1), if a band has no covered
point, or if a gap in the table hides a point the reference reaches. For `--method pe` the
reference is a PE on twice the grid density, with `--validation_samples` frequencies per band
(default 21). Midpoints do not bound every interpolation error, and a PE reference is not a
measurement; [FIDELITY.md](FIDELITY.md) has the error budgets measured so far.

**Reciprocity: one march per receiver height.** A PE march puts one source in the field and
reads every height from it, so the table can be filled from either end. `--march source`
marches once per bearing and source height and reads the receivers, as the tool always did.
`--march receiver`, the default, marches once per bearing and receiver height with the source
put at the receiver's height, and reads every source height. In the effective-sound-speed
model the field from a source at z_s to a receiver at z_r along a bearing equals the field from
z_r to z_s along the reverse bearing with the wind reversed, and that path's effective profile
c + w·h is the same, so it is the same march with the heights exchanged. A table has a couple
of receiver heights and many source heights, so it is several times fewer marches.

On the grid this holds only if a source near the ground is placed as the adjoint of the
receivers' read, and it was not. The fourth-order ground row is a half cell whose weights differ
from Numerov's rows above it, so a source within three grid points of the ground came out
wrong: one at the first node 1.34 dB high and one at the second 0.37 dB low against the exact
two-ray field, and a source on the ground (a flush microphone, in a reciprocal march) had no
self-starter at all. The source vector now divides each of those nodes by its exact weight, a
function of the mode's vertical wavenumber and of the ground's admittance applied as an
operator (`place_ground_source` in `src/pe.cpp` has the derivation). A source at 0 to 3.5
grid spacings then agrees with the receiver read at the same height to 0.003 dB at grazing and
within 0.1 dB at 45°, and with the exact field over rigid ground and grass alike
(`tests/test_pe.cpp`). Sources from four spacings up are placed exactly as before. In tables
built before this, the 5 ft row's 50-100 Hz bands were up to 0.39 dB off past the first range
(0.12 dB mean); nothing else moved by more than 0.0002 dB.

Rebuilt both ways from the frozen recipe of the 2017 table study (9 source heights 5-700 ft,
flush and 1.2 m receivers, 8 bearings, 50 Hz-4 kHz, ranges to 4000 ft every 50 ft), on 14 cores:

| run | `--march source` | `--march receiver` | receiver - source, by elevation (cells above -40 dB, max \|Δ\|) |
| --- | --- | --- | --- |
| 229355 (R44, stable morning) | 220-234 s | 91-106 s | 0.013 dB below 5°, 0.001 at 5-10°, 0.003 at 10-20°, 0.025 at 20-35°, 0.11 at 35-50° |
| 229369 (R44, unstable) | 231 s | 91 s | 0.008, 0.001, 0.003, 0.023, 0.10 |
| 284200 (B407, day 284) | 249 s | 96 s | 0.013, 0.003, 0.005, 0.026, 0.11 |

rms over all covered cells is 0.003-0.010 dB. The 35-50° differences are the 1.2 m microphone
at 63-80 Hz, where it sits on the second grid node, and the flush one's steep-angle residual
(0.035 dB); the source-marched value is the closer to the exact field there. Two other places
differ more and neither matters for a footprint. The first range (one step from the source, a
single self-starter step that still carries modes past `--max_angle`) differs by up to 0.29 dB
at 50 Hz for the 5 ft source; against the exact field there the receiver march is the closer
(6.1 dB, source march 6.5, exact 5.9). And still-air deep shadows below -40 dB re free field,
2-10% of the cells (2-4 kHz past 1 km), differ by up to 0.55 dB: the absorbing layer's floor
follows each march's domain, which reaches 250 m past the highest height the march reads, and
that domain is also why the speed-up is 2.2-2.6x rather than the 4.5x the march count promises.
Per-band numbers are in the study's `compare_v2.txt`.

**Turbulence in a table** (report only; nothing defaults to it). `--gamma_t`, `--outer_scale`,
`--seed` and `--realizations` march each frequency sample through phase screens, as the PE does
in [Turbulence, and the energy it scatters into shadows](#turbulence-and-the-energy-it-scatters-into-shadows),
`--realizations` times with distinct seeds, and each cell is the energy mean over the band's
samples and the realizations; every march has its own seeds. `--gamma_profile FILE` (CSV of z
in `--units` and γ in m^-2/3, linear between, held outside) gives a strength that varies with
height: the screens are drawn at `--gamma_t` (default the profile's largest value) and each row's
phase scaled by sqrt(γ(z)/γT), so each height has its own variance. That is an approximation to
a field whose statistics change with height: the vertical spectrum and the outer scale stay a
homogeneous field's, where near the ground the eddies shrink with height. A turbulent table is
reciprocal only in its ensemble mean, not realization by realization, so the two marches agree
to sampling noise. `--gamma_t 0` gives the still table bit for bit whatever the other options
say. `tests/test_pe_table.cpp` checks the screens' variance per height, that a profile scales
the march's screens by its square root (mean-field decay), that three times the strength puts
4.77 dB more into a deep shadow, and that the realizations average as stated.

```sh
# The frozen recipe of the 2017 table study (the harness's propagation/pe_table_study.py)
niceops_propagation_table --method pe --atmosphere RUN_atmosphere.csv \
    --database databases/RO-44_fixed_load_factor.nod \
    --source_heights 5 10 20 35 50 100 200 350 700 --receiver_heights 0 3.937 \
    --max_range 4000 --step 50 --bearings 8 --band_limits 50 4000 \
    --impedance_model variable_porosity --flow_resistance 200 -o tables/RUN.nc
# The same, turbulent: 3e-6 (HARMONOISE "moderate"), 4 realizations a frequency
    ... --gamma_t 3e-6 --realizations 4 -o tables/RUN_turbulent.nc
```

## Reference solutions

The ray models are checked against a full-wave solution, a wide-angle parabolic equation, and the parabolic equation
is itself checked against exact solutions and the standard benchmark cases for outdoor propagation.

### Checking the refraction models against a parabolic equation

`niceops_pe_compare` runs one path through all four models and writes each band's level
relative to free field: straight rays, `linear`, `stratified`, and a wide-angle parabolic
equation. The last is a full-wave solution, with the ground's impedance written into the
grid. Absorption and turbulence are left out of all four, unless `--gamma_t` puts turbulence
in the PE (below).  With `--shadow_model resolved` it writes a fifth column per band,
`resolved_<f>`: `stratified`'s rays with the transition `niceops --shadow_model resolved` takes,
with `--shadow_phase` and `--cramped_ground` as there, and `--tangent_decay`,
`--fit_linear_boundary` and `--no_boundary_law` as sensitivities (see "Low sources near grazing").
Its per-band h_f, r_s(f), r_l(f) and r_b go to stderr. `-j`/`--threads` (default: the
environment's `NICEOPS_THREADS`, else one per hardware thread) sets the threads the PE's bands
run on, as for `niceops_propagation_table`; the output does not depend on it.

The PE is marched in range by a split-step Padé method (Collins, JASA 93, 1736 (1993)): one
range step applies the one-way equation's exact propagator, replaced by a rational function
that is accurate to 1e-7 for every wave within 50° of horizontal (`--max_angle`, up to 60) and
costs one tridiagonal solve per term. A step is eight wavelengths long (`--range_step`) and
takes about ten terms. The result does not depend on the step: between 2 and 32 wavelengths it moves
by less than 0.001 dB, 60 dB into a shadow. Waves steeper than `--max_angle` are absorbed
rather than carried, so **a level on a path steeper than that from the source is written as
NaN** (within 83 m of the point under a 100 m source, at the default).

Two more things make the field, and not just the propagator, good to that angle. The march
starts from Collins' self-starter (`--starter self`): the first range step with the point
source's modal amplitudes in it, applied to a delta function on the grid, instead of a
narrow-angle Gaussian. And the grid is fourth order in height throughout: the derivative
(`--vertical_order 4`, Numerov's form), the ground condition, and the interpolation that puts
the source on the grid and takes the receivers off it (`--boundary_order 4`). Against the
exact two-ray field over soft ground at 1 kHz, from a 100 m source:

| Elevation | Gaussian starter, second order, 20 points | self-starter, fourth order, 10 points | … 20 points |
| --- | --- | --- | --- |
| 0–20° | 0.01–0.02 dB rms | 0.005 | 0.005 |
| 20–30° | 0.14 | 0.006 | 0.003 |
| 30–35° | 0.59 | 0.01 | 0.002 |
| 35–40° | 0.61 | 0.02 | 0.004 |
| 40–50° | 1.75 | 0.01–0.07 | 0.001–0.007 |

Fourth order is also what lets the grid be coarse. It is 10 points per wavelength by default
(`--ppw`; 20 under `--scheme cn` or with either order set to 2). On
`atmospheres/example_windy.csv` out to 1.5 km, 250 Hz – 2 kHz, against an 80-point run:

| | 20 points | 10 points | 6 points |
| --- | --- | --- | --- |
| second order throughout | 0.02–0.04 dB rms (0.14 worst) | 0.06–0.19 (0.64) | |
| derivative fourth, ground and interpolation second | 0.01–0.02 (0.06) | 0.05–0.09 (0.29) | |
| fourth order throughout | under 0.0005 (0.001) | 0.001–0.002 (0.011) | 0.004–0.01 (0.08) |

The middle row is why the ground condition and the interpolation had to follow the
derivative: with those at second order the whole scheme converged at second order.

The self-starter also corrects a bias at every angle. The Gaussian is normalized with the
wavenumber at the ground, so a source aloft in a refracting or windy column comes out at the
wrong strength: 0.11–0.15 dB high against the exact residue series for a 100 m source where
the sound speed is 3% lower than at the ground, and 0.10 dB on `example_windy.csv`. The
self-starter is within 0.04 dB.

`--scheme cn` is the march this replaced, Crank–Nicolson on the Padé(1,1) wide-angle equation
at ten steps per wavelength (`--range_ppw`), kept as an independent check. It is good to about
30°, and its phase error builds with range. Through `atmospheres/example_windy.csv` from a
100 m source out to 1.5 km, 125 Hz – 2 kHz, the two differ by 0.03–0.12 dB rms with the same
starter and grid order; against the exact field in still air it is Crank–Nicolson that is off,
by 0.7 dB rms at 25–30° and 2.1 dB at 30–35° at 1 kHz.

| Five bands, 125 Hz – 2 kHz, 100 m source to 1.5 km, 14 cores | PE time |
| --- | --- |
| `--scheme cn` | 25.4 s |
| `--scheme pade` (default: 50°, self-starter, fourth order, 10 points) | 0.47 s |
| … `--max_angle 60` | 0.67 s |
| … `--max_angle 35` | 0.37 s |
| … `--ppw 20` | 1.4 s |
| … `--ppw 20 --max_angle 35 --vertical_order 2 --boundary_order 2 --starter gaussian` | 0.6 s |

(The ray models' columns take another 0.4 s.) One march is 39 times faster than
Crank–Nicolson at the defaults. Bands and their frequencies run in parallel. For
Crank–Nicolson, against a
30 × 30 Crank–Nicolson reference, 10 in height is within about 0.1 dB rms for a receiver near the ground, but
0.47 dB rms (2.2 dB at worst) for one 1.2 m up at 1 kHz; 20 brings that to 0.09 dB, at twice the
cost. On the stored shadow-calibration cases the change moved the PE by 0.16–0.39 dB rms (1.1 dB
at worst) and left the calibration's edge error against it unchanged.

The PE is checked against exact solutions:

| Case | Agreement | Test |
| --- | --- | --- |
| Two rays over rigid ground | within 0.3 dB | `tests/test_pe.cpp` |
| Weyl–van der Pol field over grass | within 0.5 dB | `tests/test_pe.cpp` |
| Residue series in an upward-refraction shadow, rigid ground | within 0.02 dB at 250, 500 and 1000 Hz, down to 47 dB below free field | `tests/test_pe_shadow.cpp` |
| Benchmark case 1: direct and image rays over the benchmark's ground | within 0.2 dB at 10 Hz from 50 m (0.29 dB at 25 m, inside a wavelength of the source) and 0.08 dB at 100 Hz and 1 kHz from 20 m, to 1 km | `tests/test_benchmark.cpp` |
| Benchmark case 2: normal-mode series, downward refraction | within 0.06 dB at 100 Hz to 1 km and 1 kHz to 500 m, 0.2 dB at 10 Hz | `tests/test_benchmark.cpp` |
| Benchmark case 3: residue series, upward refraction, impedance ground | within 0.03 dB at 100 Hz and 1 kHz, to 60 dB below free field | `tests/test_benchmark.cpp` |

The shadow rows are the checks that matter for arbitrating the ray models. Their profiles
have k² linear in height, and the field in the shadow is then exactly a sum of creeping waves
built from Airy functions of complex argument. Doubling the grid resolution moves no shadow
onset and no level past 100 m by more than 0.7 dB.

```sh
niceops_pe_compare --atmosphere atmospheres/example_windy.csv --toward 270 \
    --source_height 328 --receiver_height 3.937 --bands 125 250 500 1000 --max_range 8000
```

Through `atmospheres/example_windy.csv` with a 1.2 m receiver, this is where each model falls
10 dB below free field and stays there, upwind:

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

And the rms difference from the PE beyond 400 m, downwind of the 100 m source:

| band | straight | linear | stratified |
| --- | --- | --- | --- |
| 125 Hz | 2.8 dB | 4.0 dB | 0.9 dB |
| 250 Hz | 7.1 dB | 7.6 dB | 0.4 dB |
| 500 Hz | 3.1 dB | 0.9 dB | 0.9 dB |
| 1 kHz | 8.1 dB | 7.2 dB | 2.1 dB |

What that says:

- **Downwind and in the lit region, `stratified` is right** and the others are not.
- **Upwind, the shadow begins between the two ray models' predictions.** It moves toward
  `stratified` as the frequency rises and the wavelength shrinks against the near-ground shear
  layer, which is what diffraction would do. `stratified` starts the shadow too early and
  `linear` too late, by similar amounts.

`tools/plot_pe_compare.py OUTDIR CSV:TITLE ...` plots that output: one figure per case, a
panel per band, and a summary of where each shadow begins. It runs under the build's Python,
`build/pyenv/bin/python`.

**The shadow's edge.** Over the stretch where the PE falls from −3 to −25 dB, the ray models'
level error against it (rms, dB):

| case | before | with the boundary field |
| --- | --- | --- |
| linear lapse (−0.018 /s), sources 30 and 100 m, 125 Hz – 2 kHz | 5–20 | 0.1–1.0 |
| log-law winds, `stratified`, 16 cases, 125 Hz – 2 kHz, pooled | 30.7 | 6.1 |

| log-law winds, pooled rms | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz |
| --- | --- | --- | --- | --- | --- |
| `stratified` before | 7.1 | 22.6 | 34.6 | 37.4 | 39.3 |
| with the boundary field and the per-band boundary | 6.0 | 5.2 | 5.7 | 4.7 | 8.4 |

On the linear profile, where the ray model's geometry is right, the boundary field puts it on
the PE. On log laws it needs the per-band boundary too. Scaled by the layered model's own
boundary, it falls 10–20 dB too fast at 125–250 Hz and pooled at 19.5 dB.

For a receiver on the ground (3 mm up) the same comparison gives 7.1 dB pooled, against 10.0
for `linear`, with the lapse control at 0.5 dB. By band, 125 Hz is the worst (10.2 dB) and
500 Hz the best (4.4).

Two things that looked like the treatment and are not:
- **Clipping the profile below a fixed height.** Once the height passes the microphone, the
  microphone sits in the clipped layer and shadows at every range.
- **Smoothing the profile per band.** A linear profile smoothed only over its lowest meters
  moved the lapse control by 20 dB at 2.5 km. The ray model is too sensitive to structure
  below the microphone for either to work.

The PE's deep shadow in still air, 40–80 dB down, is not to be believed either: with no
turbulence in it, nothing scatters energy in.

**Band sample counts.** PE band sample counts must be positive. Zero or negative counts are
rejected by the solver API, including for an empty batch, rather than producing invalid band
levels. Valid sample counts retain the existing midpoint quadrature.

### Turbulence, and the energy it scatters into shadows

The PE takes turbulence as random phase screens (Gilbert, Raspet and Di 1990): `--gamma_t`
(HARMONOISE strength), `--outer_scale` (von Kármán, default 10 m) and `--seed`, with each of
`--samples` frequencies its own realization, so a band level averages over the band and the
ensemble at once. The effective index's structure parameter is γT / 4, the value for which
Ostashev's coherence and HARMONOISE's both reduce to one form (Attenborough & Taherzadeh 2026,
eqs. 18 and 20). One screen per outer scale of range carries the phase that slab accumulates.
Its vertical spectrum is the 3-D one integrated across the plane, so the vertical coherence
is the 3-D medium's. Checks, `tests/test_pe_turbulence.cpp`:
- **The screens' variance:** 0.782 k² Δx (γT/4) L0^(5/3), to 0.06%, less the 4% the finite
  column leaves out.
- **The ensemble-mean field:** it falls as exp(−σ²/2) with range to within the sampling noise
  of 300 realizations, which checks the screens and the march together.
- **Single scattering:** a third of the strength lowers the deep-shadow level by 4.82 dB,
  against 10 lg 3 = 4.77.

At γT = 3e-5 (what the 2017 co-located pole and board pairs fit, in panam), the PE fills the
day-284 shadows to −3 dB re free field at a board, median, only 5–10 dB below the lit level.
`--scattering` adds that energy to the shadow-boundary field, in place of the floor:
S = a + 10 lg γT + b lg(f / 1 kHz) + c lg(R / 100 m) + d lg(z_s / 100 m) dB re free field.

| receiver | a | b | c | d | fitted on |
| --- | --- | --- | --- | --- | --- |
| flush (3 mm, or a plate) | 39.18 | −1.55 | 2.10 | 5.30 | 12 profiles, sources 50 and 120 m, 1.5 dB rms |
| 1.2 m | 36.67 | +5.61 | 5.13 | 5.30 (the flush fit's) | 6 profiles, sources 50 m |

S is held to 0 dB, R to 100 m–10 km and z_s to 10–300 m; a plate lifts it by K as it lifts
the coherent field.

**The measurements do not want that much.** Against the day-284 SEL, the upwind 5–10° bin
goes +0.28 / 3.06 dB (median / rms) with the floor, +0.49 / 3.04 at γT = 1e-6, +0.95 / 3.02 at
3e-6, +1.55 / 3.16 at 1e-5 and +2.96 / 4.00 at 3e-5. Up to about 3e-6
is consistent with them, which is HARMONOISE's "moderate" range. The pair fit's 3e-5 is ten
times that. It is three times HARMONOISE's maximum, and probably carries more than turbulence
(roughness, the pole's height), or turbulence at the pole's height, where it is strongest.
`--scattering` stays off by default, since the floor does as well on these data.

### The Attenborough et al. benchmark

Attenborough et al., "Benchmark cases for outdoor sound propagation models", JASA 97, 173
(1995), is the standard set of test cases for outdoor propagation codes. It fixes a source at
5 m, a receiver at 1 m, 10, 100 and 1000 Hz, a locally reacting ground with impedances it
tabulates, absorption at 20 °C and 70% RH, and four profiles:

| case | profile | exact solution |
| --- | --- | --- |
| 1 | uniform | direct and image rays with the spherical-wave reflection coefficient (their eq. 12) |
| 2 | c rising 0.1 /s (downward refraction) | normal-mode series (eq. 10) |
| 3 | c falling 0.1 /s (upward refraction) | residue series of creeping waves (eq. 11) |
| 4 | rising 0.1 /s to 100 m, falling to 300 m, uniform above (a duct) | none |

`niceops_benchmark` runs all four to 200 m and to 10 km with the PE and, where one exists,
the exact solution (`src/exact_solutions.hpp`). `tools/plot_benchmark.py DIR` draws them in
the paper's layout, so each figure sits beside their figs. 11–14, and prints how far the PE is
from the exact solution. A full run takes about 10 seconds.

```sh
niceops_benchmark --out figures/benchmark
build/pyenv/bin/python tools/plot_benchmark.py figures/benchmark
```

Set beside the paper's FFP and PE curves, all four cases agree in every panel. Where there is
an exact solution, the PE is within 0.3 dB of it past 20 m on every 200 m panel. On the 10 km
panels it is within 0.3 dB for case 1, and for case 3 within 0.1 dB, or 0.5 dB at 10 Hz,
where the linear profile and the exact solution's profile part company within a few Airy
lengths of the ground. Case 2 at 10 km is the exception, covered below. Case 4 has no exact
solution; it reproduces the paper's curve feature for feature, down to the 100 Hz humps at
2.5–3.2 km and 4.3–4.7 km.

Three things came out of it:

- **Their eq. (11), the upward-refraction residue series, has the ground condition wrong.**
  It writes the eigenvalue condition as Ai'(b) + i k0 l beta w Ai(b) = 0, w = exp(2 pi i/3),
  where the derivation gives Ai'(b) - i k0 l beta w² Ai(b) = 0: the admittance term is off by
  a factor exp(i pi/3). Over the benchmark's 100 Hz ground the printed form is 1.5 dB off the
  PE at 100 m and 12 dB off at 1 km. The derived form agrees with the PE to 0.02 dB down to
  -90 dB. Over rigid ground the two coincide, which is why that case never showed it.
  `src/exact_solutions.cpp` carries the derivation.
- **Case 2's normal-mode solution carries a field no PE can.** Its profile, c0/sqrt(1 - 2z/R),
  reaches infinite sound speed at z = R/2 = 1716 m, so every ray turns below that height,
  including the near-vertical ones. Those fall back within R = 3.4 km of the source, landing
  at R sin(2 × elevation), like mortar rounds. They put structure at 3.4 and 6.8 km into the
  mode sum; the paper's own normal-mode figure has it, and its FFP and PE figure does not.
  With modes steeper than 75° left out (`max_elevation`), the series and the PE agree to
  0.06 dB out to 1 km at 100 Hz. With them in, they differ by 5 dB.
- **Case 2 at 10 km depends on where the profile is cut off.** The paper says so, and it holds
  here: the PE with the linear profile stopped at 1000 m matches their curves in form, not
  interference null for null.

The inputs are theirs. The impedances are the four-parameter model's values they give, not a
NICE-OPS impedance model. The absorption is the 0.18 dB/km at 100 Hz and 5.4 dB/km at 1 kHz
they quote from ANSI S1.26-1978. ISO 9613-1 gives 0.22 and 4.98 dB/km at the same conditions.

## Performance

### How much gets evaluated

The sound exposure integral covers the window from each observer's first to its last sample
within 10 dB of its peak, dips included (the duration of the aircraft noise metrics), which on a
shipped case is a little under 60% of the track. The observer loop therefore does not march every
emission time: it finds the peak with a per-direction bound on the level, and evaluates only the
states that bound cannot rule out of the window. The run logs how much it evaluated:

```
Evaluated N of M (P%); the rest is provably outside every observer's 10 dB window.
```

It is exact, not an approximation — the output is byte-identical to evaluating every state. It
is used on the spectral path (`--spectral`, `--ground_effect`, plates, refraction, tables), where
an evaluation costs 150 ns to a microsecond; a broadband run, where an evaluation costs about 20
ns and the bound proves only about 6% of the pairs outside the window, marches every state
instead. It is skipped, with every state evaluated, for `-a`, `-i`, `--band_history`,
`--epnl_file`, `--band_frames`, `-D`, `--lden_file`, `--cnel_file`, `--lamax_file`,
`--nax_file`, `--ta_file`, `--csel_file`, `--full_history` and `--db_cutoff` as well; see
[`CLAUDE.md`](../CLAUDE.md) for why each one cannot use it.

### `--db_cutoff`

On a large observer grid most of the work goes into observers that are never loud enough to
matter. `--db_cutoff L` skips any evaluation that provably falls below `L`, bounding the level
from the loudest sphere the state blends and spherical spreading alone.

The threshold sits 10 dB below `L` and is applied to the **instantaneous level**, which makes the
result exact for any observer whose *peak level* reaches `L`: the SEL integral covers the window
from the peak down to peak − 10, so no sample in that window is ever culled.

The cells of the output are SELs, though, and an SEL exceeds its observer's peak level by 1 to 12
dB (the window's duration), so a cell can be at or above `L` in SEL while its peak is below `L`,
and its window then loses samples. On the demo case, `--db_cutoff 85` changes 157 of the 592 cells
whose SEL is at or above 85 dB, by at most 0.43 dB; a cell is unchanged to the last bit once it is
about 3 dB above the cutoff. Read the cutoff as a threshold on peak level, and leave a few dB of
margin when it is the SEL you care about.

Observers that never reach the cutoff come out as `-inf` — they contributed no energy, and those
cells will not plot. That is the trade the flag asks for; it is off by default.

`--db_cutoff` works on the spectral path through per-direction bounds on the band-limited,
weighted level (`Database::spectral_direction_bounds`). With `--ground_effect` or
`--plate_table`, each band is raised by the most the ground or the plate can add.

### Threads

The observer loop runs on every core by default. Each observer is independent, so the output is
deterministic: runs at `-j 1`, `-j 2` and `-j 32` are byte-for-byte identical.

Observers are dealt out to the workers in chunks as they finish rather than split evenly up
front, because on a machine with both performance and efficiency cores — every Apple silicon
part, and most recent x86 laptops — an even split makes the slowest core the critical path.
On an M4 Pro the even split ran the demo case no faster on all 14 threads than on its 10
performance cores; pulling chunks recovers about 10%.

### Cost of the spectral and ground paths

**The spectral and ground paths are expensive.** Acoustic loop on the Be407 case (94.3 M evaluations, 26 bands, 14 threads,
default band limits): broadband 0.34 s, `--spectral` 1.44 s, `--spectral --ground_effect`
10.9 s. So carrying the spectrum costs about 4x, and the ground effect about 32x — a complex
reflection coefficient and a Faddeeva evaluation per band per observer per state, none of
which can be hoisted because it all depends on the geometry. The Faddeeva evaluation is about
half of the ground path on its own.

Those three are from before the September 2026 optimization round and all of them are now
lower. Measured against the same build on the S-76D demo track (94.3 M evaluations, 14
threads), over a synthetic single-condition 24-band database built the way
`write_flat_spectral_database` in `tests/test_spectral_propagation.cpp` builds one:

| | before | after |
| --- | --- | --- |
| broadband | 0.32 s | 0.16 s |
| `--spectral` | 1.36 s | 1.16 s |
| `--spectral --ground_effect` | 8.95 s | 7.81 s |

A single condition makes the operating-condition blend degenerate, so those two spectral
figures are optimistic against a real multi-condition database like the Be407 one in the row
above; the ground-effect term, which is most of the last row, does not depend on the condition
count. All three outputs agree with the previous build to 1.4e-12 dB.

## Plotting

Plotting embeds Python through pybind11. The first configure creates a virtual environment at
`build/pyenv` and installs matplotlib into it, and the embedded interpreter puts that directory
at the front of `sys.path` — so a run depends on what the build made rather than on whatever
happens to be installed for the system interpreter, and needs no `PYTHONPATH`. The install
happens once per build directory and needs network access; afterwards configures are offline.

Two options control this:

| Option | Meaning |
| --- | --- |
| `-DNICEOPS_PLOTTING=OFF` | Drop plotting, and Python, pybind11 and matplotlib with it. The acoustics need none of them. `-p`, `-f` and `-F` are then refused with a message rather than silently doing nothing. |
| `-DNICEOPS_BUNDLE_PYTHON=OFF` | Keep plotting but skip the virtual environment, using the interpreter's own packages instead. |

Headless is handled: the facade selects the `Agg` backend itself.

The reference figures under `figures/` date from June 2020 and were produced with three settings
that differ from the defaults. All three are now reachable:

```sh
./build/niceops -t tracks/S-76D/197401AC.csv -v databases/S-76D.json -d databases/S-76D_M3.nod \
        --trim -6000 --grid -2000 1100 -1500 1600 -o out.csv -f out.eps
```

Matched that way the two runs land on the same 124 × 124 grid and agree cell for cell: the largest
difference anywhere is **1.05 dB**, 99% of cells are within 0.5 dB, and the peak levels differ by
0.07 dB. Contour levels now default to 65–100, the values the references were drawn with; the old
45–90 put its top level below the peak of the demo case, leaving the loudest patch unfilled.

## Changes that affect results

Several corrections have moved NICE-OPS's results. Compare old and new results only with
these in mind. The largest is the first:

- **The 4x EAA multiplier** (September 2026), below: the demo footprint rose by 9.93 dB mean.
- **`-a`** adds the ambient level in energy. It used to be added as pressure amplitude, as the
  D original did, which read up to 3 dB high ([Command-line options](#command-line-options)).
- **`-D`** is normalized to a 24-hour day. It used to write the raw sum, 49.37 dB high
  ([Several trajectories, and the day-night level](#several-trajectories-and-the-day-night-level)).
- **`--start_utc`** is required with `-D` for a WGS84 track without a clock. It used to
  default to midnight, so the whole footprint took the night penalty.
- **`--lamax_file`** writes LAFmax, the IEC 61672 fast-weighted maximum in reception time, by
  default. It used to write the instantaneous peak of the history, 0.7–0.9 dB higher at the
  median; `--lamax_weighting instantaneous` still does, to the bit. `--nax_file` counts against
  the same weighted maximum
  ([Maximum level, time above and C-weighted SEL](#maximum-level-time-above-and-c-weighted-sel)).
  A trajectory sampled more coarsely than half the weighting's time constant is refused, so
  WGS84 waypoints at the default `-s 0.5` now need `-s 0.0625` or finer for `--lamax_file` or
  `--nax_file` (or `--lamax_weighting instantaneous` for the old number).
- **`--full_history`** now governs `-D` as well as SEL. It used to leave DNL on each operation's
  10 dB-down window, 0.3–0.6 dB below the whole history at the median; the default is unchanged
  ([Several trajectories, and the day-night level](#several-trajectories-and-the-day-night-level)).
- **The point receiver's microphone height** sets the evaluation point. Point-receiver SEL
  moved by up to about ±0.8 dB ([Receivers](#receivers)).
- **Four length options** (`--grid`, `--gridstep`, `--receiver_height`, `--path_smoothing`)
  were meters before `--units` existed ([Units](#units)).
- **Stratified propagation keeps every arrival** it solves, where it kept 24
  ([Refraction through the profile's layers](#refraction-through-the-profiles-layers)).
- **The hover nose is turned by the frame** whenever the run knows one, not by a fit to the
  track that absorbed its mean crab. Only tracks that fly fast and then hover move
  ([The sphere's nose](#the-spheres-nose)).
- **`--refraction stratified` resolves a calm floor** only as finely as the run's highest band
  can. A profile with calm, isothermal air at the ground (a log law written below its roughness
  length) used to put the shadow boundary kilometers out and flush microphones 22–28 dB low;
  profiles without such a floor move by under 0.001 dB
  ([Refraction through the profile's layers](#refraction-through-the-profiles-layers)).

### The 4x EAA multiplier

**Results from before September 2026 are not comparable to results after.** `eval` used to
multiply the database's `EAA` by a constant 4.0, under a `//TODO Add ground effects` inherited
from the D original. That is not a ground effect, and it overstated atmospheric absorption
fourfold — 168 dB at 6 km where the physical value is 14. panam *writes* that `EAA` as
`SPLA - SPLAa` over a 1000 m reference and consumes it as `EAA / distance * (slant_range -
radius)`; the tool now does the same.

Removing it raised the demo footprint by **9.93 dB mean, up to 22.7 dB**, most of it far from
the track where the old factor was removing the most. `golden/197421_LAE.csv` was regenerated;
`golden/state_197421.csv` is trajectory state and unaffected.

It was measured against `--ground_effect`, the model it was standing in for, before being
removed: soft ground costs a roughly range-independent 2–4 dB relative to a rigid surface,
where the factor removed 7 dB at 250 m rising to 38 dB at 2 km. A real ground effect is bounded
by the surface impedance; a multiplier on a distance-proportional absorption is not.
