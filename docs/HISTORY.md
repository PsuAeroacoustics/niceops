# NICEOPS history

How NICEOPS got to its current form, and why old results differ from new ones. For what the
tool does now, see the [README](../README.md) and the [manual](MANUAL.md); for how accurate it
is, see [FIDELITY](FIDELITY.md).

## A step's split made opt-in (October 2026)

The jump split across a step's line (below) is taken only with `--tilted_step_split` now, and a
run's levels are again those before it: the level steps where a microphone crosses the line. Its two
constants were fitted over one step; on a second it gained less overall and was worse above the
line at 125 and 250 Hz (FIDELITY, "Across a step"), it roughly doubles the cost of the microphones
it reaches, and no validated scene moved by more than 0.03 dB with it. It removes the seam in maps.

## A step's split reaches past the edges and the floor; a dip may be split between footprints (October 2026)

The jump split across a step's line stopped where the microphone's foot left the edges or went
under the floor, and the level stepped there instead: by 1.0-2.1 dB across the line just past the
blocks' ends, and 1.8 dB where the line meets the floor. It now goes on as far as the pair's rays
reach, and its foot slides back along the line to the floor; against boundary elements by the
floor the error went from 2.0-2.4 dB to 0.3 at 125 Hz and from 1.1-1.4 to 0.65-0.78 at 250. Runs
summed in energy took the split's twin rays as added energy, up to 3.8 dB too loud near the line;
they take no split now. A dip whose low roof is two touching footprints of one height (10/8/8/10)
was left to the cascades, stepping 6.5-7.2 dB through the roofs' plane; it takes the profile of the
same dip in one piece. The stored tests, EDGE's LoD1 scenes against the FDTD, DigiCity and the
verify cases do not move.

## A dip takes its row's profile (October 2026)

A row of touching blocks with a roof lower than both its neighbours (a dip) was left to the
cascades of its roof edges, which stepped by up to 22 dB where a microphone crossed the roofs'
plane and spiked at roof height. Such a row now takes its profile field, as a row with a street
did already (THEORY §6.5, "Over a row of roofs"). Against boundary elements over a 10/8/10 m dip
the error went from 4.92 to 1.39 dB rms at 125-500 Hz, 2.70 to 1.51 at 1 kHz, and the step at the
plane to the exact field's own, 0.4 dB at most. Touching steps up or down keep their cascades.
Streets, the stored tests and DigiCity do not move; of the EDGE scenes only UE05's LoD1 stand-in
does (FIDELITY). `--no_profile_pe` gives the old field.

## A step's line no longer steps the level (October 2026)

Where a microphone crossed the line through a step's two edges (a ridge or a taller roof's edge
and the edge of the lower roof it steps down to) the point-source model's level stepped by up to
6.4 dB: below the line the pair was the receiver's aperture, above it the first edge's ray and the
cascade. The field's jump at the line is now split between its sides (THEORY §6.5, "The jump across
a step's line"), so the level is continuous there. Against boundary elements near such a line the
error went from 0.67 to 0.26 dB rms (352 cases) and from 1.21 to 0.47 (129 cases within a metre of
it); the stored step cases from 0.53 to 0.48, EDGE's terrace from 0.32 to 0.27. Far above the line
levels move by 0.04-0.08 dB rms. EDGE's LoD1 scenes and DigiCity do not move; EDGE's pitched
scenes move by 0.21 dB at most. It was made opt-in the same day (above).

## The urban model's search made faster, and a fast mode (October 2026)

Eleven exact changes to the point-source model's path search and band sum left every ray and band
level as it was, to the bit with clang and on x86 builds without FMA (CI's). A GCC build with FMA,
on arm64 or otherwise, moves the profile PE's band levels by about 4e-9 dB: its complex divides
are now written out as compiler-rt's, where libgcc's are Smith's method. The changes cut the search's time 5.4-5.8x on a pitched suburb and about 2x
on DigiCity. The largest was proving which part of an edge an image of two facades sees, where
such images used to take the whole edge. `niceops --buildings` now takes `--path_radius`, the
reach of each end's buildings, fixed at 150 m until now; below it is a fast mode, at a measured
cost (FIDELITY), and `--diffuse_radius` follows it by default.

The source-coherence pair search now switches from its z-sorted sweep to direct enumeration once
the sweep finds 20% of all pairs. Dense launch clusters made sorting and retaining candidates
slower than forming the full pair list; sparse sets still use the sweep. The threshold and its
microbenchmark are recorded in [performance notes](internals/performance.md#urban-source-coherence-pairs).

## A row's profile up to a size in wavelengths (October 2026)

A row's profile field, the parabolic equation over a row of flat-roofed blocks, used to be taken at
every band, at a cost that grows as the row's size in wavelengths squared: over a terrace of eight
rows under a low approach it cost 37 times the rays it stands for, a run there about nine times.
It now gives way to those rays between k D = 3000 and 6000 (D the row across), 1-2 kHz over that
terrace and above 5.5 kHz over a single street, and is not marched beyond; and its column stops at
what the march's end needs rather than reaching up to the source, which moves the field by 0.024 dB
rms. Together: 5.6 times the rays, a run about twice one without. Against the rays the profile moved
the terrace's A-weighted levels by -0.20 dB on average and now by -0.19; its worst pair moves from
-2.1 to -1.8 dB. `niceops --buildings` also takes `--no_profile_pe`. No golden moves: none crosses
a row.

## Through a roof's plane with both ends above it (October 2026)

With the source and the microphone both above a flat roof's plane, the ray diffracted over the
roof's two edges used to take the cascade, while the microphone just below the plane took the
receiver's aperture. With the source low over the roof, the far edge sits in the near edge's
transition, which the cascade cannot carry: the level stepped by 3-8 dB as the microphone crossed
the plane (11 dB with the source 0.3° over a 2.44 m block), and above the plane the cascade was up
to 4.3 dB low against boundary elements. The pair now takes each end's mirror below the plane,
negated (THEORY "Both ends above the face's plane"). Over 528 cases through the plane behind a step
the error goes from 0.88 to 0.20 dB rms and the worst from −4.3 to −1.0 dB; the step at the plane is
0.45 dB at most. A wall between two corners takes the same rule, with both ends in front of it:
on DigiCity's field 1.2 m up a quarter of the microphones move by more than 0.1 dB with a source
12 m up, by up to 4.4 dB, and the measured microphones by 0.04 dB at most; EDGE's scores do not move.
The facade-scattering
golden moves by 1.2e-5 dB at one observer; the others do not move. A step's tilted pair keeps the
cascade with both ends above its line, where the mirror does not hold.

## A shared face's boundaries go to one side (October 2026)

A microphone exactly in a roof's plane, behind a shared face or a step, or exactly on a step's
lower-roof reflection boundary, used to take both sides' waves at once: the pair's cascade from
one side and the blocked or reflected ray from the other. It read 3-9.6 dB high in the plane
against boundary elements and spiked by 13 dB on the boundary. Each now takes the side the ray set
gives it. Over the step's 72 cases above its lower roof the error goes from 2.34 to 0.38 dB rms.
A rooftop microphone beside a roof of the same height lies in that roof's plane.

## Twice-diffracted rays off a roof at both ends (October 2026)

A twice-diffracted ray may now reflect in a roof, level or a slope, at both of its ends; it used to
take one at most. Only microphones and sources over roofs see them: on the pitched suburb's
rooftop microphones, under a drone, the band levels move by 0.078 dB at most.

## The row's parabolic equation, faster and converged (October 2026)

The parabolic equation over a row of roofs now marches by split-step Padé and takes a row to
10 kHz in 0.29 s of CPU rather than 6.7. Its domain above the roofs used to be a fixed 15
wavelengths high with a layer as thick, which reflected the near-grazing waves at high frequency:
a street's parts were 1.1 dB off their converged values at 3.15 kHz and up to 3 dB at 8 kHz. The
domain now grows with the row's length, and the roofs sit on grid nodes, where the level used to
step by up to 0.2 dB as the frequency moved a roof between two. Against boundary elements one
street goes 1.65 → 1.61 dB rms and two 2.15 → 2.10. EDGE's UE12 and UE13 move with it (FIDELITY).

## A row of roofs by a parabolic equation (October 2026)

Over a row of flat-roofed blocks with a street between, crossed near grazing, the point model
used to sum the cascades over two, three and four roof edges, which do not converge there: over
one street they were 8 dB high at 500 Hz. Those rays now give way to the row's two-dimensional
field from a parabolic equation over its cross-section, carried as four rays (THEORY "Over a row
of roofs"). Against boundary elements one street goes from 3.4 to 1.65 dB rms and two from 6.7 to
2.15, the cascades having been as much as 14 dB off; at a microphone in the roofs' plane behind
the row, where the cascades stepped 11 dB, the field is now smooth. On EDGE only UE12's and UE13's
receiver R3 moves, toward the FDTD (UE12's floor source from 3.5 dB above it to 1.0 below). The
shipped goldens do not move. `--no_profile_pe` in point_field restores the cascades, for
comparison. A row costs seconds of CPU at high frequency (FIDELITY).

## Plane-wave tables over roofs (October 2026)

`niceops_urban_table` now takes a microphone on or over a roof, which it used to refuse: the roofs
below it reflect, specularly and to the diffracted rays, a microphone on a roof or a slope takes its
reflection there, `--roofs` reads a CSV district's pitched roofs, and the roofs scatter and light
the walls' diffuse field. A corner's end wave toward a microphone, or a facade's image of one, in
its solid sector in plan, which a prefilter used to drop before the exact tests, is kept as the
point-source model keeps it. On the test district's run that moves the levels by 0.001 dB at most,
and the roofs' scattering by 0.005 dB more with `--facade_scattering`. `--no_roofs` restores the
refusal, for comparison.

## Roofs light the scattering walls (October 2026)

With scattering, a point source now lights the walls' and roofs' patches off a level roof or a
slope as well as off the walls and the floor; those paths used to light nothing. Diffuse levels
can only rise: on the gable suburb by 0.04 dB on average, at most 0.16 dB at s = 0.3.
`--no_roof_lighting` in point_field leaves them out, for comparison.

## A shared face with an end at its line (October 2026)

Where one end of a ray over a shared face lay within about 1e-5 rad of the line between the edges
and the other above it, the aperture's field lost its precision and grew, to 6-13 dB over the
direct wave; it is now taken at its limit. Only such ends move: on the gable suburb, microphones
1.5 m up that lie exactly in a slope's extended plane, by up to 10 dB, which used to depend on the
compiler (GCC against clang) and on the device (the GPU 3.4 dB off the CPU).

## Twice-diffracted rays off a roof (October 2026)

A twice-diffracted ray now reflects in a level roof or a slope between one end and its edge,
where it used to be traced only off the floor. These are the partners a once-diffracted ray's
roof reflection needs where a second edge cuts it off; without them the field stepped there by
the whole ray. Results with roofs move accordingly: on the gable suburb up to 2.4 dB at 125 Hz
(mean +0.07 dB), on DigiCity's flyover up to 3.1 dB in 1.6% of band levels. `--no_twice_roofs`
in point_field leaves them out, for comparison.

## Walls with sloped-roof reflections (October 2026)

The point model now traces specular paths through a slope between two facade chains, and
once-diffracted paths with facade reflections between an endpoint and its slope reflection.
These paths were previously absent. Their floor partners record real reflection points and
retain a downward launch for a floor-contact source. Specular source images are traced even
when diffraction is disabled. The comparison switch is `--no_slope_facades` in point_field;
`--unending_edges` also exposes the existing finite-edge comparison there.

## Lateral-offset validation fixes (October 2026)

The fold check for WGS84 tracks with velocity columns mixed recorded velocity with
acceleration from fitted positions. Changing only the recorded speed could therefore admit
an offset beyond a turn's radius or refuse a valid offset. Both derivatives now come from
the fitted positions; recorded velocity still sets the offset direction.

Across a slow interval, choosing the nearest faster sample's direction could switch the
offset sideways when the aircraft stopped and departed in another direction. Those offsets
are now refused rather than fitted into a path with artificial speed and acceleration.
Leading and trailing slow intervals still hold the available direction, and zero offsets
retain the original path.

## Origins (2020–2021)

Eric Greenwood and Rob Rau developed NICEOPS in D, starting in January 2020. It used a
CGAL-based barycentric interpolator, the SPLINTER B-spline library and matplotlib-d. The S-76D
case, with `S-76D_M3.nod`, became the default in March 2020. The AW139, B205 and B407 databases
and tracks from the joint NASA/FAA flight tests were added in January 2021, together with a
second S-76D database. These are the `_M` files in `databases/`: `AW139_M1`, `B205_M1`,
`B407_M2`, `S-76D_M3` and `S-76D_M4`, plus `B407_M3` from the same family, which was added to
the repository later. All are broadband, with discrete load factors.

The June 2020 reference figures came from that code; see
[Reference figures (June 2020)](#reference-figures-june-2020).

## The C++17 port (September 2026)

In September 2026 NICEOPS was ported to C++17 (`666046b`) and checked against the D reference.
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
consumes it as `EAA / distance * (slant_range - radius)`, with no such factor. NICEOPS does the
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

Besides the departures from D above, these corrections have moved NICEOPS's results. Compare
old and new results only with them in mind.

### Metrics

- **`-D` is normalized to a 24-hour day.** It used to write the raw sum, which read
  10 lg 86 400 = 49.37 dB high
  ([Several trajectories, and the day-night level](MANUAL.md#several-trajectories-and-the-day-night-level)).
- **`--start_utc` is required** with `-D`, `--lden_file` or `--cnel_file` for a WGS84 track
  without a `utcsec` column. It used to default to midnight, which is night, so the whole
  footprint took the 10 dB penalty.
- **A static operation needs `--start_utc`** for the day-night metrics. Its states carried a
  clock of zero, so `-D`, `--lden_file`, `--cnel_file` and the period levels dated every hold at
  UTC midnight without a word, which is night for `-D` (+10 dB); `--start_utc` was ignored for it.
  A hold now takes its clock from `--start_utc` and is refused without one. No golden or
  validation result moved: no shipped test or harness run asked for a day metric on a hold.
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
  |w| = 500, and NICEOPS followed, until panam `e978703` evaluated it by `wofz` at any |w|.
  NICEOPS does the same, pinned to panam `b758e3a`. The reflection coefficient Q moves by up to
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

- **`--refraction stratified` interpolates the layered rays between every fourth state** by
  default (`--ray_stride 4`; October 2026). Each path's rays were solved at every trajectory
  state; now at every fourth and the last, with the states between interpolated where the two
  agree and solved where they do not. Against every state solved, SEL moves by up to 0.012 dB on
  the windy demo profile and 0.057 dB in the inversion duct (99% of cells within 0.003 dB),
  LAmax by up to 0.002 dB on a validation approach, and single history samples by up to about
  1 dB where the solved history steps between states
  ([FIDELITY](FIDELITY.md#interpolated-rays---ray_stride)). `--ray_stride 1` reproduces the old
  results to the bit; the stratified goldens made before are kept as `*_exact.csv` and pin it.

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

### Urban receivers

- **The ideal canyon traces the floor's rays.** Its first tables left the floor to the ground
  model, as though the floor's reflection shared the direct ray's fate; at 15 and 30 degrees over
  a 20 m canyon they read 20-25 dB of shielding where the wave solution has 2-4 dB. Tracing them
  moved the demo canyons' low-pass microphone from −3.7 to −2.1 dB against open ground.
- **UTD replaced the CNOSSOS-EU placeholder at the roof edges**, and the table became coherent.
  The placeholder left shadows 2-4 dB too deep (−14.5 dB against the wave solution's −10.9 in
  one), and the energy sum jumped at every boundary. The demo canyons' range went from −2.1 to
  +2.1 dB to −1.0 to +2.2 dB, and a 400 Hz wave band, which had moved them by up to 1.1 dB, now
  moves them by 0.08.
- **Twice-diffracted rays whose stationary point does not converge are rejected.** The point
  was found by plain alternation, which stopped at 200 steps whether or not it had converged.
  Near a corner where a roof edge meets a vertical one it was still creeping toward the corner,
  and those iterates were kept as rays a few centimetres long between the edges, strong for
  their spreading: 0.1-0.2% of the point-source model's rays. Rejecting them moved DigiCity's
  four microphones by at most 0.063 dB, the 49-block district's levels by 0.07 dB at the 99th
  percentile and 0.37 dB at most, and 10 of DigiCity's 1427 grid points, in deep shadow west of
  the city (−58 to −63 dB re 1 m), by 1-5.6 dB.
- **The shared-face aperture went from the midpoint rule at a 24th of a wavelength to Filon's
  rule**, which moved the point-source model's levels by 0.015 dB at most. With the rejection
  above, the table builder's seven test microphones moved by at most 0.024 dB.
- **A receiver exactly on a shadow or reflection boundary takes the edge's term from the lit
  side.** On a boundary the coefficient's limit was taken a hair above the boundary's angle
  whichever side that was, which is the shadow side for the incident wave's boundary and face 0's
  reflection; the path tests keep a grazing ray there, so the two added, 1 + 1/2 where the field
  is 1/2, about 9 dB high. Only receivers exactly on a boundary moved, which round-numbered
  geometry puts many on: the 49-block district's 20 test microphones by 0.015 dB at the median and
  4.0 dB at most (each now nearer the mean of its neighbours 1 cm away: 0.034 dB at the median,
  from 0.13), DigiCity's grid by 1e-8 dB, the urban tables (which step off exact boundaries) not
  at all.
- **The point-source model absorbs the air**, ISO 9613-1 along each ray in each sub-band, and
  `niceops_point_field` needs `--humidity` (or `--no_air_absorption`). It had no absorption at
  all. At 70% RH and 20 °C, DigiCity's 22 microphones with the loudspeaker at hover spot 5 fall
  by up to 0.22 dB at 1 kHz, 0.86 dB at 4 kHz and 1.8 dB at 6.3 kHz; on the 49-block district,
  microphones 240-290 m from a drone 30 m up fall by 1.2-1.9 dB at 1 kHz and 5-9 dB at 4 kHz.
- **A twice-diffracted ray over edges that are not parallel spreads by the second edge's own
  curvature.** Its spreading was the parallel edges' $1/\sqrt{s_1 s_2 s_3 (s_1 + s_2 + s_3)}$
  whatever the edges' angle. The wave reaching the second edge is curved about the first edge
  across it and about the source along it, and past the second edge it spreads by its curvature
  along that edge, which for perpendicular edges is the first edge's, $s_2$, not $s_1 + s_2$
  ([THEORY §6.5](THEORY.md#65-urban-receivers)). Over a screen's top and a perpendicular
  screen's edge the ray was 1.45 dB high against the edge-source integral, and is within 0.007 dB
  of it. DigiCity's four microphones moved by at most 0.12 dB.
- **`niceops_point_field --edge_source`** takes the diffracted rays from the edge-source integral,
  which keeps the edges' ends. Without it nothing changes.
- **`niceops_point_field --scattering` and `--scattering_spectrum`** let walls reflect part of what
  reaches them diffusely, and read a building's `scattering` column, which the tool used to ignore.
  With every coefficient 0 nothing changes.
- **`niceops_fdtd_field`** is new: the full-wave reference for a point source among buildings.
- **`niceops_point_field`'s roofs reflect**, as the floor does, on specular rays and on either side
  of a once-diffracted one. Only what stands above a roof hears it, directly or by a higher edge:
  over EDGE's rows of unequal houses the drone's levels moved by up to 3.4 dB in a band (UE07),
  and nothing at a microphone below every roof moved. `--no_roofs` leaves them out.
- **`niceops --buildings`** is new: every microphone hears the source's own rays among the
  buildings, each reading the sphere along its launch, in place of the open-ground model.
- **`--buildings` adds rays leaving the source apart in energy** (`--ray_coherence cone`, the
  default), except those within `--coherence_angle` (1°) of each other. It first summed every ray
  with its phase, as if the rotor radiated in phase in every direction, which
  `--ray_coherence coherent` keeps to the bit. The DigiCity hovers moved by 0.77, 0.44 and
  0.22 dB rms at 15, 30 and 60 m, and the demo approach's pass SEL by at most 0.08 dB. A database
  may record its vehicle's coherence (panam's `directional_coherence`), which the run then takes
  (`measured`), with the main rotor's delay where it gives one. In between, the measured model
  was the default, with the delay switched on by the vehicle file's blade count and sense; a delay
  is an assumption about the source's phase, so it is now the database's to make. The level's
  spread over the unknown phases was first studied with every group of rays within a launch-angle
  window taking a uniformly random phase (behind environment switches, now gone), which put a wide
  spread on broadband bands that add in energy; `--phase_interval` takes the measured spread
  instead.
- **The UTD transition function comes from a table below X = 50**, fitted to the Faddeeva form it
  replaced at Chebyshev nodes; it agrees with that form to 7e-12, its own error, and the demo
  approach's `--buildings` levels moved by under 1e-9 dB.
- **The shared faces are sampled** (`niceops --buildings` and `niceops_point_field`): the uniform
  form's integral only where the cascade cannot stand in for it, at wavenumbers two thirds of an
  octave apart, the cascade carrying it elsewhere. It was taken at every wavenumber, which
  `--exact_shared_faces` keeps to the bit. The demo approach's states moved by 0.002 dB at most and
  its SELs by under 0.0001; `niceops_point_field`'s test district by up to 0.1 dB in a band.
- **`niceops --buildings` samples the states** (`--state_tolerance`, 0.02 dB by default): the
  buildings' transfer is evaluated at some states and interpolated between them, the source exact
  at every state. Every state was evaluated, which `--state_tolerance 0` keeps to the bit. The demo
  approach's SELs moved by 0.003 dB at most and its LAmax by 0.007, for a third of the pairs.
- **The GPU is the default wherever there is one**: `niceops --buildings` (new there),
  `niceops_point_field`, `niceops_urban_table --wave_below` and `niceops_fdtd_field`, each falling
  back to the CPU where the GPU does not take the case; they were the CPU's unless `--gpu`.
  `--no-gpu` keeps the CPU, to the bit. On a Mac their results move by what single precision moves
  them: the demo approach's `--buildings` SELs by 0.001 dB at most, `niceops_point_field`'s band
  levels by 0.05. The GPU now samples the shared faces as the CPU does.
- **The demo approach's `--buildings` run** took 70 s at first; tabulating the transition function
  and sharing a shared face's nodes across wavenumbers took it to 37 s, sampling the shared faces
  to 10 s, sampling the states to 3.3 s, and the GPU to 1.7 s. Its levels moved only as the
  entries above say.
- **`--buildings` puts the microphones where the observer file does**, raised by its height column,
  as any run without a ground model. It had kept the ground model on, so an observer file without
  a height column put every microphone 1.2 m up (about 3 dB at street level against the observer's
  own point), and `--no-ground_effect` was refused with it. A database without spectra is now
  refused; it ran with every level `-inf` and exited cleanly. Observer files with a height column,
  as the tests' and the demo's, did not move.
- **A roof under a reflected leg is found whatever stands far off.** Its candidates came from the
  plan box of the leg's unfolded line, so whether a real bounce was found depended on the plan
  grid's cell size, set by every building in the district. **A bounce on the join of two level
  roofs is one ray**; both buildings took it, 1.4 dB loud. **Twice-diffracted rays are reciprocal
  at `--max_edge_reflections 0`.** **Walls against part of another's are split at the junction**:
  a side shared in part kept its facade and corners, and rays diffracted at the junction and ran
  along the seam between the buildings; no shipped or validation district has one. The demo
  approach's `--buildings` SELs did not move.
- **`--edge_source`: angles about an edge are clamped onto its air, and a partner's part reaching
  a shared vertex ends on it.** A district turned by 17-61° moved bands by up to 1 dB, and a
  partner stopping 1e-20 short of a vertex added its whole field to a ray, 10-16 dB.
  On EDGE's scale models the edge source's error against the FDTD fell from 2.73 to 1.47 dB rms
  for UE06's drone and from 0.92 to 0.23 for UE14's, where it had been 2.6 dB high at one band, and
  rose from 3.81 to 4.03 for UE07's floor source; against the measurements, read swapped, it went
  from 4.02 to 3.98 dB rms over all band levels. UTD did not move.
- **`--fade`**: a leg leaving an L-shaped building is tested against its other wing, and a ray
  across a building is not let through by keeping only the three nearest screens.
- **`niceops_point_field`'s roofs scatter** under `--scattering`, with their building's
  coefficients, and a building with any wall within `--diffuse_radius` scatters as a whole. A roof
  took its first wall's ρ(1 − s) and nothing returned the rest: microphones above the roofs rose
  by 1.2 dB on average and up to 1.9 dB (7.4 in a band), windows by up to 0.2 and the street by
  0.04 (FIDELITY). DigiCity's microphone 7 did not move: the loudspeaker stands below the roofs.
- **The GPU's functions are the precise ones**, on macOS 15 and later too, where Metal had
  defaulted to the fast ones whatever the math mode: the demo approach's GPU SELs went from within
  5e-4 dB of the CPU's to within 1e-6. A cone's sum below zero in a sub-band is taken as zero, as
  on the CPU.
- **`niceops_fdtd_field`'s transfer divides by the march's own time difference**, not by -iω, and
  the air's delay weighting goes on the source as on the pressure. Levels rose by 20 lg(x/sin x),
  x = ωΔt/2: up to 0.05 dB at the default 8 cells a wavelength, 0.001-0.045 dB in EDGE's runs;
  with `--humidity` the free-field test's 0.013 dB bias went to 0.002. Under `--mirror` the
  default duration covers the whole space, where it covered the half gridded (up to 1.55 dB at
  1 kHz in a run that had stopped early). A grid too large for the GPU marches on the CPU instead
  of failing.
- **The phase interval's quantile past a shape of 1e4** is Wilson and Hilferty's; the series had
  stopped short there and put the lower end on the level, within 0.001 dB.
- **The cone's phase interval takes a database's recorded spread and tonal fractions**, as the
  measured model's did; it took the defaults whatever the record said. No shipped database has a
  record, so nothing moved.
- **Buildings sharing a side are judged along it**, where one or both roofs are pitched: the
  shared wall stands, and the edge diffracts, only where one roof is above the other, cut wherever
  either roof breaks or the two cross. They had been compared by their highest points, so a gable
  against a taller flat roof kept a wall whose bottom was above its top, and the flat roof lost
  its edge on that side. A mansard's level deck now reflects and scatters, and a ridge is an edge
  only where its two planes are the roof on either side. The facade prefilter samples a wall
  between its bottom and its top; it had sampled up to half its top, so a wall raised more than
  half way by a neighbour never reflected, flat roofs too. UE06's drone against the FDTD went
  from 2.12 to 2.03 dB rms by UTD and from 1.47 to 1.20 by the edge source, and UE05's floor
  source against the measurements from −0.71 to −0.47 dB mean (3.09 to 3.14 rms).
- **The hybrid relief's grating rays reach the microphone straight from the relief**, or off the
  floor. They could reflect in up to two other walls first, which the diffuse field also
  re-emits, so that energy arrived twice. The source's image tree now reaches as deep as the
  grating rays need, so `--no_diffraction` no longer drops them; an order's share is the mean
  over the band's sub-bands; a floor bounce keeps its phase, and a plate its reflected term.
  Beyond a microphone's diffuse reach a Lambert relief reflects whole, and a relief building's
  roof no longer scatters. On EDGE's reliefs the hybrid went from 1.39, 1.76 and 1.98 dB rms
  (UE24-26) to 1.25, 1.70 and 1.80, and DigiCity's microphone 7 from −1.97 / 2.72 to −1.91 / 2.64.
- **A pitched roof's slopes reflect** (`--roofs`), as mirrors in their planes, with the floor at
  either end; and a source or microphone may lie on a slope, taking its reflection there (a plate
  on a slope is refused, as on a wall). Once-diffracted rays take a slope's reflection on either
  side of the edge, and slopes scatter, with tilted diffuse patches. They reflected and scattered
  nothing, so the field stepped across each slope edge's reflection boundary and across every
  boundary another edge cast on a slope's reflection (8 dB in the test), a scattering slope's
  rays lost its share s, and a point on a slope was refused. EDGE's pitched scenes do not move:
  every microphone there is low in a house's shadow, which no slope reflection reaches. In a
  suburb of 64 gable houses heard from a drone 60 m up, 197 of 810 garden band levels moved by
  more than 0.1 dB, from −3.7 to +1.8 dB.
- **`--edge_source`'s twice-diffracted rays decide each wavenumber alone**, and its end ratios let
  the other edge's point slide. The choice between the cascade-times-ratios and the integral was
  made at a run's least wavenumber, so a 1-2 kHz band moved 0.7-1.8 dB with whether 125 Hz was in
  the run; the ratios held the other edge's stationary point fixed, 1-3 dB off on coupled edges
  with one end near; and a ray with both edges' ends near, or near a partner's boundary, now takes
  the integral. DigiCity moved by up to 0.09 dB in a band (now within 0.001 of the slow integral)
  and the demo district by up to 0.7 (now 0.06); DigiCity's run takes 17 s where it took 10.5.
- **`niceops_urban_table` takes a microphone flush on a wall**: the wall's reflection at the
  microphone, `1 + sqrt(rho)`, on every ray and the first-order reflections the run would have
  evaluated, which the table now holds for it. One on or over a roof is refused by name. Before,
  such a microphone was refused as inside the building or took no doubling, by which side of the
  wall's line it rounded to.
- **`niceops_point_field --flow_resistance`'s diffuse floor** takes each band's absorption at its
  middle sub-band, as the run does; it took it at the band's centre (a few thousandths of a dB).
- **`--export_sphere` carries `rotor_rpm`**, as it carried the coherence record: a run on an
  exported sphere took the rotor's rate from the vehicle file instead, which also moved the band
  up to which the measured model's rotor delay holds.
- **Two edges across a step take the uniform form**: a ridge or a taller roof's far edge above a
  step down to a lower roof's edge, where the faces below the line between the edges meet beneath
  it. They took the cascade of two coefficients, which against boundary elements was 4.0-5.0 dB
  low from the drone 1.25 m behind EDGE's pent-roofed terrace (UE02) and 0.96 dB rms, 0.58 dB low,
  over a full-scale step; now 0.15-0.53 dB there and 0.71 rms, unbiased. EDGE's UE02 drone row,
  measured minus predicted, went from +3.89 dB (4.30 rms) to +1.93 (3.33), UE01's from +1.90
  (3.53) to +1.41 (3.34), and the pitched scenes' 554 band levels from +0.72 (3.08) to +0.47
  (2.96); the flat-roofed scenes did not move. DigiCity's microphone 7 moved by under 0.1 dB.
  Under `--edge_source` such a pair takes its transition from the uniform form too: the
  integral, which it took before, was 1-3 dB low there.
- **UTD takes an edge's ends** (THEORY 6.5, "An edge's ends"). Its edges were infinite: a
  diffracted ray whose Keller point left the edge at its end was dropped, so the field stepped where
  a building's corner cut an edge off (3.8 dB in the test's corner), and inside the edge the end's
  wave was missing (up to 7 dB against the edge-source integral past a box's corner, now 0.86).
  Each wedge term now loses each end's Fresnel tail, and a ray past an end is kept within 30
  radians of the run's lowest band. With it, an end that is not one stopped being one: an edge cut
  in two at a T-junction's seam goes on as the next piece, and an eave running square into a
  taller neighbour's wall that covers it goes on as its image (`Edge::beyond`); a low roof against
  such a wall had missed image theory by 1.3 dB as an end, and by 1.4 under `--edge_source`, whose
  floor-only continuation had the same false end. The goldens moved: `urban_buildings` and
  `urban_buildings_plate` by up to 0.041 dB, `urban_buildings_diffuse` 0.061,
  `buildings_scattering` 0.014. On EDGE, UTD against the FDTD went from 1.95 to 1.61 dB rms over
  the 13 runs (UE14's floor row 1.22 to 0.84, its drone row 0.68 to 0.45, UE06's floor row 2.77 to
  1.33 with R1 from +4.7 to −1.1 dB mean), the edge source from 1.35 to 1.32; against the
  measurements, UTD from 3.82 to 3.74 dB rms, the edge source 3.73 to 3.72.

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
- **Databases are compared by canonical path.** `-d databases/X.nod -d ./databases/X.nod` used
  to count as two databases, so `--provenance_file` and `--screen_database` refused the run; they
  now accept it, and the file loads once. A run that returns to a database after another one
  used to reload it, and now holds it until its last trajectory. Results are unchanged.

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
source, without a hover correction. Against each aircraft's approaches (NICEOPS minus the
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
Nothing else changed: the recipe, panam `eeb040b` and NICEOPS `e13577e`'s ray tool are the
2026-10-04 build's, and every straight-ray sphere came out bit-identical. Before the fix RO-44
scored +0.12 / 1.27 and RO-66 +0.86 / 1.57. RO-44's gain is below 20 degrees (0–5 degrees, rms
2.36 → 1.94; 5–10 degrees, 1.41 → 1.23); RO-66 is unchanged. NICEOPS predictions were never
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
