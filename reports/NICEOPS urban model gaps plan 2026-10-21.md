# Plan: the urban model's remaining gaps (2026-10-21)

The review rounds of 2026-10-20 and 2026-10-21 closed every bug they found (commits 29ff0aa to 776df25 on
`gpu-metal`). What they left are five gaps in what the point-source model and the plane-wave tables represent. This
plan takes each in turn: what is wrong and how we know, the approach, how it will be checked, what it costs, and
when it is done. The [next-steps report](NICEOPS%20Urban%20Propagation%20Next%20Steps.md) links here.

## Order

| Step | Gap | Why here |
| --- | --- | --- |
| 0 | EDGE's pitched drone rows read low (diagnosis) | Cheap, and its answer sets how much step 2 matters. **Done** (below) |
| 1a | Two edges across a step take the cascade | Found by step 0: most of UE02's bias. **Done** (below) |
| 1 | No vertex (corner) term in UTD | Affects every district, flat or pitched; the most visible artifact. **Done** (below) |
| 1b | Three or more edges in each other's transitions (a street at grazing, a roof's near edge with a step) | Found by step 0; the largest error left in the BEM checks. **Done** for rows with a street (below) |
| 2 | Missing slope and roof paths | Needs step 1's end treatment for its own edges' ends. **Done**, the rest measured and closed (below) |
| 3 | Roofs and slopes do not light the diffuse field | Reuses step 2's roof and slope images. **Done** |
| 4 | The plane-wave tables have no roofs | Independent of 1-3; can run alongside. **Done**, two leftovers deferred |
| 5 | Re-score and document | After each step, and once at the end. **Done** through step 2's close |

Each step lands as its own commit (or a few), with both suites passing, FIDELITY re-measured where it moves, and
HISTORY and the internals updated, as this branch's work has.

## Next steps (updated 2026-10-22)

Steps 0-5 are done or closed on measurement. What the work since has established, and what it points
at, in order. Each item says what was measured, what would be done, how it would be checked and what
it costs. Nothing here has been started.

**Progress (2026-10-22, later).** Item 1 measured and closed: `--fade` takes the suburb's sweep
from 10 steps to 9 (largest 3.04 → 2.23 dB) and DigiCity's 20 lines from 23 to 21 (median 1.6 →
1.8, largest 5.1 either way): seven of DigiCity's 21 are made partners whose own legs graze another
edge, so the handover moves rather than goes. Not worth a default, or an extension to twice rays,
on this evidence. Item 2 re-diagnosed: no partner was missing. A microphone exactly in a roof's
plane, or on the lower roof's reflection boundary, took both sides' waves; giving each boundary to
the ray set's side took the step's 72 cases above its lower roof from 2.34 to 0.38 dB rms (FIDELITY,
"Across a step"). What is left of it is the switch's own step (2.5-4.3 dB across a roof's plane with
the source above it, 0.9 across the reflection boundary), which a blend of aperture and cascade
would remove, the GPU's aperture kernel included.

**Item 2 built for flat roofs (2026-10-22).** The plane's step was not the switch alone: the
cascade itself steps there by up to 12 dB, the far edge sitting in the near edge's transition with the
source low over the roof, so a blend toward it would have kept the wrong level above the plane. With
both ends above a flat roof's plane the pair now takes each end's mirror below it, negated, the
end nearer the plane governing (THEORY "Both ends above the face's plane"). Against 528 new BEM
cases through the plane behind the step: 0.88 -> 0.20 dB rms, worst -4.3 -> -1.0; the step at the
plane 8.4 -> 0.45 dB at most; far above the plane (336 cases, sources to 40 degrees) it and the cascade
agree. The GPU takes it unchanged (the sign rides in the mirrored aperture's node weights), at no
measurable cost. Walls between corners take it too: a quarter of DigiCity's field 1.2 m up moves by
over 0.1 dB with a source 12 m up (to 4.4 dB), its measured microphones by 0.04 dB at most; EDGE's
scores do not move (UE20 and UE27 by 0.38 dB at most in a band). Left: a
step's tilted pair with the source just above its line steps 2.2-2.7 dB where the microphone crosses
the line; there the direct ray and both single diffractions are at their boundaries at once, the
mirror is 12 dB off, and neither the cascade nor the mirror fixes it. Also left: chains of three and
four edges over a row of roofs still step 3.9 dB in the roofs' plane where no profile is taken; the
reflection boundary's 0.94 dB is unchanged (where it meets the plane, 10 m behind, 1.05-4.16 ->
0.86-1.16 dB).

**Item 3 closed (2026-10-22, not pursued).** Built (round trips across each street, faded out above
k w = 120, the walls' reflection per band) and measured against new BEM: it held at single frequencies
for streets 0.5-3 times as deep as wide (1.19 -> 0.58, 0.90 -> 0.59, 1.39 -> 0.67 dB rms). But BEM at 16
frequencies a band showed the march alone already within 0.22 / 0.40 / 0.66 dB rms of the band levels
(one street, two, a 20 m street), at 125 Hz -0.31 to +0.18: the "bias" was where single frequencies
fall among the street cavity's resonances. Band-averaged the trips gained 0.15 dB at most, and at a
run's three samples a band their fine structure made the levels worse (0.49 / 0.72 / 0.93 against
0.38 / 0.54 / 0.82). Reverted, recorded in FIDELITY ("Over a row of roofs") and internals; patch in job
tmp/item3. The user chose not to pursue it further for now.

**Item 3 measured in two dimensions (2026-10-22, not built).** A two-way pass over each street --
the field below the next roof reflected off its wall and marched back, below this roof reflected
forward again, the part over the next roof added -- is what the bias lacked. Against the BEM, rms:
one street 1.61 -> 0.58 dB at 24 round trips, two streets 2.10 -> 0.69, the step 1.18 -> 0.41, the
steep source 1.03 -> 0.29; the 125 Hz bias -2.3 -> -0.9. But a street is a weakly damped cavity at
high frequency: at 1 kHz the partial sums swing 0.25-0.55 dB rms with the trip count (forward only
0.28), and each trip costs two marches over the full column. Trips only at low frequency (16, the
last half tapered, faded out between a street 60 and 120 radians wide) keep most of the gain (one
street 0.72, two 0.87, step 0.64, steep 0.46, 1 kHz unchanged) at a cost still to measure.
Prototype: job tmp/twoway/tw.cpp. Decision for the user: build that, and whether a frequency fade
is acceptable.

**1. Partners beyond the orders, by `--fade` (first).** The largest defect left in the field.
- *Evidence.* Both continuity sweeps put most of their steps in one class: a ray whose partner at a
  boundary is beyond the traced orders. DigiCity (80 lines, 79,460 positions): 355 of the 434 steps
  above 1 dB, a diffracted ray's leg meeting another edge or its reflection leaving its wall (265),
  and a specular ray's leg meeting an edge whose partner needs three or four reflections (90); that
  sweep predates three roof edges, reflections between edges and two reflections at each end, which
  took 20 of its lines from 31 steps to 23, so these counts are high. The
  gable suburb (8 lines, 882 positions): all 10, once-diffracted rays with three or four wall
  reflections and twice-diffracted rays with walls. Raising `--max_edge_reflections` to 4 makes it
  worse (490 steps, 10.1 dB at most), and UE06's R1 (+4.7 dB, a reflection 3 cm inside a wall's end
  kept whole) is the same failure. The research review's "truncate together" (Calamia 2009) is the
  remedy.
- *What exists.* `--fade` makes that partner from the ray it partners, for specular and
  once-diffracted rays, in `niceops_point_field` only, CPU only, off by default. Its effect on the
  sweeps has not been recorded.
- *Approach.* Measure `--fade` on both sweeps first (steps removed, cost). Then extend it to
  twice-diffracted rays (the suburb's main case), bring it to the run (`niceops --buildings`) and the
  GPU, and decide whether it becomes the default.
- *Check.* The sweeps' step counts; the BEM and image-theory tests unchanged; UE06 R1 and EDGE
  against the FDTD; cost on the demo approach and DigiCity's flyover.
- *Cost.* Moderate; the GPU path is most of it. The sweeps take 20 minutes (suburb) to an hour or two
  (DigiCity) each.

**2. A roof between two edges.** The step tables' remaining error.
- *Evidence.* FIDELITY, "Across a step": microphones above the lower roof are 2.55 dB rms against the
  boundary elements (10.8 dB at worst), with either form, nearly all of it at the lower roof's
  reflection boundary for the first edge's wave, where the ray set steps by up to 6 dB. The partner
  that would take over is a twice-diffracted ray reflected in the lower roof between its edges; a
  reflection between two edges is a facade or the floor today. The row of roofs' profile march does
  not help touching blocks (a step up: 2.04 dB rms in two dimensions against the cascade's 1.99).
- *Approach.* A level roof or a slope as a reflection between the edges, the second edge mirrored in
  its plane, as the floor between two edges already is.
- *Check.* The stored step tables (`tests/test_point_paths.cpp`), and continuity across that
  boundary.
- *Cost.* Small to moderate.

**3. The row of roofs' low-frequency bias.** (Closed 2026-10-22: not a bias in band levels; see above.)
- *Evidence.* The profile field is low against the boundary elements by 2.2-2.4 dB at 125 Hz over one
  street, 2.2-3.6 over two, and about 1 dB at 250-500 Hz: the street's reverberation between its
  walls, which a one-way march cannot carry back. Taking the walls' rays back in as cascades
  overshot by 3 dB.
- *Approach.* Measure in two dimensions before building: the walls' rays at a share, or a reverberant
  term from the ideal canyon's solution (§6.6), against the stored BEM cases.
- *Cost.* Moderate. Risk: no form may hold across street widths.

**4. The row of roofs' cost in a run.** No shipped case exercises it: the demo goldens, DigiCity and
EDGE barely reach it. A row costs 0.29 s to 10 kHz for each microphone and state that sees one, and on
the GPU path the profile factors are still formed on the CPU. A terraced district with a low approach,
run with the profile and without, would say whether a run can afford it and what it moves. Needs a run
flag to switch it off (`--no_profile_pe` is point_field's only).
- *Measured (2026-10-22), by sampling, not a full run* (a full one ran over two hours on this laptop,
  which sleeps). `niceops --buildings` now takes `--no_profile_pe`; the shipped test district reaches no
  row (byte-identical, a CTest). A synthetic terraced district under the demo approach: 8 flat-roofed
  rows 30 ft deep and high, 600 ft long, 40 ft streets, 264 observers 4 ft up in the streets and
  behind. 22% of (observer, state) pairs take a profile. Per such pair, one core, 25 Hz-10 kHz: 3.18 s
  with the profile, 0.086 s with the cascades (37x) -- the march spans 8 rows, about 160 m, and
  often both ways -- so a run there costs about 9x. What it moves (30 pairs): A-weighted with a flat
  source spectrum, profile less cascades -0.20 dB mean, -0.01 median, -2.1 at worst; by band -0.1
  to -0.5 dB mean, -5.7 at worst (200 Hz). The march's cost goes as f^2: on one 8-row cross-section
  the bands to 1 kHz took 3.6% of it, to 2 kHz 8.2%, to 4 kHz 21.6%. Its column also reaches the
  source's height (11.8 s one way with the source 122 m up), where the start needs only the wave's
  angles. Decision for the user: which lever, if any (a frequency cap with a fade to the cascades;
  the column's height; profile off by default in the run).
- *Lever 1 built (2026-10-22, the user's choice).* The profile fades out by the row's size in
  wavelengths, k D from 3000 to 6000 (D the row across), not by frequency: a fixed 2-4 kHz cap was
  cheaper here (0.77 s a pair) but would strip EDGE's scale-model rows, which sit at kilohertz. On
  the same 30 pairs: 3.18 -> 0.98 s a pair (37x -> 12x the cascades, the run about 9x -> 3.3x), the
  A-weighted effect -0.20 -> -0.19 dB mean, worst -2.13 -> -1.79. The rays it stands for are kept
  for above the fade. Lever 2 (the column's height) is still open if a run needs more.

**5. Deferred, with the evidence for deferring.**
- *Step 1, part 2 (the double transition).* A strip under 0.5 m wide where a face's plane passes
  through a vertex; worst 3.68 dB at 125 Hz on the test box; the render does not show it. If it is
  ever wanted, hand those rays to the edge-source integral (less work) rather than derive the
  generalized-Fresnel coefficient.
- *Step 2's three remaining paths* (an end's slope before its walls on a diffracted ray, walls on a
  twice-diffracted ray's roof side, roofs on thrice and chain rays). No step from them in the
  suburb's sweep; the one measured kind of their family, a roof at both ends, moved levels by 0.08 dB
  at most.
- *Step 4's leftovers* (in the tables, a slope with facades between it and the microphone; the plane
  wave's lighting off a slope). They touch only tables over pitched roofs; measure on the suburb's
  table before building.
- *Step 1b for touching blocks and pitched steps.* No method found that beats the cascade there;
  item 2 addresses the measured main error over a lower roof.
- *Frequency sampling of the profile, and Numerov in height.* Built, measured and dropped (internals).
- *Data, not the model* (hypotheses, not findings): from the floor source EDGE's measurements read
  4-6 dB above the FDTD at UE12's and UE13's R3; UE17's R1 stays open; R3 and R4 are free-field
  capsules where R1 and R2 are pressure ones, which may explain their high-frequency trends.

**6. Re-score** after items 1-3, as step 5 did: EDGE, DigiCity's microphone 7 and accuracy figure,
the demo approach's goldens, the suburb with its sweep, and the 3-D render.

## Step 0: why EDGE's pitched drone rows read low

**What is wrong.** Measured minus predicted, from the drone source (S3, 1.5 m up, aimed 30° down): UE01 +1.9 dB
mean (3.5 rms), UE02 +3.9 (4.3), UE03 +0.7, UE04 +1.3, UE05 +0.3, UE15 +0.3, UE17 +1.1 (4.7 rms). The floor-source
rows of the same scenes are within ±0.5 dB on average. The model is low from above, not from the floor. The slope
reflections added on 2026-10-20 do not reach these microphones: every one is low in a house's shadow, and the drone is
below the plane of the slope that faces them.

**Approach.** Diagnosis, before any model change:

- Break each drone row down by receiver and band. A bias in every band points at the source or its calibration; one
  that grows with frequency points at geometry or missing paths; one at particular receivers points at paths.
- Run the FDTD on the pitched scenes in the bands where its staircase is fine. At 3 mm cells the staircase is 5 dB
  loud at 6.3 kHz; find the band below which it agrees with the boundary elements on the long gable, and use it as
  the reference there. FDTD minus model then separates the model's error from the data's.
- Check the drone loudspeaker's directivity off its axis (the receivers are 40-70° off it), and the per-band
  calibration taken from the LoD1 scenes' line-of-sight pairs, which never saw S3 at those angles.
- List the rays at the worst receivers, and the paths the model does not trace that a full-wave field would carry
  there: a slope's reflection off a wall, a ray diffracted twice with a roof or slope bounce.

**Done when** the bias is assigned to the data (source, calibration, geometry), to specific missing paths, or to
neither, with the numbers in FIDELITY. If it is missing paths, step 2 is ordered by which ones.

**Result (2026-10-21).** Assigned, and not to missing paths:

- **UE02's drone bias is the model's**: the twice-diffracted ray over the ridge and the deck's edge, which share no face
  (a 0.12 m step between them), took the cascade. The deck's edge lies 1.2° inside the ridge's shadow boundary. Against
  the exact 2-D field (boundary elements, converged to 0.1 dB) the cascade was 4.0-5.0 dB low at R2's distance, where
  the measurements read 4.3-4.9 high; UE01's step explains 0.2-1.6 dB of its +1.9. UE03-05 share faces and were
  already uniform. The fix is step 1a.
- **R3 and R4 are G.R.A.S. 40BF free-field capsules, R1 and R2 B&K 4136 pressure ones** (EDGE's equipment table).
  Above 8 kHz R4 reads up to 10-13 dB high and R3 5 dB low in every pitched scene, UE03 included, where the 2-D model
  is exact; the calibration is per loudspeaker, not per receiver. These may be the capsules' response off axis, or
  their positions: the data do not say which.
- **UE15 against the FDTD** (3 and 2 mm cells, extrapolated to zero cell for the stepped slopes): R3 and R4 within
  0.4 dB of the model, R1 and R2 within ±2.3; at R1 the FDTD is 0.45 dB above the model on average where the
  measurements are 2.25 above, so about 1.8 dB of R1's excess lies outside the full-wave field. R1 reads 1-2 dB high
  in most pitched scenes. It may be the receiver.
- **UE17's R1** (+7 dB, rising to +12 at 20 kHz, a 1/5-scale detached house) is too large for the FDTD at these
  cells, and stays open.
- **Found on the way: a street at grazing.** Two 10 m blocks 10 m apart, a source 1.1° above the roofs: the model is
  up to 8 dB high at 500 Hz, and its cascades over two, three and four edges do not converge (dropping the
  thrice-diffracted rays moves the mean from +1.1 to −3.8 dB). That is step 1b.

Slope diffraction (Holm's higher-order series, rederived) was tried for the transitions and diverged from the second
order; the Kirchhoff aperture, generalized, is what worked.

## Step 1a: two edges across a step

**Done (2026-10-21).** The shared face's aperture (`double_wedge`) generalized to a pair whose faces below the line
between the edges are tilted off it and meet beneath it, a closed step, with the line touching both wedges from above.
Against boundary elements: UE02's drone row at 1.25 m from 4.0-5.0 dB low to within 0.53, UE01 from rms 0.66 to 0.23,
a full-scale step's 229 cases from rms 0.96 (0.58 low) to 0.71 (unbiased). What it leaves: the level steps 0.5-2 dB
across the line through the two edges (the cascade 0-1.5), worst for a source grazing the taller roof, which the
three edges over that roof would put right (step 1b). FIDELITY has the tables ("A point source among buildings").

## Step 1b: three or more edges in transition

**What is wrong.** Where three or more roof edges lie in each other's transition regions -- a row of roofs at
grazing, a step whose taller roof the source grazes -- the cascades of two, three and four edges each carry a share
of the field and do not converge. The street at grazing is 8 dB high at 500 Hz.

**Approach.** Measure first, as in step 1's part 2. Candidates: the aperture applied in turn at each intermediate edge
(a multiple Kirchhoff integral over the taut string's corners, as Vogler's multiple knife-edge integral does for
screens; costly per pair, but only where the chain is in transition); Tzaras and Saunders's heuristic for the
transition functions' distance parameters; and the edge source's higher orders. The boundary-element tool now takes
several blocks (`--blocks`), so streets and rows can be solved exactly.

**Done (2026-10-22) for rows with a street.** The row's field is marched, a wide-angle parabolic equation over its
cross-section from the first block's front face, finished through the last block's back corner by that wedge's
UTD Green's function (THEORY "Over a row of roofs"), prototyped in two dimensions against the boundary elements
before it went into the point model. In the point model it replaces every ray over the row with
four, the parts from the source or its floor image to the microphone or its image, each its own length, weighed in
by elevation (40-25°), by the line's clearance of the blocks' ends, and by each end's angle off the vertical over
its outer corner. Against boundary elements, rms: one street 1.65 dB (cascades 3.40), at 1 kHz 0.33 (3.49), a
source 16.7° up 1.13 (2.19), a step across a street 1.15 (1.86), two streets 2.15 (6.70). What is left is a low
bias at 125 Hz, 2.3-3.6 dB, the streets' reverberation the march cannot carry back (at single
frequencies; averaged over a band it is not there, next steps, item 3). Touching blocks (no street)
keep the cascades, which hold as well there; so does UE02's pitched step. EDGE: only UE12's and UE13's receiver R3,
behind both terraces, takes it, and moves toward the FDTD (UE12's floor source from 3.5 dB above it to 1.0 below;
the row's rms 2.30 → 1.65 dB), away from the measurements, which read 4-6 dB above the FDTD there (harness
RESULTS). Cost: 0.29 s a row to 10 kHz on one core after the optimization that followed (split-step Padé),
from 6.7; that work also found the high-frequency domain unconverged (up to 3 dB at 8 kHz) and fixed it.

## Step 1: a vertex term for UTD

**What is wrong.** UTD diffracts at each edge as if it were infinite. A diffracted ray whose Keller point leaves the
edge at its end is dropped, so the field steps where a building's corner cuts the edge off. Two roof edges meeting at
a corner, each with its Keller point within centimetres of it, give the one-frame, 4 dB dips in the DigiCity 3-D
render. About 1% of the render's surface cells sit 3 dB or more off their neighbours in a frame, most beside corners.
The edge-source integral handles ends exactly, but it is CPU only and slow: DigiCity's four microphones take 17 s
with it, against UTD's 1.3.

**The reference.** For rigid wedges, the edge-source integral (Svensson, Fred and Vanderkooy 1999) summed over a
corner's edges is the corner's first-order diffraction exactly: each edge's integral stops at the corner, and the end
contributions are the vertex term. So the asymptotic form of those end contributions *is* a vertex coefficient, and we
already have the exact integral to check it against. The UTD vertex coefficient of Albani, Capolino, Carluccio and
Maci (IEEE Trans. Antennas Propag. 57(12), 3911-3925, 2009) does the same for perfectly conducting faceted
structures, with transition functions built on generalized Fresnel integrals. Its form near the double transition
(the Keller point at the end while the observer is near the edge's own shadow or reflection boundary) is the one to
follow there; it is electromagnetic, so it needs re-deriving for rigid (Neumann) acoustic faces.

**Approach**, in two parts:

1. **Ends as a Fresnel transition.** Keep a UTD ray whose Keller point lies past an end, within a reach (as the edge
   source already does, `edge_source_reach`), and multiply every once-diffracted UTD ray by each edge's end factor: the
   asymptotic form of the end ratio `end_ratios` computes from integrals, a complex Fresnel integral of the end's
   signed path excess, which is 1 far inside the edge, 1/2 at the end and falls off past it. The excess comes from the
   ray's own geometry, so this costs a few flops a ray. On a twice-diffracted ray each edge takes its own, with the
   coupled curvature `end_ratios` uses (A - C²/B). This removes the step wherever one wedge's coefficient is smooth
   along the edge, which is everywhere but the double transition.
2. **The double transition.** Where the end and the edge's shadow or reflection boundary meet, the product above is
   not uniform. Measure how far off it is against the integral first. If it matters, adopt the generalized Fresnel
   transition of Albani et al., re-derived for the rigid wedge; if only marginally, hand those rays to the integral,
   as `--edge_source` already selects by transition parameter.

The GPU: `GpuWedge` gains the two ends' excess per unit wavenumber, and the kernel evaluates the Fresnel factor per
problem (a rational approximation, as `faddeeva` is). The plane-wave tables take the same factor at their edges' ends.

**Validation.**

- Against the edge-source integral, the exact reference for finite rigid edges: sweeps over single corners, roof
  corners where three edges meet, and DigiCity's containers, per band. Target: within 0.5 dB everywhere outside the
  double transition, and the double transition reported separately.
- Continuity: a microphone moved across an edge's end-boundary in micrometre steps; the step falls from the
  diffracted ray's whole level to under 0.05 dB.
- FDTD on EDGE's small-house scene UE14, where ends matter: against the FDTD its drone row is 0.68 dB rms by UTD and
  0.22 by the edge source, the floor row 1.22 and 0.74. The end factor should close most of that gap.
- The DigiCity 3-D render re-run: the share of cells 3 dB off their neighbours should fall from about 1%.

**Cost and risk.** Moderate. The factor is closed form; the work is its derivation for each ray kind (once, twice,
thrice and chains, the shared faces, the edge source's own fade) and the GPU kernel. The risk is the double transition,
which is why it is measured before it is built.

**Result (2026-10-21).** Part 1 is in (THEORY 6.5, "An edge's ends"; FIDELITY, "An edge's ends in UTD"):

- Each of D's four terms loses each end's Fresnel tail, weighted by 1/sqrt(1 + k Delta / X_i) so that the term
  carrying a ray's own boundary keeps its jump (one factor on D broke those boundaries by 0.29 dB). A ray past an
  end is kept within 30 radians of the run's lowest band, through the end, with the Keller point's coefficient; past
  an end a source or microphone beyond a face's plane takes that face's coefficient weighed by its depth. Twice and
  chain rays take each end with the other edges' points free on their lines; shared faces take the ratio; the GPU
  and the plane-wave tables take the same.
- Found on the way: **where an edge goes on.** A T-junction's seam cut edges into collinear pieces, each with false
  ends, and an eave running into a taller wall is no end where the wall covers its wedge (image theory), but is one
  on a flush terrace, where the wall stands only behind the eave's facade. Taking those right closed the image-theory
  test to 0.02 dB (it was 1.3 by UTD and 1.4 by the edge source).
- Against the edge source past a box's corner: 0.18 dB rms, 0.86 at worst (1.4 and 7.1 with the edges unending).
  Continuity at an end and across a face's plane past it: 1e-5 dB per micrometre.
- EDGE against the FDTD: UTD 1.95 -> 1.61 dB rms over the 13 runs; UE14 floor 1.22 -> 0.84, drone 0.68 -> 0.45 (the
  plan's check); UE06 floor 2.77 -> 1.33. Against the measurements UTD 3.82 -> 3.74.
- **Part 2, the double transition, measured:** within 0.5 m of a face's plane through the vertex the ends are no
  better than none (worst 3.68 dB against 3.13 at 125 Hz on the box). It is a strip a fraction of a metre wide, and
  the render below no longer shows it; left as a known limitation (FIDELITY) until a measurement says it matters.
- **The DigiCity render** (177 frames, 11,250 surface cells), against its run of 2026-10-20: cells 3 dB or more off
  their neighbours' median 0.33% -> 0.17% of cell-frames; one-frame dips of 3 dB 0.025% -> 0.010%, and the roof
  cells ever showing one 14.4% -> 0.17% (the corner dips). The plan's "about 1%" was a different count, which could
  not be reconstructed; these are the harness's own (job script, both fields scored the same way). Microphone 7 is
  unchanged (-4.1 dB mean, 5.3 rms).
- Cost: DigiCity's 22 microphones from 14 source positions, 17 bands, one thread: 6.0 -> 9.0 s, 30% more rays.

## Step 2: the slope and roof paths that are missing

**Slope–wall milestone (2026-10-21).** Specular wall chains on either side of one slope, and
single-edge paths with end -> walls -> slope -> edge on either side, are implemented. The
source's specular tree is independent of the diffraction switch; floor partners fold into real
coordinates. Oblique image theory, reciprocal paths, floor-contact launch, pruning agreement
and a screened reflection boundary are checked; measured results belong in FIDELITY.
Committed ee2644e.

**Twice-diffracted roof milestone (2026-10-21).** A twice-diffracted ray takes a level roof or a
slope (after the floor or not) between one bare end and its edge (`twice_roofs`,
`--no_twice_roofs`). Where a second edge cuts off a once-diffracted roof ray the field now steps
by 1-2% of that ray (UTD's own error at the second edge; the same without a roof) instead of
100%. Suburb: 236 of 810 band levels > 0.1 dB, max 2.4 dB, no cost; DigiCity flyover: 1.6% of
band levels > 0.1 dB, max 3.1 dB, +6% time. FIDELITY has the numbers.

Still open in this step: end -> slope -> walls -> edge; walls on a bounce side of a twice ray;
bounces at both ends; roofs on thrice and chain rays. Step 3's roof-reflected diffuse lighting
remains separate.

**Measured, and closed (2026-10-22).** Bounces at both ends are in: the source's roof bounces
against the microphone's. They move the suburb's rooftop microphones by 0.078 dB at most, for 5%
more time, and the street microphones not at all. To see whether the other three matter, a
continuity sweep of the suburb (8 lines 2.5 cm apart, 882 positions, the drone and a street source)
found 10 steps above 1 dB: two from once-diffracted rays with three or four wall reflections and
eight on one shadowed line from twice-diffracted rays with walls, the partners DigiCity's sweep
also found missing (more reflections, or a third edge). None is a ray off a roof or a slope. The
three paths left are energy the model does not carry, not a step in its field, and on this
evidence small; they stay out until a measurement says otherwise.

Found on the way, not step 2's, and fixed: with a microphone exactly in a slope's extended plane,
the shared face's aperture had one end on the line between the edges and lost its precision
(GCC against clang up to 10 dB; GPU against CPU 3.4 dB). It is now taken at its limit; the
compilers agree 1 µm off the plane and the GPU is within 0.04 dB. Left: at the exact point the
tracer admits grazing rays that exist nowhere else (up to 6 dB off the field 1 µm away).

**Original scope.** Three kinds of path were missing:

- a slope's reflection together with a facade's (wall then slope, or slope then wall);
- a twice- or thrice-diffracted ray with a roof or slope bounce on a side (level roofs have the same gap: "Not done:
  twice-diffracted rays off a roof" in the internals);
- a once-diffracted ray with a slope on a side that also has facades.

None of these makes the field step, because none is ever traced; each is energy that is simply missing. They matter
most in streets and courtyards of pitched houses, and possibly behind EDGE's pitched rows (step 0).

**Approach.** The tilted mirror does not commute with the walls' vertical ones, so the level roofs' shortcut (mirror
in z anywhere along the facade chain) does not carry over. It does not need to: put the slope *between* the two ends'
facade chains. The source's image through its chain A and the microphone's through its chain B see the slope
directly. The reflection is then the tilted mirror between those two images, checked by tracing each chain back from
the reflection point (`folded_path_open`). That is the structure the once-diffracted roof sides already use. For
diffracted rays the slope goes on the side of an edge between that end's chain and the edge, as step 1's partners do,
now with facades on that side. The twice-diffracted search gains roof and slope sides at either end the same way. The
images are already traced; the cost is the pairs, which a plan-view test prunes. A slope's building must lie within
both images' beams, the windows `image_tree` already carries.

**Validation.**

- Image-theory exactness: a wall at right angles to a long slope reflects as the two mirrors' composite image, to
  0.05 dB.
- Continuity: a microphone across each new path's boundaries (a wall's end, the slope's edge), stepped.
- The pruning test (`pruning finds every ray among pitched houses`) extended to the new kinds: the pruned search
  finds what the unpruned one does.
- EDGE's pitched scenes, the gable suburb, and FIDELITY's tables re-scored.

**Cost and risk.** Moderate to large, mostly in the pair search's cost among many pitched houses. The suburb's 110 s
is the budget to watch.

## Step 3: roofs and slopes light the diffuse field

**Done (2026-10-21).** A point source's lighting paths take one level roof (anywhere along the
patch's chain) or slope (next to the source), carrying its rho (1 - s) and counted among the path's
reflections (`DiffuseOptions::roof_lighting`, `--no_roof_lighting`). Image theory holds to 1e-9 for
a wall over a flat roof and over a slope. Suburb, s = 0.3: +0.04 dB mean, +0.16 at most, +1.6% time.
The plane wave's lighting (`irradiance`) waits for step 4's roof traces.

**What is wrong.** A wall or roof patch is lit only by paths through facade reflections and the floor
(`point_irradiance`, `irradiance`). A wave that reflects in a roof or a slope on its way to a patch lights nothing, so
a drone above a district under-lights the walls its roofs reflect onto. That applies to level roofs as much as to
slopes.

**Approach.** Give each patch's lighting paths one roof or slope reflection, as the rays have: the patch's images in
the roof's level plane or the slope's tilted one, with their facade chains as in step 2, and the same plan-view
pruning. Each such path carries the roof's or slope's ρ(1 - s), as a facade's does.

**Validation.** Lambert's law for a wall lit only off a roof (analytic: a large lit wall beside a large roof, the
source above the roof); the forward solve against the reverse; and energy bookkeeping on a closed courtyard, where
the total diffuse energy should rise by what the roof reflections carry.

**Cost and risk.** Small to moderate. The risk is the lighting setup's time on large districts, which
`prepare_point_lighting` already caches per patch.

## Step 4: roofs in the plane-wave tables

**Done (2026-10-21).** The tables reflect in the level roofs and slopes below a microphone
(`trace_off_roof`; a slope next to the microphone; diffracted rays to its mirror in a roof), a
microphone on a roof or a slope takes its reflection there, `--roofs` and `--no_roofs` in the
builder, and the roofs scatter and light the plane wave's diffuse field off a level roof. Against a
point source 1e8 m off: 0.002 dB over a flat roof and a slope, on and over each; flush against
0.1 mm above, 0.01 dB on rigid roofs. Left: a slope with facades between it and the microphone,
and the plane wave's lighting off a slope.

**What is wrong.** `niceops_urban_table` traces walls and the floor only, and it refuses a receiver on or over a
roof. It also takes no `--roofs`. A rooftop or window-level microphone therefore needs `niceops --buildings`, which is
slower for long tracks.

**Approach.** In `BuildingReceiver`, let the specular traces reflect in level roofs and slopes (`first_hit` meeting
a roof from above, the trace turning there), and the diffracted images include the roofs' mirrors, as the point model
does. Then accept receivers on and over roofs (a roof's doubling through `surface_factor`, as the walls have), and
`--roofs` in the builder. Scattering roofs come with it (`DiffuseOptions::roofs`).

**Validation.** The table and the point-source model meet when the source is far away: a rooftop receiver's table
entry against `niceops_point_field` with a source a kilometre off, both ways round, to the band-averaging noise. A
receiver flush on a roof against one a hair above it. The canyon and long-building tests unchanged.

**Cost and risk.** Moderate. The table builder's speed on large districts is the thing to keep.

## Step 5: re-score and document

**Done for steps 1-4 (2026-10-21).** EDGE: only UE07's drone rows move after step 2, by 0.04-0.08
dB rms the wrong way, from the twice-diffracted roof rays (harness f7e3fe7); steps 3 and 4 do not
reach it (no scattering, no tables). DigiCity microphone 7 and the accuracy figure: unchanged (the
loudspeaker stands below the container roofs; harness 1ae0bbc). The demo approach's `--buildings`
SELs: within 1.2e-6 dB (the shared face's limit); the table runs' goldens 0.001 dB (step 4's
prefilter) and 0.005 with scattering. The gable suburb: step 2 up to 2.4 dB at 125 Hz in near-nulls,
step 3 +0.04 dB mean with scattering. The 3-D render: 0.85% of cell-frames move by more than 0.1 dB,
at most 1.7 dB, its continuity shares unchanged. FIDELITY and HISTORY carry each. Step 1b, for rows
with a street, followed on 2026-10-22 (above).

**Done for step 1b (2026-10-22).** EDGE: only UE12's and UE13's receiver R3 moves, toward the FDTD
(harness 0fa41a3, 599817f). DigiCity microphone 7 and the accuracy figure: unchanged. The gable
suburb: unchanged by construction, its roofs all pitched. The demo approach's goldens: unchanged.
The 3-D render: no cell-frame moves, at the same cost (27,000 s of CPU).

After each step: EDGE (FDTD and measured, flat and pitched), DigiCity microphone 7 and the accuracy figure, the demo
approach's `--buildings` SELs, the gable suburb, and the 3-D render. FIDELITY's tables are re-measured where they
move, and HISTORY says what moved and why. The next-steps report is updated as each gap closes.

## Decisions for you

Settled: step 0's FDTD time was spent; step 1a went before step 1; step 1b took the profile march
(measured against the boundary elements before it was built); step 2 closed on measurement.

Open:

- **Whether `--fade` becomes the default** (next step 1), once the sweeps say what it removes and what
  it costs. It changes every run's levels near those boundaries.
- ~~The low-frequency bias's treatment~~ (next step 3): closed; there is no bias in band levels.
- **Step 1, part 2**: deferring is the recommendation; if wanted, the edge-source integral for those
  rays rather than a new coefficient.
