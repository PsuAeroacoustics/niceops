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
| 1 | No vertex (corner) term in UTD | Affects every district, flat or pitched; the most visible artifact |
| 1b | Three or more edges in each other's transitions (a street at grazing, a roof's near edge with a step) | Found by step 0; the largest error left in the BEM checks |
| 2 | Missing slope and roof paths | Needs step 1's end treatment for its own edges' ends |
| 3 | Roofs and slopes do not light the diffuse field | Reuses step 2's roof and slope images |
| 4 | The plane-wave tables have no roofs | Independent of 1-3; can run alongside |
| 5 | Re-score and document | After each step, and once at the end |

Each step lands as its own commit (or a few), with both suites passing, FIDELITY re-measured where it moves, and
HISTORY and the internals updated, as this branch's work has.

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

## Step 2: the slope and roof paths that are missing

**What is wrong.** Three kinds of path are not traced at all:

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

After each step: EDGE (FDTD and measured, flat and pitched), DigiCity microphone 7 and the accuracy figure, the demo
approach's `--buildings` SELs, the gable suburb, and the 3-D render. FIDELITY's tables are re-measured where they
move, and HISTORY says what moved and why. The next-steps report is updated as each gap closes.

## Decisions for you

- **Step 1, part 2**: if the double transition matters, should it take the generalized-Fresnel coefficient (more
  work, closed form, GPU-able) or hand those rays to the edge-source integral (less work, CPU only, slower)? I would
  measure first and recommend then.
- **Step 0's FDTD runs** take GPU time (minutes to tens of minutes a scene at 3 mm). Fine to spend?
- **Order**: steps 1 and 4 are independent; if rooftop receivers in tables are wanted soon, step 4 can go first.
- **Step 1a was taken before step 1**, as step 0 found it behind the largest measured bias; step 1 (the vertex term)
  is next unless step 1b is wanted first.
- **Step 1b**: which multi-edge form, measured on streets and rows against boundary elements before any is built.
