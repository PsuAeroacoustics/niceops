# Grazing incidence: what curved rays, rebuilt spheres and the measured hovers showed

Refraction explains most of the grazing-incidence propagation error, and what is left is not modeled propagation. With each run's own profile, curved rays take the 0–2° board/pole transfer-function error from 7.6 to 0.8–1.2 dB rms with the ground model unchanged, and the parabolic equation (PE) agrees with the ray models at every geometry checked. Two faults in `--refraction stratified` were found and fixed on the way: rays bounced along a surface duct were never searched, and a low source in a weak upwind gradient got no shadow boundary at all. The B407 spheres were rebuilt three times. The shipped database now depropagates along curved rays, refuses to divide out a ground board's near-grazing null, and grids with a floored, azimuth-stretched Shepard weight. Against the day-286 approach SEL it improves rms from 1.42 to 1.37 dB under `stratified` (1.31 to 1.26 under the station atmosphere). The residuals that remain sit near the horizon and are not explained by propagation as modeled, nor by the boards, the gridding or the database blend: low sources at short range read 3–4 dB loud, fast high flight 1.5–2.8 dB quiet at long range. The 1.2 m poles carry the same error as the boards beside them, so the cause is shared: the source's levels near the horizon, or something along the path the models leave out.

This report follows the plan in `Grazing incidence soft ground propagation.md` (30 September). Its first two ranked steps, the full-wave check and the curved-ray pair rerun, are done; the profile ensemble and the reciprocal PE table were not needed for what was found. All numbers are from runs made on 1–2 October 2026 on branch `stratified-ducts` (NICE-OPS, the validation harness, and the harness's public-lineage panam clone); none of it is pushed. `H/` below is the validation harness (`~/Code/Python/niceops-validation`).

## Propagation

### Bounced arrivals in a duct

`--refraction stratified` found one eigenray per path. In the day-286 inversion the air aloft is faster than everything between ground and source, so beyond the reach of every direct ray the energy arrives by bouncing along the ground, and the model reported a near-silent shadow there. It now sums every arrival in energy, bounced ones included, each earlier bounce costing the plane-wave |R|²; past every direct ray's reach the strongest bounce becomes the primary (NICE-OPS `e1617d0`).

| Case | Before | After |
| --- | --- | --- |
| Day-286 inversion, far third of the range, against the PE | 8–24 dB under | mostly within ±3 dB |
| Downwind regression, far third | −4 to −6.5 dB | within 1 dB |
| Run time, `example_windy.csv` / `example_inversion.csv` (demo grid, B407, `--ground_effect`) | 23.6 / 57.5 s | 37 / 86 s |

Three changes bought the cost down from +107% to +55%: splitting each bracket at the range discontinuities where 1/p crosses a record local maximum of c (the inversion's direct solve went from 9.6 to 2.8 µs a pair); sharing one 32-sample grid of launch parameters across all four ray families; and searching only bounced rays past the hand-over to the linear fit. Not fixed: a flush microphone under warm, calm ground air (run 286288) is still about 21 dB low at 8,000–13,000 ft, where rays turned back from aloft turn up again within the lowest meter and the PE diffracts through that layer.

### Board/pole pairs with curved rays

The pairs compare each ground board with the 1.2 m pole beside it, T = L_board − L_pole + G_pole, so the source cancels and only the receivers' ground interaction is left. Taking the ray geometry from each run's measured profile (`niceops_ray_geometry`, NICE-OPS `cdd1b9e`, `2dfb18b`):

| Elevation | Frames | Straight rays | Linear | Stratified | Another run's profile |
| --- | --- | --- | --- | --- | --- |
| 0–2° | 4,861 | 7.6 dB rms | 0.8 | 1.2 | 2.4–3.5 |
| 2–5° | | 2.2 | 1.0 | 1.0 | |
| 5–10° | | 1.3 | 1.1 | 1.0 | |

Per band, straight rays were 5–16 dB off at 200 Hz–1.25 kHz; curved rays leave ±1–2 dB. Variable porosity σe = 200 needs no refit: the "harder ground" fits of σ 3,200–6,400 at 0–2° were compensating for refraction. The shuffled-profile control shows it is each run's own atmosphere that helps, not curvature in general (`H/propagation/pairs_refraction.py`).

### A low source in a weak upwind gradient

At the slow end of the day-286 approaches `stratified` was 6.5 dB over straight rays at 63–250 Hz, at microphones 150–450 ft from a source 15–50 ft up, where the PE matched straight rays within 0.5 dB. Two faults (NICE-OPS `90cd576`, `0316225`):

- **Faster air aloft hid the boundary.** Any air above the source faster than the ground was read as "rays from aloft reach the receiver", so no ground-grazing boundary was set. They land only past their skip distance, kilometers out; short of it the receiver is in an ordinary shadow. With no boundary there was no transition and no hand-over to the linear fit, and the shadow took the layered shadow ray, focused 6.5 dB by the log law's lowest meters, plus the whole ground term: 12 dB over free field at 63 Hz where the PE gives 4. The boundary now stands (`boundary_under_aloft`) and is taken back where a ray from aloft or along a duct does arrive, which leaves the duct's handling exactly as it was.
- **A low source leaves the transition no room.** When a band's boundary lies within two diffraction lengths of the source, the shadow-transition table has no lit zone to hand over from, and the source is under the 5 l its table was validated to; V put the whole path 2–5 dB under the PE. Such a band keeps the ray model to its boundary and blends to V over the next diffraction length, without the creeping-wave loss in the blend (which put a 30 m source in a linear lapse 6 dB under the PE near its boundary).

| | Before | After | Straight rays |
| --- | --- | --- | --- |
| Day-286 frames, 63–630 Hz, sources under 100 ft and 20 kt | +8.9 dB | +4.4 dB | +3.8 dB |
| Day-286 frames, 2–5° within 2,000 ft | +4.7 dB | +3.1 dB | +2.4 dB |
| Day-286 PE matrix past 1,000 ft, `stratified` − PE, 200–630 / 800–2,000 Hz | 9.5 / 11.7 dB rms | 7.5 / 9.2 dB rms | |

The 24 day-284 PE cases and the steady-run frames do not move (day 284 overall +0.09 → +0.14 dB). The development cases (windy, log law, lapse) move by at most 0.1 dB rms.

### The PE at the day-283 low-pass geometry

At a fixed emission angle (3°, 5°, 8°) on three day-283 lapse profiles, sources 60–300 ft, flush and 1.2 m receivers, both ray models stay within ±1 dB of the PE at every source height, with no trend in height (`H/propagation/day283_pe/`). Modeled propagation does not make a low source look quieter than a high one at the same angle.

## Scoring

Two biases in the frame statistics were removed before the residuals below were read:

| Bias | Mechanism | Size |
| --- | --- | --- |
| Measured-peak gate | Keeping frames within 10 dB of the measured peak selects on the measurement, so quiet-by-chance frames drop out at long range | about 1 dB at long range |
| Decibel medians | A 0.5 s third-octave frame fluctuates by σ = 2–4 dB, and its median in dB sits about 0.115 σ² below the energy mean the SEL integrates | 0.6–1.9 dB |

Frame scores are energy means with no gate on the measured level (`band_frames_summary.energy_error`, `energy_slope`). Several conclusions changed on rescoring; an older number in a project note that predates this should be treated as provisional.

## The spheres

### Depropagating along curved rays

Straight-ray depropagation is good to about 0.5 dB at a 10° rim, ±1.5 dB at 5–10° and up to 10 dB below 5° (`H/propagation/depropagation_refraction.py`). panam's `depropagate_hemisphere` now takes a `ray_model`: `niceops_ray_geometry` returns per pair the launch depression (where the sample is filed), the tube's range (spreading), the arc length (absorption) and the arrival angle (the board's response). panam depends on NICE-OPS only through that executable and a CSV; the default stays straight rays.

### Ground-board correction near grazing

A curved ray reaching a flush board within half a degree of grazing had the board's modeled response divided out, and that response falls to −10 dB at 2°, −16 at 1°, −26 at 0.3° and −41 at 0.05° (4 kHz). Behind four day-283 low descents the shipped spheres carried 3–4 kHz levels of 113 dB where straight rays give about 80, and an adaptive grid exposed up to 123 dB. panam now drops a band needing more than 15 dB of board correction (`max_response_correction_db`, panam `ea119c2`), as it already dropped a band needing more than 30 dB of absorption correction.

### Gridding the rim

The hemisphere grid was a 10° elevation step with Shepard interpolation over a fixed 25° radius. It is big enough to fill overhead, where samples are up to 18° apart, and far too big for the rim, where they are under a degree apart: on low descents the samples fall about 6 dB from 8–10° to 2–4° below the horizon, and the grid read 2.4–2.9 dB loud at 2–4°. The written sphere and NICE-OPS's bilinear reading add nothing further (`H/spheres/rim_samples.py`).

An adaptive grid (panam `4a8debd`) gives each node a radius from its own samples: R = max(1.3 d₈, 1.3 d₃ mics, 0.5 s), capped at 60°, d₈ the distance to the 8th sample, d₃ mics the distance within which three microphones contribute, s the sample's angular resolution. A node no radius under the cap can bracket with three microphones retries with fewer before it is a gap (panam `ea119c2`); every sphere with gaps now misses fewer directions than the fixed grid did (286318: 54% against 68%).

That grid let each node follow the nearest microphone's line of samples, because Shepard's weight grows without bound at a sample, so each microphone's own offset (calibration, installation, local ground) printed as ripples along its track. Per-microphone flat corrections were rejected: there is no concrete evidence of a specific bias, and they would fit noise. The lines run nose to tail, so the ripples lie across azimuth, while the rim's steep change is in elevation. A neighborhood stretched in azimuth averages across the lines and keeps that gradient, but only once the weight stops interpolating (panam `83f989f`):

| Grid, low descent 283145 | Azimuthal ripple | Offset spread left in the samples | Grid − samples at 2–4° | at 4–6° |
| --- | --- | --- | --- | --- |
| Fixed 25° radius | 0.52 dB rms | 0.43 dB | +0.95 dB | +0.60 dB |
| Adaptive Shepard | 0.72 | 0.36 | −0.24 | +0.03 |
| Shepard, azimuth reach ×5 | 0.65 | 0.39 | −0.13 | −0.02 |
| Biweight, ×5 | 0.48 | 0.50 | +0.10 | +0.26 |
| **Shepard, ×5, distance floored at the radius (shipped)** | **0.52** | **0.46** | **−0.03** | **+0.16** |

At 10–40° every variant stays within 0.25 dB of its samples ahead, abeam and behind, so the stretch does not smear real azimuthal directivity.

### Databases against the day-286 approaches

| Database (NICE-OPS commit) | SEL rms, station / `stratified` | at 5–10° peak elevation | Frames, under 100 ft and 20 kt | Frames, 2–5° at 4,000–8,000 ft |
| --- | --- | --- | --- | --- |
| Curved rays, fixed 25° grid (`56ce39b`) | 1.31 / 1.42 dB | 2.10 / 2.28 dB | +4.45 dB | −0.13 dB |
| Adaptive Shepard, board cap (`041f39b`) | 1.28 / 1.39 | 2.02 / 2.20 | +3.90 | −0.55 |
| Biweight, ×5 (not shipped) | 1.28 / 1.40 | 1.92 / 2.07 | +3.91 | −0.80 |
| **Floored Shepard, ×5 (`5d761b8`, shipped)** | **1.26 / 1.37** | **1.89 / 2.03** | +3.89 | −0.81 |

SEL is paired over 1,112 microphone-runs with the 10 dB-down window. The sharper rim helps at short range and costs a few tenths at long range: the rim is a short-range against long-range trade-off, and the next section says why. `self_test` passes on every build.

### Two rim experiments that failed

- **Rim from all ranges.** Taking samples below 14° from up to 10,000 ft instead of 2,000 made the long-range forward shortfall worse (−2.2 → −2.6 dB at 4–6°), put the near shallow angles 9 dB loud and took SEL rms to 1.84 dB (panam `ce0b16e`). Far grazing samples depropagate too loud, which the 2,000 ft cap protects against.
- **Height against range.** Passes flown high imply a forward rim 1–3.5 dB louder than low passes of the same card at the same angle. Fitted over 4–45° with a constant per card and angle, neither a height term (+0.02 ± 0.22 dB per 100 ft on the clean level group) nor a range term (−0.34 ± 0.92 dB/km) is resolved (`H/spheres/altitude_pairs.py`, `height_vs_range.py`).

## Source-side residuals

### Fast high flight at long range

Fast (80 kt and up) high (400–2,000 ft) states heard 2–10° below the horizon from far ahead are 1–2 dB quiet, flat from 2,000 to 12,000 ft and across 63–630 Hz, the same under straight and curved rays: an emission-angle effect, not a range or propagation one. Ruled out for it:

- **The database blend.** The approaches' fast segment (cards A33, A35, A40, A42) blends real ~80 kt spheres from days 283 and 285, under 1% on the synthetic extended-angle copies, whose forward rims agree within about 1 dB (`H/spheres/approach_blend.py`).
- **Day 286 alone.** Day-286 spheres are louder at the forward rim than other days' at the same condition (+0.47 ± 0.25 dB over a (μ, γ) surface; +1.17 ± 0.42 on L4, the only card flown every day), shrinking with elevation, the shortfall's shape. But a database of day-286 spheres only moves the shortfall by −0.2, +0.3 and +0.7 dB at 4–6, 6–8 and 8–10°: it helps the steeper rim, not the 2–6° where the shortfall is largest.
- **Gridding, the range cap and the output stage** (above).

SEL by speed overhead, energy means, predicted − measured: 40–60 kt −0.3 dB, 60–80 kt +0.7, 80–100 kt −0.4, 100–140 kt −1.1.

### The hover

NICE-OPS pointed the sphere's nose along the velocity, which at hover is a few tenths of a knot of drift in any direction. Below 5 kt the nose is now the heading (the AC format's column 8), turned into the local frame by an offset inferred from the track's fast states or from `--frame_bearing` (NICE-OPS `9fc714a`). No shipped track and no day-286 approach goes that slow (the approaches bottom out at 13 kt), so neither moves.

The database's hover is synthesised from the slowest level pass averaged fore to aft. The day-283 hovers (H1–H12, never used before) were flown over the center of a 500 ft ring of 13 boards, low (4–30 ft, heard at 0.5–3° below the horizon) and high (165–235 ft, 18–25°), at four headings (`H/hover/`). Against them, with the nose on the heading:

| Correction, measured − predicted at 18–25° | Within 45° of the nose | Sides and aft |
| --- | --- | --- |
| 63 Hz | −8 to −18 dB | −7 to −19 dB |
| 315–1,250 Hz | −3 to +7 dB | +4 to +13 dB |
| 1.6–5 kHz | −4 to +3 dB | −1 to +6 dB |

The two passes of each card agree within 0.4–1.7 dB. A 40 kt pass averaged fore to aft is quietest exactly where a hover is loudest. panam's `build_empirical_database(hover_correction=…)` (panam `c0e8639`) blends the correction into the synthesised hover: as measured from 15°, none at 10° and below, fading to its azimuth mean from 30° to the pole. The high hovers then score within ±1.7 dB in every band and azimuth bin, against −8 to +13 before. The low hovers and the approaches' slow end do not move: they are heard 1–5° below the horizon, on the rim the correction leaves alone. It was built on the earlier spheres and is **not shipped**: above 30° it is extrapolated. The flight card (`B407_Flight_Cards/FLT283_B407_static.xlsx`) names every card:
- H1 and H12 are ground idle, and H2 and H11 flight idle. All four are on the ground, not flight conditions, so their errors are expected: H1/H12 about 16 dB loud in every band, H2/H11 about 11 dB quiet above 800 Hz. They belong neither in the correction nor in its scoring.
- H3–H6 are hover in ground effect (the low hovers) and H7–H10 hover out of ground effect (the high hovers). The correction is built from the latter.
- The card's relative direction (90 = nose toward mic 7, increasing clockwise) agrees with the AC heading on all 20 runs. The heading reads a median of about 9° clockwise of the card on the nose-to-mic-7 runs, the INS offset seen in flight (+4.5° on the B407), so the heading-based nose is confirmed.
- The card's "average temperature (F)" column reads 22–35, which is °C. Its humidity of 42–54% is the hovers' own; the station configurations should be checked against it.

### The day-283 low passes

The day-283 descents (D4, D24, D14, D12) cross the array 55–112 ft up, with profiles. At a fixed emission angle forward, the same steady pass implies a source 3–5 dB quieter from under 75 ft than from 200–400 ft, under straight and curved rays alike (`H/spheres/low_pass.py`). Separating height from range relies on the spread of angles, and the two methods disagree: the depropagated samples put the trend on height (+1.1 to +1.4 dB per 100 ft at 200–630 Hz, no range term), the frames mostly on range (−1.5 dB/km). It is not a site effect (per-microphone constants move the slopes by at most 0.15 dB per 100 ft), and per the PE above it is not modeled propagation.

### Poles against boards

The three 1.2 m poles (50, 51, 52) stand beside boards 34, 36 and 38 at x = −3,759 ft. Scored frame by frame beside them, as point receivers at their fitted heights over the site ground with the pair fit's turbulence (`H/propagation/band_frames.py --poles`, `pole_board.py`; 462,000 paired frame-bands from the 25 approaches and 52 steady runs), pole minus board is within ±0.5 dB in nearly every cell, near grazing included:

| Energy-mean error, board / pole | 63–160 Hz | 200–630 Hz | 800–2,000 Hz |
| --- | --- | --- | --- |
| 0–2°, 2,000–4,000 ft | +2.7 / +2.8 dB | +3.7 / +3.3 dB | +1.2 / +1.7 dB |
| 2–5°, 2,000–4,000 ft | −0.4 / −0.5 | −1.0 / −1.4 | −3.6 / −4.1 |
| 10–20°, under 2,000 ft | −0.2 / +0.0 | −0.7 / −1.1 | −0.1 / −0.1 |

The near-horizon error is common to both receivers.

## What is ruled out, and what remains

| Candidate | Verdict | Evidence |
| --- | --- | --- |
| Refraction at 0–2° (straight rays) | Was the cause; fixed | Pairs 7.6 → 0.8–1.2 dB rms |
| Ground model (VP σe = 200) | Stands | Pairs with curved rays; no refit needed |
| `stratified` duct and low-source faults | Fixed | PE comparisons above |
| Modeled propagation for the remaining residuals | Ruled out | PE within ±1 dB at fixed angle, all heights |
| Scoring artifacts (gate, dB medians) | Removed | 1–2 dB |
| Board correction near grazing | Fixed (cap) | Aft 3–4 kHz spikes gone |
| Rim gridding | Fixed | ±0.2 dB of the samples; ripple below the fixed grid's |
| Database blend, range cap, output stage | Ruled out | |
| Day-286 source offset | Partial | +0.5 to +1.2 dB at the steeper rim only |
| Microphone-line site effect | Ruled out | Per-microphone constants |
| Boards' absolute response near grazing | Ruled out | Poles carry the same error |
| Hover sphere at 18–25° | Corrected, not shipped | ±1.7 dB in-sample |
| **Source levels near the horizon, or an unmodeled path effect** | **Open** | Low/near sources 3–4 dB loud; fast/far 1.5–2.8 dB quiet at 2–6° |

The two open residuals have opposite signs and sit in the same angular band: near the horizon a source at short range and low height comes out too loud, one at long range and high too quiet. That is the signature the rim trade-off kept showing, and it is what a rim built from short-range samples would carry if the attenuation near the ground grew with range or height faster than any model here predicts. Turbulence, ground that varies between source and microphone, and terrain are the candidates outside the PE's physics; the source's own near-horizon directivity, which every sphere is built from the same short-range samples, is the other.

## Next

1. **Separate source from path.** Compare microphone lines across the array at matched emission angle and range: an error that moves with the line is the path; one that stays with the emission angle is the source.
2. **The hover.** Decide whether to ship the hover correction re-applied to the shipped spheres (`spheres/hover_sphere.py --spheres spheres/out/min2_stratified_shepfloor5/Be407`). The cards are now known (above). The low hovers are in ground effect, a rotor condition the synthesised hover does not represent, so their +3 to +7.5 dB at 200–630 Hz need not be the rim alone.
3. **The flush microphone under warm, calm ground air** (286288), 21 dB low at 8,000–13,000 ft, which the ray models cannot reach.
4. **panam upstream.** The ray-model depropagation, the board-correction cap, the adaptive grid and the hover correction are on the harness's public-lineage panam clone only; they are generic code and would go to panam master, not beta26.

## Source or path: the microphone lines (2 October, evening)

Next step 1 is done. Both open residuals follow the emission direction and flight state, not the microphone line. The lines do carry offsets of their own, about ±0.5–2 dB, but those offsets are not the residual.

**Method** (`H/propagation/line_vs_angle.py`, output `line_vs_angle.txt`; shipped database `5d761b8`, `stratified` with plates, boards only):
- A cell fixes the source: run, band, emission elevation (1°), sphere phi (10°) and theta (5°), range (500 ft), and on an approach ground speed (10 kt). Every line that hears a cell hears the same source from the same direction at the same distance.
- A line's effect is its energy-mean error less the mean over the lines in the cell.
- The noise level comes from shuffling line labels among the frames of each cell.

**The lines carry real offsets.** Line-effect rms is about twice the shuffled value in every band and elevation bin: 1.1–2.8 dB against 0.7–1.7. The pattern repeats across the 52 steady runs and the 25 approaches. At 200–630 Hz and 5–10°, 8 of 9 lines have the same sign in both sets:

| Line x (ft) | +608 | +425 | 0 | −1002 | −1905 | −2808 | −3759 | −5115 | −6087 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Steady runs | +1.6 | −0.7 | +0.6 | −0.4 | −0.1 | +0.5 | +0.1 | −0.7 | +0.6 |
| Approaches | +1.8 | +0.0 | +0.1 | −0.4 | −0.2 | +0.5 | +0.1 | −1.1 | +1.1 |

The offsets shrink above the rim (+608: +0.9 to +1.2 dB at 10–20°; −5115: −0.3 to −0.6). An installation offset would not depend on angle, so this points at the ground along each line's paths more than at the microphones. The offsets sum to zero across lines by construction, so they cannot carry a residual whose sign is common to every line.

**Fast, high flight heard far ahead** (≥80 kt, 400–2,000 ft, beyond 4,000 ft, 2–6°, 63–630 Hz):
- **Approaches:** −2.34 dB overall. Every forward line shows it, from −3.4 dB (+425) to −1.1 dB (−2808).
- **Steady passes, same geometry and lines:** +0.1 dB overall.
- **The same runs heard near and steep** (under 2,000 ft, 10–30° ahead): within ±0.7 dB on every approach run.
- **By run:** the shortfall is in all ten fast approach runs (−1.2 to −4.1 dB). Among steady runs it appears only in descents 283144 and 283145 (85 kt, −2.0 and −3.7 dB) and day 284's 1,000 ft passes 284209 and 284210 (−1.6 dB). The other seven fast steady runs are −0.6 to +3.0 dB.

The shortfall is common to every line, absent at the same geometry on most steady passes, and present in particular runs. That makes it a property of the source state, not of the path.

**The slow end of the approach heard near the horizon** (under 100 ft, under 40 kt, within 2,000 ft, 1–5°): +2.95 dB overall, and its sign follows the direction of emission.
- Lines hearing it ahead: +3.3 and +4.0 dB.
- Line hearing it abeam: +2.3 dB.
- Lines hearing it abeam and aft, or aft: −0.5 and −1.8 dB.

Paths ahead and behind cross the same kind of ground at the same range. A sign that changes from fore to aft is the source's directivity at that state: decelerating toward a hover, which no sphere was flown at.

**What this leaves.** Both residuals are on the source side.
- **Slow end of the approach:** this is the missing decelerating, near-hover condition. The hover correction does not reach it, because it leaves 1–5° alone.
- **Fast, high flight far ahead:** the forward rim of a fast state heard from far off is 1.5–4 dB louder than the rim its own near, low samples imply. That holds for the approaches and for some steady runs, notably the low descents. The near samples set the rim when the aircraft is low; the far frames hear it high. Lines cannot separate a height dependence of the source from a gain that grows with range near the ground. Frequency can: the shortfall is flat across 63–630 Hz, which argues against scattering and refraction and for the source.

**Suggested next:** a forward-rim term in source height, for example rim samples weighted toward the height of the state being predicted. Alternatively, accept the forward rim as a ±2 dB uncertainty for fast flight.

The per-line offsets could be corrected, but they are ±1 dB at the rim and they average out of footprint-scale results. A correction is worth fitting only if site maps show surface differences along the +608 and −5115 lines.

## Code, data and commits

| Repository | Commit | What |
| --- | --- | --- |
| NICE-OPS | `e1617d0` | Every arrival under `--refraction stratified`, bounced ones in ducts |
| NICE-OPS | `cdd1b9e`, `2dfb18b` | `niceops_ray_geometry`: refracted-ray geometry for a batch of pairs |
| NICE-OPS | `56ce39b` | B407 databases depropagated along curved rays |
| NICE-OPS | `90cd576`, `0316225` | Boundary short of the skip distance; room for the transition |
| NICE-OPS | `9fc714a` | Sphere nose on the heading at hover |
| NICE-OPS | `041f39b` | B407 databases: adaptive grid, board correction capped |
| NICE-OPS | `5d761b8` | B407 databases: floored Shepard, azimuth reach ×5 (shipped) |
| panam | `23b8bcf`, `8143c63` | `ray_model` in depropagation; samples carry range and height |
| panam | `4a8debd`, `ce0b16e` | Adaptive radius; `rim_range` |
| panam | `c0e8639` | `hover_correction` in `build_empirical_database` |
| panam | `ea119c2` | `max_response_correction_db`; relaxed microphone count |
| panam | `200cb20`, `83f989f` | Biweight kernel and azimuth aspect; `shepard_floor` |
| harness | `0bbdee2` … `058862e` | Scripts: `propagation/` (pairs, frames, PE matrices, poles), `spheres/` (rebuilds, rim, blends, low passes, plots), `hover/` |

Sphere builds are in `H/spheres/out/` (`min2_stratified_shepfloor5` is the shipped one; `hover_corrected` holds the hover correction on the earlier spheres). Frame sets are in `H/propagation/band_frames_*`, SEL runs in `H/propagation/out_*`.

## Figures

Lower hemisphere seen from below (nadir at the center, horizon on the rim, nose up), one condition under every build recipe tried, and rim cuts ahead, abeam and behind (`H/spheres/plot_methods.py`):

- Low descent 283145 (D4, 80 kt, −3°, 63 ft over the array): [maps](Grazing%20incidence%20findings%202026-10-02%20-%20figures/283145_maps.png), [rim](Grazing%20incidence%20findings%202026-10-02%20-%20figures/283145_rim.png). The aft spike in the uncapped curved-ray builds and the rim's steep fall-off are both visible.
- Level pass 283107 (L3, 80 kt, 248 ft): [maps](Grazing%20incidence%20findings%202026-10-02%20-%20figures/283107_maps.png), [rim](Grazing%20incidence%20findings%202026-10-02%20-%20figures/283107_rim.png).
- The same samples gridded four ways, dBA contours, showing the microphone-line ripples and their removal (`H/spheres/plot_smoothing.py`): [smoothing](Grazing%20incidence%20findings%202026-10-02%20-%20figures/smoothing_dba.png).

## Addendum: source against path, by microphone line

Within one run the same emission direction at the same range reaches several of the array's cross-track lines from one source, so a constant per (run, 1° elevation, 30° azimuth, 500 ft range) cell plus a constant per line separates what a line's path or site adds from what the source and geometry share (`H/propagation/line_effects.py`). The line constants spread 0.2–0.3 dB at 10–40° below the horizon and 0.7–1.7 dB at 2–5°, on the day-286 approaches and on the steady runs alike: near the horizon there is a site-dependent path component. Its one consistent piece is microphones 42–44 (x = −5,115 ft, y = 0 to −800), measured 2–6 dB louder than modeled at 200–2,000 Hz near grazing in both sets, which is what a patch attenuating less than the modeled ground would give.

The two large residuals are not that. The fast high far-forward shortfall is on every line (−1.3 to −4.1 dB), larger on the lines furthest ahead of the aircraft: broad, and growing with distance. The low slow excess sits on the two lines at the approach end (x = 0 and +608, +3.3 to +4.6 dB at 63–630 Hz) and not on the lines behind, where direction and line are confounded; where they overlap both matter (ahead at 200–630 Hz, +4.7 dB at x = +608 against +1.9 at x = 0; at 63–160 Hz ahead and abeam +2.5 to +4.3, behind −0.5). A uniform path effect, such as ground everywhere harder than modeled near grazing, cannot be told from the source this way; a second aircraft over the same site, or the same aircraft over different ground, can.

## Addendum: a second aircraft (B206L3, day 279)

Today's recipe does not generalize as an overall improvement, and the near-horizon pattern does. The B206L3's 23 day-279 approaches were scored against two databases rebuilt without any approach card: the September recipe (straight rays, 10° rim, 25° grid) and today's (curved rays, 2° rim, floored Shepard ×5, board correction capped).

| B206L3 SEL rms, station / `stratified` (906 mic-runs) | All | 5–10° | 10–20° | 20–40° | 40–90° |
| --- | --- | --- | --- | --- | --- |
| September recipe | 1.48 / 1.47 dB | 2.36 / 2.40 | 2.23 / 2.19 | 1.33 / 1.31 | 1.09 / 1.08 |
| Today's recipe | 1.50 / 1.49 dB | 2.11 / 2.16 | 2.29 / 2.27 | 1.38 / 1.37 | 1.13 / 1.11 |
| Shipped database (`7a46a72`, contaminated) | 1.41 / 1.40 dB | 2.28 / 2.34 | 2.18 / 2.14 | 1.24 / 1.23 | 1.01 / 0.99 |

Today's recipe gains at 5–10° (0.25 dB rms), as on the B407, and loses about 0.05 elsewhere. In the frames the B206L3 shows the B407's near-horizon signature at smaller size: short range at 2–10° +1.6 to +2.1 dB, low states +1.5 to +1.9, long range at 2–5° −0.7 to −1.7; the sharper rim moves it the same way (short range better, long range more negative). Two aircraft, two days, one site and one method share it.

**The shipped B206L3 database contains 12 of its 23 validation approaches.** The B206L3's approach cards carry explicit zeros for bank and acceleration, so `is_steady_flight_card` counts all 23 steady, and the 12 with an 8 s steady stretch were built into spheres (`7a46a72`, 27 September). The R66 has 2 such cards; the B407, AS350B3, EC130B4 and R44 have none. `H/spheres/rebuild.py --exclude-conditions A` rebuilds without them; the in-sample advantage is about 0.07–0.1 dB rms.

On the B206L3's decelerating approaches the line breakdown cannot separate site from source: at a fixed emission cell each line hears a different stage of the approach, so the line constants form a gradient along the track (x = 0 and +608 −1 to −3 dB, x = −6,087 +1.2 to +5) at every elevation. The B407's steady passes, whose state does not change along the track, are the clean evidence for the site effect at x = −5,115.

## Addendum: run 286288 is an elevated duct, not a warm layer (2 October, night)

The flush-microphone failure listed under Next item 3 was investigated with the PE, profile variants and the day-286 frames (harness `propagation/micro_shadow.py`, `micro_shadow_frames.py`; NICE-OPS `5d761b8`, shipped B407 database).

**The profile.** Toward 270 the effective sound speed is highest at the ground, 334.13 m/s. It falls 0.21 m/s in the lowest foot and 0.42 m/s by 20 ft, to a minimum near 230 ft, and rises to 333.92 m/s at the inversion's top near 330 ft. The near-ground fall is the log-law headwind assumed below the LIDAR's 46 ft gate (z0 = 0.1 ft); the temperature *rises* 0.18 K over the lowest 10 ft. So it is an elevated duct. Its returning rays turn back up 1–20 ft above the ground, and ray theory delivers them to neither receiver.

**Against the PE** (100 ft source, mean of model − PE over 8,000–13,000 ft):

| Band | PE, dB re free field (flush) | `stratified` − PE, flush | `stratified` − PE, 1.2 m | straight − PE, flush |
| --- | --- | --- | --- | --- |
| 63 Hz | +0.7 | −19.3 | −10.9 | −3.0 |
| 250 Hz | −4.9 | −21.2 | −15.7 | −6.1 |
| 1 kHz | −22.3 | −14.2 | −11.9 | +6.2 |
| 2 kHz | −28.4 | −13.2 | −13.4 | +9.5 |

The 1.2 m microphone fails too; it is lit only at 12–13 kft.

**Two mechanisms, separated by profile variants.** Wind below the gate was held at its 1 ft, 10 ft or 46 ft value, or set calm.
- *Tunneling.* The PE barely notices the near-ground layer. Holding the wind flat below 46 ft moves it by under 1 dB at 63–250 Hz and by about 6 dB at 1 kHz, 10,000–13,000 ft. A WKB evanescent factor, exp(−2ω∫κ dz) with κ = √(p² − 1/c(z)²), gives 0.1–0.4 dB at 63 Hz and 1–7 dB at 1 kHz for rays turning 1–9 ft up, and tens of dB for rays turning above 20 ft. So the PE's return is carried by the rays that turn lowest, and the leak is the turning-point continuation the methods report proposed.
- *The skip zone.* With the shear removed, the duct's rays reach the ground, but only from about 12 kft. Short of that no ray arrives, the boundary is infinite, and `stratified` falls to −55 to −80 dB where the PE is −10 to −48 dB. The PE's return starts 1–2 kft earlier, filled from the caustic. With the measured profile the shadow transition's floor happens to sit near the PE in the skip zone, and the return is what fails.

**Against the measurement.** Frames for all 25 approaches were produced under straight and layered rays with the same build and database. Where `stratified` predicts at least 10 dB under straight rays (46,677 of 3.77 million frame-bands, in every run, a median of 0.4–2° and 4–9 kft from sources under 300 ft), the energy-mean errors are:

| Band | Straight rays | `stratified` |
| --- | --- | --- |
| 50–160 Hz | +3.8 dB | −10.5 dB |
| 200–630 Hz | −0.1 | −15.5 |
| 800–2,000 Hz | −7.1 | −20.6 |

The SNR gates change these by under 1.5 dB. The data sit with the PE: the ray model's deep shadows and skip zones under day 286's profiles are 10–20 dB deeper than measured.

**It does not move the SEL.** None of these frames falls within 10 dB of its microphone-band's measured peak, except 11 windows at 800–2,000 Hz. Approach exposure is set near and overhead, and per-band exposure errors are the same under both ray models to 0.02 dB. The gap matters for a low source heard far off: hover, low transit, take-off and landing, and footprint contours far from the track.

**What a fix needs.**
1. A tunneled family in `stratified`: rays whose lower turning point lies above the receiver, landed at their turning range with the per-band WKB (or Airy-uniform) factor. The factor is per band, so it has to reach the band loop as an extra loss, not as ray geometry.
2. A caustic and skip-zone treatment: an Airy fill short of the duct's first landing, of the same kind `ShadowTransition` gives a ground-grazing boundary.
3. Alternatively, a PE table for low sources (`--propagation_table`), which carries both and is already plumbed.

## Addendum: three treatments of the duct's micro-shadow, compared (3 October)

All three are built, off by default, uncommitted on `stratified-ducts`: `--tunneling` and `--skip_fill` under `--refraction stratified` (Airy tails below a turning point and into a caustic's dark side, `layered_wave_arrivals` in `src/stratified.cpp`), and a PE `--propagation_table` per run, which now runs beside `--refraction` and reads plate microphones. Harness scripts: `propagation/micro_shadow_fixes.py` against the stored PE, and `micro_shadow_frames.py` against the measured frames. The tables are in `propagation/micro_shadow/tables/`.

**Against the PE** (100 ft source, flush, 2,000–13,000 ft, five day-286 profiles; `stratified` − PE, mean / rms):

| | 50–160 Hz | 200–630 Hz | 800–2,000 Hz |
| --- | --- | --- | --- |
| neither | −5.4 / 9.8 | −6.1 / 11.1 | −3.6 / 9.9 |
| `--tunneling` | −3.4 / 7.3 | −3.8 / 8.5 | −1.4 / 7.3 |
| `--skip_fill` | −4.5 / 8.9 | −5.6 / 10.5 | −2.7 / 9.2 |
| both | −1.7 / 5.1 | −2.3 / 6.5 | +0.3 / 6.7 |

Sources from 330 ft up do not move: there the source sits at or above the duct's top, and no family turns back from aloft. On 286288's profile with the shear below 46 ft removed, both together take 50–630 Hz from −13 to −21 dB to −2 to −4. 800–2,000 Hz stays 7 dB low, at the skip zone's near edge, where a duct's direct family ends without a caustic.

**Against the measurement**, on the frames where `stratified` predicts at least 10 dB under straight rays (energy mean, predicted − measured):

| | straight | `stratified` | `--tunneling` | `--skip_fill` | both | PE table |
| --- | --- | --- | --- | --- | --- | --- |
| 25 runs, 50–160 Hz | +3.8 | −10.5 | −7.3 | −3.6 | −2.1 | |
| 25 runs, 200–630 Hz | −0.1 | −15.5 | −10.7 | −7.2 | −5.0 | |
| 25 runs, 800–2,000 Hz | −7.6 | −20.6 | −20.6 | −20.7 | −18.9 | |
| 5 runs, sources under 150 ft at 3–6 kft, 50–630 Hz | +1.0 to +1.6 | −11.5 to −13.3 | | | −2.5 to −4.5 | +0.3 to +2.5 |
| 5 runs, other frames beyond 4,000 ft under 5°, 50–160 / 200–630 Hz | +2.3 / +1.4 | +4.5 / +4.2 | +4.5 / +4.2 | +4.8 / +4.4 | +4.8 / +4.4 | +2.4 / +1.5 |

- **The PE table is the most accurate where it applies.** It also removes `stratified`'s +2 to +3 dB over-prediction on the far frames that are not in deficit.
- **Run 286276's far frames are wrong for every model but plain `stratified`.** These are 150–330 ft sources at 6–20 kft, and the table reads 13.5 dB too loud there, straight rays 8.8. The cause is outside the PE: turbulence, the profile's duct, or the source.
- **None of the options moves the SEL.** Within 10 dB of each microphone's measured peak, the exposure error is the same under all six configurations to 0.02 dB (rms 2.1 / 1.4 / 1.7 dB by band group on the five runs).

**Cost on a day-286 run**, timed on a loaded machine and so approximate: `--tunneling` about +20%, `--skip_fill` 2–4×, and the table about the ray model's own. Each table (8 bearings × 8 source heights to 330 ft, 50 Hz–2 kHz, 16 kft) takes 2.5 minutes. Read between its nodes, the table is within about 1 dB of a direct PE run at 63–250 Hz but 5 dB rms overall in this duct, almost all at 1–2 kHz.

## Addendum: turbulence does not explain run 286276's far frames (3 October)

There, 250–300 ft sources heard at 1.9° about 8,000 ft ahead, the PE table read 15–20 dB above the measurement at 125–800 Hz. The PE was rerun on the run's profile with von Kármán screens, γT from 1e-6 to 1e-4, outer scales of 10 and 50 m, 40 realizations a band (harness `propagation/micro_shadow_turbulence.py`). Against the still field marched the same way, turbulence *raises* 250 Hz–1.25 kHz by 1–5 dB at every strength, filling the ground's interference. It lowers 63–125 Hz by at most 2.5–3 dB, and only at 1e-4. So it cannot close a gap of the opposite sign.

The measured field has neither model's shape. Relative to straight rays it falls smoothly from −3.6 dB at 5 kft to −11.6 at 7–7.5 kft, then recovers only to −6 to −7 at 9–10.5 kft. The profile's inversion puts the PE 4–9 dB *above* straight rays over the same span. The layered rays have a sharp skip zone at 8–8.5 kft and a caustic at 9–10 kft, 15 dB above the measurement there.

It is not one run. On the seven mid-morning runs with such frames (286274–286280, sources 150–400 ft at 6–10 kft, 200–630 Hz, SNR-gated), the measurement is −6 to −8 dB against straight rays in five of them, while their profiles imply +3 to +5. The correlation between the two is −0.19. The profile's refraction does not show in these data, and the deficit does not follow it. That points at something fixed to the geometry: ground or obstructions along the far path, or the source's level 1–3° below the horizon. Turbulence and a different propagation method do not fix it.

## Addendum: the far-frame deficit is not the terrain (3 October)

The source was public data: the USGS 3DEP 1 m lidar DEM, resampled to 2 m over the approach corridor (6 × 2 km, harness `propagation/terrain_paths.py`, DEM in `propagation/terrain/`). It reproduces the reference point's ground to 0.1–0.3 ft. Each far frame's source was placed on the track at its emission time, and the line of sight to the microphone was compared with the terrain every 2 m.

- **The corridor is a lakebed.** It has 9 m of relief in all, and the ground rises about 3 ft between the array and a source 8,000 ft out. A ray arriving at 1.9° is already 33 ft up 1,000 ft from the microphone.
- **No path is obstructed.** On every path the worst Fresnel parameter is under 0 (median −0.19 at 6–10 kft, −0.32 at 3–6 kft). The worst point sits within about 40 ft of the board, which is the ground interaction the ground model already carries.
- **Nothing is built in the way.** OpenStreetMap maps no buildings in the corridor. The paved runway 09/27 (46 m wide) lies 2,600–4,500 ft north of the track and crosses none of these paths.
- **The deficit varies run to run, which no fixed feature can do.** Over the same terrain, the same microphones and the same card, 286277 reads +3.8 dB against straight rays and 286278 −6.6. Over the seven runs the far frames span −8.4 to +3.8 dB, while 3–6 kft frames from the same runs and source heights stay within −1.4 to +1.6.

Ruled out: terrain, mapped structures and paved surfaces (above), turbulence (previous addendum). Left: a time-varying atmosphere the single per-run profile does not capture, such as horizontal inhomogeneity along a 2–3 km path or the profile's timing, and the source at 1–3° below the horizon.

## Addendum: the atmosphere hypothesis, tested (3 October)

The far-frame deficit on runs 286274–286280 (150–400 ft sources heard 6–10 kft ahead) was tested against the atmosphere in five ways. Harness scripts: `propagation/micro_shadow_ensemble.py`; the per-frame splits were run inline on `micro_shadow/frames_*`.

1. **The profiles are interpolated, and the real inversion was stronger.** Over these runs the tethered balloon sat at about 440 ft and read 38–45 °F, swinging 2–4 °C in minutes. The runs' 0–500 ft profiles hold about 38 °F there, and their ground temperatures climb smoothly, so they come from ascents before and after. Against the three ground stations, the real-time inversion was +6.4 to +10.8 °F, 2.5–5.7 °F stronger than the profiles on 286277–286280. The stations differ among themselves by about 6 °F (SWS2 colder). The deficit follows none of these: 286277 and 286278 have similar inversions and read +3.8 and −6.6 dB. Since a stronger inversion means more downward refraction, this goes the wrong way.
2. **The deficit is flat in frequency.** On 286274–286276 it is −7 to −9 dB across 50–630 Hz and −4 to −7 dB at 800–2,000. A refraction shadow or a near-ground layer deepens with frequency. (286277 reads +17 dB at 800–2,000 Hz, which looks like contamination.)
3. **Same emission time, near microphones fine.** For the aircraft at x = −9,000 to −6,000 ft, frames 3–6 kft away read +0.1 to +1.6 dB against straight rays, while those 6–10 kft away read −3.5 to −5.4. So the deficit is not the flight state at emission. It follows the path or the emission angle.
4. **An ensemble over the inversion's measured spread smears the caustics but does not lower the level.** The profile above the ground was warmed or cooled by −2 to +6 °F at 440 ft, tapering to 0 at the ground. The ensemble-mean PE comes out at +3 to +6.5 dB against straight rays across 5–11 kft, where the measurement is −3.6 to −11.6 dB. Even the least-refracting member is 6 dB above the measurement at 7–8.5 kft.
5. **A warm surface layer the profile misses does not do it either.** Adding +2 to +6 °F at the ground over 20–100 ft moves the PE by under 3 dB at 200–630 Hz, leaving it 6–15 dB above the measurement. Only the layered rays throw a shadow there (+6 °F: −6 to −20 dB), because they resolve a layer the waves pass through.

A correction to the framing. Straight rays are a poor yardstick near grazing: at 0.5–1.5° and 2–4 kft the measurement is 13 dB above them. Under `stratified`, the error by emission elevation and range is within ±3 dB in most cells, the worst well-populated one −4.2. So the "−6 to −8 dB" was partly the yardstick. Under both Airy options (`--tunneling --skip_fill`), the 8–14 kft cells at 0.5–2.5° are 6–7 dB too loud, where plain `stratified` is within 0–4 dB. The deficit-frame score, which selects on `stratified`'s own shadows, did not show that.

**Verdict:** a range-independent atmosphere does not explain the far frames. That covers the measured inversion, its measured fluctuation, a missing warm surface layer and turbulence (previous addenda). Not tested: horizontal inhomogeneity along a 2–3 km path, which a range-independent PE cannot represent. Left: the source's forward rim at 1.5–2.5° from a 100 kt, 270 ft state, and run-to-run source variation, since consecutive runs on one card differ by 10 dB.

## Addendum: the rim against range, all ranges (3 October)

The rim is built from samples within 2,000 ft. To see what lifting the cap would bring in, each pass was depropagated out to 12,000 ft and down to 1° (harness `spheres/rim_range_fit.py`). Samples were then fitted as level = c(speed and glide-path group, 1° elevation, 30° azimuth) + r(run) + g(range bin), with g relative to 1,000–2,000 ft and standard errors clustered by run. r(run) absorbs everything constant over a pass, including a level pass's height, so on level passes g is a clean range test.

| 5–10° below the horizon, ahead; g at 4,500–6,500 ft re 1,000–2,000 | 63–160 Hz | 200–630 Hz |
| --- | --- | --- |
| B407, all 115 steady runs (curved rays where profiled) | +1.1 ± 0.1 | +1.2 ± 0.2 |
| B407, straight rays | +1.0 ± 0.2 | +1.1 ± 0.2 |
| B407, descents only (85) | +1.2 ± 0.2 | +1.3 ± 0.2 |
| B407, level passes only (28) | −0.2 ± 0.3 | +1.2 ± 0.7 |
| B206L3, all 111 (straight rays) | +0.8 ± 0.2 | +1.2 ± 0.3 |
| B206L3, descents only (76) | +0.8 ± 0.3 | +1.3 ± 0.3 |
| B206L3, level passes only (29; 3,000–4,500 ft) | +0.2 ± 0.1 | +0.5 ± 0.2 |

Under 1,000 ft the 200–630 Hz term is −0.9 to −2 dB at 1–10°. Together these are the frames' near-horizon signature (short range over-predicted, long range under-predicted), found in the depropagated samples themselves.

**It is the descents, not the path.** Both aircraft's range term is carried by the descents. On level passes, where a far sample and a near one at the same angle differ only in path, it is absent at 63–160 Hz and about +0.5 dB at 200–630 Hz. Within a descent at one emission angle, range and height move together, and the far forward samples come from the pass's early part. That part is a different state: over the B407 descents' windows, height falls from 741 to 272 ft between the first and last thirds, speed by a median 3 kt, and 22 of 84 runs change glide path by more than 1°. This is the source-state reading the line analysis reached, from independent data.

**Consequence:** a rim fitted jointly with a range loss would carry the descents' state drift into every sphere as a propagation term. It is not worth building, and the 2,000 ft cap stays. What would address the residual is filing each sample under its own flight condition (μ and γ at emission time), not the run's mean, so a descent's early, faster and steeper part is not averaged into its sphere.

**Aft, a ray-model artifact beyond the cap:** with curved rays, aft samples past 4,500 ft depropagate +11 to +17 dB too loud at 5–10° (−9 dB at 10–20°). Under straight rays the same samples are within ±1.8 dB at 63–630 Hz. Depropagating through the ray model's too-deep shadows is the cause, which is a second reason to keep the cap with curved rays.

**Approach cards:** panam `bbf8224` refuses approach cards outright. The B206L3 databases are rebuilt without them (NICE-OPS `019a0d5`, September recipe; SEL rms 1.48 / 1.47 dB out of sample). The R66 rebuild is waiting on its data fetch.
