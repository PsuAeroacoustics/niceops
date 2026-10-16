# Independent checks: NORAH2 roll sign, and AAM v3 interpolation/Doppler/SEL details

This note re-checks two claims in the existing notes, working only from local documents and data. Nothing in the NORAH2 or AAM packages was executed. All analysis was done with Python scripts that read the shipped text files. The scripts are in `/Users/evg5332/.claude/jobs/a090ec0c/tmp/rollcheck/`:
- `fit_roll.py`, `inout.py`, `clamp.py`: NORAH frame fits;
- `phigrad.py`: sensitivity of hemisphere level to φ;
- `pages.py`: AAM page lookup.

Package roots:
- NORAH2: `/Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/`, written `NORAH2/` below.
- AAM: `.../Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf`, written "TR". Page numbers are the printed numbers (PDF page − 16).

---

## A. NORAH2 roll sign

### Verdict
**Confirmed: the shipped departure banks the hemisphere toward the outside of both turns.**

I tried every refutation route below and none survives. The defect is in the multi-event module's track-to-trajectory step, not in SingleEvt. SingleEvt applies ROLL consistently. The multi-event module, however, gives each **track segment** a fixed ROLL sign:
- an `R` segment always gets −atan(V²/gR);
- an `L` segment always gets +atan(V²/gR).

The sign does not depend on the direction the segment is flown. A segment flown one way is a right turn and flown the other way is a left turn, yet the approach and the departure carry the same sign on it. **Whatever the sign convention, one of the two must therefore be wrong.** The geometry shows which: the departure banks outward and the approach inward.

The earlier note's "unexplained" 0.15° rms / 1.1° max φ residual on the approach is also explained (see the clamping finding below): it is not a frame error.

### Evidence 1: the trajectory rows themselves (no assumption about SingleEvt)
Track NE (`NORAH2/ARP/Tracks.csv`):

| idSeg | type | details |
| --- | --- | --- |
| 0 | O | initial heading 45° |
| 1 | S | 2000 m, slip 90° |
| 2 | R | radius 2000 m, 45° |
| 3 | S | 4000 m |
| 4 | L | radius 1500 m, 90° |
| 5 | S | 4000 m |

Headings are clockwise from Y_M (north) (Manual §4 p.16, §5.1.5 p.23). In every `.inp`, the HEADING column equals the track bearing atan2(ΔX, ΔY) computed from consecutive positions (difference 0.0–0.1° in turns). The exception is the slip-90° leg, where HEADING − track = +90°. That matches the manual's "rotation angle = horizontal angle + slip".

Turn sense below comes from the sign of the cross product of velocity and centripetal acceleration, both computed from positions. It does not use the R/L labels or the ROLL column.

| File | Segment | Heading change | Physical turn | ROLL | atan(V²/gR) |
| --- | --- | --- | --- | --- | --- |
| `Case2/SE/R22_H1_DEP_STD1_NE.inp` | R (t = 72.5–127.5 s) | 45° → 90° | **right** | **−2.34°** | 2.34° (55 kt, R = 2000 m) |
| same | L (t = 269–352.5 s) | 90° → 0° | **left** | **+3.12°** | 3.12° (55 kt, R = 1500 m) |
| `Case2/SE/R22_H1_APP_STD1_NE.inp` (same track, reversed) | L reversed (t = 428–511 s) | 180° → 270° | **right** | **+3.12°** | |
| same | R reversed (t = 652.5–707 s) | 270° → 225° | **left** | **−2.34°** (−2.17° at 53 kt) | |

`A600_H1_APP_STD1_NE`, `R22_H1_APP_STD2_NE`, `A600_H1_APP_STD2_NE` and `CH7_H1_APP_STD1_NE` repeat the approach pattern: the reversed-L right turn has positive ROLL and the reversed-R left turn negative.

**The same physical turn sense carries opposite ROLL signs in the departure and the approach.** No single sign convention can make both right. This holds before asking how SingleEvt uses ROLL.

The documented convention also points at the departure. D1.5d §A.3.3 eqs. 18–20 (`e132005.pdf`, p.16–17) define:
- Θ = atan2(ΔX, ΔY), clockwise from north;
- K = ΔΘ/ΔS;
- Φ = atan(K V_g²/g).

Under these, **Φ > 0 in a right (clockwise) turn**. The approach agrees with eq. 20 and the departure has the opposite sign. Neither the Manual nor D1.5d states the ROLL sign in words. Eq. 20 is the only written convention, and the departure violates it.

### Evidence 2: how SingleEvt applies ROLL (fit to `.his` emission angles)
Method:
- The `.his` rows give the aircraft position at emission (XAC_M, YAC_M, ZAC_M, tem) and the observer (XMICM, YMICM, ZMIC + HMIC).
- I rebuilt the hemisphere frame: x along HEADING in the horizontal plane, y starboard, z down, then rolled about x.
- θ = acos(d·x) and φ = atan2(d·y, d·z) were compared with THETAOBS and PHIOBSAC.

Results, over rolled samples only:

| File | Hypothesis | rms Δφ | max Δφ |
| --- | --- | --- | --- |
| Case2 R22 DEP | +ROLL = starboard wing down | **0.003°** | 0.005° |
| | +ROLL = port wing down | 5.67° | 6.25° |
| | roll ignored | 2.83° | 3.13° (= the roll) |
| Case4 R22 DEP | +ROLL = starboard down | **0.003°** | 0.005° |
| Case2 R22 APP / A600 APP | +ROLL = starboard down | 0.35°* | 1.12°* |
| | +ROLL = port down | 5.52° | 6.25° |
| Case4 R22 APP STD2 | +ROLL = starboard down | 0.20°* | 0.75°* |

\*All of this residual is the clamping described below. Excluding the clamped rows, the residual is ≤ 0.006° in every file.

Across all rows, θ fits to 0.003° rms. On unrolled rows, φ fits to 0.003–0.004° rms (R22 DEP n = 1343, R22 APP n = 1286). That pins the φ sign convention independently of roll: positive φ is starboard, as D1.5d eq. 11 says.

So **SingleEvt applies +ROLL as starboard-wing-down (right bank)**, the aviation convention, and is consistent with eq. 20. Note that each `.his` file holds one microphone, at (−300, −400, −3.8) m.

### Evidence 3: a convention-free version of the same test
`inout.py` drops ROLL's sign altogether:
- It takes only |ROLL|.
- It takes the turn centre from the trajectory's own curvature, as the direction of the finite-difference acceleration.
- It tilts the rotor "up" axis either toward the centre (the physical coordinated-turn bank) or away from it.

Results, excluding polar θ < 20° and clamped rows:

| File | Bank toward centre | Bank away from centre |
| --- | --- | --- |
| R22 DEP (Case2) | 5.66° rms | **0.026° rms** (max 0.10°) |
| R22 DEP (Case4) | 5.69° | **0.027°** |
| R22 APP STD1 (Case2) | **0.032°** | 5.28° |
| R22 APP STD2 (Case4) | **0.027°** | 4.85° |

The departure's emission angles are those of a hemisphere **banked away from the turn centre, in both its right and its left turn**. The approaches bank toward it. This uses no sign convention for ROLL, heading or φ, only geometry, so it is the most direct statement of the bug.

### Refutation attempts and why each fails
1. **"ROLL's sign convention might be positive = port down."**
   - SingleEvt's own outputs reject it (5.7° rms against 0.003°).
   - Evidence 3 does not use the convention at all.
2. **"HEADING might be counter-clockwise, or the turn direction misread."**
   - HEADING equals atan2(ΔX, ΔY) from the positions, and turn sense was computed from positions.
   - Z is up: the departure climbs from 0 to 200 m. X is east: heading 90° = +X, per Manual §5.1.5. With Z up, (X, Y, Z) is right-handed, so clockwise-from-above really is a right turn.
3. **"Maybe the approach is flown backwards, keeping the departure heading."** No. Approach HEADING follows the reversed track (180° → 270° → 225°). The only slip is the final 90° leg, the mirror of the departure's first leg.
4. **"The R22 DEP turns are in FOV mode at 55 kt; maybe a different rule applies to level-flight hemispheres."**
   - D1.5d §A.3.4 says level-flight hemispheres are "tilted using the bank angle" in turns, with no sign exception.
   - The approach turns are also FOV/level (mode 1 at 300 m) and bank inward.
5. **"The manual warns that tracks are defined 'as seen in a departure'"** (Manual §5.1.5 p.24). That is the point: the departure is the case the definition is meant for, and it is the one that comes out wrong. The approach is right only because reversing the segment order flips the turn sense while the stored sign is kept. Two sign errors cancel there.
6. **"Maybe outward bank is a physical choice, e.g., a skidding turn."** The magnitude is exactly atan(V²/gR), the coordinated-turn bank. That only makes sense with lift tilted toward the centre.
7. **"Maybe the radar-trajectory examples show the intended convention."** `ARP/Radar_Trajectories.csv` has Roll = 0 on every row (348 rows), so it cannot adjudicate. There is no other departure with turns in the package: Case2_Operations has only R22 DEP STD1.

### How large the error is
- The departure hemisphere is looked up at φ off by 2|ROLL|: 4.7° in the right turn and 6.2° in the left.
- Sensitivity, A-weighted (`phigrad.py`, from 10°-neighbour differences in filled bins):

  | Hemisphere | median \|dL_A/dφ\| | 90th percentile |
  | --- | --- | --- |
  | `R22_Flyover_55kts_0deg.hem` | 0.04 dB/° | 0.14 dB/° |
  | `R22_Takeoff_54kts_7deg.hem` | 0.06 dB/° | 0.18 dB/° |

- For a 5.6° φ error this is about 0.2–0.3 dB typical and 0.8–1.0 dB at the 90th percentile, along the turns only.
- The error is small for these gentle demo turns (bank ≤ 3°). It scales with bank angle, so a 30° bank (60° φ error) would be gross.
- This is a sensitivity estimate. NORAH was not rerun.

### New finding: SingleEvt clamps φ at ±90° for directions above the hemisphere
- Wherever the geometric φ in the rolled frame exceeds 90°, which is a direction above the hemisphere's equator, `.his` reports PHIOBSAC = ±90.00 exactly, keeps the sign and leaves THETAOBS unchanged. This held in every case: 58/58 rows in R22 APP STD1 and CH7 APP, 50/50 in R22 APP STD2 and A600 APP STD2 (`clamp.py`).
- So upward directions are mapped onto the equator at the same θ. They are neither discarded nor mirrored.
- In these examples it happens only when the bank lifts a distant, near-grazing observer above the equator.
- **This accounts for the whole residual that norah2.md §4 left unexplained.**
- It is NORAH's answer to the "upper hemisphere" question, and differs from NICEOPS, which uses a full sphere with 180 − φ completion.

### Relevance to NICEOPS
- NICEOPS takes the sphere frame from the TPP normal (`SphereFrame`: nose, normal, starboard; see CLAUDE.md, hot-path notes). That normal is the load-factor direction, from the trajectory's own acceleration.
- So bank follows the curvature physically, and NICEOPS cannot get a turn-sense sign wrong this way. There is no separately stored ROLL whose sign could disagree with the path.
- The NORAH bug is a warning for any interface where NICEOPS might *accept* an externally supplied roll or bank column: state the sign convention in words, and cross-check it against path curvature on load.

### Corrections to norah2.md
- §4's "possible roll-sign problem (needs confirmation)" is **confirmed**. The defect is located in the multi-event track-to-trajectory conversion: R → negative, L → positive, irrespective of the direction flown. SingleEvt is consistent.
- §4's "0.15° rms (1.1° max) φ residual, unexplained" is **explained by φ clamping at ±90°**.

---

## B. AAM v3 details

### B1. Interpolation across flight conditions: dB or energy?
**AAM v3.1 (2021 TR): linear in dB, stated explicitly three times.**
- "The outcome of the search routine is the identification of up to four sound spheres, which are linearly interpolated in decibels based on speed and flight path angle" (TR p.33 §2.5.1.1).
- "The interpolation procedure remains, namely a linear interpolation on decibels ... first two sets of two sources interpolated to the desired airspeed and then these resultant spheres are linearly interpolated to the desired flight path angle" (TR p.34).
- Fixed-wing: spheres "are linearly interpolated if so indicated, in decibels, based on power setting and thrust vector angle" (TR p.36), and "the interpolation procedures remain a linear interpolation on decibels" (TR p.37).
- The FRAME/load-factor path (up to 8 spheres, TR p.41–42) states no domain.

**The RNM 2002 paper says energy.** I fetched Conner & Page, Heli Japan 2002 (NTRS 20030003730) and extracted its text to `/Users/evg5332/.claude/jobs/a090ec0c/tmp/rnm2002.txt`. Its two-column layout is garbled in extraction, but the sentence reads: source characteristics "extracted from the hemisphere files [by a linear interpolation] in the energy domain at the required [airspeed] ... followed by a second linear interpolation in the energy domain on the vehicle flight path angle".

The AAM TR carries what is plainly the same sentence (TR p.31: "...extracted from the sphere files using a linear interpolation at the required airspeed, followed by a second linear interpolation on the vehicle flight path angle"), with **"in the energy domain" removed**. On p.33–34 it then says "decibels" and "remains", which implies this has always been the method.

**Resolution:**
- For AAM v3, the documented method is **linear in dB**, so aam_v3.md is right.
- other_tools.md is right that the 2002 RNM paper says energy. Its line 30, "RNM / AAM ... interpolated linearly in energy", should not be applied to AAM v3.
- What the code did in 2002 or does now cannot be settled from documents. Either RNM changed domain at some version, or one document is wrong. The AAM TR is the later document and the more specific one (it repeats "decibels" and covers fixed-wing too), so take dB for AAM v3.
- For NICEOPS comparisons: NICEOPS blends in dB (barycentric over (μ, γ, CT)), which matches AAM v3 as documented.

### B2. Levels between sphere grid nodes (θ, φ)
**Not documented anywhere in the TR.** An exhaustive search covered:
- all 173 occurrences of "interpolat";
- "bilinear", "nearest", "closest", "neighbo", "grid node", "mesh point", "resolution", "directivity", "look up", "emission angle", "sphere location".

Every within-sphere interpolation statement concerns *building* spheres in ART:
- ZGRID LaPlace/spline gridding of scattered control points (TR p.146; GRIDOPTIONS p.154);
- pole averaging over φ and linear fill of the adjacent θ columns (TR p.146);
- energy averaging inside NPOLEGRIDAREAS (TR p.147);
- ART3's former "linear interpolation" (TR p.146).

None of it concerns how AAM *reads* a sphere during propagation.

What the documents do show:
- AAM evaluates at **continuous** angles. The `.TH` "Theta"/"Phi" columns are the "Sound sphere location for source level (degrees in sphere coordinate system)" (TR p.84–86). The tutorial's Figure 5-29 lists Theta = 3.01°, Phi = 33.42° on a 5° grid (TR p.203). Whether the level is then taken from the nearest node or interpolated (and if so, in dB or energy) is not shown: that example uses the 200 dB OMNI sphere, where every option gives 200.
- A node set to **1E35** is "not propagated" (TR p.145). How an interpolation stencil that touches a 1E35 node behaves is unstated.
- The only "linear (in energy)" interpolation in the TR is the 3D Visualizer's interpolation in *time* (TR p.176). It is unrelated.
- Gap: unchanged from aam_v3.md §4. It needs the code, or the RNM manual (Page, Wilmer & Plotkin 2008, WR 08-04).

### B3. Sign of u in f′ = f/(1 + u/c)
- **u is not defined anywhere in the TR.** Eq. 2 (TR p.11 §2.2.1) is followed directly by "In general, unless data will always be applied...".
- The sign follows from the sentence before it: "The net result of the Doppler shift is to increase frequencies for sources that are approaching a receiver, and to decrease frequencies for sources that are moving away".
- f′ > f needs 1 + u/c < 1, so **u is negative when approaching and positive when receding**. u is therefore the rate of change of source–receiver distance: u = −V·r̂, with r̂ the unit vector from source to receiver. Equivalently, u = −V cos θ_rel, where θ_rel is the angle between velocity and the ray.
- This is an inference from the prose, not a definition.
- Other Doppler facts in the TR:
  - The `.TIG` output's "Doppler" column is the "Doppler frequency multiplier; non zero if the .NC files are reflective of a stationary source" (TR p.91).
  - Moving-frame spheres (DOPPLERSTATE 1) are corrected by "a difference based on Equation (2) applied to the flight test and the modeled airspeed" (TR p.11). Whether that is applied per sphere before blending or to the blend is still unstated.
  - In the tutorial `.TH` (DOPPLERSTATE 1), FrqDop = FreqNC exactly (TR p.203).

### B4. SEL quadrature
- Emission points sit on the trajectory at the segmentation time step, 2 s by default (TIMESPACING, TR p.88).
- Sound "is 'binned' at the receiver locations based on the arrival time. Only after the propagation has been completed, are the receiver sounds interpolated onto a uniform time mesh and added" (TR p.10).
- SEL is defined only conceptually: the 1-second constant level with the same energy (TR p.50).
- Integration covers "the entire noise event" (TR p.51, contrasting it with SENEL). DBDOWN10 "forces all integrated metric calculations to only include points within 10 dB of the maximum level" (TR p.64).
  - That wording reads as a **threshold on samples**, not a contiguous window. NICEOPS's `peak_range` is a contiguous first-crossing window. The two differ whenever the history dips below and recovers.
- **Not stated anywhere:**
  - the quadrature rule (rectangle or trapezoid);
  - the step of the uniform receiver mesh;
  - whether the arrival-time interpolation is linear in dB or in energy.

  PEARLS (macro/micro steps, TR p.72) changes the sampling but not the rule. This confirms the gap recorded in aam_v3.md.

### B5. "Duration adjustment"
- It appears **once**, in Table 2-8 (TR p.34): "FPA selections occur before SPD selections. (A duration adjustment will be made by AAM in the course of its normal calculation procedure)."
- There is no formula, keyword or other mention. "duration" elsewhere refers to run-up, blasting, static events, EPNL's duration correction and PEARLS.
- **Inference** (labelled as such): the parenthesis is explaining why speed is the *soft* criterion. Selecting a sphere at the wrong speed does not get a 10 log(V_ref/V)-style exposure correction of the kind integrated-database models (NOISEMAP/Omega, INM) need. In AAM, exposure duration comes out of the time-stepped simulation along the actual trajectory at the actual speed, "in the course of its normal calculation procedure".
- So no separate adjustment term probably exists. The sphere's speed affects only its spectrum and directivity.
- This fits the TR's own explanation of why scaling an integrated metric is "no longer possible" for 3-D spheres (TR p.37, §2.5.2).

### B6. Upper-hemisphere directions with a lower-hemisphere sphere
**Not documented.** What the TR does say:
- A file with PHI −90…+90 "represent[s] hemispheres, with only the lower portion populated" (TR p.214). Table 4-37 lists PHI as −90° to +90° (TR p.163–164). Full spheres (−180…180) are also legal, e.g. `OMNI_200.NC` (TR p.204–205).
- ART regions with no data are 1E35, and AAM "will not propagate the source point" there (TR p.145).
- WRAPAROUND "takes a hemisphere (lower half) ... and mirror images the values on the top half of the sphere ... so that a fully defined sphere will be populated with reasonable source levels instead of 1E35 (no data) points" (TR p.161 §4.1.2.38; keyword summary table TR p.150).
  - "Mirror image" onto the top half is most naturally the reflection through the rotor/wing plane, φ → 180° − φ, which keeps port and starboard. The TR does not give the formula.
  - That is exactly the ambiguity behind NICEOPS's legacy −φ hover bug.
- EXTEND45 / EXTENDX copy the values at φ = 45° (or a chosen φ) out to the upper extent of the SPHERE domain (TR p.153).

**Inference:** a direction outside a hemisphere file's φ extent is treated like 1E35, i.e. not propagated. That is why ART offers WRAPAROUND and EXTEND to fill the top before AAM sees it. The TR never says so for a direction outside the stored axis range, as opposed to a 1E35 node inside it. Clamping, as NORAH does, is also possible.

Contrast:
- NORAH2 clamps to the equator (A, above).
- NICEOPS completes the sphere at load with 180 − φ and repairs legacy files that used −φ.

### Summary of corrections to the other notes
- **aam_v3.md:**
  - "linearly in dB" stands.
  - Add that the RNM 2002 paper said "energy domain" and the AAM TR silently dropped those words.
  - u's sign is inferable from prose (u > 0 receding), not defined.
  - DBDOWN10's wording implies a sample threshold, not a contiguous window.
  - "Duration adjustment" is most plausibly a remark that time-stepping handles duration, not an extra term.
- **other_tools.md:** line 30 ("RNM / AAM ... interpolated linearly in energy") is right for the 2002 RNM paper and wrong for the AAM v3 TR. Split it.
- **norah2.md:** the roll bug is confirmed and located (multi-event R/L → fixed sign). The φ residual is explained by clamping at ±90°.
