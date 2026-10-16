# Experiment: blending source spheres in dB, in energy, or picking the nearest

The literature does not give this number, so this note measures it. Across operating conditions, NICEOPS blends sphere levels **linearly in dB** with barycentric Delaunay weights over (μ, γ, CT). Here that blend is compared with three alternatives:

- a barycentric blend **in energy**, 10 log10 Σ w 10^(L/10);
- a **NORAH2-style blend**: inverse-distance weights over the same simplex, applied to energies;
- **nearest-sphere selection**, the sparse-data fallback that AAM/RNM and NORAH2 use.

Each is measured three ways:

- as a source-level difference per sphere direction;
- as a change in the SEL footprint;
- as accuracy, in a leave-one-out test against stored spheres.

All of it ran on the shipped databases and tracks. Nothing in the repository was modified.

## Bottom line

- **Energy blending is louder, but only a little where it matters.** It never comes out below the dB blend, as Jensen's inequality requires.
  - Over footprint directions of real trajectory states, the median gap is **0.05–0.15 dB**, the 95th percentile **0.5–1.2 dB**, and the maximum **3.6–6.1 dB**.
  - In the SEL footprint it adds **+0.25 dB on average**: 95th percentile **+0.5–0.6 dB**, maximum **+2.4–2.9 dB**. Fewer than 1% of grid cells move by more than 1 dB.
  - The gap is set entirely by the spread of the blended levels: about 0.03 dB at a 1–2 dB spread, 0.2 dB at 3–5 dB, 0.5 dB at 5–8 dB, and 1–2.7 dB at 8–20 dB.
- **Nearest-sphere selection is the larger departure, and it is worse.**
  - Per direction it differs from the dB blend by 2–3.4 dB at the 95th percentile and up to ±9–12 dB.
  - In SEL it moves 2–11% of footprint cells by more than 1 dB, with errors of −3.6 to +4.1 dB.
  - In leave-one-out it is the least accurate method on every database: rms 6–18% worse than the dB blend. The dB blend beats it on 78–88% of held-out spheres.
- **In leave-one-out, dB and energy are statistically indistinguishable.** On three of four databases dB has the lower rms by 1–3%. On B407, energy is lower (1.54 vs 1.76 dB), but only because it damps a single held-out sphere with a deep null.
  - In any case, the interpolation error (~1.1–1.8 dB rms over footprint directions) is roughly ten times the dB–energy difference.
  - The choice of blending domain is therefore a **systematic ~+0.25 dB bias question, not an accuracy question**.
- **NICEOPS cannot switch to an energy blend just by changing the formula.** The condition set is a degenerate prism, so the Delaunay triangulation is not unique.
  - The dB blend is invariant to which valid triangulation is chosen; the max difference measured was 0.000 dB.
  - The energy blend is not invariant: 95th percentile 0.10–0.19 dB, maximum 0.8–1.5 dB between two valid triangulations.
- **Correction to the task framing:** AAM/RNM blend **in dB**, not in energy. Up to four spheres are "linearly interpolated in decibels" ([aam_v3.md §2](aam_v3.md), TR p.33–34). Among the tools reviewed, only NORAH2 blends in energy ([norah2.md §3](norah2.md), D1.5d eqs. 3–10). NICEOPS and AAM share the dB convention; NORAH2 is the outlier.

## Method

### Data

| Database | Groups | Flight conditions (μ, γ) | Measured* | Tracks used |
| --- | --- | --- | --- | --- |
| `S-76D_M3.nod` | 1312 | 164 | 83 | 30 S-76D AC tracks (incl. demo 197421) |
| `AW139_M1.nod` | 808 | 101 | 52 | 41 AW139 tracks (per-direction stats only) |
| `B205_M1.nod` | 880 | 110 | 61 | 13 B205 tracks |
| `B407_M2.nod` | 1248 | 156 | 107 | none shipped; 20,000 random states (below) |

\*"Measured" excludes the padding conditions: μ = 0 (three per database) and γ = −24° / +35°.

Four things about these databases shape the result.

**The CT axis is synthetic.** Each flight condition has 8 copies at CT/CT_ref = 0, 1, 1.1, 1.2, 1.5, 2, 3 and 5.
- In every direction, each copy is the reference sphere plus exactly 20 log10(CT/CT_ref) (std across directions 0.000 dB).
- The CT = 0 copy equals the reference.
- So directivity varies only with (μ, γ). Blending along CT only interpolates a known scalar law.

**The hull is padded.**
- 78/81 (S-76D) and 46/49 (the others) of the padding conditions at γ = −24°/+35° are exact copies of a measured sphere at the same μ.
- A state outside the measured γ range therefore gets roughly constant extrapolation.
- No trajectory state fell outside the normalisation box, so no state was clamped.

**All four databases are pre-versioned (`database_version` absent).**
- The legacy hover-directivity repair was applied to the 24 μ = 0 spheres of each, as `Database::Database` does.

**Each database has 36 × 37 directions:** φ every 10° from −180° to 170°, and θ every 5° from 0° to 180°.

### States

- Each track was run through the existing build (`/Users/evg5332/Code/D/NICEOPS/build/niceops -S`, broadband, defaults).
- The query was (μ, γ_e = `fpa_e` (QSAM default), CT), as `process_trajectory` passes it.
- Per-direction statistics use every 20th state: 8,729 S-76D, 7,870 AW139 and 2,281 B205 states.
- B407 has no shipped track. Its 20,000 states are uniform over the measured (μ, γ) box, with CT/CT_ref drawn from 0.9 to 2.
- Across trajectory states, CT/CT_ref runs from 0.75 to 1.23 (1st to 99th percentile), so trajectories barely use the CT axis.

### Weights

- **Normalisation and triangulation.** The script copies `Database::Database` / `ConditionInterpolator::build`:
  - min–max normalisation of μ, γ and CT;
  - a boundary-duplicated vertex set extended by 0.1;
  - `scipy.spatial.Delaunay` over that set;
  - clamping;
  - the ≥ 1 single-vertex rule;
  - nearest vertex outside the triangulation.
- **Validation.** A numpy re-implementation of the broadband path takes these weights through to SEL:
  - sphere frame, φ/θ lookup, bilinear interpolation on the sphere grid, spreading, EAA, 10 dB-down window, reception-time weight.
  - It reproduces the NICEOPS LAE output to **≤ 2.4e-13 dB on all 43 S-76D and B205 tracks** at the default 120 × 120 grid.
  - The NICEOPS run of the demo case is byte-identical to `golden/197421_LAE.csv`.
  - So the dB-blend baseline here *is* NICEOPS, even though scipy's triangulation is not NICEOPS's. See the triangulation section for why that does not matter in dB.

### Variants, all at the sphere grid nodes

- **(a) dB:** Σ w L. This is NICEOPS.
- **(b) energy:** 10 log10 Σ w 10^(L/10), with the same weights.
- **(c) nearest:** the nearest stored condition in normalised (μ, γ, CT).
- **(c2) nearest FC + CT law:** the nearest (μ, γ) sphere plus the database's exact 20 log10(CT/CT_ref). This isolates the choice of directivity.
- **(d) NORAH2-style:** inverse normalised distance to the vertices of the same simplex, applied to energies.
  - This is D1.5d eq. 3 transplanted to NICEOPS's 3-D simplex and normalisation.
  - It does not use NORAH's 2-D (V, γ) space or its F_fc = 2 γ scaling.

Only the dBA blend differs between variants. EAA is always blended with the (a) weights.

### Direction sets

- **all:** all stored directions.
- **lower hemisphere:** |φ| ≤ 90°.
- **footprint:** |φ| ≤ 60°, 20° ≤ θ ≤ 160°.

## Results

### 1. Source level per direction (variant − dB blend, dB), footprint directions

| Database (states) | energy: median / mean / p95 / p99 / max | NORAH-IDW energy: median / p95 / max (min) | nearest: median / \|·\| p95 / range | nearest FC + CT: \|·\| p95 | spread of blended levels: median / p95 |
| --- | --- | --- | --- | --- | --- |
| S-76D (8,729 traj.) | 0.10 / 0.22 / 0.83 / 1.56 / 5.27 | 0.16 / 1.48 / 7.48 (−1.99) | 0.00 / 2.91 / −9.1…+10.4 | 2.91 | 2.70 / 7.30 |
| AW139 (7,870 traj.) | 0.05 / 0.12 / 0.50 / 1.01 / 3.55 | 0.07 / 0.75 / 6.89 (−2.13) | 0.04 / 1.98 / −8.7…+9.5 | 1.97 | 2.03 / 5.81 |
| B205 (2,281 traj.) | 0.15 / 0.32 / 1.21 / 2.37 / 6.14 | 0.15 / 1.55 / 8.94 (−2.45) | −0.15 / 3.44 / −9.7…+11.9 | 3.42 | 3.36 / 8.18 |
| B407 (20,000 random) | 0.14 / 0.23 / 0.74 / 1.38 / 48.5† | 0.15 / 1.12 / 53.4† (−2.97) | −0.03 / 2.65 / −66.9…+30.5† | 2.47 | 3.32 / 7.47 |

†B407's extremes come from one measured sphere (μ = 0.245, γ = −0.09°) that has a 26 dB minimum in some direction. Its blended spread then reaches 75 dB. This is a database artefact, not a property of the blend.

How the results vary by direction set and by state:

- Over all stored directions the energy gap is slightly larger: p95 0.60–1.34 dB, max 4.0–6.7 dB (B407 aside). The lower hemisphere matches "all" to within 0.05 dB.
- The worst footprint direction in each state has a median energy–dB gap of 1.0–1.9 dB (p95 2.4–4.2 dB). For nearest it is 2.3–3.5 dB (p95 7.3–9.2 dB), and for NORAH-IDW 1.2–2.2 dB (p95 3.3–4.9 dB).

The energy–dB gap depends only on the spread of the blended levels. On S-76D footprint directions, binned by spread (max − min over conditions with weight):

| spread (dB) | share of samples | median gap | p95 gap | max gap |
| --- | --- | --- | --- | --- |
| 0–1 | 11% | 0.009 | 0.021 | 0.029 |
| 1–2 | 24% | 0.032 | 0.080 | 0.114 |
| 2–3 | 21% | 0.091 | 0.183 | 0.257 |
| 3–5 | 25% | 0.221 | 0.486 | 0.701 |
| 5–8 | 16% | 0.531 | 1.140 | 1.759 |
| 8–12 | 3% | 1.11 | 2.18 | 3.65 |
| 12–20 | 0.3% | 2.65 | 4.35 | 5.27 |

The other databases agree within a few hundredths of a dB per bin, because the gap is a pure function of the levels and weights. As an analytic check, the worst-case gap for two spheres D dB apart is 0.03 / 0.26 / 0.71 / 2.69 / 5.62 / 9.18 dB at D = 1 / 3 / 5 / 10 / 15 / 20 dB. That maximum falls at a weight of 0.48–0.21 on the louder sphere.

### 2. SEL footprint (variant − dB blend, dB)

Each track uses the full default grid: 120 × 120 observers at 25 m spacing, ±1500 m.

| Set | Variant | mean | median | p5 | p95 | min | max | \|Δ\| > 1 dB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S-76D, 30 tracks (432,000 cells) | energy | +0.26 | +0.23 | +0.09 | +0.52 | −1.23 | +2.41 | 0.4% |
| | NORAH-IDW | +0.28 | +0.25 | +0.03 | +0.64 | −1.30 | +2.46 | 0.8% |
| | nearest | +0.19 | +0.19 | −0.73 | +1.08 | −3.19 | +4.04 | 9.8% |
| | nearest FC + CT | +0.25 | +0.24 | −0.67 | +1.14 | −3.59 | +4.09 | 10.9% |
| B205, 13 tracks (187,200 cells) | energy | +0.25 | +0.21 | +0.09 | +0.59 | −2.27 | +2.93 | 0.4% |
| | NORAH-IDW | +0.23 | +0.21 | +0.03 | +0.51 | −3.33 | +2.04 | 0.4% |
| | nearest | −0.07 | −0.04 | −0.74 | +0.49 | −2.36 | +4.00 | 1.9% |
| | nearest FC + CT | −0.05 | −0.02 | −0.73 | +0.51 | −2.79 | +4.01 | 1.8% |

Further detail from these runs:

- **Contoured region.** Restricting to cells within 20 dB of each track's peak (138k S-76D and 49k B205 cells) changes little:
  - energy +0.23 / +0.30 dB mean, p95 0.50 / 0.63 dB, max 1.78 / 2.12 dB;
  - nearest |Δ| p95 1.17 / 0.93 dB, with 8.4% / 3.7% of cells over 1 dB.
- **Per-track means.** The energy offset is +0.18 to +0.42 dB on every track. The nearest-sphere offset swings from −0.36 to +0.54 dB depending on the track: it is not a consistent bias but noise that depends on which sphere happens to be nearest.
- **Demo track 197421.**
  - energy: +0.26 dB mean, p95 +0.51, max +1.75;
  - nearest: p5 −1.54, p95 +1.12, max +3.12, with 22% of cells over 1 dB;
  - NORAH-IDW: +0.25 dB mean.
- **Why SEL can fall under an energy blend.** The energy blend never lowers an instantaneous level, yet SEL can drop (to −1.2 / −2.3 dB). Raising the peak raises the 10 dB-down threshold and shortens the window. The window's endpoints moved for 80% (S-76D) and 44% (B205) of observers. This is the same knife-edge that CLAUDE.md describes for decimation.
- **Instantaneous levels.** Over all state–observer pairs of the demo track, energy − dB has median +0.22 dB, p95 +0.85 dB, max +5.6 dB. Nearest − dB ranges from −8.9 to +9.2 dB. The SEL integral averages most of this out.

### 3. Leave-one-out accuracy

Method:
- Hold out each measured flight condition, together with all its CT copies and any padding copies of the same sphere.
- Rebuild the triangulation from the rest.
- Predict the held-out sphere at its own (μ, γ, CT_ref), then compare over footprint directions.

The errors also contain the held-out sphere's own measurement scatter, so they bound the interpolation error from above.

| Database (held out) | dB rms / bias / p95 \|e\| | energy rms / bias | NORAH-IDW rms / bias | nearest rms / bias / p95 \|e\| | energy beats dB | dB beats nearest |
| --- | --- | --- | --- | --- | --- | --- |
| S-76D (82) | **1.84** / −0.02 / 3.83 | 1.86 / +0.16 | 1.89 / +0.24 | 2.15 / +0.04 / 4.60 | 46% of spheres | 80% |
| AW139 (51) | **1.12** / 0.00 / 2.26 | 1.14 / +0.08 | 1.12 / +0.11 | 1.32 / +0.07 / 2.67 | 33% | 88% |
| B205 (60) | **1.58** / 0.00 / 3.37 | 1.63 / +0.16 | 1.58 / +0.20 | 1.67 / −0.11 / 3.53 | 42% | 78% |
| B407 (106) | 1.76 / −0.05 / 2.61 | 1.54 / +0.08 | **1.53** / +0.12 | 2.04 / −0.05 / 2.96 | 37% | 85% |

How to read the table:

- **dB is unbiased; the energy blends are not.** The dB blend has essentially zero bias. The energy blends are biased high by +0.08 to +0.24 dB, which is the same offset that appears in the footprint.
- **The rms values are close.** They agree to within 0.02–0.05 dB, except on B407. There, the deep-null sphere (max |e| ≈ 65 dB under every method) dominates the rms, and energy's damping of the null is what wins. By median per-sphere rms, dB is better on B407 too: 1.08 vs 1.13 dB.
- **Nearest** is unbiased, but it is the noisiest method on every database.

### 4. Triangulation dependence and the CT axis

**Triangulation dependence.**
- A second, equally valid triangulation was built with Qhull `QJ` instead of the default `Qt`. It shares the vertex set of the enclosing simplex in only 23–28% of trajectory states.
- The **dB blend is identical** under both: max |Δ| 0.000 dB on all states. The levels are additive, L = D(μ, γ) + g(CT), the conditions form a prism, and any barycentric weights reproduce both the (μ, γ) marginal and the CT coordinate, so a linear blend cannot tell triangulations apart.
- The **energy blend is not identical**: p95 0.19 / 0.10 / 0.18 dB and max 0.78 / 0.96 / 1.51 dB (S-76D / AW139 / B205).
- This is also why scipy's triangulation (12,330 simplices for S-76D, against NICEOPS's logged 10,049, including zero-volume slivers) still reproduces NICEOPS to round-off.

**The CT axis on its own.** The stored copies follow 20 log10 CT exactly. In LF, a dB blend lies below that law (the chord of a concave function), and an energy blend 10 log10(Σ w CT²) lies above it by the same amount:

| interval | LF 1.0–1.1 | 1.2–1.5 | 1.5–2 | 2–3 | 3–5 |
| --- | --- | --- | --- | --- | --- |
| largest deviation from the 20 log10 CT law, either blend | 0.01 dB | 0.05 | 0.09 | 0.18 | 0.28 |

The exact law sits midway between the two blends. At trajectory load factors this contributes ≤ 0.02 dB (S-76D) and ≤ 0.11 dB (AW139, B205) to the energy–dB gap. Effectively all of the gap comes from directivity differences between (μ, γ) spheres.

## What this says about NICEOPS

1. **Keep dB, but say why.** NICEOPS's dB blend agrees with AAM's convention and is unbiased in leave-one-out. It is also the only blend that is well defined on the current degenerate Delaunay set. Nowhere does the code or its documentation justify the choice, and the comparison with NORAH2 deserves a sentence there:
   - the domain choice costs about 0.25 dB of systematic offset in SEL against an energy-blending tool such as NORAH2;
   - up to ~2–3 dB in individual cells;
   - about 1 dB in single directions at the ~5 dB spreads that are common between neighbouring flight conditions.
2. **An energy option would need a different scheme.** If an energy blend (for NORAH2 parity) were ever wanted, a formula swap is not enough, because the result would depend on Qhull/Bowyer–Watson tie-breaking by up to ~1.5 dB. It would need a tensor-product scheme: a 2-D (μ, γ) triangulation crossed with linear interpolation in CT. Since the CT copies are the exact 20 log10 law, CT could instead be applied analytically, as the `fixed_load_factor` path already does.
3. **Nearest-sphere fallbacks cost more than the blending domain does.** NICEOPS falls back to nearest only for databases with fewer than dim + 1 conditions, or outside the triangulation. That is the right call: nearest selection adds ±1–4 dB of track-dependent SEL noise and is the least accurate method in leave-one-out.
4. **Interpolation error dominates both.** Leave-one-out gives ~1.1–1.8 dB rms over footprint directions on every database. Database density and measurement scatter matter about ten times more than dB versus energy.
5. **Incidental observations.** These were noticed while reproducing the path and were not part of the question.
   - **Padding means extrapolation.** A third to a half of each database's flight conditions are padding copies at γ = −24° / +35°, so the blend outside the measured γ range (S-76D: −8.7° to +14.4°) is constant extrapolation of the edge spheres. That is a builder (panam) property worth stating in the review.
   - **φ folds onto the lower hemisphere.** `phi_degrees()` computes atan(across / |up|) (`src/sphere_database.cpp`, `phi_degrees`), so the lookup only ever reaches |φ| ≤ 90°. An observer above the tip-path plane reads the mirrored lower-hemisphere direction, never the repaired upper surface. That is harmless for ground observers below an elevated source, but it may matter for climbs and for receivers near or above the source. It should be checked separately; this experiment does not test it.

## Caveats

- **Sphere-grid nodes only.** Part 1 compares blends at the sphere-grid nodes. In NICEOPS the blend happens at the four corners of each sphere cell, which are then interpolated bilinearly in dB. The energy variants here change only the condition blend, not the bilinear step. NORAH2 also interpolates angles in energy, which this experiment does not include.
- **NORAH2 analogue.** The "NORAH-IDW" variant uses NORAH's weighting rule on NICEOPS's 3-D simplex and normalisation. It is not a reimplementation of NORAH2's 2-D (V, γ) triangulation with F_fc = 2.
- **What nearest means here.** "Nearest" is Euclidean distance in NICEOPS's normalised (μ, γ, CT) box, or in (μ, γ) for variant c2. AAM selects by speed and FPA with dFPA/dSPD gates and a higher-speed tie rule, so its actual choices can differ.
- **Footprint scope.** Footprints are broadband only (`--no-spectral`) and use dBA spheres. They were run for S-76D and B205. AW139 has per-direction statistics but no footprint run, and B407 has random states only.
- **What leave-one-out can and cannot show.** It removes a sphere and its copies, so it tests a sparser database than the real one. Its errors also include each held-out sphere's own measurement noise. It ranks the methods; it does not measure absolute accuracy at full density.
- **Triangulations.** scipy's triangulation differs from NICEOPS's Bowyer–Watson. The dB results are exact regardless, as shown above. The energy and IDW numbers depend on which triangulation was used, by the amounts given in section 4.
- **Legacy databases.** All four databases are legacy: no `database_version`, and the hover repair applied. Their CT copies are synthetic, so "CT interpolation" here tests only the 20 log10 law, not measured thrust dependence.

## Reproduction

All scripts are in `/Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend/` and run with `/Users/evg5332/Code/Python/niceops-validation/.venv/bin/python` (numpy 2.5, scipy 1.18, netCDF4 1.7).

| Script | What it does |
| --- | --- |
| `nod.py` | `.nod` loader with the legacy hover repair |
| `inspect_db.py` | caches the grids to `*.npz` |
| `blend.py` | NICEOPS weights and the alternative blends |
| `states.sh`, `batch.sh` | NICEOPS `-S` / LAE runs, written to `states/` and `fp/` |
| `direction_stats.py` | part 1 (`direction_stats.out`, `.json`) |
| `footprint.py`, `aggregate.py` | part 2 (`fp/*/summary.txt`, `agg_S-76D.out`, `agg_B205.out`) |
| `loo.py` | part 3 (`loo.out`, `loo.json`) |
| `extra_checks.py` | part 4 (`extra_checks.out`) |
| `inspect_ct.py`, `pads.py`, `check_b407.py` | database-structure checks |
