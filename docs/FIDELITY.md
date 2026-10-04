# Fidelity status and validation

Updated 2026-10-04. This page distinguishes regression agreement, numerical
accuracy against a reference solver, and accuracy against measurements.
Passing a golden comparison establishes stable behavior; it does not establish
physical accuracy. Historical grazing-incidence reports contain successive
hypotheses and rescoring, so use their latest addenda when interpreting them.

## Current changes

| Change | What it establishes | Remaining limit |
| --- | --- | --- |
| Measured B407 hover candidate | High-hover correction, with absolute held-out scores and fresh forward predictions | Experimental; measurements cover two elevation rings, and the low ring is in ground effect |
| Emission-time sample export and regridding | Every sample carries its own ground speed, advance ratio and flight-path angle; narrow condition bins replace run-wide labels in experimental spheres | The pilot uses two descent runs; a complete training rebuild and held-out flight validation remain necessary |
| Propagation-table midpoint validation | Between-node predictions are checked against direct solves with a specified error budget; unsupported points are counted | Cell midpoints alone cannot bound all interpolation error, and a PE reference is not a measurement |
| Source-axis diagnostics | Each run reports how many trajectory states were clamped, and how many received analytic load-factor scaling | Axis ranges include synthetic extensions; this is not a measured angular-coverage or convex-hull test |
| Common frame scoring | Models are compared on the same measurements, with geometry bins and explicit coverage counts | SNR and flyover-floor gates still require sensitivity checks for a new study |

The shipped databases and propagation defaults have not been replaced. The
experimental source artifacts are under `fidelity_runs/`, which is ignored by
git. The tools and their validation checks are tracked.

## Hover: measured improvement and a rejected extension

The source is the B407 slow level pass with its Doppler frequency shift removed,
corrected against the day-283 flying hovers. Ground/flight idle cards are excluded.
The Fourier order is fixed at 3 before scoring. Training fits include
sum-to-zero offsets per heading card. Validation reports both absolute errors
and errors after removing the held-out offset; the absolute score retains the
wind-related variation a prediction cannot know.

For 315–1,250 Hz, the fixed-order, leave-one-heading-out bin RMS is:

| Ring | Uncorrected de-Dopplerized baseline | Corrected | Decision |
| --- | --- | --- | --- |
| High, out of ground effect | 7.38 dB | 3.51 dB | Retain the high-elevation correction as a candidate |
| Low, in ground effect | 6.41 dB | 6.90 dB | Do not apply this low-ring Fourier correction by default |

Leaving one repeat pass out instead gives 2.87 dB high and 6.14 dB low.
These are held-out **binned correction** scores using cached forward predictions,
not fresh predictions from separately built databases for every fold. Fresh
predictions from the all-training candidate are a separate, in-sample check.

Fresh forward runs used the same current GCC binary for the shipped and
candidate databases: 16 flying hovers and 25 held-out day-286 approaches.
For the high-hover frames at 17.5–25°, the 315–1,250 Hz energy-mean error changes
from **−6.81 to −1.83 dB**. This is the all-training check, not the held-out RMS
above. On the approaches, pooled band-group energy-mean errors move by at most
0.013 dB. Per-run geometry bins can move more: the largest movement is 1.14 dB
on run 286267 at 2–5° and 50–160 Hz, reducing error from +3.36 to +2.22 dB.
Thus the pooled result should not be read as every frame being unchanged.

Neither comparison introduces candidate prediction gaps. The hover set has
280,683 measurement frames passing the SNR gate: 1,245 lack a shipped prediction,
158 lack a candidate prediction, and the common score excludes those 1,245.
The approach common score excludes 91 existing missing predictions.

Local review artifacts:

- `fidelity_runs/hover_high_only_coverage_2026-10-04/`: final experimental
  database, input hashes, correction arrays, fresh frames and score CSVs.
- `fidelity_runs/shipped_current/`: the matching current-binary baseline frames.
- `fidelity_runs/emission_descents_final_2026-10-04/`: aligned descent samples,
  exact recipe controls, and the upstream emission-time patch.
- `fidelity_runs/emission_state_spheres_final_2026-10-04/`: nine experimental
  condition spheres and their coverage counts.

The default candidate applies the high-ring correction fully at 17.5–25° below
the horizon, tapering to zero over 5° beyond each edge. It preserves the
synthesis elsewhere. This is conservative about unmeasured directions and can
leave a localized elevation ridge. `--extent two_ring` retains the broader
experimental construction for comparison, but its low-ring generalization
does not support making it the default.

De-Dopplerization left 45 previously supported source cells empty at the
recorded spectrum's edge. The candidate explicitly restores the original
flight-source values at those cells before hover synthesis and records that
fallback in its manifest. A missing measurement must not silently become a
new confident prediction or disappear from the score.

Outstanding physical questions include the 63 Hz disagreement between the
forward-correction and direct-depropagation routes, separating wind effects
from source directivity, and the unmeasured elevations between the two rings.
No numerical uncertainty curve is supplied for these unmeasured directions.

## Emission-time flight conditions

The exporter preserves the existing recipe, including its 2,000 ft path cap,
and rejects approach cards and runs outside the original training manifest.
It adds emission times to panam's scattered samples in memory, then obtains
speed and flight-path angle from the tracking record at those times. Advance
ratio uses the ground-speed and rotor-tip-speed convention of the existing
database builder. Physical microphone identifiers remain attached.

The pilot controls reproduce the original spheres with **0 dB maximum
difference**, including identical finite-data masks:

| Descent | Run mean speed | Sample speed span | Run mean path angle | Sample angle span |
| --- | --- | --- | --- | --- |
| 283146 | 80.22 kt | 76.69–84.66 kt | −8.12° | −9.60° to −6.54° |
| 283148 | 64.46 kt | 61.46–68.20 kt | −5.92° | −7.62° to −3.92° |

The regridder assigns each sample to one 5 kt × 2° condition bin, grids in
linear energy, and labels the resulting sphere by its samples' mean condition.
It requires at least 200 samples and three physical microphones per retained
condition, and preserves missing angular nodes. Nine pilot spheres were built;
three have angular gaps. These are training experiments, not evidence of
improved held-out exposure predictions. They demonstrate why a complete rebuild
must assess coverage as well as fit residuals.

## Reproducing source work

These optional tools require NumPy and the local `niceops-validation` harness
with its panam dependencies, cached frame measurements and hydrated recordings.
Use its Python environment. Set `PYTHONDONTWRITEBYTECODE=1` and put
`MPLCONFIGDIR` and `NUMBA_CACHE_DIR` in writable scratch directories when
running in a restricted environment. All output directories must be fresh.

```sh
python -B tools/build_hover_candidate.py \
  --harness /path/to/niceops-validation \
  --out fidelity_runs/hover_candidate

python -B tools/export_emission_samples.py \
  --harness /path/to/niceops-validation \
  --out fidelity_runs/emission_samples --runs 283146 283148

python -B tools/regrid_emission_samples.py \
  --harness /path/to/niceops-validation \
  --samples fidelity_runs/emission_samples \
  --out fidelity_runs/condition_spheres
```

The hover candidate writes a self-contained `.nod`, correction arrays,
held-out scores, a fallback source, and a manifest with input hashes. Its flight
sphere links are input references; the resulting `.nod` contains the spectra.
The exporter writes state-labelled samples, control comparisons, and a guarded
patch that adds emission time to upstream scattered output. It does not modify
the external checkout or the recordings.

For fresh forward predictions, point the harness at this build with
`NICEOPS_SRC` and `NICEOPS_BIN`, then run its `propagation/band_frames.py` with
an absolute output path, `--database` pointing to the candidate, `--config
station_plate --day 283 --tables hover` for hovers, and `--config stratified_plate
--day 286` for approaches. Run the shipped database through the same current
binary and configurations. Historical cached geometries can differ from the
current code, and the scorer deliberately rejects a mismatched comparison.

```sh
python -B tools/score_fidelity_frames.py \
  --baseline fidelity_runs/shipped/hover \
  --candidate fidelity_runs/hover_candidate/frames \
  --output fidelity_runs/hover_candidate/hover_scores.csv --hover
```

The scorer writes per-run and pooled energy-mean errors by frequency, elevation
and range, plus a separate coverage CSV. Hovers use the ambient SNR gate only;
approaches also use the recording-floor gate. It never selects on a predicted
shadow, the measured peak, or the model's residual. New candidate coverage gaps
return exit code 2; existing gaps remain visible in the report.

## Propagation-table error budgets

Append validation options to a normal table-generation command:

```sh
./build/niceops_propagation_table --method pe \
  --atmosphere atmospheres/example_inversion.csv --units ft \
  --bands 125 250 --source_heights 100 110 --receiver_heights 3.937 \
  --max_range 4000 --step 10 --bearings 1 \
  --validation_csv fidelity_runs/table_midpoints.csv \
  --validation_samples 21 --validation_max_error 1 \
  -o fidelity_runs/table.nc
```

Validation evaluates the midpoint of every source-height, receiver-height and
range interval, and the midpoint of every periodic bearing interval. Singleton
axes keep their sole coordinate. Ground-model tables use direct ground-model
evaluation. PE tables use direct marches at these geometries, twice the
generation height-grid density, half the output-range spacing, and the
specified reference frequency count (21 by default).

The CSV holds table and reference levels, signed errors and statuses in SI
coordinates. Stderr reports bias, RMS, maximum absolute error and gap counts
per band. The default error budget is 1 dB. Exit code 2 means a finite reference
falls in a table gap, the budget is exceeded, or a band has no covered checks.
The table remains available for inspection. Reference values outside the PE's
angle domain are reported separately and are not scored.

This comparison includes both interpolation error and differences in numerical
resolution/band quadrature. Increase the reference frequency count again to
check convergence before treating it as an oracle. Passing at cell midpoints
does not guarantee accuracy everywhere or continuity at the fallback boundary.
Use reported failing coordinates to add source/receiver-height nodes, reduce
range spacing, or increase bearing density. The footprint lookup continues to
interpolate in dB; these changes do not alter that established behavior.

## Validation required before replacing a source database

Keep complete runs out of training. Score the held-out measurements by source
height, speed, emission angle, range, frequency and microphone installation,
alongside per-band SEL and total exposure. Report coverage, per-run variation
and clustered uncertainty; pooled energy means can hide opposing errors.
Repeat SNR-gate sensitivity checks. A propagation improvement on frames selected
by the old model's own shadows must also pass the other geometry regimes.

The latest local findings show that tunneling and skip-zone fill improve some
duct-shadow frames but can overpredict other far frames. They remain optional.
A range-dependent source fit would confound descent state drift with path loss;
the source work above preserves the existing path cap instead.

## Checks

`ctest --test-dir build-gcc --output-on-failure -j4` passed all 115 tests on
2026-10-04, including the new midpoint acceptance/rejection check and source
clamping diagnostics. The original golden outputs still pass.

```sh
python -B tests/test_fidelity_tools.py
```

The optional Python checks verify emission-clock alignment, guarded source
patching, separation of heading offsets, absolute held-out scoring, correction
extent, condition-bin bounds, frame alignment and common prediction coverage.
