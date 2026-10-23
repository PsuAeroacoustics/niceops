# Time-domain synthesis of received signals in NICE-OPS: an outline

*2026-10-23. This is a plan, not a result. The inventory was taken read-only from:*

- *NICE-OPS `gpu-metal` (3589b90);*
- *panam `directional-coherence` (f28426f, unpushed);*
- *the local PSU-WOPWOP (`wopwop3 2`, `development` b0feb06) and RNM/AAM (`AAM_Rau`) trees.*

*Numbers are labelled as **measured** (by a prior study), **checked** (read or computed for this outline) or **estimate**. Effort figures are estimates of focused engineering time.*

---

## 1. Goal

**Received signal.** This means acoustic pressure p(t) in Pa at a named receiver, for one trajectory, on absolute reception time. Each receiver has a type: bare point, pole microphone, or flush plate. Each sits in a NICE-OPS environment: open ground, terrain, refraction, or `--buildings`.

**Receivers, not grids.** Synthesis is for `-b` microphone lists only, and grids are refused.
- The demo track is 6551 × 0.02 s = 131 s. At 51.2 kHz that is 6.7 M samples, or 27 MB of float32 per receiver (checked).
- The 14,400-point demo grid would come to about 390 GB (estimate).

**Output format.** Follow the 2017 `*_pascal.nc` convention so that measured and synthesized records go through the same analysis code.
- **What a measured file holds** (checked with `ncdump -h` on `283421_01`):
  - one microphone per file (`nmics = 1`) and `float pressure(icount, nmics)` in Pa;
  - scalar globals `sample_rate`, `start_time`, `number_samples`, `mic_number`, `Latitude`, `Longitude`, `X`, `Y`, `Z`.
- **Files.** Write one file per microphone with the same attributes. A multi-microphone variant would need per-microphone coordinate variables, and would be a different format.
- **Sample rate.** This is a run option. By default it is the compared microphone's own `sample_rate`, read per file. The cached B407 hovers are 51,200 Hz (checked). panam's THEORY says the 2017 set also has 25.6 and 25.0 kHz by run and channel.
- **Start time.** `start_time` is an integer number of seconds in the files seen. Write the synthesized value to full precision. Measured records then need sub-sample alignment (§7.2).
- **Optional WAV.** float32 at 1.0 = 1 Pa, at the native rate, with an optional resample to 48 kHz for listening.

**What p(t) enables:**
- LAmax (F/S), SEL, TA and NAx taken from the waveform rather than from band energy × `reception_factor`.
- Metrics that need a waveform: loudness, sharpness, tonality, roughness, fluctuation strength, and impulsiveness (crest factor, kurtosis, pulse peak-to-RMS, `bvi_harmonic_ratio` in `panam/acoustics/impulsiveness.py`).
- Spectrogram and modulation comparisons against the 2017 records.
- Listening tests.
- Binaural signals for listening tests, later: each path already carries its arrival direction, so convolving it with a published HRTF set is a small addition to the renderer.

**Offline first.** Real-time limits such as low filter orders and wavefront clustering do not apply. The artefact lessons from real-time work still do.

**All in NICE-OPS.** The renderer is ours, not an external auralizer's. The milestones' main check is that the synthesized signal, analysed back to bands, reproduces NICE-OPS's own levels model by model, and a third-party renderer would put a second implementation of delays, filters and interpolation in the way of it. Real-time engines also take shortcuts this work does not need, such as low-order filters, which fail at the rotor harmonics (§4.2(d)). The code added is modest: delay lines, short per-path filters, an oscillator bank and filtered noise.

**Out of scope:** interactive, head-tracked or loudspeaker-array playback, which would need a real-time engine (revisit only if that becomes a goal), and nonlinear propagation.

---

## 2. Architecture

### 2.1 Synthesize at the source, then gain, delay and filter per path
This is the NASA Auralization Framework GTF pattern (Rizzi & Sahai review §5.3) and VA's outdoor renderer (Stienen 2022). The alternative is to synthesize at the receiver from received spectra. It is rejected for three reasons:

- **NICE-OPS already has the per-path quantities** that receiver-based synthesis has to reconstruct:
  - `RayPath` (`src/propagation.hpp`);
  - eigenray `travel_time` (`refraction.cpp`, `stratified.cpp`);
  - urban `PointRay` lengths and `point_ray_amplitude(ray, …, k)` (`point_paths.hpp`).
- **The receiver-based paper states the limit.** Aumann, Rizzi and Heath: mixed spectra cannot be split back into direct and reflected parts over finite-impedance ground, or for elevated receivers, whose two paths have different emission angles. Their conclusion is that a propagation program supplying spectra per path would resolve this. NICE-OPS can be that program.
- **Third-octave energy has already lost the direct/image phase**, and that phase is what makes the audible ground comb.

### 2.2 Who does what
| Component | Owner | Product |
|---|---|---|
| Source signal per emission direction and condition | panam | Harmonic plus residual sphere in the `.nod` (§3) |
| Paths: delay, launch and arrival direction, spreading, interaction transfers | NICE-OPS | Path tracks (§2.3) |
| Rendering: variable delay lines, time-varying filters, correlated noise fields | NICE-OPS (new) | p(t) netCDF and WAV |
| Comparison with measurement | niceops-validation harness | Derived tables only (§7.4) |

### 2.3 Path tracks as the renderer's internal product
A record per (state, receiver, path) does not scale. Urban would be 6551 states × 96 microphones × about 10³ rays, roughly 6×10⁸ records (estimate). So the unit is a **path track**, one per path key per receiver:

- **Key:** the interaction sequence (facades, edges, floor bounces, profile, kind), or the stratified family and bounce count.
- **Lifetime:** birth and death state indices.
- **Geometry on the state grid** (or at anchor states, with spline coefficients):
  - delay τ;
  - spreading;
  - absorption length;
  - launch direction in the sphere frame;
  - arrival direction.
- **Transfer:** the complex interaction transfer on a sparse log-k grid, with the bulk delay removed (§4.2).
- **Coherence:** the model's own coherence assumption for this path (§3.5).
- **Format:** netCDF-4 when written.
- **Specular and edge paths:** once a key is known, the path length can be recomputed analytically at every state without re-tracing. That fixes the delay-sampling problem caused by the urban model's 32-state stride (§5).

Path tracks are internal to the renderer. They are written to a file only for debugging and the closure tests, so their format can change freely.

---

## 3. Source side (panam)

### 3.1 What exists
**Spheres today:**
- `.nod` spheres hold third-octave magnitude only.
- The 2017 spectral databases carry 31 bands, 10 Hz to 10 kHz.
- Locally, `B407_fixed_load_factor.nod` has a `frequency` dimension. `B407_M3.nod` and `B407_ambient_gated.nod` do not (checked with `ncdump -h`).
- The AAM hemisphere is PHI 19 × THETA 37 (10° × 5°, −90..90) at r_ref = 100 ft. As NICE-OPS loads it, that becomes 36 × 37 = 1332 channels, with the upper half mirrored at coverage 0.
- `depropagate_hemisphere` drops phase at the spectrogram.
- Flight spheres are in **received** frequency, and NICE-OPS applies no Doppler (THEORY §11).

**Building blocks:**
- `geometry.dedopplerize`: retarded-time read at t_e + r/c.
  - Scales by r/radius only (checked), so its output keeps the received amplitude, convective amplification included.
  - Uses two linear `np.interp` calls.
  - Ignores wind, absorption, the board and the hub lever arm.
- `vold_kalman_filter`: complex envelopes per order.
- `order_tracking`, `ridge_track`, `estimate_clock_offset`.
- `spectra.tone_aware_band_power(return_parts=True)`, and `tone_split`, which writes separate AAM spheres and not `.nod` fields.
- `axisymmetric_bem.scattering`: complex P_d and P_r at any frequency. Only band tables are cached.
- `synthetic_flyover.synthesize_pass`: a closed-loop test source.
- `directional_coherence.py`: unpushed branch only. The harness submodule (d727b96) predates it.

**Absent:**
- time-synchronous averaging;
- angle-domain resampling;
- any rotor-azimuth reference;
- complex harmonic amplitudes per direction;
- wind-aware de-Dopplerization;
- complex receiver deconvolution.

**Linear-interpolation loss in `dedopplerize`.** At a fractional offset d the magnitude is |1 − d + d e^{−iω}|, which is |cos(πf/fs)| at the half-sample worst case (computed):

| fs | f | Worst case | Mean over offsets |
|---|---|---|---|
| 25.6 kHz | 5 kHz | −1.75 dB | −1.08 dB |
| 25.6 kHz | 10 kHz | −9.45 dB | −3.88 dB |
| 51.2 kHz | 10 kHz | −1.75 dB | −1.08 dB |

### 3.2 What the 2017 data constrain
**No rotor reference.** There is no tach and no once-per-rev signal anywhere in the B407 data: not in the AC tracks, the VN-200 raw files, the flight cards or `B407FullRefList.csv`. The rate has to be found acoustically. The comb search gives main BPF ≈ 27.6 Hz and tail BPF ≈ 83.4 Hz (measured).

**Across directions** (harness `coherence/README.md`):
- The fundamental is 0.94–0.99 coherent between any two directions, and its phase turns with shaft azimuth (measured).
- The rotor delay dψ/Ω helps up to about 4.5 BPF (124 Hz) and hurts from about 160 Hz. For example, the flyover 250 Hz band at 1–5° gives Re +0.64 without it and −0.16 with it (measured).
- The residual mean coherency is g = exp(−(fθ/850)^1.6). Its g = 0.5 angle is about 5° at 135 Hz and 1.4° at 500 Hz.
- Windowed |ρ|² stays at 0.4–0.6 within 10° at 500 Hz while the mean is low (measured). The relative phase is deterministic within a window and wanders between windows; the directions are not independent.
- **Consequence:** phase between directions cannot be interpolated on a 10° × 5° grid above roughly 120–150 Hz.

**Within one direction.** Phase across harmonics (pulse shape, BVI, crest factor) is a different quantity, and each microphone measures it directly. Store it separately from cross-direction phase.

**Timing (estimates):**
- Upos = 1 m is about 3 ms, which is 30° at 27.6 Hz and 130° at 125 Hz.
- The hub-to-antenna lever arm is unknown for the B407. 1–1.5 m would give about ±4 ms by direction. This is a hypothesis to calibrate.
- Wind adds about r·W·n/c², roughly 13 ms at 300 m in 5 m/s. Day 284's median wind was 11 kt.

**Dwell (estimates):**
- At 100 kt and 250 ft, one microphone crosses a 5° overhead cell in about 0.9 revolution.
- Hovers give about 55 revolutions in 8 s, about 17 dB of incoherent-ambient suppression.
- In-flight averages therefore need pooling across microphones and runs, which needs a shared Ψ.
- A drift of 0.1% in Nr over 8 s moves the BPF phase by about 80°.

### 3.3 Sphere designs
| | 1: Phase-locked harmonics + residual statistics | 2: Per-direction filters | 3: Depropagated grains |
|---|---|---|---|
| Holds | Complex shaft-order amplitudes (main, tail) relative to Ψ(t_e); harmonic coherence and revolution count; residual PSD; residual blade-rate AM | Magnitude, realised as minimum-phase FIR per direction | De-Dopplered pressure snippets per cell, Ψ at start |
| Reproduces | Tones, within-direction waveform shape, swish | Broadband spectrum. Minimum phase gives the wrong crest factor and kurtosis for BVI | Everything measured |
| Needs | Ψ reference, complex board deconvolution, wind travel times | Nothing new; can be built at load from `amplitude` | Grain alignment; ambient and board baked in |
| Size (estimate) | ≈2.2 MB per source condition at K_main 160, K_tail 48; ≈100–150 MB per database with zlib | Small | ≈6–24 GB, a separate file |
| Distributable | Above the pooling thresholds of §7.4 | Yes | **No.** Near-raw; local only |

**Recommendation.** Use Design 1 for the periodic part. Use Design 2 to shape its residual, and as the fallback for coverage-0 cells and level-only databases. Keep Design 3 as a local reference and listening tool only.

**Design 1 details to fix early:**
- **Shaft orders, not only BPF multiples,** so that blade dissimilarity is kept.
- **Store amplitudes de-rotated,** c′_k = c_k e^{−iksψ_dir}, and let the reader re-rotate. The azimuth phase then lives in the sphere and **never** also in a delay line (§3.5).
- **Reader rule:**
  - complex interpolation for k·f_shaft ≤ `coherence_rotation_limit_Hz`;
  - above that limit, interpolate the magnitude and take the phase relative to the fundamental from the nearest node, so the pulse shape survives.
- **Blending conditions.**
  - Blend in shaft order. Blend magnitude and aligned phase separately; a complex blend of unaligned offsets cancels.
  - Ω also differs between conditions. Synthesize Ψ(t_e) as the running integral of a blended Ω(t_e), continuous across states. The c_k come from the blend; Ψ does not restart at each state.
- **Load factor.** Store each waveform once per source condition. Load-factor groups refer to it through an index; they do not copy it. The +20 lg n scaling is magnitude-only, and crude for harmonic shape.

### 3.4 Predicted source: PSU-WOPWOP
- **How.** Use an `&ObserverIn` spherical grid `attachedTo` the hub, with periodic loading and tMin..tMax spanning one revolution.
- **Output.** `pressure.fn` is Plot3D stream binary, little-endian float32. The header is (nbObsDim1, nbObsDim2, nt, nbSources+1), and time is function 1 (`observerContainer.f90:WriteTimeHistoryFunctionFile`). `complexPressureFileName` gives the harmonic equivalent.
- **What it covers.** Deterministic thickness and loading only. Broadband (BPM, Pegg) comes out as band spectra.
- **It is not flyover-ready.** An attached observer is in the wind-tunnel frame: its frequency is emission-referenced, its amplitude and direction are not. The converter has to:
  - remap the attached-frame (reception-time) direction onto the emission direction that `hemigen` and NICE-OPS index by. The two differ by Mach aberration, of order M rad at the sides;
  - apply the Doppler-factor amplitude conversion to the measured-sphere convention (§3.6).
  Test the converter against a WOPWOP run with a stationary ground observer.
- **Prior art.**
  - `wopwop2aam.d` (Rob Rau; on `origin/dub_build`, not on gpu-metal) reads `spl_spectrum.fn` into AAM `AMPLITUDE` and drops phase.
  - RNM/AAM PT spheres (`Inventory.F`, `GetPT.F`, `ProPT.F`, `AddPT.F`) are the existing amplitude-plus-phase format. `GetPT.F` interpolates in p² and linear phase.
- **Cost.** A NICE-OPS-compatible sweep is about 164 flight conditions × 8 CT ≈ 1300 trimmed CHARM load sets (estimate). The sweep costs more than the reader does.

### 3.5 Coherence mapped onto signals
`src/source_coherence.hpp` is a second-order statistic: by default the cone is g = 1 within 0.75° and tapers to 0 at 1.25°. The measured form is exp(−(fθ/θ₀)^p). For signals it becomes a generator.

**Rule 1: the renderer reproduces each model's own coherence assumption,** so that energy closure holds model by model.
- Open-ground direct and image paths are coherent apart from T(f). The source cone is not applied there today. A 1.2 m microphone at short range has a direct–image launch difference above 1°, so applying the cone would delete the comb that G_b models.
- Stratified arrivals are summed incoherently today.
- Changing either is a model change. It is scored on the energy path first.

**Rule 2: under `--buildings`, the source is a direction-correlated random field, not discrete coherence groups.**
- Group ids are ill-posed: cone membership is not transitive, and paths cross the taper as the aircraft moves.
- Per band, draw signals whose cross-spectra between launch directions follow g(fθ) in the mean and the spread terms in the window statistic. Two ways to do it:
  - a common component times √g plus a per-direction component times √(1−g), smooth over direction;
  - a slow random relative phase with its decorrelation time matched to the 16/f-second window statistic.
- Slowly moving short-term interference between nearby paths is what the measurements describe. It is not an artefact to suppress.

**Rule 3: the rotor delay is a harmonic phase, not a delay-line delay.**
- dψ/Ω (`SourceCoherence::delay_length`) applies only to harmonic components with k·f_shaft ≤ the rotation limit.
- It never applies to broadband or to higher harmonics, where the data say it hurts.
- With Design 1 it is already in the re-rotation of c′_k and must not be applied again.

### 3.6 Frame and amplitude convention (must be agreed jointly)
**Frequency axis: emitted.** The delay line produces Doppler and the time compression, so a received-frequency source double-counts Doppler.
- Add a per-variable frame flag. The group-level `DOPPLER_SHIFT_REMOVED` is not enough.
- Proposed root attributes: `waveform_format_version`, `waveform_frequency_frame = "emitted"`, `waveform_phase_reference`.

**Amplitude: received, convective amplification included, for measured spheres.**
- That is what `dedopplerize` already produces.
- A delay line y(t) = x(t − τ(t)) keeps amplitude sample by sample and stretches duration by 1 + ṙ/c. Its energy is therefore ΣP_i T(1 + ṙ/c), which is exactly how `reception_factor` (`src/util/sel.hpp:88`) weights NICE-OPS's SEL (checked).
- **The renderer applies spreading and delay and no moving-source amplitude factor.**
- Predicted spheres are converted to this convention before they are written (§3.4). They are never compensated in the renderer.

**De-Dopplering band or PSD data** must conserve band power: S_e(f_e) = D·S_r(D f_e), with D = 1/(1 − M_r). Relabelling the frequency axis alone does not conserve power.

### 3.7 `.nod` change (Design 1)
- Add a `/waveform/source<N>` subgroup per source condition, holding:
  - `main_harmonics(PHI, THETA, k_main, ri)` and `tail_harmonics(…)`;
  - `harmonic_coherence` and `harmonic_n_rev`;
  - `residual_psd(PHI, THETA, band)` and `residual_modulation(PHI, THETA, band, m, ri)`;
  - shaft rates and their spread, and blade counts.
- Complex values use a trailing `ri = 2` dimension.
- Each existing group gets an `i4` attribute `waveform_source`.
- Lower hemisphere only. The mirrored half has no phase, so the reader falls back to Design 2 there.
- Readers that do not know `waveform_format_version` ignore it.
- **Cost risk.** netCDF group metadata already dominates load time, and nested subgroups add to the `nc_open` cost. Do the single-group layout fix in panam first, or at the same time.

---

## 4. Propagation side (NICE-OPS)

### 4.1 Model by model
| Model | Today | Delay | Phase internally | Time-domain treatment |
|---|---|---|---|---|
| Broadband EAA (`Database::eval`) | One dBA per pair | r/c | None | **Refuse.** Nothing to synthesize |
| Spectral + ISO 9613-1 (`eval_spectrum`, `atmosphere.cpp`) | Band dB | Straight range / c_surface | α(f) at any f | Delay line, 1/r, minimum-phase absorption FIR. With wind in `--atmosphere`, τ uses the effective sound speed along the path even without refraction (≈13 ms at 300 m, 5 m/s, estimate) |
| Ground, point receiver (`ground_effect.cpp`) | Chessell band G_b with T(f) | τ = (R2 − R1)/c, explicit | `spherical_reflection_coefficient` gives complex Q at any κ | Separate **image path**: its own launch direction, 1/R2, and an FIR of Q(f) with R2/c removed. Use a **causal** impedance (Miki, or a fitted rational/multipole model), not Delany–Bazley; check that h_Q(t) is zero before onset. Drop G_b |
| Plate (`plate_receiver.cpp`) | mean ‖P_d + Q P_r‖² | mm | P_d, P_r complex, 5 per band | Short FIR from the complex response interpolated in log f, or recomputed from `axisymmetric_bem.scattering` on a finer grid |
| Terrain (`terrain.cpp`) | CNOSSOS dB; Kirchhoff ‖F(v)‖² | δ computed, unused | Fresnel F(v) analytic | Kirchhoff: complex F(v), delay from δ. CNOSSOS: minimum phase, band-smoothed |
| `--refraction linear` | Band dB; `ShadowTransition` V(ξ, s) power | `DirectEigenray.travel_time`; the reflected ray has none | Residue series complex, stored as power | Eigenray delays. Add `ReflectedEigenray::travel_time`. V is used as one band-smoothed path (see rule (c) below) |
| `--refraction stratified` | Band dB; arrivals incoherent; caustic cap +10 dB | Per-`Arrival` `travel_time`, launch, arrival, spreading, bounces | — | One tap per arrival. **Add a caustic count** (KMAH): each caustic is −π/2, i.e. a Hilbert transform of that arrival. The +10 dB cap is an energy device |
| PE (`pe.cpp`) | Band dB; tables magnitude only | ψ phase (`solve_pe` keeps it; `solve_pe_band` does not) | ψ per frequency | Cost goes as f², so no dense grid. Magnitude as a band-smoothed path with the ray delay. Tables could later store complex ratio or group delay |
| Urban point (`point_paths`, `urban_point`, `src/gpu/*`) | Band energy (pair sum or GPU quadrature) | Every `PointRay::length` | `point_ray_amplitude(s)` complex at any k: UTD, floor Q, plate, shared faces | Already an echogram. A_i(k) on the sparse log-k grid, FIR per ray, tap at L_i/c. GPU: read back the complex A_i(k) that `ray_amplitudes` already forms |
| Edge source BTM (`edge_source.cpp`, CPU) | Complex spectra via Filon in path length | Exact onset | Any k | Keep in the frequency domain on the sparse grid, like every other path. The partner corrections subtract UTD per k and have no time-domain counterpart |
| Row-of-roofs PE (`profile_pe.cpp`) | Complex `ProfileParts` per sub-band k | Relative to taut string | Sparse in k | Interpolate magnitude and unwrapped phase in log k |
| Urban tables, canyon wave | dB, ρ | First-order images computed, not exposed | 2-D FDTD pulse, discarded | Low priority |
| Diffuse facades (`diffuse_facades.cpp`) | Steady radiosity energy | None | None | Decorrelated noise tail with band envelopes; needs time-dependent radiosity or a decay estimate per band |
| Turbulence scattering (`turbulence_scattering.cpp`, PE `ShadowScattering`) | Energy, dB re free field | None | None | Band-limited noise decorrelated from the rays, at band energy, a short spread after the direct delay, with slow modulation |
| FDTD (`fdtd.cpp`) | Native p(t) (`--recording`) | Exact | Exact | Oracle only |

### 4.2 Rules that cut across models
**(a) Sparse transfer grid, not FFT bins.**
- Once L_i/c is taken out, the remaining transfers are smooth in log f: Q, UTD, P_d/P_r, F(v), and the shared-face and profile phase relative to the taut string.
- The existing 3 sub-bands per third octave (`sub_bands()`, about 90 points to 10 kHz) are enough to design a 64–256-tap FIR per path per update (estimate).
- Dense bins are about 20–200× the cost (estimate) and are kept as a validation option.
- **Interface:** `evaluate_transfer(path, k[])` plus FIR design, with `point_ray_amplitudes` as the template.

**(b) Delay-free filters.** Every per-path FIR is designed with its onset delay removed. All delay, including the fractional part, lives in the delay-line read.
- Estimate the impulse-response length per interaction from the spread of path-length differences.
- Set the design length from that estimate and window the result. This guards against wrap-around from UTD transition functions, BTM tails and shared-face apertures.

**(c) No double counting.**
- Band-averaged energy stand-ins are never applied on top of the explicit paths that replace them. This covers Chessell's sinc average, T(f), the urban pair sum and incoherent arrival sums.
- PE relative levels and V(ξ, s) already contain the ground reflection. For each receiver and path, use either:
  - the explicit direct + image pair with Q(f), or
  - one path filtered by the PE or V band magnitude, labelled band-smoothed (no comb).
  Never use both.

**(d) Minimum phase only where it is valid.**
- **Absorption: valid.** It is causal, and dispersion is negligible. At kilometre range α·L spans more than 60 dB, so zero-pad the cepstral FFT, floor the magnitude (e.g. −120 dB), and check the result against direct convolution.
- **Interference spectra: not minimum phase in general.** 1 + (R1/R2)Q e^{−iωτ} is non-minimum-phase wherever |Q R1/R2| > 1, which happens near grazing over soft ground. A magnitude-only term that contains ground interference gets rule (c), not a minimum-phase FIR.
- **Avoid low-order IIR below about 80 Hz.** It under-fits there (Stienen §5.4), which is where the BPF harmonics are.

**(e) Frequency coverage.**
- The databases stop at 10 kHz. Above the top band, define an explicit roll-off (for example −12 dB/octave), which is a decision for the user.
- Tones below about 100 Hz need the harmonic source of §3, not bands.

---

## 5. Time variation

| Quantity | Today | For synthesis |
|---|---|---|
| States | 0.02 s on shipped tracks | Adequate for magnitudes (NAF hop 11.6 ms) |
| Stratified | Every 4th state, `lerp_scaled` | Linear interpolation puts slope kinks in τ. Respline, or solve every state |
| Urban | Every 32nd state + bisection in dB; no phase, no identity | Path tracks with analytic length per state (§2.3), or `--state_tolerance 0` with key matching |
| Delay | — | Per audio sample |

**Delays.**
- Fit τ(t_e) per path track as a **smoothing** spline (C² or better), with the tolerance set from Upos and a stated bandwidth limit on Doppler variation.
- Raw AC positions enter the run directly, since resampling never coarsens, and δf/f = δv_r/c turns position noise into FM.
- Acceptance check: the residual FM on a synthetic tone.
- A linear ramp per state would step the Doppler at 50 Hz (Schäfer et al. 2024).

**Fractional-delay read.**
- Use windowed sinc or Lagrange of order ≥ 3, with oversampling where tones alias. See Randall & Andrew (2025) on aliasing and imaging from time-varying delay.

**Where directivity is applied.** Each path has its own launch direction, so it has its own source spectrum or c_k. Two placements:
- **Emission time (default).** Synthesize the per-direction source signal, then apply the per-path delay. There is no shared buffer across paths with different directivity.
- **After the read, in reception time.** Then the directivity's frequency axis must be scaled by that path's own 1/(1 − M_r,i). Without the scaling, band edges are misplaced by up to about 15% at M_r = ±0.15 (estimate).
- A shared source buffer with several read cursors (VA's SIMO) is valid only for paths that share a launch direction to within the directivity's resolution.

**Filters.**
- Update at the state rate between delay-aligned responses.
- Linear tap interpolation between minimum-phase FIRs is not minimum phase in between. That is harmless when neighbouring responses are close; otherwise crossfade, or interpolate the magnitude and rebuild.
- Prefer per-wavefront time-variant processing to crossfading whole impulse responses (Stienen).

**Broadband noise.**
- Do **not** use random-phase OLA blocks at a 20 ms hop. The 10–20 Hz third octaves are 2.3–4.6 Hz wide and need blocks of at least about 1 s.
- Instead, run continuous white noise in emission time through a fixed third-octave (or finer) filter bank, with per-sample-interpolated band gains.
- Krishnamurthy & Rizzi's 2 s blocks with a 11.6 ms parameter hop are the cited alternative; they blur fast changes.

**Birth and death.**
- Paths appear and vanish at:
  - shadow boundaries;
  - reflection and visibility boundaries;
  - new stratified arrivals;
  - urban re-traces.
- Fade over a stated window, never switch hard.
- Where a model has a continuous transition, use it rather than a time crossfade: the shadow transition's ξ = −2 → −1 blend, and the edge-source/UTD zone continuity.
- **Measure before designing.** Log path keys at consecutive evaluated states on one DigiCity microphone and count births and deaths per second. Then choose between crossfades and every-state evaluation. Sample a few seconds of track only.

**Full march.**
- Synthesis joins `needs_full_history` (`niceops.cpp:2834`), so there is no windowed peak search, `--db_cutoff` or culling.
- Day-night penalties stay in emission time; synthesis is in reception time throughout.

---

## 6. Milestones

**M0. Ground-path export.** About 1 week (estimate).
- Deliverable: per state and receiver, write the run's band transfer together with τ = (R2 − R1)/c, R1, R2, the grazing angle and Q at band centre.
- Check: a standalone re-implementation of Chessell's G_b from those fields reproduces `--band_history` to 1e-9 dB. This tests the decomposition, not the physics.

**M1a. Direct-path audio, harness prototype, no C++ changes.** 1–2 weeks (estimate), in `niceops-validation`.
- Input: `--band_history --no-ground_effect` on `B407_fixed_load_factor.nod`, the local spectral B407 database, for `-b` microphones. This file already has `emission_time`, `reception_time`, range, angles and unweighted bands for every state (`write_band_history_file`, `niceops.cpp:268`).
- Render: filter-bank noise, a windowed-sinc delay read from `reception_time`, and output in `*_pascal.nc` layout.
- Frequency: the spheres are received-frequency, so each path gets its own buffer with band edges de-Dopplered by its own (1 + ṙ_i/c) and power-conserving density scaling. This is labelled approximate, because the bands were measured over shifted edges.
- Checks:
  - whole-event SEL from ∫p² equals `--full_history` SEL within ±0.1 dB, as the mean over ≥ 10 seeds;
  - per-band SEL within 2σ, with σ ≈ 4.34/√(BT) dB;
  - no frame-by-frame check in low bands. At 25 Hz in a 0.5 s frame, σ ≈ 2.5 dB.

**M1b. Ground image and plate in C++.** 3–4 weeks (estimate).
- Deliverable: the explicit image path with a causal Q(f) FIR, plus the plate P_d + Q P_r FIR. The 2017 ground boards are flush plates.
- Checks:
  - the band-averaged comb equals the **exact** band integral ∫|1 + Q(R1/R2)e^{−i2πfτ}|² df from `spherical_reflection_coefficient` on a fine grid, within the band-filter-shape error;
  - Chessell's own error against that integral is reported separately, so that it is not mistaken for a synthesis bug (estimate: tenths of a dB in low bands);
  - with turbulence off only, unless an Arntzen & Simons T/R filter pair is included.

**M1c (optional). Tones from `tone_split`.** Synthesize them at emitted k·f_shaft, with Ω from the vehicle JSON (`Main Rotor.tip speed`/`radius`) or the coherence record, and amplitude from tonal band power. Needs a panam writer change and a reader change.

**M2. First measured comparison.** 1–2 weeks (estimate).
- Compare at held-out microphones only (§7.2), on flush-plate microphones only.
- Hovers first, because they are the only data on disk. They test the source and not the delay line.
- Flight runs are fetched run by run through `spheres/fetch.py`.
- Compare band frames, LAmax and SEL above the 10 dB SNR gate. A stochastic band source cannot match tonal spectrogram structure.
- Pass: no worse than NICE-OPS's own band-level error at the same microphones (from the harness tables).

**M3-spike. Virtual tach.** 2 weeks (estimate). This is the go/no-go for Design 1.
- Build Ψ(t_e) acoustically from de-Dopplerized emission-time signals, refined with the VKF.
- Check it against logged RPM in `example_data/F100_20_P#1[FA1FZC-3.h5].mat` (public `origin/master`).
- Measure the stability of Ψ on B407 hovers.
- If Ψ cannot be held under about 30° at BPF (a proposed threshold), Design 1 falls back to per-direction magnitude plus within-direction phase relative to the fundamental.
- Written fresh on public master (§7.4).

**M3. Design 1 builder (panam).** 2–3 months (estimate).
- `dedopplerize` rewrite as its own item (1–2 weeks): windowed sinc, wind travel times from `uniform_wind_ray_model`, per-channel rates, lever arm. Closed-loop test on an extended `synthesize_pass`.
- Complex board table and per-harmonic deconvolution. Exclude pole microphones; cap the correction near nulls.
- Builder on hovers, then level passes; then the `.nod` writer, `validate_nod`, and `docs/file_formats.md`.
- Prerequisites: merge `directional-coherence` and bump the harness submodule.
- Checks: held-out harmonic magnitude, within-direction phase relative to the fundamental, and impulsiveness metrics frame by frame (thresholds §8).

**M4. Harmonic synthesis in NICE-OPS.** 4–6 weeks (estimate).
- Oscillator bank on the running integral Ψ(t_e(t)) with re-rotated c_k(launch).
- Residual through the filter bank with blade-rate AM.
- Correlated source field under `--buildings` (§3.5).
- Convention test, two cases on the demo track:
  - the measured-convention sphere with no renderer factor reproduces `--full_history` SEL within ±0.1 dB;
  - a deliberately wrong configuration (an extra 1/(1 − M_r) factor) fails.

**M5. WOPWOP reader plus one condition.** 3–4 weeks (estimate).
- `pressure.fn` / complex-pressure reader, aberration remap and amplitude conversion (§3.4), checked against a stationary-observer WOPWOP run.
- The (μ, γ, CT) sweep is a separate item and is costed separately.

**M6. Refraction and terrain.** 3–4 weeks (estimate).
- `ReflectedEigenray::travel_time`, stratified arrival tracks with caustic counts, Kirchhoff screen paths, and band-smoothed PE and V paths.
- Check: magnitude and phase at named geometries against single-frequency `solve_pe` ψ. The band-smoothed PE path is checked for magnitude only.

**M7. Urban.** 2–3 months (estimate).
- First the birth/death measurement (§5).
- Then `PointRay` path tracks with A_i(k) on the sparse grid: CPU first, then GPU read-back.
- Check against `niceops_fdtd_field --recording` on a named small scene, within the FDTD's resolved band.
- A time-domain BTM (Svensson's sample-bin integration of the level-curve density) is a later option, blocked on time-domain partner corrections.

**M8. Stochastic terms.** 4–6 weeks (estimate).
- Diffuse facade tail, turbulence-scattering noise, and per-path phase and log-amplitude sequences (Rietdijk et al.).
- The sequences have:
  - a frozen-turbulence correlation time from the path's transverse speed through the air (aircraft speed by relative position, plus wind);
  - a **cross-path structure function** for paths that share air (direct and image near the source), not independent draws.
- Check: ensemble-mean energy within the 2σ seed spread, and the mean comb reproduces Chessell's T(f).

---

## 7. Validation

### 7.1 Self-consistency (every milestone)
- **Energy closure is to a tolerance, not to 1e-9.**
  - A time-varying sinc read aliases and images in proportion to dτ/dt and the kernel length. State an interpolator error bound at the top band.
  - Random parts close in expectation. Write the tolerance as a multiple of σ/√N.
- **The 1e-9 checks are export checks only** (M0).
- **Each model keeps its own coherence assumption** (§3.5), so closure is per model.
- **CPU and GPU.** Closure tests run on CPU. A separate GPU-vs-CPU test bounds the read-back difference in dB per band of ∫p². GPU amplitudes are `complex<float>` with float-float phases, about 3e-7 rad.

### 7.2 Against 2017 records
- **Home.** Scripts and derived tables live in `~/Code/Python/niceops-validation`. Aligned measured snippets and ambient stay outside git.
- **Not circular.**
  - The spheres were built from these runs. Rebuild with `panam.synthesis.leave_mics_out` (`synthesis.py:451`), or exclude the harness's designated held-out runs.
  - Compare only at the excluded microphones and runs, with one fixed list for M2–M4.
- **Alignment.**
  - `start_time` is good to 1 s only. Estimate the sub-sample offset per run and microphone with `estimate_clock_offset`, and record it with its uncertainty.
  - Check that `utcsec` and `start_time` use the same time scale. GPS and UTC differed by 18 s in 2017; this is a hypothesis that has not been checked.
- **Ambient.** Compare only bands and frames with SNR ≥ 10 dB, matching `band_snr_gate_db`. Compute impulsiveness on band-limited signals with the same gate applied to both records.
- **Phase.**
  - Compare magnitude, and phase relative to the fundamental within a direction.
  - Absolute or cross-direction phase only below `coherence_rotation_limit_Hz` (124.1 Hz for B407), and only after the virtual tach passes M3-spike.
- **Comparisons:**
  - band frames (`--band_frames` against a band analysis of p(t));
  - LAmax and SEL;
  - harmonic magnitude and relative phase;
  - blade-rate modulation depth;
  - crest factor, kurtosis and `bvi_harmonic_ratio`.

### 7.3 Psychoacoustics and listening
- Compare loudness, sharpness, tonality, roughness and fluctuation strength against the recordings. Then run a similarity or plausibility test (Pieren et al.; Rietdijk).
- Krishnamurthy & Rizzi report that phase modulation carries impulsiveness: RMSE 21.5% with measured phase modulation against 54.4% without. These figures are as reported in the full text and have not been re-checked; confirm them before using them as targets.
- Listening tests need IRB approval at Penn State. Plan the lead time.

### 7.4 Data boundary
- Raw `*_pascal.nc` pressure, measured clips and spectrograms, ambient recordings and Design 3 grains are raw and stay local.
- **Design 1 coefficients are publishable only** where `harmonic_n_rev` ≥ N and at least M microphones or runs are pooled. Cells below that threshold are local-only, and so is any audio that depends on them. In flight a single cell can be one revolution of one microphone, which is near-raw.
- The coherence record JSON is derived and fine to use.
- The virtual tach and harmonic separation are written fresh on public panam master. The private beta26 separation workflow and its demo audio are not a source, and are not ported or quoted. Generic `vold_kalman_filter.py` on master is fine.

---

## 8. Decisions needed, with recommendations

1. **Amplitude convention.** Recommendation: measured spheres carry received amplitude on an emitted-frequency axis; the renderer applies no moving-source factor; WOPWOP spheres are converted before writing (§3.6). Also decide whether WOPWOP and measured spheres may share one database. Recommendation: yes, only after the M5 converter test passes.
2. **Sphere design and K_main.** Recommendation: Design 1 + Design 2, Design 3 local only. K_main 160 covers BVI at about 100–150 MB per database (estimate); 80 halves it.
3. **Phase above about 124 Hz between directions.** Recommendation: not deterministic. Use magnitude plus the correlated field of §3.5, or WOPWOP where deterministic phase is wanted.
4. **Coherence in the renderer.** Recommendation: reproduce each model's current assumption; any change goes through the energy model first.
5. **M1 Doppler with received-frequency spheres.** Recommendation: per-path band-edge de-Doppler with power-conserving scaling, labelled approximate. The alternative is a static render.
6. **Output sample rate and the ceiling above 10 kHz.** Recommendation: the compared microphone's own rate, with −12 dB/octave above the top band.
7. **Held-out runs and microphones.** Recommendation: the harness's existing held-out list, fixed for M2–M4.
8. **Numeric thresholds** (proposals, not measured):
   - harmonic magnitude within X dB for k ≤ 10 BPF in hover and ≤ 4–6 BPF in flight;
   - within-direction relative phase within Y° rms;
   - impulsiveness metrics within Z%;
   - virtual-tach Ψ error under about 30° at BPF.
9. **Publishability thresholds N and M** for Design 1 cells.
10. **Fetch budget** for cloud-only flight data. B407 acoustics are 86.8 GB against about 48 GB free on 2026-09-27 (measured then). Recommendation: per-run fetches of held-out runs only.
11. **Turbulence.** Recommendation: per-path sequences with cross-path structure (M8). Until then, close with turbulence off. Leaving turbulence out is not neutral for annoyance (Lincke et al. 2025).
12. **Urban time variation.** Recommendation: measure birth and death rates first (§5), then choose.
13. **IRB** for listening tests.
14. **Phase modulation (BVI jitter).** No time-extensible model exists in the literature. Recommendation: flag it as a known gap and measure the statistics in M3 as research.

---

## 9. References

### Literature
- Rizzi & Sahai, auralization review: https://ntrs.nasa.gov/api/citations/20200002351/downloads/20200002351.pdf
- Aumann, Rizzi, Heath, receiver-based NAF: https://ntrs.nasa.gov/api/citations/20200002559/downloads/20200002559.pdf
- Aumann et al., NAF architecture (InterNoise 2015): https://ntrs.nasa.gov/citations/20160006430
- Krishnamurthy & Rizzi, AM helicopter flyover: https://ntrs.nasa.gov/citations/20200002058
- Pera et al., AS350 synthesis (SciTech 2018): https://ntrs.nasa.gov/citations/20180002090
- Krishnamurthy, Tuttle, Rizzi, F1A in NAF: https://ntrs.nasa.gov/citations/20200002468 ; plugin: https://ntrs.nasa.gov/citations/20200011447
- Boucher et al. 2023: https://ntrs.nasa.gov/api/citations/20230002310/downloads/TUSQ_v8_final.pdf
- Stienen, *Real-Time Auralisation of Outdoor Sound Propagation* (RWTH 2022): https://www.logos-verlag.eu/ebooks/OA/978-3-8325-5629-7.pdf
- VACore outdoor renderer: https://git.rwth-aachen.de/ita/VACore/-/blob/9ff1749209faaeb1c7b474601e6be39124b0992b/src/Rendering/VAOutdoorNoiseAudioRenderer.h
- Schäfer, Fatela, Vorländer 2024: https://acta-acustica.edpsciences.org/10.1051/aacus/2023070 ; ART: https://acta-acustica.edpsciences.org/10.1051/aacus/2021018
- Arntzen & Simons, NLR-TP-2015-282: https://reports.nlr.nl/bitstream/10921/1014/1/TP-2015-282.pdf ; Arntzen thesis: https://repository.tudelft.nl/file/File_f1fa2900-30d2-4fc8-92c3-24e99b4e90a9
- Rietdijk, Forssén, Heutschi 2017, turbulence sequences: https://research.chalmers.se/en/publication/250825
- Lincke & Pieren 2024: https://acta-acustica.edpsciences.org/articles/aacus/abs/2024/01/aacus240030/aacus240030.html ; Lincke, Kawai, Pieren 2025: https://acta-acustica.edpsciences.org/10.1051/aacus/2025053
- Mascarenhas, Cotté, Doaré 2022: https://perso.ensta-paris.fr/~cotte/files/Mascarenhas_Cotte_al_jasa2022.pdf
- Damiano & van Waterschoot, pyroadacoustics: https://dafx2020.mdw.ac.at/proceedings/papers/DAFx20in22_paper_9.pdf
- Randall & Andrew 2025: https://acta-acustica.edpsciences.org/articles/aacus/abs/2025/01/aacus240026/aacus240026.html
- Carpentier, SMC 2024: https://smcnetwork.org/smc2024/papers/SMC2024_paper_id118.pdf
- J.O. Smith, Doppler simulation: https://ccrma.stanford.edu/~jos/pasp/Doppler_Simulation.html
- Erraji, Stienen, Vorländer 2021: https://acta-acustica.edpsciences.org/10.1051/aacus/2021010
- Palenda & Vorländer 2026: https://acta-acustica.edpsciences.org/articles/aacus/abs/2026/01/aacus250149/aacus250149.html
- Dreier & Vorländer 2024: https://acta-acustica.edpsciences.org/10.1051/aacus/2024026
- Pieren et al.: https://elib.dlr.de/135813 ; Rizzi et al., AIAA 2020-2582: https://elib.dlr.de/129935
- Svensson, Fred, Vanderkooy, JASA 106:2331 (1999).
- Dragna & Blanc-Benon, time-domain impedance models (causal ground for item 4.1).

### Code
**NICE-OPS (`~/Code/D/NICEOPS-gpu`, 3589b90):**
- `src/util/sel.hpp` (`reception_factor`);
- `src/niceops.cpp` (`write_band_history_file`, `needs_full_history`);
- `src/propagation.hpp`;
- `src/sphere_database.cpp`;
- `src/atmosphere.cpp`, `src/atmosphere_profile.cpp`;
- `src/ground_effect.cpp`;
- `src/plate_receiver.cpp`;
- `src/terrain.cpp`;
- `src/refraction.cpp`;
- `src/stratified.cpp` (`Arrival`, `lerp_scaled`);
- `src/pe.cpp`, `src/profile_pe.cpp`;
- `src/point_paths.{hpp,cpp}`, `src/urban_point.cpp`, `src/gpu/band_energy.hpp`, `src/gpu/point_rays.hpp`;
- `src/edge_source.cpp`;
- `src/diffuse_facades.cpp`, `src/turbulence_scattering.cpp`, `src/source_coherence.hpp`;
- `src/fdtd.cpp`;
- `src/band_frames.cpp`.

**panam (`~/Code/Python/panam`, `directional-coherence`):**
- `geometry.dedopplerize`, `hemigen`;
- `depropagation.depropagate_hemisphere`;
- `noise_abatement_2017.build_sphere`;
- `vold_kalman_filter`, `order_tracking`;
- `spectra.tone_aware_band_power`;
- `axisymmetric_bem.scattering`;
- `synthetic_flyover.synthesize_pass`;
- `synthesis.leave_mics_out`;
- `directional_coherence`;
- `acoustics.impulsiveness`;
- `database.build_empirical_database`.

**PSU-WOPWOP (`~/Documents/Code/Fortran/wopwop3 2`):**
- `observerContainer.f90` (`WriteTimeHistoryFunctionFile`, `CalcObsContComplexPressure`);
- `binaryIO.f90`;
- `wopwop3_manual/Output.tex`.

**RNM/AAM (`AAM_Rau`):** `Inventory.F`, `GetPT.F`, `ProPT.F`, `AddPT.F`.

**Converter:** `wopwop2aam.d` (NICE-OPS `origin/dub_build`).
