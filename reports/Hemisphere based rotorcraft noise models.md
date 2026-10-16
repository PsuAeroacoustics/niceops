# NICEOPS leads on propagation, lags on sources

> **Status.** This is a snapshot of NICEOPS `stratified-ducts` at `e472b6f`, written on 2026-10-03. Most of the ranked shortcomings below have since been addressed:
>
> - PR #8 made ground effect the default, refused a missing humidity, normalized DNL to 24 hours and refused vehicle mismatches.
> - PR #9 started reading the upper hemisphere, made the seam periodic, read the Doppler flag, and added static operations, NAx and coverage.
> - PR #15 added `--terrain`, EPNL, time-weighted LAmax, Lden/CNEL, `--blend` policies, `--nose` and `--wind`.
>
> Still open: a single compact source, no Doppler shift applied, and no ambient-pressure scaling. `docs/MANUAL.md` describes the current behavior. The validation harness's `paper/model_review.md` (2026-10-05) supersedes the tool comparison here. It covers the build side as well, and corrects the SELENE and HELENA entries: SELENE does not floor its shadow zones, and HELENA blends four (V, γ) hemispheres by inverse distance.

Every research-grade rotorcraft footprint tool reviewed here does the same core thing NICEOPS does. That covers AAM v3.1, NORAH2, RNM and its FRAME-QS descendant, DLR SELENE, HELENA, Airbus CAROT/MOTUS, ANOPP2/AMAT and the 2023–26 eVTOL work. Each tool stores a fixed-radius sphere or hemisphere of band levels for each steady operating condition, picks or blends spheres along the trajectory, rotates them into the flight frame, and propagates along rays. The tools differ in four things: which state variables index the spheres, the domain they blend in, how much of the sphere they trust, and how much propagation physics sits downstream. AEDT/HNM is the exception, using noise-distance curves at three lateral angles.

On the first and last of those four, NICEOPS is at or beyond the published state of the art. Its (μ, γ_e, C_T) Delaunay blend is the indexing the maneuver literature recommends. No footprint tool found combines its spectral ground model, layered-eigenray refraction, residue-series shadow transitions and a benchmarked PE. The closest is the TU Delft ray model used by Gennaretti et al. (IJA 2017): layered arcs, ray-tube spreading, and Arntzen's shadow-zone model with ground-wave, diffraction and turbulence-scattering terms. Its internal harness result, 1.26 dB rms SEL over 1,112 held-out microphone-runs of B407 approaches, is a multi-run error distribution of a kind no sphere tool has published. It is not yet published itself.

Its weaknesses sit elsewhere:

- **Source coverage at the envelope edges.** The hover is synthesised, steep descents are copies of level spheres, and the rim is zero-gated. This is where its largest measured residuals live: +3.9 dB at the low, slow end of approaches and −1.5 to −2.8 dB on the fast forward rim.
- **Defaults that are not the validated configuration.** A default run has no ground model and assumes 20% relative humidity.
- **No terrain or screening.** AAM and NORAH2 both model it.
- **Thin outputs.** NICEOPS reports A-weighted SEL only, and its "DNL" is 49.4 dB above a standard DNL.
- **One whole-vehicle source.** It does not read the Doppler flag, and it has no static operating modes.
- **Documentation drift.** CLAUDE.md and the README now misstate the SEL window, the spectral default, the assert policy and the map projection.

Blending in dB rather than energy, the obvious point of difference from NORAH2, matters least. Our own measurement puts it at +0.25 dB mean in SEL. That is a quarter to a seventh of the 1.1–1.8 dB rms leave-one-out interpolation error, and the dB–energy difference in that test (0.02–0.05 dB rms, apart from B407's single deep-null sphere) is a small fraction of it. dB is the only blend that is well defined on NICEOPS's degenerate condition set.

All NICEOPS code references below are to the `stratified-ducts` branch at `e472b6f`. Where the code-verified note ([verification_niceops.md](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/verification_niceops.md)) disagrees with the inventory, the verified reading is used. The inventory was written at `e472b6f`; the note was run at `fb9fad5`. Of the files it checks, only `niceops.cpp`, `trajectory.*`, `ground_effect.*`, `refraction.*` and `stratified.*` differ between the two, and the line numbers cited here have been re-checked at `e472b6f`. AAM page numbers are the printed TR page numbers; NORAH2 page numbers are the printed page footers of each document.

## One lookup table, five lineages, and a regulatory retreat from interpolation

The US lineage runs from RNM (NASA Langley and Wyle, late 1990s) to AAM, which Volpe has owned since spring 2020 ([AAM TR p.ii][aamAnn]).

**RNM.** The original RNM indexed hemispheres by airspeed, flight-path angle and pylon angle. It took the *nearest* γ and pylon angle, then interpolated airspeed only, "assum[ing] a linear variation in pressure squared". Within a hemisphere it used a Coons-surface interpolation in (θ, φ), and it stepped the aircraft at 0.5 s ([Lucas & Marcolini, c. 1997](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)). By 2002 the same sentence described two successive linear interpolations "in the energy domain", first in airspeed and then in flight-path angle ([Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)).

**AAM v3.1.** The current AAM Technical Reference repeats that sentence with "in the energy domain" removed ([TR p.31][aam31]). Two pages later it says twice (and again for fixed-wing on p.36–37) that up to four spheres are "linearly interpolated in decibels", in speed and then in FPA, and that this procedure "remains" ([TR p.33–34][aam33]). The documented AAM v3 method is therefore **dB-linear, like NICEOPS**. Either the domain changed at some RNM version or one of the two documents is wrong; the TR is the later and more specific source. The search has hard gates: 3° in FPA on approach, 5° on departure and 20 kt in speed. FPA may disqualify spheres but speed may not. Ties at the speed midpoint go to the faster sphere ([TR p.33–34][aam33]), and a sparse database falls back to the nearest sphere ([TR p.31][aam31]).

AAM adds three things on top of the RNM search:

- **Q-SAM1/Q-SAM2.** These substitute an effective flight-path angle. Q-SAM2 derives the thrust direction n̂ = (a + g)/|a + g| and the roll from a waypoint spline ([TR p.39–40][aam39]).
- **FRAME.** The FRAME database type adds "effective weight" (W·LF) as a third axis, with up to eight spheres per lookup ([TR p.41–42][aam41]).
- **NCSPEC.** This names a sphere explicitly per trajectory point. It is the recommended mode for AAM vehicles, because automatic lookup "is limited to the NC file header data" ([TR p.30–31][aam30]).

AAM does not account for wind in source lookup: ground speed stands in for airspeed ([TR p.226][aam226]).

**FRAME-QS.** This is the research end of the US lineage, and the closest published analogue to NICEOPS's state handling. It holds about 500 FRAME-predicted spheres over (V, γ, n_z) and blends them by natural-neighbour weights on a 3-D Delaunay triangulation. Its γ_e is the angle between the velocity and the inertial tip-path plane, sin γ_e = v·n̂/V, with n̂ the inertial TPP normal (the direction of a + g). On an aggressive Bell 430 pitch-up it caught the transient OASPL rise "within 2 dB in all three directions"; holding C_T fixed missed about 3.5 dB, close to 20 log10 n_z ([Greenwood, Rau, May & Hobbs 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).

**European research tools.** These converged on the main-rotor aerodynamic state:

- **DLR SELENE** stores hemispheres by μ, α_TPP and C_T on a Delaunay mesh with a nearest-neighbour fallback outside the envelope ([Olsman & Lummer](https://dspace-erf.nlr.nl/bitstreams/61aeceb7-6f3b-486f-8d7e-d335cd4db426/download)). That is almost exactly NICEOPS's `interpolator.cpp`.
- **MANOEUVRES** re-indexed its maps from (V, γ, W) to (μ, C_T, α_TPP) ([Trainelli et al. 2015](https://dspace-erf.nlr.nl/server/api/core/bitstreams/132a1b61-32a5-4197-9f86-3fe20f153d5b/content)).
- **Gennaretti et al.** concluded that "at least three parameters, including disk loading, are required to adequately characterize noise hemispheres" ([Gennaretti et al., AST 2016](https://consensus.app/papers/details/6d70fa64cca3505e90368ff77076782d/?utm_source=claude_desktop)).
- **HELENA**'s indexing is not documented in anything retrievable; a "three closest flight conditions" rule is unverified. It propagates with SAE ARP 866A absorption, Zaporozhets–Tokarev ground reflection and ray-traced refraction ([ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)).
- **Airbus CAROT/MOTUS** splits measured data into main-rotor tonal, anti-torque tonal and residual hemispheres, stores the tonal ones de-Dopplerized, and applies simulation-derived "masks" for design changes ([Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)). Its parameterisation and interpolation are unpublished.

**The regulatory line went the other way.** NORAH1 forbade interpolation outright: "No interpolation between flight conditions is allowed", and the vertical profile is "dictated by the available hemispheres" ([van Oosten et al., Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)).

NORAH2 (EASA, 2023–24) added a 2-D Delaunay blend over normalised (V, γ). γ is multiplied by F_fc = 2, "empirically determined to give minimum interpolation error", with no error figure given. The three triangle vertices are combined by an **inverse-distance-weighted energy average**, not barycentric weights ([D1.5d §A.3.1 p.12–13](https://www.easa.europa.eu/en/downloads/132005/en)).

ECAC Doc 32 (1st edition, endorsed 19 May 2026) is built on NORAH but **did not adopt the blend**: "there is no agreed method for the interpolation between hemispheres", so users assign a hemisphere per constant-speed, constant-angle segment. It disclaims single-event accuracy ([ECAC Doc 32](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf)). NICEOPS's barycentric blending is therefore outside the endorsed European method. That absence of an agreed method is ECAC's whole stated reason. No published measurement of interpolation error exists either, and the experiment below supplies one for NICEOPS's databases.

**The regulatory US tools have no sphere.** HNM, inherited by INM and AEDT's helicopter mode, stores SEL-versus-distance at eight distances from 200 to 10,000 ft. These cover three moving modes at left 45°, centre and right 45°, plus static-mode LAmax curves with an eight-term azimuthal Fourier directivity ([HNM methodology, 1988](https://rosap.ntl.bts.gov/view/dot/10018)). AEDT/INM has 16 helicopter modes: 8 dynamic and 8 static or vertical (ground and flight idle, HIGE, HOGE, and vertical ascent/descent in and out of ground effect) ([ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)).

NASA's 2025 recommendations say the ±45° extrapolation "may not accurately reflect" lateral noise and that the fourth-power noise fraction is "only marginally applicable" to helicopters and eVTOL ([Letica & Rizzi, NASA/TM-20250008619](https://ntrs.nasa.gov/citations/20250008619)). The scale of the problem is known. Inside 45° slant, AEDT/INM and AAM/HELENA differ by about ±1 dBA SEL on the MD-902. Toward the rotor plane the difference is "3-5 dBA SEL" ([ACRP 02-44 App. B](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)).

**eVTOL work keeps the pattern.** AARON computes lower-hemisphere F1A spheres that AMAT or AAM fly ([Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf); [Lopes 2022](https://ntrs.nasa.gov/citations/20220018599)). The MIT ASCENT 084 hemisphere library reports two useful results ([ASCENT 084 2025](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-084-2025-Annual-Report.pdf)):

- Hemispheres 13× coarser move ground levels by only 1–3 dB on average.
- Tilt-rotor transition is a gap that no three-parameter index covers.

Averages hide local error. In ACRP 02-44, cutting a Bell 412 sphere's lateral spacing to 30–45° moved SEL(A) by up to 2.9 dB and Lmax by up to 9.4 dB at a 6,300 ft sideline, and by 0.0 dB under the track ([ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)). NICEOPS's 10° φ × 5° θ grid is three to four times finer laterally, and the harness puts its bilinear reading within about 0.3 dB of the depropagated samples at 4–20° elevation ([niceops_inventory.md §1](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/niceops_inventory.md)).

## The tools side by side

The table condenses the comparison. "AAM" is v3.1.0 as documented in the May 2021 TR. "NORAH2" is the D1.5d method plus the V2.0.74 public package, with ECAC Doc 32 noted where it differs. Cells marked *(ours)* are findings from our own inspection of shipped files, not from the documents.

| Aspect | NICEOPS (stratified-ducts, e472b6f) | AAM v3.1 | NORAH2 / ECAC Doc 32 | RNM / FRAME-QS | SELENE / HELENA / CAROT | AEDT / HNM |
| --- | --- | --- | --- | --- | --- | --- |
| **Source format and grid** | netCDF `.nod`, one group per condition. 30.48 m sphere, 36 φ (10°) × 37 θ (5°). dBA + EAA, plus optional 31 one-third-octave bands, 10 Hz–10 kHz (the AS350, B407, Be206L3, EC130, RO-44 and RO-66 `_fixed_load_factor` databases have them; S-76D_M3 does not). Full sphere stored. | One NetCDF file per source per condition. Any regular (θ, φ) grid, "typically" 5° ([p.144][aam144]). Broadband 1/3-OB, narrowband, 1/12-OB or pure tones with phase. Radius per file. 1E35 = no data ([p.145][aam145], [p.163–164][aam163]) | ASCII `.hem` (the HELENA format). 19 × 19 lower hemisphere, 10° bins, 60 m. 31 bands, with absorption to 60 m at ICAO reference conditions included. 38–67% of bins filled, −999 elsewhere *(ours)*. Hover is an azimuth-only table ([D1.5d p.14–15, App. A](https://www.easa.europa.eu/en/downloads/132005/en)) | RNM: broadband and pure-tone hemispheres, 5° typical. FRAME-QS: about 500 predicted spheres | SELENE: fixed-radius hemispheres. HELENA: hemispheres or "noise carpets". CAROT: separate MR tonal, TR tonal and residual hemispheres | No sphere. NPD at 8 distances × 3 lateral angles. Static modes: an azimuth-only Fourier fit (HNM), or 15° steps at 200 ft (AEDT) |
| **Operating-state index** | (μ, γ_e or γ, C_T); `-q` (γ_e) is the default. C_T = LF·W/(ρV_tip²πR²). On `fixed_load_factor` databases: (μ, γ) + 20 lg(C_T/C_T,ref) above the reference, uncapped | (V, γ) for helicopters, plus nacelle angle (nearest) for tiltrotors. FRAME mode adds W·LF ([p.30][aam30], [p.41][aam41]) | (V, γ) only, using ground-relative 3-D speed. Hover, idle and taxi are discrete modes with fixed dB offsets ([D1.5d p.11–19](https://www.easa.europa.eu/en/downloads/132005/en)) | RNM: (V, γ, pylon). FRAME-QS: (V, γ_e, n_z) | SELENE: (μ, α_TPP, C_T). CAROT: unpublished | 16 discrete modes; speed enters through a 10 lg(V_ref/V) duration term |
| **Cross-condition blend** | Barycentric over 4 spheres on an exact-predicate 3-D Delaunay mesh, **in dB**. Hull padded by 10%. Nearest condition outside the hull. No distance gates | Up to 4 spheres, **in dB**, speed then FPA. Gates: dFPA 3°/5°, dSPD 20 kt ([p.33–34][aam33]). Nearest fallback ([p.31][aam31]) | 3 vertices, IDW (1/δ) **in energy**, over normalised (V, 2γ). Nearest outside the hull ([D1.5d p.12–13](https://www.easa.europa.eu/en/downloads/132005/en)). ECAC Doc 32: none | RNM 1997: nearest γ, then **p²**-linear in V. FRAME-QS: natural neighbour | SELENE: Delaunay, nearest outside. HELENA: unverified | Mode substitution |
| **Orientation** | TPP normal = load-factor unit vector. Nose = (air) velocity projected into the TPP; below 5 kt, the INS heading instead, as a hard switch (`9fc714a`). A track with no heading column, or with too few states over 20 kt to relate the heading to the frame and no `--frame_bearing`, falls back to the drift, with a warning; before `9fc714a` every hover did | Body axes from user-supplied roll, pitch (AOA + γ) and yaw. Legacy small-angle, zero-roll ANGLES routine below thresholds. MAST TILT. Q-SAM2 derives the TPP normal and roll ([p.24–25][aam24], [p.39–40][aam39]). Roll + = port down, opposite to NORAH2's ([p.20][aam20]) | Heading + slip, rolled by bank. **No pitch** (verified to 0.003° *(ours)*). Departure turns bank outward, a sign bug *(ours)*. Build side (NORAH1): attitude ignored, x-axis along track | FRAME-QS: rotated by the inertial TPP normal | SELENE: aircraft attitude angles. HELENA: TPP angles from the EUROPA performance code | Track-aligned only |
| **Within-sphere interpolation** | Periodic bilinear **in dB**; per band on the spectral path | Not in the TR. ACRP: "linear on (phi, theta) for each one-third octave band", domain unstated | Bilinear **in energy**. Empty bins take the nearest filled bin (great circle) as a constant ([D1.5d p.15–16](https://www.easa.europa.eu/en/downloads/132005/en)) | RNM: Coons surface | Undocumented in every retrievable source (TUM eVTOL work: linear, after Doppler) | Linear in φ between the 3 NPD directions; ±45° held beyond |
| **Upper hemisphere** | The stored upper half is panam's mirror (φ → 180 − φ). The lookup uses \|up\|, so no propagated level ever reads the upper half (only `--self_test` and the sphere exporters touch it) | Full spheres are legal (OMNI_200 spans ±180°). Hemisphere files populate only φ −90…+90 ([p.214][aam214]); how AAM treats directions above them is undocumented (1E35 nodes inside a file are not propagated, [p.145][aam145]). WRAPAROUND mirrors the lower half ([p.161][aam161], [p.204][aam204]) | None. SingleEvt clamps φ to ±90° at the same θ *(ours)* | FRAME: computed; the halves differ | HELENA: lower hemisphere (the .hem format NORAH2 inherited). DLR: lower, from a ground disk to 300 m. CAROT: undocumented | None: ±45° lateral curves held beyond 45° |
| **Doppler** | No shift applied; spheres keep the received shift (panam writes `DOPPLER_SHIFT_REMOVED = 0` in the per-run `.nc` files; the `.nod` has no flag). The flag is **never read**. Reception-time weight 1 + ṙ/c | DOPPLERSTATE flag. f′ = f/(1 + u/c), with a difference correction from test speed to modelled speed. Ground model at the shifted frequencies ([p.10–11][aam10]) | None. Reception time t_e + r/c, with c = 346.1 m/s ([D1.5d p.20](https://www.easa.europa.eu/en/downloads/132005/en)) | RNM: Doppler shift in propagation | SELENE: shift plus amplification. CAROT: re-Dopplerizes its de-Dopplerized tonal spheres | Implicit in the NPD |
| **Time step and reception time** | Trajectory's own step (0.02 s on shipped tracks); reception-time weight 1 + ṙ/c | 2 s default segmentation (1 s typical, 0.2–0.5 s for ART); emissions binned by arrival time, then interpolated to a uniform mesh ([p.10][aam10]) | Multi-event DTIME 0.5 s; SingleEvt uses every input point; t_r = t_e + r/c with fixed c = 346.1 m/s ([D1.5d p.20](https://www.easa.europa.eu/en/downloads/132005/en)) | RNM: 0.5 s | SELENE: emission time solved from t − τ − \|R\|/c₀ = 0, Doppler amplification included | Segment-based; no time step |
| **Absorption** | ISO 9613-1 per band, with an optional layered T/RH/p profile. Broadband path: the database EAA, linear in range, built at 20 °C / 20% RH | ANSI S1.26 (R2004), layered, tabulated at 1000 ft. SAE-AIR-1845 option ([p.13–14][aam13]) | SAE ARP 5534 with Rickley band conversion; uniform atmosphere ([D1.5d p.20–22](https://www.easa.europa.eu/en/downloads/132005/en)) | ANSI S1.26-1978, altitude-averaged | HELENA: SAE ARP 866A | In the NPD, plus SAE-AIR-1845 |
| **Ground** | **Off by default.** Chien–Soroka with Chessell boundary loss and band averaging. Four impedance models, roughness, turbulence coherence. GR1425 plate BEM; point receivers per height. Flat ground at z = 0 | Chien–Soroka/Chessell/Daigle with flow resistivity. Incoherent per band, coherent for pure tones. Chessell decoherence ([p.10–11][aam10]) | Chien–Soroka + Delany–Bazley with sinc coherence. Default grass, 1.2 m microphones. Mean ground plane over terrain ([D1.5d p.22–27](https://www.easa.europa.eu/en/downloads/132005/en)) | Two-ray, finite impedance | HELENA: Zaporozhets–Tokarev | Empirical lateral-attenuation lump |
| **Terrain and shielding** | **None.** Observer z moves only the direct path | Rasmussen GTD wedge/screen, per-cell impedance, ISO 9613-2 foliage. Cannot be combined with wind ([p.10–12][aam10]) | DEM, CNOSSOS-EU screening (0–25 dB). Building obstacles are in the method but not supported by the prototype ([Final report §3.2](https://www.easa.europa.eu/en/downloads/139729/en)) | RNM v3+: GTD | CAROT: topography | Slant-range adjustment only |
| **Refraction, wind, turbulence** | `linear` (equivalent-linear arcs) or `stratified` (layered eigenrays, every arrival summed in energy, including a duct's bounced rays). Shadow-transition tables, creeping-wave loss, PE propagation tables (usable beside either ray model); scattering and the duct's Airy tails (`--tunnelling`, `--skip_fill`) opt-in | RNMwea: precomputed omnidirectional tables over flat uniform ground, single event, Berry–Daigle shadows ([p.14–16][aam14]) | **None** ("not covered") | RNM 7: stratified with winds. APET: GFPE tables | SELENE: ray tracing with a level floor in shadows. HELENA: ray tracing | None |
| **Metrics and SEL window** | A-weighted SEL only. 10 dB-down window from first to last crossing. `-D` is a penalised energy sum, not divided by 86 400 s | Lmax, SEL (A and C), DNL, CNEL, PNLT, EPNL, NA, TA, TVL. Full-history SEL by default; DBDOWN10 (a sample threshold) is an option ([p.50–52][aam50], [p.64][aam64]) | LAmax, SELA, PNLTM, EPNL, NAx, Lden. Full history by default (ISEL10DB = 0) ([Manual Annex A p.46][nman]) | RNM: SEL, Lmax, EPNL, DNL and others | 10 dB-down SEL (HELENA, SELENE) | SEL, LAmax, EPNL, DNL |
| **Maneuvers** | Q-SAM γ_e, C_T from load factor, TPP tilt. No transient or wake term | QSAM1: γ + v̇/g. QSAM2. FRAME: 8 spheres ([p.39–42][aam39]) | Turns tilt the *level* hemisphere by the bank angle, even in descent. No acceleration | FRAME-QS: within 2 dB on a pull-up | MANOEUVRES: (μ, C_T, α_TPP) | Decelerating and accelerating modes (B/C/E/F) |
| **Multiple sources** | One whole-vehicle sphere at the track point. The tail-rotor JSON block is unread | Up to 10 broadband + 10 narrowband + 10 pure-tone spheres at body offsets; tones summed coherently ([p.30][aam30], [p.40–41][aam40]) | One | RNM: multiple | CAROT: MR, TR and residual (+ engine) | One |
| **Sphere generation** | panam: steady runs over flush boards. Ambient gate at 10 dB with subtraction. Plate BEM divided out (capped at 15 dB). ISO absorption undone (capped at 30 dB), with curved rays where profiles exist. Adaptive Shepard in power on a 2° working grid (weight floored, reach ×5 in azimuth), written at 10° φ × 5° θ; 2° rim. Mirrored upper half, synthetic hover and γ, load-factor copies | ART: fly OMNI_200, subtract every propagation term per band. SNR gate at 5 dB, no subtraction. ZGRID. 1E35 holes. High-frequency absorption cap at 5 kHz ([p.145–161][aam145]) | Measured: 17 inverted plate microphones. Last-good-band repair. Scaled to 60 m, binned at 10°, repeat runs merged in energy ([Final report §2.2.3](https://www.easa.europa.eu/en/downloads/139729/en)) | FRAME: physics model calibrated to measurement | DLR: back-propagated to a 300 m sphere | AAM-simulated certification events → NPDs |
| **Validation** | Internal, unpublished: B407 day-286 approaches, 1,112 mic-runs, **1.26 dB rms** SEL (station atmosphere), 1.37 (stratified). Windy day 284: about 1.7. S-76D: "within 2 dBA" | Qualitative ("agree well"); no statistics ([p.3][aam3], [p.12][aam12]) | Against certification: ≤ 1.8 EPNdB. Against Nord2000: ±1 dB. No operational comparison ([D3](https://www.easa.europa.eu/en/downloads/139664/en)) | RNM XV-15: 0.5 dB, but self-consistency only. FRAME-QS: 2 dB on 2 maneuvers | SELENE: "good agreement". HELENA: no numbers | Model-to-model only: ±1 dBA SEL against AAM/HELENA inside 45°, 3–5 dBA toward the rotor plane (MD-902, ACRP 02-44). 1.7 dB SELA under-track against AAM for a UAM quadrotor. A 40 dBA contour at over twice AMAT's lateral extent (NASA TM 2025) |

Sources for the NICEOPS column:

- Grid, blend, fold, Doppler, metrics, refusals and window: verified in code ([verification_niceops.md §§1–10](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/verification_niceops.md)).
- Validation: the harness README ([niceops-validation README](/Users/evg5332/Code/Python/niceops-validation/README.md)) and the grazing-incidence findings ([Grazing incidence findings](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md)).
- S-76D result: [Pascioni et al., VFS 76, 2020](https://ntrs.nasa.gov/api/citations/20205008035/downloads/Pascioni_VFS2020.pdf).
- The sphere-generation recipe and the day-286 approach figures (1.26 / 1.37 dB SEL rms, +3.9 dB at the slow end) are for the B407 database of `5d761b8`, the build `e472b6f` ships. That build replaced `0decafa` (straight-ray depropagation, a fixed 10° Shepard grid with a 25° radius, no cap on the board correction; the 2° floor is unchanged) by way of `56ce39b` and `041f39b`. The options study (§1), the day-284 windy-day figures and the four footprints (§11) are harness README results measured at earlier NICEOPS commits (`d4c868d`, `70ce85e`, `d852613`) on the straight-ray build or its day-284-held-out variant, and have not been rerun on `5d761b8`.

### The NORAH2 rows rest partly on reverse-engineering

The NORAH2 frame findings come from refitting the `.his` emission angles of the shipped example runs ([NORAH2 package, ARP cases](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/>)). Heading plus roll reproduces THETAOBS and PHIOBSAC to 0.003° rms. Adding a pitch tilt by the vertical angle breaks the fit by 3–54° in φ. So SingleEvt applies **no pitch and no γ tilt**, consistent with D1.5d's statement that pitch "is included implicitly" ([D1.5d §A.3.3 p.16](https://www.easa.europa.eu/en/downloads/132005/en)).

The departure example banks its hemisphere **away from** the turn centre in both turns (0.026° rms against 5.7° for inward bank), while the approaches bank inward. The cause is the multi-event module: it gives each `R` track segment −atan(V²/gR) and each `L` segment +atan(V²/gR) regardless of the direction in which the segment is flown. One of the two senses must therefore be wrong, and D1.5d eq. 20 says it is the departure. The φ error is 2|bank|. For the example's 2–3° banks that costs about 0.2–0.3 dB typical and 0.8–1.0 dB at the 90th percentile, along the turns only. This is a sensitivity estimate; NORAH2 was not rerun. The same refit shows that SingleEvt clamps any direction above the equator to φ = ±90° at the same θ. Scripts: [norah2_frame](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/scripts/norah2_frame/).

NICEOPS cannot make this class of error, because its bank is the load-factor direction derived from the path itself. It is a warning, though, for any future input column that carries roll explicitly.

## Exchanging spheres with AAM and NORAH2

The three formats agree on the angle grid and disagree on almost everything stored with it:

| Convention | NICEOPS / panam | AAM v3.1 (ART) | NORAH2 (.hem) |
| --- | --- | --- | --- |
| θ, φ | θ from nose; φ = 0 down, + starboard | Same (ANOPP) ([p.29][aam29]); the RNM 1997 paper had + port | Same (D1.5d eq. 11) |
| Roll sign | None stored (bank from the path) | + = left (port) down ([p.20][aam20]) | + = starboard down (SingleEvt fit *(ours)*) |
| Heading | Nose from (air) velocity, or below 5 kt from the INS heading (column 8), turned into the local frame by an offset inferred from the track or by `--frame_bearing` | N = 90°, E = 0° ([p.20][aam20]) | Clockwise from north |
| Radius | 30.48 m | Per file; ART default 100 ft | 60 m (other radii rescaled on load) |
| Absorption in stored level | None (lossless) | None (lossless) | To 60 m at ICAO 298.15 K / 70% RH (FREEFIELD = 2) |
| Ground reflection in stored level | None (plate BEM divided out) | None (removed by AAM's own ground model) | None for flight; included in the hover table (FREEFIELD = 0) |
| Doppler | Kept; per-run `.nc` files carry `DOPPLER_SHIFT_REMOVED = 0`, the `.nod` carries no flag | DOPPLERSTATE 1 = moving; NetCDF flag 1 = removed; Bullwink 1 = stationary | Kept, no flag |
| Speed field | μ from ground speed | SPEED labelled airspeed, used as ground speed ([p.226][aam226]) | ACSPEED labelled IAS; lookup by ground-relative 3-D speed |
| No data | Floor −1000 dB, dropped by `blend_present` | 1E35, not propagated | −999, nearest-bin constant fill |
| Opposite rotor rotation | Not handled | Not handled | φ → −φ for mirrored types |

Sources: [aam_v3.md §§1, 3, 5](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/aam_v3.md), [norah2.md §§2–3](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/norah2.md), [norah2_aam_checks.md §A](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/norah2_aam_checks.md), [hemisphere_generation.md §5](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/hemisphere_generation.md).

panam already writes both formats (`write_aam_hemisphere_netcdf`, `write_norah2_hemisphere`, `write_norah2_triangulation`) and reads `.hem` (`load_norah2_hemisphere`). NICEOPS exports AAM spheres (`--export_aam`) but imports neither format. Two hazards follow:

- **Ground speed in an IAS field.** panam writes the run's mean ground speed into NORAH2's ACSPEED field (`noise_abatement_2017.py:490, 993`). That is harmless inside NORAH2, which looks up by ground speed, but wrong for anyone who reads the field as labelled.
- **Lost error cancellation.** An ART sphere had its ground reflection removed with AAM's own Chien–Soroka model and flow resistivity. Propagating it with NICEOPS's plate or ground models loses the error cancellation AAM gets ([aam_v3.md §11](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/aam_v3.md)).

## Blending in dB costs a quarter decibel; picking the nearest sphere costs up to four

The literature has no number for the blending domain. NORAH2's F_fc = 2 was tuned "to give minimum interpolation error" with no figure published, and ECAC cites the absence of an agreed method. So this review measured it directly.

All numbers in this section are **our own computation**. We re-implemented NICEOPS's broadband path in numpy, copying its normalisation, the 10% boundary duplication, the Delaunay weights, clamping and the nearest fallback. The re-implementation reproduces the NICEOPS LAE output to ≤ 2.4e-13 dB on all 43 S-76D and B205 tracks at the default 120 × 120 grid. We then swapped only the condition blend for three alternatives:

- an energy blend with the same weights;
- a NORAH-style inverse-distance energy blend over the same simplex;
- the nearest stored condition.

The data were the shipped `S-76D_M3`, `AW139_M1`, `B205_M1` and `B407_M2` databases, with the S-76D, AW139 and B205 tracks; B407 has no shipped track and was evaluated at 20,000 random states ([experiment scripts](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/scripts/blending/)). The scripts, and those for the NORAH2 roll refit above, are committed with the research notes together with their summary outputs. The per-track state CSVs and intermediate arrays are not included; `scripts/blending/states.sh` regenerates the states. The scripts could later move to `niceops-validation`.

The comparison ran against an existing NICEOPS build, not one made at `e472b6f`. The broadband path it reproduces does not differ between `fb9fad5` and `e472b6f`: `src/sphere_database.cpp`, `src/interpolator.cpp`, `src/delaunay.cpp`, `src/predicates.cpp` and `src/util/sel.hpp` are unchanged between `fb9fad5` and `e472b6f`, and so are those four databases (`B407_M2` is not the B407 build of `5d761b8`, which changed only `B407_fixed_load_factor.nod` and `B407_ambient_gated.nod`) and the tracks used. The one change in `src/niceops.cpp` and `src/trajectory.cpp` that reaches a broadband run, `9fc714a`'s heading nose, applies only below 5 kt. No shipped track goes below 12.5 kt.

### The size of the effect

**The energy–dB gap is set by the spread of the levels being blended, and is small where it matters.** Over footprint directions of real trajectory states (|φ| ≤ 60°, 20° ≤ θ ≤ 160°), the energy blend exceeds the dB blend as follows:

| Database | Median | p95 | Maximum |
| --- | --- | --- | --- |
| S-76D | 0.10 dB | 0.83 dB | 5.3 dB |
| AW139 | 0.05 dB | 0.50 dB | 3.6 dB |
| B205 | 0.15 dB | 1.21 dB | 6.1 dB |

Jensen's inequality fixes the sign. The magnitude tracks the spread of the blended levels: about 0.03 dB at a 1–2 dB spread, 0.22 dB at 3–5 dB, 0.53 dB at 5–8 dB, and 1–2.7 dB at 8–20 dB.

**In the SEL footprint the integral averages most of it away.**

| Footprint set | Variant | Mean | p5 / p95 | Min / max | Cells beyond ±1 dB |
| --- | --- | --- | --- | --- | --- |
| S-76D, 30 tracks, 432,000 cells | energy − dB | **+0.26** | +0.09 / +0.52 | −1.23 / +2.41 | 0.4% |
| | NORAH-style IDW − dB | +0.28 | +0.03 / +0.64 | −1.30 / +2.46 | 0.8% |
| | nearest − dB | +0.19 | −0.73 / +1.08 | −3.19 / +4.04 | **9.8%** |
| B205, 13 tracks, 187,200 cells | energy − dB | **+0.25** | +0.09 / +0.59 | −2.27 / +2.93 | 0.4% |
| | nearest − dB | −0.07 | −0.74 / +0.49 | −2.36 / +4.00 | 1.9% |

The nearest-sphere offset is not a bias. It swings from −0.36 to +0.54 dB from track to track, depending on which sphere happens to be nearest. The energy offset is +0.18 to +0.42 dB on every track.

SEL can *fall* under an energy blend, by up to 2.3 dB, even though the instantaneous level never does. A louder peak raises the 10 dB-down threshold and pulls the window in; the window endpoints moved for 80% of S-76D observers. The window still moves under the current first-to-last rule, because the threshold moves with the peak. CLAUDE.md's decimation table already covers that rule in its outermost-crossing column (1.13–1.54 dB maximum cell error), so decimation stays unsafe under it; only the table's first-crossing column and the prose around it describe the superseded rule.

### Accuracy against held-out spheres

Leave-one-out tests whether a blend can predict a measured sphere it has not seen. Each measured flight condition is held out, along with its C_T and padding copies, and predicted from the rest. The errors are taken over footprint directions:

| Database (held out) | dB rms / bias | energy rms / bias | IDW-energy rms / bias | nearest rms | dB beats nearest |
| --- | --- | --- | --- | --- | --- |
| S-76D (82) | **1.84** / −0.02 | 1.86 / +0.16 | 1.89 / +0.24 | 2.15 | 80% |
| AW139 (51) | **1.12** / 0.00 | 1.14 / +0.08 | 1.12 / +0.11 | 1.32 | 88% |
| B205 (60) | **1.58** / 0.00 | 1.63 / +0.16 | 1.58 / +0.20 | 1.67 | 78% |
| B407_M2 (106) | 1.76 / −0.05 | 1.54 / +0.08 | **1.53** / +0.12 | 2.04 | 85% |

The dB blend is unbiased and the energy blends run high by 0.08–0.24 dB. The rms values are indistinguishable except on B407. There, one held-out sphere with a 26 dB null dominates, and energy wins by damping it; by median per-sphere rms, dB wins on B407 too (1.08 against 1.13). Interpolation error itself, about 1.1–1.8 dB rms, is roughly ten times the dB–energy difference.

**The blending domain is therefore a systematic offset question against NORAH2 (+0.25 dB from the condition blend), not an accuracy question.** Database density and measurement scatter matter an order of magnitude more. Tuinstra and Sijtsma reach the analytical counterpart: dB averaging sits 0.03, 0.15 and about 0.6 dB below the energy mean for ±1, ±2 and ±6 dB scatter. They also show that energy-domain directivity fits give "erroneous results", including negative energies, under strong directivity ([Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)).

This covers the condition blend alone. NORAH2 also interpolates within each hemisphere in energy, on 10° × 10° cells against NICEOPS's 10° × 5°, and fills 33–62% of each hemisphere from the nearest bin. Neither step was included here. Jensen's inequality makes the energy bilinear step add to the offset, and the fill can move it either way, so the full NICEOPS–NORAH2 offset is unmeasured. The NORAH-style variant also used NICEOPS's 3-D simplex, not NORAH2's 2-D (V, 2γ) mesh.

### Why NICEOPS cannot simply switch to energy

The shipped conditions form a degenerate prism. Every flight condition carries eight synthetic C_T copies at 0, 1, 1.1, 1.2, 1.5, 2, 3 and 5 × C_T,ref. Those above 1 × are exactly the reference sphere plus 20 lg(C_T/C_T,ref), and the 0 × copy repeats the reference. The Delaunay triangulation of that set is therefore not unique.

A dB blend is **invariant** to which valid triangulation is chosen: the maximum difference measured between Qhull's `Qt` and `QJ` meshes was 0.000 dB. The levels are additive in (μ, γ) and C_T, and any barycentric weights reproduce both marginals. An energy blend is not invariant. It moves by up to 0.8–1.5 dB between two equally valid meshes (p95 0.10–0.19 dB), which would make the result depend on Bowyer–Watson tie-breaking.

The CT axis itself contributes at most 0.02–0.11 dB to the gap at trajectory load factors. The exact 20 lg C_T law sits midway between the dB chord and the energy chord, within 0.28 dB even over 3–5 g.

An energy option for NORAH2 parity would therefore need a different scheme: a 2-D (μ, γ) triangulation crossed with linear or analytic C_T, as the `fixed_load_factor` path already does. A formula swap alone is not enough.

Two incidental findings from the experiment matter more than the domain:

- **Padding copies.** A third to half of each database's flight conditions are copies of a measured sphere at γ = −24° and +35°. Outside the measured γ range (−8.7° to +14.4° for S-76D), the blend is therefore constant extrapolation of the edge spheres.
- **Nearest fallback.** Falling back to the nearest condition, as AAM does for sparse databases and NORAH2 does outside its hull, is the least accurate method on every database. NICEOPS reserves it for queries outside the hull, which is the right call.

## Ranked shortcomings, with evidence and fixes

The ranking weighs measured or well-bounded dB impact against how many runs each item touches. Items 1–4 change real predictions by a decibel or more today. Items 5–8 limit what NICEOPS can be used for, though item 7's nose switch also moves predictions at the hover end. Items 9–11 are about trust in the code and the record.

| # | Shortcoming | Expected impact | Who does it better | Fix effort |
| --- | --- | --- | --- | --- |
| 1 | The validated configuration is not the default | About 1 dB rms SEL (calm day), more in wind; median 1.2–3.0 dB from the atmosphere defaults | AAM, NORAH2 (ground model on by default) | Small (defaults and CTest pins) |
| 2 | Source coverage at the envelope edges | +3.9 dB at the slow end of approaches; −1.5 to −2.8 dB on the fast forward rim; hover bands 4–19 dB off | NORAH2 method (dedicated hover arrays; the public package ships placeholders), HNM (measured HIGE over hard and soft ground), FRAME (physics fill), ART (honest 1E35 holes) | Large (panam coverage mask, censored gate, hover data) |
| 3 | No terrain, screening or non-flat ground | Up to ±10 dB under screening; 35 dB in terrain shadow for straight rays | AAM (GTD), NORAH2 (CNOSSOS), sonX, UYGUR | Large (DEM, mean ground plane, screening) |
| 4 | Unrecorded reference conditions; silent mismatches; no ambient-pressure scaling | 0.75 dB per 10% weight error; ~3 dB thickness noise per 4000 ft of site altitude (scaled from 4.7–5.4 dB over 7000 ft); load-factor law unvalidated above ~2 g, uniform in direction | NORAH2 metadata, FRAME non-dimensional scaling | Small to medium (root attributes, JSON check, p_amb scaling) |
| 5 | Metrics: SEL only; non-standard DNL; window-only SEL | +49.4 dB on "DNL"; no LAmax/EPNL comparison possible | AAM, NORAH2 | LAmax/NAx small, EPNL medium |
| 6 | One compact source; Doppler flag ignored | Large for design studies (main rotor alone "severely" underpredicts); ~2 dB on low-frequency tones fore and aft at 100 kt if the flag is misread (our arithmetic), ≤ 0.5 dB when conventions match (ASCENT 038) | AAM (30 spheres), CAROT/MOTUS | Small for the flag, large for multiple sources |
| 7 | No static operations; the hover nose switches hard from velocity to heading at 5 kt | 2.2–4.4 dB steps at the switch for a 30–90° nose rotation (library-level; how often real tracks present such angles is unmeasured) | NORAH2, HNM/AEDT, AAM | Small (smooth nose blend), medium (static records) |
| 8 | Upper hemisphere folded by \|up\| | 0 dB on today's mirrored databases; blocks genuine upper data | AAM (full spheres), FRAME | Small (atan2), but `DirectionBounds` must then bound both halves |
| 9 | Silent numerical failure paths | Up to 4.4 dB (inward spreading), +6.02 dB (coherent ambient), 1.8 dB (spline ends), 0.43 dB (`--db_cutoff`) | — | Small (merge the branch) |
| 10 | Documentation drift | Misleads anyone reasoning from CLAUDE.md | — | Small |
| 11 | Validation breadth and UQ | Error bars not propagated; B407 at one flat site dominates, a second aircraft scored, a six-aircraft set not yet summarised; no cross-tool sphere comparison | Nobody yet; this is open | Medium |

### 1. The validated configuration is several flags away from the default

The defaults do not match the configuration that was validated:

- **Spectral is the default but only sometimes runs.** The code defaults to spectral (`src/niceops.cpp:546-552, 2251, 2357-2361`) and falls back to broadband silently when the database has no spectra. The demo database `S-76D_M3.nod` has none, so the default demo and the regression golden are broadband. The same command line therefore runs different physics depending on which `.nod` it is pointed at.
- **The ground model is off.** No ground effect is applied by default.
- **The atmosphere defaults are silent.** Humidity defaults to 20%, and pressure is derived from density·R·T. The code review measured these silent defaults at 1.2–3.0 dB median ([code review X3](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md)).

The harness options study shows what the defaults cost. Replacing the flat +6 dB board correction with the soft-ground reflection model is "the largest single step on both days": on the spectral path with station air, rms SEL falls from 2.98 to 1.38 dB on calm day 286 and from 4.18 to 2.43 dB on windy day 284 (options study at NICEOPS `d852613`, on an earlier B407 build). The broadband EAA path with the flat correction sits at 3.40 / 4.34 dB ([harness README](/Users/evg5332/Code/Python/niceops-validation/README.md)). The broadband EAA is linear in range and anchored at 1 km, while true A-weighted absorption is concave in range. With the atmosphere matched, spectral minus broadband therefore runs from −0.75 dB at 0–250 m to +1.76 dB at 1.5–1.75 km.

Both competitors default to a ground model. NORAH2 assumes flat soft ground (σ = 200 kPa·s/m², 1.2 m microphones) unless told otherwise ([Manual Annex B p.48][nman]). AAM applies Chien–Soroka with a flow resistivity from SETUP PARA at a 5 ft grid height ([TR p.75][aam75]). A default NICEOPS footprint is a free-field prediction, which no community receiver experiences.

**Fix:**

- Make `--ground_effect` with a 1.2 m point receiver the default whenever spectra exist.
- Refuse to run spectral-default without an explicit atmosphere (or a station file), rather than assume 20% RH.
- Print the physics actually used on stdout, not only in `log/niceops.log`.
- Keep the goldens by pinning their configuration explicitly in the CTest command lines.

### 2. Source coverage is weakest exactly where the residuals are largest

The open residuals on the day-286 approaches are source-side by elimination:

- **The low, slow end.** States under 100 ft and 20 kt run **+3.9 dB** (energy mean) under stratified rays with the shipped database. On the previous build, stratified gave +4.4 against +3.8 for straight rays, so the excess is not propagation.
- **The fast forward rim.** States at 80 kt or more and 400–2,000 ft, heard 2–10° below the horizon from far ahead, run **−1.5 to −2.8 dB** at 2–6° (−2.3 dB overall beyond 4,000 ft), larger on the lines furthest ahead, and the same under every ray model ([Grazing incidence findings](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md)).

Four builder properties are candidate causes of these residuals; none has yet been tested in isolation:

- **The rim gate zeroes data.** panam's ambient gate *zeroes* bins below the SNR threshold, so rim cells can only read low. The long-range forward error depends on the rim choice. On the earlier straight-ray build (`0decafa`), 660–1000 ft sources read −1.8 dB with its 2° rim against −0.3 dB with a 10° rim, even though the 2° rim won overall. Curved-ray depropagation (`56ce39b`) lifted them to −1.4 dB, and the shipped `5d761b8` grid costs a few tenths more at long range (2–5° at 4,000–8,000 ft: −0.13 dB on `56ce39b`, −0.81 on `5d761b8`) ([forward-line note](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/forward-line-error-source.md); [Grazing incidence findings](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md)).
- **The hover is synthesised.** It comes from the slowest near-level sphere, averaged fore to aft in energy. Against the day-283 HOGE hovers it is 4–13 dB too quiet at 315–1250 Hz to the sides and aft, and 7–19 dB too loud at 63 Hz. A measured correction that brings it within ±1.7 dB in-sample exists but is not shipped.
- **Steep conditions are copies.** The γ = −24° and +35° conditions are copies of level spheres, so steep descents blend toward level-flight directivity. That is the BVI regime where D1.5d's own R22 data show the largest sensitivity ([D1.5d §A.3.1 p.11](https://www.easa.europa.eu/en/downloads/132005/en)).
- **The hover has no variants.** There is no IGE/OGE distinction and no hard/soft-ground variant.

Other tools handle the envelope edges more honestly or more physically:

- **ART** writes 1E35 into unmeasured cells, which AAM refuses to propagate ([TR p.145][aam145]). Its authors explicitly abandoned constant fill as something that "may or may not be acceptable" ([Conner, Burley & Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- **NORAH2** measures hover on a dedicated quarter-circle array 150 m out and carries hover, idle and taxi as their own modes ([Final report §4](https://www.easa.europa.eu/en/downloads/139729/en)). The public V2.0.74 package does not yet realise this: its only hover table is a B412 soft-ground HIGE table filed under R22, which already includes the ground reflection (FREEFIELD = 0). Its idle is a copy of a flyover sphere. Its offsets (+8/−10/−2 dB) differ from D1.5d's (+12/−12/−2.5) *(ours)*.
- **FRAME** fills unmeasured cells from a calibrated physics model, which is what FRAME-QS's load-factor sweep rests on ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- **RNM's** best-documented failure was of exactly this kind: up to 3 dB SEL over-prediction on the advancing side of a CH-146 approach, attributed to extrapolated hemisphere data ([Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)).
- **Korean UAM work** found that three ways of filling null hemisphere cells moved a vertiport impact range between 200 and 700 m ([Lee & Haan 2025](https://dael.euracoustics.org/confs/fa2025/data/articles/000610.pdf)).

**Fix:**

- Have panam write a per-cell, per-band coverage mask (measured, gated-low, filled) into each `.nod` group, and have NICEOPS report, per observer, the share of SEL energy drawn from unmeasured cells.
- Replace the zeroing gate with a censored estimate: keep the gated value as an upper bound and fill from neighbours or a FRAME-type fit.
- Ship the measured hover correction as its own condition group.
- Use FRAME-calibrated spheres, or the 2017 Volpe hover-array data (recorded but not distributed by NASA, per [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf)), in place of the synthetic hover and the γ copies.
- Measure each database against ECAC Doc 32's density guidance (descent angles no more than 3° apart at four speeds; level flight at 0.9 VH, +10, −15 and −30 kt) and log which trajectory states fall outside the measured γ range. For S-76D that is anything steeper than −8.7°.

### 3. Flat ground at z = 0 is the largest physics gap against AAM and NORAH2

**What NICEOPS assumes.** NICEOPS's ground models take the source height as `source.z` and the horizontal ground distance (`src/propagation.hpp:68-81`). An observer's own z moves only the direct path, and propagation tables likewise assume flat, uniform ground.

**What the competitors model.**

- **AAM** classifies each source–receiver terrain slice (flat, concave, convex, wedge or screen) and applies Rasmussen's GTD with per-cell impedance and foliage. It cannot combine this with wind ([TR p.10–12][aam10]).
- **NORAH2** builds a CNOSSOS-EU mean ground plane per path, with screening bounded at 25 dB ([D1.5d §A.4.4–A.4.5 p.24–33](https://www.easa.europa.eu/en/downloads/132005/en)). Against Nord2000 it agrees to ±1 dB generally, but to only ±10 dB where screening exceeds about 20 dB ([EASA D3](https://www.easa.europa.eu/en/downloads/139664/en)).
- **Gaussian beams** over terrain cut straight-ray error "by 5 dB in the illuminated zones and 35 dB in the terrain shadow zones" ([Yunus et al. 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)).
- **sonX** reproduces grazing incidence to about ±1 dB out to 2.8 km with mapped ground cover and buildings ([Wunderli et al. 2022](https://consensus.app/papers/details/9b9419377e9e578d8b9d3ddff58f3bd9/?utm_source=claude_desktop)).

**Why the validation cannot see it.** The Amedee lakebed is flat, so the harness never exercises this gap. Every hillside or heliport study would.

**Fix.** Two stages, cheapest first:

- Add a DEM and a per-path least-squares mean ground plane, feeding the existing ground models with equivalent heights and incidence angle, plus CNOSSOS or Deygout screening per band. This is NORAH2's method, and it slots into `RayPath` without touching the source side.
- Later, extend `--propagation_table` with a terrain-following PE, as APET does with GFPE ([APET](https://software.nasa.gov/software/LAR-17761-1)).

### 4. Reference conditions are unrecorded and mismatches pass silently

Several reference conditions are either not recorded in the database or not checked against it.

**Unrecorded build atmosphere.** EAA is built at 293.15 K, 101.325 kPa and 20% RH, but the `.nod` files record no build atmosphere ([code review X3](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md)).

**Unchecked vehicle file.** Only the tip speed is compared against the database, and only in a log line. A 10% weight error moves the footprint by 0.75 dB median without warning (X2).

**The fixed_load_factor correction.** On `fixed_load_factor` databases, the level correction is 20 lg(C_T/C_T,ref), with C_T computed from the JSON and C_T,ref read from the database (`src/sphere_database.cpp:1221-1243`). This equals 20 lg LF only when the JSON reproduces the database reference. A heavier weight or a hot-day density therefore adds level at *every* state, uncapped. On the C_T-axis databases the clamp is whatever C_T maximum the file holds: 5 g for panam's list. At LF = 10 that gives +13.98 dB on M3 against +20 dB on `fixed_load_factor`. The 20 lg n_z law itself is supported at moderate load factor: FRAME-QS's 3.5 dB at 1.5 g matches it.

**The load-factor law has no directivity.** Every C_T copy above the reference is the reference sphere plus exactly 20 lg(C_T/C_T,ref) in every direction and band (the 0 × copy repeats it). FRAME-QS's support is in OASPL, and Q-SAM without thrust data is limited to load-factor increases under 10% ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)). FRAME-QS's own database stops at 2 g, and Greenwood, Sim & Boyd saw low-frequency loading noise fall 6 dB in some directions as C_W rose 24%. Above about 2 g the copies, and the analytic correction on `fixed_load_factor` databases, are an extrapolation.

**Ground speed, and no altitude scaling.** The spheres are labelled with ground speed. The `speed_reference` attribute exists in NICEOPS but panam does not write it, so every database is ground-referenced. More seriously, the levels carry no ambient-pressure or hover-tip-Mach scaling. At constant IAS and weight, the AS350 rose 4.7–5.4 dB in thickness noise from sea level to 7000 ft. Data collapse only when C_W and M_He are matched and amplitude is scaled with ambient pressure; "no simple empirical model can be used to 'correct' the data" ([Watts, Greenwood & Stephenson 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf); [Greenwood, Sim & Boyd](https://ntrs.nasa.gov/api/citations/20160009056/downloads/20160009056.pdf)). NICEOPS's μ and C_T are non-dimensional, but its tip speed is a database constant. Running the Amedee B407 database (about 4000 ft) at a sea-level heliport is therefore an unflagged extrapolation. Scaling the AS350 result linearly gives a rough size: about 3 dB in in-plane thickness noise between 4000 ft and sea level (our interpolation), with BVI able to move several dB either way. The 20 lg(p_amb/p_ref) correction in the fix is valid only once C_W and M_He are matched, so it should be applied with an M_He check, not on its own.

**What NORAH2 records.** Its hemispheres carry TAMB, RELHUM and PAMB, IAS (ACSPEED), γ, pitch, roll, rotor RPM and the wind components in every file ([D1.5d App. A p.34–35](https://www.easa.europa.eu/en/downloads/132005/en)). Its RPM field is corrupt in about 15 shipped files *(ours)*.

**Fix:**

- Have panam write the build atmosphere, the speed reference, the rotor RPM and the site's pressure altitude as root attributes.
- Have NICEOPS refuse a vehicle file whose W, R, V_tip or ρ reproduces C_T,ref worse than about 2%, unless overridden.
- Add 20 lg(p_amb/p_ref) amplitude scaling and an M_He axis or correction, and say in the output which database conditions were extrapolated.
- Cap the analytic load-factor correction at the largest load factor the database's source supports (2 g for FRAME-type data), say so in the log, and replace the copies with FRAME-calibrated C_T sweeps where available.

### 5. One metric, a non-standard "DNL", and a window competitors do not use by default

**SEL is the only metric.** The only per-observer metric is A-weighted SEL. A grep for LAmax, EPNL, PNL, NOY, tone correction, loudness, C-weighting, time-above and Leq finds only a comment.

**`-D` is not a DNL.** It sums night-penalised exposures without subtracting 10 lg 86 400 = 49.37 dB (`src/util/sel.hpp:175-183`; `src/niceops.cpp:2874-2881`), and each operation contributes only its 10 dB-down window. The README says this outright: "With one daytime operation the DNL equals the SEL exactly". Calling it DNL is a naming defect.

**`-m` averages in the wrong domain.** It writes the arithmetic mean of per-observer SEL in dB, not an energy mean.

**What the competitors report.** AAM reports Lmax, SEL (A and C), DNL, CNEL, NEF, PNL/PNLT/EPNL, WECPNL, NA, TA and TVL ([TR p.50–52][aam50]). NORAH2 reports LAE, LAmax, PNLTM, EPNL and NAx, with user-defined Lden/Ldn weightings ([Manual §5.1.3 p.21–22][nman]).

**The window differs too.** Both competitors integrate the **whole** history by default: DBDOWN10 is opt-in in AAM ([TR p.64][aam64]) and ISEL10DB = 0 is NORAH2's land-use default. A like-for-like comparison therefore needs a flag on one side. Even with it on, AAM's documented rule is a threshold on samples ("only include points within 10 dB of the maximum"), not NICEOPS's contiguous first-to-last interval, so dips are excluded there and included here. NORAH2 does not state how ISEL10DB = 1 defines its period. SELENE's definition matches NICEOPS's. The size of the full-history versus 10 dB-down difference on NICEOPS footprints has not been measured.

**Fix:**

- Compute LA(t) on the reception-time axis that `reception_factor` already defines, and emit LAmax, a full-event SEL option and NAx.
- For a like-for-like comparison, add an option that integrates only the samples above peak − 10 dB, and another that integrates the full history.
- Implement EPNL from the per-band spectral path. The band set reaches 10 Hz–10 kHz, which covers the 24 PNL bands. Before trusting it, report per band and cell how much of the sphere above 2.5 kHz survived the 10 dB ambient gate. The 2017 spheres are SNR-limited there, and PNL and the tone correction are weighted toward 1–5 kHz. NORAH2 repairs masked high bands by last-good-band extrapolation (ICAO Doc 9501); panam would need an equivalent before EPNL is meaningful.
- Rename `-D` to "night-weighted cumulative SEL" or normalise it by 86 400 s, and make `-m` an energy mean.

### 6. A single compact source, and a Doppler flag that is never read

NICEOPS places one whole-vehicle sphere at the track point (`src/niceops.cpp:778`). The JSON's `Tail Rotor` block is never read. Measured empirical spheres fold the components together correctly at the measured conditions. The single sphere fails in three cases:

- when components scale differently off-database;
- when tones must be summed coherently;
- when a design change is wanted.

The evidence for why this matters:

- Main rotor alone "will produce results that severely underpredict integrated metrics" ([Boyd et al. 2017](https://ntrs.nasa.gov/citations/20170002590)).
- The Bell 206B3 main and tail rotors have strongly different directivities ([Greenwood & Schmitz 2014](https://pure.psu.edu/en/publications/separation-of-main-and-tail-rotor-noise-from-ground-based-acousti/)).
- Vold–Kalman separation now splits main rotor, tail rotor and broadband residual into separate hemispheres ([Rachaprolu & Greenwood 2024](https://rosap.ntl.bts.gov/view/dot/72623/dot_72623_DS1.pdf)).

AAM carries up to ten broadband, ten narrowband and ten pure-tone spheres at body offsets. Its tones are summed coherently with absolute phase ([TR p.30][aam30], [p.40–41][aam40]). CAROT/MOTUS separates main rotor, tail rotor and residual, and applies its design masks per component.

The Doppler convention is right for panam's spheres. They keep the received shift (`DOPPLER_SHIFT_REMOVED = 0` in the per-run `.nc` files; the assembled `.nod` carries no flag), so NICEOPS correctly applies no shift. But NICEOPS never *reads* the flag (`src/sphere_database.cpp` writes it only in the AAM exporter). An imported AAM or ART sphere marked de-Dopplerized would therefore be propagated with stationary-frame frequencies through absorption, ground and A-weighting.

There is a trap here too. AAM's DOPPLERSTATE keyword (1 = moving) is coded inversely to the NetCDF variable `DOPPLER_SHIFT_REMOVED` (1, by its name, = stationary), and Bullwink's input (0 = moving, 1 = stationary) follows the file's sense, not the keyword's ([TR p.11][aam10], [p.165][aam165]). Blending reception-referenced spheres across speed also blends different shifts with no correction. AAM instead applies a test-to-model speed difference ([TR p.11][aam10]). When conventions are applied consistently, de-Dopplerized processing and moving-observer prediction agree "typically within 0.5 dB or less" ([ASCENT 038 2025](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-038-2025-Annual-Report.pdf)).

A misread flag costs more. By our arithmetic, at 100 kt (M ≈ 0.15) the fore/aft shift f/(1 − M cos θ) is about 0.7 of a third-octave band. Where A-weighting falls 3 dB per band (63–125 Hz, the main-rotor harmonics), mis-shifting moves a tonal band's A-weighted contribution by about 2 dB; above 500 Hz it is under 1 dB. panam can already build a de-Dopplerized hover (`hover_source` with `remove_doppler=True`) beside received-frame flight spheres, and the `.nod` records the convention nowhere.

**Fix:**

- Have panam write `DOPPLER_SHIFT_REMOVED` per `.nod` group, and have NICEOPS read it and refuse mixed or set flags. Better, apply f′ = f/(1 − M_r) per band before absorption, ground and weighting.
- Allow several spheres per condition group, each with a body-frame offset and a type (broadband or tonal), summed incoherently at first. This is the hook for MR/TR separation and for per-rotor eVTOL spheres. The dimension-agnostic `Delaunay` already makes a fourth index axis, such as nacelle tilt, structurally cheap.

### 7. No static operations, and a hard switch to the hover nose

**No static operations.** Hover must be flown as trajectory states, so a 30 s HIGE at the pad needs a synthetic track. HNM/AEDT, NORAH2 and AAM all treat static operations as first-class:

- AEDT has 8 static and vertical modes ([ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)).
- NORAH2 has IGH, OGH, FIDL, RIDL and TAX ([Manual §5.1.6 p.26][nman]).

**A hard switch at 5 kt.** The sphere nose is the (air) velocity projected into the TPP (`src/sphere_database.cpp:678-688`, fed from `src/niceops.cpp:763-782`). Below 5 kt it is the INS heading instead (`nose_direction`, `src/trajectory.cpp:828-835`), turned into the local frame by an offset inferred from the track's states over 20 kt, or by `--frame_bearing` (`heading_offset`, `:799-826`). A track with neither falls back to the drift, with a warning. Before `9fc714a` every hover nose followed the drift, a few tenths of a knot in any direction. The switch is hard: the code review found that a nose rotated by 30–90° between the two branches steps the level by 2.2–4.4 dB (X4, library-level; how often real hovers and approach ends present such angles is unmeasured). NORAH2's hover tables are azimuth-only, and HNM's static directivity is an eight-term Fourier series in azimuth. NICEOPS's full (φ, θ) hover sphere is richer in kind; it is let down only by the data behind it and by the hard switch at 5 kt.

**Fix:**

- Add a static-operation record (position, heading, mode, duration) that expands to stationary states internally.
- Blend the nose from velocity to heading over about 3–8 kt with a smooth weight, rather than switching.

### 8. The upper hemisphere is stored but never propagated

**The fold.** `phi_degrees` computes atan(across/|up|) (`src/sphere_database.cpp:711-723`). Every propagation lookup goes through it, so a direction above the TPP reads its mirror image below. For panam databases this loses nothing: panam fills the upper half by reflection through the horizontal plane, which in steady flight is the same reflection, exact to 0.0000 dB over 703 cells.

**Where it would matter.** The fold matters wherever genuine upper-hemisphere data could:

- **Receivers above the source.** Scattered observers carry their own z, and the default point-receiver height is 1.2 m.
- **TPP tilt in maneuvers.** A 0.1 g deceleration tilts the normal about 5.7° aft, so forward observers at depressions under about 5.7°, the lower half of the 2–10° band where the validation residuals concentrate, sit above the TPP.
- **Refracted rays.** Under downward refraction, long-range rays are launched above the horizontal.

**What the literature says.** FRAME OASPL spheres differ above and below the horizon, and "models which assume noise radiation to be symmetric about the rotor tip-path-plane may lead to inaccurate" results ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)). AAM accepts full spheres. NORAH2 clamps to the equator *(ours)*. Elevated crane microphones at Eglin did measure above the TPP ([Conner, Burley & Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).

**Fix.** Use φ = atan2(across, −up) over the full ±180° and read the stored upper half. This is a no-op on every mirrored database and enables FRAME- or crane-derived upper data. Keep `repair_legacy_hover_directivity` as is.

### 9. Silent numerical failure paths

The code review reproduced several defects. All of them leave the exit status at 0 or put the explanation only in the log:

- **NaN and garbage inputs.** NaN or −inf results come from a flat time column, a NaN coordinate or a unit-slipped atmosphere, with 4–10 dB shifts and exit status 0. Errors appear only in the log.
- **Thread exceptions.** A worker exception at `-j ≥ 2` calls `std::terminate`.
- **Inward spreading.** Spreading is unbounded inward: 20 lg(r/R) has no floor, and inside the sphere EAA·(r − R) adds level. That affects 28 demo observers by up to 4.4 dB.
- **Coherent ambient.** The `-a` ambient is mixed in *coherently*, adding +6.02 dB at equal levels.
- **Bounds near the nose axis.** `DirectionBounds` is exceeded by up to 0.7 dB within 1° of the nose axis, so the windowed path is not strictly exact there.
- **Spline end conditions.** Natural end conditions force zero acceleration at geodetic track ends. That is worth up to 1.8 dB, and 7.5–11° of TPP-normal error, over the last seconds.
- **Asserts.** The assert guards are compiled out in the documented Release build (below).
- **`--db_cutoff` is not exact per cell.** 157 of 592 cells at or above 85 dB change, by up to 0.43 dB, although the README says every such cell is "unchanged to the last bit"; only the per-peak guarantee holds.
- **Multi-track outputs.** `-i` and `--band_history` overwrite each other across `-t` tracks, a metric geodetic track loses its frame, and DNL on a geodetic track with no clock comes out +10 dB.

Most of these are fixed on the unmerged `review-fixes-2026-10-02` branch ([code review](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md)). The fix is to merge that branch, make ambient addition incoherent, and add NaN guards to `golden_compare`.

### 10. Documentation contradicts the code in seven verified places (plus two unverified)

The reader of CLAUDE.md is the next maintainer, so stale load-bearing claims there are a defect in their own right. The conflicts, the last two taken from the code review and not re-verified here:

| Topic | CLAUDE.md / README says | Code at e472b6f does |
| --- | --- | --- |
| Spectral default | "`--spectral` is opt-in" (CLAUDE.md); opt-in, and broadband is the "default" column (README :168, :286) | Defaults to spectral and falls back to broadband when the database has no spectra (`niceops.cpp:546-552, 2357-2361`) |
| SEL window | `peak_range` "walks out from the peak and stops at the **first** sample below peak − 10". Names `expand_window()`. Cites the "~0.7% of observers" fragility | First-to-last crossing, dips included (`util/sel.hpp:22-60`). `expand_window` no longer exists; the fragility argument describes a superseded rule |
| Windowed path with spectral | Skipped under `--spectral`; `--db_cutoff` refused with `--spectral` | Both run on spectral (`niceops.cpp:1499-1500, 1369-1372, 1446-1453`) |
| Spectral databases | "No database in `databases/` carries spectra" | Six do: the AS350, B407, Be206L3, EC130, RO-44 and RO-66 `_fixed_load_factor.nod` files each carry 31 bands |
| Database tracking | README :257-258: "Databases are not tracked in this repository" | `git ls-files databases` lists 24 files |
| Asserts | "Never add `-DNDEBUG` to Release"; asserts live as a port contract | CMake's default `CMAKE_CXX_FLAGS_RELEASE` is `-O3 -DNDEBUG`, and nothing removes it: `CMakeLists.txt:16-24` defaults to Release and has only a comment saying not to add it |
| Projection | Transverse Mercator (`geodetic.hpp:8`, CLAUDE.md, README) | Ellipsoidal American polyconic (`geodetic.cpp:68-85`) |
| Evaluated fraction | "70% of states evaluated", "1.2x faster" | About 94% evaluated per the code review (not re-verified here) |
| `--db_cutoff` exactness | README: every cell at or above the cutoff "unchanged to the last bit" | 157 of 592 cells at or above 85 dB change, by up to 0.43 dB, per the code review (not re-verified here) |

Sources: [verification_niceops.md §§1, 3, 6, A–E](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/verification_niceops.md); [niceops_inventory.md §6](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/niceops_inventory.md).

Three stale code strings sit alongside: `niceops.cpp:2699, 2705` still name "--no-spectral or --db_cutoff" as broadband routes, and the `--propagation_table` help (`niceops.cpp:2467-2472`) still says "not with --refraction or --plate_table", although since `e472b6f` a table runs beside either ray model and plates read it. Our own inventory note also gave the load-factor law as "+20 log10 LF"; the code applies 20 lg(C_T/C_T,ref) (item 4).

**Fix.** Rewrite the "SEL integral is oversampled" and "Walking out from the peak" sections of CLAUDE.md for the first-to-last rule, and drop the decimation table's first-crossing column, since its outermost-crossing column already measures `ten_db_down_interval`. Then remove the `NDEBUG` flag explicitly in CMakeLists.txt (or retract the claim) and correct the rest of the table.

### 11. Validation is the strongest in the field, but narrow and without error bars

**What exists.** The harness holds the only multi-run, per-observer SEL error distribution any sphere tool has, though it is internal and not yet published. Every other quantified sphere result predicts the run that built the sphere: the XV-15 RNM case ("within 0.5 dB", [Lucas & Marcolini](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)) and the R44 table in ASCENT 038 (within 1.2 dB, [ASCENT 038 2018](https://ascent.aero/documents/2019/07/ascent-project-038-2018-annual-report.pdf)).

Beyond the B407 headline:

- **A second aircraft.** The B206L3 scores 1.50 / 1.49 dB rms SEL (station / stratified) with the A cards that went into its database excluded.
- **Six aircraft, unsummarised.** The harness also tracks a six-aircraft set of 139 approach runs (`approach/`), not yet summarised per aircraft.
- **Footprints.** Four measured footprints come out at +0.5/1.3, 0.0/1.9, +0.7/1.7 and −0.8/1.3 dB (median/rms, stratified + plate; runs 286270, 286277, 284205, 284211). The flat +6 dB correction over-predicts their sidelines by 3–8 dB.

Sources: [niceops_inventory.md §10](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/niceops_inventory.md); [harness README](/Users/evg5332/Code/Python/niceops-validation/README.md).

**Its limits.**

- **One aircraft and one site dominate.** The headline rests on B407 approaches at one flat site, where the approaches are held out of the sphere build (steady runs only) but share aircraft, site and campaign.
- **Contamination.** Twelve Be206L3 approach cards were built into the shipped Be206L3 database, worth about 0.07–0.1 dB rms.
- **No LAmax or turn validation.** No LAmax validation exists, and the 2017 turns and accelerating turns are unused.

The turns are large effects. Measured LA,max differences between paired turns reach 7.5–13.3 dBA ([Stephenson et al. 2020](https://ntrs.nasa.gov/api/citations/20205003220/downloads/VFS76_Paper_v7.pdf)), and those differences are exactly what Q-SAM with C_T claims to capture.

**Error bars.** No hemisphere footprint tool has propagated an error bar. The three ingredients each have published magnitudes:

- source repeatability, ±0.5 dBA in guided level flight ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)) to 7.5 dB at a single XV-15 microphone ([Boyd 2003](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf));
- database processing, about 2 dB of refraction contamination near the rim under uniform-atmosphere depropagation ([Munt et al. 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full));
- weather, up to 4 dBA in the illuminated zone ([Fuerkaiti et al. 2021](https://consensus.app/papers/details/b05b0c351cd55240bdf8f87bd18ec91f/?utm_source=claude_desktop)).

**Fix:**

- Score the 2017 turn and accelerating-turn pairs on LAmax and SEL.
- Publish per-aircraft statistics for the existing six-aircraft, 139-run approach set, then rebuild one aircraft's database with every approach card excluded.
- Use the repeat runs plus the measured profiles to put a Monte-Carlo band on a footprint. That would be a first for the field.

**No cross-tool calibration yet.** No published head-to-head comparison of hemispheres for the same vehicle and condition exists ([hemisphere_generation.md §2](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/hemisphere_generation.md)), and NICEOPS has the material for one:

- Compare panam's B407 sphere for run 109 with the AAM tutorial's ART sphere for the same run ([TR p.199–207][aam199]), and the shipped 2017 databases with the Watts et al. ART hemispheres, cell by cell over 30–2500 Hz. This isolates builder differences: 10 dB gate with subtraction against a 5 dB gate without, plate BEM against −6 dB, curved against straight depropagation.
- Load NORAH2's R44, R66 and AS350 hemispheres with `load_norah2_hemisphere`, undo their 60 m ICAO absorption and rescale them to 30.48 m. Compare them with RO-44, RO-66 and AS350 over the |φ| ≤ 60° cells both cover, and fly one NICEOPS track with each database.
- Following ANOPP2's ASOFT, reduce a NICEOPS sphere to the AEDT 0/±45° slices and fly it, to quantify what AEDT-style directivity costs on the day-286 approaches.

## Where NICEOPS is ahead of the state of the art

**State indexing and orientation.** NICEOPS's state handling is what the maneuver literature converged on, and the regulatory tools have not adopted it:

- **Indexing.** (μ, γ_e, C_T), with C_T taken from the trajectory's own load factor and the sphere oriented by the inertial TPP normal, is Q-SAM2 plus FRAME-QS's load-factor axis. It uses the same Delaunay-plus-nearest-fallback design as SELENE, built here on exact predicates. AAM needs the FRAME database type to reach three axes, and NORAH2 stays at (V, γ). The `--no-qsam` ablation alone moves the demo by 0.10 dB median and 1.7 dB at worst.
- **Orientation.** Because the bank comes from path curvature, NICEOPS cannot reproduce NORAH2's turn-sense sign error.
- **Air reference.** Its `--speed_reference air` machinery anticipates the crosswind result that BVI is unchanged in crabbed flight but changes by about 11 dB in sideslip ([Greenwood & Sim](https://ntrs.nasa.gov/api/citations/20130013862/downloads/20130013862.pdf)). No other tool documents an air-referenced database convention; AAM and NORAH2 both say they neglect wind in lookup.

**Propagation.** No footprint tool found combines all of the following, each with reference checks:

- per-band ISO 9613-1 absorption through layered profiles;
- a Chien–Soroka ground pinned to panam at 1e-9 dB with four impedance models;
- a BEM ground-plate receiver;
- equivalent-linear and layered-eigenray refraction, with every arrival in a duct summed;
- residue-series shadow transitions and creeping-wave loss;
- a split-step Padé PE within 0.3 dB of the Attenborough et al. benchmark exact solutions, used to build propagation tables.

The comparators fall short in specific ways:

- SELENE floors the shadow wherever the eigenray search fails.
- AAM's refraction is an omnidirectional pre-tabulation over flat ground, and cannot be combined with terrain.
- NORAH2 has no refraction at all.
- APET's GFPE tables are the only other full-wave coupling found, and nothing found uses them as a validation oracle.
- The TU Delft/Arntzen model (Gennaretti et al. 2017) fills shadows with ground-wave, diffraction and turbulence-scattered energy, the counterpart of NICEOPS's ShadowTransition plus `--scattering`; no accuracy figure is published for it ([Gennaretti et al., IJA 2017](https://repository.tudelft.nl/file/File_8a12a003-8adb-4fca-b37d-03c081592697)).

The measured payoff on the windy day (held-out database, harness runs that predate the `5d761b8` rebuild): the 5–10° bin went from +4.5/5.8 dB under straight rays to +0.8/3.0 under stratified rays, and overall rms from 2.36 to 1.70 ([harness README](/Users/evg5332/Code/Python/niceops-validation/README.md)). On co-located board/pole pairs at 0–2° (4,861 frames), straight rays are 7.6 dB rms against 0.8 for linear and 1.2 for stratified ([Grazing incidence findings](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md)).

Its known limits (the README and the code review, re-checked at `5d761b8`; the sample count and caps are unchanged at `e472b6f`): the stratified 32-sample arrival search is unconverged in ducts (4.3% of pairs over 1 dB against a 2048-sample reference, about 5% against 128 per the README; up to 10.5 dB in single cells). Per-ray focusing caps let the summed gain reach 17 dB. An elevated duct turning just above the ground (run 286288) comes out 11–21 dB under the PE by default. On the measured day-286 frames where it falls short, the opt-in Airy tails (`--tunnelling`, `--skip_fill`) take the error from −10.5 and −15.5 dB to −2.1 and −5.0 dB at 50–160 and 200–630 Hz, but barely touch 800–2000 Hz (−20.6 to −18.9 dB). Over all far frames, though, the two options over-predict at 8–14 kft and 0.5–2.5° by 6–7 dB, where plain stratified is within 0–4 dB (harness, uncommitted), so the deficit-frame score alone does not make them a default. The shadow-boundary law cannot get below about 6 dB rms pooled at the edge, and linear refraction can be 10 dB off the PE near the ground ([niceops_inventory.md §11](../research_notes/Hemisphere%20based%20rotorcraft%20noise%20models/niceops_inventory.md)). Rim cells were themselves depropagated through the measured atmosphere, so scoring long-range refraction against the same campaign is partly circular ([Munt et al. 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full)). panam depropagates to 2000 ft (610 m), past Munt's recommended r < 450 m.

**Database building.** panam's pipeline addresses three contamination mechanisms that the literature names but other builders leave alone:

- **Refraction baked into the sphere.** Munt et al. estimate it at about 2 dB near the rotor plane under uniform-atmosphere depropagation; panam depropagates along curved rays where profiles exist.
- **Ground-board response.** panam divides out the plate's BEM response instead of a flat −6 dB. ARP 4055 is specified only within 30° of normal incidence, which every rim cell violates ([SAE ARP4055](https://saemobilus.sae.org/standards/arp4055-ground-plane-microphone-configuration-propeller-driven-light-aircraft-noise-measurement)).
- **Noise floor.** panam subtracts the ambient instead of only gating it. AAM's default gating biases near-floor levels high by up to about 1.2 dB at its 5 dB threshold, by our arithmetic: 10 lg(1 + 10^−0.5).

**Time resolution, kinematics and speed.** NICEOPS integrates at the trajectory's 0.02 s, against RNM's 0.5 s and AAM's 1–2 s, and applies a reception-time Jacobian worth 0.4 dB on the check set. The demo case's 94.3 million observer-state pairs go through the acoustic loop in about 0.16 s of wall time on 14 threads (about 1 s for the whole run, most of it database load), which makes the propagation options affordable at footprint scale.

**Validation scale.** As item 11 says, its 1.26 dB rms over 1,112 held-out approach microphone-runs is unmatched in the literature, and would be the field's first such distribution once written up.

## Conclusion

The comparison changes where effort pays off. The design arguments that dominate the hemisphere literature are already settled in NICEOPS's favour or are second-order:

- dB versus energy in the condition blend is worth 0.25 dB;
- barycentric versus natural-neighbour weights, and (V, γ) versus (μ, C_T, α_TPP) indexing, are questions NICEOPS already answers the recommended way;
- no other published propagation chain combines what NICEOPS's does, though its own limits in ducts and at shadow edges are documented.

What separates NICEOPS from being the reference tool is mostly outside the code: the sphere database at its edges, and the gap between the validated and the default configuration. The two largest measured residuals, the +3.9 dB slow end and the −1.5 to −2.8 dB fast forward rim, are both source-side by elimination: neither goes away under any ray model. Their causes are not yet isolated. Section 2 names the builder candidates: the zeroing gate, the synthesised hover, the γ copies and the hover's lack of variants. The database blend is ruled out. Day 286's spheres are +0.5 to +1.2 dB louder at the forward rim, but a day-286-only database moves only the 6–10° rim and leaves the 2–6° shortfall in place. Each can be tested cheaply: rebuild with the measured hover correction and rescore the under-100 ft, under-20 kt frames; rebuild with a censored rather than zeroing gate and rescore the 2–10° forward frames. Until those runs exist, the coverage mask is a diagnostic, not a fix.

On the regulatory side, ECAC's retreat from interpolation shows the field lacks exactly the evidence NICEOPS can now supply: a measured interpolation error (about 1.1–1.8 dB rms in leave-one-out), a quantified domain effect, and a held-out multi-run validation. Written up with LAmax and terrain added, that evidence would make NICEOPS the strongest available argument for interpolating hemispheres at all.

[aamAnn]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=3
[aam3]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=19
[aam10]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=26
[aam12]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=28
[aam13]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=29
[aam14]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=30
[aam20]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=36
[aam24]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=40
[aam29]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=45
[aam30]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=46
[aam31]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=47
[aam33]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=49
[aam39]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=55
[aam40]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=56
[aam41]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=57
[aam50]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=66
[aam64]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=80
[aam75]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=91
[aam144]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=160
[aam145]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=161
[aam161]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=177
[aam163]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=179
[aam165]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=181
[aam199]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=215
[aam204]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=220
[aam214]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=230
[aam226]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=242
[nman]: <file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>
