# NORAH2 (NOise of Rotorcraft Assessed by a Hemisphere-approach, version 2)

Scope note for the report writer. The primary source is the public **NORAH2 User Manual**, EASA report SC03-D1.2, by A. Diez, H. Olsen and L. Meliveo, report date 31 May 2023, PDF built 11 Nov 2024, shipped in the `NORAH2_V2.0.74_public` package. The manual is a software and file-format manual. It does not describe the physics: it defers that to the methodology reports (its refs [2] and [3], p.44). To answer the method questions, three public EASA documents were also read:
- **D1.5d** "Rotorcraft noise modelling guidance", Olsen, Tuinstra and van Oosten, 31 Jan 2024, contract EASA.2020.FC.06. This is the method specification that NORAH2 implements.
- The **Final report** "NOISE: Final report", Tuinstra, Olsen and van Oosten, NLR, April 2024.
- The **webinar** slides of the final dissemination event, 13 Mar 2024.

The shipped data and example runs in the package were also inspected: 204 `.hem` files, 11 `*_Triangulation.int` files, and the `ARP/Case1-4` inputs and outputs. Executables were not run. Some conventions were reverse-engineered from the shipped `.his` outputs; these are marked as inferences.

Source keys used in links:
- **Manual**: the local PDF. Page numbers are the printed "PAGE n" footers.
- **D1.5d**, **Final report**, **Webinar**: EASA downloads. Page numbers are printed footers, or slide numbers for the webinar.
- **Package**: local files under `NORAH2_V2.0.74_public/`.

## 1. Name, developers, history, intended use

### Takeaway
NORAH means "NOise of Rotorcraft Assessed by a Hemisphere-approach". It is a **European** model: built for the EU (DG MOVE, then EASA) by NLR (Netherlands) with SINTEF (Norway) and Anotec Engineering (Spain). It is not part of the NASA, Volpe, FAA, RNM, HNM or AEDT lineage. Its purpose is to add rotorcraft to strategic noise mapping and land-use planning (the END, ECAC Doc 29). NORAH1 (2017) used fixed hemispheres only. NORAH2 (2023-24) added four things:
- interpolation between hemispheres;
- hover, idle, taxi and turns;
- topography with CNOSSOS-EU-style screening;
- radar-trajectory input.

### Cited Findings
- The acronym is defined: NORAH1/NORAH2 = "NOise of Rotorcraft Assessed by a Hemisphere-approach version 1/2" — [Manual, Abbreviations p.7](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Customer, contract and authors.**
  - The document is EASA report SC03-D1.2 under contract EASA.2020.FC06.SC03, by Anja Diez, Herold Olsen and Luis Meliveo, dated 31 May 2023 — [Manual, cover p.2-3](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
  - It describes itself as the "functional requirements and architecture for the NORAH2 prototype software" — [Manual, Summary p.4](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **History.**
  - NORAH1 was developed under the EASA-funded project MOVE-C2-2014-269 and released in 2017, by NLR, SINTEF and Anotec.
  - In 2020 EASA launched a 4-year framework contract, "Environmental Research – Rotorcraft Noise", to extend NORAH1.
  - The updated method was presented to ECAC AIRMOD, and documented in "D1.5 Draft rotorcraft noise modelling guidance" — [Manual §1 p.8](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- The four NORAH1-to-NORAH2 extensions, as the manual lists them:
  1. arbitrary trajectories via interpolation between hemispheres of different flight conditions;
  2. turns, hover, idle and taxi;
  3. a ground-absorption model for varying topography;
  4. screening by topography.

  — [Manual §1 p.8](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- The software is self-described as a prototype "not optimized for efficiency", with a command-line interface and ASCII inputs and outputs — [Manual §1 p.8](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Rationale.** Existing guidance does not cover rotorcraft:
  - ICAO Doc 9911, ECAC Doc 29 and Directive 2015/996 are fixed-wing only, although helicopters fall within the END.
  - DG MOVE tendered MOVE-C2-2014-269 (2015–2017) to fill that gap.
  - That contract measured eight small and medium types and produced the NORAH prototype.

  — [Final report, Summary p.3](https://www.easa.europa.eu/en/downloads/139729/en).
- The method is meant to complement ECAC Doc 29, "suited to support strategic noise mapping activities" — [D1.5d §A p.6](https://www.easa.europa.eu/en/downloads/132005/en).
- **Comparison with AEDT.** The webinar criticises AEDT's helicopter treatment on three counts:
  - its types do not match the European fleet;
  - its noise-distance data have "limited directivity consideration";
  - its noise data are independent of climb/descent angle or speed.

  — [Webinar slide 9](https://www.easa.europa.eu/en/downloads/139577/en).
- **Related models.** NORAH sits alongside HELENA (Airbus/FP6) and AAM (Page et al. 2009) as "next-generation" hemisphere models. The hemisphere file format is described as "used in the NORAH and HELENA" models — [Final report §2.1.2 p.11](https://www.easa.europa.eu/en/downloads/139729/en); [D1.5d App. A p.34](https://www.easa.europa.eu/en/downloads/132005/en).
- NLR describes NORAH2 as led by an NLR consortium with Anotec and SINTEF supporting EASA from 2020. It describes the hemisphere approach as suitable for drones and UAM as well — [NLR news](https://nlr.org/news/the-norah2-project-helps-predict-noise-nuisance-of-rotorcraft).
- **Versions present in the package.**
  - Folder name: `NORAH2_V2.0.74_public`.
  - `NORAH2.exe`: 230 MB, dated Jun 2024. `SingleEvt.exe`: 1.7 MB, dated Jun 2024.
  - The example logs print "SingleEvt V2.0 $Rev: 60 $" and "NORAH V2.0 $Rev: 63 $".
  - The manual's example outputs show "$Rev: 75 $" and "V0.2 $Rev: 51 $".

  — [Package bin/, ARP/Case2/log.txt](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/>); [Manual p.35-36](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Implementation.**
  - The single-event module is Fortran. The multi-event module is Python, run from the command line.
  - Both are Windows executables (XP to 10).

  — [Manual §1.1 p.8-9; §6.1 p.39](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
  - The Final report says Fortran 95. It also says the software is accompanied by training videos on the EASA website — [Final report §2.3 p.15; §5 p.34](https://www.easa.europa.eu/en/downloads/139729/en).

### Inferences
- The task's premise of a "NASA/Volpe/FAA ... successor lineage" is wrong. Nothing in these documents links NORAH to RNM, HNM or AEDT, other than AEDT being named as a model with inadequate helicopter directivity.
- The intellectual lineage is HELENA plus ICAO Annex 16 Chapter 8 measurement practice.
- `RunCase.bat` has commented calls to `src\py\MultiEvt_py3x.py` and `src\IVF\Debug\SingleEvt.exe`, and no `src/` ships. This suggests two things:
  - the 230 MB `NORAH2.exe` is a frozen Python interpreter bundle;
  - the Fortran is Intel Visual Fortran.

  — [Package ARP/Case1/RunCase.bat](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case1/RunCase.bat>).
- The manual is dated May 2023 but describes the 2024 package. It contains stale cross-references: "chapter 6" for file formats (actually chapter 5), and "App B.1" for sigma (actually Annex D).

### Gaps
- No public source code was found, so the implementation can only be checked against the specification and the example outputs.
- The full D3.2 (2017) NORAH1 method report and the NORAH 2.0 validation report (Olsen et al. 2024, SC03 D3) were not obtained.
- Version history between Rev 51 and Rev 75 is not documented.

## 2. Source hemisphere format: geometry, radius, angular grid, spectra, metadata

### Takeaway
- **Geometry and spectra.** A forward-flight hemisphere is a 19 × 19 table on a 10° grid: polar angle θ from 0 to 180°, and lateral (azimuth) angle φ from −90 to +90°. Each cell holds 31 one-third-octave bands, 10 Hz to 10 kHz, as **SPL at a 60 m reference distance**. The stored levels include atmospheric absorption to 60 m at ICAO reference conditions.
- **Coverage.** Only the **lower** hemisphere is stored (φ = 0 is straight down). Data are typically present only for φ within ±60° and for θ within the measured 10-dB-down span; empty bins are −999.
- **Hover.** Hover hemispheres use a different, 1-D format: directivity versus azimuth about the vertical only.
- **Metadata.** Each file records the measured flight condition: IAS, flight-path angle, pitch, roll, rotor RPM, wind, and atmosphere.

### Cited Findings
- **D1.5d definition of the hemisphere.**
  - Levels are at a fixed reference distance of **60 m** and include atmospheric absorption under ICAO reference conditions (101325 Pa, 298.15 K, 70% RH).
  - Bands are one-third octave, from band 10 (10 Hz) to band 40 (10 kHz).
  - The grid is binned in **10°** intervals of azimuth φ and polar angle θ.

  — [D1.5d §A.3.2 p.14](https://www.easa.europa.eu/en/downloads/132005/en).
- **Body-axis definition** (D1.5d eq. 11):
  - x = r_h cos θ; y = r_h sin θ sin φ; z = r_h sin θ cos φ;
  - −90° ≤ φ ≤ 90°, 0° ≤ θ ≤ 180°;
  - negative φ is port, positive φ is starboard;
  - θ < 90° is forward emission, θ > 90° is rearward.

  — [D1.5d §A.3.2 p.14-15](https://www.easa.europa.eu/en/downloads/132005/en).
- **Required coverage** (eq. 12): −60° ≤ φ ≤ 60°, and θ_t1 ≤ θ ≤ θ_t2, the polar angles at the 10-dB-down instants. Larger angles "require complex measurement setups" — [D1.5d §A.3.2 p.15](https://www.easa.europa.eu/en/downloads/132005/en).
- **File format** (shared with HELENA):
  - a title line;
  - N table constants: POLDIST (reference distance), FREEFIELD (content flag), NOVALUE (−999), and TAMB/RELHUM/PAMB (the hemisphere's atmospheric conditions);
  - measurement metadata averaged over runs: Tm, RHm, Pm at 10 m, RmOmega (rotor RPM), ACSPEED (indicated airspeed, kt), GAMM (path angle), PITCH, ROLL, TW/CW/HW (total, cross and head wind, kt);
  - axis definitions: THETAOBSAC with 19 values and PHIOBSAC with 19 values;
  - NPARAD and NFREQ = 31 with the band centres;
  - then one block per PHIOBSAC, with rows of θ by band.

  — [D1.5d App. A p.34-35](https://www.easa.europa.eu/en/downloads/132005/en).
- **Shipped database contents.**
  - 204 `.hem` files for 11 types: A109 (21), AS350 (17), B412 (15), Cabri (17), EC120 (28), EC135 (19), R22 (21), R44 (17), R66 (20), S300 (13) and S92 (16).
  - 203 of them use POLDIST 60, FREEFIELD 2 ("Atmospheric absorption included in hemisphere") and the 19 × 19 THETAOBSAC/PHIOBSAC grid with 31 bands from 10 Hz to 10 kHz.
  - TAMB is 298.1 K and RELHUM 70% in every forward-flight file.

  — [Package Hemispheres/](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/>) (parsed by script).
- **Fill fraction.** Between 37.7% and 67.0% of the 361 (φ, θ) bins carry data (mean 54.5%). Examples:
  - R22 approach at 53 kt and −6°: data at φ = −60…+60 (with a few cells at −80, −70 and +70) and θ = 20…160.
  - EC135 flyover at 102 kt: φ = −60…+60 and θ = 10–30…150–160.

  — [Package Hemispheres/R22_Approach_53kts_6deg.hem, EC135_Flyover_102kts_0deg.hem](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/>).
- **The one hover-format file**, `R22_Ingroundhover_0kts_0deg_test.hem`:
  - Its internal title is `B412_HIGE_SoftGrnd`.
  - It has POLDIST 70, and FREEFIELD 0, meaning "GR & Atmospheric absorption included in hemisphere".
  - It has a single axis, PHIOBSEC, with 13 values from −180 to 180° in 30° steps, by 31 bands. There is no elevation dependence.

  — [Package Hemispheres/R22_Ingroundhover_0kts_0deg_test.hem](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/R22_Ingroundhover_0kts_0deg_test.hem>).
- At load, SingleEvt logs "Current Polar distance[m]: 70 / New Polar distance[m]: 60" for that file. It rescales to the common 60 m reference — [Package ARP/Case4/log.txt](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case4/log.txt>).
- `R22_Fullrpmidle_0kts_0deg_test.hem` is a full 19 × 19 hemisphere whose internal title and metadata are those of `R22_Flyover_55kts_0deg` (ACSPEED 55.3). It is a placeholder copy of a flyover hemisphere — [Package Hemispheres/](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/>).
- **Hover, idle and taxi directivity** is "defined as function of polar angle only". Here the hover "polar angle" is the horizontal angle, PHIOBSEC, "Polar angle used in hover hemispheres" — [Final report §3.3.2 p.22](https://www.easa.europa.eu/en/downloads/139729/en); [Manual §5.4.2 p.35](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- NORAH2 added "a specific hemisphere format ... for special operations such as hover-in-ground-effect, hover-out-of-ground-effect, taxi and idle" — [Final report §5.1 p.35](https://www.easa.europa.eu/en/downloads/139729/en).
- **Output band range.** The manual's feature list gives 1/3-octave output from 10 Hz to 10 kHz. IFR1 and IFRN select the band range (defaults 10 and 40). That choice also changes SEL(A), and PNLT needs at least bands 17–40 (50 Hz–10 kHz) — [Manual §2.2 p.10; Annex B p.48-49](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- The lowest band was lowered from 50 Hz to 10 Hz to capture the main-rotor blade passage frequency (BPF) — [Final report §2.4 p.16](https://www.easa.europa.eu/en/downloads/139729/en).

### Inferences
- **Hemisphere versus sphere.** The file format and grid are a genuine lower hemisphere. Upward emission (φ beyond ±90°) does not exist, and populated cells end around ±60–70° lateral. That is enough for ground observers from an aircraft at altitude, but grazing paths fall back on the nearest-bin fill (§5 below).
  - NICE-OPS uses a full sphere (per its CLAUDE.md, phi spans the upper surface after the legacy-hover repair). Its sphere radius is set by the database, not fixed at 60 m.
  - The channel angle convention looks the same in substance. NORAH φ = atan2(lateral, down) with 0 = down, and θ = angle from the nose. NICE-OPS's phi = atan(across/|up|) with 0 = down, and theta from the nose (per CLAUDE.md).
- **Content per cell.** NORAH stores full band spectra in every cell. NICE-OPS's default path stores broadband dBA plus an EAA factor; its `--spectral` path carries bands.
- **The RmOmega field is inconsistent across files.**
  - Most values are about 100 (that is, % RPM; D1.5d lists the unit as "%" while the files say "rpm").
  - Several are negative (−94 to −926) or −999. Files affected: A109 at 60 kt/12°, 73 kt/6°, 113 kt and 126 kt; EC120 at 53 kt/11° and 67 kt; EC135 at 57 kt/5°; R66 at 80 kt/5.5°, 105 kt and 57 kt; four takeoff files at −999.
  - RPM is metadata only and does not enter the calculation, so this is a data-quality signal rather than a model input.

### Gaps
- What the FREEFIELD codes mean beyond 0 and 2 (for example, whether 1 is free field without absorption).
- How a FREEFIELD = 0 hover table, which includes the ground reflection, is combined with NORAH's own ground model. Is GR added again or skipped? This could not be determined.
- How the 70 m to 60 m rescale is done: presumably 20·log10(70/60) plus absorption.

## 3. Parameterisation by operating state, and selection or interpolation between hemispheres

### Takeaway
- **Parameters.** The operating state is just two numbers: **airspeed V** and **flight-path angle γ**. No thrust or CT, advance ratio, acceleration, weight or RPM.
- **Interpolation.** NORAH2 normalises them, builds a 2-D Delaunay triangulation (precomputed per type in `*_Triangulation.int`), and blends the three corner hemispheres. The weights are **inverse normalised distance, applied to energies**, not barycentric weights. Outside the convex hull it takes the nearest hemisphere.
- **Special modes.** Hover, idle and taxi are discrete modes with fixed dB offsets. Turns reuse level-flight (or climb) hemispheres tilted by bank angle.
- **No effective flight path angle.** There is no QSAM-style effective γ: γ is the geometric path angle.

### Cited Findings
- The flight condition "is characterised by both airspeed and rate of flight angle". BVI makes descent the most sensitive region. In the NORAH1 R22 data, LAmax varies by up to 10 dB over the envelope but generally within 1 dB with rate of climb — [D1.5d §A.3.1 p.11](https://www.easa.europa.eu/en/downloads/132005/en).
- **Recommended database density:**
  - descent angles at most 3° apart, at 4 velocities;
  - climbs at for example 3, 6 and 9° at Vy, plus the maximum climb angle;
  - level flight at 0.9 VH, plus +10 kt (or VH), −15 kt and −30 kt.

  — [D1.5d §A.3.1 p.11](https://www.easa.europa.eu/en/downloads/132005/en).
- **Interpolation algorithm** (D1.5d eqs. 3–10):
  1. Normalise γ and V by the database range (max − min). Multiply γ by a flight-condition scaling factor F_fc = 2, "empirically determined to give minimum interpolation error".
  2. Build the Delaunay triangulation of the normalised (γ̄_j, V̄_j). This can be a lookup table.
  3. If the target lies inside a triangle T_k, compute δ_j = sqrt((γ̄ − γ̄_j)² + (V̄ − V̄_j)²) and blend: L̂ = 10·log10( Σ_{j∈T_k} 10^(L̃_j/10)/δ_j  /  Σ_{j∈T_k} 1/δ_j ).
  4. Otherwise take the nearest neighbour, argmin δ_j.

  — [D1.5d §A.3.1 p.12-13](https://www.easa.europa.eu/en/downloads/132005/en).
- **Triangulation file.** `Heli_triangulation.int` is part of the database and "should not be edited". It holds three namelist-headed tables:
  - the hemisphere list (iHem, HEMSpeed in kt, HEMAngle in deg, file name);
  - the Delaunay triangles (iHem1–3);
  - default corrections for out-of-ground hover, reduced-RPM idle and full-RPM idle.

  — [Manual §5.3.4 p.33-34](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Shipped triangulations.**
  - R22: 20 hemispheres in 27 triangles, including a 0 kt / 0° in-ground-hover entry.
  - The others: A109 21/35, AS350 17/26, B412 15/23, Cabri 17/28, EC120 28/47, EC135 19/30, R44 17/26, R66 20/31, S300 13/19, S92 16/24 (hemispheres/triangles). None of these has a hover entry.
  - Every file carries the same correction table: Outoffgroundhover +8 dB, Reducedrpmidle −10 dB, Fullrpmidle −2 dB.
  - Speeds and angles in the `.int` are the measured averages (e.g. 49.9 kt / −7°), not the nominal file-name values.

  — [Package Hemispheres/*.int](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/Hemispheres/>).
- **Hover, idle and taxi** (D1.5d §A.3.5). Measured data are preferred. Otherwise there are three fallback approaches. The third derives the other modes from in-ground hover with constant directivity:
  - HOGE = HIGE + 12 dB;
  - reduced-RPM ("ground") idle = HIGE − 12 dB;
  - full-RPM ("flight") idle = HIGE − 2.5 dB.

  These offsets came from inverted ground-plate microphones, and EMPA found they may not hold for other setups. With no hover data, a "hover hemisphere" may be built from 2 sideline microphones in the slowest level flight, corrected to a 150 m circle, with constant directivity in φ. Taxi uses in-ground hover (wheeled or skid) and full-RPM idle hemispheres — [D1.5d §A.3.5 p.17-19](https://www.easa.europa.eu/en/downloads/132005/en).
- **Modes in NORAH2:** START, DEP, FIDL (full-RPM idle), RIDL (reduced-RPM idle), TAX, IGH, OGH, FOV and APP. Each trajectory point carries an iMode — [Manual §5.1.6 p.26; §5.3.1 p.30](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- Profile steps give speed and vertical angle; "the hemisphere will be interpolated from other hemispheres based on the given values" — [Manual §5.1.6 p.26](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Turns** use level-flight hemispheres tilted by bank angle, with the closest-matching speed:
  - in level turns, the level hemisphere;
  - in take-off turns, the take-off hemisphere;
  - in descending turns, the **level-flight** hemisphere.

  — [D1.5d §A.3.4 p.17](https://www.easa.europa.eu/en/downloads/132005/en).
  - The Final report instead says "closest matching velocity and flight trajectory angle", with interpolation selecting automatically — [Final report §3.3.1 p.22](https://www.easa.europa.eu/en/downloads/139729/en).
- **Speed for radar tracks.** The 3-D speed V_A = sqrt(ΔX² + ΔY² + ΔZ²)/Δt (not the ground speed V_g) "should be used to find the best matching hemisphere". The resulting velocity is "with respect to ground surface, which neglects wind" — [D1.5d §A.3.3 p.16-17](https://www.easa.europa.eu/en/downloads/132005/en).
- **Speed in the files.** Hemisphere metadata store the *indicated* airspeed (ACSPEED, kt) — [D1.5d App. A p.34-35](https://www.easa.europa.eu/en/downloads/132005/en). Radar input asks for the "indicator speed" — [Manual §5.2.2 p.27](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Substitution.** A type without its own data uses a class reference's hemispheres plus ΔL_EPNL, the difference in certification levels: take-off for climb, overflight for level, approach for descent. Mirrored main/tail-rotor configurations (e.g. [A600] in the R22 class) reverse the azimuth (φ → −φ) — [D1.5d §A.2.2 p.9](https://www.easa.europa.eu/en/downloads/132005/en); [Manual Annex C p.52](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Substitute.csv** carries per-(Heli, HeliHem, Mode, speed range, angle range) ddB values; the last matching row wins. A helicopter may substitute only one HeliHem type — [Manual §5.1.4 p.22-23](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Out of scope in NORAH1:** accelerated and decelerated flight, transitions between steady conditions, turns and hover. Further research was recommended — [Final report §2.4 p.16](https://www.easa.europa.eu/en/downloads/139729/en).
- **Measurement basis.** Flight tests ran at ≥ 90% MTOW — [Final report §2.2.2 p.13](https://www.easa.europa.eu/en/downloads/139729/en). Rotor speed was held within ±1% and airspeed within ±5 kt — [Webinar slide 37](https://www.easa.europa.eu/en/downloads/139577/en).

### Inferences
- **The implementation matches the D1.5d normalisation.** In the example `.his` output (`ARP/Case3/SE/A600_H1_APP_STD2_NE.his`):
  - a column `TriOrNear` reports positive values (triangle IDs, e.g. 9 = R22 hemispheres 1, 3 and 17) and negative values (−4, −2, −20);
  - the negative values are evidently "nearest hemisphere #"; the manual says only "Triangle ... or nearest hemisphere", p.35;
  - for a 50 kt / −5° step outside the R22 hull, the code picked hemisphere 4 (53.3 kt / −5°) over hemisphere 1 (49.9 kt / −7°);
  - that choice is what the γ-doubling (F_fc = 2) normalisation predicts; raw (kt, deg) distance would pick hemisphere 1.
- **Comparison with NICE-OPS.**
  - NORAH: 2-D (V, γ), inverse-distance energy blend of 3 corners, nearest-neighbour extrapolation.
  - NICE-OPS: 3-D Delaunay over (μ, γ, CT) with barycentric weights over 4 corners, per its CLAUDE.md.
  - NORAH's inverse-distance blend within a triangle is not a partition of unity in the barycentric sense, so it is not linear-exact. It is discontinuous in gradient across edges, and singular (correctly limiting to the vertex) as δ → 0.
  - NORAH has no CT or load-factor parameter. Thrust effects (weight, load factor in turns) enter only through the bank tilt and the measured condition.
- **The shipped hover/idle offsets (+8, −10, −2 dB) differ from the D1.5d Approach 3 values (+12, −12, −2.5 dB).** In the public package, R22 hover uses a B412 soft-ground HIGE table as a "test" stand-in. No other type has a hover hemisphere, so hover or idle for other types would fall back on the nearest forward-flight hemisphere.
- **No QSAM-like treatment.** No source maps acceleration or deceleration into an effective path angle. Hemispheres are steady-state only, and the measured pitch is "included implicitly" (§4).

### Gaps
- Whether the 1/δ weights use the normalised distances exactly as written (eq. 8). Also whether hemisphere bins that are −999 in one corner are gap-filled before blending: presumably yes, per eq. 14, but this is not stated.
- How OGH, idle and taxi pick between the hover hemisphere, the triangulation and the correction table in software.

## 4. Hemisphere orientation in the flight frame

### Takeaway
- **Body axes as specified.** D1.5d defines the hemisphere in "aircraft body axes" but says pitch attitude is implicit in the hemisphere and not to be applied.
- **What the code actually does.** SingleEvt orients the hemisphere by a **horizontal rotation angle** (track heading plus slip) and **tilts it by the roll/bank angle**. It applies **no pitch and no flight-path-angle tilt**; this was verified against the shipped outputs to about 0.003°.
- **Possible roll-sign problem.** In the shipped departure examples, the roll sign appears to bank the hemisphere toward the outside of turns.

### Cited Findings
- "Rotorcraft pitch attitude is included implicitly in the hemispheres and shall not specified explicitly." The bank angle "allows to accurately model the emission pattern" — [D1.5d §A.3.3 p.16](https://www.easa.europa.eu/en/downloads/132005/en).
- **Bank from a radar track:**
  - Φ = atan(K V_g²/g), with curvature K = ΔΘ/ΔS;
  - heading Θ = atan2(ΔX, ΔY);
  - path angle γ = acos(ΔZ/ΔS), as printed. As written this is dimensionally odd; asin or atan was presumably meant.

  — [D1.5d §A.3.3 eqs. 18-21 p.17](https://www.easa.europa.eu/en/downloads/132005/en).
- **Angle definitions:**
  - horizontal angle: clockwise from the Y_M (north) axis to the track;
  - slip: from track to rotorcraft orientation, positive clockwise;
  - rotation angle = horizontal angle + slip, "used to orient the hemisphere within the microphone coordinate system";
  - backward departure means slip = 180°.

  — [Manual §4 p.16-17](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- `USE_ROLL` (default 1) makes the multi-event module estimate roll as atan(V²/gR) for track arcs. It is not available for radar operations, which must supply roll — [Manual Annex B p.49; §5.1.5 p.24; §5.2.2 p.27](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- The single-event trajectory row is TEM, X, Y, Z, iMode, SPEED, VANG, HEADING, ROLL, ddB; **there is no pitch column** — [Manual §5.3.1 p.30-31](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>); [Package ARP/Case1/SE/A600_H1_APP_STD2_NE.inp](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case1/SE/A600_H1_APP_STD2_NE.inp>).
- **Reverse-engineering check (own analysis).** The `.inp` trajectory rows were joined with the `.his` outputs (THETAOBS, PHIOBSAC), using the observer at microphone height (ZMIC + HMIC). θ = acos(d·f) and φ = atan2(d·s, d·down) were recomputed under several frame hypotheses. The fit was evaluated over 333–1621 samples per file in `ARP/Case2`/`Case3` (`R22_H1_DEP_STD1_NE`, `R22_H1_APP_STD1_NE`, `A600_H1_APP_STD2_NE`, `A600_APP_RT1`):
  - **Heading plus roll:** rms error 0.003° in both θ and φ.
  - Adding a pitch tilt by VANG: errors up to 0.9–4.4° in θ and 3–54° in φ.
  - Heading alone: φ off by exactly the roll (3.12°) in turns.

  — [Package ARP/Case2, Case3 SE files](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/>).
- **Sign checks from the same data.**
  - On a straight leg with slip 90°, PHIOBSEC = bearing-to-observer − heading = −91.6° for an observer to port.
  - PHIOBSAC = −83.8°, which equals −(90° − 6.2° depression).
  - Both confirm negative φ = port, consistent with D1.5d.
- **Wind** is recorded in hemisphere metadata (TW/CW/HW) but is not an input to NORAH2 calculations; speed is "with respect to ground ... neglects wind" — [D1.5d §A.3.3 p.16](https://www.easa.europa.eu/en/downloads/132005/en).

### Inferences
- The effective NORAH frame is x along the rotation angle (heading + slip) in the **horizontal plane**, rolled by ROLL about that axis. Body pitch and γ are not applied, because the hemisphere was measured in that steady condition.
  - Compare NICE-OPS: a `SphereFrame` from the velocity and the TPP normal (nose, normal, starboard), with an optional "qsam" effective γ (per its CLAUDE.md).
  - The NORAH approach is equivalent to assuming the measured hemisphere's attitude recurs whenever (V, γ) recurs. NICE-OPS instead re-derives orientation from the trajectory's TPP normal.
- **Roll sign, possible bug (needs confirmation before relying on it).** In the shipped `R22_H1_DEP_STD1_NE` departure:
  - ROLL = −2.34° on the right turn (heading 45° to 90°) and +3.12° on the left turn (heading 90° to 0°);
  - the fitted rotation applies positive ROLL as starboard-wing-down;
  - so the departure hemisphere is banked toward the **outside** of each turn.

  In the mirrored approach (`R22_H1_APP_STD1_NE`), the same segments are flown in the opposite sense and keep the same ROLL signs, which there bank the hemisphere *inward*, as is physical.

  This fits the manual's caveat that track L/R orientation is defined "as seen in a departure" (§5.1.5 p.24). The roll sign appears to be computed from that departure-sense L/R and not reversed for approaches, or vice versa.
  - The analysis rests on my reading of the code's φ sign, which was confirmed independently on a no-roll straight leg.
  - The approach file still shows a 0.15° rms (1.1° max) φ residual under the best hypothesis, unexplained.

### Gaps
- Whether the hemispheres were rotated to compensate for the measured pitch and roll during processing (the metadata record PITCH and ROLL). Nothing says they were.
- Whether hover hemispheres follow the heading with slip, as is implied by PHIOBSEC = bearing − heading, which matches the data.

## 5. Angular interpolation within a hemisphere

### Takeaway
NORAH uses **bilinear interpolation in energy** (10^(L/10)) over the 10° (φ, θ) cells. Empty bins are filled by **constant-value extrapolation from the nearest filled bin**, nearest by great-circle angle, with ties averaged energetically. Interpolation between flight conditions is also done in energy.

### Cited Findings
- Eq. 13 is the bilinear form of the four corner energies, 10^(L_h(f_c)_{m,n}/10) etc., weighted by (φ_{m+1} − φ)/Δφ, (φ − φ_m)/Δφ and the same for θ, then 10·log10 — [D1.5d §A.3.2 p.15](https://www.easa.europa.eu/en/downloads/132005/en).
- **Outside the populated range:** "constant value extrapolation should be applied from the nearest filled data bin, filling any gaps". The nearest bin minimises ρ_{m,n} = acos(x(φ, θ)·x_{m,n}), using unit vectors from eq. 11. "In case there are multiple closest bins, the energetic average of the closest bins is taken" — [D1.5d §A.3.2 eqs. 14-15 p.15-16](https://www.easa.europa.eu/en/downloads/132005/en).

### Inferences
- Because roughly 45% of every hemisphere is −999, grazing and near-axis emission angles (φ beyond ±60°, θ near 0 or 180°) get **constant extrapolation** of the outermost measured bin, with no physics of the unmeasured directivity.
  - This matters for distant lateral observers, where a 300 m altitude source at 3 km is at φ ≈ 84°.
  - NICE-OPS has a full sphere from prediction, so no such fill is needed.
- **NORAH interpolates in energy; NICE-OPS blends dBA in dB, linearly** (per its CLAUDE.md `Pair` and `cell_stencil` descriptions). Energy averaging biases toward the louder corner. The two will differ most at steep directivity gradients and between dissimilar conditions.

### Gaps
- Whether the gap-fill is applied per band or per bin across all bands. The text implies per bin.
- How the θ = 0/180° poles (where φ is degenerate) are handled.

## 6. Doppler, amplitude and emission time

### Takeaway
NORAH evaluates each emission point at its emission time and moves the result to the **reception time t_r = t_e + r/c**, with c fixed at 346.1 m/s. Time-integrated metrics are integrated over reception time. **No Doppler frequency shift or convective amplitude correction is mentioned** anywhere, and no de-Dopplerization of measured data is described. Hemispheres are built from ground microphones in steady flight, so they presumably retain the in-flight Doppler shift.

### Cited Findings
- "t_r = t_e + r/c ... At ICAO reference conditions ... c = 346.1 m/s. To evaluate time integrated metrics, e.g. SEL or EPNL, it is required to express the predicted noise levels as a function of recorded time" — [D1.5d §A.4.1 eq. 22 p.20](https://www.easa.europa.eu/en/downloads/132005/en).
- The `.his` output lists, per sample: aircraft XAC/YAC/ZAC, `tem` (emission time), `trec` (reception time), DOBS (distance), THETAOBS, PHIOBSAC, PHIOBSEC, LA and PNLT — [Manual §5.4.2 p.35](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- In a shipped `.his`, DOBS = 3903.27 m with tem 703.50 s and trec 714.78 s, giving an implied c of 346.0 m/s. The sound speed is constant, apparently independent of TAMB — [Package ARP/Case3/SE/A600_H1_APP_STD2_NE.his](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case3/SE/A600_H1_APP_STD2_NE.his>).
- **Hemisphere processing steps:**
  - corrections for cable, microphone response, free field and windscreen;
  - last-good-band (LGB) handling, scaling to 60 m, and an equal-energy check above 8 kHz;
  - 10° binning, then merging of repeat runs.

  No Doppler step is listed — [Final report §2.2.3 p.14-15](https://www.easa.europa.eu/en/downloads/139729/en).

### Inferences
- The emission angles are evidently the emission-time (retarded) angles. The stored hemispheres therefore represent the field as heard from those angles in that steady flight, Doppler-shifted and convectively amplified, which is the same convention NICE-OPS uses ("spheres keep Doppler").
- Re-applying a Doppler shift in NORAH would double-count it; neither model does.
- Substitution between speeds (interpolation between hemispheres at different V) does not correct for the Doppler difference.

### Gaps
- No public statement on whether NORAH holds the Doppler convention deliberately.
- How the reception-time series is resampled before integration, given that reception times are non-uniform.

## 7. Propagation: spreading, absorption, atmosphere, ground, terrain, refraction, shielding

### Takeaway
The level at the observer is ΔL_p = ΔL_s + ΔL_a + ΔL_g + ΔL_d, per band:
- **ΔL_s, spherical spreading:** −20·log10(r/60).
- **ΔL_a, absorption:** SAE ARP5534 pure-tone coefficients converted to bands by Rickley's method, applied as −α(r − 60).
- **ΔL_g, ground:** a Chien–Soroka/Delany–Bazley spherical-wave ground reflection with a band-averaging coherence term.
- **ΔL_d, diffraction:** CNOSSOS-EU-based diffraction and mean-ground-plane ground effect for DEM terrain.

The atmosphere is uniform (no gradients). Refraction, wind and turbulence are **not** modelled, and there is no empirical lateral attenuation.

### Cited Findings
- **Decomposition.** ΔL_p = ΔL_s + ΔL_a + ΔL_g + ΔL_d, with ΔL_s = −20·log10(r/r_h) and r_h = 60 m — [D1.5d §A.4.1 eqs. 23-24 p.20](https://www.easa.europa.eu/en/downloads/132005/en).
- **The manual's `.sub` file** shows the per-band breakdown SPLmic = SPLrRef + ATMABS + GR + 20logDobs:
  - SPLrRef is the free-field level at 60 m, with the hemisphere's built-in absorption to 60 m "eliminated" at the hemisphere's ambient conditions;
  - ATMABS is the absorption from rotorcraft to observer by SAE ARP 5534 with Rickley's band correction;
  - the GR row is DIFFR/GR.
  - In the example, 20log(Dobs/60) = −17.5 at DOBS 449.6 m, and GR is about +5–6 dB at low frequency for a 1.2 m microphone.

  — [Manual §5.4.2 p.36](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Absorption.**
  - D1.5d eq. 27 gives the pure-tone ARP5534 coefficient with f_rO = 630.7 Hz and f_rN = 75692 Hz (the text swaps the O/N labels), at ICAO reference conditions.
  - Table 4 tabulates the band dB/km values: 0.0 below 50 Hz, 6.3 at 1 kHz, 101.0 at 10 kHz.
  - ΔL_a = −α(f)(r − r_h).

  — [D1.5d §A.4.2 p.20-22](https://www.easa.europa.eu/en/downloads/132005/en).
- **Atmosphere inputs:** TAMB (default 298 K, "No vertical gradient"), RHAMB (70%) and PAMB (101325 Pa) are scalar run options — [Manual Annex A p.46; Annex B p.48](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Refraction** "is not covered by the present method" because it needs temperature and velocity profiles — [D1.5d §A.2.3 p.10](https://www.easa.europa.eu/en/downloads/132005/en). Atmospheric refraction is listed as a possible future module — [Manual §2.3 p.11](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Ground model.**
  - The field is p = e^{ikr1}/r1 + Q e^{ikr2}/r2.
  - Treating Q as constant within a band: ΔL_g = 10·log{1 + (r1/r2)²|Q|² + 2(r1/r2)|Q|·I}.
  - I = [sin(0.727 f_c ΔR/c)/(0.727 f_c ΔR/c)]·cos(6.325 f_c ΔR/c + ψ). As frequency rises, I → 0 (incoherent sum).
  - Q = R_p + (1 − R_p)F(d), with R_p = (Z_s cos ξ − 1)/(Z_s cos ξ + 1).
  - Boundary loss F(d) = 1 + i d√π e^{−d²} erfc(−i d), with numerical distance d = ((1 + i)/2) √(k r2) (1/Z_s + cos ξ).
  - The erfc is evaluated by Abramowitz–Stegun and Matta–Reichel (App. B).

  — [D1.5d §A.4.3 eqs. 28-34 p.22-23; App. B p.36-37](https://www.easa.europa.eu/en/downloads/132005/en).
- **Impedance:**
  - Delany–Bazley, Z_s = 1 + 0.0511(f/σ)^−0.754 + i·0.0768(f/σ)^−0.732, locally reacting;
  - recommended σ is concrete 65·10⁶ Pa·s/m² for cities and grass 200·10³ for rural areas;
  - Table 5 gives the CNOSSOS-EU classes A–H, from 12.5 to 200 000 kPa·s/m².

  — [D1.5d §A.4.3 eq. 35, Table 5 p.23-24](https://www.easa.europa.eu/en/downloads/132005/en).
  - The manual's Annex D instead lists Encyclopedia of Acoustics values in kPa·s/m²: concrete, painted 200,000; concrete 65,000; grass 200; snow 29. The default SIGMA = 200000 Pa·s/m² is grass — [Manual Annex D p.54; Annex B p.48](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Varying terrain:**
  - each source–receiver vertical section is represented by a least-squares mean ground plane (MGP), as in CNOSSOS-EU;
  - the ground term uses equivalent heights and the incidence angle relative to the MGP;
  - source and receiver heights are floored at 0.1 m above the MGP;
  - mixed ground uses a log-averaged σ weighted by segment length.

  — [D1.5d §A.4.4 eqs. 36-41 p.24-27](https://www.easa.europa.eu/en/downloads/132005/en).
- **Screening and diffraction (CNOSSOS-EU):**
  - path-length difference δ; ΔL_d = 10·C_h·log10(3 + 40 C″ δ/λ) when 40 C″ δ/λ ≥ −2, with C_h = min(f_m h_0/250, 1);
  - C″ for multiple diffraction;
  - bounds 0–25 dB;
  - source-side and receiver-side ground terms come from image-source and image-receiver diffractions (eqs. 45–47);
  - the "rubber-band" convex path over DEM points, and only terrain points above the line of sight count as diffraction points.

  — [D1.5d §A.4.5 p.27-33; App. D p.42-43](https://www.easa.europa.eu/en/downloads/132005/en).
- **Implementation notes from the Final report:**
  - topography screening is applied only if the line of sight is obstructed, giving at least 5 dB;
  - negative path differences are ignored for terrain;
  - screening by thin constructed obstacles needs explicit obstacle input, which "is currently not supported in NORAH2".

  — [Final report §3.2 p.20-21](https://www.easa.europa.eu/en/downloads/139729/en).
- **Topography and sigma inputs** are ESRI ARC/INFO ASCII grids. Without topography the ground is "flat, soft ground" — [Manual §3.1 p.13; §5.3.2-5.3.3 p.32-33](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- Example Case3 says that without terrain files "Diffraction effects will not be calculated, only absorption and Ground reflections" — [Package ARP/Case3/Case3_Run_Options.txt](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case3/Case3_Run_Options.txt>).
- **Ground-plate microphones.** Functional requirement 15 is "a correction for ground plate microphones". Per-observer HMIC can be as small as 0.01 m in the example irregular grid — [Manual §2.1 p.10](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>); [Package ARP/grids/irregularGridXYZHmicSigma.txt](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/grids/irregularGridXYZHmicSigma.txt>).
- **Urban use** would need vector building input and ray-path tracking, which are not in the prototype — [Final report §6.2 p.38-39](https://www.easa.europa.eu/en/downloads/139729/en).

### Inferences
- **Close match to NICE-OPS's spectral ground path.** NORAH's ground model is the same physics as NICE-OPS's `--spectral --ground_effect` (the panam `ega()` port): a Chien–Soroka spherical reflection, Delany–Bazley impedance, a Faddeeva/erfc boundary loss, and a sinc band-averaging term. A point-by-point comparison should check the coherence constants (0.727 and 6.325 f_c ΔR/c) and the boundary-loss evaluation (NORAH uses the A&S and Matta–Reichel approximations at about 1e-6; NICE-OPS uses Weideman's).
- **Absorption.** NORAH's absorption is ARP5534 with Rickley band conversion. NICE-OPS uses ISO 9613-1 per band (spectral path), or a fitted EAA (broadband).
- **What NORAH lacks** relative to NICE-OPS's newer propagation work: refraction, layered atmospheres, wind, turbulence, shadow zones and PE. NORAH has terrain diffraction and MGP ground, which NICE-OPS does not (NICE-OPS assumes flat ground, per its CLAUDE.md).
- **Lateral attenuation** is handled by the explicit ground model and directivity, not an SAE-AIR-5662-style empirical lateral attenuation. None is mentioned.

### Gaps
- Whether ATMABS is evaluated at the user's TAMB/RHAMB/PAMB or fixed at ICAO reference conditions. D1.5d says reference conditions; the manual exposes ambient inputs.
- Whether NORAH's I term includes turbulence decoherence (none is stated).
- How a FREEFIELD = 0 hover table interacts with GR.

## 8. Receivers, metrics, time integration and outputs

### Takeaway
- **Receivers:** points on a regular or irregular grid, each with its own height above ground (default 1.2 m), flow resistivity and terrain elevation.
- **Metrics:** LAmax, SELA, LAE, PNLTM, EPNL (with 10-dB-down flags), NAx, and user-defined weighted Leq metrics such as Lden.
- **SEL window: the full time history by default**, not 10-dB-down. The 10-dB-down window is a certification option.
- **Time step:** that of the trajectory (default DTIME 0.5 s).

### Cited Findings
- **ISEL10DB** sets the SELA integration time: 0 = full time history (land-use planning, the default); 1 = the 10-dB-down period (certification) — [Manual Annex A p.46; Annex B p.48](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **I10db flag** in the `.onl`: 3 means no 10-dB-down in either PNLT or LA; 2 = none in LA only; 1 = none in PNLT only; 0 = both found — [Manual §5.4.1 p.35](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Per-observer metric columns** in `.onl`/`.his`: LAE, LAMAX, PNLTM, EPNL, SELA, T1pnlt and T2pnlt (the PNLT 10-dB-down times) — [Manual §5.4.1-5.4.2 p.35](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **DTIME** (default 0.5 s) is used only by the multi-event module to discretise the flight path. "SingleEvent uses all the points in the input file (it does not interpolate)".
  - Guidance: 0.5 s for low or fast flights; 1 s for high or slow ones.
  - Too large a step gives "corrugated" contours; halving DTIME roughly doubles run time.

  — [Manual Annex B p.50](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Culling options:**
  - DMAX (default 3000 m): grid points farther than DMAX from every flight-path point are not computed and are output as −999.
  - GRDEXT (default off; at least 2000 m recommended): deletes flight-path points farther than GRDEXT from the grid.

  — [Manual Annex B p.50-51](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Cumulative metrics:** Leq_w = 10·log10((t0/T0) Σ w_i 10^{0.1 L_sei}), with L_se being SELA. Metrics.csv defines SELA, LAeq24h, Lden (evening weight 3.162, night 10), Ldn, Ln and Ld; only A-weighted metrics can be defined. LAmax is built in. NAx thresholds are set in options (e.g. NA65 NA70) — [Manual §5.1.3 p.21-22](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Requirement 14:** output in SEL, EPNdB, LAmax and NAX, several metrics per run — [Manual §2.1 p.10](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- NAx and EPNL are new in NORAH2 — [Final report §5.1 p.34](https://www.easa.europa.eu/en/downloads/139729/en).
- **Observers:**
  - The OBSFILE has XMIC_M, YMIC_M, optionally ZMIC_M, HMIC and SIGMA. Alternatively a regular grid is set by X0/Y0/LX/LY/DX/DY.
  - Defaults are HMIC 1.2 m, SIGMA 200000 and Z0_M 0.
  - Precedence for sigma: sigma file, then observer file, then default. The sigma used is echoed in output.
  - An irregular grid gives no NMPlot map.

  — [Manual §5.1.2 p.19-21; Annex A p.47](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Output files:**
  - single-event: `.onl` (overall levels); optional `.his` (time history), `.spc` (spectra versus time) and `.sub` (per-band breakdown), selectable per observer or all;
  - multi-event: `Case_Metric.onl` (matrix), `.irr` (irregular grids), `Case.ONL` (all metrics) and `Case_Metric.grd` (NMPlot, feet by default, UNITSNMPLOT=METR for metres);
  - ZEROVALUE floors the multi-event output at 0 dB(A), or −999 for NORAH1 behaviour.

  — [Manual §5.4-5.5 p.34-38; Annex B p.49](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- Single-event noise is the 1/3-octave level versus reception time, from which "any required time integrated metric may be calculated". Multi-event and contours follow ECAC Doc 29 Vol. 2 chapters 5–6 — [D1.5d §A.5 p.33](https://www.easa.europa.eu/en/downloads/132005/en).

### Inferences
- **LAE appears to be the event-average LAeq.**
  - In the shipped `.his`, SELA − LAE = 29.9 dB, which implies T ≈ 980 s. That is about the reception-time span of the 1031 s trajectory.
  - In the single-point hover case (`ARP/Case4/SE/R22_DEP_RT2`), LAE = LAMAX = SELA = 103.64 dB.
  - So LAE ≈ LAeq over the event and SELA = LAE + 10·log10(T/1 s). The manual does not define LAE.
- **Contrast with NICE-OPS.** NICE-OPS's SEL uses the 10-dB-down window walked out from the peak (its "first crossing" rule is called out as a fragility). NORAH's land-use default integrates the **whole** history, so it has no window knife-edge. It does, however, depend on trajectory length and on the DMAX/GRDEXT truncation.
  - When comparing footprints with NORAH, set ISEL10DB = 1 or account for the difference.
- **DTIME = 0.5 s** in NORAH versus the 0.02 s shipped tracks in NICE-OPS. NORAH's guidance treats the step as an accuracy/cost trade, consistent with NICE-OPS's finding that SEL is heavily oversampled at 50 Hz.

### Gaps
- Exact EPNL implementation details: tone correction, and whether integration is on half-second PNLT.
- Whether SELA integrates on reception time with uniform or variable weights.
- The documents contain no DNL or CNEL metric names beyond Ldn (the US DNL equivalent is Ldn with night weight 10).

## 9. Multiple sources and vehicles, flight-track input, noise-abatement features

### Takeaway
NORAH2 has three trajectory inputs:
- tracks & profiles, Doc 29 style: segments, procedure steps and modes;
- smoothed radar trajectories, in X/Y or lat/lon via an EPSG projection, carrying speed, vertical angle, roll, slip and mode;
- direct single-event `.inp` 4-D trajectories.

The multi-event module runs one single-event per operation and energy-sums them with day/evening/night weights. There is one source per vehicle (the hemisphere), no separate rotor sources, and no built-in procedure optimisation.

### Cited Findings
- The three trajectory input options are tracks&profiles, radar trajectories and NORAH trajectories (`.inp`). The `.inp` can be produced by third-party software — [Manual §1.1 p.9](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Tracks:**
  - segments of type O (initial point with heading), S (straight, length) or R/L (arc, radius and turned angle), plus a slip angle;
  - the last segment must be straight;
  - tracks serve both approach and departure, but L/R is defined "as seen in a departure".

  — [Manual §5.1.5 p.23-24](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Profiles:**
  - steps per Heli/Operation (APP or DEP)/Profile;
  - each step gives one of time, altitude or distance (the others computed), plus MODE, final speed (kt) and a constant vertical angle;
  - a START step is required;
  - overflights are made as APP or DEP with a virtual helipad.

  — [Manual §5.1.6 p.24-26; §8.1 p.43](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Radar trajectories:** RTID, t, lat/lon *or* X/Y, Z, speed (kt), Vang, roll, slip and mode. They "need to be smoothed before they are included in NORAH2 to avoid unrealistic jumps, especially in height" — [Manual §5.2.2 p.27](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>); [Package ARP/Radar_Trajectories.csv](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Radar_Trajectories.csv>).
- **Radar track quality:** 5–10 m noise and timestamps rounded to 1 s give 10–20 m error. Spline smoothing is recommended, resampled at for example 1 s, with V, γ and bank derived by finite differences — [D1.5d §A.3.3 p.16-17](https://www.easa.europa.eu/en/downloads/132005/en).
- **Coordinates:**
  - the topography frame (X, Y, Z above sea level) and the microphone frame (X_M, Y_M, Z_M) differ by translation only;
  - the origin is given in X/Y or in lat/lon with an EPSG/ESRI PROJ code, which must match the DEM's;
  - headings are clockwise from Y_M (north).
  - The example log shows EPSG:27394 converting the lat/lon origin to X/Y.

  — [Manual §4 p.16; §5.1.1 p.18-19](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>); [Package ARP/Case2/log.txt](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/Case2/log.txt>).
- **Operations:** Case_Operations.csv (Heli, StartPoint, OP, Profile, Track, D/E/N counts, fractional allowed) and Radar_Operations.csv (Heli, RTID, OP, D/E/N). Single events are merged into totals on the grid — [Manual §5.2 p.26-28; §3.1 p.13](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Flight paths** are built per ECAC Doc 29 Vol. 1 §4.8 as procedural steps from available hemisphere conditions — [D1.5d §A.3.3 p.16](https://www.easa.europa.eu/en/downloads/132005/en).

### Inferences
- **No noise-abatement design features.** Nothing documents optimisation or abatement-procedure tools. Abatement is assessed only by running alternative profiles.
- **The radar path uses the ground-referenced 3-D speed and no wind.** That contrasts with NICE-OPS's `--speed_reference air`.
- **No geodetic tracks for tracks & profiles.** The lat/lon route exists only for radar input, and conversion uses a projected CRS (pyproj-style EPSG codes) rather than a local transverse Mercator.

### Gaps
- Whether multiple helicopters in one event, or simultaneous sources, are supported. The design is one event, one vehicle.

## 10. Stated limitations, assumptions, validation and known issues

### Takeaway
NORAH2 is explicitly a prototype. Its limitations:
- steady-state hemispheres only;
- (V, γ) parameterisation only;
- no refraction, wind or turbulence;
- measured directivity only to ±60° lateral;
- class-based substitution not validated;
- hover, taxi and turn guidance not validated;
- urban building input not supported.

Validation so far covers data quality and implementation consistency, not single- or multi-event predictive accuracy:
- certification-level checks: maximum 1.6 dB, s.d. 0.4 dB;
- an S-92 certification re-simulation: maximum 1.8 dB, s.d. 1.2 dB;
- a Nord2000 comparison: typically ±1 dB;
- an independent MATLAB implementation.

### Cited Findings
- The prototype "comes with certain limitations", is "not optimized for efficiency", and has a simple CLI with ASCII I/O — [Manual §1 p.8](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **Fleet coverage:** about 350 types are represented by 8–15 classes, on the assumption that similar types have comparable directivity — [Manual Annex C p.52](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>). The extension raised coverage from about 70% to 80–81% of the European fleet (Cabri G2, A109, S-92 added) — [Final report §4 p.23; §6.1 p.38](https://www.easa.europa.eu/en/downloads/139729/en).
- **Validation against certification levels**, from tripod microphones during the new campaigns, gives maximum 1.6 dB and s.d. 0.4 dB:
  - Cabri: +0.1;
  - A109: 0.0, 1.0 and −1.6;
  - S-92: 0.7, 1.5 and 0.6.

  An S-92 certification run simulated from GPS tracks gives maximum 1.8 dB and s.d. 1.2 dB — [Final report §5.2.1 p.36](https://www.easa.europa.eu/en/downloads/139729/en).
- **Nord2000 comparison:** an R22 departure near a ridge with mixed ground. Agreement is generally ±1 dB, with ±10 dB deviations under strong screening and 2–3 dB for long ground-to-ground paths — [Final report §5.2.2 p.36-37](https://www.easa.europa.eu/en/downloads/139729/en); [Webinar slides 61-62](https://www.easa.europa.eu/en/downloads/139577/en).
- **Independent MATLAB re-implementation** of the screening and ground module from the guidance found ambiguities and "a couple of bugs", which were fixed — [Final report §5.2.3 p.37](https://www.easa.europa.eu/en/downloads/139729/en).
- **Recommended further validation:**
  - the substitution method;
  - taxi, hover and turn accuracy;
  - single-event comparison with measurement including topography and screening;
  - multi-event comparison against monitoring stations.

  The low-frequency "thumping" BPF content is not captured by existing metrics — [Final report §6.2 p.38-39](https://www.easa.europa.eu/en/downloads/139729/en).
- **Hover offsets** derived from inverted ground-plate microphones "may not be valid for other microphone setups (e.g., 4 m height over hard or soft ground)" (an EMPA finding) — [D1.5d §A.3.5 p.19](https://www.easa.europa.eu/en/downloads/132005/en).
- **Refraction** is not covered — [D1.5d §A.2.3 p.10](https://www.easa.europa.eu/en/downloads/132005/en).
- **DMAX trade-off:** small values (1000 m) can affect the high SELA contours; 3000 m is recommended — [Manual Annex B Fig. B.2 p.51](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).
- **CSV and tab handling:** stray tabs can make NORAH2 fail to read files, and space-separated CSV is not recognised — [Manual §5 p.18](<file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/doc/NORAH2%20User%20Manual%20(public).pdf>).

### Inferences
- **Known issues visible in the public package** (own inspection):
  - hover and idle hemispheres are "test" placeholders: a B412 HIGE table filed under R22, and a copy of the R22 55-kt flyover filed as R22 full-RPM idle;
  - only R22 has a hover entry in its triangulation;
  - the shipped hover/idle offsets (+8, −10, −2 dB) disagree with D1.5d (+12, −12, −2.5 dB);
  - RPM metadata is corrupt in about 15 files;
  - the roll sign looks inverted in departure turns (§4);
  - there are internal documentation inconsistencies: Substitutes.csv uses mode "HOV", which is not in the mode list, and sigma units are confused (Pa·s/m² versus kPa·s/m²; "Pa/(s.m2)").
- **None of the validation compares NORAH ground predictions with independent single-event flight measurements** at non-certification microphones. The certification check validates the measurement data. The S-92 check is the closest to a model validation, and it is at the certification microphones the hemispheres were partly built from.

### Gaps
- The NORAH 2.0 validation report (Olsen et al. 2024, SC03 D3) and the NORAH1 validation (Knepper & Olsen 2021, SC01 D2.5) were not obtained. They may hold single-event comparisons.

## 11. How hemispheres are generated

### Takeaway
Hemispheres are **empirical, from flight tests**. The setup follows ICAO Annex 16 Chapter 8 practice:
- a 420 m line of 17 inverted microphones on 40 cm ground plates, plus 3 microphones on 1.2 m tripods at the certification positions;
- DGPS and an on-board recorder;
- 0.1 s one-third-octave spectra from 10 Hz to 10 kHz.

Processing then:
1. corrects each spectrum to free field;
2. applies last-good-band masking repair;
3. scales to 60 m with ICAO reference absorption;
4. bins into 10° cells per run;
5. merges the repeat runs that meet speed and path tolerances.

Hover is measured on a quarter-circle of microphones at 30° spacing, 150 m from the hover spot. There is no prediction-based database and no stated de-Dopplerization.

### Cited Findings
- **NORAH1 campaigns:** R22, R44, R66, EC120, EC135, AS350, B412 and S300 at Flevoland, Bestkazerne and Marugán. 17 ground-plate microphones on a 420 m line, plus 3 tripods at 1.2 m. Coverage is polar angles within the SEL 10-dB-down interval and ±60° lateral. 1/3-octave spectra every 0.1 s from 10 Hz to 10 kHz, with a 10 m weather mast — [Final report §2.2.1 p.13](https://www.easa.europa.eu/en/downloads/139729/en).
- **On board:** a carry-on flight data recorder with DGPS, attitude, and video of the instrument panel processed live to extract IAS, RPM, VS and ALT. All tests at ≥ 90% MTOW. Over 170 test conditions, 800 runs and 60 flight hours — [Final report §2.2.2 p.13-14](https://www.easa.europa.eu/en/downloads/139729/en).
- **Processing:**
  - corrections for cable length, microphone frequency response, "free field" and windscreen;
  - last good band per ICAO Annex 16, with masked bands replaced by scaled LGB levels;
  - scaling to 60 m;
  - an equal-energy check above 8 kHz: if the next band rises by more than 3 dB, it is replaced by the current band;
  - 10° binning per run, then merging of repeat runs;
  - acceptance requires 90% of speed samples within ±5 kt and 90% of position samples within γ ± γ/12.

  — [Final report §2.2.3 p.14-15](https://www.easa.europa.eu/en/downloads/139729/en).
- **Masked-band repair method.** D1.5d says masked bands above 2 kHz are reconstructed by the ICAO Doc 9501 frequency-extrapolation method at 60 m, assuming a flat (equal-energy) spectrum after the last good band — [D1.5d §A.3.2 p.14](https://www.easa.europa.eu/en/downloads/132005/en).
- **NORAH2 campaigns:**
  - Cabri G2 at Mollerussa (Apr 2022): a 420 m line of 17 inverted ground-plate microphones plus 3 tripods, and a **quarter-circle hover array at 30° intervals** with both plate and tripod microphones.
  - A109 at Bestkazerne/De Peel (Sep 2023): hover spot 150 m from the array centre.
  - S-92 at Stavanger (Sep 2022): 12 microphones (9 on plates, 3 on tripods). Hemisphere construction "takes into account the actual microphone location".

  — [Final report §4.1.2 p.24; §4.2.2 p.27-28; §4.3.2 p.30-31](https://www.easa.europa.eu/en/downloads/139729/en); [Webinar slide 39](https://www.easa.europa.eu/en/downloads/139577/en).
- **Hover test heights:** HIGE 1.5 m, HOGE 30 m, flight and ground idle at 0 m — [Webinar slide ~48, test-condition table](https://www.easa.europa.eu/en/downloads/139577/en).
- **Test matrix:** speeds Vy, 0.66 Vy and 1.33 Vy at descent angles of 3, 6, 9 and 12°, at 110/120 m and 150/160 m heights, with at least 36 valid runs. Hover follows CAEP WG1 guidelines (CAEP12_WG1_3_IP03) — [Webinar slides 36-37](https://www.easa.europa.eu/en/downloads/139577/en).
- An earlier NLR ERF 2012 paper describes "Measurement of helicopter noise hemispheres utilizing a 100m vertical array" (Tuinstra & Sijtsma). It is cited as the way to reach higher lateral angles — [D1.5d refs [16] p.44](https://www.easa.europa.eu/en/downloads/132005/en).

### Inferences
- **Ground reflection removal** is done by measuring with **inverted microphones on ground plates** and applying the "free field" correction. That is the standard ground-plate pressure-doubling (about −6 dB) correction, which is my inference; the correction value is not stated. So there is no explicit image-source subtraction.
- The 1.2 m tripod microphones are used only for quality control against certification levels.
- **No de-Dopplerization or retarded-time back-propagation detail is given.** Scaling to 60 m presumably uses spherical spreading and ARP5534 absorption at the measured conditions, then re-adds absorption at reference conditions; that would explain FREEFIELD = 2.
- **Contrast with NICE-OPS's databases**, which are built by panam from predictions or other source spheres. NORAH's are purely measured, which explains the partial angular coverage and the absence of thrust as a parameter.

### Gaps
- The exact free-field (ground-plate) correction values.
- The emission-angle reconstruction (retarded-time geometry with DGPS) and the averaging rules across runs (energy versus dB) are only summarised. They are presumably detailed in the NORAH1 D3.2/D4 reports and in van Oosten et al. 2024, which were not obtained.
