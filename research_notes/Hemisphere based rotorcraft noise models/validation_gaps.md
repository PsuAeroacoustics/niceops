# Closing the validation gaps: hemisphere-based rotorcraft noise models

Scope: this note follows up the Gaps sections of `validation_limitations.md`. For each finding it says whether the number is a **measured** comparison against flight data, a **model-to-model** or **simulated** check, or an author **claim** with no number behind it. Full texts were read for:
- Watts et al. 2019 (NASA/TM-2019-220264)
- the EASA NORAH2 modelling guidance (D1.5d) and the final webinar
- ECAC Doc 32 (1st ed., 2026)
- Tuinstra and Sijtsma (NLR vertical array)
- Munt et al. 2001 (DERA/QinetiQ)
- Tuinstra and van Oosten (ERF 2018)
- Lucas and Marcolini (RNM)
- the ASCENT Project 038 annual reports for 2018 and 2025

The HELENA, CH-146/RNM and UQ papers were available as abstracts only. The NORAH2 D3 validation deliverable (Olsen et al. 2024) is already summarised in `validation_limitations.md` and is not repeated here, beyond what its webinar adds.

## 1. What does the 2017 Noise Abatement flight-test data report (Watts et al. 2019) contain?

### Takeaway
The report is NASA/TM-2019-220264, March 2019, by Watts, Greenwood, Smith and Stephenson. It is a data-description report, not a validation study: it reports no model-vs-measurement comparison and no repeatability statistics.

It documents:
- six light helicopters: R44, R66, B206L3, B407, AS350B3 (utility configuration, +75% drag) and EC130B4;
- 24 flight days between 15 Aug and 26 Oct 2017, about 1,450 test points;
- a 52-position flush ground-board array;
- tethersonde, temperature-string and LIDAR weather data;
- 50 Hz INS/GPS state data;
- per-run and averaged ART hemispheres at a 100 ft radius over 30–2,500 Hz A-weighted, flagged `DOPPLER_SHIFT_REMOVED`.

The B407 days NICE-OPS uses are flights 283–286 (10–13 Oct 2017) at Amedee AAF.

### Cited Findings
- **Identity and stated objectives** — [Watts, Greenwood, Smith, Stephenson, "Noise Abatement Flight Test Data Report," NASA/TM-2019-220264, Mar 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - The test was a joint NASA/FAA/US Army flight test.
  - Objectives: (1) source noise "to enable modeling of the tested aircraft in empirical models, like AAM, and semiempirical models, like FRAME"; (2) maneuvering-flight data for Fly Neighborly guidance; (3) low-noise approach data "for validation of theoretical models, such as the coupled PSU-WOPWOP/CHARM prediction method."
- **Aircraft (Table 1)** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):

  | Aircraft | TOGW (lbf) | MR blades | MR rotation | MR tip speed (ft/s) | Tail rotor |
  | --- | --- | --- | --- | --- | --- |
  | R44 | 2029 | 2 | CCW | 705 | conventional |
  | R66 | 2041 | 2 | CCW | 705 | conventional |
  | 206L3 | 3462 | 2 | CCW | 744 | conventional |
  | 407 | 3824 | 4 | CCW | 757 | conventional |
  | AS350B3 | 4247 | 3 | CW | 723.5 | conventional |
  | EC130B4 | 4529 | 3 | CW | 705 | Fenestron |

  - The AS350B3 flew in a utility configuration, with a work platform and cowlings that raised airframe drag by 75% over stock. It is therefore not a stock-airframe source.
- **Sites and dates (Table 13)** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - R44 and R66 flew at Eglin AFB Range C62 (15–25 Aug 2017), with static tests at Wagner Field.
  - B206L3, B407, AS350B3 and EC130B4 flew at Amedee AAF, Sierra Army Depot, CA (4–26 Oct 2017).
  - B407 flights: 283 (10/10/2017, 75 test points), 284 (10/11, 58), 285 (10/12, 48) and 286 (10/13, 58).
  - Test points per aircraft (summed from Table 13): R44 281, R66 222, B206L3 232, B407 239, AS350B3 191, EC130B4 287, for a total of about 1,452.
- **Microphone array** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - 52 positions were "essentially the same" at both sites. A "Turn Array" variant dropped mics 2–5 and added 53–56.
  - Primary sensor: G.R.A.S. 67AX, a 1/2-inch microphone flush-mounted in a 15-inch ground board, offset from centre ("based on SAE ARP4055, but adapted from inverted microphones to flush mounted").
  - Three 46AE microphones sat on 4 ft tripods at −45°, 0° and 45° to emulate certification placement.
  - Sample rates: 25,000 S/s (16-bit WAMS) and 25,500 S/s (24-bit WAARN).
  - Hover and idle used a 500 ft-radius semicircle of 13 inverted ground-board microphones, recorded by Volpe and *not distributed by NASA*.
- **Weather** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - A tethered sonde was held at the primary 500 ft AGL altitude, with a temperature string every 10 ft logged once a minute.
  - A ZephIR 300 LIDAR (13 heights to 900 ft) sat 2,300 ft before the PAPI. The distributed LIDAR file gives 10-minute averages at 984, 696, 495, 348, 197, 125 and 46 ft AGL.
  - Two to five 4 ft ground weather stations were spread among the array.
  - Where the balloon was below 500 ft, temperature above it was extrapolated at the standard lapse rate (−3.564 °F/1000 ft). The distributed temperature profiles are therefore partly synthetic aloft.
- **Vehicle state** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - The NASA ANTS unit uses a VectorNav VN-200 Kalman-filtered INS/GPS logged at 50 Hz: position, velocity, acceleration, attitude and rates.
  - No accuracy figure for ANTS is given in the report.
- **Test matrix** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - Level flight: 40 KIAS to VH, mostly at 200 ft AGL, with two certification speeds at 500 ft.
  - Descents: 3–12° from 60 KIAS to Vy, PAPI-guided, captured about 10,000 ft out.
  - Climbs: Vy ± 10 KIAS.
  - Turns: 20–45° bank, five techniques (constant speed, constant torque, accelerating roll-in, turn from climb/descent, decelerating roll-in).
  - Approach conditions were designed on site from the day-1 descent data, so there were no predefined approach codes. The approaches ended in a 15 ft hover taxi.
- **Hemisphere processing** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - Processing used ART within RNM, with spectra taken "typically every 0.5 seconds" and "depropagated to a hemisphere of a 100 foot radius, centered at the aircraft GPS receiver. Atmospheric and speed effects are accounted for."
  - Level-flight hemispheres were averaged over ground-speed groupings.
  - Descent hemispheres were *not* averaged: scatter in the flown condition (IAS plus visual PAPI) left "no obvious groupings."
  - Plots cover 30–2,500 Hz A-weighted, a band chosen "to minimize the effect of background noise because these are quiet aircraft and the distances … are large."
  - Sideslip descents were flown for the AS350B3 and EC130B4.
- **Hemisphere file schema (Table 24)** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf):
  - Fields: BB/NB/PT switches, `DOPPLER SHIFT REMOVED` (0/1), `RADIUS` (ft), target `FLIGHT PATH ANGLE`, `SPEED` (*ground* speed, kts), `PHI`, `THETA`, `FREQUENCY`, `AMPLITUDE`.
  - Weight, mast tilt and pylon fields are "not used for this test."
- **Distribution, conclusions and stated validation intent:**
  - The data are "publicly available upon request." The conclusions state that "maneuvering flight can lead to significant increases in noise for all aircraft tested" and that low-noise approaches were demonstrated for all six. — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf)
  - FRAME "must be calibrated using measured steady-flight acoustic data for each type of helicopter," which is why the steady matrix was flown. — [Watts et al. 2019, §2](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf)

### Inferences
- The report confirms the 2017 data were designed partly *for* model validation. Its only validation customer named in the text is the physics-based PSU-WOPWOP/CHARM system, so a sphere-model validation on these data (NICE-OPS) would be new use rather than a replication.
- Three properties of the distributed hemispheres matter to any sphere tool built on them:
  1. Speed is **ground speed**, and the condition is a *target* FPA. This is the database-side counterpart of NICE-OPS's `--speed_reference` and its air-referenced wind correction.
  2. The radius is 100 ft (30.5 m). NORAH and ECAC use 60 m.
  3. Spectra stop at 2,500 Hz, so any A-weighted level from these spheres omits content above that.
- Because the descent hemispheres are unaveraged single runs, any blend of them carries single-run variability. The ±0.5 to 7.5 dB repeatability floors in the earlier note apply in full.

### Gaps
- The report gives no INS/GPS accuracy, no microphone calibration uncertainty, no ground-board frequency response, and no run-to-run repeatability statistics. Those must come from other papers or from the data themselves.
- The cited ART depropagation reference ("Ref. 13") was not retrieved. Its treatment of ground reflection and of refraction is therefore not confirmed from the primary source.

## 2. What validation exists for the European hemisphere tools (NORAH1, NORAH2, HELENA, SELENE, HELIACT, ECAC Doc 32)?

### Takeaway
No European hemisphere tool has a published, quantified footprint validation against in-service operations:
- NORAH1 and NORAH2 are validated against certification levels and Nord2000 only.
- HELENA's EC130/EC135 comparison is described in its abstract without numbers.
- DLR SELENE's EC135 validation is a qualitative "good agreement."

The new regulatory method, **ECAC Doc 32 (1st ed., endorsed 19 May 2026)**, is built on NORAH but deliberately *drops* hemisphere interpolation ("no agreed method"). It forbids nothing, but disclaims single-event accuracy.

### Cited Findings
- **ECAC Doc 32 on interpolation and scope** — [ECAC.CEAC Doc 32, "Report on Standard Method of Computing Rotorcraft Noise Contours," 1st ed., endorsed by DGCA/166 on 19 May 2026](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf):
  - "At present there is no agreed method for the interpolation between hemispheres to reflect flight conditions not directly covered by specific hemispheres. Consequently, users will need to select the most appropriate hemisphere to assign to each flight path segment."
  - The method "currently only supports flight paths consisting of segments of constant speed and constant climb/descent angles" matching an available hemisphere. Interpolation "is expected in a future edition."
  - It "is primarily designed to estimate long-term average noise exposure"; for single movements "the uncertainties will be greater, and it should not be relied upon for that purpose."
  - "Refraction … is not covered by the present method."
- **ECAC Doc 32 hemisphere conventions** — [ECAC Doc 32](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf):
  - 60 m reference radius, 10° bins, one-third-octave bands from 10 Hz to 10 kHz, and absorption at the ICAO reference atmosphere.
  - Minimum required coverage is −60° ≤ φ ≤ 60° and the 10-dB-down polar interval.
  - Within a hemisphere, bilinear interpolation is done **in energy** (10^(L/10)). Beyond the data, constant-value extrapolation is taken from the nearest filled bin by great-circle angle, with an energetic average when several bins tie.
  - Recommended database density: descent-angle intervals of at most 3° at four velocities, and level flight at 0.9 VH, +10 kt, −15 kt and −30 kt.
  - Over the R22 flight envelope, levels vary "up to 10 dB … [and] generally … within 1 dB as a function of rate of climb."
  - Radar track input is accurate to 5–10 m; one-second rounding adds 10–20 m of position error.
- **NORAH2 interpolation (EASA D1.5d)** — [EASA NOISE-SC03 D1.5d "Rotorcraft noise modelling guidance"](https://www.easa.europa.eu/en/downloads/132005/en):
  - Speed and FPA are normalised by the database range, and the FPA is multiplied by a scaling factor F_fc = 2 "empirically determined to give minimum interpolation error"; no error value is reported.
  - A Delaunay triangulation is built over the conditions. Inside the convex hull, the three vertex hemispheres are combined as an **inverse-distance-weighted energy average**: L̂ = 10 log10(Σ 10^(L_j/10)/δ_j ÷ Σ 1/δ_j). This is not barycentric.
  - Outside the hull, the nearest neighbour is used.
  - Hemispheres are looked up by the 3-D velocity V_A, "not V_g." Velocity is ground-relative and "neglects wind."
  - Turns tilt the closest level-flight hemisphere by the bank angle. For descending turns this is again a *level-flight* hemisphere.
  - Hover offsets (+12, −12 and −2.5 dB) are given with the caveat that "EMPA found that these values may not be valid for other microphone setups (e.g., 4m height over hard or soft ground)."
- **NORAH2 final webinar** — [EASA NORAH2 final dissemination webinar, 13 Mar 2024](https://www.easa.europa.eu/en/downloads/139577/en):
  - It frames NORAH1's no-interpolation design as trading "Position error" and "Hemisphere condition mismatch" for accuracy at the database points.
  - Validation is listed as: source vs certification (measurement and simulation), method by "peer review" and comparison with Nord2000, and peer implementation.
  - It recommends to "broaden validation of the method (… NORAH predictions vs measured noise in an operational environment)."
  - No interpolation-error figure appears.
- **NORAH1 hemisphere construction** — [Tuinstra & van Oosten, "The development of a European helicopter noise model," ERF 2018](https://dspace-erf.nlr.nl/items/c70f8916-c3ad-4c15-8bb4-fb7d67a88225/full):
  - Eight types; more than 170 test conditions, 800 runs and 60 flight hours.
  - "Helicopter attitude is not considered": the body x-axis is assumed aligned with the local track axis.
  - Repeat runs are merged by **energy summation** per bin, with a ground-plate pressure-doubling correction G.
  - Runs are rejected unless 90% of velocity samples lie within ±5 kt of nominal and 90% of positions within γ ± γ/12.
  - The paper contains no measured-vs-predicted validation numbers.
- **HELENA, FRIENDCOPTER (Gervais et al.)** — [Gervais et al., "Validation of EC130 and EC135 Environmental Impact Assessment Using HELENA," AHS 66th Forum, 2010 (abstract only)](https://consensus.app/papers/details/9ec15701cb4b5ca5ac9fd637a70797e0/?utm_source=claude_desktop):
  - Steady single-segment procedures were compared with measurements on both aircraft, "in terms of commonly used integrated noise metrics, as well as detailed time-histories."
  - HELENA was compared with an NPD model.
  - A multi-segment procedure was computed "preliminar[ily]" against maneuvering flight data.
  - No dB values are in the abstract.
- **DLR SELENE, EC135** — [Guntzer, Spiegel, Lummer, "Genetic Optimizations of EC135 Noise Abatement Flight Procedures Using an Aeroacoustic Database," 2009 (DLR elib)](https://elib.dlr.de/63145):
  - Optimised procedures were flown in 2008, and "a good agreement between simulated and measured noise levels is observed."
  - The only number is a *predicted* reduction of 10 dB SEL between 3 and 1 km before landing. This is a claim; no error statistic is given.
- **HELIACT (DERA/QinetiQ, UK)** — [Munt, Browne, Pidd, Williams, "A measurement and prediction method for determining helicopter noise contours," 27th ERF, Moscow, 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full):
  - A UK hemisphere tool built from vertical beamforming arrays (75 m and 200 m from the track, 28 microphones, 5° resolution over 60° ≤ |φ| ≤ 90°).
  - The authors note that for such tools "the accuracy has not been fully validated."

### Inferences
- The European regulatory line has moved *backwards* on interpolation. NORAH2 (2024) introduced energy-domain inverse-distance blending; ECAC Doc 32 (2026) did not adopt it. A tool that interpolates (NICE-OPS: barycentric blending over μ, γ_e, CT) is therefore outside the endorsed European method. The absence of any published interpolation-error measurement is the reason ECAC gives.
- The NORAH lineage differs from NICE-OPS in ways that matter for comparison:
  - NORAH ignores attitude in hemisphere construction; NICE-OPS uses a TPP-normal sphere frame.
  - NORAH uses ground-relative velocity; NICE-OPS uses a ground or air reference according to the database.
  - NORAH tilts *level* hemispheres for turns, including descending turns; NICE-OPS uses Q-SAM with CT and γ_e.
  - NORAH interpolates in energy; NICE-OPS blends dBA linearly.

### Gaps
- HELENA's EC130/EC135 dB results remain unretrieved; only the abstract is available.
- No NORAH1 validation against operational measurements was found. A search for Zurich and Frankfurt operational comparisons returned nothing.
- A search snippet stated that "variations in the measured sound exposure level from one approach to another are on the order of 3–7.5 dB(A)" for EC135 tests, with geometric accuracy of about 6 m horizontal and 12 m vertical. The originating document could not be identified among the results (DLR elib 98261, SAE F-0070-2014-9418), so the figure is unattributed and should not be cited without retrieval.
- The NORAH2 webinar was cited by its EASA download number (139577), inferred from the downloaded file. The landing page itself was not opened.

## 3. What quantified sphere-model (RNM/AAM) and physics-model comparisons against US flight test exist?

### Takeaway
The quantified US comparisons fall into four groups:
- **Early RNM self-consistency.** XV-15 level flyovers came "generally within 0.5 dB," with excursions blamed on the two-point trajectory input. The CH-146 abstract claims agreement against an independent NASA array, without numbers.
- **PSU physics-based predictions on the 2017 data.** These are within about 2 dBA of peak SEL for five of six aircraft in an 80 kt, 6° BVI descent. The EC130 came out 5–6 dB low, because the Fenestron shroud was not modelled.
- **A few-point check of a single-run sphere** at three microphones: within 1.2 dB, with trajectory and meteorology changing the prediction by at most 0.2 dB.
- **Joby flyovers and the Eglin 2007 campaign.** For Joby, no AAM-vs-measured footprint comparison was found. For Eglin 2007, no validation numbers were found.

### Cited Findings
- **RNM original validation (XV-15)** — [Lucas & Marcolini, "Rotorcraft Noise Model," Wyle/NASA LaRC (NTRS 20040110391; c. 1997)](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf):
  - For XV-15 constant-speed level flyovers at 70, 90, 110 and 140 kt (60° nacelle), SEL vs sideline: "Generally, the noise predictions are within 0.5 dB, but there are several instances where the sound levels are greater."
  - The authors attribute the discrepancies to the trajectory, which was modelled by two endpoints only.
  - The hemispheres came from the same test via ART, so this is closer to a self-consistency check than independent validation.
- **RNM interpolation scheme** — [Lucas & Marcolini](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf):
  - Within a hemisphere, "Coons" surface interpolation in θ, φ.
  - Between hemispheres, airspeed only, assuming "a linear variation in pressure squared with airspeed," i.e. energy-domain interpolation.
  - FPA and pylon angle are taken as nearest, not interpolated, and a user keyword can switch the speed interpolation to nearest as well.
  - Propagation adds spreading, ANSI S1.26 absorption, ground, wind and topography terms.
- **RNM/ART2, CH-146 (Bell 412), Moose Jaw 1998 (NATO CCMS)** — [Page et al., "Rotorcraft noise model (RNM) and acoustic repropagation technique (ART2) validation and application using CH-146 (Bell 412) measurement data," JASA 2000 (abstract)](https://consensus.app/papers/details/822703a2f7cb5b5abb8234068399f832/?utm_source=claude_desktop):
  - USAF-array data built the hemispheres, which were propagated back to the measurement locations.
  - "Agreement with the original USAF measurements demonstrated consistency of the ART/RNM process. Agreement with the independent NASA measurements provided validation."
  - No dB values are given in the abstract. This is the only located case of hemispheres validated against an *independent* array.
- **Crane-suspended vertical arrays (CH-46E, CH-53E, AH-1)** — [Lucas, "Acoustic directivity measurements of the CH-46E, CH-53E, and AH-1A helicopters," JASA 1998 (abstract)](https://consensus.app/papers/details/0fe21bb4d4f85cdd9f0fbb24c0386f6f/?utm_source=claude_desktop):
  - Two 180 ft crane arrays plus a ground array formed a "U," with the rotor plane flown roughly in line with the crane tips.
  - This is the RNM-era answer to near-rotor-plane coverage.
- **PSU-WOPWOP/CHARM vs the 2017 data (not sphere-based; same data NICE-OPS uses)** — [ASCENT Project 038 2018 Annual Report, Task 10](https://ascent.aero/documents/2019/07/ascent-project-038-2018-annual-report.pdf):
  - Case: 80 kt, 6° descent ("strong BVI").
  - "Peak SEL values and the contour shapes agree within about 2 dBA or less for all of the aircraft, except for the EC130," which was "lower than the measured data by 5 to 6 dB," attributed to the unmodelled Fenestron shroud.
  - The rotation-direction contour skew (right for CCW rotors, left for CW) was reproduced.
  - Post-processing was identical for prediction and data.
  - The comparison was made *before* the broadband ground-reflection fix. With that fix, the R66 (broadband-dominated) moved much closer, and the B407 became "slightly overpredicted."
- **Sphere-based few-point check, R44 run 227164 (3°, 60 kt)** — [ASCENT 038 2018, Table 11.1](https://ascent.aero/documents/2019/07/ascent-project-038-2018-annual-report.pdf):

  | Mic | Measured | Predicted (four ideal/actual trajectory × std/actual meteo variants) |
  | --- | --- | --- |
  | 19 | 95.50 | 95.6–95.8 |
  | 23 | 83.79 | 83.1–83.3 |
  | 15 | 82.90 | 81.7–81.9 |

  - Every variant used "Sphere #164", the hemisphere of the same run, so this is a self-consistency check.
  - Swapping the ideal trajectory for the actual one, or standard meteorology for actual, moved predictions by at most 0.2 dB.
  - The table is marked "preliminary," and which microphones built the sphere is "TBD."
- **PSU later work** — [ASCENT Project 038 2025 Annual Report](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-038-2025-Annual-Report.pdf):
  - Early comparisons were "within several dBA of the flight test data."
  - A Pegg broadband-model error (C̄_L overestimated about 1.5x) caused 1.67 to 13.39 dBA overprediction depending on C̄_L.
  - EC130 HOGE hover vs 2017 measurements: "within 2–4 dB at several microphone locations" after modelling the Fenestron shroud.
  - Neglecting shroud scattering can give "discrepancies of up to 10 dBA" (citing Falissard et al. 2011).
  - Preliminary PSU-to-AAM coupling was implemented in 2022.
- **2017 RNM/ART use of the data** — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf): the hemispheres were produced with RNM's ART, but the report contains no RNM footprint prediction compared with the measured approach or turn contours.

### Inferences
- Every quantified "within 0.5 dB" sphere result found (XV-15, the R44 table) predicts the *same run* that built the sphere. Those results measure processing round-trip consistency, not predictive skill. Independent-run validation of sphere models (CH-146) exists only as an unquantified abstract claim.
- The PSU results give a benchmark for the same 2017 data a sphere tool uses: about 2 dBA peak-SEL agreement for a steady BVI descent, from physics alone. A held-out NICE-OPS approach error of about +1.5 dB median (internal, see the earlier note) is in the same range. NICE-OPS is evaluated on decelerating approaches, a harder case.
- Trajectory and meteorology moved predictions by at most 0.2 dB at near-field mics in the R44 table. This suggests source-state selection dominates near the track. Refraction matters mainly at the long ranges where NICE-OPS's rim-dependent error sits.

### Gaps
- No quantified Eglin 2007 (Acoustics Week III) or Eglin 2010 sphere-model validation was found. The Eglin 2007 document (govdocs mirror) could not be downloaded.
- No AAM prediction of the Joby measured flyovers against the measured footprint was found. The Joby literature reports hemispheres and repeatability only.
- Volpe's 2019 in-situ Fly Neighborly report and the Page/Rapoza/Jacobs Inter-Noise 2021 paper both returned HTTP 403 and were not read.
- The FAA Heliport Noise Model (HNM) is NPD-based rather than hemisphere-based. Its validation was not pursued.

## 4. Quantified evidence on specific sphere-processing error sources

### Takeaway
Several of the gaps named earlier now have numbers, mostly from analysis or simulation rather than flight-test validation:
- **dB vs energy averaging.** Averaging in dB biases the estimate low by 0.03, 0.15 and about 0.6 dB for ±1, ±2 and ±6 dB scatter (NLR, analytical). An energy-domain least-squares fit gave grossly wrong results under strong directivity.
- **Ground-reflection removal at an elevated microphone.** A standard reflection correction errs by nearly 10 dB on rotor tones (DERA, simulated), against ±1 dB for a beamforming vertical array.
- **Refraction during de-propagation.** It costs about 2 dB at 40° < θ < 140° near the rotor plane (φ = 80°), with "very large randomly varying errors" outside that range. This motivated an r < 450 m limit.
- **Doppler.** De-Dopplerised flight data and moving-observer predictions agree "typically within 0.5 dB."
- **Flush ground-board measurement.** It is only specified for incidence within 30° of normal (ARP4055), which the rim of every hemisphere violates.

### Cited Findings
- **SPL-based vs energy-based averaging and fitting** — [Tuinstra & Sijtsma (NLR), "Measurement of helicopter noise hemispheres utilizing a 100 m vertical array," ERF paper (AH-64 trials; 93-mic vertical array + 13 ground mics)](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download):
  - Direct dB averaging converges not to the energy mean but to a bias E3 ≈ 2.171 (σ_Z/p²_rms)².
  - The resulting "highest accuracy" is 0.03 dB for weak (±1 dB), 0.15 dB for medium (±2 dB) and about 0.6 dB for strong (±6 dB) scatter. *Analytical/statistical result, not flight validation.*
  - With fewer than 50 averages per hemisphere point (typical in flyovers), the finite-sample energy estimator is no better, so "a regression model based on SPL is justified."
  - Energy-based least-squares directivity fits "produce erroneous results" at 63 Hz, including "negative acoustic energies," because a global least-squares fit minimises absolute rather than relative error under strong directivity.
- **Ground-reflection removal for elevated microphones** — [Munt et al., 27th ERF 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full):
  - Setup: a simulated hover at 35 m AGL, 220 m from the array, main-rotor thickness tones only, over perfectly reflecting ground.
  - Standard reflection-correction theory applied to a single elevated microphone gave a direct-field error "nearly 10dB."
  - A beamforming vertical array recovered the direct field to "better than ±1 dB." The low-frequency array was better than ±2 dB over 20–160 Hz.
  - For a highly elevated microphone, direct and reflected waves leave the source at different times and angles, so there are "two unknown source levels … the problem is undetermined."
- **NLR beamforming result** — [Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download):
  - Delay-and-sum over six-microphone clusters removes the ground reflection well up to 200 Hz. There are "discrepancies" from 200 to 500 Hz, and the method was not used above 500 Hz (1 m spacing).
  - The authors recommend avoiding ground microphones for lateral angles beyond 60°.
- **Refraction contaminating de-propagation** — [Munt et al. 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full):
  - Setup: MAPE parabolic-equation runs (100 Hz, 5 m/s wind, 100 m source over grass, turbulence-averaged).
  - Assuming a uniform atmosphere during de-propagation gives "errors of the order of 2dB … in the range 40° < θ < 140° and very large randomly varying errors … outside this range" at φ = 80° (the microphone 573 m to the side).
  - At φ = 60° the uniform assumption holds only for 25° < θ < 155°.
  - Upwind and downwind attenuation differ by at least 15 dB at the θ extremes.
  - Conclusion: propagation distance should be limited to r < 450 m at the edge of certification meteorology. *Model-based, not measured.*
- **Far-field distance and resolution trade-off** — [Munt et al. 2001](https://dspace-erf.nlr.nl/items/d80bbadf-b1e8-4e15-8c4e-e91d683b4814/full):
  - Measurements assume a point source with R ≫ vehicle size.
  - At 120 kt and 5° angular resolution, a 1 Hz narrowband needs r ≥ 686 m, while 3.5 Hz needs r ≥ 196 m. Long ranges improve tonal resolution but invite refraction error.
  - No quantification of a near-field error from hemisphere radius was given.
- **Doppler convention** — [ASCENT Project 038 2025 Annual Report (summarising the 2020 report)](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-038-2025-Annual-Report.pdf):
  - "De-Dopplerization (used in flight test data processing) and moving observers (used in noise predictions to eliminate Doppler effects) were demonstrated to be effectively equivalent (typically within 0.5 dB or less)." *Measured-vs-predicted processing comparison.*
  - MOTUS (Airbus) builds its correction database from "de-Dopplerized noise hemispheres." — [same report](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-038-2025-Annual-Report.pdf)
- **Ground-plane microphone validity** — [SAE ARP4055, "Ground-Plane Microphone Configuration for Propeller-Driven Light-Aircraft Noise Measurement"](https://saemobilus.sae.org/standards/arp4055-ground-plane-microphone-configuration-propeller-driven-light-aircraft-noise-measurement):
  - Equivalent free-field levels are obtainable from 50 Hz to at least 10 kHz "for sound incidence angles within 30 degrees of the normal to the ground."
  - The 2017 array used flush mics "based on ARP4055, but adapted from inverted microphones to flush mounted." — [Watts et al. 2019](https://aeroacoustics.psu.edu/files/2022/09/Watts.TM2019-compressed.pdf)
- **Upper-hemisphere and rim coverage in practice:**
  - NORAH/ECAC require only |φ| ≤ 60° plus the 10-dB-down polar band and fill the remainder by nearest-bin constant extrapolation. — [ECAC Doc 32](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf)
  - "Measurement of higher lateral angles or polar angles approaching 0°/180° is possible but requires complex measurement setups." — [ECAC Doc 32](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf)
  - Opposite-rotation types in a class are handled by mirroring azimuth (the "±" in the class formula). — [EASA D1.5d](https://www.easa.europa.eu/en/downloads/132005/en)
- **Sparse-array interpolation on the sphere** — [Greenwood, "A Physics-Based Approach to Characterizing Helicopter External Noise Radiation from Ground-Based Noise Measurements," Univ. of Maryland thesis, 2008 (abstract)](https://drum.lib.umd.edu/items/71b7b329-0778-4cd5-8b65-e5c1d09cfe6c/full):
  - "Current methods of interpolation of data on RNM acoustic radiation spheres are found to lead to high levels of inaccuracy when using sparse microphone arrays."
  - Radial basis functions are "shown to lead to significantly improved accuracy."
  - No numbers are in the abstract.
- **Ground reflections in hemisphere back-propagation (drones)** — [Herold et al., AIAA/CEAS Aeroacoustics 2024 (abstract)](https://consensus.app/papers/details/ec9f283cd5ad5631bfe32ccc9a95a44f/?utm_source=claude_desktop): back-propagation "followed by energetic averaging," with "the influence of ground reflections" named as a current limitation (unquantified).
- **Weather sensitivity of hemisphere-based footprints (eVTOL, model-to-model)** — [Fuerkaiti et al., JASA 2021 (abstract)](https://consensus.app/papers/details/b05b0c351cd55240bdf8f87bd18ec91f/?utm_source=claude_desktop): vertical temperature and wind gradients cause a "mismatch up to 4 dBA in the illuminated zone and a significant drop in the refractive shadow zone" relative to a homogeneous atmosphere.

### Inferences
- **dB-domain blending in NICE-OPS.** NICE-OPS blends the broadband dBA linearly in dB across the four Delaunay-bracketing conditions and across bilinear sphere corners. RNM, NORAH2 and ECAC all interpolate in energy. Tuinstra and Sijtsma's statistics suggest the difference is small where neighbouring values differ by a couple of dB: averaging in dB sits below the energy mean by roughly 0.15 dB at ±2 dB spread.
  - Their scatter figures (±1, ±2, ±6 dB) are reused here as a stand-in for the spread between blended neighbours, with weights instead of equal averages. That reuse is mine, not theirs.
  - Under the same reading, the bias grows toward 0.6 dB or more where the blended neighbours differ by about 6 dB or more, as across a BVI boundary or between a hover sphere and a forward sphere.
  - That direction (dB blending lower than energy) is opposite in sign to the NICE-OPS slow-end over-prediction, so it cannot explain it.
- **Ground boards.** The flush ground-board, single-height arrangement of the 2017 and 2019 arrays avoids the "nearly 10 dB" elevated-microphone reflection problem at near-normal incidence. ARP4055's 30° validity cone does not cover the near-horizon rim data, which are exactly the long-range emission angles where NICE-OPS's internal error depends on how the rim is treated.
- **Refraction in the database itself.** Munt's r < 450 m and θ-range limits imply that hemisphere data near the rim are themselves refraction-contaminated, at about 2 dB inside the limits and worse outside, because ART de-propagates through a uniform (or averaged) atmosphere. A propagation model applied downstream (NICE-OPS `--refraction`) cannot undo refraction already baked into the database. Validating long-range refraction against those same rim data is partly circular.
- **Doppler.** De-Dopplerised flight hemispheres used with a moving-source, non-Doppler-shifting propagator are consistent to about 0.5 dB, per the PSU equivalence result. That supports NICE-OPS's convention of not applying a frequency shift, given the database's `DOPPLER_SHIFT_REMOVED` flag. The internal memory note records that the panam spheres keep Doppler, so the flag in each file is the thing to check, not the convention.

### Gaps
- No flight-test measurement of the dB- vs energy-domain *interpolation between operating conditions* was found, in any tool. NORAH2's F_fc = 2 was "empirically determined to give minimum interpolation error" without a published error.
- No quantification of hemisphere radius (30 m vs 60 m vs larger) or near-field error for helicopter spheres was found. Munt gives only the far-field/resolution trade-off.
- No quantified error from mirroring (opposite-rotation azimuth flip, or lower-to-upper completion) was found anywhere.
- Greenwood's 2008 RBF-vs-current-method accuracy numbers are in the thesis body, which was not read.

## 5. Uncertainty quantification of hemisphere-based footprints

### Takeaway
No Monte-Carlo or Bayesian UQ of a hemisphere-based footprint against flight test was found. The UQ literature that exists sits elsewhere:
- **Physics-model inputs:** wake-model uncertainty gives up to 3.5 dB at BVI-dominated observers; gust uncertainty in hover loading noise has also been studied.
- **Fixed-wing system-noise simulation:** single-event UQ.
- **Model-to-model weather sensitivity:** up to 4 dBA in the illuminated zone.

The regulatory method (ECAC Doc 32) explicitly declines single-event accuracy claims.

### Cited Findings
- **Tip-vortex core model uncertainty** — [Vouros et al., JAHS 2021 (abstract)](https://consensus.app/papers/details/966005587576505089908f24ba8ea695/?utm_source=claude_desktop): uncertainty in the core model propagates to "up to 25% and 3.50 dB" at BVI-dominated observers when no calibration data exist (polynomial chaos; simulation only).
- **Gust-induced hover loading noise** — [Gill et al., JAHS 2026 (abstract)](https://consensus.app/papers/details/cd857453e2695f0bbde455e0f17590a5/?utm_source=claude_desktop): UQ by polynomial chaos expansion, validated against Monte Carlo. Uncertainty is largest on the rotor axis, and gust length dominates at higher harmonics.
- **Single-event UQ including the trajectory (fixed-wing)** — [Bertsch et al., J. Aircraft 2025 (abstract)](https://consensus.app/papers/details/c6162eb544785edb94de1b0568ec0301/?utm_source=claude_desktop): a fixed-wing (supersonic transport) single-event study that extends UQ "to fully incorporate the simulation of the flight trajectory." It is the closest methodological template found; no rotorcraft equivalent was found.
- **Regulatory stance** — [ECAC Doc 32](https://www.ecac-ceac.org/images/documents/ECAC-CEAC-DOC_32-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_ROTORCRAFT_NOISE_CONTOURS.pdf):
  - The method is for long-term averages; for single events "the uncertainties will be greater, and it should not be relied upon."
  - It quantifies only input-track uncertainty (5–10 m radar, 10–20 m with time rounding).
- **Repeat-run acceptance windows as an implicit source-variability bound** — [Tuinstra & van Oosten, ERF 2018](https://dspace-erf.nlr.nl/items/c70f8916-c3ad-4c15-8bb4-fb7d67a88225/full): NORAH merges runs only when 90% of samples lie within ±5 kt and γ ± γ/12. No resulting dB spread is reported.
- **Snapshot-array design for capturing within-run variability** — [Houston et al., "Snapshot Array Design Considerations for Rotorcraft Noise Characterization," VFS 81st Forum 2025 (abstract)](https://consensus.app/papers/details/be5ccb59d81053e1bf27b6d584953c9c/?utm_source=claude_desktop); [Stephenson et al., "Overview of a Comprehensive MD530F Acoustic Flight Test," VFS 81, 2025 (abstract)](https://consensus.app/papers/details/e5a0b6439a705c3eb34dcbae797542c3/?utm_source=claude_desktop):
  - The July 2024 Army/NASA MD530F test used a uniform-coverage snapshot array, embedding lower-count subsets to study channel count.
  - Its stated motivation: vehicles "may not maintain a steady aerodynamic state during flight, making source hemispheres measured using traditional linear arrays unreliable."
  - Conditions were chosen by one-factor-at-a-time and three design-of-experiment approaches, which is an emerging statistical basis for database UQ.
- **Weather (model-to-model)** — [Fuerkaiti et al. 2021](https://consensus.app/papers/details/b05b0c351cd55240bdf8f87bd18ec91f/?utm_source=claude_desktop): up to 4 dBA in the illuminated zone.

### Inferences
- A defensible footprint error bar for a sphere tool would have to combine three pieces, none of which any published sphere tool has propagated:
  1. measured source repeatability (±0.5 dB at best in level flight, up to 7.5 dB for BVI);
  2. database processing bias (about 0.1–0.6 dB from dB averaging, about 2 dB of refraction contamination near the rim);
  3. propagation and weather (up to 4 dBA illuminated, more in shadow).
- NICE-OPS already holds the ingredients for a first such UQ: repeat runs in the 2017 set, measured weather profiles, and multiple propagation models. That would be novel relative to the literature found.
- The 2024 MD530F design-of-experiments matrix is the first sign of a statistically designed sphere database. A NICE-OPS blending-error study could use it to separate interpolation error from source variability.

### Gaps
- No Monte-Carlo or Bayesian UQ of a hemisphere-based helicopter or eVTOL footprint against flight test was found, in Consensus or in web search.
- No published study combines source-state uncertainty (from wind, through air- versus ground-referenced speed and the effective flight-path angle) with propagation uncertainty.
- The full texts of the 2025 VFS snapshot-array papers were not read, so their channel-count error numbers are unknown.
