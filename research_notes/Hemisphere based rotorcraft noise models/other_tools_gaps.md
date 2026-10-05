# Closing the gaps in other_tools.md: HNM/AEDT, sonAIR, HELENA and the European tools, RNM 7/APET/AMAT, and recent eVTOL footprint tools

Scope note for the report writer:

- **Purpose:** this file adds to `other_tools.md` and does not replace it. It answers that file's Gaps sections where a source could be found, and says so where none could.
- **Primary documents read in full:**
  - The FAA HNM methodology report (1988). The PDF is a scan; I OCR'd it, so equations are paraphrased rather than quoted.
  - ACRP 02-44 *Guidance for Helicopter Community Noise Prediction* (Sept 2015).
  - Lucas & Marcolini's RNM paper (c. 1997).
  - The ANOPP2 AMAT slides (Dec 2022).
  - Letica & Rizzi, NASA/TM-20250008619 (Sept 2025).
  - The ASCENT Project 084 2025 annual report (Apr 2026).
  - A TU Delft MSc thesis on sonAIR (2023).
  - The Empa F-35A presentation (Dec 2023).
  - The Swiss BAFU programme recommendation (updated Feb 2026).
  - Gennaretti et al. (IJA 2017).
  - Moses et al. (DAGA 2025).
  - Guntzer et al. (ERF 2023), already in other_tools.md; re-read for CAROT.
- **Still blocked:**
  - Empa pages and the sonAIR documentation PDF (WAF 403, also with a browser user agent).
  - MDPI (403).
  - The AEDT Technical Manual itself. It was not found online; AEDT facts below come from NASA's 2025 TM, which cites AEDT 3g.
- **Marking:** **[snippet only]** means the fact rests on a search-engine or abstract summary, not on a document I read. **[OCR]** means it was read from my OCR of a scanned page.

## 1. FAA HNM and the INM/AEDT helicopter mode: directivity, static and dynamic modes, lateral attenuation, ground, and how the NPD data are made

### Takeaway
HNM, which INM and AEDT's helicopter mode inherited, represents a helicopter with no hemisphere at all. It uses:

- **Moving modes:** SEL-versus-distance tables at eight distances (200 to 10,000 ft), for three moving modes in three lateral directions (left 45 deg, centre, right 45 deg).
- **Static modes:** LAmax-versus-distance tables for four static modes, each for hard and soft ground. A static mode's azimuthal directivity is a least-squares Fourier fit in azimuth (five cosine and three sine terms), derived from hover-in-ground-effect measurements.
- **Adjustments:** empirical terms for the advancing-tip Mach number, a 1/V duration term, a helicopter-specific low-angle lateral term, and a 90 deg-dipole "noise fraction".

AEDT 3g keeps this structure. It has 16 operational modes (11 of them static), the in-ground-effect height is hard-coded at 1.5 rotor diameters, and the static directivity is referenced to 200 ft. NASA's 2025 recommendations say that its extrapolation beyond +/-45 deg and its fourth-power noise fraction are wrong for helicopters and eVTOL.

### Cited Findings

**HNM methodology (FAA, draft report, 30 April 1988)** — all from [HNM methodology draft report, ROSA-P dot/10018 (PDF dot_10018_DS1.pdf)](https://rosap.ntl.bts.gov/view/dot/10018) **[OCR]**

- **HNM Version 1 contents:**
  - A database of 21 helicopters.
  - "Noise curves for hard and soft surface static modes and for left, right and center in-flight modes".
  - Algorithms "for hard and soft surfaces, duration, directivity, source noise, lateral effects, and noise fractions".
  - Its output metric is DNL.
- **Static modes:**
  - The four stationary modes are ground idle, flight idle, hover in ground effect (HIGE) and hover out of ground effect (HOGE).
  - Each is a table of maximum A-weighted level at eight distances, 200 to 10,000 ft, with separate tables for hard ground (pavement) and soft ground (grass).
  - "The data base also contains coefficients to define the acoustical emission directivity pattern for each static mode."
- **Moving modes:**
  - Takeoff, landing and level flyover, stored as SEL at the same eight distances.
  - Level flyover is stored at one airspeed; other speeds come from "a second order regression equation on main rotor tip speed".
- **Data provenance:**
  - FAA flight tests of 1976, 1978 and 1980, plus CERL TR N-38. They covered 6 deg approaches, takeoffs and 500 ft level flyovers, with some speed and altitude variations.
  - Reduction followed FAR-36 practice: 24 one-third-octave bands every 0.5 s.
  - SEL was integrated over the 10 dB-down duration and corrected to 77 F / 70% RH.
  - Where signal-to-noise was below 5 dB above 1.25 kHz, the spectrum was extrapolated at -2 dB per third octave.
- **In-flight lateral directivity:**
  - The left and right curves are "normalized to an elevation angle of 45 degrees", and the centre curve to 90 deg.
  - Between 45 deg and overhead the level is "determined by interpolation between the reference elevation angles" (linear in elevation angle).
- **Lateral (low-angle) effect:**
  - Between 45 and 15 deg elevation a helicopter-specific linear adjustment in elevation angle is added to the 45 deg SEL.
  - Below 15 deg a separate linear function of elevation angle is used.
  - The adjustment combines "excess ground attenuation and source directivity that are helicopter specific", and "unlike fixed-wing aircraft, the adjustment is not necessarily an attenuation".
- **Source-noise correction:** used for level flyover only, as a second-order polynomial in the change of advancing-tip Mach number M_adv = (1.689 V + Omega R)/c, with helicopter-specific coefficients.
- **Speed-duration term:** 10 lg(V_ref/V).
- **Noise fraction:**
  - Turns and finite segments use a noise fraction "based on a 90 degree dipole model" of the time history. The database assumes a straight, infinitely long flyby.
  - It does not apply to static operations.
- **Static directivity:**
  - It is "an adjustment to the sound level as a function of the helicopter azimuth angle which is measured clockwise from the nose".
  - The function is "a least square error fit using five cosine terms and three sine terms", that is, Sum_{n=0..4} a_n cos(n theta) + Sum_{m=1..3} b_m sin(m theta), from the report's Figure 4.
  - Coefficients are given per helicopter for hard and soft surfaces in its Appendix C.
  - It should be used only if the helicopter is oriented consistently all year; otherwise "omnidirectional patterns should be utilized".
- **Static duration and conversion:**
  - Duration multiplies the energy (+10 lg t).
  - When HIGE is used for taxi, LAmax is converted to SEL by a dipole function of slant range and speed.
- **Answer to the earlier gap (number of azimuths):** HNM stores no discrete azimuths. Hover directivity is a continuous eight-coefficient Fourier series in azimuth only, with no polar dependence.

**ACRP 02-44 final report (Page, Hobbs, May, Boeker; Volpe and Wyle; September 2015)** — [ACRP 02-44 Final Report, TRB](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)

- The INM rotorcraft NPD database "contains 3 directivity directions: 45o right, directly undertrack and 45o left". Propagation accounts for "geometric spreading, air absorption and finite ground impedance". Terrain enters only as an adjustment of slant range for ground altitude.
- AEDT/INM's 16 helicopter operational modes:
  - **Dynamic:** A (approach), D (departure), L (level flyover), B and C (decelerating approaches), E and F (accelerating departures), T (taxi).
  - **Static:** G (ground idle), H (flight idle), I and J (HIGE, HOGE), V and W (vertical ascent IGE/OGE), Y and Z (vertical descent IGE/OGE).
- "There is no fore-aft directivity in INM, except for static operations (hover and idle)." AEDT/INM "interpolates the NPD data linearly on phi".
- **How NPDs are made today:**
  - AAM simulates FAR-36-like events: a 6 deg / 63 kt approach, a 4 deg / 60 kt departure, and level flight at 492 ft / 129 kt and at 1000 ft and 5000 ft / 120 kt.
  - Conditions are 59 F, 70% RH, no wind and soft ground, with receivers at 4 ft (centre and 500 ft either side).
  - The output goes through LCorrect, the SAE-AIR-1845 Appendix B Simplified Adjustment Procedure, to give SEL, Lmax, EPNL and PNLTmax NPDs.
  - Spectral classes are assigned from the spectrum at PNLTmax. Speed coefficients are a least-squares second-order regression of PNLTmax on speed.
  - The original HNM-era helicopter NPDs used Volpe's 1979-85 methodology instead.
- **Static directivity in AEDT/INM:** Lmax and PNLTmax are given "at points 200 feet from the center of a modeled helicopter static operation every 15 degrees radially around the helicopter". From these come 0 deg NPDs and "a relative directivity adjustment for four different static operational modes": flight idle, ground idle, HIGE and HOGE. Without measured data, a similar helicopter's directivity is substituted.
- **Database counts (Table 3-1):**
  - INM/AEDT: 26 rotorcraft types and 125 recommended substitutions.
  - AAM: 37 types.
  - HELENA: EC135 and A109 hemispheres and EC130 "carpets".
  - AAM spheres are generally 5 deg fore/aft by 10 deg lateral.
- **Recommendations:**
  - Include third-octave bands from 10 Hz, because of low-frequency rotor noise and C-weighted metrics.
  - Model lateral directivity with enough angular fidelity for directional BVI, and enough lateral extent for roll.
- **Lateral decimation study, Bell 412, AAM:**
  - Decimating the sphere's lateral resolution to 30 deg or 45 deg spacing changed Lmax by up to 9.4 dB and SEL(A) by up to 2.9 dB at the 6,300 ft starboard point at 128 kt.
  - At 78 kt the changes were up to 5.3 dB Lmax and 2.7 dB SEL(A) at 6,300 ft port.
  - Under the track the change was 0.0 dB.
- **Appendix B, omnidirectional source, level flight at 120 kt:**
  - AAM and HELENA SEL versus lateral distance "agree fairly well".
  - INM departs from both with lateral distance at 492 ft altitude, attributed to ground-effect modelling.
  - With the real MD-902 sphere, INM's three-direction directivity differs from AAM/HELENA by about +/-1 dBA inside 45 deg, and by "3-5 dBA SEL" outside 45 deg towards the rotor plane.
- **Weather in segment models:** feasible only for a horizontally stratified atmosphere over flat ground, as a one-time pre-computation by source elevation, distance and bearing. Three propagation points per segment (CPA and the two ends) for segments up to 2000 ft.

**NASA's AEDT recommendations, 2025 (Letica & Rizzi)** — [NASA/TM-20250008619, Sept 2025](https://ntrs.nasa.gov/citations/20250008619) (PDF read)

- The reference is AEDT 3g (Technical Manual DOT-VNTSC-FAA-24-08, Aug 2024).
- "The 'in-ground effect' altitude should be user-specified and not hard-coded to 1.5x main rotor diameter (AEDT Technical Manual Eq. 3-74)".
- A helicopter NPD standard "along the lines of SAE AIR-1845A is needed". It should cover a "user-defined reference distance (currently set to 200 ft) for the static directivity adjustment" and the "cone angle used for hover condition measurement".
- "The extrapolation used in the AEDT lateral directivity adjustment may not accurately reflect the noise at lateral receptors with azimuthal angles > 45°."
- In the quadrotor case, the outermost fixed-wing AEDT 40 dBA contour lay at over double the lateral distance given by AMAT simulation with the full hemisphere.
- The noise fraction's "fourth-power model currently used is only marginally applicable to helicopters and AAM vehicles". To avoid an "incorrect noise fraction adjustment", all helicopter SEL NPDs in the study were corrected to a 160 kt reference speed.
- "Very few helicopters in the database provide LAmx, LPNTmx, or LEPN data for dynamic operational mode procedural steps." AEDT's maximum-level approximation "is almost certainly based on fixed-wing aircraft data".
- **Static-step workarounds:**
  - Static steps had to be suppressed for the comparisons by cutting their duration to 0.5 s and their NPD level to -140 dBA.
  - Fixed-wing mode has no hover. The noise fraction would control the contour shape, which "is inappropriate for this flight condition".
- **Recommendation for point-to-point operations:** a hybrid; helicopter type for overflights (lateral directivity), fixed-wing fixed-point profiles for departures and approaches (no limit on NPD states).

**Lateral attenuation and ground**
- In AEDT, helicopter exposure includes a lateral-attenuation term LA_ADJ alongside LD_ADJ (from Letica & Rizzi 2024, already cited in other_tools.md). The 2025 TM refers to "other AEDT adjustments (e.g., lateral attenuation adjustment)" and does not say which formula helicopter mode uses. — [NASA/TM-20250008619](https://ntrs.nasa.gov/citations/20250008619)

### Inferences
- **Hover directivity: NICE-OPS is ahead of HNM/AEDT in kind.**
  - HNM/AEDT has one LAmax curve per static mode, an azimuth-only Fourier or 15 deg directivity normalised at 200 ft, and a duration multiplier.
  - NICE-OPS has a full (phi, theta) sphere and orients hover by INS heading below 5 kt.
  - But NICE-OPS has no *static operation* as such. Hover must be flown as trajectory states, so a "30 s HIGE at the pad" needs a synthetic track. HNM/AEDT, NORAH2 and AAM all treat this as a first-class mode. This is a usability gap for heliport studies.
- **NICE-OPS hover data is synthesised.** It comes from the slowest near-level sphere, averaged fore to aft (niceops_inventory.md). HNM's HIGE directivity came from dedicated hover measurements over hard and soft ground. NICE-OPS has neither IGE/OGE distinction nor hard/soft hover data.
- **Lateral attenuation is the clearest regulatory contrast.**
  - HNM's below-45 deg term is an empirical lump of ground effect and directivity, and can be a gain.
  - NICE-OPS computes ground and directivity separately (`--ground_effect`, plates, refraction).
  - This favours NICE-OPS physically. But it means NICE-OPS will not reproduce AEDT numbers off-track even with identical source data, which matters for any regulatory cross-check.
- **NICE-OPS reports only A-weighted SEL** (niceops_inventory.md). That falls short of HNM/AEDT/ACRP practice, which carries LAmax, PNLT/EPNL and C-weighted metrics, and of ACRP's 10 Hz lower band.

### Gaps
- The AEDT 3g Technical Manual itself was not retrieved. The exact equations used for helicopter lateral attenuation (whether SAE-AIR-5662 or the HNM low-angle coefficients) and for the static directivity in AEDT are therefore unconfirmed. The 15 deg/200 ft data format is from ACRP 02-44; the 1.5 D IGE height and 200 ft reference are from NASA's TM.
- The HNM methodology was OCR'd from a scan. The equation forms (low-angle adjustment, Mach polynomial, dipole conversion) are reconstructed from garbled text and should be checked against the PDF before being quoted as formulas.
- The HNM v2.2 User's Guide ([ROSA-P dot/8828](https://rosap.ntl.bts.gov/view/dot/8828)) was not opened, so changes between v1 and v2.2 are unknown.

## 2. Empa sonAIR: does it model helicopters, and with what source and propagation models?

### Takeaway
Yes, but only in simplified form. sonAIR's peer-reviewed description lists "simplified models for a wide variety of propeller driven aircraft and helicopters", alongside its detailed turbofan models. No helicopter-specific sonAIR emission paper, parameterisation or validation was found. The published validation (about 20,000 flights) covers only 22 turbofan models.

In Swiss military noise calculations, helicopters sit in a combined "Prop/Heli" category and are grouped where no emission model exists. For helicopter Lmax at heliports, the Swiss regulatory recommendation is still FLULA2, not sonAIR.

The propagation core, sonX, is the strongest part. It is physics-based, in one-third octaves from 25 Hz to 5 kHz, with:

- Ground effect over mapped land cover.
- Vertical gradients of wind, temperature and humidity.
- Foliage.
- Building shielding and reflections.
- Forest and cliff echoes.

It reproduces grazing-incidence (lateral attenuation) cases to about +/-1 dB without an empirical lateral term.

### Cited Findings
- sonAIR "Several sound emission models are available, from detailed models for turbofan powered aircraft, which describe airframe and engine noise separately with a three dimensional sound directivity pattern, to simplified models for a wide variety of propeller driven aircraft and helicopters". Propagation "pays particular attention to air attenuation through a stratified atmosphere". The chain runs in one-third octaves from 25 Hz to 5 kHz. **[snippet only: abstract]** — [Wunderli et al., "sonAIR – a GIS-Integrated Spectral Aircraft Noise Simulation Tool…", Acta Acustica united with Acustica, 2018](https://consensus.app/papers/details/aa9083c93f1b5e0daa674c3fd47b3004/?utm_source=claude_desktop)
- **The validation is turbofan only.** About 20,000 single flights were simulated with "the 22 currently available sonAIR emission models of turbofan aircraft", at 22 positions 1.8-53 km from Zurich and Geneva. Mean differences were "well below ±1 dB", with standard deviations of ±1.7 or ±2.4 dB depending on model type. **[snippet only: abstract]** — [Jäger et al., "Validation of the sonAIR aircraft noise simulation model", Noise Mapping 2021](https://consensus.app/papers/details/9d45c3ebc42155b0a5f66d58f92b382a/?utm_source=claude_desktop)
- **Emission-model form** (turbofans; read in full):
  - Airframe is L_0(log Ma, log rho, flaps, gear, speed brakes, procedure) plus a polar radiation term.
  - Engine is L_0(Ma, N1) plus polar and azimuthal radiation terms that depend on N1.
  - The model is a regression on measured fly-bys, "each aircraft (sub)type needs its own model", "as a consequence, a limited number of aircraft types have been added".
  - Propagation is either "direct" or sonX. sonX "accounts for atmospheric absorption, ground effect, foliage attenuation, and the influence of vertical gradients of wind, temperature, and relative humidity".
  - — [Boelhouwer, "Assessing the performance of the sonAIR aircraft noise model in predicting noise levels at Schiphol Airport", TU Delft MSc thesis, May 2023](https://repository.tudelft.nl/file/File_6cf754ab-afef-4a83-b126-1ae6bb26b43f)
- **Swiss military practice (Dec 2023 presentation, read in full):**
  - The fleet mix is "Jets, Props, Helis", with "Gruppierung bei fehlenden Emissionsmodellen" (grouping where emission models are missing). Movement tables list a combined "Prop/Heli" category: 3,550 to 20,160 movements a year depending on base.
  - sonAIR's stated advances over FLULA2:
    - Emission separated from propagation.
    - 3-D directivity.
    - Explicit Doppler, hence speed dependence.
    - A detailed propagation model with ground cover and "physikalische Modellierung von Situationen mit streifendem Schalleinfall" (physical modelling of grazing incidence).
    - Calculation in third octaves.
  - Local differences against FLULA2 come from Doppler (louder behind the runway threshold), ground reflection, lateral directivity in banked flight, and terrain shielding on slopes.
  - — [Schäffer (Empa), "Lärmberechnungen des F-35A", Bern, 5-7 Dec 2023](https://www.newsd.admin.ch/newsd/message/attachments/85145.pdf)
- **Swiss regulatory status (read in full):**
  - BAFU recommends FLULA2 (Version004) and IMMPAC 3.3 for "die Berechnung der Beurteilungspegel Lmax von Helikoptern auf Helikopterflugplätzen" (helicopter Lmax at heliports).
  - A supplement dated February 2026 recommends sonAIR for military airfields (LSV Annex 8).
  - — [BAFU, "Empfohlene Lärmberechnungsprogramme für die Fluglärmberechnung in der Schweiz"](https://www.bafu.admin.ch/dam/de/sd-web/ORQNEak1MGUx/empfohlene_laermberechnungsprogrammefuerdiefluglaermberechnungin.pdf)
- **Grazing incidence:** sonAIR's physics-based propagation gives "deviations on the order of only ±1 dB at all measurement locations" in grazing-incidence situations up to 2.8 km. The best-practice models need an empirical lateral attenuation term. **[snippet only: abstract]** — [Wunderli et al., JASA 2022](https://consensus.app/papers/details/9b9419377e9e578d8b9d3ddff58f3bd9/?utm_source=claude_desktop)
- **Cost of the extra physics:**
  - Meteorological effects "reduce noise exposure close to the airport during daytime due to acoustic shadow zones, but double computational time".
  - Buildings improve local estimates but are "impractical for large-scale noise mapping".
  - **[snippet only: abstract]**
  - — [Schalcher et al., Acta Acustica 2025](https://consensus.app/papers/details/10f22b9b8c4d5f1aac18a3ee8bbdb464/?utm_source=claude_desktop)
- Empa's multicopter emission method (Wunderli, Schäffer et al., IJERPH 20(1):96, 2023) exists, but MDPI returned 403. Its form and whether sonAIR uses it are unverified. — [doi:10.3390/ijerph20010096](https://doi.org/10.3390/ijerph20010096)

### Inferences
- sonAIR cannot be cited as a helicopter hemisphere tool. Its helicopter handling is "simplified", probably at FLULA2-like fidelity, and no paper says otherwise. Its propagation core is the relevant comparison for NICE-OPS:
  - sonX has mapped ground cover, foliage, building shielding and reflection, and cliff echoes, which NICE-OPS lacks.
  - NICE-OPS has ducted, multi-arrival stratified rays, shadow-boundary diffraction and a PE reference, which no sonAIR paper found describes.
- sonAIR is regression-parameterised (Mach, N1, density, configuration) rather than table-of-hemispheres. That is a third design pattern besides NPD and hemisphere tables. It suits a fleet measured at an airport, not a rotorcraft whose emission depends on (mu, CT, alpha_TPP).

### Gaps
- No helicopter-specific sonAIR emission model (its parameters, directivity dimensionality or types) could be found; Empa's documentation is WAF-blocked. The sonAIR documentation PDF ([dokumentation-sonair-2020](https://www.empa.ch/documents/d/s509/dokumentation-sonair-2020)) is the likely source.
- No helicopter validation for sonAIR was found.
- **A search-engine pitfall to avoid:** one summary attributed "off by a maximum of 2 dB for the three helicopters (R66, AS350, Bell 407)" to sonAIR. Those numbers belong to Penn State's prediction system (Botre et al.), not sonAIR.

## 3. HELENA, CAROT/MOTUS, SELENE and the older HELISHAPE/HELINOISE work: parameterisation, interpolation, propagation and validation

### Takeaway
New facts:

- **HELENA:**
  - It is co-owned by AgustaWestland, Eurocopter, EADS, Turbomeca, Anotec, DLR, CIRA and NLR.
  - It takes sources as hemispheres or "noise carpets".
  - It propagates with SAE ARP 866A absorption, Zaporozhets & Tokarev ground reflection, and refraction by ray tracing (Tuinstra 2007).
  - It agrees "fairly well" with AAM on an omnidirectional benchmark.
  - Its databases were EC135 and A109 hemispheres and EC130 carpets.
- **The TU Delft ray model** used in Gennaretti et al.'s quasi-steady hemisphere work has layered constant-gradient arcs, ICAO absorption, Delany-Bazley ground and ray-tube spreading. It also has Arntzen's shadow-zone model, which includes "scattering due to turbulence". That is the first hemisphere footprint propagation found with a turbulence term.

Still unknown:

- HELENA's flight-condition axes and interpolation scheme. An unattributed snippet says "three closest flight conditions".
- CAROT's parameterisation and interpolation.
- SELENE's absorption standard.

### Cited Findings
- **HELENA, ACRP 02-44 description:**
  - "HELENA explicitly models spherical spreading, atmospheric attenuation, and source directivity. The noise source can be represented either in hemispheres or in noise carpets."
  - The data come from CFD or measurement. Measurements are converted to hemispheres in a post-processing step that "generally" includes averaging.
  - "Propagation modeling is done explicitly and the model includes the effects of spherical spreading, atmospheric attenuation [SAE ARP 866A…], ground reflection [Zaporozhets & Tokarev, 2002], and Sound refraction by the atmosphere by ray tracing [Tuinstra, 2007]."
  - It is co-owned by "Agusta-Westland, Eurocopter, EADS, Turbo Meca, Anotec, DLR, CIRA and NLR".
  - Its user manuals are Meliveo 2010, HELENA 2.4 and "Helena PV 2.0".
  - — [ACRP 02-44 Final Report, Sept 2015](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- HELENA's operational modes are "Landing, Takeoff, Level flight, Static Hover". Its spectra are a "1/3-octave band directivity matrix". Its database held "EC135, A109 hemispheres, EC130 Carpets", plus Clean Sky single-engine light and twin light/medium/heavy generic types. "A performance modeling capability within HELENA also has capability to compute the noise impacts" of decelerations. — [ACRP 02-44, Table 3-1 and Sec. 4.2](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- HELENA validation (Gervais et al., AHS Forum 66, 2010), from the abstract only: HELENA is "designed to provide high accuracy rotorcraft noise footprint predictions based on a database of source data which contains directivity and frequency information as a function of flight condition". It was compared with measured EC130 and EC135 single-segment procedures in integrated metrics and time histories, against an NPD-based model, and on one multi-segment manoeuvre. **[snippet only: abstract; no numbers available]** — [Gervais et al. 2010](https://consensus.app/papers/details/9ec15701cb4b5ca5ac9fd637a70797e0/?utm_source=claude_desktop)
- FRIENDCOPTER's own summary: HELENA "takes into account the complexity of noise generation specific to rotorcraft". Noise abatement procedures gave "up to 50 % noise footprint reduction, more than 10 dB reduction at specific observer locations" on three aircraft. These are procedure-benefit numbers, not model accuracy. — [CORDIS, FRIENDCOPTER (502773) reporting](https://cordis.europa.eu/project/id/502773/reporting)
- **CAROT / MOTUS:**
  - The H130 database (FRIENDCOPTER flight tests) is shown as mean A-weighted hemisphere level against steady flight condition, on axes that read as airspeed and vertical speed.
  - The BVI rise is at negative vertical speed. Level flight rises with speed, through tip Mach number.
  - CAROT runs over flat homogeneous ground or "realistic scenarii" with topography, demography and traffic.
  - The paper states no interpolation scheme, absorption or ground model.
  - — [Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1) (re-read)
- CAROT's flight-procedure side is Dieumegard's CIFRE thesis at ENAC with Airbus Helicopters, which won the 3AF thesis prize. **[snippet only]** — [ENAC alumni news](https://alumni.enac.fr/en/news/prix-de-la-these-3af-l-enac-a-l-honneur-944)
- **TU Delft ray model in a hemisphere footprint study (read in full):**
  - The atmosphere is "divided in a number of layers, and in each of the layers a constant speed-of-sound gradient is assumed", with circular arcs per layer.
  - The emission direction is where the launch ray crosses the hemisphere.
  - Spreading is corrected for focusing by the spacing of neighbouring rays.
  - Absorption uses the ICAO method.
  - Ground uses a plane-wave reflection coefficient with Delany-Bazley impedance.
  - "Based on the approach developed by Arntzen et al., the strong decrease of the sound level at the transition between the illuminated and the shadow zone is determined, as well as the noise levels penetrating the shadow zones due to ground waves, diffraction and scattering due to turbulence."
  - Hemisphere indexing compared (mu, gamma) with (mu, CT, alpha_TPP).
  - — [Gennaretti, Bernardini, Serafini, Anobile, Hartjes, "Helicopter noise footprint prediction in unsteady maneuvers", Int. J. Aeroacoustics 16(3) 2017](https://repository.tudelft.nl/file/File_8a12a003-8adb-4fca-b37d-03c081592697)
- Gennaretti et al. (AST 2016): "at least three parameters, including disk loading, are required to adequately characterize noise hemispheres". "Similarity of kinematic parameters alone may yield steady flight acoustic predictions poorly correlated with those obtained for unsteady maneuvers." **[snippet only: abstract]** — [Gennaretti et al., "Numerical characterization of helicopter noise hemispheres", Aerosp. Sci. Technol. 2016](https://consensus.app/papers/details/6d70fa64cca3505e90368ff77076782d/?utm_source=claude_desktop)
- DLR's pre-SELENE footprint model (Yin 2004): the "HEMISPHERE" propagation code includes "Doppler effect correction, spherical spreading, atmospheric absorption and ground reflection". Its source is UPM-Mantic/APSIM, or measured arrays reverse-propagated. The absorption standard is not stated. **[snippet only: catalogue record; not open access]** — [Yin, "Towards Noise Abatement Flight Procedure Design: DLR Rotorcraft Noise Ground Footprints Model and its Validation", 2004, elib.dlr.de/12976](https://elib.dlr.de/12976/)
- **HELISHAPE and HELINOISE** were Brite-EuRam rotor aerodynamics and aeroacoustics programmes, producing, for example, the HERNOP-2 time-domain FW-H code. They are source-prediction projects, not footprint tools. **[snippet only]** — [DLR elib 36440](https://elib.dlr.de/36440/); [CNR IRIS record](https://iris.cnr.it/handle/20.500.14243/210502)

### Inferences
- HELENA's propagation (ARP 866A absorption, impedance ground, ray-traced refraction) is broadly what NICE-OPS's straight and linear modes do. NICE-OPS uses ISO 9613-1 and adds ducting, shadow-boundary diffraction and PE checks.
- **The shadow-zone level is where NICE-OPS should check itself against the literature.** The TU Delft/Arntzen model puts ground-wave, diffraction and turbulence-scattered energy into shadows. NICE-OPS's `--scattering` (a fit to PE phase screens) is the analogous term, and is off by default.
- Gennaretti 2016 supports NICE-OPS's inclusion of CT, a disk-loading parameter. Indexing by (mu, gamma) alone, which RNM/AAM's (V, gamma) amounts to, is the weaker choice for manoeuvres.

### Gaps
- **HELENA's interpolation axes and scheme.** A search summary said HELENA "interpolates the corresponding noise hemisphere spectra from the three closest flight conditions", then applies Doppler, absorption and ground reflection. The originating document could not be identified, so this is **unverified**. The Gervais 2010 paper and the Meliveo user manuals were not obtainable.
- No HELENA accuracy statistics (bias or rms) were found.
- CAROT's flight-condition axes, interpolation, absorption and ground models remain undocumented in anything retrievable. Dieumegard's thesis is the likely source but was not located online.
- SELENE's absorption standard and ground-impedance model are still unconfirmed.
- No HELISHAPE/HELINOISE footprint tool was found. Those programmes appear to have been source-prediction only.

## 4. RNM 7 and APET, ANOPP2's AMAT and AARON, and whether "ARNM" exists

### Takeaway
RNM's core method is now documented in its original paper:

- **Within a hemisphere:** a four-way Coons-surface interpolation in (theta, phi).
- **Between hemispheres:** nearest flight path angle and pylon angle, then linear interpolation in pressure squared on airspeed, which can be switched to nearest.
- **Time step:** half a second.
- **Propagation:** straight rays with an altitude-averaged absorption coefficient, and a complex-impedance two-ray ground.

The 2016 tech brief says RNM 7 adds propagation "through a horizontally stratified atmosphere over uniform terrain with winds", but nothing retrievable describes how. APET adds tables from a ray tracer and GFPE.

AMAT is confirmed as a command-line wrapper of ANOPP2's straight-ray propagation module, with:

- Doppler, spreading and absorption: lossless, SAE ARP 866A, ANSI S1.26-1978, ANSI S1.26-2014 or ICAO 1978.
- Hard, soft (flow resistance) or no ground.
- An optional layered atmosphere.

AARON couples CAMRAD II with F1A tonal noise and an ANOPP2 self-noise module for broadband. No ANOPP2 tool called "ARNM" exists in any source found.

### Cited Findings
- **RNM, original paper:**
  - Sources are broadband third-octave hemispheres and narrowband pure-tone hemispheres with phase, up to ten of each.
  - The propagation loop advances the aircraft "at half second intervals".
  - Euler angles come from linear interpolation between track points.
  - Within a hemisphere, "a four-way interpolation for theta and phi is made using a 'Coons' surface interpolation technique".
  - Between hemispheres, the search finds those "with flight path and pylon angles that are nearest to the specified conditions", and then "two sound hemispheres which will be used to interpolate sound levels based on airspeed". That interpolation "assumes a linear variation in pressure squared with airspeed", and a keyword switches it to nearest.
  - "The propagation model this program uses assumes that the ray paths are straight lines, there is no wind present". Temperature gradients enter only through weighted-average sound speed and absorption tabulated every 1,000 ft.
  - Ground is "a finite impedance … provided by the user" (direct plus reflected ray).
  - Pure tones are summed coherently, then added to the bands in pressure squared.
  - — [Lucas & Marcolini, "Rotorcraft Noise Model", Wyle/NASA Langley, c. 1997 (NTRS 20040110391)](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)
- RNM 7: "RNM calculates the effects of sound propagation over varying ground terrain for acoustic sources using geometrical theory of diffraction algorithms, and through a horizontally stratified atmosphere over uniform terrain with winds". Its developers were Page, Wilmer and Plotkin of Wyle for NASA Langley. **[snippet only]** — [NASA Tech Briefs LAR-17753-1](https://www.techbriefs.com/component/content/article/24556-lar-17753-1); [NASA software catalog LAR-17753-1](https://software.nasa.gov/software/LAR-17753-1)
- APET is "a modified version of the Rotorcraft Noise Model (RNM, Version 7) coupled with the newly developed Spectral Attenuation Method (SAM) that contains the Ray Tracing Program (RTP) and two versions of a Greens Function Parabolic Equation (GFPE) code". **[snippet only]** — [NASA software catalog APET LAR-17761-1](https://software.nasa.gov/software/LAR-17761-1)
- ACRP 02-44 on AAM, RNM's Volpe successor:
  - AAM computes "at one second intervals" typically.
  - It "can optionally account for varying ground terrain or atmospheric gradient effects", and computes nonlinear propagation for high-level sources.
  - "The AAM algorithm interpolates linearly on (phi, theta) for each one-third octave band before propagating."
  - "Generally noise spheres in AAM are defined using 5-degree fore/aft and 10-degree lateral spacing."
  - — [ACRP 02-44, Sec. 3.2 and 4.1](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **AMAT (ANOPP2 Mission Analysis Tool), read in full:**
  - It is a "Command line tool that wraps the ANOPP2 Straight Ray Propagation Internal Functional Module (ASRPIFM)".
  - It accounts for "Doppler shift, spherical spreading, atmospheric absorption, ground effects".
  - Its absorption options are "lossless … SAE ARP 866a, ANSI S1 26 1978, ANSI S1 26 2014, ICAO6 1978".
  - Its ground options are "Hard/Soft/No Ground, Flow Resistance".
  - "Layered atmosphere (not required for NPD generation)."
  - Sources may be one-third-octave SPL, narrowband or PSD.
  - It runs in level flyover, approach, takeoff and manoeuvre modes, and generates AEDT NPDs.
  - — [Lopes, "ANOPP2 Mission Analysis Tool (AMAT)", UNWG SG1, 7 Dec 2022, NTRS 20220018599](https://ntrs.nasa.gov/citations/20220018599)
- **AARON and the ANOPP2 tool chain:**
  - AARON is "Comprehensive analysis using CAMRAD II for vehicle trim and lifting-line motion and forces". It uses the Farassat 1A module (AF1AIFM) for tonal noise and the "ANOPP2 Self-Noise Internal Functional Module (ASNIFM)" for broadband, with optional CFD-coupled aerodynamics.
  - The other tools in the chain:
    - ARMIT attaches metadata to hemispheres.
    - ASOFT filters a hemisphere to selected azimuths, for example 0 deg for "fixed-wing mode" or -45/0/45 deg for "helicopter mode", for V&V.
    - AMPIT is MPI parallelism.
    - ASOCIT converts between AAM (netCDF), ANOPP2 restart and ANOPP ACD formats.
  - — [Lopes 2022](https://ntrs.nasa.gov/citations/20220018599)
- AMAT "is a time-marching simulation program that uses the full source noise hemisphere, akin to the Volpe Advanced Acoustic Model (AAM)". NASA states it is not "the 'truth'". — [Letica & Rizzi, NASA/TM-20250008619](https://ntrs.nasa.gov/citations/20250008619)

### Inferences
- **On interpolation, NICE-OPS differs from RNM in three ways:**
  - Its across-condition blend is barycentric in dB over a 3-D Delaunay mesh. RNM's is nearest-neighbour in two axes, then linear in pressure squared in airspeed only.
  - Its within-sphere interpolation is bilinear in dB. RNM's is a Coons patch.
  - Its time step is the trajectory's own 0.02 s, against RNM's 0.5 s and AAM's typical 1 s.
- Neither interpolation domain (energy or dB) has been shown to be better for footprints. NLR's argument for dB applies to fitting. This is a sensible item for NICE-OPS to test, for example on the self-test residual.
- The ASOFT "azimuth filter" V&V pattern would be easy to copy in NICE-OPS: reduce a sphere to the AEDT 0/+/-45 deg slices and fly it. It would give a direct, apples-to-apples comparison against AEDT-style NPD predictions on the day-284/286 data.
- **Absorption standards.** AMAT offers five and RNM uses ANSI S1.26; NICE-OPS uses ISO 9613-1 only, which is numerically close to the ANSI S1.26 family. Adding SAE ARP 866A would matter only for reproducing legacy NPD and certification corrections, not for physics.

### Gaps
- No document was found describing RNM 7's stratified-atmosphere and wind implementation (ray tracing or an effective-sound-speed correction), or validating it.
- The APET tables' parameterisation (source height, range, bearing?) and validation are not documented in anything retrievable.
- AMAT's ground model (which impedance model, with flow resistance as input) and whether it applies convective amplification as well as frequency Doppler are not stated in the slides.
- No source mentions an ANOPP2 "ARNM"; it should be treated as non-existent unless the brief's author has a reference.

## 5. Newer (2023-2026) eVTOL/UAM hemisphere footprint tools and the NORAH2 extension

### Takeaway
2023-2026 work keeps the hemisphere pattern and pushes it to eVTOL. The groups and tools found:

- **NASA:** AARON/AMAT feeding AEDT NPDs.
- **MIT, ASCENT Project 084:** an ANOPP2 hemisphere library for tilt-VTOL, lift-plus-cruise and STOL vehicles. It reports 1-3 dB agreement with industry data in steady flight, and that 13x coarser hemispheres move ground levels by only 1-3 dB on average. It flags transition as a modelling gap AEDT cannot represent.
- **TUM:** a CFD/FW-H hemisphere workflow, explicitly without ground reflection or absorption.
- **Korean, Chinese and ROM-based variants.**
- **NORAH2 (released January 2025):**
  - Adds three helicopters, taking fleet coverage above 80%.
  - Adds building and topography screening.
  - Adds guidance for taxi, idle, hover and turns.
  - Its sources say the approach is "suitable" for drones and UAM, but none found says eVTOL hemispheres are in its database.

The recurring eVTOL complication is a source that depends on more state than (V, gamma): rotor tilt, powerplant state, crossflow and interactions. NICE-OPS's fixed 3-D (mu, gamma, CT) axes share that limit.

### Cited Findings
- **NORAH2 (NLR news, 13 June 2024):**
  - "Three medium- to heavy weight helicopters were added to the noise database", giving types "representative of the bulk of the helicopter fleet (more than 80%) operating in Europe".
  - "The propagation modelling method was extended with screening effects from buildings and topography."
  - "Specific modelling guidance was introduced to enable inclusion of taxi, idling, hover and turns operations."
  - "The hemisphere approach is suitable for modelling the noise of helicopters but also for future transport modes such as drones and urban air mobility (UAM)."
  - The consortium was NLR with Anotec Engineering and SINTEF.
  - — [NLR, "The NORAH2 project helps predict noise nuisance of rotorcraft"](https://www.nlr.org/newsroom/nieuws/the-norah2-project-helps-predict-noise-nuisance-of-rotorcraft/)
- A search summary named the three added types as Cabri G2, A109 and S-92 (70% to 80%). **[snippet only; not seen in a document I read]**
- EASA released the NORAH2 software prototype and hemisphere database for free download in January 2025, as a Windows prototype needing no installation. **[snippet only]** — [AIN, 7 Jan 2025](https://ainonline.com/aviation-news/aerospace/2025-01-07/easa-releases-rotorcraft-noise-software-and-database)
- **ASCENT Project 084 (MIT, FAA; 2025 annual report, Apr 2026):**
  - "Modeling noise using a simple uniform inflow model results in reasonably close results compared to the validation data (within 1-3 dB) for steady state conditions for the tilt rotor, lift plus cruise, and STOL vehicle."
  - "The noise output of the tilt-VTOL vehicle is found to be particularly sensitive to the powerplant condition during transition. Given that standard noise modeling tools like the AEDT are suited for a helicopter or traditional airplane models, but not something in between, this represents a notable gap."
  - "Reducing the resolution by approximately 13 times resulted in differences on the ground of only 1 – 3 dB on average", with roughly 40-fold faster hemisphere generation.
  - The maximum-noise contour was sensitive to "increased, un-modeled crossflow when flying with rotors at nonstandard angles".
  - About 15 dB more noise over 45 deg of rotor tilt; about 8 dB average increase from rotor-wing interaction.
  - — [ASCENT Project 084 2025 Annual Report](https://s3.wp.wsu.edu/uploads/sites/2479/2026/04/ASCENT-Project-084-2025-Annual-Report.pdf)
- **NASA RVLT's final AEDT recommendations** for UAM: user-defined operational modes (for example decelerating climb), no restrictions on step sequence, and lateral directivity beyond left-centre-right. See section 1. — [NASA/TM-20250008619, Sept 2025](https://ntrs.nasa.gov/citations/20250008619)
- **TUM UAM footprint workflow (DAGA 2025):**
  - Tilt-rotor cruise hemispheres come from CFD and FW-H, tonal only, at 100 m. Hemispheres are placed at waypoints, with "linear interpolation across the noise hemisphere" along the emission vector and Doppler applied before interpolation.
  - "Background noise, ground reflection, and atmospheric absorption were not considered in this analysis."
  - A time-step-to-grid "CR" ratio is used to pick trajectory resolution.
  - — [Moses, Schmähl, Hornung, "Evaluating the Ground Noise Footprint of Urban Air Mobility Vehicles", DAS|DAGA 2025](https://pub.dega-akustik.de/DAS-DAGA_2025/files/upload/paper/299.pdf)
- **Alternatives to hemisphere tables:**
  - An "acoustic modal database" for AS350 footprints runs in one-fifth of RNM's time for 1800 observers. — [Shao et al., JASA 2024](https://consensus.app/papers/details/334b6cd535385d9f96e72d5844829c9b/?utm_source=claude_desktop)
  - A URANS-trained reduced-order model predicts AH-1G hemispheres to under 0.47 dB mean OASPL error. — [Kwon et al., AIAA Aviation 2023](https://consensus.app/papers/details/9b6d63a7de7c5b618fa452e0a1046c42/?utm_source=claude_desktop)
  - **[snippet only: abstracts]**
- The Korean vertiport studies (Lee & Haan; Jo & Wie, Forum Acusticum 2025) are already covered in other_tools.md and are not repeated here.

### Inferences
- **What NICE-OPS shares with eVTOL tools, and what it lacks.** Like every eVTOL tool found, NICE-OPS fixes the sphere index at three operating parameters.
  - Tilt-rotor transition needs at least one more axis (nacelle or tilt angle, as RNM has had since the 1990s). So do powerplant state and crossflow.
  - NICE-OPS's `Delaunay` is dimension-agnostic, per CLAUDE.md, which makes a fourth axis structurally cheap. The `.nod` schema and panam would need the extra coordinate.
- ASCENT's finding that 13x coarser hemispheres move the ground by 1-3 dB on average is a useful bound for NICE-OPS's own grid-resolution sensitivity. It contrasts with ACRP's 3-9 dB local errors when only the *lateral* resolution is cut to 30-45 deg. So averages hide large local errors near the rotor plane, which is where NICE-OPS's grazing-pair work lives.
- NORAH2's added screening (buildings and topography) and its static and turn guidance are exactly what NICE-OPS lacks: no terrain or obstacles, and no static operations. NICE-OPS's atmospheric physics is well beyond NORAH2's.

### Gaps
- No source confirms that NORAH2's released database contains any eVTOL or drone hemisphere, or that a NORAH2 eVTOL extension project exists. Only the "suitable for" statement was found.
- The MDPI *Drones* 2025 survey ("Noise Prediction and Mitigation for UAS and eVTOL Aircraft") could not be opened (403).
- No 2023-2026 tool was found that indexes hemispheres by more than three parameters, or that interpolates eVTOL hemispheres across transition. ASCENT 084 names this as an open gap.
- EASA's eVTOL noise guidance (beyond the NORAH2 material in other_tools.md and norah2.md) was not located in this pass.
