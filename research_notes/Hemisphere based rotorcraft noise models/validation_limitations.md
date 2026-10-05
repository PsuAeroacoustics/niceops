# Validation and documented shortcomings of hemisphere-based (source noise sphere) rotorcraft and eVTOL noise prediction

Scope note: this covers published comparisons of sphere/hemisphere-based tools (RNM, AAM, NORAH/NORAH2, HELENA, FRAME/FRAME-QS/Q-SAM, NICE-OPS) against flight test or certification data, plus the measured variability that bounds what any such comparison can show. Where only an abstract or a search snippet was readable this is said explicitly. Several primary documents (Watts et al. 2019 Noise Abatement data report, the HELENA EC130/EC135 validation paper, the Joby CAA flyover paper, Botre's full JAHS paper) could not be read in full; numbers attributed to them come from abstracts.

## 1. What accuracy do hemisphere-based tools achieve against flight test, and which conditions are worst?

### Takeaway
Published, quantified sphere-vs-flight-test errors are scarce and mostly small-sample: roughly 2 dBA SEL "with few exceptions" (NICE-OPS, one S-76D decelerating approach), 2-4 dB SEL (a PSU physics-based system, five helicopters), within 2 dB OASPL transients (FRAME-QS, two Bell 430 maneuvers), and 1.6-1.8 EPNdB for NORAH2 against S-92 certification. No published study found reports a held-out, multi-run SEL error distribution (mean/std over many observers and approaches). Worst conditions are consistently BVI-dominated descents and transient maneuvers, large lateral distance / near-rotor-plane emission angles, near-horizon fore/aft directions, and hover.

### Cited Findings
- **NICE-OPS (approach, S-76D):** simulated vs measured ground contours for an S-76D −7.5° decelerating approach from 60 KIAS showed "good agreement ... within 2 dBA with few exceptions" and "higher uncertainty at large lateral distances." This is a slide-level claim from one condition, not a statistical validation. — [Pascioni, Greenwood, Watts, Stephenson, Smith, "Medium-Sized Helicopter Noise Abatement Flight Test," VFS 76th Forum, Oct 2020 (NTRS slides)](https://ntrs.nasa.gov/api/citations/20205008035/downloads/Pascioni_VFS2020.pdf)
- The same work states the reason a sphere database alone is insufficient for approaches: the measured array covered "only a small portion of the area impacted," and "even small longitudinal accelerations can have a substantial influence on the noise generation processes during approach." NICE-OPS needed heuristics (operating-envelope closure; loading-noise scaling ΔSPL = 20 log10(CT/CW)) because decelerating states "can't measure in steady flight." — [Pascioni et al. 2020 (abstract)](https://consensus.app/papers/details/ffeba51f6be759c9ba58e81dc73c0ccc/?utm_source=claude_desktop); [slides](https://ntrs.nasa.gov/api/citations/20205008035/downloads/Pascioni_VFS2020.pdf)
- **FRAME-QS (Q-SAM applied to FRAME spheres), Bell 430:** for an aggressive pitch-up (~10 deg/s, peak 1.5 g, from a 6° descent at 80 kt, ~330 ft) the variable-CT model captured "the transient increase in noise within 2 dB in all three directions" (OASPL, observers ~10° below the horizon, ahead and 30° either side). The constant-thrust (CT = CW) model missed the increase; the difference is close to the 20 log10(nz) ≈ 3.5 dB heuristic. The "conventional" model (no sphere reorientation, no effective-FPA shift) showed no increase at all. Validation was "in low wind conditions" for two maneuvers only. — [Greenwood, Rau, May, Hobbs, "A Maneuvering Flight Noise Model for Helicopter Mission Planning," AHS 71st Forum 2015 (NTRS 20160005979); later JAHS 65(2), 2020](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005979.pdf)
- **FRAME-D (semi-empirical BVI), Bell 206B maneuvers:** for a fast cyclic pull-up, peak amplitudes and pulse shapes were "well captured at three microphone locations." The slower cyclic maneuver was captured "not as well ... due to the variability of the measured BVI noise." For a roll transient, "the predicted levels are too high towards the advancing side, suggesting that more improvements of the modeling of turning flight are probably required." Qualitative (time histories), no dB statistics. — [Greenwood, Schmitz, Sickenberger, "A Semi-Empirical Noise Modeling Method for Helicopter Maneuvering Flight Operations" (NTRS 20120008939; JAHS 60(2), 2015)](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf)
- **PSU physics-based prediction system (not sphere-based, but validated on the same class of data):** SELdBA, OASPL and LA "matched the trends and levels within 2–4 dB of measured flight-test data" for five commercial helicopters in level flight, descents, level and decelerating turns, and descending turns. It "missed blade-vortex-interaction (BVI) noise directivity during some complex maneuvers." The 6° steady descent was louder than any complex procedure. — [Botre et al., JAHS 2026 (abstract)](https://consensus.app/papers/details/d322dbd85eeb57d6a3aadb168aa1b9bc/?utm_source=claude_desktop); [Botre, PhD dissertation, Penn State 2020](https://etda.libraries.psu.edu/catalog/18150mcb355); [Botre et al., VFS 75th Forum 2019 (abstract: "within a few SELdBA")](https://consensus.app/papers/details/b473d549e7bf520090851342146d8bb1/?utm_source=claude_desktop)
- **NORAH2 vs certification (S-92):** NORAH2 predicted 96.0 / 99.0 / 97.4 EPNdB for take-off / flyover / approach against certified 94.6 / 97.2 / 97.5, i.e. +1.6 / +1.8 / −0.1 EPNdB. The hemisphere measurements themselves reproduced certification to within −0.6 to +1.6 dB across Cabri G2, A109 and S-92 (e.g. A109 approach 91.7 vs 90.1 EPNdB). — [Olsen, Tuinstra, van Oosten (NLR/SINTEF/Anotec), "D3 NORAH 2.0 validation," EASA Research Project NOISE-SC03, 31 Jan 2024](https://www.easa.europa.eu/en/downloads/139664/en)
- **NORAH2 has not been validated against measured in-service operations.** The deliverable states that "comparisons with independent measurements of real flight operations could not be made" because of budget priorities. Second-by-second recordings of S-92 normal operations (July 2022) were collected for that purpose but "the final analyses ... has not been made." It lists in-depth tests of hemisphere interpolation, turns and hover/taxi against measurements as future work. — [EASA D3 NORAH 2.0 validation, 2024](https://www.easa.europa.eu/en/downloads/139664/en)
- **NORAH2 vs Nord2000 (model-to-model, propagation only):** over an artificial 2×2 km flat-plus-mountain case (R22 take-off with a 60° turn), differences were "generally in the order of ±1 dB." They were up to ±10 dB where terrain screening exceeded ~20 dB, and Nord2000 was 2–3 dB higher behind the helicopter over flat soft ground, which the authors attribute to Nord2000's turbulence reducing the ground effect. Their conclusion is that NORAH2 "may slightly underestimate the noise for long distance ground-to-ground propagation." — [EASA D3, 2024](https://www.easa.europa.eu/en/downloads/139664/en)
- **RNM (Eglin 2003/2005, MD600N):** RNM was run seven times, once with each of seven "housekeeping" (110 kt, 150 ft) source semi-spheres, and the predicted SEL varied by ≈ ±0.5 dBA at all sideline distances to ±3000 ft. That figure is the floor from source/measurement repeatability, not a prediction error. — [Conner, Burley, Smith, "Flight Acoustic Testing and Data Acquisition for the Rotor Noise Model (RNM)," AHS 62nd Forum 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
- **RNM, XV-15:** measured and RNM-predicted footprints for a 3°-to-9° segmented approach, plus a centerline LA time-history comparison, were shown as figures, without tabulated dB error. — [Conner & Page, "A Tool for Low Noise Procedures Design and Community Noise Impact Assessment: The Rotorcraft Noise Model (RNM)," AHS 2002](https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-2002-ahsitsm-dac.pdf)
- **Worst condition: near the rotor plane / large lateral distance.** Model-to-model, INM (no lateral directivity) vs AAM/HELENA differed by about ±1 dBA SEL inside the 45° slant angle, but "for distances outside of the 45° slant angle (up closer to the rotor plane) the differences can be significant, for the MD-902 on the order of 3-5 dBA SEL," worsening with bank angle or canyons. — [Page et al., ACRP 02-44 "Guidance for Helicopter Community Noise Prediction," Final Report, Sept 2015, App. B](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **Worst condition: near-horizon fore/aft (eVTOL).** Computational predictions of the Joby aircraft at 60 and 100 kt exceeded the NASA-measured hemisphere over its whole measured portion, "with the largest differences directly forward and aft of the aircraft near the horizon, with up to a 12 dB(A) discrepancy." Source: search snippet of the abstract only; the full paper was not accessible. This is a physics-model vs hemisphere comparison, not a hemisphere-model validation. — ["Flyover Noise Computations of the Joby Aviation Aircraft" (ResearchGate, 2023)](https://www.researchgate.net/publication/371996452_Flyover_Noise_Computations_of_the_Joby_Aviation_Aircraft)
- **Worst condition: maneuvers (the quasi-static-blind case).** Measured peak LA,max differences between paired turns: Bell 206 descending vs level left 35° turn +8.2 dBA; Bell 407 descending vs level turns +7.5 dBA; Bell 205 accelerating vs decelerating right turns 13.3 dBA; Bell 205 decelerating vs accelerating left turn +5 dBA; S-76D 0.9 dBA. These are the magnitudes a steady-sphere model gets wrong if it does not shift the effective state. — [Stephenson, Watts, Greenwood, Pascioni, "Development and Validation of Generic Maneuvering Flight Noise Abatement Guidance for Helicopters," VFS 76th Forum 2020 (JAHS 67(1), 2022)](https://ntrs.nasa.gov/api/citations/20205003220/downloads/VFS76_Paper_v7.pdf)
- **Worst condition: hover / near-vertiport (eVTOL).** For Joby hover, levels at a fixed observer varied 2–5 dB within a single run. Back-to-back runs differed by nearly 5 dB at θ = 155° (HIGE), which the authors call "unexpected." Repeated approach/departure states were "often too dissimilar to truly compare as repeated conditions." — [Pascioni, Watts, Houston, Lind, Stephenson, Bain, "Acoustic Flight Test of the Joby Aviation Advanced Air Mobility Prototype Vehicle," AIAA/CEAS Aeroacoustics 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
- **Unpublished internal NICE-OPS harness results on the 2017 B407 day-286 approaches** (held-out, ~25 approaches; for comparison only, not literature):
  - Forward-line SEL over-prediction is about +1.5 dB median, concentrated at the low, slow end of the approach (+2 dB below 100 ft, +3.2 dB at 0–20 kt). With energy means, the slow end is +1.5 to +3.9 dB.
  - At long range (>4000 ft, 2–10° elevation) the error depends on the database rim: −0.7 dB (400–660 ft sources) and −1.8 dB (660–1000 ft) with the shipped 2° rim, against +0.25 and −0.3 dB with a 10° rim. With energy means it is about −1.5 dB.
  - [local memory note: forward-line-error-source.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/forward-line-error-source.md)

### Inferences
- The literature's "within 2 dB" / "2–4 dB" figures are single-case or few-case claims, mostly read off contour or time-history plots. No published sphere-tool validation found gives a distribution of SEL errors over many held-out approaches and observers. A NICE-OPS validation on the 2017 held-out set (approaches, turns, accelerating turns) with per-observer statistics would appear to be the most complete such comparison published, if written up.
- The consistent pattern across tools is that error grows toward the rotor plane (large lateral or long range, near-horizon) and in hover or low-speed deceleration. These are the regions where sphere coverage is thinnest and state sensitivity highest.
- NORAH2, the European regulatory tool, has formally been validated only against certification points and another model, not against in-service footprints. Certification EPNL agreement (±1.8 dB) says little about footprint accuracy at low elevation.

### Gaps
- Watts et al. 2019, "Noise Abatement Flight Test Data Report" (NASA TM; [ResearchGate listing](https://www.researchgate.net/publication/333320599_Noise_Abatement_Flight_Test_Data_Report)), was not readable here, so its exact dates, sites and run counts for the 2017 campaign are unconfirmed from the literature.
- No quantified predicted-vs-measured comparison for the 2017 light-helicopter campaign (B407, AS350, EC130, R44, R66, B206) was found in the open literature. The Fly Neighborly AAM sphere products ([Page et al., Inter-Noise 2021](https://consensus.app/papers/details/3ce8695f4b3c530586c69e924047f80f/?utm_source=claude_desktop)) and the AAM SEL scaling study ([Houston et al., JASA 2024 abstract](https://consensus.app/papers/details/9a01710b79fc5d2898822cbd46231053/?utm_source=claude_desktop)) use the spheres but do not report validation errors in their abstracts.
- HELENA validation against EC130/EC135 measurements exists ([ResearchGate listing](https://www.researchgate.net/publication/289803872_Validation_of_EC130_and_EC135_environmental_impact_assessment_using_HELENA)), but its dB results could not be retrieved.
- The Eglin 2007/2010 campaigns, the FAA/Volpe HNM validations, FRIENDCOPTER/HELINOVI and DLR comparisons were not found with quantified sphere-model errors in accessible sources.

## 2. Attributed error sources

### Takeaway
The documented, quantified error sources are, roughly in order of measured size:
1. Operating-state mismatch: the quasi-static assumption, acceleration, load factor and turns, at 5–13 dBA in LA,max.
2. Run-to-run source variability: 1–7.5 dB depending on condition.
3. Ambient-condition and wind effects on the source: up to 5.4 dB from altitude, up to 9–11 dB on BVI from wind.
4. Orientation and tracking used when building the sphere: ±1 dBA for a benign rotor, larger for directional sources.
5. Sphere coverage and extrapolation near the rotor plane: 3–5 dBA SEL model-to-model at low elevation.
6. Propagation over soft ground and long range: 2–3 dB, up to ±10 dB under screening.

Hemisphere radius, de-propagation ground-reflection removal and flush-microphone corrections are named as concerns but are rarely quantified in the validation literature.

### Cited Findings
- **Quasi-static / operating-state mismatch.**
  - Steady-sphere models that do not re-select the sphere for acceleration and load factor miss transient increases: about 3.5 dB from load factor alone in a 1.5 g pull-up. — [Greenwood et al. 2015/2020](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005979.pdf)
  - Q-SAM was originally limited to "moderate accelerations and decelerations below 0.1 g." Adding CT as a database parameter extended it to turns. — [Sickenberger, PhD dissertation, Univ. of Maryland, 2013](https://api.drum.lib.umd.edu/server/api/core/bitstreams/c9b00512-faff-4491-9c02-4a162fb6b6b4/content)
  - Numerical work concludes "at least three parameters, including disk loading, are required to adequately characterize noise hemispheres," and that kinematic similarity alone "may yield steady flight acoustic predictions poorly correlated with those obtained for unsteady maneuvers." — [Gennaretti et al., "Numerical characterization of helicopter noise hemispheres," Aerosp. Sci. Technol. 2016 (abstract)](https://consensus.app/papers/details/6d70fa64cca3505e90368ff77076782d/?utm_source=claude_desktop)
- **No interpolation between flight conditions (NORAH).** NORAH's method says: "No interpolation between flight conditions is allowed. The recommended practice is to reconstruct flight paths using existing flight conditions in the database." The vertical profile is "dictated by the available hemispheres." — [van Oosten, Meliveo, Tuinstra, Olsen, "The new EU helicopter noise model NORAH," Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
- RNM originally used "the closest matching flight condition available hemisphere, without interpolation" for some parameters. — [Conner & Page 2002](https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-2002-ahsitsm-dac.pdf)
- **Sphere coverage gaps and extrapolation.**
  - NORAH requires only azimuth ±60° and the polar 10-dB-down interval to be measured. "For directivity angles beyond the measured intervals, constant value extrapolation should be applied from the nearest filled data bin." Bins are 10°, the reference radius is 60 m and bands run 10 Hz–10 kHz. — [van Oosten et al., Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
  - In RNM/ART semi-spheres, the perimeter regions with "no acoustic data available" are visible. — [Conner, Burley, Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
  - "While it is relatively straight forward to interpolate between source hemispheres based on the flight condition, extrapolation becomes a challenging task," which motivated FRAME. — [Sickenberger 2013](https://api.drum.lib.umd.edu/server/api/core/bitstreams/c9b00512-faff-4491-9c02-4a162fb6b6b4/content)
  - Joby array geometry meant "it was not possible to acquire data for emission angles that correspond to large lateral distances." — [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
- **Lateral directivity and bank angle.** ACRP 02-44 ranks the primary weakness of AEDT/INM rotorcraft modeling as source modeling: "spectral content, lateral directivity and operational sensitivity." Its recommendation 2 calls for lateral source characteristics "with sufficient a) angular fidelity to capture directional BVI noise and b) lateral extent to account for changes in vehicle roll angle." It notes that a 30° lateral angle captured the BVI "hot spot" better than 45°. — [ACRP 02-44, 2015](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **Vehicle and rotor orientation in sphere construction.** Building the MD600N 6° approach sphere from measured rather than estimated attitude changed the RNM footprint by "always between ±1 dBA" (a benign six-bladed rotor with suppressed BVI). For the XV-15, with intense BVI, the sensitivity "can be significant." The authors recommend direct measurement of tip-path-plane orientation (flapping) because fuselage attitude is only a proxy. — [Conner, Burley, Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
- **Ambient conditions.** At constant indicated airspeed and gross weight, in-plane thickness-noise pulse amplitude rose 4.7 dB (moderate speed) and 5.4 dB (high speed) between sea level and 7000 ft (AS350 SD1; Salton Sea, Amedee 4000 ft, Sweetwater 7000 ft; 2014–15; 1510 data points). Nondimensional flight conditions (CW, MHe) collapsed these, and BVI also changed with altitude. — [Watts, Greenwood, Stephenson, "Measurement and Characterization of Helicopter Noise at Different Altitudes," AHS 72nd Forum 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf); theory in [Greenwood & Schmitz, J. Aircraft 51(1), 2014] and [Greenwood, Sim, Boyd, JAHS 69, 012006 (2024)](https://s3.amazonaws.com/na-st01.ext.exlibrisgroup.com/01ALLIANCE_WSU/storage/alma/06/A4/F0/CB/0C/A7/2A/3D/D0/71/E0/C9/37/84/12/7B/77-jahs2274_pub.pdf) (link is a time-limited signed URL; cite by DOI/JAHS)
- **Wind as a source-state error, not only a propagation error.**
  - In an 80 KIAS descent with a 15 kt crosswind, holding the fuselage on the inertial path cuts BVI by about 11 dB, while crabbing to zero sideslip leaves it essentially unchanged. A headwind gradient (15 kt at 300 ft to 0 at ground) cuts peak BVI about 9 dB, against about 1 dB for an equal crosswind gradient.
  - Directivity rotates toward the wind, more at 50 KIAS than 130 KIAS. Lower-harmonic noise shifts only about 7°.
  - The authors recommend crosswind-oriented test runs for database fidelity.
  - [Greenwood & Sim, "The Effects of Crosswind Flight on Rotor Harmonic Noise Radiation," J. Aircraft 55(5), 2018 (PMC summary)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6662637/)
- **Run-to-run source variability (measurement-side "error").**
  - XV-15 repeated level flights: 106.7, 106.6, 99.6 and 107.1 dB at one microphone 10 s before overhead, a 7.5 dB spread. A nominal prediction looked 3 dB under or 4 dB over depending on which run was "blindly" chosen. — [Boyd, "Challenges in Rotorcraft Acoustic Flight Prediction and Validation," AIAA 2003-0709](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf)
  - Snapshot arrays "captured significant variation in acoustic emissions throughout a single run and between multiple runs of similar conditions." — [Stephenson et al., CEAS Aeronaut. J. 2023 (abstract)](https://consensus.app/papers/details/bb3903fd8fcf582dabbbdecda80000a6/?utm_source=claude_desktop)
- **Propagation.**
  - RNM uses straight rays without wind; Boyd notes that echoes from canyon walls or buildings are not included. — [Boyd 2003](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf); [Conner & Page 2002](https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-2002-ahsitsm-dac.pdf)
  - NORAH uses Chien–Soroka ground reflection with Delany–Bazley impedance and no refraction. — [Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
  - The NORAH2 validators say refraction "have not been addressed in this project" and that weather effects "will have a very strong influence on the effective ground attenuation for long distance propagation over flat soft ground." Turbulence (in Nord2000) raises levels 2–3 dB there. — [EASA D3 2024](https://www.easa.europa.eu/en/downloads/139664/en)
  - Over-water vs over-ground propagation at 4–6° elevation changed SEL by about 2 dB (fixed-wing data, DCA). — [ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **Low-frequency content and metric choice.** Bands down to 10 Hz are required for rotorcraft. Omitting them does not change SEL(A), LAmax or EPNL but does change C-weighted and unweighted metrics. — [ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **Directivity flattening in integrated tools (eVTOL).** Calibrated to match LmaxA at 1000 ft, AEDT and AAM still differed by 1.7 dB SELA directly under the track, "attributable to polar spectral directivity that accumulates in the SELA calculation." — [Rizzi, Page et al., "Comparison of two community noise models applied to a NASA UAM concept vehicle," Inter-Noise 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)
- **Hemisphere radius / near field and Doppler / convective amplification.**
  - NORAH uses 60 m, matched to ICAO Doc 9501's frequency-extrapolation distance. — [Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
  - Joby's hemisphere radius is 100 ft (30 m). — [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
  - Convective amplification is "a small effect (≈1 dB or less)" for the Joby flyover hemispheres and was not corrected. — [Pascioni, Thai, Bain, "Propeller Source Noise Separation from Flight Test Measurements of the Joby Aviation Aircraft," Aeroacoustics 2024](https://ntrs.nasa.gov/api/citations/20240004305/downloads/Aeroacoustics2024_NASApaper_Pascioni_v4.pdf)
- **Ground-plate / flush-microphone correction.** The Joby processing subtracts a uniform 6.02 dB (pressure doubling) "independent of frequency and emission angle." — [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
- **Ambient / SNR.** Joby data were removed where vehicle levels were less than 3 dB above ambient, and a frequency cutoff was set where the vehicle crossed ambient + 3 dB. — [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)

### Inferences
- Source-state errors (deceleration, load factor, wind-relative state) are an order of magnitude larger than most propagation differences at high elevation. Once the state is right, the floor is set by run-to-run variability.
- Two items are flagged in the literature without quantified effect on footprints: (a) the flat 6 dB flush-plate correction, which is accurate only at high frequency and near-normal incidence; and (b) constant-value extrapolation toward the rotor plane (NORAH) or a fixed rim. For a tool whose long-range error depends on the database rim (internal NICE-OPS result above), this is a measurable gap.
- The crosswind result means a ground-referenced sphere database used in wind mixes a source error of up to about 10 dB (BVI) with the propagation error. The air- vs ground-speed reference of the database therefore matters as much as refraction.

### Gaps
- No published quantification found of the error introduced by removing ground reflections during de-propagation (inverse ground-effect) for elevated microphones, as opposed to ground boards.
- No published quantification found of hemisphere-radius / near-field error for helicopter spheres, i.e. 60 m vs 30 m vs larger.
- No quantified INS/GPS tracking-error propagation into sphere-model SEL found, beyond RNM's attitude-sensitivity study.
- No published study isolates turbulence as an error source in sphere-based footprint validation against measurement; the NORAH2/Nord2000 result is model-to-model.

## 3. Specific campaigns

### Takeaway
The US government campaigns produced the main sphere databases:
- Eglin 2003/2005 for RNM.
- Bell 430 maneuver test, 2011, for FRAME-QS.
- 2014–15 altitude test.
- 2017 light-helicopter Noise Abatement test at Amedee.
- 2019 medium-sized test at Coyle Field.
- 2021 Joby test (NASA AAM National Campaign).
- A later Army/NASA/Navy snapshot-array test.

Published validations from them are mostly guidance-oriented (Fly Neighborly) rather than model-error studies. In Europe, NORAH/NORAH2 (DG-MOVE and EASA) measured 8 + 3 helicopter types for hemispheres but validated against certification and Nord2000 only.

### Cited Findings
- **2017 light helicopters and 2019 medium helicopters.** "Ten distinct vehicles were tested at three different test ranges"; nine are documented. Twelve turning-condition pairs were evaluated on peak band-limited (50–2500 Hz) LA,max as a BVI surrogate. Light: R44, R66, Bell 206, Bell 407, AS350, EC130. Medium: Bell 205, S-76D, AW139 (MH-65 also flown). The 2019 medium-weight campaign was at Coyle Field. — [Stephenson et al., VFS 76 (2020)](https://ntrs.nasa.gov/api/citations/20205003220/downloads/VFS76_Paper_v7.pdf); [Pascioni et al. 2020](https://ntrs.nasa.gov/api/citations/20205008035/downloads/Pascioni_VFS2020.pdf)
- 2017 data report: Watts et al., "Noise Abatement Flight Test Data Report" (NASA TM, 2019). — [ResearchGate listing](https://www.researchgate.net/publication/333320599_Noise_Abatement_Flight_Test_Data_Report) (not read)
- The resulting Fly Neighborly guidance covers the AS350, AS365, AW139, Bell 205, 206 and 407, R44, R66 and S-76D, modeled in AAM with flight-test spheres. — [Page et al., Inter-Noise 2021 (abstract)](https://consensus.app/papers/details/3ce8695f4b3c530586c69e924047f80f/?utm_source=claude_desktop)
- **Altitude campaign (Sep 2014–Feb 2015).** AS350 SD1 and EH-60L at Salton Sea (0 ft), Amedee (4000 ft) and Sweetwater (7000 ft); 165 flight hours, 1510 data points. One stated purpose was to "establish the statistical variability in acoustic flight testing of helicopters." — [Watts, Greenwood, Stephenson 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf)
- **Eglin 2003 ("Acoustics Week") and 2005 (Heavy Lift Rotorcraft Acoustics).** Benchmark semi-sphere databases for RNM/ART. No testing was done with ground winds or crosswind components above 10 kt. Repeated MD600N housekeeping runs gave ±0.5 dBA SEL. — [Conner, Burley, Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
- **Bell 430 maneuver acoustic flight test** (Watts et al., NASA TM-2014-218266): the data behind FRAME-QS near-horizon validation. — cited in [Watts et al. 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf) and [Greenwood et al. 2015](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005979.pdf)
- **Joby / NASA AAM National Campaign (2021, Big Sur; published 2022).** 58 channels, more than 100 test points (31 unique conditions), 100 ft hemisphere radius, winds averaging under 3 kt.
  - Flyover repeatability was 1–2 dB (60 kt: two flyovers on different days; 100 kt: six flyovers over two days).
  - Port/starboard symmetry was within 2 dB.
  - Hover varied 2–5 dB within a run, with back-to-back disparities.
  - Wing-borne 103 kt was more than 6 dB quieter than 61 kt semi-thrust-borne; 45.2 dBA scaled to 500 m.
  - [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
- **Snapshot arrays (US Army/NASA/Navy).** Hemispheres from snapshot arrays, put through land-use software, matched traditional linear-array SEL within 0.1 dBA at the sideline but differed by 1.6 dBA on the centerline. The authors estimate at least 51 microphones are needed without changes to the land-use software, and call for a semi-empirical interpolation in place of linear frequency weighting. — [Stephenson et al., CEAS Aeronaut. J. 2023 (abstract)](https://consensus.app/papers/details/bb3903fd8fcf582dabbbdecda80000a6/?utm_source=claude_desktop)
- **NORAH (DG-MOVE).** Hemispheres for 8 types representing about 70% of the European fleet, with 15 acoustic classes. Missing types are handled by a class-level offset from certification levels, mirrored for opposite rotor rotation. — [Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
- **NORAH2 (EASA, released Dec 2024).** New hemispheres for the Cabri G2, A109 and S-92. — [EASA D3 2024](https://www.easa.europa.eu/en/downloads/139664/en); [EASA project page](https://www.easa.europa.eu/en/research-projects/environmental-research-rotorcraft-noise)
- A NORAH2 urban demonstration for Cannes exists. — [Aerospace 13(1):37, 2026](https://doi.org/10.3390/aerospace13010037) (not read in full)
- **HELENA / FRIENDCOPTER (Airbus Helicopters).** A hemisphere/noise-carpet tool, compared with measurements on the EC130 and EC135 in steady single-segment procedures; numbers not retrieved. — [ACRP 02-44](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf); [ResearchGate listing](https://www.researchgate.net/publication/289803872_Validation_of_EC130_and_EC135_environmental_impact_assessment_using_HELENA)
- Its successor workflow, MOTUS, corrects measured hemisphere databases with simulation to represent new technology. — [Guntzer et al., VFS 80th Forum 2024](https://consensus.app/papers/details/c87e294065ce55bc914a9d7f9c7bfb6e/?utm_source=claude_desktop); [CEAS 2025](https://link.springer.com/article/10.1007/s13272-025-00876-5)
- **Other hemisphere-from-flight-test work.** An 84-microphone KNAC array built hemispheres for an electric unmanned helicopter from more than 30 level flights, combining flight directions and altitudes "under low wind conditions." — [Cho et al., AIAA/CEAS Aeroacoustics 2024 (abstract)](https://consensus.app/papers/details/473208349f6d5063b5d07f94777f8431/?utm_source=claude_desktop)

### Inferences
- The 2017 campaign's held-out maneuvering and approach runs appear unused in the literature as a quantitative sphere-model validation set. Published use is guidance tables and AAM operational plots. This is an open niche, and NICE-OPS's harness is positioned to fill it.

### Gaps
- No quantified results found for the Eglin 2007/2010 tests, the FAA/Volpe HNM, FRIENDCOPTER/HELINOVI validation numbers, the Airbus NAP validation, or DLR measurement comparisons; these were either not located or paywalled.
- Uber/NASA eVTOL hemisphere measurement: not found.

## 4. Uncertainty quantification, run-to-run variance, weather sensitivity

### Takeaway
Measured repeatability of steady-flight source levels is about ±0.5 dBA SEL (well-guided, low-wind level flight) to 1–2 dB (eVTOL flyovers). It degrades to 2–5 dB in hover and up to 7.5 dB at single microphones under BVI-prone or poorly controlled states. Ambient conditions and wind shift source levels by about 5 dB (altitude) and up to 9–11 dB (BVI in wind), unless the database is nondimensional and air-referenced. No formal Monte-Carlo UQ of a sphere-based footprint against flight test was found.

### Cited Findings
- **±0.5 dBA SEL** across seven MD600N housekeeping-run spheres (110 kt, 150 ft, two test days), attributed to wind-driven control and state changes; cockpit guidance was used and winds were capped at 10 kt. — [Conner, Burley, Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
- **7.5 dB spread** among four nominally identical XV-15 level flights at one microphone. "Conclusions drawn from a prediction would depend greatly on which of the four repeated runs were chosen." — [Boyd 2003](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf)
- **Joby:** 1–2 dB flyover repeatability across days; 2–5 dB within a hover run; about 5 dB between back-to-back HIGE runs. HIGE altitude differed by under 2 ft (5.7 vs 7.3 ft AGL), and one propeller's RPM differed by about 5%, which alone should give under 1 dB. — [Pascioni et al. 2022](https://ntrs.nasa.gov/api/citations/20220006729/downloads/Aeroacoustics2022_Pascioni_STRIVES5.pdf)
- **Nondimensional flight conditions reduce both the mean shift and the spread.** "The deviation is less for the non-dimensional flight conditions than it is for the constant IAS conditions" (AS350 SD1 and EH-60L medians, means and standard deviations, Tables A14, A15 and A18). — [Watts, Greenwood, Stephenson 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf)
- **Wind / crosswind:** 11 dB BVI change with piloting technique in a 15 kt crosswind, 9 dB with a headwind gradient. — [Greenwood & Sim 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6662637/)
- Boyd attributes the greater variability of flight test, relative to wind tunnels, to "real world factors" that are difficult to control. — [Boyd 2003](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf)
- **Within-run variability seen by snapshot arrays:** "significant variation in acoustic emissions throughout a single run and between multiple runs of similar conditions." — [Stephenson et al. 2023](https://consensus.app/papers/details/bb3903fd8fcf582dabbbdecda80000a6/?utm_source=claude_desktop)
- A machine-learning classifier of "noisy" operating regions was trained on more than 1500 steady conditions for ten helicopters and generalised to helicopters held out of training. It addresses the operating-state dimension of uncertainty, not dB error bars. — [Greenwood, "Estimating Helicopter Noise Abatement Information with Machine Learning," AHS 74th Forum 2018](https://ntrs.nasa.gov/api/citations/20180006308/downloads/20180006308.pdf)

### Inferences
- A sphere-model validation claiming better than about 1 dB SEL agreement for single runs is below the measured repeatability floor in all but the best-controlled level flights. Approaches and hover should be judged against a 2–5 dB envelope, or better, against the spread of repeated runs.
- Weather sensitivity enters twice, through the source state (air- vs ground-relative μ, γ, CT) and through propagation. Published UQ treats them separately; no study combining them was found.

### Gaps
- No published Monte-Carlo or Bayesian UQ of a hemisphere-based footprint against flight test found (searched Greenwood/Rau, NASA, PSU). Greenwood's CV lists ambient-condition and ML work but no explicit sphere UQ paper.
- No quantified sensitivity of sphere-model SEL to measured temperature/wind profiles (refraction) against flight test found in the hemisphere literature.

## 5. Recommended practices and proposed improvements in the literature

### Takeaway
The literature converges on these improvements:
1. Parameterise spheres nondimensionally (μ, CT or CW, effective flight-path angle, MH) and select them quasi-statically with load factor (Q-SAM / FRAME-QS).
2. Fly the database air-referenced and crosswind-oriented.
3. Widen lateral and near-rotor-plane coverage and measure tip-path-plane orientation.
4. Use semi-empirical or physics-based extrapolation instead of constant-value fill (FRAME, MOTUS, a semi-empirical snapshot interpolation).
5. Add refraction and turbulence where long-range or soft-ground propagation matters.
6. Validate interpolation, turns and hover against operational measurements, which no regulatory tool has yet done.

### Cited Findings
- **Q-SAM with CT:** "Accurate noise predictions must account for the significant changes in the rotor tip-path-plane orientation and load factor that occur during maneuvering flight," and FRAME can precompute spheres for load factors "difficult to measure using more empirical source noise characterization approaches." — [Greenwood et al. 2015/2020](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160005979.pdf)
- **Nondimensional databases scale across ambient conditions.** — [Watts et al. 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf); [Greenwood, Sim, Boyd JAHS 2024](https://s3.amazonaws.com/na-st01.ext.exlibrisgroup.com/01ALLIANCE_WSU/storage/alma/06/A4/F0/CB/0C/A7/2A/3D/D0/71/E0/C9/37/84/12/7B/77-jahs2274_pub.pdf)
- **Crosswind-oriented runs and crabbing** for database fidelity. — [Greenwood & Sim 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6662637/)
- **Measure the tip-path plane, not just fuselage attitude**, when building spheres. — [Conner, Burley, Smith 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)
- **ACRP 02-44 recommendations:**
  - Lateral directivity with BVI-resolving angular fidelity and roll-angle extent.
  - One-third-octave bands from 10 Hz.
  - A "hybrid analytical-empirical database development" approach, which received "general agreement in principle" from Boeing, UMD, NASA, Sikorsky and Airbus.
  - [ACRP 02-44, 2015](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf)
- **NORAH2 validators** recommend refraction modelling for long-distance soft-ground cases (noting Nord2000 is about 100x slower) and in-depth tests of hemisphere interpolation, turns and hover/taxi against measurement. — [EASA D3 2024](https://www.easa.europa.eu/en/downloads/139664/en)
- **NORAH 2018** notes that because the method is one-third-octave and separates source from propagation, "more sophisticated methods (e.g. raytracing or shielding) can easily be implemented." — [Euronoise 2018](https://www.euronoise2018.eu/docs/papers/35_Euronoise2018.pdf)
- **Snapshot arrays** for capturing non-steady directivity, with at least 51 microphones and a semi-empirical interpolation in place of linear frequency weighting. — [Stephenson et al. 2023](https://consensus.app/papers/details/bb3903fd8fcf582dabbbdecda80000a6/?utm_source=claude_desktop)
- **Alternatives to tabulated hemispheres:** acoustic-mode / near-field holography representations report under 0.3 dB average error against numerical reference and one-fifth the cost of the hemisphere method for 100,000 observers (simulation only). — [Shao et al., Chinese J. Aeronautics 2024 (abstract)](https://consensus.app/papers/details/a147875b4ae45d24bec92d956e44d975/?utm_source=claude_desktop)
- **Boyd:** add wind and urban echoes to RNM. Flight comparisons require examining time histories and vehicle state, not blind run selection. — [Boyd 2003](https://ntrs.nasa.gov/api/citations/20030002351/downloads/20030002351.pdf)

### Inferences
- Against this state of the art, a sphere-based tool with the following would already exceed what the published regulatory tools (NORAH2, AEDT) document:
  - nondimensional (μ, γ_e, CT) blending,
  - Q-SAM with load factor,
  - air-referenced conditions,
  - refraction (linear, stratified, PE) and turbulence-scattering options,
  - point-receiver and plate ground models,
  - held-out flight-test validation.
- The residual gaps the literature points to and has not solved are:
  - (a) source coverage at the low-speed / hover end of approaches and near the rotor-plane rim,
  - (b) run-to-run variability as an explicit error bar,
  - (c) validation of interpolation and turns against operations.
  The first two match the internal NICE-OPS day-286 findings: slow-end over-prediction and rim-dependent long-range error.

### Gaps
- Spherical-harmonic or other physics-based extrapolation applied to measured helicopter hemispheres, with flight-test validation, was not found within this search budget.
- No published coupling of measured hemispheres with PE/ray propagation validated against rotorcraft flight test was found, other than the NORAH2/Nord2000 model-to-model comparison and (unpublished) NICE-OPS work.
