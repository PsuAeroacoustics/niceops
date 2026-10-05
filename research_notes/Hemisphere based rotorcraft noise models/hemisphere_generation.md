# Generation, Parameterization and Interpolation of Rotorcraft (and eVTOL/UAM) Source Noise Hemispheres for Footprint Simulation

Scope note: these notes cover how the published literature acquires or predicts source-noise hemispheres, what it uses to index them by operating state, how it interpolates them, which frame and Doppler conventions it uses, and where it reports the method failing. Most primary sources were read in full text (NTRS PDFs, UMD dissertations, an FAA/ASCENT report). A few are cited from abstracts or publisher listings only, and those are marked. Eric Greenwood (the user) is an author on many of these; his bibliography was checked against the PSU Aeroacoustics Lab publication list ([PSU Aeroacoustics publications](https://aeroacoustics.psu.edu/publications/)).

## 1. Flight-test hemisphere acquisition: arrays, de-propagation, tracking, SNR handling, radius choice, coverage gaps

### Takeaway
The standard method flies one steady condition through a linear "goal-post" ground array perpendicular to the track. Received spectra are "frozen" to emission angles and de-propagated to a fixed-radius sphere, correcting for spherical spreading, absorption, and ground reflection (−6 dB pressure doubling for ground boards). Recent NASA/Army work replaces this with 2-D "snapshot" arrays that capture all angles at one instant. Coverage within about 10° of the horizon (the tip-path plane) and above it is chronically missing, and that is where the measured hemisphere is least reliable.

### Cited Findings
**Classic linear-array method (RNM/ART)**
- Conner, Burley & Smith describe the RNM/ART (Acoustic Repropagation Technique) process ([Conner, Burley & Smith, "Flight Acoustic Testing and Data Acquisition for the Rotor Noise Model (RNM)," AHS 62nd Annual Forum, Phoenix, May 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)):
  - The aircraft flies at a constant operating condition through a linear microphone array perpendicular to the ground track. Spectra are computed at fixed intervals ("typically every 0.5 seconds").
  - Each spectrum is related to the aircraft position relative to each microphone, giving levels against emission angle, which are de-propagated to "a semi-sphere of selected radius."
  - RNM's propagation physics (spherical spreading, atmospheric absorption, ground reflection, Doppler shifts) is used in reverse.
  - The typical sphere grid has 5° angular resolution.
- Eglin 2003/2005 array geometry from the same paper ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)):
  - The 2005 test used a 22-microphone array: two cranes each suspending four microphones (highest at 175 ft AGL, cranes 800 ft apart), 12 ground-board microphones between them, and 1200 ft sideline microphones.
  - Elevated "goal-post" microphones were added specifically to improve in-plane (near-horizon) coverage.
  - Flying 25 ft below the highest crane microphones gave data above the tip-path plane (φ > ±90°).
- On coverage gaps, Conner et al. report ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)):
  - For a 90 kt, 6° descent at 250 ft, "a significant number of semi-sphere grid points lie outside the region where measured data was acquired, especially near the rotor tip-path-plane to the port side."
  - Earlier ART versions filled those cells by holding the level constant from the nearest measured angle up to the tip-path plane. The authors judged that this "may or may not be acceptable."
  - The current ART "does not attempt to extrapolate." It writes a sentinel value of 1×10^35, which tells RNM not to use that part of the sphere.
- ART gives "a user selectable blend of Laplace and spline interpolation schemes" to fill the regular sphere grid from scattered emission points. A wider range of contributing points lets sparse data populate the sphere, "but with a reduced fidelity" ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- Repeatability: the SEL variation between "housekeeping" repeat runs was about ±0.5 dBA at all sideline distances, attributed to wind variation ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- Joby/NASA eVTOL hemispheres (2021-22) were processed at a 100 ft radius covering emission directions below the aircraft. Repeated flight profiles gave "many readings ... within 1 to 2 dBA of each other" (secondary trade-press source quoting Pascioni) ([Vertical Mag, "Findings from NASA and Joby's first published eVTOL acoustic test data"](https://verticalmag.com/features/nasa-first-published-evtol-acoustic-test-data-set-stage-future/)).

**Snapshot arrays (2022–2025)**
- Houston, Stephenson, Pascioni & Stutz motivate the snapshot array ([Houston, Stephenson, Pascioni & Stutz, "Snapshot Array Design Considerations for Rotorcraft Noise Characterization," VFS 81st Annual Forum, Virginia Beach, May 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)):
  - UAM/FVL vehicles "may not maintain a steady aerodynamic state during flight, making source hemispheres measured using traditional linear arrays unreliable or difficult to interpret."
  - Independently controlled propulsors "allow for non-unique trim states for the same operating condition," citing Joby data.
  - Their solution is to measure all emission angles simultaneously in a 2-D array.
- Array design (MD530F, July 2024, same paper):
  - Microphone positions come from a Lebedev-quadrature distribution projected onto the ground from 61 m (200 ft), with 25-, 36-, 37-, 48-, 61- and 73-element subsets embedded.
  - The base distribution has "no points within 10° of the horizon plane." A ring of 12 microphones at 7.5° from the horizon and six beneath the vehicle were added.
  - Ground boards are 400 mm. GPS timing syncs the stations. Tracking uses the Army/NASA second-generation Aircraft Navigation and Tracking System.
- De-propagation (same paper):
  - SPL_hemi = SPL_mic + A_spread + A_atm + A_ground. Hard ground boards are modelled as simple pressure doubling, A_ground = −6.02 dB.
  - Straight-line propagation is used because the distances are short (≈91–457 m).
  - Each snapshot uses the 0.5 s of data centred on the array-centre passage. Pressures sampled at 50 kHz are windowed around each microphone's computed reception time t_i,0, and spectra are formed with Welch's method (Hamming, 50% overlap).
- SNR floor handling (same paper):
  - At high frequency near the horizon, adding back absorption "can pose a problem when the levels added back are greater than the original signal." Some runs were modelled as having more than 20 dB of attenuation at 3 kHz over 457 m.
  - Fix: above 250 Hz, each spectrum is compared with an ambient spectrum raised by 6 dB. Above the intersection frequency A_atm is set to zero, though spreading is still applied.
- Sphere storage (same paper): radius 30.48 m (100 ft), "centered about the vehicle, typically at the main rotor hub for conventional helicopters" ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- Accuracy against channel count (leave-one-out cross-validation, 75 kt level flight; [Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)):
  - Point errors range from ≈0 to more than 14 dB at low resolution.
  - Median error falls from 1.75 dB (IQR 0.87–3.26) with two microphones inside a 30° radius of influence to 0.31 dB (IQR 0.14–0.46) with 25.
  - Adding the 12-microphone horizon ring to the 25-element array cut the median from 1.68 to 1.48 dB.
  - The 61- and 73-element arrays gave medians of 0.71 and 0.89 dB.
  - The largest errors are "to the rear of the aircraft near the horizon."
  - Main-rotor BPF (36–44 Hz) errors are higher and more uncertain than broadband (500 Hz–5 kHz).
  - The 61-element array was found adequate to "qualitatively" match linear-array hemispheres in steady level and descending flight.
- An earlier 21-microphone snapshot array, with positions chosen by hand, was flown in 2022 ([Stephenson & Houston, "Rotorcraft Source Noise Characterization via Acoustic Snapshot Array: Development and Evaluation," CEAS Aeronautical Journal 15, 659–670, 2023, DOI 10.1007/s13272-023-00672-z](https://link.springer.com/article/10.1007/s13272-023-00672-z); cited via Houston et al. 2025, full text not accessed).

**Vertical arrays and beamforming (NLR)**
- NLR measured an AH-64 with a 100 m vertical array of 93 microphones plus 13 ground microphones ([Tuinstra & Sijtsma, "Measurement of Helicopter Noise Hemispheres Utilizing a 100m Vertical Array," European Rotorcraft Forum paper, NLR](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)):
  - Delay-and-sum beamforming focused on the helicopter (up to 500 Hz) filtered ground reflections out of the vertical-array signals.
  - Their design rules: R ≫ helicopter size so that it is a point source; R < ≈450 m to limit meteorological error; and avoid ground microphones for lateral angles beyond about 60°, using a vertical array instead.
- Directivity was then fitted by regression. Two-variable Taylor (polynomial) series were better at low order, Fourier series better at following rapid directivity changes at high order. Expected-error measures were calibrated to "weak (±1dB), medium (±2dB) and strong (±6dB) scatter" ([Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)).

**Upper hemisphere and near-horizon data**
- In Greenwood, Sim & Boyd's AS350 hemispheres (Lambert projection), "Data within 10° of the plane of the horizon are not available because the acoustic measurements are made from the ground" ([Greenwood, Sim & Boyd, "The Effects of Ambient Conditions on Helicopter Harmonic Noise Radiation: Theory and Experiment," AHS 72nd Forum 2016; JAHS 69(1), 2024, DOI 10.4050/JAHS.69.012005](https://ntrs.nasa.gov/api/citations/20160009056/downloads/20160009056.pdf)).
- Schmitz et al. combined a helicopter-mounted in-flight microphone boom with a ground linear array on a Bell 206B ([Schmitz, Greenwood, Sickenberger, Gopalan, Sim, Conner, Moralez & Decker, "Measurement and Characterization of Helicopter Noise in Steady-State and Maneuvering Flight," AHS 63rd Forum, Virginia Beach, May 2007](https://ntrs.nasa.gov/api/citations/20070021458/downloads/20070021458.pdf)):
  - The test used DGPS/INU, direct TPP measurement, a pursuit guidance display, and testing "only in very low winds."
  - The authors conclude that microphones close to the helicopter can estimate far-field-radiating sources "if the microphones are positioned correctly relative to the noise source."

**Source separation and de-Dopplerization as hemisphere pre-processing**
- Greenwood & Schmitz use time-domain de-Dopplerization ([Greenwood & Schmitz, "Separation of Main and Tail Rotor Noise from Ground-Based Acoustic Measurements," J. Aircraft 51(2):464–472, 2014, DOI 10.2514/1.C032046](https://pure.psu.edu/en/publications/separation-of-main-and-tail-rotor-noise-from-ground-based-acousti/); ERF 2009 conference version):
  - Ground signals are transformed to "virtual in-flight microphones traveling with the helicopter."
  - They are then periodically averaged with main and tail rotor blade passages, separating the two rotors' harmonic noise. Bell 206B3 data show strongly different directivities for the two rotors.
- Rachaprolu & Greenwood replace harmonic averaging with a Vold–Kalman order-tracking filter after de-Dopplerization ([Rachaprolu & Greenwood, "Helicopter Noise Source Separation Using an Order Tracking Filter," JAHS 69(1), 2024, DOI 10.4050/JAHS.69.012006](https://rosap.ntl.bts.gov/view/dot/72623/dot_72623_DS1.pdf)):
  - This preserves phase and amplitude and works even for time-varying sources, including BVI.
  - Bell 430 data were split into main rotor, tail rotor and broadband residual, and depropagated to separate hemispheres per component.
  - The de-Dopplerization uses higher-order interpolation than Greenwood & Schmitz "to limit high-frequency distortion."
- Pascioni, Thai & Bain applied the Vold–Kalman filter to Joby eVTOL flight data, using per-propeller time-varying RPM and a single microphone channel ([Pascioni, Thai & Bain, "Propeller Source Noise Separation from Flight Test Measurements of the Joby Aviation Aircraft," 30th AIAA/CEAS Aeroacoustics Conference, Rome, June 2024](https://ntrs.nasa.gov/citations/20240004305); abstract only):
  - This produced per-propeller noise hemispheres and ranked propulsors "without the use of phased arrays."

### Inferences
- Every measured database inherits a structural hole near and above the TPP horizon. The literature handles it three ways:
  - Flag the cells and refuse to use them (ART's 10^35 sentinel).
  - Fill them by holding a level or extrapolating (old ART), which the authors themselves questioned.
  - Fill them from a physics model (FRAME; see section 3).
- A footprint tool should know which cells of its sphere were measured. Low-elevation, long-range observers, which are exactly where footprints extend, are served by the least-certain cells.
- Absorption correction at long de-propagation ranges can manufacture high-frequency energy out of the noise floor. Houston et al.'s ambient+6 dB cutoff is one published guard; a sphere builder without one can inflate high-frequency levels on near-horizon cells.
- Reported uncertainty floors are roughly ±0.5 dBA SEL for repeat runs (Conner 2006), 1–2 dBA for eVTOL repeats (Joby, trade press), and sub-dB to more than 3 dB median interpolation error on the sphere depending on array density (Houston 2025). Footprint differences below about 1 dB should not be over-interpreted.

### Gaps
- Not extracted: the snapshot-array CEAS paper's own quantitative snapshot-vs-linear comparison (full text was paywalled), and the Joby NASA paper's hemisphere-processing details (abstract and trade press only).
- No primary source found that quantifies the effect of tracking error (position or attitude) in degrees on hemisphere levels. Only Conner et al.'s orientation-method sensitivity is reported (section 5).
- The year and ERF number of the Tuinstra & Sijtsma vertical-array paper could not be confirmed. It cites 2004 work, so it is 2004 or later.
- Snider (Bell) hemisphere-acquisition papers were not located. Snider appears as a co-author of the Bell 430 maneuver test (Watts, Snider, Greenwood & Baden, AHS 68th Forum 2012; NASA TM-2014-218266 data report) per the [PSU list](https://aeroacoustics.psu.edu/publications/); those were not read.

## 2. Prediction-based hemispheres (CAMRAD II / ANOPP2 / AARON / PAS / FW-H)

### Takeaway
Predicted hemispheres are generated the same way everywhere: a comprehensive-analysis trim (CAMRAD II or similar) gives blade loads, and an FW-H Formulation 1A solver evaluates periodic thickness and loading noise at a sphere of observers at fixed radius. The result is usually lower-hemisphere only and tonal only. Broadband noise is added empirically, if at all. Tail-rotor content and broadband noise turn out to be essential for integrated metrics.

### Cited Findings
- Boyd, Greenwood, Watts & Lopes describe the CAMRAD II + ANOPP2 chain for an AS350 ([Boyd, Greenwood, Watts & Lopes, "Examination of a Rotorcraft Noise Prediction Method and Comparison to Flight Test Data," NASA/TM-2017-219370, Jan 2017](https://ntrs.nasa.gov/citations/20170002590)):
  - CAMRAD II trims; ANOPP2-F1A computes pressure time histories "for a set of observers on a hemispherical surface under the rotorcraft"; these are converted to 1/3-octave spectra.
  - An optional broadband hemisphere from a scaled Pegg model (main rotor only) is added.
  - The combined hemisphere is "flown" along a straight path with ANOPP2-AnoppGround, which accounts for "hemisphere motion ('Doppler effect'), atmospheric absorption and attenuation, and spherical spreading," at 0.5 s output intervals.
- Findings from the same report ([Boyd et al. 2017](https://ntrs.nasa.gov/citations/20170002590)):
  - Main rotor alone "will produce results that severely underpredict integrated metrics." Including tail-rotor frequency content "is essential."
  - With inputs calibrated on level flight, the descent case's predicted A-weighted metrics were 4–5% low because the measured descent had more broadband content.
  - It is "unclear whether this calibration method is valid for various flight conditions of the same vehicle or for vehicles of other size and weight classes."
- For the NASA RVLT quadrotor, the ANOPP2 AARON tool (F1A, periodic loading and thickness noise "under a quasistatic operating condition") produced source spheres at fixed radius from the CG ([Rizzi, Page & Cheng (inferred authorship: Volpe filename "InterNoise-2021_Rizzi-Page-Cheng"), "Comparison of two community noise models applied to a NASA urban air mobility concept vehicle," InterNoise 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)):
  - The mesh is 10° (θ = 0° nose to 180° tail; φ from −90° port to +90° starboard). Only the lower hemisphere was computed.
  - Output was 1/3- and 1/12-octave bands in NetCDF for Volpe AAM. This is the "Gen 1.2" database.
- Directivity simplification test (same paper): spheres were calibrated to AEDT NPD L_Amax at 1000 ft. An omnidirectional sphere matched L_Amax under-track but SEL_A differed by 1.7 dB, attributed to polar directivity accumulating in SEL. With full 3-D spheres, AEDT increasingly over-predicted at lateral positions (to 5000 ft sideline) because it lacks azimuthal directivity ([IN21](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)).
- Pascioni, Rizzi & Schiller store, per propeller, the magnitude and phase of Fourier-transformed FW-H pressure on source hemispheres at R_s = 12R and 5° resolution, truncated above the 5th BPF ([Pascioni, Rizzi & Schiller, "Noise Reduction Potential of Phase Control for Distributed Propulsion Vehicles," AIAA SciTech 2019 (venue per the AIAA-member footnotes and NTRS record; not printed in the PDF text)](https://ntrs.nasa.gov/api/citations/20200002397/downloads/20200002397.pdf)):
  - Complex pressures are interpolated at the emission point on each hemisphere, then corrected for propagation time and spreading, and summed coherently.
- Krishnamurthy, Rizzi & Cheng use per-rotor source noise hemispheres of one-revolution pressure histories at evenly spaced points ([Krishnamurthy, Rizzi & Cheng, "Prediction-Based Auralization of a Multirotor Urban Air Mobility Vehicle," AIAA SciTech 2021](https://ntrs.nasa.gov/api/citations/20205010694/downloads/20205010694_ST21_Paper.pdf)):
  - "The center of the source noise hemisphere must be centered around each rotor. Rotor displacement from vehicle pitch and roll must be accounted for in the hemisphere centering."
  - Isolated-source runs were needed so that ANOPP2 would not use one hemisphere for all rotors.
- Krishnamurthy, Tuttle & Rizzi compute an F1A solution directly at the emission point tracked across the hemisphere, rather than interpolating a precomputed periodic hemisphere, so that aperiodic (unsteady) rotor noise can be auralized ([Krishnamurthy, Tuttle & Rizzi, "Auralization of Unsteady Rotor Noise using a Solution to the Ffowcs Williams-Hawkings Equation," VFS 75th Annual Forum, 2019](https://ntrs.nasa.gov/api/citations/20200002468/downloads/20200002468.pdf)).
- The FRAME-QS database of "nearly 500 predicted acoustic spheres" (40–130 kt, γ from +3° to −18°, n_z 0.75–2.0) is itself prediction-based, from a measurement-calibrated physics model ([Greenwood, Rau, May & Hobbs, AHS 71st Forum 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).

### Inferences
- Predicted spheres solve the coverage problem (any angle, including the upper hemisphere, can be computed) but usually miss broadband, tail-rotor, interactional and BVI-wake fidelity. Measured spheres have the opposite profile. Hybrid approaches such as FRAME and calibrated CAMRAD exist for this reason.
- Prediction-based multirotor work stores complex pressure (amplitude and phase) per rotor so that rotors can be summed coherently. A single dB-only vehicle sphere cannot represent a change in relative rotor phase (section 3).

### Gaps
- No CHARM-generated or PSU-WOPWOP-generated hemisphere database paper was found in this session with quantified validation. Those workflows exist in the literature but were not verified here.
- No published head-to-head comparison was found of a prediction-based against a measured hemisphere for the same vehicle and condition, with dB statistics over the full sphere.

## 3. Parameterizing the operating state: airspeed/μ, flight-path angle, C_T, Q-SAM effective flight-path angle, FRAME, non-dimensional scaling, transients, multirotor/eVTOL, rotor phase

### Takeaway
Empirical tools (RNM) index spheres by dimensional airspeed and flight-path angle, plus nacelle angle for tiltrotors. Q-SAM maps accelerating and maneuvering states onto steady ones through an effective flight-path angle that gives the same TPP angle of attack. FRAME shows the governing set is non-dimensional (C_T or C_W, μ, M_H or effective tip Mach, inflow λ). Data indexed dimensionally (IAS, weight) do not transfer across ambient conditions. Transients are well served quasi-statically if load factor/C_T is carried. For multirotor eVTOL, non-unique trim and relative rotor phase break the one-condition-one-sphere premise.

### Cited Findings
**Dimensional indexing (RNM)**
- RNM sound-hemisphere files carry the flight condition as attributes. For tiltrotors there are three independent variables (flight path angle, pylon angle, airspeed); "for conventional helicopters, only flight path angle and airspeed would be used." There may be separate broadband and narrowband (pure-tone) hemispheres per condition ([Lucas & Marcolini, "Rotorcraft Noise Model," Wyle/NASA LaRC, c. 1997, NTRS 20040110391](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)).

**Q-SAM**
- Gopalan's dissertation establishes the quasi-static equivalence between longitudinal decelerating flight and steady descent ([Gopalan, "Quasi-Static Acoustic Mapping of Helicopter Blade-Vortex Interaction Noise," Ph.D. dissertation, Univ. of Maryland, 2004 (advisor F. H. Schmitz)](https://api.drum.lib.umd.edu/server/api/core/bitstreams/744acd6f-a007-4f7e-9e6f-1971df854a2b/content)):
  - "The entire operating state of the helicopter is shown to be similar during equivalent flight conditions at the same flight velocity."
  - It implements storage-based acoustic mapping extended to slowly maneuvering longitudinal flight for BVI, and shows multi-segment decelerating approaches to be effective for noise abatement.
- Q-SAM formulation as stated by Greenwood et al. ([Greenwood, Rau, May & Hobbs, "A Maneuvering Flight Noise Model for Helicopter Mission Planning," AHS 71st Forum, May 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf); journal version Greenwood & Rau, JAHS 65(2), 2020, DOI 10.4050/JAHS.65.022007, per [PSU list](https://aeroacoustics.psu.edu/publications/)):
  - The acoustic state is equivalent when TPP angle of attack, airspeed, RPM and load factor match.
  - To first order the TPP angle of attack depends on D/W, γ and a_lon/g. At constant airspeed this gives the effective flight-path angle γ_e = γ + a_lon/g.
  - Generalized to 3-D: sin γ_e = (v · n̂_a)/V, where n̂_a is the "inertial tip-path-plane" unit normal, i.e. the direction of the net specific force (a + g).
  - Sphere attitude: û = (v − (v·n̂_a)n̂_a)/V, sin Δθ = u_z, and a roll expression (their Eqs. 6–9).
- Q-SAM had been validated "in several experiments" and was integrated into RNM as RNM-QS ([Greenwood, Schmitz & Sickenberger, "A Semi-Empirical Noise Modeling Method for Helicopter Maneuvering Flight Operations," AHS 68th Forum, Fort Worth, May 2012; JAHS 60(2), 2015, DOI 10.4050/JAHS.60.022007](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf); [Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- Its premise is "that only slow and moderate perturbations (including maneuvering flight) are considered," which Schmitz et al. argue is "not very restrictive" because such maneuvers are the safest and most likely flown ([Schmitz et al. 2007](https://ntrs.nasa.gov/api/citations/20070021458/downloads/20070021458.pdf)).
- Steady turning flight: without data at different load factors, Q-SAM is limited to "small increases (<10%) in the load factor" ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- Experimental support from the same Bell 206B program ([Schmitz et al. 2007](https://ntrs.nasa.gov/api/citations/20070021458/downloads/20070021458.pdf)):
  - Acceleration during a −7.5° glideslope reduced BVI "on average, by up to 9 dB."
  - Measured TPP angles varied linearly with V², flight path angle and acceleration (standard deviation < 1.5°, with a constant ≈2° bias attributed to the air-data boom).

**FRAME**
- FRAME applies parameter identification to "a non-dimensional assumed analytical model of the rotor noise sources" ([Greenwood, Schmitz & Sickenberger 2012](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf); [Greenwood, "Fundamental Rotorcraft Acoustic Modeling from Experiments (FRAME)," Ph.D. dissertation, Univ. of Maryland, 2011](https://aeroacoustics.psu.edu/publications/); [Greenwood & Schmitz, "A Parameter Identification Method for Helicopter Noise Source Identification and Physics-Based Semiempirical Modeling," AHS 66th Forum 2010; JAHS 63(3), 2018, DOI 10.4050/JAHS.63.032001](https://aeroacoustics.psu.edu/publications/)):
  - Governing parameters: C_T, μ, M_H and λ_0.
  - Sources: a modified Beddoes prescribed wake with identified inflow and vortex-strength parameters feeding indicial aerodynamics and a compact-chord acoustic solver for BVI, plus thickness and loading noise.
- FRAME-D extends this to transients with a dynamic wake (vortex vertical position from an ODE of freestream plus induced inflow) and needs no extra modelling parameters, so "existing steady-flight wind tunnel and flight test data" can build maneuver models ([Greenwood, Schmitz & Sickenberger 2012](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf)).
- Validation on a Bell 206B ([Greenwood, Schmitz & Sickenberger 2012](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf)):
  - For a 1.8 g cyclic pull-up, peak amplitudes and pulse shapes were "well captured at three microphone locations."
  - Roll transients were "predicted fairly well" but "too high towards the advancing side."
  - No "super BVI" (vortex-bundling-driven levels beyond steady descent) was observed in the measured maneuvers.

**Non-dimensional scaling and ambient conditions**
- Greenwood & Schmitz on altitude effects ([Greenwood & Schmitz, "The Effects of Ambient Conditions on Helicopter Rotor Source Noise Modeling," AHS 67th Forum 2011; J. Aircraft 51(1):90–103, 2014, DOI 10.2514/1.C032045](https://ntrs.nasa.gov/api/citations/20110011597/downloads/20110011597.pdf)):
  - At constant indicated airspeed, true airspeed and μ rise with altitude while the speed of sound falls, giving "large changes to BVI noise levels" and more intense thickness noise.
  - They conclude "existing measurement-based empirically derived helicopter rotor noise source models may give incorrect noise estimates when they are used at conditions where data were not measured."
- AS350 measurements at three sites from sea level to 7000 ft ([Greenwood, Sim & Boyd 2016/2024](https://ntrs.nasa.gov/api/citations/20160009056/downloads/20160009056.pdf)):
  - At the same IAS and weight:
    - In-plane thickness noise rose about 5 dB over 7000 ft and scales with the effective hover tip Mach number M_He (Gopalan's definition).
    - Low-frequency loading noise fell by 6 dB in some directions as C_W rose 24%.
    - BVI fell 4 dBA in its peak direction.
    - Wake skew χ stays constant while μ changes.
  - Radiated noise "is shown to scale with ambient pressure when the flight condition ... is defined in nondimensional terms," and "no simple empirical model can be used to 'correct' the data."
  - FRAME predictions out of plane were "generally within 3 dB of measured levels."
- Watts, Greenwood & Stephenson give the companion test campaign ([Watts, Greenwood & Stephenson, "Measurement and Characterization of Helicopter Noise at Different Altitudes," AHS 72nd Forum, 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf)):
  - Amplitude differences up to 5.4 dB over 7000 ft of altitude at constant IAS and weight.
  - Amplitudes scale with ambient pressure when C_W and the effective hover tip Mach number are matched.
  - The campaign also established statistical run-to-run variability.
- Effects of weight and of true versus indicated airspeed on BVI in steady descent were studied separately ([Stephenson & Greenwood, AHS 71st Forum, 2015](https://aeroacoustics.psu.edu/publications/); listing only, not read).

**Wind and air-relative state**
- Greenwood & Sim modelled crosswinds on a Bell 430 ([Greenwood & Sim, "The Effects of Crosswind Flight on Rotor Harmonic Noise Radiation," AHS 69th Forum 2013; J. Aircraft 55(5):2137–2148, 2018, DOI 10.2514/1.C034634](https://ntrs.nasa.gov/api/citations/20130013862/downloads/20130013862.pdf)):
  - With the fuselage held along the inertial path (sideslip), TPP angle of attack fell about 2° and peak BVISPL was about 11 dB lower at 80 KIAS, −7.5°.
  - When "crabbed" for zero aerodynamic sideslip, "the BVI noise state is unchanged."
  - A headwind gradient (15 kt at 300 ft to 0 at the ground) cut peak BVI about 9 dB. The same gradient as a crosswind (crabbed) cut it only about 1 dB.
  - "Steady head or tail winds have no effect on the rotor's acoustic state," but gradients do.
  - In all cases, directivity rotates toward the wind relative to ground observers.
  - Recommendation: fly test tracks across the prevailing wind.

**Transients and the quasi-steady assumption**
- Sickenberger's first-principles longitudinal maneuver model (Bell 206B-3) "demonstrated that Quasi-Static Acoustic Mapping (Q-SAM) methods can be used to reliably predict the noise radiated during transient maneuvers" and that blade flapping "does not significantly affect the acoustics" ([Sickenberger, "Modeling Helicopter Near-Horizon Harmonic Noise due to Transient Maneuvers," Ph.D. dissertation, Univ. of Maryland, 2013 (title page 2012)](https://drum.lib.umd.edu/items/255f1df7-59bb-496b-a486-2e767895fb6a)):
  - Thickness and low-frequency loading contributions depend on TPP orientation relative to the observer; BVI depends on inflow and TPP orientation.
- FRAME-QS validation on a Bell 430 ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)):
  - An aggressive pitch-up (nearly 10°/s from 80 kt, 6° descent, peak 1.5 g) was observed about 10° below the horizon ahead and 30° either side.
  - With variable C_T, FRAME-QS captured the transient OASPL rise "within 2 dB in all three directions."
  - Holding C_T = C_W missed the large rise from load factor; the difference was ≈ 20 log10 n_z ≈ 3.5 dB.
  - A "Conventional" fixed-sphere model showed no maneuver increase at all.
  - The TPP held its orientation relative to the fuselage within 1° (optical TPP tracking).
- Greenwood, Schmitz & Sickenberger state that it "may not be possible to model the increase in BVI noise radiation due to transient maneuvering flight to the accuracy desired using the quasi-static (Q-SAM) method alone," because of the wake's dynamic motion. FRAME-D addresses this ([Greenwood, Schmitz & Sickenberger 2012](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf)).
- Stephenson et al.'s turning-flight program covers nine vehicles' data (ten tested, at three ranges) and twelve turning pairs ([Stephenson, Watts, Greenwood & Pascioni, "Development and Validation of Generic Maneuvering Flight Noise Abatement Guidance for Helicopters," VFS 76th Forum 2020; JAHS 67(1), 2022, DOI 10.4050/JAHS.67.012012](https://ntrs.nasa.gov/api/citations/20205003220/downloads/VFS76_Paper_v7.pdf)):
  - Peak band-limited (50–2500 Hz) A-weighted SPL was the BVI surrogate.
  - Turn-initiation altitude offsets were corrected. The authors derived six validated guidance principles.

**Multirotor / eVTOL**
- Non-unique trim: independently controlled propulsors allow different trim states for the same operating condition (Joby data), which motivated snapshot arrays ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- Rotor phase effects:
  - With equal RPM, "random phasing, which occurs in most small UAVs, can produce up to 20 dB uncertainty in the tonal sound pressure level at a given ground observer" ([Pascioni & Rizzi, "Tonal Noise Prediction of a Distributed Propulsion Unmanned Aerial Vehicle," AIAA paper, c. 2018 (NTRS 20190000886; venue not confirmed in text)](https://ntrs.nasa.gov/api/citations/20190000886/downloads/20190000886.pdf)).
  - A 6 dB BPF reduction from phase control requires the rotation-rate deviation to stay within about 0.5% ([Pascioni, Rizzi & Schiller 2019](https://ntrs.nasa.gov/api/citations/20200002397/downloads/20200002397.pdf)).
- Penn State work on eVTOL flight-state variance, rotor phase synchronization hardware, and wind effects on multirotor noise exists, but only listings were seen: Valente, Greenwood & Johnson, "An Experimental Investigation of eVTOL Flight State Variance on Noise," VFS 79th Forum, 2023 ([ASCENT listing](https://ascent.aero/publication/an-experimental-investigation-of-evtol-flight-state-variance-on-noise/)); Valente, Greenwood & Johnson, rotor phase synchronization, VFS 80th Forum, 2024; ElSharkawy et al., QuietDrones 2024 ([PSU list](https://aeroacoustics.psu.edu/publications/)).
- The UAM noise gap assessment notes that AAM "permit[s] the user to define explicitly an appropriate noise sphere for each trajectory point" and that automatic sphere selection "could be improved to allow for UAM specific operational modes" ([Rizzi, Huff, Boyd et al., "Urban Air Mobility Noise: Current Practice, Gaps, and Recommendations," NASA/TP-2020-5007433, Oct 2020](https://ntrs.nasa.gov/api/citations/20205007433/downloads/NASA-TP-2020-5007433.pdf)).

### Inferences
- A hemisphere database indexed by (V, γ) alone is the RNM baseline and loses load factor and ambient effects. The literature's progression is:
  - (V, γ) → (V, γ_e) for acceleration (Q-SAM).
  - → (V, γ_e, n_z) or C_T for maneuvers (FRAME-QS).
  - → non-dimensional (μ, C_T/C_W, M_H or effective M_tip, λ or γ_e) for transfer across altitude and temperature (Greenwood & Schmitz 2011/2014; Greenwood, Sim & Boyd 2016/2024).
- Indexing by (μ, γ_e, C_T), as NICE-OPS does, matches the published recommendation. The source in which these spheres are documented (FRAME-QS) indexes by (V, γ, n_z).
- The 3.5 dB result says a footprint tool that evaluates C_T from load factor during maneuvers recovers most of the transient rise. Holding C_T fixed misses it entirely.
- Crosswind results imply the sphere should be selected and oriented by the air-relative velocity (crabbed flight leaves the acoustic state unchanged; directivity rotates with the wind). Wind *gradients* change trim and are not captured by any steady-wind correction.
- For multirotor eVTOL, a single dB sphere per condition is a lossy representation whenever rotor phase is uncontrolled. Variability of up to 20 dB tonal at a point is possible in principle, and averaging-type spheres smooth it out.

### Gaps
- The original Q-SAM journal papers by Gopalan and Schmitz were not read. One JAHS 52(1), 2007 paper on quasi-static flight path control of BVI in steady wind was located on [Ingenta](https://www.ingentaconnect.com/content/ahs/jahs/2007/00000052/00000001/art00002) and [ResearchGate](https://www.researchgate.net/publication/239417941_Quasi-Static_Flight_Path_Control_of_Helicopter_Blade_Vortex_Interaction_Noise_in_a_Steady_Wind), but its author list was not verified.
- No quantitative eVTOL "hemisphere per condition" validation (dB error) was found. The Valente et al. 2023 variance numbers were not accessed.
- No published wake-skew-based sphere indexing scheme was found, beyond FRAME's use of inflow.

## 4. Interpolation across operating conditions and within the sphere; spherical harmonics; tonal vs broadband; amplitude vs dB

### Takeaway
Across conditions, published schemes range from nearest-neighbour (RNM in flight-path angle; AAM NCSPEC holds a sphere constant between points), through linear in pressure-squared in one variable (RNM in airspeed), to natural-neighbour interpolation over a 3-D Delaunay triangulation (FRAME-QS in V, γ, n_z). Within a sphere, scattered measured points are gridded by Laplace/spline blends (ART), Coons surfaces (RNM), or modified quadratic Shepard on pressure-squared (NASA snapshot work). Regression onto Taylor, Fourier or spherical-harmonic bases is an alternative. Coherent multirotor work interpolates complex amplitude and phase.

### Cited Findings
- RNM's two-step lookup ([Lucas & Marcolini, c. 1997](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)):
  - First, "a four-way interpolation for theta and phi ... using a 'Coons' surface interpolation technique."
  - Second, interpolation between two hemispheres' airspeeds "assumes a linear variation in pressure squared with airspeed," which the user can override to take the nearest sphere.
  - Hemispheres are first selected as those "with flight path and pylon angles that are nearest to the specified conditions" (nearest-neighbour in those variables).
- FRAME-QS: "Spheres for intermediate conditions not in the database are computed as needed using natural neighbor interpolation over the three-dimensional (V, γ, n_z) Delaunay triangulation" of ≈500 spheres ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- AAM NCSPEC mode "utilizes a constant noise source between simulation points" and "does not allow interpolation between different spheres." Users insert guard points to stop AAM interpolating speed ([Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)).
- Snapshot hemisphere gridding: "The interpolation is done on the pressure-squared acoustic data for each frequency band using an inverse distance weighting scheme, described as Modified Quadratic Shepard" (Franke & Nielson 1991), with a 30° radius of influence. The authors note that both the Shepard choice and the 30° radius are assumptions ("arbitrary," though enlarging the radius had no noticeable effect) ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- ART grids scattered emission-angle data with a user-chosen blend of Laplace and spline interpolation ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- Regression bases ([Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)):
  - Taylor (polynomial) fits are better at low order; Fourier fits follow rapid changes at high order.
  - An energy-based least-squares fit produced "erroneous results" at 63 Hz under strong directivity, because minimising absolute error over the whole hemisphere lets local (relative) error grow where levels are low.
- Spherical harmonics: Mobley's dissertation extends EMPA's spherical-harmonic helicopter directivity models, adding removal of contaminated data, averaging of harmonic representations, and "methods for interpolating between spherical harmonic source representations" of different conditions. It is validated with dipole tests and UH-1 data from White Sands ([Mobley, "Further development of 3-D rotary-wing acoustic directivity using a spherical harmonic representation," Ph.D. dissertation, Univ. of Dayton, 2012](https://ecommons.udayton.edu/graduate_theses/471); abstract only).
- Coherent sources: complex pressure magnitude and phase per propeller are interpolated on each hemisphere ("Interpolated pressure amplitude ... Interpolated phase on source hemisphere") before phase-correct summation ([Pascioni, Rizzi & Schiller 2019](https://ntrs.nasa.gov/api/citations/20200002397/downloads/20200002397.pdf)).
- Tonal vs broadband:
  - RNM carries separate pure-tone (narrowband) and broadband (1/3-octave) hemispheres per condition ([Lucas & Marcolini](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf); [Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
  - Main rotor harmonic and broadband have different directivity, and BPF-band sphere interpolation errors exceed broadband ones ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).

### Inferences
- The measured-data community interpolates energy (p²), not dB (Houston 2025 within-sphere; RNM across airspeed). Interpolating dB is a geometric mean of energies and biases low between dissimilar spheres. A tool that blends dBA values (as NICE-OPS's grid does) will differ from a p² blend wherever neighbouring spheres differ by several dB. The size of that difference was not found quantified in the literature.
- Delaunay/barycentric or natural-neighbour interpolation across conditions has precedent in FRAME-QS. RNM's nearest-in-γ approach produces step changes at segment boundaries that a footprint will show as artifacts.
- Interpolating between spheres of different flight speed also blends different Doppler shifts and different convective amplification when spheres keep the received spectrum (section 6).

### Gaps
- No published study was found that compares dB-domain against p²-domain interpolation across operating conditions, or barycentric against natural-neighbour, with footprint dB statistics.
- Full text of the EMPA spherical-harmonic papers and Mobley's accuracy numbers were not accessed.

## 5. Sphere frame/orientation conventions, source centre, multiple-source decomposition

### Takeaway
Conventions vary:
- RNM: body-axis spheres at the CG.
- FRAME and the Maryland group: TPP-fixed spheres, or horizon-fixed spheres for steady data.
- NASA snapshot: main-rotor hub.
- ANOPP2/AARON: CG.
- Multirotor auralization: one sphere per rotor, centred on that rotor.

Orientation errors matter for directive (BVI-rich) vehicles, and measured roll, pitch and heading or direct TPP measurement is recommended. Main/tail rotor separation into separate hemispheres is demonstrated and needed for integrated metrics.

### Cited Findings
- RNM sphere coordinates: origin at the noise source CG, Z-axis along the fuselage. θ runs from 0 at the bow to 180° at the stern; φ is 0 directly beneath, +90° port, −90° starboard. The netCDF format carries a "Radius" and an "XYZ" sphere centre ([Lucas & Marcolini](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)).
- Orientation in ART and its sensitivity ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)):
  - ART "assumes that the vehicle fuselage orientation is the same as the propulsor orientation." Without attitude data it estimates orientation from the track, with roll assumed zero.
  - MD600N (benign directivity): measured versus estimated orientation changed the footprint by within ±1 dBA.
  - XV-15 at 110 kt, 80° nacelle, 6° approach (intense BVI): the difference was much larger, with pitch dominant. The delta-SEL contour scale in their Fig. 24 spans several dB; exact values could not be extracted from the text layer.
  - The authors call for "direct measurement of rotor tip-path-plane orientation" via flapping.
- In FRAME-QS the steady spheres are rotated by the computed change in the inertial TPP normal n̂_a. Pitch and roll changes are derived from v and n̂_a ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- Upper and lower hemispheres differ: FRAME OASPL spheres (100 kt level) show different magnitude and directivity above and below the horizon, "consistent with the findings of Sickenberger," so "models which assume noise radiation to be symmetric about the rotor tip-path-plane may lead to inaccurate" results ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- Horizon-fixed hemispheres were used to present steady turning-flight BVI directivity ([Schmitz et al. 2007](https://ntrs.nasa.gov/api/citations/20070021458/downloads/20070021458.pdf)).
- Centre choices in published work:
  - NASA snapshot: main rotor hub, 30.48 m ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
  - AARON: aircraft CG ([Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)).
  - Per-rotor hemispheres centred on each rotor, with pitch and roll displacement accounted for ([Krishnamurthy, Rizzi & Cheng 2021](https://ntrs.nasa.gov/api/citations/20205010694/downloads/20205010694_ST21_Paper.pdf)).
  - Per-propeller at 12R ([Pascioni, Rizzi & Schiller 2019](https://ntrs.nasa.gov/api/citations/20200002397/downloads/20200002397.pdf)).
  - A 30.48 m reference distance, 2.86 main-rotor diameters ([Pera, Rizzi, Krishnamurthy, Fuller & Christian, "A Recording-Based Method for Auralization of Rotorcraft Flyover Noise," AIAA paper, 2018](https://ntrs.nasa.gov/api/citations/20180002090/downloads/20180002090.pdf)).
  - FRAME-QS validation levels normalised to 10 rotor radii from the hub ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- Multiple sources:
  - Main and tail rotor separation with separate directivities ([Greenwood & Schmitz 2014](https://pure.psu.edu/en/publications/separation-of-main-and-tail-rotor-noise-from-ground-based-acousti/); [Rachaprolu & Greenwood 2024](https://rosap.ntl.bts.gov/view/dot/72623/dot_72623_DS1.pdf)).
  - Tail rotor content is essential for PNLT and EPNL ([Boyd et al. 2017](https://ntrs.nasa.gov/citations/20170002590)).
  - RNM can hold hemispheres "centered at different points on the aircraft" (e.g. main and/or tail rotor) ([Lucas & Marcolini](https://ntrs.nasa.gov/api/citations/20040110391/downloads/20040110391.pdf)).

### Inferences
- A tool that mirrors or symmetrises a lower-hemisphere sphere about the TPP to fill the upper half contradicts published FRAME/Sickenberger results showing the halves differ. Mirroring is at best a placeholder for coverage.
- Sphere centre and radius only matter if observers come close (within a few sphere radii) or sources are spatially separated (multirotor, tail rotor at about 1 diameter aft). At typical footprint ranges a CG/hub offset is negligible geometrically. But a radius smaller than the source extent makes the compact-source assumption invalid near the vehicle.
- For hover and low speed, the published orientation conventions (fuselage or TPP) require a heading. A velocity-aligned frame is undefined there. The literature consistently orients by attitude or TPP rather than by velocity.

### Gaps
- No direct quantitative study was found of TPP-frame versus body-frame versus velocity-frame sphere orientation error on footprints, beyond Conner et al.'s measured-versus-estimated attitude comparison.

## 6. Doppler: emission-referenced vs reception-referenced spheres

### Takeaway
Two conventions coexist:
- **Emission-referenced spheres:** the frequency shift is removed, and the propagation code reapplies Doppler for the trajectory. This is the convention of RNM/ART, ANOPP2-AnoppGround, de-Dopplerized virtual-microphone methods, and predicted spheres.
- **Reception-referenced spheres:** the measured (shifted) spectrum is kept for that emission direction and flight speed. This is valid only if the propagation code then applies no further shift.

Mixing them double-applies or omits the shift.

### Cited Findings
- RNM propagation corrects for "spherical spreading, atmospheric absorption, ground reflection and attenuation, Doppler shifts, and the difference in phase between the direct and reflected rays." ART uses the same physics in reverse to build spheres ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- ANOPP2-AnoppGround flies the predicted hemisphere and accounts for "hemisphere motion ('Doppler effect')" ([Boyd et al. 2017](https://ntrs.nasa.gov/citations/20170002590)).
- Time-domain de-Dopplerization transforms ground signals to an observer "moving with the aircraft and removing the Doppler frequency shift from the signal." Amplitudes can be normalised to the virtual microphone using the known emission-time distance ([Rachaprolu & Greenwood 2024](https://rosap.ntl.bts.gov/view/dot/72623/dot_72623_DS1.pdf); [Greenwood & Schmitz 2014](https://pure.psu.edu/en/publications/separation-of-main-and-tail-rotor-noise-from-ground-based-acousti/)).
- Recording-based auralization de-Dopplerizes ground recordings to a virtual microphone at 30.48 m "moving with the source" before segmenting by emission angle ([Pera et al. 2018](https://ntrs.nasa.gov/api/citations/20180002090/downloads/20180002090.pdf)).
- The snapshot de-propagation equation (spreading, absorption, ground board) contains no Doppler term. Its 0.5 s spectra are formed at the received frequencies ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- Project-internal context (not literature; [local memory note](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/sphere-doppler-convention.md)): panam-built `.nod`/sphere files carry `DOPPLER_SHIFT_REMOVED = 0`. They are reception-referenced, so NICE-OPS applies only the reception-time energy factor and no frequency shift.

### Inferences
- Snapshot and Houston-style spectral hemispheres, like panam's, appear to keep the received spectrum. RNM/ART, ANOPP2 and the de-Dopplerized virtual-microphone methods remove it. A database should carry an explicit flag, and a footprint tool must branch on it.
- A reception-referenced sphere bakes in one flight speed's shift per direction. Interpolating spheres across airspeed blends shifts, and the change in A-weighting and absorption with the shift is approximated. For broadband dBA this is small; for tonal or band-limited metrics near the nose or tail at high speed it is larger. No published quantification was found.
- Separately from frequency, the convective (1 − M cos θ) amplitude and time-compression factors are physically part of the received field. A de-Dopplerized emission sphere still contains convective amplification for its own speed, since it was measured from a moving source. The reception-time energy factor (dt_obs/dτ) is a propagation-side term in either convention.

### Gaps
- No published paper was found that explicitly documents the AAM/RNM netCDF "Doppler removed" attribute or compares footprints built from emission- versus reception-referenced spheres. The AAM Technical Reference (Page et al., Volpe, 2020, cited by Houston et al.) would be the place to confirm, and was not accessed.

## 7. Known shortcomings and validation of hemisphere / Q-SAM / FRAME approaches against measurements

### Takeaway
Steady-condition hemispheres reproduce steady flyovers to about 1–3 dB. They fail when:
- the state is not steady or unique (maneuvers with load factor, multirotor non-unique trim, random rotor phase);
- the conditions differ from the test conditions (altitude/temperature at constant IAS, wind gradients);
- the coverage is missing (near and above the horizon);
- the source assumptions break (main rotor only, no broadband, symmetric about the TPP).

Quasi-static mapping with load factor/C_T and non-dimensional indexing recovers transients to within about 2 dB in the one published FRAME-QS comparison.

### Cited Findings
- FRAME-QS (variable C_T) captured an aggressive pitch-up's transient OASPL rise "within 2 dB in all three directions." A constant-sphere model missed it entirely, and a constant-C_T model missed about 3.5 dB ([Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)).
- FRAME across ambient conditions: out-of-plane predictions "generally within 3 dB" of measured. RNM/ART hemispheres at the same dimensional condition differ markedly between sea level and 7000 ft ([Greenwood, Sim & Boyd 2016/2024](https://ntrs.nasa.gov/api/citations/20160009056/downloads/20160009056.pdf)).
- Conventional empirical models "may give incorrect noise estimates when they are used at conditions where data were not measured" ([Greenwood & Schmitz 2011/2014](https://ntrs.nasa.gov/api/citations/20110011597/downloads/20110011597.pdf)). Altitude alone changed levels by up to 5.4 dB at constant IAS and weight ([Watts, Greenwood & Stephenson 2016](https://ntrs.nasa.gov/api/citations/20160009057/downloads/20160009057.pdf)).
- Wind: headwind gradients changed peak BVI by about 9 dB and sideslip by about 11 dB ([Greenwood & Sim 2013/2018](https://ntrs.nasa.gov/api/citations/20130013862/downloads/20130013862.pdf)).
- Quasi-static BVI may miss dynamic wake motion in aggressive transients. FRAME-D roll predictions were too high on the advancing side ([Greenwood, Schmitz & Sickenberger 2012](https://ntrs.nasa.gov/api/citations/20120008939/downloads/20120008939.pdf)). Sickenberger nonetheless found Q-SAM reliable for the longitudinal transients studied ([Sickenberger 2013](https://drum.lib.umd.edu/items/255f1df7-59bb-496b-a486-2e767895fb6a)).
- Sphere interpolation error from array sparsity reaches more than 14 dB at points with few nearby microphones. Medians are 0.3–1.8 dB depending on density, and worst at the rear near the horizon ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- Linear-array hemispheres are "unreliable or difficult to interpret" for UAM/FVL vehicles that do not hold a steady aerodynamic state ([Houston et al. 2025](https://ntrs.nasa.gov/api/citations/20250002951/downloads/VFS81_Snapshot_Array_04142025_v2.pdf)).
- Uncontrolled rotor phase allows up to 20 dB tonal uncertainty at a point ([Pascioni & Rizzi](https://ntrs.nasa.gov/api/citations/20190000886/downloads/20190000886.pdf)).
- Simplified directivity: omnidirectional spheres give a 1.7 dB SEL_A error under-track, and lateral errors grow without azimuthal directivity ([Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)).
- Prediction-based hemispheres: main rotor only "severely" underpredicts integrated metrics, and the descent A-weighted metrics were 4–5% low without extra broadband ([Boyd et al. 2017](https://ntrs.nasa.gov/citations/20170002590)).
- Orientation method: up to ±1 dBA for a benign vehicle, and substantially more for a BVI-rich tiltrotor ([Conner et al. 2006](https://ntrs.nasa.gov/api/citations/20060020242/downloads/20060020242.pdf)).
- European practice: an NLR-led method for European land-use planning uses measured hemispheres "for a range of conditions within the flight envelope" and clusters the fleet into representative types ([Tuinstra, Stevens, van Oosten & Olsen, "Towards a European helicopter noise calculation method," 44th ERF, paper 46, 2018](https://dspace-erf.nlr.nl/items/22c30b45-c166-4027-ac73-3395d624262b); abstract only).
- Airbus HELENA uses hemisphere databases and has validation papers for EC130 and EC135 ([HELENA validation paper listing](https://www.researchgate.net/publication/289803872_Validation_of_EC130_and_EC135_environmental_impact_assessment_using_HELENA)); not read, no dB figures obtained.

### Inferences
What a hemisphere-based footprint tool should do, by convergence across these sources:
1. Index spheres by non-dimensional state (μ, C_T or C_W, γ_e, M_tip) rather than IAS and weight.
2. Compute γ_e and C_T from the trajectory's acceleration and load factor (Q-SAM/FRAME-QS).
3. Orient by the TPP or attitude, using the air-relative velocity.
4. Interpolate in energy.
5. Flag unmeasured cells rather than extrapolating.
6. Keep tonal and broadband (and main and tail rotor) separable.
7. Know whether the spheres are Doppler-removed.

Where the concept breaks down, by convergence:
- aggressive transients with wake dynamics;
- multirotor phase and non-unique trim;
- near-horizon and upper-hemisphere radiation;
- wind gradients during descent;
- any operation outside the measured envelope when the database is dimensional.

### Gaps
- No published footprint-level (SEL contour) validation of a hemisphere tool against a dense measured ground grid with dB error statistics was found in this session, beyond Conner et al.'s repeat-run spread and the point-microphone comparisons above. HELENA and European-method validation numbers were not accessed.
- No published measurement of how much hover or low-speed heading versus drift orientation changes a footprint was found.
