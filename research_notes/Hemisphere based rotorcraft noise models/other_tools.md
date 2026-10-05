# Hemisphere / directivity-sphere rotorcraft noise tools other than AAM and NORAH2

Scope note for the report writer:

- **Covered:**
  - RNM, with its RNM-QS, FRAME-QS and APET descendants.
  - FAA HNM and the INM/AEDT helicopter mode, plus the ECAC Doc 29 / ICAO Doc 9911 segmentation approach, as a contrast.
  - NASA ANOPP2/AMAT and the NASA Auralization Framework.
  - Empa sonAIR, FLULA2 and sonX.
  - DLR SELENE, and the DLR/Eurocopter "few aerodynamic parameters" database work.
  - HELENA (FRIENDCOPTER; used in Clean Sky PhoeniX).
  - Airbus Helicopters CAROT / MOTUS.
  - The Clean Sky MANOEUVRES quasi-steady work (Politecnico di Milano, Roma Tre, AgustaWestland).
  - NLR hemisphere measurement.
  - JAXA and Penn State.
  - TU Delft eVTOL work (UYGUR Gaussian beams).
  - Korean UAM hemisphere work (KARI, Hanseo, Chungbuk).
- **Contrast only:** AAM and NORAH2 are covered by other researchers.
- **Not retrieved:** several primary documents returned 403 or WAF blocks: Empa, Volpe, ROSA-P, Acta Acustica, and the INM/AEDT technical manuals.
- **Marking:** **[snippet only]** means a fact rests on a search-engine summary, not on a document I read.
- **Search quota:** it ran out before the HNM, AEDT technical-manual and sonAIR-helicopter details could be closed (see Gaps).

## 1. Source representation and how it depends on operating state (and interpolation between states)

### Takeaway
Almost every research-grade rotorcraft footprint tool uses the same core idea as NICE-OPS. That includes RNM, SELENE, HELENA, CAROT, FRAME-QS, AMAT, MANOEUVRES, JAXA and the Korean UAM tools. Each stores a fixed-radius hemisphere or sphere of 1/3-octave levels per steady operating condition (some add narrowband or tonal data). A small set of flight parameters selects or interpolates it.

They differ mainly in which parameters index the database:

- **RNM / AAM:** (airspeed, flight path angle, nacelle angle), interpolated linearly in energy.
- **SELENE:** (mu, alpha_TPP, CT) on a Delaunay triangulation, with nearest-neighbour outside the envelope. This is the closest published analogue to NICE-OPS's (mu, gamma, CT) barycentric scheme.
- **FRAME-QS:** (V, gamma, n_z), natural-neighbour interpolation over a Delaunay triangulation.
- **HNM / INM / AEDT (regulatory):** noise-distance curves at three lateral directions only (centreline and +/-45 deg), for a handful of fixed modes.

### Cited Findings

**RNM (NASA Langley / Wyle; Conner & Page)**
- RNM's vehicle sources are "sound hemispheres". They can come from theory, wind tunnel, flight test or a combination. They hold broadband 1/3-octave levels or pure tones with phase. Points on a hemisphere are "a fixed radius and two spherical angles", stored as netCDF. — [Conner & Page, "A Tool for Low Noise Procedures Design and Community Noise Impact Assessment: The Rotorcraft Noise Model (RNM)", Heli Japan 2002, NTRS 20030003730](https://ntrs.nasa.gov/citations/20030003730)
- Each hemisphere holds one quasi-steady flight condition with three independent variables: airspeed, flight path angle, and nacelle pylon angle (90 deg for a conventional helicopter). A vehicle may carry up to ten broadband and ten pure-tone hemispheres at different points on the aircraft (engine, main rotor, tail rotor, MR/TR wake interaction). — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- "Hemisphere selection" is a hierarchical search over nacelle angle, then flight path angle, then speed. With enough hemispheres, RNM interpolates **linearly in the energy domain (pressure squared)**, first in airspeed and then in flight path angle. With sparse hemispheres it takes the closest flight condition without interpolating. — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- Hemispheres come from flight test via the Acoustic Repropagation Technique (ART2), shipped with RNM. It needs synchronized steady flight over a microphone array; the data are back-propagated to a fixed-radius hemisphere. — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- The emission direction is where the straight source-receiver line crosses the hemisphere. The v1-era geometry "assumes a zero roll angle and makes small angle approximations for all angles except yaw and heading" when computing the hemisphere Euler angles. — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM version 7 represents each source "by a sphere of spectral data at a reference distance". It accepts broadband 1/3-octave, narrow band (user-defined bandwidth and spacing) or pure-tone data. — [NASA Tech Briefs, RNM (LAR-17753-1), May 2016](https://www.techbriefs.com/component/content/article/24556-lar-17753-1); [NASA software catalog LAR-17753-1](https://software.nasa.gov/software/LAR-17753-1)

**RNM-QS / Q-SAM and FRAME-QS (UMD, NASA Langley)**
- Q-SAM (Gopalan & Schmitz, UMD) makes a first-order equivalence between steady and accelerating or decelerating flight with the same tip-path-plane angle of attack. Integrated into RNM as RNM-QS, it extends steady straight-line spheres to moderate longitudinal accelerations. **[snippet only]** — [Gopalan, "Quasi-Static Acoustic Mapping of Helicopter Blade-Vortex Interaction Noise", UMD thesis](https://drum.lib.umd.edu/items/cb681d36-4777-4d14-80aa-b8a1dace6f1e/full)
- FRAME-QS (Greenwood, Rau, May, Hobbs, AHS 71, 2015) applies Q-SAM to a FRAME-generated database of "nearly 500 predicted acoustic spheres". The database covers 40-130 kt, flight path angles +3 to -18 deg and load factors 0.75-2.0. Intermediate conditions use **natural-neighbour interpolation over the 3-D (V, gamma, n_z) Delaunay triangulation**. — [Greenwood et al., "A Maneuvering Flight Noise Model for Helicopter Mission Planning", NTRS 20160005979](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)
- The same paper defines an effective flight path angle, gamma_e. It is the angle between the spline-derived inertial velocity and the "inertial tip-path-plane" normal tilted by inertial acceleration. Maneuvers thereby map onto an equivalent steady TPP angle of attack. The paper notes that RNM and AAM "will not compute the change in orientation of the steady-state flight acoustic spheres due to maneuvering flight"; the user may specify orientation manually. — [Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)
- In the FRAME approach, spheres are generated over advance ratio mu, thrust coefficient CT and TPP angle of attack. Flight-test hemisphere databases typically vary airspeed and flight path angle. — [Greenwood & Sim, "The Effects of Crosswind Flight on Rotor Harmonic Noise Radiation", NTRS 20130013862](https://ntrs.nasa.gov/api/citations/20130013862/downloads/20130013862.pdf)

**FAA HNM / INM / AEDT helicopter mode (contrast: noise-distance curves, not hemispheres)**
- In helicopter mode, AEDT uses manufacturers' "noise-operating condition-distance curves", with the operating condition set directly by procedural step. The reason is that helicopter mode has no performance model and no correlating parameter (fixed-wing mode uses corrected net thrust). — [Letica & Rizzi, "On the modeling of urban air mobility vehicle takeoff and landing operations in the FAA AEDT", NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- Helicopter NPD data for modes L (level), D (departure) and A (approach) are given at sideline elevation angles phi = 0 and +/-45 deg. The lateral directivity adjustment LD_ADJ interpolates between them and uses the +/-45 deg data beyond. Other modes are substituted by AEDT rules (e.g. accelerating departure Mode E becomes Mode D). This "severely limits the number of unique operational states". — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- AEDT helicopter exposure is E = E_NPD + NF_ADJ + DUR_ADJ - LA_ADJ + MN_ADJ + LD_ADJ, where:
  - MN_ADJ is the advancing-tip-Mach-number adjustment.
  - NF_ADJ, the noise fraction, uses a fourth-power 90 deg dipole model.
  - DUR_ADJ is a log ratio of segment speed to the mode's reference speed.
  - — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- The INM 7.0 Technical Manual added helicopter adjustments for lateral directivity, static directivity and static duration. **[snippet only]** — [INM 7.0 Technical Manual, ROSA-P dot/12188](https://rosap.ntl.bts.gov/view/dot/12188)
- HNM v2.2 is Volpe's model for the FAA "for determining the impact of helicopter noise in the vicinity of terminal operations". **[snippet only; document not opened]** — [HNM Version 2.2 User's Guide, ROSA-P dot/8828](https://rosap.ntl.bts.gov/view/dot/8828)
- ACRP (Wyle, Volpe, NLR, KB Environmental) found that improving AEDT for helicopters should focus first on source modelling: spectral content, lateral directivity and operational sensitivity. It recommended lateral source characteristics with enough angular resolution to capture directional BVI. **[snippet only; report pages are image-only online]** — [ACRP, "Recommended Community Noise Model Enhancements to Improve Prediction of Helicopter Activity Impacts", 2015, doi:10.17226/22079](https://nap.nationalacademies.org/read/22079)
- ICAO Doc 9911, ECAC Doc 29 and Directive 2015/996 are "limited to fixed-wing aircraft, thus not covering rotorcraft noise". — [EASA, Research Project NOISE final report](https://www.easa.europa.eu/en/downloads/139729/en)
- ECAC Doc 29 5th edition (endorsed by DGCA/166, 19 May 2026): the ANP database excludes "data for military aircraft and helicopters". — [ECAC Doc 29 5th Ed., Vol. 2 Technical Guide](https://ecac-ceac.org/images/documents/ECAC-CEAC-DOC_29_5th_Edition-REPORT_ON_STANDARD_METHOD_OF_COMPUTING_NOISE_CONTOURS_AROUND_CIVIL_AIRPORTS-Volume_2-TECHNICAL_GUIDE.pdf)
- Contrast, the NORAH guidance: hemispheres are binned at 10 deg in azimuth and polar angle, explicitly unlike "the Noise Power Distance methods for fixed wing". — [EASA NOISE D1.5d Rotorcraft noise modelling guidance, Jan 2024](https://www.easa.europa.eu/en/downloads/132005/en)

**NASA ANOPP2 / AMAT / AARON (UAM)**
- AMAT (ANOPP2 Mission Analysis Tool) is a time-marching simulation that flies arbitrary paths with source noise hemispheres computed by AARON, ANOPP2's Aeroacoustic ROtor Noise tool. AMAT also generated the NPD data used as AEDT inputs. — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- For the NASA RVLT quadrotor, AARON (Farassat F1A, periodic loading and thickness noise) computed pressure time histories on a sphere of observers at a fixed radius from the CG. Only the lower hemisphere was computed. Each sphere, one per operating state, went to netCDF in 1/3- and 1/12-octave bands for AAM (the "Gen 1.2" database). — [Rizzi, Page & Cheng, InterNoise 2021, NTRS 20210015134](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)
- The same study derived omnidirectional and axisymmetric spheres from the full 3-D sphere and calibrated all of them (+0.21 dB) to a 90 kt / 1000 ft reference. The AEDT NPDs came from the 0 deg-azimuth slice only, under a fixed-point fixed-wing profile. — [Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)
- ANOPP2 is "NASA's aeroacoustic toolset", built to combine acoustic approaches of varying fidelity for conventional and unconventional aircraft. — [Rizzi, Lopes & Burley, "NASA's Aeroacoustic Tools and Methods for Analysis of Aircraft Noise", NTRS 20160006703](https://ntrs.nasa.gov/citations/20160006703)
- The ANOPP2 tool names I found are AARON, AMAT, AANNT (artificial neural network tool) and the AFFIFMs (Farassat Formulations Internal Functional Modules). The names "AFFAL", "ASNIM" and "ARNM" in the brief were **not found** in any source. — search summary over [NTRS 20160014858](https://ntrs.nasa.gov/citations/20160014858) and [ANOPP2 v1.5.0 software page](https://software.nasa.gov/software/LAR-18567-1) **[snippet only]**

**DLR SELENE (and the DLR/Eurocopter "few aerodynamic parameters" database)**
- SELENE = "Sound Exposure Level starting from Evaluation of Noise Emissions". Its database of fixed-radius "noise directivity hemispheres" came from dedicated EC135 ACT/FHS and BO105 flight experiments. Hemispheres are **stored by main-rotor advance ratio mu, tip-path-plane angle alpha_TPP and thrust coefficient CT**. "An unstructured mesh is then generated by use of a Delaunay triangulation" to interpolate unmeasured conditions; "outside the measured flight envelope, a nearest neighbour search is used". — [Olsman & Lummer, "Influence of wind on the noise footprint of a helicopter landing", ERF (c. 2010-12; year not confirmed)](https://dspace-erf.nlr.nl/bitstreams/61aeceb7-6f3b-486f-8d7e-d335cd4db426/download)
- The stated rationale: alpha_TPP relates to vortex miss distance, mu to the epicycloidal shape of the tip-vortex paths, and CT to blade loading. Together they govern main-rotor BVI in steady and quasi-steady flight. The aircraft's attitude angles orient the interpolated hemisphere in space. — [Olsman & Lummer](https://dspace-erf.nlr.nl/bitstreams/61aeceb7-6f3b-486f-8d7e-d335cd4db426/download)
- DLR describes SELENE as using "a database of noise hemispheres for the description of noise directivity and level as a function of the flight condition", from flight experiments or numerical simulation. It is licence-only. — [DLR SELENE software page](https://www.dlr.de/en/as/research-and-transfer/software-solutions/acoustics/software-selene)
- Le Duc, Spiegel, Guntzer, Lummer & Götz (ODAS 2008) modelled ground footprints in arbitrary maneuvers using a 2004 flight-test aeroacoustic database. Their premise is that "the instantaneous noise emission of the whole helicopter can be characterized by the instantaneous value of few aerodynamic parameters". — [elib.dlr.de/58498](https://elib.dlr.de/58498)
- DLR's EC135 work (PAVE and FRIENDCOPTER):
  - Noise was found to be "mainly dependant on the main rotor advance ratio RmMu, its thrust coefficient RmCT and its tip path plane angle of attack PmAlpha".
  - Directivity spheres were back-propagated from up to 43 microphones on an 800 m disk to a 300 m-radius sphere.
  - Engine torque tracks BVI, and Fenestron noise at low torque can be removed with sideslip.
  - — [Spiegel, Guntzer, Le Duc & Buchholz, ERF 34, 2008](https://dspace-erf.nlr.nl/bitstreams/16abfa07-16f9-4711-adbf-2a87e3cbd01b/download)

**HELENA (FRIENDCOPTER; used in Clean Sky PhoeniX)**
- HELENA (Helicopter Environmental Noise Analysis) was developed in the EU FRIENDCOPTER project and computes footprints from experimental or numerical (CFD) noise databases. It was the noise module of Clean Sky's PhoeniX platform, alongside the EUROPA performance code and the GSP engine code. EUROPA supplies each segment's position, altitude and tip-path-plane angles over time. — [Smith et al., "Achieving rotorcraft noise and emissions reduction for 'Clean Sky'", ERF 2014](https://dspace-erf.nlr.nl/bitstreams/c920dfb8-f809-440a-9387-af4cf81fcdd9/download)
- A FRIENDCOPTER data dictionary, written mainly at DLR (22 coordinate systems, 400 variables), was adopted "for collective development of the HELENA code". — [Spiegel et al. 2008](https://dspace-erf.nlr.nl/bitstreams/16abfa07-16f9-4711-adbf-2a87e3cbd01b/download)
- HELENA was validated for EC130 and EC135 environmental impact assessment (Gervais, Gareton, Dummel, Heger, AHS 66, 2010; cited, not read). — [reference list in Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)
- HELENA was reportedly developed and tested for the EC135 and "showed improvement over NPD tables". **[snippet only; source document not identified]**

**Airbus Helicopters CAROT / MOTUS**
- CAROT (Compute the Acoustics of a Rotorcraft Over Terrain) is Airbus Helicopters' in-house propagation software, used "on a daily basis". Its hemisphere databases come from field measurement. The H130 database, characterised in FRIENDCOPTER, gives mean A-weighted hemisphere level against steady flight condition: a BVI rise at negative vertical speed, and a rise with speed in level flight. — [Guntzer, Caillet, Cariou, Pinacho, Dieumegard, Roca León, "A Comprehensive Helicopter Acoustic Modeling Tool based on Simulation and Experiment", ERF 2023 (VFS 2023 Cheeseman best paper)](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)
- MOTUS:
  - The ROSI time-domain tool separates the signal into main-rotor tonal, anti-torque tonal and residual hemispheres. Engine noise can also be split out when bench data exist, giving four.
  - Tonal hemispheres are stored **de-dopplerized**, in a "wind-tunnel" frame with a fixed frequency vector.
  - Design changes are applied as "masks" at three fidelity levels: per direction and frequency, a sphere-mean spectrum, or a constant shift.
  - Masks come from simulation (HOST/CAMRAD II, free wake, FLOWer CFD, Kirchhoff/FW-H) or from experiment.
  - — [Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)

**Clean Sky MANOEUVRES (PoliMi, Roma Tre, AgustaWestland)**
- The MANOEUVRES Pilot Acoustic Indicator uses a pre-computed database of "static acoustic maps", reparameterised from (V, gamma, W) to (mu, CT, alpha_TPP) triplets. A rotor-state (blade-attitude) measurement system feeds the quasi-steady estimate. — [Trainelli et al., "Project MANOEUVRES", ERF 41, 2015](https://dspace-erf.nlr.nl/server/api/core/bitstreams/132a1b61-32a5-4197-9f86-3fe20f153d5b/content)
- Gennaretti et al. compared two quasi-steady indexings against a fully unsteady FW-H (Farassat 1A) maneuver solution on a hub-centred hemisphere. Approach A indexes by (mu, gamma), Approach B by (mu, CT, alpha_TPP).
  - B is better for OASPL when inertial loads change disk loading.
  - B is better for BVISPL when inertial loads change rotor attitude.
  - Neither is "fully adequate to replace unsteady aeroacoustic solvers" in arbitrary maneuvers.
  - — [Gennaretti, Serafini, Molica Colella, Bernardini, "Simulation of helicopter noise in maneuvering flight", ERF (c. 2014)](https://dspace-erf.nlr.nl/server/api/core/bitstreams/bfa55d24-072c-4cf5-a4b8-56f3b5770836/content)

**NLR hemisphere measurement**
- NLR measured AH-64 hemispheres with:
  - A 96 m vertical array of 93 microphones at 1 m spacing, with delay-and-sum beamforming in six-microphone clusters to remove the ground reflection.
  - A 13-microphone ground array.
  - A hemisphere that is a **regression model**, not a table: per-band polynomial (Taylor) or truncated Fourier series in (polar, lateral) angle, fitted by least squares.
  - Fitting **in SPL (dB) rather than p_rms^2**, because energy fits under strong directivity give large local errors and even negative energies.
  - — [Tuinstra & Sijtsma, "Measurement of helicopter noise hemispheres utilizing a 100 m vertical array", NLR/ERF (year not confirmed)](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)

**Empa sonAIR / FLULA2 (Switzerland)**
- sonAIR (Empa, developed 2012-2016, continuously updated) has an emission model in 1/3-octave bands from 25 Hz to 5 kHz, with 3-D directivity and engine/airframe separation. Emission depends on flight parameters: air density, fan speed N1, landing gear and high-lift devices. In the time-step simulation, "sound emissions are read and interpolated from a lookup table". — [Zellmann, "Aircraft noise calculations with sonAIR", CSA workshop, Oct 2019](https://kth.se/polopoly_fs/1.995593.1592494492!/Presentation_CSAworkshop_Zellmann.pdf)
- sonAIR source models range from detailed turbofan models, with airframe and engine directivity separated in 3-D, down to "simplified models for a wide variety of propeller driven aircraft and helicopters". **[snippet only; Empa site blocked]** — [Empa sonAIR page](https://www.empa.ch/web/s509/sonair)
- FLULA2 moves a point source with directional emission stepwise along the trajectory. Its database is measured (Zurich campaigns of 1991 and 1996). Its directivity is described as 2-D (longitudinal), against sonAIR's 3-D. **[snippet only]** — [Empa FLULA2 page](https://www.empa.ch/web/s509/flula2)
- DLR's PANAM, a fixed-wing system-noise module unrelated to the user's database builder, predicts emission of novel aircraft concepts. It has been combined with sonAIR for scenarios (Delfs et al. 2018, MDPI Energies). — [Zellmann 2019](https://kth.se/polopoly_fs/1.995593.1592494492!/Presentation_CSAworkshop_Zellmann.pdf)

**eVTOL / UAM tools and auralization**
- TU Delft / Politecnico di Torino build a Noise Hemisphere Database by sampling FW-H on a sphere or hemisphere around the vehicle. The flow comes from PowerFLOW LB/VLES or a low-fidelity solver, and the UYGUR Gaussian-beam tracer propagates it. — [Yunus, Casalino, Avallone, Ragni, "Efficient prediction of urban air mobility noise in a vertiport environment", Aerosp. Sci. Technol. 139 (2023) 108410](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)
- Korea (Chungbuk): a hemisphere was built from KARI's OPPAV flight test with 80 ground microphones at 60 m altitude (Goheung NATC), using a program called "SPHPERE-N" [sic]. Three variants of it, differing in how null cells are handled, gave different vertiport impact ranges (200-700 m horizontally). — [Lee & Haan, "Noise Impact Range Assessment Using UAM Noise Hemisphere", Forum Acusticum 2025, doi:10.61782/fa.2025.0610](https://dael.euracoustics.org/confs/fa2025/data/articles/000610.pdf)
- Hanseo University with KARI built a **time-domain** noise hemisphere database for auralization, by back-propagating microphone pressure time histories from the K-UAM Grand Challenge flight data to the hemisphere surface. They note that existing hemisphere methods work in the frequency domain, "complicating accurate auralization". — [Jo & Wie, "Prediction of environmental noise from rotorcraft using a noise hemisphere with time-domain data", Forum Acusticum 2025, doi:10.61782/fa.2025.0693](https://dael.euracoustics.org/confs/fa2025/data/articles/000693.pdf)
- NASA Auralization Framework (NAF): DLLs with a plugin architecture (the Advanced Plugin Libraries). It synthesises sound at the source from a source-noise definition and propagates it in the time domain. A receiver-based mode instead interpolates frequency-domain propagated spectra in time and handles the direct/ground-reflected path difference. — [Aumann, Rizzi, Heath, "Receiver-based auralization of broadband aircraft flyover noise using the NAF", NTRS 20200002559](https://ntrs.nasa.gov/api/citations/20200002559/downloads/20200002559.pdf)
- A NAF synthesis plugin computes rotor loading and thickness noise sample by sample (Farassat 1A) from blade loading, motion and geometry. It applies to steady periodic, unsteady periodic and aperiodic (transition) cases, so it needs no hemisphere interpolation. — [Krishnamurthy, Tuttle, Rizzi, "A Synthesis Plug-in for Steady and Unsteady Loading and Thickness Noise Auralization", NTRS 20200011447](https://ntrs.nasa.gov/api/citations/20200011447/downloads/20200011447.pdf)

**JAXA, Penn State**
- JAXA:
  - Measured source characteristics with the MuPAL-epsilon research helicopter (2005).
  - Measured propagation with microphones hung from two balloons at 50-200 m (October 2006).
  - Flew an onboard real-time route and flight-condition optimiser at Taiki-cho, Hokkaido.
  - — [JAXA "Sora to Sora" No. 20 (c. 2007, inferred from text)](https://www.aero.jaxa.jp/eng/about/publications/magazine/pdf/en_sorasora_20.pdf)
- Penn State's "comprehensive noise prediction system" (Botre, Brentner, Horn, Wachspress, VFS 75, 2019) was validated on level flight, descent, left and right level turns, a decelerating turn and a descending turn. It matched SEL contours "within a few SELdBA". I read the abstract only, so the source representation is undetermined. — [PSU Pure record](https://pure.psu.edu/en/publications/validation-of-helicopter-noise-prediction-system-with-flight-data/)

### Inferences
- The European research tools (SELENE, the HELENA/DLR work, MANOEUVRES) and FRAME have converged on the main-rotor aerodynamic state for indexing hemispheres: (mu, CT, alpha_TPP), or close equivalents such as gamma_e and n_z. RNM/AAM and the regulatory tools stay with kinematic parameters (airspeed, flight path angle). NICE-OPS's (mu, gamma or gamma_e, CT) sits between the two and is closest to SELENE and FRAME-QS.
- SELENE's design (Delaunay triangulation, nearest-neighbour fallback outside the hull) matches NICE-OPS's `interpolator.cpp` almost exactly. The others differ:
  - RNM interpolates linearly in energy, one axis at a time.
  - FRAME-QS uses natural-neighbour weights rather than barycentric ones.
- Whether a tool interpolates in dB or in energy is rarely stated. RNM says energy, and NLR argues that dB is the better domain for *fitting* hemispheres.
- The AEDT helicopter path cannot represent lateral (BVI) directivity beyond +/-45 deg, or more than about three flight modes. A hemisphere tool and AEDT will therefore disagree off-track by construction, as NASA's comparisons show.

### Gaps
- The full HNM technical manual and the AEDT/INM helicopter chapters could not be retrieved (403). So the details of HNM's static (hover/idle) directivity, such as the number of azimuths, and of AEDT's "dynamic operation modes" are unverified.
- sonAIR's helicopter emission model could not be confirmed: whether it is 3-D, which flight parameters it uses, and which helicopter types it covers. The Empa pages and documentation are WAF-blocked, and no helicopter-specific sonAIR paper was found.
- HELENA's exact database parameters and interpolation scheme were not found in a primary source. That DLR's FRIENDCOPTER work used (mu, CT, alpha_TPP) is confirmed; that HELENA uses them is inferred.
- The paper read does not state CAROT's flight-condition parameterisation (airspeed and vertical speed are implied by the H130 figure) or its interpolation. The primary CAROT paper (Dieumegard et al., VFS 78, 2022) was not retrieved.
- No ONERA- or Leonardo-owned ground-footprint tool was found by name. MANOEUVRES (with AgustaWestland as partner) is the only Leonardo-linked item. HELISHAPE/HELINOISE and "CAMPOS" were not checked before the search quota ran out.
- The ANOPP2 tool names "AFFAL", "ASNIM" and "ARNM" in the brief could not be verified and may not exist.
- No KARI-owned footprint tool name was found, only university/KARI joint hemisphere papers.

## 2. Propagation physics: absorption standard, ground model, terrain, refraction/wind, turbulence, shielding/urban

### Takeaway
Most hemisphere tools propagate along straight rays, with spherical spreading, standard air absorption and a two-ray impedance-ground model. RNM, for example, uses ANSI S1.26 absorption, a Chien-Soroka/Chessell ground and, from v3, GTD terrain. Refraction is a bolt-on where it exists:

- **SELENE:** ray tracing.
- **RNM 7:** a stratified atmosphere with winds.
- **TU Delft UYGUR:** Gaussian beams.
- **NASA APET:** SAM/RTP/GFPE attenuation tables.

Shadow zones are handled crudely in the ray tools; SELENE floors the level. Urban shielding appears only in sonX (reflection and scattering), in Gaussian beams over terrain (UYGUR), and in the Korean ISO 9613-2-based building studies. Footprint tools essentially omit turbulence; it appears only in auralization (Airbus GENEPASS).

### Cited Findings
- RNM (2002) atmosphere:
  - Sound speed comes from the temperature profile by the ideal gas law; the default is the U.S. Standard Atmosphere 1976.
  - Air absorption follows **ANSI S1.26-1978**, using the humidity profile.
  - Weighted-average sound speed and absorption are tabulated at 1000 ft intervals and read without interpolation.
  - Atmospheric loss is the slant distance times the altitude-averaged coefficient.
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM (2002) level equation: L(r) = L(r0) + A_spread + A_atm + A_grd + A_wind + A_topo.
  - The ground term is a complex impedance after Chien & Soroka and Chessell, with Daigle's corrections, evaluated at Doppler-shifted frequencies. The validation used a flow resistivity of 300 cgs Rayls.
  - Ray paths are straight, with no wind; wind was "planned".
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM v3.0 added terrain, using Rasmussen's GTD extensions to finite impedance.
  - A 2-D slice is classified as flat, concave, convex, a wedge/screen with one flat, or a wedge/screen with two flats.
  - Shielding "transition[s] continuously from a Pierce wedge to a Maekawa screen", applied per 1/3-octave band.
  - Echoes from canyons or surfaces behind the receiver are not treated.
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM 7 adds propagation "through a horizontally stratified atmosphere over uniform terrain with winds", alongside GTD over varying terrain. — [NASA Tech Briefs RNM (LAR-17753-1), 2016](https://www.techbriefs.com/component/content/article/24556-lar-17753-1)
- APET (Acoustic Propagation and Emulation Toolset): SAM builds a 2-D attenuation table over terrain with a Ray Tracing Program and two GFPE codes ("ray bending, spreading, atmospheric absorption, ground effects, and variable weather"). A modified RNM 7 consumes the table, and AUDIB computes detectability. — [NASA Tech Briefs APET (LAR-17754-1), 2017](https://www.techbriefs.com/component/content/article/27046-lar-17754-1)
- SELENE propagation:
  - It includes spherical spreading, Doppler shift, wind, atmospheric absorption and ground reflection.
  - Its ray-tracing module (the high-frequency limit of the linearised Euler equations) handles wind and temperature gradients.
  - Where the eigenray search fails, it assumes a shadow zone and sets SPL to a low value rather than minus infinity.
  - Caustics near the shadow boundary produce irregular contours.
  - Under the track, contours change mainly because wind changes the flight condition; far from it, mainly because of propagation.
  - — [Olsman & Lummer, ERF](https://dspace-erf.nlr.nl/bitstreams/61aeceb7-6f3b-486f-8d7e-d335cd4db426/download)
- DLR's own description: SELENE propagates "either propagation in a still atmosphere or by ray tracing, to incorporate atmospheric effects such as wind and/or temperature gradients". — [DLR SELENE page](https://www.dlr.de/en/as/research-and-transfer/software-solutions/acoustics/software-selene)
- HELENA's propagation models are "specifically tailored for rotorcraft noise" and include propagation distance, atmospheric absorption and ground reflection. — [Smith et al., ERF 2014](https://dspace-erf.nlr.nl/bitstreams/c920dfb8-f809-440a-9387-af4cf81fcdd9/download)
- CAROT runs over idealized flat homogeneous terrain or over "realistic scenarii" with local topography and demography. In the auralization chain, GENEPASS "accounts for the ground effect and the atmospheric turbulence in the synthesis of the signal". — [Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)
- sonX works in 1/3-octave bands from 25 Hz to 5 kHz. It computes direct sound in a homogeneous atmosphere, with three optional refinements:
  - METEO: a meteorological correction.
  - REFLECT: specular reflections and scattering from buildings and walls.
  - FOREST: diffuse reflections from forest edges and cliffs.
  - Attenuations are read from a database, or computed directly for simple situations.
  - — [Zellmann 2019](https://kth.se/polopoly_fs/1.995593.1592494492!/Presentation_CSAworkshop_Zellmann.pdf)
- sonX's ground effect is an analytical spherical-wave solution extended to finite segments and variable ground. Meteorology is handled by ray tracing, and coherent reflections from rigid surfaces by a Kirchhoff-Helmholtz integral. **[snippet only]** — [Empa sonX page](https://empa.ch/web/s509/sonx)
- For the Geneva 2017 annual scenario, sonAIR's physics-based ground model gives less smooth contours than FLULA2 but equivalent results on average. **[snippet only]** — [Acta Acustica 2025, aacus250029](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2025/01/aacus250029/aacus250029.html)
- NASA AMAT propagates "using straight ray propagation". In the AEDT-versus-AAM UAM comparison, AAM's ground model "considers direct and reflected rays and can exhibit significant 'dips'", unlike AEDT's. — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615); [Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)
- UYGUR:
  - It traces rays from the aircraft centre and solves the dynamic ray-tracing equations for spreading and travel time.
  - It sums the Gaussian beams that pass near each receiver.
  - It handles 3-D terrain, multiple reflections and wind profiles.
  - Against straight rays it reduces prediction error "by 5 dB in the illuminated zones and 35 dB in the terrain shadow zones".
  - — [Yunus et al. 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)
- The Korean vertiport studies:
  - Lee & Haan propagate per **ISO 9613-2**, with a 5 m source radius, and draw horizontal and vertical (façade) maps around buildings in Seoul. — [Lee & Haan 2025](https://dael.euracoustics.org/confs/fa2025/data/articles/000610.pdf)
  - Jo & Wie cite ISO 9613-2:2024 and Daigle et al. for ground effects, applied to STFT narrowband data. — [Jo & Wie 2025](https://dael.euracoustics.org/confs/fa2025/data/articles/000693.pdf)
- JAXA's balloon measurements found that vertical changes in wind speed cause wide variation in propagation, while the no-wind case "agree[s] well with the predicted values". The optimiser plans curved routes "in consideration of noise propagation changes with wind". — [JAXA Sora to Sora No. 20](https://www.aero.jaxa.jp/eng/about/publications/magazine/pdf/en_sorasora_20.pdf)
- Contrast: NORAH2's propagation was benchmarked against NORD2000. — [EASA NOISE final report](https://www.easa.europa.eu/en/downloads/139729/en)

### Inferences
- NICE-OPS's straight / linear / stratified / PE hierarchy goes further than any footprint tool found here. Only two pair a full-wave solution with a hemisphere footprint tool: APET (GFPE tables feeding RNM) and NICE-OPS's `--propagation_table`. APET's tables are aimed at terrain and detectability.
- SELENE's published shadow handling (a floor wherever the eigenray search fails) is the same as NICE-OPS's behaviour before ShadowTransition. No tool reviewed publishes a creeping-wave, residue-series or boundary-transition treatment.
- The absorption standards differ: ANSI S1.26-1978 in RNM and the NORAH benchmark, ISO 9613-1/-2 in the Korean studies and NICE-OPS. The choice can move long-range high-frequency levels by around 1 dB/km, which is worth flagging in any cross-tool comparison.

### Gaps
- The sources read do not state the absorption standard or ground-impedance model used by SELENE, HELENA or CAROT.
- Whether sonAIR's helicopter models use the full sonX chain (METEO, REFLECT) is unconfirmed.
- AMAT's ground and absorption options are not documented in what was read.
- No source was found on turbulence (coherence loss, scattering into shadows) in any hemisphere footprint tool, outside auralization.

## 3. Doppler / emission-time handling, maneuvers, time step and metrics

### Takeaway
The tools march along the trajectory in source time and place each contribution at its arrival time. RNM bins contributions by arrival time; NAF and Hanseo synthesise in the time domain. Frequency Doppler is applied in:

- **RNM:** in its ground model.
- **SELENE.**
- **CAROT:** by re-dopplerizing de-dopplerized tonal hemispheres.

Maneuvers are handled quasi-statically: a steady hemisphere is chosen by (mu, CT, alpha_TPP) or gamma_e and oriented by attitude. Studies agree this beats matching path slope but does not equal unsteady FW-H. The common footprint metric is SEL integrated within 10 dB of the peak (HELENA, SELENE, the AAM comparisons), as in NICE-OPS.

### Cited Findings
- RNM "is a source time based model, which uses the flight track based source time directly, 'binning' sound at the receiver locations based on arrival time".
  - Only after propagation are the receiver sounds interpolated onto a uniform time mesh and summed.
  - Broadband arrival time uses the vehicle centre (incoherent).
  - Pure tones keep absolute phase from the first trajectory point, so multiple hemispheres interfere coherently.
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM's multiple-track mode defaults to a maximum track-point spacing of 2.0 s, paired with a community grid of 500 ft or coarser at up to 120 kt. Segment interpolation assumes constant acceleration. — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM metrics:
  - Single tracks: SPL, SPL(A), PNLT, SEL, SEL(C), SEL(A), Lmax(A), EPNL.
  - Multiple tracks: DNL, CNEL, NEF, WECPNL, Leq.
  - Output goes to GIS, Tecplot and NMPLOT.
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- SELENE:
  - The flight path is built from waypoints with Bézier splines, with flyability checks (height-velocity diagram, vortex ring state, fuselage pitch).
  - The emission time tau solves t - tau - |R|/c0 = 0, and a Doppler amplification factor is included.
  - SEL runs from the first to the last point within 10 dB of the peak.
  - Wind also changes the flight condition, through the flight-mechanics code.
  - — [Olsman & Lummer](https://dspace-erf.nlr.nl/bitstreams/61aeceb7-6f3b-486f-8d7e-d335cd4db426/download)
- HELENA: "t2 and t1 are defined as the instances at which the A-weighted Sound Pressure Level decreases 10dB below the maximum level". — [Smith et al., ERF 2014](https://dspace-erf.nlr.nl/bitstreams/c920dfb8-f809-440a-9387-af4cf81fcdd9/download)
- CAROT/MOTUS:
  - Tonal hemispheres are stored de-dopplerized and re-dopplerized "according to the computed emission-to-reception vectors over time".
  - Contributors with different frequency vectors are summed as uncorrelated.
  - The output feeds GENEPASS for auralization.
  - — [Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)
- NLR's hemisphere acquisition solves for emission time iteratively, t_m = t_e + |X_S(t_e) - X_mic|/c, using 0.2 s FFT blocks with 50% overlap. — [Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)
- Quasi-steady maneuvers:
  - Gennaretti et al.: indexing by (mu, CT, alpha_TPP) beats (mu, gamma), but neither replaces unsteady FW-H. — [Gennaretti et al.](https://dspace-erf.nlr.nl/server/api/core/bitstreams/bfa55d24-072c-4cf5-a4b8-56f3b5770836/content)
  - FRAME-QS: spheres are selected and oriented by gamma_e and n_z through transient pitch-up and roll maneuvers. — [Greenwood et al. 2015](https://ntrs.nasa.gov/api/citations/20160005979/downloads/20160005979.pdf)
- DLR flight tests:
  - Maneuvers have strong effects on dBA snapshots.
  - Accelerations change CT, so pilots must fly a "tunnel" smoothly.
  - The pilot guidance: avoid long periods under 50 kt, and fly faster than 90 kt where possible.
  - — [Spiegel et al., ERF 2008](https://dspace-erf.nlr.nl/bitstreams/16abfa07-16f9-4711-adbf-2a87e3cbd01b/download)
- Crosswind:
  - BVI changes markedly when the fuselage is held along the inertial path, but not when the aircraft crabs to zero aerodynamic sideslip.
  - Wind gradients matter less crosswind than headwind.
  - — [Greenwood & Sim, NTRS 20130013862](https://ntrs.nasa.gov/api/citations/20130013862/downloads/20130013862.pdf)
- AEDT helicopter mode treats approach and departure as constant-state segments, with noise-fraction and duration adjustments. Fixed-wing mode adds an extra noise-fraction term beyond the last approach segment, which helicopter mode lacks. For the same UAM approach this gives "bulb"-shaped contours in one mode and "wedge"-shaped ones in the other. — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- NAF synthesises sample by sample at the instantaneous emission angle for transitional (unsteady) eVTOL flight, so it avoids hemisphere interpolation. — [Krishnamurthy, Tuttle, Rizzi, NTRS 20200011447](https://ntrs.nasa.gov/api/citations/20200011447/downloads/20200011447.pdf)

### Inferences
- The 10-dB-down SEL window is not peculiar to NICE-OPS or the D original. HELENA and SELENE document it explicitly, and the NASA comparisons use it in AAM, so NICE-OPS's first-crossing fragility (CLAUDE.md) probably exists in those tools too. SELENE defines t1 as the *last* point within 10 dB, which reads like the outermost-crossing rule that NICE-OPS measured as an alternative.
- None of the tools reviewed documents a wind-relative ("air") flight condition as a database convention. SELENE handles it implicitly through its flight-mechanics code, and Greenwood & Sim show that the trim/crab choice matters. This supports NICE-OPS's explicit `speed_reference` attribute.

### Gaps
- The sources read do not state the time-step defaults of SELENE, HELENA, CAROT or AMAT.
- I found nothing on whether RNM 7 or AMAT apply amplitude Doppler (convective amplification) or only the frequency shift.

## 4. Validation and known limitations

### Takeaway
Published validation is mostly qualitative or point-wise, within about 2-3 dB SEL:

- **RNM:** the NATO Bell 412 trials.
- **JAXA:** within 2 dB at five points.
- **PSU:** "within a few SELdBA".
- **SELENE and HELENA:** procedures verified in flight.

The recurring limitations:

- Extrapolated hemisphere regions near the tip-path plane.
- The quasi-steady assumption in maneuvers.
- Straight-ray propagation with no wind or shadow physics.
- For the regulatory tools, the coarse three-angle lateral directivity.

### Cited Findings
- RNM was compared with the NATO/CCMS Bell 412SP (CH-146) trials at CFB Moose Jaw, June 1998: a 6 deg / 103 kt approach on a 250 ft grid, with Euler angles set to zero because attitude was not measured.
  - Agreement was excellent under the track and on the retreating side.
  - RNM over-predicted the advancing side by up to 3 dB SEL.
  - The time history was over-predicted early and late, near the TPP.
  - Both errors were attributed to **extrapolated hemisphere data**, beyond the lateral reach of the ART2 measurement array.
  - — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- RNM was adopted by DoD and NATO as the standard tool for military rotorcraft environmental impact assessment, and is distributed with NOISEMAP/BASEOPS. — [Conner & Page 2002](https://ntrs.nasa.gov/citations/20030003730)
- DLR flight validation:
  - The SELENE-optimised approach procedure (Clean Sky) "was successfully verified by flight test with noise measurements on the ground". — [DLR SELENE page](https://www.dlr.de/en/as/research-and-transfer/software-solutions/acoustics/software-selene)
  - DLR/Eurocopter noise-abatement procedures gave measured reductions "close to 10" dB. — [Spiegel et al. 2008](https://dspace-erf.nlr.nl/bitstreams/16abfa07-16f9-4711-adbf-2a87e3cbd01b/download)
- The DLR "Hemisphere" method was validated with BO105 data from the 2001 RONAP flight test at Cochstedt, including 6 deg descent and 12 deg takeoff footprints. **[snippet only]** — search summary of [ERF/TU Delft records](https://repository.tudelft.nl/file/File_24b06680-4b05-4fbb-bff1-ff8cabf02979)
- JAXA: the optimised route was about 4 dB quieter than a simulated ILS approach at five ground points, with "error between the predicted and measured noise ... within 2 dB". — [JAXA Sora to Sora No. 20](https://www.aero.jaxa.jp/eng/about/publications/magazine/pdf/en_sorasora_20.pdf)
- PSU: SEL contours "within a few SELdBA" across level, descent and turning maneuvers. — [Botre et al. 2019](https://pure.psu.edu/en/publications/validation-of-helicopter-noise-prediction-system-with-flight-data/)
- Airbus MOTUS limitations:
  - Without source separation, a total-noise hemisphere database "cannot identify the variations of the several noise contributors".
  - Contributors have so far been extracted only for "flight conditions of interest".
  - Absolute levels from simulation remain hard to get right, which is why the masks use simulated *differences*.
  - — [Guntzer et al., ERF 2023](https://dspace-erf.nlr.nl/items/4699deee-c7f6-4304-b8a7-37a2d80d44d1)
- AEDT versus AAM for a UAM quadrotor:
  - With the spheres calibrated to match LAmax at 1000 ft under the track, SELA still differed by 1.7 dB, because polar directivity accumulates in SEL.
  - Off-track differences come from AEDT's 0 deg-azimuth-only NPDs and its simpler ground model.
  - — [Rizzi, Page & Cheng 2021](https://ntrs.nasa.gov/api/citations/20210015134/downloads/IN21_AAM_v12.pdf)
- AEDT's helicopter lateral directivity adjustment "creates large differences to the sides of the flight path" compared with fixed-wing mode. With axisymmetric helicopter NPDs, the overflight difference becomes negligible. — [Letica & Rizzi, NTRS 20240003615](https://ntrs.nasa.gov/citations/20240003615)
- Hemisphere processing:
  - Energy-domain fitting can yield negative energies under strong directivity.
  - Treating an 18 m AH-64 as a point source at 90 m standoff gives about +/-5.5 deg polar-angle uncertainty.
  - — [Tuinstra & Sijtsma](https://dspace-erf.nlr.nl/bitstreams/10f99d37-2225-43e9-9bf1-ac584ac6cdd7/download)
- Korean UAM: how null hemisphere cells were filled changed the predicted impact range. Coverage and extrapolation of measured hemispheres is therefore a first-order uncertainty. — [Lee & Haan 2025](https://dael.euracoustics.org/confs/fa2025/data/articles/000610.pdf)
- Against Gaussian beams, straight-ray propagation in vertiport and terrain settings errs by about 5 dB in lit zones and about 35 dB in terrain shadow. — [Yunus et al. 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)

### Inferences
- The best-documented failure mode across tools is on the source side: extrapolated or missing hemisphere regions near the TPP and laterally beyond the array (RNM/NATO, the Korean null cells). This matches NICE-OPS's own finding that forward-line errors are dominated by the approach's end and by low, slow states rather than by propagation.
- No tool reviewed publishes a quantitative comparison against a full-wave reference (PE or exact solutions) of the kind NICE-OPS holds itself to. In what was found, APET's GFPE builds tables; it is not used as a validation oracle.

### Gaps
- No quantitative validation statistics (bias or rms over many microphones) were found for SELENE, HELENA, CAROT or sonAIR helicopters.
- The HELENA validation paper (Gervais et al., AHS 66, 2010) and the CAROT paper (Dieumegard et al., VFS 78, 2022) were identified from reference lists but not read.
- Empa's helicopter validation, if any, could not be accessed.
