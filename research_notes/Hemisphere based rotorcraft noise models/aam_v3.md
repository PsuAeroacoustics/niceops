# Advanced Acoustic Model (AAM) v3.1.0: how it represents sources and computes ground noise

Primary source: *Advanced Acoustic Model (AAM) Technical Reference and User's Guide*, May 2021 (report date 25 May 2021, PDF created 27 May 2021), J.A. Page, A. Rapoza, A. Oberg, A. Hastings, G. Baker, M. Shumway, Volpe National Transportation Systems Center. It is 250 PDF pages, and the printed page numbers are the PDF page minus 16. Secondary source: `ReadMe-AAM-v3.1.0.txt` from the same distribution. Every citation below gives the printed page number ("TR p.N") with a section number where one exists. Links open the local PDF at the right page. Note on provenance: development moved from Wyle/NASA to Volpe in spring 2020 ([TR p.ii Announcements][pAnn]). "Wyle/Blue Ridge" is therefore the code's lineage. The v3.1.0 document is Volpe's.

## 1. Source sphere format: geometry, angle convention, spectral content, stored quantities, metadata, QuietSphere

### Takeaway
An AAM source is a set of NetCDF files, one per noise source and per flight condition. Each file holds levels in dB on a regular (theta, phi) grid at a fixed radius, `AMPLITUDE(PHI, THETA, FREQUENCY)`, plus a phase array for pure tones. Theta is the polar angle from the nose (0 deg) to the tail (180 deg). Phi is the roll angle about the body x-axis: 0 deg straight down, +90 deg starboard, -90 deg port. This is the ANOPP convention. Phi from -90 to +90 makes a lower hemisphere; phi from -180 to 180 makes a full sphere. Spectra are 1/3-octave bands 10 Hz–10 kHz, 121 1/12-octave bands, arbitrary narrowbands, or pure tones with phase. No PSD form is documented. QuietSphere is only a utility that shifts every 1/3-OB level of an existing sphere by a constant number of dB.

### Cited Findings
- Source characteristics are "defined by sets of sound spheres (or hemispheres), each sphere being centered on a noise source of the aircraft, or as a single compact source at the center". They may be 1/3-octave (broadband), narrowband, 1/12-octave, or "pure tone sound pressure levels and phase" ([TR p.1 §1.1][p1]).
- "A sound sphere is stored as a single computer file containing the acoustic data for a single flight condition". All data reach AAM in NetCDF ([TR p.8–9 §2][p8]). "Sound sphere (or hemisphere) files each contain one noise source for a single flight condition" ([TR p.30 §2.5][p30]; repeated at [TR p.163 §4.1.4][p163]).
- Points on the sphere are "described in terms of a fixed radius and two spherical angles". The origin is at the centre of the noise source and the X-axis is aligned with the body-fixed system. Theta is zero toward the nose and 180 deg toward the tail. "Phi is zero directly underneath the vehicle, positive 90° on the starboard (right) side, and minus 90° on the port (left) side." The notation is "consistent with that used by ... ANOPP" ([TR p.29 §2.5][p29]; same wording in ART [TR p.147 §4.1.1][p147]). The sphere origin is customarily at the CG for fixed-wing and tiltrotor aircraft and "at the rotor hub for helicopters" ([TR p.24 §2.4.1][p24]).
- NetCDF variables (Table 4-37):
  - `BB`, `NB`, `PT`: type flags. A sphere holds exactly one type, so BB + NB + PT = 1.
  - `DOPPLER_SHIFT_REMOVED`: 0/1.
  - `EMPTY_WEIGHT`, `FUEL_WEIGHT`, `LOAD_WEIGHT`.
  - `RADIUS` (ft or m).
  - `FLIGHT_PATH_ANGLE`: positive = departure, negative = arrival.
  - `PYLON_ANGLE`: 0 deg = airplane mode, 90 deg = helicopter mode.
  - `SPEED` (knots, m/s or mph).
  - `AOA`.
  - `XYZ`: sphere centre in body-fixed coordinates.
  - `PHI`: listed as "-90° to +90°".
  - `THETA`: 0–180 deg.
  - `FREQUENCY` (Hz): band centre or tone frequency.
  - `AMPLITUDE(Phi, Theta, Frequency)` (dB).
  - `PHASE(Phi, Theta, Frequency)`: pure tone only, "relative to reference time".
  - `POWER`, `THRUSTANGLE`, `ACCONFIG`, `INTERPOLATIONSTATE`: fixed-wing.

  ([TR p.163–164 §4.1.4][p163]).
- Storage order: frequency varies fastest, then theta, with phi slowest (Fortran order) ([TR p.168 §4.7][p168]).
- Full spheres do exist even though Table 4-37 says phi runs -90 to +90. The auto-generated `OMNI_200.NC` has PHI = -180, -175, …, 180 (73 values), THETA = 0…180 in 5 deg steps (37 values), 31 1/3-OB frequencies 10 Hz–10 kHz, RADIUS = 100 ft, all levels 200 dB, DOPPLER_SHIFT_REMOVED = 0 ([TR p.204–205 §5.8.1][p204]). It is called a "200 dB omni directional full sphere with 5 degree grid resolution" ([TR p.204][p204]). ART's SPHERE keyword allows "phi values range from -180 to 180 deg" ([TR p.159 §4.1.2.31][p159]).
- Hemispheres: the distributed eVTOL files use 5 deg spacing with phi -90…+90 and theta 0…180. "These two source files represent hemispheres, with only the lower portion populated" ([TR p.214 §5.12.1][p214]). ART's WRAPAROUND "takes a hemisphere (lower half) ... and mirror images the values on the top half of the sphere" ([TR p.161 §4.1.2.38][p161]).
- Typical resolution: "Typically, the data are provided in 5-degree increments on the sphere" ([TR p.144 §4.1][p144]). The SPHERE.EXE example is a 37 × 37 grid ([TR p.175 §4.17][p175]).
- Radius:
  - ART's default is 100 ft ([TR p.60–61, ART keyword][p60]).
  - The distributed eVTOL spheres use 1000 ft ([TR p.214][p214]).
  - Figure 2-17 shows a fixed-wing sphere at a "Reference Distance of 1000 Ft." ([TR p.29][p29]).
  - Bullwink rescales input data from the .DAT radius to the requested output radius "according to pure geometric spreading" ([TR p.166 Table 4-41 note 1][p166]).
  - Spheres combined in an interpolation must share a radius, or AAM stops with "Radius mismatch" ([TR p.137 §3.6.3][p137]).
- Spectral formats:
  - 1/3-OB band numbers must lie in 10–40, i.e. 10 Hz–10 kHz ([TR p.89 USERBAND][p89]; [TR p.139 errors][p139]).
  - 1/12-OB uses exactly 121 bands from 10 Hz to 10 kHz ([TR p.162 Table 4-36][p162]; band list [TR p.170–171 §4.11.3][p170]).
  - Narrowband allows user-defined bands and spacing, "with arbitrary bands permitted" ([TR p.2][p2]). Up to 1024 frequencies are allowed in Bullwink ([TR p.166][p166]). Sparse NB files carrying only selected tones (for example main- and tail-rotor harmonics) are allowed ([TR p.61 AUDIB][p61]).
  - Pure tones carry level and phase ([TR p.28–29 §2.5][p28]).
- Missing data: regions with no valid data are set to 1E35, and "AAM will not propagate the source point to the receiver location" ([TR p.145 §4.1][p145]).
- Per vehicle: up to 1000 spheres, files `vvvvvnnn.nc` (5-character vehicle ID plus a 3-digit index) ([TR p.76 SETUP PARA][p76]; [TR p.145][p145]). There can be up to 10 broadband, 10 narrowband and 10 pure-tone spheres per flight condition at different body-fixed centres ([TR p.30][p30]; [TR p.163][p163]).
- QuietSphere "reads in ASCII data, allows the user to add or subtract a fixed dB value in all directions and across frequency bands". It writes an ASCII NetCDF-style file for NCGen, under a new 5-character vehicle name ([TR p.172 §4.13][p172]). The distribution table says it "Takes a NCDump format ASCII file (1/3 OB only) and adds x dB to all bands" ([TR p.46 Table 3-3][p46]). The v3.1.0 ReadMe lists an "Updated QuietSphere utility for more flexible input formats" ([ReadMe][readme]).

### Inferences
- The angle convention (theta from the nose, phi a roll angle about the body x-axis with 0 = down) is the same as NICE-OPS's: phi = 0 straight down, source spheres on [-90, 90], upper surface at 180 - phi. The `DOPPLER_SHIFT_REMOVED` attribute name used by panam/NICE-OPS is AAM's. This suggests panam's .nod layout descends from the AAM/ART NetCDF convention. That is an inference from the matching names, not something the document says.
- AAM stores one condition per file. NICE-OPS stores many conditions as groups in one `.nod` file. AAM has no in-file representation of a condition grid: the "database" is the directory listing, sorted by header values at run time (see §2).
- QuietSphere is a what-if utility (for example, a notional vehicle N dB quieter). It is not part of the propagation chain and plays no role in sphere generation or selection.

### Gaps
- The document never says whether `AMPLITUDE` is a band SPL at the sphere radius in a lossless field. "Free field loss-less sound level at a distance r0" in Eq. 1 ([TR p.10][p10]) implies it is, but the file format is never stated to be atmospheric-absorption-free. ART's depropagation removes absorption, which also implies it (§11).
- No PSD storage is documented. Narrowband levels are per band, and the bandwidth is passed separately through the AUDIB keyword ([TR p.61–62][p61]).
- How a full sphere's phi range from -180 to 180 is indexed at lookup is not stated.

## 2. Parameterization by operating state; inter-sphere selection and interpolation; maneuvers

### Takeaway
Rotorcraft spheres are indexed by **nacelle (pylon) tilt angle, flight path angle and airspeed**. Conventional helicopters use FPA and speed only. Lookup is hierarchical:
- vehicle type, sphere type and sphere centre must match exactly;
- the nacelle angle is snapped to the closest available;
- up to four spheres are chosen and blended **linearly in dB**, first in speed at two bracketing FPAs, then in FPA. This is a two-pass, bilinear-style blend.

There are hard limits on how far apart spheres may be: dFPA 3 deg on approach, 5 deg on departure, and dSPD 20 kt. When the database is sparse, AAM falls back to the nearest condition. Thrust or CT, RPM, weight and load factor play no part, with two exceptions: Q-SAM1/Q-SAM2 substitute an effective flight path angle, and FRAME databases add load factor (through an "effective weight") as a third dimension with up to 8 spheres. Fixed-wing spheres use power setting and thrust-vector angle, gated by configuration state and an INTERPOLATIONSTATE flag. NCSPEC bypasses all of this and names the file per trajectory point.

### Cited Findings
- Header data read per file:
  - sphere type;
  - flight path and nacelle tilt angles (rotorcraft);
  - thrust vector angle and thrust setting (fixed-wing);
  - airspeed (rotorcraft);
  - radius;
  - Doppler state.

  "This information is tabulated and sorted" ([TR p.9 §2.1][p9]).
- Rotorcraft: "The three independent variables describing the quasi-steady flight condition for a rotorcraft are the flight path angle, the nacelle tilt angle, and the airspeed. For conventional helicopters, only flight path angle and airspeed are needed" ([TR p.30 §2.5][p30]; [TR p.163 §4.1.4][p163]).
- History: RNM v1 "considered only airspeed when linearly interpolating between sound spheres" [Lucas 1998]. The present version does "a two-dimensional interpolation considering both airspeed and flight path angle", or takes explicit per-segment files with NCSPEC ([TR p.31 §2.5.1.1][p31]).
- "At each point along the flight trajectory the source sound characteristics are extracted from the sphere files using a linear interpolation at the required airspeed, followed by a second linear interpolation on the vehicle flight path angle." For sparsely populated databases "the prior RNM methodology of selecting the closest matching flight condition available sphere without interpolation is utilized" ([TR p.31][p31]).
- Since RNM v7, AAM "rounds the .NC file's elevation angles to the nearest 0.5-degree increment before sorting". The reason given is that ART assigns a nominal condition by averaging trajectory parameters ([TR p.31][p31]).
- Sort order (Table 2-5): vehicle type; XYZ centre; sound type; nacelle tilt (0 to 90 deg); flight path angle (approach to departure); speed in knots ([TR p.31][p31]).
- Search criteria (Table 2-6):
  - vehicle type, file type and sphere location must match, with no interpolation;
  - nacelle angle: "No interpolation, closest chosen";
  - FPA: "Linear interpolation, dFPA deg max";
  - speed: "Linear interpolation, dSPD kts max".

  If any of the first three does not match, execution halts ([TR p.33][p33]).
- "Only a single nacelle tilt angle is permitted for all four sources". On an exact tie the larger nacelle angle is taken ([TR p.33][p33]).
- Interpolation limits (Table 2-7): approach FPA 3.0 deg, departure FPA 5.0 deg, speed 20.0 kt ([TR p.33][p33]).
- "Up to four sound spheres ... are linearly interpolated in decibels based on speed and flight path angle" (Figure 2-20: I1, I2 at FPAng1 and I3, I4 at FPAng2, each pair interpolated to SPD_desired, then to FPA_desired) ([TR p.33–34][p33]).
- Selection rules (Table 2-8):
  - FPA selection precedes speed selection and "MAY disqualify .NC files".
  - Speed selection "may NOT disqualify .NC files (i.e., in the absence of flight data within the specified dSPD criteria, ANY available sphere speed may be chosen)".
  - FPAng1 ≤ FPA_desired ≤ FPAng2.
  - On a tie at the speed midpoint "the largest speed source file will be selected (generally a conservative result)".
  - Table 2-8 also says "(A duration adjustment will be made by AAM in the course of its normal calculation procedure)" ([TR p.34][p34]).
- Errors: spheres that share the same nacelle, FPA and speed trigger "Same profile data found in multiple spheres" ([TR p.137–138][p137]). Interpolated spheres must agree in band set, NB frequencies and PT frequencies, with PT frequencies "within 3% of one another" ([TR p.134–139][p134]).
- Fixed-wing selection:
  - Spheres are indexed by power setting in its own units (%N1, EPR, ...), with a thrust-vector-angle limit of dTVA 20 deg.
  - Interpolation happens only within the same configuration (clean, flaps, flaps + gear, run-up, afterburner variants, F-35B Mode 4, hover; integer codes in Table 2-12) and only when INTERPOLATIONSTATE = 1.
  - There is no extrapolation beyond the data.
  - Interpolation order: power first, then TVA.

  ([TR p.34–38 §2.5.1.2, §2.5.2][p34]).
- Fixed-wing altitude thrust correction: 10^(0.00002(1000−Z)) above 1000 ft MSL, about −2 dB per 10,000 ft ([TR p.40 §2.5.4][p40]).
- NCSPEC: "explicitly define which noise sphere to use at each trajectory point" (single-event mode only). It is recommended for Advanced Air Mobility vehicles because automatic lookup "is limited to the NC file header data" ([TR p.30–31][p30]; [TR p.69][p69]). With NCSPEC, all files must match in phi, theta and frequency dimensions ([TR p.138][p138]).
- Q-SAM1 (Gopalan 2004; Gopalan & Schmitz 2003): γ_eff = γ + v̇/g (Eq. 13). The text says "a 0.1 g acceleration ... corresponds to a 0.1 radians (5.7 degrees) climb angle". The physical link is that the tip-path plane during slow acceleration matches steady flight at γ_eff, more nose-up in deceleration. QSAM1 makes AAM "select and rotate the helicopter source characteristics using the effective flight path angle" ([TR p.39 §2.5.3.1][p39]). It is rotorcraft-only ([TR p.73][p73]).
- Q-SAM2 [Greenwood 2015]:
  - Each axis of the waypoints is fit with a cubic spline that misses no waypoint by more than 1 ft per axis.
  - a⃗ = a_x î + a_y ĵ + (a_z + g) k̂ (Eq. 14); thrust direction n̂ = a⃗/|a⃗| (Eq. 15).
  - γ_eff = sin⁻¹(n̂·v⃗/|v⃗|) (Eq. 16).
  - Roll from the projection of velocity into the inertial tip-path plane: û = (v⃗ − (v⃗·n̂)n̂)/|v⃗| (Eq. 17), tan φ = −((k̂×û)·n̂)/(n̂·(û×k̂×û)) (Eq. 18).

  ([TR p.39–40 §2.5.3.2][p39]; equations checked against the page image). Note: the reference list has no "Greenwood 2015" entry.
- FRAME (Greenwood's "Fundamental Rotorcraft Acoustics Modeling from Experiments"):
  - The sphere database is "evenly spaced in flight path angle and speed, and includes load factor as a new dimension".
  - A cubic spline gives accelerations, from which AAM "determines roll and pitch angles and selects up to 8 spheres for interpolation".
  - The lookup variable is "effective weight" = base weight (QSAM2 keyword input, lbs) × load factor.
  - Search is FPA first, then speed, then load factor. This only works when spheres are matched on the FPA and speed grid.

  ([TR p.41–42 §2.5.7][p41]; [TR p.74 QSAM2][p74]).
- Airspeed versus ground speed: "At present, AAM does not account for the wind when determining airspeed and performing a source data look-up. Winds are taken into account by the RNMwea propagation module" ([TR p.226 Glossary "Flight Profile"][p226]; also "Ground speed", [TR p.227][p227]). The multiple-operations rotorcraft profile specifies "Indicated airspeed (kts)" ([TR p.115 Table 3-64][p115]).

### Inferences
- AAM's rotorcraft condition space is (nacelle, γ, V). Thrust and CT are implicit: they follow from γ and V at a fixed weight. NICE-OPS uses (μ, γ, CT) with a Delaunay-barycentric blend over scattered conditions. AAM blends speed then FPA, with fall-backs to nearest, so it implicitly assumes a roughly grid-like (γ, V) layout. NICE-OPS's blend handles arbitrary scatter. Both blend linearly in dB.
- The dFPA/dSPD limits make AAM's output discontinuous in the operating condition. When a condition is too far from its neighbours it snaps to the nearest sphere, and on ties it picks the higher speed. NICE-OPS has no such distance gating.
- Q-SAM1/Q-SAM2 put manoeuvre effects into the lookup variable (γ_eff) and the sphere orientation. This is close to NICE-OPS's `fpa_e` (flight path angle relative to the TPP normal, `qsam` flag). NICE-OPS also carries CT through load factor. FRAME's "effective weight" plays the same role as NICE-OPS's CT = LF·W/(ρ V_tip² πR²) dimension.
- Weights are stored in each sphere, but outside FRAME/QSAM2 they are not used for lookup (see Gaps).

### Gaps
- The document does not say whether Q-SAM1's γ_eff is clipped, or how v̇ is computed for Q-SAM1 (finite differences on the segmented track?).
- The meaning and size of the "duration adjustment" in Table 2-8 is not explained.
- Internal inconsistency: Table 2-8's tie rule says "Approaches (FPA>0.0) Use the higher value of FPA; Departures (FPA<0.0) Use the lower". Elsewhere FPA is positive for climb/departure (Table 2-5 sorts "lowest, approach to highest, departure"; Table 4-37 says "positive = departure, negative = arrival") ([TR p.31, 34, 163][p34]). The signs in Table 2-8 look swapped. The code behaviour cannot be determined from the document.
- Thrust, CT, RPM, gross weight and density altitude are not used as rotorcraft lookup variables (only FRAME's effective weight). No rotorcraft analogue of the fixed-wing altitude correction is documented.

## 3. Sphere orientation in the flight frame

### Takeaway
The sphere is fixed to **body axes**, not to the tip-path plane or the velocity vector. The body frame is reached from a flight-path-fixed frame (heading β, flight path angle γ) by roll, pitch (= AOA + FPA) and yaw, which the user supplies per trajectory point. For rotorcraft, a constant MAST TILT can rotate the sphere nose-up, and Q-SAM1/2 rotate it using γ_eff and the derived roll. The legacy angle routine (ANGLES, unchanged since RNM 1.0) assumes zero roll and small angles. An iterative Newton Euler-angle solver switches in automatically above small thresholds, or when multiple sources are more than 10 ft apart.

### Cited Findings
- Orientation conventions: flight path angle (+ = climb), heading (N = 90 deg, E = 0 deg), pitch (+ = nose up), roll (+ = left down), yaw (+ = facing left) ([TR p.20 §2.3][p20]). Pitch is the "Summation of the angle of attack and the flight path angle". AOA is the "angle inscribed between the flight path and the X-axis of the body-fixed coordinate system" ([TR p.226–227 Glossary][p226]). In multiple-operations profiles, AOA is "the angle between the climb angle and the pitch angle", where "the pitch angle is the angle between the horizontal and the thrust vector" ([TR p.115 Table 3-64][p115]).
- Trajectory interpolation: kinematics assume constant acceleration between inputs. "Roll, pitch and yaw are interpolated linearly". Engine state is a step function per segment. Thrust-vector angle is interpolated linearly ([TR p.21 §2.3][p21]).
- θ is found from the dot product of the aircraft unit direction vector A⃗_c = cosτ cosξ î + cosτ sinξ ĵ + sinτ k̂ with the unit source-to-receiver vector, θ = cos⁻¹(A⃗_c·S⃗) (Eqs. 4–7). φ = cos⁻¹[(h_s/cosτ + S_dis cosθ tanτ)/(S_dis sinθ)] (Eq. 8), "where h_s is the height of the source above the ground level" ([TR p.24–25 §2.4.1][p24]; equation checked against the page image).
- "This geometric implementation (subroutine ANGLES) remains unchanged since RNM Version 1.0. It assumes a zero roll angle and makes small angle approximations for all angles except yaw and heading". A Newton-iteration Euler solver (after Gopalan 2001) is used automatically when |yaw| > 2 deg, |roll| > 5 deg, |pitch| > 2 deg, Σ > 3 deg, EULER ON, or when there is more than one sphere of a type with spacing above 10 ft. These thresholds "cannot be changed by the user" ([TR p.25 §2.4.1, Table 2-4][p25]). EULER OFF gives "the geometric calculations as performed by RNM4 and earlier" ([TR p.65][p65]).
- Sphere centres: they are body-fixed offsets (x, y, z)_body. They are converted to earth-fixed coordinates by a chain of nine 4×4 transforms (translate, rotate by −β, −γ, then ψ, θ, φ, back by γ, β, translate). "The model rigorously accounts for a full complement of vehicle orientation angles (roll, angle of attack, and yaw)" for source positions ([TR p.26–28 §2.4.2, Eqs. 9–12][p26]).
- MAST TILT: "Constant tilt angle uniformly added to source (.NC) orientation (+angle, nose up)", rotorcraft only. With ABSANGLES it is added to the user's AOA ([TR p.58, 69][p69]). ABSANGLES: "the geometry module bypasses any coordinate rotation and directly utilizes the user input earth-based coordinate system source orientation", and is allowed only for single-source vehicles ([TR p.59][p59]).
- Q-SAM1 rotates the source with γ_eff ([TR p.39][p39]). Q-SAM2 derives roll from the inertial tip-path-plane projection (Eqs. 17–18) ([TR p.40][p40]).

### Inferences
- AAM's sphere frame is (body x-axis, body y-axis), driven by user-supplied AOA, roll and yaw. NICE-OPS builds its frame from the TPP normal (the load-factor direction) and the velocity: `sphere_frame(velocity, normal)`, with the heading as the nose at hover. Q-SAM2's n̂ = (a + g)/|a + g| is exactly a load-factor TPP normal, so with QSAM2 the two models orient the sphere the same way. In plain AAM, the sphere tilts only through the AOA input.
- Eq. 8 uses the source height above ground h_s rather than the actual receiver vector. This suggests the first-order ANGLES path assumes a receiver at ground level on flat terrain. NICE-OPS computes φ from three dot products against an orthonormal frame, which is exact for any receiver height.

### Gaps
- τ and ξ are drawn in Figure 2-15 (τ the body-axis elevation, ε another angle) but never defined in the text. Presumably τ = pitch (γ + AOA) and ξ = heading plus yaw.
- There is no statement on how rotor tip-path-plane flapping or longitudinal cyclic is represented other than through AOA, MAST TILT or Q-SAM.

## 4. Within-sphere angular interpolation

### Takeaway
The technical reference does not say how AAM reads a level at an arbitrary (θ, φ) between grid nodes: nearest node, bilinear, or anything else. It documents only the ART-side gridding of scattered depropagated points onto the regular mesh (ZGRID, a Laplace/spline blend) and the pole treatments.

### Cited Findings
- AAM evaluates "the location on the source sphere intersected by a straight line between the source center and receiver" ([TR p.24 §2.4.1][p24]). Nothing more is said about lookup between nodes.
- Gridding inside ART: depropagated control points "are passed as random points into the ZGRID algorithm where the sphere mesh is populated". ZGRID gives "a user-selectable blend of LaPlace or Spline interpolation". LaPlace is recommended (CAY = 0), and NRNG sets how far each control point reaches ([TR p.146 §4.1.1][p146]; [TR p.145–146][p145]; GRIDOPTIONS [TR p.154][p154]).
- Poles: "At θ = 0 degrees and 180 degrees, where all φ represent the same physical point ... the average over all φ values is applied". "Values at θ adjacent to the next two columns of the 'nose' (bow) and 'tail' (stern) are filled via linear interpolation" ([TR p.146][p146]).

### Inferences
- From the 5 deg grid and the absence of any mention, it is likely (unverified) that AAM interpolates between grid nodes in some form, probably bilinear in dB, as RNM did. This needs confirmation from the code or the RNM manual (Page, Wilmer, Plotkin 2008, WR 08-04).
- ZGRID smooths the measured directivity on the sphere before propagation. A NICE-OPS user reading an AAM-made sphere inherits that smoothing.

### Gaps
- The interpolation method, its wrap-around behaviour at φ = ±180 deg, and how it treats a 1E35 neighbour are all undocumented here.

## 5. Doppler handling, convective amplification, emission and arrival time

### Takeaway
A sphere carries a flag for whether its spectra are in the moving (ground-observed) frame or the stationary aircraft frame (de-Dopplerized). AAM applies the frequency shift f′ = f/(1 + u/c). For moving-frame (flight-test) data it applies a *difference* correction between the test speed and the modelled speed. The ground model is evaluated at Doppler-shifted frequencies. No convective amplitude factor is mentioned anywhere. Arrival time comes from a straight-line path at the path-averaged sound speed. No retarded-time iteration is needed, because sound is "binned" at receivers by arrival time and then interpolated onto a uniform time mesh.

### Cited Findings
- Eq. 2: f′ = f·[1/(1 + u/c)]. "Unless data will always be applied for the same vehicle speed as it was collected, a Doppler correction will be needed. This can be accomplished within the AAM noise model, with the model using a difference based on Equation (2) applied to the flight test and the modeled airspeed." AAM can accept "de-dopplerized" data from aeroacoustic codes or wind tunnels, "particularly useful when there are tonal characteristics, as AAM outputs will report both the moving and stationary frequencies". The flag DOPPLERSTATE means 0 = stationary (aircraft system), 1 = moving (ground-based system) ([TR p.11 §2.2.1][p11]; page image checked).
- Ground reflection "uses the Doppler shifted frequencies based on the speed of the vehicle and its direction relative to the receiver" ([TR p.10 §2.2][p10]).
- In TERRAIN mode, "Doppler effect on frequency for stationary .NC file source acoustics not included" ([TR p.88 Table 3-39 note 1][p88]). With WEATHER (curved rays), "Acoustic source data must not have the Doppler shift removed" ([TR p.91][p91]).
- Doppler-shifted NB and PT frequencies that fall outside the 1/3-OB window are not binned (warnings ADDNB/ADDPT) ([TR p.141–142][p141]).
- ART's DOPPLERSTATE "is only set in the header and does NOT affect the processing of the measurement data" ([TR p.153 §4.1.2.6][p153]). AAM's ART keyword option STATIONARY makes AAM "compute the Doppler frequency shifts for itemization in the .TH file". ART then uses that file to "de-Dopplerize the microphone signal to aid in the generation of stationary source .NC files" ([TR p.61][p61]; [TR p.85][p85]). The B407 tutorial sphere was written with DOPPLERSTATE 1 (moving) ([TR p.206][p206]).
- Arrival time: "the arrival time is calculated based on a straight-line propagation path from the sound sphere center to the receiver using the average speed of sound on the path". For broadband and narrowband sources the vehicle centre is used ([TR p.40 §2.5.5][p40]). "Sound is 'binned' at the receiver locations based on the arrival time. Only after the propagation has been completed, are the receiver sounds interpolated onto a uniform time mesh and added" ([TR p.10][p10]). Absolute time and phase tracking "was implemented in order to eliminate the need for the iterative schemes found in early versions of the flight track module" ([TR p.41][p41]).

### Inferences
- NICE-OPS keeps spheres in the moving frame and applies no shift. That corresponds to AAM's DOPPLERSTATE = 1 used at the test speed. AAM additionally corrects for a mismatch between the modelled speed and the sphere's speed. NICE-OPS blends between speed-indexed spheres instead.
- Since AAM interpolates between spheres of different speeds, it is unclear whether the speed-difference correction is applied per sphere before blending or to the blend. The document does not say.
- Mixed conventions: Bullwink's input line 5 is "Doppler shift status (0=moving source, 1= stationary source)" ([TR p.165 Table 4-40][p165]). That is the opposite coding to DOPPLERSTATE (0 = stationary, 1 = moving). The NetCDF variable is `DOPPLER_SHIFT_REMOVED`, where 1 should mean stationary. So the integer in the file and the DOPPLERSTATE keyword are inverses of each other. Anyone writing AAM-compatible files from panam should store `DOPPLER_SHIFT_REMOVED`, not DOPPLERSTATE.

### Gaps
- u in Eq. 2 is undefined: its sign convention (presumably the source velocity component along the receiver direction, positive receding) and whether wind enters are not stated.
- No convective amplification (1 − M cosθ)^-n is documented, nor any dwell-time or retarded-time weighting in the integration. The only related remark is Table 2-8's unexplained "duration adjustment".

## 6. Propagation physics

### Takeaway
AAM's ground level is the straight-ray budget of Eq. 1: L = L(r₀) + A_spread + A_atm + A_grd + A_topo + A_wind + A_foli + A_nlin. Spreading is point-source 20 log(r/r₀). Absorption follows ANSI S1.26 (R2004) through a layered straight-ray atmosphere, with an SAE-AIR-1845 option. Over flat ground the excess ground attenuation uses Chien–Soroka/Chessell with Daigle's corrections and flow-resistivity impedance: incoherent over 1/3-OB bands, coherent for pure tones, with Chessell turbulent decoherence. Terrain uses Rasmussen's geometrical theory of diffraction with wedge/screen shielding and spatially varying impedance. Refraction and wind come from a separate pre-computed ray-tracing table (RNMwea), only over flat uniform terrain and never together with terrain. Foliage follows ISO 9613-2 (1/3-OB only), and Burgers-based nonlinear gain is for jets. Echoes, and a separate "lateral attenuation" term, are not modelled.

### Cited Findings
- Eq. 1 and its terms ([TR p.9–10 §2.2][p9]):
  - L(r₀) is the "Free field loss-less sound level at a distance r₀ from the source".
  - A_spread is "Geometrical spherical spreading loss (point source)".
  - A_atm is "ANSI/ISO atmospheric absorption standard [ANSI, 2004]".
  - A_grd is Chien & Soroka 1975 and Chessell 1977 with Daigle 1979 corrections, using Doppler-shifted frequencies.
  - A_topo is Rasmussen 1984. It "asymptote[s] to Maekawa's thin screen theory". "Echo effects ... are not currently treated". "Topography cannot currently be modeled simultaneously with A_wind".
  - Individual sources "are propagated independently".
- Spreading in practice: in the B407 .TH example, slant range 6837.41 ft and sphere radius 100 ft give Spread = −36.70 dB. That equals 20 log₁₀(6837.41/100) ([TR p.203 Figs. 5-28, 5-29][p203]; arithmetic is mine).
- Atmosphere:
  - The straight-ray module uses a "horizontally stratified temperature and pressure profile (no winds)", defaulting to the US Standard Atmosphere 1976 with constant 70% RH.
  - Sound speed and absorption coefficients are tabulated "at 1,000-foot intervals" as "weighted average[s]" by ANSI S1.26 (R2004).
  - Claimed accuracy is ±10% for 50–10,000 Hz, 0–40 °C and 10–100% RH.
  - "Absorption coefficients are computed based on the center band frequencies when Narrow Band, 1/12 Octave bands or pure tones/phase is in use."

  ([TR p.13–14 §2.2.3][p13]). ATMOS allows up to 32 levels and is incompatible with WEATHER. Above the profile's top the highest level is extrapolated ([TR p.61][p61]).
- SAE-AIR-1845 (ABS1845): coefficients are updated to CAEP7 values, extrapolated toward 0 below 50 Hz, isothermal 59 °F and 70% RH, 1/3-OB only, single-event only. Any user atmosphere is ignored ([TR p.14][p14]; [TR p.59][p59]).
- Flat-ground model: "broadband and narrowband signals are summed incoherently over one-third octave bands while pure-tone signals are added considering the single frequency phase shift due to the reflection". Chessell turbulent decoherence has "typical" value 2·10⁻⁴ rad·s/√ft, and 0 disables it ([TR p.11 §2.2.2][p11]; page image checked). SETUP PARA recommends 0.0004 ([TR p.75 Table 3-25 note 8][p75]). Flow-resistivity table (kPa·s/m²): snow 30, forest floor 50, grass 225, roadside dirt 650, packed sand 1650, hard dirt 3000, exposed dirt/rock 6000, concrete 10,000, asphalt 50,000, water 1,000,000 ([TR p.11 Table 2-1][p11]; [TR p.76 Table 3-26][p76]).
- Terrain:
  - A source-to-receiver "slice" through the terrain is classified as 2 points (flat), 3 points concave or convex, 4 points (wedge or screen with one flat), or 5 points (with two flats).
  - Shielding "transition[s] continuously from a Pierce wedge ... to a Maekawa screen". It is applied per 1/3-OB band.
  - Multiple reflections in valleys, ground impedance (per-cell .IMP file), turbulent decoherence and foliage are included. It reduces to the flat-earth model on flat ground.

  ([TR p.11–12 §2.2.2][p11]; TERRAIN keyword [TR p.88][p88]). The flow resistivity in SETUP PARA is ignored in TERRAIN mode ([TR p.75–76][p75]).
- Curved rays (RNMwea, since RNM 4.0):
  - Run as a pre-processor for "an omnidirectional point source over flat, uniform terrain", over a range of source heights. AAM requires at least two heights to interpolate ([TR p.174][p174]).
  - Rays are launched to reach a regular mesh, using PCBoom eikonal routines for a horizontally stratified medium.
  - Receiver-plane crossings, reflection points, incidence angles and times are gridded by Laplacian interpolation (PLOT88), then fed to a phase/incidence form of Chessell's Eq. 3.
  - Impedance is Delany & Bazley.
  - Downwind multiple reflections use a third receiver-plane crossing and Embleton's method; "Multiple downwind paths are not addressed".
  - Shadow-zone diffraction follows Berry & Daigle 1988, and absorption ANSI S1.26-2014.
  - Results go to a .RAZ table looked up by geometry, "only ... for single event analyses over uniform flat terrain".

  ([TR p.14–16 §2.2.4][p14]; [TR p.172–175 §4.16][p172]). WEATHER is limited to broadband and narrowband spheres ([TR p.91][p91]). 3D_WEATHER reads an external 3-D time-dependent RAZ file from NASA's SAM code (APET) ([TR p.92][p92]).
- Turbulence: beyond Chessell decoherence of direct and reflected rays, the text only notes qualitatively that "Turbulence can also 'fill in' a shadow zone by scattering sound into it" in describing the physics ([TR p.14][p14]).
- Foliage: ISO 9613-2, for dense foliage only, with d_f = d₁ + d₂ along 15 deg arcs of radius 5 km. Octave-band attenuations are applied to their 1/3-OB constituents. It is not applied when the vehicle is above the trees, and is 1/3-OB only ([TR p.17–18 §2.2.5][p17]).
- Nonlinear propagation (jets) uses Gee's generalized Burgers equation, either from a waveform embedded in the sphere (direction-independent gain) or from tables of gain against directivity, distance, temperature, RH and frequency ([TR p.1][p1]; [TR p.18–19 §2.2.6][p18]).
- Computational culling:
  - A track is skipped for a receiver if the slant range at the point of closest approach exceeds the cutoff (default 30,000 ft).
  - Individual points are skipped beyond R_max = R_min·10^(Δβ/20), where Δβ (default 40 dB, so 100×) is the "spherical loss factor".
  - CALCALL disables both checks.

  ([TR p.75–76][p75]).

### Inferences
- AAM's broadband path is fully spectral: every band is propagated with its own absorption and ground effect. This corresponds to NICE-OPS's `--spectral --ground_effect` path, not to its default broadband EAA path. AAM has no equivalent of the broadband dBA plus fitted EAA.
- AAM's refraction is a pre-tabulated omnidirectional-source correction on a ground mesh, indexed by source height, and is incompatible with terrain. NICE-OPS's `--refraction linear|stratified` traces per pair, with ground models, shadow transitions and duct arrivals. That is substantially richer, and AAM has nothing comparable to the PE reference.
- No lateral-attenuation term (such as SAE AIR 5662) appears; ground and terrain physics carry that role. No receiver-side plate or board model is mentioned.

### Gaps
- Whether 1/3-OB absorption is evaluated at the band centre or integrated over the band. The document states band centres only for NB, 1/12-OB and PT.
- The A_spread reference is not stated as an equation. Only the .TH example (above) shows it as 20 log(r/R_sphere).
- Two ANSI S1.26 editions are cited (R2004 for straight rays, 2014 for RNMwea), and the document does not reconcile them.
- No spectral broadening, no turbulence scattering model, and no urban or building reflections beyond the thin screens used for terrain features.

## 7. Receiver modeling

### Takeaway
Receivers are a rectangular grid or up to 500 points of interest (POI). A grid sits at one height above ground level, normally 5 ft, with spacing of at least 100 ft and at least 5 points per axis. Each POI carries its own z AGL. Nested subgrids (new in v3.1.0) double the resolution locally. HOTSPOT gives a radial mesh from a single vehicle position. Receiver height enters the ground and terrain models. There is no plate or microphone-mounting model on the prediction side; ART applies incidence-angle corrections to measured data only.

### Cited Findings
- Grid spacing ≥ 100 ft. Ground elevation is MSL. "All points on the grid must be at the same elevation AGL. Normally the grid elevation is set for 5 feet AGL" ([TR p.75 SETUP PARA][p75]). Each direction needs at least 5 points ([TR p.54–55 §3.5.1.1][p54]).
- POI: up to 500 points, x and y in grid coordinates, z in ft AGL (from the .ELV terrain in TERRAIN mode) ([TR p.73][p73]).
- SUBGRIDS (v3.1.0): each subgrid has twice the resolution of its parent. NI and NJ must be odd and at least 3. Shared points are not recomputed ([TR p.86–88][p86]; [ReadMe][readme]).
- HOTSPOT: one vehicle state; points of interest on a radial grid with a user-set receiver height. The output .HS file is structured like a .TIG file, for detection studies ([TR p.66–67][p66]).
- Terrain, impedance and foliage grids "must ... encompass the entire computational domain", and must extend far enough that levels are at least 10 dB down at their edges ([TR p.53 §3.4][p53]).
- Multiple-operations profiles: altitudes below 5 ft AGL are raised to 5 ft "For compatibility and consistency with NMAP" ([TR p.114][p114]).

### Inferences
- AAM's grid receivers all share one height, as NICE-OPS's grids effectively do. NICE-OPS additionally supports per-observer heights, plates and point receivers in one run. AAM's POIs allow a per-point z but no receiver-type model.

### Gaps
- How receiver height enters Chien–Soroka (point receiver at z AGL is assumed).
- Whether grid points above sloping terrain lean with it. The document says only "at a specified height above the local terrain" ([TR p.9][p9]).

## 8. Time integration and metrics

### Takeaway
AAM is a time-stepping simulation:
- the trajectory is segmented at a default 2.0 s (TIMESPACING overrides it; turns are sub-segmented, by default every 10 deg);
- each emission is propagated and binned by arrival time, then interpolated onto a uniform receiver time mesh and summed;
- metrics are computed from the resulting spectral time histories.

SEL integrates the *whole* time history by default; DBDOWN10 restricts it to within 10 dB of the maximum. Narrowband, 1/12-OB and pure-tone signals are binned or rolled up into 1/3-OB bands for the standard metrics. PNL, PNLT, EPNL and WECPNL need 1/3-OB sources. PEARLS mixes coarse "macro" steps with bursts of fine "micro" steps for time-varying loudness.

### Cited Findings
- Default 2.0 s spacing, overridable in single-event mode ([TR p.9][p9]; [TR p.20 §2.3][p20]; [TR p.88 TIMESPACING][p88]). Turns are approximated by an even number of straight segments with a sub-angle of at most 30 deg (Eldred & Miller), "This may result in smaller time steps" ([TR p.22–23][p22]). TURNSPACING defaults to 10 deg ([TR p.89][p89]). The ART tutorial recommends a nominal 0.5 s and its example uses 0.2 s ([TR p.199–200 §5.8][p199]).
- Binning by arrival time, then interpolation onto a uniform time mesh and summation ([TR p.10][p10]).
- Metrics:
  - Lmax, SEL (overall, C- and A-weighted), DNL, CNEL, NEF, PNL, PNLT, EPNL, WECPNL, DPRIME, NA, TA, TVL and PA ([TR p.3][p3]; [TR p.50–52 §3.3][p50]).
  - CNEL uses SEL, not SENEL, because "SENEL is not currently supported" ([TR p.51][p51]).
  - PNL is sampled "twice per second" over 24 bands from 50 Hz to 10 kHz. PNL, PNLT, EPNL and WECPNL are computed "only if the noise spheres contain one-third octave band (OTOB) data" ([TR p.51–52][p51]).
  - TVL comes from the TVLInputs post-processor (Moore–Glasberg/Schlittenlacher model) at 20 ms or at AAM's time resolution ([TR p.52][p52]; [TR p.177][p177]).
- "DBDOWN10 forces all integrated metric calculations to only include points within 10 dB of the maximum level. The default without this keyword is to include all frequencies for the full the time history" ([TR p.64][p64]). USERBAND restricts the metric bands to a subset of 10–40 ([TR p.89][p89]).
- NB and PT signals are binned into 1/3-OB bands by ADDNB/ADDPT ([TR p.141–142][p141]). For 1/12-OB grids, "the 1/12 OB data has been rolled up into the 1/3 OBs for computation of the standard metrics" ([TR p.225][p225]).
- Outputs:
  - .GRD (NMPlot) and .PLT grids.
  - .POI/.TH/.NB/.PT spectral time histories with per-step propagation terms: Spread, AirAb, gOut, Nlc, FoliAttn, θ, φ, slant range, emission and arrival times.
  - .TIG for video and audibility.
  - NCTHIST NetCDF time histories (needed for 1/12-OB detail), with spectral levels against FREQUENCY and TIME.

  ([TR p.44 Table 3-1][p44]; [TR p.69–70][p69]; [TR p.84–86][p84]).
- PEARLS: macro steps of typically 0.2–1.0 s, micro steps of 10⁻² to 10⁻⁴ s, a set pearl width and interval, with integer-multiple rules. The total point limit is "4000 points" ([TR p.72–73][p72]). The error-message section gives 8,000 ([TR p.137][p137]).
- TA needs a time spacing under 2 s, or AAM warns ([TR p.142][p142]).

### Inferences
- NICE-OPS integrates SEL only over the 10-dB-down window (the first-crossing `peak_range`), which corresponds to AAM's opt-in DBDOWN10, not to its default. AAM's default sums the full history. A like-for-like comparison of the two codes needs DBDOWN10 on in AAM, or NICE-OPS's window widened. AAM's DBDOWN10 window definition (first or outermost crossing) is not stated.
- AAM's 2 s default step is far coarser than NICE-OPS's 0.02 s tracks. The document's footnote argues the coarse step suits 500 ft community grids ([TR p.20 footnote][p20]).

### Gaps
- The SEL quadrature rule (rectangle or trapezoid on the uniform receiver mesh), the uniform mesh's step, and the interpolation used for it (energy-linear?) are not given for AAM itself. Only the visualizer is described as interpolating linearly in energy ([TR p.176][p176]).

## 9. Multiple sources and source offsets

### Takeaway
A vehicle can carry up to 10 broadband, 10 narrowband and 10 pure-tone spheres at body-fixed offsets, for example main rotor, tail rotor, engine and MR/TR wake interaction. Each is propagated independently. Broadband and narrowband sum incoherently and take arrival time from the vehicle centre. Pure tones sum coherently, with absolute phase tracked from the start of the track. All sources of a vehicle share one selected nacelle angle and condition set. Formation flight (multiple vehicles) is also supported.

### Cited Findings
- "For each set of flight conditions there can be up to ten broadband and ten pure-tone sound spheres centered at different points on the vehicle (i.e., engine, main rotor (MR), tail rotor (TR), and MR/TR wake interference)" ([TR p.30][p30]). The count is 10 BB, 10 NB and 10 PT in SETUP PARA ([TR p.75–76][p75]). Helicopter and fixed-wing sources "cannot be intermixed for a given vehicle" ([TR p.30][p30]).
- Broadband and narrowband arrival time comes from the vehicle centre, "since phase is not a concern with broadband (incoherent) noise signals". This assumes source separation is small relative to path length. For pure tones, phasing and small path differences are "critical". Phase is referenced to the first trajectory time, p(t) = Σ Aₙ e^{i(ωₙt + θₙ)}, and |p|² = Re² + Im² ([TR p.40–41 §2.5.5, Eqs. 20–21][p40]).
- Sphere positions are transformed through the full roll, AOA and yaw rotation ([TR p.28][p28]). Multiple spheres more than 10 ft apart trigger the iterative Euler solver ([TR p.25][p25]).
- "Multiple vehicles, such as those in formation flight, may also be modeled" ([TR p.3][p3]).

### Inferences
- NICE-OPS uses one compact sphere per condition, centred at the hub. AAM's multi-sphere option allows separate MR/TR/engine components and coherent tonal interference. Neither the shipped CH-146 nor the Bo105 sets are documented as using more than one sphere, and the B407 example uses one broadband sphere at (0, 0, 0) ([TR p.206][p206]).

### Gaps
- Whether interpolation across conditions is done per component sphere independently (it presumably is, since search is keyed on matching XYZ centre).

## 10. Stated limitations, assumptions, validation, known issues and version differences

### Takeaway
The document lists many functional restrictions:
- no echoes;
- terrain and wind cannot be combined;
- curved rays only over flat terrain in single-event mode;
- a legacy small-angle, zero-roll angle routine;
- no wind in source lookup (ground speed is treated as airspeed);
- no extrapolation in fixed-wing power;
- several metrics and features 1/3-OB-only;
- 1E35 holes are not propagated.

Validation is cited only qualitatively ("agree well") against NATO F-16 terrain tests and RNM/NMSim comparisons. No error statistics are given. The document also contradicts itself in a handful of places.

### Cited Findings
- Limitations:
  - Echo effects are not treated ([TR p.3][p3]; [TR p.10][p10]).
  - Topography and wind are not simultaneous ([TR p.10][p10]).
  - Curved rays are flat-terrain, single-event, BB/NB only ([TR p.14][p14]; [TR p.91][p91]).
  - The ANGLES assumptions are as in §3 ([TR p.25][p25]).
  - No wind in airspeed or source lookup ([TR p.226–227][p226]).
  - Track dispersion exists only in multiple-operations mode and is Gaussian across a width of ±3σ ([TR p.41 §2.5.6][p41]; [TR p.115][p115]).
  - Foliage, PNL-family metrics and ABS1845 are 1/3-OB only ([TR p.17, 51–52, 59][p17]).
  - Fixed-wing extrapolation is not allowed ([TR p.37–38][p37]).
  - AAM v3.0.0 does not accommodate multi-POI NCTHIST files ([TR p.70][p70]).
  - Windows-only (64-bit), built with Intel Fortran. HIGHMEM and MEDMEM can crash on allocation ([TR p.42][p42]; [TR p.45][p45]).
- Validation claims:
  - Algorithms derive from ANOPP and "have been validated through many years of testing".
  - RNM, applied to the 1994 NATO F-16 Narvik departure and mountain flyover, "agree[s] well" with measurement [Lee et al. 1996; Page 2002]. NMSim was validated there and at Grand Canyon [Miller 2003].
  - "The calculation and aircraft type capabilities of AAM are a superset of those in RNM and NMSim".

  ([TR p.3][p3]; [TR p.12][p12]). No quantitative error metrics are given anywhere in the document.
- Version history:
  - RNM (late 1990s, NASA/Wyle, for the XV-15) gained helicopters, then ART; the Navy funded terrain in the early 2000s.
  - AAM added fixed-wing under SERDP SI-1304.
  - The Quarry Noise Model and foliage were added under ODOT.
  - 1/12-OB and PEARLS came from Uber Elevate.

  ([TR p.4–5 §1.2][p4]). NASA transferred AAM to Volpe in spring 2020 ([pAnn]).
- v3.1.0 changes ([ReadMe][readme]):
  - SUBGRIDS for multi-operations and single-event modes (with the beta BaseOps);
  - mixed BB/NB custom spheres with improved NC inventory, letting BB and NB files coexist regardless of SETUP PARA;
  - improved checking of NC files "for compliance with AAM requirements, principally for UAM/eVTOL single event analyses";
  - PEARLS fixes "including fixes to the array buffer sizes";
  - standardised error messages;
  - an updated QuietSphere;
  - a revised manual.

  The ReadMe is signed "Juliet Page, Volpe ... 27 May 2010", which looks like a typo for 2021.
- Earlier-version items noted in the document:
  - RNM v1 interpolated in airspeed only ([TR p.31][p31]).
  - RNM 4.0 introduced curved rays ([TR p.14][p14]).
  - RNM v7 began rounding FPA to 0.5 deg ([TR p.31][p31]).
  - EULER OFF reproduces RNM4 geometry ([TR p.65][p65]).
  - ART3 (linear interpolation) was replaced by ART5's ZGRID and gained a high-frequency absorption cap ([TR p.146][p146]).
  - Older ART extrapolated missing regions; current ART leaves 1E35 holes ([TR p.145][p145]).
  - The memory model dates from v2.0 ([TR p.66][p66]).
- Internal inconsistencies:
  - PEARLS point limit: 4000 ([TR p.72–73][p72]) against 8000 ([TR p.137][p137]).
  - Turbulence parameter: 2·10⁻⁴ "typical" ([TR p.11][p11]) against 0.0004 "recommended" ([TR p.75][p75]).
  - Table 4-37 says phi runs "-90° to +90°" while full spheres use ±180 ([TR p.163][p163] against [p204]).
  - Bullwink's Doppler coding is the inverse of DOPPLERSTATE ([TR p.165][p165] against [TR p.11][p11]).
  - The Table 2-8 approach/departure sign conflict (§2).
  - The glossary defines body X "from the bow to the stern" ([TR p.226][p226]), while θ = 0 is toward the nose.
  - The flight-path-fixed glossary has "Y is positive on the port side" ([TR p.227][p227]) while body Y is starboard.

### Inferences
- For a NICE-OPS comparison, the most consequential AAM assumptions are:
  1. lookup by ground-referenced speed;
  2. FPA/speed gating with nearest-sphere fall-back;
  3. the full-history SEL default;
  4. spectral (not broadband-EAA) propagation;
  5. ART's ZGRID smoothing and 1E35 holes in measured spheres.

### Gaps
- No quantitative validation of AAM v3 rotorcraft predictions against measurement appears in this document.

## 11. Sphere generation: ART depropagation, ground-reflection removal, ambient handling, analytical spheres

### Takeaway
Measured spheres are made in two steps:
1. AAM runs the measured trajectory with a 200 dB omnidirectional sphere (`OMNI_200.NC`) and writes a .TH file with every propagation term (spreading, absorption, ground or terrain, weather), θ, φ and arrival time per microphone and step.
2. ART interpolates each microphone's spectral time history to those arrival times, subtracts the propagation terms band by band, keeps only bands above a signal-to-noise threshold (no ambient subtraction by default), and grids the scattered (θ, φ) control points onto the sphere mesh with ZGRID.

Ground reflection is removed exactly as AAM's own ground model predicts it, using per-microphone flow resistivity, terrain or curved rays. Undefined regions stay at 1E35 unless the user mirrors, extends or smooths them. Analytical or synthetic spheres are written with Bullwink.

### Cited Findings
- "Depropagation": AAM computes spreading, absorption, ground effect and source–receiver geometry; "ART subtracts these from measured noise levels yielding the source characteristics ... levels on a spherical set of (θ, φ) are obtained" ([TR p.3 §1.1][p3]). "ART contains no propagation physics except in the case of fixed-wing aircraft with non-linear effects" ([TR p.145][p145]).
- Procedure ([TR p.146–147 §4.1.1][p146]):
  - An omnidirectional source is propagated and the .TH file written.
  - "The measured microphone data is interpolated to the exact arrival time."
  - Corrections "are applied spectrally to the microphone data", only for "microphone data that exceed the ambient level by the specified signal-to-noise ratio".
  - The control points go into ZGRID.
  - Poles are averaged over φ, and the adjacent θ columns linearly interpolated.
  - North Pole (in-plane forward, vertical-array) microphones are treated like the lateral array, with energy averaging inside the NPGRIDAREAS regions.
- OMNI_200: AAM auto-creates a 200 dB omni sphere (radius 100 ft by default) when the ART keyword is used, and replaces it if the spectral content differs ([TR p.60–61][p60]). The tutorial calls it "a 200 dB omni directional full sphere with 5 degree grid resolution" ([TR p.204][p204]).
- Ambient: SNRATIO is a minimum signal-to-noise in cB, default 50 cB (5 dB). "There is NO subtraction of ambient noise from the signal" ([TR p.159 §4.1.2.29][p159]).
  - USESYNSPL instead applies OMEGA 5 logic: unchanged if at least 16 dB above ambient; energy-subtract the ambient when 6–16 dB above; delete when 6 dB or less. The whole spectrum is dropped if more than 14 bands are deleted. Low bands are filled flat, high bands by the last slope or −6 dB/band (whichever is steeper), interior gaps linearly ([TR p.161 §4.1.2.37][p161]).
  - MINAMBIENT uses the minimum per band as the ambient ([TR p.155][p155]).
- High-frequency absorption cap (default since ART 5): the absorption at 5 kHz (band 37) is used for all higher bands "to prevent unrealistically high depropagated source levels". NOABSCAP disables it ([TR p.146][p146]). ABSOFF turns off absorption removal entirely ([TR p.150][p150]). OMEGASLOPE caps levels above a band to a set dB/band slope ([TR p.158][p158]).
- Ground removal: each microphone in the ART MICROPHONE keyword carries a flow resistivity ("1E6 = hard, 200 = soft") ([TR p.155 Table 4-17][p155]). The AAM run may use "terrain, non-uniform ground impedance, and curved rays" to "correct measurement data" ([TR p.145][p145]). In the B407 example, the AAM run uses TERRAIN with flat100k.elv/softer100k.imp, flow resistivity 200, turbulence 0.0004, the measured ATMOS, TIMETRAJ GPS track with AOA and roll, and TIMESPACING 0.2 s ([TR p.200 Fig. 5-26][p200]). The .TH example shows gOut of +5.23 dB at 10 Hz falling to −17.6 dB at 1 kHz for a ground microphone at 6837 ft slant range ([TR p.203][p203]).
- Microphone response: INCIDENCE and ADJUSTBYINCIDENCE apply free-field corrections by angle of incidence for vertically oriented, skyward-facing microphones ([TR p.151, 154][p151]). MICADJUST adds per-band level corrections ([TR p.155][p155]).
- Sphere controls:
  - SPHERE sets θ/φ min, max and step, and the maximum slant range to depropagate ([TR p.159][p159]).
  - SYMMETRICAL averages port and starboard in dB ([TR p.159][p159]).
  - EXTEND45/EXTENDX copy the φ = 45 deg (or user-chosen) values to the top of the domain ([TR p.153][p153]).
  - WRAPAROUND mirrors the lower hemisphere upward ([TR p.161][p161]).
  - TAPER, NOSETAILEXT and PHISMOOTH (N-point smoothing in φ per θ) are further options ([TR p.157–160][p157]).
  - CTLOVERRIDE combines depropagated data from several runs (1/3-OB only) ([TR p.152][p152]).
- Assigned condition: ART assigns the average airspeed and FPA of the trimmed trajectory to the header (the B407 log gives 102.52 kt, 0.286 deg). SPEED and FPA keywords override this ([TR p.207 Fig. 5-34][p207]; [TR p.153, 159][p153]). The trajectory "should be pre-trimmed to only include those portions of the run which are considered 'on condition'". Its time spacing "ha[s] a strong influence on the grid spacing and resolution on the lower portion of the sphere" ([TR p.160 THFILE][p160]).
- The tutorial data are the B407 run 109 from the Amedee 2017 tests (condition L4, 100 kt, FPA 0, entry altitude 200 ft). It used 49 microphones, 31 1/3-OB bands, a 5 deg sphere from θ 0–180 and φ −90 to 90, slant max 3000 ft, NRNG 7, radius 100 ft and DOPPLERSTATE 1 ([TR p.199–207 §5.8][p199]). The flight test uses a "U"-shaped array of ground and elevated microphones; without elevated microphones "That limits the θ measurement range" ([TR p.198][p198]).
- 1/12-OB measured data: `ProcessNCpascalto12.py` FFTs pressure in 0.5 s Hamming-windowed blocks at 2 Hz resolution, then groups into 121 1/12-OB bands ([TR p.170 §4.11][p170]).
- Analytical spheres: "may be obtained from flight tests, wind tunnel measurements, or theoretical predictions" ([TR p.1][p1]). Examples:
  - Bullwink writes BB/NB/PT NetCDF from ASCII ([TR p.164–166 §4.3][p164]).
  - The eVTOL spheres are ANOPP-PAS-based omni 1/12-OB ([TR p.213][p213]).
  - The B737 samples were made "from publically available noise-power-distance data with assumed source spectral and directivity information" ([TR p.47 Table 3-4][p47]).
  - Shipped rotorcraft are the Bo105 (20 spheres) and CH-146 (Bell 407 variant, 41 spheres), both from ART ([TR p.47][p47]).
  - FRAME produces sphere databases from flight-test-derived models ([TR p.41][p41]).

### Inferences
- ART's ground removal uses the same Chien–Soroka/Chessell (or Rasmussen) model and the same nominal flow resistivity that AAM later applies forward. Errors in the ground model therefore partly cancel when an ART sphere is re-propagated with matching settings in AAM. Spheres propagated in a different code (such as NICE-OPS's panam-ported EGA, or a board/plate receiver model) will not get that cancellation, and the mismatch goes straight into predicted levels.
- AAM's ambient handling is gating, not subtraction (unless USESYNSPL). Depropagated levels near the noise floor are biased high by up to about 1.2 dB at the 5 dB threshold. That figure is my arithmetic: 10 log(1 + 10^-0.5) ≈ 1.19 dB. Holes are left as 1E35 instead of being filled.
- ART sets a sphere's (V, γ) to the averaged on-condition values, and AAM rounds γ to 0.5 deg. A sphere set from many runs therefore forms a scattered (V, γ) cloud, not a grid, which is exactly the case AAM's bracketing search handles worst (§2).

### Gaps
- ZGRID's algorithm details (beyond "LaPlace or Spline", NRNG and CAY) are in an external PLOT88 manual, not in this document.
- The document gives no treatment of combining multiple microphone heights to remove ground interference (for example, using elevated microphones to separate direct and reflected paths). Ground removal is purely model-based.

[pAnn]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=3
[readme]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/ReadMe-AAM-v3.1.0.txt
[p1]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=17
[p2]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=18
[p3]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=19
[p4]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=20
[p8]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=24
[p9]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=25
[p10]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=26
[p11]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=27
[p12]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=28
[p13]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=29
[p14]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=30
[p16]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=32
[p17]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=33
[p18]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=34
[p20]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=36
[p21]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=37
[p22]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=38
[p24]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=40
[p25]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=41
[p26]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=42
[p28]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=44
[p29]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=45
[p30]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=46
[p31]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=47
[p33]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=49
[p34]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=50
[p37]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=53
[p39]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=55
[p40]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=56
[p41]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=57
[p42]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=58
[p44]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=60
[p45]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=61
[p46]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=62
[p47]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=63
[p50]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=66
[p51]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=67
[p52]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=68
[p53]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=69
[p54]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=70
[p59]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=75
[p60]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=76
[p61]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=77
[p64]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=80
[p65]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=81
[p66]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=82
[p69]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=85
[p70]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=86
[p72]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=88
[p73]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=89
[p74]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=90
[p75]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=91
[p76]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=92
[p84]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=100
[p85]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=101
[p86]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=102
[p88]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=104
[p89]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=105
[p91]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=107
[p92]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=108
[p114]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=130
[p115]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=131
[p134]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=150
[p137]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=153
[p139]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=155
[p142]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=158
[p144]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=160
[p145]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=161
[p146]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=162
[p147]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=163
[p150]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=166
[p151]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=167
[p152]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=168
[p153]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=169
[p154]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=170
[p155]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=171
[p157]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=173
[p158]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=174
[p159]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=175
[p160]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=176
[p161]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=177
[p162]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=178
[p163]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=179
[p164]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=180
[p165]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=181
[p166]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=182
[p168]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=184
[p170]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=186
[p172]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=188
[p174]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=190
[p175]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=191
[p176]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=192
[p177]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=193
[p198]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=214
[p199]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=215
[p200]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=216
[p203]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=219
[p204]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=220
[p206]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=222
[p207]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=223
[p213]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=229
[p214]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=230
[p225]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=241
[p226]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=242
[p227]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=243
[p141]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=157
[p138]: file:///Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/AAMv3/AAM-V3.10Beta-Distribution-May-2021/AAM_v3_TechnicalReferenceAndUsersGuide-27May2021.pdf#page=154
