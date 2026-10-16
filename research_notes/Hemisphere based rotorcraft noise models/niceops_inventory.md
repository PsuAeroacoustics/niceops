# NICEOPS and panam: what the model actually does (code-verified inventory)

Read-only inventory, compiled 2026-10-02 from code and documents. Roots used in the citations below:
**N** = `/Users/evg5332/Code/D/NICEOPS` (checkout on branch `stratified-ducts`, HEAD `e472b6f`; the review-fix branch `review-fixes-2026-10-02` in `/Users/evg5332/Code/D/NICEOPS-review-fixes` is NOT merged, so behaviour described is HEAD's);
**P** = panam as read from the validation harness's public-lineage clone `/Users/evg5332/Code/Python/niceops-validation/external/panam` (branch `stratified-ducts`, `83f989f`). This is the lineage that built the shipped B407 databases. The iCloud checkout `~/Documents/Code/Python/panam` was not read: it is in iCloud (eviction risk) and appears to carry the private `beta26` work, which must not be summarised;
**H** = validation harness `/Users/evg5332/Code/Python/niceops-validation`;
**M** = user memory `/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/`;
**R** = `N/reports/` (code review and grazing-incidence reports, untracked).
Where a doc claim disagrees with code it is flagged "doc/code conflict".

## 1. Source representation (sphere format, grid, conventions, hemisphere handling)

### Takeaway
NICEOPS reads a netCDF `.nod` database: one group per operating condition, each a 30.48 m (100 ft) radius sphere on a regular 36 (phi, 10 deg) x 37 (theta, 5 deg) grid carrying broadband dBA plus an EAA attenuation rate, and optionally a 31-band one-third-octave spectrum (10 Hz to 10 kHz). Only the lower hemisphere is measured. panam mirrors it to the upper surface, but NICEOPS's lookup folds every direction into the lower hemisphere anyway (|up|), so the upper half is never read.

### Cited Findings
- Each `.nod` group holds the scalars `advance_ratio`, `flight_path_angle`, `thrust_coefficient`, `radius` and `rotor_scale`, plus per-channel `dBA` and `EAA` and optional `frequency` + `amplitude(PHI,THETA,frequency)`. Root holds `same_grid` (required: `same_grid == false` is refused, the interpolation path is unimplemented), `database_version`, `fixed_load_factor`, `shared_grid_and_frequency`, rotor radius, tip speed and weight — [N/src/sphere_database.cpp:764-781, 856-1000](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); confirmed by `ncdump -h` of `N/databases/B407_fixed_load_factor.nod` and `S-76D_M3.nod`.
- The sphere radius is radius x rotor_scale and must be constant across spheres. Every shipped database checked (B407, Be206L3, AS350, S-76D) has 30.48 m, from panam's `r_ref_ft = 100` — [N/src/sphere_database.cpp:876-883, 1140-1143](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [P/noise_abatement_2017.py:97](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py); netCDF summary script run on `N/databases/*.nod`.
- Grid: 36 unique phi (-180 to 170, 10 deg step) x 37 theta (0 to 180, 5 deg step) = 1332 channels (B407 and the other 2017 aircraft) or 1406 (S-76D_M3). B407_fixed_load_factor: 179 spheres, 31 bands, 10 Hz to 10 kHz — netCDF read of `N/databases/B407_fixed_load_factor.nod`.
- Angle convention (AAM/RNM "ART"): theta is measured from the nose (projected velocity in the tip-path plane), 0 ahead and 180 aft. phi is the roll angle about the nose axis, with phi = 0 straight down along minus the TPP normal and phi > 0 starboard. The code computes phi = atan(across / |up|), which is confined to [-90, 90] — [N/src/sphere_database.cpp:678-734](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); the code review confirms "phi > 0 = starboard matches panam and the B407 advancing-side lobe" — [R/NICEOPS code review 2026-10-02.md, sec. 10 D](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Upper hemisphere: panam's source spheres cover only phi in [-90, 90] (lower half). `mirror_phi_to_upper_surface` reflects through the horizontal plane (phi -> 180 - phi) — [P/flight_acoustics.py:3073-3136](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py). NICEOPS's `phi_degrees` uses `std::fabs(up)`, so a direction above the TPP is read at its mirror image below. The reviewer confirmed the stored upper half is an exact mirror (0.0000 dB over 703 cells), "so the |up| fold loses nothing" — [N/src/sphere_database.cpp:711-724](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [R/code review sec. 10 D](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Legacy repair: databases with `database_version < 1` had hover spheres completed with `-phi` instead of `180 - phi`. `repair_legacy_hover_directivity` relabels the mirrored half (phi_correct = wrap(180 + phi_stored)) and refuses any sphere that does not match the broken shape. In S-76D_M3 that is the 24 mu = 0 spheres — [N/src/sphere_database.cpp:260-327, 920](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [N/CLAUDE.md "Legacy databases"](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md).
- Two propagation paths: broadband (dBA - spreading - EAA·(r - R)/1000) and spectral (per band: L_b - spreading - alpha_b(r - R), then A-weighted energy sum). **Spectral is the default** (`bool spectral = true`, flag `--spectral,!--no-spectral`) and falls back to broadband when the database has no spectra — [N/src/niceops.cpp:551, 1010-1034, 2357-2361](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp). Doc/code conflict: README's table and CLAUDE.md call `--spectral` opt-in; the code review flags this as stale (item X5) — [N/README.md:168, 279-291](/Users/evg5332/Code/D/NICEOPS/README.md); [R/code review sec. 10 D X5](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Interpolation is in dB in both paths: barycentric across conditions, then bilinear across the four grid corners. A corner or condition with no measurement ("empty", floored at `no_energy_floor_db`) drops out of the blend rather than pulling it toward silence — [N/src/sphere_database.cpp:1281-1320, 1395-1450](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp).
- Shipped databases (all git-tracked in `N/databases/`): AS350, AW139, B205, B407, Be206L3, EC130, RO-44, RO-66 and S-76D, in `*_fixed_load_factor.nod` (spectral, CT analytic) and `*_ambient_gated.nod` / `*_M*.nod` (broadband, CT stored) variants. Doc/code conflict: README says "Databases are not tracked in this repository", but `git ls-files databases` lists them — [N/README.md:258-259](/Users/evg5332/Code/D/NICEOPS/README.md).
- Source sphere data on OneDrive `Data/Sphere_Database` (listed only): Be206L3 94 .nc (+ LowNR 9), Be407 113, Be430 87, EC130 107, MD-902 84, Mi-8 85, R-44 59, R-66 82, SH-60B 16; AS350 in subfolders Amedee 257, B3 95, SaltonSea 189, Sweetwater 166; plus `Be407_spectral.nod`. `Data/Sphere_Database_Rebuilt` holds the panam rebuilds: AS350 96, Be206L3 122, Be407 118, EC130 38, RO-44 50, RO-66 89 — directory listing of `/Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Data/Sphere_Database*`.

### Inferences
- In hemisphere-tool terms NICEOPS is a lower-hemisphere model, like AAM/NORAH2. The full-sphere storage is cosmetic, because the lookup never reads the upper half. An observer above the TPP (elevated receiver, steep bank, terrain above the aircraft) gets the mirrored lower-hemisphere level.
- The 10 deg phi step means that abeam, the elevation resolution near the horizon is 10 deg. The harness findings name the AAM output grid plus bilinear reading as a coarsening stage at the rim, though measured to be within about 0.3 dB of the samples at 4-20 deg (see Q10).

### Gaps
- No database in the repo records `speed_reference`, and none of the netCDF headers checked carries global attributes. The build atmosphere behind EAA is not recorded in the file either (code review X3).

## 2. Operating-state parameterisation and interpolation across conditions

### Takeaway
Conditions are (mu, gamma, CT). mu = |v|/V_tip uses total speed; gamma is by default the QSAM "effective" path angle relative to the TPP; CT = LF·W/(rho·V_tip²·pi·R²). Each state is blended by barycentric weights from an exact-predicate Delaunay triangulation of the normalised conditions. Queries are clamped to the envelope, with boundary points pushed out 10% so the edges stay inside the hull, and fall back to the nearest condition outside the hull. Manoeuvres enter only through the load factor (TPP tilt plus CT or a 20 log10 LF level offset) and fpa_e. Below 5 kt the sphere nose follows the INS heading. The default speed reference is the ground, unless the database says otherwise.

### Cited Findings
- Per state: LF = |n|, where n is acceleration in g with gravity added; CT = LF·W/(rho·V_tip²·pi·R²); mu = |u|/V_tip, u being the air velocity, which equals ground velocity unless a wind is set — [N/src/trajectory.cpp:458-472](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp).
- fpa = asin(u_z/|u|); fpa_e = asin(û·n̂), the path angle relative to the TPP whose normal is the load-factor direction. Both are 0 at true hover — [N/src/trajectory.cpp:429-452](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp). `-q/--qsam` (default on) selects fpa_e for the lookup — [N/src/niceops.cpp:2322, 775](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp). Reviewer: `--no-qsam` moves the demo by 0.10 dB median, 1.7 dB max — [R/code review sec. 10 D](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Normalisation and clamping: each axis is normalised by its database min/max and clamped to [0, 1] (`clamp01`) — [N/src/sphere_database.cpp:1206-1250, 1160-1175](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp).
- Triangulation: `ConditionInterpolator` builds a Delaunay mesh in 2 or 3 dimensions (Bowyer-Watson on Shewchuk exact predicates, `src/delaunay.cpp`, `src/predicates.cpp`). Every condition on a face of the unit cube gets a copy pushed to -0.1 or 1.1 (`extend_domain_by = 0.1`). A query that lands exactly on a condition takes that sphere whole; outside the hull it takes the nearest condition ("rather than extrapolating past the envelope"). Fewer than dim+1 points means nearest only — [N/src/interpolator.cpp:31-152](/Users/evg5332/Code/D/NICEOPS/src/interpolator.cpp); [N/src/interpolator.hpp:46-52](/Users/evg5332/Code/D/NICEOPS/src/interpolator.hpp).
- Load factor, two policies:
  - (a) `fixed_load_factor` databases: only the LF = 1 sphere is stored, the triangulation is 2-D over (mu, gamma), and a level correction of +20 log10(CT/CT_ref) is added only above the reference, with no upper cap.
  - (b) other databases: panam materialises copies at LF = [0, 1, 1.1, 1.2, 1.5, 2, 3, 5], each +20 log10(LF), with LF = 0 carrying the 1 g level, and CT is a real third axis that clamps at 5 g.
  - Reviewer: LF 10 and 20 give +13.98 dB on M3 (clamped) against +20.0 and +26.0 on fixed_load_factor (unbounded) — [N/src/sphere_database.cpp:1220-1245](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [P/flight_acoustics.py:3452-3500](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [H/spheres/rebuild.py:66-76](/Users/evg5332/Code/Python/niceops-validation/spheres/rebuild.py); [R/code review X6](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Synthetic envelope: panam re-emits every near-level sphere at gamma = -24 and +35 deg ("on the assumption that directivity at a given airspeed carries over to a steeper flight path"), and writes the hover at gamma {-24, 0, 35}. Shipped databases span gamma -24 to 35 — [P/flight_acoustics.py:3184-3190, 3328-3350](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); netCDF read (B407 gamma -24.0 to 35.0, 115 unique). Memory: "The gamma -24/+35 database entries are copies of the level sphere at each mu" — [M/propagation-roadmap.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/propagation-roadmap.md).
- Hover (mu = 0) is synthesised from the slowest near-level sphere, averaged fore to aft in energy, with SPLA and EAA recomputed from the averaged spectrum. An optional `hover_correction` callable exists but is not shipped — [P/flight_acoustics.py:3299-3326](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [R/Grazing incidence findings 2026-10-02.md, "The hover"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
- Heading at hover: below `heading_nose_speed` = 5 kt the nose is the AC track's column-8 heading, turned into the local frame by an offset fitted from states over 20 kt (circular mean) or from `--frame_bearing`. With neither, the nose follows the drift, with a warning. The switch is a hard threshold — [N/src/trajectory.cpp:799-835](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp); [N/src/trajectory.hpp:80](/Users/evg5332/Code/D/NICEOPS/src/trajectory.hpp); [N/src/niceops.cpp:763-797](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp). Reviewer: a 30-90 deg nose jump at the threshold gives 2.2-4.4 dB steps (X4).
- Speed reference: `--speed_reference ground|air`, defaulting to the database's root attribute `speed_reference` (absent means ground). Air subtracts the `--atmosphere` wind at the state's height from mu, fpa, fpa_e and the nose, while LF and CT stay inertial. panam does not write the attribute yet, so every database is ground-referenced — [N/src/trajectory.cpp:761-777](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp); [N/src/sphere_database.cpp:103-140](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [N/README.md:372-393](/Users/evg5332/Code/D/NICEOPS/README.md). panam labels spheres with ground speed (`advance_ratio = 0.514444·speed_knots/V_tip`), and "the tracking data carries no airspeed" — [P/flight_acoustics.py:3466-3468](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [P/noise_abatement_2017.py:793-796](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py).
- Vertical-speed sign: the AC format is inconsistent (the 2017 tracks store vz positive down). NICEOPS detects this by correlating vz with dz/dt and negates it — [N/src/trajectory.cpp:474-524](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp).

### Inferences
- There is no explicit acceleration or transient term (for example, BVI onset during a deceleration or manoeuvre). Acceleration acts only through the TPP tilt, CT or the 20 log10 LF offset, and fpa_e. That is the QSAM approach.
- Steep-descent directivity beyond the measured gamma range is level-flight directivity relabelled. A steep approach blends toward copies of level spheres, not toward measured descent spheres.
- Tip-speed variation (rotor RPM) is not a parameter. The database's tip speed is a constant, and a JSON mismatch is only logged (code review X2).

### Gaps
- I did not verify how panam's per-run `flight_condition` defines gamma (ground track, presumably) beyond the docstring.

## 3. Sphere orientation and within-sphere interpolation

### Takeaway
The sphere frame is the TPP normal (the unit load-factor vector, so in steady flight it is vertical) with the nose taken as the velocity projected into the TPP and starboard as nose x normal. The angles come from three dot products, and the level from periodic bilinear interpolation in dB on the (phi, theta) grid.

### Cited Findings
- `sphere_frame(velocity, normal)`: normal = n̂ (load factor); nose = v̂ projected onto the TPP and normalised; starboard = nose x normal — [N/src/sphere_database.cpp:678-689](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp). The normal passed in is the state's (nx, ny, nz) load-factor unit vector — [N/src/niceops.cpp:773-776](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- theta = acos(nose·d); phi = atan(starboard·d / |normal·d|), with the 0/0 pole case caught — [N/src/sphere_database.cpp:711-734, 2062-2079](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp).
- Grid lookup: `GridAxis::locate` (uniform-axis arithmetic, periodic wrap in phi) plus `cell_stencil`, giving four corner weights shared by the broadband and spectral paths — [N/src/sphere_database.cpp:544-640, 454-476](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [N/CLAUDE.md "Why the hot path is shaped..."](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md).
- The measurement frame matches the playback frame only in steady flight. panam's `hemigen` files samples by elevation relative to the horizon and azimuth relative to the ground-track heading — [P/flight_acoustics.py:867-924](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py). Memory notes "both sides use a level frame in steady flight" — [M/propagation-roadmap.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/propagation-roadmap.md).
- Known edge: `sphere_frame` returns a zero nose for exactly zero velocity, or for velocity parallel to the TPP normal (with no heading available) — [R/code review sec. 3](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).

### Inferences
- In manoeuvres the playback sphere tilts with the TPP while the spheres were built in a level frame. This is the standard QSAM "rotate the sphere with the TPP" assumption.

### Gaps
- No interpolation in log-frequency or on spectral shape: each band is interpolated independently in dB.

## 4. Doppler, emission and reception time, convective amplification

### Takeaway
NICEOPS applies no frequency shift and no explicit convective-amplification factor. The spheres are built from moving-source measurements with `DOPPLER_SHIFT_REMOVED = 0`, so the received shift and the amplitude effects for each sphere's own speed are baked in. The time history is sampled in emission time and integrated in reception time, weighted by (1 + (dr/dt)/c).

### Cited Findings
- `reception_factor` = max(1 + (dr/dt)/c, 1e-3) per sample, with dr/dt from central differences of range (or of c·travel_time under refraction) — [N/src/util/sel.hpp:76-97](/Users/evg5332/Code/D/NICEOPS/src/util/sel.hpp); [N/src/niceops.cpp:1886, 1927](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp). Reviewer: omitting it would move SEL by 0.397 dB on the check set, and the demo by +0.078 dB median, +0.41 max — [R/code review sec. 10 D](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Sound speed for reception time comes from the vehicle-file atmosphere, or from the surface level of `--atmosphere` — [N/src/niceops.cpp:977-999](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- Doppler convention: `Be407102.nc` has `DOPPLER_SHIFT_REMOVED = 0`, and panam never calls `dedopplerize` in the build, so NICEOPS correctly applies no shift. If a database ever sets it to 1, the shift would have to be applied before absorption, ground and A-weighting. Blending spheres of different speeds blends different shifts — [M/sphere-doppler-convention.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/sphere-doppler-convention.md).
- panam's `hemigen` computes `mach_r`, but `depropagate_hemisphere` discards it (`az_deg, el_deg, r_geom, t_obs, _ = hemigen(...)`). Depropagation applies no Doppler-amplitude or convective correction, only ambient gating, receiver response, spreading and absorption — [P/flight_acoustics.py:1252, 1358-1395](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py).

### Inferences
- Convective amplification is therefore implicit and tied to each sphere's measured speed. It is correct at the speeds measured and blended linearly in dB between them. An air-referenced database would also need the shift and amplification defined against airspeed.

### Gaps
- Whether the source-time/observer-time mapping in panam uses straight-line retarded time only (it does in `hemigen`, t_obs = t + r/c), or uses ray travel time when `ray_model` is set, was not traced in code (the docstring says the ray model supplies reception time).

## 5. Propagation physics

### Takeaway
The default chain is spherical spreading from the 30.48 m sphere, plus ISO 9613-1 absorption per band (spectral, the default) or the database's broadband EAA rate (`--no-spectral`). No ground effect is applied by default. Opt-in layers: stratified absorption; Chien-Soroka ground effect with four impedance models, roughness and turbulence coherence; a flush ground-plate (GR1425 BEM) or point-receiver model; curved rays (linear or layered eigenrays) with shadow zones, creeping-wave decay, a shadow-boundary diffraction table, duct multi-bounce arrivals, Airy tunnelling and skip-zone fills; turbulent scattering into shadows; and external PE propagation tables. Ground is flat at z = 0 everywhere. There is no terrain, shielding or urban model.

### Cited Findings
- **Spreading**: 20 log10(r/R_sphere), extrapolated inward without a floor at HEAD. Inside the sphere the level grows and r = 0 gives +inf; 28 demo observers are affected by up to 4.4 dB. The review branch floors it at 0.01 R — [N/src/sphere_database.cpp:1343-1344, 1415](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [R/code review sec. 4, sec. 11](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- **Broadband EAA**: EAA·(r - R)/1000, EAA being SPLA - SPLA_after_1000 m of ISO 9613-1. panam computes it at 293.15 K, 101.325 kPa, 20% RH by default. The D original's 4x multiplier was removed (the footprint rose 9.93 dB mean) — [N/src/sphere_database.cpp:1322-1344](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [P/flight_acoustics.py:3054-3070, 3825-3860](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [N/README.md:261-277](/Users/evg5332/Code/D/NICEOPS/README.md). Reviewer: EAA is linear in range, anchored at 1 km, while the true A-weighted loss is concave. With matched atmosphere, spectral minus broadband runs from -0.75 dB (0-250 m) to +1.76 dB (1.5-1.75 km) (X5).
- **ISO 9613-1 absorption** (spectral path), ported from panam's `iso_9613_1_1993.py` and pinned in tests. The atmosphere comes from the vehicle JSON's `Atmosphere` section or `--temperature/--pressure/--humidity`. Humidity defaults to 20% and pressure is derived from density·R·T. Reviewer: these silent defaults are worth 1.2-3.0 dB median (X3) — [N/src/atmosphere.cpp](/Users/evg5332/Code/D/NICEOPS/src/atmosphere.cpp); [N/src/trajectory.cpp:389-424](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp); [N/README.md:293-317](/Users/evg5332/Code/D/NICEOPS/README.md).
- **Layered atmosphere** `--atmosphere` (CSV of T, RH, optional p and wind against height): per band, the mean absorption coefficient over the heights the path spans, times the path length. The profile must reach the aircraft (no extrapolation upward), and the bottom level is held to the ground. Spectral only. Wind is used only by refraction and `--speed_reference air` — [N/src/atmosphere_profile.cpp](/Users/evg5332/Code/D/NICEOPS/src/atmosphere_profile.cpp); [N/README.md:319-370](/Users/evg5332/Code/D/NICEOPS/README.md).
- **Ground effect** `--ground_effect` (off by default; needs the spectral path): Chien-Soroka spherical-wave reflection with Chessell boundary-loss correction and band averaging, ported from panam `ega()` and pinned to 1e-9 dB. It is a gain (+6 dB over rigid ground). Options:
  - impedance `delany_bazley` (default), `miki`, `variable_porosity` (sigma_e, alpha_e) or `delany_bazley_layer`;
  - `--flow_resistance`, default 225 kPa·s/m²;
  - turbulence `chessell` (rad·s/√m) or `harmonoise` (`--gamma_t`);
  - `--roughness`, applied as exp(-2(k·sigma·cos theta)²).
  - Pure-tone mode is deliberately not ported.
  - [N/src/ground_effect.hpp:1-90, 241-248](/Users/evg5332/Code/D/NICEOPS/src/ground_effect.hpp); [N/README.md:176-185, 476-506](/Users/evg5332/Code/D/NICEOPS/README.md). The site fit for Amedee is variable porosity sigma_e about 200 — [P/noise_abatement_2017.py:81](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py).
- **Receivers**: a flush microphone (height 0) with `--plate_table` uses panam's axisymmetric BEM of a GRAS GR1425 plate (400 mm, 8 mm, microphone at 3/4 R per SAE ARP 4055), p/p_free = P_d + Q·P_r. Shipped tables: `plates/gr1425_db225.nc` and `gr1425_vp200.nc`. Any other height, with `--ground_effect`, is a point receiver; the default height is 1.2 m and an observer file's fifth column sets it per observer. Per-trajectory ground tables are interpolated, or evaluated exactly with `--exact_ground` — [N/README.md:963-994, 1126-1201](/Users/evg5332/Code/D/NICEOPS/README.md); [N/src/niceops.cpp:2255](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- **Ground geometry is flat at z = 0**: `straight_ray` gives the ground model source_height = source.z and the horizontal ground distance. An observer's own z shifts only the direct path — [N/src/propagation.hpp:68-81](/Users/evg5332/Code/D/NICEOPS/src/propagation.hpp). `--ground_elevation` is one constant datum subtracted from geodetic altitude — [N/README.md:148, 1250](/Users/evg5332/Code/D/NICEOPS/README.md). Propagation tables also "assume one stratified atmosphere over flat, uniform ground"; ground elevation is listed as not done — [N/README.md:1010-1013](/Users/evg5332/Code/D/NICEOPS/README.md); [N/CLAUDE.md propagation_table "Not done"](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md).
- **Refraction** `--refraction none|linear|stratified` (default none; needs `--atmosphere`):
  - *linear*: a per-path least-squares line through the effective sound speed c(z) + w(z)·h (Harmonoise/Nord2000 equivalent-linear), with closed-form circular-arc direct and reflected eigenrays. The ray gives the sphere launch direction, the spreading r/cos(launch), absorption along the arc, the ground geometry and the travel time.
  - *stratified*: piecewise-linear layers and eigenrays over four p-families, every arrival summed in energy, duct multi-bounce arrivals with plane-wave |R|² per bounce, 32 launch samples, and focusing capped at +10 dB per ray. Hand-over to the linear fit occurs past 0.9 of the layered boundary.
  - [N/src/refraction.cpp](/Users/evg5332/Code/D/NICEOPS/src/refraction.cpp); [N/src/stratified.cpp:790](/Users/evg5332/Code/D/NICEOPS/src/stratified.cpp); [N/src/stratified.hpp:46, 182](/Users/evg5332/Code/D/NICEOPS/src/stratified.hpp); [N/README.md:395-475, 508-705](/Users/evg5332/Code/D/NICEOPS/README.md).
  - Limitations (README): the effective-sound-speed approximation overstates wind by 1/cos(elevation), and under linear "the source does not feel the wind".
- **Shadow zones**: the grazing ray's level, less the least-attenuated creeping-wave mode (over rigid ground or the ground's own impedance), capped by `--shadow_floor` (default 20 dB). Near the boundary (point and plate receivers) a residue-series table V(xi, s), built at a reference source height of 20 l, blends from 2F before the boundary. The per-band boundary law has fitted constants: 0.363 at 1.2 m, 0.75 flush, and fallback multiples 1.60/1.90·(f/1 kHz)^-0.08 — [N/README.md:413-437, 581-705](/Users/evg5332/Code/D/NICEOPS/README.md); [N/CLAUDE.md ShadowTransition](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md).
- **Micro-shadow and skip-zone options** (stratified only, off by default): `--tunnelling` (Airy tail below the turning point) and `--skip_fill` (caustic Airy tail into the skip zone). "None moves the SEL" (under 0.02 dB) — [N/README.md:1038-1091](/Users/evg5332/Code/D/NICEOPS/README.md).
- **Turbulence**: in the ground model, as coherence (Chessell or HARMONOISE gamma_T). Into shadows, `--scattering gamma_T` adds an empirical S = a + 10 lg gamma_T + b lg f + c lg R + d lg z_s dB re free field, fitted to the turbulent PE. It is off by default because the measurements want at most about 3e-6 — [N/README.md:869-905](/Users/evg5332/Code/D/NICEOPS/README.md).
- **PE**: a wide-angle split-step Padé PE (Collins), with self-starter, fourth-order grid and impedance ground, used as a reference (`niceops_pe_compare`, `niceops_benchmark`) and to build `--propagation_table`. A table gives each band's level re free field by (bearing, source height, receiver height, range, band) for one stratified atmosphere over flat uniform ground; uncovered pairs fall back to the run's own model. Benchmark agreement: within 0.3 dB of the exact solutions for Attenborough et al. 1995 cases 1-3 — [N/README.md:706-797, 907-961, 996-1124](/Users/evg5332/Code/D/NICEOPS/README.md); [N/src/propagation_table.hpp](/Users/evg5332/Code/D/NICEOPS/src/propagation_table.hpp).
- **Terrain, shielding and urban: not implemented.** A grep of `src/` for urban/terrain/building/shield/barrier finds no implementation. The roadmap's step 7 ("per-receiver urban reciprocity tables") is planned only — [N/NICEOPS Atmospheric and Urban Propagation Options.md, "Steps"](/Users/evg5332/Code/D/NICEOPS/NICEOPS%20Atmospheric%20and%20Urban%20Propagation%20Options.md); [M/propagation-roadmap.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/propagation-roadmap.md).
- Ambient: `-a` mixes an ambient level into every sample *coherently* (20 log10 of summed pressure amplitudes, +6.02 dB at equal levels). This is faithful to the D original and was left for a maintainer decision — [R/code review sec. 4](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).

### Inferences
- Against AAM/NORAH2-class tools, NICEOPS's distinctive capabilities are:
  - the spectral path with per-band ISO absorption plus a finite-impedance spherical-wave ground model with a microphone-plate BEM;
  - curved-ray refraction with diffraction-corrected shadow edges and duct arrivals;
  - an interface for PE propagation tables.
- Its missing pieces are terrain (even sloped flat ground), shielding, lateral attenuation or engine-installation effects, and any urban or building treatment.

### Gaps
- Not verified: whether `--refraction` with the broadband (non-spectral) path is refused or supported. README says `--atmosphere` is "spectral path only".

## 6. Receivers and metrics

### Takeaway
Observers are a regular grid (default 120 x 120, ±1500 m at 25 m, end-exclusive) or a scattered CSV with optional per-microphone heights. The only per-observer metric is A-weighted SEL. Integration covers the window from the first to the last sample within 10 dB of the peak, in reception time. `-D` is a night-penalised (+10 dB, 22:00-07:00) energy sum over trajectories, **not normalised by 86 400 s**, so it is not a standard DNL. There is no LAmax map output, and no PNL/EPNL/tone correction, C-weighting or loudness. Time histories, band histories and a summary metrics file are available.

### Cited Findings
- Grid: default 120x120 spanning -1500 to 1475 m, D-style end-exclusive `iota`; `--grid`, `--gridstep` and `--units ft|m` override it — [N/src/niceops.cpp:56-66, 486-503](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp); [N/README.md:68, 164](/Users/evg5332/Code/D/NICEOPS/README.md).
- Scattered observers: `x,y,z`, `id,x,y,z` or `id,x,y,z,height` (feet by default) — [N/README.md:68, 974-988](/Users/evg5332/Code/D/NICEOPS/README.md).
- SEL window: `ten_db_down_interval` runs from the first to the last sample at or above peak - 10 dB, dips included ("the EPNL duration of 14 CFR 36 A36.4.5"). SEL = 10 log10(T·Σ 10^(L/10)·w_i) — [N/src/util/sel.hpp:22-61, 140-166](/Users/evg5332/Code/D/NICEOPS/src/util/sel.hpp). Doc/code conflict: CLAUDE.md still describes the D-style first-crossing walk from the peak, and the code review flags it stale (the change is commit `9d6e986`) — [N/CLAUDE.md "The SEL integral is oversampled..."](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md); [R/code review sec. 1](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Integration step: the trajectory's own mean sample interval (0.02 s on shipped tracks), not `-s` — [N/src/niceops.cpp:729-751](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- DNL: `day_night_level` adds 10 dB to night samples (22:00-07:00 local, via `--utc_offset`) before the window is found. Energies are summed over all `-t` and written as 10 log10(Σ). There is no 86 400 s normalisation: "With one daytime operation the DNL equals the SEL exactly" — [N/src/util/sel.hpp:127-186](/Users/evg5332/Code/D/NICEOPS/src/util/sel.hpp); [N/src/niceops.cpp:1970-1978, 2684-2687, 2874-2881](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp); [N/README.md:209-221](/Users/evg5332/Code/D/NICEOPS/README.md). Reviewer: this is 49.4 dB above an ANSI/FAA DNL for the same event (X5).
- Other outputs:
  - `-i`: per-observer level time histories;
  - `--band_history`: per-state band levels with emission/reception time, range, elevation, phi, theta and level_dBA;
  - `-m`: average and maximum over observers plus the fraction of observers at or above 65 dB;
  - `-S`: trajectory state;
  - `--powerfile`: shaft power, if the database has it (none does);
  - `--latlon`; `--export_sphere` / `--export_aam`.
  - [N/README.md:72-91](/Users/evg5332/Code/D/NICEOPS/README.md); [N/src/niceops.cpp:237-257](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- Weighting: A-weighting only (IEC 61672 pole-zero form). `--band_limits` (default 20-20000 Hz) restricts the bands summed — [N/src/weighting.hpp:23-63](/Users/evg5332/Code/D/NICEOPS/src/weighting.hpp); grep for PNL/NOY/tone/loudness in `src/` found nothing.
- `--db_cutoff` culls provably sub-threshold evaluations. The review found its README "exact" claim wrong per cell (157 of 592 cells at or above 85 dB change, by up to 0.43 dB), while the per-peak-level guarantee holds — [R/code review sec. 3](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).

### Inferences
- A comparison with AAM or NORAH2 should note that NICEOPS reports SEL (LAE) only. LAmax and time-above could be post-processed from `-i` histories but are not built in. Its "DNL" is a penalised cumulative SEL.

### Gaps
- No built-in Leq/time-above metrics found.

## 7. Trajectory input

### Takeaway
Two formats are accepted. The 27-column NASA/FAA AC state CSV carries positions, velocities and accelerations; a header-named WGS84 track (time, lat, lon, alt, optional vE/vN/vUp) is projected locally and smoothed by a spline with an RMS tolerance (default 10 ft). Velocity and load factor come from spline derivatives. Resampling or smoothing of AC tracks (`adapt()`) is disabled. `--trim` crops to the lowest or slowest point. Several trajectories per run are supported, each with its own database and vehicle file.

### Cited Findings
- Format is decided by width: 27 columns is AC (feet); anything else is read by header names — [N/README.md:62-70, 1229-1260](/Users/evg5332/Code/D/NICEOPS/README.md); [N/src/trajectory.cpp:526-560](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp).
- Geodetic path: a smoothing spline per axis with `--path_smoothing` RMS tolerance (default 10 ft). Velocity is the first derivative; load factor is the second derivative in g, with +1 on z — [N/src/trajectory.cpp:290-347](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp). Error table: 0.085 dB median at the default, 11.6 dB at exact interpolation — [N/README.md:1275-1303](/Users/evg5332/Code/D/NICEOPS/README.md).
- Doc/code conflict: the projection is described as transverse Mercator in code and README ([N/src/geodetic.hpp:8](/Users/evg5332/Code/D/NICEOPS/src/geodetic.hpp)), while the code review says it is the American polyconic and fixed the docs on its branch — [M/code-review-2026-10-02.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/code-review-2026-10-02.md). Not adjudicated here.
- Reviewer: the natural end conditions of the smoothing spline force zero acceleration at the track ends, worth up to 1.8 dB on a geodetic demo footprint and 7.5-11 deg of TPP-normal error over the last seconds — [R/code review sec. 4](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- `adapt()` (uniform resampling with smoothing) is commented out in `process_trajectory` — [N/src/niceops.cpp:710-719](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- `--trim X` (`watts_trim`) crops from x = X to whichever comes first, the minimum z or the minimum mu. It is off by default — [N/src/trajectory.cpp:779-797](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp).
- Multiple trajectories: `-t` is repeatable, and `-d`/`-v`/`-o`/`-S` take one value or one per track. A bad track is logged and skipped. Plotting is skipped. At HEAD, `-i`/`--band_history` overwrite each other across tracks (fixed on the review branch) — [N/README.md:209-226](/Users/evg5332/Code/D/NICEOPS/README.md); [R/code review item 5](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Shipped tracks: AW139 41 files, B205 13, S-76D 30 (`N/tracks/`). The AC tracks' local frame runs +x along bearing 210 deg (shipped tracks) or 270 deg (2017 test) — [N/README.md:1322-1329](/Users/evg5332/Code/D/NICEOPS/README.md); [H/spheres/rebuild.py FRAME_BEARING](/Users/evg5332/Code/Python/niceops-validation/spheres/rebuild.py).

### Inferences
- NICEOPS consumes trajectories with their kinematics; it does not generate procedures. There is no flight-dynamics or trim model beyond deriving LF and CT from accelerations, and the vehicle JSON's drag is unread (reviewer).

### Gaps
- None material.

## 8. Multiple sources and source position

### Takeaway
There is one source per aircraft: a single total-vehicle sphere (main rotor, tail rotor, engine and everything else together, as measured), placed at the trajectory point. There is no separate tail-rotor or engine source, no offset of the source from the tracked point, and no rotor-shaft tilt separate from the TPP-normal assumption.

### Cited Findings
- `source_state` uses `Vec3(state.x, state.y, state.z)`, the track position, as the sphere centre — [N/src/niceops.cpp:773-776](/Users/evg5332/Code/D/NICEOPS/src/niceops.cpp).
- The vehicle file carries only `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` and `Atmosphere` — [N/src/trajectory.cpp:385-424](/Users/evg5332/Code/D/NICEOPS/src/trajectory.cpp). A grep of `src/` for tail rotor or source offset finds nothing.
- Spheres are depropagated from ground-board measurements of the whole aircraft — [P/noise_abatement_2017.py:754-990](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py). The hover decomposition shows that main-rotor and tail-rotor tones are not separable below about 300 Hz in the measured data — [M/grazing-pairs-refraction.md "Hover decomposition"](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md).

### Inferences
- Source-separation effects (tail rotor near the ground, different directivity frames per source) are folded into the single empirical sphere.

### Gaps
- None.

## 9. Database generation in panam

### Takeaway
panam builds each source sphere from steady-flight ground-board recordings in four steps:
1. Emission-time geometry along straight lines (or, newly, curved rays).
2. 0.5 s third-octave frames, ambient-gated per bin against a measured ambient run (10 dB SNR) with ambient subtracted, and the plate's BEM response divided out (capped).
3. Spreading undone to 100 ft and ISO 9613-1 absorption undone (capped at 30 dB) at the run's station atmosphere.
4. Samples gridded by Shepard IDW (now adaptive, on a 2 deg grid), restricted to emission angles at least 2 deg below the horizon and ranges up to 2000 ft, then written to the AAM 10 deg x 5 deg lower-hemisphere grid.

`build_empirical_database` then mirrors the upper hemisphere, synthesises the hover, adds synthetic gamma copies and load-factor copies, and computes dBA and EAA (at 1000 m, 20 deg C, 20% RH). There is no Doppler removal and no convective correction.

### Cited Findings
- `depropagate_hemisphere` defaults: `r_ref = 100` ft, `freq_range = (0, 2000)`, `window_time = 0.5`, `azi_step = elv_step = 10`, `rmax = 25`, `apply_absorption_deprop = False`, `min_elevation_deg = 0`, `band_snr_gate_db = 3`; ambient gating only when an ambient reference is passed — [P/flight_acoustics.py:927-967](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py).
- The 2017 `build_sphere` defaults override these:
  - `band_snr_gate_db = 10`; gate against the measured ambient run (`gate_ambient = True`); `max_absorption_correction_db = 30`; `apply_absorption_deprop = True`;
  - `min_elevation_deg = 10`; `max_propagation_range_ft = 2000`;
  - `board_correction = 'plate_bem'`: divide by the plate's BEM response on SITE_GROUND. A flat 0.5 pressure (-6 dB) was used before 2026-09-26;
  - atmosphere from the run's ground station, fallback 293.15 K / 101.325 kPa / 20%;
  - `third_octave = True` on the legacy band centres and grid.
  - [P/noise_abatement_2017.py:754-947](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py).
- Harness rebuild settings for the shipped databases: `COMMON = dict(band_snr_gate_db=10, max_absorption_correction_db=30, max_propagation_range_ft=2000, point_stride=10, steady_only=True, board_correction='plate_bem')`, `LOAD_FACTORS = [0, 1, 1.1, 1.2, 1.5, 2, 3, 5]`, `EXTENDED_FPA = (-24, 35)`, frame bearing 270 — [H/spheres/rebuild.py:58-100](/Users/evg5332/Code/Python/niceops-validation/spheres/rebuild.py).
- The shipped B407 recipe (NICEOPS `5d761b8`) is:
  - a 2 deg rim (`min_elevation_deg = 2`);
  - curved-ray (stratified) depropagation for the 51-52 runs that have profiles, straight rays otherwise;
  - the board correction capped (`max_response_correction_db = 15`);
  - an adaptive Shepard grid on 2 deg (aspect 5, `shepard_floor` 1.0).
  - [M/grazing-pairs-refraction.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md); [R/Grazing incidence findings 2026-10-02.md, "Databases against the day-286 approaches"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
- Depropagation core order: ambient gate and subtraction, then the receiver response, then spreading ×(r/r_ref)², then absorption ×10^(alpha(r - r_ref)/10) with the cap — [P/flight_acoustics.py:1358-1395](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py). Gridding averages in linear power (`weights @ power_pts`) — [P/flight_acoustics.py:~1615-1640](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py).
- Adaptive interpolation defaults: k = 8, min_mics = 3, kappa = 1.3, resolution_factor = 0.5, max_radius 60 deg. The radius is R = max(1.3·d_8, 1.3·d_3mics, 0.5·s). Nodes beyond the cap are NaN gaps, later filled by the AAM writer's resampling — [P/flight_acoustics.py:4697-4710](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [M/grazing-pairs-refraction.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md).
- Database assembly: `extract_SPL` gives SPLA and EAA (EAA = SPLA - SPLA after `distance` = 1000 m of absorption at 293.15 K, 101.325 kPa, 20% RH). Optional `infreqs` band limit (the harness rebuild sets none, so B407 carries 10 Hz-10 kHz). Then `_complete_sphere` mirror, hover from the slowest near-level sphere, extended gammas, load-factor copies, and `add_sphere_group` (mu from ground speed, CT = LF·weight coefficient, dBA + 20 log10 LF) — [P/flight_acoustics.py:3054-3500, 3825-3860](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py).
- Doppler: not removed (`DOPPLER_SHIFT_REMOVED = 0`) — [M/sphere-doppler-convention.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/sphere-doppler-convention.md). panam can also write NORAH2 `.hem` files from the same hemisphere (`write_norah2_hemisphere`) — [P/flight_acoustics.py:2462](/Users/evg5332/Code/Python/niceops-validation/external/panam/flight_acoustics.py); [P/noise_abatement_2017.py:961-975](/Users/evg5332/Code/Python/niceops-validation/external/panam/noise_abatement_2017.py).
- Contamination history: before ambient gating, top bands were amplified ambient (Be407 aft pole 216 dB at 10 kHz; onset at 3150 Hz aft). `Be407_M2` and `M3` differ only by a 2000 Hz versus about 3500 Hz cut — [M/nod-database-generator.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/nod-database-generator.md). The current build's peak check refuses spheres with any band above 194 dB — [H/spheres/rebuild.py MAX_PEAK_DB](/Users/evg5332/Code/Python/niceops-validation/spheres/rebuild.py). The shipped B407_fixed_load_factor maximum dBA is 114.5 (netCDF read).
- Known builder caveats:
  - The gate *zeroes* bins below SNR, so the rim can only read low — [M/grazing-pairs-refraction.md "Methods research round 3"](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md).
  - Fore/aft averaging in the synthetic hover erases a measured 6-8 dB fore-aft asymmetry at 315-1250 Hz. A measured hover correction exists but is not shipped — [R/Grazing incidence findings, "The hover"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
  - Validation contamination: 12 of the B206L3 approach cards were built into the shipped Be206L3 database (about 0.07-0.1 dB rms advantage) — [M/grazing-pairs-refraction.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md).

### Inferences
- The panam spheres are "measured-at-the-ground, depropagated" source models, like AAM/RNM empirical spheres. They differ in the explicit ambient gate, the plate-BEM board correction, the optional curved-ray depropagation and the adaptive gridding. The sphere radius (100 ft) and the AAM-compatible grid make them exportable to AAM and NORAH2 formats.

### Gaps
- Not checked: whether the curved-ray, cap and adaptive changes are on panam `master`. The report says they live only on the harness's public-lineage clone and still need upstreaming.

## 10. Validation status against the 2017 Noise Abatement (Amedee) data

### Takeaway
B407 approach SEL (day 286, 25 approaches, 1,112 mic-runs, 10 dB-down window like for like) is predicted at **1.26 dB rms with the station atmosphere and 1.37 dB rms with stratified refraction** using the shipped database. The 5-10 deg peak-elevation bin is the worst (1.89 / 2.03 dB rms). On the windy day 284 (held-out spheres), stratified refraction plus plates brings overall rms to about 1.7 dB (straight rays 2.36 after removing cut records). The open residuals are source-side near the horizon: the low, slow approach end runs +3 to +4 dB, and fast, high flight at long range forward runs -1.5 to -2.8 dB at 2-6 deg.

### Cited Findings
- Shipped-database SEL scores against the day-286 approaches (1,112 mic-runs), station / stratified:

  | Database | SEL rms | 5-10 deg | Frames under 100 ft and 20 kt |
  | --- | --- | --- | --- |
  | `56ce39b` | 1.31 / 1.42 | 2.10 / 2.28 | +4.45 |
  | `041f39b` | 1.28 / 1.39 | | +3.90 |
  | `5d761b8` (shipped) | 1.26 / 1.37 | 1.89 / 2.03 | +3.89 |

  — [R/Grazing incidence findings 2026-10-02.md "Databases against the day-286 approaches"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
- Earlier milestones (harness README):
  - Rebuild with 10 deg -> 5 deg rim: the 5-10 deg bin goes +2.42/3.23 -> +1.21/1.89 and overall rms 1.39 -> 1.25. A 2 deg rim adds 0.07 dB.
  - Like-for-like 10 dB-down measured SEL: the 2 deg database is +1.43/2.05 at 5-10 deg and +0.27/1.30 overall.
  - [H/README.md "Findings so far"](/Users/evg5332/Code/Python/niceops-validation/README.md).
- Windy day 284 (16 steady runs, 738 mic-runs, held-out database):
  - 5-10 deg: straight +4.5/5.8, linear +3.7/4.5, stratified +0.3/3.6;
  - overall rms 2.36 / 2.02 / 1.85 after removing nine cut records; later 1.85 -> 1.70 after the shadow fixes;
  - upwind 5-10 deg: +0.28/3.06.
  - [H/README.md](/Users/evg5332/Code/Python/niceops-validation/README.md); [N/README.md:687-696](/Users/evg5332/Code/D/NICEOPS/README.md).
- Option-by-option overall rms (calm 286 / windy 284), from `options_study.py`:

  | Configuration | rms (dB) |
  | --- | --- |
  | broadband EAA + flat +6 dB | 3.40 / 4.34 |
  | spectral, default air | 2.32 / 3.55 |
  | soft-ground reflection | 1.38 / 2.43 |
  | plate BEM | 1.30 / 2.36 |
  | linear | 1.42 / 2.02 |
  | stratified | 1.45 / 1.70 |
  | + scattering 3e-6 | 1.46 / 1.69 |

  "The ground correction is the largest single step on both days" — [H/README.md "What each modelling option is worth"](/Users/evg5332/Code/Python/niceops-validation/README.md).
- Footprints, four runs (median / rms, stratified + plate):

  | Run | Day | Result (dB) |
  | --- | --- | --- |
  | 286270 | calm | +0.5 / 1.3 |
  | 286277 | calm | 0.0 / 1.9 |
  | 284205 | windy | +0.7 / 1.7 |
  | 284211 | windy | -0.8 / 1.3 |

  The flat correction over-predicts sidelines by 3-8 dB — [H/README.md footprints table](/Users/evg5332/Code/Python/niceops-validation/README.md).
- Board/pole transfer pairs at 0-2 deg (4,861 frames): straight rays 7.6 dB rms, linear 0.8, stratified 1.2, with VP sigma_e = 200 unchanged — [R/Grazing incidence findings "Board/pole pairs"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
- PE agreement: stratified matches the PE within 0.2 dB (energy mean) for sources at or above 660 ft on all five day-286 profiles. Downwind stratified is 0.4-2.1 dB rms from the PE, and upwind the PE's shadow onset lies between linear and stratified — [M/propagation-roadmap.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/propagation-roadmap.md); [N/README.md:803-833](/Users/evg5332/Code/D/NICEOPS/README.md).
- Residuals (energy means, day 286):
  - low, slow approach end (under 100 ft, under 20 kt): +3.9 dB, against +3.8 for straight rays, so source-side;
  - fast (at least 80 kt), high (400-2,000 ft) states heard 2-10 deg below the horizon from far ahead: -1 to -2 dB, flat with range, the same under all ray models;
  - SEL by speed: 40-60 kt -0.3, 60-80 +0.7, 80-100 -0.4, 100-140 -1.1 dB.
  - [R/Grazing incidence findings "Source-side residuals"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md); [M/forward-line-error-source.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/forward-line-error-source.md).
- Hover: the synthesised hover against the day-283 HOGE hovers at 18-25 deg is 4-13 dB too quiet at 315-1250 Hz to the sides and aft, and 7-19 dB too loud at 63 Hz. A measured correction brings it within ±1.7 dB in-sample but is not shipped — [R/Grazing incidence findings "The hover"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md).
- Second aircraft B206L3: SEL rms 1.50/1.49 (station/stratified) with today's recipe, against 1.48/1.47 with September's, in a clean comparison excluding contaminated A cards — [M/grazing-pairs-refraction.md "Second aircraft"](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/grazing-pairs-refraction.md).
- An earlier six-aircraft approach validation exists (139 runs, plate receiver, measured SEL with ambient removed) in `H/approach/` — [M/propagation-roadmap.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/propagation-roadmap.md).

### Inferences
- Validation is concentrated on one aircraft (B407) at one site (flat lakebed, Amedee), largely in-sample for the steady runs and one calm approach day. Errors of about 1.3 dB rms SEL are the headline, with near-horizon (grazing) emission the weak spot.

### Gaps
- No per-aircraft summary of the six-aircraft approach validation was read (the tables are in `H/approach/`; not opened here).
- No validation of LAmax or of other metrics, since none is produced.

## 11. Known open issues and limitations (docs, memory, code review)

### Takeaway
The arithmetic core has been independently re-implemented and agrees to about 1e-13 dB. The open issues are physics scope, source-side residuals near the horizon, silent defaults and input validation, a few documented invariants that fail, and stale documentation. Many robustness fixes exist on an unmerged branch.

### Cited Findings
- Physics scope not implemented:
  - terrain (flat ground, constant elevation);
  - shielding and urban;
  - multiple sources;
  - same_grid = false databases;
  - adapt() resampling;
  - the panam half of the air speed reference;
  - shaft power data (no database carries it);
  - PE-table ground elevation and refracted-direction sphere lookup.
  - [N/CLAUDE.md](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md); [N/src/sphere_database.cpp:779-781](/Users/evg5332/Code/D/NICEOPS/src/sphere_database.cpp); [M/nod-database-generator.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/nod-database-generator.md).
- Physics limitations acknowledged in the docs:
  - linear refraction uses one gradient per path, can be up to 10 dB off the PE near the ground (van Maercke & Salomons 2009), and the source does not feel the wind;
  - the effective-sound-speed approximation overstates wind by 1/cos(elevation);
  - the stratified 32-sample arrival search is unconverged in ducts (4.3% of pairs over 1 dB against 2048 samples; up to 10.5 dB in single cells);
  - the per-ray focusing cap lets summed gain reach 17 dB;
  - an elevated duct with near-ground turning (run 286288) is 11-21 dB under the PE;
  - the shadow boundary law cannot get below about 6 dB rms pooled at the shadow edge;
  - the PE band average uses edge weights (0.16-0.47 dB rms; fix not applied);
  - `--scattering` strength is uncertain.
  - [N/README.md:453-474, 557-568, 684-686](/Users/evg5332/Code/D/NICEOPS/README.md); [R/code review sec. 10 E S1, S5; sec. 4](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Source and database side:
  - the rim below the minimum elevation is a flat fill (2 deg rim shipped);
  - the forward rim is 1.5-2.8 dB quiet for fast, high flight;
  - the low, slow end is +3 to +4 dB;
  - the synthetic hover is unrepresentative;
  - steep gamma uses level-sphere copies;
  - the database is built on ground speed;
  - spheres blend test days (a day-286 forward-rim offset of +0.5 to +1.2 dB);
  - the load-factor law is 20 log10(LF), unbounded on fixed_load_factor databases.
  - [R/Grazing incidence findings "What is ruled out, and what remains"](/Users/evg5332/Code/D/NICEOPS/reports/Grazing%20incidence%20findings%202026-10-02.md); [R/code review X6](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Code-review top defects reproduced at the review baseline (fixed on unmerged `review-fixes-2026-10-02` unless noted):
  - `golden_compare` accepts NaN/inf;
  - bad input (flat time column, NaN coordinate, unit-slipped atmosphere) gives NaN/-inf or 4-10 dB shifts with exit 0;
  - errors only in the log;
  - a worker exception calls `std::terminate` at -j >= 2;
  - multi-track `-i` overwrite;
  - a metric geodetic track loses its frame;
  - DNL on a clockless geodetic track is +10 dB;
  - `DirectionBounds` is exceeded by up to 0.7 dB within 1 deg of the nose axis;
  - asserts compiled out despite the docs (`-DNDEBUG`);
  - inward spreading unbounded (not fixed in output terms: floored at 0.01 R only);
  - coherent `-a` ambient (not fixed).
  - [M/code-review-2026-10-02.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/code-review-2026-10-02.md); [R/code review sec. 2, 11](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).
- Vehicle-file consistency: only the tip speed is compared with the database, and only in a log line. A 10% weight error moves the footprint 0.75 dB median silently (X2).
- Stale docs at HEAD:
  - CLAUDE.md and README describe the first-crossing SEL window and "70% evaluated, 1.2x"; the code uses first-to-last crossing and evaluates about 94%;
  - "spectral opt-in" (the code defaults to on);
  - "databases not tracked";
  - projection naming (disputed).
  - [R/code review sec. 1](/Users/evg5332/Code/D/NICEOPS/reports/NICEOPS%20code%20review%202026-10-02.md).

### Inferences
- For a point-by-point comparison, treat README and CLAUDE.md statements about the SEL window, the spectral default and evaluation fractions as superseded by the code paths cited above.

### Gaps
- Panam-side assumptions have not been code-reviewed (the review lists them as not done) — [M/code-review-2026-10-02.md](/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/code-review-2026-10-02.md).
- Real 2017 hover or approach-end tracks have not been run through the code-review checks.
