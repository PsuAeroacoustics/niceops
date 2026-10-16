# NICEOPS at near-grazing incidence over soft ground: what is modeled, and the measured residuals

Compiled 2026-09-30 from local code, logs and notes only. No computation was rerun; every number is quoted from the file cited.

Path abbreviations used in the citations:

- `N/` = `/Users/evg5332/Code/D/NICEOPS/` (NICEOPS working tree, branch `band-history`, HEAD 84462d7)
- `H/` = `/Users/evg5332/Code/Python/niceops-validation/` (validation harness; NICEOPS and panam pinned as submodules)
- `P/` = `H/external/panam/` (panam, pinned at public master becbf63)
- `M/` = `/Users/evg5332/.claude/projects/-Users-evg5332-Code-D-NICEOPS/memory/` (project memory notes)

Sign convention throughout: error = NICEOPS − measured, so positive is over-prediction. "median / rms" pairs are in dB. Elevation bins are the **peak elevation of the aircraft seen from the microphone** in the SEL tables, and the **emission elevation** in the frame tables.

---

## 1. Ground model as implemented (point receiver, `--ground_effect`)

### Takeaway
The point-receiver ground term is a Chien–Soroka/Weyl–van der Pol spherical-wave reflection over a **locally reacting** surface admittance, evaluated at **one center frequency per third-octave band** with Chessell's band-averaging factors, and fed an "equivalent straight geometry" when rays are curved. The validation runs use variable porosity σe = 200 kPa·s/m², αe = 0, with no turbulence or roughness on the flush path.

### Cited Findings

**Formulation**
- Ported from panam's `ega()`; sources named in the header are Chien & Soroka 1975, Delany & Bazley 1970, Chessell 1977, Daigle, Embleton & Piercy 1979, and Wyle WR 85-19. Header states "Assumes level terrain, no wind, and straight rays." — `N/src/ground_effect.hpp` lines 3–16
- Plane-wave coefficient: R_p = (cos θ − β)/(cos θ + β), where cos θ = (h_s + h_r)/R_image is the cosine from the surface normal (the variable is named `cos_grazing`) and β = 1/Z is the normalized admittance (stored in a variable named "impedance"). — `N/src/ground_effect.cpp` lines 29–40, 208; `N/src/ground_effect.hpp` lines 117–119
- Numerical distance: √w = sqrt(i κ / (1 + β cos θ)) · (cos θ + β), with κ = π f R_image / c (= k R₂/2). — `N/src/ground_effect.cpp` lines 47–81, 236
- Boundary-loss term: F(w) = 1 + i √π · √w · w_F(√w), where w_F is the Faddeeva function (panam's form is 1 + i√(πw) e^(−w)(1 − erf(−i√w))). Spherical coefficient Q = R_p + F(w)(1 − R_p). — `N/src/ground_effect.cpp` lines 86–105
- **Cut-off:** for |w| > 500 the boundary-loss term is set to zero, so Q = R_p (plane wave). panam uses the same cut. — `N/src/ground_effect.cpp` lines 20–22, 83–86
- |Q| exceeds 1 toward grazing (the surface wave), reaching about 2.02 just inside |w| = 500, i.e. +9.6 dB, not +6 dB. This is the bound the peak search uses. — `N/src/ground_effect.cpp` lines 336–344
- An asymptotic boundary loss for large |w| was tried and rejected: error plateaus at 1.6e-5 at |w| ≥ 50 and 1.9e-6 at |w| ≥ 100, which cannot hold the 1e-9 dB pinning to panam. — `N/CLAUDE.md` ("The asymptotic boundary loss was tried and rejected")
- Faddeeva: Weideman 32-term rational approximation, max relative error 3e-13 over Re z ∈ [0.09, 22.3], Im z ∈ [−1.2, 15.6]; the lower-half-plane reflection is safe only in that domain. Fast paths (asymptotic for |z| ≥ 6, Taylor table below) agree within 1e-11. — `N/CLAUDE.md` (`src/faddeeva` entry); `N/README.md` "Ground-plane microphones"
- The model is pinned to panam's `ega()` at 1e-9 dB and to `ground_plane.pole_level` (15 cases, 1e-6 dB). — `N/CLAUDE.md`; `P/docs/ground_plane_corrections.md` "In NICEOPS"

**Band averaging (one frequency per band)**
- The band level is Chessell's eqs. 20–21: magnitude = 1 + |Q|²/ρ² + (2/ρ)·[Re Q cos(6.325159 f τ) − Im Q sin(6.325159 f τ)] · T · sinc(0.727477 f τ), with ρ = R_image/R_direct, τ = (R_image − R_direct)/c, T the turbulence factor. Q, β and f are all taken **at the band center only**; the band average enters only through the sinc factor and the modified phase constant. — `N/src/ground_effect.cpp` lines 16–18, 250–306
- "Pure-tone mode is deliberately not ported." — `N/src/ground_effect.hpp` lines 230–234
- Contrast: the plate receiver uses five sub-frequencies per band (section 2), and the PE reference uses `--samples` frequencies per band. — `N/src/plate_receiver.hpp` lines 25–28; `N/README.md` "Turbulence, and the energy it scatters into shadows"

**Impedance models** (`--impedance_model`, `--flow_resistance`, `--porosity_rate`, `--layer_depth`)

| model | formula as coded | parameters |
| --- | --- | --- |
| `delany_bazley` (default) | Z = 1 + 9.08 X^−0.75 + i 11.9 X^−0.73, X = f/σ | σ, default 225 ("grassy field") |
| `miki` | Z = 1 + 5.50 X^−0.632 + i 8.43 X^−0.632 | σ |
| `variable_porosity` | Z = (1+i) sqrt(1000 σe/(π γ ρ0 f)) + i c0 αe/(8 π γ f), γ = 1.4, ρ0 = 1.2 | σe, αe (1/m, default 0) |
| `delany_bazley_layer` | Z = Z_c coth(−i k d), DB Z_c and k | σ, d (default 0.05 m) |

— `N/src/ground_effect.cpp` lines 114–170; `N/src/ground_effect.hpp` lines 40–74, 86–115
- A non-passive result (Re β ≤ 0) is refused; a thin low-resistivity DB layer does this at low frequency (0.03 m at σ 100 below about 60 Hz). — `N/src/ground_effect.cpp` lines 159–166; `P/docs/ground_plane_corrections.md` "Validation"
- Header comment: "Delany-Bazley is poor for grassland (their Table 2); on the 2017 Amedee data variable porosity with sigma_e ~ 200 fitted best." — `N/src/ground_effect.hpp` lines 64–66

**Turbulence and roughness**
- Chessell: T = exp(−(C/2 · f · √R)²), C = `--turbulence` (0 to 16e-4 rad·s/√m). HARMONOISE: T = exp(−(3/8)·0.364·k²·γ_p^(5/3)·R·γT), γ_p = h_s h_r/(h_s + h_r), `--turbulence_model harmonoise --gamma_t`. — `N/src/ground_effect.hpp` lines 76–84; `N/src/ground_effect.cpp` lines 212–228
- HARMONOISE coherence is exactly 1 when h_r = 0 (γ_p = 0; the code also skips it unless h_s + h_r > 0 and γT > 0). — `N/src/ground_effect.cpp` lines 222–228; `P/docs/ground_plane_corrections.md` "Validation" table ("a receiver on the ground gives 1")
- Roughness: coherent reflection × exp(−2 (k σ_h cos θ)²) (Kirchhoff/Ament); "the scattered part is taken as lost. It acts near normal incidence and vanishes toward grazing." `--roughness`, default 0. — `N/src/ground_effect.hpp` lines 107–110; `N/src/ground_effect.cpp` lines 229–234
- Turbulence and roughness are not applied to the plate receiver ("its reflection is millimeters from its direct path"). — `N/README.md` "Ground-plane microphones"; `N/src/plate_receiver.hpp` lines 30–31

**Receiver height 0**
- Header: "Zero puts the microphone on the surface, where the direct and reflected paths coincide and the model tends to its +6 dB limit" (rigid ground). — `N/src/ground_effect.hpp` lines 96–99
- In the code the sinc factor is applied only when the path delay is strictly positive ("at zero delay the two paths are the same path and the ratio is 0/0"). — `N/src/ground_effect.cpp` lines 247, 285–292
- In a run, a microphone of height 0 with `--plate_table` is a plate; any other height with `--ground_effect` is a point receiver at that height. — `N/README.md` "Receivers"; `N/src/niceops.cpp` line 849
- The harness's board-free comparison puts the microphone 1 mm up (`FLUSH_FT = 0.001 / 0.3048`), described as "1 mm up, which NICEOPS's point receiver needs" / "whose point receiver cannot sit on the plane". — `H/propagation/footprints.py` lines 48–53; `H/propagation/installation_validation.py` docstring; `H` git log 408ca6a

**Curved-ray geometry into the ground model**
- "Ground models get an equivalent straight geometry. A point receiver gets h' = L2 sin(grazing) − h_r and d' = L2 cos(grazing) from the reflected ray; a plate gets the direct ray's arrival angle and length. The ground models themselves are untouched." — `N/CLAUDE.md` (`src/refraction` entry)
- The reflected ray goes to the microphone's height above the plane z = 0; the ground is the plane z = 0 whatever the observer's own z. — `N/src/refraction.hpp` lines 189–196
- One direct and one reflected ray; "The extra ground-bouncing arrivals of strong downward refraction at long range are not included." — `N/README.md` "Refraction"
- The roadmap records a caveat on this approach: Nguyen et al. 2022, citing Ostashev & Wilson, "argue the spherical-wave coefficient's assumptions weaken under refraction"; Nord2000 is validated this way to 1 km and acceptable to 3 km (Ahearn et al. 2012). — `N/NICEOPS Atmospheric and Urban Propagation Options.md` line 100
- The sound speed for the ground impedance and reception time is the profile's surface level. — `N/README.md` "Stratified atmosphere"

**Tabulation**
- By default the point-receiver term is interpolated from a per-trajectory `GroundTable` over s = sqrt(h_s/r) (641 nodes, "s resolves grazing incidence, where the model changes fastest") and ln r (40 per decade); level histories agree with the exact model to 0.001 dB. `--exact_ground` bypasses it. — `N/src/ground_effect.hpp` lines 169–187

**Parameters used in the validation runs**
- Every harness run with a ground uses `--ground_effect --impedance_model variable_porosity --flow_resistance 200`, plus `--plate_table plates/gr1425_vp200.nc` for the plate; no `--porosity_rate`, `--gamma_t`, `--turbulence` or `--roughness` is passed. — `H/propagation/profile_validation.py` lines 58–59; `H/propagation/footprints.py` lines 49–50
- panam's site ground: `SITE_GROUND = dict(model='variable_porosity', sigma_e=200.0, alpha_e=0.0)`. — `P/noise_abatement_2017.py` line 81
- `--shadow_floor` default 20 dB; `--scattering` default 0 (off); `--receiver_height` default 3.937 ft (1.2 m). — `N/src/niceops.cpp` lines 2101–2116

### Inferences
- **Local reaction throughout.** Every impedance model enters as one angle-independent surface admittance β per band; the hard-backed layer is evaluated as a normal-incidence surface impedance. A search of `N/src` and `N/README.md` for "extended reaction" finds nothing; the only matches for "locally reacting" are the PE header (`N/src/pe.hpp` line 6) and the benchmark description. No extended-reaction or slit-pore/Wilson model exists in the code.
- At h_r = 0 the coded expression reduces to |1 + Q|² (ρ = 1, τ = 0, sinc skipped), so a flush point receiver in soil is governed entirely by Q at grazing, i.e. by β and the numerical distance. This is consistent with the memory note that "the level goes as about −40 log|β|" (`M/propagation-roadmap.md`, methods survey).
- Because Q is taken at the band center, any rapid variation of Q across a band near grazing is not averaged in the point model; the plate path found a linear-in-frequency Q inadequate (up to 0.8 dB above −15 dB, 2.9 dB in grazing nulls) and uses a quadratic (section 2).
- The CLI help for `--ground_effect` reads "(the default)" and the header says the model "is wired into the spectral path (the default since the 2017 approach validation)", while the README presents `--ground_effect` as an option to add. The harness always passes it explicitly, so the ambiguity does not affect the validation numbers.

### Gaps
- No measured ground impedance exists for the site (section 5), so σe = 200, αe = 0 is a fit to co-located pole/board pairs, not a measurement.
- No file states the validity range in elevation of the Weyl–van der Pol form as used here; the only quantified check is in panam against the Sommerfeld integral at plate scales (off by up to 6% at kR 0.3–30, 2% at kR > 10; agreement to 0.005 in Q at 500 ft). — `P/docs/ground_plane_corrections.md` "Validation"

---

## 2. Plate (flush ground-board) receiver (`--plate_table`)

### Takeaway
The plate receiver applies a precomputed axisymmetric-BEM solution for the GRAS GR1425 plate (400 mm, 8 mm thick, microphone at 3/4 R) as p/p_free = P_d + Q·P_r, with Q recomputed per geometry from the same ground model, averaged in energy over five sub-frequencies per band. Below about 1 kHz it collapses onto the soft-ground reflection, so at grazing it inherits the ground model's behavior; it matches the co-located pairs to 1.3–1.5 dB at 2–10° but fails at 0–2° (6.8 dB).

### Cited Findings

**What it computes and its inputs**
- p/p_free = P_d + Q P_r; P_d and P_r depend on frequency, source elevation and azimuth (far field of the plate); Q is the spherical-wave coefficient at the geometry. The table records the ground it was built on and a mismatch is refused. — `N/src/plate_receiver.hpp` lines 3–23
- Band average: |P_d + Q P_r|² at five sub-frequencies per third-octave band, averaged in energy. P_d, P_r interpolated bilinearly (real and imaginary) in elevation and azimuth; "below the table's lowest elevation they are held." — `N/src/plate_receiver.hpp` lines 25–28
- Run path: Q taken as quadratic in log frequency through each band's edges and center (2n + 1 evaluations). Within 0.07 dB of the exact form wherever the level is above −15 dB re free field; up to ~1.5 dB in the deep near-grazing cancellations (−30 to −70 dB). A linear Q was rejected: up to 0.8 dB above −15 dB and 2.9 dB in the grazing nulls, "because P_d − P_r nearly cancels there and any error in Q shows." — `N/src/plate_receiver.hpp` lines 76–82; `P/docs/ground_plane_corrections.md` "In NICEOPS"
- Shipped table `plates/gr1425_vp200.nc` (ncdump header): 31 bands × 5 sub × 80 elevations × 36 azimuths; `radius_m = 0.2`, `thickness_m = 0.008`, `edge_thickness_m = 0.0025`, `taper_length_m = 0.02`, `mic_x_m = 0`, `mic_y_m = 0.15`, `sound_speed_mps = 342.9`, `ground_model = variable_porosity`, `flow_resistance = 200`, `porosity_rate = 0`. Elevation nodes start at 0.05° and are log-spaced up to 10°, then roughly 2° steps to 90°. — `N/plates/gr1425_vp200.nc` (ncdump -h / -v elevation)
- A second table, `gr1425_db225.nc`, is for Delany–Bazley 225. — `N/README.md` "Ground-plane microphones"
- Azimuth: horizontal propagation direction in the plate frame; microphone offset "outboard" (away from the track) by default, `--plate_mic_azimuth`. — `N/src/plate_receiver.hpp` lines 155–164
- Q is taken at the plate top (`thickness_m`, "the plate top above the ground, where Q is taken"). — `N/src/plate_receiver.hpp` line 37
- The BEM itself (panam `axisymmetric_bem.py`): rigid body (flat top, conical taper, vertical rim) on a soft half-space, exterior Neumann boundary integral with the exact half-space Green's function (complex-image form, Ochmann/Taraldsen), azimuthal modes |m| ≤ ka + 10, ~60 generator segments. — `P/docs/ground_plane_corrections.md` "Plate thickness", "Exact ground Green's function", "Axisymmetric, compiled BEM"

**Behavior versus elevation near grazing**
- "Below about 1 kHz the plate is too small to matter, and every model but flat and De Jong collapses onto the soft-ground reflection. Above 2 kHz at 10° and up the plate restores the doubling that soft soil does not." — `H/README.md` "Ground-board models against the co-located pairs"
- K, the plate over bare ground, near grazing: within 0.1 dB from 0.5° to 10°; about 0 dB up to 1 kHz, +0.9 dB at 2 kHz, +2.6 at 4 kHz, +7 at 8 kHz; varies with azimuth by ±1 dB at 1–4 kHz and ±2.5 dB at 8 kHz. — `N/src/refraction.hpp` lines 343–348; `N/README.md` "Refraction through the profile's layers"
- Plate thickness is second order: real plate re a 0.3 mm plate is within about ±1 dB everywhere; at 1–10° elevation −0.6/+0.3 dB at 1 kHz and +0.8/+0.4 at 2 kHz (toward/away from the microphone side, exact kernel). "It cannot explain the multi-dB 0-2 deg conflict." — `P/docs/ground_plane_corrections.md` "Thickness effect with the exact kernel"
- The real plate on rigid ground is not exactly 2: −0.35/+0.2 dB at 500 Hz and about ±1.2 dB at 2 kHz. — `P/docs/ground_plane_corrections.md` "Axisymmetric, compiled BEM"
- Measured deviation of the board from +6 dB (first pass, pooled bands, DB σ 100 pole model): −10.6, −6.6, −3.7, −2.0 dB at 0–2, 2–5, 5–10, 10–20°. With the refitted pole model the measured level near grazing rises by up to about 8 dB and the ideal rigid plane's 0–2° error falls from 12.7 to 3.5–5 dB. — `P/docs/ground_plane_corrections.md` "First comparison, all runs" and "Path 4"
- In SEL the plate and the flush-in-soil reflection differ by 0.1 dB ("A-weighted helicopter SEL sits below 2 kHz"). — `H/README.md` "Footprints, four runs by five models"

**Validation**

| component | reference | result |
| --- | --- | --- |
| exact half-space Green's function | Sommerfeld integral by quadrature, 15 cases, 3 grounds, kR 0.3–30 | 2e-9 relative |
| axisymmetric BEM | exact series, rigid sphere (hemisphere on rigid ground), ka 0.5, 2, 5 | max 2.7e-4 relative (0.002 dB) |
| three plate formulations in the thin limit | each other, 250–1000 Hz, 5–70° | within 0.25 dB |
| real plate, axisymmetric vs 3-D surface | each other, 500 Hz and 2 kHz, 10° and 45° | within 0.3 dB |
| NICEOPS `level_exact` | panam `board_level` | 1e-6 dB |

— `P/docs/ground_plane_corrections.md` "Validation", "In NICEOPS"; `N/README.md` "Ground-plane microphones"
- "The combined soft-ground plate problem has no closed form"; published BEM correction curves (Kingan et al., Blandeau et al.) "were not accessible". — `P/docs/ground_plane_corrections.md` "Validation"

**Against the co-located pole/board pairs** (four aircraft, 771 runs, 68,210 frames; every model on variable porosity σe 200; rms of the per-band median error, 50 Hz–2.5 kHz, dB)

| model | 0–2° | 2–5° | 5–10° | 10–20° | 20–40° | 40–90° |
| --- | --- | --- | --- | --- | --- | --- |
| flat +6 dB | 6.15 | 5.81 | 3.91 | 1.89 | 1.13 | 0.86 |
| soft-ground reflection, flush in the soil, no board | 6.74 | 1.59 | 1.37 | 0.61 | 0.60 | 0.81 |
| Fresnel zone, strip | 6.71 | 1.54 | 1.31 | 0.44 | 0.42 | 0.66 |
| Fresnel zone, 400 mm disc | 6.79 | 1.66 | 1.38 | 0.53 | 0.36 | 0.63 |
| De Jong (nMID) edge diffraction | 6.48 | 2.73 | 1.95 | 0.76 | 0.56 | 0.72 |
| BEM, thin disc | 6.74 | 1.41 | 1.24 | 0.34 | 0.32 | 0.56 |
| BEM, real plate | 6.76 | 1.48 | 1.27 | 0.40 | 0.20 | 0.36 |

— `H/README.md` "Ground-board models against the co-located pairs"; run log `H/propagation/ground_models.log` (68,210 frames; pair gain offsets 50: +0.94, 51: −0.73, 52: −0.48 dB)
- Caveat stated there: "The pole's model was fitted on half these runs at 30° and up, so the steep bins are partly in-sample." — `H/README.md`
- The measured transfer function near grazing rests on the pole's modeled ground reflection (−10 to −14 dB at 125–500 Hz below 2°). — `P/docs/ground_plane_corrections.md` "Literature review", implication 4

### Inferences
- Since the plate equals the soft-ground reflection below ~1 kHz, the plate model adds nothing new at grazing in the bands that carry A-weighted helicopter SEL; whatever is wrong at 0–2° is in Q (ground, refraction, or the reference), not in the plate geometry. The panam log reaches the same conclusion ("That points past the plate, to refraction or a range-dependent ground").
- The plate table's Q uses the plate-top height (8 mm) and no turbulence/roughness, so the HARMONOISE γT = 3e-5 fitted for the pole has no effect on any flush-microphone prediction.

### Gaps
- The test report (Watts et al., NASA TM 2019) says a 15 inch (381 mm) board; the user confirmed GRAS GR1425 (400 mm) lying on the ground. The note flags the report as wrong; the discrepancy rests on the user's statement. — `P/docs/ground_plane_corrections.md` "The pole microphones' response" and "Plate thickness"
- Taper length (20 mm) is assumed, not on the GRAS drawing (checked at 10 and 30 mm: up to ~0.6 dB at high frequency). — `P/docs/ground_plane_corrections.md` "Plate thickness"
- About +2 dB of excess at 2–5 kHz at 40–90° in the pair data remains after the 46AE incidence response; roughness σ_h = 8 mm (later 10–15 mm) reduces it. Not a grazing issue but it shows the pole reference is imperfect. — `P/docs/ground_plane_corrections.md` "The pole microphones' response", "Orientation, windscreen, roughness"

---

## 3. Refraction, shadow and scattering near grazing; atmosphere below the lowest measured level

### Takeaway
Refraction is by eigenrays (`--refraction linear`: one equivalent linear profile per path; `stratified`: piecewise-linear layers, four ray families), keeping **one direct and one reflected arrival**, with a +10 dB caustic cap, a hand-over to the linear fit's ray at 0.9 of the layered shadow boundary, a tabulated residue-series transition through the shadow boundary with a per-band effective boundary, and either a 20 dB shadow floor or a fitted scattering term. In the harness the wind below the lowest LIDAR gate (46 ft) is an assumed log law with z0 = 0.1 ft; the measurements support γT ≤ ~3e-6 for scattering, a tenth of the pole-pair fit.

### Cited Findings

**Arrivals, caps, hand-over**
- `linear`: "One direct and one reflected ray. The extra ground-bouncing arrivals of strong downward refraction at long range are not included." "One gradient per path. A log-law profile's steep lowest meters are averaged into it." Equivalent-linear models up to 10 dB off a PE near the ground at 500–2000 Hz (van Maercke & Salomons 2009). The fit is least squares over the heights the path spans, widened to at least 10 m. — `N/README.md` "Refraction"
- The effective-sound-speed approximation "overstates the wind by 1/cos(elevation)". — `N/README.md` "Refraction"
- `stratified`: "One arrival per pair. Non-monotonic profiles give caustics and several arrivals; the solve returns the first family's root. Focusing is capped at `max_focusing_gain_db` (+10)." — `N/CLAUDE.md` (`src/stratified` entry); `N/src/stratified.hpp` line 46
- Effective sound speed is linear between the profile's own levels and "held constant below the lowest and above the highest". — `N/src/stratified.hpp` lines 8–10; `N/README.md` "Refraction through the profile's layers"
- Hand-over: past 0.9 r_s (`handover_fraction`) the lit ray is the linear fit's, kept past the fit's own boundary too. Short of r_s the layered ray focuses into a caustic (+10 dB) the PE does not show. 0.9 gave 0.6 dB rms on day 284 over [0.5, 1] r_s (1.4 at r_s) and 1.7 on log laws (3.3). — `N/CLAUDE.md`; `N/README.md` "The measured profiles of the 2017 test's windy day"; `N/src/stratified.hpp` line 117
- A microphone below the plane is solved as if on it (ground boards 3–10 cm under the survey datum gave a ray grazing under the plane, tube 46× its range, 33 dB low). — `N/CLAUDE.md`
- Ground-grazing boundary rule: `fastest_above(z_s) > fastest(z_r, z_s)`; day-284 mornings have 330.5 m/s aloft against 339 at the ground. — `N/CLAUDE.md`

**Shadow loss and transition**
- Shadow loss: grazing ray's level less the least-attenuated creeping-wave mode, capped at `--shadow_floor` (default 20 dB). Over rigid ground about 24 dB/km at 1 kHz for a gradient of 0.02 /s; with `--ground_effect` the rate runs over the ground's admittance. Ordinary ground attenuates up to 12% less than rigid at small |q| and up to 2.3× more once soft on the wave's scale (benchmark 1 kHz ground at 0.1 /s: 146 dB/km vs rigid 70). Bands with Re q > 0 (arg β < −60°) take the rigid rate. Rate held to the residue series' slope to 0.02–0.3%. — `N/README.md` "Refraction"; `N/CLAUDE.md` (`ShadowAttenuation`)
- `ShadowTransition`: per band and microphone height, a table V(ξ, s) of the residue series relative to the ray's own spreading, ξ = (r − r_b)/F, F = 2 k0 l² (about 100–400 m). Ray model below ξ = −2, blend to −1, V beyond; floored at V(0) − floor, capped at 20 log10 2. Built at a reference source height 20 l. Ray theory's ground term at the grazing geometry was 10–25 dB low over grass and ~5 dB high over rigid ground near the boundary. — `N/CLAUDE.md`; `N/src/refraction.hpp` lines 296–326
- Per-band boundary through layers: r_b = r_s + w (r_l − r_s), w = x/(x + 0.363) at 1.2 m, 0.75 for flush; where the fit casts no shadow, 1.60 (f/1 kHz)^−0.08 · r_s flush and 1.90 (f/1 kHz)^−0.08 · r_s at 1.2 m. — `N/src/refraction.hpp` lines 452–457; `N/CLAUDE.md`
- Plates: field at the plate's top (8 mm) plus K per band and azimuth taken at 2° elevation. — `N/CLAUDE.md`; `N/src/refraction.hpp` lines 343–348
- "The law cannot be made much better by fitting": per-case optimum 2.7 dB pooled; every smooth law tried sits at 6–9 dB. "Do not spend another day on it without turbulence in the reference." — `N/CLAUDE.md`

**Shadow-edge accuracy against the (still-air) PE** (rms from the PE's last level above −3 dB down to −25 dB)

| set | before | now |
| --- | --- | --- |
| day-284 profiles, 18 cases (sources 50, 120, 250 m), flush | 12.9 | 8.1 |
| log laws, 16 cases, 1.2 m | 6.2 | 6.0 |
| log laws, 16 cases, flush | 7.3 | 7.6 |
| linear lapse (−0.018 /s), sources 30 and 100 m, 125 Hz–2 kHz | 5–20 | 0.1–1.0 |

— `N/README.md` "Refraction through the profile's layers" and "Checking the refraction models against a parabolic equation". Rerun at 20 points/wavelength: 8.14/8.16, 6.04/6.00, 7.60/7.57. — `H` git log 085b6c1
- Flush (3 mm) log-law comparison by band: 125 Hz worst (10.2 dB), 500 Hz best (4.4). — `N/README.md`
- Lit-region/downwind: `stratified` within 0.4–2.1 dB rms of the PE beyond 400 m (100 m source, `example_windy.csv`: 0.9, 0.4, 0.9, 2.1 dB at 125, 250, 500, 1000 Hz); straight 2.8–8.1, linear 0.9–7.6. — `N/README.md` "Checking the refraction models against a parabolic equation"
- Upwind shadow onset (−10 dB re free field, 100 m source, 1.2 m receiver): PE 1.3, 1.07, 0.91, 0.84 km at 125–1000 Hz; stratified 1.2, 0.53, 0.53, 0.53; linear 1.9, 1.4, 1.3, 1.4. — `N/README.md`
- "Without `--gamma_t` the deep shadow, 40–90 dB down, is not physical." — `N/CLAUDE.md`

**Scattering (`--scattering`, off by default)**
- S = a + 10 lg γT + b lg(f/1 kHz) + c lg(R/100 m) + d lg(z_s/100 m) dB re free field. Flush (3 mm or plate): a 39.18, b −1.55, c 2.10, d 5.30 (12 profiles, sources 50 and 120 m, 1.5 dB rms). 1.2 m: a 36.67, b +5.61, c 5.13, d 5.30 (the flush fit's; 6 profiles, sources 50 m). S held to ≤ 0 dB, R 100 m–10 km, z_s 10–300 m. — `N/README.md` "Turbulence, and the energy it scatters into shadows"
- PE turbulence: von Kármán phase screens, C² = γT/4, outer scale default 10 m, uniform in height, one screen per outer scale of range. Single scattering linear in γT confirmed (×3 → 4.82 dB, theory 4.77). — `N/CLAUDE.md`; `N/README.md`
- Evidence on strength, day-284 upwind 5–10° bin (median / rms): floor +0.28 / 3.06; γT 1e-6 +0.49 / 3.04; 3e-6 +0.95 / 3.02; 1e-5 +1.55 / 3.16; 3e-5 +2.96 / 4.00. "Up to about 3e-6 is consistent with them, which is HARMONOISE's 'moderate' range." At 3e-5 the PE fills day-284 shadows to −3 dB re free field at a board (median). — `N/README.md`; `H/README.md`
- Moving `--shadow_floor` from 20 to 50 dB "changed nothing" for the day-284 upwind errors. — `H/README.md` "Three faults in `stratified`'s shadow"
- Center-line flyover frames, day 284 (one-bearing comparison): on 284211 (45 m) the still-air PE puts the shadow onset within about 50 m of the measured one in every band 250 Hz–2 kHz; only turbulence (γT 3e-6) gets the depth, −17 to −18 dB against −12 to −17 measured, where still air falls below −35. Measured minus turbulent PE: +1.3 to +2.6 dB mean, 5.5–6.3 dB rms. On 284153 (76 m) there is no shadow, "but the measured level falls 4-6 dB/km faster than every model predicts." — `H` git log 085b6c1; `H/propagation/flyover.log` (1,835 and 1,705 frame-bands from 8 microphones)

**How the harness builds the profile below the lowest measured level**
- Temperature: tethered balloon, 0–500 ft AGL every 10 ft (°F), held at the 500 ft value above. Wind: ZephIR 300 at 46, 125, 197, 348, 495, 696, 984 ft AGL, taken as knots, direction "from", true; interpolated in east/north components, held above 984 ft, "and below 46 ft brought to zero at the ground by a log law with z0 = 0.1 ft". Humidity: ground stations', same at every height. Pressure: stations', hydrostatic above. — `H/propagation/profile_validation.py` docstring lines 14–22
- Code: levels include 0, 1, 3, 10, 20 ft plus the temperature and LIDAR heights and a 5,000 ft top; below the lowest valid gate, `scale = log(z/Z0_FT)/log(h0/Z0_FT)` applied to that gate's u, v (so direction is held constant below 46 ft); `Z0_FT = 0.1`, `TOP_FT = 5000.0`. A run needs at least three valid LIDAR gates. — `H/propagation/profile_validation.py` lines 60–62, 132–160
- NICEOPS then holds the lowest level's values below it and refuses trajectories above the top. — `N/README.md` "Stratified atmosphere"
- Roadmap: "The biggest risk is not the method but the input: the near-ground gradient controls the grazing geometry, and it is rarely known well." — `N/NICEOPS Atmospheric and Urban Propagation Options.md` line 168

### Inferences
- The shadow boundary under `stratified` is set by the speed at the ground against the speed at the source, and in the harness the ground-level wind is not measured: it is the assumed log law between 0 and 46 ft. The calibration README's finding that the upwind limiting ray "is fixed by c(z_s)/c(0) alone" (`N/CLAUDE.md`, "The log-law finding") means the unmeasured lowest 14 m directly sets the layered shadow onset.
- The 1.2 m scattering constants come from one source height (50 m) and borrow the flush fit's height term; they are less constrained than the flush ones.
- The scattering term and the shadow transition need `--refraction` and `--ground_effect`; straight-ray runs have neither, so the calm-day straight-ray residuals contain no turbulence treatment for flush microphones at all.

### Gaps
- No turbulence measurement (C_T², C_v²) is recorded anywhere in the files read; γT is only inferred (pole-pair fit 3e-5; SEL consistency ≤ 3e-6).
- Multiple arrivals under downward refraction (the calm-day inversion case) are not modeled and no file quantifies what they would add.
- The per-band boundary law was calibrated against a PE without turbulence.

---

## 4. Quantified residuals against measurements at low elevation, long range, low frequency, or flush microphones

### Takeaway
On calm day 286 with the 2° spheres and like-for-like 10 dB-down SEL, the 5–10° bin is still +1.4 to +2.0 dB high depending on the model (62 mic-runs), while frame-level analysis shows a range trend of −2.4 dB/km at fixed direction, with NICEOPS about 2 dB too quiet beyond 4,000 ft at 2–10° and 63–630 Hz decaying 1–2 dB/km more slowly in the measurement than predicted. On windy day 284, `stratified` brings 5–10° to about +0.3 to +0.8 / 3.0–3.6 dB. At 0–2° the flush microphone reads 5–10 dB above every ground-aware model from 500 Hz up.

### Cited Findings

**4a. Day 286 (calm, 1–4 kt), 25 B407 approaches, 1,112 mic-runs, plate + VP 200, SEL by peak elevation (median / rms)**

Like-for-like measured SEL (`--measured standard`, 10 dB-down first-to-last crossing). No mic-run peaks below 5° (n = 0 in 0–5°).

| database | model | 5–10° (n=62) | 10–20° (n=211) | 20–40° (n=253) | 40–90° (n=586) | all |
| --- | --- | --- | --- | --- | --- | --- |
| sep25 (original, min elev 10°) | station | +2.83 / 3.54 | +0.72 / 1.45 | +0.24 / 0.92 | +0.35 / 1.35 | +0.44 / 1.50 |
| min10 (control rebuild) | station | +2.75 / 3.49 | +0.69 / 1.44 | +0.24 / 0.92 | +0.33 / 1.35 | +0.41 / 1.49 |
| min5 | station | +1.51 / 2.11 | +0.33 / 1.28 | +0.16 / 0.92 | +0.31 / 1.35 | +0.29 / 1.31 |
| min2 (promoted) | station | +1.43 / 2.05 | +0.32 / 1.28 | +0.10 / 0.92 | +0.29 / 1.35 | +0.27 / 1.30 |
| min2 | profile (absorption through profile) | +1.45 / 2.05 | +0.45 / 1.33 | +0.22 / 1.03 | +0.40 / 1.47 | +0.40 / 1.39 |
| min2 | linear | +1.57 / 2.09 | +0.72 / 1.50 | +0.22 / 1.04 | +0.42 / 1.48 | +0.47 / 1.44 |
| min2 | stratified | +1.63 / 2.15 | +0.59 / 1.42 | +0.23 / 0.97 | +0.39 / 1.45 | +0.43 / 1.40 |

— `H/propagation/database_effect_B407_standard.txt`
- After the plate boundary field and sub-plane fix (`out_plate`), 5–10°: linear +1.25 / 1.82, stratified +1.98 / 2.62 (overall 1.42 and 1.46). — `H/propagation/day286_B407_plate.txt`. After the three shadow fixes, stratified 5–10° rms 2.62 → 2.55; overall 1.46 → 1.45. — `H/README.md`; `N/README.md`
- Against the whole-reception-window measured SEL (the earlier default), the same runs read lower: min10 station +2.42 / 3.23, min5 +1.21 / 1.89, min2 +1.14 / 1.83 at 5–10°; overall 1.39, 1.25, 1.25. The whole window holds 0.4–1.2 dB more than the 10 dB-down SEL, most at 5–10°. — `H/propagation/database_effect_B407.txt`; `H/README.md` "The default measured SEL flatters low-elevation over-prediction"
- Earliest run (17 approach runs, 756 mic-runs, original spheres): baseline median −0.95, rms 1.60; station −0.01, rms 1.41; linear rms 1.47; stratified 1.44; "+2.7 dB over-prediction at 5–10° peak elevation" against "+0.3 at 10–20°". The −0.95 baseline bias was the 20% humidity default. — `M/propagation-roadmap.md` "Measured validation (2026-09-27)"
- All-aircraft, flat-correction era (6,118 mic-runs, standard 10 dB-down SEL): +4.8 dB (median) below 5°, +2.9 dB at 5–10°, within about 1 dB above 10°. — `P/docs/ground_plane_corrections.md` "Why this matters"
- Footprints on two calm runs (flat / EGA / plate / linear+plate / stratified+plate): 286270 +2.0/2.6, +0.8/1.3, +0.6/1.3, +0.6/1.3, +0.5/1.3; 286277 +2.2/3.5, −0.1/1.6, −0.1/1.6, 0.0/1.6, 0.0/1.9. "The flat correction over-predicts every sideline by 3–8 dB." — `H/README.md`

**4b. Day 286 frame-level residuals** (25 approaches, `stratified` + plate + scattering 3e-6, 0.5 s A-weighted window frames, 65,634 frames; median error by emission elevation and range, n in brackets)

Forward lines (x > −500 ft; the lines at x = 0 and +608 ft):

| emission elevation | 0–1000 ft | 1000–2000 ft | 2000–4000 ft | 4000–8000 ft | 8000–20000 ft |
| --- | --- | --- | --- | --- | --- |
| 2–5° | +2.84 (1400) | +0.55 (501) | +0.35 (929) | −3.06 (98) | |
| 5–10° | +1.64 (807) | +1.50 (1672) | +0.75 (1326) | −1.31 (1021) | −4.63 (116) |
| 10–20° | +1.24 (52) | +2.01 (323) | +0.79 (474) | −0.54 (231) | −13.28 (111) |

Other lines:

| emission elevation | 0–1000 ft | 1000–2000 ft | 2000–4000 ft | 4000–8000 ft | 8000–20000 ft |
| --- | --- | --- | --- | --- | --- |
| 2–5° | +0.27 (49) | −0.08 (649) | +0.11 (1578) | −1.95 (185) | |
| 5–10° | −0.48 (308) | +0.02 (3480) | +0.13 (1498) | −1.91 (859) | −0.15 (109) |
| 10–20° | −0.10 (2647) | +0.07 (5937) | −0.04 (4660) | +0.54 (1502) | |
| 20–40° | +0.15 (5655) | +0.51 (7869) | +1.02 (3618) | +0.65 (773) | |
| 40–90° | +0.60 (6680) | +1.19 (4667) | +1.52 (819) | | |

— `H/propagation/band_frames/aweighted.txt`
- Range slope at fixed direction cell (elevation × sphere θ, 10° × 30°) and run: all −2.40 ± 0.02 dB/km (n = 65,634); forward lines −2.69 ± 0.04 (n = 9,970); other lines −1.99 ± 0.03 (n = 55,664). Standard errors treat frames as independent, "so it is optimistic by the frames' correlation within a mic". — `H/propagation/band_frames/aweighted.txt`; `H/propagation/band_frames_summary.py` docstring
- "beyond 4,000 ft at 2–10° NICEOPS is about 2 dB too quiet on every line. Per band, 63–630 Hz decays 1–2 dB/km more slowly than predicted ... Above 3 kHz the measured frames at range are noise-limited (per-band slopes of −13 to −39 dB/km at 4–5 kHz)." "The state is matched by height, which is exact only for a monotonic descent." — `H/README.md` "The forward-line over-prediction is the end of the approach"
- Forward lines carry about 30% of the squared SEL error on both days: +1.65 dB median on day 286 against +0.34 on the other lines; the frames carry +1.46. By height at emission (forward / other): 0–100 ft +2.04 / +0.07; 100–200 ft +1.23 / −0.18; 200–400 ft +1.15 / +0.10; 400–800 ft −0.02 / +0.41; 800–3000 ft −2.27 / +0.84. By ground speed: 0–20 kt +3.19 / −4.68; 20–40 kt +0.97 / −0.33; 40–60 kt +1.23 / +0.02; 60–80 kt +1.29 / +1.42; 80–200 kt −0.89 / +0.04. — `H/propagation/band_frames/aweighted.txt`; `H/README.md`
- By sphere θ at emission, 5–20° elevation: θ 0–30 forward +0.43 / other −0.47; 30–60 +1.78 / +0.96; 60–90 −0.02 / +0.20; 90–120 other −0.27; 120–150 −0.47; 150–180 −1.08. — `H/propagation/band_frames/aweighted.txt`
- Frame gating: measured band ≥ 10 dB over the mic's ambient and ≥ 3 dB over the record's own floor, within 10 dB of that mic-band's loudest frame, emission elevation ≥ 2°. — `H/propagation/band_frames_summary.py` docstring
- Memory note on interpretation: "The long-range under-prediction on the calm day is separate; the first suspect is the calm-day inversion refracting sound downward." The user rejected raising panam's 2,000 ft sphere range cap. — `M/forward-line-error-source.md`

**4c. Day 284 (windy, 9–15 kt; LIDAR median 11.0 kt), 16 B407 steady runs, held-out spheres (`min2_no284`), plate, standard SEL, 738 mic-runs**

| stage | model | 0–5° (n=2) | 5–10° (n=71) | 10–20° (n=131) | 20–40° (n=155) | 40–90° (n=379) | all |
| --- | --- | --- | --- | --- | --- | --- | --- |
| in-sample spheres (`out_284_full`) | station (straight) | +8.88 / 8.90 | +3.36 / 4.48 | +0.80 / 1.89 | +0.15 / 2.21 | −0.09 / 1.59 | +0.22 / 2.26 |
| in-sample | linear | +4.53 / 4.57 | +2.26 / 3.22 | +0.88 / 1.84 | +0.24 / 2.25 | −0.09 / 1.59 | +0.20 / 2.00 |
| held out (`out_284_no284`) | station | +11.31 / 11.31 | +4.52 / 5.76 | +0.94 / 2.09 | +0.21 / 2.36 | −0.27 / 1.78 | +0.22 / 2.67 |
| held out | linear | +5.73 / 5.73 | +3.53 / 4.37 | +1.02 / 2.09 | +0.36 / 2.40 | −0.29 / 1.78 | +0.20 / 2.36 |
| held out | stratified, before fixes | −25.41 / 25.42 | −0.68 / 11.17 | +0.99 / 2.21 | +0.44 / 2.40 | −0.27 / 1.79 | −0.01 / 4.18 |
| held out, sub-plane fix only | stratified | +0.95 / 1.06 | +0.01 / 4.51 | +1.25 / 2.27 | +0.44 / 2.40 | −0.27 / 1.79 | +0.11 / 2.39 |
| held out, + plate boundary field (`_fixed`/`_all`) | stratified | +1.01 / 1.14 | +0.26 / 3.58 | +1.14 / 2.27 | +0.44 / 2.40 | −0.27 / 1.79 | +0.14 / 2.24 |
| held out, + plate boundary field | linear | +6.49 / 6.50 | +3.70 / 4.50 | +0.97 / 2.07 | +0.36 / 2.40 | −0.29 / 1.77 | +0.21 / 2.38 |

— `H/propagation/day284_B407.txt`, `day284_B407_attribution.txt`, `day284_B407_fixed.txt`
- Holding day 284 out of the spheres raises the 5–10° error by about 1.2 dB (straight +3.36 → +4.52), "which is the size of that absorption" of refraction by in-sample spheres. — `H/README.md`; `H/propagation/day284_B407.txt`
- The day puts 0.8% of evaluations in shadow under `linear` and 16.6% under `stratified`. — `H/README.md`
- Later (2026-09-29, three shadow fixes + flush boundary law), `stratified` 5–10°: +0.26 / 3.58 → +0.79 / 2.99; overall 1.85 → 1.70 (cut records excluded, 729 mic-runs). Upwind only: 5–10° (65) +0.45 / 3.72 → +0.28 / 3.06; 10–20° (80) +0.99 / 2.24 → +0.71 / 1.87. — `H/README.md`; `N/README.md` "Refraction through the profile's layers"
- Excluding the nine cut day-284 mic-runs: overall rms 2.67 / 2.38 / 2.24 → 2.36 / 2.02 / 1.85 (straight / linear / stratified); the 5–10° bin does not move. — `H/README.md`
- Structure of the day-284 upwind errors before the fixes: "sideline mics 2–9 dB low under `stratified`, far mics 4–9 dB high under every model." — `H/README.md`
- Footprints: 284205 D34 (windy) flat +2.9/6.4, EGA +1.5/3.8, plate +1.5/3.7, linear +1.5/2.8, stratified +0.7/1.7; 284211 L4 +2.3/4.5, 0.0/1.9, +0.1/1.8, −0.4/1.2, −0.8/1.3. — `H/README.md`

**4d. Modeling-option ladder** (every B407 mic-run, overall rms, calm 286 / windy 284, dB; NICEOPS d852613)

| option | 286 | 284 |
| --- | --- | --- |
| broadband, EAA, flat +6 dB | 3.40 | 4.34 |
| spectral, default air (7 °C, 20% RH) | 2.32 | 3.55 |
| spectral, station air | 2.98 | 4.18 |
| soft-ground reflection (flush in soil, no board) | 1.38 | 2.43 |
| ground-board plate (BEM) | 1.30 | 2.36 |
| absorption through measured profile, straight rays | 1.39 | 2.37 |
| refraction, linear fit | 1.42 | 2.02 |
| refraction, layered | 1.45 | 1.70 |
| + scattering γT 3e-6 | 1.46 | 1.69 |

— `H/README.md` "What each modeling option is worth" (`H/propagation/options.log`, `options_ega.log` list the 41 runs). The ground correction takes the share of mic-runs within 2 dB from 54 to 91% on day 286.

**4e. Flush microphone versus pole at 0–2°** (co-located pairs; not an SEL residual but the direct test of the ground model at grazing)
- Every ground-aware model scores 6.5–6.8 dB rms at 0–2° (table in section 2). "Nothing works at 0–2°, where the measurement sits 5–10 dB above every ground-aware model from 500 Hz up. That regime wants a far harder ground than every other elevation, which points to refraction." — `H/README.md`
- Band-by-band (first pass): soft-ground and Fresnel models match only below about 300 Hz at 0–2° and below about 500 Hz at 2–5°; above that they predict −10 to −25 dB while the measurement "bottoms out near −10 dB around 300–600 Hz and rises again", residuals +10 to +15 dB by 1–2 kHz. The disc boundary integral "recovers from −17 to about −13 dB at 1–2 kHz where the data reach about −2 dB." Below 10° the SNR gate leaves few frames above 1 kHz (hundreds against thousands below 300 Hz). — `P/docs/ground_plane_corrections.md` "First comparison" and "Path 3"
- Held-out-run ground profile (thin disc, rms dB):

| ground | 0–2° | 2–5° | 5–10° | 10–20° | 20–30° | 30–40° | 40–90° |
| --- | --- | --- | --- | --- | --- | --- | --- |
| before (DB 100, no turbulence) | 6.18 | 4.21 | 2.89 | 1.59 | 0.78 | 0.79 | 0.75 |
| DB 200 | 8.29 | 2.74 | 1.82 | 0.81 | 0.44 | 0.51 | 0.55 |
| DB 400 | 6.49 | 1.66 | 1.47 | 0.50 | 0.32 | 0.41 | 0.53 |
| DB 800 | 4.62 | 2.06 | 1.74 | 0.58 | 0.40 | 0.53 | 0.66 |
| DB 3200 | 1.89 | 3.32 | 2.24 | 0.77 | 0.54 | 0.74 | 0.86 |
| VP σe 100 | 7.76 | 2.06 | 1.37 | 0.61 | 0.30 | 0.38 | 0.55 |
| **VP σe 200** | 6.50 | **1.28** | **1.17** | **0.37** | **0.27** | 0.43 | 0.60 |
| VP σe 400 | 5.18 | 1.71 | 1.52 | 0.47 | 0.35 | 0.53 | 0.69 |
| VP σe 1600 | 2.75 | 2.98 | 2.10 | 0.70 | 0.48 | 0.69 | 0.83 |
| real plate, VP 200 | 6.26 | 1.33 | 1.25 | 0.41 | 0.23 | 0.35 | 0.37 |
| real plate, VP 150 | 6.85 | 1.45 | 1.21 | 0.45 | 0.21 | 0.33 | 0.38 |
| real plate, DB 400 | 6.36 | 1.66 | 1.53 | 0.51 | 0.28 | 0.40 | 0.43 |

— `P/docs/ground_plane_corrections.md` "Pinning the ground" and "The ground profile redone with the real plate"
- "0-2 deg wants a far harder ground (sigma 3200-6400: 1.6-1.9 dB) than every other elevation, so it conflicts with them. Those frames are the lowest and longest-range." — `P/docs/ground_plane_corrections.md` "Pinning the ground"
- Steep-angle fit (≥ 30°, all models 1.39 dB mean abs residual): DB σ 1300, Miki σ 1700, VP σe 380 (αe → 0), DB layer σ 1300 depth 0.11 m, γT 3.0–3.2e-5; log10 σe = 2.61 ± 0.05, log10 γT = −4.57 ± 0.04. Adding the 46AE response moved σ from 1300 to about 5300 (DB) and 380 to 1900 (VP) with the loss unchanged (1.39 → 1.40), at a cost at 0–2° (2.9 → 4.5 dB): "The flow resistivity is not identified by the steep-angle fit." — `P/docs/ground_plane_corrections.md` "Path 4", "The pole microphones' response"
- Out-of-sample below 30° with those steep-fit grounds: DB 1300 + turbulence 2.91 at 0–2°, 2.66 at 2–5°; VP 380 + turbulence 4.69, 1.66; DB layer 0.11 m 2.90, 2.67. "Variable porosity is best from 2 to 40 deg. Delany-Bazley (or the layer) is best at 0-2 deg." — `P/docs/ground_plane_corrections.md` "Path 4"
- Installation effect in the sphere build (B407 run 283115, plate minus flat): +0.2 to +0.4 dB below 100 Hz; +0.8 to +2.6 dB at 0.5–2 kHz, largest toward the sidelines; within ±0.2 dB above 5 kHz. — `P/docs/ground_plane_corrections.md` "In the sphere builds"

**What has been ruled out** (as stated in the files)
- For the 5–10° over-prediction: refraction and absorption; the hover sphere (flooring speed at 37 kt changes nothing); QSAM (turning it off is worse); a frame mismatch. — `M/propagation-roadmap.md`
- The filled rim below the spheres' elevation floor explains about half: 10° → 5° floor takes 5–10° from +2.42 / 3.23 to +1.21 / 1.89 (whole-window measure); 2° gains only 0.07 dB more, "so the remaining +1.1 dB is not the filled rim." — `H/README.md`
- For the forward-line excess: long-range propagation and the sphere range cap (the excess is at short range, low and slow; the long-range residual has the opposite sign). — `H/README.md`; `M/forward-line-error-source.md`
- For the 0–2° flush residual: plate thickness (second order, ±1 dB), plate model choice (all models within 0.3 dB of each other there), microphone position. — `P/docs/ground_plane_corrections.md`
- For the day-284 upwind errors: the shadow floor (20 → 50 dB changed nothing). — `H/README.md`
- A Doppler frequency shift is correctly not applied (spheres have `DOPPLER_SHIFT_REMOVED = 0`). — `M/sphere-doppler-convention.md`

### Inferences
- The remaining calm-day 5–10° SEL excess (+1.4 to +2.0 dB) and the long-range deficit (−2 dB beyond 4,000 ft) have opposite signs and coexist in the same elevation band; the frame tables show the sign flips with range (positive inside 4,000 ft on forward lines, negative beyond on all lines). An SEL-by-peak-elevation statistic therefore mixes a short-range source/directivity excess with a long-range propagation deficit.
- The measured-decays-slower finding on calm day 286 (63–630 Hz) and the measured-falls-faster finding on windy run 284153 (4–6 dB/km faster than every model) point in opposite directions on different days, which is what a refraction-related (meteorology-dependent) cause would do and what a fixed ground-impedance error would not. This is an inference; no file states it.
- The 0–5° SEL bin is effectively unpopulated (n = 0 on day 286, n = 2 on day 284), so SEL validation says nothing below 5° peak elevation; the only evidence there is the pair transfer functions and the 2–5° frame rows.
- The frames were produced with NICEOPS adaf7c0 per `H/propagation/band_frames/manifest.txt`, whereas `H/README.md` attributes the finding to 84462d7; 84462d7 only added `--band_history`, so the propagation physics should be the same, but the manifest and README disagree on the commit.

### Gaps
- The per-band range-slope table (output of `band_frames_summary.py`) is not stored as text; only the README's summary ("63–630 Hz decays 1–2 dB/km more slowly", "−13 to −39 dB/km at 4–5 kHz") and the A-weighted table survive. The 2.5 GB frame set was not reprocessed for this note.
- No text file holds the final day-284 by-elevation table after the 2026-09-29 fixes (+0.79 / 2.99); it is quoted only in the two READMEs. `day284_B407_fixed.txt` is the previous stage (+0.26 / 3.58).
- No band-resolved residual by elevation exists for the SEL tables; band-level evidence is frames (day 286 approaches) and the pair transfer functions only.
- Other aircraft (AS350B3, B206L3, EC130B4, R44, R66) have measured tables (139 approach runs) but no propagation-option results were found in the files read.

---

## 5. Test-site ground and meteorological information available

### Takeaway
There is no measured ground impedance and no soil characterisation beyond "high-desert site"; the ground is known only through the fit to three co-located pole/board pairs. Meteorology is a ZephIR 300 LIDAR with lowest gate 46 ft, a tethered-balloon temperature string 0–500 ft every 10 ft, ground stations, and balloon sondes; nothing measures wind below 46 ft or turbulence.

### Cited Findings
- "The site is Amedee Army Airfield, CA, a high-desert site. Its ground impedance is unknown." — `P/docs/ground_plane_corrections.md` "The 2017 NASA hardware"
- Fitted site ground: variable porosity σe ≈ 200 (≈ 400 for Delany–Bazley), γT ≈ 3e-5, roughness 10–15 mm (8 mm in the first fit), from pairs 50/34 (y = +500 ft), 51/36 (y = 0), 52/38 (y = −500 ft), all at x = −3,760 ft. — `P/docs/ground_plane_corrections.md` "The 2017 NASA hardware", "Pinning the ground", "The ground profile redone with the real plate"
- Centerline pair differs: per-pair DB fit 120–180 kPa·s/m² at the sideline pairs, about 60 at the centerline; centerline residual stays worst (about 5.7 dB vs 3.6–4.2). "That is consistent with something other than ground impedance at the centerline." Pair 50/34 sits about 1.5 dB above the other two, "a calibration or gain difference". — `P/docs/ground_plane_corrections.md` "First look", "Fixed flow resistance", "First comparison"
- Pole height not recorded; fitted 4.08, 4.05, 4.08 ft (±0.005 bootstrap); report says "4 foot tripods". Poles are GRAS 46AE, diaphragm plane along the flight track (Part 36), foam windscreen of unrecorded size. — `P/docs/ground_plane_corrections.md` "Path 4", "The pole microphones' response", "Orientation, windscreen, roughness"
- Boards: GRAS 67AX in GR1425 (400 mm, 8 mm thick, 2.5 mm rim), laid on top of the ground, microphone 150 mm from center, outboard of the track; all ground microphones were flush 67AX including those labeled `invgb7`. — `P/docs/ground_plane_corrections.md` "The 2017 plates", "In the sphere builds"
- Ground boards sit 3–10 cm below the survey datum. — `N/CLAUDE.md`; `H/README.md`
- Array: lines across the track at X ≈ +608, 0, −1002, −1905, −2808, −3759, −5115, −6087 (ft), out to ±1500 laterally; about 52 microphones. — `M/noise-abatement-2017-data.md`
- LIDAR: ZephIR 300; Vh, Vv, dir at 46, 125, 197, 348, 495, 696, 984 ft AGL; 99999 = dead; logged 16–85 s after run start. Units taken as knots (confirmed two ways: day 284 LIDAR reads 1.05× the balloon in kt; R66 day 236 crab 7–9° vs knots-based prediction 8–15°, m/s would give 16–30°). Direction taken as "from", true; 180° ambiguity checked against aircraft drift (7/7 straight runs on day 284). — `M/noise-abatement-2017-data.md`; `H/README.md` "Conventions worth knowing"
- Temperature: `B407_Temperature_Data.csv`, per run, °F at 0–500 ft AGL in 10 ft steps (the balloon's temperature string). — `M/noise-abatement-2017-data.md`; `P/docs/ground_plane_corrections.md` "Pinning the ground"
- Ground stations: air temperature (°F), pressure (kPa), humidity (%); station vanes report a stale direction in calm. `Balloon_Sonde/` exists; its wind direction is off by a steady ~117° on day 284. — `M/noise-abatement-2017-data.md`
- B407 steady runs by day: 283, 284, 285, 286 → 29, 16, 48, 26 runs with LIDAR medians 3.1, 11.0, 3.2, 1.5 kt. Only one day-286 steady run has a LIDAR row. B407 approach (validation) runs are all day 286. — `M/noise-abatement-2017-data.md`
- Day 286: calm, 1–4 kt, inversion up to +7.6 K over 500 ft. — `M/propagation-roadmap.md`
- Day-284 mornings: upwind shear near the ground under an inversion, 330.5 m/s aloft against 339 at the ground. — `H/README.md`; `N/CLAUDE.md`
- Test report `Watts.TM2019.pdf` is in the data folder; both copies were unreadable on 2026-09-26 (OneDrive). — `M/noise-abatement-2017-data.md`

### Inferences
- Refinements that need a measured impedance spectrum, a layer depth, porosity, or a range-dependent ground cannot be fed from data on hand; only a re-fit of the three pairs (all at one x station) is possible, and the fit is known to trade flow resistivity against turbulence and against the pole's high-frequency response.
- Refinements that need the wind or temperature gradient in the lowest 14 m (the shadow-controlling layer) rest on the assumed log law (z0 = 0.1 ft) and on the 10 ft temperature spacing; refinements needing turbulence strength, stability (Obukhov length) or friction velocity have no direct measurement.

### Gaps
- Soil type, moisture, vegetation cover, surface roughness measurements: none found in any file read.
- Whether humidity varies with height (sondes exist but the harness uses the stations' value at every height).
- Whether the Watts report contains ground or micrometeorological detail beyond what is quoted; it could not be read.

---

## 6. Source-side limits at low elevation

### Takeaway
The 2017 spheres are built from steady runs with a 2,000 ft cap on the propagation path, a minimum-elevation floor (10° originally, 2° in the promoted rebuild), the plate BEM response divided out on the fitted site ground, and straight rays in a homogeneous atmosphere; directions below the floor are inverse-distance filled. Lowering the floor from 10° to 5° removed about half of the 5–10° over-prediction; 5° → 2° gained 0.07 dB.

### Cited Findings
- Builder defaults: `min_elevation_deg=10.0`, `max_propagation_range_ft=2000.0`, `board_correction='plate_bem'`. Rebuild settings: "a 2000 ft cap on the propagation path, point_stride 10, steady runs", `board_correction='plate_bem'`. — `P/noise_abatement_2017.py` lines 755–760; `H/spheres/rebuild.py` lines 6, 52–53
- The range cap's stated purpose (CLI help): beyond it the "atmosphere is a poor model (default 2000)". — `P/noise_abatement_2017.py` lines 1129–1132
- Depropagation assumes straight rays in a homogeneous atmosphere, "so test-day refraction is baked into the spheres." — `M/propagation-roadmap.md`
- Plate correction in depropagation: `depropagate_hemisphere(..., receiver_response_db=...)` divides each bin by the plate's band-averaged response per emission point, after the ambient gate and before spreading; response is `axisymmetric_bem.board_level` on the GR1425 over `SITE_GROUND` (VP σe 200), mic outboard, Q at the run's sound speed. — `P/docs/ground_plane_corrections.md` "In the sphere builds"
- The rim below the floor is inverse-distance (shepIDW) filled "flat from the data above", too loud; in-sample steady descents show forward mics +2.9 dB median, up to +11.7 at 7.5–9°; level flight at the same mics is exact; the excess is broadband 63 Hz–4 kHz. — `M/propagation-roadmap.md`; `H/README.md` "Low-elevation over-prediction"
- Rebuild comparisons (118 spheres, dBA, by elevation below the horizon):

| comparison | 0–5° | 5–10° | 10–20° | 20–40° | 40–90° |
| --- | --- | --- | --- | --- | --- |
| min10 control − sep25 original | −0.02 [−0.10, +0.04] max 0.32 | −0.02 [−0.08, +0.03] | −0.02 [−0.08, +0.04] | −0.02 | −0.02 |
| min5 − min10 | −0.95 [−4.51, +0.36] max 12.95 | −0.37 [−2.58, +0.40] max 13.23 | −0.07 [−1.07, +0.14] | −0.00 | +0.00 |
| min2 − min5 | +0.00 [−2.66, +0.81] max 9.90 | +0.00 [−0.55, +0.04] max 7.70 | +0.00 [−0.08, +0.00] | +0.00 | +0.00 |

(median [5–95%]). min5 − min10 at 0–10°, descent fore: median −1.01, 5th percentile −5.82; descent aft −1.09; level fore −0.15; climb 0.00. — `H/spheres/out/compare_sep25_min10.txt`, `compare_min10_min5.txt`, `compare_min5_min2.txt`
- Effect on day-286 SEL at 5–10° (station, standard measure): sep25 +2.83 / 3.54; min10 +2.75 / 3.49; min5 +1.51 / 2.11; min2 +1.43 / 2.05. Change from sep25 in that bin: min5 −1.57 [−3.1, −0.3], min2 −1.65 [−3.2, −0.4]. — `H/propagation/database_effect_B407_standard.txt`; run directories `H/propagation/out_min10`, `out_min5`, `out_min2`, `out_db_sep25` (each 25 runs; `out_min2/manifest.txt`: niceops f217339, panam becbf63)
- The 2° build is the promoted one (`databases/B407_{fixed_load_factor,ambient_gated}.nod`); "It is marginally ahead of 5° in that bin under all four models and level elsewhere, but its 2–5° data lean harder on the plate and ground models near grazing." — `H/README.md`
- The held-out database for day 284 is `spheres/out/min2_no284`. One day-284 run (284203) failed to build ("No valid hemisphere samples were generated"). — `H/README.md`; `H/spheres/out/min2.log`
- Spheres' reference: database's `EAA` built over 1,000 m; NICEOPS database sphere grid is 5° in direction; "At kilometer ranges [the launch angle] is a few degrees from the straight line, which is a large fraction of the database's 5° grid." — `N/CLAUDE.md`; `N/README.md` "Refraction"
- Other source-side facts: γ −24/+35 database entries are copies of the level sphere at each μ; the μ = 0 hover sphere is quieter (88 dB max vs 92–98) and differs from every flight sphere. — `M/propagation-roadmap.md`
- High-band contamination: top source bands are ambient floor amplified by depropagation; for Be407 from 3,150 Hz at the aft pole; onset moves aft (clean to 10 kHz at θ ≤ 120°, 6.3 kHz at 150°, 3,150 Hz at 165–180°). — `M/nod-database-generator.md`
- Speed basis: spheres take μ, γ and azimuth reference from ground velocity; the B407 spheres include 16 steady runs from day 284 flown in 9–15 kt. — `H/README.md` "Speed basis"

### Inferences
- With a 2° floor and a 2,000 ft cap, a sphere direction at 2–5° below the horizon is populated only by emission points within 2,000 ft slant range at source heights of roughly 70–175 ft (2,000 ft × sin 2–5°), i.e. by low passes close to microphones, and its level depends on the plate/ground model at grazing used in depropagation (VP 200, straight rays). Any error of the grazing ground model is therefore present twice, once in the sphere (divided out) and once in prediction (multiplied in); at matched geometry they cancel, at longer range than the build they do not. This arithmetic is mine, not from a file.
- The remaining +1.1 to +1.4 dB at 5–10° after the 2° build coincides with the forward-line short-range excess (section 4b), which the files attribute to the end of the approach (hover sphere and hover-to-slowest-flight blend), not to the rim.

### Gaps
- No experiment with a longer range cap exists (the user rejected it); no sphere build with refracted depropagation exists (roadmap intent only).
- No file quantifies how many samples populate the 2–5° directions in the min2 build.

---

## 7. Already-identified candidate causes of the grazing-incidence residuals, with evidence

### Takeaway
The files name six candidates. Two are established (the sphere rim; the end-of-approach source coverage). The rest, for the 0–2° flush residual and the long-range low-frequency deficit, are open, with refraction of the direct and reflected rays near grazing and a ground that differs from the fitted one (near the receiver or along the path) the two explicitly named alternatives.

### Cited Findings

| # | candidate | which residual | evidence for | evidence against / limits |
| --- | --- | --- | --- | --- |
| 1 | Sphere rim filled below `min_elevation_deg = 10` | +2.7 dB at 5–10° (day 286) | 5° floor: +2.42/3.23 → +1.21/1.89; in-sample descents +2.9 median, up to +11.7 forward; flat rim visible in every B407 sphere | 2° floor gains only 0.07 dB; "+1.1 dB is not the filled rim" |
| 2 | End-of-approach source coverage (hover sphere, hover-to-slowest-flight blend, forward arc) | forward-line +1.5 to +1.65 dB | excess at short range, height < 400 ft, 0–20 kt (+3.19); other lines right at same height and speed; θ 30–60 forward +1.78 | earlier test: flooring speed at 37 kt "changes nothing" for the 5–10° bin; state matched by height, exact only for monotonic descent |
| 3 | Refraction near grazing, not modeled in the pole/board reference or mis-modeled in the run | 0–2° flush residual (5–10 dB above models from 500 Hz); −2 dB beyond 4,000 ft at 2–10°; 63–630 Hz decaying 1–2 dB/km too fast in the model | 0–2° "wants a far harder ground (σ 3200–6400)", which "is how a refraction-free model would try to imitate" differential bending at a 1.2 m pole and a flush board; day 286 had an inversion (+7.6 K over 500 ft), "first suspect is the calm-day inversion refracting sound downward"; downward-refraction multi-bounce arrivals are not in the model | on day 286 `linear` and `stratified` do not improve the SEL (rms 1.42–1.46 vs 1.30–1.39) and make 5–10° slightly worse; on windy 284153 (no shadow) the measurement falls 4–6 dB/km *faster* than every model; near-ground wind below 46 ft is assumed |
| 4 | Ground differs from the fitted VP σe 200: harder near the receiver, range-dependent, or wrong model class | 0–2° flush residual | DB 3200 gives 1.89 dB at 0–2°; VP 1600 2.75; "A ground that changes along the path is the other candidate"; methods survey: "the 0–2° flush residual may be the ground near the receiver, since the level goes as about −40 log|β|"; centerline pair fits a different σ (≈ 60 vs 120–180); steep-angle fits give σ 1300–5300 | a hard ground conflicts with every other elevation (DB 3200: 3.32 at 2–5°, 2.24 at 5–10°); no measured impedance to decide; Delany–Bazley known poor for grassland but the site is desert |
| 5 | The pole reference itself (pole ground reflection, turbulence, 46AE response, roughness) | 0–2° pair residual | the measured transfer at grazing "rests on the pole's modeled ground reflection (−10 to −14 dB at 125–500 Hz below 2°)"; refitting the pole model raised the measured near-grazing level by up to 8 dB; γT 3e-5 is 10× what the SEL supports, "probably carries more than turbulence (roughness, the pole's height)" | after refit the 0–2° error is still 6.3–6.8 dB with the pinned ground; applies to the pair comparison only, not to the SEL/frame residuals at flush boards |
| 6 | Turbulent scattering / decorrelation | day-284 shadow depth; possibly grazing nulls | turbulent PE (γT 3e-6) gets shadow depth −17 to −18 vs −12 to −17 measured where still air is below −35 | SEL consistent only up to ~3e-6; on day 286 scattering moves overall rms 1.45 → 1.46; plate path has no coherence loss by construction (γ_p = 0) |

Sources: row 1 `H/README.md`, `H/propagation/database_effect_B407.txt`, `M/propagation-roadmap.md`; row 2 `H/README.md`, `H/propagation/band_frames/aweighted.txt`, `M/forward-line-error-source.md`; row 3 `P/docs/ground_plane_corrections.md` "Pinning the ground", `M/forward-line-error-source.md`, `M/propagation-roadmap.md`, `H/propagation/day286_B407_plate.txt`, `H` git log 085b6c1, `N/README.md` "Refraction"; row 4 `P/docs/ground_plane_corrections.md`, `M/propagation-roadmap.md` (methods survey, 2026-09-30); row 5 `P/docs/ground_plane_corrections.md` "Literature review" and "Path 4", `N/README.md` "Turbulence"; row 6 `H` git log 085b6c1, `N/README.md`, `H/README.md`.

Further items the files list as ruled out or second order for grazing: plate thickness (±1 dB), plate model family (Fresnel vs BEM within 0.1–0.2 dB up to 20°), microphone position on the plate, the shadow floor value, QSAM, frame mismatch, sphere range cap, Doppler shift. — `P/docs/ground_plane_corrections.md`; `H/README.md`; `M/propagation-roadmap.md`; `M/sphere-doppler-convention.md`

Directions the files themselves record (not recommendations of this note): the 2026-09-30 methods survey's "leading conclusions" were reciprocal per-run PE excess-attenuation tables for low bands and low elevations (Collins split-step Padé for a 5–10× larger range step; precedent Cotté 2019), effective-sound-speed error for elevated sources (Colas 2023), NWP profiles with MOST fill-in, and stability-dependent turbulence; "The user wanted more research before writing the document." — `M/propagation-roadmap.md`

### Inferences
- Candidates 3 and 4 are not separable with the data on hand: both predict a flush microphone reading higher than the soft-ground model at 0–2°, and the site has neither a measured impedance nor wind/temperature gradients in the lowest 14 m. The files say as much ("refraction or a range-dependent ground").
- Implementation features that bear on the low-frequency, long-range deficit and that no file has tested against it: Q evaluated at a single band-center frequency in the point model; the |w| > 500 cut to the plane-wave coefficient; local reaction for all impedance models; αe fixed at 0; a single reflected arrival under downward refraction; ground term taken through an equivalent straight geometry. These are properties of the code (section 1), listed here as untested, not as identified causes.
- The long-range deficit was measured under `stratified` + plate + scattering 3e-6; no frame-level comparison under straight rays or `linear` was found, so how much of the −2.4 dB/km trend is introduced or removed by the refraction model on the calm day is not known from these files.

### Gaps
- No file tests an alternative impedance model (slit-pore, Wilson relaxation, two-parameter, extended reaction) or a non-zero αe against the 0–2° or long-range data.
- No file tests multi-bounce downward-refraction arrivals or a PE at the day-286 long-range geometries (PE cases stored are log laws, a linear lapse, `example_windy`, and 18 day-284 upwind cases — `H/propagation/calibration/README.md`, `H/propagation/pe_runs/cases.txt`).
- No file separates the long-range deficit by microphone line beyond forward/other, or by upwind/downwind on day 286.
