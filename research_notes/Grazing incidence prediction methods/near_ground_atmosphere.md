# The lowest ~15 m of the surface layer in grazing-incidence propagation: profile models, sensitivity, and methods for what rays miss

Scope: elevated sources (helicopter, 10–700 m), receivers 0–1.2 m (flush boards and 1.2 m poles), 0.3–4 km, 50 Hz–5 kHz, 0–10° elevation. This builds on `research_notes/Grazing incidence soft ground propagation/` (refraction_turbulence_ground_theory.md §6, numerical_methods_grazing.md §3 and §6). It does not repeat what those notes established: lin-log profiles from MOST, MOST's breakdown in strong inversions, the 5–15 dB cost of linearising (van Maercke & Salomons), the Harmonoise classes, the 2 dB versus 8–10 dB mean/event rms (Wilson et al. 2015), Waxler's quiet height, HPM, and the 5–7 dB of the Hart et al. surrogates.

Source tags: **[FULL]** full text read; **[ABS]** abstract or publisher summary only; **[SN]** search-engine snippet only, so weakest; **[OWN CALC]** my own derivation or arithmetic; **[RECALLED]** a standard result from memory that I did not re-verify in a source this session.

---

## 1. Near-surface profile models beyond MOST for stable, very stable, calm and morning-transition conditions, and for a warm surface layer under an inversion aloft. Can their parameters be estimated from the instruments available?

### Takeaway
There is no agreed flux–profile law for the very stable surface layer. The stability functions in use split three ways once ζ = z/L exceeds about 1. Extrapolating Businger–Dyer gives φ = 1 + 5ζ, a *linear* profile. Levelling functions (Cheng–Brutsaert, local scaling) give φ → about 7, a *steep log* profile. SHEBA gives φ_m ∝ ζ^(1/3), an intermediate power law. Each puts a different gradient in the lowest 15 m, which is exactly the layer below the lidar. The practical route from the available data has three steps. First, a bulk Richardson number from 2 m temperature, the lidar's first gate and the sonde, inverted non-iteratively (Gryanik et al. 2020). Second, once Ri_b passes about 0.2–0.25, give up on similarity and use a jet-shaped profile (Banta: shear below the jet nose near-linear, clustering about 0.1 s⁻¹ when turbulent). Third, for the morning transition and the warm-ground-under-inversion case, use a two-layer construction: a shallow unstable surface layer topped by the remnant inversion. That construction is not a published acoustic recipe, and all three choices should be carried as an ensemble, not picked.

### Cited Findings
**Stable and very stable stability functions**
- CASES-99 (60 m tower, Kansas, Oct 1999): "The gradient stability functions were found to be linear, namely φm = 1 + 5.8ζ and φh = 1 + 5.4ζ up to a limit of … ζ = 0.8". Beyond a transition, "both functions were observed gradually to approach a constant, with a value of approximately 7". A two-parameter form covers neutral to very stable. **[ABS]** — [Cheng & Brutsaert 2005, BLM](https://consensus.app/papers/details/0a60909858e05356a253dd61da80e33a/?utm_source=claude_desktop)
- SHEBA (20 m tower, five levels, 11 months over Arctic pack ice): "In the limit of very strong stability, φm follows a ζ^1/3 dependence, whereas φh initially increases with increasing ζ, reaches a maximum at ζ ≈ 10, and then tends to level off". New φm and φh cover "the entire range … from neutral to very stable". **[ABS]** — [Grachev et al. 2007, BLM](https://consensus.app/papers/details/e721d71e25a75635aa364a43c3f6b1df/?utm_source=claude_desktop)
- Local scaling (SABLES98; sonics at 5.8, 13.5 and 32 m; 14 thermocouples): φm and φh increase with ζ "up to a value of ζ≈1–2, above which the values remain constant". As a result "the turbulent mixing is underestimated when linear similarity functions (Businger et al., 1971) are used". **[ABS]** — [Yagüe et al. 2006, NPG](https://consensus.app/papers/details/0f7651c4ddac552a8bacc8fde5ba25b9/?utm_source=claude_desktop)
- Regimes in terms of local z/Λ*: nearly neutral (< 0.02), weakly stable (0.02–0.6), very stable (0.6–50), extremely stable (> 50). In the very stable regime the dimensionless gradients go as (z/Λ*)^(3/5). "The existence of scaling laws in 'extremely-stable' conditions is doubtful." Pr falls from 0.9 to about 0.7, and steady turbulence needs "Ri < 0.7". **[ABS]** — [Sorbjan & Grachev 2010, BLM](https://consensus.app/papers/details/634ae590dc1d54dd9484d14e7257bbc9/?utm_source=claude_desktop). In Ri terms: nearly neutral < 0.02, weakly stable 0.02–0.12, very stable 0.12–0.7, extremely stable > 0.7, where turbulence "must be non-stationary, i.e. decaying or sporadic". **[ABS]** — [Sorbjan 2010, QJRMS](https://consensus.app/papers/details/12605594f65f5caf9acde6b0a229b648/?utm_source=claude_desktop)
- A 2021 review of the four main formulations: "None of them turned out to be fully consistent with the literature", and "a convincing description of flux–profile relationships still needs to be found". **[ABS]** — [Casasanta et al. 2021, Atmosphere](https://consensus.app/papers/details/ce9c136c44bc5e33bfac2677bc0a4f34/?utm_source=claude_desktop)
- Functions derived from the MYNN PBL scheme: φm = 1 + 4.8ζ, φh = 0.74 + 6.0ζ, giving critical bulk and flux Ri of 0.26 and 0.21 as ζ → ∞. **[ABS]** — [Nakanishi et al. 2022, JMSJ](https://consensus.app/papers/details/81f166aa27b756619d83b949940bdfb1/?utm_source=claude_desktop)
- A non-iterative parametrization: "transfer coefficients result as a function of the bulk Richardson number (Rib) and roughness parameters". It agrees with SHEBA for 0 < Ri_b < 0.5 and with Businger–Dyer for Ri_b < 0.08 (0.1 in the journal version). **[ABS]** — [Gryanik, Lüpkes, Grachev & Sidorenko 2020, JAS](https://consensus.app/papers/details/48450b6196d45082bf2c08b956210d18/?utm_source=claude_desktop); [preprint](https://consensus.app/papers/details/63fbf1f6e34359d8b9a764291a9d69d5/?utm_source=claude_desktop)
- Beljaars & Holtslag (1991, J. Appl. Meteorol. 30, 327–341) give stable-case ψ functions with an exponential term in exp(−0.35ζ). One quoted variant is "− 0.7ζ + 0.75(ζ − 14.3 exp(−0.35ζ) + 10.7)". **[SN]** — [search summary of arXiv:0911.2286 and others](https://arxiv.org/pdf/0911.2286). (The exact constants should be checked against the paper before coding.)

**Low-level jets and the sub-jet layer**
- CASES-99 HRDL lidar: "The vertical profile of the mean along-wind component U(z) was often nearly linear below this jet maximum (or 'nose')". Its departure from linearity comes "within a few meters of the surface, where ∂U/∂z (and also ∂θ/∂z) becomes very large". For TKE > ~0.1 m² s⁻², sub-jet shear U_X/Z_X "cluster[s] around a shear value of 0.1 s⁻¹", which suggests the nose adjusts to keep that shear. TKE was small for jet Ri_J > ~0.4 and rose below it. θ(z) was "often roughly linear in the SBL between the surface layer and the nocturnal inversion top". **[FULL]** — [Banta, Pichugina & Newsom 2003, JAS 60, 2549](https://meteor.geol.iastate.edu/windresearch/resources/llj.banta.tke.and.llj.maximum.2003.pdf)
- A five-parameter "log–jet" fit, U(z) = (u*/κ) ln(z/z0) + U_jet(z; U_m, z_m, S), was fitted to WindCube 100S lidar profiles from 80 m up (fits with R² ≥ 0.90 kept) and used to compare lidar with ERA5 and NORA3. **[FULL abstract + method summary]** — [Bui, Bakhoday-Paskyabi & Reuder 2026, WES 11, 2307](https://wes.copernicus.org/articles/11/2307/2026/)

**Morning transition**
- Cabauw and Flatland observations: "the morning transition is largely driven by entrainment rather than by the surface heat flux. … the heat that warms the layer between the surface and 200 m, eroding the surface inversion, comes almost entirely from above". "The small positive surface heat flux appears to act primarily to relax the stability." Stronger winds shorten the time from flux crossover to onset. **[FULL extended abstract]** — [Angevine, Klein Baltink & Bosveld (AMS preprint of BLM 101, 209, 2001)](https://ams.confex.com/ams/pdfpapers/43198.pdf)
- Superadiabatic layers are "generally shallow and located near the earth's surface", with "temperature gradients of more than 50 K/km" close to the surface, commonly under clear skies, light winds and dry soils. **[SN, low-quality aggregator]** — [weather-prediction aggregator page](https://theweatherprediction.com/habyhints/31/)
- In a conventionally neutral boundary layer (neutral surface layer under a capping inversion) the log law underestimates the wind. PE runs (CNPE, near-ground sources) show that this leads to a phase error that grows with distance and to under-prediction of "the energy refracted downward". **[ABS]** — [Libianchi et al. 2023, JASA](https://consensus.app/papers/details/876f5b20e90c5f1682c0043ca24a52da/?utm_source=claude_desktop)

**Ray-model equivalent profiles**
- Equivalent near-ground wind profiles for an axisymmetric eigenray tracer were set "by correlating predictions with measurements" in downward and upward refraction. The paper gives guidelines for the wind profile, but the abstract carries no numbers. **[SN/ABS]** — [Prospathopoulos & Voutsinas 2007, JASA 122, 1391](https://www.openarchives.gr/aggregator-openarchives/edm/ntua/000011-123456789_17924)

### Inferences
- **[OWN CALC] What each very-stable law does to the shape of the lowest 15 m.** With dU/dz = (u*/κz) φ(ζ):
  - Linear φ = 1 + βζ gives U = (u*/κ)[ln(z/z0) + βz/L]. Once z ≫ L/β this is *linear*, with gradient βu*/(κL), uniform in height.
  - Levelling φ → φ∞ ≈ 7 gives dU/dz = 7u*/(κz). That is still a *log* profile, seven times steeper, so most of the shear (and, for θ, most of the temperature step) sits in the lowest meters.
  - SHEBA's φm ∝ ζ^(1/3) gives dU/dz ∝ z^(−2/3), in between.
  
  For a given lidar wind at 14 m (46 ft) and a given 2 m temperature, the three laws spread the same total Δc_eff over 0–14 m very differently. That spread is what reaches the acoustics; see §2 and §3.
- **[RECALLED]** The other forms needed:
  - Businger–Dyer (stable): φm = φh = 1 + 5ζ, ψ = −5ζ; Businger's original 1 + 4.7ζ with Pr0 = 0.74.
  - Cheng–Brutsaert: φ = 1 + a[ζ + ζ^b(1+ζ^b)^((1−b)/b)]/[ζ + (1+ζ^b)^(1/b)], with (a, b) = (6.1, 2.5) for momentum and (5.3, 1.1) for heat.
  - Grachev 2007: φm = 1 + a_m ζ (1+ζ)^(1/3)/(1 + b_m ζ) with a_m = 6.5, b_m = 1.3; φh = 1 + (a_h ζ + b_h ζ²)/(1 + c_h ζ + ζ²) with a_h = b_h = 5, c_h = 3.
  - Unstable and morning surface layer: free-convection scaling, in which dθ/dz ∝ z^(−4/3).
  
  Check every coefficient against the papers before coding.
- **Estimating parameters from what the 2017 tests had** (inference; not a published acoustic method):
  1. Bulk Ri_b between 2 m (surface station T, with wind from the 10 m mast) and the lidar's 14 m gate (wind) plus the sonde (T). Invert with Gryanik et al.'s non-iterative fit to get ζ and hence u* and θ*. z0 = 0.1 ft is already assumed. The scalar roughness z0h is not observable from these instruments and over desert can differ from z0 by orders of magnitude (gap).
  2. If Ri_b < ~0.1, the Businger–Dyer, Cheng–Brutsaert and Grachev profiles agree closely, and a single profile is defensible.
  3. If 0.1 < Ri_b < ~0.25, carry all three as an ensemble.
  4. If Ri_b exceeds ~0.25 (above the MYNN critical value and inside Sorbjan's "very stable"), drop similarity. Instead, interpolate θ between the 2 m value and the sonde, and U between a small surface value and the lidar. Use a Banta-style linear sub-jet profile with a log foot. Carry the depth of the log foot as an uncertain parameter.
- **Warm ground layer under an inversion aloft** (day 286, run 286288). This is a morning-transition structure: a shallow layer heated from the surface, superadiabatic in its lowest meter or so, under a remnant nocturnal inversion that Angevine et al. say is eroded mainly from above. MOST in its unstable form can describe the shallow layer only if a flux exists. With calm wind it is close to free convection, where the gradient goes as z^(−4/3), concentrated at the bottom. A defensible profile is therefore two-layer:
  - (i) A surface lapse layer of depth d (0.3–3 m) and temperature excess ΔT, from the surface-station 2 m T, the sonde's lowest point, and if possible a skin-temperature estimate.
  - (ii) The measured inversion above it.
  
  d and ΔT are the two parameters that matter (§3). Neither is constrained by the lidar.
- The Bui et al. log–jet form is an off-the-shelf five-parameter fit that matches what the brief has: lidar from 14 m up, with z0 fixed. It extrapolates downward through its log term, so it inherits the neutral-log assumption below the first gate and does not fix the very-stable shape problem above.

### Gaps
- No source found gives measured temperature profiles in the lowest 0–2 m over desert or airfield ground just after sunrise under a remnant inversion (ΔT, d). The magnitudes in §3 are therefore parametric. The only figure found ("more than 50 K/km") is an aggregator snippet.
- Beljaars–Holtslag constants: snippet only. The Cheng–Brutsaert and Grachev coefficients are recalled, not re-read.
- I found no published method that blends a lidar first gate, a surface station and a sonde into an acoustic profile for the lowest 15 m. The procedure above is a suggestion.
- Scalar roughness z0h over desert, which controls θ* from a skin–air difference, was not researched.

---

## 2. How sensitive is propagation at 0–10° elevation (30–300 m sources, 1–4 km) to these choices? What do the acoustic studies say, and how are profile ensembles and UQ used?

### Takeaway
Whether a straight-line arrival at a flush microphone survives the bottom layer is set by the layer's total effective-sound-speed step, not its depth. The critical arrival angle is θc ≈ √(2Δc_eff/c) [OWN CALC]: 1.8° for 0.3 K, 3.3° for 1 K and 5.8° for 3 K of warm-ground excess. That covers every 30–100 m source beyond about 1 km and 300 m sources beyond about 3–5 km. How many decibels this costs depends on the layer's thickness relative to the diffraction length (§3), so at low and mid frequencies the lowest meter matters much less than geometric rays imply. The acoustic literature agrees that profile *shape* (log versus linear, CNBL versus log) changes elevated-source levels by several dB. It also agrees that inherent prediction scatter even with excellent met is 6–8 dB for events, and that ensembles should be propagated by LHS or polynomial chaos with a handful to a few hundred samples.

### Cited Findings
- Elevated sources upwind: nonlinear (e.g. logarithmic) c_eff profiles produce "local SPL increases of up to 7 dB", a caustic that "typically emerges at intermediate distances—potentially extending to a few kilometers—before transitioning into the shadow zone". Source height and the vertical gradients of c_eff are the key factors, and estimation formulas are given. **[ABS]** — [Gharbi et al. 2026, Acta Acustica](https://consensus.app/papers/details/89f6fa0ca10555559ddad6abe66d18f4/?utm_source=claude_desktop). The companion paper says the increase needs "an elevated source location, a non-linear (e.g. logarithmic) wind or temperature profile", and uses "two fundamentally different modeling approaches". **[ABS]** — [Baumann, Gharbi, Juliust & Schady 2025/2026, Acta Acustica](https://consensus.app/papers/details/b976cee37b6056ed83c70dd054eb71b8/?utm_source=claude_desktop)
- eVTOL footprint with a ray model and FW-H hemisphere sources, the same architecture as NICEOPS: weather "can contribute to mismatch up to 4 dBA in the illuminated zone and a significant drop in the refractive shadow zone". **[ABS]** — [Fuerkaiti et al. 2021, JASA](https://consensus.app/papers/details/b05b0c351cd55240bdf8f87bd18ec91f/?utm_source=claude_desktop)
- ASCENT Project 40 extended Wilson et al. (2014) UQ "to the geometry of aircraft noise propagation", using a simplified Nord2000 for speed, to give the SPL effect of individual met uncertainties (temperature profile, wind profile) and of their combination at terminal-area altitudes. **[ABS, no numbers]** — [Patankar et al. 2019, JASA (abstract)](https://consensus.app/papers/details/0695aca775da5a87abc4eb5735f03aa6/?utm_source=claude_desktop)
- Wilson et al. 2008 (NCEJ 56, 465): actual mean profiles and MOST profiles give similar predictive errors. A windless neutral atmosphere gives "large positive bias errors in upward refraction". Skill improves with bandwidth, especially downwind. **[ABS]** — [Wilson, Lewis, Weatherly & Andreas 2008](https://ince.publisher.ingentaconnect.com/contentone/ince/ncej/2008/00000056/00000006/art00007;jsessionid=19wl5h4dupv75.x-ic-live-01)
- In Wilson's 2024 retrospective, even "with accurate meteorological measurements from a location and time very close to the propagation path, sound level predictions have inherent random errors of about 6-8 dB". Errors "diminish when mean meteorological measurements (over an interval of, say, a half hour) are used to predict mean sound levels". The finding came from LES ground truth combined with acoustic modeling. **[FULL]** — [Wilson 2024, Acoustics Today 20(4), 78](https://acousticstoday.org/wp-content/uploads/2024/12/sound-predicting-wilson.pdf)
- UQ methods. Wilson et al. 2014 (JASA 136, 1013) is the reference framework. In a later study (2016), importance sampling "converges most rapidly when uncertainties in atmospheric and ground properties dominate", compared with MC, LHS, and stratification by met class. The atmosphere there is from surface-layer similarity and the field from PE. **[SN]** — [sonicfield library record](https://sonicfield.org/library/uncertainty-quantification-turbulent-scattering-and-random-sampling-in-outdoor-s)
- Lloyd's-mirror UQ demonstration:
  - LHS reaches MC's 10,000-sample accuracy for the mean TL "using only 100 computations", and generalized polynomial chaos (GPC) "may achieve the required accuracy with only 5 samples".
  - GPC outperforms MC and LHS for small and medium variability, but "MCS and LHS perform better for large variability".
  - A 1% sound-speed error can move near-ground TL at 50 m "by as much as 15 dB".
  - The GPC coefficients give a cheap surrogate; the authors suggest substituting such a surrogate inside "a model of aircraft noise near the ground".
  
  **[FULL]** — [Martinelli, Wilson, Wixom & Pettit 2023, Acoustics Today 19(2), 28](https://acousticstoday.org/wp-content/uploads/2023/06/Uncertainty-in-Acoustical-Modeling-Sheri-L.-Martinelli-D.-Keith-Wilson-Andrew-S.-Wixom-and-Chris-L.-Pettit-1.pdf)
- CNBL versus log-law wind (near-ground sources, CNPE): the underestimated wind gives a phase error growing with distance and an underestimate of downward-refracted energy. **[ABS]** — [Libianchi et al. 2023](https://consensus.app/papers/details/876f5b20e90c5f1682c0043ca24a52da/?utm_source=claude_desktop)
- Effective sound speed approximation (ESSA) in WAPE: "large phase error (several dB) on monochromatic calculations with a point source", but "small phase error (<1 dB)" for an extended broadband wind-turbine source. **[FULL]** — [Kayser et al. 2023, JASA 153, 1846](https://perso.ensta-paris.fr/~cotte/files/Kayser_al_jasa2023.pdf)
- JAXA's GF-PE aircraft model was tested against balloon-suspended elevated sources over a runway and grass (Taiki, 2012–2013), with air-to-ground impulse responses under various met conditions. **[ABS, no numbers]** — [Yokota et al. 2014](https://consensus.app/papers/details/6525eb26644f55d0aec1c4815df1b818/?utm_source=claude_desktop)
- A measured analogue at a desert site: F-18 sonic-boom lateral-cutoff tests at Edwards AFB found "considerable noise 1 – 2 nm (1.9 – 3.7 km) beyond numerical predictions" (PCBoom, geometric acoustics). During "strong temperature inversions" near the ground, levels were "80 dB at 6.6 nm (12 km)" where "<60 dB [was] expected at 8 nm (15 km)". **[FULL slides]** — [Cliatt et al. 2015, NASA FaINT, ISNA/Sonic Boom Forum](https://ntrs.nasa.gov/api/citations/20160001338/downloads/20160001338.pdf)

### Inferences
- **[OWN CALC] Geometric sensitivity.** Snell's law across a near-ground layer with effective-sound-speed excess Δc at the ground (warm ground, or upwind shear) gives θ_ground² ≈ θ_top² − 2Δc/c. Rays arriving at the layer top below θc = √(2Δc/c) turn before reaching z = 0, *whatever the layer depth*.

  | near-ground excess | Δc (m/s) | θc |
  |---|---|---|
  | ΔT = 0.3 K | 0.18 | 1.8° |
  | ΔT = 1 K | 0.59 | 3.3° |
  | ΔT = 3 K | 1.8 | 5.8° |
  | upwind Δu cos φ = 0.5 m/s | 0.5 | 3.1° |
  | upwind Δu cos φ = 1.3 m/s | 1.3 | 5.0° |

  (c ≈ 343 m/s, T ≈ 293 K, dc/dT ≈ 0.59 m/s/K.)

  Straight-line elevations, for comparison:

  | source height | 1 km | 2 km | 4 km |
  |---|---|---|---|
  | 30 m | 1.7° | 0.9° | 0.4° |
  | 100 m | 5.7° | 2.9° | 1.4° |
  | 300 m | 16.7° | 8.5° | 4.3° |

  So ~1 K of warm ground air blocks geometric arrivals for 30 m sources at all ranges ≥ 0.5 km, for 100 m sources beyond ~1.7 km, and for 300 m sources beyond ~5 km. Downwind of the same layer, or under an inversion, the step has the opposite sign and only flattens the arrival.
- **[OWN CALC] What the assumed log law contributes below the first gate.** With u(14 m) = 3 m/s and z0 = 0.1 ft (0.03 m): u* ≈ 0.20 m/s, u(1 m) ≈ 1.7 m/s, so Δu(1→14 m) ≈ 1.3 m/s. Up- or downwind that equals a ~2 K temperature step (θc ≈ 5°). A wrong z0, or a very-stable law that piles the shear lower (§1), changes this by tens of percent. That is enough to move θc across the 1–3° arrival angles of the far microphones. Under calm wind the shear term vanishes and temperature dominates.
- The 1.2 m pole and the flush board respond differently. A warm layer thinner than 1.2 m leaves the pole microphone above the turning region, so only the board sits in the micro-shadow. That is consistent with run 286288 being a flush-microphone failure. The lowest meter is essentially invisible to a 1.2 m receiver except through the ground-reflected path.
- Given the 6–8 dB inherent event scatter (Wilson 2024) and the large-variability regime where MC/LHS beat GPC (Martinelli 2023), a near-ground ensemble for NICEOPS should be small and low-dimensional. Candidate dimensions: (ΔT or Δc of the surface layer, its depth d, the stable-law choice or ζ, z0). Propagate it with LHS (tens of samples) or with tensor GPC on 2–3 parameters (5–25 samples). The NICEOPS ray models cost milliseconds per pair, so any of these is affordable.
- Gharbi et al.'s upwind caustic (+7 dB, log versus linear) is the opposite-signed sibling of the NICEOPS log-law shadow-boundary finding (CLAUDE.md, stratified): profile *curvature* near the ground relocates caustics and shadow boundaries for elevated sources by kilometers. Both say that linearising the bottom 15 m is not safe for a helicopter footprint.

### Gaps
- No source found quantifies, for elevated sources, the dB sensitivity to the lowest 1–10 m separately from the rest of the profile. The angle table above is arithmetic; the dB consequence is in §3.
- Gharbi et al. and Baumann et al.: abstracts only (Acta Acustica full text returned 403). Their estimation formulas and source heights were not read.
- Patankar et al. (ASCENT 40) and Wilson et al. 2014: no numbers retrieved.
- Wilson 2003 (JASA 113, 750; transitional height) was not re-read. It is cited in the earlier notes.

---

## 3. A warm surface layer under an inversion aloft: the flush-receiver micro-shadow. How does the literature treat it, and what recovers it?

### Takeaway
I found no acoustic paper on this exact configuration, a sub-meter warm layer under an inversion aloft with a flush receiver. Two bodies of work bear on it: tunneling and leakage through evanescent layers (ducts, ocean acoustics), and shadow-zone diffraction and scattering (Cheinet; Gilbert et al.). Both say the ray-theory silence is unphysical and that low frequencies leak through thin barriers. [OWN CALC] The leak has a closed form. For a linear warm layer of depth d, the grazing field at the ground is reduced by about exp(−(2/3)(d/l)^(3/2)), where l = (R_c/2k²)^(1/3) is the diffraction length at the layer's curvature R_c = c/|dc/dz|. That is about 1 dB at 100 Hz and 11 dB at 1 kHz for a 3 K, 1 m layer. Compare the ~21 dB deficit the layered ray model shows at 8,000–13,000 ft. The cheap fix is a uniform (Airy) turning-point continuation of the dipping-ray family below its turning height. PE remains the arbiter, and turbulence scattering adds a floor above the tunneled field.

### Cited Findings
- "Ray theory gives no useful information about sound levels below acoustic surface ducts, predicting instead a region of silence (shadow zone), whose lack of physical reality has been proved repeatedly in practice." **[SN]** — [Roebuck, "An heuristic study of the propagation of sound to long ranges below and within acoustic surface ducts", Proc. IOA](https://www.ioa.org.uk/system/files/publications/I%20ROEBUCK%20AN%20HEURISTIC%20STUDY%20OF%20THE%20PROPAGATION%20OF%20SOUND%20TO%20LONG%20RANGES%20BELOW%20AND%20WITHIN%20ACOUSTIC%20SURFACE%20DUCTS.pdf)
- Tunneling "is a leakage mechanism which diminishes with increasing frequency and with increasing distance of the source/receiver from the duct boundary". Energy leaks "through paths typically neglected by ray models", and at low frequencies it "can be surprisingly strong so that a receiver in a ray-theory shadow zone is actually well insonified". **[SN]** — [NJIT, "Anomalous parabolic equation results for propagation in leaky surface ducts"](https://digitalcommons.njit.edu/fac_pubs/17145)
- Tunneling rays from an evanescent region, and the smooth transition from total reflection to transmission, "can be described in a ray tracing approximation by extending the treatment to complex rays". **[SN]** — ["Acoustic tunneling", Proc. Roy. Soc. Edinburgh A](https://resolve.cambridge.org/core/journals/proceedings-of-the-royal-society-of-edinburgh-section-a-mathematics/article/abs/acoustic-tunneling/B24D266928714E6A7B0D57DF8BF288C2)
- Near-surface refractive shadows, with stratification and turbulence derived together from standard met through flux–profile relations and a second-order-moment PE: shadow levels "decrease with adverse wind as a consequence of enhanced mean upward refraction, and increase with the absolute value of the surface heat flux as a consequence of enhanced turbulence scattering". **[ABS]** — [Cheinet 2012, JASA (search summary)](https://data-ww3.ifremer.fr/BIB/Cheinet_JASA2012.pdf) (PDF returned 403)
- "Two mechanisms contribute to sound in the shadow zone—diffraction and turbulent scattering—with scattering due to turbulence being the predominant mechanism in the short wavelength limit." **[SN]** — [search summary of Gilbert, Raspet & Di / NASA 19910007390](https://ntrs.nasa.gov/citations/19910007390)
- LEE (FDTD) has been evaluated for strong and saturated turbulence on long 2-D harmonic paths, the regime where PE's one-way and narrow-scatter assumptions are stressed. **[SN]** — [Ehrhardt, Cheinet et al. 2013, sonicfield record](https://sonicfield.org/library/evaluating-a-linearized-euler-equations-model-for-strong-turbulence-effects-on-s)
- Gaussian beams improve "the accuracy in caustics and shadow zones" over ray fields. Mo et al. reproduce the Attenborough benchmark (h_s = 5 m, h_r = 1 m, 10/100/1000 Hz, downward, upward and mixed) and note that "previous work based on Gaussian beams was reported to have difficulties with upward refractive conditions". Agreement is shown graphically; no rms is given. **[FULL]** — [Mo, Yeh, Lin & Manocha 2017, JASA 141, 2289](https://gamma.cs.unc.edu/PAPERS/OutdoorSoundPropagation17.pdf), citing Gabillet et al. 1993 (JASA 93, 3105)
- The desert sonic-boom data in §2 (Cliatt et al. 2015) are a measured case of large levels well beyond the geometric-acoustics limit, with a near-ground inversion present. **[FULL slides]** — [NASA FaINT](https://ntrs.nasa.gov/api/citations/20160001338/downloads/20160001338.pdf)

### Inferences
- **[OWN CALC] Tunneling loss through a linear warm layer.**
  - Take a layer of depth d with c decreasing linearly by Δc from the ground to its top. A ray reaching the top at grazing angle θ has horizontal wavenumber k_x = k(d) cos θ. Below the turning height the field is evanescent, with vertical decay rate q(z) = √(k_x² − k(z)²) ≈ k√(2Δc/c·(1 − z/d) − θ²).
  - The amplitude reduction from the turning height to the ground is exp(−X), with X = ∫ q dz = (2/3) k d (a − θ²)^(3/2)/a and a = 2Δc/c.
  - At grazing (θ → 0), X = (2/3) k d √(2Δc/c). Writing Δc/c = d/R_c, this is exactly **X = (2/3)(d/l)^(3/2)** with l = (R_c/2k²)^(1/3), the same diffraction length that NICEOPS's ShadowTransition already uses. The loss is 8.69X dB ≈ **5.8 (d/l)^(3/2) dB**. It is also the asymptotic decay of Ai, so the WKB and exact Airy answers agree away from the turning point.
  - For ΔT = 3 K over d = 1 m (dc/dz ≈ −1.8 s⁻¹, R_c ≈ 195 m), l = 3.1 m at 100 Hz and 0.66 m at 1 kHz:

    | frequency | X | loss |
    |---|---|---|
    | 100 Hz | 0.12 | 1.1 dB |
    | 250 Hz | 0.31 | 2.7 dB |
    | 500 Hz | 0.62 | 5.4 dB |
    | 1 kHz | 1.24 | 10.7 dB |
    | 2 kHz | 2.47 | 21 dB |
    | 4 kHz | 4.94 | 43 dB |

  - For ΔT = 1 K over 1 m, scale X by 1/√3: 100 Hz 0.6 dB, 1 kHz 6.2 dB.
  - Halving d at fixed ΔT halves X (X ∝ d √ΔT at fixed frequency).
  - Arrivals steeper than grazing at the layer top tunnel *less*: X falls by the factor (1 − θ²/θc²)^(3/2).
- So the ray model's ~21 dB deficit at 8,000–13,000 ft is consistent with a model that sets the field to zero, or to a deep creeping-wave floor, where the physical field has lost only a few dB below ~500 Hz. The PE diffracting through the layer is the expected behavior, not a PE artifact. If the measured spectrum's deficit relative to the PE grows with frequency, roughly as f^(2/3), that would confirm the mechanism.
- **Cheap treatments, in order of cost** (inferences):
  1. **Turning-point continuation of the dipping family.** For a receiver below a ray's turning height, multiply the ray amplitude at the turning point by the Airy factor Ai(ζ_r)/Ai(0)-type continuation, or by exp(−∫q dz) away from the turning point, along the layered profile. The phase is unchanged (the evanescent field carries the turning ray's horizontal phase). This is the standard uniform caustic treatment applied to the lower caustic of the family. It adds one quadrature per pair over a layer NICEOPS already represents, and it reduces continuously to the lit ray as the receiver reaches the turning height.
  2. **Ground interaction below the layer.** The flush board sits on z = 0 inside the evanescent zone. The ground reflection there is an evanescent-wave reflection: the plane-wave coefficient evaluated at the complex angle cos θ = k_x/k(0) > 1. Weyl–van der Pol with that complex angle is the consistent continuation, and its correctness should be checked against the PE.
  3. **A turbulence floor.** Cheinet's result (shadow levels rise with |surface heat flux|) means that a warm, convective surface layer *also* scatters into its own micro-shadow. NICEOPS's ShadowScattering term (fitted γ_t) could be applied to the micro-shadow too, but the defensible γ_t is ≲ 3e-6 (CLAUDE.md), and an energy sum with the tunneled field may make it irrelevant at low frequency.
  4. **PE as the reference** for the case and for calibrating (1). It already agrees with the rays elsewhere within ±1 dB.
- A layer that is acoustically thin (d ≪ l) can simply be *removed* from the ray profile, merged into the first layer above it, for the bands where X < ~0.1 (under 1 dB). That gives a frequency-dependent "acoustic smoothing" of the bottom layer: below the diffraction length, the ray model should see a profile smoothed over l. This is a cheaper but cruder alternative to (1).

### Gaps
- No measured or PE results were found for this exact configuration (sub-meter warm layer under an inversion, elevated source, flush receiver). The tunneling estimate is mine and needs checking against the NICEOPS PE, ideally on the 286288 profile itself.
- Cheinet 2012's numbers (shadow-level dB against wind and heat flux) were not retrievable (403). The Roebuck, NJIT and Edinburgh tunneling sources are snippets.
- Whether Weyl–van der Pol at a complex (evanescent) incidence angle is valid for a flush receiver over σe = 200 kPa·s/m² ground was not established in any source.

---

## 4. Hybrid and fast methods for footprint-scale use: PE-based tables, ray–PE coupling, Gaussian beams, normal modes, interpolation over profile parameters, ML surrogates. What accuracy and cost are reported?

### Takeaway
For a footprint code the defensible options are physics-based and tabulated: PE excess-attenuation tables over a few fictive source heights (Cotté's method, as used by Colas et al.), interpolation in range and height, and low-dimensional profile-parameter interpolation or GPC surrogates. Full-wave LEE confirms the PE to 0.5 dB on flat ground but costs about 1000 times as much (1200 CPU-h against 1 CPU-h). Gaussian beams and complex rays are the ray-side fixes for shadows and caustics, with no published rms for atmospheric grazing. Neural surrogates have not reached sub-dB accuracy: PINNs for linearised Euler take hours of training for a 200 m domain, and data-driven PE surrogates remain at several dB.

### Cited Findings
- **LEE against PE for wind turbines.**
  - Flat ground: "The difference between the DL obtained with the three methods [LEE, WAPEvec, WAPEtopo] does not exceed 0.5 dB for x > 1 km".
  - Hill (31° slope): WAPEtopo is about 3 dB out per frequency, but OASPL is within 1 dB.
  - LEE cost: "approximately 1200 central processing unit (CPU) hours, which is equivalent to 3 days on a 16-core machine", using a 300 m × 300 m moving frame for a 3000 m × 300 m domain ("reduces the computational cost by one order of magnitude"). The PE at λ/10 "requires 1 CPU hour".
  - To limit cost, the Cotté (2019) method runs propagation only for "a number of fictive source heights distributed over a vertical line", which is a tabulation in source height.
  
  **[FULL]** — [Colas, Emmanuelli, Dragna, Blanc-Benon, Cotté & Stevens 2023, JASA 154, 1413](https://perso.ensta-paris.fr/~cotte/files/Colas_Emmanuelli_al_jasa2023%20Wind%20turbine%20sound%20propagation.pdf)
- Validation of a CNPE wind-turbine model against measurements at 178, 535 and 845 m: good agreement in OASPL and third-octave spectra except "in strong upwind situations", where levels at larger distances are overestimated "because turbulence is not considered". **[SN]** — [WES 8, 639 (2023)](https://wes.copernicus.org/articles/8/639/2023/)
- PE–Gaussian-beam hybrids: the analytic-ray-curve plus Gaussian-beam field runs "at close-to-interactive speed" for 3-D scenes and matches the Attenborough benchmark graphically. **[FULL]** — [Mo et al. 2017](https://gamma.cs.unc.edu/PAPERS/OutdoorSoundPropagation17.pdf)
- Neural PDE solvers for outdoor acoustics (DLR):
  - 1-D linearised Euler: a SIREN-PINN reaches MAE 1.85e-3 Pa after "1 h 2 min" of training; ELM-ODE reaches 7.40e-3 Pa in 2.98 s.
  - 2-D (±100 m domain): PINN "7h 36m 52s", ELM-ODE 1362 s, with MAE 6.2e-3 and 2.9e-3 at t = 0.2 s.
  
  **[FULL]** — [Juliust, Baumann, Schady, Gharbi & Dietrich 2025, Forum Acusticum](https://dael.euracoustics.org/confs/fa2025/data/articles/000602.pdf)
- The first PINN for ABL sound propagation (Pettit & Wilson, POMA 2020) adds physics penalties to the data loss to regularize a surrogate trained on less data. **[SN, no accuracy found]** — [search summary](https://labs.sonicfield.org/library/a-physics-informed-neural-network-for-sound-propagation-in-the-atmospheric-bound)
- CNNs and recurrent networks trained on simulated range-dependent wind and temperature fields predict infrasound transmission loss with "an error of 4 dB" and uncertainty estimates, at negligible cost compared with PE. **[SN]** — [CTBTO SnT 2024/2025 contributions](https://conferences.ctbto.org/event/27/contributions/5423/)
- GPC surrogates: see §2 (Martinelli 2023: 5 GPC samples against 100 LHS against 10,000 MC for small or medium variability).
- Earlier notes already hold Kayser 2022's kriging metamodel (WindTUNE), Pettit & Wilson's POD and cluster-weighted models, Hart et al.'s 5–7 dB against data, the HPM (PE + FFP + ray), and Waxler's normal modes. They are not repeated here.

### Inferences
- NICEOPS's split-step Padé PE (0.47 s for a five-band `pe_compare` case) is already cheaper than the PE costs in this literature by orders of magnitude, so the table approach Colas et al. use is within reach per run. The existing `--propagation_table` format (bearing × source height × receiver height × range × band) *is* the Cotté fictive-source-height method.
- Where to spend the PE: on a near-ground ensemble. Tabulate per run over 3–5 samples of the bottom-layer parameters (ΔT, d, or the stable-law choice). Then either interpolate in those parameters (they enter smoothly through X and θc; §2–3) or report the band. That costs 3–5 tables per atmosphere and gives the "ensemble of near-surface gradients scored as a band" recommended earlier.
- **Ray–PE coupling** (inference, no source read): trace rays from the source down to a height H (say 20–50 m) at a range r0 close to the receivers. Then hand the field on that vertical aperture to a PE as a starter for the last few hundred meters. That would confine the PE to the near-ground slab where diffraction and tunneling matter, while the elevated, wide-angle part stays with rays and the PE's 50–60° angle limit is never approached. Underwater acoustics uses similar ray/beam-to-PE handoffs, but I found no atmospheric implementation in this session.
- ML surrogates are not competitive for NICEOPS: training costs hours for toy domains, and the errors are several dB. The physics-based fast paths (turning-point continuation, PE tables, GPC over 2–3 parameters) are each cheaper and more accurate.

### Gaps
- No accuracy figures were found for Gaussian beams against PE at near-grazing over impedance ground (Gabillet et al. 1993 was not read).
- No published ray–PE vertical-aperture coupling for atmospheric acoustics was found.
- Pettit & Wilson's PINN accuracy is unknown, and the CTBTO 4 dB is for infrasound at 100+ km, so it does not transfer.
- Nord2000's meteorological approximations (equivalent linear profile, heuristic shadow) were not read: the FORCE PDF link returned HTML. Their accuracy against PE was not obtained.

---

## 5. Can NWP or reanalysis (HRRR, WRF LES, ERA5) fill the lowest layer, and what acoustic accuracy is reported?

### Takeaway
NWP and reanalysis are weakest in exactly this layer and regime. HRRR fails to form strong nocturnal inversions and under-represents low-level shear and static stability, and ERA5's largest near-surface errors are at night. Mesoscale-to-acoustic coupling has worked for low frequencies and short ranges (MIUU + OASES to 1.4 km, ≲100 Hz; SUBMESO + MW-WAPE; Meso-NH down to 1 m). No study found reports a dB accuracy for helicopter-type geometries. For the bottom 15 m, NWP is useful at most as a prior on the *class* (stable or transition), not as a profile.

### Cited Findings
- HRRR shows "large warm biases on clear-sky days resulting from the model's inability to form strong nocturnal inversions", attributed in part to its 3 km grid and to poorly resolved drainage flows. **[SN]** — [WRF workshop 2019 slides](https://www2.mmm.ucar.edu/wrf/site_linked_files/workshops/slides/2019/8.6.pdf)
- Against WFIP3 coastal lidar and profilers (1 year, 3 sites), HRRR captured seasonal and diurnal cycles, but "low-level horizontal wind shear and static stability were too weak in the model". **[SN]** — [egusphere-2026-97 preprint](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-97/egusphere-2026-97-manuscript-version3.pdf)
- ERA5: "near-surface layers, especially during nighttime, exhibited larger biases and errors … mainly due to reanalysis limitations in capturing surface cooling, temperature inversion". **[SN]** — [Remote Sensing 17, 393 (2025)](https://mdpi-res.com/d_attachment/remotesensing/remotesensing-17-00393/article_deploy/remotesensing-17-00393.pdf)
- ERA5 and NORA3 low-level jets were compared with FINO1 lidar through a five-parameter log–jet fit, with bias correction as the aim (bias correction is in the title; magnitudes not read). **[ABS]** — [Bui et al. 2026, WES](https://wes.copernicus.org/articles/11/2307/2026/)
- The MIUU mesoscale model simulated wind and temperature profiles over hilly terrain in eastern Norway as input to OASES. Out to 1.4 km and at ≲100 Hz, "acoustic predictions based on both measured and simulated atmospheric profiles agreed well with measured single-frequency sound levels in most cases". **[SN]** — [Hole & Mohr, "Forecasts of sound propagation by a mesoscale atmospheric model"](https://sonicfield.org/library/forecasts-of-sound-propagation-by-a-mesoscale-atmospheric-model). FFI reports note that weather-driven variability at long range is "tens of dB". **[SN]** — [FFI](https://ffi-publikasjoner.archive.knowledgearc.net/handle/20.500.12242/1055)
- SUBMESO micrometeorology was coupled to a mean-wind WAPE over moderately complex terrain. **[ABS]** — [Lihoreau, Gauvreau, Bérengier, Blanc-Benon & Calmet 2006, JASA 120, 110](https://hal.archives-ouvertes.fr/hal-00274779)
- Meso-NH run "at very fine scales (up to 1m)", with a drag-force canopy, was validated under stable, unstable and neutral conditions as input to TLM acoustics. **[SN]** — [Aumond et al. (Météo-France)](https://bibliotheque.meteo.fr/pub/DOC00025081-application-the-transmission-line-matrix-method-fo.html)
- A freeway-noise forecast used effective-sound-speed profiles from a 1 km NWP with near-ground vertical resolution "finer than 20 m". "Vertical nesting can improve the prediction of finer structures in near-ground temperature and velocity profiles, such as morning temperature inversions and low level jet-like features". WRF vertical nesting was used below 400 m AGL. **[SN]** — [ASU, "Simulating meteorological profiles to study noise propagation from freeways"](https://asu.elsevierpure.com/en/publications/simulating-meteorological-profiles-to-study-noise-propagation-fro); [Ovenden et al. 2015 (UCL)](https://discovery-pp.ucl.ac.uk/id/eprint/1461305/2/Ovenden.1461305_2015.pdf)
- LES-plus-PE studies (Wilson et al. 2008, via Wilson 2024) used LES as "ground truth" and found 6–8 dB inherent random errors for event predictions even with nearby accurate profiles. **[FULL]** — [Wilson 2024](https://acousticstoday.org/wp-content/uploads/2024/12/sound-predicting-wilson.pdf)

### Inferences
- With a lidar from 14 m and a sonde, the 2017 tests already have better upper-profile data than any NWP. NWP cannot help below the first gate, because its own lowest levels (typically ~10–20 m) and its stable-layer physics are what fails. At most, HRRR or ERA5 surface fluxes could help classify a run (stable, transition, convective) and set a prior on ζ for the §1 ensemble.
- The desert or airfield setting adds a further reason for caution: HRRR's documented failures are tied to unresolved drainage flows and terrain, and the morning-transition timing depends on entrainment from above (Angevine). Coarse NWP smooths both.

### Gaps
- No acoustic study was found that reports dB accuracy from HRRR, ERA5 or WRF profiles for elevated sources at 1–4 km. The Hole & Mohr and ASU results are snippets.
- ERA5 and HRRR near-surface bias magnitudes (K, m/s) at desert sites were not retrieved.
