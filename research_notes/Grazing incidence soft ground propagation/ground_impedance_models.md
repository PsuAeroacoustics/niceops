# Ground impedance models and ground characterisation at near-grazing incidence over soft soil

**Provenance key used throughout**

- **[FT]** full text read: Attenborough & Taherzadeh 2026 (Appl. Acoust. 242, 111067); Howe 1984 NASA CR-3862; Kruse dissertation (Oldenburg, summary volume); Nguyen et al. 2022 (arXiv 2205.13695); Nord2000 Road User's Guide; Nord2000 update paper (Acoustics 2011); Attenborough 2022 JASA "Reflections" note on Embleton et al.
- **[ABS]** abstract only (Crossref/PubMed): everything else, including Attenborough, Bashir & Taherzadeh 2011, Dragna, Attenborough & Blanc-Benon 2015, Li et al. 1998, Di & Gilbert 1993, Thomasson 1976, Hothersall & Harriott 1995. The full texts of these were paywalled or blocked.
- **[SNIP]** search-engine summary only (weakest): Dragna & Blanc-Benon 2014, Guillaume et al. 2015, Faure et al. 2017.
- **[OWN]** computed by me for this note (Python; script at the end of "Inferences" of the first question). Not from any source.

Sign convention everywhere: exp(-i omega t), Z = R + iX with X > 0 for porous ground, beta = 1/Z.

Not obtained at all: *Predicting Outdoor Sound* (2007/2021), Salomons' book, the Nord2000 and Harmonoise technical reports proper, ANSI S1.18 and NORDTEST ACOU 104 texts, CNOSSOS report. Anything about those is second-hand via the sources above.

---

## 1. Which impedance models are recommended for grassland/soil at low frequency?

### Takeaway
The Open University/Lyon line of work (2011, 2015, 2026) consistently says: drop one-parameter Delany–Bazley (and Miki as a layer) because they are physically inadmissible and fit soft ground worst; use a two/three-parameter physically admissible model (slit-pore, phenomenological/Zwikker–Kosten, or Wilson relaxation, with T = 1/sqrt(porosity)), hard-backed where the data demand it. The 2026 paper settles on slit-pore/Wilson with a class table; variable porosity is only "marginally" better than Delany–Bazley on grass and needs negative porosity rates at about a quarter of grass sites.

### Cited Findings

**Formulas (all [FT], Attenborough & Taherzadeh 2026, eqs. 3–12)** — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)

- Spherical-wave reflection (eq. 3): Q = Rp + (1 − Rp)F(w); F(w) = 1 + i sqrt(pi) w exp(−w²) erfc(−iw); w = sqrt(i k0 R2 / 2)(cos θ + β); Rp = (cos θ − β)/(cos θ + β). θ from the normal.
- Delany–Bazley (eq. 4): Zc = 1 + 9.08 (f/Rs)^−0.754 + i 11.9 (f/Rs)^−0.732; k = (ω/c0)[1 + 10.8 (f/Rs)^−0.70 + i 10.3 (f/Rs)^−0.595].
- Hard-backed layer (eq. 5): Z(d) = Zc coth(−i k d).
- Miki three-parameter (eq. 6): Zc = (sqrt(T)/Ω)[1 + 0.070 (ω/2πRs)^−0.632 + i 0.107 (ω/2πRs)^−0.632]; k = (ω sqrt(T)/c0)[1 + 0.109 (ω/2πRs)^−0.618 + i 0.160 (ω/2πRs)^−0.618].
- Zwikker–Kosten: Zc = (sqrt(T)/Ω) sqrt(1 + iΩRs/(ω ρ0 T)), k = Ω Zc ω/c0 (as printed).
- Phenomenological (modified ZK, isothermal pores): Zc = [1/(Ω sqrt(γ))] sqrt(T + i Rs Ω/(ω ρ0)), k = Ω sqrt(γ) Zc ω / c0... (printed as "k = Ω γ^(1/2) Zc ω/c0").
- Hamet–Bérengier: Zc = (1/Ω)(K/γ)^(1/2) {1 + [(γ−1)/γ](1/F0)}^(−1/2) Fμ^(1/2), Fμ = 1 + iωμ/ω, F0 = 1 + iω0/ω, ωμ = (Rs/ρ0)(Ω/K), ω0 = ωμ(K/N_Pr). Recommended in HARMONOISE for porous asphalt only.
- Variable porosity (eq. 10): Z = a(1 + i) sqrt(Rse/f) + i b αe / f, a = 1/sqrt(π ρ0 γ), b = c0/(8πγ), Rse = Rs/T, αe = 4/de (de an effective depth). This matches the expression already in `src/ground_effect.cpp`.
- Slit pore (eq. 11): Zc = [1/(ρ0 c0)] sqrt[(T/Ω²) ρ(ω)/C(ω)], k = ω sqrt[T ρ(ω) C(ω)], ρ(ω) = ρ0/H(λ), C(ω) = (γ P0)^−1 [γ − (γ−1) H(λ sqrt(N_Pr))], λ² = 3 ω T ρ0/(Ω Rs), H(λ) = 1 − tanh(λ sqrt(−i))/(λ sqrt(−i)).
- Wilson relaxation (eq. 12): Z = (sqrt(T)/Ω)[(1 + (γ−1)/sqrt(1 − iωτe))(1 − 1/sqrt(1 − iωτv))]^(−1/2); k = (ω sqrt(T)/c0)[(1 + (γ−1)/sqrt(1 − iωτe)) / (1 − 1/sqrt(1 − iωτv))]^(1/2); τv = 2 ρ0 T/(Ω Rs), τe = N_Pr τv.
- "Acceptable fits to short range data have been obtained by assuming T = 1/sqrt(Ω)", which cuts a hard-backed layer to three parameters, two if depth does not matter.

**What the three papers conclude**

- [FT] AT2026's objections to the HARMONOISE/CNOSSOS class table: one-parameter model fits low-flow-resistivity surfaces (classes A–C) relatively poorly; hard-backed fits still worse than alternatives; single σ per class ignores in-class variation; the model is inadmissible because "(a) it predicts negative real parts of complex density and surface impedance at low frequencies, and (b) it is non causal"; no account of moisture; no account of roughness. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Average fitting errors to short-range level-difference spectra (AT2026 Table 2, summarizing the 2011 paper). 26 grasslands: Delany–Bazley 9.3, Zwikker–Kosten and slit pore 8.7, variable porosity 6.7. 12 grasslands better fitted as layers: DB layer 7.5, ZK layer 7.3, slit-pore layer 7.0, variable porosity 7.6. 5 pine forest floors: DB 13.8, ZK 5.0, slit 5.0, variable porosity 10.8. 3 beech floors: DB layer 22.6, VP 24.8, ZK layer 13.8, slit layer 11.6. 4 gravel/sand: DB layer 15.2, ZK layer 9.4, slit layer 9.4, VP 9.8. (Error measure is eq. 15, a deviation normalized by the measurement standard deviation; the column is headed "dB".) — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] "Although the variable porosity model is marginally more accurate than the one parameter model for estimating the surface impedance of some grass-covered surfaces, best fits ... for 6 of the 26 grassland sites and for a bare soil site, have been obtained with negative values of the porosity rate or of 'effective depth'", and comparable fits come from slit-pore and ZK there. Negative porosity rates were also needed for dry and wet sand in an impedance tube and for several forest floors. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Miki three-parameter: similar to the others at low σ/high porosity, but at higher σ and lower porosity "its predictions of the real part of the surface impedance differ significantly from those of other models" (Fig. 2: 100 kPa s m⁻², Ω = 0.4, d = 3 cm). Wilson relaxation is "very similar to those of the slit pore and phenomenological models" and is the one carried forward. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [ABS] 2011 paper: 42 short-range level-difference spectra; "For all but eight sites, slit-pore, phenomenological and variable porosity models yield lower fitting errors than those given by the widely used one-parameter semi-empirical model. Data for 12 of 26 grassland sites and for three beech wood sites are fitted better by hard-backed-layer models." Fitted parameters for forest floors, porous asphalt and gravel are consistent with non-acoustic measurements. — [Attenborough, Bashir & Taherzadeh 2011](https://doi.org/10.1121/1.3569740)
- [ABS] 2015 paper: one-parameter models "are not physically admissible. Even when corrected to satisfy a passivity condition in respect of surface impedance they do not satisfy the condition that the real part of complex density must be greater than zero." Variable porosity is better for many grasslands, two-parameter slit pore for forest floors and gravel; slit pore is shown to be admissible. Crucially for this project: "differences between frequency domain predictions at longer ranges using physically admissible models rather than one parameter models are significantly greater than those resulting from short range spatial variability and comparable with seasonal variability over grassland." — [Dragna, Attenborough & Blanc-Benon 2015](https://pubmed.ncbi.nlm.nih.gov/26520321/)
- [SNIP] Dragna & Blanc-Benon 2014: admissibility means the model extends to the whole complex frequency plane and its inverse transform is real, causal and passive; they propose a "modified Miki" model that is admissible. — [Dragna & Blanc-Benon 2014](https://doi.org/10.3813/aaa.918719)
- [ABS] Taraldsen & Jonasson 2011: a two-parameter Darcy-law model (effective flow resistivity + effective layer depth) fits field data 200 Hz–2.5 kHz better than Delany–Bazley and is causal; "particularly relevant for extrapolation of the models to lower frequencies", the harmonised methods going down to the 25 Hz band with Delany–Bazley as default. — [Taraldsen & Jonasson 2011](https://doi.org/10.1121/1.3500694)
- [ABS] Attenborough 1992: "Further evidence is produced of the unreliability that would be associated with extrapolation of fits based on the single-parameter semiempirical model ... outside the frequency range of the fitted data." — [Attenborough 1992](https://doi.org/10.1121/1.404251)
- [FT] Kruse: Delany–Bazley relations "are considered valid for 0.01 < f/σ < 1; this makes the use of these relationships for low frequencies and high flow resistivities commonly occurring in natural grounds questionable." (σ = 200 kPa s m⁻² at 63 Hz is f/σ = 3e-4 in those units, far outside.) A positive αe gives X > R, like a hard-backed layer; R > X at low frequency indicates "a compacted / encrusted upper layer". — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [ABS] Raspet & Sabatier 1996 is the origin of the corrected exponential-porosity model (earlier versions gave non-physical results for porosity decreasing with depth). — [Raspet & Sabatier 1996](https://doi.org/10.1121/1.414498)

**Typical parameter values**

- [FT] AT2026 Table 4, proposed classes (porosity; effective flow resistivity kPa s m⁻²; effective layer thickness m):
  1. Railway ballast, large gravel: 0.4 ± 0.1; 0.02 ± 0.01; 0.2 ± 0.1
  2. Porous asphalt: 0.2 ± 0.1; 8 ± 6; 0.04
  3. Snow: 0.6 ± 0.2; 15 ± 10; 0.15 ± 0.1
  4. Forest floors, coarse sand: 0.5 ± 0.1; 50 ± 40; 0.06 ± 0.04
  5. Uncompacted grass (pasture, meadow): 0.25 ± 0.15; 250 ± 100; 0.05
  6. Fine sand (dry): 0.3 ± 0.1; 300 ± 150; 0.05
  7. Compacted soil or "maintained" grass (lawns, verges, sports fields; compacted by mowing, trafficking, rolling): 0.2 ± 0.1; 1500 ± 700; infinite
  8. Compacted gravel road, hot rolled asphalt, sealed concrete: ≤ 0.1; ≥ 20,000; infinite
  — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] AT2026 Table 3, individual fits (porosity; σ; depth): pine forest 0.43; 26.7; ∞. Beech floor 0.35; 31.0; 0.070. Long grass 0.65; 51.0; 0.046. Lawn 0.32; 413; ∞. Sports field 0.22; 663.8; ∞. Arable 0.5; 2250; ∞ ("suggests that it might have had a surface crust due to weathering"). Snow 0.774; 9.6; 0.05 (measured, not fitted). — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] HARMONOISE/CNOSSOS/Nord2000 classes A–H: 12.5, 31.5, 80, 200, 500, 2000, 20,000, 200,000 kPa s m⁻² (C = "uncompacted, loose ground (turf, grass, loose soil)", D = "normal uncompacted ground (forest floors, pasture field)", E = "compacted field and gravel (compacted lawns, park area)"). — [Nord2000 Road User's Guide, Table 10](https://egra.cedex.es/EGRA-ingles/I-Documentacion/National_Methods/Users_Guide_Nord2000_Road.pdf)
- [ABS] Embleton, Piercy & Daigle 1983 (Delany–Bazley effective σ): new snow 10–30, grass 150–300, mature asphalt 30,000 cgs rayls (= kPa s m⁻²); the viscous/thermal boundary layer caps even sealed concrete at 1e5–1e6. — [Embleton et al. 1983](https://doi.org/10.1121/1.390029)
- [FT] Howe fits "grass covered ground" with R = 1 + 0.051 (κ/f)^0.75, X = 0.077 (κ/f)^0.73, κ = 230,000 MKS, giving |1/Z| = 0.0065 at 10 Hz and 0.035 at 100 Hz. — [Howe 1984, NASA CR-3862](https://ntrs.nasa.gov/api/citations/19850008280/downloads/19850008280.pdf)

### Inferences

- **[OWN] Model choice alone is worth 2–4.5 dB at 2° elevation in the 250–1000 Hz bands, in the direction of the observed deficit.** Level re free field (dB, third-octave averaged) at a flush receiver, range 1.5 km, homogeneous ground, full Weyl–van der Pol:

  | elev | model (σ kPa s m⁻² / Ω) | 63 | 125 | 250 | 500 | 1000 Hz |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1° | Delany–Bazley 200 | 2.0 | −4.6 | −8.3 | −11.8 | −15.4 |
  | 1° | slit pore 250 / 0.25 | 1.9 | −2.6 | −5.7 | −8.0 | −10.4 |
  | 2° | Delany–Bazley 80 | −1.1 | −4.4 | −7.7 | −11.0 | −14.3 |
  | 2° | Delany–Bazley 200 | 2.6 | −0.4 | −3.3 | −6.5 | −9.9 |
  | 2° | Delany–Bazley 500 | 4.5 | 2.7 | 0.4 | −2.3 | −5.4 |
  | 2° | Miki 200 (T = Ω = 1) | 0.1 | −2.9 | −5.6 | −8.5 | −11.4 |
  | 2° | slit pore 250 / 0.25 (class 5 center) | 2.6 | 0.6 | −1.3 | −3.3 | −5.4 |
  | 2° | same, hard-backed 5 cm | 5.5 | 1.7 | −1.8 | −3.9 | −5.4 |
  | 2° | slit pore 150 / 0.4 | 0.3 | −2.0 | −4.1 | −6.4 | −8.9 |
  | 2° | slit pore 350 / 0.1 | 4.3 | 3.0 | 1.6 | 0.1 | −1.2 |
  | 2° | slit pore 1500 / 0.2 (class 7) | 4.9 | 4.1 | 3.0 | 1.8 | 0.3 |
  | 2° | variable porosity σe 100, αe 50 /m | 2.2 | −0.7 | −3.1 | −5.6 | −8.2 |
  | 5° | Delany–Bazley 200 | 4.3 | 3.2 | 1.5 | −0.8 | −3.5 |
  | 5° | Miki 200 | 3.3 | 1.8 | 0.0 | −2.3 | −4.8 |
  | 5° | slit pore 250 / 0.25 | 4.3 | 3.6 | 2.6 | 1.3 | −0.2 |
  | 5° | slit pore 1500 / 0.2 | 5.4 | 5.1 | 4.7 | 4.2 | 3.4 |
  | 10° | Delany–Bazley 200 | 5.1 | 4.5 | 3.6 | 2.2 | 0.2 |
  | 10° | slit pore 250 / 0.25 | 5.1 | 4.7 | 4.2 | 3.5 | 2.5 |

  The phenomenological model at 250 / 0.25 is within 0.1 dB of slit pore everywhere. The variable-porosity parameters are illustrative, not taken from a published fit. At equal nominal flow resistivity the slit-pore impedance is larger than Delany–Bazley's from 125 Hz up (e.g. 14.4 + 13.2i against 8.7 + 10.1i at 250 Hz) because of the 1/Ω factor, so the ground is "harder" and the flush level higher.
- **[OWN] For an elevated source and a flush receiver the ground term depends on elevation angle, not on range.** Repeating the table at 3 km gives the same numbers to 0.1 dB at 2–10° (only the 63 Hz band at 1° differs, by up to 2 dB). Numerical distance is already large (|w| = 1.7–26 at 2°, 63–1000 Hz, DB 200, 1.5 km), so the field is the plane-wave one: level ≈ 20 log|2 sin ψ/(sin ψ + β)|, which I checked against the full formula to 0.1 dB wherever |w| > 3.
- **[OWN] Consequences of that.** (a) The flush-microphone sensitivity is at most −20 log|β|, not −40 log|β|; the fourth-power law belongs to source and receiver both on the ground. Doubling σ (slit pore, Ω = 0.25) raises the flush level by 1.5–3.0 dB at 1°, 1.1–2.2 dB at 2°, 0.5–1.5 dB at 5°, 0.3–1.0 dB at 10°. (b) A measured decay that is 1–2 dB/km slower than predicted cannot come from a homogeneous ground term at fixed elevation. With a fixed source height, elevation falls with range, and then the slope is d(level)/dψ: from 5° to 2° the 250 Hz band drops 4.8 dB with DB 200, 3.9 dB with slit 250/0.25 and 1.7 dB with slit 1500/0.2. A harder effective ground than modeled flattens the decay in exactly the 125–630 Hz bands.
- The model families divide cleanly for this use: semi-infinite slit-pore/phenomenological/Wilson have R ≈ X at low frequency; Delany–Bazley, Miki, variable porosity with positive αe, and any hard-backed layer have X > R. Which side the site is on is what decides the low-frequency behavior, and AT2026 says roughly half of grass sites sit on each side (12 of 26 better as layers).
- The [OWN] script is `calc.py`/`calc2.py` in the session scratchpad (`.../scratchpad/gi/`); it implements AT2026 eqs. 3–12 as transcribed above and nothing else.

### Gaps
- Full text of the 2011 and 2015 JASA papers was not accessible, so per-site fitted values (the 42-site table), the variable-porosity parameter ranges for grass, and the 2015 paper's longer-range dB differences are not quoted. Only its abstract's qualitative statement is.
- Horoshenkov's Padé/log-normal pore-size model: nothing retrieved. No source found applying it to outdoor ground at low frequency.
- Attenborough's original three- and four-parameter models: formulas not retrieved; AT2026 does not carry them forward, preferring models needing at most three parameters.
- Nothing found specific to dry lake bed (playa) or desert soil impedance. The "arable, surface crust, 2250 kPa s m⁻²" fit is the nearest analogue, and is my analogy, not a source's.

---

## 2. Extended reaction versus local reaction

### Takeaway
Local reaction is adequate for natural ground at σ of order 100 kPa s m⁻² and above; it fails for snow, ballast, porous asphalt, forest litter and thin low-σ layers. For grass/soil the error near grazing is under 0.1 dB, so this is not a candidate for the 2 dB deficit.

### Cited Findings
- [FT] Extended-reaction effective admittance, to be used in place of β in the same sound-field model (Kruse eqs. 35–37, ψ the grazing angle, n the complex index of refraction k/k0, m = ρ0/ρm = 1/(n Zc)): semi-infinite βe = m sqrt(n² − cos²ψ); hard-backed layer of depth d: βe = −i m sqrt(n² − cos²ψ) tan(k0 d sqrt(n² − cos²ψ)) (sign as for exp(−iωt); the equation was partly garbled in text extraction). Normal incidence of the layer form reduces to Ze = i Zc cot(k d). — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [FT] Local reaction can be expected for near-normal incidence, a high index of refraction (high flow resistivity), or high loss in the medium. "Most grounds have flow resistivities of 100 krayl/m and above ... It is therefore, in general, adequate to assume local reaction." Kruse's Figs. 14–15: at 10 kPa s m⁻² the semi-infinite extended/local difference grows slowly with frequency, and a 10 cm hard-backed layer differs markedly below about 500 Hz; at 100 kPa s m⁻² "the effect of the extended reaction is negligible". — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [ABS] Li, Waters-Fuller & Attenborough 1998: "heuristic modifications of the classical analytical approximation for the field due to a point source above an impedance plane are shown to be adequate for predicting sound propagation from a point source close to layered porous media with extended reaction. The main deficiencies are for extremely low flow resistivities and at close ranges." Multiple layering explains earlier difficulties over snow. — [Li et al. 1998](https://doi.org/10.1121/1.423307)
- [ABS] Li & Tao 2014: the plain WvdP form "is inadequate for predicting sound fields in the vicinity of non-locally reacting surfaces; a correction term is often required"; with it there is a singularity when source is directly above the receiver, which they remove. — [Li & Tao 2014](https://doi.org/10.1121/1.4836275)
- [ABS] Attenborough, Hayek & Lawther 1980: exact solution over a porous half-space as a branch-line integral plus a pole residue; steepest descent with pole subtraction gives the approximate solution for an "externally reacting" boundary, and "the precise nature of the further approximations that follow the assumption of a locally reacting (or impedance) boundary is made explicit". — [Attenborough et al. 1980](https://doi.org/10.1121/1.385074)
- [ABS] Nicolas, Berry & Daigle 1985 (snow, 5–50 cm, ranges to 15 m): layer theory agrees above about 200–500 Hz; "below these frequencies the measurements suggest that the layer of snow behaves as semi-infinite porous ground." — [Nicolas et al. 1985](https://doi.org/10.1121/1.391902)
- [FT] Railway ballast (σ about 5 kPa s m⁻²): "sufficiently low that allowance should be made for external reaction ... but Fig. 3 shows that reasonable agreement is obtained by assuming local reaction" at 7.5 m. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [ABS] 2011 paper: over ballast, "if extended reaction is taken into account, the hard-backed-layer version of the slit-pore model gives the most reasonable parameter values." — [Attenborough et al. 2011](https://doi.org/10.1121/1.3569740)
- [ABS] Porous road surface: local reaction "is not valid"; using the grazing-incidence effective impedance in a single-domain BEM is much more accurate than the normal-incidence impedance. — [Anfosso-Lédée et al. 2007](https://doi.org/10.1121/1.2749407)

### Inferences
- **[OWN] Magnitude at the geometry of interest** (βe substituted into both Rp and w; flush receiver and 1.2 m pole; 1–10° elevation; 63–1000 Hz third octaves): slit pore 250 kPa s m⁻² / Ω 0.25 semi-infinite: ≤ 0.1 dB everywhere, 0.2 dB at worst at the pole's 500 Hz dip. Hard-backed 6 cm layer at 50 / 0.5 (forest-floor class): ≤ 0.4 dB flush, up to 0.8 dB at the pole's interference dip at 250 Hz. 15 cm layer at 15 / 0.6 (snow class): ≤ 1.0 dB flush, up to 1.2 dB. So extended reaction matters at the 0.5–1 dB level only for class 3–4 ground, and not at all for grass or soil.
- The effect has no systematic sign in those numbers.

### Gaps
- No quantitative published number for the local/extended difference in dB at km range was found; the 1998 paper's figures were not accessible.
- The exact form of Li et al.'s "heuristic modification" (whether a correction term beyond substituting βe is needed) was not read. Li & Tao's abstract says a correction term "is often required".

---

## 3. Surface waves over soft and layered ground

### Takeaway
A surface wave exists when the surface is sufficiently "springlike" (X > R, as for thin hard-backed layers or rough ground); outdoors over natural ground it has been confirmed only over thin snow. The Faddeeva form of F(w) contains the pole term exactly, provided the complex root and half-plane are handled correctly, and for an elevated source at ≥ 2° over grass it is not excited at all.

### Cited Findings
- [FT] "Embleton et al identify a component in the solution for near-grazing total sound field above a finite impedance plane as a cylindrically spreading airborne surface wave ... Since its attenuation depends on the real part of ground impedance, the only observation of the surface wave outdoors over a naturally occurring surface has been over a thin snow layer." It is "related to that created over a rough surface near-grazing incidence". — [Attenborough 2022, JASA 151, R11](https://pubs.aip.org/asa/jasa/article-abstract/151/5/R11/2839449/Pioneering-study-of-outdoor-sound-propagation)
- [ABS] Albert 2003: pistol shots 1 m above snow, 60 m range, microphones 0.1–4.75 m. "A strong, low frequency 'tail' following the initial arrival"; exponential decay with height, measured α = 0.28 m⁻¹; confirmed against "the explicit surface wave pole term using residue theory". — [Albert 2003](https://doi.org/10.1121/1.1559191)
- [ABS] Raspet & Baird 1989: the term "is a true surface wave that can exist independent of the body wave in air." — [Raspet & Baird 1989](https://doi.org/10.1121/1.397588)
- [ABS] Daigle, Stinson & Havelock 1996: over a lattice of cavities the surface wave is a distinct arrival with phase speed below c and exponential decay with height; pulses computed from the measured impedance agree with measurement. — [Daigle et al. 1996](https://doi.org/10.1121/1.415386)
- [ABS] Thomasson 1976: the exact solution is "a single integral along a steepest-descent contour and a Hankel function"; Wenzel "found that in a special case expansions of his and Ingard's solution differ by a 'surface wave' term." — [Thomasson 1976](https://doi.org/10.1121/1.380943)
- [ABS] Stinson 1995, on the Chien–Soroka formula: "the correct choice of a complex root is specified. This choice is made by ensuring that the sound pressure be a bounded and continuous function of the geometrical parameters and the surface impedance." — [Stinson 1995](https://doi.org/10.1121/1.413381)
- [ABS] Ochmann 2004: Sommerfeld's real line-image form "converges only for masslike impedances and cannot be used for the more important case of ... springlike surface impedances"; the singular part becomes a Hankel function, "which shows that surface waves are contained in the whole solution". — [Ochmann 2004](https://doi.org/10.1121/1.1819504)
- [ABS] Bashir, Taherzadeh & Attenborough 2013: closely spaced strips act as a hard-backed slit-pore layer and wider spacings as periodic roughness; both generate surface waves of comparable magnitude. — [Bashir et al. 2013](https://doi.org/10.1121/1.4824846)
- [FT] Howe: with source and receiver on rough ground at low frequency, surface irregularities generate a surface wave that delays the onset of the ground shadow "to such an extent that the field exceeds that obtaining for a smooth, rigid plane. At larger distances the surface wave is absorbed by the ground." He cites Piercy et al. that beyond a mile the main contribution is "from surface waves of frequencies f ≈ 200 Hz" (the symbol before 200 is lost in the scan). — [Howe 1984, NASA CR-3862](https://ntrs.nasa.gov/api/citations/19850008280/downloads/19850008280.pdf)

### Inferences
- **[OWN] Condition in this formulation.** With w = sqrt(i k0 R2/2)(sin ψ + β), Im w < 0 requires −Im β − Re β > sin ψ. Only then does the reflection identity w(z) = 2 exp(−z²) − w(−z) bring in the 2 exp(−w²) residue, which is the surface wave. Margins −Im β − Re β that I computed:
  - Delany–Bazley 200: +0.004 to +0.008 (63–500 Hz). Semi-infinite slit pore / phenomenological: negative at all frequencies (never). Variable porosity (100, 50 /m): +0.008 to +0.012.
  - Slit 250/0.25 hard-backed 5 cm: +0.015, +0.020, +0.012 at 63, 125, 250 Hz.
  - Slit 50/0.5 hard-backed 6 cm: +0.042, +0.068, +0.078 at 63, 125, 250 Hz.
  - Against sin ψ = 0.017, 0.035, 0.087, 0.174 at 1°, 2°, 5°, 10°.

  So over grass-like ground the pole is reached only below about 1° and only for a hard-backed layer; over forest-floor-like layers it is reached up to about 2–4° in the 63–250 Hz bands. Where it is reached the flush level exceeds pressure doubling: I get +6.6 dB at 63 Hz, 1°, 1.5 km for the 5 cm grass layer, and +9 to +10 dB at 31.5–63 Hz, 0.5°, for the soft layer.
- **[OWN] The standard formulation captures it correctly for local reaction.** Against direct numerical evaluation of the exact complex-image integral (next question) the Weyl–van der Pol result including those surface-wave cases agrees to better than 0.01 dB. NICEOPS's `faddeeva` documents the lower-half-plane reflection as safe only in its tested domain (Im z ≥ −1.2); thin soft layers at very low elevation are where Im w goes negative, so that domain limit is the thing to check, not the physics.
- For extended reaction and multilayer ground the pole is no longer at cos θ = −β with a constant β; substituting βe(ψ) evaluated at the specular angle puts the pole only approximately. This is my reading of why Li et al. call the modification "heuristic" and find it worst at extremely low σ and short range; I have not seen their derivation.

### Gaps
- No full-text source for pole location in layered/extended-reaction ground was read. Attenborough–Hayek–Lawther 1980 is the primary reference (abstract only).
- No outdoor observation of a surface wave over grass or soil from an elevated source was found; the literature is ground-to-ground.

---

## 4. Validity of Weyl–van der Pol / Chien–Soroka, and exact alternatives

### Takeaway
For locally reacting ground at these ranges (k r of 10³–10⁵) the approximation is not a source of error: I find ≤ 0.01 dB against the exact image integral for grass-like impedance and ≤ 0.15 dB even for snow-like layers. The documented failures are at k r of order 1 with very small |Z|, and for non-locally-reacting ground.

### Cited Findings
- [ABS] Di & Gilbert 1993: "By representing the plane-wave reflection coefficient as the Laplace transform of an image source distribution, a well-behaved image integral, instead of the usual Sommerfeld integral, is obtained. The approach is valid for both locally and extended reacting surfaces. For a locally reacting ground surface, the image integral is an especially simple, rapidly convergent integral." An upper integration limit is given that works "with realistic values of ground impedance". — [Di & Gilbert 1993](https://doi.org/10.1121/1.405435)
- [ABS] Gu et al. 2025: Di & Gilbert's formulation "is an efficient and accurate method"; its upper limit had "not been examined in the context of an extended reacting reflecting plane", and they give a method to choose it. — [Gu et al. 2025](https://pubmed.ncbi.nlm.nih.gov/39791994/)
- [ABS] Thomasson 1976: exact solution as one steepest-descent integral plus a Hankel function that "could easily be integrated numerically". — [Thomasson 1976](https://doi.org/10.1121/1.380943)
- [ABS] Nobile & Hayek 1985: asymptotic series with recursion for higher terms, "extremely accurate, even for relatively low values of the parameter kR. Not surprisingly, the plane-wave solution often shows major deviations from the exact integral solution." — [Nobile & Hayek 1985](https://doi.org/10.1121/1.392902)
- [FT] Kruse, summarizing Nobile & Hayek and later comparisons: "For not too small surface impedances (e.g. Z = (2+3i)), it will be practically identical to the exact solution for k·r > 0.1, even for near grazing incidence. Only for very low impedances (e.g. Z = (0.3+0.5i)) one has to anticipate differences ... for k·r < 3." "If dealing with typical ground impedances which are high at low frequencies, the approximate solutions are usually indistinguishable from those of the more exact expressions." — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [ABS] Attenborough, Hayek & Lawther 1980: "the existing solutions for a locally reacting boundary are compared, and with one exception they are shown to be identical." — [Attenborough et al. 1980](https://doi.org/10.1121/1.385074)
- [ABS] Ochmann 2004: complex equivalent sources give a half-space Green's function valid for both masslike and springlike impedance, where the real line-image form fails for springlike. — [Ochmann 2004](https://doi.org/10.1121/1.1819504)
- [SNIP] Taraldsen 2005, "The complex image method" (Wave Motion 43, 91–97): a simple derivation, particularly simple for an impedance boundary and valid for reflection from an arbitrary layered medium. No abstract was available. — [Taraldsen 2005](https://doi.org/10.1016/j.wavemoti.2005.07.001)
- [ABS] Kawai, Hidaka & Nakajima 1982, J. Sound Vib. 83, 125–138, exists as cited; no abstract retrieved. — [Kawai et al. 1982](https://doi.org/10.1016/s0022-460x(82)80080-1)

### Inferences
- **[OWN] Direct check.** I evaluated p_refl = exp(ikR2)/R2 − 2kβ ∫₀^∞ exp(−kβq) exp(ikR(q))/R(q) dq, R(q) = sqrt(r² + (H + iq)²), which follows from writing 1/(cos θ + β) as a Laplace integral (the Di & Gilbert / complex-image construction; convergent for Re β > 0), by trapezoid quadrature, and compared it with Q from eq. 3. Flush receiver, r = 300, 1500, 5000 m, ψ = 0.5°, 2°, 10°, f = 31.5, 63, 250, 630 Hz:
  - Delany–Bazley 200; slit 250/0.25 hard-backed 5 cm; the same with 1 cm roughness: |difference| ≤ 0.011 dB in all 108 cases.
  - Slit 50/0.5 hard-backed 6 cm: ≤ 0.005 dB.
  - Slit 15/0.6 hard-backed 15 cm (snow-like, local reaction): ≤ 0.03 dB except 0.08 and 0.15 dB at 63 Hz, 300 m, 0.5° and 2°.
- The approximation's error scales with 1/(k R2), and k R2 ≥ 350 here. There is nothing to gain from an exact evaluation on the broadband or spectral path for local reaction. It would matter only if a multilayer extended-reaction ground were introduced.

### Gaps
- Published dB comparisons of WvdP against exact solutions at long range and grazing (as opposed to short range) were not found; the comparison above is mine.
- Kawai 1982 and Taraldsen 2005 contents are known only by title.

---

## 5. Roughness at grazing: effective admittance

### Takeaway
Sub-wavelength roughness is not a loss of coherence at grazing; it changes the effective admittance: a term proportional to k0 times roughness volume per unit area is added to the imaginary part (boss models), or equivalently the real part of the impedance drops while the imaginary part stays. That makes the ground effectively softer and more surface-wave-prone, the opposite of a Kirchhoff factor that vanishes at grazing.

### Cited Findings
- [FT] "Whether the surface is porous or not, a mean roughness height much smaller than the incident sound wavelengths, alters the effective surface impedance near grazing incidence since a high proportion of the incident sound energy is scattered, coherently, in the direction of specular reflection. The predicted effect ... is to reduce real part of the effective surface impedance while leaving the imaginary part relatively unaffected." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Approximate formula for rms roughness height ⟨H⟩ (eq. 16, as I read it from a garbled layout): Zr ≈ Zs − (⟨H⟩σ/(γ ρ0 c0)) (2/ν − 1), ν = 1 + (2/3)π⟨H⟩, with Re(Zr) = 0 if the correction exceeds Re(Zs). "The influence of small-scale roughness is predicted to be largest when the ground surface is acoustically hard or has a high flow resistivity." Proposed for mean roughness heights below 0.1 m. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Alternative: "choosing a lower flow resistivity value within the range for the appropriate ground surface class. The effective flow resistivity of the surface of a field was reduced by a factor of 16 after roughness heights of between 0.07 and 0.2 m were created by ploughing" (citing Chambers & Sabatier 2002). — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] "If the roughness size is comparable with the wavelengths of interest (≥ 0.3 m), then non-specular (incoherent) scattering dominates ... A rough hard ground with ⟨H⟩ ≥ 0.3 m can be assumed to increase the sound level by a frequency-independent 3 dB." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Howe's coherent boundary condition for hemispherical bosses of radius a on a locally reacting plane of impedance Z (eqs. 2.27–2.28): ∂φ/∂x2 = −i k0 η* φ, where η* is the flat-ground admittance 1/Z, scaled by a factor of order (1 + const·σ), plus a term −i·(const)·σ k0 a, σ being the fractional area covered; for semicylinders (eq. 2.31) the same structure with the roughness term multiplied by cos² of the angle between propagation and the normal to the ridges (eq. 3.14), so ridges parallel to propagation act as smooth ground of impedance Z/(1 + const·σ). The numerical constants are illegible in the scan. Valid for σ ≪ 1 and wavelengths large against the roughness; "expected to be adequate at frequencies less than about 1000 Hz". — [Howe 1984, NASA CR-3862](https://ntrs.nasa.gov/api/citations/19850008280/downloads/19850008280.pdf)
- [FT] Howe: for grass at f < 450 Hz, |1/Z| < 0.1, so "the contributions to η* from Z and the surface irregularities can be of the same order of magnitude." When k0 a is large the roughness has "a significant effect on propagation". — [Howe 1984](https://ntrs.nasa.gov/api/citations/19850008280/downloads/19850008280.pdf)
- [ABS] Attenborough & Taherzadeh 1995: modified Howe and Tolstoy theories agree with each other and tolerably with data over hard rough boundaries; the modified Tolstoy form, which lets the bosses have their own impedance, is better over rough finite-impedance surfaces. — [Attenborough & Taherzadeh 1995](https://doi.org/10.1121/1.414454)
- [ABS] Attenborough & Waters-Fuller 2000: extended boss theory (incoherent scatter, randomness, finite impedance) agrees with grazing-angle impedance measurements over rough sand. "The rough surface impedance model may be used, as an alternative to layering, to explain measurements where the reactance exceeds the resistance at low frequencies", and it predicts features measured "with an elevated source and vertically separated receivers at a range of 30 m over cultivated ground surfaces." — [Attenborough & Waters-Fuller 2000](https://doi.org/10.1121/1.1288940)
- [ABS] Boulanger et al. 1998: Twersky boss effective impedance plus a heuristic grating correction is "in tolerable agreement with the data except for close and random packing". — [Boulanger et al. 1998](https://doi.org/10.1121/1.424358)
- [ABS] Boulanger, Attenborough & Qin 2005: Twersky's formulation extended heuristically to porous roughness; the modified Tolstoy and modified Twersky theories "give different predictions for relatively large bosses"; modified Twersky is good over sand with semi-elliptical roughness and over uncultivated soil. — [Boulanger et al. 2005](https://doi.org/10.1121/1.1850211)
- [ABS] Nichols et al. 2011: roughness size recovered from complex excess attenuation with a simplified boss model to within 4%, up to 32% packing density. — [Nichols et al. 2011](https://doi.org/10.1121/1.3556897)
- [SNIP] Faure et al. 2017: an effective impedance built from the roughness spectrum of the ground, validated against 1/10-scale 2-D rough surfaces, for use in time-domain codes. — [Faure et al. 2017](https://doi.org/10.1016/j.apacoust.2016.11.019)
- [ABS] Kayser, Gauvreau & Ecotière 2019 (PE with effective-impedance roughness, Morris screening over σ, roughness height, correlation length): "variations of ground roughness induce non-negligible differences in sound pressure levels regarding the ground absorption, even for high height sound source, i.e., nongrazing incidence." — [Kayser et al. 2019](https://doi.org/10.1121/1.5131652)
- [ABS] But Kayser et al. 2020, full emission + propagation chain, downwind: "the relative air humidity and the ground roughness parameters appear to be negligible on sound pressure level sensitivity", refraction and ground absorption dominating. — [Kayser et al. 2020](https://doi.org/10.1121/10.0002872)
- [FT] Nord2000 treats roughness differently: an rms height class per segment, N (0 m, heights within ±0.25 m), S (0.25 m, ±0.5 m), M (0.5 m, ±1 m), L; it feeds a coherence coefficient. "The method for including ground roughness in Nord2000 Road has not yet been validated, and it is recommended to use class N." — [Nord2000 Road User's Guide, Table 12](https://egra.cedex.es/EGRA-ingles/I-Documentacion/National_Methods/Users_Guide_Nord2000_Road.pdf)

### Inferences
- **[OWN] Magnitude with AT2026 eq. 16 on class-5 grass** (slit pore 250/0.25). ⟨H⟩ = 1 cm takes Re Z from 14.4 to 10.2 at 250 Hz and from 8.1 to 4.0 at 1 kHz; 3 cm takes it to 2.9 at 250 Hz and to zero from 500 Hz. Change in flush-receiver level, 1.5 km:

  | elev | ⟨H⟩ | 63 | 125 | 250 | 500 | 1000 Hz |
  | --- | --- | --- | --- | --- | --- | --- |
  | 1° | 1 cm | 0.0 | −0.8 | −1.0 | −1.5 | −2.5 |
  | 1° | 3 cm | +0.3 | −2.9 | −2.2 | −2.5 | −3.6 |
  | 2° | 1 cm | 0.0 | −0.3 | −0.5 | −1.0 | −2.0 |
  | 2° | 3 cm | 0.0 | −0.6 | −0.9 | −1.3 | −2.6 |
  | 5° | 1 cm | 0.0 | −0.1 | −0.2 | −0.5 | −1.2 |
  | 5° | 3 cm | +0.1 | 0.0 | +0.2 | +0.5 | −0.7 |

  At a 1.2 m pole the 500 Hz dip deepens by 3–9 dB at 1–2°.
- **Direction.** Small-scale roughness modeled this way makes the prediction quieter at 2° and below, by 0.5–2.5 dB in 250–1000 Hz. It therefore does not explain a model that is already too quiet: adding an effective-admittance roughness term would widen the gap, not close it. The sign is against small-scale roughness as the cause.
- The Kirchhoff (Ament) factor and the effective-admittance picture apply to different scales. At grazing the Rayleigh parameter k ⟨H⟩ sin ψ is tiny (0.003 at 250 Hz, 2°, 3 cm), so the Kirchhoff factor is unity to four digits while the effective-admittance shift is 1–2 dB. The two are not alternatives; the second is what is missing.
- Howe's statement that the roughness and impedance contributions are comparable for grass below 450 Hz applies to bosses with k0 a of order 0.1–0.5 (a = 5–25 cm at 100 Hz), i.e. tussocks, furrows and hummocks, not blade-scale texture.

### Gaps
- The boss-model constants (Tolstoy/Twersky/Howe) could not be transcribed reliably; the formulas above give structure only. *Predicting Outdoor Sound* ch. 5 has them and was not accessible.
- I could not confirm the exact typography of AT2026 eq. 16 (whether σ is flow resistivity; units of ⟨H⟩ in ν); the numbers above assume σ is flow resistivity in Pa s m⁻² and ⟨H⟩ in meters.
- Faure et al. 2015 (the Internoise/Euronoise paper named in the brief) was not found; the 2017 Applied Acoustics paper is its successor. No full text.
- No measurement of roughness-induced level change from an elevated source at km range was found.

---

## 6. Spatial and temporal variability of ground impedance

### Takeaway
Fitted flow resistivity over grass changes by a factor 1.2–1.8 between summer and winter, and by up to 7.5 with soil moisture in one study; the direction is always wetter = harder. By my calculation a factor of two in σ is 1–2 dB at 2° and 0.5–1.5 dB at 5°, so seasonal and moisture state is the same size as the discrepancy being chased.

### Cited Findings
- [FT] "Several studies have shown that the impedance of ground surfaces increases with moisture content. Consequently, the flow resistivity obtained by fitting level difference spectra tends to be higher in Winter than in Summer. Also, for some ground surfaces, the fitted effective layer depth is less in Winter than in Summer." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] AT2026 Table 5, wet/dry (winter/summer) ratios of mean fitted values. Flow resistivity: compacted grass 1.78, uncompacted grass 1.21, fine sand 1.1–3.5, coarse sand 5.0–9.8, ash forest 1.49, lime 3.34, cherry 1.40, oak 1.81. Layer thickness: compacted grass 0.38, uncompacted grass 0.97, ash 0.51, lime 0.69, cherry 1.03, oak 1.12. The grass values are from Guillaume et al. 2015; sand is 50% saturation against dry. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] Recommended handling: "estimates of attenuation due to ground effect should be made by using higher values of flow resistivity, lower values of porosity and lower values of layer depth, in the ranges appropriate to classes 4 to 7." The class ranges (e.g. 250 ± 100 for pasture) "allow for the potential variations due to differences in degree of compaction, moisture content and roughness." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [ABS] Cramond & Don 1987 (impulse method; grassland, barren earth, sand; up to 40% moisture by weight): effective flow resistivity of the top surface "may change from 200 cgs rayl at 10% moisture to 1500 cgs rayl at 35% moisture content. Results taken over grassland and sand follow the same trend if normalized to their saturation value." Wetting also layers the soil, with position-dependent resonances near 5 kHz. — [Cramond & Don 1987](https://doi.org/10.1121/1.395565)
- [ABS] Horoshenkov & Mohamed 2006: "a relatively small change (e.g., from 0% to 11% by volume) in the degree of water saturation can result in a much larger change (e.g., twofold) in the acoustic surface admittance"; neither the one- nor the two-parameter model fits partly saturated sand accurately. — [Horoshenkov & Mohamed 2006](https://doi.org/10.1121/1.2338288)
- [SNIP] Guillaume et al. 2015: two-microphone narrow-band fitting of flow resistivity and effective thickness on three grounds in winter and summer, analyzed for spatial and temporal variability and the consequence for sound levels. Numbers not accessible beyond those in AT2026 Table 5. — [Guillaume et al. 2015](https://doi.org/10.1016/j.apacoust.2015.01.024)
- [ABS] Dragna et al. 2015: the difference between admissible and one-parameter models at longer range is "significantly greater than those resulting from short range spatial variability and comparable with seasonal variability over grassland." — [Dragna et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26520321/)
- [ABS] Horoshenkov, Khan & Benkreira 2013: plants on soil behave as an equivalent fluid (Miki three-parameter) whose flow resistivity and tortuosity follow leaf area density and leaf angle; 50–1600 Hz impedance tube. — [Horoshenkov et al. 2013](https://doi.org/10.1121/1.4798671)
- [FT] Nguyen et al. 2022 (ray model, 80 m source, 1.5 m receiver, downwind, nine impedance models × four wind-profile models): impedance-model choice is more than 75% of total A-weighted model uncertainty out to 0.5 km, under 25% at 1 km, and the wind profile is over 95% beyond 2.5 km; the impedance share falls more slowly for C- and unweighted levels, "indicating that impedance models were likely to have a greater effect on the low frequency content." Total spread exceeds 10 dBA beyond 3.5 km. Going from σ = 12.5 to 200,000 kPa s m⁻² changes the *uncertainty* by only 1.7 dB. — [Nguyen et al. 2022](https://arxiv.org/pdf/2205.13695)
- [ABS] Wilson, Pettit, Ostashev & Vecherin 2014 treat uncertain ground and atmosphere by Latin hypercube sampling with importance sampling; first example is "the Lloyd's mirror effect for an elevated sound source in a homogeneous atmosphere". No dB figures in the abstract. — [Wilson et al. 2014](https://doi.org/10.1121/1.4890644)
- [FT] Parry's thesis (Monte Carlo on simple engineering models): with geometric uncertainty present, ground hardness "is found not be greatly influential", with little difference between 500 and 20,000 kPa s m⁻². This concerns near-ground geometries and statistics of excess attenuation, not an elevated source at grazing. — Parry 2020, PhD thesis, University of Sheffield, "Statistically assessing forward-wise and inverse uncertainties in outdoor acoustics" (local copy; no URL retrieved)
- [ABS] Soh et al. 2010: impedance measured directly from the experiment's own signals at several points along a 400 m path over an agricultural field, consistent with earlier measurements there. — [Soh et al. 2010](https://doi.org/10.1121/1.3501116)

### Inferences
- **[OWN]** Using the sensitivity from question 1 (flush receiver; +1.1 to +2.2 dB per doubling of σ at 2°, +0.5 to +1.5 dB at 5°, +1.5 to +3.0 dB at 1°): the compacted-grass winter/summer ratio of 1.78 is about 1–1.8 dB at 2°; Cramond & Don's factor 7.5 is about 4–6 dB at 2° and 2–4 dB at 5°. AT2026's class-5 range (150–350 kPa s m⁻², Ω 0.1–0.4) spans −4.1 to +1.6 dB at 250 Hz and −8.9 to −1.2 dB at 1 kHz at 2° (rows "150/0.4" and "350/0.1" in the table above), i.e. a 5–8 dB band from within-class uncertainty alone.
- The level is as sensitive to porosity as to flow resistivity in the slit-pore family, because Z scales as 1/Ω. A site fitted with Delany–Bazley carries no porosity information at all.
- Nguyen et al.'s finding that impedance stops mattering beyond 1 km rests on a downwind ray model with multiple ground reflections and a 1.5 m receiver; it does not transfer to a single-bounce, near-neutral, flush-microphone geometry, where my numbers show the ground term is range-independent and stays 1–8 dB wide.

### Gaps
- Guillaume et al. 2015's spatial-variability statistics (standard deviation of σ across a field, in kPa s m⁻² and in dB of level) were not obtained.
- No study was found giving the uncertainty of km-range levels from an *elevated* source due to impedance uncertainty; the figures above are my own propagation of published parameter ranges.
- Nothing on dry lake bed surfaces (crusting, desiccation cracks, moisture cycles).

---

## 7. In-situ measurement and whether short-range fits transfer

### Takeaway
The standard methods fit level differences between two microphones 1–1.75 m from a source; they are trustworthy from roughly 400 Hz upward and poor below 300–400 Hz because the level difference there barely depends on impedance. The 63–630 Hz bands where the discrepancy sits are therefore mostly extrapolation of a model, which is exactly where the choice of model matters most.

### Cited Findings
- [FT] "Short-range measurements ... are of the difference in sound level spectra between vertically separated microphones located at 1.75 m (NORDTEST) or 1 m (ANSI) from a point source." Level difference LD = EA(M1) − EA(M2) with Q from eq. 3. The NORDTEST data are third-octave; five repeats per site. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [ABS/web listing] ANSI/ASA S1.18-2018: real and imaginary impedance "from in-situ measurements of the sound pressure levels at two vertically separated microphones using specified geometries"; either fit an impedance model to templates or deduce impedance directly; "applicable to all nominally flat, commonly occurring surfaces including grassland or snow covered ground." — [ANSI/ASA S1.18-2018](https://webstore.ansi.org/standards/asa/ansiasas1182018)
- [FT] Kruse tested the ANSI S1.18 two-microphone method on several grounds: "the method already achieves reproducible and credible results, but only for frequencies above approx. 400 Hz and not for acoustically hard grounds." — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [FT] Why: "the level difference below approx. 400 Hz will not depend much on the flow resistivity of the ground ... Conversely, this means that minor errors in the level difference lead to sizeable errors in the deduced surface impedance." "Even small errors in the transfer function in the range of accuracy of high quality measurement systems cause considerable errors in the ground impedance", and microphone position errors are not negligible. — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [FT] With optimised geometries and switched microphone positions, "reliable measurements above 100 Hz are made possible"; "contrary to the recommendations of ANSI S1.18, the application of large loudspeakers has to be advised." — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [FT] "The frequency of the ground dip is most sensitive to changes in the surface impedance for near grazing incidence", hence heights of 0.2–0.4 m at 1.75 m. Kruse's geometry A is 24° grazing. — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [ABS] Taraldsen & Jonasson 2011: their field impedance measurements cover "200 Hz to 2.5 kHz", while the harmonised methods go "down to the 25 Hz third octave band". — [Taraldsen & Jonasson 2011](https://doi.org/10.1121/1.3500694)
- [FT] Kruse on extended reaction in the inversion: deducing material properties from an effective impedance "will require measurements at different angles of incidence." — [Kruse dissertation](http://www.akustik.uni-oldenburg.de/literatur/Kruse/Dissertation_Kruse_Ohne_Artikel.pdf)
- [ABS] Attenborough & Waters-Fuller 2000 used an elevated source and vertically separated receivers at 30 m over cultivated ground, i.e. a longer-baseline variant, and needed the rough-surface model to explain it. — [Attenborough & Waters-Fuller 2000](https://doi.org/10.1121/1.1288940)

### Inferences
- A two-microphone fit at 1–1.75 m samples a patch of ground under a meter across, at 15–30° grazing, above about 300–400 Hz. The flight-test geometry samples hundreds of meters (next question), at 2–10°, from 63 Hz. Three separate extrapolations are involved: in frequency (model-dependent; see the 2–4.5 dB spread between models at equal σ), in angle (roughness and extended reaction both grow toward grazing), and in area (spatial variability).
- The frequency extrapolation is the one with a known bias. A Delany–Bazley σ fitted to the 400 Hz–4 kHz dip reproduces the mid-frequency impedance and then, below that, gives Re Z too small relative to admissible models of the same ground — the 2015 abstract's point — which lowers the predicted grazing level in the low bands.
- Because the flush-microphone level at 2–5° is itself a clean function of β (the plane-wave formula holds to 0.1 dB), the flight-test or a stationary elevated-source measurement over the array is a better impedance probe at 63–630 Hz than any short-range fit. This is an inference from my numbers, not a published recommendation.

### Gaps
- ANSI S1.18 and NORDTEST ACOU 104 texts were not read. From memory only, and not to be cited: NORDTEST uses source 0.5 m, microphones 0.5 and 0.2 m, 1.75 m range, 200 Hz–2.5 kHz; ANSI geometries A/B use separations of 1.75 and 1.0 m with heights of 0.2–0.5 m. These need checking against the standards.
- No paper was found that tests explicitly whether short-range fitted parameters predict km-range grazing levels from an elevated source.

---

## 8. Impedance discontinuities, mixed ground and the Fresnel zone

### Takeaway
For a flush receiver the controlling ground is a long thin ellipse starting at the microphone and stretching toward the source: with a quarter-wavelength criterion, about 900, 410 and 130 m at 63, 250 and 1000 Hz at 2° elevation, shrinking to 85, 22 and 6 m at 10°, and only a few meters wide. It is never "the ground at the microphone" at low frequency and low elevation.

### Cited Findings
- [FT] Harmonoise: "A Fresnel weight for a ground segment is basically equal to the fraction of the Fresnel ellipse that is covered by the segment"; the ellipse is the intersection with the ground of a Fresnel ellipsoid around the image source and receiver. — [Salomons et al., DAGA 2009](https://pub.dega-akustik.de/NAG_DAGA_2009/data/articles/000250.pdf)
- [FT] HARMONOISE weights excess attenuation by Fresnel-zone area following Hothersall & Harriott: A_ground = Σ w_i EA_i (eq. 24). Boulanger et al. showed better agreement with laboratory data when the weighting is on the pressure ratio before the logarithm, A_ground = −10 log(Σ w_i P_i) (eq. 26); AT2026 propose this, and it agrees better with 2-D BEM in their strip cases. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] "The HARMONOISE scheme uses a wavelength fraction of 1/16 when defining the Fresnel ellipsoid ... even though reference [26] on which it is based recommends a wavelength fraction of 1/3 or 1/4. Fig. 12 shows that a wavelength fraction of 1/3 results in spectra somewhat closer to those predicted by BEM, so use of a wavelength fraction of 1/3 is proposed." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] "A Fresnel zone method cannot predict the spectral oscillations associated with diffraction effects", but it is more satisfactory than the CNOSSOS G-factor method even after band averaging and turbulence. — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [FT] AT2026 test cases (road source 0.01–0.75 m, receiver 1.5 m at 50 m, soft strip 7–10 m or 7–20 m): the two weightings differ most when the specular point is in the middle of the strip (75.8 against 72.6 dB LAeq for one source, BEM 75.5); overall A-weighted differences are small because road noise peaks at 1–2 kHz, and "the difference would be more pronounced for a source that its main energy in the 500 Hz band." — [AT2026](https://doi.org/10.1016/j.apacoust.2025.111067)
- [ABS] Hothersall & Harriott 1995: compared a boundary integral solution, a de Jong-type semi-empirical diffraction formula and a Fresnel-zone method for discontinuities perpendicular and parallel to the path and for strips. "The de Jong approach proved flexible in application to all the boundary conditions considered and agreed well with the boundary integral equation and the experimental results except in one case of propagation over a strip when the grazing angle was very small." — [Hothersall & Harriott 1995](https://doi.org/10.1121/1.412136)
- [FT] de Jong, Moerkerken & van der Toorn 1983: a diffraction term at the admittance discontinuity (road surface to grass verge) added to the two homogeneous solutions; including it "increases the agreement between measured and calculated curves considerably" at 7.5–15 m; grassland flow resistance assessed by fitting, 40–200 kN s m⁻⁴ in their cases. — [de Jong et al. 1983](https://doi.org/10.1016/0022-460X(83)90941-0) (J. Sound Vib. 86, 23–46; read from the local copy in ~/Desktop/papers)
- [ABS] Lam & Monazzam 2006 extend the de Jong semi-empirical formulation to multiple impedance discontinuities. — [Lam & Monazzam 2006](https://doi.org/10.1121/1.2216905)
- [FT] Nord2000: "The ground effect ... is calculated for each type of ground to be found inside the Fresnel-zone and the resulting ground effect is calculated as a weighted average taking into account the fraction of the Fresnel-zone covered by each type of ground surface." Validation: standard deviation 1 dB to 400 m, 2 dB at 600–1000 m over flat ground against PE. — [Kragh, Acoustics 2011](https://acoustics.asn.au/conference_proceedings/AAS2011/papers/p110.pdf)
- [SNIP] Nord2000 defines the ellipsoid by |SP| + |PR| − |SR| = F_λ λ with foci at source (image) and receiver, and uses a rectangular zone of 4/π times the ellipse area. — [Nord2000 comprehensive model report, Plovsing](https://forcetechnology.com/-/media/force-technology-media/pdf-files/projects/nord2000/nord2000-comprehensive-outdoor-sound-propagation.pdf)

### Inferences
- **[OWN] Size of the zone.** Ground points P with |S'P| + |PR| − |S'R| ≤ Fλ, S' the image source, flat ground, range 1.5 km, c = 343 m/s. Extent from the receiver toward the source / maximum half-width, in meters:

  | receiver | F | elev | 63 Hz | 250 Hz | 1000 Hz |
  | --- | --- | --- | --- | --- | --- |
  | flush | 1/4 | 2° | 898 / 24.7 | 410 / 8.4 | 129 / 2.3 |
  | flush | 1/4 | 5° | 289 / 14.0 | 85 / 3.8 | 22 / 1.0 |
  | flush | 1/4 | 10° | 85 / 7.6 | 22 / 2.0 | 5.6 / 0.5 |
  | flush | 1/16 | 2° | 407 / 8.3 | 129 / 2.3 | 34 / 0.6 |
  | flush | 1/16 | 5° | 84 / 3.8 | 22 / 1.0 | 5.6 / 0.2 |
  | flush | 1/16 | 10° | 22 / 1.9 | 5.6 / 0.5 | 1.4 / 0.1 |
  | 1.2 m pole | 1/4 | 2° | 909 / 25.2 | 444 / 9.2 | 179 / 3.2 |
  | 1.2 m pole | 1/4 | 5° | 307 / 15.0 | 108 / 4.8 | 45 / 1.8 |
  | 1.2 m pole | 1/4 | 10° | 96 / 8.7 | 34 / 2.9 | 16 / 1.2 |

  Behind a flush receiver the zone extends only Fλ/(1 + cos ψ) ≈ Fλ/2 (0.7 m at 63 Hz for F = 1/4). For the pole the specular point is 33.6, 13.6 and 6.8 m in front of the mast at 2°, 5° and 10°, and at higher frequency the zone no longer reaches the mast.
- **[OWN] Closed form** for a flush receiver when the zone is short against the range: length toward the source ≈ Fλ/(1 − cos ψ) ≈ 2Fλ/ψ² (ψ in radians). It over-estimates once the length approaches the range (2.2 km against the exact 0.9 km for 63 Hz at 2°).
- Consequences for a test array: (a) a ground board or plate, a graded pad, a road or a taxiway within tens to hundreds of meters up-path of a microphone is inside the zone for the low bands at low elevation; (b) each frequency band and each elevation averages a different patch, so a single impedance per microphone cannot be right for all bands unless the site is uniform over about a kilometer; (c) the area weighting is on pressure squared per AT2026, so a hard patch covering a fraction w of the zone raises the level toward the hard-ground value roughly in proportion to w, which is the direction of the observed deficit.
- Hothersall & Harriott's one failure of the de Jong approach was a strip at very small grazing angle, so at 2° neither simple method is validated; the choice of F (1/16 to 1/3) changes the zone length by a factor of up to 5.

### Gaps
- The Nord2000 and Harmonoise technical reports were not obtained, so their exact F_λ values (Nord2000 is reported elsewhere to use a frequency-dependent fraction) and transverse-size rules are not quoted. The only F values cited are those AT2026 attributes to HARMONOISE (1/16) and to Hothersall & Harriott (1/3 or 1/4).
- Rasmussen's work on impedance discontinuities was not retrieved.
- No published Fresnel-zone sizes for an airborne source and ground receiver were found; the table is mine.
