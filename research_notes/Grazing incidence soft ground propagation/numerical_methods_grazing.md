# Numerical methods for near-grazing propagation over soft ground (reference and production-by-table)

Provenance key used throughout:
**[FT]** = I read the full text. **[AB]** = abstract only. **[SN]** = search-engine snippet or secondary account only. **[README]/[LICENSE]** = repository file read.
"Own estimate" = my arithmetic or reasoning, not a published number.

Access note: AIP, ScienceDirect, Taylor & Francis, DTIC and PubMed web pages returned 403 or cookie walls. Full texts read were Nyborg et al. 2023 (DTU Orbit), Colas et al. 2023 (arXiv), Cotté 2019 (author copy) and the ncpaprop LICENSE. Everything else is abstract-level (via NCBI E-utilities) or snippet-level. The Consensus search quota was exhausted, so coverage of Salomons' GFPE guidance, Gilbert & Di's starters, FFP cost, Rosenbaum's HPM thresholds, moment-equation PEs and the Kayser kriging accuracy numbers is thin; those are listed under Gaps rather than filled from memory.

---

## 1. PE variants: angle limit, range step, impedance boundary, accuracy, cost

### Takeaway
In the wind-turbine literature (the closest analogue: 60–240 m sources, 2 m receivers, 0.2–3 km, grass) a Crank–Nicolson Padé(1,1) PE at λ/8–λ/10 is the de facto reference and agrees with a linearized-Euler solver to 0.5 dB beyond 1 km for a near-ground receiver. The two published ways to go faster are split-step Padé(2,2) at Δx = 2λ (Cotté) and Fourier split-step (GFPE/BTPE) at Δr = 2.5–10λ (Nyborg), and the Fourier route shows visible interference-pattern shifts over hard ground at those steps.

### Cited Findings

**Crank–Nicolson wide-angle PE (CNPE / GTPE / WAPE)**
- [FT] DTU's GTPE (CN wide-angle, Sack & West 1995 terrain form, Fortran/MPI): "limiting elevation angle of 30°" attributed to the second-order Gaussian starter; grid from a convergence study Δz = Δr = λ/8 at all frequencies; absorbing layer 50λ thick; computational height 500 m, with 1 km and 2 km tested and showing no effect; effective sound speed. Because of the angle limit with a 100 m source, "ground distances shorter than 200 m are considered invalid" (that is, elevation above about 27°). — [Nyborg et al. 2023, JASA 154:1299](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] ECL/ENSTA WAPE (two variants, second-order FD in z, Crank–Nicolson, Thomas algorithm, Salomons second-order starter, PML top, impedance ground folded into the tridiagonal coefficients): converged at Δx = Δz = λ/10 on flat ground; a 31° hill needed λ/50. Cost: 1 CPU hour for all 35 frequencies (50 Hz–1 kHz, 3 km × 300 m) on flat ground, 30 CPU hours for the hill case. — [Colas et al. 2023, JASA 154:1413 (arXiv)](https://arxiv.org/pdf/2309.09669)
- [FT] Accuracy against a linearized-Euler (LEE) reference, flat ground, variable-porosity ground with σ = 50 kN s m⁻⁴, source 58–142 m, receiver 2 m: "The difference between the ΔL obtained with the three methods does not exceed 0.5 dB for x > 1 km" at the three frequencies plotted (up to 1 kHz). Frequency-averaged error over the whole (x, z) plane: above 3 dB close to the source (PE angle limit), below 0.5 dB farther out in a still atmosphere. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] Effective-sound-speed error: with a 10 m/s hub-height ABL, the effective-sound-speed WAPE's frequency-averaged error rises above 3 dB across much of the domain, but "close to the ground the error is still relatively small"; the vector WAPE (Ostashev et al. 2020 Sec. VI.A) shows no such growth. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [README] Univ. Gustave Eiffel's public WAPE implements Ostashev et al. 2020 with Padé(1,1) and Crank–Nicolson, recommends Δx = Δz = λ/10, and warns of instabilities with oversized steps or too thin an absorbing layer. — [Universite-Gustave-Eiffel/WAPE](https://github.com/Universite-Gustave-Eiffel/WAPE)
- [FT] Field validation of a CNPE with an elevated source over grass (hub 119 m, microphones at 1.7 m, 178/535/845 m, three incoherent point sources): overall level within 2 dB outside upwind; third-octave mean absolute differences 1.6–3.2 dB; more than 20 dB off at 845 m in strong upwind because turbulence was not modeled. — [Könecke et al. 2023, Wind Energ. Sci. 8:639](https://wes.copernicus.org/articles/8/639/2023/)

**Split-step Padé (Collins 1993) in the atmosphere**
- [FT] Cotté uses split-step Padé (N,N) after Collins (1993) and Dallois et al.; "accurate results are obtained in a typical wind turbine configuration with N = 2 and mesh sizes Δx = 2λ and Δz = λ/10", ground impedance condition at z = 0, absorbing layer on top, effective sound speed. It is stated (citing Cotté 2018, JSV 422:343) to be "more computationally effective than a classical wide-angle PE based on a Padé (1,1) approximation". "The angular validity increases with the order N of the development and depends on the mesh size Δx chosen (Dallois et al., 2002)." — [Cotté 2019, JASA 145:1363](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [FT] Cost with that scheme: "a set of PE calculations for the 49 frequencies takes 8 min to run on one core" (Xeon X5650, 2.66 GHz) for one source height and one azimuth; results are used out to about 1.2 km. — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)

**Fourier split-step: GFPE and Beilis–Tappert (BTPE)**
- [FT] BTPE (split-step Fourier with a terrain phase screen, after Parakkal et al. 2012; "many similarities to the GFPE"): Δz = λ/10 throughout; Δr = 2.5λ at and below 50 Hz, 5λ for 63–100 Hz, 10λ at 125 Hz and above; fourth-degree polynomial starter; GFPE-style absorbing layer; domain 500 m to 2 km range and 1000 m to 5 km. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] BTPE step-size convergence (σ = 2000 kPa s/m², Gaussian hill, source 100 m, receiver 2 m): "average convergence within ±0.5 dB between Δr = 0.1λ and 2.5λ for 50 Hz and between Δr = 2.5λ and 10λ for 500 Hz", but the interference pattern is displaced "around 50 m" for 10λ versus 0.1λ at 50 Hz and "around 10 m" at 500 Hz. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] GTPE versus BTPE versus JEPE, flat ground, source 100 m, receiver 2 m: at 50 Hz "a slight offset of 1 dB but generally agree well"; the offset grows with flow resistivity; at 500 Hz in downward refraction the interference minima are phase-shifted by an amount that grows with flow resistivity. Over terrain with σ = 2000 kPa s/m² the pointwise difference reaches "up to 20 dB at 50 Hz" at the dips. The authors attribute this to grid spacing and to dispersion error at large elevation angles, and note the disagreement is "diffused" when several incoherent sources are used. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Upward refraction: the PE models "start to deviate when entering the shadow zone"; the authors call the no-turbulence shadow levels "unrealistically low" and the differences there numerical noise. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [AB] GFPE surface-wave pole: "The original formulation of the Green's function parabolic equation (GFPE) can have numerical accuracy problems for large normalized surface impedances"; the improved form with the surface-wave pole subtracted "is shown to be accurate for surface impedances varying over 2 orders of magnitude, with the largest having a magnitude exceeding 1000", and is slightly faster because the surface-wave term is no longer computed separately. — [Gilbert 2015, JASA 137:EL8](https://pubmed.ncbi.nlm.nih.gov/25618104/)
- [SN] Cooper & Swanson give GFPE parameter guidelines by comparing with analytical solutions; their acceptance criterion is "within 3 dB of the analytical solution", and they note PE use is "generally limited to experts because of their dependence on numerous algorithm parameters". — [Cooper & Swanson 2007, Appl. Acoust. 68:390](https://www.sciencedirect.com/science/article/abs/pii/S0003682X06000624)
- [AB] Rotated GFPE over a valley slope, road source: predictions within 3 dBA of measurements to 1 km in windless conditions, with "typical parameters of this method ... optimized for computational efficiency". — [Van Renterghem, Botteldooren & Lercher 2007, JASA 121:2522](https://pubmed.ncbi.nlm.nih.gov/17550151/)

**Moving-medium PEs without the effective-sound-speed approximation**
- [AB] "The narrow-angle parabolic equation (NAPE) with the effective sound speed approximation (ESSA) ... is valid only for angles less than 20°"; the low-Mach equations and "the WAPE with the ESSA are applicable only under rather restrictive conditions on the medium velocity"; an "effective numerical implementation of the WAPE for arbitrary Mach numbers in the Padé (1,1) approximation" is given. — [Ostashev, Wilson & Muhlestein 2020, JASA 147:3969](https://pubmed.ncbi.nlm.nih.gov/32611146/)
- [AB] A split-step *spectral* algorithm for an extra-wide-angle PE "accounts for the propagation angles up to 90°", valid for small refractive-index variations, and is "particularly useful" in moving media because derivative terms become wave vectors. — [Ostashev, Muhlestein & Wilson 2019, JASA 145:1031](https://pubmed.ncbi.nlm.nih.gov/30823829/) (an erratum exists: [JASA 153:2944](https://pubmed.ncbi.nlm.nih.gov/37195626/))
- [AB] Phase-preserving WAPE/NAPE: "omitting the second-order terms ... can result in large cumulative phase errors for long propagation ranges"; implementation needs "minimal modifications to existing codes" and is more accurate than the ESSA versions. — [Ostashev, Colas, Dragna & Wilson 2024, JASA 155:1086](https://pubmed.ncbi.nlm.nih.gov/38341733/)

**Other**
- [FT] JEPE (Karasalo & Sundström 1996): vertical discretization to a differential-algebraic system, fourth-order strongly damped marching, adaptive range step from a local error estimate, boundary-fitted curvilinear grid. In Nyborg's flat-ground case it does not show the near-field interferences the other two do, which the authors attribute to the BTPE and GTPE starting conditions. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)

### Inferences
- The in-house CN Padé(1,1) at 20 points per wavelength in height is finer than any grid quoted above (λ/8 to λ/10). The published evidence says its remaining error at grazing is not the propagator: Colas's 0.5 dB agreement with LEE is at λ/10.
- For the 0–10° regime the angle limit is not binding for any variant. What differs between variants is the ground condition and the range step. Finite-difference schemes put the impedance into the matrix and are robust for soft ground; Fourier schemes handle it through a reflection-coefficient integral plus a surface-wave term, which is where the GFPE's soft-ground difficulty lives (Gilbert 2015 exists because of it).
- The effective-sound-speed error matters aloft and at range, and is small near the ground (Colas). Under reciprocity the *source* becomes the elevated point, so ESSA error at 30–300 m height is relevant to a reciprocal table even though the receiver is on the ground. Ostashev 2024's phase-preserving form is described as a minimal change to an existing CN code and is the cheap fix.
- Own estimate of cost per frequency for a 5 km range and 500 m useful height, using NICEOPS's own figures (2.7 ns per point per step, Δz = λ/20, Δr = λ/10, top margin max(250 m, 100λ)):

  | f | height points | range steps | CN Padé(1,1), one core |
  | --- | --- | --- | --- |
  | 100 Hz | 4.9e3 | 1.5e4 | 0.2 s |
  | 1 kHz | 4.4e4 | 1.5e5 | 18 s |
  | 4 kHz | 1.8e5 | 5.9e5 | 280 s |

  Cost goes as f². With Cotté's split-step Padé(2,2) at Δx = 2λ there are 20 times fewer steps and two tridiagonal solves per step, so roughly 10 times cheaper at equal Δz (own estimate; Cotté gives no ratio, only "more computationally effective"). Cotté's own 8 min for 49 frequencies on a 2010 Xeon core is consistent with that order.
- A Fourier split-step at Δr = 10λ takes 100 times fewer steps than the CN march, but each step is an FFT pair over the column and Nyborg's 10–50 m interference shifts show the accuracy price. That matters less for a band- and time-integrated SEL than for a single-frequency comparison, but it disqualifies it as the reference.

### Gaps
- Salomons' own GFPE step-size and near-ground guidance (Salomons 1998 JASA 104:100; *Computational Atmospheric Acoustics* 2001) was not retrieved. My unverified recollection is that the GFPE range step is limited by the refraction phase factor and by receiver height near the ground, with steps of order 5–50λ quoted for high receivers and much smaller steps needed for low receivers over soft ground in strong gradients. Treat as unverified.
- Cooper & Swanson's actual recommended parameter values were behind a paywall; only their 3 dB criterion is confirmed.
- Collins 1993, Sack & West 1995 and Dallois et al. 2001/2002 were not read directly; angle limits for Padé(n,n) as a function of n and Δx are therefore unquoted. Kayser's WAPE thesis was not retrieved.
- No published per-frequency run times for a 5 km × 500 m atmospheric domain were found for any variant; the table above is my arithmetic.
- The generalized terrain PE's known sensitivity to ground-derivative errors is only mentioned in passing by Colas (citing Sack & West); no quantitative accuracy for flat soft ground specifically.

---

## 2. Reciprocity for an elevated moving source

### Takeaway
Every elevated-moving-source PE study I could read runs the PE *forward* from a handful of fictive source heights and interpolates to the nearest; none uses a reciprocal run from the receiver. The theoretical basis for the reciprocal run (flow reversal) is established, and the forward studies supply the numbers needed to judge it: ΔH of 10–15 m between source heights for under 1 dB error, and one to five frequencies per third-octave band.

### Cited Findings
- [FT] Cotté's moving-monopole (MM) model: the propagation term is ΔL(ω, x) relative to free field, computed by PE for Nh fictive source heights between Hmin = 35 m and Hmax = 125 m, and applied to each blade segment by nearest-neighbor lookup (maximum height error ΔH/2). The free-field source level with its directivity multiplies ΔL. — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [FT] Error versus number of source heights, maximum over 500–1200 m at 2 m height, against a 19-height reference (ΔH = 5 m). OASPL: 3 heights (ΔH = 45 m) 0.4–3.6 dB depending on azimuth; 7 heights (ΔH = 15 m) 0.3–1.2 dB; 10 heights (ΔH = 10 m) 0.3–0.8 dB. The worst azimuths are crosswind-to-upwind (τ = 110°) where a shadow boundary is present; downwind is 0.3–0.4 dB even with 3 heights. A single point source is off by up to 11.0 dB upwind because it puts the shadow "too close". — [Cotté 2019, Table I](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [FT] Factorizing directivity out of the PE is validated there: the MM model (monopole ΔL times directional free-field level) and the Amiet-PE model (directivity built into the starter by back-propagation) are both "in excellent agreement" with the analytical solution in a homogeneous atmosphere with seven or more heights, and MM is about 24 times cheaper (ratio M·Nβ/Nh at Nh = 10). Amiet-PE cost: more than 30 h per azimuth; MM: 24 min (Nh = 3) to 2.5 h (Nh = 19). — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [FT] The back-propagation starter: a known field (direct plus image with spherical-wave reflection coefficient Q) on a vertical line at x_is is marched backwards to x = 0 in a homogeneous atmosphere over rigid ground, then marched forward in the real medium; used "in order to preserve the angular validity of the split-step Padé (2,2) method". — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [FT] Colas et al. adopt the same fictive-source scheme: "seven fictive source heights are sufficient ... for the flat case"; "for the case with a hill, 30 source heights are required" because "changing the source height by a few meters can affect the position of the focusing zone at the ground of more than 100 m". — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] Frequency sampling per third-octave band in Colas et al., Table I: one frequency per band for 50–160 Hz; two at 200 and 250 Hz; three at 315 Hz; four at 400, 500 and 630 Hz; five at 800 Hz and 1 kHz — 35 frequencies for 50 Hz–1 kHz. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] Cotté uses 49 frequencies for his third-octave spectra. Nyborg et al. and Könecke et al. compute at band center frequencies only. — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf); [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Single-frequency, single-source results are fragile at interference dips: model-to-model differences of up to 20 dB at 50 Hz over hard ground, and "the agreement of the spectra obtained by the PE models might vary significantly with the chosen receiver positions"; distributed incoherent sources "average out the interference effects". — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Azimuth: Cotté builds c_eff(z) = c(z) + U(z) cos τ per propagation direction and steps τ by 10° for directivity plots. Kayser et al. (as summarized by Nyborg) found "the angle between the direction to the receiver and the wind direction as well as the ground flow resistivity have significant influences". — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf); [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [SN] Flow reversal theorem: symmetry "with respect to interchange of the source and receiver positions and the simultaneous reversal of the ambient flow ... an extension for moving fluid of the reciprocity principle valid in quiescent media". — [Godin 1997, Wave Motion 25:143](https://www.sciencedirect.com/science/article/abs/pii/S0165212596000376)
- [FT] Doppler/convection on the source enters Colas's LEE normalization as k_eff = k0/(1 + M) at the source height — a reminder that a moving-medium source term is not the still-air one. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)

### Inferences
- **Why reciprocity pays here and not for wind turbines.** A turbine's sources span 90 m of height and there are few receivers, so 7–10 forward runs are cheap. NICEOPS has sources over 30–300 m (270 m of height: 19–28 forward runs at Cotté's ΔH = 10–15 m) but only two receiver heights (0 and 1.2 m). Over flat, range-independent ground and a stratified atmosphere the field depends only on (range, z_s, z_r, azimuth to the wind, f), so one reciprocal run per (z_r, azimuth, f) yields ΔL at every source height and range at full grid resolution, with no source-height interpolation error at all. The number of receivers (tens or 15,000) and of trajectory states does not enter the run count.
- **Under the effective-sound-speed approximation reciprocity is exact and needs no theorem**: the medium is motionless with c_eff(z) for that azimuth, so swap source and receiver and keep the *same* c_eff profile. That is identical to reversing the wind and propagating the other way. With a vector-wind PE it is Godin's flow reversal and the wind must be reversed explicitly.
- **Receiver exactly on the ground becomes source exactly on the ground.** A Gaussian starter centerd at z = 0 has half its width below the boundary, and its image term is weighted by a plane-wave reflection coefficient that is poorly defined at zero height over soft ground. Two ways round it that need no literature: (a) validate the z = 0 starter against Weyl–van der Pol with the source on the ground (the in-house tests already hold the PE to it for an elevated source); (b) use Cotté's back-propagation idea — build the exact homogeneous-air field of a ground source on a vertical line a few tens of meters out from the spherical-wave reflection formula and start the march there, which also removes the starter's angle limit.
- **Flush plate microphones** are not a point at z = 0; the existing plate table (K per band and azimuth) would multiply the reciprocal bare-ground field in the same way it does for the shadow transition.
- **Directivity coupling**: at 0–10° the direct and ground-reflected paths leave the source within about twice the elevation angle of each other, so applying the sphere level at the eigenray launch direction times a monopole ΔL is the MM factorization Cotté validated. It fails only where the sphere's directivity changes by decibels over a degree or two.
- **Table resolution** (own estimate from Cotté's Table I): nearest-neighbor at ΔH = 10 m costs up to about 1 dB near a shadow boundary; a reciprocal table stored at, say, 2–5 m in height with linear interpolation in dB should be well under that. Azimuth every 10–20° in U cos τ; interpolating in cos τ rather than τ follows the physics.
- **Frequency sampling**: Colas's one to five frequencies per band is the only explicit published rule found. For a moving source the SEL time integral sweeps the interference pattern past the receiver, which is a second averaging the wind-turbine papers do not have; center-frequency-only may be tolerable for SEL but not for a band time history.
- Own estimate of reciprocal table cost with the existing CN PE (5 km, 300 m of source height plus 250 m margin, Colas's frequency sampling extended to five per band above 1 kHz): about 3.5 core-minutes per (azimuth, receiver height) for bands to 1 kHz and about 100 core-minutes to 5 kHz. For 18 azimuths and 2 receiver heights: about 2 core-hours to 1 kHz, 60 core-hours to 5 kHz. Bands above 1–2 kHz are better left to the ray model, since absorption removes them at range and the f² cost is all there.

### Gaps
- I found no paper that uses a reciprocal PE (source at the ground receiver, field read aloft) for an elevated moving source. Cotté 2019, Colas 2023 and Nyborg 2023 (all read in full) run forward. Barlas et al. 2017 and Kayser's papers were not read in full, so I cannot rule it out there. The absence is a finding of this search, not a proof.
- Gilbert & Di's exact (Laplace-transform) starter and Salomons' starter with an impedance image term were not retrieved; nothing quotable on the z = 0 source over impedance ground.
- No published dB errors for interpolating PE tables across range, frequency or azimuth; only source height (Cotté) is quantified.
- Godin 1997 was seen as a snippet only; whether flow reversal holds exactly for the specific one-way wide-angle moving-medium PEs (rather than the full wave equation) was not checked.

---

## 3. Alternatives at grazing: FFP, normal modes, hybrid stitching

### Takeaway
The only hybrid built for aircraft noise is the FAA/Volpe Hybrid Propagation Model (PE + FFP + straight ray), which exists because of the PE's elevation limit. Normal modes are established for nocturnal ducting over impedance ground, where a lossy surface mode and a "quiet height" near 2 m are physical features any low-receiver model should reproduce. I could not retrieve quantitative cost or accuracy for either.

### Cited Findings
- [SN] "HPM is a composite of three propagation methods: a parabolic equation (PE) model; a fast field program (FFP); and a straight ray-trace model", developed with Penn State for FAA's AEDT; documented in Volpe reports DOT-VNTSC-FAA-12-05 (Assessment of the Hybrid Propagation Model, vols. 1–2). — [Volpe/ROSAP report listing](https://rosap.ntl.bts.gov/view/dot/9629/dot_9629_DS1.pdf)
- [AB] The hybrid PE–FFP model "overcomes the elevation angle limitations associated with using a pure parabolic equation approach", treats uneven terrain, mixed impedance and range-dependent meteorology, with "emphasis on low frequency noise" and only "preliminary results" on turning vertical-plane fields into contour maps. — [Rosenbaum et al., Penn State record](https://pure.psu.edu/en/publications/enhanced-sound-propagation-modeling-of-aviation-noise-using-a-hyb/)
- [AB] A highway study likewise couples "a parabolic equation model ... to an analytical Green's function solution close to the source field". — [Ovenden, Shaffer & Fernando 2009, JASA 126:25](https://pubmed.ncbi.nlm.nih.gov/19603859/)
- [SN] Waxler's modal model treats "propagation of sound over an impedance ground plane in a stratified atmosphere which is downward refracting near the ground but upward refracting at high altitudes", with a ground impedance having "both real and imaginary parts so that the ground is lossy as well as compliant". — [Waxler 2004, JASA 115:1437](https://pubs.aip.org/asa/jasa/article-abstract/115/4/1437/546564/Modal-expansions-for-sound-propagation-in-the)
- [FT, lay-language paper] Nocturnal duct experiment, 100–500 Hz sweeps, 1.8 km, 4 m tower: a surface mode that is the slowest and most attenuated; a "quiet height" at about 2 m with levels "about 10 dB" below those on the ground, appearing "about a half a kilometer from the source" once the surface mode has died; attenuation "increases dramatically with increasing frequency". — [Waxler, Talmadge et al., ASA 152nd meeting](https://acoustics.org/pressroom/httpdocs/152nd/wax.html)
- [README] ncpaprop's Modess is a normal-mode expansion for a single tone in a stratified atmosphere under the effective-sound-speed approximation; the package targets "narrowband propagation from 0.1 to 10 Hz". — [ncpaprop-release](https://github.com/chetzer-ncpa/ncpaprop-release)

### Inferences
- Waxler's quiet height bears directly on NICEOPS's two receiver heights: in a surface inversion a 1.2 m microphone and a flush one can differ by around 10 dB beyond 0.5 km at a few hundred hertz. A ground-effect model built on a single reflected ray cannot produce that; the PE and the modes can.
- A mode sum is naturally "reciprocal": p is a sum of ψ_n(z_s) ψ_n(z_r) H0(k_n r), so one eigenproblem per (f, azimuth) gives every source height, both receiver heights and all ranges. Own estimate of trapped-mode count for a duct of depth H and speed excess Δc: N is about (2 H f / c)·sqrt(2Δc/c); for H = 100 m, Δc = 5 m/s that is about 10 modes at 100 Hz and about 100 at 1 kHz. At low frequency in an inversion that is cheaper than a PE march and exact in angle. The in-house experience with the Attenborough case 2 series (steep modes beyond the PE's reach, surface-wave pole, continuation failures past |q| = 6) says the root-finding, not the summation, is the cost.
- Modes do not help upwind (the residue series converges only in the shadow) or in neutral conditions; they are a special-purpose exact answer for low-frequency ducting, not a general table generator.
- An FFP is range-independent by construction, which matches the flat-ground table assumption, and has no angle limit, so it is the natural partner above 30° as HPM uses it. At 0–10° it offers nothing the PE does not, and turbulence cannot be added to it.

### Gaps
- The HPM switching angles, cost and validation numbers: the Volpe PDFs would not download (HTML landing pages only) and the Penn State abstract gives none. Boeker & Rosenbaum's "Intelligent switching between different noise propagation algorithms" report exists on ROSAP but was not read.
- No FFP cost or accuracy figures were retrieved (Raspet et al. 1985, Lee et al. 1986, Wilson 1993, CERL-FFP, FFLAGS). Nothing citable on FFP pole handling over soft ground at grazing.
- Raspet, Baird & Wu 1992 and Waxler's later impulse papers were not read; mode counts and attenuation rates above are my estimate and a lay summary respectively.
- Whether ncpaprop's Modess/ePape accept a finite ground impedance was not confirmed from the README; the manual PDF was not read.

---

## 4. Time-domain and other full-wave methods

### Takeaway
Linearized-Euler finite-difference solvers with a time-domain impedance condition and a moving frame are a working 2D reference to 3 km and 1 kHz, at about 1200 CPU hours per run against 1 CPU hour for the PE, and they confirm the PE to 0.5 dB at a 2 m receiver. They are a reference for validating the PE, not a production or table-building tool. No open implementation of that class was found.

### Cited Findings
- [FT] LEE solver (ECL): 2D conservative LEE on a curvilinear mesh, optimised fourth-order 11-point centerd differences (Bogey & Bailly), fourth-order six-step Runge–Kutta (Berland), selective filtering, PML top. Ground: broadband time-domain impedance condition (Rienstra 2006; Troian et al. 2017) with the impedance as a rational function of frequency so the convolution becomes auxiliary ODEs; "not implemented on the impedance but on the reflection coefficient for stability reasons". Ground model: variable porosity, 50 kN s m⁻⁴, 100 m⁻¹. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] Moving frame: a 300 m × 300 m window follows the pulse through a 3000 m × 300 m domain, which "reduces the computational cost by one order of magnitude". For accuracy to 1 kHz and 3 km: grid 0.05 m (λ/6.8 at 1 kHz), CFL 0.5, 127,000 time steps, 36 × 10⁶ points, "approximately 1200 ... CPU hours, which is equivalent to 3 days on a 16-core machine" per source height. One run gives all frequencies. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] A side effect of the moving frame near the source: the reflected wave cannot reach the top of the window, giving "a zone without any interferences and a value of ΔL close to 0" for x < 100 m. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] The Gaussian pulse source itself costs accuracy at the top of the band: the difference between the pulse's spectrum and an ideal monopole "reaches around 1.5 dB at 1 kHz" in their appendix check. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [AB] First demonstration of the class: FDTD LEE with a time-domain impedance boundary condition, compared with frequency-domain PE in upward and downward refraction. — [Cotté & Blanc-Benon 2009, JASA 125:EL202](https://pubmed.ncbi.nlm.nih.gov/19425623/); the search snippet for the companion work gives the range as 500 m ([Cotté et al., TDBC paper](https://www.researchgate.net/publication/215514751_Time-Domain_Impedance_Boundary_Conditions_for_Simulations_of_Outdoor_Sound_Propagation)).
- [Extracted by an automated reader from the arXiv HTML; treat with care] The AIAA 2022 precursor (hill case) reports a 0.5 m grid, CFL 0.8, "over 130000 iterations" for 5 km, "a 1 km propagation simulation takes around 240 CPU hours, and a full simulation takes around 20000 CPU hours". The same extraction gave a maximum frequency of "11 kHz", which is not possible on a 0.5 m grid and is probably a misread; I did not verify these figures in the text. — [Colas et al. 2022, AIAA 2022-2923 (arXiv)](https://arxiv.org/html/2207.13406v1)
- [AB] LEE solvers have also been evaluated for strong-turbulence propagation (intensity-variance saturation). — [Ehrhardt, Cheinet, Juvé & Blanc-Benon 2013, JASA 133:1922](https://pubmed.ncbi.nlm.nih.gov/23556562/)
- [SN] openPSTD (pseudo-spectral time domain, 2D, GUI-driven) is GPL. — [openPSTD](https://github.com/openPSTD/openPSTD)
- [SN] A review of efficient FDTD for outdoor sound exists (hybrid FDTD–PE, moving frames). — [Van Renterghem 2014, Int. J. Aeroacoustics 13:385](https://journals.sagepub.com/doi/10.1260/1475-472X.13.5-6.385)

### Inferences
- The LEE cost scales with f³ in 2D for fixed range (points per area times time steps), against f² for the PE. At "a few hundred Hz" and 5 km a 2D LEE is tens of CPU hours (own estimate: Colas's 1200 h scaled by (300/1000)³ × (5/3) is about 54 h) — feasible as a one-off reference, particularly for a vector-wind check of the effective-sound-speed error aloft that the reciprocal table would inherit.
- What the LEE adds over the existing references is narrow: no angle limit, a true moving medium, and back-scatter. For flat ground at 0–10° the in-house exact solutions (Weyl–van der Pol, residue series, mode series) are stronger references because they are exact, and Colas shows the LEE agreeing with the PE rather than correcting it.
- Two of the LEE's own error sources (the pulse's 1.5 dB spectral error at the band top; the impedance fitted as a rational function) are comparable to or larger than the PE's error at grazing, so it is not automatically the more accurate of the two.

### Gaps
- No open LEE/FDTD code with a time-domain impedance condition and moving frame was found. The ECL solver is not public as far as I could see. openPSTD is GPL and oriented to urban scenes at short range; I did not confirm whether it supports refraction or frequency-dependent impedance.
- Ostashev & Wilson's time-domain boundary conditions for soil models (relaxation/Zwikker–Kosten), Dragna's long-range papers, discontinuous Galerkin and BEM-with-refraction were not retrieved; nothing citable on their km-range practicality.
- SPARROW, DTU WindSTAR-Pro, and transmission-line-matrix (TLM) codes from Univ. Gustave Eiffel: not found or not checked.

---

## 5. Turbulence in PE at grazing

### Takeaway
Still-air PE in an upwind shadow is wrong by tens of decibels against measurement, so a table meant for upwind use needs turbulence or an explicit floor. The published practice found here is brute-force phase screens with tens of realizations (64 in Hart et al.); the moment-screen method propagates the second moment directly and removes the realizations, at the cost of a matrix-valued march.

### Cited Findings
- [FT] CNPE without turbulence versus measured wind-turbine noise: "exceeds 20 dB at 845 m" in strong upwind conditions, while other directions agree to 1.6–3.2 dB per band. — [Könecke et al. 2023](https://wes.copernicus.org/articles/8/639/2023/)
- [FT] Nyborg et al. leave turbulence out deliberately because it "increases the computational time and complexity significantly", and state that shadow-zone transmission losses without it "are unrealistically low". — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Cotté 2019 also omits turbulence and notes that the level "might be underestimated when a shadow zone is present". — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [AB] Hart et al.: CNPE datasets of 5000 Latin-hypercube samples, one with "single realizations of turbulence", one with "average TL fields for 64 realizations of turbulence"; against an 8 km propagation experiment, "predictions obtained through a single realization of turbulence agree better with the experimental observations". — [Hart, Wilson, Pettit & Nykaza 2021, JASA 149:4384](https://pubmed.ncbi.nlm.nih.gov/34241466/)
- [SN] Moment-screen method: "direct numerical solution of a parabolic equation (PE) for the second moment of the sound field in a refracting medium with random scattering ... determines the mean-square sound pressure without requiring generation of random realizations"; "the second-moment matrix is factored into components that are independently propagated with a conventional PE algorithm". — [Wilson, Ostashev & Lewis 2009, Waves Random Complex Media 19:369](https://www.tandfonline.com/doi/abs/10.1080/17455030802647462)
- [README] The public WAPE code exposes turbulence through a single strength parameter. — [Universite-Gustave-Eiffel/WAPE](https://github.com/Universite-Gustave-Eiffel/WAPE)

### Inferences
- A reciprocal table multiplies cleanly by realizations: own estimate, 64 realizations of the to-1-kHz table is about 130 core-hours for 18 azimuths and 2 heights, versus 2 core-hours still. Because the quantity wanted is a band-averaged, time-integrated energy, each frequency within a band can be its own realization (as the in-house `solve_pe_band` already does), so the effective sample count is frequencies-per-band times realizations and far fewer than 64 marches per frequency should be needed.
- Turbulence only matters to the table where the still-air field is well below the scattered field, that is, in the upwind shadow and at interference minima. An alternative to realizations is to build the still-air table and add the existing `ShadowScattering` fit in energy, which is what the code does now; the table would replace the ray/transition model, not the scattering term.
- Reciprocity under turbulence holds realization by realization for a frozen medium, and the ensemble mean-square is symmetric if the turbulence statistics are horizontally homogeneous; height-dependent turbulence strength is unaffected by the swap.
- Hart's finding that single realizations match observations better than 64-realization means is about comparing with *individual* measurements; for a predicted SEL footprint the ensemble mean is the right target.

### Gaps
- Gilbert, Raspet & Di 1990; Chevret, Blanc-Benon & Juvé 1996; Cheinet et al. 2012; Salomons on screen spacing, realization counts and height-dependent strength: not retrieved. No citable number for "realizations needed for a band-averaged level to ±x dB".
- The moment-screen method's cost (rank of the factorization, how it scales with frequency and range) was not available from the abstract. "Cheng et al." was not found.
- No published comparison of mutual-coherence-function PE against phase screens at grazing over soft ground was found.

---

## 6. Surrogates and emulators of PE excess attenuation

### Takeaway
Published PE surrogates are built for uncertainty quantification over meteorological and ground parameters, not as drop-in propagation kernels, and the only accuracy figure I could confirm is coarse: 5–7 dB median absolute error against field data. Nothing found supports replacing a per-run PE table with a learned emulator at the sub-decibel level NICEOPS works to.

### Cited Findings
- [AB] Pettit & Wilson: single-frequency point source; atmosphere by Monin–Obukhov similarity, ground by a relaxation model, propagation by PE; parameters as independent random variables sampled by Latin hypercube; proper orthogonal decomposition of the field plus cluster-weighted models gives "compact representations of the sound pressure random field" and full-field sensitivities. — [Pettit & Wilson 2007, JASA 122:1374](https://pubmed.ncbi.nlm.nih.gov/17927400/)
- [AB] Hart et al.: CNPE training data, 5000 Latin-hypercube samples per dataset; ensemble decision trees, neural networks and cluster-weighted models; "regression predictions have 5-7 dB in median absolute error" against an experiment out to 8 km; "surrogate data quality depends on an accurate characterization of refractive and scattering conditions". — [Hart et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34241466/)
- [AB] Kayser et al.: quasi-Monte Carlo sampling of environmental parameters feeding an Amiet emission model coupled to a PE; "a metamodel of the global (emission-propagation) wind turbine noise model was built using the kriging interpolation technique to drastically reduce calculation time"; downwind conditions only. — [Kayser, Gauvreau, Écotière & Mallet 2022, JASA 151:390](https://pubmed.ncbi.nlm.nih.gov/35105016/)
- [AB] The Morris screening behind it: "parameters involved in atmospheric refraction and in ground absorption have a significant influence"; "relative air humidity and the ground roughness parameters appear to be negligible". — [Kayser, Cotté, Écotière & Gauvreau 2020, JASA 148:3623](https://pubmed.ncbi.nlm.nih.gov/33379890/)
- [AB] For ground parameters alone (flow resistivity, roughness height, correlation length), roughness variations "induce non-negligible differences ... even for high height sound source". — [Kayser, Gauvreau & Écotière 2019, JASA 146:3222](https://pubmed.ncbi.nlm.nih.gov/31795674/)
- [SN] The kriging metamodel is deployed as the online tool WindTUNE. — [WindTUNE paper](https://www.researchgate.net/publication/371957559_WindTUNE_a_new_tool_for_modeling_wind_farm_noise_uncertainties)
- [SN] Operator-learning surrogates of 2D PE solvers exist in underwater acoustics. — [DeepONet PE surrogate, JASA 154:3276](https://pubs.aip.org/asa/jasa/article-abstract/154/5/3276/2922093/Training-a-deep-operator-network-as-a-surrogate)

### Inferences
- These surrogates take a low-dimensional parametrization of the atmosphere (MOST: friction velocity, heat flux or Obukhov length, roughness) plus flow resistivity. NICEOPS takes a measured profile file, which is not in that parameter space; a surrogate would have to be trained on profile shapes, or the profile projected onto MOST parameters with a loss the validation work already shows is large (the log-law shadow boundary finding).
- Pettit & Wilson's POD idea is the useful part for a table: the ΔL(r, z_s) field over azimuth and frequency is highly compressible, so a handful of POD modes could interpolate across azimuth and frequency more safely than pointwise interpolation through interference patterns. That is an inference, not something they report for this use.
- Given a reciprocal table costs core-hours per run at worst, a per-run PE is affordable and carries no training-distribution risk; a surrogate is only attractive if many runs share an atmosphere class.

### Gaps
- Kayser 2022's metamodel accuracy (kriging error in dB, number of training simulations), and any Kayser 2023 PCE work: not in the abstract and full text not accessible.
- Hart et al. 2016 and the surrogate-versus-CNPE (rather than versus-experiment) error: DTIC copy blocked. The 5–7 dB figure mixes surrogate error with model-input error and should not be read as emulator accuracy.
- "Pinto" metamodels: not found.
- Training costs (CPU hours for the 5000-sample CNPE sets) not available.

---

## 7. Open-source implementations and licenses

### Takeaway
Two usable external cross-checks exist: Univ. Gustave Eiffel's MATLAB WAPE (EUPL-1.2, vector-wind Padé(1,1) with impedance ground and turbulence) and ncpaprop (permissive, but built for infrasound at 0.1–10 Hz). Neither can be copied into an MIT codebase without care; both can be run alongside.

### Cited Findings
- [README] WAPE: MATLAB, no toolboxes; implements Ostashev et al. 2020 (no effective-sound-speed assumption, arbitrary Mach number); ground impedance from airflow resistivity and layer thickness; turbulence strength parameter; Padé(1,1) Crank–Nicolson; license "EUPL-1.2". — [Universite-Gustave-Eiffel/WAPE](https://github.com/Universite-Gustave-Eiffel/WAPE)
- [LICENSE] ncpaprop 2.1.0, © 2016 University of Mississippi: permission "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies", conditions being retention of the notice, prominent notice of modifications, no endorsement, and passing on a copy of the license; some included code under a Triad/LANL permissive license. — [ncpaprop LICENSE](https://raw.githubusercontent.com/chetzer-ncpa/ncpaprop-release/master/LICENSE)
- [README] ncpaprop: C/C++ with FFTW, GSL, PETSc and SLEPc; modules include Modess, WMod, ModBB and ePape; aimed at "narrowband propagation from 0.1 to 10 Hz and broadband propagation from 0.05 Hz to 2 or 3 Hz". — [ncpaprop-release](https://github.com/chetzer-ncpa/ncpaprop-release)
- [SN] openPSTD: GPL. — [openPSTD](https://github.com/openPSTD/openPSTD)
- [FT] DTU's GTPE is Fortran/MPI and KTH's JEPE is in-house; neither is described as public. Nyborg's BTPE is a MATLAB re-implementation of Parakkal et al. 2012. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)

### Inferences
- EUPL-1.2 is a copyleft license; treat WAPE like GPL code — run externally as a benchmark, do not port its source. Its value is specific: it is the one public code that does the vector-wind WAPE, which is the check the effective-sound-speed reciprocal table most needs.
- ncpaprop's license is close to NCSA/BSD with a modification-notice clause and would be compatible with MIT redistribution under its conditions, but its frequency regime and ground treatment make it a weak cross-check for 20 Hz–5 kHz over grass.

### Gaps
- pyram/RAM (Collins' split-step Padé, underwater), PyPE/pywaveprop, SPARROW, WindSTAR-Pro: not checked; licenses unknown.
- WAPE's README states no validation; I did not run it.

---

## 8. Benchmark cases beyond Attenborough et al. 1995

### Takeaway
The best new benchmark for this regime is a measured one with public data: Könecke et al. 2023 (119 m source, grass, 1.7 m microphones to 845 m, met mast to 100 m, CC BY). Nyborg et al. 2023 supplies a clean numerical inter-model case with a 100 m source and three ground classes, and Colas et al. 2023 a LEE-verified one.

### Cited Findings
- [FT] Attenborough et al. 1995 and the later numerical studies built on it "are limited to a source height of 5–10 m, which is naturally not applicable to wind turbines". — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Nyborg numerical cases: source 100 m, receiver 2 m, ranges to 5 km; c_eff(z) = c0 + b ln(z/z0 + 1) with c0 = 340 m/s, z0 = 0.1 m, b = +1, 0, −1; ground by Attenborough's four-parameter model with σ = 12.5, 250, 2000 kPa s/m² (pore shape factor 0.75, grain shape factor 0.5, porosity 0.3); frequencies 50 and 500 Hz versus range, and third-octave centers 20 Hz–2 kHz at 500 and 1500 m; flat, four Gaussian hills, undulating terrain; three independent PEs (GTPE, BTPE, JEPE). Impedances were computed once and fed to all models. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Nyborg loudspeaker campaign: loudspeaker on a nacelle at 109 m, seven microphones from 158 m to 978 m, measured flow resistivity 250 kPa s/m² to 270 m and 500 beyond, met mast (temperature gradient 3–105 m, power-law wind), 2 h of 2-min spectra, 100 Hz–2 kHz. Models depart from measurement above 800 Hz, attributed to loudspeaker directivity; differencing two microphones removes it. — [Nyborg et al. 2023](https://backend.orbit.dtu.dk/ws/files/334925908/1299_1_10.0020562.pdf)
- [FT] Könecke et al.: multi-MW turbine, hub 119 m, rotor 114 m, flat grassland in northern Germany; microphones at 1.70 m at 178, 535 and 845 m; 100 m mast (29, 57, 76, 100 m), 10-min averages; data "provided and available for research purposes" at doi:10.25835/0012136, article under CC BY 4.0. — [Könecke et al. 2023](https://wes.copernicus.org/articles/8/639/2023/)
- [FT] Colas et al. flat case: hub 100 m (sources 58–142 m), neutral LES boundary layer with u* = 0.512 m/s (10 m/s at hub), wake ducting producing a focus on the ground at 1.7 km, variable-porosity ground 50 kN s m⁻⁴ / 100 m⁻¹, 50 Hz–1 kHz, to 3 km, receiver 2 m; LEE and two PEs agree to 0.5 dB beyond 1 km. — [Colas et al. 2023](https://arxiv.org/pdf/2309.09669)
- [FT] Cotté 2019 case: hub 80 m, log profile with 8 m/s at hub, receiver 2 m, 0.2–1.2 km, all azimuths; upwind shadow from about 400 m. — [Cotté 2019](https://perso.ensta-paris.fr/~cotte/files/Cotte_jasa145_2019%20Extended%20source%20models%20for%20wind%20turbine%20noise%20propagation.pdf)
- [SN] Harmonoise reference model: a combination of PE, BEM and ray methods for "geometrically relatively simple but physically complex" road and rail situations, with a statistical atmosphere. — [Defrance et al. 2007, Acta Acustica 93:213 (TNO record)](https://repository.tno.nl/SingleDoc?docId=22334)
- [SN] IEA Wind Task 39 initiated a wind-turbine noise code benchmark in its first phase, proceeding "from wind turbine noise generation to airborne noise propagation over large distances". — [IEA Wind Task 39](https://iea-wind.org/task39/)

### Inferences
- The elevation angles in these cases: 100 m source at 0.5–5 km is 11° down to 1.1°; Könecke's three stations are 34°, 12.5° and 8°. Only the 535 and 845 m stations are inside the 0–10° regime, and only just. The Nyborg numerical case to 5 km is the one that reaches the 1–5° range NICEOPS cares about.
- Nyborg's log profile with b = ±1 is the same family as the in-house `example_windy` finding (shadow boundary fixed by c(z_s)/c(0)); running `pe_compare` on it with σ = 250 kPa s/m² and reading the published BTPE/GTPE/JEPE curves is a direct external check of both the PE and the per-band boundary law, at 50 and 500 Hz.
- All the wind-turbine cases use receivers at 1.7–2 m. None has a flush microphone, so the on-ground reciprocal starter has no external benchmark; the in-house Weyl–van der Pol and residue-series tests remain the only reference for z_r = 0.

### Gaps
- IEA Task 39 propagation-benchmark case definitions and results were not found in citable form. Harmonoise benchmark cases (heights, ranges, impedances) were not read; they target low road/rail sources and are probably outside the elevated-source regime.
- Whether Nyborg's loudspeaker data are public was not stated in the paper.
- No NATO or ISO benchmark with an elevated source over soft ground was found.
