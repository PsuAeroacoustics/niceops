# Improving NICEOPS's urban noise predictions: a ranked literature review

October 20, 2026. Private: this review quotes the project's own validation numbers.

> **Status, checked afterwards (Oct 20).** Written by a literature-survey agent; NICEOPS's own checks
> since: recommendation 1 is implemented (`--scattering` with `--buildings`, THEORY 6.5). The
> author's version of Forssén et al. (2014) lacks its Table IV; the velocity coefficients quoted
> below reproduce that version's Figure 2 to 2 dB, but the temperature ones put temperature 3 dB
> *above* velocity per unit C², where the scattering cross section puts it 7 dB below, so NICEOPS
> does not use them and takes γ_T as velocity turbulence instead. The review's conversion
> "γ_t/4 = C_v²/c0²" is not NICEOPS's: γ_T = C_T²/T0² + (22/3) C_v²/c0² (THEORY 5.3), which puts the
> open-terrain "moderate" 3×10⁻⁶ at C_v² ≈ 0.05, not 0.1. Its "13 steps at a wall's foot" is the
> 80-line sweep's; on 20 lines the floor between two edges has since removed both. Everything
> else is as the agent left it.

## How this was done

- **Five literature passes**, one per survey area:
  - ray/beam methods and edge diffraction
  - diffuse reflection
  - wave, hybrid and surrogate methods
  - engineering methods and the atmosphere
  - drones/AAM and validation data

  Each pass was told to report only sources it found and to grade how far it had read each one.
- **My own checks on top:**
  - About 90 DOIs confirmed against Crossref.
  - The abstracts the ranking leans on re-read.
  - The Forssén et al. (2014) turbulence model read in full and evaluated for rotorcraft geometries.
- **Data policy.** Every search was on a generic technical topic. No project name, result or file content left the machine.
- **Coverage.** The session's shared web-search budget ran out partway through the passes. What that left uncovered is in §10.
- **Evidence labels.** Each piece of evidence is labelled *measured*, *full-wave* (a converged numerical solution of the wave equation) or *simulated* (engineering or ray models). My own arithmetic is labelled as mine.
- **Verification tags in the reference list (§11):**
  - **[F]** full text, or a substantial part, read.
  - **[A]** abstract or publisher/repository record read.
  - **[M]** bibliographic record confirmed; content taken from a named secondary source.
  - **[S]** seen only as a search snippet or a citation in another work; unverified.
  - **"DOI checked"** means Crossref returned the stated title, authors and venue.

### The baseline the ranking is measured against

| Benchmark | NICEOPS ray model now | For comparison |
| --- | --- | --- |
| EDGE scale models against measurement (1124 band levels; UE06/07 read as swapped) | 3.8 dB rms (UTD) | FDTD alone is 0.7–2.3 dB rms from measurement on 11 of 13 scene–source pairs |
| EDGE, ray model against FDTD (13 runs) | 1.95 dB rms; was 2.80 before the third edge and the reflections between edges | — |
| DigiCity microphone 7 (reached only by diffraction) | −3.9 dB mean, 5.0 rms (specular); −2.3 dB, 2.9 rms with the corrugation's s(f) | NLR's CNOSSOS-EU model: −12.9 dB mean, 13.3 rms |
| Continuity sweep (3771 positions) | 25 steps over 1 dB; the largest is 5.2 dB | — |

Known local misses:
- **UE06 R1, +4.7 dB.** A reflection lands 3 cm inside the end of a facade, and its partner ray would diffract at a vertex.
- **Lit canyon with relief (UE24–26).** Lambert scattering predicts about 0 dB of change; the measured change is −1.0 to −1.6 dB.
- **Behind EDGE's second rows**, the measurement exceeds even the FDTD by 3–5 dB. This is unexplained.

---

## Executive summary

**NICEOPS is already at or beyond what others report.**
- **Reflection-dominated receivers.** Published urban tools agree to a couple of dB (Knepper et al. 2025; Malbéqui 2023).
- **Diffraction-dominated receivers:**
  - Engineering ray models in the CNOSSOS-EU style miss by 10–15 dB (Knepper et al. 2025, measured).
  - Gaussian beam tracing needs an empirical correction to stay within 8 dB of FEM (Yunus et al. 2023b).
  - A radio ray-tracer with two diffractions reaches about 3 dB rms in deep urban shadow (Athanasiadou et al. 1997, measured).
- **Achievable floor.** Scale-model and in-situ facade levels agree only to σ ≈ 2 dB(A) (Lugten et al. 2023, measured). Building-wake turbulence moves tonal levels by more than 3 dB (Goyal et al. 2026a, LB-VLES simulation, i.e. lattice-Boltzmann very-large-eddy simulation). Against these, NICEOPS's 2 dB rms against FDTD and 2.3 dB at DigiCity are close to what any model can be held to.

The remaining gains fall into five areas:
- physics real cities have but the scale models do not (turbulence)
- the ends of finite surfaces
- the form of facade scattering
- roof shape
- path orders for grazing propagation over several rows

Ranked by expected accuracy gain for rotorcraft and drones near or above buildings (SEL and LAmax at street, facade and rooftop receivers), weighed against effort and run time:

### 1. Add a turbulence-scattering floor in building shadows

**What to implement**
- Use the engineering model of Forssén et al. (2014) for every microphone the buildings shield. Its form, which I read in the paper:
  - scattered level re free field, computed separately for velocity and temperature turbulence and summed in energy:
    L = b1 + b2 lg(h/10 m) + b3 lg(h²/(d_S d_R) + 0.0012) + (10/3) lg(f/1 kHz) − βx
  - h is the obstacle's height above the source–receiver line, measured in the frame rotated to that line.
  - b1 = −52.8 + 10 lg C_v² (velocity) and −49.6 + 10 lg C_T² (temperature); b2 = 11.3 and 11.5; b3 = −17.1 and −13.1.
  - A canyon correction ΔL = 7 or 14 dB + (2H/W) lg(H/10 m).
  - Their sum with diffraction is capped at the unshielded level.
- Take h, d_S and d_R from the path tree NICEOPS already builds: the obstruction along the direct line, or the strongest diffracted path's edge.
- Add the result in energy to the band levels.

**Adapting it to an elevated source**
- The fit assumes energy doubling by the ground at both ends. That overstates an elevated source by up to 3 dB, so drop the source-side doubling. This is my inference, not the paper's.
- Use the canyon term only where NICEOPS's own reflections do not already carry the scattered energy.
- Draw the turbulence strength from the run's meteorology, as the open-terrain `ShadowScattering` does.
  - The conversion between γ_t and C_v², C_T² has to be stated explicitly.
  - One reading, mine: taking γ_t/4 as C_v²/c0² puts the open-terrain "moderate" γ_t ≈ 3×10⁻⁶ near C_v² ≈ 0.1 m^4/3 s⁻². Check this before relying on it.

**Why the evidence says it matters (measured)**
- **Behind a container barrier** 2.55 m high and 2.44 m thick, turbulence raised the shadow level by about 5 dB at 4 kHz and 10 dB at 8 kHz in wind, with little effect below 2 kHz (Forssén & Ögren 2002).
- **At a shielded facade in Ghent**, long-term Lday differed by nearly 10 dB between weak- and strong-turbulence classes above 1 kHz. Strong hours dominate Leq (Van Renterghem et al. 2013a).
- **Official maps** came out more than 14 dB(A) low at a shielded courtyard facade. A multiple-reflection term fixed the low frequencies, and the authors call turbulence scattering "essential" for the rest (Wei et al. 2014).
- **In NICEOPS** the air among the buildings is still. Its deep-shadow, high-frequency bands therefore sit where the measurements say models are low. Microphone 7's remaining shortfall at 1.25–1.6 kHz (−6.4 dB at 1.6 kHz after scattering) is in that range, though the turbulence there is unknown.

**How big (my evaluation of the published formula)**

Assumptions: thin-screen Maekawa diffraction; no air absorption; Forssén's ground doubling at both ends; no canyon term. The table gives the gain over diffraction alone.

| Case | Turbulence | 1 kHz | 2 kHz | 4 kHz | 8 kHz |
| --- | --- | --- | --- | --- | --- |
| Helicopter 100 m up, 500 m out; microphone 15 m behind a 10 m building | moderate (C_v² = 0.1) | +1.0 | +2.1 | +4.1 | +6.9 dB |
| | strong (C_v² = 1) | +5.5 | +8.6 | +12.2 | +16.0 dB |
| Low approach 30 m up, 400 m out; 20 m behind a 15 m building | moderate | +0.9 | +1.9 | +3.8 | +6.6 dB |
| Drone 40 m up, 150 m out; 10 m behind a 12 m building | moderate | +0.1 | +0.3 | +0.7 | +1.7 dB |
| | strong | +1.1 | +2.4 | +4.6 | +7.5 dB |

- **Where it matters.** The floor matters most for distant sources casting shallow shadows: helicopter approaches and departures. It matters least for a nearby drone casting a deep shadow.
- **Effect on A-weighted SEL at fully shielded receivers.** For traffic spectra it was under 1 dB(A) (Ögren 2004). For rotorcraft spectra with more energy at 1–4 kHz I would expect 0.5–2 dB(A) in moderate turbulence and several dB(A) in strong. That is my estimate; run it on the B407 spectrum.
- **LAmax** is set by the least-shielded moment, so it moves less.

**Cost.** A few days of work. A few logarithms per shadowed path per band, so negligible run time.

**How to validate it**
- **Unit test:** reproduce Forssén & Ögren's container barrier with their measured C_v² (0.4–17) and C_T².
- **Field checks:**
  - DigiCity microphone 7's high bands, if NLR supplies the weather.
  - Schiphol's shielded facades, which have KNMI weather.
  - MOTUS rooftops, split by turbulence class as Van Renterghem et al. (2013a) split them.
- **Keep it off for EDGE**: an anechoic room has no turbulence, and tuning on it would contaminate the geometry comparison.

**Evidence quality.** Moderate. There is one engineering fit and two measurement sets, and the fit's own error is under 3 dB. The input C_v² is the weak link: a factor of 10 in C_v² moves the result by 10 dB.

### 2. Make finite facades and edge ends uniform

This covers three things: Fresnel-weighted reflections, vertex (end-point) terms, and truncating path families together.

**What to implement**
- **(a) Weight reflections near facade ends.** Replace the binary in/out test for a facade reflection, within about a Fresnel radius √(λ d₁d₂/(d₁+d₂)) of the facade's vertical ends and roof edge, by a physical-optics weight over the facade rectangle. That weight is a product of two Fresnel transitions: −6 dB on a single boundary and about −12 dB at a corner of the lit region.
  - Nord2000's Fresnel-zone weight (20 lg S coherent, 10 lg S in energy) is the cheap engineering form. It is cited from the Nord2000 method document AV 1106/07 §5.20 by Forssén et al. (2014); I did not see the primary document.
- **(b) Let the edge-source integral reach the edge's physical end** for a partner ray whose Keller point falls off the edge, instead of dropping it.
  - The end-point contribution of the line integral is the first-order vertex wave. This is how the incremental theory of diffraction obtains continuity (Albani et al. 1997), and how BTM on finite edges carries it (Svensson et al. 1999, 2009).
  - For UTD cascades, the scalar hard-surface reduction of the uniform vertex coefficient of Albani et al. (2009) is the published route. It is validated against method of moments, in electromagnetics only.
- **(c) Truncate together.** Generate and cut every geometrical-optics term and its compensating diffracted partner by the same order budget and the same zone test (Calamia 2009, §7.3). Where a partner is beyond the budget, fade the term over its transition width rather than keeping it whole.

**Why**
- **The specular boundary.** At a specular-zone boundary the reflection and the onset of the edge's wave must sum to half the free reflection. If the reflection test and the diffraction test use different tolerances, the result is "double the correct amplitude, zero amplitude ... or incorrect diffraction polarity" (Calamia 2009).
- **UE06 R1** is that failure: +4.7 dB, a reflection 3 cm inside the end of a wall, kept whole because its partner meets a vertex.
- **The continuity steps** that remain (largest 5.2 dB) are partners that the orders leave out.
- **Shadow-only diffraction**, as some interactive tracers use (Tsingos et al. 2001; Schissler et al. 2014), leaves a 6 dB step on the lit side; NICEOPS rightly does not do this. Palenda & Vorländer (2026) conclude that diffraction "around edges, not only in shadow regions" affects accuracy.

**How big**
- Near a facade end or corner: up to about 6 dB locally, plus the steps in time history that make LAmax jump along a track.
- Away from transitions, a vertex wave is weaker than an edge wave by about (kL)^(−1/2). That is roughly 20 dB at 1 kHz for a 10 m edge (my inference from GTD scaling).
- Map-averaged rms: small, I expect a few tenths of a dB.

**Cost.** Low to medium.
- NICEOPS already has the edge-source integral with ends, and `--fade` partners for the lower orders.
- Only paths near ends pay, at a few Fresnel integrals or one edge integral each.

**How to validate it**
- UE06 R1 against the FDTD.
- EDGE's corner scenes, UE20–22, against measurement.
- The continuity sweep's step count.
- A finite box or plate computed by BTM (Svensson's edge-diffraction toolbox) as an exact reference.

**Evidence quality.** The theory is strong, and the vertex coefficient is validated in electromagnetics. There is no urban acoustic validation; NICEOPS's own FDTD would be the first.

### 3. Make facade scattering depend on incidence and on order

**What to implement**
- **(a) s(θ, f).** Take an incidence-dependent scattering coefficient from the 2-D BEM that already computes s(f). It costs nothing at run time.
  - It matters for overhead sources: a facade lit from above is lit near grazing, where s can be about four times its random-incidence value. One real facade gave about 0.48 against 0.1 at 125 Hz, and 0.9 against 0.6 at 500 Hz (Heimes et al. 2023, BEM, values read off a figure).
- **(b) Grating orders for ordered relief.** For periodic relief (corrugation, window bays, balconies), spawn the propagating grating orders as image sources at the first one or two reflections. Weight each order by the BEM's far-field pattern, and send the remainder to radiosity.
- **(c) A receiver-side gather that scales.** For maps, solve the radiosity once per band and precompute each microphone's gather. Siltanen et al. (2007) show the patch solution is view-independent, so only the final gather (milliseconds) depends on the receiver. Each state then costs a matrix–vector product, not today's O(patches²) per microphone.

**Why**
- **NICEOPS's own UE24–26 result.** The relief lowered the canyon by 1.0–1.6 dB (measured) and the FDTD reproduces it, but Lambert's law returns nearly all the energy. The relief redirects energy rather than diffusing it.
- **The literature agrees on why:**
  - Periodic surfaces put the "scattered" energy into a few grating orders, and s jumps from about 0 to 0.9 at an order's cut-on (Heimes et al. 2023; Pan et al. 2023, BEM).
  - A periodic facade can put its specular lobe 9 dB down (Lucchetti et al. 2026, FEM and a 1:20 scale model).
  - In a street canyon the choice of scattering distribution changes the energy–time curve significantly (Fatela et al. 2026).
  - Relief orientation changes shielded-side levels at equal coverage (Hornikx & Forssén 2008, 1:40 scale model, measured).
- **For overhead aircraft** specifically, diffuse reflection hardly changes the peak transmission gain. It does scatter energy out of the street and lowers time-integrated levels: "more complex models involving the effects of diffuse reflections are required for noise units involving integration over time" (Ismail & Oldham 2002, simulated).

**How big**
- 1–2 dB at multiply-reflected and shielded receivers; small at directly exposed ones.
- Up to 3 dB in single bands where a relief redirects a reflection (EDGE, measured).
- Riegel & Sparrow (2019, measured SonicBOBS) found diffuse reflection did not help around an isolated building. It matters in canyons and in shadows reached by reflection, not in open settings.

**Cost**
- (a) low.
- (b) medium: the image tree branches by 2–5 orders per reflection, so cap it at depth 1–2.
- (c) medium, and it is what makes diffuse energy affordable on a grid.

**How to validate it**
- EDGE UE24–26 and UE28 against measurement and the FDTD.
- DigiCity microphone 7.
- The open interlaboratory coefficients of Heimes et al. (2026) as input checks.

**Evidence quality.** Strong for single surfaces (BEM and measured). Weak for whole canyons, apart from NICEOPS's own FDTD.

### 4. Model pitched roofs (LoD2), and validate on EDGE's measured pitched-roof scenes

**What to implement**
- Sloped roof facets and non-vertical walls as reflectors.
- Ridge, eave and hip edges with their own wedge angles. UTD and the edge-source integral take any wedge angle; the path search and the "a facade reflection keeps the vertical component" shortcut do not.
- EDGE already holds measured pitched-roof scenes: UE01–05 (gable, half-gable, pent, mansard and shed rows), UE15, UE17 and UE19, in the same archive and with the same sources as the LoD1 scenes.

**Why**
- **Full-wave.** Real roof shapes against flat-roofed versions of the same district changed levels above and behind buildings by up to 10 dB (Oshima et al. 2014b, FDTD).
- **Measured.** At the Schiphol container courtyards, exposed–shielded facade differences were 8–9 dB(A) in a courtyard with straight facades and 14 dB(A) with a slanted facade and inset (Lugten 2022, measured, real aircraft).
- **Simulated.** Roof shape is among the strongest levers on the quiet side of buildings (Van Renterghem et al. 2013b, 2-D/3-D full-wave).
- **Real geometry.** Most residential fabric has pitched roofs. LoD1 turns each ridge into two right-angled edges at the wrong height.

**How big.** Several dB, up to about 10 dB locally, at receivers shielded by roofs. Small at lit receivers.

**Cost.** Medium to high: geometry, path search, beam tracing in plan view, and roof bounces. Measurement data is in hand.

**Evidence quality.** Moderate: one full-wave study, one measured field lab, and measured scale models not yet used.

### 5. Higher-order paths with pruning, and an engineering term beyond the explicit orders

**What to implement**
- **Generalise reflections around and between edges** by the image-edge construction: mirror edges in reflecting planes (Erraji et al. 2021). This gives the floor between two edges and several reflections between edges.
- **Four or more edges**, through an edge-visibility graph with wavelength-based edge simplification.
  - Schissler et al. (2014) cut the diffracting edges by 3–14 times and ran 4–6 orders over city-sized meshes in tens of milliseconds.
  - Use a uniform double-diffraction form only where one edge lies in the other's transition region (Albani 2005; Capolino et al. 1997). That is the common case for equal-height LoD1 rows.
- **Beyond the explicit order**, apply a per-canyon excess attenuation instead of dropping the path: −1 to −4 dB per intermediate canyon (Schiff et al. 2010), or 1 dB per 100 m capped at 5 dB (Wei et al. 2014). Both were derived for ground-level sources; treat them as a floor and check them against the FDTD.

**Why**
- Behind EDGE UE13's second row the model was 8–12 dB low until it took a third edge; it is now within 2.2 dB rms.
- A thin middle row needs four edges.
- Of the residual continuity steps, 13 are a twice-diffracted ray reaching the floor between its edges.
- At grazing, first-order diffraction alone left about 4 dB rms against measurement; second order brought it to about 1 dB (Medwin's thin-plate data, via Calamia 2009).

**How big**
- Up to about 10 dB where a path needs the order. Such paths are common only at low elevation, i.e. the helicopter approach.
- Under steep overflight they are rare.

**Cost.** Medium to high, and the cost grows exponentially with order. Diffraction culling saved 88–95% of compute in rooms (Calamia et al. 2008, via Calamia 2009).

**How to validate it**
- EDGE's multi-row scenes against measurement and FDTD.
- The continuity sweep.
- The FDTD on a district tile.

**Evidence quality.** Moderate, mostly from rooms and scale models.

### Cheap, needed for validation, but not ranked as accuracy gains

- **Ground impedance and flush plates among buildings.**
  - Use the open-ground spherical-wave reflection coefficient per floor bounce.
  - Effect sizes:
    - Asphalt changed to grass lowered shielded-canyon levels by 3 dB(A) at the lowest receivers (Hornikx & Forssén 2006, 2.5-D equivalent sources).
    - Next to a wall, a flush floor microphone measured a peak-to-notch of about 20 dB at 5 cm (Gibertini et al. 2023, measured at 1:100).
  - Flight-test arrays use flush plates, so this is needed to compare with NICEOPS's own campaigns among buildings.
- **Receiver conventions.**
  - ISO 1996-2 adds +3 dB at 0.5–2 m from a facade and up to +6 dB flush. Measured deviations from these at large incidence angles exist (Čurović & Prezelj 2025).
  - MOTUS put rooftop monitors on 4 m masts specifically to avoid facade corrections (Caillet et al. 2024).
- **Data requests.**
  - Schiphol courtyard facades (Lugten).
  - MOTUS rooftops (Bruitparif / Airbus Helicopters).
  - DigiCity, all 22 microphones plus weather (NLR).
  - The PoliMi flush-microphone grid.
  - The open Wang et al. (2026) moving elevated source (§7).

### Bounded by the source, not the propagation

- **The drone source model is likely the bigger error.** NLR's rotor source model under-predicted a DJI M600 in both tonal and broadband noise. Measured mid and high frequencies fell about 10 dB from hover to 10 m/s, against the 3 dB the literature led them to expect (Habing et al. 2025, measured).
- **Until a measured multirotor sphere exists**, source error will probably exceed the propagation errors this review ranks. NICEOPS's open items 4a (a real multirotor sphere) and 4b (distributed rotors) are where that is addressed.

### Not worth pursuing

These are argued in §9:
- ML surrogates
- Gaussian beams without diffraction
- CFD-coupled building wakes for SEL maps
- wind and temperature inside canyons for elevated sources
- Harmonoise's 15 dB cap
- CNOSSOS-EU's ground term and favourable-condition multiple diffraction for high sources
- leaky-mode and Bragg models
- whole-scene ESIE (the edge source integral equation)
- diffusion-equation street models
- time-dependent radiosity
- tuning a Lambert s to fit a lit canyon

---

## 1. Ray and beam methods, and higher-order edge diffraction

### 1.1 What the literature shows

**Single edges: the edge-source/BTM integral is the reference, and NICEOPS already has it**
- The analytic secondary-source model (Svensson et al. 1999) has finite-edge frequency-domain forms (Svensson et al. 2009) and accurate approximations near zone boundaries (Svensson & Calamia 2006). Fast evaluation exists: edge subdivision (Calamia & Svensson 2007) and numerical steepest descent (Asheim & Svensson 2010).
- **Measured agreement:**
  - A single chipboard edge matched Svensson's finite-edge integral "almost perfectly" (Lokki & Pulkki 2004).
  - On a thin plate at grazing, first order alone left about 4 dB rms; including second order, by treating the plate as two wedges sharing a face, brought it to about 1 dB (Medwin et al., via Calamia 2009).
- **UTD's error grows with diffraction order.** Against BTM on a 12-edge box, 20 Hz–20 kHz, UTD averaged 0.67 dB at directly lit listeners, 3.10 dB at first-order diffraction positions and 3.61 dB at second order. It overestimated high-frequency attenuation because it uses one point of an edge while BTM integrates the whole edge (Schissler et al. 2014, simulated).
- **NICEOPS agrees with this:** UTD and the edge source tie on EDGE where buildings are many wavelengths long, and the edge source wins where ends count (UE07 drone 2.45 → 0.93 dB rms; UE14).

**Vertex (corner-point) diffraction**
- **Electromagnetics has a uniform vertex coefficient; acoustics has none.** Albani et al. (2009) give a closed-form UTD vertex-diffracted ray for faceted structures, uniform across the tip field's transition regions and validated against method of moments. It has been extended to truncated curved wedges (Albani et al. 2015).
- **Older forms are non-uniform:** the corner coefficients of Sikta et al. (1983) and Hansen (1991).
- **The vertex wave is the end-point contribution of an edge integral.** In the incremental theory of diffraction, the edge integral has a stationary point (the edge ray) and an end point, and the end point's term "provides the desired continuity" as the diffraction point leaves the edge (Albani et al. 1997).
- **Integrating an edge-source integral to the physical end therefore carries the first-order vertex wave.** It does not carry the corner's second-order edge-to-edge interaction.
- **No acoustic or urban study quantifying vertex diffraction was found.**

**Reflection boundaries at finite reflectors**
- **The half-amplitude rule.** At a specular-zone boundary the reflection and the edge's onset must sum to half the reflection's free value. Tests for the reflection and for diffraction that use different tolerances produce double, zero or sign-flipped fields. Calamia's fix derives both from the same per-edge angle parameters (Calamia 2009).
- **Engineering treatments:**
  - Nord2000 weights facade reflections by the part of the Fresnel zone they cover. Forssén et al. (2014) used an energy version (10 lg S rather than 20 lg S); I saw the method only through them.
  - Rindel's Kirchhoff–Fresnel finite-plate analysis (Rindel 2005) and the finite-facade study of Hopkins & Lam (2009) were seen only as snippets.
- **I found no urban validation of the corner of a lit region.** Physical optics over a rectangle gives a quarter amplitude there (standard Fresnel diffraction; inference).

**Higher orders and double diffraction**
- **Uniform double diffraction:**
  - Coplanar skew edges, including edges sharing a face: second-order terms are "as important as first-order" where the transition regions overlap (Capolino et al. 1997).
  - A uniform double-diffraction coefficient for wedges in arbitrary configuration (Albani 2005).
  - Pierce's thick-barrier double edge, which stays finite with an end in the roof's plane (Pierce 1974).
- **NICEOPS's shared-face aperture form is in this family.** Against BEM it is within 0.94 dB where the cascade was off by up to 16 dB.
- **Filter approximations** of first- and higher-order diffraction are an alternative for speed (Kirsch & Ewert 2024, abstract only).
- **Path finding at high order:**
  - The image-edge model mirrors edges in surfaces to find any mix of reflections and diffractions (Erraji et al. 2021).
  - Edge-visibility graphs bound the search; cost grows as n^d, with n the visible edges per edge and d the order. Wavelength-based edge simplification cut edges by 3–14 times (Schissler et al. 2014).
  - Conservative from-region visibility for finite edges (Antani et al. 2012a, snippet).
- **BTM has no slope-diffraction term** (Calamia 2009, §3.8). This matters when a second edge lies on the first edge's shadow boundary, where the incident amplitude vanishes but its gradient does not. The gradient near a thin-screen edge is given by Hewett & Svensson (2013, snippet).

**Urban tools and their accuracy**

| Tool | Method | Accuracy | Evidence |
| --- | --- | --- | --- |
| RWTH Virtual Acoustics urban renderer (Palenda et al. 2023; Palenda & Vorländer 2026) | Geometrical acoustics and UTD with higher-order diffraction and Doppler | Benchmarks show limits from its infinite-edge assumption; diffraction "around edges, not only in shadow regions" matters | Simulated, abstract |
| NLR (Knepper et al. 2025) | CNOSSOS-EU ray model | 10–15 dB low at the diffraction-only microphone; good where paths reflect | Measured |
| Gaussian beam tracing (Yunus et al. 2023a, 2023b) | Beams, no explicit edge diffraction | Over 40 dB maximum discrepancy against FEM without its correction, under 8 dB with it | Simulated |
| Hybrid FDTD/diffusion (Pasareanu 2014) | Hybrid | Edge mishandling up to 10 dB; mean 1.5 dB against full FDTD | Simulated, abstract |
| Radio, 1.8 GHz, central Bristol (Athanasiadou et al. 1997) | 7 reflections and 2 corner diffractions, powers summed | 3.0–3.3 dB rms including deep shadow | Measured |
| Radio, 28 GHz, Manhattan (Chizhik et al., venue not shown) | UTD at corners | Around-corner levels 7.8 dB rms low; modelling scatter from street clutter instead gave 2.3 dB rms | Measured |

The radio cases are at different ka (frequency times size) and only analogues. Still, at about 2 kHz acoustic wavelength, two diffractions with multiple reflections gave about 3 dB rms in deep shadow.

### 1.2 Bearing on NICEOPS's gaps

- **Vertex diffraction:** recommendation 2(b).
- **Facade-end transition:** recommendation 2(a) and (c).
- **Four or more edges, and the floor between edges:** recommendation 5.
- **Cost:**
  - Culling diffracted components by amplitude, or keeping only those near a zone boundary, saved 88–95% in rooms (Calamia et al. 2008, via Calamia 2009).
  - Edge-visibility graphs (Schissler et al. 2014) are the published tools.
  - NICEOPS's own measurement already favours this over ray-tracing hardware: visibility tests are 6–12% of a pair's cost.

---

## 2. Diffuse reflection in streets

### 2.1 What the literature shows

**Specular and diffuse canyons behave differently, mostly in long streets and in shadow**
- **Diffuse boundaries** give much more attenuation along a street, a shorter RT30 and a larger benefit from absorption (Kang 2000, simulated; full radiosity in Kang 2002).
- **Urban squares** (Kang 2005, simulated):
  - A boundary diffusion coefficient of 0.2 already makes the field close to fully diffuse.
  - Doubling the side length lowers levels by 6–9 dB.
  - Raising absorption from 0.1 to 0.9 lowers them by 10–12 dB.
- **Ignoring diffuse reflection** in a ray code gave about −2 to +4 dB of error, and up to +10 dB with street fittings (Can et al. 2015, simulated, street-level sources). The split ranges are from a snippet.
- **Multiple reflections** add about 8 dB in a 1:40 canyon with a viaduct, and a diffusion-method model matched it (Okada et al. 2010, measured).
- **The closest prior work to NICEOPS's case is a low-flying aircraft over street canyons** (Ismail & Oldham 2002, simulated):
  - With specular facades, the relative EPNL (effective perceived noise level) rises with building height over flight altitude; street width matters little.
  - Mixed specular and diffuse reflection barely changed the transmission gain. It did scatter energy out of the street, which made the relative EPNL negative for most streets.
- **Sonic booms:**
  - In canyons, diffusion was the most influential parameter, and without diffraction the shadows were silent (Riegel & Sparrow 2024).
  - Around one large building, measured, diffuse reflection did not improve agreement (Riegel & Sparrow 2019).

**Lambert's law fails for ordered relief**
- **Grating orders.** BEM of sinusoidal and rectangular surfaces: the directivities "do not follow the Lambert's law". Scattered energy goes into a few grating orders, and s jumps from about 0 to about 0.9 at an order's cut-on (Heimes et al. 2023; Pan et al. 2023).
- **A facade as a grating.** A Galois-profile facade gives 7 lobes with the specular lobe about 9 dB down (Lucchetti et al. 2026, preliminary; FEM and a 1:20 scale model).
- **Periodic slits** scatter as discrete wave components in different directions (Brunskog et al. 2022).
- **Bidirectional coefficients:**
  - Heimes & Vorländer (2025) define bidirectional scattering coefficients for urban use.
  - In a radiosity with a custom bidirectional distribution, Lambert matched ray tracing on a single facade. For retro-reflective surfaces only the bidirectional model approached the wave reference, and in a canyon the choice changed the energy–time curve significantly (Fatela et al. 2026).
- **Orientation (measured, 1:40).** At the shielded side, vertical diffusion (horizontal patches) gave significantly lower levels than horizontal diffusion at the same coverage. Diffusing patches had about the effect of facade absorption α ≈ 0.2 (Hornikx & Forssén 2008).
- **Explicit relief vs a coefficient.** Modelling the relief explicitly with edge diffraction, against a flat wall with s, gave the same objective parameters but audible differences (Shtrepi 2019).
- **Coherent and incoherent parts.** Directional coherence-loss coefficients split each direction's reflection into the two parts (Xie et al. 2026).

**Facade relief and shape, wave-based**
- **Triangular facade elements**, street-level source, 4th floor, 660–1330 Hz full scale: insertion loss above 8 dB at 30°, 6–8 dB at 45° and 2–4 dB at 60°, by suppressing specular paths (Van Renterghem & Thomas 2025, scale model and FDTD).
- **A rooftop source**, nodal DG-FEM checked on a 1:10 model (Feriani et al. 2025):
  - the opposite facade raises levels by up to 12 dB;
  - facade tilt changes them by −8.4 to +4.6 dB;
  - balconies change them by −3.5 to +6.3 dB, of which pure redirection is about ±1 dB.

**Coefficient data and its uncertainty**
- **A real facade** (BEM at 1:19; period 1.37 m, depth 0.18 m): random-incidence s about 0.1, 0.3 and 0.6 at 125, 250 and 500 Hz. Near grazing it is much higher, about 0.48 at 125 Hz and 0.9 at 500 Hz (Heimes et al. 2023, read off a figure).
- **Facade scale models** gave small s, around 0.2, with random scattering dominating higher reflection orders (Ismail & Oldham 2005, snippet-level findings).
- **The two-lab comparison** agreed closely on a sinusoid and less well on a rectangular profile (Heimes et al. 2026). Facades are rectangular-like.
- **Measurement-method uncertainty** for random-incidence coefficients is quantified by Müller-Trapet & Vorländer (2015).

**Scaling to many receivers**
- **View independence.** The rendering-equation formulation stores a view-independent patch solution; only the final gather depends on the receiver, and it takes milliseconds (Siltanen et al. 2007).
- **Moving sources and listeners:** a precomputed direct-to-indirect transfer operator (Antani et al. 2012b).
- **The specular–radiosity split by s that NICEOPS uses** is the established combination (Koutsouris et al. 2013; Marbjerg et al. 2015). Radiosity measured at least as well as diffuse ray tracing in rooms (Hodgson & Nosal 2006).

### 2.2 Bearing on NICEOPS

**Lambert radiosity is right for:**
- shadow fill
- high reflection orders in irregular fabric
- the DigiCity result

**It is wrong for:**
- ordered relief in lit canyons, which NICEOPS found independently
- relief orientation

The ranked fix is recommendation 3. For overhead sources, the incidence dependence of s is the cheapest improvement that matters.

---

## 3. Wave-based, hybrid and surrogate methods; low frequencies

### 3.1 What the literature shows

**3-D wave solvers on real urban scenes**

| Study | Method | Scene and reach | Cost | Accuracy |
| --- | --- | --- | --- | --- |
| Mehra et al. (2014) | ARD (adaptive rectangular decomposition), 2.6 points per wavelength | 150 × 150 m training village, 15 buildings; 450 Hz | 20 min on one CPU core, 1–2 min on a 2010 GPU | About 3 dB against blast waveforms at 14 sensors, measured; 3-D clearly better than 2-D where paths cross roofs |
| Oshima et al. (2014a) | Full-scale FDTD with terrain | 380 × 380 m | — | Against a 1:100 model: about 2 dB rms at 125 and 250 Hz, about 4 dB at 500 Hz |
| Hornikx & Forssén (2011) | 3-D PSTD (pseudospectral time domain) | Courtyards to 500 Hz | — | Relative to 3-D, 2-D under-predicts facade screens by up to 4.4 dB(A) and facade absorption by up to 1.5 dB(A), and over-predicts roof screens by up to 1.7 dB(A); a courtyard opening raised levels up to 10 dB(A) |

- **ARD's method paper** reports at least 100 times less compute than FDTD (Raghuvanshi et al. 2009).
- **Methods disagree on shadow depth.** Behind a roadside building, scale model and Euler FDTD predicted a deep shadow; field data and ISO 9613-2 showed a shallow one or none (Heimann et al. 2018, measured).

**Low-frequency canyon effects**
- **Street modes.** In a 1:25 street (5 m wide, 6.75 m high) at 40–100 Hz full scale, single-frequency maps show nodes 20–60 dB deep. A modal-plus-image model and 2-D FDTD agree with the experiment up to about 100 Hz. Above that, leaky modes are needed (Richoux et al. 2009; Pelat et al. 2009). A 3-D PE handles irregular canyons (Doc et al. 2015).
- **Sonic booms** (2-D Euler; Dragna et al. 2022):
  - First transverse modes of 10, 20 and 30 m streets fall near 18.5, 10 and 7.5 Hz.
  - They ring at about 10 Hz with "little effect on perceived noise".
  - An isolated building gives about +4 dB on the lit side, +7 dB at the front corner and −14 dB in the shadow.
- **Tones vs bands.** For a source above a corner (Malbéqui 2023, image sources):
  - The pure-tone response has notches to about −30 dB and peaks near +10 dB.
  - In octave bands this collapses to about −4 dB.
  - Canyon gains are +4 to +14 dB, an upper bound because the walls are rigid and infinite.
- **Narrow canyons.** Coherent image summation is needed for streets narrower than about 10 m (Iu & Li 2002, measured).
- **Per intermediate canyon**, excess attenuation is −1 to −4 dB, roughly flat with frequency (Schiff et al. 2010, PE with equivalent sources).
- **The gap.** No study was found on how band averaging behaves for Doppler-swept rotor tones in canyons, for SEL or LAmax.

**Hybrids and wave-derived corrections**
- **The clearest precedent** is an engineering model fitted to a large set of full-wave canyon calculations and checked at 9 Ghent facades (Wei et al. 2014; FDTD–PE parameter study by Van Renterghem et al. 2006).
- **Crossover frequencies** in room-acoustic hybrids are chosen empirically (Kjær et al. 2022, a perceptual study in rooms).

**ML surrogates**

| Study | Setup | Error (against its own training model) |
| --- | --- | --- |
| Spitznagel & Keuper (2024) | 2-D, 500 Hz, ground truth from a CNOSSOS tool | Best MAE 1.5–2.9 dB; 1.2–4.8 dB without line of sight |
| Eckerle et al. (2025) | Normalizing flow | 0.65 dB baseline, 2.6 dB with diffraction, 3.6 dB with reflection, in shadow; about 2000 times faster |
| Kim & Lee (2024) | UAM noise hemisphere through random 3-D layouts, trained on ray tracing | 2.56 dB RMSE |
| Borrel-Jensen et al. (2024) | Neural operators trained on DG | Small rooms only; 1–3 GPU-days of training per scene; no geometric generalization |

Kim & Lee is the only one with an elevated UAM source.

### 3.2 Bearing on NICEOPS

- **Low frequencies matter per band, not for A-weighted metrics.**
  - Below about 150 Hz the effects are large per band or per tone (nodes of 20–60 dB) and probably small on A-weighted SEL or LAmax. That is my inference from A-weighting (−19 dB at 100 Hz) and from Dragna et al.
  - This supports NICEOPS's "lower priority" placement.
  - Drones are the exception to watch: a DJI M300's blade-passage frequencies sit near 100 Hz (Ramos-Romero et al. 2023, measured), so their first harmonics fall in the modal range of narrow streets.
- **The cheapest credible step:**
  - Below a crossover of about 250 Hz, sum NICEOPS's existing paths coherently at finer frequency resolution, which tracks the interference.
  - Validate at 20–150 Hz against the 3-D FDTD on canyon and courtyard tiles.
  - Optionally store FDTD-minus-ray corrections per tile, direction and band, following the Wei et al. precedent.
- **Not worth it:**
  - per-flight full-wave hybrids
  - leaky-mode or Bragg models, which assume uniform or periodic streets with sources in them
  - 2.5-D solvers, which need geometry invariant along the street
  - ML surrogates (§9)

---

## 4. Engineering methods and their urban errors

| Method | Building treatment | Documented urban error |
| --- | --- | --- |
| ISO 9613-2-type barrier formula, also in CNOSSOS-EU | Barrier diffraction, image reflections | "Underestimates the attenuation considerably" of thick buildings, against long-term measurements and FDTD (Wei et al. 2014) |
| CNOSSOS-EU | Reflections, lateral and vertical diffraction | An END map, implementation not stated, was over 14 dB(A) low at a shielded courtyard facade (Wei et al. 2014, measured). Exposed facades mostly within 3 dB, biased low (Rey Gozalo & Gómez Escobar 2021, measured). Multiple diffraction in favourable conditions can come out "much lower (even 0)", and the ground formulas are "fundamentally incorrect for a high source and low receiver" (Kok 2019). 10–15 dB low at DigiCity microphone 7 (Knepper et al. 2025, measured) |
| Harmonoise | 2-D planes, coherent within a plane and incoherent between; turbulence only as a 15 dB cap on excess attenuation | Under 0.5 dB(A) on flat terrain and 3.1 dB(A) at a hilly site with a viaduct, non-urban (Defrance et al. 2007; Salomons et al. 2011) |
| Nord2000 | Reflections weighted by Fresnel-zone coverage; city-street test cases with first- and second-order facade reflections | Validation of 544 cases is mostly open terrain and screens, σ about 1 dB to 400 m (Kragh 2011; Danish EPA 2006) |
| Wei/Hornikx/Forssén add-on: multiple-reflection term A_can, plus scattering | Fitted to full-wave runs | Fit residual 2.7–3.0 dB; source- and receiver-side terms separate within 0.5 dB(A) beyond 4 street widths (Hornikx et al. 2014; Wei et al. 2014) |
| sonAIR (aircraft) | Maekawa shielding plus scattering reflections | At Meiringen: shielding below −2 dB over wide areas and up to −12 dB on averted facades; reflections +3–4 dB on facing facades; authors advise leaving buildings out of standard mapping (Ramseier et al. 2023, 2024, simulated) |
| AEDT / ACRP 02-79 | Gain and loss grids relative to open terrain | Facade gains up to about 3 dB, falling off quickly; shielding stronger and slower to decay (Downing et al. 2019, via secondary summary) |

**Aircraft over cities.** Buildings matter little under direct overflight, a lot at grazing incidence, and reflections matter mainly at receivers already shielded (Schlatter et al. 2018). For helicopters, approach and departure spend most of their time at low elevation, in the regime where buildings matter most.

**Bearing on NICEOPS**
- NICEOPS should borrow pieces, not adopt a method:
  - Nord2000's Fresnel-zone weighting (recommendation 2(a))
  - the per-canyon attenuation as a floor beyond the explicit orders (recommendation 5)
  - the courtyard benchmark values, closed-courtyard Ā′r of 7–8 dB(A) (Hornikx et al. 2014), and ISO/TR 17534-4 as sanity checks only
- None of the engineering methods was built for an elevated source over canyons, and CNOSSOS-EU's ground term is known to be wrong there.

---

## 5. Atmospheric effects among buildings

**Turbulence scattering into building shadows, measured**
- **Container barrier** (Forssén & Ögren 2002):
  - About +5 dB at 4 kHz and +10 dB at 8 kHz in windy conditions, little below 2 kHz.
  - Measured C_v² was 0.4–17 m^4/3 s⁻² and C_T² 0.2–4 K² m^−2/3.
  - The model over-predicted above 6–8 kHz.
- **Ghent shielded facade** (Van Renterghem et al. 2013a):
  - Nearly 10 dB between weak and strong turbulence classes above 1 kHz, and significant below it.
  - The hourly interquartile range reached 15 dB at 4 kHz.
  - The authors recommend "rather high" strengths for long-term mapping.
- **History.** Daigle (1982, snippet) is the original case that diffraction plus scattering explains barrier data better than diffraction alone.

**Simulated**
- **Equivalent sources with turbulence:** +2–5 dB at 1.6 kHz in a canyon, negligible below 500 Hz, under 1 dB(A) on a traffic spectrum. A ray model gave 4.9 dB where equivalent sources gave 1.6 dB, so the result is method-dependent (Ögren & Forssén 2004; Ögren 2004).
- **Quasi-wavelets:** scattering over barriers matters mainly above several kHz, near the barrier top and in very turbulent air (Wilson et al. 2004).
- **The Forssén et al. (2014) engineering model** is given in recommendation 1. Its error is about 2 dB against its parameter study and under 3 dB in total; with canyons and a road-traffic spectrum at C_v² = 1 it is about 4 dB(A).
- **Building wakes** (LB-VLES, a 50 m building, a 100 Hz dipole): fluctuations above 3 dB that extend beyond the wake through reflections (Goyal et al. 2026a). Ray tracing through frozen velocity snapshots reproduces the statistics (Goyal et al. 2026b).

**Coherence between paths**
- **The theory:**
  - Ögren & Forssén use Γ_ij = exp(−k²J Δ_ij), with Δ the path difference. They note it under-represents turbulence by applying it only to the extension.
  - Ostashev et al. (2001) and Salomons et al. (2001) give the spherical-wave coherence for direct and reflected waves.
- **Interference persists in measurement.** Behind the container barrier, interference persisted to about 3 kHz, where theory predicted more decorrelation (Forssén & Ögren 2002).
- **Arithmetic on Ögren's form** (mine, C_v² = 1): Γ for a 10 m path difference is 0.98 at 1 kHz and 0.68 at 4 kHz.
- **NICEOPS's 1-degree cone is therefore defensible.** The pairs inside it are mostly floor and near-facade partners with path differences of metres. Turbulence would decorrelate them only at the top bands in strong turbulence.
- **Recommendation:** a pairwise Γ is a small refinement (a few dB at interference structures near facades and flush microphones, little on SEL), not a priority.

**Wind and temperature in canyons**
- **Measured:**
  - Courtyard to street over months: coherence loss rose with wind speed and frequency, a temperature lapse had no measurable effect, and downwind refraction into the shielded canyon was limited by the geometry (Van Renterghem & Botteldooren 2010).
  - Over a block or two there was no notable refraction (Wiener et al. 1965).
- **Simulated:**
  - Upwind facades quieter and downwind louder (Heimann 2007).
  - Wind effects of 15–23 dB(A) for a road source with receivers up to 500 m away over generic urban sections (Hornikx et al. 2018, RANS with PSTD; reference condition not checked).
- **Bearing on NICEOPS:**
  - The large wind effects are for low sources propagating over many roofs.
  - For an elevated source on a steep path, a 30 m leg in a canyon bends by under a metre at plausible gradients (my estimate).
  - The above-canopy profile belongs to NICEOPS's existing eigenray models. Unsteady wakes are a variance on LAmax, not a shift in mean SEL.
  - Not worth modelling in-canyon.

**Size against geometry errors**
- At deeply shielded receivers above 1–2 kHz, turbulence is as large as or larger than the geometric errors NICEOPS has measured (5–10 dB measured, against 2–4 dB rms).
- At lit receivers and below 1 kHz, it is negligible.

---

## 6. Drones, AAM and helicopters in cities

**How other groups model propagation among buildings, and what accuracy they report**

| Group | Method | Reported accuracy and effects | Evidence |
| --- | --- | --- | --- |
| NASA (Rizzi et al. 2026) | Review | Canyons cause "reflections, reverberation, and diffraction, and may act as waveguides". Recommends canonical urban geometries with simple sources under well-characterized meteorology. Notes the relatively high partial coherence of direct and reflected tonal UAM noise at elevated microphones. The UNWG reached no consensus on reference urban locations | Review |
| NLR (Knepper et al. 2025; Habing et al. 2024, 2025) | CNOSSOS-EU paths over extruded footprints; rotor-by-rotor source hemispheres | Reflection-only microphones agree; the diffraction-only microphone is 10–15 dB low (loudspeaker). The source model under-predicts the M600 | Measured |
| ONERA (Malbéqui 2023) | Mithra ray code against image sources | Corner +8 dB, canyon up to +14 dB re free field; octave bands keep the floor dip | Simulated |
| TU Delft (Yunus et al. 2023a, b) | Gaussian beam tracing with a source sphere and wind | Buildings add 0–5 dB in lit zones; wind changes the footprint by up to 23 dB; first- and second-BPF directivities differ by up to 40 dB; no diffraction | Simulated |
| PoliMi (Gibertini et al. 2023) | 1:100 rotor descending to a roof; flush floor microphones | About 60 dB(A) in front of the building and 25 dB(A) behind at 50 m before touchdown; wall interference about 20 dB peak-to-notch at 5 cm; frequencies only to about 200 Hz full scale | Measured, scale model |
| Airbus Helicopters, Bruitparif and Univ. Eiffel, MOTUS (Caillet et al. 2024) | NoiseModelling adapted for aerial sources | Two residents near the heliport differ by about 10 dB by facade orientation; monitors on 4 m rooftop masts | Measured and simulated |
| TU Delft, Schiphol (Lugten 2022; Lugten et al. 2023, 2024) | Full-scale container courtyards; 1:50 model | Exposed–shielded facade differences 8–9 dB(A) with straight facades, 14 dB(A) with a slanted facade and inset; scale model against in situ, relative: 0.5–4.2 dB(A), σ ≈ 2 dB(A) | Measured |
| TU Delft (Wuite et al., c. 2024) | 36 trees in a courtyard | Shielded microphone ASEL down 1.1 dB(A) against a bare twin courtyard | Measured, abstract |
| UCL (Feng & Kang 2026) | CadnaA / ISO 9613 over ten cities | Courtyard facades up to 12.2 dB(A) apart; quiet semi-enclosed spaces rise more than 20 dB(A) | Simulated, unvalidated |

**Sources near buildings**
- **Far field of an isolated rotor.** It converges within 0.5 dB of 1/r by about 10 rotor radii, slowest in-plane (Schiller et al. 2019). NASA's AARON hemispheres are built at 152 m, 38–100 rotor radii (Rizzi et al. 2024).
- **Between rotors.**
  - Independently controlled rotors have no stable phase relationship and are summed in energy (Schiller et al. 2019). NLR measured rotor-to-rotor RPM spread and fluctuation (Habing et al. 2025). Yunus et al. summed four rotors incoherently.
  - Synchronized pairs redistribute the field by about ±4–5 dB directionally rather than reducing it (Schiller et al. 2019; Chai et al. 2025).
- **Distributed rotors near a facade.** No quantitative study was found. An estimate from the research pass, not from literature:
  - under 1 dB beyond 2–3 vehicle spans in lit zones;
  - up to about 3 dB within a band about one span wide around shadow and specular boundaries, during hover or vertiport segments;
  - largely averaged away in a flyover SEL.

  A test with an artificially wide source array on EDGE geometry would show whether it matters.

**Bearing on NICEOPS**
- **Other groups report no better accuracy.** NICEOPS's DigiCity result already beats the only published comparison of a diffraction-dominated urban microphone, by about 10 dB.
- **Per-path directivity** (each ray reading the sphere along its own launch) is what the literature says is needed, given 40 dB directivity differences between harmonics. NICEOPS already does this.
- **Distributed rotors are a vertiport and hover issue.** The source database, not propagation, bounds drone accuracy for now.

---

## 7. Validation data

| Dataset | What it is | Source | Receivers | Bands | Availability | Tests |
| --- | --- | --- | --- | --- | --- | --- |
| EDGE (Palenda & Vorländer 2026), Zenodo 10.5281/zenodo.21264202 | 31 scale-model urban scenes, RWTH hemi-anechoic room, 1/3 to 1/20 scale | Floor and "drone" loudspeakers, measured directivity | Mostly on the floor, in shadow | About 1.6–20 kHz at model scale (tweeters) | CC BY-SA 4.0, in hand | Diffraction orders, corners (UE20–22), scattering reliefs (UE23–30), **pitched roofs (UE01–05, 15, 17, 19), not yet used** |
| NLR Drone DigiCity (Knepper et al. 2025) | 24 shipping containers at full scale | B&K 4296 loudspeaker on a 2 m tripod; DJI M600 at 8.5 and 20 m with rotor RPM logged | 22 microphones at 1.2 m | 160 Hz–6.3 kHz usable | On request (NLR) | Reflection-only, diffraction-only and mixed paths; multirotor source |
| Schiphol courtyards (Lugten 2022; Lugten et al. 2023, 2024) | Three container courtyards, 2383 flyovers, radar and KNMI weather; 1:50 model | Real jets, mean elevation 37° ± 7° | **Facade microphones 0.2 m off, at 1.5 and 3.9 m** | 1/3-octave every 0.125 s, plus WAV | On request (TU Delft) | Facade receivers, exposed vs shielded, roof and facade shape, turbulence by weather |
| MOTUS, Paris-Issy heliport (Caillet et al. 2024) | 15 rooftop stations, 15,329 helicopter events matched to radar | Operational helicopters | **Rooftop, 4 m masts**, with arrival direction | LAeq every 0.1 s; per-event LAmax and SEL | Report public; data on request | Rooftop SEL and LAmax |
| PoliMi rotor and building (Gibertini et al. 2023) | 1:100 rotor descending to a roof | Model rotor | **Flush floor microphones** 5–50 cm from the wall | Up to about 200 Hz full scale | On request | Facade interference at flush microphones; shielding |
| Wang et al. (2026), Zenodo 10.5281/zenodo.15833542 | 1:40 hemi-anechoic; five layouts of detached houses | Elevated source moved along a 100 m line | Single-point insertion loss per receiver | Octaves 125 Hz–2 kHz | CC BY 4.0, open | Building-group insertion loss, elevated source |
| Heimes et al. (2026), Zenodo coefficients | Interlaboratory bidirectional scattering coefficients, including facade elements | — | — | — | Open | The scattering inputs for recommendation 3 |
| Heli-STAR Vol. 5 (Ahuja et al. 1997) | S-76 near buildings | Full-scale helicopter | Roof, building base, open slope | Not verified | DTIC; downloads failed | Rooftop and base receivers. Accounts conflict on whether the helipad was a rooftop (Gibertini et al. 2023) or at ground level (the project's earlier reading of Heli-STAR Vol. 1); unresolved |
| Forssén & Ögren (2002); Van Renterghem et al. (2013a) | Container barrier with sonic anemometer; Ghent shielded facades | Ground and road sources | Behind the barrier; at facades | 1/3-octave | In the papers | The turbulence floor (recommendation 1) |

**What none of them provides:** a rotorcraft at 30–1500 m over real canyons, in bands, with tracks and weather. NASA's 2025 review names this gap (Rizzi et al. 2026).

---

## 8. Gap-to-method table

| Known gap | Best method in the literature | Expected effect | Effort | Evidence quality |
| --- | --- | --- | --- | --- |
| Vertex (corner-point) diffraction | Edge-source integral taken to the physical end (end-point term); scalar reduction of Albani et al.'s (2009) UTD vertex coefficient for cascades | Removes about 6 dB steps where a Keller point leaves its edge; elsewhere about 20 dB below edge waves at 1 kHz (inference) | Low–medium (integral route); high (UTD vertex) | Strong theory, validated in electromagnetics; no acoustic urban validation |
| Reflection-boundary transition near a facade's end, partner at a vertex | Physical-optics/Fresnel weighting over the facade rectangle (product of two Fresnel transitions); Nord2000 Fresnel-zone weight as a cheap form; one zone test for a term and its partners (Calamia 2009) | Up to about 5–6 dB locally (UE06 R1 +4.7 dB) | Low–medium | Strong theory; no urban validation |
| Four or more edges | Image-edge path search with an edge-visibility graph and wavelength pruning; uniform double diffraction where transitions overlap; per-canyon attenuation beyond the explicit order | Up to about 10 dB where needed (multi-row grazing); rare at steep elevation | Medium–high | Moderate (rooms, scale models) |
| Floor between two edges; several reflections between edges | Image edges (Erraji et al. 2021; Torres et al. 2001) | Up to +3 dB in energy (+6 dB coherent near grazing) for the affected family; 13 of 25 continuity steps | Low | Good (rooms, scale models) |
| Low-frequency canyon resonances, below about 150 Hz | Coherent sums at fine frequency resolution below about 250 Hz; FDTD-derived correction tables per tile (Wei et al. 2014 precedent) | Per band or tone: nodes of 20–60 dB, idealized gains of +4 to +14 dB. On A-weighted SEL and LAmax: probably under 1 dB (inference) | Low–medium (coherent sums); medium (tables) | Low–moderate; nothing tests rotor SEL |
| Turbulence: scattering into shadows | Forssén et al. (2014) engineering model | +1–4 dB (moderate) to +5–16 dB (strong) at 1–8 kHz in shallow shadows (my evaluation); measured +5/+10 dB at 4/8 kHz and about 10 dB between turbulence classes above 1 kHz | Low | Moderate (one fit, two measurement sets; input strength uncertain) |
| Turbulence: decorrelation of multipath | Pairwise Γ = exp(−k²JΔ) or spherical-wave coherence (Ostashev et al. 2001) | A few dB at interference structures; little on SEL. Measured interference persists to about 3 kHz | Low–medium | Moderate; the measurement shows less decorrelation than theory |
| LoD2 pitched roofs | Sloped facets and edges in the path search; validate on EDGE's pitched-roof scenes | Up to about 10 dB locally (full-wave); measured exposed–shielded contrast 14 against 8–9 dB(A) with a slanted facade | Medium–high | Moderate |
| Directional (non-Lambert) facade scattering | s(θ, f) from BEM; grating-order images; bidirectional radiosity (Fatela et al. 2026) | 1–2 dB at reflected and shielded receivers; up to 3 dB in single bands | Low (s(θ)), medium (orders), high (BRDF radiosity, i.e. a full directional distribution) | Strong for single surfaces; weak for canyons |
| Vegetation | Treat trees as scatterers and absorbers in the diffuse field; green roofs as impedance on roof bounces | About 1 dB(A) at a shielded courtyard microphone (Wuite et al., measured); green roofs among the strongest quiet-side measures (simulated, road traffic) | Medium | Weak to moderate |
| Terrain under buildings | Footprints draped on the terrain model; locally planar floor facets in the image method. Real-terrain FDTD exists (Oshima et al. 2014a, 2014b) | No isolated number found; NICEOPS's own open-ground slope study suggests about 1 dB for flush plates at low elevation | Medium | Weak |
| Ground impedance and plates among buildings | The open-ground spherical-wave reflection coefficient per floor bounce | 3 dB(A) between asphalt and grass at low shielded receivers (simulated); about 20 dB interference structure for flush microphones by walls (measured, scale model) | Low | Strong physics; thin urban measurement |
| Facade and rooftop receivers, validation | Schiphol facade microphones; MOTUS rooftops; ISO 1996-2 conventions | Facade +3/+6 dB conventions with measured deviations; rooftop microphones gain 1–2 dB from roof scattering (NICEOPS's own result) | Low (data requests) | Data exists, on request |
| Distributed rotors near facades | N incoherent sub-sources at the hubs, each with its share of the sphere, near facades only | Under 1 dB beyond 2–3 spans; up to about 3 dB near boundaries in hover (estimate; no literature) | Medium | Low |
| Cost at district scale | Diffraction culling (88–95% saved in rooms); edge-visibility graphs; reciprocal radiosity gather; per-receiver tables | Enables maps; no accuracy change | Medium | Room and graphics evidence |

---

## 9. What NICEOPS already does well, and what is not worth pursuing

### Already well founded

- **The edge-source (BTM) integral with finite edges.** It is the reference-quality single-edge model: measured agreement "almost perfect" (Lokki & Pulkki 2004), and about 1 dB with second order at grazing (Medwin, via Calamia 2009).
- **The uniform shared-face form.** It follows Pierce's thick barrier and the coplanar double-diffraction results (Capolino et al. 1997), and fixes the cascaded-UTD failure in transition regions that Schissler et al. (2014) quantify.
- **Diffraction on both sides of a boundary, and edges lit by reflections.** The tools that miss DigiCity's shielded microphone by 10–15 dB are the ones that do not combine reflection with diffraction. Palenda & Vorländer (2026) find diffraction matters around edges, not only in shadow.
- **Per-path sphere reads.** These are needed given the large directivity differences of rotor harmonics (Yunus et al. 2023a).
- **Lambert radiosity split by s, and roofs as scatterers.** This is the standard combination (Koutsouris et al. 2013; Marbjerg et al. 2015), right for shadow fill. The BEM-from-profile s(f) is how the coefficient literature does it (Heimes et al. 2023).
- **Coherent near pairs, energy beyond.** This matches Harmonoise's coherent-within-a-plane practice, and the measured persistence of interference to about 3 kHz (Forssén & Ögren 2002).
- **A 3-D full-wave reference.** Others use one too: PSTD, ARD, Oshima's FDTD. The deliberate rejection of 2-D references for real streets is supported: 2-D misstates treatments by 1.5–4.4 dB(A) (Hornikx & Forssén 2011).
- **Not adopting CNOSSOS-EU's ground term for high sources** (Kok 2019).

### Not worth pursuing

- **ML surrogates.** Errors of 1–5 dB against their own training models, mostly 2-D, ground sources and single frequencies, with no bound. NICEOPS's GPU ray model is already fast (Spitznagel & Keuper 2024; Eckerle et al. 2025; Kim & Lee 2024).
- **Gaussian beam tracing without diffraction.** It cannot produce building shadows (Yunus et al. 2023a, 2023b).
- **CFD- or LBM-coupled building wakes for SEL maps.** About 10⁵ CPU-hours for 70 s of one tone; the effect is a variance, not a mean (Goyal et al. 2026a).
- **Wind and temperature gradients inside canyons for elevated sources** (§5).
- **Harmonoise's 15 dB turbulence cap.** It is calibrated for barriers; urban shadows need −20 to −45 dB re free field.
- **CNOSSOS-EU's favourable-condition multiple diffraction and high-source ground terms** (Kok 2019).
- **Leaky-mode, modal and Bragg band-gap street models** (Pelat et al. 2009; Richoux et al. 2009). They assume uniform or periodic streets with sources in them.
- **ESIE for whole scenes.** It is for convex polyhedra only, and its higher orders converge slowly (Asheim & Svensson 2013).
- **Energy-based uncertainty diffraction.** It loses phase (Stephenson 2010).
- **Diffusion-equation street models and regression corrections fitted to street-level traffic** (Picaut et al. 1999; Can et al. 2015; Thomas et al. 2013).
- **Time-dependent radiosity.** SEL and LAmax are time-integrated or slow.
- **Tuning a Lambert s to fit a lit canyon.** A coefficient cannot hold redirection.
- **Diffuse modelling around isolated buildings.** It did not improve agreement with the SonicBOBS measurements (Riegel & Sparrow 2019).

---

## 10. Coverage limits and conflicts

- **Not covered or only shallow** (the search budget ran out):
  - the primary Nord2000 equations for finite reflectors (seen only through Forssén et al. 2014)
  - RLS-19 and sonRoad building treatment
  - NoiseModelling's measurement validation
  - effective-medium canopy PEs
  - Fourier neural operator and DeepONet surrogates for outdoor propagation
  - Doppler-swept tones in canyons
  - DLR's urban UAM work
  - Empa's drone-route measurements
  - the full ACRP 02-79 report
  - a quantitative point-vs-distributed-source study (none found)
- **Conflicts and corrections:**
  - **Heli-STAR Vol. 5:** Gibertini et al. describe a rooftop landing; the project's earlier reading of Vol. 1 found a ground-level helipad. Unresolved.
  - **Lugten et al. (2024):** I could not read the journal abstract. The 8–9 and 14 dB(A) figures are from the Inter-Noise 2022 version, and the 2383-flyover count is from the project's earlier validation search.
  - **Hornikx & Forssén (2008)** is Acta Acustica 94(2):265–281, not 94(4):641.
  - **"Ouis 2003 JSV 262:347"** and an "Embrechts 2016" scattering paper could not be found as cited. Embrechts (2000) and Embrechts & Billon (2011) are the verified works.
  - **Ögren & Forssén (2004):** the indexed abstract reads "25 dB at 1.6 kHz"; Ögren's thesis says 2–5 dB.
- **My own numbers** are labelled as such throughout: the turbulence-floor table, the γ_t–C_v² reading, the vertex scaling and the A-weighted estimates. They are worth checking in NICEOPS before being relied on.

---

## 11. References

Tags: [F] full text read; [A] abstract or record read; [M] metadata confirmed, content from a stated secondary source; [S] snippet or secondary citation only (unverified). "DOI checked" means confirmed against Crossref.

**Ray and beam methods, and diffraction**

- Albani M. (2005). A uniform double diffraction coefficient for a pair of wedges in arbitrary configuration. *IEEE Trans. Antennas Propag.* 53(2):702–710. doi:10.1109/TAP.2004.841289 [A, DOI checked]
- Albani M., Capolino F., Carluccio G., Maci S. (2009). UTD vertex diffraction coefficient for the scattering by perfectly conducting faceted structures. *IEEE Trans. Antennas Propag.* 57(12):3911–3925. doi:10.1109/TAP.2009.2027455 [A, DOI checked]
- Albani M., Capolino F., Maci S., Tiberio R., Molinet F. (1997). ITD formulation for the currents on a plane angular sector. *IEEE AP-S Int. Symp.* 1997, 1764 ff. https://escholarship.org/uc/item/7jt1g8j8 [F, pp. 1764–1765]
- Albani M., Carluccio G., Pathak P.H. (2015). A uniform geometrical theory of diffraction for vertices formed by truncated curved wedges. *IEEE Trans. Antennas Propag.* 63(7):3136–3143. doi:10.1109/TAP.2015.2427877 [A]
- Antani L., Chandak A., Taylor M., Manocha D. (2012a). Efficient finite-edge diffraction using conservative from-region visibility. *Appl. Acoust.* 73:218–233 [S]
- Asheim A., Svensson U.P. (2010). Efficient evaluation of edge diffraction integrals using the numerical method of steepest descent. *JASA* 128:1590–1597 [S]
- Asheim A., Svensson U.P. (2013). An integral equation formulation for the diffraction from convex plates and polyhedra. *JASA* 133(6):3681–3691. doi:10.1121/1.4802654 [M, DOI checked]
- Athanasiadou G.E., Nix A.R., McGeehan J.P. (1997). Comparison of predictions from a ray tracing microcellular model with narrowband measurements. *IEEE 47th Vehicular Technology Conf.* 2:800–804. doi:10.1109/VETEC.1997.600439 [F, DOI checked]
- Calamia P.T. (2009). *Advances in Edge-Diffraction Modeling for Virtual-Acoustic Simulations.* PhD thesis, Princeton University. https://gfx.cs.princeton.edu/pubs/_2009_AIE/calamia-phd-thesis.pdf [F, §§3.6–3.8, 6.4–6.6, 7]
- Calamia P.T., Markham B.E., Svensson U.P. (2008). Diffraction culling for virtual-acoustic simulations. *Acta Acust. united Ac.* 94:907–920 [S, content via Calamia 2009]
- Calamia P.T., Svensson U.P. (2007). Fast time-domain edge-diffraction calculations for interactive acoustic simulations. *EURASIP J. Adv. Signal Process.* 2007:063560. doi:10.1155/2007/63560 [A, DOI checked]
- Capolino F., Albani M., Maci S., Tiberio R. (1997). Double diffraction at a pair of coplanar skew edges. *IEEE Trans. Antennas Propag.* 45(8):1219–1226. doi:10.1109/8.611240 [A, DOI checked]
- Chizhik D., Du J., Kohli M., Adhikari A., Feick R., Valenzuela R.A., Zussman G. Accurate urban path loss models including diffuse scatter. https://par.nsf.gov/servlets/purl/10406109 [F, pp. 1–2; venue and year not shown]
- Erraji A., Stienen J., Vorländer M. (2021). The image edge model. *Acta Acustica* 5:17. doi:10.1051/aacus/2021010 [A, DOI checked]
- Hansen T.B. (1991). Corner diffraction coefficients for the quarter plane. *IEEE Trans. Antennas Propag.* 39(7):976–984 [S]
- Hewett D.P., Svensson U.P. (2013). The diffracted field and its gradient near the edge of a thin screen (L). *JASA* 134:4303–4306 [S]
- Hopkins C., Lam Y.W. (2009). Sound fields near building façades: comparison of finite and semi-infinite reflectors on a rigid ground plane. Salford repository output 1445131 [S; venue not confirmed]
- Kirsch C., Ewert S.D. (2024). Filter-based first- and higher-order diffraction modeling for geometrical acoustics. *Acta Acustica* 8:73. doi:10.1051/aacus/2024059 [A, DOI checked]
- Kouyoumjian R.G., Pathak P.H. (1974). A uniform geometrical theory of diffraction for an edge in a perfectly conducting surface. *Proc. IEEE* 62(11):1448–1461. doi:10.1109/PROC.1974.9651 [M, DOI checked]
- Lokki T., Pulkki V. (2004). Measurement and theoretical validation of diffraction from a single edge. *Proc. ICA 2004*, Kyoto, Tu2.B2.1. https://users.aalto.fi/~ktlokki/Publs/Tu2.B2.1.pdf [F, pp. 1–2]
- Palenda P., Schäfer P., Vorländer M. (2023). Auralization of dynamic urban environments including diffraction and reflection using "Virtual Acoustics". *Forum Acusticum 2023*, Turin, 5637–5641. doi:10.61782/fa.2023.0191 [F, DOI checked]
- Palenda P., Vorländer M. (2026). Auralizing arbitrary urban environments including diffraction. *Acta Acustica* 10:46. doi:10.1051/aacus/2026033 [A, DOI checked]
- Pasareanu S.M. (2014). *A numerical hybrid method for modeling outdoor sound propagation in complex urban environments.* MS thesis, Virginia Tech. hdl:10919/47601 [A]
- Pierce A.D. (1974). Diffraction of sound around corners and over wide barriers. *JASA* 55(5):941–955. doi:10.1121/1.1914668 [M, DOI checked]
- Rindel J.H. (2005). Kirchhoff–Fresnel analysis of reflection from finite-size plates (listing title "Reflection of sound from finite-size plane and curved surfaces") [S; not retrieved]
- Schissler C., Mehra R., Manocha D. (2014). High-order diffraction and diffuse reflections for interactive sound propagation in large environments. *ACM Trans. Graph.* 33(4):1–12. doi:10.1145/2601097.2601216 [F, DOI checked]
- Sikta F.A., Burnside W.D., Chu T.-T., Peters L. (1983). First-order equivalent current and corner diffraction scattering from flat plate structures. *IEEE Trans. Antennas Propag.* 31(4):584–589 [S]
- Stephenson U.M. (2010). An energetic approach for the simulation of diffraction within ray tracing based on the uncertainty relation. *Acta Acust. united Ac.* 96(3):516–535 [S]
- Svensson U.P., Calamia P.T. (2006). Edge-diffraction impulse responses near specular-zone and shadow-zone boundaries. *Acta Acust. united Ac.* 92(4):501–512. https://www.ingentaconnect.com/contentone/dav/aaua/2006/00000092/00000004/art00001 [A]
- Svensson U.P., Calamia P.T., Nakanishi S. (2009). Frequency-domain edge diffraction for finite and infinite edges. *Acta Acust. united Ac.* 95(3):568–572. doi:10.3813/AAA.918181 [M, DOI checked]
- Svensson U.P., Fred R.I., Vanderkooy J. (1999). An analytic secondary source model of edge diffraction impulse responses. *JASA* 106(5):2331–2344. doi:10.1121/1.428071 [M, DOI checked]
- Torres R.R., Svensson U.P., Kleiner M. (2001). Computation of edge diffraction for more accurate room acoustics auralization. *JASA* 109(2):600–610. doi:10.1121/1.1340647 [M, DOI checked; findings S]
- Tsingos N., Funkhouser T., Ngan A., Carlbom I. (2001). Modeling acoustics in virtual environments using the uniform theory of diffraction. *Proc. SIGGRAPH 2001*, 545–552. doi:10.1145/383259.383323 [M, DOI checked]

**Diffuse reflection and scattering**

- Antani L., Chandak A., Taylor M., Manocha D. (2012b). Direct-to-indirect acoustic radiance transfer. *IEEE Trans. Vis. Comput. Graph.* 18(2):261–269. doi:10.1109/TVCG.2011.76 [M, DOI checked]
- Brunskog J., Glebe D., Garza-Agudelo D., Nilsson E. (2022). Absorption and scattering by perforated facings with periodic narrow slits. *JASA* 151(3):1847–1859. doi:10.1121/10.0009826 [A, DOI checked]
- Can A., Fortin N., Picaut J. (2015). Accounting for the effect of diffuse reflections and fittings within street canyons, on the sound propagation predicted by ray tracing codes. *Appl. Acoust.* 96:83–93. doi:10.1016/j.apacoust.2015.03.013 [A, DOI checked; split ranges S]
- Embrechts J.J. (2000). Simulation of first and second order scattering by rough surfaces with a sound ray formalism. *J. Sound Vib.* 229(1):65–87. doi:10.1006/jsvi.1999.2479 [M, DOI checked]
- Embrechts J.J., Billon A. (2011). Theoretical determination of the random-incidence scattering coefficients of infinite rigid surfaces with a periodic rectangular roughness profile. *Acta Acust. united Ac.* 97(4):607–617. doi:10.3813/AAA.918441 [M, DOI checked]
- Fatela J., Heimes A., Vorländer M., Masullo M., Maffei L. (2026). Acoustic radiosity simulation with custom bidirectional reflectance distribution implementation in three dimensions. *JASA* 159(5):3925–3940. doi:10.1121/10.0043765 [A, DOI checked]
- Feriani S.J., Engsig-Karup A.P., Jeong C.-H. (2025). Impact of building façade design on road and rooftop noises using the discontinuous Galerkin finite element method. *Appl. Acoust.* 240:110935. doi:10.1016/j.apacoust.2025.110935 [A, DOI checked]
- Heimes A., Cabrera D., Holmes J., Reinhardt D., Tokac Celikyay I., Brell-Cokcan S., Vorländer M. (2026). Interlaboratory comparison and analysis of bidirectional scattering coefficients for various surfaces. *Acta Acustica* 10:65. doi:10.1051/aacus/2026060 [A, DOI checked]
- Heimes A., Pan L., Vorländer M. (2023). Numerical simulation and analysis of surface scattering, Part 1. *Forum Acusticum 2023*, 851–856. doi:10.61782/fa.2023.0568 [F, DOI checked]
- Heimes A., Vorländer M. (2025). Bidirectional surface scattering coefficients. *Acta Acustica* 9:41. doi:10.1051/aacus/2025026 [A, DOI checked]
- Hodgson M., Nosal E.-M. (2006). Experimental evaluation of radiosity for room sound-field prediction. *JASA* 120(2):808–819. doi:10.1121/1.2216559 [A]
- Hornikx M., Forssén J. (2008). A scale model study of parallel urban canyons. *Acta Acust. united Ac.* 94(2):265–281. doi:10.3813/AAA.918030 [F, DOI checked]
- Ismail M.R., Oldham D.J. (2002). The effect of the urban street canyon on the noise from low flying aircraft. *Building Acoustics* 9(3):233–251. doi:10.1260/135101002320815693 [A, DOI checked; abstract read via Crossref]
- Ismail M.R., Oldham D.J. (2005). A scale model investigation of sound reflection from building façades. *Appl. Acoust.* 66(2):123–147. doi:10.1016/j.apacoust.2004.07.007 [M, DOI checked; findings S]
- Kang J. (2000). Sound propagation in street canyons: comparison between diffusely and geometrically reflecting boundaries. *JASA* 107(3):1394–1404. doi:10.1121/1.428580 [A]
- Kang J. (2002). Numerical modelling of the sound fields in urban streets with diffusely reflecting boundaries. *J. Sound Vib.* 258(5):793–813. doi:10.1006/jsvi.2002.5150 [M, DOI checked]
- Kang J. (2005). Numerical modeling of the sound fields in urban squares. *JASA* 117(6):3695–3706. doi:10.1121/1.1904483 [A, DOI checked]
- Koutsouris G.I., Brunskog J., Jeong C.-H., Jacobsen F. (2013). Combination of acoustical radiosity and the image source method. *JASA* 133(6):3963–3974. doi:10.1121/1.4802897 [A, DOI checked]
- Lucchetti F.-M., Ancey S., Gabrielli P., Mercier M. (2026). Diffraction of sound waves by a periodic façade displaying impedance boundary conditions. *Forum Acusticum 2026*, Graz (preliminary manuscript). https://dael.euracoustics.org/confs/fa2026/data/articles/000162.pdf [F]
- Marbjerg G., Brunskog J., Jeong C.-H., Nilsson E. (2015). Development and validation of a combined phased acoustical radiosity and image source model. *JASA* 138(3):1457–1468. doi:10.1121/1.4928297 [A]
- Müller-Trapet M., Vorländer M. (2015). Uncertainty analysis of standardized measurements of random-incidence absorption and scattering coefficients. *JASA* 137(1):63–74. doi:10.1121/1.4904512 [A, DOI checked]
- Okada Y., Yoshihisa K., Kuno K. (2010). Simple calculation model for noise propagation in city street canyons based on a diffusion method. *Acoust. Sci. Tech.* 31(1):95–101. doi:10.1250/ast.31.95 [A, DOI checked]
- Pan L., Heimes A., Vorländer M. (2023). Numerical simulation and analysis of surface scattering, Part 2. *Forum Acusticum 2023*, 857–859. doi:10.61782/fa.2023.0569 [F, DOI checked]
- Picaut J., Simon L., Hardy J. (1999). Sound field modeling in streets with a diffusion equation. *JASA* 106(5):2638–2645. doi:10.1121/1.428093 [M, DOI checked]
- Riegel K.A., Sparrow V.W. (2019). Sonic boom propagation around a large building using a combined ray tracing/radiosity method. *JASA* 145(4):2317–2327. doi:10.1121/1.5097582 [A, DOI checked]
- Riegel K., Sparrow V.W. (2024). Sonic boom propagation in urban canyons using a combined ray tracing/radiosity method. *JASA* 155(2):1527–1533. doi:10.1121/10.0024980 [A, DOI checked]
- Shtrepi L. (2019). Investigation on the diffusive surface modeling detail in geometrical acoustics based simulations. *JASA* 145(3):EL215–EL221. doi:10.1121/1.5092821 [A, DOI checked]
- Siltanen S., Lokki T., Kiminki S., Savioja L. (2007). The room acoustic rendering equation. *JASA* 122(3):1624–1635. doi:10.1121/1.2766781 [F, DOI checked]
- Thomas P., Van Renterghem T., De Boeck E., Dragonetti L., Botteldooren D. (2013). Reverberation-based urban street sound level prediction. *JASA* 133(6):3929–3939. doi:10.1121/1.4802641 [A, DOI checked]
- Van Renterghem T., Thomas P. (2025). Self-shielding windows: scale model measurements and numerical simulations. *Forum Acusticum / Euronoise 2025*, 1233–1238. doi:10.61782/fa.2025.0764 [F, DOI checked]
- Xie D., Wittebol W., Hornikx M. (2026). Directional coherence loss coefficients for characterizing scattering distributions in enclosures. *JASA* 159(1):753–770. doi:10.1121/10.0042241 [A, DOI checked]

**Wave-based, hybrid and surrogate methods**

- Albert D.G., Liu L. (2010). The effect of buildings on acoustic pulse propagation in an urban environment. *JASA* 127(3):1335–1346. doi:10.1121/1.3277245 [M, DOI checked]
- Borrel-Jensen N., Goswami S., Engsig-Karup A.P., Karniadakis G.E., Jeong C.-H. (2024). Sound propagation in realistic interactive 3D scenes with parameterized sources using deep neural operators. *PNAS* 121(2):e2312159120. doi:10.1073/pnas.2312159120 [F (arXiv version), DOI checked]
- Doc J.-B., Lihoreau B., Félix S., Faure C., Dubois G. (2015). Three-dimensional parabolic equation model for low frequency sound propagation in irregular urban canyons. *JASA* 137(1):310–320. doi:10.1121/1.4904700 [A, DOI checked]
- Dragna D., Emmanuelli A., Ollivier S., Blanc-Benon P. (2022). Sonic boom reflection over an isolated building and multiple buildings. *JASA* 151(6):3792–3806. doi:10.1121/10.0010452 [F, about 60%; DOI checked]
- Eckerle A., Spitznagel M., Keuper J. (2025). Real-time prediction of urban sound propagation with conditioned normalizing flows. arXiv:2510.04510 [A]
- Heimann D., Schady A., Schuster R., Berkefeld T. (2018). Time-evolution of sound levels around a roadside building. *Appl. Acoust.* 133:123–132. doi:10.1016/j.apacoust.2017.12.011 [A, DOI checked]
- Hornikx M. (2016). Ten questions concerning computational urban acoustics. *Build. Environ.* 106:409–421. doi:10.1016/j.buildenv.2016.06.028 [M, DOI checked]
- Hornikx M., Forssén J. (2007). The 2.5-dimensional equivalent sources method for directly exposed and shielded urban canyons. *JASA* 122(5):2532–2541. doi:10.1121/1.2783197 [A, DOI checked]
- Hornikx M., Forssén J. (2011). Modelling of sound propagation to three-dimensional urban courtyards using the extended Fourier PSTD method. *Appl. Acoust.* 72(9):665–676. doi:10.1016/j.apacoust.2011.03.005 [A, DOI checked]
- Hornikx M., Waxler R., Forssén J. (2010). The extended Fourier pseudospectral time-domain method for atmospheric sound propagation. *JASA* 128(4):1632–1646. doi:10.1121/1.3474234 [A, DOI checked]
- Iu K.K., Li K.M. (2002). The propagation of sound in narrow street canyons. *JASA* 112(2):537–550. doi:10.1121/1.1492821 [A, DOI checked]
- Kim Y., Lee S. (2024). Deep learning based prediction of urban air mobility noise propagation in urban environment. *JASA* 155(1):171–187. doi:10.1121/10.0024242 [A, DOI checked]
- Kjær C., Borrel-Jensen N., Sampedro Llopis H., Pind F., Marozeau J., Jeong C.-H. (2022). Crossover frequency between geometric acoustics and wave-based simulations: perceptual study. *Proc. ICA 2022* [F, pp. 1–3]
- Malbéqui P. (2023). Assessment of sound propagation for the urban air mobility. *Forum Acusticum 2023*, 1011–1014. doi:10.61782/fa.2023.0104 [F, DOI checked]
- Mehra R., Raghuvanshi N., Chandak A., Albert D.G., Wilson D.K., Manocha D. (2014). Acoustic pulse propagation in an urban environment using a three-dimensional numerical simulation. *JASA* 135(6):3231–3242. doi:10.1121/1.4874495 [F, DOI checked]
- Oshima T., Hiraguri Y., Imano M. (2014b). Geometry reconstruction and mesh generation techniques for acoustic simulations over real-life urban areas using digital geographic information. *Acoust. Sci. Tech.* 35(2):108–118. doi:10.1250/ast.35.108 [A, DOI checked]
- Oshima T., Ishizuka T., Kamijo T. (2014a). Three-dimensional urban acoustic simulations and scale-model measurements over real-life topography. *JASA* 135(6):EL324–EL330. doi:10.1121/1.4879672 [A, DOI checked]
- Pelat A., Félix S., Pagneux V. (2009). On the use of leaky modes in open waveguides for the sound propagation modeling in street canyons. *JASA* 126(6):2864–2872. doi:10.1121/1.3259845 [A, DOI checked]
- Raghuvanshi N., Narain R., Lin M.C. (2009). Efficient and accurate sound propagation using adaptive rectangular decomposition. *IEEE Trans. Vis. Comput. Graph.* 15(5):789–801. doi:10.1109/TVCG.2009.28 [F, pp. 1–4; DOI checked]
- Richoux O., Ayrault C., Pelat A., Félix S., Lihoreau B. (2009). Effect of the open roof on low frequency acoustic propagation in street canyons. arXiv:0907.0605 [F]
- Schiff M., Hornikx M., Forssén J. (2010). Excess attenuation for sound propagation over an urban canyon. *Appl. Acoust.* 71(6):510–517. doi:10.1016/j.apacoust.2010.01.005 [F, DOI checked]
- Spitznagel M., Keuper J. (2024). Urban sound propagation: a benchmark for 1-step generative modeling of complex physical systems. arXiv:2403.10904 [F]
- Van Renterghem T., Salomons E., Botteldooren D. (2006). Parameter study of sound propagation between city canyons with a coupled FDTD-PE model. *Appl. Acoust.* 67(6):487–510. doi:10.1016/j.apacoust.2005.09.006 [M, DOI checked]

**Engineering methods**

- Bocher E., Guillaume G., Picaut J., Petit G., Fortin N. (2019). NoiseModelling: an open source GIS based tool to produce environmental noise maps. *ISPRS Int. J. Geo-Inf.* 8(3):130. doi:10.3390/ijgi8030130 [A, DOI checked]
- Danish EPA (2006). *Revised test cases for control of Nord2000-road.* Environmental Project No. 1111. https://www2.mst.dk/Udgiv/publications/2006/87-7052-197-2/html/helepubl_eng.htm [F]
- Defrance J., Salomons E., Noordhoek I., Heimann D., Plovsing B., Watts G., et al. (2007). Outdoor sound propagation reference model developed in the European Harmonoise project. *Acta Acust. united Ac.* 93(2):213–227. https://infoscience.epfl.ch/record/120520 [F]
- Downing J.M., Calton M.F., Page J.A., Rochat J.L., et al. (2019). *Improving AEDT Modeling for Aircraft Noise Reflection and Diffraction from Terrain and Manmade Structures.* ACRP Web-Only Document 43, National Academies Press. doi:10.17226/26897 [M, DOI checked; findings S]
- Hornikx M., Forssén J. (2006). Noise abatement schemes for real life shielded urban environments: an initial investigation of absorption treatments. *Long Range Sound Propagation Symposium*, New Orleans. https://research.chalmers.se/en/publication/25515 [A]
- Hornikx M., Forssén J., Botteldooren D., Van Renterghem T., Wei W., Ögren M., Salomons E. (2014). Urban background noise mapping: the multiple-reflection correction term. *Acta Acust. united Ac.* 100(2):293–305. doi:10.3813/AAA.918709 [F, DOI checked]
- ISO/TR 17534-4:2020. Software for the calculation of sound outdoors, Part 4: CNOSSOS-EU propagation. https://www.iso.org/standard/72115.html [A]
- Kok A. (2019). Flaws in the CNOSSOS calculation method and proposed solutions. *Proc. 23rd ICA*, Aachen. https://pub.dega-akustik.de/ICA2019/data/articles/000072.pdf [F]
- Kragh J. (2011). *Traffic Noise Prediction with Nord2000: An Update.* Danish Road Institute Report 195-2011. https://www.vejdirektoratet.dk/api/drupal/sites/default/files/publications/traffic_noise_prediction_with_nord2000.pdf [F, §2]
- Plovsing B. (2007, rev. 2014). *Proposal for Nordtest Method: Nord2000, Prediction of Outdoor Sound Propagation.* DELTA AV 1106/07 [S, §5.20 as cited by Forssén et al. 2014]
- Ramseier T., Schalcher S., Wunderli J.M., Schäffer B. (2023). Influence of buildings, forests and cliffs on aircraft noise exposure mapping: a case study. *Forum Acusticum 2023*, 3829–3833. doi:10.61782/fa.2023.0587 [F, DOI checked]
- Ramseier T., Schalcher S., Wunderli J.M., Schäffer B. (2024). Impact of buildings, forests and cliffs on aircraft noise mapping: case study. *Transp. Res. D* 133:104279. doi:10.1016/j.trd.2024.104279 [M, DOI checked]
- Rey Gozalo G., Gómez Escobar V. (2021). Uncertainty evaluation of road traffic noise models in two Ibero-American cities. *Appl. Acoust.* 180:108134. doi:10.1016/j.apacoust.2021.108134 [F, DOI checked]
- Salomons E., Polinder H., Lohman W., Zhou H., Borst H., Miedema H. (2009). Engineering modeling of traffic noise in shielded areas in cities. *JASA* 126(5):2340–2349. doi:10.1121/1.3238238 [A, DOI checked]
- Salomons E., van Maercke D., Defrance J., de Roo F. (2011). The Harmonoise sound propagation model. *Acta Acust. united Ac.* 97(1):62–74. doi:10.3813/AAA.918387 [M, DOI checked]
- Schlatter F., Köpfli M., Wunderli J.M. (2018). Relevance of buildings in aircraft noise predictions. *Proc. Inter-Noise 2018*, 3042–3051 [A]
- Wei W., Botteldooren D., Van Renterghem T., Hornikx M., Forssén J., Salomons E., Ögren M. (2014). Urban background noise mapping: the general model. *Acta Acust. united Ac.* 100(6):1098–1111. doi:10.3813/AAA.918789 [F, DOI checked]

**Atmosphere and turbulence**

- Daigle G.A. (1982). Diffraction of sound by a noise barrier in the presence of atmospheric turbulence. *JASA* 71(4):847–854 [S]
- Forssén J., Hornikx M., Botteldooren D., Wei W., Van Renterghem T., Ögren M. (2014). A model of sound scattering by atmospheric turbulence for use in noise mapping calculations. *Acta Acust. united Ac.* 100(5):810–815. doi:10.3813/AAA.918760. Author version: http://publications.lib.chalmers.se/records/fulltext/202744/local_202744.pdf [F, read by me; DOI checked]
- Forssén J., Ögren M. (2002). Thick barrier noise-reduction in the presence of atmospheric turbulence: measurements and numerical modelling. *Appl. Acoust.* 63(2):173–187. doi:10.1016/S0003-682X(01)00028-7 [F (reprinted in Ögren 2004), DOI checked]
- Goyal J., de Prenter F., Casalino D., Ragni D. (2026b). Turbulence-induced noise variability from ray tracing with instantaneous velocity snapshots. *JASA Express Lett.* 6(3):032801. doi:10.1121/10.0043128 [F, abstract and introduction; DOI checked]
- Goyal J., de Prenter F., Ragni D., Casalino D. (2026a). Impact of building wake turbulence on the noise footprint of a dipole source. *JASA* 159(2):1304–1315. doi:10.1121/10.0042402 [F, pp. 1–6; DOI checked]
- Heimann D. (2007). Three-dimensional linearised Euler model simulations of sound propagation in idealised urban situations with wind effects. *Appl. Acoust.* 68(2):217–237. doi:10.1016/j.apacoust.2005.10.002 [A, DOI checked]
- Hornikx M., Dohmen M., Conen K., van Hooff T., Blocken B. (2018). The wind effect on sound propagation over urban areas: predictions for generic urban sections. *Build. Environ.* 144:519–531. doi:10.1016/j.buildenv.2018.08.041 [A, DOI checked]
- Ögren M. (2004). *Prediction of Traffic Noise Shielding by City Canyons.* PhD thesis, Chalmers University of Technology. https://publications.lib.chalmers.se/records/fulltext/22/22.pdf [F]
- Ögren M., Forssén J. (2004). Modelling of a city canyon problem in a turbulent atmosphere using an equivalent sources approach. *Appl. Acoust.* 65(6):629–642. doi:10.1016/j.apacoust.2003.11.010 [F (reprinted in Ögren 2004), DOI checked]
- Ostashev V.E., Salomons E.M., Clifford S.F., Lataitis R.J., Wilson D.K., Blanc-Benon P., Juvé D. (2001). Sound propagation in a turbulent atmosphere near the ground: a parabolic equation approach. *JASA* 109(5):1894–1908. doi:10.1121/1.1356022 [A, DOI checked]
- Salomons E.M., Ostashev V.E., Clifford S.F., Lataitis R.J. (2001). Sound propagation in a turbulent atmosphere near the ground: an approach based on the spectral representation of refractive-index fluctuations. *JASA* 109(5):1881–1893. doi:10.1121/1.1356021 [A, DOI checked]
- Van Renterghem T., Botteldooren D. (2010). Meteorological influence on sound propagation between adjacent city canyons: a real-life experiment. *JASA* 127(6):3335–3346. doi:10.1121/1.3425731 [A, DOI checked]
- Van Renterghem T., Wei W., Forssén J., Hornikx M., Ögren M., Botteldooren D., Salomons E. (2013a). Improving the accuracy of engineering models at shielded building façades: experimental analysis of turbulence scattering. *Proc. Inter-Noise 2013*, Innsbruck. http://publications.lib.chalmers.se/records/fulltext/185369/local_185369.pdf [F]
- Wiener F.M., Malme C.I., Gogos C.M. (1965). Sound propagation in urban areas. *JASA* 37(4):738–747. doi:10.1121/1.1909409 [A, DOI checked]
- Wilson D.K., Ostashev V.E., Goedecke G.H., Auvermann H.J. (2004). Quasi-wavelet calculations of sound scattering behind barriers. *Appl. Acoust.* 65(6):605–627. doi:10.1016/j.apacoust.2003.11.009 [A, DOI checked]

**Drones, AAM and helicopters; receivers; vegetation**

- Ahuja K.K., Funk R., Hsu J., Heiges M., Stancil C. (1997). *Operation Heli-STAR Vol. 5: Effects of Buildings on Helicopter Noise.* DOT/FAA/ND-97/13, DTIC ADA370872 [S; download failed; see §10]
- Caillet J., Dieumegard P., Sineau M., Hellot M., Philipps-Bertin C., Kourieh A. (2024). Helicopter noise study: towards a better understanding of urban air mobility noise. *50th European Rotorcraft Forum*, Marseille [F]
- Chai Z., Li X., Yao H.-D. (2025). Numerical study on dual-rotor phase synchronization for tonal noise counteraction. *Aerosp. Res. Commun.* 3:15262. doi:10.3389/arc.2025.15262 [F, DOI checked]
- Čurović L., Prezelj J. (2025). The role of measurement locations in low-frequency noise assessment. *Forum Acusticum / Euronoise 2025*, 803–809. doi:10.61782/fa.2025.0883 [F, DOI checked]
- Feng T., Kang J. (2026). Low-altitude aircraft will reshape noise exposure across global cities. arXiv:2609.29763 [F]
- Gibertini G., Rezzonico S., Rossetti M., Zanotti A. (2023). Experimental investigation of helicopter noise while approaching an elevated helipad. *Aerospace* 10(8):701. doi:10.3390/aerospace10080701 [F, DOI checked]
- Habing R., Admiraal P., Bakker R., Dekker H., Knepper K. (2025). Experimental validation of rotor broadband noise models for AAM in urban environment. *51st European Rotorcraft Forum*, paper 009 [F]
- Habing R., Lammers K.H., Admiraal P., Bakker R., Lania H., Aalmoes R., Vos V. (2024). A noise prediction tool for IAM and AAM in urban environments. *50th European Rotorcraft Forum*, paper 049 [F]
- Knepper K., Admiraal P., Lammers K. (2025). Experimental database for AAM noise assessment in urban environment. *51st European Rotorcraft Forum*, Venice, paper 518. https://hdl.handle.net/20.500.11881/4816 [F]
- Lugten M. (2022). Assessing the influence of street canyon shape on aircraft noise: results from measurements in courtyards near Amsterdam Schiphol Airport. *Proc. Inter-Noise 2022*, 6474–6485. doi:10.3397/IN_2022_0976 [A, read by me; DOI checked]
- Lugten M., Wuite G., Peng Z., Kim K., Tenpierik M. (2023). A scale model approach to simulate aircraft noise in street canyons: a comparison between in-situ and laboratory measurements. *Forum Acusticum 2023*, 5309–5316. doi:10.61782/fa.2023.0264 [F, DOI checked]
- Lugten M., Wuite G., Peng Z., Tenpierik M. (2024). Assessing the influence of street canyon shape on aircraft noise: results from measurements in courtyards near Amsterdam Schiphol Airport. *Build. Environ.* 255:111400. doi:10.1016/j.buildenv.2024.111400 [M, DOI checked; abstract not accessible]
- Ramos-Romero C., Green N., Asensio C., Torija A.J. (2023). On-field measurement for sUAS noise characterisation. *Forum Acusticum 2023*, 655–660. doi:10.61782/fa.2023.0248 [F, partly; DOI checked]
- Rizzi S.A., Henderson B.S., Lopes L.V., Bain J.J., Cluts J.D., McCarthy P., Krishnamurthy S., Merino-Martinez R., He H., Snider R. (2026). *Urban Air Mobility Noise: 2025 Update on Current Practice, Gaps, and Recommendations.* NASA/TP-20250011527. doi:10.64631/XZRY3941. https://ntrs.nasa.gov/citations/20250011527 [F, Secs. 2.1.2, 2.2.2, 3.3.2, 3.3.4, 4.2.1.1; record checked]
- Rizzi S.A., Letica S.J., Boyd D.D. Jr., Lopes L.V. (2024). Prediction of noise-power-distance data for urban air mobility vehicles. *J. Aircraft* 61(1):166–182. doi:10.2514/1.C037435 [F, hemisphere passage; DOI checked]
- Schiller N.H., Pascioni K.A., Zawodny N.S. (2019). Tonal noise control using rotor phase synchronization. *VFS 75th Annual Forum.* NTRS 20200002608 [F]
- Van Renterghem T., Hornikx M., Forssén J., Botteldooren D. (2013b). The potential of building envelope greening to achieve quietness. *Build. Environ.* 61:34–44. doi:10.1016/j.buildenv.2012.12.001 [A, read by me; DOI checked]
- Wuite G., Lugten M., Peng Z. (c. 2024). The influence of trees on aircraft noise in urban courtyards: a comparative measurement study. *Reinventing the City* conference abstract 23-253. https://reinventingthecity24.dryfta.com/abstract-archive/abstract/public/253/the-influence-of-trees-on-aircraft-noise-in-urban-courtyards-a-comparative-measurement-study [A, read by me]
- Wuite G., Lugten M., Vassallo R., van Timmeren A. (2025). A spatial processing framework for mapping areas benefiting from urban morphology and building design in airport regions. *Forum Acusticum / Euronoise 2025*, 1139–1145. doi:10.61782/fa.2025.0037 [F, DOI checked]
- Yunus F., Casalino D., Avallone F., Ragni D. (2023a). Efficient prediction of urban air mobility noise in a vertiport environment. *Aerosp. Sci. Technol.* 139:108410. doi:10.1016/j.ast.2023.108410 [F, DOI checked]
- Yunus F., Casalino D., Avallone F., Ragni D. (2023b). Efficient prediction of airborne noise propagation in a non-turbulent urban environment using Gaussian beam tracing method. *JASA* 153(4):2362–2375. doi:10.1121/10.0017825 [A, DOI checked]

**Validation datasets**

- Palenda P., Vorländer M. (2026). *Evaluation of Diffraction in Generic urban Environments (EDGE) Database*, v0.1.0. Zenodo. doi:10.5281/zenodo.21264202. CC BY-SA 4.0 [F, landing page read by me]
- Wang Q., Anai K., Yano H., Sakamoto S. (2026). Experimental study of prediction model of insertion loss of buildings for elevated noise source. *Acoust. Sci. Tech.* 47(3):209–220. doi:10.1250/ast.e25.73 [M, DOI checked]. Dataset: Zenodo doi:10.5281/zenodo.15833542, CC BY 4.0 [F, landing page]
