> **Survey of 2026-10-11.** It was done without sight of the earlier
> [`Urban canyon propagation methods.md`](Urban%20canyon%20propagation%20methods.md) (2026-10-04)
> and reaches a compatible design. It adds explicit first-order facade images ahead of the
> table, line-of-sight-probability class tables and a tone test. It does not know the existing
> `--propagation_table` seam (`extra_loss`, per-band `path_gains`).
>
> Seven research passes (geometric, wave, PE/hybrid, engineering, aircraft/UAM, acceleration,
> physics/validation), each checked by a skeptic, then a synthesis, a critic and a revision.
> The per-pass findings and verdicts are in
> `research_notes/Urban propagation survey 2026-10-11/`.
>
> The staged plan built on it is the Claude Doc
> [NICEOPS: Urban Propagation Plan](https://claude.ai/code/artifact/4875c3e1-f4c9-4bae-9f51-2e9fb8f4d94a).

# Urban propagation models for NICEOPS: street canyons, courtyards and multipath from an elevated, moving rotorcraft source

*Synthesis for the NICEOPS lead, built from six research passes, each checked by a skeptic, then revised after an external review. Claims the skeptics marked fabricated or wrong have been removed or corrected. Numbers that appear only in papers whose full text nobody retrieved are flagged "(unverified)". Inline numbers refer to the reference list in §8.*

---

## 1. Executive summary: recommended architecture

- **Split the problem at the receiver, and split it by reflection order.**
  - **Above the canopy, nothing changes.** Straight, linear or stratified eigenrays, ISO 9613-1 or layered absorption, spreading and the refractive shadow treatment all stay as they are.
  - **At an urban receiver, a few paths are evaluated explicitly at runtime:** the direct path, its floor reflection (the existing spherical-wave ground model), and the **first-order facade images**. Each gets its own sphere lookup, Doppler factor, spreading, absorption and ground term.
  - **Everything else goes into a tabulated residual** G_res(band, arrival azimuth, arrival elevation [, range]): higher-order reflections, diffraction, diffuse scattering and turbulence. It is stored as a **delta relative to the existing open-ground model** and multiplies the band energies in the spectral loop.
  - **Why split it this way.** A table keyed only on arrival direction assumes every image path takes the direct path's sphere value and Doppler. That fails close in: at 30–100 m slant range from 10–30 m buildings, the images leave the rotor 10–40° apart. Rotor directivity (BVI lobes, the in-plane thickness lobe) changes by several dB over 10–20°. Overhead, the direct path can be receding while a facade image is still approaching. The source sphere is not known when the table is built, so no range axis can recover this. Explicit first-order paths remove the problem for the dominant reflections. The residual loses its directional detail anyway.
  - **Precedent.** The ACRP 02-79 "blended method" for AEDT has the same shape: building gain/loss grids relative to open terrain, precomputed by a separate model and added to the integrated prediction [10]. The residual table is the roadmap's "per-receiver urban reciprocity table". Structurally it is the plate receiver's K(f, azimuth) table generalised to the full sky.
- **Build the residual offline by reciprocity, with a hybrid solver.** Put the source at the microphone and record what leaves the urban canopy, binned by direction:
  - At low frequency, use a 3-D wave method (PSTD, ARD, FDTD, DG or LBM).
  - At high frequency, use image sources plus **edge diffraction** and diffuse facade scattering [18]. Use uniform single-edge forms (UTD or Pierce wedge), Pierce's thick barrier or second-order BTM for double edges such as thick roofs and parapets, and BTM as the exact check [35, 38, 42].
  - **Set the crossover from a cost budget, not a literature limit.** [45] and [49] reached about 450–500 Hz on 2011–2014 hardware. Other 3-D work goes further: distributed-memory FDTD over real topography [78], DG-FEM [9], LBM [32]. Time-domain cost grows as **f⁴** (f³ cells times f steps).
  - Cross-fade the two solvers in an overlap band. The overlap doubles as a validation.
- **Diffraction is not optional.** Ray and radiosity models without diffraction give silent shadow zones (0 PLdB against 94 PLdB lit [16]). NLR's ray tracer based on CNOSSOS-EU uses 2-D building footprints plus heights. It underpredicted diffraction-dominated microphones by 10–15 dB against Drone DigiCity measurements, while matching reflection-dominated ones well [5].
- **Store band energy for the residual, but check the tones first.**
  - **Tones.** Rotor noise is tonal. Main-rotor harmonics start around 15–25 Hz, BVI harmonics run to several hundred Hz, and the tail-rotor fundamental is roughly 50–120 Hz. A band average of a coherent comb is the wrong quantity for a single tone. The Doppler sweep over a pass (about 1.6 third-octaves) may justify it for SEL, but that has to be shown by a Tier-1 test (§5). If the test fails, store finer resolution below about 300 Hz.
  - **Coherence.** Sum coherently only pairs with path differences under about 1–2 m (floor, near-facade). Apply the turbulence factor exp(−σ²(1 − ρ)) to them. Sum facade flutter in energy.
  - **Turbulence floor.** In shielded directions, add a turbulence-scattering floor. Its size is uncertain:
    - The 2-D equivalent-sources study gives **2–5 dB at 1.6 kHz**, almost nothing below 500 Hz, and at most 1 dB on a traffic-spectrum A-weighted level [28, 68]. Its own author expects the method to underestimate.
    - The stronger evidence that the floor matters is long-term facade measurements [30] and the Wei et al. fitted model [21].
    - The existing PE phase screens and the `ShadowScattering` fitting approach are the tools for calibrating it.
- **Keep the windowed SEL exact.**
  - **The bound.** Per receiver, the bound is the energy sum of each explicit path's bound, each formed from `DirectionBounds` at that path's own geometry plus that path's reflection and ground maxima, plus the residual table's maximum.
  - **The open-sky limit is exactly 0 dB by construction.** With buildings removed, there are no facade images and the residual is zero, so an urban receiver reproduces today's open-ground receiver to the bit.
- **Tier the tables, because storage is the binding constraint.**
  - **Grid and size.** Resolve elevation where shielding switches on: 1.5° over 0–30° and 5° above, which gives 32 rows by 72 azimuths. At 27 bands that is 62 kB per receiver as uint8 dB in 0.25 dB steps (249 kB as float32), or **0.9 GB for a 14,400-point grid**.
  - **Wide-area grids** use per-*class* tables. Weight each class's lit and shielded tables by a **line-of-sight probability against elevation**, taken from radio-propagation statistics (ITU-R P.1410, Al-Hourani [69, 70]). Key the classes on frontal area index, H/W and the distance to the first building row [14].
  - **Validation microphones and inversion** use per-receiver tables.
  - **Prototype.** Start with the analytic ideal-canyon image lattice [19] before any solver work.
- **Do not adopt a standardised engineering method wholesale.**
  - CNOSSOS-EU's ground term is "fundamentally incorrect for a high source and low receiver" [26]. That was said of industrial sources, but the geometry carries over.
  - ISO 9613-2 and Harmonoise building terms were derived for near-ground sources.
  - **Nord2000 is the exception for high sources.** It handles arbitrary source heights with curved rays and Fresnel-zone-weighted ground and reflections, and it has been used for aircraft. Vancouver YVR chose it under ASCENT 040, and the FAA compared it with INM/AEDT [82].
  - Borrow individual pieces only: Nord2000's **Fresnel-zone weighting for finite facades**, the Pierce wedge, Heutschi-style reflections [25], and the template of a fitted canyon term [20, 21].
- **Keep inversion out of canyons.** Courtyard position-to-position spreads of 6–15 dB(A) [2] mean any microphone used to build source spheres in panam should stay in open terrain. The explicit first-order paths give panam a forward model of the dominant contamination, and the residual flags the directions it cannot model. Using any of this for de-propagation is a research risk (§7).

---

## 2. The physics that matters for an elevated source and ground or canyon receivers

### 2.1 Geometry decides almost everything

The aircraft-noise literature agrees on this. Under direct overflight, buildings have a **minor** effect. At **grazing incidence**, barrier shielding is substantial, and reflections matter mostly at receivers that are already shielded (Zurich departures, sonAIR [12]). Other findings:

- **ACRP 02-79** sized the blended method's region of influence from where shielding changes SEL by more than 0.5 dBA and from arrivals below **152 m AGL**. That gives 2,621 m along the extended centreline and 2,200 m lateral [10].
- **ACRP 02-79's summary**, as relayed by the external review and not checked here: building effects are localised, and most models showed reflection gains of up to 3 dB near a facade, falling off quickly with distance [10].
- **At the Schiphol container field lab**, aircraft orientation alone explained R² = 0.58 of measured shielding variance [2].
- **Hao and Kang** measured flyovers at 200–400 ft over 20 urban sites of 250 × 250 m [14]:
  - Site-to-site spread in mean attenuation reached about 4.6 dB at 3150 Hz within 1000 m of the track.
  - Altitude mattered little (about 2 dB).
  - Attenuation correlated mainly with **building frontal area index** and **distance from the track to the first building row**.
  - Those are the natural keys for class tables.

That is bad news for helicopters specifically. Approach, departure and urban-heliport operations at 30–300 m spend most of their time at the low elevation angles where buildings matter most.

Two sub-regimes for a receiver in a canyon or courtyard:

- **Lit (line of sight to the source through the canyon opening).** The field is the direct path, the ground reflection, and a finite set of facade images. The number of images that escape through the open top grows as elevation falls and as H/W rises. The first-order images are the ones the runtime evaluates explicitly. Reported magnitudes:
  - Gaussian beam tracing of an eVTOL among building blocks: 0–5 dB above flat terrain [6].
  - sonAIR at Meiringen: reflections of up to 3–4 dB locally on facing facades [11]. That case is an F/A-18 in an alpine valley, the metric is the annual rating level Lr, and the paper recommends *against* including buildings in standard mapping.
  - Lyon (1974), the classical canyon paper that treats sources above the street, reportedly calls aircraft noise "slightly amplified" by channelling [15] (quote unverified).
  - Feriani et al., DG-FEM with a rooftop source, validated on a 1:10 model [9]:
    - The **opposite building's reflections add up to 12 dB in a narrow canyon** at a facade receiver.
    - Facade tilting changes levels by −8.4 to +4.6 dB.
    - Balconies matter little for the rooftop source (about ±1 dB), but range from −3.5 to +6.3 dB for road sources, with low-frequency resonances.
- **Shielded (rooftop between source and receiver).** The field is rooftop diffraction, then canyon multiple reflections, then turbulence scattering. Reported magnitudes:
  - Schiphol courtyards, measured: 9.5 dB(A) of shielding for line-of-sight flight points against 5.8 dB(A) for non-line-of-sight ones, with distribution peaks at 13 and 5 dB(A), and up to **14.6 dB(A)** between the exposed and shielded sides [2].
  - The same group's 1:50 scale model reproduced relative levels to 0.5–4.2 dB(A) [3].
  - The 8–9 dB(A) spreads with straight facades and about 14 dB(A) with a slanted facade and inset come from the InterNoise 2022 version of [4]. They are not confirmed in the journal text.
  - ISO-9613/CadnaA district runs for 100 m eVTOL routes reportedly predict rises of more than 20 dB(A) in semi-enclosed courtyards that are currently quiet [13] (unverified in the text). That would be a modelling result with elevated line sources and hard ground.

ONERA's canonical image-source cases for UAM [1] (about +8–9 dB near a corner, up to +14 dB in a canyon) are often quoted. But **every case puts the source at 10–30 m beside 300 m walls**, below the roofline. They show channelling for a source inside the canyon. They do not transfer to a receiver in a canyon under a source well above the roofs.

### 2.2 Facade reflections: +6, +3 dB, and reverberant build-up

- **Facade doubling.** At a facade, coherent pressure doubling gives +6 dB for a flush microphone.
- **Decorrelation away from the wall.** Direct and reflected paths decorrelate within a band once their delay exceeds roughly the inverse bandwidth. The level then tends to the energetic +3 dB, with the familiar comb at low frequency. For a rotor tone there is no band to average over, so this changes with the source; see §2.8.
- **The floor reflection, per image path.** In a canyon the floor reflection interacts with every facade image, so the ground model has to be applied **per image path**. At runtime that means each explicit facade image carries its own spherical-wave ground term. For overhead incidence the grazing angles are steep, so the reflection is nearly specular.
- **Asphalt** is hard but not rigid, with an effective flow resistivity of about 12,000–30,000 kPa·s·m⁻² (numerically equal to cgs rayls/cm). FHWA TNM uses 20,000 for pavement.
- **Reverberant build-up** from multiple facade reflections is what standard maps miss at shielded facades, and it belongs to the residual table:
  - Wei et al. fitted a canyon multiple-reflection term to an extensive set of full-wave runs. They validated it at 9 Ghent facades, 3 of them in inner yards, with strongly improved shielded predictions [21]. The widely quoted "more than 14 dBA END underestimate" and the fit statistics are unverified.
  - Salomons et al. give an earlier engineering model for shielded city areas [22].
  - Hornikx et al. split the multiple-reflection term into source-side and receiver-side parts [20]. The receiver-side part is the one an aircraft needs, but it is valid only for separations of at least 4 street widths and was derived for traffic sources.
- **Finite facades.** A facade a few Fresnel zones wide does not reflect like an infinite plane at 30–100 m range. Nord2000's Fresnel-zone weighting [82] is a cheap, range-dependent correction to the explicit image amplitudes. Ögren's thesis [68] applies the same approach in a canyon and finds it close to the equivalent-sources solution, though smoother near resonances.

### 2.3 Diffuse versus specular facades

Diffuse and specular canyon boundaries behave very differently:

- Kang showed that diffuse (radiosity) and specular (image) canyon boundaries give **considerably different** along-street attenuation and RT30 [18].
- In a ray-plus-radiosity sonic-boom study, diffuse reflections had the largest effect on canyon signatures [16].
- Measured data are thin. ISO 17497-1 gives only a random-incidence scattering coefficient. ISO 17497-2 defines a directional diffusion coefficient, but measured facade data at steep incidence are essentially absent.

Energy methods that can build the diffuse part of the residual table:

- **Diffusion-equation street models** [24, 76].
- **Transmission-line matrix (TLM)** time-domain modelling with impedance boundaries [77].
- **Acoustic radiance transfer** [75]. This generalises radiosity to arbitrary directional reflection, which makes it the natural builder for a direction-binned diffuse table with a per-facade s(f).

Practical consequence: the builder needs a scattering coefficient s(f) per facade class. The sensitivity of **SEL** (not just RT) to s is unquantified for elevated sources. That is an early study to do with the prototype.

### 2.4 Rooftop diffraction

Rooftop diffraction is the dominant path for shielded receivers and the regime where existing urban tools fail [5]. Usable formulations:

- **Pierce's rigid wedge and thick (three-sided) barrier asymptotics** [35]. Maekawa's equivalent thin screen "may lead to sizeable errors" for thick or double edges.
- **Hadden–Pierce** for arbitrary point-source position [36].
- **UTD** [37]. For single hard wedges, frequency-domain forms are "fairly accurate even at relatively low frequency" against exact Biot–Tolstoy [42]. Cascaded UTD fails when the second edge lies in the first edge's transition zone. That is common for thick roofs and parapets under an elevated source; use the Pierce thick-barrier solution or second-order BTM there.
- **The exact BTM line-integral form** [38, 39], with fast edge subdivision [40], and the Edge Source Integral Equation (ESIE) for all higher orders on convex polyhedra [41].
- **Energy-based uncertainty diffraction (UBD)**, a ray-tracing add-on that agrees with the secondary-edge-source model to under 1 dB for most screen and slit cases [43]. Validation beyond screens is thin.

How many diffraction orders a closed courtyard under an elevated source needs is not quantified anywhere found.

### 2.5 Coherence, turbulence and motion

- **Turbulent scattering into building shadows rises with frequency, but its size is disputed:**
  - The canyon equivalent-sources model gives **2–5 dB at 1.6 kHz**, less at lower frequencies, almost nothing below 500 Hz, and at most 1 dB on an A-weighted traffic level [28, 68]. Its author expects the method to underestimate somewhat. (The indexed abstract of [28] reads "25 dB at 1.6 kHz". The author's thesis, which discusses the same paper, says "in the order of 2–5 dB", so the abstract has evidently lost an en-dash.)
  - Long-term measurements found scattering significant at shielded facades even below 1 kHz, with strong-turbulence periods dominating Leq [30]. This, with [21], is the main evidence that the floor matters.
  - The engineering scattering model has single- and double-canyon corrections [29]; its coefficients are unverified.
  - Building-wake turbulence produces localised level variability of more than 3 dB for a tonal dipole source [32]. No frequency is stated in the abstract. A companion letter appeared in JASA Express Letters 6(3), 032801 (2026).
- **Coherence between paths.**
  - **Bandwidth.** A crude heuristic says two paths decorrelate within a third-octave once the path difference exceeds about c/(0.23 f). For facade flutter in a 20 m canyon (2W = 40 m) that happens above about 37 Hz, and above about 75 Hz at 60° incidence.
  - **Short pairs.** Floor and near-facade reflections, with path differences of 1–2 m, stay coherent to about 1 kHz.
  - **Turbulence factor.** Turbulence multiplies the cross-term by exp(−σ²(1 − ρ)), where ρ is the transverse correlation over the separation between the paths. Plain exp(−σ²) overstates decorrelation for nearby pairs.
  - **Motion.** A moving source sweeps the interference pattern across the receiver, so static geometry with fine-frequency coherent sums makes interference look more persistent than it is.
  - **Rule:** sum coherently only pairs with Δr ≲ 1–2 m, which the existing ground model already does for the floor. Sum facade flutter in energy.
- **Doppler and directivity per path.**
  - **The excursion.** At M ≈ 0.18 the source Doppler factor 1/(1 − M cos θ) runs from 0.85 to 1.22 across a pass, which is larger than a third-octave.
  - **The far field.** In the far field the spread between paths at any instant is small.
  - **Close in it is not.** Each image leaves the rotor in a different direction, spread by roughly the canyon extent over R: 10–40° at 30–100 m for 10–30 m buildings. Where a facade's normal is parallel to the track, the image's velocity is reversed. So directly overhead at low altitude, the direct path can be receding while an image is still approaching.
  - **Consequence.** This is why the first-order images are evaluated explicitly, each with its own emission direction and Doppler factor. The residual can stay quasi-static.
  - **Per-path time-variant Doppler** beyond that matters for auralization or narrowband work [65].

### 2.6 Meteorology in canyons

- **Large urban wind effects** (+15 to +23 dB(A) from downward refraction over generic urban sections [33]) are reported for **ground sources propagating across many blocks**. They do not transfer to steep paths from overhead.
- **3-D street layouts.** Heimann's general finding is quieter upwind and louder downwind facades, except that upwind facades are slightly louder in parallel streets [34].
- **Measurements:**
  - A months-long courtyard-to-street experiment exists [31]. Its reported weak temperature-gradient effect, small refraction and strong, wind-dependent coherence loss were not checked against the text.
  - Wiener, Malme and Gogos measured siren propagation in Boston from both high-mounted and street sources, and saw no appreciable refraction over a block or two [80].
- **The defensible split for NICEOPS.** The mean-profile refraction stays in the existing eigenray solvers above the canopy, and the canopy is treated as still air with elevated turbulence. That is an approximation. In a moving medium, reciprocity strictly gives way to the flow-reversal theorem.

### 2.7 Low frequency

- **Canyons are leaky open waveguides at low frequency** [53, 54, 56].
- **Resonances.** Sonic-boom simulations show canyon resonances lengthen signals but change perceived level only modestly for wide streets, more as streets narrow [52].
- **Facade openings and 2-D models** [45]:
  - A facade opening can raise courtyard levels by up to 10 dB(A) relative to propagation over the roof.
  - 2-D models get mitigation measures wrong. They underpredict facade screens by up to 4.4 dB(A) and facade absorption by up to 1.5 dB(A), and overpredict roof screens by up to 1.7 dB(A).
  - 3-D beat 2-D against measurement, especially over rooftops [49].
- **Where geometric acoustics stops.** Geometric acoustics with diffraction add-ons is defensible where canyon and building dimensions are several wavelengths. λ = 3.4 m at 100 Hz, so a 10 m canyon is only about 3λ. The modal, interference-dominated bands below a few hundred Hz are where a wave solver earns its cost. Those are also the bands where rotor tones sit (§2.8).

### 2.8 Rotor tones

- **The tones.** Main-rotor harmonics begin around 15–25 Hz, BVI harmonics extend to several hundred Hz, and the tail-rotor fundamental is roughly 50–120 Hz.
- **What a tone sees.** For a tone, the canyon's coherent response at the Doppler-shifted frequency is what matters, not its band average. Over a pass, the Doppler excursion (about 1.6 third-octaves) sweeps each harmonic across the comb features, and SEL integrates over that sweep. That may make band energy adequate for SEL, but it is an assumption until tested (§5, Tier 1).
- **What NICEOPS knows.** NICEOPS spheres carry band levels, so the tonal content within a band is unknown at runtime. Even a narrowband table would need a tonal source model to use.

---

## 3. Comparison of model families

| Family | Physics captured | Frequency validity | Cost | Maturity / open source | Fit to ~1e8 evaluations |
|---|---|---|---|---|---|
| Explicit first-order images at runtime (direct + floor + facade images) | Per-path directivity, Doppler, spherical-wave ground; finite-facade Fresnel weighting [82] | All bands (specular part) | ~1 µs per path per evaluation on the spectral+ground path, so ~3–5 µs per urban evaluation | Reuses the existing sphere lookup and `ground_effect.cpp` | **Runtime core** for urban receivers |
| Analytic ideal-canyon image lattice [19] | Specular facades + ground, coherent or incoherent | ≳100–200 Hz; no diffraction, open ends or finite height | µs per receiver per band | Textbook; trivial to write | Prototype table builder, and a unit test for any tracer |
| Image sources + edge diffraction (UTD/Pierce; thick barrier or BTM for double edges; BTM/ESIE reference) [35–42] | Reflection–diffraction–reflection paths, rooftop shielding | UTD fair to low kr for single hard wedges [42]; cascaded UTD fails in transition zones; BTM exact for rigid wedges | UTD ~0.1–1 µs per path per band; BTM ms per path | Mature theory; Svensson's toolbox exists (license not checked) | **Core GA builder** for the residual |
| 3-D ray/beam tracing (incl. Gaussian beams, GPU, UBD diffraction) [6, 7, 43, 63, 64] | Many-order specular (± diffuse); beams carry refraction | ≳ a few hundred Hz; diffraction must be added | GBT vertiport 55.7 s against 3.7 s for straight rays [6]; GPU GBS stage 790× but **~88× end to end** [63] | Mature commercial (Mithra, SoundPLAN); research GBT; Pachyderm used on Schiphol courtyards [79] | Offline builder only |
| Energy models: radiosity, diffusion equation, TLM, acoustic radiance transfer [16, 18, 24, 75, 76, 77] | Diffuse facade energy, reverberation; ART adds directional reflection | High frequency, diffuse field | Seconds to minutes per geometry | I-Simpa/SPPS (GPL-3) | Diffuse part of the high-frequency residual; ART fits direction-binned output |
| Engineering canyon terms (Wei/Hornikx A_can [20, 21], Salomons [22], Schiff [23], Heutschi [25]) | Fitted multiple-reflection, shielding, roof shape | Octaves 63 Hz–4 kHz | ns | Published fits, **all for street-level sources** | Usable in the loop but needs refitting for elevated incidence |
| Line-of-sight probability statistics (radio: ITU-R P.1410, Al-Hourani, P.2108) [69–71] | P(LoS) against elevation from built-up fraction, density and height distribution | Geometry only, frequency-independent | ns | Standardised (ITU-R) | **Weights lit/shielded class tables** for wide-area grids |
| Radio roof-to-street diffraction (Walfisch–Bertoni, Ikegami) [72, 73] | Closed-form incidence-dependent rooftop-to-street loss | Derived at cm wavelengths; diffraction terms must be re-derived at acoustic kr | ns | Mature in radio | Template for an analytic class-table term |
| Standard methods (CNOSSOS, ISO 9613-2, Harmonoise) [26, 27, 10] | 2.5-D vertical-plane diffraction, image reflections | Octave/third-octave | Cheap formulas; path finding dominates | NoiseModelling (GPL-3) [27] | Cross-check only; CNOSSOS's ground term invalid for high sources |
| Nord2000 [82] | Curved rays, Fresnel-zone-weighted ground and finite reflectors, arbitrary source height | Third-octave | Cheap formulas | Used for aircraft (YVR, FAA comparison) | Borrow the Fresnel-zone weighting; cross-check |
| 3-D FDTD / PSTD / DG / ARD / LBM [9, 32, 44, 45, 48, 49, 78] | Everything, incl. modes, impedance, wind (FDTD/PSTD/LBM) | Set by budget; cost ∝ f⁴. 100×100×50 m block: FDTD (10 ppw) 1.2e10 cells at 1 kHz; PSTD (3 ppw, practical with PML and source smoothing) 4e7 at 500 Hz, 3e8 at 1 kHz, 2e10 at 4 kHz; ARD 450 Hz in minutes on a 2014 desktop [49] | openPSTD (GPL-3) [48] | **Low-frequency residual builder** and reference |
| 2.5-D ESM / BEM [46, 47, 61] | Exact Helmholtz in long canyons; plane-wave sweep reuses one factorisation | Full band in 2-D; homogeneous air | Cheap per incidence after factorisation | Research codes | Builder for long uniform streets |
| PE hybrids (FDTD-PE, PE-ESM, GFPE-Kirchhoff, 3-D ADI PE) [51, 23, 55] | Grazing over-roof propagation with wind | PE wedge ~30–35° for Padé(1,1) | Seconds to hours per geometry | Research | Low-elevation reference only; steep incidence outside the PE wedge |
| Modal / leaky-mode street models [53–55] | Low-frequency canyon modes, open-roof loss | Few modes (low frequency) | Small | Research | Niche low-frequency complement |
| ML surrogates [66, 67] | Learned field maps | Error grows with frequency and reflection order | ~real time per map | Research, 2-D, ground sources | Low priority; no bounds, poor for inversion |
| Precomputed reciprocal tables (runtime pattern) [59, 60] | Whatever the builder captured | Builder-limited | One gather per band, ~0.1 µs per evaluation (estimate) | Shipped in Microsoft Project Acoustics (closed source; latest release seen Nov 2022) | **The runtime form for the residual** |

---

## 4. Recommended approach in detail, mapped onto NICEOPS

### 4.1 Where it sits in the code

For each (state, observer), `spectral_level_at_observer` already receives a `RayPath`, or one per eigenray: the sphere-lookup direction, the range for spreading and absorption, and the ground model's geometry. For an **urban receiver**, the change is:

1. **The direct path is unchanged**, including its floor reflection. It covers spreading, per-band absorption (uniform, or `LayeredAbsorption`), the refractive shadow loss, and the existing ground term (`excess_ground_attenuation`, `GroundTable`, `PlateReceiver`). It also keeps the scattering floor of `ShadowTransition` and `ShadowScattering` where those apply.
2. **Each visible first-order facade image becomes one more `RayPath`.** The image receiver is the microphone mirrored in the facade plane, with a visibility test against the facade's extent. Each image path gets:
   - its own sphere lookup on its own launch direction;
   - its own Doppler factor;
   - its own spreading and absorption;
   - a facade reflection factor per band, with Fresnel-zone weighting for a finite facade [82];
   - its own spherical-wave ground term.
   For a receiver within about 1–2 m of the facade, that facade's image combines coherently with the direct path, as a two-ray term with facade impedance. Every other facade image adds in energy. A receiver needs only its nearest 2–4 facade planes, stored per receiver.
3. **The residual `UrbanReceiverTable::gain(band, az_arr, el_arr[, r])`** multiplies the band energies. It is defined relative to the open-ground model and covers higher-order reflections, diffraction, diffuse energy and the turbulence floor. It is indexed by the **arrival direction at the microphone**: the reversed launch direction for straight rays, the eigenray's arrival angle for refraction. The plate path already uses that quantity. With several eigenrays, each gets its own lookup and the results add in energy.
4. **Only flagged receivers pay.** Open-field observers keep today's path. Building one table per distinct receiver mirrors `GroundTable` per distinct height. The observer loop is observer-outer, so a receiver's 62 kB uint8 table stays in L2 across its state walk. An urban evaluation costs about one spectral+ground evaluation per explicit path, so 3–5× today's ~1 µs.

The urban term belongs on the **spectral path only**, as `--atmosphere` does: the broadband EAA path has no bands to weight. A broadband gain built from a reference rotor spectrum is possible, but it would need its own justification.

**Peak-search bound.** Per receiver, the bound is the energy sum of:

- each explicit path's bound, from `DirectionBounds` at that path's own image geometry plus that path's maximum reflection factor and the ground's `max_power_ratio`;
- the residual table's maximum gain.

That keeps `windowed_sound_exposure_level` exact. A check like `refraction_*_windowed_is_exact` should ship with it. Facade images will loosen the bound for facade receivers. Shielding is a loss and needs no bound.

### 4.2 What gets precomputed, keyed on what, and how big

- **Key:** receiver index (or class id) × third-octave band × arrival azimuth × arrival elevation. An optional **range axis** (3–4 nodes, e.g. 30, 60, 120 m and ∞) covers low states. With the first-order paths explicit, the residual's dependence on source directivity and range should be weaker. Nobody has tested that for 30–100 m passes near 10–30 m buildings (§7).
- **Resolution:** non-uniform in elevation, because shielding switches on at the roofline elevation arctan(H/d) and buildings matter most at low elevation.
  - **Rows:** 1.5° steps over 0–30° (20 rows) and 5° steps over 30–90° (12 rows), by 72 azimuths at 5°.
  - **Shadow edges:** keep an explicit shadow mask, or interpolate in energy, so one cell never straddles a 10–20 dB shadow step.
  - **Azimuth:** the residual is smoother in azimuth than the full field, so a 10° step may suffice, which would halve storage.
- **Size:** 2,304 cells × 27 bands = 62,208 values per receiver.
  - **Per receiver:** 62 kB as uint8 dB in 0.25 dB steps (a 64 dB range), or 249 kB as float32.
  - **Per 14,400-point grid:** 0.9 GB uint8, or 3.6 GB float32. Hence two tiers:
  - **Class tables** for wide-area grids. Classify each grid receiver from LoD1 or LoD2 geometry by H/W, offset across the street (x/W), street azimuth (stored relative to the street axis and rotated at lookup), frontal area index, distance to the first building row [14], and facade class (absorption, scattering). Store separate lit and shielded tables. Weight them per cell by P(LoS | elevation, morphology) from ITU-R P.1410 or the Al-Hourani sigmoid [69, 70]. Tens of classes means megabytes. The Schiphol group already uses elevation-angle shielding maps this way [74].
  - **Per-receiver tables** for scattered microphone sets, validation sites and inversion microphones. A 96-microphone array is about 6 MB as uint8.
- **Builder output:** band energy, with coherent fine-frequency sums only for path pairs with Δr ≲ 1–2 m, each weighted by exp(−σ²(1 − ρ)). Everything else is summed in energy. In shielded directions, add a gamma_t-scaled scattering floor in energy. Fit it with the existing PE phase-screen machinery, as `ShadowScattering` was, and check it against [28, 30].
- **Invalidation:** store a hash of the geometry patch and the builder version with each table.

### 4.3 The builder

- **Reciprocal by construction.** Source at the microphone. Record the outgoing field on a hemisphere or box above the roofs, then transform to far-field directions (Kirchhoff–Helmholtz) for wave solvers, or bin escape directions for GA.
  - **What reciprocity costs.** Reciprocal GA does **not** save traces. A forward trace from one state reaches every receiver, so forward costs about N_states and reciprocal about N_receivers, which is the larger number on the demo grid.
  - **What it buys.** The table is **trajectory-independent**, reusable across tracks, altitudes and vehicles, and directly indexable by direction.
  - **Wave solvers need it.** A forward domain would have to contain the aircraft.
- **Subtract what the runtime does explicitly.**
  - **In GA this is exact:** drop the direct, floor and first-order facade paths from the residual.
  - **In a wave solver** the total field contains them coherently, so an energy subtraction can go negative in interference nulls.
  - **Below the crossover**, the simplest choice is that the table holds the whole urban delta and the runtime's facade images are switched off. Rotor harmonics' directivity is broader at low frequency, which makes the separability error smaller there. That has to be checked by the separability test.
- **Low band.** A 3-D wave solver on a ~100 m neighbourhood, up to whatever crossover the budget allows.
  - **Cost:** PSTD at a practical 3 points per wavelength [44] is about 4e7 points at 500 Hz for a 100 × 100 × 50 m block. Each octave up multiplies the work by 16.
  - **Precedent:** ARD reached 450 Hz in minutes on a 2014 desktop [49].
  - **Licensing:** openPSTD is GPL-3 [48], so run it out of process, or write a small in-house MIT FDTD/PSTD for rigid and impedance boxes.
  - **Scale:** per receiver this is affordable for class representatives and for tens to hundreds of microphones, not for 14,400 individual grid points.
- **High band.** Write it in-house (MIT, CPU BVH), with optional GPU acceleration as a build option like `NICEOPS_PLOTTING`. It needs:
  - image sources from second order up, each with a visibility test through the canyon opening;
  - Weyl–van der Pol ground reflection per image (reuse `ground_effect.cpp`);
  - edge diffraction for roof and corner edges: UTD or Pierce for single edges, the Pierce thick barrier or second-order BTM for double edges and parapets;
  - a diffuse fraction s(f), handled by acoustic radiance transfer [75] or ray energy exchange.
- **Overlap band.** Both solvers run there. Their disagreement is a standing validation metric.

### 4.4 Composition with refraction and the shadow models

Atmosphere above the roofs, urban paths and residual below, joined at the microphone by arrival direction. Three interfaces need care:

1. **Do not double-count scattering.** If the refractive shadow floor and the building-shadow floor both apply, combine them once (energy sum or max), never stack them.
2. **The refraction near-field.** The direct path keeps its ground term, so `ShadowTransition`'s ground-term fade (`ground_weight`) applies to it unchanged. Facade image paths need a decision: either share the direct path's fade or take their own transition from their own geometry. The residual table, built in still air, takes none.
3. **In-canyon wind is neglected.** This is defensible for steep paths and questionable for distant, low-elevation states.

### 4.5 Inverse use in panam

The explicit first-order paths give panam a forward model of the dominant urban contamination, path by path, with each path's own source direction. Three things still argue against de-propagating through canyons:

- Shadow and null entries in the residual carry a large dynamic range.
- Diffuse contributions carry no directional information.
- Turbulence makes the shielded-direction level stochastic.

Recommendation: keep inversion microphones in open terrain. Compute the urban model for each flight-test microphone anyway, and **exclude source directions where the residual |G_res| exceeds about 1 dB, or where facade image paths contribute more than about 1 dB**. That turns the urban model into a data-quality filter for sphere building rather than a deconvolution operator. Modelling first-order contamination explicitly, rather than excluding it, is a later option once the separability tests pass.

### 4.6 Staged plan

1. **Prototype on the lattice.** Build the analytic ideal-canyon lattice, split into explicit first-order paths plus residual, as `UrbanReceiverTable`. Add the bound plumbing and the open-sky identity test.
2. **Tone test.** Run the Tier-1 tone-sweep test (§5) to decide band energy against finer resolution below about 300 Hz.
3. **Diffraction and diffuse energy.** Add edge diffraction (including double edges) and a diffuse fraction, then measure the SEL sensitivity to s(f).
4. **Wave low band.** Add a 3-D wave low band for 3–4 canonical geometries, then cross-fade.
5. **Real district.** Build class tables from LoD2 for a real district, with LoS-probability weighting and the turbulence floor.
6. **Range axis**, if §5's separability test demands it.

---

## 5. Full-wave reference strategy

This is the urban counterpart of what the CN-PE and the Attenborough benchmark do for refraction. Three tiers:

**Tier 0: analytic and semi-analytic (in-tree unit tests).**
- Ideal two-wall canyon image lattice, coherent and incoherent [19].
- Single rigid wedge: BTM impulse response [38, 39] against the UTD/Pierce frequency-domain coefficient [35, 37, 42].
- Thick barrier (double edge) [35], against cascaded UTD to show where the cascade fails.
- **Open-sky identity:** with buildings removed, an urban receiver reproduces the open-ground receiver to the bit.
- **Explicit-path split:** runtime explicit images plus residual must reproduce the full lattice for an omnidirectional source at any range. A directional synthetic source shows the separability gain of the split against a direction-only table.
- The GA tracer must reproduce the lattice exactly with diffraction off.

**Tier 1: 2-D cross-canyon cuts.**
- **Solver.** A small 2-D FDTD/PSTD, or openPSTD out of process, for plane-wave incidence at elevation β onto canyons with H/W = 0.5, 1, 2 and rigid or impedance facades. Compare GA plus diffraction in the overlap band.
- **Tone sweep.** Pass rotor harmonics (main-rotor, BVI, tail-rotor) through the coherent canyon response, Doppler-swept as in a real pass, and integrate SEL. Compare against the band-energy table. If the error exceeds about 1 dB, store finer resolution below about 300 Hz.
- **CN-PE for low-elevation over-roof slices only.** It is one-way and limited to about 30–35°, so it cannot carry facade backscatter. A directional source can be introduced as in [58]. Its role is to check grazing per-canyon terms before anyone applies one to helicopters, such as Schiff's excess attenuation of about 1–4 dB per intermediate canyon between ground-level canyons [23].

**Tier 2: 3-D reciprocal wave runs to the budgeted crossover** (PSTD/ARD/FDTD), on:
- a straight canyon,
- an intersection or corner,
- a closed courtyard with and without a facade opening [45],
- one LoD2 block from a candidate measurement site.

Tests to run:
- **(a) Separability.** The table plus explicit paths against point-source runs at 30, 60, 120 and 300 m, using a **directional source** with realistic lobes (a sphere from a shipped `.nod`). This sets whether a range axis is needed and whether the low-band full-delta choice in §4.3 holds.
- **(b) GA plus diffraction against wave** in the overlap band.
- **(c) Band energy against a coherent field.**
- **(d) Turbulence**, via random realisations or phase screens. Fit the scattering floor and check it against the 2–5 dB at 1.6 kHz of [28, 68] and the facade measurements of [30].

**Benchmarks to reuse:**
- The CRREL training-village pulse data, compared against 2-D FDTD [50] and 3-D ARD [49].
- The 1:40 parallel-canyon scale-model study [47].
- The 3-D FDTD and scale-model measurements over real-life topography [78].
- The Trikootam et al. 2017 urban benchmark [57]. That is only a one-page meeting abstract; the case definitions would have to be obtained from the authors before it can serve as a suite.

---

## 6. Validation data and a proposed measurement

**Existing (best first):**

- **ACRP 02-79** [10]. Its empirical-dataset survey lists a legacy **"NYC Urban Helicopter"** dataset, which would be helicopter data, urban, and already paired with model comparisons. Availability needs checking with Volpe or Blue Ridge Research. The same project's **566 measured SEL events** at LGB, El Segundo, Kittyhawk and Playa del Oro give the realistic engineering-correction floor (verified against the project presentation):
  - Measured minus AEDT was −1.9 ± 4.8 dB.
  - It became +0.9 ± 3.9 dB with TNM building offsets and +0.1 ± 4.4 dB with ISO 9613-2.
  - Both corrections **over-corrected** shielded departures: −5.8 → +2.0 dB with TNM and +2.1 dB with ISO 9613-2.
- **Schiphol container field lab** (TU Delft): radar-matched fixed-wing flyovers into three courtyards at about 37° mean elevation, with shielding statistics [2] and a validated 1:50 scale-model protocol [3, 4]. Kim et al. ran Pachyderm against these measurements [79], the closest precedent for stage 3 of the plan (title and results not retrieved).
- **NLR Drone DigiCity**: 22 microphones at 1.2 m among 24 containers, an omni source (100–400 Hz tones plus broadband) and an M600 drone at 8.5 and 20 m. It separates reflection-dominated from diffraction-dominated positions [5], which makes it ideal for the diffraction stage. Its tonal content also suits the tone test.
- **Hao and Kang**: flyovers at 200–400 ft over 20 morphologically characterised sites [14]. This is the data for choosing and checking class-table keys.
- **Heli-STAR Vol. 5, "Effects of Buildings on Helicopter Noise"** (Ahuja, Funk, Hsu, Heiges, Stancil; GTRI), DOT/FAA/ND-97/13, September 1997 [81]. Full-scale helicopter data near buildings. The aircraft (reportedly an S-76 on rooftop approaches) and the microphone count (reportedly about 100) are unverified.
- **The Polimi elevated-helipad scale model** [17].

**Proposed campaign.** No dataset combines a full-scale helicopter, real canyons and a controlled elevation angle.

- **Site.** One street canyon (H/W ≈ 1) and one courtyard, documented in LoD2 or LiDAR.
- **Microphones:**
  - an open-terrain reference microphone (1.2 m and a ground plate) to anchor the source through NICEOPS's existing open-field path;
  - canyon floor at centreline and at W/4;
  - a facade microphone, flush and at 2 m;
  - courtyard positions on the exposed and shielded sides;
  - one rooftop microphone above the canopy.
  - Record narrowband (raw time series), not only bands, so the tone test can be repeated on real data.
- **Flight matrix.** The S-76 or Bell 407, for which `.nod` databases exist. Sweep elevation angle at the canyon (overhead, about 45°, 30°, 15°) × azimuth relative to the street axis (along and across), with 2–3 repeats each for turbulence statistics. Include low passes (30–100 m slant range) to test separability directly.
- **Supporting measurements:**
  - GPS/INS tracking;
  - a roof-level sonic anemometer and a surface temperature profile, for gamma_t and the above-canopy profile;
  - a reciprocal **direct measurement of the urban transfer**: a drone-lifted or crane-lifted omni loudspeaker at known directions and 100–300 m range, recorded by the canyon microphones. This measures the table itself, independent of rotor directivity.
- **Optional.** A 1:50 scale model of the same site in an anechoic room, following [3], for dense direction sampling.

Keep this data out of sphere building, or use it only through the direction-exclusion filter in §4.5.

---

## 7. Risks and open questions

1. **Separability at low altitude.** A direction-only table is exact only as range/H → ∞. At 30–100 m, three things depend on the source's position and directivity, not just the arrival direction: the image path differences, the Fresnel zones on finite facades, and each image's emission direction and Doppler. Explicit first-order paths remove this for the dominant reflections. The residual still assumes it, and the low-band full-delta choice reintroduces it below the crossover. Tier 2 test (a) must settle both.
2. **The plane-wave ground term.** The explicit paths keep the spherical-wave ground term at runtime. A residual built from plane-wave incidence still drops spherical-wave and ground-wave terms in its higher-order paths. Those matter near grazing inside a canyon and for low states, which is another reason for the range axis.
3. **Turbulence in shadows.** Its size is uncertain:
   - The 2-D model gives 2–5 dB at 1.6 kHz [28, 68], while facade measurements show strong-turbulence periods dominating Leq [30].
   - If scattering sets shielded levels above about 1 kHz, deterministic geometry detail there is wasted and the gamma_t choice dominates.
   - NICEOPS already learned with `--scattering` that the pair-fit gamma_t was far too strong. The same calibration problem recurs here.
4. **Diffuse fraction and facade detail.**
   - s(f) for real facades at steep incidence is essentially unmeasured.
   - LoD1 against LoD2 (roof edges, slanted facades) may matter by several dB for overhead sources: about 14 dB(A) with a slanted facade in the InterNoise version of [4].
   - Effects tuned against street noise, such as sky-facing windows and balconies, may change sign for an overhead source. Feriani's balconies are ±1 dB for a rooftop source but −3.5 to +6.3 dB for a road source [9].
5. **Tones, and coherent against band energy below about 250 Hz.** Canyon modes and floor impedance interact coherently, and rotor harmonics sit in exactly those bands. Two things are open: whether Doppler sweep makes band energy adequate for SEL, and whether finer or complex storage is needed. The sphere database carries bands only, so a narrowband table would also need a tonal source model.
6. **Bound tightness.** Facade images and large residual gains may loosen the peak-search bound enough to erase the windowing saving at facade receivers. That costs performance, not correctness.
7. **Storage, classing and LoS weighting error.** Two things are unknown: whether class tables hold about 1 dB against per-receiver tables, and whether the LoS-probability weighting holds for SEL over a whole pass. Consecutive states see correlated obstruction, which a per-elevation probability ignores.
8. **Runtime cost of explicit paths.** 3–5× the spectral+ground evaluation cost at urban receivers. That is acceptable for microphone sets, but a fully urban 14,400-point grid would pay it everywhere.
9. **Licensing.** openPSTD, NoiseModelling and I-Simpa are GPL-3. Run them out of process only; never link them into MIT NICEOPS.
10. **In-canyon wind and reciprocity.** The still-air canyon assumption rests on short canyon-to-canyon experiments with ground sources [31], plus the 1965 Boston siren work [80]. Applying it to steep helicopter paths is an extrapolation.
11. **Literature gaps to close before citing:**
    - Hornikx's "Ten questions concerning computational urban acoustics" (Building and Environment, ~2016).
    - The Pierce & Hadden (1978) finite-impedance wedge.
    - Casalino et al.'s LBM urban UAM work beyond [32].
    - NASA RVLT urban propagation.
    - Walfisch–Bertoni and Ikegami roof-to-street diffraction [72, 73], re-derived at acoustic kr.
    - Full texts for the unverified items: the Wei et al. 2014 statistics [21], the Lugten 2022/2024 spreads [4], Feng & Kang's courtyard figure [13], Lyon's quote [15], the Heli-STAR details [81], Kim et al. [79], and the ACRP "up to 3 dB near a facade" summary [10].

---

## 8. References

*Aircraft, rotorcraft, UAM and urban measurements*
1. Malbéqui, P. (2023). Assessment of sound propagation for the urban air mobility. Forum Acusticum 2023, Turin, pp. 1011–1014. https://doi.org/10.61782/fa.2023.0104
2. Wuite, G., Lugten, M., Tenpierik, M. (2023). Determining causes of variance in ground-level aircraft noise. Forum Acusticum 2023. https://doi.org/10.61782/fa.2023.0438
3. Lugten, M., Wuite, G., Peng, Z., Kim, K., Tenpierik, M. (2023). A scale model approach to simulate aircraft noise in street canyons. Forum Acusticum 2023. https://doi.org/10.61782/fa.2023.0264
4. Lugten, M., Wuite, G., Peng, Z., Tenpierik, M. (2024). Assessing the influence of street canyon shape on aircraft noise: results from measurements in courtyards near Amsterdam Schiphol Airport. Building and Environment 255, 111400. https://doi.org/10.1016/j.buildenv.2024.111400
5. Knepper, K., Admiraal, P., Lammers, K. (2025). Experimental database for AAM noise assessment in urban environment. 51st European Rotorcraft Forum. https://dspace-erf.nlr.nl/items/fb7cf7fd-6323-4a2e-a628-2deda217add0
6. Yunus, F., Casalino, D., Avallone, F., Ragni, D. (2023). Efficient prediction of urban air mobility noise in a vertiport environment. Aerospace Sci. Technol. 139, 108410. https://doi.org/10.1016/j.ast.2023.108410
7. Yunus, F. et al. (2023). Efficient prediction of airborne noise propagation in a non-turbulent urban environment using Gaussian beam tracing. J. Acoust. Soc. Am. 153(4), 2362–2375. https://doi.org/10.1121/10.0017825
8. Bian, H., Tan, Q., Zhong, S., Zhang, X. (2021). Assessment of UAM and drone noise impact on the environment based on virtual flights. Aerospace Sci. Technol. 118, 106996. https://doi.org/10.1016/j.ast.2021.106996
9. Feriani, S. J., Engsig-Karup, A. P., Jeong, C.-H. (2025). Impact of building façade design on road and rooftop noises using the discontinuous Galerkin finite element method. Applied Acoustics 240, 110935 (volume and DOI not confirmed). https://doi.org/10.1016/j.apacoust.2025.110935 ; https://orbit.dtu.dk/en/publications/impact-of-building-fa%C3%A7ade-design-on-road-and-rooftop-noises-using/
10. National Academies (2019). Improving AEDT Modeling for Aircraft Noise Reflection and Diffraction from Terrain and Manmade Structures. ACRP Web-Only Document 43. https://doi.org/10.17226/26897 ; Downing, Calton, Page, Rochat, ACRP 02-79 presentation (source of the verified statistics, dataset list and region-of-influence criteria): https://rosap.ntl.bts.gov/view/dot/57382/dot_57382_DS1.pdf
11. Ramseier, T., Schalcher, S., Wunderli, J. M., Schäffer, B. (2023). Influence of buildings, forests and cliffs on aircraft noise exposure mapping. Forum Acusticum 2023. https://doi.org/10.61782/fa.2023.0587
12. Schlatter, F., Köpfli, M., Wunderli, J. M. (2018). Relevance of buildings in aircraft noise predictions. Proc. Inter-Noise 2018, 3042–3051. https://ince.publisher.ingentaconnect.com/content/ince/incecp/2018/00000258/00000004/art00007
13. Feng, T., Kang, J. (2026). Low-altitude aircraft will reshape noise exposure across global cities. arXiv:2609.29763. https://arxiv.org/abs/2609.29763
14. Hao, Y., Kang, J. (2014). Influence of mesoscale urban morphology on the spatial noise attenuation of flyover aircrafts. Applied Acoustics 84, 73–82. https://doi.org/10.1016/j.apacoust.2013.12.001 ; https://discovery.ucl.ac.uk/id/eprint/10058528/
15. Lyon, R. H. (1974). Role of multiple reflections and reverberation in urban noise propagation. J. Acoust. Soc. Am. 55(3), 493–503. https://trid.trb.org/View/138963
16. Riegel, K., Sparrow, V. W. (2024). Sonic boom propagation in urban canyons using a combined ray tracing/radiosity method. J. Acoust. Soc. Am. 155(2), 1527–1533. https://doi.org/10.1121/10.0024980
17. Gibertini, G., Rezzonico, S., Rossetti, M., Zanotti, A. (2023). Experimental investigation of helicopter noise while approaching an elevated helipad. Aerospace 10(8), 701. https://doi.org/10.3390/aerospace10080701

*Canyon physics and engineering terms*

18. Kang, J. (2000). Sound propagation in street canyons: comparison between diffusely and geometrically reflecting boundaries. J. Acoust. Soc. Am. 107(3), 1394–1404. https://doi.org/10.1121/1.428580
19. Li, K. M., Lai, C. Y. C. (2009). A note on noise propagation in street canyons. J. Acoust. Soc. Am. 126(2), 644–655. https://doi.org/10.1121/1.3158599
20. Hornikx, M., Forssén, J., Botteldooren, D., Van Renterghem, T., Wei, W., Ögren, M., Salomons, E. (2014). Urban background noise mapping: the multiple-reflection correction term. Acta Acust. united Acust. 100(2), 293–305. https://doi.org/10.3813/AAA.918709
21. Wei, W., Botteldooren, D., Van Renterghem, T., Hornikx, M., Forssén, J., Salomons, E., Ögren, M. (2014). Urban background noise mapping: the general model. Acta Acust. united Acust. 100(6), 1098–1111; erratum 101, 204 (2015). https://doi.org/10.3813/AAA.918789
22. Salomons, E., Polinder, H., Lohman, W., Zhou, H., Borst, H., Miedema, H. (2009). Engineering modeling of traffic noise in shielded areas in cities. J. Acoust. Soc. Am. 126(5), 2340–2349. https://doi.org/10.1121/1.3238238
23. Schiff, M., Hornikx, M., Forssén, J. (2010). Excess attenuation for sound propagation over an urban canyon. Applied Acoustics 71(6), 510–517. https://doi.org/10.1016/j.apacoust.2010.01.005
24. Okada, Y., Yoshihisa, K., Kuno, K. (2010). Simple calculation model for noise propagation in city street canyons based on a diffusion method. Acoust. Sci. Tech. 31(1), 95. https://www.jstage.jst.go.jp/article/ast/31/1/31_1_95/_article
25. Heutschi, K. (1995). A simple method to evaluate the increase of traffic noise emission level due to buildings, for a long straight street. Applied Acoustics 44(3), 259–274. https://trid.trb.org/View/425170
26. Kok, A. (2019). Flaws in the Cnossos calculation method and proposed solutions. Proc. 23rd ICA, Aachen. https://pub.dega-akustik.de/ICA2019/data/articles/000072.pdf
27. Bocher, E., Guillaume, G., Picaut, J., Petit, G., Fortin, N. (2019). NoiseModelling: an open source GIS based tool to produce environmental noise maps. ISPRS Int. J. Geo-Inf. 8(3), 130. https://doi.org/10.3390/ijgi8030130

*Turbulence and meteorology*

28. Ögren, M., Forssén, J. (2004). Modelling of a city canyon problem in a turbulent atmosphere using an equivalent sources approach. Applied Acoustics 65(6), 629–642. https://doi.org/10.1016/j.apacoust.2003.11.010 ; https://research.chalmers.se/en/publication/305 (the indexed abstract's "25 dB" is "2–5 dB" per [68])
29. Forssén, J., Hornikx, M., Botteldooren, D., Wei, W., Van Renterghem, T., Ögren, M. (2014). A model of sound scattering by atmospheric turbulence for use in noise mapping calculations. Acta Acust. united Acust. 100(5), 810–815. https://doi.org/10.3813/AAA.918760
30. Van Renterghem, T., Wei, W., Forssén, J., Hornikx, M., Ögren, M., Botteldooren, D., Salomons, E. (2013). Improving the accuracy of engineering models at shielded building facades: experimental analysis of turbulence scattering. Inter-Noise 2013. https://research.chalmers.se/en/publication/185369
31. Van Renterghem, T., Botteldooren, D. (2010). Meteorological influence on sound propagation between adjacent city canyons: a real-life experiment. J. Acoust. Soc. Am. 127(6), 3335–3346. https://doi.org/10.1121/1.3425731
32. Goyal, J., de Prenter, F., Ragni, D., Casalino, D. (2026). Impact of building wake turbulence on the noise footprint of a dipole source. J. Acoust. Soc. Am. 159, 1304–1315. https://doi.org/10.1121/10.0042402 ; companion: JASA Express Letters 6(3), 032801 (2026)
33. Hornikx, M., Dohmen, M., Conen, K., van Hooff, T., Blocken, B. (2018). The wind effect on sound propagation over urban areas: predictions for generic urban sections. Building and Environment 144, 519–531. https://researchportal.hw.ac.uk/en/publications/the-wind-effect-on-sound-propagation-over-urban-areas-predictions/
34. Heimann, D. (2007). Three-dimensional linearised Euler model simulations of sound propagation in idealised urban situations with wind effects. Applied Acoustics 68, 217–237. https://doi.org/10.1016/j.apacoust.2005.10.002

*Edge diffraction*

35. Pierce, A. D. (1974). Diffraction of sound around corners and over wide barriers. J. Acoust. Soc. Am. 55(5), 941–955. https://doi.org/10.1121/1.1914668
36. Hadden, W. J., Pierce, A. D. (1981). Sound diffraction around screens and wedges for arbitrary point source locations. J. Acoust. Soc. Am. 69(5), 1266–1276. https://doi.org/10.1121/1.385809
37. Kouyoumjian, R. G., Pathak, P. H. (1974). A uniform geometrical theory of diffraction for an edge in a perfectly conducting surface. Proc. IEEE 62(11), 1448–1461.
38. Svensson, U. P., Fred, R. I., Vanderkooy, J. (1999). An analytic secondary source model of edge diffraction impulse responses. J. Acoust. Soc. Am. 106, 2331–2344.
39. Svensson, U. P., Calamia, P. T. (2006). Edge-diffraction impulse responses near specular-zone and shadow-zone boundaries. Acta Acust. united Acust. 92(4). https://www.ingentaconnect.com/contentone/dav/aaua/2006/00000092/00000004/art00001
40. Calamia, P. T., Svensson, U. P. (2007). Fast time-domain edge-diffraction calculations for interactive acoustic simulations. EURASIP J. Adv. Signal Process. 63560. https://doi.org/10.1155/2007/63560
41. Asheim, A., Svensson, U. P. (2013). An integral equation formulation for the diffraction from convex plates and polyhedra. J. Acoust. Soc. Am. 133, 3681–3691. https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/2592097
42. Ouis, D. (2003). Noise attenuation by a hard wedge-shaped barrier. J. Sound Vib. 262, 347–364. https://doi.org/10.1016/S0022-460X(02)01134-3
43. Stephenson, U. M. (2010). An energetic approach for the simulation of diffraction within ray tracing based on the uncertainty relation. Acta Acust. united Acust. 96(3), 516–535.

*Wave-based and hybrid methods*

44. Hornikx, M., Waxler, R., Forssén, J. (2010). The extended Fourier pseudospectral time-domain method for atmospheric sound propagation. J. Acoust. Soc. Am. 128(4), 1632–1646. https://doi.org/10.1121/1.3474234
45. Hornikx, M., Forssén, J. (2011). Modelling of sound propagation to three-dimensional urban courtyards using the extended Fourier PSTD method. Applied Acoustics 72(9), 665–676. https://doi.org/10.1016/j.apacoust.2011.03.005
46. Hornikx, M., Forssén, J. (2007). The 2.5-dimensional equivalent sources method for directly exposed and shielded urban canyons. J. Acoust. Soc. Am. 122(5), 2532–2541. https://doi.org/10.1121/1.2783197
47. Hornikx, M., Forssén, J. (2008). A scale model study of parallel urban canyons. Acta Acust. united Acust. 94(4), 641. https://doi.org/10.3813/AAA.918074
48. Hornikx, M., Krijnen, T., van Harten, L. (2016). openPSTD: the open source pseudospectral time-domain method for acoustic propagation. Comput. Phys. Commun. https://doi.org/10.1016/j.cpc.2016.02.029 ; https://github.com/openPSTD/openPSTD
49. Mehra, R., Raghuvanshi, N., Chandak, A., Albert, D. G., Wilson, D. K., Manocha, D. (2014). Acoustic pulse propagation in an urban environment using a three-dimensional numerical simulation. J. Acoust. Soc. Am. 135(6), 3231–3242. https://doi.org/10.1121/1.4874495
50. Albert, D. G., Liu, L. (2010). The effect of buildings on acoustic pulse propagation in an urban environment. J. Acoust. Soc. Am. 127(3), 1335–1346. https://doi.org/10.1121/1.3277245
51. Van Renterghem, T., Salomons, E., Botteldooren, D. (2006). Parameter study of sound propagation between city canyons with a coupled FDTD-PE model. Applied Acoustics 67(6), 487–510. https://doi.org/10.1016/j.apacoust.2005.09.006
52. Dragna, D., Emmanuelli, A., Ollivier, S., Blanc-Benon, P. (2022). Sonic boom reflection over an isolated building and multiple buildings. J. Acoust. Soc. Am. 151(6), 3792–3806. https://doi.org/10.1121/10.0010452
53. Pelat, A., Félix, S., Pagneux, V. (2009). On the use of leaky modes in open waveguides for the sound propagation modeling in street canyons. J. Acoust. Soc. Am. 126(6), 2864–2872. https://doi.org/10.1121/1.3259845
54. Pelat, A., Félix, S., Pagneux, V. (2011). A coupled modal-finite element method for the wave propagation modeling in irregular open waveguides. J. Acoust. Soc. Am. 129(3), 1240–1249. https://doi.org/10.1121/1.3531928
55. Doc, J.-B., Lihoreau, B., Félix, S., Faure, C., Dubois, G. (2015). Three-dimensional parabolic equation model for low frequency sound propagation in irregular urban canyons. J. Acoust. Soc. Am. 137(1), 310–320. https://doi.org/10.1121/1.4904700
56. Richoux, O., Ayrault, C., Pelat, A., Félix, S., Lihoreau, B. (2009). Effect of the open roof on low frequency acoustic propagation in street canyons. arXiv:0907.0605. https://arxiv.org/abs/0907.0605
57. Trikootam, S. C., Hornikx, M., Fernandes, G. C., van Hooff, T., Blocken, B., Van Renterghem, T. (2017). Benchmark cases for urban sound propagation (meeting abstract). J. Acoust. Soc. Am. 141(5), 4035. https://doi.org/10.1121/1.4989306
58. Rosenbaum, J. E., Atchley, A. A., Sparrow, V. W. (2016). Inclusion of elevation angle source directivity in the parabolic equation. Noise Control Eng. J. 64(2), 181–194. https://doi.org/10.3397/1/376370

*Precomputation, acceleration and moving sources*

59. Raghuvanshi, N., Snyder, J. (2014). Parametric wave field coding for precomputed sound propagation. ACM Trans. Graph. 33(4). https://doi.org/10.1145/2601097.2601184
60. Raghuvanshi, N., Snyder, J. (2018). Parametric directional coding for precomputed sound propagation. ACM Trans. Graph. 37(4). https://doi.org/10.1145/3197517.3201339
61. Mehra, R., Raghuvanshi, N., Antani, L., Chandak, A., Curtis, S., Manocha, D. (2013). Wave-based sound propagation in large open scenes using an equivalent source formulation. ACM Trans. Graph. 32(2). https://www.microsoft.com/en-us/research/publication/wave-based-sound-propagation-large-open-scenes-using-equivalent-source-formulation/
62. Ben-Hur, Z., Alon, D. L., Mehra, R., Rafaely, B. (2019). Efficient representation and sparse sampling of head-related transfer functions using phase correction based on ear alignment. IEEE/ACM TASLP. https://cris.bgu.ac.il/en/publications/efficient-representation-and-sparse-sampling-of-head-related-tran/
63. Sheng, Z., Duan, L., Jiang, H. (2025). Accelerating Gaussian beam tracing method with dynamic parallelism on GPUs. arXiv:2501.13382. https://arxiv.org/abs/2501.13382
64. Gbikpi-Benissan, G., Magoulès, F. Beam-tracing domain decomposition method for urban acoustic pollution. arXiv:1912.03998. https://arxiv.org/abs/1912.03998
65. Schäfer, P., Reich, L., Vorländer, M. (2022). Linking atmospheric and urban auralization models. Acta Acustica 6, 28. https://doi.org/10.1051/aacus/2022021
66. Spitznagel, M., Keuper, J. (2024). Urban Sound Propagation: a benchmark for 1-step generative modeling of complex physical systems. arXiv:2403.10904. https://arxiv.org/abs/2403.10904
67. Eckerle, A., Spitznagel, M., Keuper, J. (2025). Real-time prediction of urban sound propagation with conditioned normalizing flows. arXiv:2510.04510. https://arxiv.org/abs/2510.04510

*Added in revision*

68. Ögren, M. (2004). Prediction of traffic noise shielding by city canyons. PhD thesis, Chalmers University of Technology, Report F04-01. https://publications.lib.chalmers.se/records/fulltext/22/22.pdf (Ch. 5: turbulence raises the level "in the order of 2–5 dB" at 1.6 kHz, almost nothing below 500 Hz, at most 1 dB A-weighted for a traffic spectrum)
69. ITU-R Recommendation P.1410. Propagation data and prediction methods required for the design of terrestrial broadband radio access systems (built-up parameters α, β, γ; line-of-sight probability).
70. Al-Hourani, A., Kandeepan, S., Lardner, S. (2014). Optimal LAP altitude for maximum coverage. IEEE Wireless Commun. Lett. 3(6), 569–572.
71. ITU-R Recommendation P.2108. Prediction of clutter loss; background in ITU-R Report P.2402-1 (2026). https://www.itu.int/dms_pub/itu-r/opb/rep/R-REP-P.2402-1-2026-PDF-E.pdf
72. Walfisch, J., Bertoni, H. L. (1988). A theoretical model of UHF propagation in urban environments. IEEE Trans. Antennas Propag. 36(12), 1788–1796.
73. Ikegami, F., Yoshida, S., Takeuchi, T., Umehira, M. (1984). Propagation factors controlling mean field strength on urban streets. IEEE Trans. Antennas Propag. 32(8), 822–829.
74. Wuite, G., Lugten, M. et al. (2025). Elevation-angle shielding maps for aircraft noise (title not retrieved). Forum Acusticum 2025. https://doi.org/10.61782/fa.2025.0037
75. Siltanen, S., Lokki, T., Kiminki, S., Savioja, L. (2007). The room acoustic rendering equation. J. Acoust. Soc. Am. 122(3), 1624–1635.
76. Picaut, J. (2002). Numerical modeling of urban sound fields by a diffusion process. Applied Acoustics 63(9), 965–991.
77. Guillaume, G., Picaut, J., Dutilleux, G., Gauvreau, B. (2011). Time-domain impedance formulation for transmission line matrix modelling of outdoor sound propagation. J. Sound Vib. 330(26), 6467–6481.
78. Oshima, T. et al. (2014). Three-dimensional urban acoustic simulations and scale-model measurements over real-life topography. J. Acoust. Soc. Am. (bibliographic details not retrieved). https://kaken.nii.ac.jp/grant/KAKENHI-PROJECT-23360255
79. Kim, K., Lugten, M., Wuite, G., Peng, Z., Tenpierik, M. (2023). Geometric-acoustics (Pachyderm) simulation of aircraft noise in Schiphol courtyards (title not retrieved). Proc. Building Simulation 2023 (IBPSA), pp. 1763–1770.
80. Wiener, F. M., Malme, C. I., Gogos, C. M. (1965). Sound propagation in urban areas. J. Acoust. Soc. Am. 37(4), 738–747.
81. Ahuja, K. K., Funk, R., Hsu, J., Heiges, M., Stancil, C. (1997). Heli-STAR Vol. 5: Effects of Buildings on Helicopter Noise. DOT/FAA/ND-97/13. https://rosap.ntl.bts.gov/view/dot/42032
82. Nord2000 and its aircraft use: FORCE Technology, Nord2000 validation of the propagation model, AV 1117/06, https://forcetechnology.com/-/media/force-technology-media/pdf-files/projects/nord2000/nord2000-validation-of-the-propagation-model-av-1117-06.pdf ; ASCENT Project 040 annual report (2017), https://ascent.aero/documents/2018/06/ascent-040-2017-annual-report.pdf/ ; Ahearn, M. et al. (2012), FAA comparison of Nord2000 with INM/AEDT, https://rosap.ntl.bts.gov/view/dot/9687

*[8] and [62] are kept in the list for stable numbering and are no longer cited in the text.*

---

## Changes from review

- **Turbulence magnitude (A1):** I checked this against the sources. The indexed abstract of [28] does read "25 dB at 1.6 kHz". The author's own PhD thesis (added as [68]) says "in the order of 2–5 dB" at 1.6 kHz and at most 1 dB A-weighted, so the abstract has evidently lost an en-dash. I corrected the figure in §1, §2.5, §4.2, §5(d) and §7.3. The case for a turbulence floor now rests mainly on [30] and [21].
- **Separability (A2, E):** I restructured the recommendation by reflection order. The direct path, the floor reflection and the first-order facade images are now evaluated explicitly at runtime, each with its own sphere lookup, Doppler and ground term. Only the higher-order, diffracted, diffuse and turbulence energy is tabulated. This changed §1, §2.5, §4.1–4.6, the §5 tests and §7.1. I also added a design choice for the low band, where a wave-solver residual cannot be cleanly subtracted.
- **Open-sky target (A3):** The table is now a delta relative to the open-ground model, with the spherical-wave ground term kept at runtime. The open-sky limit is exactly 0 dB by construction, and the Tier-0 test is now an identity.
- **Tones (A4):** I added §2.8, a tone-sweep test in Tier 1, a stage in the plan, narrowband recording in the campaign, and §7.5.
- **Elevation grid and storage (A5):** The grid is now non-uniform (1.5° below 30°, 5° above) with a shadow mask, stored as uint8 in 0.25 dB steps. Storage is recomputed: 62 kB per receiver, 0.9 GB per 14,400-point grid, 6 MB for 96 microphones. I removed the HRTF analogy (C17).
- **Wave-method limits (A6):** The solver crossover is now set by a cost budget. I added f⁴ scaling, PSTD at a practical 3 points per wavelength with recomputed cell counts, and other 3-D work [78, 9, 32].
- **Nord2000 (A7):** It is now separated from CNOSSOS as valid for high sources and used for aircraft. Its Fresnel-zone weighting for finite facades is adopted for the explicit images [82].
- **New model families (B8–B12):**
  - line-of-sight-probability weighting of class tables (ITU-R P.1410, Al-Hourani, P.2108) and Walfisch–Bertoni/Ikegami [69–74];
  - energy models: diffusion, TLM, acoustic radiance transfer [75–77];
  - Hao & Kang now cited, and used for the class keys [14];
  - Kim et al.'s Pachyderm validation [79];
  - Wiener et al. 1965 [80].
- **Second-order physics (C13–C18):**
  - The coherence factor is now exp(−σ²(1 − ρ)).
  - Sums are coherent only for pairs with Δr ≲ 1–2 m.
  - Added the cascaded-UTD caveat, with the thick barrier or second-order BTM for double edges, and added Kouyoumjian–Pathak as [37].
  - Asphalt flow resistivity is now in kPa·s·m⁻², with TNM's 20,000.
  - The bound now includes each explicit path's own bound.
- **Citations (D19–D24):**
  - Removed the use of placeholder [25]. Heutschi is now [25], split out of [20].
  - Moved [43] out of the radiosity row.
  - Cited [22], [24] and [58]. [8] and [62] are noted as uncited.
  - Added pages and the erratum to [21].
  - Added Feriani's road-source balcony range, and flagged its volume and DOI.
  - Heli-STAR is now cited by authors and title, with the aircraft and microphone count flagged [81].
  - Named NLR's CNOSSOS-based ray tracer in [5].
  - Reworded Schiff as "excess attenuation of about 1–4 dB".
  - Flagged as unverified: [32] (no frequency stated), [13], [15], and Project Acoustics' status.
- **Where I did not follow the review:** It asked for ACRP's 152 m threshold and statistics to be marked unverified. I extracted the ACRP 02-79 presentation text and confirmed the 566 events, the −1.9 ± 4.8 / +0.9 ± 3.9 / +0.1 ± 4.4 dB statistics, the −5.8 → +2.0 dB shielded-departure over-correction, the "NYC Urban Helicopter" dataset and the 152 m AGL criterion, so they stay as verified. I reworded the 152 m claim to match the slide: it sizes the blended method's region of influence. The review's "up to 3 dB near a facade" summary was not found in the slides and is attributed to the review.
