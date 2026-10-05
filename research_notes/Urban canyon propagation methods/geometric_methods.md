# Geometric (ray-based) methods for urban street-canyon and building propagation, for an elevated rotorcraft source

Scope: image-source and image-edge methods with edge diffraction, ray, beam and Gaussian-beam tracing, how facade scattering and coherence are handled, coupling to curved-ray atmospheres, cost, and what the virtual-acoustics groups use. The target is NICE-OPS: 10^8 to 10^9 source-receiver pairs per run, energy per third-octave band, about 10 Hz to 10 kHz. Its current per-pair budget is about 19 ns broadband, 150 ns spectral and about 1 µs with the ground model ([NICE-OPS CLAUDE.md](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md)).

Source note: several primary papers were paywalled or blocked by the publisher (403). These are Palenda & Vorländer 2026, Wei, Van Renterghem & Botteldooren 2015 (JASA), Kang 2000, Casalino et al. 2011, and the Yu et al. 2019 MDPI paper. For these the notes use the abstract or repository text only, and say so where they do. Full texts were read for the following:
- Schäfer, Reich & Vorländer 2022
- Schissler, Mehra & Manocha 2014
- Schissler, Mückl & Calamia 2021
- Schissler & Manocha 2016
- Van Renterghem, Salomons & Botteldooren 2006
- Wei et al. 2014 (Acta Acustica)
- Wei, Cooke & Havie-Clark 2015 (IOA)
- Stephenson (DAGA 2012)
- Yunus et al. 2023
- Ramseier et al. (Empa, Forum Acusticum 2023)
- Antani et al. (UNC TR10-011)
- Sheng et al. 2025 (arXiv)

---

## 1. Image-source method with edge diffraction: BTM, UTD, Pierce, multi-edge approximations, Maekawa. What diffraction order is needed, and what are the published errors?

### Takeaway
Image sources plus edge diffraction is the most validated deterministic geometric family for buildings.
- **Accuracy:** BTM is the accurate reference. UTD, the cheap asymptotic, is within about 1 dB of BTM in the lit zone and about 3-4 dB at first- and second-order diffracted receivers, and it overestimates high-frequency attenuation. Energy-based UTD at high order overestimates the diffracted field, by up to 20 dB at 84 Hz when every path is summed.
- **Orders needed:** first-order diffraction over the roof edge plus a few facade reflections dominates the shielded canyon. Wei et al. 2015 (IOA) found that levels saturate after about 4 reflections when horizontal reflections dominate. In the work cited here, second-order diffraction mattered mainly for receivers that are themselves shadowed.
- **Cost:** a full image-edge path search costs milliseconds to about a second per source-receiver pair. That is 10^3 to 10^7 times the NICE-OPS budget.

### Cited Findings
**BTM and UTD in general**
- BTM (Biot-Tolstoy-Medwin) "is the most accurate, but it suffers from high computational complexity and hence is difficult to apply in real time". VDaT (volumetric diffraction and transmission) estimates the BTM amplitude response and path length by spatial sampling. A GPU version handles "thousands of direct and specular reflection path segments in small-to-medium-size scenes" in real time — [Pisha et al. 2020, via SonicField library](https://labs.sonicfield.org/library/approximate-diffraction-modeling-for-real-time-sound-propagation-simulation)
- UTD assumes an infinite edge and "works well for models with very large edges in outdoor scenes". BTM "has been shown to be more accurate than UTD at low frequencies, where diffraction effects play an important role" — [Antani et al., Fast Geometric Sound Propagation with Finite-Edge Diffraction, UNC TR10-011](https://gamma.cs.unc.edu/BTM/BTM.pdf)

**UTD against BTM and FDTD**
- UTD against offline BTM for a rectangular obstacle (12 edges, one source): average errors were 0.67 dB (direct-sound receiver), 3.10 dB (first-order diffraction receiver) and 3.61 dB (second-order diffraction receiver). "The UTD diffraction model overestimates the attenuation of high frequencies" — Schissler, Mehra & Manocha 2014, "High-Order Diffraction and Diffuse Reflections for Interactive Sound Propagation in Large Environments", ACM TOG 33(4) 39, read in full text ([UNC GAMMA group page](https://gamma.umd.edu/subjects/sound-synthesis/))
- Diffraction-path cost in a graph-based search grows as O(n^d), with n the visible neighbours per edge and d the maximum order. The 2014 method reached 5th-6th order diffraction "at real-time rates" in its benchmarks, but per-frame cost "varies considerably with the source and listener positions" — Schissler et al. 2014 (same paper)
- UTD against FDTD (Schissler, Mückl & Calamia 2021, "Fast Diffraction Pathfinding for Dynamic Sound Propagation"):
  - On a 20-sided absorptive cylinder with diffraction order of at least 8, the two differ by less than 5 dB "for most positions in the scene". "UTD overestimates the level of the diffracted sound, particularly for areas with high-order diffraction."
  - On a fence-like scene, rendering all UTD diffraction paths "overestimated by about 20dB for the 84Hz frequency band". Rendering only the loudest path matched FDTD much better at low frequency.
  - Sources: [Schissler, Mückl & Calamia 2021, ACM TOG 40(4) 138](https://research.facebook.com/file/4521255027934795/Fast-Diffraction-Pathfinding-for-Dynamic-Sound-Propagation.pdf); [Consensus entry](https://consensus.app/papers/details/6d3ef5f57a1f53fbae263e17b5495bab/?utm_source=claude_desktop)
- Wavelength-dependent mesh simplification is needed to keep UTD sensible: "the UTD edge diffraction algorithm tends to be more accurate for longer edges; the presence of a high number of small edges can result in inaccurate results" — Schissler et al. 2014

**BTM against measurement**
- Remillieux et al. (2012) combined image sources (specular, rigid surfaces) with an extension of BTM. They validated it against a BEM solution and against field data from explosive sonic-boom simulants (0.78 kg detonation cord) around a rectangular building at Virginia Tech, finding "a good overall agreement". They then applied it to a cluster of six buildings with elevated source positions — [Remillieux et al. 2012, Applied Acoustics, via Consensus](https://consensus.app/papers/details/6b8196c9fc385b87b763b6f9e06e177c/?utm_source=claude_desktop). The abstract gives no dB error figure.

**Filter approximations to diffraction**
- Universal diffraction filter approximation (UDFA): half-order low-pass filters give the singly diffracted field from infinite wedges "with errors below 0.1 dB". Higher-order diffraction uses a heuristic extension — [Ewert et al. 2023, JASA](https://consensus.app/papers/details/608d7e1b23015cd19fc0ea8425b6dfd1/?utm_source=claude_desktop)
- The 2024 extension of UDFA to double edges: "For double diffraction at a three-sided barrier, errors are considerably decreased in comparison to a state-of-the-art sequential approach" — [Kirsch et al. 2024, Acta Acustica](https://consensus.app/papers/details/5c0a9bfd950e53bd9ff0a237614e22c5/?utm_source=claude_desktop)

**Which orders matter**
- From BTM-based room auralization: "second-order diffraction to nonshadowed receivers can often be neglected" — [Torres, Svensson & Kleiner 2001, JASA](https://consensus.app/papers/details/43078de916a65cf6b45e6775e670ef02/?utm_source=claude_desktop)
- Reflection order, comparing ISO 9613-2 and Pierce's solution behind buildings with Cadna/A ray/beam tracing:
  - Where vertical diffraction dominates, few reflections suffice.
  - Where horizontal reflections dominate, "the sound pressure level saturates after four reflections … at least four reflections are recommended".
  - ISO 9613-2 and Pierce's model can differ by "more than 10 dB" near a facade.
  - Source: [Wei, Cooke & Havie-Clark 2015, Proc. IOA 37(2)](https://www.apexacoustics.co.uk/wp-content/uploads/2016/02/Acoustics2015-How-many-reflections-should-I-use.pdf)
  - Caveat: the sources there are at ground level (industrial), not elevated.
- Wei, Van Renterghem & Botteldooren 2015 (JASA 138(5) 2744) is a simplified analytical model for shielded urban locations. It covers multi-edge diffraction over successive buildings and multiple reflections between parallel facades. It is built on Pierce's diffraction theory with the Fresnel integral approximated by trigonometric functions, and validated against Pierce's theory and a full-wave method for idealized configurations — [JASA abstract](https://pubs.aip.org/asa/jasa/article-abstract/138/5/2744/690076/Simplified-analytical-model-for-sound-level); [abstract text via SonicField](https://labs.sonicfield.org/library/simplified-analytical-model-for-sound-level-prediction-at-shielded-urban-locatio). The dB error figures were not accessible (paywall).

**Wei et al. 2014, "Urban background noise mapping: the general model"**
- Shielding is written as A_diffr = −10 log10(10^(−0.1 A_bar) + 10^(−0.1 A_can)). A_bar is the single path over the roof. A_can is the path with at least one reflection in the source and/or receiver canyon, which "quickly dominates A_bar".
- A_can,flat is an image-source analytic form fitted to 2D FDTD. Its standard deviation from the simulations is 2.7 dB.
- The roof-shape correction A_can,roof is −5 dB with both canyons present and −2.5 dB with one canyon.
- "Using the ISO standard to calculate A_bar underestimates the attenuation considerably" compared with in-situ measurements and FDTD.
- At an enclosed yard the official END noise map underpredicted the measured Lday by more than 14 dBA (position 9). The canyon term fixes low frequencies but not high frequencies, where a turbulence-scattering term "becomes therefore essential".
- Source: [Wei et al. 2014, Acta Acustica united with Acustica 100(6) 1098, author copy](https://publications.lib.chalmers.se/records/fulltext/206833/local_206833.pdf)

**Multi-edge approximations (Deygout, Epstein-Peterson)**
- In radio propagation, compared with the rigorous Vogler method: Bullington is the least accurate; Epstein-Peterson, Deygout and Giovaneli are "of similar accuracy". Deygout "has been reported to overestimate the diffraction loss". Epstein-Peterson error "increases when the obstruction objects are closely spaced" — [NTIA, A Comparative Analysis of Multiple Knife-Edge Diffraction Methods](https://its.ntia.gov/media/33472/dem_n-08.pdf); [Formulation of Multiple Diffraction by Trees and Buildings, PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882151/)

**Maekawa's detour law against an exact model**
- For double diffraction (r_S = r_R = d = 10λ), Stephenson's uncertainty-based diffraction (UBD) and a modified Maekawa detour law (MDL, with an equivalent thin screen after Pierce) were compared with Svensson's secondary-source model (SSM, the exact edge-diffraction toolbox). In the shadow zone both overestimate the diffracted energy, MDL more so. UBD is "by far better than the MDL".
- UBD against SSM: under 3 dB over a wide distance range, rising linearly to 5 dB at 100λ (in regions about −40 dB down).
- Source: [Stephenson, DAGA 2012, "Combining Ray Tracing with Higher Order Diffraction"](https://pub.dega-akustik.de/DAGA_2012/data/articles/000063.pdf)

### Inferences
- For an elevated source over a roof edge into a canyon, the dominant energy path is one roof-edge diffraction plus 0 to about 4 facade reflections in the receiver canyon (Wei 2014's A_can). That is a low-order image-edge problem. The high diffraction orders that game engines chase (5-10) matter for street-level sources wrapping around blocks, much less for a source above the roofline. For a low vertiport or UAM source near or below roof height, horizontal multi-edge paths around corners return, and the order requirement goes up.
- UTD's systematic errors (about 3 dB per diffraction against BTM, overestimated HF attenuation, energy-sum blow-up at high order and low frequency) point the same way. A Pierce/UTD-class single-edge coefficient is adequate for third-octave energy if the diffraction order is capped at 1-2 and only the dominant diffraction path per edge sequence is kept. The "loudest path" result in Schissler 2021 supports this.
- A full deterministic image-edge search per pair is unaffordable at NICE-OPS scale (Section 6). Any adoption has to amortize the path search across source states.

### Gaps
- Wei, Van Renterghem & Botteldooren 2015 (JASA) dB error statistics and the orders used: paywalled, abstract only.
- No acoustics-specific quantitative comparison of Deygout or Epstein-Peterson against BTM/BEM for building rows was found. The evidence above is from radio propagation.
- Hadden & Pierce's double-wedge and wide-barrier accuracy, and the classic Pierce (1974) wedge solution's error against BTM at low kr: no source fetched in this pass.
- Casalino et al. 2011 ("Simulation of Helicopter Community Noise in Complex Urban Geometry", AIAA J. — [AIAA](https://arc.aiaa.org/doi/10.2514/1.J050774)) carried no abstract in the index, so its method and validation are unknown here. Note that the 2011 AIAA helicopter-urban paper is by Casalino et al. (CIRA), not Remillieux. Remillieux et al. is the 2012 Applied Acoustics building and sonic-boom paper (with a 2011 JASA abstract). The request's "Remillieux 2011" probably conflates the two.

---

## 2. Ray, beam and Gaussian-beam tracing: diffuse reflections, diffraction inside ray tracers, convergence in canyons

### Takeaway
Stochastic ray and particle tracing handles diffuse facades naturally but converges slowly. In published comparisons 1000 rays gave errors of about 6.7 dB where 10,000 rays served as the reference.

Variance reduction cuts the ray count by about 10x for similar accuracy (2.3 dB average error at 1000 rays): radiosity-like caching, temporal reuse, bidirectional or A* search over precomputed edge-visibility graphs.

Gaussian beam tracing (GBT) removes caustic and ray-density artifacts and handles refraction. The published UAM/vertiport GBT models omit edge diffraction, which leaves building shadows 20-35 dB too quiet. Stephenson's uncertainty-relation diffraction is the main way to put diffraction into a particle tracer, at about 2-5 dB from exact.

### Cited Findings
**Path tracing against cached diffuse reflections (Schissler et al. 2014)**
- Diffuse reflections up to 10 orders. Path tracing with 1000 rays had an average error of 6.69 dB against 10,000-ray path tracing. The cached diffuse method at 1000 rays had 2.27 dB average error, at about 10.6 ms against 99.2 ms for 10,000-ray path tracing.
- Cache entries below −60 dB were discarded. Prior work (RESound) took about 250-500 ms for 3 orders of diffuse reflection with 200K rays on 60-280K-triangle models.
- Source: Schissler, Mehra & Manocha 2014, ACM TOG 33(4) 39, full text

**City-scale timings (Schissler et al. 2014, 4-core CPU)**

| Scene | Triangles | Sources | Specular / diffuse / diffraction order | Time per frame | Per source |
| --- | --- | --- | --- | --- | --- |
| Large city | 254,903 | 14 | 10 / 6 / 4 | 71.0 ms | about 5 ms |
| Small city | 89,792 | 8 | 10 / 6 / 6 | 19.3 ms | about 2.4 ms |

- Visibility-graph precomputation: 178.5 s and 23.8 MB for the large city.
- Source: Schissler et al. 2014, full text

**Multi-source city (Schissler & Manocha 2016)**
- City, 206,976 triangles, 50 sources clustered to 24, 10 specular and 50 diffuse orders: 47.2 ms propagation per frame, 27,711 paths. It uses listener-based backward tracing with "acoustic reciprocity for spherical sources" — [Schissler & Manocha 2016, ACM TOG, "Interactive Sound Propagation and Rendering for Large Multi-Source Scenes"](http://gamma.cs.unc.edu/MULTISOURCE/paper.pdf)

**Diffraction pathfinding (Schissler et al. 2021)**
- Silhouette-edge preprocessing, diffraction "flag" geometry, bidirectional path tracing and A* over the edge-visibility graph. It is 2.7-568x faster per source than the 2014 method for 10th-order diffraction, on a single thread.
- 100 rays with A* and cache took 3.6 ms, against 860 ms for 10,000 rays with A* and cache, in one scene.
- Source: [Schissler et al. 2021](https://research.facebook.com/file/4521255027934795/Fast-Diffraction-Pathfinding-for-Dynamic-Sound-Propagation.pdf)

**Uncertainty-relation (Heisenberg) diffraction**
- Stephenson's method diffracts sound particles with a probability that grows the closer they pass an edge. For double diffraction it is within 3 dB of Svensson's secondary-source model over a wide range, deviating up to 5 dB at 100λ in deep shadow, and is far better than Maekawa — [Stephenson DAGA 2012](https://pub.dega-akustik.de/DAGA_2012/data/articles/000063.pdf); [Stephenson 2010, Acta Acustica 96, via ResearchGate](https://www.researchgate.net/publication/233521923_An_Energetic_Approach_for_the_Simulation_of_Diffraction_within_Ray_Tracing_Based_on_the_Uncertainty_Relation)

**Particle tracing for diffuse facades**
- Can et al. (2015) ran 32,175 sound-particle-tracing simulations to quantify the error of neglecting facade diffusion in canyons (receiver at 1.5 m). The findings are in Section 3 — [Can et al. 2015, Applied Acoustics](https://consensus.app/papers/details/5b57a1256924590689a24ba3a6213ec9/?utm_source=claude_desktop)

**Gaussian beam tracing for UAM (Yunus et al. 2023, Aerospace Sci. Tech.)**
- The model:
  - Noise sphere input, 3D terrain and buildings, wind.
  - Reflections only: "the present GBT method does not account for the diffracted wavefield behind the building blocks". Building dimensions were chosen 1.75x the first-BPF wavelength "to minimize the diffraction effects".
- What buildings did to the footprint:
  - +0-5 dB in illuminated zones from multiple reflections.
  - Shadow differences over 20 dB, which the authors acknowledge would be smaller with diffraction.
  - Against a straight-ray propagator, GBT "reduces the prediction error by 5 dB in the illuminated zones and 35 dB in the terrain shadow zones".
- Cost: 55.7 s per GBT case against 3.7 s for straight rays, on a 36-core Xeon Gold 6140.
- Source: [Yunus et al. 2023, author copy](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf); [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1270963823003073)

**Other Gaussian-beam work**
- MANIA-H is a GBT framework with moving sources, multiple urban reflections and refraction. A "Virtual Beam Extending" fix handles beam truncation at reflecting surfaces. It was validated against conventional-helicopter flight-test SPL time histories and demonstrated over downtown New York City "with practical computational cost". No numbers are in the abstract — [Yang et al. 2025, AIAA Aviation](https://consensus.app/papers/details/58b4673d283651ee90e8c5a204973d59/?utm_source=claude_desktop)
- A CUDA GBT with dynamic parallelism is "up to 790 times faster in city environment and 188 times faster in open plane environment" than the original CPU program. The Gaussian beam summation (gathering all beams at each receiver) dominates as the number of rays grows — [Sheng et al. 2025, arXiv:2501.13382](https://arxiv.org/html/2501.13382v1)

**Ray tracing with radiosity, and stochastic ray tracing, for sonic booms**
- Riegel et al. 2024: urban canyons modeled without diffraction show "strong shadow zones where no sound was detected as this model does not include diffraction". Diffuse reflection had the largest effect on the received signatures — [Riegel et al. 2024, JASA 155(2) 1527](https://consensus.app/papers/details/41a56437547f599480b280dc8cec544c/?utm_source=claude_desktop)
- Stochastic ray tracing was investigated as a lower-memory replacement for radiosity — [Marini et al. 2024, JASA abstract](https://consensus.app/papers/details/d31263252d0550f3b441cf2698dd8337/?utm_source=claude_desktop)

**Eikonal / jet-marching alternative**
- Eikonal and transport equations are solved on a tetrahedral mesh, with UTD for edges, as a non-ray-shooting geometric-optics alternative for building geometries — [Potter & Cameron 2022, arXiv](https://consensus.app/papers/details/572e40390b5e5e1da6075c48828762a0/?utm_source=claude_desktop)

### Inferences
- The two prominent GBT UAM models omit edge diffraction (Yunus explicitly). Every published geometric model that omits diffraction (Yunus GBT, Riegel radiosity) reports physically wrong, near-silent shadows behind buildings. For an elevated source the roof-edge diffracted path is first-order essential. A tracer without diffraction is disqualified for NICE-OPS's shielded receivers.
- Stochastic tracing's Monte Carlo noise (several dB at 1000 rays) is incompatible with NICE-OPS's bit-reproducible, golden-tested design and its peak-search bounds. A deterministic path enumeration (image-edge), or a deterministic beam method, fits the codebase's culture better than per-pair stochastic sampling. Stochastic methods may still be acceptable inside an offline per-receiver precomputation that is then frozen.

### Gaps
- No published convergence study of ray count against canyon level for an *elevated* source was found. All the convergence numbers are for games and rooms, or for street-level sources.
- Phonon tracing (Bertram et al.) for outdoor scenes: no validation found in this pass.
- Yunus et al.'s GBT accuracy against FEM is said to come from their earlier papers, which were not fetched. The number of receivers behind the 55.7 s was not extracted.

---

## 3. Specular or diffuse facades: how much it changes canyon levels, and realistic scattering coefficients

### Takeaway
Facade diffusion is a first-order effect in canyons with many reflections, and a small one for single-reflection or isolated-building cases.
- **Canyons, ray models:** it shifts levels by about −2 to +4 dB, and up to +10 dB once street fittings are included.
- **Shielded canyons, 2D full-wave:** it adds about 10 dB of shielding above about 500 Hz.
- **Isolated building:** it is negligible (Riegel 2019).

Measured or estimated facade scattering coefficients are small to moderate. Onaga et al. put absorption plus scattering at 0.1-0.25 for real streets. Ismail & Oldham found scattering "small and not very sensitive" to irregularity. A search summary gives 0.07-0.74 across five measured streets. Diffuse energy grows with reflection order, so specular is a fair model for the first 1-2 reflections and diffuse dominates after that.

### Cited Findings
**Kang 2000**
- Radiosity (diffuse boundaries) against image sources (geometric boundaries). With geometric boundaries, "the sound attenuation along the length becomes significantly less; the RT30 is considerably longer; and the extra attenuation caused by air or vegetation absorption is reduced". There are "considerable differences between the sound fields resulting from the two kinds of boundaries" — [Kang 2000, JASA 107, 1394](https://consensus.app/papers/details/3d8adc693cac522ca61dd2bc02538b6b/?utm_source=claude_desktop). Full-text dB values not accessed.

**Ismail & Oldham 2005**
- Scale-model facades with realistic irregularity: "the measured value of the scattering coefficient is found to be small and not very sensitive to the degree of surface irregularity". "A progression of energy from a specular reflection field to a diffuse reflection field for successive orders of reflections is observed". For higher orders, "the dominant mechanism … is via random scattering" — [Ismail & Oldham 2005, Applied Acoustics](https://consensus.app/papers/details/0887509e867b5b7bb49df13cd4cb6fb5/?utm_source=claude_desktop)

**Onaga & Rindel 2007**
- Image method plus radiosity. Scattering raises SPL at short distances and lowers it at long distances. The estimated "sum of absorption coefficient and scattering coefficient of facades of actual urban streets range from 0.1 to 0.25" (from measured reverberation times) — [Onaga & Rindel 2007, Applied Acoustics](https://consensus.app/papers/details/7b1aedaadcc350d2965ee9abba3a0afb/?utm_source=claude_desktop)

**Yu et al. 2019 (IMCRI)**
- An iterative image-source and radiosity hybrid ("image patches") validated against full-scale measurements in five streets with different scattering. Agreement in SPL and RT was good, and better than existing energy models for RT — [Yu et al. 2019, Applied Sciences 9(14) 2803](https://consensus.app/papers/details/47b4e8dbba0a5b53873b28eefedacf38/?utm_source=claude_desktop)
- The search-engine summary of related work states that "the estimated scattering coefficient in five measured streets varies from 0.07 to 0.74" — [MDPI page, full text 403](https://www.mdpi.com/2076-3417/9/14/2803). Taken from a snippet, not verified in full text.

**Can et al. 2015**
- 32,175 particle-tracing simulations. "Diffuse reflections can lead to reduction of 2 dB to an increase of 4 dB of sound pressure levels in the absence of fittings, and can lead to an increase of 10 dB … in the presence of fittings". Regressions are provided to correct specular ray codes — [Can et al. 2015, Applied Acoustics](https://consensus.app/papers/details/5b57a1256924590689a24ba3a6213ec9/?utm_source=claude_desktop)

**Van Renterghem et al. 2006 (2D coupled FDTD-PE, canyon to canyon)**
- "Diffusely reflecting façades and balconies lead to an important increase in shielding compared to flat façades." "Starting from about 500 Hz, a gain in shielding of about 10 dB" for a profiled facade (scattering simulated in 2D only).
- "Rigid façades yield significantly lower shielding compared to partly reflecting façades".
- Source: [Van Renterghem, Salomons & Botteldooren 2006, Applied Acoustics 67, 487](https://users.ugent.be/~tvrenter/publicaties/streetcanyonsFDTDPE.pdf)

**Riegel et al. (sonic booms)**
- 2019, isolated building: including diffuse reflections "does not have a dramatic effect … and does not improve the comparisons with the measurements … not worth the additional computational cost for an isolated building" (validated against NASA SonicBOBS) — [Riegel et al. 2019, JASA 145(4)](https://consensus.app/papers/details/e8ffde29ba555a1b9de471bd8a38b223/?utm_source=claude_desktop)
- 2024, canyons: "Diffuse reflections was the parameter that had the greatest impact". More reflections led to "a dramatic reduction in the PLdB" — [Riegel et al. 2024](https://consensus.app/papers/details/41a56437547f599480b280dc8cec544c/?utm_source=claude_desktop)

**Angle-dependent scattering**
- The ISO 17497-1 scattering coefficient "is the most uncertain standard metric applied in geometrical acoustics". RWTH is testing angle-dependent scattering (sound steering, retro-reflection) for building facades — [Heimes et al. 2023, JASA abstract](https://consensus.app/papers/details/2af3cf63851057fd92b528cf73e6b49f/?utm_source=claude_desktop)
- In a street canyon, bidirectional (BRDF-like) scattering gives "significant differences in the energy time curve" compared with random-incidence scattering. For a single facade, random-incidence radiosity agrees with ray tracing — [Fatela et al. 2026, JASA (CC BY)](https://consensus.app/papers/details/332fc91cc5a25ec98574bf9d780f2ec0/?utm_source=claude_desktop)

**Measured canyon impulse responses**
- About 100 inner-city street canyons in Ghent were measured (Thomas et al. 2013). They are used to build frequency-dependent stochastic envelope tails on top of an image-source model — [Vos et al. 2024, Inter-Noise](https://consensus.app/papers/details/2de2cf9110d755518c57fb4a10182a27/?utm_source=claude_desktop)

**Facade structure as periodic scattering (radio analogue)**
- Facade scattering from nearly periodic windows and masonry "takes the form of diffraction orders, rather than diffuse scattering from a randomly rough surface" — [Pongsilamanee & Bertoni 2004, IEEE TAP](https://consensus.app/papers/details/b592c4a263b75e79a4dae528899f7fe3/?utm_source=claude_desktop)

**Empa sonAIR (aircraft noise mapping)**
- The REFLECT module models facade scattering analytically (Heutschi's model). It considers energy exchange up to third order (source-reflector, reflector-reflector, reflector-receiver), with visibility checks only and no diffraction.
- Results near an airfield:
  - Building shielding cut the rating level by more than 2 dB over large areas, and by up to 12 dB on facades facing away from the runway.
  - Reflections raised it only locally, by up to 3-4 dB on facades facing the runway.
- Source: [Ramseier et al. 2023, Forum Acusticum, DOI 10.61782/fa.2023.0587](https://dael.euracoustics.org/landing_pages/aaua/51980.html)
- Heutschi (2009, Acta Acustica 95) holds that facade reflections are "more likely specular or coherent at low frequencies and increasingly diffuse or incoherent at higher frequencies". The degree of incoherence was estimated for typical facades numerically — [Heutschi, IOA proceedings, as summarized by search](https://www.ioa.org.uk/system/files/proceedings/k_heutschi_coherent_and_incoherent_reflections_at_building_facades.pdf). The PDF link returned 404 when fetched, so this claim is from the search summary.

### Inferences
- **What NICE-OPS should implement:** specular image paths (with a Pierce/UTD roof-edge coefficient) for orders 0-2, plus a diffuse energy term for higher orders. Wei 2014's A_can plus turbulence scattering, Can 2015's regressions, and IMCRI-style radiosity are three versions of that hybrid. A fully specular canyon model will overpredict levels along long streets (Kang) and underpredict shielding at high frequency (Van Renterghem 2006).
- **Default facade parameters:** scattering s ≈ 0.1-0.3, frequency dependent and rising with frequency, with absorption α ≈ 0.05-0.1 for masonry or glass. These defaults are inferred from the cited ranges. Sensitivity to s is a few dB in canyons, so the uncertainty in s is comparable to the method error of the propagation model itself.
- **Isolated buildings and sparse suburbs** (the Empa and Riegel 2019 regime): specular reflection plus shielding captures most of the effect.

### Gaps
- Kang 2000's numerical dB differences (full text not accessed).
- The 0.07-0.74 street scattering range is from a search snippet and should be confirmed from Yu et al. 2019 or its cited source.
- No measured, frequency-resolved facade scattering coefficients (for example ISO 17497-1 on real facades) were found in this pass.

---

## 4. Coherent or incoherent summation: when phase matters (rotor harmonics, narrow canyons)

### Takeaway
For third-octave energy at mid and high frequency, incoherent (energy) summation is standard and appropriate. Facade irregularity, turbulence and band averaging destroy the interference pattern, and coherent models overpredict levels at shielded positions.

Phase matters at low frequency in narrow canyons. There the field is modal: an image-plus-modal model agrees with experiment only up to about 100 Hz at full scale, and above that, leaky-mode behavior is needed. For a coherent multi-path model, small path-delay errors (about 1 ms) visibly shift comb patterns.

Rotor harmonics at 10-30 Hz have wavelengths (11-34 m) comparable to canyon widths. That is outside the geometric regime altogether.

### Cited Findings
- Coherent ray summation in narrow canyons (under 10 m wide) is needed to capture the "interference effects due to multiple reflections from building facades and ground surfaces", and was validated in an anechoic chamber and outdoors — [Iu & Li 2002, JASA](https://consensus.app/papers/details/700595b2db9656eb97e3cd41a2c57e19/?utm_source=claude_desktop)
- Incoherent image sources replaced by rows of line sources give an exact analytic sound-energy solution, validated indoors and outdoors — [Li et al. 2009, "A note on noise propagation in street canyons", JASA](https://consensus.app/papers/details/9b88f9f425595ff8b702f2f16883471b/?utm_source=claude_desktop)
- "The use of a coherent line source for shielded urban canyon observer positions leads mostly to an overprediction of levels and can yield erroneous results for noise abatement schemes." The 2.5D equivalent-sources solution needs a "very fine frequency resolution … due to the multiple reflections in the canyons" — [Hornikx & Forssén 2007, JASA](https://consensus.app/papers/details/856d51f760d4589bb057e7e8eaa20024/?utm_source=claude_desktop)
- Low-frequency street canyon (scale model, 2D FDTD, modal decomposition plus 2D image sources): "quite good agreement until f = 100 Hz at full scale. For higher frequency … the notion of leaky modes must be introduced" — [Richoux, Ayrault, Pelat & Félix 2009, Applied Acoustics](https://consensus.app/papers/details/b81b9d038b6d51d3890b30087f72b956/?utm_source=claude_desktop)
- For the open top of a street, energetic (ray) methods assume it is perfectly absorbing. At low frequency the reflection off the geometric discontinuity at the opening "is still an open question" — [Pelat et al., Noise in the Built Environment 2010](https://consensus.app/papers/details/a04d68c7d8c95f328357cc9dc4144412/?utm_source=claude_desktop)
- Sonic boom over urban areas (2D Euler): canyon resonances produce low-frequency oscillations whose amplitude, frequency and decay vary from canyon to canyon. Despite this, "the median of the perceived noise levels mostly differs by less than 1 dB from the value obtained for flat ground" — [Dragna et al. 2022, JASA](https://consensus.app/papers/details/f4929c18a0735403ba75590ddb7cdc76/?utm_source=claude_desktop); [Dragna et al. 2022, isolated and multiple buildings](https://consensus.app/papers/details/41cdca47c0b951c78544c534fa336c1f/?utm_source=claude_desktop)
- Coupling a curved atmospheric path to an urban image-edge model (Schäfer et al. 2022):
  - Path-delay errors up to 1.2 ms (about 0.01% of the delay) change the interference pattern of many-path transfer functions. The maximum deviation was "approx. 0.5 dB at 4 kHz", larger at higher frequencies (for example 1.5 dB at 8 kHz, against about 75 dB of air attenuation).
  - Mean psychoacoustic deviation was 0.35 sone in loudness and 0.001 acum in sharpness.
  - Source: [Schäfer, Reich & Vorländer 2022, Acta Acustica 6, 28](https://d-nb.info/1262974569/34)
- Wei et al. 2014 subtract 3 dB from facade-mounted measurements "to account for the coherence between the direct and reflected waves near the façades for very low frequencies (< 250 Hz)" — [Wei et al. 2014](https://publications.lib.chalmers.se/records/fulltext/206833/local_206833.pdf)

### Inferences
- **NICE-OPS should stay incoherent (energy) for the urban paths in third-octave bands from about 100-200 Hz up.** It already works in energy, and the literature (Hornikx & Forssén, Heutschi, Wei) says coherent summation over-predicts at shielded urban positions once facade roughness and turbulence are present.
- **At the facade itself**, a microphone or receiver on the facade sees a coherent doubling below about 250 Hz. A +3 dB pressure-doubling rule is simpler and safer than tracking phase.
- **Below about 100 Hz**, which covers the rotor fundamentals and first harmonics at 10-30 Hz, a geometric urban model is not physically justified. Canyon widths are 0.3-3 wavelengths, the field is modal and leaky, and diffraction fills shadows almost completely.
  - Option one: apply no urban correction below a crossover, treating buildings as transparent to energy at 10-50 Hz.
  - Option two: use a wave-based reduced model (2.5D equivalent sources or modal) for the lowest bands.
  - This crossover is an inference from Richoux's 100 Hz agreement limit and wavelength arguments, not a published rule.

### Gaps
- No published comparison of coherent against incoherent canyon summation for an *elevated, moving* tonal source (BPF harmonics) was found.
- Heutschi's quantitative incoherence-against-frequency curves for typical facades were not accessible.

---

## 5. Coupling with atmospheric refraction: virtual-source coupling and other hybrids

### Takeaway
The established geometric coupling is Schäfer, Reich & Vorländer's virtual source method (VSM).
1. Find the curved free-field eigenray from source to receiver (the Atmospheric Ray Tracing framework, ART).
2. Replace the source by a virtual point source along the arrival direction, at distance c_urban × travel time.
3. Run a straight-ray image-edge model from that virtual source.
4. Re-normalize spreading, air absorption and launch direction to the curved path.

Errors were at most 2.5° in arrival elevation, under about 1 dB in spreading, and 1.2 ms in delay, worst when the source is low and far. This is the closest published analogue to bolting an urban model onto NICE-OPS's existing stratified eigenrays.

Gaussian beam tracers (Yunus 2023, MANIA-H 2025) couple refraction and urban reflection natively, but without diffraction.

### Cited Findings
**How VSM works**
- The ART framework gives curved eigenrays in a stratified, moving atmosphere. The image edge model assumes a homogeneous, still medium and straight paths.
- The virtual source lies on the line cast from the receiver opposite the atmospheric path's incident direction, at d_atmos = c_urban · τ_atmos. c_urban is taken "at receiver altitude".
- Delay was chosen over spreading as the defining parameter "as it is used to generate Doppler shifts". Both choices were tested and "lead to similar results".
- Source: [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**Corrections applied to the urban paths**
- Spreading: g'_spread = (g_spread,atmos / (1/d_atmos)) · 1/l_path. For diffracted paths the UTD distance factor sqrt(s'/(s(s'+s))) is used instead of 1/d.
- Air absorption: A_air,path(f) = α_air(f) · (l_path − d_atmos) + A_air,atmos(f), with α evaluated at receiver altitude.
- Launch direction: the virtual source is rotated by the elevation and azimuth offsets between the atmospheric and straight launch directions, to preserve directivity lookup.
- Source: [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**Proof of concept**
- Setup:
  - ISA temperature profile with a logarithmic wind profile, aircraft flying upwind on a take-off.
  - Street canyon between two 10×10×255 m buildings.
  - Image-edge orders: combined order up to 3, reflection order up to 2, diffraction order up to 2, with back-face culling and occlusion detection.
- Errors against new curved eigenrays traced to each first interaction point:
  - Arrival elevation up to 2.5°.
  - Spreading near 0 dB for most positions, approaching the 1 dB JND early in the trajectory.
  - Delay up to 1.2 ms.
  - Error grows with distance and is largest when the launch elevation is closest to horizontal.
- The virtual source sat up to 360 m and 460 m from the aircraft at the start and end of the flyover.
- Source: [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**Stated assumptions and limitations**
- The assumption is that "atmospheric propagation effects do not change significantly among the different urban paths", since only one free-field eigenray is used.
- Buildings would in reality modify the wind near the ground, so "the assumption of a stratified atmosphere would not hold anymore, at least in the vicinity of the ground".
- Rendering artifacts appear when dominant urban paths switch on or off, "likely caused by only considering the phase of the diffraction filters indirectly".
- Source: [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**Open-source components**
- ARTMatlab and pigeon (image edge model) in the RWTH GitLab ITAGeometricalAcoustics repository. The paper is CC BY 4.0 — [Schäfer et al. 2022](https://d-nb.info/1262974569/34); [Acta Acustica landing page](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2022/01/aacus220008/aacus220008.html)
- The image edge model (Erraji, Stienen & Vorländer 2021) is a deterministic, open-source path algorithm producing reflection and diffraction paths with metadata (angles, materials, wedge geometry) for later acoustic modeling — [Erraji et al. 2021, Acta Acustica, DOAJ](https://doaj.org/article/31fed7f4cc1545faa5439d7dfb545ac9)

**Palenda & Vorländer 2026**
- An open-source auralization tool for "arbitrary urban environments" built on geometrical acoustics and UTD. It covers "mixed reflections and diffraction effects, including higher-order diffraction and Doppler effects", in shadow and non-shadow regions.
- About 28 benchmark scenarios: wedges, blocks, building corners, curved surfaces, canyon-like structures, static, linear, circular and 3D trajectories.
- Its stated limitation is the infinite-edge assumption of UTD. It found that "diffraction and higher-order diffraction effects" significantly influence accuracy.
- Supplementary code and data are Windows-only and licensed CC BY-NC-SA 4.0.
- Sources: [Zenodo record 19112715](https://zenodo.org/records/19112715); [Acta Acustica 2026](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2026/01/aacus250149/aacus250149.html)

**Gaussian beam hybrids**
- GBT propagates a noise-sphere source through an "inhomogeneous anisotropic atmosphere" with 3D buildings. Wind changed the footprint by up to 23 dB and raised vertiport-pad levels by up to 30 dB in one case, but there was no diffraction — [Yunus et al. 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)
- MANIA-H combines GBT with refraction by sound-speed and wind variation and multiple urban reflections — [Yang et al. 2025](https://consensus.app/papers/details/58b4673d283651ee90e8c5a204973d59/?utm_source=claude_desktop)

**Stratified ray tracing for urban sonic booms**
- Riegel et al. 2021 exploit the z-stratification of the atmosphere to also slice the building collision tests ("only the sections of the building that are within the striation are checked"), which cut computation time — [Riegel et al. 2021, JASA abstract](https://consensus.app/papers/details/ad2fe445ec41558180ed0ea0e62d7f96/?utm_source=claude_desktop)

**Full-wave reference for meteorology in canyons**
- In 2D FDTD-PE canyon-to-canyon results, downwind propagation decreases shielding by "near 10 dB" (u* = 1 m/s), with refraction the main mechanism. Upwind, turbulent scattering dominates. "The combination of effects … is in general not a simple addition of the separate effects" — [Van Renterghem et al. 2006](https://users.ugent.be/~tvrenter/publicaties/streetcanyonsFDTDPE.pdf)

### Inferences
- **VSM maps directly onto NICE-OPS.** NICE-OPS already computes, per pair, the arrival direction, travel time, tube spreading and per-band absorption of each stratified eigenray (`RayPath`, `Arrival` in `src/propagation.hpp` and `src/refraction.hpp` — [CLAUDE.md](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md)). An urban term could be an excess attenuation in dB per band, a function of the arrival direction at the receiver and of local geometry. The atmospheric eigenray then supplies everything else.
  - This is cheaper than VSM proper, because it skips the virtual-source position and the per-path absorption re-normalization. The extra absorption of urban detours of a few tens of meters is negligible below about 4 kHz.
  - The error budget Schäfer measured (2.5° in arrival elevation, under 1 dB) bounds what a direction-only interface costs.
- **Where the interface breaks:** VSM's "one atmospheric eigenray per receiver" assumption fails exactly where NICE-OPS's grazing and shadow-zone modeling is most delicate, at low elevation (Schäfer: largest error near-horizontal). In a refractive shadow, the arrival is a creeping or diffracted wave with no well-defined direction. The urban term there has to fall back to something direction-agnostic, for example a diffuse-canyon correction.
- **Multiple atmospheric arrivals:** NICE-OPS keeps several (ducted and bounced). Each would carry its own urban transfer by its own arrival direction, summed in energy, which is consistent with NICE-OPS's existing incoherent arrival sum.

### Gaps
- Palenda & Vorländer 2026: the quantitative benchmark errors, the orders used and the runtimes were not readable (publisher 403, and ResearchGate 403).
- The license of the RWTH ITAGeometricalAcoustics, pigeon and ART code was not verified: the GitLab page is behind an Anubis bot wall. Check the repository's LICENSE file directly before any reuse.
- No validation of VSM-type coupling against a full-wave (PE or FDTD) reference with both refraction and buildings was found. Schäfer validated only against re-traced curved rays, not against wave solutions or measurements.

---

## 6. Runtime per source-receiver pair and per receiver; GPU, Embree and OptiX implementations

### Takeaway
Published geometric urban methods cost between about 0.1 ms and 1 s per source-receiver evaluation:
- **Image-edge (Schäfer, Matlab and C++):** 0.8-1.2 s per pair.
- **Game-engine UTD path search with diffuse cache:** 2-5 ms per source per frame on 4 cores, and down to about 0.2-3.6 ms single-thread with A* graph search.
- **GBT:** tens of seconds per case on 36 cores. GPU GBT is up to 790x faster.

NICE-OPS's existing pair evaluation is 19 ns to 1 µs. A per-pair urban path search at 10^8 to 10^9 pairs is therefore 10^3 to 10^6 times over budget: roughly 10^5 to 10^9 core-seconds per run. Only formulations that do the geometric search once per receiver (or per receiver and source direction bin) and reuse it across thousands of source states are viable.

### Cited Findings
**Image-edge and VSM (Schäfer et al. 2022)**
- 0.8-1.2 s per single source/receiver position for ART + VSM + image edge model on a workstation, combined order 3.
- 83 s of audio (3,900 blocks) took about 1 h of precomputation.
- "A lot of computation time is spent transferring data" between Matlab and C++.
- Source: [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**Game-engine propagation (UNC and Meta)**

| Method and scene | Orders | Time | Hardware |
| --- | --- | --- | --- |
| Schissler 2014, large city, 254,903 triangles, 14 sources | spec 10, diffuse 6, diffraction 4 | 71.0 ms per frame (about 5 ms per source) | 4-core CPU |
| Schissler 2014, small city, 8 sources | | 19.3 ms per frame | 4-core CPU |
| Schissler & Manocha 2016, city, 206,976 triangles, 50 sources in 24 clusters | 10 spec, 50 diffuse | 47.2 ms per frame | |
| Schissler et al. 2021, diffraction search | 10th order | 2.7-568x faster per source than 2014; 0.2-3.6 ms at 100 rays | single thread |

- Schissler 2014 also reports 178.5 s of visibility-graph precomputation for the large city.
- Schissler 2021 preprocessing times by scene: Sun Temple (4.7M triangles) 52.9 s simplification and 16,175 s graph; Warehouse 11.3 s and 83.6 s; Sibenik 0.76 s and 8.07 s.
- Sources: Schissler et al. 2014 (full text); [Schissler & Manocha 2016](http://gamma.cs.unc.edu/MULTISOURCE/paper.pdf); [Schissler et al. 2021](https://research.facebook.com/file/4521255027934795/Fast-Diffraction-Pathfinding-for-Dynamic-Sound-Propagation.pdf)

**Gaussian beams**
- 55.7 s per GBT case against 3.7 s for straight rays (Xeon Gold 6140, 36 cores). The GBT cost "remains constant" with wind and with frequency (the same sphere resolution) — [Yunus et al. 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf)
- CUDA GBT is up to 790x faster than CPU in a city and 188x in free field. Beam summation at receivers dominates cost as the ray count grows — [Sheng et al. 2025](https://arxiv.org/html/2501.13382v1)

**Ray-throughput ceilings**
- An OptiX-based ray tracer traced "0.57 Giga rays … in 1 second", 338x a kd-tree CPU tracer. RTX hardware is quoted at up to about 6-10 Grays/s. Note: this is from an EM shooting-and-bouncing-rays paper, via search summary — [PLOS ONE 2021, OptiX SBR](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0253743)
- Embree is an "open source ray tracing framework for x86 CPUs", fully open "to avoid licensing restrictions" (Wald et al. 2014, ACM TOG 33(4), "Embree: A Kernel Framework for Efficient CPU Ray Tracing"; read from a local copy). The license is Apache 2.0 per general knowledge, not confirmed in the fetched text.
- GPU guided multi-view ray tracing computes "early specular paths and first order diffraction" interactively — [Taylor et al., Guided Multiview Ray Tracing for Fast Auralization](http://gamma.cs.unc.edu/Sound/Guided/Guided.pdf)

**Engineering urban mapping**
- sonAIR (Empa): a single flight's propagation takes "a few minutes" with BASIC (direct sound, shielding by terrain and buildings). Adding building reflections (REFLECT, third-order energy exchanges, visibility only) adds about 10-20 minutes. Adding forest reflections brings it to about 2 hours (i7 workstation, 32 GB) — [Ramseier et al. 2023](https://dael.euracoustics.org/landing_pages/aaua/51980.html)
- NoiseModelling (CNOSSOS-EU) uses a "rubber-band like" path finder with an R-tree. It ignores rays going backwards to source or receiver: "if a receiver is located inside a U-shaped building, only diffractions on horizontal edges will be taken into account" — [NoiseModelling docs](https://noisemodelling.readthedocs.io/en/latest/Numerical_Model.html)

**NICE-OPS reference cost**
- One evaluation is 19.4 ns (broadband). Spectral is about 150 ns, `--ground_effect` about 1 µs. The default run is 94.3M evaluations — [CLAUDE.md](/Users/evg5332/Code/D/NICEOPS/CLAUDE.md)

### Inferences
**Order-of-magnitude budget for 10^9 pairs**

| Method | Per-pair cost | Total | Wall time |
| --- | --- | --- | --- |
| Schäfer-style | 1 s | 10^9 core-s | infeasible |
| UNC/Meta game engines | about 3 ms | 3×10^6 core-s | about 1 day on 32 cores |
| Ray launch only, tracing 10^3 rays per pair at 1 Grays/s | about 1 µs | 10^3 s | |

The last row shows that raw intersection throughput is not the bottleneck. Path validation, image-source bookkeeping and per-band diffraction coefficients are.

**Exploiting NICE-OPS's structure**
- For a source 50-1500 m away, the urban transfer at a receiver depends on the incidence direction (and frequency band) much more than on the exact source position. The far-field and plane-wave assumption is the same one Schäfer invoked: aircraft are "very far away from the urban environment if compared to the dimensions of the buildings".
- Precompute per receiver a table U_r(band, θ, φ) of urban excess attenuation over arrival directions, with a few hundred to a few thousand direction bins.
- That takes 10^5 receivers × 10^3 directions = 10^8 geometric solves once per city, reusable across runs and trajectories, and each solve can be an expensive BTM/UTD image-edge solve.
- The hot loop then adds one table lookup, about 10-50 ns, in the same style as NICE-OPS's existing `GridAxis` and propagation-table lookups.
- This is an inference, and so is its limit: it fails for low vertiport and UAM sources within about 10 building heights of the receiver, where curvature of the incident wavefront and source position within the canyon matter. Those cases need a per-pair (or per-state-cluster) solve.

**Reciprocity**
- Receiver-side tracing amortizes over all source states. Schissler & Manocha 2016 use backward (listener-based) tracing for exactly this reason.

### Gaps
- No published per-receiver direction-table ("urban HRTF") implementation for aircraft noise was found. Its accuracy against per-pair image-edge solutions is untested.
- No acoustic-specific Embree or OptiX outdoor benchmark with rays/s *and* validated levels was found. Throughput numbers come from graphics and EM.
- Palenda 2026 runtimes are unknown (paywall).

---

## 7. What auralization and virtual-acoustics groups (RWTH ITA VA, Pachyderm, RAVEN, GSound) use outdoors, what their validations showed, and the code licenses

### Takeaway
RWTH ITA (VA, ART, pigeon/image edge model, the Palenda 2026 tool) is the only group building specifically for aircraft over urban areas. Its stack is deterministic image-edge with UTD, coupled to curved-ray atmospheres via the virtual source method. Validation so far is against re-traced rays and simple benchmarks, not full-wave or field data in cities.

The UNC/Meta line (GSound, Schissler) is UTD plus diffuse path caches, validated against BTM (0.7-3.6 dB) and FDTD (under 5 dB mostly, but 20 dB errors at 84 Hz with all paths). RAVEN (RWTH, closed) and Pachyderm (GPL) are room-acoustics-first.

License picture for an MIT C++ codebase:
- **BSD:** EDtoolbox. Linkable or portable.
- **GPL v3:** old EDB toolbox, Pachyderm, NoiseModelling. Benchmark-only.
- **CC BY-NC-SA 4.0:** Palenda 2026 supplementary. Non-commercial share-alike, so benchmark-only.
- **All rights reserved:** GSound.
- **Not verified:** RWTH ITAGeometricalAcoustics.

### Cited Findings
**RWTH Virtual Acoustics (VA)**
- "A real-time auralization framework for scientific research … open source and fully controllable", initiated around 2012 by Frank Wefers at ITA/IHTA — [Zenodo VA record](https://zenodo.org/records/13788752)
- VA rendered the Schäfer 2022 flyover. The paper calls urban reflections and diffraction "of great importance": the perceived aircraft position flipped sides because the direct path was blocked and the reflection off the opposite building dominated — [Schäfer et al. 2022](https://d-nb.info/1262974569/34)

**RWTH image-edge model (Erraji et al. 2021) and the urban tool (Palenda & Vorländer 2026)**
- Deterministic and open source, with UTD diffraction.
- Palenda 2026 is open source with supplementary material under **CC BY-NC-SA 4.0**. It reports the infinite-edge UTD assumption as a limitation.
- Sources: [Erraji et al. 2021](https://doaj.org/article/31fed7f4cc1545faa5439d7dfb545ac9); [Zenodo 19112715](https://zenodo.org/records/19112715)
- An urban sound auralization and visualization framework was case-studied at IHTApark (RWTH) — [Sustainability 2022, 14(4) 2026](https://doi.org/10.3390/su14042026). Not read in detail.

**Edge-diffraction toolboxes (Svensson, NTNU)**
- Current EDtoolbox: "external, convex Neumann scattering problems" in the time or frequency domain. Higher-order diffraction in the frequency-domain version, lower-order only in the time-domain version. **BSD License** — [EDtoolbox GitHub](https://github.com/upsvensson/Edge-diffraction-Matlab-toolbox)
- The older EDB toolbox (EDB1, last released 2012) is under the **GNU GPL v3 or later** — [EDB toolbox files](https://ulfps.folk.ntnu.no/software/edbtoolbox/edbtoolboxversionhistory.txt); [EDB1 source mirror](https://code.soundsoftware.ac.uk/hg/human-echolocation/file/tip/private/EDB1wedge2nd.m). The license statement is from the search summary of these pages.

**Other packages**
- **Pachyderm** (Rhino/Grasshopper plug-in, ORASE nonprofit): open source, **GPL** — [food4rhino](https://www.food4rhino.com/en/app/pachyderm-acoustical-simulation); [Ladybug forum](https://discourse.ladybug.tools/t/acoustic-simulations/442)
- **GSound** (Schissler, UNC): "all rights reserved" by Carl Schissler/UNC. The pygsound Python wrapper and the derivative GSound-SIR are described as Apache 2.0 — [pygsound on PyPI](https://pypi.org/project/pygsound/0.3); [arXiv 2503.17866](https://www.arxiv.org/abs/2503.17866?context=cs). From the search summary. Verify which parts of the wrapped C++ core the permissive license actually covers.
- **NoiseModelling** (CNOSSOS-EU, UMRAE/Lab-STICC): **GPL v3** — [NoiseModelling license](https://noisemodelling.readthedocs.io/en/latest/License.html)

**Validations by family**

| Family | Reference | Result |
| --- | --- | --- |
| UTD (GSound line) | BTM | 0.67 / 3.10 / 3.61 dB |
| UTD | FDTD | under 5 dB mostly; +20 dB at 84 Hz when all paths are summed |
| Image source + BTM (Remillieux) | BEM and field explosive tests | "good overall agreement" |
| Ray tracing + radiosity (Riegel) | NASA SonicBOBS | "general agreement" |
| GBT (Yunus) | FEM, in earlier papers | |
| GBT (MANIA-H) | helicopter flight test | |
| Engineering canyon model (Wei 2014) | 9 Ghent facades, 90-day measurements | END maps under by more than 14 dBA at a courtyard; model strongly improved |

Sources: Schissler 2014 (full text); [Schissler 2021](https://research.facebook.com/file/4521255027934795/Fast-Diffraction-Pathfinding-for-Dynamic-Sound-Propagation.pdf); [Remillieux 2012](https://consensus.app/papers/details/6b8196c9fc385b87b763b6f9e06e177c/?utm_source=claude_desktop); [Riegel 2019](https://consensus.app/papers/details/e8ffde29ba555a1b9de471bd8a38b223/?utm_source=claude_desktop); [Yunus 2023](https://iris.polito.it/bitstream/11583/2978683/1/2023YunusAST.pdf); [Yang 2025](https://consensus.app/papers/details/58b4673d283651ee90e8c5a204973d59/?utm_source=claude_desktop); [Wei 2014](https://publications.lib.chalmers.se/records/fulltext/206833/local_206833.pdf)

**High-fidelity references used in UAM**
- FW-H coupled to BEM accelerated by an adaptive black-box fast multipole method, for a VTOL in an urban area — [Barbarino et al. 2021, Aerospace Sci. Tech.](https://consensus.app/papers/details/84e3abd66836550f8c5ee1615122bda7/?utm_source=claude_desktop)
- FW-H coupled to FEM with chimera grids for a helicopter in a realistic urban environment — [Barbarino et al. 2010](https://consensus.app/papers/details/0fe7ec3f03895a30ac25fce6838a928a/?utm_source=claude_desktop)
- These are candidate references for validating any geometric model's building terms.

### Inferences
**Ranking for NICE-OPS** (inference, weighing the evidence above)
1. **Best fit: deterministic image-edge with a single-edge diffraction coefficient (Pierce/UTD-class, or UDFA-style filters), low order.** Up to 1-2 diffractions and about 3-4 specular facade reflections, plus a diffuse higher-order term (Wei/Can/IMCRI style). Computed per receiver as a direction-indexed urban transfer table, and coupled to the existing stratified eigenrays VSM-style through the arrival direction. Accuracy class: about 1-3 dB lit and first-order shadow against BTM, and about 3 dB standard deviation against FDTD for the canyon term (Wei 2014). This is deterministic, so it suits golden tests and the windowed peak-search bounds.
2. **BTM** as the offline reference used to calibrate or verify (1), via BSD EDtoolbox. Too costly for the run itself, fine for a precompute or a test oracle.
3. **Stochastic ray or particle tracing with diffuse reflections:** handles scattering naturally, but Monte Carlo noise of several dB at affordable ray counts, and it needs UTD or UBD bolted on for diffraction. Acceptable only inside an offline per-receiver precompute.
4. **Gaussian beam tracing:** attractive because it natively unifies refraction and reflection and is GPU-friendly. But the published UAM implementations lack diffraction, producing 20-35 dB shadow errors, and cost tens of seconds per case on 36 cores. NICE-OPS's existing eigenray machinery already covers GBT's refraction advantage.
5. **Maekawa or ISO 9613-2 style charts:** cheapest, but documented errors over 10 dB at shielded urban positions (Wei 2014; Wei et al. 2015 IOA). Use them only as a fallback.

**Regardless of method**
- Below about 100 Hz, which covers the rotor fundamentals, geometric urban models are outside their validity, and the low bands need a wave-based or no-correction treatment.
- Turbulence scattering into shielded canyons is essential at high frequency (Wei 2014; Van Renterghem 2006). NICE-OPS already has a turbulence-scattering model for refractive shadows (`ShadowScattering`) that could be extended.

### Gaps
- No head-to-head study of image-edge, ray tracing and GBT on the *same* urban scene with an elevated source against a wave or measurement reference was found. The ranking above is synthesized across studies with different scenes.
- RAVEN's outdoor validation and license (RWTH, proprietary as generally known) was not researched in this pass.
- Steam Audio (Valve; Embree or Radeon Rays based) was not researched. It is believed to be Apache 2.0 open source since 2024, but that is unverified here.
- I-Simpa/SPPS (the sound-particle code family used for studies like Can et al. 2015) license was not verified.
- Licenses of the RWTH ITAGeometricalAcoustics, ART and pigeon repositories were not verified (bot wall).
