# Acceleration, precomputation and emulation of urban (building-aware) propagation at footprint scale

Context for the reader: NICE-OPS evaluates ~10^4 observers x ~10^4 states (~10^8 source-receiver pairs) at
~20-150 ns per pair on the spectral path, integrates SEL over the 10 dB-down window, and culls states with
per-direction upper bounds. Two candidate architectures: (1) per-receiver, direction-indexed transfer tables
built by reciprocity (trace outward from each receiver) for sources well above rooftops; (2) per-pair 3-D
tracing (Embree/OptiX) for near-field cases such as vertiports. Numbers below are from the cited sources;
anything computed here is marked as an inference.

## 1. Reciprocity-based precomputation (receiver-side transfer functions, Project Acoustics/Triton, sonAIR, Gao et al. 2023)

### Takeaway
Every production system that makes building-aware propagation cheap at runtime does so by precomputing at the
*constrained* end (listener/receiver probes) and storing smooth, perceptually or energetically reduced
parameters rather than raw responses: Microsoft's parametric coding stores ~100-120 KB per listener probe
(80-160 MB per game map) and decodes in ~100 us per source, and Empa's sonAIR precomputes band attenuations
at virtual source positions on an altitude-adaptive cell grid and interpolates by aircraft position. Neither
stores a pure direction-only table; that is a NICE-OPS-specific simplification justified only in the far field
(Section 2).

### Cited Findings
- **Parametric wave field coding (Raghuvanshi & Snyder, ACM TOG 2014).** The acoustic field is "a chaotic 7D
  function of time and the positions of source and listener"; coding impulse responses as a few perceptual
  parameters from the energy decay makes the parameter fields "spatially smooth", compressible with a lossless
  PNG-like scheme, removing two of the seven dimensions; entire game maps fit "within 100MB" and run-time
  decoding takes "100 us per source". — [ACM DL](https://dl.acm.org/doi/10.1145/2601097.2601184)
- **Parametric directional coding (Raghuvanshi & Snyder, ACM TOG 37(4) 2018) — the closest analogue of a
  receiver-side direction-aware table.** — [paper PDF](https://www.microsoft.com/en-us/research/wp-content/uploads/2018/10/paramd.pdf)
  - Uses *reciprocal* simulation: probes are placed at the (constrained) runtime listener positions, "typical
    spacing of 3-4m", and each probe's 3-D wave simulation (90 m x 90 m x 30 m volume) yields the response for
    all runtime *source* positions in that volume. Arrival direction at the listener is obtained by simulating
    three dipole sources at the probe (centered differences) plus one monopole: four simulations per probe.
  - Cost: at 1 kHz a simulation has ~120 million cells; the raw field over 0.5 s would be 5.5 TB ("30 hours just
    for disk I/O"); their streaming encoder does it in 5 hours with 40 GB RAM and no disk. Probes "can be
    computed in parallel with each one taking 5-6 hours at 1kHz and 20 minutes at 500Hz on a single 8-core
    machine". Grid spacing 12.75 cm at 1 kHz.
  - Storage: Citadel 780 probes = 80 MB; HighRise 1350 probes = 160 MB (inference: ~100-120 KB per probe).
    Parameter fields over source position are smoothed and subsampled on a 1.5 m grid, quantized, delta-coded
    and zlib-compressed.
  - Encoded parameters per (source cell, listener probe): onset delay, initial loudness, primary arrival
    direction, reflection delay, six directional reflection loudnesses, decay time. Quanta
    {2 ms, 2 dB, (6.0 deg, 2.8 deg), 2 ms, 3 dB, 3 (log1.05 steps)}.
  - Directional reflection energy uses six world-axis cosine-squared lobes: "this directional basis forms a
    partition of unity which preserves overall energy, and does not ring to the opposite hemisphere like
    low-order spherical harmonics." Arrival direction is encoded as (elevation, azimuth) with shortest-arc
    deltas to avoid 0/2pi jumps; near the source (<2 m, unoccluded) it is replaced by line of sight because the
    flux direction is singular there.
  - Initial-arrival loudness is aggregated over a 10 ms window and reflections over 80 ms; the perceptual fusion
    windows (1 ms precedence) set what is stored.
- **Reciprocal Latent Fields (Seute, Vasudev, Richan, Buffoni; SIGGRAPH 2026 conf.)**: a volumetric grid of
  trainable latent embeddings decoded by a *symmetric* function (reciprocity built in) to predict the acoustic
  parameters of precomputed impulse responses; "reduces the memory footprint by several orders of magnitude"
  versus wave-coding baselines; perceptually indistinguishable from ground truth in a listening test. No
  ratios or runtimes in the abstract. — [arXiv 2602.06937](https://arxiv.org/abs/2602.06937)
- **sonAIR (Empa; official Swiss aircraft-noise model from 2026).** — [Empa sonAIR](https://www.empa.ch/web/s509/sonair)
  - Attenuation is *precomputed* with the physical model sonX "at cell corners along flight routes, where
    virtual source positions are placed", and the precomputed values are "weighted by the aircraft position
    within the cell, with cell sizes dynamically adjusted based on altitude and horizontal distance to
    receivers"; each sonAIR version carries its own attenuation database (BASIC or METEO module, with or without
    building shielding). (Search-result summary of the 2025 sensitivity study; full text was behind a captcha.)
    — [Acta Acustica 2025, aacus250029](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2025/01/aacus250029/aacus250029.html)
  - "the vectorized building dataset is exclusively used to calculate shielding effects, making it relevant
    solely for grazing sound incidence." — [Empa sonAIR](https://www.empa.ch/web/s509/sonair)
  - Building reflections (module REFLECT) are "implemented in a simplified way. Instead of calculating
    diffraction, only a visibility check is performed. Energy exchanges up to the third order" (source-reflector,
    reflector-reflector, reflector-receiver); scattering (not coherent) solution used. Cost on an i7/32 GB
    workstation: "a few minutes" per flight for BASIC with or without buildings; building reflections add
    "about 10-20 minutes"; adding forest reflections makes it "about two hours". — [Ramseier et al., FA2023](https://dael.euracoustics.org/confs/fa2023/data/articles/000587.pdf)
- **Gao et al. 2023, "Developing virtual acoustic terrain for Urban Air Mobility trajectory planning"
  (Transportation Research Part D).** Combines acoustic ray tracing with reciprocity to convert noise
  constraints at urban receivers into 3-D exclusion zones aircraft should avoid; demonstrated on two Chicago
  downtown cases; the point is to avoid recomputing footprints during trajectory planning. — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1361920923001918)
  (abstract-level only; full text 403).
- **Other reciprocity uses**: backward tracing from the listener gives sub-linear scaling in source count and a
  5x speedup over forward tracing in a multi-source system — [Schissler & Manocha, multi-source](http://gamma.cs.unc.edu/MULTISOURCE/paper.pdf);
  "Versa" (2025) creates virtual emitter samples by swapping emitter/listener poses and notes that
  emitter/listener gain patterns complicate naive reciprocity — [arXiv 2510.20602](https://arxiv.org/abs/2510.20602).
- **Heutschi (Acta Acustica u. Acustica 95(4), 2009)**: reflection gain in street gorges tabulated from "numerous
  ray tracing simulations" into look-up tables parameterized by street width, facade height, receiver-lane
  distance and gap width/position — an engineering precedent for tabulated urban gain. — [Chalmers record citing it](https://research.chalmers.se/en/publication/194530)
  (description from search summary).

### Inferences
- The Microsoft pipeline *is* the reciprocity architecture NICE-OPS option (1) proposes, but it stores the
  field over *source position* within 45 m of the probe because game sources are near. For a helicopter
  hundreds of meters away the 3-D source-position field collapses to (direction, range) and, if range enters
  only through spreading (Section 2), to direction alone -- which is what makes a 2-D table per receiver
  sufficient.
- Storage estimate for option (1) (inference, not from a source): an upper-hemisphere direction grid at ~3 deg
  is ~2,300 directions; x 24 third-octave bands x 14,400 receivers x 2 bytes (int16 at 0.01-0.1 dB) = ~1.6 GB;
  at 5 deg (~825 directions) ~0.57 GB; 1 byte at 0.25 dB halves it. Per receiver it is ~110 KB at 3 deg, i.e. an
  L2-resident table for the duration of one observer's state loop, matching NICE-OPS's observer-outer loop.
  This is the same order as Microsoft's ~100-120 KB per probe.
- Microsoft's choice of six cosine-squared lobes (partition of unity, no ringing) over low-order spherical
  harmonics is directly relevant: low-order SH of a band-energy table can go negative and ring across a
  shadow boundary, which would break a culling upper bound; non-negative bases (lobes, piecewise-linear on a
  grid) keep the max-over-cell bound trivial.
- sonAIR's altitude-adaptive cell grid over source positions is the aviation-proven alternative when
  direction-only is not valid (lower altitudes, nearer receivers): store attenuation at virtual source
  positions, coarser with altitude and horizontal distance.

### Gaps
- Could not verify the author list "Gao, Porcayo, Clarke" or any numbers (ray counts, grid, runtime) for the
  virtual acoustic terrain paper; the ScienceDirect and HKUST pages were blocked.
- sonAIR's actual cell-size law, interpolation error and database size were not obtainable (Acta Acustica 2025
  full text behind captcha; sonAIR 2018 paper on ResearchGate blocked).
- No source found for Mehra et al. 2013 (equivalent-source method for large open scenes) storage/runtime numbers;
  only its design (per-object transfer functions, analytic inter-object interactions) — [MSR page](https://www.microsoft.com/en-us/research/publication/wave-based-sound-propagation-large-open-scenes-using-equivalent-source-formulation/).
- No published *quantitative* study of how many incidence directions an outdoor building-aware transfer
  function needs was found (HRTF literature gives the analogous result; see Section 6).

## 2. Validity of far-field / plane-wave incidence for an aircraft above a city

### Takeaway
No published quantitative criterion for "far enough for a direction-only urban lookup" was found. The usable
criteria have to be derived from geometry: the error of using the direct path's spreading for a reflected path
is ~8.7 x (2d/r) dB for a reflector at distance d from the receiver and source range r, and the shift of a
diffraction shadow boundary scales with (neighbourhood size / r). Both say the lookup is good to ~0.5 dB when
r exceeds ~35x the distance to the farthest significant reflector or edge, which holds for en-route and
approach overflights but fails at vertiports.

### Cited Findings
- Empirically, building effects on aircraft noise are concentrated at grazing incidence: "Buildings have a minor
  effect in case of direct fly-overs; however, this is not the case for situations with grazing sound
  incidence, where barrier effects can lead to a substantial reduction". — [Acta Acustica 2025 (sonAIR sensitivity study)](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2025/01/aacus250029/aacus250029.html)
- Specular street-canyon modelling of low-flying aircraft: relative EPNL increase "in proportion to the ratio of
  building height to flying altitude", while "changing the street width has only a minor effect". — [Ismail & Oldham 2002](https://journals.sagepub.com/doi/10.1260/135101002320815693)
- In the Microsoft encoder the flux-derived arrival direction becomes singular as source and listener meet and
  is replaced by the line of sight within 2 m. — [paramd PDF](https://www.microsoft.com/en-us/research/wp-content/uploads/2018/10/paramd.pdf)
- ISO 16283-3 facade sound-insulation testing substitutes a loudspeaker at 45 deg incidence for road, rail or
  aircraft sources (Annex D covers aircraft) -- i.e. standard practice treats aircraft noise on a facade as a
  single incidence direction. (Loudspeaker distance not obtained.) — [ISO 16283-3 summary](https://standards.globalspec.com/std/9961017/iso-16283-3)

### Inferences
- **Spreading error.** By the triangle inequality a reflected path is never shorter than the direct one. A
  facade at distance d from the receiver lengthens the path by at most 2d, so a direction-only table that
  applies the direct spreading to a reflected arrival over-predicts it by at most 20 log10(1 + 2d/r):
  0.56 dB at d/r = 1/30, 1.1 dB at d = 20 m and r = 300 m, 0.35 dB at r = 1 km. Reflected arrivals typically
  carry less energy than the direct, so the total-level error is smaller.
- **Direction error across the neighbourhood.** Scatterers within radius a of the receiver see the source
  from directions differing by ~a/r rad (5.7 deg for a = 50 m, r = 500 m). The table must therefore be built
  for, and smoothed over, at least that angular scale, or a rooftop shadow boundary lands in the wrong place.
- **Diffraction.** For an edge at r_r from the receiver and source at r_s, the path-length difference (and so
  the Fresnel number/Maekawa attenuation) departs from its plane-wave limit by a relative amount of order
  r_r/r_s; same criterion.
- **Practical criterion (proposal).** Direction-only is acceptable when r_source >= ~30-35 x (distance from the
  receiver to the farthest building edge or facade that contributes within ~10 dB). For ground observers in
  a city with ~20-50 m street scales that is r >= ~0.6-1.7 km for 0.5 dB. Below it, add range as a second
  table axis (log-spaced, sonAIR-style), or fall back to per-pair tracing (option 2). A (direction,
  log-range) table with ~4-8 range nodes would cover most approach geometry.

### Gaps
- No paper was found that quantifies plane-wave-versus-point-source error for building scattering of aircraft
  noise as a function of altitude; the criterion above is derived, not measured, and should be checked with a
  per-pair reference (image sources + UTD/BTM) on a city block before it is relied on.

## 3. GPU and BVH ray-tracing throughput for acoustics

### Takeaway
Raw ray casting is not the bottleneck: CPUs do 10^7-10^8 incoherent rays/s per socket (2014 hardware) and RTX
GPUs 10^8-10^9; the expensive part is *path validation* (image sources and diffraction-edge search), which in
the best interactive systems costs 0.14-7.5 ms per source-listener pair in city-scale scenes. At 10^8 pairs
that is 10^4-10^6 s, so per-pair 3-D tracing is viable only for small near-field regions (vertiports) or with
receiver-side precomputation.

### Cited Findings
- **Embree (Wald et al., ACM TOG 2014; Apache-2.0).** Full path tracing including shading: 9-117 Mrays/s across
  five scenes on Haswell, 2x SandyBridge and Xeon Phi; primary rays 27-464 Mrays/s. A functionally similar
  OptiX path tracer on a GTX Titan ran 36-80 Mrays/s on the same scenes; Karras et al. report 350 Mrays/s
  path tracing on a GTX Titan. Single-ray SIMD kernels give 2x (Haswell) to 8.6x (Xeon Phi) frame-time speedup
  over scalar traversal. — [Embree paper](https://www.sci.utah.edu/~wald/Publications/2014/embree/embree.pdf);
  license — [GitHub RenderKit/embree](https://github.com/RenderKit/embree)
- **RTX hardware.** On an RTX 2070 (Vulkan/DX12) measured 92-806 Mrays/s primary and ~50-635 Mrays/s for
  secondary (random) rays across Sponza (66K tris), Crytek Sponza (262K) and Hairballs (224M). — [Examination of the Nvidia RTX, CEUR Vol-2485](https://ceur-ws.org/Vol-2485/paper3.pdf).
  Turing RT cores are quoted at "10+ Giga-Rays per second" vs ~1.1 for high-end Pascal (search summary of
  the same literature). — [search source list incl. CEUR paper](https://ceur-ws.org/Vol-2485/paper3.pdf)
- **OptiX for acoustics**: an OptiX ray tracer for underwater sound was "up to 310 times faster" than the CPU
  Fortran code Bellhop, though the RT-core share of the gain could not be isolated. — [Ulmstedt & Stalberg thesis](https://www.diva-portal.org/smash/get/diva2:1352170/FULLTEXT01.pdf)
- **GSound (Schissler & Manocha, AES 41st)**: with 1000 visibility rays, 4 reflection orders and 1 diffraction
  order the system found 1,566 paths (indoor), 11,642 (outdoor) and 76,088 (cathedral). — [GSound](https://www.semanticscholar.org/paper/GSound:-Interactive-Sound-Propagation-for-Games-Schissler-Manocha/dff531a59e38071602c0478d23cc06e1986d4129)
- **Multi-source backward tracing (Schissler & Manocha, TOG)**, 3.5 GHz 4-core CPU, 4-wide SIMD BVH: City scene
  206,976 triangles, 50 sources in 24 clusters, propagation 47.2 ms per frame producing 27,711 paths; Sibenik
  75K tris/18 sources 49.2 ms (208,832 paths); Tradeshow 200 sources 182.7 ms (871,982 paths). — [MULTISOURCE paper](http://gamma.cs.unc.edu/MULTISOURCE/paper.pdf)
- **High-order diffraction + diffuse reflections (Schissler, Mehra, Manocha, TOG 2014)**: precomputed global
  edge-visibility graph plus wavelength-dependent mesh simplification (edges reduced to ~30-90% of the
  original). Small City 89,792 tris / 64,853 edges -> 25,656 after simplification, visibility graph 31.5 s,
  5.18 MB; 8 sources, 10 specular, 6 diffuse, 6 diffraction orders: 19.3 ms/frame. Large City 254,903 tris /
  222,680 edges -> 70,605, graph 299.5 s, 23.8 MB; 14 sources with 4 diffraction orders: 71.0 ms/frame. 4-core
  i7-4770K, 8 threads. Accuracy comparable to path tracing with 10,000 rays and 10 reflection orders. — [HIGHDIFF paper](https://gamma.cs.unc.edu/HIGHDIFF/paper.pdf)
- **Fast diffraction pathfinding (Schissler et al., SIGGRAPH 2021)**: silhouette-edge heuristic + edge
  visibility graph + A*; diffraction "using 0.14ms - 7.5ms per source", preprocessing under a minute,
  10th-order diffraction "2.7 to 568 times faster per source" than prior work. — [paper PDF](https://research.facebook.com/file/4521255027934795/Fast-Diffraction-Pathfinding-for-Dynamic-Sound-Propagation.pdf)
- **Steam Audio** (occlusion, transmission, reflections, pathing; probe-based baking with layers for reverb,
  static-source reflections and pathing) is fully open source under **Apache-2.0** since v4.5.2 (Feb 2024). — [Phoronix](https://phoronix.com/news/Steam-Audio-SDK-Fully-Open); [FOSS Force](https://fossforce.com/2024/02/valve-releases-steam-audio-as-open-source/)
- **Gaussian beam tracing for UAM at vertiports (Yunus, Casalino, Avallone, Ragni; Aerospace Sci. Tech. 139,
  108410, 2023)**: replaces straight rays to include complex directivity, 3-D terrain and wind; buildings raise
  ground-level noise by 5 dB in illuminated zones and create shadows; vs straight rays "reduces the prediction
  error by 5 dB in the illuminated zones and 35 dB in the terrain shadow zones". No runtime given in the abstract. — [TU Delft portal](https://research.tudelft.nl/en/publications/efficient-prediction-of-urban-air-mobility-noise-in-a-vertiport-e/)

### Inferences
- Cost per pair, option (2): taking ~1 ms per source-listener pair (between the 0.14-7.5 ms diffraction figures
  and the ~1-5 ms/source of the city frames) gives 10^8 pairs ~ 10^5 s (~28 h) on a 2014-2021 4-core CPU,
  three to four orders of magnitude above NICE-OPS's current ~2-15 s of loop CPU. GPU throughput might buy 10-100x
  but not 10^4x. Option (2) therefore has to be restricted to a small observer set (e.g. a vertiport
  neighbourhood of 10^2-10^3 receivers) or to the few states inside each SEL window, using NICE-OPS's existing
  bounds to cull.
- Option (1) moves the same ray work to precompute: per receiver, trace ~10^4-10^5 rays outward (2,300
  directions x tens of rays per direction for edge/reflection sampling) -- at 10^7-10^8 rays/s that is
  ~1-10 ms per receiver, ~15-150 s for 14,400 receivers before path validation; with the 1-7 ms validation per
  direction it grows to minutes-hours, still a one-off per city and receiver set.
- License map for an MIT host: Embree (Apache-2.0), Steam Audio (Apache-2.0) are compatible; NoiseModelling
  (GPL-3.0, Section 5) is benchmark-only; OptiX is a proprietary NVIDIA SDK (not verified here; treat as an
  optional, separately licensed backend).

### Gaps
- No published throughput figure for a *city-block image-source + edge-diffraction* solver expressed per
  source-receiver pair for outdoor aircraft geometry (high source, ground receiver) was found; all timings
  above are game/VR frames with sources near the listener.
- No measured numbers for NVIDIA VRWorks Audio were found.
- GSound's license and current availability were not determined.

## 4. Surrogate / ML models for building shielding and urban noise maps

### Takeaway
Image-to-image surrogates (CNN, GAN, normalizing flows) reproduce 2-D urban noise maps at 0.3-2.6 dB RMSE or
0.65 dB MAE with 1,800-2,000x speedups, but they are trained on fixed 2-D scenes and source parameterizations
and give no error bound; the only surrogate with a *certified* worst-case bound (monotone NNs, 1.2-5.6 dB) does
not model buildings. For a culling-based exact pipeline a learned surrogate is a poor fit except as a
candidate generator.

### Cited Findings
- UAM noise in an urban environment: modified CNN trained/evaluated on 45,000 2-D noise maps; RMSE 2.56 dB
  against ray tracing; >1,800x faster. — [ResearchGate record (abstract via search)](https://www.researchgate.net/publication/377186372_Deep_learning_based_prediction_of_urban_air_mobility_noise_propagation_in_urban_environment)
- Conditioned normalizing flows (Full-Glow), Eckerle, Spitznagel, Keuper (2025): 256x256 urban sound-pressure
  maps in real time on one RTX 4090 (102 ms per map per search summary), ">2000 times" faster than the
  reference solver; baseline NLoS MAE 0.65 dB; up to 24% better NLoS accuracy than prior deep models on
  diffraction and reflection datasets. — [arXiv 2510.04510](https://arxiv.org/abs/2510.04510)
- GAN surrogate for traffic noise mapping (Yang, Han, Lu, Zhang 2024): RMSE 0.3024 dB(A), SSIM 0.8528;
  integrated into Grasshopper for early design. — [arXiv 2405.13227](https://arxiv.org/abs/2405.13227)
- CNN for accelerated traffic noise mapping: RMSE 1.8 dB, 91.3% of predictions within 3 dB (search summary;
  not fetched). — [Scientific Reports](https://www.nature.com/articles/s41598-026-46367-3)
- Certified, physics-informed monotone NNs for eVTOL noise (Park, Elfar, Fleming, Shoukry 2025): monotone in
  speed, RPM, altitude, distance; azimuth split into 22 sectors at 1 dB tolerance; trained on 226,935
  FLOWUnsteady + PSU-WOPWOP samples on a 4 km x 4 km grid at 100 m; certified worst-case bounds 1.23-5.62 dB
  per sector with active sampling (8.83-19.54 dB uniform sampling). No building/urban propagation. — [arXiv 2509.20306](https://arxiv.org/html/2509.20306v1)
- Neural-operator surrogates for outdoor acoustics with parameterized wind and terrain profiles exist (JASA
  160(3), 2653); details not obtained (403). — [JASA](https://pubs.aip.org/asa/jasa/article/160/3/2653/3405543/Neural-operator-surrogates-for-outdoor-acoustics)
- A multi-edge-diffraction simplified analytical model for shielded urban locations deviates from full-wave
  simulations by <2 dB for 125-1000 Hz octaves (search summary). — [PubMed 26627751](https://pubmed.ncbi.nlm.nih.gov/26627751/)

### Inferences
- The reported errors (0.65-2.6 dB) are of the same order as the effects being modelled (3-5 dB building
  reflection gains, Section 7), and are averages over maps; NICE-OPS's exact windowed search needs bounds that
  are never exceeded, which none of the urban surrogates provide.
- The useful idea to borrow is the certified-monotone construction: a per-receiver table with a stored
  max-over-cell value is already a certified bound, at zero training cost.

### Gaps
- The UAM CNN paper's authors, journal, input resolution and frequency treatment could not be verified (page
  blocked).
- No surrogate was found that handles a moving elevated source over 3-D buildings with per-band output.

## 5. Engineering simplifications in noise mapping (CNOSSOS-EU, Nord2000, NoiseModeling, sonAIR)

### Takeaway
Strategic noise mapping makes millions of receivers tractable by cutting path search, not by precomputation:
reflection order 1, lateral diffraction dropped for extended sources (explicitly including aircraft), 10 m
grids at 4 m height, tiling, and "efficiency settings" accepted when the 95% confidence interval stays within
1 dB of the full calculation. Aircraft mapping still ignores buildings altogether.

### Cited Findings
- CNOSSOS-EU propagation is ray-based over all source-receiver paths; reflected paths from vertical obstacles
  follow specular reflection in the horizontal plane; laterally diffracted paths from vertical edges "usually
  have negligible contributions" for extended sources ("road, railway and aircraft") "and can therefore be
  omitted" (JRC draft, via search summary). — [CNOSSOS-EU JRC draft](https://sicaweb.cedex.es/wp-content/uploads/2021/08/CNOSSOS-JRCdraft.pdf)
- Irish EPA guidance (May 2025): lateral (vertical-edge) diffraction under CNOSSOS-EU "only applies to
  industrial sources"; for multiple reflections between two walls "at least the first reflection shall be
  considered"; Guidance Note 35 recommends reflection order = 1 for strategic mapping, accepting
  underestimation "behind buildings or within enclosed courtyards", because raising it to 2 or more
  "increase[s] calculation times by an order of magnitude"; 2-3 for localized studies. Grid receivers on a
  regular 10 x 10 m grid; assessment receivers 4 m above ground; facade receivers per segment (Case 2a).
  Efficiency settings (reflection search radius, minimum source-receiver distance, reflection count, tiling)
  should keep the 95% confidence interval within 1.0 dB of the base case, tested on a 5 x 5 km model with a
  central 1 x 1 km calculation area. Aircraft noise under the Directive uses ECAC Doc 29 4th ed. — [EPA Ireland guidance Part 2](https://www.epa.ie/publications/monitoring--assessment/noise/Part-2-Calculation-Methodology-(May-2025).pdf)
- NoiseModelling (Univ. Gustave Eiffel; **GPL-3.0**): implements CNOSSOS-EU (Directive 2015/996) for road and
  rail, configurable reflection order, H2GIS spatial database, third/octave bands 50-10,000 Hz. — [GitHub](https://github.com/Universite-Gustave-Eiffel/NoiseModelling).
  A production run used 4,500,417 receivers and 1,304,158,043 source/receiver pairs; line sources segmented
  equidistantly to speed ray tracing (search summary of the v3 overview paper). — [NoiseModelling v3 overview](https://www.researchgate.net/publication/343904801_Overview_of_the_Noise_Modelling_open-source_software_version_3_and_its_applications)
- Facade convention: a microphone 0.5-2 m in front of a reflecting surface takes a -3 dB correction to recover
  the incident level (search summary of CNOSSOS/ISO 1996 sources). — [CNOSSOS-EU JRC draft](https://sicaweb.cedex.es/wp-content/uploads/2021/08/CNOSSOS-JRCdraft.pdf)
- Nord2000: screens on flat ground predicted ~1 dB above reference methods (138 downwind thin-screen cases,
  25-400 m). — [Nord2000 validation AV 1117/06](https://forcetechnology.com/-/media/force-technology-media/pdf-files/projects/nord2000/nord2000-validation-of-the-propagation-model-av-1117-06.pdf)
- Aircraft mapping: "neither shielding by buildings nor reflections by buildings, forests or cliffs are
  considered in noise mapping"; Empa recommends detailed ground properties for general mapping and buildings
  only for specific cases, because buildings add "substantial computational effort" and give "highly
  variable (i.e., non-smooth) noise contours". — [Acta Acustica 2025](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2025/01/aacus250029/aacus250029.html); [Ramseier et al. FA2023](https://dael.euracoustics.org/confs/fa2023/data/articles/000587.pdf)

### Inferences
- The engineering methods' main trick (first-order reflections only, no lateral diffraction for aircraft,
  energy addition) maps naturally onto a per-receiver direction table: first-order specular images plus
  vertical-plane edge diffraction is a small, enumerable set of paths per receiver per direction.
- The 1 dB / 95% CI acceptance criterion is a reasonable external benchmark for whatever NICE-OPS's table
  discretization costs.

### Gaps
- Could not obtain the END Annex I text on "incident sound" at the most exposed facade (and the Directive's
  own +3 dB statement) from a primary source in this session.
- Nord2000's treatment of urban reflections and its runtime per path were not found.

## 6. Interpolation over source position and angular smoothness of urban transfer functions

### Takeaway
Energy-domain (band-averaged, incoherent) transfer quantities are smooth enough to subsample coarsely --
Microsoft subsamples over source position at 1.5 m and quantizes direction at 6 x 2.8 deg; sonAIR interpolates
attenuation between cell corners -- whereas the coherent field is hopeless (spherical-harmonic order ~kr). The
sharpest feature is a diffraction shadow boundary, whose angular transition width is about sqrt(lambda/d_edge)
radians, so a ~2-3 deg direction grid resolves edges more than ~20 m away up to ~1 kHz.

### Cited Findings
- "These raw parameter fields are quite smooth, affording spatial coherence and good compression"; they are
  smoothed and subsampled on a 1.5 m grid, and direction quantized at (6.0 deg, 2.8 deg). — [paramd PDF](https://www.microsoft.com/en-us/research/wp-content/uploads/2018/10/paramd.pdf)
- Parametric coding exists precisely because the raw IR field is "chaotic" and hard to interpolate; perceptual
  parameters are smooth. — [ACM DL 2014](https://dl.acm.org/doi/10.1145/2601097.2601184)
- sonAIR: attenuation at virtual source positions on cell corners, weighted by the aircraft's position in the
  cell, cell size adapted to altitude and horizontal distance. — [Acta Acustica 2025](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2025/01/aacus250029/aacus250029.html)
- HRTF analogue: a "well-defined finite number of spherical harmonics" captures the spatial variation (the
  spatial dimensionality), and an optimal sampling scheme follows from it (Zhang et al., IEEE TASLP 20(2) 2012). — [ANU record](https://researchportalplus.anu.edu.au/en/publications/on-high-resolution-head-related-transfer-function-measurements-an/)
- Spherical-harmonic fits of sparse directivity data reached <1 dB mean error with 89% fewer measurement points
  (search summary). — [search result, BYU](https://physics.byu.edu/docs/publication/3493)
- Certified eVTOL surrogate partitions azimuth into 22 sectors at a 1 dB "bounded azimuthal variation"
  tolerance (free-field source directivity, not buildings). — [arXiv 2509.20306](https://arxiv.org/html/2509.20306v1)
- Empirical non-smoothness: building effects make aircraft-noise contours "highly variable (i.e., non-smooth)". — [Ramseier et al. FA2023](https://dael.euracoustics.org/confs/fa2023/data/articles/000587.pdf)

### Inferences
- Coherent transfer: the standard HRTF result is truncation order N ~ k a for a scatterer region of radius a
  (not quoted in the fetched abstract). For a = 20 m at 1 kHz (k = 18.3 /m) that is N ~ 370, ~1.4 x 10^5
  coefficients per band: storing complex transfer functions is impractical; store band *energies*.
- Shadow-boundary width: the Fresnel transition at an edge a distance d_edge from the receiver spans about
  sqrt(lambda / d_edge) rad in source direction: 7.5 deg at 1 kHz for d_edge = 20 m, 24 deg at 100 Hz, 1.7 deg
  at 4 kHz for 100 m. A 2-3 deg grid with bilinear interpolation is adequate for ~<=1 kHz and near edges;
  high bands behind distant edges will be under-resolved (step-like), which affects accuracy but not the bound
  if the stored per-cell maximum is used for culling.
- Moving source: because NICE-OPS integrates energy over the SEL window, the integrand needs only band
  energies at each state; interpolation in direction between states is not needed (each state looks up its own
  direction). Temporal smoothness matters only for the 10 dB-down window detection, which is already robust to
  steps by the "outermost crossing" rule.

### Gaps
- No measured or simulated study of urban building-scattering band-energy smoothness versus incidence direction
  for aircraft was found; the Fresnel-width argument is derived.

## 7. Bounding urban gain for culling

### Takeaway
Measured and modelled building *gains* for aircraft noise are modest -- +3-5 dB from facade reflections, up to
~7-9 dB in illuminated regions for sonic booms -- while losses are large (-12 dB shielding, 8-14 dB facade
contrasts). A rigorous per-receiver bound is available for free: every reflected path is at least as long as
the direct one, so K specular paths add at most 10 log10(1 + K) dB in energy (20 log10(1 + K) coherent), and a
precomputed table's per-direction maximum is an exact bound for the table model.

### Cited Findings
- sonAIR, Meiringen airbase: building shielding reduced the annual rating level by >2 dB over large areas and
  "up to 12 dB directly on the building facades averted from the runway"; building reflections increased it
  "locally, close to the facades, by up to 3-4 dB" on the sides facing the runway. — [Ramseier et al. FA2023](https://dael.euracoustics.org/confs/fa2023/data/articles/000587.pdf)
- Vertiport Gaussian beam study: buildings increased ground-level noise by 5 dB in illuminated zones. — [Yunus et al. 2023](https://research.tudelft.nl/en/publications/efficient-prediction-of-urban-air-mobility-noise-in-a-vertiport-e/)
- Sonic boom over urban areas: perceived levels increase "by at most 7-9 dB in illuminated regions of buildings
  compared to the flat ground case" (JASA 152(6) 3323, 2022; search summary). — [JASA](https://pubs.aip.org/asa/jasa/article/152/6/3323/2839538/Sonic-boom-reflection-over-urban-areas)
- Schiphol courtyard experiment (shipping-container courtyards, ten microphones): differences between positions
  of 8-9 dB(A) for straight facades, up to 14 dB(A) with a slanted facade and inset. — [Lugten, Internoise 2022](https://research.tudelft.nl/en/publications/assessing-the-influence-of-street-canyon-shape-on-aircraft-noise-/).
  Mean difference between facades with and without line of sight: -1.3 to 5.0 dBA for arrivals, 8.7 to 13.6
  dBA for departures (Building and Environment 2024 version, search summary). — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0360132324002427)
- Street canyons: reflections between opposite buildings can raise levels by up to 12 dB for a rooftop
  source; diffusing fittings up to 10 dB in unfavourable configurations (search summary). — [search result list](https://www.researchgate.net/publication/285673582_The_effect_of_street_canyon_design_on_traffic_noise_exposure_along_roads)
- A facade microphone at 0.5-2 m takes -3 dB to recover incident level (energy doubling by the facade). — [CNOSSOS-EU JRC draft](https://sicaweb.cedex.es/wp-content/uploads/2021/08/CNOSSOS-JRCdraft.pdf)

### Inferences
- **Path-count bound.** Reflected path length >= direct distance (triangle inequality), so each image source's
  spreading factor is <= the direct path's; with reflection coefficients |R| <= 1 and diffraction factors <= 1,
  K additional paths give at most +10 log10(1 + K) dB incoherent (+20 log10(1 + K) coherent) over free field
  at the same straight-line range. First-order reflections from two canyon walls plus ground: K = 3 -> +6 dB
  incoherent. This holds for any source position, so it can be a single per-receiver constant.
- **Table bound.** If the table stores per (receiver, direction cell, band) the maximum gain over the cell
  (or the grid node values with bilinear interpolation, which cannot exceed the max node), then
  max_over_directions(table) per band is an exact bound for the table model, and can be fed to NICE-OPS's
  existing bound machinery exactly as the propagation table's per-band maximum already is
  (`path_gains` in the bound setup). A per-receiver-and-direction-bin max would keep DirectionBounds' tightness.
- Because measured gains are small (3-5 dB) relative to losses (up to 12-14 dB), the urban model will mostly
  *lower* levels; the culling bound should take max(table gain, 0) per bin so shielded states are not
  under-bounded when a reflection rescues them, as the sonAIR authors note reflections matter most "in shielded
  situations without direct sound incidence".

### Gaps
- No published per-receiver maximum-gain statistics for 3-D city models with elevated sources were found.
