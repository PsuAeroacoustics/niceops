# Validation data, benchmarks, inputs and open-source code for an urban street-canyon propagation model (elevated aircraft sources)

Scope: what exists (as of Oct 2026) to build and validate a per-third-octave, energy-based urban canyon propagation model in NICEOPS (MIT; GPL code only as an external executable for benchmarking). Licenses were checked against the repository or data portal itself wherever marked "verified".

## 1. Measured datasets of outdoor sound in canyons/courtyards from elevated or moving sources

### Takeaway
The single most relevant measured dataset (Lugten et al., 2383 real aircraft flyovers into three full-scale container courtyards at Schiphol) is published but I found no public data deposit. The best *downloadable* urban measurement set is RWTH Aachen's new EDGE scale-model database (Zenodo, CC BY-SA 4.0, Aug 2026). It has CAD geometry and measured responses for about 31 generic urban scenes, but its sources are not aircraft-like elevated moving sources. Measured helicopter, drone and UAM data inside real street canyons is essentially absent; NASA's 2025 UAM noise review lists it as an open gap.

### Cited Findings

**Lugten et al. – full-scale courtyard field lab, Schiphol (aircraft flyovers)**
- The test site is three adjacent courtyards enclosed by stacked shipping containers. The stacking differs per courtyard (slanted roofs, building insets). Sound levels for **2383 aircraft flyovers** were analysed over five months at **ten microphone positions**, for arrivals and departures. Front/rear differences range from **11 dB(A)** (straight facades) to **14 dB(A)** (slanted roof plus overhang). Published as "Assessing the influence of street canyon shape on aircraft noise: Results from measurements in courtyards near Amsterdam Schiphol Airport", *Building and Environment* vol. 255, May 2024 — [ScienceDirect S0360132324002427](https://www.sciencedirect.com/science/article/abs/pii/S0360132324002427). An Inter-Noise 2023 version also exists — [Ingenta](https://ince.publisher.ingentaconnect.com/content/ince/incecp/2023/00000265/00000001/art00055;jsessionid=2bvoeo67ef7hs.x-ic-live-03)
- The field lab ("Urban Comfort Lab") was built in 2021 near the Kaagbaan, Hoofddorp, from 120 shipping containers forming three courtyards — [Geluid 2023-4](https://mp.nl/sites/default/files/2023-12/GELUID_2023_4.pdf); [AMS Institute](https://www.ams-institute.org/urban-challenges/resilient-cities/noise/); [TU Delft story](https://www.tudelft.nl/en/stories/articles/building-differently-to-combat-aircraft-noise)
- Earlier Lugten, Steemers & Kang pilot: 45 flyovers recorded in front of and behind real buildings near Schiphol. Mean shielding was about **11 dB(A) for landings and 14 dB(A) for take-offs**. In octave bands it was 0–7 dB below 125 Hz and 8–14 dB above, for landings — [Cambridge repository PDF](https://api.repository.cam.ac.uk/server/api/core/bitstreams/044505b4-6042-47fe-a02f-a6ac4a5f70b9/content)
- That paper cites two older studies of street-canyon effects on aircraft/rotorcraft noise:
  - Donavan, "Model study of the propagation of sound from V/STOL aircraft into urban environs", MSc thesis, MIT, 1973
  - Ismail & Oldham, "The effect of the urban street canyon on the noise from low flying aircraft", *Building Acoustics* 9(3):233–251, 2002

  Both are cited in the reference list of the [Cambridge repository PDF](https://api.repository.cam.ac.uk/server/api/core/bitstreams/044505b4-6042-47fe-a02f-a6ac4a5f70b9/content).
- Data availability: searches of 4TU.ResearchData (API query "Lugten" returned nothing), Zenodo and the TU Delft portal found **no public deposit** of the flyover data — [4TU search API](https://data.4tu.nl/v2/articles/search); the ScienceDirect page returned 403, so the article's data statement could not be read.

**RWTH Aachen – EDGE database (Evaluation of Diffraction in Generic urban Environments)** (verified)
- Palenda & Vorländer, IHTA, RWTH Aachen. Zenodo DOI 10.5281/zenodo.21264202, v0.1.0, published 19 Aug 2026, license **CC BY-SA 4.0**. Contents:
  - "scale model measurements of prototypical urban environments intended for the evaluation of sound propagation models with a focus on diffraction", measured in IHTA's anechoic chamber
  - 3D CAD models of each scene
  - source descriptions as SOFA/WAV directivity (S1, S2, S3)
  - surface descriptions: carpet and MDF absorption JSON, plus two scatterer geometries in STEP/FreeCAD
  - equipment and measurement CSVs
  - ZIP 2.1 GB, PDF documentation 240 MB

  About 31 scene folders (UE00 onward), each with metadata, geometry, measurements and pictures. Funded by DFG project 499742946 — [Zenodo record](https://zenodo.org/records/21264202); README and CITATION.cff read via the [Zenodo API](https://zenodo.org/api/records/21264202)
- The companion paper presents an open-source UTD-based urban auralization tool evaluated against "benchmark scenarios": Palenda & Vorländer, "Auralizing arbitrary urban environments including diffraction", *Acta Acustica* 2026, DOI 10.1051/aacus/2026033. Its supplementary results and code are on Zenodo 10.5281/zenodo.19112715 under **CC BY-NC-SA 4.0** — [Zenodo 19112715](https://zenodo.org/records/19112715)

**TU Eindhoven – wind-tunnel scale-model database (moving atmosphere)**
- Cosnefroy, Qin, Ricci, Gillmeier, van Hooff & Hornikx, "A database of scale model measurements for urban sound propagation in a moving atmosphere and first comparisons with simulations", ICA 2022. Setup:
  - 1:40 scale in TU/e's atmospheric boundary-layer wind tunnel
  - urban canyons and solar-panel arrays
  - up to 100 m full-scale range and 2 kHz full-scale

  The portal gives no license or data-repository link; only the paper PDF is posted — [TU/e research portal](https://research.tue.nl/en/publications/a-database-of-scale-model-measurements-for-urban-sound-propagatio/)

**Classic scale-model canyon experiments (papers only, no data deposits found)**
- Hornikx, Forssén & Kropp, "Scale model measurements and numerical modelling for directly exposed and shielded urban street canyons", Euronoise 2006 (Tampere) — [Chalmers](https://research.chalmers.se/en/publication/22994)
- Related Chalmers/TU/e studies:
  - "A scale model study of parallel urban street canyons" — [TU/e](https://research.tue.nl/en/publications/a-scale-model-study-of-parallel-urban-street-canyons/)
  - "Sound propagation to two-dimensional shielded urban areas: a numerical and scale model study of façade treatments" — [Chalmers](https://research.chalmers.se/en/publication/22993)
  - 2.5D equivalent-sources method for exposed and shielded canyons, JASA 122(5):2532 — [JASA](https://pubs.aip.org/asa/jasa/article-abstract/122/5/2532/409/The-2-5-dimensional-equivalent-sources-method-for?redirectedFrom=fulltext)
- Ismail & Oldham, "A scale model investigation of sound reflection from building façades", *Applied Acoustics* (2005). Facades act as plane specular reflectors with localized diffusers (windows, balconies). Canyon reverberation time decreases with increasing facade irregularity — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0003682X0400115X); [IOA proceedings PDF](https://www.ioa.org.uk/system/files/proceedings/mr_ismail_and_dj_oldham_a_scale_model_investigation_of_the_role_of_diffuse_reflections_in_urban_noise_propagation.pdf)
- Picaut: full-scale measurements in the first 50 m of rue Kervégan, Nantes (210 m long pedestrian street, about 18 m high, 7.90 m wide), plus a diffusion-equation model of street sound fields — [ASA lay paper](https://acoustics.org/pressroom/httpdocs/149th/picaut.htm); [JSV diffusion model](https://www.sciencedirect.com/science/article/abs/pii/S0022460X02951503)
- "A scale model experiment for the study of sound propagation in urban areas", *Applied Acoustics* (2001) — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0003682X00000281)

**Sonic boom / urban**
- SonicBAT (F/A-18 off Florida, plus Armstrong and KSC) was about turbulence effects on booms, used to validate turbulence models built into PCBoom. It was **not** an urban/building campaign — [NTRS 20200002482](https://ntrs.nasa.gov/api/citations/20200002482/downloads/20200002482.pdf)
- NASA measured conventional and low-amplitude N-waves outside two instrumented one-storey houses in June 2006 and July 2007, with multiple microphones per house. The data show the "building spiking" effect from diffraction around finite buildings — [sonicfield summary](https://labs.sonicfield.org/library/sonic-booms-spectral-analysis-and-diffraction-by-buildings); [Penn State](https://pure.psu.edu/en/publications/diffraction-of-sonic-booms-around-buildings-resulting-in-the-buil/)
- QSF18 (Galveston, 2018) measured at four stations across the city for community response — [sonicfield/BYU QSF18](https://labs.sonicfield.org/library/making-sonic-boom-measurements-during-adverse-weather-summarizing-byus-participa)
- Numerical urban-boom studies report street-canyon resonances whose amplitude depends on street width and boom spectrum — [sonicfield](https://sonicfield.org/library/reflection-of-sonic-boom-over-realistic-urban-areas); [Inter-Noise 2023](https://ince.publisher.ingentaconnect.com/contentone/ince/incecp/2023/00000265/00000001/art00017;jsessionid=1xsu1w5u2c2c4.x-ic-live-03)

**Drones / UAM**
- NASA/TP-20250011527 ("Urban Air Mobility Noise: 2025 Update", May 2026) on the gap:
  - "Measurements within such environments are critical for developing and validating multipath propagation tools … but can have large uncertainty bands"
  - "Defining canonical urban geometries and testing in those environments, with simple sources under controlled or well-characterized meteorological conditions, may expedite development" (Sec. 3.3.2)
  - Vertiport materials' "impedance properties are less known"

  — [NTRS PDF](https://ntrs.nasa.gov/api/citations/20250011527/downloads/NASA-TP-20250011527_v1.pdf)
- Field studies of hovering drones near buildings exist but are perception or level studies, not propagation-validation sets:
  - DJI Mavic 2 Pro / Matrice 200 hovering about 5 m from facades at 22–23 m, in residential courtyards — [Forests 17(8):907, 2026](https://doi.org/10.3390/f17080907)
  - Hovering-UAV soundscape study — [arXiv 1912.00087](https://arxiv.org/pdf/1912.00087)
  - "Maneuvering drone noise investigation in a reverberant environment" — [JASA 157(6):4468, 2025](https://pubs.aip.org/asa/jasa/article/157/6/4468/3350490/Maneuvering-drone-noise-investigation-in-a)
- UaVirBASE is a public UAV recording dataset for source localization (128 UAV recordings, 14.6 GB, at varying altitudes, distances and angles) — [Applied Sciences 15(10):5378](https://www.mdpi.com/2076-3417/15/10/5378)

**Helicopter, urban**: I found no publicly downloadable helicopter-in-street-canyon measurement set. The classic reference is NASA's 1978 urban helicopter noise requirements paper — [NTRS 19780024877](https://ntrs.nasa.gov/api/citations/19780024877/downloads/19780024877.pdf).

### Inferences
- For NICEOPS, validation needs three levels:
  - EDGE (downloadable, with geometry) for the diffraction and reflection kernels under controlled conditions
  - Lugten's flyover statistics (front/rear and per-position dB(A), octave-band shielding) as the only real-aircraft urban comparison, from the published figures and tables unless the authors share data
  - The user's own 2017 flight-test data, which is open-field, so it does not cover canyons
- EDGE's CC BY-SA 4.0 is fine for *using* the data in validation. Redistributing derived data files (e.g., checked-in goldens built from EDGE geometry) would carry share-alike obligations on those files, not on NICEOPS code. That is a judgment; check with counsel if goldens will ship.
- The Palenda tool's supplementary package is CC BY-NC-SA. Its runnable components ship as Windows wheels (pigeon_exe, pelican_exe, dovelet). It is usable as an external benchmark executable, not as linkable code.

### Gaps
- Lugten 2024 data availability statement not read (ScienceDirect 403); data presumably on request.
- EDGE scale factor, frequency range, source heights and scene list were not extracted: they are in the 240 MB PDF, which was not downloaded.
- The 2006/2007 NASA house boom data and QSF18 data were not found as public downloads.
- The NASA UAM 2025 update references no drone- or eVTOL-in-canyon public dataset.

## 2. Benchmark cases for urban/outdoor propagation used to compare codes

### Takeaway
The established urban code-comparison benchmark is Hornikx et al. (2017). It uses 2D cross-sections of typical urban geometries (screened and line-of-sight), with reference solutions from openPSTD (Fourier PSTD) and FDTD solving the linearized Euler equations. It is the urban analogue of Attenborough 1995, but I found no public deposit of its reference data. Newer candidates are EDGE (measured, 3D) and the Palenda & Vorländer stress-test scenarios.

### Cited Findings
- "Benchmark cases for urban sound propagation" (Hornikx, Trikootam, Fernandes, van Hooff, Blocken, Van Renterghem), JASA 141(5_Supplement):4035, 2017:
  - "multiple benchmark cases … two-dimensional cross-sections of typical urban geometries", source "either geometrically screened or … in the line of sight"
  - includes multiple reflections, edge diffraction and meteorology (refraction/scattering)
  - reference results from Fourier-PSTD (openPSTD v2.0) and FDTD, both LEE

  — [JASA abstract](https://pubs.aip.org/asa/jasa/article/141/5_Supplement/4035/625703/Benchmark-cases-for-urban-sound-propagation); [sonicfield summary](https://sonicfield.org/library/benchmark-cases-for-urban-sound-propagation); [ResearchGate](https://www.researchgate.net/publication/317553786_Benchmark_cases_for_urban_sound_propagation)
- Related general outdoor benchmark: "Benchmark cases for outdoor sound propagation models" (TNO) — [TNO resolver](https://resolver.tno.nl/uuid:aa12dd89-a86d-4cab-a63d-1c5665603428). It is also used to verify lattice-Boltzmann codes — [PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0147206)
- "A Platform for Benchmark Cases in Computational Acoustics" (Hornikx et al.) proposed a shared repository for such cases — [ResearchGate](https://www.researchgate.net/publication/277994015_A_Platform_for_Benchmark_Cases_in_Computational_Acoustics)
- Chalmers "Excess attenuation for sound propagation over an urban canyon" gives canonical single-canyon excess attenuation results — [Chalmers 122579](https://research.chalmers.se/en/publication/122579)
- Palenda & Vorländer (2026) evaluate their UTD tool with "benchmark scenarios intended to stress-test its capabilities". They find limits from the infinite-edge assumption, and that higher-order diffraction matters "not only in shadow regions but also around edges" — [Zenodo 19112715](https://zenodo.org/records/19112715)
- A Zenodo "Urban Sound Propagation" set has 25,000 OSM-derived 500 m × 500 m layouts in 10 cities, with NoiseModelling v4.0 maps (baseline, reflection, diffraction variants). It is anonymous and **CC BY-NC-ND 4.0** — a surrogate-model training set, simulated, not measured — [Zenodo 10609793](https://zenodo.org/records/10609793)

### Inferences
- A practical benchmark ladder for NICEOPS:
  1. Analytic wedge and thin-screen (Biot–Tolstoy/Medwin or BTM via EDtoolbox) for diffraction kernels
  2. Hornikx 2017 2D cross-sections, rerun with openPSTD as an external GPL executable, since the reference data do not appear to be published
  3. EDGE measured 3D scenes
  4. Lugten field statistics
- The NoiseModelling-generated 25k set is only useful for cross-model sanity checks; it is not a reference solution.

### Gaps
- I did not find whether the 2017 benchmark geometries and results were deposited (4TU or supplementary); they likely require the paper or authors.
- I found no Euronoise or Inter-Noise urban round-robin with published data, and no explicit building analogue of Attenborough 1995 beyond Hornikx 2017.

## 3. Facade acoustic properties (absorption, impedance, scattering/diffusion) and how they are measured

### Takeaway
Facade absorption is small (about 0.01–0.1 broadband), so the dominant uncertain parameter is scattering. ISO 17497-1 (random-incidence scattering coefficient, reverberation room, full or scale) and ISO 17497-2 (directional diffusion coefficient, free field) are the standards. The newest openly licensed data is RWTH's bidirectional scattering coefficients (CC BY 4.0, SOFA-compatible). A 2026 Acta Acustica interlaboratory comparison covers bidirectional coefficients.

### Cited Findings
- ISO 17497-1:2004 (Amd 1:2014) measures the random-incidence scattering coefficient from surface roughness in a reverberation room, "either in full scale or on a physical scale model". The result describes how far the reflection deviates from specular — [ISO 17497-1](https://www.iso.org/standard/31397.html); [SIS Amd 1:2014](https://www.sis.se/en/produkter/construction-materials-and-building/protection-of-and-in-buildings/acoustics-in-building-sound-insulation/iso1749712004amd12014/)
- ISO 17497-2:2012 measures the directional diffusion coefficient in a free field: the uniformity of the reflected polar distribution — [ISO 17497-2](https://www.iso.org/standard/55293.html)
- Heimes & Vorländer, "Bidirectional surface scattering coefficients", *Acta Acustica* 2025, DOI 10.1051/aacus/2025026:
  - a methodology for directional (incidence × reflection angle) scattering coefficients from measured and simulated reflected-pressure distributions on finite samples
  - validated on a sinusoidal test surface
  - proposes a **SOFA-compatible storage format** for propagation models

  Data, simulations and Python scripts are on Zenodo 10.5281/zenodo.15673498, **CC BY 4.0** (about 11 GB) — [Zenodo 15673498](https://zenodo.org/records/15673498)
- "Interlaboratory Comparison and Analysis of Bidirectional Scattering Coefficients for Various Surfaces", *Acta Acustica* 2026 — [Acta Acustica](https://acta-acustica.edpsciences.org/articles/aacus/full_html/2026/01/aacus250212/aacus250212.html). The page returned 403, so its content and data statement were not read.
- Cox et al., "Characterizing diffusive surfaces using scattering and diffusion coefficients" — [ResearchGate](https://www.researchgate.net/publication/238498422_Characterizing_diffusive_surfaces_using_scattering_and_diffusion_coefficients)
- Ismail & Oldham: real facades (masonry with windows, doors, balconies) behave as specular planes with localized diffusers, and more irregularity shortens canyon reverberation — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0003682X0400115X)
- A hybrid model studies the effect of facade scattering on street-canyon propagation — [Applied Sciences 9(14):2803](https://mdpi.com/2076-3417/9/14/2803/htm)
- Balcony configuration and facade material in canyons — [ASCE J. Archit. Eng.](https://ascelibrary.org/doi/10.1061/JAEIED.AEENG-1625)
- HOSANNA (EU FP7, to April 2013) studied green facades and roofs. It assumes non-vegetated facades with a **broadband absorption coefficient of 0.1**, and finds green walls on all facades give about **2–3 dB** in typical streets — [CORDIS](https://cordis.europa.eu/article/id/92989-green-defences-against-noise); [HOSANNA brochure](https://www.hosanna.bartvanderaa.com/includes/upload/DELIVERABLES/HSNNA_SUMMARY_BROCHURE_JANUARY_2013.pdf)
- Typical tabulated absorption (from a secondary chart, room-acoustics data):
  - painted brick 0.01–0.03, 125 Hz–4 kHz
  - 6 mm glass 0.18 at 125 Hz falling to 0.02 at 4 kHz

  — [Sengpiel table](https://sengpielaudio.com/calculator-RT60Coeff.htm), per the search summary (secondary source)
- In-situ facade or material absorption methods:
  - two-microphone with ambient noise — [Applied Acoustics](https://www.sciencedirect.com/science/article/abs/pii/S0003682X04002026)
  - MLS reflection method — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/0003682X93900322)
  - discrete complex image source method — [arXiv 2404.11399](https://arxiv.org/pdf/2404.11399)
  - existing building fabrics — [ResearchGate](https://www.researchgate.net/publication/288565652_In_situ_measurement_of_the_sound_absorption_characteristics_of_existing_building_fabrics)
- EDGE ships its own surface descriptions (carpet and MDF absorption JSON, two scatterer geometries), so its validation can use the same material data as the measurement — [Zenodo 21264202](https://zenodo.org/records/21264202)

### Inferences
- For an energy-based third-octave model, a facade can be parametrized as absorption α(f) (about 0.02–0.1) plus a scattering coefficient s(f). Glazing gives more low-frequency absorption from panel resonance; balconies and recesses contribute most to s.
- The Heimes & Vorländer SOFA format is a good, CC BY input format if NICEOPS ever moves beyond Lambert-plus-specular.
- Since facade α barely moves levels, most of the effort should go to geometry and scattering.

### Gaps
- I found no open database of *measured* scattering coefficients for real full-scale facades (brick, curtain wall, balconies) with a permissive license.
- CNOSSOS-EU default facade reflection values were not retrieved.

## 4. Building geometry sources, quality, height derivation, licenses

### Takeaway
Use the national LoD2 models where they exist:
- 3D BAG (Netherlands): CC BY 4.0
- Berlin LoD2: dl-de/zero-2.0
- NYC 3D model: public, license unstated

Elsewhere, Overture or OSM footprints (ODbL), with heights from USGS 3DEP lidar (public domain) in the US. Imagery-derived heights (Microsoft) are measurably poorer: 3.6–4.3 m RMSE against lidar. Google Open Buildings polygons have no heights; its 2.5D temporal product gives 4 m rasters only for the Global South.

### Cited Findings
- **Overture Maps buildings**: theme license **ODbL**. Sources:
  - OpenStreetMap (ODbL)
  - Esri Community Maps (CC BY 4.0)
  - Microsoft Global ML Building Footprints (listed by Overture as ODbL)
  - Google Open Buildings (CC BY 4.0)
  - USGS 3DEP (heights)
  - East Asian Building Data (CC BY 4.0)
  - BTN 2024 Spain (CC BY 4.0)

  — [Overture attribution docs](https://docs.overturemaps.org/attribution/). Overture says about 785 M outlines and some heights from 3DEP lidar — [Overture release notes](https://overturemaps.org/overture-january-2024-release-notes/). Its building-heights whitepaper — [Linux Foundation PDF](https://project.linuxfoundation.org/hubfs/Overture%20Maps/Building%20Heights%20Whitepaper_041423.pdf). The OvertureMaps/data GitHub repo (tooling) is MIT (verified via GitHub API).
- **Microsoft GlobalMLBuildingFootprints** (verified): README says 1.4 B buildings from Bing Maps imagery 2014–2024, "freely available … under **CDLA Permissive 2.0**". The repo LICENSE file says "CLDA Permissive 2.0" (sic). Recent releases add height estimates derived from Maxar and Vexcel imagery, e.g., 1.2 M US heights in Feb 2026, and 3.4 M for Belgium, Netherlands and Germany in Jan 2026 — [GitHub repo](https://github.com/microsoft/GlobalMLBuildingFootprints). **Conflict**: Overture's attribution page lists this source as ODbL — [Overture](https://docs.overturemaps.org/attribution/).
- **Height accuracy against lidar** (USGS; Liu et al., accepted in *PE&RS* 2026). Footprint heights were compared with 3DEP-derived heights (3DEP vertical accuracy about 10–20 cm RMSEz):
  - USA Structures: RMSE **2.4 m** (Centennial CO) and **3.1 m** (Houston); R² 0.87 and 0.86
  - Microsoft GBF: RMSE **3.6 m** and **4.3 m**; R² 0.54 and 0.60; slope 0.53 and 0.72, i.e., systematic underestimation of tall or complex buildings
  - Method: lidar classified with Point Transformer; height = median roof-point elevation − median nearby ground elevation

  — [CaGIS 2026 abstract](https://cartogis.org/docs/conferences/CaGIS_2026/abstracts/research/Liu_et_al_research_abstract_CaGIS_2026.pdf)
- **USGS 3DEP**: public-domain lidar point clouds for CONUS, available on AWS — [AWS Open Data registry](https://registry.opendata.aws/usgs-lidar/). The minimum classification is ground, non-ground and unclassified; USGS is building a deep-learning classifier (3DEPPCC) to add buildings and vegetation — [USGS 3DEPPCC](https://www.usgs.gov/centers/community-for-data-integration-cdi/science/3deppcc-automated-dl-based-point-cloud)
- **Google Open Buildings v3** polygons cover 58 M km² of Africa, South and South-East Asia, Latin America and the Caribbean. Attributes are area, confidence (0.65–1.0) and geometry, with **no height**. Dual-licensed **CC BY 4.0 or ODbL** at the user's choice — [Google Open Buildings](https://sites.research.google/gr/open-buildings/)
- **Google Open Buildings 2.5D Temporal**: annual 2016–2023 rasters of building presence, fractional count and **height** at 4 m effective resolution from Sentinel-2. Same regions, same dual license — [Earth Engine catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1?hl=en)
- **3D BAG (Netherlands)**: about 10 M buildings, default LoD2.2, generated from BAG footprints and AHN lidar. CityJSON, GeoPackage, OGC API and WMS. **CC BY 4.0** — [4TU dataset](https://data.4tu.nl/datasets/94a8f656-f7c9-4cf6-a2e7-c4c823b3db92/1); [OSM wiki on the 3D BAG import](https://wiki.openstreetmap.org/wiki/Organised_Editing/Activities/TomTom/Netherlands_-_3D_BAG_Building_Heights_Import). Tooling (verified via GitHub API):
  - 3DBAG/3dbag-pipeline: Apache-2.0
  - 3DBAG/roofer (LoD2 reconstruction): **GPL-3.0**
- **Kadaster "3D Geluid" (Netherlands)**: a national dataset made *for noise modelling*:
  - buildings (BAG footprints plus height from AHN and aerial imagery)
  - ground-factor polygons ("bodemvlakken")
  - terrain TIN and contours
  - annual releases 2020–2024, 1.7–3.4 GB national zips

  The license was not stated on the API page — [PDOK 3D Geluid API](https://api.pdok.nl/kadaster/3d-geluid/ogc/v1/collections)
- **Berlin LoD2**: about 550,000 LoD2 buildings over 890 km², CityGML with textures, license **Datenlizenz Deutschland – Zero – 2.0** — [Berlin GDI metadata](https://gdi.berlin.de/geonetwork/srv/api/records/3c7c49af-00a4-3bcd-bc00-20e7f0f1b7bf); [Cesium community post](https://community.cesium.com/t/3d-city-model-of-berlin-published-under-open-data-license/2167)
- **NYC 3D Model by Community District**: every NYC building in 2014, from DoITT's 2014 aerial survey, with roof structure and extra detail for iconic buildings. Format is .3dm (Rhino). Designated "public" with no license text in the catalog — [data.gov](https://catalog.data.gov/dataset/nyc-3d-model-by-community-district)
- A UAV-derived building height and floor-count dataset exists for one campus (UT Dallas, 2025) — [PMC12704084](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12704084/)

### Inferences
- ODbL is a database license. Using OSM or Overture geometry as *input* to a run creates no obligation on NICEOPS's MIT code. Redistributing a derived building database (e.g., a shipped test city) would require ODbL share-alike on that database. CC BY 4.0 (3D BAG, Google) and CDLA-Permissive-2.0 only need attribution. dl-de/zero needs nothing. US 3DEP is public domain.
- For canyon acoustics, height error matters through canyon aspect ratio H/W and through rooftop diffraction-edge height. Imagery-derived heights at 3–4 m RMSE with underestimation bias would systematically bias shielding. Lidar-derived (3DEP, AHN) or official LoD2 heights are preferred.
- OSM `height` and `building:levels` tags are patchy; LoD1 extrusion from levels × 3 m is the usual fallback. (Common practice, not sourced here.)

### Gaps
- The quantitative quality of OSM height tags and Overture heights was not found.
- The 3D Geluid license and CityGML LoD2 licenses for other cities (Hamburg, Vienna, Helsinki, Japan PLATEAU) were not checked.

## 5. Open-source code: ray tracing, edge diffraction, urban propagation engines (language, status, license)

### Takeaway
For an MIT C++17 code, the linkable options are:
- Embree (Apache-2.0, active)
- Steam Audio (Apache-2.0)
- Svensson's EDtoolbox (BSD, MATLAB, active; usable as a reference or port source)
- Geometry libraries: Clipper2 (BSL-1.0), earcut.hpp (ISC), Boost.Geometry (BSL-1.0), CDT (MPL-2.0)

GPL tools (NoiseModelling, I-Simpa, openPSTD, Pachyderm) are benchmark executables only. GSound is non-commercial and OptiX is proprietary. RAVEN has no public repository I could find. CGAL is GPL/LGPL per package, so avoid it or use only its LGPL parts carefully.

### Cited Findings
License and activity from the GitHub API unless noted (`pushed_at` = last push):

| Code | Language | License (verified) | Last push / status | Notes |
|---|---|---|---|---|
| Intel Embree (RenderKit/embree) | C++ | Apache-2.0 | 2026-10-01 | CPU ray-tracing kernels |
| Valve Steam Audio (ValveSoftware/steam-audio) | C++ | Apache-2.0 | 2026-03-25 | Geometric acoustics, reflections and pathing, game-oriented |
| EDtoolbox (upsvensson/Edge-diffraction-Matlab-toolbox) | MATLAB | BSD (LICENSE.md: "offered under the BSD license", © 2018 P. Svensson; text has only the two redistribution clauses, i.e. BSD-2-style, not BSD-3) | 2026-09-25 | BTM edge diffraction, frequency and time domain |
| EDB1 toolbox (older Svensson) | MATLAB | GPL v3+ (file headers) | legacy | [soundsoftware mirror](https://code.soundsoftware.ac.uk/hg/human-echolocation/file/tip/private/EDB1diff2ISES.m) |
| NoiseModelling (Universite-Gustave-Eiffel) | Java | GPL-3.0 | 2026-10-02 | CNOSSOS-EU / ISO 9613 urban mapping |
| I-Simpa (Universite-Gustave-Eiffel) | C++ | GPL-3.0 | 2026-01-18 | Ray/particle tracing, room and urban |
| openPSTD | C++ (GUI); earlier Python/Blender | GPL-3.0 | 2019-02-12 (inactive) | Fourier PSTD; the 2017 urban benchmark reference |
| Pachyderm Acoustic (PachydermAcoustic_Universal) | C# | GPL v3+ (General_Public_License.txt, © 2008–2017 A. van der Harten); Pachyderm_Noise_Control also GPL-3.0 | 2026-10-03 | Rhino/Grasshopper hybrid GA and FDTD |
| GSound / pygsound (GAMMA-UMD) | C++/Python | **UNC non-commercial**: executable only, no copying or redistribution | 2023-01-17 | Not usable |
| NVIDIA OptiX (NVIDIA/optix-dev headers) | C++/CUDA | NVIDIA SDK license, proprietary, NVIDIA GPUs only. Sec. 2.1 bars "avionics, navigation … Critical Applications" without a separate agreement and bars use that would subject it to an Open Source License | 2025-12-12 | Avoid |
| pyroomacoustics (LCAV) | Python | MIT | 2026-07-17 | Room ISM/ray tracing; reference only |
| CGAL | C++ | Per-file GPL-3+ or LGPL-3+; commercial license from GeometryFactory | 2026-10-02 | GPL packages incompatible with linking into MIT distribution |
| Clipper2 (AngusJohnson) | C++ | BSL-1.0 | 2026-04-20 | Polygon boolean and offset |
| earcut.hpp (mapbox) | C++ | ISC | 2026-09-28 | Polygon triangulation (roof caps for extrusion) |
| Boost.Geometry | C++ | BSL-1.0 | 2026-08-17 | Polygons, R-trees |
| CDT (artem-ogre) | C++ | MPL-2.0 | 2026-09-21 | Constrained Delaunay |
| TTL (SINTEF) | C++ | AGPL-3.0 | 2025-11-25 | Avoid |
| 3DBAG roofer | C++ | GPL-3.0 | 2026-07-30 | LoD2 reconstruction from lidar; executable only |

Sources:
- GitHub pages for each repo:
  - [embree](https://github.com/RenderKit/embree)
  - [steam-audio](https://github.com/ValveSoftware/steam-audio)
  - [EDtoolbox](https://github.com/upsvensson/Edge-diffraction-Matlab-toolbox)
  - [NoiseModelling](https://github.com/Ifsttar/NoiseModelling)
  - [I-Simpa](https://i-simpa.univ-gustave-eiffel.fr/presentation-1/licenses)
  - [openPSTD Mendeley](https://data.mendeley.com/datasets/kr876wdm92/1)
  - [GSound](https://github.com/GAMMA-UMD/pygsound)
  - [CGAL](https://github.com/CGAL/cgal)
  - [Clipper2](https://github.com/AngusJohnson/Clipper2)
  - [earcut.hpp](https://github.com/mapbox/earcut.hpp)
  - [CDT](https://github.com/artem-ogre/CDT)
- NoiseModelling GPL v3 is also stated in its docs — [NoiseModelling license](https://noisemodelling.readthedocs.io/en/latest/License.html)
- A separate pseudospectral TD code, FRidh/pstd, is BSD-3-clause per search — [github.com/FRidh/pstd](https://github.com/FRidh/pstd)

Other:
- **RAVEN (RWTH)**: C++, hybrid image-source plus ray tracing with diffraction and scattering, MATLAB interface. I found no public repository or open-source release statement; RWTH GitLab blocked automated queries — [FA2025 paper](https://dael.euracoustics.org/confs/fa2025/data/articles/000143.pdf)
- **IHTA urban auralization tool (Palenda & Vorländer 2026)**: UTD-based, higher-order diffraction and Doppler, described as open source. The released package is Windows wheels (pigeon_exe, pelican_exe, dovelet, pynamic) inside a CC BY-NC-SA 4.0 Zenodo archive, also on the IHTA Package Registry. A separate permissive code license was not found — [Zenodo 19112715](https://zenodo.org/records/19112715)
- **Mendeley openPSTD**: "open source implementation of the Fourier PSTD method … geared towards applications in the built environment"; EU project OPENPSTD — [Mendeley Data](https://data.mendeley.com/datasets/kr876wdm92/1); [CORDIS](https://cordis.europa.eu/article/id/174940-opensource-software-for-sound-propagation)

### Inferences
- A clean permissive stack is feasible:
  - Embree for ray-geometry queries
  - footprint extrusion via Clipper2 (cleanup) plus earcut.hpp (roof caps)
  - an in-house UTD or BTM diffraction kernel, ported from the formulas (EDtoolbox's BSD license also permits porting its code with attribution)
  - Steam Audio only if its pathing helps, which is unlikely for outdoor energy-band use
- For cross-checks, run NoiseModelling (CNOSSOS ray paths), I-Simpa and openPSTD as external executables. That matches the "GPL only as external executable" rule.
- The EDtoolbox license is BSD with two clauses, not BSD-3 as the task assumed. That is more permissive, not less.

### Gaps
- I found no permissive, maintained C/C++ UTD library for outdoor acoustics. UTD implementations seen are in GPL tools, the IHTA wheels (license unclear) or MATLAB (EDtoolbox is BTM, not UTD).
- RAVEN availability, and the Palenda tool's separate code license (if any, e.g. on the IHTA Package Registry), remain unconfirmed.
- The GSound commercial-license terms are unknown (contact UNC).
