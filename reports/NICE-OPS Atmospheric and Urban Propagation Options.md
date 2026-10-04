# NICE-OPS: Atmospheric and Urban Propagation Options

Sep 26, 2026 · @E G

## Summary

Use tabulated curved-ray eigenrays for the atmosphere, per-receiver reciprocity tables for cities, and PE/FFP only as a reference. Gaussian beams are not the right core.

- **Atmosphere:** eigenray tables per source height, built on the full moving-medium Snell invariants, feeding the existing ground model.
- **Cities:** per-receiver urban transfer tables indexed by arrival direction, built by tracing outward from each receiver (reciprocity).
- **PE/FFP:** an in-house GFPE and FFP to validate shadow edges and inversions, and optionally to build offline tables at low elevation angles.

The deciding factor is the ground. NICE-OPS already has a spherical-wave, finite-impedance ground model with boundary loss and turbulence coherence. No published Gaussian beam code keeps any of that: all use a plane-wave reflection coefficient or rigid ground.

The refraction study behind AAM found that for elevated sources the main effect is the changed ground-incidence geometry, not focusing ([Plotkin et al. 2013](https://rosap.ntl.bts.gov/view/dot/12035)). Curved-ray geometry fed into the existing ground model keeps exactly that.

Start with layered absorption and a wind-corrected source state, then build the eigenray tables. The last section lays out the steps.

## What shaped the choice

Four properties of NICE-OPS decide which methods fit: one propagation seam, a directivity lookup that refraction moves, shadow zones inside community-scale grids, and an atmosphere that can be tabulated.

### One seam, five quantities

All propagation passes through `spectral_level_at_observer` (`src/sphere_database.cpp`) and `reception_factor` (`src/util/sel.hpp`). Per source–receiver pair they use:

1. the emission direction, which `lookup_angles` turns into sphere angles;
2. spreading, `20 log10(r / r_sphere)`;
3. per-band absorption, `alpha_b * (r - r_sphere)`;
4. ground geometry, source height and ground distance, for the ground model, `GroundTable` or the plate receiver;
5. travel time, as range over c, in the reception-time weighting.

Any method has to supply these five per eigenray, about 94 million times per run on the demo case.

### Refraction moves the directivity lookup

For a linear effective sound-speed profile, rays are circular arcs of radius R = c0 / |dc/dz|. The launch angle then departs from the straight-line angle by about:

```latex
\Delta\theta \approx \frac{d}{2R} = \frac{d\,|dc/dz|}{2c_0}
```

Exact circular-arc launch-angle shift, source at 150 m, receiver at 1.2 m, downward refraction:

| dc/dz (1/s) | d = 0.5 km | d = 1 km | d = 2 km | d = 4 km |
| --- | --- | --- | --- | --- |
| 0.005 | 0.2° | 0.4° | 0.8° | 1.7° |
| 0.01 | 0.4° | 0.8° | 1.7° | 3.3° |
| 0.02 | 0.8° | 1.7° | 3.3° | 6.6° |
| 0.05 | 2.1° | 4.1° | 8.2° | 16.1° |

The S-76D database grid is 5° in theta and 10° in phi. A strong gradient at long range therefore moves the lookup by a full cell or more, and it lands on the shallow emission angles where rotor directivity changes fastest.

### Shadow zones fall inside community-scale grids

Under upward refraction the shadow begins at a horizontal range of:

```latex
x_s = \sqrt{2R}\left(\sqrt{h_s} + \sqrt{h_r}\right)
```

| dc/dz magnitude (1/s) | R (km) | h\_s = 15 m | h\_s = 50 m | h\_s = 150 m | h\_s = 300 m |
| --- | --- | --- | --- | --- | --- |
| 0.005 | 68.6 | 1.84 km | 3.02 km | 4.94 km | 6.82 km |
| 0.01 | 34.3 | 1.30 km | 2.14 km | 3.49 km | 4.82 km |
| 0.02 | 17.1 | 0.92 km | 1.51 km | 2.47 km | 3.41 km |
| 0.05 | 6.9 | 0.58 km | 0.96 km | 1.56 km | 2.16 km |
| 0.1 | 3.4 | 0.41 km | 0.68 km | 1.11 km | 1.53 km |

The demo grid spans about ±450 m, so it rarely reaches a shadow. Low approach segments and multi-kilometer community grids do.

### The atmosphere can be tabulated

Given a horizontally stratified atmosphere, flat ground and one receiver height, propagation depends only on source height, range and bearing. That is the structure `GroundTable` already exploits for (h\_s, d).

Reciprocity makes it cheaper still. By the flow-reversal theorem ([Godin 1997](<https://doi.org/10.1016/S0165-2125(96)00037-6>)), one fan traced upward from the receiver height, with the wind reversed, fills in every source height at once.

The tables and shadow ranges above are my own calculations from the linear-profile formulas. Real near-ground profiles are log-linear, with steeper gradients low down.

## The options

Six options, from cheapest to most general. Options 1, 2 and 5 are the recommended path; 3 and 4 are alternatives weighed and set aside; 6 is for later.

### 1. Layered absorption and wind-corrected source state

The cheapest meteorological improvement: integrate absorption through the profile, and use airspeed rather than ground speed for the source state.

- **Absorption.** Integrate the ISO 9613-1 coefficient through T, RH and p profiles. Precompute a cumulative F\_b(z) per band; path absorption is then S \[F\_b(z\_s) - F\_b(z\_r)\] / |z\_s - z\_r|. That costs O(1) per band, is exact for straight rays and first-order for curved ones.
- **Size of the effect.** Absorption along the profile, against homogeneous ground conditions, shifted LAE by -0.6 to +1.3 dB(A) at 500 m ([Zellmann & Wunderli 2014](https://www.acoustics.asn.au/conference_proceedings/INTERNOISE2014/papers/p557.pdf)).
- **Existing practice.** [14 CFR 36 Appendix A](https://www.law.cornell.edu/cfr/text/14/appendix-A_to_part_36) requires layering when alpha at 3150 Hz departs from its 10 m value by more than 0.5 dB per 100 m, in equal layers of at least 30 m.
- **Source state.** Wind changes airspeed, and so mu and C\_T. DLR found EC135 contours below the track dominated by this source-side change, and those farther out by propagation (Olsman & Lummer, ERF 2012, paper 046).

### 2. Stratified eigenray tables (recommended core)

Precompute ray tables per source height, look up the eigenrays per pair, and pass their geometry to the existing ground model.

- **Tables.** One per source-height node, every 10 to 25 m, on a polar (range x bearing) grid, stored per ray branch. Each entry holds the launch wave normal, grazing angle, travel time, ray-tube spreading, a caustic count and a validity mask.
- **Height spacing.** [Cotté 2019](https://doi.org/10.1121/1.5093307) found 10 source heights give 0.3 dB(A) worst-case error and 3 heights give 2.4 dB(A).
- **Full equations cost nothing extra.** In a stratified medium the horizontal slowness (s\_x, s\_y) is invariant (Pierce, *Acoustics*, ch. 8). The effective-sound-speed shortcut overstates the wind term by 1/cos(theta): 15% at 30° and double at 60°. Steep angles are exactly what an elevated source produces.
- **Ground.** Pass the reflected ray's grazing angle and the travel-time path difference to the existing ground model, as Nord2000 and [L'Espérance et al. 1992](<https://doi.org/10.1016/0003-682X(92)90022-K>) do. Nord2000 is validated this way to 1 km and acceptable to 3 km ([Ahearn et al. 2012](https://rosap.ntl.bts.gov/view/dot/9687)). Caveat: [Nguyen et al. 2022](https://arxiv.org/abs/2205.13695), citing Ostashev & Wilson, argue the spherical-wave coefficient's assumptions weaken under refraction.
- **Shadow zones.** Use the first residue mode plus a turbulence floor. Deep in the shadow the level is set by turbulent scattering, typically 20 to 30 dB below free field ([Gilbert, Raspet & Di 1990](https://doi.org/10.1121/1.399088)). [Forssén et al. 2014](https://doi.org/10.3813/AAA.918760) give an engineering single-scattering model.
- **Caustics.** Cap them with an Airy-type correction. Left uncapped they produce discontinuities above 10 dB ([Salomons 1998](https://doi.org/10.1121/1.423966)).
- **Cost.** Tables build in well under a second; a lookup should cost 50 to 100 ns per pair (my estimate). Searching eigenrays pair by pair costs about 5 µs each even with analytic layers ([Tuinstra 2014](https://doi.org/10.1260/1475-472X.13.5-6.337)), so it serves as an oracle, not in production.
- **Precedent.** [Plotkin et al. 2013](https://rosap.ntl.bts.gov/view/dot/12035) note that a stratified atmosphere over flat ground allows one-time precomputation by source elevation, distance and bearing.
- **Test oracles.** For c = c0 + az, rays are circular arcs with a closed-form eigenray and travel time. [infraGA](https://github.com/LANL-Seismoacoustics/infraGA) ([Blom & Waxler 2017](https://doi.org/10.1121/1.4980096)) traces the full 3D moving-medium equations and can run externally.

### 3. Gaussian beam tracing (considered, not recommended as the core)

Gaussian beams handle caustics and multipath well, but they are weakest on exactly the terms that dominate a helicopter footprint.

- **Strengths.** Caustics stay finite, no eigenray search is needed, and multipath in ducts and under downward refraction comes out naturally.
- **Ground.** BELLHOP applies a plane-wave reflection coefficient at the central ray, with no beam shift. Every atmospheric beam code found uses plane-wave reflection or rigid ground.
- **Shadow zones.** The shadow field depends explicitly on the beam parameter ([White et al. 1987](https://academic.oup.com/gji/article/89/2/579/632948)). Errors reach 5 dB in shadow and interference zones ([Bao et al. 2025](https://doi.org/10.1098/rspa.2024.0485)).
- **No width rule.** There is still no rule for the initial beam width. [Porter 2019](https://doi.org/10.1121/1.5125262) says this is why BELLHOP defaults to geometric "hat" beams, which are ray tubes, in effect option 2.
- **Urban.** Summing beams on the terrain surface gave errors up to 40 dB, below 8 dB only after an empirically fitted offset ([Yunus et al. 2023](https://doi.org/10.1121/10.0017825)). There was no edge diffraction; shadows were filled by beam width, with error growing with frequency.
- **Cost.** Coherent summation needs on the order of 10^7 beams per source position per frequency at 3 km and 2.5 kHz (my estimate from BELLHOP's own rule). bellhopcuda's documentation calls receivers spread across the domain its slow case.
- **Licensing.** BELLHOP and bellhopcuda are GPL-3. The only permissive Gaussian beam code found is [AcousticRayTracers.jl](https://github.com/org-arl/AcousticRayTracers.jl), a 2½-D Julia tracer.
- **Where it would earn a place.** A fully 3D inhomogeneous atmosphere, such as terrain-driven flow or LES fields, where per-height tables fail. Even there it would need the existing ground model grafted on.

### 4. Parabolic equation and FFP (reference, possibly offline tables)

PE is accurate in the grazing band and invalid on the steep near-track paths that dominate SEL, so it belongs in validation.

- **Angle limits.** Wide-angle and terrain PEs hold to about 30 to 35° elevation, and to about 60° with analytic starters (Salomons, *Computational Atmospheric Acoustics*, 2001). For a 100 m source, a 35° PE is valid only beyond about 143 m of horizontal range.
- **No launch angles.** PE returns a field, not rays, so the directivity lookup has nothing to key on.
- **Precedent.** Penn State/Volpe's HPM stitched GTPE below 35°, FFP from 35 to 48° and straight rays above 48° ([Boeker & Rosenbaum 2012](https://rosap.ntl.bts.gov/view/dot/9639); [Rosenbaum 2011](https://etda.libraries.psu.edu/catalog/11494)). It took more than 2 days per case, against under a minute for rays.
- **Best use.** An in-house GFPE ([Gilbert & Di 1993](https://doi.org/10.1121/1.407454)) with surface-pole subtraction ([Gilbert 2015](https://doi.org/10.1121/1.4902423)) and an exact starter ([Gilbert & Di 2007](https://doi.org/10.1121/1.2719044)). No open implementation was found, but the algorithm is compact. It would validate shadow edges, inversions and the turbulence floor.
- **Optional tables.** Reciprocal, band-averaged excess attenuation below about 25 to 30° elevation, blended into option 2, with directivity looked up at the ray launch angle. Estimated cost is 2 to 5 core-hours per meteorological profile.
- **Related work.** Coyle's 2012 Penn State M.S. thesis applied GFPE with weather and terrain (seen only as a listing).

### 5. Urban: per-receiver transfer tables by reciprocity (recommended for overflight)

For a source well above the rooftops, precompute each receiver's response to sound arriving from each direction, and look it up by the eigenray's arrival direction.

- **Build.** From each receiver, trace outward over the arrival hemisphere: image sources to order 2, Pierce or UTD double-edge roof diffraction, and an incoherent Lambert term. Include the street-surface reflection so the ground model does not count it twice.
- **Hot-loop cost.** One lookup per band, keyed by the arrival direction option 2 already supplies.
- **Size.** About 1 GB dense (14,400 receivers x 648 directions x 26 bands, float32). A parametric form is about 100 times smaller: a few paths per direction, each a path difference and a reflection count, with band values computed on the fly.
- **Closest precedents.** [Schäfer, Reich & Vorländer 2022](https://doi.org/10.1051/aacus/2022021) coupled a curved eigenray to an image-edge model through a virtual source. Launch-angle error was under 0.1° and spreading error about 1 dB, at 0.8 to 1.2 s per pair in MATLAB. [Gao, Porcayo & Clarke 2023](https://doi.org/10.1016/j.trd.2023.103794) used receiver-side reciprocity for UAM.
- **Production precedent.** [sonAIR](https://doi.org/10.1051/aacus/2025035) precomputes building attenuation on an airspace grid, at about 2 weeks per airport. Aircraft more than 15° above the obstacle horizon get a simple on-the-fly model. No published per-receiver, direction-indexed table for aircraft noise was found.
- **Magnitudes.** Rear facades near Schiphol were shielded by about 11 dB(A) for landings at 50 to 100 m and about 14 dB(A) for climb-outs ([Lugten et al.](https://api.repository.cam.ac.uk/server/api/core/bitstreams/044505b4-6042-47fe-a02f-a6ac4a5f70b9/content)). Reflections added 0 to 5 dB in illuminated zones for a hovering eVTOL ([Yunus et al. 2023](https://doi.org/10.1016/j.ast.2023.108410)).
- **Largest modeling uncertainty.** Specular versus diffuse facades can flip the sign of the canyon effect on EPNL ([Ismail & Oldham 2002](https://doi.org/10.1260/135101002320815693)).
- **Where it breaks.** Sources within about 10 times the neighborhood size: drones in canyons, vertiports, low passes. The reflection-path spreading error is about 20 log10(1 + dL/R), 3 dB at R = 100 m for dL = 40 m.
- **Validation.** EDtoolbox's BTM edge diffraction (BSD-3), openPSTD run externally (GPL), and the RWTH EDGE scale-model database.

### 6. Full 3D per-pair tracing (later, for near-field cases)

The general answer when the source is close to or inside the city: trace each source–receiver pair through the geometry.

- **Tools.** [Intel Embree](https://github.com/RenderKit/embree) (Apache-2.0) for BVH ray queries on CPU; OptiX on NVIDIA GPUs.
- **Method.** Schäfer-style virtual sources per pair, or full image-edge path finding with UTD.
- **Scope.** Vertiport and drone scenarios, not whole footprints: it costs microseconds to milliseconds per pair.

## Comparison

Only the eigenray tables deliver exact launch angles and keep the existing ground model at a hot-loop cost comparable to today's.

| Option | Refraction and ducting | Shadow zones | Impedance ground | Urban | Launch angle for directivity | Hot-loop cost per pair |
| --- | --- | --- | --- | --- | --- | --- |
| 1. Layered absorption | No (absorption only) | No | Unchanged | No | Straight line | Negligible |
| 2. Eigenray tables | Yes, per branch | Heuristic floor | Keeps the existing model | Through option 5 | Exact | \~50–100 ns (estimate) |
| 3. Gaussian beams | Yes, caustics smooth | Beam-width artifact | Plane-wave only | Specular, no diffraction | Smeared across beams | High |
| 4. PE / FFP | Yes, best at grazing | Yes, with turbulence | Yes | No | None | Offline only |
| 5. Urban tables | Uses option 2 | Diffraction included | Street surface included | Yes, source far away | Uses option 2 | One lookup per band |
| 6. Per-pair 3D | Yes | Depends on method | Must be grafted on | Yes, including near field | Exact | µs to ms |

## Risks

The biggest risk is not the method but the input: the near-ground gradient controls the grazing geometry, and it is rarely known well.

- **Input uncertainty dominates.** Equivalent-linear profiles deviated from PE by up to 10 dB at 500 to 2000 Hz ([van Maercke & Salomons 2009](https://pub.dega-akustik.de/NAG_DAGA_2009/data/articles/000250.pdf)). Scenario classes or ensembles are easier to defend than one deterministic profile.
- **The culling bounds stop being bounds.** `DirectionBounds` and the spectral bounds assume straight-line spreading and the straight-line direction. Downward refraction and caustics add gain, so each table cell needs a maximum propagation gain and each urban receiver a maximum urban gain. The SEL window itself is safe: since `9d6e986` it runs first-to-last, so a shadow dip cannot truncate it.
- **Sphere angle convention.** In wind, the wave normal and the ray direction differ by up to 1.7° at 10 m/s. NICE-OPS looks up the emission-time straight-line direction, which equals the wave normal only in still air. panam's hemigen function confirms the spheres use that same still-air, emission-time geometry, with azimuth taken from the ground-track heading. Under refraction the lookup key is therefore the eigenray's launch wave normal.
- **Table discontinuities.** Shadow boundaries, caustics and changes in branch count make naive interpolation fail. Each branch needs a validity mask.
- **Low rotor harmonics.** At 10 to 30 Hz, wavelengths are 11 to 34 m, at the edge of ray validity near the ground and around buildings.
- **Facade scattering.** The specular-versus-diffuse assumption can flip the sign of the urban effect, and no rotorcraft data settle it.
- **Stratification limits.** The stratified assumption fails over terrain, near buildings and across fronts.

## Licensing

The permissive tools cover what production needs; every GPL code is usable only as an external benchmark executable, never linked.

| Code or data | License | Role here |
| --- | --- | --- |
| [infraGA](https://github.com/LANL-Seismoacoustics/infraGA) (LANL) | MIT-style; modified versions must be marked | Ray-tracing oracle; could be adapted |
| [Intel Embree](https://github.com/RenderKit/embree) | Apache-2.0 | BVH ray queries for urban tracing |
| EDtoolbox (Svensson) | BSD-3 | Edge-diffraction (BTM) reference |
| [AcousticRayTracers.jl](https://github.com/org-arl/AcousticRayTracers.jl) | MIT | Only permissive Gaussian beam code (2½-D Julia) |
| Steam Audio (Valve) | Apache-2.0 | Reflection and occlusion; no refraction |
| pyram; PEmodel; simplePE.jl | BSD-3; MIT; MIT | PE starting points |
| NCPAprop | Permissive code, but needs FFTW and GSL (GPL) | External executable only |
| BELLHOP / Acoustics Toolbox; bellhopcuda | GPL-3 | External benchmark only |
| NoiseModeling; I-Simpa; openPSTD | GPL-3 | External benchmark only |
| WAPE (Univ. Gustave Eiffel) | EUPL-1.2 | Reference only |
| cTraceo | CC BY-NC-SA 3.0 | Unusable |
| OSM buildings | ODbL | Footprints; heights patchy |
| Microsoft GlobalMLBuildingFootprints | CDLA-Permissive-2.0 | Footprints with ML heights |
| USGS 3DEP lidar | US public domain | Building and terrain heights |

No open GFPE implementation was found. The algorithm is compact enough to write in-house.

## Where to start

Start with a refactor that changes no output: move the five propagation quantities behind one per-eigenray record. Then add the profile input and layered absorption. Make a closed-form linear-profile eigenray the first refracting model.

That order follows what each step unlocks. The refactor makes every later model a drop-in; the analytic model is both an early result and the oracle the general tables are tested against.

### Settle the conventions first

- **Sphere angles: confirmed.** panam's `hemigen` (`flight_acoustics.py`) builds the spheres from the straight line between the emission-time position and the microphone, in still air, with azimuth measured from the ground-track heading. Under refraction the lookup key becomes the launch wave normal; the heading stays as it is.
- **Advance ratio: open.** NICE-OPS takes mu from the track's ground-frame speed (`derive_flight_parameters` in `src/trajectory.cpp`). panam's `read_vehicle_data` takes the condition speeds from the test log's indicated-airspeed column, though a later variable calls them ground speed. Wind separates the two, so decide which the database means before wind enters.
- **Depropagation must match propagation.** panam depropagates along straight rays through a homogeneous atmosphere, so any test-day refraction is baked into the spheres. Once NICE-OPS propagates through profiles, panam should depropagate with the same model and the test-day profiles.

### Steps

1. **Refactor the seam, output unchanged.**
   - One per-eigenray record: launch direction, spreading, per-band path absorption, travel time and ground geometry. The straight-line model fills it exactly as today.
   - Choose the model outside the observer loop, as a template parameter, so the inner loop pays no dispatch. The windowed path already works this way with its level accessor.
   - Express ground geometry as an equivalent straight-ray (h\_s, d) that reproduces the grazing angle and path difference. `GroundTable`, the exact ground model and the plate receiver then work unchanged.
   - Done when the goldens pass unchanged at 1e-9 and the acoustic loop is no slower.
2. **Profile input and layered absorption.**
   - `--atmosphere <file>` with self-describing columns z, T, RH, p, u, v in the NCPAprop style, lengths in `--units`, scalars such as z0 in the header. Interpolate with a C² spline; `SmoothingSpline` already is one.
   - Precompute F\_b(z) per band once per run. Path absorption becomes S \[F\_b(z\_s) - F\_b(z\_r)\] / |z\_s - z\_r|, falling back to alpha S when the heights are equal.
   - Done when a uniform profile reproduces the goldens, F\_b matches numerical quadrature, and a layered case matches the 14 CFR 36 procedure.
3. **Wind-corrected source state,** if the advance-ratio question settles on airspeed. Subtract the profile wind at the source altitude from the track velocity before `derive_flight_parameters`.
4. **Closed-form linear-profile eigenrays.**
   - For c = c0 + az, the direct ray, the n-bounce reflected rays, launch angle, travel time and spreading all have closed forms. Fit the equivalent linear gradient per bearing from the profile, as Harmonoise and Nord2000 do.
   - It is exact for its profile, cheap enough to run per pair, and exercises the whole pipeline: directivity at the launch angle, curved-ray ground geometry, reception time and bounds.
   - Done when it reduces to the straight-line result as the gradient goes to zero and reproduces the launch-angle and shadow-onset tables above.
5. **Propagation-gain bounds.** Extend `DirectionBounds` and the spectral bounds with a maximum propagation gain. Check the windowed run against a full march, as ctest already does for the plate receiver.
6. **General stratified eigenray tables.**
   - Reciprocal fans traced up from the receiver height with the wind reversed, using the full moving-medium invariants. Per-branch validity masks, an Airy cap at caustics, and a residue-mode-plus-turbulence floor in shadows.
   - Done when they match step 4 on linear profiles to table resolution, and infraGA, run externally, on log-linear profiles with wind.
7. **Reference, then urban.** An in-house GFPE to calibrate the shadow floor, then the urban reciprocity tables.

The first deliverable is steps 1 and 2 together. The refactor lands with the goldens untouched, and layered absorption is the first measurable change, worth up to about 1 dB(A) at 500 m.

## Sources

The key sources behind this survey, grouped by topic. **Read** means the full text or code was read; **Abstract** means only the abstract, metadata or a secondary account.

### Gaussian beams

- [Porter 2019](https://doi.org/10.1121/1.5125262), JASA 146:2016. BELLHOP retrospective: no beam-width rule, why geometric beams are the default. **Read**
- [Porter & Bucker 1987](https://doi.org/10.1121/1.395269), JASA 82:1349. Gaussian beams for ocean acoustics. **Abstract**
- [Červený, Popov & Pšenčík 1982](https://doi.org/10.1111/j.1365-246X.1982.tb06394.x), GJRAS 70:109. Gaussian beam summation. **Abstract**
- [Gabillet et al. 1993](https://doi.org/10.1121/1.405722), JASA 93:3105. First atmospheric Gaussian beams; 2D, off by up to 10 dB in strong wind refraction (per Yunus). **Abstract**
- [White et al. 1987](https://academic.oup.com/gji/article/89/2/579/632948), GJRAS 89:579. Shadow-zone fields depend on the beam parameter. **Abstract**
- [Bao et al. 2025](https://doi.org/10.1098/rspa.2024.0485), Proc. R. Soc. A 481. Moving-media dynamic ray tracing; up to 5 dB error in shadows; rigid ground only. **Read**
- [Yunus et al. 2023](https://doi.org/10.1121/10.0017825), JASA 153:2362, and [Yunus 2023 thesis](https://repository.tudelft.nl/record/uuid:72d10b7a-6790-41fc-9b15-26f9cccdb77f). 3D beams with wind and buildings; 40 dB surface artifact; no diffraction. **Read**
- [Nguyen et al. 2022](https://arxiv.org/abs/2205.13695). BELLHOP adapted to wind-farm propagation; plane-wave ground reflection by choice. **Read**
- [Mo et al. 2017](https://doi.org/10.1121/1.4977005), JASA 141:2289. Ray-curve tracer plus Gaussian beams on the Attenborough benchmark. **Read**
- [Zhang, Duan & Jiang 2025](https://doi.org/10.1016/j.cpc.2025.109722), CPC 315. GPU Gaussian beams in a city: 88x overall speedup. **Read**
- [bellhopcuda](https://github.com/A-New-BellHope/bellhopcuda). GPU port; GPL-3; spread receivers are its slow case. **Read**

### Rays in moving, stratified atmospheres

- [Pierce, *Acoustics*](https://doi.org/10.1007/978-3-030-11214-1), ch. 8–9. Moving-media ray equations, Blokhintzev invariant, caustics, shadow zones. **Read**
- [Ostashev & Wilson 2015](https://doi.org/10.1201/b18922), *Acoustics in Moving Inhomogeneous Media*. General theory. **Abstract**
- [Blom & Waxler 2017](https://doi.org/10.1121/1.4980096), JASA 141:2681, and [infraGA](https://github.com/LANL-Seismoacoustics/infraGA). Eigenray search and Jacobian amplitudes. **Read**
- [Schäfer & Vorländer 2021](https://doi.org/10.1051/aacus/2021018), Acta Acust. 5:26. ART stratified ray tracer, full equations with advection. **Read**
- [Tuinstra 2014](https://doi.org/10.1260/1475-472X.13.5-6.337), IJA 13:337. Analytic-layer eigenrays at about 5 µs each. **Read**
- [L'Espérance et al. 1992](<https://doi.org/10.1016/0003-682X(92)90022-K>), Appl. Acoust. 37:111. Heuristic curved-ray model with partial coherence. **Abstract**
- [Salomons 1998](https://doi.org/10.1121/1.423966), JASA 104:3259. Caustics in log profiles; >10 dB discontinuities without diffraction. **Abstract**
- [van Maercke & Salomons 2009](https://pub.dega-akustik.de/NAG_DAGA_2009/data/articles/000250.pdf), NAG/DAGA. Harmonoise curved ground against PE: up to 10 dB. **Read**
- [Plotkin et al. 2013](https://rosap.ntl.bts.gov/view/dot/12035), Wyle WR 13-01. AAM/RNMwea refraction; ground geometry dominates; tabulation. **Read**
- [Page et al. 2015](https://onlinepubs.trb.org/onlinepubs/acrp/docs/ACRP02-44_FR.pdf), ACRP 02-44. Refraction in airport noise models; AEDT gap. **Read**
- [Heath & McAninch 2006](https://ntrs.nasa.gov/citations/20060004780), AIAA 2006-411. NASA stratified ray code; critical emission angle. **Read**
- Olsman & Lummer 2012, ERF paper 046. EC135: wind's effect on source state versus on propagation. **Read**
- [Berton 2021](https://ntrs.nasa.gov/citations/20210011612), AIAA Aviation. Lateral attenuation fading by 30–40° elevation. **Read**
- [Ahearn et al. 2012](https://rosap.ntl.bts.gov/view/dot/9687), DOT/FAA/AEE/2012-5. Nord2000 for aircraft; validated to 1 km. **Read**
- [Zellmann & Wunderli 2014](https://www.acoustics.asn.au/conference_proceedings/INTERNOISE2014/papers/p557.pdf), Inter-Noise. Layered absorption: -0.6 to +1.3 dB(A). **Read**
- [14 CFR 36 Appendix A](https://www.law.cornell.edu/cfr/text/14/appendix-A_to_part_36). Layered absorption rule. **Read**

### Turbulence and shadow zones

- [Gilbert, Raspet & Di 1990](https://doi.org/10.1121/1.399088), JASA 87:2428. Turbulent scattering into shadow zones. **Abstract**
- [Forssén et al. 2014](https://doi.org/10.3813/AAA.918760), Acta Acust. 100:810. Engineering single-scattering model. **Read**
- [Daigle, Piercy & Embleton 1983](https://doi.org/10.1121/1.390152), JASA 74:1505. Coherence of direct and reflected paths. **Abstract**
- [Berry & Daigle 1988](https://doi.org/10.1121/1.396385), JASA 83:2047. Residue series against controlled experiments. **Abstract**

### Parabolic equation and FFP

- [Gilbert & Di 1993](https://doi.org/10.1121/1.407454), JASA 94:2343. The Green's-function PE. **Abstract**
- [Gilbert 2015](https://doi.org/10.1121/1.4902423), JASA-EL 137. Surface-pole subtraction. **Abstract**
- [Gilbert & Di 2007](https://doi.org/10.1121/1.2719044), JASA-EL 121. Exact starter for any source height. **Abstract**
- [Salomons 2001](https://doi.org/10.1007/978-94-010-0660-6), *Computational Atmospheric Acoustics*. Rays, PE, FFP. **Abstract**
- [Attenborough et al. 1995](https://doi.org/10.1121/1.412302), JASA 97:173. PE/FFP/ray benchmark cases. **Abstract**
- [Cotté 2019](https://doi.org/10.1121/1.5093307), JASA 145:1363. PE tabulated by source height: 10 heights, 0.3 dB(A). **Read**
- [Colas et al. 2023](https://doi.org/10.1121/10.0020834), JASA 154:1413. WAPE against LEE; cost figures. **Read**
- [Nyborg et al. 2023](https://doi.org/10.1121/10.0020562), JASA 154:1299. PE angle limits and grid requirements. **Read**
- [Boeker & Rosenbaum 2012](https://rosap.ntl.bts.gov/view/dot/9639), Inter-Noise, and [Rosenbaum 2011](https://etda.libraries.psu.edu/catalog/11494), Penn State PhD. HPM hybrid: over 2 days per case. **Read** / **Abstract**
- [Godin 1997](<https://doi.org/10.1016/S0165-2125(96)00037-6>), Wave Motion 25:143. Flow-reversal reciprocity. **Abstract**
- [NCPAprop](https://github.com/chetzer-ncpa/ncpaprop-release). Normal modes and PE; license and dependencies. **Read**
- Coyle 2012, Penn State M.S. thesis. GFPE with weather and terrain. Listing only.

### Urban propagation

- [Schäfer, Reich & Vorländer 2022](https://doi.org/10.1051/aacus/2022021), Acta Acust. Curved eigenray coupled to an image-edge model. **Read**
- [Schalcher et al. 2025](https://doi.org/10.1051/aacus/2025035), Acta Acust. sonAIR/sonX buildings and meteorology; cost. **Read**
- [Gao, Porcayo & Clarke 2023](https://doi.org/10.1016/j.trd.2023.103794), TRD. Receiver-side reciprocity for UAM. **Abstract**
- [Yunus et al. 2023](https://doi.org/10.1016/j.ast.2023.108410), AST. eVTOL at a vertiport: 0–5 dB from reflections. **Read**
- [Lugten, Steemers & Kang](https://api.repository.cam.ac.uk/server/api/core/bitstreams/044505b4-6042-47fe-a02f-a6ac4a5f70b9/content). Schiphol facade shielding, 11–14 dB(A). **Read**
- [Lugten et al. 2024](https://doi.org/10.1016/j.buildenv.2024.111400), Building and Environment. 2383 flyovers in a courtyard field lab. **Abstract**
- [Ismail & Oldham 2002](https://doi.org/10.1260/135101002320815693), Building Acoustics. Low-flying aircraft in canyons; diffuse reflection flips the sign. **Abstract**
- [Kang 2000](https://doi.org/10.1121/1.428580), JASA. Specular versus diffuse facades. **Abstract**
- [Riegel & Sparrow 2024](https://doi.org/10.1121/10.0024980), JASA. Sonic booms in canyons; diffusion dominates. **Abstract**
- [Dragna et al. 2022](https://doi.org/10.1121/10.0016442), JASA. Sonic booms in urban sections; \~8000 CPU-h per case. **Read**
- [Remillieux et al. 2011](https://doi.org/10.1121/1.3654486), JASA. Image sources plus BTM, validated against measurements. **Abstract**
- [Wei, Van Renterghem & Botteldooren 2015](https://doi.org/10.1121/1.4932585), JASA. Pierce multi-edge diffraction within 2 dB of full-wave. **Abstract**
- [Palenda & Vorländer 2026](https://doi.org/10.1051/aacus/2026033), Acta Acust. UTD higher-order artifacts. **Read**
- [Schissler, Mehra & Manocha 2014](https://doi.org/10.1145/2601097.2601216), ACM TOG. UTD against BTM: 3–4 dB at first and second order. **Read**
- [Van Renterghem & Botteldooren 2010](https://doi.org/10.1121/1.3425731), JASA. Limited refraction into shielded canyons. **Abstract**
- [Rizzi et al. 2020](https://ntrs.nasa.gov/citations/20205007433), NASA/TP-2020-5007433. UAM noise gaps, including urban propagation. **Read**
