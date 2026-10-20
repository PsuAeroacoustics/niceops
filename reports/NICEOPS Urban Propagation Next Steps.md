# NICEOPS: Urban Propagation Next Steps

Oct 19, 2026

## Summary

The urban model now runs a point source among LoD1 buildings inside a trajectory run (`--buildings`), with
reflections from facades, roofs and the floor, diffuse scattering by the facades and roofs, and diffraction at one or two edges by UTD or the edge-source
integral. A full-wave FDTD reference checks it, and EDGE scale-model data and a DigiCity reconstruction validate it.
Rays leaving the source apart add in energy beyond a 1° cone, unless the database records the vehicle's coherence.

What follows is ordered by how much each step would improve predictions for drones, helicopters and AAM among
buildings, judged against what EDGE, DigiCity and the FDTD have shown. Recommended order:

1. **Facade scattering in the run.** Done (Oct 19): `--facade_scattering` (1b).
2. **A third edge, and reflections between edges.** In progress (1a): three roof edges, a reflection between
   two edges and two at each end of a twice-diffracted ray are in; vertex partners and higher orders remain.

The rest follow as data and need decide.

**The research review** ([Urban propagation research review 2026-10-20](Urban%20propagation%20research%20review%202026-10-20.md))
ranks what the literature says would improve the predictions most, for rotorcraft and drones among buildings:

1. **Turbulence scattering into building shadows** (Forssén et al. 2014). Done (Oct 20): `--scattering` γ_T with
   `--buildings`, off by default; unvalidated, the strength being the weak input. It matters where a microphone
   stays shielded through an event.
2. **Finite facades and edge ends treated uniformly**: reflections near a facade's end weighted by their Fresnel
   zone, the vertex term, and a ray whose partner is beyond the orders faded rather than kept whole (1a's
   leftover, UE06's R1). Begun (Oct 20): the edge source carries the vertex for a ray that reflects between two
   edges (UE06 2.63 → 1.01 dB rms); by UTD it remains.
3. **Facade scattering by incidence and grating order**: s(θ, f) from the BEM, grating orders for periodic relief
   (UE24-26), and a receiver-side gather that makes the diffuse field affordable on a grid. Done (Oct 20) but the
   gather: `--relief`/`--facade_relief` takes a relief's grating orders from the BEM (`profile_scattering.py
   --orders`), its specular share at each ray's incidence, the other orders as grating rays, and by default what it
   scatters onto other walls diffusely at each path's incidence. Against EDGE's reliefs the orders match the FDTD
   (UE24-26 rms 1.24, 1.59, 1.63 against its 1.25, 1.60, 1.60; Lambert's s(f) 1.48, 1.85, 2.25); the default
   hybrid beats s(f) everywhere, DigiCity's microphone 7 included (−2.0 / 2.7 against −2.3 / 2.9).
4. **Pitched roofs (LoD2)**, validated on EDGE's measured pitched-roof scenes (UE01-05, 15, 17, 19). Done (Oct
   20): `--roofs` gives a CSV district's buildings roof planes (the least of them, capped at the top), with walls
   up to the roof, eaves, rakes, ridges and hips as wedges of their own angles, and rays over four roof edges or
   more traced as chains (`max_chain`). Against EDGE's pitched scenes the planes are 3.1 dB rms where the same
   houses flat are 4.5 at the mean height and 5.8 at the ridge (the mansard's floor source 14.3 → 4.5 with
   chains); against boundary elements a gable is 0.17 dB rms. Open: slopes reflect no specular ray; over a
   mansard's nearly flat wedges the cascade is 3.9 dB rms off the BEM, which a cascade of transitions (Tzaras and
   Saunders) or the edge source's higher orders would address; the FDTD's staircase is 5 dB loud on a slope at
   EDGE's 3 mm cells, so it is no reference there; GeoJSON districts take no roofs yet (OSM's `roof:shape`).
5. **Higher orders with pruning** (1a), and an engineering term for paths beyond the explicit orders.

For drones the review expects the source model to outweigh these until a measured multirotor sphere exists (4a).

## 1. Paths the ray model misses

### 1a. A third edge, and reflections between edges

**Evidence.**

- Behind a second row of buildings (EDGE UE13, E1 R3) the model is 8.5–13 dB low below 4 kHz. The path goes over
  three edges, and the model traces at most two.
- A twice-diffracted ray carries only one reflection on each side, so the reflected paths of a canyon are missing
  between edges. EDGE UE06's floor source is about 5 dB off for this reason.
- In the DigiCity continuity sweep (`urban/continuity.py` in the harness) there are 434 steps over 1 dB in 79,000
  positions. The reduced sweep has 30, all of them untraced higher-order partners, none at an edge's end.

**Approach.** Extend the path search in `src/point_paths` to three edges and to reflections between edges. The number
of paths grows quickly with order, so this needs pruning like `src/edge_visibility`'s. The partners such paths make at
boundaries need the same continuity treatment the second order has.

**Check.** UE13 against measurement and FDTD, UE06's floor source, and the continuity sweep's step count.

**Done so far (Oct 19-20).** Three changes to `src/point_paths`, all on by default:

- **Two reflections at each end of a twice-diffracted ray** (`max_twice_reflections`, was one): the street
  between two rows after the roof. UE13's receivers between the rows went from −3.6 dB mean, 4.5 rms to −2.0,
  3.2 against the FDTD.
- **One reflection between two edges** (`max_between_reflections`), in a facade or the floor: the second edge
  mirrored in it. UE06's floor source went from 5.3 to 2.8 dB rms against the FDTD; its R2-R4 are within 1.6.
- **Three roof edges** (`triple_diffraction`, `roofs`; `all` for any three): UE13's receiver behind the second
  row went from 8.2 to 2.2 dB rms. Every edge cost eight times the demo run for no level that moved.

Rays against the FDTD over EDGE's 13 runs: 2.80 → 1.95 dB rms. Against the measurements (1124 band levels, UE06
and UE07 read swapped): 4.01 → 3.82 dB rms by UTD. DigiCity microphone 7: −4.4 / 5.5 → −4.1 / 5.2 dB (specular), with the source calibrated on three other
microphones that the new rays reach too.
The continuity sweep (20 lines, 3771 positions): 31 → 23 steps over 1 dB, the largest 9.5 → 4.9 dB. The demo
`--buildings` run costs 80% more on the CPU (3.5 → 6.3 s). FIDELITY and THEORY 6.5 have the details.

**Still open:**

- **A reflection near a facade's end whose partner would be a vertex, by UTD.** UE06's R1 is 4.7 dB high by
  UTD: its ray reflects 3 cm inside a tall house's wall, where the field is about half, and the corner that
  would diffract the partner meets the low roof's front edge. With `--edge_source` the ray now takes the
  integral unfolded about the wall, whose end over the front edge is the vertex's wave: R1 −0.1 dB, the scene
  1.0 dB rms, EDGE's 13 runs 1.35 against UTD's 1.95. A Fresnel weight on the reflection would not do it (its
  zone is a few centimetres there: −0.4 dB at 8 kHz). Making it the default needs the edge source on the GPU
  and cheaper, or a vertex term in UTD (Albani et al. 2009). Cheaper since (Oct 20): looser tolerances, the
  cascade times the edges' end ratios out of both transitions, and a pair's integrals over threads take
  DigiCity's four microphones from 145 to 10.5 s and the demo district's 266 pairs from 497 to 19 s, within 0.012
  dB A-weighted; EDGE's 13 runs take 3.6 s. `niceops --buildings --edge_source` runs it in a trajectory. Still
  0.8 s of CPU a pair over the demo district (UTD 4 ms): what is left is near-transition rays at low frequency and
  shared faces, for which a uniform asymptotic form (Capolino and Albani's generalized Fresnel integral) is the
  step; and its largest departures from UTD, +20 dB deep in a building's shadow by a corner's waves, are
  unchecked.
- **Vertex diffraction** in general, more than one reflection between edges,
  reflections on a ray of three edges, and the shared face's aperture on a ray of three (it takes the cascade).
- **Four edges or more**, which a thin middle row needs (each of its roof edges counts).

### 1b. Facade scattering in the run

**Evidence.** `niceops_point_field --scattering_spectrum` already carries the diffuse field.

- **DigiCity microphone 7:** with the containers' corrugation (36 mm on a 278 mm pitch, coefficient from
  `tools/profile_scattering.py`), it goes from −4.1 dB mean / 5.2 dB rms to −2.3 / 2.9.
- **EDGE canyons with relief** (UE24–26 against UE23): measured −1.0 to −1.6 dB, which a specular-only model misses.

**Done (Oct 19).** `niceops --buildings` takes `--facade_scattering` (or each building's scattering column),
`--facade_scattering_spectrum`, `--patch_size` and `--diffuse_radius`. `UrbanDistrict::levels` builds
`DiffuseFacades` per microphone, gives the rays rho (1 - s) for the walls and roofs each field covers, and adds the
diffuse energy band by band, on the CPU and the GPU. Each lighting path reads the sphere along its own launch
(`DiffuseFacades::SourcePower`).

- **Checks.** For an omnidirectional source the run is `niceops_point_field`'s to 1e-9 dB (the GPU to 1e-6), so
  DigiCity microphone 7 with the corrugation (−2.3 dB mean, 2.9 rms) holds in a run. A directional source lights a wall as
  Lambert's law says to 0.5%. The demo district's states sampled agree with every state to 0.003 dB. A run with
  nothing in reach is the specular run byte for byte. FIDELITY has the demo district's SELs (s = 0.3: −0.22 to
  +1.03 dB).
- **Cost.** 0.2 s on the GPU and 0.4 s on the CPU over the demo district's seven microphones. Each microphone's
  field is O(patches²) to build, so a grid of thousands of microphones is out of reach this way.
- **Still open.** A spectrum per building: now one spectrum multiplies each building's coefficient.

## 2. Ground and receivers, matched to how measurements are made

### 2a. Ground impedance and plate microphones under `--buildings`

**The gap.** The floor among the buildings is one energy absorption coefficient. `--ground_effect`, `--plate_table`
and `--ground_classes` are refused.

**Why it matters.** Flight-test arrays use flush ground boards, and real streets mix pavement, grass and parks.

**Approach.** Give the floor reflection in the ray sum the spherical-wave reflection coefficient the open-ground
model already has, per band and per ray, from the reflected ray's grazing angle and path. Plate receivers then follow
as in the open-ground run. A floor of mixed classes needs each floor bounce's class, from `--ground_classes`.

**Check.** One far building with the impedance ground equals the open-ground run with the same ground (as
`buildings_far_is_open_without_ground` does for a rigid floor).

**Done (Oct 20).** `--buildings` takes `--ground_effect`, `--ground_classes` and `--plate_table`. Each ray records
its floor bounces (leg and grazing angle), and each bounce takes the open-ground model's Q at the sub-band's middle,
on the CPU and the GPU. A flush microphone on a plate takes P_d on the direct ray and Q P_r on its floor partner.
`niceops_point_field --flow_resistance` gives the same floor.

- **Checks.** Over open ground with Q held at the band's centre, the ray sum is the open-ground model to 1e-5 dB;
  Q's change across a band adds up to 0.12 dB. With one far building, the run is the open-ground run to 0.06 dB of
  SEL (4 ft microphones on grass), 0.14 dB (on concrete) and 0.003 dB (flush microphones on plates); the difference
  is the floor ray reading the sphere along its own launch. The GPU matches the CPU to 1e-4 dB.
- **Effect.** On the test district a grass floor takes 0.8 to 4.0 dB off the rigid floor's SELs, the flush
  microphones most; 5% more time.
- **Still open.** Each bounce takes the microphone's ground class, as the open-ground model does. A floor that
  changes along the street (a park beside the road) needs a map of the ground, which `--ground_classes` does not
  give. The ground model's turbulence is not taken among buildings.

### 2b. Facade and rooftop microphones, checked against data

**State.** The run takes microphones where the observer file puts them, at any height outside the buildings, roofs
included: a microphone 4 ft above a 45 ft roof, and one 0.66 ft off a wall at 20 ft, both run. Neither has been
compared with a measurement. In `niceops_point_field` the roofs scatter too, which raises a rooftop microphone by
1-2 dB over walls scattering alone (FIDELITY, "Scattering walls (point source)").

**Why it matters.** Urban exposure is assessed at facades, not at 1.2 m above the street.

**Data.**
- **Schiphol courtyards** (TU Delft and AMS): microphones 0.2 m off the facades.
- **MOTUS** (Bruitparif, Paris-Issy heliport): rooftop sites, and the closest data to the application.

## 3. Cost, the barrier to district-scale maps

**The numbers.** The demo approach (6,551 states × 7 microphones) takes 1.2 s on the GPU (3.2 s on 14 threads), against
0.3 s without buildings; an evaluated pair costs about 3 ms of CPU, so a district grid of 10,000 receivers is still far off.

**Done (Oct 19).** 70 s → 1.2 s on the demo approach (`docs/internals/ground.md`, FIDELITY):
- the UTD transition function from a table, and a shared face's node geometry worked out once for every wavenumber
  (exact; 70 → 37 s);
- the shared faces sampled: the integral only where the cascade cannot stand in for it, two thirds of an octave
  apart, the cascade times the interpolated ratio elsewhere (within 0.1 dB per ray; 37 → 10 s);
- the states sampled: the buildings' transfer at adaptively chosen states, the source exact at every state (SEL within
  0.003 dB, LAmax within 0.007; 10 → 3.3 s).
Measured and left out: windowing the coherence sum by the sinc's reach (4%).
- the GPU in the run, by default wherever there is one (the shared faces sampled there too, the cone's sum in the
  pair-sum kernel): 3.3 → 1.2 s. The measured coherence and `--phase_interval` stay on the CPU.

**Ray-tracing hardware does not help yet.** RT cores accelerate ray–triangle tests, which here are the path search's
visibility tests: 1–2% of a `--buildings` run before the sampling above, about 6% of an evaluated pair after it, and
12% of a pair on the 49-block district. It could pay later for scenes
of many triangles (LoD2, terrain, trees) or a stochastic diffuse field (`docs/internals/gpu.md`).

In order of payoff:

1. **The peak search, back.** It is off because a coherent sum could exceed the straight path's bound. Under the
   default (energy beyond a cone) a bound is much easier to state: the energy sum is bounded by the sum of the rays'
   bounds, and the cone's pairs add at most their own cross terms.
2. **Per-receiver reciprocity tables** (the roadmap's urban step). Trace each microphone's tree of paths once and
   reuse it for every source position.
3. **The plane-wave table where the source is far** from the buildings. Decided: worth it only where it saves a lot
   and the source is well clear. The switch needs a criterion, such as distance from the nearest building over that
   building's height, measured against the point model.

**Check.** Each must reproduce the CPU point model to the existing tolerances: the peak search exactly (as `*_windowed_is_exact` tests do), the tables to a stated dB.

## 4. Sources for drones and AAM

### 4a. A real multirotor sphere database

Every urban study so far used the B407's hover sphere as a stand-in. Two ways to get a drone database:
- PANAM's hexacopter recording;
- NLR's DJI M600 flights at DigiCity (8.5 and 20 m, logged position and rotor speeds), if NLR shares them.

Either would also give the first multirotor coherence record (`directional_coherence.fit_coherence`). The
hexacopter's tonal fraction stays 0.23–0.40 from 125 to 315 Hz, where the default law gives 0.05–0.20.

### 4b. Distributed rotors near buildings

**The problem.** A large AAM vehicle's rotors span 10–15 m. Hovering 10–20 m from a facade, a single point source
misplaces the reflections.

**Approach.** Place each rotor as its own source, mutually incoherent and each with its share of the sphere. Also
check the rays shorter than the sphere's radius, which get no spreading gain inside it, in the same close-range
cases.

## 5. Lower priority

- **Low-frequency wave effects.** Canyon resonances below about 150 Hz, where the rotor tones sit. The wave solver
  handles ideal canyons only (`--wave_below`). The FDTD could supply low-band corrections per scene.
- **Terrain under buildings, vegetation.** `--buildings` stands everything on the plane z = 0 and refuses
  `--terrain`. Pitched roofs are in (`--roofs`, item 4 above).
- **Atmosphere in streets.** Refraction is minor at urban ranges. The decorrelation of multipath by turbulence is
  partly stood in for by the cone.

## Validation still open

- **NLR Drone DigiCity:** the measured data and the drone flights, if NLR shares them. The user is contacting NLR;
  what to ask for is in the harness's `urban/README.md`. The microphone 7 comparison now rests on spectra digitized
  from the paper.
- **EDGE UE06 and UE07:** the two measurement files *may be* swapped against their geometry. The FDTD of each
  geometry matches the other's measured changes to about 1 dB. Other explanations: the models measured differed from
  those documented, or another experimental fault. A question for the authors, if wanted
  (pascal.palenda@akustik.rwth-aachen.de).
- **Schiphol courtyards and MOTUS:** data requests, for 2b. The cases, what each tests and whom to ask are in
  [Urban propagation validation cases 2026-10-18](Urban%20propagation%20validation%20cases%202026-10-18.md).

## Decisions that bound this work

- **Coherence.** Rays leaving the source apart add in energy beyond a 1° cone (`--coherence_angle`). The measured
  model, and the main rotor's delay with it, come only from a database's record, which PANAM's
  `directional_coherence` writes. Do not assume the source's phase.
- **The plane-wave table** is for sources well clear of the buildings, and only where it saves a lot.
- **Gaussian beams** were set aside as the core in the propagation roadmap.

## Where things are

| | Location |
| --- | --- |
| Point model | `src/point_paths`, `src/building_paths`, `src/edge_source`, `src/edge_visibility` |
| In a run | `src/urban_point`, `src/source_coherence` |
| Full-wave reference | `src/fdtd`, `tools/fdtd_field.cpp` |
| Diffuse facades | `src/diffuse_facades` |
| Developer notes | `docs/internals/ground.md` |
| Measurements | `docs/FIDELITY.md` |
| Harness: DigiCity, EDGE, continuity sweep | `urban/digicity/`, `urban/edge/`, `urban/continuity.py` |
| Harness: coherence | `coherence/` |
| Validation leads | [Urban propagation validation cases 2026-10-18](Urban%20propagation%20validation%20cases%202026-10-18.md) |
