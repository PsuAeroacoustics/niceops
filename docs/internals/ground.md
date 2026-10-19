# Receivers, ground effect, terrain and weighting

The user-facing description is in the manual's [Ground effect](../MANUAL.md#ground-effect),
[Receivers](../MANUAL.md#receivers),
[Ground-plane microphones](../MANUAL.md#ground-plane-microphones) and
[Terrain](../MANUAL.md#terrain). The equations are in [THEORY.md](../THEORY.md) §5 and §6. The
tables' interpolation accuracy is in
[FIDELITY.md](../FIDELITY.md#ground-and-plate-tables), and the measured effect of ground and
terrain under [Ground effect](../FIDELITY.md#ground-effect) and
[Terrain at Amedee](../FIDELITY.md#terrain-at-amedee). How the ground path was profiled, and the
optimizations measured and rejected there, are in [performance.md](performance.md).

## Receivers

Every level is evaluated at the microphone, the observer's point raised by its height above the
ground (`process_trajectory`). Which receiver a microphone gets is a user rule, given in the
manual. What to keep:

- **One ground model and one `GroundTable` per distinct (receiver height, ground class)**, built
  once per trajectory, not per observer. An observer file's fifth column lets one run mix plates
  and poles, and its sixth puts each microphone over a `--ground_classes` class. Keys are kept in
  the order first met, so with one class the tables, logs and the refraction setup's
  `point_heights` are what they were with heights alone.
- **The peak search's ground bound is the larger of the plate's and the point receiver's**, which
  does not depend on height, so a mixed run needs no per-observer bound. The point receiver's is
  `point_ground_gain_bounds` over every class point receivers stand on, the elementwise maximum
  of each class's `ground_gain_bounds`.
- **`mixed_receivers_are_their_subsets`** checks a mixed run (`tests/data/mixed_receivers.csv`)
  against runs over each kind. It runs with `--exact_ground`, because a table's extent follows its
  own run's observers, so tables built for different observer sets differ at the last digits.
- **Ground classes** (`src/ground_classes.{hpp,cpp}`, read as in
  [inputs.md](inputs.md#the-ground-classes-file)). Each class is the run's `GroundModel` with
  the file's overrides, and each has its admittance (`class_impedance`, one `ground_impedance()`
  call per class per trajectory). The hot loop picks the observer's admittance once per observer
  (`impedance_here`) and passes it where the run's was passed, on flat ground and over terrain, so
  an evaluation costs nothing more. What to keep:
  - **One class is the run without classes, byte for byte**: `class_models` is
    `{options.ground}`, `class_impedance[0]` is the one call made without classes, and
    `ground_class_single_is_scalar_{point,plate,terrain}` (and the explicit-columns twin) compare
    the files. A class file's numbers in those tests are whole, so the reader's and CLI11's
    parsing cannot differ by an ULP.
  - **`ground_impedance_bands` is the first class's**, and the refraction models and a table's
    plate offset take it. That is right only because `--refraction` and `--propagation_table`
    are refused with classes (an assert in the refraction setup says so). Per-class shadow
    transitions and bounce losses belong with ground along the path, not under the microphone.
  - **Flush microphones are on one class**, refused otherwise, and the plate is bound over it
    (`ground_model` points at it). The check is `PlateReceiver::ground_mismatch`, so `bound()`'s
    message on a run without classes is unchanged.
  - **Rejected:** a `DirectionBounds` table per class (more cache traffic for a bound that moves
    little, since each class's bound is (1 + max|Q|)² and |Q| nears 1 at grazing for every
    surface; revisit only if the measured evaluated fraction moves: half concrete moved it from
    98.48% to 98.57%, [performance.md](performance.md#terrain-ground-classes-and-screening));
    scaling the bound per observer (impossible, since the ground gain is folded into the
    A-weighted direction table); substituting a lone class for `options.ground` (it would enable
    refraction over a ground the user did not name); ids as row indexes (they move when rows are
    reordered, and do not match land-cover codes).
  - **A later blend over surfaces** (Fresnel-zone weighting along each path) stays under the same
    bound only if it blends power ratios convexly; otherwise the bound must be re-derived.
  - Tests: `ground_classes_mixed_are_their_subsets` (`--exact_ground`, as
    `mixed_receivers_are_their_subsets`), `ground_classes_tabled_matches_exact`,
    `ground_classes{,_terrain}_impedance_is_the_observers`,
    `ground_classes{,_terrain}_windowed_is_exact`, and the refusals.
- The direct path ends at the microphone, not on the ground. Moving it there changed
  point-receiver results; see [HISTORY](../HISTORY.md#changes-that-affect-results).

## Urban receivers (`--urban_table`)

`src/urban.{hpp,cpp}` holds the facades, the first-order image test, the ideal canyon under
plane-wave incidence, and `UrbanTable`; `tools/urban_table.cpp` builds tables; `urban_level` and
`urban_bound_geometry` in `process_trajectory` are the run's side. The design and its staged plan
are the Claude Doc "NICEOPS: Urban Propagation Plan"; the literature behind it is in
`reports/Urban canyon propagation methods.md` and `reports/Urban propagation survey 2026-10-11.md`.
What to keep:

- **A table of 0 dB with no facades is the run without one, byte for byte.** The direct path goes
  through `spectral_level_along` with the table's loss as `extra_loss`, and subtracting 0.0 is
  exact; with no valid reflection `urban_level` returns that level as it came, without an
  energy round trip; the bound multiplies by `max(1, ...)`. `urban_table{,_plate}_open_is_open_ground`
  compare the files with `compare_files`.
- **First-order reflections are paths, not table entries.** Each is
  `straight_ray(source, mirrored microphone)`: the sphere is read along the line leaving the
  source, so it gets its own directivity. Putting them in the table, keyed on the arrival
  direction alone, roughly doubles the error at 30 to 120 m (FIDELITY, Urban receivers).
- **An image's ground model is evaluated exactly** (no `GroundTable`): a table's extent follows
  the run's observers, and a mirrored microphone lies outside it. The vertical mirror keeps
  heights, so the ground geometry of the unfolded path is exact.
- **An image's separation is mirrored back** before the plate reads its azimuth: the sound
  crosses the plate travelling the mirrored way. `urban_table_plate_golden` fails without it.
- **The bound** multiplies `path_gains` by the table's largest gain, and the direct path and the
  images share one bound geometry, built by the eigenray rule (largest numerator over the
  harmonic sum of squared ranges). The image term has no CTest that fails without it (see
  [testing.md](testing.md#mutation-guarded-tests)); keep it, since a strongly directional sphere
  can make a reflection louder than the direct path's bin bound.
- **The roof edges diffract by UTD, and the table sums every ray coherently** (`coherent_ratio`,
  `canyon_rays` in `urban.cpp`, summed by `band_averaged_energy`, which `src/building_paths`
  shares; the kernel is `src/wedge`). Continuity across a
  shadow or reflection boundary is a matter of phase -- the diffracted ray makes up the jump of
  the specular ray it shadows only when the two add as complex amplitudes -- so the specular
  families and the diffracted rays share one sum, band-averaged pair by pair (the exact sinc of
  each path difference) over three sub-bands. In energy the field jumped at every boundary. What
  to keep:
  - **The run's own rays are left out of the sum, not subtracted from it.** Subtracting their
    band energy from the total left the cross terms with the rest, which went negative in
    interference nulls; the run adds its reflections in energy, so their interference is not
    the table's to carry.
  - **Overhead is guarded:** cos(pi/2) is 6e-17 in floating point, and a ray that steep leaves
    its trace 1e17 m up, taking its phase with it (13.8 dB wrong at 63 Hz before the guard).
  - **Pairs whose sinc is below 1/30 are skipped**, the rays sorted by length; that and computing
    each ray's phase once per sub-band took a seven-receiver build from 139 to 40 s of CPU,
    and moved the specular-limit test's mean by 4e-5 dB. Directions are threaded.
  - **Each pair's sine comes from the rays' own**, sin(h dL) = sin(h L_i) cos(h L_j) - cos(h L_i)
    sin(h L_j), so a pair costs a multiply-add and a divide, two pairs at a time through GCC/Clang
    vector types (the compiler judges a loop of double divides not worth vectorizing; 0.5-0.7 ns
    a pair against 1.6 on an M4 Pro). **Below h dL = 0.1 the sinc is its series instead**, and that
    is load-bearing: the difference of products carries the absolute error of h L itself, 1e-12
    at 10 kHz a kilometre from the frame's origin, and a district table moved by 2.4e-4 dB in one
    direction until the near pairs took the series. With it every table is within 2e-10 dB of
    the per-pair `sin`, and `the band-averaged sum is each pair's own sinc` holds it there
    (without the series it fails by sixteen times its tolerance).
  - **The kernel's boundary tie:** on a shadow or reflection boundary the geometrical-optics test
    and the coefficient's limit must take the same side (`hard_wedge_field` moves theta off a
    boundary once); in floating point theta + theta0 can land a hair under pi.
  - **The transition function's asymptotic series runs to its sixth term** past X = 50, meeting
    the Faddeeva form within 1e-8.
  - `--no_diffraction` keeps the specular rays in energy, which the separability study and the
    four-family tests compare against specular lattices.
- **The canyon table's specular part comes from four traces, not a lattice sum.** A plane wave
  reaches a point by at most one ray per cross-canyon direction and per way of arriving, from
  above or off the floor, so `canyon_ratio` traces four rays back from the microphone and its
  floor image; `specular_paths` (brute force, with the walls carried below the floor for the
  floor's rays) is the test reference only.
- **The floor is the table's, not the ground model's, inside a canyon.** The first version left
  the floor reflection to the ground model, as if it shared the direct ray's fate. The wave
  solver showed otherwise: from the floor image, 2 z_r lower, a ray clears the roofs after a
  different number of reflections, and at low elevations it is often the only one open. At 15
  and 30 degrees over a 20 m canyon the table read 20-25 dB of shielding where the wave solution
  has 2-4 dB. The table is therefore relative to open ground's two rays, 1 + rho_f.
- **A first-order reflection's floor partner is checked at run time**
  (`first_order_floor_partner`, against `urban_walls_below`): the reflection is evaluated with
  the ground model only when its partner is open too, and in free field otherwise, and the
  table holds the partner exactly when the run does not. Leaving the check out overcounts by
  3 dB in the directions where the two differ; both pinned goldens fail without it.
- **Receivers are matched by position** (1 mm), because the observer reader keeps no ids. The
  table records its frame bearing, and the run refuses another: azimuths are compass bearings
  and the geometry is local.
- **Stage-1 refusals**: `--terrain`, `--refraction`, `--propagation_table`, `--db_cutoff`.
  Lifting refraction means looking the table up along each eigenray's arrival direction;
  `--db_cutoff` needs `urban_bound_geometry` in `is_culled`.
- Tests: `[urban]` in `tests/test_urban.cpp` (two-family count and first-order images against
  `specular_paths`, the coherent field's continuity, specular limit and shadow, the table's
  lookup, bound and file), the
  hidden `[.urban_separability]` study, and the `urban_table_*` CTests (open-sky identity,
  windowed against marched for canyons, plates and the +25 dB bump table, the pinned goldens and
  the refusals).

### LoD1 buildings (`src/buildings`, `src/building_paths`, `src/diffuse_facades`, `src/geojson_buildings`, `src/plan_grid.hpp`)

`niceops_urban_table --buildings` builds the table among footprints extruded to flat roofs:
`Buildings` holds the geometry and the visibility test, `BuildingReceiver` one microphone's rays and
its ratio, `building_gains` the table. It is the ideal canyon's model for any arrangement, and two
long buildings reproduce `canyon_ratio` to 1e-9 dB (FIDELITY, Buildings). What to keep:

- **The visibility test is exact for prisms, so there is no ray-tracing library.** A segment is
  blocked when a stretch of it longer than `Buildings::contact` (1 micrometre) lies strictly inside
  a footprint below its roof: breakpoints at every side it meets and at the two heights, a
  midpoint test between them. Touching is not blocking -- a ray ending on a wall where it
  reflects, leaving an edge where it diffracts, or grazing a wall -- which is what lets every
  path be tested with the same call. Embree, which the plan named, is left for LoD2 meshes.
- **The visibility queries walk a plan-view grid** (`PlanGrid`): `blocked` the buildings' boxes,
  `first_hit` the facades', cell by cell along the line, each item tested once per query (a
  per-thread stamp) and `first_hit` stopping at the cell that holds its best hit. The grid only
  prunes; `the plan grid finds what trying everything finds` holds both to the exhaustive answer
  on 4000 random queries. It is what makes a real district's thousands of buildings affordable.
- **The floor is a mirror world**, `below_floor`: every building from `-top` to `top`, the same
  trick as the canyon's walls carried below the floor, so a floor-reflected ray is one straight
  line from the microphone's floor image. Below the floor a rising ray can enter a mirrored
  building through its mirrored roof -- the incident wave coming down through a real one -- which
  the facade search cannot see, so `trace_specular` tests every leg there whole. A microphone
  lower than the roofs around it has its floor image above their mirrored roofs, so its rays did
  not move (every golden held); a high patch's floor image is below some of them.
- **Facades run clockwise around their building**, so that `facade_side` is positive outside, as
  `first_order_image` and the run expect.
- **Specular rays are traced, not enumerated.** A vertical facade turns only a ray's horizontal
  component, so a ray that leaves along u reached the microphone along one of u's reflections in
  the facades' orientations: those candidates (`max_turns` changes of orientation, at most
  `max_departures`, breadth first) are traced deterministically through up to `max_order`
  bounces and kept when they leave along u. Parallel facades' mirrors cancel in pairs, so a grid
  has four candidates however often a ray bounces up a street; at 1.5 degrees a ray needs about
  fifty bounces to climb out of a 25 m street, which is why `max_order` is 64. A rising ray
  enters a building only through a wall, so the first facade it meets is its first obstruction.
- **Diffraction points are on Keller's cone, and rays whose point falls off the edge are
  dropped.** The edge's own faces are already in its coefficient, so a diffracted ray never
  reflects first in a facade containing its edge; `diffracted_path_open` refuses that by its
  in-front test without a special case. Reentrant corners do not diffract.
- **Oblique incidence is the 2-D wedge at k sin(beta) and r = s sin(beta)**, over
  sqrt(s sin(beta)); two tests hold a roof edge and a vertical corner to `hard_wedge_field`.
- **The images come from plan-view beam tracing, breadth first.** Each image carries the part of
  its last facade through which the microphone can be reached; a further facade counts only
  where it lies in that beam, and each image keeps the edges whose plan crosses it. The test
  `the beams find every diffracted ray the unpruned search does` compares the rays with
  `BuildingPathOptions::beams` off. Breadth first, so that a search stopped at `max_images`
  has dropped the highest orders rather than whole branches; depth first, a stopped search
  read 0.36 dB low in energy on the FIDELITY district. The facades' mutual visibility is a
  sampled prefilter (20 points on each face, in the mirror world), a heuristic that a gap
  narrower than its samples could fool; every ray is still tested exactly. It is worked out
  when first asked, after the beam's own test, with the samples in the same order as before:
  most pairs are never asked, and a point source's two trees per pair went from 6.3 to 2.1 ms.
- **Everything that does not depend on the direction is built once per microphone**: the
  frames, the orientations, the edges within `radius`, the images and their edges, and each
  (image, edge, floor)'s foot on the edge (`BuildingReceiver::Foot`): its distance along and
  across, the range of cot(beta) that puts a Keller point on the edge, and whether the image
  can lie in the wedge's air at all, which does not depend on where along the edge it is seen
  from. A direction costs the traces, the Keller points and their tests; 85% of the
  (image, edge, floor) triples miss the edge, and they now cost two comparisons. The prefilters
  are a hair generous and the exact tests stand behind them, so the tables did not change by a
  bit.
- **A ray's wedge angles are worked out once** (`WedgeAngles`, `wedge.hpp`): the coefficient's
  four cotangents and cosines depend on the angles alone, and a ray is evaluated at 93
  sub-bands. `hard_wedge_coefficient(n, theta, theta0, k, r)` is the split form to the bit. With
  the pair sum and the foot prefilters this took a builder receiver from 1.75 to 0.77 s on the
  49-block benchmark district (`niceops_urban_table`, 72 by 33 directions, 31 bands).
- **The run's facades are the table's exclusions.** The tool writes the facades within
  `--facade_radius` that face the microphone, and passes the same list as `run_facades`, so the
  table leaves out exactly the first-order reflections (and floor partners) the run evaluates.
  Passing an empty list instead fails both `urban_buildings` goldens.
- **Grazing and overhead are clamped as in the canyon**: elevations under 1e-3 rad are taken at
  it, and a direction within 1e-9 of vertical is vertical.

- **GeoJSON footprints are projected by the run's own transform** (`src/geojson_buildings`,
  `LocalFrame::track_aligned` with the run's origin and bearing), so a building means the same
  place as a geodetic track. Rings lose their closing vertex and any vertex within a millimetre
  of the last; inner rings are courtyards. Every ring keeps the solid on its left (outer
  counterclockwise, holes clockwise), so one code path makes the facades, edges and visibility of
  all of them. A feature raised off the ground is left out, not extruded from the floor.
- **Shared sides are found exactly** (the same two vertices, run the other way, to a micrometre),
  as OpenStreetMap terraces share nodes. A side shared only in part keeps its facade and edges;
  `Buildings::buried` then keeps an edge point or a patch that lies inside the neighbour out of
  the rays. The lower roof's edge along a taller wall would otherwise diffract as a convex
  wedge, and it is reentrant.
- **Round-off is the enemy of aligned geometry, three ways**, each found by reading the demo
  district back from GeoJSON, where 1e-8 m of round-off made walls that were parallel or square
  1e-10 rad off:
  - Facades whose normals agree within the orientation tolerance share one orientation, for the
    departures and for the traces' reflections both (`Buildings::reflect`). Taken one by one, the
    near-parallel walls multiplied the candidate departures to their cap, and a seven-receiver
    build that takes two seconds was stopped after ten minutes unfinished. The tool's default is
    1 degree, not the 1e-5 rad round-off needs, because real streets are not exactly parallel
    either: a ray bouncing up a street whose walls diverge by half a degree turns by a degree a
    pair of reflections, its departure is then a composition dozens of mirrors deep, and with
    exact orientations those long chains are never candidates (`a street half a degree out of
    parallel keeps its long bounce chains`: none exactly, 28 reflections within a degree). The departures are also made only
    in the orientations of the facades within the radius (`orientations_near`).
  - Departures within 1e-6 are one, and a trace escapes "along u" within 1e-7. Two orientations
    a hair off square compose into a rotation by round-off, not the identity, and the direct ray
    was found three times. A path turns its departure rigidly, so no ray can be found twice.
  - The table is evaluated 1e-6 rad clockwise of its azimuths (`building_gains`): a direction
    exactly parallel to a wall or through a corner is a knife edge of the model, and a grid
    district samples those lines.
- **Scattering walls are radiosity, solved backwards** (`DiffuseFacades`, THEORY §6.5). The
  form factors and the microphone's gathering depend only on the geometry and the band, so
  w = S (I - F^T D)^-1 g is solved once per receiver and band, and a direction costs only its
  irradiance. `energy_forward` solves the other way round, per direction, as the test's
  reference. What to keep:
  - **Every specular reflection scatters its share.** A patch's irradiance is gathered from every
    specular arrival, traced from the patch with the microphone's own tracer
    (`specular_departures`, `trace_specular`) and carrying rho (1 - s) per reflection. Lighting the
    walls with the plane wave's first hit alone loses the share of every later reflection: 2-3 dB
    at 4.5 to 13.5 degrees on the FIDELITY district.
  - **Diffuse energy stays diffuse, with the whole rho.** A field arriving from every direction
    reflects into every direction, specularly or not.
  - **The microphone's two paths from a patch interfere** band-averaged as open ground's rays do,
    so a flush plate gets the same doubling from the diffuse field as from the direct sound;
    between patches the floor path adds in energy.
  - **With s = 0 nothing changes, byte for byte**: the tool builds no diffuse field, and
    rho (1 - 0) is rho. The run's first-order reflections carry rho (1 - s);
    `urban_buildings_diffuse_golden` fails if they carry rho.
  - Patches under three of their sizes from a receiving point are summed over 4 by 4 parts;
    each patch is lit, and seen, from its centre, so shadow edges are resolved to a patch.
  - **A point source lights the same patches** (`point_irradiance`, `niceops_point_field
    --scattering`): the source's specular paths to each patch are the patch's own images
    (`image_tree` from its centre, to `--max_reflections`) checked by `folded_path_open`, exactly
    as the microphone's specular rays are found, each cos(incidence) / L^2 times its rho (1 - s)
    and the floor's. An image whose first reflection is the patch's own wall is skipped: that is
    the wall reflecting what arrives, not a second arrival. Only the active patches (weights not
    negligible) are lit, in parallel over patches. The weights are built once per microphone for
    every source, and so are the patches' images (`prepare_point_lighting`): 20 positions of
    DigiCity's loudspeaker at its four microphones take 2.2 s on one thread, against 1.0 s with
    specular walls, where tracing the images afresh for each source took 20 s. The lighting paths
    and the gathering take the air's absorption (`DiffuseOptions::air`); the patches' exchange
    does not.
  - **Only walls in reach scatter.** A wall beyond `--diffuse_radius` of a microphone has no
    patches to return its share, so for that microphone it keeps its whole rho: the rays' reflections
    are per microphone (`reflection_for`), and `point_field_scattering_out_of_reach_is_specular`
    holds a run with none in reach to the specular run. `point_energy` is tested against Lambert's law for a source and microphone
    together in front of a wall, rho s / 2 d^2 (0.2% at patches of d / 5).
  - **The coefficient is the relief's, band by band** (`--scattering_spectrum`, from
    `tools/profile_scattering.py`): a constant s is a poor stand-in, the relief scattering nothing
    where it is shallow against the wavelength and nearly everything above. A building's own
    coefficient then weights the spectrum, so one run holds plain and relieved buildings.

### A point source among buildings (`src/point_paths`, `src/double_wedge`)

`PointSourceField` is the tables' ray model for a source near the buildings (`niceops_point_field`,
THEORY §6.5). It shares `image_tree` and `folded_path_open` with `BuildingReceiver`: the
microphone's images and the source's are the same beam-traced trees, and every leg is tested
exactly. What to keep:

- **Double diffraction over a shared face is an aperture integral, not the cascade.** Measured
  against boundary elements (`tools/thick_barrier_bem.py`, FIDELITY), the cascade of two UTD
  coefficients is 1-4 dB low over a container stack with the source 18° below the roof and up to
  11 dB low at 1.6°, and that does not fall away with frequency. The first edge's wave reaching
  the second is a transition field, and no choice of distance parameters makes a ray of it:
  - Kouyoumjian and Pathak's own parameters (`s1 s2/(s1 + s2)`, `s2 s3/(s2 + s3)`, the
    cascade) are the ones that compensate the second edge's jump exactly, and still miss.
  - Whole-ray parameters, each edge's taken from the whole ray before and after it
    (`s1 (s2 + s3)/S`, `(s1 + s2) s3/S`), with the ordinary spreading gave a
    frequency-independent step of 0.5 to 16 dB at the roof plane, and up to 24 dB high against
    the BEM: the spreading has to change with them.
  - The first edge's parameter from the whole ray and the second's from the cascade (a hybrid)
    is up to 15 dB high.
  The aperture integral holds the transition itself; its end point is the cascade, and the
  Keller-ratio correction `c` puts the first edge's large-angle coefficient right. Slope
  diffraction is no help here: over a hard face the incident field's normal derivative is zero.
- **Above the face's plane, subtract the aperture's own unbounded field, not the UTD ray.**
  Subtracting the first edge's UTD ray (its diffraction and the source's wave) was 0.1 dB from the
  BEM near the plane but 10 dB high for microphones well above it: the correction `c` is right
  only toward the face, so the aperture's picture of the first edge's wave in other directions
  differs from the UTD ray by as much as the second edge's whole effect. The same integral with
  the face continuing carries the same error, and the difference is what the second edge does.
- **The two evaluations are weighted by depression angle, not averaged.** Each is within a few
  tenths of the other; a plain average switches to one of them when an end crosses the plane, a
  0.3-1 dB step. The angle weights go smoothly to the evaluation that remains, and swap with the
  ends, so reciprocity holds to round-off.
- **The 2-D field converts to 3-D per ray, by sin(beta)/sqrt(l_t).** At beta = 90° that is exact
  for every kind of ray, which is what lets a 2-D BEM test the whole 3-D implementation
  (`over a long block the point-source field is the boundary-element solution's`); it also
  reduces the shared-face form to the 3-D cascade at its end point.
- **The test's reference values come from `tools/thick_barrier_bem.py`** at 60 elements a
  wavelength, converged to 0.05 dB on those points. Rerunning it reproduces them exactly.
- **`source_reflections`** exists for the far-field test: the plane-wave builder leaves out edges
  lit by a reflection or the floor, and the comparison needs the same ray set.
- **Do not range-for over `PointSourceField(...).rays()`.** The temporary dies before the loop
  runs (C++17 does not extend it through a member call); two tests did, and read garbage.
- `tools/building_input.hpp` is the district and microphone reading both urban tools share; a
  change to the CSV or GeoJSON rules goes there once.

**Where its time goes, and what not to undo.** On the 49-block benchmark district (a drone 30 m
up, 185 street microphones, 24 bands) a pair took 1.55 s of CPU and takes 0.04 s, 88% of it the
bands; DigiCity's 1427-point grid (17 bands) went from 0.33 to 0.043 s a pair. In order of what
it bought:

- **Each end proves which stretch of each edge it cannot see** (`visible_stretch`,
  `edge_visibility.cpp`), once for all its partners (`PathEnd::View`), and both diffracted
  searches drop a stationary point outside the stretch before checking its legs; the
  twice-diffracted search's reach tests bound the point by the stretch rather than the edge, so
  most such pairs are never solved. For an image to one reflection the path from the edge's
  point q unfolds into the segment from q to the end's mirror; a convex building, shrunk by
  1e-4 m, blocks it for one interval of q (a polygon clipped in two parameters, u along the
  segment and w = u t), and the stretch runs from the first point the intervals of the buildings
  near leave uncovered to the last. On the district a source's search over its 185 microphones
  costs 3.3 s on one thread without it and 1.1 s with it; on a generated city of 900 blocks, a
  source 40 m up over 40 street microphones, 10.7 and 2.2 s. The rays are the same to the byte
  there, on the district and on DigiCity's grid. What it rests on:
  - **It errs one way only.** A point outside the stretch is one `folded_path_open` refuses:
    `a point of an edge proved hidden is hidden` samples every 2 cm and finds every point it sees
    inside the stretch, and `pruning the diffracted paths loses no ray` holds the rays the same
    with pruning off. Footprints that are not convex, or have courtyards, are left out, which
    only widens a stretch; images of two facades or more are not tried.
  - **The union and the stretch are what pay.** Proving whole edges hidden one building at a time
    prunes only 7-25% of them: a source 30 m up sees some of all but 8% of the edges within
    150 m, and a street microphone some of 45-72%. Covered by every building's interval at once,
    per image and floor, the test's edges that no sample sees are proved hidden 97% of the time,
    and the stretch excludes 84% of the unseen points on the rest.
  - **The buildings come from the plan grid** (`Buildings::near_box`), not a scan of all of
    them, which a city of thousands would make cost more than it saves.

- **Twice-diffracted pairs are ruled out before they are solved**, which removes 90% of them
  before the stretches above. Whether the source's image lights the first edge, and the
  microphone's the second, does not depend on where along the edge (`PathEnd::View::lit`, a hair
  generous); and a stationary point lies between its image's foot and the other edge's feet, at
  the fraction off_p / (off_p + off_q) of the way, which bounds it from the edges' spans
  (`keller_can_reach`, the cheap `keller_between` first). The search runs over the lit edges'
  lists, not every edge with a test. `PointPathOptions::prune` switches these and the stretches
  off, and `pruning the diffracted paths loses no ray` holds the rays the same to the bit with
  and without (a slack of -0.5 m in the reach test fails it, as does the stretch pulled in by
  0.5 m); the ray dumps of DigiCity's grid and the district were byte-identical too.
- **The fixed point is found by Steffensen's method**, in closed form for parallel edges
  (`parallel_stationary`; the map is affine there), and a pair that does not converge in 100
  steps is no ray. The plain alternation it replaced took 16 steps on average and stopped at
  200 whether or not it had converged, and near a corner where a roof edge meets a vertical one
  it was still creeping toward the corner -- 0.36 mm a step, 4 cm short -- when it stopped. Those
  iterates were accepted as rays with s2 of a few centimetres, and their spreading made them
  strong (what rejecting them moved is in
  [HISTORY](../HISTORY.md#urban-receivers)). A corner's own diffraction is left out of the
  model, and these were not it. Steffensen's step falls back to the map itself where it does not
  bring the residual down, and reuses the map it has already taken there.
- **The checks after a solve run cheap to dear**: the angles, then `buried`, then the legs. The
  first leg is blocked for half of what converges.
- **The shared-face aperture is integrated by Filon's rule** (`double_wedge.cpp`), 6 to 24 times
  fewer integrand evaluations than the midpoint rule at a 24th of a wavelength and ten times its
  accuracy ([FIDELITY](../FIDELITY.md#a-point-source-among-buildings) has both against a 192nd of
  a wavelength). Four things in it are load-bearing:
  - **A segment's end nodes sit a hair inside it.** The wedge coefficient jumps at a cut, and on
    the cut it takes one side's value; a node exactly there put the other side's value into a
    panel's quadratic, up to 28 dB off. A microphone in a roof's plane (theta = pi) puts both
    boundaries at the edge, t = 0, so the first node is nudged too (6.5 dB off without).
    `the aperture's Filon rule is the fine midpoint rule's` takes both cases, and fails by 6.8 dB
    without the nudge.
  - **The panels grade from the edge and the cuts only**, at 0.15 of the distance plus 1% of
    min(w, s1); the taper's start and the aperture's end are panel ends but not features.
    Grading from them as well cost twice the evaluations for nothing.
  - **The moments choose series or recurrence by the panel's phase in magnitude.** Every path
    length that counts grows along the aperture, so the phase is positive in every case met; a
    panel along which one shrank would otherwise take the series far past the radian it holds
    to. The Filon test's obtuse corners (n = 1.25) are what take the wave off the second
    wedge's other face, which right-angled faces never reach.
  - **The end nodes' hair matters at high k.** A geometrical-optics wave's stationary point lies
    on its cut, where its amplitude per unit length goes as 1/sqrt(u − u_c); how much of it the
    first panel misses depends on how far in the node sits. 1e-7 of each wavenumber's own short
    segment is fine; 1e-7 of a 200 m one is not, at 10 kHz, which is why
    `shared_face_geometry` (the `--gpu` path's aperture, worked out once for every wavenumber)
    takes a hair of the shortest aperture instead ([gpu.md](gpu.md)).
  Its second wedge's coefficient comes from a `WedgeSweep`: theta is fixed along the aperture,
  so each term's angle is theta's part plus theta0's, by the angle-addition formulas, two sincos
  calls a node instead of twelve calls. Within 1e-8 rad of a term's pole the sum of products
  cancels, and that term takes `wedge_angles`' own path; elsewhere the two agree to 1e-9
  (`the wedge sweep is the coefficient to round-off`).
- **The pair sum's vectorized difference form** (above, under the canyon): with a thousand rays
  a pair it is still a quarter of the district's time.
- **`band_levels` sorts indices**, not a copy of the rays and their paths.
- **Whether a direction is in a wedge's air is decided by signs** (`Edge::air_by_sign`): its
  components along `air` and `air1`, the two faces' normals toward the air, settle it wherever it
  lies 4e-7 rad or more from a face -- wider than onto_face's 1e-7 snap, the widest slack any
  caller gives the angle -- and the angle is taken only inside that margin, and for the rays kept.
  Widen a caller's slack past it and the signs would decide what the angle would not. The ray
  dumps did not change by a byte.
- **Each edge pair's geometry is worked out once a field** (`EdgePair`, two `Side`s, in a table
  hashed by the pair and cleared by generation, so that its memory follows the pairs met rather
  than the square of the edges near): each segment's span about the other's line, which the
  pre-checks read, and its points as that line sees them, affine in their parameter. A pair was
  met eight times on average, once per image pair. The stationary point's iteration reads the
  source's and microphone's places from their ends (`PathEnd::View`), so `keller_point` takes
  one square root an edge, not two, and a once-diffracted ray's Keller point takes none: both
  ends' places are the views', worked out once an end rather than once a pair of images.
- **A moving source traces each end once** (`PathEnd`, `--sources`). A field's image trees depend
  on the pair only through which near edges each image keeps, so an end traced with the edges
  near every end it will meet (its `partners`, `edges_near` over their points, ascending), each
  pair then keeping those near it, gives the pair's own lists in the pair's own order. A pair
  whose edges an end was not traced with is refused, as are ends in the wrong roles or traced with
  other options, and `a field from ends traced once is the pair's own` checks all of it. The
  single-pair constructor builds its two ends the same way, so every field goes by this road.
  Images deeper than the edge order keep no edges (`image_tree`'s `edge_order`): nothing diffracts
  through them. That order is `max_edge_reflections`, but at least one wherever double
  diffraction is on, because the twice-diffracted search takes the microphone's images to one
  reflection whatever `max_edge_reflections` is. Without the floor of one, a run at 0 loses two
  thirds of its twice-diffracted rays; `twice-diffracted rays reflect at the microphone whatever
  the edge order` holds it.
  The ends are cheap -- about 5% of a pair, 2 ms against 42 ms of CPU on the district -- so
  what `--sources` buys over a run per source is mostly being one process, with every worker
  busy to the end. The tool runs the sources a block at a time (64, or four per worker) to bound
  what it holds: the block's ends, levels and, for `--rays`, ray text, written out before the
  next. Within a block the pairs go microphone by microphone, which keeps the workers busy
  however few microphones there are, and each worker keeps the microphone's end while its pairs
  last; across a block boundary it is traced again, which the 3% makes noise.
- **On a boundary the coefficient takes the side its ray is kept on** (`term_angles` in
  `wedge.cpp`): the incident wave's where |beta| < pi, face 0's reflection's where beta < pi, face
  n's where beta > (2n - 1) pi. The path tests count a path through an edge as open, so a ray
  grazing the edge exactly is kept, and its edge's term must be the lit side's -1/2. A hair above
  every boundary, whatever its side, put 1 + 1/2 on two of the four: 9 dB on a corner's shadow
  boundary. `the field is continuous across a corner's shadow and reflection boundaries` holds a
  tall box's corner to 0.5 dB between neighbours 2.5 cm apart. `hard_wedge_field` instead moves
  its direction off a boundary once, for both of its tests.
- **Continuity is the model's own test of its ray set.** Every boundary needs its partner traced
  (THEORY 6.5), and where the orders run out the field steps. The validation harness's
  `urban/continuity.py` sweeps a city for steps and classifies each by the ray that appears or
  vanishes and why; [FIDELITY](../FIDELITY.md#a-point-source-among-buildings) has DigiCity's.
  Raising `--max_edge_reflections` moves the gap up an order rather than closing it.
- **The air's absorption is a factor on each ray's amplitude** (`sub_band_absorption`, Np/m at
  each sub-band's middle frequency), taken where `sub_band_problems` forms the amplitudes and
  only when one is given, so a run without it is the same bits as a build without absorption (the ray
  dumps and levels of DigiCity's grid and the district were compared). On the GPU each problem
  carries its coefficient and each ray its length (`BandProblem::absorption`, `GpuRay::length`),
  and both amplitude kernels multiply by exp(-a L). `--rays` gives each ray's level at the band
  centre, so its absorption is the band centre's. `the air absorbs each ray along its length`
  holds the levels to the closed form for a lone ray.
- **A ray's reflections are taken once a band** (`reflected_amplitude`, both models): the square
  roots of its facades' and the floor's reflection, which no sub-band changes; the coefficient
  is the only per-wavenumber part. `ray_amplitude` and `point_ray_amplitude` compose the two,
  so the table builder and the point-source model share one spelling of each.
- With the first two (and `folded_path_open`'s scratch kept per thread, and single
  diffraction's edges matched by one pass along both sorted lists) a pair's path search fell by
  about a third and the whole pair by 11%. The bands are now two thirds of a pair, the shared
  faces' integrals half of that and the pair sum most of the rest. The Filon panel's moments
  below a radian are a fixed polynomial in phase^2 (`filon_series`), not a series summed term by
  term, which took the shared-face integral down by a further third.
- **The transition function is a table below X = 50** (`utd_transition`, `wedge.cpp`). Its
  argument is real, so F / sqrt(X) is w along a line in the upper half-plane, smooth and slowly
  varying: a polynomial of degree 10 in sqrt(X) on each of 32 intervals of [0, sqrt(50)], fitted to
  the Faddeeva form at each interval's Chebyshev nodes the first time it is wanted. It agrees with
  that form to 7e-12, which is the form's own error along the line (SciPy's `wofz` puts it at 5e-12
  below sqrt(X) = 1 and near 7), and costs 3 ns a call against 15. The Faddeeva form and the
  transition around it had been a third of a `--buildings` run; the demo approach went from 70 to
  54 s. It is built on first use, not at load, because the Faddeeva form reads a constant that
  `faddeeva.cpp` initializes at load, in an order across translation units the language leaves
  open. `the tabulated transition function is the Faddeeva form` holds it at 2e-11; a degree of 5
  fails it at 8e-9.
- **A shared face's nodes are worked out once for all of a ray's sub-bands** (`shared_face_spectrum`;
  `Aperture` and `NodeTable` in `double_wedge.cpp`). Each wavenumber keeps its own nodes, ending at
  its own taper and extent, but they are graded from the same edge and cuts, so most of them fall
  on points, to the bit, that another wavenumber has taken already. A node's geometry -- its path
  lengths, rates and weights, and the second wedge's angles and distance parameter there -- does
  not depend on the wavenumber, and is found by its t in a table instead of being worked out again.
  The arithmetic is the same either way, so the spectrum is `shared_face_field` at each wavenumber
  to the bit (`the shared face's spectrum is its field at each wavenumber`, 72 sub-bands over the
  Filon test's faces), and the demo approach's levels did not change by a byte: 54 to 37 s. A
  node's wave counts follow from the cuts that bound it, so they agree wherever its t does; the
  table compares them anyway, and taking the comparison out passes the test.

Measured and left alone:

- **The GPU's aperture on the CPU.** `shared_face_geometry` works one set of nodes out for every
  wavenumber, ending panels at every one's taper and extent (it is more accurate, 0.006 dB from the
  fine midpoint rule against 0.017). On the CPU the demo approach took 70 s with it against 54 s
  without: each wavenumber then integrates over the other wavenumbers' panel ends as well, and the
  extra nodes cost more than the geometry they share. Sharing the geometry of each wavenumber's
  own nodes (above) keeps the CPU's integral and is faster than either.
- **Phases measured from the reference path**, as the GPU's are. The profile showed `sincos`
  reaching into a second routine, which looked like the reduction for phases of 1e4-1e5 radians;
  measured from the reference, the phases are small and the run took the same time (37.5 s against
  37.3). The calls are many, not slow.
- **Dropping weak rays does not pay.** Leaving out of each band the weakest rays, together under a
  share of its energy (ranked at the middle sub-band, whose amplitudes the kept rays reuse), saved
  1.3x at a share of 1e-4, 1.5x at 1e-3 and 1.7-1.8x at 1e-2 on the district and DigiCity, for
  level changes of 0.09-0.21, 0.72-0.84 and 1.4-3.5 dB at worst (0.04-0.06, 0.17-0.20 and
  0.69-0.72 dB at the 99th percentile). Two things cap it. A ray with a thousandth of the energy
  can still be a near-twin of a strong one -- the same path off the floor, a hair longer -- and
  its cross term is then sqrt of its share, not its share: ranking by energy plus that overlap
  (the rays within pi / h of it) holds the change to 0.006-0.14 dB, but then little is dropped
  and the saving is 1.1-1.2x. (Those savings were measured with the twice-diffracted search the
  larger cost; with the stretches above the bands are most of a pair again.) Nor can an energy
  test prune that search cheaply: those rays are relatively strongest in the lowest band, where a
  60 dB floor keeps most of them, and bounding a ray's coefficient costs what the visibility
  tests it would skip cost. A shared face's ray cannot be ranked by its
  cascade either: on the district its uniform form is up to 53 dB stronger, a source exactly in
  a roof's plane.
- **Nor does the band average's integral form on the CPU.** Where every pair of rays is within
  the pair sum's reach (the bands up to about 200 Hz on the district), Gauss-Legendre over the
  sub-band, as `--gpu` takes it, would agree to 1e-12 dB, but the vectorized pair sum takes 0.42 ns
  a pair there -- 0.54 s for 1.3 billion pairs over the district's 185 microphones -- and the
  integral would take 49 million sincos calls for the same problems.
- **Neither do cheaper ray sets.** Twice-diffracted rays off the floor at both ends carry 0.1% of
  any band's energy but are near-twins of others: leaving them out moved levels by 1.5-2.7 dB.
  Without double diffraction the district runs 13x faster and some microphones receive nothing;
  those rays carry up to 41% of a band elsewhere. A 100 m `--radius` moves levels by up to 42 dB.

- **Roofs reflect, one horizontal surface to a ray** (`PointPathOptions::roofs`,
  `folded_roof_path_open`, `Buildings::on_roof`). Among vertical walls a ray that leaves a roof
  keeps rising, so the floor and the roofs exclude each other and a roof is the floor at its own
  height for that ray; the facades' mirrors commute with it. What to keep:
  - **The bounce is tested on the roof to round-off, not by `inside()`**, which keeps 1 µm off a
    footprint's sides. The roof's ray has to appear exactly where its edge's diffracted ray changes
    sign, at the roof's reflection boundary; tested by `inside()` it stopped short and the field
    stepped by 0.19 dB 1 µm off the boundary.
  - **The edge must stand above the roof.** A roof's own edges carry its reflection in their wedge
    coefficients, the roof being one of the wedge's faces; a ray off the roof to them would count it
    twice.
  - **A roof counts as one of a side's `max_edge_reflections`**, and its candidates come from
    `near_box` over the plan box from the image to the edge (or from the source to the
    microphone's image), sorted so the rays come in a fixed order.
  - **Roof rays skip the visibility proof**: their views take the whole edge, since the proof knows
    facades' apertures and not a roof's outline. The exact path tests decide.
  - **They take UTD under `--edge_source`** (`edge_source_factors` returns nothing for them): the
    roof's outline bounds the stretch of the edge its image sees, which the integral does not know.
  - **`point_ray_launch` reads `source_plane`**, which floor rays set to 0 as well. Telling a floor
    bounce by the source image's sign could not tell a roof, whose image may lie above z = 0.
  - Cost: 4% on `--buildings`' district run (60.0 to 62.5 s).
  - Not done: twice-diffracted rays off a roof; roofs in the plane-wave tables, where a street
    microphone would see a roof's ray only through an edge lit by a reflected wave, which they leave
    out; and a reflection between two edges, which UE06's floor source needs (FIDELITY).
  - Tests: `a roof reflects as the floor does`, `the field is continuous where a roof's edge hides
    its reflection`, `on a wide roof the field is the direct ray and the roof's`, and `rays off a
    roof are reciprocal and leave along minus their length's gradient`.
- **`point_ray_launch` is minus the length's gradient**, and `a point ray leaves the source along
  minus its length's gradient` checks it so, by nudging the source, over every kind of ray; it fails
  without either floor rule. It is what a directional source is read along.

### In a run (`--buildings`, `src/urban_point`)

`UrbanDistrict::levels` gives every state's level at every microphone before the observer loop,
which marches them. What to keep:

- **The district is opaque to `niceops.cpp`.** `buildings.hpp`'s `Ring` is not `contour.hpp`'s, and
  both live in `niceops`, so the run's translation unit sees `UrbanDistrict` (a pimpl) and never
  `Buildings`. `building_input.hpp` moved from `tools/` to `src/` for it.
- **Sources go in blocks, their trees built in parallel, and pairs microphone by microphone**, as
  `niceops_point_field --sources` takes them.
- **Each ray carries the sphere's band amplitude along its launch** (`Database::spectrum_along`, at
  1 m), with the sphere's conventions: no spreading gain and no absorption inside its radius. That
  is what makes `buildings_far_is_open_without_ground` hold -- one building 1000 km off and a floor
  absorbing everything give the run without a ground model, at 1e-9 (in fact to the printed digit).
- **No peak search.** A coherent sum of a source's rays can exceed the straight path's bound by
  20 lg N, so every state is evaluated, as under the options that need a full history.
- Cost: 9.4 ms of CPU a pair; the demo track's 6,551 states at seven microphones take 37 s on 14
  threads (70 s before the transition table and the shared nodes above). Sampled, the shared
  faces' integrals are still most of it: roughly a fifth the transition function's table, a
  quarter the Filon panels with their sines and cosines, a fifth the rest of the aperture's loop
  and 6% the wedge coefficient; the ray sum is 7% and the path search 1%. What it changes is in
  FIDELITY.
- Tests: `buildings_far_is_open_without_ground`, the refusals `buildings_refuses_*`, and the
  launch test above.

**How its rays add** (`src/source_coherence`, `partially_coherent_energies` in
`src/urban_point.cpp`; THEORY 6.5). The model came out of the harness's `coherence/`, which says how
each part was measured; what to keep in the code:

- **The rotor's delay rides on the rays, not the pairs.** A ray carries its effective length
  L + (c s / Omega) psi, so the pair's path difference and delay are one difference, and the sine and
  cosine of the half width times it come from the rays' own (`band_averaged_energy`'s trick). The
  pair's wrap to +-pi is then a shift by one turn's length on the pairs whose raw azimuth difference
  passes pi, whose phase, sine and cosine are per sub-band constants. A pair needs no trigonometry.
  Doing it per pair (a polar and a sin per pair per sub-band) was most of the sum's time.
- **The decay is taken at the band's centre, once per pair per band**, as it was measured (on
  band coherency), not per sub-band, and a pair once below e^-14 is dropped for the bands above:
  its decay only grows with frequency.
- **No sort.** An early version sorted the pairs by separation to break early; at 724 rays that is
  a sort of 262,000 pairs per field, and the shrinking live list does the job.
- **The chord, not the angle.** |d_i - d_j| in degrees needs no acos, and a decay in the chord is a
  positive-definite kernel on the sphere, which a decay in the great-circle angle with an exponent
  above 1 is not; they differ by hundredths of a degree where the decay is not yet complete.
- **Buffers per thread, bound once.** The pair arrays live in a `thread_local` `PairScratch`, bound
  to references at the top: each use of a `thread_local` itself is a lookup on macOS, and in the
  pair loops that cost 15% more than the allocations it saved.
- **The interval's groups.** Rays with mean coherency 0.99 or more are one group (single linkage,
  per band), and each pair of groups adds half the square of its summed cross terms, each weighted
  by the root of its pair's spread. Treating every pair of rays as independent instead understates
  the variance wherever rays travel together: a flush microphone's direct ray and the floor's, with
  a facade's reflection and its floor partner, give four pairs that deviate as one, a factor of 4.
  A pair of rays each alone in its group is its own slot, outside the hash map.
- **Coherent mode keeps the old code path** (`PointSourceField::band_levels`), so its output is
  the earlier coherent sum to the bit; `partially_coherent_energies` under `RaySum::coherent`
  agrees with it to 1e-9 dB (`test_source_coherence.cpp`).
- **The default is the cone, not a fit** (`cone_coherence()`). Weighted by A-weighted power, the
  flight data cannot tell any cone up to 2.5° from the energy sum, and the energy sum assumes
  nothing of the source; the cone is there for the one emission (shared launches, diffraction
  partners). Its mean does not depend on the band, so it is formed once per field and the live
  list is the pairs inside the cone: the cone costs what the coherent sum does.
- **The rotor's delay comes from the database alone** (`measured_coherence()`). A vehicle file's
  blade count and sense were briefly taken as switching it on; but a delay is an assumption about
  the source's phase, so only a database that records it, where someone fitted it, applies it.
  The record's fits were made with its rotor's delay removed or not, and the defaults for what it
  leaves out follow the same.
- **Rejected** (measured in the harness): the delay at every frequency (worse than none above about
  160 Hz); an effective rotation scaled by kappa (0.75 in the passes, 1 in hover, within 0.03 of
  kappa = 1 in the score); the hexacopter's rotors as mutually incoherent sources on a ring
  (worse than energy for any radius); and a uniformly random phase beyond a window for the
  interval, which put a wide spread on broadband bands that add in energy.
- Cost on the DigiCity hovers: the cone the coherent sum's, the measured sum 1.4 times it; with
  `--phase_interval` 2.7 and 3.1 times (most of it the per-band spread and group slots, run over
  every pair).
- Tests: `[source_coherence]` (the delay's sign and wrap; coherency 1 and spread 0 at no
  separation; rays leaving together staying coherent; the measured sum within the coherent and the
  energy sums sub-band by sub-band; the variance against its formula for two rays; the record's
  precedence; the record's round trip through an exported sphere and its refusals), and
  `buildings_far_phase_*_is_the_level` (one ray: both ends are the level). panam's
  `test_directional_coherence.py` pins the estimator's sign against a synthetic turning pattern,
  and its `test_niceops_loads_databases.py` that NICEOPS reads the record and refuses a bad one.

### The edge-source integral (`src/edge_source`, `--edge_source`)

The Biot-Tolstoy-Medwin solution as Svensson's line of secondary sources (THEORY 6.5), in the
frequency domain, for the point-source model's diffracted rays. What to keep:

- **First order is integrated in eta, not z**, with the peak closed. Near a boundary a
  directivity term is a peak of width psi at eta = 0, and its integral is the closed arctangent;
  the rest is smooth. The two sides of the apex share eta's integrand. The denominator is
  2 sinh^2(nu eta/2) + 2 sin^2(nu phi/2): written cosh - cos it cancels at the apex within 1.5e-8
  of a boundary, where cos rounds to 1, and returned inf there.
- **eta(z) and z(eta) are written so nothing cancels.** cosh eta - 1 =
  (r_S b + r_R a)^2 / (r_S r_R (m l - a b + r_S r_R)), with r_S b + r_R a = (r_S + r_R)(z - apex) and,
  where a b > 0, m l - a b by its conjugate. Far out along an edge m l - a b is a difference of two
  numbers of 1e30 for a source 1e7 m off, which put 1e-3 dB on the plane-wave check. The edge's point
  at eta comes from a quadratic whose discriminant factors as (8 d r_S r_R sinh eta)^2.
- **Past a dozen radians of phase the single-wavenumber field goes along steepest-descent paths**
  (`descent`, 20 Gauss-Laguerre points). Every singularity of that integrand -- beta's poles, the
  branch points at h = 0 and i pi -- has R(h) real and at most R0, so it lies at least the window's
  phase from the Laguerre variable's axis however near a boundary the terms are.
- **For many wavenumbers, Filon's rule of degree four in R** (`edge_source_spectrum`): the
  integrand there does not depend on k. Panels geometric at 1.25 from the apex, the first one
  included (at 1.5 the error was 4.6e-4, from the 1/sqrt the integrand has at the apex); moments by
  the recurrence from a radian of phase and a Horner series below (36 terms summed per call took
  432 us a ray for 72 wavenumbers; 82 us now, against 880 us for 72 descent calls).
- **Second order is integrated over its level curves** (`LevelCurves`, `twice_diffracted_spectrum`):
  G(Lambda) does not depend on k, so one set of curves serves every sub-band. Rejected before it:
  a nested integral per wavenumber in a smooth window, whose taper needs 40 to 80 radians to be
  within 2e-3 to 4e-5 and so costs milliseconds a sub-band; steepest descent in two dimensions,
  which needs beta1 continued into complex z2 and risks crossing its poles; and a semi-cascade
  (the exact first order at the second edge, UTD after), which fails when the edges are near.
  Load-bearing in it:
  - **Chebyshev points of the first kind.** At u = 0 the curve is a point and G only a limit; a
    Lobatto node there returned 0, the interpolant never converged, and a pair took two minutes.
  - **u G is interpolated, not G.** Past a boundary point G goes as 1/sqrt(Lambda - B).
  - **The breaks include the boundary points** (`boundary_points`): where a directivity term's
    pole meets the real line on its edge's ridge, G steps. Without them the perpendicular pair
    was 1-11% off.
  - **A piece is judged against the whole transform's integral of |G|, estimated first from three
    points a segment.** Judged against its own segment, a small segment where G crosses zero was
    held to 1e-6 of almost nothing: one DigiCity pair took 3.9 s, 24 ms now.
  - **What is held still near a break, and how near two breaks are one**, is a hundredth of a
    radian at the largest wavenumber or a thousandth of the directivity peak's own scale
    (psi^2 r/4), whichever is smaller, floored at 1e-10 of the least path. A millionth of the span
    -- 4 mm on the long block's 4000 m edges -- smeared its grazing roof's transition and put the
    block 4.5 dB off at 1 kHz; a hundredth of a radian alone put it 0.62 dB rms off. Two breaks a
    floor-standing pair puts 1.5e-7 apart cost 26 pieces and 1.2 s until merged (5 ms).
  - **From the least path the pieces start geometric in u, a factor of eight apiece**, because
    what G does near it is on the scale of how near; one piece over the whole segment never
    sampled the small u where the long block's structure lies. Only the first segment needs it:
    applied at every break, the fine hold and the geometric pieces took DigiCity's four
    microphones from 61 s to 137 s on one thread for nothing.
  - **Directivities by angle addition** (`HalfAngle`): one sincos for each angle that varies along
    the curve, the four terms' half-angles from it; for parallel edges both are built once.
- **The plane factor is a half over a shared face.** Without it the integral over a block's roof
  is the cascade without UTD's grazing half, to 0.01-0.1 dB at k = 20: twice the field. With it the
  long block is within 0.10 dB rms of the BEM.
- **A twice-diffracted ray over edges that are not parallel takes its partners' corrections**
  (`twice_diffracted_partners`, `localized_twice_diffracted_*`). A boundary plane of one edge
  crosses the other, and there the second order carries half a once-diffracted field that only
  the cut-off first order cancels. In a window about the stationary point the raw integral drifted
  from -68.94 to -70.44 dB as the window grew from 40 to 160 radians on the perpendicular pair;
  with the corrections it holds -68.944 to -68.948 and is the corrected cascade's within 0.02 dB.
  The subtracted partner is the model's own ray for it (`PartnerModel`): with `--edge_source`
  the edge-source integral over the whole edge.
- **A point on the other edge's own line has no angle about it.** Where a roof edge meets the
  vertical corner below it, `where_along` skips the shared vertex and a partner's part must have
  length; without both the vertex read as a wave's boundary, and its stretch of no length gave
  every band of 20 of RWTH EDGE's source-receiver pairs a NaN. A factor that is not finite falls back to the
  cascade (zero for a ray past an end) as a last guard.
- **Vertical edges on the floor go through it** (`through_floor` in `point_paths.cpp`). The
  floor's images split one doubled edge between a ray and its mirror, and the tracer keeps the one
  above the floor; integrated over [0, L] each kept ray had a false end at the floor, a metre or
  so from every street-level Keller point, which moved DigiCity's microphone 7 by -1.5 dB where
  the edges' real ends move it by +0.6. `with the edge source a corner on the floor is one edge
  with its image` holds the kept pair to the doubled edge's field.
- **Once-diffracted rays past an end are kept** within `edge_source_reach` (100 radians at the
  lowest band) and faded between 50 and 100 radians of the end's excess, so the tracer's cut is
  no step; a point below the floor is its mirror's and is not. Twice-diffracted rays past an
  end are not kept: it would loosen the pair prune and cost a level-curve integral per kept
  ray, and none of the 30 steps above 1 dB in a 3771-position DigiCity sweep was one. Every
  step there was a partner the tracer has no edge, reflection or floor leg for.
- **Weights, not switches**: once-diffracted rays always take the integral, which matches UTD on
  long edges to p90 0.01 dB; twice-diffracted ones blend into the cascade in each wavenumber's
  transition parameter (0.3 to 1) and path excess (100 to 1000 radians).
- **Do not gate by a ray's share of the energy.** Leaving the integral to rays above 45 dB down
  moved DigiCity's microphone 7 by up to 0.48 dB: hundreds of twice-diffracted rays each under
  that share carry half a dB of end effects between them. (The near-twin effect above is the same
  lesson for dropping rays outright.)
- **The cascade's spreading for edges that are not parallel** is the second edge's own curvature
  along it (`twice_diffracted_ray`); its distance parameters stay the point-source ones, because
  the astigmatic one the curvature gives the second edge broke `the point-source field is
  reciprocal` by 0.04 dB.
- **CPU only.** `--edge_source` with `--gpu` is refused.
- **Where its time goes.** DigiCity's four microphones (17 bands, one thread) take 58 s against
  0.22 s for UTD, nearly all of it the level curves' integrand: two `atan2`, two `sinh` and the
  `asinh` of the two etas, and the directivities' half-angles. `std::nearbyint` in the
  directivity build and `std::hypot` in the geometry together cost 5% of the run.

Guard tests (`tests/test_edge_source.cpp`, and two in `test_point_paths.cpp`): the plane-wave
limits against Sommerfeld and the n = 3/2 series; the descent against the direct integral; UTD at
high k; reciprocity; continuity across a boundary and an end; the spectrum against the field; the
level curves against the nested integral on parallel, shared-face and perpendicular pairs; the
cascade against the localized integral at 0, 45 and 90 degrees; the long block through the field;
the floor's mirror pairs; a corner's adjacent edges.

### The canyon wave solver (`src/canyon_wave`)

A 2-D staggered pressure-velocity FDTD of a flush-roof canyon, run by reciprocity, for the urban
tables' low band (THEORY §6.6). What to keep:

- **The far field comes from the opening, and that needs the roofs rigid and flush.** The
  Neumann Green's function of z > H turns the opening's vertical velocity into the exact far field
  only because dp/dz vanishes on the rest of the plane. Unequal walls, or absorbing roofs, need
  another transform (a closed contour, or the half-space representation over the roofs too).
- **The recording is the vertical velocity on the face row at z = H**, the opening's columns
  only, and is divided by the source's own transform bin by bin (`V/Q`), so the pulse's shape and
  amplitude drop out. Bin 0 is left at zero: the pulse has no mean, so there is nothing to divide
  by. The transform's sign is exp(+i w t), the ground model's convention.
- **Normalization is checked over open ground**, not assumed: the same run and transform over a
  long recording line reproduce 4 cos^2(k z sin phi) within 0.1 dB from 45 to 135 degrees. That
  test is what would catch a lost factor of 2 in the Green's function or the source.
- **Walls absorb by default (0.05).** A rigid canyon rings for seconds in its cross modes
  (n c / 2W), which leak out of the opening only weakly; the bands were still moving by 0.8 dB
  between 1 and 4 s of a rigid run, and settle by 2 s at 0.05. The wall boundary is
  semi-implicit over the half cell (stable for any Z >= 0); an explicit v = p/Z there is not.
- **Runs stop at 30 dB of decay or 4 s.** Bands moved by under 0.01 dB beyond about -28 dB.
- **Cost** (M4 Pro): a 20 m canyon to 630 Hz at 10 points per wavelength is 275,000 cells and
  40,000 steps, 15 s on one thread and 8 s on eight. Threading per half step costs a thread
  start each time, so the table builder runs receivers in parallel instead. Each doubling of
  the top frequency costs eight times as much.
- **What follows the march is threaded too, and has to be.** Each recording point's transform
  (`spectra`) and each azimuth of a canyon's gains (`canyon_wave_gains`) go to the run's
  threads, the same bits at any count. Serially the gains alone are 11 s of a 20 m canyon's
  79 s to 707 Hz, and fourteen such canyons to 354 Hz take 65 s on 14 threads against 34 s
  threaded. With `--gpu` the march itself is the GPU's ([gpu.md](gpu.md)).
- **15 points per wavelength for tables**: 15 against 20 agree within 0.3 dB, 10 within 1.4.
- **The far-field sum steps its phase by recurrence** over the evenly spaced points: one complex
  multiply per point instead of a sine and a cosine. It is the table build's inner loop.
- **Oblique incidence is the reduced frequency f |u_perp|** (`canyon_wave_ratio`), exact for rigid
  walls; a reduced band under 5 Hz is open ground.
- **The builder's wave band** (`--wave_below`): each distinct canyon is solved once, canyons in
  parallel on one thread each, to the top of the last blended band. `blend_wave_band` mixes in
  energy and scales each facade's reflection by 1 - w: the wave solution holds the first-order
  reflections, so leaving them in would count them twice. Canyons with unequal walls are warned
  about and left geometric (the far-field transform needs flush roofs).
- Tests: `[canyon_wave]` (open ground against the image solution, independence of the grid's
  margins), and the hidden studies `[.canyon_wave_study]` (convergence, the geometric limit) and
  `[.canyon_tone_study]` (rotor tones against band tables, the stage-2 gate), whose results are
  in FIDELITY under Canyon wave solver.

### The full-wave reference (`src/fdtd`, `niceops_fdtd_field`)

A 3-D staggered pressure-velocity FDTD (THEORY §9.4) of a point source among rigid solids: the
reference the point-source model is held to on scenes no exact solution covers, such as EDGE's
scale models. It is not part of any prediction. What to keep:

- **Fourth order is the default, and its walls are mirrors.** The staggered 9/8 and -1/24
  differences reach two cells; beyond a closed face the pressure takes its mirror value and the
  velocity its negative (`pressure_difference`, `velocity_difference`, `fdtd_dp4`, `fdtd_dv4`).
  With that rule the two differences are exact negative transposes of each other on any
  arrangement of open and closed faces -- checked on random 1-D layouts, |G + D^T| = 0 -- so the
  march conserves the discrete energy and is stable to C = 6/(7 sqrt 3), whatever the solids.
  Do not "simplify" a ghost to zero or to the near value: that breaks the transpose.
- **The warp is what makes fourth order pay.** Its phase error is fast by nearly the same amount
  in every direction (0.45-0.61% at 8 cells a wavelength), so reading each frequency where the
  grid carries the physical wavenumber, averaged over directions (`grid_frequency`), leaves only
  the spread. Yee's error runs from 0 on the diagonals to its worst on the axes; warping it made
  the 30-wavelength test worse, so it is not warped.
- **Fourth order's sources and receivers are tricubic** over 4 by 4 by 4 cells, folded across
  the floor and the plane of symmetry, and trilinear where a solid is among the 64. Trilinear
  spreading alone cost second order 0.04-0.09 dB at 25 cells a wavelength, an error that grows as
  the square of the cell.
- **Voxels by the cell's centre, faces by flags.** One byte a cell: air, and whether each of its
  three low faces is open (air both sides) or a solid's (air one side). A box whose faces fall on
  multiples of the cell is exact, which is why the grid snaps its origin to a multiple of the cell;
  EDGE's LoD1 coordinates are whole millimetres, so 3 and 4 mm cells grid them exactly.
- **The PML is convolutional and lives in slabs** (Roden and Gedney, `Axis::slab`): one auxiliary
  field per axis and per face or cell, only within that axis's slab. A split-field PML, as the
  canyon solver has, would store three pressures in every cell of a 3-D grid.
- **The floor and the plane of symmetry are rigid faces, not PML.** `--mirror` halves the grid;
  a source on the plane injects half its strength, its image the other half.
- **Sources and receivers spread over their eight cells by trilinear weights, renormalised over
  the air among them.** The spreading filters what a point would radiate: the levels come out
  0.01-0.09 dB low at 25 cells a wavelength, rising with frequency (`tests/test_fdtd.cpp`).
- **Absorbing walls are the canyon solver's half-cell impedance** (stable for any Z), on the faces
  flagged as a solid's, from the air cell's side. The floor stays rigid. At fourth order two
  things go with it, and both are needed (`velocity_difference`, `wall_pressure`, `fdtd_dv4`,
  `fdtd_wall_p4`):
  - **Beyond a closed face the velocity is reflected about that face's own value**, 2 v_face - v,
    not negated. A rigid face's is zero, so nothing changes there; an absorbing one's is not, and
    negating gave its velocity 13/12 of its weight over the two cells' divergences. The wall then
    absorbed 7% too much (0.0537 for 0.05), a fraction of a dB in one reflection but more in a
    reverberant street.
  - **The face's update reads 25/24 of its cell's pressure less 1/24 of the next**, the transpose
    of those weights, so the march still dissipates for any impedance.
  A gated reflection at normal incidence measures it (FIDELITY); the unit test's direct-plus-
  reflected levels cannot see 7% of 0.05.
- **The air is applied after the march, by delay**: exp(-a c (t - t_peak)) on every sample is
  exp(-a L) on every arrival, which is how the ray model absorbs. Putting the loss in the march
  would need a frequency-dependent term per cell.
- **Truncation is judged band by band**: a band's level with and without the run's last tenth. The
  raw tail energy is no guide; it is mostly the pulse's content below the bands, which leaves the
  PML slowly.
- **Rigid canyons ring.** Between parallel rigid walls the image sources' energies fall only as
  1/n^2 and the tail as 1/n: in EDGE's UE13 the bands between the rows still moved 1.6-2.1 dB in
  the last tenth of a 36 ms run. MDF absorbs about 0.05 (EDGE's `mdf.json`); with that and 90 ms
  they settle.
- **The GPU march is four dispatches a step** (velocity, pressure, injection, recording), 256
  steps a batch, two batches in flight; `nk_u8` entered the kernel dialect for the flags. It
  matches the CPU's float march to 1e-4 of the peak. The injection covers every cell of the
  source's stencil, 64 for a tricubic one: dispatched for 32, it injected half the source of every
  fourth-order run without a plane of symmetry (which folds the stencil to 32 cells), 6 dB low.
  The GPU test runs both with the plane and without it for that reason. Throughput on the M4 Pro: 4.2-4.4e9
  cell-steps a second (600 M cells, 7158 steps in 1033 s). Memory is 17 bytes a cell.
- **Fourth order costs what second order does per cell-step on the GPU**: 4.2e9 either way on
  a 512 by 512 by 256 box. A step moves about 50 bytes a cell (the flags twice, p, and the three
  velocities read and written), so 4.2e9 is 210 GB/s against the M4 Pro's 273: the march is
  bound by memory, and the wider stencil's extra reads come from cache. Everything fourth order
  gains is therefore fewer cells and steps. On EDGE's UE00, order 2 at 3 mm took 505 s (floor
  source, 356 M cells) and 1453 s (drone, 900 M); order 4 at 6 mm takes 40 and 116 s, 12.5 times
  less, and is the closer of the two to the converged level ([FIDELITY](../FIDELITY.md)).
- **The CPU skips the general path where it can** (`FDTD_CLEAN`): a cell whose neighbours two
  deep along every axis are open air, outside every PML slab, takes the plain differences in a
  vectorisable loop over its run of the row (`clean_velocity`, `clean_pressure`), 1.4-1.6 times
  faster. The GPU kernels ignore the flag: being bound by memory, they have no arithmetic to save.
- **Each step marches only its light cones' box** (`FdtdSetup::cull`, `cone_boxes`): the cells
  the source's stencil can have reached by the step's end, and from which a receiver's stencil
  can still be reached by the last, both growing at the scheme's fastest group speed plus 16
  cells. A cell outside the first holds zeros it would only have added to; one outside the
  second holds values that can no longer arrive. Three things to keep:
  - **The speed is the dispersion relation's, not c** (`fdtd_fastest_group_speed`). Fourth
    order's leapfrog outruns sound: 1.01 c along an axis at 7.5 cells a wavelength, and 1.23 c
    near the Nyquist wavenumber on a diagonal. The pulse carries nothing up there, but a frozen
    cell is a sharp edge in a full field, which excites every wavenumber, so the receivers' cone
    has to shrink at the fastest. Second order's is c, at k -> 0.
  - **The margin is the dispersive front's tail.** Against full marches (a block and a wall
    along a strip, 310 and 620 cells long), no margin left 1e-5 (order 4) to 2e-4 (order 2) of
    each receiver's peak, and 4 cells reached the float floor of a few 1e-7 (1e-6 on the GPU) at
    either order and either length; band levels moved by under 1e-5 dB. 16 is four times that,
    and costs 2% of the cell-steps.
  - **What it saves depends on the run's length.** It cuts only the first and last domain
    crossings, so a run of the default length (twice the box's diagonal) on UE00 marches 70% of
    its cell-steps, 28.5 and 81 s instead of 40 and 116, its levels within 5e-6 dB of the full
    march's. UE13's 90 ms canyon run would keep about 89%, by the boxes alone.
  `culling the march to its light cones` compares culled and full marches at 4e-6 of each peak,
  which no margin, or a receivers' cone shrinking at half the speed, fails; `the FDTD's fastest
  waves` pins the speed itself, which a margin this wide would otherwise hide.
- **Rejected: marching the GPU a column to a thread.** Carrying a column's pressures in
  registers, so each is read once, ran at 3.0e9 cell-steps a second with the columns in rows and
  2.66e9 in 16 by 16 tiles, against the cell kernel's 4.2e9; a fast path for clean cells added
  3%. With the cell kernel already at the memory's bandwidth there is nothing for reuse to save,
  and a column's serial walk hides less latency.
- **Rejected: stopping a run when its recordings fall quiet.** A receiver deep in a shadow can
  still be gaining 1e-5 of the emitted energy, and some arrivals are late (UE00's ends diffract
  at 17 ms), so any threshold is a guess at what is still to come. `truncation_db` reports what
  the last tenth moved instead.
- **Rejected: a thinner PML.** With receivers 5 cells inside it, 8 to 20 cells all leave 5-6e-3
  of the peak (0.07-0.13 dB in a band) and 6 cells three times that; the residue is the low
  frequencies, for which 20 cells is less than a wavelength, not the depth. Ten cells would save
  11% of UE00's cells for nearly twice the worst band's error.
- **Rejected**: a 2.5-D boundary-element solve (a 2-D BEM per wavenumber along the rows, then the
  transform). It is exact for rows of constant section but EDGE's stepped rows (UE06/07) vary
  along them, and a canyon's near-trapped modes need a fine wavenumber grid. A scalar 7-point
  scheme halves the memory but has no simple PML in second-order form, and a sponge needs a
  wavelength of depth.
- **Keep test geometry on the grid.** A face between cell faces is moved to the nearest one: a
  concave corner with both planes 4 mm off the grid, at 8 mm cells, was 0.56 dB out at 500-800 Hz,
  and 0.01 dB once they were on it.
- Tests: `[fdtd]` -- the two rays over a floor (and the plane of symmetry, and the air applied by
  delay), the exact wedge (`edge_source_field` over the whole line plus geometrical optics), a
  concave corner against its three images, an impedance wall's plane-wave reflection, a point a
  fraction of a cell inside a solid refused, fourth order's phase over twenty wavelengths, the
  culled march against the full one, the fastest group speed, and the GPU against the CPU.

`tools/profile_scattering.py` gives a relief's scattering coefficient for `--scattering_spectrum`
from the same 2-D boundary elements as `thick_barrier_bem.py`: two closed boxes of one outline, the
relief on one, CHIEF points inside both, the correlation of their far fields over the front
half-plane. The boxes' edges scatter alike and leave s alone.

## Weighting

`src/weighting.{hpp,cpp}` holds A-weighting and band-limited energy integration.

- **The A curve is panam's `dBAw`**, the IEC 61672 pole-zero form rather than a table, so it
  evaluates at whatever band centers a database carries. It is pinned in
  `tests/test_weighting.cpp`. In this form the curve is **not** identically 0 dB at 1 kHz but
  −1.49e-6 dB; a test asserting exact 0 there will fail.
- **`a_weighting_gains()` folds the weighting and the band limit into one per-band multiply**, with
  excluded bands exactly zero.
- **`weighted_level()` returns −∞ for silence, not 0 dB**, so an excluded spectrum cannot survive an
  energy sum as a quiet-but-real level.
- `c_weighting()` / `c_weighting_gains()` are the C curve for `--csel_file`, from the same
  constants; see [metrics.md](metrics.md#time-weighting-day-metrics-c-sel).

## Faddeeva function

`src/faddeeva.{hpp,cpp}` computes w(z) = exp(−z²) erfc(−iz). It exists for the ground model's
boundary-loss term, which panam writes as `exp(-w)(1 - erf(-i sqrt w))` using SciPy's complex
`erf`. With z = √w that *is* w(z), and evaluating it as written overflows. The UTD transition
function's table (`wedge.cpp`) is fitted to it, on the line arg z = 45°, the first time a wedge
coefficient is wanted.

- **Three methods, chosen by the argument** (`faddeeva()`):
  - |z| ≥ 6 with Im z ≥ 0, or x² − y² ≥ 36: a 12-term asymptotic series, 1e-11 relative.
  - |z| < 6 with Re z ≥ 0 and −0.25 ≤ Im z ≤ 6: a sixth-order Taylor series about the nearest node
    of a 0.05 grid (15,246 nodes, filled once at start-up from the Weideman series), 3.6e-12
    relative, losing digits further below the real axis.
  - Everything else: Weideman's 32-term rational approximation (`faddeeva_reference` is that series
    everywhere).

  The domains were sampled from real runs (480,000 arguments from point and plate receivers on the
  demo case) and each method measured against SciPy's `wofz` on them.
- **The Weideman series runs as four interleaved Horner chains over t⁴** rather than one over t.
  31 dependent complex multiply-adds is a latency chain, not a throughput problem, and it ran at
  about a quarter of the machine's floating-point rate. Splitting by the power of t modulo four
  gives chains a quarter as deep that issue together (the gain is in
  [performance.md](performance.md#the-spectral-and-ground-paths)). The coefficient table is left
  exactly as Weideman's construction produced it and indexed by
  stride, so it stays checkable against the source.
- **It is written in real arithmetic.** Every `std::complex` multiply carries a branch to recover
  NaN/infinity cases that cannot arise here. The two trailing divides are folded into one:
  `2p/d^2 + c/d == (2p + c d)/d^2`.
- **Accuracy.** Max relative error 3e-13 over the domain the ground model reached while the
  boundary loss was cut at |w| = 500 (Re z ∈ [0.09, 22.3], Im z ∈ [−1.2, 15.6]), against `wofz`;
  past it the tests hold it to `wofz` at 1e-11 out to |z| = 700.
- **The lower-half-plane reflection `w(z) = 2exp(-z²) - w(-z)` is safe** because a passive ground
  keeps z within 45 degrees of the real axis wherever it is below it, so exp(−z²) is at most 1 in
  size.

## Ground effect

`src/ground_effect.{hpp,cpp}` is the excess ground attenuation, ported from panam's `ega()` and
pinned against it to 1e-9 dB in `tests/test_ground_effect.cpp`. The impedance, turbulence and
roughness models follow Attenborough & Taherzadeh (2026), as panam's `ground_plane` module does,
and are pinned against `ground_plane.pole_level` there (1e-6 dB). A change to any of them is a
change to panam first.

- **It returns a gain, not a loss**: +6 dB over rigid ground with the receiver on the surface
  (and no turbulence).
- **`ground_impedance()` is hoisted out of the observer loop**, and the allocation-free overload is
  what the loop calls. Returning a vector there allocated 94 million times and cost 20% of the run.
- **The band loop in `excess_ground_attenuation` is real arithmetic**, with its band-invariant
  branches (boundary loss on/off, the sinc, turbulence) hoisted so the body stays straight-line.
  Bands outside the metric's limits are not computed at all: `first_band`/`last_band` bound a
  contiguous run, which is always what the limits select because band centers are sorted. The
  details, and what not to undo, are in [performance.md](performance.md).
- **The boundary-loss term has no cutoff.** It is evaluated by the Faddeeva function at any |w|,
  as panam's `ega()` does, and pinned to panam b758e3a (`ega` 1e-9 dB, `pole_level` and
  `board_level` 1e-6). The cut at |w| = 500 it replaced is in
  [HISTORY](../HISTORY.md#changes-that-affect-results). Two things to keep:
  - **The gain bound past the edge is analytic** (`reflection_coefficient_bounds`). Past
    |w| = 500, |Q| ≤ |R_p| + |1 − R_p| (1.003/(2|w|) + |P|), with P = 2i√π z exp(−z²) the surface
    wave, whose size peaks at |z|² = 1/(2 cos 2θ), θ the parameter's direction, which does not
    depend on range. For soils that peak is far inside the edge. For a nearly lossless stiff
    surface (a thin layer's band, θ a hair under −45°) the wave carries on: |Q| = 151 at
    |w| = 1250 at grazing. The sweep stops at the edge and the formula takes over; `a surface wave
    carried past |w| = 500 ...` holds it to a sweep out to κ = 1e10.
  - **A mass-like admittance (Im β > 0) keeps NICEOPS's branch, not panam's.** Its parameter is
    in the second quadrant. panam squares it and takes the principal root, the negative, whose w()
    grows as exp(|Re w|) and overflows by |w| ~ 1000. Such a surface has no surface wave, and the
    parameter itself gives none. A thin layer gives this in some bands.

### The ground table (`GroundTable`)

By default the point receiver interpolates a per-trajectory table instead of evaluating the model
at every observer and state; `--exact_ground` bypasses it.

- **Nodes in s = √(h_s/r) and ln r**, r = hypot(h_s, d): 641 angle nodes and 40 per decade by
  default. s resolves grazing incidence, where the model changes fastest. Values are exact at the
  nodes and interpolated bilinearly, as a power ratio, in between. Ranges outside [r_min, r_max]
  are clamped to the table's edge.
- **The interpolant is a convex combination of exact values**, so `ground_gain_bounds()` bounds it
  as it bounds the model, and the peak search and `--db_cutoff` stay exact.
- **Why it exists:** the model costs about 20 ns per band and the observer loop reaches it 94
  million times on the demo case, where it was 5.5 of 7.1 s. A lookup is four reads per band. The
  tables are built in well under 0.1 s.
- **The table holds one receiver height.** A path with its own receiver height (terrain, below)
  evaluates the model exactly.

## Plate receiver (`--plate_table`)

`src/plate_receiver.{hpp,cpp}` is a flush ground-plane microphone on a GR1425 plate, from panam's
axisymmetric boundary-element tables (`axisymmetric_bem.write_netcdf()`): p/p_free = P_d + Q P_r,
with P_d and P_r tabulated over band, sub-frequency, elevation and azimuth, and Q the ground's
spherical-wave coefficient at the geometry. The file comment in `plate_receiver.hpp` has the netCDF
layout.

- **The table records its ground, and `load()`/`bound()` refuse a mismatch** with the run's ground
  model and band centers. P_d and P_r depend on the ground.
- **Q excludes roughness and turbulence** (`spherical_reflection_coefficient`): the plate's image
  path is millimeters longer than the direct one, so neither decorrelates it.
- **P_d and P_r are interpolated bilinearly** (real and imaginary parts) in elevation and azimuth,
  held below the lowest and above the highest elevation, with azimuth periodic.
- **Three evaluations exist.** `level_exact()` takes Q at every sub-frequency, as panam's
  `board_level` does, and matches it to 1e-6 dB in `tests/test_plate_receiver.cpp`. `level()` and
  `power_ratio()`, which the run calls, take Q quadratic in log frequency through each band's
  edges and center (2n + 1 evaluations for n bands, not 5n), about 2.5x cheaper for the whole
  run. Its error against `level_exact` is in
  [FIDELITY](../FIDELITY.md#ground-and-plate-tables).
- **`tabulate()` tables those Q values** over the image path's grazing cosine (as √cos) and log
  length, 161 angles and 20 per decade, coarser than `GroundTable` because Q is smooth and the
  plate's interference is in P_d and P_r. `level_exact()` never uses it.
- **Orientation.** The tables put the microphone on +y of the plate. `plate_frame_azimuth` turns
  the propagation direction by the microphone's offset direction, and it is the one function both
  the straight path (`Database::spectral_level_along`) and the refracted path use.
  `--plate_mic_azimuth` defaults to `outboard` (`outboard_azimuth`: away from the track's nearest
  state, or to the track's right directly under it).
- **`PlateReceiver::gain_bounds`** bounds the plate for the peak search and `--db_cutoff`;
  `plate_windowed_is_exact` checks that skipped states change nothing.

The acoustic loop's cost with each receiver, tabled and exact, is in
[performance.md](performance.md#the-spectral-and-ground-paths).

## Terrain (`--terrain`)

`src/terrain.{hpp,cpp}`, and `terrain_ray` in `src/propagation.hpp`. A bare-earth grid (netCDF,
heights in meters; the validation harness's `terrain/dem_to_netcdf.py` makes one from USGS 3DEP).
Grid observers stand on it, and each pair's ground geometry is taken over the least-squares mean
ground line along its ground projection, observer to the source's footprint (CNOSSOS-EU's mean
plane, NORAH2's practice over terrain). Spreading, absorption and the sphere lookup
keep the true 3-D path. The geometry is in [THEORY.md](../THEORY.md) §6.3.

### Things to preserve

- **Level terrain is the flat run, byte for byte** (`terrain_flat_*_is_flat`: point, plate and
  scattered). That rests on three spellings:
  - `Terrain::elevation` and `TerrainRaster::height` interpolate from a corner (`a + f (b - a)`),
    because the product form's four weights do not sum to exactly 1 and a constant grid came out
    1.9e-13 ft off;
  - `plane_geometry` divides by √(1 + b²) rather than multiplying by its reciprocal, so a level
    line divides by exactly 1;
  - the test grid is written as `-((-1250 * 0.3048) - (54.2651 * 0.3048))`, the height
    `--terrain_offset -1250` and `--ground_elevation 54.2651` take back to zero in the program's
    own arithmetic.
- **The ground's z is H + `--terrain_offset` − `--ground_elevation`**, the same
  `--ground_elevation` that is a WGS84 track's z = 0. NICEOPS converts no datums. Leaving out
  either is a plausible wrong answer, so the run checks what it can (the list is in the manual's
  [Terrain](../MANUAL.md#terrain)). The AC track's own lat/lon check uses `Trajectory::recorded`,
  kept as read, before any trim. Every refusal has a CTest on its message.
- **A point receiver's height above each line is its own** (`RayPath::receiver_height`, NaN for
  every other path), so `spectral_level_impl` evaluates the ground model exactly when it differs
  from the model's, *before* the `GroundTable` branch, which holds one height. Do not read the
  table at the pole's height instead: a height error of δ moves the path difference by 2δ sin
  (grazing), a radian at 1 kHz for 0.3 m at 30°. A plate lies along the line under its microphone
  and keeps its table (source height h_s − h_r). `a path's own receiver height is the ground
  model's, past its table` holds the branch.
- **Heights below the line are raised to zero**, CNOSSOS-EU's rule: a receiver in a dip, or a
  source behind rising ground, which is screened only with `--terrain_screening` (below). The
  run counts them and warns.
- **A negative distance along the line is taken by its size, not clamped to zero.** A nearly
  overhead source over ground falling toward it has its foot behind the microphone's. The image
  construction depends on d², so |d| keeps the model's direct path equal to the real range (zero
  shortened it), and `terrain_ray` turns `separation` round so a plate is crossed from the side the
  arrival actually comes from. `a source whose foot lies behind the microphone's keeps the real
  range` holds both.
- **`terrain_hills_point_golden` and `terrain_hills_plate_golden` pin the hills runs** at 1e-6
  (`golden/terrain_hills_*.csv`, from the C++ build). The windowed/marched pair compares one binary
  with itself, so a call-site mutation that handed the plate's microphone to `terrain_ray` as a
  point receiver moved cells by 0.01 dB and passed every other terrain test. Regenerate them only
  with a deliberate change, and say so.
- **The bounds need nothing new** without screening (with it, see
  [Screening](#screening---terrain_screening)). `ground_gain_bounds` and
  `PlateReceiver::gain_bounds` hold for every geometry and receiver height, and the bound test's geometry is the true 3-D path, as the
  level's spreading is. `terrain_hills_point_windowed_is_exact` and
  `terrain_hills_plate_windowed_is_exact` hold it. The tables' reach (`reach()` in
  `process_trajectory`) extends below z = 0 under terrain; node placement does not depend on the
  reach, so the flat case is unmoved.
- **Refused with `--refraction`.** The ray models reflect off and shadow onto z = 0; a tilted
  plane per path handed to them would be another model. That refusal is the only one; the one
  with `--propagation_table` had no CTest, and `terrain_table_refuses_refraction` holds that
  `--refraction` still wins when both are given.
- **A propagation table reads the same path.** `straight_level` finds each pair's `terrain_ray`
  path once, when a ground model or a table will read it, and hands it to `table_level` and then,
  for a pair the table misses, to the ground model. The plane counts as taken, and a raised
  height as raised, only for a pair one of them read: a pair the table misses with no ground
  model is free field. Level terrain stays byte-identical through the same three spellings
  (`terrain_flat_table_*_is_flat`).
  - **The bearing is the true separation's; the plate's azimuth is the turned one's.** Read at
    |d′| on a fixed bearing, an anisotropic table is continuous as d′ passes through zero; the
    turned separation's bearing jumps 180° there. The plate is still crossed from the side the
    arrival comes from. `a table read over terrain is continuous through d' = 0` holds both.
  - **A point receiver outside the table's rows falls back and is counted, never clamped**, for
    the reason in the bullet above: a height error moves the path difference. The run warns with
    the heights seen, and up front for a table with one receiver row.
  - **The bound needs nothing new.** The table term is its largest node (times the largest K
    node for plates), which no interpolation exceeds, whichever cell a pair reads, and the
    pairs that fall back are covered by the ground bounds above
    (`terrain_hills_table_*_windowed_is_exact`, the bump table among them).
- **Observers' ground is read from the grid itself in `run()`**, once, not from the raster, which
  is a resampling of it.

### Screening (`--terrain_screening`)

`src/screening.{hpp,cpp}` holds the CNOSSOS-EU kernel and the knife edge,
`GroundPlaneCache::may_block`, `blocked` and `fit_range` their terrain, and `screened_level` in
`process_trajectory` the CNOSSOS call (the knife edge's is below). The equations are
in [THEORY §6.4](../THEORY.md#64-terrain-screening). What to keep:

- **An unblocked pair takes the untouched call.** `straight_level` asks `blocked()` first and
  otherwise goes on exactly as without screening, so level ground, a uniform slope and the hills
  (where nothing blocks) are the unscreened runs byte for byte
  (`terrain_flat_*_screened_is_flat`, `terrain_tilted_*_is_unscreened`,
  `terrain_hills_unblocked_is_unscreened`, each with a log check that nothing was screened). A
  pair with no ground model and no table is never given a `terrain_ray` path: it stays
  `spectral_level_at_observer`'s.
- **The onset is strict.** A sample exactly on the line of sight does not block. CNOSSOS's own
  rule screens down to δ = −λ/20, but over level ground a sample near a 1.2 m microphone has
  δ ≈ h²/2s, about 0.14 m at 5 m, under λ/20 at 63 Hz, so it would screen level terrain. (The
  Directive also asks δ > λ/4 − δ′, δ′ by the ground images, which keeps the low bands out; the
  conformance case TC06, a path clearing a plateau's edge by 0.67 m, is diffracted at 500 Hz
  and 1 kHz in the TR and unscreened here, `ISO/TR 17534-4 TC06` in `tests/test_screening.cpp`.)
- **A minimum protrusion** (`--screening_min_height`, 0.5 m): a sample blocks only when it is
  that far above the line of sight. At Sweetwater the 1 m lidar grid's texture 3–5 m from a
  flush board, 0.03–0.14 m above a path at 1–2° elevation, screened it, and the screened level came
  out up to 1.1 dB above the unscreened one: with δ near zero the 4.77 dB floor and the two side
  ground terms beat the whole path's plate term. The threshold leaves the quick test (a necessary
  condition) as it was; it moved the ridge goldens by up to 0.9 dB (points) and 1.6 dB (plates)
  in about 90 of 441 cells, all at the shadow's edge, where it was regenerated.
- **One spelling of the slope** (`sight_slope`, (z − h_m)/s) for the horizon, the line of sight
  and the scan. The quick test reads each ray's horizon at the last sample strictly inside the
  path (ceil(D/Δ) − 1, never the footprint's), and with no turn reads ray a only, as `line_from`
  does.
- **The profile is one per pair**, blended in bearing with `line_from`'s own turn fraction, and
  the microphone is at its height above the raster's ground at its foot, the base of the running
  sums, so a flush microphone is never under its own profile. An observer file's z kept without
  `--observers_on_terrain` is where the microphone is evaluated, so the horizon is taken from it
  (`set_observer_at`), held no lower than the ground at the foot; taking the height above the
  raster instead screened such a microphone from where it is not (`terrain_ridge_screened_file_z`:
  3000 ft up, nothing is screened). A source under the profile at its
  footprint (states down to 10 m under the terrain are let through) is raised onto it.
- **No ground for the whole path on a screened pair** (`ground = nullptr`, `table = nullptr`):
  the ground enters only through the two side terms, so it counts once. A propagation table's
  level is not read either, since its ground is the whole path's; that is why the combination
  needs a ground model (refused otherwise) and the log counts the pairs that left the table.
- **F is clamped to 1** for the bound, and the run counts how often it binds. The bound's term is
  `screen_gain_bound(G_S, G_R)` = max(1, G_S) max(1, G_R)/3, with G_S over every class in use,
  a plate's included (the source side is always the point model), and G_R the point and plate
  bounds. Without a ground model it is 1/3 and changes nothing. The 3 is `cnossos_floor`,
  written once. The run logs how many bands the term raises above the ground term: the
  windowed-is-exact tests pass without it, so that line is what `terrain_ridge_screened_bound`
  checks.
- **Side lines come from the running sums** (`fit_range`, re-based to k0), so a screened pair
  fits nothing sample by sample; raw heights and the horizon cost about 16 bytes per sample per
  ray per thread, only with the option.
- **Rejected:** Deygout on the resampled profile (every sample of a rounded crest is a hull
  vertex, so a path just grazing it took about 13.6 dB, depending on `--terrain_step`); blending
  the two rays' losses in dB (a pair partly screened by a profile not its own, and a duplicate of
  the four-way ground selection on half-blocked pairs); a persistent per-ray hull (about 20 bytes
  per sample, and only useful with the dB blend); the sphere looked up toward the first edge
  (breaks the direction bounds' bins); a fade-in below δ = 0 (breaks the level-terrain identity).
- Tests: `terrain_ridge_screened_{point,plate,noground,classes,table}_windowed_is_exact` behind
  the ridge fixture (`integrated/terrain_ridge.nc`, 25 m high, 250 m to the track's side), each
  with a log check that pairs were screened; the point and plate runs pinned at 1e-6
  (`terrain_ridge_screened_{point,plate}_golden`); the refusals; and the Catch2 `[screening]` cases,
  including the bound sweep and the half-plane comparison.

#### The knife edge (`--terrain_screening kirchhoff`)

`knife_edge`, `knife_edge_loss` and the band-loss table in `src/screening.cpp`. In
`straight_level` a blocked pair computes its loss (`knife`) and goes on down the unscreened call,
which takes it off whichever level that reaches. What to keep:

- **The same blocked pairs as CNOSSOS.** It reuses `blocked()`, `--screening_min_height` and
  `rubber_band`, so where nothing is blocked the runs are the unscreened ones byte for byte
  (`terrain_flat_*_kirchhoff_is_flat`, `terrain_tilted_plate_kirchhoff_is_unscreened`,
  `terrain_hills_kirchhoff_is_unscreened`, each with a log check that nothing was screened).
- **The loss goes on every unscreened level**: the ground model over the mean plane
  (`spectral_level_along`'s extra loss), the table (`table_level`'s `added_loss`), and with
  neither the free field along the straight ray, which replaces `spectral_level_at_observer` for
  a blocked pair. Each call site has its golden; the table's and the free field's were added
  because dropping the loss there passed every windowed/marched pair.
- **Keep the whole path's ground, not CNOSSOS's two side terms.** Behind the Sweetwater hill
  `cnossos` came 29 to 48 dB below the unscreened level per blocked state from 500 Hz to 2 kHz
  (median over 11 runs), past its 25 dB cap on Δdif, where the measured shielding was 4–7 dB: the
  rest is its side terms in place of the path's ground.
- **No bound term.** |F(v)|² ≤ 1/4 for v ≥ 0 and falls monotonically, so the band loss is at
  least 10 lg 4 and a blocked pair is below its unscreened level, which the bounds already cover.
  `screen_gain_bound` is skipped, and behind the ridge the run evaluates the unscreened run's
  fraction (98.44%, against 99.26% for `cnossos`). A lit-side extension (v < 0, where |F|² reaches
  1.37) would need a term.
- **The edge** is the higher of the two horizon rays where they cross the aperture plane, which
  stands at the sample deepest above the line of sight (the first, on a tie). δ is summed as
  H²/(√(d² + H²) + d), which keeps its digits for H ≪ d. `the equivalent edge is the lowest point
  both ends see` checks it against a bisection on a peak, a hill and 200 random profiles.
- **The band loss is a table in ln x**, x = δ f_c/c, which is all it depends on: built on first
  use (a function-local static, so thread-safe), every 1/32 from 10⁻⁸ to 10⁶ (1033 nodes, about
  2 ms), with each node's exact slope for cubic Hermite. It is within 10⁻⁹ dB of the quadrature
  for x from 10⁻⁶ to 10³ and 6·10⁻⁸ dB at its top end; beyond its ends the quadrature is
  evaluated directly.
- **Rejected:** the eight-point quadrature per band at run time. Eight Fresnel integrals per band,
  each a series or a continued fraction, came to 50–64 µs per screened evaluation for δ of 0.1 to
  1 m, and the full ridge grid took 74 s against 8 s unscreened; the table takes 0.12 µs for 31
  bands.
- Tests: `terrain_ridge_kirchhoff_{point,plate,noground,table}_windowed_is_exact`, each with a log
  check that pairs were screened (the table's that they kept the table); all four runs pinned
  (`terrain_ridge_kirchhoff_*_golden`, the table's at 1e-5); `terrain_ridge_kirchhoff_is_not_cnossos`;
  and the Catch2 `[kirchhoff]` cases: the Fresnel integrals against SciPy, |F|² at grazing and in
  the shadow, the half-plane comparison, the band loss's floor and order, the table against the
  quadrature and the quadrature against Simpson, and the edge against bisection.
- **Cost** behind the ridge on the full demo grid (point receivers on variable-porosity ground of
  200 kPa·s/m², 8 threads, one run each on a shared M4 Pro at load 5 to 10, so read the ratios):
  59.1 s of CPU unscreened, 65.7 s `kirchhoff`, 77.5 s `cnossos`; peak RSS 239, 250 and 249 MiB.
  The knife edge's 6.6 s over 8.92 million screened evaluations is about 0.7 µs each, mostly the
  blocked test, the rubber band and the edge.

### How the plane is found (`GroundPlaneCache`)

A direct fit per pair is about 300 samples for a 3 km path, a microsecond, 94 million times on the
demo grid. Instead:

1. **Per trajectory**, the grid is resampled onto the local frame every `--terrain_step` (default:
   the grid's own spacing, at least 5 m) as a `TerrainRaster`: an inverse projection per node, not
   per sample. Its margin is two steps plus the diagonal times the bearing step, which is how far a
   bracketing ray reaches.
2. **Per observer**, profiles run along 1440 fixed bearings (0.25°, 13 m apart at 3 km), each
   sampled every step, kept as running sums of z and k·z. The least-squares line over any whole
   number of steps then costs a few flops (`line_from_sums`). A bearing is sampled lazily, only when
   a pair first needs it and only as far as that pair needs. Heights are summed relative to the
   observer's, or a 1000 m site cancels in `sum_kz - n/2 sum_z`.
3. **Per pair**, the line is interpolated linearly in bearing and in length between the four fits
   around it, which keeps it continuous. A path shorter than one step takes the first step's fit
   along its bearing, not a one-sample "fit" with zero slope, which was wrong for every nearly
   overhead state on a slope.

Along a bearing, at a whole number of steps, this is the direct least-squares fit
(`tests/test_terrain.cpp`, 1e-9 in slope). On a uniform plane it is exact up to the cosine's
curvature between bearings: a 0.05 slope errs by 1e-7, or 0.4 mm at 3 km
(2e-7 in slope at worst along a ray at a whole step).
On rolling synthetic terrain with ±15 m hills it stays within 0.15 m of the direct fit at both ends
of every path tried. The mutation that read one bearing twice fails two cases in
`tests/test_terrain.cpp`.

Cost on the full demo grid (S-76D track, B407 database, 14,400 observers × 6,551 states, 14
threads, `--humidity 20`), CPU time:

| | level ground (no `--terrain`) | level `--terrain` | hilly `--terrain` |
| --- | --- | --- | --- |
| `--ground_effect`, 1.2 m point receivers | 32 s | 35 s | 127 s |
| `--plate_table` | 71 s | 75 s | 75 s |

The planes cost about 3 s of CPU, roughly 35 ns a pair. Point receivers over real terrain cost 4x,
because each pair's receiver height differs from the one height the `GroundTable` was built for,
so the model is evaluated exactly. On level terrain the height is the pole's and the table is read.
Plates keep their table either way. For the 2017 arrays (about 40 microphones) none of this is
noticeable.

With propagation tables and screening, on the same grid at 8 threads on an M4 Pro (asserts on,
variable-porosity ground of 200 kPa·s/m², best of three; the machine was shared, load 10 to 31,
so read the ratios rather than the seconds). The table is the five-row ground-model test table
(`tables/ground_model_terrain_test.nc`: rows 0, 1.2 m, 8, 16 and 25 ft, source heights and ranges
every 20 ft):

| run | CPU | acoustic loop | evaluated | peak RSS |
| --- | ---: | ---: | ---: | ---: |
| level ground (no `--terrain`), 1.2 m point receivers | 23.2 s | 2.90 s | 98.48% | 231 MiB |
| hills, point receivers | 93.4 s | 11.8 s | 98.90% | 240 MiB |
| hills, point receivers, `--propagation_table` | 34.2 s | 4.32 s | 98.89% | 311 MiB |
| hills, `--plate_table` | 61.3 s | 7.72 s | 98.91% | 275 MiB |
| hills, `--plate_table --propagation_table` | 35.4 s | 4.46 s | 98.92% | 343 MiB |
| hills, point receivers, `--terrain_screening cnossos` | 98.4 s | 12.4 s | 99.63% | 249 MiB |
| ridge, point receivers | 61.8 s | 7.84 s | 98.44% | 240 MiB |
| ridge, point receivers, `--terrain_screening cnossos` | 80.0 s | 10.2 s | 99.26% | 250 MiB |

- **A table over terrain is cheaper than the ground model it replaces**: the point receivers'
  exact evaluation (4.0x the level run) becomes a lookup at their height above the plane (1.5x).
  The table covered 99.999% of the evaluations; the rest took the ground model. The extra 70 MiB
  is the table itself. Its accuracy at this row spacing is in
  [FIDELITY](../FIDELITY.md#propagation-tables).
- **Screening with nothing to screen** (the hills, byte-identical to the unscreened run) cost
  5% in CPU: the quick test per pair, and the bound's screened term, up to 4.1 dB above the
  ground term in 28 of 31 bands, which evaluated 0.7 points more of the pairs.
- **Behind the ridge** 9.27 million evaluations were screened (9.9% of those made, at 6057 of
  14,400 observers; 8.95 million at 6048 with the minimum protrusion added since), 3.51 million of them (38%) reached the 25 dB cap in some band, and F was
  held at 1 on neither side. The 18 s of extra CPU, the pairs the looser bound adds included, is
  about 2 µs per screened evaluation, as at one thread on the 441-observer grid (0.55 s for
  289,159).
- **The screening cache** holds about 16 bytes more per sample per ray. At these fixtures' 42.8 m
  step a path of 4.4 km is about 100 samples, about 2.3 MB per thread with every ray sampled:
  the peak footprint rose 3.4 MiB at one thread and 7.5 MiB at eight (441 observers), and 7.1
  MiB at eight on the full grid. At a 5 m step and 3 km, 600 samples, it is about 14 MB per
  thread.

### Terrain tests

- *Level terrain reproduces the run without it, byte for byte*: point receivers, plates, and
  scattered observers stood on it (`terrain_flat_*_is_flat`).
- *A uniformly tilted plane gives the closed-form h_s, h_r and d*, nearly overhead included, for
  point receivers and plates, and the rest of the path record stays the straight line's. The
  plane's own geometry is a rotation, hand-checked for rising and falling ground
  (`tests/test_terrain.cpp`).
- *The hills runs are pinned* (`terrain_hills_point_golden`, `terrain_hills_plate_golden`), and
  the windowed runs over hills equal the march (`terrain_hills_*_windowed_is_exact`).
- *The terrain reaches the ground model, not just the observers.* Microphones stood on the hills
  through an observer file, without `--terrain`, give different levels, and the same pair with no
  ground model agrees (`terrain_reaches_the_ground_model`,
  `terrain_without_ground_places_observers`). A level that reads the `GroundTable` at the wrong
  height is caught in `tests/test_spectral_propagation.cpp`.
- *A propagation table over terrain* is the flat read on level terrain, byte for byte
  (`terrain_flat_table_*_is_flat`), keeps the windowed search exact over the hills
  (`terrain_hills_table_*_windowed_is_exact`), counts its pairs toward the raised-height warning
  with no ground model (`terrain_hills_table_warns_unscreened`), warns of receiver rows
  (`terrain_table_warns_one_row`, `terrain_table_warns_receiver_rows`), and warns of a steep
  plane for a refracting table but not for a ground-model one
  (`terrain_steep_refracting_table_warns`, `terrain_steep_ground_model_table_is_quiet`, on the
  5.7° plane `integrated/terrain_steep.nc`). The geometry is checked in `tests/test_terrain.cpp`
  (`a table read ...`). A ground-model table over the hills is the terrain ground model to its
  interpolation (`terrain_hills_table_point_is_the_ground_model`, from 20 to 250 Hz over rows
  every half foot, `tables/ground_model_terrain_low_test.nc`; and
  `terrain_hills_table_plate_is_the_plate_model`); the measured error is in
  [FIDELITY](../FIDELITY.md#propagation-tables).
- *Every refusal has a CTest* that passes on its message.

What terrain changes at Amedee, and the slope sensitivity, are in
[FIDELITY](../FIDELITY.md#terrain-at-amedee).
