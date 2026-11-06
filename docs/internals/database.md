# The sphere database and condition blending

`sphere_database`, `predicates`, `delaunay`, `interpolator` and `condition_selection`, along with
legacy-database repair, the broadband and spectral paths, and where the databases come from. The
equations of blending and of reading a sphere are in [THEORY §2](../THEORY.md); how a user
chooses among them is in the manual's [Noise databases](../MANUAL.md#noise-databases).

## Modules

### `src/sphere_database.{hpp,cpp}`

The `Database` class. It loads a netCDF noise database (`.nod`): one group per operating
condition (μ advance ratio, γ flight path angle, C_T thrust coefficient), each holding
scattered (φ, θ) channels of dBA level and EAA excess-attenuation factor.

- **`eval()` computes the ground level:** geometry → sphere angles (`lookup_angles`),
  barycentric interpolation between operating conditions (`--blend` offers two other
  selections; see `condition_selection` below), periodic bilinear interpolation on the sphere
  grid, then spherical spreading plus excess attenuation, `EAA * (range - r_sphere) / 1000`.
  EAA is the total A-weighted attenuation over the 1000 m reference distance the database was
  built at, which is exactly how panam writes and consumes it.
- Ground effect: see `--ground_effect`, on by default wherever the spectral path runs.
- Only `same_grid == true` databases are supported (the interpolation path is unimplemented).
- Free functions (`bracket`, `periodic_bilinear_interp`, `lookup_angles`, the geometry helpers)
  are exported in the header so the unit tests can reach them.
- **Only φ has a period (360).** A θ outside [θ_first, θ_last] goes through `wrap_bracket` to
  the cell (last, first), whose fraction then falls outside [0, 1]: a mild extrapolation that
  leans on the θ_first node. That is D's `periodic_bilinear_interp`, which `GridAxis::locate`
  reproduces; leave it unless the goldens are being regenerated.

**The loader's guards:**

- **Scalars and arrays are read as `double`** whatever the file stores (`nc_get_var` copies the
  external type's bytes unconverted, so a float32 `.nod` would load as a radius of 0 and NaN
  footprints with exit 0), and a "scalar" reads element 0 only.
- **A sphere layout that leaves grid cells without a channel** (they would read as 0 dBA), or
  drops channels that match no angle of sphere 0's axes, is warned about once per distinct
  layout, the first three in full. None of the shipped databases trips it; `B407_M3`'s 934
  zero channels are written zeros, not unreached cells.
- **`condition_origin` is read per group** (`Database::origin`, `origin_weights`): panam's
  `measured`, `extended_flight_path_angle` or `synthesized_hover`. Absent is measured, so the
  shipped databases that predate it, the goldens' among them, run as before; any other value is
  refused. The run's warning counts states whose blend puts less than half its weight on measured
  spheres, other than states on the synthesized hover itself (weight 0.99 or more), which a
  hovering trajectory is meant to use. The threshold is a majority, not fitted. On the 2019 S-76D
  build (58 measured spheres, 62 extended, 3 synthesized hovers) it warns on 21 of the 31
  approaches, whose flare and climb-out over the array lie between the measured conditions (21
  to 1164 states each; 12% of the approaches' blend weight is extended), and on none of the 39
  descents and 3 level runs the build was made from (3% and 1%). The per-observer exposure share
  in `--provenance_file` is the finer measure.
- **The netCDF fill value is refused** in dBA, EAA and amplitude (`refuse_fill_values`). A group
  panam never finished writing holds 9.97e36 everywhere; it is finite, so the non-finite floor
  passed it, and it read as a 1e37 dB level (the September 2026 files, HISTORY).
- **The suspect-cell screen** (`screen_suspect_cells`) runs once the grid and spectra are filled.
  It warns by default and changes nothing, so the default path is untouched; `--mask_suspect_cells`
  floors the entries to `no_energy_floor_db`. Each condition is screened on its own, spread over up
  to eight threads and joined in condition order, which the warning and the CSV list them in; on
  one thread the band medians took 15% of a spectral database's load. Things to keep:
  - **The rules and margins were measured** ([FIDELITY](../FIDELITY.md#database-screen)): a dBA
    cell against the median of its sphere's reached cells (40 dB), a band against the median of
    its own cell's bands (45 dB, excess only). The sound files reach 29.1 and 35.9 dB.
  - **Unreached cells are left out.** They hold the grid's initial {0, 0}, which the layout
    warning reports; taken as measured they would be flagged as zeros. `reached_cells` records
    which entries a channel filled.
  - **Rejected: the median of the eight grid neighbors.** It flagged sound cells beside a strip of
    zeros (half their neighbors zeros), and a second pass without the flagged neighbors still
    left a grid edge cut off behind a strip. A strip wider than half a neighborhood is its own
    median, which hid most of the S-76D columns.
  - **Rejected: the same band in neighboring cells.** Gating makes sound spectra jump 100 dB and
    more between cells, so no margin separates good files from NASA's 8–10 kHz cells.
  - **A masked dBA cell drops its bands too.** The spectral path reads only `spectra_`, and a
    spectrum raised or lowered as a whole has no band out of line with its own cell.
  - **Masking recomputes `sphere_max_` and `max_level_`** from the masked grid. `cull_range` and
    the bounds rest on them, and a kept 1450 dBA maximum keeps every observer in range. The
    bounds themselves are built from the grid on demand, so they follow the mask. The
    `cull_range` check in `a dBA cell far from the rest of its sphere is reported, and dropped
    when asked` is a mutation guard: without the recompute the rest of the suite passed.
  - `--self_test` re-reads the file, so it refuses `--mask_suspect_cells`.
- **`sphere_frame` with no velocity, or velocity along the thrust axis,** takes +x's projection
  for its nose, instead of a zero frame that looks every observer up at (φ 0, θ 90).
- **The sphere writers check `nc_close`'s status**, which is where HDF5's buffering reports a
  full disk.

**Exported spheres:**

- **An exported sphere carries the root the loader interprets**
  (`Database::Sphere::speed_reference`, `tip_speed`, `rotor_radius`, `vehicle_weight_newtons`,
  written at the root where the source had them), so it reloads as air- or ground-referenced as
  its source was.
- `--export_sphere` and `--export_aam` refuse an output that is the database being read (they
  write with `NC_CLOBBER`).
- `interpolate_sphere` warns when the requested condition is outside the envelope (`blend_at`
  clamps it, and the file is labeled with what was asked for).
- `an exported sphere reloads as the sphere it was, cell for cell` and the AAM readback test pin
  the payload. `an AAM export is pitched by its flight path angle` pins the frame's sign (a
  level export misses its nose by 10 sin 8° dB), and the readback test AAM's variable order. They are mutation guards: swapped dBA/EAA, or a dropped load-factor correction,
  passed the rest of the suite (see [testing.md](testing.md#mutation-guarded-tests)).

### `src/predicates.{hpp,cpp}`

Exact orientation and in-sphere tests via Shewchuk expansions, with a filtered double fast
path.

- **The exactness is not optional.** The operating conditions are a prism (164 flight
  conditions × 8 thrust coefficients), so whole families of points are *exactly* cospherical,
  and a rounded in-sphere test answers the same question two different ways, which leaves
  Bowyer-Watson's cavity non-star-shaped.
- **The recursion is allocation-free.** `exact_determinant_recurse` takes its matrix as an
  **array of pointers** and builds its minors the same way: the entries are expansions, which
  are heap-allocated vectors, and copying `(size-1)^2` of them at every level of a factorial
  recursion put `malloc` and `free` at the top of the database-load profile. Callers use the
  `std::vector<Expansion>` overload of `exact_determinant`, which adapts to the pointer form
  once. It asserts `size <= MAX_PREDICATE_SIZE`, so the fixed stack buffers the recursion uses
  guard themselves rather than relying on each entry point's `fatal()`. The exact path is taken
  often here (mass degeneracy is the point), so it is worth keeping allocation-free.
- Assumes normal (non-subnormal) coordinates; conditions are normalized into [0, 1] first.

### `src/delaunay.{hpp,cpp}`

Delaunay triangulation in arbitrary dimension, Bowyer-Watson on those predicates.

- Dimension is a constructor argument, not a template parameter, so a fourth operating
  parameter costs nothing here.
- A simplex whose circumsphere passes exactly through the new point is left standing: a point
  *on* a circumsphere does not violate the Delaunay condition, and leaving it is what keeps the
  choice consistent under mass degeneracy.

### `src/interpolator.{hpp,cpp}`

Barycentric weights over the operating-condition space, on top of `delaunay.hpp`. It replaced
CGAL's `Delaunay_triangulation_3`, and Boost, GMP and MPFR went with it. Databases with fewer than
`dim + 1` conditions (a single exported sphere) skip the triangulation and take the nearest
condition.

### `src/condition_selection.{hpp,cpp}`

`--blend nearest|aam` and `--load_factor_correction one_sided|two_sided`, carried into
`Database::blend_at` (and `source_state`, `interpolate_sphere`) as a `BlendOptions` whose
defaults are the original behavior. These are comparison options, not improvements: scored
against the 2017 data, both alternatives are worse on SEL than barycentric
([FIDELITY](../FIDELITY.md#source-options-blend-and-load-factor-law)). The AAM rule set itself is
in [THEORY §2.3](../THEORY.md). Things to keep:

- **The default path is byte-identical.** `blend_barycentric_is_golden` names
  `--blend barycentric --load_factor_correction one_sided` and compares with the golden, and the
  barycentric runs reproduced the fidelity paper's pinned 6bbd6bb outputs digit for digit.
- **`BlendOptions` is at namespace scope**, with `Database::BlendOptions` an alias. Nested in
  `Database`, its member initializers made `= {}` default arguments of `Database`'s own members
  ill-formed ("default member initializer needed within definition of enclosing class").
- **`nearest` searches the normalized space** the triangulation uses (`normalized_conditions_`),
  with the query clamped. In raw units the path angle decides every distance. Its padding slots
  repeat the chosen index at weight 0, not sphere 0 as a degenerate barycentric blend does, so
  the bounds group it with itself.
- **`aam` is built from the raw conditions, before the constructor normalizes them in place**
  (`aam_layers_`): one `AamSelector` per stored C_T (one in all on a fixed_load_factor
  database, where the law carries C_T), over advance ratio and path angle.
  - Speeds stay in advance ratio. The 20 kt gate is converted per query with the database's own
    tip speed, the vehicle's only when the database records none, and a run with neither is
    refused.
  - Layers group C_T equal to 1e-9 relative, which every shipped load-factor database (8 values)
    meets. One whose C_T is scattered per sphere would make every layer one sphere and pin each
    state to the sphere nearest its C_T, so the run warns on single-sphere layers
    (`aam_single_sphere_layers`).
- **The AAM rules, as read**, each have a mutation-guarded test in
  `tests/test_condition_selection.cpp` or `[blend]` in `test_sphere_database.cpp`:
  - sphere angles rounded to 0.5 degrees (halves away from zero), the state's not;
  - gates are maxima, so a pair exactly as wide still interpolates; 3 degrees for a descending
    state, 5 otherwise;
  - path-angle ties by Table 2-8's *labels* (approach higher, departure lower), not its
    contradictory signs; speed ties to the faster;
  - past a gate, the nearest row, then the nearest speed in it (hierarchical, never a combined
    distance); two spheres in a row at one speed take the nearer stored angle;
  - thrust layers tie to the higher C_T.

  Do not "simplify" the fallback to a Euclidean nearest: AAM selects path angle first.
- **Two-sided puts its floor inside `blend_at`** (`min_load_factor`), so a C_T at or below zero
  (negative g) never reaches the logarithm; `process_trajectory` only counts the states held. A
  negative `level_correction_db` is still one uniform shift of the blend, which is why the
  direction bounds (`group_blends`' per-state scale) and `cull_range` stay exact; the
  `blend_*_windowed_is_exact` CTests run two-sided to hold that. Two-sided on a database whose
  load factors are spheres is refused, as is `--min_load_factor` without two-sided.

### Tests of the blend options

- **Unit tests of each AAM rule** on small synthetic condition sets
  (`tests/test_condition_selection.cpp`): bracketing, both path-angle gates, a pair exactly as
  wide as a gate, ties in path angle and speed, the fallback to the nearest row and then the
  nearest speed, rounding, and duplicates within a row.
- **Database-level tests** (`tests/test_sphere_database.cpp`, tag `[blend]`): nearest searches
  the normalized space; aam converts its gate with the tip speed, refuses a run with none, and
  picks the nearest thrust layer; the direction bounds hold under every policy.
- **CTests:** naming the default reproduces the golden footprint; nearest and aam each move it
  and differ from each other; the windowed peak search stays exact against a march under each.

## Legacy databases (`database_version` < 1)

Databases carry a root-level `database_version`. Absent or 0 means the file predates the
versioned builder, and its **hover spheres have broken directivity**: the upper surface was
completed with `-phi` instead of `180 - phi`. Since the source spheres span only φ in [-90, 90]
(φ = 0 straight down), negating is a left/right flip *within* the lower half that maps the range
onto itself, so the mirrored channels wear lower-hemisphere labels and the upper surface is
missing. Forward-flight spheres in those files were completed correctly; only hover is affected
(in `S-76D_M3.nod`, now in the repository's history only, exactly the 24 spheres at μ = 0).

Only the labels are wrong, so `Database::Database` repairs them on load via
`repair_legacy_hover_directivity()`: for the mirrored half, `phi_correct = wrap(180 +
phi_stored)`. It refuses to touch anything that does not match the broken shape (it must be
confined to [-90, 90], have an even channel count, and have its halves actually mirrored) and
logs how many spheres it fixed. `self_test` applies the same repair, or it would compare a grid
built from repaired labels against queries at the raw ones.

Two consequences worth knowing:

- **Channel angles are read per sphere**, not once from sphere 0. Reusing sphere 0's array
  happened to give the right answer for `S-76D_M3.nod` (the repaired hover labels are
  byte-identical to sphere 0's), but only because every forward sphere there shares one layout.
  It would silently misfile data in any database with heterogeneous layouts.
- **The repair leaves the demo footprint bit-identical**, for the same reason, so `golden/` did
  not need regenerating.

Regenerating the database with the current builder (panam, `DATABASE_FORMAT_VERSION = 1`) makes
the fix-up unnecessary; the loader then trusts the file as written.

## Format 2: frames, lossless spheres, density and pure-tone sources

The root attributes are read in any version;
the sources only in format 2. Things to keep:

- **The frame is chosen once, in `source_state`** (`sphere_frame_kind_`), so every path that
  reads a sphere (the hot loop, the bounds, buildings, the GPU) gets the right frame through
  `SourceState::frame`. `flight_path_sphere_frame` takes the nose as given, climb included,
  which is what `nose_direction` already returns: the heading at hover, with its weighted climb.
  Degenerate cases fall back to `sphere_frame`'s whole frame, not a half-built one. A level
  flight-path database matches the same data labeled level to 1e-9, not to the bit: the two
  frames normalize in a different order.
- **Lossless absorption is one member, `absorption_start_`** (the sphere radius, or 0), used by
  `eval` and `eval_spectrum`. For the default it is the radius, so the expressions are the old
  ones to the bit. Inside the sphere a lossless level is the surface's absorbed over the radius,
  in both paths, which is why `eval_spectrum` floors the path at the radius there. The urban
  model's rays absorb over their whole length and add back the part inside the sphere; for a
  lossless sphere they add back nothing (`inside_sphere` in `urban_point.cpp`). The bounds stay
  bounds: absorbing from the hub only lowers levels, and a gaining EAA already declines to bound.
- **`--condition_density atmosphere` in the vehicle file's own air takes the density that file
  states**, not one recomputed from the pressure the file derived from it, so the
  `condition_density_same_air_is_vehicle` CTest is byte-identical. The ratio is applied only when
  it is not 1. Only niceops forms states, so `state_condition` takes the ratio rather than a
  profile.
- **The sources are gridded over the sorted unique angles**, as a band amplitude is; format 2 has
  no unsorted legacy to permute. Each source must be in every condition, keep its harmonics,
  `phase_reference` and `sync_group`, and have every tone an integer order of its shaft rate
  (relative 1e-5). Stored as complex float, pressure relative to the reference, in
  `[(condition * cells + cell) * tones + tone]`, the layout `tone_directivity` reads.
- **The phase-resolution check** (`widest_phase_steps`) flags a tone at load if any condition has
  a step of more than 90° between neighbors (phi, across the seam when phi spans the circle, and
  theta) both within 20 dB of its peak. It warns, marks the tone energy-only, and refuses
  nothing. It judges the phase as it will be interpolated: unwound about
  the shaft where that is the smoother (below).
- **A shaft-order tone is interpolated unwound** (`ToneSource::unwind_order`,
  `cell_shaft_azimuth_`, `shaft_azimuth`). Harmonic m of a rotor winds its phase m times around
  the shaft, which is what a coarse grid cannot follow: on PSU-WOPWOP's 15° tandem the widest
  neighbor step ran 65-180° for orders 3-66. Multiplied by exp(-i sigma m alpha), sigma the sign
  of the shaft rate and alpha the cell's azimuth about the sphere frame's up axis toward
  starboard, the same spheres step at most 5-17° through order 36. The convention is pinned on
  that tandem: PSU-WOPWOP's rotors turn at +30.05 rad/s about +z (counterclockwise seen from
  above), their stored phase winds as exp(+i m alpha), and the other sign makes every order worse.
  `tone_directivity` sums the unwound corners under the same weights and multiplies the result by
  exp(+i sigma m alpha) at the direction itself, so the nodes give back the stored phase and the
  level is untouched. Each tone takes the unwound phase only where its widest step is the smaller
  (a strongly mixed harmonic, or a fixture whose phase does not wind, stays wound). What it buys
  on the tandem is in [FIDELITY](../FIDELITY.md#pure-tones-against-psu-wopwop). Tests: `a shaft-order
  tone is interpolated with its phase unwound about the shaft` (both senses, exact between nodes,
  a non-winding phase left wound).
- **Only niceops's run evaluates tones** (below), so the constructor refuses `band_content =
  broadband` with tones unless `DatabaseLoadOptions::tones_evaluated`, and warns for `total`.
  Putting the refusal in the constructor covers every other tool that loads a database; `run()`
  sets the flag.
- **A replayed rotor turns at its logged rate** (`--rotor_history`, `ToneEvaluator::Group::replayed`):
  the replayed sources are one coherent group exempt from the shared-rate check a sync group
  takes, since free-running rotors differ by 0.2-2 %; each takes its own logged azimuth and rate
  at every state. Test: `replayed rotors interfere at their logged azimuths, and sound at their
  logged rates` (its unequal-rate pair).
- **`tone_speed_exponent` scales both parts of a tone** (`StateTones::amplitude_scales`): the
  phase-locked pressure by (|Omega_run| / |Omega_db|)^nu and the incoherent part's power by its
  square. The speed test carries an incoherent part so that dropping its scale fails
  (`a tone's pressure follows its rotor's speed where the database gives tone_speed_exponent`).
- **An incoherent part is level only, added in energy per source** (`ToneSource::incoherent_db`):
  it takes the tone's own arrivals with no directivity phase or azimuth, so its direct and
  reflected arrivals still interfere, and it is never summed into a sync group. Its -inf is no
  data and is filled from the neighboring cells, as a tone's is: a writer must give a quiet
  direction a finite floor. PANAM's measured writer once wrote a non-negative fit's zeros as
  -inf, and the fill put energy beside them (a 5 dB spike through the off-receiver ground ratio).
- **Multirotor axes share the helicopter's three slots** (`ConditionAxes`, THEORY 2.1):
  `condition(i)` is (x, z, M_h) there, `blend_at` takes the multirotor axes in its own branch, and the
  helicopter code is untouched, branch for branch, so its blends are the old ones to the bit.
  `active_axes_` names the axes the database varies; only those are normalized and triangulated
  (a database at one weight and air holds one M_h, and blends over (x, z), as a
  fixed_load_factor one does over (mu, gamma)), and the blend pads the unused barycentric slots
  with weight 0 on the first sphere, as before. A query off a constant axis sets that axis's
  clamped flag. `--blend aam` is refused in `blend_at` and in `run()`.

### The multirotor schema (format 2, `condition_axes = multirotor`)

What a builder writes, and the loader checks:

| Where | Name | Meaning |
| --- | --- | --- |
| root attribute | `condition_axes` | `multirotor` |
| root attribute | `rotor_count` | N, a whole number |
| root attribute | `rotor_radius_m` | R, m |
| root attribute | `hover_thrust_coefficient` | C_T,h per rotor, T / (rho pi R^2 (Omega R)^2), at hover |
| root attribute | `thrust_coefficient_density` | `run`: the axes were formed with each condition's own air; anything else (or absent) is refused |
| root variable | `fixed_load_factor` | 0 or absent; set, the database is refused |
| condition group, scalar | `inplane_speed_ratio` | x >= 0, the in-plane air speed over v_h |
| condition group, scalar | `axial_speed_ratio` | z, the air speed along the thrust over v_h, positive in climb |
| condition group, scalar | `hover_tip_mach` | M_h = Omega_h R / c, > 0 |
| condition group, scalars (optional) | `airspeed`, `climb_rate` | the near-field term's flight Mach vector (D = 1 without them, with a warning) |
| condition group, scalars (optional) | `air_density`, `weight_n` | informative; not read |
| source subgroup, scalar | `speed_ratio` | the rotor's \|shaft_rate\| over Omega_h at that condition, > 0 |
| source subgroup, scalar | `shaft_rate` | rad/s, signed: > 0 counterclockwise seen from above (the sphere frame's up axis); one sense at every condition |

The helicopter scalars (`advance_ratio`, `flight_path_angle`, `thrust_coefficient`) are not read.
Every other format 2 rule holds (sources, tones as whole shaft orders of `shaft_rate`, phase
references, `band_content`). At a state, `multirotor_condition` (state_condition.cpp) forms
v_h, Omega_h, x, z and M_h from the vehicle file's weight times the load factor (held at 0.1
below it), the density `--condition_density` names (the run's air by default for these
databases) and the run's sound speed at the state's height; `run()` writes Omega_h into the
state's `Blend::hover_rotor_speed`, and `tone_shaft_rate` and `tone_frequencies` scale the
blended speed ratios by it. `check_vehicle` compares only the rotor radius: the axes are
normalized by weight, and no tip speed indexes them. Tests: `a multirotor state's x, z and M_h
are momentum similarity's` (hover, level, climb, descent, a 1.5 g turn, a tilted disk, wind, the
load-factor floor against hand values), `a multirotor database blends over the axes it varies; a
node is its condition`, `a multirotor database refuses what it cannot place a state with`, and
the CTests `multirotor_*`.

## Pure tones in a run (`src/tones.hpp`, THEORY 4.4)

The tones of a `band_content = broadband` database, per arrival,
with phase. Things to keep:

- **Nothing on the default path.** Every tone term is behind `tone_run` in `process_trajectory`:
  a database without tones takes the old calls, and the goldens are byte-identical. A tone run
  marches every state (`expand_from_peak` is off): the bounds are on the bands alone, and
  `--db_cutoff` is refused for the same reason.
- **`Database::tone_directivity` uses `eval_spectrum`'s cell and weights** (`cell_stencil`,
  `blend_present`), so a tone and its bands are read at the same point. The level is blended in
  dB, the phase is the argument of the complex sum under the same weights: no unwrapping, and no
  dip between neighbors of different phase.
- **One reception time per state** (`tone_arrival`): every arrival's emission time is solved
  against the state's center's direct lag, so hubs and image arrivals compare. That is what gives
  PSU-WOPWOP's far-field e^(+i omega n.d/(c(1 - M.n))) between hubs (`two hubs' relative phase`)
  without any far-field approximation in the code.
- **Q enters conjugated.** `spherical_reflection_coefficient` is written for e^(-i omega t); the
  tones are in e^(+i omega t) like the spheres. The two-ray test (`a tone over the ground is the
  exact two-ray sum`) holds the sign.
- **One combining form for a source and a group** (`combine` in `evaluate`): a group of one adds
  `0 + d`, which is `d`, so `--tone_coherence coherent` and `energy` agree to the bit for one
  source (`one source sums the same coherently or in energy, to the bit`; CTest
  `tones_one_coherent_is_energy`). The turbulence coherence is taken from the state's center for
  the same reason.
- **The band fallback reads the band path's own levels** (`with_tones` in `process_trajectory`):
  with every term and in straight-ray free field, so every path (plates, tables, terrain,
  screening, refraction, urban) gives its correction without its own tone code. The lambdas
  already write band levels when given a buffer; a tone run gives them one even when no history
  is kept. A band the database leaves empty there keeps the free field.
- **Database sync groups must share a shaft rate at every condition**; a run-time group is held at
  its members' mean instead. Only `rotor_azimuth` sources join a coherent group.
- **The A-weighted level weights each tone at its own frequency**; outputs built from the band
  levels afterwards (the C-weighted SEL, EPNL's PNLT, band histories and frames) see the tone's
  power filed in its band and weight it as the band. `add_tone_power` writes that power into the
  history's bands.
- **The ground per frequency uses an admittance table** (`AdmittanceTable`, 96 points per octave),
  not the band admittances: a tone's frequency is its own, Doppler-shifted. The ground table
  (`GroundTable`) is the band model's and is not read for tones.
- **The default phasing is each rotor's recorded azimuth** (`azimuth`, first sphere): the phases
  are stored locked to it (Â = A e^(-i m psi_ref)), so evaluating them at it gives back the
  phasing they were predicted with. Zero would put every rotor of a tandem at one blade azimuth,
  which is not the vehicle that was predicted. A group whose members' recorded azimuths differ
  between spheres by other than a common turn warns: one default cannot hold both.
- **Wind is solved in the frame moving with the air** (`tone_arrival`'s `wind`): the hub at
  V - W, the microphone received at x - W T. The wind is the profile's mean over the heights from
  the ground to the state, the straight ray's own on the level ground the per-frequency path
  runs on. Calm air takes the old arithmetic in its own branch, to the bit (`a tone in a uniform
  wind`), so a run without a wind profile is unchanged.
- **A predicted sphere's floor cells must stay finite.** The loader keeps a non-finite level as
  an empty entry, which `tone_directivity` drops and renormalizes over, as the bands do for
  unmeasured cells: a -inf written where a prediction is merely quiet fills a directivity null
  with its neighbors' level. panam's PT build keeps PSU-WOPWOP's floor level there; on the tandem,
  writing -inf put 15-33 dB into bands made of such tones.
- **Verified against PSU-WOPWOP** in the validation harness (`tones/`): the database's cells
  through `tone_directivity` to float precision; the coherent sum of a tandem's two rotor spheres
  against PSU-WOPWOP's whole-vehicle sphere, converging as the range grows; and ground-microphone
  band levels of a flyover against PSU-WOPWOP's direct prediction, with and without the near-field
  term ([FIDELITY](../FIDELITY.md#pure-tones-against-psu-wopwop)).
- Tested in `[tones]` (`test_sphere_database.cpp`) and the `tones_*` CTests, on databases
  `make_tone_database` writes from B407's (two rotors at 110 and 220 Hz, 95 dB). Deferred: tones on
  refracted and terrain paths per frequency, a bound for the peak search.
- **The near-field term** (`--tone_near_field`, THEORY 4.4) is `ToneSource::near_field`, Lap[D A]
  per condition, cell and tone, computed at load by `compute_tone_near_fields` with
  `sphere_laplacian` (`src/sphere_laplacian.cpp`) and interpolated in `tone_directivity` with the
  amplitude's own weights and unwinding. Keep:
  - **Differentiate whichever of f and its unwound B varies less at each node** (relative differences
    to its four grid neighbors; a tie unwinds). Unwinding everywhere fails twice: order +-1, whose B
    is sin(beta) with a cone point on the shaft (its error grew as the grid was refined, 0.05 at 15°,
    0.26 at 2.5°), and a tone whose content near the shaft does not wind, where B winds instead.
    Choosing by the phase step alone, with a 20 dB floor, picked the wound form near the shaft for
    pure modes (both scored zero there) and lost an order of magnitude.
  - **A lower hemisphere is extended evenly across the horizontal plane**, the half Database reads
    as its mirror image: one-sided fits at its edge were first order and ran to 0.3 relative.
  - **Interior nodes take three-point differences** in (theta, phi); only the grid's poles, an open
    edge and a theta edge without a pole take the tangent-plane quadratic fit (two grid steps).
  - **The cap**, twice the tone's largest |A| over the sphere, only bites within a wavelength or two
    of the hub; the counted condition, the term over half of |A| there, is the expansion's limit.
  - **D's speed** is the condition's own (airspeed, or advance ratio times the root's tip speed),
    over c = 340.29 m/s; a database without it gets D = 1 and a warning.
  - Tested in `[near_field]` (`test_near_field.cpp`): spherical harmonics on 15°, 5° and 2.5° grids,
    full and lower hemispheres, tolerance 0.25 (l h)^2 + 0.1 h; and a port of the rotating source of
    PSU-WOPWOP's near-field check (exact convected Green's function) at R = 2-16 c/Ω, hover and
    M 0.2, shaft orders 3 and 6: the corrected error falls as about 1/R^2, never exceeds the far
    field's, and is under 0.3 of it where (m + 1) c / (2 |Ω| R) <= 1/2. `--tone_near_field off`
    takes the far-field arithmetic unchanged, so its results are those of the far field alone, to the bit.

## Broadband and spectral propagation

- **There is no 4x EAA multiplier.** `eval` applies EAA exactly as panam wrote it. The D
  original's factor of 4, why it was removed, and what that did to the goldens are in
  [HISTORY](../HISTORY.md#departures-from-d).
- **Two propagation paths, and they do not agree.** The default carries one broadband dBA and
  attenuates it with the database's `EAA` factor; `--spectral` carries the per-band spectrum and
  attenuates each band at its own ISO 9613-1 coefficient. `Database::eval_spectrum` mirrors
  `eval`'s bracketing and bilinear cell exactly, so a difference between them is the physics and
  never the interpolation. `--spectral` is the default whenever the database carries spectra
  (and falls back to broadband when it does not), so the goldens, made from `S-76D_load_factors.nod`,
  which has none, are the broadband path. On B407, spectral minus `--no-spectral` is median
  -1.2 dB (-3.2 to +2.3) on the demo track: the humidity and temperature of the run against the
  ones the database's EAA was built at, and EAA being a linear-in-range rate anchored at 1 km
  where the true loss is concave.
- **`spectra_` is condition-innermost** (`[(cell·n_bands + band)·n_conditions + condition]`,
  as its declaration in `sphere_database.hpp` says), so the per-band loop reads the conditions of one band
  contiguously, on one cache line. The netCDF file stores them the other way round, and the
  loader transposes once.
- `--db_cutoff` works on the spectral path through `spectral_direction_bounds`, the
  band-limited level's own bound (see [performance.md](performance.md)).

## Where databases come from

Databases are built by panam, not here; a wrong result may be bad database data rather than bad
code. Which databases ship, and which carry spectra, is in the manual's
[Shipped databases](../MANUAL.md#shipped-databases). Their rebuild history, including the
frame-bearing rebuild, the EC130 standard-day spheres and the fill-value bug, is in
[HISTORY](../HISTORY.md#the-2017-noise-abatement-databases); their scores are in
[FIDELITY](../FIDELITY.md#shipped-databases-approach-scores-per-aircraft).

### Source spheres

The AAM/RNM source spheres live in the group's OneDrive under `Data/Sphere_Database`, one
directory of `.nc` files per aircraft, named as the source names them (AS350, Be206L3, Be407,
Be430, EC130, MD-902, Mi-8, R-44, R-66, SH-60B), each with its own `vehicle.cfg` and
`vehicle.json`. `Be407_spectral.nod` sits beside them: 26 bands, 10–3150 Hz, built from `Be407/`
as below. The built `.nod` files in `databases/` are tracked; the source spheres are far too
large and stay in OneDrive.

### Building a database with spectra

Use panam's `build_empirical_database(..., include_spectra=True)`, which carries the per-band
amplitudes through from the source spheres. Band-limit it (`infreqs=(10.0, 3500.0)`), because
the top few bands of a measured sphere are usually noise-floor artifacts.

### The 2017 Noise Abatement databases

The AS350, B407, Be206L3, EC130, RO-44 and RO-66 databases are built from the 2017 Noise
Abatement flight test by panam (`noise_abatement_2017.build_all`, then
`build_empirical_database`). The validation harness keeps the recipes (`spheres/rebuild.py`,
`spheres/hover_sphere.py`, `spheres/release_database.py`, `spheres/rebuild_database.py`).

- **What each file records at its root:** `speed_reference`, the build atmosphere, and
  `hover_source` when the hover was synthesized from a sphere outside the build; per group,
  `DOPPLER_SHIFT_REMOVED` and `coverage`.
- **B407 flight spheres** are depropagated along curved rays through each run's own atmosphere,
  with a 2-degree emission floor, an adaptive 2-degree grid and tone-aware third-octave bands. A
  rotor harmonic near a band edge is filed whole in its own band instead of being split by whole
  FFT bins.
- **The B407 hover** is synthesized from the slowest level pass with the Doppler shift taken
  out, then corrected by the measured hovers. Two rings were measured: 18–25 degrees below the
  horizon (out of ground effect) and 0.5–3 degrees (in ground effect).
  - Each ring is a Fourier series in azimuth from the heading, with a sum-to-zero offset per
    heading for the wind.
  - The two are blended linearly in elevation from 3 to 17.5 degrees. Below the high ring each
    azimuthal order fades toward the nadir as cos^m.
  - Nothing measures 3–17 degrees or the region under the aircraft, so the hover there is an
    interpolation. The low ring is the noisier of the two (5 dB held-out rms against 2.8).
- **AS350, Be206L3, EC130, RO-44 and RO-66** take the B407's flight recipe and its
  de-Dopplerized hover source, but no hover correction. Only the B407's hovers are tracked: the
  AS350 flew the same hover cards on day 292, but the dataset has no tracking file for them.
  Their hover is the slowest level pass averaged fore to aft (AS350 run 289151, Be206L3 277103,
  EC130 298295, RO-44 229395, RO-66 234520). No approach card is in any of them, because each
  aircraft's approaches are its validation set.
- **Sphere builds need each site's frame bearing.** The harness's stratified rays must use the
  track frame's own +x bearing, not one bearing for every aircraft: the Eglin layout, where
  RO-44 and RO-66 flew, runs along 140°. Predictions are unaffected by this, since they fit
  each track's own bearing.
- **The dataset spells the EC130 station folder `EC130B4_Ground Stations`**, where panam looks
  for `EC130B4_Ground_Stations`; a build that misses it depropagates at a standard day instead
  of the stations' atmosphere.
