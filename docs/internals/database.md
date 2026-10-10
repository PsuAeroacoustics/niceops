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
  layout, the first three in full. None of the 21 shipped databases trips it; `B407_M3`'s 934
  zero channels are written zeros, not unreached cells.
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
  the payload. They are mutation guards: swapped dBA/EAA, or a dropped load-factor correction,
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
(in `S-76D_M3.nod`, exactly the 24 spheres at μ = 0).

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

## Broadband and spectral propagation

- **There is no 4x EAA multiplier.** `eval` applies EAA exactly as panam wrote it. The D
  original's factor of 4, why it was removed, and what that did to the goldens are in
  [HISTORY](../HISTORY.md#departures-from-d).
- **Two propagation paths, and they do not agree.** The default carries one broadband dBA and
  attenuates it with the database's `EAA` factor; `--spectral` carries the per-band spectrum and
  attenuates each band at its own ISO 9613-1 coefficient. `Database::eval_spectrum` mirrors
  `eval`'s bracketing and bilinear cell exactly, so a difference between them is the physics and
  never the interpolation. `--spectral` is the default whenever the database carries spectra
  (and falls back to broadband when it does not), so the goldens, made from `S-76D_M3.nod`,
  which has none, are the broadband path. On B407, spectral minus `--no-spectral` is median
  -1.2 dB (-3.2 to +2.3) on the demo track: the humidity and temperature of the run against the
  ones the database's EAA was built at, and EAA being a linear-in-range rate anchored at 1 km
  where the true loss is concave.
- **`spectra_` is condition-innermost** (`[(cell·n_bands + band)·n_conditions + condition]`,
  `sphere_database.hpp:826`), so the per-band loop reads the conditions of one band
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
