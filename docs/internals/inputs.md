# Units and input checking

What `--units` means to a user, and the refusals a user sees, are in the manual
([Units](../MANUAL.md#units), [Errors and refusals](../MANUAL.md#errors-and-refusals)). This
file is what to preserve in the code behind them, and the rules a new reader must follow.

## `--units`

### The seam

Internal computation is metric (`ft2m` and `LengthUnit` in `src/util/constants.hpp`). The
conversion happens once, in `run()`, where the command line is parsed, and nothing downstream of
that seam sees anything else. The seam converts exactly these options, each only when given:
`--grid`, `--gridstep`, `--trim`, `--ground_elevation`, `--receiver_height`, `--path_smoothing`,
`--roughness`, `--layer_depth`, `--roughness_length`, `--surface_layer_depth`,
`--terrain_offset`, `--terrain_step`, `--screening_min_height`, `--patch_size`,
`--diffuse_radius` and `--path_radius`, an operations file's `offset` column, and a ground
classes file's `layer_depth` and `roughness` columns (`ground_classes_from`). The anemometer
height in `--wind` is a length in `--units` too. Speeds are not lengths: `--wind`'s speed names
its unit (`12kt`, `6.2mps`, `14mph`), as an `--atmosphere` file's wind columns must, and a bare
number is refused rather than read in `--units`. A new length option joins the list above, or `--units m` silently reads it as
feet.

### Things to preserve

- **A default is a physical quantity, not a number in the current unit.** Defaults are written
  in meters and skip the conversion, which is why the seam tests CLI11's `count()` rather than
  scaling everything. Scale the defaults too and `--units m` shrinks the observer grid to a
  third of its extent.
- **`--units` declares what the files already are; it converts nothing.** The shipped tracks are
  feet, so `--units m` on the demo case reads feet as meters. Nothing can detect that. Where a
  column header names its own unit (`alt_ft`, `alt_m`) it is checked against `--units` and a
  disagreement fails the file; a bare `alt` is always accepted.
- **`write_state` stays metric** whatever `--units` says: `golden/state_197421.csv` is the
  D-derived oracle and encodes the internal values.
- **`fps2g` is a whole constant**, not composed as `unit.to_meters * mps2g`, because the product
  of two constants is not the same double as two multiplies, and the state golden is pinned at
  1e-9. Verified end to end: a metric twin of the demo track under `--units m` reproduces the
  golden footprint to 1.5e-12 dB.

How the options stood before `--units` existed is in
[HISTORY](../HISTORY.md#changes-that-affect-results).

## Input checking and reporting

A bad input must fail, or warn, with its cause named; it must never give a plausible footprint
with exit status 0. Know these before adding a reader, and follow the same rules in it.

### Reporting

- **Errors and warnings reach stderr** as well as `log/niceops.log` (`Logger::write_line`); `-L`
  still decides what is written at all. A skipped trajectory is `Skipping <file>: <why>` on the
  terminal and exit status 1.
- **An exception in an observer-loop worker is carried to the caller** by `parallel_chunks`
  (otherwise it terminates the process at any thread count above one), so `input_error` skips
  one trajectory as `run()`'s handler intends.

### What the readers refuse

- **Trajectories.** An AC track's time must increase and every column the states are built from
  must be finite. A geodetic track's time, altitude and velocity must be finite, and a unit
  declared in a geodetic header (`Alt (ft)`, `Climb Rate (fpm)`, `Time (ms)`) is checked against
  `--units` (`check_declared_units`). `-D` on a geodetic track with no clock needs
  `--start_utc`, rather than defaulting to UTC midnight and so to the night penalty.
- **CSV layout.** A text line *after* the first row of data is an error (`CsvReport`; text above
  it is a header). Observer files keep one width, and an observer file that yields no observers
  is an error. A sixth observer column is the microphone's ground class, which must be a whole
  number from 0 to 2^31 − 1 and is checked against `--ground_classes` in `run()`; it needs the
  height column, so no width is ambiguous. `--atmosphere` rejects non-finite cells. `--trim`
  past the crop is diagnosed.
- **Padding.** `read_numeric_csv` accepts fields padded on either side of the comma.
  `arrays/coyle_approach.csv`, the only scattered-observer set in the repository, pads *before*
  the comma, and all 84 shipped AC tracks pad on the left, which `std::stod` always skips. Keep
  both working; how the first was once lost is in
  [HISTORY](../HISTORY.md#changes-that-affect-results).
- **Atmosphere values are range-checked** (`temperature_problem` and friends in
  `atmosphere.hpp`, the bounds `--atmosphere` already uses) from the command line and the
  vehicle file. The spectral path's log says where each value came from, and refuses a spectral
  run whose humidity is not given from any source. It is not assumed to be 20 percent: the CTest
  lines that relied on that pin `--humidity 20`.
- **Vehicle against database.** The vehicle file's tip speed, weight and rotor radius, and for a
  `fixed_load_factor` database its C_T at 1 g, are compared with the database's
  (`process_trajectory`). A difference over 1 percent
  warns, and over 2 percent is refused unless `--accept_vehicle_mismatch` is given. The weight
  is compared at 1 percent because panam builds with g = 9.82 and NICEOPS uses 9.80665. The
  vehicle file's keys are listed in the manual's [Vehicle file](../MANUAL.md#vehicle-file).
- **Databases.** A database must vary where the triangulation needs it (advance ratio, flight
  path angle, and C_T unless `fixed_load_factor`) and may not hold two spheres at one condition;
  either would otherwise give a silent NaN or a hang. Test fixtures with one thrust coefficient
  carry the flag.
- **Several trajectories in one run** derive their output names as `-o` does, and two
  trajectories asked to write one file are refused. The naming is in the manual's
  [Several trajectories](../MANUAL.md#several-trajectories-and-the-day-night-level).
- **The operations file** (`read_operations`, `src/operations.cpp`) has its own reader, because
  `read_labeled_csv` drops every row that is not wholly numeric and a row here holds paths and
  clock times. It is described below.
- **`--atmosphere` warns, and does not refuse,** when a profile of several levels has none at or
  below the ground (its lowest is held down), and when a level above the ground is calm in a
  profile that has wind elsewhere: the thin layer of
  [atmosphere.md](atmosphere.md#the-surface-wind-profile---wind). Refusing would be wrong, because
  `atmospheres/example_windy_calm_floor.csv` exists to run that case. A profile calm throughout
  (`example_inversion.csv`) is still air and does not warn. When `<file>.provenance.json` (the
  builder's sidecar) sits beside an `--atmosphere`, the run manifest lists it as
  `atmosphere_provenance`.
- **`--terrain` checks the datum and the frame against the track** where the track can say (an
  AC file's lat/lon and alt beside its x/y/z), refuses states under the ground and grids that
  miss the scene, and warns when scattered observers' z disagrees with the terrain; see the
  terrain section of [ground.md](ground.md). The grid's heights must be meters.

### The operations file

`--operations` is read by `read_operations` and turned into the per-trajectory lists by
`expand_operations` (`niceops.cpp`), which runs right after the command line is parsed, before
the manifest and every check that reads those lists. Things to keep:

- **The reader is pure.** It checks syntax and that the named files exist, and converts nothing:
  lengths stay in `--units` (the offset and `wind_height` header units are checked against
  `--units` in `run()`, where the unit is settled), and the clock stays local seconds until
  `expand_operations` shifts it by `--utc_offset`.
- **Every refusal names the line** (`<file> line N`), and a missing file names the row's name
  and the path as resolved. The reader keeps line numbers for that, and the run's skip message is
  `Skipping operation <name> (<file> line N): ...`.
- **Rows are never dropped.** A blank line or a `#` line is skipped anywhere; every other line
  after the header is a row, must have the header's field count (a trailing comma is one more,
  empty, field), and is parsed in full: `6.2 ops` is refused, not read as 6.2. A quote is refused
  rather than half-supported.
- **Columns are matched by `LabeledCsv::normalize`**, so `local_time` is `localtime`. An unknown
  column is refused with the known names listed, so a misspelled `cnt` cannot be a file of counts
  of 1. A `start_utc` column is refused with a pointer to `local_time`.
- **Every header unit is checked**, since normalize drops it from the name. Humidity,
  temperature and pressure name their fixed units, the offset and `wind_height` are checked
  against `--units`, `wind_direction` may say degrees, and any label on another column is
  refused: `count (per year)` read as daily counts would be 25.6 dB high and exit 0.
- **Paths resolve against the operations file's folder** and are made `lexically_normal`, so
  the log, the `Loading database` line, the manifest's deduplication (on role and path) and the
  database cache's keys see one spelling. Absolute paths stand.
- **Names are file names**: `[A-Za-z0-9_.-]`, not starting with a dot, and unique ignoring case
  (some file systems ignore it). Two rows with one name are refused, never suffixed, so a name
  never moves when a row is inserted.
- **Per-trajectory options are defaults, not lists.** `expand_operations` refuses more than one
  `-d`, `-v`, `-o`, `-S`, `--atmosphere`, `--wind`, `--humidity`, `--temperature`, `--pressure`,
  `--start_utc`, `--csel_file`, `--epnl_file` or `--band_frames`. A database or vehicle comes from
  the row or from an explicit `-d`/`-v` (`count() > 0`), never from the demo default. An air
  column that is present must be filled in every row unless its option was given once; the
  list-level air checks then run unchanged on the filled lists.
- **`local_time` goes through the clock**, not a weight: `GeodeticInput::clock_replaces_track`
  makes `start_utc + (t - t0)` the clock of an AC, WGS84 or static track whatever clock it
  carries, so the penalties, the period levels and a flight across 22:00 are all exact.
- **Period counts never touch the clock.** `expand_operations` passes them through
  `OperationRun::periods` to `PathJob::period_counts`; the reader refuses them beside `count` or
  `local_time`, and `run()` refuses them with `--laeq_file`. The weighting is in
  [metrics.md](metrics.md).
- **Lateral offsets stay in `--units` until the seam**, which scales `OperationRun::offsets`
  with the other lengths; the job hands each to `GeodeticInput::lateral_offset`. What the reader
  does with it is in [trajectory.md](trajectory.md#lateral-subtracks).

### The ground classes file

`--ground_classes` is read by `read_ground_classes` (`src/ground_classes.cpp`), which follows the
operations file's rules: it is pure, every refusal names the line, rows are never dropped, a
quote is refused, and columns are matched by `LabeledCsv::normalize` with an unknown one
refused and the known ones listed. It has its own reader for the same reason: its rows hold
impedance model names. Beyond those rules:

- **Ids are labels, not row indexes**: whole numbers from 0 to 2^31 − 1 written as digits only
  (`1.0`, `-1` and `1e3` are refused), unique, in any order. A land-cover raster's codes can
  then be a later file's ids.
- **What a row leaves out is the command line's.** `ground_classes_from` starts every class from
  the run's `GroundModel` after the seam, so a row of only an id is that model field for field,
  and receiver height, turbulence and the boundary-loss switch always come from it. A turbulence
  column is refused with the reason: turbulence is the air's.
- **Units.** `layer_depth` and `roughness` are lengths in `--units`, converted by
  `ground_classes_from` with the same single multiply the seam uses, and a header unit on them
  is checked against `--units`. The flow resistance may declare kPa and the porosity rate 1/m
  (a material rate whatever `--units` says, as on the command line); any other label is refused.
- **Values** carry the limits the models need: a positive flow resistance and layer depth, a
  non-negative porosity rate and roughness. A non-passive admittance can only be found once the
  database's bands are known, in `process_trajectory`.
- **The observer file and the classes name each other** (`run()`): a sixth column without
  `--ground_classes`, `--ground_classes` without one (or with grid observers), and a microphone
  on an id the file does not define are each refused with the row named.

### The profile builder's column reader

`niceops_build_atmosphere` reads its columns with `read_column_source`
(`src/atmosphere_builder.cpp`), which finds them with `find_profile_columns`, the same names and
unit rules as `--atmosphere`, so a header one accepts the other does. Beyond `--atmosphere`'s
rules (non-finite cells, text after data, a height unit that disagrees with `--units`, a wind
without its unit, two levels at one height):

- **Fill values are explicit.** A cell equal to a `--missing` value is dropped for its variable
  at its level and counted in the log and the provenance. The reader never guesses a sentinel
  from its size.
- **Any other implausible value is refused** with the file and data row, unless
  `--drop_implausible`, which drops and counts it. The ranges are `--atmosphere`'s, except that a
  relative humidity up to 110 is read and left to the build, which refuses one above 100 or clips
  it (`--clip_humidity`).
- **Every variable is optional**, so a file carrying only wind or only temperature reads; one
  carrying none of them is refused. A variable two files carry is refused unless one is named.
- **Lengths given on the command line** (`--z0`, `--top`, `--levels`, `--surface_layer_depth`,
  `--ground_elevation`, `--gates`) are converted once, in the tool's `main`, when given.

### Tests

- **CTests of refusals and warnings pass on the message text** (`PASS_REGULAR_EXPRESSION`),
  which replaces the exit-status check, so a usage error or a missing file cannot satisfy them.
- `golden_compare` itself is tested on NaN and infinity.
