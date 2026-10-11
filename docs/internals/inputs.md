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
`--terrain_offset` and `--terrain_step`. The anemometer height in `--wind` is a length in
`--units` too. Speeds are not lengths: `--wind`'s speed names its unit (`12kt`, `6.2mps`,
`14mph`), as an `--atmosphere` file's wind columns must, and a bare number is refused rather than
read in `--units`. A new length option joins the list above, or `--units m` silently reads it as
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
  is an error. `--atmosphere` rejects non-finite cells. `--trim` past the crop is diagnosed.
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
  is compared at 1 percent because panam builds with g = 9.82 and NICE-OPS uses 9.80665. The
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
- **Period counts and lateral offsets are read and refused in one place** in `run()`, right after
  `expand_operations`, until the run uses them.

### Tests

- **CTests of refusals and warnings pass on the message text** (`PASS_REGULAR_EXPRESSION`),
  which replaces the exit-status check, so a usage error or a missing file cannot satisfy them.
- `golden_compare` itself is tested on NaN and infinity.
