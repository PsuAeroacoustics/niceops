# Units and input checking

## `--units`

Default demo case (no arguments): trajectory `tracks/S-76D/197421AC.csv`, database `databases/S-76D_M3.nod`, vehicle config `databases/S-76D.json`, output `integrated/197421_LAE.csv`. For the full option set see `docs/MANUAL.md`, which documents every flag, or the CLI11 setup in `src/niceops.cpp`.

**`--units` picks the length unit for the whole interface** -- `ft` (default) or `m` -- covering
every length on the command line (`--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
`--receiver_height`, `--path_smoothing`, `--roughness_length`, `--surface_layer_depth`, the
anemometer height in `--wind`), every file carrying a length, and the figure axes. Speeds are not
lengths: `--wind`'s speed names its unit (`12kt`, `6.2mps`, `14mph`), as an `--atmosphere`
file's wind columns must, and a bare number is refused rather than read in `--units`.
Internal computation is metric (`ft2m` and `LengthUnit` in `src/util/constants.hpp`); the
conversion is one seam in `run()`.  Four of those options used to be metric while every file was
feet, which is the inconsistency this removes.

Three things are load-bearing here:

- **A default is a physical quantity, not a number in the current unit.**  Defaults are written
  in meters and skip the conversion, which is why the seam tests CLI11's `count()` rather than
  scaling everything.  Scale the defaults too and `--units m` shrinks the observer grid to a
  third of its extent.
- **`--units` declares what the files already are; it converts nothing.**  The shipped tracks are
  feet, so `--units m` on the demo case reads feet as meters.  Nothing can detect that.  Where a
  column header names its own unit (`alt_ft`, `alt_m`) it is checked against `--units` and a
  disagreement fails the file; a bare `alt` is always accepted.
- **`write_state` stays metric** whatever `--units` says: `golden/state_197421.csv` is the
  D-derived oracle and encodes the internal values.

`fps2g` is chosen as a whole constant rather than composed as `unit.to_meters * mps2g`, because
the product of two constants is not the same double as two multiplies and the state golden is
pinned at 1e-9.  Verified end to end: a metric twin of the demo track under `--units m`
reproduces the golden footprint to 1.5e-12 dB. Trajectory CSVs are the 27-column NASA/FAA AC state format. Vehicle JSON keys: `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg, converted to Newtons), `Atmosphere.density`.

## Input checking and reporting

A bad input used to give a plausible-looking wrong footprint with exit status 0 and nothing on
screen; these now fail, or warn, with the cause named. Know them before adding a reader:

- **Errors and warnings reach stderr** as well as `log/niceops.log` (`Logger::write_line`); `-L`
  still decides what is written at all. A skipped trajectory is `Skipping <file>: <why>` on the
  terminal and exit status 1.
- **An exception in an observer-loop worker is carried to the caller** by `parallel_chunks` (it
  terminated the process at any thread count above one), so `input_error` skips one trajectory
  as `run()`'s handler intends.
- **Readers refuse what they used to swallow.** An AC track's time must increase and every
  column the states are built from must be finite; a geodetic track's time, altitude and
  velocity must be finite; a unit declared in a geodetic header (`Alt (ft)`, `Climb Rate (fpm)`,
  `Time (ms)`) is checked against `--units` (`check_declared_units`); a text line *after* the
  first row of data is an error (`CsvReport`; text above it is a header); observer files keep
  one width; `--atmosphere` rejects non-finite cells. `--trim` past the crop is diagnosed, and
  `-D` on a geodetic track with no clock needs `--start_utc` (it defaulted to UTC midnight, i.e.
  the night penalty, silently).
- **Atmosphere values are range-checked** (`temperature_problem` and friends in
  `atmosphere.hpp`, the bounds `--atmosphere` already used) from the command line and the vehicle
  file. The spectral path's log says where each value came from, and refuses a spectral run whose
  humidity is not given from any source. It is not assumed to be 20 percent:
  the CTest lines that used the old default pin `--humidity 20`. The vehicle file's tip speed,
  weight, rotor radius and C_T are compared with the database's, and a mismatch over 2 percent is
  refused unless `--accept_vehicle_mismatch` is given.
- **A database must vary where the triangulation needs it** (advance ratio, flight path angle,
  and CT unless `fixed_load_factor`) and may not hold two spheres at one condition; both were
  silent NaN or hang before. Test fixtures with one thrust coefficient carry the flag.
- **Several trajectories in one run** name `-i` / `--band_history` files for their track and
  derive `-S`, `--csel_file`, `--epnl_file`, `--band_frames` and `--write_atmosphere` names as
  `-o` does (`<track>_LCE.csv`, `<track>_EPNL.csv`, `<track>_band_frames.nc`,
  `<track>_atmosphere.csv`), and two trajectories asked to write one file are refused.
- **`--terrain` checks the datum and the frame against the track** where the track can say
  (an AC file's lat/lon and alt beside its x/y/z), refuses states under the ground and grids
  that miss the scene, and warns when scattered observers' z disagrees with the terrain; see the
  `src/terrain` entry in [ground.md](ground.md). The grid's heights must be meters.
- **CTests that test these pass on the message text** (`PASS_REGULAR_EXPRESSION`), which
  replaces the exit-status check, so a usage error or a missing file cannot satisfy them;
  `golden_compare` itself is tested on NaN and infinity.

- **`arrays/coyle_approach.csv` is the only scattered-observer set in the repo**, and it is a
  working `-b` input again: every one of its rows was being dropped, because its fields are
  padded *before* the comma and `read_numeric_csv` counted trailing whitespace as junk. The run
  then evaluated zero observers and exited 0. An observer file that yields nothing is now an
  error. All 84 shipped AC tracks pad on the *left*, which `std::stod` always skipped, so the
  goldens never saw this.
