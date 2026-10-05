# The D port: replacements and fidelity

## D standard-library replacements

The port swapped each D dependency for a small local equivalent; when a behavior looks odd, check the D-side original named in the header comment.

| D original | C++ replacement |
| --- | --- |
| `std.getopt` | CLI11 |
| `std.json` | nlohmann_json |
| `std.csv` `csvReader!(double)` | `src/util/csv.hpp` (`read_numeric_csv`, plus `read_labeled_csv` for the header-named formats, which have no D original) |
| `std.experimental.logger` | `src/util/logger.hpp` (`log`/`info`/`trace`/`fatal`) |
| `std.numeric.findRoot` | `src/util/rootfind.hpp` (bisection) |
| `std.math.approxEqual` | `approx_equal()` in `sphere_database.hpp` |
| `unittest` blocks | Catch2 cases in `tests/` |
| `scope(exit) nc_close` | `NcGuard` RAII in `sphere_database.cpp` |
| matplotlib-d | pybind11-embedded matplotlib |
| SPLINTER (vendored B-spline) | `src/smoothing_spline.{hpp,cpp}` |
| `netcdf.d` | the netCDF C API directly |

## Code Conventions

## Porting pitfalls already hit (do not reintroduce)

The port was validated against the D reference and these four bugs were found and fixed. Each was a plausible-looking mistranslation:

- **`iota` is EXCLUSIVE of its end value.** D's `iota(-1500.0, 1500.0, 25.0)` yields **120** values ending at **1475**, not 121 ending at 1500 (verified by compiling D). This governs the observer grid and the contour levels. Comments in the port previously asserted the opposite.
- **Don't sort `phis_scattered`/`thetas_scattered` in place.** D builds the unique axis values from a *copy* (`.dup.sort().uniq().array()`). Sorting in place collapses the per-channel arrays that the grid fill pairs index-wise with `dBAs_scattered`, leaving most cells zero. Symptoms: `dba_max` too low, database load suspiciously fast (~0.3 s instead of ~2.9 s), `--self_test` returning exactly 0, levels tens of dB low.
- **D default-initializes `double` to NaN, not 0.** `ambient_level` relies on this: the `isNaN` guard means D adds no ambient unless `-a` is given. Initializing it to `0.0` adds a spurious 0 dB floor. (The mix itself is no longer D's: `-a` adds the ambient in
  *energy* (`level_with_ambient`, `sel.hpp`), where D summed pressure amplitudes, which added two
  uncorrelated levels coherently -- +6.02 dB over two equal levels, up to 3 dB high -- and was never
  tested. A departure from D by the port-fidelity bar: measured against the physics, not a reading.)
- **`approx_equal` mirrors D's `approxEqual`**: relative term against the *second* argument, defaults 1e-2 relative / 1e-5 absolute. It is used for database channel matching, so tightening it silently drops channels.

General rule: when a D construct looks obvious, check it against Phobos semantics rather than intuition — and prefer verifying with a throwaway `ldc2` program, which is how the `iota` question was settled.

- **Port fidelity outranks correctness, with one deliberate exception.** Several `NOTE:`
  comments mark places that faithfully reproduce D quirks; leave them alone unless you are also
  regenerating the golden outputs, and say so explicitly. The exception is the 4x EAA multiplier,
  removed above: it was a physics error rather than a quirk, it was measured against the model it
  stood in for before being removed, and `golden/197421_LAE.csv` was regenerated from the C++.
  That is the bar for departing from D — a measured case that the D behavior is wrong, not a
  reading that it looks wrong.
  - Contour defaults are 65/100, not c74bab3's 45/90: the older pair puts its top level below the demo case's 93.6 dB peak, leaving it unfilled, and the shipped reference figures use 65/100. Figures only -- the goldens are CSV.
  - `watts_trim` and the observer grid extent are now CLI options (`--trim`, `--grid`), both defaulting to the port's behavior so the goldens are untouched. With `--trim -6000 --grid -2000 1100 -1500 1600` the demo case reproduces the June 2020 reference figures to within 1.05 dB on every cell.
  - D's `iota` is EXCLUSIVE of the stop value; observer grids and contour levels reconstruct that off-by-one deliberately (`iota_exclusive` in `src/niceops.cpp`).
  - `write_levels` emits one extra trailing newline because the D original did.
  - `level_along_bearing` uses `heading[0]` twice where the second component wants `heading[1]`; that is the D original's latent bug, kept.
  - `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does not run on fatal paths.
