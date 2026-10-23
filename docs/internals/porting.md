# The D port: replacements and fidelity

## D standard-library replacements

The port swapped each D dependency for a small local equivalent. When a behavior looks odd,
check the D-side original named in the header comment.

| D original | C++ replacement |
| --- | --- |
| `std.getopt` | CLI11 |
| `std.json` | nlohmann_json |
| `std.csv` `csvReader!(double)` | `src/util/csv.hpp` (`read_numeric_csv`, plus `read_labeled_csv` for the header-named formats, which have no D original) |
| `std.experimental.logger` | `src/util/logger.hpp` (`log`/`info`/`trace`/`fatal`) |
| `std.numeric.findRoot` | `src/util/rootfind.hpp` (bisection) |
| `std.math.approxEqual` | `approx_equal()` in `sphere_database.hpp` |
| `std.range.iota` (floating point) | `iota_exclusive()` in `src/util/iota.hpp` |
| `unittest` blocks | Catch2 cases in `tests/` |
| `scope(exit) nc_close` | `NcGuard` RAII in `sphere_database.cpp` |
| matplotlib-d | pybind11-embedded matplotlib ([build.md](build.md#plotting)) |
| SPLINTER (vendored B-spline) | `src/smoothing_spline.{hpp,cpp}` |
| `netcdf.d` | the netCDF C API directly |

## Porting pitfalls already hit (do not reintroduce)

The port was validated against the D reference, and these four bugs were found and fixed. Each
was a plausible-looking mistranslation.

- **`iota` is EXCLUSIVE of its end value.** D's `iota(-1500.0, 1500.0, 25.0)` yields **120**
  values ending at **1475**, not 121 ending at 1500 (verified by compiling D). This governs the
  observer grid and the contour levels, which reproduce the off-by-one deliberately through
  `iota_exclusive` (`src/util/iota.hpp`, called from `src/niceops.cpp`; `tests/test_iota.cpp`).
  Comments in the port once asserted the opposite.
- **Don't sort `phis_scattered`/`thetas_scattered` in place.** D builds the unique axis values
  from a *copy* (`.dup.sort().uniq().array()`). Sorting in place collapses the per-channel arrays
  that the grid fill pairs index-wise with `dBAs_scattered`, leaving most cells zero. Symptoms:
  `dba_max` too low, the database loading suspiciously fast (~0.3 s instead of ~2.9 s),
  `--self_test` returning exactly 0, and levels tens of dB low.
- **D default-initializes `double` to NaN, not 0.** `ambient_level` relies on this: the `isNaN`
  guard means D adds no ambient unless `-a` is given. Initializing it to `0.0` adds a spurious
  0 dB floor. The mix itself is no longer D's: `-a` adds the ambient in *energy*
  (`level_with_ambient`, `sel.hpp`), where D summed pressure amplitudes (see
  [HISTORY](../HISTORY.md#departures-from-d)).
- **`approx_equal` mirrors D's `approxEqual`**: the relative term is against the *second*
  argument, with defaults 1e-2 relative and 1e-5 absolute. It is used for database channel
  matching, so tightening it silently drops channels.

General rule: when a D construct looks obvious, check it against Phobos semantics rather than
intuition, and prefer verifying with a throwaway `ldc2` program, which is how the `iota`
question was settled.

## Port fidelity and its exceptions

**Port fidelity outranks correctness, with three deliberate exceptions.** Several `NOTE:`
comments mark places that faithfully reproduce D quirks; leave them alone unless you are also
regenerating the golden outputs, and say so explicitly. The exceptions are the 4x EAA
multiplier, the energy ambient mix and the spreading cap at the sphere surface (see
[HISTORY](../HISTORY.md#departures-from-d)). Each was a physics error rather than a quirk, and each was
measured against the physics before it changed; `golden/197421_LAE.csv` was regenerated from the
C++ build ([testing.md](testing.md#golden-files)). That is the bar for departing from D: a
measured case that the D behavior is wrong, not a reading that it looks wrong.

## Quirks kept on purpose

- `write_levels` emits one extra trailing newline because the D original did.
- `level_along_bearing` (`sphere_database.cpp`) uses `heading[0]` twice where the second
  component wants `heading[1]`. That is the D original's latent bug, kept.
- The `-F` state figures are named `state_CT_.png` and so on, with the separator before the
  extension, as the D original named them ([build.md](build.md#plotting)).
- θ is not periodic, but the sphere lookup reads it through the D original's
  `periodic_bilinear_interp`: in a database whose θ nodes stop short of 0° or 180°, a direction
  past the last node falls in the cell that joins the last node to the first, and is read by a
  slight linear extrapolation from them (the cell bracketing in `sphere_database.cpp` reproduces
  the knot sort). Kept for port fidelity.
- `fatal()` is `[[noreturn]]` and aborts, matching D's logger; RAII cleanup intentionally does
  not run on fatal paths.

## Departures that leave the goldens alone

Two port-time changes touch only figures and options: the contour levels' default (65–100 dB
rather than D's 45–90) and the `--trim` and `--grid` options that replaced `watts_trim` and the
fixed grid extent, both defaulting to the port's behavior. Why, and the command that reproduces
the June 2020 reference figures, are in [HISTORY](../HISTORY.md#the-c17-port-september-2026)
and [Reference figures](../HISTORY.md#reference-figures-june-2020).
