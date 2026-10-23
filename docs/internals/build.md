# Build and CI: the details

The rationale behind the build rules in CLAUDE.md. The basic commands and the CMake options are
in the README; the GCC recipe and the fresh-clone recipe are in CLAUDE.md, the one place they
live.

## Link-time optimization

Link-time optimization is on by default (`-DNICEOPS_LTO=OFF` turns it off), where
`check_ipo_supported` says the toolchain has it. The acoustic loop and the database live in
different translation units, so without LTO each of the ~100 million evaluations of a default
run is an out-of-line call across that boundary. What it is worth, and why, is measured in
[performance.md](performance.md). It leaves the golden comparison passing unchanged, and
GCC 16 and Apple Clang produce identical output with it.

## Asserts

**Live when testing, compiled out in shipped binaries.** `-DNICEOPS_ASSERTS=ON` strips CMake's
`-DNDEBUG` from the optimized configurations and defines `NICEOPS_ASSERTS_REQUIRED`;
`tests/test_build.cpp` fails if such a build still has NDEBUG. CI and the GCC recipe in
CLAUDE.md pass the option. A plain `cmake -DCMAKE_BUILD_TYPE=Release` build is a shipped build,
so its asserts are off (about 4% of the demo's CPU).

The D original ran with its asserts live, and the `assert(std::isfinite(...))` guards in
`sphere_database.cpp` (the sphere level, the excess-attenuation factor, the excess attenuation
and the spreading loss in `Database::eval`, among others) are what the suite exercises, so test
with the option on. Every assert is
side-effect free, so a failing assert is the only place the two builds can differ.

## Command Line Tools without Xcode

With Command Line Tools only, pass the compiler and SDK flags the README gives. No other special
flags are needed.

## GCC and libstdc++

**CI builds on GCC/libstdc++; a macOS clang build will not catch what that finds.** Two classes
of bug are invisible on macOS and fatal on Linux:

- headers that libc++ drags in transitively but libstdc++ does not (`<algorithm>` is the usual
  one);
- static-library link cycles, which ld64 resolves and GNU ld does not (see
  [One static library](#one-static-library)).

So before pushing anything structural, build the tree with Homebrew GCC (the recipe in
CLAUDE.md) rather than waiting for a CI round-trip.

The goldens pass under GCC too, so a compiler change does not move the 1e-9 comparison despite
the different FMA contraction.

**Do not assert bit-identity against a libm call on a literal argument** in a test. The literal
gets constant-folded at build time (GCC folds through MPFR, correctly rounded) while the
function under test reaches libm at run time, and the two differ by an ULP. `acosc` was tested
that way and only GCC noticed.

## Untracked output directories

`state/`, `integrated/`, `log/` and `figures/` are not tracked, so a fresh checkout has none of
them and every writer must create its own parent directory (`ensure_parent_directory` in
`src/util/paths.hpp`). `Trajectory::write_state` once opened its ofstream directly, and so
aborted on a clean CI checkout while passing on any working tree that had run the tool before.
A working tree is therefore *not* a faithful test of a checkout: verify against a clone (the
recipe is in CLAUDE.md).

What such a failure looks like: `fatal()` logs its reason and calls `abort()`, and CTest
reports `Subprocess aborted`. The reason goes to `log/niceops.log` and, since warnings and
errors are echoed to stderr as well (`Logger::write_line`), to the terminal. CI also prints the
log's tail and its fatal and error lines when a job fails, and uploads `log/` as an artifact.

## One static library

**Keep the C++ sources in one static library.** `interpolator.cpp` was once its own
`interpolator` target, a D-era seam that walled off CGAL's heavy headers. Once CGAL went, that
split was a *cycle*: `interpolator.cpp` calls into `Delaunay`, which lives in `niceops_core`,
while `niceops_core` linked `interpolator`. GNU ld walks archives in a single pass and fails
(`undefined reference to niceops::Delaunay::barycentric`); ld64 resolves it either way, **so
this breaks Linux only and a macOS build will not show it**. If you add a target, check that
its dependencies point one direction.

## GPU backend

`NICEOPS_GPU` picks the backend of `src/gpu`: `AUTO` (the default) is Metal on macOS when the
compiler is Clang and none elsewhere, since GCC cannot build the Objective-C++ that drives Metal;
`METAL` and `OFF` force it. The Metal backend enables CMake's `OBJCXX` language for
`gpu_metal.mm` alone (with ARC) and links the Metal and Foundation frameworks. No shader toolchain
is involved: the kernels are compiled at run time from source that `cmake/embed_kernels.cmake`
embeds as strings, rewriting the generated file only when a kernel changes. A new kernel source
goes in `NICEOPS_KERNEL_SOURCES`. Without a backend, `gpu_disabled.cpp` stands in, and a run that
asks for the GPU is refused. See [gpu.md](gpu.md).

## Plotting

Plotting (`-p`, `-f`, `-F`) is the only thing that needs Python. `src/plotting.hpp` is a small
facade in namespace `plt` (`figure`, `contourf`, `tricontourf`, `plot`, labels, `colorbar`,
`savefig`, `show`, `clear`) that replaces the D original's matplotlib-d binding.

- **Embedding.** `src/plotting.cpp` embeds a Python interpreter through pybind11 and calls
  matplotlib directly. It needs `<pybind11/stl.h>` for the `std::vector` conversions: without
  it every call that takes a vector compiles clean and throws `pybind11::cast_error` at run
  time, which is how the whole plotting path once shipped broken (see
  [testing.md](testing.md#the-plotting-gap)).
- **Lazy creation.** The interpreter is started, and `matplotlib.pyplot` imported, on the first
  plotting call (a function-local static), so runs that do not plot never pay the import cost.
- **The bundled environment.** With `NICEOPS_BUNDLE_PYTHON` on (the default), the first
  configure creates a virtual environment in `build/pyenv` and pip-installs matplotlib into it,
  reaching the network only when matplotlib is missing. The build passes the environment's
  site-packages directory as `NICEOPS_VENV_SITE_PACKAGES`, and the embedded interpreter puts it
  at the front of `sys.path`, so a run depends on what the build made rather than on the system
  interpreter. If the environment or the install fails, CMake only warns and the run falls back
  to the system interpreter's packages; CI therefore checks `build/pyenv` explicitly.
- **The `Agg` backend** is selected before pyplot is imported. It is load-bearing for licensing
  as well as for headless runs: matplotlib's Qt backends would pull in PyQt, which is GPL or
  commercial (see [Licensing](#licensing)).
- **Builds without plotting.** `-DNICEOPS_PLOTTING=OFF` drops Python, pybind11 and matplotlib
  from the build and compiles `src/plotting_disabled.cpp` instead: every entry point exists and
  only warns, and `available()` returns false. `run()` refuses `-p`, `-f` and `-F` up front in
  such a build, before anything runs, so the stubs are never reached; they exist so that a
  future caller that forgets the check degrades quietly instead of failing to link. CI's
  "without plotting" job is the only build of that file, so it does not rot unseen.
- **Figure names.** The `-F` state figures are written as `base + "_" + name + "_" + ext`, so
  the separator lands before the extension (`state_CT_.png`, `state_mu_.png`,
  `state_gamma_.png`, `state_z_.png`). That is the D original's naming, which the reference
  figures carry; keep it.
- **The flight track** is drawn over the footprint. It was commented out in the particular
  `niceops.d` the port was derived from, but it is live in that file's ancestor and on the D
  original's master, and every reference figure has it.
- **Contour levels** come from `iota_exclusive(cmin, cmax + cstep, cstep)`, with defaults 65,
  100 and 5 dB (see [porting.md](porting.md) for why not the D original's pair).

No CTest exercises plotting; CI covers it with a separate figure step
([testing.md](testing.md#the-plotting-gap)).

## The run manifest and build information

`--manifest` (`src/manifest.{hpp,cpp}`, `util/sha256.hpp`) records the build, the options, and
the SHA-256 of every input and output. Things to keep:

- **The git state is read at every build, not at configure.** The `niceops_build_info` target runs
  `cmake/build_info.cmake`, which writes `build/generated/niceops_build_info.hpp` only when its
  text changes, so an unchanged tree recompiles nothing and a commit after configuring is still
  seen. A tree with no git records `unknown`. Only `src/build_info.cpp` includes the header.
- **Outputs are found, not listed.** Every writer passes its file through
  `ensure_parent_directory` (`util/paths.hpp`), which records it once a `RunManifest` exists.
  A writer that skips it is missing from the manifest as well as unsafe in a fresh clone. The log
  file goes through it too, and the manifest names it without hashing it: runs append to it.
- **Warnings and errors are counted whatever `-L` says**: `Logger::write` formats and notes
  anything at warn or above before the level test. Below warn, nothing changed.
- **Options are recorded as given, not resolved.** CLI11's `config_to_str(true)` writes an empty
  value for every default it did not capture, which reads as if the default were empty;
  capturing them all would change every `--help` line. The version pins the defaults.
- **Nothing is opened twice for it.** A database's root attributes come from the loader, which has
  the file open (`Database::root_attributes`); a second `nc_open` cost about 0.55 s, since
  netCDF-4 reads every group's metadata at open. Only a database that never loaded is opened at
  write. The inputs are hashed in the background from the start of the run
  (`RunManifest::start_hashing`), so a 67 MB database's 0.2 s overlaps the computation.
- **`main()` writes the manifest**, after `run()` returns or throws, so the status it records is
  the process's. A manifest path that is an input is refused before the run, before the request
  is handed to `main()`, so a refused one leaves nothing to write.

## Licensing

The project is MIT (`LICENSE`). `THIRD_PARTY_NOTICES.md` carries the verbatim notices for
netCDF, CLI11, nlohmann/json, pybind11 and Catch2, all BSD, MIT or BSL and none copyleft.
Regenerate it from the fetched sources under `build/_deps` if a dependency is added, changed or
dropped. The plotting packages in `build/pyenv` are invoked rather than linked and are covered
by a summary paragraph rather than reproduced.

The `Agg` backend choice in `plotting.cpp` is load-bearing for this: matplotlib's Qt backends
would pull in PyQt, which is GPL or commercial.
