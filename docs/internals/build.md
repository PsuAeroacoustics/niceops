# Build and CI: the details

The rationale behind the build rules in CLAUDE.md.

Link-time optimization is on by default (`-DNICEOPS_LTO=OFF` turns it off), where
`check_ipo_supported` says the toolchain has it. It is worth 18% of the acoustic loop,
because that loop and the database live in different translation units; it leaves the
golden comparison passing unchanged, and GCC 16 and Apple Clang produce identical output
with it.

- **Asserts: live when testing, compiled out in shipped binaries.** `-DNICEOPS_ASSERTS=ON` strips CMake's `-DNDEBUG` from the optimized configurations and defines `NICEOPS_ASSERTS_REQUIRED`; `tests/test_build.cpp` fails if such a build still has NDEBUG. CI and the test recipe below pass the option. A plain `cmake -DCMAKE_BUILD_TYPE=Release` build is a shipped build, so its asserts are off (about 4% of the demo's CPU). The D original ran with asserts live, and `sphere_database.cpp`'s `assert(isFinite(...))` guards are what the suite exercises, so test with the option on. Every assert is side-effect free, so a failing assert is the only place the two builds can differ.
- If Xcode proper isn't installed (Command Line Tools only), pass `-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`. No special flags are needed; the vendored SPLINTER library that used to require them is gone.
- **CI builds on GCC/libstdc++; a macOS clang build will not catch what that finds.**
  Two classes of bug are invisible here and fatal there: headers libc++ drags in
  transitively but libstdc++ does not (`<algorithm>` is the usual one), and static-library
  link cycles, which ld64 resolves and GNU ld does not. Before pushing anything structural,
  build the tree with Homebrew GCC rather than waiting for a CI round-trip:

  ```sh
  cmake -S . -B build-gcc -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON \
        -DCMAKE_CXX_COMPILER=$(brew --prefix)/bin/g++-16 -DNICEOPS_PLOTTING=OFF
  cmake --build build-gcc -j8 && ctest --test-dir build-gcc --output-on-failure
  ```

  The goldens pass under GCC too, so a compiler change does not move the 1e-9 comparison
  despite the different FMA contraction. **Do not assert bit-identity against a libm call
  on a literal argument** in a test: the literal gets constant-folded at build time (GCC
  folds via MPFR, correctly rounded) while the function under test reaches libm at run
  time, and the two differ by an ULP. `acosc` was written that way and only GCC noticed.

- **`state/`, `integrated/`, `log/` and `figures/` are not tracked**, so a fresh checkout
  has none of them and every writer must create its own parent directory
  (`ensure_parent_directory` in `src/util/paths.hpp`). `Trajectory::write_state` once
  opened its ofstream directly and so aborted on a clean CI checkout while passing on any
  working tree that had run the tool before. A working tree is therefore *not* a faithful
  test of a checkout — verify against a clone:

  ```sh
  git clone --no-hardlinks . /tmp/fresh && cd /tmp/fresh
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_ASSERTS=ON && cmake --build build -j8
  ctest --test-dir build --output-on-failure
  ```

  Note what such a failure looks like: `fatal()` logs its reason and calls `abort()`, and
  CTest reports `Subprocess aborted`. The log line goes to `log/niceops.log` and, since
  warnings and errors are echoed to stderr as well (`Logger::write_line`), to the terminal --
  it used to reach the log alone, which is why CI still dumps the log on failure.
- **Keep the C++ sources in one static library.** `interpolator.cpp` used to be its own
  `interpolator` target, a D-era seam that walled off CGAL's heavy headers. Once CGAL went,
  that split was a *cycle*: `interpolator.cpp` calls into `Delaunay`, which lives in
  `niceops_core`, while `niceops_core` linked `interpolator`. GNU ld walks archives in a
  single pass and fails (`undefined reference to niceops::Delaunay::barycentric`); ld64
  resolves it either way, **so this breaks Linux only and a macOS build will not show it**.
  If you add a target, check that its dependencies point one direction.

## Licensing

MIT (`LICENSE`), cherry-picked from master, where it had been added after this branch forked.
`THIRD_PARTY_NOTICES.md` carries the verbatim notices for netCDF, CLI11, nlohmann/json,
pybind11 and Catch2 -- all BSD/MIT/BSL, none copyleft. Regenerate it from the fetched sources
under `build/_deps` if a dependency is added, changed or dropped; the plotting packages in
`build/pyenv` are invoked rather than linked and are covered by a summary paragraph rather than
reproduced.

The `Agg` backend choice in `plotting.cpp` is load-bearing for this: matplotlib's Qt backends
would pull in PyQt, which is GPL or commercial.
