---
name: regenerate-goldens
description: Rebuild golden/state_197421.csv from the D reference implementation, or build the D oracle to adjudicate a C++/D discrepancy. golden/197421_LAE.csv and the other goldens come from the C++ build.
---

# Regenerating the D-derived golden

`golden/state_197421.csv` (trajectory state, checked by `state_golden`) is the only golden that
comes from the **D reference implementation**, which was removed from the tree once the port was
verified. The C++ port matches it to **1e-9** (actual max relative difference ~1e-14, a few ULPs
of LDC-vs-clang round-off).

`golden/197421_LAE.csv` and every other golden come from the C++ build: regenerate one by running
the current `niceops` on its case and copying the output over it (the table is in
`docs/internals/testing.md#golden-files`). A D build would reproduce what the C++ deliberately
departed from.

D can adjudicate any path except the deliberate departures: the EAA multiplier, the energy
ambient mix, and the spreading cap at the sphere surface (see
docs/HISTORY.md#departures-from-d).

Only regenerate the state golden when a numerics change was deliberate, and say so in the
commit. Port fidelity outranks correctness here, so a golden mismatch usually means the port is
wrong, not the golden.

## Build the D original in a scratch worktree

Do **not** do this in the working tree:

```sh
# c74bab3 is the last commit that still contains the D sources and the D-capable
# CMakeLists; 9ccbf0a is the D-era build (c74bab3's CMakeLists is C++-only).
git worktree add --detach /tmp/dref 9ccbf0a
git show c74bab3:niceops.d > /tmp/dref/niceops.d   # see below: 9ccbf0a's copy is NOT the oracle
cmake -S /tmp/dref -B /tmp/drefbuild -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_D_COMPILER=$(which ldc2) -DNCDF=/opt/homebrew/lib/libnetcdf.dylib \
      -DCMAKE_CXX_FLAGS=-I/opt/homebrew/include/eigen3 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build /tmp/drefbuild --target niceops
```

## Three things will bite you

- **Take `niceops.d` from `c74bab3`, not from the D-era commit.** The port was derived
  from a later revision of it: contour defaults 45/90 (the D-era copy has 65/100), and
  `watts_trim()`/`adapt()` both commented out. The older copy runs trimming and spline
  adaptation, computes a different trajectory, and yields a golden the port can never match.
- **D prints only 6 significant figures** by default, so a stock D build agrees with the
  port to just 5e-6. For a 1e-9 golden, patch the writers to `%.17g` (`write_levels` in the
  `niceops.d` copied from `c74bab3`, `write_state` in `trajectory.d`, which is the same in both
  commits) — formatting only, no numerics.
- The D CMakeLists hardcodes `/usr/local` paths for netCDF and Eigen; override as above on
  Apple Silicon. D also crashes if `state/` does not already exist (it does not create the
  directory, unlike the port).
