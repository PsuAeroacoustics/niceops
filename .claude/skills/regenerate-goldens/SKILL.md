---
name: regenerate-goldens
description: Rebuild golden/197421_LAE.csv and golden/state_197421.csv from the D reference implementation. Use when the goldens must be regenerated because a deliberate numerics change broke regression_golden or state_golden, or when the D oracle is needed to adjudicate a C++/D discrepancy.
---

# Regenerating the golden reference files

> **`golden/197421_LAE.csv` no longer comes from D.** Removing the spurious 4x EAA multiplier
> deliberately departed from the D original, so a D build now reproduces the bug rather than the
> reference. Regenerate that file by running the current `niceops` on the demo case and copying
> `integrated/197421_LAE.csv` over it. Everything below still applies to
> `golden/state_197421.csv`, which is trajectory state and untouched by propagation, and to using
> D to adjudicate a discrepancy in any path other than the EAA attenuation.

`golden/state_197421.csv` comes from the **D reference
implementation**, which was removed from the tree once the port was verified. The
C++ port matches them to **1e-9** (actual max relative difference ~1e-14, a few ULPs
of LDC-vs-clang round-off).

Only regenerate them when a numerics change was deliberate — port fidelity outranks
correctness here, so a golden mismatch usually means the port is wrong, not the golden.

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
  port to just 5e-6. For a 1e-9 golden, patch the writers to `%.17g` (`write_levels` in
  `niceops.d`, `write_state` in `trajectory.d`, both recovered from `c74bab3`) — formatting
  only, no numerics.
- The D CMakeLists hardcodes `/usr/local` paths for netCDF and Eigen; override as above on
  Apple Silicon. D also crashes if `state/` does not already exist (it does not create the
  directory, unlike the port).
