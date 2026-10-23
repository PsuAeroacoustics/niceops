# Developer guide

This directory is the developer documentation for NICEOPS. Each file covers one area of the
code: what to preserve there, why it is built the way it is, the designs that were tried and
rejected, what it costs, and the tests that guard it. **Read the matching file before you change
a module.** Much of what looks over-engineered was measured, and many of the obvious
simplifications were tried and failed.

## How the docs are divided

| File | Reader | What goes there |
| --- | --- | --- |
| [`README.md`](../../README.md) | Anyone arriving | What NICEOPS is, building, a first run, the units trap, an index of the docs. Kept short. |
| [`docs/MANUAL.md`](../MANUAL.md) | End users | How to use each option and file: what it does, when to use it, its limits and cost, and the refusals. No equations beyond one-liners, no measurements, no code or test names, no history. |
| [`docs/THEORY.md`](../THEORY.md) | Users and reviewers who need the equations | The models as run: equations, conventions and references. A change to a model's equations updates it in the same commit. |
| [`docs/FIDELITY.md`](../FIDELITY.md) | Anyone deciding how far to trust a model | Every measured and validation result, stated once; the other docs link to it. |
| [`docs/HISTORY.md`](../HISTORY.md) | Anyone comparing with older results | The D original, the port, the departures from D, and every change that moved results. "Used to" and "before" wording belongs only here. |
| `docs/internals/` (this directory) | Developers | Things to preserve, rationale, rejected designs, measured costs, guard tests: one file per area (see the [module map](#module-map)), among them [`urban.md`](urban.md) for propagation among buildings and [`gpu.md`](gpu.md) for the GPU backends. |
| [`CLAUDE.md`](../../CLAUDE.md) | AI coding agents | The rules and invariants, pointing here for detail. |

Dates, commit hashes, harness paths and session narrative stay out of MANUAL and THEORY.

## Pipeline

The pipeline in `run()`/`process_trajectory` (`src/niceops.cpp`):

1. Load the `Database` and the `Trajectory`.
2. Build a `Database::SourceState` once per trajectory state (frame, blend).
3. For each observer and state, evaluate the level at the microphone, either broadband dBA
   with EAA or per band with ISO 9613-1 absorption.
4. Integrate in reception time over the outermost 10-dB-down window (`src/util/sel.hpp`).
5. Write the CSVs and optional figures.

A default run makes about 94 million evaluations, so anything left in the inner loop costs
about 100 million times over.

## Module map

Sources are in `src/`; a name without an extension means its `.cpp` and `.hpp`.

| Area | Sources | Notes |
| --- | --- | --- |
| Command line, the observer loop, metrics; hot loop, SEL window, windowed peak search, load time, rejected optimizations | `niceops.cpp`, `sphere_database.cpp`, `util/sel.hpp`, `util/parallel.hpp` | [performance.md](performance.md) |
| Sphere database, Delaunay and exact predicates, the sphere lookup, `--blend` (including `nearest`/`aam`), legacy repair, broadband vs spectral | `sphere_database`, `predicates`, `delaunay`, `interpolator`, `condition_selection` | [database.md](database.md) |
| ISO 9613-1, the `RayPath` seam, `--atmosphere`, `--wind`, the profile builder | `atmosphere`, `propagation.hpp`, `atmosphere_profile`, `atmosphere_builder` | [atmosphere.md](atmosphere.md) |
| `--refraction linear` / `stratified`, shadow zones, `ShadowTransition`, ducts, `--tunneling` | `refraction`, `stratified`, `airy` | [refraction.md](refraction.md) |
| Parabolic equation (split-step Padé, Crank–Nicolson), benchmark exact solutions, `--propagation_table` | `pe`, `pade_propagator`, `exact_solutions`, `propagation_table`, `pe_table`, `util/lanes.hpp`, `tools/` | [pe.md](pe.md) |
| Receivers (plates, poles, `--ground_classes`), ground effect, Faddeeva, weighting, `--terrain`, `--terrain_screening` | `ground_effect`, `ground_classes`, `faddeeva`, `weighting`, `terrain`, `screening`, `plate_receiver`, `terrain_ray` in `propagation.hpp` | [ground.md](ground.md) |
| Propagation among buildings: `--urban_table`, `--buildings` and how its rays add, the edge-source integral, facade relief, the roof-row profile PE, turbulence scattering into shadows, the canyon solver and the FDTD reference (`niceops_fdtd_field`) | `urban`, `wedge`, `double_wedge`, `edge_source`, `canyon_wave`, `fdtd`, `buildings`, `building_paths`, `point_paths`, `edge_visibility`, `diffuse_facades`, `facade_relief`, `grating_paths`, `profile_pe`, `turbulence_scattering`, `geojson_buildings`, `building_input.hpp`, `urban_point`, `source_coherence`, `plan_grid.hpp` | [urban.md](urban.md) |
| EPNL, reception-time band frames, LAmax/NAx/TA, DNL/Lden/CNEL, C-SEL | `pnl`, `band_frames`, `util/time_weighting.hpp` | [metrics.md](metrics.md) |
| AC and geodetic tracks, smoothing spline, hover nose, lateral subtracks, `--nose`, `--speed_reference`, `niceops_procedure` | `trajectory`, `geodetic`, `smoothing_spline`, `util/csv.hpp`, `procedure`, `state_condition`, `tools/procedure.cpp` | [trajectory.md](trajectory.md) |
| `--units`, input checking, the operations file (`--operations`) | `run()` in `niceops.cpp`, the readers, `operations` | [inputs.md](inputs.md) |
| Goldens, `self_test`, mutation-guarded tests | `golden/`, `tests/` | [testing.md](testing.md) |
| Replacements for D and its libraries, port quirks kept on purpose | `util/iota.hpp`, `util/rootfind.hpp`, `util/logger.hpp`, `vec3.hpp` | [porting.md](porting.md) |
| Contour polygons (`--contour_file`) | `contour` | [metrics.md](metrics.md#contours) |
| Gridded levels as CF netCDF (`--raster_file`) | `raster` | [metrics.md](metrics.md#raster) |
| Energy by sphere (`--provenance_spheres`) | `sphere_energy` | [metrics.md](metrics.md#source-provenance) |
| Build, link layout, shared helpers, the run manifest | `CMakeLists.txt`, `util/paths.hpp`, `util/netcdf.hpp`, `util/netcdf_write.hpp`, `util/constants.hpp`, `util/hermite.hpp`, `manifest`, `build_info`, `util/sha256.hpp`, `cmake/build_info.cmake` | [build.md](build.md) |
| Plotting (pybind11-embedded matplotlib, created lazily, `Agg` backend) | `plotting`, `plotting_disabled.cpp` | [build.md](build.md#plotting) |
| GPU backends, the kernel dialect, band energies and point rays on the GPU (`niceops --buildings`, `niceops_point_field`), the canyon solver's and the FDTD's marches | `gpu/gpu.hpp`, `gpu/gpu_metal.mm`, `gpu/gpu_disabled.cpp`, `gpu/band_energy`, `gpu/point_rays`, `gpu/band_stream`, `gpu/kernels/`, `cmake/embed_kernels.cmake` | [gpu.md](gpu.md) |

## Repository layout

- `src/`: the implementation. Every `.cpp` lives here and is listed in `CMakeLists.txt`.
- `tests/`: Catch2 unit tests, the `golden_compare` tool, the CTest helper scripts and
  programs, and test inputs in `tests/data/`.
- `golden/`: reference outputs. The regression tests compare against these at a relative
  tolerance of 1e-9; [testing.md](testing.md#golden-files) says where each comes from.
- `tools/`: developer programs; see [Companion tools](#companion-tools).
- `databases/`, `tracks/`, `arrays/`, `atmospheres/`, `plates/`: shipped data; the user's view
  is in [MANUAL.md](../MANUAL.md#input-files).
- `*batch.sh`: batch scripts from 2020. Each one runs every shipped track of one aircraft; the
  `*nqbatch.sh` variants run with `--qsam=false`.
- `docs/`: the documentation above. `.claude/skills/` holds the `regenerate-goldens` skill.
- `reports/`: snapshots of literature reviews (grazing incidence, hemisphere models, urban
  propagation), the propagation roadmap and code reviews. They are not maintained
  documentation.
- `research_notes/`: research behind those reports, and the source-model fidelity studies
  ([Fidelity studies 2026-10-05](../../research_notes/Fidelity%20studies%202026-10-05.md)).
  Also not maintained documentation.
- `state/`, `integrated/`, `log/`, `figures/`, `tables/`: run outputs. These are not tracked,
  and each writer creates what it needs. `tables/` is not shipped data: it is where the CTests
  write the propagation tables they then run (`ground_model_test.nc` from
  `niceops_propagation_table`, `bump_test.nc` and `bump_absorbed_test.nc` from the hidden Catch2
  case `[.bump_table]`), and where the manual's examples put the tables they build (`*.nc` is
  ignored, except `plates/*.nc`).

## Companion tools

The C++ tools in `tools/` are built alongside `niceops` and link `niceops_core`:

| Executable | Source | Purpose |
| --- | --- | --- |
| `niceops_propagation_table` | `tools/propagation_table.cpp` | Builds `--propagation_table` files, by the PE (`--method pe`, the default) or by straight rays with the ground model (`--method ground_model`). |
| `niceops_urban_table` | `tools/urban_table.cpp` | Builds `--urban_table` files for microphones in ideal street canyons or among LoD1 buildings. |
| `niceops_point_field` | `tools/point_field.cpp` | The field of a point source among LoD1 and LoD2 (pitched-roof) buildings at microphones, band by band. |
| `niceops_fdtd_field` | `tools/fdtd_field.cpp` | The full-wave (FDTD) field of a point source among rigid buildings, with the same inputs and output as `niceops_point_field`: the reference that tool is held to. |
| `niceops_pe_compare` | `tools/pe_compare.cpp` | The refraction models against the parabolic equation on one path. |
| `niceops_pe_field` | `tools/pe_field.cpp` | The PE's whole field on one vertical plane, for plotting. |
| `niceops_ray_geometry` | `tools/ray_geometry.cpp` | The refracted ray's geometry for a batch of source–receiver pairs, so ground models outside `niceops` can be evaluated on the same rays. |
| `niceops_benchmark` | `tools/benchmark.cpp` | The Attenborough et al. (1995) benchmark cases by the PE and, where one exists, the exact solution. |
| `niceops_build_atmosphere` | `tools/build_atmosphere.cpp` | `--atmosphere` profiles from measured columns (`build`), and a site roughness length from wind profiles (`fit_z0`); see [atmosphere.md](atmosphere.md#the-profile-builder-niceops_build_atmosphere). |

`niceops_procedure` (`tools/procedure.cpp`) is described with the trajectories in
[trajectory.md](trajectory.md). The manual covers the eight that users run:
[niceops_procedure](../MANUAL.md#niceops_procedure),
[niceops_propagation_table](../MANUAL.md#niceops_propagation_table),
[niceops_build_atmosphere](../MANUAL.md#niceops_build_atmosphere),
[niceops_urban_table](../MANUAL.md#niceops_urban_table),
[niceops_point_field](../MANUAL.md#niceops_point_field),
[niceops_fdtd_field](../MANUAL.md#niceops_fdtd_field),
[niceops_pe_compare](../MANUAL.md#niceops_pe_compare),
[niceops_benchmark](../MANUAL.md#niceops_benchmark).

The Python scripts:

- `plot_pe_compare.py` plots `niceops_pe_compare` output, and `plot_benchmark.py` plots
  `niceops_benchmark` output in the layout of the Attenborough et al. figures.
- `pnl_reference.py` generates the PNL/EPNL reference values that `tests/test_pnl.cpp` embeds,
  from panam's own functions (it needs a panam checkout).
- `thick_barrier_bem.py` solves a rigid block on rigid ground by boundary elements, for a line
  source: the reference values `tests/test_point_paths.cpp`, `tests/test_pitched_roofs.cpp` and
  `tests/test_profile_pe.cpp` embed (it needs NumPy and SciPy).
- `profile_scattering.py` computes a facade relief's scattering coefficient and grating orders
  per band by the same boundary elements: the files `--facade_relief` reads
  ([urban.md](urban.md)).
- `urban_scattering_study.py` measures how far street microphones' SEL depends on the facades'
  scattering coefficient, over synthetic districts built with `niceops_urban_table`.
- `score_fidelity_frames.py`, `build_hover_candidate.py`, `export_emission_samples.py` and
  `regrid_emission_samples.py` are experimental source-model studies. They need NumPy, and all
  but `score_fidelity_frames.py` read the validation harness (`--harness`); each writes only to
  the output it is given, and none of them builds a production database.
  `tests/test_fidelity_tools.py` checks their numerics, outside CTest. How they were run, and what
  they found, is in the
  [fidelity studies note](../../research_notes/Fidelity%20studies%202026-10-05.md).

## Propagation roadmap

The plan for propagation (refraction first, then urban) is the Claude Doc "NICEOPS:
Atmospheric and Urban Propagation Options".
[`reports/`](../../reports/NICEOPS%20Atmospheric%20and%20Urban%20Propagation%20Options.md) holds
a snapshot of it. The internals files name roadmap items by what they are (the
linear refraction model, the stratified rays, the PE, the propagation tables) and link here
rather than citing step numbers.

The urban model's next steps -- the paths it misses, ground and receivers, cost, drone and AAM
sources, and the validation still open -- are in
[`reports/NICEOPS Urban Propagation Next Steps.md`](../../reports/NICEOPS%20Urban%20Propagation%20Next%20Steps.md).
