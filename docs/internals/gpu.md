# GPU

The GPU paths: the backend interface, the kernel dialect, how `niceops_point_field --gpu` forms
the rays' amplitudes and sums the bands, and how `niceops_urban_table --gpu` marches the canyon
wave solver. Accuracy is in [FIDELITY](../FIDELITY.md#a-point-source-among-buildings) and
[FIDELITY](../FIDELITY.md#canyon-wave-solver).

## Layout

| Where | What |
| --- | --- |
| `src/gpu/gpu.hpp` | The backend-neutral interface: `Device`, shared `Buffer`s, compiled `Kernel`s, `Batch`es of dispatches. No backend type appears in it. |
| `src/gpu/gpu_metal.mm` | The Metal backend, in Objective-C++. |
| `src/gpu/gpu_disabled.cpp` | The backend of a build without one: `Device::get()` is null. |
| `src/gpu/band_energy` | `BandEnergies`: a batch of ray sets on the GPU, either `RaySet`s (amplitudes from the host) or `PointRaySet`s (rays whose amplitudes the GPU forms). `plan_nodes` picks each problem's form. |
| `src/gpu/point_rays` | `pack_point_rays`: a `PointSourceField`'s rays as the GPU's records. |
| `src/gpu/kernels/` | Kernel sources in the dialect: `prelude_metal.h`; `common.h` (complex floats, the group sum, `nk_phase`); `wedge.h` (Faddeeva, the UTD transition, the wedge coefficient, Filon's panel); the `*_types.h` records shared with the host; `band_energy.kern`, `point_rays.kern`, `canyon_wave.kern`, `fdtd.kern`; `wedge_check.kern`, for the tests. |
| `src/canyon_wave.cpp` (`run_gpu`) | The canyon solver's march on the GPU, beside `run()`. |
| `src/fdtd.cpp` (`march_gpu`) | The 3-D FDTD reference's march on the GPU, beside `march_cpu` ([ground.md](ground.md#the-full-wave-reference-srcfdtd-niceops_fdtd_field)). |
| `cmake/embed_kernels.cmake` | Embeds every file of `NICEOPS_KERNEL_SOURCES` as a string, for `gpu::kernel_source()`. |

A batch dispatches, in order: `ray_amplitudes` (rays not over a shared face, one problem a group),
`face_amplitudes` (each shared face's ray at each problem, 32 threads splitting its Filon
panels), `band_phasors` (each ray's phase at the band's middle), then `band_pair_sum` and
`band_quadrature` over their problems.

## Things to preserve

- **The CPU is the reference.** Every GPU path is opt-in, is tested against the CPU at a stated
  tolerance rather than bit for bit, and never feeds a golden. The CPU path stays byte-identical
  whether or not the build has a backend.
- **Kernels are written once, in the dialect.** A kernel names only what a prelude defines:
  `NK_KERNEL`, `NK_IN`/`NK_OUT`/`NK_VALUE` with their slots, `NK_GROUP_ARGS` (which carries its own
  leading comma), `NK_GROUP`/`NK_THREAD`/`NK_THREADS`, the address-space words for helper
  functions' pointers (`NK_GLOBAL_IN`, `NK_LOCAL`, `NK_SHARED`), `NK_SHARED_ARRAY`, `NK_BARRIER`,
  `NK_CONSTANT` for tables, `nk_sincos`, and `nk_u8` for byte-wide flags. No standard library, exceptions, recursion or `double`;
  complex numbers are `nk_complex` from `common.h`, not a backend's vector type. A record the
  host fills goes in a `*_types.h` of fixed-width integers and floats, included by both sides, so
  the layouts agree. A CUDA or HIP backend needs a `prelude_cuda.h` mapping the same names
  (`extern "C" __global__`, `blockIdx.x`, `__shared__`, `__constant__`, `sincosf`, and an
  `NK_GROUP_ARGS` that is empty) and a `gpu_cuda.cpp` behind `gpu.hpp`, compiling the same sources
  with NVRTC; arguments bind by position, which the slots already are.
- **No shader toolchain.** Metal compiles the sources at run time (`newLibraryWithSource`), which
  needs the framework only; Command Line Tools without Xcode have no offline `metal` compiler. A
  program is the prelude, `common.h`, then the sources a `Device::kernel` call names, in order;
  each program and kernel is compiled once and cached.
- **Safe math** (`MTLMathModeSafe`): the Kahan compensation in the pair sum and the float-float
  phases depend on the compiler not reassociating them.
- **What float cannot hold is formed in double on the host.**
  - Phases: a length and a wavenumber each travel as a pair of floats, their product is formed
    exactly (`fma`) and reduced modulo 2 pi before its sine and cosine (`nk_phase`): about
    3e-7 rad of error, where `k L` in float loses 0.03 rad at 10 kHz and 2 km. A Filon panel's
    phase `k u0` goes the same way; in float it lost 6e-4 rad over tens of metres.
  - A wedge's angles: each term's cotangent and the cosine in its transition function's
    argument, from `wedge_angles` or `WedgeSweep`. Near a boundary the cotangent has a pole and
    the cosine a zero; as floats each keeps its relative precision, where the angles themselves
    could not place the term so finely.
  - A Filon panel's lengths: its start's, and its other nodes' distances from it, so that a panel
    of micrometres beside a cut keeps them.
- **The shared face's aperture is worked out once for the run's wavenumbers**
  (`shared_face_geometry`): nodes, path lengths, rates, amplitudes and the second wedge's angles,
  which do not depend on k, once per ray, where the CPU path works them out again at each of the
  72 sub-bands. The nodes run to the lowest wavenumber's extent, graded as `aperture()` grades
  them, with a panel end at every wavenumber's taper start and extent: without them the Filon
  test's faces come out 0.019 dB from the fine midpoint rule at worst, three times what
  [FIDELITY](../FIDELITY.md#a-point-source-among-buildings) gives with them. **The end
  nodes sit a hair of the shortest aperture inside each segment, not a hair of the segment.** A
  geometrical-optics wave's stationary point lies on its cut, where its amplitude per unit path
  length goes as 1/sqrt(u − u_c), so how much of that the first panel misses depends on how far
  in the first node sits: 1e-7 of a 200 m segment (2e-5 m, the 50 Hz aperture's) moved a 10 kHz
  field by 0.06 dB. The CPU path's own apertures are short at high k, which is why it never
  showed there.
- **The band average takes the cheaper of two forms** (`plan_nodes`): the pair sum, as the CPU's,
  or the Gauss–Legendre integral of |S|² over the sub-band (THEORY 6.5), whenever its nodes times
  the rays are fewer than the pairs within reach. The node count is
  `h Δ / 2 + 4 (h Δ)^(1/3) + 6`, Δ the spread of lengths: on the district's problems it is within
  2e-12 dB of every pair summed. A problem past 4096 nodes takes the pairs.
- **`niceops_point_field --gpu` streams.** Its workers trace pairs without pause and push each
  pair's packed rays (`pack_point_rays`) onto a bounded queue; one GPU thread takes them a batch
  at a time, keeps two batches in flight and writes the levels. The bound holds the workers back if
  the GPU falls behind, and a GPU failure wakes them to stop. The ray dumps come from the CPU.
- **CI never builds the Metal backend.** It runs on Linux, where `NICEOPS_GPU=AUTO` is `OFF`.
  `point_field_gpu_refused` covers the disabled backend there; everything else here is tested
  only on a Mac (`[gpu]`, `point_field_gpu_*`).

## The canyon wave solver

- **One dispatch a step** (`canyon_step`). A group owns a 32 by 8 tile of cells, loads the old
  pressure over the tile and a cell around it, forms the new velocity on every face its
  pressures need -- its own faces, and again its neighbours' on its right and top -- writes its
  own faces, then its cells' new pressure, with `run()`'s rules for the PML, the impedance
  walls, the recording and the source term. Every field is double-buffered, old in and new
  out, so a group reads its neighbours' old faces while they write their new ones. Two
  dispatches a step (velocity, then pressure) took 27.6 s on the 20 m canyon to 707 Hz; one
  takes 23.1 s.
- **The stopping rule is `run()`'s, judged late.** Each step on a 32nd leaves one partial energy
  per tile; the host encodes 1024 steps a batch, keeps two batches in flight, sums the partials
  of each finished batch in order and stops where `run()` would. A batch past that step runs on
  harmlessly, since the field never depends on when the run stops, and the recording is cut
  there. The step the test canyon stopped at was the CPU's.
- **What follows the march is the CPU's and is threaded** (ground.md): with it serial, the GPU
  was busy 10.8 s of a 23 s solve. The table builder marches two canyons at a time on the GPU,
  so that one's transforms overlap the next one's march, and the two share the run's `-j`
  between them (`CanyonWaveSetup::threads`): the solver never takes every core of its own accord.
- **`--gpu` is refused before any work** wherever it would do nothing (`niceops_urban_table`
  without `--wave_below`) or cannot run (`gpu::Device::require`, which both tools call, so the
  refusal they print is one text).
- **A dispatch checks its group size against the pipeline's limit** and refuses an unallocated
  buffer: the kernels' group sums and tiles assume the size they are given, and `ensure`
  allocates at least one element, because a dispatch binds every buffer it names, used or not: a
  batch with no wedges (`--no_diffraction`, an open scene) still binds the wedge buffer.

## Measured

On the 49-block district, an M4 Pro (10 performance and 4 efficiency cores, 20 GPU cores).

- **The pair sum alone** (40 microphones, one source 30 m up, 24 bands of three sub-bands, 1.1
  billion pairs within reach): the CPU's vectorized sum on 14 threads takes 51 ms
  (2.2e10 pairs/s, about 0.6 ns a pair a thread). The same algorithm on the GPU takes 14 ms;
  taking each problem by the cheaper form takes 4.3 ms, 12 times the CPU.
- **The 60-point track over 185 microphones**: 30 s on the CPU, 3.9 s with `--gpu`. The GPU is
  busy 3.3 s of it and sets the pace: the shared faces' panels (`face_amplitudes`) take about
  half of that, the other rays' amplitudes 0.4 s and the band sums the rest. The workers spend
  72% of their time on the path search and 25% packing the rays, two thirds of that the shared
  faces' aperture geometry (415,000 faces on this run, 286 nodes each).

- **The canyon wave solver**: the 20 m canyon to 707 Hz at 15 points per wavelength (1.08
  million cells, 66,764 steps) marches in 10.8 s of GPU time, about 6.7e9 cell updates a second,
  against about 67 s on one core. Fourteen such canyons to 354 Hz take 22 s with `--gpu`
  against 34 s on 14 threads, the CPU's own being bandwidth-bound: 14 threads ran them only 2.6
  times as fast as one.

Measured and left alone:

- **Staging a problem's rays in threadgroup memory** made the pair sum slower (5.6e10 against
  8.2e10 pairs/s): 32 KB a group cut how many groups a core keeps in flight, and the rays were
  already cached.
- **The quadrature threshold**: taking the integral where nodes times rays is under 4, 1, 0.5,
  0.25 or 0.1 times the pairs gave 5.5, 4.3, 4.5, 5.4 and 7.2 ms.
- **Batches the workers wait for**: tracing 64 pairs, then handing them to the GPU, ran as slowly
  as the CPU alone (43.6 s against 42.7 s): the cores idled at each batch's end and through the
  single-threaded hand-off. Streaming took it to 35 s.
- **The CPU's amplitudes, uploaded**: forming them on the CPU and sending them over in float
  (35 s) left the workers the amplitudes' 40% of their time; the GPU forms them in a fraction of
  its 3 s.
- **Fewer aperture nodes.** Most of a shared face's nodes are the panel ends at each of the 72
  sub-bands' taper start and extent. Over 60 of the track's faces, at every sub-band, against
  the midpoint rule at a 192nd of a wavelength: with them, 352 nodes a face and 0.009 dB at
  worst; without, 155 nodes and 0.042 dB; without but graded at 0.05 of the distance rather than
  0.15, 445 nodes and 0.017 dB; with only every third sub-band's, 196 nodes and 0.046 dB. The
  ends are the cheapest way to that accuracy.

## Not on the GPU

- **The path search** (`PointSourceField::trace`): branchy, and most of the workers' time.
- **Ray-tracing hardware** (RT cores; Metal's intersection API, in hardware on M3 and later;
  OptiX; Vulkan's ray queries) accelerates what a renderer spends its time on: walking a bounding
  volume hierarchy and testing rays against triangles. Here that is the visibility tests inside
  the path search, and the path search is 1-2% of a `--buildings` run (the demo approach among
  the test district's eleven buildings, measured by sampling) and 12% of a pair on the 49-block
  district; the rest is each ray's wave arithmetic at each sub-band, which no ray-tracing unit
  does. With `niceops_point_field --gpu` the search is most of the workers' time, but the GPU
  sets that run's pace. Two places it could pay later: scenes of many triangles (pitched LoD2
  roofs, terrain meshes, trees), where the exact test for extruded footprints no longer applies;
  and a diffuse field traced stochastically, as game-audio engines use it, for scattering beyond
  the first facade. Either needs the path tests' contact rules (a segment ending on a wall,
  grazing one, leaving an edge) in float: a conservative filter, with the double test deciding
  what lies near a surface.
- **The edge-source integral** (`--edge_source`, [ground.md](ground.md#the-edge-source-integral-srcedge_source---edge_source)),
  refused with `--gpu`. Its factors are the same at every sub-band from one set of nodes, so the
  natural split is the host forming each ray's amplitudes and the GPU summing them as a `RaySet`
  (amplitudes from the host); that is not done.
- **Apple GPUs have no double precision**, so the PE and the stratified eigenrays, which need
  it, stay on the CPU for Metal. For the PE that was measured, not assumed: in float it holds
  its accuracy only at a range step of a wavelength or less, four times the arithmetic, on a
  march too sequential to fill a GPU ([pe.md](pe.md#grid-and-accuracy)).
