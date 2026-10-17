# PE, the benchmark's exact solutions, and propagation tables

The parabolic equation (PE) and the benchmark's exact solutions are the full-wave references of
the [propagation roadmap](README.md#propagation-roadmap): they check the ray models, and the PE
builds propagation-loss tables. The PE never runs inside a footprint. The residue series does,
at setup: it tabulates the shadow-transition field
([refraction.md](refraction.md#shadowtransition-the-level-near-a-shadow-boundary)) and, under
`--shadow_model resolved`, the cramped bands' field. The scheme and the benchmark are in
[THEORY §9](../THEORY.md); the comparisons with exact solutions, with the benchmark and with the
ray models are in [FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions). This
file is what to preserve.

Sources: `src/pe.{hpp,cpp}`, `src/pade_propagator.{hpp,cpp}`, `src/pe_table.{hpp,cpp}`,
`src/exact_solutions.{hpp,cpp}`, `src/airy.{hpp,cpp}`, `src/propagation_table.{hpp,cpp}`,
and the tools `tools/pe_compare.cpp` (`niceops_pe_compare`), `tools/benchmark.cpp`
(`niceops_benchmark`) and `tools/propagation_table.cpp` (`niceops_propagation_table`).

## Scheme and starter

- **A finite-difference PE, not a Green's-function PE (GFPE).** The impedance condition
  dpsi/dz = -i k beta psi goes straight into the grid, so no surface-wave pole has to be handled
  by hand. The default fourth-order ground row carries psi_0 as an unknown (below); only the
  second-order row (`boundary_order` 2, and Crank–Nicolson) eliminates it by a one-sided
  difference. The convention is exp(-i omega t), as in `ground_effect.cpp`, whose admittance it
  takes. When speed became the question the answer was again not GFPE, which would have
  brought that pole back and splits refraction from diffraction (a step limit exactly where the
  gradient is strong), but a longer step for the same grid: the split-step Padé march.
- **The default march is split-step Padé** (`PeScheme::split_step_pade`):
  psi(r + dr) = [c_0 + sum_j c_j (Q - p_j)^-1] psi(r), a rational stand-in for
  exp(i k_a dr (sqrt(1 + Q) - 1)). One tridiagonal solve per term, about ten terms for an
  eight-wavelength step, where Crank–Nicolson takes eighty steps: 39x per march at the
  defaults. Timings are under [Crank–Nicolson](#cranknicolson-the-independent-check).
- **The coefficients are fitted at run time, not tabulated and not Collins' Taylor Padé.** The
  step has to divide `output_step`, so k_a dr is different for every frequency. r is the type
  (n, n) interpolant of the exact propagator at Chebyshev nodes on the interval of x the run
  needs: elevation to `max_angle` (50 degrees), shifted by the profile's range of k(z), plus
  room for the surface wave. A Taylor Padé is exact at x = 0 and degrades outward, which is the
  wrong shape for that.
- **Two of the nodes sit among the evanescent modes (x = -1.5 and -6), and they are not
  optional.** Interpolating a function of modulus one gives an all-pass r, |r| = 1 along the
  whole axis, whose phase is right only on the interval. Steeper modes then travel upward at a
  twentieth of their true speed, linger, and come back as several dB of noise (3–5 dB rms on the
  rigid still-air test, which is how it was found). Pinned out there, r falls away below the
  interval and they die in a few steps. Two is the number: three or four overshoot, |r| > 1 by
  1e-4 to more than 1 between the interval and the first node. This is Collins' stability
  constraint.
- **Nothing in the construction guarantees stability, so it is measured.** Every fit is
  checked on dense samples: error on the interval, |r| <= 1 over the whole spectrum including
  the absorbing layer's imaginary parts, every pole below the axis. The fewest terms that pass
  are used; if none up to 16 does, the step is shortened. Do not remove the check to save the
  milliseconds.
- **Results are never interpolated in range.** The step is `output_step` over a whole number,
  so the march lands on every output. Over eight wavelengths psi turns through radians and
  linear interpolation is meaningless.
- **A path steeper than `max_angle` from the source is NaN**, because those modes are absorbed
  and the number would be hundreds of dB down. Crank–Nicolson reported a plausible-looking level
  there with the wrong phase. `max_angle` stops at 60 degrees (the interval is clamped at
  x = -0.75, short of the branch point at -1).
- **Wide-angle accuracy took three things, and the propagator was the least of them.** With
  `max_angle` alone raised from 35 to 50, the error against the exact two-ray field at
  30–50 degrees did not move, and doubling the grid with the Gaussian starter barely moved it
  either. What it took:
  - *The self-starter* (`PeStarter::self`, Collins 1992): the first step is the propagator times
    (1 + Q)^(-1/4), applied to a delta function shared between grid points,
    psi = sqrt(2 pi i / k_a) (1 + Q)^(-1/4) exp(...) delta. It is fitted and checked exactly as the
    propagator is (`PadeRequest::starter`), except that its modulus is bounded by twice its own
    value at the interval's low end rather than by one: it rises toward the branch point and is
    applied once.
  - *Fourth order in height* (`vertical_order` 4, Numerov). (Q - p) y = psi becomes
    [a delta^2 + (1 + delta^2/12) D] y = (1 + delta^2/12) psi; in y' = c y, with c = a + D/12 the
    column's off-diagonal, the off-diagonals are one again and the one-number-a-row
    factorization still serves. It costs a weight per term per row and a third more time. Second
    order at 20 points a wavelength is a radian of phase out by 35 degrees at 1 kHz from a 100 m
    source.
  - *The grid's edges at fourth order* ([below](#grid-and-accuracy)).
- **The starter's first step is taken in parts if its fit is refused.** Its fit is harder than
  the propagator's (the amplitude factor rises toward the branch point) and fails at long steps
  on the coarse grid. A silent fallback to the Gaussian showed as a 0.007 dB offset that
  depended on the range step, so the starter is taken over a half or a quarter of the step and
  that length's own propagator over the rest. `a long first step is still self-started` pins
  it.
- **The Gaussian starter is biased for any source where k differs from k at the ground**,
  because it is normalized with k_a, and reads high for a source aloft in a refracting column
  (by how much: [FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions)). Do not
  take an offset of that size against a stored run made with the Gaussian for a regression. The
  Gaussian
  runs only for a source within a grid point of the ground on the second-order ground row
  (which has no psi_0 to share the delta with), or when the starter's fit is refused even at a
  quarter of the step; `PeResult::self_started` says which ran.
- **Normalization:** Salomons' second-order Gaussian starter times sqrt(i k_a) makes
  p = psi exp(i k_a r)/sqrt(r) a unit point source. `tests/test_pe.cpp` holds it to the exact
  two-ray field (rigid) and to Weyl–van der Pol through NICEOPS's own
  `spherical_reflection_coefficient` (grass), which checks the ground condition's sign too.
- **The kernel pairs terms into two-lane vectors by hand** (`Lanes`, as `Pair` in
  `sphere_database.cpp`), with the term count a template parameter. From scalar source the
  compiler ran one term at a time: 11.7 ns a row against 5.9 (second order, eight terms). Each
  solve stores one complex number per row (e_r = a / pivot; the off-diagonals are all a, so the
  forward and back coefficients coincide), which halves what a step streams through the cache.
- **Not twisted.** Eight interleaved terms are throughput-bound (16 cycles of arithmetic a row
  against 16 of latency), so Crank–Nicolson's two-ended factorization would buy nothing here.

## Grid and accuracy

- **The ground condition and the interpolation are fourth order too** (`boundary_order` 4),
  and that is what makes ten points a wavelength the default. With the operator alone at fourth
  order the grid still converged at second (the convergence table is in
  [FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions)).

  How:
  - *The ground row comes from the equation, not from a wider stencil*, which would leave the
    tridiagonal. With F = y'' = k_a^2 (psi - D y) and y' = lambda y at z = 0,
    y_1 - (1 + dz lambda) y_0 = dz^2 (7 F_0 + 6 F_1 - F_2) / 24 + O(dz^5). It is a row in y_0,
    y_1, y_2 that is linear in the pole, like Numerov's rows, so every solve is still a
    resolvent of one operator, on the grid that includes psi_0.
  - *psi_0 is therefore carried, as `psi[m]`*, not derived from psi_1 and psi_2. Each step
    rebuilds it from the two rows above and the old field, a few operations per term; the
    screens turn it with the rest; the Gaussian starter sets it.
  - *Row 1's right side differs per term* (it loses c_0 H R, R = 7 psi_0 + 6 psi_1 - psi_2), so
    the forward sweep takes row 0 a term at a time before the paired loop.
  - *Cubic interpolation at both ends of the path.* The receivers read four nodes, and the
    self-starter's delta function is spread over four with the same weights: by reciprocity a
    source placed to fourth order is a receiver read to fourth order. Within three nodes of the
    ground the scheme is not symmetric, and `place_ground_source` (pe.cpp) divides out each
    mode's exact ratio there instead, so a source at any height, node 0 included, self-starts
    (see [Reciprocity](#reciprocity-one-march-per-receiver-height)).
- **`points_per_wavelength` 0 means 10 where all of that holds and 20 where anything is second
  order** (Crank–Nicolson, `vertical_order` 2, `boundary_order` 2). Do not run a second-order
  configuration at 10: it converges an order of magnitude worse
  ([FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions)).
- **The result does not depend on the step**
  ([FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions)), past the first few
  steps. Within them (the modes beyond `max_angle` still
  leaving) it is 0.1–0.4 dB. `tests/test_pade_propagator.cpp` holds that, agreement with
  Crank–Nicolson at 40 steps a wavelength (on its own second-order grid and Gaussian starter),
  the exact field to 50 degrees, the elevated source's strength, and fourth-order convergence in
  height (16x from doubling the grid, against 4x with `boundary_order` 2).
- **The domain margin is at least 250 m, and the absorbing layer 100 wavelengths.** In a shadow
  the near-ground field is a weak remainder while energy leaks upward past it. With
  60 wavelengths of margin (25 m at 1 kHz), reflection off the layer set a floor: the PE came
  out 3.8 dB high at -44 dB against the residue series.
- **Results land on exact multiples of `output_step`.** Crank–Nicolson's range step follows the
  wavelength, so it interpolates between steps; without that each frequency had its own grid,
  and band averaging read past the end of one. The split-step lands on them.
- **The band average takes the midpoints of `samples` equal slices of the third-octave band**
  (`solve_pe_bands`), one sample being the center. Sample counts must be positive: zero or
  negative counts are refused by the solver API, an empty batch included. Placing the first
  and last sample on the band edges with a full 1/N weight each read off the true band mean
  ([FIDELITY](../FIDELITY.md#shadow-transition-calibration)), and the stored PE tables and the
  shadow-edge constants were fitted against that. They stand within their
  tolerances (every PE-pinned test passes), but anything fitted to PE band levels should be
  re-scored against fresh PE runs (`--no_pe` re-scores the ray models against stored ones
  without re-solving).
- **`solve_pe_bands` pools every (band, frequency) march, highest frequency first**, because
  cost goes as f^2 and a band at a time left the 2 kHz band on a 2-thread share. A five-band
  `pe_compare` case went from 194 s to 14 s that way under Crank–Nicolson. Summing is in
  frequency order afterwards, so any thread count gives the same bits, under either scheme.

## Crank–Nicolson, the independent check

`PeScheme::crank_nicolson` (`--scheme cn`) is the march the split-step replaced, kept as an
independent check: Padé(1,1) at ten steps a wavelength, good to about 30 degrees. How closely
it agrees with the split-step, and its own error against the exact still-air field, are in
[FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions). PE runs stored before the
split-step reproduce to
their printed digits under `--scheme cn`.

- **Its left-hand tridiagonal matrix is factored once per run, and the sweep is shaped for
  latency** (2.7 ns per point per step, from 9.0). The Thomas sweep is a recurrence; each row
  waited on two dependent complex multiplies. Three things, all needed:
  - The factorization is twisted (burn at both ends): rows below the middle eliminate
    downward, rows above upward. Each sweep then carries two independent chains.
  - The pivot is folded into the stored coefficients, so each row is one complex multiply-add.
  - That multiply-add is written as two explicit `std::fma` per part (`minus_mul`). Left to
    itself the compiler contracts inside the multiply and leaves a subtract at the end, which
    puts 11 cycles of latency on the chain where 8 will do; that alone was 5.1 -> 2.7 ns.
    Only where the target has a fused multiply-add (`__FMA__`, AArch64): on x86-64 without
    `-mfma`, which is how CI's GCC builds, each `std::fma` is a call into libm on the chain, so
    there it is the plain products, and the results differ from an FMA build by round-off.

  Output matches the plain Thomas sweep to about 1e-9 dB (the same arithmetic reordered).
- **Crank–Nicolson's grid is 20 points per wavelength in height and 10 in range.** Ten in
  height was too coarse for a receiver 1.2 m up; 20 costs twice as much
  ([FIDELITY](../FIDELITY.md#parabolic-equation-against-exact-solutions)). Rerunning the stored
  calibration cases at 20 did not move the shadow-edge calibration fitted against them, so the
  fits stand. The turbulent mean-field test
  runs at 10 in height with 1000 realizations; at 20 its tolerance was too tight for its sample
  size. Halving Crank–Nicolson's range steps costs most cases only 0.02–0.03 dB, but fails
  `test_benchmark` by 0.2–0.4 dB at long range. Do not coarsen its range step to buy speed;
  that is what the split-step is for.
- **Cost.** PE time for a five-band `pe_compare` case (125 Hz–2 kHz, 100 m source to 1.5 km,
  14 cores; the ray models' columns take another 0.4 s):

  | configuration | PE time |
  | --- | --- |
  | `--scheme cn` | 25.4 s |
  | `--scheme pade` (default: 50°, self-starter, fourth order, 10 points) | 0.47 s |
  | ... `--max_angle 60` | 0.67 s |
  | ... `--max_angle 35` | 0.37 s |
  | ... `--ppw 20` | 1.4 s |
  | ... `--ppw 20 --max_angle 35 --vertical_order 2 --boundary_order 2 --starter gaussian` | 0.6 s |

  The benchmark tool takes about 10 s under the split-step, 2.5 min under Crank–Nicolson.

## Turbulence and scattering

- **Turbulence is phase screens** (`PeSetup::gamma_t`, `outer_scale`, `screen_spacing`,
  `seed`; `pe_compare --gamma_t`). von Kármán in the effective index with C^2 = gamma_t / 4
  (HARMONOISE's gamma_t, the form Ostashev's and HARMONOISE's coherences both reduce to),
  uniform in height, one screen per outer scale of range, drawn by FFT (`ScreenSource`). The
  screen's vertical spectrum is the 3-D spectrum integrated across the plane, so the vertical
  coherence is 3-D's. `solve_pe_band` makes each frequency its own realization. Still air is
  untouched: gamma_t 0 skips it all. `PeSetup::turbulence_z`/`turbulence_gamma` give a strength
  that varies with height (each row's phase scaled by sqrt(gamma(z)/gamma_t)), and the table
  tool runs several realizations per frequency (below).
- **Turbulence caps the step at the screen spacing.** So compare a turbulent run with a still
  one marched the same way: the field at a point whose reflected ray is steeper than
  `max_angle` depends on the step. `pe_turbulent_mean_field_decays` takes its still field with
  `gamma_t = 1e-30` for that reason; with 0 the two differed by 8% before any turbulence.
- **`tests/test_pe_turbulence.cpp` checks:**
  - the screens' variance, 0.782 k^2 dx (gamma_t/4) L0^(5/3), to 0.06%, less the 4% the finite
    column leaves out;
  - the ensemble-mean field's exp(-sigma^2/2) decay with range, within the sampling noise of
    300 realizations, which checks the screens and the march together;
  - single scattering: three times the strength puts 10 lg 3 more into a deep shadow;
  - `ShadowScattering` against the PE.
- **`ShadowScattering` (`--scattering`) is fitted to it**, a + 10 lg gamma_t + b lg f +
  c lg R + d lg z_s dB re free field, added in energy to the transition's field where it runs,
  in place of the floor ([THEORY §8.1](../THEORY.md) has the coefficients and clamps). Single
  scattering is linear in gamma_t, and the PE agrees. Two things to know before trusting a
  strength:
  - **The pair fit's gamma_t is far too strong for scattering**, and the measured SEL wants
    about a tenth of it ([FIDELITY](../FIDELITY.md#turbulent-scattering-strength)). It is off
    by default for that reason.
  - The 1.2 m constants come from one source height; the height term is the flush fit's.
- **Without `--gamma_t` the deep shadow, 40–90 dB down, is not physical**: with no turbulence,
  nothing scatters energy in. Absorption is optional (`PeSetup::absorption`, Np/m, the
  imaginary part of k so it acts along the actual path) and off in `pe_compare`, because the
  ray-model comparison leaves it out of all four models.

## Exact solutions and the benchmark

- **The shadow is checked against an exact solution** (`tests/test_pe_shadow.cpp`). For
  k^2 = k0^2 (1 + 2z/R) over rigid ground the field is a residue series of Airy creeping waves.
  The complex Airy function it needs is `src/airy.{hpp,cpp}`: series for |z| <= 5, asymptotic
  beyond, with the connection formula near the negative axis. Five is where both are about
  1e-7; near the positive real axis the series cancels exp(4/3 z^1.5) of its digits away. The
  reference values came from Bessel-integral quadrature, not memory, which had Ai(-10) wrong.
- **`airy_asymptotic` and `hankel1_0`'s large-argument branch are in real arithmetic**, for
  the reason `faddeeva` is (the library's complex divide, pow, exp and abs are calls with
  NaN-recovery branches that cannot fire here): `z^(3/2)` and `z^(1/4)` are `z sqrt(z)` and
  `sqrt(sqrt(z))` on the principal branch, the divides are by a reciprocal formed once, and the
  series stop on squared magnitudes. 2–3x per call (H0 209 -> 67 ns, Ai 151 -> 73), and a
  refraction run with a point receiver, whose ShadowTransition table is most of it, went from
  4.9 to 1.6 s. They agree with what they replaced to 1e-16 (H0; 2e-13 where exp(-Im z)
  underflows) and 4e-13 (Ai, where the new one is the closer to a long-double evaluation).
- **The benchmark** (`src/exact_solutions.{hpp,cpp}`, `tools/benchmark.cpp`,
  `tests/test_benchmark.cpp`) is Attenborough et al., JASA 97, 173 (1995); the PDF is not in
  the repository. Its exact solutions are direct plus image ray (case 1), the downward
  normal-mode series (case 2) and the upward residue series (case 3), all over impedance ground
  with absorption as Im k0.
- **The mode sums are derived in `exact_solutions.cpp`, not copied, because the paper's
  eq. (11) is wrong** (the corrected condition and its effect are in
  [THEORY §9.2](../THEORY.md)). Eq. (10) is right up to an overall sign. Over rigid ground the
  two forms coincide, so `test_pe_shadow` could never have caught it. Do not "fix" the code back
  to the paper.
- **Case 2's exact profile, c0/sqrt(1 - 2z/R), turns every ray below z = R/2.** Near-vertical
  ones fall back within R = 3.4 km. The PE cannot carry them (30–45° and an absorbing top), and
  neither can the benchmark's linear profile. `downward_pressure`'s `max_elevation` drops modes
  steeper than the cut. The test uses 75°: paths reaching 1 km leave below 8.4° or above 81.6°.
  At 10 Hz only about 60 modes lie below 75°, too few for a clean cut within 6 wavelengths, so
  the test uses the full series there. At 1 kHz past 500 m the PE's own grid is the limit,
  0.8 dB at 10 points a wavelength and 0.3 at 20.
- **`impedance_airy_zeros` is built for |q| <= 6 and throws past it.** It follows each zero of
  Ai' from q = 0, from a first-order seed once |a'| > 4|q|² + 10, and checks the count with the
  argument principle on a box, which caught a zero lost to the continuation in development. The
  benchmark reaches |q| = 4.9. At |q| = 8.5 the continuation loses its way, and at 10.8 the
  box's edge overflows a double. The low zero at 1 kHz is the ground's surface wave, near q²,
  and the continuation reaches it. `first_impedance_airy_zero` follows the first zero alone,
  which stays well apart from its neighbors, to any |q| (checked to 1e5 against
  a_1 - 1/q - a_1/(3q^3)): what the shadow loss needs. **Past |q| = 6 with Re q <= 0** (every
  upward problem over ordinary ground) the zeros are seeded instead from
  zeta = (n - 3/4) pi - atan(q / sqrt(x)), b = -x, and polished by Newton; checked to step by pi
  in zeta. That removed the limit for `upward_pressure`. Re q > 0 past 6 (downward, and
  surface-wave grounds) still throws.
- **A mode sum returns NaN where it cannot be trusted:** terms canceling by more than nine
  orders (the residue series near the source, whose terms grow before decaying), or a tail not
  yet dead at the 12000-mode cap. Callers skip NaN rather than tighten the cap.
- **`hankel1_0` switches from the power series to the asymptotic expansion at
  |z| = 12 - 3 Im z / |z|.** The series loses exp(|z| + Im z) of its digits, and the expansion's
  smallest term is exp(-2|z|). A switch on |z| alone put H0(3+6i) off by 1e-6. Its test values
  are Mathematica's, as are the root-finder's.
- **The benchmark's inputs are its own:** the tabulated impedances (Z = 38.79+38.41i,
  12.81+11.62i, 5.96+2.46i), not a NICEOPS impedance model, and the absorption it quotes
  (0.18 and 5.4 dB/km, ANSI S1.26-1978; ISO 9613-1 gives 0.22 and 4.98). Its profiles stop at
  1000 m, the height of its contour plots. Only case 2 at 10 km cares, and the paper notes the
  same dependence.
- **The 10 Hz benchmark legs start at 50 m, not 20, and case 2's 100 Hz leg at 25.** The first
  is inside a wavelength and a half of the source: the converged march is 0.29 dB over the exact
  field at 25 m, and Crank–Nicolson's own error happened to bring it inside 0.2. The second runs
  at ten points a wavelength and closes on the series as 1/r, 0.100 dB under at 20 m against a
  0.1 margin. No other tolerance moved.

## `niceops_pe_compare`

It runs one path through straight rays, `linear`, `stratified` and the PE (and
`resolved_<f>` under `--shadow_model resolved`, with `--tangent_decay`,
`--fit_linear_boundary` and `--no_boundary_law` as sensitivities; see
[refraction.md](refraction.md#resolved-shadow---shadow_model-resolved)). What the comparison
shows is in [FIDELITY](../FIDELITY.md#refraction-models-against-the-pe).

- **`--no_pe`** gives the ray models' columns alone, on the PE's range grid, for re-scoring a
  model change against stored PE runs without re-solving them.
- **Two things that looked like the upwind treatment and are not:**
  - *Clipping the profile below a fixed height.* Once the height passes the microphone, the
    microphone sits in the clipped layer and shadows at every range.
  - *Smoothing the profile per band.* A linear profile smoothed only over its lowest meters
    moved the lapse control by 20 dB at 2.5 km.

  The ray model is too sensitive to structure below the microphone for either to work; the
  per-band boundary law is the treatment
  ([refraction.md](refraction.md#shadowtransition-the-level-near-a-shadow-boundary)).

## Propagation tables

`--propagation_table`: a table of each band's level re free field at the straight-line distance,
by (bearing, source height, receiver height, range, band), for one stratified atmosphere over
uniform ground. `niceops_propagation_table` builds it from the PE or from the straight-ray ground
model; the format is the specification, so other methods can write it. The composition is in
[THEORY §9.3](../THEORY.md); interpolation and reciprocity errors are in
[FIDELITY](../FIDELITY.md#propagation-tables).

- **It enters through the existing seam.** A covered pair is `spectral_level_along` on the
  straight ray with no ground model and `extra_loss` = minus the table: sphere by the
  straight-line direction, spreading over the straight-line range, absorption as the run applies
  it (zeros if the table `includes_absorption`). Nothing else in the evaluation knows the table
  exists.
- **Uncovered pairs fall back, per pair, to the run's own model**: `--refraction` if on, else
  `--ground_effect` if on, else free field. A pair is uncovered if any corner it would
  interpolate is NaN (the PE past its angle), or it is outside an axis. Zero-weight corners are
  not read, so a point on the grid beside a NaN is still covered. Microphone heights outside the
  table are refused up front, not dropped to the fallback.
- **The peak search's bound takes max(ground model's gain, table's largest value) per band**
  (`path_gains` in the bound setup), since a pair may take either.
  `propagation_table_windowed_is_exact` holds it; the PE table was checked the same way by hand.
- **`lookup` finds an even axis's cell by arithmetic and corrects it against the knots**
  (`AxisLocator`), so it is the cell a binary search gives, bit for bit; the fmod on the bearing
  is skipped inside [0, 360). -20% CPU on a table run, output byte-identical.
  `lookup finds the same cell on even and uneven axes` sweeps knots, midpoints and a hair either
  side of them against a plain search. The four-band lookup and the rejected height cache are in
  [performance.md](performance.md#measured-and-rejected-optimizations).
- **`BandMap` is per run, not in the table**, because the table is loaded once
  (`RunOptions::propagation_table`, shared) and each trajectory can name its own database.
- **Bearings are compass, toward.** The run turns the local separation back to east/north by
  inverting `CompassToLocal` at `--frame_bearing` (90, +x east, when unset).
- **Microphone height comes from `--receiver_height` under a table**, as under
  `--ground_effect`; a run without either once put microphones on the ground.
- **The plumbing check is the ground-model table** (`propagation_table_is_the_ground_model`).
  It reproduces `--ground_effect --exact_ground` to the table's interpolation, and the error
  falls fourfold per halving of the height spacing, which is what says it is interpolation.
- **Plates read the table at its lowest receiver height plus K**, the plate over bare ground
  there per band, azimuth and elevation on a 1 km sphere (`table_plate_offset` in niceops.cpp),
  as ShadowTransition places a plate, bilinear in azimuth and elevation. Things to keep:
  - **K depends on elevation.** It was taken at 2 degrees for every pair; through the ground-model
    table's flush row that put the demo footprint up to 0.34 dB off
    ([FIDELITY](../FIDELITY.md#propagation-tables)).
  - **A raised lowest row is warned about.** Its error grows with the row's height (0.13 dB at
    1.2 m, 1.3 dB at 30 m) and is not the grid: a twelve times finer K grid leaves it unchanged.
  - **The elevation step follows the bare term's fringes**: eight steps per λ/(2h) in the sine of
    the elevation at the highest band, clamped to 0.05–0.5 degrees (0.5 for a flush row, about
    900 elevations for 1.2 m at 10 kHz).
  - **The bound multiplies the table's largest gain by the largest K node**, which bounds the
    interpolant. No current CTest needs it: through a ground-model table, table plus K is the
    plate's level, which the plate's own gain bound covers. A PE table, whose fringes need not
    line up with bare ground's, can need it.
  - **Flush microphones are not checked against the table's receiver heights**, since they read
    its lowest row. Checking them refused every plate run whose table had no row at 0; the
    warning above replaces that. **With `--refraction` the table takes the pairs it covers and
  the ray model the rest** (`table_level` first in `refracted_level`); that run evaluates every
  state, since its bounds were not built for both geometries.
- **Turbulence and roughness are baked into the table, not applied on it.** A ground-model table
  takes `--turbulence`, `--gamma_t` (HARMONOISE) and `--roughness` into the ground model it
  tabulates, and records them (`record_coherence`: `turbulence_model`, `turbulence` or `gamma_t`,
  `roughness_m`, at full precision); a turbulent PE table records `phase_screens` and its
  `gamma_t` the same way. The run reads them back (`table_coherence`) and warns when its own
  settings differ, since table pairs ignore them.
  Applying the run's coherence on top of a table would need the direct and reflected fields
  apart (a still PE table holds only their sum, a turbulent one has it already), which is a new
  file format; building the table with the run's settings needs none, and a ground-model table
  so built reproduces `--ground_effect` with them (`propagation_table_carries_turbulence`,
  `propagation_table_carries_harmonoise`). A table that records nothing is still air, which every
  table built before the record was unless it says `gamma_t` with no `turbulence_model`; that one
  wrote it to six decimals, so it is compared at that precision, and "0.000000" counts as
  turbulent. The check runs either
  way round: a still run on a turbulent table is warned too, since its uncovered pairs and the
  plate offset's bare ground are still (`a run's turbulence and roughness are held to its
  table's, either way round`).
- **Thread count is invisible.** `-j` (or `NICEOPS_THREADS`, which only the PE tools read) gives
  the same bits at any count, turbulent tables included, because each march's seeds are fixed
  before the work is dealt out (`pe_tools_thread_count_is_invisible`).
- **Over terrain the table is read in each pair's mean-plane frame**, from the same
  `terrain_ray` path the ground model takes (`table_read`; the rules are in
  [ground.md](ground.md#things-to-preserve)). Flush microphones read the lowest row at the
  source's height above the microphone, and the K grid is unchanged, because it is the plate in
  the plane's own frame. A refracting table so read assumes a terrain-following stratification;
  `table_refracts` is anything but a ground-model table, unless its provenance records the
  atmosphere `none`, and only such a table warns past 5 degrees of tilt.
- **Not done:** a sphere lookup along a refracted launch direction (the table cannot say what it
  is).

### Reciprocity: one march per receiver height

`src/pe_table.{hpp,cpp}`. `--march source` marches once per bearing and source height and reads
every receiver height. `--march receiver`, the default, marches once per bearing and receiver
height with the source at the receiver's height, and reads every source height. In the
effective-sound-speed model the field from z_s to z_r along a bearing equals the field from z_r
to z_s along the reverse bearing with the wind reversed, and the effective profile c + w.h is
the same, so it is the same march with the heights exchanged. A table has a couple of receiver
heights and many source heights, so it is several times fewer marches; the domain, which reaches
250 m past the highest height a march reads, is why the speed-up is about 2.5x rather than the
march count's ratio.

- **On the grid this holds only if a source near the ground is placed as the adjoint of the
  receivers' read.** The fourth-order ground row is a half cell whose weights differ from
  Numerov's rows above it, so a source within three grid points of the ground came out wrong
  (one at the first node 1.34 dB high and one at the second 0.37 dB low against the exact
  two-ray field), and a source on the ground (a flush microphone, in a reciprocal march) had no
  self-starter at all. `place_ground_source` divides each of those nodes by its exact weight, a
  function of the mode's vertical wavenumber and of the ground's admittance applied as an
  operator; pe.cpp has the derivation. A source at 0 to 3.5 grid spacings then agrees with the
  receiver read at the same height to 0.003 dB at grazing and within 0.06 dB at 45°, and with
  the exact field over rigid ground and grass alike (`tests/test_pe.cpp`). Sources from four
  spacings up are placed as they always were. Tables built before this differ slightly; see
  [HISTORY](../HISTORY.md#changes-that-affect-results).
- **A turbulent table is reciprocal only in its ensemble mean**, not realization by
  realization, so the two marches agree to sampling noise. Every march has its own seeds.

### Turbulent tables

`--gamma_t`, `--outer_scale`, `--seed` and `--realizations` march each frequency sample through
phase screens `--realizations` times with distinct seeds; each cell is the energy mean over the
band's samples and the realizations. `--gamma_profile FILE` (CSV of z in `--units` and gamma in
m^-2/3, linear between, held outside) gives a strength that varies with height: the screens are
drawn at `--gamma_t` (default the profile's largest value) and each row's phase is scaled by
sqrt(gamma(z)/gamma_t). That approximates a field whose statistics change with height: the
vertical spectrum and the outer scale stay a homogeneous field's, where near the ground the
eddies shrink with height. `--gamma_t 0` gives the still table bit for bit whatever the other
options say. `tests/test_pe_table.cpp` checks the screens' variance per height, that a profile
scales the march's screens by its square root (mean-field decay), that three times the strength
puts 10 lg 3 more into a deep shadow, and that the realizations average as stated.

### Table file format

The file is netCDF-4, in SI whatever `--units` says. `src/propagation_table.hpp` is the
specification; `PropagationTable::load` refuses a file that does not match it.

| Item | Kind | Content |
| --- | --- | --- |
| `bearing` | coordinate variable on its own dimension, `units` = degrees | compass bearing the sound travels toward, uniform from 0 over 360; one value means isotropic |
| `source_height` | coordinate, m | above the ground, finite, ascending |
| `receiver_height` | coordinate, m | above the ground, finite, ascending |
| `range` | coordinate, m | horizontal, finite, ascending |
| `band` | coordinate, Hz | band centers, ascending; the run's bands within `--band_limits` are matched to 0.5% and must all be there |
| `relative_level` | double, dimensioned (bearing, source_height, receiver_height, range, band) | dB re free field at the straight-line distance; NaN where the method has no value |
| `propagation_table_version` | global, one integer | 1 (`PropagationTable::format_version`); a missing or different version is refused |
| `includes_absorption` | global, one integer | nonzero when the values already include atmospheric absorption |
| `method` | global text | how it was built (`pe`, `ground_model`, or another writer's name) |
| `bearing_convention` | global text | written for readers; not read |
| `provenance_*` | global text | free-form build record, read back with the table; `provenance_turbulence_model`, `provenance_turbulence`, `provenance_gamma_t` and `provenance_roughness_m` are read as the table's turbulence and roughness (`table_coherence`) |

Values are interpolated linearly in dB along each axis.
