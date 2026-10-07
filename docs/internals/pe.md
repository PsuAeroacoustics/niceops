# PE, the benchmark's exact solutions, and propagation tables

- `src/pe.{hpp,cpp}`, `tools/pe_compare.cpp` — roadmap step 7a: a wide-angle parabolic
  equation as the full-wave reference, and `niceops_pe_compare`, which runs one path through
  straight, linear, stratified and PE. Load-bearing:
  - **A finite-difference PE, not the GFPE the roadmap named.** The impedance condition
    dpsi/dz = -i k beta psi goes straight into the grid (psi_0 eliminated by a one-sided
    difference), so no surface-wave pole has to be handled by hand. The convention is
    exp(-i omega t), as in `ground_effect.cpp`, whose admittance it takes. When speed became
    the question the answer was again not GFPE, which would have brought that pole back and
    splits refraction from diffraction (a step limit exactly where the gradient is strong),
    but a longer step for the same grid: the split-step Padé march below.
  - **The default march is split-step Padé** (`PeScheme::split_step_pade`,
    `src/pade_propagator.{hpp,cpp}`): psi(r + dr) = [c_0 + sum_j c_j (Q - p_j)^-1] psi(r), a
    rational stand-in for exp(i k_a dr (sqrt(1 + Q) - 1)). One tridiagonal solve per term,
    about ten terms for an eight-wavelength step, where Crank–Nicolson takes eighty
    steps. 39x per march at the defaults; a five-band `pe_compare` case's PE 25.4 s ->
    0.47 s, the benchmark tool 2.5 min -> 10 s. Things to keep:
    - **The coefficients are fitted at run time, not tabulated and not Collins' Taylor
      Padé.** The step has to divide `output_step`, so k_a dr is different for every
      frequency. r is the type (n, n) interpolant of the exact propagator at Chebyshev nodes
      on the interval of x the run needs: elevation to `max_angle` (50 degrees), shifted by
      the profile's range of k(z), plus room for the surface wave. A Taylor Padé is exact
      at x = 0 and degrades outward, which is the wrong shape for that.
    - **Two of the nodes sit among the evanescent modes (x = -1.5 and -6), and they are
      not optional.** Interpolating a function of modulus one gives an all-pass r, |r| = 1
      along the whole axis, whose phase is right only on the interval. Steeper modes then
      travel upward at a twentieth of their true speed, linger, and come back as several dB
      of noise (3-5 dB rms on the rigid still-air test, which is how it was found). Pinned
      out there r falls away below the interval and they die in a few steps. Two is the
      number: three or four overshoot, |r| > 1 by 1e-4 to more than 1 between the interval
      and the first node. This is Collins' stability constraint.
    - **Nothing in the construction guarantees stability, so it is measured.** Every fit is
      checked on dense samples: error on the interval, |r| <= 1 over the whole spectrum
      including the absorbing layer's imaginary parts, every pole below the axis. The fewest
      terms that pass are used; if none up to 16 does, the step is shortened. Do not remove
      the check to save the milliseconds.
    - **Results are never interpolated in range.** The step is `output_step` over a whole
      number, so the march lands on every output. Over eight wavelengths psi turns through
      radians and linear interpolation is meaningless.
    - **A path steeper than `max_angle` from the source is NaN**, because those modes are
      absorbed and the number would be hundreds of dB down. Crank–Nicolson reported a
      plausible-looking level there with the wrong phase. `max_angle` stops at 60 degrees
      (the interval is clamped at x = -0.75, short of the branch point at -1).
    - **Wide-angle accuracy took three things, and the propagator was the least of
      them.** With `max_angle` alone raised from 35 to 50 the error against the exact
      two-ray field at 30-50 degrees did not move (0.6-1.7 dB rms at 1 kHz). Doubling the
      grid with the Gaussian starter barely moved it either. What it took:
      - *The self-starter* (`PeStarter::self`, Collins 1992): the first step is the
        propagator times (1 + Q)^(-1/4), applied to a delta function shared between two
        grid points, psi = sqrt(2 pi i / k_a) (1 + Q)^(-1/4) exp(...) delta. It is fitted
        and checked exactly as the propagator is (`PadeRequest::starter`), except that
        its modulus is bounded by twice its own value at the interval's low end rather
        than by one: it rises toward the branch point and is applied once.
      - *Fourth order in height* (`vertical_order` 4, Numerov). (Q - p) y = psi becomes
        [a delta^2 + (1 + delta^2/12) D] y = (1 + delta^2/12) psi; in y' = c y, with
        c = a + D/12 the column's off-diagonal, the off-diagonals are one again and the
        one-number-a-row factorization still serves. It costs a weight per term per row
        and a third more time. Second order at 20 points a wavelength is a radian of
        phase out by 35 degrees at 1 kHz from a 100 m source.
      - Together, and with the grid's edges at fourth order (next): 0.07 dB rms or
        better at every elevation to 60 degrees at ten points a wavelength, 0.01 to 50 at
        twenty.
    - **The ground condition and the interpolation are fourth order too**
      (`boundary_order` 4), and that is what makes ten points a wavelength the default.
      With the operator alone at fourth order the grid still converged at second. On the
      windy profile against 80 points: second order throughout 0.02-0.04 dB rms at 20
      points, fourth-order operator with second-order edges 0.01-0.02, fourth throughout
      under 0.0005, and 0.001-0.002 at 10 (0.004-0.01 at 6). How:
      - *The ground row comes from the equation, not from a wider stencil*, which would
        leave the tridiagonal. With F = y'' = k_a^2 (psi - D y) and y' = lambda y at z = 0,
        y_1 - (1 + dz lambda) y_0 = dz^2 (7 F_0 + 6 F_1 - F_2) / 24 + O(dz^5). It is a row
        in y_0, y_1, y_2 that is linear in the pole, like Numerov's rows, so every solve is
        still a resolvent of one operator -- on the grid that includes psi_0.
      - *psi_0 is therefore carried, as `psi[m]`*, not derived from psi_1 and psi_2. Each
        step rebuilds it from the two rows above and the old field, a few operations per
        term; the screens turn it with the rest; the Gaussian starter sets it.
      - *Row 1's right side differs per term* (it loses c_0 H R, R = 7 psi_0 + 6 psi_1 -
        psi_2), so the forward sweep takes row 0 a term at a time before the paired loop.
      - *Cubic interpolation at both ends of the path.* The receivers read four nodes,
        and the self-starter's delta function is spread over four with the same weights:
        by reciprocity a source placed to second order is a receiver read to second order.
    - **`points_per_wavelength` 0 means 10 where all of that holds and 20 where anything
      is second order** (Crank–Nicolson, `vertical_order` 2, `boundary_order` 2). Do not
      run a second-order configuration at 10: it is 0.06-0.19 dB rms and 0.6 at worst.
    - **The starter's first step is taken in parts if its fit is refused.** Its fit is
      harder than the propagator's (the amplitude factor rises toward the branch point)
      and fails at long steps on the coarse grid. It used to fall back to the Gaussian
      silently, which showed as a 0.007 dB offset that depended on the range step; now it
      takes the starter over a half or a quarter of the step and that length's own
      propagator over the rest. `a long first step is still self-started` pins it.
    - **The Gaussian starter is biased for any source where k differs from k at the
      ground**, because it is normalized with k_a: 0.11-0.15 dB high at every range for a
      100 m source in the benchmark's upward column (against the residue series), 0.10 dB
      on `example_windy.csv`. Every PE level from before the self-starter carries that,
      so do not take a 0.1 dB offset against a stored run for a regression. The
      self-starter is within 0.04 dB. A source within a grid point of the ground, or a
      starter fit that is refused, falls back to the Gaussian; `PeResult::self_started`
      says which ran.
    - **Compare a turbulent run with a still one marched the same way.** Turbulence caps
      the step at the screen spacing, and the field at a point whose reflected ray is
      steeper than `max_angle` depends on the step. `pe_turbulent_mean_field_decays`
      takes its still field with `gamma_t = 1e-30` for that reason; with 0 the two
      differed by 8% before any turbulence.
    - **The result does not depend on the step**: 2 to 32 wavelengths moves it by under
      0.001 dB, 60 dB into a shadow, past the first few steps. Within them (the modes
      beyond `max_angle` still leaving) it is 0.1-0.4 dB. `tests/test_pade_propagator.cpp`
      holds that, agreement with Crank–Nicolson at 40 steps a wavelength (on its own
      second-order grid and Gaussian starter), the exact field to 50 degrees, the elevated
      source's strength, and fourth-order convergence in height (16x from doubling the
      grid, against 4x with `boundary_order` 2).
    - **The kernel pairs terms into two-lane vectors by hand** (`Lanes`, as `Pair` in
      `sphere_database.cpp`), with the term count a template parameter. From scalar source
      the compiler ran one term at a time: 11.7 ns a row against 5.9 (second order, eight
      terms). Each solve stores one
      complex number per row (e_r = a / pivot; the off-diagonals are all a, so the forward
      and back coefficients coincide), which halves what a step streams through the cache.
    - **Not twisted.** Eight interleaved terms are throughput-bound (16 cycles of
      arithmetic a row against 16 of latency), so the two-ended factorization below would
      buy nothing here.
    - **Turbulence caps the step at the screen spacing.**
    - **The 10 Hz benchmark legs start at 50 m, not 20, and case 2's 100 Hz leg at 25.**
      The first is inside a wavelength and a half of the source: the converged march is
      0.29 dB over the exact field at 25 m, and Crank–Nicolson's own error happened to
      bring it inside 0.2. The second runs at ten points a wavelength and closes on the
      series as 1/r, 0.100 dB under at 20 m against a 0.1 margin. No other tolerance
      moved.
  - **`PeScheme::crank_nicolson` is the march this replaced, kept as an independent
    check**: Padé(1,1) at ten steps a wavelength, good to about 30 degrees. It
    differs from the split-step by 0.03-0.12 dB rms on `example_windy.csv` to 1.5 km, and
    against the exact still-air field the error is its own (0.7 dB rms at 25-30 degrees,
    2.1 at 30-35, where the split-step has 0.1 and 0.6). Its left-hand tridiagonal matrix
    is factored once per run. It is latency-bound, and shaped for that (2.7 ns per point
    per step, from 9.0). The Thomas sweep is a recurrence; each row waited on two dependent
    complex multiplies. Three things, all needed:
    - The factorization is twisted (burn at both ends): rows below the middle eliminate
      downward, rows above upward. Each sweep then carries two independent chains.
    - The pivot is folded into the stored coefficients, so each row is one complex
      multiply-add.
    - That multiply-add is written as two explicit `std::fma` per part (`minus_mul`).
      Left to itself the compiler contracts inside the multiply and leaves a subtract at
      the end, which puts 11 cycles of latency on the chain where 8 will do; that alone
      was 5.1 -> 2.7 ns.

    Output matches the plain Thomas sweep to ~1e-9 dB (the same arithmetic reordered), and
    PE runs stored before the split-step reproduce to their printed digits under
    `--scheme cn`. `solve_pe_bands` pools every (band, frequency) march, highest frequency
    first, because cost goes as f^2 and a band at a time left the 2 kHz band on a 2-thread
    share. A five-band `pe_compare` case went from 194 s to 14 s that way. Summing is in
    frequency order afterwards, so any thread count gives the same bits, under either
    scheme.
  - **Crank–Nicolson's grid is 20 points per wavelength in height and 10 in range.**
    (The split-step's is 10 in height; above.) Against a
    30 x 30 reference, 10 in height was ~0.1 dB rms near the ground but 0.47 dB rms (2.2 at
    worst) for a 1.2 m receiver at 1 kHz; 20 is 0.09 dB, for twice the cost. Rerunning the
    stored calibration cases at 20 moved the PE by 0.16-0.39 dB rms (1.1 at worst) and did
    not move the shadow-edge calibration fitted against it, so the fits stand. The
    turbulent mean-field test runs at 10 in height with 1000 realizations; at 20 its
    tolerance was too tight for its sample size. Halving Crank–Nicolson's range steps costs
    most cases only 0.02-0.03 dB, but fails `test_benchmark` by 0.2-0.4 dB at long range. Do
    not coarsen its range step to buy speed; that is what the split-step is for.
  - **The band average takes the midpoints of `samples` equal slices of the third-octave band**
    (`solve_pe_bands`), one sample being the center. It once put the first and last sample on the band
    edges with a full 1/N weight each, which read 0.16 to 0.47 dB rms off the true band mean downwind
    (bias +0.1 to +0.35 dB) and 0.4 to 1.5 upwind, and the stored PE tables and the shadow-edge
    constants below were fitted against that. They stand within their tolerances (every PE-pinned test
    passes), but anything fitted to PE band levels should be re-scored against fresh PE runs
    (`--no_pe` re-scores the ray models against stored ones without re-solving).
  - **Normalization:** Salomons' second-order Gaussian starter times sqrt(i k_a) makes
    p = psi exp(i k_a r)/sqrt(r) a unit point source. `tests/test_pe.cpp` holds it to the
    exact two-ray field (0.3 dB, rigid) and to Weyl–van der Pol through NICE-OPS's own
    `spherical_reflection_coefficient` (0.5 dB, grass), which checks the ground condition's
    sign too.
  - **Results land on exact multiples of `output_step`.** Crank–Nicolson's range step
    follows the wavelength, so it interpolates between steps; without that each frequency
    had its own grid, and band averaging read past the end of one. The split-step lands on
    them (above).
  - **The shadow is checked against an exact solution** (`tests/test_pe_shadow.cpp`). For
    k^2 = k0^2 (1 + 2z/R) over rigid ground the field is a residue series of Airy creeping
    waves, and the PE matches it within 0.02 dB down to -47 dB at 250 Hz to 1 kHz. The complex
    Airy function it needs is `src/airy.{hpp,cpp}`: series for |z| <= 5, asymptotic beyond,
    with the connection formula near the negative axis. Five is where both are ~1e-7; near
    the positive real axis the series cancels exp(4/3 z^1.5) of its digits away. The
    reference values came from Bessel-integral quadrature, not memory, which had Ai(-10)
    wrong. **`airy_asymptotic` and `hankel1_0`'s large-argument branch are in real arithmetic**,
    for the reason `faddeeva` is (the library's complex divide, pow, exp and abs are calls with
    NaN-recovery branches that cannot fire here): `z^(3/2)` and `z^(1/4)` are `z sqrt(z)` and
    `sqrt(sqrt(z))` on the principal branch, the divides are by a reciprocal formed once, and the
    series stop on squared magnitudes. 2-3x per call (H0 209 -> 67 ns, Ai 151 -> 73), and a
    refraction run with a point receiver, whose ShadowTransition table is most of it, went from
    4.9 to 1.6 s. They agree with what they replaced to 1e-16 (H0; 2e-13 where exp(-Im z)
    underflows) and 4e-13 (Ai, where the new one is the closer to a long-double evaluation).
  - **The domain margin is at least 250 m, and the absorbing layer 100 wavelengths.** In a
    shadow the near-ground field is a weak remainder while energy leaks upward past it.
    With 60 wavelengths of margin (25 m at 1 kHz), reflection off the layer set a floor: the
    PE came out 3.8 dB high at -44 dB against the residue series.
  - **Turbulence is phase screens** (`PeSetup::gamma_t`, `outer_scale`, `screen_spacing`,
    `seed`; `pe_compare --gamma_t`). von Karman in the effective index with C^2 = gamma_t / 4
    -- HARMONOISE's gamma_t, the form Ostashev's and HARMONOISE's coherences both reduce to
    -- uniform in height, one screen per outer scale of range, drawn by FFT (`ScreenSource`).
    The screen's vertical spectrum is the 3-D spectrum integrated across the plane, so the
    vertical coherence is 3-D's. `solve_pe_band` makes each frequency its own realization.
    `tests/test_pe_turbulence.cpp` checks the screens' variance (0.06%), the mean field's
    exp(-sigma^2/2) decay (within sampling noise), and ShadowScattering against the PE.
    Still air is untouched: gamma_t 0 skips it all.
  - **`ShadowScattering` (`--scattering`) is fitted to it**, a + 10 lg gamma_t + b lg f +
    c lg R + d lg z_s dB re free field, added in energy to the transition's field where it
    runs, with the floor off. Single scattering is linear in gamma_t and the PE agrees (x3 ->
    4.82 dB, theory 4.77). Two things to know before trusting a strength:
    - **The pair fit's gamma_t = 3e-5 is far too strong for scattering.** It fills the
      day-284 shadows to -3 dB re free field; the measured SEL wants <= ~3e-6
      (HARMONOISE's "moderate"), and does as well with the floor. It is off by default for
      that reason.
    - The 1.2 m constants come from one source height; the height term is the flush fit's.
  - Without `--gamma_t` the deep shadow, 40–90 dB down, is not physical. Absorption is optional
    (`PeSetup::absorption`, Np/m, the imaginary part of k so it acts along the actual path)
    and off in `pe_compare`, because the ray-model comparison leaves it out of all four.
  - **Findings** (`docs/MANUAL.md` has the tables). Downwind, `stratified` is within 0.4–2.1 dB rms of
    the PE and the others are 1–8 dB off. Upwind, the PE's shadow onset lies between
    `stratified` (too early) and `linear` (too late), moving toward `stratified` with
    frequency; the per-band boundary law is that treatment.
  - **`--no_pe`** gives the ray models' columns alone, on the PE's range grid, for
    re-scoring a model change against stored PE runs without re-solving them.
- `src/exact_solutions.{hpp,cpp}`, `tools/benchmark.cpp`, `tests/test_benchmark.cpp` — step 7b:
  the Attenborough et al. benchmark (JASA 97, 173 (1995); the PDF is not in the repo) and its
  exact solutions: direct plus image ray (case 1), the downward normal-mode series (case 2)
  and the upward residue series (case 3), all over impedance ground with absorption as
  Im k0. The PE reproduces the paper's figs. 11–14 in every panel. `docs/MANUAL.md` has the
  tables. Load-bearing:
  - **The mode sums are derived in `exact_solutions.cpp`, not copied, because the paper's
    eq. (11) is wrong.** Its upward ground condition, Ai'(b) + i k0 l beta w Ai(b) = 0, should
    read Ai'(b) - i k0 l beta w² Ai(b) = 0, a factor exp(i pi/3) in the admittance term. The
    printed form misses the PE by 1.5–12 dB at 100 Hz, and the derived one agrees to 0.02 dB.
    Eq. (10) is right up to an overall sign. Over rigid ground the two forms coincide, so
    `test_pe_shadow` could never have caught it. Do not "fix" the code back to the paper.
  - **Case 2's exact profile, c0/sqrt(1 - 2z/R), turns every ray below z = R/2.** Near-vertical
    ones fall back within R = 3.4 km like mortar rounds. The PE cannot carry them (30–45° and
    an absorbing top), and neither can the benchmark's linear profile. Uncut, they put 5 dB
    between the series and the PE by 1 km at 100 Hz, and features at R and 2R that the paper's
    own normal-mode figure shows and its FFP/PE figure lacks. `downward_pressure`'s
    `max_elevation` drops modes steeper than the cut. The test uses 75°: paths reaching 1 km
    leave below 8.4° or above 81.6°. At 10 Hz only ~60 modes lie below 75°, too few for a
    clean cut within 6 wavelengths, so the test uses the full series there. At 1 kHz past
    500 m the PE's own grid is the limit, 0.8 dB at 10 points a wavelength and 0.3 at 20.
  - **`impedance_airy_zeros` is built for |q| <= 6 and throws past it.** It follows each zero
    of Ai' from q = 0, from a first-order seed once |a'| > 4|q|² + 10. It checks the count
    with the argument principle on a box, which caught a zero lost to the continuation in
    development. The benchmark reaches |q| = 4.9. At |q| = 8.5 the continuation loses its
    way, and at 10.8 the box's edge overflows a double. The low zero at 1 kHz is the ground's
    surface wave, near q², and the continuation reaches it. `first_impedance_airy_zero`
    follows the first zero alone, which stays well apart from its neighbors, to any |q|
    (checked to 1e5 against a_1 - 1/q - a_1/(3q^3)): what the shadow loss needs.
    **Past |q| = 6 with Re q <= 0** (every upward problem over ordinary ground) the zeros
    are seeded instead from zeta = (n - 3/4) pi - atan(q / sqrt(x)), b = -x, and polished by
    Newton; checked to step by pi in zeta. That removed the limit for `upward_pressure`.
    Re q > 0 past 6 (downward, and surface-wave grounds) still throws.
  - **A mode sum returns NaN where it cannot be trusted:** terms canceling by more than nine
    orders (the residue series near the source, whose terms grow before decaying), or a tail
    not yet dead at the 12000-mode cap. Callers skip NaN rather than tighten the cap.
  - **`hankel1_0` switches from the power series to the asymptotic expansion at
    |z| = 12 - 3 Im z / |z|.** The series loses exp(|z| + Im z) of its digits, and the
    expansion's smallest term is exp(-2|z|). A switch on |z| alone put H0(3+6i) off by 1e-6.
    Its test values are Mathematica's, as are the root-finder's.
  - **The benchmark's inputs are its own:** the tabulated impedances (Z = 38.79+38.41i,
    12.81+11.62i, 5.96+2.46i), not a NICE-OPS impedance model, and the absorption it quotes
    (0.18 and 5.4 dB/km, ANSI S1.26-1978; ISO 9613-1 gives 0.22 and 4.98). Its profiles stop
    at 1000 m, the height of its contour plots. Only case 2 at 10 km cares, and the paper
    notes the same dependence. A full run is ~10 s (2.5 min under Crank–Nicolson).
- `src/propagation_table.{hpp,cpp}`, `tools/propagation_table.cpp` -- `--propagation_table`: a
  table of each band's level re free field at the straight-line distance, by (bearing, source
  height, receiver height, range, band), for one stratified atmosphere over uniform ground.
  Built by `niceops_propagation_table` from the PE or from the straight-ray ground model; the
  format (netCDF, `propagation_table_version` 1) is the specification, so other methods can
  write it. Things to keep:
  - **It enters through the existing seam.** A covered pair is `spectral_level_along` on the
    straight ray with no ground model and `extra_loss` = minus the table: sphere by the
    straight-line direction, spreading over the straight-line range, absorption as the run
    applies it (zeros if the table `includes_absorption`). Nothing else in the evaluation knows
    the table exists.
  - **Uncovered pairs fall back, per pair, to the run's own model**: `--ground_effect` if on,
    free field if not. A pair is uncovered if any corner it would interpolate is NaN (the PE
    past its angle), or it is outside an axis. Zero-weight corners are not read, so a point
    on the grid beside a NaN is still covered. Microphone heights outside the table are
    refused up front, not dropped to the fallback.
  - **The peak search's bound takes max(ground model's gain, table's largest value) per band**
    (`path_gains` in the bound setup), since a pair may take either.
    `propagation_table_windowed_is_exact` holds it; the PE table was checked the same way by
    hand.
  - **`lookup` finds an even axis's cell by arithmetic and corrects it against the knots**
    (`AxisLocator`), so it is the cell a binary search gives, bit for bit; the fmod on the
    bearing is skipped inside [0, 360). -20% CPU on a table run, output byte-identical.
    `lookup finds the same cell on even and uneven axes` sweeps knots, midpoints and a hair
    either side of them against a plain search. Caching source and receiver height cells
    in the caller was also tried, together with skipping bearings for isotropic tables.
    After the four-band optimization, a 961-observer table run showed no measurable gain
    (2.953 versus 2.952 s median acoustic time); the additional API and caches were removed.
  - **`BandMap` is per run, not in the table**, because the table is loaded once
    (`RunOptions::propagation_table`, shared) and each trajectory can name its own database.
  - **Lookup accumulates four bands together**, sharing each corner's address and weight.
    Each band's additions stay in corner order, including the initial addition to zero;
    this enables parallel arithmetic without reassociating a sum. Reordered band columns,
    partial ranges, the scalar tail, and NaN gaps are covered by the table tests. Zero-weight
    corners remain unread, and excluded output bands remain untouched.
  - **Bearings are compass, toward.** The run turns the local separation back to east/north by
    inverting `CompassToLocal` at `--frame_bearing` (90, +x east, when unset).
  - **Microphone height comes from `--receiver_height` under a table**, as under
    `--ground_effect`; before, a run without either put microphones on the ground.
  - **The plumbing check is the ground-model table.** It reproduces `--ground_effect
    --exact_ground` to the table's interpolation, 0.19 dB at worst at 10 ft spacing, 0.02 at
    5. The mean error falls fourfold per halving, which is what says it is interpolation.
  - **Plates read the table at its lowest receiver height plus K**, the plate over bare
    ground there per band and azimuth at 2 degrees (`table_plate_offset` in niceops.cpp), as
    ShadowTransition places a plate. **With `--refraction` the table takes the pairs it
    covers and the ray model the rest** (`table_level` first in `refracted_level`); that
    run evaluates every state, since its bounds were not built for both geometries.
  - **Not done:** ground elevation (flat ground at z = 0; `--terrain` refuses a table), and a
    sphere lookup along a refracted launch direction (the table cannot say what it is).
