# NICE-OPS history

How NICE-OPS got to its current form, and why old results differ from new ones. For what the
tool does now, see the [README](../README.md) and the [manual](MANUAL.md).

## Origins (2020–2021)

Eric Greenwood and Rob Rau developed NICE-OPS in D, starting in January 2020. It used a
CGAL-based barycentric interpolator, the SPLINTER B-spline library and matplotlib-d. The S-76D
case became the default in March 2020. The AW139, B205 and B407 databases and tracks from the
joint NASA/FAA flight tests were added in January 2021. They are the `_M1`, `_M2`, `_M3` and
`_M4` files in `databases/`, all broadband, with discrete load factors.

The June 2020 reference figures came from that code. They are kept in a local `figures/`
directory, which is not tracked. The C++ build reproduces them to within 1.05 dB on every cell when given the settings they were drawn with (see the
manual's [Plotting](MANUAL.md#plotting)).

## The C++17 port (September 2026)

In September 2026 NICE-OPS was ported to C++17 (`666046b`) and checked against the D reference.
On the demo case the port matched D to 1e-9 relative, with an actual maximum difference of about
1e-14, a few ULPs of LDC-vs-clang floating point. Each D dependency was replaced with a small
local equivalent. [`internals/porting.md`](internals/porting.md) lists the replacements, the four
mistranslations found during validation, and the D quirks kept on purpose.

The D sources were removed once the port was verified complete. They remain in history at
`c74bab3` on `origin/no-stk`. Check that commit out when the D oracle is needed again, for
example through the `regenerate-goldens` skill.

## Departures from D

Port fidelity is the rule: a D behavior is changed only when a measured case shows it is wrong.
Two changes have met that bar:

- **The 4x EAA multiplier.** D multiplied the database's excess-attenuation factor by 4, under a
  `//TODO Add ground effects`. Removing it raised the demo footprint by 9.93 dB mean (up to
  22.7 dB), and `golden/197421_LAE.csv` was regenerated from the C++ build.
- **The ambient mix.** D added the `-a` ambient level as pressure amplitude, which read up to
  3 dB high. It is now added in energy.

Since then the SEL golden has been regenerated once more, when spreading was capped at the
sphere surface. `golden/state_197421.csv` is still D-derived, so the trajectory state matches D.
[`internals/testing.md`](internals/testing.md) gives the provenance of each golden.

## Changes that affect results

**Results from before September 2026 are not comparable with results from after.** The
manual's [Changes that affect results](MANUAL.md#changes-that-affect-results) lists every
correction that moved NICE-OPS's numbers and how far each moved them. Besides the EAA
multiplier, it covers the DNL normalization, the LAFmax time weighting, the microphone height,
`--units` and the stratified-refraction fixes.

## The 2017 Noise Abatement databases (2026)

The AS350, B407, Be206L3, EC130, RO-44 and RO-66 databases were built by panam from the 2017
Noise Abatement flight test. They have been rebuilt several times, most recently on
2026-10-04 (the AS350 on 2026-10-05). The manual's
[Where the 2017 databases come from](MANUAL.md#where-the-2017-databases-come-from) records the
current recipe, what each rebuild fixed, and the approach-validation statistics before and
after. The AW139, B205 and S-76D `_fixed_load_factor` files were converted from those aircraft's
discrete-load-factor `_M*` databases. They agree with their sources to 1e-6 dB at load factors
from 0.5 to 3.0.
