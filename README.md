# NICE-OPS

**N**oise **I**nformed **C**ommunity **E**nvironment **O**perations **P**lanning **S**ystem — a
helicopter noise footprint tool. It predicts ground-level dBA noise exposure over an observer
grid (or a scattered observer set) from a flight trajectory and a rotor-noise "sphere" database.

NICE-OPS was developed by Eric Greenwood and Rob Rau. It began life in D; this is the C++17
implementation, which reproduces the D original's results on the demo case to within a few
units in the last place.

## Building

Requirements: a C++17 compiler and netCDF. Plotting additionally needs Python 3 with
development headers; matplotlib itself is installed by the build, not by you. On Ubuntu:

```sh
sudo apt-get install libnetcdf-dev python3-dev
```

Plotting can be dropped entirely, which removes the Python requirement with it:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNICEOPS_PLOTTING=OFF
```

CLI11, nlohmann_json, Catch2 and pybind11 are fetched at configure time, so the first configure
needs network access.

Build out of source, always with explicit `-S`/`-B`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

This configuration ships: `assert()` is compiled out, as CMake's Release normally has it. To test
the way CI does, keep the asserts live with `-DNICEOPS_ASSERTS=ON`.

On macOS with Command Line Tools but no full Xcode, add
`-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)`.

Link-time optimization is on wherever the toolchain supports it, since the acoustic loop and
the database it calls live in different translation units; `-DNICEOPS_LTO=OFF` turns it off if
a build needs to be faster to produce than to run.

Run the tests with `ctest --test-dir build -C Release`.

## Running

With no arguments, NICE-OPS runs a demo case: the trajectory `tracks/S-76D/197421AC.csv`
against the database `databases/S-76D_M3.nod` and vehicle `databases/S-76D.json`, writing
`integrated/197421_LAE.csv`.

```sh
./build/niceops                                     # the demo case
./build/niceops -t my_track.csv -d my_db.nod -v my_vehicle.json -o out.csv
./build/niceops --help                              # the full option list
```

Paths are relative to the working directory, so run from the repository root when using the
shipped data.

### Inputs

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV. Column 0 is UTC seconds, column 1 elapsed time; positions and velocities are in **feet**. |
| Trajectory (`-t`), WGS84 | A CSV whose header names `time`, `lat`, `lon` and `alt` columns — see [Geodetic coordinates](#geodetic-coordinates). Any column order, and optional `vE`/`vN`/`vUp`. |
| Trajectory (`-t`), static | A CSV whose header names `duration`, `x`, `y`, `z` and `heading`: each row is a hold at one point for `duration` seconds, nose on the compass `heading`. Every hold shares one position, in `--units`; a move is a trajectory. Sampled at 50 Hz, at 1 g with zero velocity. |
| Database (`-d`) | netCDF `.nod`: one group per operating condition (advance ratio μ, flight path angle γ, thrust coefficient C<sub>T</sub>), each holding scattered (φ, θ) channels of dBA level and EAA excess-attenuation factor. |
| Vehicle (`-v`) | JSON with `Main Rotor.radius`, `Main Rotor.tip speed`, `Vehicle.weight` (mass in kg) and `Atmosphere.density`. |
| Observers (`-b`) | CSV of scattered observer positions, in feet: `x,y,z`, `id,x,y,z`, or `id,x,y,z,height` with each microphone's height above the ground (see [Receivers](#receivers)). The first line is always taken as a header. Without it, a 120×120 grid spanning −1500…1475 m is used. |

Input and output files use **feet**; the computation is metric internally.

### Outputs

- `-o` — sound exposure level (SEL) per observer. Gridded output writes the x and y axes on the
  first two lines and then one row of levels per y.
- `-D` — day-night level, summed over every trajectory in the run (see below).
- `-i` — per-observer level time histories.
- `--band_history` — per-observer band time histories, for comparing with measured spectra
  frame by frame. One file per observer, `PREFIX_i<N>.csv`, with a row per trajectory state:
  `emission_time`, `reception_time` (s), `range` (the ray tube's under `--refraction`),
  `elevation` (the emission direction below horizontal, deg), the sphere's `phi` and `theta`,
  `level_dBA` (the A-weighted, band-limited level `-i` writes), then `band_<Hz>`: each band's
  unweighted level at the microphone, dB re 20 µPa, with spreading, absorption, the ground or
  plate and any shadow loss applied. A-weighted and summed in energy the bands give
  `level_dBA`. Bands outside `--band_limits` and culled states are left empty; `-a`'s ambient
  goes into `level_dBA` only. Needs the spectral path.
- `-m` — appends the energy-mean SEL, maximum and above-threshold counts for the case.
- `-S` — the derived trajectory state (t, position, velocity, load factor, C<sub>T</sub>, μ, flight path angles).
- `--powerfile` — per-state shaft power, if the database carries any.
- `--latlon` — writes a `_latlon` companion beside every level file above, as latitude, longitude, level.
- `--export_sphere` — a whole interpolated noise sphere, written as a one-group `.nod`.

## Units

`--units` selects the length unit for the whole interface: `ft` (the default) or `m`. It covers

- every length on the command line — `--grid`, `--gridstep`, `--trim`, `--ground_elevation`,
  `--receiver_height`, `--path_smoothing`;
- every file carrying a length — trajectory input, observer lists, level output, time histories;
- the axis labels on the figures.

The computation is metric throughout. The conversion happens once, where the command line is
parsed, and nothing downstream of that seam sees anything else.

**`--units` says what your files already are; it does not convert them.** The shipped tracks and
observer sets are in feet, so running them with `--units m` reads feet as meters and flies the
aircraft three times too far. Nothing can detect that, so the flag is a declaration you have to
get right. Where a file names its own unit in a column header — `alt_ft`, `alt_m` — the name is
checked against `--units` and a disagreement fails the file rather than being silently believed;
a bare `alt` claims nothing and is always accepted.

Two things stay put whatever you select:

- **Defaults are physical quantities, not numbers in the current unit.** The observer grid
  defaults to ±1500 m in steps of 25 m whichever unit is in force — it reads as ±4921.26 ft or
  ±1500 m, the same ground either way. If defaults scaled with the flag, `--units m` would
  quietly shrink the footprint to a third of its area.
- **The `-S` state file is metric.** It is a dump of the internal state, and
  `golden/state_197421.csv` — the D-derived oracle — encodes it that way.

Non-length quantities keep the units of the standards they come from, because converting them
would make them harder to check against the source: temperature in K, pressure and flow
resistance in kPa (ISO 9613-1 is written in kPa), frequency in Hz, angles in degrees, time in
seconds, levels in dB.

> Before this flag existed, `--grid`, `--gridstep`, `--receiver_height` and `--path_smoothing`
> were meters while `--trim` and `--ground_elevation` were feet and every file was feet. Scripts
> written against that need those four converting, or `--units m` adding.

## Command line arguments
## Command line arguments

### Inputs and outputs

| Option | Meaning |
| --- | --- |
| `-d`, `--database` | netCDF noise database. Once, or once per trajectory. |
| `-t`, `--trajectoryfile` | CSV trajectory. Repeat for several flight paths. |
| `-v`, `--vehicle` | JSON vehicle definition. Once, or once per trajectory. |
| `-o`, `--output` | SEL output file. Once per trajectory (default `integrated/197421_LAE.csv`); with several tracks and one name or none, each output is named after its track. Two trajectories asked to write one file are refused. |
| `-b`, `--observer` | Scattered observer file; switches off the default grid. |
| `-S`, `--state` | Save the derived trajectory state. With several tracks and one name, each is saved as `<track>_state.csv` beside it. |
| `-i`, `--history` | Prefix for per-observer time history files, `<prefix>_i<N>.csv`. With several tracks each file is also named for its track, `<prefix>_<track>_i<N>.csv`, so they do not overwrite one another. |
| `--band_history` | Prefix for per-observer band time history files (see Outputs), named as `-i` names its files. Needs the spectral path. |
| `-m`, `--metricfile` | Append overall noise metrics to this file: the case name, the energy mean of the per-observer SELs (10 lg of the mean of 10^(L/10), not a mean of decibels), their maximum, and the share above 65 dB. |
| `--origin` | WGS84 reference for the local frame: latitude longitude, degrees. Defaults to the first position of the first geodetic track. |
| `--frame_bearing` | Compass bearing of the local +x axis. Default: x east, y north. The shipped AC tracks use 210. |
| `--ground_elevation` | Elevation of the observer plane above mean sea level, in `--units`. Subtracted from a WGS84 track's altitude column. |
| `--path_smoothing` | How far, RMS per axis in `--units`, the path fitted through a WGS84 track may sit from the given positions. Default 10 ft. Zero interpolates exactly, which on measured positions makes the load factor track coordinate rounding. |
| `--start_utc` | Seconds since UTC midnight at *t* = 0 of a WGS84 track, for the day-night window. Required with `-D` when the track has no `utcsec` column (it used to default to midnight, which is night, so the whole footprint took the 10 dB penalty). |
| `--latlon` | Also write each level file as latitude, longitude, level. |
| `-l`, `--logfile` | Log file (default `log/niceops.log`). |
| `-L`, `--loglevel` | `trace`, `info`, `warn`, `error` or `fatal` (default `info`). |

### Physics and metrics

| Option | Meaning |
| --- | --- |
| `-a`, `--ambient` | Ambient level in dB, added to every sample's level in energy (two equal levels give +3 dB). Off unless given. It used to be added as pressure amplitude, as the D original did, which read up to 3 dB high. |
| `-q`, `--qsam` / `--no-qsam` | Use the flight path angle relative to the tip-path-plane normal rather than the vertical. On by default. |
| `-s`, `--step` | Resampling step for trajectory adaptation. **Not** the SEL integration step, which comes from the trajectory's own sampling. |
| `--units` | Length unit for the command line and for every file carrying a length: `ft` (default) or `m`. See **Units**. |
| `--trim` | Crop each trajectory from this x, in `--units`, to its lowest/slowest point. Off by default; it changes which part of the flight is integrated. |
| `--grid`, `--gridstep` | Observer grid bounds and spacing, in `--units`. Default ±1500 m (±4921.26 ft), step 25 m (82.021 ft) — a fixed extent either way. |
| `-D`, `--dnlfile` | Write the day-night level, summed over every trajectory: the energy average over a 24-hour day, with the 22:00–07:00 night penalized. A daytime operation reads its SEL less 49.37 dB. |
| `--lamax_file` | Write each observer's LAmax: the highest A-weighted level over the operation, dB. Forces the whole history, so the peak is exact. |
| `--nax_file`, `--nax_threshold` | Write each observer's NAx: the number of operations whose own LAmax reaches the threshold, dBA (required with `--nax_file`: there is no default, so the study's threshold is always stated). The count adds across trajectories. EPNL is not implemented. |
| `--coverage_file` | Write, per observer, the energy-weighted share of its levels that came from cells the database did not measure. Needs a database that records `coverage` (none of the shipped ones does yet). Straight rays only: refused with `--refraction` or `-a`. |
| `--full_history` | Integrate each SEL over the whole history rather than the event's 10 dB-down window. Slower, since no window is found. |
| `--utc_offset` | Hours from UTC to local time, for the DNL night window. |
| `--powerfile` | Write per-state shaft power, when the database has it. |
| `--spectral`, `--no-spectral` | On by default when the database carries spectra: propagate its per-band spectra to the ground with ISO 9613-1 absorption. `--no-spectral` takes one broadband level with the database's EAA factor instead (databases without spectra always do). |
| `--band_limits` | Frequency range integrated into the metric, Hz: `low high`. Default `20 20000`. |
| `--humidity`, `--temperature`, `--pressure` | Atmosphere for `--spectral`, overriding the vehicle file. Percent, K, kPa. The spectral path refuses to run without a humidity (no `--humidity`, no vehicle-file value, no `--atmosphere`): absorption depends strongly on it, so there is no silent 20 %. |
| `--atmosphere` | A stratified atmosphere: a CSV of temperature, humidity, and optionally pressure and wind, against height. Absorption is integrated through it. Spectral path only. See below. |
| `--refraction` | `none` (default); `linear`, curved rays through the `--atmosphere` profile's effective sound speed fitted to a straight line per path; or `stratified`, through the profile's own layers. Both give shadow zones under upward refraction. See below. |
| `--shadow_floor` | Most a refractive shadow zone attenuates, dB. Default 20. |
| `--tunneling` | `--refraction stratified`: also take the rays turned back from aloft that turn again just above the microphone, which the wave reaches through the layer below them (an Airy tail per band). Off by default; see "Micro-shadows and skip zones" below. |
| `--skip_fill` | `--refraction stratified`: fill the skip zone short of where rays from aloft first land with the caustic's Airy tail per band. Off by default; see below. |
| `--scattering` | Fill refraction shadows with the energy turbulence scatters into them, at HARMONOISE strength γT, instead of `--shadow_floor`. Default 0 (off). Needs `--refraction` and `--ground_effect`. |
| `--speed_reference` | Which velocity sets the flight condition (μ, the path angle and the sphere's nose): `ground`, or `air`, the ground velocity less the `--atmosphere` profile's wind at the aircraft. Default: what the database records it was built with, ground if it records nothing. See below. |
| `--ground_effect`, `--no-ground_effect` | Add excess ground attenuation (Chien–Soroka over a Delany–Bazley impedance) to the spectral path. On by default wherever the spectral path runs on a database that carries spectra; `--no-ground_effect` turns it off. Asking for it with `--no-spectral` is an error. |
| `--max_load_factor` | Largest load factor, in g, that a `fixed_load_factor` database's load-factor correction is taken to (default 2). Above it the correction is held, and the run says how many states it held. |
| `--accept_vehicle_mismatch` | A vehicle file that disagrees with the database by more than 2 percent in tip speed, rotor radius, weight or thrust coefficient at 1 g is refused. Between 1 and 2 percent it is a warning; this flag turns the refusal into a warning. |
| `--flow_resistance` | Ground flow resistance, kPa·s/m². 30 snow, 225 grass (default), 650 soil, 3000 dirt, 10000 concrete, 50000 asphalt. |
| `--receiver_height`, `--turbulence` | Point receiver's microphone height (in `--units`, default 1.2 m; an observer file's fifth column overrides it per observer) and turbulence (rad·s/√m, 0–16e-4) for `--ground_effect`; the height also places the microphone for `--propagation_table`. |
| `--impedance_model` | Ground impedance: `delany_bazley` (default), `miki`, `variable_porosity` (with `--porosity_rate`, 1/m) or `delany_bazley_layer` (with `--layer_depth`, in `--units`). `--flow_resistance` is the (effective) flow resistivity. |
| `--turbulence_model`, `--gamma_t` | `chessell` (default, `--turbulence`) or `harmonoise`, whose coherence takes `--gamma_t` (~1e-6 moderate, ~1e-5 strong). |
| `--roughness` | RMS roughness height, in `--units`: the coherent reflection falls by exp(−2(kσ cos θ)²). Default 0. |
| `--plate_table` | Observers are flush microphones on a rigid ground plate (GRAS 67AX in a GR1425 plate). Takes a netCDF table from panam; needs `--ground_effect` over the ground the table was built for. See below. |
| `--propagation_table` | A propagation loss table (netCDF, from `niceops_propagation_table`): each band's level relative to free field by bearing, source height, receiver height and range, for one stratified atmosphere over uniform ground. Replaces the ground and ray models wherever it covers a pair; pairs it does not cover take the run's own model, `--refraction` included. A plate microphone reads it at its lowest receiver height plus the plate's offset over bare ground. Spectral path. See below. |
| `--exact_ground` | Evaluate the ground model (or the plate's reflection coefficients) at every observer and state instead of interpolating per-trajectory tables of it: 3-4x slower, for checking the tables. |
| `--plate_mic_azimuth` | Which way the microphone is offset on its plate: `outboard` (default, away from the track) or degrees in the local frame from +x toward +y. |

### Performance

| Option | Meaning |
| --- | --- |
| `-j`, `--threads` | Worker threads for the observer loop. 0 (default) uses one per hardware thread; 1 runs serially. |
| `--db_cutoff` | Skip evaluations that provably fall below this level. See the caveat below. |

### Spheres and diagnostics

| Option | Meaning |
| --- | --- |
| `--export_sphere` | Interpolate a whole noise sphere and write it as a one-group `.nod` database. |
| `--export_aam` | The same sphere in AAM/RNM format. Needs a database with spectra. |
| `--sphere_condition` | The operating condition for it: advance ratio, flight path angle, C<sub>T</sub>. |
| `--self_test` | Validate every sphere in the database against its own stored levels. |
| `-p`, `--plotting` | Plot the footprint with the flight track over it; `--cmin`/`--cmax`/`--cstep` set the contour levels. |
| `-f`, `--figfile`, `-F`, `--statefigfile` | Save the footprint and trajectory-state figures to files. |

## Notes on particular features

### Several trajectories, and the day-night level

`-t` may be repeated. `-d`, `-v`, `-o` and `-S` then take either one value, shared by every
trajectory, or one per trajectory.

```sh
./build/niceops -t morning.csv -t evening.csv -t night.csv -D integrated/dnl.csv
```

DNL is a property of a set of operations rather than of one flight: each trajectory's SEL-style
integral gains a 10 dB penalty on samples flown at night, and the results are summed as energy.
Night is 22:00–07:00 local; `--utc_offset` shifts the trajectory's UTC seconds into local time.
The sum is normalized to a 24-hour day, as the DNL is defined: a single daytime operation
gives the SEL less 10 lg 86 400 = 49.37 dB. (Before the normalization, `-D` wrote the raw sum,
which read 49.37 dB high.)

Plotting is skipped for multi-trajectory runs, since one figure file cannot hold several paths.

A trajectory that cannot be read is logged and skipped rather than ending the run, so one bad
file in a batch costs that file alone. The exit status is nonzero if anything was skipped.

### `--db_cutoff`

On a large observer grid most of the work goes into observers that are never loud enough to
matter. `--db_cutoff L` skips any evaluation that provably falls below `L`, bounding the level
from the loudest sphere the state blends and spherical spreading alone.

The threshold sits 10 dB below `L` and is applied to the **instantaneous level**, which makes the
result exact for any observer whose *peak level* reaches `L`: the SEL integral covers the window
from the peak down to peak − 10, so no sample in that window is ever culled.

The cells of the output are SELs, though, and an SEL exceeds its observer's peak level by 1 to 12
dB (the window's duration), so a cell can be at or above `L` in SEL while its peak is below `L`,
and its window then loses samples. On the demo case, `--db_cutoff 85` changes 157 of the 592 cells
whose SEL is at or above 85 dB, by at most 0.43 dB; a cell is unchanged to the last bit once it is
about 3 dB above the cutoff. Read the cutoff as a threshold on peak level, and leave a few dB of
margin when it is the SEL you care about.

Observers that never reach the cutoff come out as `-inf` — they contributed no energy, and those
cells will not plot. That is the trade the flag asks for; it is off by default.

### Spectra

A database may carry per-band spectra alongside the overall A-weighted level: a
`frequency` dimension and an `amplitude(phi, theta, frequency)` array per sphere. When it
does, `--export_aam` can write an interpolated sphere back out in the AAM/RNM hemisphere
format the source spheres arrive in. Without spectra the flag refuses, since
`AMPLITUDE(PHI, THETA, FREQUENCY)` is the whole payload of that format.

Six of the databases in `databases/` carry spectra: the AS350, B407, Be206L3, EC130, RO-44 and
RO-66 `_fixed_load_factor` files. The S-76D demo database does not, so the demo runs the
broadband path. To build one with spectra, use panam's
`build_empirical_database(..., include_spectra=True)`, which now carries the per-band
amplitudes through from the source spheres. Band-limit it — `infreqs=(10.0, 3500.0)` —
because the top few bands of a measured sphere are usually noise-floor artifacts.

The AAM/RNM source spheres live in the group's OneDrive under `Data/Sphere_Database`, one
directory of `.nc` files per aircraft (AS350, Be206L3, Be407, Be430, EC130, MD-902, Mi-8, R-44,
R-66, SH-60B), each with its own `vehicle.cfg` and `vehicle.json`. `Be407_spectral.nod` sits
beside them — 26 bands, 10–3150 Hz, built from `Be407/` as above. The built `.nod` files in
`databases/` are tracked; the source spheres are far too large and stay in OneDrive.

#### Where the 2017 databases come from

The AS350, B407, Be206L3, EC130, RO-44 and RO-66 databases are built from the 2017 Noise
Abatement flight test by panam (`noise_abatement_2017.build_all`, then
`build_empirical_database`). The validation harness keeps the recipes (`spheres/rebuild.py`,
`spheres/hover_sphere.py`, `spheres/rebuild_database.py`). Each file records its build at the
root: `speed_reference`, the build atmosphere, and `DOPPLER_SHIFT_REMOVED` and `coverage` per
group.

- **B407** is the reference aircraft and the one validated most. Its spheres are depropagated
  along curved rays through each run's own atmosphere, with a 2-degree emission floor, an
  adaptive 2-degree grid and tone-aware third-octave bands. A rotor harmonic near a band edge is
  filed whole in its own band instead of being split by whole FFT bins.
- **The B407 hover** is synthesized from the slowest level pass with the Doppler shift taken
  out, then corrected by the measured hovers. Two rings were measured: 18–25 degrees below the
  horizon (out of ground effect) and 0.5–3 degrees (in ground effect). Each ring is a Fourier
  series in azimuth from the heading, with a sum-to-zero offset per heading for the wind. The
  two are blended linearly in elevation from 3 to 17.5 degrees. Below the high ring each
  azimuthal order fades toward the nadir as cos^m.
  - Nothing measures 3–17 degrees or the region under the aircraft, so treat the hover there
    as an interpolation.
  - The low ring is the noisier of the two (5 dB held-out rms against 2.8).
- **The other aircraft** keep the September recipe: straight rays, a 10-degree floor and a
  25-degree grid. AS350, EC130 and RO-44 still hold the September spheres. Be206L3 and RO-66
  are built without the approach cards, which are the validation set.

A database written before panam c53588f (2026-09-26) may have its last group unwritten. Every
entry is then the netCDF fill value (9.97e36), which the loader reads as a finite 1e37 dB
level. The September AS350, EC130, RO-44 and RO-66 files had this, and spectral runs that
reached the fastest +35-degree condition gave an infinite SEL.

#### A note on levels changing

**Results from before September 2026 are not comparable to results after.** `eval` used to
multiply the database's `EAA` by a constant 4.0, under a `//TODO Add ground effects` inherited
from the D original. That is not a ground effect, and it overstated atmospheric absorption
fourfold — 168 dB at 6 km where the physical value is 14. panam *writes* that `EAA` as
`SPLA - SPLAa` over a 1000 m reference and consumes it as `EAA / distance * (slant_range -
radius)`; the tool now does the same.

Removing it raised the demo footprint by **9.93 dB mean, up to 22.7 dB**, most of it far from
the track where the old factor was removing the most. `golden/197421_LAE.csv` was regenerated;
`golden/state_197421.csv` is trajectory state and unaffected.

It was measured against `--ground_effect`, the model it was standing in for, before being
removed: soft ground costs a roughly range-independent 2–4 dB relative to a rigid surface,
where the factor removed 7 dB at 250 m rising to 38 dB at 2 km. A real ground effect is bounded
by the surface impedance; a multiplier on a distance-proportional absorption is not.

#### Propagating them to the ground

The spectral path, which is the default whenever the database carries spectra (`--no-spectral`
turns it off), carries the spectrum itself to the observer instead of a single number. Each
band is spread by the same `20 log10(r/r_sphere)` and absorbed at its own ISO 9613-1
coefficient over `r - r_sphere`, and only then are the bands A-weighted and summed into a
level. The two paths answer different questions and will not agree:

| | `--no-spectral` | spectral (the default with spectra) |
| --- | --- | --- |
| Attenuation | database's `EAA`, one dB/km factor | ISO 9613-1 α(f), per band |
| Fitted for | the spectrum the source happened to have | the atmosphere you state |
| Atmosphere | none — `EAA` is baked in | temperature, pressure, humidity |
| Needs | any database | one carrying spectra |

The absorption model is ported from panam's `panam_acoustics/iso_9613_1_1993.py` and pinned
against it in `tests/test_atmosphere.cpp`; the A-weighting is panam's `dBAw`, pinned in
`tests/test_weighting.cpp`.

**Band-limit deliberately.** Summing whatever bands a database happens to carry quietly
makes the answer a property of that database's bandwidth, so two databases are not
comparable. `--band_limits low high` states the range; bands outside it contribute exactly
nothing, and the run logs how many of the database's bands survived. A limit that excludes
every band is an error rather than silence.

The atmosphere comes from the vehicle JSON's `Atmosphere` section:

```json
"Atmosphere": {
    "density": "1.225",
    "temperature": "300.71",
    "pressure": "105.74",        // optional; else derived as rho*R*T
    "relative humidity": "20"    // optional; a spectral run is refused without one
}
```

Only `density` is required, as before. Pressure defaults to the one implied by that file's
own density and temperature rather than to a standard day — `S-76D.json` pairs 1.225 kg/m³
with 300.71 K, which is not an ISA point, so assuming 101.325 kPa would contradict its own
numbers. `--humidity`, `--temperature` and `--pressure` override all of it.

#### Stratified atmosphere

`--atmosphere file.csv` replaces that single atmosphere with a profile, and integrates each
band's absorption through it. Along a straight path, what a band loses is the ISO 9613-1
coefficient's mean over the heights the path spans, times the path length. That is the
fine-layered form of what certification does with 30 m layers (14 CFR 36 Appendix A), and a
unit test checks the two agree.

On the demo track with the B407 database, `atmospheres/isa.csv` moves 90% of the footprint by
-0.10 to -0.005 dB (mean -0.04 dB, extremes -0.57 and +1.11 dB), against a homogeneous run at
ISA sea-level conditions. On a grid it costs nothing measurable, because every observer stands
at one height and each state's mean is formed once. Scattered observers at different heights
pay a subtraction and a multiply per band per evaluation.

Columns are found by name, in any order, one row per level; rows may run top down:

| Quantity | Column names | Unit |
| --- | --- | --- |
| Height above the observers | `z`, `agl` | `--units`; `z_ft` / `z_m` must agree with it |
| or altitude above sea level | `alt`, `altitude`, `height` | `--units`, less `--ground_elevation` |
| Temperature | `T_K`, `T_C` | named; a bare `T` is told from the values |
| Relative humidity | `RH`, `humidity` | percent |
| Pressure (optional) | `p_kPa`, `p_hPa`, `p_Pa` | named; a bare `p` is told from the values |
| Wind (optional) | `wind_speed_mps` or `wind_speed_kt`, with `wind_direction` | direction is the compass bearing it blows **from**, degrees |
| or wind components | `u_mps` and `v_mps`, or `u_kt` and `v_kt` | toward east and toward north |

```csv
z_ft,T_K,RH,p_kPa
0,288.150,20,101.325
1000,286.169,20,97.717
2000,284.188,20,94.213
```

Without a pressure column, pressure is integrated hydrostatically up from the surface
pressure, which is the vehicle file's or `--pressure`'s. A single row is a uniform
atmosphere. `atmospheres/isa.csv` is the ISA to 10,000 ft at 20% RH, and
`atmospheres/B407_uniform.csv` holds `B407.json`'s own atmosphere. With the latter a run
reproduces the homogeneous one, which a CTest case checks.

The rules that keep it from guessing:

- **The profile has to reach the aircraft.** A trajectory climbing above the top level is
  refused, not extrapolated. Below the lowest level, that level is held down to the ground.
- **Surface properties stay at the surface.** The profile's level at the ground supplies the
  sound speed for the reception-time integral and for the ground impedance.
- **`--temperature` and `--humidity` cannot be combined with it.** Nor can `--pressure`
  unless the file gives no pressure, in which case it is the surface pressure.
- **Only the spectral path uses it.** The broadband path's `EAA` is the database's own
  fitted attenuation, and there is nothing per band to integrate.
- **Wind columns must name their unit.** Nothing in the values tells m/s from knots. Wind
  matters only to `--refraction` and to `--speed_reference air`; absorption does not depend
  on it. The log names every column it ignored.

#### Speed reference

A database's flight conditions come from some velocity: μ = |v| / V_tip, the path angle
asin(v_z / |v|), and the sphere's azimuth reference, which is v projected into the tip-path
plane. Every database so far took v over the ground. `--speed_reference air` takes it
through the air instead: the ground velocity less the profile's wind at each state's height,
turned into the local frame by `--frame_bearing`. The load factor, and so CT, stays inertial,
because the rotor's thrust balances the aircraft's acceleration over the ground whatever
the air does.

- **The database says which it is,** in a root text attribute `speed_reference`, `ground`
  or `air`, and NICE-OPS follows it by default. A file without one is ground, which is
  every database so far; panam does not write the attribute yet.
  Giving `--speed_reference` against the database is logged as an override.
- **Air needs a wind.** With no wind columns in `--atmosphere`, an air-referenced run is
  refused rather than taken as calm. `--speed_reference ground` says so explicitly.
- **It is after `--trim`,** which cuts an approach at its slowest point over the ground.
- **Calm air changes nothing**, to the bit, which a CTest case checks.

What it moves, on the 2017 Noise Abatement approaches: 1–3% in μ and about 1° of nose on
calm days. On a 15 kt day, 5–10% in μ and 8–15° of nose. Matching a ground-built database with
an air-referenced trajectory mixes the two, so change both together.

#### Refraction

`--refraction linear` bends the rays through the profile. For each path it fits a straight
line to the effective sound speed, c(z) + w(z)·h, where w is the wind and h the horizontal
direction of propagation. The fit is least squares over the heights the path spans, widened
to at least 10 m. This is the equivalent linear profile of Harmonoise and Nord2000. In a
linear profile every ray is an arc of a circle, so the eigenrays have closed forms, and the
path supplies:

- **The sphere's direction:** the direct ray's launch angle. At kilometer ranges it is a few
  degrees from the straight line, which is a large fraction of the database's 5° grid.
- **Spreading:** from the ray tube's area, which in a linear profile is exactly
  r / cos(launch angle).
- **Absorption:** along the arc.
- **Ground geometry:** a point receiver's ground model gets the reflected ray's grazing angle
  and length; a plate microphone gets the direct ray's arrival angle. Both ground models are
  otherwise unchanged.
- **Reception time:** from the ray's travel time.
- **Shadow zones:** under upward refraction, receivers beyond the ray that grazes the ground
  get the grazing ray's level, decaying with distance into the shadow. The decay is the least
  attenuated creeping-wave mode. With `--ground_effect` it runs over the ground's own
  impedance; without, over rigid ground, about 24 dB/km at 1 kHz for a gradient of 0.02 /s.
  Real ground is not always the deeper shadow: to first order it attenuates up to 12% less
  than rigid at low frequency and weak gradients, and up to 2.3 times as much once the ground
  counts as soft on the wave's scale. The benchmark's 1 kHz ground, at 0.1 /s, is 146 dB/km
  against rigid's 70. The rate is held to the slope of the exact residue series deep in the
  shadow to 0.3%. It is capped at `--shadow_floor`, where energy scattered in by turbulence
  takes over; `--scattering` puts that energy in directly instead (below). A ground reactive enough to carry a surface wave into the shadow (a porosity
  rate of tens per meter, a thin hard-backed layer, low bands) takes the rigid rate in those
  bands, and the log names them.
- **Near the shadow boundary**, with a point receiver or a flush plate microphone, the level
  comes from the exact field instead. Within a couple of diffraction lengths of the boundary
  (F = 2k₀ℓ², about 100–400 m) ray theory fails on both sides: its ground term, taken at
  the grazing ray's geometry, sinks toward the grazing ground wave, 10–25 dB below the exact
  field over grass, and over rigid ground it doubles where the exact field does not.
  NICE-OPS tabulates the residue series of the path's linear profile per band, relative to
  the ray's own spreading, as a function of distance past the boundary in diffraction
  lengths. That function barely depends on the source's height, which is Fock's result, so
  one table per band and microphone height serves every path. The ray model stands until
  2F before the boundary, blends into the table by 1F before it, and the table (floored 20 dB
  below its boundary level) takes over from there. It follows the exact field to 0.1–1 dB
  for sources from 30 to 300 m. With `--refraction stratified` each band gets its own
  boundary first, as the next section describes.

`atmospheres/example_windy.csv` is an illustrative profile, not a measurement: a neutral
log-law westerly (5.2 m/s at 10 m, 8.5 m/s at 1000 ft) over the ISA lapse. On the demo track
with the B407 database and `--frame_bearing 210`, compared with the same profile under
`--refraction none`:

| | mean | upwind half | downwind half | worst cell | acoustic loop |
| --- | --- | --- | --- | --- | --- |
| no ground effect | -0.77 dB | -1.51 dB | -0.04 dB | -7.85 dB | 2.3 s (1.6 s without) |
| `--ground_effect` | -0.56 dB | -0.80 dB | -0.32 dB | -3.81 dB | 3.6 s (2.1 s without) |
| `--plate_table` | | | | | 6.1 s (5.4 s without) |

Most of the upwind loss comes from the low states near the end of the approach. From 10 m up,
the shadow begins about a kilometer upwind.

What to know before relying on it:

- **Wind directions are compass bearings**, turned into the local frame by `--frame_bearing`,
  whose default is +x east. The shipped AC tracks run +x along 210°, so without
  `--frame_bearing 210` their wind is read 120° off. The log says so whenever wind is present
  and the bearing was not given.
- **One gradient per path.** A log-law profile's steep lowest meters are averaged into it.
  Equivalent-linear models have been found up to 10 dB off a parabolic-equation solution near
  the ground at 500–2000 Hz (van Maercke & Salomons 2009). The general stratified model
  planned as step 6 removes this approximation.
- **The effective-sound-speed approximation** overstates the wind by 1/cos(elevation). That is
  worst on steep paths, where refraction matters least.
- **One direct and one reflected ray.** The extra ground-bouncing arrivals of strong downward
  refraction at long range are not included.
- **The source does not feel the wind.** Its advance ratio still comes from the track's ground
  speed.
- **The ray ends at the microphone,** as the straight path does (see [Receivers](#receivers)),
  so uniform air gives the straight result to 1e-9 for every receiver, scattered microphones
  on sloping ground included.
- **The peak search and `--db_cutoff` still work.** They bound a refracted state with the
  same direction tables as a straight one, taken along the ray's launch direction at its ray
  tube's range. A CTest case checks that the windowed run equals a full march.

#### Ground effect

`--ground_effect` adds the interference between the direct ray and the one reflected off the
ground — Chien–Soroka over a Delany–Bazley impedance, with Chessell's boundary-loss correction
and band averaging, ported from panam's `ega()` and pinned against it in
`tests/test_ground_effect.cpp`.

**It is a gain, not a loss**, despite the name: over rigid ground with the receiver on the
surface it tends to +6 dB, the pressure doubling of a perfect reflector. It goes sharply
negative only in the interference dip, where the reflected path arrives out of phase. On the
Be407 demo case:

| | mean | max |
| --- | --- | --- |
| spectral, no ground | 65.10 | 96.38 |
| + concrete (10000) | 68.45 (**+3.35**) | 99.00 |
| + grass (225) | 64.73 (**−0.36**) | 97.89 |

Hard ground is 3.7 dB louder than soft, and grass roughly cancels its own reflection gain
against the dip. `--flow_resistance` spans five orders of magnitude and the answer is sensitive
to it, so name the surface rather than guessing.

`--turbulence` decorrelates the two paths and fills the nulls in; leaving it at zero gives
unrealistically deep interference dips.

The impedance, turbulence and roughness models follow Attenborough & Taherzadeh, *Prediction
of outdoor ground effect*, Applied Acoustics 242 (2026) 111067, as panam's `ground_plane`
module does. They are pinned against `ground_plane.pole_level` in
`tests/test_ground_effect.cpp`. Delany–Bazley is poor for grassland. On the 2017 Amedee pole
microphones, variable porosity with σe ≈ 200 fitted best, with HARMONOISE γT ≈ 3e-5 and
10–15 mm roughness (panam `docs/ground_plane_corrections.md`).

#### Refraction through the profile's layers

`--refraction stratified` drops the single straight-line fit. The effective sound speed is
taken as linear *between the profile's own levels*, and held constant below the lowest and
above the highest. Inside each layer a ray is still a circular arc, so crossing a layer has a
closed form, and a ray is a sum over its layers.

Eigenrays are solved per source–observer pair, over four families of rays:

- launched down, arriving from above;
- launched down, dipping below the microphone and turning back up to it;
- launched up, turning aloft, then arriving from above;
- launched up, turning aloft, then dipping.

The primary arrival is the first family's root, found from the linear model's answer. A
receiver at or above the source is left to the linear model.

The range need not be monotonic along a family, so the pair can have more arrivals than that.
A source at a local maximum of the sound speed is one way. Its shallowest rays skim the slow
layer and travel far, its steepest turn high and travel far, and the arrivals lie between. The
other families are sampled in launch angle, 32 times plus 8 for each bounce the path can hold
(at most 512), and every crossing of the range is solved for.
Each arrival is a level of its own, through its own launch direction on the sphere and its own
ray tube, and they are summed in energy.

**Surface ducts.** When the air above the source outruns all of the air below it, a ray that
reaches the ground turns aloft and comes back down. A morning inversion does this. A ray that
bounces m times before its last descent travels X_0(p) + m L(p), where L is one hop from the
ground to its turning height and back. Those bounced rays are found the same way, from the same
samples. Each earlier bounce costs the plane-wave reflection loss, and the last is the ground
model's, as for the primary. Beyond the reach of every direct ray, about 2.4 km for a 30 m
source under day 286's morning inversion, the bounced rays are all there is. The strongest
becomes the primary. `atmospheres/example_inversion.csv` is that profile, from run 286276,
without its 1–3 kt of wind: 3.5 K of inversion capped near 350 ft.

Without them, `stratified` was 8–24 dB under the PE past 4,000 ft on day 286's profiles for
100–330 ft sources, and up to 54 dB under it at single ranges. The table compares the energy
mean over 50 Hz – 2 kHz against the PE, per range bin, for a flush receiver:

| day 286, source 100 ft | 1,000–4,000 ft | 4,000–8,000 ft | 8,000–13,000 ft |
| --- | --- | --- | --- |
| run 286264 | +0.4 | −1.1 (was −12.2) | −0.5 (was −24.3) |
| run 286270 | 0.0 | +0.8 (was −7.2) | −2.5 (was −8.1) |
| run 286276 | −0.7 | −1.5 (was −10.7) | −1.1 (was −9.2) |
| run 286282 | +4.4 | +0.2 | −3.1 (was −17.1) |

Downwind on `example_windy.csv` and a log law, at 100 ft, the far third of a 10,000 ft path goes
from 4.0–6.5 dB under the PE to within 1 dB of it. Upwind, the linear lapse and the day-284
profiles do not move.

Two limits remain. Within one pair of samples, two close arrivals at a caustic can still be
missed. Against a 2048-sample reference (12 bearings, sources 5–300 m, 50 m to 10 km), 2.9% of
pairs under `example_inversion.csv` are still more than 1 dB out, and 0.01% under
`example_windy.csv`; a fixed 32 samples left 11.3% and 1.05%. Scaling the count with the hops
cost 4–8% of a stratified run.
Day 286's last profile (run 286288) is an elevated duct. Toward 270 the air is fastest at the
ground, because the log-law headwind assumed below the LIDAR's first gate (46 ft) takes
0.2 m/s off in the lowest foot, and a ray turned back from the inversion near 330 ft turns up
again 1–20 ft above the ground. For a 100 ft source it reaches neither a flush nor a 1.2 m
microphone at 8,000–13,000 ft, which stay on the shadow's floor 11–21 dB under the PE, while
the wave tunnels through the layer (flattening it moves the PE by under 1 dB to 250 Hz and
about 6 dB at 1 kHz). Short of the duct's return lies a skip zone that no ray reaches, which
the PE fills from the caustic 1–2 kft early. Measured day-286 frames agree with the PE: where
`stratified` predicts 10 dB or more under straight rays, it is 10–20 dB too quiet. Those
frames are far from each microphone's peak and do not move the SEL.

It is checked two ways. On a linear profile it reproduces the linear model's closed forms,
the shadow boundary included. On layered profiles its rays match numerically traced ones to a
micrometer. Everything else follows `--refraction linear`: the ray tube, the ground geometry,
reception time, the shadow's decay and floor, and the peak search.

What the layers change is near the ground. A log-law wind is zero at the ground, so for
propagation upwind the fastest air is at the surface. Snell's law then fixes the limiting
ray's launch angle from the speeds at the source and at the ground alone, whatever lies
between. For a source at 100 m in `atmospheres/example_windy.csv`, the layered model's
upwind shadow begins 535 m out, where a fitted line puts it past 1.6 km.

**The wave does not take that boundary at face value.** Rays follow the thin near-ground
shear layer, but the wave resolves it only once its diffraction length is short against the
distance between the two boundaries. At 125 Hz the parabolic equation below shows no shadow
until about a kilometer, and even at 2 kHz it has not reached the layered boundary. So with a
point receiver each band gets an effective boundary between the two,
r_b = r_s + w (r_l − r_s), where w = x / (x + 0.363) and x is the band's diffraction length
at the linear fit's curvature over r_l − r_s. The boundary field is then scaled by the
curvature that makes r_b. Past the layered boundary the lit level is the linear fit's ray:
- The layered model's near-grazing rays there have no reflected partner.
- On a log law its dipping family reaches tens of meters past its own end.

The law's one constant was fitted on twelve log laws (u* 0.2–0.6 m/s, z0 0.01–0.1 m,
sources 30 and 100 m, 125 Hz – 2 kHz) and checked on `example_windy.csv`, which was held out:
6.8 dB rms over the shadow edge against the PE there, against 18.9 dB with the layered
boundary as it stands. On the demo grid with a point receiver, `stratified` then averaged
3.1 dB below `linear` (5th percentile −12.4 dB), where it was 6.1 (−19); with the fixes below,
3.3 (−13.2). Runs without
`--ground_effect` keep the layered boundary and the creeping-wave loss.

**Flush plate microphones** get the same treatment:
- **The field:** the boundary field for a receiver at the plate's top (8 mm for the GR1425).
- **The plate over bare ground, K:** added on top. Near grazing K barely depends on elevation
  (within 0.1 dB from 0.5° to 10°). It is about 0 dB up to 1 kHz, +0.9 dB at 2 kHz, +2.6 at
  4 kHz and +7 at 8 kHz, and it varies with azimuth by ±1 to ±2.5 dB, so it is tabulated per
  band and azimuth at 2° and applied at each pair's own azimuth.
- **At ground level** the boundary law holds up. Against the PE for a receiver 3 mm up, over
  16 log-law cases, the shadow-edge error is 7.1 dB rms, against 10.0 for `linear`. The
  constant that suits both heights together is 0.39, only 0.1 dB better than 0.363.

Two faults in the layered model surfaced with flush microphones, both fixed:
- **Microphones just below the plane.** Ground boards sit a few centimeters below the survey
  datum. The layers hold their lowest speed below z = 0 and the direct-ray solve knew no
  ground, so it sent a ray grazing along under the plane. On the 2017 test's windy day that
  was 33 dB low. Such a microphone is now solved for on the plane.
- **Receivers a few millimeters up.** The ray tube past the boundary came out 0/0. It is now
  spread as a linear profile's limiting ray.

**The measured profiles of the 2017 test's windy day** (LIDAR winds and balloon temperatures,
morning of day 284) exposed three more, all in how `stratified` meets its shadow. They were
found by running the PE on the geometries of the microphones the model got most wrong.
- **An inversion aloft hid the boundary.** Those mornings have an upwind shear near the
  ground under an inversion, so the air a hundred meters up is a little faster than at the
  source. The layered model took any faster air above the source to mean rays from aloft
  reach the receiver, and set no boundary. With none, there was no transition and no
  hand-over to the linear fit: the shadow came in as a 10–35 dB cliff at the layered
  boundary, which the PE puts 30–50% further out. A ray turned down aloft turns up again
  wherever the air is faster still, here in the shear, meters above a board. Rays from aloft
  reach the receiver only if the air aloft outruns everything between the two, which is the
  test now.
- **The layered ray focuses short of its boundary.** The steep lowest layer that sets the
  boundary also bends the rays that just clear it into a caustic, a spike of up to +10 dB the
  PE does not show. From 0.9 of the layered boundary on, the level is the linear fit's ray.
  Against the PE over the last half of the way to the boundary, that is 0.6 dB rms on the
  day-284 profiles (1.4 dB handing over at the boundary itself) and 1.7 on log laws (3.3).
- **Past the fit's own boundary the level stepped.** The model went back to the layered
  shadow ray, whose spreading differs, by 2 dB on a log law and 10 on day 284. The fit's
  limiting ray now carries on, with its spreading continuous.

**The day-286 approaches' hovering end** (sources 15–50 ft, under 20 kt) showed two more.
`stratified` was 6.5 dB over straight rays there on the measured frames, 63–250 Hz, at
microphones 150–450 ft away; the PE on those paths agreed with straight rays.
- **Faster air aloft still hid the boundary near the source.** The test above is right about
  whether rays from aloft reach the receiver at all, but they come down only past their skip
  distance, kilometers out. Short of it the boundary now stands, and is taken back where a ray
  from aloft or along a duct does arrive. Without it the shadow took the layered shadow ray,
  focused 6.5 dB by the log law's lowest meters, plus the whole ground term: 12 dB over free
  field at 63 Hz where the PE gives 4.
- **A low source leaves the transition no room.** When a band's boundary is nearer the
  source than two diffraction lengths, the table has no lit zone to hand over from, and the
  source is under the 5 l the table was checked to. Such a band keeps the ray model to its
  boundary and blends to the table over the next diffraction length.

Together, on the measured frames:

| 63–630 Hz energy-mean error | before | after | straight rays |
| --- | --- | --- | --- |
| sources under 100 ft and 20 kt | +8.9 dB | +4.4 dB | +3.8 dB |
| within 2,000 ft, 2–5° | +4.7 dB | +3.1 dB | +2.4 dB |
| every other height, speed, range or elevation cell | | within 0.2 dB of before | |

The 4 dB left at the hover is in straight rays too, so it is the source. On the day-286 PE
matrix past 1,000 ft, `stratified` against the PE went from 9.5 to 7.5 dB rms at 200–630 Hz
and 11.7 to 9.2 at 800–2,000 Hz (5.4 to 5.1 at 50–160).

**The boundary law, refitted for flush microphones.** With those fixed, the law put flush
microphones' boundaries about twice too far out. Fitted per case against the PE, the best
boundary is about 0.9 of a diffraction length at the fit's curvature past the layered one on
day 284 and 1.2 on flush log laws, where the law gave 1.7–2.7. A microphone 1.2 m up keeps
0.363 (refitted with the day-284 profiles added: 0.387). A flush one takes 0.75, between the
flush log laws' own optimum of 0.47 and day 284's 0.94. Where the linear fit casts no shadow
at all, as over a strong inversion, the law has nothing to interpolate toward. There each
band takes 1.60 (f / 1 kHz)^−0.08 of the layered boundary for a flush microphone, and 1.90
(f / 1 kHz)^−0.08 at 1.2 m, the per-band optima pooled over the calibration set.

Shadow-edge rms against the PE (from its last level above −3 dB to −25 dB), dB:

| | before | now |
| --- | --- | --- |
| day-284 profiles, 18 cases (sources 50, 120, 250 m), flush | 12.9 | 8.1 |
| log laws, 16 cases, 1.2 m | 6.2 | 6.0 |
| log laws, 16 cases, flush | 7.3 | 7.6 |

"Before" for day 284 already has the first fix; without it the cliff made those shadows
worse still. The best any boundary per case can do is 2.7 dB pooled, and no smooth law found
gets below about 6: the optimum is sharp where the shadow wall is steep, at high frequency.

Against the measured day-284 SEL, held-out spheres, `stratified` on the plate (729 mic-runs,
SNR ≥ 10 dB, records that miss the flyover left out), median / rms in dB:

| | before | now |
| --- | --- | --- |
| upwind, peak elevation 5–10° (65) | +0.45 / 3.72 | +0.28 / 3.06 |
| upwind, 10–20° (80) | +0.99 / 2.24 | +0.71 / 1.87 |
| all | +0.09 / 1.85 | +0.09 / 1.70 |

The calm day 286 barely moves (overall 1.46 → 1.45 dB rms): its winds cast few shadows.

Costs, with a point receiver: the demo grid's acoustic loop took 12.0 s when `stratified` went
in (3.3 s linear, 2.1 s without refraction). A 96-microphone scattered array like
`arrays/coyle_approach.csv` takes 35 ms. Searching for every arrival costs about half as much
again. A whole demo-grid run (B407, `--ground_effect`, `--frame_bearing 210`, 14 threads) goes
from 23.6 to 37 s on `example_windy.csv` and from 57.5 to 86 s on `example_inversion.csv`. The
search is skipped past the hand-over to the linear fit, unless a duct gives bounced rays to
look for. Focusing at a caustic is capped at +10 dB.

#### Checking the refraction models against a parabolic equation

`niceops_pe_compare` runs one path through all four models and writes each band's level
relative to free field: straight rays, `linear`, `stratified`, and a wide-angle parabolic
equation. The last is a full-wave solution, with the ground's impedance written into the
grid. Absorption and turbulence are left out of all four, unless `--gamma_t` puts turbulence
in the PE (below).

The PE is marched in range by a split-step Padé method (Collins, JASA 93, 1736 (1993)): one
range step applies the one-way equation's exact propagator, replaced by a rational function
that is accurate to 1e-7 for every wave within 50° of horizontal (`--max_angle`, up to 60) and
costs one tridiagonal solve per term. A step is eight wavelengths long (`--range_step`) and
takes about ten terms. The result does not depend on the step: between 2 and 32 wavelengths it moves
by less than 0.001 dB, 60 dB into a shadow. Waves steeper than `--max_angle` are absorbed
rather than carried, so **a level on a path steeper than that from the source is written as
NaN** (within 83 m of the point under a 100 m source, at the default).

Two more things make the field, and not just the propagator, good to that angle. The march
starts from Collins' self-starter (`--starter self`): the first range step with the point
source's modal amplitudes in it, applied to a delta function on the grid, instead of a
narrow-angle Gaussian. And the grid is fourth order in height throughout: the derivative
(`--vertical_order 4`, Numerov's form), the ground condition, and the interpolation that puts
the source on the grid and takes the receivers off it (`--boundary_order 4`). Against the
exact two-ray field over soft ground at 1 kHz, from a 100 m source:

| Elevation | Gaussian starter, second order, 20 points | self-starter, fourth order, 10 points | … 20 points |
| --- | --- | --- | --- |
| 0–20° | 0.01–0.02 dB rms | 0.005 | 0.005 |
| 20–30° | 0.14 | 0.006 | 0.003 |
| 30–35° | 0.59 | 0.01 | 0.002 |
| 35–40° | 0.61 | 0.02 | 0.004 |
| 40–50° | 1.75 | 0.01–0.07 | 0.001–0.007 |

Fourth order is also what lets the grid be coarse. It is 10 points per wavelength by default
(`--ppw`; 20 under `--scheme cn` or with either order set to 2). On
`atmospheres/example_windy.csv` out to 1.5 km, 250 Hz – 2 kHz, against an 80-point run:

| | 20 points | 10 points | 6 points |
| --- | --- | --- | --- |
| second order throughout | 0.02–0.04 dB rms (0.14 worst) | 0.06–0.19 (0.64) | |
| derivative fourth, ground and interpolation second | 0.01–0.02 (0.06) | 0.05–0.09 (0.29) | |
| fourth order throughout | under 0.0005 (0.001) | 0.001–0.002 (0.011) | 0.004–0.01 (0.08) |

The middle row is why the ground condition and the interpolation had to follow the
derivative: with those at second order the whole scheme converged at second order.

The self-starter also corrects a bias at every angle. The Gaussian is normalized with the
wavenumber at the ground, so a source aloft in a refracting or windy column comes out at the
wrong strength: 0.11–0.15 dB high against the exact residue series for a 100 m source where
the sound speed is 3% lower than at the ground, and 0.10 dB on `example_windy.csv`. The
self-starter is within 0.04 dB.

`--scheme cn` is the march this replaced, Crank–Nicolson on the Padé(1,1) wide-angle equation
at ten steps per wavelength (`--range_ppw`), kept as an independent check. It is good to about
30°, and its phase error builds with range. Through `atmospheres/example_windy.csv` from a
100 m source out to 1.5 km, 125 Hz – 2 kHz, the two differ by 0.03–0.12 dB rms with the same
starter and grid order; against the exact field in still air it is Crank–Nicolson that is off,
by 0.7 dB rms at 25–30° and 2.1 dB at 30–35° at 1 kHz.

| Five bands, 125 Hz – 2 kHz, 100 m source to 1.5 km, 14 cores | PE time |
| --- | --- |
| `--scheme cn` | 25.4 s |
| `--scheme pade` (default: 50°, self-starter, fourth order, 10 points) | 0.47 s |
| … `--max_angle 60` | 0.67 s |
| … `--max_angle 35` | 0.37 s |
| … `--ppw 20` | 1.4 s |
| … `--ppw 20 --max_angle 35 --vertical_order 2 --boundary_order 2 --starter gaussian` | 0.6 s |

(The ray models' columns take another 0.4 s.) One march is 39 times faster than
Crank–Nicolson at the defaults. Bands and their frequencies run in parallel. For
Crank–Nicolson, against a
30 × 30 Crank–Nicolson reference, 10 in height is within about 0.1 dB rms for a receiver near the ground, but
0.47 dB rms (2.2 dB at worst) for one 1.2 m up at 1 kHz; 20 brings that to 0.09 dB, at twice the
cost. On the stored shadow-calibration cases the change moved the PE by 0.16–0.39 dB rms (1.1 dB
at worst) and left the calibration's edge error against it unchanged.

The PE is checked against exact solutions:

| Case | Agreement | Test |
| --- | --- | --- |
| Two rays over rigid ground | within 0.3 dB | `tests/test_pe.cpp` |
| Weyl–van der Pol field over grass | within 0.5 dB | `tests/test_pe.cpp` |
| Residue series in an upward-refraction shadow, rigid ground | within 0.02 dB at 250, 500 and 1000 Hz, down to 47 dB below free field | `tests/test_pe_shadow.cpp` |
| Benchmark case 1: direct and image rays over the benchmark's ground | within 0.2 dB at 10 Hz from 50 m (0.29 dB at 25 m, inside a wavelength of the source) and 0.08 dB at 100 Hz and 1 kHz from 20 m, to 1 km | `tests/test_benchmark.cpp` |
| Benchmark case 2: normal-mode series, downward refraction | within 0.06 dB at 100 Hz to 1 km and 1 kHz to 500 m, 0.2 dB at 10 Hz | `tests/test_benchmark.cpp` |
| Benchmark case 3: residue series, upward refraction, impedance ground | within 0.03 dB at 100 Hz and 1 kHz, to 60 dB below free field | `tests/test_benchmark.cpp` |

The shadow rows are the checks that matter for arbitrating the ray models. Their profiles
have k² linear in height, and the field in the shadow is then exactly a sum of creeping waves
built from Airy functions of complex argument. Doubling the grid resolution moves no shadow
onset and no level past 100 m by more than 0.7 dB.

```sh
niceops_pe_compare --atmosphere atmospheres/example_windy.csv --toward 270 \
    --source_height 328 --receiver_height 3.937 --bands 125 250 500 1000 --max_range 8000
```

Through `atmospheres/example_windy.csv` with a 1.2 m receiver, this is where each model falls
10 dB below free field and stays there, upwind:

| source | band | straight | linear | stratified | PE |
| --- | --- | --- | --- | --- | --- |
| 100 m | 125 Hz | never | 1.9 km | 1.2 km | 1.3 km |
| | 250 Hz | never | 1.4 km | 0.53 km | 1.07 km |
| | 500 Hz | never | 1.3 km | 0.53 km | 0.91 km |
| | 1 kHz | never | 1.4 km | 0.53 km | 0.84 km |
| 30 m | 125 Hz | never | 0.91 km | 0.72 km | 0.65 km |
| | 250 Hz | never | 0.46 km | 0.38 km | 0.50 km |
| | 500 Hz | 1.2 km | 0.42 km | 0.19 km | 0.38 km |
| | 1 kHz | never | 0.50 km | 0.23 km | 0.34 km |

And the rms difference from the PE beyond 400 m, downwind of the 100 m source:

| band | straight | linear | stratified |
| --- | --- | --- | --- |
| 125 Hz | 2.8 dB | 4.0 dB | 0.9 dB |
| 250 Hz | 7.1 dB | 7.6 dB | 0.4 dB |
| 500 Hz | 3.1 dB | 0.9 dB | 0.9 dB |
| 1 kHz | 8.1 dB | 7.2 dB | 2.1 dB |

What that says:

- **Downwind and in the lit region, `stratified` is right** and the others are not.
- **Upwind, the shadow begins between the two ray models' predictions.** It moves toward
  `stratified` as the frequency rises and the wavelength shrinks against the near-ground shear
  layer, which is what diffraction would do. `stratified` starts the shadow too early and
  `linear` too late, by similar amounts.

`tools/plot_pe_compare.py OUTDIR CSV:TITLE ...` plots that output: one figure per case, a
panel per band, and a summary of where each shadow begins. It runs under the build's Python,
`build/pyenv/bin/python`.

**The shadow's edge.** Over the stretch where the PE falls from −3 to −25 dB, the ray models'
level error against it (rms, dB):

| case | before | with the boundary field |
| --- | --- | --- |
| linear lapse (−0.018 /s), sources 30 and 100 m, 125 Hz – 2 kHz | 5–20 | 0.1–1.0 |
| log-law winds, `stratified`, 16 cases, 125 Hz – 2 kHz, pooled | 30.7 | 6.1 |

| log-law winds, pooled rms | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz |
| --- | --- | --- | --- | --- | --- |
| `stratified` before | 7.1 | 22.6 | 34.6 | 37.4 | 39.3 |
| with the boundary field and the per-band boundary | 6.0 | 5.2 | 5.7 | 4.7 | 8.4 |

On the linear profile, where the ray model's geometry is right, the boundary field puts it on
the PE. On log laws it needs the per-band boundary too. Scaled by the layered model's own
boundary, it falls 10–20 dB too fast at 125–250 Hz and pooled at 19.5 dB.

For a receiver on the ground (3 mm up) the same comparison gives 7.1 dB pooled, against 10.0
for `linear`, with the lapse control at 0.5 dB. By band, 125 Hz is the worst (10.2 dB) and
500 Hz the best (4.4).

Two things that looked like the treatment and are not:
- **Clipping the profile below a fixed height.** Once the height passes the microphone, the
  microphone sits in the clipped layer and shadows at every range.
- **Smoothing the profile per band.** A linear profile smoothed only over its lowest meters
  moved the lapse control by 20 dB at 2.5 km. The ray model is too sensitive to structure
  below the microphone for either to work.

The PE's deep shadow in still air, 40–80 dB down, is not to be believed either: with no
turbulence in it, nothing scatters energy in.

#### Turbulence, and the energy it scatters into shadows

The PE takes turbulence as random phase screens (Gilbert, Raspet and Di 1990): `--gamma_t`
(HARMONOISE strength), `--outer_scale` (von Kármán, default 10 m) and `--seed`, with each of
`--samples` frequencies its own realization, so a band level averages over the band and the
ensemble at once. The effective index's structure parameter is γT / 4, the value for which
Ostashev's coherence and HARMONOISE's both reduce to one form (Attenborough & Taherzadeh 2026,
eqs. 18 and 20). One screen per outer scale of range carries the phase that slab accumulates.
Its vertical spectrum is the 3-D one integrated across the plane, so the vertical coherence
is the 3-D medium's. Checks, `tests/test_pe_turbulence.cpp`:
- **The screens' variance:** 0.782 k² Δx (γT/4) L0^(5/3), to 0.06%, less the 4% the finite
  column leaves out.
- **The ensemble-mean field:** it falls as exp(−σ²/2) with range to within the sampling noise
  of 300 realizations, which checks the screens and the march together.
- **Single scattering:** a third of the strength lowers the deep-shadow level by 4.82 dB,
  against 10 lg 3 = 4.77.

At γT = 3e-5 (what the 2017 co-located pole and board pairs fit, in panam), the PE fills the
day-284 shadows to −3 dB re free field at a board, median, only 5–10 dB below the lit level.
`--scattering` adds that energy to the shadow-boundary field, in place of the floor:
S = a + 10 lg γT + b lg(f / 1 kHz) + c lg(R / 100 m) + d lg(z_s / 100 m) dB re free field.

| receiver | a | b | c | d | fitted on |
| --- | --- | --- | --- | --- | --- |
| flush (3 mm, or a plate) | 39.18 | −1.55 | 2.10 | 5.30 | 12 profiles, sources 50 and 120 m, 1.5 dB rms |
| 1.2 m | 36.67 | +5.61 | 5.13 | 5.30 (the flush fit's) | 6 profiles, sources 50 m |

S is held to 0 dB, R to 100 m–10 km and z_s to 10–300 m; a plate lifts it by K as it lifts
the coherent field.

**The measurements do not want that much.** Against the day-284 SEL, the upwind 5–10° bin
goes +0.28 / 3.06 dB (median / rms) with the floor, +0.49 / 3.04 at γT = 1e-6, +0.95 / 3.02 at
3e-6, +1.55 / 3.16 at 1e-5 and +2.96 / 4.00 at 3e-5. Up to about 3e-6
is consistent with them, which is HARMONOISE's "moderate" range. The pair fit's 3e-5 is ten
times that. It is three times HARMONOISE's maximum, and probably carries more than turbulence
(roughness, the pole's height), or turbulence at the pole's height, where it is strongest.
`--scattering` stays off by default, since the floor does as well on these data.

#### The Attenborough et al. benchmark

Attenborough et al., "Benchmark cases for outdoor sound propagation models", JASA 97, 173
(1995), is the standard set of test cases for outdoor propagation codes. It fixes a source at
5 m, a receiver at 1 m, 10, 100 and 1000 Hz, a locally reacting ground with impedances it
tabulates, absorption at 20 °C and 70% RH, and four profiles:

| case | profile | exact solution |
| --- | --- | --- |
| 1 | uniform | direct and image rays with the spherical-wave reflection coefficient (their eq. 12) |
| 2 | c rising 0.1 /s (downward refraction) | normal-mode series (eq. 10) |
| 3 | c falling 0.1 /s (upward refraction) | residue series of creeping waves (eq. 11) |
| 4 | rising 0.1 /s to 100 m, falling to 300 m, uniform above (a duct) | none |

`niceops_benchmark` runs all four to 200 m and to 10 km with the PE and, where one exists,
the exact solution (`src/exact_solutions.hpp`). `tools/plot_benchmark.py DIR` draws them in
the paper's layout, so each figure sits beside their figs. 11–14, and prints how far the PE is
from the exact solution. A full run takes about 10 seconds.

```sh
niceops_benchmark --out figures/benchmark
build/pyenv/bin/python tools/plot_benchmark.py figures/benchmark
```

Set beside the paper's FFP and PE curves, all four cases agree in every panel. Where there is
an exact solution, the PE is within 0.3 dB of it past 20 m on every 200 m panel. On the 10 km
panels it is within 0.3 dB for case 1, and for case 3 within 0.1 dB, or 0.5 dB at 10 Hz,
where the linear profile and the exact solution's profile part company within a few Airy
lengths of the ground. Case 2 at 10 km is the exception, covered below. Case 4 has no exact
solution; it reproduces the paper's curve feature for feature, down to the 100 Hz humps at
2.5–3.2 km and 4.3–4.7 km.

Three things came out of it:

- **Their eq. (11), the upward-refraction residue series, has the ground condition wrong.**
  It writes the eigenvalue condition as Ai'(b) + i k0 l beta w Ai(b) = 0, w = exp(2 pi i/3),
  where the derivation gives Ai'(b) - i k0 l beta w² Ai(b) = 0: the admittance term is off by
  a factor exp(i pi/3). Over the benchmark's 100 Hz ground the printed form is 1.5 dB off the
  PE at 100 m and 12 dB off at 1 km. The derived form agrees with the PE to 0.02 dB down to
  -90 dB. Over rigid ground the two coincide, which is why that case never showed it.
  `src/exact_solutions.cpp` carries the derivation.
- **Case 2's normal-mode solution carries a field no PE can.** Its profile, c0/sqrt(1 - 2z/R),
  reaches infinite sound speed at z = R/2 = 1716 m, so every ray turns below that height,
  including the near-vertical ones. Those fall back within R = 3.4 km of the source, landing
  at R sin(2 × elevation), like mortar rounds. They put structure at 3.4 and 6.8 km into the
  mode sum; the paper's own normal-mode figure has it, and its FFP and PE figure does not.
  With modes steeper than 75° left out (`max_elevation`), the series and the PE agree to
  0.06 dB out to 1 km at 100 Hz. With them in, they differ by 5 dB.
- **Case 2 at 10 km depends on where the profile is cut off.** The paper says so, and it holds
  here: the PE with the linear profile stopped at 1000 m matches their curves in form, not
  interference null for null.

The inputs are theirs. The impedances are the four-parameter model's values they give, not a
NICE-OPS impedance model. The absorption is the 0.18 dB/km at 100 Hz and 5.4 dB/km at 1 kHz
they quote from ANSI S1.26-1978. ISO 9613-1 gives 0.22 and 4.98 dB/km at the same conditions.

#### Receivers

Each observer's level is evaluated at its microphone: the observer's point raised by the
microphone's height above the ground. What the microphone is decides the ground model:

| Microphone | Ground model | How to ask for it |
| --- | --- | --- |
| flush with the ground, height 0 | the plate (below) | `--plate_table` |
| on a pole, height h | the point receiver at h | `--ground_effect`, with h from `--receiver_height` or the observer file |
| anything, without `--ground_effect` | none: the free field at the microphone | |

With a four-column observer file, every observer takes the run's receiver: all on plates
with `--plate_table`, otherwise all at `--receiver_height`. A fifth column gives each
microphone's own height, in `--units`, so one run can hold an array of flush microphones with
a few on poles, as the test campaigns do. With `--plate_table`, height 0 is a plate and
anything else is a point receiver at that height:

```csv
id,x_ft,y_ft,z_ft,height_ft
2,-788.3,-1085.1,-1.1,0
10,-405.9,1940.8,-10.1,3.937
```

The log reports how many observers took each receiver. `tests/data/mixed_receivers.csv` is a
worked example, and a CTest case checks that it gives the same levels as separate runs over
each kind of microphone.

Before the microphone height set the evaluation point, the point receiver's direct path ended
on the ground even though its ground model put the microphone 1.2 m up. Correcting that moved
point-receiver SEL on the demo grid by up to about ±0.8 dB (median 0). Almost all of it comes
from the low states at the end of the approach, where 1.2 m is a sizeable fraction of the
source's height. Plate and free-field results did not change.

#### Propagation loss tables

`--propagation_table` takes what the path does to each band from a table computed elsewhere,
so that a full-wave method -- the parabolic equation above, or something else later -- can
carry a footprint. The table holds each band's level relative to a point source's free field
at the straight-line distance: the ground, refraction, shadows and diffraction, everything
beyond spherical spreading. The run then takes, per band,

    L = L_sphere(straight-line direction) - spreading(straight-line range) - absorption + table

which is the composition `--ground_effect` uses, with the table where the ground model and
the ray model were. Atmospheric absorption stays the run's (homogeneous, or layered through
`--atmosphere`) unless the table says it includes it.

It assumes one stratified atmosphere over flat, uniform ground, for now. The level then
depends on where the observer is only through the compass bearing and horizontal range from
the source and the two heights above the ground, which are the table's axes:

| Axis | |
| --- | --- |
| `bearing` | compass bearing the sound travels toward, degrees, uniform from 0 over 360; one value means the same in every direction |
| `source_height`, `receiver_height` | above the ground, m, ascending |
| `range` | horizontal, m, ascending |
| `band` | band centers, Hz; the run's bands within `--band_limits` must all be there |

Values are interpolated linearly in dB along each axis. A NaN means the method has no value
there. A pair the table does not cover, by a NaN among the values it would interpolate or by
falling outside an axis, takes the run's own propagation instead: `--refraction` if it is
on, `--ground_effect` if that is, free field if not. A flush microphone on a plate
(`--plate_table`) reads the table at its lowest receiver height and adds K, the plate's level
over bare ground there, per band and azimuth at 2 degrees, as the shadow transition places a
plate. With `--refraction` too, the run evaluates every state rather than windowing from the
peak, since the bounds were not built for the two geometries together. The log says how many did. Receivers are the exception: a microphone
height outside the table's is refused, since that table was built for another run. Bearings
are compass bearings, so a run in a rotated frame needs `--frame_bearing`, as with wind.

The file is netCDF, SI whatever `--units` says: a `relative_level(bearing, source_height,
receiver_height, range, band)` array, one coordinate variable per axis, and global attributes
`propagation_table_version` (1), `method`, `includes_absorption` and `provenance_*`.
Anything that can fill that array can be read; `src/propagation_table.hpp` is the
specification.

### Micro-shadows and skip zones

Under day 286's last profile (run 286288) the air toward 270 is fastest at the ground, because
the log-law headwind assumed below the lidar's 46 ft gate takes 0.2 m/s off in the lowest
foot, and the inversion near 330 ft makes an elevated duct above it. A 100 ft source's sound
turned back by the inversion turns up again 1-20 ft above the ground, so ray theory delivers it
to neither a flush nor a 1.2 m microphone; and short of where it first comes down, about
12,000 ft out, no ray arrives at all. The parabolic equation hardly notices the near-ground
layer (flattening it moves the PE by under 1 dB to 250 Hz) and fills the skip zone from the
caustic. Three treatments were built and compared (harness `propagation/micro_shadow*.py`):

- `--tunneling` adds the rays that turn within reach of the microphone, each with the Airy
  tail below its turning point, Ai^2((3/2 omega I)^(2/3)), I the integral of
  sqrt(p^2 - 1/c^2) from the microphone to the turning height. Near the family's caustic the
  pair of rays is replaced by the caustic's band-averaged Airy field.
- `--skip_fill` adds, for each family from aloft whose range has an interior extremum with the
  receiver on its dark side, the caustic's Airy tail, I = C Ai^2(omega^(2/3) (X_c - r)
  (2 / X'')^(1/3)), C matched to the two rays' energy on the lit side.
- `--propagation_table` built from the PE for each run's profile, beside `--refraction
  stratified` for the pairs it does not cover.

| 100 ft source, flush, 2,000-13,000 ft; `stratified` - PE, mean / rms dB | 50-160 Hz | 200-630 Hz | 800-2,000 Hz |
| --- | --- | --- | --- |
| Day-286 PE matrix (five profiles), neither | -5.4 / 9.8 | -6.1 / 11.1 | -3.6 / 9.9 |
| `--tunneling` | -3.4 / 7.3 | -3.8 / 8.5 | -1.4 / 7.3 |
| `--skip_fill` | -4.5 / 8.9 | -5.6 / 10.5 | -2.7 / 9.2 |
| both | -1.7 / 5.1 | -2.3 / 6.5 | +0.3 / 6.7 |

Against measured day-286 frames where `stratified` predicts at least 10 dB under straight rays
(1.4% of frame-bands to 2 kHz; sources under 300 ft, 0.4-2 degrees, 4-9 kft), energy-mean
predicted - measured:

| | straight | `stratified` | `--tunneling` | `--skip_fill` | both | PE table (5 runs) |
| --- | --- | --- | --- | --- | --- | --- |
| 25 runs, 50-160 Hz | +3.8 | -10.5 | -7.3 | -3.6 | -2.1 | |
| 25 runs, 200-630 Hz | -0.1 | -15.5 | -10.7 | -7.2 | -5.0 | |
| 25 runs, 800-2,000 Hz | -7.6 | -20.6 | -20.6 | -20.7 | -18.9 | |
| 5 runs, sources under 150 ft at 3-6 kft, 50-630 Hz | +1.0 to +1.6 | -11.5 to -13.3 | | | -2.5 to -4.5 | +0.3 to +2.5 |

None moves the SEL: the exposure within 10 dB of each microphone's measured peak is unchanged
to 0.02 dB under every option. All three are off by default. What they leave:

- The tails above 800 Hz. The tunneled field at the ground and the skip zone's edge short of
  the direct rays' reach (a duct's direct family ends without a caustic, and the PE diffracts
  past it) are still 10-20 dB low there.
- On run 286276's far frames (150-330 ft sources, 6-20 kft) every model but plain `stratified`
  is too loud, the PE table by 13.5 dB and straight rays by 8.8: something the PE leaves out.
- The table is read across 45 degree bearings and eight source heights; off those nodes it is
  within about 1 dB of a direct PE run at 63-250 Hz but 5 dB rms overall in this duct, almost
  all at 1-2 kHz, where the return's caustic moves between nodes. Each run's table took
  2.5 minutes to build.
- Cost on a day-286 run: `--tunneling` about +20%, `--skip_fill` 2-4x, a table about the
  ray model's own.

`niceops_propagation_table` writes one, by either of two methods:

```sh
# The parabolic equation through a profile, 36 bearings where it has wind
niceops_propagation_table --method pe --atmosphere atmospheres/example_windy.csv \
    --database databases/B407_fixed_load_factor.nod --band_limits 20 2000 \
    --source_heights 40 90 140 190 240 --receiver_heights 3.937 \
    --max_range 15000 --step 10 -o tables/windy.nc

# The straight-ray ground model --ground_effect uses: the same footprint, through a table
niceops_propagation_table --method ground_model --database databases/B407_fixed_load_factor.nod \
    --sound_speed 335.636 --flow_resistance 225 --source_heights 40 50 60 ... 1130 \
    --receiver_heights 3.937 --max_range 15000 --step 10 -o tables/ground_model.nc
```

The second is the check on the plumbing. On the demo case its footprint matches
`--ground_effect --exact_ground` to -0.002 dB mean and 0.19 dB at worst with 10 ft spacing,
and to 0.02 dB at worst with 5 ft spacing. The ground effect itself moves cells by up to
4.3 dB. CTest runs that comparison (`propagation_table_is_the_ground_model`), and checks that
the peak search stays exact with a table that adds gain (`propagation_table_windowed_is_exact`):
the search's bound takes each band's largest table value.

With `--method pe` through `atmospheres/isa.csv`, 21 bands from 20 Hz to 2 kHz, 23 source
heights every 50 ft and ranges to 15,000 ft every 10 ft, the table takes 103 s on 14 cores. The
footprint's acoustic loop goes from 0.06 to 0.08 s. The table covered 99.7% of evaluations; the
rest were steeper than the PE's 50 degrees and took the ground model. Against the ray model's
footprint it reads -0.32 dB mean (-0.9 to 0 dB, 5th to 95th percentile; -2.8 dB at worst), the
standard lapse rate's mild upward refraction thinning the field at long range.

Interpolating between source heights in dB is where the table's error is. Over the ground model it
shrinks fourfold as the spacing halves, but the PE's band-averaged field is only as smooth as
its band. Build tables a few wavelengths apart at the highest band, or densely at the lowest
heights, where the interference changes fastest.

#### Ground-plane microphones

A flush ground-plane microphone sits on a plate: rigid close in, soft ground beyond. The point
receiver can model one surface or the other, not both. `--plate_table` models the plate
itself, using panam's axisymmetric boundary-element solution for a GRAS GR1425 plate (400 mm,
8 mm thick, lying on the ground) with the microphone at 3/4 R per SAE ARP 4055. The solution is
stored as two complex transfer functions of the direct and image waves:

    p / p_free = P_d + Q P_r

P_d and P_r are tabulated over frequency (five sub-frequencies per band), elevation and
azimuth. Q is the ground's spherical-wave reflection coefficient, which NICE-OPS computes for
each geometry. The table depends on the ground, so it is built for one ground and NICE-OPS
refuses any other. Two tables are shipped in `plates/`:

| table | ground |
| --- | --- |
| `gr1425_db225.nc` | `--flow_resistance 225` (Delany–Bazley grass, the default) |
| `gr1425_vp200.nc` | `--impedance_model variable_porosity --flow_resistance 200` (the fitted Amedee ground) |

Build another with panam's `axisymmetric_bem.table(bands, …)` and `write_netcdf`. A table at
the databases' 31 bands takes about 20 s on 14 cores. Turbulence and roughness do not apply to
the plate, because its reflection is millimeters from its direct path.

Q is taken as quadratic in log frequency across each band, through the band's edges and
center. This is within 0.07 dB of evaluating Q at every sub-frequency wherever the level is
above −15 dB. The exact per-sub-frequency form (`PlateReceiver::level_exact`) matches panam's
`board_level` to 1e-6 dB in `tests/test_plate_receiver.cpp`. Its cost on the full demo grid
(S-76D track, B407 fixed-load-factor database, 94.3 M evaluations, 14 threads):

| | acoustic loop |
| --- | --- |
| spectral, no ground | 1.6 s |
| `--ground_effect` | 2.1 s (7.1 s with `--exact_ground`) |
| `--ground_effect --plate_table` | 5.4 s (13.3 s with `--exact_ground`) |

By default both receivers interpolate per-trajectory tables instead of evaluating the model at
every observer and state. The tables are built exactly at their nodes in well under 0.1 s.
- **Point receiver:** the whole ground correction, over source elevation and log range. Level
  histories agree with the exact model within 0.001 dB.
- **Plate receiver:** its reflection coefficients. SEL agrees within 0.0005 dB rms.

The interpolated values never exceed the model's bounds, so the peak search and `--db_cutoff`
stay exact. `--exact_ground` evaluates the model directly. Either way the Faddeeva function uses
fast paths: an asymptotic series for |z| ≥ 6 and a Taylor table below, which agree with
Weideman's series within 1e-11.

The peak search and `--db_cutoff` bound the plate with `PlateReceiver::gain_bounds`, so both
keep working. The `plate_windowed_is_exact` ctest checks that skipped states change nothing.

**It is expensive.** Acoustic loop on the Be407 case (94.3 M evaluations, 26 bands, 14 threads,
default band limits): broadband 0.34 s, `--spectral` 1.44 s, `--spectral --ground_effect`
10.9 s. So carrying the spectrum costs about 4x, and the ground effect about 32x — a complex
reflection coefficient and a Faddeeva evaluation per band per observer per state, none of
which can be hoisted because it all depends on the geometry. The Faddeeva evaluation is about
half of the ground path on its own.

Those three are from before the September 2026 optimization round and all of them are now
lower. Measured against the same build on the S-76D demo track (94.3 M evaluations, 14
threads), over a synthetic single-condition 24-band database built the way
`write_flat_spectral_database` in `tests/test_spectral_propagation.cpp` builds one:

| | before | after |
| --- | --- | --- |
| broadband | 0.32 s | 0.16 s |
| `--spectral` | 1.36 s | 1.16 s |
| `--spectral --ground_effect` | 8.95 s | 7.81 s |

A single condition makes the operating-condition blend degenerate, so those two spectral
figures are optimistic against a real multi-condition database like the Be407 one in the row
above; the ground-effect term, which is most of the last row, does not depend on the condition
count. All three outputs agree with the previous build to 1.4e-12 dB.

`--db_cutoff` works on the spectral path through per-direction bounds on the band-limited,
weighted level (`Database::spectral_direction_bounds`). With `--ground_effect` or
`--plate_table`, each band is raised by the most the ground or the plate can add.

### How much gets evaluated

The sound exposure integral covers the window from each observer's first to its last sample
within 10 dB of its peak, dips included (the duration of the aircraft noise metrics), which on a
shipped case is a little under 60% of the track. The observer loop therefore does not march every
emission time: it finds the peak with a per-direction bound on the level, and evaluates only the
states that bound cannot rule out of the window. The run logs how much it evaluated:

```
Evaluated N of M (P%); the rest is provably outside every observer's 10 dB window.
```

It is exact, not an approximation — the output is byte-identical to evaluating every state. It
is used on the spectral path (`--spectral`, `--ground_effect`, plates, refraction, tables), where
an evaluation costs 150 ns to a microsecond; a broadband run, where an evaluation costs about 20
ns and the bound proves only about 6% of the pairs outside the window, marches every state
instead. It is skipped, with every state evaluated, for `-a`, `-i`, `--band_history`, `-D` and
`--db_cutoff` as well; see CLAUDE.md for why each one cannot use it.

### Threads

The observer loop runs on every core by default. Each observer is independent, so the output is
deterministic: runs at `-j 1`, `-j 2` and `-j 32` are byte-for-byte identical.

Observers are dealt out to the workers in chunks as they finish rather than split evenly up
front, because on a machine with both performance and efficiency cores — every Apple silicon
part, and most recent x86 laptops — an even split makes the slowest core the critical path.
On an M4 Pro the even split ran the demo case no faster on all 14 threads than on its 10
performance cores; pulling chunks recovers about 10%.

### Geodetic coordinates

A trajectory can be given in WGS84 coordinates — latitude, longitude and altitude — instead of
the 27-column AC state format. That is what a mission planner has: a route, not a flight-test
record.

```
time,lat,lon,alt
0.0,39.82914323,-74.40996141,1134.60
0.5,39.82903000,-74.41007000,1133.80
...
```

Both formats are read by `-t`, and which one a file is is decided by its width: 27 columns is
the positional AC format, anything else is read by the names in its header. Column order does
not matter and unrecognized headers are an error rather than a guess.

| Column | Meaning |
| --- | --- |
| `time` | Seconds from the start of the track. `utcsec` (seconds since UTC midnight) serves instead, and sets the clock for `-D`. |
| `lat`, `lon` | Degrees. |
| `alt` | **Feet above mean sea level.** `--ground_elevation` gives the elevation of the ground the observers stand on, which is subtracted to get height above them. |
| `vE`, `vN`, `vUp` | Optional, ft/s. Used when all three are present; otherwise velocity is differentiated from the path. |

Names are matched case-insensitively, ignoring punctuation and a parenthesized unit, so
`Latitude`, `LAT` and `lat (deg)` are one column. `alt_ft` is recognized and `alt_m` is
deliberately *not*: a file naming meters fails to load rather than being read as feet.

**Where the frame sits.** `--origin <lat> <lon>` fixes the local x/y origin; without it, the
first position of the first geodetic track is used. `--frame_bearing` turns the axes so +x runs
along a compass bearing — the default leaves them east/north, which is what a map wants. One
frame serves the whole run, because every trajectory in it shares the observers.

**What is derived, and how well.** A WGS84 track gives position only, so the velocity and the
load factor behind C<sub>T</sub> come from the first and second derivatives of a smoothing
spline through the path. `--path_smoothing` says how far, in feet RMS **per axis**, that fit
may sit from the given positions. It is an acceptable deviation, not an abstract weight: the
solver finds the smoothest curve whose RMS residual meets the bound, so asking for 10 ft means
the fit will use 10 ft of slack, not merely stay under it.

The default is 10 ft. Zero interpolates exactly, and exact interpolation is the wrong thing for
both kinds of input: measured positions make the second derivative track coordinate rounding,
and sparse waypoints are a statement of intent rather than a constraint to be met to the inch —
flying precisely through them manufactures aggressive maneuvering, which is itself a noise
source.

Feeding the demo track's own lat/lon/alt columns back in, against the footprint computed from
its full inertial state:

| Input | Median | 99th pct | Max | Cells over 1 dB (of 14400) |
| --- | --- | --- | --- | --- |
| 50 Hz positions, default `--path_smoothing 10` | 0.073 dB | 0.72 | 1.42 | 25 |
| 50 Hz positions, `--path_smoothing 0.0656` (= 0.02 m) | 0.031 dB | 0.36 | 1.44 | 6 |
| 50 Hz positions, `--path_smoothing 0` (exact) | 8.0 dB | 9.57 | 10.58 | 14400 |
| 50 Hz positions and velocity columns, default | 0.015 dB | 0.36 | 1.42 | 4 |
| Waypoints 10 s apart, resampled to 0.1 s, default | 0.121 dB | 0.86 | 2.10 | 76 |

The projection itself is not the lossy part: it reproduces that file's own x/y to under a meter
over 3.2 km, and altitude passes through exactly. What costs is differentiating twice.

The default is not the best row in the table — a dense 50 Hz track does better at 0.0656 ft,
which is tuned to *that* file's coordinate resolution. It is the best *default*, because it is
the only setting that is good on every input shape at once: 0.073 dB on dense positions,
0.121 dB on sparse waypoints, 0.015 dB when velocities are given. Exact interpolation is
uniformly catastrophic on measured input, and a value tuned for one sampling rate is wrong for
another. Tune it down if you know your coordinate resolution; leave it alone otherwise.

**The ends of a dense track are fitted freely.** A natural smoothing spline is pinned to zero
second derivative at the first and last knot, which is a statement that the aircraft is not
accelerating there. A track that starts or stops mid-manoeuvre is therefore fitted with a
velocity and a load factor that bend toward a straight line over the last few seconds: on a
constant-acceleration test path (1.5 m/s²) the natural fit's acceleration is zero at the end,
0.50 four seconds in and 1.14 eight seconds in, where the freed fit holds it within 11%
throughout. For positions sampled at
0.5 s or finer, with a nonzero `--path_smoothing`, the fit is therefore made over the track
extended by 10 s at each end along a local quadratic (fitted to the last 3 s), so that the
natural condition lands on the extension and not on the data. The extension is not used
afterwards. Sparser tracks, and `--path_smoothing 0`, keep the natural ends; there the
interpolant is the statement of intent, and a quadratic fitted to a few waypoints is not. On
the table above this takes the default row from 0.109 dB median and 84 cells over 1 dB to the
values shown, and leaves the exact and waypoint rows as they were. The shipped AC tracks carry
their own state and never go through this.

Sparse waypoints are resampled onto the `-s` grid first — never coarsened — and their load
factor is an artifact of the interpolant, which is logged.

> The `--path_smoothing 0` row was, until recently, reported as 0.196 dB. That measured a bug:
> `SmoothingSpline` short-circuited at λ = 0 and returned zero curvature, so the load factor
> came out as exactly 1 g on every state and the rounding amplification never happened. With
> the interpolating spline carrying its own curvature, exact interpolation costs what the
> surrounding text always said it would.

**Note on velocity columns.** `--path_smoothing` is a distance and applies to the fitted
position only. A file carrying `ve`/`vn`/`vu` has those interpolated exactly; they are
measurements of the quantity wanted, not something to be re-derived.

**Output.** `--latlon` writes a companion beside every level file (`out.csv` →
`out_latlon.csv`) holding latitude, longitude and level, one row per observer, for dropping a
footprint onto a map. It needs a frame: either a geodetic track, which carries one, or
`--origin`. A grid regular in the local frame is not regular in latitude and longitude, so this
file is always scattered, never gridded.

**The transform.** `src/geodetic.hpp` is the American polyconic projection about the reference
point (it agrees with a transverse Mercator to 3 mm at 10 km, and reproduces the shipped
tracks' own x/y), with altitude carried through as z — a map projection with a height attached rather than a true ENU
frame. Over the few kilometers a footprint spans the difference is far below the modeling
error. `from_east_north_up` applies the frame's rotation without its translation, for
velocities; `to_geodetic` inverts the projection by Newton iteration, and is what `--latlon`
uses.

The shipped AC tracks are the case where the axes are not east/north: their x/y run along a
compass bearing of **210°**, the approach heading, not along the wandering heading in the
file's own column 8. Their ground plane sits 54.2651 ft below their altitude column. So
```sh
./build/niceops -t track.csv --origin 39.8104276 -74.4263086 --frame_bearing 210 \
                --ground_elevation 54.2651
```
is the demo case restated in WGS84, and is what the first row of the table above measures
(`--path_smoothing 0.0656` gives the second).

### Plotting

Plotting embeds Python through pybind11. The first configure creates a virtual environment at
`build/pyenv` and installs matplotlib into it, and the embedded interpreter puts that directory
at the front of `sys.path` — so a run depends on what the build made rather than on whatever
happens to be installed for the system interpreter, and needs no `PYTHONPATH`. The install
happens once per build directory and needs network access; afterwards configures are offline.

Two options control this:

| Option | Meaning |
| --- | --- |
| `-DNICEOPS_PLOTTING=OFF` | Drop plotting, and Python, pybind11 and matplotlib with it. The acoustics need none of them. `-p`, `-f` and `-F` are then refused with a message rather than silently doing nothing. |
| `-DNICEOPS_BUNDLE_PYTHON=OFF` | Keep plotting but skip the virtual environment, using the interpreter's own packages instead. |

Headless is handled: the facade selects the `Agg` backend itself.

The reference figures under `figures/` date from June 2020 and were produced with three settings
that differ from the defaults. All three are now reachable:

```sh
./build/niceops -t tracks/S-76D/197401AC.csv -v databases/S-76D.json -d databases/S-76D_M3.nod \
        --trim -6000 --grid -2000 1100 -1500 1600 -o out.csv -f out.eps
```

Matched that way the two runs land on the same 124 × 124 grid and agree cell for cell: the largest
difference anywhere is **1.05 dB**, 99% of cells are within 0.5 dB, and the peak levels differ by
0.07 dB. Contour levels now default to 65–100, the values the references were drawn with; the old
45–90 put its top level below the peak of the demo case, leaving the loudest patch unfilled.

## License

NICE-OPS is under the MIT License; see `LICENSE`. Every third-party component it
builds against is permissive — BSD, MIT or the Boost Software License — and none carries a
reciprocal or source-disclosure obligation. `THIRD_PARTY_NOTICES.md` reproduces their notices,
which a binary distribution needs to carry with it.

## Repository layout

- `src/` — the implementation. Every `.cpp` lives here.
- `tests/` — Catch2 unit tests and the `golden_compare` tool.
- `golden/` — reference outputs captured from the D implementation; the regression tests compare
  against these at 1e-9 relative tolerance.
- `databases/`, `tracks/` — shipped noise databases and flight tracks.
- `atmospheres/` — atmosphere profiles for `--atmosphere`.
- `tools/` — `pe_compare.cpp`, the refraction models against the parabolic equation, and
  `benchmark.cpp`, the Attenborough et al. benchmark cases; with a plotting script for each.

`CLAUDE.md` carries the developer-facing detail: the port's history, the shape of the hot path,
and the pitfalls worth knowing before changing the numerics.

### Arrival completeness and PE sample counts

Stratified propagation retains all solved ray and wave arrivals. It previously
kept only the 24 arrivals with the smallest spreading ranges, which discarded
energy without accounting for source directivity or ground losses. The 10 dB
focusing cap applies per ray; distinct arrivals add in energy and can exceed it
in total. Deep ducts can therefore take longer to evaluate and predict more
energy. Launch-angle sampling and the bounce search remain approximate.

PE band sample counts must be positive. Zero or negative counts are rejected by
the solver API, including for an empty batch, rather than producing invalid band
levels. Valid sample counts retain the existing midpoint quadrature.
