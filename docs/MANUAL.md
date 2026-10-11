# NICE-OPS manual

This manual is the user's reference: how to run NICE-OPS, every option, the input and output
formats, which model to use when, and what each one costs and refuses. It states each model in a
paragraph or two. [THEORY.md](THEORY.md) has the equations and references,
[FIDELITY.md](FIDELITY.md) the evidence of how well each model agrees with exact solutions and
flight-test measurements, and [HISTORY.md](HISTORY.md) the changes that moved results. To build
the code, start with the [README](../README.md); developers start at
[internals/README.md](internals/README.md).

## Contents

- [Overview: how a prediction is made](#overview-how-a-prediction-is-made)
- [Running NICE-OPS](#running-nice-ops): [Quick start](#quick-start) · [Command-line
  options](#command-line-options) · [Units](#units) · [Errors and refusals](#errors-and-refusals)
- [Inputs](#inputs): [Input files](#input-files) · [Trajectories](#trajectories) · [Vehicle
  file](#vehicle-file) · [Observers](#observers)
- [Noise databases](#noise-databases): [Shipped databases](#shipped-databases) · [Spectra](#spectra) · [Screening a database](#screening-a-database)
  · [Choosing spheres: blending and `--blend`](#choosing-spheres-blending-and---blend) · [The
  load-factor law below 1 g](#the-load-factor-law-below-1-g) · [Speed reference](#speed-reference) ·
  [The sphere's nose](#the-spheres-nose) · [Reading the sphere's
  angles](#reading-the-spheres-angles)
- [Propagation](#propagation): [Broadband and spectral paths](#broadband-and-spectral-paths) ·
  [Stratified atmosphere](#stratified-atmosphere) · [A profile from the surface
  wind](#a-profile-from-the-surface-wind) · [Ground effect](#ground-effect) ·
  [Receivers](#receivers) · [Ground-plane microphones](#ground-plane-microphones) ·
  [Terrain](#terrain) · [Refraction](#refraction) · [Refraction through the profile's
  layers](#refraction-through-the-profiles-layers) · [Low sources near grazing: the resolved
  shadow](#low-sources-near-grazing-the-resolved-shadow) · [Micro-shadows and skip
  zones](#micro-shadows-and-skip-zones) · [Propagation loss tables](#propagation-loss-tables)
- [Metrics and output files](#metrics-and-output-files): [Output files](#output-files) · [Run records and option files](#run-records-and-option-files) · [Exposure
  (SEL) and the 10 dB-down window](#exposure-sel-and-the-10-db-down-window) · [Several trajectories,
  and the day-night level](#several-trajectories-and-the-day-night-level) · [Maximum level, time
  above and C-weighted SEL](#maximum-level-time-above-and-c-weighted-sel) · [EPNL and reception-time
  frames](#epnl-and-reception-time-frames)
- [Companion tools](#companion-tools): [`niceops_procedure`](#niceops_procedure) · [`niceops_propagation_table`](#niceops_propagation_table) ·
  [`niceops_pe_compare`](#niceops_pe_compare) · [`niceops_benchmark`](#niceops_benchmark) · [Plot
  scripts](#plot-scripts)
- [Performance](#performance): [How much gets evaluated](#how-much-gets-evaluated) ·
  [`--db_cutoff`](#--db_cutoff) · [Threads](#threads) · [Relative cost](#relative-cost)
- [Plotting](#plotting)
- [Comparing with older results](#comparing-with-older-results)

## Overview: how a prediction is made

NICE-OPS is a source–path–receiver model, run for every trajectory state and every observer: on
the demo case, 6,551 states by 14,400 observers.

- **Source.** A database of noise spheres around the aircraft, one per operating condition. Each
  state is reduced to the advance ratio μ, the flight path angle γ (to the tip-path plane under
  `--qsam`, the default, or to the horizontal) and the thrust coefficient C<sub>T</sub>. The
  velocity is taken in the database's recorded speed reference (the ground if it records none)
  unless `--speed_reference` overrides it ([Speed reference](#speed-reference)). The condition
  picks the spheres, and the direction to the observer the level on them
  ([THEORY §2](THEORY.md#2-the-source), [§3](THEORY.md#3-trajectories)).
- **Path.** The broadband path uses the database's own attenuation factor (EAA). The spectral
  path, the default for a database with spectra, carries each band with ISO 9613-1 absorption,
  the ground effect, and optionally a stratified atmosphere, refraction or a propagation table
  ([THEORY §4–§9](THEORY.md#4-spreading-and-atmospheric-absorption)).
- **Receiver.** Free field, a point receiver above the ground, or a flush microphone on a ground
  plate ([Receivers](#receivers), [THEORY §6](THEORY.md#6-receivers-and-terrain)).
- **Metric.** Each observer's history is integrated in reception time into an exposure, or read
  for its maximum, time above a threshold or EPNL; several trajectories combine into day-night
  metrics ([THEORY §10](THEORY.md#10-noise-metrics)).

## Running NICE-OPS

### Quick start

With no arguments NICE-OPS runs the demo case: `tracks/S-76D/197421AC.csv` against
`databases/S-76D_M3.nod` and `databases/S-76D.json`, writing `integrated/197421_LAE.csv`. Run from
the repository root, since paths are relative; output directories are created as needed.

```sh
./build/niceops                                     # the demo case (broadband)
./build/niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json \
                --humidity 70 -o integrated/197421_B407.csv   # spectral, with ground effect
```

A database with spectra runs the spectral path with the ground effect by default, which needs a
relative humidity (vehicle file, `--humidity` or `--atmosphere`). Each run prints a one-line
`Physics:` summary of the path, ground model and atmosphere it used.

### Command-line options

`niceops --help` prints the same list. Lengths are in `--units` (feet by default); a default
is a physical quantity, given here in meters with its value in feet.

#### Inputs and outputs

| Option | Meaning |
| --- | --- |
| `-t`, `--trajectoryfile` | Trajectory CSV ([Trajectories](#trajectories)). Repeat for several flight paths. Default the demo track. |
| `-d`, `--database` | netCDF noise database. Once, or once per trajectory. Default `databases/S-76D_M3.nod`. |
| `-v`, `--vehicle` | JSON vehicle file ([Vehicle file](#vehicle-file)). Once, or once per trajectory. Default `databases/S-76D.json`. |
| `-o`, `--output` | SEL output file, once per trajectory. Default `integrated/197421_LAE.csv`. With several tracks and one name (or none), each is named after its track. |
| `-b`, `--observer` | Scattered observer file ([Observers](#observers)); replaces the default grid. |
| `--grid` | Observer grid bounds: `xmin xmax ymin ymax`; the grid starts at each minimum and stops short of each maximum. Default −1500 … 1475 m in 25 m steps (120 values per axis); −4921.26 … 4839.24 ft. |
| `--gridstep` | Observer grid spacing. Default 25 m (82.021 ft). |
| `-S`, `--state` | Save the derived trajectory state, always metric. With several tracks and one name, `<track>_state.csv`. |
| `-i`, `--history` | Prefix for per-observer level histories, `<prefix>_i<N>.csv` (`<prefix>_<track>_i<N>.csv` with several tracks). |
| `--band_history` | Prefix for per-observer band histories, named as `-i` names its files. Needs the spectral path. |
| `-m`, `--metricfile` | Append the case's overall metrics to this file ([Output files](#output-files)). |
| `--powerfile` | Write per-state shaft power, when the database carries it. |
| `--latlon` | Also write each level file as latitude, longitude, level. Needs `--origin` unless the trajectory is geodetic. |
| `--origin` | WGS84 origin of the local frame: `latitude longitude`, degrees. Default the first position of the first geodetic track. |
| `--frame_bearing` | Compass bearing of the local +x axis, degrees. Default x east, y north. The shipped AC tracks use 210. Turns every compass quantity (wind, headings, table bearings) into the frame. |
| `--ground_elevation` | Elevation of the local z = 0 on the trajectory's vertical datum. Subtracted from a WGS84 track's altitude. Default 0. |
| `--path_smoothing` | How far, rms per axis, the path fitted through a WGS84 track may sit from its positions. Default 10 ft (3.048 m). Zero interpolates exactly. |
| `-s`, `--step` | Resampling step, s (default 0.5), for a WGS84 track whose sampling is uneven or coarser than this. A finer track is never coarsened, and an AC track is used as sampled. Not the SEL integration step. |
| `--trim` | Crop each trajectory from this x to its lowest or slowest point. Off by default; it changes which part of the flight is integrated. |
| `--start_utc` | Seconds since UTC midnight at *t* = 0 of a WGS84 track, for the day-night periods. Required with `-D`, `--lden_file`, `--cnel_file` or a period level when the track has no `utcsec` column. Once, or once per trajectory. |
| `--units` | `ft` (default) or `m`: the length unit of the command line and of every file ([Units](#units)). |
| `-l`, `--logfile` | Log file. Default `log/niceops.log`. |
| `--manifest` | Write a JSON record of the run ([Run records and option files](#run-records-and-option-files)). |
| `--config` | Read options from a TOML file; the command line overrides it. |
| `--version` | Print the build's version (its git state) and exit. |
| `-L`, `--loglevel` | `trace`, `info`, `warn`, `error` or `fatal`. Default `info`. |

#### Source and database

| Option | Meaning |
| --- | --- |
| `-q`, `--qsam` / `--no-qsam` | Take the flight path angle relative to the tip-path plane (Q-SAM), rather than to the horizontal. On by default. |
| `--blend` | How a state picks its spheres: `barycentric` (default), `nearest` or `aam` ([Choosing spheres](#choosing-spheres-blending-and---blend)). |
| `--max_load_factor` | Largest load factor, g, that a `fixed_load_factor` database's load-factor law is taken to. Default 2. Above it the correction is held, and the run counts the states. |
| `--load_factor_correction` | `one_sided` (default): the law applies above 1 g only. `two_sided`: below 1 g too ([The load-factor law below 1 g](#the-load-factor-law-below-1-g)). |
| `--min_load_factor` | With `two_sided`, the smallest load factor, g, the law is taken to. Default 0.5; at most 1. Refused without `two_sided`. |
| `--speed_reference` | `ground` or `air`: which velocity sets μ, the path angle and the `velocity` nose. Default: what the database records, ground if nothing ([Speed reference](#speed-reference)). |
| `--nose` | `velocity`, `ground_velocity` or `heading`. Default: what the database's `azimuth_reference` records, `velocity` if nothing ([The sphere's nose](#the-spheres-nose)). |
| `--accept_vehicle_mismatch` | Warn, rather than refuse, when the vehicle file disagrees with the database by more than 2 percent ([Vehicle file](#vehicle-file)). |
| `--screen_database` | Write every entry the load-time screen finds suspect to this CSV ([Screening a database](#screening-a-database)). One database per run. |
| `--mask_suspect_cells` | Drop the suspect entries, as gating drops a direction, instead of only warning. Not with `--self_test`. |

#### Atmosphere and wind

| Option | Meaning |
| --- | --- |
| `--spectral` / `--no-spectral` | On by default when the database carries spectra: propagate each band with ISO 9613-1 absorption. `--no-spectral` takes one broadband level with the database's EAA factor; a database without spectra always does. Asking for `--spectral`, or for any option that needs bands (`--ground_effect`, `--atmosphere`, `--wind`, `--propagation_table`, `--band_history`, `--csel_file`, `--epnl_file`, `--band_frames`), with a database that has no spectra is an error. |
| `--band_limits` | Frequency range summed into the metric, Hz: `low high`. Default `20 20000`. |
| `--humidity`, `--temperature`, `--pressure` | Relative humidity (percent), temperature (K) and pressure (kPa) for the spectral path, overriding the vehicle file. Not `--temperature` or `--humidity` with `--atmosphere`. `--pressure` only when the file has no pressure column; it is then the surface pressure.  Once, or once per trajectory. |
| `--atmosphere` | A stratified atmosphere: CSV of temperature, humidity, and optionally pressure and wind, against height ([Stratified atmosphere](#stratified-atmosphere)). Spectral path only.  Once, or once per trajectory. |
| `--wind` | Build the profile from one surface wind instead: `SPEED,DIRECTION[,HEIGHT]` ([A profile from the surface wind](#a-profile-from-the-surface-wind)). Not with `--atmosphere`. Spectral path only. The four rows below and `--write_atmosphere` are refused without it.  Once, or once per trajectory. |
| `--roughness_length` | `--wind`: aerodynamic roughness length z₀. Default 0.03 m (0.0984 ft), open flat terrain. Not `--roughness`. |
| `--surface_layer_depth` | `--wind`: height to which the wind follows its law; held above. Default 200 m (656.2 ft). |
| `--stability` | `--wind`: `neutral` (default), or a TA Luft class `I`, `II`, `III/1`, `III/2`, `IV`, `V` (Pasquill `F` to `A`). The classes are unvalidated, and a run that takes one warns. |
| `--profile_levels` | `--wind`: profile levels per factor of ten in height. Default 16. Results do not depend on it. |
| `--write_atmosphere` | `--wind`: write each trajectory's built profile in the `--atmosphere` format. With several tracks and one name, `<track>_atmosphere.csv`. |

#### Refraction and shadows

| Option | Meaning |
| --- | --- |
| `--refraction` | `none` (default), `linear` or `stratified`: curved rays through the `--atmosphere` or `--wind` profile ([Refraction](#refraction)). Needs a profile. |
| `--shadow_floor` | The most a refractive shadow attenuates, dB. Default 20. Past it, energy scattered in by turbulence is taken to dominate. |
| `--scattering` | Fill refraction shadows with the energy turbulence scatters into them, at HARMONOISE strength γ<sub>T</sub>, in place of the `--shadow_floor` hold. Default 0 (off). Needs `--refraction` and the ground effect. |
| `--shadow_model` | `--refraction stratified` with the ground effect: `legacy` (default) or `resolved` ([Low sources near grazing](#low-sources-near-grazing-the-resolved-shadow)). |
| `--shadow_phase` | `resolved`: the vertical phase, radians, a band must gather to resolve the profile near the ground. Default 1. |
| `--cramped_ground` | `resolved`: how a band whose boundary is within two diffraction lengths of the source is filled: `series` (default), `straight` or `ray`. |
| `--tunneling` | `--refraction stratified`: add the rays that turn back just above the microphone, through an Airy tail per band. Off by default ([Micro-shadows](#micro-shadows-and-skip-zones)). |
| `--skip_fill` | `--refraction stratified`: fill the skip zone short of where rays from aloft first land with the caustic's Airy tail. Off by default. |

#### Ground, receivers and terrain

| Option | Meaning |
| --- | --- |
| `--ground_effect` / `--no-ground_effect` | The ground reflection, on each band ([Ground effect](#ground-effect)). On by default wherever the spectral path runs. Asking for it with `--no-spectral` is an error. |
| `--flow_resistance` | Ground flow resistivity, kPa·s/m². Default 225 (grass). 30 snow, 650 soil, 3000 dirt, 10000 concrete, 50000 asphalt. |
| `--impedance_model` | `delany_bazley` (default), `miki`, `variable_porosity` or `delany_bazley_layer`. |
| `--porosity_rate` | `variable_porosity`: the porosity's rate of change with depth, 1/m whatever `--units` says. Default 0. |
| `--layer_depth` | `delany_bazley_layer`: depth of the hard-backed layer. Default 0.05 m. |
| `--roughness` | RMS roughness height of the ground; the coherent reflection falls by exp(−2(kσ cos θ)²). Default 0. |
| `--turbulence_model` | `chessell` (default, with `--turbulence`) or `harmonoise` (with `--gamma_t`): how turbulence decorrelates the direct and reflected paths. |
| `--turbulence` | Chessell's turbulence parameter, rad·s/√m, typically 0 to 16e-4. Default 0. |
| `--gamma_t` | HARMONOISE turbulence strength γ<sub>T</sub>, m<sup>−2/3</sup>: about 1e-6 moderate, 1e-5 strong. Default 0. |
| `--receiver_height` | Height of a point receiver's microphone. Default 1.2 m (3.937 ft). An observer file's fifth column overrides it per observer. Also places the microphone for `--propagation_table`. |
| `--plate_table` | Observers are flush microphones on a ground plate: a netCDF table built for the run's ground ([Ground-plane microphones](#ground-plane-microphones)). Needs the ground effect. |
| `--plate_mic_azimuth` | Which way the microphone is offset on its plate: `outboard` (default, away from the track) or degrees in the local frame from +x toward +y. |
| `--exact_ground` | Evaluate the ground model (or the plate's reflection coefficients) at every observer and state rather than from per-trajectory tables. Slower; for checking the tables. |
| `--propagation_table` | A propagation loss table (netCDF) that replaces the ground and ray models wherever it covers a pair ([Propagation loss tables](#propagation-loss-tables)). Spectral path only. |
| `--terrain` | A bare-earth elevation grid (netCDF, for example USGS 3DEP) ([Terrain](#terrain)). Needs the run's frame: `--origin` (with `--frame_bearing` for rotated axes; the AC tracks use 210), or a geodetic track. Not with `--refraction` or `--propagation_table`. The three rows below are refused without it. |
| `--terrain_offset` | Added to every terrain height to put the grid on the trajectory's vertical datum. Default 0. |
| `--terrain_step` | Spacing of the terrain samples each path's mean ground plane is fitted to. Default the grid's own spacing, at least 5 m. |
| `--observers_on_terrain` | With `--terrain` and `-b`: take each observer's ground height from the terrain rather than the file's z. |

#### Metrics

| Option | Meaning |
| --- | --- |
| `-a`, `--ambient` | Ambient level, dB, added to every sample's level in energy. Off unless given. Refused with `--csel_file`, `--epnl_file`, `--band_frames` and `--coverage_file`. |
| `--full_history` | Integrate every exposure (SEL, C-SEL, and each operation's share of DNL, Lden, CNEL and the period levels) over the whole history rather than the 10 dB-down window. |
| `-D`, `--dnlfile` | Write the day-night level, summed over every trajectory ([Several trajectories](#several-trajectories-and-the-day-night-level)). |
| `--lden_file` | Write the day-evening-night level (EU Directive 2002/49/EC). |
| `--lday_file`, `--levening_file`, `--lnight_file` | Write Lden's period levels: the energy heard 07:00–19:00, 19:00–23:00 or 23:00–07:00, over that period's length ([Several trajectories, and the day-night level](#several-trajectories-and-the-day-night-level)). |
| `--laeq_file` | Write LAeq,T over the local clock window `--laeq_window START END` (hours; `7 23` is LAeq,16h, `23 7` wraps). Needs the window, which has no default. |
| `--cnel_file` | Write the community noise equivalent level (California Title 21). |
| `--utc_offset` | Hours from UTC to local time, for the day-night periods. Default 0. |
| `--lamax_file` | Write each observer's LAmax over every operation, time-weighted per `--lamax_weighting`. |
| `--lamax_weighting` | `fast` (default), `slow`, `leq0.5` or `instantaneous` ([Maximum level](#maximum-level-time-above-and-c-weighted-sel)). Refused without `--lamax_file`, `--nax_file` or `--ta_file`. |
| `--nax_file`, `--nax_threshold` | Write each observer's NAx: the number of operations whose own LAmax reaches the threshold, dBA. The threshold has no default and is required. |
| `--ta_file`, `--ta_threshold` | Write each observer's time above: seconds of reception time the weighted level exceeds the threshold, dBA, summed over operations. The threshold has no default and is required. |
| `--csel_file` | Write the C-weighted SEL (LCE) per observer. One name, or one per trajectory (`<track>_LCE.csv`). Spectral path only; refused with `--db_cutoff` or `-a`. |
| `--epnl_file` | Write each observer's EPNL, PNLTM, duration and tone corrections (14 CFR 36 A36.4). Once, or once per trajectory (`<track>_EPNL.csv`). Needs the spectral path and the 24 bands 50 Hz–10 kHz. Refused with `-a`. See [EPNL](#epnl-and-reception-time-frames). |
| `--epnl_masking` / `--no-epnl_masking` | Exclude tones in bands below their noy threshold from the tone correction. On by default; `--no-epnl_masking` is the procedure as written. |
| `--band_frames` | Write every observer's band levels in 0.5 s reception-time frames, as netCDF. Once, or once per trajectory (`<track>_band_frames.nc`). Needs the spectral path. Refused with `-a`. |
| `--coverage_file` | Write, per observer, the energy-weighted share of its level that came from sphere cells the database did not measure. Needs a database that records `coverage`. Straight rays only: refused with `--refraction` or `-a`. |
| `--provenance_spheres` | With `--provenance_file`, also write every sphere's share of every observer's SEL energy to a netCDF file ([Output files](#output-files)). |
| `--provenance_file` | Write, per observer, where its SEL energy came from: the shares from states clamped on each database axis or held at a load-factor limit, and the sphere that supplied the most ([Output files](#output-files)). One database per run; refused with `-a`. |

#### Threads, culling and diagnostics

| Option | Meaning |
| --- | --- |
| `-j`, `--threads` | Worker threads for the observer loop. 0 (default) uses one per hardware thread; 1 runs serially. The output is the same at any count. |
| `--db_cutoff` | Skip evaluations that provably fall below this level, dB. Off by default; read [the caveat](#--db_cutoff). |
| `--export_sphere` | Interpolate a whole noise sphere at `--sphere_condition` and write it as a one-group `.nod` database. |
| `--export_aam` | The same sphere in AAM/RNM format. Needs a database with spectra. |
| `--sphere_condition` | The condition to export at: advance ratio, flight path angle, C<sub>T</sub>. `--blend` and `--load_factor_correction` apply. |
| `--self_test` | Check every channel of the database against what the run's evaluation returns for it, then exit. |
| `-p`, `--plotting` | Plot the footprint with the flight track over it ([Plotting](#plotting)). |
| `-f`, `--figfile` | Save the footprint figure to this file. |
| `-F`, `--statefigfile` | Save the trajectory-state figures to this file. |
| `--cmin`, `--cmax`, `--cstep` | Contour levels of the footprint figure and of `--contour_file`, dB. Defaults 65, 100 and 5. |
| `--contour_file` | Write contour polygons and their areas as GeoJSON ([Output files](#output-files)). Needs the observer grid and a frame. |
| `--raster_file` | Write every gridded level of the run to one CF netCDF a GIS opens ([Output files](#output-files)). Needs the observer grid. |
| `--contour_metric` | What `--contour_file` contours: `sel` (default; one trajectory), `dnl`, `lden`, `cnel`, `lday`, `levening`, `lnight`, `laeq` or `lamax`, whose own output must be asked for too. |

### Units

`--units` selects the length unit of the whole interface, `ft` (the default) or `m`: every length
on the command line (`--grid`, `--gridstep`, `--trim`, `--ground_elevation`, `--receiver_height`,
`--path_smoothing`, `--terrain_offset`, `--terrain_step`, `--roughness`, `--layer_depth`,
`--roughness_length`, `--surface_layer_depth` and the anemometer height in `--wind`), every file
carrying a length (trajectories, observers, `--atmosphere` heights, level output, histories),
and the figure axes.

**`--units` says what your files already are; it does not convert them.** The shipped tracks
and observer sets are in feet, so running them with `--units m` reads feet as meters and flies
the aircraft three times too far, which nothing can detect. Where a file names its unit in a
column header (`alt_ft`, `alt_m`, `z_ft`, `Alt (ft)`), the name is checked against `--units` and
a disagreement fails the file; a bare `alt` claims nothing and is accepted.

- **Defaults are physical quantities**, the same ground in either unit: the observer grid is
  −1500 … 1475 m, which reads as −4921.26 … 4839.24 ft. Scaled with the flag, the default grid
  would change size by a factor of 3.28 between `--units ft` and `--units m`.
- **The `-S` state file is metric**, a dump of the internal state.
- **Speeds name their own unit** (`--wind 12kt,270`, `wind_speed_kt`); a bare number is refused.
- **Other quantities keep the units of their standards**: K, kPa (pressure and flow
  resistivity, as ISO 9613-1 is written), Hz, degrees, seconds, dB. The `--terrain` grid and a
  propagation table are in meters whatever `--units` says.

### Errors and refusals

A bad input must fail or warn with its cause named; it should never give a plausible footprint
with exit status 0. Errors and warnings go to stderr as well as to the log file (`-L` decides
what is written at all). A trajectory that cannot be read, or whose inputs are refused, is
skipped with `Skipping <file>: <why>` and the run goes on with the others; the exit status is 1
if anything was skipped. A refusal that concerns the whole run (an option, the observers, the
first database) stops it with `Stopping: <why>` and exit status 2.

**Files.** An AC track's time must increase and its state columns be finite. A WGS84 track's
time, altitude and velocity must be finite, its column names recognized, and any unit its header
declares must agree with what is read: lengths in `--units`, speeds in `--units` per second,
times in seconds (`Alt (ft)` under `--units m`, `Climb Rate (fpm)` and `Time (ms)` are
refused). A text line after the first row of data is an error (above it, a header). A track whose states
go below the ground plane, with no `--terrain`, is warned about (beyond the larger of 1 m and a
WGS84 track's `--path_smoothing`, and after `--trim`): its altitude datum and `--ground_elevation`
disagree, or an AC track's z is not height above the ground. A `--trim` that leaves fewer than
two states is refused, naming where the track's crop ends. An observer file must keep one
width and yield at least one observer. An `--atmosphere` file must be finite, reach every state,
and name its wind units. A database must vary in μ, γ and (unless `fixed_load_factor`)
C<sub>T</sub>, with no two spheres at one condition, and no group may hold the netCDF fill value
(a group the builder never wrote). Air must be plausible: a temperature in °C,
a pressure in Pa or a humidity fraction is refused. The vehicle file must agree with the
database ([Vehicle file](#vehicle-file)).

**Combinations** that cannot be served are refused before any level is computed: those on the
command line alone stop the run, and those that depend on a trajectory or its database (no
humidity, no clock, a coarse track, `--speed_reference air` without wind) skip that trajectory.
The main ones:

| Refused | Why |
| --- | --- |
| A spectral run with no humidity from any source | Absorption depends strongly on it; there is no silent default. |
| `--spectral`, or any option that needs bands, on a database without spectra | There are no bands to carry. |
| `--ground_effect`, `--atmosphere`, `--wind`, `--propagation_table`, `--band_history`, `--csel_file`, `--epnl_file`, `--band_frames` with `--no-spectral` | The broadband path carries one A-weighted level. |
| `--csel_file`, `--epnl_file`, `--band_frames`, `--coverage_file` with `-a` | The ambient is an A-weighted level with no spectrum. |
| `--provenance_file` with `-a` | The ambient's energy came from no state or sphere. |
| `--provenance_spheres` without `--provenance_file` | It is a companion to that file, and takes its limits. |
| `--csel_file` with `--db_cutoff` | The cutoff is proven on the A-weighted level. |
| `--refraction` without `--atmosphere` or `--wind`; `--atmosphere` with `--wind` | One profile, stated once. |
| `--shadow_model resolved` without `--refraction stratified` and the ground effect; `--scattering` without `--refraction` and the ground effect | They act through the shadow transition, which needs a ground model. |
| `--tunneling`, `--skip_fill` without `--refraction stratified` | They are stratified-ray corrections. |
| `--terrain` with `--refraction` or `--propagation_table`, or without a frame | Those models are built over flat ground; terrain needs to know where x/y lie. |
| `--plate_table` without the ground effect, or over a different ground | The table is built for one ground. |
| `--latlon`, `--contour_file`, `--nose heading` without a frame; `--contour_file` or `--raster_file` with `-b`; `--contour_file` with `--contour_metric sel` over several trajectories | Nothing to turn compass directions or invert positions with; scattered observers have no grid to contour; one SEL per trajectory. |
| `-D`, `--lden_file`, `--cnel_file` or a period level on a track with no clock | Day or night would be made up. |
| `--laeq_file` without `--laeq_window`, or an empty window | There is no default period. |
| A `--lamax_weighting` on a track sampled coarser than half its time constant | The "maximum" would be the largest sample. |
| `--speed_reference air` without a wind | Calm air is not assumed. |
| Two trajectories writing one file; a list of `-d`, `-v`, `-o` that is neither one nor one per track | Outputs would overwrite one another. |
| `--mask_suspect_cells` with `--self_test`; `--screen_database` or `--provenance_file` with more than one database; `--screen_database` or `--manifest` naming an input | The self test checks the file as written; one list would overwrite the other, and a sphere index names a sphere of one database; the input would be overwritten. |
| `-p`, `-f`, `-F` in a build without plotting | Rather than silently doing nothing. |

## Inputs

### Input files

| What | Format |
| --- | --- |
| Trajectory (`-t`) | 27-column NASA/FAA AC state CSV, a WGS84 CSV, or a list of static holds ([Trajectories](#trajectories)). Lengths in `--units`. |
| Database (`-d`) | netCDF `.nod`: one group per operating condition (μ, γ, C<sub>T</sub>), each a sphere of scattered (φ, θ) channels of A-weighted level and EAA, and optionally spectra ([Noise databases](#noise-databases)). |
| Vehicle (`-v`) | JSON ([Vehicle file](#vehicle-file)). |
| Observers (`-b`) | CSV of observer positions in `--units` ([Observers](#observers)). Without it, the default grid. |
| Atmosphere (`--atmosphere`) | CSV profile ([Stratified atmosphere](#stratified-atmosphere)). |
| Terrain (`--terrain`) | netCDF: 1-D `latitude` and `longitude` (degrees, strictly increasing) and `elevation(latitude, longitude)` in meters whatever `--units` says; an optional global `vertical_datum` text attribute. Missing nodes are NaN, `_FillValue` or `missing_value` ([Terrain](#terrain)). |
| Plate table (`--plate_table`) | netCDF from panam ([Ground-plane microphones](#ground-plane-microphones)). |
| Propagation table (`--propagation_table`) | netCDF ([Propagation loss tables](#propagation-loss-tables)). |

### Trajectories

`-t` reads three kinds of CSV, told apart by their shape: 27 columns is the positional AC
format, and anything else is read by the names in its header.

**AC state files** are the NASA/FAA flight-test format of the shipped tracks. Column 0 is UTC
seconds and column 1 elapsed time; lengths are in `--units`. Each state carries its own
velocity and acceleration, and the track is used as sampled. The shipped tracks' x/y run along a
compass bearing of **210°**, so they need `--frame_bearing 210` wherever a compass direction
enters (wind, headings, terrain, tables). Their ground plane is 54.2651 ft below their altitude
column.

**WGS84 tracks** give what a mission planner has: a route, not a flight-test record.

| Column | Meaning |
| --- | --- |
| `time` | Seconds from the start. `utcsec` (seconds since UTC midnight) serves instead and sets the clock for the day-night metrics. |
| `lat`, `lon` | Degrees. |
| `alt` | Altitude in `--units`, on the datum of `--ground_elevation`; z is the altitude less `--ground_elevation`. |
| `vE`, `vN`, `vUp` | Optional velocity, `--units` per second, used when all three are present. |
| `heading` | Optional compass heading of the nose, degrees from true north (`hdg` serves too); a unit other than degrees, a header naming a magnetic heading, and a second heading column are refused. It sets the hover nose and `--nose heading` ([The sphere's nose](#the-spheres-nose)). Resampled linearly, through north the short way. |

Names are matched case-insensitively, ignoring punctuation and a parenthesized unit (`Latitude`,
`LAT` and `lat (deg)` are one column), in any order. `alt_ft` is read under `--units ft` and
`alt_m` under `--units m`; a column naming the other unit is refused.

- **The frame.** `--origin <lat> <lon>` fixes the local origin (default: the first position of
  the first geodetic track), and `--frame_bearing` turns +x onto a compass bearing (default
  east/north). One frame serves the whole run.
- **Velocity and load factor** come from the first and second derivatives of a smoothing spline
  through the positions. Velocity columns, when given, are interpolated exactly instead.
- **`--path_smoothing`** is how far, rms per axis, the fit may sit from the positions (default
  10 ft); the fit uses that allowance. Zero is wrong for both kinds of input: on measured
  positions the load factor tracks coordinate rounding, and through sparse waypoints it
  manufactures maneuvering. The default works on dense positions, sparse waypoints and velocity
  columns alike ([FIDELITY](FIDELITY.md#trajectory-smoothing)); lower it only if you know your
  coordinate resolution.
- **Resampling.** A track uneven or coarser than `-s` (default 0.5 s) is resampled onto it, never
  coarsened. The load factor of sparse waypoints is an artifact of the interpolant, and the log
  says so.

The demo case restated in WGS84 is

```sh
./build/niceops -t track.csv --origin 39.8104276 -74.4263086 --frame_bearing 210 \
                --ground_elevation 54.2651
```

**Static operations.** A CSV whose header names `duration`, `x`, `y`, `z` and `heading` is a list
of holds at one position (in `--units`), each `duration` seconds long with the nose on the compass
`heading`. The states are sampled at 50 Hz, at 1 g and zero velocity.

### Vehicle file

JSON; values may be numbers or numeric strings. Other keys (the tail rotor, blade counts, drag)
are ignored. For example:

```json
{ "Main Rotor": { "radius": "6.705", "tip speed": "220.0" },
  "Vehicle": { "weight": "5000" },
  "Atmosphere": { "density": "1.225", "temperature": "300.71",
                  "pressure": "105.74", "relative humidity": "20" } }
```

| Key | Meaning |
| --- | --- |
| `Main Rotor.radius`, `Main Rotor.tip speed` | m and m/s. Required. |
| `Vehicle.weight` | Mass, kg. Required. |
| `Atmosphere.density` | kg/m³, for C<sub>T</sub>. Required. |
| `Atmosphere.temperature` | K. Optional; 293.15 K if absent. |
| `Atmosphere.pressure` | kPa. Optional; if absent, ρRT from the file's own density and temperature, not a standard day, since a file's density and temperature need not be an ISA pair. |
| `Atmosphere.relative humidity` | Percent. Optional, but a spectral run needs a humidity from some source. |

C<sub>T</sub> always takes `density`. The temperature, pressure and humidity are the spectral
path's air, which `--temperature`, `--pressure` and `--humidity` override.

**Checked against the database.** The tip speed, rotor radius and weight, and for a
`fixed_load_factor` database the thrust coefficient at 1 g, index the database's conditions, so
they are compared with what the database records. Within 1 percent is silent, 1–2 percent warns,
and past 2 percent the run is refused unless `--accept_vehicle_mismatch` is given. Vehicle files
ship for the AW139, B205, B407 and S-76D only; the AS350, Be206L3, EC130, RO-44 and RO-66
databases need one you supply.

### Observers

Without `-b` the observers are a grid, by default −1500 … 1475 m in 25 m steps on each axis (120
values, 14,400 observers), set by `--grid` and `--gridstep`. Grid observers stand at z = 0, or on
the terrain. `-b FILE` reads scattered observers instead, in `--units`: `x,y,z`, `id,x,y,z`, or
`id,x,y,z,height`, the fifth column being each microphone's height above the ground
([Receivers](#receivers)). The first line is always a header, and is discarded even when it
holds numbers. `tests/data/mixed_receivers.csv` is an example (coordinates shortened here):

```csv
id,x_ft,y_ft,z_ft,height_ft
2,-788.3,-1085.1,-1.1,0
10,-405.9,1940.8,-10.1,3.937
```

## Noise databases

A database (`.nod`, netCDF) holds one group per operating condition: advance ratio μ, flight
path angle γ and thrust coefficient C<sub>T</sub>. Each group stores one sphere around the rotor:
scattered (φ, θ) channels of the overall A-weighted level and of the EAA excess-attenuation
factor, and optionally per-band spectra. A `fixed_load_factor` database stores 1 g spheres only
and carries the load factor by an analytic law. Databases are built by panam, from measured or
computed hemispheres. [THEORY §2](THEORY.md#2-the-source) describes how they are read.

### Shipped databases

| Files in `databases/` | Contents |
| --- | --- |
| `AS350_`, `B407_`, `Be206L3_`, `EC130_`, `RO-44_`, `RO-66_fixed_load_factor.nod` | **Spectral**, 10 Hz–10 kHz, built by panam from the 2017 Noise Abatement flight test. The current, validated databases. |
| `AS350_`, `B407_`, `Be206L3_`, `EC130_`, `RO-44_`, `RO-66_ambient_gated.nod` | Broadband versions of the same builds. |
| `AW139_`, `B205_`, `S-76D_fixed_load_factor.nod` | Broadband, converted from the `_M*` files below. |
| `AW139_M1`, `B205_M1`, `B407_M2`, `B407_M3`, `S-76D_M3`, `S-76D_M4.nod` | Broadband, with discrete load factors, from joint NASA/FAA flight tests. `S-76D_M3.nod` is the demo database. |
| `AW139.json`, `B205.json`, `B407.json`, `S-76D.json` | Vehicle files. |

Of the 2017 databases, **B407 is the reference** and the most validated; its spheres are
depropagated along curved rays through each run's own atmosphere. Its hover is synthesized from
the slowest level pass and corrected by hovers measured 0.5–3° and 18–25° below the horizon, so
treat the hover from 3° to 17.5° and under the aircraft as an interpolation. The other five take
the B407's recipe without a hover correction. No approach was used to build any of them: the
approaches are the validation set ([FIDELITY](FIDELITY.md#shipped-databases-approach-scores-per-aircraft)).
Each file records its build at the root (`speed_reference`, the build atmosphere, `hover_source`)
and `DOPPLER_SHIFT_REMOVED` and `coverage` per group. The spheres keep the measured,
Doppler-shifted spectra, and NICE-OPS applies no frequency shift. The rebuilds are in
[HISTORY](HISTORY.md#the-2017-noise-abatement-databases).

### Spectra

A database may carry per-band spectra beside the overall A-weighted level: a `frequency`
dimension and an `amplitude(PHI, THETA, frequency)` array per sphere. It then runs the spectral
path by default, and serves every option that works on bands. `--export_aam` writes a sphere in
the AAM/RNM hemisphere format, whose payload is `AMPLITUDE(PHI, THETA, FREQUENCY)`, so it needs
spectra. Of the shipped databases, the six `_fixed_load_factor` files of the 2017 aircraft carry
spectra; the S-76D demo database does not, so the demo case runs broadband.

### Screening a database

Every database is screened as it loads, for entries no sound database holds, which would
otherwise be blended as real levels:

- **a dBA cell far from the rest of its sphere**: more than 40 dB from the median of that sphere's
  cells, whether a spike or a near-silent cell;
- **a band far above its cell's spectrum**: more than 45 dB above the median of that cell's bands.
  A band far *below* the rest is what ambient gating leaves, and is not flagged.

The screen warns once per kind, with the count, the spheres affected and the worst entry, and
otherwise changes nothing. `--screen_database FILE` lists every entry it found: sphere, its
condition, φ, θ, band (blank for the overall level), the level and the median it was judged
against. `--mask_suspect_cells` drops them, so the blends read those entries as unmeasured and take
the other conditions' level there. Masking is a judgment about data: rebuilding the database is
the fix, and a masked run should say so.

Four shipped databases warn: `B407_M2`, `B407_M3`, `S-76D_M4` and `S-76D_fixed_load_factor`. The
demo database and the 2017 databases do not. A sphere more than half defective has a defective
median, and the screen cannot judge it. The counts are in
[FIDELITY](FIDELITY.md#database-screen). The screen costs well under a second at load and nothing
afterward.

### Choosing spheres: blending and `--blend`

By default (`--blend barycentric`) the database's conditions are normalized into the unit cube and
triangulated, and a state is blended in dB from the spheres at the corners of the tetrahedron that
holds it. A condition outside an axis's range is clamped to that axis's end, and the envelope is
padded so that the clamped edge falls inside the hull; a query still outside takes the nearest
vertex. Each run logs how many states were clamped on each axis and how many took the load-factor
correction, and warns if any were clamped. Neither count says which directions were measured.

Two other policies exist for comparing with other codes and with regulatory practice. Neither
predicts better ([FIDELITY](FIDELITY.md#source-options-blend-and-load-factor-law)); they show how
much of a difference between two models is the choice of sphere. **`nearest`** takes the one
stored condition nearest the state, whole, in the normalized space: a trajectory tool's version of
ECAC Doc 32's one hemisphere per segment. **`aam`** is AAM v3.1's rotorcraft search: up to four
spheres bracketing the state in speed (within 20 kt) and then path angle (within 3° descending,
5° otherwise), or the nearest sphere where nothing brackets it; the run logs how many states
blended four, three, two and one sphere. The rules are in [THEORY §2.3](THEORY.md#23-blending-the-spheres).

### The load-factor law below 1 g

A `fixed_load_factor` database carries the load factor by the law 20 lg(C<sub>T</sub> /
C<sub>T,1g</sub>) dB, up to `--max_load_factor` (default 2 g). By default the law applies above
1 g only (`one_sided`), as panam applies it, so the sub-1 g states of an approach's flare and
pushovers read the 1 g sphere unscaled. `--load_factor_correction two_sided` applies it below 1 g
too, held below `--min_load_factor` (default 0.5 g), where the logarithm would take the level tens
of dB down; the run warns with the number of states held. On a database whose load factors are
spheres there is no law, and `two_sided` is refused. The evidence for the law below 1 g is modest
and mixed ([FIDELITY](FIDELITY.md#source-options-blend-and-load-factor-law)), so the default
stays one-sided.

### Speed reference

The flight condition (μ, the path angle and, under `--nose velocity`, the nose) comes from a
velocity. `--speed_reference ground` takes it over the ground; `air` takes the ground velocity
less the profile's wind at each state's height. The load factor, and so C<sub>T</sub>, stays
inertial either way.

- **The database says which it was built with**, in a root attribute `speed_reference`, and the
  run follows it; a file without one is ground, and the 2017 databases record `ground`. Giving
  `--speed_reference` against the database is logged as an override.
- **Air needs a wind.** With no wind in `--atmosphere` and no `--wind`, an air-referenced run is
  refused rather than taken as calm. A `--wind` of zero is calm stated.
- **It is applied after `--trim`**, and calm air changes nothing, to the bit.

Matching a ground-built database with an air-referenced run mixes the two references, so change
both together. On a 15 kt day air reference moves μ by 5–10% and the nose by 8–15°
([FIDELITY](FIDELITY.md#speed-reference)).

### The sphere's nose

The sphere's azimuth is measured from its nose, so where the nose points decides which side of
the sphere each observer reads. `--nose` chooses it:

- **`velocity`**: the velocity `--speed_reference` selects, projected into the tip-path plane.
  At hover the velocity is drift, so below 3 kt the nose is the heading, turning into the
  velocity by 8 kt.
- **`ground_velocity`**: the ground velocity whatever `--speed_reference` says, with the same
  hover turn, so μ can be taken through the air while the nose stays on the ground track.
- **`heading`**: the compass heading at every speed. In a crosswind a helicopter flies crabbed
  into the wind while its track runs straight, and its directivity turns with the airframe. Only
  the nose changes; μ, the path angles and C<sub>T</sub> still come from the velocity.

The default is what the database's root attribute `azimuth_reference` records: `heading` gives
`heading`, and `track` or nothing gives `velocity` (of the shipped databases only the AS350's
record it, as `track`); any other value is refused. An explicit
`--nose` wins, with a warning when the file records the other; `--export_sphere` carries the
attribute on. **The frame must be known; it is never inferred**: with no frame (an AC track
without `--frame_bearing`, and no `--origin`) `--nose heading` is refused. A track with no heading
(a WGS84 track without a `heading` column) is refused under `--nose heading`, and falls back to the
velocity with a warning when only the database asks for heading mode; a state without a heading is refused,
naming its time. The log gives the mean crab (heading less ground track, positive with the nose
right of the track) and its rms above 8 kt, and warns past 45°, which is a frame turned the wrong
way rather than a crosswind.

On the windy day of the 2017 test, where the aircraft crabbed 4–19° into the wind, `--nose
heading` lowered the SEL rms by about 0.15 dB, and on the calm day it changed almost nothing
([FIDELITY](FIDELITY.md#the-spheres-nose); the hover rule is in
[THEORY §2.4](THEORY.md#24-reading-the-sphere)).

### Reading the sphere's angles

The direction from the source to the observer gives two angles on the sphere, which
`--band_history` writes as `phi` and `theta`:

- θ, from the nose: 0° straight ahead, 180° straight behind.
- φ, around the nose axis: 0° straight down along the tip-path-plane normal, +90° to starboard
  and −90° to port in the tip-path plane, ±180° straight up.

The level is interpolated bilinearly between the sphere's nodes. Unless every sphere stores its
upper half, a direction above the tip-path plane reads its mirror image below.

## Propagation

### Broadband and spectral paths

The sphere gives the level at its surface, radius R<sub>s</sub>, and it is carried over the range
r ≥ R<sub>s</sub> to the microphone. The **broadband** path (a database without spectra, or
`--no-spectral`) applies spherical spreading and the database's EAA, its A-weighted attenuation
per 1000 m; no ground, atmosphere or refraction model applies. The **spectral** path spreads each
band the same way, absorbs it at its own ISO 9613-1 coefficient, applies the ground, receiver,
refraction and table models per band, then A-weights the bands within `--band_limits` and sums
them in energy. The two answer different questions and will not agree:

| | `--no-spectral` | spectral (the default with spectra) |
| --- | --- | --- |
| Attenuation | the database's EAA, one dB/km factor | ISO 9613-1 α(f), per band |
| Fitted for | the spectrum the source happened to have | the atmosphere you state |
| Atmosphere | none: EAA is baked in | temperature, pressure, humidity |
| Needs | any database | one carrying spectra |

**Band-limit deliberately.** Summing whatever bands a database carries makes the answer a
property of its bandwidth. `--band_limits low high` states the range; the run logs how many bands
survived, and a range with none is an error. **The air** is the vehicle file's
([Vehicle file](#vehicle-file)), overridden by `--humidity`, `--temperature` and `--pressure`,
or replaced by a profile.

### Stratified atmosphere

`--atmosphere file.csv` replaces the single atmosphere with a profile, and integrates each band's
absorption through it: the coefficient's mean over the heights the path spans, times the path
length ([THEORY §4.3](THEORY.md#43-atmospheric-absorption)). Columns are found by name, in any
order, one row per level; rows may run top down:

| Quantity | Column names | Unit |
| --- | --- | --- |
| Height above the observers | `z`, `agl` | `--units`; `z_ft` / `z_m` must agree with it |
| or altitude above sea level | `alt`, `altitude`, `height` | `--units`, less `--ground_elevation` |
| Temperature | `T_K`, `T_C` | named; a bare `T` is told from the values |
| Relative humidity | `RH`, `humidity` | percent |
| Pressure (optional) | `p_kPa`, `p_hPa`, `p_Pa` | named; a bare `p` is told from the values |
| Wind (optional) | `wind_speed_mps` or `wind_speed_kt`, with `wind_direction` | compass bearing it blows **from**, degrees |
| or wind components | `u_mps` and `v_mps`, or `u_kt` and `v_kt` | toward east and toward north |

```csv
z_ft,T_K,RH,p_kPa
0,288.150,20,101.325
1000,286.169,20,97.717
```

- **Pressure.** Without a column, it is integrated hydrostatically up from the surface pressure
  (the vehicle file's or `--pressure`'s). A single row is a uniform atmosphere.
- **The profile has to reach the aircraft.** A trajectory above the top level is refused. Below
  the lowest level that level is held to the ground, so give the profile a level at the ground.
- **Surface properties stay at the surface**: the ground level supplies the sound speed for the
  reception time and the ground impedance.
- **Not with `--temperature` or `--humidity`**; `--pressure` only when the file has no pressure.
- **Wind columns must name their unit.** Wind matters only to `--refraction` and
  `--speed_reference air`. The log names every column it ignored.

`atmospheres/` holds `isa.csv` (the ISA to 10,000 ft at 20% RH), `example_windy.csv` (an
illustrative log-law westerly), `example_inversion.csv` (a morning inversion),
`example_windy_calm_floor.csv` (the windy profile over a 0.1 ft calm floor) and `B407_uniform.csv`
(`B407.json`'s own atmosphere, which reproduces the homogeneous run).

### A profile from the surface wind

A measured profile (LIDAR, a balloon, a met tower) is rare; a surface wind, from a METAR or a
station beside the microphones, is not. `--wind SPEED,DIRECTION[,HEIGHT]` builds the profile from
that one observation in place of an `--atmosphere` file: the practical route to refraction on a
windy day without a sounding.

```sh
niceops ... --humidity 25 --wind 12kt,270,6.56 --refraction stratified --frame_bearing 210
```

- **The speed names its unit**: `12kt` (or `kts`, `knots`), `6.2mps` (`m/s`), `14mph`.
- **The direction is the compass bearing the wind blows from**, degrees true: a METAR's 27012KT
  is `12kt,270`. `--frame_bearing` turns it into the frame.
- **The anemometer height** is in `--units`; omitted, it is 10 m, the WMO standard exposure.

The wind follows the neutral log law through the observation, over `--roughness_length` (default
0.03 m, open flat terrain such as an airfield; snow is about 0.005 m, low crops 0.1 m, farmland
with hedges 0.25 m, suburbs 1 m), up to `--surface_layer_depth` (default 200 m), and is held
above, in one direction at every height. The surface air is the run's own (vehicle file,
`--temperature`, `--humidity`, `--pressure`, which may be combined with `--wind`). Above the
ground the air lapses dry-adiabatically, the pressure is hydrostatic and the water vapor is held,
so the relative humidity rises with height, capped at 100 percent with a warning
([THEORY §7](THEORY.md#7-the-stratified-moving-atmosphere)). The profile is sampled at
`--profile_levels` per factor of ten (default 16) and passes through the observation exactly; no
level above the ground is calm.

`--write_atmosphere FILE` writes the profile each trajectory used in the `--atmosphere` format
(z in `--units`, `T_K`, `RH`, `p_kPa`, `u_mps`, `v_mps`), which fed back reproduces the run. A
calm observation (`0kt`) gives a profile with the temperature lapse alone. Completely dry air
(`--humidity 0`) does not replay, since `--atmosphere` takes a humidity never above 1 percent for
a fraction, and the run warns. On the 2017 windy day, `--wind` from the ground stations'
anemometers recovered all of the improvement over straight rays that the measured LIDAR and
balloon profile gave, and on the calm day it did no harm ([FIDELITY](FIDELITY.md#wind-profile-and-stability)).

**Stability classes are unvalidated.** `--stability` takes the TA Luft classes `I` to `V`
(Pasquill `F` to `A`), each a power-law wind and a temperature gradient in the surface layer
([THEORY §7](THEORY.md#7-the-stratified-moving-atmosphere)). None has been shown to help, one
day scored after the fact cannot choose a class, and a run that takes one warns. Use `neutral`
unless the class is known from measurements.

### Ground effect

The ground effect is the interference between the direct ray and the one reflected off the
ground: Chien–Soroka's spherical-wave reflection coefficient with Chessell's band averaging, as
panam computes it ([THEORY §5](THEORY.md#5-ground-effect)). It is **on by default** wherever the
spectral path runs; `--no-ground_effect` turns it off.

**It is a gain, not a loss.** Over rigid ground with the receiver on the surface it tends to
+6 dB, the pressure doubling of a perfect reflector, and it goes sharply negative only in the
interference dip. On the S-76D demo track with the B407 spectral database, concrete reads several
dB louder than grass, which roughly cancels its own reflection gain against the dip
([FIDELITY](FIDELITY.md#ground-effect)). `--flow_resistance` spans five orders of magnitude and
the answer is sensitive to it, so name the surface rather than guessing.

- **Impedance.** Delany–Bazley by default; `--impedance_model` selects Miki, variable porosity
  (`--porosity_rate`) or a hard-backed layer (`--layer_depth`). Delany–Bazley is poor for
  grassland: the 2017 Amedee pole microphones fit variable porosity with σ<sub>e</sub> ≈ 200.
- **Turbulence** decorrelates the two paths and fills the nulls; at zero the dips are
  unrealistically deep. Chessell's form takes `--turbulence`, HARMONOISE's `--gamma_t`.
  `--roughness` lowers the coherent reflection.
- **Tables.** The ground gain is interpolated from a per-trajectory table, within about 0.001 dB
  of the exact model; `--exact_ground` evaluates it directly, several times slower.

### Receivers

Each observer's level is evaluated at its microphone, the observer's point raised by the
microphone's height. What the microphone is decides the ground model:

| Microphone | Ground model | How to ask for it |
| --- | --- | --- |
| flush with the ground, height 0 | the plate | `--plate_table` |
| on a pole, height h | the point receiver at h | the default on the spectral path, with h from `--receiver_height` or the observer file |
| anything, with `--no-ground_effect` or on the broadband path | none: the free field at the microphone | |

With a four-column observer file every observer takes the run's receiver. A fifth column gives
each microphone's own height, so one run can hold flush microphones with a few on poles: under
`--plate_table`, height 0 is a plate and anything else a point receiver. The log reports how many
observers took each.

### Ground-plane microphones

A flush ground-plane microphone sits on a plate: rigid close in, soft ground beyond, which the
point receiver cannot model. `--plate_table` models the plate from panam's axisymmetric
boundary-element solution for a GRAS GR1425 plate (400 mm, 8 mm thick, on the ground) with the
microphone at 3/4 of the radius, per SAE ARP 4055 ([THEORY §6.2](THEORY.md#62-ground-plane-microphones)).
A table is built for one ground, and the run refuses any other:

| Table in `plates/` | Ground |
| --- | --- |
| `gr1425_db225.nc` | `--flow_resistance 225` (Delany–Bazley grass, the default) |
| `gr1425_vp200.nc` | `--impedance_model variable_porosity --flow_resistance 200` (the fitted Amedee ground) |

Build another with panam's `axisymmetric_bem.table(bands, …)` and `write_netcdf` (about 20 s on
14 cores at 31 bands). Turbulence and roughness do not apply to the plate, whose reflection is
millimeters from its direct path. `--plate_mic_azimuth` says which way the microphone is offset
(default outboard, away from the track). The plate's reflection coefficients are interpolated
from per-trajectory tables, within 0.0005 dB rms in SEL ([FIDELITY](FIDELITY.md#ground-and-plate-tables));
`--exact_ground` evaluates them directly.

### Terrain

Without `--terrain` the ground is the plane z = 0 everywhere. `--terrain FILE` gives it a shape:
grid observers stand on the terrain (with the microphone at `--receiver_height` or on its plate
above it), and each source–receiver pair's ground term is taken over its **mean ground plane**,
the least-squares line through the terrain heights along the path's ground projection. That is
the mean plane of CNOSSOS-EU, which NORAH2 and ECAC Doc 32 use over terrain
([THEORY §6.3](THEORY.md#63-terrain)). Spreading, absorption and the sphere lookup still take the
real three-dimensional path. A level terrain file reproduces the run without it, bit for bit.

**Datums.** The local z = 0 is `--ground_elevation` on the trajectory's vertical datum, so

    z_ground(x, y) = H(latitude, longitude) + --terrain_offset − --ground_elevation

with H the grid's height, interpolated bilinearly. NICE-OPS does not convert between datums:
`--terrain_offset` states the difference, and the log repeats it. At the 2017 Amedee site the
3DEP grid is NAVD88 while the flight altitudes and surveyed microphone heights are ellipsoidal,
23.6 m (77.43 ft) apart, and the AC tracks put their z = 0 at 3926.74 ft ellipsoidal:

```sh
niceops -t 286264AC.csv -b 286264_observers.csv ... \
    --origin <ref_lat> <ref_lon> --frame_bearing <bearing of the track's +x> \
    --ground_elevation 3926.74 --terrain amedee_3dep.nc --terrain_offset -77.43
```

**What the run checks**, since a datum or frame mistake gives a plausible-looking wrong answer.
The grid must cover every observer and every state's footprint, with a margin and no missing
node (the refusal names the extent needed), in meters. An AC track's own latitudes and longitudes
must lie within 30 m rms of where the run's frame puts its x/y (a warning above 3 m), and where
its alt − z is one number, `--ground_elevation` must match it to 1 m. A state more than 10 m below
the terrain is refused as a datum mismatch; up to 10 m warns. Scattered observers keep their
file's z, with a warning where it differs from the terrain by more than 1 m;
`--observers_on_terrain` stands them on the terrain instead.

**What is not modeled.** Terrain does not screen or diffract: a ridge between source and receiver
puts the source below the mean plane, and its height is set to zero, so the level there is too
high. The run warns when it happens, with how many observers had such paths (a receiver in a dip
raises the same flag). Each path's ground is one plane. `--refraction` and `--propagation_table` are refused, because they are built over flat
ground, and on the broadband path `--terrain` only places the observers.

**Cost and effect.** Plates keep their tables over terrain; point receivers over real terrain are
evaluated exactly, about 4 times the CPU, which a test array of tens of microphones does not
notice. At Amedee terrain moves flush-plate SEL by hundredths of a dB; at a sloped site it can
move it by a dB at low elevation ([FIDELITY](FIDELITY.md#terrain-at-amedee),
[sensitivity](FIDELITY.md#sensitivity-studies)).

### Refraction

`--refraction linear` bends the rays through the `--atmosphere` or `--wind` profile. For each path
it fits a straight line to the effective sound speed (the sound speed plus the wind along the
path) over the heights the path spans, the equivalent linear profile of Harmonoise and Nord2000.
Every ray is then a circular arc with closed-form eigenrays, and the path gives the sphere's
direction, the spreading, absorption along the arc, the ground geometry, the reception time and
shadow zones ([THEORY §8.1](THEORY.md#81-linear-sound-speed-profile-circular-rays)).

Under upward refraction, receivers past the ray that grazes the ground are in a shadow. The level
decays along the least attenuated creeping wave and is held at most `--shadow_floor` (default
20 dB) below the boundary level, past which turbulent scattering is taken to dominate;
`--scattering` puts the scattered energy in directly. With the ground effect on, the level near
the boundary comes from a tabulated exact field, since ray theory fails within a couple of
diffraction lengths of it.

What to know before relying on it:

- **Wind directions are compass bearings**, turned by `--frame_bearing` (default +x east). The
  shipped AC tracks need `--frame_bearing 210`, or their wind is read 120° off; the log says so.
- **One gradient per path.** A log-law profile's steep lowest meters are averaged in, and
  equivalent-linear models have been found up to 10 dB off a parabolic equation near the ground
  at 500–2000 Hz. `stratified` removes this approximation.
- **The effective-sound-speed approximation** overstates the wind by 1/cos(elevation), worst on
  steep paths, where refraction matters least.
- **One direct and one reflected ray**: the extra bounced arrivals of strong downward refraction
  are left to `stratified`.
- **The source feels the wind only through `--speed_reference`**: refraction bends the path, not
  the operating condition.
- **Uniform air gives the straight result** to 1e-9, and the peak search and `--db_cutoff` stay
  exact.

How the models compare with the parabolic equation and the 2017 flights is in
[FIDELITY](FIDELITY.md#refraction-models-against-the-pe).

### Refraction through the profile's layers

`--refraction stratified` drops the single fit. The effective sound speed is linear between the
profile's own levels, so a ray is an arc in each layer, and eigenrays are solved per pair over
several families: launched down or up, arriving directly or after dipping below the microphone,
the ground-reflected ray and, in a surface duct, rays that bounce before their last descent
([THEORY §8.2](THEORY.md#82-stratified-profile-ray-tracing-through-layers)). Every arrival is
kept and summed in energy, each through its own direction on the sphere. Use it whenever the
profile has structure near the ground (a log-law wind above all), and downwind, where it is the
model that agrees with the parabolic equation.

- **Ducts.** Under a morning inversion, past every direct ray, the strongest bounced ray becomes
  the primary arrival.
- **Give the profile a level at the ground.** The speed is held below the lowest level, so a
  profile starting at a LIDAR's first gate (46 ft) has a 14 m layer of constant air, which the
  rays follow for kilometers. Calm, isothermal air at the ground is resolved only as finely as the
  run's highest band can, so it does not throw the shadow boundary out.
- **Shadows.** With the ground effect on, each band's boundary lies between the layered rays' and
  the linear fit's, and from 0.9 of the layered boundary on the level is the linear fit's ray.
  Without it the layered boundary and the creeping-wave loss stand.
- **Limits.** A receiver at or above the source takes the linear model's ray. Focusing at a
  caustic is capped at +10 dB per ray. Two close arrivals at a caustic can be missed, and an
  elevated duct can leave a micro-shadow the wave does not have ([below](#micro-shadows-and-skip-zones)).
- **Cost.** On the demo grid with a point receiver, about five times `linear`'s
  acoustic loop ([Relative cost](#relative-cost)), and more in a duct, where the arrivals
  multiply.

### Low sources near grazing: the resolved shadow

`--shadow_model resolved` is an opt-in change to how `--refraction stratified` places and fills a
ground-grazing shadow, for low sources under light wind (landings, say), where the default
`legacy` shadow is several dB deeper than a parabolic equation's. Each band sees only the profile
it resolves: below the height where it gathers `--shadow_phase` radians of vertical phase
(default 1), the profile is replaced by its tangent, and the band's boundary comes from that. A
band whose boundary is within two diffraction lengths of the source takes the residue series at
the source's height (`--cramped_ground series`, the default), or the ground model along the
straight path (`straight`) or the ray's geometry (`ray`) ([THEORY §8](THEORY.md#8-refraction)).

It improves low R44 approaches markedly but over-reaches in strong wind with a 1.2 m receiver,
failing the windy B407 day's do-no-harm test, which is why it is opt-in
([FIDELITY](FIDELITY.md#resolved-shadow-against-the-pe)). Treat other `--shadow_phase` values as
sensitivities. Scattered observers are looked up by their height above the ground. Its series
table costs 2–3 s of CPU per band per receiver height (5–8 s for 24 bands on 14 threads).

### Micro-shadows and skip zones

Under an elevated duct with the fastest air at the ground (a log-law headwind below an inversion),
sound turned back aloft turns up again a few feet above the ground and reaches neither a flush nor
a 1.2 m microphone, and short of where it first comes down no ray arrives. The wave fills both from
the caustic. Two opt-in corrections exist
([THEORY §8.2](THEORY.md#82-stratified-profile-ray-tracing-through-layers)): `--tunneling` adds
the rays that turn within reach of the microphone, with the Airy tail below their turning point
(about +20% cost), and `--skip_fill` adds a caustic's Airy tail on the dark side of where rays from
aloft land (2–4 times). They bring such deep frames closer to the parabolic equation and the
measurements, but leave the tails above 800 Hz 10–20 dB low, can over-predict other far frames,
and do not move the SEL ([FIDELITY](FIDELITY.md#micro-shadows-and-skip-zones)). A table built from
the parabolic equation for the day's profile is the alternative.

### Propagation loss tables

`--propagation_table` takes what the path does to each band from a table computed elsewhere, so
that a full-wave method such as the parabolic equation
([`niceops_propagation_table`](#niceops_propagation_table)) can carry a footprint. The table
holds each band's level relative to free field at the straight-line distance: the ground,
refraction, shadows and diffraction. The run takes the sphere's level in the straight-line
direction, less spherical spreading and absorption, plus the table
([THEORY §9.3](THEORY.md#93-propagation-loss-tables)). Absorption stays the run's unless the
table declares that it includes it. A table assumes one stratified atmosphere over flat, uniform
ground, so its axes are `bearing` (compass bearing the sound travels toward, uniform over 360°;
one value means every direction), `source_height` and `receiver_height` (above the ground, m),
`range` (horizontal, m) and `band` (Hz; every run band within `--band_limits` must be there).

- **Interpolation** is linear in dB along each axis.
- **Coverage.** A pair the table does not cover, by a NaN among the values it would interpolate
  or by falling outside an axis, takes the run's own model: `--refraction` if on, the ground
  effect if on, free field otherwise.
- **Receivers.** A point microphone whose height is outside the table's is refused. A flush plate
  microphone reads the lowest receiver height plus the plate's level over bare ground there, at
  the pair's azimuth and elevation. For plates, give the table a receiver row at 0: over a raised
  row the bare-ground term is a comb in elevation, which the plate's offset has to follow.
- **Bearings are compass bearings**, so a rotated frame needs `--frame_bearing`.
- **With `--refraction` too**, every state is evaluated rather than windowed from the peak.
- **The file** is netCDF, SI: `relative_level(bearing, source_height, receiver_height, range,
  band)`, a coordinate variable per axis, and the attributes `propagation_table_version` (1),
  `method`, `includes_absorption` and `provenance_*`. Anything that fills it can be read
  ([internals/pe.md](internals/pe.md) has the specification).

Interpolating between source heights is where a table's error is: build tables a few wavelengths
apart at the highest band, or densely at the lowest heights, and check them with
`--validation_csv`. Reading a table costs little.

## Metrics and output files

### Output files

Each metric file (`-D`, `--lden_file`, `--cnel_file`, `--lamax_file`, `--nax_file`,
`--ta_file`, `--csel_file`, `--coverage_file`) is laid out as `-o` is, and is described with its
metric below. `--coverage_file` writes, per observer, the energy-weighted share of its level that
came from sphere cells the database did not measure. With spectra, each cell counts by the share
of its A-weighted band energy at the sphere that was measured, so a band the builder's SNR gate
dropped counts by the energy it carries; coverage recorded per channel counts a cell as
measured or not. The
others:

- `-o`: the SEL per observer. A grid writes the x and y axes on the first two lines, then one row
  of levels per y; scattered observers give one row each.
- `-m`: appends a line for the case: its name, the energy mean of the per-observer SELs (10 lg of
  the mean of 10^(L/10), not a mean of decibels), their maximum, the fraction of observers at or
  above 65 dB, and the area at or above 65 dB in the run's units squared, contoured as
  `--contour_file` contours it (NaN for scattered observers).
- `-i`: per-observer level histories, one file per observer.
- `--band_history`: per-observer band histories, one file per observer with a row per state:
  `emission_time`, `reception_time` (s), `range` (the ray tube's under `--refraction`),
  `elevation` (of the emission direction below horizontal, deg), the sphere's `phi` and `theta`,
  `level_dBA` (what `-i` writes), then `band_<Hz>`, each band's unweighted level at the microphone
  (dB re 20 µPa) after spreading, absorption, ground or plate and any shadow loss. Bands outside
  `--band_limits` and culled states are empty; the `-a` ambient goes into `level_dBA` only.
- `-S`: the derived state (time, position, velocity, load factor, C<sub>T</sub>, μ, path
  angles), metric whatever `--units` says.
- `--provenance_file`: one row per observer, grid or scattered: `x`, `y`, `exposure_db` (the SEL
  for one operation, the energy sum of the SELs for several), then the shares of that exposure
  from states whose advance ratio, path angle or C<sub>T</sub> lay outside the database and were
  clamped to its edge, and from states whose load factor was held at `--max_load_factor` or
  `--min_load_factor`; then the sphere that supplied the most of it (its index in the database),
  its share, and its stored μ, γ and C<sub>T</sub>. Each state counts as the SEL counts it (inside
  the 10 dB-down window, or the whole history under `--full_history`), and its energy is shared
  among its spheres by their blend weights, so a sphere with nothing measured toward the observer
  is still credited. The run evaluates every state, as `--coverage_file` does, whose share is
  weighted the same way. An observer that heard nothing is written as NaN.
- `--provenance_spheres`: netCDF, the whole sphere breakdown behind `--provenance_file`'s top
  sphere. Per observer, its `x` and `y` (meters) and `exposure_db`; per sphere, its stored
  `advance_ratio`, `flight_path_angle` (deg) and `thrust_coefficient`, in the database's group
  order (the index `top_sphere` gives); and `energy_share(observer, sphere)`, each sphere's
  share of the observer's exposure, summing to one (NaN for an observer that heard nothing).
  Summing the shares over a group of spheres gives that group's share, for any grouping: by
  speed, by flight day, measured or synthetic. Mostly zeros, it is compressed.
- `--contour_file`: GeoJSON (WGS84 longitude and latitude, RFC 7946), one MultiPolygon feature
  per level from `--cmin` to `--cmax` in steps of `--cstep`, covering where `--contour_metric` is
  at or above it, holes included, with the level and its area (m², km², square miles). It works
  in a build without plotting and opens in any GIS. Contours are drawn by linear interpolation
  between grid nodes, so a contour's position is as good as the grid's spacing; the area of a
  smooth footprint converges with the square of it. A culled observer counts as below every
  level, and a region reaching the grid's edge is closed along it, so widen `--grid` until the
  outermost level closes inside it. It needs the run's frame (`--origin`, or a geodetic track)
  and the observer grid.
- `--raster_file`: a CF-1.8 netCDF of the grid, with `x` and `y` in meters in the run's frame
  and, when the run has a frame, each node's WGS84 `latitude` and `longitude`, which GDAL, QGIS
  and xarray read as the grid's position. It holds one variable for each gridded level the run
  writes: `sel` (one trajectory only), `dnl`, `lden`, `cnel`, `lday`, `levening`, `lnight`,
  `laeq`, `lamax`, `nax`, `time_above` and `unmeasured_share` (`--coverage_file`), each present
  only if its own file is asked for, with its values. An observer that heard nothing is the
  fill value, NaN, where the text files write `-inf`.
- `--latlon`: a `_latlon` companion beside every level file (`out.csv` → `out_latlon.csv`), one
  row per observer of latitude, longitude and level, always scattered.
- `--epnl_file`, `--band_frames` ([EPNL](#epnl-and-reception-time-frames)); `--write_atmosphere`
  ([A profile from the surface wind](#a-profile-from-the-surface-wind)); `--powerfile`;
  `--export_sphere` and `--export_aam` (a whole sphere, as a one-group `.nod` or in AAM/RNM
  format).

### Run records and option files

`--manifest FILE` writes one JSON file saying what produced the run's outputs:

- the build: its version (`git describe`, `-dirty` when the tree had uncommitted changes, or
  `unknown` when built outside git), compiler, build type, and whether asserts and plotting were
  on (`--version` prints the version alone);
- the command line, and the options given on it or in a `--config` file (the rest are the
  build's defaults);
- every input (trajectories, databases, vehicle, observer, atmosphere, terrain and table files)
  and every output, figures included, each with its size and SHA-256, which `shasum -a 256`
  reproduces; the log file is named only, since runs append to it;
- each database's root attributes and scalar variables, which record how it was built;
- the number of warnings and errors, the first 20 of them, the exit status and the run time.

It is written after the run, including one that stops on a refused input, so the exit status it
records is the one the process returns. A command line that does not parse (an unknown option,
an unreadable `--config`) stops before there is a run to record, and writes none. A run writing more than 1000 files (per-observer histories) lists the
first 1000 and counts the rest. Hashing reads every input once more: about 0.3 s for the largest shipped
database (72 MB). A manifest that names one of the run's inputs or its log file is refused before the run; one
that names an output is refused after it, keeping the output.

`--config FILE` reads options from a TOML file, one `name = value` per line, using the long
option name without its dashes. Flags take `true`; a negated flag is written as its own name.
A key that names no option is refused. The command line overrides the file. The demo case of [Quick start](#quick-start) as a file:

```toml
no-spectral = true
no-ground_effect = true
```

### Exposure (SEL) and the 10 dB-down window

The level history is sampled at the trajectory's emission times, but a microphone accumulates
exposure in reception time: each sample lasts T(1 + ṙ/c), compressed while the source closes and
stretched while it recedes, with T the trajectory's own sample interval (not `-s`)
([THEORY §10.2](THEORY.md#102-sound-exposure-level)). The sum runs over the 10 dB-down window,
from the first to the last sample within 10 dB of the peak, dips between them included: the
duration the aircraft noise metrics use (14 CFR 36 A36.4.5). `--full_history` sums the whole
history instead. Use the window to match a measured SEL's convention, and `--full_history` to sit
beside AAM or NORAH2, which integrate everything. `-a` adds the ambient to every sample before
the window is found.

### Several trajectories, and the day-night level

`-t` may be repeated. `-d`, `-v`, `-o`, `-S`, `--csel_file`, `--epnl_file` and `--band_frames`
then take one value shared by every trajectory, or one per trajectory, and so does each
trajectory's air: `--atmosphere`, `--wind`, `--humidity`, `--temperature`, `--pressure` and
`--start_utc`, so a day's operations can each fly in their own weather and at their own time.
The profile's shape (`--roughness_length`, `--stability`, `--surface_layer_depth`,
`--profile_levels`) is the run's. One output name for
several tracks becomes one per track (`<track>_state.csv`, `<track>_LCE.csv`, `<track>_EPNL.csv`,
`<track>_band_frames.nc`, `<track>_atmosphere.csv`), and `-i` and `--band_history` files are
named for their track. Plotting is skipped.

```sh
./build/niceops -t morning.csv -t evening.csv -t night.csv -D integrated/dnl.csv
```

DNL is a property of a set of operations: each trajectory's exposure gains a penalty on samples
heard in the evening or at night, the results are summed as energy, and the sum is normalized to a
24-hour day, so one daytime operation gives its SEL less 49.37 dB
([THEORY §10.4](THEORY.md#104-daynight-metrics)). Each penalty applies to every sample emitted in
its period, local time (the propagation delay is not added):

| Metric | Option | Evening | Night | Definition |
| --- | --- | --- | --- | --- |
| DNL | `-D` | — | 22:00–07:00, +10 dB | 14 CFR Part 150 |
| Lden | `--lden_file` | 19:00–23:00, +5 dB | 23:00–07:00, +10 dB | EU Directive 2002/49/EC, Annex I |
| CNEL | `--cnel_file` | 19:00–22:00, +10 lg 3 = 4.77 dB | 22:00–07:00, +10 dB | California Code of Regulations, Title 21, §5001 |

A period's opening hour belongs to it (07:00 is day, 19:00 evening); a member state's shifted
evening is not offered.

**The period levels** (`--lday_file`, `--levening_file`, `--lnight_file`, and `--laeq_file` over
`--laeq_window`) are each the energy heard in one period over that period's length, with no
penalty: Lnight, for example, is the EU's night indicator and the basis of the WHO's 2018 night
guideline. Samples outside the period carry no energy, and an operation's 10 dB-down window is
found among those inside. An observer that heard nothing in the period is written as `-inf`.
Lden is its three periods recombined as the directive writes it, exactly for operations that
lie wholly in one period, and under `--full_history` for any. The clock is the track's `utcsec` column (column 0 of an AC track) or
`--start_utc`, shifted by `--utc_offset`, and a track without one is refused. By default each operation contributes its own
penalized 10 dB-down window; `--full_history` integrates its whole history, as AAM and NORAH2 do,
which reads a few tenths of a dB higher ([FIDELITY](FIDELITY.md#day-night-window-choice)). Use it
for a DNL meant to sit beside their contours.

### Maximum level, time above and C-weighted SEL

**LAmax (`--lamax_file`)** is a sound level meter's largest reading through a time weighting,
which `--lamax_weighting` picks ([THEORY §10.3](THEORY.md#103-maximum-level-time-above-and-number-above)):

| Value | Maximum | Coarsest trajectory step |
| --- | --- | --- |
| `fast` (default) | LAFmax, IEC 61672-1 F (τ = 0.125 s) | 0.0625 s |
| `slow` | LASmax, IEC 61672-1 S (τ = 1 s) | 0.5 s |
| `leq0.5` | the largest running 0.5 s Leq: ECAC Doc 32's LASmax surrogate, and NORAH2's | 0.25 s |
| `instantaneous` | the largest sample, with no time weighting | any |

The weighting is evaluated exactly in reception time, with nothing resampled. Over several
operations LAmax is the largest of theirs. A weighting is refused on a track whose mean sample
interval exceeds the table's last column, where the "maximum" would be the largest sample under
another name. The 2017 AC tracks (0.02 s) pass; WGS84 waypoints need `-s 0.0625` or finer under
`fast`; a coarse AC track needs a finer file.

**NAx (`--nax_file`)** counts the operations whose own weighted maximum reaches `--nax_threshold`.
**Time above (`--ta_file`)** is the seconds of reception time the weighted level exceeds
`--ta_threshold` (strictly), summed over operations; a crossing within a sample is credited when
it happens, and under `instantaneous` each sample counts whole. Neither threshold has a default.

**C-weighted SEL (`--csel_file`)** is the exposure integral over the C-weighted level, from the
same band levels at the microphone and the same `--band_limits` as the A-weighted one. It finds
its own 10 dB-down window on the C-weighted history (or uses the whole history under
`--full_history`). C-weighted energy sits mostly at low
frequency, so **the band limits matter**: the default 20 Hz limit leaves out a database's 10–16 Hz
bands. Refused on the broadband path, with `--db_cutoff` and with `-a`.

### EPNL and reception-time frames

`--epnl_file FILE` writes each observer's EPNL, the certification metric of 14 CFR 36 Appendix A,
A36.4 (identical to ICAO Annex 16 Vol. I, Appendix 2); `--band_frames FILE` writes the
one-third-octave levels it is computed from. Both need the spectral path, are refused with `-a`,
and give each trajectory its own file. The procedure follows panam's term for term
([THEORY §10.5](THEORY.md#105-effective-perceived-noise-level)).

- **Records.** Each state's band energy is spread over the stretch of reception time it arrives
  in and summed into 0.5 s frames, so the frames hold the history's whole exposure. Frame *k*
  covers [0.5 *k*, 0.5 (*k* + 1)) s on the trajectory's clock, shared by every observer.
- **Bands.** All 24 bands from 50 Hz to 10 kHz are needed, each matched to the database band
  within 5 percent of its nominal center. A database lacking one, or `--band_limits` leaving one
  out, is refused with the bands named. The shipped spectral databases carry 10 Hz–10 kHz.
- **Levels** are those at the microphone the run models, with no frequency shift: a flush plate
  carries its pressure doubling, as the 2017 measurements do, where 14 CFR 36 prescribes the
  1.2 m microphone of the default point receiver.
- **Masked bands** (`--epnl_masking`, the default). A predicted band falls as far as absorption
  takes it, and read literally the plunge reads as tones. A band below its noy threshold
  therefore enters the tone correction at that threshold and carries no tone, as tones in masked
  bands are excluded from a measured spectrum. `--no-epnl_masking` is the procedure as written.
- **Not included**: the corrections A36 makes to a *measured* spectrum (instrumentation,
  background, reference conditions, pseudo-tones from ground reflections). Every predicted
  audible tone counts, ground-effect ripple included.

**`--epnl_file`** has one row per observer: `observer` (the `-i` index), `x`, `y`, `z` in
`--units`, `EPNL` (EPNdB), `PNLTM` (TPNdB), `duration_correction` (D; EPNL = PNLTM + D),
`tone_correction` (PNLTM − PNL at its record), `tone_band` (Hz, `nan` if none),
`bandshare_adjusted` (1 when the band-sharing adjustment applied), `PNLTM_time`,
`duration_first`, `duration_last` (record centers, s, reception time) and `clipped` (1 when the
duration reaches the first or last record heard, so the event may not be whole). An observer
with no finite PNLT is written as `-inf,-inf,nan,...`.

**`--band_frames`** is netCDF-4: `level(observer, frame, band)`, float, dB re 20 µPa, unweighted;
`time(frame)`, the frame centers; `band_center(band)`, Hz, within `--band_limits`; the
microphones' `x`, `y`, `z` and `microphone_height` in `--units`; and for a grid `grid_x` and
`grid_y` (observer = *iy* × size(`grid_x`) + *ix*). Levels are NaN outside the frames an observer
received and −inf in a band that received nothing. Global attributes record the trajectory, the
database, `frame_length_seconds` (0.5) and `band_frames_version` (1). On the demo grid it is about
318 MB.

**Cost.** Both evaluate every state. On the demo grid with B407, `--epnl_file` adds about a third
to the run; `--band_frames` more than doubles it and holds about 0.9 GB.

## Companion tools

The build also makes three tools for checking the propagation models and building tables, and
one for building tracks. They share the `--units` convention (feet by default) and the profile
format.

### `niceops_procedure`

Flies a procedure into a track: a start and a list of segments, integrated into a WGS84 track with
velocity and heading columns, which `niceops` reads as it is. Use it where you would otherwise
write sparse waypoints, whose load factor is an artifact of the interpolant
([FIDELITY](FIDELITY.md#trajectory-smoothing)); a flown procedure carries its own.

```sh
./build/niceops_procedure approach.json -o tracks/approach.csv
./build/niceops -t tracks/approach.csv --origin 39.8291 -74.41 --frame_bearing 210 ...
```

The procedure is JSON. Every length and speed names its unit in its key (`_ft` or `_m`; `_kt`,
`_mps` or `_fps`; rates `_kt_per_s` or `_mps2`; angles `_deg`; times `_s`), and a bare number is
refused:

```json
{"start": {"lat": 39.8291, "lon": -74.41, "alt_ft": 1000, "speed_kt": 100, "course_deg": 210},
 "segments": [
   {"type": "level", "duration_s": 20},
   {"type": "descent", "angle_deg": 6, "until_alt_ft": 700},
   {"type": "speed", "to_speed_kt": 40, "rate_kt_per_s": 1.5, "angle_deg": -3},
   {"type": "turn", "bank_deg": -15, "turn_deg": 30},
   {"type": "speed", "to_speed_kt": 0, "rate_kt_per_s": 1},
   {"type": "hover", "duration_s": 15}]}
```

| Segment | Keys |
| --- | --- |
| `level` | `duration_s` or `distance` |
| `climb`, `descent` | `angle` (a descent's is taken downward) and `duration_s` or `until_alt` |
| `speed` | `to_speed`, optional `rate` (default 0.5 m/s²) and `angle` (default level) |
| `turn` | `bank` (signed, positive right, under 60°) and `turn`, the course change |
| `hover` | `duration_s`; it must start at rest |

- **The start** is a position, `alt` on the track's datum (what `niceops` takes
  `--ground_elevation` off), a ground speed and a compass course. An optional `ground_alt` (default
  0) is the ground on that datum: a segment that goes below it is refused.
- **Transitions are flown**: the path angle changes at a limited pull (`limits.pull_g`, default
  0.1 g), the bank at a limited roll rate (`limits.roll_rate_deg_per_s`, default 10), and speed at
  the segment's rate, so a join puts no spike in the load factor. A climb overshoots its altitude
  by the leveling-off. A turn rolls in, holds and rolls out within itself, ending level-winged on
  the asked course to within a degree; one too short to reach full bank banks only as far as it
  needs. A target altitude already passed, and any key the format does not know, are refused.
- **The end** is flown on to the next written step, so the last state is written.
- **Calm air**: the heading is the course and the ground speed the airspeed. `step_s` (default
  0.2) is the spacing of the written states.
- Whether the flown conditions stay inside a database's envelope is what `--provenance_file`
  reports.

### `niceops_propagation_table`

Builds a table for `--propagation_table` ([Propagation loss tables](#propagation-loss-tables)) by
one of two methods. `--method pe` (the default) marches the wide-angle parabolic equation through
the profile's effective sound speed along each bearing, band-averaged over `--samples`
frequencies (default 5); levels steeper than `--max_angle` (default 50°, at most 60°) are NaN, and
those pairs take the run's own model. `--method ground_model` tabulates the straight-ray ground
model of `--ground_effect` on one bearing: a footprint through it matches `--ground_effect
--exact_ground` up to interpolation, which makes it the check on a table's plumbing.

```sh
niceops_propagation_table --method pe --atmosphere atmospheres/example_windy.csv \
    --database databases/B407_fixed_load_factor.nod --band_limits 20 2000 \
    --source_heights 40 90 140 190 240 --receiver_heights 3.937 \
    --max_range 15000 --step 10 -o tables/windy.nc
```

| Option | Meaning |
| --- | --- |
| `--atmosphere` | The profile, as for `niceops`. Required for `pe`; gives the sound speed for `ground_model`. |
| `--sound_speed` | Sound speed at the ground, m/s, for the impedance and the ground model. Default the profile's at the ground, else 343.2. Match the run's, which `niceops` logs. |
| `--database`, `--bands`, `--band_limits` | Band centers from a database, or given, and the range kept. |
| `--source_heights`, `--receiver_heights`, `--max_range`, `--step` | The table's heights above the ground, farthest range and range spacing, in `--units`. Required. |
| `--bearings` | Bearings evenly over 360°. Default 36, or 1 when the profile has no wind. |
| `--flow_resistance`, `--impedance_model`, `--porosity_rate`, `--layer_depth` | The ground, as for `niceops`, except that `--layer_depth` is in meters whatever `--units` says. |
| `--ppw` | `pe`: grid points per wavelength in height. Default the PE's own. |
| `--march` | `pe`: `receiver` (default) marches once per receiver height and reads every source height by reciprocity, 2–3 times faster; `source` marches once per source height. They agree to about 0.01 dB. |
| `--gamma_t`, `--outer_scale`, `--realizations`, `--seed`, `--gamma_profile` | `pe`: turbulence by phase screens, for reference only. Defaults still air, 10 m, 1 realization, seed 1. `--gamma_profile` is a CSV of z (in `--units`) and γ. |
| `--validation_csv`, `--validation_samples`, `--validation_max_error` | Check the table between its nodes (below). Defaults 21 samples and 1 dB. |
| `-j`, `--threads` | Worker threads. Default the environment's `NICEOPS_THREADS`, else one per hardware thread. The table is the same bits at any count. A `NICEOPS_THREADS` that is not a positive whole number is refused. |
| `--units`, `-o` | `ft` (default) or `m`; the table to write (SI). |

A PE table through the ISA, 21 bands to 2 kHz, 23 source heights and ranges to 15,000 ft every
10 ft took 103 s on 14 cores.

**Checking a table between its nodes.** `--validation_csv FILE` solves the method again at the
midpoint of every source-height, receiver-height, range and bearing interval (an axis with one
value keeps it): a PE table against marches at twice the height-grid density, half the range
spacing and `--validation_samples` frequencies per band, a ground-model table against the model
itself. The CSV has a row per point and band (table, reference, signed error, status, SI
coordinates); each band's bias, rms, largest error and gap count go to stderr. The tool exits
with status 2, keeping the table, if a covered point is off by more than `--validation_max_error`
dB, a gap hides a finite reference, or a band has no covered point; points past the PE's angle
limit are reported but not scored. The error includes the reference's own resolution (raise
`--validation_samples` to check convergence), a PE is not a measurement, and passing at midpoints
bounds neither the error everywhere nor the hand-over to the run's own model. Use the failing
coordinates to add height nodes, shorten the range step or add bearings. Measured table errors
are in [FIDELITY](FIDELITY.md#propagation-tables).

### `niceops_pe_compare`

Runs one path through every model and writes each band's level relative to free field: straight
rays, `linear`, `stratified`, and the wide-angle parabolic equation (PE), a full-wave solution
with the ground's impedance in the grid. Absorption is left out, and turbulence too unless
`--gamma_t` puts it in the PE. With `--shadow_model resolved` it adds a `resolved_<f>` column per
band, with `--shadow_phase` and `--cramped_ground` as in `niceops`.

```sh
niceops_pe_compare --atmosphere atmospheres/example_windy.csv --toward 270 \
    --source_height 328 --receiver_height 3.937 --bands 125 250 500 1000 --max_range 8000
```

The PE is a split-step Padé march with Collins' self-starter and a fourth-order grid
([THEORY §9.1](THEORY.md#91-the-wide-angle-parabolic-equation)). `--max_angle` (default 50°, at
most 60°) limits it: a level steeper than that from the source is NaN. `--range_step`, `--ppw`
(default 10 points per wavelength, 20 under `--scheme cn` or a second-order grid), `--starter`,
`--vertical_order` and `--boundary_order` tune the grid, and `--scheme cn` is the older
Crank–Nicolson march, good to about 30°, kept as a check (`--range_ppw` range steps per
wavelength, default 10). `--step` sets the output range spacing
(default a hundredth of `--max_range`), and `--no_pe` runs the ray models alone. Turbulence
enters as phase screens: `--gamma_t`, `--outer_scale` (default 10 m), `--seed`, and `--samples`
frequencies per band (default 5), each its own realization. `-j` defaults as for
`niceops_propagation_table`. Five bands to 1.5 km take about half a second on 14 cores.

The PE agrees with exact solutions to a few hundredths of a dB
([FIDELITY](FIDELITY.md#parabolic-equation-against-exact-solutions)), but **its deep shadow in
still air, 40–80 dB down, is not to be believed**: with no turbulence, nothing scatters energy
in. The 2017 measurements support HARMONOISE's moderate strength, γ<sub>T</sub> up to about
3e-6, for `--scattering` ([FIDELITY](FIDELITY.md#turbulent-scattering-strength)).

### `niceops_benchmark`

Runs the four cases of Attenborough et al., "Benchmark cases for outdoor sound propagation
models", JASA 97, 173 (1995), to 200 m and to 10 km, with the PE and, where one exists, the exact
solution ([THEORY §9.2](THEORY.md#92-analytic-reference-solutions)). `--cases` picks some of
the four (default all), and `--out` the output directory. It takes about ten seconds:

```sh
niceops_benchmark --out figures/benchmark
build/pyenv/bin/python tools/plot_benchmark.py figures/benchmark
```

The plot script draws the results in the paper's layout and prints how far the PE is from the
exact solution ([FIDELITY](FIDELITY.md#the-attenborough-et-al-benchmark)).

### Plot scripts

The scripts in `tools/` run under the build's Python, `build/pyenv/bin/python`:

- `tools/plot_pe_compare.py OUTDIR CSV:TITLE ...` plots `niceops_pe_compare` output: one figure
  per case, a panel per band, and a summary of where each shadow begins.
- `tools/plot_benchmark.py DIR` plots `niceops_benchmark` output beside the paper's figures.

## Performance

### How much gets evaluated

Exposure covers each observer's 10 dB-down window, a little under 60% of a shipped track. So on the
spectral path the observer loop finds the peak with a per-direction bound on the level and
evaluates only the states the bound cannot rule out, logging `Evaluated N of M (P%); the rest is
provably outside every observer's 10 dB window.` The search is exact: the output is identical to
evaluating every state. A broadband run, where an evaluation is cheap, marches every state, as does
any run with `-a`, `-i`, `--band_history`, `--epnl_file`, `--band_frames`, `-D`, `--lden_file`,
`--cnel_file`, `--lamax_file`, `--nax_file`, `--ta_file`, `--csel_file`, `--coverage_file`,
`--full_history` or `--db_cutoff` ([internals/performance.md](internals/performance.md) says why).

### `--db_cutoff`

`--db_cutoff L` skips any evaluation that provably falls more than 10 dB below `L`, bounding the
level from the loudest sphere the state blends, spherical spreading and the most the ground or a
plate can add. That is exact for any observer whose *peak level* reaches `L`. An SEL exceeds its
peak level by 1 to 12 dB, though, so a cell at or above `L` in SEL can have its window cut: on
the demo case `--db_cutoff 85` changes about a quarter of the cells at or above 85 dB, by at most
0.43 dB, and none about 3 dB above the cutoff. Read it as a threshold on peak level, with a few dB
of margin for SEL. Observers that never reach it come out as `-inf` and do not plot.

### Threads

The observer loop runs on every hardware thread by default (`-j 0`), and the output is
byte-for-byte the same at any count. Work is dealt out in chunks, so machines with performance
and efficiency cores are used fully. `niceops` does not read `NICEOPS_THREADS`; the companion
tools do.

### Relative cost

The acoustic loop on the full demo grid (S-76D track, 94.3 million observer–state pairs,
8 threads), with the B407 spectral database and `--humidity 20`:

| Run | Acoustic loop |
| --- | --- |
| spectral, `--no-ground_effect` | 2.5 s |
| spectral with the ground effect (the default) | 2.9 s; 11.6 s with `--exact_ground` |
| `--plate_table` | 7.1 s; 21.0 s with `--exact_ground` |
| `--refraction linear` (`example_windy.csv`), ground effect | 7.8 s |
| `--refraction stratified` (`example_windy.csv`), ground effect | 37 s, more in ducts |

A broadband run is several times cheaper than spectral. On top of these:
`--shadow_model resolved` spends 5–8 s on setup; `--tunneling` adds about 20% and `--skip_fill`
2–4 times to a stratified run; point receivers over real terrain take about 4 times the CPU; a
propagation table is cheap to read but takes minutes to build; `--epnl_file` and `--band_frames`
evaluate every state.

## Plotting

`-p` shows the footprint with the flight track over it, `-f FILE` saves it, and `-F FILE` saves
the trajectory-state figures. Contours default to 65–100 dB in 5 dB steps (`--cmin`, `--cmax`,
`--cstep`). Multi-trajectory runs skip plotting. Plotting embeds Python, with matplotlib in the
build's own environment (`build/pyenv`), and selects the `Agg` backend itself, so it works
headless. A build without plotting refuses `-p`, `-f` and `-F`; the build options are in the
[README](../README.md).

## Comparing with older results

Several deliberate changes have moved NICE-OPS's results, the largest being the removal of a 4x
multiplier on the database's EAA. Before comparing old and new results, read
[HISTORY.md](HISTORY.md#changes-that-affect-results).
