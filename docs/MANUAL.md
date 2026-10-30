# NICEOPS manual

This manual is the user's reference: how to run NICEOPS, every option, the input and output
formats, which model to use when, and what each one costs and refuses. It states each model in a
paragraph or two. [THEORY.md](THEORY.md) has the equations and references,
[FIDELITY.md](FIDELITY.md) the evidence of how well each model agrees with exact solutions and
flight-test measurements, and [HISTORY.md](HISTORY.md) the changes that moved results. To build
the code, start with the [README](../README.md); developers start at
[internals/README.md](internals/README.md).

## Contents

- [Overview: how a prediction is made](#overview-how-a-prediction-is-made)
- [Running NICEOPS](#running-niceops): [Quick start](#quick-start) · [Command-line
  options](#command-line-options) · [Units](#units) · [Errors and refusals](#errors-and-refusals)
- [Inputs](#inputs): [Input files](#input-files) · [Trajectories](#trajectories) · [Vehicle
  file](#vehicle-file) · [Observers](#observers)
- [Noise databases](#noise-databases): [Source coherence](#source-coherence) · [Shipped databases](#shipped-databases) · [Spectra](#spectra) · [Screening a database](#screening-a-database)
  · [Choosing spheres: blending and `--blend`](#choosing-spheres-blending-and---blend) · [The
  load-factor law below 1 g](#the-load-factor-law-below-1-g) · [Speed reference](#speed-reference) ·
  [The sphere's nose](#the-spheres-nose) · [Reading the sphere's
  angles](#reading-the-spheres-angles)
- [Propagation](#propagation): [Broadband and spectral paths](#broadband-and-spectral-paths) ·
  [Stratified atmosphere](#stratified-atmosphere) · [A profile from the surface
  wind](#a-profile-from-the-surface-wind) · [Ground effect](#ground-effect) ·
  [Mixed ground](#mixed-ground) · [Receivers](#receivers) · [Ground-plane microphones](#ground-plane-microphones) ·
  [Terrain](#terrain) · [Refraction](#refraction) · [Refraction through the profile's
  layers](#refraction-through-the-profiles-layers) · [Low sources near grazing: the resolved
  shadow](#low-sources-near-grazing-the-resolved-shadow) · [Micro-shadows and skip
  zones](#micro-shadows-and-skip-zones) · [Propagation loss tables](#propagation-loss-tables) ·
  [Urban receivers](#urban-receivers) · [Buildings around the
  microphones](#buildings-around-the-microphones)
- [Metrics and output files](#metrics-and-output-files): [Output files](#output-files) · [Run records and option files](#run-records-and-option-files) · [Exposure
  (SEL) and the 10 dB-down window](#exposure-sel-and-the-10-db-down-window) · [Several trajectories,
  and the day-night level](#several-trajectories-and-the-day-night-level) · [An operations
  file](#an-operations-file) · [Maximum level, time
  above and C-weighted SEL](#maximum-level-time-above-and-c-weighted-sel) · [EPNL and reception-time
  frames](#epnl-and-reception-time-frames)
- [Companion tools](#companion-tools): [`niceops_procedure`](#niceops_procedure) · [`niceops_propagation_table`](#niceops_propagation_table) ·
  [`niceops_build_atmosphere`](#niceops_build_atmosphere) ·
  [`niceops_urban_table`](#niceops_urban_table) · [`niceops_point_field`](#niceops_point_field) ·
  [`niceops_fdtd_field`](#niceops_fdtd_field) · [`niceops_pe_compare`](#niceops_pe_compare) · [`niceops_benchmark`](#niceops_benchmark) · [Plot
  scripts](#plot-scripts)
- [Performance](#performance): [How much gets evaluated](#how-much-gets-evaluated) ·
  [`--db_cutoff`](#--db_cutoff) · [Threads](#threads) · [Relative cost](#relative-cost)
- [Plotting](#plotting)
- [Comparing with older results](#comparing-with-older-results)

## Overview: how a prediction is made

NICEOPS is a source–path–receiver model, run for every trajectory state and every observer: on
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

## Running NICEOPS

### Quick start

With no arguments NICEOPS runs the demo case: `tracks/S-76D/197421AC.csv` against
`databases/S-76D_load_factors.nod` and `databases/S-76D.json`, writing `integrated/197421_LAE.csv`. Run from
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
| `--operations` | A CSV of operations, one row per track flown, with its count and optionally its database, vehicle, clock, lateral offset and air ([An operations file](#an-operations-file)). Replaces `-t`. |
| `-d`, `--database` | netCDF noise database. Once, or once per trajectory. Default `databases/S-76D_load_factors.nod`. |
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
| `--path_smoothing` | How far, rms per axis, the path fitted through a WGS84 track may sit from its positions. Default 3.048 m (10 ft). Zero interpolates exactly. |
| `-s`, `--step` | Resampling step, s (default 0.5), for a WGS84 track whose sampling is uneven or coarser than this. A finer track is never coarsened, and an AC track is used as sampled. Not the SEL integration step. |
| `--trim` | Crop each trajectory from this x to its lowest or slowest point. Off by default; it changes which part of the flight is integrated. |
| `--start_utc` | Seconds since UTC midnight at *t* = 0 of a WGS84 track or static operation, for the day-night periods. Required with `-D`, `--lden_file`, `--cnel_file` or a period level when the track has no `utcsec` column (a static operation never has one), unless an operations file dates the row (`local_time`) or gives period counts; ignored when the track has one. Once, or once per trajectory. |
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
| `--scattering` | Fill refraction shadows with the energy turbulence scatters into them, at HARMONOISE strength γ<sub>T</sub>, in place of the `--shadow_floor` hold; or, with `--buildings`, the buildings' shadows ([Turbulence in the buildings' shadows](#turbulence-in-the-buildings-shadows)). Default 0 (off). Needs `--refraction` and the ground effect, or `--buildings`. |
| `--shadow_model` | `--refraction stratified` with the ground effect: `legacy` (default) or `resolved` ([Low sources near grazing](#low-sources-near-grazing-the-resolved-shadow)). |
| `--shadow_phase` | `resolved`: the vertical phase, radians, a band must gather to resolve the profile near the ground. Default 1. |
| `--cramped_ground` | `resolved`: how a band whose boundary is within two diffraction lengths of the source is filled: `series` (default), `straight` or `ray`. |
| `--tunneling` | `--refraction stratified`: add the rays that turn back just above the microphone, through an Airy tail per band. Off by default ([Micro-shadows](#micro-shadows-and-skip-zones)). |
| `--skip_fill` | `--refraction stratified`: fill the skip zone short of where rays from aloft first land with the caustic's Airy tail. Off by default. |
| `--ray_stride` | `--refraction stratified`: solve each path's layered rays only every this many trajectory states and interpolate the states between where the two agree ([Refraction through the profile's layers](#refraction-through-the-profiles-layers)). Default 4; 1 solves every state. |

#### Ground, receivers and terrain

| Option | Meaning |
| --- | --- |
| `--ground_effect` / `--no-ground_effect` | The ground reflection, on each band ([Ground effect](#ground-effect)). On by default wherever the spectral path runs, except with `--buildings`, where it is off unless asked for. Asking for it with `--no-spectral` is an error. |
| `--flow_resistance` | Ground flow resistivity, kPa·s/m². Default 225 (grass). 30 snow, 650 soil, 3000 dirt, 10000 concrete, 50000 asphalt. |
| `--impedance_model` | `delany_bazley` (default), `miki`, `variable_porosity` or `delany_bazley_layer`. |
| `--porosity_rate` | `variable_porosity`: the porosity's rate of change with depth, 1/m whatever `--units` says. Default 0. |
| `--layer_depth` | `delany_bazley_layer`: depth of the hard-backed layer. Default 0.05 m (0.164 ft). |
| `--roughness` | RMS roughness height of the ground; the coherent reflection falls by exp(−2(kσ cos θ)²). Default 0. |
| `--turbulence_model` | `chessell` (default, with `--turbulence`) or `harmonoise` (with `--gamma_t`): how turbulence decorrelates the direct and reflected paths. |
| `--turbulence` | Chessell's turbulence parameter, rad·s/√m, typically 0 to 16e-4. Default 0. |
| `--gamma_t` | HARMONOISE turbulence strength γ<sub>T</sub>, m<sup>−2/3</sup>: about 1e-6 moderate, 1e-5 strong. Default 0. |
| `--receiver_height` | Height of a point receiver's microphone. Default 1.2 m (3.937 ft). An observer file's fifth column overrides it per observer. Also places the microphone for `--propagation_table`. |
| `--ground_classes` | A CSV of ground surfaces keyed by an integer id, which an observer file's sixth column picks per microphone ([Mixed ground](#mixed-ground)). Needs `-b` with that column. Not with `--refraction` or `--propagation_table`. |
| `--plate_table` | Observers are flush microphones on a ground plate: a netCDF table built for the run's ground ([Ground-plane microphones](#ground-plane-microphones)). Needs the ground effect. |
| `--plate_mic_azimuth` | Which way the microphone is offset on its plate: `outboard` (default, away from the track) or degrees in the local frame from +x toward +y. |
| `--exact_ground` | Evaluate the ground model (or the plate's reflection coefficients) at every observer and state rather than from per-trajectory tables. Slower; for checking the tables. |
| `--propagation_table` | A propagation loss table (netCDF) that replaces the ground and ray models wherever it covers a pair ([Propagation loss tables](#propagation-loss-tables)). Spectral path only. |
| `--urban_table` | Urban receivers: for each microphone it holds, what the street around it does to the direct path, and its facades' first-order reflections ([Urban receivers](#urban-receivers)). Built by `niceops_urban_table`. Spectral path only; not with `--terrain`, `--refraction`, `--propagation_table` or `--db_cutoff`. |
| `--buildings`, `--wall_absorption`, `--floor_absorption` | Buildings around the microphones: every microphone hears the source's own rays among them, in place of the open-ground model ([Buildings around the microphones](#buildings-around-the-microphones)); the floor rigid, absorbing, or with `--ground_effect` the ground model's at every bounce. Spectral path only, with a database carrying spectra. |
| `--facade_scattering`, `--facade_scattering_spectrum` | With `--buildings`: the share of what the walls and roofs reflect that they scatter diffusely, where a building gives none (default 0, or 1 with a spectrum or `--facade_relief`), and a CSV of `band_hz` and `s` that multiplies it band by band. |
| `--facade_relief`, `--relief_scattering`, `--grating_reflections` | With `--buildings`: walls with a periodic relief, reflecting into its grating orders ([Walls with a periodic relief](#walls-with-a-periodic-relief)): one table of orders or more; where what they scatter beyond their specular order goes, `hybrid` (default), `orders` or `lambert`; and a grating ray's reflections besides its relief (default 2). Not with `--facade_scattering_spectrum`. |
| `--roofs` | With a CSV `--buildings`: pitched roofs (LoD2), one row per roof plane ([Pitched roofs](#pitched-roofs)). |
| `--edge_source` | With `--buildings`: diffraction by the edge-source integral, on the CPU, as `niceops_point_field --edge_source` takes it. Refused with `--gpu`. |
| `--patch_size`, `--diffuse_radius` | With scattering facades: the diffuse field's patches' size at most (default 4 m (13.123 ft)), and the horizontal distance from a microphone within which a building with a wall there scatters (default the rays' own reach, `--path_radius`). |
| `--ray_coherence` | With `--buildings`: how rays leaving the source in different directions add: `cone`, `measured`, `coherent` or `energy`. Default `measured` where the database records its vehicle's coherence, `cone` otherwise ([How the rays add](#how-the-rays-add)). |
| `--coherence_angle` | With `--buildings`: the cone's half-angle, degrees, within which rays leaving the source add with their phases. Default 1. Given alone, it selects `cone`; refused with any other `--ray_coherence`. |
| `--phase_interval` | With `--buildings`, under `cone` or `measured`: also write the SEL at each end of this percent interval over the source's unknown phases, as `<output>_phase_lower` and `<output>_phase_upper` ([The phase interval](#the-phase-interval)). |
| `--state_tolerance` | With `--buildings`: dB by which the buildings' effect, taken linearly between the states evaluated, may miss a state's A-weighted level. Default 0.02; 0 evaluates every state ([Not every state is evaluated](#not-every-state-is-evaluated)). |
| `--exact_shared_faces` | With `--buildings`: take the integral over a flat roof between its two edges, or down a step between two roofs' edges, at every wavenumber, not only where the cascade of the two edges cannot stand in for it. The reference; much slower. |
| `--no_profile_pe` | With `--buildings`: over a row of flat-roofed blocks with a street (or a roof lower than both its neighbors), crossed near grazing, the cascades' rays rather than the row's profile marched by a parabolic equation, for comparison, or as a fast mode ([Speed and fast modes](#speed-and-fast-modes)). |
| `--tilted_step_split` | With `--buildings`: across the line through a ridge or roof edge and the edge of the lower roof it steps down to, split the field's jump either side rather than take it as the rays leave it, which steps where the microphone crosses the line. It removes that seam from maps behind a step, and about doubles those microphones' cost. Not taken under `--ray_coherence energy`, which has no jump to split. Off by default: fitted on one step, and on another worse above the line ([FIDELITY](FIDELITY.md#a-point-source-among-buildings), "Across a step"). |
| `--path_radius` | With `--buildings`: the horizontal distance from a source or microphone within which buildings' edges and walls enter its rays. Default 150 m (492.13 ft); less is a fast mode, at a cost in dB that depends on the layout ([Speed and fast modes](#speed-and-fast-modes)). |
| `--gpu`, `--no-gpu` | With `--buildings`: form and sum each state's rays on the GPU (Metal, on a Mac), the default wherever there is one and the rays add `cone`, `coherent` or `energy` without `--phase_interval`; `measured`, the interval and `--edge_source` are the CPU's. `--gpu` asks for it and is refused where it cannot run; `--no-gpu` keeps the CPU ([On the GPU](#on-the-gpu)). |
| `--terrain` | A bare-earth elevation grid (netCDF, for example USGS 3DEP) ([Terrain](#terrain)). Needs the run's frame: `--origin` (with `--frame_bearing` for rotated axes; the AC tracks use 210), or a geodetic track. Not with `--refraction`. The three rows below are refused without it. |
| `--terrain_offset` | Added to every terrain height to put the grid on the trajectory's vertical datum. Default 0. |
| `--terrain_step` | Spacing of the terrain samples each path's mean ground plane is fitted to. Default the grid's own spacing, at least 5 m. |
| `--observers_on_terrain` | With `--terrain` and `-b`: take each observer's ground height from the terrain rather than the file's z. |
| `--terrain_screening` | With `--terrain`: `none` (default), `cnossos` or `kirchhoff`, two ways of screening a path the terrain rises above ([Terrain](#terrain)). Spectral path only. |
| `--screening_min_height` | With `--terrain_screening`: how far the terrain must rise above a path's line of sight to screen it, in `--units`. Default 0.5 m (1.640 ft), above a terrain grid's texture. |

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
| `--nax_file`, `--nax_threshold` | Write each observer's NAx: the number of operations whose own LAmax reaches the threshold, dBA, each trajectory counting as many operations as it stands for. The threshold has no default and is required. |
| `--ta_file`, `--ta_threshold` | Write each observer's time above: seconds of reception time the weighted level exceeds the threshold, dBA, summed over operations (each trajectory times its count). The threshold has no default and is required. |
| `--csel_file` | Write the C-weighted SEL (LCE) per observer. One name, or one per trajectory (`<track>_LCE.csv`). Spectral path only; refused with `--db_cutoff` or `-a`. |
| `--epnl_file` | Write each observer's EPNL, PNLTM, duration and tone corrections (14 CFR 36 A36.4). Once, or once per trajectory (`<track>_EPNL.csv`). Needs the spectral path and the 24 bands 50 Hz–10 kHz. Refused with `-a`. See [EPNL](#epnl-and-reception-time-frames). |
| `--epnl_masking` / `--no-epnl_masking` | Exclude tones in bands below their noy threshold from the tone correction. On by default; `--no-epnl_masking` is the procedure as written. |
| `--band_frames` | Write every observer's band levels in 0.5 s reception-time frames, as netCDF. Once, or once per trajectory (`<track>_band_frames.nc`). Needs the spectral path. Refused with `-a`. |
| `--coverage_file` | Write, per observer, the energy-weighted share of its level that came from sphere cells the database did not measure. Needs a database that records `coverage`. Straight rays only: refused with `--buildings`, `--refraction` or `-a`. |
| `--provenance_spheres` | With `--provenance_file`, also write every sphere's share of every observer's SEL energy to a netCDF file ([Output files](#output-files)). |
| `--provenance_file` | Write, per observer, where its SEL energy came from: the shares from states clamped on each database axis or held at a load-factor limit, and the sphere that supplied the most ([Output files](#output-files)). One database per run; refused with `-a`. |

#### Threads, culling and diagnostics

| Option | Meaning |
| --- | --- |
| `-j`, `--threads` | Worker threads for the observer loop. 0 (default) uses one per hardware thread; 1 runs serially. The output is the same at any count. |
| `--db_cutoff` | Skip evaluations that provably fall below this level, dB. Off by default; read [the caveat](#--db_cutoff). |
| `--export_sphere` | Interpolate a whole noise sphere at `--sphere_condition` and write it as a one-group `.nod` database. |
| `--export_aam` | The same sphere in AAM/RNM format, in AAM's flight-path frame (fly the AAM deck at AOA 0). Needs a database with spectra. |
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
`--path_smoothing`, `--terrain_offset`, `--terrain_step`, `--screening_min_height`, `--roughness`,
`--layer_depth`, `--roughness_length`, `--surface_layer_depth`, `--patch_size`, `--diffuse_radius`,
`--path_radius` and the anemometer height in `--wind`), every file carrying a length
(trajectories, observers, `--atmosphere` heights, a CSV `--buildings` and `--roofs`, level output,
histories),
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
| `--tunneling`, `--skip_fill` or `--ray_stride` without `--refraction stratified` | They are stratified-ray options. |
| `--terrain` with `--refraction`, or without a frame | The ray models are built over flat ground; terrain needs to know where x/y lie. |
| `--terrain_screening` without `--terrain`, or with `--no-spectral` or a database without spectra; `--terrain_screening cnossos` with `--propagation_table` and `--no-ground_effect` | The edges are the terrain's; the diffraction is per band; a path `cnossos` screens leaves the table for the ground model. |
| `--plate_table` without the ground effect, or over a different ground | The table is built for one ground. |
| `--ground_classes` with grid observers, without an observer file's ground column, or naming a class the file does not define; the column without `--ground_classes` | Each names the other. |
| `--ground_classes` with `--refraction`, `--propagation_table`, `--no-ground_effect`, `--no-spectral` or a database without spectra; flush microphones on more than one class, or on one the plate table was not built for | The ray models need the ground along the path, a table and a plate table each hold one ground, and the classes are a ground for the ground effect. |
| `--latlon`, `--contour_file`, `--nose heading` without a frame; `--contour_file` or `--raster_file` with `-b`; `--contour_file` with `--contour_metric sel` over several trajectories | Nothing to turn compass directions or invert positions with; scattered observers have no grid to contour; one SEL per trajectory. |
| `-D`, `--lden_file`, `--cnel_file` or a period level on a track with no clock, unless an operations file gives period counts | Day or night would be made up. |
| `--laeq_file` without `--laeq_window`, or an empty window; `--nax_file` or `--ta_file` without its threshold | There is no default period or threshold. |
| An option that shapes another without it: `--min_load_factor` without `two_sided`; `--shadow_phase` or `--cramped_ground` without `--shadow_model resolved`; `--lamax_weighting` without `--lamax_file`, `--nax_file` or `--ta_file`; `--roughness_length`, `--surface_layer_depth`, `--stability`, `--profile_levels` or `--write_atmosphere` without `--wind`; `--terrain_offset`, `--terrain_step` or `--observers_on_terrain` without `--terrain`; `--screening_min_height` without `--terrain_screening` | Each shapes an option that is not given. |
| A `--lamax_weighting` on a track sampled coarser than half its time constant | The "maximum" would be the largest sample. |
| `--speed_reference air` without a wind | Calm air is not assumed. |
| Two trajectories writing one file; a list of `-d`, `-v`, `-o` that is neither one nor one per track | Outputs would overwrite one another. |
| `-t` with `--operations`; a per-trajectory option given more than once with it; a row with no database or vehicle and no `-d` or `-v` | The file names the operations; the demo aircraft never fills a row. The file's own refusals are in [An operations file](#an-operations-file). |
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
| Operations (`--operations`) | CSV, one row per track flown ([An operations file](#an-operations-file)). |
| Atmosphere (`--atmosphere`) | CSV profile ([Stratified atmosphere](#stratified-atmosphere)). |
| Ground classes (`--ground_classes`) | CSV of surfaces keyed by class id ([Mixed ground](#mixed-ground)). |
| Terrain (`--terrain`) | netCDF: 1-D `latitude` and `longitude` (degrees, strictly increasing) and `elevation(latitude, longitude)` in meters whatever `--units` says; an optional global `vertical_datum` text attribute. Missing nodes are NaN, `_FillValue` or `missing_value` ([Terrain](#terrain)). |
| Plate table (`--plate_table`) | netCDF from panam ([Ground-plane microphones](#ground-plane-microphones)). |
| Propagation table (`--propagation_table`) | netCDF ([Propagation loss tables](#propagation-loss-tables)). |
| Urban table (`--urban_table`) | netCDF from `niceops_urban_table` ([Urban receivers](#urban-receivers)). |
| Buildings (`--buildings`, `--roofs`) | CSV or GeoJSON footprints, and a CSV of roof planes ([Buildings around the microphones](#buildings-around-the-microphones)). |

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
| `heading` | Optional compass heading of the nose, degrees from true north (`hdg` serves too); a unit other than degrees, a header naming a magnetic heading, and a second heading column are refused, wherever the header says it (`heading (rad)` and `heading_rad` alike). It sets the hover nose and `--nose heading` ([The sphere's nose](#the-spheres-nose)). Resampled linearly, through north the short way. |

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
`heading`. The states are sampled at 50 Hz, at 1 g and zero velocity. A static operation has no
clock of its own, so `-D`, `--lden_file`, `--cnel_file` and the period levels refuse it unless
`--start_utc`, or an operations file's `local_time`, dates its start.

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
the terrain. `-b FILE` reads scattered observers instead, in `--units`: `x,y,z`, `id,x,y,z`,
`id,x,y,z,height` or `id,x,y,z,height,ground`, the fifth column being each microphone's height
above the ground ([Receivers](#receivers)) and the sixth the ground class it stands on
([Mixed ground](#mixed-ground)). The first line is always a header, and is discarded even when it
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

### Source coherence

A database may record its vehicle's coherence between emission directions, which `--buildings`
weights its rays by, as optional root attributes written by panam (its file-format
documentation lists them): the main rotor's blades and sense of rotation, the band up to
which its delay holds, the decay of the mean coherency and of the coherence with frequency times
angle, and each band's coherence between distant directions. A database that records any of them
is taken with `--ray_coherence measured` unless the option says otherwise, the rotorcraft
measurements' values filling what it leaves out; the rotor's delay comes from the database alone.
Under `--ray_coherence cone` the record's spread and tonal fractions still set `--phase_interval`'s
width.
One that records none is taken with the cone. The rotor turns at the database's `rotor_rpm` where
it gives one, and at the vehicle file's tip speed over radius otherwise. `--export_sphere` carries
the record and `rotor_rpm` on, so a sphere cut from a database is taken as the database is. A value NICEOPS cannot
use, such as a band count that is not the database's, is refused when the database loads.

### Shipped databases

| Files in `databases/` | Contents |
| --- | --- |
| `AS350_`, `B407_`, `Be206L3_`, `EC130_`, `RO-44_`, `RO-66_fixed_load_factor.nod` | **Spectral**, 10 Hz–10 kHz, built by panam from the 2017 Noise Abatement flight test. The current, validated databases. |
| `S-76D_fixed_load_factor.nod` | **Spectral**, built by panam from the 2019 Noise Abatement flight test at Coyle Field: level flight at 45–155 kt, descents to −10.7° and climbs at +11° to +22°. Its approaches' flare and climb-out lie between these conditions, and a run there warns ([Choosing spheres](#choosing-spheres-blending-and---blend)). |
| `S-76D_load_factors.nod` | Broadband, with discrete load factors, from the same spheres. The demo database, which the goldens use. |
| `AW139.json`, `B205.json`, `B407.json`, `S-76D.json` | Vehicle files. |

Of the 2017 databases, **B407 is the reference** and the most validated; its spheres are
depropagated along curved rays through each run's own atmosphere. Its hover is synthesized from
the slowest level pass and corrected by hovers measured 0.5–3° and 18–25° below the horizon, so
treat the hover from 3° to 17.5° and under the aircraft as an interpolation. The other five take
the B407's recipe without a hover correction. No approach was used to build any of them: the
approaches are the validation set ([FIDELITY](FIDELITY.md#shipped-databases-approach-scores-per-aircraft)).
Each file records its build at the root (`speed_reference`, the build atmosphere, `hover_source`)
and `DOPPLER_SHIFT_REMOVED` and `coverage` per group. The spheres keep the measured,
Doppler-shifted spectra, and NICEOPS applies no frequency shift. The rebuilds are in
[HISTORY](HISTORY.md#the-2017-noise-abatement-databases).

### Spectra

A database may carry per-band spectra beside the overall A-weighted level: a `frequency`
dimension and an `amplitude(PHI, THETA, frequency)` array per sphere. It then runs the spectral
path by default, and serves every option that works on bands. `--export_aam` writes a sphere in
the AAM/RNM hemisphere format, whose payload is `AMPLITUDE(PHI, THETA, FREQUENCY)`, so it needs
spectra. AAM fixes its sphere to the body, pitched by the deck's angle of attack plus the flight
path angle, so the export is pitched nose-up by the condition's flight path angle from NICEOPS's
level sphere, and records `sphere_frame = flight_path`: an AAM deck flown at AOA 0 then hears
what NICEOPS hears. A database is level; one that records another `sphere_frame` is refused.
Third-octave bands (10 Hz to 10 kHz) are written as AAM's broadband sphere type (`BB`), any
other band set as narrowband (`NB`). The export is a lower hemisphere, `PHI` -90° to 90°, as
AAM's spheres are, with an empty band written as AAM's missing level, -999. Of the shipped databases, the six `_fixed_load_factor`
files of the 2017 aircraft carry spectra; the S-76D demo database does not, so the demo case
runs broadband.

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

No shipped database warns. A sphere more than half defective has a defective
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

The axis ranges include spheres a database builder made rather than measured: panam extends the
measured path angles to steeper and shallower ones, and synthesizes a hover, and records which in
each sphere's `condition_origin`. A state can therefore sit inside every range and still be far
from any measured condition. Each run logs the mean blend weight its states put on measured,
extended and synthesized spheres, and warns, with their number and time span, about states that
take less than half their weight from measured spheres (a state on the hover itself excepted: a
hovering trajectory is meant to use it). `--provenance_file` gives each observer's share of its
exposure from the extended and the synthesized spheres. A database without `condition_origin`
counts every sphere as measured.

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
  the lowest level that level is held to the ground, with a warning, so give the profile a level
  at the ground ([`niceops_build_atmosphere`](#niceops_build_atmosphere) always writes one).
- **No calm level above the ground under wind.** `--refraction stratified` takes one as a thin
  layer that can put levels 20 dB or more low, and the run warns of it. A profile calm at every
  level is still air and does not warn.
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
(`--humidity 0`) does not replay, since `--atmosphere` refuses a humidity column that never
exceeds 1 percent as a likely fraction, and the run warns. On the 2017 windy day, `--wind` from the ground stations'
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
dB louder than grass, and grass roughly cancels its own reflection gain against the dip
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

### Mixed ground

When microphones stand on different surfaces, for example a concrete pad among grass,
`--ground_classes FILE` gives each surface as a row of a CSV, and an observer file's sixth column
(`id,x,y,z,height,ground`) names the class of each microphone:

```csv
class,name,impedance_model,flow_resistance
0,grass,,
21,pad,delany_bazley,10000
```

- **Columns.** `class` is required: an integer id from 0 to 2147483647, in any order (land-cover
  codes work). The others are optional: `name` (for the log), `impedance_model`,
  `flow_resistance` (kPa·s/m²), `porosity_rate` (1/m), `layer_depth` and `roughness` (in
  `--units`). A column left out, or a cell left empty, takes the command line's value, so a file
  lists only what differs; above, class 0 is whatever the command line says.
- **Turbulence is the air's**, so it stays run-wide (`--turbulence`, `--gamma_t`), and a
  turbulence column is refused. The receiver height comes from the observer file.
- **Flush microphones** on a plate (`--plate_table`) must all be on one class, the ground the
  plate table was built for.
- **The limit.** Each microphone is evaluated as if the whole ground were its class. Use a class
  where the patch extends well beyond the microphone toward the aircraft: tens of meters for a
  1.2 m pole at grazing incidence, less for a high source. A flush plate needs only a patch
  larger than the plate. A small pad under a pole microphone is given the full effect of a
  concrete plain, which overstates it at long range.
- **A land-cover class is not an impedance.** Ground mapped as built-up or paved is seldom a
  concrete plain: measured at grazing, such a patch behaved like ground of about 1000–2000
  kPa·s/m², and concrete's 10,000 or more overstated it by several dB
  ([FIDELITY](FIDELITY.md#ground-classes)). Give a hard class a moderate flow resistance unless
  the surface is known to be sealed pavement.
- **What it does not combine with.** `--refraction` (its terms depend on the ground along the
  path), `--propagation_table` (a table holds one ground), grid observers (they have no class
  column), `--no-ground_effect` and the broadband path. `--terrain` works.
- **Cost.** One ground table per microphone height and class, about 10 MB each; above 16 the run
  warns and `--exact_ground` is the alternative. Classes defined but not used are only logged.

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
`--plate_table`, height 0 is a plate and anything else a point receiver. A sixth column puts each
microphone over its own ground class ([Mixed ground](#mixed-ground)). The log reports how many
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
the mean plane of CNOSSOS-EU, which NORAH2 uses over terrain
([THEORY §6.3](THEORY.md#63-terrain)). Spreading, absorption and the sphere lookup still take the
real three-dimensional path. A level terrain file reproduces the run without it, bit for bit.

**Datums.** The local z = 0 is `--ground_elevation` on the trajectory's vertical datum, so

    z_ground(x, y) = H(latitude, longitude) + --terrain_offset − --ground_elevation

with H the grid's height, interpolated bilinearly. NICEOPS does not convert between datums:
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

**Screening.** By default terrain does not screen or diffract: a ridge between source and
receiver puts the source below the mean plane, and its height is set to zero, so the level there
is too high. The run warns when it happens, with how many observers had such paths (a receiver in
a dip raises the same flag). `--terrain_screening` screens instead, by one of two models, on the
same paths: those on which the terrain rises more than `--screening_min_height` above the line of
sight. Every other path is evaluated as without it, so level ground, a uniform slope, or hills no
path crosses give the unscreened run exactly. Use it where ridges, bluffs or berms stand between
low flight and the receivers.

- `cnossos` takes the CNOSSOS-EU diffraction over the terrain's convex outline between source
  and microphone, with the ground effect on each side of it in place of the whole path's, capped
  at 25 dB ([THEORY §6.4](THEORY.md#64-terrain-screening)). It is the method European noise
  mapping prescribes; choose it to match that.
- `kirchhoff` keeps the level the path would have unscreened, with its ground over the mean
  plane, its plate or its table, and takes off the Fresnel-Kirchhoff loss of one knife edge,
  standing where the terrain rises highest above the line of sight, as high as both ends see
  over the terrain there ([THEORY §6.4](THEORY.md#64-terrain-screening)). Behind a hill it
  shields less than `cnossos` above the low bands, which is closer to what was measured, and it
  costs little more than the unscreened run.

What both leave out, and where they step:

- Only the vertical plane through source and microphone is taken: no paths around the ends of an
  obstacle, and no buildings or vegetation. The rays are straight. Behind an isolated hill, which
  sound also reaches around its sides, both models overstate the shielding above the low bands
  ([FIDELITY](FIDELITY.md#terrain-screening)); use them for long ridges, bluffs and berms.
- Crests are found at `--terrain_step`, so a sharp ridge is slightly under-screened at the
  default step; a finer step costs about proportionally more.
- As a source drops behind a crest the loss switches on at 4.8 dB (`cnossos`) or 6.0 dB
  (`kirchhoff`), not from zero, so the level steps; that shows in LAmax and time above more than
  in SEL.
- `kirchhoff` stands one edge for all the crests on a path, so two ridges well apart shield less
  than they would one after the other. `cnossos` diffracts over every crest of the outline.
- With `--propagation_table`, a path `cnossos` screens is not read from the table but takes the
  screening with the ground model, which must then be on; unblocked paths keep the table. The two
  physics meet at the edge of the screened region, and the log counts the paths that left the
  table. `kirchhoff` keeps the table on every path.
- Their accuracy, against an exact solution and behind one measured hill, is in
  [FIDELITY](FIDELITY.md#terrain-screening).

Each path's ground is otherwise one plane. `--refraction` is refused, because the ray models are
built over flat ground, and on the broadband path `--terrain` only places the observers. A
`--propagation_table` is read over each path's mean plane
([Propagation loss tables](#propagation-loss-tables)).

**Cost and effect.** Plates keep their tables over terrain; point receivers over real terrain are
evaluated exactly, about 4 times the CPU, which a test array of tens of microphones does not
notice; a ground-model `--propagation_table` with rows spanning their heights above the planes
reads them instead, at about 1.5 times. Screening adds a test per path and keeps about twice the
terrain samples in memory per thread; `cnossos` adds two exact ground evaluations per screened
path, and `kirchhoff` a table lookup per band. At
Amedee terrain moves flush-plate SEL by hundredths of a dB; at a sloped site it can move it by a dB
at low elevation ([FIDELITY](FIDELITY.md#terrain-at-amedee),
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
- **Interpolated states (`--ray_stride`, 4 by default).** The rays to a microphone change little
  from one trajectory state to the next, so they are solved at every fourth state and the states
  between are interpolated. A state whose two solved neighbors differ in kind (lit and shadowed, a
  different set of arrivals, a caustic or boundary passing between them) is solved itself, and
  the level at every state is still formed from its own geometry. It roughly halves a stratified
  run against solving every state, for SEL and LAmax within thousandths to hundredths of a
  decibel; a single sample of a time history can differ by up to a decibel where the solved
  history itself steps from one state to the next
  ([FIDELITY](FIDELITY.md#interpolated-rays---ray_stride)). `--ray_stride 1` solves every state:
  to compare ray models state by state, or to reproduce a stratified result that solved every
  state. Larger strides gain little more, since what is left is the rest of the evaluation.

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
- **Turbulence and roughness come from the table.** A covered pair carries what the table was
  built with, and the run's `--turbulence`, `--gamma_t` and `--roughness` reach only the pairs
  that fall back. Build the table with the run's (`niceops_propagation_table --turbulence`,
  `--gamma_t`, `--roughness`); a run whose settings differ from what the table records is warned.
- **With `--refraction` too**, every state is evaluated rather than windowed from the peak.
- **Over terrain** each pair reads the table at its heights above, and range along, its path's
  mean ground plane ([Terrain](#terrain)), as the ground model does there. A point microphone's
  height above the plane is its own, not its pole's: build the table with receiver rows from 0 to
  above the tallest microphone plus the terrain's relief. A pair whose microphone falls outside
  the rows takes the run's own model, and the run warns with the heights seen; a table with one
  receiver row is warned about before the run. Reading between rows is where a point receiver's
  error is over terrain: the ground's interference fringes in receiver height are finer than a
  practical row spacing above a few hundred hertz. A flush plate reads the lowest row, as on flat
  ground, which needs no extra rows. A table that refracts (a PE table, or any not built by the
  ground model with no atmosphere) is taken to have its atmosphere follow each path's mean plane;
  the run warns when pairs read it over a plane tilted more than 5 degrees. With
  `--terrain_screening cnossos`, a path the terrain blocks takes the screening and the ground
  model instead of the table; `kirchhoff` reads the table and takes its knife edge's loss off it
  ([Terrain](#terrain)).
- **The file** is netCDF, SI: `relative_level(bearing, source_height, receiver_height, range,
  band)`, a coordinate variable per axis, and the attributes `propagation_table_version` (1),
  `method`, `includes_absorption` and `provenance_*`. Anything that fills it can be read
  ([internals/pe.md](internals/pe.md) has the specification).

Interpolating between source heights is where a table's error is: build tables a few wavelengths
apart at the highest band, or densely at the lowest heights, and check them with
`niceops_propagation_table --validation_csv`. Reading a table costs little.

### Urban receivers

`--urban_table` puts microphones in streets. For each microphone it holds, the table gives what
the buildings around it do to the sound arriving from each direction, and the facades that
reflect it ([THEORY §6.5](THEORY.md#65-urban-receivers)). The run evaluates the direct path as it
always does and adds the table's gain, in dB per band, at the source's compass bearing and
elevation from the microphone. That gain carries the shielding where buildings block the direct
path and the street's further reflections. Each first-order facade reflection is evaluated as a
path of its own, with the sphere read in the direction the reflected sound leaves the rotor, and
the paths are summed in energy. Microphones without a table are evaluated as without one.

- **Building one.** [`niceops_urban_table`](#niceops_urban_table) builds tables two ways. For
  microphones in ideal street canyons, two long parallel walls, it takes the street's bearing,
  width, the microphone's offset and the walls' heights; a row of width 0 is an open microphone,
  0 dB with no facades. Above a few hundred hertz that table comes from rays; below, with
  `--wave_below`, from the wave equation, which a street only a few wavelengths wide needs. For
  microphones among real buildings it takes their footprints and roof heights (LoD1, with pitched
  roofs from `--roofs`) and traces the rays among them.
- **Receivers are matched by position**, to a millimeter: a receiver is a microphone (an
  observer raised by its height) in the run's local frame. Build the table from the run's own
  observer file. A table none of whose receivers is in the run is refused; receivers the run does
  not have are warned about.
- **The frame.** Azimuths are compass bearings and the geometry is local, so the run's
  `--frame_bearing` must be the one the table was built with (east when neither gives one).
- **Coverage.** A direction outside the table's elevations, or a NaN, takes the open-ground model
  with no urban terms at all; the log counts these.
- **Interpolation** is bilinear in dB, which smooths the steps at shadow edges over one cell.
- **The street's floor** is the table's, not the ground model's: inside a canyon the floor's
  reflection often arrives where the direct sound is blocked. The table assumes hard pavement
  unless built with `--floor_absorption`; put urban microphones on a hard ground class
  ([Mixed ground](#mixed-ground)) so the open-ground prediction it scales is a street's.
- **What the models hold and leave out.** The walls and floor reflect specularly and the
  building edges diffract, all summed with their phases so that the field stays continuous where a
  reflection or the direct sound is cut off. Among buildings, the roof edges and the buildings'
  convex corners diffract, and a diffracted ray may reflect in the facades on its way to the
  microphone, and with `--facade_scattering` the walls also scatter part of what they reflect
  diffusely, which fills the shadows that specular reflection leaves. Roofs and slopes reflect
  only below the microphone. Diffraction at two edges in turn, an ideal street's ends and
  turbulence in the shadows are not modeled, nor are the low-frequency wave effects unless the
  table has a wave band (ideal canyons only).
- **Not with** `--terrain`, `--refraction`, `--propagation_table` or `--db_cutoff`.

A microphone in a street costs about one extra evaluation per reflecting facade.

### Buildings around the microphones

`--buildings district.csv` evaluates every microphone with the point-source model of
[`niceops_point_field`](#niceops_point_field): the rotor's own rays among LoD1 buildings,
reflected in the facades, the roofs and the floor and diffracted at their edges once, at two in
turn (reflecting in a facade between them or not) or at three roof edges, each reading the sphere
in the direction it leaves the rotor, added as [below](#how-the-rays-add)
([THEORY §6.5](THEORY.md#65-urban-receivers)). Where `--urban_table` takes the source as a distant
plane wave arriving from its direction, this follows the source among the buildings. Use it for a
source near them -- a drone among the roofs, a slow approach to a vertiport -- and for metrics of a
moment, such as LAmax and the numbers above a level, where the table's error does not average out
over a pass as it does in an exposure. Each evaluation costs far more than the open-ground model,
so this is for a street's microphones, not a grid.

**Not with** `--urban_table`, `--terrain`, `--refraction`, `--propagation_table`, `--atmosphere`,
`--wind`, `--db_cutoff` or `--coverage_file`.

#### The district, the microphones and the database

- **The district** is `niceops_point_field`'s: a CSV of footprint vertices (`building`, `x`, `y`,
  `top`, optionally `absorption`, `scattering` and `relief`) in the run's `--units` and frame, or
  GeoJSON with `--origin` and `--frame_bearing`. A building's coefficients apply to its walls and
  its roof, and `--wall_absorption` and `--facade_scattering` to those that give none.
- **Where a microphone stands.** At its observer's point, raised by the observer file's height
  column where it has one, as in any run without a ground model; with a ground model
  ([The floor](#the-floor)), as in the open-ground run, at `--receiver_height` unless the height
  column says otherwise. That puts it on the floor, on a roof, on a wall, or on a pole above any of
  them. An observer whose `z` is 0, or whose point lies on a roof or a wall of the district, stands
  on the district's floor (z = 0), and a microphone on a roof or a wall takes that surface's
  reflection as one on the floor takes the floor's (a microphone flush on a flat roof hears the
  roof's doubling). Any other observer's floor is at its own `z`, as in the open-ground run: the
  buildings stand on that floor for it, and the source is taken as high above it as the track says
  above z = 0. With `--plate_table` a flush microphone sits on its plate on the floor or a flat
  roof, the roof's reflection in place of the ground's; one on a pitched roof's slope takes the
  slope's reflection, as on a flat roof.
- **Refused**, naming the observer or the state: a microphone below its floor, inside a building
  or on a building's edge, or flush on a wall or a slope with a plate (a plate lies level); and a
  source state below a microphone's floor or inside a building
  ([THEORY §6.5](THEORY.md#65-urban-receivers)).
- **The database** must carry per-band spectra, which each ray reads along its own launch; one
  without them is refused.

#### Pitched roofs

`--roofs roofs.csv` (LoD2, with a CSV district) has one row per roof plane: `building` (an id of
the district's), `x`, `y` and `z` (a point the plane passes through, in the run's `--units` and
frame) and `slope_x` and `slope_y` (its rise per unit run along x and along y). A building given
no plane keeps its flat roof.

- **The roof** is the least of its building's planes, capped at its `top`: a top below the
  planes' highest point is a level deck there (a mansard's), and one at or above it changes
  nothing. A least of planes is convex, so a roof with a valley or a step is several buildings
  sharing sides. Along a shared side the wall stands, and the edge diffracts, only where one roof
  is above the other, judged along the side wherever the two roofs' profiles cross or break.
- **What it does.** Its walls stand up to the roof, and its eaves, rakes, ridges and hips diffract
  as wedges of their own angles; every path the flat roofs take applies. A level deck reflects and
  scatters as a flat roof does. A slope reflects specularly, as a mirror in its plane with the
  floor at either end, with walls before or after it. A once-diffracted ray may reflect in walls
  and then a slope on its way from either end to the edge, and a twice-diffracted ray in a slope
  or a level roof between one end and its edge. Slopes scatter diffusely as the walls do, with
  their building's coefficients ([THEORY §6.5](THEORY.md#65-urban-receivers)).
- **How well.** On EDGE's measured pitched scenes the roof planes do better than the same houses
  as flat blocks, at the mean height or at the ridge ([FIDELITY](FIDELITY.md#pitched-roofs)).

#### The floor

- **Without a ground model** the floor is rigid by default, or absorbs `--floor_absorption` of the
  energy at every angle and frequency.
- **With `--ground_effect`** (or `--ground_classes` or `--plate_table`, which name a ground) it is
  the ground model's instead: each ray's every bounce on the floor -- the floor's own ray, and a
  diffracted ray's bounces before, between and after its edges -- takes the spherical-wave
  reflection coefficient Q the open-ground model takes, at that bounce's grazing angle, its leg as
  the image path, in each third of each band ([THEORY §6.5](THEORY.md#65-urban-receivers)). With
  one far building the run is close to the open-ground run: the floor's ray reads the sphere along
  its own launch, and Q changes across a band. A ground floor costs a little more time than a
  rigid one; how much it moves the SELs, and how close it comes to open ground, is in
  [FIDELITY](FIDELITY.md#the-point-source-model-in-a-run). `--floor_absorption` with a ground
  model is refused.
- **Each microphone's floor** is its ground class (`--ground_classes`), the ground at the
  microphone taken for every bounce, as the open-ground model takes it.
- **A flush microphone** (height 0) with `--plate_table` is on its plate: a ray reaching it
  directly takes the plate's P<sub>d</sub>, and the floor's ray beside it Q P<sub>r</sub>, at the
  elevation and azimuth the ray arrives from, which is the open-ground plate's level where only
  those two arrive.
- **Turbulence.** The ground model's `--turbulence` and `--gamma_t` are refused: how the rays'
  phases relate is the source's ([How the rays add](#how-the-rays-add)), and turbulence fills the
  shadows through `--scattering` ([below](#turbulence-in-the-buildings-shadows)). The diffuse
  facades and that fill carry energy alone, and take the ground's diffuse-field reflection, band by
  band.

#### Facades that scatter

A wall or roof whose scattering coefficient s is above 0 reflects specularly only ρ(1 − s) of what
reaches it and scatters the rest of ρ diffusely, by Lambert's law, as `niceops_point_field
--scattering` takes it ([THEORY §6.5](THEORY.md#65-urban-receivers)).

- **The diffuse field.** Each microphone has its own, over the buildings with a wall within
  `--diffuse_radius` of it, walls and roofs cut into patches of at most `--patch_size`. The source
  lights each patch by its specular paths there, off the walls, the floor and one roof or slope,
  every path reading the sphere along its own launch as a ray does, and the field adds to the rays
  in energy, band by band. A building beyond the radius reflects its whole ρ.
- **By band.** `--facade_scattering_spectrum` gives s band by band, as
  `tools/profile_scattering.py` writes it, times each building's coefficient; it must reach every
  band the metric sums (`--band_limits`).
- **Cost.** Each microphone's field is built once for the run, so the cost grows with the
  microphones and the patches, not with the states
  ([FIDELITY](FIDELITY.md#the-point-source-model-in-a-run) has the times).
- **How well.** For a corrugated container's s(f), on NLR's DigiCity, see
  [FIDELITY](FIDELITY.md#scattering-walls-point-source).

#### Walls with a periodic relief

Ribs, a sawtooth in plan or a container's corrugation do not scatter by Lambert's law: they reflect
into grating orders, each on the cone of the wave's elevation
([THEORY §6.5](THEORY.md#65-urban-receivers)).

- **The tables.** `tools/profile_scattering.py --orders` writes a relief's orders by boundary
  elements from its profile. `--facade_relief` takes one table or more, a building carrying the one
  its `relief` column names (from 1; 0 for none), or without the column the first wherever it
  scatters (`--facade_scattering` or its `scattering` column, 1 by default with this). Every wall
  of the building carries it, its +along the wall's base direction, which matters only for a
  relief that is not symmetric. Not with `--facade_scattering_spectrum`.
- **The rays.** Specular rays take the relief's specular order at the incidence they meet it at;
  rays that diffract reflect in it whole. A grating ray reflects in `--grating_reflections` walls
  besides its relief (default 2).
- **The rest** goes by `--relief_scattering`: `hybrid` (the default), each wall's first scattering
  as grating rays straight to the microphone and what it scatters onto other walls diffusely;
  `orders`, the grating rays alone; `lambert`, all of it diffusely, by Lambert's law at each
  lighting path's incidence. Against EDGE's measured reliefs and DigiCity the hybrid beats a
  relief's random-incidence s(f) everywhere, the orders alone match the FDTD in a canyon and fall
  short in a shadow, and Lambert's law at incidence the reverse
  ([FIDELITY](FIDELITY.md#scattering-walls-point-source)).
- **Other walls.** With a relief in the run, a wall that carries none reflects specularly whatever
  its scattering coefficient: the coefficient then only says which buildings carry the relief.

#### Turbulence in the buildings' shadows

`--scattering` γ<sub>T</sub>, off by default. Where the buildings shield a microphone from the
source, turbulence scatters energy into the shadow, most at high frequency: Forssén et al.'s
engineering model for a thin screen, the roof top farthest above the line from source to
microphone ([THEORY §6.5](THEORY.md#65-urban-receivers)), added in energy and held to open
ground's level. A microphone in sight of the source is untouched. It matters where a microphone
stays shielded through an event -- a courtyard, a street behind a row under a low approach -- and
hardly on a pass that each microphone sees at some point
([FIDELITY](FIDELITY.md#the-point-source-model-in-a-run)). The strength is the weak input:
γ<sub>T</sub> of 10<sup>-6</sup> (moderate) to 10<sup>-5</sup> (strong) is
C<sub>v</sub><sup>2</sup> of 0.016 to 0.16 m<sup>4/3</sup> s<sup>-2</sup>, where measurements
beside a barrier in wind have reported 0.4 to 17, and a factor of ten in it is 10 dB of scattered
level.

#### How the rays add

The spheres carry levels, not phases, so how coherent the source is between two rays' launch
directions is a model ([THEORY §6.5](THEORY.md#65-urban-receivers)). The run's log says which it
used.

- **`cone`** (the default for a database that records no coherence) assumes nothing of the
  source: rays add in energy, except those leaving within `--coherence_angle` (1°) of each other,
  which are one emission -- a flush microphone's direct ray and the floor's, a diffracted ray and
  the ray it stands in for at a shadow boundary -- and add with their phases, so the field stays
  continuous and the pressure doubles at a flush microphone. A microphone above the floor whose
  reflection leaves more than the angle from its direct ray takes the two in energy, without their
  interference.
- **`measured`.** What the database records comes first. A database may record its vehicle's
  coherence ([Source coherence](#source-coherence)); the run then takes it: each pair's cross term
  weighted by a coherency that falls with frequency times the angle between their launches, and,
  where the database gives the main rotor's blades and sense of rotation, the delay a pattern
  turning with the rotor puts between them in the bands up to 4.5 blade-passage frequencies.
  `--ray_coherence measured` takes the rotorcraft measurements' coherence without the rotor where
  the database records none.
- **`coherent`** adds every ray with its phase, as if the rotor radiated in phase in every
  direction; **`energy`** adds every pair in energy, which breaks the field's continuity at shadow
  boundaries.

#### The phase interval

`--phase_interval 90` writes, beside the output, the SEL with every state at the 5th percentile of
its level over the phases the model leaves unknown (`<output>_phase_lower`), and at the 95th
(`<output>_phase_upper`). It goes with `cone` or `measured`. Names that collide with another output
are refused before the acoustic calculations. The states are put at the same end together, since
one vehicle's phases do not change from moment to moment, so for a pass the interval is wide rather
than narrow. It is widest at low frequency and close to the rotor, where the source's tones may
keep a stable but unknown phase between directions that the energy sum leaves out; broadband rays
far apart leave nothing to be unsure of. How far that coherence reaches between directions, and
each band's tonal share of it, come from the database's record where it gives them, under the cone
too; the rotorcraft measurements' values otherwise. The log says which. It costs more than the
level alone.

#### Not every state is evaluated

The buildings' effect on each band -- the band's level less the direct ray's in free field --
changes slowly along a track, while the source's own level changes from one state to the next. So
per microphone the rays are evaluated at states 32 apart, and then wherever the effect taken
linearly between two evaluated states misses the level halfway between by more than
`--state_tolerance` (0.02 dB by default), until it does not. The other states take the source
exactly, from the sphere at their own state along the direct ray, plus the effect interpolated. A
state where the source makes a band loud whose effect is known at neither evaluated state around
it is evaluated outright. `--state_tolerance 0` evaluates every state. There is no peak search.

#### Diffraction over shared faces, and by the edge-source integral

- **Flat roofs and walls between two corners**, and a step from one roof's edge down to a lower
  roof's, are crossed by the uniform form's integral only where the cascade of the two edges'
  coefficients cannot stand in for it -- the lower bands, and sources or microphones near the line
  between the edges -- at wavenumbers two thirds of an octave apart, and the cascade carries it in
  between and beyond. `--exact_shared_faces` takes the integral at every wavenumber.
- **Across a step's line.** `--tilted_step_split` splits the field's jump either side of the line
  through a ridge or roof edge and the edge of the lower roof it steps down to (see the option
  table); off by default.
- **The edge-source integral** (`--edge_source`), as `niceops_point_field --edge_source` takes it
  ([niceops_point_field](#niceops_point_field)), carries the edges' ends and the transitions the
  cascade misses, on the CPU (with `--gpu` refused). It costs far more a pair than UTD, so a
  trajectory's few thousand pairs take minutes, not seconds
  ([FIDELITY](FIDELITY.md#a-point-source-among-buildings)).

#### On the GPU

Where the machine has one, each evaluated state's rays are formed and summed there, in single
precision, while the CPU traces the next states' rays; the levels agree with the CPU's closely
([FIDELITY](FIDELITY.md#the-point-source-model-in-a-run)). `--no-gpu` keeps them on the CPU, and
the measured coherence, `--phase_interval` and `--edge_source` always do: `--gpu` with
`--ray_coherence measured`, `--phase_interval` or `--edge_source` is refused, and with a database
whose recorded coherence selects the measured model the run stops at that track.

#### Speed and fast modes

Two options trade dB for time, and neither is the default. `--path_radius` below its 150 m leaves
out the buildings farther from each end, which pays where they spread out and costs nothing where
everything is near. `--no_profile_pe` takes the cascades' rays over a row of roofs instead of the
row's march, most of the time wherever such rows are crossed. What each costs is measured in
[FIDELITY](FIDELITY.md#the-point-source-model-in-a-run).

## Metrics and output files

### Output files

Each metric file (`-D`, `--lden_file`, `--cnel_file`, `--lday_file`, `--levening_file`,
`--lnight_file`, `--laeq_file`, `--lamax_file`, `--nax_file`, `--ta_file`, `--csel_file`,
`--coverage_file`) is laid out as `-o` is, and is described with its
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
  its share, and its stored μ, γ and C<sub>T</sub>; last, the shares from spheres the database
  marks as extended to another flight path angle and as a synthesized hover (`condition_origin`;
  zero on a database that does not record it). Each state counts as the SEL counts it (inside
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
- `--phase_interval`: `_phase_lower` and `_phase_upper` companions beside the output
  ([Buildings around the microphones](#buildings-around-the-microphones)), laid out as it is.
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
- every input (the operations file, trajectories, databases, vehicle, observer, atmosphere, terrain and table files)
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
named for their track. Plotting is skipped. Each distinct database is loaded once, however its
trajectories are ordered, and two spellings of one path are one database.

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
`--start_utc` (a static operation's only clock), shifted by `--utc_offset`, and a track without one is refused unless an operations file gives period counts ([Period counts](#an-operations-file)). By default each operation contributes its own
penalized 10 dB-down window; `--full_history` integrates its whole history, as AAM and NORAH2 do,
which reads a few tenths of a dB higher ([FIDELITY](FIDELITY.md#day-night-window-choice)). Use it
for a DNL meant to sit beside their contours.

### An operations file

`--operations FILE` describes a day's traffic as a table: each row is a track and how many times
an average day flies it, and may name its own database, vehicle, clock, lateral offset and air. It replaces `-t`
(the two are refused together), and is the way to model a fleet mix or a schedule without a
command line per operation.

```text
# Runway 21 departures and arrivals, an average annual day
name,track,database,vehicle,count,local_time,wind_speed,wind_direction
dep21,../tracks/dep21.csv,../databases/B407_fixed_load_factor.nod,../databases/B407.json,6,10:00,8kt,200
dep21_night,../tracks/dep21.csv,../databases/B407_fixed_load_factor.nod,../databases/B407.json,0.5,23:30,3kt,180
arr03,../tracks/arr03.csv,../databases/AS350_fixed_load_factor.nod,../databases/AS350.json,4,14:00,8kt,200
```

The file is plain comma-separated text. The first line that is not blank or a `#` comment is the
header; every row after it must have as many fields, and fields are trimmed. Quoting is not
supported, so a field cannot hold a comma or a quote. Column names match whatever their case and
punctuation (`Local Time`, `local_time` and `LOCAL_TIME` are one column). A unit in parentheses
or brackets is checked where a column has one, and nothing is converted: a unit on a column that has none,
such as `count (per year)`, is refused, and `wind_direction` may say only `(deg)`.

| Column | Meaning | When the cell or column is empty |
| --- | --- | --- |
| `track` | The trajectory file: an AC track, a WGS84 track or a static operation. Required. | Refused. |
| `name` | The row's name: what its files are named after, and the `-m` case name. Letters, digits, `_`, `.` and `-`, not starting with a dot. | The track's file name, less its extension. |
| `count` | Operations per average day. Any positive number, fractions included. | 1. |
| `day`, `evening`, `night` | Operations per average day in each period, in place of `count` and `local_time` (below). Each zero or more, and a row's three add up to more than zero. | 0. |
| `local_time` | `HH:MM` or `HH:MM:SS` (00:00 to 23:59:59), local time by `--utc_offset`, at the track's first state. It replaces any clock the track carries. | The track's own clock, or `--start_utc`. |
| `database`, `vehicle` | As `-d` and `-v`. | The one `-d` or `-v` given; refused without it. |
| `output` | This row's SEL file, as `-o`. | Named after the row (below). |
| `atmosphere` | As `--atmosphere`. | See the air, below. |
| `wind_speed`, `wind_direction`, `wind_height` | Joined as `--wind`'s `speed,direction[,height]`: the speed names its unit (`8kt`), the direction is a compass bearing, the height is in `--units`. The first two go together. | See the air, below. |
| `humidity` (%), `temperature` (K), `pressure` (kPa) | As the options of those names. | See the air, below. |
| `offset` | A lateral subtrack: how far to the left of the direction of travel the row flies its track, in `--units` (negative is to the right). | 0, the track itself. |

**Paths** in the file are taken from the operations file's own folder, so a study folder moves as
a whole; absolute paths are used as given. Each resolved path is in the log, and a track,
database, vehicle or atmosphere file that does not exist is refused before anything is computed,
naming the row and the path.

**Per-trajectory options.** `-d`, `-v`, `--atmosphere`, `--wind`, `--humidity`, `--temperature`,
`--pressure` and `--start_utc` may each be given once with `--operations`, as the value for the
rows that leave that column empty; given twice they are refused, since per-operation values go
in the file. The demo's S-76D database and vehicle never fill a row.

**The air.** An air column that is present must be filled in every row, unless its option was
given once, which then fills the empty cells. The run's own checks then apply to the rows as
they would to the command line: rows mixing an atmosphere file with a wind, for example, are
refused.

**Counts.** A row's count multiplies its contribution to every run total after its exposure is
formed: DNL, Lden, CNEL and the period levels add count times the operation's energy, NAx adds
the count when the operation's maximum reaches the threshold, and the time above adds count
times the operation's seconds ([THEORY §10.3](THEORY.md#103-maximum-level-time-above-and-number-above),
[§10.4](THEORY.md#104-daynight-metrics)). LAmax, the largest level, is not raised by a count.
A count of 3 gives what three `-t` of the same track give, and a count of 10 adds 10 dB to the
day metrics. The coverage and provenance reports weight each row the same way.

**The clock.** `local_time` dates the row's first state, and every later state by the track's
own time; the day-night penalties then fall where they would for a flight at that hour, across
22:00 or midnight too. A row without `local_time` keeps its track's clock (an AC track's first
column, a WGS84 track's `utcsec`), or `--start_utc`, and a track with neither is refused by the
day metrics. Give each operation its own row and `local_time`: a row's count flies it that many
times at the same hour.

**Period counts.** A file may instead give each row's operations as `day` (07:00 to 19:00),
`evening` (19:00 to 23:00) and `night` (23:00 to 07:00) counts, as an AEDT or INM study's fleet
tables do. They take the clock's place: each operation is heard as its plain SEL, and each
metric weights a period's count as it weights an operation flown then. A night operation counts
ten times in DNL, Lden and CNEL. An evening operation counts three times in CNEL, about 3.16
times (+5 dB) in Lden, and once in DNL, which has no evening: DNL counts it as day, as AEDT and
INM do. DNL and CNEL weight every evening operation as one flown 19:00 to 22:00, though by the
clock 22:00 to 23:00 is their night. `--lday_file`, `--levening_file` and `--lnight_file` each
take only their own period's count. NAx, the time above, the coverage and the provenance take the
row's three counts added up. A track's own clock is not read for the
day metrics (the log says so), and a track with no clock runs. A file has either period columns
or `count` and `local_time`, never both; `--laeq_file`, whose window is a clock window, is
refused with period counts. Use `local_time` where the hour matters within a period, or when a
flight straddles a period's edge: period counts weight the whole operation by one period.

**Subtracks.** A route is rarely flown on its centerline, and a study spreads its operations over
subtracks either side of it: one row per subtrack, each with its share of the count and its
`offset`. The offset is measured square to the track's horizontal direction of travel, so on a
turn the subtrack flies a tighter turn on the inside and a wider one on the outside.

- A WGS84 track is displaced before it is smoothed, so its speed, load factor and thrust are the
  subtrack's own: on the inside of a turn it flies slower and with less bank. Velocity columns
  in the file are the centerline's, so with an offset they only set the direction, and the
  velocity comes from the displaced path (the log says so).
- An AC track is moved sideways state by state, keeping the speeds and accelerations it
  recorded, which is right on a straight leg and approximate in a turn.
- Where the track is slower than 3 kt (a hover's drift) the direction is held from the nearest
  faster point. If the directions before and after a slow interval differ, the offset is
  refused: switching sides would make the displaced path jump. Give that subtrack's path
  explicitly or use zero offset. A track that is never faster than 3 kt is refused, and so is
  an offset on a static operation: move the hold's x and y instead.
- An offset toward the inside of a turn must be under 90% of the turn's radius; one that is not
  would fold the path back on itself, and is refused naming the time.
- `--trim` applies to the displaced track, and the outputs, the `-S` state included, are the
  subtrack's.

```text
name,track,count,offset (ft)
dep21_c,../tracks/dep21.csv,6,0
dep21_l,../tracks/dep21.csv,2,500
dep21_r,../tracks/dep21.csv,2,-500
```

**Outputs.** Each row's SEL file, and its `-S`, `--csel_file`, `--epnl_file`, `--band_frames` and
`--write_atmosphere` files, are one operation's, whatever its count. With several rows each is
`<name>_LAE.csv` (and `<name>_state.csv` and so on) in the folder of the name given to the option
(`integrated/` for the default `-o`); a row's `output` cell names its own. A one-row file keeps
the `-o` it was given. The `-i` and `--band_history` files carry the row's name, as they carry
the track's for several `-t`. These options take one name with `--operations`. Two rows may not
share a name, compared ignoring case: give such rows a `name` column.

**Refusals.** An unknown or repeated column, a file with no `track` column or no rows, a row with
the wrong number of fields, a quote, `count` or `local_time` beside period columns, period
counts that are negative or add up to nothing, `--laeq_file` with period counts, a count that is
not a positive number (`0`, `-2`, `nan`,
`6.2 ops`; comment a row out with `#` to leave it out), a `local_time` that is not a clock time
(`25:00`, `7pm`), a `start_utc` column (use `local_time`), a temperature, pressure or humidity
column declaring another unit, an `offset` or `wind_height` column declaring a length unit other
than `--units`, a wind speed without its unit, and a name outside the allowed characters. A
row is skipped, and the run exits 1, for an offset on a static operation, on a track never faster
than 3 kt, or toward the inside of a turn by 90% of its radius or more.

**Cost.** One run per row, as for a `-t` per row. Each distinct database is loaded once,
however its rows are ordered. Rows run in file order.

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

The build also makes tools for building tracks (`niceops_procedure`), profiles
(`niceops_build_atmosphere`) and tables (`niceops_propagation_table`, `niceops_urban_table`), for
the field among buildings (`niceops_point_field`, and `niceops_fdtd_field` as its full-wave
reference), and for checking the propagation models (`niceops_pe_compare`, `niceops_benchmark`).
They share the `--units` convention (feet by default) and the profile format.

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
| `speed` | `to_speed`, optional `rate` (default `limits.acceleration`) and `angle` (default level) |
| `turn` | `bank` (signed, positive right, under 60°) and `turn`, the course change |
| `hover` | `duration_s`; it must start at rest |

- **The start** is a position, `alt` on the track's datum (what `niceops` takes
  `--ground_elevation` off), a ground speed and a compass course. An optional `ground_alt` (default
  0) is the ground on that datum: a segment that goes below it is refused.
- **Transitions are flown**: the path angle changes at a limited pull (`limits.pull_g`, default
  0.1 g), the bank at a limited roll rate (`limits.roll_rate_deg_per_s`, default 10), and speed at
  the segment's rate (`limits.acceleration_mps2` or `_kt_per_s` where a segment gives none,
  default 0.5 m/s²), so a join puts no spike in the load factor. A climb overshoots its altitude
  by the leveling-off. A turn rolls in, holds and rolls out within itself, ending level-winged on
  the asked course to within a degree; one too short to reach full bank banks only as far as it
  needs. A target altitude already passed, and any key the format does not know, are refused.
- **The end** is flown on to the next written step, so the last state is written.
- **Calm air**: the heading is the course and the ground speed the airspeed. `step_s` (default
  0.2) is the spacing of the written states.
- **`--check DATABASE --vehicle VEHICLE`** reads the written track back as niceops reads it and
  reports the states the database would clamp on each axis (advance ratio, path angle, C<sub>T</sub>)
  or hold at `--max_load_factor` (default 2 g, as in niceops), with the stretches of time they
  span, so a procedure can be kept inside its database before any acoustics run. `--no-qsam`
  matches niceops's `--no-qsam`. A vehicle file that disagrees with the database is refused as
  niceops refuses it, and `--accept_vehicle_mismatch` maps the states anyway. The procedure is calm, so the database's speed reference makes
  no difference. How much of the footprint those states carry is what `--provenance_file`
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
| `--gamma_t`, `--outer_scale`, `--realizations`, `--seed`, `--gamma_profile` | `pe`: turbulence by phase screens, for reference only. Defaults still air, 10 m, 1 realization, seed 1. `--gamma_profile` is a CSV of z (in `--units`) and γ. `ground_model`: `--gamma_t` alone, the HARMONOISE coherence of the direct and reflected paths, as `niceops --turbulence_model harmonoise --gamma_t`. |
| `--turbulence`, `--roughness` | `ground_model`: Chessell's coherence (rad·s/√m) and the ground's RMS roughness (in `--units`), as for `niceops`. Default 0. Not with `--gamma_t`, and not for `pe`. |
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

### `niceops_build_atmosphere`

Builds an `--atmosphere` profile from measured columns, and fits a site's roughness length to
wind profiles. Use `build` for anything measured or modeled in height: a met tower, a sodar or
LIDAR, a temperature string or tethered balloon, a radiosonde, or a weather model's column already
extracted to CSV. Use `fit_z0` for a site with a LIDAR or tower whose roughness length nobody
measured.

```sh
niceops_build_atmosphere build --column lidar.csv --column string.csv --missing 9999 \
    --humidity 62.5 --z0 0.1 --top 5000 -o atmospheres/run.csv
niceops --atmosphere atmospheres/run.csv --pressure 101.2 --refraction stratified ...

niceops_build_atmosphere fit_z0 runs/*.csv --missing 9999 -o site_z0.json
niceops_build_atmosphere build ... --z0_from site_z0.json -o atmospheres/run.csv
```

**Inputs to `build`.** Each `--column FILE` gives every variable it carries, under the column
names of [Stratified atmosphere](#stratified-atmosphere). Any variable may be missing, and the
humidity may also be a dew point (`Td_C`, `Td_K`) or a specific humidity (`q_kgkg`, `q_gkg`, with
the pressure in the same rows). A variable two files carry is refused; `--wind_from`,
`--temperature_from`, `--humidity_from` and `--pressure_from` name the file to take it from.
`--temperature`, `--humidity`, `--pressure` and `--wind SPEED,DIRECTION[,HEIGHT]` (as niceops's
`--wind`) are a surface station: its air is the ground level's, and its wind a level at the
anemometer. Lengths are in `--units`, feet by default, on the command line and in the files; an
altitude column has `--ground_elevation`, the observers' elevation, taken off it.

- **`--missing V ...`** lists cell values that are fill, not data (9999, 99999, a LIDAR's
  1.944e4). They are dropped for their variable at their level, and counted. Any other
  implausible value is refused with its file and row; `--drop_implausible` drops those too. A
  wind column left with no value is refused, since it would build still air; a temperature,
  humidity or pressure column left with none is warned of and ignored.
- **Instrument cleaning stays with the caller**: a LIDAR's 180° direction ambiguity, a balloon
  that was not aloft, or a GRIB2 file. The tool takes clean columns.

**The profile.**

- **The ground level** takes the station's value where one is given. Otherwise it takes the
  measurements interpolated to the ground; levels below it (an instrument standing lower than
  the microphones) count only for that. Otherwise the lowest value is held down, with a warning.
  The wind at the ground is calm.
- **Below the lowest wind level**, `--fill log` (the default) carries that level's wind down by
  the `--wind` log law, which makes the roughness length the only thing the fill adds. `--fill
  hold` holds the level's wind down; `--fill none` takes it linearly to calm at the ground. The
  fill's levels start at ten roughness lengths (or the lowest wind level, if lower), so no level
  above the ground is calm. A temperature, humidity or pressure measured lower is a level too,
  except at or below the roughness length, where the law is calm: there it sets the values
  around it without being a level, with a warning.
- **`--z0` or `--z0_from` is required.** Nothing measures the roughness length, and the wind
  near the ground depends on it, so a profile never carries a default. A weather model's z0 is a
  grid cell's, not the site's.
- **Between measurements** temperature and humidity are linear in height, the wind linear in its
  components, and pressure in its logarithm. **Above** each variable's highest measurement it is
  held; `--above lapse` lapses the temperature dry-adiabatically instead, and `--humidity_aloft
  held_concentration` keeps the water vapor and lets the relative humidity follow the cooling air.
  `--surface_layer_depth` carries the law up from a single wind level, as `--wind` does; with
  several levels it is refused. `--stability` applies only to a temperature measured at the ground
  alone, with `--above lapse`.
- **Levels.** By default they are geometric, `--levels_per_decade` (16) a decade, through every
  measured level. `--levels 0,1,3,10,20` instead gives the levels to use, to which every
  measured level, the ground and the top are added. Give the levels a profile was run with when
  reproducing it: the stratified ray search follows the layers.
- **`--top` is required**: a run refuses a track above the profile.
- **Pressure** is written when measured, or with `--write_pressure` or `held_concentration` from
  `--pressure`, hydrostatically. Otherwise the file has none and the run integrates its own.

**What it refuses**, naming the cause: a calm wind measured above the ground (`--calm_levels
drop` drops it), a roughness length at or above the lowest wind level, a `--levels` height at or
below the roughness length under the log fill, pressure rising with
height, a relative humidity above 100 (`--clip_humidity` clips it), two levels of one variable at
one height, a variable two files offer, a speed without its unit, and any output `--atmosphere`
would refuse. The written file is read back through `--atmosphere`'s reader before the tool exits.

**Provenance.** Beside the profile, `OUT.provenance.json` (or `--provenance`) records the build,
the command line, each input's digest, the rules, the roughness length and where it came from,
which heights each variable was measured over and which stretches of wind were filled, the counts
of dropped values, and every warning. The CSV itself stays plain.

**`fit_z0`** takes one column file per run, named for the run, and fits one roughness length to
every run's gates between `--gates LOW,HIGH` (default 46,197 ft), with the log law's slope free
per run ([THEORY §7](THEORY.md#7-the-stratified-moving-atmosphere)). A run enters only if it is
near neutral: `--max_richardson` (0.05) on the bulk Richardson number from the temperatures in the
same file, at least `--min_speed` (3mps) at the highest gate, and at most `--max_veer` (20°) of
turning. A run without temperatures at both gate heights is left out unless `--assume_neutral`.
The output JSON holds z0 with its 95% bootstrap interval over runs (`--draws`, `--seed`), the runs
used and each one left out with its reason; `--table` writes a one-row CSV. Fewer than
`--min_runs` (5) runs are refused. The fit is ill-conditioned (a 1% error in the speed ratio
moves z0 by a factor of about 1.4), so pool runs, and quote the interval with z0. A low gate
biased by its probe volume would bias z0.

Both commands take milliseconds to a second.

### `niceops_urban_table`

Builds an `--urban_table` ([Urban receivers](#urban-receivers)), either for microphones in ideal
street canyons (`--canyons`) or for microphones among LoD1 buildings (`--buildings` with
`--receivers`, below). Each file is a CSV read by column name. Lengths are in `--units`; a column
may name its unit (`x_ft`, `height_m`), which must then be `--units`.

```sh
niceops_urban_table --canyons streets.csv --database databases/B407_fixed_load_factor.nod \
    --frame_bearing 210 -o tables/streets.nc
```

| Column | Meaning |
| --- | --- |
| `x`, `y`, `z`, `height` | The observer and its microphone height, as the run's observer file gives them. An `id` column is ignored. |
| `axis` | Compass bearing of the street, degrees, either way along it. |
| `width` | Between the walls; 0 for an open microphone, when the rest of the row is not read. |
| `offset` | The microphone's distance from wall A, the wall on the left looking along `axis`. |
| `height_a`, `height_b` | The walls' tops above the ground. |
| `absorption_a`, `absorption_b` | Optional: each wall's energy absorption coefficient, every band. |

| Option | Meaning |
| --- | --- |
| `--canyons` | The CSV above. Either this or `--buildings`. |
| `--buildings`, `--receivers` | The district and its microphones, below. Each needs the other. |
| `--roofs` | With a CSV `--buildings`: pitched roofs, as `niceops --roofs` takes them. |
| `--no_roofs` | With `--buildings`: the roofs reflect nothing, and a microphone on or over one is refused. For comparison. |
| `--database`, `--bands` | Band centers from a database, or given. |
| `--azimuths` | Azimuths evenly over 360°. Default 72. |
| `--elevations` | Elevations, degrees. Default 0 to 30 by 1.5, then 35 to 90 by 5. |
| `--frame_bearing` | The compass bearing of the run's local +x axis. Default 90 (east). |
| `--absorption` | The walls' absorption where a row gives none. Default 0, rigid. |
| `--floor_absorption` | The street's absorption, for the rays its floor reflects. Default 0, rigid pavement. |
| `--facade_length` | Each wall's length for the run's reflections, centered on the microphone. Default 10 km, effectively infinite. |
| `--sound_speed` | For the diffracted path's wavelengths, m/s. Default 343.2. |
| `--no_diffraction` | Leave out the roof edges' diffraction and sum the specular rays in energy: the field then jumps wherever a ray is cut off. For comparison only. |
| `--wave_below` | Solve each canyon with walls of equal height by the wave equation and use its table for the bands centered at or below this frequency, Hz. Default 0, none. |
| `--wave_blend` | Bands above `--wave_below` over which the wave table fades into the geometric one. Default 2. |
| `--wave_ppw` | The wave solver's grid points per wavelength. Default 15. |
| `--gpu`, `--no-gpu` | March the wave equation on the GPU (Metal, on a Mac), as below: the default wherever there is one and `--wave_below` is given. `--gpu` asks for it and is refused without one, or without `--wave_below`; `--no-gpu` keeps the CPU. |
| `-j`, `--threads` | Canyons solved at once, or a receiver's directions shared out. Default the environment's `NICEOPS_THREADS`, else one per hardware thread. |
| `--units`, `-o` | `ft` (default) or `m`; the table to write (SI). |

**The wave band.** Below a few hundred hertz a street is a few wavelengths wide, and its
resonances and diffraction are not rays. `--wave_below` solves each canyon whose walls are of
equal height by the wave equation in its cross-section, once per distinct canyon, and takes that
solution for the low bands, which already hold the first-order reflections, so the run's own
reflections fade out of those bands as the wave table fades in. Canyons with walls of unequal
height keep the geometric table, with a warning naming the row. Give the walls a little
absorption: a rigid canyon rings for seconds, and a warning says when a canyon had not decayed
when its solution stopped. Each doubling of the frequency costs about eight times as much: with
the default blend, `--wave_below 400` solves each canyon to 707 Hz, about a minute at 15 points per
wavelength, and canyons are solved several at a time. On a machine with a GPU the wave equation is
marched there instead, in single precision, two canyons at a time (`--no-gpu` keeps the CPU);
[FIDELITY](FIDELITY.md#canyon-wave-solver) has how much faster that is and how closely its tables
agree with the CPU's.

**Buildings.** `--buildings` gives the buildings as footprints extruded to flat roofs, either as a
CSV with one row per footprint vertex (below) or as GeoJSON; `--receivers` gives the microphones,
with the run's observer-file columns (`x`, `y`, `z`, `height`; build from the run's own observer
file). The buildings stand on the ground at each receiver's `z`. A microphone flush on a wall (within
a micrometer of it) takes the wall's reflection there, as one on the floor takes the floor's: a ray
from in front arrives doubled, `1 + sqrt(rho)` in amplitude, and the table then holds the
first-order reflections the run would otherwise evaluate itself, which take it too. A microphone
over a roof hears it as `niceops --buildings` does, specularly and off its edges, and one flush on
a roof or a slope takes its reflection there as one on a wall does; `--roofs` gives a CSV
district's pitched roofs as `niceops --roofs` takes them, and `--no_roofs` leaves the roofs'
reflections out and refuses a microphone on or over one, for comparison. One on a building's
edge, or inside one, is refused.

```sh
niceops_urban_table --buildings district.csv --receivers observers.csv \
    --database databases/B407_fixed_load_factor.nod --frame_bearing 210 -o tables/district.nc
```

| Column of `--buildings` | Meaning |
| --- | --- |
| `building` | An id. A building's rows come together, its vertices in order around it, either way. |
| `x`, `y` | The vertex, in the run's local frame. |
| `top` | The roof above the ground, the same on each of the building's rows. |
| `absorption` | Optional: its walls' energy absorption coefficient, every band, the same on each of its rows. Default `--absorption`. |
| `scattering` | Optional: the share of what its walls reflect that they scatter diffusely, every band, the same on each of its rows. Default `--facade_scattering`. |

**GeoJSON** (a `.geojson` or `.json` file, as an OpenStreetMap export or a city's LoD1 footprints
give it): each Polygon or MultiPolygon feature is a footprint in WGS84 longitude and latitude, its
inner rings courtyards. It is projected into the run's frame, so give the run's own `--origin` and
`--frame_bearing` -- for the shipped AC tracks, the origin their files were made with and 210.
Heights are meters, from the properties `height` or `building:height` (a number, or text such as
`12.5 m` or `40 ft`), else `building:levels` times `--level_height`; buildings with none are
refused unless `--default_height` is given. A feature raised off the ground (`min_height` or
`building:min_level` above zero: a bridge, an overhang) is left out and counted. Optional numeric
properties `absorption` and `scattering` set the walls'. Buildings that share a side exactly, as
terraced houses in OpenStreetMap do, stand against each other there: the shared wall is inside the
taller one, which keeps only what stands above the lower roof.

| Option | Meaning |
| --- | --- |
| `--origin` | With GeoJSON: the run's WGS84 origin, `latitude longitude`. Required. |
| `--level_height` | With GeoJSON: meters per story, for buildings that give only `building:levels`. Default 3. |
| `--default_height` | With GeoJSON: the height, meters, of buildings that give none. Default: refuse them. |
| `--max_order` | Facade reflections on a ray from the source, at most. Default 64: bouncing up a street is cheap, and grazing directions need it. |
| `--max_diffracted_order` | Facade reflections on a ray after an edge diffracts it, at most. Default 6. |
| `--radius` | The edges, and the facades a diffracted ray reflects in, within this of a microphone. Default 150 m. |
| `--orientation_tolerance` | Walls whose orientations agree within this, degrees, reflect alike: a ray then bounces up a street a little out of parallel as up the canyon it nearly is, and a district of many orientations costs far less. Default 1; 3 for an irregular old town. The run's own first-order reflections stay exact. |
| `--facade_radius` | The facades within this of a microphone, facing it, are written for the run's first-order reflections. Default 100 m. |
| `--facade_scattering` | The share of what a wall reflects that it scatters diffusely, where a building gives none. Default 0: specular walls. |
| `--patch_size` | With scattering walls: the walls are cut into patches no larger than this for the diffuse field. Default 4 m. |
| `--diffuse_radius` | With scattering walls: the walls within this of a microphone scatter. Default 100 m. |

**Scattering walls.** Real facades -- windows, balconies, ledges -- scatter part of what they
reflect, and that diffuse energy is what reaches a microphone the direct sound and its
reflections do not. With a scattering coefficient s the walls reflect rho (1 - s) specularly,
in the table's rays and in the run's own first-order reflections, and rho s diffusely
([THEORY §6.5](THEORY.md#65-urban-receivers)). Measured facade data are scarce: ISO 17497-1
gives a random-incidence scattering coefficient, and facades measured at steep incidence hardly
exist, so s is a modeling choice rather than a value to look up. It matters for a particular
microphone more than for a district: across synthetic streets a modest s moves the median
microphone's SEL little but many microphones' by more, most in deep streets
([FIDELITY](FIDELITY.md#buildings)). The default leaves the walls specular, and tables built with
it are unchanged. The diffuse field adds several seconds to a receiver's build on one thread.

A microphone inside a footprint is refused, as is `--wave_below`. A table for a street
microphone among a hundred blocks takes a few seconds per receiver at the defaults; the time
grows about 1.7-fold with each further order of diffracted reflection
([FIDELITY](FIDELITY.md#buildings)), and a warning says if a receiver's search was cut short.

### `niceops_point_field`

The field of a point source among LoD1 buildings, at microphones: for a source near the buildings
-- a loudspeaker in a street, a drone among the roofs -- where an `--urban_table`'s plane wave will
not do. It traces the source's own rays, reflected in the facades, the roofs and the floor, diffracted at
the roof edges and the corners, once, at two in turn (with a facade between them or not) or at three
roof edges in turn, and writes each band's level re the
source's level at 1 m, so a measured source spectrum at 1 m adds to it directly.

```sh
niceops_point_field --buildings district.csv --source 12.5 30 2 --microphones mics.csv \
    --bands 160 200 250 315 400 500 630 800 1000 --humidity 70 --units m -o levels.csv
```

`--buildings` is `niceops_urban_table`'s, CSV or GeoJSON (with `--origin` and `--frame_bearing`,
and `--level_height` and `--default_height` as there), and the buildings stand on a flat floor.
`--roofs` gives a CSV district pitched roofs, as [`niceops --roofs`](#pitched-roofs) takes them, and a ray over four roof edges
or more -- a mansard heard from low beside it, several rows of roofs -- is traced up to
`--max_chain` of them. Roofs reflect as the floor does, with their building's
absorption: a level roof replaces the floor, and reaches an end or an edge above it. A slope
may reflect with walls and the floor on either side, or next to the edge of a once-diffracted
ray. A twice-diffracted ray may reflect in a level roof or a slope between either end, or both, and its edge.
Over a row of two or more flat-roofed blocks with a street between, or touching with a roof lower
than both its neighbors (a dip; touching blocks of one height count as one roof), crossed near grazing, where
the cascades over three and four roof edges do not settle, the rays over the row give way to its
two-dimensional field from a parabolic equation over its cross-section ([THEORY §6.5](THEORY.md),
"Over a row of roofs"), its floor bounces at the ends taking the floor as any ray's do. It is taken
up to a size in wavelengths, k times the row's width across from 3000, gone by 6000: 1-2 kHz over a
terrace 160 m across, above 5.5 kHz over one street.
`--no_slope_facades` omits paths combining slopes and walls for comparison,
`--no_twice_roofs` the twice-diffracted rays off a roof, `--no_roof_lighting` the scattering
patches' lighting off a roof, and `--no_profile_pe` the row's profile field;
`--tilted_step_split` adds the jump split across a step's line, which is off by default;
`--unending_edges` omits UTD's finite-edge correction. `--microphones` is a CSV with `x`, `y` and `height`
columns, the height above the floor: a microphone (or a source) at a flat roof's height over it, or
on a wall, lies on that surface and takes its reflection there, as one on the floor takes the
floor's, and one on a pitched roof's slope the slope's; one on a building's edge is refused. Each output row is a microphone's `x`, `y` and `height`, then one level per band, dB, the
rays summed coherently and averaged over the band, each absorbed by the air along its length
(ISO 9613-1, as the spectral path's). A free-field level at the microphone would be −20 lg r for
its distance r in meters, less that absorption.

The absorption depends strongly on the humidity, so a run without `--humidity` is refused rather
than given one nobody chose; `--no_air_absorption` leaves the air out instead. The temperature
and pressure default to the standard's reference, 293.15 K and 101.325 kPa. Over a few tens of
meters the absorption is small below a few kilohertz; over hundreds it takes most of the top
bands.

For a source that moves -- a drone's states along its track -- `--sources` takes a CSV of
positions with the same three columns, and each output row then starts with the source's row in
it, the sources one after another over every microphone. The levels are those of one run per
source, to the bit, and the run is quicker than those runs. The output is written as the run goes,
a block of sources at a time, so a run that fails partway exits with an error and leaves the blocks
it finished.

On a machine with a GPU (Metal, on a Mac) each ray's amplitude at every band is formed and the bands
summed there, while the CPU traces the next pairs' rays; `--no-gpu` keeps the CPU, and
`--edge_source` and `--fade` always do. The rays are the same, and the levels agree closely with
the CPU's ([FIDELITY](FIDELITY.md#a-point-source-among-buildings) has the comparison): in the low
bands, where every pair of rays interferes, the GPU averages over the band exactly,
while the CPU leaves out the pairs that the averaging has all but cancelled. The GPU works in
single precision. `--gpu` asks for it, and is refused by a build without a GPU backend, a
machine without a GPU, or with `--edge_source` or `--fade`.

`--edge_source` takes the diffracted rays' fields from the edge-source integral ([THEORY
§6.5](THEORY.md#65-urban-receivers)) instead of UTD's coefficients. The integral runs over each
edge's whole length, so a building's corners and the ends of its roof edges diffract too, and a ray
fades rather than stops as its diffraction point leaves its edge. Every once-diffracted ray takes it
but one off a roof, which takes UTD's coefficient. A twice-diffracted ray takes it where the cascade
of two coefficients cannot be trusted -- near either edge's shadow or reflection boundary (but not
down a step, where the uniform form holds the transition better), near an edge's end, over a shared
face -- and the cascade elsewhere, the two blended smoothly band by band; out of both transitions,
the cascade times each edge's end ratio. It is several times UTD's cost on a few microphones and
far more over a large district at low frequency, and CPU only: `--edge_source` is refused with
`--gpu` or `--fade`, and keeps a default run on the CPU. Its accuracy and
cost are in [FIDELITY](FIDELITY.md#a-point-source-among-buildings).

Walls with relief -- corrugated cladding, balconies, ribs, window bays -- send part of what they
reflect in every direction rather than specularly. A wall's scattering coefficient `s`, from a
building's `scattering` column or `--scattering` for those that give none, takes that share of its
reflection off the rays and returns it diffusely, by Lambert's law: the buildings with a wall
within `--diffuse_radius` of each microphone are cut into patches of at most `--patch_size`, walls
and roofs, each lit by the source's own reflections and the floor's, and the energy they scatter
among themselves and to the microphone adds to the rays' in each band
([THEORY §6.5](THEORY.md#65-urban-receivers)). A roof scatters with its building's coefficient, as
it absorbs with it. A shadowed microphone gains most: the patches it can see may be lit when it is
not. A microphone above the roofs gains from theirs, which a street's does not see. Where the source and the microphone see each other along a
street, a relief mostly redirects the reflections it changes, which Lambert's law does not hold,
and the coefficient does little good there
([FIDELITY](FIDELITY.md#scattering-walls-point-source)). The coefficient depends strongly on frequency and on the
relief's depth against the wavelength: `tools/profile_scattering.py` computes a relief's from its
profile, band by band, and `--scattering_spectrum` takes that file. Each band's coefficient is then
the spectrum's, interpolated in log frequency, times the building's own value, which defaults to 1:
give 1 to the buildings with that relief and 0 to plain ones. A run whose bands the spectrum does
not cover is refused. The patches cost time and memory as the square of their number, so keep
`--patch_size` near a tenth of the street's width. Diffraction lights no patch.

| Option | Meaning |
| --- | --- |
| `--buildings` | The district. |
| `--source` | The source's `x`, `y` and height above the floor. Refused inside a building. |
| `--sources` | Instead of `--source`, a CSV of source positions, `x`, `y` and `height` by name. |
| `--microphones` | The CSV above. A microphone inside a building is refused. |
| `--bands`, `--database` | Band centers, given or from a database. |
| `--absorption` | The walls' energy absorption where a building gives none. Default 0, rigid. |
| `--floor_absorption` | The floor's. Default 0, rigid. |
| `--relief`, `--relief_scattering`, `--grating_reflections` | Walls with a periodic relief, as `niceops --facade_relief` takes them (above): the reliefs' grating-order tables, where what they scatter beyond their specular order goes (`hybrid`, `orders` or `lambert`), and a grating ray's other reflections. Not with `--scattering_spectrum`. |
| `--flow_resistance`, `--impedance_model` | A floor of finite impedance instead, as `niceops --buildings --ground_effect` takes it: each bounce takes the ground's Q. The flow resistance in kPa s/m²; the model `delany_bazley` (default) or `miki`, the two that need nothing else. Not with `--floor_absorption`. |
| `--scattering` | The walls' and roofs' scattering coefficient where a building gives none, 0 to 1. Default 0, specular; 1 with `--scattering_spectrum`. |
| `--scattering_spectrum` | A CSV of `band_hz` and `s`: the coefficient band by band, each building's weighting it, as above. |
| `--patch_size` | With scattering walls, the diffuse field's patches' size at most. Default 4 m. |
| `--diffuse_radius` | With scattering walls, a building with a wall within this horizontal distance of a microphone scatters, walls and roof; those beyond reflect specularly. Default `--radius`. |
| `--max_reflections` | Facade reflections on a ray that diffracts nowhere, and on a path that lights a patch. Default 4. |
| `--max_edge_reflections` | Facade reflections on each side of a single diffraction. Default 2. |
| `--max_twice_reflections` | Facade reflections on each side of a double diffraction. Default 2. |
| `--max_between_reflections` | Reflections, in a facade or the floor, between the two edges of a double diffraction, 0 or 1. Default 1. |
| `--triple_diffraction` | Rays diffracted at three edges: `roofs` (default; three roof edges, the rows a street's sound crosses), `all` (any three, every chain of corners round a block too: far more rays for little), or `none`. |
| `--max_thrice_reflections` | Facade reflections on each side of a triple diffraction. Default 0. |
| `--roofs` | Pitched roofs for a CSV district, one row per plane, as `niceops --roofs` takes them. |
| `--max_chain` | Rays diffracted at more than three roof edges, up to this many: the corners of the taut string over the buildings' cross-section between the source and the microphone, the floor at either end or not. Default 8, at most 14; 3 for none. |
| `--radius` | Facades and edges within this horizontal distance of the source or a microphone. Default 150 m. |
| `--no_diffraction`, `--no_double_diffraction` | Leave out every diffracted ray, or those at two edges or more. |
| `--no_roofs` | Roofs do not reflect, for comparison. |
| `--no_slope_facades`, `--no_twice_roofs`, `--no_roof_lighting` | Leave out, for comparison: the rays off both a slope and a wall; a level roof's or slope's reflection on a twice-diffracted ray; with scattering, the patches' lighting off a level roof or a slope. |
| `--no_profile_pe` | Over a row of roofs crossed near grazing, the cascades' rays instead of the row's profile field, as above. |
| `--tilted_step_split` | Split the field's jump across a step's line, as `niceops --tilted_step_split`. Off by default. |
| `--unending_edges` | UTD's edges without their ends, each taken as its whole infinite line, for comparison. |
| `--no_prune` | Solve and check every diffracted path rather than ruling most out first: the same rays, slower. For checking. |
| `--fade` | Where a boundary leaves the field to an order the run does not trace -- a third reflection about an edge, a reflection between two edges -- make that ray from the one it partners, so the level does not step there (specular and once-diffracted rays). Not with `--edge_source` or `--gpu`. |
| `--fade_reach` | With `--fade`, the largest detour by an edge a partner is made for. Default: where the lowest band's Fresnel parameter reaches 3, at most 0.5 m. |
| `--fade_incident` | With `--fade`, the made partner's coefficient from the incident wave's terms alone. For comparison. |
| `--cascade` | Over a shared face (a roof from edge to edge, a wall from corner to corner) or a step (a roof's edge down to a lower roof's), the cascade of two edge coefficients instead of the uniform form. For comparison only. |
| `--exact_shared_faces` | Over a shared face or a step, the uniform form's integral at every wavenumber rather than only where the cascade cannot stand in for it, as `niceops --buildings --exact_shared_faces` does. The reference; much slower. |
| `--edge_source` | Diffraction by the edge-source integral, which keeps the edges' ends, as above. Not with `--gpu`. |
| `--origin`, `--frame_bearing`, `--level_height`, `--default_height` | With a GeoJSON district, as `niceops_urban_table`'s; `--frame_bearing` defaults to 90 (east). |
| `--orientation_tolerance` | As `niceops_urban_table`'s. Default 1 degree. |
| `--sound_speed` | m/s. Default 343.2. |
| `--humidity` | Relative humidity, percent, for the air's absorption. Required, unless `--no_air_absorption`. |
| `--temperature`, `--pressure` | For the absorption: K, default 293.15; kPa, default 101.325. |
| `--no_air_absorption` | Leave the air's absorption out. Not with the three above. |
| `--rays` | Also write every ray: its source's row (with `--sources`), its microphone, kind, reflections, length, its level at each band center (no band averaging; absorbed at the band center's coefficient), and its path from the source as `x y z;x y z;...`, for drawing. A leg off the floor is drawn unfolded below it, so only its `x` and `y` are the ray's. A relief's grating rays follow, of kind `grating order m`, each with its level in its own band and −300 in the others; a ray meeting a relief carries its specular share. A row's profile field is four rays of kind `profile field`, one for each of its parts, `off the floor at the source` or `at the microphone` naming the images, each drawn along its taut string over the roofs' corners. |
| `--gpu`, `--no-gpu` | Form the rays' amplitudes and sum the bands on the GPU, as above: the default wherever there is one. `--gpu` asks for it; `--no-gpu` keeps the CPU. |
| `-j`, `--threads` | Pairs of source and microphone evaluated at once. Default the environment's `NICEOPS_THREADS`, else one per hardware thread. |
| `--units`, `-o` | `ft` (default) or `m`, for every length in and out; the file to write. |

A twice-diffracted ray takes no roof bounce; with both the source and the microphone above a
shared face, its double diffraction is the cascade (the integral, with `--edge_source`). The field is
continuous across a boundary wherever the diffracted ray that takes over there is one of the
orders traced. Where it would need more reflections about its edge than `--max_edge_reflections`,
a third edge, or a reflection between two edges, the level steps over a few centimeters; a map
shows these as thin lines, and a microphone that sits on one is uncertain by the step. Accuracy
against an exact solution is in [FIDELITY](FIDELITY.md#a-point-source-among-buildings).

### `niceops_fdtd_field`

The full-wave field of a point source among rigid buildings, at microphones: the wave equation
marched on a grid in three dimensions ([THEORY §9.4](THEORY.md#94-a-full-wave-field-among-buildings)),
the reference `niceops_point_field` is checked against where no exact solution exists. It takes the
same `--buildings`, `--source`, `--microphones` and `--bands` and writes the same rows, each band's
level re the source's at 1 m, so the two subtract. It is meant for small scenes: a scale model, or a
street at low frequency.

```sh
niceops_fdtd_field --buildings case/buildings.csv --source -2 0 1.5 --microphones mics.csv \
    --bands 1600 2000 2500 3150 --units m --mirror --gpu -o levels.csv
```

The domain is the box around the buildings, the source and the microphones, `--margin` wider on
every side and above (default a wavelength at the lowest band's lower edge), then an absorbing
layer. Buildings cut by the box carry on through that layer. The differences are fourth order
(`--order 2` gives Yee's second order, kept for comparison), and the cell (`--cell`) defaults to
an eighth of a wavelength at the top band's upper edge (a twelfth at second order). Fourth order
is the more accurate of the two at fewer cells, over long paths most of all. The cost goes as the
inverse fourth power of the cell, and memory as the cube, 17 bytes a cell. Each step marches only
the cells the pulse can have reached and from which a microphone can still hear before the run
ends; the line the run writes to the standard error gives the share of the cell-steps that
leaves. `--mirror` takes the plane `y` = the source's `y` as one of symmetry and grids half the
space: right only where the buildings are symmetric about it. A microphone on the far side is
taken at its mirror point.

The walls and solids are rigid unless `--wall_absorption` gives their absorption at normal
incidence, the same at every frequency; the floor stays rigid. The air is lossless in the march, and `--humidity`
(with `--temperature` and `--pressure`) applies its absorption to every arrival afterwards by its
delay, as `niceops_point_field` absorbs each ray along its length. `--solids` adds prisms, rigid
or absorbing as the walls are,
from a CSV of `solid`, `x`, `y`, `bottom` and `top` (each solid's vertices in order around it), for
geometry LoD1 buildings cannot hold, such as a facade's relief. `--roofs` gives a CSV district
pitched roofs, as `niceops --roofs` takes them; a cell is solid where its center lies under the
roof, so a slope is a staircase of cells, which diffracts as a smooth slope only well below the
step: behind a gable the staircase leaves the top bands several dB loud
([FIDELITY](FIDELITY.md#pitched-roofs)), so check a sloped roof at a finer `--cell` than a flat one
needs.

The run lasts `--duration` (default twice the time sound takes across the domain's diagonal, the
whole space's under `--mirror`), and
a microphone whose bands still move by more than 0.02 dB over the run's last tenth is named on the
standard error: lengthen the run. Between parallel rigid walls the reflections die away slowly;
give such walls their absorption. On a machine with a GPU (Metal, on a Mac) the march is the
GPU's, in single precision as the CPU's march, wherever the GPU can hold the grid; a grid too large
for it marches on the CPU, and the line on the standard error says so. `--gpu` asks for the GPU and
is refused where there is none or it cannot hold the grid; `--no-gpu` keeps the CPU, which is for
small grids. Check a result by repeating it with a finer `--cell`.

| Option | Meaning |
| --- | --- |
| `--buildings`, `--source`, `--microphones`, `--bands`, `--units`, `-o` | As `niceops_point_field`'s. |
| `--solids` | More prisms, as above. |
| `--roofs` | Pitched roofs for a CSV district, as above. |
| `--cell`, `--points_per_wavelength` | The grid's cell, or how many make a wavelength at the top band's upper edge. Default 8 (12 with `--order 2`). |
| `--order` | The differences' order, 4 (default) or 2. |
| `--margin` | Air beyond everything each way before the absorbing layer. |
| `--pml` | The absorbing layer's thickness, cells. Default 20. |
| `--mirror` | Grid half the space, the plane `y` = the source's `y` one of symmetry. |
| `--wall_absorption` | The walls' (and solids') absorption at normal incidence. Default 0, rigid. |
| `--humidity`, `--temperature`, `--pressure` | The air's absorption, applied by delay. Default none. |
| `--duration` | Simulated time, s. |
| `--courant` | The time step as a fraction of the scheme's limit. Default 0.99. |
| `--pulse_frequency` | Where the pulse's spectrum peaks, Hz. Default the top band's upper edge over 2.25. |
| `--sound_speed` | m/s. Default 343.2. |
| `--spectrum` | Also write the transfer function across the bands, every quarter of the inverse duration: microphone, frequency, real and imaginary parts. |
| `--recording` | Also write the pulse and every microphone's pressure at every step. |
| `--gpu`, `--no-gpu`, `-j` | March on the GPU (the default wherever there is one that holds the grid) or keep the CPU; CPU threads. |

### `niceops_pe_compare`

Runs one path through every model and writes each band's level relative to free field, as CSV
on standard output: a `range_<unit>` column, then per band `straight_<f>`, `linear_<f>`,
`stratified_<f>` and `pe_<f>`, the last the wide-angle parabolic equation (PE), a full-wave
solution with the ground's impedance in the grid. Absorption is left out, and turbulence too
unless `--gamma_t` puts it in the PE. The ground is `--flow_resistance`, `--impedance_model` and
`--porosity_rate`, as for `niceops`. `--tunneling` and `--skip_fill` act on `stratified` as in
`niceops`. With `--shadow_model resolved` it adds a `resolved_<f>` column per band, with
`--shadow_phase` and `--cramped_ground` as in `niceops`; `--fit_linear_boundary`,
`--tangent_decay` and `--no_boundary_law` vary that model's parts, as sensitivities.

```sh
niceops_pe_compare --atmosphere atmospheres/example_windy.csv --toward 270 \
    --source_height 328 --receiver_height 3.937 --bands 125 250 500 1000 --max_range 8000
```

The PE is a split-step Padé march with Collins' self-starter and a fourth-order grid
([THEORY §9.1](THEORY.md#91-the-wide-angle-parabolic-equation)). `--max_angle` (default 50°, at
most 60°) limits it: a level steeper than that from the source is NaN. `--range_step`, `--ppw`
(default 10 points per wavelength, 20 under `--scheme cn` or a second-order grid), `--starter`,
`--vertical_order` and `--boundary_order` tune the grid, and `--scheme cn` is a
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
`--cnel_file`, a period level, `--lamax_file`, `--nax_file`, `--ta_file`, `--csel_file`,
`--coverage_file`, `--provenance_file`, `--full_history`, `--db_cutoff`, or `--propagation_table`
with `--refraction` ([internals/performance.md](internals/performance.md) says why). `--buildings`
samples states its own way ([Not every state is evaluated](#not-every-state-is-evaluated)).

### `--db_cutoff`

`--db_cutoff L` skips any evaluation that provably falls more than 10 dB below `L`, bounding the
level from the loudest sphere the state blends, spherical spreading and the most the ground or a
plate can add. That is exact for any observer whose *peak level* reaches `L`. An SEL exceeds its
peak level by 1 to 12 dB, though, so a cell at or above `L` in SEL can have its window cut: on
the demo case a quarter of the cells at or above the cutoff move by a fraction of a dB
([FIDELITY](FIDELITY.md#--db_cutoff-and-sel)). Read it as a threshold on peak level, with a few dB
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
| `--refraction stratified` (`example_windy.csv`), ground effect | 18 s; 34 s with `--ray_stride 1`; more in ducts |

A broadband run is several times cheaper than spectral. On top of these: `--shadow_model resolved`
and the refraction models' shadow transitions build their tables only once a ray nears a shadow
boundary, up to a second at 8 threads and several at `-j 1`; `--tunneling` adds about 20% and
`--skip_fill` 2–4 times to a stratified run; point receivers over real terrain take about 4 times
the CPU; a propagation table is cheap to read but takes minutes to build; `--epnl_file` and
`--band_frames` evaluate every state.

## Plotting

`-p` shows the footprint with the flight track over it, `-f FILE` saves it, and `-F FILE` saves
the trajectory-state figures. Contours default to 65–100 dB in 5 dB steps (`--cmin`, `--cmax`,
`--cstep`). Multi-trajectory runs skip plotting. Plotting embeds Python, with matplotlib in the
build's own environment (`build/pyenv`), and selects the `Agg` backend itself, so it works
headless. A build without plotting refuses `-p`, `-f` and `-F`; the build options are in the
[README](../README.md).

## Comparing with older results

Results from different versions of NICEOPS can differ by deliberate changes to the models, the
largest of them the removal of [a 4x multiplier on the database's EAA](HISTORY.md#the-4x-eaa-multiplier). Read
[HISTORY.md](HISTORY.md#changes-that-affect-results) before comparing results across versions.
