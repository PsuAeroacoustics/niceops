# NICE-OPS theory

This document states the models NICE-OPS implements, as equations, with the references they
come from. It describes the code as it runs: where NICE-OPS departs from a cited reference, the
departure is stated. [The manual](MANUAL.md) covers the options that select each model and when
to use them. [FIDELITY.md](FIDELITY.md) holds the evidence: how each model was checked, and how
well it agrees with exact solutions and with measurements. The
[developer guide](internals/README.md) covers the implementation.

References are numbered in order of first citation and listed at the end.

## Contents

1. [Overview and notation](#1-overview-and-notation)
2. [The source](#2-the-source)
3. [Trajectories](#3-trajectories)
4. [Spreading and atmospheric absorption](#4-spreading-and-atmospheric-absorption)
5. [Ground effect](#5-ground-effect)
6. [Receivers and terrain](#6-receivers-and-terrain)
7. [The stratified, moving atmosphere](#7-the-stratified-moving-atmosphere)
8. [Refraction](#8-refraction)
9. [The parabolic equation and reference solutions](#9-the-parabolic-equation-and-reference-solutions)
10. [Noise metrics](#10-noise-metrics)
11. [What is not modeled](#11-what-is-not-modeled)
12. [References](#references)

## 1. Overview and notation

NICE-OPS is a source–path–receiver model. For each state of a flight trajectory it forms a
source directivity from a database of measured noise "spheres" [1]. For each observer it
carries the level along a path through the atmosphere to the ground. It then integrates each
observer's level history into exposure metrics. The approach follows the maneuvering-flight
noise model of Greenwood and Rau [2]. Every observer is paired with every state; by default a pair is skipped only where a
bound proves it below the observer's 10 dB-down window, which leaves the result unchanged. The demo
case has 6,551 states and 14,400 observers.

| Symbol | Meaning |
| --- | --- |
| $f$, $\omega = 2\pi f$ | frequency, angular frequency |
| $c$, $k = \omega/c$ | sound speed, wavenumber |
| $L$ | sound pressure level, dB re 20 µPa |
| $\lg$ | base-10 logarithm |
| $r$ | slant range from source to observer |
| $R_s$ | the database's sphere radius |
| $h_s$, $h_r$ | source and receiver heights above the ground plane |
| $d$ | horizontal distance |
| $\beta = 1/Z$ | normalized surface admittance ($Z$ the normalized impedance) |

The time dependence is $e^{-i\omega t}$, so a passive ground has $\mathrm{Re}\,\beta > 0$. SI
units are used throughout. The `--units` option changes only the interface.

## 2. The source

### 2.1 The operating condition

Each trajectory state is reduced to three non-dimensional parameters, following the
quasi-static acoustic mapping (Q-SAM) of Schmitz, Greenwood and co-workers [3, 1, 4]. The
load-factor vector $\mathbf{n}$ is the aircraft's acceleration less gravity, in units of
$g = 9.80665$ m/s²:

$$
\mathbf{n} = \frac{\mathbf{a}}{g} + \hat{\mathbf{z}}, \qquad n = \lvert\mathbf{n}\rvert, \qquad \hat{\mathbf{n}} = \mathbf{n}/n .
$$

$\hat{\mathbf{n}}$ stands in for the main rotor's tip-path-plane normal. With the velocity
$\mathbf{u}$, the tip speed $V_\mathrm{tip}$, rotor radius $R$, weight $W$ and air density $\rho$
from the vehicle file:

$$
\mu = \frac{\lvert\mathbf{u}\rvert}{V_\mathrm{tip}}, \qquad
C_T = \frac{n\,W}{\rho\,V_\mathrm{tip}^2\,\pi R^2}, \qquad
\gamma = \arcsin\frac{\mathbf{u}\cdot\hat{\mathbf{n}}}{\lvert\mathbf{u}\rvert}\ \ (\text{Q-SAM, default}), \qquad
\gamma = \arcsin\frac{u_z}{\lvert\mathbf{u}\rvert}\ \ (\texttt{--no-qsam}).
$$

The Q-SAM path angle is the velocity's angle to the tip-path plane. It is the angle at which a
steady flight condition presents the same rotor inflow, which is what lets steady-flight spheres
represent a maneuvering state. At $\lvert\mathbf{u}\rvert = 0$ both path angles are taken as zero
(§3). $\rho$ is the vehicle file's density in every case, including when `--temperature` or
`--pressure` sets the propagation atmosphere.

**The speed reference.** $\mathbf{u}$ is taken in the reference the database records it was
built with, over the ground if it records none; `--speed_reference` overrides it, and an
air-referenced run needs a wind profile. Through the air, $\mathbf{u}$ is the ground velocity
less the profile's wind at the state's height, turned into the local frame by the compass
rotation of §3. The same velocity sets the sphere's nose (§2.4). The load factor, and so $C_T$,
stays inertial: the rotor's thrust balances the aircraft's acceleration over the ground whatever
the air does.

### 2.2 The sphere database

A database (`.nod`, netCDF) holds one sphere per operating condition $(\mu, \gamma, C_T)$. Each
sphere gives the A-weighted level $L_\mathrm{sph}(\phi, \theta)$, an excess-attenuation factor
$\mathrm{EAA}(\phi, \theta)$, and optionally one-third-octave band levels
$L_{\mathrm{sph},b}(\phi, \theta)$, all on the sphere of radius $R_s$ around the rotor hub. The
spheres are built by PANAM; the shipped spectral databases come from the 2017 Noise Abatement
flight test [5, 6]. Measured levels are propagated back from the microphones to $R_s$
("depropagated") through each run's own atmosphere. The received spectra keep their Doppler
shift (§11).

### 2.3 Blending the spheres

The conditions are normalized into the unit cube,
$\tilde{x} = (x - x_\min)/(x_\max - x_\min)$ for each of $\mu$, $\gamma$ and $C_T$, and
triangulated by an incremental Delaunay triangulation [7, 8]. The triangulation is built on
exact orientation and in-sphere predicates [9], with the error-free sum of [10]. Before
triangulating, every condition on a face of the unit cube is duplicated, with each coordinate
at 0 or 1 moved to −0.1 or 1.1. The copy carries the same sphere, so the clamped edge of the
envelope falls inside the hull. A state's normalized condition $\mathbf{q}$, clamped to the
cube, lies in one simplex with vertices $\mathbf{P}_v$. Its barycentric weights are

$$
w_v = \frac{\left\lvert \det\left(\mathbf{P}_k - \mathbf{q}\right)_{k \ne v} \right\rvert}
           {\sum_{v'} \left\lvert \det\left(\mathbf{P}_k - \mathbf{q}\right)_{k \ne v'} \right\rvert},
$$

and the sphere levels are blended linearly **in decibels**:

$$
L(\phi,\theta) = \sum_v w_v\, L_v(\phi,\theta), \qquad
\mathrm{EAA}(\phi,\theta) = \sum_v w_v\, \mathrm{EAA}_v(\phi,\theta).
$$

A raw weight of 1 or more on a vertex takes that sphere whole. A sphere that has no data in a
direction (stored at or below −999 dB) is dropped from that direction's blend, and the remaining
weights are renormalized. A query still outside the hull takes the nearest (extended) vertex.
The blend depends only on the operating condition, so it is formed once per state.

**Alternatives.** Two other ways to choose the spheres exist for comparison with other models.
`--blend nearest` takes the stored condition nearest the state, whole, at every state. Nearest
is measured in the normalized space the triangulation is built in, with the state clamped into
it first, and a tie goes to the database's earlier sphere. This is the closest a trajectory
model can come to the one-hemisphere-per-segment practice of ECAC Doc 32 [11].

`--blend aam` is the rotorcraft sphere search of the Advanced Acoustic Model, v3.1
[12, §2.5.1.1, pp. 31–34]. AAM has no thrust axis. On a database whose load factors are
spheres, each state searches only the layer of spheres whose $C_T$ is nearest its own (a tie
takes the higher). A `fixed_load_factor` database has one layer, and the load-factor law below
carries $C_T$ as it does under the barycentric blend. Within the layer:

1. The spheres are sorted into rows of one flight path angle, each angle rounded to the nearest
   0.5° first (halves away from zero). A row's angle is its rounded one; the state's is not
   rounded.
2. The rows bracketing the state's path angle are interpolated between if they are at most 3°
   apart on approach ($\gamma < 0$) or 5° apart on departure ($\gamma \ge 0$) [12, Table 2-7].
   A wider pair is not interpolated: the nearer row is taken alone, as is the end row for a
   state past either end. A pair exactly as wide as a limit is still interpolated.
3. Within each row, the speeds $\mu V_\mathrm{tip}$ bracketing the state's are interpolated
   between if they are at most 20 kt apart. Otherwise the nearest speed in the row is taken,
   however far: speed does not disqualify a sphere. A state exactly halfway takes the faster
   sphere.
4. The two rows' speed blends are blended in path angle: up to four spheres, linearly in dB,
   speed first and then path angle [12, Fig. 2-20].

The 20 kt gate is converted to advance ratio with the database's tip speed, or the vehicle
file's for a database that records none. The path angle is the one the run uses (§2.1). The
reference leaves three cases open, and NICE-OPS settles them so:

- **Ties in path angle** go toward level flight: an approach takes the higher (less steep)
  angle, a departure the lower. This follows the labels of [12, Table 2-8], whose signs
  contradict the same document's convention that positive angles climb.
- **A state past every gate** takes the nearest row in path angle, then the nearest speed in
  it, which is the "closest matching flight condition" fallback of [12, p. 31]. Path angle
  comes first; it is not the sphere nearest in a combined distance.
- **Two spheres in one row at one speed** (their angles round alike): the one whose own angle
  is nearer the state's.

Two consequences of these rules follow. The path-angle gate changes at level flight, so a state
between rows 2° either side of level is interpolated while it climbs (5° gate) and takes the
nearer row alone once it descends (3° gate): the output steps there, as it does wherever a
bracketing pair crosses a gate. And the shipped databases pad their envelope with copies of
edge spheres at $\gamma = -24°$ and $35°$, so that the triangulation covers steep states. To
the AAM search those copies are ordinary rows, wider than either gate from the outermost
measured row, so a state more than halfway from that row to a padding row takes the padding
row alone, which repeats edge spheres.

AAM's FRAME option (load factor as a searched axis), its duration adjustment and its
nacelle-angle axis are not reproduced.

**Load factor.** A `fixed_load_factor` database stores only the 1 g spheres and is triangulated
over $(\mu, \gamma)$. The load factor enters through the thrust law

$$
\Delta L_n = 20 \lg \frac{C_T}{C_{T,\mathrm{ref}}}, \qquad C_T > C_{T,\mathrm{ref}},
$$

added to every level and band. $C_T$ is capped at `--max_load_factor` (default 2) times
$C_{T,\mathrm{ref}}$. Below 1 g the law is not applied by default, as PANAM applies it.
`--load_factor_correction two_sided` applies it below 1 g too, down to `--min_load_factor`
(default 0.5), below which the correction is held: near zero load factor, or below it where
$C_T$ is negative, the logarithm would fail.

### 2.4 Reading the sphere

Each sphere is fixed to the rotor. Its axes are $\hat{\mathbf{n}}$, the nose $\hat{\mathbf{e}}$
(a horizontal direction $\mathbf{v}$, below, projected into the tip-path plane) and starboard
$\hat{\mathbf{s}} = \hat{\mathbf{e}} \times \hat{\mathbf{n}}$. When the projection vanishes (no
velocity, or a velocity along $\hat{\mathbf{n}}$), the nose falls back to $+x$ (or $+y$, if
$\hat{\mathbf{n}}$ lies along $x$) projected into the tip-path plane.

**The nose.** By default the database's `azimuth_reference` decides: a database built on the
heading takes `--nose heading` (below), and any other takes the velocity of §2.1
(`--nose velocity`). `--nose ground_velocity` takes the ground velocity whatever the speed
reference. Near hover the
velocity is drift, so when the track carries a heading and the run knows its frame, the nose
turns toward the heading at low speed. With $b$ the frame's bearing (`--frame_bearing`, or 90°
for a geodetic frame settled without one), compass heading $h$ points at the local angle
$\psi_h = b - h$, since compass bearings turn clockwise. With the horizontal speed
$\lvert u_h\rvert$ and the smoothstep weight

$$
w = s^2(3 - 2s), \qquad s = \frac{\lvert u_h\rvert - 3\ \mathrm{kt}}{5\ \mathrm{kt}}\ \text{clamped to}\ [0, 1],
$$

the horizontal angle of $\mathbf{v}$ turns from $\psi_h$ toward the velocity's own by the
fraction $w$, the short way round, and its vertical part is $w\,u_z/\lvert u_h\rvert$. Below 3 kt
the nose is the heading alone, and from 8 kt up it is the velocity. Under `--nose heading` the
horizontal direction is $\psi_h$ at every speed, with the same climb $w\,u_z/\lvert u_h\rvert$, so
a heading along the track gives the velocity's nose; in a crosswind it is turned by the crab
angle. When the run has no frame, the low-speed turn uses instead the circular mean of heading
against track over the states above 20 kt (at least 50 of them), and `--nose heading` is
refused. Without a heading, or without either offset, the nose is the velocity at every speed.
$\mu$, $\gamma$ and $C_T$ do not depend on the nose.

**The angles.** For the unit vector $\hat{\mathbf{d}}$ from source to observer:

$$
\theta = \arccos(\hat{\mathbf{e}} \cdot \hat{\mathbf{d}}), \qquad
\phi = \arctan\!\left(\frac{\hat{\mathbf{s}} \cdot \hat{\mathbf{d}}}{\lvert \hat{\mathbf{n}} \cdot \hat{\mathbf{d}} \rvert}\right) .
$$

This form of $\phi$ holds below the tip-path plane. $\theta$ is 0° ahead and 180° behind.
$\phi$ is 0° straight down along $-\hat{\mathbf{n}}$, +90° to starboard, −90° to port and,
above the plane, ±180° straight up. A direction above the tip-path plane reads its mirror
image below unless every sphere stores its upper half. The level in a direction is interpolated
bilinearly between grid nodes, periodically in $\phi$. $\theta$ has no period: past the last
$\theta$ node the cell joining the last and first nodes extrapolates slightly, a quirk kept
from the D original. The blend in §2.3 is applied at each of the four grid corners first. A
corner with no data after the blend is dropped from the bilinear interpolation in the same way,
and the other corners' weights are renormalized.

## 3. Trajectories

**AC state tracks** (the 27-column NASA/FAA format) carry positions, velocities and
accelerations, and are used as sampled. The format does not agree with itself about the sign of
the vertical speed, so it is decided from the track: with $\dot z_i$ the central difference of
the heights, if both $v_z$ and $\dot z$ have an rms of at least 0.15 m/s and the uncentered
cosine $\sum v_{z,i}\dot z_i / \sqrt{\sum v_{z,i}^2 \sum \dot z_i^2}$ is −0.5 or less, $v_z$ is
negated and the path angles are recomputed. $\mu$ and the load factor do not depend on the sign.

**Static operations** are holds at one point, each for a given duration with the nose on a
compass heading. They are sampled every 0.02 s at 1 g with zero velocity, so $\mu = 0$. At
zero velocity (true hover, flown or static) both path angles are zero: the arcsine of §2.1 is
undefined there, and a stationary aircraft neither climbs nor descends.

**WGS84 tracks** carry positions, and optionally east, north and up velocities. Without
velocities, velocities and accelerations come from a cubic smoothing spline [13] that
minimizes

$$
\sum_i \left(y_i - f(t_i)\right)^2 + \lambda \int f''(t)^2\, dt
$$

per axis. $\lambda$ is chosen so that $\sum_i (y_i - f(t_i))^2 = N\sigma^2$, with $\sigma$ =
`--path_smoothing` (default 10 ft rms per axis, in `--units`) and $N$ the number of points;
$\sigma = 0$ interpolates. The velocity is the spline's first derivative, and the load factor
of §2.1 comes from its second. A dense track (mean spacing 0.5 s or less, $\sigma > 0$) is
padded at each end with 10 s of points on a quadratic fitted to the last 3 s, so the spline's
natural end conditions, which pin the acceleration to zero, do not bend the ends of the track;
the sum and $N$ then include the padded points. If the file gives the velocities, they are
rotated into the frame, interpolated exactly, and differentiated once for the acceleration.
States are resampled at min(`-s`, the file's mean step); a track is never coarsened, and a
uniform track already at or finer than `-s` keeps its own samples.

**Coordinates.** Latitude and longitude are mapped to a local plane by the American polyconic
projection about a reference point [14] on the WGS84 ellipsoid [15]
($a = 6\,378\,137$ m, $f = 1/298.257223563$). $z$ is the altitude less `--ground_elevation`
(the observer plane's elevation), otherwise untransformed. The plane is then rotated so that
$+x$ points along `--frame_bearing` $b$:

$$
x = E \sin b + N \cos b, \qquad y = -E \cos b + N \sin b .
$$

The same rotation, without the translation, turns velocities and compass bearings (wind
directions, headings) into the local frame. Without `--frame_bearing` the axes are east and
north ($b = 90°$).

**Lateral subtracks.** An operations file's `offset` $d$ (positive to the left) moves a track
square to its horizontal direction of travel, as dispersed subtracks do in airport noise models
[57, 58]. With $\mathbf v = (v_x, v_y)$ the horizontal velocity and $\mathbf a = (a_x, a_y)$
the horizontal acceleration, the left normal and the signed curvature (positive turning left)
are

$$
\hat{\mathbf n} = \frac{(-v_y,\ v_x)}{|\mathbf v|}, \qquad
\kappa = \frac{v_x a_y - v_y a_x}{|\mathbf v|^3},
$$

and the displaced path $\mathbf r + d\,\hat{\mathbf n}$ moves at $(1 - d\kappa)$ times the
original's speed, on a radius $1/\kappa - d$. An offset with $1 - d\kappa \le 0.1$ at any
point, 90% of the radius or more toward the inside of a turn, is refused. Where
$|\mathbf v|$ is below 3 kt, $\hat{\mathbf n}$ is taken from the nearest point that is faster.
A slow interval bounded by different travel directions is refused: switching held normals
would introduce a discontinuity in the displaced path. Normals within $10^{-6}$ in Euclidean
norm are treated as the same direction to allow numerical roundoff.

- **WGS84 tracks.** Both derivatives used to calculate $\kappa$ come from a preliminary fit
  of the raw samples at the same $\sigma$. The file's velocity columns, when present, set
  only $\hat{\mathbf n}$ and the 3 kt threshold, not the curvature. Each raw sample is moved
  by $d\,\hat{\mathbf n}$, and the displaced samples are fitted as above.
  The velocity and load factor are the displaced path's, and the file's velocities are not
  used for them.
- **AC tracks.** Each state is moved by $d\,\hat{\mathbf n}$ from its own recorded velocity,
  and keeps its recorded velocity and acceleration: exact on a straight leg, while in a turn
  the speed and load factor stay the centerline's. The $\kappa$ of the fold test comes from a
  fit of the recorded $x$ and $y$ at $\sigma = 10$ ft, not from the recorded acceleration,
  whose sample-to-sample scatter exceeds a turn's centripetal term.

## 4. Spreading and atmospheric absorption

### 4.1 The broadband path

A database without spectra, or a run with `--no-spectral`, propagates the A-weighted level
directly:

$$
L = L_\mathrm{sph}(\phi,\theta) - 20 \lg\frac{\max(r, R_s)}{R_s} -
  \mathrm{EAA}(\phi,\theta)\,\frac{\max(r, R_s) - R_s}{1000\ \mathrm{m}} + \Delta L_n .
$$

EAA is the A-weighted atmospheric attenuation over 1000 m in the database's build atmosphere.
PANAM writes it as the difference between the A-weighted level with and without absorption over
that distance. Spreading is capped at the sphere surface: inside $R_s$ the sphere level is
taken as is.

### 4.2 The spectral path

When the database carries spectra (the default for such a database), each band $b$ is
propagated separately:

$$
L_b = L_{\mathrm{sph},b}(\phi,\theta) - 20 \lg\frac{\max(r, R_s)}{R_s} -
  \alpha_b \max(\ell - R_s, 0) + \Delta L_n - \Delta L_{\mathrm{extra},b} + 10 \lg G_b .
$$

$r$ is the spreading range, $\ell$ the absorbed path length ($\ell = r$ on a straight ray), $G_b$
the ground or receiver power ratio (§5, §6), and $\Delta L_{\mathrm{extra},b}$ any shadow-zone
loss or table correction (§8, §9.3). The bands are then weighted and summed in energy (§10.1).

### 4.3 Atmospheric absorption

The absorption coefficient is that of ISO 9613-1 [16], evaluated at each band's center
frequency. With $T_0 = 293.15$ K, $T_{01} = 273.16$ K, $p_r = 101.325$ kPa, ambient temperature
$T$, pressure $p_a$ and relative humidity $h_r$ in percent:

$$
\frac{p_\mathrm{sat}}{p_r} = 10^{\,-6.8346\,(T_{01}/T)^{1.261} + 4.6151}, \qquad
h = h_r \frac{p_\mathrm{sat}}{p_a},
$$

$$
f_{rO} = \frac{p_a}{p_r}\left(24 + 4.04\times10^4\, h\, \frac{0.02 + h}{0.391 + h}\right), \qquad
f_{rN} = \frac{p_a}{p_r}\left(\frac{T}{T_0}\right)^{-1/2}
         \left(9 + 280\, h\, e^{-4.170\left[(T/T_0)^{-1/3} - 1\right]}\right),
$$

$$
\alpha = 8.686\, f^2 \left[ 1.84\times10^{-11} \frac{p_r}{p_a}\left(\frac{T}{T_0}\right)^{1/2} +
  \left(\frac{T}{T_0}\right)^{-5/2}\left(
   \frac{0.01275\, e^{-2239.1/T}}{f_{rO} + f^2/f_{rO}} +
  \frac{0.1068\, e^{-3352.0/T}}{f_{rN} + f^2/f_{rN}} \right) \right] \ \mathrm{dB/m}.
$$

The sound speed is $c = 343.2\sqrt{T/T_0}$ m/s. The spectral path refuses to run without a
relative humidity: absorption depends strongly on it, so no default is assumed.

**Layered absorption.** Through a stratified atmosphere (§7), $\alpha_b(z)$ is tabulated at 1 m
intervals and integrated exactly for a piecewise-linear coefficient,
$F_b(z) = \int_0^z \alpha_b\, dz'$. The path takes the height-mean coefficient

$$
\bar{\alpha}_b = \frac{F_b(z_s) - F_b(z_o)}{z_s - z_o}
$$

times its length. This is the fine-layer limit of the 30 m layers of 14 CFR 36, Appendix A
[17], and the approach of the sonAIR model [18]. Two approximations come with it. The mean is
taken over the whole height span, including the part inside the sphere radius that the path
length excludes. And on a curved ray the mean is still taken over the straight span
$[z_o, z_s]$, a first-order approximation; the ray's own arc length is what it multiplies.
Pressure in the profile, when not given, is hydrostatic for dry air:
$p(z) = p(0) \exp\!\left(-\frac{g}{R_d}\int_0^z \frac{dz'}{T}\right)$, with
$R_d = 287.05$ J/(kg K), from the surface pressure of the vehicle file or `--pressure`.

## 5. Ground effect

The ground model applies to the spectral path and is on by default there. It returns a gain,
the ratio of the mean-square pressure over the ground to the free-field mean-square pressure.
It is ported from PANAM's `ega()`, which follows the formulation of Stusnick, Plotkin and
Sutherland [19], and is pinned to it at $10^{-9}$ dB. The model assumes a point
source over flat, locally reacting ground, with straight rays.

### 5.1 Ground impedance

The ground is locally reacting, with a normalized impedance $Z$ given by one of four models.
$X = f/\sigma$, with $f$ in Hz and the flow resistivity $\sigma$ in kPa·s/m² (`--flow_resistance`,
default 225, grass):

| Model | Impedance |
| --- | --- |
| Delany–Bazley [20] (default) | $Z = 1 + 9.08 X^{-0.75} + 11.9 i X^{-0.73}$ |
| Miki [21] | $Z = 1 + 5.50 X^{-0.632} + 8.43 i X^{-0.632}$ |
| Variable porosity [22] | $Z = (1+i)\sqrt{\dfrac{\sigma_e}{\pi \gamma \rho_0 f}} + \dfrac{i c\, \alpha_e}{8 \pi \gamma f}$ |
| Hard-backed layer of depth $d_l$ [22] | $Z = Z_c \coth(-i k_c d_l)$ |

In the variable-porosity model $\sigma_e$ is the effective flow resistivity in Pa·s/m², 1000
times `--flow_resistance`, which is given in kPa·s/m² for every model. $\alpha_e$ is the rate of
change of porosity with depth (`--porosity_rate`), $\gamma = 1.4$ and $\rho_0 = 1.2$ kg/m³. For
the layer, $Z_c = 1 + 9.08X^{-0.754} + 11.9iX^{-0.732}$ and
$k_c = k\,(1 + 10.8X^{-0.70} + 10.3iX^{-0.595})$. The models use $\beta = 1/Z$, and a ground with
$\mathrm{Re}\,\beta \le 0$ is refused as non-passive.

**The convention for $\alpha_e$.** The variable-porosity term $c\,\alpha_e/(8\pi\gamma f)$ is what a
first-order (WKB) solution gives for porosity $\Omega_0 e^{-\alpha z}$ with isothermal
compressibility, taking $\alpha_e = \alpha/\Omega_0$. The form usually quoted,
$0.436(1+i)\sqrt{\sigma_e/f} + 19.48\,i\alpha_e/f$, has twice the coefficient, so its $\alpha_e$
is defined differently: halve an $\alpha_e$ fitted in that form before passing it to
`--porosity_rate`. The printed eq. (10) of [22] has not been checked. The 2017 site ground uses
$\alpha_e = 0$.

### 5.2 The spherical-wave reflection coefficient

For a source at height $h_s$ and receiver at $h_r$, a horizontal distance $d$ apart, the direct
and image paths are $R_1 = \sqrt{d^2 + (h_s - h_r)^2}$ and $R_2 = \sqrt{d^2 + (h_s + h_r)^2}$. With
$\cos\vartheta = (h_s + h_r)/R_2$ measured from the normal:

$$
R_p = \frac{\cos\vartheta - \beta}{\cos\vartheta + \beta}, \qquad
Q = R_p + (1 - R_p)\, F(w),
$$

$$
p_w = \sqrt{\frac{i k R_2/2}{1 + \beta\cos\vartheta}}\,(\cos\vartheta + \beta), \qquad
F = 1 + i\sqrt{\pi}\, p_w\, \mathrm{w}(p_w),
$$

where $w = p_w^2$ is the numerical distance and
$\mathrm{w}(z) = e^{-z^2}\mathrm{erfc}(-iz)$ is the Faddeeva function. This is the
Chien–Soroka form of the reflection coefficient [23, 24, 22]. Writing $F$
through $\mathrm{w}(p_w)$, rather than $e^{-w}\mathrm{erfc}(-i\sqrt{w})$, keeps it stable at
any $\lvert w\rvert$. It also keeps the branch of $p_w$ itself, so a mass-like ground
($\mathrm{Im}\,\beta > 0$) does not acquire a spurious growing surface wave.

$\lvert Q\rvert$ is not bounded by 1. Toward grazing incidence the surface wave carried by
$F$ takes it to about 2.02 for an ordinary soil, a gain of up to $(1 + \lvert Q\rvert)^2$,
+9.6 dB, rather than the +6 dB of pressure doubling.

**The Faddeeva function** is computed by Weideman's rational series with $N = 32$ terms
[25]:

$$
\mathrm{w}(z) = \frac{2 P(Z)}{(L - iz)^2} + \frac{1}{\sqrt{\pi}(L - iz)}, \qquad
Z = \frac{L + iz}{L - iz}, \qquad L = \sqrt{N/\sqrt{2}},
$$

with $P$ a polynomial of degree 31. The lower half-plane uses
$\mathrm{w}(z) = 2e^{-z^2} - \mathrm{w}(-z)$. For speed, the hot path uses the asymptotic series
for $\lvert z\rvert \ge 6$ in the upper half-plane or where $x^2 - y^2 \ge 36$ ($z = x + iy$), a
sixth-order Taylor expansion about nodes 0.05 apart for $\lvert z\rvert < 6$ with
$\mathrm{Re}\,z \ge 0$ and $-0.25 \le \mathrm{Im}\,z \le 6$, and Weideman's series for every
other argument. Both fast paths agree with the series to about $10^{-11}$, and the series agrees
with SciPy's `wofz` to $3\times10^{-13}$.

**Roughness.** An rms roughness height $\sigma_r$ (`--roughness`, default 0) reduces the
coherent reflection by the Kirchhoff (Ament) factor [26]:
$Q \leftarrow Q\, e^{-2(k\sigma_r\cos\vartheta)^2}$.

### 5.3 Band averaging and turbulence

A one-third-octave band is not a pure tone. NICE-OPS uses Chessell's band-averaged form
[27], with the corrections noted by Daigle, Embleton and Piercy [28], with the path delay
$\tau = (R_2 - R_1)/c$ and range ratio $\varrho = R_2/R_1$:

$$
G = 1 + \frac{\lvert Q\rvert^2}{\varrho^2} +
  \frac{2}{\varrho}\left[\mathrm{Re}\,Q\cos(\Omega f\tau) - \mathrm{Im}\,Q\sin(\Omega f\tau)\right]
    T(f)\, \frac{\sin(\Delta f\tau)}{\Delta f\tau}, \qquad \Delta L_g = 10\lg G,
$$

with $\Omega = \pi(2^{1/6} + 2^{-1/6}) = 6.325159$ and $\Delta = \pi(2^{1/6} - 2^{-1/6}) = 0.727477$.
These are the interference evaluated at the band's arithmetic center, averaged over its width.
Over a rigid surface with the receiver on the ground, and with no turbulence ($C = 0$, the
default) or under HARMONOISE, $G \to 4$ (+6 dB). Chessell's factor depends on range alone, so
it decorrelates even coincident paths.

$T(f)$ is the coherence of the direct and reflected waves, reduced by atmospheric turbulence.
It uses the direct range $R_1$ in both models:

- **Chessell** (default) [27]: $T = \exp\left[-\left(\tfrac12 C f \sqrt{R_1}\right)^2\right]$,
  with $C$ (`--turbulence`, rad·s/√m) zero by default.
- **HARMONOISE** [29, 22]:
  $T = \exp\left[-\tfrac38\, B\, k^2 \gamma_p^{5/3} R_1\, \gamma_T\right]$, with $B = 0.364$,
  $\gamma_p = h_s h_r/(h_s + h_r)$ and the turbulence strength $\gamma_T$ below.

The HARMONOISE turbulence strength combines the temperature and wind-velocity structure
parameters $C_T^2$ and $C_v^2$ [30]:

$$
\gamma_T = \frac{C_T^2}{T_0^2} + \frac{22}{3}\frac{C_v^2}{c_0^2}\quad (\mathrm{m}^{-2/3}).
$$

### 5.4 Tabulated ground gain

The point-receiver gain of §5.2–§5.3 depends on the geometry only through the source height
and the horizontal distance, the receiver height being fixed for a table. By default it is not
evaluated for every pair. Per trajectory it is computed exactly at the nodes of a grid in
$s = \sqrt{h_s/r}$, $r = \sqrt{h_s^2 + d^2}$, and $\ln r$ (641 values of $s$, 40 nodes per
decade of range), and interpolated bilinearly between them as a power ratio. $s$ resolves
grazing incidence, where the model changes fastest. Ranges outside the table are clamped to its
edge. The interpolant is a convex combination of exact values, so it never leaves the model's
bounds. `--exact_ground` evaluates the model directly. The accuracy of the table is in
[FIDELITY](FIDELITY.md#ground-and-plate-tables).

### 5.5 Ground classes

`--ground_classes` gives a few surfaces $c$, each with its own impedance model, flow resistance,
porosity rate, layer depth and roughness, and an observer file names the class each microphone
stands on. A microphone on class $c$ takes $\beta_c(f) = 1/Z_c(f)$ (§5.1) and its roughness in
§5.2–§5.3, with every other term unchanged. Turbulence is a property of the air and stays the
same for every class. Each table of §5.4 holds one receiver height over one class.

The class is taken as the whole ground under every path to that microphone: a locally reacting
uniform plane. That is right when the patch around the microphone covers the specular point and
its first Fresnel zone. Over flat ground the specular point lies

$$
x_s = \frac{d\,h_r}{h_s + h_r}
$$

from the microphone: about 8 m for a 1.2 m pole with the source 150 m up and 1 km away, and
38 m with the source 30 m up, while the Fresnel zone reaches well past $x_s$ at grazing
incidence and low frequency. A flush plate microphone ($h_r = 0$) reflects at the microphone, so
a patch larger than the plate qualifies. Weighting the surfaces over the Fresnel zone, as
Nord2000 and CNOSSOS-EU do [59, 32, 33], is not done.

For the peak search (§1) the point receiver's bound, $(1 + \max|Q|)^2$ per band, is the
largest over the classes the run's point receivers stand on. It does not depend on the receiver
height, so one bound still serves every observer.

## 6. Receivers and terrain

### 6.1 The point receiver

The default receiver is a point `--receiver_height` above the ground (1.2 m, as 14 CFR 36
prescribes [17]), or each microphone's own height from an observer file. It is §5's model
with $h_r$ set to that height, over the run's ground or over the microphone's ground class
(§5.5).

### 6.2 Ground-plane microphones

`--plate_table` models a microphone mounted in a circular plate on the ground: the GRAS 67AX in a
GR1425 plate, 0.2 m in radius, with the microphone 3/4 of the radius from the center, as in
SAE ARP 4055 [31]. The plate's response is tabulated by PANAM's axisymmetric boundary
element method as complex transfer functions of the direct and image waves, $P_d$ and $P_r$,
over frequency, elevation and azimuth. The band gain averages the table's $N_\mathrm{sub}$
sub-frequencies $f_j$ across the band (five in the GR1425 tables):

$$
G_b = \frac{1}{N_\mathrm{sub}} \sum_{j=1}^{N_\mathrm{sub}} \left\lvert P_d(f_j) + Q(f_j)\, P_r(f_j) \right\rvert^2 ,
$$

with $Q$ the spherical-wave coefficient of §5.2 *without* the roughness factor, evaluated at
the plate's top surface. Neither turbulence model (§5.3) applies, because the plate's image
path is millimeters longer than the direct path. This is the reference form. A run takes $Q$
at each band's edges and center and interpolates it quadratically in $\lg f$ to the $f_j$
(by default from a per-trajectory table of $Q$ over the grazing cosine and the image path's
length). That is within 0.07 dB of the reference above −15 dB re free field, and up to about
1.5 dB off in deep near-grazing nulls. `--exact_ground` evaluates those values of $Q$ directly
instead of from the table.

$P_d$ and $P_r$ are interpolated bilinearly, real and imaginary parts, in elevation and
azimuth. They are held below the table's lowest elevation and above its highest, and periodic
in azimuth. The azimuth is that of the horizontal propagation direction, source to microphone,
in the plate's frame, whose $+y$ is the microphone's offset from the plate's center. The offset
points along a fixed azimuth or, by default, outboard: away from the flight track, along the
horizontal line from the track's nearest state to the observer, or to the track's right for an
observer directly under it. The table records the ground it was computed over, and a run with
a different ground is refused.

### 6.3 Terrain

With `--terrain`, each source–observer pair takes its own mean ground plane: a least-squares
line $z = a + b s$ fitted to the terrain heights along the pair's horizontal path, with $s$ the
horizontal distance from the observer. This is the mean-plane construction of CNOSSOS-EU
[32, 33], also used over terrain in NORAH2 [34] and ECAC Doc 32 [11]. The ground model then runs
in the rotated frame of that plane:

$$
h_s' = \frac{z_s - (a + bD)}{\sqrt{1 + b^2}}, \qquad
h_r' = \frac{z_m - a}{\sqrt{1 + b^2}}, \qquad
d' = \frac{D + b(z_s - z_m)}{\sqrt{1 + b^2}},
$$

with $D$ the horizontal distance and $z_m$ the microphone's height. This is a rotation of the
plane, so the source–microphone distance is unchanged. A negative equivalent height is raised
to zero, as in CNOSSOS-EU. For a flush plate microphone the plate lies along the line under the
microphone: the ground model takes the source height $h_s' - h_r'$ (raised to zero if negative),
and the receiver keeps its plate. A negative $d'$ (a nearly overhead source over ground falling
toward it) is used by its magnitude, which keeps the ground model's direct path equal to the
true range, and the arrival azimuth the plate reads is reversed.

The per-pair line is not fitted pair by pair. It is interpolated linearly in bearing and in
length from per-observer fits along 1440 bearing rays (0.25°), each sampled one raster step
apart, which is exact on a uniformly sloping plane. For an elevated source the specular
reflection point lies near the receiver, at $d'h_r'/(h_s' + h_r')$, so a line fitted over the
whole path also weighs terrain the reflection never sees; that is the standard practice, and on
smooth ground the two agree. Spreading, absorption and the sphere
direction still use the true straight line. On level terrain at $z = 0$ every quantity reduces
to the flat-ground geometry. A propagation table (§9.3) is read in the same rotated frame.
Terrain screening and diffraction are not modeled (§11).

## 7. The stratified, moving atmosphere

An `--atmosphere` profile gives temperature, relative humidity and optionally pressure and wind
against height. Values are linear in height between levels, pressure is log-linear, and the
lowest level is held below the profile. Sound in moving air is treated in the
effective-sound-speed approximation [30, 35]:

$$
c_\mathrm{eff}(z) = c\left(T(z)\right) + \hat{\mathbf{h}} \cdot \mathbf{v}_w(z),
$$

with $\hat{\mathbf{h}}$ the horizontal propagation direction and $\mathbf{v}_w$ the horizontal
wind. The approximation overstates the wind's effect on steep paths by a factor
$1/\cos(\text{elevation})$, which matters little because steep paths refract little.

**A profile from one surface wind** (`--wind`). The wind follows the neutral logarithmic law
[36] through the observation $u_\mathrm{ref}$ at height $z_\mathrm{ref}$,

$$
u(z) = u_\mathrm{ref}\,\frac{\ln(z/z_0)}{\ln(z_\mathrm{ref}/z_0)}, \qquad z_0 < z \le \delta,
$$

held constant above the surface-layer depth $\delta$ (`--surface_layer_depth`, default 200 m)
and zero at or below $z_0$. Its direction is the same at every height. The log law describes
only the neutral surface layer, about a tenth of a boundary layer one to two kilometers deep.
Above it the wind turns and grows slowly toward the geostrophic wind, which a single surface
observation says nothing about, so the profile holds the surface layer's top value rather than
extrapolate the law or model the turning. The roughness length
$z_0$ (`--roughness_length`) defaults to 0.03 m, the Davenport–Wieringa class "open flat terrain" [37, 38]. The
anemometer height $z_\mathrm{ref}$ defaults to the standard 10 m.

As run, the law is sampled at levels. The ground, $z = 0$, is calm. The first level above it is
at $\min(10 z_0, z_\mathrm{ref})$, and from there levels are spaced geometrically,
`--profile_levels` (default 16) per decade, to a top at the larger of $\delta$ and 5% above the
highest state or microphone; $z_\mathrm{ref}$, $\delta$ and the top are levels themselves. The
profile is linear between levels and from zero at the ground to the first level, so no level
above the ground is calm.

The temperature lapse is dry adiabatic (−0.0098 K/m), which makes the layer neutral. Water
vapor holds its surface molar concentration, so the relative humidity follows the cooling air,
capped at 100%. Pressure is hydrostatic (§4.3). The surface temperature, humidity and pressure
are the run's own: the vehicle file's, or `--temperature`, `--humidity` and `--pressure`.

**Stability classes** (`--stability`, unvalidated). Each class replaces the log law with a power
law $u = u_\mathrm{ref}(z/z_\mathrm{ref})^m$, held above $\delta$ as the log law is, and the lapse
with a class temperature gradient. The class gradient holds up to $\delta$, and every class
lapses dry-adiabatically above it. The exponents are the Klug–Manier classes of TA Luft [39],
and the gradients are mid-range values of the $\Delta T$ classes of U.S. NRC Regulatory
Guide 1.23 [40]:

| Class (TA Luft) | Pasquill | Character | $m$ | $dT/dz$ below $\delta$ |
| --- | --- | --- | --- | --- |
| I | F | very stable | 0.42 | +2.75 K/100 m |
| II | E | stable | 0.35 | +0.5 K/100 m |
| III/1 | D | neutral to stable | 0.22 | −1.0 K/100 m |
| III/2 | C | neutral to unstable | 0.20 | −1.6 K/100 m |
| IV | B | unstable | 0.15 | −1.8 K/100 m |
| V | A | very unstable | 0.09 | −2.0 K/100 m |

How the scheme scores is in [FIDELITY](FIDELITY.md#wind-profile-and-stability).

## 8. Refraction

### 8.1 Linear sound-speed profile: circular rays

`--refraction linear` fits each source–observer pair's $c_\mathrm{eff}(z)$ over the pair's
height span (at least 10 m) with a straight line $c(z) = c_0 + a z$ by least squares. This is the
equivalent-linear approach of HARMONOISE and Nord2000 [41, 42]. The rays are
then arcs of circles of radius $\mathcal{L} = c_0/\lvert a\rvert$ centered at height
$-c_0/a$. For the direct ray, with $\kappa(z) = z(2\mathcal{L} + \sigma z)$ and
$\sigma = \mathrm{sign}\,a$, the arc's center lies a horizontal distance
$x_c = [r^2 + \sigma(\kappa_r - \kappa_s)]/(2r)$ from the source. The ray-tube spreading range
is

$$
R_\mathrm{eq}^2 = r \left\lvert \frac{\partial r}{\partial \theta_s} \right\rvert \frac{\sin\lvert\theta_r\rvert}{\cos\theta_s}
 \;\;\Rightarrow\;\; R_\mathrm{eq} = \frac{r}{\cos\theta_s}\ \text{on a circle},
$$

and the travel time is
$t = \mathrm{artanh}\left(r\rho_c/(h_c^2 + x_c r)\right)/\lvert a\rvert$, with
$h_c = \mathcal{L} + \sigma z_s$ and $\rho_c = \sqrt{x_c^2 + h_c^2}$.
The level is read from the sphere at the ray's launch direction, and absorbed over the arc
length (§4.3). The reflected ray is found by iteration on its launch angle, and a point
receiver's ground model is evaluated at the equivalent straight geometry with the ray's grazing
angle $\psi$ and path length; a plate microphone takes the direct ray's arrival angle. Compass
wind directions are turned into the frame by `--frame_bearing` (default 90°, $+x$ east). In
uniform air the rays are the straight lines of §4.

**Shadow zone.** In an upward-refracting atmosphere ($a < 0$) no ray reaches beyond

$$
r_b = t(z_s) + t(z_r), \qquad t(z) = \sqrt{z(2\mathcal{L} - z)} \approx \sqrt{2\mathcal{L} z}.
$$

Beyond it, sound arrives as a creeping wave diffracted along the ground. The level starts from
the limiting ray's: its spreading range is continued in proportion to range,
$R_\mathrm{eq} = r\mathcal{L}/h_c$ ($h_c = \mathcal{L} - z_s$), and its travel time at $c_0$. It decays at the rate of the
first mode [43, §9-5; 44]. Over rigid ground, with $\tau_1 = 1.01879$ the first zero of
$\mathrm{Ai}'(-x)$:

$$
\alpha_\mathrm{sh} = \frac{20}{\ln 10}\,\tau_1 \sin 60^\circ\, \frac{(\pi f)^{1/3} \lvert a\rvert^{2/3}}{c_0}\ \ \mathrm{dB/m}.
$$

Over an impedance ground, the root $b_1$ of $\mathrm{Ai}'(b) + q\,\mathrm{Ai}(b) = 0$, with
$q = -i s\, e^{4\pi i/3} \beta$ and $s = (\pi f/\lvert a\rvert)^{1/3}$, replaces it:
$\alpha_\mathrm{sh} = 8.686\, \mathrm{Im}(b_1 e^{-2\pi i/3})\, \pi f/(c_0 s^2)$. Bands whose
admittance admits a surface wave ($\arg\beta < -60°$) cannot follow that root and take the
rigid-ground rate. Real ground is not always the deeper shadow. To first order the impedance
rate is up to 12% below the rigid one at low frequency and weak gradients, and up to 2.3 times
it once the ground is soft on the wave's scale: at 1 kHz in the benchmark's case 3 (§9.2, $c$
falling 0.1 /s) it is 146 dB/km against rigid ground's 70. The loss $\min(\alpha_\mathrm{sh}(r - r_b), L_\mathrm{floor})$, with
`--shadow_floor` $L_\mathrm{floor}$ (default 20 dB), is subtracted. The floor stands for the
energy turbulence scatters into the shadow.

**The shadow transition.** The ray picture is wrong near the boundary, where the field
diffracts smoothly from lit to shadow. Its ground term, taken at the grazing ray's straight
equivalent geometry, sinks toward the grazing ground wave over soft ground and doubles over
rigid ground, where the exact field does neither, and near the ground the equivalent-linear
method can depart from a PE by up to 10 dB [45]. When a ground model is on, NICE-OPS replaces
the ray model there by a tabulated field.

For each band, and each microphone height, it tabulates $V(\xi, s)$: the level of the
residue-series field of a profile with $k^2$ linear in height over the impedance plane (§9.2),
relative to the ray's spreading alone. Its arguments are the Fock-scaled distance past the
boundary $\xi = (r - r_b)/F$, with the diffraction length $F = 2s^2/k_0$, and
$s = k_0 \ell$, $\ell = (c_0/(2\lvert a\rvert k_0^2))^{1/3}$. Relative to the ray's spreading it
hardly depends on the source's height, which is Fock's result [44], so one table serves every
path. The table is built with the source at $\min(20, 0.1 s^2)\,\ell$, over
$-2 \le \xi \le 6$ in steps of 0.1 and $5 \le s \le 3000$ at 8 nodes per decade. $V$ is capped at
$20\lg 2$, the most two coherent rays give, held at $V(0) - L_\mathrm{floor}$ deep in the
shadow, and extrapolated at its last slope past $\xi = 6$.

The level is the ray model's, ground term included, up to $\xi = -2$, and the ray's spreading
plus $V$ from $\xi = -1$ on, where $V$ replaces both the ground term and the creeping-wave
loss. Between them the two are blended linearly in $\xi$. A cramped band, one whose boundary
lies within two diffraction lengths of the source ($r_b < 2F$), keeps the ray model up to its
boundary and blends over $0 \le \xi \le 1$ instead: the table has no lit zone to hand over from.
Bands with $s < 5$ keep the ray model and the creeping-wave loss. A flush plate microphone takes
the field for a receiver at the plate's top, plus a tabulated offset $K$(band, azimuth): the
plate's level over bare ground near grazing (computed at 2° elevation, where it barely depends
on elevation), every 10° of azimuth and interpolated linearly to the pair's own. Without a
ground model there is no transition, and the creeping-wave loss above applies.

**Each band's boundary.** Through layers (§8.2) a pair has two boundaries, the layered rays'
$r_s$ and the linear fit's $r_l$. A wave resolves the thin fast layer near the ground that sets
$r_s$ only once its diffraction length is short against the gap between them, so each band takes

$$
r_b = r_s + w\,(r_l - r_s), \qquad w = \frac{x}{x + C}, \qquad x = \frac{F_l}{r_l - r_s},
$$

with $F_l$ the band's diffraction length at the fit's curvature, $C = 0.363$ for a raised
receiver and $C = 0.75$ for a flush one (below 5 cm, or a plate). Where the fit casts no shadow
($r_l$ infinite, as over a strong inversion), $r_b = m\,(f/1\ \mathrm{kHz})^{-0.08}\, r_s$, with
$m = 1.60$ flush and $1.90$ otherwise. The band is then treated as a linear profile whose
boundary is $r_b$, gradient $\lvert a\rvert = c_0(\sqrt{2z_s} + \sqrt{2z_r})^2/r_b^2$, for
$s$, $\xi$ and the creeping-wave rate. With $r_s = r_l$, a linear profile, it is the one
boundary. The law and its constants are a fit to the parabolic equation particular to NICE-OPS;
the calibration is in [FIDELITY](FIDELITY.md#shadow-transition-calibration).

**Turbulent scattering** (`--scattering`, off by default). Turbulence scatters energy into the
shadow. With `--scattering` $\gamma_T$, the scattered level

$$
S = \min\!\left(0,\ a_0 + 10\lg\gamma_T + b_0 \lg\frac{f}{1\ \mathrm{kHz}} +
  c_0' \lg\frac{R}{100\ \mathrm{m}} + d_0 \lg\frac{z_s}{100\ \mathrm{m}}\right)
$$

relative to free field, with $R$ the slant range, is added in energy to the coherent field $V$,
and replaces the `--shadow_floor` hold. $R$ is clamped to [100, 10 000] m and $z_s$ to
[10, 300] m, and the energy sum is capped at the table's maximum. A plate adds $K$ to $S$ as it
does to $V$. The coefficients are NICE-OPS's own regression on parabolic-equation runs with
von Kármán phase screens (§9.1) through measured upwind profiles, not a published formula. The
flush coefficients, $(a_0, b_0, c_0', d_0) = (39.18, -1.55, 2.10, 5.30)$, apply below 5 cm or
on a plate, and $(36.67, 5.61, 5.13, 5.30)$, fitted for a receiver 1.2 m up, to every other
receiver. Scattering acts through the shadow transition, so it needs a ground model. How
strong a $\gamma_T$ the measurements support is in
[FIDELITY](FIDELITY.md#turbulent-scattering-strength).

### 8.2 Stratified profile: ray tracing through layers

`--refraction stratified` traces rays through $c_\mathrm{eff}(z)$, taken as piecewise linear
between the profile's levels and constant below the lowest and above the highest. Snell's law
fixes the ray parameter $p = \cos\theta/c$. Across a layer of gradient $g$ from $c_0$ to $c_1$,
with $q = \sqrt{1 - p^2c^2}$, the ray advances

$$
\Delta x = \frac{p\,(c_0 + c_1)\,\Delta z}{q_0 + q_1}, \qquad
\Delta t = \frac{1}{g} \ln\frac{c_1 (1 + q_0)}{c_0 (1 + q_1)} .
$$

The layered solve applies when the source is above the receiver ($z_s > z_r$) and the
horizontal distance is at least 1 mm; otherwise the pair takes the linear model's ray.

**Eigenrays.** The eigenrays are the roots of $X(p) = r$ over four direct families: launched up
or down, arriving on the first crossing of $z_r$, or on the second after dipping below $z_r$ and
turning back up without touching the ground. The primary arrival is the downward family's
root, found from the linear model's answer, and each root is solved by the Illinois variant of
regula falsi. Range need not be monotonic in $p$ along a family, so the other families are
sampled in launch angle, 32 times plus 8 for each bounce the path can hold (at most 512), and
every crossing of $r$ is solved for. The reflected eigenray has its own two families (launched
down or up). Rays launched up reach the receiver from aloft only if the air aloft outruns all
the air between source and receiver, and only past their skip distance; short of it the
ground-grazing boundary stands.

The spreading of each arrival comes from the ray tube,
$R_\mathrm{eq}^2 = r\lvert dX/d\theta_s\rvert\sin\lvert\theta_r\rvert/\cos\theta_s$, with any
focusing gain capped at 10 dB per ray. In the shadow the limiting ray's tube is carried out in
proportion to range and its travel time continued at $c(0)$, as in §8.1.

**Ducts.** When the air above the source outruns all of the air below it, a ray that reaches
the ground turns aloft and comes back down. Rays that bounce $m$ times arrive at
$X_m = X_0 + m L_\mathrm{hop}$, with $L_\mathrm{hop}$ one hop from the ground to the turning
height and back, and are found from the same samples. Each intermediate bounce costs the
plane-wave loss $-10\lg\lvert R_p(\psi)\rvert^2$, and the last bounce takes the ground model.
Past the reach of every direct ray, the strongest bounced arrival becomes the primary. Every
solved arrival is kept, each through its own launch direction on the sphere and its own ray
tube, and they are summed **incoherently**, $L = 10\lg\sum_i 10^{L_i/10}$.

**The resolution limit.** A ray whose invariant is $1/c_\mathrm{max}$, with $c_\mathrm{max}$ the
fastest air on the path, runs horizontally where the air is that fast. In flat or nearly flat
air (a calm floor at the ground, for instance) ray theory then sends it, and the shadow
boundary with it, out without limit. A wave confined to a height $h$ diffracts out of it within
about $k h^2$. So the limiting ray's turning region is taken as the heights within one radian of
vertical phase of the fastest air, $k\int\sin\theta\,dz = 1$, and the ray's horizontal travel
across that region, of height $h$, is capped at $1.5\,k h^2$, with $k$ at the run's highest band.
In a linear gradient ray theory travels $\tfrac43 k h^2$ there, so the cap never binds on a
resolved turning. The cap is not applied where the fastest air is the source's own and there is
faster air above it (downward refraction). The phase is integrated beyond the path's ends,
through the profile, when the path ends before the radian is gathered.

**The hand-over and the shadow.** When a shadow transition serves the receiver (a ground model
is on) and the linear fit's boundary $r_l$ lies beyond the layered boundary $r_s$, the pair takes
the linear fit's ray from $0.9\,r_s$ onward, lit or in the fit's shadow, carrying both
boundaries for §8.1's per-band boundary. Near its boundary the layered ray focuses into a
caustic the wave does not show. Without a transition there is no hand-over, and the layered
shadow ray takes the creeping-wave loss at the fit's gradient. A boundary set by faster air
aloft is withdrawn where a ray from aloft or along a duct does arrive.

**Wave corrections.** By default a caustic's focusing gain is capped at 10 dB. Two opt-in
Airy corrections exist. `--tunneling` adds rays that turn back at a height $z_t$ just above the
receiver, which the wave reaches as the evanescent tail $\mathrm{Ai}^2(x)$,
$x = (\tfrac32\omega I)^{2/3}$, $I = \int_{z_r}^{z_t}\sqrt{p^2 - 1/c^2}\,dz$; near the family's
caustic the pair of rays is replaced by the caustic's band-averaged Airy field. `--skip_fill` adds
a caustic's tail on the dark side of an interior extremum of a family's range,
$I = C\,\mathrm{Ai}^2\!\left(\omega^{2/3} d\,(2/X'')^{1/3}\right)$, with $d$ the distance short of
the caustic and $C$ matched to the two rays' energy on the lit side. Each is averaged over 7
sub-frequencies of the third octave and normalized to the ray tube at the 10 kHz reference
frequency, so it is never a gain over the tube [43].

### 8.3 The resolved shadow

`--shadow_model resolved`, with `--refraction stratified` and a ground model, changes how each
band's shadow boundary is set; the default, `legacy`, is §8.1–§8.2. Ray theory fixes a
ground-grazing boundary from $c(z_s)/c(0)$ alone, so a layer near the ground thinner than any
band resolves can set it. Under this option each band sees only the profile it resolves. For
the band at $f$ it finds the height $h_f$ where

$$
k_f \int_0^{h_f} \sqrt{\frac{2\max\left(c(0) - c(z), 0\right)}{c(0)}}\, dz = \Phi ,
$$

with $\Phi$ = `--shadow_phase` (default 1 rad, the radian of the resolution limit in §8.2).
Below $h_f$ the effective sound speed is replaced by its tangent there, the slope of the layer
containing $h_f$. The band's layered boundary $r_s(f)$ comes from that profile, and its
linear-fit boundary $r_l(f)$ is the path's fit corrected by the least-squares line through the
tangent's change to the profile. Both go through the boundary law of §8.1. A band whose
integral never reaches $\Phi$ resolves no upward refraction near the ground. On a linear
profile

$$
h_f = \left(\frac{3\Phi}{2k\sqrt{2\lvert a\rvert/c_0}}\right)^{2/3},
$$

1.31 creeping-wave layer thicknesses $\ell$ at one radian, the resolved profile is the profile
itself, and nothing changes. The decay keeps the boundary's curvature: the creeping-wave rate,
$s$ and the diffraction length are those of the linear profile whose boundary is $r_b(f)$.

A cramped band ($r_b < 2F$, which is $\sqrt{y_s} + \sqrt{y_r} < 2$ in Airy lengths) takes, under
`--cramped_ground series` (the default), the residue series of §9.2 at the source's own height,
with the curvature whose boundary is the band's. It is tabulated per band over $\xi$, $s$ and
$\sqrt{y_s}$ ($0$ to $2$ every 0.4; $\xi$ every 0.2). Near the source, where the series is not
used, the band falls back to the ray model. `--cramped_ground straight` takes the ground model
along the straight path through the lit side and the blend instead, and `ray` the ray model as
`legacy` does.

The per-band boundaries depend only on bearing, source height and receiver height, so a run
tabulates them once per trajectory, every 2° of bearing and at nodes even in $\sqrt{z_s}$, and
interpolates in $1/r_b$. Scattered observers are looked up at their receiver height above the
ground. The comparison with the PE and with flight data is in
[FIDELITY](FIDELITY.md#resolved-shadow-against-the-pe).

## 9. The parabolic equation and reference solutions

The parabolic equation (PE) never runs inside a footprint. It checks the refraction models
(`niceops_pe_compare`, `niceops_benchmark`) and builds propagation-loss tables (§9.3). The
analytic residue series of §9.2 does run in a footprint, at setup: it tabulates the
shadow-transition field of §8.1 and the cramped bands of §8.3.

### 9.1 The wide-angle parabolic equation

With $\psi = p\sqrt{r}\,e^{-ik_a r}$ and reference wavenumber $k_a = \omega/c(0)$, the one-way
equation is

$$
\frac{\partial \psi}{\partial r} = i k_a\left(\sqrt{1 + \mathcal{Q}} - 1\right)\psi, \qquad
\mathcal{Q} = \frac{1}{k_a^2}\left(\frac{\partial^2}{\partial z^2} + k^2(z) - k_a^2\right),
$$

with $k(z) = \omega/c_\mathrm{eff}(z) + i\alpha$ and the impedance condition
$\partial\psi/\partial z = -ik\beta\psi$ at the ground. It is marched in range by the split-step
Padé method of Collins [46]: $\exp[i k_a \Delta r(\sqrt{1+\mathcal{Q}} - 1)]$ is replaced
by a rational function $r_0 + \sum_j c_j/(\mathcal{Q} - p_j)$, so a step costs one tridiagonal
solve per term. **Departure:** the coefficients are not the Taylor–Padé coefficients. They come
from a rational interpolant fitted at run time on Chebyshev nodes over the propagation angles
needed (`--max_angle`, default 50°, up to 60°), with the fewest terms that meet a $10^{-7}$
tolerance. Collins's stability constraint enters as two further nodes among the evanescent
modes, at $\mathcal{Q} = -1.5$ and $-6$, where the interpolant is pinned to zero, so modes
steeper than the accurate interval decay within a few steps instead of lingering. Waves steeper
than `--max_angle` are therefore absorbed rather than carried, and a level on a path steeper
than that from the source is reported as NaN. The default range step is eight wavelengths.

The depth operator is discretized to fourth order (Numerov), with a fourth-order impedance row
at the ground and fourth-order interpolation that places the source on the grid and reads the
receivers off it. The grid is 10 points per wavelength (20 with a second-order option or
Crank–Nicolson). A Crank–Nicolson scheme, the Padé(1,1) case [35, App. F], is available as an
independent check. The field starts from Collins's self-starter [47], the first range step
applied to a delta function at the source, or from the Gaussian starter of [35]. An artificial
absorbing layer 100 wavelengths thick closes the top of the domain.

**Turbulence** enters as random phase screens [48; 35, App. I], one every outer scale
$L_0$ (default 10 m) of range. Each screen's phase spectrum is that of a von Kármán spectrum
with structure parameter $\gamma_T/4$, the value for which Ostashev's coherence and HARMONOISE's
reduce to one form [22]:

$$
S_\varphi(\kappa_z) = 2\pi k_a^2\, \Delta x\; 0.033\, \frac{\gamma_T}{4}\;
 \frac{\sqrt{\pi}\,\Gamma(4/3)}{\Gamma(11/6)}\;\left(\kappa_z^2 + L_0^{-2}\right)^{-4/3} ,
$$

the 3-D spectrum integrated across the plane, so the vertical coherence is the 3-D medium's.
The screen's phase variance is $0.782\,k^2\Delta x\,(\gamma_T/4)\,L_0^{5/3}$. A height-dependent
strength $\gamma(z)$ scales each row's phase by $\sqrt{\gamma(z)/\gamma_T}$, keeping the vertical
spectrum of a homogeneous field.

The relative level is $\Delta L = 20\lg(\lvert\psi\rvert R/\sqrt{r})$ relative to free field.
Band levels are energy means over sub-frequencies (five by default) at the midpoints of equal
slices of the band in log frequency. With turbulence each sub-frequency is its own realization,
and further realizations can be averaged in, so a band level averages over the band and the
ensemble at once.

**The scattering fit.** The coefficients of $S$ in §8.1 were fitted to the turbulent PE's level
deep in the shadow, at points where it is 10 dB clear of the still-air field.

### 9.2 Analytic reference solutions

These are the benchmark cases of Attenborough et al. [49]: a source at 5 m, a receiver at 1 m,
10, 100 and 1000 Hz, a locally reacting ground with the paper's four-parameter impedances, and
absorption entered as $\mathrm{Im}\,k_0$ from the values the paper quotes:

| Case | Profile | Exact solution |
| --- | --- | --- |
| 1 | uniform | direct and image rays with the spherical-wave reflection coefficient |
| 2 | $c$ rising 0.1 /s (downward refraction) | normal-mode series |
| 3 | $c$ falling 0.1 /s (upward refraction) | residue series of creeping waves |
| 4 | rising 0.1 /s to 100 m, falling to 300 m, uniform above (a duct) | none |

In homogeneous air the field is $p = e^{ikR_1}/R_1 + Q\,e^{ikR_2}/R_2$ (§5.2). In a linear
profile, $k^2(z) = k_0^2(1 \mp 2z/\mathcal{R})$, the field is a sum of normal modes,

$$
p = \pi i \sum_n \frac{u_n(z)\,u_n(z_s)}{N_n}\, H_0^{(1)}(\kappa_n r), \qquad N_n = \int_0^\infty u_n^2\, dz,
$$

with $\ell = (\mathcal{R}/(2k_0^2))^{1/3}$ and $\mathrm{w}_3 = e^{2\pi i/3}$:

- **Downward refraction:** $u_n = \mathrm{Ai}(t_n + z/\ell)$, with
  $\mathrm{Ai}'(t_n) + i k_0 \ell \beta\, \mathrm{Ai}(t_n) = 0$.
- **Upward refraction:** $u_n = \mathrm{Ai}(b_n - \mathrm{w}_3 z/\ell)$, with
  $\mathrm{Ai}'(b_n) - i k_0 \ell \beta\, \mathrm{w}_3^2\, \mathrm{Ai}(b_n) = 0$.

The Airy and Hankel functions follow the DLMF [50].

**Departure: the upward mode condition.** The condition printed in [49, Eq. 11],
$\mathrm{Ai}'(b) + i k_0 \ell \beta\,\mathrm{w}_3\,\mathrm{Ai}(b) = 0$, differs from the one derived
here by a factor $e^{i\pi/3}$ on the admittance term. The two agree over rigid ground. Over
impedance ground the printed form departs from the PE and the derived form agrees with it
([FIDELITY](FIDELITY.md#the-attenborough-et-al-benchmark)). Eq. (10) of [49] is right up to an
overall sign.

**Case 2's exact profile.** The normal-mode solution's profile, $c_0/\sqrt{1 - 2z/\mathcal{R}}$,
reaches infinite sound speed at $z = \mathcal{R}/2$, so every ray turns below that height,
near-vertical ones included. Those fall back within $\mathcal{R}$ of the source, at
$\mathcal{R}\sin(2 \times \text{elevation})$, and put structure into the mode sum that no PE,
and no linear profile, carries. The unit test's comparison with the PE leaves out modes steeper
than 75° (tapered over the last 5°); `niceops_benchmark` keeps them all. At long range
case 2 also depends on where the profile is cut off, as the paper says.

### 9.3 Propagation-loss tables

A `--propagation_table` stores $\Delta L(\text{bearing}, z_s, z_r, r, f)$ relative to free field
at the straight-line distance, over one stratified atmosphere and flat, uniform ground. It is
computed through one profile either by the PE (by default marched by reciprocity, one march per
receiver height) or by straight rays with the spherical-wave ground model
(`--method ground_model`). Reciprocity holds because, in the effective-sound-speed model, the
field from $z_s$ to $z_r$ along a bearing equals the field from $z_r$ to $z_s$ along the reverse
bearing with the wind reversed, which has the same effective profile; a source near the ground
is placed as the adjoint of the receivers' read.

A run interpolates the table multilinearly in dB and applies it in place of the ground and
refraction models, with the straight-line direction and range:

$$
L_b = L_{\mathrm{sph},b} - 20\lg\frac{\max(r, R_s)}{R_s} - \alpha_b \max(r - R_s, 0) + \Delta L_n + \Delta L .
$$

A table that declares `includes_absorption` drops the $\alpha_b$ term. A pair the table does
not cover (outside an axis, or a NaN cell) falls back to the run's own ground and refraction
model. A flush plate microphone (§6.2) reads the table at its lowest receiver height $h_0$ and
adds

$$
K(f, \phi, \varepsilon) = L_\mathrm{plate}(f, \phi, \varepsilon) - L_\mathrm{bare}(f, \varepsilon; h_0),
$$

the plate's level (§6.2) less bare ground's at $h_0$ (§5), for a source on a 1 km sphere at azimuth $\phi$ in the plate's frame and elevation $\varepsilon$, interpolated
bilinearly in $\phi$ and $\varepsilon$. Unlike the shadow transition's $K$(band, azimuth) of
§8.1, taken at 2°, it depends on elevation, because over a raised $h_0$ the bare term is a comb
in elevation that $K$ has to follow; over $h_0 = 0$ it is the plate's offset from flush.

Over terrain (§6.3) the table is read in each pair's rotated frame:
$\Delta L(\beta, h_s', h_r', |d'|, f)$ for a point receiver, and
$\Delta L(\beta, \max(0, h_s' - h_r'), h_0, |d'|, f)$ plus $K(f, \phi', \varepsilon')$ for a flush
plate, with $\varepsilon' = \arctan(\max(0, h_s' - h_r')/|d'|)$ and $\phi'$ the plate-frame
azimuth of the separation as §6.3 turns it when $d' < 0$. The bearing $\beta$ is the true
horizontal separation's, not the turned one's, so an anisotropic table stays continuous as $d'$
passes through zero. Spreading, absorption and the sphere direction keep the true straight line,
and since the rotation preserves distance the table's free-field reference is the true range
except where a height was raised to zero. A point receiver whose $h_r'$ lies outside the table's
receiver heights is not clamped onto a row; it falls back to the run's own model. Reading a
refracting table this way takes its stratification to follow each path's mean plane. CNOSSOS-EU
evaluates its favorable (downward-refracting) conditions over the mean plane in the same way
[32, 33]. For a linear profile the circular ray tangent to a plane tilted by $\alpha$ has a
radius larger by $1/\cos\alpha$ (0.14% at 3°), so the read is a small perturbation of the
vertical-gradient geometry on gentle slopes.

Turbulence and ground roughness are the table's own. A ground-model table built with
Chessell's $C$ or HARMONOISE's $\gamma_T$ carries the coherence $T(f)$ of §5, and one built with
a roughness the Kirchhoff factor; a PE table's turbulence is its phase screens (§9.1). The run's
`--turbulence`, `--gamma_t` and `--roughness` reach only the pairs that fall back.

## 10. Noise metrics

### 10.1 Frequency weighting and band summation

The A and C weightings are the analytic forms of IEC 61672-1 [51]:

$$
W_C(f) = 10\lg\frac{K_C\, f^4}{(f^2 + f_1^2)^2 (f^2 + f_4^2)^2}, \qquad
W_A(f) = W_C(f) + 10\lg\frac{K_A\, f^4}{(f^2 + f_2^2)(f^2 + f_3^2)},
$$

with $f_1 = 20.598997$ Hz, $f_2 = 107.65265$ Hz, $f_3 = 737.86223$ Hz, $f_4 = 12194.22$ Hz,
$K_C = 2.242881\times10^{16}$ and $K_A = 1.562339$ (the constants normalize each to 0 dB at
1 kHz). Each band is weighted at its center frequency, and the bands inside `--band_limits`
(default 20 Hz–20 kHz) are summed in energy:

$$
L_A = 10\lg \sum_b 10^{\left(L_b + W_A(f_b)\right)/10} .
$$

### 10.2 Sound exposure level

The level history is sampled at emission times $t_i$, a constant step $T$ apart (the
trajectory's own interval), but exposure accumulates in reception time $t_i + r_i/c$. Each
sample therefore lasts $T(1 + \dot r_i/c)$, compressed while the aircraft approaches and
stretched while it recedes:

$$
\mathrm{SEL} = 10\lg\left( \frac{T}{t_\mathrm{ref}} \sum_{i=i_1}^{i_2} \left(1 + \frac{\dot r_i}{c}\right) 10^{L_i/10} \right),
\qquad t_\mathrm{ref} = 1\ \mathrm{s},
$$

with $\dot r_i$ a central difference. The sum runs over the 10 dB-down window, from the first to
the last sample within 10 dB of the peak, dips between them included, as in 14 CFR 36
[17, A36.4.5]. `--full_history` sums every sample, here and in the day metrics (§10.4). Under
refraction $r_i$ is $c$ times the ray's travel time. An ambient level $L_\mathrm{amb}$ (`-a`) is
added in energy, $10\lg(10^{L/10} + 10^{L_\mathrm{amb}/10})$, before the window is found.

### 10.3 Maximum level, time above and number above

The signal is piecewise constant in mean square, $P_i = 10^{L_i/10}$, over the reception-time
durations above. The exponential time weighting of IEC 61672-1 [51] (fast, $\tau = 0.125$ s,
by default) is applied exactly for that input:

$$
y_i = P_i + (y_{i-1} - P_i)\, e^{-T(1 + \dot r_i/c)/\tau}, \qquad
L_{AF\mathrm{max}} = 10\lg\max_i y_i ,
$$

starting settled on the first sample. Across a sample the weighted level relaxes
monotonically, so its maximum falls at a sample's end. Slow weighting ($\tau = 1$ s) and a 0.5 s
running $L_\mathrm{eq}$, the maximum used by ECAC Doc 32 [11], are also available; the running
average is linear between the instants where either edge of its window crosses a sample
boundary. So is `instantaneous`, the largest sample with no time weighting. It is not a meter
quantity. A track whose mean step exceeds $\tau/2$ is refused under each time weighting
(0.25 s for the running $L_\mathrm{eq}$); `instantaneous` takes any step.

Time above a threshold (TA) is the reception time during which the weighted level strictly
exceeds it, solved for analytically within each sample (each sample counts whole under
`instantaneous`). The number of events above (NAx) counts the operations whose maximum reaches
the threshold. Over several operations, the maximum level is the largest of theirs and TA adds.
A trajectory $j$ standing for $N_j$ operations (an operations file's count, or its period
counts added up as in §10.4; 1 otherwise) is counted that many times:

$$
\mathrm{NA}x = \sum_j N_j \left[L_{\mathrm{max},j} \ge x\right], \qquad
\mathrm{TA} = \sum_j N_j\, \mathrm{TA}_j, \qquad
L_\mathrm{max} = \max_j L_{\mathrm{max},j} .
$$

The coverage report's two energies and the provenance's per-state energies are weighted by
$N_j$ in the same way, so their shares are $N_j$-weighted ratios over the trajectories.

### 10.4 Day–night metrics

Each operation's exposure is computed with a penalty $P$ added to every sample heard in a
penalized period, by local time at emission (the propagation delay is not added). The 10 dB-down
window is found on the penalized history, or the whole history is summed under
`--full_history`. The operations are then summed and normalized to a 24-hour day:

$$
L_\mathrm{dn} = 10\lg\left(\frac{1}{86\,400\ \mathrm{s}} \sum_j N_j\, t_\mathrm{ref}\, 10^{\mathrm{SEL}_{P,j}/10}\right),
$$

where trajectory $j$ stands for $N_j$ operations: an operations file's count, or 1. The count
multiplies the energy after the exposure is formed; it does not enter the levels, the window
or the penalty. An operations file's `local_time` sets the clock the penalties are read from.

| Metric | Evening | Night | Definition |
| --- | --- | --- | --- |
| DNL | — | 22:00–07:00, +10 dB | 14 CFR Part 150 [52] |
| Lden | 19:00–23:00, +5 dB | 23:00–07:00, +10 dB | Directive 2002/49/EC [53] |
| CNEL | 19:00–22:00, +10 lg 3 dB | 22:00–07:00, +10 dB | Cal. Code Regs. tit. 21, §5001 [54] |

A period's opening hour belongs to it. The directive's form,

$$
L_\mathrm{den} = 10\lg\left[\frac{12\cdot10^{L_d/10} + 4\cdot10^{(L_e+5)/10} + 8\cdot10^{(L_n+10)/10}}{24}\right],
$$

with each $L$ the equivalent level over its own period, is the same number as the penalized
energies summed over 86 400 s. Title 21's weighting of each evening operation three times and
each night operation ten times is the same again.

The period levels (`--lday_file`, `--levening_file`, `--lnight_file`) are those $L_d$, $L_e$ and
$L_n$, and `--laeq_file` is the same over a clock window $[t_1, t_2)$ of length $T_w$:

$$
L_{\mathrm{Aeq},T} = 10\lg\left(\frac{1}{T_w} \sum_j N_j\, t_\mathrm{ref}\, 10^{\mathrm{SEL}_{w,j}/10}\right),
$$

where $\mathrm{SEL}_w$ is the operation's exposure over its samples inside the window (its 10 dB-down
window found among them), $T_w$ is 12 h, 4 h and 8 h for the three periods, and no penalty
applies. An operation lying wholly in one period recombines into $L_\mathrm{den}$ exactly; one
that straddles a period's edge does so exactly under `--full_history`, and to the difference of
the windows otherwise.

Period counts replace the clock. A trajectory given $N_{d,j}$, $N_{e,j}$ and $N_{n,j}$
operations by day (07:00–19:00), evening (19:00–23:00) and night (23:00–07:00) contributes, to
metric $m$,

$$
w_{m,j}\, t_\mathrm{ref}\, 10^{\mathrm{SEL}_j/10}, \qquad
w_{m,j} = a_{m,d} N_{d,j} + a_{m,e} N_{e,j} + a_{m,n} N_{n,j},
$$

with the unpenalized $\mathrm{SEL}_j$ and $a_{m,p} = 10^{P_{m,p}/10}$ the weighting of an
operation flown in period $p$ (1 inside a period level's period, 0 outside):

| Metric | $a_d$ | $a_e$ | $a_n$ |
| --- | --- | --- | --- |
| DNL | 1 | 1 | 10 |
| Lden | 1 | $\sqrt{10}$ | 10 |
| CNEL | 1 | 3 | 10 |
| $L_d$ | 1 | 0 | 0 |
| $L_e$ | 0 | 1 | 0 |
| $L_n$ | 0 | 0 | 1 |

DNL and CNEL weight an evening operation as one flown 19:00–22:00, so DNL counts it as day, the
AEDT and INM convention; by the clock, 22:00–23:00 is night in both. The entries are each
metric's weighting at 13:00, 21:00 and 03:00. $L_{\mathrm{Aeq},T}$ has no periods to count. In
§10.3 (number above, time above, coverage and provenance), $N_j = N_{d,j} + N_{e,j} + N_{n,j}$.

### 10.5 Effective perceived noise level

EPNL follows 14 CFR 36, Appendix A, §A36.4 [17], which is identical to ICAO Annex 16
[55]. It is ported term for term from PANAM.

**The records.** The one-third-octave spectra are formed as linear averages over 0.5 s records
in reception time (§A36.3). Each state's band energy is spread evenly over its reception-time
duration $T(1 + \dot r_i/c)$ and split among the records it overlaps in proportion:

$$
L_k = 10\lg\left(\sum_i 10^{L_i/10}\, \frac{\text{overlap}_{i,k}}{0.5\ \mathrm{s}}\right).
$$

Record $k$ covers $[0.5k, 0.5(k+1))$ s on the trajectory's clock, so together the records hold
the whole history's exposure in every band. The first and last records are averaged over the
whole 0.5 s. All 24 bands from 50 Hz to 10 kHz are required, each matched to the database band
whose center is within 5 percent of the nominal one; bands outside that range do not enter
EPNL. For each record $k$ and band $i$:

1. Noy values $n(i,k)$ from the formulas of Table A36-3 (§A36.4.7.3): four regions bounded by
   SPL(a), SPL(b), SPL(e) and SPL(d), and 0 noy below SPL(d).
2. $\mathrm{PNL}(k) = 40 + \dfrac{10}{\lg 2}\lg\left(0.85\, n_\max + 0.15 \sum_i n(i,k)\right)$.
3. The tone correction $C(k)$ of §A36.4.3.1, from the ten-step slope procedure: slopes from
   band 4 (80 Hz has none); a slope that changes by more than 5 dB is encircled, and with it
   SPL($i$) when the slope is positive and rising, or SPL($i-1$) when it has turned from
   positive to zero or negative; encircled levels are replaced by their neighbors' mean (band
   24 by SPL(23) + s(23)); the slopes are recomputed with $s'(3) = s'(4)$ and an imaginary 25th
   band; three adjacent slopes are averaged; the background SPL″ is built from
   SPL″(3) = SPL(3); and $C$ is the largest Table A36-2 factor of
   $F = \mathrm{SPL} - \mathrm{SPL}''$ over $F \ge 1.5$ dB: $F/3 - 1/2$, $F/6$ and $10/3$ below
   500 Hz and above 5 kHz, and $2F/3 - 1$, $F/3$ and $20/3$ from 500 Hz to 5 kHz inclusive.
   Bands 1–3 never carry a tone correction. A NaN band makes PNL and $C$ NaN, and an EPNL
   history holding one is refused.
4. $\mathrm{PNLT}(k) = \mathrm{PNL}(k) + C(k)$ and its maximum PNLTM.
5. The band-sharing check of §A36.4.4.2: when $C$ at PNLTM is below its mean $\bar C$ over the
   five records centered there (over the records that exist, at either end of the history),
   $\Delta_B = \mathrm{PNL}(k_M) + \bar C - \mathrm{PNLTM}$, added to EPNL as a separate term
   as ICAO Annex 16 [55] and FAA AC 36-4 [56] carry it; otherwise $\Delta_B = 0$. It is kept
   separate because, read literally, 14 CFR 36 raises PNLTM by $\Delta_B$ and takes $D$ from the
   raised PNLTM, so the two cancel in $\mathrm{PNLTM} + D$: the adjustment could never raise
   EPNL, only lower it through a narrower duration.
6. The duration correction $D$, from the unadjusted PNLTM, over the records $k_1 \ldots k_2$
   found as follows. $k_1$ and $k_2$ start at the first and last records with
   $\mathrm{PNLT} \ge \mathrm{PNLTM} - 10$ (the longest duration when there are several peaks),
   and each moves one record outward when that record is strictly closer to
   $\mathrm{PNLTM} - 10$. Every record between them is summed, dips included (§A36.4.5.5, as
   PANAM's F026 reads it). With $T_0 = 10$ s and $\Delta t = 0.5$ s:

$$
D = 10\lg\sum_{k=k_1}^{k_2} 10^{\mathrm{PNLT}(k)/10} + 10\lg\frac{\Delta t}{T_0} - \mathrm{PNLTM}, \qquad
\mathrm{EPNL} = \mathrm{PNLTM} + D + \Delta_B ,
$$

where $10\lg(\Delta t/T_0)$ is taken as exactly −13 for $\Delta t = 0.5$ s, the constant
§A36.4.5.4 states, rather than −13.0103.
Nothing checks that the history holds the whole 10 dB-down event beyond flagging a duration
that reaches the first or last record (the `clipped` column of `--epnl_file`), whose event
may have been cut off.

**Departures.** A band of $-\infty$ (zero energy), for which the tone correction is undefined,
enters it at its noy threshold SPL(d) and carries no tone, as in PANAM. By default
(`--epnl_masking`), every band below its noy threshold is masked: it enters the procedure at
that threshold and carries no tone correction itself, following the exclusion of tones in
masked bands in FAA AC 36-4 [56]. Without masking, deeply absorbed high bands would read as
tones. PNL is unchanged by masking, since a masked band is 0 noy either way. The corrections the
regulation applies to a measured spectrum (instrumentation, background noise, reference
conditions, pseudo-tones from ground reflections) are not applied to a prediction, so every
tone in the predicted audible spectrum counts, ground-effect ripple included.

### 10.6 C-weighted sound exposure level

C-SEL is the exposure integral of §10.2 over the C-weighted level, formed from the same band
levels at the microphone as the A-weighted level, over the same `--band_limits`, with the same
reception-time weights. Its 10 dB-down window is found on the C-weighted history itself, so it
is not the A-weighted window; `--full_history` sums every sample. The ambient of `-a` is an
A-weighted level with no spectrum, so it does not enter, and a run combining the two is refused.
C-weighting gives the low bands their full weight, so `--band_limits` matters more for C-SEL
than for SEL.

## 11. What is not modeled

- **Doppler frequency shift.** The spheres carry the measured, Doppler-shifted spectra, and
  NICE-OPS applies no shift. A database whose spheres have had the shift removed is refused. The
  source's motion enters only through the reception-time factor $1 + \dot r/c$ (§10.2).
- **Diffraction and screening** by terrain or buildings (§6.3). Terrain enters only through
  each path's mean ground plane.
- **Refraction over terrain.** The ray models assume flat ground, so `--terrain` is refused
  with `--refraction`. A propagation table is read over each path's mean plane (§9.3), which
  for a refracting table assumes a stratification that follows the ground.
- **Ground reflection on the broadband path.** The broadband path has only spreading and
  EAA (§4.1).
- **Coherent summation of refracted arrivals.** Arrivals are summed in energy (§8.2).
- **Wave effects inside a footprint run** beyond the tabulated shadow transition (§8.1, §8.3),
  the opt-in Airy corrections (§8.2) and a propagation table (§9.3). The PE itself never runs
  in a footprint.
- **Turbulence and roughness at a ground-plane microphone** (§6.2), and turbulent scattering into
  a shadow unless `--scattering` is given (§8.1).
- **Coupling between wind and the source condition** beyond the speed reference (§2.1).

## References

1. Greenwood, E., "Fundamental Rotorcraft Acoustic Modeling from Experiments (FRAME)," Ph.D. Dissertation, Dept. of Aerospace Engineering, University of Maryland, College Park, MD, 2011. http://hdl.handle.net/1903/11518
2. Greenwood, E., and Rau, R., "A Maneuvering Flight Noise Model for Helicopter Mission Planning," *Journal of the American Helicopter Society*, Vol. 65, No. 2, 2020, Paper 022007. https://doi.org/10.4050/JAHS.65.022007
3. Schmitz, F. H., Greenwood, E., Sickenberger, R. D., Gopalan, G., Sim, B. W.-C., Conner, D., Moralez, E., III, and Decker, W. A., "Measurement and Characterization of Helicopter Noise in Steady-State and Maneuvering Flight," *Proceedings of the American Helicopter Society 63rd Annual Forum*, Virginia Beach, VA, May 2007. https://doi.org/10.4050/VFS-F63-000258
4. Greenwood, E., Schmitz, F. H., and Sickenberger, R. D., "A Semiempirical Noise Modeling Method for Helicopter Maneuvering Flight Operations," *Journal of the American Helicopter Society*, Vol. 60, No. 2, 2015, Paper 022007. https://doi.org/10.4050/JAHS.60.022007
5. Watts, M. E., Greenwood, E., Smith, C. D., and Stephenson, J. H., "Noise Abatement Flight Test Data Report," NASA TM-2019-220264, NASA Langley Research Center, Hampton, VA, March 2019.
6. Pascioni, K. A., Greenwood, E., Watts, M. E., Smith, C. D., and Stephenson, J. H., "Medium-Sized Helicopter Noise Abatement Flight Test," *Proceedings of the Vertical Flight Society 76th Annual Forum*, Virtual, October 2020. https://doi.org/10.4050/F-0076-2020-16497
7. Bowyer, A., "Computing Dirichlet Tessellations," *The Computer Journal*, Vol. 24, No. 2, 1981, pp. 162–166. https://doi.org/10.1093/comjnl/24.2.162
8. Watson, D. F., "Computing the n-Dimensional Delaunay Tessellation with Application to Voronoi Polytopes," *The Computer Journal*, Vol. 24, No. 2, 1981, pp. 167–172. https://doi.org/10.1093/comjnl/24.2.167
9. Shewchuk, J. R., "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates," *Discrete & Computational Geometry*, Vol. 18, No. 3, 1997, pp. 305–363. https://doi.org/10.1007/PL00009321
10. Knuth, D. E., *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, Reading, MA, 1997.
11. "Report on Standard Method of Computing Rotorcraft Noise Contours," ECAC.CEAC Doc 32, 1st ed., European Civil Aviation Conference, Neuilly-sur-Seine, France, May 2026.
12. Page, J. A., Rapoza, A., Oberg, A., Hastings, A., Baker, G., and Shumway, M., "Advanced Acoustic Model (AAM) Technical Reference and User's Guide," Software Manual, Version 3.1, U.S. Department of Transportation, Volpe National Transportation Systems Center, Cambridge, MA, May 2021.
13. Reinsch, C. H., "Smoothing by Spline Functions," *Numerische Mathematik*, Vol. 10, No. 3, 1967, pp. 177–183. https://doi.org/10.1007/BF02162161
14. Snyder, J. P., "Map Projections — A Working Manual," U.S. Geological Survey Professional Paper 1395, U.S. Government Printing Office, Washington, DC, 1987. https://doi.org/10.3133/pp1395
15. "Department of Defense World Geodetic System 1984: Its Definition and Relationships with Local Geodetic Systems," NIMA TR8350.2, 3rd ed., Amendment 1, National Imagery and Mapping Agency, July 1997 (amended January 2000).
16. "Acoustics — Attenuation of Sound during Propagation Outdoors — Part 1: Calculation of the Absorption of Sound by the Atmosphere," ISO 9613-1:1993, International Organization for Standardization, Geneva, 1993.
17. "Noise Standards: Aircraft Type and Airworthiness Certification," Title 14, Code of Federal Regulations, Part 36, Appendix A, "Aircraft Noise Measurement and Evaluation Under § 36.101." URL: https://www.ecfr.gov/current/title-14/part-36 [retrieved 8 October 2026].
18. Zellmann, C., and Wunderli, J. M., "Influence of the Atmospheric Stratification on the Sound Propagation of Single Flights," *Proceedings of INTER-NOISE 2014*, Vol. 249, No. 4, Institute of Noise Control Engineering, Melbourne, Australia, 2014, pp. 3770–3779.
19. Stusnick, E., Plotkin, K. J., and Sutherland, L. C., "Short-Range Acoustic Propagation Model," Wyle Laboratories Research Report WR 85-19, Wyle Laboratories, Arlington, VA, July 1985.
20. Delany, M. E., and Bazley, E. N., "Acoustical Properties of Fibrous Absorbent Materials," *Applied Acoustics*, Vol. 3, No. 2, 1970, pp. 105–116. https://doi.org/10.1016/0003-682X(70)90031-9
21. Miki, Y., "Acoustical Properties of Porous Materials — Modifications of Delany-Bazley Models," *Journal of the Acoustical Society of Japan (E)*, Vol. 11, No. 1, 1990, pp. 19–24. https://doi.org/10.1250/ast.11.19
22. Attenborough, K., and Taherzadeh, S., "Prediction of Outdoor Ground Effect," *Applied Acoustics*, Vol. 242, 2026, Paper 111067. https://doi.org/10.1016/j.apacoust.2025.111067
23. Chien, C. F., and Soroka, W. W., "Sound Propagation Along an Impedance Plane," *Journal of Sound and Vibration*, Vol. 43, No. 1, 1975, pp. 9–20. https://doi.org/10.1016/0022-460X(75)90200-X
24. Chien, C. F., and Soroka, W. W., "A Note on the Calculation of Sound Propagation Along an Impedance Surface," *Journal of Sound and Vibration*, Vol. 69, No. 2, 1980, pp. 340–343. https://doi.org/10.1016/0022-460X(80)90618-5
25. Weideman, J. A. C., "Computation of the Complex Error Function," *SIAM Journal on Numerical Analysis*, Vol. 31, No. 5, 1994, pp. 1497–1518. https://doi.org/10.1137/0731077
26. Ament, W. S., "Toward a Theory of Reflection by a Rough Surface," *Proceedings of the IRE*, Vol. 41, No. 1, 1953, pp. 142–146. https://doi.org/10.1109/JRPROC.1953.274171
27. Chessell, C. I., "Propagation of Noise Along a Finite Impedance Boundary," *Journal of the Acoustical Society of America*, Vol. 62, No. 4, 1977, pp. 825–834. https://doi.org/10.1121/1.381603
28. Daigle, G. A., Embleton, T. F. W., and Piercy, J. E., "Some Comments on the Literature of Propagation near Boundaries of Finite Acoustical Impedance," *Journal of the Acoustical Society of America*, Vol. 66, No. 3, 1979, pp. 918–919. https://doi.org/10.1121/1.383207
29. van Maercke, D., and Defrance, J., "Development of an Analytical Model for Outdoor Sound Propagation Within the Harmonoise Project," *Acta Acustica united with Acustica*, Vol. 93, No. 2, 2007, pp. 201–212.
30. Ostashev, V. E., and Wilson, D. K., *Acoustics in Moving Inhomogeneous Media*, 2nd ed., CRC Press, Boca Raton, FL, 2015. https://doi.org/10.1201/b18922
31. "Ground-Plane Microphone Configuration for Propeller-Driven Light-Aircraft Noise Measurement," SAE ARP4055A, SAE International, Warrendale, PA, 2020.
32. "Commission Directive (EU) 2015/996 of 19 May 2015 Establishing Common Noise Assessment Methods According to Directive 2002/49/EC of the European Parliament and of the Council," *Official Journal of the European Union*, L 168, 1 July 2015, pp. 1–823.
33. Kephalopoulos, S., Paviotti, M., and Anfosso-Lédée, F., "Common Noise Assessment Methods in Europe (CNOSSOS-EU)," JRC Reference Report EUR 25379 EN, Publications Office of the European Union, Luxembourg, 2012. https://doi.org/10.2788/31776
34. Olsen, H., Tuinstra, M., and van Oosten, N., "D1.5d Rotorcraft Noise Modelling Guidance," Research Project NOISE, EASA.2020.FC.06, Royal Netherlands Aerospace Centre (NLR), for the European Union Aviation Safety Agency, Cologne, January 2024.
35. Salomons, E. M., *Computational Atmospheric Acoustics*, Kluwer Academic Publishers, Dordrecht, 2001. https://doi.org/10.1007/978-94-010-0660-6
36. Stull, R. B., *An Introduction to Boundary Layer Meteorology*, Kluwer Academic Publishers, Dordrecht, 1988. https://doi.org/10.1007/978-94-009-3027-8
37. Wieringa, J., "Updating the Davenport Roughness Classification," *Journal of Wind Engineering and Industrial Aerodynamics*, Vol. 41, Nos. 1–3, 1992, pp. 357–368. https://doi.org/10.1016/0167-6105(92)90434-C
38. "Guide to Instruments and Methods of Observation," WMO-No. 8, 2024 ed., World Meteorological Organization, Geneva, 2024.
39. "Erste Allgemeine Verwaltungsvorschrift zum Bundes-Immissionsschutzgesetz (Technische Anleitung zur Reinhaltung der Luft — TA Luft)," 27 February 1986, *Gemeinsames Ministerialblatt*, No. 7, 28 February 1986, p. 95.
40. "Meteorological Monitoring Programs for Nuclear Power Plants," Regulatory Guide 1.23, Rev. 1, U.S. Nuclear Regulatory Commission, Office of Nuclear Regulatory Research, Washington, DC, March 2007.
41. Salomons, E., van Maercke, D., Defrance, J., and de Roo, F., "The Harmonoise Sound Propagation Model," *Acta Acustica united with Acustica*, Vol. 97, No. 1, 2011, pp. 62–74. https://doi.org/10.3813/AAA.918387
42. Plovsing, B., "Nord2000. Comprehensive Outdoor Sound Propagation Model. Part 2: Propagation in an Atmosphere with Refraction," DELTA Acoustics & Vibration Report AV 1851/00 (revised), Hørsholm, Denmark, 2006.
43. Pierce, A. D., *Acoustics: An Introduction to Its Physical Principles and Applications*, 3rd ed., ASA Press/Springer, Cham, Switzerland, 2019. https://doi.org/10.1007/978-3-030-11214-1
44. Fock, V. A., *Electromagnetic Diffraction and Propagation Problems*, International Series of Monographs on Electromagnetic Waves, Vol. 1, Pergamon Press, Oxford, 1965.
45. van Maercke, D., and Salomons, E., "The Harmonoise Sound Propagation Model: Further Developments and Comparison with Other Models," *Proceedings of NAG/DAGA 2009 International Conference on Acoustics*, Deutsche Gesellschaft für Akustik, Berlin, 2009, pp. 1622–1625.
46. Collins, M. D., "A Split-Step Padé Solution for the Parabolic Equation Method," *Journal of the Acoustical Society of America*, Vol. 93, No. 4, 1993, pp. 1736–1742. https://doi.org/10.1121/1.406739
47. Collins, M. D., "A Self-Starter for the Parabolic Equation Method," *Journal of the Acoustical Society of America*, Vol. 92, No. 4, 1992, pp. 2069–2074. https://doi.org/10.1121/1.405258
48. Gilbert, K. E., Raspet, R., and Di, X., "Calculation of Turbulence Effects in an Upward-Refracting Atmosphere," *Journal of the Acoustical Society of America*, Vol. 87, No. 6, 1990, pp. 2428–2437. https://doi.org/10.1121/1.399088
49. Attenborough, K., Taherzadeh, S., Bass, H. E., Di, X., Raspet, R., Becker, G. R., Güdesen, A., Chrestman, A., Daigle, G. A., L'Espérance, A., Gabillet, Y., Gilbert, K. E., Li, Y. L., White, M. J., Naz, P., Noble, J. M., and van Hoof, H. A. J. M., "Benchmark Cases for Outdoor Sound Propagation Models," *Journal of the Acoustical Society of America*, Vol. 97, No. 1, 1995, pp. 173–191. https://doi.org/10.1121/1.412302
50. Olver, F. W. J., Olde Daalhuis, A. B., Lozier, D. W., et al. (eds.), "NIST Digital Library of Mathematical Functions," Release 1.2.8, National Institute of Standards and Technology, 15 September 2026. URL: https://dlmf.nist.gov/ [retrieved 8 October 2026].
51. "Electroacoustics — Sound Level Meters — Part 1: Specifications," IEC 61672-1:2013, Ed. 2.0, International Electrotechnical Commission, Geneva, 2013.
52. "Airport Noise Compatibility Planning," Title 14, Code of Federal Regulations, Part 150, Appendix A, "Noise Exposure Maps." URL: https://www.ecfr.gov/current/title-14/part-150 [retrieved 8 October 2026].
53. "Directive 2002/49/EC of the European Parliament and of the Council of 25 June 2002 Relating to the Assessment and Management of Environmental Noise," *Official Journal of the European Communities*, L 189, 18 July 2002, pp. 12–25.
54. "Noise Standards: Definitions," California Code of Regulations, Title 21, Sec. 5001.
55. "Environmental Protection, Volume I — Aircraft Noise," Annex 16 to the Convention on International Civil Aviation, 8th ed., International Civil Aviation Organization, Montréal, July 2017.
56. "Noise Standards: Aircraft Type and Airworthiness Certification," Advisory Circular AC 36-4D, Federal Aviation Administration, Washington, DC, 12 October 2017.
57. "Report on Standard Method of Computing Noise Contours around Civil Airports, Volume 2: Technical Guide," ECAC.CEAC Doc 29, 4th ed., European Civil Aviation Conference, Neuilly-sur-Seine, France, December 2016.
58. "Procedure for the Calculation of Airplane Noise in the Vicinity of Airports," SAE AIR1845A, SAE International, Warrendale, PA, 2012.
59. Plovsing, B., "Nord2000. Comprehensive Outdoor Sound Propagation Model. Part 1: Propagation in an Atmosphere without Significant Refraction," DELTA Acoustics & Vibration Report AV 1849/00 (revised), Hørsholm, Denmark, 2006.
