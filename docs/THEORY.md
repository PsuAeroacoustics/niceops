# NICE-OPS theory

This document states the models NICE-OPS implements, as equations, with the references they
come from. It describes the code as it is: where NICE-OPS departs from a cited reference, the
departure is stated. [The manual](MANUAL.md) covers the options that select each model, the
evidence for the defaults, and how each model was checked. [`internals/`](internals/) covers
the implementation.

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
source directivity from a database of measured noise "spheres" [1].
For each observer it carries the level along a path through the atmosphere to the ground. It then
integrates each observer's level history into exposure metrics. The approach follows the
maneuvering-flight noise model of Greenwood and Rau [2]. Every observer is evaluated
against every state. The demo case has 6,551 states and 14,400 observers.

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
quasi-static acoustic mapping (Q-SAM) of Schmitz, Greenwood and co-workers
[3, 1, 4]. The load-factor vector $\mathbf{n}$ is the aircraft's
acceleration less gravity, in units of $g = 9.80665$ m/s²:

$$
\mathbf{n} = \frac{\mathbf{a}}{g} + \hat{\mathbf{z}}, \qquad n = \lvert\mathbf{n}\rvert, \qquad \hat{\mathbf{n}} = \mathbf{n}/n .
$$

$\hat{\mathbf{n}}$ stands in for the main rotor's tip-path-plane normal. With the velocity
$\mathbf{u}$ (over the ground by default, through the air under `--speed_reference air`), the
tip speed $V_\mathrm{tip}$, rotor radius $R$, weight $W$ and air density $\rho$ from the vehicle
file:

$$
\mu = \frac{\lvert\mathbf{u}\rvert}{V_\mathrm{tip}}, \qquad
C_T = \frac{n\,W}{\rho\,V_\mathrm{tip}^2\,\pi R^2}, \qquad
\gamma = \arcsin\frac{\mathbf{u}\cdot\hat{\mathbf{n}}}{\lvert\mathbf{u}\rvert}\ \ (\text{Q-SAM, default}), \qquad
\gamma = \arcsin\frac{u_z}{\lvert\mathbf{u}\rvert}\ \ (\texttt{--no-qsam}).
$$

The Q-SAM path angle is the velocity's angle to the tip-path plane. It is the angle at which a
steady flight condition presents the same rotor inflow, which is what lets steady-flight spheres
represent a maneuvering state. $\rho$ is the vehicle file's density in every case, including
when `--temperature` or `--pressure` sets the propagation atmosphere.

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
triangulated by an incremental Delaunay triangulation [7, 8]. The
triangulation is built on exact orientation and in-sphere predicates [9], with the
error-free sum of [10]. A state's normalized condition $\mathbf{q}$, clamped to the
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

A sphere that has no data in a direction (stored at or below −999 dB) is dropped from that
direction's blend, and the remaining weights are renormalized. A query outside the convex hull
takes the nearest condition.

Two alternatives exist for comparison with other models. `--blend nearest` takes the nearest
condition in the normalized space, as the one-hemisphere-per-segment practice of ECAC Doc 32
[11] would at every state. `--blend aam` is the bracketing search of the Advanced Acoustic
Model, v3.1 [12, §2.5.1.1]. It groups the spheres into rows of one flight path angle,
rounded to 0.5°. It interpolates between the bracketing rows only if they are at most 3° apart
on approach or 5° apart on departure. Within a row, it interpolates between bracketing speeds at
most 20 kt apart. Otherwise it takes the nearest row or speed. The result is up to four spheres,
blended linearly in dB, speed first.

**Load factor.** A `fixed_load_factor` database stores only the 1 g spheres and is triangulated
over $(\mu, \gamma)$. The load factor enters through the thrust law

$$
\Delta L_n = 20 \lg \frac{C_T}{C_{T,\mathrm{ref}}}, \qquad C_T > C_{T,\mathrm{ref}},
$$

added to every level and band. $C_T$ is capped at `--max_load_factor` (default 2) times
$C_{T,\mathrm{ref}}$. Below 1 g the law is not applied by default. `--load_factor_correction
two_sided` applies it down to `--min_load_factor` (default 0.5).

### 2.4 Reading the sphere

Each sphere is fixed to the rotor. Its axes are $\hat{\mathbf{n}}$, the nose $\hat{\mathbf{e}}$
(the velocity projected into the tip-path plane) and starboard
$\hat{\mathbf{s}} = \hat{\mathbf{e}} \times \hat{\mathbf{n}}$. Below 8 kt the nose turns toward
the aircraft's heading, with a smoothstep weight
$w = s^2(3 - 2s)$, $s = (\lvert u_h\rvert - 3\ \mathrm{kt})/(5\ \mathrm{kt})$. For the unit vector
$\hat{\mathbf{d}}$ from source to observer:

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
from the D original. The blend in §2.3 is applied at each of
the four grid corners first.

## 3. Trajectories

**AC state tracks** (the 27-column NASA/FAA format) carry positions, velocities and
accelerations, and are used as sampled.

**WGS84 tracks** carry positions only. Velocities and accelerations come from a cubic smoothing
spline [13] that minimizes

$$
\sum_i \left(y_i - f(t_i)\right)^2 + \lambda \int f''(t)^2\, dt
$$

per axis. $\lambda$ is chosen so that the residual equals a tolerance: `--path_smoothing`
(default 10 ft rms per axis) times the number of points. Dense tracks are padded at each end with
10 s of points on a quadratic fitted to the last 3 s, so the spline's natural end conditions do
not bend the ends of the track.

**Coordinates.** Latitude and longitude are mapped to a local plane by the American polyconic
projection about a reference point [14] on the WGS84 ellipsoid [15]
($a = 6\,378\,137$ m, $f = 1/298.257223563$). Altitude passes straight through. The plane is
then rotated so that $+x$ points along `--frame_bearing` $b$:

$$
x = E \sin b + N \cos b, \qquad y = -E \cos b + N \sin b .
$$

The same rotation turns compass bearings (wind directions, headings) into the local frame.

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
[17], and the approach of the sonAIR model [18]. On a curved ray the mean is still
taken over the straight span $[z_o, z_s]$, a first-order approximation. Pressure in the profile,
when not given, is hydrostatic for dry air:
$p(z) = p(0) \exp\!\left(-\frac{g}{R_d}\int_0^z \frac{dz'}{T}\right)$, with
$R_d = 287.05$ J/(kg K).

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

In the variable-porosity model $\sigma_e$ is the effective flow resistivity in Pa·s/m²,
$\alpha_e$ the rate of change of porosity with depth, $\gamma = 1.4$ and
$\rho_0 = 1.2$ kg/m³. For the layer, $Z_c = 1 + 9.08X^{-0.754} + 11.9iX^{-0.732}$ and
$k_c = k\,(1 + 10.8X^{-0.70} + 10.3iX^{-0.595})$. The models use $\beta = 1/Z$, and a ground with
$\mathrm{Re}\,\beta \le 0$ is refused as non-passive.

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

**The Faddeeva function** is computed by Weideman's rational series with $N = 32$ terms
[25]:

$$
\mathrm{w}(z) = \frac{2 P(Z)}{(L - iz)^2} + \frac{1}{\sqrt{\pi}(L - iz)}, \qquad
Z = \frac{L + iz}{L - iz}, \qquad L = \sqrt{N/\sqrt{2}},
$$

with $P$ a polynomial of degree 31. The lower half-plane uses
$\mathrm{w}(z) = 2e^{-z^2} - \mathrm{w}(-z)$. For speed, the hot path uses the asymptotic series
for $\lvert z\rvert \ge 6$ and a sixth-order Taylor expansion about tabulated nodes elsewhere.
Both agree with the series to about $10^{-11}$, and the series agrees with SciPy's `wofz` to
$3\times10^{-13}$.

**Roughness.** An rms roughness height $\sigma_r$ (`--roughness`, default 0) reduces the
coherent reflection by the Kirchhoff (Ament) factor [26]:
$Q \leftarrow Q\, e^{-2(k\sigma_r\cos\vartheta)^2}$.

### 5.3 Band averaging and turbulence

A one-third-octave band is not a pure tone. NICE-OPS uses Chessell's band-averaged form
[27], with the corrections noted by Daigle, Embleton and Piercy [28], with the path delay $\tau = (R_2 - R_1)/c$ and range ratio
$\varrho = R_2/R_1$:

$$
G = 1 + \frac{\lvert Q\rvert^2}{\varrho^2} +
  \frac{2}{\varrho}\left[\mathrm{Re}\,Q\cos(\Omega f\tau) - \mathrm{Im}\,Q\sin(\Omega f\tau)\right]
    T(f)\, \frac{\sin(\Delta f\tau)}{\Delta f\tau}, \qquad \Delta L_g = 10\lg G,
$$

with $\Omega = \pi(2^{1/6} + 2^{-1/6}) = 6.325159$ and $\Delta = \pi(2^{1/6} - 2^{-1/6}) = 0.727477$.
These are the interference evaluated at the band's arithmetic center, averaged over its width.
Over a rigid surface with the receiver on the ground, $G \to 4$ (+6 dB).

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

## 6. Receivers and terrain

### 6.1 The point receiver

The default receiver is a point `--receiver_height` above the ground (1.2 m, as 14 CFR 36
prescribes [17]), or each microphone's own height from an observer file. It is §5's model
with $h_r$ set to that height.

### 6.2 Ground-plane microphones

`--plate_table` models a microphone mounted in a circular plate on the ground: the GRAS 67AX in a
GR1425 plate, 0.2 m in radius, with the microphone 3/4 of the radius from the center, as in
SAE ARP 4055 [31]. The plate's response is tabulated by PANAM's axisymmetric boundary
element method as complex transfer functions of the direct and image waves, $P_d$ and $P_r$,
over frequency, elevation and azimuth. The band gain averages five sub-frequencies $f_j$ across
the band:

$$
G_b = \frac15 \sum_{j=1}^{5} \left\lvert P_d(f_j) + Q(f_j)\, P_r(f_j) \right\rvert^2 ,
$$

with $Q$ the reflection coefficient of §5.2 evaluated at the plate's top surface. The table
records the ground it was computed over, and a run with a different ground is refused.

### 6.3 Terrain

With `--terrain`, each source–observer pair takes its own mean ground plane: a least-squares
line $z = a + b s$ fitted to the terrain heights sampled along the pair's horizontal path. This
is the mean-plane construction of CNOSSOS-EU [32, 33], also used over
terrain in NORAH2 [34] and ECAC Doc 32 [11]. The ground model then runs in the rotated frame of
that plane:

$$
h_s' = \frac{z_s - (a + bD)}{\sqrt{1 + b^2}}, \qquad
h_r' = \frac{z_m - a}{\sqrt{1 + b^2}}, \qquad
d' = \frac{D + b(z_s - z_m)}{\sqrt{1 + b^2}},
$$

with $D$ the horizontal distance and $z_m$ the microphone's height. A negative equivalent height
is raised to zero, as in CNOSSOS-EU. Spreading, absorption and the sphere direction still use
the true straight line. Terrain screening and diffraction are not modeled.

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
[36]

$$
u(z) = u_\mathrm{ref}\,\frac{\ln(z/z_0)}{\ln(z_\mathrm{ref}/z_0)}, \qquad z_0 < z \le \delta,
$$

held constant above the surface-layer depth $\delta$ (default 200 m) and zero at or below $z_0$.
The roughness length $z_0$ defaults to 0.03 m, the Davenport–Wieringa class "open flat terrain"
[37, 38]. The anemometer height $z_\mathrm{ref}$ defaults to the standard 10 m. The
temperature lapse is dry adiabatic (−0.0098 K/m), and water vapor holds its surface molar
concentration. `--stability` replaces these with a power law $u = u_\mathrm{ref}(z/z_\mathrm{ref})^m$ and a
class temperature gradient. The exponents are the Klug–Manier classes of TA Luft [39], and
the gradients are mid-range values of the $\Delta T$ classes of U.S. NRC Regulatory Guide 1.23
[40]. The stability scheme has not been validated.

## 8. Refraction

### 8.1 Linear sound-speed profile: circular rays

`--refraction linear` fits each source–observer pair's $c_\mathrm{eff}(z)$ over the pair's
height span (at least 10 m) with a straight line $c(z) = c_0 + a z$ by least squares. This is the
equivalent-linear approach of HARMONOISE and Nord2000 [41, 42]. The rays are
then arcs of circles of radius $\mathcal{L} = c_0/\lvert a\rvert$ centered at height
$-c_0/a$. For the direct ray, with $\kappa(z) = z(2\mathcal{L} + \sigma z)$ and
$\sigma = \operatorname{sign} a$, the arc's center lies a horizontal distance
$x_c = [r^2 + \sigma(\kappa_r - \kappa_s)]/(2r)$ from the source. The ray-tube spreading range
is

$$
R_\mathrm{eq}^2 = r \left\lvert \frac{\partial r}{\partial \theta_s} \right\rvert \frac{\sin\lvert\theta_r\rvert}{\cos\theta_s}
 \;\;\Rightarrow\;\; R_\mathrm{eq} = \frac{r}{\cos\theta_s}\ \text{on a circle},
$$

and the travel time is
$t = \operatorname{artanh}\left(r\rho_c/(h_c^2 + x_c r)\right)/\lvert a\rvert$, with
$h_c = \mathcal{L} + \sigma z_s$ and $\rho_c = \sqrt{x_c^2 + h_c^2}$.
The level is read from the sphere at the ray's launch direction, and absorbed over the arc
length. The reflected ray is found by iteration on its launch angle, and the ground model is
evaluated at the equivalent straight geometry with the ray's grazing angle $\psi$ and path length.

**Shadow zone.** In an upward-refracting atmosphere ($a < 0$) no ray reaches beyond

$$
r_b = t(z_s) + t(z_r), \qquad t(z) = \sqrt{z(2\mathcal{L} - z)} \approx \sqrt{2\mathcal{L} z}.
$$

Beyond it, sound arrives as a creeping wave diffracted along the ground. It decays at the
rate of the first mode [43, §9-5; 44]. Over rigid ground, with $\tau_1 = 1.01879$ the
first zero of $\mathrm{Ai}'(-x)$:

$$
\alpha_\mathrm{sh} = \frac{20}{\ln 10}\,\tau_1 \sin 60^\circ\, \frac{(\pi f)^{1/3} \lvert a\rvert^{2/3}}{c_0}\ \ \mathrm{dB/m}.
$$

Over an impedance ground, the root $b_1$ of $\mathrm{Ai}'(b) + q\,\mathrm{Ai}(b) = 0$, with
$q = -i s\, e^{4\pi i/3} \beta$ and $s = (\pi f/\lvert a\rvert)^{1/3}$, replaces it:
$\alpha_\mathrm{sh} = 8.686\, \mathrm{Im}(b_1 e^{-2\pi i/3})\, \pi f/(c_0 s^2)$. The loss
$\min(\alpha_\mathrm{sh}(r - r_b), L_\mathrm{floor})$, with `--shadow_floor` default 20 dB,
is subtracted.

**The shadow transition.** The ray picture is wrong near the boundary, where the field
diffracts smoothly from lit to shadow. Near the ground the equivalent-linear method can depart
from a PE by up to 10 dB [45]. NICE-OPS tabulates the difference between the exact
field of a linear profile over an impedance plane (§9.2) and the ray model, as a function of the
Fock-scaled distance past the boundary $\xi = (r - r_b)/F$, with $F = 2s^2/k_0$, and of
$s = k_0 \ell$, $\ell = (c_0/(2\lvert a\rvert k_0^2))^{1/3}$. Near the boundary it blends the ray
model into that correction. Where a stratified profile places each band's boundary is a fit to
the parabolic equation particular to NICE-OPS (see the manual's
[Refraction](MANUAL.md#refraction)).

**Turbulent scattering** (`--scattering`, off by default). Turbulence scatters energy into the
shadow. With `--scattering` $\gamma_T$, the shadow loss is limited by a scattered level

$$
S = \min\!\left(0,\ a_0 + 10\lg\gamma_T + b_0 \lg\frac{f}{1\ \mathrm{kHz}} +
  c_0' \lg\frac{R}{100\ \mathrm{m}} + d_0 \lg\frac{z_s}{100\ \mathrm{m}}\right)
$$

relative to free field, added in energy to the coherent field. The coefficients are NICE-OPS's
own regression on parabolic-equation runs with von Kármán phase screens (§9.1), not a published
formula. For a flush receiver they are $(a_0, b_0, c_0', d_0) = (39.18, -1.55, 2.10, 5.30)$, and
for one 1.2 m up $(36.67, 5.61, 5.13, 5.30)$.

### 8.2 Stratified profile: ray tracing through layers

`--refraction stratified` traces rays through $c_\mathrm{eff}(z)$, taken as piecewise linear
between the profile's levels. Snell's law fixes the ray parameter $p = \cos\theta/c$. Across a
layer of gradient $g$ from $c_0$ to $c_1$, with $q = \sqrt{1 - p^2c^2}$, the ray advances

$$
\Delta x = \frac{p\,(c_0 + c_1)\,\Delta z}{q_0 + q_1}, \qquad
\Delta t = \frac{1}{g} \ln\frac{c_1 (1 + q_0)}{c_0 (1 + q_1)} .
$$

The eigenrays are the roots of $X(p) = r$ over four families: launched up or down, with or without
a ground bounce. Each is solved by the Illinois variant of regula falsi. The spreading of each
arrival comes from the ray tube, $R_\mathrm{eq}^2 = r\lvert dX/d\theta_s\rvert\sin\lvert\theta_r\rvert/\cos\theta_s$,
with any focusing gain capped at 10 dB. In a duct, rays that bounce $m$ times arrive at
$X_m = X_0 + m L_\mathrm{hop}$. Each intermediate bounce costs the plane-wave loss
$-10\lg\lvert R_p(\psi)\rvert^2$, and the last bounce takes the ground model.

Multiple arrivals are summed **incoherently**, $L = 10\lg\sum_i 10^{L_i/10}$. Near caustics and
for rays that tunnel through a turning point (`--tunneling`), the band level is the band average
of the Airy-function intensity [43]. Far into the shadow, a stratified pair hands over to
the equivalent linear model of §8.1, and its creeping wave and transition.

## 9. The parabolic equation and reference solutions

The parabolic equation (PE) and the analytic solutions do not run inside a footprint. They check
the refraction models (the `niceops_pe_compare` and `niceops_benchmark` tools), and they build
propagation-loss tables (§9.3).

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
by a rational function $r_0 + \sum_j c_j/(\mathcal{Q} - p_j)$. **Departure:** the coefficients are not the
Taylor–Padé coefficients. They come from a rational interpolant fitted at run time on Chebyshev
nodes over the propagation angles needed (up to 50°), with Collins's stability constraints, and
the fewest terms that meet a $10^{-7}$ tolerance. The depth operator is discretized to fourth
order (Numerov), at 10 points per wavelength, including an original fourth-order impedance row. A
Crank–Nicolson scheme, the Padé(1,1) case [35, App. F], is available. The field starts
from Collins's self-starter [47], or from the Gaussian starter of [35]. An
artificial absorbing layer 100 wavelengths thick closes the top of the domain.

**Turbulence** enters as random phase screens [48; 35, App. I], one every
outer scale $L_0 = 10$ m. Each screen's phase spectrum is that of a von Kármán spectrum with
structure parameter $\gamma_T/4$:

$$
S_\varphi(\kappa_z) = 2\pi k_a^2\, \Delta x\; 0.033\, \frac{\gamma_T}{4}\;
 \frac{\sqrt{\pi}\,\Gamma(4/3)}{\Gamma(11/6)}\;\left(\kappa_z^2 + L_0^{-2}\right)^{-4/3} .
$$

The relative level is $\Delta L = 20\lg(\lvert\psi\rvert R/\sqrt{r})$ relative to free field.
Band levels are energy means over five sub-frequencies and over the turbulence realizations.

### 9.2 Analytic reference solutions

These are the benchmark cases of Attenborough et al. [49]. In homogeneous air
the field is $p = e^{ikR_1}/R_1 + Q\,e^{ikR_2}/R_2$ (§5.2). In a linear profile,
$k^2(z) = k_0^2(1 \mp 2z/\mathcal{R})$, the field is a sum of normal modes,

$$
p = \pi i \sum_n \frac{u_n(z)\,u_n(z_s)}{N_n}\, H_0^{(1)}(\kappa_n r), \qquad N_n = \int_0^\infty u_n^2\, dz,
$$

with $\ell = (\mathcal{R}/(2k_0^2))^{1/3}$ and $\mathrm{w}_3 = e^{2\pi i/3}$:

- **Downward refraction:** $u_n = \mathrm{Ai}(t_n + z/\ell)$, with
  $\mathrm{Ai}'(t_n) + i k_0 \ell \beta\, \mathrm{Ai}(t_n) = 0$.
- **Upward refraction:** $u_n = \mathrm{Ai}(b_n - \mathrm{w}_3 z/\ell)$, with
  $\mathrm{Ai}'(b_n) - i k_0 \ell \beta\, \mathrm{w}_3^2\, \mathrm{Ai}(b_n) = 0$.

**Departure:** the upward-refraction mode condition printed in [49, Eq. 11],
$\mathrm{Ai}'(b) + i k_0 \ell \beta\,\mathrm{w}_3\,\mathrm{Ai}(b) = 0$, differs from the one derived
here by a factor $e^{i\pi/3}$ on the admittance term. The two agree over rigid ground. Over
impedance ground the printed form departs from the PE by 1.5–12 dB at 100 Hz, and the derived
form agrees with it to 0.02 dB. The Airy and Hankel functions follow the DLMF [50].

### 9.3 Propagation-loss tables

A `--propagation_table` stores $\Delta L(\text{bearing}, z_s, z_r, r, f)$ relative to free field
at the straight-line distance, computed by the PE through one profile. A run interpolates it
multilinearly in dB and applies it in place of the ground and refraction models:
$L_b = L_{\mathrm{sph},b} - 20\lg(r/R_s) - \alpha_b(r - R_s) + \Delta L$.

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

The level history is sampled at emission times $t_i$, a constant step $T$ apart, but exposure
accumulates in reception time $t_i + r_i/c$. Each sample therefore lasts
$T(1 + \dot r_i/c)$, compressed while the aircraft approaches and stretched while it recedes:

$$
\mathrm{SEL} = 10\lg\left( \frac{T}{t_\mathrm{ref}} \sum_{i=i_1}^{i_2} \left(1 + \frac{\dot r_i}{c}\right) 10^{L_i/10} \right),
\qquad t_\mathrm{ref} = 1\ \mathrm{s},
$$

with $\dot r_i$ a central difference. The sum runs over the 10 dB-down window, from the first to the
last sample within 10 dB of the peak, as in 14 CFR 36 [17, A36.4.5]. `--full_history` sums
every sample. Under refraction $r_i$ is $c$ times the ray's travel time. An ambient level
$L_\mathrm{amb}$ (`-a`) is added in energy, $10\lg(10^{L/10} + 10^{L_\mathrm{amb}/10})$, before
the window is found.

### 10.3 Maximum level, time above and number above

The signal is piecewise constant in mean square, $P_i = 10^{L_i/10}$, over the reception-time
durations above. The exponential time weighting of IEC 61672-1 [51] (fast, $\tau = 0.125$ s,
by default) is applied exactly for that input:

$$
y_i = P_i + (y_{i-1} - P_i)\, e^{-T(1 + \dot r_i/c)/\tau}, \qquad
L_{AF\mathrm{max}} = 10\lg\max_i y_i .
$$

Slow weighting ($\tau = 1$ s) and a 0.5 s running $L_\mathrm{eq}$, the maximum used by ECAC Doc 32
[11], are also available. A track sampled more coarsely than $\tau/2$ is refused. Time above a
threshold (TA) solves for the crossings of $y$ analytically within each sample. The number of
events above (NAx) counts the operations whose maximum reaches the threshold.

### 10.4 Day–night metrics

Each operation's exposure is computed with a penalty $P$ added to every sample heard in a
penalized period, by local time at emission. The operations are then summed and normalized to a
24-hour day:

$$
L_\mathrm{dn} = 10\lg\left(\frac{1}{86\,400\ \mathrm{s}} \sum_\mathrm{ops} t_\mathrm{ref}\, 10^{\mathrm{SEL}_P/10}\right).
$$

| Metric | Evening | Night | Definition |
| --- | --- | --- | --- |
| DNL | — | 22:00–07:00, +10 dB | 14 CFR Part 150 [52] |
| Lden | 19:00–23:00, +5 dB | 23:00–07:00, +10 dB | Directive 2002/49/EC [53] |
| CNEL | 19:00–22:00, +10 lg 3 dB | 22:00–07:00, +10 dB | Cal. Code Regs. tit. 21, §5001 [54] |

### 10.5 Effective perceived noise level

EPNL follows 14 CFR 36, Appendix A, §A36.4 [17], which is identical to ICAO Annex 16
[55]. It is ported term for term from PANAM. The one-third-octave spectra are formed as
linear averages over 0.5 s records in reception time (§A36.3). Each state's band energy is
spread over its reception-time duration and split among the records it overlaps. For each
record $k$:

1. Noy values $n(i,k)$ from the formulas of Table A36-3.
2. $\mathrm{PNL}(k) = 40 + \dfrac{10}{\lg 2}\lg\left(0.85\, n_\max + 0.15 \sum_i n(i,k)\right)$.
3. The tone correction $C(k)$ of §A36.4.3.1, from the 10-step slope procedure and Table A36-2.
4. $\mathrm{PNLT}(k) = \mathrm{PNL}(k) + C(k)$, its maximum PNLTM, and the band-sharing check of §A36.4.4.2.
5. The duration correction $D$, summed over the records within 10 dB of PNLTM, with
   $T_0 = 10$ s and $\Delta t = 0.5$ s, and EPNL:

$$
D = 10\lg\left[\frac{1}{T_0}\sum_k 10^{\mathrm{PNLT}(k)/10}\,\Delta t\right] - \mathrm{PNLTM}, \qquad
\mathrm{EPNL} = \mathrm{PNLTM} + D .
$$

**Departures.** $D$ is computed unrounded, where the regulation writes $10\lg(0.5/10)$ as −13.
By default, a band below its noy threshold is masked: it enters the procedure at that threshold
and carries no tone correction, following the exclusion of tones in masked bands in
FAA AC 36-4 [56]. Without masking, deeply absorbed high bands would read as tones. The
corrections the regulation applies to a measured spectrum (instrumentation, background noise,
reference conditions) are not applied to a prediction.

## 11. What is not modeled

- **Doppler frequency shift.** The spheres carry the measured, Doppler-shifted spectra, and
  NICE-OPS applies no shift. A database whose spheres have had the shift removed is refused. The
  source's motion enters only through the reception-time factor $1 + \dot r/c$ (§10.2).
- **Diffraction by terrain or buildings**, and screening (§6.3).
- **Ground reflection on the broadband path.** The broadband path has only spreading and
  EAA (§4.1).
- **Coherent summation of refracted arrivals.** Arrivals are summed in energy (§8.2).
- **Coupling between wind and the source condition** beyond `--speed_reference air` (§2.1).

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
