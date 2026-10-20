#pragma once

// The full-wave reference for a point source among rigid solids: the
// three-dimensional wave equation marched by a staggered pressure-velocity
// FDTD (Yee's scheme), the reference the ray model of point_paths.hpp is held
// to where no exact solution exists -- several buildings, stepped roofs,
// facade profiles (FIDELITY "Urban receivers").
//
// Pressure lives at the cells' centres and each velocity component on the
// faces normal to it.  A cell is solid or air, by whether its centre lies in
// a solid (voxels, so a box whose faces fall on the grid is exact and
// anything else is a staircase); a face with a solid on either side is rigid,
// its velocity held at zero, and so are the faces at a floor (z = lo.z) and at
// a symmetry plane (y = lo.y), which a mirror image makes exact.  Every other
// side of the domain is a convolutional PML (Roden and Gedney 2000), with
// auxiliary fields only in its slabs.  The solids run on into the PML, so a
// building cut by the domain's edge continues out of it.
//
// The source is a monopole of volume velocity Q(t), a Gaussian, injected into
// the eight cells around it by the weights that interpolate the pressure
// there; the receivers read the pressure by the same weights.  Its free field
// is rho Q'(t - r / c) / (4 pi r), so a receiver's transfer function re the
// source's field at 1 m is
//
//   H(f) = P(f) 4 pi / (-i omega rho Q(f)),
//
// |H| = 1/r in free field, the ray model's normalisation.  Yee's second-order
// differences are slow by (k h)^2 (1 - C^2) / 24 along the axes and exact
// along the diagonals at C = 1, an error that depends on direction; the
// fourth-order differences (the default) are fast by nearly the same amount
// in every direction, which the warp takes out (FdtdSetup::order).

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

// A rigid prism: a polygon in plan, x and y, less its holes, from `bottom` to
// `top`.
struct FdtdSolid
{
    std::vector<double> x, y;
    std::vector<Ring> holes;
    double bottom = 0.0, top = 0.0;
    // A pitched roof's planes, the solid below the least of them
    // (buildings.hpp); empty for a flat top.
    std::vector<RoofPlane> roof;
};

// The solids of LoD1 buildings on a floor at z = 0: each footprint to its
// roof (its holes are not solid), from far below the floor, so that the
// solids continue into a PML under a domain without one.
std::vector<FdtdSolid> fdtd_solids(const std::vector<Footprint>& footprints);

struct FdtdSetup
{
    double cell = 0.004;           // m
    // The air the domain must hold, before the PML: the box [lo, hi], m.  The
    // grid snaps lo down to a multiple of the cell, so that solids whose faces
    // lie on multiples of it are gridded exactly.
    Vec3 lo, hi;
    bool floor = true;             // a rigid floor at z = lo.z (else a PML below)
    bool mirror = false;           // a rigid plane at y = lo.y (else a PML there)
    int pml = 20;                  // cells
    double sound_speed = 343.2;    // m/s
    double density = 1.2;          // kg/m^3
    // The spatial differences' order: 2, Yee's, or 4, the staggered 9/8 and
    // -1/24 with each far value beyond a rigid face its mirror, which keeps
    // the march energy-conserving.  Fourth order's phase error is nearly the
    // same in every direction, so `warp` removes most of it: each frequency is
    // read where the grid carries the physical wavenumber
    // (FdtdRecording::transfer); second order's depends on direction, and it
    // is not warped.  Fourth order's sources and receivers are tricubic.
    int order = 4;
    bool warp = true;
    double courant = 0.99;         // of the scheme's limit, 1 / sqrt(3) for order 2, 6 / (7 sqrt(3)) for 4
    // The solids' absorption at normal incidence, as a real impedance
    // (1 + sqrt(1 - alpha)) / (1 - sqrt(1 - alpha)) times rho c, locally
    // reacting, the same at every frequency (as canyon_wave.hpp's walls).
    // The floor and the plane of symmetry stay rigid.
    double wall_absorption = 0.0;
    double duration = 0.03;        // s of simulated time
    // The pulse: rho Q' peaks in spectrum here, Hz, and is down 10 dB at
    // 2.25 times it and 30 dB at 3.25.
    double pulse_frequency = 4000.0;
    unsigned threads = 1;
    // March on the GPU (gpu/gpu.hpp), in float as the CPU's, where there is
    // one and it holds the grid, and on the CPU otherwise;
    // FdtdRecording::gpu says which.  With gpu_required the CPU is refused
    // instead.
    bool gpu = false;
    bool gpu_required = false;
    // March each step only the cells that matter: those the source can have
    // reached, and from which a receiver can still be reached before the
    // end, at the scheme's fastest group speed and a margin.  Off marches
    // every cell, which tests compare it with.
    bool cull = true;
};

struct FdtdRecording
{
    double dt = 0.0;               // s
    double cell = 0.0;             // m
    std::size_t steps = 0;
    std::size_t receivers = 0;
    // Q at the half steps (n + 1/2) dt, m^3/s, and the pressure at each
    // receiver at the steps (n + 1) dt, Pa: [n * receivers + receiver].
    std::vector<double> source;
    std::vector<float> pressure;
    double density = 1.2;
    double sound_speed = 343.2;
    double pulse_delay = 0.0;      // the time of Q's peak, s
    int order = 2;
    bool warp = true;

    // The frequency at which the grid carries the wavenumber 2 pi f / c,
    // averaged over directions: f itself without `warp`.
    double grid_frequency(double frequency) const;
    std::size_t cells = 0;         // nx ny nz
    std::size_t nx = 0, ny = 0, nz = 0;
    double seconds = 0.0;          // the march's wall time
    bool gpu = false;              // marched on the GPU
    std::string gpu_declined;      // why not, where the GPU was wanted and there is one
    double marched = 1.0;          // the share of the cell-steps marched (cull)

    // H at frequency f for a receiver (above), from the first `used` steps
    // (0 for all of them), read at grid_frequency(f).  With `attenuation`,
    // the air's, Np/m of amplitude: every sample of the pressure and of Q by
    // exp(-a c (t - t_peak)), which is the transfer at a complex frequency
    // where every arrival of path L carries exp(-a L), as the ray model
    // absorbs its rays.
    std::complex<double> transfer(std::size_t receiver, double frequency, std::size_t used = 0,
            double attenuation = 0.0) const;
    // |H|^2 averaged over [f1, f2], on a grid of `per_band` frequencies or one
    // every quarter of 1 / duration, whichever is finer.
    double mean_energy(std::size_t receiver, double f1, double f2, std::size_t used = 0, std::size_t per_band = 64,
            double attenuation = 0.0) const;
    // The band's level, dB re the source at 1 m, as band_level() takes the
    // ray model's: 10 lg of |H|^2 averaged over the third octave.  With
    // `attenuation`, the air's in each third of the band in log frequency
    // (sub_band_absorption's three).
    double band_level(std::size_t receiver, double band, std::size_t used = 0,
            const double* attenuation = nullptr) const;
    // How much the band's level moves without the recording's last tenth, dB:
    // whether the run lasted long enough for it.
    double truncation_db(std::size_t receiver, double band) const;
};

// The largest group speed of the march's waves, over every wavevector, re c,
// at order 2 or 4 and nu = c dt / h.  Second order's is 1, at k -> 0.
// Fourth order's leapfrog outruns c: 1.23 c at nu 0.49, near the Nyquist
// wavenumber along a diagonal, and 1.01 c along an axis at 7.5 points a
// wavelength.  FdtdSetup::cull's cones grow at it.
double fdtd_fastest_group_speed(int order, double nu);

// The march.  A source or receiver within half a cell of the floor or the
// symmetry plane is taken on it, in the first cell's centre (the field is
// even there); a source on the symmetry plane radiates into half the space
// with half its strength, as its image completes it.  A source off the plane
// with `mirror` is refused, and so is one in a solid.
FdtdRecording solve_fdtd(const FdtdSetup& setup, const std::vector<FdtdSolid>& solids, const Vec3& source,
        const std::vector<Vec3>& receivers);

} // namespace niceops
