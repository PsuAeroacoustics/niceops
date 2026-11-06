// The first-order near-field correction of a tone (THEORY 4.4): the sphere
// Laplacian it takes, and the corrected field against the exact field of a
// rotating source.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <functional>
#include <vector>

#include "sphere_laplacian.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

using cplx = std::complex<double>;

struct Grid
{
    std::vector<double> phis, thetas;
};

// phi over [-180, 180) (a full sphere) or [-90, 90] (a lower hemisphere), theta
// over [0, 180], both at `step` degrees.
Grid make_grid(double step, bool full)
{
    Grid g;
    for (double p = full ? -180.0 : -90.0; p < (full ? 180.0 - 1e-9 : 90.0 + 1e-9); p += step)
        g.phis.push_back(p);
    for (double t = 0.0; t <= 180.0 + 1e-9; t += step)
        g.thetas.push_back(t);
    return g;
}

// max |Lap f - lambda f| over max |lambda f|, f a degree-l spherical harmonic
// (eigenvalue -l (l + 1)) given as a function of (forward, starboard, up).
double eigen_error(const Grid& g, const std::function<cplx(double, double, double)>& f, int l, int order)
{
    std::vector<cplx> values;
    for (const double p : g.phis)
        for (const double t : g.thetas)
        {
            const auto d = grid_direction(p, t);
            values.push_back(f(d.forward, d.starboard, d.up));
        }
    const auto lap = sphere_laplacian(g.phis, g.thetas, values, order);
    const double lambda = -l * (l + 1.0);
    double worst = 0.0, scale = 0.0;
    std::size_t at = 0;
    for (std::size_t k = 0; k < values.size(); k++)
    {
        if (std::abs(lap[k] - lambda * values[k]) > worst)
            at = k;
        worst = std::max(worst, std::abs(lap[k] - lambda * values[k]));
        scale = std::max(scale, std::abs(lambda * values[k]));
    }
    if (std::getenv("NICEOPS_NEAR_FIELD_DIAG"))
        WARN("relative error " << worst / scale << ", worst at phi " << g.phis[at / g.thetas.size()] << " theta "
                << g.thetas[at % g.thetas.size()]);
    return worst / scale;
}

} // namespace

TEST_CASE("sphere_laplacian: spherical harmonics are its eigenfunctions", "[near_field]")
{
    // Degree-l harmonics (eigenvalue -l (l + 1)) of the direction (x forward, y
    // starboard, z up). (x + i y)^m winds as e^{i m alpha} about the shaft;
    // (x + i y)^m ((2m + 3) z^2 - 1) is degree m + 2; both are even across the
    // horizontal plane, as a lower hemisphere is read. (y + i z)^m + (y - i z)^m
    // winds about the nose, through the grid's poles. The relative error (max over
    // the grid of |Lap f + l(l+1) f| over max l(l+1)|f|) falls as h^2 in the
    // interior: 15 deg 0.009-0.27 for l = 1-8, 5 deg 0.0036-0.034, 2.5 deg
    // 0.0009-0.0085.
    struct Case
    {
        const char* name;
        std::function<cplx(double, double, double)> f;
        int degree, order;
    };
    for (const bool full : {true, false})
    {
        for (const double step : {15.0, 5.0, 2.5})
        {
            const Grid g = make_grid(step, full);
            for (int m = 1; m <= 6; m++)
            {
                const std::vector<Case> cases {
                    {"shaft", [m](double x, double y, double) { return std::pow(cplx(x, y), m); }, m, m},
                    {"shaft, counter", [m](double x, double y, double) { return std::pow(cplx(x, -y), m); }, m, -m},
                    {"shaft, degree m + 2",
                            [m](double x, double y, double z) {
                                return std::pow(cplx(x, y), m) * ((2.0 * m + 3.0) * z * z - 1.0);
                            },
                            m + 2, m},
                    {"nose", [m](double, double y, double z) { return std::pow(cplx(y, z), m) + std::pow(cplx(y, -z), m); },
                            m, 0},
                };
                for (const auto& c : cases)
                {
                    // A tone wound m times per turn about another axis than the one
                    // unwound: on a coarse grid only for the slowest windings.
                    if (c.order == 0 && step > 10.0 && m > 2)
                        continue;
                    const double error = eigen_error(g, c.f, c.degree, c.order);
                    // Second order in the interior, about 0.06 (l h)^2 (h the step, rad),
                    // first order at the fitted nodes (poles, near the shaft).
                    const double h = step * pi / 180.0;
                    const double tolerance = 0.25 * (c.degree * h) * (c.degree * h) + 0.1 * h;
                    INFO(c.name << ", step " << step << (full ? ", full sphere" : ", lower hemisphere") << ", m " << m
                                << ": relative error " << error);
                    CHECK(error < tolerance);
                }
            }
        }
    }
}

namespace {

// A rotating source of compact monopoles and dipoles (thrust and drag), as
// wopwop3's docs/noise-sphere-near-field/wilcox_check.py builds it: B blades of
// Nr stations each, one revolution sampled at Ntau, flight Mach M along the
// forward axis, the medium at rest. Coordinates (forward, starboard, up); time
// convention e^{i omega t}.
struct Rotor
{
    static constexpr double c = 340.29;
    double omega = 0.0;
    double mach = 0.0;
    std::vector<std::array<double, 3>> y, f;
    std::vector<double> q, tau;
    double w = 0.0;

    Rotor(double mach_flight, std::array<double, 3> axis, double radius = 6.58, int blades = 3, double tip_mach = 0.6,
            int stations = 6, int samples = 256)
        : mach(mach_flight)
    {
        omega = tip_mach * c / radius;
        const double na = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        for (auto& v : axis)
            v /= na;
        auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
            return std::array<double, 3> {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
        };
        auto e1 = cross(axis, std::array<double, 3> {0.0, 1.0, 0.0});
        const double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        for (auto& v : e1)
            v /= n1;
        const auto e2 = cross(axis, e1);
        const double period = 2.0 * pi / omega;
        w = 1.0 / samples;
        for (int b = 0; b < blades; b++)
            for (int j = 0; j < stations; j++)
            {
                const double r = radius * (0.2 + 0.8 * (j + 0.5) / stations);
                const double weight = (r / radius) * (r / radius);
                for (int s = 0; s < samples; s++)
                {
                    const double t = s * period / samples;
                    const double psi = omega * t + 2.0 * pi * b / blades;
                    const double cp = std::cos(psi), sp = std::sin(psi);
                    std::array<double, 3> pos {}, force {};
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = r * (cp * e1[i] + sp * e2[i]);
                        const double along = -sp * e1[i] + cp * e2[i];
                        force[i] = weight * axis[i] - 0.15 * weight * along;
                    }
                    y.push_back(pos);
                    f.push_back(force);
                    q.push_back(r / radius);
                    tau.push_back(t);
                }
            }
    }

    // The r -> infinity directivity A(n): R e^{ikR} p as R grows.
    cplx far(const GridDirection& d, int m) const
    {
        const double k = m * omega / c;
        const double dop = 1.0 - mach * d.forward;
        cplx sum(0.0, 0.0);
        for (std::size_t s = 0; s < y.size(); s++)
        {
            const double yn = y[s][0] * d.forward + y[s][1] * d.starboard + y[s][2] * d.up;
            const double fn = f[s][0] * d.forward + f[s][1] * d.starboard + f[s][2] * d.up;
            const cplx source = q[s] + cplx(0.0, k * fn / dop);
            sum += source * std::polar(1.0, k * yn / dop - k * c * tau[s]);
        }
        return sum * w / (4.0 * pi * dop);
    }

    // R e^{ikR} p at the emission-geometry observer x = R (n - M), from the exact
    // convected Green's function.
    cplx at_radius(const GridDirection& d, double R, int m) const
    {
        const double k = m * omega / c;
        const double b2 = 1.0 - mach * mach;
        const std::array<double, 3> x {R * (d.forward - mach), R * d.starboard, R * d.up};
        cplx sum(0.0, 0.0);
        for (std::size_t s = 0; s < y.size(); s++)
        {
            const double r0 = x[0] - y[s][0], r1 = x[1] - y[s][1], r2 = x[2] - y[s][2];
            const double perp2 = r1 * r1 + r2 * r2;
            const double rs = std::sqrt(r0 * r0 + b2 * perp2);
            const double phase = (rs + mach * r0) / b2;
            const cplx g = std::polar(1.0, -k * phase) / (4.0 * pi * rs);
            // grad Rs and grad Phi
            const double grs[3] = {r0 / rs, b2 * r1 / rs, b2 * r2 / rs};
            const double gph[3] = {(grs[0] + mach) / b2, grs[1] / b2, grs[2] / b2};
            cplx fdotg(0.0, 0.0);
            for (int i = 0; i < 3; i++)
                fdotg += f[s][i] * g * (cplx(0.0, -k * gph[i]) - grs[i] / rs);
            sum += (q[s] * g - fdotg) * std::polar(1.0, -k * c * tau[s]);
        }
        return R * std::polar(1.0, k * R) * sum * w;
    }
};

} // namespace

TEST_CASE("The near-field correction against the exact field of a rotating source", "[near_field]")
{
    // A + (i / (2 k R)) Lap[D A] against the exact field at R = 2 to 16 c / Omega:
    // the far field's error falls as 1 / R, the corrected one's as 1 / R^2 (until
    // the grid's Laplacian error, which is of order (l h)^2 of the correction).
    const Grid g = make_grid(5.0, true);
    struct Flight
    {
        const char* name;
        double mach;
        std::array<double, 3> axis;
    };
    for (const auto& flight : {Flight {"hover", 0.0, {0.0, 0.0, 1.0}}, Flight {"M 0.2", 0.2, {-0.1, 0.0, 1.0}}})
    {
        const Rotor rotor(flight.mach, flight.axis);
        const double c_over_omega = Rotor::c / rotor.omega;
        for (const int m : {3, 6})
        {
            std::vector<cplx> a, da;
            std::vector<GridDirection> dirs;
            for (const double p : g.phis)
                for (const double t : g.thetas)
                {
                    const auto d = grid_direction(p, t);
                    dirs.push_back(d);
                    a.push_back(rotor.far(d, m));
                    da.push_back((1.0 - flight.mach * d.forward) * a.back());
                }
            // The rotor turns from aft toward port seen from above (clockwise), so
            // its tones wind as e^{-i m alpha}: Database's order sign(shaft rate) m.
            const auto lap = sphere_laplacian(g.phis, g.thetas, da, -m);
            double peak = 0.0;
            for (const auto& v : a)
                peak = std::max(peak, std::abs(v));
            const double k = m * rotor.omega / Rotor::c;
            std::vector<double> far_error, corrected_error;
            for (const double r_over : {2.0, 4.0, 8.0, 16.0})
            {
                const double R = r_over * c_over_omega;
                double e0 = 0.0, e1 = 0.0;
                for (std::size_t cell = 0; cell < a.size(); cell += 7)
                {
                    if (std::abs(a[cell]) < 0.1 * peak)
                        continue;
                    const cplx exact = rotor.at_radius(dirs[cell], R, m);
                    e0 = std::max(e0, std::abs(exact - a[cell]) / peak);
                    e1 = std::max(e1, std::abs(exact - a[cell] - cplx(0.0, 1.0 / (2.0 * k * R)) * lap[cell]) / peak);
                }
                far_error.push_back(e0);
                corrected_error.push_back(e1);
                if (std::getenv("NICEOPS_NEAR_FIELD_DIAG"))
                    WARN(flight.name << ", m " << m << ", R " << r_over << " c/Omega: far-field error " << e0
                                     << ", corrected " << e1);
            }
            INFO(flight.name << ", m " << m);
            // 1 / R without the correction, close to 1 / R^2 with it (R = 4 to 16
            // c / Omega): hover, m 3, 0.49 to 0.12 against 0.096 to 0.0064.
            CHECK(far_error[1] / far_error[3] > 3.5);
            CHECK(far_error[1] / far_error[3] < 4.5);
            CHECK(corrected_error[1] / corrected_error[3] > 10.0);
            // Where the expansion holds, (m + 1) / (2 R Omega / c) <= 0.5, the
            // correction removes most of the error; it never adds to it.
            const double r_over[] = {2.0, 4.0, 8.0, 16.0};
            for (std::size_t i = 0; i < far_error.size(); i++)
            {
                CHECK(corrected_error[i] < far_error[i]);
                if ((m + 1.0) / (2.0 * r_over[i]) <= 0.5)
                    CHECK(corrected_error[i] < 0.3 * far_error[i]);
            }
        }
    }
}
