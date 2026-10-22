#include "profile_pe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "exact_solutions.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// The grid and the domain, in wavelengths: points a wavelength in height and
// range steps a wavelength; the physical domain this far above the highest of
// the roofs and the ends, then an absorbing layer this thick; the Huygens line
// this far before the last edge, at most half the last roof.  The results
// move by under 0.03 dB rms over the boundary-element cases for reaches of 15
// to 40, layers of 15 to 50 and lines 2 to 8 wavelengths back.
constexpr double points_per_wavelength = 20.0;
constexpr double steps_per_wavelength = 10.0;
constexpr double reach_wavelengths = 15.0;
constexpr double layer_wavelengths = 15.0;
constexpr double huygens_wavelengths = 4.0;

// The 2-D Green's function, (i/4) H0(k r).
complex green(double k, double r)
{
    return complex(0.0, 0.25) * hankel1_0(complex(k * r, 0.0));
}

// A rigid right-angled corner's Green's function by UTD: the corner at c,
// its face 0 straight down (the back wall) and its face n = 3/2 along -u (the
// roof).  The direct wave and each face's reflection where they reach R, and
// the corner's diffraction.
struct CornerGreen
{
    double k, cu, cz, ru, rz;
    static constexpr double n = 1.5;

    // Angles from face 0, counterclockwise: down 0, +u pi/2, up pi, -u 3pi/2.
    double angle(double u, double z) const
    {
        const double a = std::atan2(u - cu, -(z - cz));
        return a < 0.0 ? a + 2.0 * pi : a;
    }
    // R's angle and Q's.  On a shadow or reflection boundary the waves' test
    // and the coefficient's limit must take the same side; R within a hair of
    // one is moved off it, for both, as hard_wedge_field does.  A Huygens
    // node lands on one exactly at some wavenumber (its height is a whole
    // number of steps of a twentieth of a wavelength), and was 30 dB loud.
    std::array<double, 2> angles(double qu, double qz) const
    {
        const double th0 = angle(qu, qz);
        double th = angle(ru, rz);
        for (double boundary : {th0 + pi, th0 - pi, pi - th0, (2.0 * n - 1.0) * pi - th0})
            if (std::abs(th - boundary) < 1e-9)
                th = boundary + 1e-8;
        return {th, th0};
    }
    // Which geometrical-optics waves reach R from Q: direct, off face 0, off
    // face n.
    std::array<bool, 3> lit(double qu, double qz) const
    {
        const auto [th, th0] = angles(qu, qz);
        return {std::abs(th - th0) < pi, th + th0 < pi, th + th0 > (2.0 * n - 1.0) * pi};
    }
    // G(Q, R) with the waves `lit` says reach R from Q.
    complex operator()(double qu, double qz, const std::array<bool, 3>& waves) const
    {
        const auto [th, th0] = angles(qu, qz);
        const double rho0 = std::hypot(qu - cu, qz - cz), rho = std::hypot(ru - cu, rz - cz);
        complex out(0.0, 0.0);
        if (waves[0])
            out += green(k, std::hypot(ru - qu, rz - qz));
        if (waves[1])
            out += green(k, std::sqrt(rho * rho + rho0 * rho0 - 2.0 * rho * rho0 * std::cos(th + th0)));
        if (waves[2])
            out += green(k, std::sqrt(rho * rho + rho0 * rho0 - 2.0 * rho * rho0 * std::cos(2.0 * n * pi - th0 - th)));
        const double L = rho * rho0 / (rho + rho0);
        out += green(k, rho0) * hard_wedge_coefficient(n, th, th0, k, L) * std::polar(1.0, k * rho) / std::sqrt(rho);
        return out;
    }
};

// The system (1 + alpha q) x = rhs on the nodes from g up, rigid at g (its
// ghost the node above), factored once for every right-hand side: the march
// takes the same one at every step over a roof or a street's floor, and both
// its fields at each.
class Factored
{
public:
    void factor(const std::vector<complex>& a, std::size_t g, double c, complex alpha)
    {
        const std::size_t m = a.size() - g;
        g_ = g;
        alpha_ = alpha;
        upper_.assign(m, alpha * c);
        upper_[0] = alpha * 2.0 * c;
        upper_[m - 1] = 0.0;
        multiplier_.resize(m);
        inverse_.resize(m);
        complex pivot = 1.0 + alpha * (-2.0 * c + a[g]);
        inverse_[0] = 1.0 / pivot;
        multiplier_[0] = 0.0;
        for (std::size_t j = 1; j < m; j++)
        {
            multiplier_[j] = alpha * c / pivot;
            pivot = 1.0 + alpha * (-2.0 * c + a[g + j]) - multiplier_[j] * upper_[j - 1];
            inverse_[j] = 1.0 / pivot;
        }
    }
    bool is(std::size_t g, complex alpha) const
    {
        return !inverse_.empty() && g == g_ && alpha == alpha_;
    }
    // x in place of rhs.
    void solve(std::vector<complex>& rhs) const
    {
        const std::size_t m = inverse_.size();
        for (std::size_t j = 1; j < m; j++)
            rhs[j] -= multiplier_[j] * rhs[j - 1];
        rhs[m - 1] *= inverse_[m - 1];
        for (std::size_t j = m - 1; j-- > 0;)
            rhs[j] = (rhs[j] - upper_[j] * rhs[j + 1]) * inverse_[j];
    }

private:
    std::size_t g_ = 0;
    complex alpha_;
    std::vector<complex> upper_, multiplier_, inverse_;
};

} // namespace

bool profile_applies(const ProfileProblem& p)
{
    if (p.blocks.empty())
        return false;
    for (std::size_t i = 0; i < p.blocks.size(); i++)
    {
        const ProfileBlock& b = p.blocks[i];
        if (!(b.u1 > b.u0) || !(b.h > 0.0) || (i > 0 && b.u0 < p.blocks[i - 1].u1))
            return false;
    }
    return p.su < p.blocks.front().u0 && p.ru > p.blocks.back().u1 && p.sz >= 0.0 && p.rz >= 0.0;
}

ProfileParts profile_parts(const ProfileProblem& p, double k)
{
    if (!profile_applies(p))
        throw std::invalid_argument("profile_parts: the blocks must be ascending and apart, the source before them "
                                    "and the receiver behind");
    const double lambda = 2.0 * pi / k;
    const double dz = lambda / points_per_wavelength;
    double top = std::max(p.sz, p.rz);
    for (const ProfileBlock& b : p.blocks)
        top = std::max(top, b.h);
    const double z_phys = top + reach_wavelengths * lambda;
    const double z_top = z_phys + layer_wavelengths * lambda;
    const std::size_t n = static_cast<std::size_t>(std::ceil(z_top / dz)) + 1;
    std::vector<double> z(n);
    std::vector<complex> a(n, complex(0.0, 0.0));
    for (std::size_t i = 0; i < n; i++)
    {
        z[i] = static_cast<double>(i) * dz;
        if (z[i] > z_phys)
        {
            const double s = (z[i] - z_phys) / (z_top - z_phys);
            a[i] = complex(0.0, s * s);
        }
    }
    const ProfileBlock& last = p.blocks.back();
    const double u_start = p.blocks.front().u0;
    const double huygens = std::min(huygens_wavelengths * lambda, 0.5 * (last.u1 - last.u0));
    const double u_stop = last.u1 - huygens;
    auto ground = [&](double u) {
        for (const ProfileBlock& b : p.blocks)
            if (u >= b.u0 && u <= b.u1)
                return b.h;
        return 0.0;
    };
    auto node = [&](double h) { return static_cast<std::size_t>(std::lround(h / dz)); };

    // The starts: the source's wave on the first block's front face, above
    // its roof, and its floor image's.
    std::vector<complex> psi[2] = {std::vector<complex>(n, 0.0), std::vector<complex>(n, 0.0)};
    const std::size_t g0 = node(p.blocks.front().h);
    const complex back = std::polar(1.0, -k * u_start);
    for (std::size_t i = g0; i < n; i++)
    {
        psi[0][i] = green(k, std::hypot(u_start - p.su, z[i] - p.sz)) * back;
        psi[1][i] = green(k, std::hypot(u_start - p.su, z[i] + p.sz)) * back;
    }

    const double c = 1.0 / ((k * dz) * (k * dz));
    // q psi on the active nodes from g, rigid at g (its ghost the node above).
    auto apply_q = [&](std::size_t g, const std::vector<complex>& f, std::vector<complex>& out) {
        const std::size_t m = n - g;
        out.resize(m);
        for (std::size_t j = 0; j < m; j++)
        {
            const std::size_t i = g + j;
            const complex below = j == 0 ? f[i + 1] : f[i - 1];
            const complex above = j + 1 < m ? f[i + 1] : complex(0.0, 0.0);
            out[j] = c * (below - 2.0 * f[i] + above) + a[i] * f[i];
        }
    };
    Factored left;
    double u = u_start;
    const double du0 = lambda / steps_per_wavelength;
    std::vector<complex> qpsi, rhs;
    while (u < u_stop - 1e-12 * lambda)
    {
        const double du = std::min(du0, u_stop - u);
        const std::size_t g = node(ground(u + 0.5 * du));
        const complex alpha_l = (1.0 - complex(0.0, k * du)) / 4.0, alpha_r = (1.0 + complex(0.0, k * du)) / 4.0;
        if (!left.is(g, alpha_l))
            left.factor(a, g, c, alpha_l);
        for (std::vector<complex>& f : psi)
        {
            // Below the ground the field is gone (a wall stepped up into) or
            // not yet there (a street's floor stepped down to).
            std::fill(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(g), complex(0.0, 0.0));
            apply_q(g, f, qpsi);
            rhs.resize(n - g);
            for (std::size_t j = 0; j < n - g; j++)
                rhs[j] = f[g + j] + alpha_r * qpsi[j];
            left.solve(rhs);
            std::copy(rhs.begin(), rhs.end(), f.begin() + static_cast<std::ptrdiff_t>(g));
        }
        u += du;
    }

    // The Huygens line: p and dp/du = i k Q psi exp(i k u), Q = (1 + 3q/4) / (1 + q/4),
    // to the receiver and to its image.
    const std::size_t g = node(last.h);
    const std::size_t m = n - g;
    Factored quarter;
    quarter.factor(a, g, c, complex(0.25, 0.0));
    std::vector<complex> dpsi[2];
    for (int f = 0; f < 2; f++)
    {
        apply_q(g, psi[f], qpsi);
        dpsi[f].resize(m);
        for (std::size_t j = 0; j < m; j++)
            dpsi[f][j] = psi[f][g + j] + 0.75 * qpsi[j];
        quarter.solve(dpsi[f]);
    }
    const complex phase = std::polar(1.0, k * u);
    const double dl = 1e-3 * lambda;
    ProfileParts out {};
    const double extent = z_phys - last.h, t0 = 0.5 * extent;
    for (int side = 0; side < 2; side++)
    {
        const CornerGreen green_w {k, last.u1, last.h, p.ru, side == 0 ? p.rz : -p.rz};
        for (std::size_t j = 0; j < m; j++)
        {
            const double zj = z[g + j], t = zj - last.h;
            if (t > extent)
                break;
            double w = j == 0 ? 0.5 * dz : dz;
            if (t > t0)
                w *= 0.5 * (1.0 + std::cos(pi * (t - t0) / (extent - t0)));
            // dG/du by a central difference, or a one-sided one where a
            // shadow boundary passes between the neighbours: there the
            // diffraction cancels the wave's step only to O(1 / k rho), and
            // the remainder over 2 dl was a spike in the receiver's level.
            const std::array<bool, 3> waves = green_w.lit(u, zj), plus = green_w.lit(u + dl, zj),
                                      minus = green_w.lit(u - dl, zj);
            const complex gq = green_w(u, zj, waves);
            complex dgdu;
            if (plus == waves && minus == waves)
                dgdu = (green_w(u + dl, zj, plus) - green_w(u - dl, zj, minus)) / (2.0 * dl);
            else if (plus == waves)
                dgdu = (green_w(u + dl, zj, plus) - gq) / dl;
            else
                dgdu = (gq - green_w(u - dl, zj, minus)) / dl;
            for (int f = 0; f < 2; f++)
            {
                const complex pu = psi[f][g + j] * phase;
                const complex dpu = complex(0.0, k) * dpsi[f][j] * phase;
                out.part[f][side] += w * (-gq * dpu + pu * dgdu);
            }
        }
    }
    // From (i/4) H0(k r) ~ exp(i k r) / sqrt(r) times sqrt(1 / (8 pi k)) exp(i pi / 4).
    const complex scale = std::sqrt(8.0 * pi * k) * std::polar(1.0, -0.25 * pi);
    for (auto& row : out.part)
        for (complex& v : row)
            v *= scale;
    return out;
}

complex profile_field(const ProfileProblem& p, double k)
{
    const ProfileParts parts = profile_parts(p, k);
    return parts.part[0][0] + p.floor * (parts.part[1][0] + parts.part[0][1]) + p.floor * p.floor * parts.part[1][1];
}

std::vector<ProfileParts> profile_spectrum(const ProfileProblem& problem, const std::vector<double>& wavenumbers)
{
    std::vector<ProfileParts> out;
    out.reserve(wavenumbers.size());
    for (double k : wavenumbers)
        out.push_back(profile_parts(problem, k));
    return out;
}

double profile_path(const ProfileProblem& p, bool source_image, bool receiver_image,
        std::vector<std::array<double, 2>>* points)
{
    // The upper hull of the ends and the roofs' corners.
    std::vector<std::array<double, 2>> hull {{p.su, source_image ? -p.sz : p.sz}};
    auto add = [&](std::array<double, 2> q) {
        while (hull.size() >= 2)
        {
            const auto& o = hull[hull.size() - 2];
            const auto& b = hull.back();
            if ((b[0] - o[0]) * (q[1] - o[1]) - (b[1] - o[1]) * (q[0] - o[0]) < 0.0)
                break;
            hull.pop_back();
        }
        hull.push_back(q);
    };
    for (const ProfileBlock& b : p.blocks)
    {
        add({b.u0, b.h});
        add({b.u1, b.h});
    }
    add({p.ru, receiver_image ? -p.rz : p.rz});
    double length = 0.0;
    for (std::size_t i = 1; i < hull.size(); i++)
        length += std::hypot(hull[i][0] - hull[i - 1][0], hull[i][1] - hull[i - 1][1]);
    if (points)
        *points = std::move(hull);
    return length;
}

} // namespace niceops
