#include "sphere_laplacian.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "util/constants.hpp"

using namespace niceops; // NOLINT

GridDirection niceops::grid_direction(double phi_deg, double theta_deg)
{
    const double phi = phi_deg * (pi / 180.0), theta = theta_deg * (pi / 180.0);
    return {std::cos(theta), std::sin(theta) * std::sin(phi), -std::sin(theta) * std::cos(phi)};
}

namespace {

using cplx = std::complex<double>;

struct Unit
{
    double x, y, z;
};

double dot(const Unit& a, const Unit& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Unit unit_of(const GridDirection& d)
{
    return {d.forward, d.starboard, d.up};
}

// The azimuth about the up axis toward starboard (Database's shaft_azimuth).
double azimuth(const Unit& n)
{
    return std::atan2(n.y, n.x);
}

// An orthonormal tangent basis at p.
std::array<Unit, 2> tangent_basis(const Unit& p)
{
    // Any axis not along p.
    const Unit a = std::fabs(p.z) < 0.9 ? Unit {0.0, 0.0, 1.0} : Unit {1.0, 0.0, 0.0};
    Unit e1 {a.x - dot(a, p) * p.x, a.y - dot(a, p) * p.y, a.z - dot(a, p) * p.z};
    const double n1 = std::sqrt(dot(e1, e1));
    e1 = {e1.x / n1, e1.y / n1, e1.z / n1};
    const Unit e2 {p.y * e1.z - p.z * e1.y, p.z * e1.x - p.x * e1.z, p.x * e1.y - p.y * e1.x};
    return {e1, e2};
}

// Azimuthal-equidistant coordinates of q in the tangent plane at p (basis e).
std::array<double, 2> tangent_coordinates(const Unit& p, const std::array<Unit, 2>& e, const Unit& q)
{
    const double c = std::clamp(dot(p, q), -1.0, 1.0);
    const double d = std::acos(c);
    const Unit t {q.x - c * p.x, q.y - c * p.y, q.z - c * p.z};
    const double nt = std::sqrt(dot(t, t));
    if (!(nt > 0.0))
        return {0.0, 0.0};
    return {d * dot(t, e[0]) / nt, d * dot(t, e[1]) / nt};
}

// Least-squares quadratic g = c0 + c1 u + c2 v + c3 u^2 + c4 u v + c5 v^2 through
// the samples: value, gradient (c1, c2) and Laplacian 2 (c3 + c5) at the origin.
struct LocalFit
{
    cplx value, du, dv, laplacian;
};

LocalFit quadratic_fit(const std::vector<std::array<double, 2>>& uv, const std::vector<cplx>& y)
{
    std::array<std::array<double, 6>, 6> a {};
    std::array<cplx, 6> b {};
    for (std::size_t s = 0; s < uv.size(); s++)
    {
        const double u = uv[s][0], v = uv[s][1];
        const std::array<double, 6> x {1.0, u, v, u * u, u * v, v * v};
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
                a[i][j] += x[i] * x[j];
            b[i] += x[i] * y[s];
        }
    }
    // Gaussian elimination with partial pivoting.
    std::array<int, 6> row {0, 1, 2, 3, 4, 5};
    for (int col = 0; col < 6; col++)
    {
        int best = col;
        for (int r = col + 1; r < 6; r++)
            if (std::fabs(a[row[r]][col]) > std::fabs(a[row[best]][col]))
                best = r;
        std::swap(row[col], row[best]);
        const double pivot = a[row[col]][col];
        assert(pivot != 0.0);
        for (int r = col + 1; r < 6; r++)
        {
            const double factor = a[row[r]][col] / pivot;
            for (int j = col; j < 6; j++)
                a[row[r]][j] -= factor * a[row[col]][j];
            b[row[r]] -= factor * b[row[col]];
        }
    }
    std::array<cplx, 6> c {};
    for (int col = 5; col >= 0; col--)
    {
        cplx sum = b[row[col]];
        for (int j = col + 1; j < 6; j++)
            sum -= a[row[col]][j] * c[j];
        c[col] = sum / a[row[col]][col];
    }
    return {c[0], c[1], c[2], 2.0 * (c[3] + c[5])};
}

// Three-point first and second derivatives on a non-uniform stencil.
void three_point(double h1, double h2, cplx fm, cplx f0, cplx fp, cplx& d1, cplx& d2)
{
    d1 = (-h2 / (h1 * (h1 + h2))) * fm + ((h2 - h1) / (h1 * h2)) * f0 + (h1 / (h2 * (h1 + h2))) * fp;
    d2 = 2.0 * (fm / (h1 * (h1 + h2)) - f0 / (h1 * h2) + fp / (h2 * (h1 + h2)));
}

} // namespace

std::vector<cplx> niceops::sphere_laplacian(const std::vector<double>& phis_deg, const std::vector<double>& thetas_deg,
        const std::vector<cplx>& f, int order)
{
    const std::size_t n_phi = phis_deg.size(), n_theta = thetas_deg.size(), n_cells = n_phi * n_theta;
    assert(f.size() == n_cells);
    std::vector<cplx> out(n_cells, cplx(0.0, 0.0));
    if (n_phi < 3 || n_theta < 3)
        return out;
    const double deg = pi / 180.0;
    double widest_phi = 0.0, widest_theta = 0.0;
    for (std::size_t m = 1; m < n_phi; m++)
        widest_phi = std::max(widest_phi, phis_deg[m] - phis_deg[m - 1]);
    for (std::size_t n = 1; n < n_theta; n++)
        widest_theta = std::max(widest_theta, thetas_deg[n] - thetas_deg[n - 1]);
    const bool full_circle = phis_deg[0] + 360.0 - phis_deg[n_phi - 1] <= widest_phi + 1e-9;
    // A hemisphere whose open edges are the horizontal plane (phi = -90 and +90):
    // Database reads the other half as its mirror image across that plane, so the
    // Laplacian is the even extension's, with the edge's ghost neighbor the one
    // inside it (the mirror of phi = 90 + h is phi = 90 - h).
    const bool mirrored = !full_circle && std::fabs(phis_deg.front() + 90.0) < 1e-9
            && std::fabs(phis_deg.back() - 90.0) < 1e-9;
    const bool north_pole = std::fabs(thetas_deg.front()) < 1e-9;
    const bool south_pole = std::fabs(thetas_deg.back() - 180.0) < 1e-9;

    std::vector<Unit> dir(n_cells);
    std::vector<double> alpha(n_cells);
    std::vector<cplx> unwound(n_cells);
    for (std::size_t m = 0; m < n_phi; m++)
        for (std::size_t n = 0; n < n_theta; n++)
        {
            const std::size_t cell = m * n_theta + n;
            dir[cell] = unit_of(grid_direction(phis_deg[m], thetas_deg[n]));
            alpha[cell] = azimuth(dir[cell]);
            unwound[cell] = order != 0 ? f[cell] * std::polar(1.0, -order * alpha[cell]) : f[cell];
        }
    const double fit_radius = 2.05 * std::max(widest_theta, widest_phi) * deg;
    const double cos_radius = std::cos(fit_radius);
    const bool may_unwind = std::abs(order) >= 2;

    // A node's grid neighbors (theta -1 and +1 in its phi row, phi -1 and +1 in its
    // theta row, across the seam of a full circle), npos where there is none.
    constexpr std::size_t none = static_cast<std::size_t>(-1);
    auto neighbors = [&](std::size_t m, std::size_t n) {
        std::array<std::size_t, 4> q {none, none, none, none};
        if (n > 0)
            q[0] = m * n_theta + n - 1;
        if (n + 1 < n_theta)
            q[1] = m * n_theta + n + 1;
        if (m > 0 || full_circle)
            q[2] = (m == 0 ? n_phi - 1 : m - 1) * n_theta + n;
        else if (mirrored)
            q[2] = n_theta + n;
        if (m + 1 < n_phi || full_circle)
            q[3] = (m + 1 == n_phi ? 0 : m + 1) * n_theta + n;
        else if (mirrored)
            q[3] = (n_phi - 2) * n_theta + n;
        return q;
    };
    // How much a representation varies from a node to its neighbors, relative to the
    // largest of the five values: the larger, the worse three-point differences
    // resolve it there.
    auto roughness = [&](const std::vector<cplx>& v, std::size_t cell, const std::array<std::size_t, 4>& q) {
        double peak = std::abs(v[cell]), widest = 0.0;
        for (const std::size_t j : q)
            if (j != none)
            {
                peak = std::max(peak, std::abs(v[j]));
                widest = std::max(widest, std::abs(v[j] - v[cell]));
            }
        return peak > 0.0 ? widest / peak : 0.0;
    };
    // Per node, whether to differentiate the unwound B (where it varies the less of
    // the two there; a tie unwinds) or f itself.  Never for order +-1, whose B = sin(beta) has a
    // cone point on the shaft, nor on the shaft itself, where alpha is undefined.
    std::vector<char> use_unwound(n_cells, 0);
    if (may_unwind)
        for (std::size_t m = 0; m < n_phi; m++)
            for (std::size_t n = 0; n < n_theta; n++)
            {
                const std::size_t cell = m * n_theta + n;
                const Unit& p = dir[cell];
                if (p.x * p.x + p.y * p.y < 1e-12)
                    continue;
                const auto q = neighbors(m, n);
                use_unwound[cell] = roughness(unwound, cell, q) <= roughness(f, cell, q);
            }

    // The quadratic fit at node `cell` over the nodes within fit_radius, of f or of
    // B, the winding restored by the identity in the header with grad(alpha) taken
    // by differences of the analytic azimuth.
    auto local = [&](std::size_t cell) {
        const Unit& p = dir[cell];
        const auto e = tangent_basis(p);
        const double sin2 = p.x * p.x + p.y * p.y;
        const bool unwind = use_unwound[cell] != 0;
        std::vector<std::array<double, 2>> uv;
        std::vector<cplx> y;
        for (std::size_t q = 0; q < n_cells; q++)
        {
            const cplx value = unwind ? unwound[q] : f[q];
            if (dot(p, dir[q]) >= cos_radius)
            {
                uv.push_back(tangent_coordinates(p, e, dir[q]));
                y.push_back(value);
            }
            // Its mirror image across the horizontal plane, off the edge itself.
            const Unit image {dir[q].x, dir[q].y, -dir[q].z};
            if (mirrored && std::fabs(dir[q].z) > 1e-12 && dot(p, image) >= cos_radius)
            {
                uv.push_back(tangent_coordinates(p, e, image));
                y.push_back(value);
            }
        }
        const LocalFit fit = quadratic_fit(uv, y);
        if (!unwind)
            return fit.laplacian;
        const double h = 1e-6;
        std::array<double, 2> grad {};
        for (int i = 0; i < 2; i++)
        {
            const Unit& t = e[i];
            const Unit plus {p.x + h * t.x, p.y + h * t.y, p.z + h * t.z};
            const Unit minus {p.x - h * t.x, p.y - h * t.y, p.z - h * t.z};
            grad[i] = std::remainder(azimuth(plus) - azimuth(minus), 2.0 * pi) / (2.0 * h);
        }
        const double n = order;
        const cplx dot_grad = grad[0] * fit.du + grad[1] * fit.dv;
        const cplx inner = fit.laplacian + 2.0 * cplx(0.0, n) * dot_grad - n * n / sin2 * fit.value;
        return inner * std::polar(1.0, order * alpha[cell]);
    };

    for (std::size_t m = 0; m < n_phi; m++)
    {
        for (std::size_t n = 0; n < n_theta; n++)
        {
            const std::size_t cell = m * n_theta + n;
            const Unit& p = dir[cell];
            const double sin2 = p.x * p.x + p.y * p.y;
            const bool pole = (n == 0 && north_pole) || (n + 1 == n_theta && south_pole);
            const bool phi_inside = full_circle || mirrored || (m > 0 && m + 1 < n_phi);
            const bool theta_inside = n > 0 && n + 1 < n_theta;
            if (pole || !phi_inside || !theta_inside)
            {
                out[cell] = local(cell);
                continue;
            }
            // Interior: theta - 1 and theta + 1 in this phi row (a pole row's node is
            // the pole itself, a valid neighbor); phi - 1 and phi + 1, across the seam
            // of a full circle.
            const bool unwind = use_unwound[cell] != 0;
            const double theta = thetas_deg[n] * deg;
            const double h1 = (thetas_deg[n] - thetas_deg[n - 1]) * deg;
            const double h2 = (thetas_deg[n + 1] - thetas_deg[n]) * deg;
            std::size_t mm = m == 0 ? n_phi - 1 : m - 1, mp = m + 1 == n_phi ? 0 : m + 1;
            double g1 = (m == 0 ? phis_deg[0] + 360.0 - phis_deg[n_phi - 1] : phis_deg[m] - phis_deg[m - 1]) * deg;
            double g2 = (m + 1 == n_phi ? phis_deg[0] + 360.0 - phis_deg[n_phi - 1] : phis_deg[m + 1] - phis_deg[m]) * deg;
            if (mirrored && m == 0)
            {
                mm = 1;
                g1 = g2 = (phis_deg[1] - phis_deg[0]) * deg;
            }
            if (mirrored && m + 1 == n_phi)
            {
                mp = n_phi - 2;
                g1 = g2 = (phis_deg[m] - phis_deg[m - 1]) * deg;
            }
            const auto& b = unwind ? unwound : f;
            cplx bt, btt, bp, bpp;
            three_point(h1, h2, b[cell - 1], b[cell], b[cell + 1], bt, btt);
            three_point(g1, g2, b[mm * n_theta + n], b[cell], b[mp * n_theta + n], bp, bpp);
            const double st = std::sin(theta), ct = std::cos(theta);
            const cplx lap_b = btt + (ct / st) * bt + bpp / (st * st);
            if (!unwind)
            {
                out[cell] = lap_b;
                continue;
            }
            // grad(alpha).grad(B) = alpha_theta B_theta + alpha_phi B_phi / sin^2(theta),
            // with alpha = atan2(sin(theta) sin(phi), cos(theta)) and sin^2(beta) = sin2.
            const double phi = phis_deg[m] * deg;
            const double a_theta = std::sin(phi) / sin2;
            const double a_phi = st * ct * std::cos(phi) / sin2;
            const cplx dot_grad = a_theta * bt + a_phi * bp / (st * st);
            const double nn = order;
            const cplx inner = lap_b + 2.0 * cplx(0.0, nn) * dot_grad - nn * nn / sin2 * b[cell];
            out[cell] = inner * std::polar(1.0, order * alpha[cell]);
        }
    }
    return out;
}
