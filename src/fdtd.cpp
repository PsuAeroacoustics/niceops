#include "fdtd.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>

#include "gpu/gpu.hpp"
#include "gpu/kernels/fdtd_types.h"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

bool inside_ring(const std::vector<double>& x, const std::vector<double>& y, double px, double py)
{
    bool in = false;
    const std::size_t n = x.size();
    for (std::size_t a = 0, b = n - 1; a < n; b = a++)
        if ((y[a] > py) != (y[b] > py) && px < x[b] + (x[a] - x[b]) * (py - y[b]) / (y[a] - y[b]))
            in = !in;
    return in;
}

// One axis of the grid: n cells of h from `origin`, the PML's thickness at
// each end, and the CPML's a and b at the cells' centres and low faces.
struct Axis
{
    std::size_t n = 0, lo = 0, hi = 0;
    double origin = 0.0;
    std::vector<float> a_c, b_c, a_f, b_f;

    bool in_slab(std::size_t i) const
    {
        return i < lo || i >= n - hi;
    }
    std::size_t slab(std::size_t i) const
    {
        return i < lo ? i : i - (n - hi) + lo;
    }
    std::size_t slabs() const
    {
        return lo + hi;
    }
};

Axis make_axis(double air_lo, double air_hi, double h, std::size_t lo, std::size_t hi, bool snap, double dt,
        double c)
{
    Axis a;
    const double start = snap ? std::floor(air_lo / h + 1e-9) * h : air_lo;
    const auto air = static_cast<std::size_t>(std::ceil((air_hi - start) / h - 1e-9));
    a.lo = lo;
    a.hi = hi;
    a.n = lo + std::max<std::size_t>(air, 1) + hi;
    a.origin = start - static_cast<double>(lo) * h;
    // Graded as d^3; a reflection of 1e-6 in theory, and the CFS shift for
    // the slowest waves, a few hundred hertz at most.
    const double m = 3.0;
    auto coefficients = [&](double u, std::size_t thickness, float& av, float& bv) {
        av = 0.0f;
        bv = 1.0f;
        if (thickness == 0 || u <= 0.0)
            return;
        const double d = std::min(u / static_cast<double>(thickness), 1.0);
        const double sigma_max = (m + 1.0) * c * std::log(1e6) / (2.0 * static_cast<double>(thickness) * h);
        const double sigma = sigma_max * std::pow(d, m);
        const double alpha = 2.0 * pi * 100.0 * (1.0 - d);
        const double b = std::exp(-(sigma + alpha) * dt);
        av = static_cast<float>(sigma / (sigma + alpha) * (b - 1.0));
        bv = static_cast<float>(b);
    };
    a.a_c.resize(a.n);
    a.b_c.resize(a.n);
    a.a_f.resize(a.n);
    a.b_f.resize(a.n);
    const double top = static_cast<double>(a.n - hi);
    for (std::size_t i = 0; i < a.n; i++)
    {
        const double centre = static_cast<double>(i) + 0.5, face = static_cast<double>(i);
        float av, bv;
        if (centre < static_cast<double>(lo))
            coefficients(static_cast<double>(lo) - centre, lo, av, bv);
        else
            coefficients(centre - top, hi, av, bv);
        a.a_c[i] = av;
        a.b_c[i] = bv;
        if (face < static_cast<double>(lo))
            coefficients(static_cast<double>(lo) - face, lo, av, bv);
        else
            coefficients(face - top, hi, av, bv);
        a.a_f[i] = av;
        a.b_f[i] = bv;
    }
    return a;
}

// The cells around a point and their weights: tricubic over 4 by 4 by 4 where
// those cells are all air, else trilinear over the air among 2 by 2 by 2.  A
// node beyond a rigid low side -- the floor, the plane of symmetry -- is its
// mirror cell in that side, as the field is even there.
struct Stencil
{
    static constexpr std::size_t max_cells = 64;
    std::uint32_t cells[max_cells] {};
    float weights[max_cells] {};
    std::size_t count = 0;
};

// Lagrange's cubic through the nodes -1, 0, 1, 2 at u in [0, 1).
inline void cubic_weights(double u, double w[4])
{
    w[0] = -u * (u - 1.0) * (u - 2.0) / 6.0;
    w[1] = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    w[2] = -(u + 1.0) * u * (u - 2.0) / 2.0;
    w[3] = (u + 1.0) * u * (u - 1.0) / 6.0;
}

struct Grid
{
    Axis x, y, z;
    double h = 0.0;
    std::vector<std::uint8_t> flags;

    std::size_t cells() const
    {
        return x.n * y.n * z.n;
    }
    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const
    {
        return (k * y.n + j) * x.n + i;
    }

    // False where the point is in a solid or off the grid.  `cubic` asks for
    // the tricubic stencil where it fits.
    bool stencil(const Vec3& point, double scale, bool cubic, Stencil& out) const
    {
        const Axis* axes[3] = {&x, &y, &z};
        const double coordinate[3] = {point.x, point.y, point.z};
        // The cell the point lies in, which must be air.
        std::size_t home[3];
        for (int d = 0; d < 3; d++)
        {
            const double u = (coordinate[d] - axes[d]->origin) / h;
            if (!(u >= 0.0) || u >= static_cast<double>(axes[d]->n))
                return false;
            home[d] = std::min(static_cast<std::size_t>(u), axes[d]->n - 1);
        }
        if ((flags[index(home[0], home[1], home[2])] & FDTD_AIR) == 0)
            return false;
        // Per axis, the nodes (mirrored at a rigid low side) and their weights.
        auto build = [&](int span, std::size_t node[3][4], double weight[3][4]) {
            for (int d = 0; d < 3; d++)
            {
                const Axis& a = *axes[d];
                const double u = (coordinate[d] - a.origin) / h - 0.5;
                const double first = std::floor(u), f = u - first;
                double w[4];
                if (span == 4)
                    cubic_weights(f, w);
                else
                {
                    w[0] = 1.0 - f;
                    w[1] = f;
                }
                for (int m = 0; m < span; m++)
                {
                    long long index = static_cast<long long>(first) + m - (span == 4 ? 1 : 0);
                    if (index < 0)
                    {
                        if (a.lo > 0 || index < -2)
                            return false;
                        index = -1 - index;
                    }
                    if (index >= static_cast<long long>(a.n))
                        return false;
                    node[d][m] = static_cast<std::size_t>(index);
                    weight[d][m] = w[m];
                }
            }
            return true;
        };
        std::size_t node[3][4];
        double weight[3][4];
        for (int span : {4, 2})
        {
            if (span == 4 && !cubic)
                continue;
            if (!build(span, node, weight))
                continue;
            // Each node's weight, gathered by cell (mirrored nodes fold onto
            // their cells); a solid among the cubic's sends it to trilinear.
            out.count = 0;
            double total = 0.0;
            bool solid = false;
            for (int a = 0; a < span && !solid; a++)
                for (int b = 0; b < span && !solid; b++)
                    for (int e = 0; e < span; e++)
                    {
                        const std::size_t c = index(node[0][a], node[1][b], node[2][e]);
                        const double w = weight[0][a] * weight[1][b] * weight[2][e];
                        if ((flags[c] & FDTD_AIR) == 0)
                        {
                            if (span == 4)
                            {
                                solid = true;
                                break;
                            }
                            continue;
                        }
                        total += w;
                        std::size_t t = 0;
                        while (t < out.count && out.cells[t] != c)
                            t++;
                        if (t == out.count)
                        {
                            out.cells[t] = static_cast<std::uint32_t>(c);
                            out.weights[t] = 0.0f;
                            out.count++;
                        }
                        out.weights[t] += static_cast<float>(w);
                    }
            if (solid || !(total > 0.0))
                continue;
            // Trilinear among air alone is renormalised; the cubic's weights
            // already sum to one.
            const double norm = span == 2 ? scale / total : scale;
            for (std::size_t t = 0; t < out.count; t++)
                out.weights[t] = static_cast<float>(out.weights[t] * norm);
            return true;
        }
        return false;
    }
};

Grid make_grid(const FdtdSetup& s, const std::vector<FdtdSolid>& solids, double dt)
{
    Grid g;
    g.h = s.cell;
    const auto pml = static_cast<std::size_t>(std::max(0, s.pml));
    g.x = make_axis(s.lo.x, s.hi.x, s.cell, pml, pml, true, dt, s.sound_speed);
    g.y = make_axis(s.lo.y, s.hi.y, s.cell, s.mirror ? 0 : pml, pml, !s.mirror, dt, s.sound_speed);
    g.z = make_axis(s.lo.z, s.hi.z, s.cell, s.floor ? 0 : pml, pml, !s.floor, dt, s.sound_speed);
    const std::size_t nx = g.x.n, ny = g.y.n, nz = g.z.n;
    // Solid per cell, column by column.
    std::vector<std::uint8_t> solid(g.cells(), 0);
    parallel_chunks(ny, std::max(1u, s.threads), [&](std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; j++)
        {
            const double yc = g.y.origin + (static_cast<double>(j) + 0.5) * g.h;
            for (std::size_t i = 0; i < nx; i++)
            {
                const double xc = g.x.origin + (static_cast<double>(i) + 0.5) * g.h;
                for (const FdtdSolid& b : solids)
                {
                    if (!inside_ring(b.x, b.y, xc, yc))
                        continue;
                    bool hole = false;
                    for (const Ring& r : b.holes)
                        hole = hole || inside_ring(r.x, r.y, xc, yc);
                    if (hole)
                        continue;
                    for (std::size_t k = 0; k < nz; k++)
                    {
                        const double zc = g.z.origin + (static_cast<double>(k) + 0.5) * g.h;
                        if (zc > b.bottom && zc < b.top)
                            solid[g.index(i, j, k)] = 1;
                    }
                }
            }
        }
    });
    g.flags.assign(g.cells(), 0);
    parallel_chunks(nz, std::max(1u, s.threads), [&](std::size_t begin, std::size_t end) {
        for (std::size_t k = begin; k < end; k++)
            for (std::size_t j = 0; j < ny; j++)
                for (std::size_t i = 0; i < nx; i++)
                {
                    const std::size_t c = g.index(i, j, k);
                    std::uint8_t f = solid[c] ? 0 : FDTD_AIR;
                    // Each low face: open between two air cells, a solid's
                    // between air and solid.
                    const bool low[3] = {i > 0 && !solid[c - 1], j > 0 && !solid[c - nx], k > 0 && !solid[c - nx * ny]};
                    const bool has[3] = {i > 0, j > 0, k > 0};
                    const std::uint8_t open[3] = {FDTD_OPEN_X, FDTD_OPEN_Y, FDTD_OPEN_Z};
                    const std::uint8_t wall[3] = {FDTD_WALL_X, FDTD_WALL_Y, FDTD_WALL_Z};
                    for (int d = 0; d < 3; d++)
                    {
                        if (!has[d])
                            continue;
                        if (!solid[c] && low[d])
                            f |= open[d];
                        else if (!solid[c] != low[d])
                            f |= wall[d];
                    }
                    g.flags[c] = f;
                }
    });
    // The clean cells, from the flags.
    const std::uint8_t open = FDTD_AIR | FDTD_OPEN_X | FDTD_OPEN_Y | FDTD_OPEN_Z;
    parallel_chunks(nz, std::max(1u, s.threads), [&](std::size_t begin, std::size_t end) {
        for (std::size_t k = begin; k < end; k++)
        {
            if (g.z.in_slab(k) || k < 2 || k + 2 >= nz)
                continue;
            for (std::size_t j = 2; j + 2 < ny; j++)
            {
                if (g.y.in_slab(j))
                    continue;
                for (std::size_t i = 2; i + 2 < nx; i++)
                {
                    if (g.x.in_slab(i))
                        continue;
                    const std::size_t c = g.index(i, j, k);
                    bool clean = true;
                    for (std::size_t stride : {std::size_t {1}, nx, nx * ny})
                        for (std::size_t m = 0; m < 5 && clean; m++)
                            clean = (g.flags[c + m * stride - 2 * stride] & open) == open;
                    if (clean)
                        g.flags[c] |= FDTD_CLEAN;
                }
            }
        }
    });
    return g;
}

struct Fields
{
    std::vector<float> p, vx, vy, vz, psi_x, psi_y, psi_z, phi_x, phi_y, phi_z;
};

// The CPU's march, the kernels' arithmetic in the same order.
// The impedance of the solids' faces where they absorb: v' = keep v -+ push p,
// from the half cell between the air's centre and the face.
struct Walls
{
    bool on = false;
    float keep = 1.0f, push = 0.0f;
};

// A face's pressure difference along one axis, p at the cell less p below it:
// second order, or fourth (9/8 and -1/24) with each far cell beyond a closed
// face replaced by its mirror in that face, which keeps the two differences
// each other's negative transpose and the march energy-conserving.
template <int Order>
inline float pressure_difference(const std::vector<float>& p, const std::vector<std::uint8_t>& flags, std::size_t c,
        std::size_t stride, std::size_t i, std::size_t n, std::uint8_t open)
{
    const float near = p[c] - p[c - stride];
    if constexpr (Order == 2)
        return near;
    else
    {
        const float high = i + 1 < n && (flags[c + stride] & open) ? p[c + stride] : p[c];
        const float low = (flags[c - stride] & open) ? p[c - 2 * stride] : p[c - stride];
        return 1.125f * near - (high - low) * (1.0f / 24.0f);
    }
}

// A cell's velocity difference along one axis, the high face's less the low
// face's: second order, or fourth with each far face beyond a closed face its
// mirror in that face's velocity, 2 v_face - v.  A rigid face's is zero, so
// that is the negative of the mirror; an absorbing one's is not, and with it
// the face's velocity weighs exactly 1 over the cells' divergences, as at
// second order.
template <int Order>
inline float velocity_difference(const std::vector<float>& v, const std::vector<std::uint8_t>& flags, std::size_t c,
        std::size_t stride, std::size_t i, std::size_t n, std::uint8_t open)
{
    const float v0 = v[c], v1 = i + 1 < n ? v[c + stride] : 0.0f;
    if constexpr (Order == 2)
        return v1 - v0;
    else
    {
        const bool open0 = (flags[c] & open) != 0, open1 = i + 1 < n && (flags[c + stride] & open);
        const float v2 = open1 ? (i + 2 < n ? v[c + 2 * stride] : 0.0f) : 2.0f * v1 - v0;
        const float vm = open0 ? v[c - stride] : 2.0f * v0 - v1;
        return 1.125f * (v1 - v0) - (v2 - vm) * (1.0f / 24.0f);
    }
}

// The pressure an absorbing face sees, from its air side: the cell's, or at
// fourth order 25/24 of it less 1/24 of the next cell's where that face is
// open -- the weights the face's velocity carries in the two cells'
// divergences (velocity_difference), so the wall's update is their transpose
// and the march dissipates for any impedance.  `above` is the air's side.
template <int Order>
inline float wall_pressure(const std::vector<float>& p, const std::vector<std::uint8_t>& flags, std::size_t c,
        std::size_t stride, std::size_t i, std::size_t n, std::uint8_t open, bool above)
{
    const std::size_t air = above ? c : c - stride;
    if constexpr (Order == 4)
    {
        if (above && i + 1 < n && (flags[c + stride] & open))
            return (25.0f / 24.0f) * p[c] - (1.0f / 24.0f) * p[c + stride];
        if (!above && (flags[c - stride] & open))
            return (25.0f / 24.0f) * p[air] - (1.0f / 24.0f) * p[c - 2 * stride];
    }
    return p[air];
}

// How much wider than the fastest wave's reach each step's box is, cells:
// the dispersive front's tail.  Four cells take a culled march to its float
// floor against a full one at either order, on 310 and 620 cells of path.
constexpr double cull_margin = 16.0;

// The bounding box of stencils' cells, inclusive, [lo, hi] along each axis.
struct CellBounds
{
    std::size_t lo[3] {SIZE_MAX, SIZE_MAX, SIZE_MAX};
    std::size_t hi[3] {0, 0, 0};

    void take(const Grid& g, const Stencil& stencil)
    {
        for (std::size_t t = 0; t < stencil.count; t++)
        {
            const std::size_t c = stencil.cells[t], row = c / g.x.n;
            const std::size_t at[3] {c - row * g.x.n, row % g.y.n, row / g.y.n};
            for (int a = 0; a < 3; a++)
            {
                lo[a] = std::min(lo[a], at[a]);
                hi[a] = std::max(hi[a], at[a]);
            }
        }
    }
};

// Each step's box (FdtdSetup::cull): the cells the source can have reached
// by the step's end, at `speed` cells a step, and from which a receiver can
// still be reached before the last step's, `margin` cells wider each way.
// Without `cull`, or receivers, the box is the grid's (or only the source's
// reach bounds it).
std::vector<FdtdBox> cone_boxes(const Grid& g, bool cull, const Stencil& source, const std::vector<Stencil>& receivers,
        std::size_t steps, double speed, double margin)
{
    const std::size_t n[3] {g.x.n, g.y.n, g.z.n};
    CellBounds from, to;
    from.take(g, source);
    for (const Stencil& r : receivers)
        to.take(g, r);
    std::vector<FdtdBox> boxes(steps);
    for (std::size_t step = 0; step < steps; step++)
    {
        const double ahead = margin + speed * static_cast<double>(step + 1);
        const double behind = margin + speed * static_cast<double>(steps - step);
        std::size_t lo[3], hi[3];
        for (int a = 0; a < 3; a++)
        {
            double low = 0.0, high = static_cast<double>(n[a]);
            if (cull)
            {
                low = std::max(low, static_cast<double>(from.lo[a]) - ahead);
                high = std::min(high, static_cast<double>(from.hi[a]) + 1.0 + ahead);
                if (!receivers.empty())
                {
                    low = std::max(low, static_cast<double>(to.lo[a]) - behind);
                    high = std::min(high, static_cast<double>(to.hi[a]) + 1.0 + behind);
                }
            }
            lo[a] = static_cast<std::size_t>(std::floor(low));
            hi[a] = std::max(lo[a], static_cast<std::size_t>(std::max(0.0, std::ceil(high))));
        }
        FdtdBox& b = boxes[step];
        b.i = static_cast<std::uint32_t>(lo[0]);
        b.j = static_cast<std::uint32_t>(lo[1]);
        b.k = static_cast<std::uint32_t>(lo[2]);
        b.ni = static_cast<std::uint32_t>(hi[0] - lo[0]);
        b.nj = static_cast<std::uint32_t>(hi[1] - lo[1]);
        b.nk = static_cast<std::uint32_t>(hi[2] - lo[2]);
        if (b.ni == 0 || b.nj == 0 || b.nk == 0)
            b.ni = b.nj = b.nk = 0;
    }
    return boxes;
}

// The clean cells [begin, end) of a row: the plain differences, with nothing
// to look up, the same arithmetic as the general path's.
template <int Order>
void clean_velocity(const float* __restrict p, float* __restrict vx, float* __restrict vy, float* __restrict vz,
        std::size_t begin, std::size_t end, std::size_t nx, std::size_t plane, float push)
{
    for (std::size_t c = begin; c < end; c++)
    {
        if constexpr (Order == 2)
        {
            vx[c] -= push * (p[c] - p[c - 1]);
            vy[c] -= push * (p[c] - p[c - nx]);
            vz[c] -= push * (p[c] - p[c - plane]);
        }
        else
        {
            vx[c] -= push * (1.125f * (p[c] - p[c - 1]) - (p[c + 1] - p[c - 2]) * (1.0f / 24.0f));
            vy[c] -= push * (1.125f * (p[c] - p[c - nx]) - (p[c + nx] - p[c - 2 * nx]) * (1.0f / 24.0f));
            vz[c] -= push * (1.125f * (p[c] - p[c - plane]) - (p[c + plane] - p[c - 2 * plane]) * (1.0f / 24.0f));
        }
    }
}

template <int Order>
void clean_pressure(float* __restrict p, const float* __restrict vx, const float* __restrict vy,
        const float* __restrict vz, std::size_t begin, std::size_t end, std::size_t nx, std::size_t plane, float push)
{
    for (std::size_t c = begin; c < end; c++)
    {
        float dx, dy, dz;
        if constexpr (Order == 2)
        {
            dx = vx[c + 1] - vx[c];
            dy = vy[c + nx] - vy[c];
            dz = vz[c + plane] - vz[c];
        }
        else
        {
            dx = 1.125f * (vx[c + 1] - vx[c]) - (vx[c + 2] - vx[c - 1]) * (1.0f / 24.0f);
            dy = 1.125f * (vy[c + nx] - vy[c]) - (vy[c + 2 * nx] - vy[c - nx]) * (1.0f / 24.0f);
            dz = 1.125f * (vz[c + plane] - vz[c]) - (vz[c + 2 * plane] - vz[c - plane]) * (1.0f / 24.0f);
        }
        p[c] -= push * (dx + dy + dz);
    }
}

template <int Order>
void march_cpu(const Grid& g, const FdtdSetup& s, float push_v, float push_p, const Walls& walls,
        const Stencil& source, const std::vector<float>& pulse, const std::vector<Stencil>& receivers,
        const std::vector<FdtdBox>& boxes, std::vector<float>& recording)
{
    const std::size_t steps = boxes.size();
    const std::size_t nx = g.x.n, ny = g.y.n, nz = g.z.n, plane = nx * ny, n = g.cells();
    Fields f;
    f.p.assign(n, 0.0f);
    f.vx.assign(n, 0.0f);
    f.vy.assign(n, 0.0f);
    f.vz.assign(n, 0.0f);
    f.psi_x.assign(g.x.slabs() * ny * nz, 0.0f);
    f.phi_x.assign(f.psi_x.size(), 0.0f);
    f.psi_y.assign(nx * g.y.slabs() * nz, 0.0f);
    f.phi_y.assign(f.psi_y.size(), 0.0f);
    f.psi_z.assign(nx * ny * g.z.slabs(), 0.0f);
    f.phi_z.assign(f.psi_z.size(), 0.0f);
    const unsigned threads = std::max(1u, s.threads);
    const std::size_t count = receivers.size();
    const std::vector<std::uint8_t>& fl_ = g.flags;
    for (std::size_t step = 0; step < steps; step++)
    {
        const FdtdBox& b = boxes[step];
        const std::size_t i0 = b.i, i1 = b.i + b.ni, j0 = b.j, j1 = b.j + b.nj;
        parallel_chunks(b.nk, threads, [&](std::size_t begin, std::size_t end) {
            for (std::size_t k = b.k + begin; k < b.k + end; k++)
                for (std::size_t j = j0; j < j1; j++)
                {
                    const std::size_t row = k * ny + j;
                    for (std::size_t i = i0; i < i1; i++)
                    {
                        const std::size_t c = row * nx + i;
                        const std::uint8_t fl = fl_[c];
                        if (fl & FDTD_CLEAN)
                        {
                            // The run of clean cells from here.
                            std::size_t e = i + 1;
                            while (e < i1 && (fl_[row * nx + e] & FDTD_CLEAN))
                                e++;
                            clean_velocity<Order>(f.p.data(), f.vx.data(), f.vy.data(), f.vz.data(), c, row * nx + e,
                                    nx, plane, push_v);
                            i = e - 1;
                            continue;
                        }
                        if ((fl & FDTD_FACES) == 0)
                            continue;
                        const bool air = (fl & FDTD_AIR) != 0;
                        if (fl & FDTD_OPEN_X)
                        {
                            float d = pressure_difference<Order>(f.p, fl_, c, 1, i, nx, FDTD_OPEN_X);
                            if (g.x.in_slab(i))
                            {
                                const std::size_t q = row * g.x.slabs() + g.x.slab(i);
                                const float psi = g.x.b_f[i] * f.psi_x[q] + g.x.a_f[i] * d;
                                f.psi_x[q] = psi;
                                d += psi;
                            }
                            f.vx[c] -= push_v * d;
                        }
                        else if (walls.on && (fl & FDTD_WALL_X))
                        {
                            const float pw = wall_pressure<Order>(f.p, fl_, c, 1, i, nx, FDTD_OPEN_X, air);
                            f.vx[c] = air ? walls.keep * f.vx[c] - walls.push * pw : walls.keep * f.vx[c] + walls.push * pw;
                        }
                        if (fl & FDTD_OPEN_Y)
                        {
                            float d = pressure_difference<Order>(f.p, fl_, c, nx, j, ny, FDTD_OPEN_Y);
                            if (g.y.in_slab(j))
                            {
                                const std::size_t q = (k * g.y.slabs() + g.y.slab(j)) * nx + i;
                                const float psi = g.y.b_f[j] * f.psi_y[q] + g.y.a_f[j] * d;
                                f.psi_y[q] = psi;
                                d += psi;
                            }
                            f.vy[c] -= push_v * d;
                        }
                        else if (walls.on && (fl & FDTD_WALL_Y))
                        {
                            const float pw = wall_pressure<Order>(f.p, fl_, c, nx, j, ny, FDTD_OPEN_Y, air);
                            f.vy[c] = air ? walls.keep * f.vy[c] - walls.push * pw : walls.keep * f.vy[c] + walls.push * pw;
                        }
                        if (fl & FDTD_OPEN_Z)
                        {
                            float d = pressure_difference<Order>(f.p, fl_, c, plane, k, nz, FDTD_OPEN_Z);
                            if (g.z.in_slab(k))
                            {
                                const std::size_t q = g.z.slab(k) * plane + j * nx + i;
                                const float psi = g.z.b_f[k] * f.psi_z[q] + g.z.a_f[k] * d;
                                f.psi_z[q] = psi;
                                d += psi;
                            }
                            f.vz[c] -= push_v * d;
                        }
                        else if (walls.on && (fl & FDTD_WALL_Z))
                        {
                            const float pw = wall_pressure<Order>(f.p, fl_, c, plane, k, nz, FDTD_OPEN_Z, air);
                            f.vz[c] = air ? walls.keep * f.vz[c] - walls.push * pw : walls.keep * f.vz[c] + walls.push * pw;
                        }
                    }
                }
        });
        parallel_chunks(b.nk, threads, [&](std::size_t begin, std::size_t end) {
            for (std::size_t k = b.k + begin; k < b.k + end; k++)
                for (std::size_t j = j0; j < j1; j++)
                {
                    const std::size_t row = k * ny + j;
                    for (std::size_t i = i0; i < i1; i++)
                    {
                        const std::size_t c = row * nx + i;
                        if (fl_[c] & FDTD_CLEAN)
                        {
                            std::size_t e = i + 1;
                            while (e < i1 && (fl_[row * nx + e] & FDTD_CLEAN))
                                e++;
                            clean_pressure<Order>(f.p.data(), f.vx.data(), f.vy.data(), f.vz.data(), c, row * nx + e,
                                    nx, plane, push_p);
                            i = e - 1;
                            continue;
                        }
                        if ((fl_[c] & FDTD_AIR) == 0)
                            continue;
                        float dx = velocity_difference<Order>(f.vx, fl_, c, 1, i, nx, FDTD_OPEN_X);
                        if (g.x.in_slab(i))
                        {
                            const std::size_t q = row * g.x.slabs() + g.x.slab(i);
                            const float phi = g.x.b_c[i] * f.phi_x[q] + g.x.a_c[i] * dx;
                            f.phi_x[q] = phi;
                            dx += phi;
                        }
                        float dy = velocity_difference<Order>(f.vy, fl_, c, nx, j, ny, FDTD_OPEN_Y);
                        if (g.y.in_slab(j))
                        {
                            const std::size_t q = (k * g.y.slabs() + g.y.slab(j)) * nx + i;
                            const float phi = g.y.b_c[j] * f.phi_y[q] + g.y.a_c[j] * dy;
                            f.phi_y[q] = phi;
                            dy += phi;
                        }
                        float dz = velocity_difference<Order>(f.vz, fl_, c, plane, k, nz, FDTD_OPEN_Z);
                        if (g.z.in_slab(k))
                        {
                            const std::size_t q = g.z.slab(k) * plane + j * nx + i;
                            const float phi = g.z.b_c[k] * f.phi_z[q] + g.z.a_c[k] * dz;
                            f.phi_z[q] = phi;
                            dz += phi;
                        }
                        f.p[c] -= push_p * (dx + dy + dz);
                    }
                }
        });
        for (std::size_t t = 0; t < source.count; t++)
            f.p[source.cells[t]] += source.weights[t] * pulse[step];
        for (std::size_t r = 0; r < count; r++)
        {
            float v = 0.0f;
            for (std::size_t t = 0; t < receivers[r].count; t++)
                v += receivers[r].weights[t] * f.p[receivers[r].cells[t]];
            recording[step * count + r] = v;
        }
    }
}

// The same on the GPU, either order: four dispatches a step, 256 steps a
// batch, two batches in flight.
void march_gpu(const Grid& g, int order, float push_v, float push_p, const Walls& walls, const Stencil& source,
        const std::vector<float>& pulse, const std::vector<Stencil>& receivers, const std::vector<FdtdBox>& boxes,
        std::vector<float>& recording)
{
    const std::size_t steps = boxes.size();
    gpu::Device& device = gpu::Device::require("--gpu");
    const std::size_t nx = g.x.n, ny = g.y.n, nz = g.z.n, n = g.cells();
    if (n > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("FDTD: " + std::to_string(n) + " cells is more than the GPU's kernels index");
    FdtdGridGpu grid {};
    grid.nx = static_cast<std::uint32_t>(nx);
    grid.ny = static_cast<std::uint32_t>(ny);
    grid.nz = static_cast<std::uint32_t>(nz);
    grid.x0 = static_cast<std::uint32_t>(g.x.lo);
    grid.x1 = static_cast<std::uint32_t>(g.x.hi);
    grid.y0 = static_cast<std::uint32_t>(g.y.lo);
    grid.y1 = static_cast<std::uint32_t>(g.y.hi);
    grid.z0 = static_cast<std::uint32_t>(g.z.lo);
    grid.z1 = static_cast<std::uint32_t>(g.z.hi);
    grid.push_v = push_v;
    grid.push_p = push_p;
    grid.keep = walls.keep;
    grid.push_wall = walls.push;
    grid.walls = walls.on ? 1u : 0u;
    std::vector<float> coefficient;
    for (const Axis* a : {&g.x, &g.y, &g.z})
        for (const std::vector<float>* v : {&a->a_c, &a->b_c, &a->a_f, &a->b_f})
            coefficient.insert(coefficient.end(), v->begin(), v->end());

    const std::vector<std::string> program {"fdtd_types.h", "fdtd.kern"};
    const gpu::Kernel velocity = device.kernel(program, order == 4 ? "fdtd_velocity4" : "fdtd_velocity");
    const gpu::Kernel pressure = device.kernel(program, order == 4 ? "fdtd_pressure4" : "fdtd_pressure");
    const gpu::Kernel inject = device.kernel(program, "fdtd_inject");
    const gpu::Kernel record = device.kernel(program, "fdtd_record");
    auto upload = [&](const void* data, std::size_t bytes) {
        gpu::Buffer b = device.buffer(bytes);
        if (bytes > 0)
            std::memcpy(b.data(), data, bytes);
        return b;
    };
    const gpu::Buffer grid_buffer = upload(&grid, sizeof(grid));
    const gpu::Buffer flags = upload(g.flags.data(), n);
    const gpu::Buffer coefficients = upload(coefficient.data(), coefficient.size() * sizeof(float));
    const gpu::Buffer p = device.buffer(n * sizeof(float)), vx = device.buffer(n * sizeof(float)),
                      vy = device.buffer(n * sizeof(float)), vz = device.buffer(n * sizeof(float));
    const std::size_t sx = g.x.slabs() * ny * nz, sy = nx * g.y.slabs() * nz, sz = nx * ny * g.z.slabs();
    const gpu::Buffer psi_x = device.buffer(sx * sizeof(float)), phi_x = device.buffer(sx * sizeof(float));
    const gpu::Buffer psi_y = device.buffer(sy * sizeof(float)), phi_y = device.buffer(sy * sizeof(float));
    const gpu::Buffer psi_z = device.buffer(sz * sizeof(float)), phi_z = device.buffer(sz * sizeof(float));
    std::vector<std::uint32_t> source_cells(source.cells, source.cells + source.count);
    std::vector<float> source_weights(source.weights, source.weights + source.count);
    const gpu::Buffer source_cell_buffer = upload(source_cells.data(), source_cells.size() * sizeof(std::uint32_t));
    const gpu::Buffer source_weight_buffer = upload(source_weights.data(), source_weights.size() * sizeof(float));
    const gpu::Buffer pulse_buffer = upload(pulse.data(), pulse.size() * sizeof(float));
    // Each receiver's cells, padded to the longest stencil with zero weights.
    std::size_t per = 1;
    for (const Stencil& r : receivers)
        per = std::max(per, r.count);
    std::vector<std::uint32_t> receiver_cells(receivers.size() * per, 0);
    std::vector<float> receiver_weights(receivers.size() * per, 0.0f);
    for (std::size_t r = 0; r < receivers.size(); r++)
        for (std::size_t t = 0; t < receivers[r].count; t++)
        {
            receiver_cells[r * per + t] = receivers[r].cells[t];
            receiver_weights[r * per + t] = receivers[r].weights[t];
        }
    const gpu::Buffer receiver_cell_buffer
            = upload(receiver_cells.data(), receiver_cells.size() * sizeof(std::uint32_t));
    const gpu::Buffer receiver_weight_buffer = upload(receiver_weights.data(), receiver_weights.size() * sizeof(float));
    const std::size_t count = receivers.size();
    const gpu::Buffer recording_buffer = device.buffer(std::max<std::size_t>(steps * count, 1) * sizeof(float));

    const std::size_t threads = 256;
    const auto source_count = static_cast<std::uint32_t>(source_cells.size());
    const auto receiver_count = static_cast<std::uint32_t>(count);
    const std::size_t chunk = 256;
    std::deque<gpu::Batch> flight;
    for (std::size_t begin = 0; begin < steps; begin += chunk)
    {
        const std::size_t end = std::min(begin + chunk, steps);
        gpu::Batch batch = device.batch();
        for (std::size_t step = begin; step < end; step++)
        {
            const auto at = static_cast<std::uint32_t>(step);
            const FdtdBox& box = boxes[step];
            const std::size_t marched = std::size_t {box.ni} * box.nj * box.nk;
            if (marched > 0)
            {
                const std::size_t groups = (marched + threads - 1) / threads;
                batch.dispatch(velocity,
                        {gpu::Arg(grid_buffer), gpu::Arg(flags), gpu::Arg(coefficients), gpu::Arg(p), gpu::Arg(vx),
                                gpu::Arg(vy), gpu::Arg(vz), gpu::Arg(psi_x), gpu::Arg(psi_y), gpu::Arg(psi_z),
                                gpu::Arg::value(box)},
                        groups, threads);
                batch.dispatch(pressure,
                        {gpu::Arg(grid_buffer), gpu::Arg(flags), gpu::Arg(coefficients), gpu::Arg(p), gpu::Arg(vx),
                                gpu::Arg(vy), gpu::Arg(vz), gpu::Arg(phi_x), gpu::Arg(phi_y), gpu::Arg(phi_z),
                                gpu::Arg::value(box)},
                        groups, threads);
            }
            batch.dispatch(inject,
                    {gpu::Arg(p), gpu::Arg(source_cell_buffer), gpu::Arg(source_weight_buffer), gpu::Arg(pulse_buffer),
                            gpu::Arg::value(source_count), gpu::Arg::value(at)},
                    (source_cells.size() + 31) / 32, 32);
            if (count > 0)
                batch.dispatch(record,
                        {gpu::Arg(p), gpu::Arg(receiver_cell_buffer), gpu::Arg(receiver_weight_buffer),
                                gpu::Arg(recording_buffer), gpu::Arg::value(receiver_count), gpu::Arg::value(at),
                                gpu::Arg::value(static_cast<std::uint32_t>(per))},
                        (count + 63) / 64, 64);
        }
        batch.commit();
        flight.push_back(batch);
        if (flight.size() == 2)
        {
            flight.front().wait();
            flight.pop_front();
        }
    }
    while (!flight.empty())
    {
        flight.front().wait();
        flight.pop_front();
    }
    recording.assign(recording_buffer.as<float>(), recording_buffer.as<float>() + steps * count);
}

} // namespace

std::vector<FdtdSolid> fdtd_solids(const std::vector<Footprint>& footprints)
{
    std::vector<FdtdSolid> out;
    for (const Footprint& f : footprints)
    {
        FdtdSolid s;
        s.x = f.x;
        s.y = f.y;
        s.holes = f.holes;
        s.bottom = -1e9;
        s.top = f.top;
        out.push_back(s);
    }
    return out;
}

// sin(omega dt / 2) = (nu / 2) |S(k h)|, S the difference's symbol along
// each axis, so d omega / d k_i = c S_i S_i' / (|S| cos(omega dt / 2)).
// Sampled on a 128-point lattice per axis, whose spacing costs less than 1e-3
// of the maximum.
double fdtd_fastest_group_speed(int order, double nu)
{
    // Per axis, S^2 and (S S')^2, with S = 2 D the symbol scaled to k h at
    // small k.
    constexpr int n = 128;
    double square[n], flow[n];
    for (int m = 0; m < n; m++)
    {
        const double x = pi * (m + 1) / n;
        const double symbol
                = order == 4 ? (27.0 * std::sin(x / 2.0) - std::sin(1.5 * x)) / 12.0 : 2.0 * std::sin(x / 2.0);
        const double slope
                = order == 4 ? (27.0 * std::cos(x / 2.0) - 3.0 * std::cos(1.5 * x)) / 24.0 : std::cos(x / 2.0);
        square[m] = symbol * symbol;
        flow[m] = square[m] * slope * slope;
    }
    double fastest = 1.0; // k -> 0
    for (int a = 0; a < n; a++)
        for (int b = a; b < n; b++)
            for (int e = b; e < n; e++)
            {
                const double norm = square[a] + square[b] + square[e];
                fastest = std::max(fastest, (flow[a] + flow[b] + flow[e]) / (norm * (1.0 - 0.25 * nu * nu * norm)));
            }
    return std::sqrt(fastest);
}

FdtdRecording solve_fdtd(const FdtdSetup& setup, const std::vector<FdtdSolid>& solids, const Vec3& source,
        const std::vector<Vec3>& receivers)
{
    if (!(setup.cell > 0.0) || !(setup.duration > 0.0) || !(setup.courant > 0.0 && setup.courant <= 1.0))
        throw std::invalid_argument("FDTD: the cell, duration and Courant number must be positive, the last at most 1");
    if (!(setup.wall_absorption >= 0.0 && setup.wall_absorption < 1.0))
        throw std::invalid_argument("FDTD: the walls' absorption must lie within [0, 1)");
    if (setup.mirror && std::abs(source.y - setup.lo.y) > 1e-6 * setup.cell)
        throw std::invalid_argument("FDTD: a symmetry plane at y = " + std::to_string(setup.lo.y)
                + " needs the source on it, not at y = " + std::to_string(source.y));
    if (setup.order != 2 && setup.order != 4)
        throw std::invalid_argument("FDTD: the order must be 2 or 4");
    const double c = setup.sound_speed, rho = setup.density, h = setup.cell;
    // The scheme's limit: 1 / sqrt(3) over the sum of the difference's
    // coefficients, 1 or 9/8 + 1/24.
    const double dt = setup.courant * h / (c * std::sqrt(3.0) * (setup.order == 4 ? 7.0 / 6.0 : 1.0));
    const Grid g = make_grid(setup, solids, dt);

    FdtdRecording out;
    out.dt = dt;
    out.cell = h;
    out.density = rho;
    out.sound_speed = c;
    out.nx = g.x.n;
    out.ny = g.y.n;
    out.nz = g.z.n;
    out.cells = g.cells();
    out.receivers = receivers.size();
    out.steps = static_cast<std::size_t>(std::ceil(setup.duration / dt));
    out.order = setup.order;
    // Second order's phase error depends on direction, so no one warp helps it.
    out.warp = setup.warp && setup.order == 4;

    // Q(t) a Gaussian: rho Q' peaks in spectrum at omega = sqrt(2) / tau.
    const double tau = std::sqrt(2.0) / (2.0 * pi * setup.pulse_frequency);
    const double t0 = 4.5 * tau;
    out.source.resize(out.steps);
    std::vector<float> pulse(out.steps);
    for (std::size_t n = 0; n < out.steps; n++)
    {
        const double u = ((static_cast<double>(n) + 0.5) * dt - t0) / tau;
        out.source[n] = std::exp(-u * u);
        pulse[n] = static_cast<float>(out.source[n]);
    }
    // p^{n+1} += rho c^2 dt Q / h^3, by the stencil's weights.
    Stencil source_stencil;
    const double half = setup.mirror ? 0.5 : 1.0;
    const bool cubic = setup.order == 4;
    if (!g.stencil(source, half * rho * c * c * dt / (h * h * h), cubic, source_stencil))
        throw std::invalid_argument("FDTD: the source is in a solid or off the grid");
    std::vector<Stencil> receiver_stencils(receivers.size());
    for (std::size_t r = 0; r < receivers.size(); r++)
        if (!g.stencil(receivers[r], 1.0, cubic, receiver_stencils[r]))
            throw std::invalid_argument("FDTD: receiver " + std::to_string(r + 1) + " is in a solid or off the grid");

    const auto push_v = static_cast<float>(dt / (rho * h));
    const auto push_p = static_cast<float>(rho * c * c * dt / h);
    // A real impedance zeta rho c whose normal-incidence absorption is alpha:
    // rho (v' - v) / dt = -(Z (v' + v) / 2 - p) / (h / 2) at a face with air
    // below it, and the mirror of that with air above.
    Walls walls;
    if (setup.wall_absorption > 0.0)
    {
        const double root = std::sqrt(1.0 - setup.wall_absorption);
        const double z = (1.0 + root) / (1.0 - root) * rho * c;
        walls.on = true;
        walls.keep = static_cast<float>((rho / dt - z / h) / (rho / dt + z / h));
        walls.push = static_cast<float>((2.0 / h) / (rho / dt + z / h));
    }
    out.pulse_delay = t0;
    const std::vector<FdtdBox> boxes = cone_boxes(g, setup.cull, source_stencil, receiver_stencils, out.steps,
            fdtd_fastest_group_speed(setup.order, c * dt / h) * c * dt / h, cull_margin);
    double marched = 0.0;
    for (const FdtdBox& b : boxes)
        marched += static_cast<double>(b.ni) * b.nj * b.nk;
    out.marched = out.steps > 0 ? marched / (static_cast<double>(out.cells) * static_cast<double>(out.steps)) : 1.0;
    const auto start = std::chrono::steady_clock::now();
    if (setup.gpu)
        march_gpu(g, setup.order, push_v, push_p, walls, source_stencil, pulse, receiver_stencils, boxes,
                out.pressure);
    else
    {
        out.pressure.assign(out.steps * receivers.size(), 0.0f);
        if (setup.order == 4)
            march_cpu<4>(g, setup, push_v, push_p, walls, source_stencil, pulse, receiver_stencils, boxes,
                    out.pressure);
        else
            march_cpu<2>(g, setup, push_v, push_p, walls, source_stencil, pulse, receiver_stencils, boxes,
                    out.pressure);
    }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return out;
}

std::complex<double> FdtdRecording::transfer(std::size_t receiver, double frequency, std::size_t used,
        double attenuation) const
{
    const std::size_t count = used == 0 ? steps : std::min(used, steps);
    const double omega = 2.0 * pi * grid_frequency(frequency);
    const std::complex<double> turn = std::polar(1.0, omega * dt);
    // The pressure at (n + 1) dt, Q at (n + 1/2) dt; a sample's path, c times
    // its delay after the pulse's peak, loses exp(-a L).
    std::complex<double> phase_p = turn, phase_q = std::polar(1.0, 0.5 * omega * dt);
    std::complex<double> P(0.0, 0.0), Q(0.0, 0.0);
    const double loss = std::exp(-attenuation * sound_speed * dt);
    double kept = 1.0;
    for (std::size_t n = 0; n < count; n++)
    {
        const double t = (static_cast<double>(n) + 1.0) * dt;
        if (attenuation > 0.0)
            kept = t <= pulse_delay ? 1.0
                    : (n % 256 == 0 ? std::exp(-attenuation * sound_speed * (t - pulse_delay)) : kept * loss);
        P += static_cast<double>(pressure[n * receivers + receiver]) * kept * phase_p;
        Q += source[n] * phase_q;
        phase_p *= turn;
        phase_q *= turn;
        // Renormalised now and then, so the rotation does not drift.
        if (n % 256 == 255)
        {
            phase_p /= std::abs(phase_p);
            phase_q /= std::abs(phase_q);
        }
    }
    return P * (4.0 * pi) / (std::complex<double>(0.0, -omega) * density * Q);
}

double FdtdRecording::grid_frequency(double frequency) const
{
    if (!warp || !(frequency > 0.0))
        return frequency;
    // sin(omega dt / 2) = C sqrt(sum D(k_d h)^2), D the difference's symbol,
    // C = c dt / h; omega averaged over directions, by solid angle over an
    // octant on a 16 by 16 grid.
    const double kh = 2.0 * pi * frequency / sound_speed * cell, courant = sound_speed * dt / cell;
    auto symbol = [&](double x) {
        return order == 4 ? (27.0 * std::sin(0.5 * x) - std::sin(1.5 * x)) / 24.0 : std::sin(0.5 * x);
    };
    const int n = 16;
    double sum = 0.0, weights = 0.0;
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
        {
            const double theta = (a + 0.5) * 0.5 * pi / n, phi = (b + 0.5) * 0.5 * pi / n;
            const double u[3] = {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
            double s2 = 0.0;
            for (double component : u)
                s2 += std::pow(symbol(kh * component), 2);
            const double s = std::min(1.0, courant * std::sqrt(s2));
            sum += std::sin(theta) * 2.0 * std::asin(s) / dt;
            weights += std::sin(theta);
        }
    return sum / weights / (2.0 * pi);
}

double FdtdRecording::mean_energy(std::size_t receiver, double f1, double f2, std::size_t used,
        std::size_t per_band, double attenuation) const
{
    const double duration = static_cast<double>(steps) * dt;
    const auto points
            = std::max<std::size_t>(per_band, static_cast<std::size_t>(std::ceil(4.0 * (f2 - f1) * duration)));
    double sum = 0.0;
    for (std::size_t m = 0; m < points; m++)
    {
        const double f = f1 + (static_cast<double>(m) + 0.5) * (f2 - f1) / static_cast<double>(points);
        sum += std::norm(transfer(receiver, f, used, attenuation));
    }
    return sum / static_cast<double>(points);
}

double FdtdRecording::band_level(std::size_t receiver, double band, std::size_t used, const double* attenuation) const
{
    const double half = std::pow(2.0, 1.0 / 6.0), f1 = band / half, f2 = band * half;
    if (!attenuation)
        return 10.0 * std::log10(std::max(mean_energy(receiver, f1, f2, used), 1e-300));
    // Each third of the band in log frequency at its own attenuation, as the
    // ray model's sub-bands (point_paths.hpp's band_level).
    double sum = 0.0;
    for (int part = 0; part < 3; part++)
    {
        const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
        sum += mean_energy(receiver, g1, g2, used, 24, attenuation[part]) * (g2 - g1);
    }
    return 10.0 * std::log10(std::max(sum / (f2 - f1), 1e-300));
}

double FdtdRecording::truncation_db(std::size_t receiver, double band) const
{
    return band_level(receiver, band) - band_level(receiver, band, steps - steps / 10);
}

} // namespace niceops
