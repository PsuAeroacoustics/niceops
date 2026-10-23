#include "turbulence_scattering.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace niceops {

double velocity_equivalent(double gamma_t, double sound_speed)
{
    return 3.0 / 22.0 * sound_speed * sound_speed * gamma_t;
}

double scattered_level(double h, double d_s, double d_r, double frequency, double cv2)
{
    if (!(h > 0.0) || !(cv2 > 0.0) || !(frequency > 0.0))
        return -std::numeric_limits<double>::infinity();
    // A screen at either end leaves the fit's domain (d_S and d_R from 10 m);
    // a millimetre keeps the logarithm finite and the level falling.
    d_s = std::max(d_s, 1e-3);
    d_r = std::max(d_r, 1e-3);
    const double b1 = -52.8 + 10.0 * std::log10(cv2), b2 = 11.3, b3 = -17.1;
    return b1 + b2 * std::log10(h / 10.0) + b3 * std::log10(h * h / (d_s * d_r) + 0.0012)
            + 10.0 / 3.0 * std::log10(frequency / 1000.0);
}

ShadowScreen shadow_screen(const Buildings& buildings, const Vec3& source, const Vec3& receiver)
{
    ShadowScreen out;
    const Vec3 d = receiver - source;
    const double length = d.norm();
    if (!(length > 0.0))
        return out;
    const Vec3 u = d / length;
    thread_local std::vector<std::uint32_t> near;
    buildings.near_box(std::min(source.x, receiver.x), std::max(source.x, receiver.x), std::min(source.y, receiver.y),
            std::max(source.y, receiver.y), near);
    // Where the plan segment crosses the side from a to b: its parameter t
    // along the segment, both within [0, 1].
    auto crossing = [&](double ax, double ay, double bx, double by, double& t) {
        const double ex = bx - ax, ey = by - ay;
        const double det = d.x * ey - d.y * ex;
        if (std::abs(det) < 1e-12)
            return false;
        const double wx = ax - source.x, wy = ay - source.y;
        t = (wx * ey - wy * ex) / det;
        const double w = (wx * d.y - wy * d.x) / det;
        return t > 0.0 && t < 1.0 && w >= 0.0 && w <= 1.0;
    };
    auto consider = [&](double t, double top) {
        const double above = top - (source.z + t * d.z);
        if (!(above > 0.0))
            return;
        const Vec3 q = Vec3(source.x + t * d.x, source.y + t * d.y, top) - source;
        const double along = q.dot(u);
        const double h = (q - along * u).norm();
        if (h > out.h)
            out = {true, h, along, length - along};
    };
    for (const std::uint32_t b : near)
    {
        const Footprint& f = buildings.footprints()[b];
        if (!f.roof.empty())
        {
            // A pitched roof: its eaves, rakes and ridges where the line
            // crosses them, at their heights there.
            for (std::size_t e = buildings.edge_begin(b); e < buildings.edge_begin(b + 1); e++)
            {
                const Edge& g = buildings.edges()[e];
                if (g.vertical)
                    continue;
                const Vec3 end = g.start + g.length * g.direction;
                double t;
                if (crossing(g.start.x, g.start.y, end.x, end.y, t))
                    consider(t, roof_height(f, source.x + t * d.x, source.y + t * d.y));
            }
            continue;
        }
        auto ring = [&](const std::vector<double>& x, const std::vector<double>& y) {
            for (std::size_t i = 0; i < x.size(); i++)
            {
                const std::size_t j = (i + 1) % x.size();
                double t;
                if (crossing(x[i], y[i], x[j], y[j], t))
                    consider(t, f.top);
            }
        };
        ring(f.x, f.y);
        for (const Ring& hole : f.holes)
            ring(hole.x, hole.y);
    }
    return out;
}

} // namespace niceops
