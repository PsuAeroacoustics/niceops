#include "contour.hpp"

#include <fstream>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <unordered_map>

#include <nlohmann/json.hpp>

#include "geodetic.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// The grid padded by one node on every side, each pad node at its edge
// neighbor's coordinate and below every level.  The crossing between an edge
// node and its pad then sits on the edge node, so a region reaching the edge
// is closed along it, and every ring closes.
struct Padded
{
    std::size_t nx = 0, ny = 0;
    std::vector<double> x, y;
    std::vector<double> v;

    double at(std::size_t ix, std::size_t iy) const
    {
        return v[iy * nx + ix];
    }
};

Padded pad(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& values)
{
    Padded p;
    p.nx = x.size() + 2;
    p.ny = y.size() + 2;
    p.x.push_back(x.front());
    p.x.insert(p.x.end(), x.begin(), x.end());
    p.x.push_back(x.back());
    p.y.push_back(y.front());
    p.y.insert(p.y.end(), y.begin(), y.end());
    p.y.push_back(y.back());
    p.v.assign(p.nx * p.ny, -std::numeric_limits<double>::infinity());
    for (std::size_t iy = 0; iy < y.size(); iy++)
        for (std::size_t ix = 0; ix < x.size(); ix++)
        {
            const double value = values[iy * x.size() + ix];
            p.v[(iy + 1) * p.nx + ix + 1] = std::isfinite(value) ? value : -std::numeric_limits<double>::infinity();
        }
    return p;
}

// A traced ring made simple: consecutive repeats dropped (points along the
// padded edge, where a pad cell has no width), out-and-back spikes removed
// (a, b, a), and a ring that touches itself at a vertex split there into
// two.  Those arise where a node's value equals the level exactly, so that
// crossings land on the node; they leave the area alone but make a ring no
// GIS accepts.  Pieces of under three points carry no area and are dropped.
std::vector<Ring> simple_rings(Ring ring)
{
    std::vector<Ring> done;
    std::vector<Ring> pending {std::move(ring)};
    while (!pending.empty())
    {
        Ring r = std::move(pending.back());
        pending.pop_back();
        bool changed = true;
        while (changed && r.size() >= 3)
        {
            changed = false;
            for (std::size_t i = 0; i < r.size() && r.size() >= 3; i++)
            {
                const std::size_t next = (i + 1) % r.size();
                const std::size_t after = (i + 2) % r.size();
                if (r[i] == r[next])
                {
                    r.erase(r.begin() + static_cast<std::ptrdiff_t>(next));
                    changed = true;
                }
                else if (r[i] == r[after])
                {
                    // a, b, a: drop b and one a.
                    const std::size_t first = std::min(next, after), second = std::max(next, after);
                    r.erase(r.begin() + static_cast<std::ptrdiff_t>(second));
                    r.erase(r.begin() + static_cast<std::ptrdiff_t>(first));
                    changed = true;
                }
            }
        }
        if (r.size() < 3)
            continue;
        // A vertex met twice: split into the loop between the two meetings
        // and the rest.
        bool split = false;
        for (std::size_t i = 0; i < r.size() && !split; i++)
        {
            for (std::size_t j = i + 1; j < r.size() && !split; j++)
            {
                if (r[i] == r[j])
                {
                    Ring loop(r.begin() + static_cast<std::ptrdiff_t>(i), r.begin() + static_cast<std::ptrdiff_t>(j));
                    Ring rest(r.begin(), r.begin() + static_cast<std::ptrdiff_t>(i));
                    rest.insert(rest.end(), r.begin() + static_cast<std::ptrdiff_t>(j), r.end());
                    pending.push_back(std::move(loop));
                    pending.push_back(std::move(rest));
                    split = true;
                }
            }
        }
        if (!split)
            done.push_back(std::move(r));
    }
    return done;
}

} // namespace

std::vector<ContourPolygon> niceops::contour_polygons(const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& values, double level)
{
    assert(x.size() * y.size() == values.size());
    std::vector<ContourPolygon> polygons;
    if (x.empty() || y.empty())
        return polygons;
    const Padded p = pad(x, y, values);
    const auto high = [&](std::size_t ix, std::size_t iy) { return p.at(ix, iy) >= level; };

    // A crossing is named by the grid edge it lies on: 2 k for the edge from
    // node k to its +x neighbor, 2 k + 1 to its +y neighbor.
    const auto node = [&](std::size_t ix, std::size_t iy) { return iy * p.nx + ix; };
    const auto point_on = [&](std::size_t key) -> std::array<double, 2> {
        const std::size_t k = key / 2;
        const std::size_t ix = k % p.nx, iy = k / p.nx;
        const std::size_t jx = key % 2 == 0 ? ix + 1 : ix;
        const std::size_t jy = key % 2 == 0 ? iy : iy + 1;
        const double a = p.at(ix, iy), b = p.at(jx, jy);
        // Where the level falls between the two values.  Against a node below
        // every level (-inf), the crossing sits on the finite node.
        double t = 0.0;
        if (std::isinf(a))
            t = 1.0;
        else if (!std::isinf(b) && b != a)
            t = (level - a) / (b - a);
        t = std::clamp(t, 0.0, 1.0);
        return {p.x[ix] + t * (p.x[jx] - p.x[ix]), p.y[iy] + t * (p.y[jy] - p.y[iy])};
    };

    // Each cell's segments, oriented with the high region on their left, so
    // that rings around high ground run counterclockwise and those around low
    // ground inside it clockwise.  Walking a cell's boundary counterclockwise,
    // a segment starts where the walk leaves the high region and ends where it
    // next enters it; at a saddle, "next" is the following crossing when the
    // center is high (the high corners join) and the preceding one when not.
    std::unordered_map<std::size_t, std::size_t> next_of; // start crossing -> end crossing
    for (std::size_t iy = 0; iy + 1 < p.ny; iy++)
    {
        for (std::size_t ix = 0; ix + 1 < p.nx; ix++)
        {
            const bool corner[4] = {high(ix, iy), high(ix + 1, iy), high(ix + 1, iy + 1), high(ix, iy + 1)};
            // Edges in counterclockwise order: bottom, right, top, left.
            const std::size_t edge[4] = {2 * node(ix, iy), 2 * node(ix + 1, iy) + 1, 2 * node(ix, iy + 1),
                    2 * node(ix, iy) + 1};
            std::size_t crossing[4];
            bool leaves[4];
            int count = 0;
            for (int e = 0; e < 4; e++)
            {
                const bool from = corner[e], to = corner[(e + 1) % 4];
                if (from != to)
                {
                    crossing[count] = edge[e];
                    leaves[count] = from;
                    count++;
                }
            }
            if (count == 2)
            {
                const int start = leaves[0] ? 0 : 1;
                next_of[crossing[start]] = crossing[1 - start];
            }
            else if (count == 4)
            {
                double sum = 0.0;
                bool finite = true;
                for (const auto& [cx, cy] : {std::array<std::size_t, 2> {ix, iy}, {ix + 1, iy},
                             {ix + 1, iy + 1}, {ix, iy + 1}})
                {
                    finite = finite && std::isfinite(p.at(cx, cy));
                    sum += p.at(cx, cy);
                }
                const bool center_high = finite && 0.25 * sum >= level;
                for (int c = 0; c < 4; c++)
                {
                    if (!leaves[c])
                        continue;
                    const int partner = center_high ? (c + 1) % 4 : (c + 3) % 4;
                    next_of[crossing[c]] = crossing[partner];
                }
            }
        }
    }

    // Link the segments into rings.
    std::vector<Ring> outers, holes;
    while (!next_of.empty())
    {
        const std::size_t first = next_of.begin()->first;
        Ring ring;
        std::size_t at = first;
        while (true)
        {
            const auto found = next_of.find(at);
            if (found == next_of.end())
                break; // cannot happen on a padded grid; drop an open chain
            ring.push_back(point_on(at));
            const std::size_t to = found->second;
            next_of.erase(found);
            at = to;
            if (at == first)
                break;
        }
        for (auto& piece : simple_rings(std::move(ring)))
        {
            const double area = ring_area(piece);
            if (area > 0.0)
                outers.push_back(std::move(piece));
            else if (area < 0.0)
                holes.push_back(std::move(piece));
        }
    }

    // Each hole into the smallest outer ring holding it.
    const auto inside = [](const Ring& ring, const std::array<double, 2>& q) {
        bool in = false;
        for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
        {
            const auto& a = ring[i];
            const auto& b = ring[j];
            if ((a[1] > q[1]) != (b[1] > q[1])
                    && q[0] < (b[0] - a[0]) * (q[1] - a[1]) / (b[1] - a[1]) + a[0])
                in = !in;
        }
        return in;
    };
    polygons.resize(outers.size());
    for (std::size_t k = 0; k < outers.size(); k++)
        polygons[k].outer = std::move(outers[k]);
    for (auto& hole : holes)
    {
        // Placed by a point strictly inside the hole: a hair to the right of
        // the midpoint of its longest edge, a hole running clockwise.  Rings
        // from one level never cross, so the outer rings holding that point
        // are the hole's own and those around it, and the smallest is its
        // own.  The mean of its vertices was not inside a U-shaped hole, and
        // fell inside an island standing in the hole.
        std::size_t longest = 0;
        double longest_length = -1.0;
        for (std::size_t i = 0; i < hole.size(); i++)
        {
            const auto& a = hole[i];
            const auto& b = hole[(i + 1) % hole.size()];
            const double length = std::hypot(b[0] - a[0], b[1] - a[1]);
            if (length > longest_length)
            {
                longest = i;
                longest_length = length;
            }
        }
        const auto& a = hole[longest];
        const auto& b = hole[(longest + 1) % hole.size()];
        const double step = 1e-6;
        const std::array<double, 2> q {0.5 * (a[0] + b[0]) + step * (b[1] - a[1]),
                0.5 * (a[1] + b[1]) - step * (b[0] - a[0])};
        std::size_t best = polygons.size();
        double best_area = 0.0;
        for (std::size_t k = 0; k < polygons.size(); k++)
        {
            const double area = ring_area(polygons[k].outer);
            if (inside(polygons[k].outer, q) && (best == polygons.size() || area < best_area))
            {
                best = k;
                best_area = area;
            }
        }
        // Every hole lies inside some high ground; one that found none would
        // be a defect, and dropping it would overstate the area.
        assert(best < polygons.size());
        if (best < polygons.size())
            polygons[best].holes.push_back(std::move(hole));
    }
    return polygons;
}

double niceops::ring_area(const Ring& ring)
{
    double twice = 0.0;
    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
        twice += (ring[j][0] * ring[i][1]) - (ring[i][0] * ring[j][1]);
    return 0.5 * twice;
}

double niceops::polygon_area(const ContourPolygon& polygon)
{
    double area = ring_area(polygon.outer);
    for (const auto& hole : polygon.holes)
        area += ring_area(hole); // negative
    return area;
}

std::vector<double> niceops::write_contour_geojson(const std::string& path, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& values, const std::vector<double>& levels,
        const LocalFrame& frame, const std::string& metric)
{
    std::vector<double> areas;
    auto features = nlohmann::ordered_json::array();
    const auto to_lon_lat = [&](const Ring& ring) {
        auto coordinates = nlohmann::ordered_json::array();
        for (const auto& point : ring)
        {
            const auto position = frame.to_geodetic(Vec3(point[0], point[1], 0.0));
            coordinates.push_back({position.longitude, position.latitude});
        }
        // GeoJSON closes a ring by repeating its first position.
        coordinates.push_back(coordinates.front());
        return coordinates;
    };
    for (const double level : levels)
    {
        const auto polygons = contour_polygons(x, y, values, level);
        double area = 0.0;
        auto multipolygon = nlohmann::ordered_json::array();
        for (const auto& polygon : polygons)
        {
            area += polygon_area(polygon);
            auto rings = nlohmann::ordered_json::array();
            rings.push_back(to_lon_lat(polygon.outer));
            for (const auto& hole : polygon.holes)
                rings.push_back(to_lon_lat(hole));
            multipolygon.push_back(rings);
        }
        areas.push_back(area);
        nlohmann::ordered_json feature;
        feature["type"] = "Feature";
        feature["properties"] = {{"metric", metric}, {"level_db", level}, {"area_m2", area},
                {"area_km2", area * 1e-6}, {"area_sq_mi", area / 2589988.110336}};
        feature["geometry"] = {{"type", "MultiPolygon"}, {"coordinates", multipolygon}};
        features.push_back(feature);
    }
    nlohmann::ordered_json collection;
    collection["type"] = "FeatureCollection";
    collection["features"] = features;
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out)
        input_error("Cannot write --contour_file ", path, ".");
    out << collection.dump() << '\n';
    if (!out.flush())
        input_error("Could not finish writing --contour_file ", path, ".");
    return areas;
}
