#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
#include "edge_visibility.hpp"

using namespace niceops;

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

} // namespace

// visible_stretch may only say a point of an edge is hidden when it is:
// wherever folded_path_open sees a point of the edge (sampled every 2 cm and
// at its ends), the point must lie within the stretch it returns.  Over every
// image to two reflections and every edge each keeps, with and without the
// floor, from points among DigiCity's containers, a block of a district with
// an L-shaped building, a street of gabled houses, and U-, T- and L-shaped
// buildings (one pitched, one with a vertex in line with its neighbours),
// which prove through their convex pieces.  It must also prove most of the edges no sample
// sees wholly hidden, and most of the unseen samples of the rest, or it is no
// use: for the point and its images of one reflection, and for those of two,
// each on its own.
TEST_CASE("a point of an edge proved hidden is hidden", "[point_paths]")
{
    Footprint ell;
    ell.x = {150.0, 190.0, 190.0, 170.0, 170.0, 150.0};
    ell.y = {0.0, 0.0, 20.0, 20.0, 50.0, 50.0};
    ell.top = 25.0;
    const Buildings containers({box(4.09, 17.66, 6.53, 23.72, 7.77), box(14.81, 17.66, 17.25, 23.72, 7.77),
            box(1.6, 10.82, 4.04, 16.88, 7.77), box(1.6, 8.38, 13.72, 10.82, 5.18),
            box(20.21, 3.86, 22.65, 16.0, 7.77), box(13.5, -1.84, 19.56, 0.6, 7.77)});
    const Buildings blocks({box(0.0, 0.0, 60.0, 60.0, 20.0), box(80.0, 0.0, 140.0, 60.0, 35.0),
            box(0.0, 80.0, 60.0, 140.0, 10.0), box(80.0, 80.0, 140.0, 140.0, 40.0), ell});
    // Two rows of houses across a street, gabled along it, eaves at 6 m and
    // ridges at 9 m and 10.5 m.
    auto gabled = [](double x0, double y0, double x1, double y1, double eaves, double ridge) {
        Footprint f = box(x0, y0, x1, y1, ridge);
        const double slope = (ridge - eaves) / (0.5 * (y1 - y0));
        f.roof = {{0.0, slope, eaves - slope * y0}, {0.0, -slope, eaves + slope * y1}};
        return f;
    };
    const Buildings street({gabled(0.0, 0.0, 10.0, 12.0, 6.0, 9.0), gabled(12.0, 0.0, 22.0, 12.0, 6.0, 10.5),
            gabled(25.0, 0.0, 33.0, 12.0, 6.0, 9.0), gabled(0.0, 22.0, 12.0, 34.0, 6.0, 10.5),
            gabled(14.0, 22.0, 22.0, 34.0, 6.0, 9.0), gabled(24.0, 22.0, 34.0, 34.0, 6.0, 10.5)});
    auto ring = [](std::vector<double> x, std::vector<double> y, double top) {
        Footprint f;
        f.x = std::move(x);
        f.y = std::move(y);
        f.top = top;
        return f;
    };
    // An L gabled along x, ridge 12 m over y = 80, hipped at its west end.
    Footprint pitched = ring({200.0, 240.0, 240.0, 215.0, 215.0, 200.0}, {60.0, 60.0, 75.0, 75.0, 100.0, 100.0}, 12.0);
    pitched.roof = {{0.0, 0.3, 6.0 - 0.3 * 60.0}, {0.0, -0.3, 6.0 + 0.3 * 100.0}, {0.3, 0.0, 6.0 - 0.3 * 200.0}};
    const Buildings shapes({ring({200.0, 240.0, 240.0, 230.0, 230.0, 210.0, 210.0, 200.0},
                                   {0.0, 0.0, 40.0, 40.0, 10.0, 10.0, 40.0, 40.0}, 15.0),
            ring({275.0, 285.0, 285.0, 300.0, 300.0, 260.0, 260.0, 275.0}, {0.0, 0.0, 30.0, 30.0, 40.0, 40.0, 30.0, 30.0},
                    20.0),
            pitched,
            ring({300.0, 340.0, 340.0, 320.0, 320.0, 310.0, 300.0}, {60.0, 60.0, 75.0, 75.0, 100.0, 100.0, 100.0},
                    10.0)});
    struct Case
    {
        const Buildings* city;
        Vec3 point;
    };
    const std::vector<Case> cases {
        {&containers, {10.5, 20.5, 2.0}},
        {&containers, {10.5, -1.1, 1.2}},
        {&containers, {8.0, 12.0, 15.0}},
        {&containers, {-3.6, 6.4, 1.2}},
        {&blocks, {70.0, 70.0, 30.0}},
        {&blocks, {70.0, 30.0, 1.5}},
        {&blocks, {145.0, 70.0, 5.0}},
        {&blocks, {160.0, 60.0, 2.0}},
        {&street, {16.0, 17.0, 1.5}},
        {&street, {5.0, 40.0, 15.0}},
        {&street, {23.0, 6.0, 12.0}},
        {&shapes, {220.0, 30.0, 1.5}},
        {&shapes, {250.0, 50.0, 3.0}},
        {&shapes, {290.0, 10.0, 25.0}},
        {&shapes, {270.0, 80.0, 1.2}},
        {&shapes, {330.0, 90.0, 4.0}},
    };
    // Per order of the image: the point itself and its images of one
    // reflection, then those of two.
    std::size_t visible[2] = {}, unseen[2] = {}, proved[2] = {}, points_unseen[2] = {}, points_proved[2] = {};
    std::size_t wrong = 0;
    std::vector<Vec3> hits;
    for (const Case& c : cases)
    {
        const std::vector<std::size_t> near = edges_near(*c.city, c.point, 150.0);
        const ImageTree tree = image_tree(*c.city, c.point, 150.0, 2, 20000, true, near);
        for (std::size_t i = 0; i < tree.images.size(); i++)
        {
            const std::size_t order = std::max<std::size_t>(tree.images[i].facades.size(), 1);
            for (std::size_t e : tree.edges[i])
                for (int floor = 0; floor < 2; floor++)
                {
                    const Edge& edge = c.city->edges()[e];
                    double from, to;
                    const bool maybe
                            = visible_stretch(*c.city, edge, c.point, tree.images[i].facades, floor == 1, from, to);
                    bool seen = false;
                    const int samples = std::max(2, static_cast<int>(edge.length / 0.02));
                    for (int k = 0; k <= samples; k++)
                    {
                        const double t = edge.length * k / samples;
                        const Vec3 q = edge.start + t * edge.direction;
                        if (folded_path_open(*c.city, q, c.point, tree.images[i].facades, floor == 1, &hits))
                        {
                            seen = true;
                            wrong += !maybe || t < from || t > to;
                        }
                        else if (maybe)
                        {
                            points_unseen[order - 1]++;
                            points_proved[order - 1] += t < from || t > to;
                        }
                    }
                    if (seen)
                        visible[order - 1]++;
                    else
                    {
                        unseen[order - 1]++;
                        proved[order - 1] += !maybe;
                    }
                }
        }
    }
    REQUIRE(wrong == 0);
    for (int o = 0; o < 2; o++)
    {
        INFO((o == 0 ? "the point and its images of one reflection: " : "images of two reflections: ")
                << visible[o] << " seen, " << unseen[o] << " unseen, " << proved[o] << " of those proved hidden; of the points unseen on edges not proved hidden, "
                << points_proved[o] << " of " << points_unseen[o] << " proved");
        CHECK(visible[o] > 100);
        CHECK(proved[o] > unseen[o] / 2);
        CHECK(points_proved[o] > points_unseen[o] / 2);
    }
}

// A footprint that is not convex is cut into convex pieces that tile it: each
// counterclockwise and turning left throughout, its vertices the footprint's,
// their areas summing to the footprint's, and with its top and roof.  One
// that is convex, or has a courtyard, is not cut.
TEST_CASE("a footprint that is not convex is cut into convex pieces", "[buildings]")
{
    auto ring = [](std::vector<double> x, std::vector<double> y) {
        Footprint f;
        f.x = std::move(x);
        f.y = std::move(y);
        f.top = 10.0;
        return f;
    };
    auto area = [](const Footprint& f) {
        double a = 0.0;
        for (std::size_t i = 0; i < f.x.size(); i++)
        {
            const std::size_t j = (i + 1) % f.x.size();
            a += f.x[i] * f.y[j] - f.x[j] * f.y[i];
        }
        return 0.5 * a;
    };
    Footprint court = box(100.0, 100.0, 140.0, 140.0, 10.0);
    court.holes.push_back({{110.0, 110.0, 130.0, 130.0}, {110.0, 130.0, 130.0, 110.0}});
    Footprint ell = ring({0.0, 40.0, 40.0, 15.0, 15.0, 0.0}, {0.0, 0.0, 15.0, 15.0, 40.0, 40.0});
    ell.top = 20.0;
    ell.roof = {{0.0, 0.3, 6.0}, {0.3, 0.0, 6.0}};
    const Buildings city({ell, ring({50.0, 90.0, 90.0, 80.0, 80.0, 60.0, 60.0, 50.0}, {0.0, 0.0, 40.0, 40.0, 10.0, 10.0, 40.0, 40.0}),
            ring({0.0, 40.0, 40.0, 20.0, 20.0, 10.0, 0.0}, {60.0, 60.0, 75.0, 75.0, 100.0, 100.0, 100.0}),
            box(60.0, 60.0, 80.0, 80.0, 10.0), court});
    for (std::size_t b = 0; b < city.size(); b++)
    {
        const Footprint& f = city.footprints()[b];
        const std::vector<Footprint>& pieces = city.convex_pieces(b);
        INFO("building " << b << ", " << pieces.size() << " pieces");
        if (b >= 3)
        {
            CHECK(pieces.empty());
            continue;
        }
        REQUIRE(pieces.size() >= 2);
        double total = 0.0;
        for (const Footprint& p : pieces)
        {
            const std::size_t n = p.x.size();
            for (std::size_t i = 0; i < n; i++)
            {
                const std::size_t j = (i + 1) % n, k = (i + 2) % n;
                CHECK((p.x[j] - p.x[i]) * (p.y[k] - p.y[j]) - (p.y[j] - p.y[i]) * (p.x[k] - p.x[j]) >= 0.0);
                bool vertex = false;
                for (std::size_t v = 0; v < f.x.size(); v++)
                    vertex = vertex || (p.x[i] == f.x[v] && p.y[i] == f.y[v]);
                CHECK(vertex);
            }
            CHECK(p.top == f.top);
            CHECK(p.roof.size() == f.roof.size());
            total += area(p);
        }
        CHECK(std::abs(total - area(f)) < 1e-9 * area(f));
    }
}
