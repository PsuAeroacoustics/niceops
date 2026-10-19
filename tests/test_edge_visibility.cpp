#include <catch2/catch_test_macros.hpp>

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
// image to one reflection and every edge each keeps, with and without the
// floor, from points among DigiCity's containers and a block of a district
// with an L-shaped building (not convex, so never used to prove anything).
// It must also prove most of the edges no sample sees wholly hidden, and most
// of the unseen samples of the rest, or it is no use.
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
    };
    std::size_t visible = 0, unseen = 0, proved = 0, wrong = 0, points_unseen = 0, points_proved = 0;
    std::vector<Vec3> hits;
    for (const Case& c : cases)
    {
        const std::vector<std::size_t> near = edges_near(*c.city, c.point, 150.0);
        const ImageTree tree = image_tree(*c.city, c.point, 150.0, 1, 20000, true, near);
        for (std::size_t i = 0; i < tree.images.size(); i++)
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
                            points_unseen++;
                            points_proved += t < from || t > to;
                        }
                    }
                    if (seen)
                        visible++;
                    else
                    {
                        unseen++;
                        proved += !maybe;
                    }
                }
    }
    INFO(visible << " seen, " << unseen << " unseen, " << proved << " of those proved hidden; of the points unseen on "
                  "edges not proved hidden, " << points_proved << " of " << points_unseen << " proved");
    REQUIRE(visible > 100);
    REQUIRE(wrong == 0);
    REQUIRE(proved > unseen / 2);
    REQUIRE(points_proved > points_unseen / 2);
}
