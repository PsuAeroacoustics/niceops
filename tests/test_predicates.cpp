// New tests (no D original): the exact geometric predicates the Delaunay
// triangulation rests on.  The cases that matter are the degenerate ones --
// exactly collinear, exactly cospherical -- because those are what the noise
// database is full of and what floating point cannot answer.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "predicates.hpp"

using namespace niceops;

TEST_CASE("Determinant signs", "[predicates]")
{
    const double identity[4] = {1, 0, 0, 1};
    REQUIRE(exact_determinant_sign(identity, 2) == 1);

    const double swapped[4] = {0, 1, 1, 0};
    REQUIRE(exact_determinant_sign(swapped, 2) == -1);

    const double singular[4] = {1, 2, 2, 4};
    REQUIRE(exact_determinant_sign(singular, 2) == 0);

    // 3x3 with a known value of -306.
    const double three[9] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    REQUIRE(exact_determinant_sign(three, 3) == -1);

    // A singular 3x3: the third row is the sum of the first two.
    const double dependent[9] = {1, 2, 3, 4, 5, 6, 5, 7, 9};
    REQUIRE(exact_determinant_sign(dependent, 3) == 0);
}

TEST_CASE("Orientation of a simplex", "[predicates]")
{
    // 2D: counter-clockwise is positive, clockwise negative, collinear zero.
    const double ccw[6] = {0, 0, 1, 0, 0, 1};
    REQUIRE(orientation(ccw, 2) == 1);
    const double cw[6] = {0, 0, 0, 1, 1, 0};
    REQUIRE(orientation(cw, 2) == -1);
    const double flat[6] = {0, 0, 1, 1, 2, 2};
    REQUIRE(orientation(flat, 2) == 0);

    // 3D: a unit tetrahedron, and a flat one.
    const double tetra[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    REQUIRE(orientation(tetra, 3) != 0);
    const double coplanar[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0};
    REQUIRE(orientation(coplanar, 3) == 0);
}

TEST_CASE("Orientation is exact where double precision is not", "[predicates]")
{
    // Shewchuk's canonical near-degenerate case: three points on the line y = x.
    // Exactly collinear, so the true determinant is zero.
    const double collinear[6] = {0.5, 0.5, 12.0, 12.0, 24.0, 24.0};
    REQUIRE(orientation(collinear, 2) == 0);

    // Nudge one point by a single unit in the last place.  The true answer is
    // now definitively one sign, by an amount far below what evaluating the
    // determinant in double precision could resolve.
    double nudged_up[6] = {std::nextafter(0.5, 1.0), 0.5, 12.0, 12.0, 24.0, 24.0};
    double nudged_down[6] = {std::nextafter(0.5, 0.0), 0.5, 12.0, 12.0, 24.0, 24.0};
    const int up = orientation(nudged_up, 2);
    const int down = orientation(nudged_down, 2);
    REQUIRE(up != 0);
    REQUIRE(down != 0);
    REQUIRE(up == -down);
}

TEST_CASE("In-sphere on the cube's circumsphere", "[predicates]")
{
    // Four corners of the unit cube.  Their circumsphere is centered on the cube
    // center and passes through all eight corners -- so every other corner is
    // *exactly* on it, and the predicate must say so rather than guess.
    const double simplex[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

    const double on_it[3] = {1, 1, 1};
    REQUIRE(in_sphere(simplex, on_it, 3) == 0);
    const double also_on_it[3] = {1, 1, 0};
    REQUIRE(in_sphere(simplex, also_on_it, 3) == 0);
    const double and_again[3] = {0, 1, 1};
    REQUIRE(in_sphere(simplex, and_again, 3) == 0);

    // The center is inside, a distant point is outside.
    const double center[3] = {0.5, 0.5, 0.5};
    REQUIRE(in_sphere(simplex, center, 3) == 1);
    const double far[3] = {2, 2, 2};
    REQUIRE(in_sphere(simplex, far, 3) == -1);

    // Just inside and just outside the circumsphere, along a diagonal.
    const double just_in[3] = {0.999, 0.999, 0.999};
    REQUIRE(in_sphere(simplex, just_in, 3) == 1);
    const double just_out[3] = {1.001, 1.001, 1.001};
    REQUIRE(in_sphere(simplex, just_out, 3) == -1);
}

TEST_CASE("In-sphere does not depend on vertex order", "[predicates]")
{
    // Swapping two vertices reverses the simplex's orientation; the answer to
    // "is this point inside the circumsphere" must not change with it.
    const double simplex[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    const double swapped[12] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1};
    const double center[3] = {0.5, 0.5, 0.5};
    const double far[3] = {5, 5, 5};
    REQUIRE(in_sphere(simplex, center, 3) == in_sphere(swapped, center, 3));
    REQUIRE(in_sphere(simplex, far, 3) == in_sphere(swapped, far, 3));
}

TEST_CASE("In-sphere works in two dimensions", "[predicates]")
{
    // The unit square's corners are concyclic.
    const double triangle[6] = {0, 0, 1, 0, 0, 1};
    const double fourth_corner[2] = {1, 1};
    REQUIRE(in_sphere(triangle, fourth_corner, 2) == 0);
    const double center[2] = {0.5, 0.5};
    REQUIRE(in_sphere(triangle, center, 2) == 1);
    const double outside[2] = {3, 3};
    REQUIRE(in_sphere(triangle, outside, 2) == -1);
}

TEST_CASE("The prism structure of the database is exactly cospherical", "[predicates]")
{
    // The shape the real data takes: a flight condition repeated at several
    // thrust coefficients.  Two conditions at two CT levels put four points on
    // a rectangle in a plane, and any fourth is on the circumcircle of the
    // other three -- exact ties, not near ones.
    const double a[3] = {0.20, -4.0, 0.004};
    const double b[3] = {0.30, -4.0, 0.004};
    const double c[3] = {0.20, -4.0, 0.006};
    const double d[3] = {0.30, -4.0, 0.006};

    // All four are coplanar (gamma is constant), so the tetrahedron on them is
    // degenerate and the predicate reports it rather than inventing a sign.
    const double flat[12] = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2],
            d[0], d[1], d[2]};
    REQUIRE(orientation(flat, 3) == 0);
}

TEST_CASE("Predicates stay self-consistent near degeneracy", "[predicates]")
{
    // The filtered fast path must never disagree with the exact one.  Swapping
    // two vertices reverses a simplex's orientation exactly, so orientation and
    // its swap must be exact negatives -- including for configurations chosen to
    // sit right on the edge of what double precision can resolve.
    std::mt19937 rng(20260918);
    std::uniform_int_distribution<int> small(-3, 3);

    for (int trial = 0; trial < 3000; trial++)
    {
        // Points on the plane z = x + y + 1, then one coordinate nudged by a
        // few ulps, so the true determinant is tiny but definite.  The offset
        // keeps every coordinate a normal number: the exact arithmetic assumes
        // no underflow, and nudging an exact zero lands in the subnormals where
        // that assumption does not hold.
        double points[12];
        for (int i = 0; i < 4; i++)
        {
            points[i * 3 + 0] = 0.1 * small(rng);
            points[i * 3 + 1] = 0.1 * small(rng);
            points[i * 3 + 2] = points[i * 3 + 0] + points[i * 3 + 1] + 1.0;
        }
        const int nudges = small(rng);
        for (int n = 0; n < std::abs(nudges); n++)
            points[11] = std::nextafter(points[11], nudges > 0 ? 1e9 : -1e9);

        const int direct = orientation(points, 3);

        double swapped[12];
        std::copy(std::begin(points), std::end(points), std::begin(swapped));
        for (int j = 0; j < 3; j++)
            std::swap(swapped[j], swapped[3 + j]);
        const int reversed = orientation(swapped, 3);

        REQUIRE(direct == -reversed);
    }
}

namespace {

// Sign of a determinant of small integers, exactly: fraction-free (Bareiss)
// elimination in 128-bit integers, so an oracle with nothing in common with the
// expansion arithmetic under test.
int integer_determinant_sign(std::vector<std::vector<__int128>> m)
{
    const std::size_t n = m.size();
    int sign = 1;
    __int128 previous = 1;
    for (std::size_t k = 0; k + 1 < n; k++)
    {
        if (m[k][k] == 0)
        {
            std::size_t swap = k + 1;
            while (swap < n && m[swap][k] == 0)
                swap++;
            if (swap == n)
                return 0;
            std::swap(m[k], m[swap]);
            sign = -sign;
        }
        for (std::size_t i = k + 1; i < n; i++)
            for (std::size_t j = k + 1; j < n; j++)
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) / previous;
        previous = m[k][k];
    }
    const __int128 d = m[n - 1][n - 1];
    return d == 0 ? 0 : (d > 0 ? sign : -sign);
}

// The lifted determinant of `dim + 1` points about a query: rows (a - q, |a - q|^2).
int lifted_sign(const std::vector<std::vector<long>>& points, const std::vector<long>& query)
{
    const std::size_t dim = query.size();
    std::vector<std::vector<__int128>> m;
    for (const auto& a : points)
    {
        std::vector<__int128> row(dim + 1);
        __int128 norm = 0;
        for (std::size_t j = 0; j < dim; j++)
        {
            row[j] = a[j] - query[j];
            norm += row[j] * row[j];
        }
        row[dim] = norm;
        m.push_back(row);
    }
    return integer_determinant_sign(m);
}

} // namespace

TEST_CASE("in_sphere agrees with an exact integer oracle, in one to five dimensions", "[predicates]")
{
    // The oracle needs no sign convention: the lifted determinant is affine in the
    // lifted query and zero on the circumsphere, so it has one sign inside and the
    // other outside, and the simplex's own centroid (scaled to stay integer) is
    // always inside.  Cube corners and a coarse lattice make exact ties common; the
    // dimension-four and -five paths were not pinned directly before, and a wrong
    // in_sphere there built a non-Delaunay mesh the Delaunay tests could not see.
    std::mt19937_64 random(2024);
    int ties = 0, inside = 0, outside = 0;
    for (std::size_t dim = 1; dim <= 5; dim++)
    {
        const long scale = static_cast<long>(dim) + 1;
        for (int trial = 0; trial < 1500; trial++)
        {
            std::vector<std::vector<long>> points(dim + 1, std::vector<long>(dim));
            std::vector<double> flat;
            const bool cube = trial % 3 == 0;
            const long range = trial % 3 == 1 ? 2 : 4;
            do
            {
                flat.clear();
                for (auto& p : points)
                    for (auto& x : p)
                    {
                        x = cube ? static_cast<long>(random() % 2) : static_cast<long>(random() % (2 * range + 1)) - range;
                    }
            } while ([&] {
                std::vector<double> test;
                for (const auto& p : points)
                    for (long x : p)
                        test.push_back(static_cast<double>(x));
                return orientation(test.data(), dim) == 0;
            }());
            std::vector<long> query(dim);
            for (auto& x : query)
                x = cube ? static_cast<long>(random() % 2) : static_cast<long>(random() % (2 * range + 3)) - range - 1;
            // Scale by dim + 1 so the centroid is an integer point.
            std::vector<std::vector<long>> scaled = points;
            std::vector<long> centroid(dim, 0), scaled_query(dim);
            for (auto& p : scaled)
                for (std::size_t j = 0; j < dim; j++)
                {
                    centroid[j] += p[j];
                    p[j] *= scale;
                }
            for (std::size_t j = 0; j < dim; j++)
                scaled_query[j] = query[j] * scale;
            const int at_query = lifted_sign(scaled, scaled_query);
            const int at_centroid = lifted_sign(scaled, centroid);
            REQUIRE(at_centroid != 0);
            const int expected = at_query * at_centroid; // +1 inside, -1 outside, 0 on it
            std::vector<double> simplex, q;
            for (const auto& p : points)
                for (long x : p)
                    simplex.push_back(static_cast<double>(x));
            for (long x : query)
                q.push_back(static_cast<double>(x));
            INFO("dim " << dim << " trial " << trial);
            REQUIRE(in_sphere(simplex.data(), q.data(), dim) == expected);
            ties += expected == 0;
            inside += expected > 0;
            outside += expected < 0;
        }
    }
    // The cases include all three outcomes in quantity.
    CHECK(ties > 100);
    CHECK(inside > 500);
    CHECK(outside > 500);
}
