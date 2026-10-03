// New tests (no D original): the Delaunay triangulation that replaces CGAL.
// The interesting cases are the degenerate ones -- a regular lattice, and the
// prism the noise database's operating conditions actually form.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "delaunay.hpp"
#include "predicates.hpp"

using namespace niceops;

namespace {

// Check that the weights reconstruct the query point and form a partition of
// unity -- the two properties everything downstream relies on.
void check_barycentric(const Delaunay& mesh, const std::vector<double>& points,
        const std::vector<double>& query, double tolerance = 1e-9)
{
    const std::size_t dim = mesh.dimension();
    std::vector<double> weights(dim + 1, 0.0);
    std::vector<std::size_t> indices(dim + 1, 0);
    REQUIRE(mesh.barycentric(query.data(), weights.data(), indices.data()));

    double sum = 0.0;
    for (double w : weights)
    {
        REQUIRE(w >= -tolerance);
        sum += w;
    }
    REQUIRE(std::fabs(sum - 1.0) < tolerance);

    for (std::size_t j = 0; j < dim; j++)
    {
        double reconstructed = 0.0;
        for (std::size_t v = 0; v <= dim; v++)
            reconstructed += weights[v] * points[indices[v] * dim + j];
        REQUIRE(std::fabs(reconstructed - query[j]) < tolerance);
    }
}

} // namespace

TEST_CASE("A single tetrahedron", "[delaunay]")
{
    const std::vector<double> points {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    Delaunay mesh(3, 4, points.data());
    REQUIRE(mesh.point_count() == 4);
    REQUIRE(mesh.simplex_count() == 1);

    check_barycentric(mesh, points, {0.1, 0.1, 0.1});
    check_barycentric(mesh, points, {0.25, 0.25, 0.25});

    // A vertex is reproduced with all the weight on itself.
    std::vector<double> weights(4, 0.0);
    std::vector<std::size_t> indices(4, 0);
    const std::vector<double> corner {0.0, 0.0, 0.0};
    REQUIRE(mesh.barycentric(corner.data(), weights.data(), indices.data()));
    for (std::size_t v = 0; v < 4; v++)
    {
        if (indices[v] == 0)
            REQUIRE(std::fabs(weights[v] - 1.0) < 1e-12);
        else
            REQUIRE(std::fabs(weights[v]) < 1e-12);
    }

    // Outside the hull is reported rather than guessed at.
    const std::vector<double> outside {2.0, 2.0, 2.0};
    REQUIRE(!mesh.barycentric(outside.data(), weights.data(), indices.data()));
}

TEST_CASE("A regular lattice, which is maximally degenerate", "[delaunay]")
{
    // Every cube of eight lattice points is cospherical, so the triangulation
    // is genuinely non-unique.  Any consistent choice is correct; what must not
    // happen is an invalid mesh or a non-terminating insertion.
    std::vector<double> points;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
                points.push_back(i);
                points.push_back(j);
                points.push_back(k);
            }
    Delaunay mesh(3, points.size() / 3, points.data());
    REQUIRE(mesh.point_count() == 64);
    REQUIRE(mesh.simplex_count() > 0);

    // Interior queries, including ones sitting exactly on lattice planes and
    // exactly at cube centers.
    check_barycentric(mesh, points, {1.5, 1.5, 1.5}, 1e-9);
    check_barycentric(mesh, points, {1.0, 1.5, 2.0}, 1e-9);
    check_barycentric(mesh, points, {2.25, 1.75, 0.5}, 1e-9);
    check_barycentric(mesh, points, {1.0, 1.0, 1.0}, 1e-9);
}

TEST_CASE("The prism the operating conditions form", "[delaunay]")
{
    // Scattered flight conditions repeated at a regular ladder of thrust
    // coefficients: the shape of the real database, and heavily cospherical.
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> conditions;
    for (int i = 0; i < 20; i++)
    {
        conditions.push_back(unit(rng));
        conditions.push_back(unit(rng));
    }
    const std::vector<double> thrust {0.0, 0.25, 0.5, 0.75, 1.0};

    std::vector<double> points;
    for (std::size_t c = 0; c < conditions.size(); c += 2)
    {
        for (double ct : thrust)
        {
            points.push_back(conditions[c]);
            points.push_back(conditions[c + 1]);
            points.push_back(ct);
        }
    }
    Delaunay mesh(3, points.size() / 3, points.data());
    REQUIRE(mesh.point_count() == 100);
    REQUIRE(mesh.simplex_count() > 0);

    // Query at points that lie exactly on the thrust ladder, where the
    // degeneracy is worst.
    std::size_t located = 0;
    std::vector<double> weights(4, 0.0);
    std::vector<std::size_t> indices(4, 0);
    for (double ct : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        for (int t = 0; t < 8; t++)
        {
            const std::vector<double> query {0.3 + 0.05 * t, 0.4 + 0.02 * t, ct};
            if (mesh.barycentric(query.data(), weights.data(), indices.data()))
            {
                located++;
                double sum = 0.0;
                for (double w : weights)
                {
                    REQUIRE(w >= -1e-9);
                    sum += w;
                }
                REQUIRE(std::fabs(sum - 1.0) < 1e-9);
            }
        }
    }
    REQUIRE(located > 0);
}

TEST_CASE("Random points in three dimensions", "[delaunay]")
{
    std::mt19937 rng(987);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> points;
    for (int i = 0; i < 200; i++)
        for (int j = 0; j < 3; j++)
            points.push_back(unit(rng));
    Delaunay mesh(3, 200, points.data());

    // Every interior query must land in a simplex and reconstruct exactly.
    std::uniform_real_distribution<double> inner(0.35, 0.65);
    for (int trial = 0; trial < 40; trial++)
    {
        const std::vector<double> query {inner(rng), inner(rng), inner(rng)};
        check_barycentric(mesh, points, query, 1e-9);
    }
}

TEST_CASE("Two dimensions and four dimensions", "[delaunay]")
{
    // The dimension is a parameter, not a template argument, so a fourth
    // operating variable would cost nothing here.
    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (std::size_t dim : {2u, 4u})
    {
        std::vector<double> points;
        const std::size_t count = dim == 2 ? 60 : 120;
        for (std::size_t i = 0; i < count * dim; i++)
            points.push_back(unit(rng));
        Delaunay mesh(dim, count, points.data());
        REQUIRE(mesh.dimension() == dim);

        std::uniform_real_distribution<double> inner(0.4, 0.6);
        int located = 0;
        for (int trial = 0; trial < 10; trial++)
        {
            std::vector<double> query(dim);
            for (std::size_t j = 0; j < dim; j++)
                query[j] = inner(rng);
            std::vector<double> weights(dim + 1, 0.0);
            std::vector<std::size_t> indices(dim + 1, 0);
            if (mesh.barycentric(query.data(), weights.data(), indices.data()))
            {
                located++;
                double sum = 0.0;
                for (double w : weights)
                    sum += w;
                REQUIRE(std::fabs(sum - 1.0) < 1e-9);
                for (std::size_t j = 0; j < dim; j++)
                {
                    double reconstructed = 0.0;
                    for (std::size_t v = 0; v <= dim; v++)
                        reconstructed += weights[v] * points[indices[v] * dim + j];
                    REQUIRE(std::fabs(reconstructed - query[j]) < 1e-9);
                }
            }
        }
        REQUIRE(located > 0);
    }
}

TEST_CASE("Nearest vertex", "[delaunay]")
{
    const std::vector<double> points {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    Delaunay mesh(3, 4, points.data());
    const std::vector<double> near_origin {0.05, 0.05, 0.05};
    REQUIRE(mesh.nearest_vertex(near_origin.data()) == 0);
    const std::vector<double> far_x {5.0, 0.0, 0.0};
    REQUIRE(mesh.nearest_vertex(far_x.data()) == 1);
}

namespace {

// The Delaunay property itself, through the public interface: whatever simplex a
// query lands in, no input point lies strictly inside its circumsphere.  The other
// tests check that a located query reconstructs and that the weights sum to one,
// which any tiling of the hull satisfies; this is what makes it Delaunay.
void check_empty_spheres(std::size_t dim, const std::vector<double>& points, std::size_t queries, unsigned seed)
{
    const std::size_t count = points.size() / dim;
    const Delaunay mesh(dim, count, points.data());
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> query(dim), weights(dim + 1), simplex((dim + 1) * dim);
    std::vector<std::size_t> indices(dim + 1);
    std::size_t located = 0;
    for (std::size_t q = 0; q < queries; q++)
    {
        // A random convex combination of input points: always inside the hull.
        std::fill(query.begin(), query.end(), 0.0);
        double total = 0.0;
        for (std::size_t k = 0; k <= dim; k++)
        {
            const std::size_t p = random() % count;
            const double w = unit(random);
            total += w;
            for (std::size_t j = 0; j < dim; j++)
                query[j] += w * points[p * dim + j];
        }
        for (double& x : query)
            x /= total;
        if (!mesh.barycentric(query.data(), weights.data(), indices.data()))
            continue;
        located++;
        for (std::size_t v = 0; v <= dim; v++)
            for (std::size_t j = 0; j < dim; j++)
                simplex[v * dim + j] = points[indices[v] * dim + j];
        for (std::size_t p = 0; p < count; p++)
        {
            if (std::find(indices.begin(), indices.end(), p) != indices.end())
                continue;
            INFO("dim " << dim << ", query " << q << ", point " << p);
            REQUIRE(in_sphere(simplex.data(), points.data() + p * dim, dim) <= 0);
        }
    }
    REQUIRE(located > queries * 9 / 10);
}

} // namespace

TEST_CASE("every simplex has an empty circumsphere, in two to four dimensions", "[delaunay]")
{
    std::mt19937_64 random(5);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (const std::size_t dim : {std::size_t(2), std::size_t(3), std::size_t(4)})
    {
        std::vector<double> points(40 * dim);
        for (double& x : points)
            x = unit(random);
        check_empty_spheres(dim, points, 400, 100 + static_cast<unsigned>(dim));
    }
    // The prism the operating conditions form, and a lattice: whole families of
    // exactly cospherical points, where "empty" has to allow points on the sphere.
    std::vector<double> prism;
    for (int a = 0; a < 6; a++)
        for (int b = 0; b < 5; b++)
            for (int c = 0; c < 4; c++)
            {
                prism.push_back(a / 5.0);
                prism.push_back(b / 4.0);
                prism.push_back(c / 3.0);
            }
    check_empty_spheres(3, prism, 400, 7);
    std::vector<double> lattice;
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            for (int c = 0; c < 3; c++)
                for (int d = 0; d < 3; d++)
                {
                    lattice.push_back(a);
                    lattice.push_back(b);
                    lattice.push_back(c);
                    lattice.push_back(d);
                }
    check_empty_spheres(4, lattice, 300, 9);
}
