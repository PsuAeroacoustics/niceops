#include <array>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <utility>
#include <vector>
#include <numeric>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
//#include <iostream>


// Hide CGAL triangulation in an anonymous namespace
namespace {
    // Define kernel
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
    typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
    // Define 3D triangulation
    typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
    // Handy typedefs for CGAL cell handles and points
    typedef Triangulation::Cell_handle Cell_handle;
    typedef Triangulation::Point Point;
    // Instantiate triangulation
    Triangulation triangles = {};

    template<typename T>
    bool approx(T x, T y, T epsilon = 10 * std::numeric_limits<T>::epsilon()) {
        static_assert(std::is_floating_point<T>::value, "Only valid for floating point types.");
        return std::fabs(x - y) < epsilon;
    }

    template<typename T>
    T extend(T x, T dmin, T dmax) {
        return approx(x, 1.0) ? dmax : (approx(x, 0.0) ? dmin : x);
    }

    template<typename T>
    T clamp(T x) {
        return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
    }

    template<typename T>
    std::array<T, 3> cross(std::array<T, 3> a, std::array<T, 3> b) {
        return std::array<T, 3>{a[1] * b[2] - a[2] * b[1],
                                a[2] * b[0] - a[0] * b[2],
                                a[0] * b[1] - a[1] * b[0]};
    }

    template<typename T>
    T dot(std::array<T, 3> a, std::array<T, 3> b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

}

extern "C"
{
int set_interpolator3(size_t n, const double *x1, const double *x2, const double *x3, double extend_domain_by = 0.1) {
    std::vector<std::pair<Point, size_t>> points = {};
    points.reserve(n);
    auto dmax = 1.0 + extend_domain_by;
    auto dmin = 0.0 - extend_domain_by;

    for (size_t i = 0; i < n; i++) {
        auto xx1 = x1[i];
        auto xx2 = x2[i];
        auto xx3 = x3[i];
        points.push_back(std::make_pair(Point(xx1, xx2, xx3), i));
        //If point is on edge of data range, then add another extended point
        //at the domain limits with the same data index

        if (approx(xx1, 0.0) || approx(xx2, 0.0) || approx(xx3, 0.0) ||
            approx(xx1, 1.0) || approx(xx2, 1.0) || approx(xx3, 1.0)) {
            auto extended_point = Point(extend(xx1, dmin, dmax), extend(xx2, dmin, dmax), extend(xx3, dmin, dmax));
            points.push_back(std::make_pair(extended_point, i));
        }
    }
    triangles = Triangulation(points.begin(), points.end());
    return !triangles.is_valid();
}

void barycentric(double x1, double x2, double x3, double *weights, size_t *indicies) {
    // Clamp values between 0.0 and 1.0
    x1 = clamp(x1);
    x2 = clamp(x2);
    x3 = clamp(x3);

    // Store a hint about the last cell found
    static thread_local Cell_handle hint;
    auto request = Point{x1, x2, x3};
    //Find tetrahedron that encloses (x1,x2,x3)
    Triangulation::Locate_type locate_type;
    int li = 0, lj = 0;
    auto tetrahedron = triangles.locate(request, locate_type, li, lj, hint);
    assert(tetrahedron->is_valid());
    hint = tetrahedron;

    // The query coincides exactly with an existing vertex, i.e. with one of
    // the database's own operating conditions.  Return that sphere's data
    // rather than a blend of its neighbours.  Without this, hover (mu = 0)
    // sits on the convex hull boundary where locate() is ambiguous between
    // the infinite-cell branch below and the finite-cell volume ratios, and
    // Database::self_test -- which queries every sphere at its own condition
    // -- disagreed by up to 5.8 dB on all 24 hover spheres.
    // Both degenerate cases put all the weight on a single vertex and differ
    // only in which one: the query's own vertex when it lands exactly on one,
    // otherwise the nearest, when locate() walked off the hull.
    const bool on_a_vertex = locate_type == Triangulation::VERTEX;
    if (on_a_vertex || triangles.is_infinite(tetrahedron)) {
        weights[0] = 1.0;
        weights[1] = 0.0;
        weights[2] = 0.0;
        weights[3] = 0.0;
        indicies[0] = on_a_vertex ? tetrahedron->vertex(li)->info()
                                  : triangles.nearest_vertex(request)->info();
        indicies[1] = 0;
        indicies[2] = 0;
        indicies[3] = 0;
    } else {
        // Compute the barycentric weights from surrounding verticies
        auto x = std::array<Point, 4>{};
        for (size_t i = 0; i < 4; i++) {
            x[i] = tetrahedron->vertex(i)->point();
        }
        auto volume = std::array<double, 4>{};
        // NB: double, not float -- a float accumulator here capped every
        // barycentric weight at ~7 significant digits.
        auto total_volume = 0.0;
        for (size_t i = 0; i < 4; i++) {
            auto facet = triangles.triangle(Triangulation::Facet(tetrahedron, i));
            auto pr = std::array<std::array<double, 3>, 3>{};
            for (size_t j = 0; j < 3; j++) {
                pr[j][0] = facet[j].x() - request.x();
                pr[j][1] = facet[j].y() - request.y();
                pr[j][2] = facet[j].z() - request.z();
            }
            // Volume of tetra formed by point and facet
            volume[i] = std::fabs(1. / 6. * dot(pr[0], cross(pr[1], pr[2])));
            // Accumulate volume
            total_volume += volume[i];
        }
        for (size_t i = 0; i < 4; i++) {
            //Calculate barycentric weights and add indicies
            weights[i] = volume[i] / total_volume;
            indicies[i] = tetrahedron->vertex(i)->info();
        }
    }
}

double bilinear_interp(double q11, double q12, double q21, double q22,
                       double x1, double x2, double y1, double y2, double x, double y) {
    auto x2x1 = x2 - x1;
    auto y2y1 = y2 - y1;
    auto x2x = x2 - x;
    auto y2y = y2 - y;
    auto yy1 = y - y1;
    auto xx1 = x - x1;
    assert(x2x1 >= 0.0 && y2y1 >= 0.0);

    if (x2x1 < 10 * std::numeric_limits<double>::epsilon()) {
        x2x1 = 2.0;
        x2x = 1.0;
        xx1 = 1.0;
    } else if (y2y1 < 10 * std::numeric_limits<double>::epsilon()) {
        y2y1 = 2.0;
        y2y = 1.0;
        yy1 = 1.0;
    }
    return (q11 * x2x * y2y +
            q21 * xx1 * y2y +
            q12 * x2x * yy1 +
            q22 * xx1 * yy1)
           / (x2x1 * y2y1);
}
}