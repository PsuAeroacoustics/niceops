#include "interpolator.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <vector>

#include "delaunay.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// The triangulation, and the map from its vertices back to database spheres.
// A vertex is not always a sphere: the domain extension below adds duplicates
// that carry their original's index.
std::unique_ptr<Delaunay> mesh;
std::vector<std::size_t> vertex_to_point;
// Vertex coordinates are kept alongside, so the nearest-condition fallback
// works even when there is no triangulation to fall back from.
std::vector<double> vertex_coordinates;
std::size_t dimension = 0;

std::size_t nearest(const double* query, std::size_t dim)
{
    std::size_t best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    const std::size_t count = dim ? vertex_coordinates.size() / dim : 0;
    for (std::size_t i = 0; i < count; i++)
    {
        double distance = 0.0;
        for (std::size_t j = 0; j < dim; j++)
        {
            const double delta = vertex_coordinates[i * dim + j] - query[j];
            distance += delta * delta;
        }
        if (distance < best_distance)
        {
            best_distance = distance;
            best = i;
        }
    }
    return best;
}

bool approximately(double x, double y)
{
    return std::fabs(x - y) < 10.0 * std::numeric_limits<double>::epsilon();
}

} // namespace

namespace niceops {

int set_interpolator(std::size_t dim, std::size_t count, const double* points,
        double extend_domain_by)
{
    if (dim == 0 || count == 0)
        return 1;

    const double high = 1.0 + extend_domain_by;
    const double low = 0.0 - extend_domain_by;

    std::vector<double> vertices;
    vertices.reserve(count * dim * 2);
    vertex_to_point.clear();
    vertex_to_point.reserve(count * 2);

    for (std::size_t i = 0; i < count; i++)
    {
        const double* point = points + i * dim;
        vertices.insert(vertices.end(), point, point + dim);
        vertex_to_point.push_back(i);

        // A point on any face of the unit cube gets a copy pushed out to the
        // widened bound, so the edge of the operating envelope stays inside the
        // hull instead of falling off it.
        bool on_boundary = false;
        for (std::size_t j = 0; j < dim; j++)
        {
            if (approximately(point[j], 0.0) || approximately(point[j], 1.0))
            {
                on_boundary = true;
                break;
            }
        }
        if (on_boundary)
        {
            for (std::size_t j = 0; j < dim; j++)
            {
                double value = point[j];
                if (approximately(value, 1.0))
                    value = high;
                else if (approximately(value, 0.0))
                    value = low;
                vertices.push_back(value);
            }
            vertex_to_point.push_back(i);
        }
    }

    dimension = dim;
    vertex_coordinates = vertices;

    // A database can hold too few conditions to triangulate at all -- a single
    // exported sphere is one point.  That is not an error: every query simply
    // takes the nearest condition, which for one sphere is the only answer
    // there is.
    if (vertex_to_point.size() < dim + 1)
    {
        mesh.reset();
        info("Only ", count, " operating condition", count == 1 ? "" : "s",
                "; every query will take the nearest rather than a blend.");
        return 0;
    }

    mesh = std::make_unique<Delaunay>(dim, vertex_to_point.size(), vertices.data());
    info("Triangulated ", count, " operating conditions in ", dim, " dimensions (",
            vertex_to_point.size(), " vertices after extending the domain, ",
            mesh->simplex_count(), " simplices).");
    return 0;
}

void barycentric_weights(const double* query, double* weights, std::size_t* indices)
{
    const std::size_t dim = dimension;
    std::vector<double> clamped(dim);
    for (std::size_t j = 0; j < dim; j++)
        clamped[j] = std::clamp(query[j], 0.0, 1.0);

    for (std::size_t v = 0; v <= dim; v++)
    {
        weights[v] = 0.0;
        indices[v] = 0;
    }

    std::vector<double> raw(dim + 1, 0.0);
    std::vector<std::size_t> vertices(dim + 1, 0);
    if (mesh && mesh->barycentric(clamped.data(), raw.data(), vertices.data()))
    {
        // A query sitting exactly on one of the database's own conditions takes
        // that sphere's data whole rather than a blend of its neighbours.  Hover
        // at mu = 0 lands here, on the boundary of the envelope, and blending it
        // with anything disagreed with the stored levels by up to 5.8 dB.
        for (std::size_t v = 0; v <= dim; v++)
        {
            if (raw[v] >= 1.0)
            {
                weights[0] = 1.0;
                indices[0] = vertex_to_point[vertices[v]];
                return;
            }
        }
        for (std::size_t v = 0; v <= dim; v++)
        {
            weights[v] = raw[v];
            indices[v] = vertex_to_point[vertices[v]];
        }
        return;
    }

    // Outside the hull even after the domain extension: fall back to the
    // nearest condition rather than extrapolating past the envelope.
    weights[0] = 1.0;
    indices[0] = vertex_to_point.empty() ? 0 : vertex_to_point[nearest(clamped.data(), dim)];
}

} // namespace niceops
