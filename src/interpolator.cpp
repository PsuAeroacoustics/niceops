#include "interpolator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "delaunay.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

bool approximately(double x, double y)
{
    return std::fabs(x - y) < 10.0 * std::numeric_limits<double>::epsilon();
}

} // namespace

namespace niceops {

// Out of line so the header can forward-declare Delaunay rather than include it.
ConditionInterpolator::ConditionInterpolator() = default;
ConditionInterpolator::~ConditionInterpolator() = default;
ConditionInterpolator::ConditionInterpolator(ConditionInterpolator&&) noexcept = default;
ConditionInterpolator& ConditionInterpolator::operator=(ConditionInterpolator&&) noexcept = default;

void ConditionInterpolator::build(std::size_t dim, std::size_t count, const double* points,
        double extend_domain_by)
{
    if (dim == 0 || count == 0)
    {
        input_error("Cannot interpolate over ", count, " operating conditions in ", dim,
                " dimensions.");
    }

    const double high = 1.0 + extend_domain_by;
    const double low = 0.0 - extend_domain_by;

    vertex_coordinates_.clear();
    vertex_coordinates_.reserve(count * dim * 2);
    vertex_to_point_.clear();
    vertex_to_point_.reserve(count * 2);

    for (std::size_t i = 0; i < count; i++)
    {
        const double* point = points + i * dim;
        vertex_coordinates_.insert(vertex_coordinates_.end(), point, point + dim);
        vertex_to_point_.push_back(i);

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
                vertex_coordinates_.push_back(value);
            }
            vertex_to_point_.push_back(i);
        }
    }

    dimension_ = dim;

    // A database can hold too few conditions to triangulate at all -- a single
    // exported sphere is one point.  That is not an error: every query simply
    // takes the nearest condition, which for one sphere is the only answer
    // there is.
    if (vertex_to_point_.size() < dim + 1)
    {
        mesh_.reset();
        info("Only ", count, " operating condition", count == 1 ? "" : "s",
                "; every query will take the nearest rather than a blend.");
        return;
    }

    mesh_ = std::make_unique<Delaunay>(dim, vertex_to_point_.size(), vertex_coordinates_.data());
    info("Triangulated ", count, " operating conditions in ", dim, " dimensions (",
            vertex_to_point_.size(), " vertices after extending the domain, ",
            mesh_->simplex_count(), " simplices).");
}

void ConditionInterpolator::barycentric(const double* query, double* weights,
        std::size_t* indices) const
{
    const std::size_t dim = dimension_;
    // Scratch reused across calls; this runs once per trajectory state.
    thread_local std::vector<double> clamped;
    thread_local std::vector<double> raw;
    thread_local std::vector<std::size_t> vertices;
    clamped.resize(dim);
    raw.assign(dim + 1, 0.0);
    vertices.assign(dim + 1, 0);

    for (std::size_t j = 0; j < dim; j++)
        clamped[j] = std::clamp(query[j], 0.0, 1.0);

    for (std::size_t v = 0; v <= dim; v++)
    {
        weights[v] = 0.0;
        indices[v] = 0;
    }

    if (mesh_ && mesh_->barycentric(clamped.data(), raw.data(), vertices.data()))
    {
        // A query sitting exactly on one of the database's own conditions takes
        // that sphere's data whole rather than a blend of its neighbors.  Hover
        // at mu = 0 lands here, on the boundary of the envelope, and blending it
        // with anything disagreed with the stored levels by up to 5.8 dB.
        for (std::size_t v = 0; v <= dim; v++)
        {
            if (raw[v] >= 1.0)
            {
                weights[0] = 1.0;
                indices[0] = vertex_to_point_[vertices[v]];
                return;
            }
        }
        for (std::size_t v = 0; v <= dim; v++)
        {
            weights[v] = raw[v];
            indices[v] = vertex_to_point_[vertices[v]];
        }
        return;
    }

    // Outside the hull even after the domain extension: fall back to the
    // nearest condition rather than extrapolating past the envelope.
    weights[0] = 1.0;
    indices[0] = vertex_to_point_.empty() ? 0 : vertex_to_point_[nearest(clamped.data())];
}

std::size_t ConditionInterpolator::nearest(const double* query) const
{
    // Delaunay searches the same points, so defer to it when there is a mesh.
    if (mesh_)
        return mesh_->nearest_vertex(query);

    std::size_t best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    const std::size_t count = dimension_ ? vertex_coordinates_.size() / dimension_ : 0;
    for (std::size_t i = 0; i < count; i++)
    {
        double distance = 0.0;
        for (std::size_t j = 0; j < dimension_; j++)
        {
            const double delta = vertex_coordinates_[i * dimension_ + j] - query[j];
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

} // namespace niceops
