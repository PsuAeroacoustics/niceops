#include "delaunay.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "predicates.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace niceops {

Delaunay::Delaunay(std::size_t dim, std::size_t count, const double* points)
    : dim_(dim)
    , count_(count)
{
    if (dim_ == 0)
        fatal("Delaunay: dimension must be at least 1.");
    if (count_ < dim_ + 1)
        fatal("Delaunay: need at least ", dim_ + 1, " points in ", dim_, " dimensions, got ", count_, ".");

    storage_.assign(points, points + count_ * dim_);
    build_enclosing_simplex();
    for (std::size_t i = 0; i < count_; i++)
        insert(i);
}

void Delaunay::build_enclosing_simplex()
{
    // A simplex large enough to contain every input point, so that each
    // insertion has somewhere to land.  Its vertices are artificial and any
    // simplex still touching one at the end lies outside the real hull.
    std::vector<double> low(dim_, std::numeric_limits<double>::infinity());
    std::vector<double> high(dim_, -std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < count_; i++)
    {
        for (std::size_t j = 0; j < dim_; j++)
        {
            low[j] = std::min(low[j], storage_[i * dim_ + j]);
            high[j] = std::max(high[j], storage_[i * dim_ + j]);
        }
    }
    double extent = 1.0;
    for (std::size_t j = 0; j < dim_; j++)
        extent = std::max(extent, high[j] - low[j]);
    // Generous: the enclosing simplex only has to be big, and being far away
    // keeps its vertices from influencing the shape of the real triangulation.
    const double scale = 1000.0 * extent;

    // The corner simplex {x : x_j >= v0_j, sum(x_j - v0_j) <= (dim+1) * scale},
    // whose interior comfortably contains the bounding box.
    std::vector<double> origin(dim_);
    for (std::size_t j = 0; j < dim_; j++)
        origin[j] = low[j] - scale;

    for (std::size_t v = 0; v <= dim_; v++)
    {
        for (std::size_t j = 0; j < dim_; j++)
        {
            double value = origin[j];
            if (v > 0 && j == v - 1)
                value += static_cast<double>(dim_ + 1) * scale * 2.0;
            storage_.push_back(value);
        }
    }

    Simplex root;
    root.vertices.resize(dim_ + 1);
    root.neighbours.assign(dim_ + 1, -1);
    for (std::size_t v = 0; v <= dim_; v++)
        root.vertices[v] = count_ + v;
    // Store every simplex positively oriented, so the walk below can test a
    // single sign rather than comparing against each simplex's own.
    if (simplex_orientation(root) < 0)
        std::swap(root.vertices[0], root.vertices[1]);
    simplices_.push_back(std::move(root));
    hint_ = 0;
}

int Delaunay::simplex_orientation(const Simplex& simplex) const
{
    std::vector<double> points((dim_ + 1) * dim_);
    for (std::size_t v = 0; v <= dim_; v++)
    {
        const double* source = coordinates(simplex.vertices[v]);
        std::copy(source, source + dim_, points.begin() + v * dim_);
    }
    return orientation(points.data(), dim_);
}

int Delaunay::orientation_against(const Simplex& simplex, std::size_t replaced,
        const double* query) const
{
    std::vector<double> points((dim_ + 1) * dim_);
    for (std::size_t v = 0; v <= dim_; v++)
    {
        const double* source = (v == replaced) ? query : coordinates(simplex.vertices[v]);
        std::copy(source, source + dim_, points.begin() + v * dim_);
    }
    return orientation(points.data(), dim_);
}

long Delaunay::locate(const double* query) const
{
    // Visibility walk: step across whichever facet the query lies beyond, until
    // no facet has it outside.  Simplices are positively oriented, so "beyond
    // facet i" is a negative orientation with vertex i swapped for the query.
    long current = static_cast<long>(hint_);
    if (current < 0 || static_cast<std::size_t>(current) >= simplices_.size()
            || !simplices_[static_cast<std::size_t>(current)].alive)
    {
        current = -1;
        for (std::size_t i = 0; i < simplices_.size(); i++)
        {
            if (simplices_[i].alive)
            {
                current = static_cast<long>(i);
                break;
            }
        }
        if (current < 0)
            return -1;
    }

    // The bound is a safety net against a cycle; a visibility walk on a valid
    // triangulation terminates, but a bug here should not hang the program.
    const std::size_t limit = 8 * simplices_.size() + 64;
    for (std::size_t step = 0; step < limit; step++)
    {
        const Simplex& simplex = simplices_[static_cast<std::size_t>(current)];
        bool moved = false;
        for (std::size_t facet = 0; facet <= dim_; facet++)
        {
            if (orientation_against(simplex, facet, query) < 0)
            {
                const long next = simplex.neighbours[facet];
                if (next < 0)
                    return -1; // walked off the triangulation entirely
                current = next;
                moved = true;
                break;
            }
        }
        if (!moved)
        {
            hint_ = static_cast<std::size_t>(current);
            return current;
        }
    }
    return -1;
}

void Delaunay::insert(std::size_t point)
{
    const double* query = coordinates(point);
    const long start = locate(query);
    if (start < 0)
        fatal("Delaunay: could not locate a simplex for point ", point, ".");

    // Grow the cavity: every simplex whose circumsphere strictly contains the
    // new point must go.  A simplex whose circumsphere passes exactly through it
    // is left alone -- a point on a circumsphere does not violate the Delaunay
    // condition, and leaving it is what keeps the choice consistent when whole
    // families of points are cospherical.
    std::vector<std::size_t> cavity;
    std::vector<char> in_cavity(simplices_.size(), 0);
    std::vector<std::size_t> stack {static_cast<std::size_t>(start)};
    in_cavity[static_cast<std::size_t>(start)] = 1;
    cavity.push_back(static_cast<std::size_t>(start));

    std::vector<double> facet_points((dim_ + 1) * dim_);
    while (!stack.empty())
    {
        const std::size_t current = stack.back();
        stack.pop_back();
        for (std::size_t facet = 0; facet <= dim_; facet++)
        {
            const long neighbour = simplices_[current].neighbours[facet];
            if (neighbour < 0)
                continue;
            const auto index = static_cast<std::size_t>(neighbour);
            if (in_cavity[index])
                continue;
            const Simplex& candidate = simplices_[index];
            for (std::size_t v = 0; v <= dim_; v++)
            {
                const double* source = coordinates(candidate.vertices[v]);
                std::copy(source, source + dim_, facet_points.begin() + v * dim_);
            }
            if (in_sphere(facet_points.data(), query, dim_) > 0)
            {
                in_cavity[index] = 1;
                cavity.push_back(index);
                stack.push_back(index);
            }
        }
    }

    // Every facet of the cavity that faces outward becomes the base of a new
    // simplex joining it to the new point.
    struct Boundary
    {
        std::vector<std::size_t> vertices; // the facet, dim vertices
        long outside;                      // simplex on the far side, or -1
        std::size_t outside_facet;         // which of its facets looks back
    };
    std::vector<Boundary> boundary;
    for (std::size_t index : cavity)
    {
        const Simplex& simplex = simplices_[index];
        for (std::size_t facet = 0; facet <= dim_; facet++)
        {
            const long neighbour = simplex.neighbours[facet];
            if (neighbour >= 0 && in_cavity[static_cast<std::size_t>(neighbour)])
                continue;
            Boundary entry;
            entry.outside = neighbour;
            entry.outside_facet = 0;
            for (std::size_t v = 0; v <= dim_; v++)
            {
                if (v != facet)
                    entry.vertices.push_back(simplex.vertices[v]);
            }
            if (neighbour >= 0)
            {
                const Simplex& outside = simplices_[static_cast<std::size_t>(neighbour)];
                for (std::size_t v = 0; v <= dim_; v++)
                {
                    if (outside.neighbours[v] == static_cast<long>(index))
                    {
                        entry.outside_facet = v;
                        break;
                    }
                }
            }
            boundary.push_back(std::move(entry));
        }
    }

    for (std::size_t index : cavity)
    {
        simplices_[index].alive = false;
        free_list_.push_back(index);
    }

    // Create the new simplices, then stitch them to the cavity wall and to one
    // another.  Two new simplices are adjacent where they share the facet made
    // of the new point and dim-1 of a boundary facet's vertices; keying on the
    // sorted vertex set finds those pairs without any geometry.
    std::vector<std::size_t> created;
    created.reserve(boundary.size());
    for (const Boundary& entry : boundary)
    {
        Simplex simplex;
        simplex.vertices = entry.vertices;
        simplex.vertices.push_back(point);
        simplex.neighbours.assign(dim_ + 1, -1);
        simplex.alive = true;
        if (simplex_orientation(simplex) < 0)
            std::swap(simplex.vertices[0], simplex.vertices[1]);

        std::size_t slot = 0;
        if (!free_list_.empty())
        {
            slot = free_list_.back();
            free_list_.pop_back();
            simplices_[slot] = std::move(simplex);
        }
        else
        {
            slot = simplices_.size();
            simplices_.push_back(std::move(simplex));
        }
        created.push_back(slot);

        // The facet opposite the new point is the boundary facet; link it back.
        const Simplex& placed = simplices_[slot];
        std::size_t opposite_new_point = 0;
        for (std::size_t v = 0; v <= dim_; v++)
        {
            if (placed.vertices[v] == point)
            {
                opposite_new_point = v;
                break;
            }
        }
        simplices_[slot].neighbours[opposite_new_point] = entry.outside;
        if (entry.outside >= 0)
            simplices_[static_cast<std::size_t>(entry.outside)].neighbours[entry.outside_facet]
                    = static_cast<long>(slot);
    }

    std::map<std::vector<std::size_t>, std::pair<std::size_t, std::size_t>> open_facets;
    for (std::size_t slot : created)
    {
        const Simplex& simplex = simplices_[slot];
        for (std::size_t facet = 0; facet <= dim_; facet++)
        {
            if (simplex.vertices[facet] == point)
                continue; // that one is the cavity wall, already linked
            std::vector<std::size_t> key;
            for (std::size_t v = 0; v <= dim_; v++)
            {
                if (v != facet)
                    key.push_back(simplex.vertices[v]);
            }
            std::sort(key.begin(), key.end());
            auto found = open_facets.find(key);
            if (found == open_facets.end())
            {
                open_facets.emplace(std::move(key), std::make_pair(slot, facet));
            }
            else
            {
                const auto other = found->second;
                simplices_[slot].neighbours[facet] = static_cast<long>(other.first);
                simplices_[other.first].neighbours[other.second] = static_cast<long>(slot);
                open_facets.erase(found);
            }
        }
    }

    hint_ = created.empty() ? hint_ : created.front();
}

bool Delaunay::barycentric(const double* query, double* weights, std::size_t* indices) const
{
    const long found = locate(query);
    if (found < 0)
        return false;
    const Simplex& simplex = simplices_[static_cast<std::size_t>(found)];
    for (std::size_t v = 0; v <= dim_; v++)
    {
        if (is_artificial(simplex.vertices[v]))
            return false; // outside the hull of the real points
    }

    // Weight of a vertex is the volume of the simplex formed by the query and
    // the facet opposite it, over the total.  Same formulation the CGAL-based
    // interpolator used, so the numbers land in the same place.
    std::vector<double> volumes(dim_ + 1, 0.0);
    std::vector<double> facet(dim_ * dim_);
    double total = 0.0;
    for (std::size_t v = 0; v <= dim_; v++)
    {
        std::size_t row = 0;
        for (std::size_t other = 0; other <= dim_; other++)
        {
            if (other == v)
                continue;
            const double* source = coordinates(simplex.vertices[other]);
            for (std::size_t j = 0; j < dim_; j++)
                facet[row * dim_ + j] = source[j] - query[j];
            row++;
        }
        // |det| of the edge vectors is the simplex volume up to a constant
        // factor, and the constant cancels in the ratio below.
        double determinant = 1.0;
        std::vector<double> matrix(facet);
        for (std::size_t column = 0; column < dim_; column++)
        {
            std::size_t pivot_row = column;
            for (std::size_t r = column + 1; r < dim_; r++)
            {
                if (std::fabs(matrix[r * dim_ + column])
                        > std::fabs(matrix[pivot_row * dim_ + column]))
                    pivot_row = r;
            }
            const double pivot = matrix[pivot_row * dim_ + column];
            if (pivot == 0.0)
            {
                determinant = 0.0;
                break;
            }
            if (pivot_row != column)
            {
                for (std::size_t j = 0; j < dim_; j++)
                    std::swap(matrix[column * dim_ + j], matrix[pivot_row * dim_ + j]);
                determinant = -determinant;
            }
            determinant *= pivot;
            for (std::size_t r = column + 1; r < dim_; r++)
            {
                const double factor = matrix[r * dim_ + column] / pivot;
                if (factor == 0.0)
                    continue;
                for (std::size_t j = column; j < dim_; j++)
                    matrix[r * dim_ + j] -= factor * matrix[column * dim_ + j];
            }
        }
        volumes[v] = std::fabs(determinant);
        total += volumes[v];
    }

    if (!(total > 0.0))
    {
        // Degenerate cell: put all the weight on the first vertex rather than
        // divide by zero.
        for (std::size_t v = 0; v <= dim_; v++)
        {
            weights[v] = v == 0 ? 1.0 : 0.0;
            indices[v] = simplex.vertices[v];
        }
        return true;
    }
    for (std::size_t v = 0; v <= dim_; v++)
    {
        weights[v] = volumes[v] / total;
        indices[v] = simplex.vertices[v];
    }
    return true;
}

std::size_t Delaunay::nearest_vertex(const double* query) const
{
    std::size_t best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < count_; i++)
    {
        const double* point = coordinates(i);
        double distance = 0.0;
        for (std::size_t j = 0; j < dim_; j++)
        {
            const double delta = point[j] - query[j];
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

std::size_t Delaunay::simplex_count() const
{
    std::size_t total = 0;
    for (const Simplex& simplex : simplices_)
    {
        if (!simplex.alive)
            continue;
        bool real = true;
        for (std::size_t vertex : simplex.vertices)
        {
            if (is_artificial(vertex))
            {
                real = false;
                break;
            }
        }
        if (real)
            total++;
    }
    return total;
}

} // namespace niceops
