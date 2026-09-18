#include "predicates.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace niceops {

namespace {

// ---------------------------------------------------------------------------
// Shewchuk's exact expansion arithmetic.
//
// An "expansion" is a vector of doubles whose exact sum is the value being
// represented.  Every primitive below is exact: nothing is rounded away, the
// extra components simply carry the bits a single double could not hold.
// ---------------------------------------------------------------------------

using Expansion = std::vector<double>;

// a + b = x + y exactly, for any a and b.  Knuth's two-sum.
inline void two_sum(double a, double b, double& x, double& y)
{
    x = a + b;
    const double b_virtual = x - a;
    const double a_virtual = x - b_virtual;
    y = (a - a_virtual) + (b - b_virtual);
}

// a * b = x + y exactly.  The fused multiply-add gives a*b - x with a single
// rounding, which is exactly the low half.
inline void two_product(double a, double b, double& x, double& y)
{
    x = a * b;
    y = std::fma(a, b, -x);
}

void add_component(Expansion& e, double value)
{
    if (value == 0.0)
        return;
    double carry = value;
    for (double& component : e)
    {
        double sum = 0.0;
        double error = 0.0;
        two_sum(carry, component, sum, error);
        component = error;
        carry = sum;
    }
    if (carry != 0.0)
        e.push_back(carry);
    e.erase(std::remove(e.begin(), e.end(), 0.0), e.end());
}

void add_expansion(Expansion& e, const Expansion& f)
{
    for (double component : f)
        add_component(e, component);
}

// Exact difference of two doubles, as at most two components.
Expansion exact_difference(double a, double b)
{
    double high = 0.0;
    double low = 0.0;
    two_sum(a, -b, high, low);
    Expansion result;
    if (low != 0.0)
        result.push_back(low);
    if (high != 0.0)
        result.push_back(high);
    return result;
}

Expansion multiply(const Expansion& a, const Expansion& b)
{
    Expansion result;
    for (double x : a)
    {
        for (double y : b)
        {
            double high = 0.0;
            double low = 0.0;
            two_product(x, y, high, low);
            add_component(result, low);
            add_component(result, high);
        }
    }
    return result;
}

int sign_of(const Expansion& e)
{
    // Components are non-overlapping, so the sum's sign is the sign of the
    // largest one.
    double largest = 0.0;
    for (double component : e)
    {
        if (std::fabs(component) > std::fabs(largest))
            largest = component;
    }
    return largest > 0.0 ? 1 : (largest < 0.0 ? -1 : 0);
}

// ---------------------------------------------------------------------------
// Determinants
// ---------------------------------------------------------------------------

// Exact determinant of a small matrix of expansions, by Laplace expansion along
// the first row.  Factorial in `size`, which is why the filtered path below runs
// first; for the dimensions in play this is at most a few dozen terms.
Expansion exact_determinant(const std::vector<Expansion>& matrix, std::size_t size)
{
    if (size == 1)
        return matrix[0];

    Expansion total;
    std::vector<Expansion> minor((size - 1) * (size - 1));
    for (std::size_t column = 0; column < size; column++)
    {
        const Expansion& pivot = matrix[column];
        if (pivot.empty())
            continue;
        std::size_t out = 0;
        for (std::size_t row = 1; row < size; row++)
        {
            for (std::size_t c = 0; c < size; c++)
            {
                if (c != column)
                    minor[out++] = matrix[row * size + c];
            }
        }
        Expansion term = multiply(pivot, exact_determinant(minor, size - 1));
        if (column % 2 == 1)
        {
            for (double& component : term)
                component = -component;
        }
        add_expansion(total, term);
    }
    return total;
}

// Determinant and permanent (the same expansion with every sign made positive)
// of a double matrix, both by the same Laplace recursion.  The permanent bounds
// how much round-off the determinant can have accumulated, which is what makes
// the filter below rigorous rather than a guess.
void laplace(const double* matrix, std::size_t size, double& determinant, double& permanent)
{
    if (size == 1)
    {
        determinant = matrix[0];
        permanent = std::fabs(matrix[0]);
        return;
    }
    determinant = 0.0;
    permanent = 0.0;
    std::vector<double> minor((size - 1) * (size - 1));
    for (std::size_t column = 0; column < size; column++)
    {
        std::size_t out = 0;
        for (std::size_t row = 1; row < size; row++)
        {
            for (std::size_t c = 0; c < size; c++)
            {
                if (c != column)
                    minor[out++] = matrix[row * size + c];
            }
        }
        double sub_determinant = 0.0;
        double sub_permanent = 0.0;
        laplace(minor.data(), size - 1, sub_determinant, sub_permanent);
        const double value = matrix[column] * sub_determinant;
        determinant += (column % 2 == 0) ? value : -value;
        permanent += std::fabs(matrix[column]) * sub_permanent;
    }
}

// Sign of the determinant of a matrix whose entries are exact expansions.
// `approximate` holds the same matrix collapsed to doubles, for the filter.
int determinant_sign(const std::vector<Expansion>& exact, const std::vector<double>& approximate,
        std::size_t size)
{
    double determinant = 0.0;
    double permanent = 0.0;
    laplace(approximate.data(), size, determinant, permanent);

    // Round-off in a Laplace expansion of this size cannot exceed a small
    // multiple of eps times the permanent.  Outside that band the sign is
    // certain and the exact path is not needed.
    const double epsilon = 2.220446049250313e-16;
    const double bound = static_cast<double>(size + 2) * epsilon * permanent;
    if (std::fabs(determinant) > bound)
        return determinant > 0.0 ? 1 : -1;

    return sign_of(exact_determinant(exact, size));
}

// Collapse an expansion matrix to doubles for the filter.
std::vector<double> to_double(const std::vector<Expansion>& matrix)
{
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t i = 0; i < matrix.size(); i++)
    {
        double sum = 0.0;
        for (double component : matrix[i])
            sum += component;
        result[i] = sum;
    }
    return result;
}

} // namespace

int exact_determinant_sign(const double* matrix, std::size_t size)
{
    if (size == 0)
        return 0;
    std::vector<Expansion> exact(size * size);
    for (std::size_t i = 0; i < size * size; i++)
    {
        if (matrix[i] != 0.0)
            exact[i] = Expansion {matrix[i]};
    }
    std::vector<double> approximate(matrix, matrix + size * size);
    return determinant_sign(exact, approximate, size);
}

int orientation(const double* points, std::size_t dim)
{
    // Rows are the edge vectors from the last point.  The differences are formed
    // exactly -- taking them in plain double first would round away precisely the
    // information the exact determinant exists to use.
    const double* base = points + dim * dim;
    std::vector<Expansion> exact(dim * dim);
    for (std::size_t row = 0; row < dim; row++)
    {
        for (std::size_t column = 0; column < dim; column++)
            exact[row * dim + column] = exact_difference(points[row * dim + column], base[column]);
    }
    return determinant_sign(exact, to_double(exact), dim);
}

int in_sphere(const double* points, const double* query, std::size_t dim)
{
    // Lift to the paraboloid about the query point: the extra column is the
    // squared distance, so "inside the circumsphere" becomes an orientation test
    // one dimension up.  Both the differences and their squares are exact.
    const std::size_t size = dim + 1;
    std::vector<Expansion> exact(size * size);
    for (std::size_t row = 0; row < size; row++)
    {
        Expansion lifted;
        for (std::size_t column = 0; column < dim; column++)
        {
            Expansion delta = exact_difference(points[row * dim + column], query[column]);
            add_expansion(lifted, multiply(delta, delta));
            exact[row * size + column] = std::move(delta);
        }
        exact[row * size + dim] = std::move(lifted);
    }
    const int sign = determinant_sign(exact, to_double(exact), size);

    // Normalise out the simplex's own orientation, so the answer means "inside"
    // however the vertices happen to be ordered.
    const int simplex_orientation = orientation(points, dim);
    if (simplex_orientation == 0)
        return 0;
    return simplex_orientation > 0 ? sign : -sign;
}

} // namespace niceops
