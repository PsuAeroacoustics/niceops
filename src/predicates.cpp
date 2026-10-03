#include "predicates.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "util/logger.hpp"

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

// Upper bound on the matrix size these predicates are called with: an in-sphere
// test in `dim` dimensions builds a (dim + 1) x (dim + 1) matrix, so this caps
// the dimension at MAX_PREDICATE_SIZE - 1.  Checked at the entry points.
constexpr std::size_t MAX_PREDICATE_SIZE = 8;

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
    // Two components per product pair is the bound the expansion can reach
    // before cancellation shortens it, so this reserve is the last allocation
    // the loop needs rather than one per growth step.
    result.reserve(2 * a.size() * b.size());
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
//
// The matrix arrives as an array of pointers and the minors are built the same
// way, rather than as copies.  Each entry is an expansion -- a heap-allocated
// vector -- and copying (size-1)^2 of them at every level of a factorial
// recursion put malloc and free at the top of the database load profile; with
// pointers the recursion allocates nothing but the expansions it actually
// computes.  `minor` is a stack array for the same reason `laplace` uses one.
//
// Callers go through the std::vector overload below, which adapts to this once;
// the pointer form is what the recursion passes to itself.
Expansion exact_determinant_recurse(const Expansion* const* matrix, std::size_t size)
{
    assert(size <= MAX_PREDICATE_SIZE);
    if (size == 1)
        return *matrix[0];

    Expansion total;
    const Expansion* minor[MAX_PREDICATE_SIZE * MAX_PREDICATE_SIZE];
    for (std::size_t column = 0; column < size; column++)
    {
        const Expansion& pivot = *matrix[column];
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
        Expansion term = multiply(pivot, exact_determinant_recurse(minor, size - 1));
        if (column % 2 == 1)
        {
            for (double& component : term)
                component = -component;
        }
        add_expansion(total, term);
    }
    return total;
}

// Exact determinant of a matrix held as expansions, which is how all three
// predicates build theirs.  Adapting to the pointer form once, here, keeps that
// four-line dance out of every call site while leaving the recursion itself
// allocation-free.
Expansion exact_determinant(const std::vector<Expansion>& matrix, std::size_t size)
{
    assert(size <= MAX_PREDICATE_SIZE);
    const Expansion* rows[MAX_PREDICATE_SIZE * MAX_PREDICATE_SIZE];
    for (std::size_t i = 0; i < size * size; i++)
        rows[i] = &matrix[i];
    return exact_determinant_recurse(rows, size);
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
    // Stack buffer, not a vector: this recursion is factorial in `size`, so an
    // allocation per level put malloc at about a fifth of the whole database
    // load.  MAX_PREDICATE_SIZE bounds the dimensions this is ever called with.
    double minor[MAX_PREDICATE_SIZE * MAX_PREDICATE_SIZE];
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
        laplace(minor, size - 1, sub_determinant, sub_permanent);
        const double value = matrix[column] * sub_determinant;
        determinant += (column % 2 == 0) ? value : -value;
        permanent += std::fabs(matrix[column]) * sub_permanent;
    }
}

// Does the filtered double evaluation settle the sign on its own?  Returns 2
// when it cannot tell, which is not a sign any caller can mistake.
//
// The permanent bounds the round-off: a Laplace expansion of this size cannot
// accumulate more than a small multiple of eps times it.  The bound also has to
// cover the error already in the entries, since the callers form them by
// subtracting in plain double, and the determinant is multilinear so a relative
// eps on every entry moves it by about size * eps * permanent.  The constant
// below is comfortably above both terms; being generous only costs an occasional
// unnecessary trip through the exact path.
int filtered_sign(const double* approximate, std::size_t size)
{
    double determinant = 0.0;
    double permanent = 0.0;
    laplace(approximate, size, determinant, permanent);
    const double epsilon = 2.220446049250313e-16;
    const double bound = static_cast<double>(4 * size + 8) * epsilon * permanent;
    if (std::fabs(determinant) > bound)
        return determinant > 0.0 ? 1 : -1;
    return 2;
}

} // namespace

int exact_determinant_sign(const double* matrix, std::size_t size)
{
    if (size == 0)
        return 0;
    if (size > MAX_PREDICATE_SIZE)
        fatal("Predicates: matrix of size ", size, " exceeds MAX_PREDICATE_SIZE.");
    const int filtered = filtered_sign(matrix, size);
    if (filtered != 2)
        return filtered;

    std::vector<Expansion> exact(size * size);
    for (std::size_t i = 0; i < size * size; i++)
    {
        if (matrix[i] != 0.0)
            exact[i] = Expansion {matrix[i]};
    }
    return sign_of(exact_determinant(exact, size));
}

int orientation(const double* points, std::size_t dim)
{
    if (dim > MAX_PREDICATE_SIZE)
        fatal("Predicates: dimension ", dim, " exceeds MAX_PREDICATE_SIZE.");
    // Rows are the edge vectors from the last point.  The cheap path forms them
    // in plain double; only if that cannot settle the sign are they re-formed
    // exactly, because taking the difference in double rounds away precisely the
    // information the exact determinant exists to use.
    const double* base = points + dim * dim;
    thread_local std::vector<double> approximate;
    approximate.resize(dim * dim);
    for (std::size_t row = 0; row < dim; row++)
    {
        for (std::size_t column = 0; column < dim; column++)
            approximate[row * dim + column] = points[row * dim + column] - base[column];
    }
    const int filtered = filtered_sign(approximate.data(), dim);
    if (filtered != 2)
        return filtered;

    std::vector<Expansion> exact(dim * dim);
    for (std::size_t row = 0; row < dim; row++)
    {
        for (std::size_t column = 0; column < dim; column++)
            exact[row * dim + column] = exact_difference(points[row * dim + column], base[column]);
    }
    return sign_of(exact_determinant(exact, dim));
}

int in_sphere(const double* points, const double* query, std::size_t dim)
{
    // Lift to the paraboloid about the query point: the extra column is the
    // squared distance, so "inside the circumsphere" becomes an orientation test
    // one dimension up.
    const std::size_t size = dim + 1;
    if (size > MAX_PREDICATE_SIZE)
        fatal("Predicates: dimension ", dim, " exceeds MAX_PREDICATE_SIZE.");
    thread_local std::vector<double> approximate;
    approximate.resize(size * size);
    for (std::size_t row = 0; row < size; row++)
    {
        double lifted = 0.0;
        for (std::size_t column = 0; column < dim; column++)
        {
            const double delta = points[row * dim + column] - query[column];
            approximate[row * size + column] = delta;
            lifted += delta * delta;
        }
        approximate[row * size + dim] = lifted;
    }

    int sign = filtered_sign(approximate.data(), size);
    if (sign == 2)
    {
        // Exact path: both the differences and their squares are formed without
        // rounding.
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
        sign = sign_of(exact_determinant(exact, size));
    }

    // Normalize out the simplex's own orientation, so the answer means "inside"
    // however the vertices happen to be ordered.
    const int simplex_orientation = orientation(points, dim);
    if (simplex_orientation == 0)
        return 0;
    return simplex_orientation > 0 ? sign : -sign;
}

} // namespace niceops
