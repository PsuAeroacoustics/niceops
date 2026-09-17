#pragma once

// RAII wrapper over the extern "C" splinter B-spline API (splinter.cpp),
// replacing D's splines.d.  D's ~this() destructor becomes natural RAII here.

#include <cstddef>
#include <vector>

#include "util/logger.hpp"

extern "C"
{
std::size_t create_spline(std::size_t nsamples, const double* x_values,
        const double* y_values, int degree, double alpha);

void destroy_spline(std::size_t key);

double eval_spline(std::size_t key, double xi);

double der_spline(std::size_t key, double xi);

double der2_spline(std::size_t key, double xi);
}

namespace niceops {

class Spline
{
public:
    Spline(const std::vector<double>& x, const std::vector<double>& f, int degree = 3,
            double alpha = 0.3)
        : spline_key(create_spline(x.size(), x.data(), f.data(), degree, alpha))
    {
    }

    ~Spline()
    {
        destroy_spline(spline_key);
    }

    Spline(const Spline&) = delete;
    Spline& operator=(const Spline&) = delete;

    double eval(double xi, int order = 0) const
    {
        switch (order)
        {
        case 0:
            return eval_spline(spline_key, xi);
        case 1:
            return der_spline(spline_key, xi);
        case 2:
            return der2_spline(spline_key, xi);
        default:
            fatal("Invalid interpolation order = ", order);
        }
    }

    // Convenience function to apply eval to entire array
    std::vector<double> eval(const std::vector<double>& xi, int order = 0) const
    {
        std::vector<double> out(xi.size());
        for (std::size_t i = 0; i < xi.size(); i++)
            out[i] = eval(xi[i], order);
        return out;
    }

private:
    std::size_t spline_key;
};

} // namespace niceops
