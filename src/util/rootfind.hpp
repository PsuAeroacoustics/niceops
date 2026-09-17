#pragma once

// Bisection root finder replacing D's std.numeric.findRoot.  D's findRoot
// stops when the tolerance function is satisfied on the bracket endpoints
// (here: |hi - lo| <= tol), which is what this implements.

#include <cmath>
#include <functional>

namespace niceops {

inline double bisect(const std::function<double(double)>& f, double lo, double hi, double tol)
{
    double flo = f(lo);
    while (std::fabs(hi - lo) > tol)
    {
        double mid = 0.5 * (lo + hi);
        double fmid = f(mid);
        if (std::fabs(fmid) < tol)
            return mid;
        if ((flo < 0.0) == (fmid < 0.0))
        {
            lo = mid;
            flo = fmid;
        }
        else
        {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

} // namespace niceops
