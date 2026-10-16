#pragma once

namespace niceops {

// The cubic Hermite interpolant on a node interval of width h, at the
// fraction t of the way across it: y0, y1 the values at its ends and d0, d1
// their slopes per unit of the abscissa.
inline double cubic_hermite(double t, double h, double y0, double d0, double y1, double d1)
{
    const double r = 1.0 - t;
    return (1.0 + 2.0 * t) * r * r * y0 + t * r * r * h * d0 + t * t * (3.0 - 2.0 * t) * y1 - t * t * r * h * d1;
}

} // namespace niceops
