#pragma once

// D's iota(begin, end, step): evenly spaced values from `begin`, EXCLUSIVE of
// `end`.  Counted by index rather than accumulated, so the spacing cannot drift.
//
// The number of values is the ceiling of the span in steps, less a sliver for
// round-off.  An exact multiple of the step is the case the exclusive rule is
// about (-1500..1500 by 100 is 30 values, not 31), and it is the case round-off
// breaks: run() converts each bound and the step from feet to meters on its own,
// so the quotient can come out one ULP above the integer and a bare ceil() then
// adds a row and a column that the metric run of the same grid does not have.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace niceops {

inline std::size_t iota_exclusive_count(double begin, double end, double step)
{
    const double steps = (end - begin) / step;
    if (!(steps > 0.0))
        return 0;
    const double count = std::ceil(steps - 1e-9 * std::max(1.0, steps));
    return count > 0.0 ? static_cast<std::size_t>(count) : 0;
}

inline std::vector<double> iota_exclusive(double begin, double end, double step)
{
    std::vector<double> values(iota_exclusive_count(begin, end, step));
    for (std::size_t i = 0; i < values.size(); i++)
        values[i] = begin + step * static_cast<double>(i);
    return values;
}

} // namespace niceops
