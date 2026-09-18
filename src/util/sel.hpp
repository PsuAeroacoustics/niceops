#pragma once

// Sound exposure level integration from niceops.d.  peak_range() finds the
// 10-dB-down-from-peak window; sound_exposure_level() integrates that window
// into a SEL: 10*log10(T * sum(10^(0.1*L))).  Kept in a header so the unit
// tests can reach them.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace niceops {

inline std::array<std::size_t, 2> peak_range(const std::vector<double>& levels, double down = 10.0)
{
    // Find the 10-dB down points
    std::array<std::size_t, 2> index = {0, levels.size() - 1};
    auto peak_position = static_cast<std::size_t>(std::max_element(levels.begin(), levels.end())
            - levels.begin());
    const double threshold = levels[peak_position] - down;

    for (std::size_t i = peak_position; i < levels.size(); i++)
    {
        if (levels[i] < threshold)
        {
            index[1] = i - 1;
            break;
        }
    }
    for (long i = static_cast<long>(peak_position); i > 0; i--)
    {
        if (levels[i] <= threshold)
        {
            if (levels[i] == threshold)
            {
                index[0] = static_cast<std::size_t>(i);
            }
            else
            {
                index[0] = static_cast<std::size_t>(i) + 1;
            }
            break;
        }
    }
    return index;
}

// Mean spacing of a monotonically increasing time series, and how far the worst
// interval strays from it (relative).  The SEL integral weights every sample by
// this, so it must come from the data rather than from an assumed step.
inline double uniform_sample_interval(const std::vector<double>& times, double* relative_spread = nullptr)
{
    if (times.size() < 2)
    {
        if (relative_spread != nullptr)
            *relative_spread = 0.0;
        return 0.0;
    }
    const double mean = (times.back() - times.front()) / static_cast<double>(times.size() - 1);
    if (relative_spread != nullptr)
    {
        double worst = 0.0;
        if (mean > 0.0)
        {
            for (std::size_t i = 1; i < times.size(); i++)
                worst = std::max(worst, std::fabs((times[i] - times[i - 1]) - mean) / mean);
        }
        *relative_spread = worst;
    }
    return mean;
}

inline double sound_exposure_level(const std::vector<double>& levels, double T)
{
    auto index = peak_range(levels);
    double sum = 0.0;
    for (std::size_t i = index[0]; i <= index[1]; i++)
    {
        sum += std::pow(10.0, 0.1 * levels[i]);
    }
    return 10.0 * std::log10(T * sum);
}

} // namespace niceops
