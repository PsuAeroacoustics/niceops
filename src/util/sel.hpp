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

#include "logger.hpp"

namespace niceops {

inline std::array<std::size_t, 2> peak_range(const std::vector<double>& levels, double down = 10.0)
{
    // An empty time history has no peak.  Guarding here rather than letting
    // levels.size() - 1 wrap and max_element dereference end(): that crashed
    // with a segfault well downstream of the real cause, an unparseable
    // trajectory that yielded no states.
    if (levels.empty())
        fatal("peak_range: empty level history; the trajectory produced no states.");
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

// Whether a sample falls in the night period.  DNL's night runs 22:00 to 07:00
// local, so the UTC seconds from the trajectory need the local offset applied
// before the window is tested.
inline bool is_night(double utcsec, double utc_offset_hours)
{
    double local = std::fmod(utcsec + utc_offset_hours * 3600.0, 86400.0);
    if (local < 0.0)
        local += 86400.0;
    return local < 7.0 * 3600.0 || local >= 22.0 * 3600.0;
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

// Day-night level: the sound exposure level with a 10 dB penalty applied to
// every sample taken at night, as the D original computes it.  The penalty goes
// on before the peak window is found, so a night operation's window is located
// on the penalized history exactly as SEL's is on the raw one.
//
// `scratch` holds the penalized copy; it is a parameter so a caller looping over
// observers can reuse one buffer instead of allocating per observer.
inline double day_night_level(const std::vector<double>& levels, double T,
        const std::vector<char>& night, std::vector<double>& scratch)
{
    scratch.resize(levels.size());
    for (std::size_t i = 0; i < levels.size(); i++)
        scratch[i] = night[i] ? levels[i] + 10.0 : levels[i];
    return sound_exposure_level(scratch, T);
}

} // namespace niceops
