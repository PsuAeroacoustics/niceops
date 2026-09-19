#include "weighting.hpp"

#include <cmath>
#include <limits>

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace niceops {

double a_weighting(double frequency)
{
    // panam returns -300 rather than -inf at and below zero, so that a weighted
    // sum stays finite; keep that.
    if (frequency <= 0.0)
        return -300.0;

    const double f2 = frequency * frequency;
    const double f4 = f2 * f2;

    // Two pole pairs and a double pole pair, written as panam writes them so
    // the constants can be compared side by side.
    const double first = 1.562339 * f4
            / ((f2 + 107.65265 * 107.65265) * (f2 + 737.86223 * 737.86223));
    const double second = 2.242881e16 * f4
            / ((f2 + 20.598997 * 20.598997) * (f2 + 20.598997 * 20.598997)
                    * (f2 + 12194.22 * 12194.22) * (f2 + 12194.22 * 12194.22));

    return 10.0 * std::log10(first) + 10.0 * std::log10(second);
}

std::vector<double> a_weighting_gains(const std::vector<double>& frequencies, BandLimits limits)
{
    std::vector<double> gains;
    gains.reserve(frequencies.size());
    for (const double frequency : frequencies)
    {
        // Exactly zero outside the limits, so an excluded band contributes
        // nothing at all rather than something very small.
        gains.push_back(limits.contains(frequency) ? std::pow(10.0, a_weighting(frequency) / 10.0)
                                                   : 0.0);
    }
    return gains;
}

double weighted_level(const double* band_levels, const std::vector<double>& gains)
{
    // exp(x * ln10/10) rather than pow(10, x/10).  They are the same number to
    // within rounding, but pow() is a general two-argument routine that has to
    // handle a variable base; this runs 26 times per evaluation and ~94M
    // evaluations per run, so the difference is not academic.
    constexpr double ln10_over_10 = 0.23025850929940456840; // ln(10)/10
    double total = 0.0;
    for (std::size_t band = 0; band < gains.size(); band++)
    {
        if (gains[band] == 0.0)
            continue;
        total += gains[band] * std::exp(band_levels[band] * ln10_over_10);
    }
    if (total <= 0.0)
        return -std::numeric_limits<double>::infinity();
    return 10.0 * std::log10(total);
}

std::size_t bands_in_range(const std::vector<double>& frequencies, BandLimits limits)
{
    std::size_t count = 0;
    for (const double frequency : frequencies)
    {
        if (limits.contains(frequency))
            count++;
    }
    return count;
}

} // namespace niceops
