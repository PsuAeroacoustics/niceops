#include "pnl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr std::size_t nb = pnl_band_count;

// Table A36-3: constants for the mathematically formulated noy values, as
// panam tabulates them.  SPL(a) is infinite where the table gives none (bands
// 10 to 22 have no region above SPL(a)); M(c) is then never used.
constexpr std::array<double, nb> spl_a = {91.0, 85.9, 87.3, 79.9, 79.8, 76.0, 74.0, 74.9, 94.6, inf,
        inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, 44.3, 50.7};
constexpr std::array<double, nb> spl_b = {64.0, 60, 56, 53, 51, 48, 46, 44, 42, 40, 40, 40, 40, 40,
        38, 34, 32, 30, 29, 29, 30, 31, 37, 41};
constexpr std::array<double, nb> spl_c = {52.0, 51, 49, 47, 46, 45, 43, 42, 41, 40, 40, 40, 40, 40,
        38, 34, 32, 30, 29, 29, 30, 31, 34, 37};
constexpr std::array<double, nb> spl_d = {49.0, 44, 39, 34, 30, 27, 24, 21, 18, 16, 16, 16, 16, 16,
        15, 12, 9, 5, 4, 5, 6, 10, 17, 21};
constexpr std::array<double, nb> spl_e = {55.0, 51, 46, 42, 39, 36, 33, 30, 27, 25, 25, 25, 25, 25,
        23, 21, 18, 15, 14, 14, 15, 17, 23, 29};
constexpr std::array<double, nb> m_b = {0.043478, 0.040570, 0.036831, 0.036831, 0.035336, 0.033333,
        0.033333, 0.032051, 0.030675, 0.030103, 0.030103, 0.030103, 0.030103, 0.030103, 0.030103,
        0.029960, 0.029960, 0.029960, 0.029960, 0.029960, 0.029960, 0.029960, 0.042285, 0.042285};
constexpr std::array<double, nb> m_c = {0.030103, 0.030103, 0.030103, 0.030103, 0.030103, 0.030103,
        0.030103, 0.030103, 0.030103, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.029960, 0.029960};
constexpr std::array<double, nb> m_d = {0.079520, 0.068160, 0.068160, 0.059640, 0.053013, 0.053013,
        0.053013, 0.053013, 0.053013, 0.053013, 0.053013, 0.053013, 0.053013, 0.053013, 0.059640,
        0.053013, 0.053013, 0.047712, 0.047712, 0.053013, 0.053013, 0.068160, 0.079520, 0.059640};
constexpr std::array<double, nb> m_e = {0.058098, 0.058098, 0.052288, 0.047534, 0.043573, 0.043573,
        0.040221, 0.037349, 0.034859, 0.034859, 0.034859, 0.034859, 0.034859, 0.034859, 0.034859,
        0.040221, 0.037349, 0.034859, 0.034859, 0.034859, 0.034859, 0.037349, 0.037349, 0.043573};

} // namespace

namespace niceops {

PnlSpectrum noys(const PnlSpectrum& band_levels)
{
    PnlSpectrum n {};
    for (std::size_t i = 0; i < nb; i++)
    {
        const double spl = band_levels[i];
        // The four regions in panam's order of precedence; every comparison is
        // false for NaN, which therefore gets 0 noy, as in numpy.
        if (spl >= spl_a[i])
            n[i] = std::pow(10.0, m_c[i] * (spl - spl_c[i]));
        else if (spl >= spl_b[i])
            n[i] = std::pow(10.0, m_b[i] * (spl - spl_b[i]));
        else if (spl >= spl_e[i])
            n[i] = 0.3 * std::pow(10.0, m_e[i] * (spl - spl_e[i]));
        else if (spl >= spl_d[i])
            n[i] = 0.1 * std::pow(10.0, m_d[i] * (spl - spl_d[i]));
        else
            n[i] = 0.0;
    }
    return n;
}

double perceived_noise_level(const PnlSpectrum& band_levels)
{
    const PnlSpectrum n = noys(band_levels);
    double largest = n[0];
    double sum = 0.0;
    for (std::size_t i = 0; i < nb; i++)
    {
        largest = std::max(largest, n[i]);
        sum += n[i];
    }
    const double total = 0.85 * largest + 0.15 * sum;
    if (!(total > 0.0))
        return -inf;
    return 40.0 + (10.0 / std::log10(2.0)) * std::log10(total);
}

ToneCorrection tone_correction(const PnlSpectrum& spl)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // Arrays are indexed from 0 (band 1 = 50 Hz); the comments use the
    // regulation's band numbers, one higher.

    // Step 1: slopes s(i) = SPL(i) - SPL(i-1) from band 4; band 3 has none, so
    // the comparisons that would need it (step 2 at band 4) are not made.
    std::array<double, nb> s;
    s.fill(nan);
    for (std::size_t i = 3; i < nb; i++)
        s[i] = spl[i] - spl[i - 1];

    // Steps 2 and 3: where the slope changes by more than 5 dB, encircle
    // SPL(i) if the slope is positive and rising, or SPL(i-1) if it has turned
    // from positive to zero or negative.
    std::array<bool, nb> encircled {};
    for (std::size_t i = 4; i < nb; i++)
    {
        if (!(std::fabs(s[i] - s[i - 1]) > 5.0))
            continue;
        if (s[i] > 0.0 && s[i] > s[i - 1])
            encircled[i] = true;
        if (s[i] <= 0.0 && s[i - 1] > 0.0)
            encircled[i - 1] = true;
    }

    // Step 4: adjusted levels SPL': an encircled band takes the mean of its
    // neighbors, band 24 the extrapolation SPL(23) + s(23).
    PnlSpectrum adjusted = spl;
    for (std::size_t i = 3; i < nb - 1; i++)
    {
        if (encircled[i])
            adjusted[i] = 0.5 * (spl[i - 1] + spl[i + 1]);
    }
    if (encircled[nb - 1])
        adjusted[nb - 1] = spl[nb - 2] + s[nb - 2];

    // Step 5: new slopes s', with s'(3) = s'(4) and an imaginary 25th band,
    // s'(25) = s'(24).
    std::array<double, nb + 1> s1 {};
    for (std::size_t i = 3; i < nb; i++)
        s1[i] = adjusted[i] - adjusted[i - 1];
    s1[2] = s1[3];
    s1[nb] = s1[nb - 1];

    // Step 6: the mean of three adjacent slopes, bands 3 to 23.
    std::array<double, nb - 1> sbar {};
    for (std::size_t i = 2; i < nb - 1; i++)
        sbar[i] = (s1[i] + s1[i + 1] + s1[i + 2]) / 3.0;

    // Step 7: the background SPL'', from SPL''(3) = SPL(3).
    PnlSpectrum background {};
    background[2] = spl[2];
    for (std::size_t i = 3; i < nb; i++)
        background[i] = background[i - 1] + sbar[i - 1];

    // Steps 8 to 10: F = SPL - SPL'' from band 3, the Table A36-2 factor, and
    // the largest.  500 Hz belongs to the middle range.  Bands 1 and 2 are not
    // eligible, and band 3's F is zero by construction.
    ToneCorrection result;
    for (std::size_t i = 2; i < nb; i++)
    {
        const double f = spl[i] - background[i];
        const double frequency = pnl_band_centers[i];
        const bool low = frequency < 500.0 || frequency > 5000.0;
        double c = 0.0;
        if (low)
            c = f >= 20.0 ? 10.0 / 3.0 : f >= 3.0 ? f / 6.0 : f >= 1.5 ? f / 3.0 - 0.5 : 0.0;
        else
            c = f >= 20.0 ? 20.0 / 3.0 : f >= 3.0 ? f / 3.0 : f >= 1.5 ? 2.0 * f / 3.0 - 1.0 : 0.0;
        // Strictly greater: the first of equal factors wins, as argmax does.
        if (c > result.factor)
        {
            result.factor = c;
            result.band = i;
        }
    }
    return result;
}

EpnlResult effective_perceived_noise_level(const std::vector<PnlSpectrum>& history, double dt,
        bool bandshare_adjustment)
{
    if (history.empty())
        throw std::invalid_argument("effective_perceived_noise_level: empty history");
    if (!(dt > 0.0))
        throw std::invalid_argument("effective_perceived_noise_level: the step must be positive");
    const std::size_t count = history.size();
    EpnlResult r;
    r.pnlt.resize(count);
    r.pnl.resize(count);
    r.tone.resize(count);
    r.tone_band.resize(count);
    bool any_finite = false;
    for (std::size_t k = 0; k < count; k++)
    {
        r.pnl[k] = perceived_noise_level(history[k]);
        const ToneCorrection c = tone_correction(history[k]);
        r.tone[k] = c.factor;
        r.tone_band[k] = c.band;
        r.pnlt[k] = r.pnl[k] + c.factor;
        any_finite = any_finite || std::isfinite(r.pnlt[k]);
    }
    if (!any_finite)
        throw std::invalid_argument("effective_perceived_noise_level: no finite PNLT in the history");

    // PNLTM at the first maximum, NaN skipped (nanargmax).
    std::size_t peak = count;
    for (std::size_t k = 0; k < count; k++)
    {
        if (std::isnan(r.pnlt[k]))
            continue;
        if (peak == count || r.pnlt[k] > r.pnlt[peak])
            peak = k;
    }
    r.peak = peak;
    r.pnltm = r.pnlt[peak];
    if (bandshare_adjustment)
    {
        // A36.4.4.2: a tone correction at PNLTM below the mean of the five
        // records centered there suggests the tone was split between bands;
        // PNLTM is then PNL there plus that mean.  Near either end of the
        // history the mean is over the records there are.
        const std::size_t lo = peak >= 2 ? peak - 2 : 0;
        const std::size_t hi = std::min(count, peak + 3);
        double sum = 0.0;
        for (std::size_t k = lo; k < hi; k++)
            sum += r.tone[k];
        const double mean = sum / static_cast<double>(hi - lo);
        if (r.tone[peak] < mean)
        {
            r.pnltm = r.pnl[peak] + mean;
            r.bandshare_adjusted = true;
        }
    }

    // The duration: from the first to the last record at or above PNLTM - 10,
    // dips inside it included (A36.4.5).
    const double threshold = r.pnltm - 10.0;
    std::size_t first = count;
    std::size_t last = 0;
    for (std::size_t k = 0; k < count; k++)
    {
        if (r.pnlt[k] >= threshold)
        {
            if (first == count)
                first = k;
            last = k;
        }
    }
    r.first = first;
    r.last = last;
    double sum = 0.0;
    for (std::size_t k = first; k <= last; k++)
        sum += std::pow(10.0, r.pnlt[k] / 10.0);
    r.duration_correction = 10.0 * std::log10(sum * dt / 10.0) - r.pnltm;
    r.epnl = r.pnltm + r.duration_correction;
    r.clipped = first == 0 || last == count - 1;
    return r;
}

std::array<std::size_t, pnl_band_count> pnl_band_indices(const std::vector<double>& band_centers,
        std::vector<double>& missing)
{
    missing.clear();
    std::array<std::size_t, pnl_band_count> index {};
    for (std::size_t i = 0; i < nb; i++)
    {
        const double nominal = pnl_band_centers[i];
        std::size_t found = band_centers.size();
        for (std::size_t b = 0; b < band_centers.size(); b++)
        {
            if (std::fabs(band_centers[b] / nominal - 1.0) <= 0.05)
            {
                found = b;
                break;
            }
        }
        if (found == band_centers.size())
            missing.push_back(nominal);
        index[i] = found;
    }
    return index;
}

} // namespace niceops
