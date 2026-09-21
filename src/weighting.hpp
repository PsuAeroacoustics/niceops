#pragma once

// Frequency weighting and band-limited integration of a spectrum.
//
// Ported from panam's dBAw() in flight_acoustics.py -- the IEC 61672 pole-zero
// form, not a table lookup, so it evaluates at whatever band centres a database
// happens to carry.
//
// Band limiting is deliberate rather than incidental.  A database's bands cover
// the range its source data covered, which is not the range a metric is defined
// over: summing whatever bands happen to be present quietly makes the answer a
// property of the database's bandwidth.  Stating the limits makes two runs over
// different databases comparable, and makes it visible when a database cannot
// cover the range asked of it.

#include <cstddef>
#include <vector>

namespace niceops {

// A-weighting in dB at a frequency.  Zero and negative frequencies weight to
// -300 dB, as panam does, rather than producing a NaN from log10(0).
double a_weighting(double frequency);

// The frequency range a metric is integrated over.  The default spans the
// audible range, so it keeps every band of any plausible database until a run
// narrows it.
struct BandLimits
{
    double low = 20.0; // Hz
    double high = 20000.0; // Hz

    bool contains(double frequency) const
    {
        return frequency >= low && frequency <= high;
    }
};

// Per-band linear gains folding the A-weighting and the band limit into a
// single multiply: 10^(A(f)/10) for a band inside the limits, and exactly 0 for
// one outside.
//
// Precomputed because the acoustic loop evaluates this for every observer at
// every trajectory state; the weighting depends only on the band centres.
std::vector<double> a_weighting_gains(const std::vector<double>& frequencies, BandLimits limits);

// Overall level in dB from per-band levels and gains from a_weighting_gains.
//
// Returns -infinity when no band contributes -- every band outside the limits,
// or an empty spectrum.  That is a real answer (silence), not an error, but it
// is also what a band limit excluding the whole database looks like, so callers
// that can say something useful about the difference should check.
double weighted_level(const double* band_levels, const std::vector<double>& gains);

// How many of these band centres fall inside the limits.  For reporting what a
// run actually integrated over.
std::size_t bands_in_range(const std::vector<double>& frequencies, BandLimits limits);

} // namespace niceops
