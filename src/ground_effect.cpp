#include "ground_effect.hpp"

#include <cmath>
#include <complex>

#include "faddeeva.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Chessell's coefficients for the broadband (third-octave averaged) form.
constexpr double spherical_spreading_coefficient = 0.727477;
constexpr double ground_reflection_coefficient = 6.325159;

// Beyond this the boundary-loss term is negligible and its series is not worth
// evaluating; panam uses the same cut, leaving the term at zero above it.
constexpr double boundary_loss_cutoff = 500.0;

} // namespace

namespace niceops {

std::vector<std::complex<double>> ground_impedance(const GroundModel& ground,
        const std::vector<double>& frequencies)
{
    using complex = std::complex<double>;
    const complex i(0.0, 1.0);

    std::vector<complex> impedance;
    impedance.reserve(frequencies.size());
    for (const double frequency : frequencies)
    {
        const double resistance_ratio = frequency / ground.flow_resistance;
        const double inverse = std::pow(resistance_ratio, -0.73);
        impedance.push_back(1.0
                / (1.0 + 9.08 * inverse / std::pow(resistance_ratio, 0.02) + i * 11.9 * inverse));
    }
    return impedance;
}

void excess_ground_attenuation(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies,
        double* attenuation)
{
    using complex = std::complex<double>;
    const complex i(0.0, 1.0);

    const double hs = source_height;
    const double hr = ground.receiver_height;

    // Direct path, and the path from the image source below the ground plane.
    const double direct_range = std::hypot(ground_distance, hs - hr);
    const double image_range = std::hypot(ground_distance, hs + hr);
    // Cosine of the grazing angle, straight from the geometry: the arccos and
    // cos panam does round-trip through an angle for nothing.
    const double cos_grazing = image_range > 0.0 ? (hs + hr) / image_range : 1.0;
    const double path_delay = (image_range - direct_range) / sound_speed;
    const double range_ratio = direct_range > 0.0 ? image_range / direct_range : 1.0;

    for (std::size_t band = 0; band < frequencies.size(); band++)
    {
        const double frequency = frequencies[band];
        const complex impedance_ratio = impedance[band];

        // Plane-wave reflection coefficient.
        const complex plane_wave = (cos_grazing - impedance_ratio) / (cos_grazing + impedance_ratio);

        complex image_coefficient = plane_wave;
        if (ground.boundary_loss_correction)
        {
            // Numerical distance.
            const complex parameter = std::sqrt(i * M_PI * frequency * image_range / sound_speed
                                              / (1.0 + impedance_ratio * cos_grazing))
                    * (cos_grazing + impedance_ratio);
            const complex w = parameter * parameter;

            complex boundary_loss = 0.0;
            if (std::abs(w) <= boundary_loss_cutoff)
            {
                // panam writes this as 1 + i sqrt(pi w) exp(-w)(1 - erf(-i sqrt w)).
                // The trailing product is the Faddeeva function at sqrt(w),
                // which is the numerically stable way to evaluate it.
                boundary_loss = 1.0 + i * std::sqrt(M_PI * w) * faddeeva(std::sqrt(w));
            }
            image_coefficient = plane_wave + boundary_loss * (1.0 - plane_wave);
        }

        const double image_magnitude = std::abs(image_coefficient);
        const double image_phase = std::arg(image_coefficient);

        // Chessell equations 20 and 21: the band-averaged combination of the
        // direct and reflected paths, rather than the pure-tone interference.
        const double frequency_delay = frequency * path_delay;
        const double normalized = image_magnitude / range_ratio;
        const double normalized_squared = normalized * normalized;

        // Turbulence decorrelates the two paths, which fills in the nulls.
        const double turbulence_factor = ground.turbulence > 0.0
                ? std::exp(-std::pow(
                          0.5 * ground.turbulence * frequency * std::sqrt(direct_range), 2.0))
                : 1.0;

        const double ground_phase
                = ground_reflection_coefficient * frequency_delay + image_phase;
        const double spherical_phase = spherical_spreading_coefficient * frequency_delay;
        const double cosine_factor = std::cos(ground_phase) * turbulence_factor;

        double magnitude;
        if (path_delay > 0.0)
        {
            // sin(x)/x, the band-averaging factor. It is the reason the delay
            // has to be strictly positive: at zero delay the two paths are the
            // same path and the ratio is 0/0.
            magnitude = 1.0 + normalized_squared
                    + 2.0 * normalized * std::sin(spherical_phase) * cosine_factor
                            / spherical_phase;
        }
        else
        {
            magnitude = 1.0 + normalized_squared + 2.0 * normalized * cosine_factor;
        }

        if (!(magnitude > 0.0))
        {
            // Complete cancellation is not physical here, and taking the log of
            // it would poison the whole spectrum with a NaN or -inf.
            input_error("Excess ground attenuation came out non-positive (", magnitude,
                    ") at ", frequency, " Hz for a source ", hs, " m up and ", ground_distance,
                    " m away. Check the flow resistance (", ground.flow_resistance,
                    " kPa s/m^2) and receiver height (", hr, " m).");
        }
        attenuation[band] = 10.0 * std::log10(magnitude);
    }
}

std::vector<double> excess_ground_attenuation(const GroundModel& ground, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies)
{
    std::vector<double> attenuation(frequencies.size());
    excess_ground_attenuation(ground, ground_impedance(ground, frequencies), source_height,
            ground_distance, sound_speed, frequencies, attenuation.data());
    return attenuation;
}

} // namespace niceops
