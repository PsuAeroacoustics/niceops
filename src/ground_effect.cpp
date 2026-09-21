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
        double* attenuation, std::size_t first_band, std::size_t last_band)
{
    // Bands outside the metric's limits are never read, so do not compute them.
    if (last_band > frequencies.size())
        last_band = frequencies.size();

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

    // Hoisted: the turbulence factor's only per-band term is the frequency.
    const double turbulence_scale
            = ground.turbulence > 0.0 ? 0.5 * ground.turbulence * std::sqrt(direct_range) : 0.0;
    const double inverse_range_ratio = 1.0 / range_ratio;
    const double root_pi = std::sqrt(M_PI);
    const double numerical_distance_scale = M_PI * image_range / sound_speed;
    // 10*log10(x) == (10/ln 10) * ln(x), and log() is the cheaper call.
    constexpr double ten_over_ln10 = 4.34294481903251827651;

    // The loop below is real arithmetic throughout rather than std::complex.
    // Every complex multiply and divide in libc++ carries a branch to recover
    // NaN and infinity cases, and there are six of them per band here; none can
    // arise, because every quantity is bounded by construction.
    // Both of these are band-invariant; branching on them inside the loop
    // blocks the compiler from keeping the body straight-line.
    const bool apply_boundary_loss = ground.boundary_loss_correction;
    const bool apply_sinc = path_delay > 0.0;
    const bool apply_turbulence = turbulence_scale > 0.0;

    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double frequency = frequencies[band];
        const double zr = impedance[band].real();
        const double zi = impedance[band].imag();

        // Plane-wave reflection coefficient, (cos - Z) / (cos + Z).
        const double sum_r = cos_grazing + zr;
        const double sum_i = zi;
        const double sum_norm = sum_r * sum_r + sum_i * sum_i;
        const double diff_r = cos_grazing - zr;
        const double diff_i = -zi;
        const double plane_r = (diff_r * sum_r + diff_i * sum_i) / sum_norm;
        const double plane_i = (diff_i * sum_r - diff_r * sum_i) / sum_norm;

        double image_r = plane_r;
        double image_i = plane_i;

        if (apply_boundary_loss)
        {
            // Numerical distance: sqrt(i k / (1 + Z cos)) * (cos + Z), with
            // k = pi f R / a.  panam squares this to get w and then takes
            // sqrt(w) again inside the boundary loss, but sqrt(w) is exactly
            // this -- its real part is positive throughout the domain the model
            // reaches -- so the round trip is skipped.
            const double k = numerical_distance_scale * frequency;
            const double base_r = 1.0 + zr * cos_grazing;
            const double base_i = zi * cos_grazing;
            const double base_norm = base_r * base_r + base_i * base_i;
            // i*k divided by (base_r, base_i)
            const double quotient_r = k * base_i / base_norm;
            const double quotient_i = k * base_r / base_norm;

            // Algebraic complex square root: libc++'s goes through abs() and
            // arg(), costing a hypot and an atan2.  The branch is the standard
            // one that avoids cancellation.
            double root_r;
            double root_i;
            {
                const double magnitude
                        = std::sqrt(quotient_r * quotient_r + quotient_i * quotient_i);
                if (quotient_r >= 0.0)
                {
                    root_r = std::sqrt(0.5 * (magnitude + quotient_r));
                    root_i = quotient_i / (2.0 * root_r);
                }
                else
                {
                    root_i = std::copysign(std::sqrt(0.5 * (magnitude - quotient_r)), quotient_i);
                    root_r = quotient_i / (2.0 * root_i);
                }
            }

            const double param_r = root_r * sum_r - root_i * sum_i;
            const double param_i = root_r * sum_i + root_i * sum_r;

            // |w| <= 500 is |parameter|^2 <= 500, which norm gives without a root.
            double loss_r = 0.0;
            double loss_i = 0.0;
            if (param_r * param_r + param_i * param_i <= boundary_loss_cutoff)
            {
                // panam writes this as 1 + i sqrt(pi w) exp(-w)(1 - erf(-i sqrt w)).
                // sqrt(pi w) is sqrt(pi) times the parameter, and the trailing
                // product is the Faddeeva function -- the stable form.
                const auto faddeeva_value = faddeeva({param_r, param_i});
                const double fr = faddeeva_value.real();
                const double fi = faddeeva_value.imag();
                const double product_r = param_r * fr - param_i * fi;
                const double product_i = param_r * fi + param_i * fr;
                // 1 + i*root_pi*product
                loss_r = 1.0 - root_pi * product_i;
                loss_i = root_pi * product_r;
            }

            // image = plane + loss * (1 - plane)
            const double one_minus_r = 1.0 - plane_r;
            const double one_minus_i = -plane_i;
            image_r = plane_r + loss_r * one_minus_r - loss_i * one_minus_i;
            image_i = plane_i + loss_r * one_minus_i + loss_i * one_minus_r;
        }

        // Chessell equations 20 and 21: the band-averaged combination of the
        // direct and reflected paths, rather than the pure-tone interference.
        const double frequency_delay = frequency * path_delay;
        // |Q|^2 without the root, and |Q|cos(p + arg Q) expanded as
        // Re(Q)cos p - Im(Q)sin p -- the same number without the arg() (an
        // atan2) or the abs() (a hypot) that reach it the long way.
        const double image_norm = image_r * image_r + image_i * image_i;
        const double ground_phase = ground_reflection_coefficient * frequency_delay;
        const double coherent
                = image_r * std::cos(ground_phase) - image_i * std::sin(ground_phase);

        // Turbulence decorrelates the two paths, which fills in the nulls.
        const double turbulence_factor = apply_turbulence
                ? std::exp(-(turbulence_scale * frequency) * (turbulence_scale * frequency))
                : 1.0;

        double cross = 2.0 * inverse_range_ratio * coherent * turbulence_factor;
        if (apply_sinc)
        {
            // sin(x)/x, the band-averaging factor. It is the reason the delay
            // has to be strictly positive: at zero delay the two paths are the
            // same path and the ratio is 0/0.
            const double spherical_phase = spherical_spreading_coefficient * frequency_delay;
            cross *= std::sin(spherical_phase) / spherical_phase;
        }

        const double magnitude
                = 1.0 + image_norm * inverse_range_ratio * inverse_range_ratio + cross;

        if (!(magnitude > 0.0))
        {
            // Complete cancellation is not physical here, and taking the log of
            // it would poison the whole spectrum with a NaN or -inf.
            input_error("Excess ground attenuation came out non-positive (", magnitude,
                    ") at ", frequency, " Hz for a source ", hs, " m up and ", ground_distance,
                    " m away. Check the flow resistance (", ground.flow_resistance,
                    " kPa s/m^2) and receiver height (", hr, " m).");
        }
        attenuation[band] = ten_over_ln10 * std::log(magnitude);
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
