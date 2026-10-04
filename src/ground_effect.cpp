#include "ground_effect.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

#include "faddeeva.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Chessell's coefficients for the broadband (third-octave averaged) form.
constexpr double spherical_spreading_coefficient = 0.727477;
constexpr double ground_reflection_coefficient = 6.325159;

// |w| up to which reflection_coefficient_bounds() samples Q; beyond it the
// bound is analytic.  panam once cut the boundary-loss term to zero here,
// because exp(-w) and erfc evaluated apart overflow; the Faddeeva form does
// not, and both now evaluate the term at any |w| (panam e978703).
constexpr double bound_sample_edge = 500.0;

// The spherical-wave reflection coefficient Q of the image source, at grazing
// cosine `cos_grazing`, for a band whose Delany-Bazley impedance is (zr, zi).
// `kappa` is k R / 2 at the image range, pi f R / a, which sets the numerical
// distance.  One copy, shared by the model and by ground_gain_bounds(), so the
// bound is on exactly what the model computes.
inline std::array<double, 2> reflection_coefficient(double cos_grazing, double kappa, double zr,
        double zi, bool boundary_loss)
{
    const double root_pi = std::sqrt(pi);
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

    if (boundary_loss)
    {
        // Numerical distance: sqrt(i k / (1 + Z cos)) * (cos + Z), with
        // k = pi f R / a.  panam squares this to get w and then takes
        // sqrt(w) again inside the boundary loss, but sqrt(w) is exactly
        // this wherever its real part is positive, so the round trip is
        // skipped.  That is every passive ground but a mass-like one
        // (Im Z > 0 here, which a thin layer gives in some bands): there the
        // parameter is in the second quadrant and the principal root is its
        // negative, whose w() carries a surface-wave term growing as
        // exp(|Re w|) with range.  Such a ground has no surface wave, and this
        // keeps the branch without one.
        const double k = kappa;
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

        // panam writes this as 1 + i sqrt(pi w) exp(-w)(1 - erf(-i sqrt w)).
        // sqrt(pi w) is sqrt(pi) times the parameter, and the trailing product
        // is the Faddeeva function -- the stable form, at any |w|.  For a
        // ground whose parameter falls below the real axis (at grazing,
        // -Im Z > Re Z: most soils at low frequency) w() holds the surface wave,
        // 2 exp(-z^2), of size at most one since the parameter stays within
        // 45 degrees of the real axis.
        const auto faddeeva_value = faddeeva({param_r, param_i});
        const double fr = faddeeva_value.real();
        const double fi = faddeeva_value.imag();
        const double product_r = param_r * fr - param_i * fi;
        const double product_i = param_r * fi + param_i * fr;
        // 1 + i*root_pi*product
        const double loss_r = 1.0 - root_pi * product_i;
        const double loss_i = root_pi * product_r;

        // image = plane + loss * (1 - plane)
        const double one_minus_r = 1.0 - plane_r;
        const double one_minus_i = -plane_i;
        image_r = plane_r + loss_r * one_minus_r - loss_i * one_minus_i;
        image_i = plane_i + loss_r * one_minus_i + loss_i * one_minus_r;
    }
    return {image_r, image_i};
}

} // namespace

namespace niceops {

std::vector<std::complex<double>> ground_impedance(const GroundModel& ground,
        const std::vector<double>& frequencies, double sound_speed)
{
    using complex = std::complex<double>;
    const complex i(0.0, 1.0);
    constexpr double gamma = 1.4;          // ratio of specific heats
    constexpr double air_density = 1.2;    // kg/m^3

    std::vector<complex> impedance;
    impedance.reserve(frequencies.size());
    for (const double frequency : frequencies)
    {
        const double resistance_ratio = frequency / ground.flow_resistance;
        complex z;
        switch (ground.impedance_model)
        {
        case ImpedanceModel::delany_bazley:
        {
            // The original form (bit-identical to before the other models).
            const double inverse = std::pow(resistance_ratio, -0.73);
            z = 1.0 + 9.08 * inverse / std::pow(resistance_ratio, 0.02) + i * 11.9 * inverse;
            break;
        }
        case ImpedanceModel::miki:
        {
            const double inverse = std::pow(resistance_ratio, -0.632);
            z = 1.0 + 5.50 * inverse + i * 8.43 * inverse;
            break;
        }
        case ImpedanceModel::variable_porosity:
            z = (1.0 + i) * std::sqrt(1000.0 * ground.flow_resistance
                                      / (pi * gamma * air_density * frequency))
                    + i * sound_speed * ground.porosity_rate / (8.0 * pi * gamma * frequency);
            break;
        case ImpedanceModel::delany_bazley_layer:
        {
            const complex zc = 1.0 + 9.08 * std::pow(resistance_ratio, -0.754)
                    + i * 11.9 * std::pow(resistance_ratio, -0.732);
            const complex k = 2.0 * pi * frequency / sound_speed
                    * (1.0 + 10.8 * std::pow(resistance_ratio, -0.70)
                            + i * 10.3 * std::pow(resistance_ratio, -0.595));
            z = zc / std::tanh(-i * k * ground.layer_depth);
            break;
        }
        }
        const complex beta = 1.0 / z;
        if (!(beta.real() > 0.0))
        {
            input_error("The ground impedance model is not passive at ", frequency,
                    " Hz (Re(1/Z) = ", beta.real(), "). A thin, low-resistivity "
                    "Delany-Bazley layer does this at low frequency; use a deeper layer, a "
                    "higher flow resistance, or another impedance model.");
        }
        impedance.push_back(beta);
    }
    return impedance;
}

ImpedanceModel parse_impedance_model(const std::string& name)
{
    if (name == "delany_bazley") return ImpedanceModel::delany_bazley;
    if (name == "miki") return ImpedanceModel::miki;
    if (name == "variable_porosity") return ImpedanceModel::variable_porosity;
    if (name == "delany_bazley_layer") return ImpedanceModel::delany_bazley_layer;
    input_error("Unknown impedance model '", name, "': use delany_bazley, miki, "
                "variable_porosity or delany_bazley_layer.");
    return ImpedanceModel::delany_bazley;
}

TurbulenceModel parse_turbulence_model(const std::string& name)
{
    if (name == "chessell") return TurbulenceModel::chessell;
    if (name == "harmonoise") return TurbulenceModel::harmonoise;
    input_error("Unknown turbulence model '", name, "': use chessell or harmonoise.");
    return TurbulenceModel::chessell;
}

void excess_ground_attenuation(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies,
        double* attenuation, std::size_t first_band, std::size_t last_band, bool as_power_ratio)
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
    // Chessell: exp(-(scale f)^2), scale = C/2 sqrt(R).  HARMONOISE:
    // exp(-3/8 B k^2 gamma_p^(5/3) R gamma_T) = exp(-(scale f)^2) with
    // scale^2 = 3/8 B (2 pi / c)^2 gamma_p^(5/3) R gamma_T.
    double turbulence_scale = 0.0;
    if (ground.turbulence_model == TurbulenceModel::chessell)
    {
        turbulence_scale
                = ground.turbulence > 0.0 ? 0.5 * ground.turbulence * std::sqrt(direct_range) : 0.0;
    }
    else if (ground.gamma_t > 0.0 && hs + hr > 0.0)
    {
        const double gamma_p = hs * hr / (hs + hr);
        const double k_per_hz = 2.0 * pi / sound_speed;
        turbulence_scale = std::sqrt(0.375 * 0.364 * k_per_hz * k_per_hz
                                     * std::pow(gamma_p, 5.0 / 3.0) * direct_range * ground.gamma_t);
    }
    // Roughness: the coherent reflection's amplitude falls by
    // exp(-2 (k sigma cos theta)^2) = exp(-(roughness_scale f)^2).
    const double roughness_scale = ground.roughness > 0.0
            ? std::sqrt(2.0) * (2.0 * pi / sound_speed) * ground.roughness * cos_grazing
            : 0.0;
    const bool apply_roughness = roughness_scale > 0.0;
    const double inverse_range_ratio = 1.0 / range_ratio;
    const double numerical_distance_scale = pi * image_range / sound_speed;
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

        const auto image = reflection_coefficient(cos_grazing,
                numerical_distance_scale * frequency, zr, zi, apply_boundary_loss);
        double image_r = image[0];
        double image_i = image[1];
        if (apply_roughness)
        {
            const double coherent_part
                    = std::exp(-(roughness_scale * frequency) * (roughness_scale * frequency));
            image_r *= coherent_part;
            image_i *= coherent_part;
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
        attenuation[band] = as_power_ratio ? magnitude : ten_over_ln10 * std::log(magnitude);
    }
}

std::array<double, 2> spherical_reflection_coefficient(double cos_grazing, double kappa,
        std::complex<double> admittance, bool boundary_loss)
{
    return reflection_coefficient(cos_grazing, kappa, admittance.real(), admittance.imag(),
            boundary_loss);
}

void add_bounce_loss(double grazing, int bounces,
        const std::vector<std::complex<double>>& admittance, std::size_t first_band,
        std::size_t last_band, double* loss)
{
    if (bounces <= 0)
        return;
    const double sin_grazing = std::sin(grazing);
    last_band = std::min(last_band, admittance.size());
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const auto q = reflection_coefficient(sin_grazing, 0.0, admittance[band].real(),
                admittance[band].imag(), false);
        loss[band] -= bounces * 10.0 * std::log10(q[0] * q[0] + q[1] * q[1]);
    }
}

std::vector<double> ground_gain_bounds(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, std::size_t first_band,
        std::size_t last_band)
{
    auto bounds = reflection_coefficient_bounds(ground, impedance, first_band, last_band);
    for (std::size_t band = first_band; band < std::min(last_band, bounds.size()); band++)
        bounds[band] = (1.0 + bounds[band]) * (1.0 + bounds[band]);
    return bounds;
}

std::vector<double> reflection_coefficient_bounds(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, std::size_t first_band,
        std::size_t last_band)
{
    // The model's magnitude is 1 + |Q|^2/rho^2 + cross, where rho >= 1 is the
    // image-to-direct range ratio and |cross| <= 2|Q|/rho (the turbulence and
    // sin(x)/x factors are each at most 1 in size), so it is at most
    // (1 + |Q|)^2, whatever the geometry.  What is left is sup |Q| per band.
    //
    // Q depends on the geometry only through the grazing cosine and kappa =
    // pi f R / a.  Up to |w| = 500, |Q| climbs past 1 toward grazing
    // incidence -- the surface wave -- to about 2.02 at the edge for an
    // ordinary soil, which is +9.6 dB, not the +6 dB of pressure doubling.  The
    // sweep samples each grazing cosine up to that edge.  Against a sweep 40x
    // denser in both directions, over flow resistances 10 to 200,000 kPa s/m^2
    // and every band, it came within 0.08%; the 1% margin covers that.
    //
    // Beyond the edge the bound is analytic.  There Q = R_p + F (1 - R_p) with
    // |R_p| <= 1 for a passive ground (Re Z > 0), and F = 1 - S + P: S the
    // asymptotic series of the Faddeeva function, |1 - S| <= 1.003/(2|w|), and
    // P = 2i sqrt(pi) z exp(-z^2) the surface wave, present only where the
    // parameter z is below the real axis.  The parameter's direction, theta,
    // does not depend on range, so |P| = 2 sqrt(pi) |z| exp(-|z|^2 cos 2 theta)
    // peaks at |z|^2 = 1/(2 cos 2 theta), or at the edge if that is nearer.  For
    // ordinary soils the peak is far inside the edge and P is nothing past it;
    // for a nearly lossless stiff surface (a thin layer's band with theta near
    // -45 degrees) the wave travels on past the edge and the bound grows with
    // it, as it should.
    constexpr int cosines = 40;
    constexpr int kappas = 48;
    constexpr double margin = 1.01;
    last_band = std::min(last_band, impedance.size());
    std::vector<double> bounds(impedance.size(), 1.0);
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double zr = impedance[band].real();
        const double zi = impedance[band].imag();
        double largest = 1.0;
        for (int c = 0; c <= cosines; c++)
        {
            // cos = 0, then 1e-4 to 1 in log steps.
            const double cos_grazing
                    = c == 0 ? 0.0 : std::pow(10.0, -4.0 + 4.0 * (c - 1) / (cosines - 1));
            const double base = std::hypot(1.0 + zr * cos_grazing, zi * cos_grazing);
            const double sum_norm = (cos_grazing + zr) * (cos_grazing + zr) + zi * zi;
            const double edge = bound_sample_edge * base / sum_norm;
            for (int k = 0; k <= kappas; k++)
            {
                const double kappa = k == kappas
                        ? edge * (1.0 - 1e-12)
                        : edge * std::pow(10.0, -3.0 + 3.0 * k / kappas);
                const auto q = reflection_coefficient(cos_grazing, kappa, zr, zi,
                        ground.boundary_loss_correction);
                largest = std::max(largest, std::hypot(q[0], q[1]));
            }
            const auto plane = reflection_coefficient(cos_grazing, 0.0, zr, zi, false);
            const double plane_size = std::hypot(plane[0], plane[1]);
            double beyond = plane_size;
            if (ground.boundary_loss_correction)
            {
                // The parameter at kappa = 1; its direction is theta at every kappa.
                const std::complex<double> beta(zr, zi);
                const auto z = std::sqrt(std::complex<double>(0.0, 1.0) / (1.0 + beta * cos_grazing))
                        * (cos_grazing + beta);
                double wave = 0.0;
                if (z.imag() < 0.0)
                {
                    const double cos_2theta = (z.real() * z.real() - z.imag() * z.imag()) / std::norm(z);
                    if (cos_2theta > 0.0)
                    {
                        const double at = std::max(0.5 / cos_2theta, bound_sample_edge);
                        wave = 2.0 * std::sqrt(pi * at) * std::exp(-at * cos_2theta);
                    }
                    else
                    {
                        // Not passive; ground_impedance() refuses such a ground.
                        wave = std::numeric_limits<double>::infinity();
                    }
                }
                const double one_minus = std::hypot(1.0 - plane[0], plane[1]);
                beyond = plane_size + one_minus * (1.003 / (2.0 * bound_sample_edge) + wave);
            }
            largest = std::max(largest, beyond);
        }
        bounds[band] = margin * largest;
    }
    return bounds;
}

GroundTable::GroundTable(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, double sound_speed,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double r_min, double r_max, int angle_nodes, int nodes_per_decade)
{
    last_band = std::min(last_band, frequencies.size());
    first_band_ = first_band;
    last_band_ = last_band;
    bands_ = last_band > first_band ? last_band - first_band : 0;
    n_s_ = std::max(angle_nodes, 2);
    r_min = std::max(r_min, 1e-3);
    r_max = std::max(r_max, r_min * 1.0001);
    log_r_min_ = std::log(r_min);
    log_r_step_ = std::log(10.0) / std::max(nodes_per_decade, 1);
    n_r_ = std::max(2, static_cast<int>(std::ceil((std::log(r_max) - log_r_min_) / log_r_step_)) + 1);
    values_.assign(static_cast<std::size_t>(n_r_) * n_s_ * bands_, 1.0F);
    std::vector<double> ratio(frequencies.size());
    for (int ir = 0; ir < n_r_; ir++)
    {
        const double r = std::exp(log_r_min_ + ir * log_r_step_);
        for (int is = 0; is < n_s_; is++)
        {
            const double s = static_cast<double>(is) / (n_s_ - 1);
            const double sin_elevation = s * s;
            const double hs = r * sin_elevation;
            const double d = r * std::sqrt(std::max(0.0, 1.0 - sin_elevation * sin_elevation));
            excess_ground_attenuation(ground, impedance, hs, d, sound_speed, frequencies,
                    ratio.data(), first_band, last_band, true);
            std::copy(ratio.begin() + first_band, ratio.begin() + last_band,
                    values_.begin() + (static_cast<std::size_t>(ir) * n_s_ + is) * bands_);
        }
    }
}

void GroundTable::power_ratio(double source_height, double ground_distance, double* ratio) const
{
    const double r = std::hypot(source_height, ground_distance);
    const double s = r > 0.0 ? std::sqrt(std::max(0.0, source_height) / r) : 1.0;
    double fs = s * (n_s_ - 1);
    double fr = (std::log(std::max(r, 1e-300)) - log_r_min_) / log_r_step_;
    fs = std::min(std::max(fs, 0.0), static_cast<double>(n_s_ - 1));
    fr = std::min(std::max(fr, 0.0), static_cast<double>(n_r_ - 1));
    const int is = std::min(static_cast<int>(fs), n_s_ - 2);
    const int ir = std::min(static_cast<int>(fr), n_r_ - 2);
    const double ws = fs - is;
    const double wr = fr - ir;
    const double w00 = (1.0 - wr) * (1.0 - ws);
    const double w01 = (1.0 - wr) * ws;
    const double w10 = wr * (1.0 - ws);
    const double w11 = wr * ws;
    const float* a = values_.data() + (static_cast<std::size_t>(ir) * n_s_ + is) * bands_;
    const float* b = a + static_cast<std::size_t>(n_s_) * bands_;
    double* out = ratio + first_band_;
    for (std::size_t k = 0; k < bands_; k++)
        out[k] = w00 * a[k] + w01 * a[k + bands_] + w10 * b[k] + w11 * b[k + bands_];
}

std::vector<double> excess_ground_attenuation(const GroundModel& ground, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies)
{
    std::vector<double> attenuation(frequencies.size());
    excess_ground_attenuation(ground, ground_impedance(ground, frequencies, sound_speed),
            source_height, ground_distance, sound_speed, frequencies, attenuation.data());
    return attenuation;
}

} // namespace niceops
