#include "atmosphere.hpp"

#include <cmath>

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Frequency-independent terms of the standard's expression. A profile builds
// a spectrum at every height, so compute these once for the whole band array.
struct AbsorptionTerms
{
    double relaxation_nitrogen, relaxation_oxygen;
    double classical, oxygen_scale, nitrogen_scale, temperature_scale;

    explicit AbsorptionTerms(const Atmosphere& atmosphere)
        : relaxation_nitrogen(atmosphere.relaxation_frequency_nitrogen()),
          relaxation_oxygen(atmosphere.relaxation_frequency_oxygen())
    {
        const double t_ratio = atmosphere.temperature / reference_temperature;
        classical = 1.84e-11 * (reference_pressure / atmosphere.pressure) * std::sqrt(t_ratio);
        oxygen_scale = 0.01275 * std::exp(-2239.1 / atmosphere.temperature);
        nitrogen_scale = 0.1068 * std::exp(-3352.0 / atmosphere.temperature);
        temperature_scale = std::pow(t_ratio, -2.5);
    }

    // Keep scalar and batch evaluation on the same compiled expression. GCC
    // otherwise contracts them differently after inlining into their callers,
    // breaking the exact scalar/batch agreement even though only ULPs differ.
#if defined(__GNUC__) || defined(__clang__)
    [[gnu::noinline]]
#endif
    double attenuation(double frequency) const
    {
        const double f2 = frequency * frequency;
        const double oxygen = oxygen_scale / (relaxation_oxygen + f2 / relaxation_oxygen);
        const double nitrogen = nitrogen_scale / (relaxation_nitrogen + f2 / relaxation_nitrogen);
        return 8.686 * f2 * (classical + temperature_scale * (oxygen + nitrogen));
    }
};

} // namespace

namespace niceops {

double Atmosphere::sound_speed() const
{
    return 343.2 * std::sqrt(temperature / reference_temperature);
}

double Atmosphere::saturation_pressure() const
{
    return reference_pressure
            * std::pow(10.0,
                    -6.8346 * std::pow(triple_point_temperature / temperature, 1.261) + 4.6151);
}

double Atmosphere::molar_concentration_water_vapor() const
{
    return relative_humidity * saturation_pressure() / pressure;
}

double Atmosphere::relaxation_frequency_oxygen() const
{
    const double h = molar_concentration_water_vapor();
    return pressure / reference_pressure * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
}

double Atmosphere::relaxation_frequency_nitrogen() const
{
    const double h = molar_concentration_water_vapor();
    const double t_ratio = temperature / reference_temperature;
    return pressure / reference_pressure * std::pow(t_ratio, -0.5)
            * (9.0 + 280.0 * h * std::exp(-4.170 * (std::pow(t_ratio, -1.0 / 3.0) - 1.0)));
}

double Atmosphere::attenuation_coefficient(double frequency) const
{
    return AbsorptionTerms(*this).attenuation(frequency);
}

std::vector<double> absorption_coefficients(const Atmosphere& atmosphere,
        const std::vector<double>& frequencies)
{
    const AbsorptionTerms terms(atmosphere);

    std::vector<double> coefficients;
    coefficients.reserve(frequencies.size());
    for (const double frequency : frequencies)
    {
        coefficients.push_back(terms.attenuation(frequency));
    }
    return coefficients;
}

std::string temperature_problem(double kelvin)
{
    if (!std::isfinite(kelvin) || kelvin < 150.0 || kelvin > 350.0)
    {
        return "a temperature of " + std::to_string(kelvin) + " K is not air (150 to 350 K are accepted)"
                + (kelvin >= -100.0 && kelvin <= 70.0 ? "; is it in Celsius? It is read in kelvin." : ".");
    }
    return "";
}

std::string pressure_problem(double kpa)
{
    if (!std::isfinite(kpa) || !(kpa > 0.0) || kpa > 120.0)
    {
        return "a pressure of " + std::to_string(kpa) + " kPa is not plausible (up to 120 kPa are accepted)"
                + (kpa > 120.0 ? "; is it in Pa or hPa? It is read in kPa." : ".");
    }
    return "";
}

std::string humidity_problem(double percent)
{
    if (!std::isfinite(percent) || percent < 0.0 || percent > 100.0)
        return "a relative humidity of " + std::to_string(percent) + " percent is outside 0 to 100.";
    if (percent > 0.0 && percent <= 1.0)
    {
        return "a relative humidity of " + std::to_string(percent) + " percent is read in percent, not as a fraction; "
                "scale a fraction (0 to 1) by 100. (0 for completely dry air is accepted.)";
    }
    return "";
}

} // namespace niceops
