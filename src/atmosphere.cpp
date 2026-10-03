#include "atmosphere.hpp"

#include <cmath>

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// The body of the standard's absorption expression, with the relaxation
// frequencies passed in rather than recomputed.  Both entry points below go
// through this, so there is one copy of the formula to check against the
// Python.
double attenuation(double pressure, double temperature, double relaxation_nitrogen,
        double relaxation_oxygen, double frequency)
{
    const double t_ratio = temperature / reference_temperature;
    const double f2 = frequency * frequency;

    // Classical and rotational losses, then the vibrational relaxation of
    // oxygen and of nitrogen.
    const double classical = 1.84e-11 * (reference_pressure / pressure) * std::sqrt(t_ratio);
    const double oxygen = 0.01275 * std::exp(-2239.1 / temperature)
            / (relaxation_oxygen + f2 / relaxation_oxygen);
    const double nitrogen = 0.1068 * std::exp(-3352.0 / temperature)
            / (relaxation_nitrogen + f2 / relaxation_nitrogen);

    return 8.686 * f2 * (classical + std::pow(t_ratio, -2.5) * (oxygen + nitrogen));
}

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

double Atmosphere::molar_concentration_water_vapour() const
{
    return relative_humidity * saturation_pressure() / pressure;
}

double Atmosphere::relaxation_frequency_oxygen() const
{
    const double h = molar_concentration_water_vapour();
    return pressure / reference_pressure * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
}

double Atmosphere::relaxation_frequency_nitrogen() const
{
    const double h = molar_concentration_water_vapour();
    const double t_ratio = temperature / reference_temperature;
    return pressure / reference_pressure * std::pow(t_ratio, -0.5)
            * (9.0 + 280.0 * h * std::exp(-4.170 * (std::pow(t_ratio, -1.0 / 3.0) - 1.0)));
}

double Atmosphere::attenuation_coefficient(double frequency) const
{
    return attenuation(pressure, temperature, relaxation_frequency_nitrogen(),
            relaxation_frequency_oxygen(), frequency);
}

std::vector<double> absorption_coefficients(const Atmosphere& atmosphere,
        const std::vector<double>& frequencies)
{
    // Hoisted: these cost a pow() and an exp() each and do not vary with
    // frequency.
    const double relaxation_nitrogen = atmosphere.relaxation_frequency_nitrogen();
    const double relaxation_oxygen = atmosphere.relaxation_frequency_oxygen();

    std::vector<double> coefficients;
    coefficients.reserve(frequencies.size());
    for (const double frequency : frequencies)
    {
        coefficients.push_back(attenuation(atmosphere.pressure, atmosphere.temperature,
                relaxation_nitrogen, relaxation_oxygen, frequency));
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
