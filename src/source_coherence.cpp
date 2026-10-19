#include "source_coherence.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

#include "util/constants.hpp"

namespace niceops {

double SourceCoherence::mean(double frequency, double degrees) const
{
    switch (sum)
    {
    case RaySum::coherent:
        return 1.0;
    case RaySum::energy:
        return 0.0;
    case RaySum::cone:
    {
        // cos^2 from a quarter inside the cone to a quarter outside: a half at
        // its edge.
        const double inner = cone * (1.0 - cone_taper), outer = cone * (1.0 + cone_taper);
        if (degrees <= inner)
            return 1.0;
        if (degrees >= outer)
            return 0.0;
        const double c = std::cos(0.5 * pi * (degrees - inner) / (outer - inner));
        return c * c;
    }
    case RaySum::measured:
        break;
    }
    return std::exp(-std::pow(std::max(frequency * degrees, 0.0) / scale, exponent));
}

double SourceCoherence::tonal(std::size_t band, double frequency) const
{
    if (band < tonal_fraction.size())
        return std::clamp(tonal_fraction[band], 0.0, 1.0);
    return 1.0 / (1.0 + std::pow(frequency / 58.0, 1.67));
}

double SourceCoherence::coherence(double frequency, double degrees, double t) const
{
    return t + (1.0 - t) * std::exp(-std::pow(std::max(frequency * degrees, 0.0) / spread_scale, spread_exponent));
}

double SourceCoherence::spread(double frequency, double degrees, double t) const
{
    const double g = mean(frequency, degrees);
    return std::max(coherence(frequency, degrees, t) - g * g, 0.0);
}

double SourceCoherence::delay_length(const Vec3& first, const Vec3& second, const Vec3& normal, const Vec3& nose,
        double sound_speed) const
{
    if (!(rotation_rate > 0.0))
        return 0.0;
    // Counterclockwise seen from above is positive about the normal: from
    // the nose toward normal x nose.
    const Vec3 left = normal.cross(nose);
    const double psi1 = std::atan2(left.dot(first), nose.dot(first));
    const double psi2 = std::atan2(left.dot(second), nose.dot(second));
    double d = psi1 - psi2;
    if (d > pi)
        d -= 2.0 * pi;
    else if (d < -pi)
        d += 2.0 * pi;
    return sound_speed * rotation_sense * d / rotation_rate;
}

SourceCoherence cone_coherence(double degrees)
{
    SourceCoherence c;
    c.sum = RaySum::cone;
    c.cone = degrees;
    std::ostringstream text;
    text << "energy beyond a cone of " << degrees << " deg";
    c.source = text.str();
    return c;
}

SourceCoherence measured_coherence(const CoherenceRecord& record, double rotation_rate)
{
    SourceCoherence c;
    c.sum = RaySum::measured;
    const double blades = record.blades.value_or(0.0);
    const int sense = record.rotation_sense.value_or(0);
    if (rotation_rate > 0.0 && blades > 0.0 && (sense == 1 || sense == -1))
    {
        c.rotation_rate = rotation_rate;
        c.rotation_sense = sense;
        c.rotation_limit = record.rotation_limit.value_or(
                SourceCoherence::default_rotation_harmonics * blades * rotation_rate / (2.0 * pi));
        c.scale = SourceCoherence::default_scale;
        c.exponent = SourceCoherence::default_exponent;
    }
    if (record.scale)
    {
        c.scale = *record.scale;
        c.exponent = *record.exponent;
    }
    if (record.spread_scale)
    {
        c.spread_scale = *record.spread_scale;
        c.spread_exponent = *record.spread_exponent;
    }
    c.tonal_fraction = record.tonal_fraction;
    c.source = record.fitted() ? "the database's measured coherence" : "the measured defaults";
    if (record.fitted() && !record.source.empty())
        c.source += " (" + record.source + ")";
    c.source += c.rotation_rate > 0.0 ? ", the database's rotor" : ", no rotor delay";
    return c;
}

std::optional<int> parse_rotation_sense(const std::string& text)
{
    std::string t;
    for (char ch : text)
        if (!std::isspace(static_cast<unsigned char>(ch)))
            t += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (t == "counterclockwise" || t == "anticlockwise")
        return 1;
    if (t == "clockwise")
        return -1;
    return std::nullopt;
}

} // namespace niceops
