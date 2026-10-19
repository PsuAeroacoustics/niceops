#pragma once

// How the sound a rotorcraft sends in two directions adds, where two rays
// leave it apart (THEORY 6.5).  The spheres carry levels, not phases, so what
// a pair of rays shares is taken from flight-test measurements of the band
// coherency between emission directions: its mean, which the level takes,
// and its spread about that mean, which the phase interval takes.
//
// For two directions theta degrees apart, in a band at f Hz:
//
//   mean coherency   g(f theta) = exp(-(f theta / scale)^exponent), times,
//                    in bands up to rotation_limit, the delay a pattern
//                    turning with the main rotor puts between them:
//                    sense * (their azimuth difference about the shaft,
//                    wrapped to +-pi) / rotation_rate
//   coherence        gamma^2 = T + (1 - T) exp(-(f theta / spread_scale)^spread_exponent),
//                    T the band's coherence between distant directions
//   spread           max(gamma^2 - g^2, 0): the mean square of a pair's
//                    coherency about the mean, which sets the interval
//
// Measured on the 2017 B407 passes and hovers and a hexacopter (the
// validation harness's coherence/).  Above about 4.5 times the blade-passage
// frequency the rotation predicts the phase worse than none at all.

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "vec3.hpp"

namespace niceops {

// How rays leaving the source apart add: as measured (the model above), all
// with their phases at zero offset, or in energy.
enum class RaySum
{
    measured,
    coherent,
    energy
};

struct SourceCoherence
{
    RaySum sum = RaySum::measured;
    // The main rotor's pattern: rad/s (0 for none), +1 counterclockwise seen
    // from above the tip-path plane, -1 clockwise, and the band centre, Hz, up
    // to which its delay holds.
    double rotation_rate = 0.0;
    int rotation_sense = 1;
    double rotation_limit = 0.0;
    // The mean coherency's decay, Hz deg.  The defaults are the fits with
    // the rotation and without it: unremoved, the low harmonics' rotation
    // reads as a faster decay.
    double scale = 850.0;
    double exponent = 1.6;
    // The coherence's decay toward T, Hz deg.
    double spread_scale = 2120.0;
    double spread_exponent = 1.0;
    // T per database band; empty for the default law, 1 / (1 + (f / 58 Hz)^1.67).
    std::vector<double> tonal_fraction;
    // Where each value came from, for the run's log.
    std::string source = "the default";

    static constexpr double default_scale = 850.0, default_exponent = 1.6;
    static constexpr double default_scale_without_rotation = 475.0, default_exponent_without_rotation = 2.3;
    static constexpr double default_rotation_harmonics = 4.5;

    bool rotates(double band) const
    {
        return rotation_rate > 0.0 && band <= rotation_limit;
    }
    double mean(double frequency, double degrees) const;
    // T for band `band` (database order) at centre `frequency`.
    double tonal(std::size_t band, double frequency) const;
    double coherence(double frequency, double degrees, double tonal) const;
    double spread(double frequency, double degrees, double tonal) const;
    // The rotor pattern's delay between two launch directions as a length,
    // m, to add to the first ray's: sound_speed * sense * dpsi / rate, dpsi
    // the first's azimuth about `normal` less the second's, wrapped to +-pi.
    // `nose` is any direction in the plane the azimuth is measured from.
    double delay_length(const Vec3& first, const Vec3& second, const Vec3& normal, const Vec3& nose,
            double sound_speed) const;
};

// What a database records of its vehicle's coherence between emission
// directions, all optional root attributes (README, "Source coherence"):
// main_rotor_blades, main_rotor_rotation ("counterclockwise" or "clockwise",
// seen from above), coherence_rotation_limit_Hz, coherence_scale_Hz_deg with
// coherence_exponent, coherence_spread_scale_Hz_deg with
// coherence_spread_exponent, coherence_tonal_fraction (one per band) and
// coherence_source (text).  The fits are the ones the record's own rotation
// implies: with the delay removed where it gives the blades and the sense,
// without it otherwise.
struct CoherenceRecord
{
    std::optional<double> blades;
    std::optional<int> rotation_sense;
    std::optional<double> rotation_limit;
    std::optional<double> scale, exponent;
    std::optional<double> spread_scale, spread_exponent;
    std::vector<double> tonal_fraction;
    std::string source;

    // Whether it holds anything of the coherence itself (not just the rotor).
    bool fitted() const
    {
        return scale || spread_scale || rotation_limit || !tonal_fraction.empty();
    }
    bool empty() const
    {
        return !blades && !rotation_sense && !fitted() && source.empty();
    }
};

// The model from a database's record, the rotor's rate (rad/s; 0 where
// unknown) and the vehicle file's blades and sense (0 where it gives none).
// The record's rotor takes precedence; where it holds fits, its rotor is the
// one they were made with, so the vehicle file's is not mixed in.  What the
// record leaves out is the default for the rotation it ends up with.
// `source` says where each part came from.
SourceCoherence source_coherence(const CoherenceRecord& record, double rotation_rate, double vehicle_blades,
        int vehicle_sense);

// The model where nothing says otherwise.  With `rotation_rate` > 0 and the
// blade count, the main rotor's delay to 4.5 blade-passage frequencies; the
// sense is +1 counterclockwise seen from above.
SourceCoherence default_source_coherence(double rotation_rate = 0.0, double blades = 0.0, int sense = 0);

// "counterclockwise" or "clockwise" (seen from above) as +1 or -1; nullopt
// for anything else.
std::optional<int> parse_rotation_sense(const std::string& text);

} // namespace niceops
