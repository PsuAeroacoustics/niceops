#pragma once

// How the sound a rotorcraft sends in two directions adds, where two rays
// leave it apart (THEORY 6.5).  The spheres carry levels, not phases, so a pair
// of rays' coherency is a model.
//
// By default nothing is assumed of the source: rays add in energy, except
// those leaving within a cone of `cone` degrees of each other, which are one
// emission -- a flush microphone's direct ray and the floor's, a diffracted
// ray and the ray it stands in for at a boundary -- and add with their phases
// (tapered over a quarter of the angle either side, so the field stays
// continuous).  Measured on rotorcraft, nothing simple adds A-weighted levels
// better than that.
//
// Where a database records its vehicle's coherence, the measured model takes
// its place.  For two directions theta degrees apart, in a band at f Hz:
//
//   mean coherency   g(f theta) = exp(-(f theta / scale)^exponent), times,
//                    in bands up to rotation_limit, the delay a pattern
//                    turning with the main rotor puts between them:
//                    sense * (their azimuth difference about the shaft,
//                    wrapped to +-pi) / rotation_rate -- only where the
//                    record gives the rotor
//   coherence        gamma^2 = T + (1 - T) exp(-(f theta / spread_scale)^spread_exponent),
//                    T the band's coherence between distant directions
//   spread           max(gamma^2 - g^2, 0): the mean square of a pair's
//                    coherency about the mean, which sets the interval (the
//                    cone's too, with g its own)
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

// How rays leaving the source apart add: in energy beyond a cone (the
// default), as a database's record says (or the measured defaults), all with
// their phases at zero offset, or all in energy.
enum class RaySum
{
    cone,
    measured,
    coherent,
    energy
};

struct SourceCoherence
{
    RaySum sum = RaySum::cone;
    // The cone's half-angle, degrees.
    double cone = default_cone;
    // The main rotor's pattern, from a database's record: rad/s (0 for none),
    // +1 counterclockwise seen from above the tip-path plane, -1 clockwise,
    // and the band centre, Hz, up to which its delay holds.
    double rotation_rate = 0.0;
    int rotation_sense = 1;
    double rotation_limit = 0.0;
    // The measured mean coherency's decay, Hz deg: the fits with the
    // rotation removed and without (unremoved, the low harmonics' rotation
    // reads as a faster decay).
    double scale = default_scale_without_rotation;
    double exponent = default_exponent_without_rotation;
    // The coherence's decay toward T, Hz deg.
    double spread_scale = 2120.0;
    double spread_exponent = 1.0;
    // T per database band; empty for the default law, 1 / (1 + (f / 58 Hz)^1.67).
    std::vector<double> tonal_fraction;
    // Where each value came from, for the run's log: the sum's, and the
    // spread's and tonal fractions' (the phase interval's).
    std::string source = "the default";
    std::string spread_source = "the defaults";

    static constexpr double default_cone = 1.0, cone_taper = 0.25;
    static constexpr double default_scale = 850.0, default_exponent = 1.6;
    static constexpr double default_scale_without_rotation = 475.0, default_exponent_without_rotation = 2.3;
    static constexpr double default_rotation_harmonics = 4.5;

    bool rotates(double band) const
    {
        return sum == RaySum::measured && rotation_rate > 0.0 && band <= rotation_limit;
    }
    // The mean coherency of two directions `degrees` apart in a band centred
    // at `frequency`: the cone's or the measured decay, 1 coherent, 0 energy.
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
// directions, all optional root attributes (MANUAL, "Source coherence"):
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
    // Whether it says anything a run would use: the rotor or a fit.
    bool specifies() const
    {
        return blades || rotation_sense || fitted();
    }
    bool empty() const
    {
        return !specifies() && source.empty();
    }
};

// The cone: rays leaving within `degrees` of each other add with their
// phases, the rest in energy.  The phase interval's spread and tonal
// fractions are a database's record's where it gives them (they measure the
// coherence the interval is about, whatever the sum), the defaults
// otherwise.  Throws std::invalid_argument for a fit's scale without its
// exponent or the reverse.
SourceCoherence cone_coherence(double degrees = SourceCoherence::default_cone, const CoherenceRecord& record = {});

// The measured model from a database's record and the rotor's rate (rad/s;
// 0 where unknown): the rotor's delay where the record gives the blades and
// the sense, and the record's fits, with the B407's fits for what it leaves
// out -- with the delay removed where it has a rotor, without otherwise.
// An empty record gives the measured defaults without a rotor.  `source`
// says where each part came from.  Throws std::invalid_argument for a fit's
// scale without its exponent or the reverse.
SourceCoherence measured_coherence(const CoherenceRecord& record = {}, double rotation_rate = 0.0);

// "counterclockwise" or "clockwise" (seen from above) as +1 or -1; nullopt
// for anything else.
std::optional<int> parse_rotation_sense(const std::string& text);

} // namespace niceops
