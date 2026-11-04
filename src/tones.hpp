#pragma once

// A format 2 database's pure tones in a run (THEORY 4.4): every tone of every
// rotor, at each of its arrivals at an observer, with its phase.
//
// Each source is a rotor whose stored amplitude is its far-field directivity at
// the sphere radius, the phase tied to the rotor's azimuth (Database::
// ToneSource).  An arrival is one straight path from the rotor's hub to the
// observer, or to the observer's image in the ground plane.  Its phasor, in the
// e^(+i omega t) convention the spheres are written in (p ~ A e^(-ikR)/R at
// reception time t):
//
//   a = A(n) (R_s / L) 10^(-alpha(f) L_abs / 20) e^(i (arg A(n) + m psi0 + omega_e delta))
//
// with n the direction leaving the hub at its emission time, L the path length
// from there, f = f_e / (1 - M.n) the received frequency, delta the emission
// time relative to the state's own time for a common reception time, and psi0
// the rotor's azimuth offset (--rotor_azimuth).  In wind every arrival is solved
// in the frame moving with the air (tone_arrival).  The image arrival carries the
// ground's spherical-wave reflection coefficient at its own frequency, which
// enters conjugated (the ground model's Q is written for e^(-i omega t)).
//
// Sources in one coherent group add their phasors harmonic by harmonic; every
// other source adds in energy.  A source's own direct and reflected arrivals add
// with the turbulence coherence on their cross term:
//
//   P = |S_d|^2 + |S_r|^2 + 2 gamma Re(S_d conj(S_r)).
//
// The powers are filed into the database's bands, A-weighted at their own
// received frequency, so they add to the band spectrum's energy (band_content
// broadband).

#include <array>
#include <atomic>
#include <complex>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "ground_effect.hpp"
#include "sphere_database.hpp"
#include "vec3.hpp"

namespace niceops {

// --tone_coherence: sources in a sync group add coherently, or every source in
// energy (each keeping its own arrivals' interference).
enum class ToneCoherence
{
    coherent,
    energy,
};

// --tone_ground: the ground per tone frequency (each arrival's Q at its own
// frequency), or the band model's correction on free-field tones.
enum class ToneGround
{
    frequency,
    band,
};

std::optional<ToneCoherence> parse_tone_coherence(const std::string& name);
std::optional<ToneGround> parse_tone_ground(const std::string& name);

struct ToneRunOptions
{
    ToneCoherence coherence = ToneCoherence::coherent;
    ToneGround ground = ToneGround::frequency;
    // --rotor_azimuth: each named source advanced by this many degrees in its own
    // sense of rotation, which multiplies its shaft order m by e^(i m psi0).
    // Unnamed sources keep the stored rotor-locked phase as it is.
    std::map<std::string, double> rotor_azimuth_deg;
    // --sync_group: coherent groups declared for the run, which replace the
    // database's groups for their members and hold their shaft rates at the
    // group's mean.
    std::vector<std::vector<std::string>> sync_groups;
    // --tone_near_field on: each arrival's tone takes the first-order near-field
    // term, A + (i / (2 k R)) Lap[D A] (THEORY 4.4), where the database holds it
    // (DatabaseLoadOptions::tone_near_field).
    bool near_field = false;
};

// "tag=deg,tag=deg" into `out`; refuses an entry without '=' or a number.
void parse_rotor_azimuths(const std::string& text, std::map<std::string, double>& out);
// "tag,tag,..." into a group's member list; refuses an empty one.
std::vector<std::string> parse_sync_group(const std::string& text);

// One straight arrival from a hub that is at `hub` at the state's time and moves
// at `velocity`: received at `point` at the state's time plus `reception_lag`
// (the same for every arrival of the state, so their emission times compare).
struct ToneArrival
{
    // Emission time relative to the state's time, s.
    double delta = 0.0;
    // Path length from the hub at emission to the point, m.
    double length = 0.0;
    // Unit vector from the hub at emission toward the point.
    Vec3 direction;
    // The Doppler factor 1 - V.n / c: the received frequency is f_e over it.
    double doppler = 1.0;
};
// `wind`: a uniform wind (m/s), in which sound travels at c relative to the air;
// the solve is then the still-air one in the frame moving with it.  Zero, calm
// air, is the same arithmetic as without it.
ToneArrival tone_arrival(const Vec3& hub, const Vec3& velocity, const Vec3& point, double reception_lag,
        double sound_speed, const Vec3& wind = Vec3(0.0, 0.0, 0.0));

// Per band, the ratio of the received energy to the free field's (the band
// path's ground, plate, table, refraction or urban terms), applied to free-field
// tones: the band fallback.  Adds sum(weighted * ratio) to `level` in energy and,
// when `bands_out` is given, power * ratio to each band's level there.  A null
// `ratio` is 1 everywhere (the per-frequency path, whose tones already carry
// their ground).
double add_tone_power(double level, const double* power, const double* weighted, const double* ratio,
        std::size_t first_band, std::size_t last_band, double* bands_out);

class ToneEvaluator
{
public:
    struct Group
    {
        std::string name;
        std::vector<std::size_t> sources;
        // Whether the members add coherently: a sync group of rotor-locked
        // sources under --tone_coherence coherent.
        bool coherent = false;
        // A run-time group (--sync_group), whose shaft rates are held at their
        // mean.
        bool forced = false;
    };

    // `band_gains` are the run's (A-weighting and band limits, a_weighting_gains);
    // a tone counts in the A-weighted level where its band's gain is not zero.
    ToneEvaluator(const Database& database, const ToneRunOptions& options, const std::vector<double>& band_gains,
            double sound_speed);

    // What one state holds for every observer: each source's hub at the state's
    // time and its tones' emission frequencies (held at their group's mean
    // shaft rate for a run-time group), and the velocity.
    struct StateTones
    {
        std::vector<Vec3> hubs;
        // [source][tone], flattened by tone_offset().
        std::vector<double> frequencies;
        // The hubs' ground velocity, and the wind the tones travel in (the
        // profile's, averaged over the heights of the state's straight paths to
        // the ground; zero in calm air or without a profile).
        Vec3 velocity;
        Vec3 wind;
        // Each run-time group's mean shaft rate at this state, rad/s (0 for others).
        std::vector<double> group_rates;
    };
    StateTones prepare(const Database::SourceState& state, const Vec3& velocity,
            const Vec3& wind = Vec3(0.0, 0.0, 0.0)) const;

    // The ground for the per-frequency path: the model and its admittance at any
    // frequency.  Null: free field, direct arrivals only.
    struct Ground
    {
        const GroundModel* model = nullptr;
        const AdmittanceTable* admittance = nullptr;
    };

    // Every tone's power at `observer`, added into `power` (unweighted, p_ref^2)
    // and `weighted` (A-weighted at the tone's own received frequency) by band;
    // both hold the database's band count.  `absorption` is the path's dB/m per
    // band, read at each tone's received frequency by interpolation in log
    // frequency.  `ground_plane` is the ground's height under the observer.
    void evaluate(const Database::SourceState& state, const StateTones& tones, const Vec3& observer,
            const double* absorption, const Ground* ground, double ground_plane, double* power,
            double* weighted) const;

    const std::vector<Group>& groups() const
    {
        return groups_;
    }
    // The largest first-order near-field phase error of the stored spheres,
    // (m + 1) c / (2 |Omega| R), over every source, tone and condition, deg, and
    // the shaft order and source it belongs to.
    double near_field_degrees() const
    {
        return near_field_degrees_;
    }
    int near_field_order() const
    {
        return near_field_order_;
    }
    const std::string& near_field_source() const
    {
        return near_field_source_;
    }
    // Tone arrivals whose received frequency fell outside every band, so they
    // were not counted.
    std::size_t outside_bands() const
    {
        return outside_bands_.load();
    }
    // Tone arrivals that took the near-field term; of them, those where it was
    // more than half the far-field amplitude there (the first-order expansion is
    // then doubtful), and those where it was capped at twice the tone's largest
    // amplitude over the sphere.
    std::size_t near_field_evaluations() const
    {
        return near_field_evaluations_.load();
    }
    std::size_t near_field_large() const
    {
        return near_field_large_.load();
    }
    std::size_t near_field_capped() const
    {
        return near_field_capped_.load();
    }
    std::size_t tone_offset(std::size_t source) const
    {
        return tone_offsets_[source];
    }
    std::size_t band_of(double frequency) const;
    static constexpr std::size_t no_band = static_cast<std::size_t>(-1);

private:
    const Database& database_;
    ToneRunOptions options_;
    std::vector<double> band_gains_;
    double sound_speed_ = 343.2;
    std::vector<Group> groups_;
    // Per source: its group, its azimuth offset (rad), where its tones start in
    // StateTones::frequencies.
    std::vector<std::size_t> group_of_;
    std::vector<double> azimuth_;
    std::vector<std::size_t> tone_offsets_;
    std::size_t tone_total_ = 0;
    std::size_t max_tones_ = 0;
    // log2 of the band edges, edges_[b] .. edges_[b + 1] for band b.
    std::vector<double> log2_edges_;
    std::vector<double> log2_centers_;
    double near_field_degrees_ = 0.0;
    int near_field_order_ = 0;
    std::string near_field_source_;
    mutable std::atomic<std::size_t> outside_bands_ {0};
    mutable std::atomic<std::size_t> near_field_evaluations_ {0};
    mutable std::atomic<std::size_t> near_field_large_ {0};
    mutable std::atomic<std::size_t> near_field_capped_ {0};

    double absorption_at(const double* absorption, double frequency) const;
};

} // namespace niceops
