#pragma once

// A flush ground-plane microphone on a rigid plate lying on the ground -- the
// GRAS 67AX in its 400 mm GR1425 plate, at 3/4 of the radius per SAE ARP 4055
// -- as a receiver for the spectral path.
//
// The point-receiver ground effect (ground_effect.hpp) treats the microphone
// as sitting on one impedance: soft ground at a pole, or a rigid surface at
// the ground plane.  A plate is both.  At high frequency the microphone sees
// the rigid plate (pressure doubling); at grazing incidence and low frequency
// the soft ground around it sets the level.  panam's axisymmetric boundary
// element model (axisymmetric_bem.py; docs/ground_plane_corrections.md) solves
// the plate exactly for a point source over an impedance half-space and writes
// the answer as two complex transfer functions of the direct and image waves:
//
//     p / p_free = P_d + Q P_r
//
// with P_d and P_r functions of frequency and the source's elevation and
// azimuth only (the far field of the ~0.2 m plate), and Q the spherical-wave
// reflection coefficient of the ground at the geometry.  Q is recomputed here
// for each source position, from the same ground model the point receiver
// uses; P_d and P_r depend on the ground too, so the table records the ground
// it was computed on and load() refuses a mismatch.
//
// The level is band averaged as panam's board_level() does: |P_d + Q P_r|^2 at
// five sub-frequencies per third-octave band (the table's), averaged in energy.
// P_d and P_r are interpolated bilinearly, real and imaginary parts, in
// elevation and azimuth; below the table's lowest elevation they are held.
//
// Turbulence and roughness do not enter: the plate's image path is millimetres
// longer than the direct one, so neither decorrelates it.
//
// The table is netCDF, from panam's axisymmetric_bem.write_netcdf():
//   dims band, sub, elevation, azimuth
//   band_centre(band), frequency(band, sub), elevation (deg), azimuth (deg,
//   uniform from 0), P_d_real/_imag and P_r_real/_imag(band, sub, elevation, azimuth)
//   attributes thickness_m (the plate top above the ground, where Q is taken),
//   ground_model, flow_resistance, porosity_rate, layer_depth_m, and the rest
//   of the geometry for the log.
// Azimuth is that of the horizontal propagation direction, source to
// microphone, in the plate's frame, where the microphone is offset from the
// plate's centre along +y.

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <limits>
#include <vector>

#include "ground_effect.hpp"
#include "vec3.hpp"

namespace niceops {

class PlateReceiver
{
public:
    // Read a table.  Not usable until bound().
    static PlateReceiver load(const std::string& path);

    // A copy for one trajectory, checked against the database's band centres,
    // the bands [first_band, last_band) the run uses, and the run's ground
    // model.  `sound_speed` (m/s) is the trajectory's atmosphere's: it sets
    // the admittance at the sub-frequencies (variable porosity and the layer)
    // and Q's numerical distance.  The table itself is shared, not copied.
    PlateReceiver bound(const std::vector<double>& bands, const GroundModel& ground,
            double sound_speed, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1)) const;

    // Level at the microphone re free field, dB, per band in [first_band,
    // last_band) of `bands`, into `correction`.  `source_height` is above the
    // ground, `ground_distance` horizontal, both m.  `azimuth` is the
    // horizontal propagation direction in the plate's frame, radians (see the
    // file comment).
    //
    // Q is taken as quadratic in log frequency across each band, through its
    // values at the band's edges and centre (shared edges: 2n + 1 evaluations
    // for n bands, not 5n).  Against level_exact(), over random geometries on
    // the GR1425 tables: within 0.07 dB wherever the level is above -15 dB re
    // free field; up to ~1.5 dB in the deep near-grazing cancellations, some
    // -30 to -70 dB down, where it cannot move a band sum.  It is ~2.5x cheaper
    // for the whole run.
    void level(double source_height, double ground_distance, double azimuth,
            double* correction, std::size_t first_band, std::size_t last_band) const;

    // Tabulate the fast path's Q values (band edges and centres) over the
    // image path's grazing cosine and length, up to `r_max` metres, so that
    // level() and power_ratio() interpolate them instead of evaluating the
    // reflection coefficient (a Faddeeva evaluation each).  Q varies smoothly
    // with the geometry -- the plate's interference is in P_d and P_r -- so
    // the grid is coarser than GroundTable's.  level_exact() never uses it.
    void tabulate(double r_max, int angle_nodes = 161, int nodes_per_decade = 20);

    // The same as a power ratio, 10^(level/10): what the band sum multiplies
    // by, without a log per band.
    void power_ratio(double source_height, double ground_distance, double azimuth,
            double* ratio, std::size_t first_band, std::size_t last_band) const;

    // The same with Q at every sub-frequency, as panam's board_level() takes it.
    void level_exact(double source_height, double ground_distance, double azimuth,
            double* correction, std::size_t first_band, std::size_t last_band) const;

    // Per band, an upper bound on 10^(level/10) over every geometry: the most
    // the plate and ground can add.  Bands outside [first_band, last_band) are
    // left at 1.  Multiplied into the band gains, it keeps the spectral level
    // bounds (--db_cutoff and the peak search) bounds.
    std::vector<double> gain_bounds(std::size_t first_band, std::size_t last_band) const;

    double thickness() const
    {
        return thickness_;
    }
    std::string description() const;

private:
    template <bool exact, bool as_ratio = false>
    void evaluate(double source_height, double ground_distance, double azimuth,
            double* correction, std::size_t first_band, std::size_t last_band) const;

    std::vector<double> edges_;    // band edges, n_bands + 1
    std::vector<double> lagrange_; // [sub][lower edge, centre, upper edge]
    std::vector<std::complex<double>> edge_admittance_;
    std::vector<std::complex<double>> centre_admittance_;
    // Q table: [r][s][point][re, im], point 2b the lower edge of band b, 2b+1
    // its centre, 2n the top edge; s = sqrt(cos), cos the image path's
    // grazing cosine.
    std::vector<float> q_table_;
    int q_n_s_ = 0;
    int q_n_r_ = 0;
    double q_log_r_min_ = 0.0;
    double q_log_r_step_ = 1.0;
    std::size_t n_bands_ = 0;
    std::size_t n_sub_ = 0;
    std::vector<double> frequencies_; // [band][sub]
    std::vector<double> elevations_;  // degrees, ascending
    std::size_t n_azimuths_ = 0;
    double azimuth_step_ = 0.0;       // radians
    // [elevation][azimuth][band][sub][P_d re, P_d im, P_r re, P_r im]: a
    // bilinear stencil reads four contiguous runs.
    std::shared_ptr<const std::vector<double>> values_;
    std::vector<double> centres_;
    std::string ground_model_;
    double flow_resistance_ = 0.0;
    double porosity_rate_ = 0.0;
    double layer_depth_ = 0.0;
    std::vector<std::complex<double>> admittance_; // [band][sub]
    GroundModel ground_;
    double sound_speed_ = 343.2;
    double thickness_ = 0.0;
    double radius_ = 0.0;
    double mic_offset_ = 0.0;
    std::string source_;
};

// Which way the microphone is offset on its plate, in the run's local frame.
// A fixed azimuth (radians, from +x toward +y), or "outboard": away from the
// flight track, along the horizontal line from the track's nearest point to
// the observer -- the placement that suits the 2017 data (panam
// ground_plane.PLATE_MIC_OFFSET_FT).
struct PlateOrientation
{
    bool outboard = true;
    double azimuth = 0.0;
};

// The azimuth, in the plate's own frame, at which a ray arrives: the plate's tables
// put the microphone on +y, so the propagation direction (source to microphone,
// atan2 of `separation`) is turned by the microphone's actual offset direction
// `mic_azimuth` in the run's frame.  A source directly inboard of the microphone
// -- the ray travelling outward, along the offset -- is pi/2.  One function for
// both places that need it (Database::spectral_level_along and the refracted
// path in niceops.cpp), which spelled it out twice and agreed only by a comment.
inline double plate_frame_azimuth(const Vec3& separation, double mic_azimuth)
{
    return std::atan2(separation.y, separation.x) - mic_azimuth + 0.5 * 3.14159265358979323846;
}

// The microphone's offset direction for "outboard": away from the track, along the
// horizontal line from the track's nearest state to the observer.  Directly under
// the track that line has no direction, and the microphone goes to the track's
// right.  `States` is anything with a `position` (Database::SourceState).
template <typename States>
double outboard_azimuth(const Vec3& observer, const States& states)
{
    std::size_t nearest = 0;
    double nearest_distance = std::numeric_limits<double>::infinity();
    for (std::size_t it = 0; it < states.size(); it++)
    {
        const Vec3 separation = observer - states[it].position;
        const double distance = separation.x * separation.x + separation.y * separation.y;
        if (distance < nearest_distance)
        {
            nearest_distance = distance;
            nearest = it;
        }
    }
    const Vec3 separation = observer - states[nearest].position;
    if (nearest_distance > 1e-6)
        return std::atan2(separation.y, separation.x);
    const std::size_t before = nearest > 0 ? nearest - 1 : nearest;
    const std::size_t after = nearest + 1 < states.size() ? nearest + 1 : nearest;
    const Vec3 track = states[after].position - states[before].position;
    return std::atan2(track.y, track.x) - 0.5 * 3.14159265358979323846;
}

} // namespace niceops
