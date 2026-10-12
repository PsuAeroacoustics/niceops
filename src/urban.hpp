#pragma once

// Urban receivers (--urban_table): what the buildings around a microphone do to
// the sound reaching it, split by how it is evaluated.
//
//   - The direct path keeps every model the run applies today (spreading,
//     absorption, the ground) and is multiplied by the receiver's table, an
//     energy ratio per band looked up by the direction of the source as seen
//     from the microphone.  The table holds everything the explicit paths below
//     do not: where the buildings block the direct path, the diffracted field
//     over them; the specular paths of second order and higher; later, diffuse
//     facade energy and turbulent scattering.  It is defined against the
//     open-ground prediction the run already makes, so a table of 0 dB with no
//     facades is the run without buildings, to the bit.
//   - Each first-order facade reflection is a path of its own: the microphone
//     mirrored in the facade, evaluated as any straight path is -- its own
//     sphere lookup, spreading, absorption and ground term -- less the facade's
//     reflection loss.  Near the buildings a reflection leaves the rotor in a
//     direction 10-40 degrees from the direct path's, which no table keyed on
//     the arrival direction alone can carry.
//
// The table is built under plane-wave incidence (a source far away compared
// with the buildings), by reciprocity: for each direction, what reaches the
// microphone from a plane wave arriving from there, as a ratio to what open
// ground delivers.  For a vertical facade the mirror leaves heights alone, so
// the ground model's geometry stays exact on an unfolded image path.  The floor
// is not left to the ground model inside a canyon, though: its reflection
// clears the roofs after a different number of wall reflections from the
// direct ray's, and is often open where the direct ray is blocked, so the
// table traces it (canyon_ratio).
//
// This header is independent of Database, so the same geometry can serve the
// table builder (tools/urban_table.cpp) and, later, panam.

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "vec3.hpp"

namespace niceops {

// A vertical rectangular facade: the wall standing on the segment from (x0, y0)
// to (x1, y1) in the run's local frame, from height `bottom` to `top` (m, the
// run's z).  Its two faces both reflect.
struct Facade
{
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    double bottom = 0.0;
    double top = 0.0;
};

// The signed horizontal distance of `p` from the facade's plane, positive on
// the left of its base direction (x0, y0) -> (x1, y1).
double facade_side(const Facade& facade, const Vec3& p);

// `p` mirrored in the facade's plane; z is unchanged.
Vec3 mirror_in(const Facade& facade, const Vec3& p);

// A horizontal vector mirrored in the facade's plane (no translation): how a
// direction turns on reflection.
Vec3 mirror_direction(const Facade& facade, const Vec3& v);

// Whether the segment from `a` to `b` passes through the facade's rectangle:
// its ends strictly on opposite sides of the plane, and the crossing within the
// facade's length and heights (edges included).  An end on the plane does not
// cross it.
bool blocks(const Facade& facade, const Vec3& a, const Vec3& b);

// The first-order reflection of the path from `source` to `microphone` in
// facades[self].  True, with `image` the microphone mirrored in that facade,
// when the specular path exists: the two strictly on the same side of the
// facade, the reflection point within it, and neither leg blocked by another
// facade of `facades`.  The run evaluates the path as straight_ray(source,
// image), whose direction is the one leaving the source and whose length is
// the path's.
bool first_order_image(const std::vector<Facade>& facades, std::size_t self, const Vec3& source,
        const Vec3& microphone, Vec3& image);

// The facades with their walls carried down without end: their images in a
// floor, through which a floor-reflected path runs once unfolded.
std::vector<Facade> walls_below(std::vector<Facade> facades);

// Whether the first-order reflection in below[self] also reaches the
// microphone off the floor at height `floor`: the same test as
// first_order_image, for the microphone's image in the floor, among the walls
// carried below it (walls_below).  The run evaluates a reflection with the
// ground model's floor partner only when this holds, and in free field
// otherwise; a flush microphone is its own image.
bool first_order_floor_partner(const std::vector<Facade>& below, std::size_t self, const Vec3& source,
        const Vec3& microphone, double floor);

// A specular path from a source to a receiver by the image-source method: the
// facades it reflects in, in order from the source, and its unfolded length.
struct SpecularPath
{
    std::vector<std::size_t> facades;
    double length = 0.0;
};

// Every specular path of order 0 to `max_order` between `source` and
// `receiver` among `facades`, by brute force over facade sequences (no facade
// twice in a row), each checked by backtracking from the receiver: the
// reflection points within their facades and no leg blocked.  Exponential in
// the order; it is the reference the canyon model and the image paths are
// tested against, not a production path finder.
std::vector<SpecularPath> specular_paths(const std::vector<Facade>& facades, const Vec3& source,
        const Vec3& receiver, std::size_t max_order);

// An ideal street canyon in its own cross-section: two parallel walls of
// infinite length, wall A on the left of the axis and wall B on its right
// (looking along the axis's bearing), `width` apart.  The microphone stands
// `offset` from wall A and `microphone_height` above the canyon floor; the
// walls' tops are `height_a` and `height_b` above the floor.
struct IdealCanyon
{
    double width = 0.0;
    double offset = 0.0;
    double height_a = 0.0;
    double height_b = 0.0;
    double microphone_height = 0.0;
};

// A ray traced back from the microphone across the canyon's cross-section:
// leaving with horizontal component `v_c` across the canyon (positive toward
// wall B) and vertical component `u_z` (both of a unit vector), reflecting
// specularly off the walls until it passes over one.  The floor is not part of
// the trace: the ground model owns it.
struct CanyonTrace
{
    bool escaped = false;   // false: still inside after max_bounces
    int bounces = 0;
    int hits_a = 0;
    int hits_b = 0;
    // The first wall the ray reaches (0 for none, 'a' or 'b'), whether it
    // reflected there or passed over it.
    char first_wall = 0;
};

CanyonTrace trace_canyon(const IdealCanyon& canyon, double v_c, double u_z, int max_bounces = 400);

// The canyon's table entry for a plane wave arriving from the unit direction
// whose cross-canyon component is `u_c` (toward wall B positive) and vertical
// component `u_z`, per band, as an energy ratio to open ground's direct ray and
// floor reflection (THEORY 6.5):
//
//   - the rays that reach the microphone specularly.  Under plane-wave
//     incidence each of four families holds at most one: arriving from above or
//     off the floor, with the plane wave's horizontal direction (an even number
//     of wall reflections) or its mirror (an odd number).  The floor-reflected
//     rays start from the microphone's image in the floor, lower by twice its
//     height, so they clear the roofs after a different number of reflections
//     and are often open where the direct ray is not;
//   - with `diffraction`, each lit roof edge's diffracted rays (wedge.hpp, a
//     right-angle wedge) to the microphone and its images in the walls and the
//     floor, all of them summed with their phases and band-averaged pair by
//     pair, which keeps the field continuous where a ray is cut off; without,
//     the specular rays in energy;
//   - but not the first-order reflection when it is open from above, nor its
//     floor partner when both are: the run evaluates them as a path of its
//     own, with the ground model (first_order_floor_partner), and adds it in
//     energy, so their interference with the rest is left out too.
//
// `reflection_a` and `reflection_b` are the walls' energy reflection
// coefficients for bands [0, bands.size()) (third-octave centers),
// `floor_reflection` the floor's (1 for rigid pavement); `out` gets the ratio
// for each band.  A microphone above both walls, or a direction along the
// axis, gives 1.
void canyon_ratio(const IdealCanyon& canyon, double u_c, double u_z, const double* reflection_a,
        const double* reflection_b, double floor_reflection, const std::vector<double>& bands, double sound_speed,
        bool diffraction, double* out);

// The compass bearing (degrees, toward) and elevation (degrees above the
// horizontal) of a local-frame direction, for a run whose local +x points
// along compass bearing `frame_bearing`.
struct CompassDirection
{
    double azimuth = 0.0;   // [0, 360)
    double elevation = 0.0;
};

CompassDirection compass_direction(const Vec3& local, double frame_bearing);

// The lowest gain a builder writes, dB: a direction no path reaches (the
// canyon with its diffraction left out) is silent, which as a level is minus
// infinity and would read back as "no value".
inline constexpr double urban_gain_floor_db = -300.0;

// Each receiver's tables and facades (--urban_table): a netCDF file in SI
// units, built by tools/urban_table.cpp, read once per run.
//
//   azimuth     compass bearing from the microphone toward the source,
//               degrees, uniform from 0 over 360
//   elevation   degrees above the horizontal, ascending
//   band        band center, Hz; the run's bands are found among these
//   gain        [receiver][azimuth][elevation][band], dB, the energy ratio
//               above as a level; NaN where the builder has no value
//
// Values are interpolated bilinearly in dB.  A pair a receiver's table does not
// cover (an elevation outside the axis, or a NaN among the corners it would
// read) takes the run's open-ground model, without the urban terms; the run
// counts them.  Receivers are matched to the run's microphones by position.
// The local frame is the run's, so the table records the compass bearing of its
// +x axis, which the run checks against its own.
class UrbanTable
{
public:
    static constexpr int format_version = 1;

    struct Axes
    {
        std::vector<double> azimuths;   // degrees
        std::vector<double> elevations; // degrees
        std::vector<double> bands;      // Hz
    };

    struct Receiver
    {
        Vec3 microphone;               // local frame, m
        std::vector<Facade> facades;   // for the first-order images
        // [facade][band], energy reflection coefficients in [0, 1]
        std::vector<double> reflection;
    };

    // Where a run's bands sit in the table: for each run band in [first,
    // last), its column.
    struct BandMap
    {
        std::size_t first = 0;
        std::size_t last = 0;
        std::vector<std::size_t> column;
    };

    UrbanTable() = default;
    // `gains` is [receiver][azimuth][elevation][band]; throws
    // std::invalid_argument if anything is not as the header says.
    UrbanTable(Axes axes, std::vector<Receiver> receivers, std::vector<double> gains, double frame_bearing,
            std::string method, std::map<std::string, std::string> provenance = {});

    static UrbanTable load(const std::string& path);
    void save(const std::string& path) const;

    const Axes& axes() const
    {
        return axes_;
    }
    const std::vector<Receiver>& receivers() const
    {
        return receivers_;
    }
    double frame_bearing() const
    {
        return frame_bearing_;
    }
    const std::string& method() const
    {
        return method_;
    }
    const std::map<std::string, std::string>& provenance() const
    {
        return provenance_;
    }
    double gain(std::size_t receiver, std::size_t azimuth, std::size_t elevation, std::size_t band) const
    {
        return gains_[cell(receiver, azimuth, elevation) * axes_.bands.size() + band];
    }

    // The run's bands in [first_band, last_band) found in the table to 0.5%;
    // throws std::runtime_error naming the missing ones.
    BandMap match_bands(const std::vector<double>& run_bands, std::size_t first_band, std::size_t last_band) const;

    // The receiver whose microphone lies within `tolerance` (m) of
    // `microphone`, or receivers().size() for none.
    std::size_t find_receiver(const Vec3& microphone, double tolerance) const;

    // The receiver's gain for each run band in [bands.first, bands.last), dB,
    // into out[band]; false, writing nothing, where it has no value.
    bool lookup(std::size_t receiver, double azimuth_degrees, double elevation_degrees, const BandMap& bands,
            double* out) const;

    // The largest gain over every receiver and direction, per run band, as a
    // power ratio and at least 1: what the peak-search bound has to allow.
    std::vector<double> max_power_ratio(const BandMap& bands) const;

private:
    std::size_t cell(std::size_t receiver, std::size_t azimuth, std::size_t elevation) const
    {
        return (receiver * axes_.azimuths.size() + azimuth) * axes_.elevations.size() + elevation;
    }

    Axes axes_;
    std::vector<Receiver> receivers_;
    std::vector<double> gains_;
    double frame_bearing_ = 90.0;
    std::string method_;
    std::map<std::string, std::string> provenance_;
};

// The two walls of an ideal canyon as facades in the local frame, for the run's
// first-order images: each `length` long, centered on the microphone, standing
// from `floor` to `floor` plus its height.  The canyon's axis runs along
// compass bearing `axis_bearing`; the local +x axis along `frame_bearing`.
std::vector<Facade> ideal_canyon_facades(const IdealCanyon& canyon, const Vec3& microphone, double floor,
        double axis_bearing, double frame_bearing, double length);

// An ideal canyon's gains for one receiver, [azimuth][elevation][band] in dB,
// from canyon_ratio at each direction of `axes`.  `reflection_a` and
// `reflection_b` are per table band.
std::vector<double> ideal_canyon_gains(const IdealCanyon& canyon, double axis_bearing,
        const UrbanTable::Axes& axes, const std::vector<double>& reflection_a,
        const std::vector<double>& reflection_b, double floor_reflection, double sound_speed, bool diffraction,
        unsigned threads = 1);

// The wave band's weight in each table band (the low band, solved by
// src/canyon_wave): 1 for every band centered at or below `wave_below` Hz,
// falling in equal steps to 0 over the next `blend` bands, 0 above; all 0
// when wave_below is 0.
std::vector<double> wave_band_weights(const std::vector<double>& bands, double wave_below, std::size_t blend);

// Blend a receiver's wave table into its geometric one: per band, in energy,
// T = w T_wave + (1 - w) T_geometric.  The wave solution holds everything the
// street does, its first-order reflections included, so each facade's
// reflection in a band is scaled by 1 - w and the run's explicit reflections
// fade out as the wave table fades in.  `gains` and `wave` are one
// receiver's [direction][band] in dB; `reflection` is its [facade][band].
void blend_wave_band(std::vector<double>& gains, std::vector<double>& reflection, const std::vector<double>& wave,
        const std::vector<double>& weight);

} // namespace niceops
