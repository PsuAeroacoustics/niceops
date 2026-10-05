#pragma once

#include "atmosphere.hpp"
#include "util/constants.hpp"

// Port of trajectory.d.  Parses the 27-column NASA/FAA AC state CSV into
// metric State records, computing non-dimensional flight parameters:
// CT = LF*W/(rho*V_tip^2*pi*R^2), mu = v/V_tip, and flight path angles fpa
// (vertical) vs fpa_e (relative to the TPP normal).  adapt() spline-smooths
// the trajectory on a uniform dt grid; watts_trim() crops it to the region
// of interest.
//
// It also reads a simpler, labeled format in WGS84 coordinates -- latitude,
// longitude and altitude against time -- which is what a mission planner has
// rather than a flight-test state file.  That path projects the positions onto
// the run's local frame (see geodetic.hpp) and derives the velocity and load
// factor the flight parameters need from the path itself.

#include <array>
#include <functional>
#include <limits>
#include <limits>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "geodetic.hpp"

namespace niceops {

constexpr double g = 9.80665; // m/s^2
constexpr double mps2g = 1.0 / g;
// ft2m lives in util/constants.hpp, beside the unit system that uses it.
constexpr double fps2g = mps2g * ft2m;

struct State
{
    double t;
    // Seconds since UTC midnight, straight from the AC format's first column.
    // Only the day/night classification behind the DNL metric uses it, so it is
    // deliberately absent from write_state's column list.
    double utcsec = 0.0;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double nx;
    double ny;
    double nz;
    double LF;
    double CT;
    double mu;
    double fpa_e;
    double fpa;
    // The horizontal wind at the state, m/s in the local frame.  Zero unless
    // Trajectory::relative_to_air() set it, which leaves mu, both path angles
    // and the sphere's nose following air_velocity() rather than the ground
    // velocity.  Last, so that brace-initializing the fields above still works.
    double wx = 0.0;
    double wy = 0.0;
    // Compass heading of the nose, degrees from true north: the AC format's
    // column 9 (data[8]), or a static operation's heading column.  NaN where the
    // file carries none (the geodetic formats).  Only the sphere's orientation
    // reads it: at hover, and at every speed under --nose heading
    // (nose_direction).
    double heading = std::numeric_limits<double>::quiet_NaN();
};

// The velocity the rotor sees: the ground velocity less the wind, the vertical
// wind taken as zero.  The ground velocity itself while there is no wind.
inline std::array<double, 3> air_velocity(const State& s)
{
    return {s.vx - s.wx, s.vy - s.wy, s.vz};
}

// The ground velocity, whatever wind relative_to_air() has set.
inline std::array<double, 3> ground_velocity(const State& s)
{
    return {s.vx, s.vy, s.vz};
}

// Where the sphere frame's nose comes from (--nose).  `velocity` is the velocity
// the flight conditions take (air_velocity: the ground velocity, or under
// --speed_reference air the air velocity); `ground_velocity` is the ground
// velocity whatever --speed_reference says, so the nose and mu can be referred
// separately; `heading` is the state's compass heading at every speed.  The
// first two turn to the heading at hover (below); the sphere frame projects any
// of them into the tip-path plane (sphere_frame).
enum class NoseSource
{
    velocity,
    ground_velocity,
    heading,
};

// "velocity", "ground_velocity" or "heading"; nullopt for anything else.
std::optional<NoseSource> parse_nose_source(const std::string& name);
const char* to_string(NoseSource source);

// The sphere's nose turns from the heading to the velocity as the horizontal
// speed rises through this band, on a smoothstep weight.  Below the lower edge
// (3 kt) the velocity's direction is drift, not where the aircraft points: a
// hover's GPS velocity is a few tenths of a knot in any direction.  Above the
// upper edge (8 kt) the nose is the velocity.  Shipped tracks never fly under
// 10 kt.  Before this, a hard switch at 5 kt stepped the level 2.2-4.4 dB for a
// nose rotated 30-90 degrees between the two.
constexpr double heading_nose_speed = 3.0 * 1852.0 / 3600.0; // m/s
constexpr double velocity_nose_speed = 8.0 * 1852.0 / 3600.0; // m/s

// The smoothstep weight the nose gives the velocity at horizontal speed
// `speed`, m/s: 0 at or below heading_nose_speed, 1 at or above
// velocity_nose_speed.
double nose_velocity_weight(double speed);

// The local frame's angle atan2(y, x) of compass heading h is offset - h; this
// is that offset, radians.  From `frame_bearing` (degrees, NaN when the run
// knows no frame), the bearing of +x, where the offset is the bearing itself.
// Only when there is no frame is it inferred from the track, as the circular
// mean of atan2(vy, vx) + heading over its states faster than 20 kt, when it
// has 50 such, the velocity being air_velocity, or the ground velocity when
// `ground` is set.  That fit folds the mean crab -- the heading's angle off the
// track, a few degrees of INS bias and, in a crosswind, the wind correction --
// into the offset, so a hover heading read through it is off by that much;
// it is the fallback, never preferred to a known frame.  nullopt when neither
// is available or the file carries no heading.
std::optional<double> heading_offset(const std::vector<State>& states, double frame_bearing,
        bool ground = false);

// The direction the sphere's nose takes for a state.  For `velocity` and
// `ground_velocity`: that velocity above velocity_nose_speed, the heading
// (horizontal) below heading_nose_speed, and between them the horizontal angle
// turned from one to the other.  For `heading`: the heading's horizontal
// direction at every speed, carrying the air velocity's climb (vz over the
// horizontal speed) by nose_velocity_weight, so that a heading along the track
// gives exactly the direction `velocity` does, climbs and descents included.
// `heading` needs an offset and a finite heading (the caller checks; without
// them this returns the air velocity, as the others do).
std::array<double, 3> nose_direction(const State& state, std::optional<double> offset,
        NoseSource source = NoseSource::velocity);

inline std::ostream& operator<<(std::ostream& os, const State& s)
{
    return os << "State(t: " << s.t << ", x: " << s.x << ", y: " << s.y << ", z: " << s.z << ", vx: "
        << s.vx << ", vy: " << s.vy << ", vz: " << s.vz << ", nx: " << s.nx << ", ny: " << s.ny
        << ", nz: " << s.nz << ", LF: " << s.LF << ", CT: " << s.CT << ", mu: " << s.mu << ", fpa_e: "
        << s.fpa_e << ", fpa: " << s.fpa << ")";
}

struct dimensional_data
{
    double tip_speed = 0.0; // m/s
    double radius = 0.0; // m
    double weight = 0.0; // kg (converted to Newtons)
    double rho = 0.0; // kg/m^3
    // Atmospheric state for the spectral propagation path, from the same file's
    // "Atmosphere" section.  Unused by CT and mu, which need only rho.
    Atmosphere atmosphere;
    // Which of the atmosphere's values the file actually gave, so the run can say
    // which ones it filled in: humidity above all defaults to 20 percent, and
    // that moves the spectral footprint by 1 to 3 dB.
    bool temperature_given = false;
    bool pressure_given = false;
    bool humidity_given = false;
};

inline std::ostream& operator<<(std::ostream& os, const dimensional_data& v)
{
    return os << "dimensional_data(tip_speed: " << v.tip_speed << ", radius: " << v.radius << ", weight: "
        << v.weight << ", rho: " << v.rho << ")";
}

// How a geodetic trajectory file is turned into the local Cartesian states the
// rest of the tool works in.
struct GeodeticInput
{
    // The frame every geodetic position is projected onto.  One run shares one
    // frame across all of its trajectories -- the observers are common to them
    // all -- so it is settled before any file is read and handed in here.
    // Unset lets a lone trajectory center its own frame on its first position.
    std::optional<LocalFrame> frame;
    // Elevation of the observer plane above the WGS84 datum, in METERS -- the
    // command line's --ground_elevation is in the user's units and is converted
    // with every other length before it reaches here.  An
    // altitude column is height above mean sea level, while z is height above
    // the ground the observers stand on: in the shipped AC tracks the two
    // differ by a constant 54.2651 ft, the terrain under that particular
    // flight.  Leaving it at zero reads altitudes as height above the
    // observers, which is what a planner working in AGL already means.
    double ground_elevation = 0.0; // meters
    // Seconds since UTC midnight at t = 0.  Only the day-night window reads it,
    // and a file carrying its own utcsec column supplies it instead.  NaN means
    // neither did: the states then have no clock and -D refuses them.  It
    // defaulted to 0, which dates the flight at UTC midnight -- night in every
    // timezone but the one a few hours west -- so a clockless track got the 10 dB
    // night penalty silently.
    double start_utc = std::numeric_limits<double>::quiet_NaN();
    // Resampling step in seconds, used when the samples are unevenly spaced or
    // coarser than this.  It never coarsens -- a 50 Hz file keeps its own
    // sample times -- because the sound exposure integral is a Riemann sum over
    // whatever step the states end up with.
    double resample_dt = 0.5;
    // How far, in meters RMS, the fitted path may sit from the given positions.
    // Zero interpolates, which is what waypoints want; measured positions want
    // enough slack to keep the second derivative out of the noise, since that
    // second derivative is the load factor and so CT.
    double path_smoothing = 0.0;
};

// The first WGS84 position in a labeled geodetic trajectory file, or nullopt
// when `filename` is not one.  Reads the header and one row, so a run can
// settle its frame before committing to loading anything.
//
// The unit is not defaulted: which altitude spellings are positioned depends on
// it, and a caller that forgot it silently looked for the feet spellings only.
std::optional<Geodetic> first_geodetic_position(const std::string& filename, LengthUnit unit);

// Read the vehicle JSON into its dimensional parameters.  Exported because the
// sphere export needs the tip speed to turn an advance ratio back into an
// airspeed, without building a whole Trajectory to get at it.
dimensional_data configure_vehicle(const std::string& filename);

// Derive the non-dimensional flight parameters (LF, CT, mu, fpa, fpa_e) from a
// state's raw kinematics.  Expects nx/ny/nz to hold the un-normalized
// load-factor components on entry, and normalizes them in place.  Exported so
// the unit tests can reach the hover case directly.
void derive_flight_parameters(State& now, const dimensional_data& vehicle);

// Negate vz (and recompute the path angles) when a track's vz runs opposite to
// dz/dt: the 2017 NASA/FAA AC files store vertical speed positive down (bar
// EC130B4 day 298), the 2019 ones positive up.  Leaves a track with too little vertical motion to
// judge untouched.  Returns true if it corrected the sign.  Exported for tests.
bool correct_vertical_speed_sign(std::vector<State>& states, const std::string& filename);

class Trajectory
{
public:
    // `unit` is the length unit the *file* is written in -- feet unless --units
    // says otherwise.  It does not reach write_state, which is metric; see
    // LengthUnit in util/constants.hpp.
    Trajectory(const std::string& filename, const std::string& vehicle_filename,
            const GeodeticInput& geodetic = {}, LengthUnit unit = feet);

    // The frame a geodetic file was projected through, or nullopt when the
    // trajectory was already Cartesian.  Lets a run adopt the frame its first
    // geodetic trajectory settled on, and invert it for geodetic output.
    const std::optional<LocalFrame>& frame() const
    {
        return frame_;
    }

    void write_state(const std::string& state_file) const;

    // Resample onto a uniform dt grid, smoothing with a natural cubic spline.
    // Pass a non-negative `tolerance` to bound how far the smoothed path may sit
    // from the measured one (weighted sum of squared residuals, as MATLAB's
    // spaps does); leave it negative to use `smoothing` as the penalty weight
    // directly.
    void adapt(double dt, double smoothing = 0.3, double tolerance = -1.0);

    // D's Get!(varname).get(): collect one field from every state
    std::vector<double> column(const std::string& varname) const;

    // Take mu, both flight path angles and the sphere's nose from the velocity
    // through the air rather than over the ground: `wind(z)` gives the
    // horizontal wind at height z, m/s in the local frame.  The load factor,
    // and so CT, stays inertial, because the rotor's thrust balances the
    // aircraft's acceleration over the ground whatever the air is doing.
    // Returns the largest wind speed met, for the log.
    double relative_to_air(const std::function<std::array<double, 2>(double)>& wind);

    // Crop from this x, in METERS.  The caller converts; the command line's
    // --trim is in the user's units and is converted at the same seam as
    // everything else.
    void watts_trim(double xstart_meters);

    std::vector<State> states;
    dimensional_data vehicle;

private:
    std::optional<LocalFrame> frame_;
};

} // namespace niceops
