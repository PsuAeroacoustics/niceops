#pragma once

// Flight procedures flown into tracks (niceops_procedure): a start state and a
// list of segments -- level, climb or descent at a path angle, a change of
// speed, a coordinated turn, a hover -- integrated into a WGS84 track with
// velocity and heading columns, which niceops reads as it is.  Sparse
// waypoints make the load factor an artifact of the interpolant (FIDELITY,
// "Trajectory smoothing"); a flown procedure carries its own.
//
// Transitions are flown, not stepped: the path angle changes at a limited
// load factor, the bank at a limited roll rate and the speed at the segment's
// acceleration, so joining segments puts no spike in the load factor.  Calm
// air: the heading is the course, and the ground speed the airspeed.

#include <limits>
#include <string>
#include <vector>

#include "geodetic.hpp"
#include "util/constants.hpp"

namespace niceops {

struct ProcedureSegment
{
    enum class Kind
    {
        level,
        climb,  // or descent: path_angle signed, up positive
        speed,  // accelerate or decelerate to target_speed
        turn,   // coordinated, at bank, through turn
        hover,  // stationary for duration
    };
    Kind kind = Kind::level;
    double duration = std::numeric_limits<double>::quiet_NaN();       // s
    double distance = std::numeric_limits<double>::quiet_NaN();       // m, level only
    double until_altitude = std::numeric_limits<double>::quiet_NaN(); // m, climb only
    double path_angle = 0.0;                                          // rad
    double target_speed = std::numeric_limits<double>::quiet_NaN();   // m/s
    double acceleration = std::numeric_limits<double>::quiet_NaN();   // m/s^2, magnitude
    double bank = std::numeric_limits<double>::quiet_NaN();           // rad, positive turns right
    double turn = std::numeric_limits<double>::quiet_NaN();           // rad, magnitude
};

struct Procedure
{
    // The start: a WGS84 position, an altitude on the track's own datum (what
    // the alt column holds; niceops takes --ground_elevation off it), a
    // ground speed and a course.
    Geodetic start {0.0, 0.0, 0.0};
    double altitude = 0.0; // m
    double speed = 0.0;    // m/s
    double course = 0.0;   // rad, compass
    // The ground's altitude on the track's datum: a state below it is refused.
    double ground_altitude = 0.0; // m
    double step = 0.2;     // s between written states
    // How hard transitions are flown.
    double pull = 0.1;                         // g above or below 1 g, changing the path angle
    double roll_rate = 10.0 * 3.14159265358979323846 / 180.0; // rad/s
    double acceleration = 0.5;                 // m/s^2, when a speed segment gives none
    std::vector<ProcedureSegment> segments;
};

struct ProcedureState
{
    double t = 0.0;
    Geodetic position {0.0, 0.0, 0.0};
    double altitude = 0.0;       // m, the track's datum
    double east = 0.0, north = 0.0, up = 0.0; // velocity, m/s
    double heading = 0.0;        // deg, compass
};

// Reads a procedure (JSON, see the manual's Companion tools).  Every length
// and speed names its unit in its key (alt_ft, speed_kt, ...); a bare number
// is refused, as on niceops's command line.  Throws input_error naming the
// problem.
Procedure read_procedure(const std::string& path);

// Flies the procedure, returning a state every `step` seconds.  Throws
// input_error for a segment that cannot be flown (a hover not at rest, a turn
// at no speed, a climb that never reaches its altitude, a segment that takes
// the aircraft below the ground).
std::vector<ProcedureState> fly(const Procedure& procedure);

// Writes the WGS84 track niceops reads: time, lat, lon, alt, vE, vN, vUp,
// heading, lengths in `unit`.
void write_procedure_track(const std::string& path, const std::vector<ProcedureState>& states, LengthUnit unit);

} // namespace niceops
