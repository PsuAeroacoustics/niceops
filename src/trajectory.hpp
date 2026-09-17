#pragma once

// Port of trajectory.d.  Parses the 27-column NASA/FAA AC state CSV into
// metric State records, computing non-dimensional flight parameters:
// CT = LF*W/(rho*V_tip^2*pi*R^2), mu = v/V_tip, and flight path angles fpa
// (vertical) vs fpa_e (relative to the TPP normal).  adapt() spline-smooths
// the trajectory on a uniform dt grid; watts_trim() crops it to the region
// of interest.

#include <ostream>
#include <string>
#include <vector>

namespace niceops {

constexpr double g = 9.80665; // m/s^2
constexpr double mps2g = 1.0 / g;
constexpr double ft2m = 0.3048;
constexpr double fps2g = mps2g * ft2m;

struct State
{
    double t;
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
};

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
};

inline std::ostream& operator<<(std::ostream& os, const dimensional_data& v)
{
    return os << "dimensional_data(tip_speed: " << v.tip_speed << ", radius: " << v.radius << ", weight: "
        << v.weight << ", rho: " << v.rho << ")";
}

class Trajectory
{
public:
    Trajectory(const std::string& filename, const std::string& vehicle_filename);

    void write_state(const std::string& state_file) const;

    void adapt(double dt, int degree = 3, double alpha = 0.3);

    // D's Get!(varname).get(): collect one field from every state
    std::vector<double> column(const std::string& varname) const;

    void watts_trim(double xstart_ft);

    // NOTE: D zero-initializes this, so the isNaN guard in the constructor is
    // dead code (it never gets set from data[0]).  Kept for faithfulness.
    double utcstart = 0.0;
    std::vector<State> states;
    dimensional_data vehicle;
};

} // namespace niceops
