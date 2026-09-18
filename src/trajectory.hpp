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

// Derive the non-dimensional flight parameters (LF, CT, mu, fpa, fpa_e) from a
// state's raw kinematics.  Expects nx/ny/nz to hold the un-normalised
// load-factor components on entry, and normalises them in place.  Exported so
// the unit tests can reach the hover case directly.
void derive_flight_parameters(State& now, const dimensional_data& vehicle);

class Trajectory
{
public:
    Trajectory(const std::string& filename, const std::string& vehicle_filename);

    void write_state(const std::string& state_file) const;

    // Resample onto a uniform dt grid, smoothing with a natural cubic spline.
    // Pass a non-negative `tolerance` to bound how far the smoothed path may sit
    // from the measured one (weighted sum of squared residuals, as MATLAB's
    // spaps does); leave it negative to use `smoothing` as the penalty weight
    // directly.
    void adapt(double dt, double smoothing = 0.3, double tolerance = -1.0);

    // D's Get!(varname).get(): collect one field from every state
    std::vector<double> column(const std::string& varname) const;

    void watts_trim(double xstart_ft);

    std::vector<State> states;
    dimensional_data vehicle;
};

} // namespace niceops
