#pragma once

// Port of sphere_database.d.  Loads a netCDF noise database (.nod): one group
// per operating condition (mu advance ratio, gamma flight path angle, CT
// thrust coefficient), each holding scattered (phi, theta) channels of dBA
// level and EAA excess-attenuation factor.  eval() computes ground level:
// geometry -> sphere angles (lookup_angles), barycentric interpolation between
// operating conditions, periodic bilinear interpolation on the sphere grid,
// then spherical spreading + excess attenuation (EAA x range) + a fixed 4 dB
// ground effect.  Only same_grid == true databases are supported
// (interpolation path is unimplemented).

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace niceops {

// Declare "C" functions from interpolator
extern "C"
{
int set_interpolator3(std::size_t n, const double* x1, const double* x2, const double* x3);
void barycentric(double x1, double x2, double x3, double* weights, std::size_t* indicies);
double bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y);
}

// Check if points "wrap around" sphere, and swap inputs to bilinear
// interpolation if necessary
double periodic_bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y);

// Find the pair of grid indices bracketing val (wrapping around the ends)
std::array<std::size_t, 2> bracket(const std::vector<double>& a, double val);

double clamp01(double x);

// acos with its argument clamped into [-1, 1], so round-off just past the ends
// does not produce NaN.
double acosc(double x);

// Sphere angles (phi, theta in degrees) of the propagation direction
std::array<double, 2> lookup_angles(const Eigen::Vector3d& source,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal,
        const Eigen::Vector3d& observer);

// Repair the azimuth labels of a sphere written by a pre-versioned builder
// (database_version < 1), whose upper surface was completed with `-phi` instead
// of `180 - phi`.  Rewrites the mirrored half of `phis_scattered` in place and
// returns true if the sphere needed repair.  Leaves anything that does not match
// the broken shape untouched.
bool repair_legacy_hover_directivity(std::vector<double>& phis_scattered);

// Faithful port of D's std.math.approxEqual(value, reference, maxRelDiff,
// maxAbsDiff): relative term is against `reference`, defaults 1e-2 / 1e-5.
bool approx_equal(double value, double reference, double max_rel_diff = 1e-2,
        double max_abs_diff = 1e-5);

class Database
{
public:
    explicit Database(const std::string& datafile_path);

    double eval(double range, double phi, double theta, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double level_at_observer(const Eigen::Vector3d& observer,
            const Eigen::Vector3d& source, const Eigen::Vector3d& velocity,
            const Eigen::Vector3d& normal, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double level_along_bearing(double range, double bearing, double altitude,
            const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal,
            double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    double zero_range_level(double altitude, const Eigen::Vector3d& velocity,
            const Eigen::Vector3d& normal, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double find_range(double contour_level, double bearing, double altitude,
            const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal,
            double advance_ratio, double flight_path_angle,
            double thrust_coefficient, double threshold, bool zero_range_check) const;

    void self_test(const std::string& datafile_path);

    double sphere_radius() const
    {
        return sphere_radius_;
    }


private:
    std::vector<double> phis;
    std::vector<double> thetas;
    std::vector<std::vector<std::vector<double>>> dBA;
    std::vector<std::vector<std::vector<double>>> EAA;
    double sphere_radius_ = 0.0;
    double advance_ratio_min = 0.0, advance_ratio_max = 0.0;
    double flight_path_angle_min = 0.0, flight_path_angle_max = 0.0;
    double thrust_coefficient_min = 0.0, thrust_coefficient_max = 0.0;
    double max_level_ = 0.0;
};

} // namespace niceops
