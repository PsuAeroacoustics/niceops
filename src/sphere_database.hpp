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
}

// Bilinear interpolation of q(x, y) over the rectangle [x1, x2] x [y1, y2].
// It lived in interpolator.cpp behind the same "C" seam, but it is plain
// arithmetic with no CGAL in it, and across that translation-unit boundary it
// could not be inlined -- two out-of-line calls on every one of tens of
// millions of evaluations.
double bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y);

// Check if points "wrap around" sphere, and swap inputs to bilinear
// interpolation if necessary
double periodic_bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y);

// Find the pair of grid indices bracketing val (wrapping around the ends)
std::array<std::size_t, 2> bracket(const std::vector<double>& a, double val);

// One axis of the sphere grid: the sorted unique angles, plus the spacing when
// the axis turns out to be uniform (every shipped database's is).  Knowing the
// spacing turns bracket() into a multiply and a compare rather than a binary
// search -- eval brackets twice, so a default run does that ~190 million times,
// and the searches were around a third of the evaluation loop.
class GridAxis
{
public:
    void reset(std::vector<double> values);

    const std::vector<double>& values() const
    {
        return values_;
    }

    std::size_t size() const
    {
        return values_.size();
    }

    double operator[](std::size_t i) const
    {
        return values_[i];
    }

    // Identical in every case to bracket(values(), val).
    std::array<std::size_t, 2> bracket(double val) const;

private:
    std::vector<double> values_;
    double first_ = 0.0;
    // Reciprocal of the uniform spacing, or 0 when the axis is not uniform and
    // the binary search has to be used after all.
    double inv_step_ = 0.0;
};

double clamp01(double x);

// acos with its argument clamped into [-1, 1], so round-off just past the ends
// does not produce NaN.
double acosc(double x);

// The sphere's orientation at one instant: the unit TPP normal, the "nose"
// (velocity projected into the TPP) and starboard.  None of it depends on
// where the observer is, so it is built once per trajectory state and reused
// across the whole observer grid.
struct SphereFrame
{
    Eigen::Vector3d normal;    // n, unit TPP normal
    Eigen::Vector3d nose;      // vn, velocity projected into the TPP, unit
    Eigen::Vector3d starboard; // pn = vn x n
};

SphereFrame sphere_frame(const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal);

// Sphere angles (phi, theta in degrees) of an already-normalized propagation
// direction.
std::array<double, 2> lookup_angles(const SphereFrame& frame, const Eigen::Vector3d& direction);

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
    // The four database spheres bracketing one operating condition, and the
    // barycentric weights that blend them.  Computing this is by far the most
    // expensive part of an evaluation -- it is a CGAL Delaunay locate() walk --
    // and it depends only on (mu, gamma, CT), so it is hoisted out of any loop
    // that holds the operating condition fixed.
    struct Blend
    {
        std::array<double, 4> weights;
        std::array<std::size_t, 4> indices;
    };

    // Everything about one trajectory state that every observer shares: where
    // the aircraft is, how its sphere is oriented, and which database spheres
    // to blend.
    struct SourceState
    {
        Eigen::Vector3d position;
        SphereFrame frame;
        Blend blend;
    };

    explicit Database(const std::string& datafile_path);

    Blend blend_at(double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    SourceState source_state(const Eigen::Vector3d& source, const Eigen::Vector3d& velocity,
            const Eigen::Vector3d& normal, double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    double eval(double range, double phi, double theta, const Blend& blend) const;

    double eval(double range, double phi, double theta, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double level_at_observer(const Eigen::Vector3d& observer, const SourceState& state) const;

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
    // The dBA level and EAA excess-attenuation factor of one (phi, theta) grid
    // cell at one operating condition.  Kept adjacent because eval always reads
    // the pair, so one cache line serves both.
    struct Sample
    {
        double dBA;
        double EAA;
    };

    // Blend the four bracketing operating conditions at one grid cell.
    Sample blend_cell(std::size_t i_phi, std::size_t i_theta, const Blend& blend) const;

    GridAxis phis;
    GridAxis thetas;
    // Flat (phi, theta, condition) grid, indexed
    // (i_phi * thetas.size() + i_theta) * n_conditions + condition.  This was a
    // vector<vector<vector<double>>>, which put two dependent pointer loads in
    // front of every one of the 32 reads an eval makes -- the largest single
    // cost of a run once the CGAL locate had been hoisted out of the loop.
    std::vector<Sample> grid;
    std::size_t n_conditions = 0;
    double sphere_radius_ = 0.0;
    double advance_ratio_min = 0.0, advance_ratio_max = 0.0;
    double flight_path_angle_min = 0.0, flight_path_angle_max = 0.0;
    double thrust_coefficient_min = 0.0, thrust_coefficient_max = 0.0;
    double max_level_ = 0.0;
};

} // namespace niceops
