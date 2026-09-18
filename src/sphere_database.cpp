#include "sphere_database.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

#include <netcdf.h>

#include "util/logger.hpp"
#include "util/rootfind.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// M_PI is not part of the C++ standard and disappears from glibc's math.h
// under -std=c++17, so define it locally.
constexpr double PI = 3.14159265358979323846;

// RAII replacement for D's scope(exit) nc_close.  Note: intentionally does not
// run on fatal() paths (fatal aborts the process).
struct NcGuard
{
    int id;
    explicit NcGuard(int i) : id(i) {}
    ~NcGuard()
    {
        nc_close(id);
    }
};

std::vector<int> get_group_ids(int ncid)
{
    // Get number of groups (spheres) in database
    int numgroups;
    if (NC_NOERR != nc_inq_grps(ncid, &numgroups, nullptr))
    {
        fatal("Can't inquire netCDF groups.");
    }
    // Get the group ids
    std::vector<int> group_ids(numgroups);
    if (NC_NOERR != nc_inq_grps(ncid, &numgroups, group_ids.data()))
    {
        fatal("Cannot get number of netCDF groups.");
    }
    info("There are ", numgroups, " groups.");
    return group_ids;
}

std::size_t get_dim_length(const std::string& dimname, int group_id)
{
    // Get number of channels
    int dim_id;
    if (NC_NOERR != nc_inq_dimid(group_id, dimname.c_str(), &dim_id))
    {
        fatal("Unable to find ", dimname, "  dimension.");
    }
    std::size_t dim_length;
    if (NC_NOERR != nc_inq_dimlen(group_id, dim_id, &dim_length))
    {
        fatal("Unable to find length of ", dimname, "  dimension.");
    }
    return dim_length;
}

template <typename T>
int get_scalar(int ncid, const std::string& variable_name, T* value)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.c_str(), &varid);
    if (NC_NOERR != status)
    {
        log("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        status = nc_get_var(ncid, varid, value);
        if (NC_NOERR != status)
        {
            fatal("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
        }
    }
    return status;
}

// Read an integer-valued scalar, converting from whatever numeric type the
// file actually uses.  nc_get_var() reads the EXTERNAL type with no conversion,
// so reading a byte or short variable through it into an int reads the wrong
// width; nc_get_var_int() converts.  This matters because "same_grid" is stored
// as short in older databases and as byte in newer ones.  Returns `fallback`
// when the variable is absent.
int get_int_scalar(int ncid, const std::string& variable_name, int fallback)
{
    int varid;
    if (NC_NOERR != nc_inq_varid(ncid, variable_name.c_str(), &varid))
    {
        return fallback;
    }
    int value = fallback;
    if (NC_NOERR != nc_get_var_int(ncid, varid, &value))
    {
        fatal("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
    }
    return value;
}

template <typename T>
void require_scalar(int ncid, const std::string& variable_name, T* value)
{
    if (NC_NOERR != get_scalar(ncid, variable_name, value))
    {
        fatal(variable_name, " required but not found.");
    }
}

template <typename T>
int get_array(int ncid, const std::string& variable_name, std::size_t countp, std::vector<T>& values)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.c_str(), &varid);
    if (NC_NOERR != status)
    {
        info("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        std::size_t startp = 0;
        values.resize(countp);
        status = nc_get_vara(ncid, varid, &startp, &countp, values.data());
        if (NC_NOERR != status)
        {
            fatal("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
        }
    }
    return status;
}

template <typename T>
void require_array(int ncid, const std::string& variable_name, std::size_t countp, std::vector<T>& values)
{
    if (NC_NOERR != get_array(ncid, variable_name, countp, values))
    {
        fatal(variable_name, " required but not found.");
    }
}

double barycentric_interp(const std::vector<std::vector<std::vector<double>>>& metric, std::size_t i_phi,
        std::size_t i_theta, const std::array<double, 4>& weights, const std::array<std::size_t, 4>& indicies)
{
    return weights[0] * metric[i_phi][i_theta][indicies[0]] + weights[1] * metric[i_phi][i_theta][indicies[1]]
        + weights[2] * metric[i_phi][i_theta][indicies[2]] + weights[3] * metric[i_phi][i_theta][indicies[3]];
}

} // namespace

// These are declared in sphere_database.hpp inside namespace niceops, so they
// must be defined there too -- defining them at global scope instead leaves
// niceops::acosc & co. undefined (breaking the unit tests at link time) and
// makes unqualified calls here ambiguous against the using-directive above.
namespace niceops {

// Repair the azimuth labels of a sphere written by a pre-versioned builder.
//
// Source spheres cover only the lower half of the roll circle, phi in [-90, 90]
// with phi = 0 straight down; the upper half is the reflection through the
// horizontal plane, phi -> 180 - phi.  Older builders completed the sphere with
// `concatenate((phi, -phi))` instead.  Negation is a left/right flip *within*
// the lower half and, because the source range is symmetric about zero, maps
// [-90, 90] onto itself -- so the second half of the channels is a copy of the
// lower hemisphere wearing lower-hemisphere labels, and the upper surface is
// missing.  Forward-flight spheres in these files were completed correctly; only
// the hover spheres are affected.
//
// The acoustic data is fine -- only the labels are wrong -- so this is
// recoverable.  For a mirrored channel the file stores phi = -phi_source, hence
//     phi_correct = wrap(180 - phi_source) = wrap(180 + phi_stored).
// Returns true if the sphere was repaired.
bool repair_legacy_hover_directivity(std::vector<double>& phis_scattered)
{
    const std::size_t n = phis_scattered.size();
    // A correctly completed sphere spans the full circle.  One confined to
    // [-90, 90] never received its upper surface.
    for (auto phi : phis_scattered)
    {
        if (phi < -90.0 - 1e-9 || phi > 90.0 + 1e-9)
        {
            return false;
        }
    }
    // The broken form emits [source, -source], so the halves must mirror.
    // Refuse to touch anything that does not match that shape.
    if (n < 2 || n % 2 != 0)
    {
        log("Sphere azimuths confined to [-90, 90] but channel count is odd; "
            "leaving directivity untouched.");
        return false;
    }
    const std::size_t half = n / 2;
    for (std::size_t c = 0; c < half; c++)
    {
        if (std::fabs(phis_scattered[half + c] + phis_scattered[c]) > 1e-9)
        {
            log("Sphere azimuths confined to [-90, 90] but halves are not "
                "mirrored; leaving directivity untouched.");
            return false;
        }
    }
    for (std::size_t c = half; c < n; c++)
    {
        double corrected = 180.0 + phis_scattered[c];
        // Wrap into [-180, 180) to match the convention used elsewhere.
        corrected = std::fmod(corrected + 180.0, 360.0);
        if (corrected < 0.0)
        {
            corrected += 360.0;
        }
        phis_scattered[c] = corrected - 180.0;
    }
    return true;
}

double periodic_bilinear_interp(double q11, double q12, double q21, double q22, double x1, double x2,
        double y1, double y2, double x, double y)
{
    // Check if points "wrap around" sphere, and swap inputs to bilinear
    // interpolation if necessary
    if (x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(q11, q21);
        std::swap(q12, q22);
    }
    if (y1 > y2)
    {
        std::swap(y1, y2);
        std::swap(q11, q12);
        std::swap(q21, q22);
    }
    return bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, x, y);
}

std::array<std::size_t, 2> bracket(const std::vector<double>& a, double val)
{
    // Grab the range up to target value (first element strictly greater than val)
    auto it = std::upper_bound(a.begin(), a.end(), val);
    // Index is one less the length of this range
    long index = static_cast<long>(it - a.begin()) - 1;

    // Wrap around if the range of a is exceeded at either end.
    if (index < 0 || index >= static_cast<long>(a.size()) - 1)
        return {a.size() - 1, 0};
    return {static_cast<std::size_t>(index), static_cast<std::size_t>(index) + 1};
}

double clamp01(double x)
{
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

double acosc(double x)
{
    if (x > 1.0)
    {
        x = 1.0;
    }
    else if (x < -1.0)
    {
        x = -1.0;
    }
    return std::acos(x);
}

std::array<double, 2> lookup_angles(const Eigen::Vector3d& source, const Eigen::Vector3d& velocity,
        const Eigen::Vector3d& normal, const Eigen::Vector3d& observer)
{
    Eigen::Vector3d r = (observer - source).normalized();
    Eigen::Vector3d n = normal.normalized();
    // Project velocity vector onto TPP to get a "nose" vector
    Eigen::Vector3d v = velocity.normalized();
    Eigen::Vector3d vn = v - v.dot(n) * n;
    vn.normalize();
    // Sphere theta angle is relative to projection of velocity onto TPP
    double theta = acosc(vn.dot(r));
    // Define starboard as the cross product of the "nose" and the TPP normal
    Eigen::Vector3d pn = vn.cross(n);
    // Project the propagation vector onto the plane normal to the "nose"
    Eigen::Vector3d rp = r - r.dot(vn) * vn;
    rp.normalize();
    // Then, we can get sphere phi by computing the angle of the propagation
    // vector in the plane
    double phi = acosc(std::fabs(n.dot(rp)));
    // Sign is still uncertain, so check to see where we are in relation to starboard
    double phi_p = acosc(pn.dot(rp) / pn.norm());
    // If we're more than 90 degrees away, phi is defined negative
    if (phi_p > 0.5 * PI)
    {
        phi = -phi;
    }
    // Convert to degrees and return phi and theta
    phi *= 180.0 / PI;
    theta *= 180.0 / PI;
    return {phi, theta};
}

bool approx_equal(double value, double reference, double max_rel_diff, double max_abs_diff)
{
    // Faithful port of D's std.math.approxEqual: the relative term is taken
    // against `reference` (the second argument) rather than symmetrically, and
    // the defaults are maxRelDiff=1e-2 / maxAbsDiff=1e-5.  Both matter -- the
    // database grid fill in Database::Database matches scattered channel
    // angles with this, and a tighter tolerance silently drops channels,
    // lowering dba_max and corrupting every interpolated level downstream.
    if (reference == 0.0)
        return std::fabs(value) <= max_abs_diff;
    return std::fabs((value - reference) / reference) <= max_rel_diff
            || std::fabs(value - reference) <= max_abs_diff;
}

Database::Database(const std::string& datafile_path)
{
    // Open netCDF file
    int file_ncid;
    info("Opening netCDF file: ", datafile_path);
    if (NC_NOERR != nc_open(datafile_path.c_str(), NC_NOWRITE, &file_ncid))
    {
        fatal("Can't open netCDF file at ", datafile_path);
    }
    NcGuard guard(file_ncid);

    // Database format version.  Absent (or 0) means the file predates the
    // versioned builder, and its hover spheres may carry the broken directivity
    // described in repair_legacy_hover_directivity().
    const int database_version = get_int_scalar(file_ncid, "database_version", 0);
    info("Database format version: ", database_version);

    // Check for "same_grid" flag to determine if interpolation is required.
    // Missing is fatal rather than silently undefined -- it used to be read
    // into an uninitialised int, so an older or hand-made file could take
    // either branch at random.
    const int same_grid_flag = get_int_scalar(file_ncid, "same_grid", -1);
    if (same_grid_flag < 0)
    {
        fatal("Database error: required variable same_grid is missing.");
    }
    const bool same_grid = same_grid_flag != 0;
    info("same_grid is set to: ", same_grid);

    // TODO: Enable interpolation when same_grid == false
    if (!same_grid)
        fatal("same_grid == false not supported");

    auto group_ids = get_group_ids(file_ncid);

    // Allocate storage for operating conditions
    std::vector<double> advance_ratios, flight_path_angles, thrust_coefficients;
    advance_ratios.reserve(group_ids.size());
    flight_path_angles.reserve(group_ids.size());
    thrust_coefficients.reserve(group_ids.size());

    // Local operating parameters
    double advance_ratio, flight_path_angle, thrust_coefficient, radius, scale, dimensional_radius;
    std::vector<double> phis_scattered, thetas_scattered, dBAs_scattered, EAAs_scattered;
    std::size_t nchannels0;
    auto dba_max = -std::numeric_limits<double>::infinity();
    std::size_t repaired_spheres = 0;
    // Cycle through sphere groups and load data
    for (std::size_t i = 0; i < group_ids.size(); i++)
    {
        auto group_id = group_ids[i];
        require_scalar(group_id, "advance_ratio", &advance_ratio);
        require_scalar(group_id, "flight_path_angle", &flight_path_angle);
        require_scalar(group_id, "thrust_coefficient", &thrust_coefficient);
        require_scalar(group_id, "radius", &radius);
        require_scalar(group_id, "rotor_scale", &scale);
        advance_ratios.push_back(advance_ratio);
        flight_path_angles.push_back(flight_path_angle);
        thrust_coefficients.push_back(thrust_coefficient);
        trace("Sphere ", group_id, " loaded with \t μ = ", advance_ratios[i],
                ",\t γ = ", flight_path_angles[i], ",\t CT = ", thrust_coefficients[i],
                ",\t scale = ", scale, ",\t radius = ", radius);
        if (i == 0)
        {
            dimensional_radius = radius * scale;
        }
        else if (!approx_equal(radius * scale, dimensional_radius))
        {
            fatal("Database error: sphere radius*scale not constant.");
        }
        auto nchannels = get_dim_length("channels", group_id);
        if (i == 0)
        {
            // Record number of channels to check for changes
            nchannels0 = nchannels;
            info("First sphere contains ", nchannels, " channels.");
        }
        else if (nchannels != nchannels0)
        {
            fatal("Database error: number of channels is not constant.");
        }
        // Read THIS sphere's own channel angles.  They are not identical across
        // spheres: in a pre-versioned database the hover spheres label their
        // mirrored half differently from the forward-flight ones, so reusing
        // sphere 0's arrays files their levels into the wrong grid cells.
        phis_scattered.resize(nchannels);
        thetas_scattered.resize(nchannels);
        require_array(group_id, "phi", phis_scattered.size(), phis_scattered);
        trace("The original phi values are:\n", phis_scattered);
        require_array(group_id, "theta", thetas_scattered.size(), thetas_scattered);
        trace("The original theta values are:\n", thetas_scattered);
        if (database_version < 1 && repair_legacy_hover_directivity(phis_scattered))
        {
            repaired_spheres++;
        }
        if (i == 0)
        {
            // Sort and identify unique elements.  D does
            //   phis = phis_scattered.dup.sort().uniq().array();
            // i.e. on a COPY.  phis_scattered/thetas_scattered must keep all
            // `nchannels` per-channel values, because the grid fill below
            // pairs them index-wise with dBAs_scattered/EAAs_scattered.
            // Sorting them in place collapsed them to the unique values, so
            // the match loop ran over 36 entries instead of 1406, left most
            // grid cells at 0, and indexed dBAs_scattered with the wrong c.
            phis = phis_scattered;
            std::sort(phis.begin(), phis.end());
            phis.erase(std::unique(phis.begin(), phis.end()), phis.end());
            thetas = thetas_scattered;
            std::sort(thetas.begin(), thetas.end());
            thetas.erase(std::unique(thetas.begin(), thetas.end()), thetas.end());
            info("The ", phis.size(), " unique phi values are:\n", phis);
            info("The ", thetas.size(), " unique theta values are:\n", thetas);
            // Allocate noise level matrices
            dBA.assign(phis.size(),
                    std::vector<std::vector<double>>(thetas.size(),
                            std::vector<double>(group_ids.size(), 0.0)));
            EAA.assign(phis.size(),
                    std::vector<std::vector<double>>(thetas.size(),
                            std::vector<double>(group_ids.size(), 0.0)));
        }
        // Fetch scattered acoustic data for each sphere
        dBAs_scattered.resize(nchannels);
        EAAs_scattered.resize(nchannels);
        require_array(group_id, "dBA", dBAs_scattered.size(), dBAs_scattered);
        require_array(group_id, "EAA", EAAs_scattered.size(), EAAs_scattered);
        // Search through scattered data for grided values
        for (std::size_t m = 0; m < phis.size(); m++)
        {
            auto phi_target = phis[m];
            for (std::size_t n = 0; n < thetas.size(); n++)
            {
                auto theta_target = thetas[n];
                for (std::size_t c = 0; c < phis_scattered.size(); c++)
                {
                    auto phi_val = phis_scattered[c];
                    if (approx_equal(phi_target, phi_val)
                            && approx_equal(theta_target, thetas_scattered[c]))
                    {
                        // Load sphere data into grid
                        dBA[m][n][i] = dBAs_scattered[c];
                        EAA[m][n][i] = EAAs_scattered[c];
                        if (dBAs_scattered[c] > dba_max)
                        {
                            dba_max = dBAs_scattered[c];
                        }
                    }
                }
            }
        }
    }
    if (repaired_spheres > 0)
    {
        log("Repaired legacy hover directivity on ", repaired_spheres, " of ", group_ids.size(),
            " spheres.  Regenerate this database with a versioned builder to avoid the fix-up.");
    }
    max_level_ = dba_max;
    info("Max level found in database is: ", max_level_);
    // Set sphere_radius to dimensional_radius
    sphere_radius_ = dimensional_radius;
    info("Sphere radius set to ", sphere_radius_);
    // Find range of conditions
    advance_ratio_max = *std::max_element(advance_ratios.begin(), advance_ratios.end());
    advance_ratio_min = *std::min_element(advance_ratios.begin(), advance_ratios.end());
    flight_path_angle_max = *std::max_element(flight_path_angles.begin(), flight_path_angles.end());
    flight_path_angle_min = *std::min_element(flight_path_angles.begin(), flight_path_angles.end());
    thrust_coefficient_max = *std::max_element(thrust_coefficients.begin(), thrust_coefficients.end());
    thrust_coefficient_min = *std::min_element(thrust_coefficients.begin(), thrust_coefficients.end());
    info("Operating condition ranges are:\n", "μ  = [", advance_ratio_min, " , ",
            advance_ratio_max, "]\n", "γ  = [", flight_path_angle_min, " , ", flight_path_angle_max,
            "]\n", "CT = [", thrust_coefficient_min, " , ", thrust_coefficient_max, "]\n");
    // Normalize operating conditions
    for (auto& x : advance_ratios)
        x = (x - advance_ratio_min) / (advance_ratio_max - advance_ratio_min);
    trace("Normalized μ   = ", advance_ratios);
    for (auto& x : flight_path_angles)
        x = (x - flight_path_angle_min) / (flight_path_angle_max - flight_path_angle_min);
    trace("Normalized γ   = ", flight_path_angles);
    for (auto& x : thrust_coefficients)
        x = (x - thrust_coefficient_min) / (thrust_coefficient_max - thrust_coefficient_min);
    trace("Normalized CT  = ", thrust_coefficients);

    if (set_interpolator3(advance_ratios.size(), advance_ratios.data(), flight_path_angles.data(),
                thrust_coefficients.data()))
    {
        fatal("Could not set interpolator.");
    }
}

double Database::eval(double range, double phi, double theta, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    // Find the sphere grid points bounding the desired phi, theta
    auto phi_range = bracket(phis, phi);
    auto theta_range = bracket(thetas, theta);
    trace("Found phi (", phi, ") in range ", phis[phi_range[0]], " -- ", phis[phi_range[1]]);
    trace("Found theta (", theta, ") in range ", thetas[theta_range[0]], " -- ",
            thetas[theta_range[1]]);

    // Normalize flight condition and limit range to [0.0,1.0]
    advance_ratio = clamp01((advance_ratio - advance_ratio_min) / (advance_ratio_max - advance_ratio_min));
    flight_path_angle = clamp01(
            (flight_path_angle - flight_path_angle_min) / (flight_path_angle_max - flight_path_angle_min));
    thrust_coefficient = clamp01(
            (thrust_coefficient - thrust_coefficient_min) / (thrust_coefficient_max - thrust_coefficient_min));
    trace("Normalized operating condition is: ", "μ  = ", advance_ratio, ", γ  = ", flight_path_angle,
            ", CT = ", thrust_coefficient);

    // Get barycentric weights and indicies
    std::array<std::size_t, 4> indicies;
    std::array<double, 4> weights;

    barycentric(advance_ratio, flight_path_angle, thrust_coefficient, weights.data(), indicies.data());
    trace("Barycentric weights: ", weights, " and indicies ", indicies);

    // Barycentric interpolation between conditions at each gridpoint
    // Bilinear interpolation between gridpoints
    const auto p1 = phis[phi_range[0]];
    const auto p2 = phis[phi_range[1]];
    const auto t1 = thetas[theta_range[0]];
    const auto t2 = thetas[theta_range[1]];

    // Blend the four bracketing grid points of one metric.  Writing this once
    // rather than twice over eight transcribed corner indices removes the main
    // way this can go quietly wrong.
    const auto interpolate = [&](const std::vector<std::vector<std::vector<double>>>& metric) {
        return periodic_bilinear_interp(
                barycentric_interp(metric, phi_range[0], theta_range[0], weights, indicies),
                barycentric_interp(metric, phi_range[0], theta_range[1], weights, indicies),
                barycentric_interp(metric, phi_range[1], theta_range[0], weights, indicies),
                barycentric_interp(metric, phi_range[1], theta_range[1], weights, indicies),
                p1, p2, t1, t2, phi, theta);
    };
    const auto sphere_level = interpolate(dBA);
    const auto excess_factor = interpolate(EAA);
    trace("Sphere level (dB): ", sphere_level, "\tExcess Factor (dB/km): ", excess_factor);
    // Apply spherical spreading and excess atmospheric absorption to surface level
    // TODO Add ground effects
    const double ground_effect = 4.0;
    auto excess_attenuation = ground_effect * excess_factor * (range - sphere_radius()) / 1000.0;
    auto spreading_loss = 20.0 * std::log10(range / sphere_radius());
    trace("Spreading Loss (dB): ", spreading_loss, "\tExcess Attenuation (dB): ", excess_attenuation);
    assert(std::isfinite(sphere_level));
    assert(std::isfinite(excess_factor));
    assert(std::isfinite(excess_attenuation));
    assert(std::isfinite(spreading_loss));
    auto ground_level = sphere_level - spreading_loss - excess_attenuation;
    trace("Ground level: ", ground_level);
    return ground_level;
}

double Database::level_at_observer(const Eigen::Vector3d& observer, const Eigen::Vector3d& source,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    Eigen::Vector3d r = observer - source;
    auto slant_range = r.norm();
    trace("R = ", slant_range, " for r = ", r.transpose());
    auto angles = lookup_angles(source, velocity, normal, observer);
    trace("Looked up [phi,theta] = ", angles[0], ", ", angles[1]);
    assert(std::isfinite(angles[0]) && std::isfinite(angles[1]));
    auto level = eval(slant_range, angles[0], angles[1], advance_ratio, flight_path_angle,
            thrust_coefficient);
    assert(std::isfinite(level));
    return level;
}

double Database::level_along_bearing(double range, double bearing, double altitude,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    bearing *= PI / 180.0;
    // Define ground velocity, and get heading vector from it
    double ground_velocity[2] = {velocity(0), velocity(1)};
    double hnorm = std::sqrt(ground_velocity[0] * ground_velocity[0] + ground_velocity[1] * ground_velocity[1]);
    double heading[2] = {ground_velocity[0] / hnorm, ground_velocity[1] / hnorm};
    // Place observer on ground at heading rotated by bearing
    // NOTE: second component faithfully reproduces the D original, which uses
    // heading[0] twice (likely a latent bug in level_along_bearing).
    Eigen::Vector3d observer = {range * (std::cos(bearing) * heading[0] + std::sin(bearing) * heading[1]),
            range * (std::sin(-bearing) * heading[0] + std::cos(bearing) * heading[0]), 0.0};
    // Put source at altitude above ground
    Eigen::Vector3d source = {0.0, 0.0, altitude};
    return level_at_observer(observer, source, velocity, normal, advance_ratio, flight_path_angle,
            thrust_coefficient);
}

double Database::zero_range_level(double altitude, const Eigen::Vector3d& velocity,
        const Eigen::Vector3d& normal, double advance_ratio, double flight_path_angle,
        double thrust_coefficient) const
{
    Eigen::Vector3d observer = {0.0, 0.0, 0.0};
    Eigen::Vector3d source = {0.0, 0.0, altitude};
    return level_at_observer(observer, source, velocity, normal, advance_ratio, flight_path_angle,
            thrust_coefficient);
}

double Database::find_range(double contour_level, double bearing, double altitude,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient, double threshold,
        bool zero_range_check) const
{
    // Check if level directly below sphere (i.e. range = 0) is above the requested contour level.
    auto zero_level = zero_range_level(altitude, velocity, normal, advance_ratio, flight_path_angle,
            thrust_coefficient);
    trace("Zero range level is: ", zero_level);
    // Initialize range to zero
    double range = 0.0;
    // If level is greater than contour level, then search for the range associated with the contour level.
    if (!zero_range_check || zero_level > contour_level)
    {
        trace("Proceeding, zero range level above contour level.");
        // Range error
        auto residual = [this, bearing, altitude, &velocity, &normal, advance_ratio, flight_path_angle,
                        thrust_coefficient, contour_level](double x) {
            auto level = level_along_bearing(x, bearing, altitude, velocity, normal, advance_ratio,
                    flight_path_angle, thrust_coefficient);
            return level - contour_level;
        };
        // Set maximum range bound from max_level assuming spherical spreading only
        auto range_max = sphere_radius() * std::pow(10.0, 0.05 * (max_level_ - contour_level));
        trace("Max search range set to: ", range_max);
        range = bisect(residual, 0.0, range_max, threshold);
        trace("Found range: ", range);
        // Verify solution if zero_range_check is false
        if (!zero_range_check)
        {
            auto check_level = level_along_bearing(range, bearing, altitude, velocity, normal, advance_ratio,
                    flight_path_angle, thrust_coefficient);
            trace("Check level: ", check_level);
            if (check_level < (contour_level - threshold))
            {
                // Bad solution, set range back to zero
                range = 0.0;
                trace("Reset bad range to 0.0.");
            }
        }
    }
    return range;
}

void Database::self_test(const std::string& datafile_path)
{
    std::size_t invalid = 0;
    std::size_t valid = 0;
    log("Begin database self test against file: ", datafile_path);
    // Open netCDF file
    int file_ncid;
    info("Opening netCDF file: ", datafile_path);
    if (NC_NOERR != nc_open(datafile_path.c_str(), NC_NOWRITE, &file_ncid))
    {
        fatal("Can't open netCDF file at ", datafile_path);
    }
    NcGuard guard(file_ncid);

    const int self_test_version = get_int_scalar(file_ncid, "database_version", 0);
    auto group_ids = get_group_ids(file_ncid);
    // Local operating parameters
    double advance_ratio, flight_path_angle, thrust_coefficient, radius, scale, dimensional_radius;
    // Cycle through sphere groups and verify data
    for (std::size_t i = 0; i < group_ids.size(); i++)
    {
        auto group_id = group_ids[i];
        // Get the sphere parameters
        require_scalar(group_id, "advance_ratio", &advance_ratio);
        require_scalar(group_id, "flight_path_angle", &flight_path_angle);
        require_scalar(group_id, "thrust_coefficient", &thrust_coefficient);
        require_scalar(group_id, "radius", &radius);
        require_scalar(group_id, "rotor_scale", &scale);
        dimensional_radius = radius * scale;
        log("Checking sphere ", group_id, " with \t μ = ", advance_ratio, ",\t γ = ",
                flight_path_angle, ",\t CT = ", thrust_coefficient, ",\t scale = ", scale,
                ",\t radius = ", radius);

        // Get the data for this sphere
        auto nchannels = get_dim_length("channels", group_id);
        std::vector<double> theta_data(nchannels);
        require_array(group_id, "theta", nchannels, theta_data);

        std::vector<double> phi_data(nchannels);
        require_array(group_id, "phi", nchannels, phi_data);
        // Interpret the file exactly as the loader did, or this compares the
        // grid built from repaired labels against queries at the raw ones.
        if (self_test_version < 1)
        {
            repair_legacy_hover_directivity(phi_data);
        }

        std::vector<double> dBA_data(nchannels);
        require_array(group_id, "dBA", nchannels, dBA_data);

        std::vector<double> EAA_data(nchannels);
        require_array(group_id, "EAA", nchannels, EAA_data);

        // Validate the database against the known data at sphere surface
        for (std::size_t ch = 0; ch < dBA_data.size(); ch++)
        {
            auto dBA = dBA_data[ch];
            auto level = eval(dimensional_radius, phi_data[ch], theta_data[ch], advance_ratio,
                    flight_path_angle, thrust_coefficient);
            if (!approx_equal(dBA, level))
            {
                log("Surface level disagreement of D:", dBA, " versus N:", level, " for φ = ",
                        phi_data[ch], ", θ = ", theta_data[ch], ", μ  = ", advance_ratio, ",γ  =  ",
                        flight_path_angle, ", CT = ", thrust_coefficient);
                invalid++;
            }
            else
            {
                valid++;
            }
        }
    }
    if (!invalid)
    {
        log("Self test passed on all ", valid, " samples.");
    }
    else
    {
        fatal("Self test failed on ", invalid, " of ", valid + invalid, " samples.");
    }
}

} // namespace niceops
