#include "sphere_database.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <filesystem>
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

// Length of an optional dimension.  Returns the netCDF status rather than
// aborting, so a caller can treat "absent" as a normal answer.
int nc_inq_dimid_len(int group_id, const char* dimname, std::size_t* length)
{
    int dimid = 0;
    const int status = nc_inq_dimid(group_id, dimname, &dimid);
    if (status != NC_NOERR)
        return status;
    return nc_inq_dimlen(group_id, dimid, length);
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

// The half-open sub-range of the sorted grid axis `a` whose entries are
// approx_equal to `value`.  approx_equal accepts exactly
// |x - value| <= max(1e-2*|value|, 1e-5), so the accepted entries are
// contiguous; the binary search is given a slightly wider bound and every
// surviving candidate is then re-tested with approx_equal itself, so the
// accepted set is identical to what a full scan of the axis would give.
std::array<std::size_t, 2> approx_equal_range(const std::vector<double>& a, double value)
{
    const double tol = std::max(1e-2 * std::fabs(value), 1e-5) * 1.001 + 1e-9;
    auto lo = std::lower_bound(a.begin(), a.end(), value - tol);
    auto hi = std::upper_bound(a.begin(), a.end(), value + tol);
    while (lo != hi && !approx_equal(*lo, value))
        ++lo;
    while (hi != lo && !approx_equal(*(hi - 1), value))
        --hi;
    return {static_cast<std::size_t>(lo - a.begin()), static_cast<std::size_t>(hi - a.begin())};
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

namespace {

// The corner-independent half of a periodic bilinear interpolation: whether the
// cell wraps around the sphere on either axis, and the deltas the corner
// weights are built from.  eval() interpolates dBA and EAA over the very same
// cell, so this is derived once and applied to both -- the corner values are
// all that differ.
struct BilinearCell
{
    bool swap_x;
    bool swap_y;
    double x2x, y2y, xx1, yy1, area;
};

// Expects the axes already in order; periodic_bilinear_cell puts them there.
BilinearCell bilinear_cell(double x1, double x2, double y1, double y2, double x, double y)
{
    BilinearCell cell;
    cell.swap_x = false;
    cell.swap_y = false;

    auto x2x1 = x2 - x1;
    auto y2y1 = y2 - y1;
    cell.x2x = x2 - x;
    cell.y2y = y2 - y;
    cell.yy1 = y - y1;
    cell.xx1 = x - x1;
    assert(x2x1 >= 0.0 && y2y1 >= 0.0);

    if (x2x1 < 10 * std::numeric_limits<double>::epsilon())
    {
        x2x1 = 2.0;
        cell.x2x = 1.0;
        cell.xx1 = 1.0;
    }
    else if (y2y1 < 10 * std::numeric_limits<double>::epsilon())
    {
        y2y1 = 2.0;
        cell.y2y = 1.0;
        cell.yy1 = 1.0;
    }
    cell.area = x2x1 * y2y1;
    return cell;
}

BilinearCell periodic_bilinear_cell(double x1, double x2, double y1, double y2, double x, double y)
{
    // Check if points "wrap around" sphere, and swap inputs to bilinear
    // interpolation if necessary
    const bool swap_x = x1 > x2;
    if (swap_x)
        std::swap(x1, x2);
    const bool swap_y = y1 > y2;
    if (swap_y)
        std::swap(y1, y2);
    auto cell = bilinear_cell(x1, x2, y1, y2, x, y);
    cell.swap_x = swap_x;
    cell.swap_y = swap_y;
    return cell;
}

double bilinear_apply(const BilinearCell& cell, double q11, double q12, double q21, double q22)
{
    if (cell.swap_x)
    {
        std::swap(q11, q21);
        std::swap(q12, q22);
    }
    if (cell.swap_y)
    {
        std::swap(q11, q12);
        std::swap(q21, q22);
    }
    return (q11 * cell.x2x * cell.y2y + q21 * cell.xx1 * cell.y2y + q12 * cell.x2x * cell.yy1
                   + q22 * cell.xx1 * cell.yy1)
            / cell.area;
}

} // namespace

double bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y)
{
    return bilinear_apply(bilinear_cell(x1, x2, y1, y2, x, y), q11, q12, q21, q22);
}

double periodic_bilinear_interp(double q11, double q12, double q21, double q22, double x1, double x2,
        double y1, double y2, double x, double y)
{
    return bilinear_apply(periodic_bilinear_cell(x1, x2, y1, y2, x, y), q11, q12, q21, q22);
}

namespace {

// Turn an upper_bound-style index (the last entry not greater than the target,
// or -1) into the bracketing pair, wrapping around the ends.
std::array<std::size_t, 2> wrap_bracket(long index, std::size_t n)
{
    if (index < 0 || index >= static_cast<long>(n) - 1)
        return {n - 1, 0};
    return {static_cast<std::size_t>(index), static_cast<std::size_t>(index) + 1};
}

} // namespace

std::array<std::size_t, 2> bracket(const std::vector<double>& a, double val)
{
    // Grab the range up to target value (first element strictly greater than val)
    auto it = std::upper_bound(a.begin(), a.end(), val);
    // Index is one less the length of this range
    return wrap_bracket(static_cast<long>(it - a.begin()) - 1, a.size());
}

void GridAxis::reset(std::vector<double> values)
{
    values_ = std::move(values);
    first_ = values_.empty() ? 0.0 : values_.front();
    inv_step_ = 0.0;
    if (values_.size() < 2)
        return;
    const double step = (values_.back() - first_) / static_cast<double>(values_.size() - 1);
    if (!(step > 0.0))
        return;
    // Take the shortcut only if every entry sits within a tenth of a step of
    // where a perfectly uniform axis would put it, so the fix-up in bracket()
    // never has to walk more than one cell.
    for (std::size_t i = 0; i < values_.size(); i++)
    {
        if (std::fabs(values_[i] - (first_ + step * static_cast<double>(i))) > 0.1 * step)
            return;
    }
    inv_step_ = 1.0 / step;
}

std::array<std::size_t, 2> GridAxis::bracket(double val) const
{
    if (inv_step_ == 0.0)
        return niceops::bracket(values_, val);

    const long n = static_cast<long>(values_.size());
    // Jump straight at the cell the uniform spacing implies...
    const double estimate = (val - first_) * inv_step_;
    long index;
    if (!(estimate >= 0.0)) // also catches NaN, which upper_bound leaves at the low end
        index = -1;
    else if (estimate >= static_cast<double>(n))
        index = n - 1;
    else
        index = static_cast<long>(estimate);
    // ...then restore upper_bound's exact invariant.  Rounding in the estimate,
    // or an axis that is uniform only to within the tolerance above, can leave
    // it one cell out; otherwise neither loop runs.
    while (index >= 0 && values_[index] > val)
        --index;
    while (index + 1 < n && values_[index + 1] <= val)
        ++index;
    return wrap_bracket(index, values_.size());
}

double clamp01(double x)
{
    return std::clamp(x, 0.0, 1.0);
}

double acosc(double x)
{
    // Clamp first: round-off can push a dot product of unit vectors a hair
    // past +/-1, where acos would return NaN.
    return std::acos(std::clamp(x, -1.0, 1.0));
}

SphereFrame sphere_frame(const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal)
{
    SphereFrame frame;
    frame.normal = normal.normalized();
    // Project velocity vector onto TPP to get a "nose" vector
    const Eigen::Vector3d v = velocity.normalized();
    frame.nose = v - v.dot(frame.normal) * frame.normal;
    frame.nose.normalize();
    // Define starboard as the cross product of the "nose" and the TPP normal
    frame.starboard = frame.nose.cross(frame.normal);
    return frame;
}

std::array<double, 2> lookup_angles(const SphereFrame& frame, const Eigen::Vector3d& direction)
{
    // Sphere theta angle is relative to projection of velocity onto TPP
    double theta = acosc(frame.nose.dot(direction));
    // Project the propagation vector onto the plane normal to the "nose"
    Eigen::Vector3d rp = direction - direction.dot(frame.nose) * frame.nose;
    rp.normalize();
    // Then, we can get sphere phi by computing the angle of the propagation
    // vector in the plane
    double phi = acosc(std::fabs(frame.normal.dot(rp)));
    // Sign is still uncertain, so check which side of starboard we are on.
    // acos is monotonically decreasing, so acos(t) > pi/2 is exactly t < 0 --
    // the comparison does not need the angle, only its sign, which saves an
    // acos, a divide and a norm on a path that runs tens of millions of times.
    if (frame.starboard.dot(rp) < 0.0)
    {
        phi = -phi;
    }
    // Convert to degrees and return phi and theta
    phi *= 180.0 / PI;
    theta *= 180.0 / PI;
    return {phi, theta};
}

std::array<double, 2> lookup_angles(const Eigen::Vector3d& source, const Eigen::Vector3d& velocity,
        const Eigen::Vector3d& normal, const Eigen::Vector3d& observer)
{
    return lookup_angles(sphere_frame(velocity, normal), (observer - source).normalized());
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

    // Optional: builders since the versioned format record the rotor scale at
    // the root, which is enough to redimensionalize without a vehicle file.
    double tip_speed_value = 0.0;
    if (NC_NOERR == get_scalar(file_ncid, "main_rotor_tip_speed_meters_per_sec", &tip_speed_value))
    {
        tip_speed_ = tip_speed_value;
        info("Database records a main rotor tip speed of ", tip_speed_, " m/s.");
    }

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
    auto eaa_min = std::numeric_limits<double>::infinity();
    sphere_max_.assign(group_ids.size(), -std::numeric_limits<double>::infinity());
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
        // Shaft power, when the builder wrote it.  Optional: no shipped database
        // carries it, and a database without it simply has the feature off.
        // The D original (DRZ_edits) made it required, which would refuse every
        // database in this repository.
        double shaft_power_value = 0.0;
        if (NC_NOERR == get_scalar(group_id, "ShaftPower", &shaft_power_value))
        {
            if (i != 0 && shaft_power_.size() != i)
                fatal("Database error: ShaftPower present on some spheres but not others.");
            shaft_power_.push_back(shaft_power_value);
        }
        else if (!shaft_power_.empty())
        {
            fatal("Database error: ShaftPower present on some spheres but not others.");
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
            auto unique_sorted = [](std::vector<double> v) {
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                return v;
            };
            phis.reset(unique_sorted(phis_scattered));
            thetas.reset(unique_sorted(thetas_scattered));
            info("The ", phis.size(), " unique phi values are:\n", phis.values());
            info("The ", thetas.size(), " unique theta values are:\n", thetas.values());
            // Allocate the noise level grid
            n_conditions = group_ids.size();
            grid.assign(phis.size() * thetas.size() * n_conditions, Sample {0.0, 0.0});
        }
        // Fetch scattered acoustic data for each sphere
        dBAs_scattered.resize(nchannels);
        EAAs_scattered.resize(nchannels);
        require_array(group_id, "dBA", dBAs_scattered.size(), dBAs_scattered);
        require_array(group_id, "EAA", EAAs_scattered.size(), EAAs_scattered);
        // Optional per-band spectra, written gridded over the sorted unique
        // angles.  Absent from every database built before the builder learned
        // to carry them, so its absence is not an error -- but a database that
        // has them on some spheres and not others is malformed.
        std::size_t n_bands = 0;
        if (NC_NOERR == nc_inq_dimid_len(group_id, "frequency", &n_bands))
        {
            std::vector<double> bands(n_bands);
            require_array(group_id, "frequency", n_bands, bands);
            if (i == 0)
            {
                frequencies_ = bands;
                spectra_.assign(phis.size() * thetas.size() * group_ids.size() * n_bands, 0.0);
                info("Database carries ", n_bands, " frequency bands from ", bands.front(),
                        " to ", bands.back(), " Hz.");
            }
            else if (bands != frequencies_)
            {
                fatal("Database error: frequency bands differ between spheres.");
            }
            int amplitude_id = 0;
            if (NC_NOERR != nc_inq_varid(group_id, "amplitude", &amplitude_id))
                fatal("Database error: sphere has a frequency dimension but no amplitude.");
            if (NC_NOERR
                    != nc_get_var_double(group_id, amplitude_id,
                            spectra_.data() + i * phis.size() * thetas.size() * n_bands))
                fatal("Database error: could not read amplitude for sphere ", i, ".");
        }
        else if (!frequencies_.empty())
        {
            fatal("Database error: spectra present on some spheres but not others.");
        }

        // Place each channel into the grid cells it matches, rather than
        // scanning every channel for every cell.  The original nesting cost
        // n_phi * n_theta * n_channels approx_equal probes per sphere -- 2.5
        // billion over the shipped database, most of the load time -- while
        // each channel in practice lands in exactly one cell.  The acceptance
        // test is unchanged; the search only skips candidates approx_equal
        // could not accept.  Two semantics carry over: a later channel
        // overwrites an earlier one in a shared cell (channels are still
        // visited in file order), and dba_max covers matched channels only.
        for (std::size_t c = 0; c < phis_scattered.size(); c++)
        {
            const auto phi_match = approx_equal_range(phis.values(), phis_scattered[c]);
            if (phi_match[0] == phi_match[1])
                continue;
            const auto theta_match = approx_equal_range(thetas.values(), thetas_scattered[c]);
            if (theta_match[0] == theta_match[1])
                continue;
            if (dBAs_scattered[c] > dba_max)
            {
                dba_max = dBAs_scattered[c];
            }
            if (dBAs_scattered[c] > sphere_max_[i])
            {
                sphere_max_[i] = dBAs_scattered[c];
            }
            if (EAAs_scattered[c] < eaa_min)
            {
                eaa_min = EAAs_scattered[c];
            }
            for (std::size_t m = phi_match[0]; m < phi_match[1]; m++)
            {
                for (std::size_t n = theta_match[0]; n < theta_match[1]; n++)
                {
                    // Load sphere data into grid
                    grid[(m * thetas.size() + n) * n_conditions + i]
                            = {dBAs_scattered[c], EAAs_scattered[c]};
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
    min_excess_factor_ = eaa_min;
    info("Max level found in database is: ", max_level_);
    info("Lowest excess attenuation factor is: ", min_excess_factor_);
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

Database::Blend Database::blend_at(double advance_ratio, double flight_path_angle,
        double thrust_coefficient) const
{
    // Normalize flight condition and limit range to [0.0,1.0]
    advance_ratio = clamp01((advance_ratio - advance_ratio_min) / (advance_ratio_max - advance_ratio_min));
    flight_path_angle = clamp01(
            (flight_path_angle - flight_path_angle_min) / (flight_path_angle_max - flight_path_angle_min));
    thrust_coefficient = clamp01(
            (thrust_coefficient - thrust_coefficient_min) / (thrust_coefficient_max - thrust_coefficient_min));
    trace("Normalized operating condition is: ", "μ  = ", advance_ratio, ", γ  = ", flight_path_angle,
            ", CT = ", thrust_coefficient);

    // Get barycentric weights and indicies
    Blend blend;
    barycentric(advance_ratio, flight_path_angle, thrust_coefficient, blend.weights.data(),
            blend.indices.data());
    trace("Barycentric weights: ", blend.weights, " and indicies ", blend.indices);
    return blend;
}

Database::SourceState Database::source_state(const Eigen::Vector3d& source,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    return {source, sphere_frame(velocity, normal),
            blend_at(advance_ratio, flight_path_angle, thrust_coefficient)};
}

double Database::eval(double range, double phi, double theta, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    return eval(range, phi, theta, blend_at(advance_ratio, flight_path_angle, thrust_coefficient));
}

Database::Sample Database::blend_cell(std::size_t i_phi, std::size_t i_theta,
        const Blend& blend) const
{
    const Sample* cell = grid.data() + (i_phi * thetas.size() + i_theta) * n_conditions;
    const Sample& s0 = cell[blend.indices[0]];
    const Sample& s1 = cell[blend.indices[1]];
    const Sample& s2 = cell[blend.indices[2]];
    const Sample& s3 = cell[blend.indices[3]];
    const auto& w = blend.weights;
    return {w[0] * s0.dBA + w[1] * s1.dBA + w[2] * s2.dBA + w[3] * s3.dBA,
            w[0] * s0.EAA + w[1] * s1.EAA + w[2] * s2.EAA + w[3] * s3.EAA};
}

double Database::eval(double range, double phi, double theta, const Blend& blend) const
{
    // Find the sphere grid points bounding the desired phi, theta
    auto phi_range = phis.bracket(phi);
    auto theta_range = thetas.bracket(theta);
    trace("Found phi (", phi, ") in range ", phis[phi_range[0]], " -- ", phis[phi_range[1]]);
    trace("Found theta (", theta, ") in range ", thetas[theta_range[0]], " -- ",
            thetas[theta_range[1]]);

    // Barycentric interpolation between conditions at each gridpoint
    // Bilinear interpolation between gridpoints
    const auto p1 = phis[phi_range[0]];
    const auto p2 = phis[phi_range[1]];
    const auto t1 = thetas[theta_range[0]];
    const auto t2 = thetas[theta_range[1]];

    // Blend the operating conditions at each of the four bracketing grid
    // points, then bilinearly interpolate between them.  Both metrics travel
    // together: they are read from the same cells, so splitting them would
    // touch every cache line twice.
    const auto q11 = blend_cell(phi_range[0], theta_range[0], blend);
    const auto q12 = blend_cell(phi_range[0], theta_range[1], blend);
    const auto q21 = blend_cell(phi_range[1], theta_range[0], blend);
    const auto q22 = blend_cell(phi_range[1], theta_range[1], blend);
    const auto cell = periodic_bilinear_cell(p1, p2, t1, t2, phi, theta);
    const auto sphere_level = bilinear_apply(cell, q11.dBA, q12.dBA, q21.dBA, q22.dBA);
    const auto excess_factor = bilinear_apply(cell, q11.EAA, q12.EAA, q21.EAA, q22.EAA);
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

double Database::shaft_power(const Blend& blend) const
{
    if (shaft_power_.empty())
        return std::numeric_limits<double>::quiet_NaN();
    const auto& w = blend.weights;
    const auto& i = blend.indices;
    return w[0] * shaft_power_[i[0]] + w[1] * shaft_power_[i[1]] + w[2] * shaft_power_[i[2]]
            + w[3] * shaft_power_[i[3]];
}

Database::Sphere Database::interpolate_sphere(double advance_ratio, double flight_path_angle,
        double thrust_coefficient) const
{
    Sphere sphere;
    sphere.advance_ratio = advance_ratio;
    sphere.flight_path_angle = flight_path_angle;
    sphere.thrust_coefficient = thrust_coefficient;
    sphere.radius = sphere_radius_;
    sphere.shaft_power = shaft_power(blend_at(advance_ratio, flight_path_angle, thrust_coefficient));
    sphere.phis = phis.values();
    sphere.thetas = thetas.values();

    // Same blend eval() applies at one direction, carried out at every cell, so
    // a sphere exported at an operating condition the database already holds
    // comes back bit-for-bit as that sphere's own data.
    const auto blend = blend_at(advance_ratio, flight_path_angle, thrust_coefficient);
    sphere.cells.resize(phis.size() * thetas.size());
    for (std::size_t m = 0; m < phis.size(); m++)
    {
        for (std::size_t n = 0; n < thetas.size(); n++)
        {
            sphere.cells[m * thetas.size() + n] = blend_cell(m, n, blend);
        }
    }
    if (!frequencies_.empty())
    {
        // Blend the bands with the same weights as the level, in dB, which is
        // how the level itself is blended.
        sphere.frequencies = frequencies_;
        const std::size_t n_bands = frequencies_.size();
        const std::size_t cells = sphere.cells.size();
        sphere.amplitudes.assign(cells * n_bands, 0.0);
        const auto& w = blend.weights;
        for (std::size_t cell = 0; cell < sphere.cells.size(); cell++)
        {
            for (std::size_t band = 0; band < n_bands; band++)
            {
                double value = 0.0;
                for (std::size_t k = 0; k < 4; k++)
                {
                    value += w[k]
                            * spectra_[(blend.indices[k] * cells + cell) * n_bands + band];
                }
                sphere.amplitudes[cell * n_bands + band] = value;
            }
        }
    }
    info("Interpolated sphere at \u03bc = ", advance_ratio, ", \u03b3 = ", flight_path_angle,
            ", CT = ", thrust_coefficient, " over ", sphere.cells.size(), " cells.");
    return sphere;
}

double Database::cull_range(const Blend& blend, double cutoff) const
{
    const double unbounded = std::numeric_limits<double>::infinity();
    if (std::isnan(cutoff))
        return unbounded;
    // A level is sphere_level - 20*log10(r/r0) - 4*EAA*(r - r0)/1000.  Dropping
    // the attenuation term bounds the level from above only if that term cannot
    // be negative, i.e. if no excess factor in the database is.  The D original
    // used to take abs() of the factor precisely because it could be; if this
    // database is one of those, decline to bound anything.
    if (!(min_excess_factor_ >= 0.0))
        return unbounded;

    // The blended sphere cannot exceed the loudest of the spheres it blends:
    // the barycentric weights are a convex combination, and so are the bilinear
    // weights applied on top of them.
    double blended_max = -unbounded;
    for (std::size_t k = 0; k < blend.weights.size(); k++)
    {
        // Skip the padding in a degenerate blend, where all the weight sits on
        // one vertex and the other three indices are left at zero.
        if (blend.weights[k] > 0.0)
            blended_max = std::max(blended_max, sphere_max_[blend.indices[k]]);
    }
    if (!std::isfinite(blended_max))
        return unbounded;

    // Spherical spreading alone takes that peak below the cutoff at this range,
    // and excess attenuation only takes it further down.
    return sphere_radius_ * std::pow(10.0, 0.05 * (blended_max - cutoff));
}

double Database::level_at_observer(const Eigen::Vector3d& observer, const SourceState& state) const
{
    Eigen::Vector3d r = observer - state.position;
    auto slant_range = r.norm();
    trace("R = ", slant_range, " for r = ", r.transpose());
    // Reuse the slant range as the normalization: Eigen's normalized() divides
    // by this same sqrt(squaredNorm()) (returning the vector untouched when it
    // is zero), so this is the identical direction without a second subtraction
    // and square root.
    if (slant_range > 0.0)
    {
        r /= slant_range;
    }
    auto angles = lookup_angles(state.frame, r);
    trace("Looked up [phi,theta] = ", angles[0], ", ", angles[1]);
    assert(std::isfinite(angles[0]) && std::isfinite(angles[1]));
    auto level = eval(slant_range, angles[0], angles[1], state.blend);
    assert(std::isfinite(level));
    return level;
}

double Database::level_at_observer(const Eigen::Vector3d& observer, const Eigen::Vector3d& source,
        const Eigen::Vector3d& velocity, const Eigen::Vector3d& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    return level_at_observer(observer,
            source_state(source, velocity, normal, advance_ratio, flight_path_angle,
                    thrust_coefficient));
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

namespace {

// Small wrappers so the writer below reads as a list of what it writes rather
// than a wall of status checks.
void nc_check(int status, const std::string& what)
{
    if (status != NC_NOERR)
        fatal("netCDF error while ", what, ": ", nc_strerror(status));
}

int define_scalar(int group, const char* name, int dimid)
{
    int varid = 0;
    nc_check(nc_def_var(group, name, NC_DOUBLE, 1, &dimid, &varid), std::string("defining ") + name);
    return varid;
}

} // namespace

void write_sphere(const Database::Sphere& sphere, const std::string& path)
{
    const auto directory = std::filesystem::path(path).parent_path();
    if (!directory.empty())
        std::filesystem::create_directories(directory);

    int file = 0;
    nc_check(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file), "creating " + path);
    NcGuard guard(file);

    // Root: the flags Database's loader looks for.  Version 1 says the sphere
    // was written by a versioned builder, so the legacy hover repair is not
    // applied to it -- correct here, since the grid it came from was already
    // repaired on load if it needed to be.
    int scalar_dim = 0;
    nc_check(nc_def_dim(file, "scalar", 1, &scalar_dim), "defining scalar dimension");
    int same_grid_var = 0;
    nc_check(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid_var),
            "defining same_grid");
    int version_var = 0;
    nc_check(nc_def_var(file, "database_version", NC_INT, 1, &scalar_dim, &version_var),
            "defining database_version");

    int group = 0;
    nc_check(nc_def_grp(file, "sphere0", &group), "defining group sphere0");

    const std::size_t n_channels = sphere.cells.size();
    int radii_dim = 0, condition_dim = 0, channel_dim = 0;
    nc_check(nc_def_dim(group, "radii", 1, &radii_dim), "defining radii dimension");
    nc_check(nc_def_dim(group, "condition", 1, &condition_dim), "defining condition dimension");
    nc_check(nc_def_dim(group, "channels", n_channels, &channel_dim), "defining channels dimension");

    const int radius_var = define_scalar(group, "radius", radii_dim);
    const int scale_var = define_scalar(group, "rotor_scale", condition_dim);
    const int mu_var = define_scalar(group, "advance_ratio", condition_dim);
    const int gamma_var = define_scalar(group, "flight_path_angle", condition_dim);
    const int ct_var = define_scalar(group, "thrust_coefficient", condition_dim);
    const int phi_var = define_scalar(group, "phi", channel_dim);
    const int theta_var = define_scalar(group, "theta", channel_dim);
    const int dba_var = define_scalar(group, "dBA", channel_dim);
    const int eaa_var = define_scalar(group, "EAA", channel_dim);
    const bool has_spectra = !sphere.frequencies.empty();
    int band_dim = 0, band_var = 0, amplitude_var = 0;
    if (has_spectra)
    {
        nc_check(nc_def_dim(group, "frequency", sphere.frequencies.size(), &band_dim),
                "defining frequency dimension");
        band_var = define_scalar(group, "frequency", band_dim);
        int phi_dim = 0, theta_dim = 0;
        nc_check(nc_def_dim(group, "phi_grid", sphere.phis.size(), &phi_dim),
                "defining phi_grid dimension");
        nc_check(nc_def_dim(group, "theta_grid", sphere.thetas.size(), &theta_dim),
                "defining theta_grid dimension");
        const int dims[3] = {phi_dim, theta_dim, band_dim};
        nc_check(nc_def_var(group, "amplitude", NC_DOUBLE, 3, dims, &amplitude_var),
                "defining amplitude");
    }
    const bool has_power = !std::isnan(sphere.shaft_power);
    const int power_var = has_power ? define_scalar(group, "ShaftPower", condition_dim) : -1;

    nc_check(nc_enddef(file), "ending define mode");

    const short same_grid = 1;
    const int version = 1;
    nc_check(nc_put_var_short(file, same_grid_var, &same_grid), "writing same_grid");
    nc_check(nc_put_var_int(file, version_var, &version), "writing database_version");

    // The loader multiplies radius by rotor_scale, so put the whole dimensional
    // radius in the first and leave the second at unity.
    const double scale = 1.0;
    nc_check(nc_put_var_double(group, radius_var, &sphere.radius), "writing radius");
    nc_check(nc_put_var_double(group, scale_var, &scale), "writing rotor_scale");
    nc_check(nc_put_var_double(group, mu_var, &sphere.advance_ratio), "writing advance_ratio");
    nc_check(nc_put_var_double(group, gamma_var, &sphere.flight_path_angle),
            "writing flight_path_angle");
    nc_check(nc_put_var_double(group, ct_var, &sphere.thrust_coefficient),
            "writing thrust_coefficient");

    // Flatten the grid back into the format's scattered channel lists.
    std::vector<double> phi_channel(n_channels), theta_channel(n_channels);
    std::vector<double> dba(n_channels), eaa(n_channels);
    for (std::size_t m = 0; m < sphere.phis.size(); m++)
    {
        for (std::size_t n = 0; n < sphere.thetas.size(); n++)
        {
            const std::size_t c = m * sphere.thetas.size() + n;
            phi_channel[c] = sphere.phis[m];
            theta_channel[c] = sphere.thetas[n];
            dba[c] = sphere.cells[c].dBA;
            eaa[c] = sphere.cells[c].EAA;
        }
    }
    nc_check(nc_put_var_double(group, phi_var, phi_channel.data()), "writing phi");
    nc_check(nc_put_var_double(group, theta_var, theta_channel.data()), "writing theta");
    nc_check(nc_put_var_double(group, dba_var, dba.data()), "writing dBA");
    nc_check(nc_put_var_double(group, eaa_var, eaa.data()), "writing EAA");
    if (has_power)
    {
        nc_check(nc_put_var_double(group, power_var, &sphere.shaft_power), "writing ShaftPower");
    }
    if (has_spectra)
    {
        nc_check(nc_put_var_double(group, band_var, sphere.frequencies.data()),
                "writing frequency");
        nc_check(nc_put_var_double(group, amplitude_var, sphere.amplitudes.data()),
                "writing amplitude");
    }

    log("Wrote interpolated sphere to ", path, " (", n_channels, " channels).");
}

void write_aam_sphere(const Database::Sphere& sphere, const std::string& path, double tip_speed)
{
    if (sphere.frequencies.empty())
    {
        fatal("Cannot write ", path,
                " in AAM format: the database carries no spectra, and AMPLITUDE(PHI, THETA, "
                "FREQUENCY) is the whole payload of the format.");
    }
    const auto directory = std::filesystem::path(path).parent_path();
    if (!directory.empty())
        std::filesystem::create_directories(directory);

    const double m2ft = 1.0 / 0.3048;
    const double mps2kt = 1.94384;

    int file = 0;
    nc_check(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file), "creating " + path);
    NcGuard guard(file);

    // The AAM/RNM hemisphere layout, as the source spheres themselves carry it:
    // a handful of scalar descriptors, the three axes, and the spectrum.
    auto scalar_dim = [&](const char* name) {
        int dim = 0;
        nc_check(nc_def_dim(file, name, 1, &dim), std::string("defining ") + name);
        return dim;
    };
    auto scalar_var = [&](const char* name, int dim, const char* unit) {
        int var = 0;
        nc_check(nc_def_var(file, name, NC_FLOAT, 1, &dim, &var), std::string("defining ") + name);
        nc_check(nc_put_att_text(file, var, "unit", std::strlen(unit), unit),
                std::string("annotating ") + name);
        return var;
    };

    const int bb_var = scalar_var("BB", scalar_dim("BB"), "");
    const int nb_var = scalar_var("NB", scalar_dim("NB"), "");
    const int pt_var = scalar_var("PT", scalar_dim("PT"), "");
    const int dsr_var = scalar_var(
            "DOPPLER_SHIFT_REMOVED", scalar_dim("DOPPLER_SHIFT_REMOVED"), "");
    const int empty_var = scalar_var("EMPTY_WEIGHT", scalar_dim("EMPTY_WEIGHT"), "POUNDS");
    const int fuel_var = scalar_var("FUEL_WEIGHT", scalar_dim("FUEL_WEIGHT"), "POUNDS");
    const int load_var = scalar_var("LOAD_WEIGHT", scalar_dim("LOAD_WEIGHT"), "POUNDS");
    const int radius_var = scalar_var("RADIUS", scalar_dim("RADIUS"), "FEET");
    const int fpa_var = scalar_var("FLIGHT_PATH_ANGLE", scalar_dim("FLIGHT_PATH_ANGLE"), "DEGREE");
    const int pylon_var = scalar_var("PYLON_ANGLE", scalar_dim("PYLON_ANGLE"), "DEGREE");
    const int speed_var = scalar_var("SPEED", scalar_dim("SPEED"), "KNOTS");
    const int masttilt_var = scalar_var("MASTTILT", scalar_dim("MASTTILT"), "DEGREE");

    int xyz_dim = 0, phi_dim = 0, theta_dim = 0, freq_dim = 0;
    nc_check(nc_def_dim(file, "XYZ", 3, &xyz_dim), "defining XYZ dimension");
    nc_check(nc_def_dim(file, "PHI", sphere.phis.size(), &phi_dim), "defining PHI dimension");
    nc_check(nc_def_dim(file, "THETA", sphere.thetas.size(), &theta_dim),
            "defining THETA dimension");
    nc_check(nc_def_dim(file, "FREQUENCY", sphere.frequencies.size(), &freq_dim),
            "defining FREQUENCY dimension");

    const int xyz_var = scalar_var("XYZ", xyz_dim, "FEET");
    const int phi_var = scalar_var("PHI", phi_dim, "DEGREE");
    const int theta_var = scalar_var("THETA", theta_dim, "DEGREE");
    const int freq_var = scalar_var("FREQUENCY", freq_dim, "HERTZ");

    int amplitude_var = 0;
    const int amplitude_dims[3] = {phi_dim, theta_dim, freq_dim};
    nc_check(nc_def_var(file, "AMPLITUDE", NC_FLOAT, 3, amplitude_dims, &amplitude_var),
            "defining AMPLITUDE");
    nc_check(nc_put_att_text(file, amplitude_var, "unit", 7, "DECIBEL"), "annotating AMPLITUDE");

    const std::string title = "NICE-OPS interpolated sphere";
    nc_check(nc_put_att_text(file, NC_GLOBAL, "title", title.size(), title.c_str()),
            "writing title");
    nc_check(nc_enddef(file), "ending define mode");

    auto put_scalar = [&](int var, double value) {
        const float as_float = static_cast<float>(value);
        nc_check(nc_put_var_float(file, var, &as_float), "writing a scalar");
    };
    // Narrowband only: these spheres carry no broadband/tonal split.
    put_scalar(bb_var, 0.0);
    put_scalar(nb_var, 1.0);
    put_scalar(pt_var, 0.0);
    put_scalar(dsr_var, 0.0);
    // The database records no weight breakdown, and the format has no way to
    // say "unknown", so these stay zero rather than inventing a figure.
    put_scalar(empty_var, 0.0);
    put_scalar(fuel_var, 0.0);
    put_scalar(load_var, 0.0);
    put_scalar(radius_var, sphere.radius * m2ft);
    put_scalar(fpa_var, sphere.flight_path_angle);
    put_scalar(pylon_var, 90.0);
    put_scalar(speed_var, sphere.advance_ratio * tip_speed * mps2kt);
    put_scalar(masttilt_var, 0.0);

    const float origin[3] = {0.0F, 0.0F, 0.0F};
    nc_check(nc_put_var_float(file, xyz_var, origin), "writing XYZ");

    auto put_axis = [&](int var, const std::vector<double>& values) {
        std::vector<float> as_float(values.begin(), values.end());
        nc_check(nc_put_var_float(file, var, as_float.data()), "writing an axis");
    };
    put_axis(phi_var, sphere.phis);
    put_axis(theta_var, sphere.thetas);
    put_axis(freq_var, sphere.frequencies);

    std::vector<float> amplitude(sphere.amplitudes.begin(), sphere.amplitudes.end());
    nc_check(nc_put_var_float(file, amplitude_var, amplitude.data()), "writing AMPLITUDE");

    log("Wrote AAM sphere to ", path, " (", sphere.phis.size(), " x ", sphere.thetas.size(),
            " x ", sphere.frequencies.size(), ").");
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
    // Optional: builders since the versioned format record the rotor scale at
    // the root, which is enough to redimensionalize without a vehicle file.
    double tip_speed_value = 0.0;
    if (NC_NOERR == get_scalar(file_ncid, "main_rotor_tip_speed_meters_per_sec", &tip_speed_value))
    {
        tip_speed_ = tip_speed_value;
        info("Database records a main rotor tip speed of ", tip_speed_, " m/s.");
    }

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
        // Shaft power, when the builder wrote it.  Optional: no shipped database
        // carries it, and a database without it simply has the feature off.
        // The D original (DRZ_edits) made it required, which would refuse every
        // database in this repository.
        double shaft_power_value = 0.0;
        if (NC_NOERR == get_scalar(group_id, "ShaftPower", &shaft_power_value))
        {
            if (i != 0 && shaft_power_.size() != i)
                fatal("Database error: ShaftPower present on some spheres but not others.");
            shaft_power_.push_back(shaft_power_value);
        }
        else if (!shaft_power_.empty())
        {
            fatal("Database error: ShaftPower present on some spheres but not others.");
        }

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

        // Every channel of a sphere is queried at that sphere's own operating
        // condition, so the blend is the same for all of them -- resolve it
        // once instead of running a CGAL locate per channel.
        const auto blend = blend_at(advance_ratio, flight_path_angle, thrust_coefficient);

        // Validate the database against the known data at sphere surface
        for (std::size_t ch = 0; ch < dBA_data.size(); ch++)
        {
            auto dBA = dBA_data[ch];
            auto level = eval(dimensional_radius, phi_data[ch], theta_data[ch], blend);
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
