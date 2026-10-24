#include "sphere_database.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <system_error>
#include <thread>
#include <type_traits>
#include <utility>

#include <netcdf.h>

#include "manifest.hpp"
#include "propagation.hpp"
#include "util/logger.hpp"
#include "weighting.hpp"
#include "util/netcdf.hpp"
#include "util/paths.hpp"
#include "util/constants.hpp"
#include "util/rootfind.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// RAII replacement for D's scope(exit) nc_close.  Note: intentionally does not
// run on fatal() paths (fatal aborts the process).
struct NcGuard
{
    int id;
    explicit NcGuard(int i) : id(i) {}
    ~NcGuard()
    {
        if (id >= 0)
            nc_close(id);
    }
    // Closes now, and says how it went: with HDF5's buffering a failed write
    // (a full disk, a quota) surfaces only here, so a writer has to look.
    int close()
    {
        const int status = id >= 0 ? nc_close(id) : NC_NOERR;
        id = -1;
        return status;
    }
};

std::vector<int> get_group_ids(int ncid)
{
    // Get number of groups (spheres) in database
    int numgroups;
    if (NC_NOERR != nc_inq_grps(ncid, &numgroups, nullptr))
    {
        input_error("Can't inquire netCDF groups.");
    }
    // Get the group ids
    std::vector<int> group_ids(numgroups);
    if (NC_NOERR != nc_inq_grps(ncid, &numgroups, group_ids.data()))
    {
        input_error("Cannot get number of netCDF groups.");
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
        input_error("Unable to find ", dimname, "  dimension.");
    }
    std::size_t dim_length;
    if (NC_NOERR != nc_inq_dimlen(group_id, dim_id, &dim_length))
    {
        input_error("Unable to find length of ", dimname, "  dimension.");
    }
    return dim_length;
}

// `optional`: the variable may well be absent, so say nothing -- a database
// with 1312 groups and no ShaftPower otherwise wrote 1312 identical lines, 97%
// of a run's log, and buried the few that mattered.
template <typename T>
int get_scalar(int ncid, const std::string& variable_name, T* value, bool optional = false)
{
    int varid;
    auto status = nc_inq_varid(ncid, variable_name.c_str(), &varid);
    if (NC_NOERR != status)
    {
        if (!optional)
            log("Could not find variable ", variable_name, " in ncid ", ncid, ".");
    }
    else
    {
        // As double, whatever the file stores (nc_get_var copies the external
        // type's bytes unconverted: a float32 file loaded as garbage, and the
        // exit status was 0), and element 0 only (it read the whole variable, so
        // a "scalar" with a longer dimension overran the one double it was
        // handed).
        static_assert(std::is_same_v<T, double>, "scalars are read as double");
        const std::size_t first[NC_MAX_VAR_DIMS] = {};
        status = nc_get_var1_double(ncid, varid, first, value);
        if (NC_NOERR != status)
        {
            input_error("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
        }
    }
    return status;
}

// A root attribute holding one number, if the file carries it.  panam writes its
// build atmosphere and rotor speed this way; older databases carry neither.
std::optional<double> read_double_attribute(int ncid, const char* name)
{
    double value = 0.0;
    switch (netcdf::single_value(ncid, NC_GLOBAL, name, value))
    {
    case netcdf::Single::missing:
        return std::nullopt;
    case netcdf::Single::not_one_number:
        input_error("Found but could not read the database attribute ", name, " as one number.");
    case netcdf::Single::ok:
        break;
    }
    return value;
}

// A group's coverage, if it carries one: 1 where the builder measured the bin,
// 0 where the gate zeroed it or it was filled from elsewhere.  It is flattened
// in its own dimension order -- along the channels when the group has no
// spectra (as dBA is), or (phi, theta, frequency) when it does (as amplitude
// is) -- and the caller decides which it is by the length.
bool read_coverage(int group_id, std::vector<int>& coverage)
{
    int varid = 0;
    if (NC_NOERR != nc_inq_varid(group_id, "coverage", &varid))
        return false;
    const auto total = netcdf::element_count(group_id, varid);
    if (!total)
        input_error("Found but could not read the dimensions of the variable coverage.");
    coverage.resize(*total);
    if (NC_NOERR != nc_get_var_int(group_id, varid, coverage.data()))
        input_error("Found but could not read the variable coverage.");
    return true;
}

// A root text attribute, or nullopt when the file has none.  panam writes them
// with netCDF4-python, whose str attributes are NC_CHAR; NC_STRING is accepted
// too.  Anything else is refused by name rather than read as absent.
std::optional<std::string> read_root_text(int ncid, const char* name)
{
    nc_type type = NC_NAT;
    std::size_t length = 0;
    if (NC_NOERR != nc_inq_att(ncid, NC_GLOBAL, name, &type, &length))
        return std::nullopt;
    std::string text;
    if (type == NC_CHAR)
    {
        text.resize(length);
        if (NC_NOERR != nc_get_att_text(ncid, NC_GLOBAL, name, text.data()))
            input_error("Found but could not read the database attribute ", name, ".");
    }
    else if (type == NC_STRING && length == 1)
    {
        char* value = nullptr;
        if (NC_NOERR != nc_get_att_string(ncid, NC_GLOBAL, name, &value))
            input_error("Found but could not read the database attribute ", name, ".");
        text = value == nullptr ? "" : value;
        nc_free_string(1, &value);
    }
    else
    {
        input_error("The database attribute ", name, " is not text.");
    }
    // Fixed-width writers pad with NULs or blanks.
    text.erase(text.find_last_not_of(std::string(" \0", 2)) + 1);
    return text;
}

// The root attributes of a CoherenceRecord, each checked as it is read: a
// fit's scale comes with its exponent, numbers are positive, the tonal
// fraction within [0, 1] (its length is checked against the bands later).
CoherenceRecord read_coherence_record(int ncid)
{
    CoherenceRecord r;
    auto positive = [&](const char* name) -> std::optional<double> {
        const auto v = read_double_attribute(ncid, name);
        if (v && !(*v > 0.0 && std::isfinite(*v)))
            input_error("The database attribute ", name, " is ", *v, "; it must be positive.");
        return v;
    };
    r.blades = positive("main_rotor_blades");
    if (const auto text = read_root_text(ncid, "main_rotor_rotation"))
    {
        r.rotation_sense = parse_rotation_sense(*text);
        if (!r.rotation_sense)
            input_error("The database attribute main_rotor_rotation is '", *text,
                    "'; expected 'counterclockwise' or 'clockwise', seen from above.");
    }
    r.rotation_limit = positive("coherence_rotation_limit_Hz");
    r.scale = positive("coherence_scale_Hz_deg");
    r.exponent = positive("coherence_exponent");
    r.spread_scale = positive("coherence_spread_scale_Hz_deg");
    r.spread_exponent = positive("coherence_spread_exponent");
    if (r.scale.has_value() != r.exponent.has_value())
        input_error("The database gives one of coherence_scale_Hz_deg and coherence_exponent without the other.");
    if (r.spread_scale.has_value() != r.spread_exponent.has_value())
        input_error("The database gives one of coherence_spread_scale_Hz_deg and coherence_spread_exponent "
                    "without the other.");
    switch (netcdf::attribute_values(ncid, NC_GLOBAL, "coherence_tonal_fraction", r.tonal_fraction))
    {
    case netcdf::Values::missing:
        r.tonal_fraction.clear();
        break;
    case netcdf::Values::not_numbers:
        input_error("Found but could not read the database attribute coherence_tonal_fraction as numbers.");
    case netcdf::Values::ok:
        for (double v : r.tonal_fraction)
            if (!(v >= 0.0 && v <= 1.0))
                input_error("The database attribute coherence_tonal_fraction holds ", v,
                        "; a fraction is within [0, 1].");
        break;
    }
    if (const auto text = read_root_text(ncid, "coherence_source"))
        r.source = *text;
    return r;
}

// The root attribute speed_reference.  Absent is the ground.
SpeedReference read_speed_reference(int ncid)
{
    const auto text = read_root_text(ncid, "speed_reference");
    if (!text)
        return SpeedReference::ground;
    const auto reference = parse_speed_reference(*text);
    if (!reference)
    {
        input_error("The database attribute speed_reference is '", *text,
                "'; expected 'ground' or 'air'.");
    }
    return *reference;
}

// The root attribute azimuth_reference, nullopt when absent (which the
// Database reads as the track).
std::optional<AzimuthReference> read_azimuth_reference(int ncid)
{
    const auto text = read_root_text(ncid, "azimuth_reference");
    if (!text)
        return std::nullopt;
    const auto reference = parse_azimuth_reference(*text);
    if (!reference)
    {
        input_error("The database attribute azimuth_reference is '", *text,
                "'; expected 'track' or 'heading'.");
    }
    return reference;
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
    // Element 0 only, as get_scalar reads: nc_get_var_int would write every
    // value the variable holds.
    int value = fallback;
    const std::size_t first[NC_MAX_VAR_DIMS] = {};
    if (NC_NOERR != nc_get_var1_int(ncid, varid, first, &value))
    {
        input_error("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
    }
    return value;
}

template <typename T>
void require_scalar(int ncid, const std::string& variable_name, T* value)
{
    if (NC_NOERR != get_scalar(ncid, variable_name, value))
    {
        input_error(variable_name, " required but not found.");
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
        // A 1-D variable: nc_get_vara reads one start and count per dimension.
        int ndims = 0;
        if (NC_NOERR != nc_inq_varndims(ncid, varid, &ndims) || ndims != 1)
            input_error("Variable ", variable_name, " in ncid ", ncid, " must be 1-D.");
        std::size_t startp = 0;
        values.resize(countp);
        static_assert(std::is_same_v<T, double>, "arrays are read as double");
        status = nc_get_vara_double(ncid, varid, &startp, &countp, values.data());
        if (NC_NOERR != status)
        {
            input_error("Found but could not read variable ", variable_name, " in ncid ", ncid, ".");
        }
    }
    return status;
}

template <typename T>
void require_array(int ncid, const std::string& variable_name, std::size_t countp, std::vector<T>& values)
{
    if (NC_NOERR != get_array(ncid, variable_name, countp, values))
    {
        input_error(variable_name, " required but not found.");
    }
}

// Open a database, or say why it could not be opened.  "Can't open" alone
// leaves no way to tell a missing file from a corrupt one or a permissions
// problem.  Shared by the loader and the self test, which must agree on how the
// file is read.
int open_database(const std::string& datafile_path)
{
    int file_ncid = 0;
    info("Opening netCDF file: ", datafile_path);
    const int open_status = nc_open(datafile_path.c_str(), NC_NOWRITE, &file_ncid);
    if (NC_NOERR != open_status)
    {
        input_error("Can't open netCDF file at ", datafile_path, ": ", nc_strerror(open_status));
    }
    return file_ncid;
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

std::optional<SpeedReference> parse_speed_reference(const std::string& name)
{
    if (name == "ground")
        return SpeedReference::ground;
    if (name == "air")
        return SpeedReference::air;
    return std::nullopt;
}

const char* to_string(SpeedReference reference)
{
    return reference == SpeedReference::air ? "air" : "ground";
}

std::optional<AzimuthReference> parse_azimuth_reference(const std::string& name)
{
    if (name == "track")
        return AzimuthReference::track;
    if (name == "heading")
        return AzimuthReference::heading;
    return std::nullopt;
}

const char* to_string(AzimuthReference reference)
{
    return reference == AzimuthReference::heading ? "heading" : "track";
}

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

// One (dBA, EAA) pair in a single two-lane vector.  Sample is exactly this in
// memory, and every step of an evaluation -- blending the four operating
// conditions, then the four bilinear corners -- applies the same weight to both
// metrics, so carrying them together halves the instruction count of the whole
// interpolation.  Neither the arithmetic nor its order changes.
//
// Both compilers this builds under spell a two-double vector this way (clang
// defines __GNUC__ too, so the one test covers both); the fallback is exactly
// what the code said before and keeps it portable.  Either spelling indexes
// with [], so everything below this block is written once.
#if defined(__GNUC__)

typedef double Pair __attribute__((vector_size(2 * sizeof(double))));

#else

struct Pair
{
    double v[2];

    double operator[](std::size_t lane) const
    {
        return v[lane];
    }
};

inline Pair operator*(Pair p, double scale)
{
    return Pair {{p.v[0] * scale, p.v[1] * scale}};
}

inline Pair& operator+=(Pair& a, Pair b)
{
    a.v[0] += b.v[0];
    a.v[1] += b.v[1];
    return a;
}

#endif

inline Pair load_pair(const Database::Sample& s)
{
    return Pair {s.dBA, s.EAA};
}

// A level at or below this carries no energy: the load-time floor
// (no_energy_floor_db) is what -inf -- a direction ambient gating emptied, or a
// band it masked -- is stored as.
constexpr double empty_level_db = Database::no_energy_floor_db + 1.0;

inline bool is_empty(double level, double weight)
{
    return weight != 0.0 && level <= empty_level_db;
}

// A group panam never finished writing holds the netCDF fill value (9.97e36)
// in every entry.  It is finite, so the non-finite floor passes it, and it
// read as a 1e37 dB level: the September 2026 AS350, EC130, RO-44 and RO-66
// files gave an infinite SEL wherever a state reached their last condition.
// No level, attenuation factor or band amplitude is within 1e30 of it.
void refuse_fill_values(const std::string& path, std::size_t sphere, const char* variable,
        const std::vector<double>& values)
{
    const auto unwritten = std::count_if(values.begin(), values.end(),
            [](double v) { return std::isfinite(v) && std::fabs(v) > 1e30; });
    if (unwritten > 0)
        input_error("Database ", path, ": sphere ", sphere, " holds the netCDF fill value in ",
                unwritten, " of ", values.size(), " ", variable, " entries, so that group was never "
                "written. Rebuild the database (panam c53588f or later writes every group).");
}

// The median of `values`, reordering them; the mean of the middle two for an
// even count.
double median_of(std::vector<double>& values)
{
    assert(!values.empty());
    const auto middle = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
    std::nth_element(values.begin(), middle, values.end());
    if (values.size() % 2 == 1)
        return *middle;
    return 0.5 * (*middle + *std::max_element(values.begin(), middle));
}

// Blending in dB is linear in the level, so one empty entry at the floor pulls
// the whole blend toward it in proportion to its weight: a 20% weight on
// -1000 dB costs 200 dB.  Where a sphere simply has no measurement in a
// direction, that turned a neighboring measured value into silence -- the
// narrow-strip footprints and 10-25 dB holes seen against the 2017 approach
// data, worst where QSAM or sparse coverage pulls in a condition with large
// gated regions.  "Not measured" is not "silent", so empty entries drop out
// and the remaining weights are renormalized.  Only when every entry with a
// weight is empty is the result empty.
//
// When nothing is empty -- the common case -- callers take the original
// arithmetic in its original order, so every result that did not involve an
// empty entry is unchanged to the bit.
inline Pair blend_present(const Pair* p, const double* w)
{
    Pair sum {0.0, 0.0};
    double total = 0.0;
    for (std::size_t k = 0; k < 4; k++)
    {
        if (w[k] != 0.0 && p[k][0] > empty_level_db)
        {
            sum += p[k] * w[k];
            total += w[k];
        }
    }
    if (!(total > 0.0))
        return Pair {Database::no_energy_floor_db, 0.0};
    return sum * (1.0 / total);
}

// Scalar form for the spectral paths, which blend one band at a time.
inline double blend_present(const double* v, const double* w)
{
    double sum = 0.0, total = 0.0;
    for (std::size_t k = 0; k < 4; k++)
    {
        if (w[k] != 0.0 && v[k] > empty_level_db)
        {
            sum += v[k] * w[k];
            total += w[k];
        }
    }
    return total > 0.0 ? sum / total : Database::no_energy_floor_db;
}

// The four bracketing operating conditions of one grid cell, blended by their
// barycentric weights.  `cell` points at that cell's run of conditions, which
// the flat grid layout keeps contiguous.
inline Pair blend_conditions(const Database::Sample* cell, const Database::Blend& blend)
{
    const Pair p[4] = {load_pair(cell[blend.indices[0]]), load_pair(cell[blend.indices[1]]),
            load_pair(cell[blend.indices[2]]), load_pair(cell[blend.indices[3]])};
    const double* w = blend.weights.data();
    if (is_empty(p[0][0], w[0]) || is_empty(p[1][0], w[1]) || is_empty(p[2][0], w[2])
            || is_empty(p[3][0], w[3]))
        return blend_present(p, w);
    Pair blended = p[0] * w[0];
    blended += p[1] * w[1];
    blended += p[2] * w[2];
    blended += p[3] * w[3];
    return blended;
}

// The four corners of a sphere grid cell and the bilinear weight on each, in
// the fixed order (lo,lo), (lo,hi), (hi,lo), (hi,hi).  `corner` is the flat
// (phi, theta) cell index; eval() scales it by n_conditions to reach the grid
// and eval_spectrum() by n_bands to reach the spectra.
//
// Both paths derive their stencil here rather than each spelling it out.  The
// spectral path has to land on the identical geometry -- that is what makes a
// disagreement between the two the physics and never the interpolation -- and
// sharing GridAxis::locate() only got them the same cell, not the same weights.
struct CellStencil
{
    std::size_t corner[4];
    double weight[4];
};

inline CellStencil cell_stencil(const GridAxis::Cell& phi_cell, const GridAxis::Cell& theta_cell,
        std::size_t n_thetas)
{
    // The four corner weights are the products of the two axes' own, and they
    // sum to one, so the divide by the cell area that periodic_bilinear_interp()
    // ends with is already done.
    const double fp = phi_cell.fraction;
    const double ft = theta_cell.fraction;
    return {{phi_cell.lo * n_thetas + theta_cell.lo, phi_cell.lo * n_thetas + theta_cell.hi,
                    phi_cell.hi * n_thetas + theta_cell.lo,
                    phi_cell.hi * n_thetas + theta_cell.hi},
            {(1.0 - fp) * (1.0 - ft), (1.0 - fp) * ft, fp * (1.0 - ft), fp * ft}};
}

} // namespace

// The reference bilinear interpolation.  eval() no longer goes through it --
// it takes its weights from GridAxis::locate(), which is arithmetic rather than
// four knot lookups and a divide -- so this is now the *definition* the locate
// path is checked against, swept cell by cell in tests/test_sphere_database.cpp.
//
// It used to be split into a cell-deriving half and a corner-applying half so
// eval() could derive the cell once and apply it to dBA and EAA separately.
// Nothing needs that split now, so it is one function again.
double periodic_bilinear_interp(double q11, double q12, double q21, double q22, double x1, double x2,
        double y1, double y2, double x, double y, double x_period)
{
    // A cell across the seam of a periodic axis: its lower knot is x1 and its
    // upper knot is x2 one period on, so x is read into [x1, x2 + period).
    // Swapping the knots instead would interpolate the whole circle between
    // the two seam knots and extrapolate beyond them -- which is what this
    // function did before the phi seam was exercised at all.
    if (x_period > 0.0 && x1 > x2)
    {
        x2 += x_period;
        if (x < x1)
            x += x_period;
    }
    // Check if points "wrap around" sphere, and swap inputs to bilinear
    // interpolation if necessary
    else if (x1 > x2)
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

    auto x2x1 = x2 - x1;
    auto y2y1 = y2 - y1;
    auto x2x = x2 - x;
    auto y2y = y2 - y;
    auto yy1 = y - y1;
    auto xx1 = x - x1;
    assert(x2x1 >= 0.0 && y2y1 >= 0.0);

    // A degenerate axis (both knots the same point) has no direction to
    // interpolate along, so weight its two corners equally.
    if (x2x1 < 10 * std::numeric_limits<double>::epsilon())
    {
        x2x1 = 2.0;
        x2x = 1.0;
        xx1 = 1.0;
    }
    else if (y2y1 < 10 * std::numeric_limits<double>::epsilon())
    {
        y2y1 = 2.0;
        y2y = 1.0;
        yy1 = 1.0;
    }

    return (q11 * x2x * y2y + q21 * xx1 * y2y + q12 * x2x * yy1 + q22 * xx1 * yy1)
            / (x2x1 * y2y1);
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

void GridAxis::set_period(double period)
{
    period_ = period;
}

void GridAxis::reset(std::vector<double> values)
{
    values_ = std::move(values);
    first_ = values_.empty() ? 0.0 : values_.front();
    inv_step_ = 0.0;
    exactly_uniform_ = false;
    if (values_.size() < 2)
        return;
    const double step = (values_.back() - first_) / static_cast<double>(values_.size() - 1);
    if (!(step > 0.0))
        return;
    // Take the shortcut only if every entry sits within a tenth of a step of
    // where a perfectly uniform axis would put it, so the fix-up in bracket()
    // never has to walk more than one cell.  Track the worst deviation as well:
    // locate() reads the interpolation weight straight off the spacing and so
    // needs the knots to be where the spacing says, not merely near.
    double worst = 0.0;
    for (std::size_t i = 0; i < values_.size(); i++)
    {
        const double slip = std::fabs(values_[i] - (first_ + step * static_cast<double>(i)));
        if (slip > 0.1 * step)
            return;
        if (slip > worst)
            worst = slip;
    }
    inv_step_ = 1.0 / step;
    // A few ULP of the axis span.  Exact multiples of 10 and 5 degrees, which
    // is what every shipped database's axes are, reconstruct to zero slip.
    exactly_uniform_ = worst <= 8.0 * std::numeric_limits<double>::epsilon()
                    * std::fabs(values_.back() - first_);
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

GridAxis::Cell GridAxis::locate(double val) const
{
    const std::size_t n = values_.size();
    if (exactly_uniform_)
    {
        // One multiply gives both the cell and the position within it.  The
        // correction loops bracket() runs are not needed here: they only ever
        // move the index when round-off puts the estimate a hair on the wrong
        // side of a knot, and at a knot the two cells agree -- the fraction is
        // 1 at the top of one and 0 at the bottom of the next, so either
        // answer interpolates to the same value.
        const double estimate = (val - first_) * inv_step_;
        if (estimate >= 0.0 && estimate < static_cast<double>(n - 1))
        {
            const double index = std::floor(estimate);
            return {static_cast<std::size_t>(index), static_cast<std::size_t>(index) + 1,
                    estimate - index};
        }
        // Past either end: the wrap-around cell, below.
    }

    if (period_ > 0.0 && n >= 2)
    {
        // Bring the value into [first, first + period), then look at the seam:
        // the cell from the last knot up to the first knot one period on.
        const double first = values_.front();
        double u = val;
        if (!(u >= first && u < first + period_))
        {
            u = first + std::fmod(u - first, period_);
            if (u < first)
                u += period_;
        }
        if (u >= values_.back())
        {
            const double width = first + period_ - values_.back();
            return {n - 1, 0, (u - values_.back()) / width};
        }
        val = u;
    }

    const auto pair = bracket(val);
    const double lo = values_[pair[0]];
    const double hi = values_[pair[1]];
    // Reproduce periodic_bilinear_interp(): it sorts the two knots, so when the
    // pair wraps (lo > hi, the cell that spans the seam) the weights come out
    // measured from the *sorted* lower knot, which is `hi`.
    const double span = hi > lo ? hi - lo : lo - hi;
    if (span < 10 * std::numeric_limits<double>::epsilon())
        return {pair[0], pair[1], 0.5};
    return {pair[0], pair[1], hi > lo ? (val - lo) / span : (lo - val) / span};
}

// The narrowest theta bin, in degrees.  The bins are uniform in the square of
// the cosine, which makes them widest at the pole *and* at the equator and
// narrowest near 45 degrees -- so this is a scan, not "the last bin".  Taking
// the last bin would return 20.7 degrees against the true 7.2 and under-pad the
// stencil, which is the one thing the exact-skip scheme cannot survive.
double DirectionBounds::narrowest_theta_bin_degrees()
{
    constexpr int half = theta_bins / 2;
    double narrowest = 180.0;
    double previous = 90.0;
    for (int k = 1; k <= half; k++)
    {
        const double edge = std::acos(std::sqrt(static_cast<double>(k) / half)) * 180.0 / pi;
        narrowest = std::min(narrowest, previous - edge);
        previous = edge;
    }
    return narrowest;
}

namespace {

// (x - low) / (high - low), or 0 on an axis that does not vary, where it was 0/0.
// The expression is the one that was written inline, so every axis that does
// vary scales to the same bits as before.
double unit_interval(double x, double low, double high)
{
    const double span = high - low;
    return span > 0.0 ? (x - low) / span : 0.0;
}

// Widest gap between adjacent knots of an axis, in degrees.
double grid_spacing(const GridAxis& axis)
{
    double widest = 0.0;
    for (std::size_t i = 1; i < axis.size(); i++)
        widest = std::max(widest, std::fabs(axis[i] - axis[i - 1]));
    return widest;
}

} // namespace

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

SphereFrame sphere_frame(const Vec3& velocity, const Vec3& normal)
{
    SphereFrame frame;
    frame.normal = normal.normalized();
    // Project velocity vector onto TPP to get a "nose" vector
    const Vec3 v = velocity.normalized();
    frame.nose = v - v.dot(frame.normal) * frame.normal;
    if (!(frame.nose.norm() > 1e-9))
    {
        // No velocity, or one along the thrust axis (hover without a heading, a
        // pure climb): there is no direction to project, and a zero nose made both
        // frame vectors zero, so every observer was looked up at (phi 0, theta
        // 90) -- the same point of the sphere from every side.  Any direction in
        // the plane will do; take +x's projection (or +y's if the thrust axis is
        // along x), so the frame is a proper one and the same each time.
        const Vec3 axis = std::fabs(frame.normal.x) < 0.9 ? Vec3(1.0, 0.0, 0.0) : Vec3(0.0, 1.0, 0.0);
        frame.nose = axis - axis.dot(frame.normal) * frame.normal;
    }
    frame.nose.normalize();
    // Define starboard as the cross product of the "nose" and the TPP normal
    frame.starboard = frame.nose.cross(frame.normal);
    return frame;
}

namespace {

// Sphere phi, in degrees, from the propagation direction's components along the
// TPP normal and starboard.  The frame is orthonormal, so these two components
// *are* the direction's projection into the plane normal to the nose, up to a
// common positive scale the angle does not depend on -- which is why neither
// caller normalizes anything, and why this takes the raw components rather than
// a vector.
//
// What the projection used to compute is acos(|up| / hypot(up, across)),
// negated when `across` is, and that is exactly atan2(across, |up|) -- without
// the normalization, and with one inverse trig call rather than a magnitude and
// a sign test.  Written as atan of the ratio rather than atan2, which costs
// twice as much for quadrant and sign handling that a non-negative second
// argument does not need: the ratio overflows to an infinity of the right sign
// when up is zero, and atan takes that to the right half of pi.
//
// Shared by lookup_angles and Database::level_from_geometry, which differ only
// in how they get their components; the 0/0 guard is delicate enough that one
// copy of it is the point.
// Above the plane (up > 0) the angle comes round the far side of the circle,
// copysign(pi, across) - atan(across / up): the same value as atan2(across,
// -up), with the one atan the rest of the function already pays for.  Below the
// plane, or when the upper half is folded onto the lower, it stays
// atan(across / |up|), so every lower-half direction takes the same arithmetic
// as before the upper half was read.
inline double phi_degrees(double up, double across, bool fold_upper)
{
    const double phi = up > 0.0 && !fold_upper
            ? std::copysign(pi, across) - std::atan(across / up)
            : std::atan(across / std::fabs(up));
    if (up == 0.0 && across == 0.0)
    {
        // Straight along the nose: the projection degenerates and phi is
        // undefined -- and 0/0 is a NaN, which has to be caught here rather
        // than propagated into the grid lookup.  normalized() leaves the zero
        // vector alone, so the projection form landed on acos(0).
        return 0.5 * pi * (180.0 / pi);
    }
    return phi * (180.0 / pi);
}

} // namespace

std::array<double, 2> lookup_angles(const SphereFrame& frame, const Vec3& direction,
        bool fold_upper)
{
    // Sphere theta angle is relative to projection of velocity onto TPP
    const double along = frame.nose.dot(direction);
    const double theta = acosc(along);
    return {phi_degrees(frame.normal.dot(direction), frame.starboard.dot(direction), fold_upper),
            theta * (180.0 / pi)};
}

std::array<double, 2> lookup_angles(const Vec3& source, const Vec3& velocity,
        const Vec3& normal, const Vec3& observer)
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

Database::Database(const std::string& datafile_path, const DatabaseLoadOptions& load)
{
    const int file_ncid = open_database(datafile_path);
    NcGuard guard(file_ncid);
    root_attributes_ = netcdf_root_attributes(file_ncid);

    // Database format version.  Absent (or 0) means the file predates the
    // versioned builder, and its hover spheres may carry the broken directivity
    // described in repair_legacy_hover_directivity().
    const int database_version = get_int_scalar(file_ncid, "database_version", 0);
    info("Database format version: ", database_version);

    // Check for "same_grid" flag to determine if interpolation is required.
    // Missing is fatal rather than silently undefined -- it used to be read
    // into an uninitialized int, so an older or hand-made file could take
    // either branch at random.
    const int same_grid_flag = get_int_scalar(file_ncid, "same_grid", -1);
    if (same_grid_flag < 0)
    {
        input_error("Database error: required variable same_grid is missing.");
    }
    const bool same_grid = same_grid_flag != 0;
    info("same_grid is set to: ", same_grid);

    // TODO: Enable interpolation when same_grid == false
    if (!same_grid)
        input_error("same_grid == false not supported");

    // Optional, absent or false in every database predating panam's
    // load_factors=None: whether every stored sphere is the same LF=1
    // reference condition, scaled to other load factors analytically rather
    // than stored.  See the member comment on fixed_load_factor_.
    fixed_load_factor_ = get_int_scalar(file_ncid, "fixed_load_factor", 0) != 0;
    info("fixed_load_factor is set to: ", fixed_load_factor_);

    // Optional, absent or false in every database predating panam's grid/band
    // deduplication: whether phi, theta and (if present) frequency are the
    // same on every sphere, in which case the builder wrote them once at the
    // root instead of once per group.  Normally true -- panam completes every
    // sphere onto a common grid and reduces every condition to the same
    // frequency bands -- but not guaranteed (mixed microphone layouts, or
    // per-file frequency resolution), so it is a flag panam checked at write
    // time, not an assumption made here.
    const bool shared_grid_and_frequency =
            get_int_scalar(file_ncid, "shared_grid_and_frequency", 0) != 0;
    info("shared_grid_and_frequency is set to: ", shared_grid_and_frequency);
    std::vector<double> shared_phis_scattered, shared_thetas_scattered;
    std::vector<double> shared_bands;
    if (shared_grid_and_frequency)
    {
        const auto shared_nchannels = get_dim_length("channels", file_ncid);
        shared_phis_scattered.resize(shared_nchannels);
        shared_thetas_scattered.resize(shared_nchannels);
        require_array(file_ncid, "phi", shared_nchannels, shared_phis_scattered);
        require_array(file_ncid, "theta", shared_nchannels, shared_thetas_scattered);
        info("Shared grid contains ", shared_nchannels, " channels.");
        std::size_t shared_n_bands = 0;
        if (NC_NOERR == nc_inq_dimid_len(file_ncid, "frequency", &shared_n_bands))
        {
            shared_bands.resize(shared_n_bands);
            require_array(file_ncid, "frequency", shared_n_bands, shared_bands);
            info("Shared frequency axis carries ", shared_n_bands, " bands from ",
                    shared_bands.front(), " to ", shared_bands.back(), " Hz.");
        }
    }

    // Optional: builders since the versioned format record the rotor scale at
    // the root, which is enough to redimensionalize without a vehicle file.
    double tip_speed_value = 0.0;
    if (NC_NOERR == get_scalar(file_ncid, "main_rotor_tip_speed_meters_per_sec", &tip_speed_value))
    {
        tip_speed_ = tip_speed_value;
        info("Database records a main rotor tip speed of ", tip_speed_, " m/s.");
    }
    double root_value = 0.0;
    if (NC_NOERR == get_scalar(file_ncid, "main_rotor_radius_meters", &root_value))
        rotor_radius_ = root_value;
    if (NC_NOERR == get_scalar(file_ncid, "vehicle_weight_newtons", &root_value))
        vehicle_weight_newtons_ = root_value;
    speed_reference_ = read_speed_reference(file_ncid);
    azimuth_reference_ = read_azimuth_reference(file_ncid);
    // NICEOPS's sphere is level: theta from the nose projected into the tip-path
    // plane, which is horizontal in steady flight, as in PANAM and NORAH2.  AAM's
    // spheres are pitched by the flight path angle, and panam rotates them to level
    // as it builds; a database that says it holds anything else is refused, not read
    // a flight path angle off.
    if (const auto frame = read_root_text(file_ncid, "sphere_frame"); frame && *frame != "level")
    {
        input_error("Database ", datafile_path, " records sphere_frame '", *frame,
                "'; NICEOPS reads level spheres only. Rebuild it with panam, which rotates "
                "flight-path spheres to level.");
    }
    // The atmosphere the broadband EAA was built in, when the builder recorded
    // all three of its values; the rotor speed the spheres were built at.
    const auto build_temperature = read_double_attribute(file_ncid, "build_temperature_K");
    const auto build_pressure = read_double_attribute(file_ncid, "build_pressure_kPa");
    const auto build_humidity = read_double_attribute(file_ncid, "build_relative_humidity_percent");
    if (build_temperature && build_pressure && build_humidity)
        build_atmosphere_ = BuildAtmosphere{*build_temperature, *build_pressure, *build_humidity};
    rotor_rpm_ = read_double_attribute(file_ncid, "rotor_rpm");
    coherence_record_ = read_coherence_record(file_ncid);
    if (!coherence_record_.empty())
        info("Database records its vehicle's coherence between emission directions",
                coherence_record_.source.empty() ? std::string(".") : " (" + coherence_record_.source + ").");
    if (build_atmosphere_)
        info("Built at ", build_atmosphere_->temperature_kelvin, " K, ",
                build_atmosphere_->pressure_kpa, " kPa, ",
                build_atmosphere_->relative_humidity_percent, " % RH.");
    info("Database flight conditions are referred to the ", to_string(speed_reference_), ".");
    info("Database sphere azimuth is measured from the ", to_string(azimuth_reference()),
            azimuth_reference_ ? "." : " (the file does not say; the track is the default).");

    auto group_ids = get_group_ids(file_ncid);
    // Every field below is filled from the loop over these groups.  With none,
    // the loop never runs and the tail of this function read an uninitialized
    // sphere radius and dereferenced max_element's end() on empty condition
    // lists -- a segfault, with nothing on stderr because the log had already
    // recorded a plausible-looking load.  Say what is wrong instead.
    if (group_ids.empty())
    {
        input_error("Database ", datafile_path, " holds no operating-condition groups; "
                "a .nod file needs one group per sphere.");
    }

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
    std::size_t groups_with_upper = 0;
    std::size_t groups_with_coverage = 0;
    bool any_coverage_variable = false;
    std::vector<int> coverage_scattered;
    // A-weighting per band, for the grid-shaped coverage's energy shares.
    std::vector<double> coverage_gains;
    // Which grid cells each channel falls in.  It depends only on the channels'
    // angle labels, and every sphere shares one layout but the legacy hover
    // spheres, so it is derived once and again only when the labels change:
    // 1312 spheres x 1406 channels x 2 searches was most of what the fill cost.
    struct ChannelCells
    {
        std::array<std::size_t, 2> phi, theta;
    };
    std::vector<ChannelCells> placement;
    std::vector<double> placed_phis, placed_thetas;
    bool placed = false;
    std::size_t layout_warnings = 0;
    // Which grid entries a channel filled; the rest hold the {0, 0} the grid
    // starts at (the layout warning reports them), which the suspect-cell
    // screen must not take for a measured 0 dBA.
    std::vector<char> reached_cells;
    // The spectral amplitude layout is a function of the same labels.
    std::vector<std::size_t> cached_cell_of;
    bool cached_cell_of_valid = false;
    // Cycle through sphere groups and load data
    // The Doppler convention of the spheres: 0 received, the frame the
    // builders write and the goldens use; 1 de-Dopplerized, stationary-frame
    // spheres that need a frequency shift before absorption and weighting.
    // Absent (-1) is what every database written so far carries.
    int doppler_removed = -1;
    for (std::size_t i = 0; i < group_ids.size(); i++)
    {
        auto group_id = group_ids[i];
        require_scalar(group_id, "advance_ratio", &advance_ratio);
        require_scalar(group_id, "flight_path_angle", &flight_path_angle);
        require_scalar(group_id, "thrust_coefficient", &thrust_coefficient);
        require_scalar(group_id, "radius", &radius);
        require_scalar(group_id, "rotor_scale", &scale);
        const int group_doppler = get_int_scalar(group_id, "DOPPLER_SHIFT_REMOVED", -1);
        if (i == 0)
            doppler_removed = group_doppler;
        else if (group_doppler != doppler_removed)
            input_error("Database error: sphere ", group_id, " has DOPPLER_SHIFT_REMOVED ",
                    group_doppler, " but the first has ", doppler_removed, ". Spheres in one "
                    "database must share a Doppler convention.");
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
            input_error("Database error: sphere radius*scale not constant.");
        }
        // Shaft power, when the builder wrote it.  Optional: no shipped database
        // carries it, and a database without it simply has the feature off.
        // The D original (DRZ_edits) made it required, which would refuse every
        // database in this repository.
        double shaft_power_value = 0.0;
        if (NC_NOERR == get_scalar(group_id, "ShaftPower", &shaft_power_value, true))
        {
            if (i != 0 && shaft_power_.size() != i)
                input_error("Database error: ShaftPower present on some spheres but not others.");
            shaft_power_.push_back(shaft_power_value);
        }
        else if (!shaft_power_.empty())
        {
            input_error("Database error: ShaftPower present on some spheres but not others.");
        }

        std::size_t nchannels;
        if (shared_grid_and_frequency)
        {
            // Already read once at the root; every sphere shares it verbatim,
            // so there is nothing here to re-read or repair.
            nchannels = shared_phis_scattered.size();
            phis_scattered = shared_phis_scattered;
            thetas_scattered = shared_thetas_scattered;
        }
        else
        {
            nchannels = get_dim_length("channels", group_id);
            // Read THIS sphere's own channel angles.  They are not identical
            // across spheres: in a pre-versioned database the hover spheres
            // label their mirrored half differently from the forward-flight
            // ones, so reusing sphere 0's arrays files their levels into the
            // wrong grid cells.
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
        }
        if (std::any_of(phis_scattered.begin(), phis_scattered.end(),
                    [](double phi) { return std::fabs(phi) > 90.0; }))
            groups_with_upper++;
        // Channel-shaped coverage is placed with the channels below.  Grid-shaped
        // coverage (with spectra) is read in the spectral block, which knows the
        // cell order it has to follow, and counts there.
        const bool group_has_coverage = read_coverage(group_id, coverage_scattered);
        const bool channel_coverage = group_has_coverage && coverage_scattered.size() == nchannels;
        any_coverage_variable = any_coverage_variable || group_has_coverage;
        if (channel_coverage)
            groups_with_coverage++;
        if (i == 0)
        {
            // Record number of channels to check for changes
            nchannels0 = nchannels;
            info("First sphere contains ", nchannels, " channels.");
        }
        else if (nchannels != nchannels0)
        {
            input_error("Database error: number of channels is not constant.");
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
            phis.set_period(360.0);
            thetas.reset(unique_sorted(thetas_scattered));
            info("The ", phis.size(), " unique phi values are:\n", phis.values());
            info("The ", thetas.size(), " unique theta values are:\n", thetas.values());
            // Allocate the noise level grid
            n_conditions = group_ids.size();
            grid.assign(phis.size() * thetas.size() * n_conditions, Sample {0.0, 0.0});
            reached_cells.assign(grid.size(), 0);
            measured_.assign(grid.size(), 0.0F);
        }
        if (!placed || placed_phis != phis_scattered || placed_thetas != thetas_scattered)
        {
            placement.resize(phis_scattered.size());
            for (std::size_t c = 0; c < phis_scattered.size(); c++)
            {
                placement[c].phi = approx_equal_range(phis.values(), phis_scattered[c]);
                placement[c].theta = approx_equal_range(thetas.values(), thetas_scattered[c]);
            }
            placed_phis = phis_scattered;
            placed_thetas = thetas_scattered;
            placed = true;
            cached_cell_of_valid = false;
            // Whether this layout fills the grid.  A cell no channel reaches stays
            // 0 dBA, which blends read as a real level; a channel that matches no
            // grid angle (the axes come from sphere 0 alone) is dropped.  Neither
            // used to be mentioned anywhere.
            std::size_t dropped = 0;
            std::vector<char> reached(phis.size() * thetas.size(), 0);
            for (const auto& channel : placement)
            {
                if (channel.phi[0] == channel.phi[1] || channel.theta[0] == channel.theta[1])
                {
                    dropped++;
                    continue;
                }
                for (std::size_t m = channel.phi[0]; m < channel.phi[1]; m++)
                    for (std::size_t n = channel.theta[0]; n < channel.theta[1]; n++)
                        reached[m * thetas.size() + n] = 1;
            }
            const std::size_t unreached = static_cast<std::size_t>(std::count(reached.begin(), reached.end(), 0));
            if (dropped > 0 || unreached > 0)
            {
                if (++layout_warnings <= 3)
                {
                    warn("Database ", datafile_path, ": sphere ", i, " leaves ", unreached, " of ", reached.size(),
                            " (phi, theta) cells without a channel, which read as 0 dBA, and drops ", dropped,
                            " of ", placement.size(), " channels that match no angle of the grid (taken from "
                            "sphere 0).");
                }
            }
        }
        // Fetch scattered acoustic data for each sphere
        dBAs_scattered.resize(nchannels);
        EAAs_scattered.resize(nchannels);
        require_array(group_id, "dBA", dBAs_scattered.size(), dBAs_scattered);
        require_array(group_id, "EAA", EAAs_scattered.size(), EAAs_scattered);
        refuse_fill_values(datafile_path, i, "dBA", dBAs_scattered);
        refuse_fill_values(datafile_path, i, "EAA", EAAs_scattered);
        // Optional per-band spectra, written gridded over the sorted unique
        // angles.  Absent from every database built before the builder learned
        // to carry them, so its absence is not an error -- but a database that
        // has them on some spheres and not others is malformed.
        std::size_t n_bands = 0;
        bool has_spectrum;
        if (shared_grid_and_frequency && !shared_bands.empty())
        {
            // Already read once at the root; every sphere shares it verbatim.
            has_spectrum = true;
            n_bands = shared_bands.size();
            if (i == 0)
            {
                frequencies_ = shared_bands;
                spectra_.assign(phis.size() * thetas.size() * group_ids.size() * n_bands, 0.0);
            }
        }
        else
        {
            has_spectrum = NC_NOERR == nc_inq_dimid_len(group_id, "frequency", &n_bands);
            if (has_spectrum)
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
                    input_error("Database error: frequency bands differ between spheres.");
                }
            }
        }
        if (has_spectrum)
        {
            int amplitude_id = 0;
            if (NC_NOERR != nc_inq_varid(group_id, "amplitude", &amplitude_id))
                input_error("Database error: sphere has a frequency dimension but no amplitude.");
            // netCDF stores amplitude(phi, theta, frequency), so it arrives
            // cell-major with the bands innermost.  Read it as it lies, then
            // scatter into the condition-innermost layout the observer loop
            // wants -- see the note on spectra_ in the header.  This costs one
            // pass over the sphere at load time and saves a cache miss per
            // condition per band on every one of ~94M evaluations.
            std::vector<double> sphere_amplitudes(phis.size() * thetas.size() * n_bands);
            // nc_get_var reads every value the variable holds, so it must hold
            // exactly this many.
            if (netcdf::element_count(group_id, amplitude_id) != sphere_amplitudes.size())
                input_error("Database error: amplitude for sphere ", i, " does not hold ", phis.size(), " x ",
                        thetas.size(), " x ", n_bands, " values (phi, theta, frequency).");
            if (NC_NOERR
                    != nc_get_var_double(group_id, amplitude_id, sphere_amplitudes.data()))
                input_error("Database error: could not read amplitude for sphere ", i, ".");
            refuse_fill_values(datafile_path, i, "amplitude", sphere_amplitudes);
            const std::size_t n_conditions_total = group_ids.size();
            // Which sorted (phi, theta) cell each row and column of the array
            // is.  The format means amplitude to be gridded over the sorted
            // unique angles, but panam before 2026-09-25 wrote its rows in the
            // order it completed the sphere (-90..90, then the mirrored upper
            // surface), which read as sorted put every spectrum in the wrong
            // direction.  When the channels are a full row-major grid -- channel
            // r * n_theta + j is row r, column j, as both writers lay them out --
            // their own angles say which cell each row is, and for a sorted file
            // that is the identity.  Otherwise keep the sorted reading.
            const std::size_t n_phi = phis.size(), n_theta = thetas.size();
            std::vector<std::size_t>& cell_of = cached_cell_of;
            if (!cached_cell_of_valid)
            {
                cell_of.resize(n_phi * n_theta);
                std::iota(cell_of.begin(), cell_of.end(), std::size_t {0});
                if (phis_scattered.size() == n_phi * n_theta)
                {
                    std::vector<std::size_t> mapped(n_phi * n_theta);
                    std::vector<bool> taken(n_phi * n_theta, false);
                    bool grid_layout = true;
                    for (std::size_t r = 0; r < n_phi && grid_layout; r++)
                    {
                        for (std::size_t j = 0; j < n_theta && grid_layout; j++)
                        {
                            const std::size_t c = r * n_theta + j;
                            const auto& pm = placement[c].phi;
                            const auto& tm = placement[c].theta;
                            if (pm[0] == pm[1] || tm[0] == tm[1]
                                    || !approx_equal(phis_scattered[c], phis_scattered[r * n_theta]))
                            {
                                grid_layout = false;
                                break;
                            }
                            const std::size_t cell = pm[0] * n_theta + tm[0];
                            if (taken[cell])
                                grid_layout = false;
                            taken[cell] = true;
                            mapped[c] = cell;
                        }
                    }
                    if (grid_layout)
                        cell_of = mapped;
                    else if (i == 0)
                        log("Sphere channels are not a full row-major grid; reading amplitude as "
                            "gridded over the sorted angles.");
                }
                cached_cell_of_valid = true;
            }
            // Grid-shaped coverage, laid out as amplitude is, so the same cell_of
            // permutation applies.  A cell takes the share of its own A-weighted
            // band energy, at the sphere, that the builder measured: a band the
            // SNR gate dropped counts by what it carries, not as the whole cell.
            // Counting any unmeasured band as the whole cell put B407 down as
            // 70-80% unmeasured on the demo grid, for 5-10 kHz bands the gate
            // drops in about half its measured cells; by energy 99.6% of its
            // lower hemisphere is measured.  A cell with no energy in any band
            // contributes nothing to a level and takes all or nothing.
            if (group_has_coverage && !channel_coverage)
            {
                if (coverage_scattered.size() != n_phi * n_theta * n_bands)
                    input_error("Database error: sphere ", i, " has a coverage variable that matches "
                            "neither its channels nor its spectral grid.");
                if (coverage_gains.size() != n_bands)
                    coverage_gains = a_weighting_gains(frequencies_, BandLimits {});
                for (std::size_t slot = 0; slot < n_phi * n_theta; slot++)
                {
                    double total = 0.0;
                    double measured = 0.0;
                    bool every_band = true;
                    for (std::size_t band = 0; band < n_bands; band++)
                    {
                        const bool band_measured = coverage_scattered[slot * n_bands + band] != 0;
                        every_band = every_band && band_measured;
                        const double a = sphere_amplitudes[slot * n_bands + band];
                        if (!std::isfinite(a))
                            continue;
                        const double energy = coverage_gains[band] * std::pow(10.0, 0.1 * a);
                        total += energy;
                        if (band_measured)
                            measured += energy;
                    }
                    measured_[cell_of[slot] * n_conditions + i] = static_cast<float>(
                            total > 0.0 ? measured / total : (every_band ? 1.0 : 0.0));
                }
                groups_with_coverage++;
            }
            for (std::size_t slot = 0; slot < n_phi * n_theta; slot++)
            {
                for (std::size_t band = 0; band < n_bands; band++)
                {
                    // A band ambient gating masked arrives as -inf (or NaN):
                    // floor it, exactly as dBA is floored above, so the blends
                    // can tell "not measured" and no NaN reaches a level.
                    const double a = sphere_amplitudes[slot * n_bands + band];
                    spectra_[(cell_of[slot] * n_bands + band) * n_conditions_total + i]
                            = std::isfinite(a) ? a : no_energy_floor_db;
                }
            }
        }
        else if (!frequencies_.empty())
        {
            input_error("Database error: spectra present on some spheres but not others.");
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
            const auto& phi_match = placement[c].phi;
            if (phi_match[0] == phi_match[1])
                continue;
            const auto& theta_match = placement[c].theta;
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
                    // Load sphere data into grid.  dBA of -inf legitimately
                    // means "no energy in this direction" (ambient gating
                    // emptied it) and blend_conditions() must still be able
                    // to read it under a zero barycentric weight -- which a
                    // sparse, coarse triangulation hits routinely, e.g. any
                    // query that lands exactly on a condition's own gamma or
                    // mu.  IEEE 754 makes 0 * -inf a NaN that poisons the
                    // whole blended sum regardless of weight, so floor it to
                    // a level far below anything physical; panam already
                    // established that swapping -inf for a very negative
                    // finite floor changes nothing where the weight is
                    // nonzero (byte-identical footprints), so this is free
                    // everywhere it matters and fixes it everywhere it does.
                    const double dBA =
                            std::isfinite(dBAs_scattered[c]) ? dBAs_scattered[c] : no_energy_floor_db;
                    const double EAA = std::isfinite(EAAs_scattered[c]) ? EAAs_scattered[c] : 0.0;
                    grid[(m * thetas.size() + n) * n_conditions + i] = {dBA, EAA};
                    reached_cells[(m * thetas.size() + n) * n_conditions + i] = 1;
                    if (channel_coverage)
                        measured_[(m * thetas.size() + n) * n_conditions + i] =
                                coverage_scattered[c] != 0 ? 1.0F : 0.0F;
                }
            }
        }
    }
    if (layout_warnings > 3)
        warn("Database ", datafile_path, ": ", layout_warnings - 3, " more sphere layouts leave cells unfilled or drop channels.");
    info("ShaftPower ", shaft_power_.empty() ? "is not in this database." : "is in this database.");
    if (repaired_spheres > 0)
    {
        log("Repaired legacy hover directivity on ", repaired_spheres, " of ", group_ids.size(),
            " spheres.  Regenerate this database with a versioned builder to avoid the fix-up.");
    }
    if (doppler_removed == 1)
        input_error("Database ", datafile_path, " has DOPPLER_SHIFT_REMOVED = 1: its spheres are "
                "stationary-frame, and NICEOPS does not yet shift them to the observer's frequencies. "
                "Propagating them unshifted would be wrong by up to about 2 dB in tonal bands; refused.");
    has_coverage_ = !group_ids.empty() && groups_with_coverage == group_ids.size();
    if (!has_coverage_)
    {
        measured_.clear();
        if (any_coverage_variable)
            warn("Database ", datafile_path, " has a coverage variable that is not used: it must "
                    "be in every group, along its channels or along its spectral grid.");
    }
    upper_hemisphere_measured_ = !group_ids.empty() && groups_with_upper == group_ids.size();
    if (!upper_hemisphere_measured_)
        info("Not every sphere stores its upper half: directions above the plane read the mirror image.");
    max_level_ = dba_max;
    min_excess_factor_ = eaa_min;
    conditions_.resize(advance_ratios.size());
    for (std::size_t i = 0; i < conditions_.size(); i++)
        conditions_[i] = {advance_ratios[i], flight_path_angles[i], thrust_coefficients[i]};
    screen_suspect_cells(datafile_path, load, reached_cells);
    info("Max level found in database is: ", max_level_);
    info("Lowest excess attenuation factor is: ", min_excess_factor_);
    // Set sphere_radius to dimensional_radius
    sphere_radius_ = dimensional_radius;
    inv_sphere_radius_ = 1.0 / sphere_radius_;
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

    // A fixed_load_factor database's CT is the same constant on every stored
    // sphere by construction (panam writes only the LF=1 reference), which
    // min/max_element already confirms: unlike same_grid, this is not
    // trusted blindly, because the whole point of the flag is that every
    // condition's thrust_coefficient is discarded from the triangulation --
    // if the file lied, load factor scaling would silently stop working
    // rather than fail loudly, so it is checked here instead.
    if (fixed_load_factor_ && !approx_equal(thrust_coefficient_min, thrust_coefficient_max))
    {
        input_error("Database error: fixed_load_factor is set but thrust_coefficient is not "
                "constant across spheres (", thrust_coefficient_min, " to ", thrust_coefficient_max,
                ").");
    }
    reference_thrust_coefficient_ = thrust_coefficient_min;
    active_condition_dims_ = fixed_load_factor_ ? 2 : condition_dims;

    // An axis that does not vary cannot be scaled to [0, 1]: (x - min) / 0 is 0/0,
    // NaN, and a triangulation built on NaN answered every query with one
    // arbitrary sphere, silently (a level-flight-only build gave 114.0 dB at four
    // conditions whose true levels were 100 to 112). The triangulation needs every
    // axis it uses to vary, so a database of several spheres that does not is
    // refused here, naming the axis; a single sphere has nothing to triangulate
    // and takes the nearest-condition path, whatever its axes do.
    if (advance_ratios.size() >= 2)
    {
        const auto require_varying = [&](const char* axis, double low, double high) {
            if (!(high > low))
            {
                input_error("Database ", datafile_path, ": the ", axis, " is ", low,
                        " in all ", advance_ratios.size(), " spheres, so the operating conditions "
                        "cannot be interpolated over it. The triangulation needs the advance "
                        "ratio and the flight path angle to vary (and the thrust coefficient too, "
                        "unless the database is fixed_load_factor).");
            }
        };
        require_varying("advance ratio", advance_ratio_min, advance_ratio_max);
        require_varying("flight path angle", flight_path_angle_min, flight_path_angle_max);
        if (!fixed_load_factor_)
            require_varying("thrust coefficient", thrust_coefficient_min, thrust_coefficient_max);
    }
    // BlendPolicy::aam's search, over the raw conditions (it brackets in speed
    // and flight path angle, not in the normalized space).  AAM has no thrust
    // axis.  On a fixed_load_factor database the law carries C_T and every
    // sphere is in the one layer; otherwise the spheres are grouped by their
    // stored C_T (a database built over load factors stores a few distinct
    // values), and a state searches the layer nearest its own.
    {
        std::vector<std::size_t> order(advance_ratios.size());
        std::iota(order.begin(), order.end(), std::size_t {0});
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return thrust_coefficients[a] < thrust_coefficients[b];
        });
        const auto same_layer = [](double a, double b) {
            return std::fabs(a - b) <= 1e-9 * std::max(std::fabs(a), std::fabs(b));
        };
        std::vector<std::vector<std::size_t>> members;
        std::vector<double> layer_ct;
        for (const std::size_t i : order)
        {
            if (members.empty() || (!fixed_load_factor_ && !same_layer(layer_ct.back(), thrust_coefficients[i])))
            {
                members.emplace_back();
                layer_ct.push_back(thrust_coefficients[i]);
            }
            members.back().push_back(i);
        }
        aam_layers_.clear();
        for (std::size_t l = 0; l < members.size(); l++)
        {
            std::sort(members[l].begin(), members[l].end());
            std::vector<double> mu, gamma;
            for (const std::size_t i : members[l])
            {
                mu.push_back(advance_ratios[i]);
                gamma.push_back(flight_path_angles[i]);
            }
            AamLayer layer {layer_ct[l], {}, members[l].size()};
            layer.selector.build(mu, gamma, members[l]);
            aam_layers_.push_back(std::move(layer));
        }
    }

    // Normalize operating conditions
    for (auto& x : advance_ratios)
        x = unit_interval(x, advance_ratio_min, advance_ratio_max);
    trace("Normalized μ   = ", advance_ratios);
    for (auto& x : flight_path_angles)
        x = unit_interval(x, flight_path_angle_min, flight_path_angle_max);
    trace("Normalized γ   = ", flight_path_angles);
    // A fixed_load_factor database has thrust_coefficient_min == max, so this
    // would divide by zero -- and there is nothing to normalize anyway, since
    // CT plays no part in a triangulation that does not include it.
    if (!fixed_load_factor_)
    {
        for (auto& x : thrust_coefficients)
            x = unit_interval(x, thrust_coefficient_min, thrust_coefficient_max);
        trace("Normalized CT  = ", thrust_coefficients);
    }

    // Pack the normalized axes into the row-major point array the
    // interpolator takes -- (mu, gamma, CT) ordinarily, (mu, gamma) alone for
    // a fixed_load_factor database, where CT carries no information and
    // including it as a third, constant coordinate would make every point
    // coplanar and the 3D triangulation degenerate.  A fourth operating
    // parameter would be another column here and nothing else -- the
    // interpolator's dimension is a constructor argument, not a fixed three.
    if (!coherence_record_.tonal_fraction.empty() && coherence_record_.tonal_fraction.size() != frequencies_.size())
        input_error("The database attribute coherence_tonal_fraction holds ", coherence_record_.tonal_fraction.size(),
                " values; it needs one per band of the database's spectra, ", frequencies_.size(), ".");
    std::vector<double> conditions(advance_ratios.size() * active_condition_dims_);
    for (std::size_t i = 0; i < advance_ratios.size(); i++)
    {
        conditions[i * active_condition_dims_ + 0] = advance_ratios[i];
        conditions[i * active_condition_dims_ + 1] = flight_path_angles[i];
        if (!fixed_load_factor_)
            conditions[i * active_condition_dims_ + 2] = thrust_coefficients[i];
    }
    normalized_conditions_ = conditions;
    // The triangulation touches no netCDF state, and nc_close takes about as long:
    // it tears down the metadata netCDF-4 read for every group at nc_open, about
    // 0.1 s for the shipped file's 1312.  So they run at once.
    auto build = [&] { interpolator_.build(active_condition_dims_, advance_ratios.size(), conditions.data()); };
    std::future<void> triangulation;
    try
    {
        triangulation = std::async(std::launch::async, build);
    }
    catch (const std::system_error&)
    {
        // No thread to be had: do them one after the other.
        guard.close();
        build();
        return;
    }
    guard.close();
    triangulation.get();
}

Database::Blend Database::blend_at(double advance_ratio, double flight_path_angle,
        double thrust_coefficient, const BlendOptions& options) const
{
    // The raw query CT, kept before it is normalized below (or, for a
    // fixed_load_factor database, not used for normalization at all): the
    // analytic load-factor correction compares it against
    // reference_thrust_coefficient_ in physical units, not normalized ones.
    const double raw_thrust_coefficient = thrust_coefficient;
    Blend blend;
    blend.advance_ratio_clamped = advance_ratio < advance_ratio_min || advance_ratio > advance_ratio_max;
    blend.flight_path_angle_clamped = flight_path_angle < flight_path_angle_min || flight_path_angle > flight_path_angle_max;
    blend.thrust_coefficient_clamped = !fixed_load_factor_
            && (thrust_coefficient < thrust_coefficient_min || thrust_coefficient > thrust_coefficient_max);

    if (options.policy == BlendPolicy::aam)
    {
        blend_aam(advance_ratio, flight_path_angle, thrust_coefficient, options, blend);
    }
    else
    {
        // Normalize flight condition and limit range to [0.0,1.0]
        advance_ratio = clamp01(unit_interval(advance_ratio, advance_ratio_min, advance_ratio_max));
        flight_path_angle = clamp01(
                unit_interval(flight_path_angle, flight_path_angle_min, flight_path_angle_max));
        if (!fixed_load_factor_)
        {
            thrust_coefficient = clamp01(
                    unit_interval(thrust_coefficient, thrust_coefficient_min, thrust_coefficient_max));
            trace("Normalized operating condition is: ", "μ  = ", advance_ratio, ", γ  = ", flight_path_angle,
                    ", CT = ", thrust_coefficient);
        }
        // CT is not part of a fixed_load_factor triangulation -- see the member
        // comment on fixed_load_factor_ -- so the query is (mu, gamma) there.
        const double condition[condition_dims] = {advance_ratio, flight_path_angle, thrust_coefficient};

        if (options.policy == BlendPolicy::nearest)
        {
            // The stored condition nearest the state in the space the
            // triangulation is built in, taken whole.  The padding repeats it,
            // rather than pointing at sphere 0 as a degenerate barycentric
            // blend does, so the bounds group it with itself alone.
            const std::size_t count = normalized_conditions_.size() / active_condition_dims_;
            const std::size_t nearest = nearest_point(active_condition_dims_, count,
                    normalized_conditions_.data(), condition);
            blend.weights.fill(0.0);
            blend.weights[0] = 1.0;
            blend.indices.fill(nearest);
        }
        else
        {
            interpolator_.barycentric(condition, blend.weights.data(), blend.indices.data());
            if (fixed_load_factor_)
            {
                // The interpolator only ever sees (mu, gamma) and writes 3 of
                // Blend's 4 slots; the 4th is padding, weight 0, matching the
                // existing convention for a degenerate blend (see cull_range's
                // comment on "Skip the padding").
                blend.weights[3] = 0.0;
                blend.indices[3] = blend.indices[0];
            }
        }
        trace("Barycentric weights: ", blend.weights, " and indicies ", blend.indices);
    }

    if (fixed_load_factor_)
    {
        // One-sided, the correction only ever adds level: panam's own
        // convention leaves load factors at or below 1 g unscaled from the
        // stored reference rather than extrapolating them down (a descent or
        // a light-thrust condition is not assumed to follow the same simple
        // scaling law a maneuvering pull-up does), and the caller holds the
        // top at --max_load_factor.  Two-sided, the same law runs below the
        // reference too, held at min_load_factor so that a state near zero g
        // (or below it) is not taken tens of dB down a logarithm.
        if (options.load_factor_correction == LoadFactorCorrection::two_sided)
        {
            const double floor = options.min_load_factor * reference_thrust_coefficient_;
            blend.level_correction_db = 20.0 * std::log10(
                    std::max(raw_thrust_coefficient, floor) / reference_thrust_coefficient_);
        }
        else
        {
            blend.level_correction_db = raw_thrust_coefficient > reference_thrust_coefficient_
                    ? 20.0 * std::log10(raw_thrust_coefficient / reference_thrust_coefficient_)
                    : 0.0;
        }
        trace("Load factor correction (dB): ", blend.level_correction_db);
    }
    return blend;
}

void Database::blend_aam(double advance_ratio, double flight_path_angle, double thrust_coefficient,
        const BlendOptions& options, Blend& blend) const
{
    // AAM's 20 kt gate in the database's own speed unit, advance ratio: the
    // spheres' mu were formed with the database's tip speed, so that is the
    // one to convert with; the caller's stands in for a database without one.
    const double tip_speed = std::isfinite(tip_speed_) ? tip_speed_ : options.tip_speed;
    if (!(tip_speed > 0.0) || !std::isfinite(tip_speed))
    {
        input_error("--blend aam needs the rotor tip speed to turn its 20 kt speed gate into advance "
                "ratio, and neither the database nor the vehicle gives one.");
    }
    const double speed_gate = options.aam_limits.speed_kt * knot / tip_speed;

    // The thrust layer nearest the state (the only one on a fixed_load_factor
    // database); a state exactly between two takes the higher thrust, as the
    // speed search takes the faster sphere.
    std::size_t layer = 0;
    for (std::size_t l = 1; l < aam_layers_.size(); l++)
    {
        const double here = std::fabs(aam_layers_[l].thrust_coefficient - thrust_coefficient);
        const double best = std::fabs(aam_layers_[layer].thrust_coefficient - thrust_coefficient);
        if (here <= best)
            layer = l;
    }
    std::array<double, 4> weights {};
    std::array<std::size_t, 4> ids {};
    const std::size_t filled = aam_layers_[layer].selector.select(advance_ratio, flight_path_angle,
            speed_gate, options.aam_limits, weights, ids);
    for (std::size_t k = 0; k < blend.weights.size(); k++)
    {
        // Unused slots repeat the first sphere at weight 0, so the bounds see
        // only the spheres actually blended.
        blend.weights[k] = k < filled ? weights[k] : 0.0;
        blend.indices[k] = k < filled ? ids[k] : ids[0];
    }
    trace("AAM weights: ", blend.weights, " and indices ", blend.indices);
}

Database::SourceState Database::source_state(const Vec3& source,
        const Vec3& velocity, const Vec3& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient, const BlendOptions& options) const
{
    return {source, sphere_frame(velocity, normal),
            blend_at(advance_ratio, flight_path_angle, thrust_coefficient, options)};
}

double Database::eval(double range, double phi, double theta, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    return eval(range, phi, theta, blend_at(advance_ratio, flight_path_angle, thrust_coefficient));
}

Database::Sample Database::blend_cell(std::size_t i_phi, std::size_t i_theta,
        const Blend& blend) const
{
    const Pair blended = blend_conditions(
            grid.data() + (i_phi * thetas.size() + i_theta) * n_conditions, blend);
    return {blended[0], blended[1]};
}

double Database::eval(double range, double phi, double theta, const Blend& blend) const
{
    // Spreading and absorption are taken from the sphere's surface outward, and
    // the level inside it is the surface's: spherical spreading extrapolated
    // inward grows as 1/r, which put 28 of the demo's observers up to 4.4 dB above
    // the sphere-surface level and diverged at r = 0 (+inf in the CSV, an abort
    // with the asserts live).  The sphere is the source's own volume, so there is
    // no propagation to account for in it.  Clamped once, so the spreading and the
    // excess attenuation below agree on the range; a no-op beyond the sphere.
    range = std::max(range, sphere_radius_);
    // Find the sphere grid cell containing the desired phi, theta, and where in
    // it they sit.  locate() gives the bracketing indices and the interpolation
    // weight together, so nothing below has to read the knot arrays back.
    const auto phi_cell = phis.locate(phi);
    const auto theta_cell = thetas.locate(theta);
    trace("Found phi (", phi, ") in range ", phis[phi_cell.lo], " -- ", phis[phi_cell.hi]);
    trace("Found theta (", theta, ") in range ", thetas[theta_cell.lo], " -- ",
            thetas[theta_cell.hi]);

    // Barycentric interpolation between conditions at each gridpoint, bilinear
    // interpolation between gridpoints.  Blend the operating conditions at each
    // of the four bracketing grid points, then bilinearly interpolate between
    // them.  Both metrics travel together: they are read from the same cells,
    // so splitting them would touch every cache line twice.
    const auto stencil = cell_stencil(phi_cell, theta_cell, thetas.size());
    const Sample* base = grid.data();
    const Pair corners[4] = {blend_conditions(base + stencil.corner[0] * n_conditions, blend),
            blend_conditions(base + stencil.corner[1] * n_conditions, blend),
            blend_conditions(base + stencil.corner[2] * n_conditions, blend),
            blend_conditions(base + stencil.corner[3] * n_conditions, blend)};
    // The same rule across the four grid corners as across the four conditions
    // (see blend_present): a corner with no measurement drops out rather than
    // dragging its neighbors toward silence, so coverage ends at the edge of
    // the last measured cell instead of fading to nothing across it.
    Pair corner;
    if (is_empty(corners[0][0], stencil.weight[0]) || is_empty(corners[1][0], stencil.weight[1])
            || is_empty(corners[2][0], stencil.weight[2]) || is_empty(corners[3][0], stencil.weight[3]))
    {
        corner = blend_present(corners, stencil.weight);
    }
    else
    {
        corner = corners[0] * stencil.weight[0];
        corner += corners[1] * stencil.weight[1];
        corner += corners[2] * stencil.weight[2];
        corner += corners[3] * stencil.weight[3];
    }
    const auto sphere_level = corner[0];
    const auto excess_factor = corner[1];
    trace("Sphere level (dB): ", sphere_level, "\tExcess Factor (dB/km): ", excess_factor);
    // Apply spherical spreading and excess atmospheric absorption to surface level.
    //
    // EAA is the total A-weighted attenuation over the 1000 m reference distance
    // the database was built at: panam writes it as SPLA - SPLAa, the level at
    // the sphere less the level after 1000 m of ISO 9613-1 absorption, and
    // consumes it as EAA / distance * (slant_range - radius).  This is that
    // formula.
    //
    // The D original multiplied it by a further 4.0 under a "TODO Add ground
    // effects" (c74bab3).  That is not a ground effect: it is a fourfold
    // overstatement of atmospheric absorption, unbounded in range -- 168 dB at
    // 6 km where the physical value is 14 -- where a real ground effect is
    // bounded by the surface impedance and roughly range-independent.  Measured
    // against the model it was standing in for, soft ground costs 2 to 4 dB
    // across the whole footprint while the factor removed 7 dB at 250 m rising
    // to 38 dB at 2 km.  See --ground_effect for the real thing.
    // Both scalings are multiplications by a stored reciprocal rather than
    // divisions: eval runs about 100 million times a run, and a double divide
    // is four times the cost of a multiply on every core this targets.  The
    // reciprocal of 1000 is not exact, so this moves the result by an ULP.
    auto excess_attenuation = excess_factor * (range - sphere_radius()) * 1e-3;
    auto spreading_loss = 20.0 * std::log10(range * inv_sphere_radius_);
    trace("Spreading Loss (dB): ", spreading_loss, "\tExcess Attenuation (dB): ", excess_attenuation);
    assert(std::isfinite(sphere_level));
    assert(std::isfinite(excess_factor));
    assert(std::isfinite(excess_attenuation));
    assert(std::isfinite(spreading_loss));
    // level_correction_db is 0 except on a fixed_load_factor database queried
    // above its reference load factor -- see blend_at().  It is a level, not
    // a loss like the two terms above it, so it adds rather than subtracts.
    auto ground_level = sphere_level - spreading_loss - excess_attenuation + blend.level_correction_db;
    trace("Ground level: ", ground_level);
    return ground_level;
}

Database::PackedSpectra Database::pack_spectra(
        const std::array<std::size_t, condition_dims + 1>& indices, std::size_t first_band,
        std::size_t last_band) const
{
    const std::size_t n_bands = frequencies_.size();
    const std::size_t n_cells = phis.size() * thetas.size();
    constexpr std::size_t n_slots = condition_dims + 1;
    PackedSpectra packed;
    packed.first_band = first_band;
    packed.last_band = std::min(last_band, n_bands);
    packed.values.resize(n_cells * n_slots * n_bands);
    packed.empty_bands.assign(n_cells * n_slots, 0);
    for (std::size_t c = 0; c < n_cells; c++)
    {
        for (std::size_t k = 0; k < n_slots; k++)
        {
            double* bands = packed.values.data() + (c * n_slots + k) * n_bands;
            for (std::size_t band = 0; band < n_bands; band++)
            {
                const double v = spectra_[(c * n_bands + band) * n_conditions + indices[k]];
                bands[band] = v;
                // The mask holds 64 bands; eval_spectrum takes the general path past that.
                if (band >= packed.first_band && band < packed.last_band && band < 64 && v <= empty_level_db)
                    packed.empty_bands[c * n_slots + k] |= std::uint64_t {1} << band;
            }
        }
    }
    return packed;
}

void Database::eval_spectrum(double range, double phi, double theta, const Blend& blend,
        const double* absorption, double* levels, std::size_t first_band,
        std::size_t last_band, const PackedSpectra* packed) const
{
    eval_spectrum(range, range, phi, theta, blend, absorption, levels, first_band, last_band,
            packed);
}

void Database::eval_spectrum(double range, double absorbed_range, double phi, double theta,
        const Blend& blend, const double* absorption, double* levels, std::size_t first_band,
        std::size_t last_band, const PackedSpectra* packed) const
{
    const std::size_t n_bands = frequencies_.size();
    assert(n_bands > 0);
    // Bands the metric excludes are never read; skip them rather than blend and
    // propagate values that weighted_level will multiply by a zero gain.
    if (last_band > n_bands)
        last_band = n_bands;

    // The same cell and the same weights as eval(): the spectral path must land
    // on the identical geometry, or a disagreement between the two could be the
    // interpolation rather than the physics.  Both derive the stencil with
    // cell_stencil(), which is what keeps them identical -- sharing locate()
    // alone would only have given them the same cell.
    const auto stencil = cell_stencil(phis.locate(phi), thetas.locate(theta), thetas.size());

    // Spreading is geometry, so it is the same for every band; absorption is
    // not, which is the whole point of propagating the spectrum.
    // Capped at the sphere's surface as eval() does: no spreading gain and no
    // negative absorption inside the source's own volume.
    const double spreading_loss = 20.0 * std::log10(std::max(range, sphere_radius()) / sphere_radius());
    const double path_length = std::max(absorbed_range - sphere_radius(), 0.0);
    const double* w = blend.weights.data();
    const double* sw = stencil.weight;
    constexpr std::size_t n_slots = condition_dims + 1;

    // Where one (cell, band, k) value is: across spectra_'s condition axis,
    // or in a packed copy's band run.
    auto value = [&](std::size_t cell, std::size_t band, std::size_t k) {
        return packed != nullptr ? packed->values[(cell * n_slots + k) * n_bands + band]
                                 : spectra_[(cell * n_bands + band) * n_conditions + blend.indices[k]];
    };
    // One band the general way, with empty entries dropped (blend_present).
    auto general_band = [&](std::size_t band) {
        double corner_level[4];
        for (std::size_t corner = 0; corner < 4; corner++)
        {
            const std::size_t cell = stencil.corner[corner];
            const double v[4] = {value(cell, band, 0), value(cell, band, 1), value(cell, band, 2),
                    value(cell, band, 3)};
            // Empty (gated) entries drop out, as in eval() -- see blend_present.
            if (is_empty(v[0], w[0]) || is_empty(v[1], w[1]) || is_empty(v[2], w[2])
                    || is_empty(v[3], w[3]))
                corner_level[corner] = blend_present(v, w);
            else
                corner_level[corner] = w[0] * v[0] + w[1] * v[1] + w[2] * v[2] + w[3] * v[3];
        }

        const double sphere_level = (is_empty(corner_level[0], sw[0])
                                            || is_empty(corner_level[1], sw[1])
                                            || is_empty(corner_level[2], sw[2])
                                            || is_empty(corner_level[3], sw[3]))
                ? blend_present(corner_level, sw)
                : sw[0] * corner_level[0] + sw[1] * corner_level[1] + sw[2] * corner_level[2]
                        + sw[3] * corner_level[3];
        // Same uniform-offset correction as eval(), applied identically to
        // every band: panam's load_factors=None scaling shifts the whole
        // spectrum together, not its shape (see blend_at()).
        levels[band] = sphere_level - spreading_loss - absorption[band] * path_length
                + blend.level_correction_db;
        assert(std::isfinite(levels[band]));
    };

    if (packed != nullptr && n_bands <= 64 && first_band >= packed->first_band
            && last_band <= packed->last_band)
    {
        // Every band the plain way first: general_band's non-empty expressions,
        // in the same order, so the same result to the bit wherever nothing is
        // empty.  With the bands contiguous and no branch, this is one loop the
        // compiler vectorizes across bands.
        const double* v[4][n_slots];
        for (std::size_t corner = 0; corner < 4; corner++)
            for (std::size_t k = 0; k < n_slots; k++)
                v[corner][k] = packed->values.data()
                        + (stencil.corner[corner] * n_slots + k) * n_bands;
        const double w0 = w[0], w1 = w[1], w2 = w[2], w3 = w[3];
        const double sw0 = sw[0], sw1 = sw[1], sw2 = sw[2], sw3 = sw[3];
        const double offset = blend.level_correction_db;
        double* __restrict out = levels;
        for (std::size_t band = first_band; band < last_band; band++)
        {
            const double c0 = w0 * v[0][0][band] + w1 * v[0][1][band] + w2 * v[0][2][band]
                    + w3 * v[0][3][band];
            const double c1 = w0 * v[1][0][band] + w1 * v[1][1][band] + w2 * v[1][2][band]
                    + w3 * v[1][3][band];
            const double c2 = w0 * v[2][0][band] + w1 * v[2][1][band] + w2 * v[2][2][band]
                    + w3 * v[2][3][band];
            const double c3 = w0 * v[3][0][band] + w1 * v[3][1][band] + w2 * v[3][2][band]
                    + w3 * v[3][3][band];
            const double sphere_level = sw0 * c0 + sw1 * c1 + sw2 * c2 + sw3 * c3;
            out[band] = sphere_level - spreading_loss - absorption[band] * path_length + offset;
        }

        // Then redo the bands where some weighted entry is empty -- EC130 gates
        // its top bands in up to half its cells, and whole directions in 4% --
        // with the empty entries dropped.  Only entries with a weight count:
        // on a fixed_load_factor database one of the four slots never has any.
        // An empty corner level can only come from empty entries, so these
        // masks also cover the second (stencil) stage.
        std::uint64_t redo = 0;
        for (std::size_t corner = 0; corner < 4; corner++)
        {
            const std::uint64_t* masks = packed->empty_bands.data() + stencil.corner[corner] * n_slots;
            for (std::size_t k = 0; k < n_slots; k++)
                redo |= w[k] != 0.0 ? masks[k] : 0;
        }
        if (redo == 0)
            return;
        // blend_present and is_empty, straight-line: an entry counts if its
        // weight is nonzero and it is above the floor, and the plain blend
        // stands when no weighted entry is empty.  Absent terms enter as
        // v * 0 -- exactly nothing -- so this matches general_band to the bit.
        auto blend4 = [](double a0, double a1, double a2, double a3, double b0, double b1,
                              double b2, double b3) {
            if (!((b0 != 0.0 && a0 <= empty_level_db) || (b1 != 0.0 && a1 <= empty_level_db)
                        || (b2 != 0.0 && a2 <= empty_level_db)
                        || (b3 != 0.0 && a3 <= empty_level_db)))
                return b0 * a0 + b1 * a1 + b2 * a2 + b3 * a3;
            const double m0 = a0 > empty_level_db ? b0 : 0.0;
            const double m1 = a1 > empty_level_db ? b1 : 0.0;
            const double m2 = a2 > empty_level_db ? b2 : 0.0;
            const double m3 = a3 > empty_level_db ? b3 : 0.0;
            double sum = 0.0, total = 0.0;
            sum += a0 * m0;
            total += m0;
            sum += a1 * m1;
            total += m1;
            sum += a2 * m2;
            total += m2;
            sum += a3 * m3;
            total += m3;
            return total > 0.0 ? sum / total : no_energy_floor_db;
        };
        while (redo != 0)
        {
            const auto band = static_cast<std::size_t>(__builtin_ctzll(redo));
            redo &= redo - 1;
            const double c0 = blend4(v[0][0][band], v[0][1][band], v[0][2][band], v[0][3][band],
                    w0, w1, w2, w3);
            const double c1 = blend4(v[1][0][band], v[1][1][band], v[1][2][band], v[1][3][band],
                    w0, w1, w2, w3);
            const double c2 = blend4(v[2][0][band], v[2][1][band], v[2][2][band], v[2][3][band],
                    w0, w1, w2, w3);
            const double c3 = blend4(v[3][0][band], v[3][1][band], v[3][2][band], v[3][3][band],
                    w0, w1, w2, w3);
            const double sphere_level = blend4(c0, c1, c2, c3, sw0, sw1, sw2, sw3);
            out[band] = sphere_level - spreading_loss - absorption[band] * path_length + offset;
        }
        return;
    }

    for (std::size_t band = first_band; band < last_band; band++)
        general_band(band);
}

void Database::spectrum_along(const SourceState& state, const Vec3& direction, double* levels,
        std::size_t first_band, std::size_t last_band, const PackedSpectra* packed) const
{
    const auto angles = lookup_angles(state.frame, direction, !upper_hemisphere_measured_);
    assert(std::isfinite(angles[0]) && std::isfinite(angles[1]));
    // At the sphere's radius nothing spreads and nothing is absorbed, so the
    // coefficients are never read; zeros keep them defined.
    thread_local std::vector<double> none;
    none.assign(frequencies_.size(), 0.0);
    eval_spectrum(sphere_radius_, sphere_radius_, angles[0], angles[1], state.blend, none.data(), levels, first_band,
            last_band, packed);
}

// The body both spectral_level_at_observer and spectral_level_along share,
// forced inline into each.  As an ordinary call from the first -- the one the
// straight-ray observer loop takes -- it cost 4% of the acoustic loop, measured
// against the body written out in place.
#if defined(__GNUC__) || defined(__clang__)
[[gnu::always_inline]]
#endif
inline double Database::spectral_level_impl(const RayPath& path, const SourceState& state,
        const double* absorption, const std::vector<double>& gains, const GroundModel* ground,
        const std::vector<std::complex<double>>* impedance, double sound_speed,
        std::size_t first_band, std::size_t last_band, const PackedSpectra* packed,
        const PlateReceiver* plate, double plate_mic_azimuth, const GroundTable* ground_table,
        const double* extra_loss, const double* ground_weight, double* band_levels) const
{
    // Reused across calls: this runs once per observer per trajectory state,
    // and the band count is fixed for the life of the database.
    thread_local std::vector<double> levels;
    levels.resize(frequencies_.size());

    // Everything the level takes from the path (propagation.hpp).
    const auto angles = lookup_angles(state.frame, path.direction, !upper_hemisphere_measured_);
    assert(std::isfinite(angles[0]) && std::isfinite(angles[1]));

    eval_spectrum(path.range, path.path_length, angles[0], angles[1], state.blend, absorption,
            levels.data(), first_band, last_band, packed);
    if (extra_loss != nullptr)
    {
        for (std::size_t band = first_band; band < std::min(last_band, levels.size()); band++)
            levels[band] -= extra_loss[band];
    }

    if (ground != nullptr)
    {
        // The ground effect is interference between the direct ray and the one
        // reflected off the surface, so it needs the geometry in the vertical
        // plane -- source height and horizontal separation -- not the slant
        // range the rest of the propagation works in.
        //
        // Scratch, like levels above: this is the innermost loop of the whole
        // run and returning a vector here allocated 94 million times.
        thread_local std::vector<double> correction;
        correction.resize(frequencies_.size());
        if (plate != nullptr)
        {
            // The propagation direction, source to microphone, in the plate's
            // frame, where the microphone is offset along +y.  This is the
            // bearing between the two, which a ray in a stratified atmosphere
            // keeps unless a crosswind deflects it.
            const double azimuth = plate_frame_azimuth(path.separation, plate_mic_azimuth);
            plate->power_ratio(path.source_height, path.ground_distance, azimuth,
                    correction.data(), first_band, last_band);
        }
        else if (!std::isnan(path.receiver_height) && path.receiver_height != ground->receiver_height)
        {
            // Over terrain (terrain_ray) the microphone's height above each
            // path's mean ground line is its own, and a GroundTable holds one
            // height, so the model is evaluated exactly at this one.
            GroundModel here = *ground;
            here.receiver_height = path.receiver_height;
            excess_ground_attenuation(here, *impedance, path.source_height,
                    path.ground_distance, sound_speed, frequencies_, correction.data(),
                    first_band, last_band, true);
        }
        else if (ground_table != nullptr)
        {
            ground_table->power_ratio(path.source_height, path.ground_distance,
                    correction.data());
        }
        else
        {
            excess_ground_attenuation(*ground, *impedance, path.source_height,
                    path.ground_distance, sound_speed, frequencies_, correction.data(),
                    first_band, last_band, true);
        }
        if (ground_weight != nullptr)
        {
            // A weight scales the correction in dB, so the power ratio's power.
            for (std::size_t band = first_band; band < std::min(last_band, correction.size()); band++)
            {
                const double w = ground_weight[band];
                if (w == 0.0)
                    correction[band] = 1.0;
                else if (w != 1.0)
                    correction[band] = std::pow(correction[band], w);
            }
        }
        if (band_levels != nullptr)
        {
            for (std::size_t band = first_band; band < std::min(last_band, levels.size()); band++)
                band_levels[band] = levels[band] + 10.0 * std::log10(correction[band]);
        }
        // As power ratios multiplied into the band sum: the same as adding the
        // correction in dB, without a log and an exp per band.  Only bands in
        // [first_band, last_band) have non-zero gains, so only they are read.
        return weighted_level(levels.data(), gains, correction.data());
    }

    if (band_levels != nullptr)
    {
        for (std::size_t band = first_band; band < std::min(last_band, levels.size()); band++)
            band_levels[band] = levels[band];
    }
    return weighted_level(levels.data(), gains);
}

double Database::spectral_level_at_observer(const Vec3& observer, const SourceState& state,
        const double* absorption, const std::vector<double>& gains, const GroundModel* ground,
        const std::vector<std::complex<double>>* impedance, double sound_speed,
        std::size_t first_band, std::size_t last_band, const PackedSpectra* packed,
        const PlateReceiver* plate, double plate_mic_azimuth, const GroundTable* ground_table,
        double* band_levels) const
{
    return spectral_level_impl(straight_ray(state.position, observer), state, absorption, gains,
            ground, impedance, sound_speed, first_band, last_band, packed, plate,
            plate_mic_azimuth, ground_table, nullptr, nullptr, band_levels);
}

double Database::spectral_level_along(const RayPath& path, const SourceState& state,
        const double* absorption, const std::vector<double>& gains, const GroundModel* ground,
        const std::vector<std::complex<double>>* impedance, double sound_speed,
        std::size_t first_band, std::size_t last_band, const PackedSpectra* packed,
        const PlateReceiver* plate, double plate_mic_azimuth, const GroundTable* ground_table,
        const double* extra_loss, const double* ground_weight, double* band_levels) const
{
    return spectral_level_impl(path, state, absorption, gains, ground, impedance, sound_speed,
            first_band, last_band, packed, plate, plate_mic_azimuth, ground_table, extra_loss,
            ground_weight, band_levels);
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
        double thrust_coefficient, const BlendOptions& options) const
{
    Sphere sphere;
    sphere.advance_ratio = advance_ratio;
    sphere.flight_path_angle = flight_path_angle;
    sphere.thrust_coefficient = thrust_coefficient;
    sphere.radius = sphere_radius_;
    sphere.speed_reference = speed_reference_;
    sphere.azimuth_reference = azimuth_reference_;
    sphere.tip_speed = tip_speed_;
    sphere.coherence = coherence_record_;
    sphere.rotor_rpm = rotor_rpm_;
    sphere.rotor_radius = rotor_radius_;
    sphere.vehicle_weight_newtons = vehicle_weight_newtons_;
    sphere.phis = phis.values();
    sphere.thetas = thetas.values();
    // blend_at() clamps the condition into the database's envelope, so a request
    // outside it (a 3 for a 0.3) gets the edge sphere while the file says what
    // was asked for.  Say so.
    {
        const bool inside = advance_ratio >= advance_ratio_min && advance_ratio <= advance_ratio_max
                && flight_path_angle >= flight_path_angle_min && flight_path_angle <= flight_path_angle_max
                && (fixed_load_factor_ || (thrust_coefficient >= thrust_coefficient_min
                        && thrust_coefficient <= thrust_coefficient_max));
        if (!inside)
        {
            warn("The requested condition (mu ", advance_ratio, ", gamma ", flight_path_angle, ", CT ",
                    thrust_coefficient, ") lies outside the database's envelope (mu ", advance_ratio_min, " to ",
                    advance_ratio_max, ", gamma ", flight_path_angle_min, " to ", flight_path_angle_max,
                    fixed_load_factor_ ? "" : ", CT ", fixed_load_factor_ ? "" : std::to_string(thrust_coefficient_min),
                    fixed_load_factor_ ? "" : " to ", fixed_load_factor_ ? "" : std::to_string(thrust_coefficient_max),
                    "): the sphere is the nearest one in it, labeled with the condition asked for.");
        }
    }

    // Same blend eval() applies at one direction, carried out at every cell, so
    // a sphere exported at an operating condition the database already holds
    // comes back bit-for-bit as that sphere's own data.  Resolving it is a
    // Delaunay locate walk, so it is done once and shared with the power lookup.
    const auto blend = blend_at(advance_ratio, flight_path_angle, thrust_coefficient, options);
    sphere.shaft_power = shaft_power(blend);
    sphere.cells.resize(phis.size() * thetas.size());
    for (std::size_t m = 0; m < phis.size(); m++)
    {
        for (std::size_t n = 0; n < thetas.size(); n++)
        {
            sphere.cells[m * thetas.size() + n] = blend_cell(m, n, blend);
            // Same correction as eval(): 0 unless this is a fixed_load_factor
            // database queried above its reference load factor.  EAA is an
            // attenuation-rate factor, not a level, so it is untouched.
            sphere.cells[m * thetas.size() + n].dBA += blend.level_correction_db;
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
                const double* conditions = spectra_.data() + (cell * n_bands + band) * n_conditions;
                const double v[4] = {conditions[blend.indices[0]], conditions[blend.indices[1]],
                        conditions[blend.indices[2]], conditions[blend.indices[3]]};
                double value;
                if (is_empty(v[0], w[0]) || is_empty(v[1], w[1]) || is_empty(v[2], w[2])
                        || is_empty(v[3], w[3]))
                {
                    // Empty entries drop out, as in eval() -- see blend_present.
                    // An entry empty in every condition stays at the floor,
                    // which AAM readers take as masked.
                    value = blend_present(v, w.data());
                    if (value <= empty_level_db)
                    {
                        sphere.amplitudes[cell * n_bands + band] = no_energy_floor_db;
                        continue;
                    }
                }
                else
                {
                    value = 0.0;
                    for (std::size_t k = 0; k < 4; k++)
                        value += w[k] * v[k];
                }
                sphere.amplitudes[cell * n_bands + band] = value + blend.level_correction_db;
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
    // If this database's excess attenuation can be a gain, decline to bound.
    if (!excess_attenuation_is_a_loss())
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
    // blend.level_correction_db is 0 except on a fixed_load_factor database
    // queried away from its reference load factor, where it shifts every level
    // of the blend by the same number of dB (see blend_at()) -- so adding it
    // here keeps this a valid upper bound rather than an approximate one,
    // whichever sign it has.
    blended_max += blend.level_correction_db;

    // Spherical spreading alone takes that peak below the cutoff at this range,
    // and excess attenuation only takes it further down.
    return sphere_radius_ * std::pow(10.0, 0.05 * (blended_max - cutoff));
}

void Database::group_blends(const std::vector<SourceState>& states, DirectionBounds& bounds,
        std::vector<std::array<std::size_t, condition_dims + 1>>& blends) const
{
    // Group the states by their blend.  The tetrahedra a trajectory visits are
    // far fewer than its states, and a table per tetrahedron is what keeps this
    // small enough to stay in cache.
    //
    // On a fixed_load_factor database the triangulation is 2D (mu, gamma only
    // -- see blend_at()), so two states can share indices while querying
    // different load factors.  Their correction is a uniform shift in dB, so
    // it goes into a per-state scale rather than the key: keyed on it, every
    // state was its own group (it varies continuously with CT), and building
    // the tables took a quarter of a broadband run.
    bounds.blend_of_state_.assign(states.size(), 0);
    bounds.state_scale_.assign(states.size(), 1.0);
    blends.clear();
    std::size_t last = 0;
    for (std::size_t i = 0; i < states.size(); i++)
    {
        const auto& indices = states[i].blend.indices;
        // Consecutive states nearly always share a tetrahedron.
        std::size_t found = blends.size();
        if (last < blends.size() && blends[last] == indices)
            found = last;
        for (std::size_t b = 0; found == blends.size() && b < blends.size(); b++)
        {
            if (blends[b] == indices)
                found = b;
        }
        if (found == blends.size())
            blends.push_back(indices);
        bounds.blend_of_state_[i] = last = found;
        const double correction = states[i].blend.level_correction_db;
        if (correction != 0.0)
            bounds.state_scale_[i] = std::pow(10.0, 0.1 * correction);
    }
}

int Database::home_bin(std::size_t i_phi, std::size_t i_theta) const
{
    const double theta = thetas[i_theta] * pi / 180.0;
    const double phi = phis[i_phi] * pi / 180.0;
    return DirectionBounds::bin(std::cos(theta), std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi), 1.0);
}

std::array<int, 2> Database::home_bins(std::size_t i_phi, std::size_t i_theta) const
{
    const int here = home_bin(i_phi, i_theta);
    // sin(180 degrees) is 1.2e-16, not 0, so the tail pole's row kept a
    // meaningful azimuth by accident; a tolerance treats both poles the same
    // rather than relying on that.
    if (std::fabs(std::sin(thetas[i_theta] * pi / 180.0)) < 1e-9)
        return {(here / DirectionBounds::azimuth_bins) * DirectionBounds::azimuth_bins,
                DirectionBounds::azimuth_bins};
    return {here, 1};
}

std::array<int, 2> Database::bin_padding() const
{
    // How many bins a single grid cell can straddle.  The bins are coarser than
    // the grid on every shipped database, but widening by the computed number
    // rather than by one keeps the bound valid whatever grid a database brings.
    return {static_cast<int>(std::ceil(
                    grid_spacing(thetas) / DirectionBounds::narrowest_theta_bin_degrees())) + 1,
            static_cast<int>(std::ceil(
                    grid_spacing(phis) / DirectionBounds::azimuth_bin_degrees())) + 1};
}

// A screen for entries no sound database holds, run once the grid and spectra
// are filled.  Two kinds of fault reach this far as finite levels the blends
// take as real: a dBA cell far above or below the rest of its sphere (the S-76D
// M4 and fixed_load_factor files' 150-1450 dBA channels, the September files'
// near-silent cells, B407_M3's written zeros), and a band far above the rest
// of its own cell's spectrum (NASA's R44 spheres' 135-170 dB at 8-10 kHz).
//
// - A dBA cell is compared with the median of its sphere's measured cells
//   (unreached cells, which hold the grid's initial 0, left out).  Comparing
//   it with its grid neighbors instead flagged sound cells along a strip of
//   zeros, whose neighbors were half zeros, and could not see a strip wider
//   than half a neighborhood, which is its own median.  A sphere more than
//   half defective has a defective median, and is beyond any such screen.
// - A band is compared with the median of its cell's non-empty bands, and
//   needs five.  Comparing a band with the same band in neighboring cells
//   does not separate: sound spectra jump 100 dB and more between cells where
//   gating thins a band.  Only an excess is suspect; a band far below its
//   cell's median is what gating leaves.
//
// It only warns unless asked to mask, so the default path is unchanged.
// Masking floors an entry to no_energy_floor_db, which the blends read as
// gated (a dBA cell with every band of its spectrum), and so only removes
// energy: the bounds, built from the grid on
// demand, stay valid.
void Database::screen_suspect_cells(const std::string& datafile_path, const DatabaseLoadOptions& load,
        const std::vector<char>& reached)
{
    struct Suspect
    {
        std::size_t condition, cell;
        std::size_t band; // n_bands for a dBA cell
        double level, reference;
    };
    const std::size_t n_phi = phis.size(), n_theta = thetas.size();
    const std::size_t n_bands = frequencies_.size();
    const std::size_t n_cond = conditions_.size();
    assert(reached.size() == grid.size());
    // Each condition is screened on its own, so they are spread over a few
    // threads (the band medians are most of a load's screen time) and joined
    // in condition order, which the warning and the CSV list them in.
    std::vector<std::vector<Suspect>> per_condition(n_cond);
    const auto screen_condition = [&](std::size_t i, std::vector<double>& around) {
        auto& found = per_condition[i];
        around.clear();
        for (std::size_t cell = 0; cell < n_phi * n_theta; cell++)
        {
            const std::size_t at = cell * n_cond + i;
            if (reached[at] != 0 && grid[at].dBA > empty_level_db)
                around.push_back(grid[at].dBA);
        }
        if (around.size() >= 3)
        {
            const double reference = median_of(around);
            for (std::size_t cell = 0; cell < n_phi * n_theta; cell++)
            {
                const std::size_t at = cell * n_cond + i;
                if (reached[at] != 0 && grid[at].dBA > empty_level_db
                        && std::fabs(grid[at].dBA - reference) > suspect_cell_margin_db)
                    found.push_back({i, cell, n_bands, grid[at].dBA, reference});
            }
        }
        if (n_bands < 5)
            return;
        for (std::size_t cell = 0; cell < n_phi * n_theta; cell++)
        {
            around.clear();
            for (std::size_t band = 0; band < n_bands; band++)
            {
                const double v = spectra_[(cell * n_bands + band) * n_cond + i];
                if (v > empty_level_db)
                    around.push_back(v);
            }
            if (around.size() < 5)
                continue;
            const double reference = median_of(around);
            for (std::size_t band = 0; band < n_bands; band++)
            {
                const double v = spectra_[(cell * n_bands + band) * n_cond + i];
                if (v > empty_level_db && v - reference > suspect_band_margin_db)
                    found.push_back({i, cell, band, v, reference});
            }
        }
    };
    const std::size_t workers = std::max<std::size_t>(1,
            std::min<std::size_t>({n_cond, 8, std::thread::hardware_concurrency()}));
    std::vector<std::future<void>> jobs;
    for (std::size_t w = 0; w < workers; w++)
    {
        jobs.push_back(std::async(std::launch::async, [&, w] {
            std::vector<double> around;
            around.reserve(std::max<std::size_t>(8, n_bands));
            for (std::size_t i = w; i < n_cond; i += workers)
                screen_condition(i, around);
        }));
    }
    for (auto& job : jobs)
        job.get();
    std::vector<Suspect> found;
    for (const auto& list : per_condition)
        found.insert(found.end(), list.begin(), list.end());

    suspect_cells_ = static_cast<std::size_t>(std::count_if(found.begin(), found.end(),
            [&](const Suspect& s) { return s.band == n_bands; }));
    suspect_bands_ = found.size() - suspect_cells_;

    const auto describe = [&](const Suspect& s) {
        std::ostringstream text;
        text << "sphere " << s.condition << " (mu " << conditions_[s.condition][0] << ", gamma "
             << conditions_[s.condition][1] << ", C_T " << conditions_[s.condition][2]
             << ") at phi " << phis[s.cell / n_theta] << ", theta " << thetas[s.cell % n_theta];
        if (s.band < n_bands)
            text << ", " << frequencies_[s.band] << " Hz";
        text << ": " << s.level << " dB against " << s.reference;
        return text.str();
    };
    const auto report = [&](bool bands, const char* one, const char* many, double margin) {
        std::size_t count = 0;
        std::vector<bool> sphere_hit(n_cond, false);
        const Suspect* worst = nullptr;
        for (const auto& s : found)
        {
            if ((s.band < n_bands) != bands)
                continue;
            count++;
            sphere_hit[s.condition] = true;
            if (!worst || std::fabs(s.level - s.reference) > std::fabs(worst->level - worst->reference))
                worst = &s;
        }
        if (count == 0)
            return;
        const auto spheres = std::count(sphere_hit.begin(), sphere_hit.end(), true);
        warn("Database ", datafile_path, ": ", count, " ", count == 1 ? one : many, " in ", spheres,
                " of ", n_cond, " spheres ", count == 1 ? "lies" : "lie", " more than ", margin, " dB from ",
                count == 1 ? "the median of its " : "the median of their ",
                bands ? "cell's spectrum" : "sphere",
                ", which no sound database does; the worst is ", describe(*worst), ". ",
                load.mask_suspect_cells ? "Dropped (--mask_suspect_cells)."
                                        : "Blended as real levels; --screen_database lists them "
                                          "and --mask_suspect_cells drops them.");
    };
    report(false, "dBA cell", "dBA cells", suspect_cell_margin_db);
    report(true, "spectral band", "spectral bands", suspect_band_margin_db);

    if (!load.screen_file.empty())
    {
        ensure_parent_directory(load.screen_file);
        std::ofstream out(load.screen_file);
        if (!out)
            input_error("Cannot write --screen_database file ", load.screen_file, ".");
        out << "# Suspect entries in " << datafile_path << ": dBA cells more than "
            << suspect_cell_margin_db << " dB from the median of their sphere, and bands more than "
            << suspect_band_margin_db << " dB above the median of their cell's spectrum.\n";
        out << "sphere,advance_ratio,flight_path_angle,thrust_coefficient,phi,theta,band_hz,"
               "level_db,reference_db,masked\n";
        out.precision(10);
        for (const auto& s : found)
        {
            out << s.condition << ',' << conditions_[s.condition][0] << ','
                << conditions_[s.condition][1] << ',' << conditions_[s.condition][2] << ','
                << phis[s.cell / n_theta] << ',' << thetas[s.cell % n_theta] << ',';
            if (s.band < n_bands)
                out << frequencies_[s.band];
            out << ',' << s.level << ',' << s.reference << ',' << (load.mask_suspect_cells ? 1 : 0)
                << '\n';
        }
        if (!out.flush())
            input_error("Could not finish writing --screen_database file ", load.screen_file, ".");
        log("Wrote ", found.size(), found.size() == 1 ? " suspect entry to " : " suspect entries to ",
                load.screen_file, ".");
    }

    if (!load.mask_suspect_cells || found.empty())
        return;
    for (const auto& s : found)
    {
        if (s.band < n_bands)
        {
            spectra_[(s.cell * n_bands + s.band) * n_cond + s.condition] = no_energy_floor_db;
            continue;
        }
        // The spectral path reads only the bands, so a dropped cell drops them
        // too: a spectrum raised or lowered as a whole has no band out of line
        // with its own cell, and would otherwise still be blended.
        grid[s.cell * n_cond + s.condition].dBA = no_energy_floor_db;
        for (std::size_t band = 0; band < n_bands; band++)
            spectra_[(s.cell * n_bands + band) * n_cond + s.condition] = no_energy_floor_db;
    }
    if (suspect_cells_ == 0)
        return;
    // The per-sphere maxima bound a blended state's level (cull_range,
    // condition_bounds); a masked cell must leave them, or a 1450 dBA spike
    // keeps every observer in range.
    max_level_ = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n_cond; i++)
    {
        sphere_max_[i] = -std::numeric_limits<double>::infinity();
        for (std::size_t cell = 0; cell < n_phi * n_theta; cell++)
        {
            const double v = grid[cell * n_cond + i].dBA;
            if (v > empty_level_db)
                sphere_max_[i] = std::max(sphere_max_[i], v);
        }
        max_level_ = std::max(max_level_, sphere_max_[i]);
    }
}

DirectionBounds Database::direction_bounds(const std::vector<SourceState>& states) const
{
    DirectionBounds bounds;
    // The bound drops the excess attenuation term, so it rests on the same
    // condition cull_range does.  Returning early leaves `bounds` unusable,
    // which is what the caller checks; nothing below it is worth allocating.
    if (!excess_attenuation_is_a_loss() || states.empty())
        return bounds;
    std::vector<std::array<std::size_t, condition_dims + 1>> blends;
    group_blends(states, bounds, blends);

    const auto [theta_pad, azimuth_pad] = bin_padding();
    const std::size_t n_thetas = thetas.size();
    bounds.numerators_.assign(blends.size() * DirectionBounds::bin_count,
            -std::numeric_limits<double>::infinity());

    for (std::size_t b = 0; b < blends.size(); b++)
    {
        double* row = bounds.numerators_.data() + b * DirectionBounds::bin_count;
        for (std::size_t m = 0; m < phis.size(); m++)
        {
            for (std::size_t n = 0; n < n_thetas; n++)
            {
                // The largest level the four blended conditions can give at
                // this cell.  Taking the maximum over the conditions rather
                // than any particular convex combination of them is what makes
                // one table serve every state on the tetrahedron.
                const Sample* cell = grid.data() + (m * n_thetas + n) * n_conditions;
                double top = -std::numeric_limits<double>::infinity();
                for (std::size_t k = 0; k < blends[b].size(); k++)
                    top = std::max(top, cell[blends[b][k]].dBA);

                const auto [first_home, home_count] = home_bins(m, n);
                for (int home_index = 0; home_index < home_count; home_index++)
                {
                    const int here = first_home + home_index;
                    const int theta_index = here / DirectionBounds::azimuth_bins;
                    const int azimuth_index = here % DirectionBounds::azimuth_bins;
                    // A direction reads this cell through a bilinear stencil that
                    // reaches a grid step either way, so every bin that stencil
                    // can be queried from has to carry the cell's value.
                    for (int dt = -theta_pad; dt <= theta_pad; dt++)
                    {
                        const int t = std::min(DirectionBounds::theta_bins - 1,
                                std::max(0, theta_index + dt));
                        for (int da = -azimuth_pad; da <= azimuth_pad; da++)
                        {
                            const int a = std::min(DirectionBounds::azimuth_bins - 1,
                                    std::max(0, azimuth_index + da));
                            double& slot = row[t * DirectionBounds::azimuth_bins + a];
                            slot = std::max(slot, top);
                        }
                    }
                }
            }
        }
    }

    // Fold the level into the energy domain, with the sphere radius, so the hot
    // path compares numerator >= energy * range^2 and needs no logarithm.
    for (double& value : bounds.numerators_)
        value = std::pow(10.0, 0.1 * value) * sphere_radius_ * sphere_radius_;
    bounds.usable_ = true;
    return bounds;
}

DirectionBounds Database::spectral_direction_bounds(const std::vector<SourceState>& states,
        const std::vector<double>& gains) const
{
    // What is being bounded: eval_spectrum's per-band level S_b, less
    // spreading and absorption, summed as L = 10 log10 sum_b g_b 10^(S_b/10).
    // Absorption only lowers each band, so dropping it leaves an upper bound,
    // and spreading is the numerator's r^2, as on the broadband path.  That
    // leaves bounding the blended spectrum's weighted level, two ways:
    //
    //   per band -- each blended S_b is a convex combination of the stored
    //   values it reads (empty entries dropped and the rest renormalized, which
    //   is still a convex combination of what remains), so it is at most their
    //   maximum, and L is increasing in every S_b.  Taking that maximum band by
    //   band over the bin can pair a tone from one cell with a tone from
    //   another, which is loose, but it holds whatever is empty.
    //
    //   per spectrum -- where nothing read is empty, the blend is one convex
    //   combination of whole stored spectra with the same weights in every
    //   band, and L, a log-sum-exp, is convex in the spectrum, so it is at most
    //   the largest L of any one stored spectrum.  Much tighter, but it needs
    //   every band the metric reads to be present in every spectrum of the
    //   bin; an empty band renormalizes that band alone and breaks the "one
    //   combination" premise.
    //
    // Each bin takes the smaller of the two it is entitled to.
    DirectionBounds bounds;
    const std::size_t n_bands = frequencies_.size();
    if (states.empty() || n_bands == 0 || gains.size() != n_bands)
        return bounds;
    std::vector<std::size_t> active;
    for (std::size_t band = 0; band < n_bands; band++)
    {
        if (gains[band] > 0.0)
            active.push_back(band);
    }
    if (active.empty())
        return bounds;
    const std::size_t n_active = active.size();

    std::vector<std::array<std::size_t, condition_dims + 1>> blends;
    group_blends(states, bounds, blends);

    const std::size_t n_thetas = thetas.size();
    const std::size_t n_cells = phis.size() * n_thetas;
    std::vector<std::array<int, 2>> home(n_cells);
    for (std::size_t m = 0; m < phis.size(); m++)
        for (std::size_t n = 0; n < n_thetas; n++)
            home[m * n_thetas + n] = home_bins(m, n);

    // Weighted level of each stored spectrum, and whether it has an empty band
    // the metric reads.  Per condition, not per group: groups share conditions,
    // and only the conditions some group uses are ever computed.
    constexpr double ln10_over_10 = 0.23025850929940456840;
    std::vector<double> spectrum_level(n_cells * n_conditions, 0.0);
    std::vector<char> spectrum_has_empty(n_cells * n_conditions, 0);
    std::vector<char> condition_done(n_conditions, 0);
    auto prepare_condition = [&](std::size_t k) {
        if (condition_done[k])
            return;
        condition_done[k] = 1;
        for (std::size_t c = 0; c < n_cells; c++)
        {
            double total = 0.0;
            bool empty = false;
            for (const std::size_t band : active)
            {
                const double v = spectra_[(c * n_bands + band) * n_conditions + k];
                empty = empty || v <= empty_level_db;
                total += gains[band] * std::exp(v * ln10_over_10);
            }
            spectrum_level[c * n_conditions + k] = 10.0 * std::log10(total);
            spectrum_has_empty[c * n_conditions + k] = empty;
        }
    };

    const auto [theta_pad, azimuth_pad] = bin_padding();
    constexpr int bins = DirectionBounds::bin_count;
    const double lowest = -std::numeric_limits<double>::infinity();
    // Per home bin (no padding yet): the band-by-band maximum, the largest
    // single-spectrum level, and whether any spectrum there had an empty band.
    std::vector<double> band_max(bins * n_active);
    std::vector<double> spectrum_max(bins);
    std::vector<char> has_empty(bins);
    std::vector<double> reach(n_active);
    bounds.numerators_.assign(blends.size() * bins, lowest);

    for (std::size_t b = 0; b < blends.size(); b++)
    {
        for (const std::size_t k : blends[b])
            prepare_condition(k);
        std::fill(band_max.begin(), band_max.end(), lowest);
        std::fill(spectrum_max.begin(), spectrum_max.end(), lowest);
        std::fill(has_empty.begin(), has_empty.end(), 0);
        for (std::size_t c = 0; c < n_cells; c++)
        {
            // One bin, or at a pole all the azimuth bins of its theta bin.
            for (int home_index = 0; home_index < home[c][1]; home_index++)
            {
                const int h = home[c][0] + home_index;
                double* maxima = band_max.data() + h * n_active;
                for (const std::size_t k : blends[b])
                {
                    spectrum_max[h] = std::max(spectrum_max[h], spectrum_level[c * n_conditions + k]);
                    has_empty[h] = has_empty[h] || spectrum_has_empty[c * n_conditions + k];
                    for (std::size_t j = 0; j < n_active; j++)
                    {
                        // Raw values, the floor included: a band empty in every
                        // entry blends to the floor, not to nothing.
                        maxima[j] = std::max(maxima[j],
                                spectra_[(c * n_bands + active[j]) * n_conditions + k]);
                    }
                }
            }
        }

        // Spread each home bin over the bins its stencil can be queried from --
        // the same reach direction_bounds() gives each cell, applied per bin.
        double* row = bounds.numerators_.data() + b * bins;
        for (int t = 0; t < DirectionBounds::theta_bins; t++)
        {
            for (int a = 0; a < DirectionBounds::azimuth_bins; a++)
            {
                std::fill(reach.begin(), reach.end(), lowest);
                double whole = lowest;
                bool empty = false;
                for (int st = std::max(0, t - theta_pad);
                        st <= std::min(DirectionBounds::theta_bins - 1, t + theta_pad); st++)
                {
                    for (int sa = std::max(0, a - azimuth_pad);
                            sa <= std::min(DirectionBounds::azimuth_bins - 1, a + azimuth_pad);
                            sa++)
                    {
                        const int h = st * DirectionBounds::azimuth_bins + sa;
                        whole = std::max(whole, spectrum_max[h]);
                        empty = empty || has_empty[h];
                        const double* maxima = band_max.data() + h * n_active;
                        for (std::size_t j = 0; j < n_active; j++)
                            reach[j] = std::max(reach[j], maxima[j]);
                    }
                }
                double total = 0.0;
                for (std::size_t j = 0; j < n_active; j++)
                    total += gains[active[j]] * std::exp(reach[j] * ln10_over_10);
                double bound = total > 0.0 ? 10.0 * std::log10(total) : lowest;
                if (!empty)
                    bound = std::min(bound, whole);
                row[t * DirectionBounds::azimuth_bins + a] = bound;
            }
        }
    }

    // The energy form, as direction_bounds().  The factor covers the rounding
    // between this log-sum-exp and weighted_level's: the bound can be exact
    // (a cell center, all weight on one condition), and a state then evaluates
    // to it within an ulp or two.
    for (double& value : bounds.numerators_)
        value = std::pow(10.0, 0.1 * value) * sphere_radius_ * sphere_radius_ * (1.0 + 1e-9);
    bounds.usable_ = true;
    return bounds;
}

double Database::unmeasured_share(const std::array<double, 2>& angles, const Blend& blend) const
{
    if (!has_coverage_)
        return 0.0;
    // The four corners and their bilinear weights, exactly as eval() takes
    // them, each read at the blended conditions.  The blend's own weights are
    // normalized so the share is a fraction whatever the rounding in them.
    const auto stencil = cell_stencil(phis.locate(angles[0]), thetas.locate(angles[1]), thetas.size());
    double total = 0.0;
    for (const double weight : blend.weights)
        total += weight;
    if (!(total > 0.0))
        return 0.0;
    double share = 0.0;
    for (std::size_t corner = 0; corner < 4; corner++)
    {
        double unmeasured = 0.0;
        for (std::size_t j = 0; j < blend.weights.size(); j++)
        {
            if (blend.weights[j] != 0.0)
                unmeasured += blend.weights[j]
                        * (1.0 - measured_[stencil.corner[corner] * n_conditions + blend.indices[j]]);
        }
        share += stencil.weight[corner] * unmeasured / total;
    }
    return share;
}

double Database::level_from_geometry(double squared_range, double along, double up,
        double across, const Blend& blend) const
{
    const double range = std::sqrt(squared_range);
    trace("R = ", range);
    // The frame is orthonormal, so `along`, `up` and `across` are the
    // displacement's components in it.  Only theta needs the direction
    // normalized; phi is an angle between two of the components and does not
    // care about their common scale.
    const double inverse = range > 0.0 ? 1.0 / range : 0.0;
    const double theta = acosc(along * inverse) * (180.0 / pi);
    const double phi = phi_degrees(up, across, !upper_hemisphere_measured_);
    trace("Looked up [phi,theta] = ", phi, ", ", theta);
    assert(std::isfinite(phi) && std::isfinite(theta));
    const auto level = eval(range, phi, theta, blend);
    assert(std::isfinite(level));
    return level;
}

double Database::level_at_observer(const Vec3& observer, const SourceState& state) const
{
    const Vec3 delta = observer - state.position;
    return level_from_geometry(delta.squaredNorm(), state.frame.nose.dot(delta),
            state.frame.normal.dot(delta), state.frame.starboard.dot(delta), state.blend);
}

double Database::level_at_observer(const Vec3& observer, const Vec3& source,
        const Vec3& velocity, const Vec3& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    return level_at_observer(observer,
            source_state(source, velocity, normal, advance_ratio, flight_path_angle,
                    thrust_coefficient));
}

double Database::level_along_bearing(double range, double bearing, double altitude,
        const Vec3& velocity, const Vec3& normal, double advance_ratio,
        double flight_path_angle, double thrust_coefficient) const
{
    bearing *= pi / 180.0;
    // Define ground velocity, and get heading vector from it
    double ground_velocity[2] = {velocity(0), velocity(1)};
    double hnorm = std::sqrt(ground_velocity[0] * ground_velocity[0] + ground_velocity[1] * ground_velocity[1]);
    double heading[2] = {ground_velocity[0] / hnorm, ground_velocity[1] / hnorm};
    // Place observer on ground at heading rotated by bearing
    // NOTE: second component faithfully reproduces the D original, which uses
    // heading[0] twice (likely a latent bug in level_along_bearing).
    Vec3 observer = {range * (std::cos(bearing) * heading[0] + std::sin(bearing) * heading[1]),
            range * (std::sin(-bearing) * heading[0] + std::cos(bearing) * heading[0]), 0.0};
    // Put source at altitude above ground
    Vec3 source = {0.0, 0.0, altitude};
    return level_at_observer(observer, source, velocity, normal, advance_ratio, flight_path_angle,
            thrust_coefficient);
}

double Database::zero_range_level(double altitude, const Vec3& velocity,
        const Vec3& normal, double advance_ratio, double flight_path_angle,
        double thrust_coefficient) const
{
    Vec3 observer = {0.0, 0.0, 0.0};
    Vec3 source = {0.0, 0.0, altitude};
    return level_at_observer(observer, source, velocity, normal, advance_ratio, flight_path_angle,
            thrust_coefficient);
}

double Database::find_range(double contour_level, double bearing, double altitude,
        const Vec3& velocity, const Vec3& normal, double advance_ratio,
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
        // Set maximum range bound from max_level assuming spherical spreading
        // only.  max_level_ is the database's raw stored maximum, which does
        // not know about a fixed_load_factor correction (see blend_at()); add
        // it here or a query above the reference load factor could have a
        // true zero-range level past this bound, and the bisection below
        // would search short and miss the real crossing.
        const auto correction = blend_at(advance_ratio, flight_path_angle, thrust_coefficient)
                .level_correction_db;
        auto range_max = sphere_radius() * std::pow(10.0, 0.05 * (max_level_ + correction - contour_level));
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
        input_error("netCDF error while ", what, ": ", nc_strerror(status));
}

int define_scalar(int group, const char* name, int dimid)
{
    int varid = 0;
    nc_check(nc_def_var(group, name, NC_DOUBLE, 1, &dimid, &varid), std::string("defining ") + name);
    return varid;
}


// The record's root attributes, in define mode.
// A root attribute of one number, where there is one; in define mode.
void write_root_number(int file, const char* name, const std::optional<double>& v)
{
    if (v)
        nc_check(nc_put_att_double(file, NC_GLOBAL, name, NC_DOUBLE, 1, &*v), std::string("writing ") + name);
}

void write_coherence_record(int file, const CoherenceRecord& r)
{
    auto number = [&](const char* name, const std::optional<double>& v) { write_root_number(file, name, v); };
    auto text = [&](const char* name, const std::string& v) {
        if (!v.empty())
            nc_check(nc_put_att_text(file, NC_GLOBAL, name, v.size(), v.c_str()), std::string("writing ") + name);
    };
    number("main_rotor_blades", r.blades);
    if (r.rotation_sense)
        text("main_rotor_rotation", *r.rotation_sense > 0 ? "counterclockwise" : "clockwise");
    number("coherence_rotation_limit_Hz", r.rotation_limit);
    number("coherence_scale_Hz_deg", r.scale);
    number("coherence_exponent", r.exponent);
    number("coherence_spread_scale_Hz_deg", r.spread_scale);
    number("coherence_spread_exponent", r.spread_exponent);
    if (!r.tonal_fraction.empty())
        nc_check(nc_put_att_double(file, NC_GLOBAL, "coherence_tonal_fraction", NC_DOUBLE, r.tonal_fraction.size(),
                         r.tonal_fraction.data()),
                "writing coherence_tonal_fraction");
    text("coherence_source", r.source);
}

// The band levels of a level sphere on its own (phi, theta) grid, but in AAM's
// frame: AAM's sphere is fixed to the body, whose pitch is the angle of attack
// plus the flight path angle, so with AOA 0 (the deck AAM's own spheres assume)
// its nose is the velocity.  Each grid direction is taken from that frame,
// pitched nose-up by the sphere's flight path angle, into the level one and read
// with eval's stencil: bilinear in dB, empty corners dropping out.  AAM pitches
// it back by the track's angle, so it hears what NICEOPS hears.  Laid out as
// Sphere::amplitudes; a level condition comes back as it is, to the bit.
std::vector<double> flight_path_frame_amplitudes(const Database::Sphere& sphere)
{
    const double gamma = sphere.flight_path_angle * (pi / 180.0);
    if (gamma == 0.0)
        return sphere.amplitudes;
    GridAxis phi_axis, theta_axis;
    phi_axis.reset(sphere.phis);
    phi_axis.set_period(360.0);
    theta_axis.reset(sphere.thetas);
    // A database that stores no upper half reads its mirror image above the
    // plane, as eval does.
    const bool fold_upper = std::none_of(sphere.phis.begin(), sphere.phis.end(),
            [](double phi) { return std::fabs(phi) > 90.0 + 1e-9; });
    // The level frame in a z-up world: nose along x, starboard nose x normal.
    const SphereFrame level {Vec3(0.0, 0.0, 1.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, -1.0, 0.0)};
    const double cos_gamma = std::cos(gamma), sin_gamma = std::sin(gamma);
    const std::size_t n_theta = sphere.thetas.size(), n_bands = sphere.frequencies.size();
    std::vector<double> rotated(sphere.amplitudes.size());
    for (std::size_t m = 0; m < sphere.phis.size(); m++)
    {
        for (std::size_t n = 0; n < n_theta; n++)
        {
            // The direction in the flight-path frame (forward, starboard, down),
            // then in the level one: its forward axis is (cos, 0, -sin) gamma.
            const double phi = sphere.phis[m] * (pi / 180.0), theta = sphere.thetas[n] * (pi / 180.0);
            const double forward = std::cos(theta), starboard = std::sin(theta) * std::sin(phi),
                         down = std::sin(theta) * std::cos(phi);
            const double level_forward = cos_gamma * forward + sin_gamma * down;
            const double level_down = -sin_gamma * forward + cos_gamma * down;
            const auto angles = lookup_angles(level, Vec3(level_forward, -starboard, -level_down), fold_upper);
            const auto p = phi_axis.locate(angles[0]);
            const auto t = theta_axis.locate(angles[1]);
            const std::size_t corner[4] = {p.lo * n_theta + t.lo, p.lo * n_theta + t.hi,
                    p.hi * n_theta + t.lo, p.hi * n_theta + t.hi};
            const double w[4] = {(1.0 - p.fraction) * (1.0 - t.fraction), (1.0 - p.fraction) * t.fraction,
                    p.fraction * (1.0 - t.fraction), p.fraction * t.fraction};
            for (std::size_t band = 0; band < n_bands; band++)
            {
                double v[4];
                bool any_empty = false;
                for (std::size_t k = 0; k < 4; k++)
                {
                    v[k] = sphere.amplitudes[corner[k] * n_bands + band];
                    any_empty = any_empty || is_empty(v[k], w[k]);
                }
                double value = 0.0;
                if (any_empty)
                {
                    value = blend_present(v, w);
                    if (value <= empty_level_db)
                        value = Database::no_energy_floor_db;
                }
                else
                {
                    for (std::size_t k = 0; k < 4; k++)
                        value += w[k] * v[k];
                }
                rotated[(m * n_theta + n) * n_bands + band] = value;
            }
        }
    }
    return rotated;
}

} // namespace

void write_sphere(const Database::Sphere& sphere, const std::string& path)
{
    ensure_parent_directory(path);

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
    // The rest of the root the loader reads, where the source had it: the
    // sphere reloads with the same speed reference and rotor, and so the same
    // wind requirement, nose and vehicle-file checks, as the database it came from.
    struct RootScalar
    {
        const char* name;
        double value;
        int var;
    };
    RootScalar root_scalars[] = {{"main_rotor_tip_speed_meters_per_sec", sphere.tip_speed, 0},
            {"main_rotor_radius_meters", sphere.rotor_radius, 0},
            {"vehicle_weight_newtons", sphere.vehicle_weight_newtons, 0}};
    for (auto& scalar : root_scalars)
    {
        if (std::isfinite(scalar.value))
            nc_check(nc_def_var(file, scalar.name, NC_DOUBLE, 1, &scalar_dim, &scalar.var),
                    std::string("defining ") + scalar.name);
    }
    const std::string reference = to_string(sphere.speed_reference);
    nc_check(nc_put_att_text(file, NC_GLOBAL, "speed_reference", reference.size(), reference.c_str()),
            "writing speed_reference");
    if (sphere.azimuth_reference)
    {
        const std::string azimuth = to_string(*sphere.azimuth_reference);
        nc_check(nc_put_att_text(file, NC_GLOBAL, "azimuth_reference", azimuth.size(),
                         azimuth.c_str()),
                "writing azimuth_reference");
    }
    const std::string frame = "level";
    nc_check(nc_put_att_text(file, NC_GLOBAL, "sphere_frame", frame.size(), frame.c_str()),
            "writing sphere_frame");

    write_coherence_record(file, sphere.coherence);
    write_root_number(file, "rotor_rpm", sphere.rotor_rpm);

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
    for (const auto& scalar : root_scalars)
    {
        if (scalar.var != 0)
            nc_check(nc_put_var_double(file, scalar.var, &scalar.value), std::string("writing ") + scalar.name);
    }

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

    nc_check(guard.close(), "closing " + path);
    log("Wrote interpolated sphere to ", path, " (", n_channels, " channels).");
}

bool aam_third_octave_bands(const std::vector<double>& frequencies)
{
    // Band n is centered at 10^(n/10) Hz.  The nominal centers (12.5, 31.5, 63, ...)
    // and base-2 ones sit within 0.031 of n in 10 lg f, a 1/12-octave set's off-band
    // members 0.25 or more away; 0.05 lies between.
    if (frequencies.empty())
        return false;
    std::array<bool, 31> seen {};
    for (const double f : frequencies)
    {
        if (!(f > 0.0))
            return false;
        const double index = 10.0 * std::log10(f);
        const long n = std::lround(index);
        if (n < 10 || n > 40 || std::fabs(index - static_cast<double>(n)) > 0.05 || seen[n - 10])
            return false;
        seen[n - 10] = true;
    }
    return true;
}

void write_aam_sphere(const Database::Sphere& sphere, const std::string& path, double tip_speed)
{
    if (sphere.frequencies.empty())
    {
        input_error("Cannot write ", path,
                " in AAM format: the database carries no spectra, and AMPLITUDE(PHI, THETA, "
                "FREQUENCY) is the whole payload of the format.");
    }
    ensure_parent_directory(path);

    const double m2ft = 1.0 / ft2m;
    const double mps2kt = 1.0 / knot;

    // Classic netCDF, as AAM's own spheres and panam's are.
    int file = 0;
    nc_check(nc_create(path.c_str(), NC_CLOBBER, &file), "creating " + path);
    NcGuard guard(file);

    // The AAM/RNM hemisphere layout, as the source spheres themselves carry it:
    // a handful of scalar descriptors, the three axes, and the spectrum.  AAM 3.1
    // reads them by position ("Variable ID mismatch" otherwise), so they are
    // declared in its order, each scalar with a dimension of its own name;
    // MASTTILT, which AAM does not read, comes after AMPLITUDE.  Scalar units are
    // blank-padded to 20 characters, the Fortran field AAM's spheres write.
    auto scalar_dim = [&](const char* name) {
        int dim = 0;
        nc_check(nc_def_dim(file, name, 1, &dim), std::string("defining ") + name);
        return dim;
    };
    auto define_var = [&](const char* name, int dim, const std::string& unit) {
        int var = 0;
        nc_check(nc_def_var(file, name, NC_FLOAT, 1, &dim, &var), std::string("defining ") + name);
        nc_check(nc_put_att_text(file, var, "unit", unit.size(), unit.c_str()),
                std::string("annotating ") + name);
        return var;
    };
    auto scalar_var = [&](const char* name, const char* unit) {
        std::string padded(unit);
        padded.resize(20, ' ');
        return define_var(name, scalar_dim(name), padded);
    };

    const int bb_var = scalar_var("BB", "");
    const int nb_var = scalar_var("NB", "");
    const int pt_var = scalar_var("PT", "");
    const int dsr_var = scalar_var("DOPPLER_SHIFT_REMOVED", "");
    const int empty_var = scalar_var("EMPTY_WEIGHT", "POUNDS");
    const int fuel_var = scalar_var("FUEL_WEIGHT", "POUNDS");
    const int load_var = scalar_var("LOAD_WEIGHT", "POUNDS");
    const int radius_var = scalar_var("RADIUS", "FEET");
    const int fpa_var = scalar_var("FLIGHT_PATH_ANGLE", "DEGREE");
    const int pylon_var = scalar_var("PYLON_ANGLE", "DEGREE");
    const int speed_var = scalar_var("SPEED", "KNOTS");

    int xyz_dim = 0, phi_dim = 0, theta_dim = 0, freq_dim = 0;
    nc_check(nc_def_dim(file, "XYZ", 3, &xyz_dim), "defining XYZ dimension");
    nc_check(nc_def_dim(file, "PHI", sphere.phis.size(), &phi_dim), "defining PHI dimension");
    nc_check(nc_def_dim(file, "THETA", sphere.thetas.size(), &theta_dim),
            "defining THETA dimension");
    nc_check(nc_def_dim(file, "FREQUENCY", sphere.frequencies.size(), &freq_dim),
            "defining FREQUENCY dimension");
    const int masttilt_dim = scalar_dim("MASTTILT");

    const int xyz_var = define_var("XYZ", xyz_dim, "FEET");
    const int phi_var = define_var("PHI", phi_dim, "DEGREE");
    const int theta_var = define_var("THETA", theta_dim, "DEGREE");
    const int freq_var = define_var("FREQUENCY", freq_dim, "HERTZ");

    int amplitude_var = 0;
    const int amplitude_dims[3] = {phi_dim, theta_dim, freq_dim};
    nc_check(nc_def_var(file, "AMPLITUDE", NC_FLOAT, 3, amplitude_dims, &amplitude_var),
            "defining AMPLITUDE");
    nc_check(nc_put_att_text(file, amplitude_var, "unit", 7, "DECIBEL"), "annotating AMPLITUDE");
    std::string masttilt_unit("DEGREE");
    masttilt_unit.resize(20, ' ');
    const int masttilt_var = define_var("MASTTILT", masttilt_dim, masttilt_unit);

    const std::string title = "NICEOPS interpolated sphere";
    nc_check(nc_put_att_text(file, NC_GLOBAL, "title", title.size(), title.c_str()),
            "writing title");
    // The frame AMPLITUDE is in (flight_path_frame_amplitudes), which panam reads
    // back to rotate the sphere to level again.
    const std::string frame = "flight_path";
    nc_check(nc_put_att_text(file, NC_GLOBAL, "sphere_frame", frame.size(), frame.c_str()),
            "writing sphere_frame");
    nc_check(nc_enddef(file), "ending define mode");

    auto put_scalar = [&](int var, double value) {
        const float as_float = static_cast<float>(value);
        nc_check(nc_put_var_float(file, var, &as_float), "writing a scalar");
    };
    // AAM's sphere types (TR Table 4-37): BB + NB + PT = 1.
    const bool broadband = aam_third_octave_bands(sphere.frequencies);
    put_scalar(bb_var, broadband ? 1.0 : 0.0);
    put_scalar(nb_var, broadband ? 0.0 : 1.0);
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

    const auto flight_path = flight_path_frame_amplitudes(sphere);
    std::vector<float> amplitude(flight_path.begin(), flight_path.end());
    nc_check(nc_put_var_float(file, amplitude_var, amplitude.data()), "writing AMPLITUDE");

    nc_check(guard.close(), "closing " + path);
    log("Wrote AAM sphere to ", path, " (", sphere.phis.size(), " x ", sphere.thetas.size(),
            " x ", sphere.frequencies.size(), ").");
}

void Database::self_test(const std::string& datafile_path)
{
    std::size_t invalid = 0;
    std::size_t valid = 0;
    log("Begin database self test against file: ", datafile_path);
    const int file_ncid = open_database(datafile_path);
    NcGuard guard(file_ncid);

    // Only the format version is re-read: the constructor has already taken
    // everything else off this same file, and re-running its side effects here
    // would corrupt the object under test rather than check it.
    const int self_test_version = get_int_scalar(file_ncid, "database_version", 0);
    // Same flag the constructor reads -- see its comment on
    // shared_grid_and_frequency for when this holds.
    const bool shared_grid_and_frequency =
            get_int_scalar(file_ncid, "shared_grid_and_frequency", 0) != 0;
    std::vector<double> shared_phi_data, shared_theta_data;
    if (shared_grid_and_frequency)
    {
        const auto shared_nchannels = get_dim_length("channels", file_ncid);
        shared_theta_data.resize(shared_nchannels);
        require_array(file_ncid, "theta", shared_nchannels, shared_theta_data);
        shared_phi_data.resize(shared_nchannels);
        require_array(file_ncid, "phi", shared_nchannels, shared_phi_data);
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
        std::vector<double> theta_data, phi_data;
        if (shared_grid_and_frequency)
        {
            theta_data = shared_theta_data;
            phi_data = shared_phi_data;
        }
        else
        {
            auto nchannels = get_dim_length("channels", group_id);
            theta_data.resize(nchannels);
            require_array(group_id, "theta", nchannels, theta_data);

            phi_data.resize(nchannels);
            require_array(group_id, "phi", nchannels, phi_data);
            // Interpret the file exactly as the loader did, or this compares
            // the grid built from repaired labels against queries at the raw
            // ones.
            if (self_test_version < 1)
            {
                repair_legacy_hover_directivity(phi_data);
            }
        }
        const auto nchannels = phi_data.size();

        std::vector<double> dBA_data(nchannels);
        require_array(group_id, "dBA", nchannels, dBA_data);

        std::vector<double> EAA_data(nchannels);
        require_array(group_id, "EAA", nchannels, EAA_data);

        // Every channel of a sphere is queried at that sphere's own operating
        // condition, so the blend is the same for all of them -- resolve it
        // once instead of running a Delaunay locate per channel.
        const auto blend = blend_at(advance_ratio, flight_path_angle, thrust_coefficient);

        // Validate the database against the known data at sphere surface
        for (std::size_t ch = 0; ch < dBA_data.size(); ch++)
        {
            // eval() can only ever read back no_energy_floor_db for a
            // direction the constructor floored, never the original -inf --
            // see its member comment -- so the raw file value gets the same
            // floor applied here before the comparison means anything.
            auto dBA = std::isfinite(dBA_data[ch]) ? dBA_data[ch] : no_energy_floor_db;
            auto level = eval(dimensional_radius, phi_data[ch], theta_data[ch], blend);
            // At a sphere's own condition and a grid node the blend and the
            // bilinear weights are [1, 0, 0, 0], so eval() returns the stored
            // double itself (measured: 0 difference over 22 million samples on every
            // shipped database).  The 1% the D original allowed -- 0.8 dB at 80 dBA
            // -- let a channel filed one grid cell off pass 60-90% of the time, which
            // is the very error this test exists to catch.  The EAA is the other half
            // of the cell and was never looked at: it comes back through eval as the
            // loss over a further 1000 m, less the spreading.
            constexpr double exact = 1e-9;
            const double far = dimensional_radius + 1000.0;
            const double eaa = std::isfinite(EAA_data[ch]) ? EAA_data[ch] : 0.0;
            const double eaa_back = level - eval(far, phi_data[ch], theta_data[ch], blend)
                    - 20.0 * std::log10(far / dimensional_radius);
            if (!approx_equal(dBA, level, exact, exact) || !approx_equal(eaa, eaa_back, exact, exact))
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
        input_error("Self test failed on ", invalid, " of ", valid + invalid, " samples.");
    }
}

} // namespace niceops
