#include "urban.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <netcdf.h>

#include "screening.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"
#include "util/netcdf.hpp"
#include "util/paths.hpp"

namespace niceops {

namespace {

// The facade's base direction and its left normal, both unit, and its length.
struct FacadeFrame
{
    double tx, ty; // along the base, (x0, y0) -> (x1, y1)
    double nx, ny; // the left normal
    double length;
};

FacadeFrame frame_of(const Facade& f)
{
    const double dx = f.x1 - f.x0;
    const double dy = f.y1 - f.y0;
    const double length = std::hypot(dx, dy);
    return {dx / length, dy / length, -dy / length, dx / length, length};
}

// Where the segment a -> b crosses the facade's plane, as the fraction of the
// way from a, given the ends' signed distances (opposite signs).
Vec3 crossing(const Vec3& a, const Vec3& b, double da, double db)
{
    const double t = da / (da - db);
    return a + (b - a) * t;
}

// Whether a point already on the facade's plane lies within the rectangle.
bool within(const Facade& f, const FacadeFrame& frame, const Vec3& q)
{
    const double s = (q.x - f.x0) * frame.tx + (q.y - f.y0) * frame.ty;
    return s >= 0.0 && s <= frame.length && q.z >= f.bottom && q.z <= f.top;
}

bool blocked_by_any(const std::vector<Facade>& facades, std::size_t skip_a, std::size_t skip_b, const Vec3& a,
        const Vec3& b)
{
    for (std::size_t k = 0; k < facades.size(); k++)
        if (k != skip_a && k != skip_b && blocks(facades[k], a, b))
            return true;
    return false;
}

constexpr std::size_t none = static_cast<std::size_t>(-1);

void extend(const std::vector<Facade>& facades, const Vec3& source, const Vec3& receiver, std::size_t max_order,
        std::vector<std::size_t>& sequence, std::vector<Vec3>& images, std::vector<SpecularPath>& out)
{
    // Validate the current sequence: images[k] is the source mirrored in
    // facades sequence[0..k], so the last image is what the receiver sees.
    // Backtrack from the receiver toward each image in turn.
    const std::size_t order = sequence.size();
    Vec3 from = receiver;
    bool valid = true;
    for (std::size_t j = order; j-- > 0;)
    {
        const Facade& f = facades[sequence[j]];
        const Vec3& image = images[j];
        const double d_from = facade_side(f, from);
        const double d_image = facade_side(f, image);
        if (!((d_from > 0.0 && d_image < 0.0) || (d_from < 0.0 && d_image > 0.0)))
        {
            valid = false;
            break;
        }
        const Vec3 q = crossing(from, image, d_from, d_image);
        if (!within(f, frame_of(f), q))
        {
            valid = false;
            break;
        }
        // The leg from `from` to q, clear of every other facade; `from` lies on
        // the facade reflected in next (or is the receiver).
        const std::size_t next = j + 1 < order ? sequence[j + 1] : none;
        if (blocked_by_any(facades, sequence[j], next, from, q))
        {
            valid = false;
            break;
        }
        from = q;
    }
    if (valid)
    {
        const std::size_t first = order > 0 ? sequence[0] : none;
        if (!blocked_by_any(facades, first, none, from, source))
            out.push_back({sequence, ((order > 0 ? images.back() : source) - receiver).norm()});
    }
    if (order == max_order)
        return;
    for (std::size_t k = 0; k < facades.size(); k++)
    {
        if (order > 0 && sequence.back() == k)
            continue;
        sequence.push_back(k);
        images.push_back(mirror_in(facades[k], order > 0 ? images.back() : source));
        extend(facades, source, receiver, max_order, sequence, images, out);
        images.pop_back();
        sequence.pop_back();
    }
}

void check_axes(const UrbanTable::Axes& axes)
{
    auto ascending = [](const std::vector<double>& v) {
        for (std::size_t k = 1; k < v.size(); k++)
            if (!(v[k] > v[k - 1]))
                return false;
        return !v.empty() && std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
    };
    if (!ascending(axes.elevations) || !ascending(axes.bands))
        throw std::invalid_argument("urban table: the elevation and band axes must be finite and ascend");
    if (axes.elevations.front() < -90.0 || axes.elevations.back() > 90.0)
        throw std::invalid_argument("urban table: elevations must lie within [-90, 90] degrees");
    if (!(axes.bands.front() > 0.0))
        throw std::invalid_argument("urban table: band centers must be positive");
    const std::size_t n = axes.azimuths.size();
    if (n == 0)
        throw std::invalid_argument("urban table: no azimuths");
    for (std::size_t k = 0; k < n; k++)
        if (std::abs(axes.azimuths[k] - 360.0 * static_cast<double>(k) / static_cast<double>(n)) > 1e-9)
            throw std::invalid_argument("urban table: azimuths must be uniform from 0 over 360 degrees");
}

void nc_check(int status, const std::string& what, const std::string& path)
{
    if (status != NC_NOERR)
        input_error("Urban table ", path, ": ", what, " (", nc_strerror(status), ").");
}

constexpr const char* facade_fields[6] = {"facade_x0", "facade_y0", "facade_x1", "facade_y1", "facade_bottom",
        "facade_top"};

} // namespace

double facade_side(const Facade& facade, const Vec3& p)
{
    const FacadeFrame f = frame_of(facade);
    return (p.x - facade.x0) * f.nx + (p.y - facade.y0) * f.ny;
}

Vec3 mirror_in(const Facade& facade, const Vec3& p)
{
    const FacadeFrame f = frame_of(facade);
    const double d = (p.x - facade.x0) * f.nx + (p.y - facade.y0) * f.ny;
    return Vec3(p.x - 2.0 * d * f.nx, p.y - 2.0 * d * f.ny, p.z);
}

Vec3 mirror_direction(const Facade& facade, const Vec3& v)
{
    const FacadeFrame f = frame_of(facade);
    const double d = v.x * f.nx + v.y * f.ny;
    return Vec3(v.x - 2.0 * d * f.nx, v.y - 2.0 * d * f.ny, v.z);
}

bool blocks(const Facade& facade, const Vec3& a, const Vec3& b)
{
    const double da = facade_side(facade, a);
    const double db = facade_side(facade, b);
    if (!((da > 0.0 && db < 0.0) || (da < 0.0 && db > 0.0)))
        return false;
    return within(facade, frame_of(facade), crossing(a, b, da, db));
}

bool first_order_image(const std::vector<Facade>& facades, std::size_t self, const Vec3& source,
        const Vec3& microphone, Vec3& image)
{
    const Facade& f = facades[self];
    const double ds = facade_side(f, source);
    const double dm = facade_side(f, microphone);
    if (!((ds > 0.0 && dm > 0.0) || (ds < 0.0 && dm < 0.0)))
        return false;
    const Vec3 mirrored = mirror_in(f, microphone);
    // The reflection point: where source -> image crosses the plane, at the
    // fraction ds / (ds + dm) of the way, the image being -dm off it.
    const Vec3 q = crossing(source, mirrored, ds, -dm);
    if (!within(f, frame_of(f), q))
        return false;
    if (blocked_by_any(facades, self, none, microphone, q) || blocked_by_any(facades, self, none, q, source))
        return false;
    image = mirrored;
    return true;
}

bool first_order_floor_partner(const std::vector<Facade>& below, std::size_t self, const Vec3& source,
        const Vec3& microphone, double floor)
{
    Vec3 image;
    return first_order_image(below, self, source, Vec3(microphone.x, microphone.y, 2.0 * floor - microphone.z),
            image);
}

std::vector<Facade> walls_below(std::vector<Facade> facades)
{
    for (auto& f : facades)
        f.bottom = -std::numeric_limits<double>::infinity();
    return facades;
}

std::vector<SpecularPath> specular_paths(const std::vector<Facade>& facades, const Vec3& source,
        const Vec3& receiver, std::size_t max_order)
{
    std::vector<SpecularPath> out;
    std::vector<std::size_t> sequence;
    std::vector<Vec3> images;
    extend(facades, source, receiver, max_order, sequence, images, out);
    return out;
}

CanyonTrace trace_canyon(const IdealCanyon& canyon, double v_c, double u_z, int max_bounces)
{
    CanyonTrace trace;
    if (v_c == 0.0)
    {
        // Along the axis: it never reaches a wall.
        trace.escaped = true;
        return trace;
    }
    double y = canyon.offset;
    double z = canyon.microphone_height;
    for (;;)
    {
        const bool toward_b = v_c > 0.0;
        const double across = toward_b ? canyon.width - y : y;
        z += u_z * across / std::abs(v_c);
        if (trace.first_wall == 0)
            trace.first_wall = toward_b ? 'b' : 'a';
        if (z > (toward_b ? canyon.height_b : canyon.height_a))
        {
            trace.escaped = true;
            return trace;
        }
        if (trace.bounces == max_bounces)
            return trace;
        trace.bounces++;
        (toward_b ? trace.hits_b : trace.hits_a)++;
        y = toward_b ? canyon.width : 0.0;
        v_c = -v_c;
    }
}

double edge_path_difference(double dy, double dz, double u_c, double u_z)
{
    return std::hypot(u_c, u_z) * std::hypot(dy, dz) - (u_c * dy + u_z * dz);
}

void canyon_ratio(const IdealCanyon& canyon, double u_c, double u_z, const double* reflection_a,
        const double* reflection_b, double floor_reflection, const std::vector<double>& bands, double sound_speed,
        bool diffraction, double* out)
{
    const std::size_t n = bands.size();
    // Along the axis the walls are never reached, and above both of them they
    // are not in the way: open ground.
    if (u_c == 0.0 || canyon.microphone_height > std::max(canyon.height_a, canyon.height_b))
    {
        std::fill(out, out + n, 1.0);
        return;
    }
    // Four families of rays reach the microphone: arriving from above or off
    // the floor, each with the plane wave's horizontal direction or its mirror.
    // Traced back, a ray off the floor leaves the microphone's image in the
    // floor; the direct families reach the source after an even number of
    // wall reflections and the mirrored ones after an odd number, and each ray
    // reflects as often as the walls make it, so each family holds one path at
    // most.  Open ground has two of them, the direct ray and the floor's.
    IdealCanyon below = canyon;
    below.microphone_height = -canyon.microphone_height;
    const CanyonTrace direct = trace_canyon(canyon, u_c, u_z);
    const CanyonTrace mirror = trace_canyon(canyon, -u_c, u_z);
    const CanyonTrace direct_floor = trace_canyon(below, u_c, u_z);
    const CanyonTrace mirror_floor = trace_canyon(below, -u_c, u_z);
    const bool open = direct.escaped && direct.bounces == 0;
    auto direct_valid = [](const CanyonTrace& t) { return t.escaped && t.bounces % 2 == 0; };
    auto mirror_valid = [](const CanyonTrace& t) { return t.escaped && t.bounces % 2 == 1; };
    // The first-order reflection is the run's to evaluate whenever it is open
    // from above, and its floor partner, through the ground model, whenever
    // that is open too (first_order_floor_partner); whatever the run does not
    // evaluate is the table's.
    const bool run_reflects = mirror_valid(mirror) && mirror.bounces == 1;
    const bool run_reflects_floor = run_reflects && mirror_valid(mirror_floor) && mirror_floor.bounces == 1;
    // The blocking edge, for the diffracted path: the top of the first wall the
    // direct ray reaches, in the cross-section from the microphone.
    double delta = 0.0;
    if (!open && diffraction)
    {
        const bool b = direct.first_wall == 'b';
        const double dy = b ? canyon.width - canyon.offset : -canyon.offset;
        const double dz = (b ? canyon.height_b : canyon.height_a) - canyon.microphone_height;
        delta = edge_path_difference(dy, dz, u_c, u_z);
    }
    const double open_ground = 1.0 + floor_reflection;
    for (std::size_t band = 0; band < n; band++)
    {
        auto energy = [&](const CanyonTrace& t) {
            return std::pow(reflection_a[band], t.hits_a) * std::pow(reflection_b[band], t.hits_b);
        };
        double sum = 0.0;
        if (direct_valid(direct))
            sum += energy(direct);
        if (direct_valid(direct_floor))
            sum += floor_reflection * energy(direct_floor);
        if (!run_reflects && mirror_valid(mirror))
            sum += energy(mirror);
        if (!run_reflects_floor && mirror_valid(mirror_floor))
            sum += floor_reflection * energy(mirror_floor);
        // Where the walls block the direct ray, the diffracted path over the
        // blocking edge, with its own floor reflection.
        if (!open && diffraction)
        {
            const double loss = std::min(cnossos_dif(delta, 0.0, 1, sound_speed / bands[band]), cnossos_cap);
            sum += open_ground * std::pow(10.0, -0.1 * loss);
        }
        out[band] = sum / open_ground;
    }
}

CompassDirection compass_direction(const Vec3& local, double frame_bearing)
{
    const double b = frame_bearing * pi / 180.0;
    const double east = local.x * std::sin(b) - local.y * std::cos(b);
    const double north = local.x * std::cos(b) + local.y * std::sin(b);
    CompassDirection d;
    d.azimuth = std::atan2(east, north) * 180.0 / pi;
    if (d.azimuth < 0.0)
        d.azimuth += 360.0;
    d.elevation = std::atan2(local.z, std::hypot(local.x, local.y)) * 180.0 / pi;
    return d;
}

UrbanTable::UrbanTable(Axes axes, std::vector<Receiver> receivers, std::vector<double> gains, double frame_bearing,
        std::string method, std::map<std::string, std::string> provenance)
    : axes_(std::move(axes))
    , receivers_(std::move(receivers))
    , gains_(std::move(gains))
    , frame_bearing_(frame_bearing)
    , method_(std::move(method))
    , provenance_(std::move(provenance))
{
    check_axes(axes_);
    if (!std::isfinite(frame_bearing_))
        throw std::invalid_argument("urban table: the frame bearing must be finite");
    const std::size_t nb = axes_.bands.size();
    for (const Receiver& r : receivers_)
    {
        if (!(std::isfinite(r.microphone.x) && std::isfinite(r.microphone.y) && std::isfinite(r.microphone.z)))
            throw std::invalid_argument("urban table: a receiver's position is not finite");
        if (r.reflection.size() != r.facades.size() * nb)
            throw std::invalid_argument("urban table: a receiver's reflection coefficients do not fill its "
                                        "facades and bands");
        for (double c : r.reflection)
            if (!(c >= 0.0 && c <= 1.0))
                throw std::invalid_argument("urban table: reflection coefficients must lie within [0, 1]");
        for (const Facade& f : r.facades)
        {
            const double length = std::hypot(f.x1 - f.x0, f.y1 - f.y0);
            if (!(std::isfinite(length) && length > 0.0 && std::isfinite(f.bottom) && std::isfinite(f.top)
                        && f.top >= f.bottom))
                throw std::invalid_argument("urban table: a facade needs a base of positive length and a top at "
                                            "or above its bottom");
        }
    }
    if (gains_.size() != receivers_.size() * axes_.azimuths.size() * axes_.elevations.size() * nb)
        throw std::invalid_argument("urban table: the gains do not fill the receivers and axes");
    for (double g : gains_)
        if (std::isinf(g))
            throw std::invalid_argument("urban table: a gain is infinite (NaN marks no value)");
}

void UrbanTable::save(const std::string& path) const
{
    ensure_parent_directory(path);
    int ncid = 0;
    nc_check(nc_create(path.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "cannot create it", path);
    std::size_t facade_count = 0;
    for (const Receiver& r : receivers_)
        facade_count += r.facades.size();
    const std::size_t nb = axes_.bands.size();

    int receiver_dim = 0, azimuth_dim = 0, elevation_dim = 0, band_dim = 0, facade_dim = 0;
    if (receivers_.empty())
        throw std::invalid_argument("urban table: no receivers to save");
    nc_check(nc_def_dim(ncid, "receiver", receivers_.size(), &receiver_dim), "defining a dimension", path);
    nc_check(nc_def_dim(ncid, "azimuth", axes_.azimuths.size(), &azimuth_dim), "defining a dimension", path);
    nc_check(nc_def_dim(ncid, "elevation", axes_.elevations.size(), &elevation_dim), "defining a dimension", path);
    nc_check(nc_def_dim(ncid, "band", nb, &band_dim), "defining a dimension", path);
    // A zero length would make the dimension unlimited; a table with no
    // facades has no facade variables at all.
    if (facade_count > 0)
        nc_check(nc_def_dim(ncid, "facade", facade_count, &facade_dim), "defining a dimension", path);

    auto define = [&](const char* name, int ndims, const int* dims, const char* units, nc_type type = NC_DOUBLE) {
        int var = 0;
        nc_check(nc_def_var(ncid, name, type, ndims, dims, &var), std::string("defining ") + name, path);
        if (units != nullptr)
            nc_put_att_text(ncid, var, "units", std::string(units).size(), units);
        return var;
    };
    const int azimuth_var = define("azimuth", 1, &azimuth_dim, "degrees");
    const int elevation_var = define("elevation", 1, &elevation_dim, "degrees");
    const int band_var = define("band", 1, &band_dim, "Hz");
    const int rx = define("receiver_x", 1, &receiver_dim, "m");
    const int ry = define("receiver_y", 1, &receiver_dim, "m");
    const int rz = define("receiver_z", 1, &receiver_dim, "m");
    const int rf = define("receiver_facades", 1, &receiver_dim, nullptr, NC_INT);
    int facade_vars[6] = {};
    int reflection_var = 0;
    if (facade_count > 0)
    {
        for (int k = 0; k < 6; k++)
            facade_vars[k] = define(facade_fields[k], 1, &facade_dim, "m");
        const int dims[2] = {facade_dim, band_dim};
        reflection_var = define("facade_reflection", 2, dims, "energy ratio");
    }
    const int gain_dims[4] = {receiver_dim, azimuth_dim, elevation_dim, band_dim};
    const int gain_var = define("gain", 4, gain_dims,
            "dB re the open-ground direct path; NaN where there is no value");

    const int version = format_version;
    nc_put_att_int(ncid, NC_GLOBAL, "urban_table_version", NC_INT, 1, &version);
    nc_put_att_double(ncid, NC_GLOBAL, "frame_bearing", NC_DOUBLE, 1, &frame_bearing_);
    nc_put_att_text(ncid, NC_GLOBAL, "method", method_.size(), method_.c_str());
    const std::string azimuth_note = "compass bearing from the microphone toward the source";
    nc_put_att_text(ncid, NC_GLOBAL, "azimuth_convention", azimuth_note.size(), azimuth_note.c_str());
    for (const auto& [key, value] : provenance_)
        nc_put_att_text(ncid, NC_GLOBAL, ("provenance_" + key).c_str(), value.size(), value.c_str());
    nc_check(nc_enddef(ncid), "ending its definition", path);

    nc_check(nc_put_var_double(ncid, azimuth_var, axes_.azimuths.data()), "writing azimuth", path);
    nc_check(nc_put_var_double(ncid, elevation_var, axes_.elevations.data()), "writing elevation", path);
    nc_check(nc_put_var_double(ncid, band_var, axes_.bands.data()), "writing band", path);
    std::vector<double> xs, ys, zs;
    std::vector<int> counts;
    std::vector<double> fields[6];
    std::vector<double> reflection;
    for (const Receiver& r : receivers_)
    {
        xs.push_back(r.microphone.x);
        ys.push_back(r.microphone.y);
        zs.push_back(r.microphone.z);
        counts.push_back(static_cast<int>(r.facades.size()));
        for (const Facade& f : r.facades)
        {
            const double values[6] = {f.x0, f.y0, f.x1, f.y1, f.bottom, f.top};
            for (int k = 0; k < 6; k++)
                fields[k].push_back(values[k]);
        }
        reflection.insert(reflection.end(), r.reflection.begin(), r.reflection.end());
    }
    nc_check(nc_put_var_double(ncid, rx, xs.data()), "writing receiver_x", path);
    nc_check(nc_put_var_double(ncid, ry, ys.data()), "writing receiver_y", path);
    nc_check(nc_put_var_double(ncid, rz, zs.data()), "writing receiver_z", path);
    nc_check(nc_put_var_int(ncid, rf, counts.data()), "writing receiver_facades", path);
    if (facade_count > 0)
    {
        for (int k = 0; k < 6; k++)
            nc_check(nc_put_var_double(ncid, facade_vars[k], fields[k].data()), "writing the facades", path);
        nc_check(nc_put_var_double(ncid, reflection_var, reflection.data()), "writing facade_reflection", path);
    }
    nc_check(nc_put_var_double(ncid, gain_var, gains_.data()), "writing gain", path);
    nc_check(nc_close(ncid), "closing it", path);
}

UrbanTable UrbanTable::load(const std::string& path)
{
    int ncid = 0;
    if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
        input_error("Cannot open urban table ", path, ".");
    struct Closer
    {
        int id;
        ~Closer()
        {
            nc_close(id);
        }
    } closer {ncid};

    int version = 0;
    switch (netcdf::single_value(ncid, NC_GLOBAL, "urban_table_version", version))
    {
    case netcdf::Single::missing:
        input_error("Urban table ", path, " has no urban_table_version attribute: it is not one, or it predates "
                                          "the format.");
    case netcdf::Single::not_one_number:
        input_error("Urban table ", path, ": urban_table_version must be one integer.");
    case netcdf::Single::ok:
        break;
    }
    if (version != format_version)
        input_error("Urban table ", path, " is format version ", version, "; this build reads ", format_version,
                ".");
    double frame_bearing = 0.0;
    if (netcdf::single_value(ncid, NC_GLOBAL, "frame_bearing", frame_bearing) != netcdf::Single::ok)
        input_error("Urban table ", path, " must record frame_bearing, the compass bearing of its local +x axis, "
                                          "as one number.");

    auto variable = [&](const char* name, const std::vector<std::string>& dims) {
        int var = 0;
        if (nc_inq_varid(ncid, name, &var) != NC_NOERR)
            input_error("Urban table ", path, " has no '", name, "' variable.");
        if (netcdf::dimension_names(ncid, var) != dims)
        {
            std::string expected;
            for (const auto& d : dims)
                expected += (expected.empty() ? "" : ", ") + d;
            input_error("Urban table ", path, ": '", name, "' must be dimensioned (", expected, ").");
        }
        auto values = netcdf::read_doubles(ncid, var);
        if (!values)
            input_error("Urban table ", path, ": cannot read '", name, "'.");
        return std::move(*values);
    };
    Axes axes;
    axes.azimuths = variable("azimuth", {"azimuth"});
    axes.elevations = variable("elevation", {"elevation"});
    axes.bands = variable("band", {"band"});
    const auto xs = variable("receiver_x", {"receiver"});
    const auto ys = variable("receiver_y", {"receiver"});
    const auto zs = variable("receiver_z", {"receiver"});
    const auto counts = variable("receiver_facades", {"receiver"});
    std::size_t facade_total = 0;
    for (double c : counts)
    {
        if (!(c >= 0.0) || c != std::floor(c))
            input_error("Urban table ", path, ": receiver_facades must hold counts.");
        facade_total += static_cast<std::size_t>(c);
    }
    std::vector<double> fields[6];
    std::vector<double> reflection;
    int facade_dim = 0;
    const bool has_facades = nc_inq_dimid(ncid, "facade", &facade_dim) == NC_NOERR;
    if (has_facades)
    {
        for (int k = 0; k < 6; k++)
            fields[k] = variable(facade_fields[k], {"facade"});
        reflection = variable("facade_reflection", {"facade", "band"});
    }
    if (fields[0].size() != facade_total)
        input_error("Urban table ", path, ": receiver_facades counts ", facade_total, " facades, and the file "
                    "holds ", fields[0].size(), ".");
    auto gains = variable("gain", {"receiver", "azimuth", "elevation", "band"});

    std::vector<Receiver> receivers(xs.size());
    const std::size_t nb = axes.bands.size();
    std::size_t next = 0;
    for (std::size_t r = 0; r < receivers.size(); r++)
    {
        receivers[r].microphone = Vec3(xs[r], ys[r], zs[r]);
        const auto count = static_cast<std::size_t>(counts[r]);
        for (std::size_t k = 0; k < count; k++, next++)
        {
            receivers[r].facades.push_back(
                    {fields[0][next], fields[1][next], fields[2][next], fields[3][next], fields[4][next], fields[5][next]});
            receivers[r].reflection.insert(receivers[r].reflection.end(), reflection.begin() + next * nb,
                    reflection.begin() + (next + 1) * nb);
        }
    }
    auto text = [&](const std::string& name) {
        std::size_t length = 0;
        if (nc_inq_attlen(ncid, NC_GLOBAL, name.c_str(), &length) != NC_NOERR)
            return std::string();
        std::string value(length, '\0');
        nc_get_att_text(ncid, NC_GLOBAL, name.c_str(), value.data());
        return value;
    };
    std::map<std::string, std::string> provenance;
    int attributes = 0;
    nc_inq_natts(ncid, &attributes);
    for (int k = 0; k < attributes; k++)
    {
        char name[NC_MAX_NAME + 1] = {};
        nc_inq_attname(ncid, NC_GLOBAL, k, name);
        const std::string key(name);
        if (key.rfind("provenance_", 0) == 0)
            provenance[key.substr(11)] = text(key);
    }
    try
    {
        return UrbanTable(std::move(axes), std::move(receivers), std::move(gains), frame_bearing, text("method"),
                std::move(provenance));
    }
    catch (const std::invalid_argument& e)
    {
        input_error("Urban table ", path, ": ", e.what(), ".");
    }
}

UrbanTable::BandMap UrbanTable::match_bands(const std::vector<double>& run_bands, std::size_t first_band,
        std::size_t last_band) const
{
    last_band = std::min(last_band, run_bands.size());
    BandMap map;
    map.first = first_band;
    map.last = last_band;
    map.column.assign(run_bands.size(), 0);
    std::string missing;
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const auto& table = axes_.bands;
        const auto nearest = std::min_element(table.begin(), table.end(), [&](double a, double b) {
            return std::abs(std::log(a / run_bands[band])) < std::abs(std::log(b / run_bands[band]));
        });
        if (std::abs(*nearest / run_bands[band] - 1.0) > 0.005)
            missing += (missing.empty() ? "" : ", ") + std::to_string(run_bands[band]);
        map.column[band] = static_cast<std::size_t>(nearest - table.begin());
    }
    if (!missing.empty())
        throw std::runtime_error("the urban table has no band at " + missing + " Hz, which the run's band limits "
                                                                               "include");
    return map;
}

std::size_t UrbanTable::find_receiver(const Vec3& microphone, double tolerance) const
{
    std::size_t best = receivers_.size();
    double nearest = tolerance;
    for (std::size_t r = 0; r < receivers_.size(); r++)
    {
        const double d = (receivers_[r].microphone - microphone).norm();
        if (d <= nearest)
        {
            nearest = d;
            best = r;
        }
    }
    return best;
}

bool UrbanTable::lookup(std::size_t receiver, double azimuth_degrees, double elevation_degrees,
        const BandMap& bands, double* out) const
{
    // Elevation: a one-point axis holds only its point.
    const auto& el = axes_.elevations;
    std::size_t e0 = 0;
    double fe = 0.0;
    if (el.size() == 1)
    {
        if (std::abs(elevation_degrees - el[0]) > 1e-6)
            return false;
    }
    else
    {
        if (!(elevation_degrees >= el.front()) || !(elevation_degrees <= el.back()))
            return false;
        const auto upper = std::upper_bound(el.begin(), el.end(), elevation_degrees);
        e0 = std::min(static_cast<std::size_t>(upper - el.begin()), el.size() - 1) - 1;
        fe = (elevation_degrees - el[e0]) / (el[e0 + 1] - el[e0]);
    }
    // Azimuth: periodic, uniform.
    const std::size_t na = axes_.azimuths.size();
    std::size_t a0 = 0, a1 = 0;
    double fa = 0.0;
    if (na > 1)
    {
        double u = azimuth_degrees >= 0.0 && azimuth_degrees < 360.0 ? azimuth_degrees
                                                                      : std::fmod(azimuth_degrees, 360.0);
        if (u < 0.0)
            u += 360.0;
        u *= static_cast<double>(na) / 360.0;
        a0 = std::min(static_cast<std::size_t>(u), na - 1);
        fa = u - static_cast<double>(a0);
        a1 = (a0 + 1) % na;
    }
    // The four corners; one of zero weight is not read, so a NaN just past a
    // direction on the grid does not spoil it.
    struct Corner
    {
        std::size_t offset;
        double weight;
    };
    Corner corners[4];
    int used = 0;
    const std::size_t azimuth_index[2] = {a0, a1};
    const double azimuth_weight[2] = {1.0 - fa, fa};
    const double elevation_weight[2] = {1.0 - fe, fe};
    const std::size_t nb = axes_.bands.size();
    for (int a = 0; a < 2; a++)
        for (int e = 0; e < 2; e++)
        {
            const double weight = azimuth_weight[a] * elevation_weight[e];
            if (weight == 0.0)
                continue;
            corners[used++] = {cell(receiver, azimuth_index[a], e0 + e) * nb, weight};
        }
    for (std::size_t band = bands.first; band < bands.last; band++)
    {
        const std::size_t column = bands.column[band];
        double sum = 0.0;
        for (int c = 0; c < used; c++)
            sum += corners[c].weight * gains_[corners[c].offset + column];
        if (!std::isfinite(sum))
            return false;
        out[band] = sum;
    }
    return true;
}

std::vector<double> UrbanTable::max_power_ratio(const BandMap& bands) const
{
    std::vector<double> ratio(bands.column.size(), 1.0);
    const std::size_t nb = axes_.bands.size();
    const std::size_t cells = gains_.size() / nb;
    for (std::size_t band = bands.first; band < bands.last; band++)
    {
        double largest = -std::numeric_limits<double>::infinity();
        for (std::size_t c = 0; c < cells; c++)
        {
            const double v = gains_[c * nb + bands.column[band]];
            if (std::isfinite(v))
                largest = std::max(largest, v);
        }
        ratio[band] = std::max(1.0, std::pow(10.0, 0.1 * largest));
    }
    return ratio;
}

std::vector<Facade> ideal_canyon_facades(const IdealCanyon& canyon, const Vec3& microphone, double floor,
        double axis_bearing, double frame_bearing, double length)
{
    // East/north of the axis and of its right, toward wall B; then local.
    const double a = axis_bearing * pi / 180.0;
    const double b = frame_bearing * pi / 180.0;
    auto local = [&](double east, double north) {
        return std::array<double, 2> {east * std::sin(b) + north * std::cos(b), -east * std::cos(b) + north * std::sin(b)};
    };
    const auto axis = local(std::sin(a), std::cos(a));
    const auto right = local(std::cos(a), -std::sin(a));
    auto wall = [&](double across, double height) {
        const double cx = microphone.x + across * right[0];
        const double cy = microphone.y + across * right[1];
        return Facade {cx - 0.5 * length * axis[0], cy - 0.5 * length * axis[1], cx + 0.5 * length * axis[0],
                cy + 0.5 * length * axis[1], floor, floor + height};
    };
    return {wall(-canyon.offset, canyon.height_a), wall(canyon.width - canyon.offset, canyon.height_b)};
}

std::vector<double> ideal_canyon_gains(const IdealCanyon& canyon, double axis_bearing,
        const UrbanTable::Axes& axes, const std::vector<double>& reflection_a,
        const std::vector<double>& reflection_b, double floor_reflection, double sound_speed, bool diffraction)
{
    const std::size_t nb = axes.bands.size();
    std::vector<double> gains(axes.azimuths.size() * axes.elevations.size() * nb);
    std::vector<double> ratio(nb);
    const double a = axis_bearing * pi / 180.0;
    for (std::size_t i = 0; i < axes.azimuths.size(); i++)
    {
        const double azimuth = axes.azimuths[i] * pi / 180.0;
        for (std::size_t j = 0; j < axes.elevations.size(); j++)
        {
            const double e = axes.elevations[j] * pi / 180.0;
            // The arrival direction's component across the canyon, toward B:
            // cos(e) sin(azimuth - axis).
            const double u_c = std::cos(e) * std::sin(azimuth - a);
            canyon_ratio(canyon, u_c, std::sin(e), reflection_a.data(), reflection_b.data(), floor_reflection, axes.bands,
                    sound_speed, diffraction, ratio.data());
            for (std::size_t band = 0; band < nb; band++)
                gains[(i * axes.elevations.size() + j) * nb + band]
                        = std::max(10.0 * std::log10(ratio[band]), urban_gain_floor_db);
        }
    }
    return gains;
}

std::vector<double> wave_band_weights(const std::vector<double>& bands, double wave_below, std::size_t blend)
{
    std::vector<double> weight(bands.size(), 0.0);
    if (!(wave_below > 0.0))
        return weight;
    std::size_t last = bands.size();
    for (std::size_t b = 0; b < bands.size(); b++)
        if (bands[b] <= wave_below * (1.0 + 1e-9))
        {
            weight[b] = 1.0;
            last = b;
        }
    if (last == bands.size())
        return weight;
    for (std::size_t j = 1; j <= blend && last + j < bands.size(); j++)
        weight[last + j] = 1.0 - static_cast<double>(j) / static_cast<double>(blend + 1);
    return weight;
}

void blend_wave_band(std::vector<double>& gains, std::vector<double>& reflection, const std::vector<double>& wave,
        const std::vector<double>& weight)
{
    const std::size_t nb = weight.size();
    if (gains.size() != wave.size() || nb == 0 || gains.size() % nb != 0 || reflection.size() % nb != 0)
        throw std::invalid_argument("blend_wave_band: the tables and weights do not match");
    for (std::size_t cell = 0; cell < gains.size() / nb; cell++)
        for (std::size_t b = 0; b < nb; b++)
        {
            const double w = weight[b];
            if (w == 0.0)
                continue;
            double& g = gains[cell * nb + b];
            const double ratio = w * std::pow(10.0, 0.1 * wave[cell * nb + b]) + (1.0 - w) * std::pow(10.0, 0.1 * g);
            g = std::max(10.0 * std::log10(ratio), urban_gain_floor_db);
        }
    for (std::size_t f = 0; f < reflection.size() / nb; f++)
        for (std::size_t b = 0; b < nb; b++)
            reflection[f * nb + b] *= 1.0 - weight[b];
}

} // namespace niceops
