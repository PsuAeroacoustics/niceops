#include "sphere_energy.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "util/logger.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

void nc_check(int status, const std::string& what, const std::string& path)
{
    if (status != NC_NOERR)
        input_error("--provenance_spheres ", path, ": netCDF error while ", what, ": ", nc_strerror(status), ".");
}

// Closes on scope exit unless closed explicitly, which is the only close whose
// status is looked at.
struct NcFile
{
    int id = -1;
    ~NcFile()
    {
        if (id >= 0)
            nc_close(id);
    }
    int close()
    {
        const int status = id >= 0 ? nc_close(id) : NC_NOERR;
        id = -1;
        return status;
    }
};

void put_text(int ncid, int var, const char* name, const std::string& value, const std::string& path)
{
    nc_check(nc_put_att_text(ncid, var, name, value.size(), value.c_str()),
            std::string("writing attribute ") + name, path);
}

} // namespace

namespace niceops {

void write_sphere_energy(const std::string& path, const SphereEnergyFile& meta,
        const std::vector<std::vector<std::pair<std::uint32_t, double>>>& energy)
{
    const std::size_t observers = energy.size();
    const std::size_t spheres = meta.conditions.size();
    assert(meta.x.size() == observers && meta.y.size() == observers && meta.exposure.size() == observers);
    ensure_parent_directory(path);
    NcFile file;
    nc_check(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file.id), "creating it", path);
    const int ncid = file.id;

    int observer_dim = -1, sphere_dim = -1;
    nc_check(nc_def_dim(ncid, "observer", observers, &observer_dim), "defining observer", path);
    nc_check(nc_def_dim(ncid, "sphere", spheres, &sphere_dim), "defining sphere", path);
    auto define = [&](const char* name, nc_type type, std::vector<int> dims, const char* units,
                          const std::string& long_name) {
        int var = -1;
        nc_check(nc_def_var(ncid, name, type, static_cast<int>(dims.size()), dims.data(), &var),
                std::string("defining ") + name, path);
        put_text(ncid, var, "units", units, path);
        put_text(ncid, var, "long_name", long_name, path);
        return var;
    };
    const int x_var = define("x", NC_DOUBLE, {observer_dim}, "m", "observer x in the run's local frame");
    const int y_var = define("y", NC_DOUBLE, {observer_dim}, "m", "observer y in the run's local frame");
    const int exposure_var = define("exposure_db", NC_DOUBLE, {observer_dim}, "dB",
            "the exposure shared out: the SEL for one operation, the energy sum of the SELs for several");
    const int mu_var = define("advance_ratio", NC_DOUBLE, {sphere_dim}, "1", "the sphere's stored advance ratio");
    const int gamma_var = define("flight_path_angle", NC_DOUBLE, {sphere_dim}, "degree",
            "the sphere's stored flight path angle");
    const int ct_var = define("thrust_coefficient", NC_DOUBLE, {sphere_dim}, "1",
            "the sphere's stored thrust coefficient");
    const int share_var = define("energy_share", NC_FLOAT, {observer_dim, sphere_dim}, "1",
            "share of the observer's energy each sphere supplied, by the states' blend weights");
    const float fill = std::numeric_limits<float>::quiet_NaN();
    nc_check(nc_def_var_fill(ncid, share_var, 0, &fill), "setting the fill value", path);
    if (observers > 0 && spheres > 0)
    {
        // A row per observer, mostly zeros: compression shrinks it to the
        // dozen spheres an observer hears.
        const std::size_t chunk[2] = {std::min<std::size_t>(observers, 16), spheres};
        nc_check(nc_def_var_chunking(ncid, share_var, NC_CHUNKED, chunk), "setting the chunking", path);
        nc_check(nc_def_var_deflate(ncid, share_var, 1, 1, 1), "setting compression", path);
    }
    put_text(ncid, NC_GLOBAL, "database", meta.database, path);
    put_text(ncid, NC_GLOBAL, "sphere_order", "sphere k is the database's k-th group, as niceops loads them", path);
    put_text(ncid, NC_GLOBAL, "observer_order",
            meta.grid ? "the grid's, x fastest, as -o writes it" : "the observer file's (-b)", path);
    nc_check(nc_enddef(ncid), "ending its definition", path);

    nc_check(nc_put_var_double(ncid, x_var, meta.x.data()), "writing x", path);
    nc_check(nc_put_var_double(ncid, y_var, meta.y.data()), "writing y", path);
    std::vector<double> exposure_db(observers);
    for (std::size_t i = 0; i < observers; i++)
        exposure_db[i] = 10.0 * std::log10(meta.exposure[i]);
    nc_check(nc_put_var_double(ncid, exposure_var, exposure_db.data()), "writing exposure_db", path);
    std::vector<double> column(spheres);
    for (int axis = 0; axis < 3; axis++)
    {
        for (std::size_t k = 0; k < spheres; k++)
            column[k] = meta.conditions[k][axis];
        nc_check(nc_put_var_double(ncid, axis == 0 ? mu_var : axis == 1 ? gamma_var : ct_var, column.data()),
                "writing the conditions", path);
    }
    std::vector<float> row(spheres);
    for (std::size_t i = 0; i < observers; i++)
    {
        const double total = meta.exposure[i];
        std::fill(row.begin(), row.end(), total > 0.0 ? 0.0f : fill);
        if (total > 0.0)
        {
            for (const auto& [sphere, value] : energy[i])
            {
                assert(sphere < spheres);
                row[sphere] = static_cast<float>(value / total);
            }
        }
        if (spheres == 0)
            continue;
        const std::size_t start[2] = {i, 0};
        const std::size_t count[2] = {1, spheres};
        nc_check(nc_put_vara_float(ncid, share_var, start, count, row.data()), "writing energy_share", path);
    }
    nc_check(file.close(), "closing it", path);
}

} // namespace niceops
