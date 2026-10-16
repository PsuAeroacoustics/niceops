#include "band_frames.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "util/logger.hpp"
#include "util/netcdf_write.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr const char* option = "--band_frames";

void nc_check(int status, const std::string& what, const std::string& path)
{
    netcdf::check_write(status, option, what, path);
}

void put_text(int ncid, int var, const char* name, const std::string& value, const std::string& path)
{
    netcdf::put_text(ncid, var, name, value, option, path);
}

} // namespace

namespace niceops {

void write_band_frames(const std::string& path, const BandFramesFile& meta,
        const std::vector<ObserverFrames>& frames)
{
    const std::size_t observers = frames.size();
    const std::size_t bands = meta.band_centers.size();
    if (meta.x.size() != observers || meta.y.size() != observers || meta.z.size() != observers
            || meta.height.size() != observers)
        input_error("--band_frames ", path, ": internal error, observer coordinates do not match the frames.");

    // One frame axis spanning every observer's frames.  An observer that heard
    // nothing contributes none, and a run where none did writes an empty axis.
    std::int64_t first = std::numeric_limits<std::int64_t>::max();
    std::int64_t end = std::numeric_limits<std::int64_t>::min();
    for (const auto& f : frames)
    {
        if (f.count == 0)
            continue;
        first = std::min(first, f.first);
        end = std::max(end, f.first + static_cast<std::int64_t>(f.count));
    }
    const std::size_t frame_count = end > first ? static_cast<std::size_t>(end - first) : 0;
    if (frame_count == 0)
        first = 0;

    ensure_parent_directory(path);
    netcdf::File file;
    nc_check(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file.id), "creating it", path);
    const int ncid = file.id;

    int observer_dim = 0, frame_dim = 0, band_dim = 0;
    nc_check(nc_def_dim(ncid, "observer", observers, &observer_dim), "defining observer", path);
    nc_check(nc_def_dim(ncid, "frame", frame_count, &frame_dim), "defining frame", path);
    nc_check(nc_def_dim(ncid, "band", bands, &band_dim), "defining band", path);

    auto define = [&](const char* name, nc_type type, std::vector<int> dims, const char* units,
                          const std::string& long_name) {
        int var = 0;
        nc_check(nc_def_var(ncid, name, type, static_cast<int>(dims.size()), dims.data(), &var),
                std::string("defining ") + name, path);
        if (units != nullptr)
            put_text(ncid, var, "units", units, path);
        put_text(ncid, var, "long_name", long_name, path);
        return var;
    };
    const int time_var = define("time", NC_DOUBLE, {frame_dim}, "s",
            "center of the reception-time frame, on the trajectory's clock");
    const int band_var = define("band_center", NC_DOUBLE, {band_dim}, "Hz",
            "one-third-octave band center, as the database gives it");
    const std::string length = meta.length_unit;
    const int x_var = define("x", NC_DOUBLE, {observer_dim}, length.c_str(), "microphone x");
    const int y_var = define("y", NC_DOUBLE, {observer_dim}, length.c_str(), "microphone y");
    const int z_var = define("z", NC_DOUBLE, {observer_dim}, length.c_str(),
            "microphone z (the observer's ground point raised by its height)");
    const int height_var = define("microphone_height", NC_DOUBLE, {observer_dim}, length.c_str(),
            "microphone height above the ground; 0 is flush (on a plate with --plate_table)");
    int grid_x_var = -1, grid_y_var = -1;
    const bool grid = !meta.grid_x.empty();
    if (grid)
    {
        int gx = 0, gy = 0;
        nc_check(nc_def_dim(ncid, "grid_x", meta.grid_x.size(), &gx), "defining grid_x", path);
        nc_check(nc_def_dim(ncid, "grid_y", meta.grid_y.size(), &gy), "defining grid_y", path);
        grid_x_var = define("grid_x", NC_DOUBLE, {gx}, length.c_str(), "observer grid x axis");
        grid_y_var = define("grid_y", NC_DOUBLE, {gy}, length.c_str(), "observer grid y axis");
    }
    const int level_var = define("level", NC_FLOAT, {observer_dim, frame_dim, band_dim}, "dB re 20 uPa",
            "unweighted band level at the microphone, linear average over the frame");
    const float fill = std::numeric_limits<float>::quiet_NaN();
    nc_check(nc_def_var_fill(ncid, level_var, 0, &fill), "setting the fill value", path);
    put_text(ncid, level_var, "comment",
            "NaN outside the frames this observer received; -inf where a frame received nothing in "
            "a band.", path);
    if (frame_count > 0 && bands > 0)
    {
        // One observer's whole record per chunk: what a per-observer reader
        // wants, and small (frames x bands x 4 bytes).  Shuffled and deflated
        // at level 1, which costs little and roughly halves a grid's file.
        const std::size_t chunk[3] = {1, frame_count, bands};
        nc_check(nc_def_var_chunking(ncid, level_var, NC_CHUNKED, chunk), "setting the chunking", path);
        nc_check(nc_def_var_deflate(ncid, level_var, 1, 1, 1), "setting compression", path);
    }

    put_text(ncid, NC_GLOBAL, "title", "NICEOPS band levels in reception-time frames", path);
    put_text(ncid, NC_GLOBAL, "trajectory", meta.trajectory, path);
    put_text(ncid, NC_GLOBAL, "database", meta.database, path);
    put_text(ncid, NC_GLOBAL, "observer_layout",
            grid ? "grid: observer = iy * size(grid_x) + ix" : "scattered, in observer-file order", path);
    put_text(ncid, NC_GLOBAL, "time_base",
            "reception time: a state emitted at t is heard at t + r/c, and its emission step is "
            "spread over the reception time it arrives in", path);
    const int version = 1;
    nc_check(nc_put_att_int(ncid, NC_GLOBAL, "band_frames_version", NC_INT, 1, &version),
            "writing band_frames_version", path);
    nc_check(nc_put_att_double(ncid, NC_GLOBAL, "frame_length_seconds", NC_DOUBLE, 1, &meta.frame_length),
            "writing frame_length_seconds", path);
    nc_check(nc_enddef(ncid), "ending its definition", path);

    std::vector<double> time(frame_count);
    for (std::size_t k = 0; k < frame_count; k++)
        time[k] = (static_cast<double>(first + static_cast<std::int64_t>(k)) + 0.5) * meta.frame_length;
    if (frame_count > 0)
        nc_check(nc_put_var_double(ncid, time_var, time.data()), "writing time", path);
    if (bands > 0)
        nc_check(nc_put_var_double(ncid, band_var, meta.band_centers.data()), "writing band_center", path);
    if (observers > 0)
    {
        nc_check(nc_put_var_double(ncid, x_var, meta.x.data()), "writing x", path);
        nc_check(nc_put_var_double(ncid, y_var, meta.y.data()), "writing y", path);
        nc_check(nc_put_var_double(ncid, z_var, meta.z.data()), "writing z", path);
        nc_check(nc_put_var_double(ncid, height_var, meta.height.data()), "writing microphone_height", path);
    }
    if (grid)
    {
        nc_check(nc_put_var_double(ncid, grid_x_var, meta.grid_x.data()), "writing grid_x", path);
        nc_check(nc_put_var_double(ncid, grid_y_var, meta.grid_y.data()), "writing grid_y", path);
    }
    // Observer by observer, each a whole row of the frame axis with its own
    // frames in place and NaN around them.
    std::vector<float> row(frame_count * bands);
    for (std::size_t i = 0; i < observers && frame_count > 0 && bands > 0; i++)
    {
        std::fill(row.begin(), row.end(), fill);
        const auto& f = frames[i];
        if (f.levels.size() != f.count * bands)
            input_error("--band_frames ", path, ": internal error, observer ", i, "'s frames hold ",
                    f.levels.size(), " levels for ", f.count, " frames of ", bands, " bands.");
        if (f.count > 0)
        {
            const std::size_t offset = static_cast<std::size_t>(f.first - first) * bands;
            std::copy(f.levels.begin(), f.levels.end(), row.begin() + static_cast<std::ptrdiff_t>(offset));
        }
        const std::size_t start[3] = {i, 0, 0};
        const std::size_t count[3] = {1, frame_count, bands};
        nc_check(nc_put_vara_float(ncid, level_var, start, count, row.data()), "writing level", path);
    }
    nc_check(file.close(), "closing it", path);
    info("Wrote ", observers, " observers' band levels in ", frame_count, " frames of ", meta.frame_length,
            " s to ", path, ".");
}

} // namespace niceops
