// The refracted ray's geometry for a batch of source-receiver pairs.
//
//   niceops_ray_geometry --atmosphere profile.csv --frame_bearing 270 --refraction stratified
//                        --ground point --paths paths.csv > rays.csv
//
// Each row of --paths is `id,dx,dy,source_height,receiver_height` (id numeric,
// optional; the row number without it): the source's
// position relative to the receiver in the run's local frame (+x along
// --frame_bearing, as niceops takes it), and both heights above the ground.
// For each, writes what niceops would hand its ground model and sphere lookup:
//
//   lit           1 if a ray reaches the receiver, 0 in a shadow
//   launch_deg    elevation at the source, degrees above the horizontal
//   arrival_deg   elevation at the receiver (negative descending)
//   spreading_db  20 lg(straight-line distance / ray-tube range): the ray's gain
//                 over spherical spreading
//   shadow        distance past the shadow boundary, in --units (0 when lit)
//   boundary      the shadow boundary's range, in --units (inf where none)
//   from_fit      1 where the lit ray is the linear fit's (past the hand-over)
//   bounces       ground bounces before the primary's last descent (a duct)
//   extras        other eigenrays (layered_arrivals), and their energy
//   extra_db      relative to the primary's, dB (-inf with none)
//   source_height, ground_distance
//                 the equivalent straight geometry the ground model takes, in
//                 --units: from the reflected ray for --ground point, from the
//                 direct ray's arrival for --ground plate
//
// The ground model, the transition near a shadow boundary and absorption are
// not applied: this is the geometry only, for evaluating ground models outside
// niceops on the same rays (the validation's board/pole pairs).  Pairs are
// independent, so they are spread over the threads.

#include <CLI/CLI.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/parallel.hpp"

using namespace niceops; // NOLINT

namespace {

struct Pair
{
    double id = 0.0;
    double dx = 0.0, dy = 0.0, z_s = 0.0, z_r = 0.0;
};

std::vector<Pair> read_pairs(const std::string& path, double to_metres)
{
    const auto table = read_labelled_csv(path);
    auto column = [&](const char* name) {
        const std::size_t k = table.column({name});
        if (k == LabelledCsv::npos)
            fatal("--paths ", path, " has no `", name, "` column");
        return k;
    };
    const std::size_t dx = column("dx"), dy = column("dy"), zs = column("source_height"),
                      zr = column("receiver_height");
    const std::size_t id = table.column({"id"});
    std::vector<Pair> out;
    for (std::size_t row = 0; row < table.rows.size(); row++)
    {
        const auto& r = table.rows[row];
        Pair p;
        p.id = id != LabelledCsv::npos ? r[id] : static_cast<double>(row);
        p.dx = r[dx] * to_metres;
        p.dy = r[dy] * to_metres;
        p.z_s = r[zs] * to_metres;
        p.z_r = r[zr] * to_metres;
        out.push_back(p);
    }
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    CLI::App app {"The refracted ray's geometry for a batch of source-receiver pairs."};
    std::string atmosphere_file, paths_file;
    std::string units = "ft";
    std::string refraction = "stratified";
    std::string ground = "point";
    double frame_bearing = 90.0;
    int threads = 0;
    app.add_option("--atmosphere", atmosphere_file, "Atmosphere profile, as for niceops.")->required();
    app.add_option("--paths", paths_file, "CSV of id,dx,dy,source_height,receiver_height.")->required();
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--frame_bearing", frame_bearing, "Compass bearing of the paths' +x axis (default 90).");
    app.add_option("--refraction", refraction, "linear or stratified (default).")
            ->check(CLI::IsMember({"linear", "stratified"}));
    app.add_option("--ground", ground, "point (default): the reflected ray's geometry; plate: the direct ray's.")
            ->check(CLI::IsMember({"point", "plate"}));
    app.add_option("-j,--threads", threads, "Threads (default: all).");
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/ray_geometry.log");
    const LengthUnit unit = units == "m" ? metres : feet;
    const auto pairs = read_pairs(paths_file, unit.to_metres);

    auto profile = read_atmosphere_profile(atmosphere_file, unit, 0.0);
    if (!profile.has_pressure())
        profile = profile.with_hydrostatic_pressure(reference_pressure);
    double z_top = profile.top();
    for (const Pair& p : pairs)
        z_top = std::max(z_top, std::max(p.z_s, p.z_r));
    z_top += 2.0 * EffectiveSoundSpeed::min_fit_span;
    const EffectiveSoundSpeed linear(profile, frame_bearing, 0.0, z_top);
    const bool layered = refraction == "stratified";
    std::unique_ptr<LayeredSoundSpeed> layers;
    if (layered)
        layers = std::make_unique<LayeredSoundSpeed>(profile, frame_bearing);
    const GroundGeometry geometry = ground == "plate" ? GroundGeometry::plate : GroundGeometry::point;

    std::vector<std::string> lines(pairs.size());
    const double inf = std::numeric_limits<double>::infinity();
    const unsigned workers = threads > 0 ? static_cast<unsigned>(threads) : default_thread_count();
    parallel_chunks(pairs.size(), workers, [&](std::size_t begin, std::size_t end) {
        std::vector<double> scratch;
        for (std::size_t k = begin; k < end; k++)
        {
            const Pair& p = pairs[k];
            const Vec3 receiver(0.0, 0.0, p.z_r);
            const Vec3 source(p.dx, p.dy, p.z_s);
            const auto at_source = linear.moments(p.z_s);
            const auto at_receiver = linear.moments(p.z_r);
            RefractedRay ray;
            if (layered)
            {
                ray = layered_direct_ray(*layers, linear, source, at_source, receiver, at_receiver, scratch, true);
                layered_finish_ground(*layers, ray, geometry, p.z_r, scratch);
            }
            else
            {
                ray = direct_ray(linear, source, at_source, receiver, at_receiver);
                finish_ground(ray, geometry, p.z_r);
            }
            const double slant = (source - receiver).norm();
            const bool curved = ray.curved;
            const double launch = curved ? ray.direct.launch : std::asin((p.z_r - p.z_s) / slant);
            const double arrival = curved ? ray.direct.arrival : launch;
            double extra = 0.0;
            for (const Arrival& a : ray.extra)
                extra += (ray.path.range / a.path.range) * (ray.path.range / a.path.range);
            std::ostringstream line;
            line << std::setprecision(10) << p.id << "," << (ray.shadow_depth > 0.0 ? 0 : 1) << ","
                 << launch * 180.0 / pi << "," << arrival * 180.0 / pi << ","
                 << 20.0 * std::log10(slant / ray.path.range) << "," << ray.shadow_depth / unit.to_metres << ","
                 << (curved ? ray.direct.boundary / unit.to_metres : inf) << "," << (ray.from_linear_fit ? 1 : 0)
                 << "," << ray.bounces << "," << ray.extra.size() << ","
                 << (extra > 0.0 ? 10.0 * std::log10(extra) : -inf) << ","
                 << ray.path.source_height / unit.to_metres << "," << ray.path.ground_distance / unit.to_metres;
            lines[k] = line.str();
        }
    });
    std::cout << "id,lit,launch_deg,arrival_deg,spreading_db,shadow,boundary,from_fit,bounces,extras,extra_db,"
                 "source_height,ground_distance\n";
    for (const auto& line : lines)
        std::cout << line << "\n";
    return 0;
}
