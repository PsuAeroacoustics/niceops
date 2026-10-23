#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "buildings.hpp"
#include "diffuse_facades.hpp"
#include "gpu/gpu.hpp"
#include "point_paths.hpp"
#include "source_coherence.hpp"
#include "sphere_database.hpp"
#include "urban_point.hpp"
#include "weighting.hpp"

using namespace niceops; // NOLINT

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

// A single-condition database whose every cell carries the same flat
// spectrum: a source the same in every direction.
std::string write_omnidirectional_database(const std::vector<double>& bands, double level, const std::string& name)
{
    Database plain("databases/S-76D_M3.nod");
    auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    sphere.frequencies = bands;
    sphere.amplitudes.assign(sphere.cells.size() * bands.size(), level);
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    write_sphere(sphere, path);
    return path;
}

} // namespace

// The lighting paths carry the source's power along their own launch: a
// source whose power falls as cos^2 off the wall's normal, beside the wall of
// "a point source lights a wall's diffuse field as Lambert's law says", gives
// rho s d^4 / pi times the integral of dA / r^8.  A power the same in every
// direction scales the unit source's field.
TEST_CASE("a directional source lights the diffuse field along each path's launch", "[diffuse][urban_point]")
{
    const double d = 1.0;
    const Buildings city({box(0.0, -5.0, 10.0, 5.0, 10.0)});
    const Vec3 point(-d, 0.0, 5.0);
    const std::vector<double> bands {1000.0};
    const std::vector<double> rho(city.facades().size(), 0.9), s(city.facades().size(), 0.5);
    DiffuseOptions options;
    options.radius = 50.0;
    options.patch_size = 0.2;
    const DiffuseFacades diffuse(city, point, rho, s, 0.0, bands, 343.2, options, 4);

    const DiffuseFacades::SourcePower toward_wall = [](const Vec3& launch, double, double* power) {
        power[0] = launch.x * launch.x;
    };
    double energy = 0.0;
    diffuse.point_energy(point, 4, {}, &energy, 4, &toward_wall);
    double integral = 0.0;
    const int n = 1000;
    const double step = 10.0 / n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            const double y = -5.0 + (i + 0.5) * step, z = -5.0 + (j + 0.5) * step;
            integral += step * step / std::pow(d * d + y * y + z * z, 4);
        }
    REQUIRE(energy == Catch::Approx(0.9 * 0.5 * std::pow(d, 4) / pi * integral).epsilon(5e-3));

    double unit = 0.0, scaled = 0.0;
    const DiffuseFacades::SourcePower constant = [](const Vec3&, double, double* power) { power[0] = 7.0; };
    diffuse.point_energy(point, 4, {}, &unit, 4);
    diffuse.point_energy(point, 4, {}, &scaled, 4, &constant);
    REQUIRE(scaled == Catch::Approx(7.0 * unit).epsilon(1e-12));
}

// A run's rays and diffuse field (UrbanDistrict::levels) against
// niceops_point_field's for a source the same in every direction: the rays
// carry rho (1 - s) of the walls each microphone's field covers and rho of the
// rest, the field adds its energy band by band, and the source's level at
// 1 m scales both.  The near building scatters; the far one is out of the
// diffuse field's reach but still reflects.
TEST_CASE("a run's diffuse field is niceops_point_field's for an omnidirectional source", "[urban_point]")
{
    const std::vector<double> bands {125.0, 500.0, 2000.0};
    const double level = 90.0;
    const Database database(write_omnidirectional_database(bands, level, "niceops_omnidirectional_urban.nod"));
    const double radius = database.sphere_radius();
    REQUIRE(radius > 0.0);

    // Two buildings, the source well beyond the sphere's radius of either, so
    // that no ray or lighting path starts inside it.
    const std::vector<Footprint> footprints {box(0.0, -10.0, 20.0, 10.0, 15.0), box(-40.0, 120.0, -10.0, 150.0, 25.0)};
    const auto district_file = (std::filesystem::temp_directory_path() / "niceops_urban_scattering.csv").string();
    {
        std::ofstream out(district_file);
        out << "building,x,y,top,absorption,scattering\n";
        for (std::size_t b = 0; b < footprints.size(); b++)
            for (std::size_t v = 0; v < footprints[b].x.size(); v++)
                out << b + 1 << "," << footprints[b].x[v] << "," << footprints[b].y[v] << "," << footprints[b].top
                    << ",0.1," << (b == 0 ? 0.4 : 0.3) << "\n";
    }
    UrbanDistrict::Input input;
    input.units = "m";
    input.to_meters = 1.0;
    input.patch_size = 2.0;
    input.diffuse_radius = 60.0;
    const auto district = UrbanDistrict::read(district_file, input);
    REQUIRE(district->scatters());

    const Vec3 source(-60.0 - radius, 40.0, 45.0 + radius);
    const Vec3 microphone(-12.0, 4.0, 1.5);
    const std::vector<Database::SourceState> states {
            database.source_state(source, Vec3(30.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), 0.15, -4.0, 0.006)};
    const std::vector<double> no_absorption(bands.size(), 0.0);
    const std::vector<double> gains = a_weighting_gains(bands, BandLimits {});
    SourceCoherence coherence;
    coherence.sum = RaySum::coherent;
    std::vector<double> levels, band_levels;
    district->levels(database, states, {microphone}, no_absorption.data(), gains, 0, bands.size(), 343.2, 2,
            coherence, levels, &band_levels, nullptr, nullptr, 0.0);

    // point_field's way, for a unit source.
    const Buildings city(footprints, 1.0 * pi / 180.0);
    const std::size_t nb = bands.size(), nf = city.facades().size();
    std::vector<double> rho(nf * nb), s(nf * nb), reflection(nf * nb);
    for (std::size_t f = 0; f < nf; f++)
        for (std::size_t b = 0; b < nb; b++)
        {
            rho[f * nb + b] = 0.9;
            s[f * nb + b] = city.facade_building(f) == 0 ? 0.4 : 0.3;
        }
    DiffuseOptions options;
    options.radius = 60.0;
    options.patch_size = 2.0;
    options.roofs = true;
    const DiffuseFacades diffuse(city, microphone, rho, s, 1.0, bands, 343.2, options, 2);
    std::size_t covered = 0;
    for (std::size_t f = 0; f < nf; f++)
    {
        covered += diffuse.scatters(f) ? 1 : 0;
        for (std::size_t b = 0; b < nb; b++)
            reflection[f * nb + b] = diffuse.scatters(f) ? rho[f * nb + b] * (1.0 - s[f * nb + b]) : rho[f * nb + b];
    }
    REQUIRE(covered > 0);
    REQUIRE(covered < nf);
    std::vector<double> energy(nb), rays(nb);
    diffuse.point_energy(source, PointPathOptions {}.max_reflections, {}, energy.data());
    PointPathOptions paths;
    paths.sample_shared_faces = true;
    PointSourceField(city, source, microphone, paths).band_levels(reflection, 1.0, bands, 343.2, rays.data());
    const double at_one_metre = level + 20.0 * std::log10(radius);
    for (std::size_t b = 0; b < nb; b++)
    {
        REQUIRE(energy[b] > 0.0);
        const double expected = 10.0 * std::log10(std::pow(10.0, 0.1 * rays[b]) + energy[b]) + at_one_metre;
        REQUIRE(band_levels[b] == Catch::Approx(expected).margin(1e-9));
        // The diffuse field is not negligible here, so the test can see it.
        REQUIRE(std::abs(expected - rays[b] - at_one_metre) > 0.01);
    }

    // On a GPU the rays are summed there and the diffuse energy added as each
    // pair's energies return: the same levels, to the GPU's round-off.
    if (gpu::Device* device = gpu::Device::get())
    {
        std::vector<double> gpu_levels, gpu_bands;
        district->levels(database, states, {microphone}, no_absorption.data(), gains, 0, bands.size(), 343.2, 2,
                coherence, gpu_levels, &gpu_bands, nullptr, nullptr, 0.0, device);
        for (std::size_t b = 0; b < nb; b++)
            REQUIRE(gpu_bands[b] == Catch::Approx(band_levels[b]).margin(1e-6));
        REQUIRE(gpu_levels[0] == Catch::Approx(levels[0]).margin(1e-6));
    }
}

// Turbulence scattering in a run: behind a long wall the microphone gains the
// energy turbulence scatters over it, more at high frequency than low, and
// never past open ground's level (the direct path's free field and the
// floor's reflection); a microphone in sight of the source is untouched, to
// the bit.
TEST_CASE("a run's shadows take the energy turbulence scatters into them", "[urban_point]")
{
    const std::vector<double> bands {125.0, 500.0, 2000.0};
    const double level = 90.0;
    const Database database(write_omnidirectional_database(bands, level, "niceops_omnidirectional_turbulence.nod"));
    const double radius = database.sphere_radius();
    const auto district_file = (std::filesystem::temp_directory_path() / "niceops_turbulence_wall.csv").string();
    {
        std::ofstream out(district_file);
        out << "building,x,y,top\n1,-1,-200,10\n1,1,-200,10\n1,1,200,10\n1,-1,200,10\n";
    }
    UrbanDistrict::Input still, turbulent;
    still.units = turbulent.units = "m";
    still.to_meters = turbulent.to_meters = 1.0;
    turbulent.turbulence = 1e-5;
    const auto quiet = UrbanDistrict::read(district_file, still), stirred = UrbanDistrict::read(district_file, turbulent);
    const Vec3 source(-200.0, 0.0, 30.0 + radius);
    const std::vector<Vec3> microphones {Vec3(15.0, 0.0, 1.5), Vec3(-15.0, 0.0, 1.5)};
    const std::vector<Database::SourceState> states {
            database.source_state(source, Vec3(30.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), 0.15, -4.0, 0.006)};
    const std::vector<double> no_absorption(bands.size(), 0.0);
    const std::vector<double> gains = a_weighting_gains(bands, BandLimits {});
    SourceCoherence coherence;
    std::vector<double> levels_quiet, bands_quiet, levels_stirred, bands_stirred;
    quiet->levels(database, states, microphones, no_absorption.data(), gains, 0, bands.size(), 343.2, 2, coherence,
            levels_quiet, &bands_quiet, nullptr, nullptr, 0.0);
    stirred->levels(database, states, microphones, no_absorption.data(), gains, 0, bands.size(), 343.2, 2, coherence,
            levels_stirred, &bands_stirred, nullptr, nullptr, 0.0);
    const double r = (microphones[0] - source).norm();
    const double open = level + 20.0 * std::log10(radius) - 20.0 * std::log10(r) + 10.0 * std::log10(2.0);
    std::vector<double> gain(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        gain[b] = bands_stirred[b] - bands_quiet[b];
        REQUIRE(gain[b] >= 0.0);
        REQUIRE(bands_stirred[b] <= open + 1e-9);
        // In sight of the source nothing changes.
        REQUIRE(bands_stirred[bands.size() + b] == bands_quiet[bands.size() + b]);
    }
    REQUIRE(gain[2] > 0.3);
    REQUIRE(gain[2] > gain[0]);
    REQUIRE(levels_stirred[1] == levels_quiet[1]);
}

// A tilted step's split is the difference of the foot fields' twin rays, which
// only cancel where the rays add with their phases.  Summed in energy the
// twins would each add, so a run takes no split there: the same levels as
// with --no_tilted_step_split, to the bit, where the coherent sum moves.
TEST_CASE("a run summed in energy takes no tilted-step split", "[urban_point][point_paths]")
{
    const std::vector<double> bands {125.0, 250.0, 500.0};
    const Database database(write_omnidirectional_database(bands, 90.0, "niceops_omnidirectional_step.nod"));
    const double radius = database.sphere_radius();
    // The step of "a tilted step's split does not step": 10 m then 7 m.
    const auto district_file = (std::filesystem::temp_directory_path() / "niceops_tilted_step.csv").string();
    {
        std::ofstream out(district_file);
        out << "building,x,y,top\n1,-400,0,10\n1,400,0,10\n1,400,10,10\n1,-400,10,10\n"
            << "2,-400,10,7\n2,400,10,7\n2,400,20,7\n2,-400,20,7\n";
    }
    UrbanDistrict::Input split, unsplit;
    split.units = unsplit.units = "m";
    split.to_meters = unsplit.to_meters = 1.0;
    unsplit.no_tilted_step_split = true;
    const auto with = UrbanDistrict::read(district_file, split), without = UrbanDistrict::read(district_file, unsplit);
    // The source 3.3 degrees above the pair's line, the microphones either
    // side of it 2 m past the lower block.
    const double a = std::atan(0.3) + 3.3 * pi / 180.0, r = 60.0 + radius;
    const Vec3 source(0.0, 10.0 - r * std::cos(a), 10.0 + r * std::sin(a));
    const std::vector<Vec3> microphones {Vec3(0.0, 22.0, 6.45), Vec3(0.0, 22.0, 6.35)};
    const std::vector<Database::SourceState> states {
            database.source_state(source, Vec3(30.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0), 0.15, -4.0, 0.006)};
    const std::vector<double> no_absorption(bands.size(), 0.0);
    const std::vector<double> gains = a_weighting_gains(bands, BandLimits {});
    for (const RaySum sum : {RaySum::energy, RaySum::coherent})
    {
        SourceCoherence coherence;
        coherence.sum = sum;
        std::vector<double> levels_with, bands_with, levels_without, bands_without;
        with->levels(database, states, microphones, no_absorption.data(), gains, 0, bands.size(), 343.2, 1,
                coherence, levels_with, &bands_with, nullptr, nullptr, 0.0);
        without->levels(database, states, microphones, no_absorption.data(), gains, 0, bands.size(), 343.2, 1,
                coherence, levels_without, &bands_without, nullptr, nullptr, 0.0);
        double moved = 0.0;
        for (std::size_t i = 0; i < bands_with.size(); i++)
            moved = std::max(moved, std::abs(bands_with[i] - bands_without[i]));
        INFO((sum == RaySum::energy ? "energy" : "coherent") << ": largest band change " << moved << " dB");
        if (sum == RaySum::energy)
            REQUIRE(bands_with == bands_without);
        else
            REQUIRE(moved > 0.1);
    }
}
