// Refraction models against the parabolic equation, on one path.
//
//   niceops_pe_compare --atmosphere profile.csv --toward 270 --source_height 328
//                      --receiver_height 3.937 --bands 125 250 500 1000
//                      --max_range 8000 [--flow_resistance 200 ...] > out.csv
//
// For a source and a receiver in the vertical plane heading along a compass
// bearing, writes each band's level relative to a point source's free field
// at the straight-line distance -- so spreading, refraction and the ground,
// and nothing about the source -- by four models:
//   straight     straight rays and the spherical-wave ground model
//   linear       --refraction linear's rays and the same ground model
//   stratified   --refraction stratified's rays and the same ground model
//   pe           the wide-angle parabolic equation (pe.hpp), band-averaged
// Absorption is left out of all four, as is turbulence: the PE has neither, so
// the ground model runs with a coherent atmosphere too.  Lengths are in
// --units (feet by default), bearings in degrees.

#include <CLI/CLI.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "ground_effect.hpp"
#include "pe.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"Refraction models against the parabolic equation, on one path."};
    std::string atmosphere_file;
    std::string units = "ft";
    double toward = 0.0;
    double source_height = 0.0;
    double receiver_height = 0.0;
    std::vector<double> bands = {125.0, 250.0, 500.0, 1000.0};
    double max_range = 0.0;
    double step = 0.0;
    double flow_resistance = flow_resistance::grassy_field;
    std::string impedance_model = "delany_bazley";
    double porosity_rate = 0.0;
    int samples = 5;
    double points_per_wavelength = 10.0;
    double range_points_per_wavelength = 10.0;
    app.add_option("--atmosphere", atmosphere_file, "Atmosphere profile, as for niceops.")->required();
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--toward", toward, "Compass bearing the sound travels along, degrees.")->required();
    app.add_option("--source_height", source_height, "Source height, in --units.")->required();
    app.add_option("--receiver_height", receiver_height, "Receiver height, in --units.")->required();
    app.add_option("--bands", bands, "One-third-octave band centres, Hz.");
    app.add_option("--max_range", max_range, "Farthest range, in --units.")->required();
    app.add_option("--step", step, "Output range spacing, in --units (default max_range / 100).");
    app.add_option("--flow_resistance", flow_resistance, "Ground flow resistance, kPa s/m^2.");
    app.add_option("--impedance_model", impedance_model, "As for niceops.");
    app.add_option("--porosity_rate", porosity_rate, "As for niceops, 1/m.");
    app.add_option("--samples", samples, "PE frequencies per band (default 5).");
    app.add_option("--ppw", points_per_wavelength, "PE grid points per wavelength in height (default 10).");
    app.add_option("--range_ppw", range_points_per_wavelength,
            "PE range steps per wavelength (default 10).");
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/pe_compare.log");
    const LengthUnit unit = units == "m" ? metres : feet;
    const double z_s = source_height * unit.to_metres;
    const double z_r = receiver_height * unit.to_metres;
    const double r_max = max_range * unit.to_metres;
    const double dr_out = (step > 0.0 ? step * unit.to_metres : r_max / 100.0);

    auto profile = read_atmosphere_profile(atmosphere_file, unit, 0.0);
    if (!profile.has_pressure())
        profile = profile.with_hydrostatic_pressure(reference_pressure);
    const double sound_speed = profile.at(0.0).sound_speed();

    // The plane's effective sound speed at the profile's levels: local frame x
    // east, y north, and the propagation direction's components in it.
    const double bearing = toward * pi / 180.0;
    const double hx = std::sin(bearing);
    const double hy = std::cos(bearing);
    const LayeredSoundSpeed layered(profile, 90.0);
    const double z_top = std::max(profile.top(), z_s) + 2.0 * EffectiveSoundSpeed::min_fit_span;
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0, z_top);
    std::vector<double> c(layered.levels());
    layered.along(hx, hy, c.data());

    GroundModel ground;
    ground.flow_resistance = flow_resistance;
    ground.impedance_model = parse_impedance_model(impedance_model);
    ground.porosity_rate = porosity_rate;
    ground.receiver_height = z_r;
    const auto admittance = ground_impedance(ground, bands, sound_speed);
    const ShadowAttenuation attenuation(bands, admittance, 0, bands.size());
    // The receiver is a point receiver: near the boundary its level comes from
    // the transition's field, as in niceops.
    const ShadowTransition transition(bands, admittance, z_r, sound_speed, 0, bands.size());

    std::vector<PeResult> pe(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        PeSetup setup;
        setup.source_height = z_s;
        setup.receiver_heights = {z_r};
        setup.max_range = r_max;
        setup.output_step = dr_out;
        setup.admittance = admittance[b];
        setup.z = layered.heights();
        setup.c = c;
        setup.points_per_wavelength = points_per_wavelength;
        setup.range_points_per_wavelength = range_points_per_wavelength;
        std::cerr << "PE at " << bands[b] << " Hz..." << std::endl;
        pe[b] = solve_pe_band(setup, bands[b], samples);
    }

    std::cout << "range_" << unit.name;
    for (const double f : bands)
        std::cout << ",straight_" << f << ",linear_" << f << ",stratified_" << f << ",pe_" << f;
    std::cout << "\n" << std::setprecision(6);

    const Vec3 source(0.0, 0.0, z_s);
    const auto at_source = linear.moments(z_s);
    const auto at_receiver = linear.moments(z_r);
    std::vector<double> scratch;
    std::vector<double> loss(bands.size()), weight(bands.size());
    for (std::size_t q = 0; q < pe[0].range.size(); q++)
    {
        const double r = pe[0].range[q];
        const Vec3 receiver(r * hx, r * hy, z_r);
        const double slant = std::hypot(r, z_s - z_r);
        auto relative = [&](const RefractedRay& ray, std::size_t b) {
            const double ground_db = excess_ground_attenuation(ground, ray.path.source_height,
                    ray.path.ground_distance, sound_speed, {bands[b]})[0];
            const double spreading_db = 20.0 * std::log10(slant / ray.path.range);
            // A weight of 0 drops the ground term, which at the surface can be
            // -infinity; 0 * -inf would be NaN.
            if (transition.apply(ray, attenuation, 20.0, 0, bands.size(), weight.data(), loss.data()))
                return spreading_db + (weight[b] > 0.0 ? weight[b] * ground_db : 0.0) - loss[b];
            double shadow_db = 0.0;
            if (ray.shadow_depth > 0.0)
            {
                shadow_loss(ray, attenuation, 20.0, 0, bands.size(), loss.data());
                shadow_db = loss[b];
            }
            return spreading_db + ground_db - shadow_db;
        };
        RefractedRay straight;
        straight.path = straight_ray(source, receiver);
        RefractedRay by_linear = direct_ray(linear, source, at_source, receiver, at_receiver);
        finish_ground(by_linear, GroundGeometry::point, z_r);
        RefractedRay by_layers = layered_direct_ray(layered, linear, source, at_source, receiver,
                at_receiver, scratch, true);
        layered_finish_ground(layered, by_layers, GroundGeometry::point, z_r, scratch);

        std::cout << r / unit.to_metres;
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            std::cout << "," << relative(straight, b) << "," << relative(by_linear, b) << ","
                      << relative(by_layers, b) << "," << pe[b].relative_level[0][q];
        }
        std::cout << "\n";
    }
    return 0;
}
