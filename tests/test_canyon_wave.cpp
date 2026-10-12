#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "canyon_wave.hpp"
#include "urban.hpp"
#include "weighting.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

const std::vector<double> third_octaves {40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0,
        500.0};

double edge(double center, int side)
{
    return center * std::pow(2.0, side / 6.0);
}

} // namespace

// The solver, its normalization and the far-field transform together: over
// flat rigid ground the transform of a horizontal recording line is the image
// solution, G = 4 cos^2(k z sin phi), band by band, from 45 to 135 degrees.
// Nearer grazing the recording line itself errs (0.2-0.8 dB at 15-30 degrees,
// however long), which the canyon's transform, over the opening alone, does
// not share.
TEST_CASE("the open-ground far field from a recording line is the image solution", "[canyon_wave]")
{
    CanyonWaveSetup s;
    s.max_frequency = 200.0;
    s.points_per_wavelength = 12.0;
    s.microphone_height = 1.5;
    s.roof_margin = 5.0;
    s.sky_margin = 5.0;
    const ApertureSpectrum a = solve_open_ground(s, 3.0, 120.0);
    REQUIRE(a.residual_db < -40.0);
    double worst = 0.0;
    for (double center : third_octaves)
    {
        if (edge(center, 1) > s.max_frequency)
            continue;
        for (double degrees = 45.0; degrees <= 135.0; degrees += 15.0)
        {
            const double ratio = a.band_ratio(edge(center, -1), edge(center, 1), degrees * pi / 180.0);
            worst = std::max(worst, std::abs(10.0 * std::log10(ratio)));
        }
    }
    INFO("largest band error " << worst << " dB");
    REQUIRE(worst < 0.15);
}

// The PML and the margins: moving the grid's edges must not move the answer.
TEST_CASE("a canyon's far field does not depend on how far the grid extends", "[canyon_wave]")
{
    CanyonWaveSetup s;
    s.width = 10.0;
    s.depth = 10.0;
    s.offset = 3.0;
    s.microphone_height = 1.5;
    s.max_frequency = 150.0;
    s.points_per_wavelength = 12.0;
    s.roof_margin = 4.0;
    s.sky_margin = 4.0;
    const ApertureSpectrum near = solve_canyon(s);
    s.roof_margin = 10.0;
    s.sky_margin = 10.0;
    const ApertureSpectrum far = solve_canyon(s);
    REQUIRE(near.residual_db < -28.0);
    double worst = 0.0;
    for (double center : {50.0, 63.0, 80.0, 100.0, 125.0})
        for (double degrees = 10.0; degrees <= 170.0; degrees += 20.0)
        {
            const double phi = degrees * pi / 180.0;
            const double a = near.band_ratio(edge(center, -1), edge(center, 1), phi);
            const double b = far.band_ratio(edge(center, -1), edge(center, 1), phi);
            worst = std::max(worst, std::abs(10.0 * std::log10(a / b)));
        }
    INFO("largest band difference " << worst << " dB");
    REQUIRE(worst < 0.1);
}

// Not run by default: the canyon's band ratios by grid resolution, and against
// the ideal canyon's specular table where the two should meet, printed.
TEST_CASE("canyon wave convergence and the geometric limit", "[.canyon_wave_study]")
{
    CanyonWaveSetup s;
    s.width = 20.0;
    s.depth = 20.0;
    s.offset = 6.0;
    s.microphone_height = 1.5;
    s.max_frequency = 630.0;
    s.threads = 8;
    std::printf("ppw  cells  steps  duration_s  residual_dB\n");
    std::vector<ApertureSpectrum> runs;
    for (double ppw : {10.0, 15.0, 20.0})
    {
        s.points_per_wavelength = ppw;
        runs.push_back(solve_canyon(s));
        std::printf("%4.0f %7zu %6zu %10.3f %12.1f\n", ppw, runs.back().cells, runs.back().steps,
                runs.back().duration, runs.back().residual_db);
    }
    const IdealCanyon canyon {s.width, s.offset, s.depth, s.depth, s.microphone_height};
        std::printf("band_Hz phi_deg  ppw10  ppw15  ppw20  geometric(incl. first order)\n");
    for (double center : {63.0, 125.0, 250.0, 500.0})
        for (double degrees : {15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 120.0, 150.0})
        {
            const double phi = degrees * pi / 180.0;
            std::printf("%7.0f %7.0f", center, degrees);
            for (const auto& r : runs)
                std::printf(" %6.2f", 10.0 * std::log10(r.band_ratio(edge(center, -1), edge(center, 1), phi)));
            double ratio = 0.0;
            const std::vector<double> band(1, center);
            // The wave run's walls absorb s.wall_absorption at normal incidence;
            // the geometric table takes the same as a specular loss.
            const std::vector<double> walls(1, 1.0 - s.wall_absorption);
            canyon_ratio(canyon, std::cos(phi), std::sin(phi), walls.data(), walls.data(), 1.0, band, s.sound_speed,
                    true, &ratio);
            // The first-order reflection, which the run adds itself: its ray
            // from above, and its floor partner when that is open too, in
            // units of open ground's two rays.
            IdealCanyon floor_image = canyon;
            floor_image.microphone_height = -canyon.microphone_height;
            const CanyonTrace mirror = trace_canyon(canyon, -std::cos(phi), std::sin(phi));
            const CanyonTrace mirror_floor = trace_canyon(floor_image, -std::cos(phi), std::sin(phi));
            if (mirror.escaped && mirror.bounces == 1)
                ratio += 0.5 * walls[0] * (mirror_floor.escaped && mirror_floor.bounces == 1 ? 2.0 : 1.0);
            std::printf(" %8.2f\n", 10.0 * std::log10(ratio));
        }
}

namespace {

// One harmonic's urban change over one pass, three ways (dB re open ground):
// the coherent canyon field; the coherent open ground times the band table
// (the table's band averaging alone); and what the run does, the band-averaged
// ground times the band table.
struct ToneErrors
{
    double truth = 0.0;
    double table = 0.0;
    double run = 0.0;
};

const std::vector<double> tone_bands {25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0,
        400.0, 500.0};

// The open ground's band average for a microphone z up, the plane wave
// arriving at sin(e) = u_z: 4 cos^2(k z u_z) over the band.
double open_band(double f_lo, double f_hi, double z, double u_z, double c)
{
    double sum = 0.0;
    const int n = 64;
    for (int k = 0; k < n; k++)
    {
        const double f = f_lo + (f_hi - f_lo) * (k + 0.5) / n;
        const double v = std::cos(2.0 * pi * f / c * z * u_z);
        sum += 4.0 * v * v;
    }
    return sum / n;
}

} // namespace

// Not run by default: does a band-energy table describe a canyon's effect on a
// rotor tone?  For each canyon class, microphone and pass, each harmonic is
// Doppler-swept through the coherent canyon response; its SEL change against
// open ground is compared with the band table's, read as the run reads it
// (bilinear in azimuth and elevation on the builder's grid).  Printed for the
// stage-2 gate (docs/FIDELITY.md).
TEST_CASE("rotor tones through a canyon against the band table", "[.canyon_tone_study]")
{
    const double c = 343.2, speed = 60.0, mach = speed / c;
    // B407-like: a four-bladed main rotor at 413 rpm, a two-bladed tail
    // rotor at 2380 rpm; harmonics kept below 460 Hz so that the Doppler
    // shift keeps them within the 500 Hz band.
    std::vector<std::pair<const char*, double>> harmonics;
    for (int m = 1; 27.53 * m <= 460.0; m++)
        harmonics.push_back({"main", 27.53 * m});
    for (int m = 1; 79.33 * m <= 460.0; m++)
        harmonics.push_back({"tail", 79.33 * m});

    UrbanTable::Axes axes;
    for (int k = 0; k < 72; k++)
        axes.azimuths.push_back(5.0 * k);
    for (int k = 0; k <= 20; k++)
        axes.elevations.push_back(1.5 * k);
    for (int k = 7; k <= 18; k++)
        axes.elevations.push_back(5.0 * k);
    axes.bands = tone_bands;

    std::printf("H/W  mic  doppler  harmonics x passes   |table err| median p90 max   |run err| median p90 max   "
                "A-wtd run err: median max\n");
    for (double depth : {10.0, 20.0, 40.0})
        for (double offset : {10.0, 5.0})
        {
            CanyonWaveSetup s;
            s.width = 20.0;
            s.depth = depth;
            s.offset = offset;
            s.microphone_height = 1.5;
            s.max_frequency = 600.0;
            s.points_per_wavelength = 12.0;
            s.threads = 8;
            const ApertureSpectrum a = solve_canyon(s);
            // The canyon's axis along compass bearing 90 (east); the local
            // frame east/north.
            const double axis = 90.0;
            const auto gains = canyon_wave_gains(a, axis, axes.azimuths, axes.elevations, axes.bands,
                    axes.bands.size(), 5.0);
            const UrbanTable table(axes, {UrbanTable::Receiver {Vec3(0.0, 0.0, s.microphone_height), {}, {}}}, gains,
                    90.0, "wave");
            const auto map = table.match_bands(axes.bands, 0, axes.bands.size());

            for (bool doppler : {true, false})
            {
                std::vector<double> table_err, run_err, weighted_err;
                for (double altitude : {60.0, 150.0, 300.0})
                    for (double heading : {0.0, 45.0, 90.0})
                        for (double lateral : {0.0, 150.0})
                        {
                            // The track: heading from the axis, offset sideways.
                            const double hd = heading * pi / 180.0;
                            const Vec3 along(std::cos(hd), std::sin(hd), 0.0);
                            const Vec3 side(-std::sin(hd), std::cos(hd), 0.0);
                            std::vector<ToneErrors> per;
                            double a_truth = 0.0, a_open_truth = 0.0, a_run = 0.0, a_open_run = 0.0;
                            for (const auto& [rotor, f_h] : harmonics)
                            {
                                (void)rotor;
                                double truth = 0.0, open_truth = 0.0, tab = 0.0, run = 0.0, open_run = 0.0;
                                for (double t = -40.0; t <= 40.0; t += 0.1)
                                {
                                    const Vec3 src = side * lateral + along * (speed * t)
                                            + Vec3(0.0, 0.0, altitude);
                                    const Vec3 mic(0.0, 0.0, s.microphone_height);
                                    const Vec3 r = src - mic;
                                    const double dist = r.norm();
                                    const Vec3 u = r / dist;
                                    const double cos_theta = -along.dot(u);
                                    const double dop = doppler ? 1.0 / (1.0 - mach * cos_theta) : 1.0;
                                    const double f = f_h * dop;
                                    // Compass: x east, y north.
                                    const double azimuth = std::atan2(u.x, u.y) * 180.0 / pi;
                                    const double elevation = std::asin(u.z) * 180.0 / pi;
                                    const double u_c = std::cos(elevation * pi / 180.0)
                                            * std::sin((azimuth - axis) * pi / 180.0);
                                    const double perp = std::hypot(u_c, u.z);
                                    const double phi = std::atan2(u.z, u_c);
                                    const double g_canyon = a.gain_at(f * perp, phi);
                                    const double g_open = a.open_gain(f * perp, phi);
                                    std::size_t band = 0;
                                    while (band + 1 < tone_bands.size()
                                            && f > tone_bands[band] * std::pow(2.0, 1.0 / 6.0))
                                        band++;
                                    double level[14] = {};
                                    REQUIRE(table.lookup(0, azimuth, elevation, map, level));
                                    const double t_band = std::pow(10.0, 0.1 * level[band]);
                                    const double open_avg = open_band(tone_bands[band] * std::pow(2.0, -1.0 / 6.0),
                                            tone_bands[band] * std::pow(2.0, 1.0 / 6.0), s.microphone_height, u.z, c);
                                    // Spherical spreading and reception time.
                                    const double w = 1.0 / (dop * dist * dist);
                                    truth += w * g_canyon;
                                    open_truth += w * g_open;
                                    tab += w * g_open * t_band;
                                    run += w * open_avg * t_band;
                                    open_run += w * open_avg;
                                }
                                const double change = 10.0 * std::log10(truth / open_truth);
                                table_err.push_back(std::abs(10.0 * std::log10(tab / open_truth) - change));
                                run_err.push_back(std::abs(10.0 * std::log10(run / open_run) - change));
                                // Equal source levels, A-weighted at the harmonic.
                                const double aw = std::pow(10.0, 0.1 * a_weighting(f_h));
                                a_truth += aw * truth;
                                a_open_truth += aw * open_truth;
                                a_run += aw * run;
                                a_open_run += aw * open_run;
                            }
                            weighted_err.push_back(std::abs(10.0 * std::log10(a_run / a_open_run)
                                    - 10.0 * std::log10(a_truth / a_open_truth)));
                        }
                auto stats = [](std::vector<double> v) {
                    std::sort(v.begin(), v.end());
                    return std::array<double, 3> {v[v.size() / 2], v[(v.size() * 9) / 10], v.back()};
                };
                const auto te = stats(table_err), re = stats(run_err), we = stats(weighted_err);
                std::printf("%3.1f  %4.0f  %7s  %4zu x %2zu              %6.2f %5.2f %6.2f        %6.2f %5.2f %6.2f       "
                            "%6.2f %6.2f\n",
                        depth / s.width, offset, doppler ? "yes" : "no", harmonics.size(), table_err.size() / harmonics.size(),
                        te[0], te[1], te[2], re[0], re[1], re[2], we[0], we[2]);
            }
        }
}
