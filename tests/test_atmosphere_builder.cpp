// The profile builder (atmosphere_builder.hpp): a surface observation alone is
// --wind's profile to the bit, and measured columns give the profile the rules
// say, with a ground level and no calm level above it.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "atmosphere_builder.hpp"

using namespace niceops;
using Catch::Matchers::ContainsSubstring;

namespace {

struct TempCsv
{
    std::string path;
    explicit TempCsv(const std::string& contents)
    {
        char pattern[] = "/tmp/niceops_builder_XXXXXX";
        const int fd = ::mkstemp(pattern);
        REQUIRE(fd != -1);
        ::close(fd);
        path = pattern;
        std::ofstream out(path);
        out << contents;
    }
    ~TempCsv()
    {
        std::remove(path.c_str());
    }
};

WindSeries gates(const std::vector<double>& z, const std::vector<double>& east,
        const std::vector<double>& north)
{
    return WindSeries {"gates.csv", z, east, north};
}

ColumnSeries series(const std::vector<double>& z, const std::vector<double>& value)
{
    return ColumnSeries {"column.csv", z, value};
}

BuildRules rules_with(double z0, double top)
{
    BuildRules rules;
    rules.roughness_length = z0;
    rules.top = top;
    return rules;
}

} // namespace

TEST_CASE("one surface observation builds --wind's profile to the bit", "[atmosphere_builder]")
{
    int compared = 0, compared_classes = 0;
    for (const Stability stability : {Stability::neutral, Stability::stable, Stability::unstable})
    {
        for (const double z0 : {0.01, 0.03, 0.1, 1.0})
        {
            // An anemometer below ten roughness lengths, and one above.
            for (const double height : {2.0 * z0, 10.0})
            {
                for (const double highest : {150.0, 900.0})
                {
                    SurfaceWind wind;
                    wind.speed = 6.17;
                    wind.direction = 250.0;
                    wind.height = height;
                    wind.roughness_length = z0;
                    wind.surface_layer_depth = 300.0;
                    wind.stability = stability;
                    wind.levels_per_decade = 13;
                    // Power laws from an anemometer near the ground carry the
                    // wind past anything plausible; --wind refuses those too.
                    if (!surface_wind_problem(wind, meters).empty())
                        continue;
                    Atmosphere air;
                    air.temperature = 291.4;
                    air.relative_humidity = 63.0;
                    air.pressure = 99.2;
                    const auto expected = synthetic_profile(wind, air, highest);

                    SurfaceAir surface;
                    surface.temperature = air.temperature;
                    surface.humidity = air.relative_humidity;
                    surface.pressure = air.pressure;
                    surface.wind = wind;
                    BuildRules rules = rules_with(z0, std::max(300.0, 1.05 * highest));
                    rules.surface_layer_depth = 300.0;
                    rules.stability = stability;
                    rules.levels_per_decade = 13;
                    rules.above = BuildRules::Above::lapse;
                    rules.humidity_aloft = BuildRules::HumidityAloft::held_concentration;
                    BuildReport report;
                    const auto built = build_profile({}, surface, rules, meters, report);
                    CAPTURE(describe(stability), z0, height, highest);
                    CHECK(built.height == expected.height);
                    CHECK(built.temperature == expected.temperature);
                    CHECK(built.relative_humidity == expected.relative_humidity);
                    CHECK(built.pressure == expected.pressure);
                    CHECK(built.wind_east == expected.wind_east);
                    CHECK(built.wind_north == expected.wind_north);
                    compared++;
                    compared_classes += stability != Stability::neutral;
                }
            }
        }
    }
    CHECK(compared >= 30);
    CHECK(compared_classes >= 14);
}

TEST_CASE("the builder fills below the lowest gate with the --wind law through it",
        "[atmosphere_builder]")
{
    MeasuredColumns measured;
    // From 270 at the 14 m gate (toward the east), veering above.
    measured.wind = gates({14.0, 38.0, 60.0}, {5.0, 6.0, 6.5}, {0.0, -1.0, -2.0});
    measured.temperature = series({0.0, 100.0}, {290.0, 289.0});
    SurfaceAir surface;
    surface.humidity = 50.0;
    BuildRules rules = rules_with(0.03, 500.0);
    BuildReport report;
    const auto profile = build_profile(measured, surface, rules, meters, report);

    REQUIRE(profile.height.front() == 0.0);
    CHECK(profile.wind_east.front() == 0.0);
    CHECK(profile.height[1] == Catch::Approx(0.3)); // ten roughness lengths
    SurfaceWind law;
    law.speed = 5.0;
    law.height = 14.0;
    law.roughness_length = 0.03;
    for (std::size_t i = 1; i < profile.height.size(); i++)
    {
        const double z = profile.height[i];
        CHECK(std::hypot(profile.wind_east[i], profile.wind_north[i]) > 0.0);
        if (z < 14.0)
        {
            CHECK(profile.wind_east[i] == Catch::Approx(law.speed_at(z)).epsilon(1e-14));
            CHECK(profile.wind_north[i] == 0.0);
        }
        else if (z >= 60.0)
        {
            // Held above the highest gate.
            CHECK(profile.wind_east[i] == 6.5);
            CHECK(profile.wind_north[i] == -2.0);
        }
    }
    // Every gate is a level, its wind exact; between gates u and v are linear.
    for (const double gate : {14.0, 38.0, 60.0})
        CHECK(std::find(profile.height.begin(), profile.height.end(), gate) != profile.height.end());
    CHECK(profile.wind_at(26.0)[0] == Catch::Approx(5.5));
    CHECK(profile.wind_at(26.0)[1] == Catch::Approx(-0.5));
    // No pressure column: the run integrates its own.
    CHECK_FALSE(profile.has_pressure());
    REQUIRE(report.filled_wind.size() == 2);
    CHECK(report.filled_wind[0] == std::pair<double, double> {0.0, 14.0});
    CHECK(report.filled_wind[1] == std::pair<double, double> {60.0, 500.0});
    // Temperature held above its highest measurement, linear below it.
    CHECK(profile.temperature.front() == 290.0); // at() needs a pressure column, and this profile has none
    CHECK(profile.temperature.back() == 289.0);

    SECTION("hold keeps the gate's wind down to the first level")
    {
        rules.fill = BuildRules::Fill::hold;
        const auto held = build_profile(measured, surface, rules, meters, report);
        CHECK(held.wind_east[1] == 5.0);
        CHECK(held.wind_east[0] == 0.0);
    }
    SECTION("none is linear from calm at the ground to the gate")
    {
        rules.fill = BuildRules::Fill::none;
        const auto linear = build_profile(measured, surface, rules, meters, report);
        CHECK(linear.height[1] == 14.0);
        CHECK(linear.wind_at(7.0)[0] == Catch::Approx(2.5));
    }
}

TEST_CASE("no level the builder writes above the ground is calm", "[atmosphere_builder]")
{
    std::mt19937_64 random(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 200; trial++)
    {
        const double z0 = std::pow(10.0, -3.0 + 3.0 * unit(random));
        std::vector<double> z, east, north;
        double level = z0 * (1.5 + 30.0 * unit(random));
        for (int k = 0; k < 1 + trial % 6; k++)
        {
            z.push_back(level);
            east.push_back(-8.0 + 16.0 * unit(random));
            north.push_back(0.1 + 8.0 * unit(random));
            level *= 1.5 + 3.0 * unit(random);
        }
        MeasuredColumns measured;
        measured.wind = gates(z, east, north);
        SurfaceAir surface;
        surface.temperature = 288.0;
        surface.humidity = 40.0;
        BuildRules rules = rules_with(z0, 2.0 * level);
        if (trial % 3 == 0)
            rules.levels_per_decade = 4;
        BuildReport report;
        const auto profile = build_profile(measured, surface, rules, meters, report);
        CAPTURE(trial, z0, z.front());
        for (std::size_t i = 1; i < profile.height.size(); i++)
            CHECK(std::hypot(profile.wind_east[i], profile.wind_north[i]) > 0.0);
    }
}

TEST_CASE("a measured calm gate is refused, or dropped; a calm column is still air",
        "[atmosphere_builder]")
{
    MeasuredColumns measured;
    measured.wind = gates({14.0, 38.0, 60.0}, {5.0, 0.0, 6.5}, {0.0, 0.0, -2.0});
    SurfaceAir surface;
    surface.temperature = 288.0;
    surface.humidity = 40.0;
    BuildRules rules = rules_with(0.03, 500.0);
    BuildReport report;
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, feet, report),
            ContainsSubstring("gives a calm wind at 124.672 ft"));
    rules.drop_calm = true;
    const auto dropped = build_profile(measured, surface, rules, meters, report);
    CHECK(report.calm_dropped == 1);
    CHECK(std::find(dropped.height.begin(), dropped.height.end(), 38.0) == dropped.height.end());

    measured.wind = gates({14.0, 38.0}, {0.0, 0.0}, {0.0, 0.0});
    BuildReport calm;
    const auto still = build_profile(measured, surface, rules, meters, calm);
    CHECK_FALSE(still.has_wind());
    CHECK(calm.warnings.size() == 1);

    // A level the log fill leaves calm, at or below the roughness length.
    measured.wind = gates({14.0}, {5.0}, {0.0});
    rules.levels = {0.02, 1.0};
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("would be calm"));
}

TEST_CASE("the ground level comes from the station, the measurements across it, or the lowest held",
        "[atmosphere_builder]")
{
    MeasuredColumns measured;
    // A balloon standing 10 m below the microphones.
    measured.temperature = series({-10.0, 20.0, 50.0}, {291.0, 288.0, 287.0});
    measured.humidity = series({20.0, 50.0}, {60.0, 55.0});
    measured.wind = gates({-10.0, 20.0, 50.0}, {1.0, 4.0, 5.0}, {1.0, 0.0, 0.0});
    BuildRules rules = rules_with(0.05, 200.0);
    BuildReport report;
    const auto profile = build_profile(measured, {}, rules, meters, report);
    CHECK(profile.height.front() == 0.0);
    CHECK(profile.temperature.front() == Catch::Approx(290.0));
    // Humidity has nothing at or below the ground: its lowest is held down.
    CHECK(profile.relative_humidity.front() == 60.0);
    CHECK(report.warnings.size() == 1);
    CHECK_THAT(report.warnings.front(), ContainsSubstring("No humidity at or below the ground"));
    // The levels below the ground set the ground's value and are gone.
    CHECK(report.below_ground_dropped == 2);
    CHECK(profile.wind_east.front() == 0.0);

    SurfaceAir station;
    station.temperature = 295.0;
    BuildReport with_station;
    const auto from_station = build_profile(measured, station, rules, meters, with_station);
    CHECK(from_station.temperature.front() == 295.0);
    CHECK(with_station.below_ground_dropped == 2);
}

TEST_CASE("an explicit level list is the levels, every measurement, the ground and the top",
        "[atmosphere_builder]")
{
    MeasuredColumns measured;
    measured.temperature = series({0.0, 3.0, 6.0}, {290.0, 289.0, 288.5});
    measured.wind = gates({14.0, 38.0}, {5.0, 6.0}, {0.0, 0.0});
    SurfaceAir surface;
    surface.humidity = 50.0;
    BuildRules rules = rules_with(0.03, 1500.0);
    rules.levels = {0.0, 0.3, 1.0, 3.0, 2000.0};
    BuildReport report;
    const auto profile = build_profile(measured, surface, rules, meters, report);
    CHECK(profile.height == std::vector<double> {0.0, 0.3, 1.0, 3.0, 6.0, 14.0, 38.0, 1500.0});
}

TEST_CASE("the builder refuses what it cannot build", "[atmosphere_builder]")
{
    MeasuredColumns measured;
    measured.wind = gates({14.0, 38.0}, {5.0, 6.0}, {0.0, 0.0});
    SurfaceAir surface;
    surface.temperature = 288.0;
    surface.humidity = 40.0;
    BuildReport report;

    BuildRules rules = rules_with(0.03, 500.0);
    rules.surface_layer_depth = 300.0;
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("would kink at the highest"));

    rules = rules_with(20.0 * 0.3048 / 0.3048, 500.0);
    rules.roughness_length = 2.9;
    measured.wind = gates({2.0, 38.0}, {5.0, 6.0}, {0.0, 0.0});
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("has to stand above the roughness length"));

    rules = rules_with(0.0, 500.0);
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("roughness length of 0 m is not plausible"));

    rules = rules_with(0.03, 500.0);
    rules.stability = Stability::stable;
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("--stability sets the temperature gradient"));

    rules = rules_with(0.03, 500.0);
    MeasuredColumns rising = measured;
    rising.pressure = series({0.0, 100.0, 200.0}, {100.0, 99.0, 99.5});
    CHECK_THROWS_WITH(build_profile(rising, surface, rules, meters, report),
            ContainsSubstring("rises with height"));

    MeasuredColumns wet = measured;
    wet.humidity = series({0.0, 100.0}, {100.0, 104.0});
    SurfaceAir dry_station;
    dry_station.temperature = 288.0;
    CHECK_THROWS_WITH(build_profile(wet, dry_station, rules, meters, report),
            ContainsSubstring("above saturation"));
    rules.clip_humidity = true;
    BuildReport clipped;
    const auto profile = build_profile(wet, dry_station, rules, meters, clipped);
    CHECK(clipped.humidity_clipped == 1);
    CHECK(profile.relative_humidity.back() == 100.0);

    rules = rules_with(0.03, 500.0);
    rules.humidity_aloft = BuildRules::HumidityAloft::held_concentration;
    CHECK_THROWS_WITH(build_profile(measured, surface, rules, meters, report),
            ContainsSubstring("needs a pressure"));

    CHECK_THROWS_WITH(build_profile(measured, {}, rules_with(0.03, 500.0), meters, report),
            ContainsSubstring("no temperature"));
}

TEST_CASE("measured pressure is log-linear between levels and hydrostatic beyond them",
        "[atmosphere_builder]")
{
    MeasuredColumns measured;
    measured.temperature = series({0.0}, {288.15});
    measured.humidity = series({0.0}, {50.0});
    measured.pressure = series({100.0, 300.0}, {100.0, 97.7});
    BuildRules rules = rules_with(0.03, 1000.0);
    rules.levels_per_decade = 8;
    BuildReport report;
    const auto profile = build_profile(measured, {}, rules, meters, report);
    REQUIRE(profile.has_pressure());
    const auto at = [&](double z) {
        const auto it = std::find(profile.height.begin(), profile.height.end(), z);
        REQUIRE(it != profile.height.end());
        return profile.pressure[static_cast<std::size_t>(it - profile.height.begin())];
    };
    CHECK(at(100.0) == 100.0);
    CHECK(at(300.0) == 97.7);
    // Isothermal air: p falls as exp(-g z / R T) outside the measurements.
    const double scale = standard_gravity / (287.05 * 288.15);
    CHECK(at(0.0) == Catch::Approx(100.0 * std::exp(100.0 * scale)).epsilon(1e-12));
    CHECK(at(1000.0) == Catch::Approx(97.7 * std::exp(-700.0 * scale)).epsilon(1e-12));
}

TEST_CASE("read_column_source reads the run's names, the humidity's other forms and fill values",
        "[atmosphere_builder]")
{
    SECTION("the wind convention is the run's")
    {
        const TempCsv file("z_m,T_C,RH,wind_speed_kt,wind_direction\n10,15,50,10,270\n"
                           "100,14,50,12,300\n");
        const auto source = read_column_source(file.path, meters, 0.0, {}, false);
        REQUIRE(source.wind);
        CHECK(source.wind->east[0] > 0.0);
        const auto run = read_atmosphere_profile(file.path, meters, 0.0);
        CHECK(source.wind->east == run.wind_east);
        CHECK(source.wind->north == run.wind_north);
        CHECK(source.temperature->value == run.temperature);
    }
    SECTION("dew point and specific humidity become relative humidity through ISO 9613-1")
    {
        const TempCsv dew("z_ft,T_C,Td_C\n0,20,10\n");
        const auto from_dew = read_column_source(dew.path, feet, 0.0, {}, false);
        auto psat = [](double t) { return 101.325 * std::pow(10.0, -6.8346 * std::pow(273.16 / t, 1.261) + 4.6151); };
        CHECK(from_dew.humidity->value[0]
                == Catch::Approx(100.0 * psat(283.15) / psat(293.15)).epsilon(1e-14));
        const TempCsv specific("z_ft,T_K,p_hPa,q_kgkg\n0,293.15,1000,0.008\n");
        const auto from_q = read_column_source(specific.path, feet, 0.0, {}, false);
        const double e = 0.008 * 100.0 / (0.622 + 0.378 * 0.008);
        CHECK(from_q.humidity->value[0] == Catch::Approx(100.0 * e / psat(293.15)).epsilon(1e-14));
        CHECK(from_q.pressure->value[0] == Catch::Approx(100.0));
    }
    SECTION("fill values are dropped per variable and counted; anything else implausible refused")
    {
        const TempCsv file("z_ft,T_C,wind_speed_kt,wind_direction\n500,9999,10,270\n46,15,99999,99999\n"
                           "125,14.5,12,280\n");
        const auto source = read_column_source(file.path, feet, 0.0, {9999.0, 99999.0}, false);
        CHECK(source.sentinels_dropped == 2);
        CHECK(source.temperature->z.size() == 2);
        CHECK(source.wind->z.size() == 2);
        // Rows come in height order.
        CHECK(source.temperature->z.front() == Catch::Approx(46.0 * 0.3048));
        CHECK_THROWS_WITH(read_column_source(file.path, feet, 0.0, {}, false),
                ContainsSubstring("data row 1 gives a temperature of 9999 C"));
        const auto dropped = read_column_source(file.path, feet, 0.0, {}, true);
        CHECK(dropped.implausible_dropped == 2);
    }
    SECTION("the reader refuses what --atmosphere refuses")
    {
        const TempCsv heightless("T_K,RH\n288,50\n");
        CHECK_THROWS_WITH(read_column_source(heightless.path, feet, 0.0, {}, false),
                ContainsSubstring("names no height column"));
        const TempCsv twice("z_ft,T_K\n0,288\n0,289\n");
        CHECK_THROWS_WITH(read_column_source(twice.path, feet, 0.0, {}, false),
                ContainsSubstring("gives two temperature levels at 0 ft"));
        const TempCsv text("z_ft,T_K\n0,288\nend\n");
        CHECK_THROWS_WITH(read_column_source(text.path, feet, 0.0, {}, false),
                ContainsSubstring("line 3 is not a row of numbers"));
        const TempCsv both("z_ft,T_C,RH,Td_C\n0,20,50,10\n");
        CHECK_THROWS_WITH(read_column_source(both.path, feet, 0.0, {}, false),
                ContainsSubstring("more than one way"));
    }
}

TEST_CASE("a built profile reads back through --atmosphere as itself", "[atmosphere_builder]")
{
    MeasuredColumns measured;
    measured.temperature = series({0.0, 30.0, 150.0}, {290.0, 291.5, 289.0});
    measured.wind = gates({14.0, 38.0, 60.0}, {5.0, 6.0, 6.5}, {0.5, -1.0, -2.0});
    SurfaceAir surface;
    surface.humidity = 70.0;
    surface.pressure = 100.5;
    BuildRules rules = rules_with(0.1 * 0.3048, 1524.0);
    rules.humidity_aloft = BuildRules::HumidityAloft::held_concentration;
    BuildReport report;
    const auto built = build_profile(measured, surface, rules, feet, report);
    char pattern[] = "/tmp/niceops_built_XXXXXX";
    const int fd = ::mkstemp(pattern);
    REQUIRE(fd != -1);
    ::close(fd);
    write_atmosphere_profile(built, pattern, feet);
    const auto read = read_atmosphere_profile(pattern, feet, 0.0);
    std::remove(pattern);
    REQUIRE(read.height.size() == built.height.size());
    for (std::size_t i = 0; i < built.height.size(); i++)
        CHECK(read.height[i] == Catch::Approx(built.height[i]).epsilon(1e-15).margin(1e-300));
    CHECK(read.temperature == built.temperature);
    CHECK(read.relative_humidity == built.relative_humidity);
    CHECK(read.pressure == built.pressure);
    CHECK(read.wind_east == built.wind_east);
    CHECK(read.wind_north == built.wind_north);
}

namespace {

// A run whose gates follow the log law exactly, from `direction`, with a
// neutral temperature column (theta constant) from the ground to 500 ft.
ColumnSource log_law_run(double z0, double friction, double direction_deg,
        const std::vector<double>& gates_ft = {46.0, 125.0, 197.0, 348.0})
{
    ColumnSource run;
    run.path = "run.csv";
    WindSeries wind;
    wind.source = run.path;
    const double from = direction_deg * pi / 180.0;
    for (const double ft : gates_ft)
    {
        const double z = ft * ft2m;
        const double speed = friction / 0.4 * std::log(z / z0);
        wind.z.push_back(z);
        wind.east.push_back(-speed * std::sin(from));
        wind.north.push_back(-speed * std::cos(from));
    }
    run.wind = wind;
    ColumnSeries temperature;
    temperature.source = run.path;
    for (double ft = 0.0; ft <= 500.0; ft += 10.0)
    {
        temperature.z.push_back(ft * ft2m);
        temperature.value.push_back(290.0 + dry_adiabatic_lapse * ft * ft2m);
    }
    run.temperature = temperature;
    return run;
}

// 1% Gaussian speed errors, the same on every standard library
// (std::normal_distribution is not): Box-Muller on the engine's own bits.
std::vector<ColumnSource> noisy_runs(std::uint64_t seed, int count)
{
    std::mt19937_64 noise(seed);
    auto uniform = [&noise] { return (static_cast<double>(noise() >> 11) + 0.5) * 0x1p-53; };
    std::vector<ColumnSource> runs;
    for (int j = 0; j < count; j++)
    {
        auto run = log_law_run(0.03, 0.2 + 0.3 * j / count, 230.0 + 40.0 * j / count);
        for (std::size_t i = 0; i < run.wind->z.size(); i++)
        {
            const double gaussian = std::sqrt(-2.0 * std::log(uniform())) * std::cos(2.0 * pi * uniform());
            const double factor = 1.0 + 0.01 * gaussian;
            run.wind->east[i] *= factor;
            run.wind->north[i] *= factor;
        }
        runs.push_back(run);
    }
    return runs;
}

std::vector<std::string> labels_for(std::size_t n)
{
    std::vector<std::string> labels;
    for (std::size_t j = 0; j < n; j++)
        labels.push_back("run" + std::to_string(j));
    return labels;
}

} // namespace

TEST_CASE("fit_roughness_length recovers the roughness length of exact log-law runs",
        "[atmosphere_builder][fit]")
{
    for (const double z0 : {0.01, 0.03, 0.1})
    {
        std::vector<ColumnSource> runs;
        for (int j = 0; j < 6; j++)
            runs.push_back(log_law_run(z0, 0.25 + 0.07 * j, 200.0 + 15.0 * j));
        FitRules rules;
        rules.draws = 200;
        const auto fit = fit_roughness_length(runs, labels_for(runs.size()), rules);
        CAPTURE(z0);
        CHECK(fit.z0 == Catch::Approx(z0).epsilon(1e-9));
        CHECK(fit.used.size() == 6);
        CHECK(fit.excluded.empty());
        CHECK_FALSE(fit.at_bound);
        CHECK(fit.residual_rms < 1e-12);
        CHECK(fit.median_run_z0 == Catch::Approx(z0).epsilon(1e-9));
        CHECK(fit.z0_low == Catch::Approx(z0).epsilon(1e-9));
        CHECK(fit.z0_high == Catch::Approx(z0).epsilon(1e-9));

        // Scaling every speed leaves z0 alone, and so does the gates' order.
        std::vector<ColumnSource> scaled = runs, reversed = runs;
        for (auto& run : scaled)
        {
            for (auto& v : run.wind->east)
                v *= 1.7;
            for (auto& v : run.wind->north)
                v *= 1.7;
        }
        for (auto& run : reversed)
        {
            std::reverse(run.wind->z.begin(), run.wind->z.end());
            std::reverse(run.wind->east.begin(), run.wind->east.end());
            std::reverse(run.wind->north.begin(), run.wind->north.end());
        }
        CHECK(fit_roughness_length(scaled, labels_for(6), rules).z0 == Catch::Approx(z0).epsilon(1e-9));
        CHECK(fit_roughness_length(reversed, labels_for(6), rules).z0 == fit.z0);
    }
}

TEST_CASE("fit_roughness_length screens out stable, veering and light-wind runs",
        "[atmosphere_builder][fit]")
{
    std::vector<ColumnSource> runs;
    for (int j = 0; j < 5; j++)
        runs.push_back(log_law_run(0.03, 0.3, 250.0));
    // Stable: theta rising 2 K over the gates.
    ColumnSource stable = log_law_run(0.03, 0.3, 250.0);
    for (std::size_t i = 0; i < stable.temperature->z.size(); i++)
        stable.temperature->value[i] += 0.04 * stable.temperature->z[i];
    // Veering: the top gate turned 40 degrees.
    ColumnSource veering = log_law_run(0.03, 0.3, 250.0);
    {
        auto& wind = *veering.wind;
        const double turn = 40.0 * pi / 180.0;
        const double e = wind.east[2], n = wind.north[2];
        wind.east[2] = e * std::cos(turn) - n * std::sin(turn);
        wind.north[2] = e * std::sin(turn) + n * std::cos(turn);
    }
    ColumnSource light = log_law_run(0.03, 0.05, 250.0);
    ColumnSource unscreened = log_law_run(0.03, 0.3, 250.0);
    unscreened.temperature.reset();
    for (const auto* run : {&stable, &veering, &light, &unscreened})
        runs.push_back(*run);
    FitRules rules;
    rules.draws = 100;
    const auto fit = fit_roughness_length(runs, labels_for(runs.size()), rules);
    REQUIRE(fit.excluded.size() == 4);
    CHECK_THAT(fit.excluded[0].second, ContainsSubstring("not neutral"));
    CHECK_THAT(fit.excluded[1].second, ContainsSubstring("turns"));
    CHECK_THAT(fit.excluded[2].second, ContainsSubstring("light wind"));
    CHECK_THAT(fit.excluded[3].second, ContainsSubstring("no temperature"));
    CHECK(fit.excluded[0].first == "run5");
    rules.assume_neutral = true;
    CHECK(fit_roughness_length(runs, labels_for(runs.size()), rules).used.size() == 6);
    rules.min_runs = 7;
    CHECK_THROWS_WITH(fit_roughness_length(runs, labels_for(runs.size()), rules),
            ContainsSubstring("fewer than the 7 a fit needs"));
}

TEST_CASE("the bootstrap's index draw is the same on every standard library",
        "[atmosphere_builder][fit]")
{
    // The engine's sequence is the standard's ([rand.predef]: the 10000th
    // value of a default-constructed mt19937_64).
    std::mt19937_64 standard;
    standard.discard(9999);
    CHECK(standard() == 9981545732273789042ull);

    // Rejection: for n = 3, 2^64 - 1 is the one value thrown back.
    struct Scripted
    {
        using result_type = std::uint64_t;
        static constexpr result_type min()
        {
            return 0;
        }
        static constexpr result_type max()
        {
            return std::numeric_limits<std::uint64_t>::max();
        }
        std::vector<result_type> values;
        std::size_t next = 0;
        result_type operator()()
        {
            return values[next++];
        }
    };
    Scripted scripted {{std::numeric_limits<std::uint64_t>::max(), 5}};
    CHECK(portable_index(scripted, 3) == 2);
    CHECK(scripted.next == 2);
    Scripted power {{std::numeric_limits<std::uint64_t>::max()}};
    CHECK(portable_index(power, 8) == 7);

    // A known sequence, pinned so that a build on another library that drew
    // differently would fail here rather than move an interval.
    std::mt19937_64 engine(1);
    std::vector<std::size_t> drawn;
    for (int k = 0; k < 10; k++)
        drawn.push_back(portable_index(engine, 10));
    CHECK(drawn == std::vector<std::size_t> {8, 2, 0, 6, 4, 9, 8, 5, 8, 4});
}

TEST_CASE("the bootstrap interval covers the roughness length of noisy runs",
        "[atmosphere_builder][fit]")
{
    const auto runs = noisy_runs(1, 20);
    FitRules rules;
    rules.draws = 1000;
    const auto fit = fit_roughness_length(runs, labels_for(runs.size()), rules);
    CHECK(fit.z0_low < fit.z0);
    CHECK(fit.z0 < fit.z0_high);
    CHECK(fit.z0_low <= 0.03);
    CHECK(0.03 <= fit.z0_high);
}

TEST_CASE("the bootstrap interval covers the truth in at least 90 of 100 seeds",
        "[.slow][atmosphere_builder][fit]")
{
    int hits = 0;
    for (std::uint64_t seed = 1; seed <= 100; seed++)
    {
        const auto runs = noisy_runs(seed, 20);
        FitRules rules;
        rules.seed = seed;
        const auto fit = fit_roughness_length(runs, labels_for(runs.size()), rules);
        hits += fit.z0_low <= 0.03 && 0.03 <= fit.z0_high;
    }
    CHECK(hits >= 90);
}
