// The stratified atmosphere (--atmosphere) and the layered absorption through
// it.
//
// The absorption is checked three independent ways: a uniform profile has to
// give back the homogeneous coefficients; a path through a lapsing profile has
// to match a brute-force quadrature of ISO 9613-1 along it; and it has to agree
// with the layered procedure of 14 CFR 36 Appendix A, which is what noise
// certification uses.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "atmosphere_profile.hpp"
#include "util/csv.hpp"

using namespace niceops;

namespace {

struct TempCsv
{
    std::string path;
    explicit TempCsv(const std::string& contents)
    {
        char pattern[] = "/tmp/niceops_profile_XXXXXX";
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

// ISA to 3 km: temperature lapsing at 6.5 K/km, pressure from the barometric
// formula, and a humidity that falls with height so every input varies.
AtmosphereProfile lapsing_profile()
{
    AtmosphereProfile profile;
    profile.source = "test";
    for (double z = 0.0; z <= 3000.0; z += 250.0)
    {
        const double t = 288.15 - 0.0065 * z;
        profile.height.push_back(z);
        profile.temperature.push_back(t);
        profile.relative_humidity.push_back(70.0 - 0.01 * z);
        profile.pressure.push_back(101.325 * std::pow(t / 288.15, standard_gravity / (0.0065 * 287.05)));
    }
    return profile;
}

const std::vector<double> third_octaves = {50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0,
        8000.0, 10000.0};

} // namespace

TEST_CASE("read_atmosphere_profile reads named units and orders by height", "[atmosphere_profile]")
{
    // A sounding top down, in Celsius and hPa, with a wind column to ignore.
    TempCsv file("z_ft, T_C, RH, p_hPa, wind_kt\n"
                 "1000, 13.0, 50, 977.2, 12\n"
                 "0, 15.0, 60, 1013.25, 5\n");
    const auto profile = read_atmosphere_profile(file.path, feet, 0.0);
    REQUIRE(profile.height.size() == 2);
    CHECK(profile.height[0] == 0.0);
    CHECK(profile.height[1] == Catch::Approx(304.8));
    CHECK(profile.temperature[0] == Catch::Approx(288.15));
    CHECK(profile.temperature[1] == Catch::Approx(286.15));
    CHECK(profile.relative_humidity[0] == 60.0);
    REQUIRE(profile.has_pressure());
    CHECK(profile.pressure[0] == Catch::Approx(101.325));
    CHECK(profile.pressure[1] == Catch::Approx(97.72));
}

TEST_CASE("read_atmosphere_profile tells bare units from the values", "[atmosphere_profile]")
{
    TempCsv kelvin_kpa("z, T, RH, p\n0, 288.15, 50, 101.325\n");
    const auto a = read_atmosphere_profile(kelvin_kpa.path, meters, 0.0);
    CHECK(a.temperature[0] == 288.15);
    CHECK(a.pressure[0] == 101.325);

    TempCsv celsius_pa("z, temperature, humidity, pressure\n0, 15, 50, 101325\n");
    const auto b = read_atmosphere_profile(celsius_pa.path, meters, 0.0);
    CHECK(b.temperature[0] == Catch::Approx(288.15));
    CHECK(b.pressure[0] == Catch::Approx(101.325));

    // "T (K)": the parenthesized unit is dropped by the header matching, so
    // this is the bare name too, and the values decide.
    TempCsv parenthesized("z (m), T (K), RH (%)\n0, 280, 20\n100, 279, 20\n");
    const auto c = read_atmosphere_profile(parenthesized.path, meters, 0.0);
    CHECK(c.temperature[1] == 279.0);
    CHECK_FALSE(c.has_pressure());
}

TEST_CASE("read_atmosphere_profile converts altitude to height above the observers",
        "[atmosphere_profile]")
{
    TempCsv file("alt_ft, T_K, RH\n1000, 285, 40\n2000, 283, 40\n");
    const double ground_elevation = 500.0 * ft2m;
    const auto profile = read_atmosphere_profile(file.path, feet, ground_elevation);
    CHECK(profile.height[0] == Catch::Approx(500.0 * ft2m));
    CHECK(profile.height[1] == Catch::Approx(1500.0 * ft2m));
}

TEST_CASE("read_atmosphere_profile refuses what it cannot read safely", "[atmosphere_profile]")
{
    // A unit in the height column that disagrees with --units.
    TempCsv meters_file("z_m, T_K, RH\n0, 288, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(meters_file.path, feet, 0.0), std::runtime_error);
    CHECK_NOTHROW(read_atmosphere_profile(meters_file.path, meters, 0.0));

    // Humidity given as a fraction.
    TempCsv fraction("z, T_K, RH\n0, 288, 0.5\n100, 287, 0.4\n");
    CHECK_THROWS_AS(read_atmosphere_profile(fraction.path, meters, 0.0), std::runtime_error);

    // Two levels at one height.
    TempCsv duplicate("z, T_K, RH\n0, 288, 50\n0, 287, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(duplicate.path, meters, 0.0), std::runtime_error);

    // Kelvin labeled as Celsius.
    TempCsv mislabeled("z, T_C, RH\n0, 288, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(mislabeled.path, meters, 0.0), std::runtime_error);

    // No humidity at all.
    TempCsv no_humidity("z, T_K\n0, 288\n");
    CHECK_THROWS_AS(read_atmosphere_profile(no_humidity.path, meters, 0.0), std::runtime_error);

    // Pressure rising with height.
    TempCsv rising("z, T_K, RH, p_kPa\n0, 288, 50, 90\n100, 287, 50, 95\n");
    CHECK_THROWS_AS(read_atmosphere_profile(rising.path, meters, 0.0), std::runtime_error);
}

TEST_CASE("AtmosphereProfile::at interpolates and holds the ends", "[atmosphere_profile]")
{
    AtmosphereProfile profile;
    profile.height = {10.0, 110.0};
    profile.temperature = {290.0, 280.0};
    profile.relative_humidity = {60.0, 40.0};
    profile.pressure = {100.0, 90.0};

    const auto middle = profile.at(60.0);
    CHECK(middle.temperature == Catch::Approx(285.0));
    CHECK(middle.relative_humidity == Catch::Approx(50.0));
    // Linear in the logarithm: the geometric mean at the midpoint.
    CHECK(middle.pressure == Catch::Approx(std::sqrt(100.0 * 90.0)));

    CHECK(profile.at(0.0).temperature == 290.0);
    CHECK(profile.at(500.0).pressure == 90.0);
}

TEST_CASE("with_hydrostatic_pressure matches the analytic atmospheres", "[atmosphere_profile]")
{
    const double g_over_r = standard_gravity / 287.05;

    SECTION("isothermal: an exponential")
    {
        AtmosphereProfile profile;
        profile.height = {0.0, 1000.0, 2000.0};
        profile.temperature = {280.0, 280.0, 280.0};
        profile.relative_humidity = {50.0, 50.0, 50.0};
        const auto filled = profile.with_hydrostatic_pressure(100.0);
        for (std::size_t i = 0; i < 3; i++)
        {
            CHECK(filled.pressure[i]
                    == Catch::Approx(100.0 * std::exp(-g_over_r * profile.height[i] / 280.0))
                               .epsilon(1e-12));
        }
    }

    SECTION("constant lapse: the barometric formula")
    {
        auto profile = lapsing_profile();
        const auto expected = profile.pressure;
        profile.pressure.clear();
        const auto filled = profile.with_hydrostatic_pressure(101.325);
        for (std::size_t i = 0; i < expected.size(); i++)
            CHECK(filled.pressure[i] == Catch::Approx(expected[i]).epsilon(1e-12));
    }

    SECTION("levels starting above the ground are integrated up to from z = 0")
    {
        AtmosphereProfile profile;
        profile.height = {100.0, 200.0};
        profile.temperature = {280.0, 280.0};
        profile.relative_humidity = {50.0, 50.0};
        const auto filled = profile.with_hydrostatic_pressure(100.0);
        CHECK(filled.pressure[0]
                == Catch::Approx(100.0 * std::exp(-g_over_r * 100.0 / 280.0)).epsilon(1e-12));
    }
}

TEST_CASE("LayeredAbsorption through a uniform profile is the homogeneous absorption",
        "[atmosphere_profile]")
{
    AtmosphereProfile profile;
    profile.height = {0.0};
    profile.temperature = {283.0};
    profile.relative_humidity = {35.0};
    profile.pressure = {95.0};
    const auto homogeneous = absorption_coefficients(profile.at(0.0), third_octaves);

    const LayeredAbsorption layered(profile, third_octaves, 0, third_octaves.size(), 0.0, 800.0);
    std::vector<double> f_source(third_octaves.size()), f_observer(third_octaves.size()),
            mean(third_octaves.size());
    for (const double z_source : {0.0005, 37.3, 150.0, 799.9})
    {
        layered.integral(z_source, f_source.data());
        layered.integral(0.0, f_observer.data());
        layered.mean_coefficients(z_source, f_source.data(), 0.0, f_observer.data(), mean.data());
        for (std::size_t band = 0; band < third_octaves.size(); band++)
            CHECK(mean[band] == Catch::Approx(homogeneous[band]).epsilon(1e-12));
    }
}

TEST_CASE("LayeredAbsorption matches a quadrature of ISO 9613-1 along the path",
        "[atmosphere_profile]")
{
    const auto profile = lapsing_profile();
    const LayeredAbsorption layered(profile, third_octaves, 0, third_octaves.size(), 0.0, 3000.0);

    // Midpoint rule at 1 cm, straight from the standard at every step.
    auto brute_force = [&](double z_a, double z_b, std::size_t band) {
        const int steps = static_cast<int>(std::ceil(std::fabs(z_a - z_b) / 0.01));
        const double h = (z_b - z_a) / steps;
        double sum = 0.0;
        for (int i = 0; i < steps; i++)
            sum += profile.at(z_a + (i + 0.5) * h).attenuation_coefficient(third_octaves[band]);
        return sum * std::fabs(h) / std::fabs(z_a - z_b);
    };

    std::vector<double> f_a(third_octaves.size()), f_b(third_octaves.size()),
            mean(third_octaves.size()), reversed(third_octaves.size());
    for (const auto& [z_a, z_b] : {std::pair {300.0, 0.0}, std::pair {1234.5, 1.2},
                                          std::pair {2999.0, 2500.25}})
    {
        layered.integral(z_a, f_a.data());
        layered.integral(z_b, f_b.data());
        layered.mean_coefficients(z_a, f_a.data(), z_b, f_b.data(), mean.data());
        layered.mean_coefficients(z_b, f_b.data(), z_a, f_a.data(), reversed.data());
        for (std::size_t band = 0; band < third_octaves.size(); band++)
        {
            CHECK(mean[band] == Catch::Approx(brute_force(z_a, z_b, band)).epsilon(1e-6));
            // Which end is the source does not matter.
            CHECK(reversed[band] == Catch::Approx(mean[band]).epsilon(1e-14));
        }
    }
}

TEST_CASE("LayeredAbsorption is continuous as the two heights meet", "[atmosphere_profile]")
{
    const auto profile = lapsing_profile();
    const LayeredAbsorption layered(profile, third_octaves, 0, third_octaves.size(), 0.0, 3000.0);
    std::vector<double> point(third_octaves.size()), f_a(third_octaves.size()),
            f_b(third_octaves.size()), just_outside(third_octaves.size()),
            just_inside(third_octaves.size());
    const double z = 812.37;
    layered.coefficients(z, point.data());
    // Either side of the millimeter where the difference quotient hands over
    // to the midpoint coefficient.
    layered.integral(z + 0.6e-3, f_a.data());
    layered.integral(z - 0.6e-3, f_b.data());
    layered.mean_coefficients(z + 0.6e-3, f_a.data(), z - 0.6e-3, f_b.data(), just_outside.data());
    layered.integral(z + 0.4e-3, f_a.data());
    layered.integral(z - 0.4e-3, f_b.data());
    layered.mean_coefficients(z + 0.4e-3, f_a.data(), z - 0.4e-3, f_b.data(), just_inside.data());
    for (std::size_t band = 0; band < third_octaves.size(); band++)
    {
        CHECK(just_outside[band] == Catch::Approx(point[band]).epsilon(1e-8));
        CHECK(just_inside[band] == Catch::Approx(point[band]).epsilon(1e-8));
    }
}

TEST_CASE("LayeredAbsorption agrees with the 14 CFR 36 Appendix A layered procedure",
        "[atmosphere_profile]")
{
    // A36.7: equal layers of at least 30 m, each taking the mean of the
    // coefficients at its top and bottom.  Both that and the table are
    // trapezoid rules on a smooth integrand, at 30 m and at 1 m, so they
    // should agree closely -- and the table is the finer of the two.
    const auto profile = lapsing_profile();
    const double z_aircraft = 300.0;
    const LayeredAbsorption layered(profile, third_octaves, 0, third_octaves.size(), 0.0,
            z_aircraft);
    std::vector<double> f_top(third_octaves.size()), f_ground(third_octaves.size()),
            mean(third_octaves.size());
    layered.integral(z_aircraft, f_top.data());
    layered.integral(0.0, f_ground.data());
    layered.mean_coefficients(z_aircraft, f_top.data(), 0.0, f_ground.data(), mean.data());

    const int layers = 10; // 30 m each
    for (std::size_t band = 0; band < third_octaves.size(); band++)
    {
        double sum = 0.0;
        for (int layer = 0; layer < layers; layer++)
        {
            const double bottom = z_aircraft * layer / layers;
            const double top = z_aircraft * (layer + 1) / layers;
            sum += 0.5
                    * (profile.at(bottom).attenuation_coefficient(third_octaves[band])
                            + profile.at(top).attenuation_coefficient(third_octaves[band]));
        }
        CHECK(mean[band] == Catch::Approx(sum / layers).epsilon(1e-4));
    }
}

TEST_CASE("CompassToLocal turns east/north into the frame whose +x runs along a compass bearing", "[atmosphere_profile]")
{
    // The propagation table is indexed by compass bearing, and the run turns the local
    // separation back with the inverse of this.  Swapping atan2's arguments in that step
    // passed every test, because the only table built had one bearing.
    for (const double bearing : {0.0, 90.0, 210.0, 333.0})
    {
        const CompassToLocal frame(bearing);
        const double b = bearing * 3.14159265358979323846 / 180.0;
        // A unit vector along the bearing (east = sin, north = cos) is the frame's +x.
        const auto along = frame({std::sin(b), std::cos(b)});
        CHECK(along[0] == Catch::Approx(1.0).margin(1e-12));
        CHECK(along[1] == Catch::Approx(0.0).margin(1e-12));
        // The frame's +y is 90 degrees to its left (counter-clockwise): bearing - 90.
        const auto left = frame({std::sin(b - 0.5 * 3.14159265358979323846), std::cos(b - 0.5 * 3.14159265358979323846)});
        CHECK(left[0] == Catch::Approx(0.0).margin(1e-12));
        CHECK(left[1] == Catch::Approx(1.0).margin(1e-12));
        // And the inverse brings both back, and a compass bearing of a local vector comes out right.
        const auto back = frame.east_north({1.0, 0.0});
        CHECK(back[0] == Catch::Approx(std::sin(b)).margin(1e-12));
        CHECK(back[1] == Catch::Approx(std::cos(b)).margin(1e-12));
        const auto round_trip = frame.east_north(frame({3.0, -4.0}));
        CHECK(round_trip[0] == Catch::Approx(3.0).margin(1e-12));
        CHECK(round_trip[1] == Catch::Approx(-4.0).margin(1e-12));
        const double toward = std::atan2(back[0], back[1]) * 180.0 / 3.14159265358979323846;
        CHECK(std::remainder(toward - bearing, 360.0) == Catch::Approx(0.0).margin(1e-9));
    }
}

// --wind: the profile built from one surface observation (synthetic_profile).

namespace {

Atmosphere surface_air()
{
    Atmosphere air;
    air.temperature = 290.0;
    air.pressure = 95.0;
    air.relative_humidity = 30.0;
    return air;
}

SurfaceWind westerly()
{
    SurfaceWind wind;
    wind.speed = 6.0;
    wind.direction = 270.0;
    wind.height = 2.0;
    return wind;
}

} // namespace

TEST_CASE("the --wind log law passes through the observation and is held above the surface layer",
        "[atmosphere_profile][wind]")
{
    SurfaceWind wind = westerly();
    wind.roughness_length = 0.03;
    CHECK(wind.speed_at(wind.height) == wind.speed);
    CHECK(wind.speed_at(wind.roughness_length) == 0.0);
    CHECK(wind.speed_at(0.0) == 0.0);
    CHECK(wind.speed_at(-1.0) == 0.0);
    // u = u_ref ln(z/z0) / ln(z_ref/z0).
    CHECK(wind.speed_at(20.0)
            == Catch::Approx(6.0 * std::log(20.0 / 0.03) / std::log(2.0 / 0.03)).epsilon(1e-14));
    CHECK(wind.speed_at(500.0) == wind.speed_at(wind.surface_layer_depth));
    CHECK(wind.speed_at(150.0) < wind.speed_at(wind.surface_layer_depth));
    CHECK(wind.temperature_gradient() == dry_adiabatic_lapse);
    CHECK(wind.first_level() == Catch::Approx(0.3));
}

TEST_CASE("the --stability classes are TA Luft's power laws, named either way",
        "[atmosphere_profile][wind]")
{
    struct Expected
    {
        const char* roman;
        const char* pasquill;
        double exponent;
        bool inversion;
    };
    // TA Luft (1986): m = 0.42, 0.35, 0.22, 0.20, 0.15, 0.09 for I to V.
    const Expected classes[] = {{"I", "F", 0.42, true}, {"II", "E", 0.35, true},
            {"III/1", "D", 0.22, false}, {"III/2", "C", 0.20, false}, {"IV", "B", 0.15, false},
            {"V", "A", 0.09, false}};
    for (const auto& expected : classes)
    {
        const Stability stability = parse_stability(expected.roman);
        CHECK(parse_stability(expected.pasquill) == stability);
        CHECK(parse_stability(std::string(" ") + expected.pasquill) == stability);
        SurfaceWind wind = westerly();
        wind.stability = stability;
        CHECK(wind.speed_at(wind.height) == Catch::Approx(wind.speed).epsilon(1e-15));
        CHECK(wind.speed_at(4.0 * wind.height) / wind.speed
                == Catch::Approx(std::pow(4.0, expected.exponent)).epsilon(1e-14));
        CHECK((wind.temperature_gradient() > 0.0) == expected.inversion);
        // Superadiabatic for the unstable classes, never past the 1.9 K/100 m
        // boundary of Pasquill B by more than class A's own range.
        CHECK(wind.temperature_gradient() >= -0.0201);
        CHECK(describe(stability).find(expected.roman) != std::string::npos);
    }
    CHECK(parse_stability("neutral") == Stability::neutral);
    CHECK(parse_stability("ii") == Stability::stable);
    CHECK(parse_stability("f") == Stability::very_stable);
    CHECK_THROWS_AS(parse_stability("G"), std::runtime_error);
    CHECK_THROWS_AS(parse_stability("stable"), std::runtime_error);
    // Every name the option offers parses.
    for (const auto& name : stability_names())
        CHECK_NOTHROW(parse_stability(name));
}

TEST_CASE("parse_surface_wind reads the speed's named unit, the direction and the height in --units",
        "[atmosphere_profile][wind]")
{
    const auto knots = parse_surface_wind("10kt,270,33", feet);
    CHECK(knots.speed == Catch::Approx(10.0 * knot).epsilon(1e-15));
    CHECK(knots.direction == 270.0);
    CHECK(knots.height == Catch::Approx(33.0 * 0.3048).epsilon(1e-15));

    const auto mps = parse_surface_wind(" 5.5 mps , 90 ", meters);
    CHECK(mps.speed == 5.5);
    CHECK(mps.direction == 90.0);
    CHECK(mps.height == 10.0); // the WMO standard exposure when omitted, whatever --units says
    CHECK(parse_surface_wind("5.5m/s,90", feet).height == 10.0);
    CHECK(parse_surface_wind("4KT,0,10", meters).height == 10.0);
    CHECK(parse_surface_wind("4KT,0,10", meters).speed == Catch::Approx(4.0 * knot).epsilon(1e-15));
    CHECK(parse_surface_wind("22.369mph,180", meters).speed == Catch::Approx(10.0).epsilon(1e-4));

    // No unit, or one that is not a speed: refused, never taken from --units.
    CHECK_THROWS_AS(parse_surface_wind("10,270", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10fps,270", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt,270,33,1", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt,west", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt,270deg", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt,270,33ft", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("10kt,270,", feet), std::runtime_error);
    CHECK_THROWS_AS(parse_surface_wind("kt,270", feet), std::runtime_error);
}

TEST_CASE("surface_wind_problem refuses what the law cannot carry", "[atmosphere_profile][wind]")
{
    CHECK(surface_wind_problem(westerly(), feet).empty());
    auto with = [](auto change) {
        SurfaceWind wind = westerly();
        change(wind);
        return surface_wind_problem(wind, meters);
    };
    CHECK_FALSE(with([](SurfaceWind& w) { w.speed = -1.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.speed = 60.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.speed = std::nan(""); }).empty());
    CHECK(with([](SurfaceWind& w) { w.speed = 0.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.direction = 361.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.direction = -10.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.roughness_length = 0.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.roughness_length = 5.0; }).empty());
    // The log law is zero at z0: an anemometer at or below it gives nothing to scale.
    CHECK_FALSE(with([](SurfaceWind& w) { w.roughness_length = 2.0; }).empty());
    CHECK(with([](SurfaceWind& w) {
        w.roughness_length = 2.0;
        w.height = 2.5;
    }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.height = 0.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.surface_layer_depth = 2.0; }).empty());
    CHECK_FALSE(with([](SurfaceWind& w) { w.levels_per_decade = 0; }).empty());
    // An anemometer a hair above z0 carries a plausible surface wind to an
    // implausible one aloft.
    CHECK_FALSE(with([](SurfaceWind& w) {
        w.roughness_length = 1.9;
        w.speed = 20.0;
    }).empty());
}

TEST_CASE("synthetic_profile never puts a calm level above the ground", "[atmosphere_profile][wind]")
{
    // --refraction stratified reads 22-28 dB low with a calm level just above
    // z = 0, so every level above the ground carries wind: every class, every
    // roughness, anemometers above and below ten roughness lengths.
    for (const auto& name : stability_names())
    {
        for (const double z0 : {0.0002, 0.03, 0.5, 2.0})
        {
            for (const double height : {2.0, 10.0, 25.0})
            {
                SurfaceWind wind = westerly();
                wind.stability = parse_stability(name);
                wind.roughness_length = z0;
                wind.height = height;
                if (!surface_wind_problem(wind, meters).empty())
                    continue;
                AtmosphereProfile profile;
                REQUIRE_NOTHROW(profile = synthetic_profile(wind, surface_air(), 450.0));
                REQUIRE(profile.has_wind());
                CHECK(profile.height.front() == 0.0);
                CHECK(profile.wind_east.front() == 0.0);
                for (std::size_t i = 1; i < profile.height.size(); i++)
                {
                    INFO(name << " z0 " << z0 << " at " << profile.height[i] << " m");
                    CHECK(profile.height[i] > profile.height[i - 1]);
                    CHECK(std::hypot(profile.wind_east[i], profile.wind_north[i]) > 0.0);
                }
            }
        }
    }
}

TEST_CASE("synthetic_profile's levels hold the observation and reach past the trajectory",
        "[atmosphere_profile][wind]")
{
    const SurfaceWind wind = westerly();
    const auto profile = synthetic_profile(wind, surface_air(), 450.0);
    const auto& z = profile.height;
    // Ten roughness lengths, then geometric at sixteen a decade.
    CHECK(z[1] == Catch::Approx(0.3));
    for (std::size_t i = 2; i < z.size(); i++)
        CHECK(z[i] / z[i - 1] <= std::pow(10.0, 1.0 / 16.0) * (1.0 + 1e-12));
    CHECK(z.size() == 1 + 1 + 52); // 0, then 0.3 m and 14 + 32 + 6 steps to 2, 200 and 472.5 m
    // The anemometer and the top of the surface layer are levels, the first with
    // the observed wind exactly.
    const auto at = [&](double height) {
        return static_cast<std::size_t>(std::find(z.begin(), z.end(), height) - z.begin());
    };
    REQUIRE(at(wind.height) < z.size());
    REQUIRE(at(wind.surface_layer_depth) < z.size());
    CHECK(profile.wind_east[at(wind.height)] == Catch::Approx(wind.speed).epsilon(1e-15));
    // The top clears the highest point by 5 percent, and never stops short of
    // the surface layer.
    CHECK(profile.top() == Catch::Approx(1.05 * 450.0));
    CHECK(synthetic_profile(wind, surface_air(), 20.0).top() == wind.surface_layer_depth);
    // Linear between levels: the interpolant follows the law to within half a
    // percent at eight levels a decade, and sixteen times closer at 32.
    for (const int levels : {8, 32})
    {
        SurfaceWind fine = wind;
        fine.levels_per_decade = levels;
        const auto refined = synthetic_profile(fine, surface_air(), 450.0);
        double worst = 0.0;
        for (double height = 0.31; height < 470.0; height *= 1.01)
        {
            worst = std::max(worst,
                    std::fabs(refined.wind_at(height)[0] / fine.speed_at(height) - 1.0));
        }
        CHECK(worst < (levels == 8 ? 0.005 : 0.0005));
    }
}

TEST_CASE("synthetic_profile's air lapses dry-adiabatically, hydrostatically, keeping its water",
        "[atmosphere_profile][wind]")
{
    const Atmosphere surface = surface_air();
    const SurfaceWind wind = westerly();
    const auto profile = synthetic_profile(wind, surface, 1000.0);
    REQUIRE(profile.has_pressure());
    CHECK(profile.temperature.front() == surface.temperature);
    CHECK(profile.pressure.front() == surface.pressure);
    CHECK(profile.relative_humidity.front() == surface.relative_humidity);
    const double concentration = surface.molar_concentration_water_vapor();
    for (std::size_t i = 0; i < profile.height.size(); i++)
    {
        const double z = profile.height[i];
        CHECK(profile.temperature[i] == Catch::Approx(surface.temperature - 0.0098 * z).epsilon(1e-14));
        // Barometric formula for a linear lapse.
        const double t = surface.temperature - 0.0098 * z;
        CHECK(profile.pressure[i]
                == Catch::Approx(surface.pressure
                        * std::pow(t / surface.temperature, standard_gravity / (0.0098 * 287.05)))
                           .epsilon(1e-12));
        // The molar concentration ISO 9613-1 takes is held; the relative humidity rises.
        CHECK(profile.at(z).molar_concentration_water_vapor()
                == Catch::Approx(concentration).epsilon(1e-12));
        if (i > 0)
            CHECK(profile.relative_humidity[i] > profile.relative_humidity[i - 1]);
    }

    // A stable class's inversion stops at the top of the surface layer; above
    // it every class lapses dry-adiabatically.
    SurfaceWind stable = wind;
    stable.stability = Stability::very_stable;
    const auto inverted = synthetic_profile(stable, surface, 1000.0);
    CHECK(inverted.at(200.0).temperature == Catch::Approx(surface.temperature + 0.0275 * 200.0));
    CHECK(inverted.at(1000.0).temperature
            == Catch::Approx(surface.temperature + 0.0275 * 200.0 - 0.0098 * 800.0));

    // Humid air saturates on its way up and is held at 100 percent.
    Atmosphere humid = surface;
    humid.relative_humidity = 95.0;
    const auto saturated = synthetic_profile(wind, humid, 2000.0);
    CHECK(saturated.relative_humidity.back() == 100.0);
    for (const double rh : saturated.relative_humidity)
        CHECK(rh <= 100.0);
}

TEST_CASE("the --wind direction is the bearing the wind blows from, as the file's is",
        "[atmosphere_profile][wind]")
{
    // From the west, it blows toward the east; from the north, toward the south.
    SurfaceWind wind = westerly();
    auto profile = synthetic_profile(wind, surface_air(), 300.0);
    CHECK(profile.wind_at(wind.height)[0] == Catch::Approx(wind.speed));
    CHECK(profile.wind_at(wind.height)[1] == Catch::Approx(0.0).margin(1e-12));
    wind.direction = 0.0;
    profile = synthetic_profile(wind, surface_air(), 300.0);
    CHECK(profile.wind_at(wind.height)[0] == Catch::Approx(0.0).margin(1e-12));
    CHECK(profile.wind_at(wind.height)[1] == Catch::Approx(-wind.speed));

    // And it is exactly what --atmosphere makes of the same wind written as
    // speed and direction, at every level.
    wind.direction = 37.0;
    profile = synthetic_profile(wind, surface_air(), 300.0);
    std::ostringstream csv;
    csv << std::setprecision(17) << "z_m,T_K,RH,wind_speed_mps,wind_direction\n";
    for (std::size_t i = 0; i < profile.height.size(); i++)
    {
        csv << profile.height[i] << ',' << profile.temperature[i] << ','
            << profile.relative_humidity[i] << ','
            << std::hypot(profile.wind_east[i], profile.wind_north[i]) << ",37\n";
    }
    TempCsv file(csv.str());
    const auto read = read_atmosphere_profile(file.path, meters, 0.0);
    REQUIRE(read.height.size() == profile.height.size());
    for (std::size_t i = 0; i < profile.height.size(); i++)
    {
        CHECK(read.wind_east[i] == Catch::Approx(profile.wind_east[i]).epsilon(1e-14).margin(1e-14));
        CHECK(read.wind_north[i] == Catch::Approx(profile.wind_north[i]).epsilon(1e-14).margin(1e-14));
    }
}

TEST_CASE("a calm --wind gives still air, with the lapse", "[atmosphere_profile][wind]")
{
    SurfaceWind calm = westerly();
    calm.speed = 0.0;
    const auto profile = synthetic_profile(calm, surface_air(), 300.0);
    CHECK_FALSE(profile.has_wind());
    CHECK(profile.temperature.back() < profile.temperature.front());
}

TEST_CASE("write_atmosphere_profile reads back as the profile it wrote", "[atmosphere_profile][wind]")
{
    const auto profile = synthetic_profile(westerly(), surface_air(), 600.0);
    for (const LengthUnit unit : {feet, meters})
    {
        TempCsv file("");
        write_atmosphere_profile(profile, file.path, unit);
        const auto read = read_atmosphere_profile(file.path, unit, 0.0);
        REQUIRE(read.height.size() == profile.height.size());
        REQUIRE(read.has_pressure());
        REQUIRE(read.has_wind());
        for (std::size_t i = 0; i < profile.height.size(); i++)
        {
            // Heights to an ULP of the unit conversion; everything else to the bit.
            CHECK(read.height[i] == Catch::Approx(profile.height[i]).epsilon(4e-16));
            CHECK(read.temperature[i] == profile.temperature[i]);
            CHECK(read.relative_humidity[i] == profile.relative_humidity[i]);
            CHECK(read.pressure[i] == profile.pressure[i]);
            CHECK(read.wind_east[i] == profile.wind_east[i]);
            CHECK(read.wind_north[i] == profile.wind_north[i]);
        }
    }
}

TEST_CASE("hydrostatic prefix agrees with independent walks from the bottom", "[atmosphere_profile]")
{
    // Nonuniform layers, near-isothermal segments, and ground below, inside,
    // or above the profile exercise both integrals and the surface reference.
    for (std::size_t count : {std::size_t(1), std::size_t(64)})
        for (double offset : {-1000.0, -100.0, 0.0, 123.0})
        {
            AtmosphereProfile p;
            for (std::size_t i = 0; i < count; i++)
            {
                p.height.push_back(offset + 7.0 * i + 0.05 * i * i);
                p.temperature.push_back(i % 4 == 0 ? 280.0 : i % 4 == 1 ? 280.0 + 1e-8
                                                                                      : 275.0 + 0.1 * i);
                p.relative_humidity.push_back(50.0);
            }
            auto segment = [](double dz, double a, double b) {
                return std::abs(b - a) < 1e-9 * a ? dz / a : dz * std::log(b / a) / (b - a);
            };
            auto walk = [&](double z) {
                if (z <= p.height.front())
                    return (z - p.height.front()) / p.temperature.front();
                double sum = 0.0;
                for (std::size_t i = 0; i + 1 < count; i++)
                {
                    if (z <= p.height[i + 1])
                    {
                        const double t = p.temperature[i] + (z - p.height[i])
                                / (p.height[i + 1] - p.height[i]) * (p.temperature[i + 1] - p.temperature[i]);
                        return sum + segment(z - p.height[i], p.temperature[i], t);
                    }
                    sum += segment(p.height[i + 1] - p.height[i], p.temperature[i], p.temperature[i + 1]);
                }
                return sum + (z - p.height.back()) / p.temperature.back();
            };
            const auto filled = p.with_hydrostatic_pressure(101.325);
            const double ground = walk(0.0);
            for (std::size_t i = 0; i < count; i++)
                REQUIRE(filled.pressure[i] == 101.325 * std::exp(-(standard_gravity / 287.05)
                        * (walk(p.height[i]) - ground)));
        }
}

TEST_CASE("profile fields share a cell without changing knots or endpoint holds", "[atmosphere_profile]")
{
    auto p = lapsing_profile();
    p.wind_east = p.relative_humidity;
    p.wind_north = p.temperature;
    auto reference = [&](const std::vector<double>& values, double z, bool logarithmic) {
        if (p.height.size() == 1 || z <= p.height.front())
            return values.front();
        if (z >= p.height.back())
            return values.back();
        const auto hi = std::upper_bound(p.height.begin(), p.height.end(), z);
        const std::size_t i = static_cast<std::size_t>(hi - p.height.begin()) - 1;
        const double f = (z - p.height[i]) / (p.height[i + 1] - p.height[i]);
        return logarithmic ? values[i] * std::pow(values[i + 1] / values[i], f)
                           : values[i] + f * (values[i + 1] - values[i]);
    };
    for (double z : {-10.0, 0.0, 0.01, 123.0, 250.0, 250.00001, 2999.99, 3000.0, 3100.0})
    {
        const auto air = p.at(z);
        const auto wind = p.wind_at(z);
        CHECK(air.temperature == reference(p.temperature, z, false));
        CHECK(air.relative_humidity == reference(p.relative_humidity, z, false));
        CHECK(air.pressure == reference(p.pressure, z, true));
        CHECK(wind[0] == reference(p.wind_east, z, false));
        CHECK(wind[1] == reference(p.wind_north, z, false));
    }
}

TEST_CASE("find_profile_columns requires temperature and humidity only for a run's profile",
        "[atmosphere_profile]")
{
    const TempCsv winds("z_ft,wind_speed_kt,wind_direction,extra\n46,10,270,1\n125,12,280,2\n");
    const LabeledCsv csv = read_labeled_csv(winds.path);
    const auto columns = find_profile_columns(csv, feet, winds.path, ProfileColumns::Require::none);
    CHECK(columns.height == 0);
    CHECK(columns.temperature == ProfileColumns::npos);
    CHECK(columns.humidity == ProfileColumns::npos);
    CHECK(columns.pressure == ProfileColumns::npos);
    CHECK(columns.speed == 1);
    CHECK(columns.direction == 2);
    CHECK(columns.speed_to_mps == knot);
    CHECK(columns.needed() == 3);
    CHECK(columns.used == std::vector<bool> {true, true, true, false});
    CHECK_THROWS_WITH(find_profile_columns(csv, feet, winds.path, ProfileColumns::Require::thermo),
            Catch::Matchers::ContainsSubstring("names no temperature column"));
    const TempCsv heightless("T_K,RH\n288,50\n");
    CHECK_THROWS_WITH(find_profile_columns(read_labeled_csv(heightless.path), feet,
                              heightless.path, ProfileColumns::Require::none, "--column"),
            Catch::Matchers::StartsWith("--column " + heightless.path + " names no height column"));
    CHECK_THROWS_WITH(read_atmosphere_profile(winds.path, feet, 0.0),
            Catch::Matchers::ContainsSubstring("--atmosphere " + winds.path
                    + " names no temperature column"));
}

TEST_CASE("geometric_levels is synthetic_profile's level placement", "[atmosphere_profile][wind]")
{
    SurfaceWind wind;
    wind.speed = 5.0;
    wind.height = 2.0;
    wind.roughness_length = 0.05;
    wind.surface_layer_depth = 150.0;
    wind.levels_per_decade = 12;
    Atmosphere surface;
    const auto profile = synthetic_profile(wind, surface, 400.0);
    auto levels = geometric_levels({420.0, 150.0, 2.0, 0.5}, 12);
    levels.insert(levels.begin(), 0.0);
    CHECK(levels == profile.height);
    // Anchors within 1e-9 of each other are one level, the lower.
    CHECK(geometric_levels({1.0, 10.0 * (1.0 + 1e-12), 10.0}, 1) == std::vector<double> {1.0, 10.0});
}
