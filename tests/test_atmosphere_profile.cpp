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

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "atmosphere_profile.hpp"

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
        profile.pressure.push_back(101.325 * std::pow(t / 288.15, 9.80665 / (0.0065 * 287.05)));
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
    const auto a = read_atmosphere_profile(kelvin_kpa.path, metres, 0.0);
    CHECK(a.temperature[0] == 288.15);
    CHECK(a.pressure[0] == 101.325);

    TempCsv celsius_pa("z, temperature, humidity, pressure\n0, 15, 50, 101325\n");
    const auto b = read_atmosphere_profile(celsius_pa.path, metres, 0.0);
    CHECK(b.temperature[0] == Catch::Approx(288.15));
    CHECK(b.pressure[0] == Catch::Approx(101.325));

    // "T (K)": the parenthesised unit is dropped by the header matching, so
    // this is the bare name too, and the values decide.
    TempCsv parenthesised("z (m), T (K), RH (%)\n0, 280, 20\n100, 279, 20\n");
    const auto c = read_atmosphere_profile(parenthesised.path, metres, 0.0);
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
    TempCsv metres_file("z_m, T_K, RH\n0, 288, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(metres_file.path, feet, 0.0), std::runtime_error);
    CHECK_NOTHROW(read_atmosphere_profile(metres_file.path, metres, 0.0));

    // Humidity given as a fraction.
    TempCsv fraction("z, T_K, RH\n0, 288, 0.5\n100, 287, 0.4\n");
    CHECK_THROWS_AS(read_atmosphere_profile(fraction.path, metres, 0.0), std::runtime_error);

    // Two levels at one height.
    TempCsv duplicate("z, T_K, RH\n0, 288, 50\n0, 287, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(duplicate.path, metres, 0.0), std::runtime_error);

    // Kelvin labelled as Celsius.
    TempCsv mislabelled("z, T_C, RH\n0, 288, 50\n");
    CHECK_THROWS_AS(read_atmosphere_profile(mislabelled.path, metres, 0.0), std::runtime_error);

    // No humidity at all.
    TempCsv no_humidity("z, T_K\n0, 288\n");
    CHECK_THROWS_AS(read_atmosphere_profile(no_humidity.path, metres, 0.0), std::runtime_error);

    // Pressure rising with height.
    TempCsv rising("z, T_K, RH, p_kPa\n0, 288, 50, 90\n100, 287, 50, 95\n");
    CHECK_THROWS_AS(read_atmosphere_profile(rising.path, metres, 0.0), std::runtime_error);
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
    const double g_over_r = 9.80665 / 287.05;

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
    // Either side of the millimetre where the difference quotient hands over
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
