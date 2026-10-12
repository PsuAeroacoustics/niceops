#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "screening.hpp"
#include "urban.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

struct Direction
{
    double azimuth = 0.0;   // compass, degrees, toward the source
    double elevation = 0.0; // degrees
};

// A local-frame unit vector toward compass bearing `azimuth` at `elevation`,
// for a run whose +x points along `frame_bearing`.
Vec3 local_unit(const Direction& d, double frame_bearing)
{
    const double a = d.azimuth * pi / 180.0;
    const double e = d.elevation * pi / 180.0;
    const double b = frame_bearing * pi / 180.0;
    const double east = std::sin(a) * std::cos(e);
    const double north = std::cos(a) * std::cos(e);
    return Vec3(east * std::sin(b) + north * std::cos(b), -east * std::cos(b) + north * std::sin(b), std::sin(e));
}

// How close the trace's wall crossings come to the walls' tops: a far source
// is not a plane wave there, so a sample that grazes a top is ambiguous.
double top_margin(const IdealCanyon& c, double v_c, double u_z)
{
    double y = c.offset, z = c.microphone_height, margin = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 400 && v_c != 0.0; k++)
    {
        const bool b = v_c > 0.0;
        z += u_z * (b ? c.width - y : y) / std::abs(v_c);
        const double h = b ? c.height_b : c.height_a;
        margin = std::min(margin, std::abs(z - h));
        if (z > h)
            break;
        y = b ? c.width : 0.0;
        v_c = -v_c;
    }
    return margin;
}

} // namespace

TEST_CASE("facade mirror, side and blocking", "[urban]")
{
    const Facade f {0.0, 0.0, 10.0, 0.0, 0.0, 5.0}; // along +x, left normal +y
    REQUIRE(facade_side(f, Vec3(3.0, 2.0, 1.0)) == Catch::Approx(2.0));
    REQUIRE(facade_side(f, Vec3(3.0, -4.0, 1.0)) == Catch::Approx(-4.0));
    const Vec3 m = mirror_in(f, Vec3(3.0, 2.0, 1.0));
    REQUIRE(m.x == Catch::Approx(3.0));
    REQUIRE(m.y == Catch::Approx(-2.0));
    REQUIRE(m.z == 1.0);
    const Vec3 back = mirror_in(f, m);
    REQUIRE(back.y == Catch::Approx(2.0));
    const Vec3 v = mirror_direction(f, Vec3(1.0, 1.0, 0.5));
    REQUIRE(v.x == Catch::Approx(1.0));
    REQUIRE(v.y == Catch::Approx(-1.0));
    REQUIRE(v.z == 0.5);

    // Through the rectangle, over it, past its end, and ending on it.
    REQUIRE(blocks(f, Vec3(5.0, 1.0, 2.0), Vec3(5.0, -1.0, 2.0)));
    REQUIRE_FALSE(blocks(f, Vec3(5.0, 1.0, 6.0), Vec3(5.0, -1.0, 6.0)));
    REQUIRE_FALSE(blocks(f, Vec3(11.0, 1.0, 2.0), Vec3(11.0, -1.0, 2.0)));
    REQUIRE_FALSE(blocks(f, Vec3(5.0, 1.0, 2.0), Vec3(5.0, 0.0, 2.0)));
    REQUIRE_FALSE(blocks(f, Vec3(5.0, 1.0, 2.0), Vec3(5.0, 3.0, 2.0)));
}

TEST_CASE("first_order_image agrees with the image-source reference", "[urban]")
{
    std::mt19937 rng(20261011);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    int compared = 0, found = 0;
    for (int trial = 0; trial < 2000; trial++)
    {
        IdealCanyon c;
        c.width = 8.0 + 30.0 * unit(rng);
        c.offset = c.width * (0.05 + 0.9 * unit(rng));
        c.height_a = 3.0 + 40.0 * unit(rng);
        c.height_b = 3.0 + 40.0 * unit(rng);
        c.microphone_height = 0.5 + 4.0 * unit(rng);
        const double axis = 360.0 * unit(rng);
        const double frame = 360.0 * unit(rng);
        const Vec3 mic(100.0 * unit(rng), 100.0 * unit(rng), 10.0 + c.microphone_height);
        const auto facades = ideal_canyon_facades(c, mic, 10.0, axis, frame, 400.0);
        const Direction d {360.0 * unit(rng), 5.0 + 80.0 * unit(rng)};
        const Vec3 source = mic + local_unit(d, frame) * (30.0 + 600.0 * unit(rng));
        const auto paths = specular_paths(facades, source, mic, 1);
        for (std::size_t k = 0; k < facades.size(); k++)
        {
            Vec3 image;
            const bool valid = first_order_image(facades, k, source, mic, image);
            const bool reference = std::any_of(paths.begin(), paths.end(),
                    [&](const SpecularPath& p) { return p.facades == std::vector<std::size_t> {k}; });
            REQUIRE(valid == reference);
            compared++;
            if (valid)
            {
                found++;
                REQUIRE((image - mirror_in(facades[k], mic)).norm() < 1e-9);
                const auto it = std::find_if(paths.begin(), paths.end(),
                        [&](const SpecularPath& p) { return p.facades == std::vector<std::size_t> {k}; });
                REQUIRE((image - source).norm() == Catch::Approx(it->length).epsilon(1e-12));
            }
        }
    }
    REQUIRE(compared == 4000);
    // Both outcomes occur, or the comparison proves little.
    REQUIRE(found > 200);
    REQUIRE(found < compared - 200);
}

// The claim the table rests on: under plane-wave incidence an ideal canyon
// delivers at most two specular paths to a point, one per cross-canyon
// direction, and trace_canyon says which.  Checked against the brute-force
// image-source search, with the source far enough to be a plane wave.
TEST_CASE("a far source reaches a canyon microphone by at most two specular paths", "[urban]")
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    int checked = 0, two = 0, shielded = 0;
    for (int trial = 0; trial < 3000; trial++)
    {
        IdealCanyon c;
        c.width = 8.0 + 30.0 * unit(rng);
        c.offset = c.width * (0.05 + 0.9 * unit(rng));
        c.height_a = 3.0 + 40.0 * unit(rng);
        c.height_b = 3.0 + 40.0 * unit(rng);
        c.microphone_height = 0.5 + 4.0 * unit(rng);
        const double axis = 360.0 * unit(rng);
        const double frame = 360.0 * unit(rng);
        const Direction d {360.0 * unit(rng), 2.0 + 85.0 * unit(rng)};
        const double a = d.azimuth * pi / 180.0, e = d.elevation * pi / 180.0, x = axis * pi / 180.0;
        const double u_c = std::cos(e) * std::sin(a - x);
        const double u_z = std::sin(e);
        // Off-axis enough that every reflection lands well inside the long
        // walls, and clear of the walls' tops.
        if (std::abs(u_c) < 0.08 || top_margin(c, u_c, u_z) < 1e-3 || top_margin(c, -u_c, u_z) < 1e-3)
            continue;
        const Vec3 mic(0.0, 0.0, c.microphone_height);
        const auto facades = ideal_canyon_facades(c, mic, 0.0, axis, frame, 4.0e5);
        const Vec3 source = mic + local_unit(d, frame) * 1.0e8;

        const CanyonTrace direct = trace_canyon(c, u_c, u_z);
        const CanyonTrace mirror = trace_canyon(c, -u_c, u_z);
        const int max_order = std::max(direct.bounces, mirror.bounces) + 3;
        if (max_order > 40)
            continue;
        const auto paths = specular_paths(facades, source, mic, static_cast<std::size_t>(max_order));
        REQUIRE(paths.size() <= 2);

        std::vector<std::size_t> expected;
        if (direct.escaped && direct.bounces % 2 == 0)
            expected.push_back(static_cast<std::size_t>(direct.bounces));
        if (mirror.escaped && mirror.bounces % 2 == 1)
            expected.push_back(static_cast<std::size_t>(mirror.bounces));
        std::vector<std::size_t> orders;
        for (const auto& p : paths)
            orders.push_back(p.facades.size());
        std::sort(expected.begin(), expected.end());
        std::sort(orders.begin(), orders.end());
        REQUIRE(orders == expected);
        checked++;
        two += paths.size() == 2;
        shielded += !(direct.escaped && direct.bounces == 0);
    }
    REQUIRE(checked > 1500);
    REQUIRE(two > 100);
    REQUIRE(shielded > 100);
}

TEST_CASE("edge_path_difference is the minimum over the edge", "[urban]")
{
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int trial = 0; trial < 200; trial++)
    {
        // Canyon coordinates: x across, y along the edge, z up.
        const double dy = -30.0 + 60.0 * unit(rng);
        const double dz = 1.0 + 30.0 * unit(rng);
        const double e = (2.0 + 80.0 * unit(rng)) * pi / 180.0;
        const double az = 2.0 * pi * unit(rng);
        const double u_c = std::cos(e) * std::cos(az), u_a = std::cos(e) * std::sin(az), u_z = std::sin(e);
        // f(t) = |E(t) - R| - u . (E(t) - R) along the edge, convex in t.
        const auto f = [&](double t) { return std::sqrt(dy * dy + t * t + dz * dz) - (u_c * dy + u_a * t + u_z * dz); };
        double lo = -1e6, hi = 1e6;
        for (int k = 0; k < 400; k++)
        {
            const double m1 = lo + (hi - lo) / 3.0, m2 = hi - (hi - lo) / 3.0;
            (f(m1) < f(m2) ? hi : lo) = f(m1) < f(m2) ? m2 : m1;
        }
        REQUIRE(edge_path_difference(dy, dz, u_c, u_z) == Catch::Approx(f(0.5 * (lo + hi))).margin(1e-7));
    }
}

TEST_CASE("canyon_ratio's limits", "[urban]")
{
    const std::vector<double> bands {125.0, 1000.0, 4000.0};
    const std::vector<double> rigid(3, 1.0), soft(3, 0.5);
    std::vector<double> out(3);
    IdealCanyon c {20.0, 6.0, 15.0, 15.0, 1.5};

    // Along the axis, or above both walls, the direct path is open and alone.
    canyon_ratio(c, 0.0, 0.3, rigid.data(), rigid.data(), bands, 343.0, true, out.data());
    REQUIRE(out == std::vector<double>(3, 1.0));
    IdealCanyon roof = c;
    roof.microphone_height = 16.0;
    canyon_ratio(roof, 0.7, 0.2, rigid.data(), rigid.data(), bands, 343.0, true, out.data());
    REQUIRE(out == std::vector<double>(3, 1.0));

    // Overhead-ish from B's side, rigid walls: the direct path is open and its
    // mirror's family reflects once, which the run evaluates itself.
    const double e = 70.0 * pi / 180.0;
    canyon_ratio(c, std::cos(e), std::sin(e), rigid.data(), rigid.data(), bands, 343.0, true, out.data());
    REQUIRE(out == std::vector<double>(3, 1.0));

    // Low over B: the direct path is blocked, so the diffracted path stands in,
    // falling with frequency, below 10 lg 3 and never under the 25 dB cap.
    const double low = 10.0 * pi / 180.0;
    canyon_ratio(c, std::cos(low), std::sin(low), soft.data(), soft.data(), bands, 343.0, true, out.data());
    std::vector<double> specular(3);
    canyon_ratio(c, std::cos(low), std::sin(low), soft.data(), soft.data(), bands, 343.0, false, specular.data());
    for (std::size_t k = 0; k < 3; k++)
    {
        const double diffracted = out[k] - specular[k];
        REQUIRE(diffracted <= 1.0 / cnossos_floor + 1e-15);
        REQUIRE(diffracted >= std::pow(10.0, -0.1 * cnossos_cap) - 1e-15);
        if (k > 0)
            REQUIRE(out[k] - specular[k] <= out[k - 1] - specular[k - 1]);
    }
    // The lowest band is above the cap, so the fall is not just the cap.
    REQUIRE(out[0] - specular[0] > std::pow(10.0, -0.1 * cnossos_cap) * 1.5);
}

// What separating the first-order paths from the table costs at short range:
// the explicit paths plus a plane-wave table, against every specular path of
// a point source.  Exact in the far field, where the table's model holds.
TEST_CASE("explicit images plus the plane-wave table converge on the image lattice", "[urban]")
{
    const IdealCanyon c {20.0, 6.0, 15.0, 18.0, 1.5};
    const double axis = 30.0, frame = 90.0, reflect = 0.8;
    const std::vector<double> bands {1000.0};
    const std::vector<double> r(1, reflect);
    const Vec3 mic(0.0, 0.0, c.microphone_height);
    const auto facades = ideal_canyon_facades(c, mic, 0.0, axis, frame, 4.0e6);
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (double range : {1.0e5, 1.0e6})
    {
        int compared = 0;
        for (int trial = 0; trial < 300; trial++)
        {
            const Direction d {360.0 * unit(rng), 20.0 + 65.0 * unit(rng)};
            const double a = d.azimuth * pi / 180.0, e = d.elevation * pi / 180.0, x = axis * pi / 180.0;
            const double u_c = std::cos(e) * std::sin(a - x), u_z = std::sin(e);
            if (std::abs(u_c) < 0.1 || top_margin(c, u_c, u_z) < 0.5 || top_margin(c, -u_c, u_z) < 0.5)
                continue;
            const Vec3 source = mic + local_unit(d, frame) * range;
            double lattice = 0.0;
            for (const auto& p : specular_paths(facades, source, mic, 30))
                lattice += std::pow(reflect, static_cast<double>(p.facades.size())) / (p.length * p.length);
            double ratio = 0.0;
            canyon_ratio(c, u_c, u_z, r.data(), r.data(), bands, 343.0, false, &ratio);
            double split = ratio / (range * range);
            for (std::size_t k = 0; k < facades.size(); k++)
            {
                Vec3 image;
                if (first_order_image(facades, k, source, mic, image))
                    split += reflect / (image - source).squaredNorm();
            }
            if (lattice == 0.0)
            {
                REQUIRE(split == 0.0);
                continue;
            }
            REQUIRE(std::abs(10.0 * std::log10(split / lattice)) < (range > 5.0e5 ? 2e-3 : 1e-2));
            compared++;
        }
        REQUIRE(compared > 100);
    }
}

// Not run by default: the same comparison at helicopter ranges, printed, for
// the separability study (docs/FIDELITY.md).
TEST_CASE("explicit images plus the table against the lattice by range", "[.urban_separability]")
{
    const IdealCanyon c {20.0, 6.0, 15.0, 18.0, 1.5};
    const double axis = 30.0, frame = 90.0, reflect = 0.8;
    const std::vector<double> bands {1000.0};
    const std::vector<double> r(1, reflect);
    const Vec3 mic(0.0, 0.0, c.microphone_height);
    const auto facades = ideal_canyon_facades(c, mic, 0.0, axis, frame, 4.0e4);
    // "split" is the run's form, the first-order images explicit; "table" puts
    // them in the table too, every path keyed on the arrival direction alone.
    std::printf("range_m  split_median  split_p95  split_max  table_median  table_p95  table_max  pairs\n");
    for (double range : {30.0, 60.0, 120.0, 300.0, 1000.0, 3000.0, 10000.0})
    {
        std::mt19937 rng(5);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        std::vector<double> errors, table_errors;
        for (int trial = 0; trial < 4000; trial++)
        {
            const Direction d {360.0 * unit(rng), 5.0 + 85.0 * unit(rng)};
            const double a = d.azimuth * pi / 180.0, e = d.elevation * pi / 180.0, x = axis * pi / 180.0;
            const double u_c = std::cos(e) * std::sin(a - x), u_z = std::sin(e);
            const Vec3 source = mic + local_unit(d, frame) * range;
            double lattice = 0.0;
            for (const auto& p : specular_paths(facades, source, mic, 24))
                lattice += std::pow(reflect, static_cast<double>(p.facades.size())) / (p.length * p.length);
            double ratio = 0.0;
            canyon_ratio(c, u_c, u_z, r.data(), r.data(), bands, 343.0, false, &ratio);
            double split = ratio / (range * range);
            for (std::size_t k = 0; k < facades.size(); k++)
            {
                Vec3 image;
                if (first_order_image(facades, k, source, mic, image))
                    split += reflect / (image - source).squaredNorm();
            }
            const CanyonTrace mirror = trace_canyon(c, -u_c, u_z);
            const double first_order = u_c != 0.0 && mirror.escaped && mirror.bounces == 1 ? reflect : 0.0;
            const double table_only = (ratio + first_order) / (range * range);
            if (lattice > 0.0 && split > 0.0)
                errors.push_back(std::abs(10.0 * std::log10(split / lattice)));
            if (lattice > 0.0 && table_only > 0.0)
                table_errors.push_back(std::abs(10.0 * std::log10(table_only / lattice)));
        }
        std::sort(errors.begin(), errors.end());
        std::sort(table_errors.begin(), table_errors.end());
        const auto at = [](const std::vector<double>& v, double q) { return v[static_cast<std::size_t>(q * (v.size() - 1))]; };
        std::printf("%7.0f  %12.3f  %9.3f  %9.3f  %12.3f  %9.3f  %9.3f  %5zu\n", range, at(errors, 0.5),
                at(errors, 0.95), errors.back(), at(table_errors, 0.5), at(table_errors, 0.95), table_errors.back(),
                errors.size());
    }
}

TEST_CASE("UrbanTable lookup, bounds and file round trip", "[urban]")
{
    UrbanTable::Axes axes;
    axes.azimuths = {0.0, 90.0, 180.0, 270.0};
    axes.elevations = {0.0, 10.0, 40.0};
    axes.bands = {500.0, 1000.0};
    UrbanTable::Receiver near {Vec3(1.0, 2.0, 1.5), {{-5.0, 0.0, 5.0, 0.0, 0.0, 12.0}}, {0.9, 0.8}};
    UrbanTable::Receiver open {Vec3(50.0, 2.0, 1.5), {}, {}};
    // gain = azimuth index + 10 * elevation index + 100 * band - 120 for
    // receiver 0, whose largest is +3 dB in the second band; receiver 1 is
    // 0 dB everywhere.
    std::vector<double> gains;
    for (int a = 0; a < 4; a++)
        for (int e = 0; e < 3; e++)
            for (int b = 0; b < 2; b++)
                gains.push_back(a + 10.0 * e + 100.0 * b - 120.0);
    gains.resize(gains.size() * 2, 0.0);
    const UrbanTable table(axes, {near, open}, gains, 210.0, "test", {{"note", "unit test"}});

    const auto map = table.match_bands({250.0, 500.0, 1000.0, 2000.0}, 1, 3);
    double out[4] = {-1.0, -1.0, -1.0, -1.0};
    REQUIRE(table.lookup(0, 90.0, 10.0, map, out));
    REQUIRE(out[1] == 1.0 + 10.0 - 120.0);
    REQUIRE(out[2] == 1.0 + 10.0 - 20.0);
    REQUIRE(out[0] == -1.0);
    REQUIRE(out[3] == -1.0);
    // Halfway in both: the mean of the four corners.
    REQUIRE(table.lookup(0, 45.0, 25.0, map, out));
    REQUIRE(out[1] == Catch::Approx(0.5 + 15.0 - 120.0));
    // Azimuth wraps from 270 back to 0.
    REQUIRE(table.lookup(0, 315.0, 0.0, map, out));
    REQUIRE(out[1] == Catch::Approx(1.5 - 120.0));
    REQUIRE(table.lookup(0, -45.0, 0.0, map, out));
    REQUIRE(out[1] == Catch::Approx(1.5 - 120.0));
    // Outside the elevation axis: no value.
    REQUIRE_FALSE(table.lookup(0, 0.0, -1.0, map, out));
    REQUIRE_FALSE(table.lookup(0, 0.0, 41.0, map, out));
    REQUIRE(table.lookup(1, 123.0, 33.0, map, out));
    REQUIRE(out[1] == 0.0);
    REQUIRE(out[2] == 0.0);

    const auto ratio = table.max_power_ratio(map);
    // Outside the run's bands, 1; a band whose gains are all at most 0 dB, 1;
    // otherwise the largest, +3 dB.
    REQUIRE(ratio[0] == 1.0);
    REQUIRE(ratio[1] == 1.0);
    REQUIRE(ratio[2] == Catch::Approx(std::pow(10.0, 0.3)).epsilon(1e-12));
    REQUIRE(ratio[3] == 1.0);

    REQUIRE(table.find_receiver(Vec3(1.0, 2.0005, 1.5), 1e-3) == 0);
    REQUIRE(table.find_receiver(Vec3(1.0, 2.01, 1.5), 1e-3) == 2);
    REQUIRE(table.find_receiver(Vec3(50.0, 2.0, 1.5), 1e-3) == 1);
    REQUIRE_THROWS_AS(table.match_bands({800.0}, 0, 1), std::runtime_error);

    const auto path = (std::filesystem::temp_directory_path() / "niceops_urban_table_test.nc").string();
    table.save(path);
    const UrbanTable back = UrbanTable::load(path);
    std::filesystem::remove(path);
    REQUIRE(back.axes().azimuths == axes.azimuths);
    REQUIRE(back.axes().elevations == axes.elevations);
    REQUIRE(back.axes().bands == axes.bands);
    REQUIRE(back.frame_bearing() == 210.0);
    REQUIRE(back.method() == "test");
    REQUIRE(back.provenance().at("note") == "unit test");
    REQUIRE(back.receivers().size() == 2);
    REQUIRE(back.receivers()[0].facades.size() == 1);
    REQUIRE(back.receivers()[0].facades[0].top == 12.0);
    REQUIRE(back.receivers()[0].reflection == std::vector<double> {0.9, 0.8});
    REQUIRE(back.receivers()[1].facades.empty());
    for (std::size_t a = 0; a < 4; a++)
        for (std::size_t e = 0; e < 3; e++)
            for (std::size_t b = 0; b < 2; b++)
                REQUIRE(back.gain(0, a, e, b) == table.gain(0, a, e, b));
}

TEST_CASE("UrbanTable refuses a malformed table", "[urban]")
{
    UrbanTable::Axes axes {{0.0, 120.0, 240.0}, {0.0, 45.0}, {1000.0}};
    const UrbanTable::Receiver r {Vec3(0.0, 0.0, 1.5), {}, {}};
    const std::vector<double> gains(3 * 2, 0.0);
    REQUIRE_NOTHROW(UrbanTable(axes, {r}, gains, 90.0, "ok"));
    auto bad_axes = axes;
    bad_axes.azimuths = {0.0, 100.0, 240.0};
    REQUIRE_THROWS_AS(UrbanTable(bad_axes, {r}, gains, 90.0, "x"), std::invalid_argument);
    bad_axes = axes;
    bad_axes.elevations = {45.0, 0.0};
    REQUIRE_THROWS_AS(UrbanTable(bad_axes, {r}, gains, 90.0, "x"), std::invalid_argument);
    REQUIRE_THROWS_AS(UrbanTable(axes, {r}, std::vector<double>(5, 0.0), 90.0, "x"), std::invalid_argument);
    auto infinite = gains;
    infinite[2] = -std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(UrbanTable(axes, {r}, infinite, 90.0, "x"), std::invalid_argument);
    UrbanTable::Receiver loud {Vec3(0.0, 0.0, 1.5), {{0.0, 0.0, 1.0, 0.0, 0.0, 5.0}}, {1.2}};
    REQUIRE_THROWS_AS(UrbanTable(axes, {loud}, gains, 90.0, "x"), std::invalid_argument);
    UrbanTable::Receiver point {Vec3(0.0, 0.0, 1.5), {{1.0, 1.0, 1.0, 1.0, 0.0, 5.0}}, {0.5}};
    REQUIRE_THROWS_AS(UrbanTable(axes, {point}, gains, 90.0, "x"), std::invalid_argument);
    REQUIRE_THROWS_AS(UrbanTable(axes, {r}, gains, std::nan(""), "x"), std::invalid_argument);
}

TEST_CASE("compass_direction inverts the local unit vector", "[urban]")
{
    for (double frame : {0.0, 90.0, 210.0, 333.0})
        for (double azimuth : {0.0, 45.0, 180.0, 300.0})
            for (double elevation : {0.0, 30.0, 85.0})
            {
                const auto d = compass_direction(local_unit({azimuth, elevation}, frame) * 250.0, frame);
                REQUIRE(std::fmod(d.azimuth - azimuth + 540.0, 360.0) - 180.0 == Catch::Approx(0.0).margin(1e-9));
                REQUIRE(d.elevation == Catch::Approx(elevation).margin(1e-9));
            }
}

// Not run by default: writes tables/urban_bump_test.nc for the CTest
// urban_table_bump_windowed_is_exact, from the canyon table the CTests build
// (tables/urban_canyon_test.nc): +25 dB at elevations 3 to 20 degrees, where
// the demo track's states are otherwise quiet, so a peak search whose bound
// left out the table's gain would skip the loudest states.
TEST_CASE("write the urban bump table", "[.urban_bump_table]")
{
    const UrbanTable canyon = UrbanTable::load("tables/urban_canyon_test.nc");
    const auto& axes = canyon.axes();
    std::vector<double> gains;
    for (std::size_t r = 0; r < canyon.receivers().size(); r++)
        for (std::size_t a = 0; a < axes.azimuths.size(); a++)
            for (std::size_t e = 0; e < axes.elevations.size(); e++)
                for (std::size_t b = 0; b < axes.bands.size(); b++)
                {
                    const double elevation = axes.elevations[e];
                    const double bump = elevation >= 3.0 && elevation <= 20.0 ? 25.0 : 0.0;
                    gains.push_back(canyon.gain(r, a, e, b) + bump);
                }
    const UrbanTable bump(axes, canyon.receivers(), gains, canyon.frame_bearing(), "ideal_canyon_bump",
            canyon.provenance());
    bump.save("tables/urban_bump_test.nc");
    REQUIRE(UrbanTable::load("tables/urban_bump_test.nc").receivers().size() == canyon.receivers().size());
}
