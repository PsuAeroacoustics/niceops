// The propagation loss table (propagation_table.hpp): its file, its lookup,
// and the edges of what it covers.  CMakeLists.txt's propagation_table_* tests
// run it end to end, against --ground_effect.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <vector>

#include "propagation_table.hpp"

using namespace niceops;

namespace {

// A table whose value is linear in every axis but bearing, and a cosine of
// bearing, so linear interpolation is exact away from bearing and checkable
// along it.
double planted(double bearing, double source, double receiver, double range, double band)
{
    return 0.01 * source - 0.5 * receiver - 0.001 * range + 0.002 * band + 3.0 * std::cos(bearing * M_PI / 180.0);
}

PropagationTable make_table(std::size_t bearings = 8)
{
    PropagationTable::Axes axes;
    for (std::size_t b = 0; b < bearings; b++)
        axes.bearings.push_back(360.0 * static_cast<double>(b) / static_cast<double>(bearings));
    axes.source_heights = {10.0, 50.0, 100.0};
    axes.receiver_heights = {0.0, 1.2, 4.0};
    axes.ranges = {10.0, 20.0, 40.0, 80.0};
    axes.bands = {100.0, 125.0, 160.0};
    std::vector<double> values;
    for (const double b : axes.bearings)
        for (const double s : axes.source_heights)
            for (const double h : axes.receiver_heights)
                for (const double r : axes.ranges)
                    for (const double f : axes.bands)
                        values.push_back(planted(b, s, h, r, f));
    return PropagationTable(axes, values, "planted", false, {{"note", "a test"}});
}

} // namespace

TEST_CASE("a propagation table is linear between its points in all but bearing", "[propagation_table]")
{
    const auto table = make_table();
    const auto bands = table.match_bands({100.0, 125.0, 160.0}, 0, 3);
    double out[3] = {};
    // On a bearing, everything else between points: exact.
    REQUIRE(table.lookup(90.0, 33.0, 0.7, 55.0, bands, out));
    for (int k = 0; k < 3; k++)
        CHECK(out[k] == Catch::Approx(planted(90.0, 33.0, 0.7, 55.0, table.axes().bands[k])).margin(1e-12));
    // Between bearings: linear between the two either side, wrapping at 360.
    REQUIRE(table.lookup(-22.5, 10.0, 0.0, 10.0, bands, out));
    const double expected = 0.5 * (planted(315.0, 10.0, 0.0, 10.0, 100.0) + planted(0.0, 10.0, 0.0, 10.0, 100.0));
    CHECK(out[0] == Catch::Approx(expected).margin(1e-12));
    REQUIRE(table.lookup(337.5, 10.0, 0.0, 10.0, bands, out));
    CHECK(out[0] == Catch::Approx(expected).margin(1e-12));
    // The ends of each axis are inside it.
    CHECK(table.lookup(0.0, 100.0, 4.0, 80.0, bands, out));
    CHECK(out[2] == Catch::Approx(planted(0.0, 100.0, 4.0, 80.0, 160.0)).margin(1e-12));
}

TEST_CASE("a propagation table covers only what it holds", "[propagation_table]")
{
    const auto table = make_table();
    const auto bands = table.match_bands({100.0, 125.0, 160.0}, 0, 3);
    double out[3] = {7.0, 7.0, 7.0};
    CHECK_FALSE(table.lookup(0.0, 5.0, 1.2, 20.0, bands, out)); // below the lowest source
    CHECK_FALSE(table.lookup(0.0, 101.0, 1.2, 20.0, bands, out)); // above the highest
    CHECK_FALSE(table.lookup(0.0, 50.0, 1.2, 9.0, bands, out)); // nearer than the first range
    CHECK_FALSE(table.lookup(0.0, 50.0, 1.2, 81.0, bands, out)); // beyond the last
    CHECK_FALSE(table.lookup(0.0, 50.0, 5.0, 20.0, bands, out)); // receiver above the table's
    CHECK(out[0] == 7.0); // and nothing is written
    CHECK(table.covers_receiver(2.0));
    CHECK_FALSE(table.covers_receiver(4.5));
}

TEST_CASE("a NaN in a propagation table is a gap, not a value", "[propagation_table]")
{
    auto axes = make_table().axes();
    std::vector<double> values;
    for (const double b : axes.bearings)
        for (const double s : axes.source_heights)
            for (const double h : axes.receiver_heights)
                for (const double r : axes.ranges)
                    for (const double f : axes.bands)
                    {
                        // The PE's way: no value within 20 m at the highest source.
                        const bool steep = s == 100.0 && r <= 20.0;
                        values.push_back(steep ? std::numeric_limits<double>::quiet_NaN() : planted(b, s, h, r, f));
                    }
    const PropagationTable table(axes, values, "gappy", false);
    const auto bands = table.match_bands({100.0, 125.0, 160.0}, 0, 3);
    double out[3] = {};
    CHECK_FALSE(table.lookup(0.0, 75.0, 1.2, 15.0, bands, out)); // a corner is a NaN
    CHECK(table.lookup(0.0, 75.0, 1.2, 50.0, bands, out)); // past it: 40 to 80 m
    // On the lower source height exactly, the NaN above has no weight.
    CHECK(table.lookup(0.0, 50.0, 1.2, 15.0, bands, out));
    CHECK(out[1] == Catch::Approx(planted(0.0, 50.0, 1.2, 15.0, 125.0)).margin(1e-12));
    // And the bound ignores it.
    const auto ratio = table.max_power_ratio(bands);
    CHECK(std::isfinite(ratio[0]));
}

TEST_CASE("one bearing means every bearing", "[propagation_table]")
{
    const auto table = make_table(1);
    CHECK(table.isotropic());
    const auto bands = table.match_bands({100.0, 125.0, 160.0}, 0, 3);
    double a[3] = {}, b[3] = {};
    REQUIRE(table.lookup(0.0, 30.0, 1.0, 30.0, bands, a));
    REQUIRE(table.lookup(213.0, 30.0, 1.0, 30.0, bands, b));
    CHECK(a[1] == b[1]);
}

TEST_CASE("a run's bands are found in the table, or the run is refused", "[propagation_table]")
{
    const auto table = make_table();
    // A database's centers rounded differently, and bands outside the limits
    // the table does not have.
    const std::vector<double> run = {80.0, 100.4, 125.0, 159.5, 200.0};
    const auto bands = table.match_bands(run, 1, 4);
    CHECK(bands.column[1] == 0);
    CHECK(bands.column[3] == 2);
    double out[5] = {-1.0, -1.0, -1.0, -1.0, -1.0};
    REQUIRE(table.lookup(0.0, 50.0, 1.2, 20.0, bands, out));
    CHECK(out[0] == -1.0);
    CHECK(out[4] == -1.0);
    CHECK(out[3] == Catch::Approx(planted(0.0, 50.0, 1.2, 20.0, 160.0)).margin(1e-12));
    CHECK_THROWS_AS(table.match_bands(run, 0, 5), std::runtime_error);
}

TEST_CASE("a propagation table survives its file", "[propagation_table]")
{
    const auto table = make_table();
    const auto path = (std::filesystem::temp_directory_path() / "niceops_propagation_table_test.nc").string();
    table.save(path);
    const auto loaded = PropagationTable::load(path);
    std::filesystem::remove(path);
    CHECK(loaded.method() == "planted");
    CHECK_FALSE(loaded.includes_absorption());
    CHECK(loaded.provenance().at("note") == "a test");
    CHECK(loaded.axes().bearings == table.axes().bearings);
    CHECK(loaded.axes().ranges == table.axes().ranges);
    for (std::size_t s = 0; s < 3; s++)
        for (std::size_t r = 0; r < 4; r++)
            CHECK(loaded.at(3, s, 1, r, 2) == table.at(3, s, 1, r, 2));
}

TEST_CASE("a propagation table refuses axes it cannot interpolate", "[propagation_table]")
{
    auto axes = make_table().axes();
    const std::size_t count = axes.bearings.size() * 3 * 3 * 4 * 3;
    auto bad = axes;
    bad.ranges = {10.0, 40.0, 20.0, 80.0};
    CHECK_THROWS_AS(PropagationTable(bad, std::vector<double>(count, 0.0), "", false), std::invalid_argument);
    bad = axes;
    bad.bearings[1] = 50.0;
    CHECK_THROWS_AS(PropagationTable(bad, std::vector<double>(count, 0.0), "", false), std::invalid_argument);
    CHECK_THROWS_AS(PropagationTable(axes, std::vector<double>(count - 1, 0.0), "", false), std::invalid_argument);
}
