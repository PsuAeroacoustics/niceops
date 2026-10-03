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

namespace {

// The cell and weight by plain search, independent of the table's own.
bool reference_cell(const std::vector<double>& axis, double x, std::size_t& lower, double& fraction)
{
    if (axis.size() == 1)
    {
        lower = 0;
        fraction = 0.0;
        return std::abs(x - axis[0]) <= 1e-6 * std::max(1.0, std::abs(axis[0]));
    }
    if (x < axis.front() || x > axis.back())
        return false;
    lower = 0;
    while (lower + 2 < axis.size() && axis[lower + 1] <= x)
        lower++;
    fraction = (x - axis[lower]) / (axis[lower + 1] - axis[lower]);
    return true;
}

} // namespace

TEST_CASE("lookup finds the same cell on even and uneven axes, at the knots and between", "[propagation_table]")
{
    // Even axes are indexed by arithmetic; uneven ones by search.  Both must give
    // the cell a search gives, including exactly at a knot and at the axis's ends,
    // where round-off in the arithmetic estimate could put it a cell off.
    for (const bool even : {true, false})
    {
        PropagationTable::Axes axes;
        for (std::size_t b = 0; b < 12; b++)
            axes.bearings.push_back(30.0 * static_cast<double>(b));
        for (int k = 0; k < 17; k++)
        {
            axes.source_heights.push_back(even ? 3.0 + 7.0 * k : 3.0 + 7.0 * k * (1.0 + 0.02 * k));
            axes.receiver_heights.push_back(even ? 0.1 * k : 0.1 * k * k);
        }
        for (int k = 0; k < 40; k++)
            axes.ranges.push_back(even ? 50.0 + 12.5 * k : 50.0 * std::pow(1.07, k));
        axes.bands = {100.0, 200.0};
        std::vector<double> values;
        for (std::size_t b = 0; b < axes.bearings.size(); b++)
            for (std::size_t s = 0; s < axes.source_heights.size(); s++)
                for (std::size_t h = 0; h < axes.receiver_heights.size(); h++)
                    for (std::size_t r = 0; r < axes.ranges.size(); r++)
                        for (std::size_t f = 0; f < 2; f++)
                            values.push_back(std::sin(0.7 * static_cast<double>(b + 3 * s + 5 * h + 7 * r + f)));
        const PropagationTable table(axes, values, "test", false);
        const auto bands = table.match_bands({100.0, 200.0}, 0, 2);
        auto at_cell = [&](std::size_t b, std::size_t s, std::size_t h, std::size_t r, std::size_t f) {
            return table.at(b % axes.bearings.size(), s, h, r, f);
        };
        // Queries: every knot of every axis, the midpoints, and a spread between.
        std::vector<double> sources, receivers, ranges;
        auto spread = [](const std::vector<double>& axis, std::vector<double>& out) {
            for (std::size_t k = 0; k < axis.size(); k++)
            {
                out.push_back(axis[k]);
                if (k + 1 < axis.size())
                {
                    out.push_back(0.5 * (axis[k] + axis[k + 1]));
                    out.push_back(axis[k] + 0.999999 * (axis[k + 1] - axis[k]));
                    out.push_back(axis[k] + 1e-7 * (axis[k + 1] - axis[k]));
                }
            }
        };
        spread(axes.source_heights, sources);
        spread(axes.receiver_heights, receivers);
        spread(axes.ranges, ranges);
        int checked = 0;
        for (const double bearing : {0.0, 17.5, 90.0, 359.0, -45.0, 725.0})
        {
            for (std::size_t i = 0; i < sources.size(); i += 3)
            {
                for (std::size_t j = 0; j < receivers.size(); j += 5)
                {
                    for (std::size_t k = 0; k < ranges.size(); k += 2)
                    {
                        std::size_t s0, h0, r0;
                        double fs, fh, fr;
                        REQUIRE(reference_cell(axes.source_heights, sources[i], s0, fs));
                        REQUIRE(reference_cell(axes.receiver_heights, receivers[j], h0, fh));
                        REQUIRE(reference_cell(axes.ranges, ranges[k], r0, fr));
                        double u = std::fmod(bearing, 360.0);
                        if (u < 0.0)
                            u += 360.0;
                        u /= 30.0;
                        const std::size_t b0 = std::min(static_cast<std::size_t>(u), std::size_t(11));
                        const double fb = u - static_cast<double>(b0);
                        double expected[2] = {0.0, 0.0};
                        for (std::size_t f = 0; f < 2; f++)
                            for (int cb = 0; cb < 2; cb++)
                                for (int cs = 0; cs < 2; cs++)
                                    for (int ch = 0; ch < 2; ch++)
                                        for (int cr = 0; cr < 2; cr++)
                                        {
                                            const double weight = (cb ? fb : 1.0 - fb) * (cs ? fs : 1.0 - fs)
                                                    * (ch ? fh : 1.0 - fh) * (cr ? fr : 1.0 - fr);
                                            if (weight != 0.0)
                                                expected[f] += weight * at_cell(b0 + cb, s0 + cs, h0 + ch, r0 + cr, f);
                                        }
                        double out[2] = {};
                        INFO("even " << even << " bearing " << bearing << " source " << sources[i] << " receiver "
                                     << receivers[j] << " range " << ranges[k]);
                        REQUIRE(table.lookup(bearing, sources[i], receivers[j], ranges[k], bands, out));
                        REQUIRE(out[0] == Catch::Approx(expected[0]).margin(1e-12));
                        REQUIRE(out[1] == Catch::Approx(expected[1]).margin(1e-12));
                        checked++;
                    }
                }
            }
        }
        CHECK(checked > 1000);
    }
}

TEST_CASE("write a bump table for the windowed-versus-marched CTests", "[.bump_table]")
{
    // Not a test: CMakeLists.txt's propagation_table_bump_build runs it (by its hidden tag) to
    // make tables/bump_test.nc from tables/ground_model_test.nc, which the CTests build first.
    // The ground model's table never exceeds +6 dB and has one bearing, so the peak search's
    // use of the table's largest gain, and the run's bearing mapping, went untested: a bound
    // that left the table's gain out differed from the marched run by 0.57 dB in 21 of 23
    // cells on a table with +25 dB beyond 1500 m, and passed every CTest.  Here: +25 dB
    // beyond 1500 m and a 6 dB cos(bearing) tilt, on eight bearings.
    const auto base = PropagationTable::load("tables/ground_model_test.nc");
    PropagationTable::Axes axes = base.axes();
    REQUIRE(axes.bearings.size() == 1);
    axes.bearings.clear();
    for (int b = 0; b < 8; b++)
        axes.bearings.push_back(45.0 * b);
    std::vector<double> values;
    for (const double bearing : axes.bearings)
        for (std::size_t s = 0; s < axes.source_heights.size(); s++)
            for (std::size_t h = 0; h < axes.receiver_heights.size(); h++)
                for (std::size_t r = 0; r < axes.ranges.size(); r++)
                    for (std::size_t f = 0; f < axes.bands.size(); f++)
                    {
                        double v = base.at(0, s, h, r, f);
                        if (std::isfinite(v))
                            v += (axes.ranges[r] > 1500.0 ? 25.0 : 0.0) + 6.0 * std::cos(bearing * M_PI / 180.0);
                        values.push_back(v);
                    }
    PropagationTable(axes, values, "bump", false, {{"note", "bump_test"}}).save("tables/bump_test.nc");
    // The same values again, saying the table carries absorption: the run then adds none of its own.
    PropagationTable(axes, values, "bump", true, {{"note", "bump_test"}}).save("tables/bump_absorbed_test.nc");
}
