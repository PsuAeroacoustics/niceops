// New tests (no D original): the WGS84 projection ported from the D
// worldtransforms.d / wgs84tocartesian.d pair.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "geodetic.hpp"
#include "util/csv.hpp"

using namespace niceops;

TEST_CASE("Projection is centred on its reference", "[geodetic]")
{
    const Geodetic reference {39.82914323, -74.40996141, 0.0};
    const LocalFrame frame(reference);

    const auto origin = frame.to_local(reference);
    REQUIRE(std::fabs(origin(0)) < 1e-6);
    REQUIRE(std::fabs(origin(1)) < 1e-6);

    // Altitude passes through untouched.
    Geodetic aloft = reference;
    aloft.altitude = 345.0;
    REQUIRE(frame.to_local(aloft)(2) == 345.0);
}

TEST_CASE("Projection has the right scale and orientation", "[geodetic]")
{
    const Geodetic reference {40.0, -75.0, 0.0};
    const LocalFrame frame(reference);

    // A tenth of a degree north is one meridian arc minute x 6: about 11.1 km,
    // and it must land on +y with no east component.
    const auto north = frame.to_local({40.1, -75.0, 0.0});
    REQUIRE(std::fabs(north(0)) < 1e-6);
    REQUIRE(std::fabs(north(1) - 11119.5) < 20.0);

    // A tenth of a degree east shrinks by cos(latitude): 11.1 km x cos(40).
    const auto east = frame.to_local({40.0, -74.9, 0.0});
    REQUIRE(east(0) > 0.0);
    REQUIRE(std::fabs(east(0) - 11119.5 * std::cos(40.0 * 3.14159265358979323846 / 180.0)) < 30.0);

    // South and west are the negatives, to the metre.
    const auto south = frame.to_local({39.9, -75.0, 0.0});
    REQUIRE(std::fabs(south(1) + north(1)) < 1.0);
    const auto west = frame.to_local({40.0, -75.1, 0.0});
    REQUIRE(std::fabs(west(0) + east(0)) < 1.0);
}

TEST_CASE("Round trip returns the original position", "[geodetic]")
{
    const LocalFrame frame({39.8, -74.4, 0.0});
    // Out to roughly 50 km in each direction, well beyond any footprint.
    for (double dlat = -0.5; dlat <= 0.5; dlat += 0.25)
    {
        for (double dlon = -0.5; dlon <= 0.5; dlon += 0.25)
        {
            const Geodetic original {39.8 + dlat, -74.4 + dlon, 1234.5};
            const auto local = frame.to_local(original);
            const auto recovered = frame.to_geodetic(local);
            // A micro-degree is about 0.1 mm; the inverse is far tighter than
            // the centimetre the D original's fixed-step search stopped at.
            REQUIRE(std::fabs(recovered.latitude - original.latitude) < 1e-9);
            REQUIRE(std::fabs(recovered.longitude - original.longitude) < 1e-9);
            REQUIRE(recovered.altitude == original.altitude);
        }
    }
}

TEST_CASE("Track-aligned frames rotate the axes", "[geodetic]")
{
    const Geodetic reference {40.0, -75.0, 0.0};
    const Geodetic north_of_it {40.1, -75.0, 0.0};

    // Bearing 90 -- +x due east -- is the plain east/north frame.
    const auto east_north = LocalFrame(reference).to_local(north_of_it);
    const auto bearing_90 = LocalFrame::track_aligned(reference, 90.0).to_local(north_of_it);
    REQUIRE(std::fabs(bearing_90(0) - east_north(0)) < 1e-6);
    REQUIRE(std::fabs(bearing_90(1) - east_north(1)) < 1e-6);

    // With +x pointing due north, a point to the north lands on +x.
    const auto bearing_0 = LocalFrame::track_aligned(reference, 0.0).to_local(north_of_it);
    REQUIRE(std::fabs(bearing_0(0) - east_north(1)) < 1e-6);
    REQUIRE(std::fabs(bearing_0(1)) < 1e-6);

    // Rotation preserves distance, whatever the bearing.
    for (double bearing = -180.0; bearing <= 180.0; bearing += 37.0)
    {
        const auto rotated = LocalFrame::track_aligned(reference, bearing).to_local(north_of_it);
        REQUIRE(std::fabs(rotated.head<2>().norm() - east_north.head<2>().norm()) < 1e-6);
    }
}

TEST_CASE("Projection matches the shipped tracks' own geometry", "[geodetic]")
{
    // The 27-column AC format carries lat/lon/alt (columns 2-4) *and* a local
    // x/y/z (columns 5-7) derived from them, which checks this port against
    // real data rather than against itself.
    //
    // That x/y is NOT east/north: it is rotated to run along the flight track
    // (by about the heading in column 8), and its origin is not recorded.  So
    // compare the quantity that survives both unknowns -- the distance between
    // samples, which no rotation or translation can change.
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(records.size() > 100);
    REQUIRE(records.front().size() == 27);

    const double ft2m = 0.3048;
    const LocalFrame frame(Geodetic {records.front()[2], records.front()[3], 0.0});

    double worst = 0.0;
    double longest = 0.0;
    // Stride well past the 1e-8 degree quantisation of the lat/lon columns, so
    // the comparison is of real displacements rather than of rounding.
    const std::size_t stride = 250;
    for (std::size_t i = stride; i < records.size(); i += stride)
    {
        const auto& previous = records[i - stride];
        const auto& current = records[i];
        const auto a = frame.to_local({previous[2], previous[3], 0.0});
        const auto b = frame.to_local({current[2], current[3], 0.0});
        const double projected = (b - a).head<2>().norm();
        const double from_file = std::hypot((current[5] - previous[5]) * ft2m,
                (current[6] - previous[6]) * ft2m);
        worst = std::max(worst, std::fabs(projected - from_file));
        longest = std::max(longest, from_file);
    }
    REQUIRE(longest > 100.0);
    REQUIRE(worst < 0.5);
}
