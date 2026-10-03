// New tests (no D original): the WGS84 projection ported from the D
// worldtransforms.d / wgs84tocartesian.d pair.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "geodetic.hpp"
#include "util/constants.hpp"
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
        REQUIRE(std::fabs(rotated.horizontal_norm() - east_north.horizontal_norm()) < 1e-6);
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
        const double projected = (b - a).horizontal_norm();
        const double from_file = std::hypot((current[5] - previous[5]) * ft2m,
                (current[6] - previous[6]) * ft2m);
        worst = std::max(worst, std::fabs(projected - from_file));
        longest = std::max(longest, from_file);
    }
    REQUIRE(longest > 100.0);
    REQUIRE(worst < 0.5);
}

TEST_CASE("The shipped tracks' frame runs along a bearing of 210 degrees", "[geodetic]")
{
    // The test above compares distances between samples, which no rotation or
    // translation can change -- so it says nothing about which way the axes
    // point.  The orientation is recoverable too, and it is a round number:
    // +x runs along a compass bearing of 210, the approach heading, not along
    // the instantaneous heading in column 8, which wanders over 25 degrees.
    // Match that and the only freedom left is the origin the file does not
    // record, which is one constant offset for the whole track.
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    const double ft2m = 0.3048;
    const Geodetic reference {records.front()[2], records.front()[3], 0.0};

    auto worst_residual = [&](double bearing) {
        const auto frame = LocalFrame::track_aligned(reference, bearing);
        // to_local is zero at the reference, so the offset is the file's own
        // coordinate for its first sample.
        const double offset_x = records.front()[5] * ft2m;
        const double offset_y = records.front()[6] * ft2m;
        double worst = 0.0;
        for (std::size_t i = 0; i < records.size(); i += 250)
        {
            const auto local = frame.to_local({records[i][2], records[i][3], 0.0});
            worst = std::max(worst, std::hypot(local(0) + offset_x - records[i][5] * ft2m,
                                            local(1) + offset_y - records[i][6] * ft2m));
        }
        return worst;
    };

    // Sub-metre over 3.2 km of track: this projection is the one that file was
    // built with.
    REQUIRE(worst_residual(210.0) < 1.0);
    // And the fit is tight enough to pin the bearing rather than merely allow
    // it -- five degrees out is hundreds of metres by the end.
    REQUIRE(worst_residual(215.0) > 100.0);
    REQUIRE(worst_residual(205.0) > 100.0);
}

TEST_CASE("An east/north/up vector takes the frame's rotation only", "[geodetic]")
{
    // A velocity has no origin to be projected about, so it rotates with the
    // axes but is not translated.  Comparing it against the projection of a
    // small displacement is what says both use the same rotation -- a sign slip
    // in one of them would put a trajectory's velocity across its own path.
    const Geodetic reference {40.0, -75.0, 0.0};
    for (const double bearing : {0.0, 37.0, 90.0, 210.0, 355.0})
    {
        const auto frame = LocalFrame::track_aligned(reference, bearing);
        const auto north_of_it = frame.to_local({40.001, -75.0, 0.0});
        // 0.001 degrees of latitude is about 111 m due north.
        const auto rotated = frame.from_east_north_up(Vec3(0.0, north_of_it.horizontal_norm(), 0.0));
        REQUIRE(std::fabs(rotated(0) - north_of_it(0)) < 1e-3);
        REQUIRE(std::fabs(rotated(1) - north_of_it(1)) < 1e-3);
        // Lengths survive a rotation, and the vertical is untouched.
        const Vec3 velocity(12.0, -34.0, 5.0);
        const auto turned = frame.from_east_north_up(velocity);
        REQUIRE(std::fabs(turned.horizontal_norm() - velocity.horizontal_norm()) < 1e-9);
        REQUIRE(turned(2) == velocity(2));
    }
    // An unrotated frame is east/north, so it leaves an ENU vector alone.
    const auto plain = LocalFrame(reference).from_east_north_up(Vec3(3.0, 4.0, 5.0));
    REQUIRE(plain(0) == Catch::Approx(3.0));
    REQUIRE(plain(1) == Catch::Approx(4.0));
    REQUIRE(plain(2) == Catch::Approx(5.0));
}

TEST_CASE("Longitude wraps across the antimeridian", "[geodetic]")
{
    // The difference of two longitudes was projected as it stood, so a reference at
    // 179.99 E and a point at 179.995 W -- 0.015 degrees apart -- came out 32,000
    // km away.  Over the shipped tracks' range nothing changes (remainder leaves a
    // difference within half a turn exactly as it is).
    const LocalFrame frame(Geodetic {10.0, 179.99, 0.0});
    const double metres_per_degree = 111319.49 * std::cos(10.0 * pi / 180.0);
    for (double longitude : {179.995, -179.995, -179.99, 179.97, -180.0, 180.0})
    {
        double expected = longitude - 179.99;
        if (expected < -180.0)
            expected += 360.0;
        const Vec3 local = frame.to_local(Geodetic {10.0, longitude, 0.0});
        INFO("longitude " << longitude);
        REQUIRE(local.x == Catch::Approx(expected * metres_per_degree).epsilon(0.01));
        REQUIRE(std::fabs(local.y) < 30.0);
    }
    // And back.
    const Vec3 across = frame.to_local(Geodetic {10.0, -179.995, 0.0});
    const Geodetic back = frame.to_geodetic(across);
    REQUIRE(std::remainder(back.longitude - (-179.995), 360.0) == Catch::Approx(0.0).margin(1e-7));
}

TEST_CASE("A position nowhere near the frame's origin is refused", "[geodetic]")
{
    // A sign-flipped --origin, or a longitude with its sign dropped, put the track
    // thousands of kilometres from a frame made for a few kilometres, and it was
    // projected anyway into a plausible-looking footprint for nowhere.
    const LocalFrame frame(Geodetic {39.8, -74.4, 0.0});
    REQUIRE_NOTHROW(frame.to_local(Geodetic {39.9, -74.3, 0.0}));
    REQUIRE_THROWS_AS(frame.to_local(Geodetic {39.8, 74.4, 0.0}), std::runtime_error);
    REQUIRE_THROWS_AS(frame.to_local(Geodetic {-39.8, -74.4, 0.0}), std::runtime_error);
    REQUIRE_THROWS_AS(frame.to_local(Geodetic {49.8, -74.4, 0.0}), std::runtime_error);
}
