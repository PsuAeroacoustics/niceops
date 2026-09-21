// New test (no D original): exercise the Trajectory CSV/JSON parsing path
// with the demo case.  CTest runs with the source root as the working
// directory, so the default relative paths resolve here.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

#include "trajectory.hpp"
#include "util/csv.hpp"

using namespace niceops;

TEST_CASE("Trajectory parses the demo case", "[trajectory]")
{
    Trajectory trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json");
    REQUIRE(trajectory.states.size() > 0);
    for (const auto& state : trajectory.states)
    {
        REQUIRE(std::isfinite(state.CT));
        REQUIRE(std::isfinite(state.mu));
        REQUIRE(std::isfinite(state.fpa));
        REQUIRE(std::isfinite(state.fpa_e));
    }
    auto x = trajectory.column("x");
    REQUIRE(x.size() == trajectory.states.size());
}

TEST_CASE("Hover states do not produce NaN flight path angles", "[trajectory]")
{
    // A stationary aircraft makes v zero, and the flight path angles used to be
    // asin of a division by it -- NaN, which then travelled through mu into the
    // sphere lookup and out into every level.  The database carries hover
    // spheres at mu = 0, so this regime is meant to be reachable.
    dimensional_data vehicle;
    vehicle.tip_speed = 200.0;
    vehicle.radius = 6.7;
    vehicle.weight = 5000.0 * 9.81;
    vehicle.rho = 1.225;

    State hover;
    hover.vx = 0.0;
    hover.vy = 0.0;
    hover.vz = 0.0;
    hover.nx = 0.0;
    hover.ny = 0.0;
    hover.nz = 1.0;
    derive_flight_parameters(hover, vehicle);
    REQUIRE(std::isfinite(hover.fpa));
    REQUIRE(std::isfinite(hover.fpa_e));
    REQUIRE(hover.fpa == 0.0);
    REQUIRE(hover.fpa_e == 0.0);
    REQUIRE(hover.mu == 0.0);
    REQUIRE(std::isfinite(hover.CT));

    // Forward flight is untouched: a climb along +x with the rotor level has
    // its flight path angle set by the vertical component alone.
    State climbing;
    climbing.vx = 50.0;
    climbing.vy = 0.0;
    climbing.vz = 50.0;
    climbing.nx = 0.0;
    climbing.ny = 0.0;
    climbing.nz = 1.0;
    derive_flight_parameters(climbing, vehicle);
    REQUIRE(std::fabs(climbing.fpa - 45.0) < 1e-9);
}

// ---------------------------------------------------------------------------
// The labelled WGS84 input format.
// ---------------------------------------------------------------------------

namespace {

// Minimal scoped temp file, as tests/test_csv.cpp uses.  The name has to end
// in .csv: the trajectory reader picks its parser by extension first, and a
// suffixless temp file fails there rather than anywhere this is testing -- which
// would leave every REQUIRE_THROWS below passing for the wrong reason.
struct TempCsv
{
    std::string path;
    explicit TempCsv(const std::string& contents)
    {
        char pattern[] = "/tmp/niceops_geodetic_XXXXXX.csv";
        const int fd = ::mkstemps(pattern, 4);
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

// The demo track's own WGS84 columns (2-4), written out under `header`.
// `stride` thins it, for the sparse-input cases.
std::string demo_track_as_geodetic(const std::string& header, std::size_t stride = 1)
{
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(records.size() > 100);
    std::ostringstream out;
    out << header << "\n" << std::setprecision(12);
    for (std::size_t i = 0; i < records.size(); i += stride)
    {
        out << records[i][1] << "," << records[i][2] << "," << records[i][3] << ","
            << records[i][4] << "\n";
    }
    return out.str();
}

// The demo track's local frame: +x along a bearing of 210 degrees from its own
// first position, with the ground 54.2651 ft below the altitude column.  Both
// numbers are properties of that file, recovered from it rather than assumed --
// see "A WGS84 track reproduces the flight it was taken from" below.
GeodeticInput demo_frame(double smoothing = 0.0)
{
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    GeodeticInput input;
    input.frame = LocalFrame::track_aligned(
            Geodetic {records.front()[2], records.front()[3], 0.0}, 210.0);
    input.ground_elevation_ft = 54.2651;
    input.path_smoothing = smoothing;
    return input;
}

} // namespace

TEST_CASE("A WGS84 track reproduces the flight it was taken from", "[trajectory][geodetic]")
{
    // The 27-column AC file carries lat/lon/alt *and* the local x/y/z derived
    // from them, so feeding its geodetic columns back in checks the whole input
    // path -- projection, fit, differentiation -- against real data rather than
    // against itself.  Its x/y origin is not recorded in the file, so the
    // comparison is up to the one constant offset that cannot be recovered.
    const TempCsv track(demo_track_as_geodetic("time,lat,lon,alt"));
    const Trajectory geodetic(track.path, "databases/S-76D.json", demo_frame(0.02));
    const Trajectory reference("tracks/S-76D/197421AC.csv", "databases/S-76D.json");

    // 0.02 s samples, already uniform and finer than the resampling step, so
    // the file's own time base is kept rather than resampled onto another.
    REQUIRE(geodetic.states.size() == reference.states.size());

    const double offset_x = reference.states.front().x - geodetic.states.front().x;
    const double offset_y = reference.states.front().y - geodetic.states.front().y;
    // Roughly 2.5 km back along the track: the file's origin is nowhere near
    // its first sample, which is why it has to be divided out.
    REQUIRE(std::hypot(offset_x, offset_y) > 1000.0);

    double worst_position = 0.0;
    double worst_height = 0.0;
    double worst_mu = 0.0;
    for (std::size_t i = 0; i < geodetic.states.size(); i++)
    {
        const auto& got = geodetic.states[i];
        const auto& want = reference.states[i];
        worst_position = std::max(worst_position,
                std::hypot(got.x + offset_x - want.x, got.y + offset_y - want.y));
        worst_height = std::max(worst_height, std::fabs(got.z - want.z));
        worst_mu = std::max(worst_mu, std::fabs(got.mu - want.mu));
    }
    // Sub-metre over a 3.2 km track, which is the projection agreeing with
    // whatever the file's own producer used.
    REQUIRE(worst_position < 1.0);
    // Smoothing is applied to all three coordinates, so height is fitted like
    // the rest rather than copied; it still lands well inside a foot.
    REQUIRE(worst_height < 0.2);
    // Speed comes from differentiating the path instead of from the file's
    // inertial columns, so it agrees closely rather than exactly.
    REQUIRE(worst_mu < 0.05);

    // With nothing to smooth, altitude passes straight through: the column is
    // above mean sea level, z is above the observers, and the difference is the
    // ground elevation the run was given.
    const Trajectory exact(track.path, "databases/S-76D.json", demo_frame(0.0));
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    for (std::size_t i = 0; i < exact.states.size(); i += 250)
    {
        REQUIRE(exact.states[i].z
                == Catch::Approx((records[i][4] - 54.2651) * ft2m).epsilon(1e-12));
    }
}

TEST_CASE("Geodetic columns are found by name, not by position", "[trajectory][geodetic]")
{
    // A planner's export is not going to have the columns in our order, or
    // spell them the way we would, so they are matched by name after
    // normalisation -- and an unrecognised header is an error, not a guess.
    const TempCsv canonical(demo_track_as_geodetic("time,lat,lon,alt", 100));
    const TempCsv dressed(demo_track_as_geodetic("Time (s),Latitude,LONGITUDE,alt_ft", 100));
    const Trajectory plain(canonical.path, "databases/S-76D.json", demo_frame());
    const Trajectory fancy(dressed.path, "databases/S-76D.json", demo_frame());

    REQUIRE(plain.states.size() == fancy.states.size());
    for (std::size_t i = 0; i < plain.states.size(); i++)
    {
        REQUIRE(plain.states[i].x == fancy.states[i].x);
        REQUIRE(plain.states[i].y == fancy.states[i].y);
        REQUIRE(plain.states[i].z == fancy.states[i].z);
    }

    const TempCsv nameless("a,b,c,d\n0,1,2,3\n1,4,5,6\n2,7,8,9\n");
    REQUIRE_THROWS(Trajectory(nameless.path, "databases/S-76D.json", demo_frame()));
}

TEST_CASE("Sparse input is resampled and dense input is left alone", "[trajectory][geodetic]")
{
    GeodeticInput input = demo_frame();
    input.resample_dt = 0.1;

    // Every 500th sample: 10 s apart, which no sound exposure integral could be
    // taken over, so the path is refitted onto the 0.1 s grid.
    const TempCsv sparse(demo_track_as_geodetic("time,lat,lon,alt", 500));
    const Trajectory densified(sparse.path, "databases/S-76D.json", input);
    REQUIRE(densified.states.size() > 1000);
    const double step = densified.states[1].t - densified.states[0].t;
    REQUIRE(std::fabs(step - 0.1) < 1e-9);

    // The full 50 Hz track is already finer than the step asked for, and
    // resampling it would only throw samples away.
    const TempCsv dense(demo_track_as_geodetic("time,lat,lon,alt"));
    const Trajectory kept(dense.path, "databases/S-76D.json", input);
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(kept.states.size() == records.size());
    REQUIRE(std::fabs((kept.states[1].t - kept.states[0].t) - 0.02) < 1e-9);
}

TEST_CASE("Velocity columns are used when the file carries them", "[trajectory][geodetic]")
{
    // East/north/up in ft/s, which the frame rotates onto its own axes.  The
    // demo track's velocity columns are already in that frame, so rotating them
    // back to east/north here and letting the reader rotate them forward again
    // has to return exactly what the file started with.
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    const double rotation = (90.0 - 210.0) * 3.14159265358979323846 / 180.0;
    const double cos_r = std::cos(rotation);
    const double sin_r = std::sin(rotation);

    std::ostringstream out;
    out << "time,lat,lon,alt,vE,vN,vUp\n" << std::setprecision(17);
    for (std::size_t i = 0; i < records.size(); i += 10)
    {
        const double vx = records[i][11];
        const double vy = records[i][13];
        out << records[i][1] << "," << records[i][2] << "," << records[i][3] << ","
            << records[i][4] << "," << vx * cos_r - vy * sin_r << "," << vx * sin_r + vy * cos_r
            << "," << records[i][15] << "\n";
    }
    const TempCsv track(out.str());
    GeodeticInput input = demo_frame();
    input.resample_dt = 0.2; // the thinned track's own step, so nothing is refitted
    const Trajectory trajectory(track.path, "databases/S-76D.json", input);

    REQUIRE(trajectory.states.size() == (records.size() + 9) / 10);
    double worst = 0.0;
    for (std::size_t i = 0; i < trajectory.states.size(); i++)
    {
        const auto& state = trajectory.states[i];
        worst = std::max({worst, std::fabs(state.vx - records[i * 10][11] * ft2m),
                std::fabs(state.vy - records[i * 10][13] * ft2m),
                std::fabs(state.vz - records[i * 10][15] * ft2m)});
    }
    // A rotation there and back, so this is round-off and nothing else.
    REQUIRE(worst < 1e-8);
}

TEST_CASE("Unusable geodetic input is rejected where it can be named", "[trajectory][geodetic]")
{
    // No time column: there is nothing to differentiate against.
    const TempCsv timeless("lat,lon,alt\n39.8,-74.4,100\n39.81,-74.41,200\n39.82,-74.42,300\n");
    REQUIRE_THROWS(Trajectory(timeless.path, "databases/S-76D.json", demo_frame()));

    // Degrees, so a latitude past the pole is a unit error rather than a place.
    const TempCsv off_planet(
            "time,lat,lon,alt\n0,139.8,-74.4,100\n1,39.81,-74.41,200\n2,39.82,-74.42,300\n");
    REQUIRE_THROWS(Trajectory(off_planet.path, "databases/S-76D.json", demo_frame()));

    // Two points fix a straight line and no curvature, so no load factor.
    const TempCsv too_short("time,lat,lon,alt\n0,39.8,-74.4,100\n1,39.81,-74.41,200\n");
    REQUIRE_THROWS(Trajectory(too_short.path, "databases/S-76D.json", demo_frame()));

    // Time has to advance, or the fit is asked to be two things at one instant.
    const TempCsv backwards(
            "time,lat,lon,alt\n0,39.8,-74.4,100\n2,39.81,-74.41,200\n1,39.82,-74.42,300\n");
    REQUIRE_THROWS(Trajectory(backwards.path, "databases/S-76D.json", demo_frame()));
}

TEST_CASE("The AC state format is still read positionally", "[trajectory][geodetic]")
{
    // The AC format's own header names lat, lon and alt columns, so a reader
    // that dispatched on the header alone would take a state file for a WGS84
    // track -- and get its position from columns whose x/y origin the file does
    // not record.  Width decides first, and first_geodetic_position agrees.
    REQUIRE_FALSE(first_geodetic_position("tracks/S-76D/197421AC.csv").has_value());

    const TempCsv track(demo_track_as_geodetic("time,lat,lon,alt", 500));
    const auto first = first_geodetic_position(track.path);
    REQUIRE(first.has_value());
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(first->latitude == Catch::Approx(records.front()[2]));
    REQUIRE(first->longitude == Catch::Approx(records.front()[3]));
}
