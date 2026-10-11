// New test (no D original): exercise the Trajectory CSV/JSON parsing path
// with the demo case.  CTest runs with the source root as the working
// directory, so the default relative paths resolve here.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

#include "atmosphere_profile.hpp"
#include "sphere_database.hpp"
#include "trajectory.hpp"
#include "util/csv.hpp"
#include "util/paths.hpp"

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
    // asin of a division by it -- NaN, which then traveled through mu into the
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
// The labeled WGS84 input format.
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
    input.ground_elevation = 54.2651 * ft2m;
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
    // Sub-meter over a 3.2 km track, which is the projection agreeing with
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
    // normalization -- and an unrecognized header is an error, not a guess.
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
    REQUIRE_FALSE(first_geodetic_position("tracks/S-76D/197421AC.csv", feet).has_value());

    const TempCsv track(demo_track_as_geodetic("time,lat,lon,alt", 500));
    const auto first = first_geodetic_position(track.path, feet);
    REQUIRE(first.has_value());
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(first->latitude == Catch::Approx(records.front()[2]));
    REQUIRE(first->longitude == Catch::Approx(records.front()[3]));
}

TEST_CASE("The run frame is found for a metric geodetic track", "[trajectory][geodetic][units]")
{
    // run() settles the frame by reading the first row before it loads anything.
    // It once did so in feet whatever --units said, so a track headed alt_m (a
    // spelling that exists only in meters) looked unpositioned: no frame, and
    // --frame_bearing silently ignored.
    const TempCsv metric(demo_track_as_geodetic("time,lat,lon,alt_m", 500));
    REQUIRE_FALSE(first_geodetic_position(metric.path, feet).has_value());
    const auto first = first_geodetic_position(metric.path, meters);
    REQUIRE(first.has_value());
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(first->latitude == Catch::Approx(records.front()[2]));
    REQUIRE(first->longitude == Catch::Approx(records.front()[3]));
}

TEST_CASE("adapt resamples onto a uniform grid and keeps the clock", "[trajectory][adapt]")
{
    // adapt() is the only smoothing/resampling path for the positional AC
    // format (--path_smoothing reaches the geodetic reader only), and nothing
    // in main calls it, so this is what keeps it honest.
    Trajectory trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json");

    // A straight line at constant speed, which a natural cubic spline
    // reproduces exactly: f'' = 0 everywhere already satisfies the natural end
    // condition, so there is no boundary artifact to allow for.
    const double speed_x = 3.0;
    const double speed_y = -4.0;
    const double start_utc = 43200.0; // noon
    trajectory.states.clear();
    for (int i = 0; i <= 100; i++)
    {
        State now {};
        now.t = 0.02 * i;
        now.utcsec = start_utc + now.t;
        now.x = 100.0 + speed_x * now.t;
        now.y = speed_y * now.t;
        now.z = 250.0;
        trajectory.states.push_back(now);
    }

    trajectory.adapt(0.1, /*smoothing=*/0.0);

    // 2 s of track at 0.1 s: 21 samples, endpoint included.
    REQUIRE(trajectory.states.size() == 21);
    for (std::size_t i = 0; i < trajectory.states.size(); i++)
    {
        const auto& now = trajectory.states[i];
        REQUIRE(now.t == Catch::Approx(0.1 * static_cast<double>(i)));
        // UTC seconds used to be dropped here: the repacked State was
        // value-initialized and nothing refilled utcsec, which moved the whole
        // flight to midnight and with it --dnlfile's day/night split.
        REQUIRE(now.utcsec == Catch::Approx(start_utc + now.t));
        REQUIRE(now.x == Catch::Approx(100.0 + speed_x * now.t));
        REQUIRE(now.vx == Catch::Approx(speed_x));
        REQUIRE(now.vy == Catch::Approx(speed_y));
        REQUIRE(now.vz == Catch::Approx(0.0).margin(1e-9));
        // Straight and level: the load factor is gravity alone.
        REQUIRE(now.LF == Catch::Approx(1.0));
        REQUIRE(std::isfinite(now.CT));
        REQUIRE(std::isfinite(now.mu));
    }
}

TEST_CASE("adapt takes the load factor from the fitted path", "[trajectory][adapt]")
{
    // The load factor is the second derivative, and it is read off the spline
    // rather than finite-differenced from the resampled positions.  A constant
    // acceleration is the case that says so: the interior of an interpolating
    // cubic through a parabola carries it exactly.
    Trajectory trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json");

    const double acceleration = 2.0; // m/s^2 along x
    trajectory.states.clear();
    for (int i = 0; i <= 200; i++)
    {
        State now {};
        now.t = 0.02 * i;
        now.x = 0.5 * acceleration * now.t * now.t;
        now.y = 0.0;
        now.z = 300.0;
        trajectory.states.push_back(now);
    }

    trajectory.adapt(0.1, /*smoothing=*/0.0);

    // Skip the two end cells, where the natural boundary condition pins the
    // curvature to zero and so cannot represent a parabola.
    const std::size_t count = trajectory.states.size();
    REQUIRE(count > 8);
    for (std::size_t i = 3; i + 3 < count; i++)
    {
        const auto& now = trajectory.states[i];
        REQUIRE(now.vx == Catch::Approx(acceleration * now.t).epsilon(1e-6));
        // nx is normalized by LF, so recover the acceleration in g from both.
        REQUIRE(now.nx * now.LF == Catch::Approx(acceleration * mps2g).epsilon(1e-6));
        REQUIRE(now.nz * now.LF == Catch::Approx(1.0).epsilon(1e-6));
    }
}

TEST_CASE("the default observer grid is bit-exact in feet", "[units]")
{
    // Every length on the command line is in feet, the grid included.  Its
    // default is +/-1500 m in steps of 25 m spelled in feet, and the goldens
    // encode the resulting observer positions -- so the two spellings have to
    // generate the identical sequence, not merely endpoints that agree.
    const double bound_ft = 4921.2598425196848;
    const double step_ft = 82.020997375328079;

    // The count comes from ceil((end - begin) / step); one ULP the wrong way
    // would give 121 points instead of 120 and move every golden.
    const double count = std::ceil((bound_ft - -bound_ft) / step_ft);
    REQUIRE(count == 120.0);

    for (int i = 0; i < 120; i++)
    {
        const double as_meters = -1500.0 + 25.0 * i;
        const double via_feet = (-bound_ft + step_ft * i) * ft2m;
        REQUIRE(via_feet == as_meters);
    }
}

TEST_CASE("path_smoothing bounds the deviation, in feet", "[units][geodetic]")
{
    // It is an acceptable deviation per axis, not an abstract weight: the fit
    // is allowed to sit that far from the given positions and will use the
    // whole allowance.  Ten feet of slack has to show up as ten feet.
    const TempCsv track(demo_track_as_geodetic("time,lat,lon,alt"));

    auto rms_offset = [&](double smoothing_ft) {
        auto input = demo_frame(smoothing_ft * ft2m);
        Trajectory fitted(track.path, "databases/S-76D.json", input);
        auto exact_input = demo_frame(0.0);
        Trajectory exact(track.path, "databases/S-76D.json", exact_input);
        REQUIRE(fitted.states.size() == exact.states.size());
        double sum = 0.0;
        for (std::size_t i = 0; i < fitted.states.size(); i++)
        {
            const double dx = fitted.states[i].x - exact.states[i].x;
            sum += dx * dx;
        }
        return std::sqrt(sum / static_cast<double>(fitted.states.size())) / ft2m;
    };

    // The exact fit passes through the given points, so the offset from it is
    // the deviation the bound allows.  Loose tolerance: the bound is on the
    // residual against the data, and x is one axis of three.
    REQUIRE(rms_offset(10.0) > 1.0);
    REQUIRE(rms_offset(10.0) < 20.0);
    // Monotone in the allowance.
    REQUIRE(rms_offset(40.0) > rms_offset(10.0));
}

TEST_CASE("--units m reads the same flight from a metric file", "[units][geodetic]")
{
    // The unit system is a statement about the interface, not about the
    // physics: the same flight written in meters and read with meters has to
    // give the same states as the feet original read as feet.
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    REQUIRE(records.size() > 100);

    std::ostringstream metric;
    metric << std::setprecision(17) << "utcsec,t,lat,lon,alt,x,y,z,h,i,j,vx,l,vy,n,vz,p,q,r,s,t2,"
                                       "ax,ay,az,x1,y1,z1\n";
    for (const auto& r : records)
    {
        REQUIRE(r.size() == 27);
        auto row = r;
        for (std::size_t c : {5u, 6u, 7u, 11u, 13u, 15u, 21u, 22u, 23u})
            row[c] *= ft2m;
        for (std::size_t c = 0; c < row.size(); c++)
            metric << (c ? "," : "") << row[c];
        metric << "\n";
    }
    const TempCsv metric_track(metric.str());

    const Trajectory in_feet("tracks/S-76D/197421AC.csv", "databases/S-76D.json", {}, feet);
    const Trajectory in_meters(metric_track.path, "databases/S-76D.json", {}, meters);

    REQUIRE(in_meters.states.size() == in_feet.states.size());
    for (std::size_t i = 0; i < in_feet.states.size(); i++)
    {
        const auto& a = in_feet.states[i];
        const auto& b = in_meters.states[i];
        // Positions and velocities are already metric internally, so they match
        // to the round-off of writing the file out and reading it back.
        REQUIRE(std::fabs(b.x - a.x) < 1e-9);
        REQUIRE(std::fabs(b.z - a.z) < 1e-9);
        REQUIRE(std::fabs(b.vx - a.vx) < 1e-9);
        // The derived flight parameters are what actually index the database.
        REQUIRE(std::fabs(b.CT - a.CT) < 1e-12);
        REQUIRE(std::fabs(b.mu - a.mu) < 1e-12);
    }
}

TEST_CASE("an altitude column may only claim the unit in force", "[units][geodetic]")
{
    // "alt_ft" under --units m would read feet as meters and fly the aircraft
    // a third as high; refusing to match is what makes that loud instead of
    // silent.  A bare "alt" claims nothing and is always taken.
    const std::string body = "0,39.81,-74.43,300\n1,39.811,-74.431,299\n2,39.812,-74.432,298\n";
    const TempCsv in_feet("time,lat,lon,alt_ft\n" + body);
    const TempCsv in_meters("time,lat,lon,alt_m\n" + body);
    const TempCsv unlabeled("time,lat,lon,alt\n" + body);

    GeodeticInput input;
    input.frame = LocalFrame(Geodetic {39.81, -74.43, 0.0});

    REQUIRE_NOTHROW(Trajectory(in_feet.path, "databases/S-76D.json", input, feet));
    REQUIRE_THROWS(Trajectory(in_feet.path, "databases/S-76D.json", input, meters));
    REQUIRE_NOTHROW(Trajectory(in_meters.path, "databases/S-76D.json", input, meters));
    REQUIRE_THROWS(Trajectory(in_meters.path, "databases/S-76D.json", input, feet));
    REQUIRE_NOTHROW(Trajectory(unlabeled.path, "databases/S-76D.json", input, feet));
    REQUIRE_NOTHROW(Trajectory(unlabeled.path, "databases/S-76D.json", input, meters));
}

namespace {
// A steady 6 degree descent at 36 m/s, 0.02 s apart, with vz stored as the
// file would: positive up (vz_positive_up) or positive down.  derive has
// already normalized the load factor, as the reader leaves it.
std::vector<State> steady_descent(bool vz_positive_up)
{
    dimensional_data vehicle;
    vehicle.tip_speed = 230.0;
    vehicle.radius = 5.3;
    vehicle.weight = 2250.0 * 9.81;
    vehicle.rho = 1.2;
    const double v = 36.0, gamma = 6.0 * pi / 180.0;
    std::vector<State> states;
    for (int i = 0; i < 500; i++)
    {
        State now {};
        now.t = 0.02 * i;
        now.x = v * std::cos(gamma) * now.t;
        now.z = 500.0 - v * std::sin(gamma) * now.t;
        now.vx = v * std::cos(gamma);
        now.vz = (vz_positive_up ? -1.0 : 1.0) * v * std::sin(gamma);
        now.nz = 1.0;
        derive_flight_parameters(now, vehicle);
        states.push_back(now);
    }
    return states;
}
} // namespace

TEST_CASE("vz stored positive down is detected and made positive up", "[trajectory][vz]")
{
    // The 2017 Noise Abatement AC files (all but EC130B4 day 298's) store vz
    // positive down; read as up, every descent was flown as a climb.
    auto states = steady_descent(false);
    REQUIRE(states[10].fpa > 0.0); // as read: a climb
    REQUIRE(correct_vertical_speed_sign(states, "descent_vz_down"));
    for (const auto& s : states)
    {
        REQUIRE(s.vz < 0.0);
        REQUIRE(std::fabs(s.fpa - (-6.0)) < 1e-9);
    }
}

TEST_CASE("vz already positive up is left alone, even at a constant rate", "[trajectory][vz]")
{
    // A steady descent has a constant dz/dt, so a centered correlation has no
    // variance to work with; the check has to decide from the uncentered one.
    auto states = steady_descent(true);
    const auto before = states;
    REQUIRE_FALSE(correct_vertical_speed_sign(states, "descent_vz_up"));
    for (std::size_t i = 0; i < states.size(); i++)
    {
        REQUIRE(states[i].vz == before[i].vz);
        REQUIRE(states[i].fpa == before[i].fpa);
    }
}

TEST_CASE("a level track gives no evidence and is left alone", "[trajectory][vz]")
{
    auto states = steady_descent(true);
    for (auto& s : states)
    {
        s.z = 300.0;
        s.vz = 0.0;
    }
    REQUIRE_FALSE(correct_vertical_speed_sign(states, "level"));
}

TEST_CASE("the shipped demo tracks store vz positive up and are read unchanged", "[trajectory][vz]")
{
    for (const auto& [track, vehicle] :
            {std::pair {"tracks/S-76D/197421AC.csv", "databases/S-76D.json"},
                    std::pair {"tracks/B205/190337AC.csv", "databases/B205.json"}})
    {
        Trajectory trajectory(track, vehicle);
        std::vector<State> copy = trajectory.states;
        REQUIRE_FALSE(correct_vertical_speed_sign(copy, track));
    }
}

// --speed_reference air: mu, the path angles and the sphere's nose follow the
// velocity through the air, the ground velocity less the wind at the state.
namespace {

std::function<std::array<double, 2>(double)> uniform_wind(double wx, double wy)
{
    return [=](double) { return std::array<double, 2> {wx, wy}; };
}

// One state flying +x at `speed` m/s, TPP level, taken from the demo track so
// that its CT and the vehicle's tip speed are real ones.
Trajectory single_state(double speed, double vz, double z = 100.0)
{
    Trajectory trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json");
    State s = trajectory.states.front();
    s.x = 0.0;
    s.y = 0.0;
    s.z = z;
    s.vx = speed;
    s.vy = 0.0;
    s.vz = vz;
    s.nx = 0.0;
    s.ny = 0.0;
    s.nz = 1.0;
    trajectory.states = {s};
    // mu and the path angles for that velocity, over the ground.
    trajectory.relative_to_air(uniform_wind(0.0, 0.0));
    return trajectory;
}

} // namespace

TEST_CASE("an air reference in calm air changes nothing, to the bit", "[trajectory][speed_reference]")
{
    const Trajectory ground("tracks/S-76D/197421AC.csv", "databases/S-76D.json");
    Trajectory air = ground;
    REQUIRE(air.relative_to_air(uniform_wind(0.0, 0.0)) == 0.0);
    for (std::size_t i = 0; i < ground.states.size(); i++)
    {
        const auto& g = ground.states[i];
        const auto& a = air.states[i];
        REQUIRE(a.mu == g.mu);
        REQUIRE(a.fpa == g.fpa);
        REQUIRE(a.fpa_e == g.fpa_e);
        REQUIRE(a.CT == g.CT);
    }
}

TEST_CASE("a headwind raises mu and flattens the path; CT stays inertial", "[trajectory][speed_reference]")
{
    auto trajectory = single_state(40.0, -3.0);
    const State before = trajectory.states.front();
    // Blowing toward -x, into the aircraft's face, at 10 m/s.
    REQUIRE(trajectory.relative_to_air(uniform_wind(-10.0, 0.0)) == 10.0);
    const State& after = trajectory.states.front();
    const double airspeed = std::hypot(50.0, 3.0);
    REQUIRE(after.mu == Catch::Approx(airspeed / trajectory.vehicle.tip_speed));
    REQUIRE(after.mu > before.mu);
    REQUIRE(after.fpa == Catch::Approx(std::asin(-3.0 / airspeed) * 180.0 / pi));
    REQUIRE(std::fabs(after.fpa) < std::fabs(before.fpa));
    // The TPP is level, so the QSAM angle is the path angle.
    REQUIRE(after.fpa_e == Catch::Approx(after.fpa));
    // The thrust balances the acceleration over the ground, and the ground
    // velocity itself is the track's.
    REQUIRE(after.CT == before.CT);
    REQUIRE(after.nz == before.nz);
    REQUIRE(after.vx == 40.0);
}

TEST_CASE("a crosswind turns the sphere's nose into the relative wind", "[trajectory][speed_reference]")
{
    auto trajectory = single_state(40.0, 0.0);
    // Blowing toward +y: the aircraft has to point toward -y to hold its track.
    trajectory.relative_to_air(uniform_wind(0.0, 5.0));
    const auto v = air_velocity(trajectory.states.front());
    REQUIRE(v == std::array<double, 3> {40.0, -5.0, 0.0});
    const auto frame = sphere_frame(Vec3(v[0], v[1], v[2]), Vec3(0.0, 0.0, 1.0));
    REQUIRE(std::atan2(frame.nose.y, frame.nose.x) == Catch::Approx(std::atan2(-5.0, 40.0)));
    REQUIRE(trajectory.states.front().mu * trajectory.vehicle.tip_speed
            == Catch::Approx(std::hypot(40.0, 5.0)));
}

TEST_CASE("the wind is the profile's at each state's height, in the track frame",
        "[trajectory][speed_reference]")
{
    // From the north, calm at the ground and 10 m/s at 100 m.
    AtmosphereProfile profile;
    profile.height = {0.0, 100.0};
    profile.temperature = {288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0};
    profile.wind_east = {0.0, 0.0};
    profile.wind_north = {0.0, -10.0};
    // The 2017 tracks' frame: +x along 270 (west), so +y along 180 (south),
    // which is where a north wind blows.
    const CompassToLocal to_local(270.0);
    auto trajectory = single_state(40.0, 0.0, 50.0);
    State high = trajectory.states.front();
    high.z = 100.0;
    trajectory.states.push_back(high);
    const double strongest = trajectory.relative_to_air(
            [&](double z) { return to_local(profile.wind_at(z)); });
    REQUIRE(strongest == Catch::Approx(10.0));
    REQUIRE(trajectory.states[0].wx == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(trajectory.states[0].wy == Catch::Approx(5.0));
    REQUIRE(trajectory.states[1].wy == Catch::Approx(10.0));
    REQUIRE(trajectory.states[1].mu * trajectory.vehicle.tip_speed
            == Catch::Approx(std::hypot(40.0, 10.0)));
}

TEST_CASE("at hover the sphere's nose is the heading, turned into the local frame", "[trajectory]")
{
    // A local frame with +x along compass 270 and +y along 180, as the 2017 AC
    // tracks are.  Fast states fly west on heading 274 (an INS heading 4 degrees
    // off the ground track), so the offset inferred from them folds that in.
    auto state = [](double vx, double vy, double heading) {
        State s {};
        s.vx = vx;
        s.vy = vy;
        s.heading = heading;
        return s;
    };
    std::vector<State> states;
    for (int i = 0; i < 60; i++)
        states.push_back(state(40.0, 0.0, 274.0));
    const auto inferred = heading_offset(states, std::numeric_limits<double>::quiet_NaN());
    REQUIRE(inferred.has_value());
    CHECK(std::remainder(*inferred * 180.0 / pi - 274.0, 360.0) == Catch::Approx(0.0).margin(1e-9));
    // A hovering state drifting, nose north: north is local -y in this frame.
    const auto nose = nose_direction(state(0.1, -0.1, 4.0), inferred);
    CHECK(nose[0] == Catch::Approx(0.0).margin(1e-12));
    CHECK(nose[1] == Catch::Approx(-1.0));
    CHECK(nose[2] == 0.0);
    // Faster than velocity_nose_speed it is the velocity, unchanged.
    const State moving = state(5.0, 1.0, 4.0);
    const auto along = nose_direction(moving, inferred);
    CHECK(along[0] == moving.vx);
    CHECK(along[1] == moving.vy);
    // Between the two the nose turns continuously from the heading to the
    // velocity: no step anywhere in the band, and the endpoints are exact.
    // The heading here is west (274), the velocity north-east, so the turn
    // is the short way round.
    const auto unit_nose = [](const std::array<double, 3>& n) {
        const double horizontal = std::hypot(n[0], n[1]);
        return std::array<double, 2> {n[0] / horizontal, n[1] / horizontal};
    };
    double previous_angle = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i <= 500; i++)
    {
        const double speed = (2.0 + 0.012 * i) * 1852.0 / 3600.0; // 2 kt to 8 kt
        const State s = state(speed * std::cos(1.0), speed * std::sin(1.0), 274.0);
        const auto n = unit_nose(nose_direction(s, inferred));
        const double angle = std::atan2(n[1], n[0]);
        if (!std::isnan(previous_angle))
            CHECK(std::fabs(std::remainder(angle - previous_angle, 2.0 * pi)) < 0.02);
        previous_angle = angle;
    }
    const auto below = unit_nose(nose_direction(state(1.0, 0.0, 274.0), inferred));
    CHECK(std::fabs(std::remainder(std::atan2(below[1], below[0]) * 180.0 / pi
                  - (*inferred * 180.0 / pi - 274.0), 360.0)) == Catch::Approx(0.0).margin(1e-6));
    // A hover-only track has no fast states: the frame bearing supplies it,
    // and without one there is none.
    std::vector<State> hover(100, state(0.1, 0.0, 90.0));
    CHECK(std::remainder(*heading_offset(hover, 270.0) * 180.0 / pi - 270.0, 360.0)
            == Catch::Approx(0.0).margin(1e-9));
    CHECK_FALSE(heading_offset(hover, std::numeric_limits<double>::quiet_NaN()).has_value());
    // East on that frame is local -x.
    const auto east = nose_direction(hover.front(), heading_offset(hover, 270.0));
    CHECK(east[0] == Catch::Approx(-1.0));
    CHECK(east[1] == Catch::Approx(0.0).margin(1e-12));
    // No heading column: the velocity.
    State bare {};
    bare.vx = 0.1;
    CHECK_FALSE(heading_offset(std::vector<State>(60, bare), 270.0).has_value());
    CHECK(nose_direction(bare, 1.0)[0] == 0.1);
}

TEST_CASE("a known frame, not the track's fit, turns the hover heading", "[trajectory][nose]")
{
    // A crabbed approach: fast states track west (+x of a 270 frame) on a
    // heading of 283, 13 degrees right of the track into a wind from the right.
    // The fit takes the mean of heading against track and so absorbs the crab:
    // read through it, a hover on 283 would point the nose along 270.  With the
    // frame known the hover nose is the heading the file gives.
    State fast {};
    fast.vx = 40.0;
    fast.heading = 283.0;
    std::vector<State> states(60, fast);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto fitted = heading_offset(states, nan);
    REQUIRE(fitted.has_value());
    CHECK(std::remainder(*fitted * 180.0 / pi - 283.0, 360.0) == Catch::Approx(0.0).margin(1e-9));
    const auto framed = heading_offset(states, 270.0);
    REQUIRE(framed.has_value());
    CHECK(*framed == 270.0 * pi / 180.0);
    State hover = fast;
    hover.vx = 0.1;
    // Compass 283 is 13 degrees clockwise of the +x bearing, so local -13 degrees.
    const auto nose = nose_direction(hover, framed);
    CHECK(std::atan2(nose[1], nose[0]) * 180.0 / pi == Catch::Approx(-13.0).epsilon(1e-12));
    // The ground flag fits against the ground velocity, not the air's.
    for (auto& s : states)
    {
        s.wx = 0.0;
        s.wy = -40.0; // an air velocity of (40, 40): 45 degrees off the ground track
    }
    CHECK(std::remainder(*heading_offset(states, nan, true) * 180.0 / pi - 283.0, 360.0)
            == Catch::Approx(0.0).margin(1e-9));
    CHECK(std::remainder(*heading_offset(states, nan) * 180.0 / pi - 328.0, 360.0)
            == Catch::Approx(0.0).margin(1e-9));
}

TEST_CASE("--nose heading points the nose along the heading, crab and all", "[trajectory][nose]")
{
    REQUIRE(parse_nose_source("velocity") == NoseSource::velocity);
    REQUIRE(parse_nose_source("ground_velocity") == NoseSource::ground_velocity);
    REQUIRE(parse_nose_source("heading") == NoseSource::heading);
    REQUIRE_FALSE(parse_nose_source("Heading").has_value());
    REQUIRE_FALSE(parse_nose_source("track").has_value());
    for (const auto source : {NoseSource::velocity, NoseSource::ground_velocity, NoseSource::heading})
        REQUIRE(parse_nose_source(to_string(source)) == source);

    // A synthetic crosswind whose crab is known.  +x east (bearing 90), the
    // aircraft tracking north at 30 m/s over the ground through a 5 m/s wind
    // blowing toward the east.  Without sideslip it points along its air
    // velocity, (-5, 30): into the wind, left of the track, on a compass heading
    // of atan2(-5, 30) = -9.46 degrees.  Compass headings turn clockwise from
    // north and the local frame counterclockwise from +x, so a sign error here
    // puts the nose 19 degrees off, on the downwind side.
    State crabbed {};
    crabbed.vy = 30.0;
    crabbed.vz = 1.5;
    crabbed.wx = 5.0;
    crabbed.heading = std::atan2(-5.0, 30.0) * 180.0 / pi + 360.0;
    const double east = 90.0 * pi / 180.0;
    const auto heading = nose_direction(crabbed, east, NoseSource::heading);
    CHECK(heading[0] / heading[1] == Catch::Approx(-5.0 / 30.0).epsilon(1e-12));
    CHECK(heading[1] > 0.0);
    // The same direction the air velocity gives (--speed_reference air), climb
    // included: the climb is the air velocity's, over its horizontal speed.
    const auto air = nose_direction(crabbed, east, NoseSource::velocity);
    CHECK(air[0] == -5.0);
    CHECK(air[1] == 30.0);
    CHECK(heading[2] == Catch::Approx(1.5 / std::hypot(5.0, 30.0)).epsilon(1e-12));
    // ground_velocity ignores the wind, whatever --speed_reference set.
    const auto ground = nose_direction(crabbed, east, NoseSource::ground_velocity);
    CHECK(ground[0] == 0.0);
    CHECK(ground[1] == 30.0);
    CHECK(ground[2] == 1.5);
    // CompassToLocal, which turns the wind, agrees on the heading's direction:
    // one convention for both.
    const double h = crabbed.heading * pi / 180.0;
    const auto local = CompassToLocal(90.0)({std::sin(h), std::cos(h)});
    CHECK(std::atan2(heading[1], heading[0]) == Catch::Approx(std::atan2(local[1], local[0])).epsilon(1e-12));
    // And in the shipped tracks' frame (+x along 210), a heading of 220 is 10
    // degrees clockwise of +x: local -10 degrees, toward -y.
    State shipped {};
    shipped.vx = 40.0;
    shipped.heading = 220.0;
    const auto turned = nose_direction(shipped, 210.0 * pi / 180.0, NoseSource::heading);
    CHECK(std::atan2(turned[1], turned[0]) * 180.0 / pi == Catch::Approx(-10.0).epsilon(1e-12));
    // With no offset or no heading it is the velocity, as the others are.
    CHECK(nose_direction(shipped, std::nullopt, NoseSource::heading)[0] == 40.0);
}

TEST_CASE("--nose heading along the track gives the velocity's sphere frame", "[trajectory][nose]")
{
    // A heading that is the track's own compass bearing reproduces the default
    // nose exactly, turns, climbs and descents included, and through the sphere
    // frame's projection into a tilted tip-path plane.  This is what lets the
    // CTest nose_heading_on_track_is_velocity hold to 1e-9.
    const double bearing = 210.0;
    for (const auto& velocity : std::vector<std::array<double, 3>> {{40.0, 0.0, 0.0},
                 {30.0, -12.0, -3.0}, {-8.0, 25.0, 4.0}, {-30.0, -30.0, 0.5}, {4.2, 0.3, -1.0}})
    {
        State s {};
        s.vx = velocity[0];
        s.vy = velocity[1];
        s.vz = velocity[2];
        s.heading = bearing - std::atan2(s.vy, s.vx) * 180.0 / pi;
        const auto v = nose_direction(s, bearing * pi / 180.0, NoseSource::velocity);
        const auto n = nose_direction(s, bearing * pi / 180.0, NoseSource::heading);
        const Vec3 normal(0.08, -0.05, 1.0);
        const auto a = sphere_frame(Vec3(v[0], v[1], v[2]), normal);
        const auto b = sphere_frame(Vec3(n[0], n[1], n[2]), normal);
        CHECK((a.nose - b.nose).norm() < 1e-12);
        CHECK((a.starboard - b.starboard).norm() < 1e-12);
    }
}

TEST_CASE("--nose heading takes the climb out at hover, continuously", "[trajectory][nose]")
{
    // Below 3 kt the velocity is drift and its climb angle means nothing, so the
    // nose is the horizontal heading; by 8 kt it carries the velocity's climb in
    // full.  Between, the same smoothstep the velocity nose turns by, with no
    // step anywhere.  The horizontal direction is the heading at every speed.
    REQUIRE(nose_velocity_weight(0.0) == 0.0);
    REQUIRE(nose_velocity_weight(heading_nose_speed) == 0.0);
    REQUIRE(nose_velocity_weight(velocity_nose_speed) == 1.0);
    REQUIRE(nose_velocity_weight(100.0) == 1.0);
    const double offset = 270.0 * pi / 180.0;
    double previous = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i <= 1000; i++)
    {
        State s {};
        s.vx = 0.01 * i * 1852.0 / 3600.0 * std::cos(0.7); // 0 to 10 kt
        s.vy = 0.01 * i * 1852.0 / 3600.0 * std::sin(0.7);
        s.vz = -0.8;
        const double speed = std::hypot(s.vx, s.vy);
        s.heading = 300.0;
        const auto n = nose_direction(s, offset, NoseSource::heading);
        CHECK(std::atan2(n[1], n[0]) * 180.0 / pi == Catch::Approx(-30.0).epsilon(1e-12));
        CHECK(std::hypot(n[0], n[1]) == Catch::Approx(1.0).epsilon(1e-14));
        // The nose's elevation, which must move continuously and end at the
        // velocity's own path angle.
        const double elevation = std::atan(n[2]);
        if (speed <= heading_nose_speed)
            CHECK(n[2] == 0.0);
        if (speed >= velocity_nose_speed)
            CHECK(elevation == Catch::Approx(std::atan2(-0.8, speed)).epsilon(1e-12));
        if (!std::isnan(previous))
            CHECK(std::fabs(elevation - previous) < 0.01);
        previous = elevation;
    }
}

namespace {

// The demo AC track as text, with `change(row, fields)` applied to every data
// row first.  The shipped tracks are the only AC-format input in the repo, so
// the validation cases are built by damaging one.
template <typename Change>
std::string demo_ac_track(Change change)
{
    std::ifstream in("tracks/S-76D/197421AC.csv");
    REQUIRE(in.good());
    std::string line, text;
    bool header = true;
    std::size_t row = 0;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        if (header)
        {
            text += line + "\n";
            header = false;
            continue;
        }
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ','))
            fields.push_back(field);
        change(row++, fields);
        for (std::size_t i = 0; i < fields.size(); i++)
            text += (i ? "," : "") + fields[i];
        text += "\n";
    }
    return text;
}

// What loading `path` says is wrong with it, or empty when it loads.
std::string load_problem(const std::string& path, LengthUnit unit = feet)
{
    try
    {
        Trajectory trajectory(path, "databases/S-76D.json", demo_frame(0.02), unit);
    }
    catch (const std::exception& problem)
    {
        return problem.what();
    }
    return "";
}

} // namespace

TEST_CASE("An AC track with a value that is not a number is refused", "[trajectory][validation]")
{
    // One NaN in a column the states are built from ran through the whole
    // acoustics and came out as NaN footprint cells, exit status 0.
    const TempCsv intact(demo_ac_track([](std::size_t, std::vector<std::string>&) {}));
    REQUIRE(load_problem(intact.path).empty());

    for (std::size_t column : {1u, 5u, 6u, 7u, 11u, 13u, 15u, 21u, 22u, 23u})
    {
        for (const char* bad : {"nan", "inf", "-inf"})
        {
            const TempCsv damaged(demo_ac_track([&](std::size_t row, std::vector<std::string>& f) {
                if (row == 3000)
                    f[column] = bad;
            }));
            const auto problem = load_problem(damaged.path);
            INFO("column " << column << " set to " << bad << ": " << problem);
            REQUIRE(problem.find("not a number") != std::string::npos);
            REQUIRE(problem.find("row 3001") != std::string::npos);
        }
    }
}

TEST_CASE("An AC track whose clock does not advance is refused", "[trajectory][validation]")
{
    // The geodetic reader already refused this.  Here a constant time column
    // gave a zero integration step (an all -inf footprint) and a descending one
    // a negative step (all NaN), both exit status 0.
    const TempCsv flat(demo_ac_track([](std::size_t, std::vector<std::string>& f) { f[1] = "0.0"; }));
    REQUIRE(load_problem(flat.path).find("Time does not increase") != std::string::npos);

    const TempCsv descending(demo_ac_track([](std::size_t row, std::vector<std::string>& f) {
        f[1] = std::to_string(-0.02 * static_cast<double>(row));
    }));
    REQUIRE(load_problem(descending.path).find("Time does not increase") != std::string::npos);

    // A repeated stamp in the middle is the same fault.
    const TempCsv stalled(demo_ac_track([](std::size_t row, std::vector<std::string>& f) {
        if (row >= 1000)
            f[1] = std::to_string(0.02 * static_cast<double>(row - 1));
    }));
    REQUIRE(load_problem(stalled.path).find("Time does not increase between data rows 1000") != std::string::npos);
}

TEST_CASE("A column header that declares a conflicting unit is refused", "[trajectory][geodetic][units]")
{
    // LabeledCsv::normalize drops "(ft)" so the name matches, which left the
    // declared unit unchecked: a foot value read as meters put the aircraft 3.28
    // times too high, and "Climb Rate (fpm)" was read as feet per second.
    const auto declares = [](const std::string& problem) {
        return problem.find("declares its") != std::string::npos;
    };

    const TempCsv feet_in_feet(demo_track_as_geodetic("time,lat,lon,Alt (ft)", 50));
    REQUIRE_FALSE(declares(load_problem(feet_in_feet.path, feet)));
    REQUIRE(declares(load_problem(feet_in_feet.path, meters)));

    const TempCsv meters_header(demo_track_as_geodetic("time,lat,lon,Altitude [m MSL]", 50));
    REQUIRE(declares(load_problem(meters_header.path, feet)));
    REQUIRE_FALSE(declares(load_problem(meters_header.path, meters)));

    // An unfamiliar unit word makes no claim and passes.
    const TempCsv unknown(demo_track_as_geodetic("time,lat,lon,Alt (msl)", 50));
    REQUIRE_FALSE(declares(load_problem(unknown.path, feet)));

    const TempCsv milliseconds(demo_track_as_geodetic("Time (ms),lat,lon,alt", 50));
    REQUIRE(declares(load_problem(milliseconds.path)));

    const std::string rows = "0,39.8,-74.4,100,10,0,-5\n1,39.8001,-74.4,99,10,0,-5\n"
                             "2,39.8002,-74.4,98,10,0,-5\n3,39.8003,-74.4,97,10,0,-5\n";
    const TempCsv fpm("time,lat,lon,alt,vE,vN,Climb Rate (fpm)\n" + rows);
    REQUIRE(declares(load_problem(fpm.path)));
    const TempCsv knots("time,lat,lon,alt,vE (kt),vN (kt),vUp\n" + rows);
    REQUIRE(declares(load_problem(knots.path)));
    const TempCsv fps("time,lat,lon,alt,vE (ft/s),vN (ft/s),vUp (ft/s)\n" + rows);
    REQUIRE_FALSE(declares(load_problem(fps.path)));
    REQUIRE(declares(load_problem(fps.path, meters)));
}

TEST_CASE("A geodetic track with a value that is not a number is refused", "[trajectory][geodetic][validation]")
{
    const std::string rows = "0,39.8,-74.4,100\n1,39.8001,-74.4,99\nNAN_ROW\n3,39.8003,-74.4,97\n";
    for (const char* bad : {"2,39.8002,-74.4,nan", "2,39.8002,-74.4,inf", "nan,39.8002,-74.4,98"})
    {
        std::string text = "time,lat,lon,alt\n" + rows;
        text.replace(text.find("NAN_ROW"), 7, bad);
        const TempCsv damaged(text);
        const auto problem = load_problem(damaged.path);
        INFO(bad << ": " << problem);
        REQUIRE(problem.find("not a number") != std::string::npos);
    }
}

TEST_CASE("A geodetic track with no clock has no time of day", "[trajectory][geodetic][validation]")
{
    // --start_utc defaulted to 0, which dates a clockless flight at UTC midnight:
    // night in most timezones, so -D added the 10 dB night penalty to the whole
    // footprint with no word said.  No clock now means NaN, which -D refuses.
    const TempCsv track(demo_track_as_geodetic("time,lat,lon,alt", 50));
    const Trajectory clockless(track.path, "databases/S-76D.json", demo_frame(0.02));
    REQUIRE_FALSE(clockless.states.empty());
    REQUIRE(std::isnan(clockless.states.front().utcsec));

    GeodeticInput dated = demo_frame(0.02);
    dated.start_utc = 43200.0;
    const Trajectory noon(track.path, "databases/S-76D.json", dated);
    REQUIRE(noon.states.front().utcsec == Catch::Approx(43200.0 + noon.states.front().t));
    REQUIRE(noon.states.back().utcsec == Catch::Approx(43200.0 + noon.states.back().t));
}

TEST_CASE("watts_trim refuses a start past the end of the crop", "[trajectory][validation]")
{
    // The crop ends at the lowest or slowest point.  A start beyond it built a
    // vector from a reversed iterator range -- "Skipping ...: vector" in the log
    // and nothing else -- and a start on it left a single state.
    const auto fresh = [] { return Trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json"); };

    Trajectory fine = fresh();
    const std::size_t before = fine.states.size();
    fine.watts_trim(-1000.0);
    REQUIRE(fine.states.size() > 1);
    REQUIRE(fine.states.size() < before);

    // The demo track ends its crop at x = 100.86 m.
    for (double start : {1000.0, 101.0, 100.9})
    {
        Trajectory past = fresh();
        std::string message;
        try
        {
            past.watts_trim(start);
        }
        catch (const std::runtime_error& problem)
        {
            message = problem.what();
        }
        INFO("trim start " << start << ": " << message);
        REQUIRE(message.find("--trim") != std::string::npos);
        REQUIRE(message.find("Start the trim earlier") != std::string::npos);
    }
}

TEST_CASE("A static operation holds the aircraft at hover, nose on each hold's heading",
        "[trajectory][static]")
{
    // Two holds at one point, in feet: 30 s then 10 s, at 50 Hz.
    const TempCsv file("duration,x,y,z,heading\n30,0,0,12,210\n10,0,0,12,240\n");
    const Trajectory hold(file.path, "databases/S-76D.json");
    REQUIRE(hold.states.size() == 2000);
    for (const auto& state : hold.states)
    {
        CHECK(state.vx == 0.0);
        CHECK(state.vy == 0.0);
        CHECK(state.vz == 0.0);
        CHECK(state.LF == Catch::Approx(1.0));
        CHECK(state.mu == 0.0);
        CHECK(state.fpa == 0.0);
        CHECK(state.fpa_e == 0.0);
        CHECK(std::isfinite(state.CT));
        CHECK(state.z == Catch::Approx(12.0 * ft2m).margin(1e-9));
    }
    CHECK(hold.states.front().heading == 210.0);
    CHECK(hold.states[1499].heading == 210.0);
    CHECK(hold.states[1500].heading == 240.0);
    CHECK(hold.states[1].t - hold.states[0].t == Catch::Approx(0.02));
    CHECK(hold.states.back().t == Catch::Approx(39.98));

    // A move between holds is a trajectory, not a static operation.
    const TempCsv moving("duration,x,y,z,heading\n5,0,0,12,0\n5,10,0,12,0\n");
    std::string message;
    try
    {
        Trajectory(moving.path, "databases/S-76D.json");
    }
    catch (const std::exception& problem)
    {
        message = problem.what();
    }
    REQUIRE(message.find("share one position") != std::string::npos);
}

TEST_CASE("write the heading tracks and database for the --nose CTests", "[.heading_tracks]")
{
    // Not a test: CMakeLists.txt's nose_heading_build runs it (by its hidden tag).
    // From the demo track it writes two copies differing only in the heading
    // column (column 9, compass degrees): integrated/heading_on_track.csv, whose
    // heading is the ground track's own compass bearing in the shipped tracks'
    // frame (+x along 210), and integrated/heading_crabbed.csv, the same turned
    // 15 degrees clockwise -- a crab of +15 at every state; and two copies
    // ending in a hover (below).  And a one-sphere
    // database from the demo database whose root says azimuth_reference =
    // heading, for the database default.
    const double bearing = 210.0;
    std::ifstream in("tracks/S-76D/197421AC.csv");
    REQUIRE(in);
    ensure_parent_directory("integrated/heading_on_track.csv");
    std::ofstream on_track("integrated/heading_on_track.csv");
    std::ofstream crabbed("integrated/heading_crabbed.csv");
    // And two that end in a hover: the last 1000 states drift at 0.3 ft/s on
    // their own compass track as heading.  Before it the one flies crabbed 15
    // degrees and the other not.  A known frame turns both hovers' headings
    // alike; the track's fit, which absorbs the first one's crab, did not.
    std::ofstream hover_crabbed("integrated/heading_hover_crabbed.csv");
    std::ofstream hover_straight("integrated/heading_hover_straight.csv");
    std::vector<std::string> lines;
    std::string line;
    REQUIRE(std::getline(in, line));
    for (auto* out : {&on_track, &crabbed, &hover_crabbed, &hover_straight})
        *out << line << '\n';
    while (std::getline(in, line))
        lines.push_back(line);
    std::size_t rows = 0;
    for (std::size_t row = 0; row < lines.size(); row++)
    {
        std::vector<std::string> fields;
        std::stringstream split(lines[row]);
        for (std::string field; std::getline(split, field, ',');)
            fields.push_back(field);
        if (fields.size() != 27)
            continue;
        const double vx = std::stod(fields[11]), vy = std::stod(fields[13]);
        const double track = bearing - std::atan2(vy, vx) * 180.0 / pi;
        const auto number = [](double value) {
            std::ostringstream text;
            text << std::fixed << std::setprecision(12) << std::setw(16) << value;
            return text.str();
        };
        const auto write = [&](std::ofstream& out, double heading) {
            fields[8] = number(heading);
            for (std::size_t k = 0; k < fields.size(); k++)
                out << (k ? "," : "") << fields[k];
            out << '\n';
        };
        write(on_track, track);
        write(crabbed, track + 15.0);
        if (row + 1000 >= lines.size())
        {
            const double scale = 0.3 / std::hypot(vx, vy);
            fields[11] = number(vx * scale);
            fields[13] = number(vy * scale);
            write(hover_crabbed, track);
            write(hover_straight, track);
        }
        else
        {
            write(hover_crabbed, track + 15.0);
            write(hover_straight, track);
        }
        rows++;
    }
    REQUIRE(rows > 1000);
    Database source("databases/S-76D_M3.nod");
    auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    sphere.azimuth_reference = AzimuthReference::heading;
    write_sphere(sphere, "integrated/heading_sphere.nod");
}

TEST_CASE("A WGS84 track's heading column is the nose's compass heading", "[trajectory][geodetic]")
{
    // The AC file's column 9 is the heading; written beside its lat/lon/alt, the
    // WGS84 reader has to hand back the same heading at every state (0.02 s
    // samples, already uniform, so nothing is resampled).
    const auto records = read_numeric_csv("tracks/S-76D/197421AC.csv", /*skip_first_record=*/true);
    std::ostringstream out;
    out << "time,lat,lon,alt,heading (deg)\n" << std::setprecision(17);
    for (const auto& r : records)
        out << r[1] << "," << r[2] << "," << r[3] << "," << r[4] << "," << r[8] << "\n";
    const TempCsv track(out.str());
    const Trajectory geodetic(track.path, "databases/S-76D.json", demo_frame(0.02));
    const Trajectory reference("tracks/S-76D/197421AC.csv", "databases/S-76D.json");
    REQUIRE(geodetic.states.size() == reference.states.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < geodetic.states.size(); i++)
        worst = std::max(worst, std::fabs(std::remainder(
                                        geodetic.states[i].heading - reference.states[i].heading, 360.0)));
    REQUIRE(worst < 1e-9);
}

TEST_CASE("A WGS84 heading is interpolated through north, not around through south",
        "[trajectory][geodetic]")
{
    // 350 then 10 degrees a second apart: halfway is north.  Interpolated as
    // plain numbers it would be 180, the aircraft facing backward.
    const TempCsv track("time,lat,lon,alt,heading\n0,39.8,-74.4,1000,350\n1,39.8001,-74.4,1000,10\n"
                        "2,39.8002,-74.4,1000,30\n");
    GeodeticInput input = demo_frame();
    input.resample_dt = 0.5;
    const Trajectory trajectory(track.path, "databases/S-76D.json", input);
    REQUIRE(trajectory.states.size() == 5);
    REQUIRE(std::fabs(std::remainder(trajectory.states[1].heading, 360.0)) < 1e-9);
    REQUIRE(trajectory.states[1].heading >= 0.0);
    REQUIRE(trajectory.states[1].heading < 360.0);
    REQUIRE(std::fabs(trajectory.states[3].heading - 20.0) < 1e-9);
}

TEST_CASE("A WGS84 heading in radians is refused", "[trajectory][geodetic]")
{
    const TempCsv radians("time,lat,lon,alt,heading (rad)\n0,39.8,-74.4,1000,3.1\n1,39.8001,-74.4,1000,3.1\n"
                          "2,39.8002,-74.4,1000,3.1\n");
    REQUIRE_THROWS_WITH(Trajectory(radians.path, "databases/S-76D.json", demo_frame()),
            Catch::Matchers::ContainsSubstring("declares its heading in rad"));
}

TEST_CASE("A WGS84 heading must be true degrees, and only one", "[trajectory][geodetic]")
{
    // Read as true, a magnetic heading turns the nose by the declination.
    const TempCsv magnetic("time,lat,lon,alt,heading (mag)\n0,39.8,-74.4,1000,30\n1,39.8001,-74.4,1000,30\n"
                           "2,39.8002,-74.4,1000,30\n");
    REQUIRE_THROWS_WITH(Trajectory(magnetic.path, "databases/S-76D.json", demo_frame()),
            Catch::Matchers::ContainsSubstring("is a magnetic heading"));
    const TempCsv turns("time,lat,lon,alt,hdg (turns)\n0,39.8,-74.4,1000,0.1\n1,39.8001,-74.4,1000,0.1\n"
                        "2,39.8002,-74.4,1000,0.1\n");
    REQUIRE_THROWS_WITH(Trajectory(turns.path, "databases/S-76D.json", demo_frame()),
            Catch::Matchers::ContainsSubstring("declares its heading in turns"));
    // Two candidates: the first match is not necessarily the true one.
    const TempCsv two("time,lat,lon,alt,hdg,true heading\n0,39.8,-74.4,1000,18,30\n1,39.8001,-74.4,1000,18,30\n"
                      "2,39.8002,-74.4,1000,18,30\n");
    REQUIRE_THROWS_WITH(Trajectory(two.path, "databases/S-76D.json", demo_frame()),
            Catch::Matchers::ContainsSubstring("has 2 heading columns"));
    // Degrees, declared or not, true or not said, are read.
    const TempCsv degrees("time,lat,lon,alt,heading (deg true)\n0,39.8,-74.4,1000,30\n1,39.8001,-74.4,1000,30\n"
                          "2,39.8002,-74.4,1000,30\n");
    REQUIRE_NOTHROW(Trajectory(degrees.path, "databases/S-76D.json", demo_frame()));
}
