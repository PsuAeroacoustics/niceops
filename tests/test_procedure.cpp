#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "procedure.hpp"
#include "trajectory.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

constexpr double knot = 1852.0 / 3600.0;

std::string write_text(const std::string& name, const std::string& text)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream(path) << text;
    return path;
}

Procedure start_at(double speed_kt, double course_deg, double alt_ft)
{
    Procedure p;
    p.start = Geodetic {39.8, -74.4, 0.0};
    p.altitude = alt_ft * ft2m;
    p.speed = speed_kt * knot;
    p.course = course_deg * deg_to_rad;
    return p;
}

} // namespace

TEST_CASE("level flight goes straight along its course at its speed", "[procedure]")
{
    auto p = start_at(100.0, 90.0, 1000.0);
    ProcedureSegment level;
    level.duration = 10.0;
    p.segments.push_back(level);
    const auto states = fly(p);
    REQUIRE(states.size() == 51); // 0 to 10 s every 0.2 s
    const LocalFrame frame(p.start);
    const Vec3 end = frame.to_local(states.back().position);
    CHECK(std::fabs(end(0) - 100.0 * knot * 10.0) < 1e-6);
    CHECK(std::fabs(end(1)) < 1e-6);
    CHECK(states.back().altitude == p.altitude);
    CHECK(std::fabs(states.back().heading - 90.0) < 1e-9);
}

TEST_CASE("a descent holds its path angle and stops at its altitude", "[procedure]")
{
    auto p = start_at(80.0, 0.0, 1500.0);
    ProcedureSegment descent;
    descent.kind = ProcedureSegment::Kind::climb;
    descent.path_angle = -6.0 * deg_to_rad;
    descent.until_altitude = 500.0 * ft2m;
    p.segments.push_back(descent);
    const auto states = fly(p);
    // In the steady part the climb rate over the speed is sin(-6 deg).
    const auto& mid = states[states.size() / 2];
    const double speed = std::sqrt(mid.east * mid.east + mid.north * mid.north + mid.up * mid.up);
    CHECK(std::fabs(mid.up / speed - std::sin(-6.0 * deg_to_rad)) < 1e-9);
    CHECK(states.back().altitude <= 500.0 * ft2m + 1.0);
    CHECK(states.back().altitude > 500.0 * ft2m - 2.0);
}

TEST_CASE("a speed segment changes speed at its rate", "[procedure]")
{
    auto p = start_at(80.0, 0.0, 500.0);
    ProcedureSegment slow;
    slow.kind = ProcedureSegment::Kind::speed;
    slow.target_speed = 40.0 * knot;
    slow.acceleration = 1.0 * knot; // 1 kt/s
    p.segments.push_back(slow);
    const auto states = fly(p);
    CHECK(std::fabs(states.back().t - 40.0) < 0.21);
    CHECK(std::fabs(states.back().north - 40.0 * knot) < 0.6);
}

TEST_CASE("niceops reads a flown turn's load factor as 1 / cos(bank)", "[procedure]")
{
    // A 30 degree coordinated turn through 180 degrees at 80 kt, written as a
    // WGS84 track and read back as niceops reads it (velocity columns,
    // differentiated for the load factor).
    auto p = start_at(80.0, 210.0, 1000.0);
    ProcedureSegment lead_in;
    lead_in.duration = 5.0;
    ProcedureSegment turn;
    turn.kind = ProcedureSegment::Kind::turn;
    turn.bank = 30.0 * deg_to_rad;
    turn.turn = 180.0 * deg_to_rad;
    ProcedureSegment lead_out;
    lead_out.duration = 5.0;
    p.segments = {lead_in, turn, lead_out};
    const auto states = fly(p);
    // Through 180 degrees, to within the roll-out's own turning.
    CHECK(std::fabs(std::remainder(states.back().heading - 30.0, 360.0)) < 3.0);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_procedure_turn.csv").string();
    write_procedure_track(path, states, feet);
    GeodeticInput input;
    input.frame = LocalFrame(p.start);
    const Trajectory trajectory(path, "databases/S-76D.json", input);
    std::filesystem::remove(path);
    // The middle of the turn: steady bank.
    const auto& s = trajectory.states[trajectory.states.size() / 2];
    CHECK(std::fabs(s.LF - 1.0 / std::cos(30.0 * deg_to_rad)) < 1e-6);
    // And the heading column came through.
    CHECK(std::isfinite(s.heading));
}

TEST_CASE("a procedure names its units, and refuses what cannot be flown", "[procedure]")
{
    const auto bare = write_text("niceops_procedure_bare.json",
            R"({"start": {"lat": 39.8, "lon": -74.4, "alt_ft": 1000, "speed": 80, "course_deg": 0},
                "segments": [{"type": "level", "duration_s": 10}]})");
    REQUIRE_THROWS_WITH(read_procedure(bare), Catch::Matchers::ContainsSubstring("a bare number is refused"));
    const auto hover = write_text("niceops_procedure_hover.json",
            R"({"start": {"lat": 39.8, "lon": -74.4, "alt_ft": 100, "speed_kt": 40, "course_deg": 0},
                "segments": [{"type": "hover", "duration_s": 10}]})");
    REQUIRE_THROWS_WITH(fly(read_procedure(hover)), Catch::Matchers::ContainsSubstring("slow to 0 first"));
    const auto both = write_text("niceops_procedure_both.json",
            R"({"start": {"lat": 39.8, "lon": -74.4, "alt_ft": 100, "alt_m": 30, "speed_kt": 40, "course_deg": 0},
                "segments": [{"type": "level", "duration_s": 10}]})");
    REQUIRE_THROWS_WITH(read_procedure(both), Catch::Matchers::ContainsSubstring("gives both alt_ft and alt_m"));
    std::filesystem::remove(bare);
    std::filesystem::remove(hover);
    std::filesystem::remove(both);
}
