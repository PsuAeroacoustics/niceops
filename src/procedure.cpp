#include "procedure.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>

#include <nlohmann/json.hpp>

#include "util/logger.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

constexpr double gravity = 9.80665; // m/s^2
constexpr double knot = 1852.0 / 3600.0; // m/s

// A key with its unit: `name_ft` or `name_m` for a length, `name_kt`, `name_mps`
// or `name_fps` for a speed, `name_deg` for an angle, `name_s` for a time.
// Exactly one spelling may be given; none returns NaN when optional.
double read_quantity(const nlohmann::json& object, const std::string& where, const std::string& name,
        const std::map<std::string, double>& units, bool required)
{
    double value = std::numeric_limits<double>::quiet_NaN();
    std::string found;
    for (const auto& [suffix, scale] : units)
    {
        const std::string key = name + "_" + suffix;
        if (!object.contains(key))
            continue;
        if (!found.empty())
            input_error(where, " gives both ", found, " and ", key, "; give one.");
        if (!object[key].is_number())
            input_error(where, ": ", key, " must be a number.");
        value = object[key].get<double>() * scale;
        found = key;
    }
    if (object.contains(name))
        input_error(where, ": ", name, " needs its unit in its name (", name, "_", units.begin()->first,
                ", ...); a bare number is refused.");
    if (required && found.empty())
        input_error(where, " needs ", name, " (", name, "_", units.begin()->first, ", ...).");
    if (!found.empty() && !std::isfinite(value))
        input_error(where, ": ", found, " is not finite.");
    return value;
}

const std::map<std::string, double> lengths = {{"ft", ft2m}, {"m", 1.0}};
const std::map<std::string, double> speeds = {{"kt", knot}, {"mps", 1.0}, {"fps", ft2m}};
const std::map<std::string, double> accelerations = {{"kt_per_s", knot}, {"mps2", 1.0}};
const std::map<std::string, double> angles = {{"deg", deg_to_rad}};
const std::map<std::string, double> times = {{"s", 1.0}};

} // namespace

Procedure niceops::read_procedure(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        input_error("Cannot open procedure ", path, ".");
    nlohmann::json doc;
    try
    {
        in >> doc;
    }
    catch (const std::exception& e)
    {
        input_error("Procedure ", path, " is not JSON: ", e.what());
    }
    Procedure p;
    if (!doc.contains("start") || !doc["start"].is_object())
        input_error("Procedure ", path, " needs a \"start\" object.");
    const auto& start = doc["start"];
    const std::string where = "Procedure " + path + ", start";
    for (const char* key : {"lat", "lon"})
        if (!start.contains(key) || !start[key].is_number())
            input_error(where, " needs ", key, " (degrees).");
    p.start = Geodetic {start["lat"].get<double>(), start["lon"].get<double>(), 0.0};
    if (!(std::fabs(p.start.latitude) <= 90.0) || !(std::fabs(p.start.longitude) <= 180.0))
        input_error(where, " is not a WGS84 position.");
    p.altitude = read_quantity(start, where, "alt", lengths, true);
    p.speed = read_quantity(start, where, "speed", speeds, true);
    p.course = read_quantity(start, where, "course", angles, true);
    const double ground = read_quantity(start, where, "ground_alt", lengths, false);
    p.ground_altitude = std::isnan(ground) ? 0.0 : ground;
    if (p.altitude < p.ground_altitude)
        input_error(where, ": alt is below ground_alt.");
    if (!(p.speed >= 0.0))
        input_error(where, ": the speed must not be negative.");
    if (doc.contains("step_s"))
        p.step = doc["step_s"].get<double>();
    if (!(p.step > 0.0))
        input_error("Procedure ", path, ": step_s must be positive.");
    if (doc.contains("limits"))
    {
        const auto& limits = doc["limits"];
        const std::string lw = "Procedure " + path + ", limits";
        if (limits.contains("pull_g"))
            p.pull = limits["pull_g"].get<double>();
        const double roll = read_quantity(limits, lw, "roll_rate", {{"deg_per_s", deg_to_rad}}, false);
        if (!std::isnan(roll))
            p.roll_rate = roll;
        const double accel = read_quantity(limits, lw, "acceleration", accelerations, false);
        if (!std::isnan(accel))
            p.acceleration = accel;
        if (!(p.pull > 0.0) || !(p.roll_rate > 0.0) || !(p.acceleration > 0.0))
            input_error(lw, ": pull_g, roll_rate and acceleration must be positive.");
    }
    if (!doc.contains("segments") || !doc["segments"].is_array() || doc["segments"].empty())
        input_error("Procedure ", path, " needs a non-empty \"segments\" list.");
    std::size_t index = 0;
    for (const auto& item : doc["segments"])
    {
        index++;
        const std::string w = "Procedure " + path + ", segment " + std::to_string(index);
        if (!item.contains("type") || !item["type"].is_string())
            input_error(w, " needs a type: level, climb, speed, turn or hover.");
        const std::string type = item["type"].get<std::string>();
        ProcedureSegment s;
        s.duration = read_quantity(item, w, "duration", times, false);
        if (type == "level")
        {
            s.kind = ProcedureSegment::Kind::level;
            s.distance = read_quantity(item, w, "distance", lengths, false);
            if (std::isnan(s.duration) == std::isnan(s.distance))
                input_error(w, " (level) needs duration_s or distance, one of them.");
        }
        else if (type == "climb" || type == "descent")
        {
            s.kind = ProcedureSegment::Kind::climb;
            s.path_angle = read_quantity(item, w, "angle", angles, true);
            if (type == "descent" && s.path_angle > 0.0)
                s.path_angle = -s.path_angle;
            s.until_altitude = read_quantity(item, w, "until_alt", lengths, false);
            if (std::isnan(s.duration) == std::isnan(s.until_altitude))
                input_error(w, " (", type, ") needs duration_s or until_alt, one of them.");
            if (s.path_angle == 0.0)
                input_error(w, " (", type, ") needs a path angle other than zero; use level.");
        }
        else if (type == "speed")
        {
            s.kind = ProcedureSegment::Kind::speed;
            s.target_speed = read_quantity(item, w, "to_speed", speeds, true);
            s.acceleration = read_quantity(item, w, "rate", accelerations, false);
            const double angle = read_quantity(item, w, "angle", angles, false);
            s.path_angle = std::isnan(angle) ? 0.0 : angle;
            if (!(s.target_speed >= 0.0))
                input_error(w, " (speed): to_speed must not be negative.");
        }
        else if (type == "turn")
        {
            s.kind = ProcedureSegment::Kind::turn;
            s.bank = read_quantity(item, w, "bank", angles, true);
            s.turn = read_quantity(item, w, "turn", angles, true);
            if (!(std::fabs(s.bank) > 0.0 && std::fabs(s.bank) < 60.0 * deg_to_rad))
                input_error(w, " (turn): bank_deg must be between 0 and 60 degrees, signed (+ right).");
            if (!(s.turn > 0.0))
                input_error(w, " (turn): turn_deg must be positive; the sign of bank_deg says which way.");
        }
        else if (type == "hover")
        {
            s.kind = ProcedureSegment::Kind::hover;
            if (std::isnan(s.duration))
                input_error(w, " (hover) needs duration_s.");
        }
        else
        {
            input_error(w, ": unknown type \"", type, "\"; level, climb, descent, speed, turn or hover.");
        }
        if (!std::isnan(s.duration) && !(s.duration > 0.0))
            input_error(w, ": duration_s must be positive.");
        p.segments.push_back(s);
    }
    return p;
}

std::vector<ProcedureState> niceops::fly(const Procedure& p)
{
    // Integrated at a fine internal step, east/north/up about the start; the
    // step divides the written one exactly.
    const double h = p.step / std::ceil(p.step / 0.01);
    const LocalFrame frame(p.start);
    // Time counts internal steps, so a long procedure's written states stay
    // evenly spaced, which spares niceops resampling them.
    std::size_t steps = 0;
    const auto every = static_cast<std::size_t>(std::llround(p.step / h));
    double t = 0.0, east = 0.0, north = 0.0, z = p.altitude;
    double speed = p.speed, course = p.course, gamma = 0.0, bank = 0.0;
    std::vector<ProcedureState> states;
    const auto write = [&]() {
        if (steps % every != 0)
            return;
        ProcedureState s;
        s.t = t;
        s.position = frame.to_geodetic(Vec3(east, north, 0.0));
        s.altitude = z;
        // + 0.0: a stationary state writes 0, not -0.
        s.east = speed * std::cos(gamma) * std::sin(course) + 0.0;
        s.north = speed * std::cos(gamma) * std::cos(course) + 0.0;
        s.up = speed * std::sin(gamma) + 0.0;
        double heading = std::fmod(course * rad_to_deg, 360.0);
        s.heading = heading < 0.0 ? heading + 360.0 : heading;
        states.push_back(s);
    };
    // Move a toward b by at most rate * h.
    const auto approach = [h](double a, double b, double rate) {
        const double change = rate * h;
        return b > a ? std::min(b, a + change) : std::max(b, a - change);
    };
    // One internal step toward the targets.
    const auto advance = [&](double target_gamma, double target_bank, double target_speed, double acceleration) {
        // The path angle turns at the limited pull: d gamma / dt = g n / V.
        const double gamma_rate = speed > 1.0 ? gravity * p.pull / speed : 2.0 * deg_to_rad;
        gamma = approach(gamma, target_gamma, gamma_rate);
        bank = approach(bank, target_bank, p.roll_rate);
        speed = approach(speed, target_speed, acceleration);
        if (speed > 1.0)
            course += gravity * std::tan(bank) / speed * h;
        east += speed * std::cos(gamma) * std::sin(course) * h;
        north += speed * std::cos(gamma) * std::cos(course) * h;
        z += speed * std::sin(gamma) * h;
        steps++;
        t = static_cast<double>(steps) * h;
        write();
    };
    write();
    std::size_t index = 0;
    for (const auto& s : p.segments)
    {
        index++;
        const std::string w = "Procedure segment " + std::to_string(index);
        const std::size_t before = states.size();
        const double start_t = t;
        // A safety net: no segment flies longer than this.
        const double limit = 4.0 * 3600.0;
        switch (s.kind)
        {
        case ProcedureSegment::Kind::level:
        {
            double travelled = 0.0;
            while (std::isnan(s.duration) ? travelled < s.distance : t < start_t + s.duration - 0.5 * h)
            {
                const double before_e = east, before_n = north;
                advance(0.0, 0.0, speed, p.acceleration);
                travelled += std::hypot(east - before_e, north - before_n);
                if (t - start_t > limit || (!std::isnan(s.distance) && speed < 0.1))
                    input_error(w, " (level) at ", speed / knot, " kt never covers its distance.");
            }
            break;
        }
        case ProcedureSegment::Kind::climb:
        {
            if (speed < 0.1)
                input_error(w, " climbs or descends from rest: give it a speed first.");
            const bool up = s.path_angle > 0.0;
            while (std::isnan(s.duration) ? (up ? z < s.until_altitude : z > s.until_altitude)
                                          : t < start_t + s.duration - 0.5 * h)
            {
                advance(s.path_angle, 0.0, speed, p.acceleration);
                if (t - start_t > limit)
                    input_error(w, " never reaches its altitude.");
            }
            break;
        }
        case ProcedureSegment::Kind::speed:
        {
            const double rate = std::isnan(s.acceleration) ? p.acceleration : s.acceleration;
            if (!(rate > 0.0))
                input_error(w, " (speed): rate must be positive.");
            while (speed != s.target_speed)
            {
                advance(s.path_angle, 0.0, s.target_speed, rate);
                if (t - start_t > limit)
                    input_error(w, " never reaches its speed.");
            }
            break;
        }
        case ProcedureSegment::Kind::turn:
        {
            if (speed < 5.0)
                input_error(w, " turns at ", speed / knot, " kt; a coordinated turn needs way on.");
            // Roll out early by the course the roll-out itself still turns.
            const double roll_time = std::fabs(s.bank) / p.roll_rate;
            const double lead = 0.5 * gravity * std::tan(std::fabs(s.bank)) / speed * roll_time;
            const double start_course = course;
            while (std::fabs(course - start_course) < std::max(s.turn - lead, 0.0))
            {
                advance(0.0, s.bank, speed, p.acceleration);
                if (t - start_t > limit)
                    input_error(w, " never completes its turn.");
            }
            break;
        }
        case ProcedureSegment::Kind::hover:
        {
            if (speed > 0.0 || std::fabs(gamma) > 0.0 || std::fabs(bank) > 0.0)
                input_error(w, " (hover) starts at ", speed / knot, " kt: slow to 0 first (a speed segment).");
            while (t < start_t + s.duration - 0.5 * h)
                advance(0.0, 0.0, 0.0, p.acceleration);
            break;
        }
        }
        for (std::size_t k = before; k < states.size(); k++)
            if (states[k].altitude < p.ground_altitude - 1e-9)
                input_error(w, " goes below the ground (ground_alt) at t = ", states[k].t, " s.");
    }
    return states;
}

void niceops::write_procedure_track(const std::string& path, const std::vector<ProcedureState>& states,
        LengthUnit unit)
{
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out)
        input_error("Cannot write ", path, ".");
    const bool feet_unit = unit.to_meters != 1.0;
    out << "time,lat,lon," << (feet_unit ? "alt_ft" : "alt_m") << ",vE,vN,vUp,heading\n";
    out << std::setprecision(12);
    double last_t = -1.0;
    for (const auto& s : states)
    {
        if (s.t <= last_t)
            continue; // the closing write can repeat the last sample
        last_t = s.t;
        out << s.t << "," << s.position.latitude << "," << s.position.longitude << ","
            << s.altitude / unit.to_meters << "," << s.east / unit.to_meters << "," << s.north / unit.to_meters
            << "," << s.up / unit.to_meters << "," << s.heading << "\n";
    }
    if (!out.flush())
        input_error("Could not finish writing ", path, ".");
}
