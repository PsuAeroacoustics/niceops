#include "trajectory.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

#include <nlohmann/json.hpp>

#include "splines.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

constexpr double PI = 3.14159265358979323846;

void diff(double dx, const std::vector<double>& f, std::vector<double>& df, std::vector<double>& ddf);

void smoother(const std::vector<double>& x, const std::vector<double>& f,
        const std::vector<double>& xi, std::vector<double>& fi, std::vector<double>& dfi,
        std::vector<double>& ddfi, int degree, double alpha)
{
    Spline spline(x, f, degree, alpha);
    // Use spline for interpolation
    fi = spline.eval(xi);
    // Use finite differences to get velocity and acceleration
    const double dx = xi[1] - xi[0];
    diff(dx, fi, dfi, ddfi);
}

void diff(double dx, const std::vector<double>& f, std::vector<double>& df, std::vector<double>& ddf)
{
    // Create storage for first derivative
    df = f;
    // Forward diff first point
    df[0] = (f[1] - f[0]) / dx;
    // Central diff interior points
    for (std::size_t i = 1; i + 1 < f.size(); i++)
        df[i] = 0.5 * (f[i + 1] - f[i - 1]) / dx;
    // Backward diff last point
    df[f.size() - 1] = (f[f.size() - 1] - f[f.size() - 2]) / dx;
    // Create storage for second derivative
    ddf = df;
    // Forward diff first point
    ddf[0] = (f[2] - 2 * f[1] + f[0]) / (dx * dx);
    // Central diff interior points
    for (std::size_t i = 1; i + 1 < f.size(); i++)
        ddf[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (dx * dx);
    // Backward diff last point
    ddf[f.size() - 1] = (f[f.size() - 1] - 2 * f[f.size() - 2] + f[f.size() - 3]) / (dx * dx);
}

// The vehicle JSONs store every value as a quoted string, and the D original
// read them with `.str` followed by `to!double`.  nlohmann's get<double>()
// throws on a string, so convert explicitly (tolerating a real JSON number).
double json_to_double(const nlohmann::json& value)
{
    if (value.is_string())
        return std::stod(value.get<std::string>());
    return value.get<double>();
}

dimensional_data configure_vehicle(const std::string& filename)
{
    info("Reading vehicle configuration: ", filename);
    std::ifstream in(filename);
    if (!in)
        fatal("Unable to open vehicle configuration file: ", filename);
    auto config = nlohmann::json::parse(in);

    dimensional_data vehicle;
    vehicle.radius = json_to_double(config["Main Rotor"]["radius"]);
    vehicle.tip_speed = json_to_double(config["Main Rotor"]["tip speed"]);
    vehicle.weight = json_to_double(config["Vehicle"]["weight"]) * g; // Actually a mass (kg), so convert to Newtons
    vehicle.rho = json_to_double(config["Atmosphere"]["density"]);
    trace("Vehicle data structure: ", vehicle);
    return vehicle;
}

} // namespace

Trajectory::Trajectory(const std::string& filename, const std::string& vehicle_filename)
{
    // Get vehicle data
    vehicle = configure_vehicle(vehicle_filename);

    // D: filename.splitter('.').retro.take(1).array()[0] -- the last '.' segment
    auto dotpos = filename.rfind('.');
    std::string extension = (dotpos == std::string::npos) ? filename : filename.substr(dotpos + 1);
    std::for_each(extension.begin(), extension.end(), [](char& c) { c = std::tolower(c); });
    info("Opening file ", filename, " with extension ", extension);
    if (extension == "csv")
    {
        trace("Loading ", filename, " as CSV.");
        // D reads this with csvReader!(double)(contents, null), which discards
        // the first record as a header regardless of content.  The demo tracks
        // all carry a text header so both paths agreed by luck; a headerless
        // trajectory would have gained a spurious first state here.
        auto records = read_numeric_csv(filename, /*skip_first_record=*/true);

        for (const auto& data : records)
        {
            if (data.size() == 27)
            {
                // Assume NASA/FAA AC state data format
                State now {};
                if (std::isnan(utcstart))
                {
                    utcstart = data[0];
                }
                // Convert to metric
                now.t = data[1];
                now.x = data[5] * ft2m;
                now.y = data[6] * ft2m;
                now.z = data[7] * ft2m;
                now.vx = data[11] * ft2m;
                now.vy = data[13] * ft2m;
                now.vz = data[15] * ft2m;
                // Convert to "g"
                now.nx = data[21] * fps2g;
                now.ny = data[22] * fps2g;
                now.nz = (data[23] * fps2g + 1.0);
                // Compute load factor
                now.LF = std::sqrt(now.nx * now.nx + now.ny * now.ny + now.nz * now.nz);
                now.nx /= now.LF;
                now.ny /= now.LF;
                now.nz /= now.LF;
                // Non-dimensionalize
                now.CT = now.LF * vehicle.weight / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed
                        * PI * vehicle.radius * vehicle.radius);
                auto v = std::sqrt(now.vx * now.vx + now.vy * now.vy + now.vz * now.vz);
                now.mu = v / vehicle.tip_speed;
                auto vn_hat = (now.vx * now.nx + now.vy * now.ny + now.vz * now.nz) / v;
                now.fpa_e = 180.0 / PI * std::asin(vn_hat);
                now.fpa = 180.0 / PI * std::asin(now.vz / v);
                states.push_back(now);
                trace("Parsed state: ", states.back());
            }
            else
            {
                fatal("Wrong number of columns.");
            }
        }
    }
    else if (extension == "nc")
    {
        fatal("netCDF files not yet supported.");
    }
    else
    {
        fatal("File extension ", extension, " not recognized.");
    }
}

void Trajectory::write_state(const std::string& state_file) const
{
    info("Writing state data to ", state_file);
    std::ofstream file(state_file);
    if (!file)
        fatal("Unable to open state file for writing: ", state_file);
    // Write header from field names (mirrors __traits(allMembers, State))
    file << "t,x,y,z,vx,vy,vz,nx,ny,nz,LF,CT,mu,fpa_e,fpa\n";
    file << std::setprecision(17);
    for (const auto& state : states)
    {
        file << state.t << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "," << state.nx << "," << state.ny << "," << state.nz << ","
            << state.LF << "," << state.CT << "," << state.mu << "," << state.fpa_e << "," << state.fpa
            << "\n";
    }
}

void Trajectory::adapt(double dt, int degree, double alpha)
{
    info("Adapting spline to dt ", dt, " with spline of degree ", degree, " and smoothing ", alpha);
    assert(dt > 0.0);

    // Unpack states
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    for (const auto& now : states)
    {
        t.push_back(now.t);
        x.push_back(now.x);
        y.push_back(now.y);
        z.push_back(now.z);
    }
    // Build interpolants
    std::vector<double> ti;
    ti.push_back(t[0]);
    double tnext = ti[0] + dt;
    while (tnext <= t[t.size() - 1])
    {
        ti.push_back(tnext);
        tnext += dt;
    }
    trace("Smoothing x.");
    std::vector<double> xi;
    std::vector<double> vxi;
    std::vector<double> axi;
    smoother(t, x, ti, xi, vxi, axi, degree, alpha);

    trace("Smoothing y.");
    std::vector<double> yi;
    std::vector<double> vyi;
    std::vector<double> ayi;
    smoother(t, y, ti, yi, vyi, ayi, degree, alpha);

    trace("Smoothing z.");
    std::vector<double> zi;
    std::vector<double> vzi;
    std::vector<double> azi;
    smoother(t, z, ti, zi, vzi, azi, degree, alpha);

    // Repack data into states
    states.clear();
    trace("States reset to ", states);
    for (std::size_t i = 0; i < ti.size(); i++)
    {
        State now {};
        now.t = ti[i];
        now.x = xi[i];
        now.y = yi[i];
        now.z = zi[i];
        now.vx = vxi[i];
        now.vy = vyi[i];
        now.vz = vzi[i];
        // Convert to "g"
        now.nx = axi[i] * mps2g;
        now.ny = ayi[i] * mps2g;
        now.nz = azi[i] * mps2g + 1.0;
        // Compute load factor
        now.LF = std::sqrt(now.nx * now.nx + now.ny * now.ny + now.nz * now.nz);
        // Normalize
        now.nx /= now.LF;
        now.ny /= now.LF;
        now.nz /= now.LF;
        // Non-dimensionalize
        now.CT = now.LF * vehicle.weight / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed * PI
                * vehicle.radius * vehicle.radius);
        auto v = std::sqrt(now.vx * now.vx + now.vy * now.vy + now.vz * now.vz);
        now.mu = v / vehicle.tip_speed;
        auto vn_hat = (now.vx * now.nx + now.vy * now.ny + now.vz * now.nz) / v;
        now.fpa_e = 180.0 / PI * std::asin(vn_hat);
        now.fpa = 180.0 / PI * std::asin(now.vz / v);
        states.push_back(now);
        trace("Adapted state: ", states.back());
    }
}

std::vector<double> Trajectory::column(const std::string& varname) const
{
    static const std::map<std::string, double State::*> fields = {{"t", &State::t}, {"x", &State::x},
        {"y", &State::y}, {"z", &State::z}, {"vx", &State::vx}, {"vy", &State::vy}, {"vz", &State::vz},
        {"nx", &State::nx}, {"ny", &State::ny}, {"nz", &State::nz}, {"LF", &State::LF}, {"CT", &State::CT},
        {"mu", &State::mu}, {"fpa_e", &State::fpa_e}, {"fpa", &State::fpa}};
    auto it = fields.find(varname);
    if (it == fields.end())
        fatal("Unknown state variable: ", varname);
    std::vector<double> x;
    x.reserve(states.size());
    for (const auto& state : states)
    {
        x.push_back(state.*(it->second));
    }
    return x;
}

void Trajectory::watts_trim(double xstart_ft)
{
    auto mu = column("mu");
    auto z = column("z");
    auto imin_z = std::min_element(z.begin(), z.end()) - z.begin();
    auto imin_v = std::min_element(mu.begin(), mu.end()) - mu.begin();
    auto ixmax = imin_z < imin_v ? imin_z : imin_v;
    auto x = column("x");
    std::size_t ixmin = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        if (x[i] > xstart_ft * ft2m)
            break;
        ixmin = i;
    }
    info("Trimming trajectory from x = ", x[ixmin], " to ", x[ixmax]);
    std::vector<State> sub(states.begin() + ixmin, states.begin() + ixmax + 1);
    states.swap(sub);
}
