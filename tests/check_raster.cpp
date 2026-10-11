// Checks what --raster_file wrote, for the CTest cases that run it
// (CMakeLists.txt, raster_*):
//
//   check_raster RASTER.nc UNIT_TO_METERS [--latlon LATLON.csv | --no-latlon] VAR=GRID.csv ...
//
//   * the file is CF-1.8, with x and y in meters that are the grid file's axes
//     (written in the run's units, UNIT_TO_METERS converting);
//   * each VAR is the grid file's levels to 1e-12 relative, a non-finite level
//     (an observer that heard nothing) as the fill value, NaN;
//   * with --latlon, the latitude and longitude are the --latlon companion's,
//     observer by observer, and every VAR names them as its coordinates;
//     with --no-latlon (a run without a frame), there are none.
//
// Exits nonzero with the first failure on stderr.

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& what)
{
    std::cerr << "check_raster: " << what << "\n";
    std::exit(1);
}

std::vector<double> numbers(const std::string& line)
{
    std::vector<double> values;
    std::stringstream stream(line);
    std::string cell;
    while (std::getline(stream, cell, ','))
    {
        if (cell.find_first_not_of(" \t\r") == std::string::npos)
            continue;
        values.push_back(std::strtod(cell.c_str(), nullptr));
    }
    return values;
}

struct Grid
{
    std::vector<double> x, y, levels;
};

// -o's layout: a header line, the x axis, the y axis, then a row per y.
Grid read_grid(const std::string& path)
{
    std::ifstream file(path);
    if (!file)
        fail("cannot open " + path);
    std::string line;
    std::getline(file, line);
    Grid grid;
    std::getline(file, line);
    grid.x = numbers(line);
    std::getline(file, line);
    grid.y = numbers(line);
    while (std::getline(file, line))
    {
        const auto row = numbers(line);
        if (row.empty())
            continue;
        if (row.size() != grid.x.size())
            fail(path + ": a row has " + std::to_string(row.size()) + " levels, not " + std::to_string(grid.x.size()));
        grid.levels.insert(grid.levels.end(), row.begin(), row.end());
    }
    if (grid.levels.size() != grid.x.size() * grid.y.size())
        fail(path + ": not one row per y");
    return grid;
}

void check(int status, const std::string& what)
{
    if (status != NC_NOERR)
        fail(what + ": " + nc_strerror(status));
}

std::vector<double> variable(int ncid, const std::string& name, std::size_t size)
{
    int var = -1;
    check(nc_inq_varid(ncid, name.c_str(), &var), "no variable " + name);
    std::vector<double> values(size);
    check(nc_get_var_double(ncid, var, values.data()), "reading " + name);
    return values;
}

std::string text_attribute(int ncid, int var, const char* name)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, var, name, &length) != NC_NOERR)
        return {};
    std::string value(length, '\0');
    check(nc_get_att_text(ncid, var, name, value.data()), std::string("reading attribute ") + name);
    return value;
}

bool close(double a, double b, double tolerance)
{
    return std::abs(a - b) <= tolerance * std::max(1.0, std::abs(b));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 4)
        fail("usage: check_raster RASTER.nc UNIT_TO_METERS [--latlon LATLON.csv | --no-latlon] VAR=GRID.csv ...");
    const std::string path = argv[1];
    const double to_meters = std::strtod(argv[2], nullptr);
    int ncid = -1;
    check(nc_open(path.c_str(), NC_NOWRITE, &ncid), "opening " + path);
    if (text_attribute(ncid, NC_GLOBAL, "Conventions") != "CF-1.8")
        fail(path + " does not say Conventions = CF-1.8");
    int x_dim = -1, y_dim = -1;
    std::size_t nx = 0, ny = 0;
    check(nc_inq_dimid(ncid, "x", &x_dim), "no x dimension");
    check(nc_inq_dimid(ncid, "y", &y_dim), "no y dimension");
    check(nc_inq_dimlen(ncid, x_dim, &nx), "x's length");
    check(nc_inq_dimlen(ncid, y_dim, &ny), "y's length");
    const auto x = variable(ncid, "x", nx);
    const auto y = variable(ncid, "y", ny);

    std::string latlon;
    bool no_latlon = false;
    int checked = 0;
    for (int k = 3; k < argc; k++)
    {
        const std::string arg = argv[k];
        if (arg == "--latlon" && k + 1 < argc)
        {
            latlon = argv[++k];
            continue;
        }
        if (arg == "--no-latlon")
        {
            no_latlon = true;
            continue;
        }
        const auto equals = arg.find('=');
        if (equals == std::string::npos)
            fail("not VAR=GRID.csv: " + arg);
        const std::string name = arg.substr(0, equals);
        const auto grid = read_grid(arg.substr(equals + 1));
        if (grid.x.size() != nx || grid.y.size() != ny)
            fail(name + ": the grid file is " + std::to_string(grid.x.size()) + " by " + std::to_string(grid.y.size())
                    + ", the raster " + std::to_string(nx) + " by " + std::to_string(ny));
        for (std::size_t i = 0; i < nx; i++)
            if (!close(x[i], grid.x[i] * to_meters, 1e-12))
                fail("x[" + std::to_string(i) + "] is not the grid file's");
        for (std::size_t i = 0; i < ny; i++)
            if (!close(y[i], grid.y[i] * to_meters, 1e-12))
                fail("y[" + std::to_string(i) + "] is not the grid file's");
        const auto values = variable(ncid, name, nx * ny);
        for (std::size_t i = 0; i < values.size(); i++)
        {
            const bool same = std::isfinite(grid.levels[i]) ? close(values[i], grid.levels[i], 1e-12)
                                                            : std::isnan(values[i]);
            if (!same)
            {
                std::ostringstream message;
                message << name << "[" << i << "] is " << values[i] << ", the grid file " << grid.levels[i];
                fail(message.str());
            }
        }
        int var = -1;
        check(nc_inq_varid(ncid, name.c_str(), &var), name);
        if (!latlon.empty() && text_attribute(ncid, var, "coordinates") != "latitude longitude")
            fail(name + " does not name latitude and longitude as its coordinates");
        checked++;
    }
    if (checked == 0)
        fail("no VAR=GRID.csv to check");

    if (!latlon.empty())
    {
        const auto latitude = variable(ncid, "latitude", nx * ny);
        const auto longitude = variable(ncid, "longitude", nx * ny);
        std::ifstream file(latlon);
        if (!file)
            fail("cannot open " + latlon);
        std::string line;
        std::getline(file, line);
        std::size_t i = 0;
        while (std::getline(file, line))
        {
            const auto row = numbers(line);
            if (row.empty())
                continue;
            if (i >= nx * ny || !close(latitude[i], row[0], 1e-12) || !close(longitude[i], row[1], 1e-12))
                fail("observer " + std::to_string(i) + "'s latitude and longitude are not the --latlon file's");
            i++;
        }
        if (i != nx * ny)
            fail("the --latlon file has " + std::to_string(i) + " observers, the raster " + std::to_string(nx * ny));
    }
    int var = -1;
    if (no_latlon && (nc_inq_varid(ncid, "latitude", &var) == NC_NOERR || nc_inq_varid(ncid, "longitude", &var) == NC_NOERR))
        fail(path + " has latitude or longitude, but the run had no frame to place the grid with");
    nc_close(ncid);
    std::cout << "check_raster: " << checked << " fields match their grid files\n";
    return 0;
}
