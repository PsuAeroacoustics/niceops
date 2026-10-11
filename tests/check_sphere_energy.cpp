// Checks what --provenance_spheres wrote against the --provenance_file of the
// same run (CMakeLists.txt, provenance_*):
//
//   check_sphere_energy SPHERES.nc PROVENANCE.csv UNIT_TO_METERS
//
//   * one observer per provenance row, at the same position, with the same
//     exposure_db;
//   * each observer's shares sum to one (to float rounding), or are all NaN
//     for an observer that heard nothing;
//   * the sphere with the largest share is the row's top_sphere, with its
//     top_sphere_share, and its stored condition is the row's.
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
#include <type_traits>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& what)
{
    std::cerr << "check_sphere_energy: " << what << "\n";
    std::exit(1);
}

void check(int status, const std::string& what)
{
    if (status != NC_NOERR)
        fail(what + ": " + nc_strerror(status));
}

std::size_t dimension(int ncid, const char* name)
{
    int dim = -1;
    std::size_t length = 0;
    check(nc_inq_dimid(ncid, name, &dim), std::string("no dimension ") + name);
    check(nc_inq_dimlen(ncid, dim, &length), name);
    return length;
}

template <typename T>
std::vector<T> variable(int ncid, const char* name, std::size_t size)
{
    int var = -1;
    check(nc_inq_varid(ncid, name, &var), std::string("no variable ") + name);
    std::vector<T> values(size);
    if constexpr (std::is_same_v<T, float>)
        check(nc_get_var_float(ncid, var, values.data()), name);
    else
        check(nc_get_var_double(ncid, var, values.data()), name);
    return values;
}

bool close(double a, double b, double tolerance)
{
    return std::abs(a - b) <= tolerance * std::max(1.0, std::abs(b));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4)
        fail("usage: check_sphere_energy SPHERES.nc PROVENANCE.csv UNIT_TO_METERS");
    const double to_meters = std::strtod(argv[3], nullptr);
    int ncid = -1;
    check(nc_open(argv[1], NC_NOWRITE, &ncid), std::string("opening ") + argv[1]);
    const std::size_t observers = dimension(ncid, "observer");
    const std::size_t spheres = dimension(ncid, "sphere");
    const auto x = variable<double>(ncid, "x", observers);
    const auto y = variable<double>(ncid, "y", observers);
    const auto exposure = variable<double>(ncid, "exposure_db", observers);
    const auto mu = variable<double>(ncid, "advance_ratio", spheres);
    const auto gamma = variable<double>(ncid, "flight_path_angle", spheres);
    const auto ct = variable<double>(ncid, "thrust_coefficient", spheres);
    const auto share = variable<float>(ncid, "energy_share", observers * spheres);
    nc_close(ncid);

    std::ifstream file(argv[2]);
    if (!file)
        fail(std::string("cannot open ") + argv[2]);
    std::string line;
    std::getline(file, line);
    std::size_t i = 0, heard = 0;
    while (std::getline(file, line))
    {
        std::vector<double> row;
        std::stringstream stream(line);
        std::string cell;
        while (std::getline(stream, cell, ','))
            row.push_back(std::strtod(cell.c_str(), nullptr));
        if (row.size() != 13)
            fail("a provenance row without 13 columns");
        if (i >= observers)
            fail("more provenance rows than observers");
        const std::string where = "observer " + std::to_string(i);
        if (!close(x[i], row[0] * to_meters, 1e-9) || !close(y[i], row[1] * to_meters, 1e-9))
            fail(where + " is not at the provenance row's position");
        if (!(std::isfinite(row[2]) ? close(exposure[i], row[2], 1e-9) : !std::isfinite(exposure[i])))
            fail(where + "'s exposure_db is not the provenance row's");
        const float* shares = &share[i * spheres];
        if (std::isnan(row[8]))
        {
            for (std::size_t k = 0; k < spheres; k++)
                if (!std::isnan(shares[k]))
                    fail(where + " heard nothing, but has a share");
            i++;
            continue;
        }
        double sum = 0.0;
        std::size_t top = 0;
        for (std::size_t k = 0; k < spheres; k++)
        {
            if (!(shares[k] >= 0.0f))
                fail(where + " has a negative or NaN share");
            sum += shares[k];
            if (shares[k] > shares[top])
                top = k;
        }
        if (!close(sum, 1.0, 1e-5))
            fail(where + "'s shares sum to " + std::to_string(sum));
        if (static_cast<double>(top) != row[8] || !close(shares[top], row[9], 1e-6))
            fail(where + "'s largest share is sphere " + std::to_string(top) + ", the provenance row's "
                    + std::to_string(row[8]));
        if (!close(mu[top], row[10], 1e-9) || !close(gamma[top], row[11], 1e-9) || !close(ct[top], row[12], 1e-9))
            fail(where + "'s top sphere's condition is not the provenance row's");
        heard++;
        i++;
    }
    if (i != observers)
        fail("the provenance file has " + std::to_string(i) + " rows, the netCDF " + std::to_string(observers)
                + " observers");
    if (heard == 0)
        fail("no observer heard anything: the check is empty");
    std::cout << "check_sphere_energy: " << heard << " of " << observers << " observers match, over " << spheres
              << " spheres\n";
    return 0;
}
