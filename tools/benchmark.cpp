// The Attenborough et al. benchmark cases for outdoor sound propagation
// (JASA 97, 173-191 (1995)), by the parabolic equation and, where one exists,
// the exact solution.
//
//   niceops_benchmark [--out figures/benchmark] [--cases 1 2 3 4] > log
//
// Writes case<N>_200m.csv and case<N>_10km.csv, each with the range and, per
// frequency, pe_<f> and exact_<f>: 20 log10 |p / p_1m|, the level re the
// direct field at 1 m, which is what their figures 11-14 plot as transmission
// loss.  An exact column is "nan" where there is no exact solution (case 4),
// or where its series cannot be summed (near the source).
//
// The benchmark, as their section I defines it:
//   source 5 m, receiver 1 m, 10, 100 and 1000 Hz, out to 200 m and to 10 km;
//   c0 = 343.23 m/s; locally reacting ground with the impedances they give
//   for their four-parameter model at each frequency (Table I, sec. I C);
//   absorption at 20 C and 70 % RH: the values they quote, 0.18 dB/km at
//   100 Hz and 5.4 dB/km at 1 kHz, and ISO 9613-1's at 10 Hz, which they
//   call negligible (it is 0.002 dB/km);
//   case 1  uniform;
//   case 2  c rising 0.1 /s (downward refraction);
//   case 3  c falling 0.1 /s (upward refraction);
//   case 4  rising 0.1 /s to 100 m, falling 0.1 /s to 300 m, uniform above.
// The profiles are linear with no stated top.  Here they stop at 1000 m,
// the height their contour plots reach, and are uniform above; that matters
// only to case 2 at 10 km, whose fine structure they note depends on where a
// calculation truncates the profile.
//
// The PE's domain is 400 m high to 200 m and 1000 m to 10 km, as their
// finite-difference PE's was (500 m at 1 kHz, which is below -80 dB before any
// path that needs more height arrives).  The exact solutions are for
// c0 / sqrt(1 -+ 2 z / R), R = c0 / 0.1 s: the forms with Airy-function modes,
// the same near the ground (exact_solutions.hpp).

#include <CLI/CLI.hpp>

#include <cmath>
#include <complex>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "exact_solutions.hpp"
#include "pe.hpp"
#include "util/constants.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT

namespace {

constexpr double c0 = 343.23;
constexpr double gradient = 0.1; // 1/s
constexpr double source_height = 5.0;
constexpr double receiver_height = 1.0;
constexpr double profile_top = 1000.0;

struct Frequency
{
    double f;
    std::complex<double> impedance; // Z, exp(-i omega t)
    double absorption_db_per_km;
};

// Their section I C, and the absorption their section IV C 1 quotes.
std::vector<Frequency> frequencies()
{
    Atmosphere air;
    air.temperature = 293.15;
    air.relative_humidity = 70.0;
    return {{10.0, {38.79, 38.41}, 1000.0 * air.attenuation_coefficient(10.0)},
            {100.0, {12.81, 11.62}, 0.18},
            {1000.0, {5.96, 2.46}, 5.4}};
}

struct Span
{
    std::string name;
    double max_range;
    double step;
    double domain_height;
};

void profile(int which, std::vector<double>& z, std::vector<double>& c)
{
    switch (which)
    {
    case 1:
        z = {0.0};
        c = {c0};
        break;
    case 2:
        z = {0.0, profile_top};
        c = {c0, c0 + gradient * profile_top};
        break;
    case 3:
        z = {0.0, profile_top};
        c = {c0, c0 - gradient * profile_top};
        break;
    default:
        z = {0.0, 100.0, 300.0};
        c = {c0, c0 + gradient * 100.0, c0 - gradient * 100.0};
        break;
    }
}

struct Column
{
    std::vector<double> pe, exact;
};

Column run(int which, const Span& span, const Frequency& f, const std::vector<double>& ranges)
{
    const double wavelength = c0 / f.f;
    PeSetup setup;
    setup.frequency = f.f;
    setup.source_height = source_height;
    setup.receiver_heights = {receiver_height};
    setup.max_range = span.max_range;
    setup.output_step = span.step;
    setup.admittance = 1.0 / f.impedance;
    setup.absorption = f.absorption_db_per_km / 1000.0 / (20.0 / std::log(10.0));
    profile(which, setup.z, setup.c);
    // At 10 Hz ten points a wavelength is 3.4 m, coarse beside a 5 m source.
    setup.points_per_wavelength = std::max(10.0, wavelength / 0.5);
    setup.domain_height = f.f >= 1000.0 && span.max_range > 1000.0 ? 500.0 : span.domain_height;
    const auto pe = solve_pe(setup);

    Column column;
    column.pe.resize(ranges.size());
    for (std::size_t q = 0; q < ranges.size(); q++)
    {
        const double slant = std::hypot(ranges[q], receiver_height - source_height);
        column.pe[q] = pe.relative_level[0][q] - 20.0 * std::log10(slant);
    }

    ExactProblem problem;
    problem.frequency = f.f;
    problem.sound_speed = c0;
    problem.absorption = setup.absorption;
    problem.admittance = setup.admittance;
    problem.source_height = source_height;
    problem.receiver_height = receiver_height;
    const double radius = c0 / gradient;
    std::vector<std::complex<double>> exact(ranges.size(),
            std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0.0));
    if (which == 1)
        for (std::size_t q = 0; q < ranges.size(); q++)
            exact[q] = homogeneous_pressure(problem, ranges[q]);
    else if (which == 2)
        exact = downward_pressure(problem, radius, ranges);
    else if (which == 3)
        exact = upward_pressure(problem, radius, ranges);
    column.exact.resize(ranges.size());
    for (std::size_t q = 0; q < ranges.size(); q++)
        column.exact[q] = 20.0 * std::log10(std::abs(exact[q]));
    return column;
}

} // namespace

int main(int argc, char** argv)
{
    CLI::App app {"The Attenborough et al. (1995) benchmark cases, by the PE and exact solutions."};
    std::string out = "figures/benchmark";
    std::vector<int> cases = {1, 2, 3, 4};
    app.add_option("--out", out, "Output directory (default figures/benchmark).");
    app.add_option("--cases", cases, "Cases to run, 1-4 (default all).")->check(CLI::Range(1, 4));
    CLI11_PARSE(app, argc, argv);

    const std::vector<Span> spans = {{"200m", 200.0, 0.5, 400.0}, {"10km", 10000.0, 5.0, 1000.0}};
    const auto bands = frequencies();

    for (const auto& span : spans)
    {
        const std::size_t count = static_cast<std::size_t>(std::llround(span.max_range / span.step));
        std::vector<double> ranges(count);
        for (std::size_t q = 0; q < count; q++)
            ranges[q] = static_cast<double>(q + 1) * span.step;
        for (const int which : cases)
        {
            std::vector<std::future<Column>> jobs;
            for (const auto& f : bands)
                jobs.push_back(std::async(std::launch::async, run, which, span, f, ranges));
            std::vector<Column> columns;
            for (auto& job : jobs)
                columns.push_back(job.get());

            const std::string path = out + "/case" + std::to_string(which) + "_" + span.name + ".csv";
            ensure_parent_directory(path);
            std::ofstream file(path);
            file << "range_m";
            for (const auto& f : bands)
                file << ",pe_" << f.f << ",exact_" << f.f;
            file << "\n" << std::setprecision(7);
            for (std::size_t q = 0; q < count; q++)
            {
                file << ranges[q];
                for (const auto& column : columns)
                    file << "," << column.pe[q] << "," << column.exact[q];
                file << "\n";
            }
            std::cout << path << std::endl;
        }
    }
    return 0;
}
