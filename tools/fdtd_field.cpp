// The full-wave field of a point source among rigid buildings at microphones
// (src/fdtd.hpp): the reference niceops_point_field is held to.  The same
// inputs and output as niceops_point_field -- each band's level re the
// source's own at 1 m, |H|^2 averaged over the third octave.
//
//   niceops_fdtd_field --buildings case/buildings.csv --source -2 0 1.5 --microphones mics.csv
//       --bands 1600 2000 2500 --units m --mirror --gpu -o levels.csv
//
// --buildings is niceops_point_field's, standing on a rigid floor at z = 0;
// --solids adds prisms (solid, x, y, bottom, top, a solid's vertices in order
// around it), for geometry that is not a building, such as a facade's relief.
// The domain is the box around the buildings, the source and the
// microphones, `--margin` wider each way (and up), then the PML.  --mirror
// takes the plane y = the source's y as a plane of symmetry -- right only
// where the solids are symmetric about it -- and grids half the space.
//
// The differences are fourth order by default (--order 2 gives Yee's), and the
// grid's cell defaults to 8 a wavelength at the top band's upper edge (12 at
// second order).  Check a run by repeating it at a finer --cell (FIDELITY
// "The full-wave reference among buildings").
// --duration defaults to twice the time sound takes across the domain's
// diagonal, the whole space's under --mirror; a microphone whose band levels move by more than 0.02 dB without
// the recording's last tenth is named on stderr.
//
// The solids are rigid, or absorb --wall_absorption at normal incidence (a
// real impedance; the floor stays rigid).  The air is lossless in the march;
// --humidity absorbs each arrival afterwards by its delay, as
// niceops_point_field absorbs each ray by its length.
//
// It marches on the GPU wherever there is one that holds the grid, and on
// the CPU otherwise (--no-gpu always); --gpu refuses the CPU.
//
// --spectrum writes H across the bands, every quarter of 1 / duration: the
// microphone's row, the frequency, Hz, and H's real and imaginary parts.

#include <CLI/CLI.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "building_input.hpp"
#include "buildings.hpp"
#include "fdtd.hpp"
#include "point_paths.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"
#include "util/parallel.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"The full-wave field of a point source among buildings, at microphones (FDTD)."};
    std::string buildings_file, solids_file, microphones_file, output, spectrum_file, recording_file;
    std::string units = "ft";
    std::vector<double> source, bands;
    double cell = 0.0, points_per_wavelength = 0.0, margin = std::numeric_limits<double>::quiet_NaN();
    int order = 4;
    double duration = 0.0, courant = 0.99, pulse_frequency = 0.0, sound_speed = 343.2;
    double wall_absorption = 0.0, humidity = std::numeric_limits<double>::quiet_NaN();
    double temperature = reference_temperature, pressure = reference_pressure;
    int pml = 20, threads = 0;
    bool mirror = false, use_gpu = true;
    app.add_option("--buildings", buildings_file, "The district: CSV of footprint vertices, or GeoJSON.")->required();
    app.add_option("--solids", solids_file, "More rigid prisms: a CSV of solid, x, y, bottom and top, each "
                                            "solid's vertices in order around it.");
    app.add_option("--source", source, "The source's x, y and height above the floor.")->expected(3)->required();
    app.add_option("--microphones", microphones_file, "The microphones: x, y and height columns.")->required();
    app.add_option("--bands", bands, "Band centres, Hz.")->required();
    app.add_option("--units", units, "Length unit of every input, option and output: ft (default) or m.")
            ->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--cell", cell, "The grid's cell, in --units (default: --points_per_wavelength at the top "
                                   "band's upper edge).");
    app.add_option("--points_per_wavelength", points_per_wavelength, "Cells a wavelength at the top band's upper "
                                                                     "edge, where --cell is not given (default 8 at "
                                                                     "--order 4, 12 at 2).");
    app.add_option("--order", order, "The spatial differences' order: 4 (default, with the frequency warp) or 2.")
            ->check(CLI::IsMember({2, 4}));
    app.add_option("--margin", margin, "Air beyond the buildings, source and microphones each way before the PML, "
                                       "in --units (default a wavelength at the lowest band's lower edge).");
    app.add_option("--pml", pml, "The PML's thickness, cells (default 20).");
    app.add_option("--duration", duration, "Simulated time, s (default twice the domain's diagonal over c).");
    app.add_option("--courant", courant, "The time step, as a fraction of the scheme's limit (default 0.99).");
    app.add_option("--pulse_frequency", pulse_frequency, "Where the pulse's spectrum peaks, Hz (default the top "
                                                         "band's upper edge over 2.25, where it is 10 dB down).");
    app.add_option("--sound_speed", sound_speed, "m/s (default 343.2).");
    app.add_option("--wall_absorption", wall_absorption, "The solids' absorption at normal incidence, a real "
                                                         "impedance at every frequency (default 0, rigid).");
    app.add_option("--humidity", humidity, "Relative humidity, percent: the air's absorption (ISO 9613-1), applied "
                                           "to each arrival by its delay. Default none.");
    app.add_option("--temperature", temperature, "With --humidity, K (default 293.15).");
    app.add_option("--pressure", pressure, "With --humidity, kPa (default 101.325).");
    app.add_flag("--mirror", mirror, "The plane y = the source's y is one of symmetry: grid half the space.");
    app.add_flag("--gpu,!--no-gpu", use_gpu,
            "March on the GPU: on by default wherever there is one that holds the grid, the CPU otherwise. --gpu "
            "asks for it and is refused where there is none or it cannot hold the grid; --no-gpu marches on the CPU.");
    app.add_option("-j,--threads", threads, "CPU threads (default: all).");
    app.add_option("-o,--output", output, "The band levels, CSV.")->required();
    app.add_option("--spectrum", spectrum_file, "H across the bands, every quarter of 1 / duration, CSV.");
    app.add_option("--recording", recording_file, "Each microphone's pressure at every step and its time, and the "
                                                    "pulse Q at its own time, half a step before, CSV.");
    CLI11_PARSE(app, argc, argv);

    try
    {
        const double to_meters = units == "ft" ? ft2m : 1.0;
        std::sort(bands.begin(), bands.end());
        const double half = std::pow(2.0, 1.0 / 6.0);
        const double f_top = bands.back() * half, f_bottom = bands.front() / half;
        if (!(points_per_wavelength > 0.0))
            points_per_wavelength = order == 4 ? 8.0 : 12.0;
        const double h = cell > 0.0 ? cell * to_meters : sound_speed / (f_top * points_per_wavelength);
        const double air = std::isnan(margin) ? sound_speed / f_bottom : margin * to_meters;

        BuildingInputOptions input;
        input.units = units;
        input.to_meters = to_meters;
        const BuildingInput district = read_buildings(buildings_file, input);
        std::vector<FdtdSolid> solids = fdtd_solids(district.footprints);
        if (!solids_file.empty())
        {
            const LabeledCsv rows = read_numeric_rows(solids_file, "solids");
            const std::size_t id = length_column(rows, solids_file, units, "solid", true, false);
            const std::size_t x = length_column(rows, solids_file, units, "x", true, true);
            const std::size_t y = length_column(rows, solids_file, units, "y", true, true);
            const std::size_t bottom = length_column(rows, solids_file, units, "bottom", true, true);
            const std::size_t top = length_column(rows, solids_file, units, "top", true, true);
            double current = std::numeric_limits<double>::quiet_NaN();
            for (const auto& r : rows.rows)
            {
                if (r[id] != current)
                {
                    solids.emplace_back();
                    solids.back().bottom = r[bottom] * to_meters;
                    solids.back().top = r[top] * to_meters;
                    current = r[id];
                }
                solids.back().x.push_back(r[x] * to_meters);
                solids.back().y.push_back(r[y] * to_meters);
            }
        }

        const Vec3 s(source[0] * to_meters, source[1] * to_meters, source[2] * to_meters);
        const LabeledCsv mics = read_numeric_rows(microphones_file, "microphones");
        const std::size_t mx = length_column(mics, microphones_file, units, "x", true, true);
        const std::size_t my = length_column(mics, microphones_file, units, "y", true, true);
        const std::size_t mh = length_column(mics, microphones_file, units, "height", true, true);
        std::vector<Vec3> points;
        for (const auto& r : mics.rows)
            points.emplace_back(r[mx] * to_meters, r[my] * to_meters, r[mh] * to_meters);
        // With a symmetry plane each microphone is taken on the source's side.
        if (mirror)
            for (Vec3& p : points)
                p.y = s.y + std::abs(p.y - s.y);

        // The box around everything.
        Vec3 lo = s, hi = s;
        auto take = [&](double x, double y, double z) {
            lo = Vec3(std::min(lo.x, x), std::min(lo.y, y), std::min(lo.z, z));
            hi = Vec3(std::max(hi.x, x), std::max(hi.y, y), std::max(hi.z, z));
        };
        for (const Vec3& p : points)
            take(p.x, p.y, p.z);
        for (const FdtdSolid& b : solids)
            for (std::size_t i = 0; i < b.x.size(); i++)
                take(b.x[i], b.y[i], std::max(0.0, std::min(b.top, 1e6)));
        FdtdSetup setup;
        setup.cell = h;
        setup.order = order;
        setup.lo = Vec3(lo.x - air, mirror ? s.y : lo.y - air, 0.0);
        setup.hi = Vec3(hi.x + air, hi.y + air, hi.z + air);
        setup.floor = true;
        setup.mirror = mirror;
        setup.pml = pml;
        setup.sound_speed = sound_speed;
        setup.courant = courant;
        // Long enough to cross the whole space twice: with the plane, the
        // half on the source's side and its image.
        Vec3 span = setup.hi - setup.lo;
        if (mirror)
            span.y *= 2.0;
        setup.duration = duration > 0.0 ? duration : 2.0 * span.norm() / sound_speed;
        setup.pulse_frequency = pulse_frequency > 0.0 ? pulse_frequency : f_top / 2.25;
        setup.threads = tool_thread_count(threads);
        // The GPU wherever there is one that holds the grid; asked for, it
        // must be there and hold it.  The flag's count takes --no-gpu too, so
        // asking is its value as well.
        setup.gpu = use_gpu;
        setup.gpu_required = use_gpu && app.count("--gpu") > 0;
        setup.wall_absorption = wall_absorption;
        // The air's absorption per sub-band, as niceops_point_field takes it.
        std::vector<double> absorption;
        if (!std::isnan(humidity))
        {
            for (const std::string& problem :
                    {temperature_problem(temperature), pressure_problem(pressure), humidity_problem(humidity)})
                if (!problem.empty())
                    throw std::runtime_error(problem);
            Atmosphere atmosphere;
            atmosphere.temperature = temperature;
            atmosphere.pressure = pressure;
            atmosphere.relative_humidity = humidity;
            absorption = sub_band_absorption(atmosphere, bands, sound_speed);
        }
        else if (app.count("--temperature") > 0 || app.count("--pressure") > 0)
            throw std::runtime_error("--temperature and --pressure go with --humidity");

        const FdtdRecording run = solve_fdtd(setup, solids, s, points);
        std::cerr << "FDTD: " << run.nx << " x " << run.ny << " x " << run.nz << " cells of " << h * 1000.0
                  << " mm (" << run.cells * 1e-6 << " M), " << run.steps << " steps of " << run.dt * 1e6 << " us, "
                  << run.seconds << " s"
                  << (run.gpu ? " on the GPU" : run.gpu_declined.empty() ? "" : " on the CPU, " + run.gpu_declined)
                  << ", "
                  << std::lround(run.marched * 100.0) << "% of the cell-steps marched" << std::endl;
        ensure_parent_directory(output);
        std::ofstream out(output);
        if (!out)
            throw std::runtime_error("cannot write " + output);
        out << "x_" << units << ",y_" << units << ",height_" << units;
        for (double b : bands)
            out << ",L" << b;
        out << '\n' << std::setprecision(10);
        // Each band's level, and how much the run's end still moved it.
        std::vector<std::vector<double>> levels(points.size(), std::vector<double>(bands.size())),
                moved(points.size(), std::vector<double>(bands.size()));
        parallel_chunks(points.size() * bands.size(), std::max(1u, setup.threads), [&](std::size_t begin, std::size_t end) {
            for (std::size_t q = begin; q < end; q++)
            {
                const std::size_t r = q / bands.size(), b = q % bands.size();
                levels[r][b] = run.band_level(r, bands[b], 0, absorption.empty() ? nullptr : &absorption[3 * b]);
                moved[r][b] = std::abs(run.truncation_db(r, bands[b]));
            }
        });
        for (std::size_t r = 0; r < points.size(); r++)
        {
            const double worst = *std::max_element(moved[r].begin(), moved[r].end());
            if (worst > 0.02)
                std::cerr << "Microphone " << r + 1 << ": the recording's last tenth moves a band by " << worst
                          << " dB; lengthen --duration." << std::endl;
        }
        for (std::size_t r = 0; r < points.size(); r++)
        {
            const auto& row = mics.rows[r];
            out << row[mx] << ',' << row[my] << ',' << row[mh];
            for (double level : levels[r])
                out << ',' << level;
            out << '\n';
        }
        if (!recording_file.empty())
        {
            ensure_parent_directory(recording_file);
            std::ofstream recording(recording_file);
            recording << "t_s,q_t_s,q";
            for (std::size_t r = 0; r < points.size(); r++)
                recording << ",p" << r + 1;
            recording << '\n' << std::setprecision(9);
            for (std::size_t n = 0; n < run.steps; n++)
            {
                recording << (static_cast<double>(n) + 1.0) * run.dt << ',' << (static_cast<double>(n) + 0.5) * run.dt
                          << ',' << run.source[n];
                for (std::size_t r = 0; r < points.size(); r++)
                    recording << ',' << run.pressure[n * points.size() + r];
                recording << '\n';
            }
        }
        if (!spectrum_file.empty())
        {
            ensure_parent_directory(spectrum_file);
            std::ofstream spectrum(spectrum_file);
            spectrum << "microphone,f_hz,re,im\n" << std::setprecision(10);
            const double step = 1.0 / (4.0 * static_cast<double>(run.steps) * run.dt);
            for (std::size_t r = 0; r < points.size(); r++)
                for (double f = f_bottom; f <= f_top; f += step)
                {
                    const std::complex<double> H = run.transfer(r, f);
                    spectrum << r + 1 << ',' << f << ',' << H.real() << ',' << H.imag() << '\n';
                }
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_fdtd_field: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
