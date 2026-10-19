#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <vector>

#include "canyon_wave.hpp"
#include "gpu/band_energy.hpp"
#include "gpu/gpu.hpp"
#include "gpu/point_rays.hpp"
#include "point_paths.hpp"
#include "urban.hpp"
#include "wedge.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

gpu::Device& device_or_skip()
{
    gpu::Device* device = gpu::Device::get();
    if (!device)
        SKIP("no GPU backend or device (backend " << gpu::Device::backend() << ")");
    return *device;
}

// Every pair, no cutoff: the quantity both of the GPU's ways approximate.
double every_pair(const std::vector<std::complex<double>>& a, const std::vector<double>& L, double k1, double k2)
{
    const double kmid = 0.5 * (k1 + k2), h = 0.5 * (k2 - k1);
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); i++)
    {
        sum += std::norm(a[i]);
        for (std::size_t j = i + 1; j < a.size(); j++)
        {
            const double d = L[j] - L[i], x = h * d;
            sum += 2.0 * std::real(a[i] * std::conj(a[j]) * std::polar(1.0, -kmid * d)) * (x < 1e-8 ? 1.0 : std::sin(x) / x);
        }
    }
    return sum;
}

} // namespace

TEST_CASE("the GPU's kernel sources are embedded and compile", "[gpu]")
{
    for (const char* name : {"prelude_metal.h", "common.h", "band_energy_types.h", "band_energy.kern"})
        REQUIRE(gpu::kernel_source(name) != nullptr);
    REQUIRE(gpu::kernel_source("no_such.kern") == nullptr);
    gpu::Device& device = device_or_skip();
    for (const char* entry : {"band_pair_sum", "band_quadrature"})
        REQUIRE(device.kernel({"band_energy_types.h", "band_energy.kern"}, entry).max_threads_per_group() >= 256);
    REQUIRE_THROWS(device.kernel({"band_energy.kern"}, "band_pair_sum")); // its types are missing
}

// Ray sets like a street microphone's (hundreds of rays over a few hundred
// metres of path, amplitudes falling as 1/L with random phases) across the
// third-octave sub-bands from 50 Hz to 10 kHz: each problem the GPU sums pair
// by pair is band_averaged_energy's, and each it integrates is every pair's.
TEST_CASE("the GPU's band energies are the CPU's", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    std::mt19937 random(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double third = std::pow(2.0, 1.0 / 6.0);
    const std::vector<double> frequencies {50.0, 125.0, 315.0, 800.0, 2000.0, 5000.0, 10000.0};
    std::vector<gpu::RaySet> sets(4);
    std::vector<std::vector<std::complex<double>>> amplitudes(sets.size());
    gpu::BandEnergies energies(device);
    for (std::size_t s = 0; s < sets.size(); s++)
    {
        gpu::RaySet& set = sets[s];
        const std::size_t n = 200 + 300 * s;
        set.length.resize(n);
        // Lengths out to 2 km, where k L reaches 4e5 at 10 kHz: the phases'
        // float-float arithmetic is what holds there.
        for (double& l : set.length)
            l = 120.0 + 400.0 * std::pow(unit(random), 2.0) + (s == 3 ? 1500.0 : 0.0);
        std::sort(set.length.begin(), set.length.end());
        auto& a = amplitudes[s];
        a.resize(n);
        for (std::size_t i = 0; i < n; i++)
            a[i] = std::polar((0.2 + unit(random)) / set.length[i], 2.0 * pi * unit(random));
        for (double f : frequencies)
        {
            const double k1 = 2.0 * pi * f / third / 343.2, k2 = k1 * std::pow(third * third, 1.0 / 3.0);
            set.add(a.data(), k1, k2);
        }
        set.plan();
        energies.add(set);
    }
    energies.submit();
    const std::vector<double>& gpu = energies.wait();
    REQUIRE(gpu.size() == sets.size() * frequencies.size());
    REQUIRE(energies.integrals() > 0);
    REQUIRE(energies.integrals() < gpu.size());
    double worst_pairs = 0.0, worst_integral = 0.0;
    for (std::size_t s = 0; s < sets.size(); s++)
        for (std::size_t p = 0; p < frequencies.size(); p++)
        {
            const gpu::RaySet& set = sets[s];
            // The CPU from the same float amplitudes the GPU was given.
            std::vector<std::complex<double>> a(set.amplitude.begin() + p * set.length.size(),
                    set.amplitude.begin() + (p + 1) * set.length.size());
            const bool integral = set.nodes[p] != 0;
            const double reference = integral ? every_pair(a, set.length, set.k1[p], set.k2[p])
                                              : band_averaged_energy(a, set.length, set.k1[p], set.k2[p]);
            double& worst = integral ? worst_integral : worst_pairs;
            worst = std::max(worst, std::abs(10.0 * std::log10(gpu[s * frequencies.size() + p] / reference)));
        }
    INFO("pairs " << worst_pairs << " dB, integrals " << worst_integral << " dB");
    REQUIRE(worst_pairs < 1e-3);
    REQUIRE(worst_integral < 1e-3);
}

namespace {

// Runs a check kernel of wedge_check.kern over `args` (stride floats each)
// and returns its complex results.
std::vector<std::complex<double>> run_check(gpu::Device& device, const char* kernel, const std::vector<float>& args,
        std::size_t stride)
{
    const std::size_t count = args.size() / stride;
    gpu::Buffer in = device.buffer(args.size() * sizeof(float)), out = device.buffer(2 * count * sizeof(float));
    std::copy(args.begin(), args.end(), in.as<float>());
    gpu::Batch batch = device.batch();
    batch.dispatch(device.kernel({"wedge.h", "wedge_check.kern"}, kernel),
            {gpu::Arg(in), gpu::Arg(out), gpu::Arg::value(static_cast<uint32_t>(count))}, (count + 63) / 64, 64);
    batch.wait();
    std::vector<std::complex<double>> result(count);
    for (std::size_t i = 0; i < count; i++)
        result[i] = {out.as<float>()[2 * i], out.as<float>()[2 * i + 1]};
    return result;
}

double relative(std::complex<double> a, std::complex<double> b)
{
    return std::abs(a - b) / std::max(std::abs(b), 1e-30);
}

} // namespace

// The kernels' UTD transition function and wedge coefficient in float
// against the host's in double, from x = 1e-8 to 1e4, and over wedges of
// n = 1 to 2 at random angles and at angles within a hair of a shadow or
// reflection boundary, where a cotangent meets its pole.
TEST_CASE("the GPU's wedge coefficients are the CPU's", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    std::vector<float> x;
    for (double e = -8.0; e <= 4.0; e += 0.01)
        x.push_back(static_cast<float>(std::pow(10.0, e)));
    for (double v : {49.999, 50.0, 50.001})
        x.push_back(static_cast<float>(v));
    const std::vector<std::complex<double>> f = run_check(device, "check_utd_transition", x, 1);
    double worst_f = 0.0;
    for (std::size_t i = 0; i < x.size(); i++)
        worst_f = std::max(worst_f, relative(f[i], utd_transition(x[i])));

    std::mt19937 random(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<float> args;
    std::vector<std::complex<double>> expected;
    for (int i = 0; i < 4000; i++)
    {
        const double n = 1.0 + unit(random);
        const double theta0 = 0.02 + (n * pi - 0.04) * unit(random);
        double theta = 0.02 + (n * pi - 0.04) * unit(random);
        // A quarter of them a hair from the incident shadow boundary.
        if (i % 4 == 0 && theta0 + pi < n * pi)
            theta = theta0 + pi + (unit(random) - 0.5) * 2e-6;
        const double k = 0.5 * std::pow(400.0, unit(random)), r = 0.1 * std::pow(1e4, unit(random));
        const WedgeAngles a = wedge_angles(n, theta, theta0);
        for (double v : a.cot)
            args.push_back(static_cast<float>(v));
        for (double v : a.c)
            args.push_back(static_cast<float>(v));
        for (double v : {n, k, r})
            args.push_back(static_cast<float>(v));
        expected.push_back(hard_wedge_coefficient(a, k, r));
    }
    const std::vector<std::complex<double>> d = run_check(device, "check_wedge_coefficient", args, 11);
    double worst_d = 0.0;
    for (std::size_t i = 0; i < expected.size(); i++)
        worst_d = std::max(worst_d, relative(d[i], expected[i]));
    INFO("F " << worst_f << ", D " << worst_d << " relative");
    REQUIRE(worst_f < 2e-6);
    REQUIRE(worst_d < 1e-4);
}

// Filon's panel in float against the same integral in double, by its exact
// moments, over panels from micrometres to tens of metres of path, either
// side of the series' switch at a radian of phase.
TEST_CASE("the GPU's Filon panel is the exact integral", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    std::mt19937 random(5);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<float> args;
    std::vector<std::complex<double>> expected;
    for (int i = 0; i < 4000; i++)
    {
        const double k = 0.5 * std::pow(400.0, unit(random));
        const double u0 = 50.0 * unit(random);
        const double d2 = 1e-6 * std::pow(1e7, unit(random)), d1 = d2 * (0.3 + 0.4 * unit(random));
        std::complex<double> g[3];
        for (auto& v : g)
            v = std::polar(0.1 + unit(random), 2.0 * pi * unit(random));
        // The quadratic through (0, g0), (d1, g1), (d2, g2), integrated
        // against exp(i k (u0 + v)) over [0, d2] by its moments in double.
        const std::complex<double> c1 = (g[1] - g[0]) / d1, c2 = ((g[2] - g[0]) / d2 - c1) / (d2 - d1);
        const std::complex<double> ik(0.0, k);
        std::complex<double> m[3];
        if (k * d2 < 0.1)
            for (int j = 0; j < 3; j++)
            {
                std::complex<double> sum = 0.0, term = 1.0;
                for (int nterm = 0; nterm < 20; nterm++)
                {
                    sum += term / static_cast<double>(nterm + j + 1);
                    term *= ik * d2 / static_cast<double>(nterm + 1);
                }
                m[j] = sum * std::pow(d2, j + 1);
            }
        else
        {
            const std::complex<double> e = std::exp(ik * d2);
            m[0] = (e - 1.0) / ik;
            m[1] = (d2 * e - m[0]) / ik;
            m[2] = (d2 * d2 * e - 2.0 * m[1]) / ik;
        }
        expected.push_back(std::exp(ik * u0) * (g[0] * m[0] + c1 * m[1] + c2 * (m[2] - d1 * m[1])));
        for (double v : {k, u0})
        {
            args.push_back(static_cast<float>(v));
            args.push_back(static_cast<float>(v - static_cast<float>(v)));
        }
        args.push_back(static_cast<float>(d1));
        args.push_back(static_cast<float>(d2));
        for (const auto& v : g)
        {
            args.push_back(static_cast<float>(v.real()));
            args.push_back(static_cast<float>(v.imag()));
        }
    }
    const std::vector<std::complex<double>> v = run_check(device, "check_filon_panel", args, 12);
    double worst = 0.0;
    for (std::size_t i = 0; i < expected.size(); i++)
        // Against the panel's size: a panel's integral can cancel to nothing.
        worst = std::max(worst, std::abs(v[i] - expected[i]) / (std::abs(expected[i]) + 1e-3 * args[12 * i + 5]));
    INFO("worst " << worst);
    REQUIRE(worst < 1e-4);
}

namespace {

Footprint block(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

} // namespace

// The GPU forming the rays' amplitudes itself (point_rays.kern: the wedge
// coefficients, and each shared face's aperture from geometry worked out once
// for every sub-band) against the same batch given the CPU's amplitudes,
// band by band, among DigiCity's containers and a district's blocks: the
// energy kernels are the same, so what differs is the amplitudes.  And the
// levels against band_levels on the CPU.
TEST_CASE("the GPU's point-source amplitudes are the CPU's", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    const Buildings containers({block(4.09, 17.66, 6.53, 23.72, 7.77), block(14.81, 17.66, 17.25, 23.72, 7.77),
            block(1.6, 10.82, 4.04, 16.88, 7.77), block(1.6, 8.38, 13.72, 10.82, 5.18),
            block(20.21, 3.86, 22.65, 16.0, 7.77), block(13.5, -1.84, 19.56, 0.6, 7.77)});
    const Buildings blocks({block(0.0, 0.0, 60.0, 60.0, 20.0), block(80.0, 0.0, 140.0, 60.0, 35.0),
            block(0.0, 80.0, 60.0, 140.0, 10.0), block(80.0, 80.0, 140.0, 140.0, 40.0)});
    struct Case
    {
        const Buildings* city;
        Vec3 source, microphone;
    };
    const std::vector<Case> cases {
        {&containers, {10.5, 20.5, 2.0}, {10.5, -1.1, 1.2}},
        {&containers, {10.5, 20.5, 2.0}, {-3.6, 6.4, 1.2}},
        {&containers, {8.0, 12.0, 15.0}, {18.0, 2.0, 1.2}},
        {&blocks, {70.0, 70.0, 30.0}, {30.0, 70.0, 1.5}},
        {&blocks, {70.0, 30.0, 5.0}, {110.0, 150.0, 1.5}},
    };
    const std::vector<double> bands {50, 100, 200, 400, 800, 1600, 3150, 6300, 10000};
    const double c = 343.2;
    // Every other case through the air, so both the GPU's absorption and
    // its absence are held to the CPU's.
    Atmosphere air;
    air.relative_humidity = 50.0;
    const std::vector<double> absorbing = sub_band_absorption(air, bands, c);
    std::vector<std::unique_ptr<PointSourceField>> fields;
    std::vector<gpu::RaySet> host(cases.size());
    std::vector<gpu::PointRaySet> point(cases.size());
    std::vector<std::vector<double>> cpu(cases.size(), std::vector<double>(bands.size()));
    gpu::BandEnergies energies(device);
    std::vector<std::size_t> host_first, point_first;
    std::size_t shared = 0;
    for (std::size_t i = 0; i < cases.size(); i++)
    {
        const Buildings& city = *cases[i].city;
        std::vector<double> reflection(city.facades().size() * bands.size(), 0.9);
        fields.push_back(std::make_unique<PointSourceField>(city, cases[i].source, cases[i].microphone));
        const PointSourceField& field = *fields.back();
        const std::vector<double> absorption = i % 2 ? absorbing : std::vector<double> {};
        field.band_levels(reflection, 1.0, bands, c, cpu[i].data(), absorption);
        field.sub_band_problems(
                reflection, 1.0, bands, c,
                [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                        const std::vector<double>& length, const SubBand& s) {
                    if (band == 0 && part == 0)
                        host[i].length = length;
                    host[i].add(amplitude.data(), s.k1, s.k2);
                },
                absorption);
        host[i].plan();
        gpu::pack_point_rays(field, reflection, 1.0, bands, c, point[i], absorption);
        shared += point[i].face_rays.size();
        host_first.push_back(energies.add(host[i]));
        point_first.push_back(energies.add(point[i]));
    }
    REQUIRE(shared > 0);
    energies.submit();
    const std::vector<double>& energy = energies.wait();
    double worst_amplitudes = 0.0, worst_cpu = 0.0;
    for (std::size_t i = 0; i < cases.size(); i++)
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const auto parts = sub_bands(bands[b], c);
            const double from_host = band_level(bands[b], parts, &energy[host_first[i] + 3 * b]);
            const double from_point = band_level(bands[b], parts, &energy[point_first[i] + 3 * b]);
            worst_amplitudes = std::max(worst_amplitudes, std::abs(from_point - from_host));
            worst_cpu = std::max(worst_cpu, std::abs(from_point - cpu[i][b]));
        }
    INFO(shared << " rays over shared faces; GPU amplitudes " << worst_amplitudes << " dB from the CPU's, levels "
                << worst_cpu << " dB from band_levels");
    REQUIRE(worst_amplitudes < 0.02);
    REQUIRE(worst_cpu < 0.1);
}

// The canyon wave solver marched on the GPU, in float, against the CPU's
// march in double: the band ratios the tables take, over a canyon with
// impedance walls of different absorption, and over open ground.
TEST_CASE("the GPU's canyon march is the CPU's", "[gpu]")
{
    device_or_skip();
    CanyonWaveSetup s;
    s.width = 12.0;
    s.depth = 10.0;
    s.offset = 3.0;
    s.microphone_height = 1.5;
    s.max_frequency = 250.0;
    s.points_per_wavelength = 12.0;
    s.roof_margin = 4.0;
    s.sky_margin = 4.0;
    s.max_duration = 1.0;
    s.wall_absorption = 0.05;
    s.wall_absorption_b = 0.2;
    const ApertureSpectrum cpu = solve_canyon(s);
    s.gpu = true;
    const ApertureSpectrum gpu = solve_canyon(s);
    REQUIRE(gpu.bins == cpu.bins);
    double worst = 0.0;
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (double f : {25.0, 50.0, 100.0, 160.0, 200.0})
        for (double phi : {0.2, 0.7, 1.2, 1.5, 2.0, 2.8})
            worst = std::max(worst, std::abs(10.0 * std::log10(gpu.band_ratio(f / half, f * half, phi)
                                                    / cpu.band_ratio(f / half, f * half, phi))));
    s.gpu = false;
    const ApertureSpectrum open_cpu = solve_open_ground(s, 6.0, 30.0);
    s.gpu = true;
    const ApertureSpectrum open_gpu = solve_open_ground(s, 6.0, 30.0);
    double worst_open = 0.0;
    for (double f : {50.0, 100.0, 200.0})
        for (double phi : {0.8, 1.2, 1.57})
            worst_open = std::max(worst_open, std::abs(10.0 * std::log10(open_gpu.gain_at(f, phi) / open_cpu.gain_at(f, phi))));
    INFO("canyon " << worst << " dB, open ground " << worst_open << " dB; steps " << cpu.steps << " and " << gpu.steps);
    REQUIRE(worst < 1e-3);
    REQUIRE(worst_open < 1e-3);
}
