#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <complex>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "canyon_wave.hpp"
#include "gpu/band_energy.hpp"
#include "gpu/band_stream.hpp"
#include "gpu/gpu.hpp"
#include "gpu/point_rays.hpp"
#include "facade_relief.hpp"
#include "ground_effect.hpp"
#include "point_paths.hpp"
#include "source_coherence.hpp"
#include "sphere_database.hpp"
#include "urban.hpp"
#include "urban_point.hpp"
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

// Street-like fields: two districts, DigiCity's containers and a district's
// blocks, sources above and among the buildings.
struct Fields
{
    Buildings containers {{block(4.09, 17.66, 6.53, 23.72, 7.77), block(14.81, 17.66, 17.25, 23.72, 7.77),
            block(1.6, 10.82, 4.04, 16.88, 7.77), block(1.6, 8.38, 13.72, 10.82, 5.18),
            block(20.21, 3.86, 22.65, 16.0, 7.77), block(13.5, -1.84, 19.56, 0.6, 7.77)}};
    Buildings blocks {{block(0.0, 0.0, 60.0, 60.0, 20.0), block(80.0, 0.0, 140.0, 60.0, 35.0),
            block(0.0, 80.0, 60.0, 140.0, 10.0), block(80.0, 80.0, 140.0, 140.0, 40.0)}};
    // A mansard house, over which a low source reaches a microphone close
    // behind by four roof edges (a chain: more wedge terms than three).
    Buildings mansard {{[] {
        Footprint f = block(0.0, -30.0, 12.0, 30.0, 9.0);
        f.roof = {{1.0, 0.0, 6.0}, {0.4, 0.0, 6.9}, {-0.4, 0.0, 11.7}, {-1.0, 0.0, 18.0}};
        return f;
    }()}};
    Buildings slopes {{[] {
        Footprint f = block(0.0, 0.0, 10.0, 20.0, 9.0);
        f.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
        return f;
    }(), block(-60.0, 25.0, 60.0, 27.0, 50.0), block(-8.0, -20.0, -6.0, 24.0, 26.0)}};
    std::vector<std::unique_ptr<PointSourceField>> fields;
    std::vector<std::vector<double>> reflection;
    // Per field, whether its microphone lies on a roof (PointRay::
    // microphone_surface), where a plate takes the roof's R.
    std::vector<bool> on_roof;

    explicit Fields(std::size_t bands, const PointPathOptions& options = {})
    {
        const std::vector<std::tuple<const Buildings*, Vec3, Vec3>> cases {
            {&containers, {10.5, 20.5, 2.0}, {10.5, -1.1, 1.2}},
            {&containers, {10.5, 20.5, 2.0}, {-3.6, 6.4, 1.2}},
            {&containers, {8.0, 12.0, 15.0}, {18.0, 2.0, 1.2}},
            {&blocks, {70.0, 70.0, 30.0}, {30.0, 70.0, 1.5}},
            {&blocks, {70.0, 30.0, 5.0}, {110.0, 150.0, 1.5}},
            // A flush microphone, for a plate.
            {&blocks, {70.0, 70.0, 30.0}, {30.0, 70.0, 0.0}},
            {&mansard, {-15.0, 0.0, 0.5}, {14.0, 1.0, 0.5}},
            {&slopes, {-20.0, 15.0, 40.0}, {42.0, 10.0, 1.5}},
            {&slopes, {42.0, 10.0, 1.5}, {-20.0, 15.0, 40.0}},
            // A flush microphone on a roof, behind its edge from a low source.
            {&blocks, {-20.0, 30.0, 5.0}, {30.0, 30.0, 20.0}},
            // Both ends above a roof's plane, the source low over it: the pair
            // takes each end's mirror (double_wedge.cpp, mirrored_weights).
            {&blocks, {-40.0, 30.0, 21.0}, {70.0, 30.0, 20.5}},
        };
        for (const auto& [city, source, microphone] : cases)
        {
            fields.push_back(std::make_unique<PointSourceField>(*city, source, microphone, options));
            reflection.emplace_back(city->facades().size() * bands, 0.9);
            on_roof.push_back(city->support(microphone).kind == Buildings::Support::Kind::roof);
        }
    }
};

} // namespace

// The GPU forming the rays' amplitudes itself (point_rays.kern: the wedge
// coefficients, and each shared face's aperture from geometry worked out once
// for every sub-band) against the same batch given the CPU's amplitudes,
// band by band, among DigiCity's containers and a district's blocks: the
// energy kernels are the same, so what differs is the amplitudes.  And the
// levels against band_levels on the CPU.  With the shared faces' integrals
// at every wavenumber, and sampled (the cascade at every problem, the
// integral at the probes, face_ratios between) as the CPU samples them.
// Every third field's floor is grass, whose bounces take Q, and the flush
// microphone's is on a plate besides.
TEST_CASE("the GPU's point-source amplitudes are the CPU's", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    for (const bool sampled : {false, true})
    {
    INFO((sampled ? "sampled" : "every wavenumber"));
    const std::vector<double> bands {50, 100, 200, 400, 800, 1600, 3150, 6300, 10000};
    PointPathOptions options;
    options.sample_shared_faces = sampled;
    const Fields f(bands.size(), options);
    const std::size_t count = f.fields.size();
    const double c = 343.2;
    // Every other case through the air, so both the GPU's absorption and
    // its absence are held to the CPU's.
    Atmosphere air;
    air.relative_humidity = 50.0;
    const std::vector<double> absorbing = sub_band_absorption(air, bands, c);
    std::vector<gpu::RaySet> host(count);
    std::vector<gpu::PointRaySet> point(count);
    std::vector<std::vector<double>> cpu(count, std::vector<double>(bands.size()));
    gpu::BandEnergies energies(device);
    std::vector<std::size_t> host_first, point_first;
    std::size_t shared = 0, thrice = 0, chained = 0, between = 0, plates = 0, roof_plates = 0, slope_walls = 0,
                mirrored = 0;
    const std::vector<double> middles = sub_band_middles(bands, c);
    FloorModel grass;
    grass.admittance = ground_impedance(GroundModel {}, middles, c);
    FloorModel plate = grass;
    plate.plate = [&](double elevation, double azimuth, std::complex<double>* direct, std::complex<double>* image) {
        for (std::size_t j = 0; j < middles.size(); j++)
        {
            direct[j] = std::polar(1.0 + 0.3 * std::sin(elevation), 0.2 * azimuth + 0.01 * j);
            image[j] = std::polar(0.9, -0.1 * j);
        }
    };
    for (std::size_t i = 0; i < count; i++)
    {
        const std::vector<double>& reflection = f.reflection[i];
        const PointSourceField& field = *f.fields[i];
        const bool flush = field.rays().front().path.back().z == 0.0;
        plates += flush;
        // A plate on the roof takes its R, the root of 0.9, at every sub-band.
        FloorModel on_roof = plate;
        on_roof.surface.assign(middles.size(), std::sqrt(0.9));
        roof_plates += f.on_roof[i];
        const FloorModel floor = flush ? plate : f.on_roof[i] ? on_roof : i % 3 == 1 ? grass : FloorModel(1.0);
        for (const PointRay& ray : field.rays())
        {
            slope_walls += !ray.roof_facades.empty() && !ray.facades.empty();
            mirrored += ray.shared && ray.shared_face.sz > 0.0 && ray.shared_face.rz > 0.0;
            thrice += ray.wedges == 3;
            chained += ray.wedges > 3;
            between += ray.between != PointRay::no_facade;
        }
        const std::vector<double> absorption = i % 2 ? absorbing : std::vector<double> {};
        // No bands is nothing to do, shared and profile rays included.
        REQUIRE_NOTHROW(field.band_levels({}, FloorModel(1.0), {}, c, nullptr));
        field.band_levels(reflection, floor, bands, c, cpu[i].data(), absorption);
        field.sub_band_problems(
                reflection, floor, bands, c,
                [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                        const std::vector<double>& length, const SubBand& s) {
                    if (band == 0 && part == 0)
                        host[i].length = length;
                    host[i].add(amplitude.data(), s.k1, s.k2);
                },
                absorption);
        host[i].plan();
        gpu::pack_point_rays(field, reflection, floor, bands, c, point[i], absorption);
        shared += point[i].face_rays.size() + point[i].samples.size();
        host_first.push_back(energies.add(host[i]));
        point_first.push_back(energies.add(point[i]));
    }
    REQUIRE(shared > 0);
    REQUIRE(mirrored > 0);
    // Rays of three wedge terms and of more, and rays with a facade between
    // two edges, take the same kernel.
    REQUIRE(thrice > 0);
    REQUIRE(chained > 0);
    REQUIRE(between > 0);
    REQUIRE(slope_walls > 0);
    REQUIRE(plates == 1);
    REQUIRE(roof_plates == 1);
    energies.submit();
    const std::vector<double>& energy = energies.wait();
    double worst_amplitudes = 0.0, worst_cpu = 0.0;
    for (std::size_t i = 0; i < count; i++)
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
}


// The cone's and energy's sums on the GPU against partially_coherent_energies
// on the same rays: a narrow cone, a wide one whose taper takes many pairs,
// and no pair at all.
TEST_CASE("the GPU's cone and energy sums are the CPU's", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    const std::vector<double> bands {50, 100, 200, 400, 800, 1600, 3150, 6300, 10000};
    const double c = 343.2;
    Atmosphere air;
    air.relative_humidity = 50.0;
    const std::vector<double> absorption = sub_band_absorption(air, bands, c);
    const Fields f(bands.size());
    struct Sum
    {
        RaySum sum;
        uint32_t gpu;
        double cone;
    };
    for (const Sum& sum : {Sum {RaySum::cone, NK_SUM_CONE, 1.0}, Sum {RaySum::cone, NK_SUM_CONE, 6.0},
                 Sum {RaySum::energy, NK_SUM_ENERGY, 1.0}})
    {
        INFO((sum.sum == RaySum::energy ? "energy" : "cone of " + std::to_string(sum.cone) + " degrees"));
        SourceCoherence coherence;
        coherence.sum = sum.sum;
        coherence.cone = sum.cone;
        gpu::BandEnergies energies(device);
        std::vector<gpu::PointRaySet> sets(f.fields.size());
        std::vector<std::size_t> first;
        std::vector<std::vector<double>> cpu(f.fields.size(), std::vector<double>(3 * bands.size()));
        for (std::size_t i = 0; i < f.fields.size(); i++)
        {
            partially_coherent_energies(*f.fields[i], f.reflection[i], 1.0, bands, 0, c, absorption, nullptr,
                    coherence, SphereFrame {}, cpu[i].data());
            gpu::PointRayPacking packing;
            packing.sum = sum.gpu;
            packing.cone = sum.cone;
            gpu::pack_point_rays(*f.fields[i], f.reflection[i], 1.0, bands, c, sets[i], absorption, packing);
            first.push_back(energies.add(sets[i]));
        }
        energies.submit();
        const std::vector<double>& energy = energies.wait();
        REQUIRE(energies.integrals() == 0);
        double worst = 0.0;
        for (std::size_t i = 0; i < f.fields.size(); i++)
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                const auto parts = sub_bands(bands[b], c);
                for (int part = 0; part < 3; part++)
                    REQUIRE(energy[first[i] + 3 * b + static_cast<std::size_t>(part)] >= 0.0);
                worst = std::max(worst, std::abs(band_level(bands[b], parts, &energy[first[i] + 3 * b])
                                                - band_level(bands[b], parts, &cpu[i][3 * b])));
            }
        INFO("worst " << worst << " dB");
        // 3e-4 dB: the GPU's float.
        REQUIRE(worst < 0.002);
    }
}

// The cone's taper is not a positive-definite kernel: three rays along a
// line of launches, each next to its neighbour (fully coherent) and the ends
// apart (in energy), with amplitudes +1, -1, +1 and the same length, sum to
// 3 - 2 - 2 < 0.  The GPU takes 0 there, as the CPU does.
TEST_CASE("the GPU's cone sum is not below zero", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    const double c = 343.2, band = 200.0;
    gpu::PointRaySet set;
    set.sum = NK_SUM_CONE;
    set.cone = 1.0f;
    const double sign[3] = {1.0, -1.0, 1.0};
    for (uint32_t i = 0; i < 3; i++)
    {
        set.length.push_back(50.0 + 1e-4 * i);
        set.rays.push_back(GpuRay {0u, 0u, i, static_cast<float>(set.length.back())});
        set.reflected.push_back(static_cast<float>(sign[i]));
        const double turn = 0.7 * i * pi / 180.0;
        for (double v : {std::cos(turn), std::sin(turn), 0.0})
            set.launch.push_back(static_cast<float>(v));
    }
    for (const SubBand& s : sub_bands(band, c))
    {
        set.k1.push_back(s.k1);
        set.k2.push_back(s.k2);
        set.absorption.push_back(0.0f);
        set.nodes.push_back(0u);
    }
    gpu::BandEnergies energies(device);
    energies.add(set);
    energies.submit();
    const std::vector<double>& energy = energies.wait();
    REQUIRE(energy.size() == 3);
    for (double e : energy)
        REQUIRE(e == 0.0);
}

// A stream of sets in batches of two, two batches in flight, gives each set
// the energies one batch of them all does, and every set back once.
TEST_CASE("the GPU's stream gives every set its own energies", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    const std::vector<double> bands {100, 400, 1600, 6300};
    const double c = 343.2;
    const Fields f(bands.size());
    const std::size_t copies = 3;
    std::vector<gpu::PointRaySet> sets(f.fields.size() * copies);
    for (std::size_t i = 0; i < sets.size(); i++)
    {
        gpu::PointRayPacking packing;
        packing.sum = i % 2 ? NK_SUM_CONE : NK_SUM_COHERENT;
        gpu::pack_point_rays(*f.fields[i % f.fields.size()], f.reflection[i % f.fields.size()], 1.0, bands, c,
                sets[i], {}, packing);
    }
    gpu::BandEnergies energies(device);
    std::vector<std::size_t> first;
    for (const gpu::PointRaySet& set : sets)
        first.push_back(energies.add(set));
    energies.submit();
    const std::vector<double> together = energies.wait();

    std::vector<std::vector<double>> streamed(sets.size());
    std::vector<int> calls(sets.size(), 0);
    {
        gpu::BandStream stream(device, 2, [&](std::size_t tag, const double* energy) {
            calls[tag]++;
            streamed[tag].assign(energy, energy + sets[tag].problems());
        });
        // Pushed in reverse, tagged by their place.
        for (std::size_t i = sets.size(); i-- > 0;)
        {
            gpu::PointRaySet copy = sets[i];
            stream.push(i, std::move(copy));
        }
        stream.finish();
    }
    for (std::size_t i = 0; i < sets.size(); i++)
    {
        INFO("set " << i);
        REQUIRE(calls[i] == 1);
        REQUIRE(streamed[i].size() == sets[i].problems());
        for (std::size_t p = 0; p < sets[i].problems(); p++)
            REQUIRE(streamed[i][p] == together[first[i] + p]);
    }
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

// A relief's specular share on the GPU: the rays' reflections complex per
// problem, as a floor's impedance makes them, the CPU's factors taken over.
TEST_CASE("the GPU takes a relief's specular share as the CPU does", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    const std::string path = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") + "/niceops_gpu_relief.csv";
    {
        std::ofstream out(path);
        out << "frequency_hz,incidence_deg,order,share,period_m\n";
        for (double f : {50.0, 20000.0})
            for (int psi = -85; psi <= 85; psi += 5)
                for (const auto& [m, share] : std::vector<std::pair<int, double>> {{-1, 0.3}, {0, 0.4}, {1, 0.3}})
                    out << f << "," << psi << "," << m << "," << share << ",0.3\n";
    }
    const FacadeRelief relief = FacadeRelief::read(path);
    const Buildings blocks {{block(0.0, 0.0, 60.0, 60.0, 20.0), block(80.0, 0.0, 140.0, 60.0, 35.0),
            block(0.0, 80.0, 60.0, 140.0, 10.0), block(80.0, 80.0, 140.0, 140.0, 40.0)}};
    std::vector<const FacadeRelief*> reliefs(blocks.facades().size(), &relief);
    PointPathOptions options;
    options.reliefs = &reliefs;
    const PointSourceField field(blocks, Vec3(70.0, 70.0, 30.0), Vec3(30.0, 70.0, 1.5), options);
    const std::vector<double> bands {250, 1000, 4000};
    const std::vector<double> reflection(blocks.facades().size() * bands.size(), 0.9);
    std::vector<double> cpu(bands.size());
    field.band_levels(reflection, 1.0, bands, 343.2, cpu.data());
    gpu::PointRaySet set;
    gpu::pack_point_rays(field, reflection, 1.0, bands, 343.2, set);
    REQUIRE(set.complex_reflections);
    gpu::BandEnergies energies(device);
    const std::size_t first = energies.add(set);
    energies.submit();
    const std::vector<double>& energy = energies.wait();
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(std::abs(band_level(bands[b], sub_bands(bands[b], 343.2), &energy[first + 3 * b]) - cpu[b]) < 0.02);
}

// A row's profile rays (PointPathOptions::profile_pe) take their field from
// the CPU, at every sub-band, as complex per-ray factors: the GPU's levels
// are band_levels', the profile weighed in wholly and in part, on a rigid
// floor, one reflecting less, grass, and a plate under a flush microphone.
TEST_CASE("the GPU takes a row's profile rays as the CPU does", "[gpu][profile_pe]")
{
    gpu::Device& device = device_or_skip();
    const Buildings street {{block(0.0, -400.0, 10.0, 400.0, 10.0), block(20.0, -400.0, 30.0, 400.0, 10.0)}};
    const std::vector<double> bands {125, 250, 500, 1000};
    const double c = 343.2;
    const std::vector<double> reflection(street.facades().size() * bands.size(), 0.9);
    const std::vector<double> middles = sub_band_middles(bands, c);
    FloorModel grass;
    grass.admittance = ground_impedance(GroundModel {}, middles, c);
    FloorModel plate = grass;
    plate.plate = [&](double elevation, double azimuth, std::complex<double>* direct, std::complex<double>* image) {
        for (std::size_t j = 0; j < middles.size(); j++)
        {
            direct[j] = std::polar(1.0 + 0.3 * std::sin(elevation), 0.2 * azimuth + 0.01 * j);
            image[j] = std::polar(0.9, -0.1 * j);
        }
    };
    std::size_t whole = 0, part = 0;
    // Grazing, and the source 33 degrees over the first roof.
    for (const Vec3& source : {Vec3(-60.0, -10.0, 12.0), Vec3(-50.0, 0.0, 42.0)})
        for (const FloorModel& floor : {FloorModel(1.0), FloorModel(0.7), grass, plate})
        {
            const Vec3 microphone(34.0, 10.0, floor.plate ? 0.0 : 1.2);
            const PointSourceField field(street, source, microphone);
            for (const PointRay& ray : field.rays())
                if (ray.profile)
                    (ray.share < 1.0 ? part : whole)++;
            std::vector<double> cpu(bands.size());
            field.band_levels(reflection, floor, bands, c, cpu.data());
            gpu::PointRaySet set;
            gpu::pack_point_rays(field, reflection, floor, bands, c, set);
            REQUIRE(set.complex_reflections);
            gpu::BandEnergies energies(device);
            const std::size_t first = energies.add(set);
            energies.submit();
            const std::vector<double>& energy = energies.wait();
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                INFO("band " << bands[b] << ", floor " << floor.reflection << (floor.impedance() ? ", grass" : "")
                             << (floor.plate ? ", plate" : ""));
                REQUIRE(std::abs(band_level(bands[b], sub_bands(bands[b], c), &energy[first + 3 * b]) - cpu[b]) < 0.02);
            }
        }
    REQUIRE(whole == 16);
    REQUIRE(part == 16);
}

// A tilted step's jump split (TiltedStepSplit) rides on the GPU as the foot
// fields' rays, each with its weight and fade per sub-band among the complex
// per-ray factors and its carry in its length: the GPU's levels are
// band_levels', with microphones 0.02 m either side of the line through the
// step's far edges and 0.5 m off it, on a rigid floor, one reflecting less
// and grass, with and without the air.
TEST_CASE("the GPU takes a tilted step's split as the CPU does", "[gpu]")
{
    gpu::Device& device = device_or_skip();
    auto long_block = [](double y0, double y1, double top) {
        Footprint b;
        b.x = {-4000.0, 4000.0, 4000.0, -4000.0};
        b.y = {y0, y0, y1, y1};
        b.top = top;
        return b;
    };
    const Buildings step({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    const std::vector<double> bands {125, 250, 500, 1000};
    const double c = 343.2;
    const std::vector<double> reflection(step.facades().size() * bands.size(), 0.9);
    FloorModel grass;
    grass.admittance = ground_impedance(GroundModel {}, sub_band_middles(bands, c), c);
    Atmosphere air;
    air.relative_humidity = 50.0;
    const std::vector<double> absorbing = sub_band_absorption(air, bands, c);
    // 3.3 degrees above the line, 60 m from the taller block's edge.
    const double t = std::atan(0.3), a = t + 3.3 * pi / 180.0;
    const Vec3 source(5.0, 10.0 - 60.0 * std::cos(a), 10.0 + 60.0 * std::sin(a));
    std::size_t split = 0;
    for (const double off : {-0.5, -0.02, 0.02, 0.5})
        for (const FloorModel& floor : {FloorModel(1.0), FloorModel(0.7), grass})
            for (int absorbed = 0; absorbed < 2; absorbed++)
            {
                const Vec3 microphone(0.0, 30.0 + off * std::sin(t), 10.0 - 0.3 * 20.0 + off * std::cos(t));
                const PointSourceField field(step, source, microphone);
                const std::vector<double> absorption = absorbed ? absorbing : std::vector<double> {};
                split += std::any_of(field.rays().begin(), field.rays().end(),
                        [](const PointRay& r) { return r.jump_weight != 0.0; });
                std::vector<double> cpu(bands.size());
                field.band_levels(reflection, floor, bands, c, cpu.data(), absorption);
                gpu::PointRaySet set;
                gpu::pack_point_rays(field, reflection, floor, bands, c, set, absorption);
                REQUIRE(set.complex_reflections);
                gpu::BandEnergies energies(device);
                const std::size_t first = energies.add(set);
                energies.submit();
                const std::vector<double>& energy = energies.wait();
                for (std::size_t b = 0; b < bands.size(); b++)
                {
                    INFO("band " << bands[b] << ", " << off << " m off the line, floor " << floor.reflection
                                 << (floor.impedance() ? ", grass" : "") << (absorbed ? ", air" : ""));
                    REQUIRE(std::abs(band_level(bands[b], sub_bands(bands[b], c), &energy[first + 3 * b]) - cpu[b])
                            < 0.02);
                }
            }
    REQUIRE(split == 24);
}
