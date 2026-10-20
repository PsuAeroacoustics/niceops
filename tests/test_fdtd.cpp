#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <complex>
#include <vector>

#include "buildings.hpp"
#include "edge_source.hpp"
#include "fdtd.hpp"
#include "gpu/gpu.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"
#include "wedge_series.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

// A third octave's level of a field known at every frequency: 10 lg of |p|^2
// averaged over the band, as FdtdRecording::band_level takes it.
template <typename Field>
double band_level_of(double band, Field field, int points = 48)
{
    const double half = std::pow(2.0, 1.0 / 6.0), f1 = band / half, f2 = band * half;
    double sum = 0.0;
    for (int m = 0; m < points; m++)
        sum += std::norm(field(f1 + (m + 0.5) * (f2 - f1) / points));
    return 10.0 * std::log10(sum / points);
}

FdtdSetup small_setup()
{
    FdtdSetup s;
    s.cell = 0.01;
    s.pml = 20;
    s.pulse_frequency = 700.0;
    s.threads = default_thread_count();
    return s;
}

} // namespace

TEST_CASE("over a rigid floor the FDTD's field is the two rays'", "[fdtd]")
{
    // 25 cells a wavelength at the top band's upper edge.  Receivers near the
    // source, far from it, on the floor, and across the symmetry plane.
    const int order = GENERATE(2, 4);
    INFO("order " << order);
    FdtdSetup s = small_setup();
    s.order = order;
    s.lo = Vec3(-0.25, -0.25, 0.0);
    s.hi = Vec3(0.6, 0.65, 0.6);
    s.duration = 0.007;
    const Vec3 source(0.0, 0.0, 0.25);
    const std::vector<Vec3> receivers {{0.4, 0.0, 0.1}, {0.3, 0.2, 0.4}, {0.0, 0.5, 0.003}, {-0.1, 0.1, 0.5}};
    const FdtdRecording run = solve_fdtd(s, {}, source, receivers);
    const double c = s.sound_speed;
    const Vec3 image(source.x, source.y, -source.z);
    for (std::size_t r = 0; r < receivers.size(); r++)
    {
        REQUIRE(std::abs(run.truncation_db(r, 500.0)) < 1e-3);
        for (double band : {500.0, 630.0, 800.0, 1000.0, 1250.0})
        {
            const double exact = band_level_of(band, [&](double f) {
                const double k = 2.0 * pi * f / c, a = (receivers[r] - source).norm(), b = (receivers[r] - image).norm();
                return std::polar(1.0 / a, k * a) + std::polar(1.0 / b, k * b);
            });
            const double error = run.band_level(r, band) - exact;
            INFO("receiver " << r + 1 << ", " << band << " Hz: " << error << " dB");
            // Fourth order 5e-4 dB at worst; taking Q's derivative as
            // -i omega would leave -0.0023 dB at 1250 Hz.
            CHECK(std::abs(error) < (order == 4 ? 0.0015 : 0.1));
        }
    }

    // The air's absorption, taken afterwards by each arrival's delay: each
    // ray loses exp(-a L) along its own length.  Fourth order 0.0018 dB at
    // worst; weighting the pressure alone, not Q, leaves 0.013.
    const double a = 0.4;
    const double attenuation[3] = {a, a, a};
    for (std::size_t r = 0; r < receivers.size(); r++)
        for (double band : {500.0, 1000.0})
        {
            const double exact = band_level_of(band, [&](double f) {
                const double k = 2.0 * pi * f / c, d1 = (receivers[r] - source).norm(), d2 = (receivers[r] - image).norm();
                return std::polar(std::exp(-a * d1) / d1, k * d1) + std::polar(std::exp(-a * d2) / d2, k * d2);
            });
            const double error = run.band_level(r, band, 0, attenuation) - exact;
            INFO("absorbed: receiver " << r + 1 << ", " << band << " Hz: " << error << " dB");
            CHECK(std::abs(error) < (order == 4 ? 0.003 : 0.1));
        }

    // The same with the plane y = 0 taken as one of symmetry, the receivers
    // off it on either side.
    s.mirror = true;
    s.lo.y = 0.0;
    std::vector<Vec3> mirrored = receivers;
    for (Vec3& p : mirrored)
        p.y = std::abs(p.y);
    const FdtdRecording half = solve_fdtd(s, {}, source, mirrored);
    for (std::size_t r = 0; r < receivers.size(); r++)
        for (double band : {500.0, 800.0, 1250.0})
            CHECK(std::abs(half.band_level(r, band) - run.band_level(r, band)) < 0.02);
    // A source off the plane of symmetry is refused.
    REQUIRE_THROWS_AS(solve_fdtd(s, {}, Vec3(0.0, 0.1, 0.25), mirrored), std::invalid_argument);
    // So is a receiver a fraction of a cell inside a solid.
    FdtdSolid box;
    box.x = {0.2, 0.3, 0.3, 0.2};
    box.y = {-0.3, -0.3, 0.3, 0.3};
    box.top = 0.2;
    REQUIRE_THROWS_AS(solve_fdtd(s, {box}, source, {Vec3(0.203, 0.0, 0.1)}), std::invalid_argument);
}

TEST_CASE("by a rigid wedge the FDTD's field is the edge-source integral's", "[fdtd]")
{
    // A right-angled wedge, its edge along y at x = 0, z = 1 m, solid below the
    // roof and beyond the wall; no floor.  The integral over the whole line is
    // its exact field (edge_source.hpp).  The edge runs through the PML at
    // y = 0.6 m, from which the waves still to come are 30 dB or more down at
    // these distances and frequencies.  Receivers deep in the shadow above the
    // roof, off the plane of the source, lit above the shadow boundary, and on
    // the source's side where the wall reflects.
    const double top = 1.0;
    Footprint block;
    block.x = {0.0, 50.0, 50.0, 0.0};
    block.y = {-50.0, -50.0, 50.0, 50.0};
    block.top = top;
    const Buildings city({block});
    const Edge* edge = nullptr;
    for (const Edge& e : city.edges())
        if (std::abs(e.direction.y) > 0.99 && std::abs(e.start.x) < 1e-9 && std::abs(e.start.z - top) < 1e-9)
            edge = &e;
    REQUIRE(edge != nullptr);

    const int order = GENERATE(2, 4);
    INFO("order " << order);
    FdtdSetup s = small_setup();
    s.order = order;
    s.floor = false;
    s.mirror = true;
    s.lo = Vec3(-0.25, 0.0, 0.75);
    s.hi = Vec3(0.3, 0.6, 1.3);
    s.duration = 0.008;
    const Vec3 source(-0.15, 0.0, 0.9);
    const std::vector<Vec3> receivers {{0.12, 0.0, 1.03}, {0.2, 0.0, 1.06}, {0.12, 0.1, 1.03}, {0.1, 0.0, 1.2},
            {-0.2, 0.0, 1.1}};
    const FdtdRecording run = solve_fdtd(s, fdtd_solids({block}), source, receivers);
    const double c = s.sound_speed;
    double square = 0.0;
    int count = 0;
    for (std::size_t r = 0; r < receivers.size(); r++)
    {
        const EdgeFrameEnds ends = edge_frame_ends(*edge, source, receivers[r]);
        for (double band : {500.0, 630.0, 800.0, 1000.0})
        {
            const double exact = band_level_of(band, [&](double f) {
                const double k = 2.0 * pi * f / c;
                return test::wedge_geometrical(edge->n, ends, k) + edge_source_field(*edge, source, receivers[r], k);
            });
            const double error = run.band_level(r, band) - exact;
            INFO("receiver " << r + 1 << ", " << band << " Hz: " << error << " dB");
            CHECK(std::abs(error) < (order == 4 ? 0.06 : 0.25));
            square += error * error;
            count++;
        }
    }
    CHECK(std::sqrt(square / count) < (order == 4 ? 0.03 : 0.12));
}

TEST_CASE("an absorbing FDTD wall reflects as its impedance says", "[fdtd]")
{
    // A wall at x = 0.3 m across the whole domain, absorbing 0.36 at normal
    // incidence: zeta = 9, R = (zeta cos - 1) / (zeta cos + 1).  No floor; the
    // source and the receivers 0.3 m or less from it, on its normal and off
    // it, where the plane wave's R holds for a wall this hard.
    const int order = GENERATE(2, 4);
    INFO("order " << order);
    FdtdSetup s = small_setup();
    s.order = order;
    s.floor = false;
    s.mirror = true;
    s.wall_absorption = 0.36;
    s.lo = Vec3(-0.35, 0.0, -0.3);
    s.hi = Vec3(0.3, 0.3, 0.3);
    s.duration = 0.006;
    FdtdSolid wall;
    wall.x = {0.3, 5.0, 5.0, 0.3};
    wall.y = {-5.0, -5.0, 5.0, 5.0};
    wall.bottom = -5.0;
    wall.top = 5.0;
    const Vec3 source(0.0, 0.0, 0.0);
    const std::vector<Vec3> receivers {{-0.2, 0.0, 0.0}, {0.1, 0.0, 0.15}, {0.0, 0.2, -0.1}};
    const FdtdRecording run = solve_fdtd(s, {wall}, source, receivers);
    const Vec3 image(0.6, 0.0, 0.0);
    const double zeta = 9.0;
    for (std::size_t r = 0; r < receivers.size(); r++)
        for (double band : {630.0, 800.0, 1000.0})
        {
            const double exact = band_level_of(band, [&](double f) {
                const double k = 2.0 * pi * f / s.sound_speed, d1 = (receivers[r] - source).norm();
                const double d2 = (receivers[r] - image).norm(), cosine = (image.x - receivers[r].x) / d2;
                const double R = (zeta * cosine - 1.0) / (zeta * cosine + 1.0);
                return std::polar(1.0 / d1, k * d1) + R * std::polar(1.0 / d2, k * d2);
            });
            const double error = run.band_level(r, band) - exact;
            INFO("receiver " << r + 1 << ", " << band << " Hz: " << error << " dB");
            CHECK(std::abs(error) < 0.15);
        }
    REQUIRE_THROWS_AS(
            [&] {
                FdtdSetup bad = s;
                bad.wall_absorption = 1.0;
                return solve_fdtd(bad, {wall}, source, receivers);
            }(),
            std::invalid_argument);
}

TEST_CASE("in a concave corner the FDTD's field is the source and its three images", "[fdtd]")
{
    // A floor at z = 0.3 m and a wall at y = 0.5 m, both through the whole
    // domain and on the grid's faces: the exact field is the source and its
    // images in each plane and in both.  No floor of the domain's own.
    const int order = GENERATE(2, 4);
    INFO("order " << order);
    FdtdSetup s = small_setup();
    s.order = order;
    s.floor = false;
    s.lo = Vec3(-0.5, -0.5, 0.3);
    s.hi = Vec3(0.5, 0.5, 0.9);
    s.duration = 0.009;
    FdtdSolid floor, wall;
    floor.x = {-50.0, 50.0, 50.0, -50.0};
    floor.y = {-50.0, -50.0, 50.0, 50.0};
    floor.bottom = -50.0;
    floor.top = 0.3;
    wall.x = floor.x;
    wall.y = {0.5, 0.5, 50.0, 50.0};
    wall.bottom = -50.0;
    wall.top = 50.0;
    const Vec3 source(-0.25, 0.1, 0.45);
    const std::vector<Vec3> receivers {{0.25, 0.1, 0.35}, {0.2, 0.35, 0.55}, {0.0, 0.45, 0.32}};
    const FdtdRecording run = solve_fdtd(s, {floor, wall}, source, receivers);
    const Vec3 images[3] = {{source.x, source.y, 0.6 - source.z}, {source.x, 1.0 - source.y, source.z},
            {source.x, 1.0 - source.y, 0.6 - source.z}};
    for (std::size_t r = 0; r < receivers.size(); r++)
        for (double band : {500.0, 630.0, 800.0})
        {
            const double exact = band_level_of(band, [&](double f) {
                const double k = 2.0 * pi * f / s.sound_speed;
                std::complex<double> p = std::polar(1.0 / (receivers[r] - source).norm(), k * (receivers[r] - source).norm());
                for (const Vec3& image : images)
                    p += std::polar(1.0 / (receivers[r] - image).norm(), k * (receivers[r] - image).norm());
                return p;
            });
            const double error = run.band_level(r, band) - exact;
            INFO("receiver " << r + 1 << ", " << band << " Hz: " << error << " dB");
            CHECK(std::abs(error) < (order == 4 ? 0.005 : 0.05));
        }
}

TEST_CASE("fourth order with the warp holds the phase over twenty wavelengths at eight cells a wavelength", "[fdtd]")
{
    // Two rays over a floor out to 2.3 m, 1-2.5 kHz bands, 8 cells a
    // wavelength at the top band's upper edge.  Over that many wavelengths
    // the second order's anisotropic phase error moves the two rays' nulls
    // (several dB at 8 cells, more than 1 dB at 12); the fourth order's is
    // nearly the same in every direction, and the warp takes it out.
    const double c = 343.2, top = 2500.0 * std::pow(2.0, 1.0 / 6.0);
    const Vec3 source(0.0, 0.0, 0.3), image(0.0, 0.0, -0.3);
    const std::vector<Vec3> receivers {{0.5, 0.0, 0.05}, {1.5, 0.2, 0.3}, {2.3, 0.0, 0.003}, {2.2, 0.3, 0.5}};
    auto run = [&](int order, double points) {
        FdtdSetup s;
        s.order = order;
        s.cell = c / (top * points);
        s.mirror = true;
        s.lo = Vec3(-0.2, 0.0, 0.0);
        s.hi = Vec3(2.5, 0.5, 0.7);
        s.duration = 0.013;
        s.pulse_frequency = top / 2.25;
        s.threads = default_thread_count();
        const FdtdRecording r = solve_fdtd(s, {}, source, receivers);
        double worst = 0.0, square = 0.0;
        int n = 0;
        for (std::size_t m = 0; m < receivers.size(); m++)
            for (double band : {1000.0, 1250.0, 1600.0, 2000.0, 2500.0})
            {
                const double exact = band_level_of(band, [&](double f) {
                    const double k = 2.0 * pi * f / c, a = (receivers[m] - source).norm(), b = (receivers[m] - image).norm();
                    return std::polar(1.0 / a, k * a) + std::polar(1.0 / b, k * b);
                });
                const double error = r.band_level(m, band) - exact;
                worst = std::max(worst, std::abs(error));
                square += error * error;
                n++;
            }
        return std::make_pair(worst, std::sqrt(square / n));
    };
    const auto fourth = run(4, 8.0), second = run(2, 8.0);
    INFO("fourth order: worst " << fourth.first << " dB, rms " << fourth.second << "; second: " << second.first << ", "
                                << second.second);
    CHECK(fourth.first < 0.5);
    CHECK(fourth.second < 0.15);
    CHECK(second.first > 2.0);
}

TEST_CASE("the FDTD's fastest waves are the dispersion relation's", "[fdtd]")
{
    // Second order never outruns c; fourth order's leapfrog does, near the
    // Nyquist wavenumber on a diagonal, which the cull's cones must allow
    // for.  Maxima of |grad omega| / c over a fine lattice, computed apart.
    const double second = 0.99 / std::sqrt(3.0), fourth = second * 6.0 / 7.0;
    CHECK(std::abs(fdtd_fastest_group_speed(2, second) - 1.0) < 1e-12);
    CHECK(std::abs(fdtd_fastest_group_speed(4, fourth) - 1.22851) < 1e-3);
    CHECK(fdtd_fastest_group_speed(4, fourth) >= 1.228);
    CHECK(std::abs(fdtd_fastest_group_speed(4, 0.3) - fdtd_fastest_group_speed(4, fourth)) > 0.05);
}

TEST_CASE("culling the march to its light cones leaves the recordings a full march's", "[fdtd]")
{
    // Along a 1.8 m strip, past a block and beside a wall: the source's
    // cone reaches the far end late, and the receivers' cones shrink to
    // nothing at the end.  Four cells of margin take a culled march to the
    // float floor of a full one, a few 1e-7 of the peak; none leaves 1e-5
    // (fourth order) to 2e-4 (second), and the margin is 16.
    Footprint block;
    block.x = {0.6, 0.8, 0.8, 0.6};
    block.y = {-0.3, -0.3, 0.3, 0.3};
    block.top = 0.2;
    Footprint wall;
    wall.x = {0.3, 1.4, 1.4, 0.3};
    wall.y = {0.25, 0.25, 0.28, 0.28};
    wall.top = 0.35;
    const int order = GENERATE(2, 4);
    INFO("order " << order);
    FdtdSetup s = small_setup();
    s.order = order;
    s.mirror = true;
    s.lo = Vec3(-0.2, 0.0, 0.0);
    s.hi = Vec3(1.6, 0.3, 0.4);
    s.pulse_frequency = 343.2 / (s.cell * (order == 4 ? 8.0 : 12.0)) / 2.25;
    s.duration = 2.0 * (s.hi - s.lo).norm() / 343.2;
    const Vec3 source(0.0, 0.0, 0.15);
    const std::vector<Vec3> receivers {{0.4, 0.05, 0.1}, {1.0, 0.0, 0.02}, {1.5, 0.15, 0.3}, {0.5, 0.2, 0.38}};
    const std::vector<FdtdSolid> solids = fdtd_solids({block, wall});
    s.cull = false;
    const FdtdRecording full = solve_fdtd(s, solids, source, receivers);
    s.cull = true;
    const FdtdRecording culled = solve_fdtd(s, solids, source, receivers);
    CHECK(full.marched == 1.0);
    CHECK(culled.marched < 0.85);
    REQUIRE(culled.pressure.size() == full.pressure.size());
    for (std::size_t r = 0; r < receivers.size(); r++)
    {
        double peak = 0.0, difference = 0.0;
        for (std::size_t n = 0; n < full.steps; n++)
        {
            const double a = full.pressure[n * receivers.size() + r], b = culled.pressure[n * receivers.size() + r];
            peak = std::max(peak, std::abs(a));
            difference = std::max(difference, std::abs(a - b));
        }
        INFO("receiver " << r + 1 << ": " << difference / peak << " of its peak");
        CHECK(difference < 4e-6 * peak);
    }
}

TEST_CASE("the FDTD marches the same on the GPU", "[fdtd][gpu]")
{
    if (!gpu::Device::get())
        SKIP("no GPU backend or device (backend " << gpu::Device::backend() << ")");
    // A block on the floor, half of it across the symmetry plane, its walls
    // absorbing, so that every kind of face and every PML slab is in play.
    // Without the plane a tricubic source has all 64 of its cells, which the
    // injection must reach: it once took 32, half the source.
    Footprint block;
    block.x = {0.1, 0.3, 0.3, 0.1};
    block.y = {-0.2, -0.2, 0.2, 0.2};
    block.top = 0.2;
    const int order = GENERATE(2, 4);
    const bool mirror = GENERATE(true, false);
    INFO("order " << order << (mirror ? ", half the space" : ", the whole space"));
    FdtdSetup s = small_setup();
    s.order = order;
    s.mirror = mirror;
    s.lo = Vec3(-0.2, mirror ? 0.0 : -0.35, 0.0);
    s.hi = Vec3(0.5, 0.35, 0.4);
    s.duration = 0.004;
    const Vec3 source(0.0, 0.0, 0.1);
    const std::vector<Vec3> receivers {{0.45, 0.0, 0.02}, {0.2, 0.3, 0.1}, {0.2, 0.1, 0.3}};
    s.wall_absorption = 0.2;
    const FdtdRecording cpu = solve_fdtd(s, fdtd_solids({block}), source, receivers);
    s.gpu = true;
    const FdtdRecording gpu_run = solve_fdtd(s, fdtd_solids({block}), source, receivers);
    REQUIRE(gpu_run.gpu);
    REQUIRE_FALSE(cpu.gpu);
    REQUIRE(gpu_run.pressure.size() == cpu.pressure.size());
    double largest = 0.0, difference = 0.0;
    for (std::size_t i = 0; i < cpu.pressure.size(); i++)
    {
        largest = std::max(largest, std::abs(static_cast<double>(cpu.pressure[i])));
        difference = std::max(difference, std::abs(static_cast<double>(cpu.pressure[i]) - gpu_run.pressure[i]));
    }
    // The same float arithmetic, but for contractions into fused multiply-adds.
    CHECK(difference < 1e-4 * largest);
    for (std::size_t r = 0; r < receivers.size(); r++)
        CHECK(std::abs(gpu_run.band_level(r, 800.0) - cpu.band_level(r, 800.0)) < 1e-3);
}
