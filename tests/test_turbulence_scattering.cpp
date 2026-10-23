#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "buildings.hpp"
#include "turbulence_scattering.hpp"

using namespace niceops; // NOLINT

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

} // namespace

// Forssén et al. (2014), eq. 4 without a canyon, velocity turbulence: 10 dB a
// decade of C_v^2 and 10/3 a decade of frequency; at h = 10 m, d_S = d_R = 40 m
// and 1 kHz, -52.8 - 17.1 lg(0.0625 + 0.0012).  Its Figure 2 (the scattering
// integral at 1 kHz, C_v^2 = 1, d_S = 40 m, no ground) runs from about -63 dB
// at h = 40 m, d_R = 10 m to about -28 dB at h = 5 m, d_R = 80 m: the fit, which
// doubles the energy at both ends, less 6 dB, is within 2 dB of both.  gamma_T
// is velocity turbulence of (3/22) c0^2 gamma_T, which gamma_T's own definition
// gives back.
TEST_CASE("the scattered level is Forssén et al.'s fit", "[turbulence_scattering]")
{
    REQUIRE(scattered_level(10.0, 40.0, 40.0, 1000.0, 1.0)
            == Catch::Approx(-52.8 - 17.1 * std::log10(0.0625 + 0.0012)).margin(1e-12));
    const double base = scattered_level(12.0, 30.0, 90.0, 1000.0, 0.1);
    REQUIRE(scattered_level(12.0, 30.0, 90.0, 10000.0, 0.1) - base == Catch::Approx(10.0 / 3.0).margin(1e-12));
    REQUIRE(scattered_level(12.0, 30.0, 90.0, 1000.0, 1.0) - base == Catch::Approx(10.0).margin(1e-12));
    REQUIRE(std::abs(scattered_level(40.0, 40.0, 10.0, 1000.0, 1.0) - 6.0 - -63.0) < 2.0);
    REQUIRE(std::abs(scattered_level(5.0, 40.0, 80.0, 1000.0, 1.0) - 6.0 - -28.0) < 2.0);
    REQUIRE(std::isinf(scattered_level(0.0, 40.0, 40.0, 1000.0, 1.0)));
    REQUIRE(std::isinf(scattered_level(10.0, 40.0, 40.0, 1000.0, 0.0)));
    const double c0 = 343.2, gamma = 3e-6;
    const double cv2 = velocity_equivalent(gamma, c0);
    REQUIRE(22.0 / 3.0 * cv2 / (c0 * c0) == Catch::Approx(gamma).epsilon(1e-12));
}

// A box between the source and the receiver, the line between them rising:
// of the two roof tops above the line where its plan crosses the box's sides,
// the near one (x = -5) stands farther above it, 8.5 m in height and 8.458 m
// across the line.  Below the line, or with no building in the way, nothing
// shields.
TEST_CASE("the shadow's screen is the roof top farthest above the line", "[turbulence_scattering]")
{
    const Buildings city({box(-5.0, -10.0, 5.0, 10.0, 12.0)});
    const Vec3 source(-20.0, 0.0, 2.0), receiver(20.0, 0.0, 6.0);
    const ShadowScreen screen = shadow_screen(city, source, receiver);
    REQUIRE(screen.shielded);
    const double length = std::hypot(40.0, 4.0), along = (15.0 * 40.0 + 10.0 * 4.0) / length;
    REQUIRE(screen.d_s == Catch::Approx(along).epsilon(1e-12));
    REQUIRE(screen.d_r == Catch::Approx(length - along).epsilon(1e-12));
    REQUIRE(screen.h == Catch::Approx(std::sqrt(15.0 * 15.0 + 10.0 * 10.0 - along * along)).epsilon(1e-12));
    REQUIRE(screen.h == Catch::Approx(8.5 * 40.0 / length).epsilon(1e-12));
    REQUIRE(!shadow_screen(city, Vec3(-20.0, 0.0, 20.0), Vec3(20.0, 0.0, 15.0)).shielded);
    REQUIRE(!shadow_screen(city, Vec3(-20.0, 30.0, 2.0), Vec3(20.0, 30.0, 2.0)).shielded);
    // A source over the roof, the receiver beside the box: the line leaves
    // the roof above it.
    REQUIRE(!shadow_screen(city, Vec3(0.0, 0.0, 20.0), Vec3(20.0, 0.0, 1.5)).shielded);
}
