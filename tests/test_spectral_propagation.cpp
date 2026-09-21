// Propagating a spectrum to the ground, band by band.
//
// None of the shipped databases carries spectra, so these build one the way
// tests/test_sphere_database.cpp does: attach a synthetic spectrum to an
// interpolated sphere, write it out, and read it back. That keeps the whole
// read/blend/propagate path under test without a database from panam.
//
// The checks are arranged so the answer is known in closed form. Absorption is
// an argument rather than something the Database computes, so a test can pass
// zeros and isolate the spreading, then pass a real coefficient and isolate the
// absorption.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <limits>
#include <vector>

#include "atmosphere.hpp"
#include "sphere_database.hpp"
#include "weighting.hpp"

using namespace niceops;

namespace {

// The band centres and the constant per-band level the fixture is built with.
const std::vector<double> kBands = {125.0, 500.0, 2000.0, 8000.0};
constexpr double kFlatLevel = 90.0;

// Write a single-condition database whose every cell carries the same flat
// spectrum, so anything the propagation does shows up undiluted by directivity.
std::string write_flat_spectral_database(const std::string& name)
{
    Database plain("databases/S-76D_M3.nod");
    auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    sphere.frequencies = kBands;
    sphere.amplitudes.assign(sphere.cells.size() * kBands.size(), kFlatLevel);

    const auto path = (std::filesystem::temp_directory_path() / name).string();
    write_sphere(sphere, path);
    return path;
}

void check_close(double got, double expected, double tolerance)
{
    REQUIRE(std::fabs(got - expected) <= tolerance);
}

} // namespace

TEST_CASE("eval_spectrum applies spreading and per-band absorption", "[spectral]")
{
    const auto path = write_flat_spectral_database("niceops_flat_spectrum.nod");
    Database db(path);
    REQUIRE(db.has_spectra());
    REQUIRE(db.band_centres() == kBands);

    const auto blend = db.blend_at(0.15, -4.0, 0.006);
    const double radius = db.sphere_radius();
    const double range = 10.0 * radius;
    const double expected_spreading = 20.0 * std::log10(range / radius);

    std::vector<double> levels(kBands.size());

    SECTION("zero absorption leaves only spherical spreading")
    {
        const std::vector<double> none(kBands.size(), 0.0);
        db.eval_spectrum(range, 0.0, 90.0, blend, none.data(), levels.data());

        // Flat source spectrum, no absorption: every band falls by the same
        // spreading loss and they stay equal.
        for (const double level : levels)
            check_close(level, kFlatLevel - expected_spreading, 1e-9);
    }

    SECTION("absorption is applied per band over the path beyond the sphere")
    {
        Atmosphere air; // panam's reference day
        const auto alpha = absorption_coefficients(air, kBands);
        db.eval_spectrum(range, 0.0, 90.0, blend, alpha.data(), levels.data());

        // The path length is measured from the sphere surface, not the source:
        // the level on the sphere already contains whatever happened inside it.
        const double path = range - radius;
        for (std::size_t band = 0; band < kBands.size(); band++)
        {
            check_close(levels[band], kFlatLevel - expected_spreading - alpha[band] * path, 1e-9);
        }

        // High frequencies must come out further down than low ones. At this
        // range that ordering is strict.
        for (std::size_t band = 1; band < kBands.size(); band++)
            REQUIRE(levels[band] < levels[band - 1]);
    }

    SECTION("absorption grows with range, spreading grows logarithmically")
    {
        Atmosphere air;
        const auto alpha = absorption_coefficients(air, kBands);
        std::vector<double> near(kBands.size());
        std::vector<double> far(kBands.size());
        db.eval_spectrum(2.0 * radius, 0.0, 90.0, blend, alpha.data(), near.data());
        db.eval_spectrum(200.0 * radius, 0.0, 90.0, blend, alpha.data(), far.data());

        // Doubling distance costs 6 dB of spreading regardless of band, so the
        // spread between the top and bottom band is absorption alone, and it
        // has to widen with range.
        const double near_spread = near.front() - near.back();
        const double far_spread = far.front() - far.back();
        REQUIRE(far_spread > near_spread);
    }

    std::filesystem::remove(path);
}

TEST_CASE("spectral_level_at_observer weights and band-limits", "[spectral]")
{
    const auto path = write_flat_spectral_database("niceops_flat_observer.nod");
    Database db(path);

    const auto state = db.source_state(Vec3(0.0, 0.0, 300.0), Vec3(60.0, 0.0, 0.0),
            Vec3(0.0, 0.0, 1.0), 0.15, -4.0, 0.006);
    const Vec3 observer(0.0, 0.0, 0.0);
    const std::vector<double> no_absorption(kBands.size(), 0.0);

    SECTION("a single retained band reduces to that band's level")
    {
        // Keep only the 500 Hz band, so the A-weighted total is that band's
        // ground level plus the A-weight there -- no summing to obscure it.
        BandLimits only_500;
        only_500.low = 400.0;
        only_500.high = 600.0;
        const auto gains = a_weighting_gains(kBands, only_500);
        REQUIRE(bands_in_range(kBands, only_500) == 1);

        const double got = db.spectral_level_at_observer(observer, state,
                no_absorption.data(), gains);

        std::vector<double> levels(kBands.size());
        const double range = (observer - state.position).norm();
        const auto angles = lookup_angles(state.frame,
                (observer - state.position) / range);
        db.eval_spectrum(range, angles[0], angles[1], state.blend, no_absorption.data(),
                levels.data());

        check_close(got, levels[1] + a_weighting(500.0), 1e-9);
    }

    SECTION("band limits only ever remove energy")
    {
        const auto wide = a_weighting_gains(kBands, BandLimits{});
        BandLimits narrow;
        narrow.low = 400.0;
        narrow.high = 2500.0;
        const auto narrowed = a_weighting_gains(kBands, narrow);

        const double all = db.spectral_level_at_observer(observer, state,
                no_absorption.data(), wide);
        const double some = db.spectral_level_at_observer(observer, state,
                no_absorption.data(), narrowed);
        REQUIRE(some < all);
    }

    SECTION("excluding every band is silence, not zero")
    {
        BandLimits nothing;
        nothing.low = 1.0e5;
        nothing.high = 2.0e5;
        const auto gains = a_weighting_gains(kBands, nothing);
        REQUIRE(db.spectral_level_at_observer(observer, state, no_absorption.data(), gains)
                == -std::numeric_limits<double>::infinity());
    }

    SECTION("absorption lowers the weighted level")
    {
        Atmosphere air;
        const auto alpha = absorption_coefficients(air, kBands);
        const auto gains = a_weighting_gains(kBands, BandLimits{});
        const double dry = db.spectral_level_at_observer(observer, state, alpha.data(), gains);
        const double none = db.spectral_level_at_observer(observer, state,
                no_absorption.data(), gains);
        REQUIRE(dry < none);
    }

    std::filesystem::remove(path);
}
