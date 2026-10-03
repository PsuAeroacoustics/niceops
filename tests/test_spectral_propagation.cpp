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
#include "ground_effect.hpp"
#include "propagation.hpp"
#include "sphere_database.hpp"
#include "weighting.hpp"

using namespace niceops;

namespace {

// The band centers and the constant per-band level the fixture is built with.
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
    REQUIRE(db.band_centers() == kBands);

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

TEST_CASE("band_levels are the bands the weighted level sums", "[spectral]")
{
    // --band_history writes these per band.  A-weighted and summed in energy
    // they must give back the level returned beside them, ground and extra
    // loss included, and asking for them must not change that level.
    const auto path = write_flat_spectral_database("niceops_flat_bands.nod");
    Database db(path);

    const auto state = db.source_state(Vec3(0.0, 0.0, 120.0), Vec3(60.0, 0.0, 0.0),
            Vec3(0.0, 0.0, 1.0), 0.15, -4.0, 0.006);
    const Vec3 observer(400.0, 150.0, 0.0);
    Atmosphere air;
    const auto alpha = absorption_coefficients(air, kBands);
    // 125 Hz is outside the limits, so first_band is 1, as run() would find it.
    BandLimits limits;
    limits.low = 300.0;
    const auto gains = a_weighting_gains(kBands, limits);
    const std::size_t first_band = 1;
    const std::size_t last_band = kBands.size();
    GroundModel ground;
    ground.receiver_height = 1.2;
    const auto impedance = ground_impedance(ground, kBands);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    auto check = [&](const std::vector<double>& bands, double level) {
        REQUIRE(std::isnan(bands[0]));
        double energy = 0.0;
        for (std::size_t band = first_band; band < last_band; band++)
            energy += gains[band] * std::pow(10.0, 0.1 * bands[band]);
        check_close(10.0 * std::log10(energy), level, 1e-9);
    };

    SECTION("without ground")
    {
        std::vector<double> bands(kBands.size(), nan);
        const double level = db.spectral_level_at_observer(observer, state, alpha.data(), gains,
                nullptr, nullptr, 343.2, first_band, last_band, nullptr, nullptr, 0.0, nullptr,
                bands.data());
        REQUIRE(level == db.spectral_level_at_observer(observer, state, alpha.data(), gains,
                nullptr, nullptr, 343.2, first_band, last_band));
        check(bands, level);
    }

    SECTION("with ground and an extra loss")
    {
        std::vector<double> bands(kBands.size(), nan);
        const std::vector<double> loss = {0.0, 3.0, 7.0, 11.0};
        const auto ray = straight_ray(state.position, observer);
        const double level = db.spectral_level_along(ray, state, alpha.data(), gains, &ground,
                &impedance, 343.2, first_band, last_band, nullptr, nullptr, 0.0, nullptr,
                loss.data(), nullptr, bands.data());
        REQUIRE(level == db.spectral_level_along(ray, state, alpha.data(), gains, &ground,
                &impedance, 343.2, first_band, last_band, nullptr, nullptr, 0.0, nullptr,
                loss.data()));
        check(bands, level);
    }

    std::filesystem::remove(path);
}

TEST_CASE("spectral_direction_bounds really bounds the spectral level", "[spectral][bounds]")
{
    // As "DirectionBounds really bounds the level" does for the broadband
    // path: the spectral peak search and --db_cutoff skip whatever this puts
    // below a threshold, so it must never be exceeded.  AS350's spectra have
    // empty (gated) bands in ~16% of cells, which is what the per-band half of
    // the bound exists for; EC130's exercise the per-spectrum half.
    for (const char* name : {"databases/AS350_fixed_load_factor.nod",
                 "databases/EC130_fixed_load_factor.nod"})
    {
        if (!std::filesystem::exists(name))
            continue;
        Database db(name);
        REQUIRE(db.has_spectra());
        const auto& bands = db.band_centers();
        const auto gains = a_weighting_gains(bands, BandLimits {});
        const std::vector<double> no_absorption(bands.size(), 0.0);

        std::vector<Database::SourceState> states;
        for (double mu : {0.0, 0.06, 0.12, 0.2})
            for (double gamma : {-9.0, -3.0, 0.0, 5.0})
                for (double ct : {0.004, 0.0075, 0.011}) // below and above the reference
                {
                    states.push_back(db.source_state(Vec3(0.0, 0.0, 60.0),
                            Vec3(40.0, 5.0, -3.0), Vec3(0.1, -0.05, 1.0), mu, gamma, ct));
                    states.push_back(db.source_state(Vec3(30.0, -50.0, 250.0),
                            Vec3(-20.0, 45.0, 2.0), Vec3(-0.08, 0.12, 1.0), mu, gamma, ct));
                }
        const auto bounds = db.spectral_direction_bounds(states, gains);
        REQUIRE(bounds.usable());
        // And with the ground effect: the gains raised by ground_gain_bounds.
        GroundModel ground;
        ground.receiver_height = 0.05;
        ground.turbulence = 8e-4;
        const auto impedance = ground_impedance(ground, bands);
        auto ground_gains = ground_gain_bounds(ground, impedance);
        for (std::size_t band = 0; band < bands.size(); band++)
            ground_gains[band] *= gains[band];
        const auto grounded = db.spectral_direction_bounds(states, ground_gains);
        REQUIRE(grounded.usable());

        const double radius = db.sphere_radius();
        std::size_t checked = 0;
        double tightest = std::numeric_limits<double>::infinity();
        for (std::size_t s = 0; s < states.size(); s++)
        {
            const auto packed = db.pack_spectra(states[s].blend.indices);
            for (double x = -1500.0; x <= 1500.0; x += 131.0)
            {
                for (double y = -1500.0; y <= 1500.0; y += 127.0)
                {
                    const Vec3 observer(x, y, 1.2);
                    const Vec3 delta = observer - states[s].position;
                    const double squared_range = delta.squaredNorm();
                    if (squared_range <= radius * radius)
                        continue;
                    const int bin = DirectionBounds::bin(states[s].frame.nose.dot(delta),
                            states[s].frame.normal.dot(delta),
                            states[s].frame.starboard.dot(delta), squared_range);
                    // No absorption: the tightest case the bound has to cover.
                    const double level = db.spectral_level_at_observer(observer, states[s],
                            no_absorption.data(), gains);
                    const double energy = std::pow(10.0, 0.1 * level) * squared_range;
                    REQUIRE(bounds.numerator(s, bin) >= energy);
                    const double with_ground = db.spectral_level_at_observer(observer, states[s],
                            no_absorption.data(), gains, &ground, &impedance, 343.2);
                    REQUIRE(grounded.numerator(s, bin)
                            >= std::pow(10.0, 0.1 * with_ground) * squared_range);
                    tightest = std::min(tightest, bounds.numerator(s, bin) / energy);
                    // The packed layout reads the same values: identical, not close.
                    REQUIRE(db.spectral_level_at_observer(observer, states[s],
                                    no_absorption.data(), gains, nullptr, nullptr, 343.2, 0,
                                    bands.size(), &packed)
                            == level);
                    checked++;
                }
            }
        }
        REQUIRE(checked > 10000);
        // Not vacuous: somewhere the bound is within 3 dB of the level.
        REQUIRE(tightest < 2.0);
    }
}

TEST_CASE("spectral_direction_bounds still bounds the level at the nose and tail axes",
        "[spectral][bounds]")
{
    // The spectral table files a cell under the bin of its center direction as the
    // broadband one does, and shared its gap: at theta = 0 or 180 that direction
    // has no azimuth, so a query a degree off the axis, in a bin the pole's value
    // was never written to, was bounded too low (0.14 dB on AS350).
    const char* name = "databases/AS350_fixed_load_factor.nod";
    if (!std::filesystem::exists(name))
        return;
    Database db(name);
    REQUIRE(db.has_spectra());
    const auto& bands = db.band_centers();
    const auto gains = a_weighting_gains(bands, BandLimits {});
    const std::vector<double> no_absorption(bands.size(), 0.0);

    std::vector<Database::SourceState> states;
    for (double mu = 0.0; mu <= 0.30001; mu += 0.02)
        for (double gamma = -12.0; gamma <= 12.0001; gamma += 3.0)
            for (double ct : {0.004, 0.0075, 0.011})
                states.push_back(db.source_state(Vec3(0.0, 0.0, 60.0), Vec3(40.0, 5.0, -3.0),
                        Vec3(0.1, -0.05, 1.0), mu, gamma, ct));
    const auto bounds = db.spectral_direction_bounds(states, gains);
    REQUIRE(bounds.usable());

    const double radius = db.sphere_radius();
    std::size_t checked = 0, exceeded = 0;
    for (std::size_t s = 0; s < states.size(); s++)
    {
        const auto& frame = states[s].frame;
        for (double theta : {0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0})
            for (double sign : {1.0, -1.0})
                for (double phi = -180.0; phi < 180.0; phi += 2.5)
                    for (double range : {150.0, 400.0, 1200.0})
                    {
                        const double t = theta * M_PI / 180.0, p = phi * M_PI / 180.0;
                        Vec3 direction = sign * std::cos(t) * frame.nose
                                + std::sin(t) * std::cos(p) * frame.normal
                                + std::sin(t) * std::sin(p) * frame.starboard;
                        direction = direction.normalized();
                        const Vec3 observer = states[s].position + range * direction;
                        const Vec3 delta = observer - states[s].position;
                        const double squared_range = delta.squaredNorm();
                        if (squared_range <= radius * radius)
                            continue;
                        const int bin = DirectionBounds::bin(frame.nose.dot(delta),
                                frame.normal.dot(delta), frame.starboard.dot(delta), squared_range);
                        const double level = db.spectral_level_at_observer(observer, states[s],
                                no_absorption.data(), gains);
                        checked++;
                        if (bounds.numerator(s, bin) < std::pow(10.0, 0.1 * level) * squared_range)
                            exceeded++;
                    }
    }
    INFO(exceeded << " of " << checked << " directions exceed the spectral bound");
    REQUIRE(checked > 2000000);
    REQUIRE(exceeded == 0);
}
