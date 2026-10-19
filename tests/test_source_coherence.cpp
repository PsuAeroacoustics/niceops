#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <vector>

#include "buildings.hpp"
#include "point_paths.hpp"
#include "source_coherence.hpp"
#include "sphere_database.hpp"
#include "urban_point.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT

namespace {

const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};
constexpr double c0 = 343.2;

Vec3 at_azimuth(double degrees, double down = 30.0)
{
    const double a = degrees * pi / 180.0, e = down * pi / 180.0;
    return Vec3(std::cos(e) * std::cos(a), std::cos(e) * std::sin(a), -std::sin(e));
}

// A level frame: the normal up, the nose along +x.
SphereFrame level_frame()
{
    return sphere_frame(Vec3(1.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
}

std::vector<double> band_levels_of(const std::vector<double>& energy)
{
    std::vector<double> out(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
        out[b] = band_level(bands[b], sub_bands(bands[b], c0), &energy[3 * b]);
    return out;
}

std::vector<double> sum(const PointSourceField& field, const SourceCoherence& coherence, double floor = 1.0,
        std::vector<double>* variance = nullptr)
{
    const std::vector<double> reflection; // no facades
    std::vector<double> energy(3 * bands.size());
    if (variance != nullptr)
        variance->assign(bands.size(), 0.0);
    partially_coherent_energies(field, reflection, floor, bands, 0, c0, {}, nullptr, coherence, level_frame(),
            energy.data(), variance != nullptr ? variance->data() : nullptr);
    return energy;
}

SourceCoherence with(RaySum how)
{
    SourceCoherence c = default_source_coherence(43.25, 4.0, 1);
    c.sum = how;
    return c;
}

} // namespace

// A pattern turning counterclockwise (seen from above) at Omega reaches a
// direction further round that way dpsi / Omega later: that ray's path is
// longer by c dpsi / Omega.  The difference is wrapped, so directions either
// side of the nose's opposite are a small angle apart, not nearly a turn.
TEST_CASE("the rotor's delay is the time its pattern takes to turn between two directions", "[source_coherence]")
{
    SourceCoherence c = default_source_coherence(40.0, 4.0, 1);
    const SphereFrame f = level_frame();
    REQUIRE(c.delay_length(at_azimuth(30.0), at_azimuth(0.0), f.normal, f.nose, c0)
            == Catch::Approx(c0 * (pi / 6.0) / 40.0));
    REQUIRE(c.delay_length(at_azimuth(0.0), at_azimuth(30.0), f.normal, f.nose, c0)
            == Catch::Approx(-c0 * (pi / 6.0) / 40.0));
    REQUIRE(c.delay_length(at_azimuth(179.0), at_azimuth(-179.0), f.normal, f.nose, c0)
            == Catch::Approx(-c0 * (2.0 * pi / 180.0) / 40.0));
    c.rotation_sense = -1;
    REQUIRE(c.delay_length(at_azimuth(30.0), at_azimuth(0.0), f.normal, f.nose, c0)
            == Catch::Approx(-c0 * (pi / 6.0) / 40.0));
    // Only the azimuth about the shaft counts.
    REQUIRE(c.delay_length(at_azimuth(30.0, 10.0), at_azimuth(30.0, 70.0), f.normal, f.nose, c0)
            == Catch::Approx(0.0).margin(1e-12));
    // The rotation's band: 4.5 blade-passage frequencies.
    REQUIRE(c.rotation_limit == Catch::Approx(4.5 * 4.0 * 40.0 / (2.0 * pi)));
    REQUIRE(c.rotates(100.0));
    REQUIRE(!c.rotates(125.0));
    // Without the blades or the sense there is no delay, and the decay is the
    // fit made without it.
    const SourceCoherence none = default_source_coherence(40.0, 4.0, 0);
    REQUIRE(none.rotation_rate == 0.0);
    REQUIRE(none.scale == SourceCoherence::default_scale_without_rotation);
    REQUIRE(none.delay_length(at_azimuth(30.0), at_azimuth(0.0), f.normal, f.nose, c0) == 0.0);
}

TEST_CASE("a pair's coherency is one at no separation, and its spread falls to the band's distant coherence",
        "[source_coherence]")
{
    const SourceCoherence c = default_source_coherence();
    for (double f : {25.0, 250.0, 2500.0})
    {
        const double t = c.tonal(99, f);
        REQUIRE(c.mean(f, 0.0) == 1.0);
        REQUIRE(c.coherence(f, 0.0, t) == 1.0);
        REQUIRE(c.spread(f, 0.0, t) == 0.0);
        REQUIRE(c.mean(f, 179.0) < 1e-6);
        // Toward the band's distant coherence as frequency times angle grows.
        REQUIRE(c.spread(f, 179.0, t) >= t);
        REQUIRE(c.spread(f, 179.0, t) < c.spread(f, 30.0, t));
        REQUIRE(c.spread(f, 1e6 / f, t) == Catch::Approx(t).margin(1e-12));
        // Never more than the coherence leaves above the mean.
        for (double degrees : {0.5, 2.0, 8.0, 30.0})
            REQUIRE(c.spread(f, degrees, t) <= 1.0 - c.mean(f, degrees) * c.mean(f, degrees) + 1e-15);
    }
    REQUIRE(c.tonal(99, 58.0) == Catch::Approx(0.5));
    SourceCoherence given = c;
    given.tonal_fraction = {0.25, 2.0};
    REQUIRE(given.tonal(0, 58.0) == 0.25);
    REQUIRE(given.tonal(1, 58.0) == 1.0);
    REQUIRE(given.tonal(2, 58.0) == Catch::Approx(0.5));
}

// The coherent sum, pair by pair here, is band_levels's, and rays that leave
// the source together -- a flush microphone's direct ray and the floor's,
// which leave toward the same point -- stay coherent under the measured
// model, the rotor's delay or not.
TEST_CASE("rays that leave the source together add with their phases", "[source_coherence]")
{
    const Buildings none;
    const PointSourceField flush(none, Vec3(0.0, 0.0, 60.0), Vec3(120.0, 30.0, 0.0));
    REQUIRE(flush.rays().size() == 2);
    std::vector<double> reference(bands.size());
    flush.band_levels({}, 1.0, bands, c0, reference.data());
    const auto coherent = band_levels_of(sum(flush, with(RaySum::coherent)));
    const auto measured = band_levels_of(sum(flush, with(RaySum::measured)));
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        REQUIRE(coherent[b] == Catch::Approx(reference[b]).margin(1e-9));
        REQUIRE(measured[b] == Catch::Approx(reference[b]).margin(1e-9));
    }
    std::vector<double> variance;
    sum(flush, with(RaySum::measured), 1.0, &variance);
    for (double v : variance)
        REQUIRE(v == 0.0);
}

// A microphone 1 m up, 300 m out from a source 100 m up: the floor's ray
// leaves 0.35 degrees below the direct one.  At 125 Hz that is close to
// coherent; at 4 kHz the two nearly add in energy.  (Not in a dip of the
// interference, where a coherency a little under one fills it by decibels.)
TEST_CASE("rays leaving apart lose their coherence as frequency times angle grows", "[source_coherence]")
{
    const Buildings none;
    const PointSourceField field(none, Vec3(0.0, 0.0, 100.0), Vec3(300.0, 0.0, 1.0));
    REQUIRE(field.rays().size() == 2);
    const auto coherent_parts = sum(field, with(RaySum::coherent)), energy_parts = sum(field, with(RaySum::energy)),
               measured_parts = sum(field, with(RaySum::measured));
    const auto coherent = band_levels_of(coherent_parts), energy = band_levels_of(energy_parts),
               measured = band_levels_of(measured_parts);
    REQUIRE(std::abs(measured[0] - coherent[0]) < 0.25);
    // At 4 kHz the pair keeps exp(-(1400 / 850)^1.6) = 0.11 of its cross term.
    REQUIRE(std::abs(measured[5] - energy[5]) < 0.15 * std::abs(coherent[5] - energy[5]));
    // One pair at one azimuth: each sub-band keeps a fraction of its cross
    // term, so lies between the two.  (A band need not: its sub-bands keep
    // different fractions of cross terms of either sign.)
    for (std::size_t p = 0; p < measured_parts.size(); p++)
    {
        const double lo = std::min(coherent_parts[p], energy_parts[p]), hi = std::max(coherent_parts[p], energy_parts[p]);
        REQUIRE(measured_parts[p] >= lo * (1.0 - 1e-12));
        REQUIRE(measured_parts[p] <= hi * (1.0 + 1e-12));
    }
}

// Two rays in different groups: the band's variance is half the square of
// their cross term, summed over its sub-bands at their widths, times their
// spread.
TEST_CASE("the band's variance is the pair's cross term squared times its spread", "[source_coherence]")
{
    const Buildings none;
    const PointSourceField field(none, Vec3(0.0, 0.0, 30.0), Vec3(40.0, 10.0, 6.0));
    REQUIRE(field.rays().size() == 2);
    const SourceCoherence model = with(RaySum::measured);
    std::vector<double> variance;
    sum(field, model, 1.0, &variance);
    std::vector<Vec3> launch;
    for (const PointRay& r : field.rays())
        launch.push_back(point_ray_launch(r));
    const double degrees = (launch[0] - launch[1]).norm() * 180.0 / pi;
    std::vector<std::complex<double>> cross(bands.size());
    std::vector<std::size_t> order;
    field.sub_band_problems(
            {}, 1.0, bands, c0,
            [&](std::size_t band, int, const std::vector<std::complex<double>>& a, const std::vector<double>& l,
                    const SubBand& s) {
                const double mid = 0.5 * (s.k1 + s.k2), half = 0.5 * (s.k2 - s.k1);
                const double f1 = bands[band] / std::pow(2.0, 1.0 / 6.0), f2 = bands[band] * std::pow(2.0, 1.0 / 6.0);
                double dl = l[0] - l[1];
                // The rays sorted by length; the delay is the first's (in ray order) over the second's.
                const SphereFrame f = level_frame();
                const double delay = model.delay_length(launch[order[0]], launch[order[1]], f.normal, f.nose, c0);
                if (model.rotates(bands[band]))
                    dl += delay;
                const double x = half * dl;
                cross[band] += (s.g2 - s.g1) / (f2 - f1) * 2.0 * std::sin(x) / x * a[0] * std::conj(a[1])
                        * std::polar(1.0, mid * dl);
            },
            {}, nullptr, &order);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const double expected
                = 0.5 * std::norm(cross[b]) * model.spread(bands[b], degrees, model.tonal(b, bands[b]));
        REQUIRE(variance[b] == Catch::Approx(expected).epsilon(1e-9));
    }
    REQUIRE(variance[0] > 0.0);
}

// The database's record first; its fits come with its own rotor, so a record
// that holds fits but no rotor keeps the vehicle file's out; what the record
// leaves out is the default for the rotation it ends with.
TEST_CASE("a database's coherence comes before the vehicle file's rotor and the defaults", "[source_coherence]")
{
    const CoherenceRecord none;
    const SourceCoherence vehicle = source_coherence(none, 40.0, 4.0, -1);
    REQUIRE(vehicle.rotation_rate == 40.0);
    REQUIRE(vehicle.rotation_sense == -1);
    REQUIRE(vehicle.scale == SourceCoherence::default_scale);
    REQUIRE(vehicle.source == "the default coherence, the vehicle file's rotor");

    CoherenceRecord rotor;
    rotor.blades = 2.0;
    rotor.rotation_sense = 1;
    const SourceCoherence database_rotor = source_coherence(rotor, 40.0, 4.0, -1);
    REQUIRE(database_rotor.rotation_sense == 1);
    REQUIRE(database_rotor.rotation_limit == Catch::Approx(4.5 * 2.0 * 40.0 / (2.0 * pi)));
    REQUIRE(database_rotor.scale == SourceCoherence::default_scale);

    CoherenceRecord fits;
    fits.scale = 600.0;
    fits.exponent = 2.0;
    fits.spread_scale = 1500.0;
    fits.spread_exponent = 0.8;
    fits.tonal_fraction = {0.9, 0.5, 0.1};
    fits.source = "a hexacopter, 2023";
    const SourceCoherence multirotor = source_coherence(fits, 40.0, 4.0, -1);
    REQUIRE(multirotor.rotation_rate == 0.0);
    REQUIRE(multirotor.scale == 600.0);
    REQUIRE(multirotor.exponent == 2.0);
    REQUIRE(multirotor.spread_scale == 1500.0);
    REQUIRE(multirotor.spread_exponent == 0.8);
    REQUIRE(multirotor.tonal(1, 1000.0) == 0.5);
    REQUIRE(multirotor.source == "the database's coherence (a hexacopter, 2023), no rotor delay");

    // Where it gives none of the fits, they are the defaults for no rotation.
    CoherenceRecord tonal_only;
    tonal_only.tonal_fraction = {0.3};
    const SourceCoherence partial = source_coherence(tonal_only, 40.0, 4.0, -1);
    REQUIRE(partial.rotation_rate == 0.0);
    REQUIRE(partial.scale == SourceCoherence::default_scale_without_rotation);

    fits.blades = 4.0;
    fits.rotation_sense = 1;
    fits.rotation_limit = 100.0;
    const SourceCoherence helicopter = source_coherence(fits, 40.0, 0.0, 0);
    REQUIRE(helicopter.rotation_rate == 40.0);
    REQUIRE(helicopter.rotation_limit == 100.0);
    REQUIRE(helicopter.scale == 600.0);
}

// Written at a database's root and read back; an exported sphere carries it,
// as it carries the rest of the root.  A tonal fraction that does not match
// the bands is refused.
TEST_CASE("a database carries its coherence through an exported sphere", "[source_coherence]")
{
    Database plain("databases/S-76D_M3.nod");
    REQUIRE(plain.coherence_record().empty());
    auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    sphere.frequencies = {125.0, 500.0, 2000.0};
    sphere.amplitudes.assign(sphere.cells.size() * sphere.frequencies.size(), 90.0);
    CoherenceRecord& r = sphere.coherence;
    r.blades = 4.0;
    r.rotation_sense = -1;
    r.rotation_limit = 90.0;
    r.scale = 700.0;
    r.exponent = 1.5;
    r.spread_scale = 1800.0;
    r.spread_exponent = 1.1;
    r.tonal_fraction = {0.8, 0.2, 0.0};
    r.source = "test";
    const auto path = (std::filesystem::temp_directory_path() / "niceops_coherence_record.nod").string();
    write_sphere(sphere, path);
    const Database back(path);
    const CoherenceRecord& b = back.coherence_record();
    REQUIRE(b.blades == r.blades);
    REQUIRE(b.rotation_sense == r.rotation_sense);
    REQUIRE(b.rotation_limit == r.rotation_limit);
    REQUIRE(b.scale == r.scale);
    REQUIRE(b.exponent == r.exponent);
    REQUIRE(b.spread_scale == r.spread_scale);
    REQUIRE(b.spread_exponent == r.spread_exponent);
    REQUIRE(b.tonal_fraction == r.tonal_fraction);
    REQUIRE(b.source == r.source);
    REQUIRE(back.interpolate_sphere(0.15, -4.0, 0.006).coherence.tonal_fraction == r.tonal_fraction);

    r.tonal_fraction = {0.5};
    write_sphere(sphere, path);
    REQUIRE_THROWS_AS(Database(path), InputError);
    r.tonal_fraction.clear();
    r.exponent.reset();
    write_sphere(sphere, path);
    REQUIRE_THROWS_AS(Database(path), InputError);
}

// A flush microphone beside a wall: the direct ray and the floor's leave
// together, and so do the wall's reflection and its floor partner.  The two
// groups' four cross pairs deviate as one, so the variance is half the square
// of their sum -- four times what treating each pair alone gives when they
// are alike.
TEST_CASE("rays that leave together deviate together", "[source_coherence]")
{
    Footprint wall;
    wall.x = {40.0, 41.0, 41.0, 40.0};
    wall.y = {-50.0, -50.0, 50.0, 50.0};
    wall.top = 200.0;
    const Buildings city({wall}, 1.0 * pi / 180.0);
    PointPathOptions options;
    options.diffraction = false;
    const PointSourceField field(city, Vec3(0.0, 0.0, 30.0), Vec3(20.0, 5.0, 0.0), options);
    REQUIRE(field.rays().size() == 4);
    const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
    // The rotor's delay in the two lowest bands, and the nose turned so that
    // the azimuth's cut lies between the direct ray's (14.0 deg) and the
    // wall's reflection's (4.8 deg): their wrapped difference is the one-turn
    // shift's.
    SourceCoherence model = with(RaySum::measured);
    model.rotation_limit = 300.0;
    const double cut = 9.4 * pi / 180.0;
    const SphereFrame f = sphere_frame(Vec3(-std::cos(cut), -std::sin(cut), 0.0), Vec3(0.0, 0.0, 1.0));
    std::vector<double> energy(3 * bands.size()), variance(bands.size());
    partially_coherent_energies(field, reflection, 1.0, bands, 0, c0, {}, nullptr, model, f, energy.data(),
            variance.data());
    std::vector<Vec3> launch;
    for (const PointRay& r : field.rays())
        launch.push_back(point_ray_launch(r));
    // Which rays leave together.
    std::vector<int> group(4);
    for (std::size_t i = 0; i < 4; i++)
        group[i] = (launch[i] - launch[0]).norm() < 1e-9 ? 0 : 1;
    REQUIRE(std::count(group.begin(), group.end(), 0) == 2);
    const Vec3 left = f.normal.cross(f.nose);
    auto azimuth = [&](const Vec3& d) { return std::atan2(left.dot(d), f.nose.dot(d)); };
    std::size_t other = 0;
    while (group[other] == 0)
        other++;
    REQUIRE(std::abs(azimuth(launch[0]) - azimuth(launch[other])) > pi);
    std::vector<std::complex<double>> together(bands.size());
    std::vector<double> apart(bands.size(), 0.0);
    std::vector<std::size_t> order;
    field.sub_band_problems(
            reflection, 1.0, bands, c0,
            [&](std::size_t band, int, const std::vector<std::complex<double>>& a, const std::vector<double>& l,
                    const SubBand& s) {
                const double mid = 0.5 * (s.k1 + s.k2), half = 0.5 * (s.k2 - s.k1);
                const double f1 = bands[band] / std::pow(2.0, 1.0 / 6.0), f2 = bands[band] * std::pow(2.0, 1.0 / 6.0);
                const double w = (s.g2 - s.g1) / (f2 - f1);
                for (std::size_t p = 0; p < 4; p++)
                    for (std::size_t q = 0; q < 4; q++)
                    {
                        const std::size_t i = order[p], j = order[q];
                        if (group[i] != 0 || group[j] != 1)
                            continue;
                        double dl = l[p] - l[q];
                        if (model.rotates(bands[band]))
                            dl += model.delay_length(launch[i], launch[j], f.normal, f.nose, c0);
                        const double x = half * dl;
                        // The band has averaged a pair this far apart away, as
                        // band_averaged_energy leaves it out.
                        if (std::abs(x) > 30.0)
                            continue;
                        const double sigma = std::sqrt(model.spread(bands[band], (launch[i] - launch[j]).norm() * 180.0 / pi,
                                model.tonal(band, bands[band])));
                        const std::complex<double> c
                                = w * 2.0 * sigma * std::sin(x) / x * a[p] * std::conj(a[q]) * std::polar(1.0, mid * dl);
                        together[band] += c;
                        apart[band] += std::norm(c);
                    }
            },
            {}, nullptr, &order);
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(variance[b] == Catch::Approx(0.5 * std::norm(together[b])).epsilon(1e-9));
    // And the mean energy, each pair's cross term at its mean coherency.
    std::vector<double> expected(3 * bands.size(), 0.0);
    field.sub_band_problems(
            reflection, 1.0, bands, c0,
            [&](std::size_t band, int part, const std::vector<std::complex<double>>& a, const std::vector<double>& l,
                    const SubBand& s) {
                const double mid = 0.5 * (s.k1 + s.k2), half = 0.5 * (s.k2 - s.k1);
                double e = 0.0;
                for (std::size_t p = 0; p < 4; p++)
                {
                    e += std::norm(a[p]);
                    for (std::size_t q = p + 1; q < 4; q++)
                    {
                        const std::size_t i = order[p], j = order[q];
                        double dl = l[p] - l[q];
                        if (model.rotates(bands[band]))
                            dl += model.delay_length(launch[i], launch[j], f.normal, f.nose, c0);
                        const double x = half * dl;
                        if (std::abs(x) > 30.0)
                            continue;
                        const double sinc = std::abs(x) < 1e-8 ? 1.0 : std::sin(x) / x;
                        e += 2.0 * model.mean(bands[band], (launch[i] - launch[j]).norm() * 180.0 / pi) * sinc
                                * (a[p] * std::conj(a[q]) * std::polar(1.0, mid * dl)).real();
                    }
                }
                expected[3 * band + static_cast<std::size_t>(part)] = e;
            },
            {}, nullptr, &order);
    for (std::size_t q = 0; q < expected.size(); q++)
        REQUIRE(energy[q] == Catch::Approx(expected[q]).epsilon(1e-9));
    // At 125 Hz the four are nearly alike, and the groups count nearly four times.
    REQUIRE(variance[0] > 3.0 * 0.5 * apart[0]);
}
