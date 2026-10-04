#include <catch2/catch_approx.hpp>
// The ground-plane receiver (plate_receiver.hpp), pinned against panam's
// axisymmetric_bem.board_level() evaluated on the same netCDF tables
// (plates/, written by axisymmetric_bem.write_netcdf), and its level bound
// checked over random geometries.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "plate_receiver.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

const std::vector<double> kBands = {10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200,
        250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
        10000};

struct Reference
{
    double source_height;
    double dx; // source minus microphone, m (panam's source_dx, source_dy)
    double dy;
    std::vector<double> level;
};

GroundModel grass()
{
    GroundModel ground;
    ground.flow_resistance = 225.0;
    return ground;
}

GroundModel amedee()
{
    GroundModel ground;
    ground.impedance_model = ImpedanceModel::variable_porosity;
    ground.flow_resistance = 200.0;
    ground.porosity_rate = 0.0;
    return ground;
}

// board_level(bands, hs, hypot(dx, dy), 340 m/s, table, source_dx=dx, source_dy=dy), panam b758e3a.
const std::vector<Reference> kGrass = {
{150.0, 12.0, -40.0, {5.982314669, 5.974559440, 5.964088689, 5.952453940, 5.938176403, 5.919795112, 5.895755728, 5.867188056, 5.829334617, 5.778462242, 5.716743088, 5.637388922, 5.524220051, 5.395938602, 5.244783422, 5.078624663, 4.942663444, 4.936503334, 5.164438593, 5.616857393, 5.945334496, 5.729833361, 5.199121493, 6.479613867, 6.965939271, 6.709760761, 6.390089769, 5.809218696, 5.967395188, 6.070952556, 5.630730367}},
{100.0, -300.0, 250.0, {5.908537732, 5.880086360, 5.842059932, 5.800784636, 5.751798267, 5.691430605, 5.616668981, 5.533166848, 5.429990891, 5.302288030, 5.161055948, 4.997543545, 4.793518086, 4.596886204, 4.402344990, 4.211094010, 3.992542293, 3.673837801, 3.133564702, 2.599338368, 2.393989718, 1.750005717, 1.458720399, 0.999880669, 0.797814872, 0.752000232, 1.013034680, 1.467031728, 2.185619359, 3.041869705, 3.948730738}},
{60.0, 500.0, -1900.0, {5.816695641, 5.719608171, 5.561891273, 5.355062691, 5.063642709, 4.642474963, 4.043669041, 3.308929593, 2.372566490, 1.277134879, 0.229895993, -0.783666678, -1.887277933, -2.916829256, -4.006345006, -5.210658230, -6.532956842, -7.805068071, -9.065624022, -10.131088293, -10.699061194, -10.793139046, -10.594427008, -10.475518845, -10.511649391, -10.630737491, -10.827591683, -9.623198723, -8.916394208, -7.462651312, -6.902129390}},
{300.0, 0.0, 5.0, {5.982690994, 5.975249901, 5.965172087, 5.953943085, 5.940132344, 5.922319481, 5.898996243, 5.871278658, 5.834615826, 5.785582440, 5.726665196, 5.652195098, 5.549355048, 5.439198499, 5.321596625, 5.216253896, 5.175353786, 5.257576901, 5.459394619, 5.539113755, 5.172469683, 4.841507007, 5.749244687, 5.784162491, 6.486306425, 6.711991228, 6.869827041, 6.384827443, 5.706586204, 6.163017400, 5.960309858}},
{40.0, 900.0, 3.0, {5.889772001, 5.828351986, 5.729352767, 5.600550601, 5.420269113, 5.160662277, 4.790113714, 4.326940504, 3.710325653, 2.927231016, 2.091930337, 1.207640708, 0.221827139, -0.666417744, -1.557268712, -2.474850672, -3.399833890, -4.257529482, -5.244557200, -6.494683015, -7.593216077, -8.239805970, -9.308068532, -9.722328842, -10.203400541, -10.204352773, -9.931654970, -9.389287425, -8.607577352, -7.681611931, -6.649292582}},
};
const std::vector<Reference> kAmedee = {
{150.0, 12.0, -40.0, {5.947704553, 5.937651160, 5.924863753, 5.911455610, 5.895866140, 5.876841968, 5.853258550, 5.826603740, 5.792896503, 5.749579906, 5.699071923, 5.636286142, 5.549363586, 5.453000091, 5.340830956, 5.217283532, 5.111868777, 5.093959738, 5.242870259, 5.578860656, 5.858774168, 5.729328748, 5.188046758, 6.144668615, 6.863951557, 6.655449555, 6.621684231, 5.911233243, 5.758439296, 6.197561723, 5.544900559}},
{100.0, -300.0, 250.0, {5.796764738, 5.756445152, 5.706517247, 5.656233010, 5.600571942, 5.536552535, 5.462562915, 5.385171703, 5.295376399, 5.191121511, 5.082855093, 4.965238426, 4.829252063, 4.710581642, 4.609127499, 4.530260638, 4.450577203, 4.294455694, 3.964778190, 3.657048782, 3.626545976, 3.209559109, 3.048134501, 2.662887322, 2.395507966, 2.115046257, 1.915813469, 1.852560500, 1.992202604, 2.411525564, 3.065070769}},
{60.0, 500.0, -1900.0, {5.475386827, 5.289752489, 5.016543993, 4.692653319, 4.280436761, 3.747531903, 3.079139134, 2.360299489, 1.556845986, 0.720277055, -0.024990849, -0.729072425, -1.487487539, -2.185503025, -2.917001206, -3.723271236, -4.613093834, -5.481895449, -6.370140922, -7.175524858, -7.684900170, -7.883948297, -7.905793271, -8.045887820, -8.426912132, -8.997975724, -9.922904285, -9.437072174, -8.797115437, -7.829121764, -7.337813896}},
{300.0, 0.0, 5.0, {5.948751471, 5.939210641, 5.927006411, 5.914149730, 5.899147470, 5.880788444, 5.857988131, 5.832209725, 5.799669715, 5.758077708, 5.710101309, 5.651590341, 5.573447018, 5.492145915, 5.407520543, 5.334057734, 5.309733448, 5.377424278, 5.539345530, 5.616730284, 5.319967674, 4.940803130, 5.581253299, 5.650068032, 6.198811207, 6.553873512, 6.930402633, 6.624912720, 5.747814750, 5.975813232, 6.107937303}},
{40.0, 900.0, 3.0, {5.669938477, 5.552583505, 5.380651381, 5.177341370, 4.918246809, 4.580570170, 4.148330949, 3.665547153, 3.092405369, 2.443303757, 1.813781503, 1.187633673, 0.513708677, -0.079049445, -0.658440617, -1.238138786, -1.805560823, -2.325769733, -2.952945580, -3.821497023, -4.621415466, -5.092436458, -6.008383400, -6.541059385, -7.274840388, -7.771529387, -8.183681551, -8.374161453, -8.263789370, -7.893000085, -7.194611574}},
};

void check_against(const std::string& path, const GroundModel& ground,
        const std::vector<Reference>& references)
{
    const auto plate = PlateReceiver::load(path).bound(kBands, ground, 340.0);
    std::vector<double> level(kBands.size());
    for (const auto& r : references)
    {
        // NICE-OPS propagates observer minus source, the reverse of panam's dx, dy.
        plate.level_exact(r.source_height, std::hypot(r.dx, r.dy), std::atan2(-r.dy, -r.dx),
                level.data(), 0, kBands.size());
        for (std::size_t band = 0; band < kBands.size(); band++)
        {
            INFO(path << " hs " << r.source_height << " dx " << r.dx << " dy " << r.dy << " band "
                      << kBands[band]);
            REQUIRE(std::abs(level[band] - r.level[band]) < 1e-6);
        }
    }
}

} // namespace

TEST_CASE("Plate receiver matches panam's board_level", "[plate_receiver]")
{
    check_against("plates/gr1425_db225.nc", grass(), kGrass);
    check_against("plates/gr1425_vp200.nc", amedee(), kAmedee);
}

TEST_CASE("Plate receiver's quadratic Q stays within 0.1 dB of the exact", "[plate_receiver]")
{
    const auto plate = PlateReceiver::load("plates/gr1425_db225.nc").bound(kBands, grass(), 343.0);
    std::mt19937_64 random(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> fast(kBands.size());
    std::vector<double> exact(kBands.size());
    double worst = 0.0;
    for (int trial = 0; trial < 20000; trial++)
    {
        const double hs = unit(random) < 0.15 ? 0.5 + 4.5 * unit(random)
                                              : 1.0 + 1000.0 * std::pow(unit(random), 3.0);
        const double d2 = 10000.0 * std::pow(unit(random), 3.0);
        const double azimuth = 2.0 * pi * unit(random);
        plate.level(hs, d2, azimuth, fast.data(), 0, kBands.size());
        plate.level_exact(hs, d2, azimuth, exact.data(), 0, kBands.size());
        for (std::size_t band = 0; band < kBands.size(); band++)
        {
            if (exact[band] > -15.0)
                worst = std::max(worst, std::abs(fast[band] - exact[band]));
        }
    }
    INFO("worst difference above -15 dB: " << worst << " dB");
    REQUIRE(worst < 0.1);
}

TEST_CASE("Plate receiver gain bounds hold over random geometries", "[plate_receiver]")
{
    for (const auto& [path, ground] : {std::pair {std::string("plates/gr1425_db225.nc"), grass()},
                 std::pair {std::string("plates/gr1425_vp200.nc"), amedee()}})
    {
        const auto plate = PlateReceiver::load(path).bound(kBands, ground, 343.0, 3, 29);
        const auto bounds = plate.gain_bounds(3, 29);
        REQUIRE(bounds[0] == 1.0);
        REQUIRE(bounds[30] == 1.0);
        std::mt19937_64 random(7);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        std::vector<double> level(kBands.size());
        double closest = 1e300; // smallest margin seen, dB
        for (int trial = 0; trial < 200000; trial++)
        {
            // Heights from grazing (0) up, ranges from overhead to 10 km.
            const double hs = trial % 7 == 0 ? 0.0 : 1000.0 * std::pow(unit(random), 3.0);
            const double d2 = trial % 11 == 0 ? 0.0 : 10000.0 * std::pow(unit(random), 3.0);
            if (hs == 0.0 && d2 == 0.0)
                continue;
            const double azimuth = 2.0 * pi * unit(random) - pi;
            for (const bool exact : {false, true})
            {
                if (exact)
                    plate.level_exact(hs, d2, azimuth, level.data(), 3, 29);
                else
                    plate.level(hs, d2, azimuth, level.data(), 3, 29);
                for (std::size_t band = 3; band < 29; band++)
                {
                    const double margin = 10.0 * std::log10(bounds[band]) - level[band];
                    closest = std::min(closest, margin);
                    REQUIRE(margin >= 0.0);
                }
            }
        }
        INFO(path << " tightest bound margin " << closest << " dB");
        // A bound, but not a vacuous one.  It takes Q's phase as adversarial
        // (|P_d| + |Q|max |P_r|, the fast path's Q widened by its quadratic's
        // weights), so it sits ~4 dB above the level reached --
        // the point receiver's (1 + |Q|)^2 is looser, 9.6 dB against +6.
        CHECK(closest < 5.0);
    }
}

TEST_CASE("Plate receiver refuses a table for another ground or other bands", "[plate_receiver]")
{
    const auto plate = PlateReceiver::load("plates/gr1425_db225.nc");
    auto other = grass();
    other.flow_resistance = 200.0;
    REQUIRE_THROWS(plate.bound(kBands, other, 343.0));
    REQUIRE_THROWS(plate.bound(kBands, amedee(), 343.0));
    auto shifted = kBands;
    shifted[10] = 101.0;
    REQUIRE_THROWS(plate.bound(shifted, grass(), 343.0));
    // A band outside the run's range may differ.
    REQUIRE_NOTHROW(plate.bound(shifted, grass(), 343.0, 11, 31));
}

TEST_CASE("Plate receiver's tabulated Q stays close to the evaluated one, inside the bound",
        "[plate_receiver]")
{
    const auto evaluated = PlateReceiver::load("plates/gr1425_vp200.nc").bound(kBands, amedee(), 343.0);
    auto tabulated = evaluated;
    tabulated.tabulate(12000.0, 161, 40);
    const auto bounds = tabulated.gain_bounds(0, kBands.size());
    std::mt19937_64 random(13);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<double> a(kBands.size()), b(kBands.size());
    double worst = 0.0, sum2 = 0.0;
    long n = 0;
    for (int trial = 0; trial < 30000; trial++)
    {
        const double r = 5.0 * std::pow(2000.0, unit(random));
        const double elevation = trial % 3 == 0 ? 0.2 * unit(random) * unit(random) : 1.5707 * unit(random);
        const double azimuth = 2.0 * pi * unit(random);
        evaluated.level(r * std::sin(elevation), r * std::cos(elevation), azimuth, a.data(), 0, kBands.size());
        tabulated.level(r * std::sin(elevation), r * std::cos(elevation), azimuth, b.data(), 0, kBands.size());
        for (std::size_t band = 0; band < kBands.size(); band++)
        {
            REQUIRE(10.0 * std::log10(bounds[band]) >= b[band]);
            if (a[band] > -15.0)
            {
                const double e = b[band] - a[band];
                worst = std::max(worst, std::abs(e));
                sum2 += e * e;
                n++;
            }
        }
    }
    // The worst cases sit where the model's own boundary-loss cutoff (|w| =
    // 500) makes Q jump, which a table smooths over.
    INFO("worst " << worst << " dB, rms " << std::sqrt(sum2 / n) << " dB");
    CHECK(std::sqrt(sum2 / n) < 0.002);
    CHECK(worst < 0.15);
}

namespace {

struct Placed
{
    Vec3 position;
};

} // namespace

TEST_CASE("the plate's azimuth is the ray's bearing turned by the microphone's offset", "[plate_receiver]")
{
    // The tables put the microphone on +y, so a ray travelling straight outward
    // along the offset arrives from the plate's inboard side: pi/2.  Nothing tested
    // this mapping; a dropped half-pi or a sign on the offset rotates every plate
    // result by up to a dB (the table is not symmetric under a quarter turn).
    constexpr double pi = 3.14159265358979323846;
    auto wrapped = [](double angle) { return std::remainder(angle, 2.0 * pi); };
    for (const double offset : {0.0, 0.7, pi / 2.0, -2.0, 3.0})
    {
        const Vec3 outward(100.0 * std::cos(offset), 100.0 * std::sin(offset), -30.0);
        CHECK(wrapped(plate_frame_azimuth(outward, offset) - 0.5 * pi) == Catch::Approx(0.0).margin(1e-12));
        // Travelling the other way, toward the track, is the opposite side.
        CHECK(std::fabs(wrapped(plate_frame_azimuth(-1.0 * outward, offset) + 0.5 * pi)) < 1e-12);
    }
    // Turning the geometry and the offset together changes nothing.
    const Vec3 separation(300.0, -120.0, -50.0);
    const double reference = plate_frame_azimuth(separation, 0.4);
    for (const double turn : {0.3, 1.9, -2.5})
    {
        const Vec3 turned(separation.x * std::cos(turn) - separation.y * std::sin(turn),
                separation.x * std::sin(turn) + separation.y * std::cos(turn), separation.z);
        CHECK(wrapped(plate_frame_azimuth(turned, 0.4 + turn) - reference) == Catch::Approx(0.0).margin(1e-12));
    }
    // It matters: the table's level is not the same a quarter turn away.
    const auto plate = PlateReceiver::load("plates/gr1425_db225.nc").bound(kBands, grass(), 340.0);
    std::vector<double> here(kBands.size()), turned(kBands.size());
    plate.level_exact(40.0, 900.0, 0.0, here.data(), 0, kBands.size());
    plate.level_exact(40.0, 900.0, 0.5 * pi, turned.data(), 0, kBands.size());
    double differ = 0.0;
    for (std::size_t band = 0; band < kBands.size(); band++)
        differ = std::max(differ, std::fabs(here[band] - turned[band]));
    CHECK(differ > 0.05);
}

TEST_CASE("outboard puts the microphone away from the track", "[plate_receiver]")
{
    constexpr double pi = 3.14159265358979323846;
    std::vector<Placed> track;
    for (int k = -10; k <= 10; k++)
        track.push_back({Vec3(10.0 * k, 0.0, 50.0)});
    // Off to either side: the direction from the track's nearest point to the observer.
    CHECK(outboard_azimuth(Vec3(0.0, 300.0, 0.0), track) == Catch::Approx(0.5 * pi).margin(1e-12));
    CHECK(outboard_azimuth(Vec3(0.0, -300.0, 0.0), track) == Catch::Approx(-0.5 * pi).margin(1e-12));
    CHECK(outboard_azimuth(Vec3(40.0, 120.0, 0.0), track) == Catch::Approx(std::atan2(120.0, 0.0)).margin(1e-12));
    // Beyond the track's end the nearest point is the end state.
    CHECK(outboard_azimuth(Vec3(400.0, 0.0, 0.0), track) == Catch::Approx(0.0).margin(1e-12));
    // Directly under it there is no direction: the track's right, which for a
    // track along +x is -y.
    CHECK(outboard_azimuth(Vec3(30.0, 0.0, 0.0), track) == Catch::Approx(-0.5 * pi).margin(1e-12));
    std::reverse(track.begin(), track.end());
    CHECK(outboard_azimuth(Vec3(30.0, 0.0, 0.0), track) == Catch::Approx(0.5 * pi).margin(1e-12));
}
