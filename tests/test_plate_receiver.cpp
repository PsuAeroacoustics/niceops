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

// board_level(bands, hs, hypot(dx, dy), 340 m/s, table, source_dx=dx, source_dy=dy).
const std::vector<Reference> kGrass = {
{150.0, 12.0, -40.0, {5.982314669, 5.974559440, 5.964088689, 5.952453940, 5.938176403, 5.919795112, 5.895755728, 5.867188056, 5.829334617, 5.778462242, 5.716743088, 5.637388922, 5.524220051, 5.395938602, 5.244783422, 5.078624663, 4.942168627, 4.936068029, 5.164064062, 5.616543974, 5.945076492, 5.729630480, 5.198984104, 6.479529887, 6.965905199, 6.709795501, 6.390205424, 5.809442055, 5.967687512, 6.071220593, 5.630780581}},
{100.0, -300.0, 250.0, {5.908537732, 5.880086360, 5.842059932, 5.800784636, 5.751798267, 5.691430605, 5.616668981, 5.533166848, 5.429990891, 5.302288030, 5.161055948, 4.997543545, 4.793518086, 4.596886204, 4.402344990, 4.211094010, 3.992542293, 3.673837801, 3.133564702, 2.599338368, 2.393517072, 1.749824331, 1.458903417, 1.000284946, 0.798358884, 0.752635957, 1.013722735, 1.467743072, 2.186313787, 3.042483213, 3.949203997}},
{60.0, 500.0, -1900.0, {5.816695641, 5.719608171, 5.561891273, 5.355062691, 5.063642709, 4.642474963, 4.043669041, 3.308929593, 2.372566490, 1.277134879, 0.229895993, -0.783666678, -1.887277933, -2.916829256, -4.006345006, -5.210658230, -6.532956842, -7.805068071, -9.065624022, -10.126693795, -10.679119135, -10.778448222, -10.583952903, -10.467733320, -10.505768295, -10.626226680, -10.824121007, -9.620496952, -8.914154388, -7.460990141, -6.900912114}},
{300.0, 0.0, 5.0, {5.982690994, 5.975249901, 5.965172087, 5.953943085, 5.940132344, 5.922319481, 5.898996243, 5.871278658, 5.834615826, 5.785582440, 5.726665196, 5.652195098, 5.549288534, 5.438879378, 5.321305882, 5.215992065, 5.175121461, 5.257371813, 5.459217273, 5.538964550, 5.172347167, 4.841412710, 5.749179055, 5.784122757, 6.486293730, 6.712008324, 6.869881406, 6.384932466, 5.706746941, 6.163134839, 5.960332778}},
{40.0, 900.0, 3.0, {5.889772001, 5.828351986, 5.729352767, 5.600550601, 5.420269113, 5.160662277, 4.790113714, 4.326940504, 3.710325653, 2.927231016, 2.091930337, 1.207640708, 0.221827139, -0.666417744, -1.557268712, -2.474850672, -3.399833890, -4.257529482, -5.244557200, -6.494683015, -7.593216077, -8.225932882, -9.294596700, -9.712121646, -10.195567020, -10.198295375, -9.926907666, -9.385437167, -8.604460633, -7.679207779, -6.647563759}},
};
const std::vector<Reference> kAmedee = {
{150.0, 12.0, -40.0, {5.947704553, 5.937651160, 5.924863753, 5.911455610, 5.895866140, 5.876841968, 5.853258550, 5.826603740, 5.792896503, 5.749579906, 5.699071923, 5.636286142, 5.549363586, 5.453000091, 5.340830956, 5.217283532, 5.111507338, 5.093647798, 5.242604668, 5.578637882, 5.858587165, 5.729174748, 5.187929626, 6.144580093, 6.863886323, 6.655416763, 6.621686655, 5.911297738, 5.758601560, 6.197878284, 5.545034275}},
{100.0, -300.0, 250.0, {5.796764738, 5.756445152, 5.706517247, 5.656233010, 5.600571942, 5.536552535, 5.462562915, 5.385171703, 5.295376399, 5.191121511, 5.082855093, 4.965238426, 4.829252063, 4.710581642, 4.609127499, 4.530260638, 4.450577203, 4.294455694, 3.964778190, 3.657048782, 3.626545976, 3.208881360, 3.047609549, 2.662639618, 2.395480106, 2.115204599, 1.916126964, 1.852998215, 1.992754395, 2.412181002, 3.065787227}},
{60.0, 500.0, -1900.0, {5.475386827, 5.289752489, 5.016543993, 4.692653319, 4.280436761, 3.747531903, 3.079139134, 2.360299489, 1.556845986, 0.720277055, -0.024990849, -0.729072425, -1.487487539, -2.185503025, -2.917001206, -3.723271236, -4.613093834, -5.481895449, -6.370140922, -7.175524858, -7.684900170, -7.873853284, -7.892165558, -8.035167419, -8.418444525, -8.991265840, -9.917563190, -9.432868853, -8.793576010, -7.826193540, -7.335287201}},
{300.0, 0.0, 5.0, {5.948751471, 5.939210641, 5.927006411, 5.914149730, 5.899147470, 5.880788444, 5.857988131, 5.832209725, 5.799669715, 5.758077708, 5.710101309, 5.651590341, 5.573392133, 5.491889156, 5.407295102, 5.333861489, 5.309564339, 5.377278020, 5.539220590, 5.616625130, 5.319879715, 4.940732088, 5.581198482, 5.650026591, 6.198783288, 6.553859081, 6.930404408, 6.624939523, 5.747899513, 5.975984458, 6.107974660}},
{40.0, 900.0, 3.0, {5.669938477, 5.552583505, 5.380651381, 5.177341370, 4.918246809, 4.580570170, 4.148330949, 3.665547153, 3.092405369, 2.443303757, 1.813781503, 1.187633673, 0.513708677, -0.079049445, -0.658440617, -1.238138786, -1.805560823, -2.325769733, -2.952945580, -3.821497023, -4.621415466, -5.092436458, -6.005803215, -6.528481400, -7.264553245, -7.763169713, -8.176891847, -8.368503099, -8.259011554, -7.888906517, -7.191103139}},
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
