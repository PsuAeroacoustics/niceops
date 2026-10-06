#pragma once

// The parabolic-equation propagation table (niceops_propagation_table --method
// pe; propagation_table.hpp is the format), as a library function so the tests
// can build one.
//
// Each cell is a band's level re free field, at a (bearing, source height,
// receiver height, range), through one effective-sound-speed profile per
// bearing.  Two ways to fill it:
//
//   source    one march per bearing and SOURCE height, every receiver height
//             read from the field.  What the tool always did.
//   receiver  one march per bearing and RECEIVER height, the source put at the
//             receiver's height and every source height read from the field:
//             reciprocity.  In the effective-sound-speed model the field from
//             z_s to z_r along a bearing is the field from z_r to z_s along the
//             reverse bearing with the wind reversed, and the effective
//             profile c + w.h is the same for both, so it is the same march
//             with the source and the receiver exchanged.  A table has a few
//             receiver heights and many source heights, so this is several
//             times fewer marches.  The default.
//
// Reciprocity holds on the grid only if the source is placed as the adjoint
// of the receivers' read.  Inside the grid it always was; within three grid
// points of the ground it is now too (pe.cpp, place_ground_source), which is
// what lets a flush microphone, z = 0, be a source.  Over impedance ground
// and through a refracting profile the two marches agree to a few thousandths
// of a dB; what is left is the domain's top, which follows the highest height
// in each march.
//
// Turbulence (gamma_t > 0) is the PE's phase screens, each table cell the
// energy mean over `samples` frequencies times `realizations` realizations.
// The ensemble mean is reciprocal as the still field is (the screens are
// statistically the same read in either direction), but one realization is
// not, so the two marches agree only to sampling noise.  Every march has its
// own seeds.

#include <complex>
#include <functional>
#include <string>
#include <vector>

namespace niceops {

enum class PeTableMarch
{
    source,
    receiver,
};

struct PeTableRequest
{
    std::vector<double> bearings; // degrees, only for progress messages
    // The effective sound speed along each bearing, at heights z.
    std::vector<double> z;
    std::vector<std::vector<double>> c_by_bearing;
    std::vector<double> source_heights; // m
    std::vector<double> receiver_heights; // m
    double max_range = 0.0; // m
    double step = 0.0; // m
    std::vector<double> bands; // Hz
    std::vector<std::complex<double>> admittance; // per band
    int samples = 5;
    double max_angle = 50.0;
    double points_per_wavelength = 0.0;
    PeTableMarch march = PeTableMarch::receiver;
    // Turbulence, as PeSetup: off at gamma_t 0.
    double gamma_t = 0.0;
    double outer_scale = 10.0;
    int realizations = 1;
    unsigned seed = 1;
    std::vector<double> turbulence_z;
    std::vector<double> turbulence_gamma;
    unsigned threads = 0;
};

// The table's values, [((bearing * ns + source) * nh + receiver) * nr + range]
// * nb + band, ranges at step, 2 step, ... max_range; NaN past max_angle.
// `progress`, if given, is called with a description of the work before it starts.
std::vector<double> build_pe_table(const PeTableRequest& request,
        const std::function<void(const std::string&)>& progress = {});

} // namespace niceops
