#pragma once

// A propagation loss table (--propagation_table): what the path between source
// and observer does to each band, computed elsewhere and looked up per pair.
//
// The quantity is the band level relative to a point source's free field at
// the straight-line distance, dB -- everything the path does beyond spherical
// spreading: the ground, refraction, shadows, diffraction.  What the source
// does (the sphere, its directivity, its spectrum) is not in it, and neither,
// unless the table says so, is atmospheric absorption, which the run applies
// as it always does.  The footprint then takes, per band,
//   L = L_sphere - spreading(straight-line range) - absorption + table,
// the same composition as --ground_effect, with the table where the ground
// model and the ray model were.
//
// It assumes what one stratified atmosphere over uniform, flat ground makes
// true: the result depends on where the observer is only through the
// horizontal range and compass bearing from the source, and on the heights of
// the two above the ground.  Five axes, then:
//   bearing          compass bearing the sound travels toward, degrees,
//                    uniform from 0 over 360; one value means the same in
//                    every direction (no wind)
//   source_height    above the ground, m, ascending
//   receiver_height  above the ground, m, ascending
//   range            horizontal, m, ascending
//   band             band center, Hz; the run's bands are found among these
// Values are interpolated linearly along each axis -- in dB, which suits a
// band-averaged level and not a single tone's interference pattern -- and a
// NaN means the method has no value there (the PE's field steeper than its
// angle limit, say).  A pair the table does not cover, by a NaN among the
// values it would interpolate or by falling outside an axis, is propagated by
// the run's own model instead (niceops.cpp counts and logs them).
//
// The file is netCDF, in SI units whatever --units says, and it is not tied
// to the PE: tools/propagation_table.cpp writes it from the PE or from the
// straight-ray ground model, and anything that can fill the array can be read.

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace niceops {

class PropagationTable
{
public:
    static constexpr int format_version = 1;

    struct Axes
    {
        std::vector<double> bearings; // degrees, compass, toward
        std::vector<double> source_heights; // m
        std::vector<double> receiver_heights; // m
        std::vector<double> ranges; // m
        std::vector<double> bands; // Hz
    };

    // How to find the cell of one ascending axis quickly: an axis whose knots
    // are evenly spaced is indexed by arithmetic and corrected against its
    // neighbors, so the cell is exactly the one a binary search would give.
    struct AxisLocator
    {
        bool uniform = false;
        double inverse_step = 0.0;
    };

    PropagationTable() = default;
    // `values` is [bearing][source][receiver][range][band]; it throws
    // std::invalid_argument if the axes are not as the header says.
    PropagationTable(Axes axes, std::vector<double> values, std::string method, bool includes_absorption,
            std::map<std::string, std::string> provenance = {});

    static PropagationTable load(const std::string& path);
    void save(const std::string& path) const;

    const Axes& axes() const
    {
        return axes_;
    }
    const std::string& method() const
    {
        return method_;
    }
    bool includes_absorption() const
    {
        return includes_absorption_;
    }
    const std::map<std::string, std::string>& provenance() const
    {
        return provenance_;
    }
    bool isotropic() const
    {
        return axes_.bearings.size() == 1;
    }
    double at(std::size_t bearing, std::size_t source, std::size_t receiver, std::size_t range,
            std::size_t band) const
    {
        return values_[index(bearing, source, receiver, range) * axes_.bands.size() + band];
    }

    // Where a run's bands are in the table: for each run band in [first,
    // last), its column.  Per run rather than in the table, because one run
    // can read several databases and the table is loaded once.
    struct BandMap
    {
        std::size_t first = 0;
        std::size_t last = 0;
        std::vector<std::size_t> column;
    };

    // Find the run's bands in the table, to 0.5%: every one in [first_band,
    // last_band) has to be there.  Throws std::runtime_error naming the
    // missing ones.
    BandMap match_bands(const std::vector<double>& run_bands, std::size_t first_band, std::size_t last_band) const;

    // Whether a receiver at this height lies within the table's.
    bool covers_receiver(double height) const;

    // The table's level for each of the run's bands in [first_band,
    // last_band), written to out[band]; false, writing nothing, if the table
    // does not cover the point.
    bool lookup(double bearing_degrees, double source_height, double receiver_height, double range,
            const BandMap& bands, double* out) const;

    // The largest value in each of the run's bands, as a power ratio: what a
    // bound on the level has to allow the table to add.  1 for a band the
    // table never raises.
    std::vector<double> max_power_ratio(const BandMap& bands) const;

private:
    std::size_t index(std::size_t bearing, std::size_t source, std::size_t receiver, std::size_t range) const
    {
        return ((bearing * axes_.source_heights.size() + source) * axes_.receiver_heights.size() + receiver)
                * axes_.ranges.size()
                + range;
    }

    Axes axes_;
    AxisLocator source_locator_, receiver_locator_, range_locator_;
    std::vector<double> values_;
    std::string method_;
    bool includes_absorption_ = false;
    std::map<std::string, std::string> provenance_;
};

// The turbulence and ground roughness a table's levels already carry, as its
// builder recorded them in its provenance.  A table pair takes them from the
// table, not from the run's --turbulence, --gamma_t or --roughness, which
// reach only the pairs the table does not cover; the run compares the two.
//   model "none"           still air (and every table that records nothing)
//   model "chessell"       the ground model's coherence, value = C, rad s/sqrt(m)
//   model "harmonoise"     the ground model's coherence, value = gamma_t, m^-2/3
//   model "phase_screens"  a turbulent PE table, value = gamma_t, m^-2/3
struct TableCoherence
{
    std::string model = "none";
    double value = 0.0;
    double roughness = 0.0; // m, the ground model's Kirchhoff roughness
    // A turbulent PE table built before the record wrote gamma_t to six
    // decimals ("0.000002"): its strength is known only that far, and that
    // text is kept to compare against.
    std::string rounded_gamma_t;
};

TableCoherence table_coherence(const std::map<std::string, std::string>& provenance);

// Whether a run with `run` would treat its covered and uncovered pairs
// alike: the same model and strength (a HARMONOISE run matches a PE table's
// phase screens of the same gamma_t) and the same roughness, either way round.
// Returns the warning's text when they differ, empty when they agree.
std::string coherence_mismatch(const TableCoherence& run, const TableCoherence& table);

// Writes a ground-model table's coherence into its provenance (keys
// turbulence_model, turbulence or gamma_t, roughness_m) at full precision.
void record_coherence(std::map<std::string, std::string>& provenance, const TableCoherence& coherence);

} // namespace niceops
