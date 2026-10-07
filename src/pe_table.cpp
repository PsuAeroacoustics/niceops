#include "pe_table.hpp"

#include <limits>
#include <sstream>
#include <stdexcept>

#include "pe.hpp"

namespace niceops {

std::vector<double> build_pe_table(const PeTableRequest& q, const std::function<void(const std::string&)>& progress)
{
    const std::size_t nbear = q.c_by_bearing.size(), ns = q.source_heights.size(), nh = q.receiver_heights.size(),
                      nb = q.bands.size();
    if (nbear == 0 || ns == 0 || nh == 0 || nb == 0 || q.admittance.size() != nb)
        throw std::invalid_argument("build_pe_table: bearings, heights, bands and one admittance per band are required");
    if (!(q.step > 0.0) || !(q.max_range >= q.step))
        throw std::invalid_argument("build_pe_table: the step must be positive and no more than the range");
    if (q.realizations < 1)
        throw std::invalid_argument("build_pe_table: realizations must be positive");
    if (q.gamma_t < 0.0 || (q.gamma_t == 0.0 && !q.turbulence_z.empty()))
        throw std::invalid_argument("build_pe_table: a turbulence profile needs a positive reference gamma_t");
    std::size_t nr = 0;
    while ((nr + 1) * q.step <= q.max_range + 1e-9)
        nr++;
    const bool turbulent = q.gamma_t > 0.0;
    // Still air has one realization whatever was asked, so a still table does
    // not depend on --realizations.
    const int realizations = turbulent ? q.realizations : 1;
    std::vector<double> values(nbear * ns * nh * nr * nb, std::numeric_limits<double>::quiet_NaN());
    auto at = [&](std::size_t b, std::size_t s, std::size_t h, std::size_t r, std::size_t band) -> double& {
        return values[(((b * ns + s) * nh + h) * nr + r) * nb + band];
    };
    const bool reciprocal = q.march == PeTableMarch::receiver;
    // A march puts its source at one of `sources` and reads all of `reads`.
    const std::vector<double>& sources = reciprocal ? q.receiver_heights : q.source_heights;
    const std::vector<double>& reads = reciprocal ? q.source_heights : q.receiver_heights;
    // Every march's (band, frequency) tasks go into one pool, highest
    // frequencies first (solve_pe_bands): a march at a time left most cores
    // idle while each march's few top-band tasks finished, about 4.5 of 14
    // busy on the frozen recipe.  The tasks are independent and each band is
    // summed in order, so the values are the same bits either way.
    std::vector<PeSetup> setups;
    std::vector<double> centers;
    setups.reserve(nbear * sources.size() * nb);
    for (std::size_t b = 0, march_index = 0; b < nbear; b++)
        for (std::size_t one = 0; one < sources.size(); one++, march_index++)
            for (std::size_t band = 0; band < nb; band++)
            {
                PeSetup setup;
                setup.source_height = sources[one];
                setup.receiver_heights = reads;
                setup.max_range = q.max_range;
                setup.output_step = q.step;
                setup.admittance = q.admittance[band];
                setup.z = q.z;
                setup.c = q.c_by_bearing[b];
                setup.max_angle = q.max_angle;
                setup.points_per_wavelength = q.points_per_wavelength;
                setup.threads = q.threads;
                if (turbulent)
                {
                    setup.gamma_t = q.gamma_t;
                    setup.outer_scale = q.outer_scale;
                    setup.turbulence_z = q.turbulence_z;
                    setup.turbulence_gamma = q.turbulence_gamma;
                    // Disjoint seeds for every march: solve_pe_bands spends
                    // samples x realizations of them from here.
                    setup.seed = q.seed
                            + static_cast<unsigned>(march_index * static_cast<std::size_t>(q.samples)
                                    * static_cast<std::size_t>(realizations));
                }
                setups.push_back(std::move(setup));
                centers.push_back(q.bands[band]);
            }
    if (progress)
    {
        std::ostringstream text;
        text << nbear * sources.size() << " marches (" << nbear << " bearings x " << sources.size()
             << (reciprocal ? " receiver" : " source") << " heights) x " << nb << " bands x " << q.samples
             << " frequencies" << (realizations > 1 ? " x " + std::to_string(realizations) + " realizations" : "")
             << " in one pool";
        progress(text.str());
    }
    const auto results = solve_pe_bands(setups, centers, q.samples, realizations);
    for (std::size_t b = 0, task = 0; b < nbear; b++)
        for (std::size_t one = 0; one < sources.size(); one++)
            for (std::size_t band = 0; band < nb; band++, task++)
                for (std::size_t k = 0; k < reads.size(); k++)
                    for (std::size_t r = 0; r < nr; r++)
                    {
                        const PeResult& result = results[task];
                        const double level = r < result.range.size() ? result.relative_level[k][r]
                                                                     : std::numeric_limits<double>::quiet_NaN();
                        if (reciprocal)
                            at(b, k, one, r, band) = level;
                        else
                            at(b, one, k, r, band) = level;
                    }
    return values;
}

} // namespace niceops
