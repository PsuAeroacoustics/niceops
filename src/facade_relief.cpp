#include "facade_relief.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>

#include "building_input.hpp"
#include "util/constants.hpp"

namespace niceops {

namespace {

// Where x falls among the ascending `knots`: the lower knot and the weight of
// the upper, held at the ends.
void bracket(const std::vector<double>& knots, double x, std::size_t& lo, double& w)
{
    if (knots.size() == 1 || x <= knots.front())
    {
        lo = 0;
        w = 0.0;
        return;
    }
    if (x >= knots.back())
    {
        lo = knots.size() - 2;
        w = 1.0;
        return;
    }
    lo = static_cast<std::size_t>(std::upper_bound(knots.begin(), knots.end(), x) - knots.begin()) - 1;
    w = (x - knots[lo]) / (knots[lo + 1] - knots[lo]);
}

} // namespace

FacadeRelief FacadeRelief::read(const std::string& path)
{
    const LabeledCsv rows = read_numeric_rows(path, "grating orders");
    auto column = [&](const char* name) {
        const std::size_t c = rows.column({name});
        if (c == LabeledCsv::npos)
            throw std::runtime_error(path + " has no " + name + " column (tools/profile_scattering.py --orders "
                                     "writes frequency_hz, incidence_deg, order, share and period_m)");
        return c;
    };
    const std::size_t fc = column("frequency_hz"), ic = column("incidence_deg"), oc = column("order"),
                      sc = column("share"), pc = column("period_m");
    FacadeRelief relief;
    relief.source_ = path;
    std::map<std::tuple<double, double, int>, double> table;
    std::vector<double> frequencies, degrees;
    for (const auto& r : rows.rows)
    {
        if (!(r[fc] > 0.0) || !(std::abs(r[ic]) < 90.0) || r[oc] != std::round(r[oc]) || !(r[sc] >= 0.0 && r[sc] <= 1.0)
                || !(r[pc] > 0.0))
            throw std::runtime_error(path + ": each frequency_hz and period_m must be positive, incidence_deg within "
                                            "(-90, 90), order a whole number and share within [0, 1]");
        if (relief.period_ == 0.0)
            relief.period_ = r[pc];
        else if (r[pc] != relief.period_)
            throw std::runtime_error(path + " gives more than one period_m");
        const int m = static_cast<int>(r[oc]);
        table[{r[fc], r[ic], m}] += r[sc];
        frequencies.push_back(r[fc]);
        degrees.push_back(r[ic]);
        relief.min_order_ = std::min(relief.min_order_, m);
        relief.max_order_ = std::max(relief.max_order_, m);
    }
    for (auto* v : {&frequencies, &degrees})
    {
        std::sort(v->begin(), v->end());
        v->erase(std::unique(v->begin(), v->end()), v->end());
    }
    relief.frequencies_ = frequencies;
    for (double d : degrees)
        relief.sines_.push_back(std::sin(d * pi / 180.0));
    const std::size_t orders = static_cast<std::size_t>(relief.max_order_ - relief.min_order_ + 1);
    relief.shares_.assign(frequencies.size() * degrees.size() * orders, 0.0);
    for (std::size_t f = 0; f < frequencies.size(); f++)
        for (std::size_t i = 0; i < degrees.size(); i++)
        {
            double sum = 0.0;
            for (int m = relief.min_order_; m <= relief.max_order_; m++)
            {
                const auto found = table.find({frequencies[f], degrees[i], m});
                const double s = found == table.end() ? 0.0 : found->second;
                relief.shares_[(f * degrees.size() + i) * orders + static_cast<std::size_t>(m - relief.min_order_)]
                        = s;
                sum += s;
            }
            if (!(std::abs(sum - 1.0) < 1e-3))
                throw std::runtime_error(path + ": the shares at " + exact_text(frequencies[f]) + " Hz and "
                        + exact_text(degrees[i]) + " degrees sum to " + exact_text(sum) + ", not 1");
        }
    return relief;
}

void FacadeRelief::orders(double sin_psi, double frequency, double sound_speed,
        std::vector<std::pair<int, double>>& out) const
{
    out.clear();
    std::size_t f0, i0;
    double wf, wi;
    bracket(frequencies_, frequency, f0, wf);
    bracket(sines_, std::clamp(sin_psi, -1.0, 1.0), i0, wi);
    const std::size_t f1 = std::min(f0 + 1, frequencies_.size() - 1), i1 = std::min(i0 + 1, sines_.size() - 1);
    const double kick = sound_speed / (frequency * period_);
    double total = 0.0;
    for (int m = min_order_; m <= max_order_; m++)
    {
        // Propagating at this frequency: |sin(theta_m)| < 1.
        if (!(std::abs(-sin_psi + m * kick) < 1.0))
            continue;
        const double s = (1.0 - wf) * ((1.0 - wi) * share(f0, i0, m) + wi * share(f0, i1, m))
                + wf * ((1.0 - wi) * share(f1, i0, m) + wi * share(f1, i1, m));
        if (s > 0.0)
        {
            out.emplace_back(m, s);
            total += s;
        }
    }
    if (total > 0.0)
        for (auto& [m, s] : out)
            s /= total;
}

double FacadeRelief::specular(double sin_psi, double frequency, double sound_speed) const
{
    thread_local std::vector<std::pair<int, double>> shares;
    orders(sin_psi, frequency, sound_speed, shares);
    for (const auto& [m, s] : shares)
        if (m == 0)
            return s;
    return 0.0;
}

double FacadeRelief::random_incidence(double frequency, double sound_speed) const
{
    // The midpoint rule over psi in (-pi/2, pi/2), weights cos(psi).
    constexpr int n = 90;
    double sum = 0.0, weights = 0.0;
    for (int j = 0; j < n; j++)
    {
        const double psi = -0.5 * pi + (j + 0.5) * pi / n;
        sum += std::cos(psi) * (1.0 - specular(std::sin(psi), frequency, sound_speed));
        weights += std::cos(psi);
    }
    return sum / weights;
}

} // namespace niceops
