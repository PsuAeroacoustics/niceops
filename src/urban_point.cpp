#include "urban_point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <unordered_map>

#include "building_input.hpp"
#include "building_paths.hpp"
#include "buildings.hpp"
#include "point_paths.hpp"
#include "source_coherence.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"
#include "weighting.hpp"

namespace niceops {

namespace {

// The regularized lower incomplete gamma P(k, x): its series below k + 1, its
// continued fraction above (Lentz).
double gamma_p(double k, double x)
{
    if (x <= 0.0)
        return 0.0;
    const double front = std::exp(-x + k * std::log(x) - std::lgamma(k));
    if (x < k + 1.0)
    {
        double term = 1.0 / k, sum = term;
        for (int n = 1; n < 1000 && std::abs(term) > 1e-16 * std::abs(sum); n++)
        {
            term *= x / (k + n);
            sum += term;
        }
        return std::min(1.0, sum * front);
    }
    double b = x + 1.0 - k, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
    for (int n = 1; n < 1000; n++)
    {
        const double an = -n * (n - k);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < 1e-300)
            d = 1e-300;
        c = b + an / c;
        if (std::abs(c) < 1e-300)
            c = 1e-300;
        d = 1.0 / d;
        const double step = d * c;
        h *= step;
        if (std::abs(step - 1.0) < 1e-15)
            break;
    }
    return std::max(0.0, 1.0 - front * h);
}

// The q-quantile of the gamma distribution with mean `mean` and variance
// `variance`, by bisection.
double gamma_quantile(double mean, double variance, double q)
{
    if (!(variance > 0.0))
        return mean;
    const double k = mean * mean / variance, theta = variance / mean;
    double lo = 0.0, hi = mean + 60.0 * std::sqrt(variance);
    for (int i = 0; i < 200; i++)
    {
        const double mid = 0.5 * (lo + hi);
        (gamma_p(k, mid / theta) < q ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// partially_coherent_energies's buffers.
struct PairScratch
{
    std::vector<Vec3> launch;
    std::vector<double> azimuth;
    std::vector<std::uint32_t> first, second;
    std::vector<double> degrees, power, spread_power;
    std::vector<std::int8_t> wrap;
    std::vector<double> mean, root_spread;
    std::vector<std::uint32_t> live;
    std::vector<std::int64_t> slot;
    std::vector<std::complex<double>> cross;
    std::vector<std::size_t> order, parent, group_size;
    std::unordered_map<std::uint64_t, std::int64_t> slots;
    std::vector<double> length, re, im, sine, cosine;
};

} // namespace

void partially_coherent_energies(const PointSourceField& field, const std::vector<double>& reflection,
        double floor_reflection, const std::vector<double>& bands, std::size_t first_band, double sound_speed,
        const std::vector<double>& air, const double* weights, const SourceCoherence& coherence,
        const SphereFrame& frame, double* energy, double* variance)
{
    const auto& rays = field.rays();
    const std::size_t n = rays.size();
    const bool measured = coherence.sum == RaySum::measured;
    // The cone's mean coherency is the same in every band.
    const bool cone = coherence.sum == RaySum::cone;
    const bool rotor = measured && coherence.rotation_rate > 0.0;
    // The rotor's delay as a length per radian of azimuth about the shaft.
    // A ray's own delay is that times its azimuth (counterclockwise seen
    // from above is positive about the normal); a pair's takes the
    // difference wrapped to +-pi, which is the rays' own difference shifted
    // by a whole turn where it is not.  So a ray carries its effective length
    // L + per_radian psi, and a pair where the wrap shifts it adds or takes
    // away one turn's length: the pairs need no trigonometry of their own.
    const double per_radian = rotor ? sound_speed * coherence.rotation_sense / coherence.rotation_rate : 0.0;
    const double turn = 2.0 * pi * per_radian;
    // The buffers persist per thread: a field has hundreds of rays, so tens
    // of thousands of pairs, and a run sums tens of thousands of fields.
    // Bound to references once: each use of a thread_local itself is a
    // lookup on macOS, which in the pair loops cost more than the buffers
    // saved.
    thread_local PairScratch scratch;
    auto& launch = scratch.launch;
    auto& azimuth = scratch.azimuth;
    launch.resize(n);
    azimuth.assign(n, 0.0);
    const Vec3 left = frame.normal.cross(frame.nose);
    for (std::size_t i = 0; i < n; i++)
    {
        launch[i] = point_ray_launch(rays[i]);
        if (rotor)
            azimuth[i] = std::atan2(left.dot(launch[i]), frame.nose.dot(launch[i]));
    }
    // Every pair: its rays, the distance between their launches as an angle
    // (the chord in degrees: within a few hundredths of a degree of the angle
    // while the decay has not run its course, and with it the decay is a
    // positive-definite kernel on the sphere), that to the decay's exponent,
    // and the turns its wrap adds (-1, 0 or +1).
    auto& first = scratch.first;
    auto& second = scratch.second;
    auto& degrees = scratch.degrees;
    auto& power = scratch.power;
    auto& spread_power = scratch.spread_power;
    auto& wrap = scratch.wrap;
    for (auto* v : {&degrees, &power, &spread_power})
        v->clear();
    first.clear();
    second.clear();
    wrap.clear();
    if (coherence.sum != RaySum::energy)
    {
        const std::size_t count = n * (n - 1) / 2;
        first.reserve(count);
        second.reserve(count);
        degrees.reserve(count);
        power.reserve(count);
        if (variance != nullptr)
            spread_power.reserve(count);
        wrap.reserve(count);
        for (std::size_t i = 0; i < n; i++)
            for (std::size_t j = i + 1; j < n; j++)
            {
                first.push_back(static_cast<std::uint32_t>(i));
                second.push_back(static_cast<std::uint32_t>(j));
                const double d = (launch[i] - launch[j]).norm() * 180.0 / pi;
                degrees.push_back(d);
                power.push_back(measured ? std::pow(d, coherence.exponent) : 0.0);
                if (variance != nullptr)
                    spread_power.push_back(coherence.spread_exponent == 1.0 ? d : std::pow(d, coherence.spread_exponent));
                const double raw = azimuth[i] - azimuth[j];
                wrap.push_back(static_cast<std::int8_t>(raw > pi ? -1 : raw < -pi ? 1 : 0));
            }
    }
    const std::size_t pairs = first.size();
    const double half_band = std::pow(2.0, 1.0 / 6.0);
    // A pair whose mean coherency is below this adds nothing worth its cost;
    // the decay only grows with the band, so it never comes back.
    constexpr double negligible = 14.0; // exp(-14) = 8e-7
    // Per band: each pair's mean coherency (at the band's centre, as it was
    // measured; the cone's is every band's) and the pairs it leaves anything
    // of; with the variance, each
    // pair's root spread and its group pair's slot (none within a group).
    auto& mean = scratch.mean;
    auto& root_spread = scratch.root_spread;
    auto& live = scratch.live;
    auto& slot = scratch.slot;
    auto& cross = scratch.cross;
    auto& order = scratch.order;
    auto& parent = scratch.parent;
    auto& group_size = scratch.group_size;
    auto& slots = scratch.slots;
    mean.assign(pairs, 0.0);
    live.resize(pairs);
    std::iota(live.begin(), live.end(), 0);
    parent.resize(n);
    // Per sub-band: each ray's phasor at the middle wavenumber over its
    // effective length, and the sine and cosine of the half width times that
    // length (band_averaged_energy's way of forming each pair's sinc).
    auto& length = scratch.length;
    auto& re = scratch.re;
    auto& im = scratch.im;
    auto& sine = scratch.sine;
    auto& cosine = scratch.cosine;
    for (auto* v : {&length, &re, &im, &sine, &cosine})
        v->resize(n);
    field.sub_band_problems(
            reflection, floor_reflection, bands, sound_speed,
            [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                    const std::vector<double>& sorted_length, const SubBand& s) {
                const double centre = bands[band];
                const bool rotates = rotor && coherence.rotates(centre);
                if (part == 0)
                {
                    const double decay = measured ? std::pow(centre / coherence.scale, coherence.exponent) : 0.0;
                    if (cone)
                    {
                        if (band == 0)
                        {
                            std::size_t kept = 0;
                            for (std::size_t q = 0; q < pairs; q++)
                            {
                                mean[q] = coherence.mean(centre, degrees[q]);
                                if (mean[q] > 0.0)
                                    live[kept++] = static_cast<std::uint32_t>(q);
                            }
                            live.resize(kept);
                        }
                    }
                    else if (variance == nullptr)
                    {
                        std::size_t kept = 0;
                        for (std::uint32_t q : live)
                        {
                            const double e = decay * power[q];
                            if (e <= negligible)
                            {
                                mean[q] = std::exp(-e);
                                live[kept++] = q;
                            }
                        }
                        live.resize(kept);
                    }
                    else
                        for (std::size_t q = 0; q < pairs; q++)
                        {
                            const double e = decay * power[q];
                            mean[q] = e <= negligible ? std::exp(-e) : 0.0;
                        }
                    if (variance != nullptr)
                    {
                        // Rays whose mean coherency is 0.99 or more, by single
                        // linkage, deviate together.
                        std::iota(parent.begin(), parent.end(), 0);
                        auto root = [&](std::size_t i) {
                            while (parent[i] != i)
                                i = parent[i] = parent[parent[i]];
                            return i;
                        };
                        for (std::size_t q = 0; q < pairs; q++)
                            if (mean[q] >= 0.99)
                                parent[root(first[q])] = root(second[q]);
                        group_size.assign(n, 0);
                        for (std::size_t i = 0; i < n; i++)
                            group_size[root(i)]++;
                        // The spread, max(gamma^2 - mean^2, 0), with gamma^2
                        // = t + (1 - t) exp(-(f theta / spread_scale)^p).
                        const double t = coherence.tonal(first_band + band, centre);
                        const double spread_decay
                                = std::pow(centre / coherence.spread_scale, coherence.spread_exponent);
                        root_spread.resize(pairs);
                        slot.assign(pairs, -1);
                        slots.clear();
                        // A pair of rays each alone in its group is its own
                        // group pair: its slot follows the shared ones'.
                        std::int64_t own = 0;
                        for (std::size_t q = 0; q < pairs; q++)
                        {
                            const std::size_t a = root(first[q]), b = root(second[q]);
                            if (a == b)
                                continue;
                            const double gamma2 = t + (1.0 - t) * std::exp(-spread_decay * spread_power[q]);
                            root_spread[q] = std::sqrt(std::max(gamma2 - mean[q] * mean[q], 0.0));
                            if (root_spread[q] == 0.0)
                                continue;
                            if (group_size[a] == 1 && group_size[b] == 1)
                            {
                                slot[q] = -2 - 2 * own++;
                                continue;
                            }
                            // Oriented by the pair's first ray's group, so a
                            // reversed pair adds its conjugate.
                            const std::uint64_t key
                                    = (static_cast<std::uint64_t>(std::min(a, b)) << 32) | std::max(a, b);
                            const auto at = slots.emplace(key, static_cast<std::int64_t>(slots.size())).first->second;
                            slot[q] = a < b ? 2 * at : 2 * at + 1;
                        }
                        const auto shared = static_cast<std::int64_t>(slots.size());
                        for (std::int64_t& v : slot)
                            if (v <= -2)
                                v = 2 * (shared + (-2 - v) / 2);
                        cross.assign(static_cast<std::size_t>(shared + own), 0.0);
                    }
                }
                const double mid = 0.5 * (s.k1 + s.k2), half = 0.5 * (s.k2 - s.k1);
                for (std::size_t r = 0; r < n; r++)
                {
                    const std::size_t i = order[r];
                    length[i] = sorted_length[r] + (rotates ? per_radian * azimuth[i] : 0.0);
                    const std::complex<double> phasor = amplitude[r] * std::polar(1.0, mid * length[i]);
                    re[i] = phasor.real();
                    im[i] = phasor.imag();
                    sine[i] = std::sin(half * length[i]);
                    cosine[i] = std::cos(half * length[i]);
                }
                // One turn's shift: its phase, and its sine and cosine at the
                // half width.
                const std::complex<double> turn_phase = rotates ? std::polar(1.0, mid * turn) : 1.0;
                const double turn_sine = rotates ? std::sin(half * turn) : 0.0;
                const double turn_cosine = rotates ? std::cos(half * turn) : 1.0;
                // A pair's cross term, 2 P_i conj(P_j) times its sinc; the
                // real part alone unless `complex`.
                auto term = [&](std::size_t q, bool complex) -> std::complex<double> {
                    const std::size_t i = first[q], j = second[q];
                    const int w = rotates ? wrap[q] : 0;
                    const double x = half * (length[i] - length[j] + w * turn);
                    if (std::abs(x) > 30.0)
                        return 0.0;
                    double sinc;
                    if (std::abs(x) < 0.1)
                    {
                        const double x2 = x * x;
                        sinc = 1.0 + x2 * (-1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 * (-1.0 / 5040.0 + x2 / 362880.0)));
                    }
                    else
                    {
                        // sin(a + w b) from the rays' sines and cosines and the turn's.
                        const double sin_d = sine[i] * cosine[j] - cosine[i] * sine[j];
                        if (w == 0)
                            sinc = sin_d / x;
                        else
                        {
                            const double cos_d = cosine[i] * cosine[j] + sine[i] * sine[j];
                            sinc = (sin_d * turn_cosine + w * cos_d * turn_sine) / x;
                        }
                    }
                    std::complex<double> p(re[i] * re[j] + im[i] * im[j], im[i] * re[j] - re[i] * im[j]);
                    if (w != 0)
                        p *= w > 0 ? turn_phase : std::conj(turn_phase);
                    return 2.0 * sinc * (complex ? p : std::complex<double>(p.real(), 0.0));
                };
                double sum = 0.0;
                for (std::size_t i = 0; i < n; i++)
                    sum += re[i] * re[i] + im[i] * im[i];
                if (variance == nullptr)
                    for (std::uint32_t q : live)
                        sum += mean[q] * term(q, false).real();
                else
                {
                    const double weight = (s.g2 - s.g1) / (centre * half_band - centre / half_band);
                    for (std::size_t q = 0; q < pairs; q++)
                    {
                        if (mean[q] == 0.0 && slot[q] < 0)
                            continue;
                        const std::complex<double> t = term(q, slot[q] >= 0);
                        sum += mean[q] * t.real();
                        if (slot[q] >= 0)
                        {
                            const std::complex<double> c = weight * root_spread[q] * t;
                            cross[static_cast<std::size_t>(slot[q] / 2)] += slot[q] % 2 == 0 ? c : std::conj(c);
                        }
                    }
                }
                energy[3 * band + static_cast<std::size_t>(part)] = std::max(sum, 0.0);
                if (variance != nullptr && part == 2)
                {
                    double v = 0.0;
                    for (const std::complex<double>& c : cross)
                        v += 0.5 * std::norm(c);
                    variance[band] = v;
                }
            },
            air, weights, &order);
}

struct UrbanDistrict::Impl
{
    BuildingInput input;
    Buildings city;
    // Each facade's energy reflection coefficient, and the floor's.
    std::vector<double> reflection;
    double floor_reflection = 1.0;
    PointPathOptions paths;
};

UrbanDistrict::UrbanDistrict() : impl_(std::make_unique<Impl>()) {}

UrbanDistrict::~UrbanDistrict() = default;

std::shared_ptr<const UrbanDistrict> UrbanDistrict::read(const std::string& file, const Input& input)
{
    BuildingInputOptions options;
    options.units = input.units;
    options.to_meters = input.to_meters;
    options.origin = input.origin;
    options.frame_bearing = input.frame_bearing;
    options.absorption = input.wall_absorption;
    std::shared_ptr<UrbanDistrict> district(new UrbanDistrict());
    Impl& d = *district->impl_;
    d.input = read_buildings(file, options);
    try
    {
        d.city = Buildings(d.input.footprints, 1.0 * pi / 180.0);
    }
    catch (const std::invalid_argument& e)
    {
        throw std::runtime_error(std::string(e.what()) + " (counting buildings from 1 in file order)");
    }
    d.reflection.resize(d.city.facades().size());
    for (std::size_t f = 0; f < d.reflection.size(); f++)
        d.reflection[f] = 1.0 - d.input.absorption[d.city.facade_building(f)];
    d.floor_reflection = 1.0 - input.floor_absorption;
    return district;
}

std::size_t UrbanDistrict::building_count() const
{
    return impl_->city.size();
}

std::size_t UrbanDistrict::facade_count() const
{
    return impl_->city.facades().size();
}

bool UrbanDistrict::scatters() const
{
    const auto& s = impl_->input.scattering;
    return std::any_of(s.begin(), s.end(), [](double v) { return v > 0.0; });
}

std::string UrbanDistrict::misplaced(const Vec3& p) const
{
    if (!(p.z >= 0.0))
        return "below the floor";
    const Buildings& city = impl_->city;
    const std::size_t b = city.inside(p.x, p.y);
    if (b < city.size() && p.z < city.footprints()[b].top)
        return "inside building " + impl_->input.ids[b];
    return "";
}

void UrbanDistrict::levels(const Database& database, const std::vector<Database::SourceState>& states,
        const std::vector<Vec3>& microphones, const double* absorption, const std::vector<double>& gains,
        std::size_t first_band, std::size_t last_band, double sound_speed, unsigned threads,
        const SourceCoherence& coherence, std::vector<double>& levels, std::vector<double>* band_levels,
        const std::vector<const Database::PackedSpectra*>* packed, PhaseInterval* interval) const
{
    const Impl& d = *impl_;
    const Buildings& city = d.city;
    const std::vector<double>& centers = database.band_centers();
    const std::size_t all = centers.size();
    last_band = std::min(last_band, all);
    first_band = std::min(first_band, last_band);
    const std::size_t nb = last_band - first_band, nm = microphones.size(), ns = states.size();
    constexpr double silence = -std::numeric_limits<double>::infinity();
    levels.assign(ns * nm, silence);
    if (band_levels != nullptr)
        band_levels->assign(ns * nm * all, silence);
    if (interval != nullptr)
    {
        interval->lower.assign(ns * nm, silence);
        interval->upper.assign(ns * nm, silence);
    }
    if (nb == 0 || nm == 0 || ns == 0)
        return;

    // The metric's bands, their facades' reflections, and the air's
    // attenuation of amplitude per metre in each of their sub-bands.
    const std::vector<double> bands(centers.begin() + static_cast<std::ptrdiff_t>(first_band),
            centers.begin() + static_cast<std::ptrdiff_t>(last_band));
    std::vector<double> reflection(d.reflection.size() * nb);
    for (std::size_t f = 0; f < d.reflection.size(); f++)
        std::fill_n(reflection.begin() + static_cast<std::ptrdiff_t>(f * nb), nb, d.reflection[f]);
    std::vector<double> nepers(nb), sub_band_air(3 * nb);
    for (std::size_t b = 0; b < nb; b++)
    {
        nepers[b] = absorption[first_band + b] * std::log(10.0) / 20.0;
        for (std::size_t part = 0; part < 3; part++)
            sub_band_air[3 * b + part] = nepers[b];
    }
    const double radius = database.sphere_radius();
    const double at_one_metre = 20.0 * std::log10(radius);
    // The interval's percentiles, as fractions.
    const double low_q = interval != nullptr ? 0.5 - interval->percent / 200.0 : 0.0, high_q = 1.0 - low_q;

    std::vector<Vec3> sources(ns);
    for (std::size_t k = 0; k < ns; k++)
        sources[k] = states[k].position;
    const std::vector<std::size_t> near_sources = edges_near(city, sources, d.paths.radius),
                                   near_microphones = edges_near(city, microphones, d.paths.radius);
    const unsigned workers = std::max(1u, threads);
    const std::size_t block = std::max<std::size_t>(64, 4 * workers);
    std::vector<std::unique_ptr<PathEnd>> source_ends;
    for (std::size_t first = 0; first < ns; first += block)
    {
        const std::size_t nk = std::min(block, ns - first);
        source_ends.clear();
        source_ends.resize(nk);
        parallel_chunks(nk, workers, [&](std::size_t begin, std::size_t end) {
            for (std::size_t k = begin; k < end; k++)
                source_ends[k] = std::make_unique<PathEnd>(city, sources[first + k], true, d.paths, near_microphones);
        });
        // Microphone by microphone, so a worker keeps the microphone's end
        // while its pairs last.
        parallel_chunks(nm * nk, workers, [&](std::size_t begin, std::size_t end) {
            std::unique_ptr<PathEnd> microphone;
            std::size_t microphone_row = nm;
            std::vector<double> weights, spectrum(all), out(nb), full(all), energy(3 * nb), variance(nb);
            for (std::size_t index = begin; index < end; index++)
            {
                const std::size_t m = index / nk, k = index % nk, state = first + k;
                if (m != microphone_row)
                {
                    microphone = std::make_unique<PathEnd>(city, microphones[m], false, d.paths, near_sources);
                    microphone_row = m;
                }
                const PointSourceField field(city, *source_ends[k], *microphone, d.paths);
                const auto& rays = field.rays();
                if (rays.empty())
                    continue;
                // Each ray's amplitude at 1 m along its launch, from the
                // sphere at its radius; no spreading gain inside the sphere,
                // and no absorption there.
                weights.assign(rays.size() * nb, 0.0);
                for (std::size_t r = 0; r < rays.size(); r++)
                {
                    database.spectrum_along(states[state], point_ray_launch(rays[r]), spectrum.data(), first_band,
                            last_band, packed != nullptr ? (*packed)[state] : nullptr);
                    const double inside = std::min(rays[r].length, radius);
                    for (std::size_t b = 0; b < nb; b++)
                    {
                        const double level = spectrum[first_band + b];
                        if (std::isfinite(level))
                            weights[r * nb + b] = std::pow(10.0, (level + at_one_metre) / 20.0) * (inside / radius)
                                    * std::exp(nepers[b] * inside);
                    }
                }
                if (coherence.sum == RaySum::coherent)
                    field.band_levels(reflection, d.floor_reflection, bands, sound_speed, out.data(), sub_band_air,
                            weights.data());
                else
                {
                    partially_coherent_energies(field, reflection, d.floor_reflection, bands, first_band, sound_speed,
                            sub_band_air, weights.data(), coherence, states[state].frame, energy.data(),
                            interval != nullptr ? variance.data() : nullptr);
                    for (std::size_t b = 0; b < nb; b++)
                        out[b] = band_level(bands[b], sub_bands(bands[b], sound_speed), &energy[3 * b]);
                }
                std::fill(full.begin(), full.end(), silence);
                std::copy(out.begin(), out.end(), full.begin() + static_cast<std::ptrdiff_t>(first_band));
                levels[state * nm + m] = weighted_level(full.data(), gains);
                if (interval != nullptr)
                {
                    // The A-weighted energy's mean and variance (the bands'
                    // phases independent) and the gamma distribution with
                    // them.
                    double mean = 0.0, spread = 0.0;
                    for (std::size_t b = 0; b < nb; b++)
                    {
                        const double g = gains[first_band + b];
                        if (std::isfinite(out[b]))
                            mean += g * std::pow(10.0, 0.1 * out[b]);
                        spread += g * g * variance[b];
                    }
                    if (mean > 0.0)
                    {
                        interval->lower[state * nm + m] = 10.0 * std::log10(gamma_quantile(mean, spread, low_q));
                        interval->upper[state * nm + m] = 10.0 * std::log10(gamma_quantile(mean, spread, high_q));
                    }
                }
                if (band_levels != nullptr)
                    std::copy(full.begin(), full.end(),
                            band_levels->begin() + static_cast<std::ptrdiff_t>((state * nm + m) * all));
            }
        });
    }
}

} // namespace niceops
