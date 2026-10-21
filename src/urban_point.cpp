#include "urban_point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <unordered_map>

#include "building_input.hpp"
#include "building_paths.hpp"
#include "diffuse_facades.hpp"
#include "facade_relief.hpp"
#include "grating_paths.hpp"
#include "ground_effect.hpp"
#include "plate_receiver.hpp"
#include "gpu/band_stream.hpp"
#include "gpu/kernels/band_energy_types.h"
#include "gpu/point_rays.hpp"
#include "buildings.hpp"
#include "point_paths.hpp"
#include "source_coherence.hpp"
#include "turbulence_scattering.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"
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


// Where an increasing cdf reaches q in [lo, hi], by bisection to round-off.
template <typename Cdf>
double invert(Cdf cdf, double lo, double hi, double q)
{
    for (int i = 0; i < 200; i++)
    {
        const double mid = 0.5 * (lo + hi);
        if (!(mid > lo && mid < hi))
            break;
        (cdf(mid) < q ? lo : hi) = mid;
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

double gamma_quantile(double mean, double variance, double q)
{
    if (!(variance > 0.0))
        return mean;
    const double k = mean * mean / variance, theta = variance / mean;
    // Past a shape of 1e4 gamma_p's series needs more than its 1000 terms
    // near the mean, and the Wilson-Hilferty cube is the quantile to 4e-7 dB
    // at the 0.5th percentile and closer nearer the middle; z is the normal
    // quantile, by bisection on erfc.
    if (k > 1e4)
    {
        const double z = invert([](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }, -40.0, 40.0, q);
        const double cube = 1.0 - 1.0 / (9.0 * k) + z / (3.0 * std::sqrt(k));
        return mean * cube * cube * cube;
    }
    return invert([&](double x) { return gamma_p(k, x / theta); }, 0.0, mean + 60.0 * std::sqrt(variance), q);
}

void partially_coherent_energies(const PointSourceField& field, const std::vector<double>& reflection,
        const FloorModel& floor, const std::vector<double>& bands, std::size_t first_band, double sound_speed,
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
            reflection, floor, bands, sound_speed,
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
    std::vector<std::pair<double, double>> scattering_spectrum;
    std::string scattering_file;
    double patch_size = 4.0, diffuse_radius = 0.0;
    double turbulence = 0.0;
    // The reliefs, each facade's (null for a plain one), how they scatter,
    // and the reflections a grating ray may take besides its order.
    std::vector<FacadeRelief> reliefs;
    std::vector<const FacadeRelief*> facade_reliefs;
    // Whether the other orders go as rays, and whether diffusely.
    bool relief_orders = true;
    bool relief_lambert = true;
    std::size_t grating_reflections = 2;
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
    options.scattering = input.scattering;
    options.roofs = input.roofs;
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
    d.paths.sample_shared_faces = !input.exact_shared_faces;
    d.paths.edge_source = input.edge_source;
    if (!input.scattering_file.empty())
        d.scattering_spectrum = read_scattering_spectrum(input.scattering_file, input.units);
    d.scattering_file = input.scattering_file;
    d.patch_size = input.patch_size;
    d.diffuse_radius = input.diffuse_radius > 0.0 ? input.diffuse_radius : d.paths.radius;
    d.turbulence = input.turbulence;
    for (const std::string& relief : input.relief_files)
        d.reliefs.push_back(FacadeRelief::read(relief));
    d.relief_orders = input.relief_scattering != "lambert";
    d.relief_lambert = input.relief_scattering != "orders";
    d.grating_reflections = input.grating_reflections;
    if (!d.reliefs.empty())
    {
        d.facade_reliefs.assign(d.city.facades().size(), nullptr);
        for (std::size_t f = 0; f < d.facade_reliefs.size(); f++)
        {
            const int which = relief_of(d.input, d.city.facade_building(f), d.reliefs.size());
            if (which > 0)
                d.facade_reliefs[f] = &d.reliefs[static_cast<std::size_t>(which - 1)];
        }
        d.paths.reliefs = &d.facade_reliefs;
        // Grating rays' source side reaches as deep as their own orders.
        if (d.relief_orders)
            d.paths.grating_reflections = d.grating_reflections;
    }
    return district;
}

bool UrbanDistrict::edge_source() const
{
    return impl_->paths.edge_source;
}

std::size_t UrbanDistrict::relief_walls() const
{
    return static_cast<std::size_t>(std::count_if(impl_->facade_reliefs.begin(), impl_->facade_reliefs.end(),
            [](const FacadeRelief* r) { return r != nullptr; }));
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

std::vector<double> UrbanDistrict::band_scattering(const std::vector<double>& bands) const
{
    if (impl_->scattering_spectrum.empty())
        return std::vector<double>(bands.size(), 1.0);
    return scattering_at(impl_->scattering_spectrum, bands, impl_->scattering_file);
}

double UrbanDistrict::patch_size() const
{
    return impl_->patch_size;
}

double UrbanDistrict::diffuse_radius() const
{
    return impl_->diffuse_radius;
}

std::string UrbanDistrict::misplaced(const Vec3& p) const
{
    return misplaced_point(impl_->city, impl_->input.ids, p);
}

bool UrbanDistrict::tilted(const Vec3& p) const
{
    const Buildings::Support::Kind kind = impl_->city.support(p).kind;
    return kind == Buildings::Support::Kind::facade || kind == Buildings::Support::Kind::slope;
}

double UrbanDistrict::floor_height(const Vec3& ground) const
{
    return std::abs(ground.z) <= Buildings::contact || impl_->city.support(ground).kind != Buildings::Support::Kind::none
            ? 0.0
            : ground.z;
}

bool UrbanDistrict::gpu_supports(const SourceCoherence& coherence, bool interval)
{
    return !interval && coherence.sum != RaySum::measured;
}

std::size_t UrbanDistrict::levels(const Database& database, const std::vector<Database::SourceState>& states,
        const std::vector<Vec3>& microphones, const double* absorption, const std::vector<double>& gains,
        std::size_t first_band, std::size_t last_band, double sound_speed, unsigned threads,
        const SourceCoherence& coherence, std::vector<double>& levels, std::vector<double>* band_levels,
        const std::vector<const Database::PackedSpectra*>* packed, PhaseInterval* interval,
        double state_tolerance, gpu::Device* device, const std::vector<MicrophoneFloor>* floors) const
{
    if (floors != nullptr && floors->size() != microphones.size())
        throw std::invalid_argument("UrbanDistrict::levels: the floors are not one per microphone");
    if (device != nullptr && !gpu_supports(coherence, interval != nullptr))
        throw std::invalid_argument("UrbanDistrict::levels: the GPU does not take the measured coherence or an interval");
    if (device != nullptr && impl_->paths.edge_source)
        throw std::invalid_argument("UrbanDistrict::levels: the edge-source integral is the CPU's");
    const Impl& d = *impl_;
    const Buildings& city = d.city;
    const std::vector<double>& centers = database.band_centers();
    const std::size_t all = centers.size();
    last_band = std::min(last_band, all);
    // The paths' options for this run's bands: with the edge source, a ray
    // whose Keller point lies past an end is kept to 100 radians of path
    // excess at the lowest band, as niceops_point_field keeps it.
    PointPathOptions paths = d.paths;
    if (paths.edge_source && first_band < last_band)
        paths.edge_source_reach = 100.0 / (2.0 * pi * centers[first_band] / sound_speed);
    // UTD's end waves fade out by end_fade_hi radians, at the lowest band last.
    if (first_band < last_band)
        paths.end_reach = end_fade_hi / (2.0 * pi * centers[first_band] / sound_speed);
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
        return 0;

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
    // Each microphone's floor: the district's energy coefficient, or a
    // ground's admittance at the sub-bands' middles (one per distinct
    // ground), with a flush microphone's plate; and the energy it reflects of
    // a diffuse field in each band, for the fields that carry energy alone.
    std::vector<FloorModel> floor_of(nm, FloorModel(d.floor_reflection));
    std::vector<std::vector<double>> floor_energy(nm, std::vector<double>(nb, d.floor_reflection));
    // The sub-bands' middle frequencies, and each one's database band.
    const std::vector<double> middles = sub_band_middles(bands, sound_speed);
    std::vector<std::size_t> problem_band;
    for (std::size_t b = 0; b < nb; b++)
        problem_band.insert(problem_band.end(), 3, first_band + b);
    if (floors != nullptr)
    {
        // Per ground, its admittance in each sub-band and its diffuse
        // reflection in each band (at the band's middle sub-band).
        struct Ground
        {
            std::vector<std::complex<double>> admittance;
            std::vector<double> diffuse;
        };
        std::map<const GroundModel*, Ground> grounds;
        for (std::size_t m = 0; m < nm; m++)
        {
            const MicrophoneFloor& f = (*floors)[m];
            if (f.ground == nullptr)
                continue;
            auto found = grounds.find(f.ground);
            if (found == grounds.end())
            {
                Ground g {ground_impedance(*f.ground, middles, sound_speed), std::vector<double>(nb)};
                for (std::size_t b = 0; b < nb; b++)
                    g.diffuse[b] = diffuse_reflection(g.admittance[3 * b + 1]);
                found = grounds.emplace(f.ground, std::move(g)).first;
            }
            FloorModel& floor = floor_of[m];
            floor.admittance = found->second.admittance;
            floor.boundary_loss = f.ground->boundary_loss_correction;
            floor.roughness = f.ground->roughness;
            floor_energy[m] = found->second.diffuse;
            if (f.plate != nullptr)
            {
                const PlateReceiver* plate = f.plate;
                const double offset = f.plate_azimuth;
                floor.plate = [plate, offset, &middles, &problem_band](double elevation, double azimuth,
                                      std::complex<double>* direct, std::complex<double>* image) {
                    // The plate's frame turns the microphone's offset onto +y
                    // (plate_frame_azimuth).
                    plate->transfer(elevation, azimuth - offset + 0.5 * pi, problem_band.data(), middles.data(),
                            middles.size(), direct, image);
                };
            }
        }
    }
    // A microphone on a wall, a level roof or a slope: that surface is its own
    // end's reflector (raised, by the facade whose coefficient it takes) -- for
    // the turbulent field's doubling at the receiver, and for a plate on a
    // level roof, R in P_d + R P_r (FloorModel::surface, below).
    std::vector<std::vector<double>> receiver_energy = floor_energy;
    std::vector<std::size_t> raised(nm, PointRay::no_surface);
    for (std::size_t m = 0; m < nm; m++)
    {
        const Buildings::Support on = city.support(microphones[m]);
        if (on.kind == Buildings::Support::Kind::none || on.kind == Buildings::Support::Kind::edge)
            continue;
        if (floor_of[m].plate && on.kind != Buildings::Support::Kind::roof)
            throw std::invalid_argument("microphone " + std::to_string(m + 1)
                    + " is flush on a wall or a sloping roof with a plate; a plate lies level, on the floor or a "
                      "level roof");
        raised[m] = on.kind == Buildings::Support::Kind::facade ? on.index : city.roof_facade(on.index);
        std::fill(receiver_energy[m].begin(), receiver_energy[m].end(),
                raised[m] == Buildings::no_facade ? 1.0 : d.reflection[raised[m]]);
    }
    const double radius = database.sphere_radius();
    const double at_one_metre = 20.0 * std::log10(radius);
    // The interval's percentiles, as fractions.
    const double low_q = interval != nullptr ? 0.5 - interval->percent / 200.0 : 0.0, high_q = 1.0 - low_q;

    std::vector<Vec3> sources(ns);
    for (std::size_t k = 0; k < ns; k++)
        sources[k] = states[k].position;
    const std::vector<std::size_t> near_sources = edges_near(city, sources, paths.radius),
                                   near_microphones = edges_near(city, microphones, paths.radius);
    const unsigned workers = std::max(1u, threads);
    const std::size_t block = std::max<std::size_t>(64, 4 * workers);

    // Walls and roofs that scatter (diffuse_facades.hpp) reflect rho of what
    // reaches them and scatter s of that.  Each microphone's diffuse field
    // covers the buildings within its reach, whose walls and roofs then give
    // the rays rho (1 - s); the rest reflect their whole rho, since nothing
    // returns what they would scatter.
    // A relief's walls give the rays their whole rho, which take its order 0
    // at their incidence; by orders they scatter nothing diffusely, and by
    // Lambert's law their random-incidence s, scaled to each path's.
    const std::vector<double> band_s = band_scattering(bands);
    std::vector<double> scatter(reflection.size(), 0.0);
    bool diffuse = false;
    const bool reliefs = !d.facade_reliefs.empty();
    auto relief_at = [&](std::size_t f) { return reliefs ? d.facade_reliefs[f] : nullptr; };
    // A relief's random-incidence s per band, once for every wall it is on.
    std::map<const FacadeRelief*, std::vector<double>> random_s;
    for (std::size_t f = 0; f < d.reflection.size(); f++)
    {
        const FacadeRelief* r = relief_at(f);
        const std::vector<double>* relief_s = nullptr;
        if (r != nullptr && d.relief_lambert)
        {
            auto found = random_s.find(r);
            if (found == random_s.end())
            {
                std::vector<double> v(nb);
                for (std::size_t b = 0; b < nb; b++)
                    v[b] = r->random_incidence(bands[b], sound_speed);
                found = random_s.emplace(r, std::move(v)).first;
            }
            relief_s = &found->second;
        }
        for (std::size_t b = 0; b < nb; b++)
        {
            scatter[f * nb + b] = relief_s != nullptr ? (*relief_s)[b]
                    : reliefs                         ? 0.0
                                                      : d.input.scattering[city.facade_building(f)] * band_s[b];
            diffuse = diffuse || scatter[f * nb + b] > 0.0;
        }
    }
    std::vector<std::unique_ptr<DiffuseFacades>> diffuse_fields(diffuse ? nm : 0);
    std::vector<std::vector<double>> microphone_reflection(diffuse ? nm : 0);
    if (diffuse)
    {
        DiffuseOptions options;
        options.radius = d.diffuse_radius;
        options.patch_size = d.patch_size;
        // The rays reflect in the roofs, so the roofs scatter too.
        options.roofs = paths.roofs;
        options.air = nepers;
        if (reliefs && d.relief_lambert)
            options.scattering_at = [&](std::size_t f, double sin_psi, double cos_elevation, std::size_t b) {
                const FacadeRelief* r = relief_at(f);
                return r != nullptr ? 1.0 - r->specular(sin_psi, bands[b] * cos_elevation, sound_speed) : -1.0;
            };
        // With the orders as rays too, those carry a relief's first scattering
        // to the microphone, and the diffuse field only what it scatters onto
        // other walls.
        if (reliefs && d.relief_lambert && d.relief_orders)
            options.first_scattering_elsewhere = [&](std::size_t f) { return relief_at(f) != nullptr; };
        std::size_t unconverged = 0;
        for (std::size_t m = 0; m < nm; m++)
        {
            diffuse_fields[m] = std::make_unique<DiffuseFacades>(city, microphones[m], reflection, scatter,
                    std::accumulate(floor_energy[m].begin(), floor_energy[m].end(), 0.0) / static_cast<double>(nb), bands,
                    sound_speed, options, workers);
            diffuse_fields[m]->prepare_point_lighting(paths.max_reflections, workers);
            unconverged += diffuse_fields[m]->converged() ? 0 : 1;
            microphone_reflection[m] = reflection;
            for (std::size_t f = 0; f < d.reflection.size(); f++)
                if (diffuse_fields[m]->scatters(f) && relief_at(f) == nullptr)
                    for (std::size_t b = 0; b < nb; b++)
                        microphone_reflection[m][f * nb + b] = reflection[f * nb + b] * (1.0 - scatter[f * nb + b]);
        }
        if (unconverged > 0)
            warn("--buildings: the diffuse field's solve did not converge at ", unconverged, " of ", nm,
                    " microphones.");
    }
    auto reflection_for = [&](std::size_t m) -> const std::vector<double>& {
        return diffuse ? microphone_reflection[m] : reflection;
    };
    // By Lambert's law alone, a relief beyond a microphone's diffuse reach
    // reflects whole, as a plain wall there does: nothing would carry what it
    // scatters.  Each microphone's rays take its own reliefs.
    std::vector<std::vector<const FacadeRelief*>> microphone_reliefs;
    std::vector<PointPathOptions> microphone_paths(nm, paths);
    if (reliefs && d.relief_lambert && !d.relief_orders && diffuse)
    {
        microphone_reliefs.assign(nm, d.facade_reliefs);
        for (std::size_t m = 0; m < nm; m++)
        {
            for (std::size_t f = 0; f < d.facade_reliefs.size(); f++)
                if (!diffuse_fields[m]->scatters(f))
                    microphone_reliefs[m][f] = nullptr;
            microphone_paths[m].reliefs = &microphone_reliefs[m];
        }
    }
    // A plate on a raised surface takes the specular R its rays' partners do.
    for (std::size_t m = 0; m < nm; m++)
        if (floor_of[m].plate && raised[m] != PointRay::no_surface)
        {
            floor_of[m].surface.resize(middles.size());
            for (std::size_t j = 0; j < middles.size(); j++)
                floor_of[m].surface[j] = raised[m] == Buildings::no_facade
                        ? 1.0
                        : std::sqrt(reflection_for(m)[raised[m] * nb + j / 3]);
        }

    // One pair: its band levels over the metric's bands, its A-weighted level
    // and its interval's ends, all -infinity where nothing arrives, and with
    // scattering walls its diffuse energy per band.
    struct Pair
    {
        std::vector<double> bands, diffuse, scattered;
        double level = silence, lower = silence, upper = silence;
    };
    struct Scratch
    {
        std::vector<double> weights, spectrum, full, energy, variance, lighting;
    };
    // Each ray's amplitude at 1 m along its launch, from the sphere at its
    // radius, into w.weights; no spreading gain inside the sphere, and no
    // absorption there.
    auto weigh = [&](const PointSourceField& field, std::size_t state, Scratch& w) {
        const auto& rays = field.rays();
        w.spectrum.resize(all);
        w.weights.assign(rays.size() * nb, 0.0);
        for (std::size_t r = 0; r < rays.size(); r++)
        {
            database.spectrum_along(states[state], point_ray_launch(rays[r]), w.spectrum.data(), first_band,
                    last_band, packed != nullptr ? (*packed)[state] : nullptr);
            const double inside = std::min(rays[r].length, radius);
            for (std::size_t b = 0; b < nb; b++)
            {
                const double level = w.spectrum[first_band + b];
                if (std::isfinite(level))
                    w.weights[r * nb + b] = std::pow(10.0, (level + at_one_metre) / 20.0) * (inside / radius)
                            * std::exp(nepers[b] * inside);
            }
        }
    };
    // The diffuse energy at microphone `m` from the source at `state`, per
    // band, into `out`: each path lighting a patch carries the sphere along
    // its launch, its power at 1 m as a ray's amplitude is (weigh).
    auto diffuse_energy = [&](std::size_t state, std::size_t m, Scratch& w, std::vector<double>& out) {
        w.lighting.resize(all);
        const DiffuseFacades::SourcePower power = [&](const Vec3& launch, double length, double* sent) {
            database.spectrum_along(states[state], launch, w.lighting.data(), first_band, last_band,
                    packed != nullptr ? (*packed)[state] : nullptr);
            const double inside = std::min(length, radius);
            for (std::size_t b = 0; b < nb; b++)
            {
                const double level = w.lighting[first_band + b];
                sent[b] = std::isfinite(level) ? std::pow(10.0, (level + at_one_metre) / 10.0) * (inside / radius)
                                * (inside / radius) * std::exp(2.0 * nepers[b] * inside)
                                               : 0.0;
            }
        };
        out.assign(nb, 0.0);
        diffuse_fields[m]->point_energy(states[state].position, paths.max_reflections, nepers, out.data(), 1,
                &power);
    };
    // With turbulence, the energy it scatters into the shadow at microphone `m`
    // per band, then the energy open ground would give there, into `out`
    // (empty where the source is in sight): the direct path's free-field
    // energy, from the sphere along the straight line as a ray's is (weigh),
    // times the scattered level, its fit's doubling by the ground at either
    // end taken as the floor's ((1 + R) / 2 each, R its energy coefficient or
    // a ground's diffuse_reflection; at a microphone on a wall, a level roof
    // or a slope, that surface's), absorbed along the line.
    const double turbulence_cv2 = velocity_equivalent(d.turbulence, sound_speed);
    auto turbulence_energy = [&](std::size_t state, std::size_t m, Scratch& w, std::vector<double>& out) {
        out.clear();
        if (!(d.turbulence > 0.0))
            return;
        const Vec3& s = states[state].position;
        const ShadowScreen screen = shadow_screen(city, s, microphones[m]);
        if (!screen.shielded)
            return;
        const Vec3 line = microphones[m] - s;
        const double r = line.norm();
        w.lighting.resize(all);
        database.spectrum_along(states[state], line / r, w.lighting.data(), first_band, last_band,
                packed != nullptr ? (*packed)[state] : nullptr);
        const double inside = std::min(r, radius);
        out.assign(2 * nb, 0.0);
        for (std::size_t b = 0; b < nb; b++)
        {
            const double level = w.lighting[first_band + b];
            if (!std::isfinite(level))
                continue;
            const double free = std::pow(10.0, (level + at_one_metre) / 10.0) * (inside / radius) * (inside / radius)
                    * std::exp(-2.0 * nepers[b] * (r - inside)) / (r * r);
            const double floor = floor_energy[m][b], receiver = receiver_energy[m][b];
            out[b] = free * 0.25 * (1.0 + floor) * (1.0 + receiver)
                    * std::pow(10.0, 0.1 * scattered_level(screen.h, screen.d_s, screen.d_r, bands[b], turbulence_cv2));
            out[nb + b] = free * (1.0 + receiver);
        }
    };
    // The scattered energy added into band levels, the total held to open
    // ground's where it is louder.
    auto add_scattered = [&](const std::vector<double>& energy, std::vector<double>& levels_db) {
        if (energy.empty())
            return;
        for (std::size_t b = 0; b < nb; b++)
        {
            const double now = std::pow(10.0, 0.1 * levels_db[b]);
            const double added = std::min(energy[b], std::max(0.0, energy[nb + b] - now));
            if (added > 0.0)
                levels_db[b] = 10.0 * std::log10(now + added);
        }
    };
    // The diffuse energy added into band levels.
    auto add_diffuse = [&](const std::vector<double>& energy, std::vector<double>& levels_db) {
        for (std::size_t b = 0; b < energy.size(); b++)
            if (energy[b] > 0.0)
                levels_db[b] = 10.0 * std::log10(std::pow(10.0, 0.1 * levels_db[b]) + energy[b]);
    };
    // With reliefs scattering by their orders, alone or with the diffuse
    // field (hybrid), the grating rays' energy at
    // microphone `m` from the source at `state` per band, added into `out`
    // (sized when empty): each ray carrying the sphere along its launch, its
    // power at 1 m as a ray's amplitude is (weigh), and each order's floor
    // variants adding with their phases (grating_band_energies).
    const bool gratings = reliefs && d.relief_orders;
    auto grating_energy = [&](const PathEnd& source, const PathEnd& microphone, std::size_t state, std::size_t m,
                                  Scratch& w, std::vector<double>& out) {
        if (!gratings)
            return;
        GratingInput g;
        g.reliefs = &d.facade_reliefs;
        g.bands = bands;
        g.sound_speed = sound_speed;
        g.reflection = &reflection_for(m);
        g.air = nepers;
        g.max_reflections = d.grating_reflections;
        // The hybrid's diffuse field carries what a relief scatters onto
        // other walls: its grating rays reach the microphone directly.
        if (d.relief_lambert)
            g.max_microphone_reflections = 0;
        g.source_reflections = paths.source_reflections;
        const FloorModel& floor = floor_of[m];
        if (floor.impedance())
            g.floor = [&](double sin_grazing, double leg, std::size_t b) {
                const auto q = spherical_reflection_coefficient(sin_grazing, pi * bands[b] * leg / sound_speed,
                        floor.admittance[3 * b + 1], floor.boundary_loss);
                return std::complex<double>(q[0], q[1]);
            };
        else
            g.floor = [&](double, double, std::size_t) { return std::complex<double>(std::sqrt(floor.reflection)); };
        // A plate at the microphone, at each band's middle sub-band.
        if (floor.plate)
        {
            g.plate = [&](double elevation, double azimuth, std::size_t b) {
                thread_local std::vector<std::complex<double>> direct, image;
                direct.resize(3 * nb);
                image.resize(3 * nb);
                floor.plate(elevation, azimuth, direct.data(), image.data());
                return std::make_pair(direct[3 * b + 1], image[3 * b + 1]);
            };
            for (std::size_t b = 0; b < nb && !floor.surface.empty(); b++)
                g.plate_surface.push_back(floor.surface[3 * b + 1]);
        }
        thread_local std::vector<GratingRay> rays;
        thread_local std::vector<double> power;
        rays.clear();
        grating_rays(city, source, microphone, g, rays);
        if (rays.empty())
            return;
        power.resize(rays.size());
        w.lighting.resize(all);
        for (std::size_t i = 0; i < rays.size(); i++)
        {
            database.spectrum_along(states[state], rays[i].launch, w.lighting.data(), first_band, last_band,
                    packed != nullptr ? (*packed)[state] : nullptr);
            const std::size_t b = rays[i].band;
            const double level = w.lighting[first_band + b];
            const double inside = std::min(rays[i].length, radius);
            power[i] = std::isfinite(level) ? std::pow(10.0, (level + at_one_metre) / 10.0) * (inside / radius)
                            * (inside / radius) * std::exp(2.0 * nepers[b] * inside)
                                            : 0.0;
        }
        if (out.empty())
            out.assign(nb, 0.0);
        grating_band_energies(rays, bands, sound_speed, power.data(), out.data());
    };
    auto evaluate = [&](const PathEnd& source, const PathEnd& microphone, std::size_t state, std::size_t m,
                            Scratch& w, Pair& out) {
        out.bands.assign(nb, silence);
        out.level = out.lower = out.upper = silence;
        // Coherent rays have no unknown phase relation, and the diffuse field
        // none either, so neither has variance.
        w.variance.assign(nb, 0.0);
        const PointSourceField field(city, source, microphone, microphone_paths[m]);
        turbulence_energy(state, m, w, out.scattered);
        // The diffuse field and the grating rays carry energy alone.
        out.diffuse.clear();
        if (diffuse)
            diffuse_energy(state, m, w, out.diffuse);
        grating_energy(source, microphone, state, m, w, out.diffuse);
        if (field.rays().empty() && out.diffuse.empty() && out.scattered.empty())
            return;
        if (!field.rays().empty())
        {
            weigh(field, state, w);
            if (coherence.sum == RaySum::coherent)
                field.band_levels(reflection_for(m), floor_of[m], bands, sound_speed, out.bands.data(),
                        sub_band_air, w.weights.data());
            else
            {
                w.energy.resize(3 * nb);
                partially_coherent_energies(field, reflection_for(m), floor_of[m], bands, first_band,
                        sound_speed, sub_band_air, w.weights.data(), coherence, states[state].frame, w.energy.data(),
                        interval != nullptr ? w.variance.data() : nullptr);
                for (std::size_t b = 0; b < nb; b++)
                    out.bands[b] = band_level(bands[b], sub_bands(bands[b], sound_speed), &w.energy[3 * b]);
            }
        }
        add_diffuse(out.diffuse, out.bands);
        add_scattered(out.scattered, out.bands);
        w.full.assign(all, silence);
        std::copy(out.bands.begin(), out.bands.end(), w.full.begin() + static_cast<std::ptrdiff_t>(first_band));
        out.level = weighted_level(w.full.data(), gains);
        if (interval != nullptr)
        {
            // The A-weighted energy's mean and variance (the bands' phases
            // independent) and the gamma distribution with them.
            double mean = 0.0, spread = 0.0;
            for (std::size_t b = 0; b < nb; b++)
            {
                const double g = gains[first_band + b];
                if (std::isfinite(out.bands[b]))
                    mean += g * std::pow(10.0, 0.1 * out.bands[b]);
                spread += g * g * w.variance[b];
            }
            if (mean > 0.0)
            {
                out.lower = 10.0 * std::log10(gamma_quantile(mean, spread, low_q));
                out.upper = 10.0 * std::log10(gamma_quantile(mean, spread, high_q));
            }
        }
    };
    // A pair's results into the outputs.
    auto store = [&](std::size_t state, std::size_t m, const Pair& p) {
        levels[state * nm + m] = p.level;
        if (interval != nullptr)
        {
            interval->lower[state * nm + m] = p.lower;
            interval->upper[state * nm + m] = p.upper;
        }
        if (band_levels != nullptr)
            std::copy(p.bands.begin(), p.bands.end(),
                    band_levels->begin() + static_cast<std::ptrdiff_t>((state * nm + m) * all + first_band));
    };
    // The pairs `todo` (state, microphone), sorted by state, into `out`
    // alike: the sources' ends a block at a time, built in parallel, and
    // the pairs microphone by microphone, so a worker keeps the
    // microphone's end while its pairs last.
    // On a GPU: the workers trace each pair, weigh its rays and pack them
    // without pause, and a stream (gpu/band_stream.hpp) sums them there and
    // turns its energies into the pair's levels.
    const uint32_t gpu_sum = coherence.sum == RaySum::cone ? NK_SUM_CONE
            : coherence.sum == RaySum::energy              ? NK_SUM_ENERGY
                                                           : NK_SUM_COHERENT;
    std::vector<std::array<SubBand, 3>> parts(nb);
    for (std::size_t b = 0; b < nb; b++)
        parts[b] = sub_bands(bands[b], sound_speed);
    auto evaluate_pairs = [&](const std::vector<std::pair<std::size_t, std::size_t>>& todo, std::vector<Pair>& out) {
        out.assign(todo.size(), Pair {});
        for (Pair& p : out)
            p.bands.assign(nb, silence);
        // stream_full first: the stream's callbacks write it, and a stream
        // dropped by an exception drains its queue through them.
        std::vector<double> stream_full(all, silence);
        std::unique_ptr<gpu::BandStream> stream;
        if (device != nullptr)
            stream = std::make_unique<gpu::BandStream>(*device, block, [&](std::size_t i, const double* energy) {
                Pair& p = out[i];
                for (std::size_t b = 0; b < nb; b++)
                    p.bands[b] = band_level(bands[b], parts[b], &energy[3 * b]);
                add_diffuse(p.diffuse, p.bands);
                add_scattered(p.scattered, p.bands);
                std::copy(p.bands.begin(), p.bands.end(), stream_full.begin() + static_cast<std::ptrdiff_t>(first_band));
                p.level = weighted_level(stream_full.data(), gains);
            });
        std::vector<std::unique_ptr<PathEnd>> source_ends;
        std::vector<std::size_t> block_states, slot(ns, 0), work;
        for (std::size_t first = 0; first < todo.size();)
        {
            // The pairs of the next `block` distinct states.
            block_states.clear();
            std::size_t last = first;
            while (last < todo.size()
                    && (block_states.size() < block || todo[last].first == block_states.back()))
            {
                if (block_states.empty() || todo[last].first != block_states.back())
                    block_states.push_back(todo[last].first);
                last++;
            }
            for (std::size_t k = 0; k < block_states.size(); k++)
                slot[block_states[k]] = k;
            source_ends.clear();
            source_ends.resize(block_states.size());
            parallel_chunks(block_states.size(), workers, [&](std::size_t begin, std::size_t end) {
                for (std::size_t k = begin; k < end; k++)
                    source_ends[k] = std::make_unique<PathEnd>(
                            city, sources[block_states[k]], true, paths, near_microphones);
            });
            work.resize(last - first);
            std::iota(work.begin(), work.end(), first);
            std::stable_sort(work.begin(), work.end(),
                    [&](std::size_t a, std::size_t b) { return todo[a].second < todo[b].second; });
            auto run = [&] {
                parallel_chunks(work.size(), workers, [&](std::size_t begin, std::size_t end) {
                    std::unique_ptr<PathEnd> microphone;
                    std::size_t microphone_row = nm;
                    Scratch w;
                    for (std::size_t i = begin; i < end; i++)
                    {
                        const auto [state, m] = todo[work[i]];
                        if (m != microphone_row)
                        {
                            microphone = std::make_unique<PathEnd>(city, microphones[m], false, paths, near_sources);
                            microphone_row = m;
                        }
                        if (!stream)
                        {
                            evaluate(*source_ends[slot[state]], *microphone, state, m, w, out[work[i]]);
                            continue;
                        }
                        const PointSourceField field(city, *source_ends[slot[state]], *microphone,
                                microphone_paths[m]);
                        // The diffuse energy before the push: the stream's
                        // callback adds it.
                        Pair& p = out[work[i]];
                        p.diffuse.clear();
                        if (diffuse)
                            diffuse_energy(state, m, w, p.diffuse);
                        grating_energy(*source_ends[slot[state]], *microphone, state, m, w, p.diffuse);
                        turbulence_energy(state, m, w, p.scattered);
                        if (field.rays().empty())
                        {
                            // The diffuse and scattered fields alone, where
                            // they reach.
                            if (!p.diffuse.empty() || !p.scattered.empty())
                            {
                                add_diffuse(p.diffuse, p.bands);
                                add_scattered(p.scattered, p.bands);
                                w.full.assign(all, silence);
                                std::copy(p.bands.begin(), p.bands.end(),
                                        w.full.begin() + static_cast<std::ptrdiff_t>(first_band));
                                p.level = weighted_level(w.full.data(), gains);
                            }
                            continue;
                        }
                        weigh(field, state, w);
                        gpu::PointRaySet set;
                        gpu::PointRayPacking packing;
                        packing.weights = w.weights.data();
                        packing.sum = gpu_sum;
                        packing.cone = coherence.cone;
                        gpu::pack_point_rays(field, reflection_for(m), floor_of[m], bands, sound_speed, set,
                                sub_band_air, packing);
                        stream->push(work[i], std::move(set));
                    }
                });
            };
            if (!stream)
                run();
            else
                try
                {
                    run();
                }
                catch (...)
                {
                    // The GPU's own failure, where it is the cause.
                    stream->finish();
                    throw;
                }
            first = last;
        }
        if (stream)
            stream->finish();
    };

    if (!(state_tolerance > 0.0) || ns < 3)
    {
        // Every state, a block at a time.
        std::vector<std::pair<std::size_t, std::size_t>> todo;
        std::vector<Pair> out;
        for (std::size_t first = 0; first < ns; first += block)
        {
            todo.clear();
            for (std::size_t state = first; state < std::min(ns, first + block); state++)
                for (std::size_t m = 0; m < nm; m++)
                    todo.emplace_back(state, m);
            evaluate_pairs(todo, out);
            for (std::size_t i = 0; i < todo.size(); i++)
                store(todo[i].first, todo[i].second, out[i]);
        }
        return ns * nm;
    }

    // Along the track, per microphone: the buildings' transfer -- each band's
    // level less the direct ray's in free field, which carries the source --
    // at states every `spacing`, then at the middle of any stretch where the
    // transfer taken linearly between its ends misses the level there by more
    // than state_tolerance, or a band within 30 dB of the loudest A-weighted
    // one by 15 times that, or where a loud band is silent at one and not the
    // other.  The other states take the direct ray's free-field levels, from
    // the sphere at their own state, plus the transfer interpolated.
    constexpr std::size_t spacing = 32;
    constexpr double loud = 30.0, band_factor = 15.0;
    // Each band's A-weighting, dB, for the loud-band tests.
    std::vector<double> gain_db(nb);
    for (std::size_t q = 0; q < nb; q++)
        gain_db[q] = 10.0 * std::log10(gains[first_band + q]);
    // A band more than `buried` under the state's loudest free-field band
    // carries nothing to any level, and its transfer means nothing: the
    // database's top bands sit near -1300 dB at some conditions and jump to
    // -300 at the next state's, and band_level floors a level at -300 dB.
    // Such a band is silent in the transfer and the estimates.
    constexpr double buried = 120.0;
    // The direct ray's band levels in free field, dB: as a lone ray's, its
    // weight over its length and the air's attenuation along it.
    auto direct = [&](std::size_t state, std::size_t m, Scratch& w, double* out) {
        const Vec3 path = microphones[m] - sources[state];
        const double r = path.norm();
        w.spectrum.resize(all);
        database.spectrum_along(states[state], path / r, w.spectrum.data(), first_band, last_band,
                packed != nullptr ? (*packed)[state] : nullptr);
        const double inside = std::min(r, radius);
        double top = silence;
        for (std::size_t b = 0; b < nb; b++)
        {
            const double level = w.spectrum[first_band + b];
            out[b] = std::isfinite(level) ? level + at_one_metre + 20.0 * std::log10(inside / radius)
                            - 20.0 * std::log10(r) - 20.0 / std::log(10.0) * nepers[b] * (r - inside)
                                          : silence;
            top = std::max(top, out[b]);
        }
        for (std::size_t b = 0; b < nb; b++)
            if (!(out[b] >= top - buried))
                out[b] = silence;
    };
    // A band's transfer is NaN where it is not known -- the band buried in
    // the free field there, or its rays' sum at band_level's floor of
    // -300 dB -- and -infinity where no ray reaches.  Between two states' at
    // fraction f: in dB where both are known, in energy where one is silent,
    // the known one where the other is not known, NaN where neither is.
    auto between = [](double a, double b, double f) {
        if (std::isfinite(a) && std::isfinite(b))
            return a + f * (b - a);
        if (std::isnan(a))
            return b;
        if (std::isnan(b))
            return a;
        const double e = (1.0 - f) * (std::isfinite(a) ? std::pow(10.0, 0.1 * a) : 0.0)
                + f * (std::isfinite(b) ? std::pow(10.0, 0.1 * b) : 0.0);
        return e > 0.0 && std::isfinite(e) ? 10.0 * std::log10(e) : silence;
    };
    // A level within a decibel of band_level's floor is at it; a state's
    // bands are heard where one is above it.
    const double band_floor = 10.0 * std::log10(band_level_floor) + 1.0;
    auto heard = [&](const double* bands) {
        return std::any_of(bands, bands + nb, [&](double v) { return v > band_floor; });
    };
    struct Known
    {
        Pair pair;
        std::vector<double> transfer;
    };
    std::vector<std::map<std::size_t, Known>> known(nm);
    auto take = [&](const std::vector<std::pair<std::size_t, std::size_t>>& todo, std::vector<Pair>& out) {
        Scratch w;
        std::vector<double> free(nb);
        for (std::size_t i = 0; i < todo.size(); i++)
        {
            const auto [state, m] = todo[i];
            direct(state, m, w, free.data());
            Known k;
            k.transfer.resize(nb);
            for (std::size_t b = 0; b < nb; b++)
            {
                const double level = out[i].bands[b];
                k.transfer[b] = !std::isfinite(free[b]) || (std::isfinite(level) && level <= band_floor)
                        ? std::numeric_limits<double>::quiet_NaN()
                        : std::isfinite(level) ? level - free[b]
                                               : silence;
            }
            k.pair = std::move(out[i]);
            known[m][state] = std::move(k);
        }
    };
    // The estimate at `state` between two known states: band levels into
    // `out` (the metric's bands), its A-weighted level returned.
    auto estimate = [&](std::size_t state, std::size_t m, const Known& a, std::size_t sa, const Known& b,
                            std::size_t sb, Scratch& w, double* out) {
        direct(state, m, w, out);
        const double f = static_cast<double>(state - sa) / static_cast<double>(sb - sa);
        // A band whose transfer is known at neither end (buried there) takes
        // the nearest band's that is known at both: the transfer varies
        // smoothly across bands, and such a band is 30 dB or more under the
        // loudest here, or the state would have been evaluated.  What is
        // still not known is silent.
        auto known_both = [&](std::size_t q) { return std::isfinite(a.transfer[q]) && std::isfinite(b.transfer[q]); };
        for (std::size_t q = 0; q < nb; q++)
        {
            std::size_t from = q;
            if (std::isnan(a.transfer[q]) && std::isnan(b.transfer[q]))
                for (std::size_t step = 1; step < nb; step++)
                {
                    if (q >= step && known_both(q - step))
                    {
                        from = q - step;
                        break;
                    }
                    if (q + step < nb && known_both(q + step))
                    {
                        from = q + step;
                        break;
                    }
                }
            const double t = between(a.transfer[from], b.transfer[from], f);
            const double v = out[q] + t;
            out[q] = std::isfinite(v) ? v : silence;
        }
        w.full.assign(all, silence);
        std::copy(out, out + nb, w.full.begin() + static_cast<std::ptrdiff_t>(first_band));
        return weighted_level(w.full.data(), gains);
    };

    std::vector<std::pair<std::size_t, std::size_t>> todo;
    std::vector<Pair> out;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> pending(nm), accepted(nm);
    for (std::size_t state = 0; state < ns; state += spacing)
        for (std::size_t m = 0; m < nm; m++)
            todo.emplace_back(state, m);
    if ((ns - 1) % spacing != 0)
        for (std::size_t m = 0; m < nm; m++)
            todo.emplace_back(ns - 1, m);
    std::size_t evaluated = 0;
    for (std::size_t m = 0; m < nm; m++)
        for (std::size_t a = 0; a + 1 < ns; a = std::min(a + spacing, ns - 1))
            pending[m].emplace_back(a, std::min(a + spacing, ns - 1));
    while (!todo.empty())
    {
        std::sort(todo.begin(), todo.end());
        evaluate_pairs(todo, out);
        evaluated += todo.size();
        take(todo, out);
        // Each pending stretch's middle is known now: split the stretch where
        // the estimate there misses, and ask for the middles of the halves.
        todo.clear();
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> next(nm);
        Scratch w;
        std::vector<double> guess(nb);
        for (std::size_t m = 0; m < nm; m++)
        {
            for (const auto& [a, b] : pending[m])
            {
                if (b - a < 2)
                {
                    accepted[m].emplace_back(a, b);
                    continue;
                }
                const std::size_t mid = a + (b - a) / 2;
                const auto known_mid = known[m].find(mid);
                if (known_mid == known[m].end())
                {
                    // Its middle is asked for in this round.
                    next[m].emplace_back(a, b);
                    todo.emplace_back(mid, m);
                    continue;
                }
                const Known& at = known_mid->second;
                const double level = estimate(mid, m, known[m].at(a), a, known[m].at(b), b, w, guess.data());
                // A middle that is silent (no ray) or floored, and an estimate
                // that is too, agree; one heard and the other not do not.
                const double* bands = at.pair.bands.data();
                const bool middle_heard = heard(bands);
                bool good = heard(guess.data()) == middle_heard
                        && (!middle_heard || std::abs(level - at.pair.level) <= state_tolerance);
                // The heard bands within `loud` of the loudest A-weighted one.
                double top = silence;
                for (std::size_t q = 0; q < nb; q++)
                    if (bands[q] > band_floor)
                        top = std::max(top, bands[q] + gain_db[q]);
                for (std::size_t q = 0; q < nb && good; q++)
                    if (bands[q] > band_floor && bands[q] + gain_db[q] >= top - loud)
                        good = std::isfinite(guess[q]) && std::abs(guess[q] - bands[q]) <= band_factor * state_tolerance;
                if (good)
                {
                    // Its middle, known, splits it for the estimates.
                    accepted[m].emplace_back(a, mid);
                    accepted[m].emplace_back(mid, b);
                }
                else
                {
                    next[m].emplace_back(a, mid);
                    next[m].emplace_back(mid, b);
                    // Each half's middle, unless the half is a step.
                    for (const auto& [lo, hi] : {std::make_pair(a, mid), std::make_pair(mid, b)})
                        if (hi - lo >= 2)
                            todo.emplace_back(lo + (hi - lo) / 2, m);
                }
            }
            // Every pending stretch of two steps or more has its middle
            // asked for; a step has no state inside to estimate.
            pending[m] = std::move(next[m]);
        }
    }

    // A state where the source makes a band loud whose transfer is known at
    // neither end of its stretch -- buried there, or reached by no ray -- is
    // evaluated outright: the database's top bands jump by 1000 dB for a
    // state now and then, and an estimate would leave such a band silent.
    {
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> unknown(nm);
        parallel_chunks(nm, workers, [&](std::size_t begin, std::size_t end) {
            Scratch w;
            std::vector<double> free(nb);
            for (std::size_t m = begin; m < end; m++)
                for (const auto& [a, b] : accepted[m])
                {
                    const Known& ka = known[m].at(a);
                    const Known& kb = known[m].at(b);
                    for (std::size_t state = a + 1; state < b; state++)
                    {
                        direct(state, m, w, free.data());
                        double top = silence;
                        for (std::size_t q = 0; q < nb; q++)
                            top = std::max(top, free[q] + gain_db[q]);
                        for (std::size_t q = 0; q < nb; q++)
                            if (free[q] + gain_db[q] >= top - loud
                                    && std::isnan(ka.transfer[q]) && std::isnan(kb.transfer[q]))
                            {
                                unknown[m].emplace_back(state, m);
                                break;
                            }
                    }
                }
        });
        todo.clear();
        for (const auto& list : unknown)
            todo.insert(todo.end(), list.begin(), list.end());
        if (!todo.empty())
        {
            std::sort(todo.begin(), todo.end());
            evaluate_pairs(todo, out);
            evaluated += todo.size();
            take(todo, out);
        }
    }

    // Every state: a known one as it was evaluated, the others estimated
    // within their accepted stretch.
    parallel_chunks(nm, workers, [&](std::size_t begin, std::size_t end) {
        Scratch w;
        std::vector<double> guess(nb);
        Pair p;
        for (std::size_t m = begin; m < end; m++)
        {
            for (const auto& [state, k] : known[m])
                store(state, m, k.pair);
            for (const auto& [a, b] : accepted[m])
            {
                const Known& ka = known[m].at(a);
                const Known& kb = known[m].at(b);
                for (std::size_t state = a + 1; state < b; state++)
                {
                    if (known[m].count(state))
                        continue;
                    const double f = static_cast<double>(state - a) / static_cast<double>(b - a);
                    p.level = estimate(state, m, ka, a, kb, b, w, guess.data());
                    p.bands.assign(guess.begin(), guess.end());
                    // The interval's ends as offsets from the level, taken alike.
                    auto offset = [&](const Known& k, double end) {
                        return std::isfinite(end) && std::isfinite(k.pair.level) ? end - k.pair.level
                                                                                : std::numeric_limits<double>::quiet_NaN();
                    };
                    const double lo_a = offset(ka, ka.pair.lower), lo_b = offset(kb, kb.pair.lower);
                    const double hi_a = offset(ka, ka.pair.upper), hi_b = offset(kb, kb.pair.upper);
                    p.lower = std::isfinite(p.level) && std::isfinite(lo_a) && std::isfinite(lo_b)
                            ? p.level + lo_a + f * (lo_b - lo_a)
                            : silence;
                    p.upper = std::isfinite(p.level) && std::isfinite(hi_a) && std::isfinite(hi_b)
                            ? p.level + hi_a + f * (hi_b - hi_a)
                            : silence;
                    store(state, m, p);
                }
            }
        }
    });
    return evaluated;
}

} // namespace niceops
