#include "tones.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>

#include "util/constants.hpp"
#include "util/logger.hpp"
#include "weighting.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace niceops {

std::optional<ToneCoherence> parse_tone_coherence(const std::string& name)
{
    if (name == "coherent")
        return ToneCoherence::coherent;
    if (name == "energy")
        return ToneCoherence::energy;
    return std::nullopt;
}

std::optional<ToneGround> parse_tone_ground(const std::string& name)
{
    if (name == "frequency")
        return ToneGround::frequency;
    if (name == "band")
        return ToneGround::band;
    return std::nullopt;
}

namespace {

std::string trimmed(const std::string& s)
{
    const auto a = s.find_first_not_of(" \t");
    if (a == std::string::npos)
        return "";
    const auto b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> comma_list(const std::string& text)
{
    std::vector<std::string> out;
    std::stringstream in(text);
    std::string item;
    while (std::getline(in, item, ','))
    {
        item = trimmed(item);
        if (!item.empty())
            out.push_back(item);
    }
    return out;
}

} // namespace

void parse_rotor_azimuths(const std::string& text, std::map<std::string, double>& out)
{
    for (const auto& item : comma_list(text))
    {
        const auto eq = item.find('=');
        if (eq == std::string::npos || eq == 0 || eq + 1 == item.size())
            input_error("--rotor_azimuth '", item, "': give each source as tag=degrees.");
        const std::string tag = trimmed(item.substr(0, eq));
        const std::string value = trimmed(item.substr(eq + 1));
        std::size_t used = 0;
        double deg = 0.0;
        try
        {
            deg = std::stod(value, &used);
        }
        catch (const std::exception&)
        {
            used = 0;
        }
        if (used != value.size() || !std::isfinite(deg))
            input_error("--rotor_azimuth '", item, "': '", value, "' is not a number of degrees.");
        if (out.count(tag) > 0)
            input_error("--rotor_azimuth names source ", tag, " twice.");
        out[tag] = deg;
    }
}

std::vector<std::string> parse_sync_group(const std::string& text)
{
    auto members = comma_list(text);
    if (members.empty())
        input_error("--sync_group '", text, "' names no source: give the group's tags, comma-separated.");
    return members;
}

ToneArrival tone_arrival(const Vec3& hub, const Vec3& velocity, const Vec3& point, double reception_lag,
        double sound_speed, const Vec3& wind)
{
    // The emission time delta solves c (T - delta) = |point - (hub + V delta)|, T the
    // reception lag: (c^2 - V^2) delta^2 - 2 (c^2 T - D.V) delta + (c^2 T^2 - D^2) = 0,
    // D = point - hub.  Of its roots the one before reception, written in the
    // form that does not cancel when delta is small against T.
    //
    // In a uniform wind W the same holds in the frame moving with the air, where
    // sound travels at c: there the hub moves at V - W and the point, at rest on
    // the ground, at -W, so it is received at point - W T.  Calm air takes the
    // first branch, the arithmetic as it always was.
    const double c = sound_speed;
    const bool calm = wind.x == 0.0 && wind.y == 0.0 && wind.z == 0.0;
    const Vec3 v = calm ? velocity : velocity - wind;
    const Vec3 received = calm ? point : point - wind * reception_lag;
    const Vec3 d0 = received - hub;
    const double a = c * c - v.squaredNorm();
    const double b = c * c * reception_lag - d0.dot(v);
    const double cc = c * c * reception_lag * reception_lag - d0.squaredNorm();
    const double disc = std::max(b * b - a * cc, 0.0);
    ToneArrival out;
    out.delta = cc / (b + std::sqrt(disc));
    out.length = c * (reception_lag - out.delta);
    const Vec3 emitted = hub + v * out.delta;
    const Vec3 toward = received - emitted;
    const double n = toward.norm();
    out.direction = n > 0.0 ? toward / n : Vec3(0.0, 0.0, -1.0);
    // f = f_e (1 - n.v_o / c) / (1 - n.v_s / c) with the air-frame velocities
    // v_s = V - W and v_o = -W; doppler is f_e / f.
    out.doppler = calm ? 1.0 - velocity.dot(out.direction) / c
                       : (1.0 - v.dot(out.direction) / c) / (1.0 + wind.dot(out.direction) / c);
    return out;
}

double add_tone_power(double level, const double* power, const double* weighted, const double* ratio,
        std::size_t first_band, std::size_t last_band, double* bands_out)
{
    double total = level > -std::numeric_limits<double>::infinity() ? std::pow(10.0, 0.1 * level) : 0.0;
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double r = ratio != nullptr ? ratio[band] : 1.0;
        total += weighted[band] * r;
        if (bands_out != nullptr && power[band] > 0.0)
            bands_out[band] = 10.0 * std::log10(std::pow(10.0, 0.1 * bands_out[band]) + power[band] * r);
    }
    return total > 0.0 ? 10.0 * std::log10(total) : -std::numeric_limits<double>::infinity();
}

ToneEvaluator::ToneEvaluator(const Database& database, const ToneRunOptions& options,
        const std::vector<double>& band_gains, double sound_speed)
        : database_(database), options_(options), band_gains_(band_gains), sound_speed_(sound_speed)
{
    const auto& sources = database.tone_sources();
    const auto& centers = database.band_centers();
    if (sources.empty())
        input_error("The pure-tone evaluation was asked of a database with no tone sources.");
    if (centers.empty())
        input_error("The database's pure tones need a band spectrum to be filed into, and it carries none.");
    std::map<std::string, std::size_t> index;
    for (std::size_t s = 0; s < sources.size(); s++)
        index[sources[s].tag] = s;
    auto tags = [&]() {
        std::ostringstream all;
        for (const auto& s : sources)
            all << (all.tellp() > 0 ? ", " : "") << s.tag;
        return all.str();
    };

    // Each rotor at the azimuth its phases were locked to (the first sphere's),
    // which reproduces the phasing it was predicted with; zero, the stored
    // rotor-locked phase itself, where the database does not record one.
    azimuth_.assign(sources.size(), 0.0);
    for (std::size_t s = 0; s < sources.size(); s++)
        if (!sources[s].azimuths_deg.empty() && std::isfinite(sources[s].azimuths_deg.front()))
            azimuth_[s] = sources[s].azimuths_deg.front() * pi / 180.0;
    for (const auto& [tag, deg] : options.rotor_azimuth_deg)
    {
        const auto it = index.find(tag);
        if (it == index.end())
            input_error("--rotor_azimuth names source ", tag, ", which the database does not carry (it has ", tags(),
                    ").");
        azimuth_[it->second] = deg * pi / 180.0;
    }

    const bool coherent = options.coherence == ToneCoherence::coherent;
    group_of_.assign(sources.size(), static_cast<std::size_t>(-1));
    // Run-time groups first: they take their members out of the database's.
    for (const auto& members : options.sync_groups)
    {
        Group g;
        g.forced = true;
        g.coherent = coherent;
        for (const auto& tag : members)
        {
            const auto it = index.find(tag);
            if (it == index.end())
                input_error("--sync_group names source ", tag, ", which the database does not carry (it has ", tags(),
                        ").");
            const std::size_t s = it->second;
            if (group_of_[s] != static_cast<std::size_t>(-1))
                input_error("--sync_group puts source ", tag, " in two groups.");
            if (sources[s].phase_reference != PhaseReference::rotor_azimuth)
                input_error("--sync_group: source ", tag, "'s phases are ", to_string(sources[s].phase_reference),
                        ", not tied to its rotor's azimuth, so it cannot add coherently with another rotor.");
            group_of_[s] = groups_.size();
            g.sources.push_back(s);
            g.name += (g.name.empty() ? "" : ",") + tag;
        }
        groups_.push_back(g);
    }
    // The database's groups, among the sources left.
    std::map<std::string, std::size_t> named;
    for (std::size_t s = 0; s < sources.size(); s++)
    {
        if (group_of_[s] != static_cast<std::size_t>(-1))
            continue;
        const auto& src = sources[s];
        if (!src.sync_group.empty() && src.phase_reference != PhaseReference::rotor_azimuth)
            warn("Database source ", src.tag, " is in sync group ", src.sync_group, " but its phases are ",
                    to_string(src.phase_reference), ", not tied to its rotor's azimuth: it adds in energy.");
        if (src.sync_group.empty() || src.phase_reference != PhaseReference::rotor_azimuth)
        {
            Group g;
            g.name = src.tag;
            g.sources = {s};
            group_of_[s] = groups_.size();
            groups_.push_back(g);
            continue;
        }
        auto found = named.find(src.sync_group);
        if (found == named.end())
        {
            Group g;
            g.name = src.sync_group;
            g.coherent = coherent;
            found = named.emplace(src.sync_group, groups_.size()).first;
            groups_.push_back(g);
        }
        group_of_[s] = found->second;
        groups_[found->second].sources.push_back(s);
    }
    // A database group must turn at one rate: its members' phases are fixed to
    // each other only if their shafts are.
    const std::size_t n_conditions = database.condition_count();
    for (const auto& g : groups_)
    {
        if (g.forced || g.sources.size() < 2)
            continue;
        for (std::size_t i = 0; i < n_conditions; i++)
        {
            const double first = std::fabs(sources[g.sources[0]].shaft_rates[i]);
            for (const std::size_t s : g.sources)
            {
                const double rate = std::fabs(sources[s].shaft_rates[i]);
                if (!(std::fabs(rate - first) <= 1e-6 * first))
                    input_error("Database sync group ", g.name, ": sources ", sources[g.sources[0]].tag, " and ",
                            sources[s].tag, " turn at ", first, " and ", rate, " rad/s at sphere ", i,
                            ". Sources held at a fixed phase must share a shaft rate.");
            }
        }
        // The default phasing is the first sphere's recorded azimuths; a group
        // whose members' recorded azimuths differ from one sphere to another by
        // other than a common turn is phased differently there, which one
        // default cannot hold.
        for (const std::size_t s : g.sources)
        {
            const auto& a = sources[s].azimuths_deg;
            const auto& a0 = sources[g.sources[0]].azimuths_deg;
            for (std::size_t i = 1; i < n_conditions && i < a.size() && i < a0.size(); i++)
            {
                const double turn = std::remainder((a[i] - a0[i]) - (a[0] - a0[0]), 360.0);
                if (std::isfinite(turn) && std::fabs(turn) > 1e-3)
                    warn("Database sync group ", g.name, ": source ", sources[s].tag, "'s recorded azimuth relative to ",
                            sources[g.sources[0]].tag, "'s is ", std::remainder(a[i] - a0[i], 360.0), " deg at sphere ",
                            i, " and ", std::remainder(a[0] - a0[0], 360.0), " deg at sphere 0; the run phases the "
                            "group as sphere 0 records unless --rotor_azimuth says otherwise.");
            }
        }
    }

    tone_offsets_.resize(sources.size());
    for (std::size_t s = 0; s < sources.size(); s++)
    {
        tone_offsets_[s] = tone_total_;
        tone_total_ += sources[s].tone_count();
        max_tones_ = std::max(max_tones_, sources[s].tone_count());
        for (std::size_t i = 0; i < n_conditions; i++)
            for (const int m : sources[s].harmonics)
            {
                const double rad = (m + 1) * sound_speed / (2.0 * std::fabs(sources[s].shaft_rates[i])
                                                                    * database.sphere_radius());
                const double deg = rad * 180.0 / pi;
                if (deg > near_field_degrees_)
                {
                    near_field_degrees_ = deg;
                    near_field_order_ = m;
                    near_field_source_ = sources[s].tag;
                }
            }
    }

    // Band edges halfway (in log frequency) between centers, and half a step out
    // beyond the first and last.
    for (const double f : centers)
        log2_centers_.push_back(std::log2(f));
    const std::size_t n = log2_centers_.size();
    const double outer = n > 1 ? 0.5 * (log2_centers_[1] - log2_centers_[0]) : 1.0 / 6.0;
    const double outer_high = n > 1 ? 0.5 * (log2_centers_[n - 1] - log2_centers_[n - 2]) : 1.0 / 6.0;
    log2_edges_.push_back(log2_centers_[0] - outer);
    for (std::size_t b = 1; b < n; b++)
        log2_edges_.push_back(0.5 * (log2_centers_[b - 1] + log2_centers_[b]));
    log2_edges_.push_back(log2_centers_[n - 1] + outer_high);
}

std::size_t ToneEvaluator::band_of(double frequency) const
{
    if (!(frequency > 0.0))
        return no_band;
    const double x = std::log2(frequency);
    if (x < log2_edges_.front() || x >= log2_edges_.back())
        return no_band;
    const auto it = std::upper_bound(log2_edges_.begin(), log2_edges_.end(), x);
    return static_cast<std::size_t>(it - log2_edges_.begin()) - 1;
}

double ToneEvaluator::absorption_at(const double* absorption, double frequency) const
{
    // The band centers' coefficients, interpolated in log frequency -- in log
    // coefficient where both neighbors absorb, which follows the f^2 of the
    // classical and relaxation terms -- and held at the ends.
    const std::size_t n = log2_centers_.size();
    const double x = std::log2(frequency);
    if (n == 1 || x <= log2_centers_.front())
        return absorption[0];
    if (x >= log2_centers_.back())
        return absorption[n - 1];
    const auto it = std::upper_bound(log2_centers_.begin(), log2_centers_.end(), x);
    const std::size_t hi = static_cast<std::size_t>(it - log2_centers_.begin());
    const std::size_t lo = hi - 1;
    const double t = (x - log2_centers_[lo]) / (log2_centers_[hi] - log2_centers_[lo]);
    const double a = absorption[lo];
    const double b = absorption[hi];
    if (a > 0.0 && b > 0.0)
        return std::exp((1.0 - t) * std::log(a) + t * std::log(b));
    return (1.0 - t) * a + t * b;
}

ToneEvaluator::StateTones ToneEvaluator::prepare(const Database::SourceState& state, const Vec3& velocity,
        const Vec3& wind) const
{
    const auto& sources = database_.tone_sources();
    StateTones out;
    out.velocity = velocity;
    out.wind = wind;
    out.hubs.resize(sources.size());
    out.frequencies.assign(tone_total_, 0.0);
    out.group_rates.assign(groups_.size(), 0.0);
    const auto& frame = state.frame;
    for (std::size_t s = 0; s < sources.size(); s++)
    {
        const auto c = database_.tone_center(s, state.blend);
        out.hubs[s] = state.position + c[0] * frame.nose + c[1] * frame.starboard + c[2] * frame.normal;
        database_.tone_frequencies(s, state.blend, out.frequencies.data() + tone_offsets_[s]);
    }
    for (std::size_t g = 0; g < groups_.size(); g++)
    {
        if (!groups_[g].forced)
            continue;
        double sum = 0.0;
        for (const std::size_t s : groups_[g].sources)
            sum += std::fabs(database_.tone_shaft_rate(s, state.blend));
        const double rate = sum / static_cast<double>(groups_[g].sources.size());
        out.group_rates[g] = rate;
        for (const std::size_t s : groups_[g].sources)
            for (std::size_t t = 0; t < sources[s].tone_count(); t++)
                out.frequencies[tone_offsets_[s] + t] = sources[s].harmonics[t] * rate / (2.0 * pi);
    }
    return out;
}

void ToneEvaluator::evaluate(const Database::SourceState& state, const StateTones& tones, const Vec3& observer,
        const double* absorption, const Ground* ground, double ground_plane, double* power, double* weighted) const
{
    const auto& sources = database_.tone_sources();
    const double c = sound_speed_;
    const double radius = database_.sphere_radius();
    const double start = database_.absorption_start();
    const bool lossless = database_.sphere_levels_lossless();
    const bool with_ground = ground != nullptr && ground->model != nullptr;
    // One reception time for every arrival: the state's center's direct lag.
    const double lag = (observer - state.position).norm() / c;
    Vec3 image = observer;
    image.z = 2.0 * ground_plane - observer.z;
    const double receiver_height = observer.z - ground_plane;
    // The turbulence coherence from the state's center, the geometry the band
    // model takes, so every source of one state shares it.
    const double center_height = state.position.z - ground_plane;
    const double center_distance = (observer - state.position).horizontal_norm();
    auto coherence = [&](double frequency) {
        return with_ground ? turbulence_coherence(*ground->model, center_height, center_distance, frequency, c) : 1.0;
    };
    auto file = [&](double p, double frequency) {
        const std::size_t band = band_of(frequency);
        if (band == no_band)
        {
            outside_bands_.fetch_add(1, std::memory_order_relaxed);
            return;
        }
        power[band] += p;
        if (band_gains_[band] != 0.0)
            weighted[band] += p * std::pow(10.0, 0.1 * a_weighting(frequency));
    };
    // |d|^2 + |r|^2 + 2 gamma Re(d conj r): one form for a source's own arrivals and a
    // coherent group's sums alike, so a group of one is the source to the bit.
    auto combine = [](std::complex<double> d, std::complex<double> r, double gamma, bool cross) {
        double p = std::norm(d) + std::norm(r);
        if (cross)
            p += 2.0 * gamma * (d * std::conj(r)).real();
        return p;
    };
    auto phasor = [&](const Database::ToneLookup& look, const ToneArrival& a, double frequency_e, double psi,
                          double& received) {
        received = frequency_e / a.doppler;
        if (look.level_db <= Database::no_energy_floor_db)
            return std::complex<double>(0.0, 0.0);
        const double alpha = absorption_at(absorption, received);
        const double absorbed = std::max(a.length - start, lossless ? radius : 0.0);
        const double amplitude = std::pow(10.0, (look.level_db - alpha * absorbed) / 20.0)
                * (radius / std::max(a.length, radius));
        return std::polar(amplitude, look.phase + psi + 2.0 * pi * frequency_e * a.delta);
    };
    auto reflect = [&](std::complex<double> p, const ToneArrival& a, const Vec3& hub, double received) {
        // The ground model's Q at the reflected arrival's own frequency and
        // geometry, conjugated into this convention, with roughness's factor.
        const double source_height = (hub + tones.velocity * a.delta).z - ground_plane;
        const double cos_grazing = a.length > 0.0 ? (source_height + receiver_height) / a.length : 1.0;
        const auto q = spherical_reflection_coefficient(cos_grazing, pi * received * a.length / c,
                ground->admittance->at(received), ground->model->boundary_loss_correction);
        const double rough = roughness_coherence(*ground->model, cos_grazing, received, c);
        return p * std::complex<double>(q[0], -q[1]) * rough;
    };

    thread_local std::vector<Database::ToneLookup> direct, reflected;
    direct.resize(max_tones_);
    reflected.resize(max_tones_);
    struct Sum
    {
        int order;
        std::complex<double> d, r;
        double frequency;
    };
    thread_local std::vector<Sum> sums;
    for (const auto& g : groups_)
    {
        sums.clear();
        for (const std::size_t s : g.sources)
        {
            const auto& src = sources[s];
            const Vec3& hub = tones.hubs[s];
            const ToneArrival d_arrival = tone_arrival(hub, tones.velocity, observer, lag, c, tones.wind);
            database_.tone_directivity(s, state, d_arrival.direction, direct.data());
            ToneArrival r_arrival;
            if (with_ground)
            {
                r_arrival = tone_arrival(hub, tones.velocity, image, lag, c, tones.wind);
                database_.tone_directivity(s, state, r_arrival.direction, reflected.data());
            }
            for (std::size_t t = 0; t < src.tone_count(); t++)
            {
                const double fe = tones.frequencies[tone_offsets_[s] + t];
                const int m = src.harmonics[t];
                const double psi = m * azimuth_[s];
                double fd = 0.0, fr = 0.0;
                const std::complex<double> pd = phasor(direct[t], d_arrival, fe, psi, fd);
                std::complex<double> pr(0.0, 0.0);
                if (with_ground)
                    pr = reflect(phasor(reflected[t], r_arrival, fe, psi, fr), r_arrival, hub, fr);
                // A tone whose phase is known adds its own arrivals coherently;
                // one without (none, or the grid does not resolve it) in energy.
                const bool phased = src.phase_reference != PhaseReference::none && !src.energy_only[t];
                if (g.coherent && phased)
                {
                    auto it = std::find_if(sums.begin(), sums.end(), [m](const Sum& x) { return x.order == m; });
                    if (it == sums.end())
                    {
                        sums.push_back({m, {0.0, 0.0}, {0.0, 0.0}, fd});
                        it = sums.end() - 1;
                    }
                    it->d += pd;
                    it->r += pr;
                    continue;
                }
                file(combine(pd, pr, coherence(fd), phased), fd);
            }
        }
        for (const Sum& x : sums)
            file(combine(x.d, x.r, coherence(x.frequency), true), x.frequency);
    }
}

} // namespace niceops
