#include "gpu/band_energy.hpp"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>

#include "gpu/kernels/band_energy_types.h"
#include "gpu/point_rays.hpp"
#include "util/constants.hpp"

namespace niceops::gpu {

namespace {

// Gauss-Legendre nodes and weights on [-1, 1], (t, w) interleaved, by Newton's
// method on P_m from Tricomi's estimate; kept once made.
const std::vector<double>& gauss_legendre(std::size_t m)
{
    static std::map<std::size_t, std::vector<double>> rules;
    static std::mutex lock;
    std::lock_guard<std::mutex> guard(lock);
    if (auto it = rules.find(m); it != rules.end())
        return it->second;
    std::vector<double> rule(2 * m);
    for (std::size_t i = 0; i < (m + 1) / 2; i++)
    {
        double z = std::cos(pi * (i + 0.75) / (m + 0.5)), slope = 0.0;
        for (int iteration = 0; iteration < 100; iteration++)
        {
            double p1 = 1.0, p2 = 0.0;
            for (std::size_t j = 1; j <= m; j++)
            {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            slope = m * (z * p1 - p2) / (z * z - 1.0);
            const double step = p1 / slope;
            z -= step;
            if (std::abs(step) < 1e-15)
                break;
        }
        const double w = 2.0 / ((1.0 - z * z) * slope * slope);
        rule[2 * i] = -z;
        rule[2 * i + 1] = w;
        rule[2 * (m - 1 - i)] = z;
        rule[2 * (m - 1 - i) + 1] = w;
    }
    return rules.emplace(m, std::move(rule)).first->second;
}

// Nodes enough that the rule integrates exp(i x t) over [-1, 1] to round-off
// for every |x| up to x_max, the integrand's widest frequency (h times the
// rays' spread of lengths): within 2e-12 dB of every pair summed, on the
// district's problems.
std::size_t nodes_for(double x_max)
{
    return static_cast<std::size_t>(std::ceil(0.5 * x_max + 4.0 * std::cbrt(std::max(x_max, 1.0)) + 6.0));
}

// A problem with more nodes than this is summed pair by pair, whatever it
// costs.
constexpr std::size_t max_nodes = 4096;

// At least one element even for none: a dispatch that runs binds every
// buffer it names, used or not, and an unallocated one cannot be bound.
template <typename T>
void ensure(Device& device, Buffer& buffer, std::size_t count)
{
    const std::size_t bytes = std::max<std::size_t>(count, 1) * sizeof(T);
    if (buffer.size() < bytes)
        buffer = device.buffer(std::max(bytes, 2 * buffer.size()));
}

} // namespace

void RaySet::add(const std::complex<double>* a, double band_k1, double band_k2)
{
    for (std::size_t i = 0; i < length.size(); i++)
        amplitude.emplace_back(static_cast<float>(a[i].real()), static_cast<float>(a[i].imag()));
    k1.push_back(band_k1);
    k2.push_back(band_k2);
}

void RaySet::plan()
{
    plan_nodes(length, k1, k2, nodes);
}

void plan_nodes(const std::vector<double>& length, const std::vector<double>& k1, const std::vector<double>& k2,
        std::vector<uint32_t>& nodes)
{
    assert(std::is_sorted(length.begin(), length.end()));
    const std::size_t n = length.size();
    nodes.assign(k1.size(), 0);
    for (std::size_t q = 0; q < k1.size(); q++)
    {
        const double h = 0.5 * (k2[q] - k1[q]);
        if (n == 0 || h <= 0.0)
            continue;
        const double reach = 30.0 / h;
        const std::size_t m = nodes_for(h * (length[n - 1] - length[0]));
        if (m > max_nodes)
            continue;
        // The pairs within reach, which the pair sum visits, counted until
        // they outnumber the nodes' work.  A node costs a sincos and a
        // complex multiply-add per ray, a pair a multiply-add and a divide: on
        // the GPU the two are about even.
        const double work = static_cast<double>(m) * static_cast<double>(n);
        // Nor more than every pair and each ray itself.
        if (work >= 0.5 * static_cast<double>(n) * static_cast<double>(n + 1))
            continue;
        double visited = 0.0;
        for (std::size_t i = 0, j = 0; i < n && visited <= work; i++)
        {
            while (j < n && length[j] - length[i] <= reach)
                j++;
            visited += static_cast<double>(j - i);
        }
        if (work < visited)
            nodes[q] = static_cast<uint32_t>(m);
    }
}

namespace {

// The kernels' sources, one program: band_energy.kern needs only its types,
// point_rays.kern the rest.
const std::vector<std::string> program {"band_energy_types.h", "point_rays_types.h", "wedge.h", "band_energy.kern",
        "point_rays.kern"};

} // namespace

BandEnergies::BandEnergies(Device& device)
    : device_(device)
    , phasors_(device.kernel(program, "band_phasors"))
    , pair_sum_(device.kernel(program, "band_pair_sum"))
    , quadrature_(device.kernel(program, "band_quadrature"))
    , ray_amplitudes_(device.kernel(program, "ray_amplitudes"))
    , face_amplitudes_(device.kernel(program, "face_amplitudes"))
    , face_ratios_(device.kernel(program, "face_ratios"))
{
}

std::size_t BandEnergies::add(const RaySet& set)
{
    if (submitted_)
        throw std::logic_error("BandEnergies::add after submit, before clear");
    if (set.nodes.size() != set.problems() || set.amplitude.size() != set.problems() * set.length.size())
        throw std::logic_error("BandEnergies::add: the set is not planned, or its amplitudes do not fill it");
    sets_.push_back({&set, nullptr});
    const std::size_t first = problem_count_;
    problem_count_ += set.problems();
    return first;
}

std::size_t BandEnergies::add(const PointRaySet& set)
{
    if (submitted_)
        throw std::logic_error("BandEnergies::add after submit, before clear");
    if (set.nodes.size() != set.problems() || set.absorption.size() != set.problems()
            || set.rays.size() != set.length.size())
        throw std::logic_error("BandEnergies::add: the set is not planned, or its rays do not match its lengths");
    sets_.push_back({nullptr, &set});
    const std::size_t first = problem_count_;
    problem_count_ += set.problems();
    return first;
}

void BandEnergies::submit()
{
    if (submitted_)
        return;
    std::size_t rays = 0, entries = 0, rule_floats = 0, wedge_count = 0, face_count = 0, panel_count = 0,
                reflected_count = 0, item_count = 0, gpu_problems = 0, sample_count = 0, probe_count = 0;
    for (const Entry& e : sets_)
    {
        const std::vector<double>& length = e.host ? e.host->length : e.point->length;
        const std::vector<uint32_t>& nodes = e.host ? e.host->nodes : e.point->nodes;
        rays += length.size();
        entries += length.size() * nodes.size();
        for (uint32_t m : nodes)
            rule_floats += 2 * m;
        if (e.point)
        {
            wedge_count += e.point->wedges.size();
            face_count += e.point->faces.size();
            panel_count += e.point->panels.size();
            reflected_count += e.point->reflected.size();
            item_count += e.point->face_rays.size() * e.point->problems() + e.point->probe_items.size();
            gpu_problems += e.point->problems();
            sample_count += e.point->samples.size();
            probe_count += e.point->probe_problems.size();
        }
    }
    if (entries > UINT32_MAX / 2)
        throw std::runtime_error("BandEnergies: too many rays in one batch");
    ensure<BandProblem>(device_, problems_, problem_count_);
    ensure<unsigned>(device_, which_pairs_, problem_count_);
    ensure<unsigned>(device_, which_quadrature_, problem_count_);
    ensure<unsigned>(device_, which_rays_, gpu_problems);
    ensure<float>(device_, amplitude_, 2 * entries);
    ensure<float>(device_, phasor_, 2 * entries);
    ensure<float>(device_, trig_, 2 * entries);
    ensure<float>(device_, length_, 2 * rays);
    ensure<float>(device_, rel_, rays);
    ensure<float>(device_, rule_, rule_floats);
    ensure<float>(device_, energy_, problem_count_);
    ensure<GpuRay>(device_, rays_, rays);
    ensure<GpuWedge>(device_, wedges_, wedge_count);
    ensure<GpuFace>(device_, faces_, face_count);
    ensure<GpuPanel>(device_, panels_, panel_count);
    ensure<float>(device_, reflected_, reflected_count);
    ensure<GpuFaceItem>(device_, face_items_, item_count);
    ensure<GpuFaceSample>(device_, samples_, sample_count);
    ensure<unsigned>(device_, probe_problems_, probe_count);
    ensure<float>(device_, probe_out_, 2 * probe_count);
    ensure<float>(device_, launch_, 3 * rays);

    BandProblem* records = problems_.as<BandProblem>();
    unsigned* pairs = which_pairs_.as<unsigned>();
    unsigned* integrals = which_quadrature_.as<unsigned>();
    unsigned* gpu_amplitudes = which_rays_.as<unsigned>();
    float* amplitude = amplitude_.as<float>();
    float* length = length_.as<float>();
    float* rel = rel_.as<float>();
    float* rule = rule_.as<float>();
    GpuRay* ray_records = rays_.as<GpuRay>();
    GpuWedge* wedges = wedges_.as<GpuWedge>();
    GpuFace* faces = faces_.as<GpuFace>();
    GpuPanel* panels = panels_.as<GpuPanel>();
    float* reflected = reflected_.as<float>();
    GpuFaceItem* items = face_items_.as<GpuFaceItem>();
    GpuFaceSample* samples = samples_.as<GpuFaceSample>();
    unsigned* probe_problems = probe_problems_.as<unsigned>();
    float* launch = launch_.as<float>();
    std::size_t q = 0, ray = 0, entry = 0, rule_at = 0, pair_count = 0, wedge_at = 0, face_at = 0, panel_at = 0,
                reflected_at = 0, item_at = 0, gpu_at = 0, sample_at = 0, probe_at = 0;
    quadrature_count_ = 0;
    for (const Entry& e : sets_)
    {
        const std::vector<double>& L = e.host ? e.host->length : e.point->length;
        const std::vector<double>& k1 = e.host ? e.host->k1 : e.point->k1;
        const std::vector<double>& k2 = e.host ? e.host->k2 : e.point->k2;
        const std::vector<uint32_t>& nodes = e.host ? e.host->nodes : e.point->nodes;
        const std::size_t n = L.size(), problems = k1.size();
        for (std::size_t i = 0; i < n; i++)
        {
            split(L[i], length[2 * (ray + i)], length[2 * (ray + i) + 1]);
            rel[ray + i] = static_cast<float>(L[i] - L[0]);
        }
        if (e.point && e.point->launch.size() == 3 * n)
            std::copy(e.point->launch.begin(), e.point->launch.end(), launch + 3 * ray);
        else
            std::fill(launch + 3 * ray, launch + 3 * (ray + n), 0.0f);
        if (e.host)
            std::copy(reinterpret_cast<const float*>(e.host->amplitude.data()),
                    reinterpret_cast<const float*>(e.host->amplitude.data()) + 2 * e.host->amplitude.size(),
                    amplitude + 2 * entry);
        else
        {
            // The set's records, their offsets moved to the batch's arrays.
            const PointRaySet& p = *e.point;
            for (std::size_t i = 0; i < n; i++)
            {
                GpuRay r = p.rays[i];
                r.term += static_cast<uint32_t>(r.kind == NK_KIND_FACE ? face_at : wedge_at);
                r.reflected += static_cast<uint32_t>(reflected_at);
                ray_records[ray + i] = r;
            }
            std::copy(p.wedges.begin(), p.wedges.end(), wedges + wedge_at);
            for (std::size_t f = 0; f < p.faces.size(); f++)
            {
                GpuFace face = p.faces[f];
                for (GpuAperture& end : face.end)
                    end.panel += static_cast<uint32_t>(panel_at);
                if (face.cascade != NK_NO_CASCADE)
                    face.cascade += static_cast<uint32_t>(wedge_at);
                if (face.ends != NK_NO_CASCADE)
                    face.ends += static_cast<uint32_t>(wedge_at);
                faces[face_at + f] = face;
            }
            // The sampled faces: their records, the probes' problems (set-local,
            // as face_ratios reads them) and the integrals at the probes.
            for (const GpuFaceSample& s : p.samples)
            {
                GpuFaceSample record = s;
                record.first = static_cast<uint32_t>(q);
                record.probe += static_cast<uint32_t>(probe_at);
                samples[sample_at++] = record;
            }
            std::copy(p.probe_problems.begin(), p.probe_problems.end(), probe_problems + probe_at);
            for (const GpuFaceItem& probe : p.probe_items)
                items[item_at++] = {static_cast<uint32_t>(q + probe.problem), probe.ray,
                        static_cast<uint32_t>(face_at + probe.face), static_cast<uint32_t>(probe_at + probe.probe)};
            probe_at += p.probe_problems.size();
            std::copy(p.panels.begin(), p.panels.end(), panels + panel_at);
            std::copy(p.reflected.begin(), p.reflected.end(), reflected + reflected_at);
        }
        for (std::size_t m = 0; m < problems; m++, q++)
        {
            const double kmid = 0.5 * (k1[m] + k2[m]), h = 0.5 * (k2[m] - k1[m]);
            BandProblem& r = records[q];
            r.rays = static_cast<uint32_t>(ray);
            r.first = static_cast<uint32_t>(entry + m * n);
            r.n = static_cast<uint32_t>(n);
            r.nodes = nodes[m];
            r.rule = static_cast<uint32_t>(rule_at / 2);
            split(kmid, r.k_hi, r.k_lo);
            split(h, r.h_hi, r.h_lo);
            r.h = static_cast<float>(h);
            r.reach = h > 0.0 ? static_cast<float>(30.0 / h) : FLT_MAX;
            r.centre = n ? static_cast<float>(0.5 * (L[n - 1] - L[0])) : 0.0f;
            r.complex_reflections = e.point && e.point->complex_reflections ? 1u : 0u;
            r.slot = static_cast<uint32_t>(r.complex_reflections ? m : m / 3);
            r.absorption = e.point ? e.point->absorption[m] : 0.0f;
            r.sum = e.point ? e.point->sum : NK_SUM_COHERENT;
            r.cone = e.point ? e.point->cone : 1.0f;
            if (r.nodes)
            {
                const std::vector<double>& rule_nodes = gauss_legendre(r.nodes);
                for (double v : rule_nodes)
                    rule[rule_at++] = static_cast<float>(v);
                integrals[quadrature_count_++] = static_cast<unsigned>(q);
            }
            else
                pairs[pair_count++] = static_cast<unsigned>(q);
            if (e.point)
            {
                gpu_amplitudes[gpu_at++] = static_cast<unsigned>(q);
                for (uint32_t face_ray : e.point->face_rays)
                    items[item_at++] = {static_cast<uint32_t>(q), face_ray,
                            static_cast<uint32_t>(face_at + e.point->rays[face_ray].term), 0};
            }
        }
        if (e.point)
        {
            wedge_at += e.point->wedges.size();
            face_at += e.point->faces.size();
            panel_at += e.point->panels.size();
            reflected_at += e.point->reflected.size();
        }
        ray += n;
        entry += n * problems;
    }

    batch_ = device_.batch();
    batch_.dispatch(ray_amplitudes_,
            {Arg(problems_), Arg(which_rays_), Arg(rays_), Arg(wedges_), Arg(reflected_), Arg(amplitude_), Arg(faces_)},
            gpu_at, 256);
    batch_.dispatch(face_amplitudes_,
            {Arg(face_items_), Arg(problems_), Arg(rays_), Arg(faces_), Arg(panels_), Arg(reflected_), Arg(amplitude_),
                    Arg(probe_out_), Arg(wedges_)},
            item_at, 32);
    const uint32_t sampled = static_cast<uint32_t>(sample_at);
    batch_.dispatch(face_ratios_,
            {Arg(samples_), Arg(problems_), Arg(probe_problems_), Arg(probe_out_), Arg(amplitude_), Arg::value(sampled)},
            (sample_at + 63) / 64, 64);
    batch_.dispatch(phasors_, {Arg(problems_), Arg(amplitude_), Arg(length_), Arg(phasor_), Arg(trig_)},
            problem_count_, 256);
    batch_.dispatch(pair_sum_,
            {Arg(problems_), Arg(which_pairs_), Arg(phasor_), Arg(trig_), Arg(rel_), Arg(energy_), Arg(launch_)},
            pair_count, 256);
    batch_.dispatch(quadrature_,
            {Arg(problems_), Arg(which_quadrature_), Arg(phasor_), Arg(rule_), Arg(rel_), Arg(energy_)},
            quadrature_count_, 64);
    batch_.commit();
    submitted_ = true;
}

const std::vector<double>& BandEnergies::wait()
{
    submit();
    batch_.wait();
    gpu_seconds_ = batch_.gpu_seconds();
    const float* energy = energy_.as<float>();
    energy_out_.assign(energy, energy + problem_count_);
    return energy_out_;
}

void BandEnergies::clear()
{
    if (submitted_)
        batch_.wait();
    submitted_ = false;
    sets_.clear();
    problem_count_ = 0;
    quadrature_count_ = 0;
}

} // namespace niceops::gpu
