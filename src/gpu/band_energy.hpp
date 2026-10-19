#pragma once

// band_averaged_energy (urban.hpp) for many ray sets at once, on the GPU.
// A problem is a set of rays sorted by length, their complex amplitudes at the
// middle of a band [k1, k2], and the band.  The GPU takes each one either pair
// by pair, as band_averaged_energy does, or as a Gauss-Legendre integral of
// |S(k)|^2 over the band, whichever costs less there: the integral wherever
// the band is narrow against the rays' spread of lengths, which is where the
// pair sum keeps every pair.  The integral takes every pair, where the pair sum
// leaves out those its sinc has damped below 1/30, so the two answers differ
// by what that cutoff leaves out (FIDELITY has how much).
//
// The amplitudes go over in float.  Each ray's phase k L is formed on the GPU
// in float-float arithmetic (the length and the wavenumber each as a pair of
// floats) and reduced modulo 2 pi before its sine and cosine are taken, so it
// carries about 3e-7 rad of error at any length; the sums are in float.

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "gpu/gpu.hpp"

namespace niceops::gpu {

// The problems on one set of rays, as a worker thread packs them: the rays'
// lengths, ascending, and for each problem its band and the rays' amplitudes.
struct RaySet
{
    std::vector<double> length;
    std::vector<std::complex<float>> amplitude; // [problem][ray]
    std::vector<double> k1, k2;                 // [problem]
    // Which way each problem goes: its Gauss-Legendre node count, or 0 for
    // the pair sum.  plan() fills it, from the lengths alone; it costs a pass
    // over the rays per problem, so it is the packing thread's to call.
    std::vector<uint32_t> nodes;

    // Adds a problem; its amplitudes are the next length.size() entries.
    void add(const std::complex<double>* a, double band_k1, double band_k2);
    void plan();
    std::size_t problems() const
    {
        return k1.size();
    }
};

// RaySet::plan's choice, for any set of lengths (ascending) and problems.
void plan_nodes(const std::vector<double>& length, const std::vector<double>& k1, const std::vector<double>& k2,
        std::vector<uint32_t>& nodes);

struct PointRaySet;

// One batch of ray sets on the GPU.  add() keeps a pointer to each set, which
// must outlive wait().
class BandEnergies
{
public:
    explicit BandEnergies(Device& device);

    // Adds a planned set and returns the index of its first problem; its
    // problems follow in order.  A RaySet brings its amplitudes; a
    // PointRaySet (point_rays.hpp) its rays, whose amplitudes the GPU forms.
    std::size_t add(const RaySet& set);
    std::size_t add(const PointRaySet& set);
    std::size_t size() const
    {
        return problem_count_;
    }
    // Copies the batch into the GPU's buffers and commits it; returns at once.
    void submit();
    // Waits for the batch and returns its energies, by problem index.
    const std::vector<double>& wait();
    // Forgets the batch, keeping the buffers for the next.
    void clear();
    // GPU time of the last wait(), s.
    double gpu_seconds() const
    {
        return gpu_seconds_;
    }
    // How many problems of the batch are taken as integrals.
    std::size_t integrals() const
    {
        return quadrature_count_;
    }

private:
    struct Entry
    {
        const RaySet* host = nullptr;
        const PointRaySet* point = nullptr;
    };
    Device& device_;
    Kernel phasors_, pair_sum_, quadrature_, ray_amplitudes_, face_amplitudes_, face_ratios_;
    std::vector<Entry> sets_;
    std::size_t problem_count_ = 0, quadrature_count_ = 0;
    Buffer problems_, which_pairs_, which_quadrature_, amplitude_, length_, phasor_, trig_, rel_, rule_, energy_;
    Buffer which_rays_, rays_, wedges_, faces_, panels_, reflected_, face_items_;
    Buffer samples_, probe_problems_, probe_out_, launch_;
    Batch batch_;
    bool submitted_ = false;
    std::vector<double> energy_out_;
    double gpu_seconds_ = 0.0;
};

} // namespace niceops::gpu
