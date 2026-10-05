#pragma once

// Band levels in fixed reception-time frames (--epnl_file, --band_frames).
//
// The observer loop evaluates one spectrum per trajectory state, in emission
// time, at the trajectory's own step (0.02 s on the shipped tracks).  A
// measured spectrum is a linear average over a fixed record in reception time
// -- 0.5 s for the one-third-octave records of 14 CFR 36 A36.3 -- so the
// prediction is put on the same footing: each state's energy is spread over
// the stretch of reception time it arrives in and summed into frames.
//
// State `it` is emitted at t and heard at t + r/c, r being the range whose rate
// of change gives reception time (reception_factor in sel.hpp; a refracted
// ray's c times its travel time).  Its emission step T arrives spread over
// T (1 + (dr/dt)/c), the same weight the SEL integral gives it, centered on the
// reception time.  That stretch is split between the frames it overlaps, in
// proportion, so a frame's energy is the integral of the received band
// pressure squared over it, and its level the linear (Leq) average:
//   L_frame = 10 lg( sum over states of 10^(L/10) x overlap / frame length ).
// Summing every frame's energy gives back exactly the whole-history exposure
// in each band (full_sound_exposure_level's integral), which is what the tests
// hold it to.
//
// Frame k covers [k F, (k + 1) F) on the trajectory's own clock, F the frame
// length, so the frames of different observers line up: a grid's frames share
// one time axis.  The first and last frames an observer hears are usually
// partly silent, and are averaged over the whole frame regardless -- silence
// is part of a record.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace niceops {

class BandFrames
{
public:
    BandFrames() = default;
    BandFrames(std::size_t bands, double frame_length) : bands_(bands), frame_length_(frame_length) {}

    // Spread the states' energy into frames.  For each state `it` < count,
    // `usable(it)` says whether it contributes (a culled state does not),
    // `reception(it)` is the reception time of its center, `duration(it)` the
    // stretch of reception time it covers, and `bands_of(it)` its band levels,
    // dB, the frames' `bands()` of them contiguous.  The frames span every
    // usable state's stretch; with none, there are none.
    template <typename Usable, typename Reception, typename Duration, typename BandsOf>
    void accumulate(std::size_t count, Usable&& usable, Reception&& reception, Duration&& duration,
            BandsOf&& bands_of)
    {
        double begin = std::numeric_limits<double>::infinity();
        double end = -std::numeric_limits<double>::infinity();
        // A stretch that is not a finite, positive interval (a NaN range, a
        // one-state track's zero step) carries nothing and must not size the
        // frames.
        auto stretch = [&](std::size_t it, double& lo, double& hi) {
            const double half = 0.5 * duration(it);
            lo = reception(it) - half;
            hi = reception(it) + half;
            return std::isfinite(lo) && std::isfinite(hi) && hi > lo;
        };
        for (std::size_t it = 0; it < count; it++)
        {
            double lo = 0.0, hi = 0.0;
            if (!usable(it) || !stretch(it, lo, hi))
                continue;
            begin = std::min(begin, lo);
            end = std::max(end, hi);
        }
        energy_.clear();
        frame_count_ = 0;
        if (!(begin < end))
            return;
        first_ = static_cast<std::int64_t>(std::floor(begin / frame_length_));
        // The last frame the stretch reaches into; an end exactly on a boundary
        // does not open the next frame.
        const std::int64_t last = std::max(first_,
                static_cast<std::int64_t>(std::ceil(end / frame_length_)) - 1);
        frame_count_ = static_cast<std::size_t>(last - first_ + 1);
        energy_.assign(frame_count_ * bands_, 0.0);
        // exp(L ln10/10) is 10^(L/10) to within rounding and cheaper than pow;
        // this runs once per band per state per observer (weighting.cpp does
        // the same).
        constexpr double ln10_over_10 = 0.23025850929940456840;
        for (std::size_t it = 0; it < count; it++)
        {
            double lo = 0.0, hi = 0.0;
            if (!usable(it) || !stretch(it, lo, hi))
                continue;
            const double* levels = bands_of(it);
            std::int64_t k = std::max(first_, static_cast<std::int64_t>(std::floor(lo / frame_length_)));
            for (; k - first_ < static_cast<std::int64_t>(frame_count_); k++)
            {
                const double frame_lo = static_cast<double>(k) * frame_length_;
                if (frame_lo >= hi)
                    break;
                const double overlap = std::min(hi, frame_lo + frame_length_) - std::max(lo, frame_lo);
                if (!(overlap > 0.0))
                    continue;
                double* out = energy_.data() + static_cast<std::size_t>(k - first_) * bands_;
                for (std::size_t band = 0; band < bands_; band++)
                    out[band] += std::exp(levels[band] * ln10_over_10) * overlap;
            }
        }
    }

    std::size_t bands() const { return bands_; }
    double frame_length() const { return frame_length_; }
    // Frame index of the first frame, on the trajectory's clock: frame k
    // covers [k F, (k + 1) F).
    std::int64_t first() const { return first_; }
    std::size_t count() const { return frame_count_; }
    // The center of frame `frame` (0 = first()), s.
    double center(std::size_t frame) const
    {
        return (static_cast<double>(first_ + static_cast<std::int64_t>(frame)) + 0.5) * frame_length_;
    }
    // Received energy, Pa^2 s / (20 uPa)^2, in one frame and band.
    double energy(std::size_t frame, std::size_t band) const { return energy_[frame * bands_ + band]; }
    // The frame's linear-average level, dB re 20 uPa; minus infinity for a
    // frame and band that received nothing.
    double level(std::size_t frame, std::size_t band) const
    {
        const double e = energy(frame, band);
        return e > 0.0 ? 10.0 * std::log10(e / frame_length_)
                       : -std::numeric_limits<double>::infinity();
    }

private:
    std::size_t bands_ = 0;
    double frame_length_ = 0.5;
    std::int64_t first_ = 0;
    std::size_t frame_count_ = 0;
    std::vector<double> energy_;
};

// One observer's frames, kept for --band_frames: its first frame index and
// count, and the levels, frame-major, as float.
struct ObserverFrames
{
    std::int64_t first = 0;
    std::size_t count = 0;
    std::vector<float> levels;
};

// What --band_frames writes besides the levels.
struct BandFramesFile
{
    std::string trajectory;
    std::string database;
    double frame_length = 0.5;
    std::vector<double> band_centers;  // Hz, the bands in `levels`
    // Microphone positions, in `length_unit`, and heights above the ground.
    std::vector<double> x, y, z, height;
    std::string length_unit;
    // A grid's axes, in `length_unit`, when the observers are one (observer =
    // iy * nx + ix); empty for a scattered set.
    std::vector<double> grid_x, grid_y;
};

// Write an (observer, frame, band) netCDF-4 file: float levels, NaN outside
// the frames each observer received, on one frame axis spanning them all.
// Throws (input_error) on any netCDF failure, closing included: HDF5 buffers
// its writes, so a full disk is reported only by nc_close.
void write_band_frames(const std::string& path, const BandFramesFile& meta,
        const std::vector<ObserverFrames>& frames);

} // namespace niceops
