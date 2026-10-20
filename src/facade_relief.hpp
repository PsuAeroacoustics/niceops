#pragma once

// A facade's relief as a grating (THEORY 6.5): a profile periodic along the
// facade and the same up it -- ribs, a sawtooth in plan, a container's
// corrugation -- that reflects a plane wave into its grating orders.  For a
// wave arriving in the profile's plane at incidence psi from the normal, the
// share of what the relief reflects that leaves in order m, at
//
//     sin(theta_m) = sin(theta_0) + m lambda / period,
//
// theta from the normal toward the facade's +along, theta_0 = -psi the
// specular direction, order 0 the specular reflection.  psi > 0 puts the
// source on the +along side, so the wave arrives travelling toward -along.  A
// wave arriving out of the profile's plane keeps the component of its
// direction up the facade and sees the profile at its horizontal wavenumber:
// it is the table's at the frequency times the cosine of its elevation, and
// each order leaves on the cone of the incident elevation.
//
// tools/profile_scattering.py --orders writes the table, by boundary
// elements: columns frequency_hz, incidence_deg, order, share and period_m.
// Shares are interpolated linearly in frequency and incidence, held beyond
// the table; an order that does not propagate at the frequency asked for is
// dropped and the others take its share, so they still sum to 1.

#include <string>
#include <utility>
#include <vector>

namespace niceops {

class FacadeRelief
{
public:
    // Throws std::runtime_error naming the file and the problem.
    static FacadeRelief read(const std::string& path);

    double period() const
    {
        return period_;
    }
    const std::string& source() const
    {
        return source_;
    }
    // Every order's share, as (order, share), at incidence sin(psi) and the
    // in-plane frequency (Hz), the wavelength there sound_speed / frequency,
    // into `out`, the specular order among them.
    void orders(double sin_psi, double frequency, double sound_speed, std::vector<std::pair<int, double>>& out) const;
    // The specular order's share.
    double specular(double sin_psi, double frequency, double sound_speed) const;
    // 1 less the specular share, averaged over incidence weighted by cos(psi)
    // (a 2-D diffuse field): the random-incidence scattering coefficient.
    double random_incidence(double frequency, double sound_speed) const;

private:
    std::string source_;
    double period_ = 0.0;
    std::vector<double> frequencies_;
    std::vector<double> sines_; // sin(psi), ascending
    int min_order_ = 0;
    int max_order_ = 0;
    // [frequency][incidence][order - min_order_]
    std::vector<double> shares_;
    double share(std::size_t f, std::size_t i, int m) const
    {
        const std::size_t orders = static_cast<std::size_t>(max_order_ - min_order_ + 1);
        return shares_[(f * sines_.size() + i) * orders + static_cast<std::size_t>(m - min_order_)];
    }
};

} // namespace niceops
