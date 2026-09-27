#include "pe.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "util/constants.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

double speed_at(const PeSetup& s, double z)
{
    if (s.z.size() == 1 || z <= s.z.front())
        return s.c.front();
    if (z >= s.z.back())
        return s.c.back();
    const auto upper = std::upper_bound(s.z.begin(), s.z.end(), z);
    const std::size_t i = static_cast<std::size_t>(upper - s.z.begin()) - 1;
    return s.c[i] + (z - s.z[i]) / (s.z[i + 1] - s.z[i]) * (s.c[i + 1] - s.c[i]);
}

} // namespace

PeResult solve_pe(const PeSetup& s)
{
    if (s.z.empty() || s.z.size() != s.c.size())
        throw std::invalid_argument("solve_pe: the sound-speed profile is empty or uneven");
    if (!(s.frequency > 0.0) || !(s.max_range > 0.0) || s.receiver_heights.empty())
        throw std::invalid_argument("solve_pe: frequency, range and receivers are required");

    const double omega = 2.0 * pi * s.frequency;
    const double c_ref = speed_at(s, 0.0);
    const double k_a = omega / c_ref;
    const double wavelength = c_ref / s.frequency;
    const double dz = wavelength / s.points_per_wavelength;
    const double dr = wavelength / s.range_points_per_wavelength;

    // The domain: everything the source and receivers need, a clear margin,
    // then an absorbing layer that takes what leaves the top.  The margin is
    // generous on purpose.  In a shadow the field near the ground is a weak
    // remainder, while energy leaks upward past it; what the layer sends back
    // down lands on that remainder as a floor.  With 60 wavelengths of margin,
    // 25 m at 1 kHz, the exact residue series (tests/test_pe_shadow.cpp) found
    // the PE 3.8 dB high at -44 dB; with 250 m it agrees to hundredths.
    const double highest = std::max(s.source_height,
            *std::max_element(s.receiver_heights.begin(), s.receiver_heights.end()));
    const double physical_top = s.domain_height > 0.0
            ? s.domain_height
            : highest + std::max({60.0 * wavelength, 0.5 * highest, 250.0});
    const double layer = s.layer_wavelengths * wavelength;
    const double top = physical_top + layer;
    const std::size_t n = static_cast<std::size_t>(std::ceil(top / dz)); // psi_n = 0

    // epsilon_j = (k_j^2 - k_a^2) / k_a^2, complex with the absorption, and
    // with the absorbing layer: a quadratic ramp to Im k ~ k_a / 10 at the top.
    std::vector<complex> epsilon(n + 1);
    for (std::size_t j = 0; j <= n; j++)
    {
        const double z = static_cast<double>(j) * dz;
        const complex k(omega / speed_at(s, z), s.absorption);
        double ramp = 0.0;
        if (z > physical_top)
        {
            const double t = (z - physical_top) / layer;
            ramp = 0.2 * t * t;
        }
        epsilon[j] = (k * k - k_a * k_a) / (k_a * k_a) + complex(0.0, ramp);
    }

    // Rows j = 1 .. n-1 (psi_0 from the ground condition, psi_n = 0).
    // Q psi_j = a (psi_{j+1} - 2 psi_j + psi_{j-1}) + epsilon_j psi_j.
    const double a = 1.0 / (k_a * k_a * dz * dz);
    const complex i(0.0, 1.0);
    const complex alpha_lhs = 0.25 - i * k_a * dr / 4.0;
    const complex alpha_rhs = 0.25 + i * k_a * dr / 4.0;
    // Ground: (-3 psi_0 + 4 psi_1 - psi_2) / (2 dz) = -i k beta psi_0, so
    // psi_0 = (4 psi_1 - psi_2) / (3 - 2 i k beta dz) = g1 psi_1 + g2 psi_2.
    const complex denominator = 3.0 - 2.0 * i * k_a * s.admittance * dz;
    const complex g1 = 4.0 / denominator;
    const complex g2 = -1.0 / denominator;

    const std::size_t m = n - 1; // unknowns psi_1 .. psi_{n-1}
    auto diagonal = [&](complex alpha, std::size_t j) { return 1.0 + alpha * (-2.0 * a + epsilon[j]); };

    // The left side does not change with range: factor it once (Thomas).
    std::vector<complex> lower(m), main(m), upper(m);
    for (std::size_t row = 0; row < m; row++)
    {
        const std::size_t j = row + 1;
        lower[row] = alpha_lhs * a;
        main[row] = diagonal(alpha_lhs, j);
        upper[row] = alpha_lhs * a;
    }
    main[0] += alpha_lhs * a * g1;
    upper[0] += alpha_lhs * a * g2;
    std::vector<complex> c_prime(m), inv_denominator(m);
    {
        complex previous = 0.0;
        for (std::size_t row = 0; row < m; row++)
        {
            const complex d = main[row] - (row > 0 ? lower[row] * previous : complex(0.0));
            inv_denominator[row] = 1.0 / d;
            previous = upper[row] * inv_denominator[row];
            c_prime[row] = previous;
        }
    }
    const complex rhs_off = alpha_rhs * a;
    std::vector<complex> rhs_main(m);
    for (std::size_t row = 0; row < m; row++)
        rhs_main[row] = diagonal(alpha_rhs, row + 1);
    const complex rhs_main_0 = rhs_main[0] + alpha_rhs * a * g1;
    const complex rhs_upper_0 = alpha_rhs * a + alpha_rhs * a * g2;

    // Salomons' second-order Gaussian starter, with its image in the ground
    // at the normal-incidence reflection coefficient.  Normalised so that
    // p = psi exp(i k_a r) / sqrt(r) is a unit point source, exp(ikR)/R.
    std::vector<complex> psi(m);
    const complex reflection = (1.0 - s.admittance) / (1.0 + s.admittance);
    const complex scale = std::sqrt(i * k_a);
    auto starter = [&](double z) {
        const double u = k_a * k_a * z * z;
        return (1.3717 - 0.3701 * u) * std::exp(-u / 3.0);
    };
    for (std::size_t row = 0; row < m; row++)
    {
        const double z = static_cast<double>(row + 1) * dz;
        psi[row] = scale * (starter(z - s.source_height) + reflection * starter(z + s.source_height));
    }

    PeResult result;
    result.relative_level.assign(s.receiver_heights.size(), {});
    const std::size_t steps = static_cast<std::size_t>(std::ceil(s.max_range / dr));
    // psi at each receiver, linear between grid points; psi_0 from the ground
    // condition.
    auto at_receiver = [&](double zr) {
        const double u = zr / dz;
        const std::size_t j = std::min(static_cast<std::size_t>(u), m - 1);
        const double t = u - static_cast<double>(j);
        auto at = [&](std::size_t jj) {
            if (jj == 0)
                return g1 * psi[0] + g2 * psi[1];
            return jj <= m ? psi[jj - 1] : complex(0.0);
        };
        return (1.0 - t) * at(j) + t * at(j + 1);
    };
    // Results at exact multiples of the output step, whatever the range step
    // (which follows the wavelength), so every frequency lands on one grid:
    // psi interpolated between the two range steps either side.
    std::vector<complex> previous_value(s.receiver_heights.size());
    for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
        previous_value[k] = at_receiver(s.receiver_heights[k]);
    std::size_t next_output = 1;
    std::vector<complex> rhs(m);
    for (std::size_t step = 1; step <= steps; step++)
    {
        // Right side: [1 + alpha_rhs Q] psi, with psi_0 from the ground
        // condition and psi_n = 0.
        rhs[0] = rhs_main_0 * psi[0] + rhs_upper_0 * psi[1];
        for (std::size_t row = 1; row + 1 < m; row++)
            rhs[row] = rhs_off * (psi[row - 1] + psi[row + 1]) + rhs_main[row] * psi[row];
        rhs[m - 1] = rhs_off * psi[m - 2] + rhs_main[m - 1] * psi[m - 1];
        // Solve.
        complex previous = 0.0;
        for (std::size_t row = 0; row < m; row++)
        {
            previous = (rhs[row] - (row > 0 ? lower[row] * previous : complex(0.0))) * inv_denominator[row];
            psi[row] = previous;
        }
        for (std::size_t row = m - 1; row-- > 0;)
            psi[row] -= c_prime[row] * psi[row + 1];

        const double r = static_cast<double>(step) * dr;
        std::vector<complex> value(s.receiver_heights.size());
        for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
            value[k] = at_receiver(s.receiver_heights[k]);
        while (next_output * s.output_step <= r + 1e-9 && next_output * s.output_step <= s.max_range + 1e-9)
        {
            const double target = next_output * s.output_step;
            const double t = (target - (r - dr)) / dr;
            result.range.push_back(target);
            for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
            {
                const complex v = (1.0 - t) * previous_value[k] + t * value[k];
                const double zr = s.receiver_heights[k];
                const double slant = std::hypot(target, zr - s.source_height);
                result.relative_level[k].push_back(
                        20.0 * std::log10(std::abs(v) * slant / std::sqrt(target) + 1e-300));
            }
            next_output++;
        }
        previous_value = value;
    }
    return result;
}

PeResult solve_pe_band(PeSetup setup, double centre, int samples)
{
    // Evenly in log frequency across the band's edges, 2^(+-1/6).  Every
    // frequency reports on the same range grid (solve_pe), so the energies add
    // point by point.
    PeResult sum;
    for (int k = 0; k < samples; k++)
    {
        const double offset = samples == 1 ? 0.0 : (static_cast<double>(k) / (samples - 1) - 0.5) / 3.0;
        setup.frequency = centre * std::pow(2.0, offset);
        const auto one = solve_pe(setup);
        if (sum.range.empty())
        {
            sum.range = one.range;
            sum.relative_level.assign(one.relative_level.size(),
                    std::vector<double>(one.range.size(), 0.0));
        }
        for (std::size_t rcv = 0; rcv < sum.relative_level.size(); rcv++)
        {
            for (std::size_t q = 0; q < std::min(sum.range.size(), one.range.size()); q++)
                sum.relative_level[rcv][q] += std::pow(10.0, 0.1 * one.relative_level[rcv][q]) / samples;
        }
    }
    for (auto& row : sum.relative_level)
        for (double& v : row)
            v = 10.0 * std::log10(v + 1e-300);
    return sum;
}

} // namespace niceops
