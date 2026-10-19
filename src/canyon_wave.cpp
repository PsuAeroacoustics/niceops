#include "canyon_wave.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <deque>
#include <stdexcept>

#include "gpu/gpu.hpp"
#include "gpu/kernels/canyon_wave_types.h"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// In-place radix-2 transform, x_j <- sum_n x_n exp(+2 pi i n j / N); N a power
// of 2.  The +i sign is the time transform for fields varying as exp(-i w t)
// (P(w) = integral p(t) exp(i w t) dt), the convention the ground model uses.
void transform(std::vector<complex>& x)
{
    const std::size_t n = x.size();
    for (std::size_t i = 1, j = 0; i < n; i++)
    {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(x[i], x[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1)
    {
        const double angle = 2.0 * pi / static_cast<double>(len);
        const complex w_len(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len)
        {
            complex w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; k++)
            {
                const complex u = x[i + k];
                const complex v = x[i + k + len / 2] * w;
                x[i + k] = u + v;
                x[i + k + len / 2] = u - v;
                w *= w_len;
            }
        }
    }
}

// The grid of one run: nx columns by nz rows of cells h wide, pressure at the
// cell centers, horizontal velocity on the faces between columns and vertical
// velocity on the faces between rows.  Row j is active over columns
// [lo[j], hi[j]); every face on the edge of an active run, and the floor of
// row 0 and the top of the last row, is rigid (its velocity stays zero).
struct Grid
{
    std::size_t nx = 0;
    std::size_t nz = 0;
    double h = 0.0;
    std::vector<std::size_t> lo, hi;
    std::size_t pml = 0;          // PML cells at the left, right and top
    std::size_t source_i = 0, source_j = 0;
    std::size_t record_row = 0;   // the face row recorded
    std::size_t record_lo = 0, record_hi = 0;
    // The region whose energy decides when to stop: rows below `energy_top`
    // outside the PML.
    std::size_t energy_top = 0;
    // Rows below `wall_rows` have impedance walls at the faces lo[j] (wall A)
    // and hi[j] (wall B), of normalized impedance zeta_a and zeta_b; 0 is
    // rigid.
    std::size_t wall_rows = 0;
    double zeta_a = 0.0;
    double zeta_b = 0.0;
};

struct Recording
{
    std::vector<float> velocity; // [step * points + point]
    std::vector<double> source;  // q at each step
    double dt = 0.0;
    std::size_t steps = 0;
    double residual_db = 0.0;
    std::size_t cells = 0;
};

Recording run(const Grid& g, const CanyonWaveSetup& s)
{
    const double c = s.sound_speed;
    const double rho = s.density;
    const double h = g.h;
    const double dt = 0.9 * h / (c * std::sqrt(2.0));
    const std::size_t nx = g.nx, nz = g.nz;

    // The PML: a quadratic profile of damping over `pml` cells, reflecting
    // 1e-5 at normal incidence in the continuum, the same for every
    // frequency.  Coefficients at the cell centers (pressure) and the faces
    // (velocity) of each column and row.
    const double depth = static_cast<double>(g.pml) * h;
    const double sigma_max = g.pml > 0 ? 3.0 * c * std::log(1e5) / (2.0 * depth) : 0.0;
    auto sigma = [&](double d) { return d <= 0.0 ? 0.0 : sigma_max * (d / depth) * (d / depth); };
    auto coefficients = [&](double sig, double& a, double& b) {
        a = (1.0 - 0.5 * sig * dt) / (1.0 + 0.5 * sig * dt);
        b = dt / (1.0 + 0.5 * sig * dt);
    };
    const double left = static_cast<double>(g.pml) * h;
    const double right = static_cast<double>(nx - g.pml) * h;
    const double top = static_cast<double>(nz - g.pml) * h;
    std::vector<double> ay_c(nx), by_c(nx), ay_f(nx + 1), by_f(nx + 1), az_c(nz), bz_c(nz), az_f(nz + 1),
            bz_f(nz + 1);
    for (std::size_t i = 0; i < nx; i++)
    {
        const double y = (static_cast<double>(i) + 0.5) * h;
        coefficients(sigma(std::max(left - y, y - right)), ay_c[i], by_c[i]);
    }
    for (std::size_t i = 0; i <= nx; i++)
    {
        const double y = static_cast<double>(i) * h;
        coefficients(sigma(std::max(left - y, y - right)), ay_f[i], by_f[i]);
    }
    for (std::size_t j = 0; j < nz; j++)
        coefficients(sigma((static_cast<double>(j) + 0.5) * h - top), az_c[j], bz_c[j]);
    for (std::size_t j = 0; j <= nz; j++)
        coefficients(sigma(static_cast<double>(j) * h - top), az_f[j], bz_f[j]);

    std::vector<double> vy((nx + 1) * nz, 0.0), vz(nx * (nz + 1), 0.0), py(nx * nz, 0.0), pz(nx * nz, 0.0);
    const double inv_rho_h = 1.0 / (rho * h);
    const double rc2_h = rho * c * c / h;
    const double inject = dt * rho * c * c / (h * h);

    // The pulse: q(t) = (t - t0) / tau exp(-(t - t0)^2 / 2 tau^2), whose
    // spectrum, w exp(-w^2 tau^2 / 2), peaks at 0.4 of the highest frequency,
    // is 0.18 of its peak there, and has no mean, so no static pressure is
    // left behind.
    const double tau = 1.0 / (2.0 * pi * 0.4 * s.max_frequency);
    const double t0 = 6.0 * tau;
    auto pulse = [&](double t) {
        const double u = (t - t0) / tau;
        return u * std::exp(-0.5 * u * u);
    };

    Recording out;
    out.dt = dt;
    const std::size_t points = g.record_hi - g.record_lo;
    for (std::size_t j = 0; j < nz; j++)
        out.cells += g.hi[j] - g.lo[j];

    const std::size_t max_steps = static_cast<std::size_t>(std::ceil(s.max_duration / dt));
    // Long enough for the pulse to cross the grid before the decay is judged.
    const double span = std::hypot(static_cast<double>(nx), static_cast<double>(nz)) * h;
    const std::size_t min_steps = static_cast<std::size_t>(std::ceil((2.0 * t0 + 2.0 * span / c) / dt));
    const double stop_ratio = std::pow(10.0, -0.1 * s.decay_db);
    double peak_energy = 0.0;
    double last_energy = 0.0;

    // An impedance wall's face velocity, from the momentum equation over the
    // half cell between the wall and the cell's center, with the wall's
    // pressure Z v (v into the wall) averaged over the step, which makes it
    // stable for any Z >= 0:
    //   v+ (rho/dt + Z/h) = v- (rho/dt - Z/h) -/+ (2/h) p_cell.
    auto keep = [&](double zeta) {
        const double z = zeta * rho * c;
        return (rho / dt - z / h) / (rho / dt + z / h);
    };
    auto push = [&](double zeta) { return (2.0 / h) / (rho / dt + zeta * rho * c / h); };
    const double keep_a = keep(g.zeta_a), push_a = push(g.zeta_a);
    const double keep_b = keep(g.zeta_b), push_b = push(g.zeta_b);
    auto update_velocity = [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t j = row_begin; j < row_end; j++)
        {
            const std::size_t lo = g.lo[j], hi = g.hi[j];
            const double* pyr = &py[j * nx];
            const double* pzr = &pz[j * nx];
            double* vyr = &vy[j * (nx + 1)];
            for (std::size_t i = lo + 1; i < hi; i++)
                vyr[i] = ay_f[i] * vyr[i] - by_f[i] * inv_rho_h * ((pyr[i] + pzr[i]) - (pyr[i - 1] + pzr[i - 1]));
            if (j < g.wall_rows)
            {
                // Wall A on the left of cell lo, into it is -y; wall B on the
                // right of cell hi - 1, into it is +y.
                if (g.zeta_a > 0.0)
                    vyr[lo] = keep_a * vyr[lo] - push_a * (pyr[lo] + pzr[lo]);
                if (g.zeta_b > 0.0)
                    vyr[hi] = keep_b * vyr[hi] + push_b * (pyr[hi - 1] + pzr[hi - 1]);
            }
            if (j == 0)
                continue;
            const std::size_t vlo = std::max(lo, g.lo[j - 1]), vhi = std::min(hi, g.hi[j - 1]);
            const double* pyb = &py[(j - 1) * nx];
            const double* pzb = &pz[(j - 1) * nx];
            double* vzr = &vz[j * nx];
            const double a = az_f[j], b = bz_f[j] * inv_rho_h;
            for (std::size_t i = vlo; i < vhi; i++)
                vzr[i] = a * vzr[i] - b * ((pyr[i] + pzr[i]) - (pyb[i] + pzb[i]));
        }
    };
    auto update_pressure = [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t j = row_begin; j < row_end; j++)
        {
            const std::size_t lo = g.lo[j], hi = g.hi[j];
            double* pyr = &py[j * nx];
            double* pzr = &pz[j * nx];
            const double* vyr = &vy[j * (nx + 1)];
            const double* vzr = &vz[j * nx];
            const double* vzu = &vz[(j + 1) * nx];
            const double a = az_c[j], b = bz_c[j] * rc2_h;
            for (std::size_t i = lo; i < hi; i++)
            {
                pyr[i] = ay_c[i] * pyr[i] - by_c[i] * rc2_h * (vyr[i + 1] - vyr[i]);
                pzr[i] = a * pzr[i] - b * (vzu[i] - vzr[i]);
            }
        }
    };
    auto rows = [&](auto&& body) {
        if (s.threads <= 1)
            body(std::size_t {0}, nz);
        else
            parallel_chunks(nz, s.threads, body);
    };

    std::size_t n = 0;
    for (; n < max_steps; n++)
    {
        rows(update_velocity);
        const double* row = &vz[g.record_row * nx + g.record_lo];
        out.velocity.insert(out.velocity.end(), row, row + points);
        const double q = pulse((static_cast<double>(n) + 0.5) * dt);
        out.source.push_back(q);
        rows(update_pressure);
        py[g.source_j * nx + g.source_i] += inject * q;

        if (n % 32 == 31)
        {
            double energy = 0.0;
            for (std::size_t j = 0; j < g.energy_top; j++)
            {
                const std::size_t lo = std::max(g.lo[j], g.pml), hi = std::min(g.hi[j], nx - g.pml);
                for (std::size_t i = lo; i < hi; i++)
                {
                    const double p = py[j * nx + i] + pz[j * nx + i];
                    energy += p * p;
                }
            }
            peak_energy = std::max(peak_energy, energy);
            last_energy = energy;
            if (n >= min_steps && energy <= stop_ratio * peak_energy)
            {
                n++;
                break;
            }
        }
    }
    out.steps = n;
    out.residual_db = peak_energy > 0.0 ? 10.0 * std::log10(std::max(last_energy, 1e-300) / peak_energy) : 0.0;
    return out;
}

// run() on the GPU, in float: the same grid, coefficients, pulse and stopping
// rule, two dispatches a step and the energy every 32 steps as partial sums
// per group.  The steps go over 1024 at a time, two such batches in flight; a
// batch past the stop runs on harmlessly (the field never depends on when the
// run stops), and the recording is cut at the step run() would stop at.
Recording run_gpu(const Grid& g, const CanyonWaveSetup& s)
{
    gpu::Device* device = gpu::Device::get();
    if (!device)
        throw std::runtime_error("canyon wave: no GPU (" + std::string(gpu::Device::backend()) + " backend)");
    const double c = s.sound_speed;
    const double rho = s.density;
    const double h = g.h;
    const double dt = 0.9 * h / (c * std::sqrt(2.0));
    const std::size_t nx = g.nx, nz = g.nz;

    const double depth = static_cast<double>(g.pml) * h;
    const double sigma_max = g.pml > 0 ? 3.0 * c * std::log(1e5) / (2.0 * depth) : 0.0;
    auto sigma = [&](double d) { return d <= 0.0 ? 0.0 : sigma_max * (d / depth) * (d / depth); };
    const double left = static_cast<double>(g.pml) * h;
    const double right = static_cast<double>(nx - g.pml) * h;
    const double top = static_cast<double>(nz - g.pml) * h;
    // ay_c, by_c, ay_f, by_f, az_c, bz_c, az_f, bz_f, as canyon_wave.kern reads them.
    std::vector<float> coefficient;
    auto put = [&](double sig) {
        coefficient.push_back(static_cast<float>((1.0 - 0.5 * sig * dt) / (1.0 + 0.5 * sig * dt)));
        return static_cast<float>(dt / (1.0 + 0.5 * sig * dt));
    };
    std::vector<float> b;
    for (std::size_t i = 0; i < nx; i++)
        b.push_back(put(sigma(std::max(left - (static_cast<double>(i) + 0.5) * h, (static_cast<double>(i) + 0.5) * h - right))));
    coefficient.insert(coefficient.end(), b.begin(), b.end());
    b.clear();
    for (std::size_t i = 0; i <= nx; i++)
        b.push_back(put(sigma(std::max(left - static_cast<double>(i) * h, static_cast<double>(i) * h - right))));
    coefficient.insert(coefficient.end(), b.begin(), b.end());
    b.clear();
    for (std::size_t j = 0; j < nz; j++)
        b.push_back(put(sigma((static_cast<double>(j) + 0.5) * h - top)));
    coefficient.insert(coefficient.end(), b.begin(), b.end());
    b.clear();
    for (std::size_t j = 0; j <= nz; j++)
        b.push_back(put(sigma(static_cast<double>(j) * h - top)));
    coefficient.insert(coefficient.end(), b.begin(), b.end());

    CanyonGridGpu grid {};
    grid.nx = static_cast<uint32_t>(nx);
    grid.nz = static_cast<uint32_t>(nz);
    grid.pml = static_cast<uint32_t>(g.pml);
    grid.wall_rows = static_cast<uint32_t>(g.wall_rows);
    grid.record_row = static_cast<uint32_t>(g.record_row);
    grid.record_lo = static_cast<uint32_t>(g.record_lo);
    grid.record_hi = static_cast<uint32_t>(g.record_hi);
    grid.source_i = static_cast<uint32_t>(g.source_i);
    grid.source_j = static_cast<uint32_t>(g.source_j);
    grid.energy_top = static_cast<uint32_t>(g.energy_top);
    grid.wall_a = g.zeta_a > 0.0;
    grid.wall_b = g.zeta_b > 0.0;
    grid.inv_rho_h = static_cast<float>(1.0 / (rho * h));
    grid.rc2_h = static_cast<float>(rho * c * c / h);
    grid.inject = static_cast<float>(dt * rho * c * c / (h * h));
    auto keep = [&](double zeta) { return (rho / dt - zeta * rho * c / h) / (rho / dt + zeta * rho * c / h); };
    auto push = [&](double zeta) { return (2.0 / h) / (rho / dt + zeta * rho * c / h); };
    grid.keep_a = static_cast<float>(keep(g.zeta_a));
    grid.push_a = static_cast<float>(push(g.zeta_a));
    grid.keep_b = static_cast<float>(keep(g.zeta_b));
    grid.push_b = static_cast<float>(push(g.zeta_b));

    const double tau = 1.0 / (2.0 * pi * 0.4 * s.max_frequency);
    const double t0 = 6.0 * tau;
    const std::size_t max_steps = static_cast<std::size_t>(std::ceil(s.max_duration / dt));
    const double span = std::hypot(static_cast<double>(nx), static_cast<double>(nz)) * h;
    const std::size_t min_steps = static_cast<std::size_t>(std::ceil((2.0 * t0 + 2.0 * span / c) / dt));
    const double stop_ratio = std::pow(10.0, -0.1 * s.decay_db);
    const std::size_t points = g.record_hi - g.record_lo;

    Recording out;
    out.dt = dt;
    for (std::size_t j = 0; j < nz; j++)
        out.cells += g.hi[j] - g.lo[j];
    out.source.resize(max_steps);
    std::vector<float> q(max_steps);
    for (std::size_t n = 0; n < max_steps; n++)
    {
        const double u = ((static_cast<double>(n) + 0.5) * dt - t0) / tau;
        out.source[n] = u * std::exp(-0.5 * u * u);
        q[n] = static_cast<float>(out.source[n]);
    }

    const std::vector<std::string> program {"canyon_wave_types.h", "canyon_wave.kern"};
    const gpu::Kernel step_kernel = device->kernel(program, "canyon_step");
    // canyon_step's tiles, 32 by 8 cells.
    const std::size_t tiles = ((nx + 31) / 32) * ((nz + 7) / 8), checks = max_steps / 32 + 1;
    gpu::Buffer grid_buffer = device->buffer(sizeof(grid));
    std::memcpy(grid_buffer.data(), &grid, sizeof(grid));
    gpu::Buffer extent = device->buffer(2 * nz * sizeof(uint32_t));
    for (std::size_t j = 0; j < nz; j++)
    {
        extent.as<uint32_t>()[j] = static_cast<uint32_t>(g.lo[j]);
        extent.as<uint32_t>()[nz + j] = static_cast<uint32_t>(g.hi[j]);
    }
    gpu::Buffer coefficients = device->buffer(coefficient.size() * sizeof(float));
    std::copy(coefficient.begin(), coefficient.end(), coefficients.as<float>());
    // Every field twice, old and new, swapped each step.
    gpu::Buffer py[2], pz[2], vy[2], vz[2];
    for (int k = 0; k < 2; k++)
    {
        py[k] = device->buffer(nx * nz * sizeof(float));
        pz[k] = device->buffer(nx * nz * sizeof(float));
        vy[k] = device->buffer((nx + 1) * nz * sizeof(float));
        vz[k] = device->buffer(nx * (nz + 1) * sizeof(float));
    }
    gpu::Buffer pulse = device->buffer(max_steps * sizeof(float));
    std::copy(q.begin(), q.end(), pulse.as<float>());
    gpu::Buffer recording = device->buffer(max_steps * points * sizeof(float));
    gpu::Buffer partial = device->buffer(checks * tiles * sizeof(float));

    const std::size_t chunk = 1024;
    struct InFlight
    {
        gpu::Batch batch;
        std::size_t begin, end;
    };
    std::deque<InFlight> flight;
    double peak_energy = 0.0, last_energy = 0.0;
    std::size_t steps = max_steps;
    bool stopped = false;
    // The checks of a finished batch, in order, by run()'s rule.
    auto judge = [&](const InFlight& f) {
        for (std::size_t n = f.begin; n < f.end && !stopped; n++)
        {
            if (n % 32 != 31)
                continue;
            double e = 0.0;
            for (std::size_t k = 0; k < tiles; k++)
                e += partial.as<float>()[(n / 32) * tiles + k];
            peak_energy = std::max(peak_energy, e);
            last_energy = e;
            if (n >= min_steps && e <= stop_ratio * peak_energy)
            {
                steps = n + 1;
                stopped = true;
            }
        }
    };
    for (std::size_t begin = 0; begin < max_steps && !stopped; begin += chunk)
    {
        const std::size_t end = std::min(begin + chunk, max_steps);
        gpu::Batch batch = device->batch();
        for (std::size_t n = begin; n < end; n++)
        {
            const int in = static_cast<int>(n % 2), out_ = 1 - in;
            const uint32_t check = n % 32 == 31 ? static_cast<uint32_t>(n / 32) : 0xffffffffu;
            batch.dispatch(step_kernel,
                    {gpu::Arg(grid_buffer), gpu::Arg(extent), gpu::Arg(coefficients), gpu::Arg(py[in]),
                            gpu::Arg(pz[in]), gpu::Arg(vy[in]), gpu::Arg(vz[in]), gpu::Arg(py[out_]),
                            gpu::Arg(pz[out_]), gpu::Arg(vy[out_]), gpu::Arg(vz[out_]), gpu::Arg(recording),
                            gpu::Arg(pulse), gpu::Arg(partial), gpu::Arg::value(static_cast<uint32_t>(n)),
                            gpu::Arg::value(check)},
                    tiles, 256);
        }
        batch.commit();
        flight.push_back({batch, begin, end});
        if (flight.size() == 2)
        {
            flight.front().batch.wait();
            judge(flight.front());
            flight.pop_front();
        }
    }
    while (!flight.empty())
    {
        flight.front().batch.wait();
        judge(flight.front());
        flight.pop_front();
    }
    out.steps = steps;
    out.residual_db = peak_energy > 0.0 ? 10.0 * std::log10(std::max(last_energy, 1e-300) / peak_energy) : 0.0;
    out.source.resize(steps);
    out.velocity.assign(recording.as<float>(), recording.as<float>() + steps * points);
    return out;
}

// The opening's spectra, per unit source volume velocity, from a run.
ApertureSpectrum spectra(const Recording& r, const Grid& g, const CanyonWaveSetup& s)
{
    ApertureSpectrum a;
    a.cell = g.h;
    a.sound_speed = s.sound_speed;
    a.duration = static_cast<double>(r.steps) * r.dt;
    a.residual_db = r.residual_db;
    a.cells = r.cells;
    a.steps = r.steps;
    const std::size_t wanted = static_cast<std::size_t>(std::ceil(1.0 / (s.frequency_step * r.dt)));
    std::size_t size = 1;
    while (size < std::max(wanted, r.steps))
        size <<= 1;
    a.frequency_step = 1.0 / (static_cast<double>(size) * r.dt);
    a.bins = static_cast<std::size_t>(std::floor(s.max_frequency / a.frequency_step)) + 1;
    const std::size_t points = g.record_hi - g.record_lo;
    a.spectrum.assign(a.bins * points, complex(0.0, 0.0));

    std::vector<complex> buffer(size);
    std::fill(buffer.begin(), buffer.end(), complex(0.0, 0.0));
    for (std::size_t n = 0; n < r.steps; n++)
        buffer[n] = r.source[n];
    transform(buffer);
    const std::vector<complex> source(buffer.begin(), buffer.begin() + static_cast<std::ptrdiff_t>(a.bins));
    // One transform per recording point, each on its own, so on any number of
    // threads the same bits: the caller's (the table builder gives a canyon
    // on the CPU one, and one on the GPU its share of the run's).
    parallel_chunks(points, std::max(1u, s.threads),
            [&](std::size_t begin, std::size_t end) {
                std::vector<complex> line(size);
                for (std::size_t p = begin; p < end; p++)
                {
                    std::fill(line.begin(), line.end(), complex(0.0, 0.0));
                    for (std::size_t n = 0; n < r.steps; n++)
                        line[n] = r.velocity[n * points + p];
                    transform(line);
                    // Bin 0 has no source to divide by (the pulse has no mean).
                    for (std::size_t bin = 1; bin < a.bins; bin++)
                        a.spectrum[bin * points + p] = line[bin] / source[bin];
                }
            });
    return a;
}

} // namespace

double ApertureSpectrum::gain(std::size_t bin, double phi) const
{
    const double k = 2.0 * pi * static_cast<double>(bin) * frequency_step / sound_speed;
    const double kc = k * std::cos(phi);
    const std::size_t points = y.size();
    const complex* v = &spectrum[bin * points];
    // The points are evenly spaced, so the phase steps by a constant factor:
    // one complex multiply per point instead of a sine and a cosine.
    const complex step = points > 1 ? std::polar(1.0, -kc * (y[1] - y[0])) : complex(1.0, 0.0);
    complex phase = std::polar(1.0, -kc * y[0]);
    complex sum(0.0, 0.0);
    for (std::size_t p = 0; p < points; p++)
    {
        sum += weight[p] * v[p] * phase;
        phase *= step;
    }
    return 4.0 * std::norm(sum);
}

double ApertureSpectrum::gain_at(double frequency, double phi) const
{
    const double x = frequency / frequency_step;
    const std::size_t b = std::min(static_cast<std::size_t>(x), bins - 2);
    const double w = std::clamp(x - static_cast<double>(b), 0.0, 1.0);
    return (1.0 - w) * gain(b, phi) + w * gain(b + 1, phi);
}

double ApertureSpectrum::open_gain(double frequency, double phi) const
{
    const double c = std::cos(2.0 * pi * frequency / sound_speed * source_height * std::sin(phi));
    return 4.0 * c * c;
}

double ApertureSpectrum::band_ratio(double f_lo, double f_hi, double phi) const
{
    const std::size_t first = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(f_lo / frequency_step)));
    const std::size_t last = std::min(bins - 1, static_cast<std::size_t>(std::floor(f_hi / frequency_step)));
    double canyon = 0.0, open = 0.0;
    if (first > last)
    {
        const std::size_t bin = std::clamp<std::size_t>(
                static_cast<std::size_t>(std::lround(0.5 * (f_lo + f_hi) / frequency_step)), 1, bins - 1);
        canyon = gain(bin, phi);
        open = open_gain(static_cast<double>(bin) * frequency_step, phi);
    }
    else
    {
        for (std::size_t bin = first; bin <= last; bin++)
        {
            canyon += gain(bin, phi);
            open += open_gain(static_cast<double>(bin) * frequency_step, phi);
        }
    }
    return canyon / open;
}

ApertureSpectrum solve_canyon(const CanyonWaveSetup& s)
{
    if (!(s.width > 0.0 && s.depth > 0.0 && s.offset > 0.0 && s.offset < s.width && s.microphone_height >= 0.0
                && s.microphone_height < s.depth))
        throw std::invalid_argument("canyon wave: the microphone must stand inside a canyon of positive width "
                                    "and depth, below the roofs");
    if (!(s.wall_absorption >= 0.0 && s.wall_absorption < 1.0)
            || !(std::isnan(s.wall_absorption_b) || (s.wall_absorption_b >= 0.0 && s.wall_absorption_b < 1.0)))
        throw std::invalid_argument("canyon wave: wall absorption must lie within [0, 1)");
    if (!(s.max_frequency > 0.0 && s.points_per_wavelength >= 6.0 && s.pml_cells >= 8))
        throw std::invalid_argument("canyon wave: give a positive frequency, at least 6 points per wavelength "
                                    "and 8 PML cells");
    const double target = s.sound_speed / s.max_frequency / s.points_per_wavelength;
    const std::size_t n_w = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(s.width / target)));
    Grid g;
    g.h = s.width / static_cast<double>(n_w);
    const std::size_t n_h = std::max<std::size_t>(1, static_cast<std::size_t>(std::lround(s.depth / g.h)));
    const std::size_t roof = static_cast<std::size_t>(std::ceil(s.roof_margin / g.h));
    const std::size_t sky = static_cast<std::size_t>(std::ceil(s.sky_margin / g.h));
    g.pml = static_cast<std::size_t>(s.pml_cells);
    g.nx = 2 * g.pml + 2 * roof + n_w;
    g.nz = n_h + sky + g.pml;
    const std::size_t wall_a = g.pml + roof;
    g.lo.resize(g.nz);
    g.hi.resize(g.nz);
    for (std::size_t j = 0; j < g.nz; j++)
    {
        g.lo[j] = j < n_h ? wall_a : 0;
        g.hi[j] = j < n_h ? wall_a + n_w : g.nx;
    }
    g.source_i = wall_a + std::min(n_w - 1, static_cast<std::size_t>(s.offset / g.h));
    g.source_j = std::min(n_h - 1, static_cast<std::size_t>(s.microphone_height / g.h));
    g.record_row = n_h;
    g.record_lo = wall_a;
    g.record_hi = wall_a + n_w;
    g.energy_top = g.nz - g.pml;
    auto zeta = [](double alpha) {
        const double r = std::sqrt(1.0 - alpha);
        return alpha > 0.0 ? (1.0 + r) / (1.0 - r) : 0.0;
    };
    g.zeta_a = zeta(s.wall_absorption);
    g.zeta_b = zeta(std::isnan(s.wall_absorption_b) ? s.wall_absorption : s.wall_absorption_b);
    g.wall_rows = n_h;

    const Recording r = s.gpu ? run_gpu(g, s) : run(g, s);
    ApertureSpectrum a = spectra(r, g, s);
    a.opening_height = static_cast<double>(n_h) * g.h;
    a.source_height = (static_cast<double>(g.source_j) + 0.5) * g.h;
    for (std::size_t i = g.record_lo; i < g.record_hi; i++)
    {
        a.y.push_back((static_cast<double>(i - wall_a) + 0.5) * g.h);
        a.weight.push_back(g.h);
    }
    return a;
}

ApertureSpectrum solve_open_ground(const CanyonWaveSetup& s, double recording_height, double line_length,
        double taper)
{
    if (!(recording_height > s.microphone_height && line_length > 0.0 && taper >= 0.0 && taper < 0.5))
        throw std::invalid_argument("open-ground wave: record above the source, along a line of positive length");
    Grid g;
    g.h = s.sound_speed / s.max_frequency / s.points_per_wavelength;
    const std::size_t line = static_cast<std::size_t>(std::ceil(line_length / g.h));
    const std::size_t margin = static_cast<std::size_t>(std::ceil(s.roof_margin / g.h));
    const std::size_t rows = static_cast<std::size_t>(std::ceil((recording_height + s.sky_margin) / g.h));
    g.pml = static_cast<std::size_t>(s.pml_cells);
    g.nx = 2 * g.pml + 2 * margin + line;
    g.nz = rows + g.pml;
    g.lo.assign(g.nz, 0);
    g.hi.assign(g.nz, g.nx);
    g.source_i = g.nx / 2;
    g.source_j = static_cast<std::size_t>(s.microphone_height / g.h);
    g.record_row = static_cast<std::size_t>(std::lround(recording_height / g.h));
    g.record_lo = g.nx / 2 - line / 2;
    g.record_hi = g.record_lo + line;
    g.energy_top = g.nz - g.pml;

    const Recording r = s.gpu ? run_gpu(g, s) : run(g, s);
    ApertureSpectrum a = spectra(r, g, s);
    a.opening_height = static_cast<double>(g.record_row) * g.h;
    a.source_height = (static_cast<double>(g.source_j) + 0.5) * g.h;
    const double ends = taper * static_cast<double>(line);
    for (std::size_t i = g.record_lo; i < g.record_hi; i++)
    {
        a.y.push_back((static_cast<double>(i) - static_cast<double>(g.source_i)) * g.h);
        // A Tukey window: a cosine taper over `taper` of the line at each end,
        // so the line's truncation does not ring.
        const double from_end = std::min(static_cast<double>(i - g.record_lo) + 0.5,
                static_cast<double>(g.record_hi - i) - 0.5);
        const double w = ends > 0.0 && from_end < ends ? 0.5 * (1.0 - std::cos(pi * from_end / ends)) : 1.0;
        a.weight.push_back(g.h * w);
    }
    return a;
}

double canyon_wave_ratio(const ApertureSpectrum& spectrum, double u_c, double u_z, double f_lo, double f_hi,
        double floor_frequency)
{
    const double perp = std::hypot(u_c, u_z);
    if (f_hi * perp < floor_frequency)
        return 1.0;
    return spectrum.band_ratio(f_lo * perp, f_hi * perp, std::atan2(u_z, u_c));
}

std::vector<double> canyon_wave_gains(const ApertureSpectrum& spectrum, double axis_bearing,
        const std::vector<double>& azimuths, const std::vector<double>& elevations, const std::vector<double>& bands,
        std::size_t band_count, double floor_frequency, unsigned threads)
{
    const std::size_t nb = bands.size();
    std::vector<double> gains(azimuths.size() * elevations.size() * nb, 0.0);
    const double a = axis_bearing * pi / 180.0;
    const double half = std::pow(2.0, 1.0 / 6.0);
    parallel_chunks(azimuths.size(), std::max(1u, threads), [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++)
            for (std::size_t j = 0; j < elevations.size(); j++)
            {
                const double e = elevations[j] * pi / 180.0;
                const double u_c = std::cos(e) * std::sin(azimuths[i] * pi / 180.0 - a);
                const double u_z = std::sin(e);
                for (std::size_t b = 0; b < std::min(band_count, nb); b++)
                {
                    const double ratio = canyon_wave_ratio(spectrum, u_c, u_z, bands[b] / half, bands[b] * half,
                            floor_frequency);
                    gains[(i * elevations.size() + j) * nb + b] = 10.0 * std::log10(ratio);
                }
            }
    });
    return gains;
}

} // namespace niceops
