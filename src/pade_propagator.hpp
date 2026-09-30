#pragma once

// The parabolic equation's range propagator as a sum of resolvents: what lets
// pe.cpp take range steps of several wavelengths.
//
// One range step dr of the one-way wave equation multiplies each vertical mode
// by exp(i sigma (sqrt(1 + x) - 1)), sigma = k_a dr, x the mode's eigenvalue of
// Q = (d^2/dz^2 + k^2 - k_a^2) / k_a^2.  Collins' split-step Pade method
// (JASA 93, 1736 (1993)) replaces that function with a rational one,
//   r(x) = c_0 + sum_j c_j / (x - p_j),
// so a step is one tridiagonal solve per term, (Q - p_j) y_j = psi, whatever
// its length.  Crank-Nicolson is the one-term case, good for a tenth of a
// wavelength; seven terms are good for eight wavelengths.
//
// The coefficients are not Collins' Taylor-series Pade approximant about
// x = 0.  That one is exact at zero and degrades outward, and here the modes
// that matter are spread over an interval: elevation angles to some maximum,
// shifted by however far refraction takes k(z) from k_a.  r is instead the
// rational function of type (n, n) that interpolates the exact propagator at
// 2n + 1 Chebyshev nodes on that interval, which is near the best uniform
// approximation over it.  Off the interval nothing is asked of its accuracy,
// only that it not amplify: the evanescent modes, down to the grid's most
// negative eigenvalue, and the absorbing layer's, whose eigenvalues carry a
// positive imaginary part.  Nothing in the construction guarantees that, so it
// is measured -- |r| over the whole spectrum the caller names, the poles all
// below the real axis -- and a fit that fails either test is rejected.

#include <complex>
#include <vector>

namespace niceops {

struct PadePropagator
{
    std::complex<double> constant; // c_0, the value at infinity
    std::vector<std::complex<double>> poles; // p_j
    std::vector<std::complex<double>> residues; // c_j
    // Measured on dense samples: the largest |r - exact| over the accurate
    // interval, and the largest |r| over the whole spectrum.
    double error = 0.0;
    double largest_modulus = 0.0;

    std::complex<double> operator()(std::complex<double> x) const;
};

// exp(i sigma (sqrt(1 + x) - 1)), on the branch that decays for evanescent
// modes (x < -1) and for absorption (Im x > 0): the cut runs down from -1.
std::complex<double> exact_propagator(std::complex<double> x, double sigma);

// (1 + x)^(-1/4) exp(i sigma (sqrt(1 + x) - 1)) on the same branch: the
// self-starter's operator (Collins, JASA 92, 2069 (1992)).  A point source's
// field one step out is this applied to a delta function at its height --
// each mode launched with the amplitude the Hankel function gives it,
// k_n^(-1/2), rather than the one a Gaussian assumes.
std::complex<double> starter_propagator(std::complex<double> x, double sigma);

struct PadeRequest
{
    double sigma = 0.0; // k_a dr
    // Fit starter_propagator rather than exact_propagator.
    bool starter = false;
    // The interval of real x on which r must match the exact propagator.
    double accurate_low = 0.0;
    double accurate_high = 0.0;
    // The rest of the spectrum, where r must only not amplify: real x down to
    // `spectrum_low`, and imaginary parts from 0 up to `spectrum_imag`.
    double spectrum_low = 0.0;
    double spectrum_imag = 0.0;
    double tolerance = 1e-7; // on |r - exact| over the accurate interval
    int max_terms = 16;
};

// The propagator with the fewest terms that meets the request.  Throws
// std::runtime_error if none up to max_terms does: the step is too long for
// the interval, and the caller should shorten it.
PadePropagator pade_propagator(const PadeRequest& request);

} // namespace niceops
