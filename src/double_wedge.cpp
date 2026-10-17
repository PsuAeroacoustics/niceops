#include "double_wedge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

double wrap(double angle)
{
    return angle < 0.0 ? angle + 2.0 * pi : angle;
}

// From the first edge, the shared face along +x; from the second, along -x.
double angle1(double x, double z)
{
    return wrap(std::atan2(z, x));
}

double angle2(double x, double z)
{
    return wrap(std::atan2(z, -x));
}

// The moments' power series below a radian of phase: m_j, the integral of
// v^j exp(i k v) over [0, s], is s^(j+1) times the sum over n of
// (i phase)^n / (n! (n + j + 1)), phase = k s.  Its even terms are the real
// part and its odd ones the imaginary, each a polynomial in phase^2 whose
// coefficients, (-1)^(n/2) / (n! (n + j + 1)), are fixed; to n = 21 they leave
// under 1e-19 of the sum.
struct FilonSeries
{
    double c[3][22];
};

constexpr FilonSeries filon_series()
{
    FilonSeries series {};
    double factorial = 1.0;
    for (int n = 0; n < 22; n++)
    {
        if (n > 0)
            factorial *= n;
        for (int j = 0; j < 3; j++)
            series.c[j][n] = ((n / 2) % 2 == 0 ? 1.0 : -1.0) / (factorial * (n + j + 1));
    }
    return series;
}

constexpr FilonSeries filon_coefficients = filon_series();

// m_j / s^(j+1) for phase under a radian, by Horner in phase^2.
complex filon_moment_series(int j, double phase)
{
    const double u = phase * phase;
    const double* c = filon_coefficients.c[j];
    double even = c[20], odd = c[21];
    for (int n = 18; n >= 0; n -= 2)
    {
        even = even * u + c[n];
        odd = odd * u + c[n + 1];
    }
    return {even, phase * odd};
}

// The integral over [u0, u2] of g(u) exp(i k u), g the quadratic through
// (u0, g0), (u1, g1) and (u2, g2): Filon's rule, exact in the oscillation
// however many wavelengths the panel spans.  Its moments are taken from their
// recurrence where the panel spans a radian or more of phase, and from their
// power series below, where the recurrence cancels.
complex filon_panel(double k, double u0, double u1, double u2, complex g0, complex g1, complex g2)
{
    const double span = u2 - u0, d1 = u1 - u0;
    const complex c1 = (g1 - g0) / d1;
    const complex c2 = ((g2 - g0) / span - c1) / (u2 - u1);
    // m0, m1, m2: the integrals of v^m exp(i k v) over [0, span].
    const double phase = k * span;
    complex m0, m1, m2;
    if (std::abs(phase) < 1.0)
    {
        m0 = span * filon_moment_series(0, phase);
        m1 = span * span * filon_moment_series(1, phase);
        m2 = span * span * span * filon_moment_series(2, phase);
    }
    else
    {
        // Dividing by i k is multiplying by -i / k.
        const complex e = std::polar(1.0, phase), by_ik(0.0, -1.0 / k);
        m0 = (e - 1.0) * by_ik;
        m1 = (span * e - m0) * by_ik;
        m2 = (span * span * e - 2.0 * m1) * by_ik;
    }
    return std::polar(1.0, k * u0) * (g0 * m0 + c1 * m1 + c2 * (m2 - d1 * m1));
}

// One of the aperture integrand's waves at a point of the aperture: its path
// length (the phase per unit k), that length's rate along the aperture, and
// its amplitude.
struct Wave
{
    double length = 0.0;
    double rate = 0.0;
    complex amplitude;
};

// The source below the face's plane and before the first edge.  The aperture
// is the half-line x = 0, z > 0, out to 30 wavelengths or three Fresnel zones,
// its far half tapered.  Each of its points sends a line source's wave to the
// second wedge, whose field -- the geometrical-optics waves past it and its
// diffraction -- the receiver sums.
//
// The integrand is four waves, each a smooth amplitude times exp(i k l), l a
// path length that grows along the aperture: the second wedge's diffracted
// wave, the incident wave, and its reflections off the shared face and off
// the wedge's other face.  Each geometrical-optics wave reaches the receiver
// or not, and the coefficient jumps, only where the second wedge's shadow and
// reflection boundaries cross the aperture, so it is cut there and between the
// cuts each wave is integrated by Filon's rule in its own path length, the
// amplitude per unit length quadratic over each panel.  The panels grow
// geometrically away from the edge and the cuts, where the amplitude changes
// fastest.  Against the midpoint rule at a 192nd of a wavelength this is
// within 0.016 dB (0.0015 dB rms) over the shared faces of DigiCity and of a
// district; the midpoint rule at a 24th of a wavelength, which it replaced,
// was within 0.027 dB (0.014-0.018 dB rms), for 6 to 24 times the evaluations.
complex aperture(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    const double w = face.width;
    const double lambda = 2.0 * pi / k;
    const double s1 = std::hypot(sx, sz);
    const double extent = std::max(3.0 * std::sqrt(lambda * std::max(s1, w)), 30.0 * lambda);
    const double theta = angle2(rx - w, rz);
    const double s3 = std::hypot(rx - w, rz);
    // The second wedge's other face, for the reflection off it.
    const double fx = -std::cos(face.n2 * pi), fz = std::sin(face.n2 * pi);
    const complex launch = std::sqrt(k / (2.0 * pi)) * std::polar(1.0, -0.25 * pi);
    // Above the face's plane the receiver sees the first edge, and the
    // aperture's waves past the second edge hold that wedge's whole field --
    // the source's own wave where the wedge does not shadow it, and the edge's
    // diffracted wave with the face's reflection of it -- which the rays of
    // the first edge carry already.  The same aperture with the face going on
    // past the second edge is that field alone, in the aperture's own
    // approximation, and is taken back out: what remains is what the second
    // edge does to it.  So above the plane the incident and face-reflected
    // waves are counted less one, the unbounded field's.
    const bool above = rz > 0.0;

    // The cuts: where theta0, the aperture point's angle at the second wedge,
    // crosses one of hard_wedge_field's boundaries.  The taper's start is a
    // panel end too, though the amplitude changes smoothly there.
    std::vector<double> cuts {0.0, 0.5 * extent, extent};
    std::vector<double> features {0.0};
    for (double boundary : {theta - pi, pi - theta, (face.n2 - 1.0) * pi, (2.0 * face.n2 - 1.0) * pi - theta})
        if (boundary > 0.0 && boundary < 0.5 * pi && w * std::tan(boundary) < extent)
        {
            cuts.push_back(w * std::tan(boundary));
            features.push_back(cuts.back());
        }
    std::sort(cuts.begin(), cuts.end());

    // The waves at t, those of `count` (+1, -1 or 0) only.  The second
    // wedge's coefficient is wanted at one theta and many theta0.
    const WedgeSweep sweep(face.n2, theta);
    auto waves = [&](double t, const int* count, Wave* out) {
        const double taper = t < 0.5 * extent ? 1.0 : 0.5 * (1.0 + std::cos(pi * (t - 0.5 * extent) / (0.5 * extent)));
        const double in = std::hypot(sx, t - sz);
        const double r = std::hypot(w, t);
        const double obliquity = 0.5 * (-sx / in + w / r);
        const complex amplitude = taper * obliquity * launch / std::sqrt(in);
        const double rate = (t - sz) / in;
        const double theta0 = std::atan2(t, w);
        out[0].length = in + r + s3;
        out[0].rate = rate + t / r;
        out[0].amplitude = amplitude * hard_wedge_coefficient(sweep.at(theta0), k, r * s3 / (r + s3)) / std::sqrt(r * s3);
        if (count[1])
        {
            const double d = std::hypot(rx, rz - t);
            out[1] = {in + d, rate + (t - rz) / d, double(count[1]) * amplitude / std::sqrt(d)};
        }
        if (count[2])
        {
            const double d = std::hypot(rx, rz + t);
            out[2] = {in + d, rate + (t + rz) / d, double(count[2]) * amplitude / std::sqrt(d)};
        }
        if (count[3])
        {
            // The aperture point mirrored in the other face's plane.
            const double along = -w * fx + t * fz;
            const double ix = w + 2.0 * along * fx + w, iz = 2.0 * along * fz - t;
            const double d = std::hypot(rx - ix, rz - iz);
            const double rate_d = -((rx - ix) * 2.0 * fz * fx + (rz - iz) * (2.0 * fz * fz - 1.0)) / d;
            out[3] = {in + d, rate + rate_d, double(count[3]) * amplitude / std::sqrt(d)};
        }
    };

    const double smallest = std::max(0.01 * std::min(w, s1), 1e-6 * extent);
    complex sum(0.0, 0.0);
    for (std::size_t c = 0; c + 1 < cuts.size(); c++)
    {
        const double a = cuts[c], b = cuts[c + 1];
        if (b - a < 1e-9 * extent)
            continue;
        // Which waves reach the receiver here (hard_wedge_field's tests).
        const double middle = std::atan2(0.5 * (a + b), w);
        const int count[4] = {1, (std::abs(theta - middle) < pi ? 1 : 0) - (above ? 1 : 0),
                (theta < pi - middle ? 1 : 0) - (above ? 1 : 0),
                middle > (face.n2 - 1.0) * pi && theta > (2.0 * face.n2 - 1.0) * pi - middle ? 1 : 0};
        // Panels a sixth to a seventh of their distance from the nearest
        // feature, and at most the distance from the second edge.
        std::vector<double> ends {a};
        for (double t = a; t < b;)
        {
            double near = std::abs(t - features[0]);
            for (double f : features)
                near = std::min(near, std::abs(t - f));
            const double h = std::min(0.15 * (near + smallest), t + w);
            double next = t + h;
            if (next > b || b - next < 0.3 * h)
                next = b;
            ends.push_back(next);
            t = next;
        }
        // Each segment's ends are taken a hair inside it, on the side of a cut
        // the segment's waves belong to.
        const double inside = 1e-7 * (b - a);
        Wave first[4], middle_wave[4], last[4];
        waves(a + inside, count, first);
        for (std::size_t p = 0; p + 1 < ends.size(); p++)
        {
            waves(0.5 * (ends[p] + ends[p + 1]), count, middle_wave);
            waves(p + 2 == ends.size() ? b - inside : ends[p + 1], count, last);
            for (int j = 0; j < 4; j++)
                if (count[j])
                {
                    // Per unit path length the amplitude is A / (dl/dt).  The
                    // lengths grow along the aperture wherever a wave counts:
                    // its stationary point, the straight path through the
                    // aperture, lies where the second wedge shadows it, or
                    // below the face's plane.
                    sum += filon_panel(k, first[j].length, middle_wave[j].length, last[j].length,
                            first[j].amplitude / first[j].rate, middle_wave[j].amplitude / middle_wave[j].rate,
                            last[j].amplitude / last[j].rate);
                }
            std::copy(last, last + 4, first);
        }
    }
    // The end point at the edge gives the aperture's own coefficient,
    // (1 + cos g) / 2 exp(i pi/4) / (sqrt(2 pi k) sin g), g the source's angle
    // below the face; the wedge's is half its Keller coefficient on the face,
    // the face doubling it.  Both grow as 1/g, and their ratio tends to 1.
    const double sin_g = -sz / s1, cos_g = -sx / s1;
    if (sin_g > 1e-6)
    {
        const complex kirchhoff = 0.5 * (1.0 + cos_g) * std::polar(1.0, 0.25 * pi) / (std::sqrt(2.0 * pi * k) * sin_g);
        const complex keller = hard_wedge_coefficient(face.n1, 0.0, angle1(sx, sz), k, 1e9);
        sum *= 0.5 * keller / kirchhoff;
    }
    return sum;
}

} // namespace

bool shared_face_applies(const SharedFace& face, double sx, double sz, double rx, double rz)
{
    return (sz < 0.0 && sx < 0.0) || (rz < 0.0 && rx > face.width);
}

std::complex<double> shared_face_field(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    // Each end below the face's plane can be the aperture's source; each
    // evaluation is weighted by its source's angle below the plane, so that
    // the one whose source nears the plane gives way to the other smoothly,
    // and the weights swap with the ends.
    const double forward = sz < 0.0 && sx < 0.0 ? std::atan2(-sz, -sx) : 0.0;
    const double reverse = rz < 0.0 && rx > face.width ? std::atan2(-rz, rx - face.width) : 0.0;
    if (forward + reverse <= 0.0)
        return {std::numeric_limits<double>::quiet_NaN(), 0.0};
    complex sum(0.0, 0.0);
    if (forward > 0.0)
        sum += forward * aperture(face, sx, sz, rx, rz, k);
    if (reverse > 0.0)
        sum += reverse * aperture({face.width, face.n2, face.n1}, face.width - rx, rz, face.width - sx, sz, k);
    return sum / (forward + reverse);
}

} // namespace niceops
