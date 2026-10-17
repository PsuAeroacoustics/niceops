#include "double_wedge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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

// The source below the face's plane and before the first edge.  The aperture
// is the half-line x = 0, z > 0, sampled at a 24th of a wavelength out to 30
// wavelengths or three Fresnel zones, its far half tapered.  Each of its
// points sends a line source's wave to the second wedge, whose field -- the
// geometrical-optics waves past it and its diffraction -- the receiver sums.
complex aperture(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    const double w = face.width;
    const double lambda = 2.0 * pi / k;
    const double s1 = std::hypot(sx, sz);
    const double extent = std::max(3.0 * std::sqrt(lambda * std::max(s1, w)), 30.0 * lambda);
    const double step = lambda / 24.0;
    const double theta = angle2(rx - w, rz);
    const double s3 = std::hypot(rx - w, rz);
    // The second wedge's other face, for the reflection off it.
    const double fx = -std::cos(face.n2 * pi), fz = std::sin(face.n2 * pi);
    const complex launch = std::sqrt(k / (2.0 * pi)) * std::polar(1.0, -0.25 * pi);
    auto wave = [&](double x, double z) {
        const double d = std::hypot(rx - x, rz - z);
        return std::polar(1.0 / std::sqrt(d), k * d);
    };
    // Above the face's plane the receiver sees the first edge, and the
    // aperture's waves past the second edge hold that wedge's whole field --
    // the source's own wave where the wedge does not shadow it, and the edge's
    // diffracted wave with the face's reflection of it -- which the rays of
    // the first edge carry already.  The same aperture with the face going on
    // past the second edge is that field alone, in the aperture's own
    // approximation, and is taken back out: what remains is what the second
    // edge does to it.
    const bool above = rz > 0.0;
    complex sum(0.0, 0.0), unbounded(0.0, 0.0);
    for (double t = 0.5 * step; t < extent; t += step)
    {
        const double taper = t < 0.5 * extent ? 1.0 : 0.5 * (1.0 + std::cos(pi * (t - 0.5 * extent) / (0.5 * extent)));
        const double in = std::hypot(sx, t - sz);
        const double r = std::hypot(w, t);
        const double obliquity = 0.5 * (-sx / in + w / r);
        const complex amplitude = taper * obliquity * launch * std::polar(1.0 / std::sqrt(in), k * in);
        const double theta0 = angle2(-w, t);
        complex field = std::polar(1.0 / std::sqrt(r * s3), k * (r + s3))
                * hard_wedge_coefficient(face.n2, theta, theta0, k, r * s3 / (r + s3));
        // The incident wave where the second wedge does not shadow it, and its
        // reflections off the shared face and off the wedge's other face
        // (hard_wedge_field's tests).
        if (std::abs(theta - theta0) < pi)
            field += wave(0.0, t);
        if (theta < pi - theta0)
            field += wave(0.0, -t);
        if (theta0 > (face.n2 - 1.0) * pi && theta > (2.0 * face.n2 - 1.0) * pi - theta0)
        {
            const double px = -w, pz = t;
            const double along = px * fx + pz * fz;
            field += wave(w + 2.0 * along * fx - px, 2.0 * along * fz - pz);
        }
        sum += amplitude * field;
        if (above)
            unbounded += amplitude * (wave(0.0, t) + wave(0.0, -t));
    }
    sum = step * (sum - unbounded);
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
