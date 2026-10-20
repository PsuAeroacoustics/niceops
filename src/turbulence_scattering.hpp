#pragma once

// Turbulence scattering into a building's shadow, for --buildings with
// --scattering (THEORY 6.5): the engineering model of Forssén, Hornikx,
// Botteldooren, Wei, Van Renterghem and Ögren, "A model of sound scattering
// by atmospheric turbulence for use in noise mapping calculations", Acta
// Acustica united with Acustica 100(5), 810-815 (2014), in its form without a
// canyon (its eq. 4):
//
//   L = b1 + b2 lg(h / 10 m) + b3 lg(h^2 / (d_S d_R) + 0.0012) + (10/3) lg(f / 1 kHz)
//
// dB re free field, h the screen's height above the source-receiver line and
// d_S, d_R the distances along that line to the foot of h.  The coefficients
// are velocity turbulence's (its Table IV): b1 = -52.8 + 10 lg C_v^2,
// b2 = 11.3, b3 = -17.1, C_v^2 in m^4/3 s^-2.  They hold the ground's
// doubling of energy at both ends; the air's attenuation is the caller's.
//
// NICEOPS's strength is HARMONOISE's gamma_T = C_T^2 / T0^2 + (22/3) C_v^2 / c0^2
// (THEORY 5.3).  Through the small angles a shadow scatters at, velocity and
// temperature turbulence scatter alike per unit gamma_T -- the cross section
// goes as C_v^2/c0^2 cos^2(theta/2) + 0.14 C_T^2/T0^2, 3/22 against 0.14 --
// so gamma_T is taken as velocity turbulence of C_v^2 = (3/22) c0^2 gamma_T,
// within 0.1 dB of either reading.  The velocity coefficients reproduce the
// paper's Figure 2 (its scattering integral at 1 kHz, C_v^2 = 1, d_S = 40 m) at
// its corners to 2 dB, less the 6 dB of the ground; its temperature
// coefficients are not used.

#include <cstddef>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

// The velocity turbulence of the same small-angle scattering as gamma_T.
double velocity_equivalent(double gamma_t, double sound_speed);

// The scattered level, dB re free field, for velocity turbulence C_v^2 at
// frequency f behind a thin screen h above the source-receiver line, d_S and
// d_R along it (m).  -infinity for no screen (h <= 0) or no turbulence.
double scattered_level(double h, double d_s, double d_r, double frequency, double cv2);

// The thin screen standing in for the buildings between a source and a
// receiver: of the roof tops above the line between them, where its plan
// crosses a footprint's side, the one farthest above it.  `shielded` where any
// stands above it; then h is its distance above the line, d_s and d_r the
// distances along the line from the source and the receiver to its foot.
struct ShadowScreen
{
    bool shielded = false;
    double h = 0.0, d_s = 0.0, d_r = 0.0;
};
ShadowScreen shadow_screen(const Buildings& buildings, const Vec3& source, const Vec3& receiver);

} // namespace niceops
