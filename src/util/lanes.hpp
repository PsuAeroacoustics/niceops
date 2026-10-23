#pragma once

// Two doubles worked on together: the GCC/Clang vector type where there is
// one, and a plain struct that says the same thing where there is not.  For
// loops the compiler will not vectorize on its own -- the PE's march across
// Pade terms (pe.cpp), and the band-averaged pair sum (urban.cpp), whose
// double divides it judges not worth it.  Either spelling indexes with [], so
// code using it is written once.

#include <cstddef>
#include <cstring>

namespace niceops {

#if defined(__GNUC__)

// Clang defines __GNUC__ too, so the one test covers both.
typedef double Lanes __attribute__((vector_size(2 * sizeof(double))));

#else

struct Lanes
{
    double v[2];

    double operator[](std::size_t lane) const
    {
        return v[lane];
    }
    double& operator[](std::size_t lane)
    {
        return v[lane];
    }
};

inline Lanes operator+(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] + b.v[0], a.v[1] + b.v[1]}};
}

inline Lanes operator-(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] - b.v[0], a.v[1] - b.v[1]}};
}

inline Lanes operator*(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] * b.v[0], a.v[1] * b.v[1]}};
}

inline Lanes operator/(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] / b.v[0], a.v[1] / b.v[1]}};
}

inline Lanes& operator+=(Lanes& a, Lanes b)
{
    return a = a + b;
}

#endif

inline Lanes both(double value)
{
    return Lanes {value, value};
}

// Two consecutive doubles, wherever they lie.
inline Lanes load_lanes(const double* p)
{
    Lanes v;
    std::memcpy(&v, p, sizeof v);
    return v;
}

} // namespace niceops
