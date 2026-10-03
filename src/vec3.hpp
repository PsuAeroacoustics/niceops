#pragma once

// A three-component vector of doubles, replacing Eigen::Vector3d.
//
// Eigen was the last library the build required beyond netCDF, and NICE-OPS
// used about a hundredth of it: positions, velocities and TPP normals, with
// dot, cross, norm and normalize.  Nothing here templates over size or defers
// evaluation, because nothing needs to.
//
// The arithmetic is written to match what Eigen did operation for operation --
// dot as x*x + y*y + z*z in that order, norm as sqrt of it, normalized as a
// division by that norm rather than a multiplication by its reciprocal -- so
// swapping the type out moves results by round-off at most.

#include <cmath>
#include <ostream>

namespace niceops {

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;

    Vec3(double x_, double y_, double z_)
        : x(x_)
        , y(y_)
        , z(z_)
    {
    }

    // Indexed access, as the call sites that read a component by number use.
    double operator()(int index) const
    {
        return index == 0 ? x : (index == 1 ? y : z);
    }

    double& operator()(int index)
    {
        return index == 0 ? x : (index == 1 ? y : z);
    }

    double dot(const Vec3& other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 cross(const Vec3& other) const
    {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }

    double squaredNorm() const
    {
        return dot(*this);
    }

    double norm() const
    {
        return std::sqrt(squaredNorm());
    }

    // Eigen's normalized() leaves a zero vector alone rather than dividing by
    // zero; level_at_observer relies on that for a coincident source.
    Vec3 normalized() const
    {
        const double squared = squaredNorm();
        if (squared > 0.0)
        {
            const double length = std::sqrt(squared);
            return {x / length, y / length, z / length};
        }
        return *this;
    }

    void normalize()
    {
        *this = normalized();
    }

    // Length in the ground plane, for the geodetic tests.
    double horizontal_norm() const
    {
        return std::sqrt(x * x + y * y);
    }

    Vec3& operator+=(const Vec3& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vec3& operator/=(double scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }
};

inline Vec3 operator+(Vec3 a, const Vec3& b)
{
    return a += b;
}

inline Vec3 operator-(Vec3 a, const Vec3& b)
{
    return a -= b;
}

inline Vec3 operator*(double scalar, const Vec3& v)
{
    return {scalar * v.x, scalar * v.y, scalar * v.z};
}

inline Vec3 operator*(const Vec3& v, double scalar)
{
    return {v.x * scalar, v.y * scalar, v.z * scalar};
}

inline Vec3 operator/(const Vec3& v, double scalar)
{
    return {v.x / scalar, v.y / scalar, v.z / scalar};
}

// Matches the layout Eigen's transpose() printed into the trace log.
inline std::ostream& operator<<(std::ostream& os, const Vec3& v)
{
    return os << v.x << " " << v.y << " " << v.z;
}

} // namespace niceops
