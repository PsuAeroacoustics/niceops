#include "geodetic.hpp"
#include <cstdio>
#include <exception>
using namespace niceops;
int main()
{
    LocalFrame f(Geodetic{10.0, 179.99, 0.0});
    for (double lon : {179.995, 180.0, -179.995, -179.99, 179.97})
    {
        try {
            Vec3 v = f.to_local(Geodetic{10.0, lon, 0.0});
            std::printf("lon %9.3f -> east %14.1f m  north %12.1f m\n", lon, v.x, v.y);
        } catch (const std::exception& e) { std::printf("lon %9.3f -> throws: %s\n", lon, e.what()); }
    }
    std::printf("(0.01 deg of longitude at 10 deg latitude is about 1095 m)\n");
}
