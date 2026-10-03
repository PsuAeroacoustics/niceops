#include "sphere_database.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
using namespace niceops;
int main(int argc, char** argv)
{
    const char* path = argc > 1 ? argv[1] : "databases/S-76D_M3.nod";
    Database db(path);
    std::vector<Database::SourceState> states;
    for (double mu = 0.0; mu <= 0.30001; mu += 0.02)
        for (double gamma = -12.0; gamma <= 12.0001; gamma += 3.0)
            for (double ct = 0.003; ct <= 0.0101; ct += 0.0014)
                states.push_back(db.source_state(Vec3(0.0, 0.0, 150.0), Vec3(60.0, 10.0, 2.0),
                        Vec3(0.1, -0.05, 1.0), mu, gamma, ct));
    const auto bounds = db.direction_bounds(states);
    const double radius = db.sphere_radius();
    std::size_t checked = 0, violations = 0;
    double worst_db = 0.0, worst_theta = 0, worst_phi = 0;
    for (std::size_t s = 0; s < states.size(); s++)
    {
        const auto& F = states[s].frame;
        for (double theta : {0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0})
            for (double sign : {1.0, -1.0})   // nose axis (theta=0) and tail axis (theta=180)
                for (double phi = -180.0; phi < 180.0; phi += 2.5)
                    for (double range : {150.0, 400.0, 1200.0})
                    {
                        const double th = theta * M_PI / 180.0, ph = phi * M_PI / 180.0;
                        Vec3 d = sign * std::cos(th) * F.nose
                                 + std::sin(th) * std::cos(ph) * F.normal
                                 + std::sin(th) * std::sin(ph) * F.starboard;
                        d = d.normalized();
                        const Vec3 observer = states[s].position + range * d;
                        const Vec3 delta = observer - states[s].position;
                        const double r2 = delta.squaredNorm();
                        if (r2 <= radius * radius) continue;
                        const int bin = DirectionBounds::bin(F.nose.dot(delta), F.normal.dot(delta),
                                F.starboard.dot(delta), r2);
                        const double level = db.level_at_observer(observer, states[s]);
                        const double lhs = bounds.numerator(s, bin);
                        const double rhs = std::pow(10.0, 0.1 * level) * r2 * (1.0 - 1e-12);
                        checked++;
                        if (lhs < rhs)
                        {
                            violations++;
                            const double excess = 10.0 * std::log10(rhs / lhs);
                            if (excess > worst_db) { worst_db = excess; worst_theta = sign > 0 ? theta : 180 - theta; worst_phi = phi; }
                        }
                    }
    }
    std::printf("%s: checked %zu directions within 5 deg of the nose/tail axes; bound exceeded in %zu; worst %.3f dB at theta=%.2f phi=%.1f\n",
            path, checked, violations, worst_db, worst_theta, worst_phi);
}
