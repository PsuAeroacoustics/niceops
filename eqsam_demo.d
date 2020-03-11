import sphere_database : Database;
import splines : Spline;
import trajectory : Trajectory;
import std.experimental.logger;

import std.math;
import std.range;
import std.algorithm;
import plt = matplotlibd.pyplot;

void main()
{
    //Set default logging to file
    sharedLog = new FileLogger("plan-ops.log");
    auto database = new Database("databases/B407_M2.nc");

    auto range = 1000.;
    auto mu = 0.15;
    auto fpa = -6.0;
    auto CT = 0.004;

    database.eval(range, 45.2, 75.0, mu, fpa, CT);
    database.eval(range, -180.0, 180.0, mu, fpa, CT);

    double[3] observer = [0.0, 0.0, 0.0];
    double[3] source = [-500, -500, 100];
    double[3] velocity = [1.0, 0.0, 0.0];
    double[3] normal = [0.0, 0.0, 1.0];
    auto level = database.level_at_observer(observer, source, velocity, normal, mu, fpa, CT);

    source = [0, -500, 500];
    level = database.level_at_observer(observer, source, velocity, normal, mu, fpa, CT);
    source = [-500, 0, 500];
    level = database.level_at_observer(observer, source, velocity, normal, mu, fpa, CT);

    auto contour_level = 65.0;
    auto bearing = 0.0;
    auto altitude = 100.0;
    auto threshold = 1.0;
    auto zero_check = true;
    range = database.find_range(contour_level, bearing, altitude, velocity,
            normal, mu, fpa, CT, threshold, zero_check);

    //Create trajectory from CSV data
    auto trajectory = new Trajectory("test.csv");

    double[] t1 = trajectory.Get!("x").get();
    double[] x1 = trajectory.Get!("y").get();

    // Trajectory adaptation parameters
    double dt = 0.25;
    int degree = 5;
    double smooth = 1e6;
    trajectory.adapt(dt, degree, smooth);

    double[] t2 = trajectory.Get!("x").get();
    double[] x2 = trajectory.Get!("y").get();
    //Plot trajectory before and after adaptation
    plt.plot(t1, x1, t2, x2);
    plt.show();
    plt.clear();
}
