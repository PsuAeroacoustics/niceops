import std.experimental.logger;

void main()
{
    import sphere_database : Database;
    import trajectory : Trajectory;
    import std.range : iota;
    import std.array;

    //Set default logging to file
    sharedLog = new FileLogger("plan-ops.log");

    auto database = new Database("databases/B407_M2.nc");

    //Create trajectory from CSV data
    auto trajectory = new Trajectory("test.csv", "vehicle.json");

    // Trajectory adaptation parameters
    double dt = 0.25;
    int degree = 5;
    double smooth = 1e6;
    trajectory.adapt(dt, degree, smooth);

    // Set up observer array
    auto nt = trajectory.states.length;
    auto x_obs = iota(-2500.0, 1100.0, 200.0).array;
    auto nx = x_obs.length;
    auto y_obs = iota(-1500.0, 1600.0, 200.0).array;
    auto ny = y_obs.length;

    double[] levels_linear;
    levels_linear.length = nt * nx * ny;

    auto it = 0;
    auto state = trajectory.states[nt / 2];

    foreach (iy, yo; y_obs)
    {
        foreach (ix, xo; x_obs)
        {
            // TODO: Calculate time series, then SEL for each observer
            //foreach (it, state; trajectory.states)
            //{
            double[3] observer = [xo, yo, 0.0];
            double[3] source = [state.x, state.y, state.z];
            double[3] velocity = [state.vx, state.vy, state.vz];
            double[3] normal = [state.nx, state.ny, state.nz];
            auto index = ix + iy * nx;
            levels_linear[index] = database.level_at_observer(observer, source,
                    velocity, normal, state.mu, state.fpa_e, state.CT);
            //}
        }
    }

    double[][] levels_plot;
    levels_plot.length = ny;
    for (auto iy = 0; iy < ny; iy++)
    {
        levels_plot[iy].length = nx;
        for (auto ix = 0; ix < nx; ix++)
        {
            levels_plot[iy][ix] = levels_linear[ix + iy * nx];
        }
    }

    immutable plotting = true;
    if (plotting)
    {
        import plt = matplotlibd.pyplot;

        double[] x = trajectory.Get!("x").get();
        double[] y = trajectory.Get!("y").get();
        double[] y1 = trajectory.Get!("mu").get();
        double[] y2 = trajectory.Get!("fpa_e").get();
        double[] y22 = trajectory.Get!("fpa").get();

        plt.figure();
        plt.plot(x, y1);
        plt.figure();
        plt.plot(x, y2, x, y22);

        plt.figure;
        plt.contourf(x_obs, y_obs, levels_plot);
        plt.plot(x, y, "k-");
        plt.colorbar();

        plt.show();
        plt.clear();
    }
}
