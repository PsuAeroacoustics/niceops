import std.experimental.logger;

void main(string[] args)
{
    import sphere_database : Database;
    import trajectory : Trajectory;
    import std.range : iota;
    import std.array;
    import std.getopt;
    import core.stdc.stdlib : exit;

    string log_file = "plan-ops.log";
    string database_file = "databases/B407_M2.nc";
    string path_file = "test.csv";
    string vehicle_file = "vehicle.json";
    string out_file = "levels.csv";
    bool plotting = true;

    auto option_results = getopt(args, "d|database", "netCDF noise database file.",
            &database_file, "l|logfile", "Output log file.", &log_file,
            "t|trajectoryfile",
            "CSV trajectory input file.", &path_file,
            "v|vehicle", "JSON vehicle definition file.", &vehicle_file, "o|output",
            "Output file name.", &out_file, "p|plotting", "Enable plotting.", &plotting);

    if (option_results.helpWanted)
    {
        defaultGetoptPrinter("NICE-R-OPS low noise planning system.", option_results.options);
        exit(0);
    }

    //Set default logging to file
    sharedLog = new FileLogger(log_file);

    auto database = new Database(database_file);

    //Create trajectory from CSV data
    auto trajectory = new Trajectory(path_file, vehicle_file);

    //Trim trajectory (per Mike)
    double xstart_ft = -6000.;
    trajectory.watts_trim(xstart_ft);

    // Trajectory adaptation parameters
    double dt = 1.0;
    int degree = 5;
    double smooth = 1e6;
    trajectory.adapt(dt, degree, smooth);

    // Set up observer array
    auto x_obs = iota(-2500.0, 1100.0, 250.0).array;
    auto nx = x_obs.length;
    auto y_obs = iota(-1500.0, 1600.0, 250.0).array;
    auto ny = y_obs.length;

    double[] levels_linear;
    levels_linear.length = nx * ny;

    // And time series data
    auto nt = trajectory.states.length;
    double[] time_history;
    time_history.length = nt;

    foreach (iy, yo; y_obs)
    {
        foreach (ix, xo; x_obs)
        {
            foreach (it, const ref state; trajectory.states)
            {
                double[3] observer = [xo, yo, 0.0];
                double[3] source = [state.x, state.y, state.z];
                double[3] velocity = [state.vx, state.vy, state.vz];
                double[3] normal = [state.nx, state.ny, state.nz];
                time_history[it] = database.level_at_observer(observer, source,
                        velocity, normal, state.mu, state.fpa_e, state.CT);
            }
            auto index = ix + iy * nx;
            levels_linear[index] = sound_exposure_level(time_history, dt);
        }
    }

    write_levels(x_obs, y_obs, levels_linear, out_file);

    if (plotting)
    {
        import plt = matplotlibd.pyplot;

        auto x = trajectory.Get!("x").get();
        auto y = trajectory.Get!("y").get();
        auto z = trajectory.Get!("z").get();

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

        plt.figure();
        plt.contourf(x_obs, y_obs, levels_plot);
        plt.plot(x, y, "k-");
        plt.axis("equal");
        plt.colorbar();

        plt.figure();
        plt.plot(x, trajectory.Get!("CT").get());

        plt.figure();
        plt.plot(x, trajectory.Get!("mu").get());

        plt.figure();
        plt.plot(x, trajectory.Get!("fpa").get(), x, trajectory.Get!("fpa_e").get());

        plt.figure();
        plt.plot(x, z);

        plt.show();
        plt.clear();
    }
}

static double sound_exposure_level(double[] levels, double T)
{
    import std.math : log10;
    import std.algorithm.iteration : sum;

    auto index = peak_range(levels);
    double[] levels_slice = levels[index[0] .. (index[1] + 1)];
    double[] p;
    p.length = levels_slice.length;
    p[] = 10.0 ^^ (0.1 * levels_slice[]);
    return 10.0 * log10(T * sum(p));
}

static size_t[2] peak_range(double[] levels, double down = 10.0)
{
    import std.algorithm.searching : maxIndex;

    // Find the 10-dB down points
    size_t[2] index = [0, levels.length - 1];
    auto peak_position = levels.maxIndex;
    auto threshold = levels[peak_position] - down;

    for (auto i = peak_position; i < levels.length; i++)
    {
        if (levels[i] < threshold)
        {
            index[1] = i - 1;
            break;
        }
    }
    for (auto i = peak_position; i > 0; i--)
    {
        if (levels[i] <= threshold)
        {
            if (levels[i] == threshold)
            {
                index[0] = i;
            }
            else
            {
                index[0] = i + 1;
            }
            break;
        }
    }
    return index;
}

unittest
{
    //One peak
    double[] levels = [0.0, 10.0, 20.0, 25.0, 20.0, 10.0, 0.0];
    auto index = peak_range(levels);
    assert(index[0] == 2);
    assert(index[1] == 4);

    //One peak, exact threshold levels
    levels = [0.0, 15.0, 20.0, 25.0, 20.0, 15.0, 0.0];
    index = peak_range(levels);
    assert(index[0] == 1);
    assert(index[1] == levels.length - 2);

    //Soft peak, no 10 dB down points
    levels = [20.0, 20.0, 20.0, 25.0, 20.0, 20.0, 20.0];
    index = peak_range(levels);
    assert(index[0] == 0);
    assert(index[1] == levels.length - 1);

    //Two peaks
    levels = [20.0, 10.0, 20.0, 25.0, 20.0, 10.0, 20.0];
    index = peak_range(levels);
    assert(index[0] == 2);
    assert(index[1] == 4);
}

static void write_levels(double[] x_obs, double[] y_obs, double[] levels, string out_file)
{
    import std.file;
    import std.stdio;

    immutable ft2m = 0.3048;

    auto nx = x_obs.length;
    auto ny = y_obs.length;

    auto file = File(out_file, "w");
    //Write header
    file.writeln("Line 1: x(ft), Line 2: y(ft),  Lines 3+ levels (x cols y rows)");
    foreach (i, x; x_obs)
    {
        file.write(x / ft2m);
        if (i < (nx - 1))
        {
            file.write(",");
        }
        else
        {
            file.write("\n");
        }
    }
    foreach (i, y; y_obs)
    {
        file.write(y / ft2m);
        if (i < (ny - 1))
        {
            file.write(",");
        }
        else
        {
            file.write("\n");
        }
    }
    foreach (i, level; levels)
    {
        file.write(level);
        if ((i + 1) % nx)
        {
            file.write(",");
        }
        else
        {
            file.write("\n");
        }
    }
    file.write("\n");
    file.close();
}
