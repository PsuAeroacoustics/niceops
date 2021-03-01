import std.experimental.logger;

import worldtransforms;

static immutable ft2m = 0.3048;

void main(string[] args)
{
    import sphere_database : Database;
    import trajectory : Trajectory;
    import std.range : iota;
    import std.array : array;
    import std.getopt : getopt, config, defaultGetoptPrinter;
    import core.stdc.stdlib : exit;
    import std.conv : to;
    import std.math : isNaN, log10;
    import std.datetime.stopwatch : StopWatch, AutoStart;

    

    //Default to demo case
    string log_file = "log/niceops.log";
    string database_file = "databases/S-76D_M3.nod";
    string path_file = "tracks/S-76D/197421AC.csv";
    string vehicle_file = "databases/S-76D.json";
    string out_file = "integrated/197421_LAE.csv";
    string state_file = "";
    bool write_state = false;
    bool plotting = false;
    string figure_file = "";
    string state_figure_file = "";
    double contour_min = 65.0;
    double contour_max = 100.0;
    double contour_step = 5.0;
    bool scattered = false;
    string observer_file = "";
    double dt = 0.5;
    double ambient_level;
    bool write_time_history = false;
    string time_history_prefix = "";
    string metric_file = "";
    bool qsam = true;
    bool self_test = false;
    //string stk_vel_file = "";
    bool stk_format = false;

    auto option_results = getopt( args, config.caseSensitive,
    "database|d", "netCDF noise database file.", &database_file,
    "logfile|l", "Output log file.", &log_file, "trajectoryfile|t",
    "CSV trajectory input file.", &path_file, "vehicle|v",
    "JSON vehicle definition file.",
    &vehicle_file, "output|o", "Output file name.", &out_file,
    "plotting|p", "Enable plotting.", &plotting, "figfile|f",
    "Filename to save figure.", &figure_file,
    "statefigfile|F", "Filename to save state figures.",
    &state_figure_file, "cmin", "Minimum contour level, dB.", &contour_min,
    "cmax", "Maximum contour level, dB.", &contour_max, "cstep",
    "Contour level step, dB.", &contour_step, "step|s",
    "Time step in seconds.", &dt, "ambient|a", "Ambient noise level, dB",
    &ambient_level, "history|i", "Time history file prefix.",
    &time_history_prefix, "metricfile|m", "File to append overall noise metrics.", &metric_file,
    "observer|b", "Observer file name.",
    &observer_file, "qsam|q", "Enable extended QSAM", &qsam, "state|S",
    "Filename to save state data", &state_file, "self_test",
    "Run self tests.", &self_test,
    "stk", "Path, observer, and aircraft file are in STK format", &stk_format);

    if (observer_file.length > 0)
    {
        scattered = true;
    }
    if (time_history_prefix.length > 0)
    {
        write_time_history = true;
    }
    if (state_file.length > 0)
    {
        write_state = true;
    }
    if (option_results.helpWanted)
    {
        defaultGetoptPrinter("Noise Informed Community Environment Operation Planning System (NICE-OPS).",
                option_results.options);
        exit(0);
    }

    //Set default logging to file
    sharedLog = new FileLogger(log_file, LogLevel.trace, CreateFolder.yes);

    // Start stopwatch
    immutable watch = StopWatch(AutoStart.yes);
    auto database = new Database(database_file);
    auto t1 = watch.peek.total!"msecs";
    log("Database loaded in ", t1 / 1000.0, " seconds.");

    if (self_test)
    {
        // TODO: Check out this literal corner case
        //2020-05-25T16:02:22.223 [all] sphere_database.d:740:self_test Surface level disagreement of D:83.6494 versus N:85.3316 for φ = -90, θ = 180, μ  = 0.164848,γ  =  -4.25762, CT = 0
        //        database.eval(database.sphere_radius, -90., 180.,0.164848, -4.25762, 0.);
        //        exit(0);
        database.self_test(database_file);
    }

    //Create trajectory from CSV data
    auto trajectory = new Trajectory( path_file, vehicle_file, stk_format);

    //Trim trajectory (per Mike)
    immutable double xstart_ft = -6000.;
    trajectory.watts_trim(xstart_ft);

    // Trajectory adaptation parameters
    immutable int degree = 5;
    immutable double smooth = 0.5;
    trajectory.adapt(dt, degree, smooth);

    // And time series data
    immutable nt = trajectory.states.length;
    double[] time_history;
    time_history.length = nt;

    double[] integrated_levels;
    double[] x_obs;
    double[] y_obs;
    double[3][] observers;
    if (scattered)
    {
        // Read in scattered observer data
        read_observer_file(observer_file, observers, stk_format, trajectory.world_transform);
        //Set length for levels and observer x, y arrays
        integrated_levels.length = observers.length;
        x_obs.reserve(observers.length);
        y_obs.reserve(observers.length);
        //Deal observers out to x_obs and y_obs arrays
        foreach (observer; observers)
        {
            x_obs ~= observer[0];
            y_obs ~= observer[1];
        }
    }
    else
    {
        // Set up observer array
        x_obs = iota(-2000.0, 1100.0, 25.0).array;
        auto nx = x_obs.length;
        y_obs = iota(-1500.0, 1600.0, 25.0).array;
        auto ny = y_obs.length;
        //Set corresponding length for levels
        integrated_levels.length = nx * ny;
        observers.reserve(nx * ny);
        foreach (iy, yo; y_obs)
        {
            foreach (ix, xo; x_obs)
            {
                observers ~= [xo, yo, 0.0];
            }
        }
    }
    foreach (index, observer; observers)
    {
        foreach (it, const ref state; trajectory.states)
        {
            double[3] source = [state.x, state.y, state.z];
            double[3] velocity = [state.vx, state.vy, state.vz];
            double[3] normal = [state.nx, state.ny, state.nz];
            if (qsam)
            {
                time_history[it] = database.level_at_observer(observer, source,
                        velocity, normal, state.mu, state.fpa_e, state.CT);
            }
            else
            {
                time_history[it] = database.level_at_observer(observer, source,
                        velocity, normal, state.mu, state.fpa, state.CT);
            }
            if (!isNaN(ambient_level))
            {
                //Add ambient noise
                time_history[it] = 20.0 * log10(
                        10.0 ^^ (time_history[it] / 20.) + 10.0 ^^ (ambient_level / 20.));
            }
        }
        integrated_levels[index] = sound_exposure_level(time_history, dt);
        if (write_time_history)
        {
            //Write out time history files
            auto time_history_filename = time_history_prefix ~ "_i" ~ to!string(index) ~ ".csv";
            write_time_history_file(time_history, dt, observer[0], observer[1],
                    time_history_filename);
        }
    }

    auto t2 = watch.peek.total!"msecs";
    log("Acoustic calculations completed in ", (t2 - t1) / 1000.0, " seconds.");
    write_levels(x_obs, y_obs, integrated_levels, out_file, scattered, stk_format);

    if (metric_file.length > 0)
    {
        import std.path : baseName, stripExtension;

        auto casename = baseName(stripExtension(path_file));
        append_metrics_file(metric_file, casename, integrated_levels);
    }

    if (write_state)
    {
        trajectory.write_state(state_file);
    }

    if (plotting || figure_file.length > 0 || state_figure_file.length > 0)
    {
        import plt = matplotlibd.pyplot;

        auto x = trajectory.Get!("x").get();
        auto y = trajectory.Get!("y").get();
        auto z = trajectory.Get!("z").get();
        import std.stdio;

        //Convert to ft for plotting
        // TODO add units throughout

        /+if(stk_format) {
            import numd.linearalgebra.matrix : Vector;
            import trajectory : rad2deg;
            foreach(idx, ref obs; observers) {
                auto lat_lon_alt = trajectory.world_transform.enu_to_lat_lon(Vector!(3, double)(obs));
                x_obs[idx] = lat_lon_alt[1]*rad2deg;
                y_obs[idx] = /+90.0 - +/(lat_lon_alt[0]*rad2deg);
            }

            foreach(idx, ref state; trajectory.states) {
                auto xyz = Vector!(3, double)(state.x, state.y, state.z);
                auto lat_lon_alt = trajectory.world_transform.enu_to_lat_lon(xyz);
                x[idx] = lat_lon_alt[1]*rad2deg;
                y[idx] = /+90.0 - +/(lat_lon_alt[0]*rad2deg);
                z[idx] = lat_lon_alt[2];
            }

            writeln(y);
        } else {+/
            x_obs[] *= 1.0 / ft2m;
            y_obs[] *= 1.0 / ft2m;
            x[] *= 1.0 / ft2m;
            y[] *= 1.0 / ft2m;
            z[] *= 1.0 / ft2m;
        //}
        

        if (plotting || figure_file.length > 0)
        {
            //Set contour levels
            auto clevels = iota(contour_min, contour_max + contour_step, contour_step).array();
            plt.figure();
            if (scattered)
            {
                plt.tricontourf(x_obs, y_obs, integrated_levels/+, clevels+/, [
                        "cmap": "YlOrRd"
                        ]);
                plt.plot(x_obs, y_obs, "ko", ["markersize": 0.1]);
            }
            else
            {
                double[][] gridded_levels;
                auto nx = x_obs.length;
                auto ny = y_obs.length;
                gridded_levels.length = ny;
                for (auto iy = 0; iy < ny; iy++)
                {
                    gridded_levels[iy].length = nx;
                    for (auto ix = 0; ix < nx; ix++)
                    {
                        gridded_levels[iy][ix] = integrated_levels[ix + iy * nx];
                    }
                }
                plt.contourf(x_obs, y_obs, gridded_levels, clevels, [
                        "cmap": "YlOrRd"
                        ]);
            }
            plt.plot(x, y, "k-");
            plt.xlabel("X, ft", ["fontsize": 14]);
            plt.ylabel("Y, ft", ["fontsize": 14]);
            plt.axis("equal");
            plt.colorbar();
            plt.tight_layout();

            if (figure_file.length > 0)
            {
                plt.savefig(figure_file);
            }
        }
        if (plotting || state_figure_file.length > 0)
        {
            import std.path;

            auto CT = trajectory.Get!("CT").get();
            auto mu = trajectory.Get!("mu").get();
            auto gamma = trajectory.Get!("fpa").get();
            auto gamma_e = trajectory.Get!("fpa_e").get();

            plt.figure();
            plt.plot(x, CT);
            plt.xlabel("X, ft", ["fontsize": 14]);
            plt.ylabel("$C_T$", ["fontsize": 14]);
            if (state_figure_file.length > 0)
            {
                plt.savefig(stripExtension(
                        state_figure_file) ~ "_CT_" ~ extension(state_figure_file));
            }

            plt.figure();
            plt.plot(x, mu);
            plt.xlabel("X, ft", ["fontsize": 14]);
            plt.ylabel("$\\mu$", ["fontsize": 14]);
            if (state_figure_file.length > 0)
            {
                plt.savefig(stripExtension(
                        state_figure_file) ~ "_mu_" ~ extension(state_figure_file));
            }

            plt.figure();
            plt.plot(x, gamma, x, gamma_e);
            plt.xlabel("X, ft", ["fontsize": 14]);
            plt.ylabel("$\\gamma$", ["fontsize": 14]);
            if (state_figure_file.length > 0)
            {
                plt.savefig(stripExtension(
                        state_figure_file) ~ "_gamma_" ~ extension(state_figure_file));
            }

            plt.figure();
            plt.plot(x, z);
            plt.xlabel("X, ft", ["fontsize": 14]);
            plt.ylabel("Z, ft", ["fontsize": 14]);
            if (state_figure_file.length > 0)
            {
                plt.savefig(stripExtension(state_figure_file) ~ "_z_" ~ extension(state_figure_file));
            }
        }
        if (plotting)
        {
            plt.show();
        }
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

static void write_levels(const double[] x_obs, const double[] y_obs,
        const double[] levels, string out_file, bool scattered, bool stk_format = false)
{
    import std.file;
    import std.stdio;

    info("Writing output to ", out_file);
    auto file = File(out_file, "w");

    if (scattered)
    {
        info("Scattered output.");
        if(!stk_format) {
            assert(levels.length == x_obs.length && levels.length == y_obs.length);
            file.writeln("x(ft), y(ft), Level(dB)");
            foreach (i, level; levels)
            {
                file.writeln(x_obs[i] / ft2m, ",", y_obs[i] / ft2m, ",", level);
            }
        } else {
            file.writeln("observer, Level(dB)");
            foreach (i, level; levels)
            {
                file.writeln(i, ",", level);
            }
        }
    }
    else
    {
        info("Gridded output.");
        auto nx = x_obs.length;
        auto ny = y_obs.length;
        assert(nx * ny == levels.length);

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
    }
    file.close();
}

static void write_time_history_file(double[] time_history, double dt, double x,
        double y, string filename)
{
    import std.stdio;

    info("Writing time history to ", filename);
    auto file = File(filename, "w");
    file.writeln("Level time history. Next row [dt (s), x (ft), y(ft)], remaining rows, dB.");
    file.writeln(dt, ",", x / ft2m, ",", y / ft2m);
    foreach (level; time_history)
    {
        file.writeln(level);
    }
}

static void append_metrics_file(string filename, string casename,
        const double[] levels, double threshold = 65.0)
{
    import std.stdio;

    auto file = File(filename, "a");
    double average = 0;
    double maximum = -double.infinity;
    size_t count = 0;
    // Iterate over range and compute average and maximum value
    foreach (level; levels)
    {
        average += level;
        maximum = level > maximum ? level : maximum;
        if (level >= threshold)
        {
            count++;
        }
    }
    average /= levels.length;
    auto area = double(count) / levels.length;
    file.writeln(casename, ",", average, ",", maximum, ",", area);
}

static void read_observer_file(string observer_file, out double[3][] observers, bool stk_format = false, WorldTransform world_transform = null)
{
    import std.conv : to;
    import std.csv : csvReader;
    import std.file : readText;
    import std.array : array, staticArray;
    import std.string : indexOf, split, isNumeric;
    import std.range : dropOne;
    import std.algorithm : filter, map, strip;
    import std.stdio : writeln;

    immutable ext = observer_file.split(".")[$-1];

    if(!stk_format)
    {
        trace("Reading observers from " ~ observer_file ~ ".");
        auto contents = readText(observer_file);
        auto records = csvReader!(double)(contents, null);
        foreach (record; records)
        {
            auto data = record.array;
            if (data.length == 3)
            {
                observers ~= [data[0] * ft2m, data[1] * ft2m, data[2] * ft2m];
            }
            else if (data.length == 4)
            {
                //Skip channel ID in first column
                observers ~= [data[1] * ft2m, data[2] * ft2m, data[3] * ft2m];
            }
            else
            {
                fatal("Observer file format unreadable.");
            }
        }
    }
    else
    {
        import std.exception : enforce;
        import std.typecons : tuple;
        import std.algorithm : sort;

        import trajectory : deg2rad;
        import numd.linearalgebra.matrix : Matrix, Vector;

        enforce(world_transform !is null, "Need a valid world transform to parse STK observer files");

        struct Layout {
            /++
                "Point Number",
                "Latitude (deg)",
                "Longitude (deg)",
                "Altitude (km)",
                "Area (km^2)",
                "Latitude Span (deg)",
                "Longitude Span (deg)"
            +/
            size_t point_num;
            double lat;
            double lon;
            double alt;
            double area;
            double lat_span;
            double lon_span;
        }

        //auto records = file_contents.csvReader!Layout(null);
        observers = 
            observer_file
                .readText
                .csvReader!Layout(null)
                .map!((a) {
                    auto lat_lon_alt = Vector!(3, double)(a.lat*deg2rad, a.lon*deg2rad, a.alt);
                    auto obs_pos_vec = world_transform.lat_long_alt_to_enu(lat_lon_alt);
                    return obs_pos_vec.mData;
                })
                .array;
        /+observers =
            observer_file
                .readText
                .split("\n")
                .map!(a => a.strip(' '))
                .filter!(a => a.length > 0 && a[0].to!string.isNumeric)
                .dropOne
                .map!((line) {
                    auto split_line = line.split(" ").filter!(b => b != "");
                    auto index = split_line.front;
                    auto lat_lon_alt_arr = line.split(" ").filter!(b => b != "").dropOne.map!(b => b.to!double).staticArray!3;

                    auto lat_lon_alt = Vector!(3, double)(lat_lon_alt_arr[0]*deg2rad, lat_lon_alt_arr[1]*deg2rad, lat_lon_alt_arr[2]*1000);

                    auto obs_pos_vec = world_transform.lat_long_alt_to_enu(lat_lon_alt);

                    return tuple(index, obs_pos_vec.mData);
                })
                .array
                .sort!((a, b) => a[0] < b[0])
                .map!(a => a[1])
                .array;
        +/

        //writeln(observers);
    }
}
