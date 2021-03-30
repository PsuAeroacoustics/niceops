import std.experimental.logger;

import trajectory;
import worldtransforms;
import sphere_database : Database;

static immutable ft2m = 0.3048;

struct AcousticResult {
    double[] sel;
    double[] dnl;
    double[][] time_history;
    bool[][] has_ambient;
    Trajectory trajectory;

    this(size_t timesteps, size_t observers, Trajectory _trajectory) {
        trajectory = _trajectory;
        sel = new double[observers];
        dnl = new double[observers];
        time_history = new double[][](observers, timesteps);
        has_ambient = new bool[][](observers, timesteps);
        foreach(ref am; has_ambient) {
            am[] = true;
        }
    }
}

auto compute_trajectory(Database database, string vehicle_file, string path_file, double[3][] observers, bool qsam, double dt, double ambient_level, double range_cutoff, WorldTransform world_transform, bool stk_format) {

    import sphere_database : Database;
    import trajectory : Trajectory;
    import std.range : iota;
    import std.array : array;
    import std.getopt : getopt, config, defaultGetoptPrinter;
    import core.stdc.stdlib : exit;
    import std.conv : to;
    import std.exception : enforce;
    import std.math : isNaN, log10, sqrt;
    import std.datetime.stopwatch : StopWatch, AutoStart;
    import std.typecons : tuple;
    // Start stopwatch
    immutable watch = StopWatch(AutoStart.yes);
    //auto database = new Database(database_file, range_cutoff);
    auto t1 = watch.peek.total!"msecs";
    //log("Database loaded in ", t1 / 1000.0, " seconds.");

    enforce(world_transform !is null, "Received invalid world transform");
    /+if (self_test)
    {
        // TODO: Check out this literal corner case
        //2020-05-25T16:02:22.223 [all] sphere_database.d:740:self_test Surface level disagreement of D:83.6494 versus N:85.3316 for φ = -90, θ = 180, μ  = 0.164848,γ  =  -4.25762, CT = 0
        //        database.eval(database.sphere_radius, -90., 180.,0.164848, -4.25762, 0.);
        //        exit(0);
        database.self_test(database_file);
    }+/

    //Create trajectory from CSV data
    auto trajectory = new Trajectory( path_file, vehicle_file, world_transform, stk_format);

    // Trajectory adaptation parameters
    immutable int degree = 5;
    immutable double smooth = 10.0;
    trajectory.adapt(dt, degree, smooth);

    // And time series data
    immutable nt = trajectory.states.length;

    bool[] day = new bool[nt];

    auto result = AcousticResult(nt, observers.length, trajectory);
    
    foreach(obs_idx; 0..observers.length) {
        if(!ambient_level.isNaN) {
            result.time_history[obs_idx][] = ambient_level;//-40000000000;//ambient_level;
        } else {
            result.time_history[obs_idx][] = 0;//-40000000000;
        }
    }
    import std.algorithm : sort, map;
    import std.range : enumerate, zip;
    
    foreach (it, const ref state; trajectory.states)
    {
        double[3] source = [state.x, state.y, state.z];
        double[3] velocity = [state.vx, state.vy, state.vz];
        double[3] normal = [state.nx, state.ny, state.nz];
        day[it] = state.day;
        
        auto sorted_observers =
            observers
                .map!(a => sqrt((a[0] - source[0])^^2.0 + (a[1] - source[1])^^2.0 + (a[2] - source[2])^^2.0 ))
                .zip(observers, iota(0, observers.length))
                .array
                .sort!"a[0] < b[0]"
                .map!(a => tuple(a[2], a[1]));

        foreach (index, observer; sorted_observers)
        {
            auto level = database.level_at_observer(observer, source,
                        velocity, normal, state.mu, qsam ? state.fpa_e : state.fpa, state.CT);

            result.time_history[index][it] = level;

            if (!isNaN(ambient_level))
            {
                //Add ambient noise
                /+result.time_history[index][it] = 20.0 * log10(
                        10.0 ^^ (result.time_history[index][it] / 20.) + 10.0 ^^ (ambient_level / 20.));+/

                if(level < ambient_level) {
                    //result.has_ambient[index][it] = true;
                    result.time_history[index][it] = ambient_level;//-40000000000;//ambient_level;
                } else {
                    result.has_ambient[index][it] = false;
                }
                
                if(level < (ambient_level - 10)) {
                    //result.has_ambient[index][it] = true;
                    trace("Killing at observer: ", index, ": ", observer, " level: ", level);
                    result.time_history[index][it] = ambient_level;//-40000000000;//ambient_level;
                    break;
                }
            } else if(level < 0) {
                //result.has_ambient[index][it] = true;
                result.time_history[index][it] = 0;//-40000000000;
            }
            
            if(level < -10 && isNaN(ambient_level)) {
                //result.has_ambient[index][it] = true;
                trace("Killing at observer: ", index, ": ", observer, " level: ", level);
                result.time_history[index][it] = 0;//-40000000000;
                break;
            }
        }
    }

    foreach (index, observer; observers)
    {
        result.sel[index] = sound_exposure_level(result.time_history[index], dt, ambient_level, result.has_ambient[index]);
        result.dnl[index] = day_night_level(result.time_history[index], dt, day);
    }

    auto t2 = watch.peek.total!"msecs";
    info("Acoustic calculations completed in ", (t2 - t1) / 1000.0, " seconds.");

    return result;
}

void main(string[] args)
{
    import trajectory : Trajectory;
    import std.range : iota;
    import std.array : array;
    import std.getopt : getopt, config, defaultGetoptPrinter, arraySep;
    import core.stdc.stdlib : exit;
    import std.conv : to;
    import std.math : isNaN, log10;
    import std.datetime.stopwatch : StopWatch, AutoStart;
    import std.algorithm : map;
    //Default to demo case
    string log_file = "log/niceops.log";
    string[] database_files = [/+"databases/S-76D_M3.nod"+/];
    string[] path_files = [/+"tracks/S-76D/197421AC.csv"+/];
    string[] vehicle_files = [/+"databases/S-76D.json"+/];
    string[] out_files = [/+"integrated/197421_LAE.csv"+/];
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

    bool stk_format = false;
    LogLevel log_level = LogLevel.info;
    bool print_version = false;
    string time_history_dir = "";
    string dnl_file = "integrated/dnl.csv";
    double range_cutoff = 5000; // meters

    arraySep = ",";
    auto option_results = getopt( args, config.caseSensitive,
        "database|d", "netCDF noise database file.", &database_files,
        "logfile|l", "Output log file.", &log_file,
        "trajectoryfile|t", "CSV trajectory input file.", &path_files,
        "vehicle|v", "JSON vehicle definition file.", &vehicle_files,
        "output|o", "Output file name.", &out_files,
        "plotting|p", "Enable plotting.", &plotting,
        "figfile|f", "Filename to save figure.", &figure_file,
        "statefigfile|F", "Filename to save state figures.", &state_figure_file,
        "cmin", "Minimum contour level, dB.", &contour_min,
        "cmax", "Maximum contour level, dB.", &contour_max,
        "cstep", "Contour level step, dB.", &contour_step,
        "step|s", "Time step in seconds.", &dt, 
        "ambient|a", "Ambient noise level, dB", &ambient_level,
        "history|i", "Time history file prefix.", &time_history_prefix,
        "metricfile|m", "File to append overall noise metrics.", &metric_file,
        "observer|b", "Observer file name.", &observer_file,
        "qsam|q", "Enable extended QSAM", &qsam,
        "state|S", "Filename to save state data", &state_file,
        "self_test", "Run self tests.", &self_test,
        "stk", "Path, observer, and aircraft file are in STK format", &stk_format,
        "loglevel|L", "Log verbosisty level (options are: trace, info, warning, error, critical, fatal) default = info", &log_level,
        "version", "Print niceops version and exit", &print_version,
        "timehistory|T", "Directory to save timehistory files to", &time_history_dir,
        "dnlfile|D", "Day night level result file", &dnl_file,
        "cutoff", "Cutoff range", &range_cutoff
    );

    if(print_version) {
        import std.stdio : writeln;
        writeln("Niceops is currently not versioned");
    }

    if(time_history_dir == "") {
        time_history_dir = "./";
    }

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
    sharedLog = new FileLogger(log_file, log_level, CreateFolder.yes);

    info("database_files.length: ", database_files.length, ", database_files: ", database_files);

    auto databases = database_files.map!(a => new Database(a, range_cutoff)).array;
    double[] x_obs;
    double[] y_obs;
    double[3][] observers;
    size_t num_observers;

    import std.exception : enforce;

    WorldTransform world_transform;
    if (scattered)
    {
        // Read in scattered observer data
        world_transform = read_observer_file(observer_file, observers, stk_format);
        
        enforce(world_transform !is null, "Failed to create world transform");
        //Set length for levels and observer x, y arrays
        //integrated_levels.length = observers.length;
        //day_night_levels.length = observers.length;
        num_observers = observers.length;
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
        //integrated_levels.length = nx * ny;
        //day_night_levels.length = nx * ny;
        observers.reserve(nx * ny);
        num_observers = nx*ny;
        foreach (iy, yo; y_obs)
        {
            foreach (ix, xo; x_obs)
            {
                observers ~= [xo, yo, 0.0];
            }
        }
    }

    enforce(path_files.length == vehicle_files.length, "Need the same number of vehicle files as flight path files");
    enforce(path_files.length == databases.length, "Need the same number of database files as flight path files");
    enforce(path_files.length == out_files.length, "Need the same number of SEL output files as flight path files");

    auto acoustic_results = new AcousticResult[path_files.length];
    import std.range : zip;
    import std.parallelism : parallel;

    foreach(idx, e; zip(path_files, databases, vehicle_files, out_files).parallel(1))
    {
        auto path_file = e[0];
        auto database = e[1];
        auto vehicle_file = e[2];
        auto out_file = e[3];

        enforce(world_transform !is null, "Passing in invalid world transform");
        acoustic_results[idx] = compute_trajectory(database, vehicle_file, path_file, observers, qsam, dt, ambient_level, range_cutoff, world_transform, stk_format);
        write_levels(x_obs, y_obs, acoustic_results[idx].sel, out_file, scattered, stk_format);
        //write_levels(x_obs, y_obs, acoustic_result.dnl, dnl_file, scattered, stk_format);

        foreach(obs_idx; 0..x_obs.length) {
            foreach(ref level; acoustic_results[idx].time_history[obs_idx]) {
                level = 20.0 * log10(10.0 ^^ (level / 20.) + 10.0 ^^ (ambient_level / 20.));
            }

            //static void write_time_history_file(double[] time_history, double dt, double x, double y, string filename)
            auto output_file = time_history_dir ~ "/" ~ time_history_prefix ~ "_"~obs_idx.to!string~"_flightpath_" ~ idx.to!string ~ ".csv";
            write_time_history_file(acoustic_results[idx].time_history[obs_idx], dt, x_obs[obs_idx], y_obs[obs_idx], output_file);
        }

        
        if (plotting || figure_file.length > 0 || state_figure_file.length > 0)
        {
            import plt = matplotlibd.pyplot;

            auto t = acoustic_results[idx].trajectory.Get!("t").get();
            auto x = acoustic_results[idx].trajectory.Get!("x").get();
            auto y = acoustic_results[idx].trajectory.Get!("y").get();
            auto z = acoustic_results[idx].trajectory.Get!("z").get();
            auto vx = acoustic_results[idx].trajectory.Get!("vx").get();
            auto vy = acoustic_results[idx].trajectory.Get!("vy").get();
            auto vz = acoustic_results[idx].trajectory.Get!("vz").get();

            import std.stdio;

            x_obs[] *= 1.0 / ft2m;
            y_obs[] *= 1.0 / ft2m;
            x[] *= 1.0 / ft2m;
            y[] *= 1.0 / ft2m;
            z[] *= 1.0 / ft2m;
            vx[] *= 1.0 / ft2m;
            vy[] *= 1.0 / ft2m;
            vz[] *= 1.0 / ft2m;

            if (plotting || figure_file.length > 0)
            {
                //Set contour levels
                auto clevels = iota(contour_min, contour_max + contour_step, contour_step).array();
                plt.figure();
                if (scattered)
                {
                    plt.tricontourf(x_obs, y_obs, acoustic_results[idx].sel/+, clevels+/, [
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
                            gridded_levels[iy][ix] = acoustic_results[idx].sel[ix + iy * nx];
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

                auto CT = acoustic_results[idx].trajectory.Get!("CT").get();
                auto mu = acoustic_results[idx].trajectory.Get!("mu").get();
                auto gamma = acoustic_results[idx].trajectory.Get!("fpa").get();
                auto gamma_e = acoustic_results[idx].trajectory.Get!("fpa_e").get();

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

                plt.figure();
                plt.plot(t, z);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("Z, ft", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_zvt_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, x);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("X, ft", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_xvt_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, y);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("Y, ft", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_yvt_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vx);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("V_x, ft/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vxvt_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vy);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("V_y, ft/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vyvt_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vz);
                plt.xlabel("t, ft", ["fontsize": 14]);
                plt.ylabel("V_z, ft/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vzvt_" ~ extension(state_figure_file));
                }

            }
            if (plotting)
            {
                plt.show();
            }
            plt.clear();
        }
    }

    double min_time = double.infinity;
    double max_time = -double.infinity;

    import std.algorithm : min, max;
    foreach(ref acoustic_result; acoustic_results) {
        min_time = min(min_time, acoustic_result.trajectory.utcstart);
        auto end_time = acoustic_result.trajectory.utcstart + acoustic_result.trajectory.states[$-1].t;
        max_time = max(max_time, end_time);
    }

    import std.math : round;
    double all_elapsed = max_time - min_time;
    long samples = round(all_elapsed/dt).to!long;
    ambient_level = 10.0^^(0.1*ambient_level);

    import std.math : log10;
    double ambient_dnl = 10.0*log10(dt*samples*ambient_level);

    double[] final_dnl = new double[observers.length];
    final_dnl[] = ambient_dnl;
    foreach(ref acoustic_result; acoustic_results) {
        final_dnl[] += acoustic_result.dnl[];
    }

    write_levels(x_obs, y_obs, final_dnl, dnl_file, scattered, stk_format);
}

static double sound_exposure_level(double[] levels, double T, double ambient_level, bool[] has_ambient)
{
    import std.math : log10;
    import std.algorithm.iteration : sum;

    auto index = peak_range(levels);
    double[] levels_slice = levels[index[0] .. (index[1] + 1)];

    foreach(idx, ref level; levels_slice) {
        if(!has_ambient[idx]) {
            level =
                20.0 * log10(
                    10.0 ^^ (level / 20.) + 10.0 ^^ (ambient_level / 20.)
                );
        }
    }

    double[] p;
    p.length = levels_slice.length;
    p[] = 10.0 ^^ (0.1 * levels_slice[]);
    return 10.0 * log10(T * sum(p));
}

static double day_night_level(double[] levels, double T, bool[] day)
{
    import std.math : log10;
    import std.algorithm;//.iteration : sum;
    import std.range;

    levels = levels.enumerate.map!(a => day[a.index] ? a.value : a.value + 10).array; // Add 10dB penalty for night ops.

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

static auto read_observer_file(string observer_file, out double[3][] observers, bool stk_format = false/+, WorldTransform world_transform = null+/)
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

        return null;
    }
    else
    {
        import std.exception : enforce;
        import std.typecons : tuple;
        import std.algorithm : min, max, sort;

        import trajectory : deg2rad;
        import numd.linearalgebra.matrix : Matrix, Vector;

        struct Layout {
            size_t point_num;
            double lat;
            double lon;
            double alt;
            double area;
            double lat_span;
            double lon_span;
        }

        import std.string : strip;
        double min_lat = double.infinity;
        double max_lat = -double.infinity;
        double min_long = double.infinity;
        double max_long = -double.infinity;
        double min_alt = double.infinity;
        double max_alt = -double.infinity;

        auto lat_lon_alt_obs = 
            observer_file
                .readText
                .strip
                .csvReader!Layout(null)
                .map!((a) {
                    min_lat = min(min_lat, a.lat*deg2rad);
                    max_lat = max(max_lat, a.lat*deg2rad);
                    min_long = min(min_long, a.lon*deg2rad);
                    max_long = max(max_long, a.lon*deg2rad);
                    min_alt = min(min_alt, a.alt);
                    max_alt = max(max_alt, a.alt);
                    return Vector!(3, double)(a.lat*deg2rad, a.lon*deg2rad, a.alt);
                })
                .array;

        trace("min_lat: ", min_lat);
        trace("max_lat: ", max_lat);
        trace("min_long: ", min_long);
        trace("max_long: ", max_long);
        trace("min_alt: ", min_alt);
        trace("max_alt: ", max_alt);

        double lat_ref = (max_lat + min_lat)/2.;
        double lon_ref = (max_long + min_long)/2.;
        double alt_ref = (max_alt + min_alt)/2.;
        auto lat_long_alt_ref = Vector!(3, double)(lat_ref, lon_ref, alt_ref*1000.0 /+ convert to meters +/);
        trace("lat_long_alt_ref: ", lat_long_alt_ref);
        auto wt = new WorldTransform(lat_long_alt_ref);

        observers = 
            lat_lon_alt_obs
                .map!(a => wt.lat_long_alt_to_enu(a).mData)
                .array;

        return wt;
    }
}
