import std.experimental.logger;

import trajectory;
import worldtransforms;
import wgs84tocartesian;
import sphere_database : Database;

static immutable ft2m = 0.3048;

struct AcousticResult {
    double[] sel;
    double[] dnl;
    double[][] time_history;
    Trajectory trajectory;
    double[] power;

    this(size_t timesteps, size_t observers, Trajectory _trajectory) {
        trajectory = _trajectory;
        sel = new double[observers];
        dnl = new double[observers];
        time_history = new double[][](observers, timesteps);
        power = new double[timesteps];
    }
}

auto compute_trajectory(Trajectory trajectory, Database database, double[3][] observers, bool qsam, double dt, double ambient_level, double db_cutoff) {

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
    auto t1 = watch.peek.total!"msecs";

    //enforce(world_transform !is null, "Received invalid world transform");
    /+if (self_test)
    {
        // TODO: Check out this literal corner case
        //2020-05-25T16:02:22.223 [all] sphere_database.d:740:self_test Surface level disagreement of D:83.6494 versus N:85.3316 for φ = -90, θ = 180, μ  = 0.164848,γ  =  -4.25762, CT = 0
        //        database.eval(database.sphere_radius, -90., 180.,0.164848, -4.25762, 0.);
        //        exit(0);
        database.self_test(database_file);
    }+/

    info("Starting flight path");
    // And time series data
    immutable nt = trajectory.states.length;

    auto result = AcousticResult(nt, observers.length, trajectory);

    foreach(obs_idx; 0..observers.length) {
        result.time_history[obs_idx][] = -double.infinity;
    }

    import std.algorithm : sort, map;
    import std.range : enumerate, zip;
    import std.typecons : Tuple, tuple;
    alias IndexTuple = Tuple!(size_t, "index", double, "range", double[3], "observer");
    auto indexed_range = new IndexTuple[observers.length];

    foreach (it, const ref state; trajectory.states)
    {
        double[3] source = [state.x, state.y, state.z];
        double[3] velocity = [state.vx, state.vy, state.vz];
        double[3] normal = [state.nx, state.ny, state.nz];
        
        foreach(idx, ref obs; observers) {
            indexed_range[idx].observer[] = obs[];
            indexed_range[idx].index = idx;
            indexed_range[idx].range = sqrt((obs[0] - source[0])^^2.0 + (obs[1] - source[1])^^2.0 + (obs[2] - source[2])^^2.0);
        }

        auto sorted_observers = 
            indexed_range
                .sort!("a.range < b.range")
                .map!(a => tuple(a.index, a.observer)); //limit the search area of observers
        

        foreach (index, observer; sorted_observers)
        {
            auto level = database.level_at_observer(observer, source,
                        velocity, normal, state.mu, qsam ? state.fpa_e : state.fpa, state.CT, state.heading);

            result.time_history[index][it] = level;

            if(!ambient_level.isNaN) {
                result.time_history[index][it] = 20.0 * log10(
                    10.0 ^^ (result.time_history[index][it] / 20.) + 10.0 ^^ (ambient_level / 20.));
            }

            if(!db_cutoff.isNaN) {
                if(level < (db_cutoff - 10)) {
                    break;
                }
            }
        }
        // D. Zachos 6.12.2023
        // query power file
        //auto power = database.eval_power(state.mu, qsam ? state.fpa_e : state.fpa, state.CT); //do we need qsam?
        //result.power[it]=power;
    }

    foreach (index, observer; observers)
    {
        result.sel[index] = result.time_history[index].sound_exposure_level(dt);
        result.dnl[index] = result.time_history[index].day_night_level(dt, trajectory.states);
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
    import std.getopt : GetOptException, GetoptResult, getopt, config, defaultGetoptPrinter, arraySep;
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
    string[] state_files = [];//"";
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
    double db_cutoff; // meters

    bool smooth = false;
    double[2] smooth_params = [3, 0.5];
    double[] _smooth_params;
    double adv;
    double ct;
    double fpa;

    bool parallel_by_observer = false;

    arraySep = ",";
    GetoptResult option_results;
    
    try {
        option_results = getopt( args, config.caseSensitive,
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
            "state|S", "Filenames to save state data", &state_files,
            "self_test", "Run self tests.", &self_test,
            "stk", "Path, observer, and aircraft file are in STK format", &stk_format,
            "loglevel|L", "Log verbosisty level (options are: trace, info, warning, error, critical, fatal) default = info. These are case sensitive.", &log_level,
            "version", "Print niceops version and exit", &print_version,
            "timehistory|T", "Directory to save timehistory files to", &time_history_dir,
            "dnlfile|D", "Day night level result file", &dnl_file,
            "cutoff|c", "Cutoff range", &db_cutoff,
            "adv", "Plot sphere with advance ratio", &adv,
            "ct", "Plot sphere with C_T", &ct,
            "fpa", "Plot sphere with flight path angle", &fpa,
            "smooth", "Smooth flight path with degree and smoothness factor. Default = 3, 0.5", &_smooth_params,
            "pbo", "Parallellize by observer instead of flight path", &parallel_by_observer
        );

        if(_smooth_params.length != 0) {
            if(_smooth_params.length != 2) {
                error("Incorrect number of smoothing parameters");
                exit(-1);
            }
            smooth = true;
            smooth_params[0] = _smooth_params[0];
            smooth_params[1] = _smooth_params[1];
        }

        if(print_version) {
            import std.stdio : writeln;
            writeln("Niceops is currently not versioned");
        }

        /+if(time_history_dir == "") {
            time_history_dir = "./";
        }+/

        if (observer_file.length > 0)
        {
            scattered = true;
        }
        if (time_history_prefix.length > 0)
        {
            write_time_history = true;
        }
        if (state_files.length > 0)
        {
            write_state = true;
        }
        if (option_results.helpWanted)
        {
            defaultGetoptPrinter("Noise Informed Community Environment Operation Planning System (NICE-OPS).",
                    option_results.options);
            exit(0);
        }
    } catch(GetOptException ex) {
        //Set default logging to file
        sharedLog = new FileLogger(log_file, log_level, CreateFolder.yes);

        fatal("Failed to parse command line options: ", ex.msg);
        exit(-1);
    } catch(Exception ex) {
        //Set default logging to file
        sharedLog = new FileLogger(log_file, log_level, CreateFolder.yes);

        fatal("Caught unknown exception while parsing command line arguments:\n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
        exit(-1);
    }

    //Set default logging to file
    sharedLog = new FileLogger(log_file, log_level, CreateFolder.yes);

    info("database_files.length: ", database_files.length, ", database_files: ", database_files);

    Database[] databases;
    try {
        import std.algorithm : canFind;

        //string[] unique_database_files;
        Database[string] unique_database_files;
        foreach(idx, database_file;  database_files) {
            if(!unique_database_files.byKey.canFind(database_file)) {
                unique_database_files[database_file] = new Database(database_file);
            }
            databases ~= unique_database_files[database_file];
        }
        //databases = database_files.map!(a => new Database(a)).array;
    } catch(Exception ex) {
        //fatal("Caught exception while opening database files.\n", ex.msg);
        fatal("Caught exception while opening database files: \n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
        exit(-2);
    }

    if(!fpa.isNaN && !ct.isNaN && !adv.isNaN) {

        /+ct = (ct - databases[0].thrust_coefficient_min)/(databases[0].thrust_coefficient_max - databases[0].thrust_coefficient_min);
        adv = (ct - databases[0].advance_ratio_min)/(databases[0].advance_ratio_max - databases[0].advance_ratio_min);
        fpa = (ct - databases[0].flight_path_angle_min)/(databases[0].flight_path_angle_max - databases[0].flight_path_angle_min);

        import std.algorithm;
        import std.math : abs, approxEqual;
        import std.range;
        import std.stdio;

        auto adv_diff = databases[0].advance_ratios.map!(a => abs(adv - a)).array.sort;
        auto min_adv_diff = adv_diff[0];
        auto adv_selection = databases[0].advance_ratios.map!(a => abs(adv - a)).enumerate.filter!(a => approxEqual(min_adv_diff, a[1]));

        writeln(adv_selection);+/

        import std.algorithm;
        import std.math;
        import std.range;
        import std.stdio;

        auto thetas = iota(90, 180, 5);
        auto phis = iota(-180, 170, 5);

        auto num_elements = phis.length*thetas.length;
        double[] x = new double[num_elements];
        double[] y = new double[num_elements];
        double[] dba = new double[num_elements];

        // D. Zachos 6.13.2023
        power = databases[0].eval_power(adv, fpa, ct);

        foreach(i, theta; thetas.enumerate) {
            auto theta_r = theta*(PI/180.0);
            foreach(j, phi; phis.enumerate) {
                auto phi_r = phi*(PI/180.0);

                auto idx = j + i*phis.length;
                
                x[idx] = databases[0].sphere_radius*cos(theta_r)*cos(phi_r);
                y[idx] = databases[0].sphere_radius*cos(theta_r)*sin(phi_r);
                //double eval(double range, double phi, double theta, double advance_ratio, double flight_path_angle, double thrust_coefficient)
                dba[idx] = databases[0].eval(databases[0].sphere_radius, phi, theta, adv, fpa, ct);
            }
        }
        import plt = matplotlibd.pyplot;
        plt.tricontourf(x, y, dba, ["cmap": "YlOrRd"]);
        plt.colorbar;
        plt.axis("equal");
        plt.show;
    }


    if(!fpa.isNaN && !ct.isNaN && !adv.isNaN) {
        return;
    }

    double[] x_obs;
    double[] y_obs;
    double[3][] observers;
    size_t num_observers;

    import std.exception : enforce;

    WorldTransform world_transform;
    double[2][] lat_lon_obs;

    try{
        if (scattered)
        {
            // Read in scattered observer data
            auto parsed_obs = read_observer_file(observer_file, observers, stk_format);
            world_transform = parsed_obs.world_transform;
            lat_lon_obs = parsed_obs.lat_lon;

            //Set length for levels and observer x, y arrays
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
    } catch(Exception ex) {
        fatal("Caught exce\
        ption while setting up observer grid: \n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
        exit(-3);
    }

    AcousticResult[] acoustic_results;
    import std.range : zip;
    import std.parallelism : parallel, totalCPUs;

    if(!parallel_by_observer) {
        acoustic_results = new AcousticResult[path_files.length];

        foreach(idx, e; zip(path_files, databases, vehicle_files, out_files, state_files).parallel(1))
        {
            try{
                auto path_file = e[0];
                auto database = e[1];
                auto vehicle_file = e[2];
                auto out_file = e[3];
                auto state_file = e[4];

                //enforce(world_transform !is null, "Passing in invalid world transform");

                auto trajectory = new Trajectory( path_file, vehicle_file, world_transform, stk_format);

                if(smooth) {
                    // Trajectory adaptation parameters
                    immutable int degree = smooth_params[0].to!int;
                    immutable double smoothing = smooth_params[1];
                    trajectory.adapt(dt, degree, smoothing);
                }

                info("Computing trajectory for flight path ", idx);

                acoustic_results[idx] = trajectory.compute_trajectory(database, observers, qsam, dt, ambient_level, db_cutoff);
                write_levels(x_obs, y_obs, acoustic_results[idx].sel, out_file, scattered, stk_format);

                if(write_state) {
                    trajectory.write_state(state_file);
                }

                if(time_history_dir != "") {
                    foreach(obs_idx; 0..x_obs.length) {
                        auto output_file = time_history_dir ~ "/" ~ time_history_prefix ~ "_"~obs_idx.to!string~"_flightpath_" ~ idx.to!string ~ ".csv";
                        write_time_history_file(acoustic_results[idx].time_history[obs_idx], dt, x_obs[obs_idx], y_obs[obs_idx], output_file);
                    }
                }
            } catch(Exception ex) {
                error("Caught unknown exception while simulating trajectory of path file ", e[0], ": \n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
            }
        }
    } else {
        import std.algorithm : each, sum;
        import std.range : chunks, enumerate;

        Trajectory[] trajectories =
            zip(path_files, vehicle_files)
                .map!((z) {
                    auto t = new Trajectory(z[0], z[1], world_transform, stk_format);
                    if(smooth) t.adapt(dt, smooth_params[0].to!int, smooth_params[1]);
                    return t;
                })
                .array;

        acoustic_results = 
            trajectories.map!(t => AcousticResult(t.states.length, observers.length, t)).array;

        immutable observers_per_core = observers.length/totalCPUs;
		immutable left_over_observers = observers.length%totalCPUs;

        immutable observers_on_core =
            iota(0, totalCPUs)
            .map!(core => core < left_over_observers ? 1 : 0) // evenly distribute the remainder accross as many cores as possible
            .map!(a => a + observers_per_core)
            .array;

        double[3][][] split_observers;
        split_observers.length = totalCPUs;
        foreach(proc; 0..totalCPUs) {
            split_observers[proc] = new double[3][observers_on_core[proc]];
            
            immutable size_t slice_start = proc == 0 ? 0 : observers_on_core[0..proc].sum;
            immutable size_t slice_end = slice_start + observers_on_core[proc];

            split_observers[proc][] = observers[slice_start..slice_end];
        }
        import std.stdio : writeln;

        foreach(proc, observer_chunk; split_observers.parallel(1)) {

            foreach(idx, database; databases)
            {
                try{
                    //enforce(world_transform !is null, "Passing in invalid world transform");
                    info("Computing trajectory for flight path ", idx);

                    auto acoustic_result = trajectories[idx].compute_trajectory(database, observer_chunk, qsam, dt, ambient_level, db_cutoff);

                    immutable size_t slice_start = proc == 0 ? 0 : observers_on_core[0..proc].sum;
                    immutable size_t slice_end = slice_start + observers_on_core[proc];

                    acoustic_results[idx].time_history[slice_start..slice_end] = acoustic_result.time_history[];
                    acoustic_results[idx].sel[slice_start..slice_end] = acoustic_result.sel[];
                    acoustic_results[idx].dnl[slice_start..slice_end] = acoustic_result.dnl[];

                    if(time_history_dir != "") {
                        foreach(obs_idx; 0..observer_chunk.length) {
                            auto output_file = time_history_dir ~ "/" ~ time_history_prefix ~ "_"~(obs_idx + slice_start).to!string~"_flightpath_" ~ idx.to!string ~ ".csv";
                            write_time_history_file(acoustic_result.time_history[obs_idx], dt, observer_chunk[obs_idx][0], observer_chunk[obs_idx][1], output_file);
                        }
                    }
                } catch(Exception ex) {
                    error("Caught unknown exception while simulating trajectory of path file: \n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
                }
            }
        }

        foreach(idx, out_file; out_files) {
            write_levels(x_obs, y_obs, acoustic_results[idx].sel, out_file, scattered, stk_format);
        }
    }

    foreach(idx; 0..path_files.length) {
        if (plotting || figure_file.length > 0 || state_figure_file.length > 0)
        {
            import plt = matplotlibd.pyplot;
            import std.array : staticArray;

            auto t = acoustic_results[idx].trajectory.Get!("t").get();
            auto x = acoustic_results[idx].trajectory.Get!("x").get();
            auto y = acoustic_results[idx].trajectory.Get!("y").get();
            auto z = acoustic_results[idx].trajectory.Get!("z").get();
            auto vx = acoustic_results[idx].trajectory.Get!("vx").get();
            auto vy = acoustic_results[idx].trajectory.Get!("vy").get();
            auto vz = acoustic_results[idx].trajectory.Get!("vz").get();

            import std.stdio;

            auto obs_pos = zip(x_obs, y_obs).map!(a => staticArray!(double, 2)([a[0], a[1]]));

            double min_x = double.infinity;
            double max_x = -double.infinity;
            double min_y = double.infinity;
            double max_y = -double.infinity;

            //foreach(pos; wgs84_obs_pos) {
            foreach(pos; obs_pos) {
                import std.algorithm : min, max;
                min_x = min(pos[0], min_x);
                max_x = max(pos[0], max_x);
                min_y = min(pos[1], min_y);
                max_y = max(pos[1], max_y);
            }

            double x_span = max_x - min_x;
            double y_span = max_y - min_y;

            if (plotting || figure_file.length > 0)
            {
                import std.algorithm : min, max;
                if(db_cutoff.isNaN) {
                    contour_min = double.infinity;
                    contour_max = -double.infinity;
                    foreach(s; acoustic_results[idx].sel) {
                        if(!s.isNaN) {
                            contour_min = min(s, contour_min);
                            contour_max = max(s, contour_max);
                        }
                    }
                } else {
                    contour_min = db_cutoff;
                    contour_max = -double.infinity;
                    foreach(s; acoustic_results[idx].sel) {
                        if(!s.isNaN) {
                            contour_max = max(s, contour_max);
                        }
                    }
                }

                if(contour_max <= contour_min) {
                    contour_max = contour_min + contour_step;
                }

                //Set contour levels
                auto clevels = iota(contour_min - contour_step, contour_max + contour_step, contour_step).array();
                plt.figure();
                if (scattered)
                {
                    plt.tricontourf(x_obs, y_obs, acoustic_results[idx].sel, clevels, [
                            "cmap": "YlOrRd"
                            ]);
                    /+plt.tricontourf(wgs84_x_obs, wgs84_y_obs, acoustic_results[idx].sel/+, clevels+/, [
                            "cmap": "YlOrRd"
                            ]);+/
                    plt.plot(x_obs, y_obs, "ko", ["markersize": 0.1]);
                    plt.xlim(["left": min_x - 0.05*x_span], ["right": max_x + 0.05*x_span]);
                    plt.ylim(["bottom": min_y - 0.05*y_span], ["top": max_y + 0.05*y_span]);
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

                /+foreach(o_idx; 0..x_obs.length) {
                    plt.annotate(o_idx.to!string, [x_obs[o_idx], y_obs[o_idx]]);
                }+/

                plt.xlabel("X, m", ["fontsize": 14]);
                plt.ylabel("Y, m", ["fontsize": 14]);
                plt.xlim(["left": min_x - 0.05*x_span], ["right": max_x + 0.05*x_span]);
                plt.ylim(["bottom": min_y - 0.05*y_span], ["top": max_y + 0.05*y_span]);
                plt.colorbar();

                if (figure_file.length > 0)
                {
                    import std.path;

                    auto fig_file_no_ext = stripExtension(figure_file);
                    auto ext = extension(figure_file);
                    auto filename = fig_file_no_ext~"_"~idx.to!string~ext;
                    plt.savefig(filename);
                }
            }
            if (plotting || state_figure_file.length > 0)
            {
                import std.path;

                auto CT = acoustic_results[idx].trajectory.Get!("CT").get();
                auto mu = acoustic_results[idx].trajectory.Get!("mu").get();
                auto gamma = acoustic_results[idx].trajectory.Get!("fpa").get();
                auto gamma_e = acoustic_results[idx].trajectory.Get!("fpa_e").get();
                auto heading = acoustic_results[idx].trajectory.Get!("heading").get();
                auto v = acoustic_results[idx].trajectory.Get!("v").get();
                auto vn_hat = acoustic_results[idx].trajectory.Get!("vn_hat").get();

                plt.figure();
                plt.plot(x, CT);
                plt.xlabel("X, m", ["fontsize": 14]);
                plt.ylabel("$C_T$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_CT_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, CT);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("$C_T$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_CT_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(x, mu);
                plt.xlabel("X, m", ["fontsize": 14]);
                plt.ylabel("$\\mu$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_mu_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, mu);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("$\\mu$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_mu_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(x, gamma, x, gamma_e);
                plt.xlabel("X, m", ["fontsize": 14]);
                plt.ylabel("$\\gamma$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_gamma_" ~ extension(state_figure_file));
                }


                plt.figure();
                plt.plot(t, gamma_e);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("$\\gamma_e$", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_gamma_e_time_" ~ idx.to!string~ extension(state_figure_file));
                }


                plt.figure();
                plt.plot(t, heading);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("heading", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(
                            state_figure_file) ~ "_heading_time_" ~ idx.to!string~ extension(state_figure_file));
                }



                plt.figure();
                plt.plot(x, z);
                plt.xlabel("X, m", ["fontsize": 14]);
                plt.ylabel("Z, m", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_z_" ~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, v, t, vn_hat);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("v, m/s", ["fontsize": 14]);
                plt.legend(["v", "vn_hat"]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, z);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("Z, m", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_zvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, x);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("X, m", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_xvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, y);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("Y, m", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_yvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vx);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("V_x, m/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vxvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vy);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("V_y, m/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vyvt_" ~ idx.to!string~ extension(state_figure_file));
                }

                plt.figure();
                plt.plot(t, vz);
                plt.xlabel("t, s", ["fontsize": 14]);
                plt.ylabel("V_z, m/s", ["fontsize": 14]);
                if (state_figure_file.length > 0)
                {
                    plt.savefig(stripExtension(state_figure_file) ~ "_vzvt_" ~ idx.to!string~ extension(state_figure_file));
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
    if(ambient_level.isNaN) {
        final_dnl[] = 0;
    } else {
        final_dnl[] = ambient_dnl;
    }
    
    foreach(ref acoustic_result; acoustic_results) {
        foreach(idx; 0..final_dnl.length) {
            final_dnl[idx] += 10.0^^(acoustic_result.dnl[idx]/10.0);
        }
    }

    foreach(idx; 0..final_dnl.length) {
        final_dnl[idx] = 10.0*log10(final_dnl[idx]);
    }

    try {
        write_levels(x_obs, y_obs, final_dnl, dnl_file, scattered, stk_format);
    } catch(Exception ex) {
        fatal("Caught exception while writing dnl file: \n", ex.msg, "\n", ex.file, ":", ex.line, "\n", ex.info);
        exit(-4);
    }


    if(plotting) {
        import plt = matplotlibd.pyplot;
        import std.array : staticArray;
        import std.algorithm : min, max;
        if(db_cutoff.isNaN) {
            contour_min = double.infinity;
            contour_max = -double.infinity;
            foreach(s; final_dnl) {
                if(!s.isNaN) {
                    contour_min = min(s, contour_min);
                    contour_max = max(s, contour_max);
                }
            }
        } else {
            contour_min = db_cutoff;
            contour_max = -double.infinity;
            foreach(s; final_dnl) {
                if(!s.isNaN) {
                    contour_max = max(s, contour_max);
                }
            }
        }

        if(contour_max <= contour_min) {
            contour_max = contour_min + contour_step;
        }

        //Set contour levels
        auto clevels = iota(contour_min - contour_step, contour_max + contour_step, contour_step).array();
        plt.figure();

        auto obs_pos = zip(x_obs, y_obs).map!(a => staticArray!(double, 2)([a[0], a[1]]));

        double min_x = double.infinity;
        double max_x = -double.infinity;
        double min_y = double.infinity;
        double max_y = -double.infinity;

        foreach(pos; obs_pos) {
            import std.algorithm : min, max;
            min_x = min(pos[0], min_x);
            max_x = max(pos[0], max_x);
            min_y = min(pos[1], min_y);
            max_y = max(pos[1], max_y);
        }

        double x_span = max_x - min_x;
        double y_span = max_y - min_y;

        plt.tricontourf(x_obs, y_obs, final_dnl, clevels, [
                "cmap": "YlOrRd"
                ]);
        plt.plot(x_obs, y_obs, "ko", ["markersize": 0.1]);
        plt.xlim(["left": min_x - 0.05*x_span], ["right": max_x + 0.05*x_span]);
        plt.ylim(["bottom": min_y - 0.05*y_span], ["top": max_y + 0.05*y_span]);


        foreach(idx, ref acoustic_result; acoustic_results) {
            auto x = acoustic_result.trajectory.Get!("x").get();
            auto y = acoustic_result.trajectory.Get!("y").get();

            plt.plot(x, y, "k-");
        }

        plt.xlabel("X, m", ["fontsize": 14]);
        plt.ylabel("Y, m", ["fontsize": 14]);
        plt.xlim(["left": min_x - 0.05*x_span], ["right": max_x + 0.05*x_span]);
        plt.ylim(["bottom": min_y - 0.05*y_span], ["top": max_y + 0.05*y_span]);
        plt.colorbar();

        plt.show();
        plt.clear();
    }
}

@nogc static double sound_exposure_level(double[] levels, double T)
{
    import std.math : log10;
    import std.algorithm.iteration : sum, map;

    auto index = peak_range(levels);
    double[] levels_slice = levels[index[0] .. (index[1] + 1)];

    immutable p_sum = levels_slice.map!(l => 10.0^^(0.1*l)).sum;
    return 10.0 * log10(T * p_sum);
}

@nogc static double day_night_level(double[] levels, double T, Trajectory.State[] state)
{
    import std.math : log10;
    import std.algorithm : map, sum;
    import std.range : zip;

    auto p_levels = zip(levels, state).map!(a => a[1].day ? a[0] : a[0] + 10); // Add 10dB penalty for night ops.

    auto index = peak_range(p_levels);

    immutable p_sum = p_levels[index[0] .. (index[1] + 1)].map!(l => 10.0^^(0.1*l)).sum;
    return 10.0 * log10(T * p_sum);
}

@nogc static size_t[2] peak_range(R)(R levels, double down = 10.0)
{
    import std.algorithm.searching : maxIndex;

    // Find the 10-dB down points
    size_t[2] index = [0, levels.length - 1];
    auto peak_position = levels.maxIndex;
    immutable threshold = levels[peak_position] - down;

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
        const double[] levels, string out_file, bool scattered, bool stk_format = false) //duplicate something like this for power in it's own file D.Zachos
{
    import std.stdio : File;
    import std.file : exists, mkdirRecurse;
    import std.path : dirName;

    // Check if directory exists;  if not, create it.
    auto directory = dirName(out_file);
    if (!exists(directory))
    {
        mkdirRecurse(directory);
    }

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
        immutable nx = x_obs.length;
        immutable ny = y_obs.length;
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
    import std.stdio : File;
    import std.file : exists, mkdirRecurse;
    import std.path : dirName;

    // Check if directory exists;  if not, create it.
    auto directory = dirName(filename);
    if (!exists(directory))
    {
        mkdirRecurse(directory);
    }

    trace("Writing time history to ", filename);
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
    import std.stdio : File;
    import std.file : exists, mkdirRecurse;
    import std.path : dirName;

    // Check if directory exists;  if not, create it.
    auto directory = dirName(filename);
    if (!exists(directory))
    {
        mkdirRecurse(directory);
    }

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
    import std.typecons : tuple;

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

        double[2][] null_lat_lon;
        WorldTransform wt;
        return tuple!("world_transform", "lat_lon")(wt, null_lat_lon);
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
                    min_lat = min(min_lat, a.lat);
                    max_lat = max(max_lat, a.lat);
                    min_long = min(min_long, a.lon);
                    max_long = max(max_long, a.lon);
                    min_alt = min(min_alt, a.alt);
                    max_alt = max(max_alt, a.alt);
                    return Vector!(3, double)(a.lat, a.lon, a.alt);
                })
                .array;

        auto lat_lon_obs = lat_lon_alt_obs.map!(a => staticArray!(double, 2)([a[0], a[1]])).array;

        trace("min_lat: ", min_lat);
        trace("max_lat: ", max_lat);
        trace("min_long: ", min_long);
        trace("max_long: ", max_long);
        trace("min_alt: ", min_alt);
        trace("max_alt: ", max_alt);

        double lat_ref = (max_lat + min_lat)/2.;
        double lon_ref = (max_long + min_long)/2.;
        double alt_ref = (max_alt + min_alt)/2.;
        auto lat_long_alt_ref = Vector!(3, double)(lat_ref, lon_ref, alt_ref);

        trace("lat_long_alt_ref: ", lat_long_alt_ref);
        auto wt = new WorldTransform(lat_long_alt_ref);

        observers = 
            lat_lon_alt_obs
                .map!(a => wt.lat_long_alt_to_enu(Vector!(3, double)(a[0], a[1], a[2]/+*1000.0+/)).mData)
                .array;

        return tuple!("world_transform", "lat_lon")(wt, lat_lon_obs);
    }
}
