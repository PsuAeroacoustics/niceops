module trajectory;

import std.csv;
import std.file : readText;
import std.algorithm.iteration : splitter;
import std.datetime;
import std.format;
import std.stdio;
import std.range : retro, take;
import std.experimental.logger;
import std.array;
import std.uni : toLower;
import std.math : isNaN, sqrt, PI, asin;
import std.string;

import splines : Spline;

import wgs84tocartesian : toCartesian;

import worldtransforms;

import numd.linearalgebra.matrix;

immutable g = 9.80665; // m/s^^2
immutable mps2g = 1.0 / g;
immutable ft2m = 0.3048;
immutable fps2g = mps2g * ft2m;
immutable deg2rad = PI/180.0;
immutable rad2deg = 180.0/PI;

auto parse_greenwood(string file_contents, Trajectory.dimensional_data vehicle) {
    import std.typecons : tuple;

    double utcstart;
    Trajectory.State[] states;

    auto records = csvReader!(double)(file_contents, null);

    foreach (record; records)
    {
        auto data = record.array;
        if (data.length == 27)
        {
            //Assume NASA/FAA AC state data format
            Trajectory.State now;
            if (isNaN(utcstart))
            {
                utcstart = data[0];
            }
            // Convert to metric
            now.t = data[1];
            now.x = data[5] * ft2m;
            now.y = data[6] * ft2m;
            now.z = data[7] * ft2m;
            now.vx = data[11] * ft2m;
            now.vy = data[13] * ft2m;
            now.vz = data[15] * ft2m;
            // Convert to "g"
            now.nx = data[21] * fps2g;
            now.ny = data[22] * fps2g;
            now.nz = (data[23] * fps2g + 1.0);
            //Compute load factor
            now.LF = sqrt(now.nx ^^ 2 + now.ny ^^ 2 + now.nz ^^ 2);
            now.nx /= now.LF;
            now.ny /= now.LF;
            now.nz /= now.LF;
            //Non-dimensionalize
            now.CT = now.LF * vehicle.weight / (
                    vehicle.rho * vehicle.tip_speed ^^ 2 * PI * vehicle.radius ^^ 2);
            auto v = sqrt(now.vx ^^ 2 + now.vy ^^ 2 + now.vz ^^ 2);
            now.mu = v / vehicle.tip_speed;
            auto vn_hat = (now.vx * now.nx + now.vy * now.ny + now.vz * now.nz) / v;
            now.vn_hat = vn_hat;
            now.fpa_e = 180.0 / PI * asin(vn_hat);
            now.fpa = 180.0 / PI * asin(now.vz / v);
            states ~= now;
            trace("Parsed state: ", states[$ - 1]);
        }
        else
        {
            fatal("Wrong number of columns.");
        }
    }

    return tuple!("states", "utcstart"/+, "world_transform"+/)(states, utcstart/+, new WorldTransform+/);
}

auto month_str_to_int(string month) {
    final switch(month) {
        case "Jan":
            return 1;
            break;
        case "Feb":
            return 2;
            break;
        case "Mar":
            return 3;
            break;
        case "Apr":
            return 4;
            break;
        case "May":
            return 5;
            break;
        case "Jun":
            return 6;
            break;
        case "Jul":
            return 7;
            break;
        case "Aug":
            return 8;
            break;
        case "Sep":
            return 9;
            break;
        case "Oct":
            return 10;
            break;
        case "Nov":
            return 11;
            break;
        case "Dec":
            return 12;
            break;
    }
}

double time_str_to_double(string timestamp) {
    import std.conv : to;
    import std.algorithm : filter;
    import std.stdio : writeln;

    auto split_datetime = timestamp.split(" ").filter!(a => a != "").array;
    auto split_time = split_datetime[3].split(":").filter!(a => a != "").array;
    auto split_seconds = split_time[2].split(".").filter!(a => a != "").array;

    auto std_format_str =
        format!("%s%02d%02dT%02d%02d%02d.%sZ")(
            split_datetime[2],
            split_datetime[1].month_str_to_int,
            split_datetime[0].to!int,
            split_time[0].to!int,
            split_time[1].to!int,
            split_seconds[0].to!int,
            split_seconds[1]
        );

    auto sys_time = SysTime.fromISOString(std_format_str);
    auto timeval = sys_time.toTimeVal;
    double fraction = timeval.tv_usec.to!double/(1000.*1000.);
    return timeval.tv_sec.to!double + fraction;
}

auto parse_stk(string file_contents, Trajectory.dimensional_data vehicle, WorldTransform world_transform) {
    import std.algorithm : filter, map;
    import std.conv : to;
    import std.exception : enforce;
    import std.range : enumerate, slide;
    import std.stdio : writeln;
    import std.typecons : tuple;

    enforce(world_transform !is null, "Invalid world transform");
    double utcstart;
    Trajectory.State[] states;

    struct Layout {
        string utc_timestamp;
        double lat;
        double lon;
        double alt;
        double lat_rate;
        double lon_rate;
        double alt_rate;
        double yaw;
        double pitch;
        double roll;
        double bank;
        double heading;
        double acc_mag;
        int day;
    }

    trace("Parsing trajectory file as csv");
    auto records = file_contents.strip.csvReader!Layout(null);
    trace("CSV parsed");
    auto start_record = records.front;

    utcstart = start_record.utc_timestamp.time_str_to_double;

    import std.array : staticArray;

    double min_lat = double.infinity;
    double max_lat = -double.infinity;
    double min_long = double.infinity;
    double max_long = -double.infinity;
    double min_alt = double.infinity;
    double max_alt = -double.infinity;

    import std.algorithm : max, min;
    states = records
        .enumerate
        .map!((record) {
            Trajectory.State now;
            now.t = record.value.utc_timestamp.time_str_to_double - utcstart;

            auto lat_long_alt_point = Vector!(3, double)(
                record.value.lat,
                record.value.lon,
                record.value.alt
            );

            min_lat = min(min_lat, lat_long_alt_point[0]);
            max_lat = max(max_lat, lat_long_alt_point[0]);
            min_long = min(min_long, lat_long_alt_point[1]);
            max_long = max(max_long, lat_long_alt_point[1]);
            min_alt = min(min_alt, lat_long_alt_point[2]);
            max_alt = max(max_alt, lat_long_alt_point[2]);

            now.lat = record.value.lat;
            now.lon = record.value.lon;

            auto point = world_transform.lat_long_alt_to_enu(lat_long_alt_point);

            now.x = point[0];
            now.y = point[1];
            now.z = point[2];

            now.day = record.value.day == 1;

            now.heading = record.value.heading*deg2rad;
            return now;
        })
        .array;

    trace("min_lat: ", min_lat);
    trace("max_lat: ", max_lat);
    trace("min_long: ", min_long);
    trace("max_long: ", max_long);
    trace("min_alt: ", min_alt);
    trace("max_alt: ", max_alt);
    // Compute accelerations and other derived quantities.
    foreach(idx, ref state; states) {
        auto gravity_vec = Vector!(3, double)(0.0, 0.0, 1.0);

        if(idx == 0) {
            immutable dt = states[idx + 1].t - state.t;

            state.vx = (states[idx + 1].x - state.x)/dt;
            state.vy = (states[idx + 1].y - state.y)/dt;
            state.vz = (states[idx + 1].z - state.z)/dt;

            state.nx = gravity_vec[0];
            state.ny = gravity_vec[1];
            state.nz = gravity_vec[2];
        } else {
            immutable dt = state.t - states[idx - 1].t;

            state.vx = (state.x - states[idx - 1].x)/dt;
            state.vy = (state.y - states[idx - 1].y)/dt;
            state.vz = (state.z - states[idx - 1].z)/dt;

            immutable x_acc = (state.vx - states[idx - 1].vx)/dt;
            immutable y_acc = (state.vy - states[idx - 1].vy)/dt;
            immutable z_acc = (state.vz - states[idx - 1].vz)/dt;

            state.nx = x_acc/g + gravity_vec[0];
            state.ny = y_acc/g + gravity_vec[1];
            state.nz = z_acc/g + gravity_vec[2];

        }

        //Compute load factor
        state.LF = sqrt(state.nx ^^ 2 + state.ny ^^ 2 + state.nz ^^ 2);
        state.nx /= state.LF;
        state.ny /= state.LF;
        state.nz /= state.LF;
        
        //Non-dimensionalize
        state.CT = state.LF * vehicle.weight / (
                vehicle.rho * vehicle.tip_speed ^^ 2 * PI * vehicle.radius ^^ 2);
        auto v = sqrt(state.vx ^^ 2 + state.vy ^^ 2 + state.vz ^^ 2);
        state.mu = v / vehicle.tip_speed;
        auto vn_hat = (state.vx * state.nx + state.vy * state.ny + state.vz * state.nz) / v;
        import std.math : isClose;

        import std.math : abs;
        if((abs(v) > 1.0e-2) && (abs(abs(vn_hat) - 1.0) > 1.0e-2)) {
            state.fpa_e = 180.0 / PI * asin(vn_hat);
            state.fpa = 180.0 / PI * asin(state.vz / v);
        } else {
            state.fpa_e = 0;
            state.fpa = 0;
        }

        state.v = v;
        state.vn_hat = vn_hat;

        trace("Parsed state: ", state);
    }

    return tuple!("states", "utcstart"/+, "world_transform"+/)(states, utcstart/+, world_transform+/);
}

class Trajectory
{
    struct State
    {
        double t;
        double x;
        double y;
        double z;
        double vx;
        double vy;
        double vz;
        double nx;
        double ny;
        double nz;
        double LF;
        double CT;
        double mu;
        double fpa_e;
        double fpa;
        bool day;
        double lat;
        double lon;
        double heading;
        double v;
        double vn_hat;
    }

    struct dimensional_data
    {
        double tip_speed; // m/s
        double radius; // m
        double weight; // kg
        double rho; // kg/m^^3
    }

    void write_state(string state_file)
    {
        import std.conv : to;
        import std.array : join;
        import std.stdio : File;

        info("Writing state data to ", state_file);
        auto file = File(state_file, "w");
        // Write header automatically from field names
        auto fields = __traits(allMembers, State);
        file.writeln([fields].join(","));
        foreach (state; states)
        {
            auto values = state.tupleof;
            file.writeln(to!(string[])([values]).join(","));
        }
    }

    double utcstart;
    State[] states;
    dimensional_data vehicle;

    private WorldTransform world_transform;

    this(string filename, string vehicle_filename, WorldTransform _world_transform, bool stk_format = false)
    {
        import std.exception : enforce;
        //enforce(_world_transform !is null, "Received invalid world transform");
        world_transform = _world_transform;
        //Get vehicle data
        vehicle = configure_vehicle(vehicle_filename);

        auto extension = filename.splitter('.').retro.take(1).array()[0];
        info("Opening file " ~ filename ~ " with extension ", extension);
        switch (toLower(extension))
        {
        case "csv":
            {
                trace("Loading " ~ filename ~ " as CSV.");
                auto contents = readText(filename);
                //auto records = csvReader!(double)(contents, null);

                if(stk_format) {
                    enforce(world_transform !is null, "Did not recieve a valid world transform");
                }

                auto parsed_data = stk_format ? contents.parse_stk(vehicle, world_transform) : contents.parse_greenwood(vehicle);

                utcstart = parsed_data.utcstart;
                states = parsed_data.states;

                break;
            }
        case "nc":
            {
                fatal("netCDF files not yet supported.");
                break;
            }
        default:
            {
                fatal("File extension " ~ extension ~ " not recognized.");
                break;
            }
        }
    }

    void adapt(double dt, int degree = 3, double alpha = 0.3)
    {
        import std.math : sqrt, PI, asin;

        info("Adapting spline to dt ", dt, " with spline of degree ", degree,
                " and smoothing ", alpha);
        assert(dt > 0.0);

        //Unpack states
        double[] t;
        double[] x;
        double[] y;
        double[] z;
        double[] lat;
        double[] lon;
        double[] heading;
        foreach (now; states)
        {
            t ~= now.t;
            x ~= now.x;
            y ~= now.y;
            z ~= now.z;
            lat ~= now.lat;
            lon ~= now.lon;
            heading ~= now.heading;
        }
        trace("t[$ - 1]: ", t[$ - 1]);
        // Build interpolants
        double[] ti;
        ti ~= t[0];
        double tnext = ti[0] + dt;
        while (tnext <= t[$ - 1])
        {
            ti ~= tnext;
            tnext += dt;
        }
        //trace("ti: ", ti);
        trace("Smoothing x.");
        double[] xi;
        double[] vxi;
        double[] axi;
        smoother(t, x, ti, xi, vxi, axi, degree, alpha);
        //trace("ti: ", ti);

        trace("Smoothing y.");
        double[] yi;
        double[] vyi;
        double[] ayi;
        smoother(t, y, ti, yi, vyi, ayi, degree, alpha);

        trace("Smoothing z.");
        double[] zi;
        double[] vzi;
        double[] azi;
        smoother(t, z, ti, zi, vzi, azi, degree, alpha);


        double[] lati;
        double[] vlati;
        double[] alati;
        smoother(t, lat, ti, lati, vlati, alati, degree, alpha);

        double[] loni;
        double[] vloni;
        double[] aloni;
        smoother(t, lon, ti, loni, vloni, aloni, degree, alpha);

        double[] headingi;
        double[] vheadingi;
        double[] aheadingi;
        smoother(t, heading, ti, headingi, vheadingi, aheadingi, degree, alpha);

        import std.algorithm : map;
        import std.array : array;
        import std.typecons : tuple;
        auto day_night = states.map!(a => tuple!("time", "day")(a.t, a.day)).array;
        size_t current_day_idx = 0;

        //Repack data into states
        states = [];
        trace("States reset to ", states);
        foreach (i, time; ti)
        {
            State now;

            if(current_day_idx < day_night.length - 1) {
                if(time >= day_night[current_day_idx + 1].time) {
                    current_day_idx++;
                }
            } else {
                current_day_idx = day_night.length - 1;
            }

            now.t = time;

            now.lat = lati[i];
            now.lon = loni[i];
            
            now.x = xi[i];
            now.y = yi[i];
            now.z = zi[i];
            now.vx = vxi[i];
            now.vy = vyi[i];
            now.vz = vzi[i];

            //auto pos_vec = Vector!(3, double)(now.x, now.y, now.z);
            //auto gravity_vec = pos_vec.normalize;
            auto gravity_vec = Vector!(3, double)(0.0, 0.0, 1.0);

            // Convert to "g"
            now.nx = axi[i] * mps2g + gravity_vec[0];
            now.ny = ayi[i] * mps2g + gravity_vec[1];
            now.nz = azi[i] * mps2g + gravity_vec[2];
            //Compute load facror
            now.LF = sqrt(now.nx ^^ 2 + now.ny ^^ 2 + now.nz ^^ 2);
            //Normalize
            now.nx /= now.LF;
            now.ny /= now.LF;
            now.nz /= now.LF;
            //Non-dimensionalize
            now.CT = now.LF * vehicle.weight / (
                    vehicle.rho * vehicle.tip_speed ^^ 2 * PI * vehicle.radius ^^ 2);
            auto v = sqrt(now.vx ^^ 2 + now.vy ^^ 2 + now.vz ^^ 2);
            now.mu = v / vehicle.tip_speed;
            auto vn_hat = (now.vx * now.nx + now.vy * now.ny + now.vz * now.nz) / v;

            import std.math : isClose;
            //if(!v.isClose(0.0, 1.0e-10)) {
            import std.math : abs;
            //if(abs(v) > 1.0e-5) {
            if((abs(v) > 1.0e-2) && (abs(abs(vn_hat) - 1.0) > 1.0e-2)) {
                now.fpa_e = 180.0 / PI * asin(vn_hat);
                now.fpa = 180.0 / PI * asin(now.vz / v);
            } else {
                now.fpa_e = 0.0;
                now.fpa = 0.0;
            }
            now.day = day_night[current_day_idx].day;
            now.heading = headingi[i];

            now.v = v;
            now.vn_hat = vn_hat;

            states ~= now;
            trace("Adapted state: ", states[$ - 1]);
        }
    }

    template Get(string varname)
    {
        double[] get()
        {
            double[] x;
            x.reserve(states.length);
            foreach (state; states)
            {
                mixin("x ~= state." ~ varname ~ ";");
            }
            return x;
        }
    }

    void watts_trim(double xstart_ft)
    {
        import std.algorithm.searching : minIndex;

        auto mu = Get!("mu").get();
        auto z = Get!("z").get();
        auto imin_z = z.minIndex;
        auto imin_v = mu.minIndex;
        auto ixmax = imin_z < imin_v ? imin_z : imin_v;
        auto x = Get!("x").get();
        size_t ixmin = 0;
        foreach (i, xi; x)
        {
            if (xi > xstart_ft * ft2m)
                break;
            ixmin = i;
        }
        info("Trimming trajectory from x = ", x[ixmin], " to ", x[ixmax]);
        states = states[ixmin .. (ixmax + 1)];
    }
}

static void smoother(in double[] x, in double[] f, in double[] xi, out double[] fi,
        out double[] dfi, out double[] ddfi, int degree = 3, double alpha = 0.3)
{
    import splines : Spline;

    auto spline = Spline(x, f, degree, alpha);
    //Use spline for interpolation
    fi = spline.eval(xi);
    //Use finite differences to get velocity and acceleration
    immutable dx = xi[1] - xi[0];
    diff(dx, fi, dfi, ddfi);
}

static pure void diff(double dx, in double[] f, out double[] df, out double[] ddf)
{
    //Create storage for first derivative
    df = f.dup;
    //Forward diff first point
    df[0] = (f[1] - f[0]) / dx;
    //Central diff interior points
    df[1 .. $ - 1] = 0.5 * (f[2 .. $] - f[0 .. $ - 2]) / dx;
    //Backward diff last point
    df[$ - 1] = (f[$ - 1] - f[$ - 2]) / dx;
    //Create storage for second derivative
    ddf = df.dup;
    //Forward diff first point
    ddf[0] = (f[2] - 2 * f[1] + f[0]) / (dx ^^ 2);
    //Central diff interior points
    ddf[1 .. $ - 1] = (f[2 .. $] - 2 * f[1 .. $ - 1] + f[0 .. $ - 2]) / (dx ^^ 2);
    //Backward diff last point
    ddf[$ - 1] = (f[$ - 1] - 2 * f[$ - 2] + f[$ - 3]) / (dx ^^ 2);
}

static Trajectory.dimensional_data configure_vehicle(string filename)
{
    import std.json : parseJSON;
    import std.file : readText;
    import std.conv : to;

    info("Reading vehicle configuration: ", filename);
    auto config = readText(filename).parseJSON;

    Trajectory.dimensional_data vehicle;
    vehicle.radius = to!double(config["Main Rotor"]["radius"].str);
    vehicle.tip_speed = to!double(config["Main Rotor"]["tip speed"].str);
    vehicle.weight = to!double(config["Vehicle"]["weight"].str) * g; //Actually a mass (kg), so convert to Newtons
    vehicle.rho = to!double(config["Atmosphere"]["density"].str);
    trace("Vehicle data structure: ", vehicle);
    return vehicle;
}
