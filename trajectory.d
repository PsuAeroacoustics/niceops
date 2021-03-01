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

    return tuple!("states", "utcstart", "world_transform")(states, utcstart, new WorldTransform);
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

auto parse_stk(string file_contents, Trajectory.dimensional_data vehicle) {
    import std.algorithm : filter, map;
    import std.conv : to;
    import std.range : enumerate, slide;
    import std.stdio : writeln;
    import std.typecons : tuple;

    double utcstart;
    Trajectory.State[] states;

    struct Layout {
        string utc_timestamp;
        /+
            "Time (UTCG)",
            "Lat (deg)",
            "Lon (deg)",
            "Alt (m)",
            "Lat Rate (deg/sec)",
            "Lon Rate (deg/sec)",
            "Alt Rate (m/sec)",
            "Yaw (deg)",
            "Pitch (deg)",
            "Roll (deg)",
            "Bank (deg)",
            "Heading (deg)",
            "Acceleration (m/sec^2)"
        +/
        double lat;
        double lon;
        double alt;
        double lat_rate;
        double lon_rate;
        double alt_rate;
        //double velocity;
        double yaw;
        double pitch;
        double roll;
        double bank;
        double heading;
        //double x_acc;
        //double y_acc;
        //double z_acc;
        double acc_mag;
    }

    auto records = file_contents.csvReader!Layout(null);
    
    auto start_record = records.front;

    utcstart = start_record.utc_timestamp.time_str_to_double;

    auto lat_long_alt_reference = Vector!(3, double)(
        start_record.lat*deg2rad,
        start_record.lon*deg2rad,
        start_record.alt
    );

    auto world_transform = new WorldTransform(lat_long_alt_reference);

    states = records
        .enumerate
        .map!((record) {
            Trajectory.State now;
            now.t = record.value.utc_timestamp.time_str_to_double - utcstart;

            auto lat_long_alt_point = Vector!(3, double)(
                record.value.lat*deg2rad,
                record.value.lon*deg2rad,
                record.value.alt
            );

            auto lat_long_alt_rate = Vector!(3, double)(
                record.value.lat_rate*deg2rad,
                record.value.lon_rate*deg2rad,
                record.value.alt_rate
            );

            auto point = world_transform.lat_long_alt_to_enu(lat_long_alt_point);

            now.x = point[0];
            now.y = point[1];
            now.z = point[2];

            auto enu_vel = world_transform.lat_long_alt_rate_to_enu_rate(lat_long_alt_point, lat_long_alt_rate);

            now.vx = enu_vel[0];
            now.vy = enu_vel[1];
            now.vz = enu_vel[2];

            return now;
        })
        .array;

    // Compute accelerations and other derived quantities.
    foreach(idx, ref state; states) {
        if(idx == 0) {
            state.nx = 0;
            state.ny = 0;
            state.nz = 1.0;
        } else {
            immutable dt = state.t - states[idx - 1].t;
            immutable x_acc = (state.vx - states[idx - 1].vx)/dt;
            immutable y_acc = (state.vy - states[idx - 1].vy)/dt;
            immutable z_acc = (state.vz - states[idx - 1].vz)/dt;

            state.nx = x_acc/g;
            state.ny = y_acc/g;
            state.nz = z_acc/g + 1.0;

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
        state.fpa_e = 180.0 / PI * asin(vn_hat);
        state.fpa = 180.0 / PI * asin(state.vz / v);

        trace("Parsed state: ", state);
    }

    return tuple!("states", "utcstart", "world_transform")(states, utcstart, world_transform);
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
        import std.conv;
        import std.array;
        import std.file;
        info( "Writing state data to ", state_file);
        auto file = File( state_file, "w");
        // Write header automatically from field names
        auto fields = __traits(allMembers, State);
        file.writeln( [ fields].join( ","));
        foreach (state; states)
        {
            auto values = state.tupleof;
            file.writeln( to!(string[])( [ values]).join( ","));
        }
    }

    double utcstart;
    State[] states;
    dimensional_data vehicle;

    WorldTransform world_transform;

    this(string filename, string vehicle_filename, bool stk_format = false)
    {
        //Get vehicle data
        vehicle = configure_vehicle( vehicle_filename);

        auto extension = filename.splitter( '.').retro.take( 1).array()[0];
        info("Opening file " ~ filename ~ " with extension ", extension);
        switch (toLower(extension))
        {
        case "csv":
            {
                trace("Loading " ~ filename ~ " as CSV.");
                auto contents = readText(filename);
                //auto records = csvReader!(double)(contents, null);

                auto parsed_data = stk_format ? contents.parse_stk(vehicle) : contents.parse_greenwood(vehicle);

                utcstart = parsed_data.utcstart;
                states = parsed_data.states;
                world_transform = parsed_data.world_transform;

                /+foreach (record; records)
                {
                    auto data = record.array;
                    if (data.length == 27)
                    {
                        //Assume NASA/FAA AC state data format
                        State now;
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
                        now.fpa_e = 180.0 / PI * asin(vn_hat);
                        now.fpa = 180.0 / PI * asin(now.vz / v);
                        states ~= now;
                        trace("Parsed state: ", states[$ - 1]);
                    }
                    else
                    {
                        fatal("Wrong number of columns.");
                    }
                }+/
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

    /+this(string aircraft_vel_filename, string flightpath_filename, string vehicle_filename)
    {
        import std.array : staticArray, array;
        import std.algorithm : filter, map, strip;
        import std.conv : to;
        import std.exception : enforce;
        import std.range : dropOne, dropExactly, enumerate, takeOne, tee, walkLength, zip;
        import std.string : isNumeric;
        import std.typecons : Tuple, tuple;
        import std.stdio : writeln;

        alias time_vel = tuple!("t", "vel");

        auto velocity_data =
            aircraft_vel_filename
                .readText
                .split("\n")
                .map!(a => a.strip(' '))
                .filter!(a => a.length > 0 && a[0].to!string.isNumeric)
                .map!((a) {
                    auto values = a.split(" ").filter!(b => b != "").map!(b => b.to!double);
                    return time_vel(values.takeOne.front, values.dropExactly(4).staticArray!3);
                });

        auto pos_data =
            flightpath_filename
                .readText
                .split("\n")
                .map!(a => a.strip(' '))
                .filter!(a => a.length > 0 && a[0].to!string.isNumeric) // Drop all lines that have no data or do not start with a number
                .dropOne // Drop first line as its the timestamp at the top.
                .map!(a =>
                    a
                    .split(" ")
                    .filter!(b => b != "")
                    .map!(b => b.strip(' '))
                    .dropExactly(4)
                    .map!(b => b.to!double).staticArray!3
                );

        enforce(velocity_data.walkLength == pos_data.walkLength, "Number of data samples in velocity file and position file do match");

        states = zip(pos_data, velocity_data)
            .map!((a) {
                State now = {
                    t: a[1].t,
                    x: a[0][0],
                    y: a[0][1],
                    z: a[0][2],
                    vx: a[1].vel[0],
                    vy: a[1].vel[1],
                    vz: a[1].vel[2]

                };
                return now;
            }).array;

        //writeln(states);
    }+/

    void adapt(double dt, int degree = 3, double alpha = 0.3)
    {
        info("Adapting spline to dt ", dt, " with spline of degree ", degree,
                " and smoothing ", alpha);
        assert(dt > 0.0);

        //Unpack states
        double[] t;
        double[] x;
        double[] y;
        double[] z;
        foreach (now; states)
        {
            t ~= now.t;
            x ~= now.x;
            y ~= now.y;
            z ~= now.z;
        }
        // Build interpolants
        double[] ti;
        ti ~= t[0];
        double tnext = ti[0] + dt;
        while (tnext <= t[$ - 1])
        {
            ti ~= tnext;
            tnext += dt;
        }
        trace("Smoothing x.");
        double[] xi;
        double[] vxi;
        double[] axi;
        smoother(t, x, ti, xi, vxi, axi, degree, alpha);

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

        //Repack data into states
        states = [];
        trace("States reset to ", states);
        foreach (i, time; ti)
        {
            State now;
            now.t = time;
            now.x = xi[i];
            now.y = yi[i];
            now.z = zi[i];
            now.vx = vxi[i];
            now.vy = vyi[i];
            now.vz = vzi[i];
            // Convert to "g"
            now.nx = axi[i] * mps2g;
            now.ny = ayi[i] * mps2g;
            now.nz = azi[i] * mps2g + 1.0;
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
            now.fpa_e = 180.0 / PI * asin(vn_hat);
            now.fpa = 180.0 / PI * asin(now.vz / v);
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
        import std.algorithm.searching;

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
    auto spline = Spline(x, f, degree, alpha);
    //Use spline for interpolation
    fi = spline.eval(xi);
    //Use finite differences to get velocity and acceleration
    double dx = xi[1] - xi[0];
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
