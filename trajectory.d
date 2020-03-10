module trajectory;

import std.csv;
import std.file : readText;
import std.algorithm.iteration : splitter;
import std.stdio;
import std.range : retro, take;
import std.experimental.logger;
import std.array;
import std.uni : toLower;
import std.math : isNaN, sqrt;
import splines : Spline;

immutable fps2g = 1.0 / 32.17405;

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
    }

    double utcstart;
    State[] states;
    this(string filename)
    {
        auto extension = filename.splitter('.').retro.take(1).array()[0];
        info("Opening file " ~ filename ~ " with extension ", extension);
        switch (toLower(extension))
        {
        case "csv":
            {
                trace("Loading " ~ filename ~ " as CSV.");
                auto contents = readText(filename);
                auto records = csvReader!(double)(contents, null);

                foreach (record; records)
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
                        now.t = data[1];
                        now.x = data[5];
                        now.y = data[6];
                        now.z = data[7];
                        now.vx = data[11];
                        now.vy = data[13];
                        now.vz = data[15];
                        // Convert to "g"
                        now.nx = data[21] * fps2g;
                        now.ny = data[22] * fps2g;
                        now.nz = (data[23] * fps2g + 1.0);
                        //Normalize
                        auto magnitude = 1.0 / sqrt(now.nx ^^ 2 + now.ny ^^ 2 + now.nz ^^ 2);
                        now.nx *= magnitude;
                        now.ny *= magnitude;
                        now.nz *= magnitude;
                        states ~= now;
                        trace("Parsed state: ", states[$ - 1]);
                    }
                    else
                    {
                        fatal("Wrong number of columns.");
                    }
                }
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
            //TODO: Transform accelerations to normal vector [check units!]
            enum g = 32.2;
            now.nx = axi[i] / g;
            now.ny = ayi[i] / g;
            now.nz = azi[i] / g + 1;
            states ~= now;
        }
    }
}

static void smoother(in double[] x, in double[] f, in double[] xi, out double[] fi,
        out double[] dfi, out double[] ddfi, int degree = 3, double alpha = 0.3)
{
    auto spline = Spline(x, f, degree, alpha);
    fi = spline.eval(xi, 0);
    dfi = spline.eval(xi, 1);
    ddfi = spline.eval(xi, 2);
}
