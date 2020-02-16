module splines;

import std.functional : memoize;
import std.algorithm.iteration : map;
import std.algorithm.searching : until;
import std.array : array;
import std.stdio;

// Smoothing spline functions:
extern (C)
{
    void cubgcv(double* x, double* f, double* df, int* n, double* y, double* c,
    int* ic, double* var, int* job, double* se, double* wk, int* ier);
    void b_spline(int n, double* x, double* y, int m, double* sx, double* sy);
}

private size_t find_index(T)(T[] x, T xi)
{
    auto range = x.until!((a, b) => a > b)( xi);
    return range.array.length - 1;
}

struct Spline
{
    double[] c;
    double[] y;
    double[] x;
    int ic;
    double varience;
    this(double[] x_, double[] f, double var = -1.0)
    {
        varience = var;
        x = x_.dup;
        int n = cast(int) x.length;
        auto df = f.dup;
        df[] = 1.0;
        y = new double[n];
        c = new double[3 * (n - 1)];
        ic = n - 1;
        int job = 0;
        double* se = null;
        double[] wk = new double[7 * (n + 2)];
        int ier = 0;
        cubgcv( x.ptr, f.ptr, df.ptr, &n, y.ptr, c.ptr, &ic, &varience,
        &job, se, wk.ptr, &ier);
    }

    double eval(double xi, int order = 0)
    {
        auto index = finder( x, xi);
        auto xi0 = x[index];
        auto c0 = c[ic * 0 + index];
        auto c1 = c[ic * 1 + index];
        auto c2 = c[ic * 2 + index];
        auto d = xi - xi0;
        assert(order >= 0 && order < 3);
        switch (order)
        {
            default:
            {
                return (c2 * d + c1) * d + c0 + y[index];
            }
            case 1:
            {
                return 2 * c2 * d + c1;
            }
            case 2:
            {
                return 2 * c2;
            }
        }
    }

    //Convienience function to apply eval to entire array
    double[] eval(double[] xi, int order = 0)
    {
        return map!(a => eval( a, order))( xi).array();
    }

    // Instantiate find_index as double, and memoize
    private:
    static auto find_d = (double[] x, double xi) => find_index( x, xi);
    alias finder = memoize!find_d;
}
