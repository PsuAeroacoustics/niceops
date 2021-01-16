module splines;

import std.algorithm.iteration : map;
import std.array : array;
import std.experimental.logger;


extern(C) @trusted @nogc nothrow
{
size_t create_spline(size_t nsamples, const double *x_values, const double *y_values, int degree, double alpha);

void destroy_spline(size_t key);

double eval_spline(size_t key, double xi);

double der_spline(size_t key, double xi);

double der2_spline(size_t key, double xi);
}

struct Spline
{
    private:
    size_t spline_key;
    public:
    this(const double[] x, const double[] f, int degree = 3, double alpha = 0.3) @safe @nogc
    {
        spline_key = create_spline(x.length,&x[0],&f[0],degree,alpha);
    }
    ~this() @safe @nogc
    {
        destroy_spline(spline_key);
    }

    double eval(double xi, int order = 0) @safe
    {
        switch(order)
        {
            case 0:
            return eval_spline(spline_key,xi);
            case 1:
            return der_spline(spline_key,xi);
            case 2:
            return der2_spline(spline_key,xi);
            default:
            fatal("Invalid interpolation order = ", order);
        }
        assert(false);
    }
    //Convienience function to apply eval to entire array
    double[] eval(const double[] xi, int order = 0) @safe
    {
        return map!(a => eval(a, order))(xi).array();
    }
}