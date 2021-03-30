#include <cstddef>
#include <unordered_map>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>

using namespace SPLINTER;


struct spline_interface {
public:
    BSpline *spline;

    spline_interface(size_t nsamples, const double *x_values, const double *y_values, int degree = 3, double alpha = 0.3) {
        // Set smoothing mode
        auto smooth_mode = BSpline::Smoothing::PSPLINE;
        // Switch off smoothing for invalid alpha
        if (alpha < 0.0) {
            smooth_mode = BSpline::Smoothing::NONE;
        }
        //Create data table and add samples to it
        DataTable samples;
        for (size_t i = 0; i < nsamples; i++) {
            samples.addSample(x_values[i], y_values[i]);
        }
        // Build B-spline into pointer
        spline = new BSpline(BSpline::Builder(samples).degree(degree).smoothing(smooth_mode).alpha(alpha).build());
    }

    ~spline_interface() {
        delete spline;
    }

    double eval(double xi) {
        return spline->eval(std::vector<double>{xi});
    }

    double derivative(double xi) {
        return spline->evalJacobian(std::vector<double>{xi})[0];
    }

    double derivative2(double xi) {
        return spline->evalHessian(std::vector<double>{xi})[0][0];
    }

    static spline_interface *
    make_spline(size_t nsamples, double *x_values, double *y_values, int degree = 3, double alpha = 0.3) {
        return new spline_interface(nsamples, x_values, y_values, degree, alpha);
    }

    static void destroy_spline(spline_interface *x_spline) {
        delete x_spline;
    }
};

#include <atomic>
#include <thread>
#include <mutex>
static std::atomic<size_t> number_of_splines = 0;
static std::unordered_map<size_t, spline_interface *> splines;
static std::mutex spline_mutex;
extern "C"
{
size_t create_spline(size_t nsamples, const double *x_values, const double *y_values, int degree, double alpha) {
    std::scoped_lock lock(spline_mutex);
    number_of_splines++;
    splines[number_of_splines] = new spline_interface(nsamples, x_values, y_values, degree, alpha);
    return number_of_splines;
}

void destroy_spline(size_t key) {
    std::scoped_lock lock(spline_mutex);
    delete splines[key];
    splines.erase(key);
}

double eval_spline(size_t key, double xi) {
    return splines[key]->eval(xi);
}

double der_spline(size_t key, double xi) {
    return splines[key]->derivative(xi);
}

double der2_spline(size_t key, double xi) {
    return splines[key]->derivative2(xi);
}
}
