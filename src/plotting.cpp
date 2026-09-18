#include "plotting.hpp"

#include <pybind11/embed.h>
// Supplies the std::vector <-> Python conversions.  Without it every call
// taking a vector -- which is most of them -- throws pybind11::cast_error at
// run time while still compiling cleanly, so the whole plotting path aborted.
#include <pybind11/stl.h>

namespace py = pybind11;

namespace plt {

namespace {

// Lazily start the embedded interpreter and import pyplot with a non-interactive
// backend.  The first call pays the import cost; non-plotting runs never do.
py::module_& pyplot()
{
    static py::scoped_interpreter guard {};
    static py::module_ module = [] {
        py::module_::import("matplotlib").attr("use")("Agg");
        return py::module_::import("matplotlib.pyplot");
    }();
    return module;
}

void label(const char* attr, const std::string& text, int fontsize)
{
    if (fontsize > 0)
    {
        pyplot().attr(attr)(text, py::arg("fontsize") = fontsize);
    }
    else
    {
        pyplot().attr(attr)(text);
    }
}

} // namespace

void figure()
{
    pyplot().attr("figure")();
}

void contourf(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<std::vector<double>>& z, const std::vector<double>& levels,
        const std::string& cmap)
{
    pyplot().attr("contourf")(x, y, z, levels, py::arg("cmap") = cmap);
}

void tricontourf(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<double>& z, const std::vector<double>& levels, const std::string& cmap)
{
    pyplot().attr("tricontourf")(x, y, z, levels, py::arg("cmap") = cmap);
}

void plot(const std::vector<double>& x, const std::vector<double>& y)
{
    pyplot().attr("plot")(x, y);
}

void plot(const std::vector<double>& x, const std::vector<double>& y, const std::string& format)
{
    pyplot().attr("plot")(x, y, format);
}

void plot(const std::vector<double>& x1, const std::vector<double>& y1,
        const std::vector<double>& x2, const std::vector<double>& y2)
{
    pyplot().attr("plot")(x1, y1, x2, y2);
}

void xlabel(const std::string& text, int fontsize)
{
    label("xlabel", text, fontsize);
}

void ylabel(const std::string& text, int fontsize)
{
    label("ylabel", text, fontsize);
}

void axis(const std::string& what)
{
    pyplot().attr("axis")(what);
}

void colorbar()
{
    pyplot().attr("colorbar")();
}

void tight_layout()
{
    pyplot().attr("tight_layout")();
}

void savefig(const std::string& filename)
{
    pyplot().attr("savefig")(filename);
}

void show()
{
    pyplot().attr("show")();
}

void clear()
{
    pyplot().attr("close")("all");
}

} // namespace plt
