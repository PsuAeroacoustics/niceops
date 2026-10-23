#pragma once

// The parabolic equation's options, shared by niceops_pe_compare and
// niceops_pe_field so that the two draw the same field from the same command
// line: the ground, the march and its grid, turbulence, and the threads its
// bands run on.  Each tool adds its own geometry -- the source, the receivers,
// the range -- to the setups these give.

#include <CLI/CLI.hpp>

#include <cstddef>
#include <string>

#include "ground_effect.hpp"
#include "pe.hpp"

namespace niceops {

struct PeOptions
{
    double flow_resistance = flow_resistance::grassy_field;
    std::string impedance_model = "delany_bazley";
    double porosity_rate = 0.0;
    int samples = 5;
    double points_per_wavelength = 0.0;
    double range_points_per_wavelength = 10.0;
    std::string scheme = "pade";
    double range_step_wavelengths = 8.0;
    double max_angle = 50.0;
    int vertical_order = 4;
    int boundary_order = 4;
    std::string starter = "self";
    double gamma_t = 0.0;
    double outer_scale = 10.0;
    unsigned seed = 1;
    int threads = 0;

    void add_to(CLI::App& app)
    {
        app.add_option("--flow_resistance", flow_resistance, "Ground flow resistance, kPa s/m^2.");
        app.add_option("--impedance_model", impedance_model, "As for niceops.");
        app.add_option("--porosity_rate", porosity_rate, "As for niceops, 1/m.");
        app.add_option("--samples", samples, "PE frequencies per band (default 5).")->check(CLI::PositiveNumber);
        app.add_option("--ppw", points_per_wavelength,
                   "PE grid points per wavelength in height (default 10 where the grid is fourth order "
                   "throughout, 20 under --scheme cn or either order 2).")
                ->check(CLI::NonNegativeNumber);
        app.add_option("--scheme", scheme,
                   "How the PE marches in range: pade (default), the split-step Pade march, or cn, "
                   "Crank-Nicolson.")
                ->check(CLI::IsMember({"pade", "cn"}));
        app.add_option("--range_step", range_step_wavelengths,
                   "--scheme pade: longest range step, in wavelengths (default 8).")
                ->check(CLI::PositiveNumber);
        app.add_option("--max_angle", max_angle,
                   "--scheme pade: elevation angle the march is accurate to, degrees (default 50, at most 60).  "
                   "Steeper paths from the source come out as NaN.")
                ->check(CLI::Range(0.1, 60.0));
        app.add_option("--vertical_order", vertical_order,
                   "--scheme pade: order of the height derivative on the grid, 4 (default) or 2.")
                ->check(CLI::IsMember({2, 4}));
        app.add_option("--boundary_order", boundary_order,
                   "--scheme pade: order of the ground condition and of the source's and receivers' "
                   "interpolation, 4 (default) or 2.")
                ->check(CLI::IsMember({2, 4}));
        app.add_option("--starter", starter,
                   "--scheme pade: self (default), the wide-angle self-starter, or gaussian.")
                ->check(CLI::IsMember({"self", "gaussian"}));
        app.add_option("--range_ppw", range_points_per_wavelength,
                   "--scheme cn: range steps per wavelength (default 10).")
                ->check(CLI::PositiveNumber);
        app.add_option("--gamma_t", gamma_t,
                "Turbulence in the PE (HARMONOISE strength; 0, the default, is none).  Each of --samples "
                "frequencies is then its own realization, so give it 20 or so.");
        app.add_option("--outer_scale", outer_scale,
                "von Karman outer scale of the PE's turbulence, m (default 10).");
        app.add_option("--seed", seed, "First realization's seed (default 1).");
        app.add_option("-j,--threads", threads,
                   "Worker threads for the PE's bands (default 0: NICEOPS_THREADS when it is set, else one "
                   "per hardware thread).  The output does not depend on it.")
                ->check(CLI::NonNegativeNumber);
    }

    // The ground they describe; the caller adds a receiver height where its
    // ground model wants one.
    GroundModel ground() const
    {
        GroundModel g;
        g.flow_resistance = flow_resistance;
        g.impedance_model = parse_impedance_model(impedance_model);
        g.porosity_rate = porosity_rate;
        return g;
    }

    // Band `band`'s march on `workers` threads, without its geometry: the
    // caller sets the source, the receivers, the range, the profile and the
    // admittance.
    PeSetup setup(std::size_t band, unsigned workers) const
    {
        PeSetup s;
        s.points_per_wavelength = points_per_wavelength;
        s.range_points_per_wavelength = range_points_per_wavelength;
        s.scheme = scheme == "cn" ? PeScheme::crank_nicolson : PeScheme::split_step_pade;
        s.range_step_wavelengths = range_step_wavelengths;
        s.max_angle = max_angle;
        s.vertical_order = vertical_order;
        s.boundary_order = boundary_order;
        s.starter = starter == "gaussian" ? PeStarter::gaussian : PeStarter::self;
        s.gamma_t = gamma_t;
        s.outer_scale = outer_scale;
        s.seed = seed + 1000u * static_cast<unsigned>(band);
        s.threads = workers;
        return s;
    }
};

} // namespace niceops
