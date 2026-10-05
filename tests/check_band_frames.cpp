// Checks what --band_frames and --epnl_file wrote, for the CTest cases that run
// them (CMakeLists.txt, epnl_*):
//
//   check_band_frames FRAMES.nc EPNL.csv [FULL_HISTORY_LAE.csv]
//
//   * the netCDF file has its dimensions, coordinates and attributes, and one
//     observer per EPNL row, at the same microphone;
//   * every observer's EPNL and PNLTM, recomputed by src/pnl.cpp from the
//     frames in the netCDF file, are what the CSV says (to the float
//     rounding of the stored levels), and EPNL = PNLTM + D;
//   * the values are plausible: EPNL between 40 and 130 EPNdB, D within
//     -15..+15 dB, the tone correction within 0..20/3 dB;
//   * given a --full_history SEL file from the same run (scattered observers),
//     the A-weighted energy in the frames is that SEL: the frames hold every
//     state's energy, in reception time, exactly as the SEL integral does.
//
// Exits nonzero with the first failure on stderr.

#include <netcdf.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "pnl.hpp"
#include "util/csv.hpp"
#include "weighting.hpp"

using namespace niceops;

namespace {

[[noreturn]] void fail(const std::string& what)
{
    std::cerr << "check_band_frames: " << what << "\n";
    std::exit(1);
}

void check(int status, const std::string& what)
{
    if (status != NC_NOERR)
        fail(what + ": " + nc_strerror(status));
}

std::size_t dimension(int ncid, const char* name)
{
    int id = 0;
    std::size_t length = 0;
    check(nc_inq_dimid(ncid, name, &id), std::string("dimension ") + name);
    check(nc_inq_dimlen(ncid, id, &length), std::string("dimension ") + name);
    return length;
}

std::vector<double> read_doubles(int ncid, const char* name, std::size_t count)
{
    int var = 0;
    check(nc_inq_varid(ncid, name, &var), std::string("variable ") + name);
    std::vector<double> values(count);
    if (count > 0)
        check(nc_get_var_double(ncid, var, values.data()), std::string("reading ") + name);
    return values;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 4)
        fail("usage: check_band_frames FRAMES.nc EPNL.csv [FULL_HISTORY_LAE.csv]");
    int ncid = 0;
    check(nc_open(argv[1], NC_NOWRITE, &ncid), std::string("opening ") + argv[1]);
    const std::size_t observers = dimension(ncid, "observer");
    const std::size_t frames = dimension(ncid, "frame");
    const std::size_t bands = dimension(ncid, "band");
    if (observers == 0 || frames == 0 || bands == 0)
        fail("an empty dimension");
    int version = 0;
    check(nc_get_att_int(ncid, NC_GLOBAL, "band_frames_version", &version), "band_frames_version");
    double frame_length = 0.0;
    check(nc_get_att_double(ncid, NC_GLOBAL, "frame_length_seconds", &frame_length), "frame_length_seconds");
    if (version != 1 || frame_length != pnl_record_seconds)
        fail("unexpected version or frame length");
    const auto time = read_doubles(ncid, "time", frames);
    for (std::size_t k = 1; k < frames; k++)
        if (std::fabs(time[k] - time[k - 1] - frame_length) > 1e-9)
            fail("the frame axis is not evenly spaced by the frame length");
    // Frame centers sit at (k + 1/2) F: frames line up on the trajectory's clock.
    const double phase = time[0] / frame_length - 0.5;
    if (std::fabs(phase - std::round(phase)) > 1e-9)
        fail("frames are not aligned on multiples of the frame length");
    const auto centers = read_doubles(ncid, "band_center", bands);
    const auto x = read_doubles(ncid, "x", observers);
    const auto y = read_doubles(ncid, "y", observers);
    const auto z = read_doubles(ncid, "z", observers);
    read_doubles(ncid, "microphone_height", observers);
    std::vector<float> level(observers * frames * bands);
    int level_var = 0;
    check(nc_inq_varid(ncid, "level", &level_var), "variable level");
    check(nc_get_var_float(ncid, level_var, level.data()), "reading level");
    nc_close(ncid);

    std::vector<double> missing;
    const auto pnl_bands = pnl_band_indices(centers, missing);
    if (!missing.empty())
        fail("the frames lack PNL bands");

    // observer,x,y,z,EPNL,PNLTM,D,C,tone band,adjusted,t_PNLTM,t1,t2,clipped
    const auto rows = read_numeric_csv(argv[2], true);
    if (rows.size() != observers)
        fail("the EPNL file has " + std::to_string(rows.size()) + " rows for " + std::to_string(observers)
                + " observers");
    for (std::size_t i = 0; i < observers; i++)
    {
        const auto& row = rows[i];
        const std::string who = "observer " + std::to_string(i) + ": ";
        if (row.size() != 14 || row[0] != static_cast<double>(i))
            fail(who + "malformed EPNL row");
        if (std::fabs(row[1] - x[i]) > 1e-6 * std::max(1.0, std::fabs(x[i]))
                || std::fabs(row[2] - y[i]) > 1e-6 * std::max(1.0, std::fabs(y[i]))
                || std::fabs(row[3] - z[i]) > 1e-6 * std::max(1.0, std::fabs(z[i])))
            fail(who + "the two files disagree on the microphone");
        // This observer's frames: the span where the levels are not NaN.
        std::vector<PnlSpectrum> history;
        bool inside = false, ended = false;
        for (std::size_t k = 0; k < frames; k++)
        {
            const float* f = &level[(i * frames + k) * bands];
            const bool present = !std::isnan(f[0]);
            if (present && ended)
                fail(who + "a gap inside the received frames");
            if (!present && inside)
                ended = true;
            inside = inside || present;
            if (!present)
                continue;
            PnlSpectrum s;
            for (std::size_t b = 0; b < pnl_band_count; b++)
                s[b] = f[pnl_bands[b]];
            history.push_back(s);
        }
        const double epnl = row[4], pnltm = row[5], d = row[6], tone = row[7];
        if (history.empty())
            fail(who + "no frames");
        const EpnlResult r = effective_perceived_noise_level(history, frame_length);
        if (std::fabs(r.epnl - epnl) > 1e-3 || std::fabs(r.pnltm - pnltm) > 1e-3)
            fail(who + "EPNL " + std::to_string(epnl) + " / PNLTM " + std::to_string(pnltm)
                    + " recompute from the frames as " + std::to_string(r.epnl) + " / "
                    + std::to_string(r.pnltm));
        if (std::fabs(epnl - (pnltm + d)) > 1e-6)
            fail(who + "EPNL is not PNLTM + D");
        if (!(epnl > 40.0 && epnl < 130.0) || !(d > -15.0 && d < 15.0)
                || !(tone >= 0.0 && tone <= 20.0 / 3.0 + 1e-9))
            fail(who + "implausible EPNL " + std::to_string(epnl) + ", D " + std::to_string(d)
                    + ", tone correction " + std::to_string(tone));
        if (!(row[11] <= row[10] && row[10] <= row[12]))
            fail(who + "PNLTM is not inside its duration");
    }

    if (argc == 4)
    {
        // x, y, SEL per observer, from a --full_history run of the same case.
        const auto sel = read_numeric_csv(argv[3], true);
        if (sel.size() != observers)
            fail("the SEL file does not have one row per observer");
        const auto gains = a_weighting_gains(centers, BandLimits {});
        for (std::size_t i = 0; i < observers; i++)
        {
            double energy = 0.0;
            for (std::size_t k = 0; k < frames; k++)
            {
                const float* f = &level[(i * frames + k) * bands];
                if (std::isnan(f[0]))
                    continue;
                for (std::size_t b = 0; b < bands; b++)
                    energy += gains[b] * std::pow(10.0, f[b] / 10.0) * frame_length;
            }
            const double framed = 10.0 * std::log10(energy);
            // Float levels: 1e-5 dB a band at most.
            if (std::fabs(framed - sel[i][2]) > 1e-4)
                fail("observer " + std::to_string(i) + ": the frames hold " + std::to_string(framed)
                        + " dB of A-weighted exposure, the whole-history SEL is " + std::to_string(sel[i][2]));
        }
    }
    std::cout << "check_band_frames: " << observers << " observers, " << frames << " frames, " << bands
              << " bands: consistent.\n";
    return 0;
}
