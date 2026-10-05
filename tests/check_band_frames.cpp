// Checks what --band_frames and --epnl_file wrote, for the CTest cases that run
// them (CMakeLists.txt, epnl_*):
//
//   check_band_frames [--no-masking] FRAMES.nc EPNL.csv [FULL_HISTORY_LAE.csv]
//   check_band_frames --history PREFIX FRAMES.nc
//   check_band_frames --compare MASKED_EPNL.csv UNMASKED_EPNL.csv
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
//     --no-masking: the EPNL file was written with --no-epnl_masking.  Without
//     it, the band each PNLTM record's tone came from must be at or above its
//     noy threshold (pnl.hpp's masking), so the tone is audible;
//   * with --compare, the same case with and without --epnl_masking: on a
//     case whose top bands fall below their thresholds, some EPNL must fall by
//     more than 0.1 dB, or the option is not reaching the tone correction.
//     (Masking can also raise an EPNL a little -- the band-sharing check and
//     the duration can pick other records, and a masked band taken at its
//     threshold moves its neighbors' background -- so a rise is reported, not
//     refused.)
//   * with --history, the frames are rebuilt from a --band_history run of the
//     same case (PREFIX_i<N>.csv): each state's band levels spread over its
//     stretch of reception time, centered on the file's own reception_time
//     column and as wide as that column's step.  The energy check above cannot
//     see where in time the energy went; this can (a state filed at its
//     emission time, or at t - r/c, lands in other frames).
//
// Exits nonzero with the first failure on stderr.

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "band_frames.hpp"
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

// One --band_history file: the header's column names, and the rows, empty
// fields as NaN.
struct HistoryFile
{
    std::vector<std::string> columns;
    std::vector<std::vector<double>> rows;
};

HistoryFile read_history(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        fail("cannot open " + path);
    HistoryFile h;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ','))
            fields.push_back(field);
        if (!line.empty() && line.back() == ',')
            fields.emplace_back();
        if (h.columns.empty())
        {
            h.columns = fields;
            continue;
        }
        if (fields.size() != h.columns.size())
            fail(path + ": a row of " + std::to_string(fields.size()) + " fields under a header of "
                    + std::to_string(h.columns.size()));
        std::vector<double> row;
        for (const auto& f : fields)
            row.push_back(f.empty() ? std::nan("") : std::stod(f));
        h.rows.push_back(row);
    }
    return h;
}

// --history PREFIX FRAMES.nc
int check_against_history(const std::string& prefix, const std::string& frames_path)
{
    int ncid = 0;
    check(nc_open(frames_path.c_str(), NC_NOWRITE, &ncid), "opening " + frames_path);
    const std::size_t observers = dimension(ncid, "observer");
    const std::size_t frames = dimension(ncid, "frame");
    const std::size_t bands = dimension(ncid, "band");
    const auto time = read_doubles(ncid, "time", frames);
    const auto centers = read_doubles(ncid, "band_center", bands);
    std::vector<float> level(observers * frames * bands);
    int level_var = 0;
    check(nc_inq_varid(ncid, "level", &level_var), "variable level");
    check(nc_get_var_float(ncid, level_var, level.data()), "reading level");
    nc_close(ncid);
    if (observers == 0 || frames == 0 || bands == 0)
        fail("an empty dimension");
    const auto first_frame = static_cast<std::int64_t>(std::llround(time[0] / pnl_record_seconds - 0.5));

    std::size_t compared = 0;
    for (std::size_t i = 0; i < observers; i++)
    {
        const std::string path = prefix + "_i" + std::to_string(i) + ".csv";
        const HistoryFile h = read_history(path);
        auto column = [&](const std::string& name) {
            for (std::size_t c = 0; c < h.columns.size(); c++)
                if (h.columns[c] == name)
                    return c;
            fail(path + " has no column " + name);
        };
        const std::size_t reception = column("reception_time");
        std::vector<std::size_t> band_column(bands);
        for (std::size_t b = 0; b < bands; b++)
        {
            std::ostringstream name;
            name << "band_" << centers[b];
            band_column[b] = column(name.str());
        }
        const std::size_t nt = h.rows.size();
        if (nt < 2)
            fail(path + ": fewer than two states");
        std::vector<double> levels(nt * bands);
        for (std::size_t it = 0; it < nt; it++)
            for (std::size_t b = 0; b < bands; b++)
                levels[it * bands + b] = h.rows[it][band_column[b]];
        // The stretch a state covers is the step of reception time around it,
        // (t_r(it+1) - t_r(it-1)) / 2, one-sided at the ends.
        auto t_r = [&](std::size_t it) { return h.rows[it][reception]; };
        BandFrames rebuilt(bands, pnl_record_seconds);
        rebuilt.accumulate(
                nt, [&](std::size_t it) { return std::isfinite(t_r(it)); }, t_r,
                [&](std::size_t it) {
                    const std::size_t lo = it > 0 ? it - 1 : it;
                    const std::size_t hi = it + 1 < nt ? it + 1 : it;
                    return (t_r(hi) - t_r(lo)) / static_cast<double>(hi - lo);
                },
                [&](std::size_t it) { return levels.data() + it * bands; });
        // The file's own frames: the span that is not NaN.
        std::size_t start = frames, stop = 0;
        for (std::size_t k = 0; k < frames; k++)
        {
            if (!std::isnan(level[(i * frames + k) * bands]))
            {
                start = std::min(start, k);
                stop = k + 1;
            }
        }
        const std::string who = "observer " + std::to_string(i) + ": ";
        if (start == frames || rebuilt.count() != stop - start
                || rebuilt.first() != first_frame + static_cast<std::int64_t>(start))
            fail(who + "the file's frames are " + std::to_string(stop - start) + " from frame "
                    + std::to_string(first_frame + static_cast<std::int64_t>(start)) + ", the history's "
                    + std::to_string(rebuilt.count()) + " from " + std::to_string(rebuilt.first()));
        // The first and last frames hold a sliver of a state, where the
        // history's ten-digit reception times move the split; the rest must
        // agree to the float and the CSV's rounding.
        for (std::size_t k = 1; k + 1 < rebuilt.count(); k++)
        {
            for (std::size_t b = 0; b < bands; b++)
            {
                const double expected = rebuilt.level(k, b);
                const double got = level[(i * frames + start + k) * bands + b];
                if (!(std::fabs(got - expected) <= 1e-3 + 1e-6 * std::fabs(expected)))
                    fail(who + "frame " + std::to_string(k) + ", " + std::to_string(centers[b])
                            + " Hz: the file has " + std::to_string(got) + " dB, the history gives "
                            + std::to_string(expected));
                compared++;
            }
        }
    }
    std::cout << "check_band_frames: " << observers << " observers' frames match their band histories ("
              << compared << " frame-bands).\n";
    return 0;
}

// --compare MASKED.csv UNMASKED.csv
int compare_masking(const std::string& masked_path, const std::string& unmasked_path)
{
    const auto masked = read_numeric_csv(masked_path, true);
    const auto unmasked = read_numeric_csv(unmasked_path, true);
    if (masked.empty() || masked.size() != unmasked.size())
        fail("the two EPNL files do not have the same, nonzero number of rows");
    double largest_drop = 0.0;
    double largest_rise = 0.0;
    for (std::size_t i = 0; i < masked.size(); i++)
    {
        const std::string who = "observer " + std::to_string(i) + ": ";
        if (masked[i].size() != 14 || unmasked[i].size() != 14 || masked[i][0] != unmasked[i][0])
            fail(who + "malformed EPNL row");
        largest_drop = std::max(largest_drop, unmasked[i][4] - masked[i][4]);
        largest_rise = std::max(largest_rise, masked[i][4] - unmasked[i][4]);
    }
    if (!(largest_drop > 0.1))
        fail("masking moved no EPNL down by more than 0.1 dB: it is not reaching the tone correction");
    std::cout << "check_band_frames: masking lowered EPNL by up to " << largest_drop << " dB and raised it by up to "
              << largest_rise << " dB.\n";
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.size() == 3 && args[0] == "--history")
        return check_against_history(args[1], args[2]);
    if (args.size() == 3 && args[0] == "--compare")
        return compare_masking(args[1], args[2]);
    ToneMasking masking = ToneMasking::below_noy_floor;
    if (!args.empty() && args[0] == "--no-masking")
    {
        masking = ToneMasking::none;
        args.erase(args.begin());
    }
    if (args.size() != 2 && args.size() != 3)
        fail("usage: check_band_frames [--no-masking] FRAMES.nc EPNL.csv [FULL_HISTORY_LAE.csv]\n"
             "       check_band_frames --history PREFIX FRAMES.nc");
    const std::string frames_path = args[0], epnl_path = args[1];
    int ncid = 0;
    check(nc_open(frames_path.c_str(), NC_NOWRITE, &ncid), "opening " + frames_path);
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
    const auto rows = read_numeric_csv(epnl_path, true);
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
        const EpnlResult r = effective_perceived_noise_level(history, frame_length, true, masking);
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
        if (masking == ToneMasking::below_noy_floor && !std::isnan(row[8]))
        {
            std::size_t tone_band = pnl_band_count;
            for (std::size_t b = 0; b < pnl_band_count; b++)
                if (pnl_band_centers[b] == row[8])
                    tone_band = b;
            if (tone_band == pnl_band_count)
                fail(who + "tone band " + std::to_string(row[8]) + " Hz is not a PNL band");
            const double at_peak = history[r.peak][tone_band];
            if (!(at_peak >= noy_threshold[tone_band]))
                fail(who + "the PNLTM record's tone is in the " + std::to_string(row[8]) + " Hz band at "
                        + std::to_string(at_peak) + " dB, below its noy threshold "
                        + std::to_string(noy_threshold[tone_band]) + " dB");
        }
    }

    if (args.size() == 3)
    {
        // x, y, SEL per observer, from a --full_history run of the same case.
        const auto sel = read_numeric_csv(args[2], true);
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
