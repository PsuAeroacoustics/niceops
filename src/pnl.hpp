#pragma once

// Perceived noise level metrics: PNL, the tone-corrected PNLT, PNLTM and EPNL,
// as 14 CFR 36 Appendix A, section A36.4, defines them (identical to ICAO
// Annex 16 Vol. I, Appendix 2).
//
// Ported from panam's flight_acoustics.py (noys, perceived_noise_level,
// tone_correction, tone_corrected_perceived_noise_level,
// effective_perceived_noise_level; panam eeb040b, with F005's treatment of
// -infinity and NaN bands), term for term, and pinned against values generated
// by running those functions (tests/test_pnl.cpp, 1e-6 dB).  panam's own
// reading of the regulation is kept where it chose:
//   * 500 Hz belongs to the 500 <= f <= 5000 Hz range of Table A36-2;
//   * bands 1 to 3 (50 to 80 Hz) never carry a tone correction -- the slopes
//     start at band 4, and band 3 is the background's anchor;
//   * the duration's limits are the records closest to PNLTM - 10 and, with
//     several peaks, those giving the longest duration (A36.4.5.5; panam
//     F026): from the outermost records at or above it, each limit moves one
//     record outward when that one is strictly closer; dips inside count;
//   * the duration correction takes 10 lg(dt / T), T = 10 s, as the -13 that
//     A36.4.5.4 states for dt = 0.5 s, and exactly for any other dt;
//   * the A36.4.4.2 band-sharing check is applied (optional, as in panam),
//     its adjustment Delta_B added to EPNL as a term of its own (panam
//     F014), PNLTM, the duration and D being the unadjusted history's: a
//     departure from ICAO Annex 16 Vol. I, Appendix 2 and AC 36-4D, which
//     add it to PNLTM (THEORY 10.5);
//   * a -infinity band (zero energy), for which the procedure is undefined,
//     enters steps 1-7 at its noy threshold SPL(d) and carries no tone (it
//     has no noys either way); taken literally it made its neighbor a full
//     tone, or, as band 3, the background's anchor, cost every tone;
//   * a NaN band is NaN noy, so PNL and the tone correction are NaN, and an
//     EPNL history holding one is refused; +infinity is refused.
// Nothing else of A36.4 is: no correction of the measured spectra to reference
// conditions (A36.9) and no treatment of pseudo-tones from ground reflection.
//
// One addition, off by default here and on in --epnl_file: ToneMasking::
// below_noy_floor, the prediction's counterpart of the certification practice
// of excluding tones in masked bands (FAA AC 36-4, Appendix 2, item 11).  A
// measured band never falls below the analyzer's noise floor; a predicted one
// falls as far as absorption takes it -- hundreds of dB below 0 dB at 10 kHz
// a kilometer out, and -1000 dB where a database band is near silent.  The
// literal procedure then reads the plunge, or an inaudible band rising out of
// it, as a tone: on the demo track over a 100 m grid, 615 of 961 observers
// took their PNLTM tone correction from a band below its noy threshold, up to
// the full 20/3 dB, and EPNL read a median 1.0 dB (up to 5.4 dB) high.
// Masked, a band below SPL(d) of Table A36-3 (the level at which it stops
// contributing noys) enters steps 1-7 at SPL(d), and carries no tone
// correction itself (panam's masked=True, since F005).  Both halves are needed: excluding
// the masked bands' tones alone leaves an audible band before a plunge with
// its background crashed (a 20/3 dB tone), and the reconstruction alone can
// make masked bands rising along SPL(d) read as one.  PNL is unchanged: such a
// band is 0 noy either way.  A spectrum with every band at or above SPL(d) is
// untouched.  EPNL can still rise a little (0.2 dB at most on that grid): the
// band-sharing check and the duration may pick other records.
//
// Band convention: the 24 one-third-octave bands 50 Hz to 10 kHz, ascending,
// at their nominal centers (pnl_band_centers).

#include <array>
#include <cstddef>
#include <vector>

namespace niceops {

inline constexpr std::size_t pnl_band_count = 24;

// The record length A36.4 prescribes: one-third-octave spectra averaged
// linearly over 0.5 s.
inline constexpr double pnl_record_seconds = 0.5;

// Nominal one-third-octave centers, Hz, 50 Hz (band 1) to 10 kHz (band 24).
inline constexpr std::array<double, pnl_band_count> pnl_band_centers = {50.0, 63.0, 80.0, 100.0,
        125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0,
        2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0};

// One record's band levels, dB re 20 uPa, in pnl_band_centers' order.
using PnlSpectrum = std::array<double, pnl_band_count>;

// Perceived noisiness of each band in noys, A36.4.7.3 with Table A36-3.  Zero
// below SPL(d), and for a band level of -infinity or NaN.
PnlSpectrum noys(const PnlSpectrum& band_levels);

// PNL, A36.4.2: 40 + (10 / lg 2) lg N, N = 0.85 n_max + 0.15 sum n.  Minus
// infinity when every band is below its noy floor.
double perceived_noise_level(const PnlSpectrum& band_levels);

// Table A36-3's SPL(d), dB: below it a band contributes no noys.
extern const PnlSpectrum noy_threshold;

// How the tone correction treats a band below its noy threshold (see above):
// `none` is the regulation's procedure as written, and panam's.
enum class ToneMasking
{
    none,
    below_noy_floor,
};

// The spectral-irregularity (tone) correction C, A36.4.3.1 steps 1-10: the
// largest of the per-band factors and the band (0-23) it came from -- band 0
// when there is no tone, as numpy's argmax gives.  A spectrum with a NaN band
// has a NaN factor and band pnl_band_count (panam's -1); +infinity throws
// std::invalid_argument.
struct ToneCorrection
{
    double factor = 0.0;
    std::size_t band = 0;
};
ToneCorrection tone_correction(const PnlSpectrum& band_levels, ToneMasking masking = ToneMasking::none);

// EPNL of a noise event from its band-level history at equal steps `dt`
// (the regulation's 0.5 s), A36.4.
struct EpnlResult
{
    double epnl = 0.0;                // EPNdB, PNLTM + D + Delta_B
    double pnltm = 0.0;               // TPNdB, after the band-sharing check (unadjusted + Delta_B)
    double pnltm_unadjusted = 0.0;    // TPNdB, the largest PNLT
    double delta_b = 0.0;             // the band-sharing adjustment, dB (0 when not applied)
    double duration_correction = 0.0; // D, dB, from the unadjusted PNLTM; EPNL = pnltm + D
    std::size_t peak = 0;             // the record PNLTM is taken at
    std::size_t first = 0;            // the duration's first and last records
    std::size_t last = 0;
    // Whether the duration touches either end of the record, so the event may
    // not have been captured whole.
    bool clipped = false;
    // Whether the band-sharing check applied Delta_B (A36.4.4.2).
    bool bandshare_adjusted = false;
    // Per-record histories.
    std::vector<double> pnlt;
    std::vector<double> pnl;
    std::vector<double> tone;
    std::vector<std::size_t> tone_band;
};

// Throws std::invalid_argument when the history is empty, holds a NaN band, or
// holds no finite PNLT (every record below the noy floor), as panam raises
// ValueError: an EPNL of nothing is not a number.  The caller decides what to report then.
EpnlResult effective_perceived_noise_level(const std::vector<PnlSpectrum>& history,
        double dt = pnl_record_seconds, bool bandshare_adjustment = true,
        ToneMasking masking = ToneMasking::none);

// Where each of the 24 PNL bands sits among a database's band centers: the
// band whose center is within 5 percent of the nominal one (nominal and exact
// base-ten or base-two centers differ by at most 1.6 percent; neighboring
// third-octave bands are 26 percent apart).  `missing` receives the nominal
// centers no band matches, and an entry of the result is then meaningless.
std::array<std::size_t, pnl_band_count> pnl_band_indices(const std::vector<double>& band_centers,
        std::vector<double>& missing);

} // namespace niceops
