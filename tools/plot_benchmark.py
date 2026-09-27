"""Plot niceops_benchmark output in the layout of Attenborough et al. (1995), figs. 11-14.

    plot_benchmark.py DIR

Reads DIR/case<N>_200m.csv and DIR/case<N>_10km.csv for each case present and
writes DIR/case<N>.png: a row per frequency, 200 m on the left and 10 km on
the right, level re the direct field at 1 m from 0 to -80 dB, as the paper
draws them.  The parabolic equation is in ink; the exact solution, where the
case has one, is dashed over it.

Also prints, per case, frequency and span, how far the PE sits from the exact
solution: rms and largest difference over the ranges both are drawn at
(above -80 dB, past 20 m), and the same with the interference nulls -- the
points more than 10 dB below the running median -- left out, since there a
hair's shift in range is tens of dB.
"""

import csv
import math
import os
import statistics
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
EXACT = "#2a78d6"
FLOOR = -80.0
EXACT_NAME = {1: "exact (direct + image ray, eq. 12)", 2: "exact (normal modes, eq. 10)",
              3: "exact (residue series, eq. 11)"}
PROFILE = {1: "uniform", 2: "downward refraction, +0.1 /s", 3: "upward refraction, -0.1 /s",
           4: "duct: +0.1 /s to 100 m, -0.1 /s to 300 m"}


def read(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    data = [[float(x) for x in r] for r in rows[1:] if r]
    cols = {h: [d[i] for d in data] for i, h in enumerate(head)}
    freqs = sorted({h.split("_", 1)[1] for h in head[1:]}, key=float)
    return cols["range_m"], freqs, cols


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=INK_2, labelsize=8, length=3, width=0.6)


def running_median(values, half=10):
    out = []
    for i in range(len(values)):
        window = [v for v in values[max(0, i - half):i + half + 1] if not math.isnan(v)]
        out.append(statistics.median(window) if window else math.nan)
    return out


def agreement(ranges, pe, exact):
    pairs = [(r, a, b) for r, a, b in zip(ranges, pe, exact)
             if r >= 20.0 and not math.isnan(b) and max(a, b) > FLOOR]
    if not pairs:
        return None
    median = running_median([b for _, _, b in pairs])
    all_d = [a - b for _, a, b in pairs]
    off_null = [a - b for (_, a, b), m in zip(pairs, median) if b > m - 10.0]
    rms = lambda d: math.sqrt(sum(x * x for x in d) / len(d))  # noqa: E731
    return (pairs[0][0], pairs[-1][0], rms(all_d), max(abs(x) for x in all_d),
            rms(off_null), max(abs(x) for x in off_null))


def band_label(f):
    f = float(f)
    return f"{f / 1000:g} kHz" if f >= 1000 else f"{f:g} Hz"


def case_figure(case, spans, out):
    freqs = spans[0][1]
    fig, axes = plt.subplots(len(freqs), 2, figsize=(10, 2.7 * len(freqs) + 1.3), facecolor=SURFACE,
                             squeeze=False)
    for col, (name, _, ranges, cols) in enumerate(spans):
        scale = 1.0 if name == "200m" else 0.001
        for row, f in enumerate(freqs):
            ax = axes[row][col]
            style(ax)
            x = [r * scale for r in ranges]
            ax.plot(x, cols[f"pe_{f}"], color=INK, linewidth=1.0, label="parabolic equation")
            exact = cols[f"exact_{f}"]
            if any(not math.isnan(v) for v in exact):
                line, = ax.plot(x, exact, color=EXACT, linewidth=1.3, label=EXACT_NAME[case])
                line.set_dashes((5, 3))
            ax.set_ylim(FLOOR, 0)
            ax.set_xlim(0, ranges[-1] * scale)
            ax.set_yticks(range(-80, 1, 20))
            ax.text(0.97, 0.9, band_label(f), transform=ax.transAxes, ha="right", va="top",
                    fontsize=10, color=INK, fontweight="bold")
            if col == 0:
                ax.set_ylabel("level re 1 m (dB)", color=INK_2, fontsize=8)
            if row == len(freqs) - 1:
                ax.set_xlabel("range (m)" if name == "200m" else "range (km)", color=INK_2, fontsize=8)
    handles, labels = [], []
    for ax in axes.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    height = fig.get_figheight()
    fig.suptitle(f"Case {case}: {PROFILE[case]}", x=0.07, y=0.995, ha="left", fontsize=13, color=INK,
                 fontweight="bold")
    fig.text(0.07, 1.0 - 0.42 / height, "Source 5 m, receiver 1 m, the benchmark's ground and "
             "absorption. Attenborough et al., JASA 97, 173 (1995).", ha="left", va="top",
             fontsize=9, color=INK_2)
    fig.legend(handles, labels, loc="upper left", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.06, 1.0 - 0.58 / height), labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 1.0 - 0.9 / height))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    folder = argv[1]
    print(f"{'case':>4} {'span':>5} {'freq':>6}  {'from':>6} {'to':>7}   {'rms':>5} {'max':>6}"
          f"   {'rms*':>5} {'max*':>6}   (* nulls left out; dB)")
    for case in (1, 2, 3, 4):
        spans = []
        for name in ("200m", "10km"):
            path = os.path.join(folder, f"case{case}_{name}.csv")
            if os.path.exists(path):
                ranges, freqs, cols = read(path)
                spans.append((name, freqs, ranges, cols))
        if len(spans) != 2:
            continue
        out = os.path.join(folder, f"case{case}.png")
        case_figure(case, spans, out)
        for name, freqs, ranges, cols in spans:
            for f in freqs:
                a = agreement(ranges, cols[f"pe_{f}"], cols[f"exact_{f}"])
                if a:
                    print(f"{case:>4} {name:>5} {band_label(f):>6}  {a[0]:6.0f} {a[1]:7.0f}   "
                          f"{a[2]:5.2f} {a[3]:6.2f}   {a[4]:5.2f} {a[5]:6.2f}")
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
