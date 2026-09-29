"""Plot niceops_pe_compare output: the refraction models against the parabolic equation.

    plot_pe_compare.py OUTDIR CASE [CASE ...]

Each CASE is CSV:TITLE, a niceops_pe_compare output and a short title for it.
Writes one figure per case -- a panel per band, level relative to free field
against range for each model -- and a summary of where each model's shadow
begins (the first range from which the level stays 10 dB below free field),
against frequency.

The parabolic equation is the reference, so it is drawn in ink; the three ray
models take the first three categorical slots, with distinct dashes so their
identity never rests on colour alone.
"""

import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
MODELS = [
    # (column prefix, label, colour, dashes, width)
    ("pe", "parabolic equation", INK, (None, None), 2.4),
    ("stratified", "stratified rays", "#2a78d6", (None, None), 1.8),
    ("linear", "linear rays", "#eb6834", (6, 3), 1.8),
    ("straight", "straight rays", "#1baf7a", (1.5, 2.5), 1.8),
]
THRESHOLD = -10.0
FLOOR = -45.0  # the plots stop here; below it only the turbulence-free PE goes


def read(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    data = [[float(x) for x in r] for r in rows[1:] if r]
    unit = head[0].split("_", 1)[1]
    to_km = 0.0003048 if unit == "ft" else 0.001
    ranges = [d[0] * to_km for d in data]
    bands = sorted({h.split("_", 1)[1] for h in head[1:]}, key=float)
    series = {}
    for b in bands:
        for prefix, *_ in MODELS:
            i = head.index(f"{prefix}_{b}")
            series[(prefix, b)] = [d[i] for d in data]
    return ranges, bands, series


def onset(ranges, levels):
    for j, v in enumerate(levels):
        if v < THRESHOLD and all(w < THRESHOLD for w in levels[j:]):
            return ranges[j]
    return None


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.6)


def band_label(b):
    f = float(b)
    return f"{f / 1000:g} kHz" if f >= 1000 else f"{f:g} Hz"


def case_figure(path, title, out):
    ranges, bands, series = read(path)
    ncols = 2
    nrows = math.ceil(len(bands) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 3.6 * nrows + 1.2), sharex=True,
                             sharey=True, facecolor=SURFACE, squeeze=False)
    lowest = min(v for values in series.values() for v in values)
    # Down to the data, in 5 dB steps, but no further than FLOOR.
    bottom = max(FLOOR, 5.0 * math.floor(lowest / 5.0) - 5.0)
    for k, b in enumerate(bands):
        ax = axes[k // ncols][k % ncols]
        style(ax)
        ax.axhline(THRESHOLD, color=INK_2, linewidth=0.8, linestyle=(0, (2, 3)))
        for prefix, label, colour, dashes, width in MODELS:
            line, = ax.plot(ranges, [max(v, bottom - 5) for v in series[(prefix, b)]], color=colour,
                            linewidth=width, label=label, solid_capstyle="round")
            if dashes[0] is not None:
                line.set_dashes(dashes)
        ax.set_title(band_label(b), loc="left", fontsize=11, color=INK, fontweight="bold")
        ax.set_ylim(bottom, 12)
        ax.set_xlim(0, ranges[-1])
        if k % ncols == 0:
            ax.set_ylabel("level re free field (dB)", color=INK_2, fontsize=9)
        if k // ncols == nrows - 1:
            ax.set_xlabel("range (km)", color=INK_2, fontsize=9)
        if k == 0:
            ax.text(ranges[-1] * 0.99, THRESHOLD + 1.2, "10 dB below free field", ha="right",
                    va="bottom", fontsize=8, color=INK_2)
    for k in range(len(bands), nrows * ncols):
        axes[k // ncols][k % ncols].set_visible(False)
    top = 1.0 - 0.95 / fig.get_figheight()  # the header's share, in figure fraction
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.suptitle(title, x=0.06, y=0.995, ha="left", fontsize=13, color=INK, fontweight="bold")
    note = ("Band-averaged; absorption and turbulence left out of all four."
            + (f" Clipped at {FLOOR:g} dB: the PE, with nothing to scatter energy in, falls to "
               f"{lowest:.0f} dB in the deep shadow." if lowest < FLOOR else ""))
    fig.text(0.06, 1.0 - 0.45 / fig.get_figheight(), note, ha="left", va="top", fontsize=9,
             color=INK_2)
    fig.legend(handles, labels, loc="upper left", ncol=4, frameon=False, fontsize=10,
               bbox_to_anchor=(0.05, 1.0 - 0.62 / fig.get_figheight()), labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, top))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def onset_figure(cases, out):
    # Only the cases with a shadow belong here: the PE's level falls for good
    # in at least half the bands.  Straight rays cast no shadow -- a crossing
    # of theirs is ground interference -- so they are left out.
    models = [m for m in MODELS if m[0] != "straight"]
    shadowed = []
    for path, title in cases:
        ranges, bands, series = read(path)
        found = {prefix: [onset(ranges, series[(prefix, b)]) for b in bands]
                 for prefix, *_ in models}
        if 2 * sum(v is not None for v in found["pe"]) >= len(bands):
            shadowed.append((title, [float(b) for b in bands], found))
    if not shadowed:
        return False
    fig, axes = plt.subplots(1, len(shadowed), figsize=(5.4 * len(shadowed), 4.6),
                             facecolor=SURFACE, squeeze=False)
    for ax, (title, freqs, found) in zip(axes[0], shadowed):
        style(ax)
        for prefix, label, colour, dashes, width in models:
            points = [(f, v) for f, v in zip(freqs, found[prefix]) if v is not None]
            if not points:
                continue
            line, = ax.plot([p[0] for p in points], [p[1] for p in points], color=colour,
                            linewidth=width, marker="o", markersize=8, markeredgecolor=SURFACE,
                            markeredgewidth=2, label=label)
            if dashes[0] is not None:
                line.set_dashes(dashes)
        ax.set_xscale("log")
        ax.set_xticks(freqs)
        ax.set_xticklabels([band_label(str(f)) for f in freqs])
        ax.minorticks_off()
        ax.set_ylim(bottom=0)
        ax.set_title(title, loc="left", fontsize=11, color=INK, fontweight="bold")
        ax.set_xlabel("band", color=INK_2, fontsize=9)
        ax.set_ylabel("shadow begins at (km)", color=INK_2, fontsize=9)
    height = fig.get_figheight()
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.suptitle("Where the shadow begins upwind", x=0.06, y=0.995, ha="left", fontsize=13,
                 color=INK, fontweight="bold")
    fig.text(0.06, 1.0 - 0.45 / height, "First range from which the level stays 10 dB below "
             "free field. Straight rays cast no shadow and are left out.", ha="left", va="top",
             fontsize=9, color=INK_2)
    fig.legend(handles, labels, loc="upper left", ncol=4, frameon=False, fontsize=10,
               bbox_to_anchor=(0.05, 1.0 - 0.62 / height), labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 1.0 - 0.95 / height))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    outdir = argv[1]
    os.makedirs(outdir, exist_ok=True)
    cases = []
    for spec in argv[2:]:
        path, title = spec.split(":", 1)
        cases.append((path, title))
        name = os.path.splitext(os.path.basename(path))[0] + ".png"
        case_figure(path, title, os.path.join(outdir, name))
        print(os.path.join(outdir, name))
    if onset_figure(cases, os.path.join(outdir, "shadow_onset.png")):
        print(os.path.join(outdir, "shadow_onset.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
