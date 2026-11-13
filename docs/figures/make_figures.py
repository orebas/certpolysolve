#!/usr/bin/env python3
"""The documentation's figures, drawn from the solver's own output (docs/WRITING.md: a picture is generated from real
output by a script kept beside it). Needs matplotlib.

    make_figures.py [PNG_DIR]        writes the SVGs here; with PNG_DIR, PNG copies there for checking

Data, in data/: the command-line tool's output for a degree-24 polynomial (degree24-*.txt), and for the pair of
roots 1 and 1 + 2^-60 under two goals (pair-clusters.txt: --clusters 1e-6; pair-accuracy.txt: the default);
time-by-degree.csv: one line per solve of the test set (degree, seconds on one core, query); compare-small-dense.csv:
the median time of one solve by degree and solver on dense polynomials with random coefficients. The commands that made
the data are in data/README.md."""
import csv
import math
import re
import statistics
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, GRID, HUE, HUE_LIGHT = "#333333", "#6b6b6b", "#dddddd", "#2e6fbf", "#9dbde6"
plt.rcParams.update({"svg.fonttype": "none", "font.family": "sans-serif", "font.size": 10, "text.color": INK,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.hashsalt": "certpolysolve-docs"})

END = re.compile(r"(-?\d+) \* 2\^\((-?\d+)\)")


def boxes(path):
    """Each box of a CLI output: (real lo, real hi, imag lo, imag hi) as Fractions, and its count."""
    out = []
    for line in (HERE / "data" / path).read_text().splitlines():
        if not line.startswith("real="):
            continue
        ends = [Fraction(int(m)) * Fraction(2) ** int(e) for m, e in END.findall(line)]
        count = int(re.search(r"count=(\d+)", line).group(1))
        out.append((ends[0], ends[1], ends[2], ends[3], count))
    return out


def save(fig, name, png_dir):
    # No date in the file and ids from a fixed salt (svg.hashsalt above): the same data gives the same bytes.
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    if png_dir:
        fig.savefig(Path(png_dir) / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def roots_figure(png_dir):
    """The 24 roots of the degree-24 polynomial, each the centre of a proved box too small to draw (the widths are
    returned for the caption). Asking for larger boxes (--clusters 0.05) gave boxes of 10^-13 and less all the
    same: the solver proves a tight box whenever that is cheap, so this picture is dots."""
    rows = boxes("degree24-roots.txt")
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.set_aspect("equal")
    ax.axhline(0, color=GRID, lw=1, zorder=0)
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, color=GRID, lw=1, ls="--", zorder=0))
    ax.scatter([float((a + b) / 2) for a, b, _, _, _ in rows], [float((c + d) / 2) for _, _, c, d, _ in rows],
               s=36, color=HUE, zorder=3)
    ax.set_xlabel("real part")
    ax.set_ylabel("imaginary part")
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-2.1, 2.1)
    save(fig, "roots-degree-24", png_dir)
    widths = [max(float(b - a), float(d - c)) for a, b, c, d, _ in rows]
    return len(rows), min(widths), max(widths)


def pair_figure(png_dir):
    """The roots 1 and 1 + 2^-60 under two goals: boxes of radius at most 1e-6, and each root alone."""
    clusters = [r for r in boxes("pair-clusters.txt") if r[0] > 0]
    accuracy = [r for r in boxes("pair-accuracy.txt") if r[0] > 0]
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.4, 2.9))
    fig.subplots_adjust(wspace=0.35)
    for ax in (left, right):
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.axhline(0, color=GRID, lw=1, zorder=0)
        ax.set_ylim(-1, 1)
    # left: the cluster box, at the scale of its radius
    (a, b, c, d, cluster_count), = clusters
    span = Fraction(5, 10**6)
    left.add_patch(plt.Rectangle((float(a - 1), -0.35), float(b - a), 0.7, fill=True, facecolor=HUE_LIGHT,
                                 edgecolor=HUE, lw=1.5))
    left.text(float((a + b) / 2 - 1), 0.55, f"one box, holds {cluster_count} roots", ha="center", color=INK, fontsize=9)
    left.set_xlim(float(-span), float(span))
    left.set_xticks([float(-span), 0, float(span)], ["−5·10⁻⁶", "0", "5·10⁻⁶"])
    left.set_xlabel("distance from 1")
    left.set_title("asked: boxes of radius at most 10⁻⁶", fontsize=10, color=INK)
    # right: the two boxes, at the scale of their distance
    offsets = [((a + b) / 2 - 1) for a, b, _, _, _ in accuracy]
    scale = Fraction(10) ** 18
    xs = [float(o * scale) for o in offsets]
    right.scatter(xs, [0, 0], s=40, color=HUE, zorder=3)
    for x, (a, b, _, _, count), y in zip(xs, accuracy, (0.45, -0.7)):
        right.text(x, y, f"{count} root, box width {float(b - a):.0e}", ha="center", color=INK, fontsize=8)
    right.set_xlim(-0.5, 1.5)
    right.set_xlabel("distance from 1, in units of 10⁻¹⁸")
    right.set_title("asked: each root alone (the default)", fontsize=10, color=INK)
    save(fig, "close-pair", png_dir)
    return cluster_count, [float(o) for o in offsets]


def time_figure(png_dir):
    """Every solve of the test set: seconds against degree, both logarithmic, with the median per band of degree."""
    rows = []
    with open(HERE / "data" / "time-by-degree.csv") as f:
        for row in csv.DictReader(f):
            rows.append((int(row["degree"]), float(row["seconds"])))
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.scatter([d for d, _ in rows], [t for _, t in rows], s=6, color=HUE, alpha=0.25, linewidths=0, zorder=2)
    bands = {}
    for d, t in rows:
        bands.setdefault(int(math.log2(d)), []).append(t)
    xs = [2 ** k * math.sqrt(2) for k in sorted(bands) if len(bands[k]) >= 5]
    ys = [statistics.median(bands[k]) for k in sorted(bands) if len(bands[k]) >= 5]
    ax.plot(xs, ys, color=INK, lw=2, marker="o", ms=5, zorder=3, label="median of each band of degree")
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_xlabel("degree")
    ax.set_ylabel("seconds, one core")
    ax.set_xticks([2 ** k for k in range(1, 19, 2)], [f"{2 ** k:,}" for k in range(1, 19, 2)])
    ax.legend(frameon=False, loc="upper left")
    save(fig, "time-by-degree", png_dir)
    return len(rows)


def comparison_figure(png_dir):
    """The comparison on dense polynomials with random coefficients in [-1, 1]: the median time of one solve against
    the degree, for CertPolySolve, MPSolve and numpy.roots, from data/compare-small-dense.csv. numpy's label gives the
    accuracy its roots reached, since it promises none."""
    rows = {}
    with open(HERE / "data" / "compare-small-dense.csv") as f:
        for row in csv.DictReader(f):
            rows.setdefault(row["arm"], []).append((int(row["degree"]), float(row["median_solve_seconds"]), float(row["median_bits"])))
    series = [("auto", "CertPolySolve, every root proved", HUE), ("mpsolve-secular", "MPSolve", "#eb6834"),
              ("numpy", "numpy.roots", "#1baf7a")]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    for arm, label, colour in series:
        points = sorted(rows[arm])
        if arm == "numpy":
            low, high = min(b for *_, b in points), max(b for *_, b in points)
            label = f"{label}, roots correct to {low:.0f} to {high:.0f} bits"
        ax.plot([d for d, *_ in points], [s for _, s, _ in points], color=colour, lw=2, marker="o", ms=4, label=label, zorder=3)
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_xlabel("degree")
    ax.set_ylabel("seconds for one solve, one core")
    ax.set_xticks([2 ** k for k in range(1, 10)], [f"{2 ** k}" for k in range(1, 10)])
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    save(fig, "compare-small-dense", png_dir)
    return sum(len(v) for v in rows.values())


if __name__ == "__main__":
    png_dir = sys.argv[1] if len(sys.argv) > 1 else None
    n, narrowest, widest = roots_figure(png_dir)
    print(f"roots-degree-24.svg: {n} roots, their boxes {narrowest:.1e} to {widest:.1e} wide")
    count, offsets = pair_figure(png_dir)
    print(f"close-pair.svg: the cluster box holds {count}; the two roots are {offsets[1] - offsets[0]:.3e} apart")
    print(f"time-by-degree.svg: {time_figure(png_dir)} solves")
    print(f"compare-small-dense.svg: {comparison_figure(png_dir)} rows")
