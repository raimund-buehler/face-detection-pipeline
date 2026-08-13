"""Render the pipeline overview used at the top of the README.

    python docs/make_pipeline_diagram.py -o docs/pipeline.svg

Two lanes: the computer-vision stage that turns video into AOI hits, and the
statistical stage that turns those into models and figures. Written as an SVG so
it stays crisp at whatever width the README is viewed at.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "path"
from matplotlib import pyplot
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

BG = "#fbfbfd"
INK = "#1f2933"
MUTED = "#7b8794"
DETECT = "#e0942f"
ANALYSE = "#3d7fbf"

BW, BH, GAP = 3.15, 1.30, 0.90
X0 = 0.55

DETECTION = ["scene video\n+ gaze samples",
             "MTCNN face detection\nper frame",
             "face / eye / mouth AOIs\ngaze & fixation hit test"]
ANALYSIS = ["per-session\ndwell time & rates",
            "merge questionnaire\n& trial data",
            "mixed-effects models\nfigures & tables"]


def box(ax, x, y, label, color, tint=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), BW, BH, boxstyle="round,pad=0.05,rounding_size=0.14",
        linewidth=1.8, edgecolor=color,
        facecolor="#f2f4f7" if tint else "#ffffff", zorder=3))
    ax.text(x + BW / 2, y + BH / 2, label, ha="center", va="center",
            fontsize=11, color=INK, zorder=4, linespacing=1.5)


def straight(ax, x1, y1, x2, y2, color):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=16, linewidth=1.7,
                                 color=color, zorder=2))


def lane(ax, items, y, color):
    for i, label in enumerate(items):
        x = X0 + i * (BW + GAP)
        box(ax, x, y, label, color, tint=(i == 0 and y > 2))
        if i:
            straight(ax, x - GAP + 0.06, y + BH / 2, x - 0.06, y + BH / 2, color)
    return X0 + (len(items) - 1) * (BW + GAP) + BW


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="docs/pipeline.svg")
    args = ap.parse_args()

    fig, ax = pyplot.subplots(figsize=(13.0, 5.0))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    y_top, y_bot = 2.95, 0.72
    end_x = lane(ax, DETECTION, y_top, DETECT)
    lane(ax, ANALYSIS, y_bot, ANALYSE)

    # Elbow from the end of the detection lane round to the start of the analysis
    # lane: out to the right, down into the gap, back along, then into the box.
    mid_y = y_top - 0.24          # just under the detection row
    first_cx = X0 + BW / 2
    ax.plot([end_x + 0.06, end_x + 0.52], [y_top + BH / 2, y_top + BH / 2],
            color=MUTED, lw=1.7, zorder=2, solid_capstyle="round")
    ax.plot([end_x + 0.52, end_x + 0.52], [y_top + BH / 2, mid_y],
            color=MUTED, lw=1.7, zorder=2, solid_capstyle="round")
    ax.plot([end_x + 0.52, first_cx], [mid_y, mid_y],
            color=MUTED, lw=1.7, zorder=2, solid_capstyle="round")
    straight(ax, first_cx, mid_y, first_cx, y_bot + BH + 0.06, MUTED)

    for y, text, color in ((y_top, "DETECTION      Python", DETECT),
                           (y_bot, "ANALYSIS      R", ANALYSE)):
        ax.text(X0, y + BH + 0.14, text, fontsize=11, color=color,
                fontweight="bold", va="bottom")

    ax.set_xlim(0, 12.85)
    ax.set_ylim(0.25, 4.85)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout(pad=0.25)
    fig.savefig(args.out, facecolor=fig.get_facecolor())
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
