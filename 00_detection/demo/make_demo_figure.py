"""Render a demonstration of what an annotated frame from this pipeline looks like.

The real pipeline draws its overlay onto the scene-camera video of an eye-tracking
session. That footage shows the participant's conversation partner and carries the
participant's gaze on top of it, so no real frame can be published. This script
reproduces the same overlay — identical colours, geometry and hit-test semantics as
``draw.py`` and ``make_frame.py`` — over a schematic stand-in face with an invented
gaze path, so the output is illustrative and contains no data of any kind.

    python 00_detection/demo/make_demo_figure.py -o docs/detection-demo.png

Only numpy, matplotlib and Pillow are needed; the detector itself is not run.
"""

import argparse
from math import atan2, degrees

import numpy as np
from matplotlib import pyplot
from matplotlib.patches import Ellipse

# Frame geometry and palette, matching make_frame.py / draw.py.
FRAME_W, FRAME_H = 1280, 720
FACE_COLOR, EYES_COLOR, MOUTH_COLOR = "orange", "red", "blue"
GAZE_POINT = (1.0, 0.7, 0.25)
GAZE_PATH = (0.0, 1.0, 0.4)
FIXATION = (0.0, 1.0, 1.0)
LABEL_X = 980


def ellipse_from_box(x, y, width, height):
    """Face AOI: centred on the detector's bounding box, 1.2x its height.

    Mirrors draw_ellipsis(..., face=True).
    """
    return (x + width / 2, y + height / 2), width, height * 1.2, 0.0


def ellipse_from_points(p_right, p_left, width, height):
    """Eye and mouth AOIs: spanned between two keypoints and rotated to match them.

    Mirrors draw_ellipsis(..., face=False).
    """
    diff = np.subtract(p_right, p_left)
    center = np.add(p_left, diff / 2)
    angle = degrees(atan2(diff[1], diff[0]))
    return tuple(center), width / 2.5, height / 5, angle


def inside(point, center, w, h, angle):
    """Hit test for a gaze or fixation sample against one AOI ellipse."""
    a = np.radians(angle)
    dx, dy = point[0] - center[0], point[1] - center[1]
    xr = dx * np.cos(a) + dy * np.sin(a)
    yr = -dx * np.sin(a) + dy * np.cos(a)
    return (xr / (w / 2)) ** 2 + (yr / (h / 2)) ** 2 <= 1.0


def draw_stand_in_face(ax, box):
    """A deliberately schematic head, so nothing here resembles a real person."""
    x, y, w, h = box
    ax.add_patch(Ellipse((x + w / 2, y + h / 2), w * 1.02, h * 1.25,
                         facecolor="#d8d2cc", edgecolor="none", zorder=1))
    ax.add_patch(Ellipse((x + w / 2, y + h * 1.46), w * 1.25, h * 0.5,
                         facecolor="#2f3640", edgecolor="none", zorder=0))
    for cx in (x + w * 0.31, x + w * 0.69):
        ax.add_patch(Ellipse((cx, y + h * 0.40), w * 0.13, h * 0.07,
                             facecolor="#3a3a3a", edgecolor="none", zorder=2))
    ax.add_patch(Ellipse((x + w / 2, y + h * 0.74), w * 0.26, h * 0.05,
                         facecolor="#8c5b5b", edgecolor="none", zorder=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="docs/detection-demo.png")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    # A plausible detection: bounding box plus the five MTCNN keypoints.
    box = (470.0, 150.0, 330.0, 380.0)
    x, y, w, h = box
    right_eye = (x + w * 0.31, y + h * 0.40)
    left_eye = (x + w * 0.69, y + h * 0.40)
    mouth_right = (x + w * 0.37, y + h * 0.74)
    mouth_left = (x + w * 0.63, y + h * 0.74)

    face_c, face_w, face_h, face_a = ellipse_from_box(x, y, w, h)
    eyes_c, eyes_w, eyes_h, eyes_a = ellipse_from_points(right_eye, left_eye, w, h)
    mouth_c, mouth_w, mouth_h, mouth_a = ellipse_from_points(mouth_right, mouth_left, w, h)

    # An invented gaze excursion running from the eye region down to the mouth, so
    # the hit-test labels below show both a miss and a hit rather than all misses.
    t = np.linspace(0, 1, 42)
    gaze_x = eyes_c[0] - 55 + (mouth_c[0] - eyes_c[0] + 55) * t + 22 * np.sin(t * 7)
    gaze_y = eyes_c[1] + (mouth_c[1] - eyes_c[1]) * t ** 1.4 + 12 * np.cos(t * 9)
    fixation = (float(mouth_c[0]), float(mouth_c[1]))

    fig, ax = pyplot.subplots(figsize=(FRAME_W / 100, FRAME_H / 100))
    ax.set_facecolor("#11151a")
    fig.patch.set_facecolor("#11151a")
    draw_stand_in_face(ax, box)

    for (c, ew, eh, a), col in ((( face_c,  face_w,  face_h,  face_a), FACE_COLOR),
                                (( eyes_c,  eyes_w,  eyes_h,  eyes_a), EYES_COLOR),
                                ((mouth_c, mouth_w, mouth_h, mouth_a), MOUTH_COLOR)):
        ax.add_patch(Ellipse(c, ew, eh, angle=a, fill=False, edgecolor=col, lw=2, zorder=4))

    ax.scatter(gaze_x, gaze_y, color=GAZE_POINT, s=100, alpha=0.2, zorder=5)
    ax.plot(gaze_x, gaze_y, color=GAZE_PATH, lw=1, zorder=5)
    ax.scatter(*fixation, color=FIXATION, s=200, alpha=0.5, facecolors="none", zorder=6)
    ax.text(fixation[0] + 12, fixation[1] + 6, "f014", color=FIXATION, fontsize=9, zorder=6)

    # The six status labels, green where the sample falls inside the AOI.
    aois = (("face", face_c, face_w, face_h, face_a),
            ("eyes", eyes_c, eyes_w, eyes_h, eyes_a),
            ("mouth", mouth_c, mouth_w, mouth_h, mouth_a))
    last_gaze = (float(gaze_x[-1]), float(gaze_y[-1]))
    for n, (name, c, ew, eh, a) in enumerate(aois):
        for prefix, pt, y0 in (("gaze", last_gaze, 40), ("fix", fixation, 160)):
            hit = inside(pt, c, ew, eh, a)
            ax.text(LABEL_X, y0 + n * 40, f"{prefix}_{name}",
                    color="green" if hit else "red", fontsize=13, zorder=7)

    ax.text(30, 40, "synthetic demonstration — no participant data",
            color="#7d8896", fontsize=11, style="italic", zorder=7)

    ax.set_xlim(0, FRAME_W)
    ax.set_ylim(FRAME_H, 0)          # image coordinates: origin top-left
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout(pad=0)
    fig.savefig(args.out, dpi=args.dpi, facecolor=fig.get_facecolor())
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
