"""Render a demonstration of what an annotated frame from this pipeline looks like.

The real pipeline draws its overlay onto the scene-camera video of an eye-tracking
session. That footage shows the participant's conversation partner, with the
participant's gaze drawn on top, so no real frame can be published. This script
produces the same overlay over a generated stand-in face: the face is detected for
real, and only the gaze path is invented.

    python 00_detection/demo/make_demo_figure.py -o docs/detection-demo.png

The AOI geometry is the same as ``draw.py``: the face ellipse spans the detector's
bounding box at 1.2x its height, and the eye and mouth ellipses span their two
keypoints with the box width and 0.3x its height, rotated to match.

The four areas are shaded as the analysis treats them — mutually exclusive and
exhaustive, following shared_transform_data.R: eyes and mouth take precedence over
the face, "face" means the rest of the face, and "background" is everything else.

By default the committed detection in ``sample_frame_detection.json`` is used, so
the figure reproduces without a detector installed. Pass ``--detect`` to re-run
detection (MTCNN as used by the pipeline, else OpenCV's YuNet).
"""

import argparse
import json
from math import atan2, degrees
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib.patches import Ellipse
import matplotlib.patheffects as pe

HERE = Path(__file__).resolve().parent
# draw.py outlines the face in orange. Here the face area is shaded rather than
# just outlined, and a warm wash over skin reads badly, so it uses a neutral
# slate instead; the eye and mouth colours are unchanged.
FACE_COLOR, EYES_COLOR, MOUTH_COLOR = "#e8edf4", "red", "blue"
FACE_RGB = (0.88, 0.91, 0.96)
EYES_RGB = (0.90, 0.15, 0.15)
MOUTH_RGB = (0.15, 0.30, 0.95)
GAZE_POINT = (1.0, 0.7, 0.25)
GAZE_PATH = (0.0, 1.0, 0.4)
FIXATION = (0.0, 1.0, 1.0)
HIT, MISS = "#31c452", "#f2545b"


def detect(image_path):
    """Locate one face and its five keypoints, preferring the pipeline's detector."""
    try:
        from mtcnn import MTCNN
        import cv2
        img = cv2.cvtColor(cv2.imread(str(image_path)), cv2.COLOR_BGR2RGB)
        res = MTCNN().detect_faces(img)[0]
        return {"box": [float(v) for v in res["box"]],
                "keypoints": {k: [float(v[0]), float(v[1])]
                              for k, v in res["keypoints"].items()},
                "score": float(res["confidence"]), "detector": "mtcnn"}
    except Exception:
        pass
    import cv2
    model = HERE / "yunet.onnx"
    if not model.exists():
        raise SystemExit("no detector available; run without --detect to use the "
                         "committed detection")
    img = cv2.imread(str(image_path))
    h, w = img.shape[:2]
    det = cv2.FaceDetectorYN.create(str(model), "", (w, h), score_threshold=0.6)
    _, faces = det.detect(img)
    f = faces[0]
    names = ["right_eye", "left_eye", "nose", "mouth_right", "mouth_left"]
    return {"box": [float(v) for v in f[:4]],
            "keypoints": {n: [float(f[4 + 2 * i]), float(f[5 + 2 * i])]
                          for i, n in enumerate(names)},
            "score": float(f[-1]), "detector": "yunet"}


def face_ellipse(box):
    x, y, w, h = box
    return (x + w / 2, y + h / 2), w, h * 1.2, 0.0


def span_ellipse(p_right, p_left, box):
    """Eye and mouth AOIs, spanned between two keypoints — mirrors draw_ellipsis."""
    _, _, w, h = box
    diff = np.subtract(p_right, p_left)
    center = np.add(p_left, diff / 2)
    return tuple(center), w, h * 0.3, degrees(atan2(diff[1], diff[0]))


def inside(px, py, ell):
    """Point-in-ellipse, vectorised so it also builds the shading masks."""
    (cx, cy), w, h, angle = ell
    a = np.radians(angle)
    dx, dy = px - cx, py - cy
    xr = dx * np.cos(a) + dy * np.sin(a)
    yr = -dx * np.sin(a) + dy * np.cos(a)
    return (xr / (w / 2)) ** 2 + (yr / (h / 2)) ** 2 <= 1.0


def classify(px, py, face, eyes, mouth):
    """The four exclusive areas, as shared_transform_data.R derives them."""
    in_face, in_eyes, in_mouth = (inside(px, py, e) for e in (face, eyes, mouth))
    face_excl = in_face & ~in_eyes & ~in_mouth
    background = ~face_excl & ~in_eyes & ~in_mouth
    return {"face": face_excl, "eyes": in_eyes,
            "mouth": in_mouth, "background": background}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="docs/detection-demo.png")
    ap.add_argument("-i", "--image", default=str(HERE / "sample_frame.jpg"))
    ap.add_argument("--detect", action="store_true",
                    help="re-run face detection instead of using the committed result")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    det_file = HERE / "sample_frame_detection.json"
    det = detect(args.image) if args.detect else json.loads(det_file.read_text())

    box, kp = det["box"], det["keypoints"]
    face = face_ellipse(box)
    eyes = span_ellipse(kp["right_eye"], kp["left_eye"], box)
    mouth = span_ellipse(kp["mouth_right"], kp["mouth_left"], box)

    img = pyplot.imread(args.image)
    fh, fw = img.shape[:2]

    # An invented scanpath: across one eye to the other, then down to the mouth.
    right_eye, left_eye = np.array(kp["right_eye"]), np.array(kp["left_eye"])
    mouth_c = np.array(mouth[0])
    waypoints = np.array([right_eye, (right_eye + left_eye) / 2, left_eye,
                          (left_eye + mouth_c) / 2, mouth_c])
    u = np.linspace(0, len(waypoints) - 1, 70)
    gx = np.interp(u, np.arange(len(waypoints)), waypoints[:, 0]) + 13 * np.sin(u * 3.1)
    gy = np.interp(u, np.arange(len(waypoints)), waypoints[:, 1]) + 11 * np.cos(u * 3.7)
    fixation = (float(mouth_c[0]), float(mouth_c[1]))

    fig, ax = pyplot.subplots(figsize=(fw / 100, fh / 100))
    ax.imshow(img)

    # Shade the exclusive areas so the partition is visible, not just implied.
    yy, xx = np.mgrid[0:fh, 0:fw]
    masks = classify(xx, yy, face, eyes, mouth)
    overlay = np.zeros((fh, fw, 4))
    for name, rgb, alpha in (("face", FACE_RGB, 0.40),
                             ("eyes", EYES_RGB, 0.26),
                             ("mouth", MOUTH_RGB, 0.26)):
        overlay[masks[name]] = (*rgb, alpha)
    ax.imshow(overlay, zorder=3)

    for ell, col in ((face, FACE_COLOR), (eyes, EYES_COLOR), (mouth, MOUTH_COLOR)):
        patch = Ellipse(ell[0], ell[1], ell[2], angle=ell[3], fill=False,
                        edgecolor=col, lw=2.2, zorder=4)
        if col == FACE_COLOR:      # keep the pale outline legible on any background
            patch.set_path_effects([pe.withStroke(linewidth=3.6, foreground="#33404d")])
        ax.add_patch(patch)

    ax.scatter(gx, gy, color=GAZE_POINT, s=90, alpha=0.22, zorder=5)
    ax.plot(gx, gy, color=GAZE_PATH, lw=1.4, zorder=5)
    ax.scatter(*fixation, s=260, alpha=0.85, facecolors="none",
               edgecolors=[FIXATION], linewidths=1.8, zorder=6)
    ax.text(fixation[0] + 16, fixation[1] + 8, "f014", color=FIXATION,
            fontsize=10, zorder=6)

    # Exactly one area per sample is green, because the areas do not overlap.
    order = ["face", "eyes", "mouth", "background"]
    label_x = fw - 190
    for prefix, pt, y0 in (("gaze", (float(gx[-1]), float(gy[-1])), 40),
                           ("fix", fixation, 210)):
        hits = classify(np.array(pt[0]), np.array(pt[1]), face, eyes, mouth)
        for n, name in enumerate(order):
            ax.text(label_x, y0 + n * 34, f"{prefix}_{name}",
                    color=HIT if bool(hits[name]) else MISS,
                    fontsize=13, fontweight="bold", zorder=7)

    ax.set_xlim(0, fw); ax.set_ylim(fh, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout(pad=0)
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight", pad_inches=0)
    print(f"wrote {args.out}  (detector: {det.get('detector', 'committed')}, "
          f"score {det.get('score', float('nan')):.3f})")


if __name__ == "__main__":
    main()
