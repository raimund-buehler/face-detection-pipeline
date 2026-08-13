"""Render a demonstration of what an annotated frame from this pipeline looks like.

The real pipeline draws its overlay onto the scene-camera video of an eye-tracking
session. That footage shows the participant's conversation partner, with the
participant's gaze drawn on top, so no real frame can be published. This script
produces the same overlay over a public-domain stock portrait instead: the face is
detected for real, and only the gaze path is invented.

    python 00_detection/demo/make_demo_figure.py -o docs/detection-demo.png

The AOI geometry is the same as ``draw.py``: the face ellipse spans the detector's
bounding box at 1.2x its height, and the eye and mouth ellipses span their two
keypoints with the box width and 0.3x its height, rotated to match.

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

HERE = Path(__file__).resolve().parent
FACE_COLOR, EYES_COLOR, MOUTH_COLOR = "orange", "red", "blue"
GAZE_POINT = (1.0, 0.7, 0.25)
GAZE_PATH = (0.0, 1.0, 0.4)
FIXATION = (0.0, 1.0, 1.0)


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


def inside(point, center, w, h, angle):
    a = np.radians(angle)
    dx, dy = point[0] - center[0], point[1] - center[1]
    xr = dx * np.cos(a) + dy * np.sin(a)
    yr = -dx * np.sin(a) + dy * np.cos(a)
    return (xr / (w / 2)) ** 2 + (yr / (h / 2)) ** 2 <= 1.0


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

    box = det["box"]
    kp = det["keypoints"]
    face_c, face_w, face_h, face_a = face_ellipse(box)
    eyes_c, eyes_w, eyes_h, eyes_a = span_ellipse(kp["right_eye"], kp["left_eye"], box)
    mouth_c, mouth_w, mouth_h, mouth_a = span_ellipse(kp["mouth_right"], kp["mouth_left"], box)

    img = pyplot.imread(args.image)
    fh, fw = img.shape[:2]

    # An invented scanpath: across one eye to the other, then down to the mouth.
    # It ends inside the mouth AOI so the labels show both a hit and a miss.
    right_eye, left_eye = np.array(kp["right_eye"]), np.array(kp["left_eye"])
    waypoints = np.array([right_eye, (right_eye + left_eye) / 2, left_eye,
                          (left_eye + np.array(mouth_c)) / 2, np.array(mouth_c)])
    u = np.linspace(0, len(waypoints) - 1, 70)
    gx = np.interp(u, np.arange(len(waypoints)), waypoints[:, 0])
    gy = np.interp(u, np.arange(len(waypoints)), waypoints[:, 1])
    gx = gx + 13 * np.sin(u * 3.1)          # sampling jitter
    gy = gy + 11 * np.cos(u * 3.7)
    fixation = (float(mouth_c[0]), float(mouth_c[1]))

    fig, ax = pyplot.subplots(figsize=(fw / 100, fh / 100))
    ax.imshow(img)

    for (c, w, h, a), col in (((face_c, face_w, face_h, face_a), FACE_COLOR),
                              ((eyes_c, eyes_w, eyes_h, eyes_a), EYES_COLOR),
                              ((mouth_c, mouth_w, mouth_h, mouth_a), MOUTH_COLOR)):
        ax.add_patch(Ellipse(c, w, h, angle=a, fill=False, edgecolor=col, lw=2.2, zorder=4))

    ax.scatter(gx, gy, color=GAZE_POINT, s=90, alpha=0.22, zorder=5)
    ax.plot(gx, gy, color=GAZE_PATH, lw=1.4, zorder=5)
    ax.scatter(*fixation, s=260, alpha=0.85, facecolors="none",
               edgecolors=[FIXATION], linewidths=1.8, zorder=6)
    ax.text(fixation[0] + 16, fixation[1] + 8, "f014", color=FIXATION,
            fontsize=10, zorder=6)

    aois = (("face", face_c, face_w, face_h, face_a),
            ("eyes", eyes_c, eyes_w, eyes_h, eyes_a),
            ("mouth", mouth_c, mouth_w, mouth_h, mouth_a))
    last_gaze = (float(gx[-1]), float(gy[-1]))
    label_x = fw - 165
    for n, (name, c, w, h, a) in enumerate(aois):
        for prefix, pt, y0 in (("gaze", last_gaze, 40), ("fix", fixation, 172)):
            hit = inside(pt, c, w, h, a)
            ax.text(label_x, y0 + n * 34, f"{prefix}_{name}",
                    color="#31c452" if hit else "#f2545b", fontsize=13,
                    fontweight="bold", zorder=7)

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
