# Detection demo

Generates `docs/detection-demo.png`, the annotated-frame figure in the top-level
README.

```
python 00_detection/demo/make_demo_figure.py -o docs/detection-demo.png
```

Frames from the study itself cannot be published — the scene camera records the
participant's conversation partner, with the participant's gaze drawn over it. So
the demo runs on a stock portrait instead. The face and its five keypoints are
detected for real; only the gaze path and the fixation are invented.

## Sample image

| | |
|---|---|
| File | `sample_frame.jpg` |
| Source | [Face portrait (Unsplash)](https://commons.wikimedia.org/wiki/File:Face_portrait_(Unsplash).jpg) on Wikimedia Commons |
| Author | William Stitt |
| Licence | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) — public domain dedication, no rights reserved |
| Changes | cropped to 16:9 and resized to 1280×720 |

The person shown is not a study participant.

## Detection

`sample_frame_detection.json` holds the committed bounding box and keypoints, so
the figure reproduces with only numpy, matplotlib and Pillow installed.

Pass `--detect` to run detection again. The script prefers MTCNN, the detector the
pipeline itself uses (`../requirements.txt`), and falls back to OpenCV's YuNet if
a `yunet.onnx` model is present in this directory. The committed detection came
from YuNet, which returns the same five keypoints as MTCNN — right eye, left eye,
nose, right mouth corner, left mouth corner — at a confidence of 0.69.

The AOI geometry is taken from `../draw.py` unchanged: the face ellipse spans the
bounding box at 1.2× its height, and the eye and mouth ellipses span their two
keypoints using the box width and 0.3× its height, rotated to the angle between
them.
