# Detection demo

Generates `docs/detection-demo.png`, the annotated-frame figure in the top-level
README.

```
python 00_detection/demo/make_demo_figure.py -o docs/detection-demo.png
```

Frames from the study itself cannot be published — the scene camera records the
participant's conversation partner, with the participant's gaze drawn over it. So
the demo runs on a stand-in image instead. The face and its five keypoints are
detected for real; only the gaze path and the fixation are invented.

## Sample image

`sample_frame.jpg` is a synthetic portrait produced by an image generator. It
depicts no real person, so there is no subject whose likeness or privacy is at
stake, and nobody holds a photographer's copyright in it. It was cropped to 16:9
and resized to 1280×720.

## Detection

`sample_frame_detection.json` holds the committed bounding box and keypoints, so
the figure reproduces with only numpy, matplotlib and Pillow installed.

Pass `--detect` to run detection again. The script prefers MTCNN, the detector the
pipeline itself uses (`../requirements.txt`), and falls back to OpenCV's YuNet if
a `yunet.onnx` model is present in this directory. The committed detection came
from YuNet, which returns the same five keypoints as MTCNN — right eye, left eye,
nose, right mouth corner, left mouth corner — at a confidence of 0.89.

The AOI geometry is taken from `../draw.py` unchanged: the face ellipse spans the
bounding box at 1.2× its height, and the eye and mouth ellipses span their two
keypoints using the box width and 0.3× its height, rotated to the angle between
them.
