from facenet_pytorch import MTCNN
import argparse
from pathlib import Path
from matplotlib import pyplot
import torch
from imutils.video import FileVideoStream
import cv2
import time
from gaze_for_frame import gaze_for_frame, fix_for_frame, blink_for_frame
from make_frame import make_frame
from filter_conf import filter_conf
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--export_dir', type=str, required=True)
args = parser.parse_args()
device = 'cuda' if torch.cuda.is_available() else 'cpu'


class FastMTCNN(object):
    """Fast MTCNN implementation."""

    def __init__(self, stride, resize=1, *args, **kwargs):
        """Constructor for FastMTCNN class.

        Arguments:
            stride (int): The detection stride. Faces will be detected every `stride` frames
                and remembered for `stride-1` frames.

        Keyword arguments:
            resize (float): Fractional frame scaling. [default: {1}]
            *args: Arguments to pass to the MTCNN constructor. See help(MTCNN).
            **kwargs: Keyword arguments to pass to the MTCNN constructor. See help(MTCNN).
        """
        self.stride = stride
        self.resize = resize
        self.mtcnn = MTCNN(*args, **kwargs)

    def __call__(self, frames):
        """Detect faces in frames using strided MTCNN."""
        if self.resize != 1:
            frames = [
                cv2.resize(f, (int(f.shape[1] * self.resize), int(f.shape[0] * self.resize)))
                for f in frames
            ]

        boxes, probs, landmarks = self.mtcnn.detect(frames[::self.stride], landmarks=True)

        keypoints = {}
        keys = range(len(frames))
        for i in keys:
            index = int(i / self.stride)
            landmark_dict = {
                "left_eye": landmarks[index][0][0],
                "right_eye": landmarks[index][0][1],
                "nose": landmarks[index][0][2],
                "mouth_left": landmarks[index][0][3],
                "mouth_right": landmarks[index][0][4]
            }
            keypoints[i] = landmark_dict

        faces = []
        for i, frame in enumerate(frames):
            box_ind = int(i / self.stride)
            if boxes[box_ind] is None:
                continue
            for box in boxes[box_ind]:
                box = [int(b) for b in box]
                faces.append(frame[box[1]:box[3], box[0]:box[2]])

        return faces, keypoints, boxes


fast_mtcnn = FastMTCNN(
    stride=4,
    resize=1,
    margin=14,
    factor=0.6,
    select_largest=True,
    keep_all=False,
    device=device
)


def run_detection(fast_mtcnn, v_cap, v_len, gaze, fix, out):
    frames = []
    frames_processed = 0
    faces_detected = 0
    batch_size = 60
    start = time.time()

    #v_cap = FileVideoStream(video).start()
    #v_len = int(v_cap.stream.get(cv2.CAP_PROP_FRAME_COUNT))

    for j in range(v_len):

        frame = v_cap.read()
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        frames.append(frame)

        if len(frames) >= batch_size or j == v_len - 1:
            faces, keypoints, boxes = fast_mtcnn(frames)

            for i, frame in enumerate(frames):
                index = int(j / 60) * 60 + i
                X, Y, data, gaze_tmstmps_frame, gaze_conf = gaze_for_frame(
                    gaze=gaze,
                    frame=frame,
                    index=index)

                X_fix, Y_fix, id_fix, conf_fix, duration_fix = fix_for_frame(
                    fix=fix,
                    frame=frame,
                    index=index
                )
                keypoint = keypoints[i]
                box = boxes[int(i / 4)][0]
                face_present, fix_present, fix_on_face, fix_on_eyes, fix_on_mouth, frame = make_frame(frame, keypoint,
                                                                                                      box, X,
                                                                                                      Y, X_fix,
                                                                                                      Y_fix, id_fix)
                blink = blink_for_frame(blinks, index + 1)
                df.iloc[index] = pd.Series({
                    "frame_id": index,
                    "face_present": face_present,
                    "fix_present": fix_present,
                    "id_fix": id_fix,
                    "conf_fix": conf_fix,
                    "duration_fix": duration_fix,
                    "pos_x_fix": X_fix,
                    "pos_y_fix": Y_fix,
                    "fix_on_face": fix_on_face,
                    "fix_on_eyes": fix_on_eyes,
                    "fix_on_mouth": fix_on_mouth,
                    "gaze_tmstmps_frame": gaze_tmstmps_frame.to_numpy(),
                    "gaze_conf": gaze_conf.to_numpy(),
                    "pos_x_gaze": X,
                    "pos_y_gaze": Y,
                    "blink": blink
                })
                # pyplot.text(980, 280, index, color="white")
                # save frame
                # path_edited = Path("{FRAME_OUTPUT_DIR}/").joinpath(
                #     f"frame{str(index + 1).rjust(6, '0')}.png")
                frame = cv2.cvtColor(frame, cv2.COLOR_RGB2BGR)
                out.write(frame)
                # pyplot.savefig(path_edited, dpi=200, bbox_inches="tight", pad_inches=0)
                print("Saved frame: ", str(index + 1).rjust(6, '0'))
                pyplot.clf()
                pyplot.cla()

            frames_processed += len(frames)
            faces_detected += len(faces)
            frames = []

            print(
                f'Frames per second: {frames_processed / (time.time() - start):.3f},',
                f'faces detected: {faces_detected}\r',
                end=''
            )


# rec_dir
recording_path = args.export_dir
rec_dir = Path(recording_path)
assert rec_dir.is_dir(), "recording path not found"

# read csv data

# gaze
gaze_path = rec_dir.joinpath("exports", "000", "gaze_positions.csv")
assert gaze_path.is_file(), "gaze_positions.csv not found"

gaze = pd.read_csv(gaze_path)
# filter for confidence
gaze, n_excl_conf = filter_conf(gaze, 0.80)
print("Excluded gaze points confidence < 0.80: ", n_excl_conf)

# filter normalized gaze points for below 0 and above 1
n_excl_bounds = len(gaze) - len(
    gaze[(gaze["norm_pos_x"] > 0) & (gaze["norm_pos_x"] < 1) & (gaze["norm_pos_y"] > 0) & (gaze["norm_pos_y"] < 1)])
gaze = gaze[(gaze["norm_pos_x"] > 0) & (gaze["norm_pos_x"] < 1) & (gaze["norm_pos_y"] > 0) & (gaze["norm_pos_y"] < 1)]
print("Excluded gaze points <0 or >1:", n_excl_bounds)

# fixations
fix_path = rec_dir.joinpath("exports", "000", "fixations.csv")
assert fix_path.is_file(), "fixations.csv not found"
fix = pd.read_csv(fix_path)

# filter for confidence
fix, n_excl_fix_conf = filter_conf(fix, 0.80)
print("Excluded fixations confidence < 0.80: ", n_excl_fix_conf)

# filter normalized fix points for below 0 and above 1
n_excl_fix_bounds = len(fix) - len(
    fix[(fix["norm_pos_x"] > 0) & (fix["norm_pos_x"] < 1) & (fix["norm_pos_y"] > 0) & (fix["norm_pos_y"] < 1)])
fix = fix[(fix["norm_pos_x"] > 0) & (fix["norm_pos_x"] < 1) & (fix["norm_pos_y"] > 0) & (fix["norm_pos_y"] < 1)]
print("Excluded fixations <0 or >1:", n_excl_fix_bounds)

# blinks

blink_path = rec_dir.joinpath("exports", "000", "blinks.csv")
assert blink_path.is_file(), "blinks.csv not found"
blinks = pd.read_csv(blink_path)
# filter for confidence
blinks, n_excl_blinks_conf = filter_conf(blinks, 0.80)
print("Excluded blinks confidence < 0.80: ", n_excl_blinks_conf)

column_names = [
    "frame_id",
    "face_present",
    "fix_present",
    "id_fix",
    "conf_fix",
    "duration_fix",
    "pos_x_fix",
    "pos_y_fix",
    "fix_on_face",
    "fix_on_eyes",
    "fix_on_mouth",
    "gaze_tmstmps_frame",
    "gaze_conf",
    "pos_x_gaze",
    "pos_y_gaze",
    "blink"]

video_path = rec_dir.joinpath("exports", "000", "world.mp4")
assert video_path.is_file(), "world.mp4 not found"
video_path = str(video_path)

save_path = str(rec_dir.joinpath("exports", "000", "world_edited.mp4"))
out = cv2.VideoWriter(save_path,
                      cv2.VideoWriter_fourcc(*'mp4v'), 30,
                      (1080, 720))

v_cap = FileVideoStream(video_path).start()
v_len = int(v_cap.stream.get(cv2.CAP_PROP_FRAME_COUNT))
df = pd.DataFrame(columns=column_names, index=range(v_len))

run_detection(fast_mtcnn, v_cap, v_len, gaze, fix, out)

result_path = rec_dir.joinpath("exports", "000", "result.csv")
df.to_csv(result_path)
out.release()
