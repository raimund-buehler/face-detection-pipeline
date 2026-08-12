
#from pathlib import Path
#import matplotlib.pyplot as plt
#import cv2


def gaze_for_frame(gaze, frame, index):
    #frame_index_path = Path(frame_path).joinpath(f"frame{str(frame_index).rjust(6, '0')}.png")
    #assert frame_index_path.is_file(), f"Can't find frame image at path: {frame_index_path}"
    # Get the array of normalized gaze points for the given index (-1 to account for different naming scheme of frame and datafile)
    gaze_points = gaze[gaze["world_index"] == index]
    gaze_points = gaze_points.sort_values(by="gaze_timestamp")
    gaze_tmstmps_frame = gaze_points["gaze_timestamp"]
    gaze_conf = gaze_points["confidence"]
    gaze_points = gaze_points[["norm_pos_x", "norm_pos_y"]]
    gaze_points = gaze_points.to_numpy()
    # Split gaze points into separate X and Y coordinate arrays
    X, Y = gaze_points[:, 0], gaze_points[:, 1]
    # Flip the fixation points
    # from the original coordinate system,
    # where the origin is at botton left,
    # to the image coordinate system,
    # where the origin is at top left
    Y = 1 - Y
    # Denormalize gaze points within the frame
    # frame_index_image = cv2.cvtColor(cv2.imread(str(frame_index_path)), cv2.COLOR_BGR2RGB)
    H, W = frame.shape[:-1]
    X, Y = X * W, Y * H
    return X, Y, gaze[gaze["world_index"] == index], gaze_tmstmps_frame, gaze_conf

def fix_for_frame(fix, frame, index):
    # frame_index_path = Path(frame_path).joinpath(f"frame{str(index).rjust(6, '0')}.png")
    fix_point = fix[(fix["start_frame_index"] <= index) & (fix["end_frame_index"] > index)]
    if not fix_point.empty:
        id = fix_point.iloc[0]["id"]
        conf = fix_point.iloc[0]["confidence"]
        dur = fix_point.iloc[0]["duration"]
        fix_point = fix_point[["norm_pos_x", "norm_pos_y"]]
        fix_point = fix_point.to_numpy()
        X, Y = fix_point[:, 0], fix_point[:, 1]
        Y = 1 - Y
        #frame_index_image = plt.imread(frame_index_path)
        H, W = frame.shape[:-1]
        X, Y = X * W, Y * H
    else:
        X = False
        Y = False
        id = False
        conf = False
        dur = False
    return X, Y, id, conf, dur

def blink_for_frame(blinks, index):
    blink = blinks[(blinks["start_frame_index"] <= index - 1) & (blinks["end_frame_index"] > index - 1)]
    if blink.empty:
        return False
    else:
        return True