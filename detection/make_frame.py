import matplotlib
from matplotlib import pyplot

#matplotlib.use('agg')
import numpy as np
from draw import draw_text, draw_ellipsis
from regions import PixCoord
import io


def make_frame(frame, keypoint, box, gaze_X, gaze_Y, X_fix, Y_fix, id_fix):
    # plot the image
    pyplot.figure(1, figsize=[10.80, 7.20])
    pyplot.imshow(frame)
    pyplot.axis("off")
    # get the context for drawing boxes
    ax = pyplot.gca()
    fig = pyplot.gcf()

    # if not result_list:
    # 	gaze_on_face = gaze_on_eyes = gaze_on_mouth = fig = ax = 0

    # draw gaze
    pyplot.scatter(gaze_X, gaze_Y, color=(1.0, 0.7, 0.25), s=100, alpha=0.2)
    pyplot.plot(gaze_X, gaze_Y, color=(0.0, 1.0, 0.4), lw=1)

    gazes = np.array(["gaze_face", "gaze_eyes", "gaze_mouth"])
    fixes = np.array(["fix_face", "fix_eyes", "fix_mouth"])

    fix_on_face = fix_on_eyes = fix_on_mouth = False
    fix_obj = [fix_on_face, fix_on_eyes, fix_on_mouth]

    # draw fix
    if not X_fix:
        fix_present = False
    else:
        fix_present = True
        pyplot.scatter(X_fix, Y_fix, color=(0.0, 1.0, 1.0), s=200, alpha=0.5, facecolors='none')
        pyplot.text(X_fix, Y_fix, id_fix, color=(0.0, 1.0, 1.0))

    if box.size == 0:
        for n in range(3):
            draw_text(False, gazes[n], 40 + n * 40)
            draw_text(False, fixes[n], 160 + n * 40)
        face_present = False
    else:
        face_present = True

    x = box[0]
    y = box[3]
    width = box[2]-x
    height = box[1]-y
    #x, y, width, height = result['box']

    # draw ellipsis around face
    ellipse_face = draw_ellipsis(x, y, width, height, ax, "orange", face=True)

    # draw ellipsis around eyes
    right_eye = keypoint['right_eye']
    left_eye = keypoint['left_eye']

    ellipse_eyes = draw_ellipsis(right_eye, left_eye, width, height, ax, "red")

    # draw ellipsis around mouth
    mouth_right = keypoint['mouth_right']
    mouth_left = keypoint['mouth_left']

    ellipse_mouth = draw_ellipsis(mouth_right, mouth_left, width, height, ax, "blue")

    # test gaze on regions
    gaze = PixCoord(gaze_X, gaze_Y)
    ellipses = np.array([ellipse_face, ellipse_eyes, ellipse_mouth])

    for n, el in enumerate(ellipses):
        if any(el.contains(gaze)):
            draw_text(True, gazes[n], 40 + n * 40)
        else:
            draw_text(False, gazes[n], 40 + n * 40)

    # test fix on regions
    if X_fix:
        fix_coord = PixCoord(X_fix, Y_fix)
        for n, el in enumerate(ellipses):
            if el.contains(fix_coord):
                draw_text(True, fixes[n], 160 + n * 40)
                fix_obj[n] = True
            else:
                draw_text(False, fixes[n], 160 + n * 40)
    else:
        for n in range(3):
            draw_text(False, fixes[n], 160 + n * 40)

    #print(pyplot.get_fignums())

    # io_buf = io.BytesIO()
    # fig.savefig(io_buf, format='raw')
    # io_buf.seek(0)
    # X = np.reshape(np.frombuffer(io_buf.getvalue(), dtype=np.uint8),
    #                      newshape=(int(fig.bbox.bounds[3]), int(fig.bbox.bounds[2]), -1))
    # io_buf.close()

    fig.tight_layout(pad=0)
    fig.canvas.draw()

    X = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    X = X.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    # X = X[:, np.where(~np.all(X == 255, axis=0))[0]]
    # X = X[:, np.where(~np.all(X == 255, axis=1))[0]]

    return face_present, fix_present, fix_obj[0], fix_obj[1], fix_obj[2], X
