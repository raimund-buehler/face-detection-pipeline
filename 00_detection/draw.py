from matplotlib import pyplot
from regions import EllipsePixelRegion, PixCoord
from astropy.coordinates import Angle
from math import atan2, degrees
import numpy as np


def draw_text(gaze, text, Y):
    try:
        if True in gaze:
            pyplot.text(980, Y, text, color="green")
        else:
            pyplot.text(980, Y, text, color="red")
    except TypeError:
        if gaze:
            pyplot.text(980, Y, text, color="green")
        else:
            pyplot.text(980, Y, text, color="red")


def draw_ellipsis(x, y, width, height, ax, color, face=False):
    if face:
        center = PixCoord(x=x + width / 2, y=y + height / 2)
        ellipse = EllipsePixelRegion(center, width, height * 1.2)
        ellipse.plot(ax=ax, color=color, fill=False)
    else:
        diff = np.subtract(x, y)
        center = np.add(y, diff / 2)
        center = PixCoord(center[0], center[1])
        angle = degrees(atan2(diff[1], diff[0]))
        angle = Angle(angle, "deg")
        ellipse = EllipsePixelRegion(center, width, height * 0.3, angle)
        ellipse.plot(ax=ax, color=color, fill=False)
    return ellipse
