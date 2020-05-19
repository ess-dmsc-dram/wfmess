import numpy as np
from .chopper import Chopper


def setup():

    choppers = dict()
    info = dict()

    choppers["WFM1"] = Chopper(frequency=70,
                               openings=np.array([83.71, 94.7, 140.49, 155.79,
                                                   193.26, 212.56, 242.32, 265.33,
                                                   287.91, 314.37, 330.3, 360.0]) + 15.0,
                               phase=47.10,
                               distance=6.6,
                               unit="deg",
                               name="WFM1")

    choppers["WFM2"] = Chopper(frequency=70,
                               openings=np.array([65.04, 76.03, 126.1, 141.4, 182.88,
                                                  202.18, 235.67, 254.97, 284.73,
                                                  307.74, 330.00, 360.0]) + 15.0,
                               phase=76.76,
                               distance=7.1,
                               unit="deg",
                               name="WFM2")

    choppers["FOL1"] = Chopper(frequency=56,
                               openings=np.array([74.6, 95.2, 139.6, 162.8, 194.3, 216.1, 245.3, 263.1, 294.8, 310.5, 347.2, 371.6]),
                               phase=62.40,
                               distance=8.8,
                               unit="deg",
                               name="Frame-overlap 1")

    choppers["FOL2"] = Chopper(frequency=28,
                               openings=np.array([98., 134.6, 154., 190.06, 206.8, 237.01, 254., 280.88, 299., 323.56, 344.65, 373.76]),
                               phase=12.27,
                               distance=15.9,
                               unit="deg",
                               name="Frame-overlap 2")

    # Number of frames
    info["nframes"] = 6

    # Length of pulse
    info["pulse_length"] = 2.86e-03

    # Position of detector
    # detector_position = 28.98 # 32.4
    info["detector_position"] = 28.42 # 32.4
    # # Monitor
    # detector_position = 25

    # Midpoint between WFM choppers which acts as new source distance for stitched data
    info["wfm_choppers_midpoint"] = 0.5 * (choppers["WFM1"].distance + choppers["WFM2"].distance)

    return {"info": info, "choppers": choppers}
