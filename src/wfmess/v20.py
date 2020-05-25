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







# def _get_frame_shifts(choppers, chopper_entries=None):
#     """
#     This function generates a list of proper shifting parameters from psc data.

#     :param input_file: An open hdf file handle
#     :param chopper_entries: Comma-separated list of entries for WFM choppers
#     :return: List of relative frame shifts in microseconds.
#     """

#     time_units = {"us": 1.0e6, "ns": 1.0e9}

#     # Read chopper timings and compute frequencies
#     freqs = []
#     try:
#         for c in chopper_entries:
#             tdc = input_file[join(c, "top_dead_center/time")]
#             timings = np.array(tdc[...], dtype=np.float64, copy=True)
#             freqs.append(np.mean(time_units[tdc.attrs["units"].decode("utf-8")] / np.ediff1d(timings)))
#     except KeyError:
#         print("Warning: chopper timings not found in file. Using default "
#               "value for WFM chopper frequencies.")
#         freqs = [70.0] * 2
#     print("WFM chopper frequencies are:", freqs)

#     # Read chopper cut-out angles
#     angles = []
#     try:
#         for c in chopper_entries:
#             angles.append(np.array(input_file[join(c, "slit_edges")][...],
#                                    dtype=np.float64, copy=True))
#     except KeyError:
#         print("Warning: chopper cut-out angles not found in file. Using "
#               "default values.")
#         # Definition V20 for wfm pulse shaping chopper 1 (closest to source)
#         # and chopper 2 (farthest from source). The sorted angles of all edges
#         # are in degrees. First entry is start angle of the first cut-out
#         # second entry is end angle of first cut-out
#         angles = [np.array([15.0, 98.71, 109.7, 155.49, 170.79, 208.26,
#                             227.56, 257.32, 280.33, 302.91, 329.37, 345.3]),
#                   np.array([15.0, 80.04, 91.03, 141.1 ,156.4 ,197.88,
#                             217.18, 250.67, 269.97, 299.73, 322.74, 345.])]

#     # factor of 0.5 * 1.0e6 (taking mean and converting to microseconds)
#     relative_shifts = (_tof_shifts(angles[0], psc_frequency=freqs[0])  + \
#                        _tof_shifts(angles[1], psc_frequency=freqs[1])) * \
#                        5.0e+05
#     return -relative_shifts
