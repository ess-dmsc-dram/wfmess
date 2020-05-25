import numpy as np
from .tools import centers_to_edges, edges_to_centers

def stitch_histogram(x=None, y=None, frame_params=None):
    """
    Stitch a N-dimensional array containing histogram data.

    :param x: The array containing the tof x-coordinate.
    :param y: The N-D array containing the neutron counts to be stitched.
    :param frame_params: A dict containing the WFM frame parameters.
    """

    nx = y.shape[-1]
    if len(x) == y.shape[-1]:
        xe = centers_to_edges(x)
        xc = x
    else:
        xc = edges_to_centers(x)
        xe = x
    xmin = xe[0]
    xmax = xe[-1]
    dx = (xmax - xmin) / float(nx)

    # Make new array for stitched image
    stitched = np.zeros_like(y, dtype=np.float)

    # Go through all original image columns and shift according to frame shifts.
    # We find to which bin the left and right pixel edges would be shifted, and
    # spread the counts evenly over these bins.
    for i in range(nx):
        ok = 0
        for j in range(len(frame_params["shifts"])):
            if xe[i] > frame_params["left_edges"][j] and xe[i] < frame_params["right_edges"][j]:
                # compute new index after time shift
                idl = int((xe[i] + frame_params["shifts"][j] - xmin) / dx)
                ok += 1
                break
        for j in range(len(frame_params["shifts"])):
            if xe[i+1] > frame_params["left_edges"][j] and xe[i+1] < frame_params["right_edges"][j]:
                # compute new index after time shift
                idr = int((xe[i+1] + frame_params["shifts"][j] - xmin) / dx)
                ok += 1
                break
        if ok == 2:
            npix = float(idr - idl + 1)
            for idx in range(idl, idr + 1):
                stitched[..., idx] += y[..., i] / npix

    return stitched
