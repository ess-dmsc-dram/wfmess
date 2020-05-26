from .frames_analytical import frames_analytical
from .frames_peakfinding import frames_peakfinding

def get_frames(data=None, instrument=None, plot=False):
    if data is not None:
        return frames_peakfinding(data=data, instrument=instrument, plot=plot)
    else:
        return frames_analytical(instrument=instrument, plot=plot)
