from .frames_analytical import frames_analytical
from .frames_peakfinding import frames_peakfinding
from .stitch_histogram import stitch_histogram
from .stitch_events import stitch_events


def get_frame_parameters(data=None, instrument=None, plot=False):
    if data is not None:
        return frames_peakfinding(data=data, instrument=instrument, plot=plot)
    else:
        return frames_analytical(instrument=instrument, plot=plot)

def stitch(x=None, y=None, files=None, frame_params=None, **kwargs):
    if y is not None:
        return stitch_histogram(x=x, y=y, frame_params=frame_params, **kwargs)
    if events is not None:
        return stitch_events(events=x, frame_params=frame_params, **kwargs)
    # if files is not None:
    #     return stitch_files(files=files, **kwargs)
