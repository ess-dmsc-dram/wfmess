from .stitch_histogram import stitch_histogram
from .stitch_events import stitch_events
from .stitch_files import stitch_files


def stitch(x=None, y=None, files=None, frame_params=None, **kwargs):
    if y is not None:
        return stitch_histogram(x=x, y=y, frame_params=frame_params, **kwargs)
    if events is not None:
        return stitch_events(events=x, frame_params=frame_params, **kwargs)
    if files is not None:
        return stitch_files(files=files, **kwargs)
