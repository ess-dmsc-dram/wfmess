from .stitch_histogram import stitch_histogram
from .stitch_events import stitch_events
from .stitch_files import stitch_files


def stitch(x=None, y=None, files=None, frames=None, **kwargs):
    if y is not None:
        return stitch_histogram(x=x, y=y, frames=frames, **kwargs)
    if events is not None:
        return stitch_events(events=x, frames=frames, **kwargs)
    if files is not None:
        return stitch_files(files=files, **kwargs)
