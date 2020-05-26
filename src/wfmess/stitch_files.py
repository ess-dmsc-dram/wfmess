"""

"""
import h5py
import numpy as np
from shutil import copyfile
from os.path import join
from . import v20
from .wfmess import get_frames


def _stitch_file(file_handle, entry=None, frames=None, plot=False):
    """
    
    """

    # Get time offsets from file
    events = {}
    events["tof"] = np.array(file_handle[entry + "event_time_offset"][...],
                           dtype=np.float64, copy=True) / 1.0e3

    # Get the data from nexus file
    events["ids"] = np.array(file_handle[entry + "event_id"][...], dtype=np.uint32, copy=True)
    events["index"] = np.array(file_handle[entry + "event_index"][...], dtype=np.uint64, copy=True)
    # event_time_zero = np.array(input_file[entry + "event_time_zero"][...], dtype=np.uint64, copy=True)

    # Stitch the data
    stitched = stitch_events(events=events, frames=frames, plot=plot) 

    # Update event_index entry
    file_handle[entry + "event_index"][...] = stitched["index"]
    # Delete old event_id and event_time_offset
    del file_handle[entry + "event_id"]
    del file_handle[entry + "event_time_offset"]
    # Create new event_id and event_time_offset
    event_id_ds = file_handle[entry].create_dataset(
        'event_id', stitched["ids"].shape, data=stitched["index"],
                                                          compression='gzip',
                                                          compression_opts=1)
    event_offset_ds = file_handle[entry].create_dataset(
        'event_time_offset', stitched["tof"].shape, data=np.array(stitched["tof"] * 1.0e3, dtype=np.uint32),
                                                          compression='gzip',
                                                          compression_opts=1)
    event_offset_ds.attrs.create('units', np.array('ns').astype('|S2'))


    return


def stitch_files(input_files=None, entries=None, plot=False, frames=None):

    if isinstance(input_files, str):
        input_files = input_files.split(",")
    elif not isinstance(input_files, list):
        input_files = [input_files]

    for f in input_files:

        ext = ".{}".format(f.split(".")[-1])
        outfile = f.replace(ext, "_stitched" + ext)

        copyfile(f, outfile)

        # Automatically find event entries if not specified
        if entries is None:
            entries = []
            key = "event_time_offset"
            with h5py.File(outfile, "r") as f:
                contents = []
                f.visit(contents.append)
            for item in contents:
                if item.endswith(key):
                    entries.append(item.replace(key, ""))
        else:
            entries = entries.split(",")

        # Loop through entries and shift event tofs
        with h5py.File(outfile, "r+") as outf:

            # Compute WFM frame shifts and boundaries from V20 setup
            if frames is None:
                v20setup = v20.setup(filename=outf)
                frames = get_frames(instrument=v20setup)
                # v20setup["info"], v20setup["choppers"])

            for e in entries:
                print("==================")
                print("Entry:", e)
                this_plot = plot
                if plot is True:
                    this_plot = outfile + "-" + entry + ".pdf"
                _stitch_file(file_handle=outf, entry=e,
                             frames=frames, plot=this_plot)

