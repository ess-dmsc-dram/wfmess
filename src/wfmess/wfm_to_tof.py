"""
  Convert WFM events to time-of-flight

  This script takes in a V20 data file which contains a WFM run that has already
  been processed by the pulse aggregator (https://github.com/ess-dmsc/pulse-aggregator)
  and converts the neutron event timings to real time-of-flight.

  Author: Neil Vaytet
  Email: neil.vaytet@esss.se
  Date: 05/2019

  Usage:

    - python convert_wfm_events_to_tof.py -i nicos_00000185.nxs -o nicos_00000185_tof.nxs

  Options (*=required):

   * -i , --input : type=str : Input file to convert
   * -o , --output : type=str : Output file name
     -e , --entries : type=str , default="entry/event_data/,entry/monitor_1/" :
                      Entries to process, separated by commas
     -p , --plot : default=False : Plot histograms?
     -g , --gauss : type=int , default=4 : Kernel width for Gaussian smoothing
     -n , --nbins : type=int , default=500 : Number of bins for histogramming
     -k , --peak-prominence : type=float , default=0.04 : Peak prominence for
                              scipy.signal.find_peaks
     -t , --inter-frame-threshold : type=float , default=0.2 : Maximum value
                                    for the signal between frames. If the
                                    signal exceeds this value, then we assume
                                    this is not a true frame gap but a
                                    fluctuation of the data inside a frame

  Example:

    - python convert_wfm_events_to_tof.py -i nicos_00000185.nxs
        -o nicos_00000185_tof.nxs -p -e entry/raw_event_data -g 6 -n 1000
        -k 0.08 -t 0.1

"""
import h5py
import numpy as np
# import argparse
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks
# from scipy.ndimage import gaussian_filter1d
from shutil import copyfile
from os.path import join
from . import v20
from .frames_analytical import frames_analytical
from .frames_peakfinding import frames_peakfinding


# def _tof_shifts(pscdata, psc_frequency=0):
#     cut_out_centre = np.reshape(pscdata,(len(pscdata)//2, 2)).mean(1)
#     cut_out_diffs = np.ediff1d(cut_out_centre)
#     return cut_out_diffs / (360.0*psc_frequency)


# def get_frame_shifts(input_file, chopper_entries=None):
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


# def make_frame_shifts(initial_shift, other_shifts=None):
#     """
#     This function constructs a list of proper shifting parameters for wavelength frames.
#     It expects the initial shift in microseconds, followed by an optional list/tuple of
#     shifts relative to the previous one. For example:

#         frame_shifts = make_frame_shifts(-6000, [-2000, -2000])
#         # results in -6000,-8000,-10000

#     The default other_shifts are the currently correct shifts for the V20 beamline at HZB.

#     :param initial_shift: Shift to apply to first wavelength frame in microseconds.
#     :param other_shifts: Iterable of shift increments with respect to the previous value
#                           (initial_shift for first value in list), also in microseconds.
#     :return: List of absolute frame shifts in microseconds.
#     """
#     frame_shift_increments = [initial_shift] + list(other_shifts)
#     frame_shifts = [sum(frame_shift_increments[:i + 1]) for i in
#                     range(len(frame_shift_increments))]

#     print("The frame_shifts are:", frame_shifts)

#     return frame_shifts


def get_frame_number(x, frames):
    """
    Find the number of the frame the input x (TOF) coordinate belongs to.
    We define here a gap_width which represents the size of the gaps between
    the frames. Any event inside a gap is discarded.
    """
    n = 0
    # Half-width of the gap between frames in microseconds
    gap_width = 800.0
    for f in frames:
        if x > f:
            n += 1
        if np.abs(x - f) < gap_width:
            n = -1
            break
    return n


def convert(input_file, entry=None, plot=False, frame_gaps=None, frame_shifts=None, nbins=5000):
    """
    This is the main function that performs the time-of-flight conversion on
    the file data.

    The options are:

     -plot : default=False : Plot histograms?
     -gauss : type=int , default=4 : Kernel width for Gaussian smoothing
     -nbins : type=int , default=500 : Number of bins for histogramming
     -peak_prominence : type=float , default=0.04 : Peak prominence for
                        scipy.signal.find_peaks
     -inter_frame_threshold : type=float , default=0.2 : Maximum value
                              for the signal between frames. If the
                              signal exceeds this value, then we assume
                              this is not a true frame gap but a
                              fluctuation of the data inside a frame

    The function first histograms the events in the file according to their
    time-of-arrival (TOA). This histogramming is simply used to find the
    boundaries between the frames via some signal processing.
    The histogrammed data is then smoothed with a gaussian filter to help the
    peak finding algorithm.

    The gaps between the frames are then found with the scipy.signal.find_peaks
    algorithm which finds the most prominent valleys in the data. Peak
    detection is often a tricky business and find_peaks will often return too
    many peaks. We go through the peaks and discard any who are lying in the
    middle of a frame (these are identified as having a y-value much higher
    than the signal background). We also discard any peaks that lie in the
    noise before the real neutron signal.

    The frame shifts are then computed from the chopper phase values.

    Finally, we go through each pulse, and then each event inside that pulse,
    find which frame that event belongs to (maybe discard it?) and shift its
    time of flight accordingly. We also re-construct the event_index array.
    If an event is discarded, the event_index is not incremented.

    To speed up the process, instead of searching through all the frames for
    each event to find which frame it belongs to and apply the correct frame
    shift, we first construct a high-resolution, uniformly spaced binning array
    in TOA, and store inside each bin the frame number that bin belongs to.
    Then, for each event we can simply get the frame number by binning the
    event's TOA into the high-resolution array.

    The real event TOFs are then saved into the output file.

    If plotting was enabled, a figure for each entry is saved as a pdf. The
    file names are the entry strings.
    """

    # Get time offsets from file
    time_offset = np.array(input_file[entry + "event_time_offset"][...],
                           dtype=np.float64, copy=True) / 1.0e3
    # Find min and max in time-of-arrival (TOA)
    xmin = np.amin(time_offset)
    xmax = np.amax(time_offset)
    nevents = len(time_offset)
    print("Read in {} events.".format(nevents))
    print("The minimum and maximum tof values are:", xmin, xmax)
    # Add padding
    dx = (xmax - xmin) / float(nbins)
    xmin -= dx
    xmax += dx
    frames_x = np.linspace(xmin, xmax, nbins + 1)

    if plot:
        # Histogram the events for plotting
        # edges = np.linspace(xmin, xmax, 513)
        y, edges = np.histogram(time_offset, bins=np.linspace(xmin, xmax, 513))
        x = 0.5 * (edges[:-1] + edges[1:])
        fig, ax = plt.subplots(2, 1, figsize=(8, 8))
        # ax[0] = fig.add_subplot(211)
        ax[0].plot(x, y, color='k', label="Raw data")
        individual_frames = [ [] for _ in range(len(frame_shifts)) ]
        for i in range(len(frame_gaps)):
            lab = "Frame boundaries" if i == 0 else None
            ax[0].axvline(frame_gaps[i], ls="dashed", color="r", label=lab)
        # xmin = np.amin(x)
        # xmax = np.amax(x)
        # dx = 0.05 * (xmax - xmin)
        # xmin -= dx
        # xmax += dx
    # # Gaussian smooth the data
    # y = gaussian_filter1d(y, gsmooth)
    # if plot:
    #     ax[0].plot(x, y, color="cyan", label="Gaussian smoothed")


#     # Minimum and maximum values in the histogram
#     ymin = np.amin(y)
#     ymax = np.amax(y)

#     # Find valleys with scipy.signal.find_peaks by inverting the y array
#     peaks, params = find_peaks(-y, prominence=peak_prominence*(ymax-ymin),
#                                distance=nbins//20)

#     # Check for unwanted peaks:
#     # 1. If a peak is found inside the noise leading the signal, it will often
#     # have a zero left base
#     to_be_removed = []
#     for p in range(len(peaks)):
#         if params["left_bases"][p] < nbins//10:
#             print("Removed peak number {} at x,y = {},{} because of a bad left base".format(p,x[peaks[p]],y[peaks[p]]))
#             to_be_removed.append(p)
#     # 2. If there are peaks in the middle of a frame, then the y value is high
#     threshold = inter_frame_threshold * (ymax - ymin) + ymin
#     for p in range(len(peaks)):
#         if y[peaks[p]] > threshold:
#             print("Removed peak number {} at x,y = {},{} because the y value exceeds the threshold {}".format(p, x[peaks[p]], y[peaks[p]], threshold))
#             to_be_removed.append(p)
#     # Actually remove the peaks
#     peaks = np.delete(peaks, to_be_removed)

#     frame_boundaries = x[peaks]
#     print("The frame_boundaries are:", frame_boundaries)
# #    print("The frame_shifts are:", frame_shifts)
#     if plot:
#         individual_frames = lists = [ [] for _ in range(len(frame_shifts)) ]
#         for i, p in enumerate(peaks):
#             lab = "Frame boundaries" if i == 0 else None
#             ax[0].axvline(x[p], ls="dashed", color="r", label=lab)
#         fig.savefig(entry.replace("/", "-")+".pdf", bbox_inches="tight")
#     if len(frame_shifts) != len(frame_boundaries) + 1:
#         print("Error: expected {} frames from chopper settings, "
#               "detected {}.".format(len(frame_shifts),
#                                     len(frame_boundaries) + 1))
#         return


    # Make a high resolution array of bins that contain the frame number they belong to
    # frames_x = np.linspace(xmin, xmax, nbins * 5)
    # print(frames_x)
    frame_id = np.array([get_frame_number(x, frame_gaps) for x in frames_x])
    frame_bin_width = frames_x[1] - frames_x[0]

    # Get the data from nexus file
    event_id = np.array(input_file[entry + "event_id"][...], dtype=np.uint32, copy=True)
    event_index = np.array(input_file[entry + "event_index"][...], dtype=np.uint64, copy=True)
    event_time_zero = np.array(input_file[entry + "event_time_zero"][...], dtype=np.uint64, copy=True)

    # Allocate new arrays
    tof = np.zeros([nevents], dtype=np.float64)
    ids = np.zeros_like(event_id)
    new_event_index = np.zeros_like(event_index)

    idx = 0
    idx_glob = 0
    event_counter = 0

    # Now go through each pulse, and then each event inside that pulse
    # Find which frame that event belongs to (maybe discard it?) and shift its
    # time of flight accordingly. If the event is discarded, the event_index is
    # not incremented
    for i in range(len(event_index) - 1):
        new_event_index[i] = idx
        for j in range(event_index[i], event_index[i+1]):
            fid = frame_id[int((time_offset[j] - xmin) / frame_bin_width)]
            # If the frame id is -1, then we are in a frame gap and event is discarded
            if fid > -1:
                ttt = time_offset[j] + frame_shifts[fid]
                if ttt > 0.0:
                    tof[idx] = ttt
                    ids[idx] = event_id[j]
                    idx += 1
                    if plot:
                        individual_frames[fid].append(ttt)

    # Discard trailing zeros
    tof = tof[:idx]
    ids = ids[:idx]

    # Update event_index entry
    input_file[entry + "event_index"][...] = new_event_index
    # Delete old event_id and event_time_offset
    del input_file[entry + "event_id"]
    del input_file[entry + "event_time_offset"]
    # Create new event_id and event_time_offset
    event_id_ds = input_file[entry].create_dataset('event_id', ids.shape, data=ids,
                                                          compression='gzip',
                                                          compression_opts=1)
    event_offset_ds = input_file[entry].create_dataset('event_time_offset', tof.shape, data=np.array(tof * 1.0e3, dtype=np.uint32),
                                                          compression='gzip',
                                                          compression_opts=1)
    event_offset_ds.attrs.create('units', np.array('ns').astype('|S2'))

    if plot:
        # ax[1] = fig.add_subplot(212)
        y, _ = np.histogram(tof, bins=edges)
        x = 0.5 * (edges[:-1] + edges[1:])
        ax[1].plot(x, y, lw=3, color="k", label="Stitched data")
        bin_width = edges[1] - edges[0]
        for i, frame in enumerate(individual_frames):
            y, e = np.histogram(frame, bins=np.arange(np.amin(frame), np.amax(frame), bin_width, dtype=np.float64))
            x = 0.5 * (e[:-1] + e[1:])
            # ax[1].plot(x, y)
            ax[1].fill_between(x, y, alpha=0.5, label="Frame {}".format(i+1))
        ax[0].grid(True, color='gray', linestyle="dotted")
        ax[0].set_axisbelow(True)
        ax[1].grid(True, color='gray', linestyle="dotted")
        ax[1].set_axisbelow(True)
        ax[0].set_xlim([xmin, xmax])
        ax[1].set_xlim([xmin, xmax])
        # dy = 0.1*(ymax - ymin)
        # ymin = 0.0
        # ymax += dy
        # ax[0].set_ylim([ymin, ymax])
        # ax[1].set_ylim([ymin, ymax])
        ax[0].set_xlabel("Raw time (microseconds)")
        ax[1].set_xlabel("Time-of-flight (microseconds)")
        ax[0].set_ylabel("Counts")
        ax[1].set_ylabel("Counts")
        ax[0].legend()
        ax[1].legend()
        fig.savefig(entry.replace("/", "-")+".pdf", bbox_inches="tight")

    return


def get_frame_parameters(data=None, instrument=None, plot=False):
    if data is not None:
        return frames_peakfinding(data=data, instrument=instrument, plot=plot)
    else:
        return frames_analytical(instrument=instrument, plot=plot)


def to_tof(input_files=None, entries=None, plot=False):


    # Compute WFM frame shifts and boundaries from TOF diagram
    v20setup = v20.setup()
    frame_boundaries, frame_gaps, frame_shifts = get_frame_parameters(v20setup)
        # v20setup["info"], v20setup["choppers"])

    if isinstance(input_files, str):
        input_files = input_files.split(",")
    elif not isinstance(input_files, list):
        input_files = [input_files]

    for f in input_files:

        ext = ".{}".format(f.split(".")[-1])
        outfile = f.replace(ext, "_tof" + ext)

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
            for e in entries:
                print("==================")
                print("Entry:", e)
                convert(input_file=outf, entry=e, plot=plot,
                               frame_gaps=frame_gaps,
                               frame_shifts=frame_shifts)

