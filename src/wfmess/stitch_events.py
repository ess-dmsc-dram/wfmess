import numpy as np
from .tools import edges_to_centers


def _get_frame_number(x, ledges, redges):
    """
    Find the number of the frame the input x (TOF) coordinate belongs to.
    Return -1 for any event not located inside a frame.
    """
    n = -1
    for i in range(len(ledges)):
        if x >= ledges[i] and x < redges[i]:
            n = i
            break
    return n

def stitch_events(events=None, frames=None, plot=False, nbins=5000):
    """
    This function takes in a list of events
    """

    # Get time offsets
    time_offset = events["tof"]

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

    if plot:
        # Histogram the events for plotting
        y, edges = np.histogram(time_offset, bins=np.linspace(xmin, xmax, 513))
        x = 0.5 * (edges[:-1] + edges[1:])
        fig, ax = plt.subplots(2, 1, figsize=(8, 8))
        ax[0].plot(x, y, color='k', label="Raw data")
        for g in frames["gaps"]:
            ax[0].axvline(x=g, color="r")


    # Make a high resolution array of bins that contain the frame number they belong to
    frames_x = edges_to_centers(np.linspace(xmin, xmax, nbins + 1))
    frame_id = np.array(
        [get_frame_number(
            x, frames["left_edges"], frames["right_edges"])
         for x in frames_x])
    frame_bin_width = frames_x[1] - frames_x[0]

    # Allocate new arrays
    stitched = {"tof": np.zeros([nevents], dtype=np.float64)}

    # Allocate optional arrays
    process_ids = False
    process_index = False
    if "ids" in events:
        event_id = events["ids"]
        stitched["ids"] = np.zeros_like(event_id)
        process_ids = True
    if "index" in events:
        event_index = events["index"]
        stitched["index"] = np.zeros_like(event_index)
        process_index = True


    idx = 0
    idx_glob = 0
    event_counter = 0

    # Now go through each pulse, and then each event inside that pulse
    # Find which frame that event belongs to (maybe discard it?) and shift its
    # time of flight accordingly. If the event is discarded, the event_index is
    # not incremented
    for i in range(len(event_index) - 1):
        if process_index:
            stitched["index"][i] = idx
        for j in range(event_index[i], event_index[i+1]):
            fid = frame_id[int((time_offset[j] - xmin) / frame_bin_width)]
            # If the frame id is -1, then we are in a frame gap and event is discarded
            if fid > -1:
                ttt = time_offset[j] + frames["shifts"][fid]
                if ttt > 0.0:
                    stitched["tof"][idx] = ttt
                    if process_ids:
                        stitched["ids"][idx] = event_id[j]
                    idx += 1
                    # if plot:
                    #     individual_frames[fid].append(ttt)

    # Discard trailing zeros
    stitched["tof"] = stitched["tof"][:idx]
    if process_ids:
        stitched["ids"] = stitched["ids"][:idx]


    if plot:
        y, _ = np.histogram(stitched["tof"], bins=edges)
        x = 0.5 * (edges[:-1] + edges[1:])
        ax[1].plot(x, y, lw=3, color="k", label="Stitched data")
        ax[0].grid(True, color='gray', linestyle="dotted")
        ax[0].set_axisbelow(True)
        ax[1].grid(True, color='gray', linestyle="dotted")
        ax[1].set_axisbelow(True)
        ax[0].set_xlim([xmin, xmax])
        ax[1].set_xlim([xmin, xmax])
        ax[0].set_xlabel("Raw time (microseconds)")
        ax[1].set_xlabel("Stitched Time-of-flight (microseconds)")
        ax[0].set_ylabel("Counts")
        ax[1].set_ylabel("Counts")
        ax[0].legend()
        ax[1].legend()
        if isinstance(plot, str):
            figname = plot
        else:
            figname = "stitched_events.pdf"
        fig.savefig(figname, bbox_inches="tight")

    return stitched
