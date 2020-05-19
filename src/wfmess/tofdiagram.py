import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def get_frame_parameters(info, choppers):
    # Set-up the V20 choppers
    # info, choppers = v20.setup()

    # Seconds to microseconds
    microseconds = 1.0e6

    # # Number of frames
    # info["nframes"] = 6

    # # Length of pulse
    # pulse_length = 2.86e-03

    # # Position of detector
    # # detector_position = 28.98 # 32.4
    # detector_position = 28.42 # 32.4
    # # # Monitor
    # # detector_position = 25

    # # Midpoint between WFM choppers which acts as new source distance for stitched data
    # wfm_choppers_midpoint = 0.5 * (choppers["WFM1"].distance + choppers["WFM2"].distance)

    # Frame colors
    colors = ['b', 'k', 'g', 'r', 'cyan', 'magenta']

    # Make figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.grid(True, color='lightgray', linestyle="dotted")
    ax.set_axisbelow(True)

    # Plot the chopper openings
    for key, ch in choppers.items():
        dist = [ch.distance, ch.distance]
        for i in range(0, len(ch.openings), 2):
            t1 = (ch.openings[i] + ch.phase) / ch.omega * microseconds
            t2 = (ch.openings[i+1] + ch.phase) / ch.omega * microseconds
            ax.plot([t1, t2], dist, color=colors[i//2])
        ax.text(t2 + (t2-t1), ch.distance, ch.name, ha="left", va="center")


    # Define and draw source pulse
    x0 = 0.0
    x1 = info["pulse_length"] * microseconds
    y0 = 0.0
    y1 = 0.0
    psize = info["detector_position"]/50.0
    rect = Rectangle((x0, y0), x1, -psize, lw=1, fc='orange', ec='k', hatch="////", zorder=10)
    ax.add_patch(rect)
    ax.text(x0, -psize, "Source pulse (2.86 ms)", ha="left", va="top", fontsize=6)

    # Now find frame boundaries and draw frames
    frame_boundaries = []
    frame_shifts = []

    for i in range(info["nframes"]):

        # Find the minimum and maximum slopes that are allowed through each frame
        slope_min = 1.0e30
        slope_max = -1.0e30
        for key, ch in choppers.items():

            # For now, ignore Wavelength band double chopper
            if len(ch.openings) == info["nframes"] * 2:

                xmin = (ch.openings[i*2] + ch.phase) / ch.omega * microseconds
                xmax = (ch.openings[i*2+1] + ch.phase) / ch.omega * microseconds
                slope1 = (ch.distance - y1) / (xmin - x1)
                slope2 = (ch.distance - y0) / (xmax - x0)

                if slope_min > slope1:
                    x2 = xmin
                    y2 = ch.distance
                    slope_min = slope1
                if slope_max < slope2:
                    x3 = xmax
                    y3 = ch.distance
                    slope_max = slope2

        # Compute line equation parameters y = a*x + b
        a1 = (y3 - y0) / (x3 - x0)
        a2 = (y2 - y1) / (x2 - x1)
        b1 = y0 - a1 * x0
        b2 = y1 - a2 * x1

        y4 = info["detector_position"]
        y5 = info["detector_position"]

        # This is the frame boundaries
        x5 = (y5 - b1)/a1
        x4 = (y4 - b2)/a2
        frame_boundaries.append([x4, x5])

        # Compute frame shifts from fastest neutrons in frame
        frame_shifts.append(-(info["wfm_choppers_midpoint"] - b2)/a2)

        ax.fill([x0, x1, x4, x5], [y0, y1, y4, y5], alpha=0.3, color=colors[i])
        ax.plot([x0, x5], [y0, y5], color=colors[i], lw=1)
        ax.plot([x1, x4], [y1, y4], color=colors[i], lw=1)
        ax.text(0.5*(x4+x5), info["detector_position"], "Frame {}".format(i+1), ha="center", va="top")


    # Plot detector location
    ax.plot([0, np.amax(frame_boundaries)], [info["detector_position"], info["detector_position"]], lw=3, color='grey')
    ax.text(0.0, info["detector_position"], "Detector", va="bottom", ha="left")
    # Plot WFM choppers mid-point
    ax.plot([0, np.amax(frame_boundaries)], [info["wfm_choppers_midpoint"], info["wfm_choppers_midpoint"]], lw=1, color='grey', ls="dashed")
    ax.text(np.amax(frame_boundaries), info["wfm_choppers_midpoint"], "WFM chopper mid-point", va="bottom", ha="right")

    # Print results
    print("The frame boundaries are:", frame_boundaries)
    frame_gaps = [0.5*(frame_boundaries[i][1]+frame_boundaries[i+1][0]) for i in range(len(frame_boundaries)-1)]
    print("The frame gaps are:", frame_gaps)
    print("The frame shifts are:", frame_shifts)

    # Save the figure
    ax.set_xlabel("Time [microseconds]")
    ax.set_ylabel("Distance [m]")
    fig.savefig("tof_diagram.pdf", bbox_inches="tight")

    return np.array(frame_boundaries), np.array(frame_gaps), np.array(frame_shifts)
