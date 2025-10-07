import matplotlib as mpl
from matplotlib import pyplot as plt


# Extract data
def plot_linear_conc(results, substrate, time_range, interval, biomass=False):
    x_positions = [results[0]["position"][compartment][0] for compartment in results[0]["position"].keys()]
    t1=time_range[0]
    tl=time_range[1]
    timepoints = [timepoint for timepoint in
                  range(t1, tl, interval)]
    substrate_conc = [
        [results[timepoint]["compartments"][compartment]["concentrations"][substrate]
         for compartment in results[0]["compartments"].keys()]
        for timepoint in timepoints
    ]

    # Set up figure and axis
    fig, ax = plt.subplots()

    # Set up colormap and normalization
    cmap = mpl.cm.plasma
    norm = mpl.colors.Normalize(vmin=min(timepoints), vmax=max(timepoints))
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required even if not used for image

    # Plot lines
    for i, y in enumerate(substrate_conc):
        t = timepoints[i]
        ax.plot(x_positions, y, color=cmap(norm(t)))

    # Labels and colorbar
    ax.set_xlabel("x position")
    if biomass:
        ax.set_ylabel(f"biomass (gDW)")
        ax.set_title(f"{substrate} Biomass")
    else:
        ax.set_ylabel(f"concentration (mM)")
        ax.set_title(f"{substrate} Concentration")
    cbar = fig.colorbar(sm, ax=ax, label="Timepoint")
    return fig, ax