import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def plot_multi_hist(data_sets, number_of_bins, labels, xlabel="Data sets", ylabel="Data values", title='', hist_range=None):
    # Computed quantities to aid plotting
    if hist_range is None:
        hist_range = (np.nanmin(data_sets), np.nanmax(data_sets))
    binned_data_sets = [
        np.histogram(d, range=hist_range, bins=number_of_bins)[0]
        for d in data_sets
    ]
    binned_maximums = np.max(binned_data_sets, axis=1)
    # x_locations = np.concatenate((np.array([0]), np.cumsum(binned_maximums)[:-1]), axis=0)
    x_locations = np.cumsum(binned_maximums)/2.0
    x_locations[1:] += np.cumsum(binned_maximums[:-1])/2.0
    # The bin_edges are the same for all histograms
    bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
    centers = 0.5 * (bin_edges + np.roll(bin_edges, 1))[1:]
    heights = np.diff(bin_edges)

    # Cycle through and plot each histogram
    fig, ax = plt.subplots()
    for x_loc, binned_data in zip(x_locations, binned_data_sets):
        lefts = x_loc - 0.5 * binned_data
        ax.barh(centers, binned_data, height=heights, left=lefts)

    ax.set_xticks(x_locations)
    ax.set_xticklabels(labels)

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.title(title)
