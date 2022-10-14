import scanpy
import numpy as np
from pickle import dump, load
import matplotlib.pyplot as plt

samples_path = '../120Samples/OriginalData/'
db_name = 'adata_all.h5ad' # obs = ['celltype_sub', 'celltype', 'subject', 'timepoint', '10Xchemistry']
'''Getting cell type populations'''
# data = scanpy.read_h5ad(samples_path+db_name)
#
# donors = sorted(set(data.obs['subject']))
# cell_types = sorted(set(data.obs['celltype']))
# timepoints = set(data.obs['timepoint'])
#
# cell_frequencies = dict()
# for t in timepoints:
#     cell_frequencies[t] = np.zeros((len(donors), len(cell_types)))
#     for d in range(len(donors)):
#         total_cells = 0
#         mask = (data.obs['timepoint'] == t) & (data.obs['subject'] == donors[d])
#         for c in range(len(cell_types)):
#             n_cells = len(data[mask & (data.obs['celltype'] == cell_types[c])])
#             cell_frequencies[t][d, c] = n_cells
#             total_cells += n_cells
#         cell_frequencies[t][d, :] /= total_cells
#
# with open(samples_path+'cell_frequencies.pickle', 'wb') as handle:
#     dump(cell_frequencies, handle)
cell_types = ['B', 'CD4T', 'CD8T', 'DC', 'DNT', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B']

def plot_multi_hist(data_sets, number_of_bins, labels, xlabel="Data sets", ylabel="Data values", title=''):
    # Computed quantities to aid plotting
    plt.figure()
    hist_range = (np.min(data_sets), np.max(data_sets))
    binned_data_sets = [
        np.histogram(d, range=hist_range, bins=number_of_bins)[0]
        for d in data_sets
    ]
    binned_maximums = np.max(binned_data_sets, axis=1)
    x_locations = np.arange(0, sum(binned_maximums), np.max(binned_maximums))

    # The bin_edges are the same for all of the histograms
    bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
    centers = 0.5 * (bin_edges + np.roll(bin_edges, 1))[:-1]
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

    # plt.show()

with open(samples_path+'cell_frequencies.pickle', 'rb') as handle:
    cell_frequencies = load(handle)
assert isinstance(cell_frequencies, dict)

log_cell_frequencies = dict()
for t in cell_frequencies.keys():
    cell_frequencies[t][np.isnan(cell_frequencies[t])] = 0
    log_cell_frequencies[t] = np.log10(cell_frequencies[t]+0.0001)

stimul_type = 'PA'
kwargs = dict(histtype='step', bins=list(np.arange(-3, 0, 0.2)))
for ctype in range(log_cell_frequencies['UT'].shape[1]):
    labels =[]
    data_sets = []
    for t in reversed(sorted(log_cell_frequencies.keys())):
        if 'UT' in t or stimul_type in t:
            data_sets.append(log_cell_frequencies[t][:, ctype])
            labels.append(t)
    plot_multi_hist(data_sets, 10,
                    labels, xlabel='', ylabel='log frequency',
                    title=cell_types[ctype])

    plt.savefig(samples_path+stimul_type+'-'+cell_types[ctype]+'-hist.png')

# plt.show()
print('hej')