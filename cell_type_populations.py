import scanpy
import numpy as np
from pickle import dump, load
import matplotlib.pyplot as plt
from plotting_tools import plot_multi_hist
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