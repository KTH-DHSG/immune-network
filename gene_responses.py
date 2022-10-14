import scanpy
import numpy as np
from pickle import dump, load
import matplotlib.pyplot as plt

samples_path = '../120Samples/OriginalData/'
db_name = 'adata_all.h5ad' # obs = ['celltype_sub', 'celltype', 'subject', 'timepoint', '10Xchemistry']

'''Getting and writing significant average expression changes'''
# data = scanpy.read_h5ad(samples_path+db_name)
#
# donors = sorted(set(data.obs['subject']))
# cell_types = sorted(set(data.obs['celltype']))
# timepoints = set(data.obs['timepoint'])
# genes = list(data.var_names)
#
# obs_info = {'subject': donors, 'cell_types': cell_types, 'timepoints': timepoints, 'genes':genes}
# with open(samples_path+'obs_info.pickle', 'wb') as handle:
#     dump(obs_info, handle)

# for t in timepoints:
#     expression_matrix = np.zeros((len(donors), len(genes)))
#     for d in range(len(donors)):
#         mask = (data.obs['timepoint'] == t) & (data.obs['subject'] == donors[d])
#         if np.sum(mask) < 100:
#             continue
#         mean_expression = np.mean(data[mask].X, axis=0)
#         expression_matrix[d, :] = mean_expression
#     with open(samples_path + 'mean_gene_expression_' + t + '.pickle', 'wb') as handle:
#         dump(expression_matrix, handle)

'''Loading short info and significant average expression changes'''


def plot_multi_hist(data_sets, number_of_bins, labels, xlabel="Data sets", ylabel="Data values", title='', log_pop=False):
    # Computed quantities to aid plotting
    plt.figure()
    hist_range = (np.min(data_sets), np.max(data_sets))
    binned_data_sets = [
        np.histogram(d, range=hist_range, bins=number_of_bins)[0]
        for d in data_sets
    ]
    if log_pop:
        binned_data_sets = np.log(binned_data_sets)
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


with open(samples_path+'obs_info.pickle', 'rb') as handle:
    obs_info = load(handle)

donors = obs_info['subject']
cell_types = obs_info['cell_types']
timepoints = obs_info['timepoints']
genes = obs_info['genes']

gene_expressions = dict()

for t in timepoints:
    with open(samples_path + 'mean_gene_expression_' + t + '.pickle', 'rb') as handle:
        gene_expressions[t] = load(handle)
    row_max = np.max(gene_expressions[t], axis=1)
    row_max[row_max == 0] = 1
    gene_expressions[t] /= row_max[:, np.newaxis]

untreated_expression = gene_expressions['UT']

threshholds = list(np.arange(0.5, 5.0, 0.2))
regs = list(np.logspace(-1, -7, 13))
n_select_genes = 10
target_timepoints = '24'
gene_file_name = target_timepoints + 'hour_response_genes.pickle'
tot_genes_up = set()
tot_genes_down = set()
for reg in regs:
    for threshhold in threshholds:
        genes = np.array(genes)
        genes_up = []
        genes_down = []
        for st in timepoints:
            if target_timepoints not in st:
                continue
            stimulated_expression = gene_expressions[st]
            diff = np.log2(stimulated_expression + reg) - np.log2(untreated_expression + reg)
            significant_upreg = diff > threshhold
            significant_downreg = -diff > threshhold
            significant_diff = significant_upreg | significant_downreg
            gene_votes = np.sum(significant_diff, axis=0)
            gene_votes_up = np.sum(significant_upreg, axis=0)
            gene_votes_down = np.sum(significant_downreg, axis=0)
            genes_up += list(genes[np.where(gene_votes_up > 100)[0]])
            genes_down += list(genes[np.where(gene_votes_down > 100)[0]])
        gene_count = len(set(genes_up + genes_down))
        if gene_count < n_select_genes:
            tot_genes_up = tot_genes_up.union(genes_up)
            tot_genes_down = tot_genes_down.union(genes_down)
            break

tot_genes = {'up': tot_genes_up, 'down': tot_genes_down}
with open(samples_path + gene_file_name, 'wb') as handle:
    dump(tot_genes, handle)
print(tot_genes_up)
print(tot_genes_down)

    # print(st)
    # print(genes[np.where(gene_votes > 100)[0]])
    # print(genes[np.where(gene_votes_up > 100)[0]])
    # print(genes[np.where(gene_votes_down > 100)[0]])
    # plot_multi_hist([gene_votes, gene_votes_up, gene_votes_down], 10,
    #                 ['total votes', 'up regulated', 'down regulated'], xlabel='', ylabel='log10 of votes',
    #                 title='Gene expression change votes in ' + st + ' response', log_pop=True)
    # plt.savefig(samples_path+'gene_vote'+st+'.png')
    # plt.hist(gene_votes, range=[10, 120])
    # plt.show()


print('hej')
