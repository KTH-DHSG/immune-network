import numpy
import scanpy
import numpy as np
from pickle import dump, load
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import defaultdict


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

samples_path = '../120Samples/OriginalData/cell_type_specific_expressions/full/'
target_path = '../120Samples/OriginalData/cell_type_specific_expressions/3h genes/'
salient_genes_file = '../120Samples/OriginalData/3hour_response_genes.pickle'

with open(salient_genes_file, 'rb') as handle:
    salient_genes = load(handle)

# salient_genes = salient_genes['up'].union(salient_genes['down'])
salient_genes = salient_genes['up']
# print(salient_genes)
with open(samples_path+'obs_info.pickle', 'rb') as handle:
    obs_info = load(handle)

donors = obs_info['subject']
cell_types = obs_info['cell_types']
timepoints = obs_info['timepoints']
genes = np.array(obs_info['genes'])
salient_genes = list(salient_genes)

with open(samples_path+'cell_frequencies.pickle', 'rb') as handle:
    cell_frequencies = load(handle)

# salient_genes_mask = np.full((len(genes,)), False)
# for s in salient_genes:
#     salient_genes_mask[np.where(genes == s)[0]] = True
# for t in timepoints:
#     big_dataframe = pd.DataFrame(cell_frequencies[t], columns=cell_types).fillna(0)
#     for ct in cell_types:
#         labels = [gene_name + '_' + ct for gene_name in salient_genes]
#         with open(samples_path + ct + '_' + t + '.pickle', 'rb') as handle:
#             data = load(handle)[:, salient_genes_mask]
#         data_frame = pd.DataFrame(data, columns=labels).fillna(0)
#         big_dataframe = pd.concat((big_dataframe, data_frame), axis=1)
#         # plot_multi_hist([np.log10(big_dataframe[g + '_' + ct]+0.00001) for g in salient_genes], 20,
#         #                 labels=[g for g in salient_genes], hist_range=[-5, 0], title=ct+' at '+t)
#         # plt.savefig(target_path + 'salient_genes_' + ct+' at '+t + '.png')
#         # plt.show()
#     with open(target_path+'salient_expressions_' + t + '.pickle', 'wb') as handle:
#         dump(big_dataframe, handle)


'''Group different responses together'''

timepoint_groups = ['UT', '3h', '24h']

tot_weight = dict()
tot_express = dict()
for tg in timepoint_groups:
    tot_weight[tg] = defaultdict(lambda: np.zeros((len(donors,))))
    tot_express[tg] = defaultdict(lambda: np.zeros((len(donors), len(salient_genes))))
    for t in timepoints:
        if tg not in t:
            continue
        with open(target_path + 'salient_expressions_' + t + '.pickle', 'rb') as handle:
            df = load(handle)
        for ct in cell_types:
            crnt_weight = df[ct]
            crnt_express = df[[gene_name + '_' + ct for gene_name in salient_genes]]
            tot_express[tg][ct] = (np.diag((tot_weight[tg][ct] + crnt_weight) ** -1)@(np.diag(tot_weight[tg][ct])@tot_express[tg][ct] + np.diag(crnt_weight)@crnt_express)).fillna(0)
            tot_weight[tg][ct] += crnt_weight

for tg in timepoint_groups:
    for ct in cell_types:
        labels = [gene_name + '_' + ct for gene_name in salient_genes]
        plot_multi_hist([np.log10(tot_express[tg][ct][g + '_' + ct]+0.00001) for g in salient_genes], 20,
                        labels=[g for g in salient_genes], hist_range=[-5, 0], title=ct + ' at ' + tg)
        plt.gcf().set_size_inches(22.5, 10.5)
        plt.savefig(target_path + 'salient_genes_' + ct + ' at ' + tg + '.png')
        # plt.show()
        plt.close()
