import numpy
import scanpy
import numpy as np
from pickle import dump, load
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import defaultdict
from plotting_tools import plot_multi_hist


samples_path = '../120Samples/OriginalData/cell_type_specific_expressions/full/'
target_path = '..   /120Samples/OriginalData/cell_type_specific_expressions/3h genes/'
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
        plt.gcf().set_size_inches(50, 20.0)
        plt.savefig(target_path + 'salient_genes_' + ct + ' at ' + tg + '.png', dpi=300)
        # plt.show()
        plt.close()
