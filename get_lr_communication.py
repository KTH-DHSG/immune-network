import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy
from copy import copy
from pickle import load, dump
import os
from dataset_lib import *

spec_type = 'mean' # mean or sum

ligand_recep_db = pd.read_csv('../NATMI/lrdbs/lrc2p.csv')
samples_path = '../covid_separated/'
sample_names = [s for s in os.listdir(samples_path) if 'h5ad' in s]

# for sample_name in sample_names:
#     sample = scanpy.read_h5ad(samples_path+sample_name)
#     gene_library = sample.to_df().columns
#     cell_types, cells_dict = separate_by_celltype(sample)
#     mean_expression_dict = dict()
#     sum_expression_dict = dict()
#     for c in cells_dict.keys():
#         mean_expression_dict[c] = cells_dict[c].to_df().mean()
#         sum_expression_dict[c] = cells_dict[c].to_df().sum()
#     cell_types_n = len(cell_types)
#     pairs = []
#     communication_spec_matrices = np.zeros((cell_types_n, cell_types_n, 0))
#     for id in ligand_recep_db.index:
#         ligand = ligand_recep_db.iloc[id]['Ligand gene symbol']
#         receptor = ligand_recep_db.iloc[id]['Receptor gene symbol']
#         if not all_genes_present([ligand, receptor], gene_library):
#             continue
#         ligand_expression = np.zeros((cell_types_n, 1))
#         receptor_expression = np.zeros((cell_types_n, 1))
#         for id in range(len(cell_types)):
#             if spec_type == 'mean':
#                 ligand_expression[id] = mean_expression_dict[cell_types[id]][ligand]
#                 receptor_expression[id] = mean_expression_dict[cell_types[id]][receptor]
#             elif spec_type == 'sum':
#                 ligand_expression[id] = sum_expression_dict[cell_types[id]][ligand]
#                 receptor_expression[id] = sum_expression_dict[cell_types[id]][receptor]
#             else:
#                 print('unkown spec_type ' + spec_type)
#
#         communication_sum = (ligand_expression@receptor_expression.T).reshape(cell_types_n, cell_types_n, 1)
#         communication_spec = communication_sum / np.sum(communication_sum)
#         communication_spec_matrices = np.concatenate((communication_spec_matrices, communication_spec), axis=2)
#         pairs.append([ligand, receptor])
#     '''saving communication matrices'''
#     communication_dict = {}
#     communication_dict['matrices'] = communication_spec_matrices
#     communication_dict['pairs'] = pairs
#     communication_dict['cell_types'] = cell_types
#
#     with open('covid_data/' + sample_name + 'specificiyty_' + spec_type + '_matrices.pickle', 'wb') as handle:
#         dump(communication_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


'''Find LR pairs associated with salient genes'''
salient_genes_file = '../120Samples/OriginalData/3hour_response_genes.pickle'

with open(salient_genes_file, 'rb') as handle:
    salient_genes = load(handle)

salient_genes = salient_genes['up']
b = np.where([f in salient_genes for f in (ligand_recep_db['Ligand gene symbol'])])
print(ligand_recep_db.loc[b])
b = np.where([f in salient_genes for f in (ligand_recep_db['Receptor gene symbol'])])
print(ligand_recep_db.loc[b])
print('hi')