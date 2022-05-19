import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy
from copy import copy
import pickle


def all_genes_present(check_list, gene_library):
    return all([gene in gene_library for gene in check_list])


def get_cell_type(dataset, cell_type):
    '''returns only the cells of cell_type from dataset'''
    mask = dataset.obs['ct_cov'] == cell_type
    return dataset[mask]


def separate_by_celltype(dataset):
    '''returns cell types present in dataset and a dictionary mapping each cell_type name to the dataset of cells in
    dataset with that type'''
    cell_types = list(set([d.obs['ct_cov'][0] for d in dataset]))
    cells_dict = dict()
    for ct in cell_types:
        specific_cells = get_cell_type(dataset, ct)
        cells_dict[ct] = copy(specific_cells)
    return cell_types, cells_dict


ligand_recep_db = pd.read_csv('../NATMI/lrdbs/lrc2p.csv')
samples_path = '../longitudinal_samples_healthy/'
sample_name = 'IGTB1290_batch0.h5ad'
sample = scanpy.read_h5ad(samples_path+sample_name)
gene_library = sample.to_df().columns
cell_types, cells_dict = separate_by_celltype(sample)
cell_types_n = len(cell_types)
pairs = []
communication_spec_matrices = np.zeros((cell_types_n, cell_types_n, 0))
for id in ligand_recep_db.index:
    ligand = ligand_recep_db.iloc[id]['Ligand gene symbol']
    receptor = ligand_recep_db.iloc[id]['Receptor gene symbol']
    if not all_genes_present([ligand, receptor], gene_library):
        continue
    ligand_expression = np.zeros((cell_types_n, 1))
    receptor_expression = np.zeros((cell_types_n, 1))
    for id in range(len(cell_types)):
        ligand_expression[id] = cells_dict[cell_types[id]].to_df()[ligand].mean()
        receptor_expression[id] = cells_dict[cell_types[id]].to_df()[receptor].mean()
    communication_sum = (ligand_expression@receptor_expression.T).reshape(cell_types_n, cell_types_n, 1)
    communication_spec = communication_sum / np.sum(communication_sum)
    communication_spec_matrices = np.concatenate((communication_spec_matrices, communication_spec), axis=2)
    pairs.append([ligand, receptor])
'''saving communication matrices'''
communication_dict = {}
communication_dict['matrices'] = communication_spec_matrices
communication_dict['pairs'] = pairs
communication_dict['cell_types'] = cell_types

with open('main/specificiyty_matrices.pickle', 'wb') as handle:
    pickle.dump(communication_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)