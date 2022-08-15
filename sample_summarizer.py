'''Turn h5ad samples to small dictionaries where each cell_type is mapped to it's count and average gene expression'''
import pandas as pd
import scanpy
from copy import copy
import pickle
import os

def all_genes_present(check_list, gene_library):
    return all([gene in gene_library for gene in check_list])


def get_cell_type(dataset, cell_type):
    '''returns only the cells of cell_type from dataset'''
    mask = dataset.obs['cell_type'] == cell_type
    return dataset[mask]


def separate_by_celltype(dataset):
    '''returns cell types present in dataset and a dictionary mapping each cell_type name to the dataset of cells in
    dataset with that type'''
    cell_types = list(set([d.obs['cell_type'][0] for d in dataset]))
    cell_types.sort()
    cells_dict = dict()
    for ct in cell_types:
        specific_cells = get_cell_type(dataset, ct)
        cells_dict[ct] = copy(specific_cells)
    return cell_types, cells_dict


samples_path = '../covid_separated/'
sample_names = [s for s in os.listdir(samples_path) if 'h5ad' in s]

for sample_name in sample_names:
    sample = scanpy.read_h5ad(samples_path+sample_name)
    gene_library = sample.to_df().columns
    cell_types, cells_dict = separate_by_celltype(sample)
    sample_summary = {}
    for c in cells_dict.keys():
        df = cells_dict[c].to_df()
        sample_summary[c] = {}
        sample_summary[c]['gene_expression'] = df.mean()
        sample_summary[c]['cell_count'] = df.shape[0]

    with open('covid_data_summary/' + sample_name + '_summary.pickle', 'wb') as handle:
        pickle.dump(sample_summary, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('../covid_data_summary/gene_list.pickle', 'wb') as handle:
    pickle.dump(sample.var_names, handle, protocol=pickle.HIGHEST_PROTOCOL)