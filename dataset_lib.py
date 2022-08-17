from copy import copy

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