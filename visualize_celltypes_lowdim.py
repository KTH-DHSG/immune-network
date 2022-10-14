import matplotlib.pyplot as plt
import pickle
import os
from dataset_lib import separate_by_celltype
import scanpy
import numpy as np

samples_folder = '../covid_separated/'
sample_files = [f for f in os.listdir(samples_folder) if 'h5ad' in f]
with open(samples_folder + "gene_list.pickle", "rb") as f:
    gene_list = pickle.load(f)

pca_folder = '../covid_celltype_pca/'
reference_patient = 'L'
patient_folder = pca_folder + 'patient_' + reference_patient + '/day_0/'
cell_types = [x[1] for x in os.walk(patient_folder)][0]
# cell_exclude_list = ['Innate T', 'NK', 'Monocytes', 'Activated CD4+ T', 'Effector Memory CD8+', 'Effector T',
#                      'Memory CD8+ T', 'Naive CD4+ T']
# cell_types = [c for c in cell_types if c not in cell_exclude_list]
# cell_types = ['Innate T']
plot_vars = [0, 2]
n_plot_points = 100
healthy_list = ['K', 'L', 'M']
include_list = ['K', 'L', 'F', 'N']
for cell_type in cell_types:
    with open(patient_folder + cell_type + "/gene_expression_pca.pkl", "rb") as pca_file:
        transformer = pickle.load(pca_file)
    n_features = transformer.components_.shape[0]
    healthy_data = np.zeros((0, n_features))
    patient_data = np.zeros((0, n_features))
    plots = []
    plot_names = []
    fig = plt.figure()
    ax = fig.add_subplot()
    for f in sample_files:
        if 'patient' not in f:
            continue
        patient = f[f.find('patient') + len('patient')]
        if patient not in include_list:
            continue
        sample = scanpy.read_h5ad(samples_folder + f)
        cell_types, cells_dict = separate_by_celltype(sample)
        gene_expressions = cells_dict[cell_type].to_df().to_numpy()
        gene_expressions_lowdim = transformer.transform(gene_expressions)
        if patient in healthy_list:
            healthy_data = np.concatenate((healthy_data, gene_expressions_lowdim), axis=0)
        else:
            patient_data = np.concatenate((patient_data, gene_expressions_lowdim), axis=0)
    healthy_plot = ax.scatter(healthy_data[:, plot_vars[0]], healthy_data[:, plot_vars[1]])
    patient_plot = ax.scatter(patient_data[:, plot_vars[0]], patient_data[:, plot_vars[1]])
    ax.legend([healthy_plot, patient_plot], ['Healthy individuals', 'covid patients'])
    ax.set_xlabel('LV' + str(plot_vars[0]))
    ax.set_ylabel('LV' + str(plot_vars[1]))
    ax.set_title(cell_type)
    plt.savefig('gene_expression_' + cell_type + '_LV' + str(plot_vars[0]) + str(plot_vars[1]) + '.png')
