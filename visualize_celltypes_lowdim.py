from typing import BinaryIO

import numpy as np
from sklearn.decomposition import SparsePCA, PCA
import matplotlib.pyplot as plt
import pickle
import os
from dataset_lib import separate_by_celltype
import scanpy
import time

samples_folder = '../covid_separated/'
sample_files = [f for f in os.listdir(samples_folder) if 'h5ad' in f]
with open(samples_folder + "gene_list.pickle", "rb") as f:
    gene_list = pickle.load(f)

pca_folder = '../covid_celltype_pca/'
reference_patient = 'M'
cell_type = 'Monocytes'
ctype_output_folder = pca_folder + 'patient_' + reference_patient + '/day_0/' + cell_type + '/'
with open(ctype_output_folder + "gene_expression_pca.pkl", "rb") as pca_file:
    transformer = pickle.load(pca_file)

plots = []
plot_names = []
fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
ax = fig.add_subplot()
for patient in ['A', 'L']:
    for f in sample_files:
        if 'patient' + patient not in f:
            continue
        # patient = f[f.find('patient') + len('patient')]
        sample = scanpy.read_h5ad(samples_folder + f)
        day = int(sample.obs['days'][0])
        cell_types, cells_dict = separate_by_celltype(sample)
        gene_expressions = cells_dict[cell_type].to_df().to_numpy()
        gene_expressions_lowdim = transformer.transform(gene_expressions)
        patient_plot = ax.scatter(gene_expressions_lowdim[:, 0], gene_expressions_lowdim[:, 1])
        plots.append(patient_plot)
        plot_names.append('Patient_' + patient + ' Day_'+str(day))
        ax.legend(plots, plot_names)

plt.show()


