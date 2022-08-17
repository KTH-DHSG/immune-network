import numpy as np
from sklearn.decomposition import SparsePCA, PCA
import matplotlib.pyplot as plt
import pickle
import os
from dataset_lib import separate_by_celltype
import scanpy
import time

folder = '../covid_separated/'
output_folder = '../covid_celltype_pca'
files = [f for f in os.listdir(folder) if 'h5ad' in f]
with open(folder + "gene_list.pickle", "rb") as f:
    gene_list = pickle.load(f)
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

ncells_for_pca = 50
ncomponents_for_pca = 3

# Generating low-dim representation through PCA for each cell type in each healthy individual
for patient in ['L', 'M', 'K']:
    for f in files:
        if 'patient' + patient not in f:
            continue
        patient_output_folder = output_folder+'/patient_'+patient
        if not os.path.exists(patient_output_folder):
            os.mkdir(patient_output_folder)
        sample = scanpy.read_h5ad(folder + f)
        gene_library = sample.to_df().columns
        day = int(sample.obs['days'][0])
        day_output_folder = patient_output_folder + '/day_' + str(day)
        if not os.path.exists(day_output_folder):
            os.mkdir(day_output_folder)
        cell_types, cells_dict = separate_by_celltype(sample)
        for ct in cell_types:
            crnt_cells = cells_dict[ct].to_df()
            if crnt_cells.shape[0] < ncells_for_pca:
                continue
            ctype_output_folder = day_output_folder + '/' + ct
            if not os.path.exists(ctype_output_folder):
                os.mkdir(ctype_output_folder)
            sample_cells = crnt_cells.sample(n=ncells_for_pca).to_numpy()
            time1 = time.time()
            transformer = SparsePCA(n_components=ncomponents_for_pca, random_state=0)
            transformer.fit(sample_cells)
            print("PCA took ", time.time() - time1, "seconds")
            with open(ctype_output_folder+"/gene_expression_pca.pkl", "wb") as pca_file:
                pickle.dump(transformer, pca_file)


