import numpy as np
import pickle
import os
import xlsxwriter


pca_folder = '../covid_celltype_pca/'
with open(pca_folder + "gene_list.pickle", "rb") as f:
    gene_list = pickle.load(f)
for patient_folder in os.scandir(pca_folder):
    if not os.path.isdir(patient_folder):
        continue
    for day_folder in os.scandir(patient_folder):
        if not os.path.isdir(day_folder):
            continue
        workbook = xlsxwriter.Workbook(day_folder.path + '/pca_genes.xlsx')
        for cell_type_folder in os.scandir(day_folder):
            if not os.path.isdir(cell_type_folder):
                continue
            worksheet = workbook.add_worksheet(cell_type_folder.name)
            with open(cell_type_folder.path + "/gene_expression_pca.pkl", "rb") as pca_file:
                transformer = pickle.load(pca_file)
            n_components = transformer.components_.shape[0]
            row = 0
            for i in range(n_components):
                worksheet.write(row, 0, 'LV'+str(i))
                row+=1
                genes_in_lv = np.where(
                    abs(transformer.components_[i, :]) > 0.2 * max(abs(transformer.components_[i, :])))
                for j in genes_in_lv[0]:
                    worksheet.write(row, 0, gene_list[j])
                    worksheet.write(row, 1, transformer.components_[i, j])
                    row += 1
        workbook.close()


# reference_patient = 'M'
# cell_type = 'Monocytes'
# ctype_output_folder = pca_folder + 'patient_' + reference_patient + '/day_0/' + cell_type + '/'
# with open(ctype_output_folder + "gene_expression_pca.pkl", "rb") as pca_file:
#     transformer = pickle.load(pca_file)
#
# plots = []
# plot_names = []
# fig = plt.figure()
# # ax = fig.add_subplot(projection='3d')
# ax = fig.add_subplot()
# for patient in ['A', 'L']:
#     for f in sample_files:
#         if 'patient' + patient not in f:
#             continue
#         # patient = f[f.find('patient') + len('patient')]
#         sample = scanpy.read_h5ad(samples_folder + f)
#         day = int(sample.obs['days'][0])
#         cell_types, cells_dict = separate_by_celltype(sample)
#         gene_expressions = cells_dict[cell_type].to_df().to_numpy()
#         gene_expressions_lowdim = transformer.transform(gene_expressions)
#         patient_plot = ax.scatter(gene_expressions_lowdim[:, 0], gene_expressions_lowdim[:, 1])
#         plots.append(patient_plot)
#         plot_names.append('Patient_' + patient + ' Day_'+str(day))
#         ax.legend(plots, plot_names)
#
# plt.show()
#
#
