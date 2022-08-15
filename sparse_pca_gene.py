import numpy as np
from sklearn.decomposition import SparsePCA, PCA
import matplotlib.pyplot as plt
import pickle
import os

folder = '../covid_data_summary/'
files = [f for f in os.listdir(folder) if 'pickle' in f]
with open(folder + "gene_list.pickle", "rb") as f:
    gene_list = pickle.load(f)

cell_types_of_interest = ['NK']
expression_type = 'sum' # mean or sum
cell_specific_points = True
days = []
sos_vectors = None
patient = ''
patients = []
tot_expressions = []
mean_expressions = []
cell_specific_expressions = []
cell_types = []
for f in files:
    if 'patient' + patient not in f:
        continue
    crnt_patient = f[f.find('patient')+len('patient')]
    patients.append(crnt_patient)
    with open(folder+f, 'rb') as handle:
        summary_dict = pickle.load(handle)
    tot_expression = None
    tot_count = 0
    for ct in summary_dict.keys():
        if ct not in cell_types_of_interest:
            continue
        if tot_expression is None:
            tot_expression = np.zeros(summary_dict[ct]['gene_expression'].shape)
            mean_expression = np.zeros(summary_dict[ct]['gene_expression'].shape)
        cell_types.append(ct)
        tot_expression += summary_dict[ct]['gene_expression'].to_numpy()*summary_dict[ct]['cell_count']
        tot_count += summary_dict[ct]['cell_count']
        cell_specific_expressions.append(summary_dict[ct]['gene_expression'].to_numpy())
    tot_expressions.append(tot_expression)
    mean_expressions.append(tot_expression/tot_count)


if cell_specific_points:
    data = np.array(cell_specific_expressions)
else:
    tot_expressions = np.array(tot_expressions)
    mean_expressions = np.array(mean_expressions)
    if expression_type == 'mean':
        data = mean_expressions
    elif expression_type == 'sum':
        data = tot_expressions

n_components = 3
transformer = SparsePCA(n_components=n_components, random_state=0)
transformer.fit(data)
# pickle.dump(transformer, open("gene_expression_pca_patients.pkl", "wb"))
# with open("gene_expression_pca_patients.pkl", "rb") as f:
#     transformer = pickle.load(f)

print("The significant genes are:")
for i in range(n_components):
    print("LV", i, ": ")
    genes_in_lv = np.where(abs(transformer.components_[i, :])>0.2*max(abs(transformer.components_[i, :])))
    print([gene_list[i] for i in genes_in_lv])
    print([transformer.components_[i, j] for j in genes_in_lv])



data_transformed = transformer.transform(data)
data_recovered = np.dot(data_transformed, transformer.components_) + transformer.mean_
print('mean error naive: ', np.linalg.norm(data - transformer.mean_))
print('mean error PCA: ', np.linalg.norm(data - data_recovered))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

if cell_specific_points:
    cell_types_set = set(cell_types)
    plots = []
    plot_names = []
    for ct in cell_types_set:
        mask = np.array(cell_types) == ct
        crnt_plot = ax.scatter(data_transformed[mask, 0], data_transformed[mask, 1], data_transformed[mask, 2])
        plots.append(crnt_plot)
        plot_names.append(ct)
    ax.legend(plots, plot_names)

else:
    patiens_set = sorted(set(patients))

    plots = []
    plot_names = []
    plotted_set = sorted([])

    # Healthy persons
    mask = np.zeros((len(patients), ))
    for p in patiens_set:
        if p not in ['K', 'L', 'M']:
            continue
        mask = np.logical_or(mask, np.array(patients) == p)

    healthy_plot = ax.scatter(data_transformed[mask, 0], data_transformed[mask, 1],
                           data_transformed[mask, 2])
    plots.append(healthy_plot)
    plot_names.append('Healthy persons')

    # Affected persons
    mask = np.zeros((len(patients), ))
    for p in patiens_set:
        if p in ['K', 'L', 'M']:
            continue
        mask = np.array(patients) == p
        affected_plot, = ax.plot(data_transformed[mask, 0], data_transformed[mask, 1],
                                   data_transformed[mask, 2], '.-')
        plots.append(affected_plot)
        plot_names.append('Patient ' + p)
    ax.legend(plots, plot_names)

ax.set_xlabel('LV0')
ax.set_ylabel('LV1')
ax.set_zlabel('LV2')

plt.show()
