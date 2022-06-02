import numpy as np
from  sklearn.decomposition import SparsePCA, PCA
import matplotlib.pyplot as plt
import pickle
import os

folder = 'covid_data_summary/'
files = [f for f in os.listdir(folder) if 'pickle' in f]

expression_type = 'sum' # mean or sum
days = []
sos_vectors = None
patient = ''
patients = []
tot_expressions = []
mean_expressions = []
for f in files:
    if 'patient' + patient not in f:
        continue
    crnt_patient = f[f.find('patient')+7]
    patients.append(crnt_patient)
    with open(folder+f, 'rb') as handle:
        summary_dict = pickle.load(handle)
    tot_expression = None
    tot_count = 0
    for ct in summary_dict.keys():
        if tot_expression is None:
            tot_expression = np.zeros(summary_dict[ct]['gene_expression'].shape)
            mean_expression = np.zeros(summary_dict[ct]['gene_expression'].shape)

        tot_expression += summary_dict[ct]['gene_expression'].to_numpy()*summary_dict[ct]['cell_count']
        tot_count += summary_dict[ct]['cell_count']
    tot_expressions.append(tot_expression)
    mean_expressions.append(tot_expression/tot_count)

tot_expressions = np.array(tot_expressions)
mean_expressions = np.array(mean_expressions)

if expression_type == 'mean':
    data = mean_expressions
elif expression_type == 'sum':
    data = tot_expressions

n_components = 3
transformer = SparsePCA(n_components=n_components, random_state=0)
# transformer = PCA(n_components=n_components, random_state=0)
transformer.fit(data)
data_transformed = transformer.transform(data)
data_recovered = np.dot(data_transformed, transformer.components_) + transformer.mean_

print('mean error naive: ', np.linalg.norm(data-transformer.mean_))

print('mean error PCA: ', np.linalg.norm(data-data_recovered))

patiens_set = sorted(set(patients))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_list = []
plotted_set = sorted([])

# Healthy persons
mask = np.zeros((len(patients), ))
for p in patiens_set:
    if p not in ['K', 'L', 'M']:
        continue
    mask = np.logical_or(mask, np.array(patients) == p)

healthy_plot_name = ax.scatter(data_transformed[mask, 0], data_transformed[mask, 1],
                       data_transformed[mask, 2])

# Affected persons
mask = np.zeros((len(patients), ))
for p in patiens_set:
    if p in ['K', 'L', 'M']:
        continue
    mask = np.logical_or(mask, np.array(patients) == p)


affected_plot_name = ax.scatter(data_transformed[mask, 0], data_transformed[mask, 1],
                           data_transformed[mask, 2])
ax.legend([healthy_plot_name, affected_plot_name], ['Healthy persons', 'Affected persons'])
plt.savefig('gene_expression_pca.png')
plt.show()
