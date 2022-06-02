import numpy as np
from  sklearn.decomposition import SparsePCA, PCA
import matplotlib.pyplot as plt
import pickle
import os

spec_type = 'mean' # mean or sum

folder = 'covid_data_summary/'
files = [f for f in os.listdir(folder) if 'pickle' in f]

spec_type = 'mean'
days = []
sos_vectors = None
patient = ''
patients = []
tot_expressions = []
for f in files:
    if 'patient' + patient not in f:
        continue
    crnt_patient = f[f.find('patient')+7]
    patients.append(crnt_patient)
    with open(folder+f, 'rb') as handle:
        summary_dict = pickle.load(handle)
    tot_expression = None
    for ct in summary_dict.keys():
        if tot_expression is None:
            tot_expression = np.zeros((summary_dict[ct]['gene_expression'].shape))
        tot_expression += summary_dict[ct]['gene_expression'].to_numpy()*summary_dict[ct]['cell_count']
    tot_expressions.append(tot_expression)

tot_expressions = np.array(tot_expressions)

n_components = 3
transformer = SparsePCA(n_components=n_components, random_state=0)
# transformer = PCA(n_components=n_components, random_state=0)
transformer.fit(tot_expressions)
tot_expressions_transformed = transformer.transform(tot_expressions)
tot_expressions_recovered = np.dot(tot_expressions_transformed, transformer.components_) + transformer.mean_

print('mean error naive: ', np.linalg.norm(tot_expressions-transformer.mean_))

print('mean error PCA: ', np.linalg.norm(tot_expressions-tot_expressions_recovered))

patiens_set = sorted(set(patients))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_list = []
plotted_set = sorted([])
for p in patiens_set:
    if p not in ['A', 'K', 'L', 'M']:
        continue
    plotted_set.append(p)
    mask = np.array(patients) == p
    plot_name = ax.scatter(tot_expressions_transformed[mask, 0], tot_expressions_transformed[mask, 1],
                           tot_expressions_transformed[mask, 2], )
    plot_list.append(plot_name)
ax.legend(plot_list, ['Patient ' + p for p in plotted_set])
plt.show()
