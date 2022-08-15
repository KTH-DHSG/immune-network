import numpy as np
from sklearn.decomposition import SparsePCA, PCA
import pickle
import matplotlib.pyplot as plt

'''getting a patient's longitudinal data'''
patient = 'A'

import os

folder = '../covid_data/'
files = [f for f in os.listdir(folder) if 'pickle' in f]

spec_type = 'mean' # mean or sum
days = []
sos_vectors = None
for f in files:
    if 'patient' + patient not in f or spec_type not in f:
        continue
    day = int(f[f.find('day')+3:f.find('day')+5])
    with open(folder+f, 'rb') as handle:
        specificity_dict = pickle.load(handle)
    if sos_vectors is None:
        sos_vectors = np.zeros((0, len(specificity_dict['cell_types'])**2))
    specificity_matrices = specificity_dict['matrices']
    sum_of_specificity = np.nansum(specificity_matrices, axis=2)

    sos_vectors = np.concatenate((sos_vectors, sum_of_specificity.reshape((1, -1))), axis=0)
    days.append(day)

days = np.array(days)
transformer = SparsePCA(n_components=3, random_state=0)
# transformer = PCA(n_components=3, random_state=0)
transformer.fit(sos_vectors)

sos_vectors_transformed = transformer.transform(sos_vectors)

sos_vectors_recovered = np.dot(sos_vectors_transformed, transformer.components_) + transformer.mean_

day_agrsort = np.argsort(days)
plt.plot(days[day_agrsort], sos_vectors_transformed[day_agrsort, :])
plt.show()
# print(sos_vectors_transformed.shape)



