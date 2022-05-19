import pickle
import matplotlib.pyplot as plt
import numpy as np

with open('specificiyty_matrices.pickle', 'rb') as handle:
    specificity_dict = pickle.load(handle)

matrices = specificity_dict['matrices']
genes = specificity_dict['matrices']
cell_types = specificity_dict['cell_types']

sum_of_specificity = np.nansum(matrices, axis=2)
normalized_sum = sum_of_specificity / np.max(sum_of_specificity)

fig, ax = plt.subplots()
im = ax.imshow(sum_of_specificity)

ax.set_xticks(np.arange(len(cell_types)), labels=cell_types)
ax.set_yticks(np.arange(len(cell_types)), labels=cell_types)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(cell_types)):
    for j in range(len(cell_types)):
        text = ax.text(j, i, "{:10.1f}".format(sum_of_specificity[i, j]),
                       ha="center", va="center", color="w")

ax.set_title("Sum of specificity in ligand-receptor communication")
# fig.tight_layout()
plt.show()