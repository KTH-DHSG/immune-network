"""
The goal is to obtain a linear discrete-time X(t+T)=Ax(t) + c model of the immune response to a certain stimulant.
Given sc-RNA data from immune cells before and after the stimulation, we first need to obtain a representation of the
state X.
In this code we obtain the state through the following steps:
1- Get significant genes in terms of differential gene expression: i.e. Genes that are consistently up-down regulated in
 each cell type;
2- (optional) Add genes that are in a Ligand-Receptor (LR) relation with the identified genes;
3- In each cell type we project the log expression of selected genes to a lower dimensional space. (Using (Sparse) PCA);
4- Decompose each PCA axis to sub-intervals and count the number of genes in each interval;
5- The (normalized) number of cells in each interval becomes the state X, used for linear identification X(t+T)=Ax(t)+c.
"""
import anndata
import scanpy
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import defaultdict
from plotting_tools import plot_multi_hist
from sklearn.decomposition import SparsePCA, PCA
from sklearn import linear_model
import random


def get_mask_anndata(dataset: anndata.AnnData, criteria: dict):
    mask = np.full((dataset.shape[0],), True)
    for d in criteria.keys():
        if type(criteria[d]) is list:
            temp_mask = np.full((dataset.shape[0],), False)
            for crit in criteria[d]:
                temp_mask = temp_mask | (dataset.obs[d] == crit)
        else:
            temp_mask = dataset.obs[d] == criteria[d]
        mask = mask & temp_mask
    return mask


def get_freqs_anndata(dataset: anndata.AnnData, criteria: dict):
    """
    Get the relative count (frequency) of rows according to one criteria
    """
    assert len(criteria) == 1
    count = []
    for d in criteria.keys():
        assert type(criteria[d]) is list
        for c in criteria[d]:
            count.append(np.sum(get_mask_anndata(dataset, {d: c})))
    count = np.array(count)
    freq = count / np.maximum(np.sum(count), 1)
    return freq


def summarize_anndata(dataset: anndata.AnnData, criteria: dict):
    """
    Gets an AnnData, 1- filters the entries according to the criteria and 2- generates a single row dataframe with the
    average value of the columns in the filtered rows and 3- adds the criteria to the dataframe and returns it.
    """
    mask = get_mask_anndata(dataset, criteria)
    df = dataset[mask].to_df().mean().to_frame().transpose()
    info_count = 0
    for d in criteria.keys():
        if type(criteria[d]) is list:
            df.insert(info_count, d, criteria[d][0])
            df.at[0, d] = criteria[d]
        else:
            df.insert(info_count, d, criteria[d])
        info_count += 1
    df.insert(info_count, 'cell_count', np.sum(mask))
    info_count += 1
    df.insert(0, 'info_columns_count', info_count+1)
    return df


def get_salient_differential_expression(df1: pd.DataFrame, df2: pd.DataFrame, threshold=1.0, required_vote_rate=0.9,
                                        min_genes=5, max_genes=10):
    """
    Returns the genes where their log2 expression is changed by at least threshold in more than required_vote_rate of
    the rows.
    """
    assert df1.shape[0] == df2.shape[0]
    assert np.all(df1.columns == df2.columns)

    row_count = df1.shape[0]
    gene_count = df1.shape[1]
    reg = np.minimum(df1[df1 > 0].min(axis=0), df2[df2 > 0].min(axis=0)).fillna(1e-10)
    df_diff = np.log2(df2+reg) - np.log2(df1+reg)
    disqualified_rows = np.sum(df_diff.fillna(0) == 0, axis=1) == gene_count  # rows with all zeros
    disqualified_rows_count = np.sum(disqualified_rows)
    qualified_rows_count = row_count - disqualified_rows_count
    if disqualified_rows_count > qualified_rows_count:
        return df1.columns[[]]

    while True:
        if threshold > 0.0:
            df_pass = df_diff > threshold
        else:
            df_pass = df_diff < threshold

        votes = np.sum(df_pass, axis=0)
        required_votes = required_vote_rate * qualified_rows_count
        salient_genes_mask = votes > required_votes
        if np.sum(salient_genes_mask) < max_genes:
            break
        threshold += 0.1*np.sign(threshold)
    if np.sum(salient_genes_mask) < min_genes:
        salient_genes_mask = votes >= np.partition(votes, -min_genes)[-min_genes]

    return df1.columns[salient_genes_mask]


def get_linear_model(df1: pd.DataFrame, df2: pd.DataFrame, alpha=0.1, test_rate=0.1):
    '''
    Getting linear model df_2 = A df_1 + c
    '''
    assert df1.shape == df2.shape
    n_rows = df1.shape[0]
    n_test = int(test_rate*n_rows)
    test_rows = random.sample(range(n_rows), n_test)
    test_mask = np.full((n_rows,), False)
    test_mask[test_rows] = True
    train_mask = ~test_mask
    train_df1 = df1[train_mask]
    train_df2 = df2[train_mask]
    test_df1 = df1[test_mask]
    test_df2 = df2[test_mask]
    sparse_model = linear_model.Lasso(alpha=alpha)
    sparse_model.fit(train_df1, train_df2)
    train_var = (train_df2 - sparse_model.predict(train_df1)).var()
    test_var = (test_df2 - sparse_model.predict(test_df1)).var()
    return sparse_model, train_var, test_var




samples_path = '../120Samples/OriginalData/'
sample_name = 'adata_small_sample.h5ad'
# sample_name = 'adata_all.h5ad'
summary_file_name = 'gene_expression_summary.xlsx'
annotated_data = scanpy.read_h5ad(samples_path + sample_name)

'''Step 1: Getting significant genes'''
# time_point1 = 'UT'
# time_point2 = ['3hPA', '3hCA', '3hMTB']
# cell_types = set(annotated_data.obs['celltype'])
# subjects = set(annotated_data.obs['subject'])
#
# 'Getting average gene expressions in each sample over all cell types'
# data_summary_all_tp1 = pd.DataFrame()
# data_summary_all_tp2 = pd.DataFrame()
# for subject in subjects:
#     criteria = {'timepoint': time_point1, 'subject': subject}
#     temp = summarize_anndata(annotated_data, criteria)
#     data_summary_all_tp1 = pd.concat((data_summary_all_tp1, temp), ignore_index=True)
#
#     criteria = {'timepoint': time_point2, 'subject': subject}
#     temp = summarize_anndata(annotated_data, criteria)
#     data_summary_all_tp2 = pd.concat((data_summary_all_tp2, temp), ignore_index=True)
# gene_start = data_summary_all_tp2['info_columns_count'][0]
#
#
# 'Getting the up-regulated genes over all cell types'
# writer = pd.ExcelWriter(samples_path + summary_file_name, engine = 'xlsxwriter')
# upregulated_genes = get_salient_differential_expression(data_summary_all_tp1.iloc[:, gene_start:],
#                                                         data_summary_all_tp2.iloc[:, gene_start:], threshold=1.0)
#
# selected_summary1 = pd.concat((data_summary_all_tp1.iloc[:, :gene_start],
#                                data_summary_all_tp1[upregulated_genes]), axis=1)
# selected_summary2 = pd.concat((data_summary_all_tp2.iloc[:, :gene_start],
#                                data_summary_all_tp2[upregulated_genes]), axis=1)
# summary = pd.concat((selected_summary1, selected_summary2), axis=0)
# summary.to_excel(writer, sheet_name='up genes - all cells')
#
#
# 'Writing the down-regulated genes over all cell types'
# downregulated_genes = get_salient_differential_expression(data_summary_all_tp1.iloc[:, gene_start:],
#                                                           data_summary_all_tp2.iloc[:, gene_start:], threshold=-1.0)
#
# selected_summary1 = pd.concat((data_summary_all_tp1.iloc[:, :gene_start],
#                                data_summary_all_tp1[downregulated_genes]), axis=1)
# selected_summary2 = pd.concat((data_summary_all_tp2.iloc[:, :gene_start],
#                                data_summary_all_tp2[downregulated_genes]), axis=1)
# summary = pd.concat((selected_summary1, selected_summary2), axis=0)
# summary.to_excel(writer, sheet_name='down genes - all cells')
#
#
# 'Getting average gene expressions in each sample in individual cell types'
# data_summary_tp1 = pd.DataFrame()
# data_summary_tp2 = pd.DataFrame()
# for ct in cell_types:
#     for subject in subjects:
#         criteria = {'timepoint': time_point1, 'subject': subject, 'celltype': ct}
#         temp = summarize_anndata(annotated_data, criteria)
#         data_summary_tp1 = pd.concat((data_summary_tp1, temp), ignore_index=True)
#
#         criteria = {'timepoint': time_point2, 'subject': subject, 'celltype': ct}
#         temp = summarize_anndata(annotated_data, criteria)
#         data_summary_tp2 = pd.concat((data_summary_tp2, temp), ignore_index=True)
#
# 'Getting the up-regulated genes for individual cell types'
# for ct in cell_types:
#     gene_start = data_summary_tp2['info_columns_count'][0]
#     ct_rows = data_summary_tp1['celltype'] == ct
#     ct_upregulated_genes = get_salient_differential_expression(data_summary_tp1[ct_rows].iloc[:, gene_start:],
#                                                                data_summary_tp2[ct_rows].iloc[:, gene_start:],
#                                                                threshold=1.0)
#     upregulated_genes = np.concatenate((upregulated_genes, ct_upregulated_genes))
#     selected_summary1 = pd.concat((data_summary_tp1[ct_rows].iloc[:, :gene_start],
#                                    data_summary_tp1[ct_rows][ct_upregulated_genes]), axis=1)
#     selected_summary2 = pd.concat((data_summary_tp2[ct_rows].iloc[:, :gene_start],
#                                    data_summary_tp2[ct_rows][ct_upregulated_genes]), axis=1)
#     summary = pd.concat((selected_summary1, selected_summary2), axis=0)
#     summary.to_excel(writer, sheet_name='up genes - ' + ct)
#
# 'Getting the down-regulated genes for individual cell types'
# for ct in cell_types:
#     gene_start = data_summary_tp2['info_columns_count'][0]
#     ct_rows = data_summary_tp1['celltype'] == ct
#     ct_downregulated_genes = get_salient_differential_expression(data_summary_tp1[ct_rows].iloc[:, gene_start:],
#                                                                data_summary_tp2[ct_rows].iloc[:, gene_start:],
#                                                                threshold=-1.0)
#     downregulated_genes = np.concatenate((downregulated_genes, ct_downregulated_genes))
#     selected_summary1 = pd.concat((data_summary_tp1[ct_rows].iloc[:, :gene_start],
#                                    data_summary_tp1[ct_rows][ct_downregulated_genes]), axis=1)
#     selected_summary2 = pd.concat((data_summary_tp2[ct_rows].iloc[:, :gene_start],
#                                    data_summary_tp2[ct_rows][ct_downregulated_genes]), axis=1)
#     summary = pd.concat((selected_summary1, selected_summary2), axis=0)
#     summary.to_excel(writer, sheet_name='down genes - ' + ct)
#
# writer.save()
# significant_genes = np.concatenate((downregulated_genes, upregulated_genes))
# with open(samples_path + 'significant_genes.pickle', 'wb') as handle:
#     pickle.dump(significant_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)

'If the genes are already generated:'
with open(samples_path + 'significant_genes.pickle', 'rb') as handle:
    salient_genes = list(set(pickle.load(handle)))
print(len(set(salient_genes)))
'''Step 5: linear regression'''
time_point1 = 'UT'
time_point2 = ['3hMTB']
cell_types = list(set(annotated_data.obs['celltype']))
subjects = list(set(annotated_data.obs['subject']))

data_in = pd.DataFrame()
data_out = pd.DataFrame()

for s in subjects:
    mask = get_mask_anndata(annotated_data, {'subject': s})
    subject_data = annotated_data[mask]

    mask = get_mask_anndata(subject_data, {'timepoint': time_point1})
    subject_data_tp1 = subject_data[mask]

    mask = get_mask_anndata(subject_data, {'timepoint': time_point2})
    subject_data_tp2 = subject_data[mask]

    freqs_tp1 = get_freqs_anndata(subject_data_tp1, {'celltype': cell_types})
    freqs_tp2 = get_freqs_anndata(subject_data_tp2, {'celltype': cell_types})
    freqs_tp1 = pd.DataFrame(freqs_tp1.reshape((1, -1)), columns=cell_types)
    freqs_tp2 = pd.DataFrame(freqs_tp2.reshape((1, -1)), columns=cell_types)

    gene_exp_tp1 = subject_data_tp1.to_df()[salient_genes].mean().to_frame().transpose()
    gene_exp_tp2 = subject_data_tp2.to_df()[salient_genes].mean().to_frame().transpose()

    tot_data_tp1 = pd.concat((freqs_tp1, gene_exp_tp1), axis=1)
    tot_data_tp2 = pd.concat((freqs_tp2, gene_exp_tp2), axis=1)
    epsilon = 1e-6
    tot_data_tp1 = np.log(tot_data_tp1 + epsilon).fillna(0)
    tot_data_tp2 = np.log(tot_data_tp2 + epsilon).fillna(0)
    data_in = pd.concat((data_in, tot_data_tp1), axis=0)
    data_out = pd.concat((data_out, tot_data_tp2), axis=0)


train_variances = pd.DataFrame()
test_variances = pd.DataFrame()

normalized_data_in = (data_in-data_in.mean())/data_in.std()
normalized_data_out = (data_out-data_out.mean())/data_out.std()
baseline_variance = normalized_data_out.var()
for i in range(20):
    model, train_variance, test_variance = get_linear_model(normalized_data_in, normalized_data_out, alpha=0.15)
    train_variances = pd.concat((train_variances, (train_variance / baseline_variance).to_frame()), axis=1)
    test_variances = pd.concat((test_variances, (test_variance / baseline_variance).to_frame()), axis=1)

print(test_variances[(test_variances.max(axis=1) < 0.7)].mean(axis=1))
reliable_predictions = np.nonzero(np.array(test_variances.max(axis=1) < 0.7).reshape((-1,)))[0]
for r in reliable_predictions:
    relevant_columns = np.argpartition(np.abs(model.coef_[r, :]), -5)[-5:]
    print(data_in.columns[r])
    print(data_in.columns[relevant_columns])
    print(model.coef_[r, relevant_columns])


# print(model.coef_)
'''Step 2: (Optional) Adding LR related genes'''

'''Step 3: Low dimensional representation'''
# font = {'weight' : 'bold',
#         'size'   : 12}
# matplotlib.rc('font', **font)
#
# for ct in cell_types:
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     print(ct)
#     tot_data = np.zeros((0, len(salient_genes)))
#     data = {}
#     plots = []
#     plot_names = []
#     for tg in timepoint_groups:
#         data[tg] = np.array(tot_express[tg][ct])
#         tot_data = np.concatenate((tot_data, data[tg]), axis=0)
#     tot_data = np.log10(tot_data + 0.00001)
#     n_components = 3
#     transformer = PCA(n_components=n_components, random_state=0)
#     transformer.fit(tot_data)
#     for tg in timepoint_groups:
#         data_transformed = transformer.transform(data[tg])
#         crnt_plt, = ax.plot(data_transformed[:, 0], data_transformed[:, 1], data_transformed[:, 2], '*')
#         plots.append(crnt_plt)
#         plot_names.append(tg)
#     ax.legend(plots, plot_names)
#     ax.set_title(ct)
#     plt.savefig(target_path + 'PCA_' + ct + '.png', dpi=300)
#
#


'''Step 4: bin counting'''
