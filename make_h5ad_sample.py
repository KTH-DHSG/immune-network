import scanpy
import random

samples_path = '../120Samples/OriginalData/'
sample_name = 'adata_all.h5ad'

target_path = '../120Samples/OriginalData/'
target_name = 'adata_small_sample.h5ad'

n_samples = 50000
annotated_data = scanpy.read_h5ad(samples_path + sample_name)
rows = random.sample(range(annotated_data.shape[0]), n_samples)

small_annotated_data = annotated_data[rows]
scanpy.write(target_path+target_name, small_annotated_data)