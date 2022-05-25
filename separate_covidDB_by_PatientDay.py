import scanpy

samples_path = '../CNP0001250/'
db_name = 'adata_CNP0001250_filt.h5ad'
data = scanpy.read_h5ad(samples_path+db_name)

patients = set(data.obs['patient'])

for p in patients:
    patient_data = data[data.obs['patient'] == p]
    days = set(patient_data.obs['days'])
    for d in days:
        patient_day_data = patient_data[patient_data.obs['days'] == d]
        condition = patient_day_data.obs['condition'][0]
        file_name = 'patient'+ p + '_day' + f'{int(d):02d}' + '_' + condition
        scanpy.write(file_name, patient_day_data)


