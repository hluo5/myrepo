import dpdata
import numpy as np

# load data of abacus/md format
data = dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar")
print("# the data contains %d frames" % len(data))

# random choose 40 index for validation_data
index_validation = np.random.choice(897, size=40, replace=False)

# other indexes are training_data
index_training = list(set(range(897)) - set(index_validation))
data_training = data.sub_system(index_training)
data_validation = data.sub_system(index_validation)

# all training data put into directory:"training_data"
data_training.to_deepmd_npy("training_data")

# all validation data put into directory:"validation_data"
data_validation.to_deepmd_npy("validation_data")

print("# the training data contains %d frames" % len(data_training))
print("# the validation data contains %d frames" % len(data_validation))
