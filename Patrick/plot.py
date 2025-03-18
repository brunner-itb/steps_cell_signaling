import steps.interface
import steps.saving as stsave
from matplotlib import pyplot as plt
import seaborn as sns
import h5py
import numpy as np
import re
import math

def traverse_datasets(hdf_file):

    """Traverse all datasets across all groups in HDF5 file."""
    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = '{}/{}'.format(prefix, key)
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    with h5py.File(hdf_file, 'r') as f:
        for (path, dset) in h5py_dataset_iterator(f):
            print(path, dset)

    return None


traverse_datasets("/home/pb/steps_cell_signaling/Patrick/saved_objects/initial_run/parallel_run.h5")
#%%
# hdf = stsave.HDF5Handler("/home/pb/steps_cell_signaling/Patrick/saved_objects/initial_run/parallel_run_1")
hdf = stsave.HDF5Handler("/home/pb/steps_cell_signaling/Patrick/saved_objects/full_run/large_model")
# with stsave.HDF5Handler("/home/pb/steps_cell_signaling/Patrick/saved_objects/initial_run/parallel_run_1") as hdf:
# results = hdf["long_run"].results
results = hdf["full_run"].results

# extract the species names for the result_selector label via regex
full_labels = [x.labels for x in results]
species_names = [re.search(r'\.(.*?)\.', s[0]).group(1) for s in full_labels if re.search(r'\.(.*?)\.', s[0])]

# Plot all results
num_species = len(species_names)
grid_size = math.ceil(math.sqrt(num_species))  # Prefer a square layout
n_rows, n_cols = grid_size, math.ceil(num_species / grid_size)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
axes = axes.flatten()  # Flatten in case of 2D array

for idx, (res, species_name) in enumerate(zip(results, species_names)):
    mean_data = np.mean(res.data[:,:,0], axis=0)
    std_data = np.std(res.data[:,:,0], axis=0)
    ax = axes[idx]
    ax.scatter(res.time[0], mean_data, label='Mean', s = 0.5)
    ax.fill_between(res.time[0], mean_data - std_data, mean_data + std_data, alpha=0.3, label='std')

    if idx >= len(results) - n_cols:
        ax.set_xlabel('Time [s]')
    # else:
    #     ax.set_xticklabels([])
    ax.set_ylabel(species_name)
    # ax.legend()

# Hide any unused subplots
for i in range(idx + 1, len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.show()
