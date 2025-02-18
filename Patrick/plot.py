import steps.interface
import steps.saving as stsave
from matplotlib import pyplot as plt
import h5py
import numpy as np

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


traverse_datasets("/home/pb/steps_cell_signaling/Patrick/saved_objects/initial_run/parallel_run.h5.h5")
#%%
with stsave.HDF5Handler("/home/pb/steps_cell_signaling/Patrick/saved_objects/initial_run/parallel_run.h5") as hdf:
    A, B, C = hdf['parallel_run'].results
    names = ["A", "B", "C"]

    # Membrane potential
    for idx, res in enumerate([A, B, C]):
        plt.figure(figsize=(10, 7))
        time = []
        for r in range(len(res.time)):
            time.append(res.time[r])
            plt.plot(res.time[r], 1e3 * res.data[r, :, 0], label=f'Compound No. {r}')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel(names[idx])
        plt.show()
