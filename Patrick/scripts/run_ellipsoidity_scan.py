import sys
from os import listdir
from os.path import isfile, join, abspath, dirname
sys.path.append(abspath(join(dirname(__file__), "../")))  # Add project root to path

from src.SimManager import SimManager
from parameters import p


"""
Initialize and run the simulation for cell signaling pathways.

This script sets up a simulation manager instance with specified parameters and species,
then loads a small model and executes a single simulation run with parallel processing.

The simulation models EGF-EGFR signaling pathway dynamics within a spherical cell mesh.

If you run a "plot_only_run" the result_selectors need to be different whether you want to plot the results on the go or 
save them. Setting "plot_only_run = True" adjusts it accordingly. But: when you do a "plot_only_run" it means there is 
no data saved, so be aware.
"""

ellipsoid_meshes_path = "/home/pb/steps_cell_signaling/Patrick/meshes_ellipsoidity/"
mesh_files = [join(ellipsoid_meshes_path, f) for f in listdir(ellipsoid_meshes_path) if isfile(join(ellipsoid_meshes_path, f))]

for mesh_idx, mesh_path in enumerate(mesh_files):
    sm = SimManager(parameters=p,
                    mesh_path = mesh_path,
                    save_path = f"/home/pb/steps_cell_signaling/Patrick/saved_objects/ellipsoidity/mesh_{mesh_idx}/result",  #without the .h5 suffix, but full file path please
                    parallel = True,
                    runname = "ellipsoidity",
                    plot_only_run = False,
                    replace = True) # whether an already existing file should be overwritten or not. Might throw an error if there is an already existing one and this is set to false.
    sm.load_model(type="small")
    sm.run(replicats = 30)
