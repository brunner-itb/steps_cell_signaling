import sys
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

# # Get the replicat_id from the command-line argument
# if len(sys.argv) < 2:
#     replicat_id = None
# else:
#     replicat_id = int(sys.argv[1])

sm = SimManager(parameters=p,
                mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/elipsoid_4.5.inp",
                save_path ="/home/pb/steps_cell_signaling/Patrick/saved_objects/full_run/large_model",  #without the .h5 suffix, but full file path please
                parallel = True,
                runname = "full_run",
                plot_only_run = False,
                replace = True) # whether an already existing file should be overwritten or not. Might throw an error if there is an already existing one and this is set to false.

sm.load_model(type="large")

sm.run(replicats = 1)
