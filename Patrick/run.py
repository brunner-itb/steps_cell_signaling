from src.SimManager import SimManager
from parameters import p

# List of molecular species involved in the simulation
species_names = ["EGF", "EGFR", "Xa", "XA", "EGF_EGFR", "EGF_EGFR2",
                 "EGF_EGFRp2", "EGF_EGFRp2_GAP", "GAP"]

"""
Initialize and run the simulation for cell signaling pathways.

This script sets up a simulation manager instance with specified parameters and species,
then loads a small model and executes a single simulation run with parallel processing.

The simulation models EGF-EGFR signaling pathway dynamics within a spherical cell mesh.

If you run a "plot_only_run" the result_selectors need to be different whether you want to plot the results on the go or 
save them. Setting "plot_only_run = True" adjusts it accordingly. But: when you do a "plot_only_run" it means there is 
no data saved, so be aware.
"""

sm = SimManager(parameters=p,
                species_names=species_names,
                mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/elipsoid_4.5.inp",
                save_file ="saved_objects/initial_run/parallel_run",
                parallel = True,
                runname = "parallel_run")

sm.load_model(type="small",
              plot_only_run=True)

sm.run(run_id=0)
