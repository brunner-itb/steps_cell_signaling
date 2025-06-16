import sys
from src.SimManager import SimManager
import logging
from parameters import p
from src.Utilities import get_repo_path

"""
Initialize and run the intermediate simulation for cell signaling pathways.

This script sets up a simulation manager instance with specified parameters and species,
then loads the intermediate model and executes a single simulation run.

The simulation models a simplified version of the EGF-EGFR signaling pathway dynamics 
within a spherical cell mesh, focusing on core ERK pathway reactions.
"""

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Starting the script...")

try:
    # Retrieve the base repository path
    logging.info("Retrieving the repository path...")
    base_path = get_repo_path()
    logging.info(f"Repository path retrieved: {base_path}")
   
    # Initialize the simulation manager
    logging.info("Initializing the simulation manager...")
    sm = SimManager(parameters=p,
                    mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes_ellipsoidity/ellipsoidity_0.0.inp",
                    save_path=f"{base_path}Patrick/saved_objects/testing/intermediate_model_test",  # Full file path without the .h5 suffix
                    parallel = True,
                    runname = "test",
                    plot_only_run = False,
                    replace = True) # whether an already existing file should be overwritten or not
    logging.info("Simulation manager initialized.")

    # Load the model
    logging.info("Loading the intermediate model...")
    sm.load_model(type="intermediate", mesh_scale=1)
    logging.info("Model loaded successfully.")

    # Run the simulation
    logging.info("Running the simulation...")
    sm.run(replicats = 1)
    logging.info("Simulation completed successfully.")

except Exception as e:
    logging.error(f"An error occurred: {e}", exc_info=True) 