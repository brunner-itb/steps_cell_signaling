import sys
import logging
from os import listdir
from os.path import isfile, join, abspath, dirname
import re
import os
# sys.path.append(abspath(join(dirname(__file__), "../")))  # Add project root to path

from Patrick.src.SimManager import SimManager
from Patrick.scripts.ellipsoidity_scan.parameters import p
from Patrick.src.Utilities import get_repo_path

"""
Initialize and run the simulation for cell signaling pathways with different ellipsoidity meshes.

This script sets up multiple simulation manager instances with specified parameters and species,
then loads a small model and executes simulation runs for each ellipsoidity mesh.

The simulation models EGF-EGFR signaling pathway dynamics within ellipsoidal cell meshes of varying shapes.

If you run a "plot_only_run" the result_selectors need to be different whether you want to plot the results on the go or 
save them. Setting "plot_only_run = True" adjusts it accordingly. But: when you do a "plot_only_run" it means there is 
no data saved, so be aware.
"""

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logging.info("Starting the ellipsoidity scan script...")

try:
    # Retrieve the base repository path
    logging.info("Retrieving the repository path...")
    base_path = get_repo_path()
    logging.info(f"Repository path retrieved: {base_path}")

    # Get all ellipsoid mesh files
    ellipsoid_meshes_path = f"{base_path}Patrick/meshes_ellipsoidity/"
    mesh_files = [join(ellipsoid_meshes_path, f) for f in listdir(ellipsoid_meshes_path) 
                 if isfile(join(ellipsoid_meshes_path, f))]
    
    logging.info(f"Found {len(mesh_files)} mesh files to process")

    # Helper function to extract ellipsoidity from mesh filename
    def get_ellipsoidity_from_mesh_name(mesh_path):
        """Extract ellipsoidity value from mesh filename."""
        filename = os.path.basename(mesh_path)
        match = re.search(r'ellipsoidity_(\d+\.?\d*)', filename)
        if match:
            return match.group(1)
        return None

    # Run simulation for each mesh
    for mesh_idx, mesh_path in enumerate(mesh_files[::5]):
        ellipsoidity = get_ellipsoidity_from_mesh_name(mesh_path)
        if ellipsoidity is None:
            logging.warning(f"Could not extract ellipsoidity from {mesh_path}, using mesh_idx instead.")
            mesh_name = f"mesh_{mesh_idx}"
        else:
            mesh_name = f"mesh_{ellipsoidity}"
        logging.info(f"Processing mesh {mesh_idx + 1}/{len(mesh_files)}: {mesh_path}")
        
        sm = SimManager(parameters=p,
                       mesh_path=mesh_path,
                       save_path=f"{base_path}Patrick/saved_objects/ellipsoidity_large_model/{mesh_name}/result",
                       parallel=True,
                       runname="ellipsoidity",
                       plot_only_run=False,
                       replace=True)
        
        logging.info("Loading the model...")
        sm.load_model(type="large_new", mesh_scale=1)
        logging.info("Model loaded successfully.")

        logging.info("Running the simulation...")
        sm.run(replicats=20)
        logging.info(f"Simulation completed for mesh {mesh_idx + 1}")

    logging.info("All simulations completed successfully.")

except Exception as e:
    logging.error(f"An error occurred: {e}", exc_info=True)
