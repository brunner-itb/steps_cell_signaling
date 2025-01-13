import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import warnings
import numpy as np
import os
import time


class SimManager:
    def __init__(self, parameters, species_names, endt, mesh_path, save_file="saved_objects/initial_run/initial_run.h5"):
        """
        Initialize the simulation model.

        Args:
            parameters (dict): Dictionary containing parameters like diffusion constants,
                               reaction rates, etc.
            endt (float): The time for which the simulation will run.
            mesh_path (str): Path to the mesh file for the model geometry.
            save_dir (str): Directory where results will be saved.
        """
        self.parameters = parameters
        self.species_names = species_names
        self.end_time = endt
        self.mesh_path = mesh_path
        self.save_file = save_file
        self.simulation = None
        self.result_selector = None
        self.result_selectors = None
        self.cell_tets = None
        self.nuc_tets = None

        self._setup_environment()

    def _setup_environment(self):
        """
        Set up directories and environment.
        """
        if not os.path.exists(self.save_file):
            os.makedirs(self.save_file)

    def load_model(self, type):
        if type == "small":
            from src.Model_small import create_model
            self.simulation, self.result_selector = create_model(self.parameters, self.species_names, self.mesh_path)
        elif type == "large":
            warnings.warn("The 'large' model type is not yet implemented", UserWarning)
        else:
            warnings.warn(f"The '{type}' model type is not yet implemented", UserWarning)

    def run(self, run_id=0):
        """
        Run the simulation with the initialized parameters.
        Save to HDF.

        Args:
            run_id (int): Identifier for the current simulation run.
        """
        with stsave.HDF5Handler(self.save_file) as hdf:
            self.simulation.toDB(hdf, f'initial_run', run_id=run_id)
            self.simulation.newRun()
            self.simulation.exo.EGF.Count = 4e4
            self.simulation.cell_surface.EGFR.Count = 7.8e4
            self.simulation.cyt.GAP.Count = 2.3e4
            # self.simulation.cyt.X.Count = 4.1e4
            # self.simulation.nuc_mem.Xa.DiffusionActive = True

            start_time = time.time()
            self.simulation.run(self.endt)
            end_time = time.time()

            print(f"Run {run_id} completed in {end_time - start_time:.2f} seconds.")

