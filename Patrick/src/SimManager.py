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
    def __init__(self, parameters, species_names, mesh_path, save_file="saved_objects/initial_run/initial_run.h5", parallel=False, runname: str = "initial_run"):
        """
        Manages the configuration, setup, and execution of biochemical simulations.

        The `SimManager` class provides a high-level interface for initializing,
        configuring, and running simulations based on user-defined parameters.
        It is responsible for loading simulation models, setting up the required environment,
        and executing the simulation runs while saving the results in an HDF5 file.

        Attributes:
            parameters (dict): Dictionary containing the simulation parameters like
                               diffusion constants, reaction rates, etc.
            species_names (list): List of molecular species names involved in the simulation.
            endtime (int): Total simulation runtime, in seconds.
            mesh_path (str): File path to the geometry mesh used in the simulation.
            save_file (str): Path where simulation results will be stored.
            simulation (object): The initialized simulation instance.

        Methods:
            __init__(parameters, species_names, endt, mesh_path, save_file):
                Initializes the simulation manager and sets up the environment.
            _setup_environment():
                Ensures the necessary directory structure exists for saving simulation results.
            load_model(type):
                Loads a simulation model based on the provided type ("small" or "large").
            run(run_id=0):
                Executes the simulation and saves the results to an HDF5 file.
        """

        self.parameters = parameters
        self.species_names = species_names
        self.endtime = self.parameters["endtime"]
        self.mesh_path = mesh_path
        self.save_file = save_file
        self.simulation = None
        self.result_selector = None
        self.mesh = None
        self.cell_tets = None
        self.nuc_tets = None
        self.parallel = parallel
        self.runname = runname

        self._setup_environment()

    def _setup_environment(self):
        """
        Set up directories and environment.
        """

        save_dir = os.path.dirname(self.save_file)  # Get the parent directory of the save file

        if not os.path.exists(save_dir):  # Check if the directory exists
            os.makedirs(save_dir)  # Create it if it doesn't exist


    def load_model(self, type, plot_only_run=False):
        self.plot_only_run = plot_only_run
        if type == "small":
            from src.Model_small import create_model
            self.simulation, self.result_selector, self.mesh = create_model(self.parameters, self.species_names, self.mesh_path, plot_only_run)
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
        if self.plot_only_run == False:
            with stsave.HDF5Handler(self.save_file) as hdf:
                self.simulation.toDB(hdf,self.runname, run_id=run_id)
                self.simulation.newRun()
                self.simulation.exo.EGF.Count = 4e4
                self.simulation.cell_surface.EGFR.Count = 7.8e4
                self.simulation.cyt.GAP.Count = 2.3e4
                # self.simulation.cyt.X.Count = 4.1e4
                # self.simulation.nuc_mem.Xa.DiffusionActive = True

                start_time = time.time()
                self.simulation.run(self.endtime)
                end_time = time.time()

                print(f"Run {run_id} completed in {end_time - start_time:.2f} seconds.")
        else:
            from src.InteractivePlotting import interactive_plots
            # TODO: assert that non parallel run? Test it first
            self.simulation.newRun()
            self.simulation.exo.EGF.Count = 4e4
            self.simulation.cell_surface.EGFR.Count = 7.8e4
            self.simulation.cyt.GAP.Count = 2.3e4

            self.result_selector = stsave.ResultSelector(self.simulation)
            SimControl = interactive_plots(self)
            SimControl.run()

