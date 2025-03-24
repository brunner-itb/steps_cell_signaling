This repo attempts to streamline the modeling of cell signaling via the STEPS package 
(https://steps.sourceforge.net/STEPS/default.php) and is very much a work in progress.

Before running, add your repo path to the config dict in the `get_repo_path()` function 
found in `Utilities.py`. 
To get your username and hostname, use `socket.gethostname()` and `getpass.getuser()` 
in a Python console. 

To run it in parallel use ``mpirun -n N python3 run.py`` with N as the number of CPU cores
you want to use.

Be aware of the options in ``run.py`` and adjust them accordingly.

If the model runs way too fast and there are no kinetics visible in the results, it is 
possible that the meshes chosen are too low in resolution. The standard resolution is
``mesh_size_min=0.166e-6, mesh_size_max=0.4e-6``.

To implement a new model, define it in its own file (see `src/Model_small.py` for an 
example) and add it to the SimManagers ``load_model`` function. Please avoid the
`.Create()` convenience function whenever possible, it is not safe for testing. If you
run a model and in the meantime change the file that describes the model it is possible
for the solver to fail. This happens because the `.Create()` function reads the file and
the line where it was called during runtime, checks in how many containers it wants to be
cast and creates the corresponding amount of objects. Don't do this. Please. 

