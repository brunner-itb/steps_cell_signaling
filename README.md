This repo attempts to streamline the modeling of cell signaling via the STEPS package 
(https://steps.sourceforge.net/STEPS/default.php) and is very much a work in progress.
1. **Dynamic Path Resolution**:
cd into steps_cell_signaling 
run `pwd`to get your path that you should add to your Python path with  `export PYTHONPATH=yourpwdpath` example `export PYTHONPATH=/Users/evelynstangneth/Signaling_repo/mnt/steps_cell_signaling`

2. **Dynamic Path Resolution**:
create a `config.json` in the following file stucture

3. **File Structure**:
   - Ensure the following directory structure is in place:
     ```
     Signaling_repo/
     └── mnt/
         └── steps_cell_signaling/
             ├── Patrick/
             │   ├── src/
             │   │   ├── __init__.py
             │   │   ├── SimManager.py
             │   │   ├── Utilities.py
             │   ├── meshes_ellipsoidity/
             │   │   ├── example.inp
             │   ├── saved_objects/
             │   │   └── testing/
             │   ├── run.py
             ├── config.json
     ```

Before running, add your repo path to the `config.json` in the 
 with the following form:
```json
{
    "users": {
        "username": {
            "hostname": "your pwd python path"
        }
            }
}
```
used in `Utilities.py`. 
To get your hostname and username, use `socket.gethostname()` and `getpass.getuser()` 
in a Python console. 

4. **Run the run.py**:
To run it in parallel use ``mpirun -n N python3 run.py`` with N as the number of CPU cores
you want to use. You can see the number of CPU cores by running htop.

Be aware of the options in ``run.py`` and adjust them accordingly.


To implement a new model, define it in its own file (see `src/Model_small.py` for an 
example) and add it to the SimManagers ``load_model`` function. Also provide an .xls file that specifies
the species, their initial conditions and diffusion constants. In `src/Utilities.py -> dataframe_cleanup()` 
this file is processed. Please adjust this according to the .xls you provide, it might need extension or
might not be necessary at all.

Please avoid the `.Create()` convenience function whenever possible, it is not safe for testing. If you
run a model and in the meantime change the file that describes the model it is possible
for the solver to fail. This happens because the `.Create()` function reads the file and
the line where it was called during runtime, checks in how many containers it wants to be
cast in and creates the corresponding amount of objects. Don't do this. Please. 

<h3> Sim Manager </h3>

The SimManager class is intended to organize all of the data necessary to run a simulation and provide a
structure for future expansion of the model. Right now all of the variables are available and writeable 
(no setter and getter functions), so be cautious if you change them on the fly.

<h3> Model file </h3>

The Model files (see `src/Model_small.py` for an example) contain the loading of the mesh, all of the important 
reactions and what should be saved to the dataframe in the end. The meshes are crucial, hard to create and
need to follow obscure STEPS rules. When providing your own mesh, take great care to make sure all of the
mesh compartments and patches are set properly, the surfaces are what you expect and the DiffusionBarrier works, 
if applicable.

If you use the "TetOpSplit" solver, you need to partition the mesh first. There are multiple partitioning 
algorithms/methods available, the easiest is the LinearMeshPartition. For the rest check the STEPS documentation.
If you use the LinearMeshPartition, take care that you split the mesh in a direction that makes sense. You
want to minimize the overlap between partitions for efficiency sake. E.g. if you have an ellipsoid with the
largest radius in the x direction, it makes sense to split the mesh in the x direction and not the y direction.

When deciding which results to save be aware that just saving everything has a computational cost, especially
on very large systems. Additionally, if you want to look at the concentrations in ParaView, make sure to save
the concentration (e.g. `rs.TETS(example_tets).MySpecies.Conc`) and not the count. If you do not save any .Conc
at all, the .xmf files will only show the mesh. For a more detailed explanation how to use ParaView check the STEPS
documentation.

<h3> Mesh Processor </h3>

This file contains multiple functions related to meshing. The `fix_surface_holes()` function patches the holes
often found on the cell meshes provided by the Bakal lab by just patching a flat surface over them.
The `create_ellipsoid_surface()` function creates an ellipsoid from points, then curves and then surface loops
using the Gmsh geo kernel. `create_full_mesh()` then uses this to create ellipsoids. These functions use the geo
kernel as the occ kernel is not able to load the cell meshes from the Bakal lab, so it will not be useful in the
future. 

Regarding the Bakal lab meshes, to create a 3 layer mesh with extracellular volume, cytosol and nucleus the
`extrude_exo_from_cytosole_create_3D_mesh()` function can be used as a starting point. This is a work in progress!
It extrudes the extracellular space from the mesh surface and creates a 3D mesh from this, but is lacking the nucleus
at the moment.

<h3> Interactive Plotting </h3>

This is used for the "plot_only_runs" and is dictated by STEPS. Here a SimControl instance is created that is able
to dynamically plot and control the simulation during serial runs. You can define which species should be plotted
in which window and which mesh compartment should be visible as well. This requires pyqt5 to work.








