import numpy as np
import sys
from os.path import isfile, join, abspath, dirname
sys.path.append(abspath(join(dirname(__file__), "../../")))  # Add project root to path
from Patrick.src.MeshProcessor import create_full_mesh



ellipsoidity = np.linspace(0,1,11)
for ellip in ellipsoidity:
   # pathname server from local /Users/evelynstangneth/Signaling_repo/mnt/steps_cell_signaling/Output/meshes_ellipsoidity
   # pathname server from server "/home/evelyn/shared_files/signaling_repo/steps_cell_signaling/Output/meshes_ellipsoidity/ellipsoidity_{ellip}_TEST.inp"
    output_file = f"/Users/evelynstangneth/Signaling_repo/mnt/steps_cell_signaling/Output/meshes_ellipsoidity/ellipsoidity_{ellip}_TEST.inp"
    create_full_mesh(output_file, ellip, mesh_size_min=0.166e-6 * 4, mesh_size_max=0.4e-6 * 4)
