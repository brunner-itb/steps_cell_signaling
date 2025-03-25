import numpy as np
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from Patrick.src.MeshProcessor import create_full_mesh
from Patrick.src.Utilities import get_repo_path

base_path = get_repo_path()

ellipsoidity = np.linspace(0,1,11)
for ellip in ellipsoidity:
    output_file = base_path + f"Patrick/meshes_ellipsoidity/ellipsoidity_{np.round(ellip,3)}.inp"
    print(output_file)
    create_full_mesh(output_file, ellip, mesh_size_min=0.166e-6, mesh_size_max=0.4e-6)