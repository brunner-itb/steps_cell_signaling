import numpy as np
from Patrick.src.MeshProcessor import create_full_mesh

# ellipsoidity = np.linspace(0,1,11)
# for ellip in ellipsoidity:
ellip = 0
output_file = f"/home/pb/steps_cell_signaling/Patrick/meshes_ellipsoidity/ellipsoidity_{ellip}_TEST.inp"
create_full_mesh(output_file, ellip, mesh_size_min=0.166e-6 * 4, mesh_size_max=0.4e-6 * 4)