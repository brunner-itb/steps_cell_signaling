import meshio
import gmsh
input_file = "/Users/evelynstangneth/Signaling_repo/mnt/steps_cell_signaling/Output/meshes_ellipsoidity/ellipsoidity_0.0_TEST.inp"

mesh = meshio.read(input_file)
mesh.write(input_file[:-4] + ".stl")

gmsh.initialize()

gmsh.model.add("stl_mesh")

# 1. Merge the STL file (loads it into Gmsh)
gmsh.merge(input_file[:-4] + ".stl")
gmsh.fltk.run()
gmsh.finalize()