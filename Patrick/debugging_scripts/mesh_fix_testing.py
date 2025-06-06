import trimesh
import pymeshfix
import gmsh
import meshio

# Load the STL file using Trimesh
file_path = "/home/pb/Downloads/bio_data_meshes/0098_0155_accelerator_20210315_bakal01_erk_main_21-03-15_12-37-27.off"
save_file = file_path[:-4] + ".stl"
mesh = trimesh.load(file_path)

# Ensure the mesh is valid and not empty
if mesh.is_empty or mesh.faces.shape[0] == 0:
    raise ValueError("The mesh is empty or has no faces!")

# Use pymeshfix with the vertex and face arrays
fixer = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
fixer.repair()
fixer.save(save_file)

gmsh.initialize()
gmsh.model.add("stl_mesh")


gmsh.merge(save_file)
gmsh.fltk.run()
# Finalize Gmsh
gmsh.finalize()


#%%

# Load and inspect the contents of the provided .inp file

inp_file_path = "/home/pb/steps_cell_signaling/Patrick/meshes_ellipsoidity/ellipsoidity_0_TEST.inp"

# Read the file to understand its structure
with open(inp_file_path, "r") as f:
    inp_contents = f.readlines()

# Extract relevant sections from the .inp file

# Identify important sections
node_section = []
element_section = []
surface_section = []
current_section = None

for line in inp_contents:
    line = line.strip()

    if line.startswith("*NODE"):
        current_section = "node"
        continue
    elif line.startswith("*ELEMENT"):
        current_section = "element"
        continue
    elif line.startswith("*SURFACE") or "*ELSET" in line:
        current_section = "surface"
        continue
    elif line.startswith("*"):
        current_section = None  # Ignore other sections

    # Store content in the appropriate list
    if current_section == "node":
        node_section.append(line)
    elif current_section == "element":
        element_section.append(line)
    elif current_section == "surface":
        surface_section.append(line)

# Display some extracted data for verification
len(node_section), len(element_section), len(surface_section)

# Extract surface elements and check for potential issues

# Dictionary to store surface elements by their ELSET name
surface_elements = {}
current_surface = None

for line in inp_contents:
    line = line.strip()

    if line.startswith("*ELEMENT") and "type=CPS3" in line:
        # New surface section detected
        current_surface = line.split("ELSET=")[-1] if "ELSET=" in line else None
        surface_elements[current_surface] = []
        continue
    elif line.startswith("*"):
        current_surface = None  # End of surface section

    # Store element data if within a surface section
    if current_surface and line:
        surface_elements[current_surface].append(line)

# Count the number of surface elements
num_surface_elements = {surf: len(elements) for surf, elements in surface_elements.items()}

# Check for duplicate surface triangles (same nodes used more than once)
triangle_counts = {}
for elements in surface_elements.values():
    for element in elements:
        nodes = tuple(sorted(map(int, element.split(",")[1:])))  # Sort nodes for uniqueness
        triangle_counts[nodes] = triangle_counts.get(nodes, 0) + 1

# Find duplicated triangles
duplicated_triangles = {tri: count for tri, count in triangle_counts.items() if count > 1}

len(surface_elements), num_surface_elements, len(duplicated_triangles)


# Extract tetrahedral elements (volume elements)
tet_elements = []
current_tet_section = False

for line in inp_contents:
    line = line.strip()

    if line.startswith("*ELEMENT") and "type=C3D4" in line:
        current_tet_section = True  # Start of tetrahedral elements
        continue
    elif line.startswith("*"):
        current_tet_section = False  # End of tetrahedral section

    if current_tet_section and line:
        tet_elements.append(tuple(map(int, line.split(",")[1:])))

# Check if any tetrahedral element shares all three nodes with a surface triangle
tetrahedral_set = set(tet_elements)
invalid_tets = []

for tri in triangle_counts.keys():
    for tet in tetrahedral_set:
        if not set(tri).issubset(set(tet)):  # Triangle fully contained in tetrahedron
            invalid_tets.append(tri)
            break  # No need to check further for this triangle

print(len(tet_elements), len(invalid_tets))


#%%
import meshio
import gmsh
input_file = "/home/pb/steps_cell_signaling/Patrick/meshes/mini_sph.inp"

mesh = meshio.read(input_file)
mesh.write(input_file[:-4] + ".stl")

gmsh.initialize()

gmsh.model.add("stl_mesh")

# 1. Merge the STL file (loads it into Gmsh)
gmsh.merge(input_file[:-4] + ".stl")
gmsh.fltk.run()
gmsh.finalize()

#%%


