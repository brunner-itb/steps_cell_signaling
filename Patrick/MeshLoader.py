import os
import steps.interface
import steps.geom as stgeom

def process_mesh(mesh_path, scale):
    """
    Process a mesh file, checking if it's '.inp' or '.xml' format, and handle it accordingly.

    Args:
        mesh_path (str): Path to the mesh file.

    Returns:
        stgeom.TetMesh: Loaded mesh object.
    """
    assert os.path.isfile(mesh_path), f"Error: The file '{mesh_path}' does not exist. Please check the path."

    file_ext = os.path.splitext(mesh_path)[-1].lower()
    file_name = os.path.splitext(mesh_path)[0]

    if file_ext == ".inp":
        # Check if an equivalent .xml file exists
        xml_path = file_name + ".xml"

        if os.path.isfile(xml_path):
            print(f"Warning: A .xml version of the file '{mesh_path}' exists. Loading '{xml_path}' instead.")
            return stgeom.TetMesh.Load(file_name, scale=scale)
        else:
            # Load .inp and save as .xml
            print(f"Loading .inp file: '{mesh_path}' and saving it as .xml.")
            mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale)
            stgeom.TetMesh.Save(mesh, file_name)
            return mesh

    elif file_ext == ".xml":
        # Directly load .xml file
        print(f"Loading .xml file: '{mesh_path}'.")
        return stgeom.TetMesh.Load(mesh_path, scale=scale)

    else:
        raise ValueError(f"Unsupported file format: '{file_ext}'. Only '.inp' and '.xml' files are supported.")