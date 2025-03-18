import os
import steps.interface
import steps.geom as stgeom
import numpy as np
import trimesh
import pymeshfix
import gmsh


def fix_surface_holes(import_file, output_file):
    '''
    Fixes potential holes in the surface of a mesh.
    This was written as the meshes provided by the collaboration partners typically contain cut-off areas which result in
    holes in the surface, but meshing programs typically require watertight surface meshes for 3D meshing.
    '''

    # # Load the STL file using Trimesh, e.g.
    # import_file = "/home/pb/Downloads/bio_data_meshes/0022_0001_accelerator_20210315_bakal01_erk_main_21-03-15_12-37-27.off"
    # output_file = "/home/pb/Downloads/bio_data_meshes/fixed_test.stl"
    mesh = trimesh.load(import_file)

    # Ensure the mesh is valid and not empty
    if mesh.is_empty or mesh.faces.shape[0] == 0:
        raise ValueError("The mesh is empty or has no faces!")

    # Use pymeshfix with the vertex and face arrays
    fixer = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
    fixer.repair()
    fixer.save(output_file)



def create_ellipsoid_surface(center_x, center_y, center_z, r_a, r_b, r_c, mesh_size=0.05):
    """
    Creates an ellipsoid surface in Gmsh using the Python API.

    Args:
        center_x (float): X-coordinate of the ellipsoid center.
        center_y (float): Y-coordinate of the ellipsoid center.
        center_z (float): Z-coordinate of the ellipsoid center.
        r_a (float): Radius along the X-axis.
        r_b (float): Radius along the Y-axis.
        r_c (float): Radius along the Z-axis.
        mesh_size (float): Mesh size for the points.

    Returns:
        int: The tag of the created surface loop.
    """

    # Define the points
    p_center = gmsh.model.geo.addPoint(center_x, center_y, center_z, mesh_size)
    p_ax = gmsh.model.geo.addPoint(center_x + r_a, center_y, center_z, mesh_size)
    p_ay = gmsh.model.geo.addPoint(center_x, center_y + r_b, center_z, mesh_size)
    p_az = gmsh.model.geo.addPoint(center_x, center_y, center_z + r_c, mesh_size)
    p_nx = gmsh.model.geo.addPoint(center_x - r_a, center_y, center_z, mesh_size)
    p_ny = gmsh.model.geo.addPoint(center_x, center_y - r_b, center_z, mesh_size)
    p_nz = gmsh.model.geo.addPoint(center_x, center_y, center_z - r_c, mesh_size)

    # Define the ellipses
    l1 = gmsh.model.geo.addEllipseArc(p_ax, p_center, p_nz, p_nz)
    l2 = gmsh.model.geo.addEllipseArc(p_nz, p_center, p_nx, p_nx)
    l3 = gmsh.model.geo.addEllipseArc(p_nx, p_center, p_az, p_az)
    l4 = gmsh.model.geo.addEllipseArc(p_az, p_center, p_ax, p_ax)
    l5 = gmsh.model.geo.addEllipseArc(p_ax, p_center, p_ay, p_ay)
    l6 = gmsh.model.geo.addEllipseArc(p_ay, p_center, p_nx, p_nx)
    l7 = gmsh.model.geo.addEllipseArc(p_nx, p_center, p_ny, p_ny)
    l8 = gmsh.model.geo.addEllipseArc(p_ny, p_center, p_ax, p_ax)
    l9 = gmsh.model.geo.addEllipseArc(p_nz, p_center, p_ay, p_ay)
    l10 = gmsh.model.geo.addEllipseArc(p_ay, p_center, p_az, p_az)
    l11 = gmsh.model.geo.addEllipseArc(p_az, p_center, p_ny, p_ny)
    l12 = gmsh.model.geo.addEllipseArc(p_ny, p_center, p_nz, p_nz)

    # Define the line loops
    ll1 = gmsh.model.geo.addCurveLoop([l5, l10, l4])
    ll2 = gmsh.model.geo.addCurveLoop([l9, -l5, l1])
    ll3 = gmsh.model.geo.addCurveLoop([-l10, l6, l3])
    ll4 = gmsh.model.geo.addCurveLoop([-l6, -l9, l2])
    ll5 = gmsh.model.geo.addCurveLoop([l8, -l4, l11])
    ll6 = gmsh.model.geo.addCurveLoop([l12, -l8, -l1])
    ll7 = gmsh.model.geo.addCurveLoop([-l11, -l3, l7])
    ll8 = gmsh.model.geo.addCurveLoop([-l2, -l7, -l12])
    #
    # Create surfaces from the curve loops
    s1 = gmsh.model.geo.addSurfaceFilling([ll1])
    s2 = gmsh.model.geo.addSurfaceFilling([ll2])
    s3 = gmsh.model.geo.addSurfaceFilling([ll3])
    s4 = gmsh.model.geo.addSurfaceFilling([ll4])
    s5 = gmsh.model.geo.addSurfaceFilling([ll5])
    s6 = gmsh.model.geo.addSurfaceFilling([ll6])
    s7 = gmsh.model.geo.addSurfaceFilling([ll7])
    s8 = gmsh.model.geo.addSurfaceFilling([ll8])

    # Define the surface loop
    sl1 = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7, s8])

    gmsh.model.geo.synchronize()
    return sl1


def extrude_cytosol_include_nucleus(input_file, output_file, height):

    """
    Extrudes the cytosol and includes the nucleus in a 3D mesh from an input STL file using Gmsh.

    This function takes an STL file representing a biological surface mesh, typically the output from
    fix_surface_holes(), extrudes an extracellular layer, and generates a 3D volumetric mesh that includes the nucleus.
    The output must be saved in Gmsh's `.msh` format.

    Parameters:
    -----------
    input_file : str
        Path to the input STL file (typically the output from a surface hole-fixing process).
    output_file : str
        Path to save the output 3D mesh in `.msh` format.
    height : float
        The extrusion height (negative value) to create an extracellular layer.


    Notes:
    ------
    - The function assumes the STL file represents a closed surface.
    - The output mesh is stored in Gmsh 2.2 ASCII format for compatibility with STEPS.

    Example Usage:
    --------------
    ```python
    extrude_cytosol_include_nucleus("input.stl", "output.msh", height=5.0)
    ```
    """

    assert output_file[-4:] == ".msh", "We need to write to Gmsh 2.2 .msh file for STEPS to be able to read it, please choose the output_file accordingly"
    gmsh.initialize()
    gmsh.model.add("stl_mesh")

    # 1. Merge the STL file (loads it into Gmsh)
    gmsh.merge(input_file)

    # 2. Identify Surfaces from STL
    stl_surfaces = gmsh.model.getEntities(2)
    if not stl_surfaces:
        print("No surfaces found in the STL file, quitting...")
        gmsh.finalize()

    # 3. Create Surface Loop for STL (assuming it's a closed surface)
    stl_surface_tags = [s[1] for s in stl_surfaces]
    stl_surface_loop_tag = gmsh.model.geo.addSurfaceLoop(stl_surface_tags)
    stl_volume_tag = gmsh.model.geo.addVolume([stl_surface_loop_tag], 1)
    gmsh.model.geo.synchronize()

    # 4. Extrude the extracellular space from the cell surface mesh
    extruded_surface_tag, extruded_volume_tag = gmsh.model.geo.extrudeBoundaryLayer(stl_surfaces, heights=[-height], recombine=True)
    gmsh.model.geo.synchronize()

    # 5. Define mesh settings
    gmsh.option.setNumber("Mesh.Algorithm", 1)  # Change meshing algorithm if needed
    gmsh.option.setNumber("Mesh.MeshSizeMin", 3.0)  # Adjust as needed
    gmsh.option.setNumber("Mesh.MeshSizeMax", 8.0)  # Adjust as needed

    # 6. Generate the Mesh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # Save as Gmsh 2.2 ASCII format
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_file)
    gmsh.fltk.run()
    # Finalize Gmsh
    gmsh.finalize()
    print(f"3D mesh saved to {output_file}")