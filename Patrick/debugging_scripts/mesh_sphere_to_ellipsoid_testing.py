import gmsh
import sys
import numpy as np


def create_ellipsoid_surface(center_x, center_y, center_z, r_a, r_b, r_c, mesh_size=0.05):
    """
    Creates an ellipsoid in Gmsh using the Python API.

    Args:
        center_x (float): X-coordinate of the ellipsoid center.
        center_y (float): Y-coordinate of the ellipsoid center.
        center_z (float): Z-coordinate of the ellipsoid center.
        r_a (float): Radius along the X-axis.
        r_b (float): Radius along the Y-axis.
        r_c (float): Radius along the Z-axis.
        mesh_size (float): Mesh size for the points.

    Returns:
        int: The tag of the created volume.
    """
    # gmsh.model.add("ellipsoid")

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

    # Define the volume
    # v1 = gmsh.model.geo.addVolume([sl1])

    gmsh.model.geo.synchronize()
    return sl1


def generate_ellipsoid_radii(ellipsoidity, volume=1):
    """
    Generate an ellipsoid, increasing ellipsoidity transitions from a sphere ( = 0) to an ellipsoid ( > 0)
    while maintaining a specified constant volume.

    Parameters:
    ellipsoidity : float
        Controls how elongated the ellipsoid becomes (higher values mean more elongation).
    volume : float, optional
        The desired volume of the ellipsoid (default is 1).

    Returns:
    List of tuples (r_a, r_b, r_c) representing ellipsoid radii.
    """

    # Initial radius for a sphere with given volume
    r_0 = (3 * volume / (4 * np.pi)) ** (1 / 3)

    factor = 1 + ellipsoidity

    # Define radii while keeping volume constant
    r_a = r_0 * factor
    r_b = r_0 / np.sqrt(factor)
    r_c = r_0 / np.sqrt(factor)

    # Ensure volume remains constant
    computed_volume = (4 / 3) * np.pi * r_a * r_b * r_c
    assert np.isclose(computed_volume, volume, atol=1e-6), "Volume constraint violated!"
    return r_a, r_b, r_c


def create_full_mesh(
        output_file,
        ellipsoidity,
        cell_volume=1676e-18,
        nucleus_volume_ratio=0.2,
        extracellular_volume_offset=3e-6,
        mesh_algorithm=1,
        mesh_size_min=0.166e-6,
        mesh_size_max=0.4e-6,
):
    """
    Creates a 3D mesh of a biological cell, including nucleus, cytosol, and extracellular space,
    while maintaining specified volume constraints and ellipsoidity.

    Parameters:
    output_file (str): Path to save the generated mesh file.
    ellipsoidity (float): Controls the elongation of the cell structure.
    cell_volume (float, optional): Total volume of the cell in cubic meters (default: 1676e-18 m³).
    nucleus_volume_ratio (float, optional): Ratio of nucleus volume to total cell volume (default: 0.2).
    extracellular_volume_offset (float, optional): Offset added to extracellular radii (default: 3e-6 m).
    mesh_algorithm (int, optional): Gmsh meshing algorithm (default: 1).
    mesh_size_min (float, optional): Minimum mesh element size (default: 0.166e-6 m).
    mesh_size_max (float, optional): Maximum mesh element size (default: 0.4e-6 m).
    """

    gmsh.initialize()
    gmsh.model.add(f"mesh_ellipsoidity_{ellipsoidity}")

    # Generate nucleus ellipsoid
    rx, ry, rz = generate_ellipsoid_radii(ellipsoidity, cell_volume * nucleus_volume_ratio)
    nucleus_surface_tag = create_ellipsoid_surface(0, 0, 0, rx, ry, rz, 0.05)
    nucleus_volume_tag = gmsh.model.geo.addVolume([nucleus_surface_tag])

    # Generate cytosol ellipsoid
    rx, ry, rz = generate_ellipsoid_radii(ellipsoidity, cell_volume)
    cytosol_surface_tag = create_ellipsoid_surface(0, 0, 0, rx, ry, rz, 0.05)
    cytosol_volume_tag = gmsh.model.geo.addVolume([cytosol_surface_tag])

    # Generate extracellular space by extending cytosol radii
    extracellular_surface_tag = create_ellipsoid_surface(
        0, 0, 0,
        rx + extracellular_volume_offset,
        ry + extracellular_volume_offset,
        rz + extracellular_volume_offset,
        0.05
    )
    extracellular_volume_tag = gmsh.model.geo.addVolume([extracellular_surface_tag])

    # Define mesh settings
    gmsh.option.setNumber("Mesh.Algorithm", mesh_algorithm)  # Choose meshing algorithm
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min * 4)  # Set minimum mesh size
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max * 4)  # Set maximum mesh size

    gmsh.model.geo.synchronize()

    # Generate and save the mesh
    gmsh.model.mesh.generate(3)
    if output_file[-4:] != ".inp":
        print()
        print("Careful, suggested output file format compatible with STEPS and the rest of this software is .inp")
        print()
    gmsh.write(output_file)

    # Uncomment for visualization/debugging
    # gmsh.fltk.run()

    gmsh.finalize()
    
create_full_mesh("/home/pb/Downloads/bio_data_meshes/fixed_test.inp", 0)

#%%
# # 2. Define the Sphere
# ellipsoid_surface_tag = create_ellipsoid_surface(center[0]/10, center[1]/10, center[2]/10, radii[0], radii[1], radii[2], 0.05)
#
# #
# # Volume of the sphere (enclosing the STL)
# ellipsoid_volume_tag = gmsh.model.geo.addVolume([ellipsoid_surface_tag])
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.embed(2, [stl_surface_loop_tag], 3, ellipsoid_volume_tag)




#
# # Define mesh settings
# gmsh.option.setNumber("Mesh.Algorithm", 1)  # Change meshing algorithm if needed
# gmsh.option.setNumber("Mesh.MeshSizeMin", 3.0)  # Adjust as needed
# gmsh.option.setNumber("Mesh.MeshSizeMax", 8.0)  # Adjust as needed
#
# # 6. Synchronize the GEO Model
# gmsh.model.geo.synchronize()
# # 7. Generate the Mesh
# gmsh.model.mesh.generate(3)
#
# # Save as Gmsh 2.2 ASCII format
# # gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
# gmsh_output = "/home/pb/Downloads/bio_data_meshes/fixed_test.inp"
# gmsh.write(gmsh_output)
# # gmsh.fltk.run()
# # Finalize Gmsh
# gmsh.finalize()
# print(f"3D mesh saved to {gmsh_output}")















# buffered_output = "/home/pb/Downloads/bio_data_meshes/buffered_test.msh"
# new_file = "/home/pb/Downloads/bio_data_meshes/new_test.msh"
# generate_buffer_zone(new_file, buffered_output, d=10.0)
# def generate_buffer_zone(mesh_file, output_file, d):
#     gmsh.initialize(sys.argv)
#     gmsh.model.add("buffer_zone")
#
#     # Load the existing 3D mesh
#     gmsh.merge(mesh_file)
#
#     # Get bounding box of the original mesh
#     nodes = np.array(gmsh.model.mesh.getNodes()[1]).reshape(-1, 3)
#     if nodes.size == 0:
#         print("No nodes found in the mesh file.")
#         gmsh.finalize()
#         return
#
#     min_coords = np.min(nodes, axis=0)
#     max_coords = np.max(nodes, axis=0)
#
#     # Compute ellipsoid radii (expanding by distance d)
#     center = (min_coords + max_coords) / 2
#     radii = (max_coords - min_coords) / 2 + d
#
#     # Create an outer ellipsoid (buffer zone) as a surface
#     # outer_sphere = gmsh.model.occ.addSphere(center[0], center[1], center[2], radii[0], radii[1], radii[2])
#     outer_sphere = gmsh.model.occ.addSphere(center[0], center[1], center[2], 2)
#
#     gmsh.model.occ.synchronize()
#
#     # Define mesh settings
#     gmsh.option.setNumber("Mesh.Algorithm", 1)  # Change meshing algorithm if needed
#     gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)  # Adjust as needed
#     gmsh.option.setNumber("Mesh.MeshSizeMax", 5.0)  # Adjust as needed
#
#     gmsh.model.mesh.generate(3)  # Generate full 3D volume mesh
#
#     # Save with partitions
#     gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
#     gmsh.write(output_file)
#
#     gmsh.finalize()
#     print(f"Buffer zone mesh with 3D volume added to {output_file}.")
#
#
# #
# #
# # import steps.geom as stgeom
# #
# # mesh = stgeom.TetMesh.LoadGmsh(gmsh_output, scale=1)
