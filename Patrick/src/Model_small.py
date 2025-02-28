import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
from src.Utilities import molar_to_molecules
import numpy as np
import os


def initialize_sphere_mesh(mesh_path, scale, nucleus_volume, volume_system, extracellular_volume, cell_surface):
    # Load mesh and compartments
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale)

    with mesh:
        # Zelle
        cell_tets = stgeom.TetList(mesh.tetGroups["Volume2"])

        # Zellkern
        nuc_tets = stgeom.TetList(
            (tet for tet in cell_tets if np.linalg.norm(tet.center) < 5e-6))

        mem_tris = cell_tets.surface
        mem_tet = stgeom.TetList()
        for tri in mem_tris:
            for tet in tri.tetNeighbs:
                mem_tet.append(tet)
        mem_tet = stgeom.TetList([tet for tet in mem_tet if tet in cell_tets])

        # Create compartments
        # nucleus
        nuc = stgeom.Compartment.Create(nuc_tets, nucleus_volume)
        # Cytoplasm
        cyt = stgeom.Compartment.Create(cell_tets - nuc_tets, volume_system)
        # Extracellular volume
        exo = stgeom.Compartment.Create(mesh.tetGroups["Volume3"], extracellular_volume)
        # Cell membrane
        cell_surface = stgeom.Patch.Create(cell_tets.surface, cyt, exo, cell_surface)

    return mesh, cell_tets

def initialize_ellipsoid_mesh(mesh_path, scale, nucleus_volume, volume_system, extracellular_volume, cell_surface):
    # Load mesh and compartments
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale)

    with mesh:

        # LISTEN

        # Extracellular space
        exo_tets = stgeom.TetList(mesh.tetGroups["Volume1"])

        # Zelle
        cell_tets = stgeom.TetList(mesh.tetGroups["Volume2"])

        # Zellkern
        nuc_tets = stgeom.TetList(mesh.tetGroups["Volume3"])

        # TETS im Cyt, die an die Membran grenzen
        mem_tris = cell_tets.surface
        mem_tet = stgeom.TetList()
        for tri in mem_tris:
            for tet in tri.tetNeighbs:
                mem_tet.append(tet)
        mem_tet = stgeom.TetList([tet for tet in mem_tet if tet in cell_tets])

        # Hälfte Cyto
        half_cyt = stgeom.TetList(tet for tet in cell_tets if tet.center.y > 0)
        # TRIS vom Durchschnitt
        slice_cyt = stgeom.TriList(half_cyt.surface & (cell_tets - half_cyt).surface)
        # die Tets am Querschnitt
        slice_tet_cyt = stgeom.TetList()
        triInds_cyt = []
        for tri in slice_cyt:
            for tet in tri.tetNeighbs:
                slice_tet_cyt.append(tet)
                triInds_cyt.append(tri.idx)

        # Häfte Nuc
        half_nuc = stgeom.TetList(tet for tet in nuc_tets if tet.center.y > 0)
        # TRIS vom Durchscnitt
        slice_nuc = stgeom.TriList(half_nuc.surface & (nuc_tets - half_nuc).surface)
        # die Tets am Querschnitt
        slice_tet_nuc = stgeom.TetList()
        triInds_nuc = []
        for tri in slice_nuc:
            for tet in tri.tetNeighbs:
                slice_tet_nuc.append(tet)
                triInds_nuc.append(tri.idx)


        # COMPARTMENTS
        # Zellkern
        nuc = stgeom.Compartment.Create(nuc_tets, nucleus_volume)

        # Cytoplasma
        cyt = stgeom.Compartment.Create(cell_tets, volume_system)

        # Zelläüßeres
        exo = stgeom.Compartment.Create(exo_tets, extracellular_volume)

        # Zellmembran
        cell_surface = stgeom.Patch.Create(cyt.surface & exo.surface, cyt, exo, cell_surface)

    return mesh, cell_tets


def create_model(p, species_names, mesh_path, plot_only_run):
    """
    Creates a STEPS simulation model based on parameters, species, and mesh geometry.

    Parameters:
    -----------
    p : dict
        Simulation parameters including diffusion constants (`p["DC"]`), reaction rates,
        and time step (`p["time step"]`).

    species_names : list of str
        Names of chemical species to be defined in the model.

    mesh_path : str
        Path to the Abaqus mesh file defining system geometry, scaled to micrometers.

    Returns:
    --------
    simulation : steps.sim.Simulation
        Configured STEPS simulation object.

    rs : steps.saving.ResultSelector
        Object for selecting and storing simulation results.

    Raises:
    -------
    AssertionError:
        If the provided `mesh_path` is invalid.

    Notes:
    ------
    - Defines species diffusion, surface and volumetric reactions (e.g., EGF-EGFR binding).
    - Partitions the mesh into compartments: nucleus, cytoplasm, extracellular space,
      and membrane patch.
    - Sets up a stochastic simulation engine with result selectors for species counts.
    """
    mdl = stmodel.Model()
    r = stmodel.ReactionManager()
    species_dict = {}

    # Create volume and surface systems
    with mdl:
        volume_system = stmodel.VolumeSystem.Create()
        nucleus_volume = stmodel.VolumeSystem.Create()
        extracellular_volume = stmodel.VolumeSystem.Create()
        cell_surface = stmodel.SurfaceSystem.Create()

        # Create a dictionary to hold the created species
        for sp_name in species_names:
            species_dict[sp_name] = stmodel.Species(name=sp_name)

        # Extracellular volume
        with extracellular_volume:
            stmodel.Diffusion(species_dict["EGF"], p["DC"])

        # Load mesh and compartments
        mesh, cell_tets = initialize_ellipsoid_mesh(mesh_path,
                               scale=10 ** -6,
                               nucleus_volume=nucleus_volume,
                               volume_system=volume_system,
                               extracellular_volume=extracellular_volume,
                               cell_surface=cell_surface)
        system_volume = mesh.Vol
        # Surface system (cell membrane)
        with cell_surface:
            species_dict["EGFR"].s + species_dict["EGF"].o < r[1] > species_dict["EGF_EGFR"].s
            species_dict["EGF_EGFR"].s + species_dict["EGF_EGFR"].s < r[2] > species_dict["EGF_EGFR2"].s
            species_dict["EGF_EGFR2"].s < r[3] > species_dict["EGF_EGFRp2"].s
            species_dict["EGF_EGFRp2"].s + species_dict["GAP"].i < r[5] > species_dict["EGF_EGFRp2_GAP"].s
            # species_dict["EGF_EGFRp2_GAP"].s < r[0] > species_dict["EGF_EGFRp2_GAP"].s
            # r[0].K = p["k[0]"], 0.1
            r[1].K = 3e7, 38e-4
            r[2].K = 1e7, 0.1
            r[3].K = 1, 0.01
            r[5].K = 1e6, 0.2

            stmodel.Diffusion(species_dict["EGFR"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFR"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFR2"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFRp2"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP"], p["DC"])


    # Initialize RNG and Simulation
    rng = strng.RNG("mt19937", 512, 2903)
    partition = stgeom.LinearMeshPartition(mesh, 1, 1, stsim.MPI.nhosts)
    simulation = stsim.Simulation("TetOpSplit", mdl, mesh, rng, False, partition)


    # Define which results to save and how they are processed
    rs = None
    if plot_only_run == False:
        rs = stsave.ResultSelector(simulation)
        result_selectors = {
            "EGF_EGFR": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFR.Count),
            "EGF_EGFRp2": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2.Count),
            "EGF_EGFRp2_GAP": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2_GAP.Count),
        }
        # Schedule saving
        for key, sel in result_selectors.items():
            simulation.toSave(sel, dt=p["time step"])

    return simulation, rs, mesh