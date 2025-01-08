import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import numpy as np
import os



def create_model(p, species_names, mesh_path):
    """
    Create the model, compartments, and simulation based on p.
    """
    c1 = 1
    c2 = 1
    mdl = stmodel.Model()
    r = stmodel.ReactionManager()
    species_dict = {}

    # Create volume and surface systems
    with mdl:
        volume_system = stmodel.VolumeSystem.Create()
        extracellular_volume = stmodel.VolumeSystem.Create()
        nucleus_volume = stmodel.VolumeSystem.Create()
        cell_surface = stmodel.SurfaceSystem.Create()

        # Create a dictionary to hold the created species
        for sp_name in species_names:
            species_dict[sp_name] = stmodel.Species(name=sp_name)


        # # Volume system
        # with volume_system:
        #     Xa > r[666] > X
        #     r[666].K = p["k[666]"]
        #     stmodel.Diffusion(X, species_dict["DCX"])
        #     stmodel.Diffusion(Xa, species_dict["DCX"])

        # Extracellular volume
        with extracellular_volume:
            stmodel.Diffusion(species_dict["EGF"], p["DC"])

        # Nucleus
        # with nucleus_volume:
        #     Xa > r[777] > XA
        #     r[777].K = 1e8
        #     stmodel.Diffusion(Xa, species_dict["DCX"])

        # Surface system (cell membrane)
        with cell_surface:
            species_dict["EGFR"].s + species_dict["EGF"].o < r[1] > species_dict["EGF_EGFR"].s
            species_dict["EGF_EGFR"].s + species_dict["EGF_EGFR"].s < r[2] > species_dict["EGF_EGFR2"].s
            species_dict["EGF_EGFR2"].s < r[3] > species_dict["EGF_EGFRp2"].s
            species_dict["EGF_EGFRp2"].s + species_dict["GAP"].i < r[5] > species_dict["EGF_EGFRp2_GAP"].s
            species_dict["EGF_EGFRp2_GAP"].s < r[0] > species_dict["EGF_EGFRp2_GAP"].s
            r[0].K = p["k[0]"], 0.1
            r[1].K = 3e7 * c1, 38e-4 * c1
            r[2].K = 1e7 * c2, 0.1 * c2
            r[3].K = 1, 0.01
            r[5].K = 1e6 * c1, 0.2 * c1

            stmodel.Diffusion(species_dict["EGFR"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFR"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFR2"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFRp2"], p["DC"])
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP"], p["DC"])

    # Load mesh and compartments
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=10 ** -6)

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


    # Initialize RNG and Simulation
    rng = strng.RNG("mt19937", 512, 2903)
    simulation = stsim.Simulation("Tetexact", mdl, mesh, rng)


    # Define results
    rs = stsave.ResultSelector(simulation)
    result_selectors = {
        "EGF_EGFR": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFR.Count),
        "EGF_EGFRp2": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2.Count),
        "EGF_EGFRp2_GAP": rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2_GAP.Count),
    }

    # Schedule saving
    for key, sel in result_selectors.items():
        simulation.toSave(sel, dt=p["time step"])

    return simulation