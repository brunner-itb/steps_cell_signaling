import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
from src.Utilities import molar_to_molecules, nostdout
import os
import time
import sys
import random
import numpy as np
import pandas as pd

def initialize_ellipsoid_mesh(mesh_path, scale, nucleus_volume, cytosol_volume, extracellular_volume, cell_surface_system):
    # Load mesh and compartments
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale)

    with mesh:
        # LISTEN
        # Extracellular space
        exo_tets = stgeom.TetList(mesh.tetGroups["Volume1"])
        # Zelle
        cytosol_tets = stgeom.TetList(mesh.tetGroups["Volume2"])
        # Zellkern
        nuc_tets = stgeom.TetList(mesh.tetGroups["Volume3"])

        # COMPARTMENTS
        # Zellkern
        nuc = stgeom.Compartment(nuc_tets, nucleus_volume, name="nuc")
        # Cytoplasma
        cyt = stgeom.Compartment(cytosol_tets, cytosol_volume, name="cyt")
        # Zelläußeres
        exo = stgeom.Compartment(exo_tets, extracellular_volume, name="exo")
        # Zellmembran
        cell_surface = stgeom.Patch(cyt.surface & exo.surface, cyt, exo, cell_surface_system, name="cell_surface")
        # DIFFUSIONS BARRIERE
        nuc_mem = stgeom.DiffBoundary(nuc.surface, name="nuc_mem")
    return mesh, exo_tets, cytosol_tets, nuc_tets


def create_model(model_dataframe, p, species_names, mesh_path, mesh_scale, plot_only_run):
    mdl = stmodel.Model()
    r = stmodel.ReactionManager()
    data_big_model_mini_sph_df = model_dataframe
    species_dict = {}

    # Volume in L
    V = 1.414e-17 * 1000 #dm^3
    fac = 1/V
    
    # Create volume and surface systems
    with mdl:
        volume_system = stmodel.VolumeSystem.Create()
        nucleus_volume = stmodel.VolumeSystem.Create()
        extracellular_volume = stmodel.VolumeSystem.Create()
        cell_surface = stmodel.SurfaceSystem.Create()

        # Create a dictionary to hold the created species
        for sp_name in species_names:
            species_dict[sp_name] = stmodel.Species(name=sp_name)

        with volume_system:
            # Simplified reactions in cytoplasm
            species_dict["ERK"] + species_dict["MEKpp"] < r[1] > species_dict["ERK_MEKpp"]
            r[1].K = 1.1e5 * fac, 0.033
            species_dict["ERK_MEKpp"] > r[2] > species_dict["ERKp"] + species_dict["MEKpp"]
            r[2].K = 16
            species_dict["ERKp"] + species_dict["MEKpp"] < r[3] > species_dict["ERKp_MEKpp"]
            r[3].K = 1.1e5 * fac, 0.033
            species_dict["ERKp_MEKpp"] > r[4] > species_dict["ERKpp"] + species_dict["MEKpp"]
            r[4].K = 5.7
            species_dict["ERKpp"] + species_dict["P3"] < r[5] > species_dict["ERKpp_P3"]
            r[5].K = 1.45e7 * fac, 0.6
            species_dict["ERKpp_P3"] > r[6] > species_dict["ERKp"] + species_dict["P3"]
            r[6].K = 0.27

            # Cytoplasm diffusion
            for index, row in data_big_model_mini_sph_df.iterrows():
                species_name = row['Species']
                cyt_dc = row['cyt DC']
                if pd.notna(cyt_dc):
                    stmodel.Diffusion(species_dict[species_name], cyt_dc)

        with extracellular_volume:
            for index, row in data_big_model_mini_sph_df.iterrows():
                species_name = row['Species']
                exo_dc = row['exo Volume DC']
                if pd.notna(exo_dc):
                    stmodel.Diffusion(species_dict[species_name], exo_dc)

        with cell_surface:
            # Simplified membrane reactions
            species_dict["EGFR"].s + species_dict["EGF"].o < r[7] > species_dict["EGF_EGFR"].s
            r[7].K = 3e7 * fac, 38e-4
            species_dict["EGF_EGFR"].s + species_dict["EGF_EGFR"].s < r[8] > species_dict["EGF_EGFR2"].s
            r[8].K = 1e7 * fac, 0.1
            species_dict["EGF_EGFR2"].s < r[9] > species_dict["EGF_EGFRp2"].s
            r[9].K = 1, 0.01
            species_dict["EGF_EGFRp2"].s + species_dict["GAP"].i < r[10] > species_dict["EGF_EGFRp2_GAP"].s
            r[10].K = 1e6 * fac, 0.2

            # Cell surface diffusion
            for index, row in data_big_model_mini_sph_df.iterrows():
                species_name = row['Species']
                cell_surface_dc = row['cell_surface DC']
                if pd.notna(cell_surface_dc):
                    stmodel.Diffusion(species_dict[species_name], cell_surface_dc)

    # Load mesh and compartments
    mesh, exo_tets, cytosol_tets, nuc_tets = initialize_ellipsoid_mesh(mesh_path,
                                                                    scale=mesh_scale,
                                                                    nucleus_volume=nucleus_volume,
                                                                    cytosol_volume=volume_system,
                                                                    extracellular_volume=extracellular_volume,
                                                                    cell_surface_system=cell_surface)
    system_volume = mesh.Vol

    ratio_v = mesh.exo.Vol / 1286e-18
    ratio_mini = mesh.cyt.Vol / 1766e-18

    # Initialize RNG and Simulation
    seed = random.randint(1,6000)
    rng = strng.RNG("mt19937", 512, seed)
    partition = stgeom.LinearMeshPartition(mesh, 1, 1, stsim.MPI.nhosts)
    sim = stsim.Simulation("TetOpSplit", mdl, mesh, rng, False, partition)

    # Define results
    rs = None
    if plot_only_run == False:
        rs = stsave.ResultSelector(sim)
        result_selectors = []

        for index, row in data_big_model_mini_sph_df.iterrows():
            species_name = row['Species']
            result_selector = row['resultsselector']
            if result_selector == 'TRIS':
                result_selectors.append((species_name, rs.TRIS(cytosol_tets.surface)))
            elif result_selector == 'TETS':
                result_selectors.append((species_name, rs.TETS(cytosol_tets)))

        for s,r in result_selectors:
            print(f"Setting result selector {s}/{len(result_selectors)}")
            rs_path = rs.SUM(getattr(r, s).Count)
            sim.toSave(rs_path, dt = p["time step"])

    return sim, rs, mesh 