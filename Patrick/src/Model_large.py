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
            #i
            species_dict["EGFRi"] + species_dict["EGFi"] < r[10] > species_dict["EGF_EGFRi"]
            r[10].K = 1.4e5 * fac, 0.011
            species_dict["EGF_EGFRi"] + species_dict["EGF_EGFRi"] < r[11] > species_dict["EGF_EGFR2i"]
            r[11].K = 1e7 * fac, 0.1
            species_dict["EGF_EGFR2i"] < r[12] > species_dict["EGF_EGFRp2i"]
            r[12].K = 1 , 0.01  #1/s
            species_dict["EGF_EGFR2i"] + species_dict["GAP"] < r[14] > species_dict["EGF_EGFRp2_GAPi"]
            r[14].K = 1e6 * fac, 0.2  # 1/Ms
            species_dict["Prot"] > r[15] > species_dict["Proti"]
            r[15].K = 1e4  #1/s
            species_dict["EGFRi"] > r[60] > None
            r[60].K = 6.67e-4 # 1/s
            species_dict["EGFi"] > r[61] > None
            r[61].K = 1.67e-4 #1/s
            species_dict["EGF_EGFRp2i"] > r[62] > None
            r[62].K = 6.67e-4 #1/s

            #Grb2
            species_dict["Grb2_Sos"] < r[35] > species_dict["Grb2"] + species_dict["Sos"]
            r[35].K = 0.0015 , 4.5e6 * fac # 1/s

            #Shc
            species_dict["Shcp_Grb2_Sos"] < r[33] > species_dict["Shcp"] + species_dict["Grb2_Sos"]
            r[33].K = 0.2, 2.1e7 * fac # 1/Ms
            species_dict["Shcp"] > r[36] > species_dict["Shc"]
            r[36].K = 340 # nM?
            species_dict["Shcp"] + species_dict["Grb2"] < r[38] > species_dict["Shcp_Grb2"]
            r[38].K = 3e7 * fac, 0.055  # 1/Ms
            species_dict["Shcp_Grb2"] + species_dict["Sos"] < r[40] > species_dict["Shcp_Grb2_Sos"]
            r[40].K = 3e7 * fac, 0.064 # 1/Ms

            #Raf
            species_dict["Raf"] + species_dict["Ras_GTP"] < r[28] > species_dict["Raf_Ras_GTP"]
            r[28].K = 1e6 * fac, 0.0053 # 1/MS
            species_dict["Raf_Ras_GTP"] < r[29] > species_dict["Rafp"] + species_dict["Ras_GTPp"]
            r[29].K = 1 , 7e5 * fac # 1/Ms
            species_dict["Rafp"] + species_dict["P1"] < r[42] > species_dict["Rafp_P1"]
            r[42].K = 7.17e7 * fac, 0.2  # 1/Ms
            species_dict["Rafp_P1"] > r[43] > species_dict["Raf"] + species_dict["P1"]
            r[43].K = 1  # 1/Ms
            species_dict["MEK"] + species_dict["Rafp"] < r[44] > species_dict["MEK_Rafp"]
            r[44].K = 1.11e7 * fac, 0.01833   # 1/Ms
            species_dict["MEK_Rafp"] > r[45] > species_dict["MEKp"] + species_dict["Rafp"]
            r[45].K = 3.5  # 1/Ms
            species_dict["MEKp"] + species_dict["Rafp"] < r[46] > species_dict["MEKp_Rafp"]
            r[46].K = 1.11e7 * fac, 0.01833 # 1/Ms
            species_dict["MEKp_Rafp"] > r[47] > species_dict["MEKpp"] + species_dict["Rafp"]
            r[47].K = 2.9  # 1/Ms
            species_dict["MEKpp"] + species_dict["P2"] < r[48] > species_dict["MEKpp_P2"]
            r[48].K = 1.43e7 * fac, 0.8   # 1/Ms
            species_dict["MEKpp_P2"] > r[49] > species_dict["MEKp"] + species_dict["P2"]
            r[49].K = 0.058   # 1/Ms
            species_dict["MEKp"] + species_dict["P2"] < r[50] > species_dict["MEKp_P2"]
            r[50].K = 2.5e5 * fac, 0.5   # 1/Ms
            species_dict["MEKp_P2"] > r[51] > species_dict["MEK"] + species_dict["P2"]
            r[51].K = 0.058  # 1/Ms
            species_dict["ERK"] + species_dict["MEKpp"] < r[52] > species_dict["ERK_MEKpp"]
            r[52].K = 1.1e5 * fac, 0.033  # 1/Ms
            species_dict["ERK_MEKpp"] > r[53] > species_dict["ERKp"] + species_dict["MEKpp"]
            r[53].K = 16  # 1/Ms
            species_dict["ERKp"] + species_dict["MEKpp"] < r[54] > species_dict["ERKp_MEKpp"]
            r[54].K = 1.1e5 * fac, 0.033   # 1/Ms
            species_dict["ERKp_MEKpp"] > r[55] > species_dict["ERKpp"] + species_dict["MEKpp"]
            r[55].K = 5.7  # 1/Ms
            species_dict["ERKpp"] + species_dict["P3"] < r[56] > species_dict["ERKpp_P3"]
            r[56].K = 1.45e7 * fac, 0.6   # 1/Ms
            species_dict["ERKpp_P3"] > r[57] > species_dict["ERKp"] + species_dict["P3"]
            r[57].K = 0.27  # 1/Ms
            species_dict["ERKp"] + species_dict["P3"] < r[58] > species_dict["ERKp_P3"]
            r[58].K = 5e6 * fac, 0.5   # 1/Ms
            species_dict["ERKp_P3"] > r[59] > species_dict["ERKp"] + species_dict["P3"]
            r[59].K = 0.3   # 1/Ms

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
            # Membrane reactions
            species_dict["EGFR"].s + species_dict["EGF"].o < r[1] > species_dict["EGF_EGFR"].s
            r[1].K = 3e7 * fac, 38e-4   # 1/Ms
            species_dict["EGF_EGFR"].s + species_dict["EGF_EGFR"].s < r[2] > species_dict["EGF_EGFR2"].s
            r[2].K = 1e7 * fac, 0.1   # 1/Ms
            species_dict["EGF_EGFR2"].s < r[3] > species_dict["EGF_EGFRp2"].s
            r[3].K = 1 , 0.01  # 1/s
            species_dict["EGF_EGFRp2"].s + species_dict["GAP"].i < r[8] > species_dict["EGF_EGFRp2_GAP"].s
            r[8].K = 1e6 * fac, 0.2  # 1/Ms 1e6, 0.2

            #i
            species_dict["EGF_EGFRp2_GAP_Grb2"].s + species_dict["Prot"].i < r[4] > species_dict["EGF_EGFRp2_GAP_Grb2_Prot"].s
            r[4].K = 1.73e-7 * fac, 1.66e-3 # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2_Prot"].s > r[5] > species_dict["EGF_EGFRp2_GAP_Grb2i"].i + species_dict["Proti"].i
            r[5].K = 0.03 # 1/Ms
            species_dict["EGFR"].s < r[6] > species_dict["EGFRi"].i
            r[6].K = 5e-5 , 5e-3  # 1/s
            species_dict["EGF_EGFR2"].s > r[7] > species_dict["EGF_EGFR2i"].i
            r[7].K =5e-5 # 1/s
            species_dict["EGF_EGFRp2_GAP"].s > r[9] > species_dict["EGF_EGFRp2_GAPi"].i
            r[9].K = 5e-5 # 1/s

            #Grb
            species_dict["EGF_EGFRp2_GAP"].s + species_dict["Grb2"].i < r[16] >  species_dict["EGF_EGFRp2_GAP_Grb2"].s
            r[16].K = 1e7 * fac, 0.055  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2"].s + species_dict["Sos"].i < r[17] > species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s
            r[17].K = 1e7 * fac, 0.06  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s + species_dict["Ras_GDP"].i < r[18] > species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GDP"].s
            r[18].K = 1.5e7 * fac, 1.3   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GDP"].s < r[19] > species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s + species_dict["Ras_GTP"].i
            r[19].K = 0.5, 1e5 * fac  # 1/Ms
            species_dict["Ras_GTP"].i + species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s < r[20] > species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GTP"].s
            r[20].K = 2.1e6 * fac, 0.4   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GTP"].s < r[21] > species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s + species_dict["Ras_GDP"].i
            r[21].K = 0.023 , 2.2e5 * fac  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Grb2_Sos"].s < r[34] > species_dict["EGF_EGFRp2_GAP"].s + species_dict["Grb2_Sos"].i
            r[34].K = 0.03, 4.5e6 * fac  # 1/Ms

            #Shc
            species_dict["EGF_EGFRp2_GAP"].s + species_dict["Shc"].i < r[22] > species_dict["EGF_EGFRp2_GAP_Shc"].s
            r[22].K = 2.1e7 * fac, 0.1   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shc"].s < r[23] > species_dict["EGF_EGFRp2_GAP_Shcp"].s
            r[23].K = 6 , 0.6   # 1/s
            species_dict["EGF_EGFRp2_GAP_Shcp"].s + species_dict["Grb2"].i < r[24] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2"].s
            r[24].K = 1e7 * fac, 0.55   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2"].s + species_dict["Sos"].i < r[25] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s
            r[25].K = 1e7 * fac, 0.0214  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s + species_dict["Ras_GDP"].i < r[26] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GDP"].s
            r[26].K = 1.5e7 * fac, 1.3   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GDP"].s < r[27] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s + species_dict["Ras_GTP"].i
            r[27].K = 0.5 , 1e5 * fac  # 1/Ms
            species_dict["Ras_GTPp"].i + species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s < r[30] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GTP"].s
            r[30].K = 7.9e6 * fac, 1.3   # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GTP"].s < r[31] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s + species_dict["Ras_GDP"].i
            r[31].K = 0.023 , 2.2e5 * fac  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s < r[32] > species_dict["EGF_EGFRp2_GAP"].s + species_dict["Shcp_Grb2_Sos"].i
            r[32].K = 0.1, 2.4e5 * fac  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp"].s < r [37] > species_dict["EGF_EGFRp2_GAP"].s + species_dict["Shcp"].i
            r[37].K = 0.3 , 9e5 * fac  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp_Grb2"].s < r[39] > species_dict["EGF_EGFRp2_GAP"].s + species_dict["Shcp_Grb2"].i
            r[39].K = 0.3 , 9e5 * fac  # 1/Ms
            species_dict["EGF_EGFRp2_GAP_Shcp"].s + species_dict["Grb2_Sos"].i < r[41] > species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"].s
            r[41].K = 3e7 * fac , 0.0429  # 1/Ms

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
        N = 1
        for s,r in result_selectors:
            print(f"Setting result selector {N}/{len(result_selectors)}: {s}")
            N += 1
            rs_path = rs.SUM(getattr(r, s).Count)
            sim.toSave(rs_path, dt = p["time step"])

    return sim, rs, mesh 