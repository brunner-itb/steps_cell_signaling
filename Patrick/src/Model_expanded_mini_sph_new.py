import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
from src.Utilities import molar_to_molecules
import os
import time
import sys
import random
import numpy as np
import pandas as pd

def initialize_mesh(mesh_path, scale, nucleus_volume, volume_system, extracellular_volume, cell_surface):
    # Load mesh and compartments
    # Sphere with diameter 12 und max.size 0.4
    assert os.path.isfile(mesh_path), "mesh_path does not exist. Please check the path and try again."
    mesh = stgeom.TetMesh.LoadAbaqus(mesh_path, scale=scale)

    with mesh:
        #for Sphere
        #Volume1 --> exo
        #Volume2 --> cyt
        #Volume3 --> nuc
        #Exo /Extracellular volume
        exo_tets = stgeom.TetList(mesh.tetGroups["Volume1"])
        #Zelle /Cytoplasm
        cell_tets = stgeom.TetList(mesh.tetGroups["Volume2"])
        #Zellkern
        #nuc_tets = stgeom.TetList(mesh.tetGroups["Volume3"])

    # Create compartments
        # Zellkern
        #nuc = stgeom.Compartment.Create(nuc_tets, vsys_nuc)
        # Cytoplasm
        cyt = stgeom.Compartment.Create(cell_tets, volume_system)
        # Extracellular volume
        exo = stgeom.Compartment.Create(exo_tets, extracellular_volume)
        # Cellmembrane
        cell_surface = stgeom.Patch.Create(cyt.surface & exo.surface, cyt, exo, cell_surface)
    #Create Diffusion Barrier 
        #Nucleus membrane
        #nuc_mem = stgeom.DiffBoundary.Create(nuc.surface)
    return mesh, cell_tets, cyt, exo


def create_model(p, p_name, factor, endt, ii,species_names, mesh_path, plot_only_run):

    #Factors for Reaction Rates
    c1 = 1
    c2 = 1
    mdl =stmodel.Model()
    r = stmodel.ReactionManager()
    data_big_model_mini_sph_df = pd.read_excel(p["big_model_mini_sph_df_path"])
    species_dict = {}

    # Diffusionskonstanten
    #DCR = p["DCR"]#2.5e-14 # m^2/s
    #p = p["DC"]#7e-12  # m^2/s
    #Volumen in L
    V = 1.414e-17 * 1000 #dm^3
    fac = 1/V
    
    # Create volume and surface systems
    with mdl:
        volume_system = stmodel.VolumeSystem.Create()
        nucleus_volume = stmodel.VolumeSystem.Create()
        extracellular_volume = stmodel.VolumeSystem.Create()
        cell_surface = stmodel.SurfaceSystem.Create()

        ##vsys = VolumeSystem.Create() #volume_system
       ## exo_vsys = VolumeSystem.Create() #extracellular_volume
        #vsys_nuc = VolumeSystem.Create() #nucleus_volume
        ##ssys = SurfaceSystem.Create()  #cell_surface
 # todo : im parameter dict erstellen make two dicts for model und for simualation
        #EGF
       # EGF, EGFR, EGF_EGFR, EGF_EGFR2, EGF_EGFRp2, GAP, EGF_EGFRp2_GAP  = Species.Create()
        #Grb2
        #Grb2, EGF_EGFRp2_GAP_Grb2, Sos, EGF_EGFRp2_GAP_Grb2_Sos, EGF_EGFRp2_GAP_Grb2_Sos_Ras_GDP, Ras_GDP, Ras_GTPp, Ras_GTP, EGF_EGFRp2_GAP_Grb2_Sos_Ras_GTP, Grb2_Sos  = Species.Create()
        #Shc
        #Shc, EGF_EGFRp2_GAP_Shc, EGF_EGFRp2_GAP_Shcp, EGF_EGFRp2_GAP_Shcp_Grb2, EGF_EGFRp2_GAP_Shcp_Grb2_Sos, EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GDP, Shcp_Grb2_Sos, Shcp, EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GTP, Shcp_Grb2 = Species.Create()
        #Raf
        #Raf, Raf_Ras_GTP, Rafp, P1, Rafp_P1, MEK, MEK_Rafp, MEKp, MEKp_Rafp, MEKpp, P2, MEKpp_P2, MEKp_P2, ERK, ERK_MEKpp, ERKp, ERKp_MEKpp, ERKpp, ERKpp_P3, P3, ERKp_P3 = Species.Create()
        #i
       # Prot, Proti, EGF_EGFRp2_GAP_Grb2i, EGF_EGFRp2_GAP_Grb2_Prot, EGFRi, EGF_EGFRi, EGF_EGFR2i, EGF_EGFRp2i, EGF_EGFRp2_GAPi, EGFi = Species.Create()
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
            # None > r[13] > EGFR.s
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
            species_dict["ERK"] + species_dict["MEKpp"] < r[52] > species_dict["ERKp_MEKpp"]
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

            #Cytoplasm diffusion
            for index, row in data_big_model_mini_sph_df.iterrows():
                species_name = row['Species']
                cyt_dc = row['cyt DC']
                if pd.notna(cyt_dc):
                    stmodel.Diffusion(species_dict[species_name], cyt_dc)

            '''stmodel.Diffusion(species_dict["GAP"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Grb2"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Sos"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Ras_GTP"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Ras_GTPp"], p["DC"] / 4)
            stmodel.Diffusion(species_dict["Ras_GDP"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Shc"], p["DC"]/4)
            stmodel.Diffusion(species_dict["Shcp"], p["DC"]/4)
            stmodel.Diffusion(species_dict["ERK"], p["DC"])
            stmodel.Diffusion(species_dict["MEK"], p["DC"])
            stmodel.Diffusion(species_dict["ERKp"], p["DC"])
            stmodel.Diffusion(species_dict["ERKpp"], p["DC"])
            stmodel.Diffusion(species_dict["MEKp"], p["DC"])
            stmodel.Diffusion(species_dict["MEKpp"], p["DC"])
            stmodel.Diffusion(species_dict["P1"], 2 * p["DC"])
            stmodel.Diffusion(species_dict["P2"], 2 * p["DC"])
            stmodel.Diffusion(species_dict["P3"], 2 * p["DC"])
            stmodel.Diffusion(species_dict["Raf"], p["DC"] / 4)
            stmodel.Diffusion(species_dict["Rafp"], p["DC"] / 4)

            # i
            stmodel.Diffusion(species_dict["Prot"], p["DC"] / 10)
            stmodel.Diffusion(species_dict["Proti"], p["DC"] / 10)
            stmodel.Diffusion(species_dict["EGFRi"], p["DC"] / 10)
            stmodel.Diffusion(species_dict["EGF_EGFRi"], p["DC"] / 20)
            stmodel.Diffusion(species_dict["EGF_EGFR2i"], p["DC"] / 40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2i"], p["DC"] / 40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAPi"], p["DC"] / 50)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2i"], p["DC"] / 70)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2_Prot"], p["DC"] / 70)

            stmodel.Diffusion(species_dict["Grb2_Sos"], p["DC"]/8)
            stmodel.Diffusion(species_dict["Shcp_Grb2_Sos"], p["DC"]/16)
            stmodel.Diffusion(species_dict["Shcp_Grb2"], p["DC"]/8)
            stmodel.Diffusion(species_dict["ERKp_P3"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["ERKpp_P3"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["MEKp_P2"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["MEKpp_P2"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["ERK_MEKpp"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["ERKp_MEKpp"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["Rafp_P1"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["Raf_Ras_GTP"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["MEK_Rafp"], p["DC"] / 8)
            stmodel.Diffusion(species_dict["MEKp_Rafp"], p["DC"] / 8)'''


        with extracellular_volume:
            for index, row in data_big_model_mini_sph_df.iterrows():
                species_name = row['Species']
                exo_dc = row['exo Volume DC']
                if pd.notna(exo_dc):
                    stmodel.Diffusion(species_dict[species_name], exo_dc)
            #stmodel.Diffusion(species_dict["EGF"], p["DC"]/10)

        #with vsys_nuc: # means with nucles_volume:

        with cell_surface:

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
            '''stmodel.Diffusion(species_dict["EGF"], p["DC"]/10)
            stmodel.Diffusion(species_dict["EGF_EGFR"], p["DC"]/20)
            stmodel.Diffusion(species_dict["EGF_EGFR2"], p["DC"]/40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2"], p["DC"]/40)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP"], p["DC"]/50)


            #Grb2
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2"], p["DC"]/60)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2_Sos"], p["DC"]/70)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GTP"], p["DC"] / 70)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Grb2_Sos_Ras_GDP"], p["DC"] / 70)

            #Shc
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shc"], p["DC"] / 60)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shcp"], p["DC"] / 60)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shcp_Grb2"], p["DC"] / 70)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos"], p["DC"] / 80)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GDP"], p["DC"] / 80)
            stmodel.Diffusion(species_dict["EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GTP"], p["DC"] / 80)'''
    # -------Mesh--------------------------
    # if os.getcwd() != '/home/anna3171/annas-magnificent-steps-project':

    #    os.chdir(r"/home/anna3171/annas-magnificent-steps-project")

    ## Kugel mit Radius 12 und max.size 0.4
    #mesh = stgeom.TetMesh.LoadAbaqus('meshes/mini_sph.inp', scale = 10**(-6))

    #für Sphere
    #Volume1 --> exo
    #Volume2 --> cyt
    #Volume3 --> nuc

   # with mesh:

    # LISTEN

        #Exo /Extracellular volume
        #exo_tets = stgeom.TetList(mesh.tetGroups["Volume1"])
        #Zelle /Cytoplasm
       # cyt_tets = stgeom.TetList(mesh.tetGroups["Volume2"])
        #Zellkern
        #nuc_tets = TetList(mesh.tetGroups["Volume3"])


    # Create compartments

        # Zellkern
        #nuc = Compartment.Create(nuc_tets, vsys_nuc)

        # Cytoplasma
        #cyt = Compartment.Create(cyt_tets, vsys)

        # Zelläüßeres
        #exo = Compartment.Create(exo_tets, exo_vsys)

        # Zellmembran
        #cell_surface = Patch.Create(cyt.surface & exo.surface, cyt, exo, ssys)

    #DIFFUSIONS BARRIERE

        #Zellkernmembran
        #nuc_mem = DiffBoundary.Create(nuc.surface)

   
    #-----------Load mesh and compartments
    mesh, cell_tets = initialize_mesh(mesh_path,
                               scale=10 ** -6,
                               nucleus_volume=nucleus_volume,
                               volume_system=volume_system,
                               extracellular_volume=extracellular_volume,
                               cell_surface=cell_surface)
    system_volume = mesh.Vol

    ratio_v = mesh.exo.Vol / 1286e-18
    ratio_mini = mesh.cyt.Vol / 1766e-18 #added so all meshes have same volume

    # ----------Start Simulation Initilize RNG and Simulation
    seed = random.randint(1,6000)
    rng = strng.RNG("mt19937", 512, seed)
    partition = stgeom.LinearMeshPartition(mesh, 1, 1, stsim.MPI.nhosts)
    sim = stsim.Simulation("TetOpSplit", mdl, mesh, rng, False, partition) #simulation

    # Define results
    rs = None
    if plot_only_run == False:
        rs = stsave.ResultSelector(sim)

    # resultsselector
        result_selectors = []

        for index, row in data_big_model_mini_sph_df.iterrows():
            species_name = row['Species']
            result_selector = row['resultsselector']
            if result_selector == 'TRIS':
                result_selectors.append((species_name, rs.TRIS(cell_surface.tris)))
            elif result_selector == 'TETS':
                result_selectors.append((species_name, rs.TETS(mesh.cyt.tets)))


        '''species = [
            #TRIS
            ("EGFR", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFR", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFR2", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Grb2", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Grb2_Sos", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Grb2_Sos_Ras_GDP", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Grb2_Sos_Ras_GTP", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shc", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shcp", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shcp_Grb2", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shcp_Grb2_Sos", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GDP", rs.TRIS(cell_surface.tris)),
            ("EGF_EGFRp2_GAP_Shcp_Grb2_Sos_Ras_GTP", rs.TRIS(cell_surface.tris)),
            #TETS
            ("EGF", rs.TETS(mesh.exo.tets)),
            ("GAP", rs.TETS(mesh.cyt.tets)),
            ("Shc", rs.TETS(mesh.cyt.tets)),
            ("Grb2", rs.TETS(mesh.cyt.tets)),
            ("Sos", rs.TETS(mesh.cyt.tets)),
            ("Ras_GDP", rs.TETS(mesh.cyt.tets)),
            ("Grb2_Sos", rs.TETS(mesh.cyt.tets)),
            ("Shcp_Grb2_Sos", rs.TETS(mesh.cyt.tets)),
            ("Shcp_Grb2", rs.TETS(mesh.cyt.tets)),
            ("Shcp", rs.TETS(mesh.cyt.tets)),
            ("Ras_GTP", rs.TETS(mesh.cyt.tets)),
            ("Ras_GTPp", rs.TETS(mesh.cyt.tets)),
            ("Raf", rs.TETS(mesh.cyt.tets)),
            ("Raf_Ras_GTP", rs.TETS(mesh.cyt.tets)),
            ("Rafp", rs.TETS(mesh.cyt.tets)),
            ("P1", rs.TETS(mesh.cyt.tets)),
            ("Rafp_P1", rs.TETS(mesh.cyt.tets)),
            ("MEK", rs.TETS(mesh.cyt.tets)),
            ("MEK_Rafp", rs.TETS(mesh.cyt.tets)),
            ("MEKp_Rafp", rs.TETS(mesh.cyt.tets)),
            ("MEKp", rs.TETS(mesh.cyt.tets)),
            ("MEKpp", rs.TETS(mesh.cyt.tets)),
            ("P2", rs.TETS(mesh.cyt.tets)),
            ("MEKpp_P2", rs.TETS(mesh.cyt.tets)),
            ("MEKp_P2", rs.TETS(mesh.cyt.tets)),
            ("ERK", rs.TETS(mesh.cyt.tets)),
            ("ERK_MEKpp", rs.TETS(mesh.cyt.tets)),
            ("ERKp", rs.TETS(mesh.cyt.tets)),
            ("ERKp_MEKpp", rs.TETS(mesh.cyt.tets)),
            ("ERKpp", rs.TETS(mesh.cyt.tets)),
            ("ERKpp_P3", rs.TETS(mesh.cyt.tets)),
            ("P3", rs.TETS(mesh.cyt.tets)),
            ("ERKp_P3", rs.TETS(mesh.cyt.tets)),
            ("Prot", rs.TETS(mesh.cyt.tets)),
            ("Proti", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFRp2_GAP_Grb2i", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFRp2_GAP_Grb2_Prot", rs.TETS(mesh.cyt.tets)),
            ("EGFRi", rs.TETS(mesh.cyt.tets)),
            ("EGFi", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFRi", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFR2i", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFRp2i", rs.TETS(mesh.cyt.tets)),
            ("EGF_EGFRp2_GAPi", rs.TETS(mesh.cyt.tets))
        ]'''

        for s,r in result_selectors:

            rs_path = rs.SUM(getattr(r, s).Count)
            print(rs_path)
            sim.toSave(rs_path, dt = p["time step"]) #keep
        
        # i = ii#sys.argv[1] #get rid of
            #rs_path.toFile(f"saved objects/Einzelne Runs/{s}_mini_sph_{p_name}_{factor}_{i}.dat") #get rid of

        
        '''time_interval = np.arange(0,endt, 1)
        time_interval = np.delete(time_interval, 0)

        for r in range(1):
            # ToDo: pack initial rates also in Dataframe

            sim.newRun()
            sim.exo.EGF.Count = p["EGF0"] * ratio_mini
            #sim.TETS(tips_one).EGF.Count = p["EGF0"]/len(tips_one)
            print("EGF in exo: " + str(sim.exo.EGF.Count))
            sim.cell_surface.EGFR.Count = 5e4 * ratio_mini  #5e4
            print("EGFR on membrane: " + str(sim.cell_surface.EGFR.Count))
            sim.cyt.GAP.Count = 1.2e4 * ratio_mini #1.2e4
            print("GAP in cytoplasm: " + str(sim.cyt.GAP.Count))
            sim.cyt.Grb2.Count = 5.10e4 * ratio_mini #5.10e4
            print("Grb2 in cyt: " + str(sim.cyt.Grb2.Count))
            sim.cyt.Sos.Count = 6.63e4 * ratio_mini #6.63e4
            print("Sos in cyt: " + str(sim.cyt.Sos.Count))
            sim.cyt.Ras_GDP.Count = 1.14e7 * ratio_mini #1.14e7
            print("Ras_GDP in cyt: " + str(sim.cyt.Ras_GDP.Count))
            sim.cyt.Raf.Count = 4e4 * ratio_mini
            print("Raf in cyt: " + str(sim.cyt.Raf.Count))
            sim.cyt.Shc.Count = 1.01e6 * ratio_mini #1.01e6
            print("Shc in cyt: " + str(sim.cyt.Shc.Count))
            sim.cyt.MEK.Count = 2.20e7 * ratio_mini #2.20e7
            print("MEK in cyt: " + str(sim.cyt.MEK.Count))
            sim.cyt.ERK.Count = 2.10e7 * ratio_mini #2.10e7
            print("ERK in cyt: " + str(sim.cyt.ERK.Count))
            sim.cyt.P1.Count = 4e4 * ratio_mini #4e4
            print("P1 in cyt: " + str(sim.cyt.P1.Count))
            sim.cyt.P2.Count = 4e4 * ratio_mini #4e4
            print("P2 in cyt: " + str(sim.cyt.P2.Count))
            sim.cyt.P3.Count = 1e7 * ratio_mini #1e7
            print("P3 in cyt: " + str(sim.cyt.P3.Count))
            sim.cyt.Prot.Count = 8.10e4 * ratio_mini
            print("Prot in cyt: " + str(sim.cyt.Prot.Count))

    # everything from here into run file
            start = time.time()
            for t in time_interval:
                sim.run(t)
                print(f"simulated for {t}s")
                end_i = time.time()
                print("It took: " + str((end_i - start) / 60 / 60) + "h")
            end = time.time()
            print("Durchlaufzeit: " + str((end - start)/60/60) + "h")'''

        return sim, rs, mesh

#standard parameterssqueue
p = {"DC" : 4e-12, "time step" : 0.01, "EGF0": 1e4}

#FOR PARAMETERS
#factor = sys.argv[2]
#factor = factor.replace(",", ".")
#factor = float(factor)

#parameter = sys.argv[3]
#endt = 10

#value = p[parameter] * factor

#create_model(p | {parameter: value}, parameter, factor, endt)

#FOR NORMALcrea
#parameter = "expanded"
#endt = 30
#i = 10
#print(f"Simulation {i}")
#create_model(p, parameter, 1, endt, i)'''
