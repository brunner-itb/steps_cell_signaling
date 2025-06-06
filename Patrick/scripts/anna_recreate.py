import steps.interface
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *
import numpy as np
import os
import time
import sys
import random

def create_model(p, p_name, factor, endt):

    mdl = Model()
    r = ReactionManager()

    # Diffusionskonstanten

    #DCR = p["DCR"]#2.5e-14 # m^2/s
    DC = p["DC"]#7e-12  # m^2/s


    #Faktor für Reaktionsraten
    c1 = 1
    c2 = 1

    with mdl:
        vsys = VolumeSystem.Create()
        exo_vsys = VolumeSystem.Create()
        vsys_nuc = VolumeSystem.Create()
        ssys = SurfaceSystem.Create()

        EGF, EGFR, EGF_EGFR, EGF_EGFR2, EGF_EGFRp2, GAP, EGF_EGFRp2_GAP, ERK, ERKp, ERKpp, P3  = Species.Create()

        with vsys:
            # ERK Deaktivierung
            ERKp + P3 > r[6] > ERK + P3
            r[6].K = p["k[666]"]
            #Cytoplasm diffusion
            Diffusion(ERK, DC)
            Diffusion(ERKp, DC)
            Diffusion(P3, 2*DC)
            Diffusion(GAP, DC/4)

        with exo_vsys:
            Diffusion(EGF, DC/10)

        with vsys_nuc:
            ERKp > r[7] > ERKpp
            r[7].K = 1e8
            Diffusion(ERKpp, DC)

        with ssys:
            EGFR.s + EGF.o < r[1] > EGF_EGFR.s
            EGF_EGFR.s + EGF_EGFR.s < r[2] > EGF_EGFR2.s
            EGF_EGFR2.s < r[3] > EGF_EGFRp2.s
            EGF_EGFRp2.s + GAP.i < r[4] > EGF_EGFRp2_GAP.s
            # None > r[13] > EGFR.s
            EGF_EGFRp2_GAP.s + ERK.i < r[5] > EGF_EGFRp2_GAP.s + ERKp.i
            r[1].K = 3e7, 38e-4   # 1/Ms
            r[2].K = 1e7 * 100, 0.1   # 1/Ms
            r[3].K = 1 * 100, 0.01  # 1/s
            r[4].K = 1e6 * 100 , 0.2   # 1/Ms 1e6, 0.2
            r[5].K = p["k[0]"], 0.1 #1e8 * c1, 0.1 * c1
            Diffusion(EGFR, DC/10)
            Diffusion(EGF_EGFR, DC/20)
            Diffusion(EGF_EGFR2, DC/40)
            Diffusion(EGF_EGFRp2, DC/40)
            Diffusion(EGF_EGFRp2_GAP, DC/50)

    # Mesh
    bc = p["bc"]
    ## Kugel mit Radius 12 und max.size 0.4
    mesh = TetMesh.LoadAbaqus(f'/home/pb/steps_cell_signaling/Patrick/meshes/elipsoid_4.5.inp', scale = 10**(-6))
    #Volume1 --> exo
    #Volume2 --> cyt
    #Volume3 --> nuc

    with mesh:

        # LISTEN

        # Exo
        exo_tets = TetList(mesh.tetGroups["Volume1"])

        # Zelle
        cyt_tets = TetList(mesh.tetGroups["Volume2"])

        # Zellkern
        nuc_tets = TetList(mesh.tetGroups["Volume3"])


        # COMPARTMENTS

        # Zellkern
        nuc = Compartment.Create(nuc_tets, vsys_nuc)

        # Cytoplasma
        cyt = Compartment.Create(cyt_tets, vsys)

        # Zelläüßeres
        exo = Compartment.Create(exo_tets, exo_vsys)

        # Zellmembran
        cell_surface = Patch.Create(cyt.surface & exo.surface, cyt, exo, ssys)

        # DIFFUSIONS BARRIERE

        # Zellkernmembran
        nuc_mem = DiffBoundary.Create(nuc.surface)

    ratio_v = exo.Vol / 1286e-18
    ratio_o = cell_surface.Area / 706e-12

    seed = random.randint(1, 6000)

    rng = RNG("mt19937", 512, seed)
    # sim = Simulation('Tetexact', mdl, mesh, rng)
    partition = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)
    sim = Simulation("TetOpSplit", mdl, mesh, rng, False, partition)
    rs = ResultSelector(sim)

    EGF_Count = rs.SUM(rs.TETS(cyt_tets).EGF.Count)
    EGF_EGFR_Count = rs.SUM(rs.TRIS(cell_surface.tris).EGF_EGFR.Count)
    EGF_EGFR2_Count = rs.SUM(rs.TRIS(cell_surface.tris).EGF_EGFR2.Count)
    EGF_EGFRp2_Count = rs.SUM(rs.TRIS(cell_surface.tris).EGF_EGFRp2.Count)
    EGF_EGFRp2_GAP_Count = rs.SUM(rs.TRIS(cell_surface.tris).EGF_EGFRp2_GAP.Count)
    ERKp_Count = rs.SUM(rs.TETS(cyt_tets).ERKp.Count)#rs.TETS(cell_tets).Xa.Count
    #ERKp_Origin = rs.TETS(mem_tet).Xa.Count
    ERKpp_Count = rs.SUM(rs.TETS(nuc_tets).ERKpp.Count)

    sim.toSave(EGF_Count, dt=p["time step"])
    sim.toSave(EGF_EGFR_Count, dt=p["time step"])
    sim.toSave(EGF_EGFR2_Count, dt=p["time step"])
    sim.toSave(EGF_EGFRp2_Count, dt=p["time step"])
    sim.toSave(EGF_EGFRp2_GAP_Count, dt=p["time step"])
    sim.toSave(ERKpp_Count, dt=p["time step"])


    for r in range(1):
        options = dict(compression="gzip",
                       compression_opts=5)  # compress the output files to save space. Larger opts = more compression
        with XDMFHandler("/home/pb/steps_cell_signaling/Patrick/saved_objects/testing", hdf5DatasetKwArgs=options) as hdf:
            sim.toDB(hdf, uid="test")
            sim.newRun()
            sim.exo.EGF.Count = p["EGF0"] * ratio_v
            print("EGF in exo: " + str(sim.exo.EGF.Count))
            sim.cell_surface.EGFR.Count = 5e4 * ratio_o
            print("EGFR on membrane: " + str(sim.cell_surface.EGFR.Count))
            sim.cyt.GAP.Count = 1.2e4
            print("GAP in cytoplasm: " + str(sim.cyt.GAP.Count))
            sim.cyt.ERK.Count = 2.1e4 #6.3e7
            sim.cyt.P3.Count = 1e3#1.5e5
            sim.nuc_mem.ERKp.DiffusionActive = True
            start = time.time()
            sim.run(endt)
            end = time.time()
        print("Durchlaufzeit: " + str(end-start))

    return


p = {"DC" : 4e-12, "k[666]" : 0.3, "k[0]" : 1e8, "time step" : 0.25, "EGF0": 1e4, "bc": 4.5}


#FOR PARAMETERS
#factor = sys.argv[2]
#factor = factor.replace(",", ".")
#factor = float(factor)

#parameter = sys.argv[3]
#endt = 15
#value = p[parameter] * factor

#create_model(p | {parameter: value}, parameter, factor, endt)

#FOR NORMAL
parameter = "tet_top"
endt = 1
print(f"Simulation of eli1, parameter: {parameter}, time: {endt}")

create_model(p, parameter, 1, endt)