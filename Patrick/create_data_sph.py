import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import steps.visual as stvis
from matplotlib import pyplot as plt
import numpy as np
import math
import os
import time
import sys
import pandas as pd

def create_model(p, p_name, run_id, endt):

    mdl = stmodel.Model()
    r = stmodel.ReactionManager()

    # Diffusionskonstanten

    DCR = p["DCR"]#2.5e-14 # m^2/s
    DCX = p["DCX"]#7e-12  # m^2/s


    #Faktor für Reaktionsraten
    c1 = 1
    c2 = 1

    with mdl:
        volume_system = stmodel.VolumeSystem.Create()
        extracellular_volume = stmodel.VolumeSystem.Create()
        nucleus_volume = stmodel.VolumeSystem.Create()
        cell_surface = stmodel.SurfaceSystem.Create()

        EGF, EGFR, EGF_EGFR, EGF_EGFR2, EGF_EGFRp2, GAP, EGF_EGFRp2_GAP, X, Xa, XA  = stmodel.Species.Create()

        with volume_system:
            # Xa Deaktivierung
            Xa > r[666] > X
            r[666].K = p["k[666]"]#6.67e-6
            stmodel.Diffusion(X, DCX)
            stmodel.Diffusion(Xa, DCX)

        with extracellular_volume:
            stmodel.Diffusion(EGF, DCX)
        with nucleus_volume:
            Xa > r[777] > XA
            r[777].K = 1e8
            stmodel.Diffusion(Xa, DCX)

        with cell_surface:
            EGFR.s + EGF.o < r[1] > EGF_EGFR.s
            EGF_EGFR.s + EGF_EGFR.s < r[2] > EGF_EGFR2.s
            EGF_EGFR2.s < r[3] > EGF_EGFRp2.s
            EGF_EGFRp2.s + GAP.i < r[5] > EGF_EGFRp2_GAP.s
            # None > r[13] > EGFR.s
            EGF_EGFRp2_GAP.s + X.i < r[0] > EGF_EGFRp2_GAP.s + Xa.i
            r[0].K = p["k[0]"], 0.1 #1e8 * c1, 0.1 * c1
            r[1].K = 3e7 * c1, 38e-4 * c1  # 1/Ms
            r[2].K = 1e7 * c2, 0.1 * c2  # 1/Ms
            r[3].K = 1, 0.01  # 1/s, eigentlich 1,0.01 aber r[2] ist laut dem Paper auch 0.01
            r[5].K = 1e6 * c1, 0.2 * c1  # 1/Ms 1e6, 0.2
            # r[13].K = 2.17 # Receptors/s
            stmodel.Diffusion(EGFR, DCR)
            stmodel.Diffusion(EGF_EGFR, DCR)
            stmodel.Diffusion(EGF_EGFR2, DCR)
            stmodel.Diffusion(EGF_EGFRp2, DCR)
            stmodel.Diffusion(EGF_EGFRp2_GAP, DCR)

    ## Kugel mit Radius 12 und max.size 0.4
    mesh = stgeom.TetMesh.LoadAbaqus('/home/pb/steps_cell_signaling/meshes/kugel04.inp', scale = 10**(-6))
    #Volume3 --> exo
    #Volume2 --> cell


    with mesh:

    # LISTEN

        #Zelle
        cell_tets = stgeom.TetList(mesh.tetGroups["Volume2"])

        #Zellkern
        nuc_tets = stgeom.TetList(
            (tet for tet in cell_tets if np.linalg.norm(tet.center) < 5e-6))

        #TETS im Cyt, die an die Membran grenzen
        mem_tris = cell_tets.surface
        mem_tet = stgeom.TetList()
        for tri in mem_tris:
           for tet in tri.tetNeighbs:
                mem_tet.append(tet)
        mem_tet = stgeom.TetList([tet for tet in mem_tet if tet in cell_tets])

    # COMPARTMENTS

        # Zellkern
        nuc = stgeom.Compartment.Create(nuc_tets, nucleus_volume)

        # Cytoplasma
        cyt = stgeom.Compartment.Create(cell_tets - nuc_tets, volume_system)

        # Zelläüßeres
        exo = stgeom.Compartment.Create(mesh.tetGroups["Volume3"], extracellular_volume)

        # Zellmembran
        cell_surface = stgeom.Patch.Create(cell_tets.surface, cyt, exo, cell_surface)

    #DIFFUSIONS BARRIERE

        #Zellkernmembran
        nuc_mem = stgeom.DiffBoundary.Create(nuc.surface & cyt.surface)


    rng = strng.RNG("mt19937", 512, 2903)
    sim = stsim.Simulation('Tetexact', mdl, mesh, rng)
    rs = stsave.ResultSelector(sim)

    EGF_EGFR_Count = rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFR.Count)
    EGF_EGFRp2_Count = rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2.Count)
    EGF_EGFRp2_GAP_Count = rs.SUM(rs.TRIS(cell_tets.surface).EGF_EGFRp2_GAP.Count)
    Xa_Count = rs.SUM(rs.TETS(cell_tets).Xa.Count)#rs.TETS(cell_tets).Xa.Count
    #Xa_Origin = rs.TETS(mem_tet).Xa.Count
    nuc_Count = rs.SUM(rs.TETS(nuc_tets).XA.Count)

    #meta data
    #Xa_Origin.metaData['distToCenter'] = [np.linalg.norm(tet.center) for tet in mem_tet]
    #Xa_Count.metaData['Vol'] = [tet.Vol for tet in cell_tets]
    #Xa_Count.metaData["center_x"] = [tet.center.x for tet in cell_tets]
    #Xa_Count.metaData["center_y"] = [tet.center.y for tet in cell_tets]
    #Xa_Count.metaData["center_z"] = [tet.center.z for tet in cell_tets]

    sim.toSave(EGF_EGFR_Count, dt = p["time step"])
    sim.toSave(EGF_EGFRp2_Count, dt = p["time step"])
    sim.toSave(EGF_EGFRp2_GAP_Count, dt = p["time step"])
    sim.toSave(Xa_Count, dt = p["time step"])
    #sim.toSave(Xa_Origin, dt = p["time step"])
    sim.toSave(nuc_Count, dt = p["time step"])

    if not os.path.isdir('saved_objects/initial_run/'):
        if not os.path.isdir("saved_objects"):
            os.mkdir("saved_objects")
        else:
            os.mkdir("saved_objects/initial_run/")
    EGF_EGFR_Count.toFile(f"saved_objects/initial_run/EGF-EGFR_sph_{p_name}_{run_id}.dat")
    EGF_EGFRp2_Count.toFile(f"saved_objects/initial_run/EGF-EGFRp2_sph_{p_name}_{run_id}.dat")
    EGF_EGFRp2_GAP_Count.toFile(f"saved_objects/initial_run/EGF-EGFRp2-GAP_sph_{p_name}_{run_id}.dat")
    Xa_Count.toFile(f"saved_objects/initial_run/Xa_Count_sph_{p_name}_{run_id}.dat")
    #Xa_Origin.toFile(f"saved_objects/initial_run/Xa_Origin_sph_{p_name}_{run_id}_{i}.dat")
    nuc_Count.toFile(f"saved_objects/initial_run/nuc_Count_sph_{p_name}_{run_id}.dat")

    for r in range(1):

        sim.newRun()
        sim.exo.EGF.Count = 4e4
        sim.cell_surface.EGFR.Count = 7.8e4
        sim.cyt.GAP.Count = 2.3e4
        sim.cyt.X.Count = 4.1e4 #4.1e7
        sim.nuc_mem.Xa.DiffusionActive = True
        start = time.time()
        sim.run(endt)
        end = time.time()
        print(end-start)

    return

#standard parameterssqueue
p = {"DCR" : 2.5e-14, "DCX" : 7e-12, "k[666]" : 0.02, "k[0]" : 1e8, "time step" : 0.01}

def run_sim(p_name, p_value, run_id, endt):

    #for p_value in values:
    create_model(p | {p_name: p_value},p_name, run_id, endt)

    return

run_id = 0
parameter = "k[666]"
endt = 15

value = p[parameter] * run_id

#run_sim("time step", values, 1) # (0.0001, 0.001, 0.01, 0.1, 1)
#run_sim( "DCR", value, run_id) #(2.5e-16, 2.5e-15, 2.5e-14, 2.5e-13, 2.5e-12))
run_sim(parameter, value, run_id, endt) #(7e-11, 7e-12, 7e-13))
#run_sim("k[0]", value, run_id) #(1e6, 1e7, 1e8, 1e9,1e10))
#run_sim("k[666]", value, run_id) #(6.67e-1, 6.67e-2, 6.67e-3, 6.67e-4, 6.67e-5, 6.67e-6, 6.67e-7, 6.67e-8))
