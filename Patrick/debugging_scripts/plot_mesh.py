import steps.interface
from steps.model import *
from steps.geom import *
from steps.rng import *
from steps.sim import *
from steps.saving import *
from steps.visual import *
import os
import numpy as np
import pandas as pd

#%%

mesh_path = f'/home/pb/steps_cell_signaling/Patrick/meshes/kugel_7.5.inp'
assert os.path.exists(mesh_path)

mdl = Model()
#mesh = TetMesh.LoadAbaqus('meshes/finger_0.45.inp', scale = 1.01e-6)
mesh = TetMesh.LoadAbaqus(mesh_path, scale = 10**(-6))
#mesh = TetMesh.LoadAbaqus(f'meshes/kugel_7.5-ext.inp', scale = 10**(-6))
#type = "Elipsoid"
#mesh = TetMesh.LoadAbaqus("meshes/hexa2.inp", scale = 1e-6)

with mdl:
    vsys = VolumeSystem.Create()
    exo_vsys = VolumeSystem.Create()
    vsys_nuc = VolumeSystem.Create()
    ssys = SurfaceSystem.Create()

with mesh:

    # EXO
    exo_tets = TetList(mesh.tetGroups["Volume1"])

    # Zelle
    cyt_tets = TetList(mesh.tetGroups["Volume2"])

    # Zellkern
    nuc_tets = TetList(mesh.tetGroups["Volume3"])

    # Zellkern
    nuc = Compartment.Create(nuc_tets, vsys_nuc)

    # Cytoplasma
    cyt = Compartment.Create(cyt_tets, vsys)

    # Zelläüßeres

    exo = Compartment.Create(exo_tets, exo_vsys)

    cell_surface = Patch.Create(cyt.surface & exo.surface, cyt, exo, ssys)

    nuc_mem = DiffBoundary.Create(nuc.surface)

# %%
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plotTriangles(ax, tris, color):
    ax.add_collection(Poly3DCollection(
        [tri.verts for tri in tris],
        facecolor=color,
        edgecolors='black',
        linewidth=0.1
    ))

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')

plotTriangles(ax, cyt.surface, (0.1, 0.2, 0.5, 0.09))
plotTriangles(ax, exo.surface, (0.1, 0.2, 0.5, 0.09))

ax.view_init(elev=0, azim=0, roll=0)
ax.set_xlim(mesh.bbox.min.x, mesh.bbox.max.x)
ax.set_ylim(mesh.bbox.min.y, mesh.bbox.max.y)
ax.set_zlim(mesh.bbox.min.z, mesh.bbox.max.z)
ax.set_xlabel('x position [m]')
ax.set_ylabel('y position [m]')
ax.set_zlabel('z position [m]')
ax.set_aspect('equal')
#plt.savefig(f"Plots/surface/3D_hexa2.pdf")
plt.show()