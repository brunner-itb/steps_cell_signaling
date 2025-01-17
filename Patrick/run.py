from src.SimManager import SimManager
from parameters import p

species_names = ["EGF", "EGFR", "Xa", "XA", "EGF_EGFR", "EGF_EGFR2",
                 "EGF_EGFRp2", "EGF_EGFRp2_GAP", "GAP"]

sm = SimManager(parameters=p,
                species_names=species_names,
                endt = 1,
                mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/kugel04.inp",
                save_file="saved_objects/initial_run/initial_run.h5")
sm.load_model(type="small")

import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import steps.visual as stvis
sc = stvis.SimControl(end_time = 1, upd_interval = 0.01)

with sc:
    # Plots
    simWin = stvis.SimDisplay.Create('My title')
    with simWin:
        stvis.ElementDisplay(sm.result_selector.cyt, color=(int(1), int(1), int(1), int(1)))

sm.run(run_id=0)
