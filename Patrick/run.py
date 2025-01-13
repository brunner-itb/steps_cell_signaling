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
# sm.run(run_id=0)

import steps.interface
import steps.model as stmodel
import steps.geom as stgeom
import steps.rng as strng
import steps.sim as stsim
import steps.saving as stsave
import steps.visual as stvis
sc = stvis.SimControl(end_time = 0, upd_interval = 0.0001)