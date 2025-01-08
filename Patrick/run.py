from SimManager import SimManager
from parameters import p

species_names = ["EGF", "EGFR", "Xa", "XA", "EGF_EGFR", "EGF_EGFR2",
                 "EGF_EGFRp2", "EGF_EGFRp2_GAP", "GAP"]

sm = SimManager(parameters=p,
                species_names=species_names,
                endt = 1,
                mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/kugel04.inp",
                save_file="saved_objects/initial_run/initial_run.h5")
sm.load_model()
sm.run(run_id=0)