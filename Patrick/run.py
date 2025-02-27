from src.SimManager import SimManager
from parameters import p

species_names = ["EGF", "EGFR", "Xa", "XA", "EGF_EGFR", "EGF_EGFR2",
                 "EGF_EGFRp2", "EGF_EGFRp2_GAP", "GAP"]

sm = SimManager(parameters=p,
                species_names=species_names,
                mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/kugel04.inp",
                save_file ="saved_objects/initial_run/parallel_run",
                parallel = True,
                runname = "parallel_run")
# the result_selectors need to be different whether you want to plot the results on the go or save them. This adjusts
# it accordingly. But: when you do a "plot_only_run" it means there is no data saved, so be aware.
sm.load_model(type="small",
              plot_only_run=True)

sm.run(run_id=0)
