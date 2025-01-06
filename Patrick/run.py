from SimManager import SimManager
from parameters import p

species_names = ["EGF", "EGFR", "Xa", "XA", "EGF_EGFR", "EGF_EGFR2",
                 "EGF_EGFRp2", "EGF_EGFRp2_GAP", "GAP"]

sm = SimManager(parameters=p, species_names=species_names, end_time = 1, mesh_path = "/home/pb/steps_cell_signaling/Patrick/meshes/kugel04.inp")