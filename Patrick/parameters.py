import pandas as pd

big_model_mini_sph_file_path = "steps_cell_signaling/Patrick/src/data_big_model_mini_sph.xls"
data_big_model_mini_sph_df = pd.read_excel(big_model_mini_sph_file_path)

p = {
    "time step" : 0.01,
    "endtime" : 5,
    "DC" : 4e-12,
    "k[0]" : 1e8,
    "k[666]" : 0.3,
    "bc": 4.5,
    # initial conditions small model
    "EGF_0": 1e4,
    "EGFR_0": 5e4,
    "GAP_0": 2.3e4,
    "ERK_0": 2.1e4,
    "P3_0": 1e3,
    # initial conditions large model
    "big_model_mini_sph_df_path": "steps_cell_signaling/Patrick/src/data_big_model_mini_sph.xls"
}

