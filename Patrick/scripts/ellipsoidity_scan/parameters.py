import pandas as pd

p = {
    "time step" : 1,
    "endtime" : 1500,
    "DC" : 4e-12,
    "k[0]" : 1e8,
    "k[666]" : 0.3,
    "bc": 4.5,
    # initial conditions for the model, currently this has no effect as it is hardcoded into the SimManager
    # TODO: implement this relative path
    "initial_conditions_path": "Patrick/data_big_model_mini_sph.xls"
}

