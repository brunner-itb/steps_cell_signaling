import numpy as np
import pandas as pd
import inspect

from steps.API_2.sim import SolverCallError

species_names = ["EGF", "EGFR", "EGF_EGFR", "EGF_EGFR2", "EGF_EGFRp2", "GAP", "EGF_EGFRp2_GAP", "ERK", "ERKp", "ERKpp",
                 "P3"]


df = pd.read_excel("/home/pb/steps_cell_signaling/Patrick/data_big_model_mini_sph.xls")

def dataframe_cleanup(df, columns=["Species"]):
    for col in columns:
        if col == "Species":
            df[col] = df[col].str.replace("'", "", regex=False)
            df[col] = df[col].str.replace(" ", "", regex=False)
    return df

df = dataframe_cleanup(df, ["Species"])


init_count_columns = df.filter(like='init count').columns
compartments_to_init_map = {col.split('init count')[0].strip(): col for col in init_count_columns}

# compartments matched to column
for compartment in compartments_to_init_map.keys():
    # print(df[compartments_to_init_map[c.name]])
    for s_idx, species in enumerate(self.species_names):
        initial_value = df.loc[df.Species == species][compartments_to_init_map[compartment]].values
        initial_value = 0 if np.isnan(initial_value) else initial_value
        try:
            getattr(getattr(self.simulation, compartment), species).Count = initial_value
            value = getattr(getattr(self.simulation, compartment), species).Count
            print(compartment, species, value)
        except steps.API_2.sim.SolverCallError:
            pass




# #%%
#
# import gmsh
# import sys
# from os import listdir
# from os.path import isfile, join
#
#
# input_path = "/home/pb/steps_cell_signaling/Patrick/bio_data_meshes/Blebbistatin/"
# onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]
#
#
# working_files = "/home/pb/steps_cell_signaling/Patrick/bio_data_meshes/" + "working_files.txt"
#
# for input_file in onlyfiles[:100]:
#     try:
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output
#
#
#         gmsh.merge(input_path + input_file)
#
#         # 1. Get all surfaces
#         surfaces = [e[1] for e in gmsh.model.getEntities(2)]
#         loop_tag = gmsh.model.geo.addSurfaceLoop(surfaces)
#         volume_tag = gmsh.model.geo.addVolume([loop_tag])
#
#
#         # 6. Generate mesh
#         gmsh.model.geo.synchronize()
#         gmsh.model.mesh.generate(3)
#
#         # 7. Write to output file
#         # gmsh.write(output_file)
#         with open(working_files, 'a') as output:
#             output.write(input_file + "\n")
#         gmsh.finalize()
#
#     except:
#         gmsh.finalize()
#         pass
# print(working_files)
#
# # #%%
# #
# # file = "0097_0080_accelerator_20210315_bakal01_erk_main_21-03-15_12-37-27.off"
# #
# # gmsh.initialize()
# # gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output
# #
# # gmsh.merge(input_path + file)
# #
# # # 1. Get all surfaces
# # surfaces = [e[1] for e in gmsh.model.getEntities(2)]
# # loop_tag = gmsh.model.geo.addSurfaceLoop(surfaces)
# # volume_tag = gmsh.model.geo.addVolume([loop_tag])
# #
# # # 6. Generate mesh
# # gmsh.model.geo.synchronize()
# # gmsh.model.mesh.generate(3)
# #
# # # 7. Write to output file
# # gmsh.write(output_file)
# # gmsh.finalize()
