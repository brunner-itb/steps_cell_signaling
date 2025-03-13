import steps.visual as stvis
import matplotlib.colors as mpl_c
# import pyqtgraph as pg

def interactive_plots(sm):
    # potentially set the default background to a custom one (white = w). Couldn't find any other way to change it for single plots only.
    # Very unreliable, sometimes turns the whole plot that color. Avoid
    # pg.setConfigOption("background", "w")
    sc = stvis.SimControl(end_time=sm.parameters["endtime"], upd_interval=sm.parameters["time step"])
    with sc:
        # 3D displays
        # nucleus, cell_surface, extracellular, EGF_d, GAP_d, ERKp_cyto, ERKp_nuc, GAP_EGF, ERKp_cyto_ERKp_nuc = stvis.SimDisplay.Create(
        #     'nucleus', 'cell surface', 'extracellular space', 'EGF', 'GAP', 'ERKp cytosol', "ERKp nucleus", 'GAP EGF', 'ERKp cytosol and nucleus'
        # )

        # Static mesh elements
        # with extracellular:
        #     stvis.ElementDisplay(sm.result_selector.exo, color=mpl_c.to_rgba("grey", alpha=0.05))
        # with cell_surface:
        #     stvis.ElementDisplay(sm.result_selector.cell_surface, color=mpl_c.to_rgba("grey", alpha=0.05))
        # with nucleus:
        #     stvis.ElementDisplay(sm.result_selector.nuc, color=mpl_c.to_rgba("grey", alpha=0.05))

        # Dynamic elements, colors from https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=5
        # with EGF_d:
        #     stvis.ElementDisplay(sm.result_selector.exo.EGF, color=mpl_c.to_rgba("#e41a1c", alpha=1), spec_size=0.015) #red
        # with GAP_d:
        #     stvis.ElementDisplay(sm.result_selector.cyt.GAP, color=mpl_c.to_rgba("#377eb8", alpha=1), spec_size=0.015) #blue
        # with ERKp_cyto:
        #     stvis.ElementDisplay(sm.result_selector.cyt.ERKp, color=mpl_c.to_rgba("#984ea3", alpha=1), spec_size=0.015) #purple
        # with ERKp_nuc:
        #     stvis.ElementDisplay(sm.result_selector.nuc.ERKp, color=mpl_c.to_rgba("#ff7f00", alpha=1), spec_size=0.015) #orange

        GAP_EGF, ERKp_cyto_ERKp_nuc = stvis.SimDisplay.Create('GAP EGF', 'ERKp cytosol and nucleus')

        with ERKp_cyto_ERKp_nuc:
            stvis.ElementDisplay(sm.result_selector.exo, color=mpl_c.to_rgba("grey", alpha=0.05))
            stvis.ElementDisplay(sm.result_selector.cell_surface, color=mpl_c.to_rgba("grey", alpha=0.05))
            stvis.ElementDisplay(sm.result_selector.nuc, color=mpl_c.to_rgba("grey", alpha=0.05))

            stvis.ElementDisplay(sm.result_selector.cyt.ERKp, color=mpl_c.to_rgba("#984ea3", alpha=1),
                                 spec_size=0.015)  # purple
            stvis.ElementDisplay(sm.result_selector.nuc.ERKp, color=mpl_c.to_rgba("#ff7f00", alpha=1),
                                 spec_size=0.015)  # orange

        with GAP_EGF:
            stvis.ElementDisplay(sm.result_selector.exo, color=mpl_c.to_rgba("grey", alpha=0.05))
            stvis.ElementDisplay(sm.result_selector.cell_surface, color=mpl_c.to_rgba("grey", alpha=0.05))
            stvis.ElementDisplay(sm.result_selector.nuc, color=mpl_c.to_rgba("grey", alpha=0.05))

            stvis.ElementDisplay(sm.result_selector.exo.EGF, color=mpl_c.to_rgba("#e41a1c", alpha=1),
                                 spec_size=0.015) #red
            stvis.ElementDisplay(sm.result_selector.cyt.GAP, color=mpl_c.to_rgba("#377eb8", alpha=1),
                                 spec_size=0.015) #blue

    return sc