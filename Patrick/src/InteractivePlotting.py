import steps.visual as stvis
import matplotlib.colors as mpl_c
# import pyqtgraph as pg

def interactive_plots(sm):
    # potentially set the default background to a custom one (white = w). Couldn't find any other way to change it for single plots only.
    # Very unreliable, sometimes turns the whole plot that color. Avoid
    # pg.setConfigOption("background", "w")
    sc = stvis.SimControl(end_time=sm.parameters["endtime"], upd_interval=sm.parameters["time step"])
    with sc:
        # # Plots
        # plots_d = stvis.PlotDisplay('EGF exo')
        # with plots_d:
        #     stvis.TimePlot(sm.result_selector.exo.EGF.Conc,
        #                    title='EGF_exo',
        #                    pen=(255, 0.647 * 255, 0),
        #                    data_size=1000,
        #                    y_range=[0, 15e-6],
        #                    y_label=('Concentration', 'M')
        #                    )

        # 3D displays
        EGF_d, cyt_d, full_d = stvis.SimDisplay.Create('EGF', 'cytosol', 'Full view')

        with EGF_d:
            # Static mesh element
            # stvis.ElementDisplay(sm.result_selector.exo, color=mpl_c.to_rgba("grey", alpha=0.05))
            # Dynamic element
            stvis.ElementDisplay(sm.result_selector.exo.EGF, color=mpl_c.to_rgba("red", alpha=1), spec_size=0.01)
        with cyt_d:
            stvis.ElementDisplay(sm.result_selector.cell_surface, color=mpl_c.to_rgba("grey", alpha=0.05))

        full_d.merge(EGF_d, cyt_d)
    return sc