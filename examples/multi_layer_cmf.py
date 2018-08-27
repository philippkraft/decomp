"""
An example with multiple decomp layers, connected with water flux modelled with CMF

requires: cmf
"""

import cmf
from decomp.cmfconnector import CmfConnector
import datetime
import numpy as np
import xarray as xr


class MultiLayerDECOMPModel:

    starttime = datetime.datetime(2018, 1, 1)
    dt = cmf.h
    depth = np.cumsum(np.arange(0.005, 30 * 0.005, 0.005))

    def __init__(self):
        p = cmf.project('N DOC')
        self.project = p

        self.decompcell = self.cell_setup()


        self.outlet = self.make_boundaries()
        self.initialize()
        self.decompcell.depose_litter(10000, 0)


    def cell_setup(self):
        """
        Creates a cmf cell with all its layers and DECOMP pools
        :return: decomp.cmfconnector.CmfConnector
        """
        p = self.project
        c = p.NewCell(0, area=1000, with_surfacewater=True)
        vgm = cmf.VanGenuchtenMualem()
        vgm.fit_w0()
        for d in self.depth:
            vgm.Ksat = 10 * 2 ** -d
            c.add_layer(d, vgm)
        c.install_connection(cmf.Richards)
        c.install_connection(cmf.GreenAmptInfiltration)
        return CmfConnector(c, 5)

    def make_boundaries(self):
        """
        Creates the boundary conditions like outlet and rain
        :return:
        """
        p = self.project
        c = p[0]
        outlet = p.NewOutlet('GW', c.x, c.y, c.z - c.soildepth)
        cmf.FreeDrainagePercolation(c.layers[-1], outlet)
        rainfall = cmf.timeseries.from_sequence(self.starttime, cmf.day, [25, 0, 0, 0, 0, 0, 0] * 200)
        p.rainfall_stations.add('Heavy rain once a week', rainfall, (0, 0, 0))
        print(cmf.describe(p.rainfall_stations))
        p.use_nearest_rainfall()

        return outlet

    def initialize(self):
        """
        Set initial conditions
        :return:
        """
        c = self.decompcell.cmf_cell
        c.saturated_depth = 0.5
        self.decompcell.Cpool = 2000 * 2**(-c.layers.lower_boundary * 10)
        self.decompcell.depose_root(10 * 2**-c.layers.lower_boundary)

    def make_integrator(self):
        """
        Create the solver
        :return:
        """
        winteg = cmf.CVodeIntegrator(1e-9)
        sinteg = cmf.HeunIntegrator()
        integ = cmf.SoluteWaterIntegrator(self.project.solutes, winteg, sinteg, self.project)
        return integ


    def make_result_table(self, for_time):
        """
        Creates a xarray.Dataset to store the results
        :param for_time: duration of model run

        :return:
        """
        time = [t.AsPython()
                for t in cmf.timerange(cmf.AsCMFtime(self.starttime),
                                       cmf.AsCMFtime(self.starttime + for_time),
                                       self.dt)]
        depth = self.depth
        layer_template = np.zeros((len(time), len(depth))) * np.NaN
        ds = xr.Dataset(
            {'N': (('time', 'depth'), layer_template.copy()),
             'Corg': (('time', 'depth'), layer_template.copy()),
             'DOC': (('time', 'depth'), layer_template.copy()),
             'Temp': (('time', 'depth'), layer_template.copy()),
             'wetness': (('time', 'depth'), layer_template.copy()),
             },
            {'time': time, 'depth': list(depth)}
        )
        return ds

    def run(self, for_time):
        p = self.project
        N, DOC = p.solutes
        integ = self.make_integrator()

        result = self.make_result_table(for_time)
        for i, t in enumerate(integ.run(self.starttime, self.starttime + for_time, self.dt)):
            T = p[0].get_weather(t).T
            c = self.decompcell.cmf_cell
            self.decompcell.run(T, self.dt / cmf.day)
            result.Temp[i] = self.decompcell.T_profile
            result.N[i] = [l[N].conc for l in c.layers]
            result.DOC[i] = [l[DOC].conc for l in c.layers]
            result.wetness[i] = [l.wetness for l in c.layers]
            result.Corg[i] = self.decompcell.Cpool
            print("%15s dt=%15s N=%g DOC=%g q=%g" %
                  (t, integ.dt, self.outlet.conc(t, N),
                   self.outlet.conc(t, DOC), self.outlet.waterbalance(t)))
        return result

    def plot_results(self, a):
        """
        Plots the result dataset
        :param a: The result xarray.Dataset
        :return:
        """
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(nrows=len(a.data_vars), sharex='all', sharey='all')
        for ax, var in zip(axes, a.data_vars):
            data = a[var]
            plt.sca(ax)
            data.plot(x='time', cmap=plt.cm.viridis_r, yincrease=False, robust=True)
        plt.show()



if __name__ == '__main__':
    model = MultiLayerDECOMPModel()
    a = model.run(for_time=datetime.timedelta(days=365))
    model.plot_results(a)


