"""
An example with multiple decomp layers, connected with water flux modelled with CMF

requires: cmf
"""

import cmf
from decomp.cmfconnector import CmfConnector
import datetime
import numpy as np
import xarray as xr

def xyz(i):
    return i * 10, 0, cmf.boltzmann(i, 5, 1)


class MultiLayerDECOMPModel:

    cellcount = 1
    starttime = datetime.datetime(2018, 1, 1)
    dt = cmf.h
    depth = np.cumsum(np.arange(0.005, 30 * 0.005, 0.005))

    def __init__(self):
        p = cmf.project('N DOC')
        self.project = p
        self.decompcells = {}
        for i in range(self.cellcount):
            c, decompcell = self.cell_setup()
            self.decompcells[c] = decompcell


        cmf.connect_cells_with_flux(p.cells, cmf.Richards_lateral)
        self.outlet = self.make_boundaries()
        self.initialize()
        for c in p[:10]:
            self.decompcells[c].depose_litter(10000, 0)


    def cell_setup(self):
        p = self.project
        i = len(p)
        c = p.NewCell(*xyz(i), area=100, with_surfacewater=True)
        vgm = cmf.VanGenuchtenMualem()
        vgm.fit_w0()
        for d in self.depth:
            vgm.Ksat = 10 * 2 ** -d
            c.add_layer(d, vgm)
        c.install_connection(cmf.Richards)
        c.install_connection(cmf.GreenAmptInfiltration)
        if i:
            c.topology.AddNeighbor(p[-2], 10)

        decompcell = CmfConnector(c, 5)
        return c, decompcell

    def make_boundaries(self):
        p = self.project
        c = p[0]
        # outlet = p.NewOutlet('Ditch', *xyz(-1))
        outlet = p.NewOutlet('GW', c.x , c.y, c.z - c.soildepth)
        #outlet.potential -= 0.5
        #p[0].connect_soil_with_node(outlet, cmf.Richards_lateral, 10, 10)
        cmf.FreeDrainagePercolation(c.layers[-1], outlet)
        rainfall = cmf.timeseries.from_sequence(self.starttime, cmf.day, [25, 0, 0, 0, 0, 0, 0] * 200)
        p.rainfall_stations.add('Heavy rain once a week', rainfall, (0, 0, 0))
        print(cmf.describe(p.rainfall_stations))
        p.use_nearest_rainfall()

        return outlet

    def initialize(self):
        for c in self.project:
            c.saturated_depth = 0.5
            self.decompcells[c].Cpool = 2000 * 2**(-c.layers.lower_boundary * 10)
            self.decompcells[c].depose_root(10 * 2**-c.layers.lower_boundary)

    def make_integrator(self):
        winteg = cmf.CVodeIntegrator(1e-9)
        sinteg = cmf.HeunIntegrator()
        integ = cmf.SoluteWaterIntegrator(self.project.solutes, winteg, sinteg, self.project)
        return integ


    def make_result_table(self, for_time, *columns):


        time = [t.AsPython()
                for t in cmf.timerange(cmf.AsCMFtime(self.starttime),
                                       cmf.AsCMFtime(self.starttime + for_time),
                                       self.dt)]
        depth = self.depth
        layer_template = np.zeros((len(time), len(self.project), len(depth))) * np.NaN
        ds = xr.Dataset(
            {'N': (('time', 'cell', 'depth'), layer_template.copy()),
             'Corg': (('time', 'cell', 'depth'), layer_template.copy()),
             'DOC': (('time', 'cell', 'depth'), layer_template.copy()),
             'Temp': (('time', 'cell', 'depth'), layer_template.copy()),
             'wetness': (('time', 'cell', 'depth'), layer_template.copy()),
             },
            {'time': time, 'cell': [c.x for c in self.project], 'depth': list(depth)}
        )
        return ds

    def run(self, for_time):
        p = self.project
        N, DOC = p.solutes
        integ = self.make_integrator()

        result = self.make_result_table(for_time, 'outlet_conc_N', 'outlet_conc_')
        for i, t in enumerate(integ.run(self.starttime, self.starttime + for_time, self.dt)):
            T = p[0].get_weather(t).T
            for c in p:
                self.decompcells[c].run(T, self.dt / cmf.day)
                result.Temp[i, c.Id] = self.decompcells[c].T_profile
                result.N[i, c.Id] = [l[N].conc for l in c.layers]
                result.DOC[i, c.Id] = [l[DOC].conc for l in c.layers]
                result.wetness[i, c.Id] = [l.wetness for l in c.layers]
                result.Corg[i, c.Id] = self.decompcells[c].Cpool
            print("%15s dt=%15s N=%g DOC=%g q=%g" %
                  (t, integ.dt, self.outlet.conc(t, N),
                   self.outlet.conc(t, DOC), self.outlet.waterbalance(t)))
        return result

    def plot_results(self, a: xr.Dataset):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(nrows= len(a.data_vars), sharex=True, sharey=True)
        for ax, var in zip(axes, a.data_vars):
            data = a[var]
            plt.sca(ax)
            data.plot(x='time', cmap=plt.cm.viridis_r, yincrease=False, robust=True)
        plt.show()



if __name__ == '__main__':
    model = MultiLayerDECOMPModel()
    a = model.run(for_time=datetime.timedelta(days=365))
    model.plot_results(a)


