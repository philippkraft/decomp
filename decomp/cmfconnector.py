# -*- coding: utf-8 -*-
""" Module for interfacing decomp++ with cmf (Versions of late Oct. 2009) """
from __future__ import division, print_function, absolute_import, unicode_literals
import decomp


class cmfconnector(object):
    """Class for creating decomp++ instances for each layer in a cmf cell
    """

    def __init__(self, cmf_cell, max_Corg_depth=1000):
        """Creates the interface for a cell from cmf
        max_Corg_depth [m] can be used to limit the number of layers
        owning decomp instances. Only layers whose upper boundary is
        less than max_Corg_depth get decomp models
        """
        c = cmf_cell
        self.cmf_cell = cmf_cell
        self.__decomplayers = [decomp.SOM() for l in c.layers
                               if l.upper_boundary < max_Corg_depth]
        self.T_depth = 2.0
        self.pH = 7.0

    def depose_litter(self, leave_mass, wood_mass):
        """Deposes leaves and wood at the first layer
        leave_mass = Fallen leaves in g/m2
        wood_mass = Fallen wood in g/m2
        """
        self.__decomplayers[0] += leave_mass * decomp.leave_litter()
        self.__decomplayers[0] += wood_mass * decomp.wood_litter()

    def plow(self, plowdepth=0.3):
        """Homogenizes the Corg content in all layers where the upper boundary
        is smaller than the plow depth
        """
        plowlayers = [
            l for l in self.cmf_cell.layers if l.upper_boundary < plowdepth - 0.01]
        sumSOM = decomp.SOM()
        for l in plowlayers:
            sumSOM += self[l]
        sumdepth = sum(l.thickness for l in plowlayers)
        for l in plowlayers:
            self[l] = sumSOM * (l.thickness / sumdepth)

    def __getitem__(self, index):
        if hasattr(index, "Position"):
            return self.__decomplayers[index.Position]
        else:
            return self.__decomplayers[index]

    def __setitem__(self, index, SOM):
        if hasattr(index, "Position"):
            self.__decomplayers[index.Position] = SOM
        else:
            self.__decomplayers[index] = SOM

    def __iter__(self):
        return iter(self.__decomplayers)

    def __getCpool(self):
        """Returns the mass of carbon stored
        """
        return [l.C for l in self.__decomplayers]

    def __setCpool(self, value):
        for i, l in enumerate(self.__decomplayers):
            if (i < len(value)):
                l[decomp.RC] = value[i]
                l.N = value[i] / 20.
            else:
                l[decomp.RC] = 0.0
    Cpool = property(__getCpool, __setCpool, "Mass of carbon per m²")


    def run(self, T, dt=1 / 24):
        """Runs the decomp model for time step dt (a float in days)
        """
        N, DOC = self.cmf_cell.project.solutes
        for i, l in enumerate(self.cmf_cell.layers):
            if i + 1 > len(self.__decomplayers):
                break
            # field capacity
            fieldcapacity = l.soil.Wetness_pF(1.8)
            # set wetness of decomp layer
            wetness = min(1, l.wetness / fieldcapacity)
            # get T damping factor
            fT = 365**(-l.upper_boundary / self.T_depth)
            # set Temperature of layer
            layerT = fT * T + (1 - fT) * self.__decomplayers[i].T
            # set DOC input
            # DOC precipitation currently disabled
            # self.__decomplayers[i][decomp.DOC] = l.Solute(DOC).state
            # Integrate the decomp
            decomp_rate = self.__decomplayers[i].integrate(dt, layerT, wetness, self.pH)
            # Update cmf
            l.Solute(N).source = decomp_rate.N
            l.Solute(DOC).source = decomp_rate[decomp.DOC]


if __name__ == "__main__":
    import cmf
    from numpy import *
    # project
    p = cmf.project('N DOC')
    N, DOC = p.solutes
    bc = cmf.BrooksCoreyRetentionCurve()
    bc.SetKsat(2.0, 0.5)

    last_cell = None
    depth = arange(0.1, 3.1, 0.1)
    for i in range(21):
        c = p.NewCell(i * 10 - 100, 0, 0, 100)
        if last_cell:
            c.topology.AddNeighbor(last_cell, 10)
        last_cell = c
        for d in depth:
            c.add_layer(d, bc)
        c.install_connection(cmf.Richards)
    decompcells = dict((c, cmfconnector(c)) for c in p)
    cmf.connect_cells_with_flux(p.cells, cmf.Richards_lateral)
    outlet = p.NewOutlet('Ditch', 2.0, 0.0, -1.0)
    p[10].connect_soil_with_node(
        outlet, cmf.Richards_lateral, 5, 10)  # ,0.0,1.)
    for c in p:
        c.saturated_depth = 0.5
        decompcells[c].Cpool = 2000 * 2**(-depth * 10)
        decompcells[c].set_root_litter(10 * 2**-depth)
    for c in p[:10]:
        decompcells[c].depose_litter(10000, 0)
    winteg = cmf.CVodeIntegrator(1e-6)
    winteg.preconditioner = 'R'
    sinteg = cmf.ImplicitEuler(1e-6)
    integ = cmf.SoluteWaterIntegrator(winteg, sinteg, p)
    rain = cmf.timeseries(integ.t, cmf.day)
    # rain.add(10)
    rain.extend([25, 0, 0, 0, 0, 0, 0] * 200)
    cmf.set_precipitation(p.cells, rain)

    def run(for_time=cmf.year):
        for t in integ.run(integ.t, integ.t + for_time, cmf.day):
            T = p[0].get_weather(t).T
            for c in p:
                decompcells[c].run(T, 1)
            print("%15s dt=%15s N=%g DOC=%g q=%g" % (t, integ.dt, outlet.conc(
                t, N), outlet.conc(t, DOC), outlet.water_balance(t)))
    from pylab import *

    def cN(): return transpose([[l.conc(N) for l in c.layers] for c in p])

    def sN(): return transpose(
        [[l.Solute(N).source for l in c.layers] for c in p])

    def dry(): return transpose([[1 - l.wetness for l in c.layers] for c in p])

    def pot(): return transpose([[l.potential for l in c.layers] for c in p])
    layers = cmf.node_list.from_sequence(cmf.get_layers(p))

    def cmfshow(a, layers, t, scale, **kwargs):
        clf()
        imshow(a, extent=(-105, 105, -3, 0), aspect='auto', **kwargs)
        colorbar()
        pos = layers.get_positions()
        f = layers.get_fluxes3d(t)
        quiver(array(pos.X), array(pos.Z), array(
            f.X), array(f.Z), scale=scale, hold=1)
        plot([outlet.Location.x], [outlet.Location.z], 'k*', ms=12, hold=1)
        annotate("outlet", [outlet.Location.x, outlet.Location.z], (30, 30),
                 textcoords='offset points', arrowprops=dict(width=2))
        axis('tight')
