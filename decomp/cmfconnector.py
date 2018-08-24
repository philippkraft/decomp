# -*- coding: utf-8 -*-
""" Module for interfacing decomp++ with cmf (Versions of late Oct. 2009) """
from __future__ import division, print_function, absolute_import, unicode_literals
import decomp
import numpy as np


class CmfConnector(object):
    """Class for creating decomp++ instances for each layer in a cmf cell
    """

    def __init__(self, cmf_cell, T_avg, max_Corg_depth=1e308):
        """
        Creates the interface for a cell from cmf
        max_Corg_depth [m] can be used to limit the number of layers
        owning decomp instances. Only layers whose upper boundary is
        less than max_Corg_depth get decomp models

        :param cmf_cell: A cmf cell with layers
        :param T_avg: The yearly average temperature in deg C
        :param max_Corg_depth: The lower boundary of Corg
        """

        c = cmf_cell
        self.cmf_cell = cmf_cell
        self.__decomplayers = [decomp.SOM() for l in c.layers
                               if l.upper_boundary < max_Corg_depth]
        self.T_profile = np.ones(c.layer_count()) * T_avg
        self.T_depth = 2.0
        self.pH = 7.0

    def depose_litter(self, leave_mass, wood_mass):
        """Deposes leaves and wood at the first layer
        leave_mass = Fallen leaves in g/m2
        wood_mass = Fallen wood in g/m2
        """
        self.__decomplayers[0] += leave_mass * decomp.leave_litter()
        self.__decomplayers[0] += wood_mass * decomp.wood_litter()

    def depose_root(self, root_mass):
        for i in range(len(self.__decomplayers)):
            self.__decomplayers[i] += root_mass[i] * decomp.root_litter()

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
            fieldcapacity = l.soil.Wetness_pF([1.8])[0]
            # set wetness of decomp layer
            wetness = min(1, l.wetness / fieldcapacity)
            # get T damping factor
            fT = 365**(-l.upper_boundary / self.T_depth)
            # set Temperature of layer
            self.T_profile[i] = fT * T + (1 - fT) * self.T_profile[i]
            # set DOC input
            # DOC precipitation currently disabled
            self.__decomplayers[i][decomp.DOC] = l[DOC].state
            # Integrate the decomp
            decomp_rate = self.__decomplayers[i].integrate(dt, self.T_profile[i], wetness, self.pH)
            # Update cmf
            l[N].source = decomp_rate.N
            l[DOC].source = decomp_rate[decomp.DOC]


