|github| |appveyor| |zenodo|

.. |github| image:: https://badge.fury.io/gh/philippkraft%2Fdecomp.svg
    :target: https://github.com/philippkraft/decomp
.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/ju30o7wbxp2p0udo?svg=true
    :target: https://ci.appveyor.com/project/philippkraft/decomp
.. |zenodo| image:: https://zenodo.org/badge/145987487.svg
   :target: https://zenodo.org/badge/latestdoi/145987487

DECOMP++
========

Model history
----------------

DECOMP is a model to calculate the mineralization rate from soil organic matter and
was developed in its original form by the Department of Chemical Engineering,
Lund University, Sweden.
The parameters and equations are based on a review by Walse, Berg and Sverdrup 
(1998 https://doi.org/10.1139/a98-001). This basis has been implemented into a
FORTRAN computer model by Wallman, Belyazid, Svensson and Sverdrup
(2006, https://doi.org/10.1016/j.envsoft.2004.09.026),
which is also a part of the forest soil chemistry and plant growth model ForSAFE by the same 
authors (2005, https://doi.org/10.1016/j.foreco.2004.10.016). The single DECOMP model and 
the version in ForSAFE differ in a way, that the ForSAFE version includes a mechanism for
Nitrogen immobilisation and release based on the N/C ratio.

This version
-------------

Philipp Kraft from the Justus Liebig University, Gießen, Germany has translated the FORTRAN code of
ForSAFE-DECOMP into C++ and changed the Nitrogen dynamics to the immobilisation / release scheme
with the formulation from the VSD model (eq. 17 in Posch & Reinds 2003, https://doi.org/10.1016/j.envsoft.2008.09.007)

The C++ version, now called DECOMP++ comes with a wrapper for the programming language Python 
and will in most cases used from there. A virtual experiment using the Python wrapper and
CMF as a transport model can be found here (Kraft et al. 2011, https://doi.org/10.5194/adgeo-27-51-2010)

The carbon model
----------------

DECOMP uses a partition of the soil organic matter (SOM) into four components:

- EDC:easily decomposable compounds
- CELL: cellulose like compounds
- LIGN: lignin like compounds
- RC: recalcitrant compounds

Each compound has its own parameters for decomposition kinetic, as a function of temperature, soil moisture
and soil solution pH. Another property of the compounds is the decomposition product: Large parts of the SOM
compounds a released to CO2 and leave the system, smaller parts a released as DOC or are transformed to another
compound. The following figure from Wallman et al. 2006 shows the composition:

.. image:: https://ars.els-cdn.com/content/image/1-s2.0-S1364815204002592-gr2.jpg
   :alt: Carbon fluxes in DECOMP

Fig. 2 from Wallman et al. 2006: The conceptual outline of the uppermost horizon in the DECOMP model. The soil is divided into an optional number
of horizons, the lower horizons having only root litter and precipitating DOC as carbon inputs. Each horizon is
treated as a continuously stirred tank reactor. Dashed lines represent release of CO2, gray lines represent DOC
and black lines solid carbon compounds.


