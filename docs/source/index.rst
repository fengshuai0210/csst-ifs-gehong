.. csst-ifs-gehong documentation master file, created by
   sphinx-quickstart on Sun Mar 12 09:05:52 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Gehong's documentation!
==================================

**GEHONG** is a Python package for modelling the data of intergral field spectrascopy mounted on the
Chinese Space Station Telescopy (CSST-IFS), which can also be used for modelling the data of
other IFS instruments. The users can feed a series of maps about the internal properties of
galaxies (e.g. stellar age/metallicity, velocity and velocity dispersion, gas-phase metallicity) to
the package, then the package can output a ``.fits`` file including a datacube for given configure of IFS instrument. 


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Usage

   usage/install.rst
   usage/example/config.rst
   usage/example/spec1d.rst
   usage/example/map2d.rst
   usage/example/cube3d.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Algorithm

   algorithm/sed/sed.rst
   algorithm/map/map.rst
   algorithm/cube/cube.rst

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: API

   api/api.rst

----

:Author: Shuai Feng
:Institute: Hebei Normal University
:Contact: sfeng@hebtu.edu.cn
:Last updated: 2025-01-31
:Version: 3.0.1
:Copyright: CSST-IFS Team

----

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`