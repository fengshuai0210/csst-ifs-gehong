GEHONG: A Python Package for Simulating IFS Spectral Data Cubes
================================================================

.. image:: /_static/gehong_logo.png
   :alt: GEHONG Logo
   :align: center
   :height: 250px



**GEHONG** (GEnerate tHe data Of iNtegral field spectrograph of Galaxy) is a Python package
developed for simulating high-resolution, three-dimensional (3D) spectroscopic data cubes
observed with ideal integral field spectrographs (IFS). It is specifically designed to support
the science preparation of the **Chinese Space Station Telescope Integral Field Spectrograph
(CSST-IFS)**, but is also applicable to other IFS instruments.

Given two-dimensional maps describing spatially resolved physical properties of target
sources, such as stellar population, surface brightness, line-of-sight velocity, and gas-phase
metallicity, **GEHONG** generates synthetic 3D datacubes by assembling position-dependent
one-dimensional spectra. These spectra can be composed of up to four components:

- stellar continuum
- ionized gas emission lines
- AGN spectra
- individual stellar spectra

**GEHONG** is capable of producing mock observations of a wide range of astrophysical systems,
including galaxies, AGNs, star clusters, and H II regions. The input parameter maps can be
derived from either analytic models, observational data, or hydrodynamic simulations.

**GEHONG** plays a key role in simulating idealized datacubes for CSST-IFS scientific pre-studies.
These datacubes can serve as inputs for:

- the CSST-IFS raw image simulator (**CSST-IFS-SIM**), which generates mock CCD images
  from ideal datacubes. `Online access <https://csst-tb.bao.ac.cn/code/csst-sims/csst_ifs_sim>`_

- the exposure time calculator (**CSST-IFS-ETC**), which estimates the signal-to-noise
  ratio of mock observations and helps optimize exposure configurations. This tool can be
  accessed both as a `Python package <https://pypi.org/project/ifs-etc/>`_ and via a
  `web-based interface <http://60.204.173.65:9091/#/tool/ifs-etc-task>`_

- other planning tools for observation simulation and survey strategy design



----

.. toctree::
   :maxdepth: 2
   :caption: 🚀 Getting Started
   :hidden:

   usage/install.rst
   usage/quickstart.rst

.. toctree::
   :maxdepth: 2
   :caption: 📘 User Guide
   :hidden:

   usage/example/config.rst
   usage/example/spec1d.rst
   usage/example/map2d.rst
   usage/example/cube3d.rst

.. toctree::
   :maxdepth: 2
   :caption: ⚙️ Algorithms & Models
   :hidden:

   algorithm/sed/sed.rst
   algorithm/map/map.rst
   algorithm/cube/cube.rst

.. toctree::
   :maxdepth: 2
   :caption: 🧪 API Reference
   :hidden:

   api/api.rst

.. toctree::
   :maxdepth: 2
   :caption: 📄 About
   :hidden:

   logo.rst
   changelog.rst
   license.rst
   cite.rst


Project Info
------------

- **Author**: Shuai Feng  
- **Institute**: Hebei Normal University  
- **Contact**: sfeng@hebtu.edu.cn  
- **Version**: 3.1.0 
- **Last Updated**: 2025-07-25 
- **License**: MIT  
- **GitHub**: https://github.com/fengshuai0210/csst-ifs-gehong  
- **Documentation**: https://csst-ifs-gehong.readthedocs.io/  
