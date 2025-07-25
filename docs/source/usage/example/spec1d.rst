.. _one-dimensional-spectrum-simulation:

Mocking One-Dimensional Spectrum
================================

The ``spec1d`` module in GEHONG simulates one-dimensional spectra of galaxies, HII regions, stars, and active galactic nuclei (AGN). It supports realistic modeling of stellar populations, ionized gas, AGN emission, and individual stars based on physically motivated or observationally calibrated inputs.

The module includes four major components:

- :ref:`stellar-population-continuum` — Simulation of stellar continuum using SSP templates or star formation histories (SFH)
- :ref:`hii-region-emission-lines` — Simulation of ionized gas emission from HII regions
- :ref:`active-galactic-nuclei-spectra` — Simulation of AGN spectra, including power-law, narrow/broad lines, and Fe II
- :ref:`single-stellar-spectra` — Simulation of individual stellar spectra from empirical or theoretical template libraries

.. toctree::
   :maxdepth: 1
   :caption: spec1d Modules

   stellar
   hii
   agn
   star
