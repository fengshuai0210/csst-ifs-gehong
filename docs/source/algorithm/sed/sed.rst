One-Dimensional Spectrum 
=========================

GEHONG provides physical models for simulating one-dimensional (1D) spectra of various astrophysical sources. These models are the building blocks for constructing spatially resolved 2D and 3D spectral datacubes.

Supported Components
--------------------

GEHONG supports the following 1D spectral components:

- **Stellar Continuum**  
  Based on user-defined star formation history (SFH) and chemical enrichment history (CEH), synthesized from SSP templates.

- **Ionized Gas Emission Lines**  
  H II region spectra modeled using photoionization grids and scaled by star formation rate (SFR) or Hα luminosity.

- **Active Galactic Nuclei (AGN)**  
  Including power-law continuum, broad-line region (BLR), narrow-line region (NLR), and Fe II templates.

- **Individual Stars**  
  High-resolution theoretical spectra generated from stellar atmosphere grids, based on Teff, logg, and [Fe/H].

Subsections
-----------

.. toctree::
   :maxdepth: 2
   :caption: Modeling of Single Spectrum

   gas.rst
   stellar.rst
   agn.rst
   star.rst
