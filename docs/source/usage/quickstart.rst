Quickstart Guide
================

This page provides a starting point for learning how to use **GEHONG** through example-driven workflows.

All tutorials are provided as **Jupyter notebooks** and hosted in the following GitHub repository:

🔗 `csst-ifs-gehong-example`  
https://github.com/fengshuai0210/csst-ifs-gehong-example

You can clone the repository with:

.. code-block:: bash

   git clone https://github.com/fengshuai0210/csst-ifs-gehong-example.git

This repository contains step-by-step examples demonstrating how to use GEHONG for different scientific goals.

Available Notebooks
-------------------

The following example notebooks are planned or under development:

- ``01_config_intro.ipynb``  
  Introduction to the `Config` object and simulation grid settings.

- ``02_map2d_example.ipynb``  
  Generate 2D parameter maps (e.g., stellar mass, velocity, SFR, metallicity).

- ``03_cube3d_starforming.ipynb``  
  Simulate a 3D spectral cube for a star-forming galaxy, using stellar and nebular emission.

- ``04_cube3d_agn.ipynb`` *(planned)*  
  Include AGN components in the cube (continuum + NLR + BLR + FeII).

- ``05_single_star_spectrum.ipynb``  
  Generate high-resolution spectra of individual stars based on atmospheric parameters.

- ``10_batch_generation_pipeline.ipynb`` *(planned)*  
  An end-to-end example of batch-simulating a population of mock galaxies for survey preparation.

Usage Instructions
------------------

Each notebook is self-contained. After installing GEHONG and setting the resource path, you can open any notebook and execute the cells directly.

GEHONG must be installed in the same Python environment where you launch Jupyter:

.. code-block:: bash

   jupyter notebook  # or jupyter lab
