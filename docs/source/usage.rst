Usage
=====

.. _installation:

Installation
------------

To use ``csst-ifs-gehong``, first install it using ``pip`` (This is not available.):

.. code-block:: console

$ pip install csst-ifs-gehong

Execution
----------------

1-D Spectrum
~~~~~~~~~~~~~~~~

The ``spec1d`` module for the modelling of single spectra, including the spectrum of single star, stellar continuum of galaxies, 
emission lines of ionized gas (such as HII region), spectra of AGN. This module also includes some tools for spectral modelling, such 
flux calibration, reddening of dust attenuation. 

Template of Emission Lines
++++++++++++++++++++++++++

The ``spec1d.EmissionLineTemplate``

Emission Line from Ionized Gas
++++++++++++++++++++++++++++++

``spec1d.IonizedGas``

Continuum of Stellar Population in Galaxies
+++++++++++++++++++++++++++++++++++++++++++

``spec1d.StellarPop``

Spectrum of Single Stellar
++++++++++++++++++++++++++

``spec1d.SingleStellar``

Spectrum of AGNs
++++++++++++++++

``spec1d.AGN``

2-D Map
~~~~~~~

``map2d``

Non-parametric Map
++++++++++++++++++

Parametric Map
++++++++++++++

Surface-brightness Map

Velocity Map

Stellar Population (Age and Metallicity) Map

3-D Cube
~~~~~~~~

``cube3d``