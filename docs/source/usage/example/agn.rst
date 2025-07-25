.. _active-galactic-nuclei-spectra:

Active Galactic Nuclei (AGN)
============================

This section describes how to simulate the spectra of active galactic nuclei (AGN) using the ``spec1d`` module in GEHONG.

AGN spectra are composed of multiple components, which can be generated individually or combined into a full AGN model.

Emission Line Template
-----------------------

The narrow-line region (NLR) template is initialized using the ``EmissionLineTemplate`` class with ``model='nlr'``.

**Usage Example**:

.. code-block:: python

    from gehong import spec1d
    nlr_temp = spec1d.EmissionLineTemplate(config, model='nlr')

AGN Component Simulations
--------------------------

Power-Law Continuum
~~~~~~~~~~~~~~~~~~~~

Simulates the AGN featureless continuum as a power law.

**Input Parameters**:

- ``m5100``: Apparent magnitude at 5100 Å
- ``alpha``: Spectral slope index
- ``vel``: Line-of-sight velocity (km/s)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    pl = spec1d.AGN_Powerlaw(config, m5100=17, alpha=-1.5, vel=10000, ebv=0.1)


.. image:: /_static/image/example_spec1d_agn_pl.png
   :align: center


Narrow-Line Region (NLR)
~~~~~~~~~~~~~~~~~~~~~~~~~

Simulates the narrow emission lines in AGN spectra.

**Input Parameters**:

- ``halpha``: Hα narrow-line flux (:math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}}`)
- ``logz``: Gas-phase metallicity [Z/H] (dex)
- ``vel``: Line-of-sight velocity (km/s)
- ``vdisp``: Velocity dispersion (km/s)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    nlr = spec1d.AGN_NLR(config, nlr_temp, halpha=300, logz=0, 
                         vel=10000, vdisp=300, ebv=0.1)


.. image:: /_static/image/example_spec1d_agn_nlr.png
   :align: center


Broad-Line Region (BLR)
~~~~~~~~~~~~~~~~~~~~~~~~

Simulates the broad Hβ emission from the BLR.

**Input Parameters**:

- ``hbeta_flux``: Hβ broad-line flux
- ``hbeta_fwhm``: Full width at half maximum (FWHM) of Hβ line (km/s)
- ``vel``: Line-of-sight velocity (km/s)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    blr = spec1d.AGN_BLR(config, hbeta_flux=150, hbeta_fwhm=4000, 
                         vel=10000, ebv=0.1)


.. image:: /_static/image/example_spec1d_agn_blr.png
   :align: center


Fe II Complex
~~~~~~~~~~~~~~~

Simulates the blended Fe II emission near Hβ.

**Input Parameters**:

- ``hbeta_broad``: Hβ broad-line flux (used to scale Fe II)
- ``r4570``: Fe II to Hβ flux ratio
- ``vel``: Line-of-sight velocity (km/s)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    fe = spec1d.AGN_FeII(config, hbeta_broad=150, r4570=0.5, 
                         vel=10000, ebv=0.1)


.. image:: /_static/image/example_spec1d_agn_fe.png
   :align: center


Combining AGN Components
-------------------------

Each component produces a spectrum on the same wavelength grid, and they can be added to form the full AGN model.

**Example**:

.. code-block:: python

    flux = pl.flux + nlr.flux + blr.flux + fe.flux

.. image:: /_static/image/example_spec1d_agn_sum.png
   :align: center

Output Attributes
------------------

- ``component.wave``: Wavelength array (Å)
- ``component.flux``: Flux array in units of :math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}}`

.. note::

   All AGN components share the same wavelength grid defined by ``config``.
