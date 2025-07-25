.. _single-stellar-spectra:

Single Stellar Spectra
=======================

This section describes how to generate spectra of individual stars using the ``spec1d.SingleStar`` module in GEHONG.

Template Library
-----------------

Single stellar templates are precomputed from high-resolution spectral libraries. Currently supported:

- ``Munari2005`` (default): Theoretical grid with wide parameter coverage
- ``XSL``: Observed empirical spectra (limited grid)

Use the ``SingleStarTemplate`` class to load the template grid:

.. code-block:: python

    from gehong import spec1d
    starlib = spec1d.SingleStarTemplate(config, template="Munari2005")

Main Parameters:

- ``config``: Simulation configuration
- ``template``: Template type, "Munari2005" or "XSL"

Generating Stellar Spectrum
----------------------------

Use the ``SingleStar`` class to generate a single stellar spectrum from atmospheric parameters.

**Input Parameters**:

- ``mag``: Apparent SDSS r-band magnitude (default: 15.0)
- ``teff``: Effective temperature in Kelvin
- ``logg``: Surface gravity (log10 cm/s²)
- ``feh``: Metallicity [Fe/H] in dex
- ``vel``: Radial velocity in km/s (default: 100)
- ``ebv``: Dust extinction E(B-V) in mag (default: 0.0)

.. code-block:: python

    star = spec1d.SingleStar(config, starlib,
                             mag=14.2, teff=5800, logg=4.4, feh=0.0,
                             vel=150, ebv=0.1)

Output Attributes:

- ``star.wave``: Wavelength array in Ångströms
- ``star.flux``: Flux array in units of :math:`10^{-17}\\ \\mathrm{erg~s^{-1}~cm^{-2}~\mathring{A}^{-1}}`

Notes
-----

- The input stellar parameters will be clipped if outside the template grid.
- Spectra are velocity-shifted, extincted, and calibrated to the input magnitude.
- For ``XSL`` templates, the resolution varies with wavelength; this is accounted for internally.

``SingleStar`` provides a straightforward way to model stellar sources in realistic simulations.
