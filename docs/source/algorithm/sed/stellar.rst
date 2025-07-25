.. _stellar-continuum:

Stellar Continuum Spectrum
==========================

In the current version, the stellar continuum spectrum is modeled using the ``spec1d.StellarContinuum`` class. This module supports multiple methods for constructing the intrinsic spectral shape, and two distinct modes for flux calibration.

Overview
--------

The stellar spectrum is constructed from single stellar population (SSP) templates. We adopt the `E-MILES <http://miles.iac.es/pages/stellar-libraries/miles-library.php>`_ library (`Vazdekis et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.3409V/abstract>`_), which provides a broad wavelength range from :math:`1680\ \mathring{\mathrm{A}}` to :math:`50000\ \mathring{\mathrm{A}}`, and dense sampling in age and metallicity.

Two options for constructing the spectral shape:

- **Composite Stellar Population (CSP)**: Assemble from a star formation history (SFH) and chemical enrichment history (CEH).
- **Single SSP**: Directly use one SSP template with given age and metallicity.

Two modes for flux calibration:

- **Magnitude-calibrated mode**: Use SDSS-r magnitude (`mag`) to normalize output spectrum (required for single SSP).
- **SFR-calibrated mode**: Use physically normalized SFH to produce absolute flux (only available for CSP).

Spectral Construction
---------------------

In CSP mode, the stellar continuum is built from a linear combination of SSP templates:

.. math::

   S(\lambda) = \sum_i w_i \cdot \mathrm{SSP}_i(\lambda)

where :math:`w_i` is the weight from SFH at lookback time :math:`t_i`, and :math:`\mathrm{SSP}_i(\lambda)` is the template spectrum for that age and metallicity.

The input SFH and CEH must be arrays with lookback time grids:

- SFH input: array of shape ``(nage, 2)``, with columns = [lookback time, relative SFR].
- CEH input: array of shape ``(nage, 2)``, with columns = [lookback time, [Fe/H]].

Only the relative SFR distribution is used in magnitude-calibrated mode; the absolute mass normalization is applied in SFR-calibrated mode.

In single SSP mode, the user provides `age` and `feh`, and the closest matching template is selected.

Velocity Dispersion Broadening
------------------------------

The E-MILES templates include intrinsic instrumental broadening, which varies with wavelength:

- UV band (:math:`\lambda < 3541\ \mathring{\mathrm{A}}`): :math:`\sigma \approx 200\ \mathrm{km\ s^{-1}}`
- Optical/NIR bands: :math:`\sigma < 100\ \mathrm{km\ s^{-1}}`

If the input stellar velocity dispersion (`vdisp`) exceeds the intrinsic value, additional broadening is applied using Gaussian convolution in log-wavelength space, following `PPXF <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C/abstract>`_.

Dust Attenuation
----------------

Dust attenuation is applied using the `Calzetti et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_ law:

.. math::

   S_{\rm dust}(\lambda) = S_{\rm nodust}(\lambda) \cdot 10^{-0.4 \cdot E(B-V) \cdot k(\lambda)}

where :math:`k(\lambda)` is the attenuation curve and :math:`E(B-V)` is the reddening input parameter (`ebv`).

Redshift and Peculiar Velocity
------------------------------

The observed-frame wavelength is computed using the cosmological redshift (`z`) and the line-of-sight peculiar velocity (`vpec`), as:

.. math::

   v_{\rm los} = (1 + z)\cdot v_{\rm pec} + z \cdot c

.. math::

   \lambda_{\rm obs} = \lambda_{\rm rest} \cdot \left(1 + \frac{v_{\rm los}}{c}\right)

This relativistic correction is consistent with the treatment used for gas emission lines.

Flux Calibration
----------------

Two options are supported:

- **Magnitude-calibrated mode**:

  - Input: SDSS-r apparent magnitude (`mag`)
  - Output flux is scaled to match `mag` after redshift and attenuation.
  - Available for both single SSP and CSP inputs.

- **SFR-calibrated mode**:

  - Input: normalized SFH + `z`
  - Output reflects physically motivated absolute flux assuming total stellar mass.
  - Only available when using SFH + CEH (not for single SSP).

Output
------

The final output is a 1D NumPy array representing the observed-frame spectrum, with:

- Spectral line broadening,
- Dust attenuation,
- Redshift and peculiar velocity correction,
- Proper flux normalization based on the chosen calibration mode.

Unit: :math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \mathring{A}^{-1}}`
