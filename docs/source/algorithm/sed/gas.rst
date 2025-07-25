.. _ionized-gas-emission-lines:

Ionized Gas Emission Lines
==========================

In the current version, the emission lines of ionized gas in galaxies are modeled exclusively as originating from star-forming HII regions. These are implemented through the ``spec1d.HII_Region`` class.

Overview
--------

The synthetic spectrum of ionized gas is constructed by combining multiple emission lines with Gaussian profiles. The base templates come from FSPS-Cloudy simulations (`Byler et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017ApJ...840...44B/abstract>`_), which tabulate emission line flux ratios under various gas metallicities.

Two flux calibration modes are supported:

- **Hα-calibrated mode**: If an observed Hα flux is provided, the spectrum is scaled accordingly.
- **SFR-calibrated mode**: If no Hα flux is given but `sfr`, `z`, and `vpec` are provided, the Hα luminosity is inferred via the `Kennicutt (1998) <https://ui.adsabs.harvard.edu/abs/1998ARA%26A..36..189K/abstract>`_ relation.

Emission Line Generation
------------------------

Each emission line is represented by a Gaussian function:

.. math::

   \mathcal{E}(\lambda) = \frac{1}{\sigma_{\rm line} \sqrt{2\pi}} \exp\left[ -\frac{(\lambda - \lambda_{\rm line})^2}{2\sigma_{\rm line}^2} \right]

where :math:`\sigma_{\rm line}` is the line width in wavelength space, computed from the velocity dispersion :math:`\sigma_v` via:

.. math::

   \sigma_{\rm line} = \frac{\sigma_v}{c} \lambda

The total spectrum is given by the weighted sum of all lines:

.. math::

   \mathcal{S}(\lambda) = \sum_i \mathcal{L}_i \mathcal{E}_i(\lambda)

where :math:`\mathcal{L}_i` is the flux of the :math:`i`th line relative to Hα, determined from the model template.

Template Selection
------------------

The emission line templates are selected based on gas-phase metallicity (`logz`). Available models include:

- ``'hii'``: HII region templates based on FSPS + Cloudy simulations.
- ``'nlr'``: Narrow-line region templates for AGN (used in other modules).

The nearest metallicity grid point is used; no interpolation is performed. The range of supported metallicity is:

.. math::

   \log Z/Z_\odot \in [-2.0, 0.5]

Inputs outside this range will be clipped with a warning.

Velocity Broadening
-------------------

Line broadening due to gas velocity dispersion is applied using a pixel-wise Gaussian convolution:

- The intrinsic resolution of the line template is deconvolved from the target dispersion.
- The convolution kernel is constructed per pixel using an adaptive window.
- This process is implemented in the `gaussian_filter1d` function.

Dust Attenuation
----------------

Dust extinction is applied using the `Calzetti et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_ attenuation law:

.. math::

   F_{\rm obs}(\lambda) = F_{\rm int}(\lambda) \cdot 10^{-0.4 \cdot E(B-V) \cdot k(\lambda)}

where :math:`k(\lambda)` is the extinction curve derived from the input `ebv` and `Rv` parameters. This reddening is applied to the entire line spectrum prior to flux calibration.

Redshift and Peculiar Velocity
------------------------------

The line-of-sight redshift of the emission spectrum is determined by either:

- the input `vel` (if `sfr` is not used), or
- the cosmological + peculiar velocity (`z`, `vpec`) if using SFR-based calibration:

.. math::

   v_{\rm los} = (1 + z) \cdot v_{\rm pec} + z \cdot c

The final spectrum is shifted accordingly using spline interpolation in wavelength space.

SFR-Calibrated Mode
-------------------

If the Hα flux is not given, the model estimates the Hα luminosity using:

.. math::

   L({\rm H}\alpha) = 7.9 \times 10^{41} \cdot {\rm SFR}~[\text{erg/s}]

This is converted to observed flux using the luminosity distance derived from redshift:

.. math::

   F_{\rm obs} = \frac{L({\rm H}\alpha)}{4\pi D_L^2 (1+z)^4}

Here, the :math:`(1+z)^4` term accounts for cosmological surface brightness dimming. The emission line spectrum is then scaled to match this predicted Hα flux.

Output Spectrum
---------------

The final synthetic spectrum includes:

- Proper flux calibration (via Hα or SFR),
- Dust attenuation,
- Gaussian broadening by velocity dispersion,
- Redshift and Doppler shift,
- Wavelength resampling onto the user-defined grid (`config.wave`).

All emission lines are modeled together and returned as a 1D array representing the final observed-frame flux.

``HII_Region`` is designed to be used in both single-spectrum and datacube simulations.
