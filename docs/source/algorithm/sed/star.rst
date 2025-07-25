.. _single-star-spectrum:

Single Star Spectrum
====================

In the current version, the spectrum of a single star is modeled using the ``spec1d.SingleStar`` class, which relies on template libraries managed by the ``SingleStarTemplate`` class.

Overview
--------

Two stellar template libraries are supported:

- ``XSL``: Empirical stellar spectra from the `X-shooter Spectral Library (XSL) <http://xsl.u-strasbg.fr/>`_, observed with VLT/X-shooter 
  (`Verro et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26A...660A..34V/abstract>`_).
- ``Munari2005``: Synthetic spectra from `Munari et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005A%26A...442.1127M/abstract>`_,
  adopting a constant instrumental resolution of 1 Å.

Each template is characterized by the stellar atmospheric parameters:

- Effective temperature :math:`T_{\rm eff}` [K],
- Surface gravity :math:`\log g` [cm/s²],
- Metallicity :math:`[{\rm Fe/H}]` [dex].

Template Selection
------------------

The stellar template that best matches the input parameters is selected through the following steps:

1. Find the nearest :math:`T_{\rm eff}` from the library.
2. Within that subset, find the nearest :math:`\log g`.
3. Within that subset, find the nearest :math:`[{\rm Fe/H}]`.

The selected template is normalized by its mean value.

Spectral Broadening
-------------------

The intrinsic resolution of each library is taken into account:

- ``XSL``: wavelength-dependent instrumental resolution with :math:`R \sim 10000`,
- ``Munari2005``: constant resolution of 1 Å across the full wavelength range.

These are converted to wavelength-dependent LSF and stored as a vector. When needed, Gaussian convolution is applied to match target velocity dispersion.

Dust Attenuation
----------------

Internal extinction is applied using the attenuation law of `Calzetti et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_:

.. math::

   F_{\rm obs}(\lambda) = F_{\rm int}(\lambda) \cdot 10^{-0.4 \cdot E(B-V) \cdot k(\lambda)}

where :math:`E(B-V)` is the user input `ebv`, and :math:`k(\lambda)` is the attenuation curve.

Redshift and Doppler Shift
--------------------------

The observed wavelength is computed using the stellar line-of-sight velocity `vel`, including relativistic correction:

.. math::

   \lambda_{\rm obs} = \lambda_{\rm rest} \cdot \left(1 + \frac{v_{\rm los}}{c} \right)

This ensures accurate transformation even for high-velocity stars.

Flux Calibration
----------------

The final flux is scaled to match the specified apparent SDSS-:math:`r` band magnitude (`mag`). The calibration uses the filter curve ``SLOAN_SDSS.r``, and the final output is expressed in :math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \mathring{A}^{-1}}`.

Output
------

The final 1D spectrum includes:

- Stellar template matched by :math:`T_{\rm eff}`, :math:`\log g`, and :math:`[{\rm Fe/H}]`,
- Gaussian spectral broadening (if needed),
- Internal dust extinction,
- Relativistic Doppler redshift from line-of-sight velocity,
- Magnitude-based flux calibration.

This class is useful for simulating isolated stars or star clusters.
