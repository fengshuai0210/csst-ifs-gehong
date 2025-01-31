.. _ionized-gas-emission-lines:

Ionized Gas Emission Lines
==========================

In current version, the emission lines of ionized gas within galaxies only consider the contribution of HII regions, 
and their generation is implemented by the ``spec1d.HII_Region`` module. 

Emission Line Generation
------------------------

The emission line generation is carried out by employing a series of Gaussian functions to represent the emission lines 
of ionized gas. In our generation, we considered 84 emission lines spanning from :math:`900\mathring{A}` to 
:math:`10500\mathring{A}`. The profile of an emission line can be described as follows:

.. math::

    \mathcal{E}(\lambda) = \frac{1}{\sigma_\text{line}\sqrt{2\pi}}\exp\left[\frac{-(\lambda - \lambda_\text{line})^2}{2\sigma_\text{line}^2}\right]

where :math:`\lambda_\text{line}` is the wavelength of the line center for the emission line, 
and :math:`\sigma_\text{line}` is the line width of emission lines which is determined by the 
velocity dispersion of ionized gas (:math:`\sigma_\text{v}`):

.. math::

    \sigma_\text{line} = \frac{\sigma_\text{v}}{c}\lambda_i

The composite spectra of ionized gas are the combination of all single emission lines:

.. math::

    \mathcal{S}(\lambda) = \sum^{N}_{i}\mathcal{L}_i\mathcal{E}_i(\lambda)

where :math:`\mathcal{L}_i` is the relative flux of :math:`i`th emission lines compared to :math:`\text{H}\alpha` flux, 
and the :math:`\mathcal{E}_i(\lambda)` is the profile of :math:`i`th emission lines.

Relative Flux Determination
---------------------------

The relative flux of the emission lines is determined through the emission line model. We adopt the emission line model 
provided by [Byler2017]_, who used ``Cloudy`` to simulate the emission line flux ratios in a star-forming region ionized 
by a young stellar cluster. In their simulation, the emission line flux ratios depend on the gas-phase metallicity, 
age of the stellar cluster, and ionization parameter. For convenience, we only use metallicity (``logz``) as the sole 
input parameter to determine the relative fluxes of each emission line, while the age of the star cluster and ionization 
parameter adopt a typical fixed value :math:`10^6\text{yrs}` and :math:`\log U = -2`.

Reddening and Redshift Effects
-------------------------------

For the emission lines of ionized gas, the treatment methods for reddening and redshift effects are the same as those for 
the stellar continuum. It is achieved by using the dust extinction (``ebv``) and the line-of-sight velocity of the gas (``vel``) 
as input parameters through Equation \ref{eq:reddening} and Equation \ref{eq:redshift}. The absolute flux of the emission 
lines of ionized gas is scaled according to the input integral flux of the :math:`\text{H}\alpha` emission line (``halpha``).