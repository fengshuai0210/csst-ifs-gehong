Spectrum of Stellar Population
==============================

`spec1d.StellarContinuum`

The spectrum of stellar pupulation is constructed by the template of single stellar population (SSP). 
In the current version, we adopt the `emiles <http://miles.iac.es/pages/stellar-libraries/miles-library.php>`_ 
library. The emiles templates have a broad wavelength coverage, spanning from :math:`1680\mathring{A}` in the ultraviolet to 
nearly :math:`50000\mathring{A}` in the near-infrared, which is greater than the instrumental detection range of CSST-IFS. 
Simultaneously, they also have a wide coverage in terms of stellar population age and metallicity. 

Intrinsic Spectral Shape
~~~~~~~~~~~~~~~~~~~~~~~~

We use the SSP template to mock the intrinsic spectral shape of the stellar continuum spectrum in a specific area 
of the galaxy. Given the high spatial resolution of CSST-IFS, the physical scale covered by one spaxel is relatively 
small (:math:`\sim 0.04 \text{kpc}` at :math:`z = 0.01`) for local galaxies. We assume that the stellar population 
within such a small area can be approximated by a single stellar population is reasonable. Thus, the general shape 
of stellar population spectra without dust attenuation is determined by the input median age (`age`) and metallicity (`feh`).

Spectral Line Broadening
~~~~~~~~~~~~~~~~~~~~~~~~

We take into account the spectral line broadening effect resulting from the irregular motion of stars, 
which can be quantified by the stellar velocity dispersion (`vdisp`). Due to the resolution limit of the 
spectrum, the emiles SSP template spectrum already has some spectral line broadening, and this broadening 
is a function of wavelength. For the ultraviolet band (:math:`\lambda < 3541\mathring{A}`), the velocity 
dispersion corresponding to its broadening can reach :math:`200\kms`, while for the optical and near-infrared 
bands, the velocity dispersion corresponding to the broadening is less than :math:`100\kms`. During the 
simulation of the galaxy spectrum, if the input velocity dispersion is smaller than the resolution of the 
template spectrum, no further broadening is performed on the template spectrum. If the input velocity dispersion
is greater than the resolution of the template spectrum, further broadening of the spectrum is carried out. 
Here, we directly adopt the method of `PPXF` [Cappellari2017]_ for spectral broadening. That is, we assume 
that the broadening caused by the velocity dispersion is a Gaussian function, and then perform convolution 
using the fast Fourier method in the natural logarithm of the wavelength coordinate.

Dust Attenuation
----------------

We only consider the intrinsic dust attenuation within the target galaxy and do not account for the foreground 
dust extinction of the Milky Way. We perform reddening on the spectrum based on the input dust extinction 
value (:math:`\text{EBV}`) to simulate the dust attenuation effect on the observed spectra. 
We adopt the dust extinction law of [Calzetti2000]_. At the rest wavelength, the intrinsic continuum 
and the reddened continuum can be expressed as follows:

.. math::

    \mathcal{S}_\text{dust}(\lambda) = \mathcal{S}_\text{nodust}(\lambda) \times e^{0.921 \ebv k(\lambda)}. 

In this equation, :math:`k(\lambda)` is the dust attenuation law from [Calzetti2000]_. 

Redshift Effect
~~~~~~~~~~~~~~~

Then, we add the redshift effect to the spectrum, which includes the cosmological redshift due 
to the expansion of the universe and the Doppler redshift caused by the motion of the galaxy. 
We use the line-of-sight velocity relative to the Galactic Standard of Rest, which is another 
input parameter (`vel`), to quantify the observed redshift. Based on the input value of the 
line-of-sight velocity (:math:`v_\text{LOS}`), we shift the spectrum at the rest wavelength 
(:math:`\lambda_\text{rest}`) along the wavelength dimension to obtain the observed spectral shape:

.. math::

    \lambda_\text{obs} = \lambda_\text{rest} (1 + \frac{v_\text{LOS}} {c})

where :math:`\lambda_\text{obs}` is the observed wavelength.

Flux Calibration
----------------

We perform flux calibration on the simulated spectra. We use the SDSS - :math:`r` band magnitude 
value as the input parameter (`mag`) and scale the redshifted spectrum obtained in the previous 
step to the observed spectrum with the flux unit of :math:`10^{-17} \text{erg/s}/\mathring{A}/\text{cm}^2`.
