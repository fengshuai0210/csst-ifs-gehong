Algorithm
=========

Model Single Spectra
--------------------

Emission Lines
~~~~~~~~~~~~~~

Continuum of Stellar Components in Galaxies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spectrum of stellar continuum is modelled by the template of single stellar population (SSP). In the current version, we adopt 
the `emiles <http://miles.iac.es/pages/stellar-libraries/miles-library.php>`_ library. Then, given the age and metallicity, we could 
obtain the profile of stellar continuum. Therefore, for the simplest case, according to the input `age` and `FeH`, the code will output
a spectra of stellar continuum. Such a spectra is not flux-callibrated. To obtain a spectrum with absolute flux (unit: 1e1-7 erg/s/A/cm^2), 
the code need a magnitude at given band (e.g. SDSS r-band) as an input parameter to calibrate the flux of spectra. 

Spectrum of AGN
~~~~~~~~~~~~~~~

The modelling of AGN spectrum comprises parts: spectrum of broad line region (BLR), spectrum of narrow line region (NLR), spectrum of 
power law, spectrum of FeII line. The final spectrum of AGN is the combination of above four parts. 

1. The spectrum of broad line region only contains the broad emission line of ionized hydrogen. We use gaussion profile to represent the profile
of each emission line. The final spectrum of BLR is the The emission line ratio is fixed and adopted 
by the 

Spectrum of Single Stellar
~~~~~~~~~~~~~~~~~~~~~~~~~~

Model Maps of Physical Parameters of Galaxies
---------------------------------------------

