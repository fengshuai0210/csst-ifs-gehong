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

* The spectrum of broad line region only contains the broad emission line of ionized hydrogen. We use gaussion profile to represent the profile
of each emission line. The final spectrum of BLR is the sum of those emission lines. In the current version, the modelled spectrum only contains
four Balmer lines (from Halpha to Hdelta). The emission line ratio is fixed and adopted as the result of 
`Ilić et al. (2006) <https://ui.adsabs.harvard.edu/abs/2006MNRAS.371.1610I/abstract>`_, which is the measurement for a Seyfert 1.5 galaxy, Mrk 817. 
The profile of BLR emission lines are determined by FWHM of Hbeta emission line, which is adopted as the input parameter in our code. The flux of 
emission lines are determined by the 

宽线区发射线的谱形由宽线Hβ的积分流量、宽线Hβ的FWHM两个参数决定。（由于宽线区的绝对流量由宽线Hβ积分流量决定，因此宽线区发射线与FeII发射线流量相关）

* The spectrum of narrow line region contains 

Spectrum of Single Stellar
~~~~~~~~~~~~~~~~~~~~~~~~~~

Model Maps of Physical Parameters of Galaxies
---------------------------------------------

