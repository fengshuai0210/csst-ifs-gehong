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

The modelling of AGN spectrum comprises parts: 

* spectrum of broad line region (BLR)
* spectrum of narrow line region (NLR)
* spectrum of power law
* spectrum of FeII line

The final spectrum of AGN is the combination of above four parts. 

The spectrum of broad line region only contains the broad emission line of ionized hydrogen. We use gaussion profile to represent the profile
of each emission line. The final spectrum of BLR is the sum of those emission lines. In the current version, the modelled spectrum only contains
four Balmer lines (from Halpha to Hdelta). The emission line ratio is fixed and adopted as the result of 
`Ilić et al. (2006) <https://ui.adsabs.harvard.edu/abs/2006MNRAS.371.1610I/abstract>`_, which is the measurement for a Seyfert 1.5 galaxy, Mrk 817. 
The profile of BLR emission lines are determined by FWHM of Hbeta emission line, which is adopted as the input parameter in our code. The flux of 
emission lines are determined by the total flux of Hbeta emission line. 

The spectrum of narrow line region contains a series of emission lines of ionized gas, which not only includes hydrogen emission lines but also 
considers those of metals. The model of emissio line flux is adopted as the result of 
`Feltre et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.3354F/abstract>`_, which is calculated by the cloudy. In optical range, this
model contains [OIII]5007, [NII]6583 ... In this model, the flux ratio of emission lines are determined by a series parameters, such as gas-phase metallicity, 
ionized paramter, metal-to-dust ratio, column density of neutral hydrogen, spectrum index of UV photon. For simplisity, we only use the 
gas-phase metallicity as the input parameter to determine the flux ratio of emission lines. The other input parameters are similar with the modelling
of HII region spectrum, such as total flux of Halpha. 

The spectum of power law is modeled by the following,

.. math::

   F_\text{PL} = F_0 \lambda^{\alpha}

where the key parameter is the slope of power law (:math:`\alpha`). The default value is :math:`-1.5`. The absolute flux of power law spectrum is determined 
by the luminosity at :math:`5100\mathring{A}`, which is another input parameter. 

The spectum of FeII lines is modeled by the template of FeII taken from `Park et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJS..258...38P/abstract>`_. 
The wavelength range of templates is from :math:`4000\mathring{A}` to :math:`5600\mathring{A}`. The flux of FeII lines is determined by the flux of broad components
of :math:`H\beta` emission line. Here, we define 

.. math::

   R4570 = \text{Fe}~\text{II}4570 / H\beta

where :math:`\text{Fe}~\text{II}4570` is the flux of FeII emission lines between :math:`4334\mathring{A}` and :math:`4684\mathring{A}`. In observation, the :math:`R4570` is 
between :math:`0.1` and :math:`1.0`. The default value is :math:`0.4`, which is the typical value for AGNs.



Model Maps of Physical Parameters of Galaxies
---------------------------------------------

