Usage
=====

.. _installation:

Installation
------------

To use ``csst-ifs-gehong``, first install it using ``pip`` (This is not available.):

.. code-block:: console
    
    $ pip install csst-ifs-gehong

Execution
----------------

1-D Spectrum
~~~~~~~~~~~~~~~~

``spec1d``: The module for the modelling of single spectra, including the spectrum of single star, stellar continuum of galaxies, 
emission lines of ionized gas (such as HII region), spectra of AGN. This module also includes some tools for spectral modelling, such 
flux calibration, reddening of dust attenuation. 

Template of Emission Lines
++++++++++++++++++++++++++

``spec1d.EmissionLineTemplate``: The module for preparing the template of emission line. 

**class** EmissionLineTemplate(self, lam_range = [500, 15000], dlam = 0.1, model = 'fsps', FWHM_inst = 0.5):

*Parameters:*
* lam_range  : the wavelength range of templates ([lower limit, upper limit]), default is [500, 15000]
* dlam       : the wavelength interval of template (unit A), default is 0.1A
* model      : the model of emission line flux ratio, including ('fsps', 'cigale', 'SDSS'), default is 'fsps'
* FWHM_inst  : the FWHM of instrument, default is 0.5


Preparing a class of emission line template. 

.. code-block:: Python

    emline_temp = spec1d.EmissionLineTemplate(model = 'fsps')

Emission Line from Ionized Gas
++++++++++++++++++++++++++++++

``spec1d.IonizedGas``

Continuum of Stellar Population in Galaxies
+++++++++++++++++++++++++++++++++++++++++++

``spec1d.StellarPop``

Spectrum of Single Stellar
++++++++++++++++++++++++++

``spec1d.SingleStellar``

Spectrum of AGNs
++++++++++++++++

``spec1d.AGN``

2-D Map
~~~~~~~

``map2d``

Non-parametric Map
++++++++++++++++++

Parametric Map
++++++++++++++

Surface-brightness Map

Velocity Map

Stellar Population (Age and Metallicity) Map

3-D Cube
~~~~~~~~

``cube3d``