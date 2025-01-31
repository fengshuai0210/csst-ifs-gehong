.. _one-dimensional-spectrum-simulation:

One-Dimensional Spectrum Simulation
========================================

This functionality is implemented by the ``spec1d`` module.

Stellar Population Continuum
--------------------------------

Stellar Population Templates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation of the stellar continuum spectrum relies on stellar population templates. 
Thus, these templates need to be imported before conducting the simulation. 
The class ``spec1d.StellarContinuumTemplate`` is used to configure the stellar population templates. 
In the current version, this module employs the single - stellar - population templates from ``emiles``.

**Main Input Parameters**:

- Simulation configuration class: ``config``

**Usage Example**:

.. code-block:: python

    from gehong import spec1d
    stellar_tem = spec1d.StellarContinuumTemplate(config)

Stellar Population Spectrum Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation of the stellar continuum spectrum is realized through the class ``spec1d.StellarContinuum``.

**Main Input Parameters**:

- SDSS - r band magnitude (``mag``): Unit is :math:`\text{mag}`, parameter range is 8 mag to 26 mag.
- Average age of the stellar population (``age``): Unit is :math:`\text{Gyr}`, parameter range is 0.06 Gyr to 17.8 Gyr.
- Average metallicity of the stellar population (``feh``): Unit is :math:`\text{dex}`, parameter range is - 2.32 to 0.22.
- Line - of - sight velocity (``vel``): Unit is :math:`\text{km s}^{-1}`, no specific parameter range.
- Line - of - sight velocity dispersion (``vdisp``): Unit is :math:`\text{km s}^{-1}`, parameter range is greater than 0 km/s.
- Dust extinction (``ebv``): Unit is :math:`\text{mag}`, parameter range is greater than or equal to 0 mag.

**Attributes**:

- Wavelength of the stellar population spectrum (``.wave``): Unit is :math:`\text{\AA}`.
- Flux of the stellar population spectrum (``.flux``): Unit is :math:`10^{-17} \text{erg/s/\AA/cm}^2`.

**Usage Example**:

.. code-block:: python

    stellar = spec1d.StellarContinuum(config, stellar_tem, mag=17, age=2, feh=-0.3, vel=30000, vdisp=150, ebv=0.2)

This code will simulate a stellar population continuum spectrum with an apparent magnitude of 17 mag, 
an average stellar population age of 2 Gyr, a metallicity of [Fe/H] = - 0.3, a line - of - sight velocity of 
30000 km/s, a velocity dispersion of 150 km/s, and a dust extinction of 0.2 mag. The simulated spectrum is as follows:

.. image:: ../../image/example_spec1d_ssp.png


Ionized Gas Emission Line Simulation
----------------------------------------

Ionized Gas Emission Line Templates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``spec1d.EmissionLineTemplate`` is used to configure the templates for ionized gas emission lines. 
A single emission line is simulated using a Gaussian model.

**Main Input Parameters**:

- Simulation configuration class: ``config``
- Emission line template (``model``): In the current version, there are two types of emission line templates: HII region template (``model = 'hii'``) and AGN narrow - line region template (``model = 'nlr'``).

**Usage Example**:

.. code-block:: python

    gas_tem = spec1d.EmissionLineTemplate(config, model = 'hii')

Ionized Gas Emission Line Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``spec1d.HII_Region`` is used to simulate the emission line spectrum of the HII region. 
The unit of the finally simulated spectral flux (``.flux``) is :math:`10^{-17} \text{erg/s/\AA/cm}^2`.

**Main Input Parameters**:

- Integrated flux of the :math:`\text{H}\alpha` emission line (``halpha``): Unit is :math:`10^{-17} \text{erg/s/cm}^2`, no specific parameter range.
- Gas - phase metallicity (``logz``): Unit is :math:`\text{dex}`, parameter range is - 2 to 0.5.
- Line - of - sight velocity (``vel``): Unit is :math:`\text{km s}^{-1}`, no specific parameter range.
- Line - of - sight velocity dispersion (``vdisp``): Unit is :math:`\text{km s}^{-1}`, parameter range is greater than 0 km/s.
- Dust extinction (``ebv``): Unit is :math:`\text{mag}`, parameter range is greater than or equal to 0 mag.

**Attributes**:

- Wavelength of the ionized gas emission line (``.wave``): Unit is :math:`\text{\AA}`.
- Flux of the ionized gas emission line (``.flux``): Unit is :math:`10^{-17} \text{erg/s/\AA/cm}^2`.

**Usage Example**:

.. code-block:: python

    gas = spec1d.HII_Region(config, gas_tem, halpha=500, logz=-0.2, vel=30000, vdisp=150, ebv=0.2)

This code will simulate an ionized gas emission line with an :math:`\text{H}\alpha` flux 
of :math:`500 \times 10^{-17} \text{erg/s/cm}^2`, a gas - phase metallicity of :math:`\log \text{Z/Z}_\odot=-0.3`, 
a line - of - sight velocity of 30000 km/s, a velocity dispersion of 150 km/s, and a dust extinction of 0.2 mag. 
The simulated spectrum is as follows:

.. image:: ../../image/example_spec1d_hii.png


Single - Star Spectrum Simulation
-------------------------------------

[Details about single - star spectrum simulation can be added here if available.]


Active Galactic Nucleus Spectrum Simulation
-----------------------------------------------

The spectrum of an active galactic nucleus consists of four components: the narrow-line region spectrum, 
the broad-line region spectrum, the iron - line spectrum, and the power - law spectrum. 
These four component spectra can be simulated either separately or together.

Power-Law Spectrum Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation of the power - law spectrum in an AGN is implemented through the class ``spec1d.AGN_Powerlaw``.

**Main Input Parameters**:

- Magnitude at 5100 Å (``m5100``): Unit is :math:`\text{mag}`, no specific parameter range.
- Power - law spectral index (``alpha``): Unit is :math:`\text{dex}`, no specific parameter range.
- Line - of - sight velocity (``vel``): Unit is :math:`\text{km s}^{-1}`, no specific parameter range.
- Dust extinction (``ebv``): Unit is :math:`\text{mag}`, parameter range is greater than or equal to 0 mag.

**Attributes**:

- Wavelength of the AGN power - law spectrum (``.wave``): Unit is :math:`\text{\AA}`.
- Flux of the AGN power - law spectrum (``.flux``): Unit is :math:`10^{-17} \text{erg/s/\AA/cm}^2`.

**Usage Example**:

.. code-block:: python

    pl = spec1d.AGN_Powerlaw(config, m5100=17, alpha=-1.5, vel=10000, ebv=0.1)

This code will simulate an AGN power - law spectrum with a magnitude of 17 mag at 5100 Å, 
a power - law spectral index of - 1.5, a line - of - sight velocity of 10000 km/s, and 
a dust extinction of 0.1 mag. The simulated spectrum is as follows:

.. image:: https://note.youdao.com/yws/res/4/WEBRESOURCE3d30cc917638c8f2b27e37249fe5d764


Narrow - Line Region Gas Emission Line Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Narrow - Line Region Gas Emission Line Templates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to the HII region spectrum simulation, first, the class ``spec1d.EmissionLineTemplate`` is used to configure the templates for ionized gas emission lines. The model should be selected as the AGN narrow - line region template (``model = 'nlr'``).

**Main Input Parameters**:

- Simulation configuration class: ``config``
- Emission line template (``model``): In the current version, there are two types of emission line templates: HII region template (``model = 'hii'``) and AGN narrow - line region template (``model = 'nlr'``).

**Usage Example**:

.. code-block:: python

    nlr_temp = spec1d.EmissionLineTemplate(config, model='nlr')

Narrow - Line Region Gas Emission Line Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simulation of the narrow - line region emission line spectrum is achieved through the class ``spec1d.AGN_NLR``.

**Main Input Parameters**:

- Simulation data configuration class: ``config``
- Narrow emission line template class: ``nlr_temp``
- Integrated flux of the :math:`\text{H}\alpha` narrow emission line (``halpha``): Unit is :math:`10^{-17} \text{erg/s/cm}^2`, no specific parameter range.
- Gas - phase metallicity (``logz``): Unit is :math:`\text{dex}`, parameter range is - 2.3 to 0.54.
- Line - of - sight velocity (``vel``): Unit is :math:`\text{km s}^{-1}`, no specific parameter range.
- Line - of - sight velocity dispersion (``vdisp``): Unit is :math:`\text{km s}^{-1}`, parameter range is greater than 0 km/s.
- Dust extinction (``ebv``): Unit is :math:`\text{mag}`, parameter range is greater than or equal to 0 mag.

**Attributes**:

- Wavelength of the AGN narrow emission line spectrum (``.wave``): Unit is :math:`\text{\AA}`.
- Flux of the AGN narrow emission line spectrum (``.flux``): Unit is :math:`10^{-17} \text{erg/s/\AA/cm}^2`.

**Usage Example**:

.. code-block:: python

    nlr = spec1d.AGN_NLR(config, nlr_temp, halpha=100, logz=0, vel=10000, vdisp=400, ebv=0