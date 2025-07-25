.. _hii-region-emission-lines:

HII Region Emission Lines
==========================

This section describes how to simulate emission-line spectra from HII regions using the ``spec1d.HII_Region`` class in GEHONG.

Emission Line Template
-----------------------

Before simulating spectra, an emission line template must be initialized using the ``EmissionLineTemplate`` class.  
For HII regions, set ``model='hii'``.

**Usage Example**:

.. code-block:: python

    from gehong import spec1d
    gas_tem = spec1d.EmissionLineTemplate(config, model='hii')


HII Region Spectrum Simulation
------------------------------

The ``HII_Region`` class supports two working modes for generating integrated emission-line spectra.

Flux-calibrated mode** (empirical Hα normalization)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Activated when ``halpha`` is provided.

**Input Parameters**:

- ``halpha``: Total Hα flux (:math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}}`)
- ``vel``: Line-of-sight velocity (km/s)
- ``vdisp``: Velocity dispersion (km/s)
- ``logz``: Gas-phase metallicity [Z/H] (dex)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    gas = spec1d.HII_Region(config, gas_tem,
                            halpha=500, vel=8000, vdisp=80,
                            logz=-0.3, ebv=0.1)


.. image:: /_static/image/example_spec1d_hii1.png
   :align: center


Physically calibrated mode (based on star formation rate)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Activated when ``halpha`` is set to ``None``, and ``sfr`` and ``z`` are provided.

**Input Parameters**:

- ``sfr``: Star formation rate (M☉/yr)
- ``z``: Redshift
- ``vpec``: Peculiar velocity (km/s)
- ``logz``: Gas-phase metallicity [Z/H] (dex)
- ``vdisp``: Velocity dispersion (km/s)
- ``ebv``: Dust extinction E(B-V)

**Usage Example**:

.. code-block:: python

    gas = spec1d.HII_Region(config, gas_tem,
                            sfr=0.5, z=0.01, vpec=300,
                            logz=-0.2, vdisp=60, ebv=0.2)


.. image:: /_static/image/example_spec1d_hii2.png
   :align: center

Output Attributes
------------------

- ``gas.wave``: 1D wavelength array (Ångströms)
- ``gas.flux``: Corresponding flux array in units of :math:`10^{-17}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}}`
