.. _section-5-three-dimensional-data-cube-simulation:

Mocking Three-Dimensional Cube
==============================

This page introduces the usage of the ``cube3d`` module to simulate integral-field spectral cubes.

Overview
--------

The module simulates a 3D datacube with dimensions (RA, DEC, Wavelength), combining both stellar continuum and ionized gas emission. The key class is:

.. code-block:: python

    from gehong import cube3d
    u = cube3d.Cube3D(config, stellar_map=stellarcontinuum, gas_map=ionizedgas)

Where:
- ``config`` is the global simulation configuration object.
- ``stellar_map`` is a ``StellarPopulationMap`` instance.
- ``gas_map`` is an ``IonizedGasMap`` instance.

StellarPopulationMap Inputs
---------------------------

The stellar population distribution is defined by ``map2d.StellarPopulationMap``. Two working modes are supported:

**(1) Surface brightness–calibrated mode (empirical):**

.. code-block:: python

    stellar_map = map2d.StellarPopulationMap(
        config,
        sb=sb_map,
        sfh=sfh_map,          # optional
        ceh=ceh_map,          # optional
        vel=vel_map,
        vdisp=vdisp_map,
        ebv=ebv_map
    )

- ``sb``: 2D surface brightness in mag/arcsec\ :sup:`2`
- ``sfh``: Scalar, 2D, or 3D star formation history
- ``ceh``: Scalar, 2D, or 3D chemical enrichment history
- ``vel``, ``vdisp``, ``ebv``: Kinematic and dust maps

**(2) Physically calibrated mode (from SFH):**

.. code-block:: python

    stellar_map = map2d.StellarPopulationMap(
        config,
        sfh=sfh_cube,
        sfh_age=sfh_age,
        ceh=ceh_cube,
        ceh_age=ceh_age,
        vdisp=vdisp_map,
        ebv=ebv_map,
        z=0.01,                 # scalar redshift
        vpec=vpec_map          # peculiar velocity map
    )

- ``sfh``: 3D array [nx, ny, nage] of SFH
- ``sfh_age``: 1D age grid
- ``z``: Scalar redshift (required)
- ``vpec``: Peculiar velocity map (km/s)

IonizedGasMap Inputs
--------------------

Ionized gas emission is defined by ``map2d.IonizedGasMap``. Two modes are also supported:

**(1) Halpha-calibrated mode:**

.. code-block:: python

    gas_map = map2d.IonizedGasMap(
        config,
        halpha=halpha_map,
        zh=zh_map,
        vel=vel_map,
        vdisp=vdisp_map,
        ebv=ebv_map
    )

- ``halpha``: 2D map of H\ :sub:`α` flux (in 10\ :sup:`-17` erg/s/cm\ :sup:`2`)
- ``zh``: Gas metallicity map
- ``vel``, ``vdisp``, ``ebv``: Kinematics and extinction

**(2) Physically calibrated mode (from SFR):**

.. code-block:: python

    gas_map = map2d.IonizedGasMap(
        config,
        sfr=sfr_map,
        zh=zh_map,
        z=0.01,               # scalar redshift
        vpec=vpec_map,
        vdisp=vdisp_map,
        ebv=ebv_map
    )

- ``sfr``: 2D star formation rate map
- ``z``: Scalar redshift (required)
- ``vpec``: Peculiar velocity map

Template Configuration
----------------------

To simulate spectra, first configure spectral templates:

.. code-block:: python

    from gehong import spec1d
    stellar_tem = spec1d.StellarContinuumTemplate(config)
    gas_tem = spec1d.EmissionLineTemplate(config, model='hii')

Generate the Cube
-----------------

Use the ``make_cube`` method of ``Cube3D`` to generate the full spectrum:

.. code-block:: python

    u.make_cube(stellar_tem=stellar_tem, hii_tem=gas_tem)

This generates a cube with shape (nx, ny, nwave) and unit of 10\ :sup:`-17` erg/s/cm\ :sup:`2`/\ Å.

Saving the Cube
---------------

You can export the cube to a FITS file with full WCS header:

.. code-block:: python

    u.savefits("mock_cube.fits")

The output will contain a primary HDU with metadata and a data HDU with the flux cube.

Inserting Single Spectra
------------------------

To insert an individual spectrum into the cube center or with an offset:

.. code-block:: python

    u.insert_spec(spec, dx_arcsec=0.0, dy_arcsec=0.0)

This is useful for injecting AGN or other localized components.
