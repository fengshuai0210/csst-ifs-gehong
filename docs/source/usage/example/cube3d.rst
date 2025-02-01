.. _section-5-three-dimensional-data-cube-simulation:

Mocking Three-Dimensional Cube
===========================================

This simulation is implemented by the ``cube3d`` module. In the following, an example will be used to introduce 
the implementation of the functions of the ``cube3d`` module.

Initialization of the Three - Dimensional Data Cube
-------------------------------------------------------

The initialization of the three - dimensional data cube is achieved through the class ``cube3d.Cube3D``. For example:

.. code-block:: python

    from gehong import cube3d
    u = cube3d.Cube3D(config, stellar_map=stellarcontinuum, gas_map=ionizedgas)

Here, ``config`` is the simulation data configuration class. ``stellarcontinuum`` and ``ionizedgas`` are classes 
that contain the two - dimensional distributions of stellar population parameters and ionized gas parameters respectively.

Two - Dimensional Distribution of Stellar Population Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``stellarcontinuum``, which contains the two - dimensional distribution of stellar population parameters, 
is created as follows:

.. code-block:: python

    stellarcontinuum = map2d.StellarPopulationMap(
        config,
        sbright=stellar_sbmap,
        logage=stellar_agemap,
        feh=stellar_fehmap,
        vel=stellar_velmap,
        vdisp=stellar_vdispmap,
        ebv=ebvmap
    )

**Input Parameters of** ``map2d.StellarPopulationMap``

- Simulation data configuration class: ``config``
- Two - dimensional array representing surface brightness distribution: ``sbright``, unit: :math:`\text{mag/arcsec}^2`
- Two - dimensional array representing the age distribution of the stellar population: ``logage``, unit: :math:`\log \text{yr}`
- Two - dimensional array representing the metallicity distribution of the stellar population: ``feh``, unit: :math:`\text{dex}`
- Two - dimensional array representing the velocity field: ``vel``, unit: :math:`\text{km/s}`
- Two - dimensional array representing the velocity dispersion field: ``vdisp``, unit: :math:`\text{km/s}`
- Two - dimensional array representing the dust extinction distribution: ``ebv``, unit: :math:`\text{mag}`


Two - Dimensional Distribution of Ionized Gas Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``ionizedgas``, which contains the two - dimensional distribution of ionized gas parameters, is created as follows:

.. code-block:: python

    ionizedgas = map2d.IonizedGasMap(
        config,
        halpha=gas_hamap,
        zh=gas_zhmap,
        vel=gas_velmap,
        vdisp=gas_vdispmap,
        ebv=ebvmap
    )


**Input Parameters of** ``map2d.IonizedGasMap``

- Simulation data configuration class: ``config``
- Two - dimensional array representing the flux distribution of the \(\text{H}\alpha\) emission line: ``halpha``, unit: :math:`10^{-17} \text{erg/s/cm}^2`
- Two - dimensional array representing the gas metallicity distribution: ``zh``, unit: :math:`\text{dex}`
- Two - dimensional array representing the velocity field: ``vel``, unit: :math:`\text{km/s}`
- Two - dimensional array representing the velocity dispersion field: ``vdisp``, unit: :math:`\text{km/s}`
- Two - dimensional array representing the dust extinction distribution: ``ebv``, unit: :math:`\text{mag}`


Configuration of Templates for Simulating Spectra
-----------------------------------------------------

Configure the stellar continuum template and the ionized gas template according to the methods in Sections 3.1 and 3.2.

.. code-block:: python

    stellar_tem = spec1d.StellarContinuumTemplate(config)
    gas_tem = spec1d.EmissionLineTemplate(config, model='hii')


Spectrum Simulation of the Three - Dimensional Data Cube
-----------------------------------------------------------

The simulation of the three - dimensional spectrum is completed through the ``.make_cube()`` method of the ``Cube3D`` class. The input parameters include the stellar continuum template and the ionized gas template.

An example is as follows:

.. code-block:: python

    u.make_cube(stellar_tem = stellar_tem, hii_tem = gas_tem)