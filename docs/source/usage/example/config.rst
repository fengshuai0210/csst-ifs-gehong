.. _simulation-data-format-configuration:

Configuration of Data Format
============================

The ``config.Config`` class defines the global configuration of spectral and spatial grids for all simulations in GEHONG. 
It determines the wavelength range, sampling step, spatial resolution, and field-of-view of simulated spectra and data cubes. 
This configuration is required by all modules, including 1D spectra, 2D maps, and 3D datacubes.

Three Usage Modes
-----------------

The class supports three configuration modes:

- ``mode='sim'``: Preset for use with ``csst-ifs-sim``, the CSST CCD image simulator.
- ``mode='etc'`` (default): Preset for use with ``csst-ifs-etc``, the exposure time calculator.
- ``mode=None``: Fully customized mode. All parameters must be specified manually.

Spectral Configuration
----------------------

The wavelength grid is defined by:

- ``wave_min`` : Starting wavelength (:math:`\mathring{\text{A}}`)
- ``wave_max`` : Ending wavelength (:math:`\mathring{\text{A}}`)
- ``dlam``     : Wavelength step size (:math:`\mathring{\text{A}}`)

Internally, these parameters define a 1D wavelength array:

.. math::

   \lambda = \text{np.arange}(\text{wave\_min}, \text{wave\_max}, \text{dlam})

Spatial Configuration
---------------------

The spatial field of view (FoV) is modeled as a regular grid of square spaxels, with:

- ``nx`` : Number of pixels along the X axis
- ``ny`` : Number of pixels along the Y axis
- ``dpix`` : Size of each spaxel (in arcseconds)

The total field of view is then:

.. math::

   \text{FoV}_x = nx \times dpix \\
   \text{FoV}_y = ny \times dpix

Preset Parameters
-----------------

The built-in presets for ``mode='etc'`` and ``mode='sim'`` are as follows:

+---------+------------+------------+--------+-----+-----+--------+
| Mode    | wave_min   | wave_max   | dlam   | nx  | ny  | dpix  |
+=========+============+============+========+=====+=====+========+
| etc     | 3500       | 10000      | 2.0    | 30  | 30  | 0.2    |
+---------+------------+------------+--------+-----+-----+--------+
| sim     | 3000       | 10500      | 1.0    | 100 | 100 | 0.1    |
+---------+------------+------------+--------+-----+-----+--------+

Usage Example
-------------

Here is an example of creating a custom configuration:

.. code-block:: python

   from gehong import config

   # Fully custom configuration
   cfg = config.Config(
       mode=None,
       wave_min=3500.0,
       wave_max=9000.0,
       dlam=1.0,
       nx=64,
       ny=64,
       dpix=0.15
   )

   print("Wavelength grid size:", len(cfg.wave))     # e.g., 5500 pixels
   print("Field of view (arcsec):", cfg.fov_x, "x", cfg.fov_y)
