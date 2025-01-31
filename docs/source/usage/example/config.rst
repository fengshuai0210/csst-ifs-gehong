.. _simulation-data-format-configuration:

Configuration of Simulation Data Format
=======================================

The class ``config.config`` is used to configure the format of the simulated spectrum, 
including the configuration in the one-dimensional wavelength direction and the two-dimensional spatial direction.

**Main Input Parameters**

- Blue - end limit of wavelength: ``wave_min``, unit: :math:`\mathring{\text{A}}`
- Red - end limit of wavelength: ``wave_max``, unit: :math:`\mathring{\text{A}}`
- Wavelength interval: ``dlam``, unit: :math:`\mathring{\text{A}}`
- Number of pixels in the x - direction: ``nx``
- Number of pixels in the y - direction: ``ny``
- Pixel size: ``dpix``, unit: :math:`\text{arcsec}`

When ``config`` is configured with the following parameters:

.. code - block:: python

    from gehong import config
    config = config.config(wave_min = 3000, wave_max = 10500, dlam = 1.5, nx = 100, ny = 100, dpix = 0.1)

The simulated one-dimensional spectrum is a one-dimensional array with a wavelength ranging 
from :math:`3000\mathring{\text{A}}` to :math:`10500\mathring{\text{A}}`, a wavelength interval 
of :math:`1.5\mathring{\text{A}}`, and containing 5000 elements. 

The simulated two-dimensional image is a two-dimensional array with a field-of-view size 
of 10 arcseconds × 10 arcseconds, a single-pixel size of 0.1 arcseconds × 0.1 arcseconds, 
and a total of 100 × 100 elements. 

The simulated three-dimensional data cube is a 100 × 100 × 5000 three-dimensional data.