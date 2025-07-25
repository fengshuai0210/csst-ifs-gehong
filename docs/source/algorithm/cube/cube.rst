Three-Dimensional Cube
=======================

The main task of three-dimensional cube generation is to arrange and integrate one-dimensional spectra according 
to the spatial positions provided by the two-dimensional map. This is achieved via the ``cube3d`` module, which 
constructs mock IFS data cubes for galaxies and other extended sources.

Core Principle
--------------

Given a set of two-dimensional physical parameter maps :math:`\mathcal{M}(x, y)`, the ``cube3d`` module generates 
a one-dimensional spectrum :math:`\mathcal{S}(\lambda)` at each spatial position :math:`(x_i, y_i)` by calling the 
appropriate spectral synthesis modules (e.g., ``StellarContinuum``, ``HII_Region``, ``AGN``). The resulting spectrum 
is then stored as the voxel :math:`\mathcal{C}(x_i, y_i, \lambda)` in the final data cube:

.. math::

   \mathcal{C}(x, y, \lambda) = \mathcal{S}(\lambda; \mathcal{M}(x, y))

By looping over all valid pixels, the full three-dimensional mock cube :math:`\mathcal{C}(x, y, \lambda)` is assembled.

The ``cube3d`` module automatically distinguishes the type of input (stellar or gas component) based on the input map classes 
(``StellarPopulationMap``, ``IonizedGasMap``, etc.), and calls the corresponding model to simulate the spectrum. The wavelength 
grid is shared across all components and defined by the configuration object.

Component Combination
----------------------

The spectrum of a galaxy typically includes both stellar continuum and nebular emission lines. The final data cube is constructed 
by adding the simulated stellar and gas components:

.. math::

   \mathcal{C}_\text{total} = \mathcal{C}_\text{stars} + \mathcal{C}_\text{gas} + \mathcal{C}_\text{AGN}

This modular design allows each component to be simulated separately, ensuring physical consistency and interpretability. The 
relative contributions of different components can be controlled via the corresponding 2D maps (e.g., Hα flux, magnitude, star 
formation rate).

Supported Use Cases
--------------------

Depending on the input maps provided, the cube simulation supports a wide range of scenarios:

- **Stellar-only galaxies**: e.g., early-type galaxies, use only ``StellarPopulationMap``.
- **Emission-line dominated sources**: e.g., HII regions, star-forming blobs, provide only ``IonizedGasMap``.
- **Full galaxy systems**: include both ``StellarPopulationMap`` and ``IonizedGasMap`` for realistic disk galaxies.
- **AGN-hosting galaxies**: add ``AGN_PhysicalModel`` or ``AGN`` component to simulate central AGN emission.

Complex Structures and Composite Sources
----------------------------------------

For more complex systems (e.g., galaxies with outflows, mergers, or superposed sources), the simulation can be 
decomposed into multiple independent components. Each component is modeled and simulated individually, and the 
resulting cubes are combined to form the final datacube:

.. math::

   \mathcal{C}_\text{final} = \sum_i \mathcal{C}_i

This approach allows users to simulate:

- AGN outflows as an additional ionized gas cube with custom velocity and metallicity;
- Merging galaxies as two distinct sources with independent dynamics and populations;
- Background or foreground sources in the same field of view.

Each individual cube is aligned on the same spatial and spectral grid and added voxel-by-voxel. This structure provides 
maximum flexibility for constructing realistic IFS datasets.

