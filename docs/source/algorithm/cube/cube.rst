Three-Dimensional Cube
=======================

The main task of three-dimensional cube generation is to arrange and integrate one-dimensional spectra according 
to the spatial positions provided by the two - dimensional map. In practice, we use the module ``cube3d`` to 
achieve this goal, whose working principle is as follows. Given the two-dimensional map :math:`\mathcal{M}(x, y)` 
of the galaxy, for the spatial position :math:`(x_i, y_i)`, we simulate the one-dimensional 
spectrum :math:`\mathcal{S}(\lambda)` based on the physical parameters :math:`\mathcal{M}(x_i, y_i)`. 
This spectrum is adopted as :math:`\mathcal{C}(x_i, y_i, \lambda)` in the three-dimensional spectrum of the galaxy. 
We perform the above process by traversing all spatial positions :math:`(x, y)`, thus realizing the mock of the 
three-dimensional spectrum :math:`\mathcal{C}(x, y, \lambda)` of the galaxy. Therefore, the input parameters of the ``cube3d`` 
module include the classes constructed by the ``map2d`` modules.

The spectrum of a galaxy usually contains both emission lines and the stellar continuum. Then, the two-dimensional maps 
used for the mock of the data cube should contain both ionized gas information (such as :math:`\ha` emission line flux) 
and stellar population information (such as magnitude). If the simulated galaxy does not include gas emission lines, 
such as early-type galaxies, then the two-dimensional maps only need to contain stellar population information. 
If simulating pure gas emission line sources, such as HII regions, only two-dimensional maps related to ionized gas 
need to be provided. 

The above is the mock for target sources with relatively simple structures. For some target sources with more complex 
structures, we need to decompose them into several simple sources. The data cube of each simple source is simulated 
separately and then combined to obtain the final data cube of the target source. For example, to simulate a galaxy 
with strong gas outflows, the galaxy should be decomposed into at least two parts: a normal galaxy and an outflowing 
ionized gas. Among them, the data cube of the normal galaxy is simulated according to the method introduced in the previous 
paragraph. The cube mock of the outflowing gas is equivalent to the mock of a pure emission line source, and only 
the emission lines need to be simulated. Then, by combining the two simulated cubes, the data cube of the galaxy with 
strong gas outflows is obtained.