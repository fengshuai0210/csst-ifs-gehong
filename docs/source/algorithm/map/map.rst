Two-Dimensional Map
=====================

To mock the IFS data of extended sources (e.g. galaxies), we need to construct a series of 2D maps of physical parameters. 
For the case of galaxies, these 2D maps can be divided into two main categories. 

The first category is used for describing the properties of the stellar components in galaxies, including the input parameters 
required to simulate the stellar continuum spectra in Section :ref:`algorithm_stellar_contimuum`, such as the distribution of surface 
brightness, the age distribution of stellar populations, and the metallicity distribution. These 2D maps are integrated 
together through the ``map2d.StellarPopulationMap`` module, whose input parameters are the two-dimensional array corresponding 
to each 2D map.

The second category is related to the properties of the gas components in galaxies, including the input parameters required 
to simulate emission lines of ionized gas in Section :ref:`ionized-gas-emission-lines`, such as the flux distribution of :math:`\text{H}\alpha` emission 
line and the distribution of gas-phase metallicity. These 2D maps related to gas components are integrated through 
the ``map2d.IonizedGasMap`` module, whose input parameters are also the two-dimensional array. 

To be able to mock IFS data that is closer to real observations, the 2D maps of physical parameters should have rich detailed 
structures, as shown in Figure :ref:`fig:case_map2d`. Therefore, it is the best choice to construct 2D maps based on real 
high-resolution observation data or numerical simulation data. As a supplement, we also provide three commonly used parametric 
models to construct the 2D maps. Those maps can be used as the two-dimensional distribution of physical parameters.

Sersic Model
----------------

The surface brightness distribution of galaxies can be simulated using a two-dimensional Sersic 
model :cite:`Sersic1963, Graham2005`. The generation of the Sersic 2D map is implemented by the 
module ``map2d.sersic_map``. 

The Sersic model is represented as follows:

.. math::

   I(x, y) = I_e \exp\left \{-b_n\left [\left (\frac {R(x, y)}{R_e}\right)^{1/n}-1\right]\right \}

In this equation, :math:`R_e` is the effective radius (``reff``) which is equal to the half-light radius, 
:math:`I_e` is the flux at the half-light radius or :math:`I_e = I(R_e)`, and the :math:`n` is the Sérsic 
index (``n``) which determines the `steepness' of the profile. :math:`n = 4` corresponds to a de Vaucouleurs’ 
profile (elliptical galaxies) while :math:`n = 1` gives the exponential profile (disk galaxies). 
The :math:`b_n` is a constant that is related to the value of :math:`n`:

.. math::

   \Gamma(2n)=2\gamma(2n, b_n)

:math:`R(x, y)` is the radius from the center that corresponds to :math:`(x, y)`. 
For inclined galaxies, the radius :math:`R(x, y)` from the center :math:`(x_0, y_0)` to :math:`(x, y)` 
is given by the following expression:

.. math::

   R(x, y) = \sqrt{R_\text{maj}^2(x, y) - \left (\frac{R_\text{min}(x, y)}{1-q}\right )^2}

and

.. math::

   R_\text{maj}(x, y) & = (x - x_0) \cos \theta + (y - y_0) \sin \theta \\
   R_\text{min}(x, y) &= -(x - x_0) \sin \theta + (y - y_0) \cos \theta

where :math:`\theta` is the position angle (``pa``), :math:`q` is the ellipticity of galaxies (``ellip``). 
In practice, we use the total magnitude of galaxies :math:`m_\text{tot}` (``mag``) to calibrate the surface 
brightness of galaxies. The corresponding total luminosity :math:`L_\text{tot}` of a galaxy that follows the 
Sérsic model can be expressed as :cite:`Ciotti1991`:

.. math::

   L_\text{tot} = 2\pi \int_0^{\infty}I(R)R \text{d} R = I_e R_e^2 2\pi n \frac {e^{b_n}}{(b_n)^{2n}} \Gamma(2n)


tanh Model
--------------

The simplest velocity map is axisymmetric, in which the rotation curve of a galaxy is radially symmetric. 
We provided the module ``map2d.tanh_map`` to generate such a velocity map. 

For the galaxies with the inclination angle of :math:`i \approx \arccos (1-q)`, the expression for the 
velocity at :math:`(x, y)` in the galaxy plane is given by:

.. math::

   V(x,y) = V(R, \theta) = V_\text{sys} + V_c(R)\,\cos\phi\,\sin i\,.

where the :math:`R` is the radius from the galaxy center, :math:`\phi` is the azimuthal angle in the galaxy plane, 
:math:`V_\text{sys}` is the systematic recession velocity of a galaxy, and :math:`V_c(R)` is the intrinsic rotational 
velocity at the radius of :math:`R` :cite:`vanderKruit1978`. The relationship between :math:`(x, y)` and 
:math:`(R, \theta)` is represented by Equation :ref:`eq:Radius`, and :math:`\phi` is expressed as:

.. math::

   \cos \phi = \frac{-(x-x_0)\,\sin \theta+(y-y_0)\,\cos \theta}{R(x, y)}

where the :math:`\theta` is the position angle of the galaxy. The rotational velocity at :math:`R` is defined as:

.. math::

   V_c(R) = V_\text{max}\tanh(R/R_t)

where :math:`V_\text{max}` is the maximum rotational velocity (``vmax``), and :math:`R_t` is the turnover 
radius (``rt``) beyond which the rotation curve becomes flat :cite:`Andersen203`.

Gredient Model
--------------

Except for surface brightness and velocity field, the 2D maps of other parameters of galaxies, such as the age and 
metallicity of stellar populations :cite:`Koleva2011, SanchezBlazquez2014`, can approximately be described by a 
simple gradient model (right panel in Figure :ref:`fig:case_map2d`). We use the module ``map2d.gred_map`` to 
construct this gradient model. The input parameters related to the case are listed in Table :ref:`tab:map2d`.

For a galaxies with an inclination angle :math:`i`, the physical parameters :math:`A` at any position :math:`(x, y)` 
can be described as follows:

.. math::

   A(x,y) = A(R, \theta) =  A_\text{eff} + \nabla_A \log \frac{R(x, y)}{R_e}.

where :math:`\nabla_A` is the gradient of parameter :math:`A` (``gred``), :math:`A_e` is the value at the effective 
radius (``aeff``). 

