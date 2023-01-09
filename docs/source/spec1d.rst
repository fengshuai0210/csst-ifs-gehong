Spec1D
=====

Emission Lines of Ionized Gas
----------------

First, you can use the ``gehong.EmissionLineTemplate()`` function to load the information of emission lines, such as the line list. 

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

Stellar Continuum of Galaxies
----------------

Spectrum of Single Stellar
----------------

Spectrum of AGN Components
----------------