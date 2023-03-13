Example
=======

Model A Spectra of AGN
----------------------

First, we model the spectra of narrow line region of AGN. 

.. code-block:: Python

    import csst-ifs-gehong.spec1d as s
    nlr_tem = s.EmissionLineTemplate(flux_table = 'nlr')