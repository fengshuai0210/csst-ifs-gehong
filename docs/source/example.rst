Example
=======

Model A Spectra of AGN
----------------------

Spectra of narrow line region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, preparing the class of emission line templates.

.. code-block:: Python

    import csst-ifs-gehong.spec1d as s
    nlr_tem = s.EmissionLineTemplate(flux_table = 'nlr')

Second, model the spectra of narrow line region. 

.. code-block:: Python

    wave = np.linspace(3000, 7000, 10000)
    g = AGN_NLR(wave0, nlr_tem, Halpha = 100, logZ = 0, vdisp = 200, vel = 1e3)
    
    plt.plot(g.wave,g.flux+1,lw=2,color='red')
    plt.xlim(3000,7000)
    plt.xlabel(r'$Wavelength(\AA)$')
    plt.ylabel(r'$Flux(erg/s/\AA/cm^2)$')
    plt.legend()
    plt.show()