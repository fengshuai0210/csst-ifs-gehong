import os
import glob
import sys
from os import path
import numpy as np
import astropy.units as u
from astropy.io import fits
from scipy.stats import norm
from scipy.interpolate import interp1d

data_path = os.path.dirname(__file__)

def readcol(filename, **kwargs):
    """
    Tries to reproduce the simplicity of the IDL procedure READCOL.
    Given a file with some columns of strings and columns of numbers, this
    function extract the columns from a file and places them in Numpy vectors
    with the proper type:

    name, mass = readcol('prova.txt', usecols=(0, 2))

    where the file prova.txt contains the following:

    ##################
    # name radius mass
    ##################
      abc   25.   36.
      cde   45.   56.
      rdh   55    57.
      qtr   75.   46.
      hdt   47.   56.
    ##################
    
    This function is a wrapper for numpy.genfromtxt() and accepts the same input.
    See the following website for the full documentation 
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html

    """
    f = np.genfromtxt(filename, dtype=None, **kwargs)

    t = type(f[0])
    if t == np.ndarray or t == np.void: # array or structured array
        f = map(np.array, zip(*f))

    # In Python 3.x all strings (e.g. name='NGC1023') are Unicode strings by defauls.
    # However genfromtxt() returns byte strings b'NGC1023' for non-numeric columns.
    # To have the same behaviour in Python 3 as in Python 2, I convert the Numpy
    # byte string 'S' type into Unicode strings, which behaves like normal strings.
    # With this change I can read the string a='NGC1023' from a text file and the
    # test a == 'NGC1023' will give True as expected.

    if sys.version >= '3':
        f = [v.astype(str) if v.dtype.char=='S' else v for v in f]

    return f

def log_rebin(lamRange, spec, oversample=False, velscale=None, flux=False):
    """
    Logarithmically rebin a spectrum, while rigorously conserving the flux.
    Basically the photons in the spectrum are simply redistributed according
    to a new grid of pixels, with non-uniform size in the spectral direction.
    
    When the flux keyword is set, this program performs an exact integration 
    of the original spectrum, assumed to be a step function within the 
    linearly-spaced pixels, onto the new logarithmically-spaced pixels. 
    The output was tested to agree with the analytic solution.

    :param lamRange: two elements vector containing the central wavelength
        of the first and last pixels in the spectrum, which is assumed
        to have constant wavelength scale! E.g. from the values in the
        standard FITS keywords: LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].
        It must be LAMRANGE[0] < LAMRANGE[1].
    :param spec: input spectrum.
    :param oversample: Oversampling can be done, not to loose spectral resolution,
        especally for extended wavelength ranges and to avoid aliasing.
        Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
    :param velscale: velocity scale in km/s per pixels. If this variable is
        not defined, then it will contain in output the velocity scale.
        If this variable is defined by the user it will be used
        to set the output number of pixels and wavelength scale.
    :param flux: (boolean) True to preserve total flux. In this case the
        log rebinning changes the pixels flux in proportion to their
        dLam so the following command will show large differences
        beween the spectral shape before and after LOG_REBIN:

           plt.plot(exp(logLam), specNew)  # Plot log-rebinned spectrum
           plt.plot(np.linspace(lamRange[0], lamRange[1], spec.size), spec)

        By defaul, when this is False, the above two lines produce
        two spectra that almost perfectly overlap each other.
    :return: [specNew, logLam, velscale]

    """
    lamRange = np.asarray(lamRange)
    assert len(lamRange) == 2, 'lamRange must contain two elements'
    assert lamRange[0] < lamRange[1], 'It must be lamRange[0] < lamRange[1]'
    s = spec.shape
    assert len(s) == 1, 'input spectrum must be a vector'
    n = s[0]
    if oversample:
        m = int(n*oversample)
    else:
        m = int(n)

    dLam = np.diff(lamRange)/(n - 1.)        # Assume constant dLam
    lim = lamRange/dLam + [-0.5, 0.5]        # All in units of dLam
    borders = np.linspace(*lim, num=n+1)     # Linearly
    logLim = np.log(lim)

    c = 299792.458                           # Speed of light in km/s
    if velscale is None:                     # Velocity scale is set by user
        velscale = np.diff(logLim)/m*c       # Only for output
    else:
        logScale = velscale/c
        m = int(np.diff(logLim)/logScale)    # Number of output pixels
        logLim[1] = logLim[0] + m*logScale

    newBorders = np.exp(np.linspace(*logLim, num=m+1)) # Logarithmically
    k = (newBorders - lim[0]).clip(0, n-1).astype(int)

    specNew = np.add.reduceat(spec, k)[:-1]  # Do analytic integral
    specNew *= np.diff(k) > 0    # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k])*spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:]*newBorders[:-1])*dLam)

    return specNew, logLam, velscale

def gaussian_filter1d(spec, sig):
    """
    Convolve a spectrum by a Gaussian with different sigma for every pixel.
    If all sigma are the same this routine produces the same output as
    scipy.ndimage.gaussian_filter1d, except for the border treatment.
    Here the first/last p pixels are filled with zeros.
    When creating a template library for SDSS data, this implementation
    is 60x faster than a naive for loop over pixels.

    :param spec: vector with the spectrum to convolve
    :param sig: vector of sigma values (in pixels) for every pixel
    :return: spec convolved with a Gaussian with dispersion sig

    """
    sig = sig.clip(0.01)  # forces zero sigmas to have 0.01 pixels
    p = int(np.ceil(np.max(3*sig)))
    m = 2*p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2

    n = spec.size
    a = np.zeros((m, n))
    for j in range(m):   # Loop over the small size of the kernel
        a[j, p:-p] = spec[j:n-m+j+1]

    gau = np.exp(-x2[:, None]/(2*sig**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel

    conv_spectrum = np.sum(a*gau, 0)

    return conv_spectrum

def calibrate(wave, flux, mag, filtername = 'SLOAN_SDSS.r'):
    """
    Calibrate the spectra according to the magnitude.

    :param wave: _description_
    :type wave: _type_
    :param flux: _description_
    :type flux: _type_
    :param mag: _description_
    :type mag: _type_
    :param filtername: _description_, defaults to './data/SLOAN_SDSS.r'
    :type filtername: str, optional
    :return: _description_
    :rtype: _type_
    """
        
    # Loading response curve
    if filtername == '5100':
        wave0 = np.linspace(3000,10000,7000)
        response0 = np.zeros(7000)
        response0[(wave0 > 5050) & (wave0 < 5150)] = 1.
    else:
        filter_file = data_path + '/data/filter/' + filtername+'.filter'
        wave0, response0 = readcol(filter_file)
    
    # Setting the response
    func = interp1d(wave0, response0)
    response = np.copy(wave)
    ind_extra = (wave > max(wave0)) | (wave < min(wave0))
    response[ind_extra] = 0
    ind_inside = (wave < max(wave0)) & (wave > min(wave0))
    response[ind_inside] = func(wave[ind_inside])
        
    # Flux map of datacube for given filter band
    preflux = np.sum(flux * response * np.mean(np.diff(wave))) / np.sum(response * np.mean(np.diff(wave)))
    
    # Real flux from magnitude for given filter
    realflux = (mag * u.STmag).to(u.erg/u.s/u.cm**2/u.AA).value

    # Normalization
    flux_ratio = realflux / preflux
    flux_calibrate = flux * flux_ratio * 1e17                        # Units: 10^-17 erg/s/A/cm^2
    
    return flux_calibrate

# ----------------
# Reddening Module

def Calzetti_Law(wave, Rv = 4.05):
    """
    Calzetti_Law

    Dust Extinction Curve of Calzetti et al. (2000)

    Parameters
    ----------
    wave : float
        Wavelength
    Rv : float, optional
        Extinction coefficient, by default 4.05

    Returns
    -------
    float
        Extinction value corresponding to the input wavelength
    """
    wave_number = 1./(wave * 1e-4)
    reddening_curve = np.zeros(len(wave))
    
    idx = np.logical_and(wave >= 1200, wave <= 6300)
    reddening_curve[idx] = 2.659 * ( -2.156 + 1.509 * wave_number[idx] - 0.198 * \
                    (wave_number[idx] ** 2)) + 0.011 * (wave_number[idx] **3 ) + Rv
                                    
    idx = np.logical_and(wave >= 6300, wave <= 22000)
    reddening_curve[idx] = 2.659 * ( -1.857 + 1.040 * wave_number[idx]) + Rv
    return reddening_curve

def reddening(wave, flux, ebv = 0.0, law = 'calzetti', Rv = 4.05):
    """
    Reddening an input spectra through a given reddening curve.

    Parameters
    ----------
    wave : float
        Wavelength of input spectra
    flux : float
        Flux of input spectra
    ebv : float, optional
        E(B-V) value, by default 0
    law : str, optional
        Extinction curve, by default 'calzetti'
    Rv : float, optional
        Extinction coefficient, by default 4.05

    Returns
    -------
    float
        Flux of spectra after reddening
    """
    curve = Calzetti_Law(wave, Rv = Rv)
    fluxNew = flux / (10. ** (0.4 * ebv * curve))
    return fluxNew

################
# Emission Lines
################

def SingleEmissinoLine(wave, line_wave, FWHM_inst):
    """
    Profile of single emission line (Gaussian profile)

    Parameters
    ----------
    wave : float
        Wavelength of 
    line_wave : float
        Wavelength of emission line at the line center
    FWHM_inst : float
        Intrinsic broadening of emission line, units: A

    Returns
    -------
    float
        Spectra of single emission line
    """
    sigma = FWHM_inst / 2.355 
    flux = norm.pdf(wave, line_wave, sigma)
    return flux

class EmissionLineTemplate():

    """
     Template for the emission lines
    """
    
    def __init__(self, instrument, lam_range = [500, 15000], dlam = 0.1, model = 'hii'):

        """
        __init__ _summary_

        Parameters
        ----------
        lam_range : list, optional
            _description_, by default [500, 15000]
        dlam : float, optional
            _description_, by default 0.1
        model : str, optional
            _description_, by default 'nlr'
        FWHM_inst : float, optional
            _description_, by default 0.5
        """
        
        self.lam_range = lam_range
        self.wave = np.arange(lam_range[0], lam_range[1], 0.1)
        self.FWHM_inst = instrument.inst_fwhm
        self.model = model
        
        # HII region model of fsps-cloudy
        if model == 'hii':
            
            # Loading emission line flux table
            flux_table_file = data_path + '/data/fsps.nebular.fits'
            line_table = fits.open(flux_table_file)

            # Emission line list
            line_list  = line_table[1].data
            line_wave  = line_list['Wave']
            line_names = line_list['Name']
        
            w = (line_wave > lam_range[0]) & (line_wave < lam_range[1])
            self.line_names = line_names[w]
            self.line_wave = line_wave[w]
            
            # Make parameter grid
            grid = line_table[2].data
            self.zh_grid   = grid['logZ']
        
        # Narrow line region model of 
        if model == 'nlr':
            
            # Loading emission line flux table
            flux_table_file = data_path + '/data/AGN.NLR.fits'
            line_table = fits.open(flux_table_file)

            # Emission line list
            line_list  = line_table[1].data
            line_wave  = line_list['Wave']
            line_names = line_list['Name']
        
            w = (line_wave > lam_range[0]) & (line_wave < lam_range[1])
            self.line_names = line_names[w]
            self.line_wave = line_wave[w]
            
            # Make parameter grid
            grid = line_table[2].data
            self.zh_grid   = grid['logZ']
        
        # Flux ratio
        flux_ratio = line_table[3].data
        self.flux_ratio = flux_ratio
        
        # Make emission line
        nline = len(line_wave)
        for i in range(nline):
            if i==0:
                emission_line  = SingleEmissinoLine(self.wave, line_wave[i], self.FWHM_inst)
                emission_lines = emission_line
            else:
                emission_line  = SingleEmissinoLine(self.wave, line_wave[i], self.FWHM_inst)
                emission_lines = np.vstack((emission_lines, emission_line))
        self.emission_lines = emission_lines.T

class HII_Region():

    """
    Class for the spectra of HII region

    Parameters
    ----------
    wave : _type_
        _description_
    temp : _type_
        _description_
    Halpha : int, optional
        _description_, by default 100
    logZ : bool, optional
        _description_, by default False
    vel : int, optional
        _description_, by default 100
    vdisp : int, optional
        _description_, by default 120
    Ebv : float, optional
        _description_, by default 0.1

    Raises
    ------
    ValueError
        _description_
    """
    
    def __init__(self, inst, temp, halpha = 100, zh = 0,
                 vel = 100, vdisp = 120, ebv = 0.1):     
        
        if (zh > -2) & (zh < 0.5):
            indz = np.argmin(np.abs(zh - temp.zh_grid))
            flux_ratio = temp.flux_ratio[indz, :]
        else:
            raise ValueError('The range of logZ is not correct!')
        
        # Make emission line spectra through adding emission lines                 
        emlines = temp.emission_lines * flux_ratio
        flux_combine = np.sum(emlines, axis = 1)
        flux_calibrate = flux_combine * halpha      # Units: erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(ebv):
            flux_dust = reddening(temp.wave, flux_calibrate, ebv = ebv)
            
        # Broadening caused by Velocity Dispersion
        velscale = 10
        lam_range = [np.min(temp.wave), np.max(temp.wave)]
        flux_logwave, logLam = log_rebin(lam_range, flux_dust, velscale=velscale)[:2]
        
        sigma_gas = vdisp / velscale                                # in pixel
        sigma_LSF = temp.FWHM_inst / (np.exp(logLam)) * 3e5 / velscale          # in pixel
        
        if sigma_gas>0: 
            sigma_dif = np.zeros(len(flux_logwave))
            idx = (sigma_gas > sigma_LSF)
            sigma_dif[idx] = np.sqrt(sigma_gas ** 2. - sigma_LSF[idx] ** 2.)
            idx = (sigma_gas <= sigma_LSF)
            sigma_dif[idx] = 0.1
            flux_broad = gaussian_filter1d(flux_logwave, sigma_dif)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = np.exp(logLam) * (1 + redshift)
        flux_red = np.interp(inst.wave, wave_r, flux_broad)
            
        self.wave = inst.wave
        self.flux = flux_red

#############
# AGN Spectra
#############

class AGN_NLR():

    """
    Class for narrow line region of AGN

    Parameters
    ----------
    wave : float
        Wavelength
    temp : class
        Class of emission line template
    Halpha : float, optional
        Total flux of Halpha emission line, units 1e-17 erg/s, by default 100
    logZ : float, optional
        Gas-phase metallicity, by default False
    vel : float, optional
        Line-of-sight velocity, by default 100
    vdisp : float, optional
        Velocity dispersion, by default 120
    Ebv : float, optional
        Dust extinction, by default 0.1

    Raises
    ------
    ValueError
        _description_
    """
    
    def __init__(self, temp, instrument, halpha = 100, zh = False,
                 vel = 100, vdisp = 120, ebv = 0.1):
        
        if (zh > -2) & (zh < 0.5):
            indz = np.argmin(np.abs(zh - temp.zh_grid))
            flux_ratio = temp.flux_ratio[indz, :]
        else:
            raise ValueError('The value of logZ is not correct!')
        
        # Make emission line spectra through adding emission lines                 
        emlines = temp.emission_lines * flux_ratio
        flux_combine = np.sum(emlines, axis = 1)
        flux_calibrate = flux_combine * halpha      # Units: 1e-17 erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(ebv):
            flux_dust = reddening(temp.wave, flux_calibrate, ebv = ebv)
            
        # Broadening caused by Velocity Dispersion
        velscale = 10
        lam_range = [np.min(temp.wave), np.max(temp.wave)]
        flux_logwave, logLam = log_rebin(lam_range, flux_dust, velscale=velscale)[:2]
        
        sigma_gas = vdisp / velscale                                # in pixel
        sigma_LSF = temp.FWHM_inst / (np.exp(logLam)) * 3e5 / velscale          # in pixel
        
        if sigma_gas>0: 
            sigma_dif = np.zeros(len(flux_logwave))
            idx = (sigma_gas > sigma_LSF)
            sigma_dif[idx] = np.sqrt(sigma_gas ** 2. - sigma_LSF[idx] ** 2.)
            idx = (sigma_gas <= sigma_LSF)
            sigma_dif[idx] = 0.1
            flux_broad = gaussian_filter1d(flux_logwave, sigma_dif)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = np.exp(logLam) * (1 + redshift)
        flux_red = np.interp(instrument.wave, wave_r, flux_broad)
            
        self.wave = instrument.wave
        self.flux = flux_red

class AGN_BLR():

    """
    Class for the broad line region of AGN
    """
    
    def __init__(self, instrument, hbeta_flux = 100, hbeta_fwhm = 2000, ebv = None, vel = 0.,
                 lam_range = [500, 15000]):
        
        wave_rest = np.arange(lam_range[0], lam_range[1], 0.1)
        
        line_names = ['Hepsilon', 'Hdelta', 'Hgamma', 'Hbeta', 'Halpha']
        line_waves = [3970.079, 4101.742, 4340.471, 4861.333, 6562.819]
        line_ratio = [0.101, 0.208, 0.405, 1.000, 2.579]              # From Ilic et al. (2006)
        
        # Make emission lines
        for i in range(len(line_names)):
            if i==0:
                emission_line  = SingleEmissinoLine(wave_rest, line_waves[i], 
                                                    hbeta_fwhm / 3e5 * line_waves[i])
                emission_lines = emission_line
            else:
                emission_line  = SingleEmissinoLine(wave_rest, line_waves[i], 
                                                    hbeta_fwhm / 3e5 * line_waves[i])
                emission_lines = np.vstack((emission_lines, emission_line))
        emlines = emission_lines.T * line_ratio
        flux_combine = np.sum(emlines, axis = 1)
        
        # Flux callibration
        flux_calibrate = flux_combine * hbeta_flux      # Units: 1e-17 erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = ebv)
        else:
            flux_dust = flux_calibrate
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave_rest * (1 + redshift)
        flux_red = np.interp(instrument.wave, wave_r, flux_dust)
            
        self.wave = instrument.wave
        self.flux = flux_red

class AGN_FeII():

    """
    Class for FeII emission lines of AGN
    """
    
    def __init__(self, instrument, hbeta_broad = 100, r4570 = 0.4, ebv = None, vel = 0.):
        
        """
         _summary_
        """
        filename = data_path + '/data/FeII.AGN.fits'
        
        # Loading FeII template
        hdulist = fits.open(filename)
        data = hdulist[1].data
        wave_rest  = data['WAVE']
        flux_model = data['FLUX']
        
        # Determine the flux of FeII
        Fe4570_temp  = 100
        Fe4570_model = hbeta_broad * r4570
        Ratio_Fe4570 = Fe4570_model / Fe4570_temp
        
        # Flux calibration
        flux_calibrate = flux_model * Ratio_Fe4570
        
        # Dust attenuation
        if np.isscalar(ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = ebv)
        else:
            flux_dust = flux_calibrate
             
        # Redshift
        redshift = vel / 3e5
        wave_r   = wave_rest * (1 + redshift)
        flux_red = np.interp(instrument.wave, wave_r, flux_dust)
        
        self.wave = instrument.wave
        self.flux = flux_red
        
class AGN_Powerlaw():
    """
    AGN_Powerlaw _summary_
    """
    
    def __init__(self, instrument, m5100 = 1000, alpha = -1.5, ebv = None, vel = 0.,):

        """
         _summary_
        """
        
        wave_rest = np.linspace(1000,20000,10000)
        flux      = wave_rest ** alpha
        
        # Flux calibration
        flux_calibrate = calibrate(wave_rest, flux, m5100, filtername='5100')
        
        # Dust attenuation
        if np.isscalar(ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = ebv)
        else:
            flux_dust = flux_calibrate
            
        # Redshift
        redshift = vel / 3e5
        wave_r   = wave_rest * (1 + redshift)
        flux_red = np.interp(instrument.wave, wave_r, flux_dust)
        
        self.wave = instrument.wave
        self.flux = flux_red
        
class AGN():
    
    def __init__(self, instrument, nlr_template, bhmass = 1e6, edd_ratio = 0.05, 
                 halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500, 
                 vel = 1000, zh = 0, ebv = 0.1, dist = 20):
        
        NLR = AGN_NLR(nlr_template, instrument, halpha = halpha_narrow, zh = zh,
                      vel = vel, vdisp = vdisp_narrow, ebv = ebv)
        if halpha_broad > 0:
            BLR = AGN_BLR(instrument, hbeta_flux = halpha_broad / 2.579, 
                          hbeta_fwhm = vdisp_broad / 2.355, ebv = ebv, vel = vel)
            
        m5100 = BHmass_to_M5100(bhmass, edd_ratio = edd_ratio, dist = dist)
        PL  = AGN_Powerlaw(instrument, m5100 = m5100, ebv = ebv, vel = vel)
        Fe  = AGN_FeII(instrument, hbeta_broad = halpha_broad / 2.579, ebv = ebv, vel = vel)
        
        self.wave = instrument.wave
        self.flux = NLR.flux + PL.flux + Fe.flux
        
        if halpha_broad > 0:
            self.flux = self.flux + BLR.flux
    
        
def BHmass_to_M5100(bhmass, edd_ratio = 0.05, dist = 21):
    """
    Caculate luminosity at 5100A according to the black hole mass

    Parameters
    ----------
    BHmass : float
        Black hole mass, unit: solar mass
    Edd_ratio : float, optional
        Eddtington ratio, by default 0.05
    Dist : int, optional
        Distance to the black hole, by default 21

    Returns
    -------
    float
        Luminosity at 5100A
    """
    
    # Calculate bolometric luminosity
    Ledd = 3e4 * bhmass
    Lbol = Ledd * edd_ratio

    # Convert bolometric luminosity to 5100A luminosity (Marconi et al. 2004)
    L5100 = Lbol / 10.9
    M5100 = 4.86 - 2.5 * np.log10(L5100)
    m5100 = M5100 + 5. * np.log10(dist * 1e5)

    return m5100

#################
# Stellar Spectra
#################

def age_metal(filename):
    """
    Extract the age and metallicity from the name of a file of
    the MILES library of Single Stellar Population models as
    downloaded from http://miles.iac.es/ as of 2016

    :param filename: string possibly including full path
        (e.g. 'miles_library/Mun1.30Zm0.40T03.9811.fits')
    :return: age (Gyr), [M/H]

    """
    s = path.basename(filename)
    age = float(s[s.find("T")+1:s.find("_iPp0.00_baseFe.fits")])
    metal = s[s.find("Z")+1:s.find("T")]
    if "m" in metal:
        metal = -float(metal[1:])
    elif "p" in metal:
        metal = float(metal[1:])
    else:
        raise ValueError("This is not a standard MILES filename")

    return age, metal

class StellarContinuumTemplate(object):

    def __init__(self, instrument, velscale = 50,
                 pathname = data_path + '/data/EMILES/Ech*_baseFe.fits', 
                 normalize = False):
        """
        Produces an array of logarithmically-binned templates by reading
        the spectra from the Single Stellar Population (SSP) library by
        Vazdekis et al. (2010, MNRAS, 404, 1639) http://miles.iac.es/.
        The code checks that the model specctra form a rectangular grid
        in age and metallicity and properly sorts them in both parameters.
        The code also returns the age and metallicity of each template
        by reading these parameters directly from the file names.
        The templates are broadened by a Gaussian with dispersion
        sigma_diff = np.sqrt(、**2 - sigma_tem**2).

        Thie script relies on the files naming convention adopted by
        the MILES library, where SSP spectra have the form below

            Mun1.30Zm0.40T03.9811.fits

        This code can be easily adapted by the users to deal with other stellar
        libraries, different IMFs or different abundances.

        :param pathname: path with wildcards returning the list files to use
            (e.g. 'miles_models/Mun1.30*.fits'). The files must form a Cartesian grid
            in age and metallicity and the procedure returns an error if they are not.
        :param velscale: desired velocity scale for the output templates library in km/s
            (e.g. 60). This is generally the same or an integer fraction of the velscale
            of the galaxy spectrum.
        :param FWHM_gal: vector or scalar of the FWHM of the instrumental resolution of
            the galaxy spectrum in Angstrom.
        :param normalize: set to True to normalize each template to mean=1.
            This is useful to compute light-weighted stellar population quantities.
        :return: The following variables are stored as attributes of the miles class:
            .templates: array has dimensions templates[npixels, n_ages, n_metals];
            .log_lam_temp: natural np.log() wavelength of every pixel npixels;
            .age_grid: (Gyr) has dimensions age_grid[n_ages, n_metals];
            .metal_grid: [M/H] has dimensions metal_grid[n_ages, n_metals].
            .n_ages: number of different ages
            .n_metal: number of different metallicities
        """
        FWHM_inst = instrument.inst_fwhm
        
        files = glob.glob(pathname)
        assert len(files) > 0, "Files not found %s" % pathname

        all = [age_metal(f) for f in files]
        all_ages, all_metals = np.array(all).T
        ages, metals = np.unique(all_ages), np.unique(all_metals)
        n_ages, n_metal = len(ages), len(metals)

        assert set(all) == set([(a, b) for a in ages for b in metals]), \
            'Ages and Metals do not form a Cartesian grid'

        # Extract the wavelength range and logarithmically rebin one spectrum
        # to the same velocity scale of the SDSS galaxy spectrum, to determine
        # the size needed for the array which will contain the template spectra.
        hdu = fits.open(files[0])
        ssp = hdu[0].data
        h2 = hdu[0].header
        lam_range_temp = h2['CRVAL1'] + np.array([0, h2['CDELT1']*(h2['NAXIS1']-1)])
        sspNew, log_lam_temp = log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]
        #wave=((np.arange(hdr['NAXIS1'])+1.0)-hdr['CRPIX1'])*hdr['CDELT1']+hdr['CRVAL1']
        
        templates = np.empty((sspNew.size, n_ages, n_metal))
        age_grid = np.empty((n_ages, n_metal))
        metal_grid = np.empty((n_ages, n_metal))

        # Convolve the whole Vazdekis library of spectral templates
        # with the quadratic difference between the galaxy and the
        # Vazdekis instrumental resolution. Logarithmically rebin
        # and store each template as a column in the array TEMPLATES.

        # Quadratic sigma difference in pixels Vazdekis --> galaxy
        # The formula below is rigorously valid if the shapes of the
        # instrumental spectral profiles are well approximated by Gaussians.
        
        # FWHM of Emiles templates
        Emile_wave = np.exp(log_lam_temp)
        Emile_FWHM = np.zeros(h2['NAXIS1'])
        Emile_FWHM[np.where(Emile_wave < 3060)] = 3.
        Emile_FWHM[np.where((Emile_wave >= 3060) & (Emile_wave < 3540))] = 3.
        Emile_FWHM[np.where((Emile_wave >= 3540) & (Emile_wave < 8950))] = 2.5
        Lwave = Emile_wave[np.where(Emile_wave >= 8950)]
        Emile_FWHM[np.where(Emile_wave >= 8950)]=60*2.35/3.e5*Lwave  # sigma=60km/s at lambda > 8950
        
        LSF = Emile_FWHM
        
        FWHM_eff=Emile_FWHM.copy()   # combined FWHM from stellar library and instrument(input)
        if np.isscalar(FWHM_inst):
            FWHM_eff[Emile_FWHM < FWHM_inst] = FWHM_inst
            LSF[Emile_FWHM < FWHM_inst] = FWHM_inst
        #else:
        #    FWHM_eff[Emile_FWHM < FWHM_inst] = FWHM_inst[Emile_FWHM < FWHM_inst]
        #    LSF[Emile_FWHM < FWHM_inst] = FWHM_inst[Emile_FWHM < FWHM_inst]
        FWHM_dif  = np.sqrt(FWHM_eff**2 - Emile_FWHM**2)
        sigma_dif = FWHM_dif/2.355/h2['CDELT1']   # Sigma difference in pixels

        # Here we make sure the spectra are sorted in both [M/H] and Age
        # along the two axes of the rectangular grid of templates.
        for j, age in enumerate(ages):
            for k, metal in enumerate(metals):
                p = all.index((age, metal))
                hdu = fits.open(files[p])
                ssp = hdu[0].data
                #if np.isscalar(FWHM_dif):
                #    ssp = ndimage.gaussian_filter1d(ssp, sigma_dif)
                #else:
                ssp = gaussian_filter1d(ssp, sigma_dif)  # convolution with variable sigma
                sspNew  = log_rebin(lam_range_temp, ssp, velscale=velscale)[0]
                #if normalize:
                #    sspNew /= np.mean(sspNew)
                templates[:, j, k] = sspNew
                age_grid[j, k]   = age
                metal_grid[j, k] = metal

        self.templates = templates/np.median(templates)  # Normalize by a scalar
        self.log_lam_temp = log_lam_temp
        self.age_grid = age_grid
        self.metal_grid = metal_grid
        self.n_ages = n_ages
        self.n_metal = n_metal
        self.lsf = log_rebin(lam_range_temp, LSF, velscale=velscale)[0]
        self.velscale = velscale
        
    def fmass_ssp(self): 

        isedpath = data_path + '/data/EMILES/model/'
        massfile = isedpath + 'out_mass_CH_PADOVA00'
        n_metal = self.n_metal
        n_ages = self.n_ages
        fage = self.age_grid[:,0]

        fMs = np.zeros((n_ages,n_metal))
        
        Metal, Age, Ms = readcol(massfile, usecols=(2, 3, 6))
        for i in range(n_metal):
            for j in range(self.n_ages):
                locmin = np.argmin(abs(Metal - self.metal_grid[j, i])) & np.argmin(abs(Age - self.age_grid[j, i]))
                fMs[j,i] = Ms[locmin]

        return fMs
    
class StellarContinuum():
    
    def __init__(self, template, instrument, mag = 15, age = 1, feh = 0, vel = 100, vdisp = 100, ebv = 0):

        """
        Modelling the spectra of stellar population
        
        Pars:
            ssp        - The class of simple stellar population
            obswave    - Wavelength, 1D-array
            mag        - Magnitude in r-band used for the calibration of spectra, scalar
            Age        - Mean age of stellar population used for determining the spectra profile, scalar (Gyr)
            FeH        - Mean FeH of stellar population used for determining the spectra profile, scalar 
            vel        - Line of sight velocity used for determining the doppler effect, scalar (km/s)
            vdisp      - Line of sight velocity dispersion used for broadening the spectra, scalar (km/s)
            Ebv        - Extinction, scalar (mag), default 0 (no dust extinction)
            dFeH       - Range of FeH when searching the best templates, default 0.2 (searching the template
                         with -0.2 < FeH - dFeh < 0.2)
                         
        Attri:
            wave       - Wavelength is the same with obswave in Pars
            flux       - Flux of best model
        """
        
        # -----------------
        # Stellar Continuum
        
        SSP_temp = template.templates

        # Select metal bins
        metals = template.metal_grid[0,:]
        minloc = np.argmin(abs(feh - metals))
        tpls = SSP_temp[:, :, minloc]
        fmass = template.fmass_ssp()[:, minloc]
        
        # Select age bins
        Ages = template.age_grid[:,0]
        minloc = np.argmin(abs(age-Ages))
        Stellar = tpls[:, minloc]
        
        wave = np.exp(template.log_lam_temp)
        
        # Broadening caused by Velocity Dispersion
        sigma_gal = vdisp / template.velscale                   # in pixel
        sigma_LSF = template.lsf / template.velscale                 # in pixel
        
        if sigma_gal>0: 
            sigma_dif = np.zeros(len(Stellar))
            idx = (sigma_gal > sigma_LSF)
            sigma_dif[idx] = np.sqrt(sigma_gal ** 2. - sigma_LSF[idx] ** 2.)
            idx = (sigma_gal <= sigma_LSF)
            sigma_dif[idx] = 0.1
            flux0 = gaussian_filter1d(Stellar, sigma_dif)
        
        # Dust Reddening
        if np.isscalar(ebv):
            flux0 = reddening(wave, flux0, ebv = ebv)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave * (1 + redshift)
        
        flux = np.interp(instrument.wave, wave_r, flux0)
        
        # Calibration
        if np.isscalar(mag):
            flux = calibrate(instrument.wave, flux, mag, filtername='SLOAN_SDSS.r')
        
        # Convert to input wavelength
        self.wave = instrument.wave
        self.flux = flux
        
#####################
# Single Star Spectra
#####################

class SingleStarTemplate():
        
    def __init__(self, inst, velscale = 20):
        
        filename = data_path + '/data/Starlib.XSL.fits'
        fwhm_inst = inst.inst_fwhm
        
        hdulist = fits.open(filename)
        lam  = hdulist[1].data['Wave']
        flux = hdulist[2].data
        par  = hdulist[3].data
        
        lam_range_temp = np.array([3500, 12000])
        TemNew, log_lam_temp = log_rebin(lam_range_temp, flux[1, :], velscale = velscale)[:2]
        
        # FWHM of XLS templates
        Temp_wave = np.exp(log_lam_temp)
        Temp_FWHM = np.zeros(len(log_lam_temp))
        Temp_FWHM[(Temp_wave < 5330)] = 13 * 2.35 / 3e5 * Temp_wave[(Temp_wave < 5330)]   # sigma = 13km/s at lambda <5330
        Temp_FWHM[(Temp_wave >= 5330) & (Temp_wave < 9440)] = 11 * 2.35 / 3e5 * Temp_wave[(Temp_wave >= 5330) & (Temp_wave < 9440)]
        # sigma = 13km/s at 5330 < lambda < 9440
        Temp_FWHM[(Temp_wave >= 9440)] = 16 * 2.35 / 3e5 * Temp_wave[(Temp_wave >= 9440)] # sigma=16km/s at lambda > 9440
        
        LSF = Temp_FWHM
        
        FWHM_eff = Temp_FWHM.copy()   # combined FWHM from stellar library and instrument(input)
        if np.isscalar(fwhm_inst):
            FWHM_eff[Temp_FWHM < fwhm_inst] = fwhm_inst
            LSF[Temp_FWHM < fwhm_inst]      = fwhm_inst
        #else:
        #    FWHM_eff[Temp_FWHM < FWHM_inst] = FWHM_inst[Temp_FWHM < FWHM_inst]
        #    LSF[Temp_FWHM < FWHM_inst]      = FWHM_inst[Temp_FWHM < FWHM_inst]
        FWHM_dif  = np.sqrt(FWHM_eff ** 2 - Temp_FWHM ** 2)
        sigma_dif = FWHM_dif / 2.355 / (lam[1] - lam[0])  # Sigma difference in pixels

        temp = np.empty((TemNew.size, par.size))
        for i in range(par.size):
            temp0 = log_rebin(lam_range_temp, flux[i, :], velscale=velscale)[0]
            #if np.isscalar(FWHM_dif):
            #    temp1 = ndimage.gaussian_filter1d(temp0, sigma_dif)
            #else:
            temp1 = gaussian_filter1d(temp0, sigma_dif)             # convolution with variable sigma
            tempNew = temp1 / np.mean(temp1)
            temp[:, i] = tempNew
            

        self.templates    = temp
        self.log_lam_temp = log_lam_temp
        self.teff_grid = par['Teff']
        self.feh_grid  = par['FeH']
        self.logg_grid = par['logg']
        self.lsf       = Temp_FWHM
        self.velscale  = velscale
        
class SingleStar():
    
    def __init__(self, template, instrument, mag = 15, teff = 10000, feh = 0, vel = 100, ebv = 0):
 
        StarTemp = template.templates
        
        # Select metal bins
        idx_feh = (np.abs(template.feh_grid - feh) < 0.5)
        tpls = StarTemp[:, idx_feh]
        
        # Select Teff bins
        teff_feH = template.teff_grid[idx_feh]
        minloc   = np.argmin(abs(teff - teff_feH))
        starspec = tpls[:, minloc]
        
        wave = np.exp(template.log_lam_temp)
        
        # Dust Reddening
        if np.isscalar(ebv):
            starspec = reddening(wave, starspec, ebv = ebv)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave * (1 + redshift)
        
        flux = np.interp(instrument.wave, wave_r, starspec)
        
        # Calibration
        if np.isscalar(mag):
            flux = calibrate(instrument.wave, flux, mag, filtername='SLOAN_SDSS.r')
        
        # Convert to input wavelength
        self.wave = instrument.wave
        self.flux = flux