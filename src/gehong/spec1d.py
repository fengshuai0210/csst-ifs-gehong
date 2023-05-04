import os
import glob
import sys
from os import path
import numpy as np
import astropy.units as u
from astropy.io import fits
from scipy.stats import norm
from scipy.interpolate import interp1d

data_path = os.getenv('GEHONG_DATA_PATH')

def readcol(filename, **kwargs):
    """
    readcol, taken from ppxf.

    Parameters
    ----------
    filename : string
        The name of input ascii file

    Returns
    -------
    float
        The value of each columns
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
    This function is taken from ppxf. 

    Parameters
    ----------
    lamRange : array
        Two elements vector containing the central wavelength
        of the first and last pixels in the spectrum
    spec : array
        Input spectrum
    oversample : bool, optional
        Oversampling can be done, not to loose spectral resolution,
        especally for extended wavelength ranges and to avoid aliasing, by default False
    velscale : float, optional
        velocity scale in km/s per pixels, by default None
    flux : bool, optional
        True to preserve total flux, by default False

    Returns
    -------
    specNew : array
        Output spectrum
    logLam : array
        Wavelength array in logarithm
    velscale : array
        velocity scale in km/s per pixels
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
    One-dimensional Gaussian convolution

    Parameters
    ----------
    spec : float array
        vector with the spectrum to convolve
    sig : float
        vector of sigma values (in pixels) for every pixel
    Returns
    -------
    float array
        Spectrum after convolution
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
    Flux calibration of spectrum

    Parameters
    ----------
    wave : float array
        Wavelength of input spectrum
    flux : float array
        Flux of input spectrum
    mag : float
        Magnitude used for flux calibration
    filtername : str, optional
        Filter band name, by default 'SLOAN_SDSS.r'

    Returns
    -------
    float array
        Spectrum after flux calibration
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
    Dust Extinction Curve of Calzetti et al. (2000)

    Parameters
    ----------
    wave : float, or float array
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
    
    idx = (wave >= 1200) & (wave < 6300)
    reddening_curve[idx] = 2.659 * ( -2.156 + 1.509 * wave_number[idx] - 0.198 * \
                    (wave_number[idx] ** 2)) + 0.011 * (wave_number[idx] **3 ) + Rv
                                    
    idx = (wave >= 6300) & (wave <= 22000)
    reddening_curve[idx] = 2.659 * ( -1.857 + 1.040 * wave_number[idx]) + Rv
    return reddening_curve

def reddening(wave, flux, ebv = 0.0, law = 'calzetti', Rv = 4.05):
    """
    Reddening an input spectra through a given reddening curve.

    Parameters
    ----------
    wave : float array
        Wavelength of input spectra
    flux : float array
        Flux of input spectra
    ebv : float, optional
        E(B-V) value, by default 0
    law : str, optional
        Extinction curve, by default 'calzetti'
    Rv : float, optional
        Extinction coefficient, by default 4.05

    Returns
    -------
    float array
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
    wave : float array
        Wavelength of spectrum
    line_wave : float
        Wavelength of emission line at the line center
    FWHM_inst : float
        Intrinsic broadening of emission line, units: A

    Returns
    -------
    float array
        Spectra of single emission line
    """
    sigma = FWHM_inst / 2.355 
    flux = norm.pdf(wave, line_wave, sigma)
    return flux

class EmissionLineTemplate():
    """
    Template for the emission lines

    Parameters
    ----------
    config : class
        The class of configuration.
    lam_range : list, optional
        Wavelength range, by default [500, 15000]
    dlam : float, optional
        Wavelength width per pixel, by default 0.1A
    model : str, optional
        Emission line model, including 'hii' for HII region and 'nlr' for narrow line region of AGN, 
        by default 'hii'
    """
    
    def __init__(self, config, lam_range = [500, 15000], dlam = 0.1, model = 'hii'):
        
        self.lam_range = lam_range
        self.wave = np.arange(lam_range[0], lam_range[1], 0.1)
        self.FWHM_inst = config.inst_fwhm
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
            self.logz_grid   = grid['logZ']
        
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
            self.logz_grid   = grid['logZ']
        
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
    config : class
        Class of configuration
    temp : class
        Class of emission line template
    halpha : float, optional
        Integral flux of Halpha emission line, by default 100 * 1e-17 erg/s/cm^2
    logz : float, optional
        Gas-phase metallicity, by default 0.0
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    vdisp : float, optional
        Velocity dispersion, by default 120.0km/s
    ebv : float, optional
        Dust extinction, by default 0.1

    Raises
    ------
    ValueError
        The value of logZ should be between -2 and 0.5.
    """
    
    def __init__(self, config, temp, halpha = 100.0, logz = 0.0,
                 vel = 100.0, vdisp = 120.0, ebv = 0.1):

        if (logz > -2) & (logz < 0.5):
            indz = np.argmin(np.abs(logz - temp.logz_grid))
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
        else:
            flux_broad = flux_logwave
            
        # Redshift
        redshift = vel / 3e5
        wave_r = np.exp(logLam) * (1 + redshift)
        flux_red = np.interp(config.wave, wave_r, flux_broad)
            
        self.wave = config.wave
        self.flux = flux_red

#############
# AGN Spectra
#############

class AGN_NLR():

    """
    Class for narrow line region of AGN

    Parameters
    ----------
    config : class
        Class of configuration
    temp : class
        Class of emission line template
    halpha : float, optional
        Integral flux of Halpha emission line, by default 100 * 1e-17 erg/s/cm^2
    logz : float, optional
        Gas-phase metallicity, by default 0.0
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    vdisp : float, optional
        Velocity dispersion, by default 120.0km/s
    ebv : float, optional
        Dust extinction, by default 0.1

    Raises
    ------
    ValueError
        The value of logZ should be between -2 and 0.5.
    """
    
    def __init__(self, config, temp, halpha = 100.0, logz = 0.0,
                 vel = 100.0, vdisp = 120.0, ebv = 0.1):
        
        if (logz > -2) & (logz < 0.5):
            indz = np.argmin(np.abs(logz - temp.logz_grid))
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
        flux_red = np.interp(config.wave, wave_r, flux_broad)
            
        self.wave = config.wave
        self.flux = flux_red

class AGN_BLR():

    """
    Class for the broad line region of AGN

    Parameters
    ----------
    config : class
        Class of configuration
    hbeta_flux : float, optional
        Integral flux of Hbeta broad line, by default 100 * 1e-17 erg/s/cm^2
    hbeta_fwhm : float, optional
        FWHM of Hbeta broad line, by default 2000.0km/s
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    ebv : float, optional
        Dust extinction, by default 0.1
    lam_range : list, optional
        Wavelength range, by default [500, 15000]
    """

    def __init__(self, config, hbeta_flux = 100.0, hbeta_fwhm = 2000.0, ebv = 0.1, 
                 vel = 0.,lam_range = [500, 15000]):

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
        flux_red = np.interp(config.wave, wave_r, flux_dust)
            
        self.wave = config.wave
        self.flux = flux_red

class AGN_FeII():
    """
    Class for FeII emission lines of AGN

    Parameters
    ----------
    config : class
        Class of configuration
    hbeta_broad : float, optional
        Integral flux of Hbeta broad line, by default 100 * 1e-17 erg/s/cm^2
    r4570 : float, optional
        Flux ratio between Fe4570 flux and Hbeta broad line flux, by default 0.4
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    ebv : float, optional
        Dust extinction, by default 0.1
    """
    
    def __init__(self, config, hbeta_broad = 100.0, r4570 = 0.4, ebv = 0.1, vel = 100.0):
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
        flux_red = np.interp(config.wave, wave_r, flux_dust)
        
        self.wave = config.wave
        self.flux = flux_red
        
class AGN_Powerlaw():

    """
    The class of power-law spectrum of AGN

    Parameters
    ----------
    config : class
        Class of configuration
    M5100 : float, optional
        Magnitude of power law spectrum between 5050A and 5150A, by default 1000.0 * 1e-17 erg/s/cm^2
    alpha : float, optional
        Index of power law, by default -1.5
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    Ebv : float, optional
        Dust extinction, by default 0.1
    """
    
    def __init__(self, config, m5100 = 1000.0, alpha = -1.5, vel = 100.0, ebv = 0.1):

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
        flux_red = np.interp(config.wave, wave_r, flux_dust)
        
        self.wave = config.wave
        self.flux = flux_red
        
class AGN():

    """
    Class of singal spectra of AGN

    Parameters
    ----------
    config : class
        Class of configuration
    nlr_template : class
        Class of emission line template
    bhmass : float, optional
        Black hole mass used for calculating the luminosity of power law spectrum at 5100A, 
        by default 1e6 solar mass
    edd_ratio : float, optional
        Eddinton ratio used for calculating the luminosity of power law spectrum at 5100A, by default 0.05
    halpha_broad : float, optional
        Integral flux of Halpha broad line, by default 100.0 * 1e-17 erg/s/cm^2
    halpha_narrow : float, optional
        Integral flux of Halpha narrow line, by default 100.0 * 1e-17 erg/s/cm^2
    vdisp_broad : float, optional
        Velocity dispersion of Halpha broad line, by default 5000.0km/s
    vdisp_narrow : float, optional
        Velocity dispersion of Halpha narrow line, by default 500.0km/s
    vel : float, optional
        Line of sight velocity, by default 1000.0km/s
    logz : float, optional
        Gas-phase metallicity of narrow line region, by default 0.0
    ebv : float, optional
        Dust extinction, by default 0.1
    dist : float, optional
        Luminosity distance of AGN, by default 20.0Mpc
    """
    def __init__(self, config, nlr_template, bhmass = 1e6, edd_ratio = 0.05, 
                 halpha_broad = 100.0, halpha_narrow = 100.0, vdisp_broad = 5000.0, vdisp_narrow = 500.0, 
                 vel = 1000.0, logz = 0.0, ebv = 0.1, dist = 20.0):
        
        NLR = AGN_NLR(config, nlr_template, halpha = halpha_narrow, logz = logz,
                      vel = vel, vdisp = vdisp_narrow, ebv = ebv)
        if halpha_broad > 0:
            BLR = AGN_BLR(config, hbeta_flux = halpha_broad / 2.579, 
                          hbeta_fwhm = vdisp_broad / 2.355, ebv = ebv, vel = vel)

        m5100 = BHmass_to_M5100(bhmass, edd_ratio = edd_ratio, dist = dist)
        PL  = AGN_Powerlaw(config, m5100 = m5100, ebv = ebv, vel = vel)
        Fe  = AGN_FeII(config, hbeta_broad = halpha_broad / 2.579, ebv = ebv, vel = vel)
        
        self.wave = config.wave
        self.flux = NLR.flux + PL.flux + Fe.flux
        
        if halpha_broad > 0:
            self.flux = self.flux + BLR.flux
    
def BHmass_to_M5100(bhmass, edd_ratio = 0.05, dist = 21.0):
    """
    Caculate magnitude at 5100A according to the black hole mass

    Parameters
    ----------
    bhmass : float
        Black hole mass, unit: solar mass
    edd_ratio : float, optional
        Eddtington ratio, by default 0.05
    dist : float, optional
        Distance to the black hole, by default 21.0Mpc

    Returns
    -------
    float
        Magnitude at 5100A
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

"""
    Extract the age and metallicity from the name of a file of
    the MILES library of Single Stellar Population models as
    downloaded from http://miles.iac.es/ as of 2016

    :param filename: string possibly including full path
        (e.g. 'miles_library/Mun1.30Zm0.40T03.9811.fits')
    :return: age (Gyr), [M/H]

    """

def age_metal(filename):
    """
    Extract the age and metallicity from the name of a file of
    the MILES library of Single Stellar Population models.

    Parameters
    ----------
    filename : string
        Full path of template files

    Returns
    -------
    age : float
        Age of SSP (Gyr)
    FeH : float
        Metallicity of SSP

    Raises
    ------
    ValueError
        This is not a standard MILES filename
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
    """
    Class of single stellar population template. 

    Parameters
    ----------
    config : class
        Class of configuration
    velscale : array
        velocity scale in km/s per pixels, by default 50.0km/s
    pathname : string, optional
        path with wildcards returning the list files to use, 
        by default data_path+'/data/EMILES/Ech*_baseFe.fits'
    normalize : bool, optional
        Set to True to normalize each template to mean=1, by default False
    """
    def __init__(self, config, velscale = 50,
                 pathname = data_path + '/data/EMILES/Ech*_baseFe.fits', 
                 normalize = False):

        FWHM_inst = config.inst_fwhm
        
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
        else:
            FWHM_eff[Emile_FWHM < FWHM_inst] = FWHM_inst[Emile_FWHM < FWHM_inst]
            LSF[Emile_FWHM < FWHM_inst] = FWHM_inst[Emile_FWHM < FWHM_inst]
        FWHM_dif  = np.sqrt(FWHM_eff**2 - Emile_FWHM**2)
        sigma_dif = FWHM_dif/2.355/h2['CDELT1']   # Sigma difference in pixels

        # Here we make sure the spectra are sorted in both [M/H] and Age
        # along the two axes of the rectangular grid of templates.
        for j, age in enumerate(ages):
            for k, metal in enumerate(metals):
                p = all.index((age, metal))
                hdu = fits.open(files[p])
                ssp = hdu[0].data
                if np.isscalar(FWHM_dif):
                    ssp = ndimage.gaussian_filter1d(ssp, sigma_dif)
                else:
                    ssp = gaussian_filter1d(ssp, sigma_dif)  # convolution with variable sigma
                sspNew  = log_rebin(lam_range_temp, ssp, velscale=velscale)[0]
                if normalize:
                    sspNew /= np.mean(sspNew)
                templates[:, j, k] = sspNew
                age_grid[j, k]   = age
                metal_grid[j, k] = metal

        self.templates = templates/np.median(templates)  # Normalize by a scalar
        self.log_lam_temp = log_lam_temp
        self.age_grid = age_grid
        self.metal_grid = metal_grid
        self.n_ages = n_ages
        self.n_metal = n_metal
        self.LSF = log_rebin(lam_range_temp, LSF, velscale=velscale)[0]
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
    """
    The class of stellar continuum
    
    Parameters
    ----------
    config : class
        Class of configuration
    template : class
        Class of single stellar population template
    mag : float, optional
        Magnitude in SDSS r-band, by default 15.0
    age : float, optional
        Median age of stellar continuum, by default 1.0Gyr
    feh : float, optional
        Metallicity of stellar continuum, by default 0.0
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    vdisp : float, optional
        Velocity dispersion, by default 120.0km/s
    ebv : float, optional
        Dust extinction, by default 0.1
    """
    def __init__(self, config, template, mag = 15.0, age = 1.0, feh = 0.0, 
                 vel = 100.0, vdisp = 100.0, ebv = 0.1):

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
        sigma_LSF = template.LSF / template.velscale                 # in pixel
        
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
        
        flux = np.interp(config.wave, wave_r, flux0)
        
        # Calibration
        if np.isscalar(mag):
            flux = calibrate(config.wave, flux, mag, filtername='SLOAN_SDSS.r')
        
        # Convert to input wavelength
        self.wave = config.wave
        self.flux = flux
        
#####################
# Single Star Spectra
#####################

class SingleStarTemplate():
    """
    Class of single stellar template

    Parameters
    ----------
    config : class
        Class of configuration
    velscale : float, option
        velocity scale in km/s per pixels, by default 20.0km/s
    """
    def __init__(self, config, velscale = 20):

        FWHM_inst = config.inst_fwhm
        filename = data_path + '/data/Starlib.XSL.fits'
        
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
        if np.isscalar(FWHM_inst):
            FWHM_eff[Temp_FWHM < FWHM_inst] = FWHM_inst
            LSF[Temp_FWHM < FWHM_inst]      = FWHM_inst
        else:
            FWHM_eff[Temp_FWHM < FWHM_inst] = FWHM_inst[Temp_FWHM < FWHM_inst]
            LSF[Temp_FWHM < FWHM_inst]      = FWHM_inst[Temp_FWHM < FWHM_inst]
        FWHM_dif  = np.sqrt(FWHM_eff ** 2 - Temp_FWHM ** 2)
        sigma_dif = FWHM_dif / 2.355 / (lam[1] - lam[0])  # Sigma difference in pixels

        temp = np.empty((TemNew.size, par.size))
        for i in range(par.size):
            temp0 = log_rebin(lam_range_temp, flux[i, :], velscale=velscale)[0]
            if np.isscalar(FWHM_dif):
                temp1 = ndimage.gaussian_filter1d(temp0, sigma_dif)
            else:
                temp1 = gaussian_filter1d(temp0, sigma_dif)             # convolution with variable sigma
            tempNew = temp1 / np.mean(temp1)
            temp[:, i] = tempNew
            

        self.templates    = temp
        self.log_lam_temp = log_lam_temp
        self.teff_grid = par['Teff']
        self.feh_grid  = par['FeH']
        self.logg_grid = par['logg']
        self.LSF       = Temp_FWHM
        self.velscale  = velscale
        
class SingleStar():

    """
    Class of single stelar spectrum

    Parameters
    ----------
    config : class
        Class of configuration
    template : class
        Class of single stellar population template
    mag : float, optional
        Magnitude in SDSS r-band, by default 15.0
    Teff : float, optional
        Effective tempreture, by default 10000.0K
    FeH : float, optional
        Metallicity of stellar, by default 0.0
    vel : float, optional
        Line of sight velocity, by default 100.0km/s
    Ebv : float, optional
        Dust extinction, by default 0.1
    """

    def __init__(self, config, template, mag = 15.0, teff = 10000.0, feh = 0.0, vel = 100.0, ebv = 0.0):
    
        StarTemp = template.templates
        
        # Select metal bins
        idx_FeH = (np.abs(template.feh_grid - feh) < 0.5)
        tpls = StarTemp[:, idx_FeH]
        
        # Select Teff bins
        Teff_FeH = template.teff_grid[idx_FeH]
        minloc   = np.argmin(abs(teff - Teff_FeH))
        starspec = tpls[:, minloc]
        
        wave = np.exp(template.log_lam_temp)
        
        # Dust Reddening
        if np.isscalar(ebv):
            starspec = reddening(wave, starspec, ebv = ebv)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave * (1 + redshift)
        
        flux = np.interp(config.wave, wave_r, starspec)
        
        # Calibration
        if np.isscalar(mag):
            flux = calibrate(config.wave, flux, mag, filtername='SLOAN_SDSS.r')
        
        # Convert to input wavelength
        self.wave = config.wave
        self.flux = flux