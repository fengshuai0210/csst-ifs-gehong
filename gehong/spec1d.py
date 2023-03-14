#import numpy as np

def calibrate(wave, flux, mag, filtername='/Users/sfeng/Fastro/Projects/CSST/MCI/SLOAN_SDSS.r'):
    """
    Calibrate the spectra according to the magnitude.

    :param wave: _description_
    :type wave: _type_
    :param flux: _description_
    :type flux: _type_
    :param mag: _description_
    :type mag: _type_
    :param filtername: _description_, defaults to '/Users/sfeng/Fastro/Projects/CSST/MCI/SLOAN_SDSS.r'
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
        filter_file = filtername+'.filter'
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
    
    def __init__(self, lam_range = [500, 15000], dlam = 0.1, model = 'nlr', FWHM_inst = 0.5):

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
        self.FWHM_inst = FWHM_inst
        self.model = model
        
        # HII region model of fsps-cloudy
        if model == 'fsps':
            
            # Loading emission line flux table
            flux_table_file = './data/fsps.nebular.fits'
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
            flux_table_file = './data/AGN.NLR.fits'
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
                emission_line  = SingleEmissinoLine(self.wave, line_wave[i], FWHM_inst)
                emission_lines = emission_line
            else:
                emission_line  = SingleEmissinoLine(self.wave, line_wave[i], FWHM_inst)
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
    
    def __init__(self, wave, temp, Halpha = 100, logZ = False,
                 vel = 100, vdisp = 120, Ebv = 0.1):     
        
        if (logZ > -2) & (logZ < 0.5):
            indz = np.argmin(np.abs(logZ - temp.logz_grid))
            flux_ratio = temp.flux_ratio[indz, :]
        else:
            raise ValueError('The range of logZ is not correct!')
        
        # Make emission line spectra through adding emission lines                 
        emlines = temp.emission_lines * flux_ratio
        flux_combine = np.sum(emlines, axis = 1)
        flux_calibrate = flux_combine * Halpha      # Units: erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(temp.wave, flux_calibrate, ebv = Ebv)
            
        # Broadening caused by Velocity Dispersion
        velscale = 10
        lam_range = [np.min(temp.wave), np.max(temp.wave)]
        flux_logwave, logLam = util.log_rebin(lam_range, flux_dust, velscale=velscale)[:2]
        
        sigma_gas = vdisp / velscale                                # in pixel
        sigma_LSF = temp.FWHM_inst / (np.exp(logLam)) * 3e5 / velscale          # in pixel
        
        if sigma_gas>0: 
            sigma_dif = np.zeros(len(flux_logwave))
            idx = (sigma_gas > sigma_LSF)
            sigma_dif[idx] = np.sqrt(sigma_gas ** 2. - sigma_LSF[idx] ** 2.)
            idx = (sigma_gas <= sigma_LSF)
            sigma_dif[idx] = 0.1
            flux_broad = util.gaussian_filter1d(flux_logwave, sigma_dif)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = np.exp(logLam) * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_broad)
            
        self.wave = wave
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
    
    def __init__(self, wave, temp, Halpha = 100, logZ = False,
                 vel = 100, vdisp = 120, Ebv = 0.1):
        
        if (logZ > -2) & (logZ < 0.5):
            indz = np.argmin(np.abs(logZ - temp.logz_grid))
            flux_ratio = temp.flux_ratio[indz, :]
        else:
            raise ValueError('The range of logZ is not correct!')
        
        # Make emission line spectra through adding emission lines                 
        emlines = temp.emission_lines * flux_ratio
        flux_combine = np.sum(emlines, axis = 1)
        flux_calibrate = flux_combine * Halpha      # Units: 1e-17 erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(temp.wave, flux_calibrate, ebv = Ebv)
            
        # Broadening caused by Velocity Dispersion
        velscale = 10
        lam_range = [np.min(temp.wave), np.max(temp.wave)]
        flux_logwave, logLam = util.log_rebin(lam_range, flux_dust, velscale=velscale)[:2]
        
        sigma_gas = vdisp / velscale                                # in pixel
        sigma_LSF = temp.FWHM_inst / (np.exp(logLam)) * 3e5 / velscale          # in pixel
        
        if sigma_gas>0: 
            sigma_dif = np.zeros(len(flux_logwave))
            idx = (sigma_gas > sigma_LSF)
            sigma_dif[idx] = np.sqrt(sigma_gas ** 2. - sigma_LSF[idx] ** 2.)
            idx = (sigma_gas <= sigma_LSF)
            sigma_dif[idx] = 0.1
            flux_broad = util.gaussian_filter1d(flux_logwave, sigma_dif)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = np.exp(logLam) * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_broad)
            
        self.wave = wave
        self.flux = flux_red

class AGN_BLR():

    """
    Class for the broad line region of AGN
    """
    
    def __init__(self, wave, Hbeta_Flux = 100, Hbeta_FWHM = 2000, Ebv = None, vel = 0.,
                 lam_range = [500, 15000]):
        
        wave_rest = np.arange(lam_range[0], lam_range[1], 0.1)
        
        line_names = ['Hepsilon', 'Hdelta', 'Hgamma', 'Hbeta', 'Halpha']
        line_waves = [3970.079, 4101.742, 4340.471, 4861.333, 6562.819]
        line_ratio = [0.101, 0.208, 0.405, 1.000, 2.579]              # From Ilic et al. (2006)
        
        # Make emission lines
        for i in range(len(line_names)):
            if i==0:
                emission_line  = SingleEmissinoLine(wave_rest, line_waves[i], 
                                                    Hbeta_FWHM / 3e5 * line_waves[i])
                emission_lines = emission_line
            else:
                emission_line  = SingleEmissinoLine(wave_rest, line_waves[i], 
                                                    Hbeta_FWHM / 3e5 * line_waves[i])
                emission_lines = np.vstack((emission_lines, emission_line))
        emlines = emission_lines.T * line_ratio
        flux_combine = np.sum(emlines, axis = 1)
        
        # Flux callibration
        flux_calibrate = flux_combine * Hbeta_Flux      # Units: 1e-17 erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = Ebv)
        else:
            flux_dust = flux_calibrate
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave_rest * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_dust)
            
        self.wave = wave
        self.flux = flux_red

class AGN_FeII():

    """
    Class for FeII emission lines of AGN
    """
    
    def __init__(self, wave, Hbeta_Broad = 100, R4570 = 0.4, Ebv = None, vel = 0.,
                 filename = './data/FeII.AGN.fits'):
        
        """
         _summary_
        """
        
        # Loading FeII template
        hdulist = fits.open(filename)
        data = hdulist[1].data
        wave_rest  = data['WAVE']
        flux_model = data['FLUX']
        
        # Determine the flux of FeII
        Fe4570_temp = 100
        Fe4570_model = Hbeta_Broad * R4570
        Ratio_Fe4570 = Fe4570_model / Fe4570_temp
        
        # Flux calibration
        flux_calibrate = flux_model * Ratio_Fe4570
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = Ebv)
        else:
            flux_dust = flux_calibrate
             
        # Redshift
        redshift = vel / 3e5
        wave_r = wave_rest * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_dust)
        
        self.wave = wave
        self.flux = flux_red
        
class AGN_Powerlaw():
    """
    AGN_Powerlaw _summary_
    """
    
    def __init__(self, wave, M5100 = 1000, alpha = -1.5, Ebv = None, vel = 0.,):

        """
         _summary_
        """
        
        wave_rest = np.linspace(1000,20000,10000)
        flux = wave_rest ** alpha
        
        # Flux calibration
        flux_calibrate = calibrate(wave_rest, flux, M5100, filtername='5100')
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(wave_rest, flux_calibrate, ebv = Ebv)
        else:
            flux_dust = flux_calibrate
            
        # Redshift
        redshift = vel / 3e5
        wave_r = wave_rest * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_dust)
        
        self.wave = wave
        self.flux = flux_red
        
def BHmass_to_M5100(BHmass, Edd_ratio = 0.05, Dist = 21):
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
    Ledd = 3e4 * BHmass
    Lbol = Ledd * Edd_ratio

    # Convert bolometric luminosity to 5100A luminosity (Marconi et al. 2004)
    L5100 = Lbol / 10.9
    M5100 = 4.86 - 2.5 * np.log10(L5100)
    m5100 = M5100 + 5. * np.log10(Dist * 1e5)

    return m5100