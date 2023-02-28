from scipy.stats import norm
import numpy as np
from astropy.io import fits

def SingleEmissinoLine(wave, line_wave, FWHM_inst):
    sigma = FWHM_inst / 2.355 
    flux = norm.pdf(wave, line_wave, sigma)
    return flux

class EmissionLineTemplate():
    
    def __init__(self, lam_range = [500, 15000], dlam = 0.1, model = 'nlr', FWHM_inst = 0.5):
        
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