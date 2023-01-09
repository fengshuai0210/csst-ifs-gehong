import numpy as np

#

def calibrate(wave, flux, mag, filtername='/Users/sfeng/Fastro/Projects/CSST/MCI/SLOAN_SDSS.r'):
    """
    Calibrate the spectra according to the magnitude.
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
# Dust attenuation

def Calzetti_Law(wave, Rv = 4.05):
    """
    Dust Extinction Curve by Calzetti et al. (2000)
    """
    wave_number = 1./(wave * 1e-4)
    reddening_curve = np.zeros(len(wave))
    
    idx = np.logical_and(wave > 1200, wave <= 6300)
    reddening_curve[idx] = 2.659 * ( -2.156 + 1.509 * wave_number[idx] - 0.198 * \
                    (wave_number[idx] ** 2)) + 0.011 * (wave_number[idx] **3 ) + Rv
                                    
    idx = np.logical_and(wave > 6300, wave <= 22000)
    reddening_curve[idx] = 2.659 * ( -1.857 + 1.040 * wave_number[idx]) + Rv
    return reddening_curve

def reddening(wave, flux, ebv = 0, law = 'calzetti', Rv = 4.05):
    """
    Redden an input spectra through a given reddening curve.
    """
    if law == 'calzetti':
        curve = Calzetti_Law(wave, Rv = Rv)
        fluxNew = flux / (10. ** (0.4 * ebv * curve))
    return fluxNew

# -------------
# Emission Line

def SingleEmissinoLine(wave, line_wave, FWHM_inst):
    width = line_wave * FWHM_inst / 1e5
    log2  = np.log(2)
    flux  = (np.exp(- 4. * log2 * (wave - line_wave) ** 2. / (width * width)) /
            (width * np.sqrt(np.pi / log2) / 2.))
    return flux

class EmissionLineTemplate():
    
    def __init__(self, lam_range = [500, 15000], flux_table = 'SDSS'):
        
        self.lam_range = lam_range
        
        if flux_table == 'cigal':
            
            flux_table_file = './data/cigal.nebular.fits'
            # Loading emission line flux table
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
            self.zgas_grid = grid['Zgas']
            self.logu_grid = grid['logU']
            self.ne_grid   = grid['ne']
            
        if flux_table == 'SDSS':
            
            flux_table_file = './data/SDSS_Emline_Table.fits'
            # Loading emission line flux table
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
            self.logOH_grid = grid['logOH']
        
        # Flux ratio
        flux_ratio = line_table[3].data
        self.flux_ratio = flux_ratio[:, w]

class EmissionLine():
    
    def __init__(self, wave, temp, Halpha = 1e-12, 
                 Zgas = False, logU = False, Ne = False, logOH = False,
                 vel = 100, vdisp = 120, Ebv = 0.1):     
        
        if (Zgas > 0) & (Zgas < 0.1) & (logU > -4.1) & (logU < -0.9) & (Ne > 5) & (Ne < 2000):
            indz = np.argmin(np.abs(Zgas - temp.zgas_grid))
            indu = np.argmin(np.abs(logU - temp.logu_grid))
            inde = np.argmin(np.abs(np.log10(Ne) - np.log10(temp.ne_grid)))
            
            flux_ratio = temp.flux_ratio[indz&indu&inde, :]
        else:
            if (logOH > 8) & (logOH < 9.2):
                indz = np.argmin(np.abs(logOH - temp.logOH_grid))
                
                flux_ratio = temp.flux_ratio[indz, :]
            else:
                raise ValueError('something went wrong')
                
        rest_wave = np.arange(temp.lam_range[0], temp.lam_range[1], 0.1)
        
        # Make emission line
        nline = len(temp.line_wave)
        for i in range(nline):
            if i==0:
                emission_line  = emline(rest_wave, temp.line_wave[i], vdisp)
                emission_lines = emission_line
            else:
                emission_line  = emline(rest_wave, temp.line_wave[i], vdisp)
                emission_lines = np.vstack((emission_lines, emission_line))
        emission_lines = emission_lines.T
        
        # Make emission line spectra                  
        emlines = emission_lines * flux_ratio
        flux_combine = np.sum(emlines, axis = 1)
        flux_calibrate = flux_combine * 1e17 * Halpha      # Units: erg/s/A/cm^2
        
        # Dust attenuation
        if np.isscalar(Ebv):
            flux_dust = reddening(rest_wave, flux_calibrate, ebv = Ebv)
            
        # Redshift
        redshift = vel / 3e5
        wave_r = rest_wave * (1 + redshift)
        flux_red = np.interp(wave, wave_r, flux_dust)
            
        self.wave = wave
        self.flux = flux_red
        #self.haflux = haflux

