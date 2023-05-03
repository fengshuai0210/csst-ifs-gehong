import numpy as np

class config():
    """
    The configuration of spectral modeling. The default value is used for ETC calculation.

    Parameters
    ----------
    wave_min : float, optional
        Minimum value of wavelength coverage, by default 3500.0A
    wave_max : float, optional
        Minimum value of wavelength coverage, by default 10000.0A
    dlam : float, optional
        Wavelength width of each spaxel, by default 2.0A
    inst_fwhm : float, optional
        Spectral resolution, by default 0.1A
    nx : int, optional
        Number of spaxel in a spatial direction, by default 30
    ny : int, optional
        Number of spaxel in a spatial direction, by default 30
    dpix : float, optional
        Pixel size in the spatial direction, by default 0.2arcsec
    """
    def __init__(self, wave_min = 3500.0, wave_max = 10000.0, 
                 dlam = 2.0, inst_fwhm = 0.1,
                 nx = 30, ny = 30, dpix = 0.2):
        self.dlam = dlam
        self.wave = np.arange(wave_min, wave_max, dlam)
        self.wave_min = wave_min
        
        self.inst_fwhm = inst_fwhm
        self.nx = nx
        self.ny = ny
        self.dpix = dpix
        self.fov_x = nx * dpix
        self.fov_y = ny * dpix