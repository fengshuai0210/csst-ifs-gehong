import numpy as np

class instrument():
    
    def __init__(self, wave_min = 3500, wave_max = 10000, 
                 dlam = 1, inst_fwhm = 2.5,
                 nx = 100, ny = 100, dpix = 0.1):

        self.dlam = dlam
        self.wave = np.arange(wave_min, wave_max, dlam)
        self.wave_min = wave_min
        
        self.inst_fwhm = inst_fwhm
        self.nx = nx
        self.ny = ny
        self.dpix = dpix
        self.fov_x = nx * dpix
        self.fov_y = ny * dpix