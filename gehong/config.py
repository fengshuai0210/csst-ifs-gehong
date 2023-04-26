import numpy as np

class config():
    """
    
    """
    def __init__(self, wave_min = 3500, wave_max = 10000, 
                 dlam = 1, inst_fwhm = 2.5,
                 nx = 100, ny = 100, dpix = 0.1):
        """
        __init__ _summary_

        Parameters
        ----------
        wave_min : int, optional
            _description_, by default 3500
        wave_max : int, optional
            _description_, by default 10000
        dlam : int, optional
            _description_, by default 1
        inst_fwhm : float, optional
            _description_, by default 2.5
        nx : int, optional
            _description_, by default 100
        ny : int, optional
            _description_, by default 100
        dpix : float, optional
            _description_, by default 0.1
        """
        self.dlam = dlam
        self.wave = np.arange(wave_min, wave_max, dlam)
        self.wave_min = wave_min
        
        self.inst_fwhm = inst_fwhm
        self.nx = nx
        self.ny = ny
        self.dpix = dpix
        self.fov_x = nx * dpix
        self.fov_y = ny * dpix