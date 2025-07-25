class Config:
    """
    Configuration class for GEHONG spectral modeling and IFS datacube simulation.

    This class defines the core parameters that control the spectral and spatial structure
    of all simulation outputs, including spectra, maps, and data cubes. It determines the
    wavelength grid, spatial sampling, and pixel scale of the simulated scene, and serves as
    the global configuration object across all GEHONG modules.

    Three usage modes are supported:
    - mode='sim': Preset configuration for use with `csst-ifs-sim` (CCD image simulation).
    - mode='etc': Preset configuration for use with `csst-ifs-etc` (exposure time calculator). [default]
    - mode=None : Fully custom configuration. All parameters must be explicitly provided.

    Parameters
    ----------
    mode : str or None, optional
        Configuration mode. One of ['sim', 'etc', None]. Default is 'etc'.
    wave_min : float, optional
        Minimum wavelength in Angstroms.
    wave_max : float, optional
        Maximum wavelength in Angstroms.
    dlam : float, optional
        Wavelength sampling step in Angstroms.
    nx : int, optional
        Number of spaxels along the X axis (image width).
    ny : int, optional
        Number of spaxels along the Y axis (image height).
    dpix : float, optional
        Spatial pixel size (spaxel scale) in arcseconds.
    """
    def __init__(self, mode='etc', wave_min=None, wave_max=None, dlam=None,
                 nx=None, ny=None, dpix=None):

        # --- Preset configurations ---
        presets = {
            "sim": {
                "wave_min": 3000.0,
                "wave_max": 10500.0,
                "dlam": 1.0,
                "nx": 100,
                "ny": 100,
                "dpix": 0.1,
            },
            "etc": {
                "wave_min": 3500.0,
                "wave_max": 10000.0,
                "dlam": 2.0,
                "nx": 30,
                "ny": 30,
                "dpix": 0.2,
            }
        }

        if mode is not None:
            if mode not in presets:
                raise ValueError(f"[Error] Unknown config mode '{mode}'. Choose from ['sim', 'etc'].")

            params = presets[mode]
            wave_min = params["wave_min"]
            wave_max = params["wave_max"]
            dlam     = params["dlam"]
            nx       = params["nx"]
            ny       = params["ny"]
            dpix     = params["dpix"]

        else:
            # Require all parameters to be provided manually
            if any(p is None for p in [wave_min, wave_max, dlam, nx, ny, dpix]):
                raise ValueError("[Error] When mode=None, all parameters must be manually specified.")

            # --- Validation ---
            if wave_min >= wave_max:
                raise ValueError("[Error] wave_min must be smaller than wave_max.")
            if dlam <= 0:
                raise ValueError("[Error] dlam must be > 0.")
            if nx <= 0 or ny <= 0:
                raise ValueError("[Error] nx and ny must be > 0.")
            if dpix <= 0:
                raise ValueError("[Error] dpix must be > 0.")

        # --- Set attributes ---
        self.wave_min = wave_min
        self.wave_max = wave_max
        self.dlam = dlam
        self.wave = np.arange(wave_min, wave_max, dlam)

        self.nx = nx
        self.ny = ny
        self.dpix = dpix
        self.fov_x = nx * dpix
        self.fov_y = ny * dpix