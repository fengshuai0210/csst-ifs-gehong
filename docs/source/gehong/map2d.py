def Sersic2D(x, y, mag=12.0, r_eff=1.0, n=2.0, ellip=0.5,
             theta=0.0, x_0=0.0, y_0=0.0, pixelscale=0.01):
    """
    Compute the 2D surface brightness distribution based on the Sérsic profile.

    Parameters
    ----------
    x : array_like
        X-coordinates of spatial positions (in pixels).
    y : array_like
        Y-coordinates of spatial positions (in pixels).
    mag : float, optional
        Total integrated magnitude of the Sérsic profile. Default is 12.0 mag.
    r_eff : float, optional
        Effective (half-light) radius in pixels. Default is 1.0.
    n : float, optional
        Sérsic index. Default is 2.0.
    ellip : float, optional
        Ellipticity, defined as 1 - b/a. Default is 0.5.
    theta : float, optional
        Position angle (PA) in degrees, measured counterclockwise from +X axis. Default is 0.0.
    x_0 : float, optional
        X-coordinate of the Sérsic model center. Default is 0.0.
    y_0 : float, optional
        Y-coordinate of the Sérsic model center. Default is 0.0.
    pixelscale : float, optional
        Pixel area in arcsec² (used for magnitude normalization). Default is 0.01.

    Returns
    -------
    sb_mag : ndarray
        Surface brightness distribution in units of AB magnitude per arcsec².
    """

    # Convert the angle in degree to angle to radians
    theta_radian = theta / 180 * np.pi

    # Produce Sersic profile
    bn = sp.gammaincinv(2. * n, 0.5)
    a, b = r_eff, (1 - ellip) * r_eff
    cos_theta, sin_theta = np.cos(theta_radian), np.sin(theta_radian)
    x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
    x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
    z = (abs(x_maj / a) ** 2 + abs(x_min / b) ** 2) ** (1 / 2)
    profile = np.exp(-bn * (z ** (1 / n) - 1))

    # Normalization
    integral = a * b * 2 * np.pi * n * np.exp(bn) / (bn ** (2 * n)) * sp.gamma(2 * n)
    prof_norm = profile / integral * pixelscale

    # Calibration
    total_flux = 10. ** ((22.5 - mag) * 0.4)
    sb_mag = 22.5 - 2.5 * np.log10(prof_norm * total_flux / pixelscale)

    return sb_mag


def VelMap2D(x, y, vmax=200.0, rt=1.0, ellip=0.5,
             theta=0.0, x_0=0.0, y_0=0.0):
    """
    Compute a 2D velocity field based on a rotating disk model.

    The rotation curve follows a hyperbolic tangent (tanh) profile:
        v(r) = vmax * tanh(r/rt)

    Parameters
    ----------
    x : array_like
        X-coordinates of spatial positions (in pixels).
    y : array_like
        Y-coordinates of spatial positions (in pixels).
    vmax : float, optional
        Asymptotic (maximum) rotational velocity in km/s. Default is 200.0.
    rt : float, optional
        Turn-over radius of the rotation curve in pixels. Default is 1.0.
    ellip : float, optional
        Ellipticity of the disk. Default is 0.5.
    theta : float, optional
        Position angle (PA) in degrees, counterclockwise from +X axis. Default is 0.0.
    x_0 : float, optional
        X-coordinate of disk center. Default is 0.0.
    y_0 : float, optional
        Y-coordinate of disk center. Default is 0.0.

    Returns
    -------
    velocity_map : ndarray
        Line-of-sight velocity field in km/s.
    """

    # Convert the angle in degree to angle to radians
    theta_radian = theta / 180 * np.pi

    # Produce tanh profile
    a, b = rt, (1 - ellip) * rt
    cos_theta, sin_theta = np.cos(theta_radian), np.sin(theta_radian)
    x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
    x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
    z = (abs(x_maj / a) ** 2 + abs(x_min / b) ** 2) ** (1 / 2)
    profile = vmax * np.tanh(z) * ((x_maj / a) / z)

    return profile


def GradMap2D(x, y, a0=10.0, r_eff=1.0, gred=-0.3, ellip=0.5,
              theta=0.0, x_0=0.0, y_0=0.0, floor_radius=0.01):
    """
    Compute a 2D map with a logarithmic radial gradient.

    The radial profile follows:
        A(r) = a0 + gred * log10(r / r_eff)

    Parameters
    ----------
    x : array_like
        X-coordinates of spatial positions (in pixels).
    y : array_like
        Y-coordinates of spatial positions (in pixels).
    a0 : float, optional
        Central value at r = r_eff. Default is 10.0 (arbitrary units).
    r_eff : float, optional
        Effective radius in pixels (normalization radius). Default is 1.0.
    gred : float, optional
        Gradient per dex in radius (dA / dlog(r)). Default is -0.3.
    ellip : float, optional
        Ellipticity of the disk. Default is 0.5.
    theta : float, optional
        Position angle in degrees (counterclockwise from +X axis). Default is 0.0.
    x_0 : float, optional
        X-coordinate of the center. Default is 0.0.
    y_0 : float, optional
        Y-coordinate of the center. Default is 0.0.
    floor_radius : float, optional
        Minimum radius for log calculation (to avoid log(0)). Default is 0.01.

    Returns
    -------
    map : ndarray
        2D distribution with log-radius gradient (arbitrary units).
    """

    # Convert PA from degrees to radians
    theta_rad = np.deg2rad(theta)

    # Axis lengths
    a, b = r_eff, (1 - ellip) * r_eff

    # Rotate coordinates to major/minor axes
    dx = x - x_0
    dy = y - y_0
    x_maj = dx * np.cos(theta_rad) + dy * np.sin(theta_rad)
    x_min = -dx * np.sin(theta_rad) + dy * np.cos(theta_rad)

    # Compute elliptical radius
    r = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)
    r = np.clip(r, floor_radius, None)  # Avoid log(0)

    # Compute profile: log-radius gradient
    profile = a0 + gred * np.log10(r)

    return profile


class Map2D:
    """
    A class to generate 2D maps of galaxy properties based on analytical models or user-defined inputs.

    This is a core utility for constructing spatial distributions (maps) of galaxy properties such as 
    surface brightness (Sérsic profile), velocity fields (tanh disk rotation), and parameter gradients 
    (e.g., age, metallicity).
    """

    def __init__(self, config):
        """
        Initialize the coordinate grids based on the provided configuration.

        Parameters
        ----------
        config : Config
            Configuration object defining the spatial and spectral layout.
        """
        self.xsamp = config.dpix
        self.ysamp = config.dpix

        # Define 1D x and y coordinates in arcseconds
        xvals = np.linspace(-(config.nx - 1) / 2.0 * self.xsamp,
                            (config.nx - 1) / 2.0 * self.xsamp, config.nx)
        yvals = np.linspace(-(config.ny - 1) / 2.0 * self.ysamp,
                            (config.ny - 1) / 2.0 * self.ysamp, config.ny)

        # Create 2D coordinate grids (shape: [ny, nx])
        x, y = np.meshgrid(xvals, yvals)

        self.nx = config.nx
        self.ny = config.ny
        self.x = x         # x-coordinate at each pixel [arcsec]
        self.y = y         # y-coordinate at each pixel [arcsec]
        self.row = xvals   # x-axis values (columns)
        self.col = yvals[::-1]  # flipped y-axis values for consistent image orientation

    def shift_rotate(self, yoff, xoff, rot):
        pa_radians = np.deg2rad(rot)
        xsh = self.x - xoff
        ysh = self.y - yoff
        xsh_rot = xsh * np.cos(pa_radians) + ysh * np.sin(pa_radians)
        ysh_rot = -xsh * np.sin(pa_radians) + ysh * np.cos(pa_radians)
        return ysh_rot, xsh_rot

    def _check_theta(self, theta):
        if abs(theta) > 180:
            theta_wrapped = ((theta + 180) % 360) - 180
            print(f"[Warning] Position angle (theta={theta}) is outside [-180, 180] range. Auto-wrapped to {theta_wrapped:.1f}°.")
            return theta_wrapped
        return theta

    def sersic_map(self, mag=12.0, r_eff=2.0, n=2.5, ellip=0.5, theta=-50.0):
        if mag < 8 or mag > 26:
            print(f"[Warning] Input mag={mag:.2f} is outside the 8–26 recommended range.")
        if r_eff <= 0:
            raise ValueError("[Error] Effective radius (r_eff) must be > 0 arcsec.")
        if n <= 0:
            raise ValueError("[Error] Sérsic index (n) must be > 0.")
        if not (0 <= ellip < 1):
            raise ValueError("[Error] Ellipticity must be in [0, 1).")

        theta = self._check_theta(theta)

        self.map = Sersic2D(self.x, self.y, mag=mag,
                            r_eff=r_eff / self.xsamp, n=n,
                            ellip=ellip, theta=theta,
                            pixelscale=self.xsamp * self.ysamp)

    def tanh_map(self, vmax=200.0, rt=2.0, ellip=0.5, theta=-50.0):
        if vmax <= 0:
            raise ValueError("[Error] Maximum rotational velocity (vmax) must be > 0 km/s.")
        if rt <= 0:
            raise ValueError("[Error] Turn-over radius (rt) must be > 0 arcsec.")
        if not (0 <= ellip < 1):
            raise ValueError("[Error] Ellipticity must be in [0, 1).")

        theta = self._check_theta(theta)

        self.map = VelMap2D(self.x, self.y, vmax=vmax, rt=rt / self.xsamp,
                            ellip=ellip, theta=theta)

    def gred_map(self, a0=10.0, r_eff=1.0, gred=-1.0, ellip=0.5, theta=0.0):
        if r_eff <= 0:
            raise ValueError("[Error] Effective radius (r_eff) must be > 0 arcsec.")
        if not (0 <= ellip < 1):
            raise ValueError(f"[Error] Ellipticity ellip={ellip} must be in [0, 1).")

        theta = self._check_theta(theta)

        self.map = GradMap2D(self.x, self.y, a0=a0, r_eff=r_eff / self.xsamp,
                             gred=gred, ellip=ellip, theta=theta)

    def load_map(self, image):
        if np.ndim(image) != 2:
            raise ValueError("[Error] Input image must be a 2D array.")
        self.map = resize(image, (self.nx, self.ny))


class StellarPopulationMap:
    """
    Container for 2D and 3D maps of stellar population properties used in spectrum synthesis.

    Supports both single-SSP and composite SFH/CEH modes, and handles observational vs physical flux calibration.
    Each spaxel's parameters are passed to the `StellarContinuum` class to generate the corresponding spectrum.

    Parameters
    ----------
    config : Config
        Global configuration object.
    sb : Map2D, optional
        Surface brightness map (mag/arcsec²), used for flux calibration in mag-mode.
    sfh : Map2D or ndarray, optional
        Either 2D map of mean stellar age (Gyr), or 3D array of shape (nx, ny, N) representing SFR(t) history.
        The corresponding time axis must be provided via `sfh_age` of shape (N,).
    sfh_age : ndarray, optional
        Array of shape (N,) representing the age grid (Gyr) used in SFH cube.
    ceh : Map2D or ndarray, optional
        Either 2D map of [Fe/H], or 3D array of shape (nx, ny, N) representing [Fe/H](t) history.
        The corresponding time axis must be provided via `ceh_age` of shape (N,).
    ceh_age : ndarray, optional
        Array of shape (N,) representing the age grid (Gyr) used in CEH cube.
    vel : Map2D, optional
        2D map of line-of-sight velocity (km/s). Required in mag-mode.
    vdisp : Map2D, optional
        2D map of stellar velocity dispersion (km/s).
    ebv : Map2D, optional
        2D map of dust extinction (E(B−V)).
    vpec : Map2D, optional
        2D map of galaxy peculiar velocity (km/s). Required in physical mode.
    z : float, optional
        Global cosmological redshift. Required in physical mode.
    """

    def __init__(self, config, sb=None, sfh=None, sfh_age=None, ceh=None, ceh_age=None,
                vel=None, vdisp=None, ebv=None, vpec=None, z=None):

        self.nx = config.nx
        self.ny = config.ny
        self.dpix = config.dpix
        self.fov_x = config.fov_x
        self.fov_y = config.fov_y
        self.z = z

        # Surface brightness and magnitude
        self.sb = sb.map if sb else None
        self.mag = self.sb - 2.5 * np.log10(self.dpix ** 2) if self.sb is not None else None

        # SFH map or cube
        self.sfh = None
        if sfh is not None:
            if hasattr(sfh, 'map'):
                self.sfh = sfh.map
            else:
                self.sfh = self._parse_3d_map(sfh, 'sfh', self.nx, self.ny, sfh_age)
        if self.sfh is None:
            raise ValueError("[Error] Missing stellar age map or SFH cube.")

        self.sfh_age = sfh_age
        if self.sfh_age is not None and np.any(self.sfh_age < 0):
            raise ValueError("[Error] sfh_age contains negative values.")

        # CEH map or cube
        self.ceh = None
        if ceh is not None:
            if hasattr(ceh, 'map'):
                self.ceh = ceh.map
            else:
                self.ceh = self._parse_3d_map(ceh, 'ceh', self.nx, self.ny, ceh_age)
        if self.ceh is None:
            raise ValueError("[Error] Missing [Fe/H] map or CEH cube.")

        self.ceh_age = ceh_age
        if self.ceh_age is not None and np.any(self.ceh_age < 0):
            raise ValueError("[Error] ceh_age contains negative values.")

        self.vel = vel.map if vel else None
        self.vpec = vpec.map if vpec else None

        # --- Mode check: magnitude-based or physical ---
        if self.mag is not None:  # mag-based mode
            if self.vel is None:
                raise ValueError("[Error] mag-based mode requires stellar velocity map (vel).")
            if self.sfh.ndim == 3 and self.sfh_age is None:
                raise ValueError("[Error] SFH cube provided, but missing sfh_age in mag-based mode.")
            if self.ceh.ndim == 3 and self.ceh_age is None:
                raise ValueError("[Error] CEH cube provided, but missing ceh_age in mag-based mode.")
        else:  # physically-calibrated mode
            if self.sfh.ndim != 3:
                raise ValueError("[Error] Physical mode requires 3D SFH cube.")
            if self.sfh_age is None:
                raise ValueError("[Error] Physical mode requires sfh_age for SFH cube.")
            if self.vpec is None or self.z is None:
                raise ValueError("[Error] Physical mode requires vpec and z.")
            if self.ceh.ndim == 3 and self.ceh_age is None:
                raise ValueError("[Error] CEH cube provided, but missing ceh_age in physical mode.")

        # vdisp check
        self.vdisp = None
        if vdisp:
            self.vdisp = vdisp.map
            bad = self.vdisp < 10
            if np.any(bad):
                print(f"[Warning] {np.sum(bad)} spaxels with vdisp < 10 km/s reset to 10.")
                self.vdisp[bad] = 10.0
        else:
            raise ValueError("[Error] Stellar velocity dispersion (vdisp) map is required.")

        # ebv check
        self.ebv = None
        if ebv:
            self.ebv = ebv.map
            bad = self.ebv < 0
            if np.any(bad):
                print(f"[Warning] {np.sum(bad)} spaxels with negative EBV reset to 0.")
                self.ebv[bad] = 0.0
        else:
            raise ValueError("[Error] Dust extinction (ebv) map is required.")

    def _parse_3d_map(self, input_array, name, nx, ny, age_axis=None):
        """
        Parse and validate a 3D stellar population property map.

        Parameters
        ----------
        input_array : ndarray
            3D array with shape (nx, ny, N).
        name : str
            Name of the parameter for error messages.
        nx, ny : int
            Expected spatial dimensions.
        age_axis : ndarray, optional
            Associated age/time array.

        Returns
        -------
        ndarray
            Validated array in shape (nx, ny, N).
        """
        if isinstance(input_array, np.ndarray):
            if input_array.ndim == 2:
                return input_array
            elif input_array.ndim == 3:
                if input_array.shape[:2] != (nx, ny):
                    raise ValueError(f"[Error] {name} shape must be (nx, ny, N), got {input_array.shape}.")
                if age_axis is not None and input_array.shape[2] != len(age_axis):
                    raise ValueError(f"[Error] Length of {name}_age does not match third dimension of {name} cube.")
                return input_array
            else:
                raise ValueError(f"[Error] {name} must be 2D or 3D numpy array.")
        else:
            raise TypeError(f"[Error] {name} must be a numpy array or Map2D.")


class IonizedGasMap:
    """
    Container for 2D maps of ionized gas properties used in nebular emission modeling.

    Supports both Hα-calibrated and SFR-calibrated modes, determined by input availability.
    Each spaxel's parameters are passed to the `HII_Region` class to generate the emission line spectrum.

    Parameters
    ----------
    config : Config
        Global configuration object.
    halpha : Map2D, optional
        Hα flux map in units of 1e-17 erg/s/cm². Required for Hα-calibrated mode.
    sfr : Map2D, optional
        Star formation rate map (M☉/yr). Required for SFR-calibrated mode.
    zh : Map2D, optional
        Gas-phase metallicity map (in log Z/Z☉). Required.
    vel : Map2D, optional
        Line-of-sight velocity map (km/s). Required in Hα mode.
    vpec : Map2D, optional
        Galaxy peculiar velocity map (km/s). Required in SFR mode.
    vdisp : Map2D, optional
        Velocity dispersion map (km/s). Required.
    ebv : Map2D, optional
        Dust extinction map E(B−V). Required.
    z : float, optional
        Global redshift. Required in SFR mode.
    """

    def __init__(self, config, halpha=None, sfr=None, zh=None,
                vel=None, vpec=None, vdisp=None, ebv=None, z=None):

        self.nx = config.nx
        self.ny = config.ny
        self.dpix = config.dpix
        self.fov_x = config.fov_x
        self.fov_y = config.fov_y
        self.z = z

        # === Required: zh ===
        if zh is None:
            raise ValueError("[Error] Gas-phase metallicity (zh) map is required.")
        self.zh = zh.map

        # === Required: vdisp ===
        if vdisp is None:
            raise ValueError("[Error] Gas velocity dispersion (vdisp) map is required.")
        self.vdisp = vdisp.map
        bad_vdisp = self.vdisp < 10
        if np.any(bad_vdisp):
            print(f"[Warning] {np.sum(bad_vdisp)} spaxels with vdisp < 10 km/s reset to 10.")
            self.vdisp[bad_vdisp] = 10.0

        # === Required: ebv ===
        if ebv is None:
            raise ValueError("[Error] Dust extinction (ebv) map is required.")
        self.ebv = ebv.map
        bad_ebv = self.ebv < 0
        if np.any(bad_ebv):
            print(f"[Warning] {np.sum(bad_ebv)} spaxels with negative EBV reset to 0.")
            self.ebv[bad_ebv] = 0.0

        # === Optional inputs ===
        self.halpha = halpha.map if halpha else None
        self.sfr = sfr.map if sfr else None
        self.vel = vel.map if vel else None
        self.vpec = vpec.map if vpec else None

        # === Mode Check: Hα vs SFR ===
        if self.halpha is not None:
            if self.vel is None:
                raise ValueError("[Error] Hα-calibrated mode requires gas velocity map (vel).")
        elif self.sfr is not None:
            if self.vpec is None:
                raise ValueError("[Error] SFR-calibrated mode requires peculiar velocity map (vpec).")
            if self.z is None:
                raise ValueError("[Error] SFR-calibrated mode requires global redshift (z).")
        else:
            raise ValueError("[Error] Either halpha map or sfr map must be provided for flux calibration.")
