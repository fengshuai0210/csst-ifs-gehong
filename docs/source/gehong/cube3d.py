class Cube3D:
    """
    Cube3D(config, stellar_map=None, gas_map=None)

    Class for constructing 3D spectral data cubes (RA × DEC × Wavelength)
    including stellar continuum and/or ionized gas emission.

    Parameters
    ----------
    config : object
        Configuration object defining spatial and spectral grids.
        Must contain attributes: nx, ny, dpix, wave, dlam, inst_fwhm.
    stellar_map : object, optional
        2D maps of stellar population properties with shape (nx, ny).
    gas_map : object, optional
        2D maps of ionized gas properties with shape (nx, ny).
    objname : str, optional
        Object name. Used for FITS header and default filename.
    ra : float, optional
        Right Ascension of cube center (in degrees). Default is 0.0 if not specified.
    dec : float, optional
        Declination of cube center (in degrees). Default is 0.0 if not specified.

    Notes
    -----
    - Internal cube shape is (nx, ny, nz), where:
        * nx: number of RA pixels (x-axis)
        * ny: number of DEC pixels (y-axis)
        * nz: number of wavelength bins
    - flux[i, j, :] corresponds to the spectrum at pixel (i_RA, j_DEC).
    - Units of output flux: 1e-17 erg/s/cm^2/Å.
    """

    def __init__(self, config, stellar_map=None, gas_map=None, objname=None,
                ra=None, dec=None):
        
        self.config = config

        # Spatial and spectral dimensions
        self.nx = config.nx  # RA pixels
        self.ny = config.ny  # DEC pixels
        self.dpix = config.dpix  # pixel size [arcsec]
        self.fov_x = config.fov_x  # [arcsec]
        self.fov_y = config.fov_y  # [arcsec]
        self.wave = config.wave  # wavelength grid [Å]
        self.nz = len(self.wave)
        self.wave0 = np.min(self.wave)  # starting wavelength
        #self.inst_fwhm = config.inst_fwhm  # instrument FWHM [Å]

        # Main datacube: shape = (RA, DEC, WAVE)
        self.flux = np.zeros((self.nx, self.ny, self.nz))

        # Input 2D parameter maps
        self.stellar_map = stellar_map
        self.gas_map = gas_map

        # Sky coordinates of cube center [deg]
        self.ra = 40.1
        self.dec = 30.0

        # Set object name
        self.objname = objname if objname is not None else 'sciobj'

    def make_cube(self, stellar_tem=None, hii_tem=None):
        """
        Generate the full spectral cube by summing stellar and ionized gas components.

        At least one of `stellar_map` or `gas_map` must be provided. If both are provided,
        their spectra are summed. If only one is present, only that component is used.

        Parameters
        ----------
        stellar_tem : object
            Template used by StellarContinuum.
        hii_tem : object
            Template used by HII_Region.

        Raises
        ------
        ValueError
            If both stellar_map and gas_map are None.
        """

        if self.stellar_map is None and self.gas_map is None:
            raise ValueError("[Cube3D.make_cube] Both stellar_map and gas_map are None. Cannot generate cube.")

        def _val(attr, i, j):
            """Return scalar, None, or indexed 2D array value."""
            if attr is None or np.isscalar(attr):
                return attr
            return attr[i, j]

        def _array2d(arr1d, age):
            """Convert 1D array + age to 2D array for StellarContinuum."""
            return np.column_stack([age, arr1d])

        def _sfh_ceh_input(arr, age, i, j):
            """Process SFH or CEH array into proper 2D format if needed."""
            if arr is None or np.isscalar(arr):
                return arr
            if arr.ndim == 2:
                return arr[i, j]
            if arr.ndim == 3:
                return _array2d(arr[i, j, :], age)
            raise ValueError(f"Unsupported shape {arr.shape} for SFH/CEH")

        smap = self.stellar_map
        gmap = self.gas_map

        for i in range(self.nx):
            for j in range(self.ny):
                flux_star = None
                flux_gas = None

                # === Stellar component ===
                if smap is not None:
                    mag = _val(smap.mag, i, j)
                    use_mag = mag is not None

                    sfh = _sfh_ceh_input(smap.sfh, smap.sfh_age, i, j)
                    ceh = _sfh_ceh_input(smap.ceh, smap.ceh_age, i, j)
                    vdisp = _val(smap.vdisp, i, j)
                    ebv = _val(smap.ebv, i, j)

                    if use_mag:
                        vel = _val(smap.vel, i, j)
                        z = None
                        vpec = None
                    else:
                        vel = 0.0  # ignored
                        z = _val(getattr(smap, 'z', None), i, j)
                        vpec = _val(getattr(smap, 'vpec', None), i, j)

                    ss = StellarContinuum(self.config, stellar_tem,
                                        mag=mag, sfh=sfh, ceh=ceh,
                                        vel=vel, vdisp=vdisp, ebv=ebv,
                                        z=z, vpec=vpec)
                    flux_star = ss.flux

                # === Gas component ===
                if gmap is not None:
                    halpha = _val(gmap.halpha, i, j)
                    use_halpha = halpha is not None

                    logz = _val(gmap.zh, i, j)
                    vdisp = _val(gmap.vdisp, i, j)
                    ebv = _val(gmap.ebv, i, j)

                    if use_halpha:
                        sfr = None
                        z = None
                        vpec = None
                        vel = _val(gmap.vel, i, j)
                    else:
                        sfr = _val(getattr(gmap, 'sfr', None), i, j)
                        z = _val(getattr(gmap, 'z', None), i, j)
                        vpec = _val(getattr(gmap, 'vpec', None), i, j)
                        vel = 0.0  # ignored in SFR mode

                    gg = HII_Region(self.config, hii_tem,
                                    halpha=halpha, sfr=sfr,
                                    z=z, vpec=vpec, logz=logz,
                                    vel=vel, vdisp=vdisp, ebv=ebv)
                    flux_gas = gg.flux

                # === Combine ===
                if flux_star is not None and flux_gas is not None:
                    self.flux[i, j, :] = flux_star + flux_gas
                elif flux_star is not None:
                    self.flux[i, j, :] = flux_star
                elif flux_gas is not None:
                    self.flux[i, j, :] = flux_gas


    def wcs_info(self):
        """
        Construct a WCS header consistent with FITS axis order.

        FITS convention:
        - Axis 1 → WAVE (fastest axis in memory)
        - Axis 2 → DEC (Y)
        - Axis 3 → RA  (X)

        Stores result in self.wcs_header.
        """
        wcs_dict = {
            'CTYPE1': 'WAVE', 'CUNIT1': 'Angstrom', 'CDELT1': self.config.dlam,
            'CRPIX1': 1, 'CRVAL1': np.min(self.wave),

            'CTYPE2': 'DEC--TAN', 'CUNIT2': 'deg',
            'CDELT2': self.dpix / 3600.,  # arcsec → deg
            'CRPIX2': np.round(self.ny / 2.), 'CRVAL2': self.dec,

            'CTYPE3': 'RA---TAN', 'CUNIT3': 'deg',
            'CDELT3': self.dpix / 3600.,
            'CRPIX3': np.round(self.nx / 2.), 'CRVAL3': self.ra,

            'BUNIT': '10**(-17)*erg/s/cm**2/Angstrom'
        }
        input_wcs = astropy.wcs.WCS(wcs_dict)
        self.wcs_header = input_wcs.to_header()


    def insert_spec(self, spec, dx_arcsec=0.0, dy_arcsec=0.0):
        """
        Add a 1D spectrum into the cube center with optional spatial offset (in arcsec).

        Parameters
        ----------
        spec : object
            Must contain `.flux` of shape (nz,)
        dx_arcsec : float
            Offset in RA direction [arcsec] (positive = east)
        dy_arcsec : float
            Offset in DEC direction [arcsec] (positive = north)

        Notes
        -----
        The spectrum is inserted at the pixel nearest to (center + offset).
        """
        # Convert arcsec to pixel offsets
        dx_pix = int(round(dx_arcsec / self.dpix))
        dy_pix = int(round(dy_arcsec / self.dpix))

        x0 = int(np.round(self.nx / 2.))
        y0 = int(np.round(self.ny / 2.))
        i = x0 + dx_pix
        j = y0 + dy_pix

        # Check bounds
        if 0 <= i < self.nx and 0 <= j < self.ny:
            self.flux[i, j, :] += spec.flux
        else:
            raise IndexError(f"[insert_spec] Offset ({dx_arcsec}\", {dy_arcsec}\") out of bounds.")


    def savefits(self, filename=None, path='./'):
        """
        Save the spectral cube into a FITS file with WCS and metadata.

        Parameters
        ----------
        filename : str, optional
            Output FITS filename (with or without .fits suffix). If None,
            defaults to "<object_name>.fits"
        path : str, optional
            Output directory. Default is current directory.
        """
        from gehong import __version__

        # Set filename from object name if not provided
        if filename is None:
            filename = f"{self.objname}.fits"
        elif not filename.endswith('.fits'):
            filename = filename + '.fits'

        # Primary header
        hdr = fits.Header()
        hdr['FILETYPE'] = 'SCICUBE'
        hdr['CODE'] = 'CSST-IFS-GEHONG'
        hdr['VERSION'] = __version__
        hdr['OBJECT'] = self.objname
        hdr['RA'] = self.ra
        hdr['DEC'] = self.dec

        hdu0 = fits.PrimaryHDU(header=hdr)

        self.wcs_info()
        hdu1 = fits.ImageHDU(self.flux, header=self.wcs_header)

        hdulist = fits.HDUList([hdu0, hdu1])
        hdulist.writeto(os.path.join(path, filename), overwrite=True)