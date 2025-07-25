data_path = 'GEHONG_DATA_PATH'

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
    if t == np.ndarray or t == np.void:           # array or structured array
        f = map(np.array, zip(*f))

    # In Python 3.x all strings (e.g. name='NGC1023') are Unicode strings by defauls.
    # However genfromtxt() returns byte strings b'NGC1023' for non-numeric columns.
    # To have the same behaviour in Python 3 as in Python 2, I convert the Numpy
    # byte string 'S' type into Unicode strings, which behaves like normal strings.
    # With this change I can read the string a='NGC1023' from a text file and the
    # test a == 'NGC1023' will give True as expected.

    if sys.version >= '3':
        f = [v.astype(str) if v.dtype.char == 'S' else v for v in f]

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
        m = int(n * oversample)
    else:
        m = int(n)

    dLam = np.diff(lamRange)[0] / (n - 1.)       # <- FIXED
    lim = lamRange / dLam + [-0.5, 0.5]          # All in units of dLam
    borders = np.linspace(*lim, num=n + 1)       # Linearly
    logLim = np.log(lim)

    c = 299792.458                               # Speed of light in km/s
    if velscale is None:                         # Velocity scale is set by user
        velscale = (np.diff(logLim)[0]) / m * c  # <- FIXED
    else:
        logScale = velscale / c
        m = int(np.diff(logLim)[0] / logScale)   # <- FIXED
        logLim[1] = logLim[0] + m * logScale

    newBorders = np.exp(np.linspace(*logLim, num=m + 1))  # Logarithmically
    k = (newBorders - lim[0]).clip(0, n - 1).astype(int)

    specNew = np.add.reduceat(spec, k)[:-1]      # Do analytic integral
    specNew *= np.diff(k) > 0                    # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k]) * spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:] * newBorders[:-1]) * dLam)

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
    p = int(np.ceil(np.max(3 * sig)))
    m = 2 * p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2

    n = spec.size
    a = np.zeros((m, n))
    for j in range(m):   # Loop over the small size of the kernel
        a[j, p:-p] = spec[j:n - m + j + 1]

    gau = np.exp(-x2[:, None] / (2 * sig**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel

    conv_spectrum = np.sum(a * gau, 0)

    return conv_spectrum


def calibrate(wave, flux, mag, filtername='SLOAN_SDSS.r'):
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
        wave0 = np.linspace(3000, 10000, 7000)
        response0 = np.zeros(7000)
        response0[(wave0 > 5050) & (wave0 < 5150)] = 1.
    else:
        filter_file = data_path + '/data/filter/' + filtername + '.filter'
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
    realflux = (mag * u.STmag).to(u.erg / u.s / u.cm**2 / u.AA).value

    # Normalization
    flux_ratio = realflux / preflux
    flux_calibrate = flux * flux_ratio * 1e17                        # Units: 10^-17 erg/s/A/cm^2

    return flux_calibrate

# ----------------
# Reddening Module


def Calzetti_Law(wave, Rv=4.05):
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
    wave_number = 1. / (wave * 1e-4)
    reddening_curve = np.zeros(len(wave))

    idx = (wave >= 1200) & (wave < 6300)
    reddening_curve[idx] = 2.659 * (-2.156 + 1.509 * wave_number[idx] - 0.198 * (wave_number[idx] ** 2) + 0.011 * (wave_number[idx] ** 3)) + Rv

    idx = (wave >= 6300) & (wave <= 22000)
    reddening_curve[idx] = 2.659 * (-1.857 + 1.040 * wave_number[idx]) + Rv
    return reddening_curve


def reddening(wave, flux, ebv=0.0, law='calzetti', Rv=4.05):
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
    curve = Calzetti_Law(wave, Rv=Rv)
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

    def __init__(self, config, lam_range=[500, 15000], dlam=0.1, model='hii'):

        self.lam_range = lam_range
        self.wave = np.arange(lam_range[0], lam_range[1], 0.1)
        self.FWHM_inst = 0.1
        self.model = model

        # HII region model of fsps-cloudy
        if model == 'hii':

            # Loading emission line flux table
            flux_table_file = data_path + '/data/fsps.nebular.fits'
            line_table = fits.open(flux_table_file)

            # Emission line list
            line_list = line_table[1].data
            line_wave = line_list['Wave']
            line_names = line_list['Name']

            w = (line_wave > lam_range[0]) & (line_wave < lam_range[1])
            self.line_names = line_names[w]
            self.line_wave = line_wave[w]

            # Make parameter grid
            grid = line_table[2].data
            self.logz_grid = grid['logZ']

        # Narrow line region model of
        if model == 'nlr':

            # Loading emission line flux table
            flux_table_file = data_path + '/data/AGN.NLR.fits'
            line_table = fits.open(flux_table_file)

            # Emission line list
            line_list = line_table[1].data
            line_wave = line_list['Wave']
            line_names = line_list['Name']

            w = (line_wave > lam_range[0]) & (line_wave < lam_range[1])
            self.line_names = line_names[w]
            self.line_wave = line_wave[w]

            # Make parameter grid
            grid = line_table[2].data
            self.logz_grid = grid['logZ']

        # Flux ratio
        flux_ratio = line_table[3].data
        self.flux_ratio = flux_ratio

        # Make emission line
        nline = len(line_wave)
        for i in range(nline):
            if i == 0:
                emission_line = SingleEmissinoLine(self.wave, line_wave[i], self.FWHM_inst)
                emission_lines = emission_line
            else:
                emission_line = SingleEmissinoLine(self.wave, line_wave[i], self.FWHM_inst)
                emission_lines = np.vstack((emission_lines, emission_line))
        self.emission_lines = emission_lines.T


class HII_Region():
    """
    Generate emission line spectra of a star-forming HII region using pre-defined emission line templates.

    Two working modes are supported for flux calibration:

    1. Hα-calibrated mode:
    If `halpha` is provided (default), the total flux of Hα (in units of 1e-17 erg/s/cm²) is used to scale the emission lines.
    In this case, the output spectrum is normalized such that the integrated Hα line matches the given value.

    2. SFR-calibrated mode:
    If `halpha` is set to None and both `sfr`, `z`, and `vpec` are provided, the Hα luminosity is estimated from the SFR
    using the Kennicutt (1998) relation:
        L(Hα) [erg/s] = 7.9 × 10⁴¹ × SFR [M☉/yr],
    which is then converted to observed flux via the luminosity distance (based on `z`), and corrected for cosmic dimming.
    This mode allows physically motivated scaling based on star formation activity.

    Input flexibility:

    - If `halpha` is not None, the model uses Hα-calibration mode and ignores `sfr`, `z`, and `vpec`.
    - If `halpha` is None, the model enters SFR-calibrated mode and requires all of: `sfr`, `z`, and `vpec`.
    - `logz` determines gas-phase metallicity and selects the closest template in the pre-computed grid.
    - The output flux includes optional dust extinction, velocity broadening, and redshift correction.

    Parameters
    ----------
    config : object
        Configuration object that defines the target output wavelength grid.
    temp : object
        Emission line template object precomputed for HII or AGN narrow-line regions.
    halpha : float or None, optional
        Observed flux of the Hα line in units of 1e-17 erg/s/cm² (default: 100.0). If None, enter SFR-calibrated mode.
    sfr : float, optional
        Star formation rate in M☉/yr, used only when `halpha` is None.
    z : float, optional
        Cosmological redshift, used for luminosity distance and cosmic dimming in SFR-calibrated mode.
    vpec : float, optional
        Peculiar velocity in km/s, used together with `z` to compute line-of-sight redshift.
    logz : float, optional
        Gas-phase metallicity (in log Z/Z☉), required to select the emission line ratios from the template (default: 0.0).
        Must be within [-2.0, 0.5], otherwise a ValueError is raised.
    vel : float, optional
        Line-of-sight velocity in km/s for redshifting the emission lines (default: 100.0).
        Only used if `vpec` is not provided.
    vdisp : float, optional
        Velocity dispersion in km/s, used to apply Gaussian broadening to emission lines (default: 120.0).
    ebv : float, optional
        Dust extinction E(B–V), applied using the Calzetti law (default: 0.1).
    """
    def __init__(self, config, temp, halpha=100.0, sfr=None, z=None, vpec=None, logz=0.0, vel=100.0, vdisp=120.0, ebv=0.1):
        self.config = config

        logz_clipped = np.clip(logz, temp.logz_grid.min(), temp.logz_grid.max())
        if logz != logz_clipped:
            print(f"[Warning] Input logz={logz:.2f} is out of range; clipped to {logz_clipped:.2f}")
        indz = np.argmin(np.abs(logz_clipped - temp.logz_grid))
        flux_ratio = temp.flux_ratio[indz, :]

        if vdisp < 0:
            print(f"[Warning] Input vdisp={vdisp:.1f} km/s is invalid. Setting vdisp=0.")
            vdisp = 0.0

        if ebv < 0:
            print(f"[Warning] Input ebv={ebv:.3f} is invalid. Setting ebv=0.")
            ebv = 0.0

        if vel < -1000:
            print(f"[Warning] Input velocity vel={vel:.1f} km/s is very negative.")

        flux_template = np.sum(temp.emission_lines * flux_ratio, axis=1)  # pure line shape

        if sfr is not None:
            if halpha is not None:
                print(f"[Warning] `halpha` is ignored since `sfr` is provided.")
            if z is None or vpec is None:
                raise ValueError("[Error] To use SFR-based calibration, both `z` and `vpec` must be provided.")
            if z < 0:
                raise ValueError(f"[Error] Input redshift z={z:.3f} is negative. Must be z ≥ 0.")

            # SFR-based calibration: first scale, then extinct
            L_Ha = sfr / 7.9e-42  # erg/s
            DL = cosmo.luminosity_distance(z).to('cm').value  # cm
            f_Ha = L_Ha / (4 * np.pi * DL**2)  # erg/s/cm^2
            f_Ha_obs = self._apply_cosmic_dimming(f_Ha, z)  # erg/s/cm^2
            flux_scaled = flux_template * f_Ha_obs * 1e17  # in 1e-17 erg/s/cm^2
            flux_dust = reddening(temp.wave, flux_scaled, ebv=ebv)
            vel = (1 + z) * vpec + z * 3e5  # km/s total LOS velocity

        else:
            # Hα-based calibration: first extinct, then scale
            flux_dust = reddening(temp.wave, flux_template, ebv=ebv)
            flux_dust *= halpha  # scale to match observed Hα flux

        flux_logwave, logLam = log_rebin([np.min(temp.wave), np.max(temp.wave)], flux_dust, velscale=10)[:2]
        flux_broad = self._apply_velocity_dispersion(flux_logwave, logLam, vdisp, temp.FWHM_inst)
        flux_red = self._apply_redshift(logLam, flux_broad, vel)

        self.wave = config.wave
        self.flux = flux_red

    def _apply_velocity_dispersion(self, flux_logwave, logLam, vdisp, FWHM_inst):
        """
        Apply Gaussian broadening based on gas velocity dispersion and instrumental LSF.
        """
        velscale = 10
        sigma_gas = vdisp / velscale
        sigma_LSF = FWHM_inst / (np.exp(logLam)) * 3e5 / velscale

        sigma_dif = np.zeros(len(flux_logwave))
        idx = sigma_gas > sigma_LSF
        sigma_dif[idx] = np.sqrt(sigma_gas ** 2. - sigma_LSF[idx] ** 2.)
        sigma_dif[~idx] = 0.1

        return gaussian_filter1d(flux_logwave, sigma_dif)

    def _apply_redshift(self, logLam, flux_broad, vel):
        """
        Apply redshift to the spectrum based on line-of-sight velocity.
        """
        wave_r = np.exp(logLam) * (1 + vel / 3e5)
        return np.interp(self.config.wave, wave_r, flux_broad)

    def _apply_cosmic_dimming(self, flux, z):
        """
        Apply cosmological surface brightness dimming (∝ 1/(1+z)^4).
        """
        if z is None or z <= 0:
            return flux
        return flux / (1 + z)**4


#############
# AGN Spectra
#############


class AGN_NLR():
    """
    Generate emission line spectra of AGN narrow-line region using template flux ratios.

    The model uses a specified Hα flux (in units of 1e-17 erg/s/cm²) to scale the emission lines.
    It applies optional dust extinction, velocity broadening, and line-of-sight redshift.

    Parameters
    ----------
    config : object
        Configuration object that defines the target output wavelength grid.
    temp : object
        Emission line template object for AGN narrow-line regions.
    halpha : float, optional
        Observed flux of the Hα line in units of 1e-17 erg/s/cm² (default: 100.0).
    logz : float, optional
        Gas-phase metallicity (in log Z/Z☉), used to select emission line ratios (default: 0.0).
        Values outside the range [-2.3, 0.54] will be clipped and a warning will be printed.
    vel : float, optional
        Line-of-sight velocity in km/s for redshifting the emission lines (default: 100.0).
    vdisp : float, optional
        Velocity dispersion in km/s, used to apply Gaussian broadening (default: 120.0).
    ebv : float, optional
        Dust extinction E(B–V), applied using the Calzetti law (default: 0.1).
    """
    def __init__(self, config, temp, halpha=100.0, logz=0.0, vel=100.0, vdisp=120.0, ebv=0.1):
        self.config = config

        if logz < -2.3:
            print(f"[Warning] Input logz={logz:.2f} < -2.3. Clipped to -2.3.")
            logz = -2.3
        elif logz > 0.54:
            print(f"[Warning] Input logz={logz:.2f} > 0.54. Clipped to 0.54.")
            logz = 0.54

        indz = np.argmin(np.abs(logz - temp.logz_grid))
        flux_ratio = temp.flux_ratio[indz, :]

        if vdisp < 0:
            print(f"[Warning] Input vdisp={vdisp:.1f} km/s is invalid. Setting vdisp=0.")
            vdisp = 0.0
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv:.3f} is invalid. Setting ebv=0.")
            ebv = 0.0
        if vel < -1000:
            print(f"[Warning] Input velocity vel={vel:.1f} km/s is very negative.")

        # Line shape and scaling
        flux_template = np.sum(temp.emission_lines * (flux_ratio / flux_ratio[6]), axis=1)
        flux_dust = reddening(temp.wave, flux_template * halpha, ebv=ebv)

        # Rebin and broaden
        flux_logwave, logLam = log_rebin([np.min(temp.wave), np.max(temp.wave)], flux_dust, velscale=10)[:2]
        flux_broad = self._apply_velocity_dispersion(flux_logwave, logLam, vdisp, temp.FWHM_inst)

        # Redshift
        flux_red = self._apply_redshift(logLam, flux_broad, vel)

        self.wave = config.wave
        self.flux = flux_red

    def _apply_velocity_dispersion(self, flux_logwave, logLam, vdisp, FWHM_inst):
        sigma_gas = vdisp / 10
        sigma_LSF = FWHM_inst / (np.exp(logLam)) * 3e5 / 10

        sigma_dif = np.zeros_like(flux_logwave)
        idx = sigma_gas > sigma_LSF
        sigma_dif[idx] = np.sqrt(sigma_gas**2 - sigma_LSF[idx]**2)
        sigma_dif[~idx] = 0.1

        return gaussian_filter1d(flux_logwave, sigma_dif)

    def _apply_redshift(self, logLam, flux_broad, vel):
        wave_r = np.exp(logLam) * (1 + vel / 3e5)
        return np.interp(self.config.wave, wave_r, flux_broad)


class AGN_BLR():
    """
    Generate emission line spectra of AGN broad-line region based on empirical Balmer line ratios.

    The model constructs a set of broad Balmer emission lines (from Hε to Hα) with fixed ratios
    and user-specified Hβ line flux and FWHM. Dust extinction and redshift are applied afterward.

    Parameters
    ----------
    config : object
        Configuration object defining the target wavelength grid (e.g., output wave sampling).
    hbeta_flux : float, optional
        Observed flux of the Hβ broad component in units of 1e-17 erg/s/cm² (default: 100.0).
    hbeta_fwhm : float, optional
        Full width at half maximum (FWHM) of the broad lines in km/s (default: 2000.0).
        Must be positive; if invalid, warning is printed and reset to 2000.
    ebv : float, optional
        Dust extinction E(B–V) applied using the Calzetti law (default: 0.1).
        If negative, will be reset to 0 with warning.
    vel : float, optional
        Line-of-sight velocity in km/s for redshifting the broad lines (default: 0.0).
    lam_range : list, optional
        Wavelength range [min, max] in Å for constructing the rest-frame line spectrum (default: [500, 15000]).
    """
    def __init__(self, config, hbeta_flux=100.0, hbeta_fwhm=2000.0, ebv=0.1,
                 vel=0.0, lam_range=[500, 15000]):

        if hbeta_fwhm <= 0:
            print(f"[Warning] Invalid FWHM input: hbeta_fwhm={hbeta_fwhm:.1f} km/s. Reset to 2000.")
            hbeta_fwhm = 2000.0
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv:.3f} is invalid. Setting ebv=0.")
            ebv = 0.0
        if vel < -10000:
            print(f"[Warning] Very large blueshift: vel={vel:.1f} km/s.")

        # --- Step 1: Build rest-frame wavelength grid ---
        wave_rest = np.arange(lam_range[0], lam_range[1], 0.1)

        # --- Step 2: Define Balmer broad lines ---
        line_names = ['Hε', 'Hδ', 'Hγ', 'Hβ', 'Hα']
        line_waves = [3970.079, 4101.742, 4340.471, 4861.333, 6562.819]
        line_ratios = [0.101, 0.208, 0.405, 1.000, 2.579]  # Ilic et al. (2006)

        # --- Step 3: Construct broadened emission line profiles ---
        emission_lines = []
        for wl in line_waves:
            sigma = hbeta_fwhm / 3e5 * wl  # Gaussian sigma in Å
            profile = SingleEmissinoLine(wave_rest, wl, sigma)
            emission_lines.append(profile)
        emission_lines = np.array(emission_lines).T  # shape: (n_wave, n_lines)

        # --- Step 4: Combine all lines with line ratios and flux scaling ---
        flux_combine = np.dot(emission_lines, line_ratios)  # weighted sum over lines
        flux_calibrated = flux_combine * hbeta_flux         # calibrated total flux

        # --- Step 5: Apply dust extinction ---
        flux_dust = reddening(wave_rest, flux_calibrated, ebv=ebv)

        # --- Step 6: Apply redshift from velocity ---
        redshift = vel / 3e5
        wave_obs = wave_rest * (1 + redshift)
        flux_obs = np.interp(config.wave, wave_obs, flux_dust)

        # --- Final output ---
        self.wave = config.wave
        self.flux = flux_obs


class AGN_FeII():
    """
    Generate Fe II pseudo-continuum emission for AGN using empirical template.

    The Fe II emission is scaled using the observed flux of broad Hβ and the empirical flux ratio
    between the Fe4570 blend and Hβ (commonly denoted R4570). The flux is then corrected for dust extinction
    and redshifted to the observed frame.

    Parameters
    ----------
    config : object
        Configuration object defining the output wavelength grid.
    hbeta_broad : float, optional
        Observed flux of the broad Hβ line in units of 1e-17 erg/s/cm² (default: 100.0).
    r4570 : float, optional
        Flux ratio of Fe4570 blend (4434–4684 Å) to broad Hβ (default: 0.4).
    ebv : float, optional
        Dust extinction E(B–V), applied using the Calzetti law (default: 0.1).
        If negative, will be reset to 0 with warning.
    vel : float, optional
        Line-of-sight velocity in km/s for redshifting the Fe II emission (default: 100.0).
    """
    def __init__(self, config, hbeta_broad=100.0, r4570=0.4, ebv=0.1, vel=100.0):
        self.config = config

        if r4570 < 0:
            print(f"[Warning] Input R4570={r4570:.2f} is invalid. Reset to 0.4.")
            r4570 = 0.4
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv:.3f} is invalid. Setting ebv=0.")
            ebv = 0.0
        if vel < -10000:
            print(f"[Warning] Very large blueshift: vel={vel:.1f} km/s.")

        # --- Step 1: Load Fe II template ---
        filename = data_path + '/data/FeII.AGN.fits'
        hdulist = fits.open(filename)
        data = hdulist[1].data
        wave_rest = data['WAVE']
        flux_template = data['FLUX']  # normalized template

        # --- Step 2: Flux calibration using R4570 ---
        flux_Fe4570_ref = 100.0  # Reference value of Fe4570 in template (unitless)
        flux_Fe4570_target = hbeta_broad * r4570
        scale_factor = flux_Fe4570_target / flux_Fe4570_ref
        flux_scaled = flux_template * scale_factor  # calibrated spectrum

        # --- Step 3: Dust extinction ---
        flux_dust = reddening(wave_rest, flux_scaled, ebv=ebv)

        # --- Step 4: Redshift ---
        wave_obs = wave_rest * (1 + vel / 3e5)
        flux_obs = np.interp(config.wave, wave_obs, flux_dust)

        # --- Final output ---
        self.wave = config.wave
        self.flux = flux_obs


class AGN_Powerlaw():
    """
    Generate a power-law continuum spectrum for AGN.

    The AGN continuum is modeled as a simple power-law:
        F(λ) ∝ λ^α
    The flux is normalized using the AB magnitude at 5100 Å (m5100),
    then corrected for dust extinction and redshift.

    Parameters
    ----------
    config : object
        Configuration object defining the output wavelength grid.
    m5100 : float, optional
        AB magnitude at rest-frame 5100 Å used for flux normalization (default: 18.0).
    alpha : float, optional
        Power-law index such that F(λ) ∝ λ^α (default: -1.5).
    vel : float, optional
        Line-of-sight velocity in km/s for redshifting the spectrum (default: 100.0).
    ebv : float, optional
        Dust extinction E(B–V), applied using the Calzetti law (default: 0.1).
        If negative, it will be reset to 0 with warning.

    Output
    ------
    self.wave : ndarray
        Observed wavelength grid (same as config.wave).
    self.flux : ndarray
        Flux density in units of erg/s/cm²/Å at each wavelength.
    """
    def __init__(self, config, m5100=18.0, alpha=-1.5, vel=100.0, ebv=0.1):
        self.config = config

        # --- Validate inputs ---
        if m5100 <= 0 or m5100 > 40:
            print(f"[Warning] Unusual m5100 input: m5100={m5100:.2f} mag")
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv:.3f} is invalid. Setting ebv=0.")
            ebv = 0.0
        if vel < -10000:
            print(f"[Warning] Very large blueshift: vel={vel:.1f} km/s.")

        # --- Step 1: Generate rest-frame power-law spectrum ---
        wave_rest = np.linspace(1000, 20000, 10000)  # in Å
        flux_raw = wave_rest ** alpha

        # --- Step 2: Flux calibration at 5100 Å region (5050–5150 Å) ---
        flux_calibrate = calibrate(wave_rest, flux_raw, m5100, filtername='5100')

        # --- Step 3: Apply dust extinction ---
        flux_dust = reddening(wave_rest, flux_calibrate, ebv=ebv)

        # --- Step 4: Apply redshift from velocity ---
        redshift = vel / 3e5
        wave_obs = wave_rest * (1 + redshift)
        flux_obs = np.interp(config.wave, wave_obs, flux_dust)

        # --- Final output ---
        self.wave = config.wave
        self.flux = flux_obs


class AGN_PhysicalModel():
    """
    Generate a physical AGN spectrum including:
        - Power-law continuum
        - Broad emission lines (BLR)
        - Narrow emission lines (NLR)
        - FeII pseudo-continuum

    Parameters
    ----------
    config : object
        Configuration object defining the target output wavelength grid.
    nlr_template : object
        Emission line template for narrow-line regions (e.g. AGN_NLR or HII).
    bhmass : float
        Black hole mass in solar masses. Must be > 0.
    edd_ratio : float
        Eddington ratio (L_bol / L_edd). Must be > 0.
    z : float
        Redshift used for luminosity distance.
    logz : float, optional
        Gas-phase metallicity (log Z/Z☉) for the NLR (default: 0.0).
    vpec : float, optional
        Peculiar velocity in km/s for redshift shift (default: 0).
    ebv : float, optional
        Dust extinction E(B–V) applied to all components (default: 0.1).
    agn_type : int, optional
        AGN type indicator (1=type 1, 2=type 2), controls whether BLR/FeII/powerlaw are included (default: 1).
    """
    def __init__(self, config, nlr_template, bhmass, edd_ratio, z,
                 logz=0.0, vpec=0.0, ebv=0.1, agn_type=1):

        if bhmass <= 0:
            raise ValueError("[Error] Black hole mass (bhmass) must be > 0 [solar masses].")
        if edd_ratio <= 0:
            raise ValueError("[Error] Eddington ratio (edd_ratio) must be > 0 [L_bol / L_edd].")
        if z < 0:
            raise ValueError("[Error] Redshift (z) must be ≥ 0.")

        self.config = config
        self.wave = config.wave.copy()
        self.flux = np.zeros_like(self.wave)

        # Determine luminosity distance from redshift
        dist = cosmo.luminosity_distance(z).to('Mpc').value

        # Compute L_bol and L_5100 in erg/s
        L_edd = 1.26e38 * bhmass              # erg/s
        L_bol = L_edd * edd_ratio             # erg/s
        L_5100 = L_bol / 10.9                 # erg/s
        M_5100 = 4.86 - 2.5 * np.log10(L_5100 / 3.828e33)  # AB mag at 10 pc assuming L_sun = 3.828e33 erg/s
        m_5100 = M_5100 + 5 * np.log10(dist * 1e5)

        # Estimate Hβ broad luminosity (erg/s) from L_5100 using Greene & Ho (2005)
        L_Hb = 1.425e42 * (L_5100 / 1e44) ** 1.133
        F_Hb = L_Hb / (4 * np.pi * (dist * 3.086e24)**2) * 1e17  # in 1e-17 erg/s/cm^2

        # Estimate FWHM from L5100 and M_BH using Greene & Ho (2005):
        # M_BH ≈ 4.4e6 × (FWHM / 1000)^2 × (L_5100 / 1e44)^0.64 => invert this
        L_5100_44 = L_5100 / 1e44
        FWHM_Hb = 1000 * np.sqrt(bhmass / (4.4e6 * L_5100_44**0.64))

        # FeII flux (empirical scaling, Shen+2011)
        FeII_Hb_ratio = min(2.0, 1.0 + 5.0 * edd_ratio)

        # Hα broad flux assuming Case B recombination ratio Hα/Hβ = 2.579
        F_Ha = F_Hb * 2.579

        # NLR Hα flux: L_Hα_NLR ∝ L_bol^0.8 (approximate scaling)
        L_Ha_NLR = 1e40 * (L_bol / 1e44)**0.8  # erg/s
        F_Ha_NLR = L_Ha_NLR / (4 * np.pi * (dist * 3.086e24)**2) * 1e17

        # NLR velocity dispersion from M–σ relation
        vdisp_nlr = 200 * (bhmass / 2.1e8)**0.177

        # Total velocity shift from redshift and peculiar velocity
        vel_obs = vpec + z * 3e5

        # --- NLR ---
        if agn_type in [1, 2]:
            nlr = AGN_NLR(config, nlr_template, halpha=F_Ha_NLR, logz=logz,
                          vel=vel_obs, vdisp=vdisp_nlr, ebv=ebv)
            self.flux += nlr.flux

        # --- BLR + FeII ---
        if agn_type == 1:
            blr = AGN_BLR(config, hbeta_flux=F_Hb, hbeta_fwhm=FWHM_Hb, vel=vel_obs, ebv=ebv)
            feii = AGN_FeII(config, hbeta_broad=F_Hb * FeII_Hb_ratio, vel=vel_obs, ebv=ebv)
            self.flux += blr.flux + feii.flux

        # --- Power-law continuum ---
        if agn_type == 1:
            pl = AGN_Powerlaw(config, m5100=m_5100, vel=vel_obs, ebv=ebv)
            self.flux += pl.flux


        # Apply surface brightness dimming
        self.flux /= (1 + z)**4


######################################
# Stellar Population Continuum Spectra
######################################


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
    age = float(s[s.find("T") + 1:s.find("_iPp0.00_baseFe.fits")])
    metal = s[s.find("Z") + 1:s.find("T")]
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
    def __init__(self, config, velscale=50,
                 pathname=data_path + '/data/EMILES/Ech*_baseFe.fits',
                 normalize=False):

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
        lam_range_temp = h2['CRVAL1'] + np.array([0, h2['CDELT1'] * (h2['NAXIS1'] - 1)])
        sspNew, log_lam_temp = log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]

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
        Emile_FWHM[np.where(Emile_wave >= 8950)] = 60 * 2.35 / 3.e5 * Lwave  # sigma=60km/s at lambda > 8950

        LSF = Emile_FWHM

        # Here we make sure the spectra are sorted in both [M/H] and Age
        # along the two axes of the rectangular grid of templates.
        for j, age in enumerate(ages):
            for k, metal in enumerate(metals):
                p = all.index((age, metal))
                hdu = fits.open(files[p])
                ssp = hdu[0].data
                sspNew = log_rebin(lam_range_temp, ssp, velscale=velscale)[0]
                if normalize:
                    sspNew /= np.mean(sspNew)
                templates[:, j, k] = sspNew
                age_grid[j, k] = age
                metal_grid[j, k] = metal

        self.templates = templates #/ np.median(templates)  # Normalize by a scalar
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
        # fage = self.age_grid[:,0]

        fMs = np.zeros((n_ages, n_metal))

        Metal, Age, Ms = readcol(massfile, usecols=(2, 3, 6))
        for i in range(n_metal):
            for j in range(self.n_ages):
                target_z = self.metal_grid[j, i]
                target_age = self.age_grid[j, i]
                dist = (np.log10(Age) - np.log10(target_age))**2 + (Metal - target_z)**2
                locmin = np.argmin(dist)
                fMs[j, i] = Ms[locmin]

        return fMs


class StellarContinuum:
    """
    Generate galaxy stellar continuum spectra using input Star Formation History (SFH)
    and Chemical Enrichment History (CEH), based on a grid of pre-computed Simple Stellar Population (SSP) templates.

    Two working modes are supported:

    1. Magnitude-calibrated mode:
    If `mag` is provided, the output spectrum is flux-calibrated using the given SDSS r-band magnitude.
    In this case, the input `sfh` only determines the *shape* of the star formation history (normalized SFH profile),
    and does not affect the flux amplitude.

    2. Physically-calibrated mode:
    If `mag` is None, the spectrum is calibrated based on the physical star formation rate (SFR) in units of M☉/yr.
    The total formed stellar mass is integrated from the input SFH and used to scale the template.
    This mode requires both `z` (cosmological redshift, for cosmic dimming) and `vpec` (peculiar velocity),
    which together define the line-of-sight redshift.

    Input flexibility:

    - Both `sfh` and `ceh` can be scalars or (N, 2) arrays.
    - Scalar inputs indicate a single SSP population and are only allowed in magnitude-calibrated mode.
    - In physically-calibrated mode, both `sfh` and `ceh` must be arrays, and the SFR values must be in physical units (M☉/yr).

    Parameters
    ----------
    config : object
        Configuration object describing the format of the output spectrum (e.g., wavelength grid).
    template : object
        Template object containing a grid of SSP spectra (wavelength × age × metallicity).
    mag : float, optional
        Observed SDSS r-band magnitude for flux calibration (default 15.0). 
        Recommended range: 8–26. If None, the model enters physically-calibrated mode.
    sfh : array-like or float, optional
        Star Formation History. Either a scalar age (Gyr), or an (N,2) array of [age (Gyr), SFR (M☉/yr)].
        In mag-mode, only the normalized SFH shape is used. In physical mode, the total SFR determines luminosity.
    ceh : array-like or float, optional
        Chemical Enrichment History. Either a scalar [Fe/H] (dex), or an (N,2) array of [age (Gyr), [Fe/H] (dex)].
    vel : float, optional
        Line-of-sight velocity in km/s, used to redshift the spectrum (default 100). Required in mag-mode.
    vdisp : float, optional
        Stellar velocity dispersion in km/s, used to apply spectral broadening (default 100).
    ebv : float, optional
        Dust extinction E(B–V), applied using the Calzetti law (default 0.1).
    z : float, optional
        Cosmological redshift. Required in physically-calibrated mode (used for luminosity distance and cosmic dimming).
    vpec : float, optional
        Galaxy peculiar velocity in km/s. Required in physically-calibrated mode.
    """
    def __init__(self, config, template, mag=15,
                 sfh=1.0, ceh=-0.3, vel=100, vdisp=100.0, ebv=0.1,
                 z=None, vpec=None):

        # Check magnitude range
        self.mag = mag
        if mag is not None:
            if mag < 8 or mag > 26:
                print(f"[Warning] Input mag={mag:.2f} is outside the 8–26 recommended range.")
            self.use_flux_calibration = mag is not None

        # Check and fix velocity dispersion
        if vdisp < 0:
            print(f"[Warning] Input vdisp={vdisp} km/s is invalid. Setting vdisp=0.")
            vdisp = 0.0

        # Check and fix extinction
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv} is invalid. Setting ebv=0.")
            ebv = 0.0

        # Check velocity (no correction but warning if very negative)
        if vel < -1000:
            print(f"[Warning] Input velocity vel={vel} km/s is very negative.")

        self.config = config
        self.template = template

        # Mode determination
        self.use_flux_calibration = mag is not None

        if self.use_flux_calibration:
            if vel is None:
                raise ValueError("In mag-based mode, 'vel' must be provided.")
            self.vel = vel
        else:
            if z is None or vpec is None:
                raise ValueError("In physical modeling mode, both 'z' and 'vpec' must be provided.")
            if z < 0:
                raise ValueError(f"Invalid redshift z={z}. Redshift must be non-negative.")
            self.z = z
            self.vel = vpec + z * 3e5  # km/s

        # Load template grids
        self.ages = template.age_grid[:, 0]  # Gyr
        self.metals = template.metal_grid[0, :]  # [Fe/H]
        fMs = template.fmass_ssp()

        # Detect input type: scalar or array
        if self.use_flux_calibration:
            self.sfh_array, is_scalar_sfh = self._validate_input(sfh, 'sfh', 0.05, 13.7, allow_negative_y=False)
            self.ceh_array, is_scalar_ceh = self._validate_input(ceh, 'ceh', 0.05, 13.7, -2.32, 0.22)
        else:
            self.sfh_array, is_scalar_sfh = self._validate_input(sfh, 'sfh', 0.05, 13.7, allow_negative_y=False)
            self.ceh_array, is_scalar_ceh = self._validate_input(ceh, 'ceh', 0.05, 13.7, -2.32, 0.22)
            if is_scalar_sfh or is_scalar_ceh:
                raise ValueError("In physical modeling mode, both SFH and CEH must be arrays with physical units.")

        if is_scalar_sfh and is_scalar_ceh:
            # Single SSP mode
            age_val = self.sfh_array
            if age_val < 0.05:
                print(f"[Warning] SFH age {age_val} Gyr < 0.05. Clipped to 0.05.")
                age_val = 0.05
            elif age_val > 13.7:
                print(f"[Warning] SFH age {age_val} Gyr > 13.7. Clipped to 13.7.")
                age_val = 13.7

            feh_val = self.ceh_array
            if feh_val < -2.32:
                print(f"[Warning] [Fe/H] {feh_val} < -2.32. Clipped to -2.32.")
                feh_val = -2.32
            elif feh_val > 0.22:
                print(f"[Warning] [Fe/H] {feh_val} > 0.22. Clipped to 0.22.")
                feh_val = 0.22

            i_age = np.argmin(np.abs(self.ages - age_val))
            i_feh = np.argmin(np.abs(self.metals - feh_val))

            Stellar = template.templates[:, i_age, i_feh]

            self.mass_weighted_age = self.ages[i_age]
            self.mass_weighted_feh = self.metals[i_feh]
            self.mass_formed = 1.0
            self.weights = None

        else:
            # Composite mode
            if is_scalar_sfh:
                self.sfh_array = np.column_stack([self.ages, np.full_like(self.ages, self.sfh_array)])
            if is_scalar_ceh:
                self.ceh_array = np.column_stack([self.ages, np.full_like(self.ages, self.ceh_array)])

            weights = self._compute_weights(fMs)
            self.weights = weights

            Stellar = np.einsum('ijk,jk->i', template.templates, weights)

            self.mass_weighted_age = np.sum(weights * self.ages[:, None]) / np.sum(weights)
            self.mass_weighted_feh = np.sum(weights * self.metals[None, :]) / np.sum(weights)
            self.mass_formed = np.sum(weights / fMs)

        # Observational effects
        log_wave = template.log_lam_temp
        rest_wave = np.exp(log_wave)

        flux0 = self._apply_velocity_dispersion(Stellar, rest_wave, vdisp, template.LSF, template.velscale)
        flux0 = self._apply_dust_extinction(rest_wave, flux0, ebv)

        obs_wave, obs_flux = self._apply_redshift(rest_wave, flux0)

        flux = self._apply_flux_calibration(obs_wave, obs_flux, config.wave)

        self.wave = config.wave
        self.flux = flux

    def _validate_input(self, input_data, name, age_min=0.05, age_max=13.7, y_min=None, y_max=None, allow_negative_y=True):
        """
        Validate and preprocess SFH or CEH input.

        Parameters
        ----------
        input_data : float or array-like
            Either a scalar value (for single SSP) or a 2D array (N, 2) of [age, value].
        name : str
            Name for logging and error messages ("sfh" or "ceh").
        age_min, age_max : float
            Allowed range for stellar population age in Gyr.
        y_min, y_max : float or None
            Optional clipping range for second column (SFR or [Fe/H]).
        allow_negative_y : bool
            Whether to allow negative values in the second column.

        Returns
        -------
        data_out : array or float
            Processed data or scalar.
        is_scalar : bool
            Whether the input was scalar.
        """
        if np.isscalar(input_data):
            return input_data, True

        data = np.asarray(input_data)
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError(f"{name.upper()} must be a scalar or (N,2) array [age, value]")

        # Check age
        if np.any(data[:, 0] < age_min) or np.any(data[:, 0] > age_max):
            print(f"[Warning] Some ages in {name.upper()} are outside range {age_min}-{age_max} Gyr.")

        # Check y values
        if not allow_negative_y and np.any(data[:, 1] < 0):
            print(f"[Warning] {name.upper()} has negative values in the second column, setting to zero.")
            data[:, 1] = np.maximum(data[:, 1], 0)

        if y_min is not None and np.any(data[:, 1] < y_min):
            print(f"[Warning] Some values in {name.upper()} < {y_min}, clipping.")
            data[:, 1] = np.maximum(data[:, 1], y_min)

        if y_max is not None and np.any(data[:, 1] > y_max):
            print(f"[Warning] Some values in {name.upper()} > {y_max}, clipping.")
            data[:, 1] = np.minimum(data[:, 1], y_max)

        return data, False

    def _compute_weights(self, mass_ratios):
        """
        Compute stellar mass weights over the age-metallicity grid.

        Parameters
        ----------
        mass_ratios : ndarray
            Mass retention fraction for each (age, metallicity) bin.

        Returns
        -------
        weight_matrix : ndarray
            Mass-weighted matrix [n_age × n_metallicity], in units of M☉.
        """
        sfh_interp = interp1d(self.sfh_array[:, 0], self.sfh_array[:, 1],
                              kind='linear', bounds_error=False, fill_value=0)
        ceh_interp = interp1d(self.ceh_array[:, 0], self.ceh_array[:, 1],
                              kind='linear', bounds_error=False, fill_value='extrapolate')

        age_edges = np.zeros(len(self.ages) + 1)
        age_edges[1:-1] = 0.5 * (self.ages[:-1] + self.ages[1:])
        age_edges[0] = self.ages[0] - (self.ages[1] - self.ages[0]) * 0.5
        age_edges[-1] = self.ages[-1] + (self.ages[-1] - self.ages[-2]) * 0.5
        dt = np.diff(age_edges)  # Gyrs
        age_bin_center = 0.5 * (age_edges[:-1] + age_edges[1:])

        sfr = sfh_interp(age_bin_center)  # in Msol/yr
        mass_formed_dt = sfr * dt * 1e9  # Msol per bin (Gyr -> yr)

        weight_matrix = np.zeros_like(mass_ratios)

        for i, (age, mass_bin) in enumerate(zip(age_bin_center, mass_formed_dt)):
            if mass_bin <= 0:
                continue
            feh = ceh_interp(age)
            j = np.argmin(np.abs(self.metals - feh))
            weight_matrix[i, j] = mass_bin * mass_ratios[i, j]  # Mformed * fraction

        total_mass = np.sum(weight_matrix)
        return weight_matrix

    def _apply_velocity_dispersion(self, flux, wave, vdisp, LSF, velscale):
        """
        Apply velocity dispersion broadening to input flux.

        Parameters
        ----------
        flux : array-like
            Input stellar flux.
        wave : array-like
            Rest-frame wavelength grid (Å).
        vdisp : float
            Galaxy stellar velocity dispersion (km/s).
        LSF : array-like
            Instrumental FWHM in Å at each wavelength.
        velscale : float
            Velocity per pixel (km/s).

        Returns
        -------
        flux_broadened : array-like
            Broadened flux.
        """
        sigma_LSF_A = LSF / 2.355
        sigma_LSF_kms = (sigma_LSF_A / wave) * 3e5  # convert to km/s
        sigma_gal_pix = vdisp / velscale
        sigma_LSF_pix = sigma_LSF_kms / velscale

        if sigma_gal_pix > 0:
            delta_sigma = np.zeros_like(flux)
            mask = sigma_gal_pix > sigma_LSF_pix
            delta_sigma[mask] = np.sqrt((sigma_gal_pix**2 - sigma_LSF_pix[mask]**2))
            delta_sigma[~mask] = 0.1
            return gaussian_filter1d(flux, delta_sigma)
        return flux

    def _apply_dust_extinction(self, wave, flux, ebv):
        """
        Apply Calzetti dust attenuation curve.

        Parameters
        ----------
        wave : array-like
            Wavelength grid (Å).
        flux : array-like
            Input flux.
        ebv : float
            E(B-V) value.

        Returns
        -------
        flux_dust : array-like
            Extincted flux.
        """
        if ebv is not None and ebv > 0:
            return reddening(wave, flux, ebv=ebv)
        return flux

    def _apply_redshift(self, wave, flux):
        """
        Shift wavelength and flux to observer frame using LOS velocity.

        Parameters
        ----------
        wave : array-like
            Rest-frame wavelength (Å).
        flux : array-like
            Flux in rest frame.

        Returns
        -------
        obs_wave : array
            Observer-frame wavelength.
        obs_flux : array
            Same flux (no energy shift).
        """
        z = self.vel / 3e5
        return wave * (1 + z), flux
    
    def _apply_flux_calibration(self, obs_wave, obs_flux, target_wave):
        """
        Calibrate output flux either by observed magnitude (mag-mode)
        or using physical luminosity and distance (physical-mode).

        Parameters
        ----------
        obs_wave : array-like
            Wavelength after redshift.
        obs_flux : array-like
            Raw output flux (Lsun/Å).
        target_wave : array-like
            Target wavelength grid.

        Returns
        -------
        flux_calibrated : array-like
            Final flux in units of 1e-17 erg/s/Å/cm².
        """
        if self.use_flux_calibration:
            # Observational mode using synthetic magnitude
            flux = np.interp(target_wave, obs_wave, obs_flux)
            flux = calibrate(target_wave, flux, self.mag, filtername='SLOAN_SDSS.r')
        else:
            # Physical mode: convert Lsun/Å to erg/s/Å/cm²
            Lsun_erg_s = 3.828e33
            flux = obs_flux * Lsun_erg_s  # erg/s/Å
            DL = cosmo.luminosity_distance(self.z).to('cm').value
            flux /= (4 * np.pi * DL**2)  # erg/s/Å/cm²
            flux /= (1 + self.z)  # bandwidth dilation effect
            flux *= 1e17  # Convert to 1e-17 erg/s/Å/cm²
            flux = np.interp(target_wave, obs_wave, flux)
        return flux
    
    def _apply_cosmic_dimming(self, flux, z):
        """
        Apply cosmological surface brightness dimming (∝ 1/(1+z)^4).

        Parameters
        ----------
        flux : array-like
            Input flux.
        z : float
            Cosmological redshift.

        Returns
        -------
        flux_diluted : array-like
            Flux after dimming correction.
        """
        if z is None or z <= 0:
            return flux
        return flux / (1 + z)**4

#####################
# Single Star Spectra
#####################


class SingleStarTemplate:
    """
    Class for handling single stellar templates from different libraries.

    Supports observed XSL templates and theoretical Munari+2005 templates.

    Parameters
    ----------
    config : object
        Configuration object containing instrument settings.
    velscale : float, optional
        Velocity scale in km/s per pixel. Default is 20.
    template : str, optional
        Template library to use: 'Munari2005' (default) or 'XSL'. Case-insensitive.
    """

    def __init__(self, config, velscale=20, template="Munari2005"):
        template = template.lower()
        data_dir = get_path()

        if template == "xsl":
            filename = path.join(data_dir, "data/Starlib.XSL.fits")
            hdulist = fits.open(filename)
            lam = hdulist[1].data['Wave']
            flux = hdulist[2].data
            par = hdulist[3].data

            lam_range_temp = np.array([lam[0], lam[-1]])
            TemNew, log_lam_temp = log_rebin(lam_range_temp, flux[1, :], velscale=velscale)[:2]

            Temp_wave = np.exp(log_lam_temp)
            Temp_FWHM = np.zeros(len(log_lam_temp))
            Temp_FWHM[Temp_wave < 5330] = 13 * 2.35 / 3e5 * Temp_wave[Temp_wave < 5330]
            Temp_FWHM[(Temp_wave >= 5330) & (Temp_wave < 9440)] = 11 * 2.35 / 3e5 * Temp_wave[(Temp_wave >= 5330) & (Temp_wave < 9440)]
            Temp_FWHM[Temp_wave >= 9440] = 16 * 2.35 / 3e5 * Temp_wave[Temp_wave >= 9440]

            valid_indices = []
            temp_list = []
            teff_list, logg_list, feh_list = [], [], []

            for i in range(par.size):
                raw_flux = flux[i, :]
                if not np.all(np.isfinite(raw_flux)) or np.all(raw_flux == 0):
                    continue
                raw_flux = np.where(raw_flux <= 0, 1e-6, raw_flux)
                temp0 = log_rebin(lam_range_temp, raw_flux, velscale=velscale)[0]
                tempNew = temp0 / np.mean(temp0)
                temp_list.append(tempNew)

                valid_indices.append(i)
                teff_list.append(par['Teff'][i])
                logg_list.append(par['logg'][i])
                feh_list.append(par['FeH'][i])

            temp = np.array(temp_list).T
            self.teff_grid = np.array(teff_list)
            self.logg_grid = np.array(logg_list)
            self.feh_grid = np.array(feh_list)

        elif template == "munari2005":
            filename = path.join(data_dir, "data/Starlib.Munari2005.fits")
            hdulist = fits.open(filename)
            lam = hdulist[1].data['Wave']
            flux = hdulist[2].data
            par = hdulist[3].data

            lam_range_temp = np.array([lam[0], lam[-1]])
            TemNew, log_lam_temp = log_rebin(lam_range_temp, flux[1, :], velscale=velscale)[:2]

            Temp_wave = np.exp(log_lam_temp)
            Temp_FWHM = np.full_like(Temp_wave, fill_value=1.0 * 2.35 / 3e5 * Temp_wave)

            valid_indices = []
            temp_list = []
            teff_list, logg_list, feh_list = [], [], []

            for i in range(par.size):
                raw_flux = flux[i, :]
                if not np.all(np.isfinite(raw_flux)) or np.all(raw_flux == 0):
                    continue
                raw_flux = np.where(raw_flux <= 0, 1e-6, raw_flux)
                temp0 = log_rebin(lam_range_temp, raw_flux, velscale=velscale)[0]
                tempNew = temp0 / np.mean(temp0)
                temp_list.append(tempNew)

                valid_indices.append(i)
                teff_list.append(par['Teff'][i])
                logg_list.append(par['logg'][i])
                feh_list.append(par['FeH'][i])

            temp = np.array(temp_list).T
            self.teff_grid = np.array(teff_list)
            self.logg_grid = np.array(logg_list)
            self.feh_grid = np.array(feh_list)

        else:
            raise ValueError(f"Unknown template: {template}. Use 'Munari2005' or 'XSL'.")

        self.templates = temp
        self.log_lam_temp = log_lam_temp
        self.LSF = Temp_FWHM
        self.velscale = velscale


class SingleStar:
    """
    Generate a single stellar spectrum using input stellar atmospheric parameters.

    This class interpolates the closest-matching stellar template from a grid defined by effective temperature (Teff),
    surface gravity (logg), and metallicity ([Fe/H]), and returns a redshifted, extincted, and magnitude-calibrated spectrum.

    Parameters
    ----------
    config : object
        Configuration object describing the wavelength grid and filter system.
    template : object
        Template object from SingleStarTemplate class, containing wavelength, flux grid, and parameter arrays.
    mag : float, optional
        Apparent SDSS r-band magnitude for flux calibration (default is 15.0). Range: 8–26.
    teff : float, optional
        Effective temperature in Kelvin. Will be clipped to template bounds if outside range.
    logg : float, optional
        Surface gravity in log10(cm/s^2). Will be clipped to template bounds if outside range.
    feh : float, optional
        Metallicity [Fe/H] in dex. Will be clipped to template bounds if outside range.
    vel : float, optional
        Line-of-sight velocity in km/s (default is 100). Used for Doppler shifting.
    ebv : float, optional
        Dust extinction E(B-V). Negative values are reset to 0. Default is 0.0.
    """

    def __init__(self, config, template, mag=15.0, teff=10000.0, logg=4.5, feh=0.0, vel=100.0, ebv=0.0):
        self.wave = config.wave

        # Check magnitude range
        if mag < 8 or mag > 26:
            print(f"[Warning] Input mag={mag:.2f} is outside the 8–26 recommended range.")
        self.mag = mag

        # Check velocity
        if vel < -1000:
            print(f"[Warning] Input velocity vel={vel} km/s is very negative.")
        self.vel = vel

        # Check extinction
        if ebv < 0:
            print(f"[Warning] Input ebv={ebv} is invalid. Setting ebv=0.")
            ebv = 0.0
        self.ebv = ebv

        # Get parameter ranges from template
        all_teff = template.teff_grid
        all_logg = template.logg_grid
        all_feh = template.feh_grid

        teff_min, teff_max = np.min(all_teff), np.max(all_teff)
        logg_min, logg_max = np.min(all_logg), np.max(all_logg)
        feh_min, feh_max = np.min(all_feh), np.max(all_feh)

        # --- Check and clip stellar parameters ---
        if teff < teff_min or teff > teff_max:
            print(f"[Warning] Input Teff={teff:.1f} K is outside template range [{teff_min:.1f}, {teff_max:.1f}] K. Clipped.")
            teff = np.clip(teff, teff_min, teff_max)
        if logg < logg_min or logg > logg_max:
            print(f"[Warning] Input logg={logg:.2f} is outside template range [{logg_min:.2f}, {logg_max:.2f}]. Clipped.")
            logg = np.clip(logg, logg_min, logg_max)
        if feh < feh_min or feh > feh_max:
            print(f"[Warning] Input FeH={feh:.2f} is outside template range [{feh_min:.2f}, {feh_max:.2f}]. Clipped.")
            feh = np.clip(feh, feh_min, feh_max)

        # --- Match closest template ---
        wave = np.exp(template.log_lam_temp)
        temp_fluxes = template.templates

        # Priority 1: match Teff
        mask_teff = (np.abs(all_teff - teff) == np.min(np.abs(all_teff - teff)))
        teff_matched_indices = np.where(mask_teff)[0]

        # Priority 2: match logg within matched teff
        logg_subset = all_logg[teff_matched_indices]
        idx_logg_in_subset = np.argmin(np.abs(logg_subset - logg))
        logg_best = logg_subset[idx_logg_in_subset]
        logg_matched_indices = teff_matched_indices[idx_logg_in_subset]

        # Priority 3: match feh within teff+logg matched
        teff_logg_mask = (all_teff == all_teff[logg_matched_indices]) & (all_logg == logg_best)
        feh_subset = all_feh[teff_logg_mask]
        idx_feh = np.argmin(np.abs(feh_subset - feh))
        final_idx = np.where(teff_logg_mask)[0][idx_feh]
        feh_best = all_feh[final_idx]

        # Warn if final match differs too much from input
        if abs(logg_best - logg) > 0.5 or abs(feh_best - feh) > 0.5:
            print(f"[Warning] No close match found for logg={logg:.2f} or FeH={feh:.2f}. Using closest available template: logg={logg_best:.2f}, FeH={feh_best:.2f}.")

        starspec = temp_fluxes[:, final_idx]

        # Dust extinction
        if self.ebv > 0:
            starspec = reddening(wave, starspec, ebv=self.ebv)

        # Redshift
        redshift = self.vel / 3e5
        wave_r = wave * (1 + redshift)

        # Interpolate to target wavelength grid
        flux = np.interp(self.wave, wave_r, starspec)

        # Calibration
        flux = calibrate(self.wave, flux, mag, filtername="SLOAN_SDSS.r")

        self.flux = flux
