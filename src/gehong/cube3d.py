import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d
from .spec1d import *
import astropy.wcs

class Cube3D():
    """
    Class of 3-dimentional spectral cube
    """
    def __init__(self, config, stellar_map = None, gas_map = None):
        
        self.config= config
        
        self.nx    = config.nx
        self.ny    = config.ny
        self.dpix  = config.dpix
        self.fov_x = config.fov_x
        self.fov_y = config.fov_y
        
        self.wave  = config.wave
        self.nz    = len(self.wave)
        self.wave0 = np.min(self.wave)
        self.inst_fwhm = config.inst_fwhm
        
        self.flux  = np.zeros((self.nx,self.ny,self.nz))
        
        self.stellar_map = stellar_map
        self.gas_map     = gas_map
        
    def make_cube(self, stellar_tem = None, hii_tem = None):
        
        for i in range(self.nx):
            for j in range(self.ny):
                if self.stellar_map is not None:
                    ss = StellarContinuum(self.config, stellar_tem, mag = self.stellar_map.mag[i,j], 
                                          age = self.stellar_map.age[i,j], feh = self.stellar_map.feh[i,j], 
                                          vel = self.stellar_map.vel[i,j], vdisp = self.stellar_map.vdisp[i,j], 
                                          ebv = self.stellar_map.ebv[i,j])
                if self.gas_map is not None:
                    gg = HII_Region(self.config, hii_tem, halpha = self.gas_map.halpha[i,j], 
                                    logz = self.gas_map.zh[i,j], vel = self.gas_map.vel[i,j], 
                                    vdisp = self.gas_map.vdisp[i,j], ebv = self.gas_map.ebv[i,j])
                    self.flux[i,j,:] = ss.flux + gg.flux
                else:
                    self.flux[i,j,:] = ss.flux
                    
    def wcs_info(self):
        wcs = fits.Header()
        wcs_dict = {'CTYPE1': 'WAVE    ', 
                    'CUNIT1': 'Angstrom', 
                    'CDELT1': self.inst.dlam, 
                    'CRPIX1': 1, 
                    'CRVAL1': np.min(self.wave), 
                    'CTYPE2': 'RA---TAN', 
                    'CUNIT2': 'deg', 
                    'CDELT2': self.dpix / 3600., 
                    'CRPIX2': np.round(self.ny / 2.), 
                    'CRVAL2': 0.5, 
                    'CTYPE3': 'DEC--TAN', 
                    'CUNIT3': 'deg', 
                    'CDELT3': self.dpix / 3600., 
                    'CRPIX3': np.round(self.nx / 2.), 
                    'CRVAL3': 1, 
                    'BUNIT' : '10**(-17)*erg/s/cm**2/Angstrom'}
        input_wcs = astropy.wcs.WCS(wcs_dict)
        self.wcs_header = input_wcs.to_header()
        
    def insert_spec(self, spec, dx = 0, dy = 0):
        x0 = np.int(np.round(self.inst.nx / 2.))
        y0 = np.int(np.round(self.inst.ny / 2.))
        self.flux[x0 + dx, y0 + dy, :] = self.flux[x0 + dx, y0 + dy, :] + spec.flux
        
    def savefits(self, filename, path = './'): 
        
        hdr = fits.Header()
        
        hdr['FILETYPE'] = 'SCICUBE'
        hdr['CODE']     = 'CSST-IFS-GEHONG'
        hdr['VERSION']  = '0.0.1'
        
        hdr['OBJECT']   = 'NGC1234'
        hdr['RA']       = 0.0
        hdr['DEC']      = 0.0
        
        hdu0 = fits.PrimaryHDU(header = hdr)
        
        self.wcs_info()
        hdr = self.wcs_header
        hdu1 = fits.ImageHDU(self.flux, header = hdr)
        
        # Output
        hdulist = fits.HDUList([hdu0, hdu1])
        hdulist.writeto(path + filename, overwrite = True)