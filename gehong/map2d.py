from __future__ import division

import scipy.special as sp
import numpy as np
from astropy.io import fits
from skimage.transform import resize

def Sersic2D(x, y, mag = 12, r_eff = 1, n = 2, ellip = 0.5, 
             theta = 0, x_0 = 0, y_0 = 0, pixelscale = 0.01):
    
    # Produce Sersic profile
    bn = sp.gammaincinv(2. * n, 0.5)
    a, b = r_eff, (1 - ellip) * r_eff
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
    x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
    z = (abs(x_maj / a) ** 2 + abs(x_min / b)  ** 2) ** (1 / 2)
    profile = np.exp(-bn * (z ** (1 / n) - 1))
    
    # Normalization
    integral = a * b * 2 * np.pi * n * np.exp(bn) / (bn ** (2 * n)) * sp.gamma(2 * n)
    prof_norm = profile / integral * pixelscale
    
    # Calibration
    total_flux = 10. ** ((22.5 - mag) * 0.4)
    sb_mag = 22.5 - 2.5 * np.log10(prof_norm * total_flux / pixelscale)
    
    return sb_mag

def VelMap2D(x, y, vmax = 200, rt = 1, ellip = 0.5, 
             theta = 0, x_0 = 0, y_0 = 0):
    
    # Produce tanh profile
    a, b = rt, (1 - ellip) * rt
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
    x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
    z = (abs(x_maj / a) ** 2 + abs(x_min / b)  ** 2) ** (1 / 2)
    profile = vmax * np.tanh(z) * ((x_maj / a) / z)
    
    return profile

def GradMap2D(x, y, a0 = 10, r_eff = 1, gred = -1, ellip = 0.5, 
              theta = 0, x_0 = 0, y_0 = 0):
    
    # Produce gradiant profile
    a, b = r_eff, (1 - ellip) * r_eff
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta
    x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta
    z = (abs(x_maj / a) ** 2 + abs(x_min / b)  ** 2) ** (1 / 2)
    profile = a0 + z * gred
    
    return profile

class Map2d(object):

    """
    Generate an x, y grid in a rectangular region, sampled with xsamp and
    ysamp spacings in the x and y directions, respectively.  Origin is at the
    upper-left corner to match numpy and matplotlib convention. astropy.io.fits
    also assumes UL origin by default.

    Parameters
    ----------
    xsamp, ysamp : float, float
        The sampling spacing in the x and y directions.
    nx, ny: int, int
        Number of samples in the x and y directions.
    """

    def __init__(self, inst):
        self.xsamp = inst.dpix
        self.ysamp = inst.dpix
        startx = -(inst.nx - 1) / 2.0 * self.xsamp
        stopx  = (inst.nx - 1) / 2.0 * self.xsamp
        starty = -(inst.ny - 1) / 2.0 * self.ysamp
        stopy  = (inst.ny - 1) / 2.0 * self.ysamp
        xvals  = np.linspace(startx, stopx, num = inst.nx)
        yvals  = np.linspace(starty, stopy, num = inst.ny)

        ones = np.ones((inst.ny, inst.nx))
        x = ones * xvals
        y = np.flipud(ones * yvals.reshape(int(inst.ny), 1))

        self.nx = inst.nx
        self.ny = inst.ny
        self.x  = x
        self.y  = y
        self.row = xvals
        # flip Y axis because we use Y increasing from bottom to top
        self.col = yvals[::-1]
    
    def sersic_map(self, mag = 12, r_eff = 2, n = 2.5, ellip = 0.5, theta = -50):
        self.mag   = mag
        self.reff  = r_eff / self.xsamp
        self.n     = n
        self.ellip = ellip
        self.theta = theta
        self.map = Sersic2D(self.x, self.y, mag = self.mag, 
                            r_eff = self.reff, n = self.n, 
                            ellip = self.ellip, theta = self.theta, 
                            pixelscale = self.xsamp * self.ysamp)
        
    def tanh_map(self, vmax = 200, rt = 2, ellip = 0.5, theta = -50):
        self.vmax  = vmax
        self.rt    = rt / self.xsamp
        self.ellip = ellip
        self.theta = theta
        self.map = VelMap2D(self.x, self.y, vmax = self.vmax, rt = self.rt, 
                            ellip = self.ellip, theta = self.theta)
    
    def gred_map(self, a0 = 10, r_eff = 1, gred = -1, ellip = 0.5, theta = 0):
        self.a0    = a0
        self.reff  = r_eff / self.xsamp
        self.gred  = gred
        self.ellip = ellip
        self.theta = theta
        self.map = GradMap2D(self.x, self.y, a0 = self.a0, r_eff = self.reff, 
                             gred = self.gred, ellip = self.ellip, theta = self.theta)
        
    def load_map(self, image):
        if np.ndim(image) == 2:
            self.map = resize(image, (self.nx, self.ny))
            
class StellarPopulationMap():
    
    def __init__(self, inst, sbright = None, logage = None, 
                 feh = None, vel = None, vdisp = None, ebv = None):
        self.nx    = inst.nx
        self.ny    = inst.ny
        self.dpix  = inst.dpix
        self.fov_x = inst.fov_x
        self.fov_y = inst.fov_y
        
        self.sbright = sbright.map
        self.logage  = logage.map
        self.feh     = feh.map
        self.vel     = vel.map
        self.vdisp   = vdisp.map
        self.ebv     = ebv.map
        
        self.mag = self.sbright + 2.5 * np.log10(self.dpix * self.dpix)
        self.age = 10 ** self.logage / 1e9
        
        self.vdisp[self.vdisp < 10] = 10
        self.ebv[self.ebv < 0]      = 0
        
class IonizedGasMap():
    
    def __init__(self, inst, halpha = None, zh = None, vel = None, vdisp = None, ebv = None):
        
        self.nx    = inst.nx
        self.ny    = inst.ny
        self.dpix  = inst.dpix
        self.fov_x = inst.fov_x
        self.fov_y = inst.fov_y
        
        self.halpha = halpha.map
        self.zh     = zh.map
        self.vel    = vel.map
        self.vdisp  = vdisp.map
        self.ebv    = ebv.map
        
        self.vdisp[self.vdisp < 10] = 10
        self.ebv[self.ebv < 0]      = 0