import numpy as np


# ----------------
# Reddening Module

def Calzetti_Law(wave, Rv = 4.05):
    """
    Dust Extinction Curve by Calzetti et al. (2000)
    """
    wave_number = 1./(wave * 1e-4)
    reddening_curve = np.zeros(len(wave))
    
    idx = np.logical_and(wave >= 1200, wave <= 6300)
    reddening_curve[idx] = 2.659 * ( -2.156 + 1.509 * wave_number[idx] - 0.198 * \
                    (wave_number[idx] ** 2)) + 0.011 * (wave_number[idx] **3 ) + Rv
                                    
    idx = np.logical_and(wave >= 6300, wave <= 22000)
    reddening_curve[idx] = 2.659 * ( -1.857 + 1.040 * wave_number[idx]) + Rv
    return reddening_curve

def reddening(wave, flux, ebv = 0, law = 'calzetti', Rv = 4.05):
    """
    Redden an input spectra through a given reddening curve.
    """
    curve = Calzetti_Law(wave, Rv = Rv)
    fluxNew = flux / (10. ** (0.4 * ebv * curve))
    return fluxNew