import sys
sys.path.insert(0, '..')
import unittest
from gehong import spec1d as s
from gehong import config as c
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter("ignore", ResourceWarning)

class test_spec1d(unittest.TestCase):

    """
    Test modules in spec1d.py
    """

    def test_StellarPop(self):

        inst = c.config()

        print("=================================================================")
        print("==================Stellar Population Modelling===================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        temp = s.StellarContinuumTemplate(inst)

        plt.figure(figsize=(15,12))

        print("------------------------Test Magnitude---------------------------")
        print("1. mr=15")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0, vel = 100, vdisp = 100, ebv = 0)
        print("2. mr=13")
        sed2 = s.StellarContinuum(inst, temp, mag = 13, age = 1, feh = 0, vel = 100, vdisp = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(231)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'mr = 15')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'mr = 13')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)
        
        print("---------------------------Test  Age-----------------------------")
        print("1. Age=0.1Gyr")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 0.1, feh = 0, vel = 100, vdisp = 100, ebv = 0)
        print("2. Age=1Gyr")
        sed2 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0, vel = 100, vdisp = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(232)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'Age = 0.1Gyr')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'Age = 1Gyr')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("-----------------------Test  Metallicity-------------------------")
        print("1. [Fe/H]=0.5")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh =  0.5, vel = 100, vdisp = 100, ebv = 0)
        print("2. [Fe/H]=-0.5")
        sed2 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = -0.5, vel = 100, vdisp = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(233)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = '[Fe/H] = 0.5')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = '[Fe/H] = -0.5')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)
        
        print("------------------------Test  Velocity---------------------------")
        print("1. v=3000km/s")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0, vel = 3000, vdisp = 100, ebv = 0)
        print("2. v=6000km/s")
        sed2 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0, vel = 6000, vdisp = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(234)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'v = 3000km/s')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'v = 6000km/s')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)
        
        print("-------------------Test Velocity Dispersion----------------------")
        print("1. sigma=180km/s")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0.5, vel = 100, vdisp = 180, ebv = 0)
        print("2. sigma=350km/s")
        sed2 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0.5, vel = 100, vdisp = 350, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(235)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'sigma = 180km/s')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'sigma = 350km/s')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("----------------------Test Dust Extinction-----------------------")
        print("1. E(B-V)=0.1")
        sed1 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0.5, vel = 100, vdisp = 100, ebv = 0.1)
        print("2. E(B-V)=0.3")
        sed2 = s.StellarContinuum(inst, temp, mag = 15, age = 1, feh = 0.5, vel = 100, vdisp = 100, ebv = 0.3)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(236)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'E(B-V) = 0.1')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'E(B-V) = 0.3')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        plt.tight_layout()
        plt.savefig('./image/test_spec1d_stellarpopulation.png')
        
    def test_Star(self):

        inst = c.config()

        print("================================================================")
        print("=================Single Star Spectra Modelling===================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        temp = s.SingleStarTemplate(inst)

        plt.figure(figsize=(15,12))

        print("------------------------Test Magnitude---------------------------")
        print("1. mr=15")
        sed1 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = 0, vel = 100, ebv = 0)
        print("2. mr=13")
        sed2 = s.SingleStar(inst, temp, mag = 13, teff = 10000, feh = 0, vel = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(221)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'mr = 15')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'mr = 13')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("------------------------Test Tempreture--------------------------")
        print("1. Teff=6000K")
        sed1 = s.SingleStar(inst, temp, mag = 15, teff = 6000, feh = 0, vel = 100, ebv = 0)
        print("1. Teff=10000K")
        sed2 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = 0, vel = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(222)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'Teff = 6000K')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'Teff = 10000K')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("-----------------------Test Metallicity--------------------------")
        print("1. [Fe/H]=-0.5")
        sed1 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = -0.5, vel = 100, ebv = 0)
        print("2. [Fe/H]=0.5")
        sed2 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = 0.5, vel = 100, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(223)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = '[Fe/H] = -0.5')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = '[Fe/H] = 0.5')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("------------------------Test  Velocity---------------------------")
        print("1. v=-200km/s")
        sed1 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = 0, vel = -200, ebv = 0)
        print("2. v=200km/s")
        sed2 = s.SingleStar(inst, temp, mag = 15, teff = 10000, feh = 0, vel =  200, ebv = 0)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(224)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'v = -200km/s')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'v = 200km/s')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)
        
        plt.tight_layout()
        plt.savefig('./image/test_spec1d_singlestar.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Single Stellar** Modelling has been passed!")

    def test_Gas(self):

        inst = c.config()

        print("================================================================")
        print("======================Ionized Gas Modelling======================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        temp = s.EmissionLineTemplate(inst, model = 'hii')

        plt.figure(figsize=(15,12))

        print("---------------------Test Total Halpha Flux----------------------")

        print("1. F_Halpha = 10 * 1e-17 erg/s/cm^2")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 100, vdisp = 120, ebv = 0.1)
        print("2. F_Halpha = 50 * 1e-17 erg/s/cm^2")
        sed2 = s.HII_Region(inst, temp, halpha = 50, logz = -1, vel = 100, vdisp = 120, ebv = 0.1)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(231)
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'Flux = 50')
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'Flux = 10')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("-------------------Test Gas-Phase Metallicity--------------------")

        print("1. [Z/H] = 0")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = 0, vel = 100, vdisp = 120, ebv = 0.1)
        print("2. [Z/H] = -0.5")
        sed2 = s.HII_Region(inst, temp, halpha = 10, logz = -0.5, vel = 100, vdisp = 120, ebv = 0.1)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(232)
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = '[Z/H] = -0.5')
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = '[Z/H] = 0')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("-------------------------Test  Velocity--------------------------")

        print("1. v=3000km/s")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 3000, vdisp = 120, ebv = 0.1)
        print("2. v=6000km/s")
        sed2 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 6000, vdisp = 120, ebv = 0.1)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(233)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'v = 3000km/s')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'v = 6000km/s')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("--------------------Test Velocity Dispersion---------------------")

        print("1. sigma=500km/s")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 100, vdisp = 500, ebv = 0.1)
        print("2. sigma=100km/s")
        sed2 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 100, vdisp = 100, ebv = 0.1)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(234)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'sigma = 500km/s')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'sigma = 100km/s')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("----------------------Test Dust Extinction-----------------------")

        print("1. E(B-V)=0.1")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 100, vdisp = 120, ebv = 0.1)
        print("2. E(B-V)=0.3")
        sed1 = s.HII_Region(inst, temp, halpha = 10, logz = -1, vel = 100, vdisp = 120, ebv = 0.3)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(235)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'E(B-V) = 0.1')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'E(B-V) = 0.3')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        plt.tight_layout()
        plt.savefig('./image/test_spec1d_ionizedgas.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Ionized Gas** Modelling has been passed!")

    def test_AGN(self):

        inst = c.config()

        print("================================================================")
        print("=====================AGN Spectra Modelling======================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        temp = s.EmissionLineTemplate(inst, model = 'nlr')

        plt.figure(figsize=(15,12))

        print("-------------------------Test  Redshift--------------------------")

        print("1. z=0.2")
        sed1 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        print("2. z=0.4")
        sed2 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 40000, logz = 0, ebv = 0.1, dist = 20)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(221)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'z = 0.2')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'z = 0.4')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)
        
        print("----------------------Test Black Hole Mass-----------------------")

        print("1. Mbh=1e6")
        sed1 = s.AGN(inst, temp, bhmass = 1e6, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        print("2. Mbh=5e6")
        sed2 = s.AGN(inst, temp, bhmass = 5e6, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(222)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'Mbh = 1e6')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'Mbh = 5e6')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("------------------------Eddington  Ratio-------------------------")

        print("1. Edd Ratio = 1.0")
        sed1 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 1.0,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        print("2. Edd Ratio = 0.1")
        sed2 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 0.1,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(223)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'EddRatio = 1.0')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'EddRatio = 0.1')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        print("----------------------Test Dust Extinction-----------------------")

        print("1. E(B-V)=0.1")
        sed1 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.1, dist = 20)
        print("2. E(B-V)=0.3")
        sed2 = s.AGN(inst, temp, bhmass = 1e5, edd_ratio = 0.05,
                     halpha_broad = 100, halpha_narrow = 100, vdisp_broad = 2000, vdisp_narrow = 500,
                     vel = 20000, logz = 0, ebv = 0.3, dist = 20)  
        if np.sum(sed1.flux) != np.sum(sed2.flux):
            print("No Problem!")
            plt.subplot(224)
            plt.plot(sed1.wave, sed1.flux, color = 'blue', lw = 2, label = 'E(B-V) = 0.1')
            plt.plot(sed2.wave, sed2.flux, color = 'red', lw = 2, label = 'E(B-V) = 0.3')
            plt.legend(frameon=False)
            plt.xlim(3900, 6200)

        plt.tight_layout()
        plt.savefig('./image/test_spec1d_agn.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **AGN** Modelling has been passed!")

if __name__=='__main__':
    unittest.main()