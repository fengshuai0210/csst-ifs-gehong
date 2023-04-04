import unittest
import temp as tem
import sed as sed
import numpy as np
import matplotlib.pyplot as plt

import warnings

warnings.simplefilter('ignore', ResourceWarning)

class test_sed(unittest.TestCase):

    """
    Test modules in sed.py
    """

    def test_StellarPop(self):

        print("================================================================")
        print("==================Stellar Population Modelling===================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        wave = np.linspace(3500,10000,7000)
        ssp = tem.emiles(FWHM_inst = 3)

        print("------------------------Test Magnitude---------------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. mr=15")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = -0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='red', label='m_r=15')

        print("2. mr=13")
        s=sed.StarPop(ssp, wave, mag = 13, Age = 1, FeH = -0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label='m_r=13')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.Mag.png')

        print("---------------------------Test  Age-----------------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. Age=0.1Gyr")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 0.1, FeH = -0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='red', label='Age=0.1Gyr')

        print("2. Age=1Gyr")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = -0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label='Age=1Gyr')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.Age.png')

        print("-----------------------Test  Metallicity-------------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. [Fe/H]=0.5")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='red', label='[Fe/H]=0.5')

        print("2. [Fe/H]=-0.5")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = -0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label='[Fe/H]=-0.5')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.FeH.png')

        print("------------------------Test  Velocity---------------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. v=3000km/s")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='red', label='v=3000km/s')

        print("2. v=6000km/s")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 6000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label='v=6000km/s')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.Vel.png')

        print("-------------------Test Velocity Dispersion----------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. sigma=180km/s")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 180, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$\sigma=180km/s$')

        print("2. sigma=350km/s")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 350, Ebv = 0)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$\sigma=350km/s$')

        plt.xlim(3700,4300)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.Sig.png')

        print("----------------------Test Dust Extinction-----------------------")

        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. E(B-V)=0.1")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 180, Ebv = 0.1)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$E(B-V)=0.1$')

        print("2. E(B-V)=0.3")
        s=sed.StarPop(ssp, wave, mag = 15, Age = 1, FeH = 0.5, vel = 3000, vdisp = 180, Ebv = 0.3)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$E(B-V)=0.3$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.title('Stellar Population Modelling')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_StellarPop.EBV.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Stellar Continumm** Modelling has been passed!")

    def test_Star(self):

        print("================================================================")
        print("=================Single Star Spectra Modelling===================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        wave = np.linspace(3700,7000,7000)
        ssp = tem.stellarlib(FWHM_inst = 3)

        print("------------------------Test Magnitude---------------------------")
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. mr=15")
        s=sed.star(ssp, wave, 15, 6000, 0, 500)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$m_r=15$')

        print("2. mr=13")
        s=sed.star(ssp, wave, 13, 6000, 0, 500)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$m_r=13$')

        plt.xlim(3700,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Single Stellar Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_SingleStar.Mag.png')

        print("------------------------Test Tempreture--------------------------")
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. Teff=6000K")
        s=sed.star(ssp, wave, 15, 6000, 0, 500)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$T_{eff}=6000K$')

        print("1. Teff=10000K")
        s=sed.star(ssp, wave, 15, 10000, 0, 500)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$T_{eff}=10000K$')

        plt.xlim(3700,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Single Stellar Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_SingleStar.Teff.png')

        print("-----------------------Test Metallicity--------------------------")
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. [Fe/H]=-0.5")
        s=sed.star(ssp, wave, 15, 6000, -0.5, 500)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$[Fe/H]=-0.5$')

        print("2. [Fe/H]=0")
        s=sed.star(ssp, wave, 15, 6000, 0, 500)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$[Fe/H]=0$')

        plt.xlim(3700,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Single Stellar Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_SingleStar.FeH.png')

        print("------------------------Test  Velocity---------------------------")
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. v=-200km/s")
        s=sed.star(ssp, wave, 15, 6000, 0, -200)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$v=-200km/s$')

        print("2. v=200km/s")
        s=sed.star(ssp, wave, 15, 6000, 0, 200)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$v=200km/s$')

        plt.xlim(3800,4200)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False)
        plt.title('Single Stellar Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_SingleStar.Vel.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Single Stellar** Modelling has been passed!")

    def test_AGN(self):

        wave = np.linspace(3500,10000,7000)
        ssp=tem.AGNlib()
        s=sed.AGN(ssp, wave, vel=0.4*300000, Dist=20)

        print("================================================================")
        print("=====================AGN Spectra Modelling======================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        wave = np.linspace(3500,10000,7000)
        ssp=tem.AGNlib()

        print("-------------------------Test  Redshift--------------------------")
        
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. z=0.2")
        s=sed.AGN(ssp, wave, vel=0.2*300000)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$z=0.2$')

        print("2. z=0.4")
        s=sed.AGN(ssp, wave, vel=0.4*300000)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$z=0.4$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('AGN SED Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_AGN.Redshift.png')

        print("----------------------Test Black Hole Mass-----------------------")
        
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. Mbh=1e6")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=1e6)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$M_{bh}=1 \times 10^6$')

        print("2. Mbh=5e6")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=5e6)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$M_{bh}=5 \times 10^6$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('AGN SED Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_AGN.BHmass.png')

        print("------------------------Eddington  Ratio-------------------------")
        
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. Edd Ratio = 1.0")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=1e6, EddRatio=1)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$Edd=1.0$')

        print("2. Edd Ratio = 0.1")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=1e6, EddRatio=0.1)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$Edd=0.1$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('AGN SED Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_AGN.Edd.png')

        print("----------------------Test Dust Extinction-----------------------")
        
        # Plot Result
        plt.figure(figsize=(8,3))

        print("1. E(B-V)=0.1")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=1e6, Ebv=0.1)
        plt.plot(s.wave, s.flux, lw=2, color='red', label=r'$E(B-V)=0.1$')

        print("2. E(B-V)=0.3")
        s=sed.AGN(ssp, wave, vel=0.2*300000, Mbh=1e6, Ebv=0.3)
        plt.plot(s.wave, s.flux, lw=2, color='blue', label=r'$E(B-V)=0.3$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('AGN SED Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_AGN.EBV.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **AGN** Modelling has been passed!")

    def test_Gas(self):

        print("================================================================")
        print("======================Ionized Gas Modelling======================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        print("---------------------Test Total Halpha Flux----------------------")

        # Plot Result
        plt.figure(figsize=(8,3))
        wave = np.linspace(3500,10000,7000)

        print("1. F_Halpha = 1e-16 erg/s/cm^2")
        s=sed.Gas(wave, 10, 0, 500, 100)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$F_{H\alpha}=1 \times 10^{16}$')

        print("2. F_Halpha = 5e-16 erg/s/cm^2")
        s=sed.Gas(wave, 50, 0, 500, 100)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$F_{H\alpha}=5 \times 10^{-16}$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Ionized Gas Emission Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_Gas.HalphaFlux.png')

        print("-------------------Test Gas-Phase Metallicity--------------------")

        # Plot Result
        plt.figure(figsize=(8,3))
        wave = np.linspace(3500,10000,7000)

        print("1. [Z/H] > -0.2")
        s=sed.Gas(wave, 200, 0, 500, 100)
        plt.plot(s.wave, s.flux, lw=2, color='red', alpha=0.5, label=r'$[Z/H]=0$')

        print("2. [Z/H] = -0.5")
        s=sed.Gas(wave, 200, np.log10(0.006/0.02), 500, 100)
        plt.plot(s.wave, s.flux, lw=2, color='green', alpha=0.5, label=r'$[Z/H]=-0.5$')

        print("2. [Z/H] < -0.7")
        s=sed.Gas(wave, 200, -1, 500, 100)
        plt.plot(s.wave, s.flux, lw=2, color='blue', alpha=0.5, label=r'$[Z/H]=-1$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Ionized Gas Emission Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_Gas.FeH.png')

        print("-------------------------Test  Velocity--------------------------")

        # Plot Result
        plt.figure(figsize=(8,3))
        wave = np.linspace(3500,10000,7000)

        print("1. v=3000km/s")
        s1=sed.Gas(wave, 100, 1, 3000, 100)
        plt.plot(s1.wave, s1.flux, lw=2, color='red', alpha=0.5, label=r'$v=3000km/s$')

        print("2. v=6000km/s")
        s2=sed.Gas(wave, 100, 1, 6000, 100)
        plt.plot(s2.wave, s2.flux, lw=2, color='blue', alpha=0.5, label=r'$v=6000km/s$')

        plt.xlim(3500,7000)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Ionized Gas Emission Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_Gas.Vel.png')

        print("--------------------Test Velocity Dispersion---------------------")

        # Plot Result
        plt.figure(figsize=(8,3))
        wave = np.linspace(3500,10000,7000)

        print("1. sigma=500km/s")
        s1=sed.Gas(wave, 200, 1, 100, 500)
        plt.plot(s1.wave, s1.flux, lw=2, color='red', alpha=0.5, label=r'$\sigma=500km/s$')

        print("2. sigma=100km/s")
        s2=sed.Gas(wave, 200, 1, 100, 100)
        plt.plot(s2.wave, s2.flux, lw=2, color='blue', alpha=0.5, label=r'$\sigma=100km/s$')

        plt.xlim(4800,5200)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Ionized Gas Emission Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_Gas.Sigma.png')

        print("----------------------Test Dust Extinction-----------------------")

        # Plot Result
        plt.figure(figsize=(8,3))
        wave = np.linspace(3500,10000,7000)

        print("1. E(B-V)=0.1")
        s1=sed.Gas(wave, 200, 1, 100, 100, Ebv=0.1)
        plt.plot(s1.wave, s1.flux, lw=2, color='red', alpha=0.5, label=r'$E(B-V)=0.1$')

        print("2. E(B-V)=0.3")
        s2=sed.Gas(wave, 200, 1, 100, 100, Ebv=0.3)
        plt.plot(s2.wave, s2.flux, lw=2, color='blue', alpha=0.5, label=r'$E(B-V)=0.3$')

        plt.xlim(3700,6800)
        plt.xlabel(r'$Wavelength (\AA)$', fontsize=12)
        plt.ylabel(r'$Flux$', fontsize=12)
        plt.legend(frameon=False, loc='upper right')
        plt.title('Ionized Gas Emission Modelling')
        plt.tight_layout()
        plt.savefig('./Figure/test_sed_Gas.EBV.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Ionized Gas** Modelling has been passed!")

if __name__=='__main__':
    unittest.main()