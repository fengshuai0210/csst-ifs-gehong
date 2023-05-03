import sys
sys.path.insert(0, '..')
import unittest
from gehong import map2d as m
from gehong import config as c
import numpy as np
import matplotlib.pyplot as plt

import warnings

class test_map2d(unittest.TestCase):

    """
    Test modules in map2d.py
    """

    def test_sbmap(self):

        inst = c.config()

        print("------------------------- Sersic Model -------------------------")
        print(" ")
        print("Pars: Total Magnitude")
        sbmap0 = m.Map2d(inst)
        sbmap0.sersic_map(mag = 12, r_eff = 2, n = 2.5, ellip = 0.5, theta = -50)
        sbmap1 = m.Map2d(inst)
        sbmap1.sersic_map(mag = 14, r_eff = 2, n = 2.5, ellip = 0.5, theta = -50)
        if np.sum(sbmap1.map) != np.sum(sbmap0.map):
            print("No Problem!")

        print("Pars: Effective Radius")
        sbmap1 = m.Map2d(inst)
        sbmap1.sersic_map(mag = 12, r_eff = 4, n = 2.5, ellip = 0.5, theta = -50)
        if np.sum(sbmap1.map) != np.sum(sbmap0.map):
            print("No Problem!")

        print("Pars: Sersic Index")
        sbmap1 = m.Map2d(inst)
        sbmap1.sersic_map(mag = 12, r_eff = 2, n = 5, ellip = 0.5, theta = -50)
        if np.sum(sbmap1.map) != np.sum(sbmap0.map):
            print("No Problem!")

        print("Pars: Ellipse")
        sbmap1 = m.Map2d(inst)
        sbmap1.sersic_map(mag = 12, r_eff = 2, n = 2.5, ellip = 0.2, theta = -50)
        if np.sum(sbmap1.map) != np.sum(sbmap0.map):
            print("No Problem!")

        print("Pars: Position Angle")
        sbmap1 = m.Map2d(inst)
        sbmap1.sersic_map(mag = 12, r_eff = 2, n = 2.5, ellip = 0.5, theta = 0)
        if np.sum(sbmap1.map) != np.sum(sbmap0.map):
            print("No Problem!")

        plt.figure(figsize=(6,5))
        plt.imshow(sbmap0.map, origin='lower', cmap='gray_r')
        plt.xticks([])
        plt.yticks([])
        plt.title('Surface Brightness Modelling')
        plt.colorbar(label='mag')
        plt.tight_layout()
        plt.savefig('./image/test_map2d_SBmap.png')

    def test_velmap(self):

        inst = c.config()

        print("------------------------ Velocity Model ------------------------")
        print(" ")
        print("Pars: Rotation Velocity")
        velmap0 = m.Map2d(inst)
        velmap0.tanh_map(vmax = 200, rt = 2, ellip = 0.5, theta = -50)
        velmap1 = m.Map2d(inst)
        velmap1.tanh_map(vmax = 100, rt = 2, ellip = 0.5, theta = -50)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Turnover Radius")
        velmap1 = m.Map2d(inst)
        velmap1.tanh_map(vmax = 200, rt = 3, ellip = 0.5, theta = -50)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Inclination Angle")
        velmap1 = m.Map2d(inst)
        velmap1.tanh_map(vmax = 200, rt = 2, ellip = 0.3, theta = -50)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Position Angle")
        velmap1 = m.Map2d(inst)
        velmap1.tanh_map(vmax = 200, rt = 2, ellip = 0.5, theta = 100)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        plt.figure(figsize=(6,5))
        plt.imshow(velmap0.map, origin='lower', cmap='RdBu_r')
        plt.xticks([])
        plt.yticks([])
        plt.title('Velocity Map Modelling')
        plt.colorbar(label='km/s')
        plt.tight_layout()
        plt.savefig('./image/test_map2d_velmap.png')

    def test_gremap(self):

        inst = c.config()

        print("------------------------ Gredient Model ------------------------")
        print(" ")
        print("Pars: Central Intensity")
        velmap0 = m.Map2d(inst)
        velmap0.gred_map(a0 = 10, r_eff = 1, gred = -1, ellip = 0.5, theta = 0)
        velmap1 = m.Map2d(inst)
        velmap1.gred_map(a0 = 20, r_eff = 1, gred = -1, ellip = 0.5, theta = 0)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Effective Radius")
        velmap1 = m.Map2d(inst)
        velmap1.gred_map(a0 = 10, r_eff = 2, gred = -1, ellip = 0.5, theta = 0)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")
        
        print("Pars: Grediant")
        velmap1 = m.Map2d(inst)
        velmap1.gred_map(a0 = 10, r_eff = 1, gred = -2, ellip = 0.5, theta = 0)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Inclination Angle")
        velmap1 = m.Map2d(inst)
        velmap1.gred_map(a0 = 10, r_eff = 1, gred = -1, ellip = 0.3, theta = 0)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        print("Pars: Position Angle")
        velmap1 = m.Map2d(inst)
        velmap1.gred_map(a0 = 10, r_eff = 1, gred = -1, ellip = 0.5, theta = 100)
        if np.sum(velmap1.map) != np.sum(velmap0.map):
            print("No Problem!")

        plt.figure(figsize=(6,5))
        plt.imshow(velmap0.map, origin='lower', cmap='RdBu_r')
        plt.xticks([])
        plt.yticks([])
        plt.title('Gredient Map Modelling')
        plt.colorbar(label='')
        plt.tight_layout()
        plt.savefig('./image/test_map2d_gredmap.png')

    def test_input_image(self):

        inst = c.config()

        print("------------------------ Feed Image------------------------")
        print(" ")
        velmap0 = m.Map2d(inst)
        velmap0.gred_map(a0 = 10, r_eff = 1, gred = -1, ellip = 0.5, theta = 0)
        velmap1 = m.Map2d(inst)
        velmap1.load_map(velmap0.map)
        if np.sum(velmap1.map) == np.sum(velmap0.map):
            print("No Problem!")

    def test_stellar_population_map(self):

        inst = c.config()

        sbmap = m.Map2d(inst)
        sbmap.sersic_map()

        velmap = m.Map2d(inst)
        velmap.tanh_map()

        vdispmap = m.Map2d(inst)
        vdispmap.gred_map()

        agemap = m.Map2d(inst)
        agemap.gred_map(a0 = 9.5, gred = -1.2)

        fehmap = m.Map2d(inst)
        fehmap.gred_map(a0 = -0.2, gred = -0.1)

        ebvmap = m.Map2d(inst)
        ebvmap.gred_map(a0 = 0.2, gred = -0.1)

        stellarcontinuum = m.StellarPopulationMap(inst, sbright = sbmap, logage = agemap, feh = fehmap, 
                                                  vel = velmap, vdisp = vdispmap, ebv = ebvmap)
        
        if np.sum(stellarcontinuum.sbright) == np.sum(sbmap.map):
            print("SBmap No Problem!")
        if np.sum(stellarcontinuum.age) == np.sum(agemap.map):
            print("Agemap No Problem!")
        if np.sum(stellarcontinuum.feh) == np.sum(fehmap.map):
            print("FeHmap No Problem!")
        if np.sum(stellarcontinuum.vel) == np.sum(velmap.map):
            print("Velmap No Problem!")
        if np.sum(stellarcontinuum.vdisp) == np.sum(vdispmap.map):
            print("VDISPmap No Problem!")
        if np.sum(stellarcontinuum.ebv) == np.sum(ebvmap.map):
            print("SBmap No Problem!")

    def test_ionized_gas_map(self):

        inst = c.config()

        halphamap = m.Map2d(inst)
        halphamap.sersic_map()

        velmap = m.Map2d(inst)
        velmap.tanh_map()

        vdispmap = m.Map2d(inst)
        vdispmap.gred_map()

        logohmap = m.Map2d(inst)
        logohmap.gred_map(a0 = 9.5, gred = -1.2)

        ebvmap = m.Map2d(inst)
        ebvmap.gred_map(a0 = 0.2, gred = -0.1)

        gas = m.IonizedGasMap(inst, halpha = halphamap, zh = logohmap, 
                              vel = velmap, vdisp = vdispmap, ebv = ebvmap)
        
        if np.sum(gas.halpha) == np.sum(halphamap.map):
            print("Halphamap No Problem!")
        if np.sum(gas.zh) == np.sum(logohmap.map):
            print("logOHmap No Problem!")
        if np.sum(gas.vel) == np.sum(velmap.map):
            print("Velmap No Problem!")
        if np.sum(gas.vdisp) == np.sum(vdispmap.map):
            print("VDISPmap No Problem!")
        if np.sum(gas.ebv) == np.sum(ebvmap.map):
            print("SBmap No Problem!")

if __name__=='__main__':
    unittest.main()