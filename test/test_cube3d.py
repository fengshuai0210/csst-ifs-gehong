import sys
sys.path.insert(0, '..')
import unittest
from gehong import spec1d as s
from gehong import map2d as m
from gehong import cube3d as b
from gehong import config as c
import numpy as np
import matplotlib.pyplot as plt

import warnings

class test_cube3d(unittest.TestCase):

    """
    Test modules in map2d.py
    """

    def test_cube(self):

        print("================================================================")
        print("=====================3D Spectrum Modelling======================")
        print(" ")
        print(">>> Start Test")
        print(" ")

        inst = c.config()
        gas_tem = s.EmissionLineTemplate(inst, model = 'hii')
        stellar_tem = s.StellarContinuumTemplate(inst)

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
        stellarcontinuum = m.StellarPopulationMap(inst, sbright = sbmap, logage = agemap, 
                                                  feh = fehmap, vel = velmap, 
                                                  vdisp = vdispmap, ebv = ebvmap)

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
        ionizedgas = m.IonizedGasMap(inst, halpha = agemap, zh = fehmap, 
                                     vel = velmap, vdisp = vdispmap, ebv = ebvmap)   

        u = b.Cube3D(inst, stellar_map = stellarcontinuum, gas_map = ionizedgas)

        u.make_cube(stellar_tem = stellar_tem, hii_tem = gas_tem)

        #u.savefits('result.fits')

        plt.figure(figsize=(16,4))
        plt.plot(u.wave,np.log10(u.flux[15,15,:]),color='blue',label=r'$v=%4.2f$'%(velmap.map[15,15]))
        plt.plot(u.wave,np.log10(u.flux[5,25,:]),color='red',label=r'$v=%4.2f$'%(velmap.map[5,25]))
        plt.xlim(3500,9500)
        plt.legend(frameon=False)
        plt.savefig('./image/test_cube3d_spec.png')

        print("-----------------------------------------------------------------")
        print(">>> Unit test for **Spec Cube** Modelling has been passed!")

if __name__=='__main__':
    unittest.main()