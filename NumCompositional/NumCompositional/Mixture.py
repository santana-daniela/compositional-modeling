from __future__ import division
from math import pow, sqrt, log, exp
import helper
from helper import R
import numpy as np
import pandas as pd


class Mixture:
    def __init__(self, molecules, zFracs):
        comps = len(molecules)
        self.comps = comps
        self.vapor = 0
        cols = ['molecule', 'z', 'x', 'y', 'mw', 'pc', 'tc', 'w', 'oa', 'ob', 'k', 'pa', 'pb', 'fl', 'fv']
        self.fluids = pd.DataFrame([[0 for _ in range(len(cols))] for _ in range(comps)], columns=cols)
        self.fluids['molecule'] = molecules
        self.fluids['z'] = zFracs
        self.pengs = {'L': [0, 0], 'V': [0, 0]}
        self.cap_pengs = {'L': [0, 0], 'V': [0, 0]}
        self.pressure = 0
        self.temperature = 0
        self.zFactors = {'L': 0, 'V': 0}
        helper.setup_fluids(self.fluids)
        self.fluids['k'] = 0.37464 + 1.54226 * self.fluids['w'] - 0.26992 * self.fluids['w'].pow(2)
        self.fluids['pb'] = self.fluids['ob'] * R * self.fluids['tc'] / self.fluids['pc']
        self.bic = np.matrix([np.zeros(comps) for _ in range(comps)])
        # hardcoded
        self.bic[0, 1] = 0.01475
        self.bic[1, 0] = self.bic[0, 1]
        self.bic[0, 2] = 0.04437
        self.bic[2, 0] = self.bic[0, 2]
        self.bic[1, 2] = 0.00845
        self.bic[2, 1] = self.bic[1, 2]

    def normalize_phases(self):
        self.fluids['x'] = self.fluids['x']/self.fluids['x'].sum()
        self.fluids['y'] = self.fluids['y'] / self.fluids['y'].sum()

    def update_mix_pengs(self):
        for l, x in [['L', 'x'], ['V', 'y']]:
            pb_mix = (self.fluids['pb']*self.fluids[x]).sum()
            sq_ax = np.sqrt(self.fluids['pa'])*self.fluids[x]
            pa_mix = (np.multiply(np.outer(sq_ax, sq_ax), (1-self.bic))).sum()
            self.pengs[l] = [pa_mix, pb_mix]
            cap_a = pa_mix * self.pressure / pow(R*self.temperature, 2)
            cap_b = pb_mix * self.pressure / R / self.temperature
            self.cap_pengs[l] = [cap_a, cap_b]

    def calc_z_factor(self, pengs):
        pa, pb = pengs
        a1 = -(1-pb)
        a2 = pa-2*pb-3*pow(pb, 2)
        a3 = -(pa*pb-pow(pb, 2)-pow(pb, 3))
        roots = helper.eos_roots(a1, a2, a3)
        roots.sort()
        return roots

    def update_fluids(self):
        pas = (1+self.fluids['k']*(1 - (self.temperature/self.fluids['tc']).pow(0.5))).pow(2)
        self.fluids['pa'] = self.fluids['oa'] * pow(R, 2) * self.fluids['tc'].pow(2) / self.fluids['pc'] * pas

    def update_z_factors(self, pressure, temperature):
        self.pressure = pressure
        self.temperature = temperature
        self.update_fluids()
        self.update_mix_pengs()
        lzf = self.calc_z_factor(self.cap_pengs['L'])[0]
        vzfs = self.calc_z_factor(self.cap_pengs['V'])
        vzf = vzfs[len(vzfs)-1]
        self.zFactors['L'] = lzf
        self.zFactors['V'] = vzf

    def update_fugacities(self):
        sq_pa = np.sqrt(self.fluids['pa'])
        for l, x in [['L', 'x'], ['V', 'y']]:
            amn_sum = np.array((np.multiply(np.outer(sq_pa, sq_pa*self.fluids[x]), (1-self.bic))).sum(axis=1).T)[0]
            z = self.zFactors[l]
            amix, bmix = self.pengs[l]
            cap_a, cap_b = self.cap_pengs[l]
            frp1 = self.fluids['pb']/bmix*(z-1)-log(z-cap_b)
            frp2 = -1/sqrt(8)*(cap_a/cap_b)*(2*amn_sum/amix-self.fluids['pb']/bmix)
            frp3 = log((z+(sqrt(2)+1)*cap_b)/(z-(sqrt(2)-1)*cap_b))
            fr = np.exp(frp1 + frp2 * frp3)
            fug = 'f' + l.lower()
            self.fluids[fug] = fr





