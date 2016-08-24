from __future__ import division

from math import pow, sqrt, log, exp

import Fluid
import helper


class Mixture:
    def __init__(self):
        self.fluids = []
        self.zFracs = []
        self.pressure = 0
        self.temperature = 0
        self.xFracs = []
        self.yFracs = []
        self.zFactors = {'L': 0, 'V': 0}
        self.fugacities = {'L': [], 'V': []}
        self.pengs = {'L': [], 'V': []}
        self.mix_pengs = {'L': [], 'V': []}

    def initialize_fluids(self, molecules, zFracs):
        self.zFracs = zFracs
        for molecule in molecules:
            fluid = Fluid.Fluid(molecule)
            helper.setup_fluid(fluid)
            self.fluids.append(fluid)
        nr = range(len(self.fluids))
        self.fugacities = {'L': [0 for _ in nr], 'V': [0 for _ in nr]}
        self.xFracs = [0 for _ in nr]
        self.yFracs = [0 for _ in nr]
        self.pengs = {'L': [0 for x in nr], 'V': [0 for y in nr]}
        self.mix_pengs = {'L': [0 for x in nr], 'V': [0 for y in nr]}
        self.bic = [[0 for _ in nr] for _ in nr]
        #hardcoded
        self.bic[0][1] = 0.01475
        self.bic[1][0] = self.bic[0][1]
        self.bic[0][2] = 0.04437
        self.bic[2][0] = self.bic[0][2]
        self.bic[1][2] = 0.00845
        self.bic[2][1] = self.bic[1][2]

    def normalize_phases(self):
        helper.normalize(self.xFracs)
        helper.normalize(self.yFracs)

    def get_pengs(self, fracs):
        peng_a = []
        pb_mix = 0
        for i in range(len(self.fluids)):
            fluid_a = self.fluids[i].peng_a
            peng_a.append(fluid_a)
            fluid_b = self.fluids[i].peng_b
            pb_mix += fluid_b*fracs[i]
        pa_mix = 0
        for m in range(len(self.fluids)):
            for n in range(len(self.fluids)):
                pa = sqrt(peng_a[m]*peng_a[n])*fracs[m]*fracs[n]*(1-self.bic[m][n])
                pa_mix += pa
        pengs = [pa_mix, pb_mix]
        return pengs

    def cap_pengs(self, pengs, pressure, temperature):
        pa, pb = pengs
        pa_cap = pa*pressure/(pow(helper.R, 2) * pow(temperature, 2))
        pb_cap = pb*pressure/temperature/helper.R
        return [pa_cap, pb_cap]

    def calc_z_factor(self, pengs):
        pa, pb = pengs
        a1 = -(1-pb)
        a2 = pa-2*pb-3*pow(pb, 2)
        a3 = -(pa*pb-pow(pb, 2)-pow(pb, 3))
        roots = helper.eos_roots(a1, a2, a3)
        roots.sort()
        return roots

    def update_z_factors(self, pressure, temperature):
        self.pressure = pressure
        self.temperature = temperature
        for fluid in self.fluids:
            fluid.prepare_eos(temperature)
        l_pengs = self.get_pengs(self.xFracs)
        self.pengs['L'] = l_pengs
        l_caps = self.cap_pengs(l_pengs, pressure, temperature)
        self.mix_pengs['L'] = l_caps
        v_pengs = self.get_pengs(self.yFracs)
        self.pengs['V'] = v_pengs
        v_caps = self.cap_pengs(v_pengs, pressure, temperature)
        self.mix_pengs['V'] = v_caps
        lzf = self.calc_z_factor(l_caps)[0]
        vzfs = self.calc_z_factor(v_caps)
        vzf = vzfs[len(vzfs)-1]
        self.zFactors['L'] = lzf
        self.zFactors['V'] = vzf

    def update_fugacities(self):
        self.calc_fugacities(self.xFracs, 'L')
        self.calc_fugacities(self.yFracs, 'V')

    def calc_fugacities(self, fracs, var):
        num = len(self.fluids)
        temp = [0 for _ in range(num)]
        for m in range(num):
            for n in range(num):
                amn = sqrt(self.fluids[m].peng_a*self.fluids[n].peng_a)*(1-self.bic[m][n])
                temp[m] += fracs[n]*amn
        z = self.zFactors[var]
        amix, bmix = self.pengs[var]
        cap_a, cap_b = self.mix_pengs[var]
        for i in range(num):
            bc = self.fluids[i].peng_b
            frp1 = (bc/bmix*(z-1)-log(z-cap_b))
            frp2 = (cap_a/cap_b)*(2*temp[i]/amix-bc/bmix)
            frp3 = log((z+(sqrt(2)+1)*cap_b)/(z-(sqrt(2)-1)*cap_b))
            fr = exp(frp1 + -1/sqrt(8)*frp2*frp3)
            self.fugacities[var][i] = fr





