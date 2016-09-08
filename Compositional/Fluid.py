from __future__ import division
from math import pow, sqrt
from helper import R


class Fluid:
    def __init__(self, molecule):
        self.molecule = molecule
        self.c_pressure = 0
        self.c_temperature = 0
        self.accentric = 0
        self.weight = 0
        self.omega_a = 0
        self.omega_b = 0
        self.peng_a = 0
        self.peng_b = 0
        self.k = 0

    def r_pressure(self, pressure):
        return pressure/self.c_pressure

    def r_temperature(self, temperature):
        return temperature/self.c_temperature

    def k_acc(self):
        if self.k == 0:
            k = 0.37464 + 1.54226*self.accentric\
                - 0.26992*pow(self.accentric, 2)
            self.k = k
        else:
            k = self.k
        return k

    def prepare_eos(self, temperature):
        k = self.k_acc()
        oa = self.omega_a
        ob = self.omega_b
        tc = self.c_temperature
        pc = self.c_pressure
        tr = self.r_temperature(temperature)
        if self.peng_b == 0:
            self.peng_b = ob*R*tc/pc
        self.peng_a = (oa*R*R*tc*tc/pc)*pow(1+k*(1-sqrt(tr)), 2)
