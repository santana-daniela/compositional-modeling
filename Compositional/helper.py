from __future__ import division
import pandas as pd
import numpy as np
import math

R = 10.73159
properties_df = pd.read_csv("molecular_properties.csv", index_col=0)


def normalize(comps):
    summed = 0
    for c in comps:
        summed += c
    for i in range(len(comps)):
        comps[i] = max(comps[i]/summed, 0.0001)


def setup_fluid(fluid):
    properties = properties_df.ix[fluid.molecule]
    fluid.weight = (properties['weight']*0.00220462)
    fluid.c_pressure = (properties['pressure']*14.6959)
    fluid.c_temperature = (properties['temperature']*1.8)
    fluid.accentric = properties['accentric']
    fluid.omega_a = properties['omega_a']
    fluid.omega_b = properties['omega_b']


def isnan(check):
    return np.isnan(check)


def eos_roots(c2, c3, c4):
    roots = np.roots([1, c2, c3, c4])
    pr = roots[~np.iscomplex(roots)]
    reals = pr[~np.isnan(pr)]

    rr = []
    for real in reals:
        if real > 0:
            rr.append(real.real)
    return rr
