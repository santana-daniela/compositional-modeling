from __future__ import division
from math import exp, pow
import numpy as np
import pandas as pd


def flash_mixture(mix, pressure, temperature):
    flashed = False
    pr = pressure / mix.fluids['pc']
    tr = temperature / mix.fluids['tc']
    kcs = (1/pr)*np.exp(5.3737*(1+mix.fluids['w'])*(1-(1/tr)))
    zcs = np.array(mix.fluids['z'])
    iters = 1000
    while not flashed:
        if func_v_mix(0, kcs, zcs) < 0:
            #liquid
            mix.vapor = 0
            mix.fluids['x'] = mix.fluids['z']
            mix.fluids['y'] = mix.fluids['z'] * kcs
        elif func_v_mix(1, kcs, zcs) > 0:
            #vapor
            mix.vapor = 1
            mix.fluids['y'] = mix.fluids['z']
            mix.fluids['x'] = mix.fluids['z'] / kcs
        else:
            identify_phases(mix, kcs)
        mix.normalize_phases()
        mix.update_z_factors(pressure, temperature)
        mix.update_fugacities()
        fl = mix.fluids['fl'] * mix.fluids['x'] * pressure
        fv = mix.fluids['fv'] * mix.fluids['y'] * pressure
        rfs = fl / fv
        error = rfs.max()
        if error < pow(10, -6) or iters < 0:
            flashed = True
        else:
            iters -= 1
            kcs = kcs * rfs


def func_v_mix(v, kcs, zcs):
    fv = (kcs-1)*zcs/((kcs-1)*v+1)
    fvm = fv.sum()
    return fvm


def der_func_v_mix(v, kcs, zcs):
    dfv = -(np.power(kcs-1, 2)*zcs)/np.power((kcs-1)*v+1, 2)
    dfvm = dfv.sum()
    return dfvm


def identify_phases(mixture, kcs):
    zcs = mixture.fluids['z']
    v = 0.5
    error = 1
    tol = pow(10, -6)
    iters = 3000
    while error > tol and iters > 0:
        fv = func_v_mix(v, kcs, zcs)
        dfv = der_func_v_mix(v, kcs, zcs)
        v2 = v - fv/dfv
        if v2 <= 0 or v2 >= 1:
            v2 = (v+1)/2
            #v2 = v2
        error = abs((v2/v)-1) #ilkay's code: abs((v2-v)/v)
        iters -= 1
        v = v2
    xf = zcs / ((kcs-1)*v+1)
    yf = kcs * zcs / ((kcs-1)*v+1)
    mixture.vapor = v
    mixture.fluids['x'] = xf
    mixture.fluids['y'] = yf
