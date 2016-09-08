from __future__ import division
from math import exp, pow


def flash_mixture(mixture, pressure, temperature):
    flashed = False
    kcs = []
    zcs = mixture.zFracs
    for fluid in mixture.fluids:
        pr = fluid.r_pressure(pressure)
        tr = fluid.r_temperature(temperature)
        w = fluid.accentric
        kc = (1/pr)*exp(5.3737*(1+w)*(1-(1/tr)))
        kcs.append(kc)
    iters = 1000
    while not flashed:
        if func_v_mix(0, kcs, zcs) < 0:
            #liquid
            for i in range(len(mixture.fluids)):
                mixture.xFracs[i] = zcs[i]
                mixture.yFracs[i] = zcs[i] * kcs[i]
        elif func_v_mix(1, kcs, zcs) > 0:
            #vapor
            for i in range(len(mixture.fluids)):
                mixture.xFracs[i] = zcs[i]
                #ilkay's code does not match with this
                mixture.yFracs[i] = zcs[i] / kcs[i]
        else:
            identify_mix(mixture, kcs)
        mixture.normalize_phases()
        mixture.update_z_factors(pressure, temperature)
        mixture.update_fugacities()
        rfs = []
        error = 0
        for i in range(len(mixture.fluids)):
            fl = mixture.fugacities['L'][i]*mixture.xFracs[i]*pressure
            fv = mixture.fugacities['V'][i]*mixture.yFracs[i]*pressure
            r = fl / fv
            rfs.append(r)
            error = max(error, abs(r-1))
        if error < pow(10, -6) or iters < 0:
            print 'final kcs:', kcs
            print 'error:', error
            flashed = True
        else:
            iters -= 1
            print iters
            for i in range(len(mixture.fluids)):
                kcs[i] = kcs[i]*rfs[i]


def func_v_mix(v, kcs, zcs):
    fvm = 0
    for i in range(len(kcs)):
        fvm += func_v(v, kcs[i], zcs[i])
    return fvm


def func_v(v, kc, zc):
    fv = (kc-1)*zc/((kc-1)*v+1)
    return fv


def der_func_v_mix(v, kcs, zcs):
    dfvm = 0
    for i in range(len(kcs)):
        dfvm += der_func_v(v, kcs[i], zcs[i])
    #for i in range(len(kcs)):
    #    dfvm = -(dfvm + (pow(kcs[i]-1, 2)*zcs[i])/pow((kcs[i]-1)*v+1, 2))
    return dfvm


def der_func_v(v, kc, zc):
    dfv = -(pow((kc-1), 2)*zc)/pow((kc-1)*v+1, 2)
    return dfv


def identify_mix(mixture, kcs):
    zcs = mixture.zFracs
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
    for i in range(len(mixture.fluids)):
        x_f = zcs[i] / ((kcs[i]-1)*v+1)
        y_f = kcs[i] * zcs[i] / ((kcs[i]-1)*v+1)
        mixture.xFracs[i] = max(x_f, 0.00001)
        mixture.yFracs[i] = max(y_f, 0.00001)
