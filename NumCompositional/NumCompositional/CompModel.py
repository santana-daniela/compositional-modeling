from __future__ import division
import Flash
import pandas as pd
from matplotlib import pyplot as plt
from Mixture import Mixture
from helper import normalize_list


psi = 1500
rankine = 300+460
molecules = ['C1', 'C4', 'C10']
zFracs = [0.5, 0.42, 0.08]


def create_composition_data(points):
    ternary_cols = ['M', 'z', 'x', 'y']
    rows = len(points)
    comp_dfs = []
    for comp in range(len(molecules)):
        comp_df = pd.DataFrame([[0 for _ in range(rows)] for _ in range(len(ternary_cols))], ternary_cols)
        comp_dfs.append(comp_df)
    count = 0
    for i, j, k in points:
        zFracs = [i+0.1,j+0.1,k+0.1]
        normalize_list(zFracs)
        mix = Mixture(molecules, zFracs)
        Flash.flash_mixture(mix, psi, rankine)
        for m in range(len(molecules)):
            comp_dfs[m][count] = mix.fluids[['molecule', 'z', 'x', 'y']].ix[m]
        count += 1
    return comp_dfs
