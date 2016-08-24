from __future__ import division
from Mixture import Mixture
from helper import normalize_list
import Flash
import matplotlib.pylab as lab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ternary
#pip install python-ternary


def project(p):
    SQRT3OVER2 = np.sqrt(3) / 2
    a, b, c = p
    x = a/2 + b
    y = SQRT3OVER2 * a
    return (x, y)

def reverse_project(x,y):
    SQRT3OVER2 = np.sqrt(3) / 2
    a = y / SQRT3OVER2
    b = x - a/2

def create_composition_data(points, scale):
    rows = len(points)
    #record interesting variables at each state
    mix_dfs = []
    #break each line of the triangle
    grid, vapor_grid = break_points(points, scale)
    #find edges of two phase region
    line_num = 0
    final_mixes = []
    for line in grid:
        pnt = 0
        mixes = solve_line(line)
        last_mix = None
        for mix in mixes:
            mix_dfs.append(mix.fluids)
            vapor_grid[line_num][pnt] = [mix.vapor, mix.fluids]
            pnt += 1
            if mix.vapor > 0 and mix.vapor < 1:
                last_mix = mix
        if last_mix is not None:
            final_mixes.append(last_mix)
        line_num += 1
    vapors = []
    c_data = []
    r, c, w = [0, 0, len(vapor_grid)]
    while len(vapors) < rows:
        vapors.append(vapor_grid[c][r][0])
        c_data.append(vapor_grid[c][r][1])
        c += 1
        if c == w:
            c = 0
            r += 1
            w -= 1
    #just a check to ensure everything was filled in
    return [vapors, c_data, final_mixes]


def break_points(points, scale):
    grid = [[] for _ in range(scale+1)]
    v_grid = [[] for _ in range(scale + 1)]
    for point in points:
        i, j, k = point
        grid[j].append(point)
        v_grid[j].append(-1)
    return [grid, v_grid]


def solve_point(point):
    i, j, k = point
    zcomp = normalize_list([j+0.01, i+0.01, k+0.01])
    print ('Solving zcomp:', zcomp)
    mix = Mixture(molecules, zcomp)
    Flash.flash_mixture(mix, psi, rankine)
    return mix


def solve_line(line):
    print ('Solving line:', line)
    tolerances = [50, 10]
    count = len(line)
    max_p = count - 1
    mixes = [None for _ in range(count)]
    #solve end points
    if count == 1:
        mix = solve_point(line[0])
        return [mix]
    mixes[0] = solve_point(line[0])
    mixes[max_p] = solve_point(line[max_p])
    simplify_line(line, mixes, tolerances)
    check = check_line(mixes)
    while not check:
        simplify_line(line, mixes, tolerances)
        check = check_line(mixes)
    return mixes


def check_line(mixes):
    for mix in mixes:
        if mix is None:
            return False
    return True


def simplify_line(line, line_mixes, tolerances):
    liq_tol, gas_tol = tolerances
    checker = -1
    last = len(line_mixes) - 1
    while checker < last:
        if line_mixes[checker] is None:
            checker += 1
            continue
        counter = checker + 1
        while checker < last and line_mixes[counter] is not None:
            counter += 1
            checker += 1
        if counter == last or checker == last:
            continue
        while counter <= last:
            if line_mixes[counter] is not None:
                break
            counter += 1
        difference = counter - checker
        v1 = line_mixes[checker].vapor
        v2 = line_mixes[counter].vapor
        generalize = False
        if v1 == 1 and v2 == 1 and difference <= gas_tol:
            generalize = True
        elif v1 == 0 and v2 == 0 and difference <= liq_tol:
            generalize = True
        if generalize:
            for i in range(checker + 1, counter):
                print ('Generalizing point:', line[i])
                line_mixes[i] = line_mixes[checker]
        else:
            mid = checker + int((counter-checker)/2)
            point = line[mid]
            mix = solve_point(point)
            line_mixes[mid] = mix
        checker = counter


def create_ternary(amt=30, heatmap=False, tie_lines=False, specific_trace=True, finish_curve=False):
    scale = 2
    points = len(list(ternary.helpers.simplex_iterator(scale)))
    while points <= amt:
        scale += 1
        points = len(list(ternary.helpers.simplex_iterator(scale)))
    scale -= 1
    points = list(ternary.helpers.simplex_iterator(scale))
    rows = len(points)
    print ('using number of points:', rows)
    print ('scale set to:', scale)
    cmap = plt.cm.get_cmap('summer')
    f, ax = plt.subplots(1, 1, figsize=(10, 8))
    figure, tax = ternary.figure(ax=ax, scale=scale)
    style = 'h'
    axes_colors = {'b': 'black', 'l': 'black', 'r': 'black'}
    mixed_data = []
    if heatmap:
        mixed_data = create_composition_data(points, scale)
        values = mixed_data[0]
        tax.heatmap(dict(zip(points, values)),
                        scale=scale,
                        cmap=cmap, vmax=1,
                        style=style, colorbar=True)

    tax.boundary(linewidth=2.0, axes_colors=axes_colors)

    tax.left_axis_label("C3", offset=0.16, color=axes_colors['l'])
    tax.right_axis_label("C1", offset=0.16, color=axes_colors['r'])
    tax.bottom_axis_label("C2", offset=-0.06, color=axes_colors['b'])

    tax.gridlines(multiple=1, linewidth=1,
                  horizontal_kwargs={'color': axes_colors['b']},
                  left_kwargs={'color': axes_colors['l']},
                  right_kwargs={'color': axes_colors['r']},
                  alpha=0.6)

    ticks = [round(i / float(scale), 1) for i in reversed(range(scale + 1))]
    tax.ticks(ticks=ticks, axis='lbr', linewidth=1, clockwise=True,
              axes_colors=axes_colors, offset=0.03)
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    #tie-lines
    prev_xs = []
    prev_ys = []
    max_phases = [0.64, 0.21, 0.15]
    print ('max phases:', max_phases)
    cnt = 0
    while max_phases[1] > 0 and cnt < 10 and tie_lines:
        pnt = (max_phases[0], max_phases[1], max_phases[2])
        tie_mix = Mixture(molecules, max_phases)
        Flash.flash_mixture(tie_mix, psi, rankine)
        diff = 0.05
        max_phases[1] -= diff
        max_phases[2] += diff
        if tie_mix.vapor == 0 or tie_mix.vapor == 1:
            print ('not in two-phase', pnt)
            continue
        cnt += 1
        fluids = tie_mix.fluids
        xs = [fluids['x'][1]*scale, fluids['x'][0]*scale, fluids['x'][2]*scale]
        ys = [fluids['y'][1]*scale, fluids['y'][0]*scale, fluids['y'][2]*scale]
        if len(prev_xs) > 0 and len(prev_ys) > 0:
            tax.plot([xs, prev_xs], linewidth=1.0, label="Two-Phase Curve", color='b')
            tax.plot([ys, prev_ys], linewidth=1.0, label="Two-Phase Curve", color='b')
        prev_xs = xs
        prev_ys = ys
        tax.plot([xs, ys], linewidth=1, label="Tie-Line", color='b')
    if specific_trace:
        z1 = specific_zFracs[1]
        z2 = specific_zFracs[0]
        z3 = specific_zFracs[2]
        specific_mix = Mixture(molecules, [z1, z2, z3])
        Flash.flash_mixture(specific_mix, psi, rankine)
        fluids = specific_mix.fluids
        zss = [fluids['z'][1] * scale, fluids['z'][0] * scale, fluids['z'][2] * scale]
        xss = [fluids['x'][1] * scale, fluids['x'][0] * scale, fluids['x'][2] * scale]
        yss = [fluids['y'][1] * scale, fluids['y'][0] * scale, fluids['y'][2] * scale]
        tax.plot([xss, zss], linewidth=2.5, label="Specific-Line", color='r')
        tax.plot([zss, yss], linewidth=2.5, label="Specific-Line", color='pink')
    if tie_lines and heatmap and finish_curve:
        prev_xs = []
        prev_ys = []
        final_mixes = mixed_data[2]
        for f_mix in final_mixes:
            print f_mix
            fluids = f_mix.fluids
            xs = [fluids['x'][1] * scale, fluids['x'][0] * scale, fluids['x'][2] * scale]
            ys = [fluids['y'][1] * scale, fluids['y'][0] * scale, fluids['y'][2] * scale]
            if len(prev_xs) > 0 and len(prev_ys) > 0:
                tax.plot([xs, prev_xs], linewidth=1.0, label="Two-Phase Curve", color='b')
                tax.plot([ys, prev_ys], linewidth=1.0, label="Two-Phase Curve", color='b')
            prev_xs = xs
            prev_ys = ys
            tax.plot([xs, ys], linewidth=1, label="Tie-Line", color='b')
    tax.show()
    return tax


storage_df = pd.DataFrame()
psi = 3000
rankine = 160+460
molecules = ['C1', 'C4', 'C10']
specific_zFracs = [ 0.1055, 0.6301, 0.2644] # i and j are switched
amount = 10000
show_heat = True
show_tie = True
show_specific = True
attempt_curve = True
tern = create_ternary(amount, show_heat, show_tie, show_specific, attempt_curve)
#storage_df.to_csv("Component_Output.csv")