from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

import ternary


def project(p):
    SQRT3OVER2 = np.sqrt(3) / 2
    a, b, c = p
    x = a/2 + b
    y = SQRT3OVER2 * a
    return (x, y)


def matplotlib_plot(values):
    scale = 3
    points = list(ternary.helpers.simplex_iterator(scale))
    print points
    while len(points) < len(values):
        scale += 1
        points = list(ternary.helpers.simplex_iterator(scale))
    scale -= 1
    points = list(ternary.helpers.simplex_iterator(scale))
    cmap = plt.cm
    f, ax = plt.subplots(1, 1, figsize=(10, 8))

    style = 'dual-triangular'
    ticks = range(scale + 1) #[#, range(scale + 2), range(scale + 1)]
    ax.set_aspect('equal')
    ax.set_title('Ternary Diagram')
    ternary.heatmap(dict(zip(points, values)),
                    scale=scale, ax=ax,
                    cmap=cmap, vmax=len(points) + 1,
                    style=style, colorbar=True)

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    plt.show()

vals = 300

values = [np.random.randint(x, vals) for x in range(vals)]
matplotlib_plot(values)
