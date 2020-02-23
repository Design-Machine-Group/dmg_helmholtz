import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np


def plot_grid_mesh(x, y, values, title=None, vmin=None, vmax=None):
    """
    x (np array): N x N array of x coordinates (pixel corners)
    y (np array): N x N array of y coordinates (pixel corners)
    values (np array): N-1 x N-1 array of values to color
    """
    if not vmin:
        vmin = values.min()
    if not vmax:
        vmax = values.max()
    levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
    cmap = plt.get_cmap('cividis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    fig, ax0 = plt.subplots(nrows=1)
    im = ax0.pcolormesh(x, y, values, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax0)
    if title:
        ax0.set_title(title)
    plt.show()

def plot_grid(values, title=None, vmin=None, vmax=None):
    """
    values (np array): N x N array of values to color 
    """
    if not vmin:
        vmin = values.min()
    if not vmax:
        vmax = values.max()

    cmap = plt.get_cmap('seismic')
    fig, ax0 = plt.subplots(nrows=1)
    im = ax0.pcolor(values, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax0)
    if title:
        ax0.set_title(title)
    plt.show()

if __name__ == '__main__':
    dx, dy = .03, .03
    y, x = np.mgrid[slice(1, 5 + dy, dy),
                    slice(1, 5 + dx, dx)]
    y = y.astype(float)
    # y[2][:] = 2.5
    zzz = np.sin(x)**10 + np.cos(10 + y * x) * np.cos(x)
    zzz = zzz[:-1, :-1]
    # plot_grid_mesh(x, y, zzz, title='something')
    plot_grid(zzz, title='something', vmin=0, vmax=2)