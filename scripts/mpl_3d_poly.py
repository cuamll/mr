#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Take a list of space-separated data files
and plot two columns of each on the same 3d axes.
"""

import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib import cm

class MplColorHelper:
    """
    Assign rgb values from a colourmap with normalisation applied.
    Stolen from stackexchange hehe
    https://stackoverflow.com/questions/26108436/
    """
    def __init__(self, cmap_name, start_val, end_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=end_val)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
        return self.scalarMap.to_rgba(val)

def polygon_under_graph(xlist, ylist):
    """
    Construct the vertex list which defines the polygon filling the space under
    the (xlist, ylist) line graph.  Assumes the xs are in ascending order.
    Stolen from the matplotlib documentation lmao
    """
    return [(xlist[0], 0.), *zip(xlist, ylist), (xlist[-1], 0.)]


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", nargs="+",
                    help="List of data files")
parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-c", "--cmap", default='viridis', help="Colour map to use")
parser.add_argument("-x", "--x_label",
                    help="String for x-axis")
parser.add_argument("-y", "--y_label",
                    help="String for y-axis")
parser.add_argument("-z", "--z_label",
                    help="String for z-axis")

args = parser.parse_args()

temps = [float(f[2:-10]) for f in args.input_file]
data = [np.loadtxt(f, skiprows=4) for f in args.input_file]
x = np.concatenate([f[:,0] for f in data])
z = np.concatenate([f[:,1] for f in data])
y = np.concatenate([temps[v] * np.ones(f[:,0].shape) for v, f in enumerate(data)])

verts = []
COL = MplColorHelper(args.cmap, min(temps), max(temps))
fc = [COL.get_rgb(t) for t in temps]

for i in range(len(data)):
    verts.append(polygon_under_graph(data[i][:,0], data[i][:,1]))

poly = PolyCollection(verts, facecolors=fc, alpha=.4)

# plot the results
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.add_collection3d(poly, zs=temps, zdir='y')
ax.scatter(x, y, z, cmap=args.cmap, c=y, s=2)

# now some hacks to make the grid nicer - no idea why these aren't exposed
pane_colour = '#FAFAFA'
grid_params = {'grid' : {'color': "#999999", 'linewidth' : 0.4, 'linestyle' : '--'}}
ax.w_xaxis.pane.set_color(pane_colour)
ax.w_xaxis.pane.set_color(pane_colour)
ax.w_xaxis.pane.set_color(pane_colour)
ax.w_xaxis._axinfo.update(grid_params)
ax.w_yaxis._axinfo.update(grid_params)
ax.w_zaxis._axinfo.update(grid_params)

plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{sfmath}')
ax.set_xlabel(args.x_label)
ax.set_ylabel(args.y_label)
ax.set_zlabel(args.z_label)

# all temps have the same q values
ax.set_xlim(min(data[0][:,0]), max(data[0][:,0]))
ax.set_ylim(min(temps), max(temps))
ax.set_zlim(0, max(np.concatenate([f[:,1] for f in data])))
plt.savefig(args.output_file)
plt.close()
