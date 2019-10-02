#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Take a list of space-separated data files
and plot two columns of each on the same 3d axes.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", nargs="+",
                    help="List of data files")
parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-l", type=int, default=0,
                    help="Plot legend - 0 or 1")
parser.add_argument("-x", "--x_label",
                    help="String for x-axis")
parser.add_argument("-y", "--y_label",
                    help="String for y-axis")

args = parser.parse_args()

temps = [float(f[2:-10]) for f in args.input_file]
print(temps)
data = [np.loadtxt(f, skiprows=4) for f in args.input_file]
x = np.concatenate([f[:,0] for f in data])
z = np.concatenate([f[:,1] for f in data])
y = np.concatenate([temps[v] * np.ones(f[:,0].shape) for v, f in enumerate(data)])

# plot the results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.grid()
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{sfmath}')
plt.legend(fontsize=6)
ax.set_xlabel(r' $ \mathbf{q} $ ')
ax.set_ylabel(r' $ T (K) $ ')
ax.set_zlabel(r' $ S^{\alpha \alpha}_{L}(\mathbf{q}_x + \mathbf{q}_y = 0) $ ')
ax.scatter(x, y, z, cmap=cm.inferno, c=y)
plt.savefig(args.output_file)
plt.close()
