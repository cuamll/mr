#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Take a list of space-separated data files and plot two columns of each.
Assumes all kinds of things, e.g. that the data will be sorted by
the first column, the output file path is valid, etc. etc.
Also does not add any kind of legend or title. Currently doesn't work
if you generate the list via a $() in bash, either. Working on that.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", nargs="+",
                    help="List of data files")
parser.add_argument("-o", "--output_file", help="Output filename")
parser.add_argument("-c1", type=int, default=0,
                    help="Column for independent variable")
parser.add_argument("-c2", type=int, default=1,
                    help="Column for dependent variable")
parser.add_argument("-l", type=int, default=1,
                    help="Plot legend - 0 or 1")
parser.add_argument("-x", "--x_label",
                    help="String for x-axis")
parser.add_argument("-y", "--y_label",
                    help="String for y-axis")

args = parser.parse_args()

# plot the results
fig, ax = plt.subplots()
plt.grid()
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{sfmath}')
plt.legend(fontsize=6)
plt.xlabel(args.x_label)
plt.ylabel(args.y_label)

for val, f in enumerate(args.input_file):


    temp = "{:4.2f}".format(float(f[2:-15]))
    data = np.loadtxt(f)
    c1 = data[:, args.c1]
    c2 = data[:, args.c2]
    color = "C{:d}".format(val)
    if (args.l):
        label = "$ T = {} $".format(temp)
    else:
        label = None
    plt.plot(c1, c2, 'o-', color=color,label = label, ms=2)

plt.legend()
plt.savefig(args.output_file)
plt.close()
