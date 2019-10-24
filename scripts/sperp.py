#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Do transverse (S^{\perp}) and longitudinal (S^{\parallel})
projections of correlation tensors S_{\alpha \beta}(\mathbf{q}).
Make PDF contour plots of the results.
'''

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def plot(arr, filename, plotargs):
    fig, ax = plt.subplots()
    im = ax.imshow(arr, **plotargs)
    plt.imshow(arr, **plotargs)
    plt.xlabel(r'$ G_x $')
    plt.ylabel(r'$ G_y $')
    plt.colorbar(im)
    fig.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("-to", "--total_file", help="S_ab^total filename")
parser.add_argument("-lo", "--longitudinal_file", help="S_ab^l filename")
parser.add_argument("-tr", "--transverse_file", help="S_ab^t filename")
parser.add_argument("-le", "--length", type=int, help="system size L")

args = parser.parse_args()
path = os.path.dirname(args.total_file) + '/python_plots/'
if not os.path.exists(path):
    os.mkdir(path)

names = [os.path.basename(f) for f in [args.total_file, args.transverse_file, args.longitudinal_file]]
abbrev = [s[5:-4] for s in names]
sp = 8
lim = int(sp * (args.length/2))

tot_raw = np.loadtxt(args.total_file)
t_raw = np.loadtxt(args.transverse_file)
l_raw = np.loadtxt(args.longitudinal_file)

s_ab_tot = tot_raw[:,2:6]
s_ab_tot = s_ab_tot.reshape((2*args.length + 1, 2*args.length + 1, 2, 2))
s_ab_t = np.zeros(s_ab_tot.shape)
s_ab_l = np.zeros(s_ab_tot.shape)

t_raw = t_raw[:,2:6].reshape((args.length + 1, args.length + 1, 2, 2))
l_raw = l_raw[:,2:6].reshape((args.length + 1, args.length + 1, 2, 2))

for i in range(args.length + 1):
    for j in range(args.length + 1):
        s_ab_t[i + args.length//2, j + args.length//2, 0, 0] = t_raw[i, j, 0, 0]
        s_ab_t[i + args.length//2, j + args.length//2, 0, 1] = t_raw[i, j, 0, 1]
        s_ab_t[i + args.length//2, j + args.length//2, 1, 0] = t_raw[i, j, 1, 0]
        s_ab_t[i + args.length//2, j + args.length//2, 1, 1] = t_raw[i, j, 1, 1]
        s_ab_l[i + args.length//2, j + args.length//2, 0, 0] = l_raw[i, j, 0, 0]
        s_ab_l[i + args.length//2, j + args.length//2, 0, 1] = l_raw[i, j, 0, 1]
        s_ab_l[i + args.length//2, j + args.length//2, 1, 0] = l_raw[i, j, 1, 0]
        s_ab_l[i + args.length//2, j + args.length//2, 1, 1] = l_raw[i, j, 1, 1]

for kx in range(-args.length, args.length + 1):
    for ky in range (-args.length, args.length + 1):
        pmx = 0
        pmy = 0
        offdiag = 1
        i = kx + args.length
        j = ky + args.length
        if kx >= (args.length//2):
            pmx = -1
            offdiag = -1 * offdiag

        if kx < -(args.length//2):
            pmx = +1
            offdiag = -1 * offdiag

        if ky >= (args.length//2):
            pmy = -1
            offdiag = -1 * offdiag

        if ky < -(args.length//2):
            pmy = +1
            offdiag = -1 * offdiag

        s_ab_t[i, j, 0, 0] = s_ab_t[i + pmx * args.length, j + pmy * args.length, 0, 0]
        s_ab_t[i, j, 1, 1] = s_ab_t[i + pmx * args.length, j + pmy * args.length, 1, 1]
        s_ab_t[i, j, 0, 1] = offdiag * s_ab_t[i + pmx * args.length, j + pmy * args.length, 0, 1]
        s_ab_t[i, j, 1, 0] = offdiag * s_ab_t[i + pmx * args.length, j + pmy * args.length, 1, 0]
        s_ab_l[i, j, 0, 0] = s_ab_l[i + pmx * args.length, j + pmy * args.length, 0, 0]
        s_ab_l[i, j, 1, 1] = s_ab_l[i + pmx * args.length, j + pmy * args.length, 1, 1]
        s_ab_l[i, j, 0, 1] = offdiag * s_ab_l[i + pmx * args.length, j + pmy * args.length, 0, 1]
        s_ab_l[i, j, 1, 0] = offdiag * s_ab_l[i + pmx * args.length, j + pmy * args.length, 1, 0]

# need to declare these before starting to assign them in the loops
s_perp_tot = np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))
s_perp_t =   np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))
s_perp_l =   np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))
s_par_tot =  np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))
s_par_t =    np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))
s_par_l =    np.zeros(( (sp*args.length) + 1, (sp*args.length) + 1 ))

# project out to however many BZs we like
for p in range(-lim, lim):
    for m in range(-lim, lim):
        i = m + lim
        j = p + lim
        kxf = m * ((2*np.pi)/(args.length))
        kyf = p * ((2*np.pi)/(args.length))
        if np.abs(kxf) < 0.001 and np.abs(kyf) < 0.001:
            kn = 0.0
        else:
            kn = 1./(kxf**2 + kyf**2)

        if np.abs(p) > 2*args.length:
            ky = np.mod(p, 2*args.length) + 2*args.length
        else:
            ky = p + 2*args.length

        if np.abs(m) > 2*args.length:
            kx = np.mod(m, 2*args.length) + 2*args.length
        else:
            kx = m + 2*args.length

        if p <= -(args.length):
            ky = p
            while ky <= -(args.length):
                ky = ky + 2*args.length

            ky = ky + (args.length)
        elif p > (args.length):
            ky = p
            while ky > (args.length):
                ky = ky - 2*args.length

            ky = ky + (args.length)
        else:
            ky = p + (args.length)

        if m <= -(args.length):
            kx = m
            while kx <= -(args.length):
                kx = kx + 2*args.length

            kx = kx + (args.length)
        elif m > (args.length):
            kx = m
            while kx > (args.length):
                kx = kx - 2*args.length

            kx = kx + (args.length)
        else:
            kx = m + (args.length)

        kx = int(kx)
        ky = int(ky)
        s_perp_tot[i, j] = (1 - kxf * kxf * kn) * s_ab_tot[kx, ky, 0, 0] + \
                        ((-1) * kxf * kyf * kn) * s_ab_tot[kx, ky, 0, 1] + \
                        ((-1) * kyf * kxf * kn) * s_ab_tot[kx, ky, 1, 0] + \
                           (1 - kyf * kyf * kn) * s_ab_tot[kx, ky, 1, 1]
        s_perp_t[i, j] =     (1 - kxf * kxf * kn) * s_ab_t[kx, ky, 0, 0] + \
                          ((-1) * kxf * kyf * kn) * s_ab_t[kx, ky, 0, 1] + \
                          ((-1) * kyf * kxf * kn) * s_ab_t[kx, ky, 1, 0] + \
                             (1 - kyf * kyf * kn) * s_ab_t[kx, ky, 1, 1]
        s_perp_l[i, j] =     (1 - kxf * kxf * kn) * s_ab_l[kx, ky, 0, 0] + \
                          ((-1) * kxf * kyf * kn) * s_ab_l[kx, ky, 0, 1] + \
                          ((-1) * kyf * kxf * kn) * s_ab_l[kx, ky, 1, 0] + \
                             (1 - kyf * kyf * kn) * s_ab_l[kx, ky, 1, 1]
        s_par_tot[i, j] =      (kxf * kxf * kn) * s_ab_tot[kx, ky, 0, 0] + \
                               (kxf * kyf * kn) * s_ab_tot[kx, ky, 0, 1] + \
                               (kyf * kxf * kn) * s_ab_tot[kx, ky, 1, 0] + \
                               (kyf * kyf * kn) * s_ab_tot[kx, ky, 1, 1]
        s_par_t[i, j] =          (kxf * kxf * kn) * s_ab_t[kx, ky, 0, 0] + \
                                 (kxf * kyf * kn) * s_ab_t[kx, ky, 0, 1] + \
                                 (kyf * kxf * kn) * s_ab_t[kx, ky, 1, 0] + \
                                 (kyf * kyf * kn) * s_ab_t[kx, ky, 1, 1]
        s_par_l[i, j] =          (kxf * kxf * kn) * s_ab_l[kx, ky, 0, 0] + \
                                 (kxf * kyf * kn) * s_ab_l[kx, ky, 0, 1] + \
                                 (kyf * kxf * kn) * s_ab_l[kx, ky, 1, 0] + \
                                 (kyf * kyf * kn) * s_ab_l[kx, ky, 1, 1]

ext = [-sp/2, sp/2, -sp/2, sp/2]
plotargs = {'interpolation': 'None', 'cmap': 'inferno', 'origin': 'lower', 'extent': ext}

for v, ar in enumerate([s_perp_tot, s_perp_t, s_perp_l]):
    plot(ar, "{0}s_perp_{1}.pdf".format(path, abbrev[v]), plotargs)

for v, ar in enumerate([s_par_tot, s_par_t, s_par_l]):
    plot(ar, "{0}s_par_{1}.pdf".format(path, abbrev[v]), plotargs)

plotargs.update(extent=[-1, 1, -1, 1])
for v, ar in enumerate([s_ab_tot, s_ab_t, s_ab_l]):
    # python reshape works the opposite way to my fortran output - row-major vs column-major
    plot(ar[:,:,0,0].T, "{0}s_xx_{1}.pdf".format(path, abbrev[v]), plotargs)
    plot(ar[:,:,0,1].T, "{0}s_xy_{1}.pdf".format(path, abbrev[v]), plotargs)
    plot(ar[:,:,1,1].T, "{0}s_yy_{1}.pdf".format(path, abbrev[v]), plotargs)
