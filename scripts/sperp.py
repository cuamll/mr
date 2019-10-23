#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Do transverse (S^{\perp}) and longitudinal (S^{\parallel})
projections of correlation tensors S_{\alpha \beta}(\mathbf{q}).
Make PDF contour plots of the results.
'''

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def Q_q(big_q):
    sml_q = big_q
    while np.abs(sml_q + 0.0001) > np.pi:
        if sml_q < 0.0:
            sml_q = sml_q + 2 * np.pi
        else:
            sml_q = sml_q - 2 * np.pi

    
    return sml_q

parser = argparse.ArgumentParser()
parser.add_argument("-to", "--total_file", help="S_ab^total filename")
parser.add_argument("-lo", "--longitudinal_file", help="S_ab^l filename")
parser.add_argument("-tr", "--transverse_file", help="S_ab^t filename")
parser.add_argument("-le", "--length", type=int, help="system size L")

args = parser.parse_args()
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

ext = [-sp, sp, -sp, sp]
plotargs = {'interpolation': 'None', 'cmap': 'inferno', 'extent': ext}

fig, ax = plt.subplots()
im = ax.imshow(s_perp_tot, **plotargs)
plt.imshow(s_perp_tot, **plotargs)
plt.colorbar(im)
plt.savefig("sperp_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_perp_t, **plotargs)
plt.imshow(s_perp_t, **plotargs)
plt.colorbar(im)
plt.savefig("sperp_t_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_perp_l, **plotargs)
plt.imshow(s_perp_l, **plotargs)
plt.colorbar(im)
plt.savefig("sperp_l_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_tot[:,:,0,0].T, **plotargs)
plt.imshow(s_ab_tot[:,:,0,0].T, **plotargs)
plt.colorbar(im)
plt.savefig("sxx_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_t[:,:,0,0].T, **plotargs)
plt.imshow(s_ab_t[:,:,0,0].T, **plotargs)
plt.colorbar(im)
plt.savefig("sxx_t_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_t[:,:,0,1].T, **plotargs)
plt.imshow(s_ab_t[:,:,0,1].T, **plotargs)
plt.colorbar(im)
plt.savefig("sxy_t_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_t[:,:,1,1].T, **plotargs)
plt.imshow(s_ab_t[:,:,1,1].T, **plotargs)
plt.colorbar(im)
plt.savefig("syy_t_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_tot[:,:,0,1].T, **plotargs)
plt.imshow(s_ab_tot[:,:,0,1].T, **plotargs)
plt.colorbar(im)
plt.savefig("sxy_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_tot[:,:,1,0].T, **plotargs)
plt.imshow(s_ab_tot[:,:,1,0].T, **plotargs)
plt.colorbar(im)
plt.savefig("syx_test.pdf", format="pdf")
plt.close()

fig, ax = plt.subplots()
im = ax.imshow(s_ab_tot[:,:,1,1].T, **plotargs)
plt.imshow(s_ab_tot[:,:,1,1].T, **plotargs)
plt.colorbar(im)
plt.savefig("syy_test.pdf", format="pdf")
plt.close()

# for f in files:
#     data = np.loadtxt(f)
#     xx = (data[:, 2][index]).reshape(Qx.shape)
#     np.savetxt("S_xx_reorder.dat", xx.flatten())
#     np.savetxt("S_xx_orig.dat", data[:,2])
#     xy = (data[:, 3][index]).reshape(Qx.shape)
#     yy = (data[:, 5][index]).reshape(Qx.shape)
#     s_perp = (1. - Qhatx * Qhatx) * xx + \
#             2. * (Qhatx * Qhaty) * xy + \
#             (1 - Qhaty * Qhaty) * yy
#     s_par = (Qhatx * Qhatx) * xx + \
#             2. * (Qhatx * Qhaty) * xy + \
#             (Qhaty * Qhaty) * yy
#     fig, ax = plt.subplots()
#     im = ax.imshow(s_perp, interpolation='None', cmap=cm.inferno, extent=ext)
#     plt.imshow(s_perp, interpolation='None', cmap=cm.inferno, extent=ext)
#     plt.colorbar(im)
#     plt.savefig("sperp_test.pdf", format="pdf")
#     plt.close()
#     plt.subplots()
#     plt.imshow(s_par, interpolation='None', cmap=cm.inferno, extent=ext)
#     plt.savefig("spar_test.pdf", format="pdf")
#     plt.close()
