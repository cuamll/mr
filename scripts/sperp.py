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
# files = [args.total_file, args.longitudinal_file, args.transverse_file]
files = [args.total_file]
bz = 1
thresh = 0.001
ext = [-bz, bz, -bz, bz]

# doing the same procedure six times: construct arrays
Q = np.linspace(-(bz)*np.pi, (bz)*np.pi, (bz * args.length) + 1, endpoint=True)
# will use these to plot later
Qx, Qy = np.meshgrid(Q, Q)
Qhatx = Qx / np.sqrt(Qx**2 + Qy**2)
Qhaty = Qy / np.sqrt(Qx**2 + Qy**2)
sml_q = np.linspace(-np.pi, np.pi, args.length + 1, endpoint=True)
Qs = np.dstack((Qx, Qy)).reshape(-1, 2)
print(Qs)
bqsq = [Q_q(q) for q in Q]
bqQs = [Q_q(q) for q in Qs]
print(len(bqQs))
# qx = Q_q(qx)
# qy = Q_q(qy)
# qx, qy = np.meshgrid(np.array(bqsq).flatten(), np.array(bqsq).flatten())
test_indices = [np.argwhere(np.abs(sml_q - i) < thresh) for i in bqsq]
indx, indy = np.meshgrid(np.array(test_indices).flatten(), np.array(test_indices).flatten())
index = (indx.flatten() + (args.length + 1) * indy.flatten()).reshape(Qx.shape)

for f in files:
    data = np.loadtxt(f)
    xx = (data[:, 2][index]).reshape(Qx.shape).T
    print(xx)
    print(data[:,2])
    xy = (data[:, 3][index]).reshape(Qx.shape).T
    yy = (data[:, 5][index]).reshape(Qx.shape).T
    s_perp = (1. - Qhatx * Qhatx) * xx + \
            2. * (Qhatx * Qhaty) * xy + \
            (1 - Qhaty * Qhaty) * yy
    s_par = (Qhatx * Qhatx) * xx + \
            2. * (Qhatx * Qhaty) * xy + \
            (Qhaty * Qhaty) * yy
    fig, ax = plt.subplots()
    im = ax.imshow(s_perp, interpolation='None', cmap=cm.inferno, extent=ext)
    plt.imshow(s_perp, interpolation='None', cmap=cm.inferno, extent=ext)
    plt.colorbar(im)
    plt.savefig("sperp_test.pdf", format="pdf")
    plt.close()
    plt.subplots()
    plt.imshow(s_par, interpolation='None', cmap=cm.inferno, extent=ext)
    plt.savefig("spar_test.pdf", format="pdf")
    plt.close()
