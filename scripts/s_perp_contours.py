#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import colormaps as cm
# currently stuck with system python 2.7 so ugly mkdir
import errno
import os
import argparse

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Directory containing field snapshots. Don't add the / at the end!")
parser.add_argument("length",type=int, help="System size L")
parser.add_argument("dpi",type=int, help="DPI for plots")
args = parser.parse_args()
indir = args.directory
length = args.length
dots = args.dpi

d = 2
sp = 8
sep = '/'
filename = 's_perp_total'
suffix = '.dat'
input_file = indir + sep + filename + suffix
outdir = indir + sep + 'contour_plots/'
mkdir_p(outdir)
# os.makedirs(outdir, exist_ok=True) # python > 3.2
s_p_raw = np.loadtxt(input_file)
kvals = s_p_raw[:,0:d]
kvals = np.unique(kvals)
intens = s_p_raw[:,d:d+1]
side = int(np.sqrt(len(intens)) + 0.01) # ensure it doesn't round down too far
intens = intens.reshape((side,side))

def g_to_index(x,y):
    return (((sp + 2*x) * (length / 2)),((sp + 2*y) * (length / 2)))

# def s_p(x, y, G_x, G_y):
#     if G_x == x and G_y == y:
#         return 0.0
#     else:
#         return (((G_x - x)*x + (G_y - y)*y)**2)/((x**2 + y**2)*((x - G_x)**2 + (y - G_y)**2))

def s_p(x, y, G_x, G_y):
    return (((G_x - x)*x + (G_y - y)*y)**2)/((x**2 + y**2)*((x - G_x)**2 + (y - G_y)**2))

# we want it empty and zero-length, so we can append to it in the loop
tot_ana_int = np.array([])

for i in range(-1*(int((sp/2) + 0.01)) + 1, (int((sp/2) + 0.01))):
    for j in range(-1*(int((sp/2) + 0.01)) + 1, (int((sp/2) + 0.01))):

        output_file = outdir + 'gx_' + str(i) + '_gy_' + str(j) + '_' + filename + '.eps'
        gx, gy = g_to_index(i,j)
        dim = int((length/2) + 0.01)
        qx = kvals[gx - dim : gx + dim + 1]
        qy = kvals[gy - dim : gy + dim + 1]
        #ana_int = np.zeros((len(qx),len(qy)))
        Qx, Qy = np.meshgrid(qx, qy)
        sim_int = intens[gx - dim : gx + dim + 1, gy - dim : gy + dim + 1]
        # sim_int = intens[gx+1+length/2:gx-length/2:-1, gy-length/2:gy+1+length/2]
        ana_int = s_p(Qx, Qy, qx[dim], qy[dim])

        # for k in range(len(qx)):
        #     for m in range(len(qy)):
        #         ana_int[k,m] = s_p(qx[k], qy[m], qx[length/2], qy[length/2])

        # for some reason sim_int comes out transposed?
        # ana_int = ana_int.T
        # sim_int = sim_int.T
        # central point will go to 0/0 for this expression; remove resulting NaN
        ana_int[np.isnan(ana_int)] = 0.0
        tot_ana_int = np.append(tot_ana_int, ana_int.flatten())

        # check the ratios
        quot = sim_int / ana_int
        quot[np.isinf(quot)] = 0.0

        # this is v. ugly but was the easiest way to get decent-looking tics
        # and colour bars that match the sizes of the subplots
        fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))
        plt.rc('text',usetex=True)
        plt.rc('font',family='sans-serif')
        ax = axes[0]
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.set_title('Simulated $ S_{\perp} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, sim_int, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        ax = axes[1]
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.set_title('Analytical $ S_{\perp} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, ana_int, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        ax = axes[2]
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.set_title('Quotient $ S_{\perp}^{simulated} / S_{\perp}^{analytic} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, quot, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        fig.tight_layout()
        fig.savefig(output_file, format='eps', dpi=dots)
        plt.close(fig)

# we can print out tot_ana_int and inspect it like s_perp_total
# just can't figure out the syntax rn
# tot_ana_int = tot_ana_int.reshape((np.sqrt(len(tot_ana_int)),np.sqrt(len(tot_ana_int))))
# fig, axes = plt.subplots()
# axes.imshow(tot_ana_int, interpolation = 'bicubic', cmap=cm.viridis)
# plt.show()
