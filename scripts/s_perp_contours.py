#!/usr/bin/env python3
import errno
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import colormaps as cm
# ignore the runtimewarning from zeros in arrays
# there's definitely a nicer way of doing this but whatever
import warnings
from utils import s_p

# currently stuck with system python 2.7 so ugly mkdir
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
irrot_input_file = indir + sep + 's_perp_irrot' + suffix
outdir = indir + sep + 'contour_plots/'
mkdir_p(outdir)
# os.makedirs(outdir, exist_ok=True) # python > 3.2
s_p_raw = np.loadtxt(input_file)
s_p_irrot_raw = np.loadtxt(irrot_input_file)
kvals = s_p_raw[:,0:d]
kvals = np.unique(kvals)
intens = s_p_raw[:,d:d+1]
irrot_intens = s_p_irrot_raw[:,d:d+1]
side = int(np.sqrt(len(intens)) + 0.01) # ensure it doesn't round down too far
intens = intens.reshape((side,side))
irrot_intens = irrot_intens.reshape((side,side))

def g_to_index(x,y):
    '''
        Basically I'm looping over Brillouin zones below and this function
        goes from the reciprocal lattice vector \vec{G} at the zone centre
        to an array index for the simulation data.
    '''
    tup = (int(0.001 + ((sp + 2*x) * (length / 2))),int(0.001 + ((sp + 2*y) * (length / 2))))
    return tup

# def s_p(x, y, G_x, G_y):
#     '''
#         This is basically the function f from Steve's correlation notes.
#         We can think of it as a kind of "projection" from the tensor 
#         (\chi / S)^{\alpha \beta}, the theoretical object, onto the
#         perpendicular component we get from neutron scattering.
#         Without superposing this function we don't see the pinch point;
#         it's the combination of the \chi/S tensor and f which gives us
#         the pinch point.

#         NB: the reason for catch warnings is that the zone centre throws
#         a warning at one index of each array, where x = G_x and y = G_y.
#         But I can't find an easy way of using meshgrid below and
#         only picking out the zone centre to skip; it complains about
#         truthiness being ambiguous for arrays. Hence, ignore the warning.
#         Bit of a hack, but whatever.

#     '''
#     # the zone centre throws a warning because x = G_x and y = G_y,
#     # but I can't find an easy way of using meshgrid below and
#     # only picking out the zone centre to skip; it complains about
#     # truthiness being ambiguous for arrays. hence, ignore the warning.
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         res = (((G_x - x)*x + (G_y - y)*y)**2)/((x**2 + y**2)*((x - G_x)**2 + (y - G_y)**2))

#     return res

# we want it empty and zero-length, so we can append to it in the loop
tot_ana_int = np.array([])

for i in range(-1*(int((sp/2) + 0.01)) + 1, (int((sp/2) + 0.01))):
    for j in range(-1*(int((sp/2) + 0.01)) + 1, (int((sp/2) + 0.01))):

        output_file = outdir + 'gx_' + str(i) + '_gy_' + str(j) + '_' + filename + '.eps'
        print("Generating contour plot for G = ({:2d}, {:2d})".format(i,j))
        gx, gy = g_to_index(i,j)
        dim = int((length/2) + 0.01)
        qx = kvals[gx - dim : gx + dim + 1]
        qy = kvals[gy - dim : gy + dim + 1]
        #ana_int = np.zeros((len(qx),len(qy)))
        Qx, Qy = np.meshgrid(qx, qy)
        sim_int = intens[gx - dim : gx + dim + 1, gy - dim : gy + dim + 1]
        irrot_sim_int = irrot_intens[gx - dim : gx + dim + 1, gy - dim : gy + dim + 1]
        # sim_int = intens[gx+1+length/2:gx-length/2:-1, gy-length/2:gy+1+length/2]
        ana_int = s_p(Qx, Qy, qx[dim], qy[dim])

        # central point will go to 0/0 for this expression; remove resulting NaN
        ana_int[np.isnan(ana_int)] = 0.0
        tot_ana_int = np.append(tot_ana_int, ana_int.flatten())

        # this was to see how close the two were
        # but because of the zone centre it's not that useful
        # quot = np.divide(sim_int, ana_int, where=ana_int!=0)

        fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(15,4))
        plt.rc('text',usetex=True)
        plt.rc('font',family='sans-serif')
        # chonk
        for chonk in range(3):
            axes[chonk].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
            axes[chonk].xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
            axes[chonk].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
            axes[chonk].yaxis.set_major_locator(tck.MultipleLocator(base=1.0))

        ax = axes[0]
        ax.set_title('Simulated $ S_{\perp}^{irrotational} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, irrot_sim_int, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        ax = axes[1]
        ax.set_title('Simulated $ S_{\perp}^{total} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, sim_int, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        ax = axes[2]
        ax.set_title('Analytical $ S_{\perp} $')
        cs = ax.contourf(Qx / np.pi, Qy / np.pi, ana_int, cmap=cm.viridis)
        fig.colorbar(cs, ax=ax)

        fig.tight_layout()
        fig.savefig(output_file, format='eps', dpi=dots)
        plt.close(fig)
