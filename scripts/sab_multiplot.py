#!/usr/bin/env python3
import os
import errno
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import colormaps as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from utils import s_p, do_plots
from mpl_toolkits.axes_grid1 import ImageGrid

# shouldn't need this anymore but using it for now
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
parser.add_argument("temperature",type=float, help="Temperature")
parser.add_argument("e_c",type=float, help="Core-energy constant")
parser.add_argument("dpi",type=int, help="DPI for plots")
args = parser.parse_args()
direc = args.directory
length = args.length
temp = args.temperature
core_energy = np.abs(args.e_c)
dim = 2
# this should be in the right neighbourhood for the fit to figure it out
dots = args.dpi
total_input_file = direc + '/s_ab_total.dat'
irrot_input_file = direc + '/s_ab_irrot.dat'
rot_input_file = direc + '/s_ab_rot.dat'

output_dir = direc + '/plots/'
total_output_file = output_dir + '/s_ab_total_mpl.eps'
irrot_output_file = output_dir + '/s_ab_irrot_mpl.eps'
rot_output_file = output_dir + '/s_ab_rot_mpl.eps'

mkdir_p(output_dir)
total_data = np.loadtxt(total_input_file)
irrot_data = np.loadtxt(irrot_input_file)
rot_data = np.loadtxt(rot_input_file)

kvals = total_data[:,0:dim]
small_q = kvals[0:(2*length)+1,1]
Qx, Qy = np.meshgrid(small_q, small_q)
qx = (Qx / np.pi)
qy = (Qy / np.pi)
qy_stack = np.stack((Qy / np.pi, Qy / np.pi))
# qx_stack = np.stack((Qx / np.pi, Qx / np.pi))
# qy_stack = np.stack((Qy / np.pi, Qy / np.pi))
# s_ab_tot = total_data[:,dim:(dim*(dim+1))].reshape((-1,dim,dim))
# s_ab_irrot = irrot_data[:,dim:(dim*(dim+1))].reshape((-1,dim,dim))
# s_ab_rot = rot_data[:,dim:(dim*(dim+1))].reshape((-1,dim,dim))

files = [total_input_file, irrot_input_file, rot_input_file]
output_files = [total_output_file, irrot_output_file, rot_output_file]
dats = [total_data, irrot_data, rot_data]
plot_titles = []
common = r" \; T = {:.4f}, \; u = {:.4f} $ ".format(temp, core_energy)
plot_t = r" $ S^{\alpha \beta}_{\text{total}}(\mathbf{q}), " + common
plot_titles.append(plot_t)
plot_t = r" $ S^{\alpha \beta}_{\text{irrot.}}(\mathbf{q}), " + common
plot_titles.append(plot_t)
plot_t = r" $ S^{\alpha \beta}_{\text{rot.}}(\mathbf{q}), " + common
plot_titles.append(plot_t)

plt.rc('text',usetex=True)
plt.rc('font',**{'family': 'sans-serif','sans-serif': ['Computer Modern']})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

# empty string evaluates to false; comment/uncomment to plot titles or not
# dotitle = 'True'
dotitle = ''

for i in range(len(files)):
    # Set up figure and image grid
    fig = plt.figure(figsize=(10, 8))
    if dotitle: fig.suptitle(plot_titles[i], fontsize=18)

    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                     nrows_ncols=(2,2),
                     axes_pad=0.48,
                     share_all=True,
                     cbar_location="right",
                     # cbar_mode="single",
                     cbar_mode="each",
                     cbar_size="7%",
                     cbar_pad=0.05,
                     )
    xx = dats[i][:,dim].reshape((len(small_q),len(small_q))).T
    xy = dats[i][:,dim+1].reshape((len(small_q),len(small_q))).T
    yx = dats[i][:,dim+2].reshape((len(small_q),len(small_q))).T
    yy = dats[i][:,dim+3].reshape((len(small_q),len(small_q))).T
    tens = [xx,xy,yx,yy]
    # print(xx.shape)
    # print(qx.shape)
    for j in range(4):
        ax = grid[j]
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        im = ax.contourf(qx, qy, tens[j], cmap=cm.inferno)
        # im = ax.pcolormesh(qx, qy, tens[j], cmap=cm.inferno, edgecolors='None')
        ax.cax.colorbar(im)
        ax.cax.tick_params(length=1, labelsize=9)
        ax.cax.toggle_label(True)

    # ax.cax.colorbar(im)
    # ax.cax.toggle_label(True)
    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    # plt.show()
    plt.savefig(output_files[i], format='eps', dpi=dots)
