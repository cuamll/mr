#!/usr/bin/env python3
'''
    sab_multiplot.py: plot the four (really only three are independent
    but it looks weird if you just plot the three) components of the
    S^{ab} correlation tensor for the lattice Coulomb gas.
    Uses matplotlib with latex to get nice labels and titles.
'''
import os
import errno
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib import cm
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
irrot_input_file = direc + '/s_ab_l.dat'
rot_input_file = direc + '/s_ab_t.dat'

output_dir = direc + '/plots/'
total_output_file = output_dir + '/s_ab_total_mpl.eps'
irrot_output_file = output_dir + '/s_ab_irrot_mpl.eps'
rot_output_file = output_dir + '/s_ab_rot_mpl.eps'
plt.rc('text',usetex=True)
# plt.rc('font',**{'family': 'sans-serif', 'size' : 18, 'sans-serif': ['Computer Modern']})
# params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
# plt.rcParams.update(params)

mkdir_p(output_dir)
total_data = np.loadtxt(total_input_file)
irrot_data = np.loadtxt(irrot_input_file)
rot_data = np.loadtxt(rot_input_file)

kvals = total_data[:,0:dim]
# s_ab_total goes out to \pm 2\pi instead of just 1st BZ
tot_q = kvals[0:(2*length)+1,1]
qx, qy = np.meshgrid(tot_q, tot_q)
qx = (qx / np.pi)
qy = (qy / np.pi)
# s_ab_t/l only have 1st BZ
# sml_q = kvals[0:(length)+1,1]
sml_q = np.linspace(-np.pi, np.pi, length + 1, endpoint=True)
sqx, sqy = np.meshgrid(sml_q, sml_q)
# sqx = (sqx / np.pi)
# sqy = (sqy / np.pi)

files = [total_input_file, irrot_input_file, rot_input_file]
output_files = [total_output_file, irrot_output_file, rot_output_file]
dats = [total_data, irrot_data, rot_data]
qs = [[qx, qy, tot_q], [sqx, sqy, sml_q], [sqx, sqy, sml_q]]
plot_titles = []
# again, I want LaTeX text output
common = r" \; T = {:.4f}, \; u = {:.4f} $ ".format(temp, core_energy)
plot_t = r" $ S^{\alpha \beta}_{\text{total}}(\mathbf{q}), " + common
plot_titles.append(plot_t)
plot_t = r" $ S^{\alpha \beta}_{\text{irrot.}}(\mathbf{q}), " + common
plot_titles.append(plot_t)
plot_t = r" $ S^{\alpha \beta}_{\text{rot.}}(\mathbf{q}), " + common
plot_titles.append(plot_t)

# empty string evaluates to false; insert something to plot titles
dotitle = ''

for i in range(len(files)):
    # Set up figure and image grid
    fig = plt.figure(figsize=(10, 8))
    if dotitle: fig.suptitle(plot_titles[i], fontsize=16)

    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                     nrows_ncols=(2,2),
                     axes_pad=(0.80,0.30),
                     share_all=True,
                     cbar_location="right",
                     # cbar_mode="single",
                     cbar_mode="each",
                     cbar_size="9%",
                     cbar_pad=0.04,
                     )
    xx = dats[i][:,   dim].reshape((len(qs[i][2]),len(qs[i][2]))).T
    xy = dats[i][:, dim+1].reshape((len(qs[i][2]),len(qs[i][2]))).T
    yx = dats[i][:, dim+2].reshape((len(qs[i][2]),len(qs[i][2]))).T
    yy = dats[i][:, dim+3].reshape((len(qs[i][2]),len(qs[i][2]))).T
    tens = [xx,xy,yx,yy]
    for j in range(4):
        ax = grid[j]
        ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        im = ax.contourf(qs[i][0], qs[i][1], tens[j], cmap=cm.inferno)
        ax.cax.colorbar(im)
        ax.cax.tick_params(length=1, labelsize=16)
        ax.cax.toggle_label(True)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    plt.savefig(output_files[i], format='eps', dpi=dots)


fig, ax = plt.subplots()
ax.set_aspect('equal')
plt.xlabel(r'$ q_{x} $')
plt.ylabel(r'$ q_{y} $')
ltrace = dats[1][:, 2].reshape((len(qs[1][2]),len(qs[1][2]))).T + dats[1][:, 5].reshape((len(qs[1][2]),len(qs[1][2]))).T
im = ax.contourf(qs[1][0], qs[1][1], ltrace, cmap=cm.inferno)
cbar = fig.colorbar(im)
# ax.cax.tick_params(length=1, labelsize=16)
fig.tight_layout()
plt.savefig(output_dir + 's_trace_l_contour.eps', format='eps', dpi=dots)
