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

plt.rc('text',usetex=True)
plt.rc('font',**{'family': 'sans-serif','sans-serif': ['Computer Modern']})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

dim = 3
sp = 8
dots = args.dpi
total_input_file = direc + '/s_perp_total.dat'
thresh = 0.0001
output_dir = direc + '/plots/'
hhl_output_file = direc + '/s_perp_total_hhl.dat'
p2_output_file = direc + '/s_perp_total_p2.dat'

for i in range((-sp*length/2),(sp*length/2 + 1)):
    for j in range((-sp*length/2),(sp*length/2 + 1)):
        for k in range((-sp*length/2),(sp*length/2 + 1)):
            # m = ((i - 1) / (sp*length + 1)**2)
            # p = np.mod((i - 1) / (sp*length+ 1), (sp*length + 1))
            # s = np.mod(i - 1, (sp*length + 1))
            # Qx = 2 * np.pi * (i - (1.0*sp*length)/2) / length
            # Qy = 2 * np.pi * (j - (1.0*sp*length)/2) / length
            # Qz = 2 * np.pi * (k - (1.0*sp*length)/2) / length
            Qx = 2 * np.pi * (i) / length
            Qy = 2 * np.pi * (j) / length
            Qz = 2 * np.pi * (k) / length

            Qnorm = 1.0/(Qx**2 + Qy**2 + Qz**2)

            while Qx <= -2.0 * np.pi:
                qx = Qx + 2.0 * np.pi
            while Qx > 2.0 * np.pi:
                qx = Qx - 2.0 * np.pi

            while Qy <= -2.0 * np.pi:
                qy = Qy + 2.0 * np.pi
            while Qy > 2.0 * np.pi:
                qy = Qy - 2.0 * np.pi

            while Qz <= -2.0 * np.pi:
                qz = Qz + 2.0 * np.pi
            while Qx > 2.0 * np.pi:
                qz = Qz - 2.0 * np.pi

            Qq[i,j,k,0] = Qx
            Qq[i,j,k,1] = Qy
            Qq[i,j,k,2] = Qz
            Qq[i,j,k,3] = Qnorm

            Qq[i,j,k,4] = qx
            Qq[i,j,k,5] = qy
            Qq[i,j,k,6] = qz

            # what we really want instead of little q is the
            # array index corresponding to little q: find that

            if (np.abs(k) > (length)*2):
                kz = np.mod(k,(2*length)) + (2*length)
            else:
                kz = k + (2*length)

            if (np.abs(j) > (2*length)):
                ky = np.mod(j,(2*length)) + (2*length)
            else:
                ky = j + (2*length)
 
            if (np.abs(i) > (2*length)):
                kx = np.mod(i,(2*length)) + (2*length)
            else:
                kx = i + (2*length)


            if (k < (-1*(2*(length/2)))):
                kz = k
                while (kz < (-1*(2*(length/2)))):
                    kz = kz + 2*(length)

                kz = kz + 2*(length/2)
            elif (k > (2*(length/2))):
                kz = k
                while (kz > 2*(length/2)):
                    kz = kz - 2*(length)

                kz = kz + 2*(length/2)
            else:
                kz = k + (2*(length/2))

            if (j < (-1*(2*(length/2)))):
                ky = j
                while (ky < (-1*(2*(length/2)))):
                    ky = ky + 2*(length)

                ky = ky + 2*(length/2)
           elif (j > (2*(length/2))):
                ky = j
                while (ky > 2*(length/2)):
                    ky = ky - 2*(length)

                ky = ky + 2*(length/2)
            else:
                ky = j + (2*(length/2))

            if (i < (-1*(2*(length/2)))):
                kx = i
                while (kx < (-1*(2*(length/2)))):
                    kx = kx + 2*(length)

                kx = kx + 2*(length/2)
           elif (i > (2*(length/2))):
                kx = i
                while (kx > 2*(length/2)):
                    kx = kx - 2*(length)

                kx = kx + 2*(length/2)
            else:
                kx = i + (2*(length/2))

            d = (1.0/np.sqrt(3))*np.array([1.,1.,1.])



            if (Qx + Qy + Qz <= thresh):
                s_perp_111[i,j,k] = ((1.0 - Qx*Qx*Qnorm) * s_ab[kx,ky,kz,0] +
                                    (-1.0 *Qx* Qy*Qnorm) * s_ab[kx,ky,kz,1] + 
                                    (-1.0 *Qx* Qz*Qnorm) * s_ab[kx,ky,kz,2] + 
                                    (-1.0 *Qy* Qx*Qnorm) * s_ab[kx,ky,kz,3] + 
                                    (1.0 - Qy* Qy*Qnorm) * s_ab[kx,ky,kz,4] + 
                                    (-1.0 *Qy* Qz*Qnorm) * s_ab[kx,ky,kz,5] + 
                                    (-1.0 *Qz* Qx*Qnorm) * s_ab[kx,ky,kz,6] + 
                                    (-1.0 *Qz* Qy*Qnorm) * s_ab[kx,ky,kz,7] + 
                                    (1.0 - Qz* Qz*Qnorm) * s_ab[kx,ky,kz,8])

                s_par_111[i,j,k] = ((Qx*Qx*Qnorm) * s_ab[kx,ky,kz,0] +
                                   (Qx*Qy*Qnorm) * s_ab[kx,ky,kz,1] + 
                                   (Qx*Qz*Qnorm) * s_ab[kx,ky,kz,2] + 
                                   (Qy*Qx*Qnorm) * s_ab[kx,ky,kz,3] + 
                                   (Qy*Qy*Qnorm) * s_ab[kx,ky,kz,4] + 
                                   (Qy*Qz*Qnorm) * s_ab[kx,ky,kz,5] + 
                                   (Qz*Qx*Qnorm) * s_ab[kx,ky,kz,6] + 
                                   (Qz*Qy*Qnorm) * s_ab[kx,ky,kz,7] + 
                                   (Qz*Qz*Qnorm) * s_ab[kx,ky,kz,8])

                


output_dir = direc + '/plots/'
hhl_output_file = output_dir + '/s_perp_total_hhl.dat'
p2_output_file = output_dir + '/s_perp_total_p2.dat'
# hhl_output_file = output_dir + '/s_perp_total_hhl.eps'
# hhbar0_output_file = output_dir + '/s_perp_total_hhbar0.dat'
# kkbar2k_output_file = output_dir + '/s_perp_total_kkbar2k.dat'

mkdir_p(output_dir)
total_data = np.loadtxt(total_input_file)

# we want the (hhl) plane first to plot as a heatmap
hhl = total_data[np.where(total_data[:,0] == total_data[:,1])]
# only want one h column
hhl = np.delete(hhl, 0, 1)

plane2 = total_data[np.where(np.abs(total_data[:,0] + total_data[:,1] + total_data[:,2]) <= thresh)]
# plane2 = np.delete(plane2, 2, 1)

fmt_arr = ['%+.10E', '%+.10E', '%+.10E']
np.savetxt(hhl_output_file, hhl, fmt_arr)
fmt_arr = ['%+.10E', '%+.10E', '%+.10E', '%+.10E']
np.savetxt(p2_output_file, plane2, fmt_arr)

print(len(hhl))
hvals = np.unique(hhl[:,0])
lvals = np.unique(hhl[:,1])
Qh, Qz = np.meshgrid(hvals, lvals)
qh = Qh / np.pi
qz = Qz / np.pi
side = int(np.sqrt(len(hhl)))
Z = hhl[:,2].reshape((side, side))

# files = [total_input_file, irrot_input_file, rot_input_file]
plot_titles = []
# common = r" \; T = {:.4f}, \; u = {:.4f} $ ".format(temp, core_energy)
# plot_t = r" $ S^{\perp}_{\text{total}}(\mathbf{Q})(hhl), " + common
# plot_titles.append(plot_t)

# do_plots(1, plot_titles, hhl_output_file, dots, qh, qz, Z)


# two lines of possible interest
hhbar0 = total_data[np.where(total_data[:,0] == -1.0 * total_data[:,1])]
hhbar0 = hhbar0[np.where(hhbar0[:,2] == 0)]
hhbar0 = np.delete(hhbar0, [1,2], 1)
kkbar2k = total_data[np.where(total_data[:,0] == -1.0 * total_data[:,1])]
kkbar2k = kkbar2k[np.where(kkbar2k[:,2] == 2.0*kkbar2k[:,0])]
kkbar2k = np.delete(kkbar2k, [1,2], 1)

fmt_arr = ['%+.10E', '%+.10E']
np.savetxt(hhbar0_output_file, hhbar0, fmt_arr)
np.savetxt(kkbar2k_output_file, kkbar2k, fmt_arr)

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
