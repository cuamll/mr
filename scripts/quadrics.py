#!/usr/bin/env python3
import os
import errno
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.lines import Line2D
import utils

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
dots = args.dpi
d = 2
bz = 2
# threshold for deciding if eigenvalues are degenerate
thresh = 0.02

input_file = direc + '/s_ab_total.dat'
output_dir = direc + '/quadrics_new/'
s_ab_output_file = output_dir + 's_ab_total_eig.dat'
chi_output_file = output_dir + 'chi_ab_total_eig.dat'
s_ab_small_output_file = output_dir + 's_ab_small_total_eig.dat'
s_ab_t_output_file = direc + '/s_ab_t.dat'
s_ab_l_output_file = direc + '/s_ab_l.dat'
mkdir_p(output_dir)

s_raw = np.loadtxt(input_file)
kvals = s_raw[:,0:d]
kvals_norm = np.zeros_like(kvals)
s_ab_tot = s_raw[:,d:(d*(d+1))]

for i in range(len(kvals)):
    if kvals[i,0] != 0.0 and kvals[i,1] != 0.0:
        kn = np.sqrt((kvals[i,0]**2 + kvals[i,1]**2))
    else:
        kn = 1.0

    kvals_norm[i] = kvals[i] / kn


# reshape to be a list of dxd matrices
s_ab_tot = s_ab_tot.reshape((-1,d,d))
# s_ab_tot = s_ab_tot / ((length**2 * bz**2) + 1
chi_tot = s_ab_tot / temp

s_ab_inv = np.zeros(s_ab_tot.shape)
chi_inv = np.zeros(s_ab_tot.shape)

s_ab_eigvals = np.zeros((len(s_ab_tot),d))
s_ab_eigvecs = np.zeros(s_ab_tot.shape)
chi_eigvals = np.zeros((len(s_ab_tot),d))
chi_eigvecs = np.zeros(s_ab_tot.shape)
# gonna do the helmholtz decomp in fourier space
s_ab_t = np.zeros(s_ab_tot.shape)
s_ab_l = np.zeros(s_ab_tot.shape)

for i in range(len(s_ab_tot)):
    # Doing eig of the inverse is actually probably unnecessary?
    eigvals_temp, eigvecs_temp = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    idx = np.argsort(eigvals_temp)
    eigvals_temp = eigvals_temp[idx]
    # eigenvalues of inverse tensor!
    eigvals_temp = 1. / eigvals_temp
    eigvecs_temp = eigvecs_temp[:,idx]
    s_ab_eigvals[i] = eigvals_temp
    s_ab_eigvecs[i] = eigvecs_temp

    if all(abs(kvals[i]) <= (0.00001 + (np.pi))):
        # now, we dot the normalised k's with each of the eigenvectors
        dots = np.zeros(len(eigvals_temp))
        for row in range(len(eigvecs_temp)):
            dots[row] = np.dot(kvals_norm[i],eigvecs_temp[:,row])

        # at q = 0 both dot products are zero, but the tensor is
        # a circle; just pick one. otherwise check the dot products
        if (kvals_norm[i,0] == 0.0 and kvals_norm[i,1] == 0.0):
            which_transverse = (0,)
            which_long = (1,)
        else:
            which_transverse = np.unravel_index(np.argmin(abs(dots)), dots.shape)
            which_long = np.unravel_index(np.argmax(abs(dots)), dots.shape)

        # the smallest dot product is the transverse one
        transverse_eigval = eigvals_temp[which_transverse]
        # the others are the longitudinal ones
        # this caused problems with value errors and stuff. removing
        # which_long = np.arange(len(eigvals_temp))!=which_transverse
        long_eigval = eigvals_temp[which_long]

        if abs(long_eigval - transverse_eigval) <= thresh:
            # print("Eigvals, q: ",eigvals_temp,long_eigval, transverse_eigval, abs(long_eigval - transverse_eigval), kvals[i])
            if (kvals_norm[i,0] == 0.0 and kvals_norm[i,1] == 0.0):
                print(which_long, which_transverse)
                k_long = np.array([[1.], [0.]])
                k_transverse = np.array([[0.], [1.]])
            else:
                # everything's an eigenvector of an identity matrix! so
                # we're fine to set k to be the longitudinal eigenvector
                k_long = np.array([[kvals_norm[i,0]], [kvals_norm[i,1]]])
                # and in general (b, -a) dot (a, b) = 0
                k_transverse = np.array([[kvals_norm[i,1]], [-1.*kvals_norm[i,0]]])

            k_long.shape = (2,1)
            k_transverse.shape = (2,1)
            eigvecs_temp[:,which_long] = k_long
            eigvecs_temp[:,which_transverse] = k_transverse

            s_ab_eigvecs[i] = eigvecs_temp

        # if dots[which_transverse] >= 0.1:
        #     pass
            # print("Eigvals, dot products, q, eigvecs: ",eigvals_temp,dots,kvals[i],eigvecs_temp)
            # raise Exception("Smallest eigenvalue is too big!")

        # now we need to construct the two components by picking
        # each eigenvalue individually.

        # for Q = 0, we don't do this: we take a straight average of the two
        # this is a different convention than I use for the real space
        # decomposition, where I put the harmonic mode with the irrot.
        if (kvals_norm[i,0] == 0.0 and kvals_norm[i,1] == 0.0):
            diag = np.diag(eigvals_temp/2.)
            s_ab_t[i] = eigvecs_temp @ diag @ np.linalg.inv(eigvecs_temp)
            s_ab_l[i] = eigvecs_temp @ diag @ np.linalg.inv(eigvecs_temp)
        else:
            # first the transverse component
            diag = np.zeros(len(dots))
            diag[which_transverse] = transverse_eigval
            diag = np.diag(diag)
            s_ab_t[i] = eigvecs_temp @ diag @ np.linalg.inv(eigvecs_temp)

            # and now the longitudinal
            diag = np.zeros(len(dots))
            diag[which_long] = long_eigval
            diag = np.diag(diag)
            s_ab_l[i] = eigvecs_temp @ diag @ np.linalg.inv(eigvecs_temp)

fmt_arr = ['%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f', '%+.8f']
sab_rot_fmt_arr = ['%+.10E', '%+.10E', '%+.10E', '%+.10E', '%+.10E', '%+.10E', '%+.10E', '%+.10E']
concat = np.concatenate((kvals,kvals_norm,s_ab_eigvals,s_ab_eigvecs.reshape((-1,d**2))),axis=1)
s_ab_t_concat = np.concatenate((kvals,s_ab_t.reshape((-1,d**2)),kvals_norm),axis=1)
s_ab_l_concat = np.concatenate((kvals,s_ab_l.reshape((-1,d**2)),kvals_norm),axis=1)
np.savetxt(s_ab_output_file, concat, fmt=fmt_arr)
np.savetxt(chi_output_file, np.concatenate((kvals,kvals_norm,chi_eigvals,chi_eigvecs.reshape((-1,d**2))),axis=1), fmt=fmt_arr)

# this should not be this hard, probably
size = int((length + 1)**2)
concat_small = np.empty(shape=(size,10))
s_ab_t_concat_small = np.empty(shape=(size,8))
s_ab_l_concat_small = np.empty(shape=(size,8))
j = 0

for i in range(len(concat)):
    if abs(concat[i,0]) <= (0.000001 + np.pi) and abs(concat[i,1]) <= (0.000001 + np.pi):
        concat_small[j] = concat[i,:]
        s_ab_t_concat_small[j] = s_ab_t_concat[i,:]
        s_ab_l_concat_small[j] = s_ab_l_concat[i,:]
        j = j + 1

concat_small = np.array(concat_small)
s_ab_t_concat_small = np.array(s_ab_t_concat_small)
s_ab_l_concat_small = np.array(s_ab_l_concat_small)
# np.savetxt(s_ab_small_output_file, concat_small, fmt=fmt_arr)
np.savetxt(s_ab_t_output_file, s_ab_t_concat_small, fmt=sab_rot_fmt_arr)
np.savetxt(s_ab_l_output_file, s_ab_l_concat_small, fmt=sab_rot_fmt_arr)


"""
Relevant peaks in the total S^{ab} tensor are at:
(\pm \pi, \pm \pi) for the charge crystal case,
(0,0) for the high-temperature conducting liquid phase,
The difficult thing is getting the right limits on the mesh
"""

xpeaks = [0, np.pi/8, np.pi/6, np.pi/4]
ypeaks = [0, np.pi/8, np.pi/3, np.pi/4]
stringpeaks = ['0_0', 'pi8_pi8', 'pi6_pi3', 'pi4_pi4']
latexpeaks = ['(0,0)','(\pi/8, \pi/8)','(\pi/6, \pi/3)', '(\pi/4, \pi/4)']
# print s_ab_inv[test,:,:]

if d == 2:
    for i in range(len(xpeaks)):
        def quadric(x, y, a, b):
            return (x**2 / a) + (y**2 / b) - 1

        # get the array index we want and also store the relevant k-value
        # as a string for pretty printing later
        cen_tuple = np.where((np.abs(kvals[:,1] - ypeaks[i]) < 0.01) & (np.abs(kvals[:,0] - xpeaks[i]) < 0.01))
        index = cen_tuple[0]
        print(xpeaks[i],ypeaks[i],index,s_ab_eigvals[index,:])
        kv = kvals[cen_tuple]
        kv_str = r' $ q = ' + latexpeaks[i] + r' $'
        # print(s_ab_eigvals[cen_tuple,:],chi_eigvals[cen_tuple,:])

        # by inspection we can see that the eigenvalues should correspond
        # to the intercepts of the contour, so they're our x and y limits.
        # also, for some reason, if we don't cast them to float,
        # linspace doesn't work at all and just prints x/ymax n times
        s_xmax = 1.1*float(np.round(np.sqrt(s_ab_eigvals[cen_tuple,0]),decimals=2))
        s_ymax = 1.1*float(np.round(np.sqrt(s_ab_eigvals[cen_tuple,1]),decimals=2))
        if (s_xmax > s_ymax):
            s_max = s_xmax
        else:
            s_max = s_ymax

        s_xlist = np.linspace(-s_max,s_max, num=150)
        s_ylist = np.linspace(-s_max,s_max, num=150)
        # s_xlist = np.linspace(-s_xmax,s_xmax)
        # s_ylist = np.linspace(-s_ymax,s_ymax)
        s_X, s_Y = np.meshgrid(s_xlist,s_ylist)
        C = quadric(s_X, s_Y, s_ab_eigvals[cen_tuple,0], s_ab_eigvals[cen_tuple,1])

        # chi_xmax = 1.1*float(np.round(np.sqrt(chi_eigvals[cen_tuple,0]),decimals=2))
        chi_xmax = 1.1*float(np.round(np.sqrt(chi_eigvals[cen_tuple,0]),decimals=2))
        chi_ymax = 1.1*float(np.round(np.sqrt(chi_eigvals[cen_tuple,1]),decimals=2))
        chi_xlist = np.linspace(-chi_xmax,chi_xmax, num=150)
        chi_ylist = np.linspace(-chi_ymax,chi_ymax, num=150)
        chi_X, chi_Y = np.meshgrid(chi_xlist,chi_ylist)
        C2 = quadric(chi_X, chi_Y, chi_eigvals[cen_tuple,0], chi_eigvals[cen_tuple,1])

        plt.rc('text',usetex=True)
        plt.rc('font',**{'family': 'sans-serif','sans-serif': ['Computer Modern']})
        fig, axes = plt.subplots(figsize=(10, 10))

        legend_elements = [Line2D([0], [0], color=utils.rd, lw=1, label=r' $ \bar{S}^{\alpha\beta}_{tot} $ ')]
        # legend_elements = [Line2D([0], [0], color=utils.blu, lw=1, label=r' $ \chi^{\alpha\beta}_{tot} $ '),
        #                    Line2D([0], [0], color=utils.rd, lw=1, label=r' $ S^{\alpha\beta}_{tot} $ ')]
        # axes.contour(chi_X, chi_Y, C2, colors=utils.blu, levels=[0])
        axes.grid()
        axes.axhline(0, color='black', lw=1.5)
        axes.axvline(0, color='black', lw=1.5)
        axes.contour(s_X, s_Y, C, colors=utils.rd, levels=[0])
        axes.legend(handles=legend_elements)

        param_title = 'Parameters: $ T $ = {:.4f}, $ \epsilon_c $ = {:.4f}'.format(temp, core_energy)
        # plt.title(r'$ \chi^{\alpha\beta}_{tot} $ and $ S^{\alpha\beta}_{tot} $ quadrics, ' + kv_str + '\n' + param_title)
        output_file = output_dir + stringpeaks[i] + '.eps'
        plt.legend()
        plt.savefig(output_file, format='eps', dpi=dots)
        plt.close()

elif d == 3:

    # for 3D
    def quadric(x,y,z,a,b,c):
        return (x**2 / a) + (y**2 / b) + (z**2 / c) - 1
    xlist = np.linspace(-1.0,1.0,200)
    ylist = np.linspace(-1.0,1.0,200)
    zlist = np.linspace(-1.0,1.0,200)
    X, Y, Z = np.meshgrid(xlist,ylist, zlist)
    C = quadric(X, Y, Z, eigvals[test,0], eigvals[test,1], eigvals[test,2])
    # need to check some mplot3d docs before finishing this!

else:
    print("d is not 2 or 3, no idea what to do here")
