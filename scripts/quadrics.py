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

input_file = direc + '/s_ab_total.dat'
output_dir = direc + '/quadrics_new/'
s_ab_output_file = output_dir + 's_ab_total_eig.dat'
chi_output_file = output_dir + 'chi_ab_total_eig.dat'
mkdir_p(output_dir)

s_raw = np.loadtxt(input_file)
kvals = s_raw[:,0:d]
s_ab_tot = s_raw[:,d:(d*(d+1))]

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

for i in range(len(s_ab_tot)):
    # following Steve, we want the inverse tensor
    # to ensure none of the principal axes blow up later
    s_ab_inv[i] = np.linalg.inv(s_ab_tot[i])
    chi_inv[i] = np.linalg.inv(chi_tot[i])
    # s_ab_eigvals[i], s_ab_eigvecs[i] = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    s_ab_eigvals_temp, s_ab_eigvecs_temp = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    chi_eigvals_temp, chi_eigvecs_temp = np.linalg.eig(np.linalg.inv(chi_tot[i]))
    idx = np.argsort(s_ab_eigvals_temp)
    s_ab_eigvals[i] = s_ab_eigvals_temp[idx]
    s_ab_eigvecs[i] = s_ab_eigvecs_temp[:,idx]
    idx = np.argsort(chi_eigvals_temp)
    chi_eigvals[i] = chi_eigvals_temp[idx]
    chi_eigvecs[i] = chi_eigvecs_temp[:,idx]
    
# print "All s_ab eigvals positive? " + str(np.all(s_ab_eigvals >= 0.0))
# print "All chi_ab eigvals positive? " + str(np.all(chi_eigvals >= 0.0))

# if not np.all(s_ab_eigvals >= 0.0):
#     print "Negative eigenvalue of s_ab_inv at",np.where(s_ab_eigvals >= 0.0)

# if not np.all(chi_eigvals >= 0.0):
#     print "Negative eigenvalue of chi_inv at",np.where(chi_eigvals >= 0.0)

np.savetxt(s_ab_output_file, np.concatenate((kvals,s_ab_eigvals,s_ab_eigvecs.reshape((-1,d**2))),axis=1))
np.savetxt(chi_output_file, np.concatenate((kvals,chi_eigvals,chi_eigvecs.reshape((-1,d**2))),axis=1))

"""
Relevant peaks in the total S^{ab} tensor are at:
(\pm \pi, \pm \pi) for the charge crystal case,
(0,0) for the high-temperature conducting liquid phase,
The difficult thing is getting the right limits on the mesh
"""
xpeaks = [0, np.pi, np.pi, np.pi/4]
ypeaks = [0, np.pi, 0, 3*np.pi/4]
stringpeaks = ['0_0', 'pi_pi', 'pi_0', 'pib4_3pib4']
latexpeaks = ['(0,0)','(\pi, \pi)','(\pi, 0)', '(\pi/4, 3\pi/4)']
# print s_ab_inv[test,:,:]

if d == 2:
    for i in range(len(xpeaks)):
        def quadric(x, y, a, b):
            return (x**2 / a) + (y**2 / b) - 1
        # def quadric(x, y, a, b, c):
        #     return x**2 / a + y**2 / b + (2 * x * y) / c - 1

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
