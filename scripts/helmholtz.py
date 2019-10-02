#!/usr/bin/env python3
'''
    helmholtz.py: 
'''
import os
import sys
import errno
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import utils
from utils import s_p, do_plots

'''
Function definitions
'''

def lor(dist, chi, kappa, gamma):
    """ 1d lorentzian """
    result = ((chi * kappa**2)/(kappa**2 + (dist)**2) + gamma)
    return result.ravel()

def two_lor(dist, chi1, kappa1, chi2, kappa2, gamma):
    """ two 1d lorentzians: both centred on x0, different chi and kappa """
    result = ((chi1 * kappa1**2)/(kappa1**2 + (dist)**2) +
              (chi2 * kappa2**2)/(kappa2**2 + (dist)**2) + gamma)
    return result.ravel()

parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Directory containing field snapshots. Don't add the / at the end!")
parser.add_argument("length",type=int, help="System size L")
parser.add_argument("temperature",type=float, help="Temperature")
parser.add_argument("e_c",type=float, help="Core-energy constant")
parser.add_argument("dpi",type=int, help="DPI for plots")
args = parser.parse_args()
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{sfmath}')
plt.rc('font',**{'family': 'sans-serif', 'size' : 14, 'sans-serif': ['Computer Modern']})

length = args.length
temp = args.temperature
core_energy = np.abs(args.e_c)
dots = args.dpi
total_input_file = args.directory + '/s_ab_total.dat'
irrot_input_file = args.directory + '/s_ab_l.dat'
rot_input_file = args.directory + '/s_ab_t.dat'
total_perp_file = args.directory + '/s_perp_total.dat'
irrot_perp_file = args.directory + '/s_perp_l.dat'
output_dir = args.directory + '/helmholtz_twolor/'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

raw_data = [np.loadtxt(total_input_file), np.loadtxt(irrot_input_file), np.loadtxt(rot_input_file)]

for f in raw_data:
    f = f[np.where(np.abs(f[:,0]) - 0.0001 <= np.pi)]
    f = f[np.where(np.abs(f[:,1]) - 0.0001 <= np.pi)]

total_red = total_data[np.where(np.abs(total_data[:,0]) - 0.0001 <= np.pi)]
irrot_red = irrot_data[np.where(np.abs(irrot_data[:,0]) - 0.0001 <= np.pi)]
rot_red = rot_data[np.where(np.abs(rot_data[:,0]) - 0.0001 <= np.pi)]

total_red = total_red[np.where(np.abs(total_red[:,1]) - 0.0001 <= np.pi)]
irrot_red = irrot_red[np.where(np.abs(irrot_red[:,1]) - 0.0001 <= np.pi)]
rot_red = rot_red[np.where(np.abs(rot_red[:,1]) - 0.0001 <= np.pi)]

# we want to ignore the central point in doing these fits
# because the q = 0 point is ascribed to the harmonic mode
big_zero_index = int((len(total_red)/2))
total_red = np.delete(total_red,(big_zero_index),axis=0)
irrot_red = np.delete(irrot_red,(big_zero_index),axis=0)
rot_red = np.delete(rot_red,(big_zero_index),axis=0)

# now, pick out the q_x = -q_y line
cut = irrot_red[np.where(irrot_red[:,0] + irrot_red[:,1] == 0.0)]
total_cut = total_red[np.where(total_red[:,0] + total_red[:,1] == 0.0)]

# and take the traces, since they're invariant
rot_trace = rot_red[:,dim] + rot_red[:,(dim * (dim + 1) - 1)]
irrot_trace = irrot_red[:,dim] + irrot_red[:,(dim * (dim + 1) - 1)]
total_trace = total_red[:,dim] + total_red[:,(dim * (dim + 1) - 1)]
cut_trace = cut[:,dim] + cut[:,(dim * (dim + 1) - 1)]
total_cut_trace = total_cut[:,dim] + total_cut[:,(dim * (dim + 1) - 1)]
print(total_cut_trace)

# we want the q_x = -q_y cut to fit to
small_zero_index = int((len(cut)/2))

# this should be a good estimate for gamma
rot_avg = np.sum(rot_trace) / len(rot_trace)
rot_trace[big_zero_index] = rot_avg

# initial estimate for flat background
gamma_init = rot_avg

# we'll want these later, when we project the fit over the whole BZ.
# note that I've already deleted the central point, so will
# add it back in later on.
kvals = np.array(total_red[:,0:2],dtype=float)

# these are rough guesses and may need tweaking; depends on parameters
chi1_init = temp
chi2_init = 0.75 * chi1_init
kappa1_init = 1.
kappa2_init = 2*np.pi/length

guess = np.array([chi1_init, kappa1_init, chi2_init, kappa2_init, gamma_init])

# we want to pull out the q values for the cut along q_x = -q_y
small_q = kvals[0:length+1,1]
Qx, Qy = np.meshgrid(small_q, small_q)
qx_stack = np.stack((Qx / np.pi, Qx / np.pi))
qy_stack = np.stack((Qy / np.pi, Qy / np.pi))


# actually do the fit
dists = np.sqrt(cut[:,0]**2 + cut[:,1]**2)
popt, pcov = curve_fit(two_lor, dists, cut_trace, guess, bounds=(0,np.inf))
perr = np.sqrt(np.diag(pcov))

dists = np.insert(dists,small_zero_index,0,axis=0)
cut_trace = np.insert(cut_trace,small_zero_index,0,axis=0)

# increase the mesh size to get decent resolution in the plots.
# also, note that since the fit goes over (-pi, pi) to (pi, -pi),
# there's a factor of sqrt(2) in the actual distances
fine_mesh = np.linspace(-np.pi*np.sqrt(2.), np.pi*np.sqrt(2.), 300, endpoint=True)

# lay the fit over the finer mesh
fitted_data = two_lor(fine_mesh, *popt)

lor1 = lor(fine_mesh, popt[0],popt[1],popt[4])
lor2 = lor(fine_mesh, popt[2],popt[3],0.)

cut = np.column_stack((fine_mesh,fitted_data,lor1,lor2))

output_file = output_dir + 'fit_params.dat'
f = open(output_file,'w')

# print(cut)
print(dists.shape, cut_trace.shape)
f.write("\nIrrot cut: Qx = -Qy, fitted, simulated, lor1, lor2\n")
f.write(np.array2string(cut))
f.write("\nSimulated data: Qx = -Qy, data\n")
qs = np.linspace(-np.pi*np.sqrt(2.), np.pi*np.sqrt(2.), length + 1, endpoint=True)
f.write(np.array2string(np.column_stack((qs,cut_trace)) ))

# again, this one's for the two lorentzians
f.write("# chi1    chi1_err    kappa1    kappa1_err    "
        "chi2    chi2_err    kappa2    kappa2_err    gamma    gamma_err"
        "\n{:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}     {:.4f}    "
        "{:.4f}    {:.4f}    {:.4f}    {:.4f}\n".format(popt[0],perr[0],
        popt[1],perr[1],popt[2],perr[2],popt[3],perr[3],popt[4],perr[4]))


# these two are latex'd titles for the plots. v ugly, i know
pat = "$ \chi_1 = {:.4f}, \kappa_1 = {:.4f} $,"
"$\chi_2 = {:.4f}, \kappa_2 = {:.4f}, "
"\gamma = {:.4f}, $".format(popt[0],popt[1],popt[2],popt[3],popt[4])
pat2 = "$ \chi_1 = {:.4f}, \kappa_1 = {:.4f} $,"
"\n$\chi_2 = {:.4f}, \kappa_2 = {:.4f}, "
"\gamma = {:.4f}, $".format(popt[0],popt[1],popt[2],popt[3],popt[4])
print(pat)

# plot the cut
# simulation data compared to the fit, plus the decomposed Lorentzian parts
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(cut[:,0] / np.pi, cut[:,1], 'o',
         color=utils.blu, ms=8, label='Simulation data')
plt.plot(cut[:,0] / np.pi, fitted_data, 'o-',
         color=utils.rd, ms=4, linewidth=2, label='Fitted data')
plt.plot(cut[:,0] / np.pi, cut[:,2], 'o-',
         color=utils.grn, ms=2, label='Lorentzian 1')
plt.plot(cut[:,0] / np.pi, cut[:,3], 'o-',
         color=utils.purp, ms=2, label='Lorentzian 2')
plt.xlabel('$ q $')
temp_str = "{:.4f}".format(temp)
plot_title = r"Cut through $ q_y = - q_x $ for $ "
"S^{\alpha \alpha}_{irrot.} $. T = " +  temp_str + '\n' + pat
plt.legend()
plt.title('')
plt.ylim(ymin=0,ymax=1.1*max(cut[:,1]))
plt.savefig(output_dir + 'sab_irrot_cut.eps', format='eps', dpi=dots)
plt.close()

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(total_cut[:,0] / np.pi, total_cut_trace, 'o',
         color=utils.rd, ms=8)
         # color=utils.blu, ms=8, label = r"$S^{\alpha \alpha}_{\text{total}}$")
plt.xlabel('$ q $')
temp_str = "{:.4f}".format(temp)
plot_title = r"Cut through $ q_y = - q_x $ for $ "
"S^{\alpha \alpha}_{irrot.} $. T = " +  temp_str + '\n' + pat
plt.legend()
plt.title('')
plt.ylim(ymin=0,ymax=1.1*max(total_cut_trace))
plt.savefig(output_dir + 'sab_total_cut.eps', format='eps', dpi=dots)
plt.close()
