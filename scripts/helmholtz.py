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
raw_data = [f[np.where(np.abs(f[:,0]) - 0.0001 <= np.pi)] for f in raw_data]
raw_data = [f[np.where(np.abs(f[:,1]) - 0.0001 <= np.pi)] for f in raw_data]
tr = [(f[:,2] + f[:,5]) for f in raw_data]
cut = [(f[np.where(f[:,0] + f[:,1] == 0.0)]) for f in raw_data]
tc = [(f[:,2] + f[:,5]) for f in cut]

# initial parameter guesses: never had a problem with them
chi1_init = temp
chi2_init = temp
kappa1_init = 1./length
kappa2_init = 1.
guess = np.array([chi1_init, kappa1_init, chi2_init, kappa2_init, temp])

# actually do the fit
dist = np.sqrt(cut[0][:,0]**2 + cut[0][:,1]**2)
popt, pcov = curve_fit(two_lor, dist, tc[1], bounds=(0,np.inf))
perr = np.sqrt(np.diag(pcov))

# increase the mesh size to get decent resolution in the plots.
# also, note that since the fit goes over (-pi, pi) to (pi, -pi),
# there's a factor of sqrt(2) in the actual distances
fine_mesh = np.linspace(-np.pi*np.sqrt(2.), np.pi*np.sqrt(2.), 300, endpoint=True)

# lay the fit over the finer mesh
fitted_data = two_lor(fine_mesh, *popt)
lor1 = lor(fine_mesh, popt[0], popt[1], popt[4])
lor2 = lor(fine_mesh, popt[2], popt[3], 0.)

output_file = output_dir + 'fit_params.dat'
f = open(output_file,'w')
f.write("\n# Irrot cut: Qx = -Qy, fitted, simulated, lor1, lor2\n")
f.write(np.array2string(
        np.column_stack((fine_mesh, fitted_data, lor1, lor2))))
f.write("\n# Simulated data: Qx = -Qy, data\n")
qs = np.linspace(-np.pi*np.sqrt(2.), np.pi*np.sqrt(2.), length + 1, endpoint=True)
f.write(np.array2string( np.column_stack((qs, tc[1])) ))

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

# plot the cut
# simulation data compared to the fit, plus the decomposed Lorentzian parts
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(cut[1][:,0] / np.pi, tc[1], 'o',
         color=utils.blu, ms=8, label='Simulation data')
plt.plot(fine_mesh / np.pi, fitted_data, 'o-',
         color=utils.rd, ms=4, linewidth=2, label='Fitted data')
plt.plot(fine_mesh / np.pi, lor1, 'o-',
         color=utils.grn, ms=2, label='Lorentzian 1')
plt.plot(fine_mesh / np.pi, lor2, 'o-',
         color=utils.purp, ms=2, label='Lorentzian 2')
plt.xlabel('$ q $')
temp_str = "{:.4f}".format(temp)
plot_title = r"Cut through $ q_y = - q_x $ for $ "
"S^{\alpha \alpha}_{irrot.} $. T = " +  temp_str + '\n' + pat
plt.legend()
plt.title('')
plt.ylim(ymin=0,ymax=1.1*max(fitted_data))
plt.savefig(output_dir + 'sab_irrot_cut_new.eps', format='eps', dpi=dots)
plt.close()

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(cut[0][:,0] / np.pi, tc[0], 'o',
         color=utils.rd, ms=8)
         # color=utils.blu, ms=8, label = r"$S^{\alpha \alpha}_{\text{total}}$")
plt.xlabel('$ q $')
temp_str = "{:.4f}".format(temp)
plot_title = r"Cut through $ q_y = - q_x $ for $ "
"S^{\alpha \alpha}_{irrot.} $. T = " +  temp_str + '\n' + pat
plt.legend()
plt.title('')
plt.ylim(ymin=0,ymax=1.1*max(tc[0]))
plt.savefig(output_dir + 'sab_total_cut_new.eps', format='eps', dpi=dots)
plt.close()
