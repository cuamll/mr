#!/usr/bin/env python3
import os
import errno
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import colormaps as cm
import utils
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from utils import s_p, do_plots

# shouldn't need this anymore but using it for now
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

np.set_printoptions(threshold=np.nan)
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
dots = args.dpi
total_input_file = direc + '/s_ab_total.dat'
irrot_input_file = direc + '/s_ab_irrot.dat'
rot_input_file = direc + '/s_ab_rot.dat'
total_perp_file = direc + '/s_perp_total.dat'
irrot_perp_file = direc + '/s_perp_irrot.dat'
output_dir = direc + '/helmholtz_twolor/'
mkdir_p(output_dir)
total_data = np.loadtxt(total_input_file)
irrot_data = np.loadtxt(irrot_input_file)
rot_data = np.loadtxt(rot_input_file)
plt.rc('text',usetex=True)
plt.rc('font',**{'family': 'sans-serif','sans-serif': ['Computer Modern']})

'''
Function definitions
'''

# still a bit confused about this
# Steve says that chi and rot_avg differ only by a factor of T
# but then says that they're both free parameters??
def irrot(dist, gamma, chi, kappa):
    """ 1d lorentzian: centred on x0, peak amplitude chi, fwhm kappa """
    # return ( ( (chi * temp) / const) / (1 + (chi) / (1 + kappa**2/(x - x0)**2)) )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        eps = 1 + (kappa**2 / dist)
        result = ( (gamma) / (1 + (chi) / (eps)) )

    return result.ravel()

def lor(dist, chi, kappa, gamma):
    """ 1d lorentzian """
    result = ((chi * kappa**2)/(kappa**2 + (dist)**2) + gamma)
    return result.ravel()

def two_lor(dist, chi1, kappa1, chi2, kappa2, gamma):
    """ two 1d lorentzians: both centred on x0, different chi and kappa """
    result = ((chi1 * kappa1**2)/(kappa1**2 + (dist)**2) + (chi2 * kappa2**2)/(kappa2**2 + (dist)**2) + gamma)
    return result.ravel()

def arctans(dist, A, B, C):
    """ can be shown to be equivalent to a distribution of Lorentzians """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = A * ((np.arctan(B * dist) - np.arctan(C * dist)) / (dist)) + gamma_init

    return result.ravel()

# these are messy, but they do the job of picking out
# the tensors for the central BZ only.
start_index = int((2 * length + 1) * (length / 2) + 1 + 0.01) - 1
end_index = int(3 * (start_index) + (length * (3/2)) + 1 + 0.01)
total_red = total_data[start_index:end_index]
irrot_red = irrot_data[start_index:end_index]
rot_red = rot_data[start_index:end_index]

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

# and take the traces, since they're invariant
rot_trace = rot_red[:,dim] + rot_red[:,(dim * (dim + 1) - 1)]
irrot_trace = irrot_red[:,dim] + irrot_red[:,(dim * (dim + 1) - 1)]
total_trace = total_red[:,dim] + total_red[:,(dim * (dim + 1) - 1)]
cut_trace = cut[:,dim] + cut[:,(dim * (dim + 1) - 1)]

# we want the q_x = -q_y cut to fit to
small_zero_index = int((len(cut)/2))

# this stems from a normalisation issue in the XY unit results.
# for Lee units results comment these four lines out.
total_trace = 0.25 * total_trace
irrot_trace = 0.25 * irrot_trace
rot_trace = 0.25 * rot_trace
cut_trace = 0.25 * cut_trace

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
# these are just to tidy up the plot call further down, very hacky
Qx, Qy = np.meshgrid(small_q, small_q)
qx_stack = np.stack((Qx / np.pi, Qx / np.pi))
qy_stack = np.stack((Qy / np.pi, Qy / np.pi))


# this actually does the fit
dists = np.sqrt(cut[:,0]**2 + cut[:,1]**2)
popt, pcov = curve_fit(two_lor, dists, cut_trace, guess, bounds=(0,np.inf))
perr = np.sqrt(np.diag(pcov))

dists = np.insert(dists,small_zero_index,0,axis=0)

# increase the mesh size to get decent resolution in the plots.
# also, note that since the fit goes over (-pi, pi) to (pi, -pi),
# there's a factor of sqrt(2) in the actual distances
fine_mesh = np.linspace(-np.pi*np.sqrt(2.), np.pi*np.sqrt(2.), 300, endpoint=True)

# lay the fit over the finer mesh
fitted_data = two_lor(fine_mesh, *popt)

tr_array = np.column_stack((Qx.flatten(),Qy.flatten(),itr.flatten()))
lor1 = lor(fine_mesh, popt[0],popt[1],popt[4])
lor2 = lor(fine_mesh, popt[2],popt[3],0.)
cut = np.column_stack((fine_mesh,fitted_data,lor1,lor2))

output_file = output_dir + 'fit_params.dat'
f = open(output_file,'w')

# print(cut)
f.write("\nIrrot cut: Qx, Qy, fitted, simulated, lor1, lor2\n")
f.write(np.array2string(cut))

# again, this one's for the two lorentzians
f.write("# chi1    chi1_err    kappa1    kappa1_err    chi2    chi2_err    kappa2    kappa2_err    gamma    gamma_err\n{:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}\n".format(popt[0],perr[0],popt[1],perr[1],popt[2],perr[2],popt[3],perr[3],popt[4],perr[4]))


# these two are latex'd titles for the plots. v ugly, i know
pat = "$ \chi_1 = {:.4f}, \kappa_1 = {:.4f} $, $\chi_2 = {:.4f}, \kappa_2 = {:.4f}, \gamma = {:.4f}, $".format(popt[0],popt[1],popt[2],popt[3],popt[4])
pat2 = "$ \chi_1 = {:.4f}, \kappa_1 = {:.4f} $,\n$\chi_2 = {:.4f}, \kappa_2 = {:.4f}, \gamma = {:.4f}, $".format(popt[0],popt[1],popt[2],popt[3],popt[4])
print(pat)

# plot the cut
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(cut[:,0] / np.pi, cut[:,1], 'o', color=utils.blu, ms=8, label='Simulation data')
plt.plot(cut[:,0] / np.pi, fitted_data, 'o-', color=utils.rd, ms=4, linewidth=2, label='Fitted data')
plt.plot(cut[:,0] / np.pi, cut[:,2], 'o-', color=utils.grn, ms=2, label='Lorentzian 1')
plt.plot(cut[:,0] / np.pi, cut[:,3], 'o-', color=utils.purp, ms=2, label='Lorentzian 2')
plt.xlabel('$ q $')
# temp_str = "{:.4f}".format(temp / (2 * np.pi))
temp_str = "{:.4f}".format(temp)
plot_title = r"Cut through $ q_y = - q_x $ for $ S^{\alpha \alpha}_{irrot.} $. T = " +  temp_str + '\n' + pat
plt.legend()
plt.title('')
plt.ylim(ymin=0,ymax=1.1*max(cut[:,1]))
plt.savefig(output_dir + 'sab_irrot_cut.eps', format='eps', dpi=dots)
plt.close()


# Now to plot the fit over the whole central BZ

# we want nans here because the neutron form factor
# is singular at the zone centre
total_trace = np.insert(total_trace,big_zero_index,np.nan,axis=0)
irrot_trace = np.insert(irrot_trace,big_zero_index,np.nan,axis=0)
rot_trace = np.insert(rot_trace,big_zero_index,np.nan,axis=0)

side = length + 1
big_dists = kvals[:,0]**2 + kvals[:,1]**2
big_dists = np.insert(big_dists,big_zero_index,0.,axis=0)
big_fit = two_lor(big_dists, *popt)
ftr = big_fit.reshape((side,side))
itr = irrot_trace.reshape((side, side))
ttr = total_trace.reshape((side, side))
rtr = rot_trace.reshape((side, side))

# plot the results
plot_titles = []
plot_t = r"Fitted $ S^{\alpha \alpha}_{irrotational}$:" + '\n' + pat
plot_titles.append(plot_t)
plot_t = r'Measured $ S^{\alpha \alpha}_{irrotational} $'
plot_titles.append(plot_t)
plot_titles = ['','']
output_file = output_dir + 'Helmholtz_sab_irrot' + '.eps'
do_plots(2, plot_titles, output_file, dots, qx_stack, qy_stack, np.stack((ftr,itr)))
ind_fit = np.unravel_index(np.nanargmax(ftr, axis=None), ftr.shape)
ind_sim = np.unravel_index(np.nanargmax(itr, axis=None), itr.shape)
print(ind_fit, ftr[ind_fit]/4, ind_sim, itr[ind_sim]/4, (ftr[ind_fit]/4)/(itr[ind_sim]/4))


# Now we've plotted the fit over the whole BZ, try simulating the projections

# try (1,1)
Gx = 1
Gy = 1
qx = small_q + (Gx * 2 * np.pi)
qy = small_q + (Gy * 2 * np.pi)
Qx, Qy = np.meshgrid(qx,qy)

d = 2
sp = 8
dim = int((length/2) + 0.01)
def g_to_index(x,y):
    '''
        Basically I'm looping over Brillouin zones below and this function
        goes from the reciprocal lattice vector \vec{G} at the zone centre
        to an array index for the simulation data.
    '''
    tup = (int(0.001 + ((sp + 2*x) * (length / 2))),int(0.001 + ((sp + 2*y) * (length / 2))))
    return tup

gx, gy = g_to_index(Gx, Gy)
total_perp_data = np.loadtxt(total_perp_file)
irrot_perp_data = np.loadtxt(irrot_perp_file)
intens = total_perp_data[:,d:d+1]
irrot_intens = irrot_perp_data[:,d:d+1]
side = int(np.sqrt(len(intens)) + 0.01) # ensure it doesn't round down too far
intens = intens.reshape((side,side))
irrot_intens = irrot_intens.reshape((side,side))
# again there's a normalisation issue
sim_int = 0.25 * intens[gx - dim : gx + dim + 1, gy - dim : gy + dim + 1]
irrot_sim_int = 0.25 * irrot_intens[gx - dim : gx + dim + 1, gy - dim : gy + dim + 1]

qx_stack = np.stack((Qx / np.pi, Qx / np.pi, Qx / np.pi))
qy_stack = np.stack((Qy / np.pi, Qy / np.pi, Qy / np.pi))

# this does f; now we need to add/subtract it in the right ways
scatt_func = s_p(Qx, Qy, Gx * 2 * np.pi, Gy * 2 * np.pi)

s_irrot_fit = ftr * (1.0 - scatt_func)
s_irrot_sim = itr * (1.0 - scatt_func)
Z_irrot = np.stack((s_irrot_fit, s_irrot_sim, irrot_sim_int))
ind_fit = np.unravel_index(np.nanargmax(s_irrot_fit, axis=None), s_irrot_fit.shape)
ind_sim1 = np.unravel_index(np.nanargmax(s_irrot_sim, axis=None), s_irrot_sim.shape)
ind_sim2 = np.unravel_index(np.argmax(irrot_sim_int, axis=None), irrot_sim_int.shape)
print(ind_fit, s_irrot_fit[ind_fit], ind_sim1, s_irrot_sim[ind_sim1], ind_sim2, irrot_sim_int[ind_sim2], irrot_sim_int[ind_sim1])

s_total_fit = (rtr * (scatt_func)) + s_irrot_fit
s_total_sim = (rtr * (scatt_func)) + s_irrot_sim
Z_total = np.stack((s_total_fit, s_total_sim, sim_int))

f.write("\n\nTotal s_perp: Qx, Qy, fitted, simulated\n")
tr_array = np.column_stack((Qx.flatten(),Qy.flatten(),s_total_fit.flatten(),sim_int.flatten()))
cut = tr_array[np.where(tr_array[:,0] + tr_array[:,1] == 4*np.pi)]
f.write(np.array2string(cut))
f.close()
# plot the cut
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
plt.plot(cut[:,0] / np.pi, cut[:,3], 'o', color=utils.blu, ms=8, label='Simulation data')
plt.plot(cut[:,0] / np.pi, cut[:,2], 'o-', color=utils.rd, ms=4, linewidth=2, label='Fitted data')
plt.xlabel('$ q $')
plt.legend()
plot_title = r"Cut through $ q_y + q_x = 4\pi $ for $ S^{\perp}_{total} $. T = " +  temp_str + '\n' + pat
plt.ylim(ymin=0,ymax=1.1*max(cut[:,3]))
plt.savefig(output_dir + 's_perp_total_cut.eps', format='eps', dpi=dots)
plt.close()

# plot of irrotational scattering
plot_titles = []
plot_t = r"Fitted $ S^{\perp}_{irrot.} $ with analytic function:" + "\n" + pat
plot_titles.append(plot_t)
plot_t = r'Measured $ S^{\perp}_{irrot.} $ with analytic function'
plot_titles.append(plot_t)
plot_t = r'Simulated $ S^{\perp}_{irrot.} $'
plot_titles.append(plot_t)
plot_titles = ['','','']
output_file = output_dir + 'Helmholtz_irrot' + '.eps'
do_plots(3, plot_titles, output_file, dots, qx_stack, qy_stack, Z_irrot)

# plot of total scattering
plot_titles = []
# plot_t = r"Fitted $ S^{\perp}_{total}: " + pat
# plot_titles.append(plot_t)
plot_t = r"Fitted $ S^{\perp}_{total} $ with analytic function:" + "\n" + pat
plot_titles.append(plot_t)
plot_t = r'Measured $ S^{\perp}_{total} $ with analytic function'
plot_titles.append(plot_t)
plot_t = r'Simulated $ S^{\perp}_{total} $'
plot_titles.append(plot_t)
plot_titles = ['','','']
output_file = output_dir + 'Helmholtz_total' + '.eps'
do_plots(3, plot_titles, output_file, dots, qx_stack, qy_stack, Z_total)
