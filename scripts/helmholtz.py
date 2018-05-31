#!/usr/bin/env python3
import os
import errno
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from scipy.optimize import curve_fit
import colormaps as cm
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
kappa_test = 2*np.pi/length
dots = args.dpi
total_input_file = direc + '/s_ab_total.dat'
irrot_input_file = direc + '/s_ab_irrot.dat'
rot_input_file = direc + '/s_ab_rot.dat'
output_dir = direc + '/helmholtz/'
mkdir_p(output_dir)
total_data = np.loadtxt(total_input_file)
irrot_data = np.loadtxt(irrot_input_file)
rot_data = np.loadtxt(rot_input_file)

# these are messy, but they do the job of picking out
# the tensors for the central BZ only.
start_index = int((2 * length + 1) * (length / 2) + 1 + 0.01) - 1
end_index = int(3 * (start_index) + (length * (3/2)) + 1 + 0.01)
total_red = total_data[start_index:end_index]
total_red = total_red[np.where(np.abs(total_red[:,1]) - 0.0001 <= np.pi)]
irrot_red = irrot_data[start_index:end_index]
irrot_red = irrot_red[np.where(np.abs(irrot_red[:,1]) - 0.0001 <= np.pi)]
rot_red = rot_data[start_index:end_index]
rot_red = rot_red[np.where(np.abs(rot_red[:,1]) - 0.0001 <= np.pi)]

# we want to move the k = 0 point to the irrotational component
# to match Steve's correlation calculations later on
zero_index = int((len(total_red)/2))
rot_red[zero_index,:] = irrot_red[zero_index,:]
irrot_red[zero_index,:] = 0.5 * total_red[zero_index,:]

kvals = np.array(total_red[:,0:2],dtype=float)
# in 3d this would need adapting, three elements to sum
rot_trace = rot_red[:,dim] + rot_red[:,(dim * (dim + 1) - 1)]
irrot_trace = irrot_red[:,dim] + irrot_red[:,(dim * (dim + 1) - 1)]
total_trace = total_red[:,dim] + total_red[:,(dim * (dim + 1) - 1)]

# there's a hole at k = 0 from moving the component across; fill it in
# NB: Steve's gamma should be equal to rot_avg here, basically
rot_avg = np.sum(rot_trace) / len(rot_trace)
rot_trace[zero_index] = rot_avg
gamma = rot_avg
chi_start = gamma / temp

# still a bit confused about this
# Steve says that chi and rot_avg differ only by a factor of T
# but then says that they're both free parameters??
def irrot(dist, gamma, chi, kappa):
    """ 1d lorentzian: centred on x0, peak amplitude chi, fwhm kappa """
    # return ( ( (chi * temp) / const) / (1 + (chi) / (1 + kappa**2/(x - x0)**2)) )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ( (gamma) / (1 + (chi) / (1 + kappa**2/(dist))) )

    return result.ravel()

# these are all the possible values of small q, i.e. once we do Q - G
small_q = kvals[0:length+1,1]
Qx, Qy = np.meshgrid(small_q, small_q)
# this does the fit over the central BZ
dists = kvals[:,0]**2 + kvals[:,1]**2
popt, pcov = curve_fit(irrot, dists, irrot_trace)
print(popt)

fitted_data = irrot(dists, *popt)

ftr = fitted_data.reshape(( int(np.sqrt(len(fitted_data))),int(np.sqrt(len(fitted_data))) ))
itr = irrot_trace.reshape(( int(np.sqrt(len(fitted_data))),int(np.sqrt(len(fitted_data))) ))
ttr = total_trace.reshape(( int(np.sqrt(len(fitted_data))),int(np.sqrt(len(fitted_data))) ))
rtr = rot_trace.reshape(( int(np.sqrt(len(fitted_data))),int(np.sqrt(len(fitted_data))) ))

# plot the results
plt.rc('text',usetex=True)
plt.rc('font',family='sans-serif')
# fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
# # chonk
# for chonk in range(2):
#     axes[chonk].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
#     axes[chonk].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].yaxis.set_major_locator(tck.MultipleLocator(base=1.0))

# ax = axes[0]
# params_title = "\gamma = {:.4f}, \chi = {:.4f}, \kappa = {:.4f} $".format(popt[0],popt[1],popt[2])
# plot_title = r"Fitted $ S^{\alpha \beta}_{irrotational}: " + params_title
# ax.set_title(plot_title)
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, ftr, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)

# ax = axes[1]
# ax.set_title(r'Measured $ S^{\alpha \beta}_{irrotational} $')
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, itr, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)
# output_file = output_dir + 'Helmholtz' + '.eps'
# fig.tight_layout()
# fig.savefig(output_file, format='eps', dpi=dots)
# plt.close(fig)

# gonna use this multiple times
pat = "\gamma = {:.4f}, \chi = {:.4f}, \kappa = {:.4f} $".format(popt[0],popt[1],popt[2])

plot_titles = []
plt = r"Fitted $ S^{\alpha \beta}_{irrotational}: " + pat
plot_titles.append(plt)
plt = r'Measured $ S^{\alpha \beta}_{irrotational} $'
plot_titles.append(plt)
output_file = output_dir + 'Helmholtz_sab_irrot' + '.eps'
do_plots(2, plot_titles, output_file, dots, np.stack((Qx, Qx)), np.stack((Qy, Qy)), np.stack((ftr,itr)))

'''
    Now we've plotted the fit with parameters, try simulating the projections
'''

# try (1,1)
Gx = 1
Gy = 1
qx = small_q + (Gx * 2 * np.pi)
qy = small_q + (Gy * 2 * np.pi)
Qx, Qy = np.meshgrid(qx,qy)

# this does f; now we need to add/subtract it in the right ways
scatt_func = s_p(Qx, Qy, Gx * 2 * np.pi, Gy * 2 * np.pi)

s_irrot_fit = ftr * (1.0 - scatt_func)
s_irrot_sim = itr * (1.0 - scatt_func)
Z_irrot = np.stack((s_irrot_fit, s_irrot_sim))

s_total_fit = (rtr * (scatt_func)) + s_irrot_fit
s_total_sim = (rtr * (scatt_func)) + s_irrot_sim
Z_total = np.stack((s_total_fit, s_total_sim))

# fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
# chonk
# for chonk in range(2):
#     axes[chonk].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
#     axes[chonk].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].yaxis.set_major_locator(tck.MultipleLocator(base=1.0))

# ax = axes[0]
# plot_titles.append(plt)
# ax.set_title(plot_title)
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, s_irrot_fit, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)

# ax = axes[1]
# ax.set_title(r'Measured $ S^{\perp}_{irrotational} $')
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, s_irrot_sim, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)
# fig.tight_layout()
# fig.savefig(output_file, format='eps', dpi=dots)
# plt.close(fig)

plot_titles = []
plt = r"Fitted $ S^{\perp}_{irrotational}: " + pat
plot_titles.append(plt)
plt = r'Measured $ S^{\perp}_{irrotational} $'
plot_titles.append(plt)
output_file = output_dir + 'Helmholtz_irrot' + '.eps'
do_plots(2, plot_titles, output_file, dots, np.stack((Qx, Qx)), np.stack((Qy, Qy)), Z_irrot)

# fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
# # chonk
# for chonk in range(2):
#     axes[chonk].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
#     axes[chonk].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
#     axes[chonk].yaxis.set_major_locator(tck.MultipleLocator(base=1.0))

# ax = axes[0]
# plot_titles = []
# pat = "\gamma = {:.4f}, \chi = {:.4f}, \kappa = {:.4f} $".format(popt[0],popt[1],popt[2])
# plt = r"Fitted $ S^{\perp}_{total}: " + pat
# plot_titles += plt

# ax.set_title(plt)
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, s_total_fit, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)

# ax = axes[1]
# ax.set_title(r'Measured $ S^{\perp}_{total} $')
# plt = r'Measured $ S^{\perp}_{total} $'
# plot_titles += plt
# cs = ax.contourf(Qx / np.pi, Qy / np.pi, s_total_sim, cmap=cm.viridis)
# fig.colorbar(cs, ax=ax)
# fig.tight_layout()
# output_file = output_dir + 'Helmholtz_total' + '.eps'
# # fig.savefig(output_file, format='eps', dpi=dots)
# plt.close(fig)

plot_titles = []

plt = r"Fitted $ S^{\perp}_{total}: " + pat
plot_titles.append(plt)
plt = r'Measured $ S^{\perp}_{total} $'
plot_titles.append(plt)
output_file = output_dir + 'Helmholtz_total' + '.eps'
do_plots(2, plot_titles, output_file, dots, np.stack((Qx, Qx)), np.stack((Qy, Qy)), Z_total)
