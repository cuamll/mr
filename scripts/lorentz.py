#!/usr/bin/env python3
import os
import errno
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from lmfit import Model
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
# this should be in the right neighbourhood for the fit to figure it out
kappa_test = 2*np.pi/length
dots = args.dpi
perp_input_file = direc + '/s_perp_total.dat'
par_input_file = direc + '/s_par_total.dat'
# output_dir = direc + '/lorentz_sperp_plus_spar/'
output_dir = direc + '/lorentz_sperp/'
mkdir_p(output_dir)
perp_data = np.loadtxt(perp_input_file)
# par_data = np.loadtxt(par_input_file)
kvals = perp_data[:,0:2]
# sum_data = perp_data[:,-1] + par_data[:,-1]
sum_data = perp_data[:,-1]
data = np.column_stack((kvals,sum_data))

def lor(x, x0, chi, kappa, bg):
    """ 1d lorentzian: centred on x0, peak amplitude chi, fwhm kappa """
    return ((chi * kappa**2)/(kappa**2 + (x - x0)**2) + bg)
    # return ((chi * kappa**2)/(kappa**2 + (x - x0)**2))

"""
Relevant peaks are at (among other places, in descending order of amplidtude):
(\pm 3 \pi, \pm 5 \pi)
(\pm 5 \pi, \pm \pi)
(\pm 3 \pi, \pm 7 \pi)
Those three are the three we'll pull out
"""

xpeaks = [np.pi, 3 * np.pi, 5 * np.pi, 3 * np.pi]
ypeaks = [np.pi, 5 * np.pi, np.pi,     7 * np.pi]
stringpeaks = ['pi_pi','3pi_5pi','5pi_pi','3pi_7pi']

for i in range(len(xpeaks)):
    # pull out the relevant line
    cen_tuple = np.where((np.abs(data[:,1] - ypeaks[i]) < 0.01) & (np.abs(data[:,0] - xpeaks[i]) < 0.01))
    centre = int(cen_tuple[0])
    small_line = data[int(centre-(length/2)):int(centre+(length/2)+1),1:]

    # initialise the model and do the fit
    gmodel = Model(lor)
    result = gmodel.fit(small_line[:,1], x=small_line[:,0], x0=ypeaks[i], chi=data[centre,2], kappa=kappa_test, bg = min(small_line[:,1]))
    # result = gmodel.fit(small_line[:,1], x=small_line[:,0], x0=ypeaks[i], chi=data[centre,2], kappa=kappa_test)
    final_ys = ((result.params['chi'].value * result.params['kappa'].value**2)/(result.params['kappa'].value**2 + (small_line[:,0] - result.params['x0'].value)**2) + result.params['bg'].value)
    ybar = np.sum(small_line[:,1])/len(small_line[:,1])
    residuals = (final_ys - ybar)
    ssreg = np.sum((final_ys - ybar)**2)
    sstot = np.sum((small_line[:,1] - ybar)**2)
    rsq = 1 - ssreg / sstot
    std_err = np.sum((small_line[:,1] - final_ys)**2)/len(final_ys)

    # write out the fit report
    output_file = output_dir + stringpeaks[i] + '.fit'
    f = open(output_file,'w')
    f.write(result.fit_report())
    f.write('\nybar = {:.6f}\n'.format(ybar))
    f.write('\nk_y          ys          fitted_ys            residuals\n')
    f.write(np.array2string(np.column_stack((small_line[:,0],small_line[:,1],final_ys,residuals))))
    f.write('\nSS_reg = {:.4f}, SS_tot = {:.4f}, R^2 = {:.4f}'.format(ssreg,sstot,rsq))
    f.write('\nS = sum((ys - final_ys)^2)/len(ys) = {:.4f}'.format(std_err))
    result.plot()

    # plot the results
    output_file = output_dir + stringpeaks[i] + '.eps'
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family': 'sans-serif', 'size' : 14, 'sans-serif': ['Computer Modern']})
    fig, ax = plt.subplots()
    peak_loc = '(' + str(int((xpeaks[i]+0.01)/np.pi)) + '$ \pi $, ' + str(int((ypeaks[i]+0.01)/np.pi)) + '$ \pi $).'
    ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
    ax.tick_params(length=1, labelsize=18)
    # plt.plot(small_line[:,0] / np.pi, result.init_fit, 'k--', label='initial fit')
    plt.plot(small_line[:,0] / np.pi, small_line[:,1], 'o', color=utils.blu, ms=8, label='Simulation')
    plt.plot(small_line[:,0] / np.pi, result.best_fit, 'o-', color=utils.rd, ms=4, linewidth=2, label='Fit')
    plt.xlabel('$ Q_x $')
    # ax = plt.axes()
    # tick_locs = [centre-np.pi,centre,centre+np.pi]
    # tick_labels = [str(int(((ypeaks[i]+0.01)/np.pi) - 1)) + '$ \pi $', str(int((ypeaks[i]+0.01)/np.pi)) + '$ \pi $', str(int(((ypeaks[i]+0.01)/np.pi) + 1)) + '$ \pi $']
    # plt.xticks(tick_locs, tick_labels)
    param_title = 'Final parameters: $ \chi $ = {:.4f}, $ \kappa $ = {:.4f}, $ \gamma $ = {:.4f}'.format(result.params['chi'].value, result.params['kappa'].value, result.params['bg'].value)
    temp_str = ' $ T $ = {:.4f}, $ \epsilon_c $ = {:.4f}'.format(temp, core_energy)
    plot_title = '$ S_{\perp}^{total}$' + peak_loc + temp_str + '\n' + param_title
    plt.legend()
    # plt.title(plot_title)
    plt.savefig(output_file, format='eps', dpi=dots)
    plt.close()
