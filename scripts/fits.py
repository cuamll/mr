#!/Users/cgray/anaconda3/bin/python
import os
import errno
import argparse
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model

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
args = parser.parse_args()
direc = args.directory
length = args.length
# this should be in the right neighbourhood for the fit to figure it out
kappa_test = 2*np.pi/length
input_dir = '/Users/cgray/code/mr/mr/out/gce_T_0.3500_e_c_-0.5000_lorentzian_test/'
input_file = direc + '/s_perp_total.dat'
output_dir = direc + '/lorentzian_fits/'
mkdir_p(output_dir)
data = np.loadtxt(input_file)

def lor(x, x0, chi, kappa, bg):
    """ 1d lorentzian: centred on x0, peak amplitude chi, fwhm kappa """
    return ((chi * kappa**2)/(kappa**2 + (x - x0)**2) + bg)

"""
Relevant peaks are at (among other places, in descending order of amplidtude):
(\pm 3 \pi, \pm 5 \pi)
(\pm 5 \pi, \pm \pi)
(\pm 3 \pi, \pm 7 \pi)
Those three are the three we'll pull out
"""

xpeaks = [3 * np.pi, 5 * np.pi, 3 * np.pi]
ypeaks = [5 * np.pi, np.pi,     7 * np.pi]
stringpeaks = ['3pi_5pi','5pi_pi','3pi_7pi']

for i in range(len(xpeaks)):
    # pull out the relevant line
    cen_tuple = np.where((np.abs(data[:,1] - ypeaks[i]) < 0.01) & (np.abs(data[:,0] - xpeaks[i]) < 0.01))
    centre = int(cen_tuple[0])
    small_line = data[int(centre-(length/2)):int(centre+(length/2)),1:]

    # initialise the model and do the fit
    gmodel = Model(lor)
    result = gmodel.fit(small_line[:,1], x=small_line[:,0], x0=ypeaks[i], chi=data[centre,2], kappa=kappa_test, bg = min(small_line[:,1]))

    # write out the fit report
    output_file = output_dir + stringpeaks[i] + '.fit'
    f = open(output_file,'w')
    f.write(result.fit_report())

    # plot the results
    output_file = output_dir + stringpeaks[i] + '.png'
    peak_loc = '(' + str(int((xpeaks[i]+0.01)/np.pi)) + '$ \pi $, ' + str(int((ypeaks[i]+0.01)/np.pi)) + '$ \pi $)'
    plt.plot(small_line[:,0], small_line[:,1], 'bo', label='data')
    plt.plot(small_line[:,0], result.init_fit, 'k--', label='initial fit')
    plt.plot(small_line[:,0], result.best_fit, 'r-', label='final fit')
    param_title = 'Final parameters: $ \chi $ = {:.4f}, $ \kappa $ = {:.4f}, background = {:.4f}'.format(result.params['chi'].value, result.params['kappa'].value, result.params['bg'].value)
    plot_title = 'Lorentzian fit to peak in $ S_{\perp}^{total} $ at ' + peak_loc + '.\n' + param_title
    plt.legend()
    plt.title(plot_title)
    plt.savefig(output_file)
    plt.close()
