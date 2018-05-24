#!/usr/bin/env python3
import os
import errno
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

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
parser.add_argument("dpi",type=int, help="DPI for plots")
args = parser.parse_args()
direc = args.directory
length = args.length
temp = args.temperature
dots = args.dpi
d = 2
bz = 2

input_file = direc + '/s_ab_total.dat'
output_dir = direc + '/quadrics/'
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
    s_ab_eigvals[i], s_ab_eigvecs[i] = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    chi_eigvals[i], chi_eigvecs[i] = np.linalg.eig(np.linalg.inv(chi_tot[i]))
    
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
(\pm \pi, \pm \pi) for the charge crystal case
(0,0) for the high-temperature conducting liquid phase
Those three are the three we'll pull out
"""
xpeaks = [np.pi, 0]
ypeaks = [np.pi, 0]
stringpeaks = ['pi_pi', '0_0']
latexpeaks = ['(\pi, \pi)','(0,0)']
# print s_ab_inv[test,:,:]
test = 400

if d == 2:
    for i in range(len(xpeaks)):
        def quadric(x, y, a, b):
            return (x**2 / a) + (y**2 / b) - 1
        # def quadric(x, y, a, b, c):
        #     return x**2 / a + y**2 / b + (2 * x * y) / c - 1

        cen_tuple = np.where((np.abs(kvals[:,1] - ypeaks[i]) < 0.01) & (np.abs(kvals[:,0] - xpeaks[i]) < 0.01))
        print(cen_tuple)

        kv = kvals[cen_tuple]
        kv_str = r' $ q = ' + latexpeaks[i] + r' $'
        # xlimit = s_ab_eigvecs[test,0,0]
        # ylimit = s_ab_eigvecs[test,1,1]
        xlimit = 10.0
        ylimit = 10.0
        # xlist = np.linspace(-1.2 * xlimit,1.2 * xlimit,400)
        # ylist = np.linspace(-1.2 * ylimit,1.2 * ylimit,400)
        xlist = np.linspace(-1*xlimit,xlimit,400)
        ylist = np.linspace(-1*ylimit,ylimit,400)
        X, Y = np.meshgrid(xlist,ylist)
        # C = quadric(X, Y, s_ab_inv[test,0,0], s_ab_inv[test,1,1], s_ab_inv[test,0,1])
        # C2 = quadric(X, Y, chi_inv[test,0,0], chi_inv[test,1,1], chi_inv[test,0,1])
        # C = quadric(X, Y, s_ab_eigvals[test,0], s_ab_eigvals[test,1])
        # C2 = quadric(X, Y, chi_eigvals[test,0], chi_eigvals[test,1])
        C = quadric(X, Y, s_ab_eigvals[cen_tuple,0], s_ab_eigvals[cen_tuple,1])
        C2 = quadric(X, Y, chi_eigvals[cen_tuple,0], chi_eigvals[cen_tuple,1])
        plt.rc('text',usetex=True)
        plt.rc('font',**{'family': 'sans-serif','sans-serif': ['Computer Modern']})
        fig, axes = plt.subplots(2, figsize=(4, 10))
        # slimit = C.max()
        # chilimit = C2.max()
        # slist = np.linspace(-1*slimit,slimit,400)
        # chilist = np.linspace(-1*chilimit,chilimit,400)
        # X, Y = np.meshgrid(slist,slist)
        # axes.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
        # axes.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
        axes[0].contour(X, Y, C, levels=[0])
        axes[0].grid()
        axes[0].arrow(0.0,0.0,s_ab_eigvecs[cen_tuple,0,0],s_ab_eigvecs[cen_tuple,1,0],color='green')
        axes[0].arrow(0.0,0.0,s_ab_eigvecs[cen_tuple,0,1],s_ab_eigvecs[cen_tuple,1,1],color='green')
        axes[0].axhline(0, color='black', lw=2)
        axes[0].axvline(0, color='black', lw=2)
        axes[0].set_title(r'$ S^{\alpha\beta}_{tot} $ quadric, ' + kv_str)
        # X, Y = np.meshgrid(chilist,chilist)
        axes[1].contour(X, Y, C2, levels=[0])
        axes[1].grid()
        axes[1].arrow(0.0,0.0,chi_eigvecs[cen_tuple,0,0],chi_eigvecs[cen_tuple,0,1],color='green')
        axes[1].arrow(0.0,0.0,chi_eigvecs[cen_tuple,1,0],chi_eigvecs[cen_tuple,1,1],color='green')
        axes[1].axhline(0, color='black', lw=2)
        axes[1].axvline(0, color='black', lw=2)
        axes[1].set_title(r'$ \chi^{\alpha\beta}_{tot} $ quadric, ' + kv_str)
        # if you wanna see it, I guess
        # plt.show()
        output_file = output_dir + stringpeaks[i] + '.eps'
        # fig, ax = plt.subplots()
        # peak_loc = '(' + str(int((xpeaks[i]+0.01)/np.pi)) + '$ \pi $, ' + str(int((ypeaks[i]+0.01)/np.pi)) + '$ \pi $).'
        # plt.plot(small_line[:,0] / np.pi, small_line[:,1], 'bo', label='data')
        # plt.plot(small_line[:,0] / np.pi, result.init_fit, 'k--', label='initial fit')
        # plt.plot(small_line[:,0] / np.pi, result.best_fit, 'r-', label='final fit')
        # param_title = 'Final parameters: $ \chi $ = {:.4f}, $ \kappa $ = {:.4f}, bg = {:.4f}'.format(result.params['chi'].value, result.params['kappa'].value, result.params['bg'].value)
        # temp_str = ' T = {:.4f} .'.format(temp)
        # plot_title = 'Lorentzian fit to peak in $ S_{\perp}^{total} $ at ' + peak_loc + temp_str + '\n' + param_title
        plt.legend()
        # plt.title(plot_title)
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
