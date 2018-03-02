import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import colormaps as cm

dir = 'T_1.650_l64_sweep/'
filename = 's_perp_total'
suffix = '.dat'
input_file = dir + filename + suffix
d = 2
length = 64
sp = 8
i = 2
j = 2
output_file = dir + 'gx_' + str(i) + '_gy_' + str(j) + filename + '.png'
s_p_raw = np.loadtxt(input_file)
kvals = s_p_raw[:,0:d]
xlist = np.linspace(kvals[0,0],kvals[-1,-1],np.sqrt(len(kvals)))
ylist = np.linspace(kvals[0,0],kvals[-1,-1],np.sqrt(len(kvals)))
X, Y = np.meshgrid(xlist, ylist)
intens = s_p_raw[:,d:d+1]
intens = intens.reshape((np.sqrt(len(intens)),np.sqrt(len(intens))))

def g_to_index(x,y):
    return (((sp + 2*x) * (length / 2)),((sp + 2*y) * (length / 2)))

def s_p(x, y, G_x, G_y):
    return (((G_x - x)*x + (G_y - y)*y)**2)/((x**2 + y**2)*((x - G_x)**2 + (y - G_y)**2))

gx, gy = g_to_index(i,j)
qx = xlist[gx-length/2:gx+1+length/2]
qy = ylist[gy-length/2:gy+1+length/2]
Qx, Qy = np.meshgrid(qx, qy)
sim_int = intens[gx-length/2:gx+1+length/2, gy-length/2:gy+1+length/2]
ana_int = s_p(Qx, Qy, qx[length/2], qx[length/2])
# central point will go to 0/0 for this expression; remove resulting NaN
ana_int[np.isnan(ana_int)] = 0.0
# check the ratios
quot = sim_int / ana_int
quot[np.isinf(quot)] = 0.0

# this is v. ugly but was the easiest way to get decent-looking tics
# and colour bars that match the sizes of the subplots
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(12,4))
ax = axes[0]
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.set_title('Simulated $ S_{\perp} $')
cs = ax.contourf(Qx / np.pi, Qy / np.pi, sim_int, cmap=cm.viridis)
fig.colorbar(cs, ax=ax)

ax = axes[1]
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.set_title('Analytical $ S_{\perp} $')
cs = ax.contourf(Qx / np.pi, Qy / np.pi, ana_int, cmap=cm.viridis)
fig.colorbar(cs, ax=ax)

ax = axes[2]
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
ax.set_title('Quotient $ S_{\perp}^{simulated} / S_{\perp}^{analytic} $')
cs = ax.contourf(Qx / np.pi, Qy / np.pi, quot, cmap=cm.viridis)
fig.colorbar(cs, ax=ax)

fig.savefig(output_file)
