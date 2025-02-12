import argparse
import numpy as np
import matplotlib
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import utils


'''

quiver.py:  Takes field snapshots which my code outputs and makes them into
            matplotlib quiver plots.

'''

parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Directory containing field snapshots. Don't add the / at the end!")
parser.add_argument("length",type=int, help="System size L")
parser.add_argument("width",type=float, help="Arrow width")
parser.add_argument("dpi",type=int, help="DPI for plots")
args = parser.parse_args()
direc = args.directory
length = args.length
arrow_width = args.width
dots = args.dpi
sep = '/'
# this should eventually be an argument to the script as well
base_fn = 'snapshot_25000_'
total_file = direc + sep + base_fn + 'total'
irrot_file = direc + sep + base_fn + 'irrot'
rot_file = direc + sep + base_fn + 'rot'
charge_input_file = direc + sep + base_fn + 'charges.dat'

c_raw = np.loadtxt(charge_input_file)
pos_x = []
pos_y = []
neg_x = []
neg_y = []

for i in range(len(c_raw)):
    if (c_raw[i,2] > 0.001):
        # print c_raw[i,0], c_raw[i,1], c_raw[i,2]
        pos_x.append(c_raw[i,0])
        pos_y.append(c_raw[i,1])
    elif (c_raw[i,2] < -0.001):
        # print c_raw[i,0], c_raw[i,1], c_raw[i,2]
        neg_x.append(c_raw[i,0])
        neg_y.append(c_raw[i,1])

files = [total_file,irrot_file,rot_file]
count = 0

for fil in files:

    input_file = fil + '.dat'
    output_file = fil + '.eps'
    count = count + 1

    raw = np.loadtxt(input_file)
    X = raw[:,0]
    Y = raw[:,1]
    U = raw[:,2]
    V = raw[:,3]

    fig, ax = plt.subplots(figsize=(10,10))
    # the scale here is because the output arrows were originally huge.
    # dividing through by 2x the maximum value seems to do the trick
    # while keeping it legible, but could be tweaked
    loc = ticker.MultipleLocator(base=length/4) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=2*np.max(raw[:,2:4]))
    # plt.plot(pos_x,pos_y,'o',color=utils.rd,markeredgecolor='r')
    # plt.plot(neg_x,neg_y,'o',color=utils.blu,markeredgecolor='b')
    if count < 3:
        plt.plot(pos_x, pos_y, 'o', markersize=48/length, color=utils.rd)
        plt.plot(neg_x, neg_y, 'o', markersize=48/length, color=utils.blu)

    fig.savefig(output_file, format='eps', dpi=dots)
    ax.clear()
    plt.close(fig)
