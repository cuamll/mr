import numpy as np
import matplotlib
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import matplotlib.pyplot as plt


'''

quiver.py:  Takes field snapshots which my code outputs and makes them into
            matplotlib quiver plots.
            To do: allow command-line variables so that this can be integrated
            with my run scripts, etc.

'''

direc = 'out/ce_T_1.000_chg_4_orig_test/'
total_file = direc + 'snapshot_100000_total'
irrot_file = direc + 'snapshot_100000_irrot'
rot_file = direc + 'snapshot_100000_rot'
charge_input_file = direc + 'snapshot_100000_charges.dat'
length = 16 
arrow_width = 0.002

c_raw = np.loadtxt(charge_input_file)
vij = c_raw[:,2].reshape((length, length))
x = np.arange(length)
y = np.arange(length)

files = [total_file,irrot_file,rot_file]

for fil in files:

    input_file = fil + '.dat'
    output_file = fil + '.png'

    # I represent the charges as a set of patches
    # which are plotted over the quiver plot.

    # Technically shouldn't have to do this in the loop since there's
    # only one charge distribution, but I found that the charges get
    # plotted in different places when reusing the patch collection?
    patches = []

    for x1 in x:
        for y1 in  y:
            # matplotlib throws the circle away if the radius is 0
            circle = Circle((x1 + 1, y1 + 1), 0.15*np.abs(vij[x1,y1]), ec='none')
            patches.append(circle)

    p = PatchCollection(patches, cmap=plt.cm.bwr)
    # These three lines set positive charges to be red and negative to be blue
    colours = 1*c_raw[:,2]
    colours_list = np.delete(colours, np.where(np.abs(colours) == 0))
    p.set_array(np.array(colours_list))
    # and this ensures we don't get ugly black borders round the charges
    p.set_edgecolor('face')

    raw = np.loadtxt(input_file)

    X = raw[:,0]
    Y = raw[:,1]
    U = raw[:,2]
    V = raw[:,3]

    fig, ax = plt.subplots()
    # the scale here is because the output arrows were originally huge.
    # dividing through by twice the maximum value seems to do the trick
    # while keeping it legible, but could be tweaked
    q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=2*np.max(raw[:,2:4]))
    ax.add_collection(p)
    # fig.colorbar(p, ax=ax)
    fig.savefig(output_file)
    plt.close(fig)
