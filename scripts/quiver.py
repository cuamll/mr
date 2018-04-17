import argparse
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

# direc = 'out/ce_T_1.500_chg_2_sign_prob/'
# direc = 'out/gce_T_1.500_e_c_-0.000_lgf_put_back_corr/'
# length = 16 
# arrow_width = 0.002

parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Directory containing field snapshots")
parser.add_argument("length",type=int, help="System size L")
parser.add_argument("width",type=int, help="Arrow width")
args = parser.parse_args()
direc = args.directory
length = args.length
arrow_width = args.width
# this should eventually be an argument to the script as well
base_fn = 'snapshot_100000_'
total_file = direc + base_fn + 'total'
irrot_file = direc + base_fn + 'irrot'
rot_file = direc + base_fn + 'rot'
charge_input_file = direc + base_fn + 'charges.dat'

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
    vals = []

    for i in range(len(c_raw)):
        red = 0.5 * (1 + c_raw[i,2])
        blue = 0.5 * (1 - c_raw[i,2])
        green = 0.0
        rgb = (red, green, blue)
        radius = np.abs(c_raw[i,2])
        x1 = c_raw[i,0]
        y1 = c_raw[i,1]
        if (np.abs(c_raw[i,2]) > 0.001):
            print x1, y1, c_raw[i,2], red, blue
            vals.append(c_raw[i,2])

        circle = Circle((x1, y1), 0.15*radius, color=rgb)
        patches.append(circle)

    # for x1 in x:
    #     for y1 in  y:
    #         # matplotlib throws the circle away if the radius is 0
    #         circle = Circle((x1 + 1, y1 + 1), 0.15*np.abs(vij[x1,y1]))
    #         patches.append(circle)

    patch_col = PatchCollection(patches, cmap=plt.cm.coolwarm)
    # patch_col.set(array=vals)
    patch_col.set_array(np.array(vals))
    # These three lines set positive charges to be red and negative to be blue
    # colours = c_raw[:,2]
    # colours_list = np.delete(colours, np.where(np.abs(colours) < 0.001))
    # col = np.array(colours_list).astype(int)
    # print "inserted array: ",col, "1s: ",(col == 1).sum(), "-1s: ",(col == -1).sum()
    # patch_col.set_array(np.array(colours_list))
    # uncomment following two lines to check the intensities of the points:
    # ret = patch_col.get_array()
    # print "returned array: ",ret, "1s: ",(ret == 1).sum(), "-1s: ",(ret == -1).sum()

    # this ensures we don't get ugly black borders round the charges
    patch_col.set_edgecolor('face')

    raw = np.loadtxt(input_file)
    X = raw[:,0]
    Y = raw[:,1]
    U = raw[:,2]
    V = raw[:,3]

    fig, ax = plt.subplots()
    # the scale here is because the output arrows were originally huge.
    # dividing through by 1.5x the maximum value seems to do the trick
    # while keeping it legible, but could be tweaked
    q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=1.5*np.max(raw[:,2:4]))

    ax.add_collection(patch_col)
    # fig.colorbar(patch_col, ax=ax)
    fig.savefig(output_file)
    ax.clear()
    plt.close(fig)
