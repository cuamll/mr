import numpy as np
import matplotlib
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import matplotlib.pyplot as plt

direc = 'out/ce_T_1.000_chg_4_orig_test/'
total_input_file = direc + 'snapshot_100000_total.dat'
irrot_input_file = direc + 'snapshot_100000_irrot.dat'
rot_input_file = direc + 'snapshot_100000_rot.dat'
charge_input_file = direc + 'snapshot_100000_charges.dat'
total_output_file = direc + 'snapshot_100000_total.png'
irrot_output_file = direc + 'snapshot_100000_irrot.png'
rot_output_file = direc + 'snapshot_100000_rot.png'
charges_output_file = direc + 'snapshot_100000_charges.png'
t_raw = np.loadtxt(total_input_file)
i_raw = np.loadtxt(irrot_input_file)
r_raw = np.loadtxt(rot_input_file)
c_raw = np.loadtxt(charge_input_file)

length = 16 
arrow_width = 0.002

vij = c_raw[:,2].reshape((length, length))
no_charges = np.sum(np.abs(vij))
x = np.arange(length)
y = np.arange(length)
patches = []

for x1 in x:
    for y1 in  y:
        # matplotlib throws the circle away if the radius is 0
        circle = Circle((x1 + 1, y1 + 1), 0.15*np.abs(vij[x1,y1]), ec='none')
        patches.append(circle)

p = PatchCollection(patches)
colours = 1*c_raw[:,2]
colours_list = np.delete(colours, np.where(np.abs(colours) == 0))
p.set_array(np.array(colours_list))

X = t_raw[:,0]
Y = t_raw[:,1]
U = t_raw[:,2]
V = t_raw[:,3]

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=2*np.max(t_raw[:,2:4]))
ax.add_collection(p)
#fig.colorbar(p, ax=ax)
# plt.axis('off')
fig.savefig(total_output_file)
plt.close(fig)


patches = []

for x1 in x:
    for y1 in  y:
        # matplotlib throws the circle away if the radius is 0
        circle = Circle((x1 + 1, y1 + 1), 0.15*np.abs(vij[x1,y1]), ec='none')
        patches.append(circle)

p = PatchCollection(patches)
colours = 1*c_raw[:,2]
colours_list = np.delete(colours, np.where(np.abs(colours) == 0))
p.set_array(np.array(colours_list))

X = i_raw[:,0]
Y = i_raw[:,1]
U = i_raw[:,2]
V = i_raw[:,3]

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=2*np.max(i_raw[:,2:4]))
ax.add_collection(p)
# plt.axis('off')
#fig.colorbar(p, ax=ax)
fig.savefig(irrot_output_file)
plt.close(fig)

patches = []

for x1 in x:
    for y1 in  y:
        # matplotlib throws the circle away if the radius is 0
        circle = Circle((x1 + 1, y1 + 1), 0.15*np.abs(vij[x1,y1]), ec='none')
        patches.append(circle)

p = PatchCollection(patches)
colours = 1*c_raw[:,2]
colours_list = np.delete(colours, np.where(np.abs(colours) == 0))
p.set_array(np.array(colours_list))

X = r_raw[:,0]
Y = r_raw[:,1]
U = r_raw[:,2]
V = r_raw[:,3]

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, pivot='middle', width=arrow_width, angles='xy', scale_units='xy', scale=2*np.max(i_raw[:,2:4]))
ax.add_collection(p)
# plt.axis('off')
#fig.colorbar(p, ax=ax)
fig.savefig(rot_output_file)
plt.close(fig)
