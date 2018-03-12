import numpy as np
import matplotlib.pyplot as plt


input_file = 'T_1.650_l64_sweep/s_ab_total.dat'
s_ab_output_file = 'T_1.650_l64_sweep/s_ab_total_eig.dat'
chi_output_file = 'T_1.650_l64_sweep/chi_ab_total_eig.dat'
s_raw = np.loadtxt(input_file)

# we want to separate out the s^{ab} tensors and keep the kvals
# this is for 2d
d = 2
bz = 2
length = 64
temperature = 1.650
kvals = s_raw[:,0:d]
s_ab_tot = s_raw[:,d:(d*(d+1))]
# reshape to be a list of dxd matrices
s_ab_tot = s_ab_tot.reshape((-1,d,d))
# s_ab_tot = s_ab_tot / ((length**2 * bz**2) + 1
chi_tot = s_ab_tot / temperature

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
    s_ab_eigvals[i], s_ab_eigvecs[i] = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    s_ab_eigvals[i], s_ab_eigvecs[i] = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))
    chi_eigvals[i], chi_eigvecs[i] = np.linalg.eig(np.linalg.inv(chi_tot[i]))

np.savetxt(s_ab_output_file, np.concatenate((kvals,s_ab_eigvals,s_ab_eigvecs.reshape((-1,d**2))),axis=1))
np.savetxt(chi_output_file, np.concatenate((kvals,chi_eigvals,chi_eigvecs.reshape((-1,d**2))),axis=1))

# plotting stuff. this is for 2D
test = 4300
print s_ab_inv[test,:,:]

if d == 2:
    # def quadric(x, y, a, b):
    #     return x**2 / a + y**2 / b - 1
    def quadric(x, y, a, b, c):
        return x**2 / a + y**2 / b + (2 * x * y) / c - 1

    kv = kvals[test]
    kv_str = r' $ q = ' + str(kv) + r' $'
    xlist = np.linspace(-3.0,3.0,400)
    ylist = np.linspace(-3.0,3.0,400)
    X, Y = np.meshgrid(xlist,ylist)
    C = quadric(X, Y, s_ab_inv[test,0,0], s_ab_inv[test,1,1], s_ab_inv[test,0,1])
    C2 = quadric(X, Y, chi_inv[test,0,0], chi_inv[test,1,1], chi_inv[test,0,1])
    # C = quadric(X, Y, s_ab_eigvals[test,0], s_ab_eigvals[test,1])
    # C2 = quadric(X, Y, chi_eigvals[test,0], chi_eigvals[test,1])
    fig, axes = plt.subplots(2, figsize=(4, 10))
    axes[0].contour(X, Y, C, levels=[0])
    axes[0].grid()
    axes[0].arrow(0.0,0.0,s_ab_eigvecs[test,0,0],s_ab_eigvecs[test,1,0],color='green')
    axes[0].arrow(0.0,0.0,s_ab_eigvecs[test,0,1],s_ab_eigvecs[test,1,1],color='green')
    axes[0].axhline(0, color='black', lw=2)
    axes[0].axvline(0, color='black', lw=2)
    axes[0].set_title(r'$ S^{\alpha\beta}_{tot} $ quadric, ' + kv_str)
    axes[1].contour(X, Y, C2, levels=[0])
    axes[1].grid()
    axes[1].arrow(0.0,0.0,chi_eigvecs[test,0,0],chi_eigvecs[test,0,1],color='green')
    axes[1].arrow(0.0,0.0,chi_eigvecs[test,1,0],chi_eigvecs[test,1,1],color='green')
    axes[1].axhline(0, color='black', lw=2)
    axes[1].axvline(0, color='black', lw=2)
    axes[1].set_title(r'$ \chi^{\alpha\beta}_{tot} $ quadric, ' + kv_str)
    # if you wanna see it, I guess
    plt.show()

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
    print "d is not 2 or 3, no idea what to do here"
