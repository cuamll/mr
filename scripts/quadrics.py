import numpy as np
import matplotlib.pyplot as plt


input_file = 'out/whatever/s_ab_total.dat'
output_file = 'out/whatever/s_ab_total_eig.dat'
s_raw = np.loadtxt(input_file)

# we want to separate out the s^{ab} tensors and keep the kvals
# this is for 2d
d = 2
kvals = s_raw[:,0:d]
s_ab_tot = s_raw[:,d:(d*(d+1))]
# reshape to be a list of dxd matrices
s_ab_tot = s_ab_tot.reshape((-1,d,d))

eigvals = np.zeros((len(s_ab_tot),d))
eigvecs = np.zeros(s_ab_tot.shape)

for i in range(len(s_ab_tot)):
    # following Steve, we want the inverse tensor
    # to ensure none of the principal axes blow up later
    eigvals[i], eigvecs[i] = np.linalg.eig(np.linalg.inv(s_ab_tot[i]))

np.savetxt(output_file, np.concatenate((kvals,eigvals,eigvecs.reshape((-1,d**2))),axis=1))

# plotting stuff. this is for 2D
test = 3000

if d == 2:
    def quadric(x,y,a,b):
        return x**2 / a + y**2 / b - 1

    xlist = np.linspace(-1.0,1.0,200)
    ylist = np.linspace(-1.0,1.0,200)
    X, Y = np.meshgrid(xlist,ylist)
    C = quadric(X, Y, eigvals[test,0], eigvals[test,1])
    plt.contour(X, Y, C, levels=[0])
    plt.grid()
    plt.axhline(0, color='black', lw=2)
    plt.axvline(0, color='black', lw=2)
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
