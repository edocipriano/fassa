import numpy as np
import vof.normals
import vof.reconstruct

stencil = np.zeros([3,3,3])

stencil[1,1,:] = 0.5
stencil[1,2,:] = 0.5
stencil[1,0,:] = 0.5

stencil[0,2,:] = 1.
stencil[0,1,:] = 1.
stencil[0,0,:] = 1.

n = vof.normals.mycs(stencil)
print("n = ", n)

cplane = vof.reconstruct.al3d(n, 0.5)
print("cplane = ", cplane)

a2 = 0.5*(0.01)*(1./50.)

x0 = np.zeros(3)
dx = np.ones(3)

if a2 < 0:
    dx[0] = -a2
if a2 > 0:
    x0[0] = 1.-a2
    dx[0] = -a2

fl = vof.reconstruct.fl3d(n, cplane, x0, dx)
print("fl = ", fl)
