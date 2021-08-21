import numpy as np
import vof.normals
import vof.reconstruct
import paris.fl3d

stencil3x3 = np.zeros([3,3,3])

stencil3x3[1,1,:] = 0.5
stencil3x3[1,2,:] = 0.5
stencil3x3[1,0,:] = 0.5

stencil3x3[2,0,:] = 1.0
stencil3x3[1,0,:] = 1.0
stencil3x3[0,0,:] = 1.0

alpha1here = 0.5

#normal = vof.normals.mycs(stencil3x3)
normal = np.array([1,0,0])
print("Normal = ", normal)

cplane = vof.reconstruct.al3d(normal, 0.25)
print("Cplane = ", cplane)


x0 = np.array([0.75,0,0])
deltax = np.array([0.25, 1, 1])
#newvol = vof.reconstruct.fl3d(normal, cplane, x0, deltax)
newvol = paris.fl3d.fl3d(normal, cplane, x0, deltax)
print("Newvol = ", newvol)
