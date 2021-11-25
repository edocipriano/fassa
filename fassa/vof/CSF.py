from fassa.Field import Field
from fassa.fieldtypes import *

class csf:
    """
    Compute the surface tension force using the classic Brackbill et al.
    Continuous Species Force model
    """
    def __init__(self, alpha, mesh):
        # Volume fraction field
        self.c    = c

    def setCurvature(self, K):
        self.K = K
        self.curvature_was_set = True

    def _computeNormals(self):
        mx = VectorField(self.mesh, fieldType="uStaggered")
        my = VectorField(self.mesh, fieldType="vStaggered")
        m  = VectorFaceField(self.mesh)
        self.n = VectorFaceField(self.mesh)

        # Bisogna forse fare 2 loop diversi per le staggered
        # Ci vuole un qualche punto in cui scrivere come fare i loop

        for i in range(1, self.mesh.nx):
            for j in range(1, self.mesh.ny+1):
                mx[i,j] = - (self.c[i+1,j]-self.c[i,j])/self.mesh.step

        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny):
                my[i,j] = - (self.c[i+1,j]-self.c[i,j])/self.mesh.step

        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny+1):
                m[i,j,0] = mx[i,j]
                m[i,j,1] = my[i,j]
                m[i,j,2] = 0.

        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny):
                m[i,j,0] = # TODO
                m[i,j,1] = # TODO

