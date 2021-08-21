import numpy as np
from Field import Field

class ProjectionMethod:
    """
        Resolve Pressure-Velocity coupling for incompressible flows
        and constant density in 2D using the Projection Method.
    """
    def __init__(self, u, v, p, mesh, time, projectionDict):
        self.u = u
        self.v = v
        self.p = p
        self.mesh = mesh
        self.time = time
        self.projectionDict = projectionDict
        self.solveHere = Field("solveHere", None, self.mesh, 1)
        self.maxIterReached = False
        self.ut = u.cells.copy()
        self.vt = v.cells.copy()
        self._setGammaField()
        self.iter = 0


    def solveHere(self, field):
        self.solveHere = field


    def _setGammaField(self):
        """
            Set the gamma field which is used for the Poisson equation
            to account for different equations near boundaries and corners
        """
        self.gamma = Field("gamma", None, self.mesh, 1./4.)

        # Internal boundaries
        for boundary in self.mesh.boundaries:
            self.gamma.cells[self.mesh.boundaries[boundary][1]] = 1./3.

        # Internal corners
        for corner in self.mesh.corners:
            self.gamma.cells[self.mesh.corners[corner][1]] = 1./2.


    def predictorStep(self, nu):

        u = self.u.cells
        v = self.v.cells
        nu = nu.cells

        # Momentum equation for temporary u-velocity
        for i in range(1, self.mesh.nx):
            for j in range(1, self.mesh.ny+1):

                ue = 0.25*( u[i+1,j] + u[i,j]   )**2
                uw = 0.25*( u[i,j]   + u[i-1,j] )**2
                un = 0.25*( u[i,j+1] + u[i,j]   )*( v[i+1,j]+v[i,j] )
                us = 0.25*( u[i,j]   + u[i,j-1] )*( v[i+1,j-1]+v[i,j-1] )

                A = (ue - uw + un - us)/self.mesh.step
                D = (nu[i,j]/self.mesh.step**2)*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4*u[i,j])

                self.ut[i,j] = u[i,j] + self.time.deltaT*(-A+D)

        # Momentum equation for temporary v-velocity
        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny):

                vn = 0.25*( v[i,j+1]   + v[i,j]   )**2
                vs = 0.25*( v[i,j]     + v[i,j-1] )**2
                ve = 0.25*( u[i,j+1]   + u[i,j]   )*( v[i+1,j]+v[i,j] )
                vw = 0.25*( u[i-1,j+1] + u[i-1,j] )*( v[i,j]+v[i-1,j] )

                A = (vn - vs + ve - vw)/self.mesh.step
                D = (nu[i,j]/self.mesh.step**2)*(v[i+1,j] + v[i-1,j] + v[i,j+1] + v[i,j-1] - 4*v[i,j])

                self.vt[i,j] = v[i,j] + self.time.deltaT*(-A+D)


    def projectionStep(self):

        p = self.p.cells
        ut = self.ut
        vt = self.vt

        maxIter  = self.projectionDict["maxIter"]
        maxError = self.projectionDict["maxError"]
        beta     = self.projectionDict["SORCoeff"]

        for iter in range(maxIter):

            for i in range(1, self.mesh.nx+1):
                for j in range(1, self.mesh.ny+1):

                    delta = p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1]
                    S = (self.mesh.step/self.time.deltaT)*(ut[i,j] - ut[i-1,j] + vt[i,j] - vt[i,j-1])
                    p[i,j] = beta*self.gamma.cells[i,j]*( delta-S )+(1-beta)*p[i,j]

            # Estimate the error
            epsilon = 0.0

            for i in range(1, self.mesh.nx+1):
                for j in range(1, self.mesh.ny+1):
                    delta = p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1]
                    S = (self.mesh.step/self.time.deltaT)*(ut[i,j] - ut[i-1,j] + vt[i,j] - vt[i,j-1])
                    epsilon = epsilon+abs( p[i,j] - self.gamma.cells[i,j]*( delta-S ) )

            epsilon = epsilon/(self.mesh.nx*self.mesh.ny)

            # Check the error
            if (epsilon <= maxError):
                break

        if iter == maxIter-1:
            self.maxIterReached = True

        self.iter = iter


    def correctVelocity(self):

        nx = self.mesh.nx
        ny = self.mesh.ny

        self.u.cells[1:nx,1:ny+1] = self.ut[1:nx,1:ny+1] - (self.time.deltaT/self.mesh.step)*(self.p.cells[2:nx+1,1:ny+1] - self.p.cells[1:nx,1:ny+1])
        self.v.cells[1:nx+1,1:ny] = self.vt[1:nx+1,1:ny] - (self.time.deltaT/self.mesh.step)*(self.p.cells[1:nx+1,2:ny+1] - self.p.cells[1:nx+1,1:ny])


    def print(self):
        if self.maxIterReached:
            print("Step: %d - Time: %f - Poisson iterations %d --> Maximum number of iterations reached\n" %(self.time.nTimeSteps, self.time.value, self.iter))
        else:
            print("Step: %d - Time: %f - Poisson iterations %d\n" %(self.time.nTimeSteps, self.time.value, self.iter))


