import numpy as np
import FieldTypes
import vof.normals
import vof.reconstruct
import paris.fl3d
import paris.al3d

class Advector:
    def __init__(self, alpha1, u, v, time):
        self.alpha1 = alpha1
        self.u = u
        self.v = v
        self.mesh = alpha1.mesh
        self.time = time
        self.setEPSC()

        # Interface normals field
        self.n = FieldTypes.VectorField(self.mesh)
        # Plane constant field
        self.cplane = FieldTypes.ScalarField(self.mesh)
        # Color function flux
        self.alphaPhi = FieldTypes.SurfaceField(self.mesh)
        # Binary color function
        self.cg = FieldTypes.ScalarField(self.mesh)

    def reconstruct(self):
        self.alpha1.updateBCs()

        # Update normals and plane constant in the same loop
        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny+1):
                stencil3x3 = self.getStencil3x3(i,j)
                self.n[i,j] = vof.normals.mycs(stencil3x3)
                mxyz = self.n[i,j,:]
                alpha1ij = self.alpha1.cells[i,j]
                self.cplane[i,j] = vof.reconstruct.al3d(mxyz, alpha1ij)


    def updateNormals(self):
        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny+1):
                stencil3x3 = self.getStencil3x3(i,j)
                self.n[i,j] = vof.normals.mycs(stencil3x3)


    def updatePlaneConstant(self):
        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny+1):
                mxyz = self.n[i,j,:]
                alpha1ij = self.alpha1.cells[i,j]
                self.cplane[i,j] = vof.reconstruct.al3d(mxyz, alpha1ij)


    def cMask(self):
        for i in range(self.mesh.nx+2):
            for j in range(self.mesh.ny+2):
                if self.alpha1.cells[i,j] > 0.5:
                    self.cg[i,j] = 1.
                else:
                    self.cg[i,j] = 0.


    def computeFluxes(self):

        # Using and operator splitting approach, every
        # direction of advection has 2 fluxes
        self.cMask()
        self.fluxes1dir("x", self.u.cells)
        self.fluxes1dir("y", self.v.cells)


    def myflux(self, direction, vel):

        mesh = self.mesh
        dt = self.time.deltaT
        alpha1 = self.alpha1.cells
        self.fluxRight = FieldTypes.ScalarField(mesh)
        self.fluxLeft  = FieldTypes.ScalarField(mesh)
        self.notEmpty = FieldTypes.BoolField(mesh, False)

        ii, jj = 0, 0
        d = -1
        if direction == "x":
            ii = 1
            dd = 0
        elif direction == "y":
            jj = 1
            dd = 1
        else:
            NotImplemented

        dxyz = mesh.step

        for j in range(1, mesh.ny+1):
            for i in range(1, mesh.nx+1):
                a2 = vel[i,j]*dt/dxyz
                a1 = vel[i-ii,j-jj]*dt/dxyz

                # Se la cella non ha liquido non fluxa niente
                if self.alpha1.cells[i,j] == 0:
                    self.fluxRight[i,j] = 0
                    self.fluxLeft[i,j] = 0
                # Se la cella è tutta liquida fluxa tutto il volumetto a2
                elif self.alpha1.cells[i,j] == 1:
                    self.fluxRight[i,j] = max(a2, 0.)    # se vDx > 0 usa il rettangolo a destra
                    self.fluxLeft[i,j]  = max(-a1, 0.)   # se vSx < 0 usa il rettangolo a sinistra
                    """
                    if a1 < 0:
                        self.fout[i,j] = -a1
                    if a2 > 0:
                        self.fout[i,j] = a2
                    """
                else:
                    stencil3x3 = self.getStencil3x3(i,j)
                    mxyz = vof.normals.mycs(stencil3x3)
                    al = paris.al3d.al3d(mxyz, self.alpha1.cells[i,j])

                    x0 = [0,0,0]
                    dx = [1,1,1]

                    #dx[dd] = a2
                    #x0[dd] = 1.-a2

                    if a1 < 0:
                        dx[dd] = -a1
                    if a2 > 0:
                        x0[dd] = 1.-a2
                        dx[dd] = a2

                    fl = paris.fl3d.fl3d(mxyz, al, x0, dx)

                    if fl != 0:
                        #self.fout[i,j] = fl # Se la cella è tagliata e fl3d trova il taglio allora fluxa la shaded area fl
                        if a1 < 0:
                            self.fluxLeft[i,j] = fl
                        if a2 > 0:
                            self.fluxRight[i,j] = fl
                    else:
                        # Se la cella è tagliata ma fl3d non trova il taglio allora dipende dalle celle vicinali:
                        # Se a destra c'è una cella liquida fluxa tutto il volumetto a2
                        if a2 > 0:
                            if self.alpha1.cells[i+ii,j+jj] == 1:
                                self.fluxRight[i,j] = a2
                            # Se a destra c'è una cella gas non fluxa liquido
                            else:
                                self.fluxRight[i,j] = 0.0
                        if a1 < 0:
                            if self.alpha1.cells[i-ii,j-jj] == 1:
                                self.fluxLeft[i,j] = -a2
                            else:
                                self.fluxLeft[i,j] = 0.0

        for j in range(1, self.mesh.ny+1):
            for i in range(1, self.mesh.nx+1):
                a2 = vel[i,j]*dt/dxyz
                a1 = vel[i-ii,j-jj]*dt/dxyz

                self.alpha1.cells[i,j] = self.alpha1.cells[i,j] - (self.fluxRight[i,j] - self.fluxLeft[i+ii,j+jj]) + \
                                (self.fluxRight[i-ii,j-jj] - self.fluxLeft[i,j]) + self.cg[i,j]*(a2-a1)


    def fluxes1dir(self, direction, vel):

        mesh = self.mesh
        dt = self.time.deltaT
        alpha1 = self.alpha1.cells
        self.vof1 = FieldTypes.ScalarField(mesh)
        self.vof2 = FieldTypes.ScalarField(mesh)
        self.vof3 = FieldTypes.ScalarField(mesh)

        ii, jj = 0, 0
        d = -1
        if direction == "x":
            ii = 1
            dd = 0
        elif direction == "y":
            jj = 1
            dd = 1
        else:
            NotImplemented

        dxyz = mesh.step

        for j in range(1, mesh.ny+1):
            for i in range(1, mesh.nx+1):
                a2 = vel[i,j]*dt/dxyz
                a1 = vel[i-ii,j-jj]*dt/dxyz
                # default: fluxes=0. (good also for c=0.)
                self.vof1[i,j] = 0.0
                self.vof2[i,j] = 0.0
                self.vof3[i,j] = 0.0
                # c = 1
                if alpha1[i,j] == 1.:
                    self.vof1[i,j] = max(-a1, 0.0)
                    self.vof3[i,j] = max(a2, 0.0)
                # 0. < c < 1.
                elif alpha1[i,j] > 0.:
                    stencil3x3 = self.getStencil3x3(i,j)
                    mxyz = vof.normals.mycs(stencil3x3)
                    alconst = paris.al3d.al3d(mxyz, self.alpha1.cells[i,j])
                    # Eulerian advection
                    x0 = [0,0,0]
                    dx = [1,1,1]
                    if a1 < 0.:
                        dx[dd] = -a1
                        #self.vof1[i,j] = paris.fl3d.fl3d(mxyz, alconst, x0, deltax)
                    if a2 > 0.:
                        x0[dd] = 1.-a2
                        dx[dd] = a2
                        #self.vof3[i,j] = paris.fl3d.fl3d(mxyz, alconst, x0, deltax)

                    fl = paris.fl3d.fl3d(mxyz, alconst, x0, dx)

                    if fl != 0:
                        #self.fout[i,j] = fl # Se la cella è tagliata e fl3d trova il taglio allora fluxa la shaded area fl
                        if a1 < 0:
                            self.vof1[i,j] = fl
                        if a2 > 0:
                            self.vof3[i,j] = fl
                    else:
                        # Se la cella è tagliata ma fl3d non trova il taglio allora dipende dalle celle vicinali:
                        # Se a destra c'è una cella liquida fluxa tutto il volumetto a2
                        if a2 > 0:
                            if self.alpha1.cells[i+ii,j+jj] == 1: # Meglio > 0.5 ?
                                self.vof3[i,j] = a2
                            # Se a destra c'è una cella gas non fluxa liquido
                            else:
                                self.vof3[i,j] = 0.0
                        if a1 < 0:
                            if self.alpha1.cells[i-ii,j-jj] == 1: # Meglio > 0.5 ?
                                self.vof1[i,j] = -a2
                            else:
                                self.vof1[i,j] = 0.0

        for j in range(1, self.mesh.ny+1):
            for i in range(1, self.mesh.nx+1):
                a2 = vel[i,j]*dt/dxyz
                a1 = vel[i-ii,j-jj]*dt/dxyz
                self.alpha1.cells[i,j] = self.alpha1.cells[i,j] - (self.vof3[i,j] - self.vof1[i+ii,j+jj]) + \
                                (self.vof3[i-ii,j-jj] - self.vof1[i,j]) + self.cg[i,j]*(a2-a1)


    def advect(self):
        #self.reconstruct()
        self.computeFluxes()    # internal reconstruction only where needed
        self.clip()
        self.alpha1.updateBCs()


    def clip(self):
        for i in range(self.mesh.nx+2):
            for j in range(self.mesh.ny+2):
                if self.alpha1.cells[i,j] < self.EPSC:
                    self.alpha1.cells[i,j] = 0.
                if self.alpha1.cells[i,j] > 1.-self.EPSC:
                    self.alpha1.cells[i,j] = 1.


    def setEPSC(self, epsc=1.e-6):
        self.EPSC = epsc


    def getStencil3x3(self, i, j):
        stencil3x3 = np.zeros([3,3,3])
        # 2D, repeated along z
        indexlist = [-1, 0, 1]
        for ii in range(3):
            for jj in range(3):
                stencil3x3[ii,jj,:] = self.alpha1.cells[i+indexlist[ii],j+indexlist[jj]]

        return stencil3x3


