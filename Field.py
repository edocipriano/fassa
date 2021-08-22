import numpy as np
from hexCMesh import hexCMesh
import BCs
import FieldTypes
import utils.write
import utils.initFields

class Field:

    def __init__(self, name, fieldDict, mesh, initValue, fieldType="centered"):
        self.name = name
        self.fieldDict = fieldDict
        self.fieldType = fieldType
        self.initValue = initValue
        self.mesh = mesh

        self.centered   = False
        self.uStaggered = False
        self.vStaggered = False

        if (self.fieldType == "centered"):
            self.centered = True
            self._setCellCenteredField()

        elif (self.fieldType == "uStaggered"):
            self.uStaggered = True
            self._setUStaggeredField()

        elif (self.fieldType == "vStaggered"):
            self.vStaggered = True
            self._setVStaggeredField()

        else:
            NotImplemented

        self.faces = FieldTypes.SurfaceField(mesh)
        self.nodes = FieldTypes.ScalarNodesField(mesh)


    def _setCellCenteredField(self):
        self.cells = np.ones([self.mesh.nx+2, self.mesh.ny+2]) * self.initValue


    def _setUStaggeredField(self):
        self.cells = np.ones([self.mesh.nx+1, self.mesh.ny+2]) * self.initValue


    def _setVStaggeredField(self):
        self.cells = np.ones([self.mesh.nx+2, self.mesh.ny+1]) * self.initValue


    def nodeInterpolation(self):
        """
            Reconstruct a field from cells-centered to vertex-centered.
            Useful for post-processing puroposes
        """
        nx = self.mesh.nx
        ny = self.mesh.ny

        if self.centered == True:
            self.nodes[0:nx+1,0:ny+1] = 0.25*(self.cells[0:nx+1,0:ny+1] + \
                                              self.cells[0:nx+1,1:ny+2] + \
                                              self.cells[1:nx+2,0:ny+1] + \
                                              self.cells[1:nx+2,1:ny+2])

        elif self.uStaggered == True:
            self.nodes[0:nx+1,0:ny+1] = 0.5*(self.cells[0:nx+1,1:ny+2] + self.cells[0:nx+1,0:ny+1])

        elif self.vStaggered == True:
            self.nodes[0:nx+1,0:ny+1] = 0.5*(self.cells[1:nx+2,0:ny+1] + self.cells[0:nx+1,0:ny+1])

        else:
            NotImplemented


    def faceInterpolation(self, i, j):
        """
            Return a dict with 4 components corresponding to the value of
            the field on the different faces of a cells ordered as: N, S, E, W

                                N
               -|---------x---------|-
             W  |                   |      Face values are computed from the
                |                   |      cells centers to the surface centers x
                |                   |      with linear interpolation with the
                x        i,j        x      neightbouring cellss
                |                   |
                |                   |
                |                   |  E
               -|---------x---------|-
                  S

            Usage: north = alpha1.faceInterpolation(i,j)["N"]
        """
        fn = 0.5*(self.cells[i,j] + self.cells[i,j+1])
        fs = 0.5*(self.cells[i,j] + self.cells[i,j-1])
        fe = 0.5*(self.cells[i,j] + self.cells[i+1,j])
        fw = 0.5*(self.cells[i,j] + self.cells[i-1,j])

        #return np.array([fn, fs, fe, fw])
        return {"N" : fn,
                "S" : fs,
                "E" : fe,
                "W" : fw}


    def interpolate(self):

        for i in range(1, self.mesh.nx+1):
            for j in range(1, self.mesh.ny+1):
                self.faces[i,j,0] = 0.5*(self.cells[i,j] + self.cells[i,j+1])
                self.faces[i,j,1] = 0.5*(self.cells[i,j] + self.cells[i,j-1])
                self.faces[i,j,2] = 0.5*(self.cells[i,j] + self.cells[i+1,j])
                self.faces[i,j,3] = 0.5*(self.cells[i,j] + self.cells[i-1,j])


    def updateBCs(self):
        """
            Change ghost cells values in order to set the Boundary Conditions

            mesh.boundaries = {"Boundary1" : (sliceOrig[0,:], sliceNbr[1,:], "alignedAxis")}
        """
        # if the field has BCs
        if self.fieldDict is not None:
            for boundary in self.mesh.boundaries:

                if self.fieldDict[boundary][0] == "fixedValue":
                    BCs.fixedValue(self, boundary)

                if self.fieldDict[boundary][0] == "zeroGradient":
                    BCs.zeroGradient(self, boundary)

                else:
                    NotImplemented
        else:
            raise Exception("Trying to update BCs for field "+self.name+" without BCs")



def onesField(mesh):
    return Field("ones", None, mesh, 1.)

def zerosField(mesh):
    return Field("zeros", None, mesh, 0.)

def nodeInterpolations(*fields):
    for field in fields:
        field.nodeInterpolation()

