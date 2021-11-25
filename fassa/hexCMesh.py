import numpy as np

class hexCMesh:

    def __init__(self, meshDict):

        # Get data from dictionary
        self.Lx = meshDict["Lx"]
        self.Ly = meshDict["Ly"]
        self.nx = meshDict["nx"]
        self.ny = meshDict["ny"]
        self.boundaries = {}
        self.corners = {}

        self._constructMesh()
        self._constructCellCenterCoordinates()

    def _constructMesh(self):

        # Compute grid step (cuboid cells)
        self.step = self.Lx/self.nx

        # Compute node positions
        self.x = np.linspace(0, self.Lx, self.nx+1)
        self.y = np.linspace(0, self.Ly, self.ny+1)
        self.z = np.array([0]) # 3rd dim NotImplemented

    def _constructCellCenterCoordinates(self):
        # 2 dimensions
        self.Cx = np.zeros(self.nx)
        self.Cy = np.zeros(self.ny)

        xcoord = self.step/2.
        for i in range(self.nx):
            self.Cx[i] = xcoord
            xcoord += self.step

        ycoord = self.step/2.
        for i in range(self.ny):
            self.Cy[i] = ycoord
            ycoord += self.step

    def setBoundary(self, name, ghostLayer, internalLayer, parallelAxis):
        """
            Define a boundary given:
            @param name: Name of the boundary
            @param ghostLayer: slice with ghost cells belonging to this boundary
            @param internalLayer: slice with layer of cells near the boundary ghost layer
            @param parallelAxis: along which axis the boundary is parallel, in 2D "x" or "y"

            Slice can be conveniently constructed using np.index_exp or np.s_
            since their syntax is more intuitive than the standard slice
        """
        self.boundaries[name] = (ghostLayer, internalLayer, parallelAxis)

    def setBasicBoundaries(self):
        """
            Set basic North, South, East, West boundaries

            @output dict{"N": (sliceNghost, sliceNnbr, "x"),
                         "S": (sliceSghost, sliceSnbr, "y"),
                         "E": (sliceEghost, sliceEnbr, "y"),
                         "W": (sliceWghost, sliceWnbr, "x")}

            From a field:
             -> the boundary ghost cells layer is called with:
                field.cell[self.mesh.boundaries[bName][0]
             -> the boundary neighbout to ghost layer is called with:
                field.cell[self.mesh.boundaries[bName][1]

        """
        self.setBoundary("N", np.index_exp[:,-1], np.index_exp[:,-2], "x")
        self.setBoundary("S", np.index_exp[:,0],  np.index_exp[:,1],  "x")
        self.setBoundary("E", np.index_exp[-1,:], np.index_exp[-2,:], "y")
        self.setBoundary("W", np.index_exp[0,:],  np.index_exp[1,:],  "y")


    def setCorner(self, name, ghostCell, internalCell):
        self.corners[name] = (ghostCell, internalCell)


    def setBasicCorners(self):
        self.setCorner("SE", np.index_exp[0,-1],  np.index_exp[1,-2])
        self.setCorner("SW", np.index_exp[0,0],   np.index_exp[1,1])
        self.setCorner("NE", np.index_exp[-1,-1], np.index_exp[-2,-2])
        self.setCorner("NW", np.index_exp[0,-1],  np.index_exp[1,-2])


def basicMesh(nx=50,ny=50):

    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : nx,
                "ny" : ny}

    return hexCMesh(meshDict)

