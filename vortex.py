import numpy as np
from hexCMesh import hexCMesh
from Time import Time
import timeStepSelector
import Field
from Field import Field, onesField, nodeInterpolations
import BCs
import utils.write
from solve.explicit.ProjectionMethod import ProjectionMethod
import fluxes
from solve.explicit.divergence import div as expdiv
from vof.Advector import Advector
import vof.utils

def main():

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 64,
                "ny" : 64}

    mesh = hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Vortex parameters
    T = 2
    CFL = 1

    # Create Time
    timeDict = {"totTime"       : 4,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-2,
                "writeSteps"    : 10}

    time = Time(timeDict, mesh)

    # Velocity fields
    u = Field("u", None, mesh, 0., fieldType="uStaggered")
    v = Field("v", None, mesh, 0, fieldType="vStaggered")

    # VOF field
    alpha1Dict = {"N" : ("zeroGradient", None),
                  "S" : ("zeroGradient", None),
                  "E" : ("zeroGradient", None),
                  "W" : ("zeroGradient", None)}

    alpha1 = Field("alpha1", alpha1Dict, mesh, 0.)

    # Read log from VOFI initialization
    vof.utils.readVofi("vofi/circle/log", alpha1)

    advector = Advector(alpha1, u, v, time)

    # Stream function
    phi = Field("phi", None, mesh, 0.)

    # Additional fields to visualise:

    # Normals
    nx = Field("nx", None, mesh, 0.)
    ny = Field("ny", None, mesh, 0.)

    # Fluxes
    fluxRight = Field("fluxRight", None, mesh, 0.)
    fluxLeft  = Field("fluxLeft",  None, mesh, 0.)

    SMALL = 1.e-12

    # Solution loop
    while time.toRun():

        # Update stream function
        for i in range(mesh.nx+1):
            for j in range(mesh.ny+1):
                phi.nodes[i,j] = 1./3.14*np.cos(3.14*time.value/T)*(np.sin(3.14*mesh.x[i])**2)*(np.sin(3.14*mesh.y[j])**2)

        # Compute velocity from the stream function
        velocityFromStreamfunction(phi, u, v)

        dtMin = timeStepSelector.alphaCo(u, v)
        time.setDeltaT(dtMin)
        print("Minimum deltaT for stabilization = %.10f"% dtMin)
        print("Current deltaT = %.10f"% time.deltaT)

        time.print()

        advector.advect()

        # Assign fields to visualise

        nx.cells[:,:] = advector.n[:,:,0]
        ny.cells[:,:] = advector.n[:,:,1]
        fluxRight.cells[:] = advector.vof3[:]
        fluxLeft.cells[:]  = advector.vof1[:]


        if time.toWrite():
            nodeInterpolations(u, v, alpha1, nx, ny, fluxRight, fluxLeft)
            filePath = "results/vortex_"+str(time.nTimeSteps)+".vtk"
            utils.write.writeVtk(filePath, mesh, u, v, nx, ny, fluxRight, fluxLeft, alpha1, phi)

    print("\nEnd\n")


def velocityFromStreamfunction(phi, u, v):
    """
    Compute staggered  velocity fields from the stream function
    """
    mesh = u.mesh

    # u = d(phi)/dy
    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny):
             u.cells[i,j] = (phi.nodes[i,j] - phi.nodes[i,j-1])/(mesh.y[j] - mesh.y[j-1])

    # v = -d(phi)/dx
    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny+1):
            v.cells[i,j] = -(phi.nodes[i,j] - phi.nodes[i-1,j])/(mesh.x[i] - mesh.x[i-1])


if __name__=="__main__":
    main()


