import numpy as np
import os
import fassa
import fassa.vof
import fassa.utils

def main():

    if os.path.exists("results") == False:
        os.system("mkdir results")

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 64,
                "ny" : 64}

    mesh = fassa.hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Vortex parameters
    T = 8
    CFL = 1

    # Create Time
    timeDict = {"totTime"       : T,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-2,
                "writeSteps"    : 10}

    time = fassa.Time(timeDict, mesh)

    # Velocity fields
    u = fassa.Field("u", None, mesh, 0., fieldType="uStaggered")
    v = fassa.Field("v", None, mesh, 0., fieldType="vStaggered")

    # VOF field
    alpha1Dict = {"N" : ("zeroGradient", None),
                  "S" : ("zeroGradient", None),
                  "E" : ("zeroGradient", None),
                  "W" : ("zeroGradient", None)}

    alpha1 = fassa.Field("alpha1", alpha1Dict, mesh, 0.)

    # Read log from VOFI initialization
    fassa.vof.utils.readVofi("../vofi/circle64/log", alpha1)

    advector = fassa.vof.Advector(alpha1, u, v, time)

    # Stream function
    phi = fassa.Field("phi", None, mesh, 0.)

    # Additional fields to visualise:

    # Normals
    nx = fassa.Field("nx", None, mesh, 0.)
    ny = fassa.Field("ny", None, mesh, 0.)

    # Solution loop
    while time.toRun():

        # Update stream function
        for i in range(mesh.nx+1):
            for j in range(mesh.ny+1):
                phi.nodes[i,j] = 1./3.14*np.cos(3.14*time.value/T)*(np.sin(3.14*mesh.x[i])**2)*(np.sin(3.14*mesh.y[j])**2)

        # Compute velocity from the stream function
        velocityFromStreamfunction(phi, u, v)

        dtMin = fassa.timeStepSelector.alphaCo(u, v)
        time.setDeltaT(dtMin)
        print("\nMinimum deltaT for stabilization = %.10f"% dtMin)
        print("Current deltaT = %.10f"% time.deltaT)

        time.print()

        advector.advect()

        # Assign fields to visualise
        nx.cells[:,:] = advector.n[:,:,0]
        ny.cells[:,:] = advector.n[:,:,1]

        if time.toWrite():
            fassa.nodeInterpolations(u, v, alpha1, nx, ny)
            filePath = "results/vortex_"+str(time.nTimeSteps)+".vtk"
            fassa.utils.write.writeVtk(filePath, mesh, u, v, nx, ny, alpha1, phi)

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


