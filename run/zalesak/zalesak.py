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
                "nx" : 128,
                "ny" : 128}

    mesh = fassa.hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Create Time
    timeDict = {"totTime"       : 2*3.14,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-2,
                "writeSteps"    : 20}

    time = fassa.Time(timeDict, mesh)

    # Velocity fields
    u = fassa.Field("u", None, mesh, 0., fieldType="uStaggered")
    v = fassa.Field("v", None, mesh, 0, fieldType="vStaggered")

    # VOF field
    alpha1Dict = {"N" : ("zeroGradient", None),
                  "S" : ("zeroGradient", None),
                  "E" : ("zeroGradient", None),
                  "W" : ("zeroGradient", None)}

    alpha1 = fassa.Field("alpha1", alpha1Dict, mesh, 0.)

    # Read log from VOFI initialization
    fassa.vof.utils.readVofi("../vofi/circle128/log", alpha1)
    alpha1.cells[61:69,50:101] = 0.

    advector = fassa.vof.Advector(alpha1, u, v, time)

    # Stream function
    phi = fassa.Field("phi", None, mesh, 0.)

    # Additional fields to visualise:
    # Normals
    nx = fassa.Field("nx", None, mesh, 0.)
    ny = fassa.Field("ny", None, mesh, 0.)

    # Solution loop
    while time.toRun():

        for i in range(1, mesh.nx):
            for j in range(1, mesh.ny):
                u.cells[i,j] = 0.5 - mesh.Cy[j]
                v.cells[i,j] = mesh.Cx[i] - 0.5

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
            filePath = "results/zalesak_"+str(time.nTimeSteps)+".vtk"
            fassa.utils.write.writeVtk(filePath, mesh, u, v, nx, ny, alpha1, phi)

    print("\nEnd\n")


if __name__=="__main__":
    main()


