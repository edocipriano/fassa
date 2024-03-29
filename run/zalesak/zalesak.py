import numpy as np
import os, math
import fassa
import fassa.vof

def main():

    if os.path.exists("results") == False:
        os.mkdir("results")

    # Create Mesh and set Boundaries
    mn = 128
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : nm,
                "ny" : nm}

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
    fassa.vof.utils.readVofi("../vofi/circle"+str(mesh.nx)+"/log", alpha1)
    alpha1.cells[math.ceil(61/128*mesh.nx):math.ceil(69/128*mesh.nx),math.ceil(50/128*mesh.nx):math.ceil(101/128*mesh.nx)] = 0.

    advector = fassa.vof.Advector(alpha1, u, v, time)

    # Stream function
    phi = fassa.Field("phi", None, mesh, 0.)

    # Additional fields to visualise:
    # Normals
    nx = fassa.Field("nx", None, mesh, 0.)
    ny = fassa.Field("ny", None, mesh, 0.)

    # Output file
    out = open("Output.out", "w")
    out.write("time[s] mass[kg]\n")
    out.close()
    liqvolume0 = volume_integral(alpha1)

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

        liqvolume = volume_integral(alpha1)
        out = open("Output.out", "a")
        out.write("%f %f\n"%(time.value, liqvolume/liqvolume0))
        out.close()

        # Assign fields to visualise
        nx.cells[:,:] = advector.n[:,:,0]
        ny.cells[:,:] = advector.n[:,:,1]

        if time.toWrite():
            fassa.nodeInterpolations(u, v, alpha1, nx, ny)
            filePath = "results/zalesak_"+str(time.nTimeSteps)+".vtk"
            fassa.writeVtk(filePath, mesh, u, v, nx, ny, alpha1, phi)

    print("\nEnd\n")


def volume_integral(alpha1):
    mesh = alpha1.mesh
    integral = 0.
    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            integral += alpha1.cells[i,j]*mesh.step**3

    return integral



if __name__=="__main__":
    main()


