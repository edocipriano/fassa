import numpy as np
import os
import fassa

def main():

    if os.path.exists("results") == False:
        os.mkdir("results")

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 50,
                "ny" : 50}

    mesh = fassa.hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Create Time
    timeDict = {"totTime"       : 5,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-4,
                "writeSteps"    : 50}

    time = fassa.Time(timeDict, mesh)

    # Velocity fields
    uDict = {"N" : ("fixedValue", 1.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    u = fassa.Field("u", uDict, mesh, 0., fieldType="uStaggered")

    vDict = {"N" : ("fixedValue", 0.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    v = fassa.Field("v", vDict, mesh, 0, fieldType="vStaggered")

    # Pressure field: using the explicit Projection Method
    # boundary conditions for pressure are not explicitely required
    pDict = {"N" : ("fixedValue", 0.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    p = fassa.Field("p", pDict, mesh, 0)

    # Physical properties
    nu = fassa.Field("nu", None, mesh, 1.e-2)

    # Create Pressure-Velocity coupling object
    projectionDict = {"maxIter"  : 500,
                      "maxError" : 1.e-4,
                      "SORCoeff" : 1.45}

    projection = fassa.ProjectionMethod(u, v, p, mesh, time, projectionDict)

    # Solution loop
    while time.toRun():

        fassa.BCs.updateBCs(u, v)

        dtMin = fassa.timeStepSelector.Co(u, v)
        time.setDeltaT(dtMin)
        print("\nMinimum deltaT for stabilization = %.10f"% dtMin)
        print("Current deltaT = %.10f"% time.deltaT)

        projection.predictorStep(nu)
        projection.projectionStep()
        projection.correctVelocity()

        projection.print()

        phi = fassa.volumetricFlux(u,v)
        div = fassa.div(phi, fassa.onesField(mesh))

        if time.toWrite():
            fassa.nodeInterpolations(u, v, p, div)
            filePath = "results/cavity_"+str(time.nTimeSteps)+".vtk"
            fassa.writeVtk(filePath, mesh, u, v, p, div)

    print("\nEnd\n")

if __name__=="__main__":
    main()


