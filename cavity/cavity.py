import numpy as np
from hexCMesh import hexCMesh
from Time import Time
from Field import Field, onesField, nodeInterpolations
import BCs
import utils.write
from solve.explicit.ProjectionMethod import ProjectionMethod
import fluxes
from solve.explicit.divergence import div as expdiv
import timeStepSelector

def main():

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 50,
                "ny" : 50}

    mesh = hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Create Time
    timeDict = {"totTime"       : 5,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-4,
                "writeSteps"    : 50}

    time = Time(timeDict, mesh)

    # Velocity fields
    uDict = {"N" : ("fixedValue", 1.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    u = Field("u", uDict, mesh, 0., fieldType="uStaggered")

    vDict = {"N" : ("fixedValue", 0.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    v = Field("v", vDict, mesh, 0, fieldType="vStaggered")

    # Pressure field: using the explicit Projection Method
    # boundary conditions for pressure are not explicitely required
    pDict = {"N" : ("fixedValue", 0.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    p = Field("p", pDict, mesh, 0)

    # Physical properties
    propertiesDict = {"nu" : 1.e-2}

    nu = Field("nu", None, mesh, propertiesDict["nu"])

    # Create Pressure-Velocity coupling object
    projectionDict = {"maxIter"  : 500,
                      "maxError" : 1.e-4,
                      "SORCoeff" : 1.45}

    projection = ProjectionMethod(u, v, p, mesh, time, projectionDict)

    # Solution loop
    while time.toRun():

        BCs.updateBCs(u, v)

        dtMin = timeStepSelector.Co(u, v)
        time.setDeltaT(dtMin)
        print("Minimum deltaT for stabilization = %.10f"% dtMin)
        print("Current deltaT = %.10f"% time.deltaT)

        projection.predictorStep(nu)
        projection.projectionStep()
        projection.correctVelocity()

        projection.print()

        phi = fluxes.volumetricFlux(u,v)
        div = expdiv(phi, onesField(mesh))

        if time.toWrite():
            nodeInterpolations(u, v, p, div)
            filePath = "results/cavity_"+str(time.nTimeSteps)+".vtk"
            utils.write.writeVtk(filePath, mesh, u, v, p, div)

    print("\nEnd\n")

if __name__=="__main__":
    main()


