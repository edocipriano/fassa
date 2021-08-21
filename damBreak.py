import numpy as np
from hexCMesh import hexCMesh
from Time import Time
import timeStepSelector
from Field import Field, onesField
import FieldTypes
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
                "nx" : 50,
                "ny" : 50}

    mesh = hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Create Time
    timeDict = {"totTime"       : 0.8,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 0.05,
                "deltaT"        : 1.e-2,
                "writeSteps"    : 1}

    time = Time(timeDict, mesh)

    # Velocity fields
    uDict = {"N" : ("zeroGradient", 0.),
             "S" : ("zeroGradient", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    u = Field("u", uDict, mesh, 0., fieldType="uStaggered")

    vDict = {"N" : ("zeroGradient", 0.),
             "S" : ("zeroGradient", 0.),
             "E" : ("fixedValue", 0),
             "W" : ("fixedValue", 0)}

    v = Field("v", vDict, mesh, 0, fieldType="vStaggered")

    # Pressure field: using the explicit Projection Method
    # boundary conditions for pressure are not explicitely required
    pDict = {"N" : ("fixedValue", 0.),
             "S" : ("fixedValue", 0.),
             "E" : ("fixedValue", 0.),
             "W" : ("fixedValue", 0.)}

    p = Field("p", pDict, mesh, 0.)

    # VOF field
    alpha1Dict = {"N" : ("zeroGradient", None),
                  "S" : ("zeroGradient", None),
                  "E" : ("zeroGradient", None),
                  "W" : ("zeroGradient", None)}

    alpha1 = Field("alpha1", alpha1Dict, mesh, 0.)
    utils.initFields.box(alpha1, [15,15], [35,35], shading=False)

    advector = Advector(alpha1, u, v, time)

    # Physical properties
    rhoL, rhoG = 10, 5
    muL, muG = 1.e-2, 1.e-5
    nuL, nuG = 1.e-2, 1e-2

    nu = Field("nu", None, mesh, 1.e-2)
    vof.utils.fieldAverage(alpha1, nu, nuL, nuG)

    # Create Pressure-Velocity coupling object
    projectionDict = {"maxIter"  : 5000,
                      "maxError" : 1.e-4,
                      "SORCoeff" : 1.45}

    projection = ProjectionMethod(u, v, p, mesh, time, projectionDict)

    # Additional fields to write
    nx = Field("nx", None, mesh, 0.)
    ny = Field("ny", None, mesh, 0.)
    fluxRight = Field("fluxRight", None, mesh, 0.)
    fluxLeft  = Field("fluxLeft",  None, mesh, 0.)

    # Solution loop
    while time.toRun():

        for j in range(mesh.ny+2):
            u.cells[1:,j] = mesh.Cx[:]/2.

        for i in range(mesh.nx+2):
            v.cells[i,1:] = mesh.Cy[:]/2.

        dtMin = timeStepSelector.Co(u, v)
        time.setDeltaT(dtMin)
        print("Minimum deltaT for stabilization = %.10f"% dtMin)
        print("Current deltaT = %.10f"% time.deltaT)

        """
        BCs.updateBCs(u, v)

        projection.predictorStep(nu)
        projection.projectionStep()
        projection.correctVelocity()
        """
        projection.print()

        advector.advect()

        nx.cells[:,:] = advector.n[:,:,0]
        ny.cells[:,:] = advector.n[:,:,1]
        
        #fluxRight.cells[:] = advector.fluxRight[:]
        #fluxLeft.cells[:]  = advector.fluxLeft[:]

        if time.toWrite():
            filePath = "results/damBreak_"+str(time.nTimeSteps)+".vtk"
            utils.write.writeVtk(filePath, mesh, u, v, p, alpha1, nu, nx, ny, fluxRight, fluxLeft)

    print("\nEnd\n")

if __name__=="__main__":
    main()


