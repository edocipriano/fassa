import numpy as np
import os
import fassa
import fassa.vof
import fassa.BCs

def main():

    if os.path.exists("results") == False:
        os.mkdir("results")

    # Create Mesh and set Boundaries
    mn = 128
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : mn,
                "ny" : mn}

    mesh = fassa.hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Vortex parameters
    T = 15
    CFL = 1.

    # Create Time
    timeDict = {"totTime"       : T,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-2,
                "writeSteps"    : 20}

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
    fassa.vof.utils.readVofi("../vofi/circle"+str(mesh.nx)+"/log", alpha1)

    advector = fassa.vof.Advector(alpha1, u, v, time)
    advector.setEPSC(1.e-12)

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

        # Update stream function
        pi = 3.1415926
        for i in range(mesh.nx+1):
            for j in range(mesh.ny+1):
              phi.nodes[i,j] = - 1.5*np.sin(2.*pi*time.value/T)*np.sin((mesh.x[i] + 0.5)*pi)*np.sin((mesh.y[j] + 0.5)*pi)/pi;

        # Compute velocity from the stream function
        velocityFromStreamfunction(phi, u, v, time, T)

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
            filePath = "results/vortex_"+str(time.nTimeSteps)+".vtk"
            fassa.writeVtk(filePath, mesh, u, v, nx, ny, alpha1, phi)

    print("\nEnd\n")

def volume_integral(alpha1):
    mesh = alpha1.mesh
    integral = 0.
    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            integral += alpha1.cells[i,j]*mesh.step**3

    return integral


def velocityFromStreamfunction(phi, u, v, time, T):
    """
    Compute staggered  velocity fields from the stream function
    """
    mesh = u.mesh

    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny+1):
            u.cells[i,j] = -(phi.nodes[i+1,j] - phi.nodes[i,j])/(mesh.x[i+1] - mesh.x[i])
            #u.cells[i,j] =  2.*np.cos(3.14*time.value/T)*np.sin(3.14*mesh.x[i])**2*np.sin(3.14*mesh.y[j])*np.cos(3.14*mesh.y[j])

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny):
            v.cells[i,j] = (phi.nodes[i,j+1] - phi.nodes[i,j])/(mesh.y[j+1] - mesh.y[j])
            #v.cells[i,j] = -2.*np.cos(3.14*time.value/T)*np.sin(3.14*mesh.y[j])**2*np.sin(3.14*mesh.x[i])*np.cos(3.14*mesh.x[i])



if __name__=="__main__":
    main()


