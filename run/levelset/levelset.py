import numpy as np
import os
import fassa
import fassa.vof

def main():

    if os.path.exists("results") == False:
        os.mkdir("results")

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 128,
                "ny" : 128}

    mesh = fassa.hexCMesh(meshDict)
    mesh.setBasicBoundaries()
    mesh.setBasicCorners()

    # Vortex parameters
    T = 8
    CFL = 1

    # Create Time
    timeDict = {"totTime"       : T,
                "writeNow"      : True,
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
    fassa.vof.utils.readVofi("../vofi/circle"+str(mesh.nx)+"/log", alpha1)

    advector = fassa.vof.Advector(alpha1, u, v, time)

    # Stream function
    phi = fassa.Field("phi", None, mesh, 0.)

    # Levelset parameters
    deltaX = mesh.step
    gamma = 0.75*deltaX
    epsilon = 5.0*deltaX
    deltaTau = 0.1*deltaX
    dimChange = 1.0

    # Levelset field
    levelset  = fassa.Field("levelset",  alpha1Dict, mesh, 0.)
    levelset0 = fassa.Field("levelset0", alpha1Dict, mesh, 0.)

    # Heaviside function
    H = fassa.Field("H", None, mesh, 0.)

    # Dirac delta function
    delta = fassa.Field("delta", None, mesh, 0.)

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

        # Advect levelset function
        #advect_levelset(levelset, u, v, time, mesh)

        # Mapping levelset
        # map_levelset(levelset0, alpha1, gamma)
        levelset0.cells = (2.*alpha1.cells-1.)*gamma

        # Solve levelset function
        levelset.cells = levelset0.cells.copy()

        for corr in range(int(epsilon/deltaTau)):
            levelset.cells += np.sign(levelset0.cells)*(1.-mag_grad(levelset))*deltaTau
            levelset.updateBCs()

        # Update Dirac function
        for i in range(mesh.nx+1):
            for j in range(mesh.ny+1):

                if (np.abs(levelset.cells[i,j]) > epsilon):
                    levelset.cells[i,j] = 0.
                else:
                    levelset.cells[i,j] = 1./(2.*epsilon)*(1.+np.cos(np.pi*levelset.cells[i,j]/epsilon))

                if (levelset.cells[i,j] < -epsilon):
                    H.cells[i,j] = 0.
                elif (epsilon < levelset.cells[i,j]):
                    H.cells[i,j] = 1.
                else:
                    H.cells[i,j] = 1./2.*(1.+levelset.cells[i,j]/epsilon+np.sin(np.pi*levelset.cells[i,j]/epsilon)/np.pi)


        # Assign fields to visualise
        nx.cells[:,:] = advector.n[:,:,0]
        ny.cells[:,:] = advector.n[:,:,1]

        if time.toWrite():
            fassa.nodeInterpolations(u, v, alpha1, nx, ny, levelset, H)
            filePath = "results/levelset_"+str(time.nTimeSteps)+".vtk"
            fassa.writeVtk(filePath, mesh, u, v, nx, ny, alpha1, phi, levelset, H)

    print("\nEnd\n")


def set_circle_levelset(levelset, centre, radius):

    mesh = levelset.mesh

    xcenter = centre[0]
    ycenter = centre[1]

    for i in range(mesh.nx+1):
        for j in range(mesh.ny+1):
            #phi.nodes[i,j] = 1./3.14*np.cos(3.14*time.value/T)*(np.sin(3.14*mesh.x[i])**2)*(np.sin(3.14*mesh.y[j])**2)
            phi.nodes[i,j] = -1.5*np.sin(np.pi*time.value/T)*np.sin((mesh.x[i] + 0.5)*np.pi)*np.sin((mesh.y[j]+0.5)*np.pi)/np.pi


def map_levelset(levelset, alpha1, gamma):
    mesh = levelset.mesh
    levelset.cells = (2.*alpha1.cells - 1.)*gamma


def advect_levelset(levelset, u, v, time, mesh):

    lso = levelset.cells.copy()

    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny):

            ue = u.cells[i,j]
            uw = u.cells[i-1,j]
            vn = v.cells[i,j]
            vs = v.cells[i,j-1]

            lse = 0.5*(levelset.cells[i+1,j] + levelset.cells[i,j])
            lsw = 0.5*(levelset.cells[i,j] + levelset.cells[i-1,j])
            lsn = 0.5*(levelset.cells[i,j+1] + levelset.cells[i,j])
            lss = 0.5*(levelset.cells[i,j] + levelset.cells[i,j-1])

            Aij = mesh.step*(ue*lse - uw*lsw + vn*lsn - vs*lss)

            levelset.cells[i,j] = lso[i,j] + time.deltaT/(mesh.step**2)*(-Aij)


def grad(field):
    mesh = field.mesh
    grad = fassa.fieldtypes.VectorField(mesh)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            grad[i,j,0] = (field[i+1,j]-field[i-1,j])/(2.*mesh.step)
            grad[i,j,1] = (field[i,j+1]-field[i,j-1])/(2.*mesh.step)
            grad[i,j,2] = 0.

    return grad


def mag_grad(field):
    mesh = field.mesh
    grad    = fassa.fieldtypes.VectorField(mesh)
    maggrad = fassa.fieldtypes.ScalarField(mesh)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            grad[i,j,0] = (field[i+1,j]-field[i-1,j])/(2.*mesh.step)
            grad[i,j,1] = (field[i,j+1]-field[i,j-1])/(2.*mesh.step)
            grad[i,j,2] = 0.
            maggrad[i,j] = (grad[i,j,0]**2 + grad[i,j,1]**2)**0.5

    return maggrad


def mag(field, mesh):
    mag = fassa.fieldtypes.ScalarField(mesh)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            mag[i,j] = (field[i,j,0]**2 + field[i,j,1]**2)**0.5

    return mag


def velocityFromStreamfunction(phi, u, v):
    """
    Compute staggered  velocity fields from the stream function
    """
    mesh = u.mesh

    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny):
            u.cells[i,j] = -(phi.nodes[i+1,j] - phi.nodes[i,j])/(mesh.x[i+1] - mesh.x[i])

    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny):
            v.cells[i,j] = (phi.nodes[i,j+1] - phi.nodes[i,j])/(mesh.y[j+1] - mesh.y[j])


if __name__=="__main__":
    main()


