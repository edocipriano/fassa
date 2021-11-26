import numpy as np
import fassa.fieldtypes

def ndS(mesh):
    """
    Compute product between face normal and face surface: n*dS
    with sign convention: positive when flux is exiting the face

                                  ^ +1*S
                                  |
                       -|---------x---------|-
                        |                   |
                        |                   |
                        |                   |
                        x-->    * i,j       x-->
                        |  -1*S             |  +1*S
                        |         ^ -1*S    |
                        |         |         |
                       -|---------x---------|-
    """
    return np.array([mesh.step, -mesh.step, mesh.step, -mesh.step])


def volumetricFlux(u, v):
    mesh = u.mesh
    phi = fassa.fieldtypes.SurfaceField(mesh)
    nds = ndS(mesh)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            phi[i,j,0] = v.cells[i,j]*nds[0]
            phi[i,j,1] = v.cells[i,j-1]*nds[1]
            phi[i,j,2] = u.cells[i,j]*nds[2]
            phi[i,j,3] = u.cells[i-1,j]*nds[3]

    return phi


def fieldFlux(u, v, fieldf):
    mesh = field.mesh
    phi = volumetricFlux(u, v)

    return phi*fieldf


def fieldFlux(phi, fieldf):
    return phi*fieldf


