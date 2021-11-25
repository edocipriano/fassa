from fassa.Field import Field
import fassa.fluxes as fluxes

def mag(f):
    """
    Compute the absolute value of a Field f
    Return a fieldtype.ScalarField
    """
    return np.abs(f)


def grad(f, mesh):
    """
    Compute the finite volume gradient of a scalar field f
    Return a vector field and a scalar field with the gradient
    of f and with the norm of the gradient f
    """
    gradvec = fassa.fieldtypes.VectorField(mesh)
    gradmag = fassa.fieldtypes.ScalarField(mesh)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            gradvec[i,j,0] = (f[i+1,j]-f[i-1,j])/(2.*mesh.step)
            gradvec[i,j,1] = (f[i,j+1]-f[i,j-1])/(2.*mesh.step)
            gradvec[i,j,2] = 0. # 2D
            gradmag[i,j] = (gradvec[i,j,0]**2 + gradvec[i,j,1]**2)**0.5

    return gradvec, gradmag

