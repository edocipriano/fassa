import numpy as np
import scipy.sparse as sparse

def laplacian(Coeff, field):
    mesh = field.mesh

    nxmat = mesh.nx**2
    nymat = mesh.ny**2

    Adense = np.zeros([nxmat, nymat])
    bdense = np.zeros([mesh.nx*mesh.ny,1])

    AEij = -2*Coeff/( (2*mesh.step)*(mesh.step) )
    AWij = -2*Coeff/( (2*mesh.step)*(mesh.step) )
    ANij = -2*Coeff/( (2*mesh.step)*(mesh.step) )
    ASij = -2*Coeff/( (2*mesh.step)*(mesh.step) )
    APij = -(AEij+AWij+ANij+ASij)

    for i in range(nxmat):
        for j in range(nymat):

            if i == j:
                Adense[i,j] = APij

            if i == j+1:
                Adense[i,j] = ANij

            if j == i+1:
                Adense[i,j] = ASij

            if i == mesh.nx+j:
                Adense[i,j] = AEij

            if j == mesh.ny+i:
                Adense[i,j] = AWij

    return sparse.csr_matrix(Adense), sparse.csr_matrix(bdense)


