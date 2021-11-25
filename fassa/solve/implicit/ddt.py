import numpy as np
import scipy.sparse as sparse

def ddt(time, field):

    mesh = field.mesh
    dt = time.deltaT

    nxmat = mesh.nx**2
    nymat = mesh.ny**2

    Adense = np.zeros([nxmat, nymat])
    bdense = np.zeros([mesh.nx*mesh.ny,1])

    for i in range(nxmat):
        for j in range(nymat):

            if i == j:
                Adense[i,j] = 1./dt

    for i in range(mesh.nx):
        for j in range(mesh.ny):
            l = (i-1)*mesh.ny + j

            bdense[l,0] = 1./dt*field.cells[i,j]

    return sparse.csr_matrix(Adense), sparse.csr_matrix(bdense)


