import numpy as np

def fieldAverage(alpha1, field, value1, value2):
    mesh = alpha1.mesh

    for i in range(mesh.nx+2):
        for j in range(mesh.ny+2):
            field.cells[i,j] = value1*alpha1.cells[i,j] + value2*(1.-alpha1.cells[i,j])


def readVofi(filename, nx, ny):

    # Cell centered field without ghost cells
    cc = np.zeros([nx, ny])
    nplines = np.zeros(nx*ny)

    f = open(filename, "r")
    lines = f.readlines()

    for count, line in enumerate(lines):
        nplines[count] = float(line)

    count = 0
    for i in range(nx):
        for j in range(ny):
            cc[i,j] = nplines[count]
            count += 1

    return cc


def readVofi(filename, field):

    nx = field.mesh.nx
    ny = field.mesh.ny

    # Cell centered field without ghost cells
    nplines = np.zeros(nx*ny)

    f = open(filename, "r")
    lines = f.readlines()

    for count, line in enumerate(lines):
        nplines[count] = float(line)

    count = 0
    for i in range(nx):
        for j in range(ny):
            field.cells[i+1,j+1] = nplines[count]
            count += 1

