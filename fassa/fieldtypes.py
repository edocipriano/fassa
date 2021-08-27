import numpy as np

def ScalarField(mesh, fieldType="centered"):

    if fieldType == "centered":
        return np.zeros([mesh.nx+2, mesh.ny+2])

    elif fieldType == "uStaggered":
        return np.zeros([mesh.nx+1, mesh.ny+2])

    elif fieldType == "vStaggered":
        return np.zeros([mesh.nx+2, mesh.ny+1])

    else:
        NotImplemented


def ScalarNodesField(mesh):
    return np.zeros([mesh.nx+1, mesh.ny+1])


def VectorField(mesh, fieldType="centered"):

    if fieldType == "centered":
        return np.zeros([mesh.nx+2, mesh.ny+2, 3])

    elif fieldType == "uStaggered":
        return np.zeros([mesh.nx+1, mesh.ny+2, 3])

    elif fieldType == "vStaggered":
        return np.zeros([mesh.nx+2, mesh.ny+1, 3])

    else:
        NotImplemented


def SurfaceField(mesh):
    return np.zeros([mesh.nx+2, mesh.ny+2, 4])


def BoolField(mesh, value):
    return np.ones([mesh.nx+2, mesh.ny+2])*value


