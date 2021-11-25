import numpy as np

SMALL = 1.e-12

def Co(u, v):
    """
    Stability condition on the velocity
    """
    mesh = u.mesh

    # Find maximum velocity
    umax = abs(u.cells.max())
    vmax = abs(v.cells.max())

    sumUoverStep = (umax+ vmax)/mesh.step + SMALL

    return 1.0/sumUoverStep


def alphaCo(u, v):
    """
    Stability condition on the velocity of interface
    advection based on the Weymouth-Yue advection scheme
    """
    mesh = u.mesh

    # Find maximum velocity
    umax = abs(u.cells.max())
    vmax = abs(v.cells.max())

    sumUoverStep = (umax+ vmax)/mesh.step + SMALL

    return 0.5/sumUoverStep
