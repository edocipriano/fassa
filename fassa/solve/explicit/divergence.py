from fassa.Field import Field
import fassa.fluxes as fluxes

def div(phi, scalar, interpolate=True):

    mesh = scalar.mesh

    if interpolate:
        scalar.interpolate()
        scalarf = scalar.faces
    else:
        scalarf = scalar

    scalarFlux = fluxes.fieldFlux(phi, scalarf)

    divergence = Field("divergence", None, mesh, 0.)

    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny+1):
            divergence.cells[i,j] = scalarFlux[i,j,:].sum()

    return divergence

