def updateBCs(*fields):
    for field in fields:
        field.updateBCs()


def fixedValue(field, boundary):
    """
    Fixed Value Boundary Condition: use linear interpolation to fix the value
    on the ghost cell which respects the fixedValue on the boundary face.

    Different situations are implemented for staggered and cell-centered fields:
    If the field is uStaggered (vx) and the boundary is aligned with the y axis,
    then the cells value of vx is already the face value and the assignment of
    the boundary condition is straight forward.
    """
    if field.uStaggered == True and field.mesh.boundaries[boundary][2] == "y":
        field.cells[field.mesh.boundaries[boundary][0]] = field.fieldDict[boundary][1]

    elif field.vStaggered == True and field.mesh.boundaries[boundary][2] == "x":
        field.cells[field.mesh.boundaries[boundary][0]] = field.fieldDict[boundary][1]

    else:
        field.cells[field.mesh.boundaries[boundary][0]] = 2*field.fieldDict[boundary][1] - \
                                                          field.cells[field.mesh.boundaries[boundary][1]]


def zeroGradient(field, boundary):
    """
    Implementation of zeroGradient Boundary Condition.
    It's the same for staggered and cell-centered fields.
    """
    field.cells[field.mesh.boundaries[boundary][0]] = field.cells[field.mesh.boundaries[boundary][1]]

