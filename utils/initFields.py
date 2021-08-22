import numpy as np

def box(field, sw, ne, shading=True, internalValue=1, externalValue=0, shadedValue=0.5):
    """
        Initialize a field with given value
        inside a box with:
        @param internalValue: value of field.cells inside the box
        @param externalValue: value of field.cells outside the box
        @param sw: South-West cell of the box
        @param ne: North-East cell of the box
    """
    mesh = field.mesh
    xSw, ySw, xNe, yNe = sw[0], sw[1], ne[0], ne[1]

    field.cells[:,:]              = externalValue
    field.cells[xSw:xNe, ySw:yNe] = internalValue

    if shading == True:
        field.cells[xSw:xNe, yNe]     = shadedValue
        field.cells[xSw, ySw:yNe]     = shadedValue
        field.cells[xNe, ySw:yNe]     = shadedValue
        field.cells[xSw:xNe, ySw]     = shadedValue
        field.cells[xNe, yNe]         = shadedValue/4
        field.cells[xSw, ySw]         = shadedValue/4
        field.cells[xSw, yNe]         = shadedValue/4
        field.cells[xNe, ySw]         = shadedValue/4

