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


# NOTE not working
def circle(field):

    mesh = field.mesh

    field.cells *= 0.0
    R = 0.25

    markerAlpha = field.cells.copy()

    for i in range(mesh.nx+2):
        for j in range(mesh.ny+2):

            deltaY = mesh.step
            x_c = mesh.Cx[i]
            y_c = mesh.Cy[j]

            x1 = x_c - deltaY/2.
            x2 = x_c + deltaY/2.
            y1 = y_c - deltaY/2.
            y2 = y_c + deltaY/2.

            v1_in = False
            v2_in = False
            v3_in = False
            v4_in = False

            if ((x1*x1+y1*y1)**0.5 - R) <= 0.:
                v1_in = True

            if ((x2*x2+y1*y1)**0.5 - R) <= 0.:
                v2_in = True

            if ((x2*x2+y2*y2)**0.5 - R) <= 0.:
                v3_in = True

            if ((x1*x1+y2*y2)**0.5 - R) <= 0.:
                v4_in = True

            if v1_in == True and v2_in == True and v3_in == True and v4_in == True:
                field.cells[i,j] = 1.

            elif v1_in == False and v2_in == False and v3_in == False and v4_in == False:
                field.cells[i,j] = 0.
            else:
                field.cells[i,j] = -3

            if (field.cells[i,j] == -3 and x_c > 0. and y_c > 0.):
                print("Cutting cell: (%d, %d)"%(i,j))

                area, xp1, xp2, yp1, yp2 = 0, 0, 0, 0, 0

                if y1 <= R:
                    xp1 = (R*R-y1*y1)**0.5
                else:
                    xp1 = 100.

                if y2 <= R:
                    xp2 = (R*R-y2*y2)**0.5
                else:
                    xp2 = 100.

                if x2 <= R:
                    yp2 = (R*R-x2*x2)**0.5
                else:
                    yp2 = 100.

                xp1_in, xp2_in, yp1_in, yp2_in = False, False, False, False

                if x1 < xp1 and x2 > xp1:
                    xp1_in = True

                if x1 < xp2 and x2 > xp2:
                    xp2_in = True

                if y1 < yp2 and y2 > yp2:
                    yp2_in = True

                if yp1_in == True and yp2_in == False and xp1_in == True and xp2_in == False:
                    integral = ((R*R)/2.)*((xp1/R)*(1.-(xp1/R)**2.)**0.5+np.arcsin(xp1/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (xp1-x1)*y1


                if yp1_in == True and yp2_in == True and xp1_in == False and xp2_in == False:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (x2-x1)*y1


                if yp1_in == False and yp2_in == False and xp1_in == True and xp2_in == True:
                    integral = ((R*R)/2.)*((xp1/R)*(1.-(xp1/R)**2.)**0.5+np.arcsin(xp1/R)) - ((R*R)/2.)*((xp2/R)*(1.-(xp2/R)**2.)**0.5+np.arcsin(xp2/R))
                    area = integral + (xp2-x1)*y2 - (xp1-x1)*y1

                if yp1_in == False and yp2_in == True and xp1_in == False and xp2_in == True:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((xp2/R)*(1.-(xp2/R)**2.)**0.5+np.arcsin(xp2/R))
                    area = integral + (xp2-x1)*y2 - (x2-x1)*y1

                if xp2 == x1 and yp2_in == True:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (x2-x1)*y1

                if xp2 == x1 and xp1_in == True:
                    integral = ((R*R)/2.)*((xp1/R)*(1.-(xp1/R)**2.)**0.5+np.arcsin(xp1/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (xp1-x1)*y1

                if xp2 == x1 and yp2 == y1:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (x2-x1)*y1

                if yp2 == y1 and yp1_in == True:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((x1/R)*(1.-(x1/R)**2.)**0.5+np.arcsin(x1/R))
                    area = integral - (x2-x1)*y1

                if yp2 == y1 and xp2_in == True:
                    integral = ((R*R)/2.)*((x2/R)*(1.-(x2/R)**2.)**0.5+np.arcsin(x2/R)) - ((R*R)/2.)*((xp2/R)*(1.-(xp2/R)**2.)**0.5+np.arcsin(xp2/R))
                    area = integral + (xp2-x1)*y2 -(x2-x1)*y1

                if x1 == xp1:
                    area = 0

                if y2 == yp2:
                    area = deltaY*deltaY

                field.cells[i,j] = area/(deltaY*deltaY)
                markerAlpha[i,j] = 1.

    # Multiply by 4
    for i in range(mesh.nx):
        for j in range(mesh.ny):

            x_c = mesh.Cx[i]
            y_c = mesh.Cy[j]

            if x_c > 0. and x_c > 0. and ((x_c*x_c + y_c*y_c)**0.5 - R*1.1) < 0:

                deltaY = mesh.step
                x_c = mesh.Cx[i]
                y_c = mesh.Cy[j]

                x1 = x_c-deltaY/2.
                x2 = x_c+deltaY/2.
                y1 = y_c-deltaY/2.
                y2 = y_c+deltaY/2.

                for ii in range(mesh.nx):
                    for jj in range(mesh.ny):

                        if (field.cells[ii,jj] == -3):
                            if (y1 < np.abs(mesh.Cy[jj])) and (np.abs(mesh.Cy[jj]) < y2) and (x1 < np.abs(mesh.Cx[ii])) and (np.abs(mesh.Cx[ii]) < x2):

                               field.cells[ii,jj] = field.cells[i,j]

