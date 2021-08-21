def fieldAverage(alpha1, field, value1, value2):
    mesh = alpha1.mesh

    for i in range(mesh.nx+2):
        for j in range(mesh.ny+2):
            field.cells[i,j] = value1*alpha1.cells[i,j] + value2*(1.-alpha1.cells[i,j])


