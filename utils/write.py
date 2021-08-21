import numpy as np

def writeVtkMesh(mesh):

    x_dim = len(mesh.x)
    y_dim = len(mesh.y)
    z_dim = 1

    nNodes = x_dim*y_dim*z_dim

    f = open("Mesh.vtk", "w")
    f.write("# vtk DataFile Version 3.0\n")
    f.write("vtk output\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_GRID\n")
    f.write("DIMENSIONS "+str(x_dim)+" "+str(y_dim)+" "+str(1)+"\n")
    f.write("POINTS "+str(nNodes)+" float\n")

    for j in range(z_dim):
        for j in range(y_dim):
            for i in range(x_dim):
                f.write(str(mesh.x[i])+" "+str(mesh.y[j])+" "+str(0.0)+"\n")

    f.close()


def writeVtk(path, mesh, *fields):

    x_dim = len(mesh.x)
    y_dim = len(mesh.y)
    z_dim = 1

    nNodes = x_dim*y_dim*z_dim
    nCells = mesh.nx*mesh.ny

    f = open(path, "w")
    f.write("# vtk DataFile Version 3.0\n")
    f.write("vtk output\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_GRID\n")
    f.write("DIMENSIONS "+str(x_dim)+" "+str(y_dim)+" "+str(1)+"\n")
    f.write("POINTS "+str(nNodes)+" float\n")

    for j in range(z_dim):
        for j in range(y_dim):
            for i in range(x_dim):
                f.write(str(mesh.x[i])+" "+str(mesh.y[j])+" "+str(0.0)+"\n")

    # For each field write cell and point data
    for field in fields:
        f.write("\n")
        f.write("CELL_DATA "+str(nCells)+"\n")
        f.write("SCALARS "+field.name+" float 1\n")
        f.write("LOOKUP_TABLE "+field.name+"_cell_table\n")
        for j in range(1, mesh.ny+1):
            for i in range(1, mesh.nx+1):
                f.write(str(field.cells[i,j])+" ")
        f.write("\n")

        f.write("POINT_DATA "+str(nNodes)+"\n")
        f.write("SCALARS "+field.name+" float 1\n")
        f.write("LOOKUP_TABLE "+field.name+"_point_table\n")
        field.nodeInterpolation()
        for j in range(y_dim):
            for i in range(x_dim):
                f.write(str(field.nodes[i,j])+" ")

    f.close()

