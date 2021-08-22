def velocityFromStreamfunction(u,v):
    """
    Compute staggered  velocity fields from the stream function
    """
    mesh = u.mesh

    # u = d(phi)/dy
    for i in range(1, mesh.nx+1):
        for j in range(1, mesh.ny):
             u.cells[i,j] = (phi.nodes[i,j] - phi.nodes[i,j-1])/(mesh.y[j] - mesh.y[j-1])

    # v = -d(phi)/dx
    for i in range(1, mesh.nx):
        for j in range(1, mesh.ny+1):
            v.cells[i,j] = -(phi.nodes[i,j] - phi.nodes[i-1,j])/(mesh.x[i] - mesh.x[i-1])
