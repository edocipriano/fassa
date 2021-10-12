import fassa
import scipy.sparse.linalg as linalg
import time

meshDict = {"Lx" : 1,
            "Ly" : 1,
            "nx" : 10,
            "ny" : 10}

mesh = fassa.hexCMesh(meshDict)

phiDict = {"N" : ("fixedValue", 1.),
           "S" : ("fixedValue", 0.),
           "E" : ("fixedValue", 0.),
           "W" : ("fixedValue", 0.)}

phi = fassa.Field("phi", phiDict, mesh, 0.)
phi.updateBCs()

start = time.time()

A, b = fassa.solve.implicit.laplacian(1., phi)
x = linalg.spsolve(A,b)

end = time.time()
print("Elapsed cpu time: ", end-start)

fassa.utils.write.writeVtk(phi)

