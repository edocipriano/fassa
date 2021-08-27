import numpy as np

class index_conversion:

    def __init__(self, mesh):
        self.mesh = mesh
        self.nx   = mesh.nx
        self.ny   = mesh.ny

    def update(self, i, j):

        self.l = (i-1)*self.ny + j

        self.P = {"grid"    : np.index_exp[i,j],
                  "storage" : np.index_exp[self.l]}

        self.N = {"grid"    : np.index_exp[i,j+1],
                  "storage" : np.index_exp[self.l+1]}

        self.S = {"grid"    : np.index_exp[i,j-1],
                  "storage" : np.index_exp[self.l-1]}

        self.E = {"grid"    : np.index_exp[i+1,j],
                  "storage" : np.index_exp[self.l+self.ny]}

        self.W = {"grid"    : np.index_exp[i-1,j],
                  "storage" : np.index_exp[self.l-self.ny]}


