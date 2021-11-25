class SurfaceTensionForce:
    """
    Base class for surface tension force models
    """

    def __init__(self, mesh):
        self.mesh = mesh

    def setCurvature(self, K):
        self.K = K

    def compute(self):


def main():
    surfaceTensionDict = {"sigma" : 0.0073,
                          "Ktype" : 0}

