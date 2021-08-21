from hexCMesh import hexCMesh

class Time:

    def __init__(self, timeDict, mesh):
        self.totTime = timeDict["totTime"]
        self.safety = timeDict["safety"]
        self.writeNow = timeDict["writeNow"]
        self.runTimeDeltaT = timeDict["runTimeDeltaT"]
        self.deltaT = timeDict["deltaT"]
        self.writeSteps = timeDict["writeSteps"]
        self.value = 0
        self.nTimeSteps = 0
        self.stopNext = False

    def toRun(self):
        """
            Return true if the current time is less than
            the final simulation time
        self.nTimeSteps += 1
        if self.value < self.totTime:
            self.update()
            return True
        else:
            return False
        """
        self.nTimeSteps += 1
        if self.writeNow == True and self.stopNext == True:
            return False
        elif self.writeNow == True and self.stopNext == False:
            self.stopNext = True
            return True
        else:
            if self.value < self.totTime:
                self.update()
                return True
            else:
                return False


    def setDeltaT(self, deltaT):
        self.deltaT = self.safety*deltaT

    def update(self):
        self.value += self.deltaT

    def toWrite(self):
        if self.writeNow:
            return True
        else:
            if self.nTimeSteps%self.writeSteps == 0 or self.nTimeSteps == 1:
                return True
            else:
                return False

def main():

    # Create Mesh and set Boundaries
    meshDict = {"Lx" : 1,
                "Ly" : 1,
                "nx" : 100,
                "ny" : 100}

    mesh = hexCMesh(meshDict)
    mesh.setBasicBoundaries()

    # Create time
    timeDict = {"totTime"       : 20,
                "writeNow"      : False,
                "runTimeDeltaT" : False,
                "safety"        : 1,
                "deltaT"        : 1.e-1,
                "writeSteps"    : 10}

    time = Time(timeDict, mesh)

    while time.toRun():
        print("Steps = %d Time = %f"%(time.nTimeSteps, time.value))
        # Resolve...

        if time.toWrite():
            #utils.write.writeVtk(fields)
            print("    Time to write!")


if __name__=="__main__":
    main()
