from fassa.hexCMesh import hexCMesh

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
        self.lastStep = False

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
            if self.value <= self.totTime:
                if self.value + self.deltaT >= self.totTime:
                  self.lastStep = True
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
            if self.nTimeSteps%self.writeSteps == 0 or self.nTimeSteps == 1 or self.lastStep == True:
                print ("Writing results")
                return True
            else:
                return False


    def print(self):
        print("Step: %d - Time: %f"%(self.nTimeSteps, self.value))

