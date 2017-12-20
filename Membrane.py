from Lipid import Lipid
import numpy as np

# Membrane Class
class Membrane(Lipid):
    def __init__(self, system, numLipids):
        self.system     = system
        self.numLipids  = numLipids
        self.size       = Lipid.initRestLength + 0.1
        self.lipid      = []                            # Contains all the lipid objects in the membrane

    def setOrientation(self, orient='monolayer', lipid1 = None, lipid2 = None):

        if orient == 'monolayer':
            # Lipids in x direction
            numX = int(np.round(np.sqrt(self.numLipids)))
            # Lipids in y direction
            numY = int(np.round(self.numLipids / numX))

            startX      = self.system.box_l[0] / 2 - numX * self.size / 2
            startY      = self.system.box_l[0] / 2 + numY * self.size / 2
            startZ      = self.system.box_l[0] / 2
            position    = np.array([startX,	startY,	startZ])

            movRight    = np.array([self.size, 0, 0])
            movDown     = np.array([0, -self.size, 0])

            for i in range(numX):
                for j in range(numY):
                    l = Lipid(self.system, position)
                    l.add()
                    self.lipid.append(l)
                    position += movRight
                position[0] = startX
                position += movDown
            print("\nCreated a monolayer of {} lipids".format(self.numLipids))

        elif orient == 'bilayer':
            # Lipids in x direction
            numX = int(np.round(np.sqrt(self.numLipids / 2)))
            # Lipids in y direction
            numY = int(np.round(self.numLipids / (2 * numX)))

            startX  = self.system.box_l[0] / 2 - numX * self.size / 2
            startY  = self.system.box_l[0] / 2 + numY * self.size / 2
            startZ  = self.system.box_l[0] / 2
            position= np.array([startX,	startY,	startZ])

            movRight= np.array([self.size, 0, 0])
            movDown = np.array([0, -self.size, 0])
            layer2  = np.array([0, 0, Lipid.initRestLength * 3.0])

            for i in range(numX):
                for j in range(numY):
                    lUp = Lipid(self.system, position)
                    lDown = Lipid(self.system, position - layer2, theta=np.pi)
                    lUp.add()
                    lDown.add()
                    self.lipid.append(lUp)
                    self.lipid.append(lDown)
                    position += movRight
                position[0] = startX
                position += movDown
            print("\nCreated a bilayer of {} lipids".format(self.numLipids))

        elif orient == 'mixedbilayer':

            if lipid2 == None:
                raise Exception("Input a second lipid instance for mixedbilayer!")

            # Lipids in x direction
            numX = int(np.round(np.sqrt(self.numLipids / 2)))
            # Lipids in y direction
            numY = int(np.round(self.numLipids / (2 * numX)))

            startX  = self.system.box_l[0] / 2 - numX * self.size / 2
            startY  = self.system.box_l[0] / 2 + numY * self.size / 2
            startZ  = self.system.box_l[0] / 2
            position= np.array([startX, startY, startZ])

            movRight= np.array([self.size, 0, 0])
            movDown = np.array([0, -self.size, 0])
            layer2  = np.array([0, 0, Lipid.initRestLength * 3.0])

            for i in range(numX):
                for j in range(numY):
                    lUp = Lipid(self.system, position, epsilon=lipid1.epsilon, sigma = lipid1.sigma ,lipidType = lipid1.type)
                    lDown = Lipid(self.system, position - layer2, theta=np.pi, epsilon=lipid2.epsilon, sigma = lipid2.sigma ,lipidType = lipid2.type)
                    self.lipid.append(lUp)
                    self.lipid.append(lDown)
                    lUp.add()
                    lDown.add()
                    position += movRight
                position[0] = startX
                position += movDown
            print("\nCreated a mixedbilayer of {} lipids".format(self.numLipids))

        elif orient == 'random':
            for i in range(self.numLipids):
                l = Lipid(self.system, midPos=np.random.random(3) * self.system.box_l)
                l.add()
                self.lipid.append(l)
            print("\nPlaced {} lipids randomly".format(self.numLipids))

        else:
            print("\nOrientation not valid!\n")