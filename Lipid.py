import numpy as np
from espressomd.interactions import FeneBond, HarmonicBond

# Lipid Class
class Lipid:
    # Takes in the position of the Middle Tail bead
    initRestLength = 0.95

    def __init__(self, system, midPos=[0.,0.,0.], epsilon=1.0, sigma=1.0, theta=0.0, phi=0.0, lipidType = {"Head": 0, "Mid": 1, "Tail": 1}):
        self.initRestLength = Lipid.initRestLength
        self.system     = system
        self.midPos     = midPos
        self.headPos    = np.zeros([3])
        self.tailPos    = np.zeros([3])
        self.epsilon    = epsilon
        self.sigma      = sigma
        self.theta      = theta
        self.phi        = phi
        self.type       = lipidType

        if theta == 0.0 and phi == 0.0:
            self.headPos = self.midPos + \
                self.initRestLength * np.array([0, 0, 1])
            self.tailPos = self.midPos + \
                self.initRestLength * np.array([0, 0, -1])

        else:
            self.setLipidOrientation()

        self.pos = {"Head": self.headPos, "Mid": midPos, "Tail": self.tailPos}
        self.partId = {}

    def add(self):
        self.addBeads()
        self.setupInternalSprings()


    def setLipidOrientation(self):
        orientVector = self.initRestLength * np.array([np.sin(self.theta) * np.cos(self.phi),
                                                       np.sin(self.theta) *
                                                       np.sin(self.phi),
                                                       np.cos(self.theta)])

        self.headPos = self.midPos + orientVector
        self.tailPos = self.midPos - orientVector

    def addBeads(self):
        for beadName in self.type:
            self.system.part.add(pos=self.pos[beadName], type=self.type[beadName])
            self.partId[beadName] = self.system.part.highest_particle_id

    # Sets up FENE Bond and Harmonic Bond
    def setupInternalSprings(self):

        # Set Up FENE
        k_bond = 30 * (self.epsilon / self.sigma**2)
        d_r_max = 1.5 * self.sigma
        fene = FeneBond(k=k_bond, d_r_max=d_r_max)

        self.system.bonded_inter.add(fene)
        self.system.part[self.partId["Head"]].add_bond((fene, self.partId["Mid"]))
        self.system.part[self.partId["Mid"]].add_bond((fene, self.partId["Tail"]))

        # Set up Harmonic Bond
        k_bend = 10 * (self.epsilon / self.sigma**2)
        r_0 = 4 * self.sigma
        harmonicBond = HarmonicBond(k=k_bend, r_0=r_0)

        self.system.bonded_inter.add(harmonicBond)
        self.system.part[self.partId["Head"]].add_bond(
            (harmonicBond, self.partId["Tail"]))



