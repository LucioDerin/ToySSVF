import numpy as np
import math as m

class Track:
    
    def __init__(self, origin=[], versor=[]):
        self.origin = np.asarray(origin)
        self.versor = np.asarray(versor)
    
    def evaluate(self,t):
        P = self.origin + self.versor*t
        return P

    def trackDistance(self,track):
        alpha = np.arccos(np.clip(np.dot(self.versor, track.versor), -1.0, 1.0))
        DP0 = track.origin - self.origin
        tstar = -(1/m.sin(alpha)**2)*(m.cos(alpha)*(np.dot(DP0,track.versor)-np.dot(DP0,track.versor)))
        sstar = m.cos(alpha)*tstar - np.dot(DP0,track.versor)
        distance = np.linalg.norm((self.evaluate(tstar)-track.evaluate(sstar)))
        return distance
    
    def pointDistance(self,P):
        d = -np.sum(self.versor*P)
        sP0 = np.sum(self.origin*self.versor)
        tstar = -(sP0+d)
        Pprime = self.evaluate(tstar)
        return np.linalg.norm((P-Pprime))

    def print(self):
        print("- Track origin's coordinates: [", self.origin[0],",",self.origin[1],",",self.origin[2],"] mm;")
        print("- Track versor's components: [", self.versor[0],",",self.versor[1],",",self.versor[2],"];")
        print(" ")