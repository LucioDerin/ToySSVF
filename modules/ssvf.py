import numpy as np
import math as m
from scipy.optimize import minimize

class ssvf:

    chi2Treshold = 0.1

    def __init__(self, dTreshold = 0.001) -> None:
        self.dTreshold = dTreshold

    def __distances(self,t,tracks):
        distances = []
        for ts in tracks:
            distances.append(t.trackDistance(ts))     
        return np.asarray(distances)

    def __chi2(self,v,tracks):
        v = np.asarray(v)
        chi2 = 0
        for t in tracks:
            chi2 += t.pointDistance(v)
        return chi2

    def __evaluateChi2(self,vfit,tracks):
        chi2 = []
        for t in tracks:
            chi2.append(t.pointDistance(vfit))
        return np.asarray(chi2)

    def vertexFinder(self,jet):
        couples = []
        tracks = jet.tracksPV + jet.tracksSV
        indices = [i for i in range(len(tracks))]

        while len(tracks)>1:
            i,t = indices.pop(0),tracks.pop(0)
            distances = self.__distances(t,tracks)
            jt = np.argmin(distances)
            j = indices[jt]
            tj = tracks.pop(jt)
            if distances[jt]<self.dTreshold:
                couples.append([t,tj])

        return couples


    def vertexFitter(self,couples):
        selectedTracks = []
        for c in couples:
            for t in c:
                selectedTracks.append(t)
        chi2 = 10000.
        SVfit0 = [0.01,0.01,10.,]
        while True:
            SVfit = minimize(self.__chi2,x0=SVfit0,args=selectedTracks).x
            chi2s = self.__evaluateChi2(SVfit,selectedTracks)
            currentChi2 = np.sum(chi2s)
            if currentChi2< self.chi2Treshold:
                break
            i = np.argmax(chi2s)
            selectedTracks.pop(i)
            if len(selectedTracks)<2:
                return -1,[-1]
        return SVfit, selectedTracks

