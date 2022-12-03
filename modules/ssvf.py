import numpy as np
import math as m
from scipy.optimize import minimize

class SSVF:
    '''
    Class that implements the vertex finding and vertex fitting of the Single Secondary Vertex Finding
    algorithm.

    Members:
    @self.dThreshold: distance threshold under which tracks are coupled in vertex finding step;
    @self.chi2Threshold: chi2 threshold value after which vertex fitting loop stops;

    Public Methods:
    @vertexFinder(jet): returns list of list of Track instances, the tracks
        coupled during the vertex finding step;
    @vertexFitter(couples): return np.array of shape (3,), the fitted SV, and a list of 
        Track instances, the tracks that have been fitted to the SV.

    Private Methods:
    @__distances(t,tracks): returns the array of distances between the track t and all of the tracks in the list;
    @__chi2(vertex,tracks): chi2 to be minimized during vertex fitting step;
    @__evaluateChi2(vfit,tracks): returns the list of the track's chi2 wrt the current fitted vertex;
    '''

    # Constructor
    def __init__(self, dThreshold = 6*0.001, chi2Threshold = 0.1) -> None:
        self.dThreshold = dThreshold
        self.chi2Threshold = chi2Threshold

    def __distances(self,t,tracks):
        '''
        Returns the array of distances between the track t and all of the tracks in the list.
        Parameters:
        @t: current track;
        @tracks: list of tracks from which evaluate the distances;
        Returns:
        @distances: np.array of the distances of track t from all of the tracks in the list;
        '''
        distances = []
        for ts in tracks:
            distances.append(t.trackDistance(ts))     
        return np.asarray(distances)

    def __chi2(self,vertex,tracks):
        '''
        Evaluates the chi2 of the current vertex.
        Parameters:
        @vertex: np.array of shape (3,), current fitted secondary vertex;
        @tracks: current list of tracks supposedly belonging to the secondary vertex;
        Returns:
        @chi2: chi2 of the fit, defined as the sum of the distances of the tracks from the vertex;
        '''
        vertex = np.asarray(vertex)
        chi2 = 0
        for t in tracks:
            chi2 += t.pointDistance(vertex)
        return chi2

    def __evaluateChi2(self,vfit,tracks):
        '''
        Evaluates the tracks' chi2 from the current vertex.
        Parameters:
        @vfit: np.array of shape (3,), current fitted secondary vertex;
        @tracks: current list of tracks supposedly belonging to the secondary vertex;
        Returns:
        @chi2: np.array, the chi2 of the tracks from the fitted vertex;
        '''
        chi2 = []
        for t in tracks:
            chi2.append(t.pointDistance(vfit))
        return np.asarray(chi2)

    def vertexFinder(self,jet):
        '''
        Implements the vertex finding step.
        Parameters:
        @jet: instance of Jet, the event on which apply SSVF;
        Returns:
        @couples: list of lists of Track, where the entry i is the i-th couple of tracks found
            by the algorithm;
        '''
        couples = []
        # Joining all the jets' track
        tracks = jet.tracksPV + jet.tracksSV + jet.tracksPileup
        # Track indices
        indices = [i for i in range(len(tracks))]
        # Coupling loop
        while len(tracks)>1:
            # Pop the current track that will be coupled
            i,t = indices.pop(0),tracks.pop(0)
            # Evaluate the distances wrt all of the other tracks
            distances = self.__distances(t,tracks)
            # Find the closest track
            jt = np.argmin(distances)
            # If the track are close enough, couple them
            if distances[jt]<self.dThreshold:
                j = indices[jt]
                tj = tracks.pop(jt)
                couples.append([t,tj])

        return couples


    def vertexFitter(self,couples):
        '''
        Implements the vertex fitting step.
        Parameters:
        @couples: couples of tracks belonging to the same vertex, 
            the result of the vertex finding step;
        Returns:
        @SVfit: np.array of shape (3,), the fitted secondary vertex;
        @selectedTracks: list of Track instances, the tracks supposed to belong to the fitted SV;
        '''

        selectedTracks = []
        # Flatten the coupled track list
        for c in couples:
            for t in c:
                selectedTracks.append(t)
        chi2 = 10000.
        # First guess on the SV: flight length on the z components and non zero x and y
        SVfit0 = [0.01,0.01,10.,]
        # Fitting loop
        while True:
            # Minimize the chi2 with the current tracks
            SVfit = minimize(self.__chi2,x0=SVfit0,args=selectedTracks).x
            # Evaluate the per-track chi2 of the current fitted vertex
            chi2s = self.__evaluateChi2(SVfit,selectedTracks)
            # Evaluate the total chi2
            currentChi2 = np.sum(chi2s)
            # If the chi2 is below the threshold, stop
            if currentChi2< self.chi2Threshold:
                break
            # Else, reject the worst track
            i = np.argmax(chi2s)
            selectedTracks.pop(i)
            # If no more couples are there, stop and return invalid results
            if len(selectedTracks)<2:
                return None,None
        return SVfit, selectedTracks