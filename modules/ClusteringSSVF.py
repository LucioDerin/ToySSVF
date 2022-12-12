import numpy as np
import math as m
from scipy.optimize import minimize
from sklearn.cluster import KMeans as KM

class ClusteringSSVF:
    '''
    Class that implements the vertex finding and vertex fitting of the Single Secondary Vertex Finding
    algorithm based on Unsupervised Learning CLustering algorithm (KMeans).

    Members:
    @self.k: number of clusters to look for;
    @self.max_iter: max iter for KMeans;
    @self.chi2Threshold: chi2 threshold value after which vertex fitting loop stops;

    Public Methods:
    @vertexFinder(jet): returns a list of tracks' belonging cluster index;
    @vertexFitter(couples): return np.array of shape (3,), the fitted SV;

    Private Methods:
    @__chi2(vertex,tracks): chi2 to be minimized during vertex fitting step;
    @__evaluateChi2(vfit,tracks): returns the list of the track's chi2 wrt the current fitted vertex;
    '''
    # Constructor
    def __init__(self, k = 3, max_iter = 300, chi2Threshold = 0.1) -> None:
        self.k = k
        self.max_iter = max_iter
        self.chi2Threshold = chi2Threshold

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
        Implements the tracks clustering step.
        Parameters:
        @jet: instance of Jet, the event on which apply SSVF;
        Returns:
        @clusters: np.array of shape (nTracks,), the belonging cluster of each track;
        '''
        jetForClustering = []
        for t in jet.tracksPV:
            track = [np.abs(t.IP),t.versor[2]]
            jetForClustering.append(track)
        for t in jet.tracksSV:
            track = [np.abs(t.IP),t.versor[2]]
            jetForClustering.append(track)
        for t in jet.tracksPileup:
            track = [np.abs(t.IP),t.versor[2]]
            jetForClustering.append(track)
        
        jetForClustering = np.asarray(jetForClustering)
        km = KM(self.k, max_iter=self.max_iter)
        km.fit(jetForClustering)
        clusters = km.labels_
        return clusters


    def vertexFitter(self,jet,clusters):
        '''
        Implements the vertex fitting step.
        Parameters:
        @jet: the jet on which the clustering has been executed;
        @clusters: tracks' belonging cluster indices;
        Returns:
        @SVfit: np.array of shape (3,), the fitted secondary vertex;
        '''

        allTracks = []
        # Flatten the track list
        for t in jet.tracksPV:
            allTracks.append(t)        
        for t in jet.tracksSV:
            allTracks.append(t)        
        for t in jet.tracksPileup:
            allTracks.append(t)

        # create clusters' track lists
        cluster0Tracks = []
        cluster1Tracks = []
        cluster2Tracks = []

        for i,t in enumerate(allTracks):
            if clusters[i] == 0:
                cluster0Tracks.append(t)
            elif clusters[i] == 1:
                cluster1Tracks.append(t)
            else:
                cluster2Tracks.append(t)  

        # First guess on the SV: flight length on the z components and non zero x and y
        SVfit0 = [0.01,0.01,10.,]
        # Fitting a vertex on each cluster
        SVfitC0 = minimize(self.__chi2,x0=SVfit0,args=cluster0Tracks).x
        SVfitC1 = minimize(self.__chi2,x0=SVfit0,args=cluster1Tracks).x
        SVfitC2 = minimize(self.__chi2,x0=SVfit0,args=cluster2Tracks).x
        clusterVertices = np.array([SVfitC0,SVfitC1,SVfitC2])

        # Evaluate the per-track per-cluster chi2 of the current fitted vertices
        chi2sC0 = self.__evaluateChi2(SVfitC0,cluster0Tracks)
        chi2sC1 = self.__evaluateChi2(SVfitC1,cluster1Tracks)
        chi2sC2 = self.__evaluateChi2(SVfitC2,cluster2Tracks)
        # Evaluate the total chi2 per-cluster
        chi2sC0 = np.sum(chi2sC0)
        chi2sC1 = np.sum(chi2sC1)
        chi2sC2 = np.sum(chi2sC2)
        chi2Clusters = np.array([chi2sC0,chi2sC1,chi2sC2])
        # If the chi2 of a cluster is below the threshold, return the respective vertex
        if np.any(chi2Clusters <= self.chi2Threshold):
            return np.array(clusterVertices[np.argmin(chi2Clusters)])
        # Else, vertex fitting failed
        else:
            return None