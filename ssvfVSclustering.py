import time as t
import numpy as np
from matplotlib import pyplot as plt
import math as m

from modules.JetGenerator import JetGenerator
from modules.SSVF import SSVF
from modules.ClusteringSSVF import ClusteringSSVF

if __name__ == "__main__":
    
    nTracks = np.linspace(2,10,6)

    ssvf = SSVF()
    clustering = ClusteringSSVF()

    distancesSSVF = []
    eDistancesSSVF = []
    timeSSVF = []
    distancesClustering = []
    eDistancesClustering = []
    timeClustering = []
    
    for n in nTracks:
        n = m.floor(n)
        print("Working.. n =",n,"...")
        generator = JetGenerator(10,nTracksPV=n,nTracksSV=n+1,nTracksPileup=n)
        jets = generator.generate()

        # SSVF
        distances = []
        start = t.time()
        for jet in jets:
            couples = ssvf.vertexFinder(jet)
            end = t.time()
            fittedVertex,_ = ssvf.vertexFitter(couples)
            if np.any(fittedVertex)==None:
                continue
            distances.append(np.linalg.norm(fittedVertex-jet.SV))
        distances = np.array(distances)
        distancesSSVF.append(np.mean(distances))
        eDistancesSSVF.append(np.std(distances)/np.sqrt(len(distances)))
        timeSSVF.append(end-start)

        # Clustering
        distances = []
        start = t.time()
        for jet in jets:
            clusters = clustering.vertexFinder(jet)
            end = t.time()
            fittedVertex = clustering.vertexFitter(jet,clusters)
            if np.any(fittedVertex)==None:
                continue
            distances.append(np.linalg.norm(fittedVertex-jet.SV))
        distances = np.array(distances)
        distancesClustering.append(np.mean(distances))
        eDistancesClustering.append(np.std(distances)/np.sqrt(len(distances)))
        timeClustering.append(end-start)

    # plot
    nTracks = 3*nTracks
    fig,(a0,a1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4, 1]},figsize=(7,6))
    a0.plot(nTracks,timeSSVF,color='orange',label="Time SSVF")
    a0.plot(nTracks,timeClustering,color='blue',label="Time C-SSVF")
    a0.set_ylabel("Execution Time [s]")
    a0.legend(loc='best')
    a01 = a0.twinx()
    a01.errorbar(nTracks,distancesSSVF,yerr=eDistancesSSVF,color='orange',linestyle='dashed',label="Error SSVF")
    a01.errorbar(nTracks,distancesClustering,yerr=eDistancesClustering,color='blue',linestyle='dashed',label='Error C-SSVF')
    a01.set_ylabel("$\sqrt{MSE}$ [mm]")
    xlim0 = a0.get_xlim()
    a01.legend(loc='best')

    y = np.array(timeSSVF)/np.array(timeClustering)
    a1.scatter(nTracks,y)
    a1.grid()
    a1.set_ylabel("$t_{SSVF}/t_{C-SSVF}$")
    a1.set_xlabel("$N_{tracks}}$")
    plt.suptitle("SSVF vs C-SSVF")
    plt.tight_layout()
    plt.show()