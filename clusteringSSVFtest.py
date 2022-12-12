from modules.JetGenerator import JetGenerator
from modules.ClusteringSSVF import ClusteringSSVF
from matplotlib import pyplot as plt
import numpy as np

if __name__=="__main__":
    
    generator = JetGenerator(100)
    jets = generator.generate()

    ssvf = ClusteringSSVF(k=3,max_iter = 10)
    distances = []

    for i,jet in enumerate(jets):
        if i%10==0:
            print("Working... ", i/len(jets), "...")
        clusters = ssvf.vertexFinder(jet)
        fittedVertex = ssvf.vertexFitter(jet,clusters)
        if np.any(fittedVertex)==None:
            continue
        distances.append(np.linalg.norm(fittedVertex-jet.SV))
    distances = np.array(distances)
    fig,(a0,a1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4, 1]},figsize=(7,6))
    titleStr = "C-SSVF test on 100 jets" 
    plt.suptitle(titleStr)
    counts, bins, _ = a0.hist(distances)
    a0.set_xlabel("$||SV_{real} - SV_{fit}||$ [mm]")
    a0.set_ylabel("Counts")
    xlimHist = a0.get_xlim()

    x = np.array(bins[:-1]) + ((bins[1]-bins[0])/2.)
    ey = np.sqrt(x)/np.sum(counts)
    a1.errorbar(x,counts/np.sum(counts),yerr=ey,fmt="o")
    a1.grid()
    a1.set_xlim(xlimHist)
    a1.set_ylabel("Counts/total")
    a1.semilogy()
    plt.tight_layout()
    plt.show()