from modules.JetGenerator import JetGenerator
from modules.SSVF import SSVF
from matplotlib import pyplot as plt
import numpy as np

if __name__=="__main__":
    
    generator = JetGenerator(100)
    jets = generator.generate()

    ssvf = SSVF()
    distances = []

    for i,jet in enumerate(jets):
        if i%10==0:
            print("Working... ", i/len(jets), "...")
        couples = ssvf.vertexFinder(jet)
        fittedVertex,_ = ssvf.vertexFitter(couples)
        if np.any(fittedVertex)==None:
            continue
        distances.append(np.linalg.norm(fittedVertex-jet.SV))
    distances = np.array(distances)
    fig,(a0,a1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4, 1]},figsize=(7,6))
    titleStr = "SSVF test on 100 jets" 
    plt.suptitle(titleStr)
    counts, bins, _ = a0.hist(distances)
    a0.set_xlabel("$||SV_{real} - SV_{fit}||$ [mm]")
    a0.set_ylabel("Counts")
    a0.semilogy()
    xlimHist = a0.get_xlim()

    x = np.array(bins[:-1]) + ((bins[1]-bins[0])/2.)
    ey = np.sqrt(counts)/np.sum(counts)
    a1.errorbar(x,counts/np.sum(counts),yerr=ey,fmt="o")
    a1.grid()
    a1.set_xlim(xlimHist)
    a1.set_ylabel("Counts/total")
    a1.semilogy()
    plt.tight_layout()
    plt.show()