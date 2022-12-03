from matplotlib import pyplot as plt
from modules.JetGenerator import JetGenerator
from modules.ssvf import ssvf

# Small test
if __name__ == "__main__":
    
    generator = JetGenerator(1)
    jets = generator.generate()

    for jet in jets:
        jet.print()
    

    vf = ssvf()
    couples = vf.vertexFinder(jets[0])
    vertex, selectedTracks = vf.vertexFitter(couples)
    jets[0].draw(couples=couples,fittedSV = vertex,fittedTracks = selectedTracks,drawJetCone=True)
    plt.tight_layout()
    plt.show()