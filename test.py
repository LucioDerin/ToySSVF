from matplotlib import pyplot as plt
from modules.JetGenerator import JetGenerator
from modules.SSVF import SSVF

# Small test
if __name__ == "__main__":
    seed = input("Seed:")
    if seed == "":
        seed = 1
    else:
        seed = int(seed)
    generator = JetGenerator(seed=seed)
    jets = generator.generate()

    for jet in jets:
        jet.print()

    vf = SSVF()
    couples = vf.vertexFinder(jets[0])
    vertex, selectedTracks = vf.vertexFitter(couples)
    jets[0].draw(couples=couples, fittedSV=vertex,
                 fittedTracks=selectedTracks, drawJetCone=True, plotly=False)
    plt.tight_layout()
    plt.show()
