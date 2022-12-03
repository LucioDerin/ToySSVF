### Toy reproduction of SSVF algorithm

Repository content:
- `.modules/Track`: Class to represent a tracks in the jet;
- `.modules/Jet`: Class to represent a jet;
- `.modules/JetGenerator`: Class that implements a simple jets Montecarlo generator;
- `.modules/ssvf`: Class that implements the SSVF algorithm;
- `./test.py`: test of the modules;
- `./notes/notes.md`: notes on the implementation of SSVF and of the tracks linear algebra;
#### Result of the SSVF algorithm
<center>
<img src="./screenshots/jetEventCone2.jpg" alt="drawing" width="700"/>
</center>

<center>
<img src="./screenshots/event2.jpg" alt="drawing" width="700"/>
</center>

Interesting event, with a track from the PV that passes near the SV and is consequently  misclassified as a SV track:
<center>
<img src="./screenshots/noise1.png" alt="drawing" width="700"/>
</center>
