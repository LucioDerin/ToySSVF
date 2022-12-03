import numpy as np
import math as m
from modules.Track import Track
from modules.Jet import Jet

class JetGenerator:
    '''
    Class that generates jet events, with PV, SV and Pileup tracks.
    Public Members:
    @self.nEvents: numbers of events to generate;
    @self.nTracksPV: mean number of tracks in the primary vertex;
    @self.nTracksSV: mean number of tracks in the secondary vertex;
    @self.nTracksPileup: mean number of tracks in pileup events;
    @self.flightLength, self.sigmaFlight: mean and stdev of the B meson flight length [mm];
    @self.sigmaTheta: stdev of the SV angle with the z axis [rad];
    @self.thetaMaxSV: maximum angle between SV and the z axis [rad];
    @self.thetaMaxJet: maximum angle of aperture of the jet [rad];
    @self.sigmaCoord: smearing on the coordinates of the tracks [mm];
    @self.nTracksRange: width of the interval in which nTracks numbers are uniformly sampled;

    Public Methods:
    @generate: returns list of nEvents jet, generates the jets;
    '''

    # Physical constants
    # average B meson flight length [mm]
    flightLength = 10
    # standard deviation of the flight length [mm]
    sigmaFlight = 0.1
    # Standard deviation of the theta angle of the SV's tracks [mm]
    sigmaTheta = 0.1
    # SV maximum theta [rad]
    thetaMaxSV = 0.2

    # jets' cone maximum aperture [rad]
    thetaMaxJet = 0.4

    # smearing on the coordinates [mm]
    sigmaCoord = 1e-3

    # width of the range of the per-vertex random number of tracks
    nTracksRange = 2

    # Constructor
    def __init__(self,nEvents=1, nTracksPV = 4, nTracksSV = 5, nTracksPileup = 3):
        self.nEvents = nEvents
        self.nTracksPV = nTracksPV
        self.nTracksSV = nTracksSV
        self.nTracksPileup = nTracksPileup

    def generate(self):
        '''
        Generates the jets.
        Returns:
        @jets: list of nEvents Jet instances, the generated jets.
        '''
        # Primary vertex is in the origin of the reference frame
        PV = np.zeros(3)
        jets = []
        generatedEvents = 0
        # Generation loop
        while generatedEvents<self.nEvents:
            # Random number of tracks per type of production
            nPV = np.random.randint(self.nTracksPV-self.nTracksRange if self.nTracksPV>self.nTracksRange else 1, self.nTracksPV+self.nTracksRange)
            nSV = np.random.randint(self.nTracksSV-self.nTracksRange if self.nTracksSV>self.nTracksRange else 1, self.nTracksSV+self.nTracksRange)
            nPileup = np.random.randint(self.nTracksPileup-self.nTracksRange if self.nTracksPileup>self.nTracksRange else 1, self.nTracksPileup+self.nTracksRange)

            # Random secondary vertex generation, spherical coordinates
            # Random distance of the SV from the origin, normally distributed
            rSV = self.flightLength  + np.random.normal(0,self.sigmaFlight)

            # SV theta, normally distributed around 0
            thetaSV = np.random.normal(0,self.sigmaTheta)
            if thetaSV > self.thetaMaxSV:
                thetaSV = 0.
            
            # SV phi, uniformly distributed in [0,2pi]
            phiSV = np.random.uniform(0,2*m.pi)

            # Conversion to cartesian coordinates
            xSV = rSV*m.sin(thetaSV)*m.cos(phiSV)
            ySV = rSV*m.sin(thetaSV)*m.sin(phiSV)
            zSV = rSV*m.cos(thetaSV)
            SV = np.asarray([xSV,ySV,zSV])

            # generating PV tracks
            tracksPV = []
            generatedTracks = 0
            while generatedTracks<nPV:
                # generating track's origin starting from the origin of the sdr (the PV)
                # and applying gaussian smearing
                origin = np.zeros(3) + np.random.normal(0,self.sigmaCoord,3)
                # Track's theta, normally distributed
                theta = np.random.normal(0,self.sigmaTheta*2)
                if theta>self.thetaMaxJet:
                    theta = self.thetaMaxJet
                # Track's phi, uniformly distributed in [0,2pi]
                phi = np.random.uniform(0,2*m.pi)
                # Building track's versor
                versor = [m.sin(theta)*m.cos(phi),m.sin(theta)*m.sin(phi),m.cos(theta)]
                # Adding track to the set
                tracksPV.append(Track(origin,versor))
                generatedTracks += 1
            
            # generating SV tracks
            tracksSV = []
            generatedTracks = 0
            while generatedTracks<nSV:
                # generating track's origin starting from the SV
                # and applying gaussian smearing
                origin = SV + np.random.normal(0,self.sigmaCoord,3)
                # Track's theta, normally distributed
                theta = np.random.normal(0,self.sigmaTheta*2)
                if theta>self.thetaMaxJet:
                    theta = self.thetaMaxJet
                # Track's phi, uniformly distributed in [0,2pi]
                phi = np.random.uniform(0,2*m.pi)
                # Building track's versor
                versor = [m.sin(theta)*m.cos(phi),m.sin(theta)*m.sin(phi),m.cos(theta)]
                # Adding track to the set
                tracksSV.append(Track(origin,versor))
                generatedTracks += 1

            # generating Pileup tracks
            tracksPileup = []
            generatedTracks = 0
            while generatedTracks<nPileup:
                # generating pileup tracks
                z = 0.5*self.flightLength
                theta = 0.5*self.thetaMaxJet
                phi = np.random.uniform(0,2*m.pi)
                # making sure pileup does not cross PV
                r = (z+1)/m.cos(theta)
                x = r*m.sin(theta)*m.cos(phi)
                y = r*m.sin(theta)*m.sin(phi)
                origin = [x,y,z]
                versor = [m.sin(theta)*m.cos(phi),m.sin(theta)*m.sin(phi),m.cos(theta)]
                # Adding track to the set
                tracksPileup.append(Track(origin,versor))
                generatedTracks += 1
            # building event and storing it
            jet = Jet(PV,SV,tracksPV,tracksSV,tracksPileup)
            jets.append(jet)
            generatedEvents += 1
        
        return jets