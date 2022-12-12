import numpy as np
import math as m

class Track:
    '''
    Class that implements the representation of a track as a straight line with an origin and a versor.
    Public Members:
    @self.origin: np.array of shape (3,), (x,y,z) coordinates of the track's origin;
    @self.versor: np.array of shape (3,), (x,y,z) components of the track's versor
        with current choice of reference system they are (cos(phi)sin(theta),sin(phi)sin(theta),cos(theta))

    Public Methods:
    @evaluate(t): returns np.array of shape (3,), evaluates the parametric representation of the line 
        where t is the value of the line's parameter; 
    @trackDistance(track): returns double, the minimum distance of the current track 
        from the track given in input;
    @pointDistance(P): returns double, the minimum distance of the track from the point P;
    @print(): prints the track's origin and versor;
    '''

    # Constructor
    def __init__(self, origin=[], versor=[]):
        '''
        Constructor of the class.
        Parameters:
        @origin: origin of the track, np.array or list of shape (3,);
        @versor: versor of the track, np.array or list of shape (3,);
        '''
        self.origin = np.asarray(origin)
        self.versor = np.asarray(versor)
        self.IP = self.pointDistance([0,0,0])
    
    def evaluate(self,t):
        '''
        Evaluates the parametric representation of the track at the parameter value t.
        Parameters:
        @t: value of the parameter of the parametric representation of the line;
        Returns:
        @P: np.array of shape (3,), the track's point at parameter = t;
        '''
        P = self.origin + self.versor*t
        return P

    def trackDistance(self,track):
        '''
        Evaluates the minimum distance between the current track and the input one.
        Parameters:
        @track: instance of Track, track from which evaluate the distance;
        Returns:
        @distance: double, the minimum distance between the tracks;
        '''
        # Angle between tracks' versors 
        alpha = np.arccos(np.clip(np.dot(self.versor, track.versor), -1.0, 1.0))
        # Handy variable (defined in the notes)
        DP0 = track.origin - self.origin
        # Values of the parameters in the tracks' parametric representation 
        # at the closest points between the two tracks
        tstar = -(1/m.sin(alpha)**2)*(m.cos(alpha)*(np.dot(DP0,track.versor)-np.dot(DP0,track.versor)))
        sstar = m.cos(alpha)*tstar - np.dot(DP0,track.versor)
        # Distance between the closest points of the two tracks
        distance = np.linalg.norm((self.evaluate(tstar)-track.evaluate(sstar)))
        return distance
    
    def pointDistance(self,P):
        '''
        Evaluates the minimum distance between the current track and the input point.
        Parameters:
        @P: np.array or list of shape(3,), point from which evaluate the distance;
        Returns:
        @distance: double, the minimum distance between the track and the point;
        '''
        P = np.asarray(P)
        # Constant of the plane containing P and orthogonal to the track
        d = -np.sum(self.versor*P)
        # Handy variable (defined in the notes)
        sP0 = np.sum(self.origin*self.versor)
        # Value of the parameter in the track's parametric representation
        # at the point of minimum distance from P
        tstar = -(sP0+d)
        # Track's closest point to P
        Pprime = self.evaluate(tstar)
        # return distance between P and Pprime
        return np.linalg.norm((P-Pprime))

    def print(self):
        '''
        Prints the track's origin and versor.
        '''
        print("- Track origin's coordinates: [", self.origin[0],",",self.origin[1],",",self.origin[2],"] mm;")
        print("- Track versor's components: [", self.versor[0],",",self.versor[1],",",self.versor[2],"];")
        print(" ")