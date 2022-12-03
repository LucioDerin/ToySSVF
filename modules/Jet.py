import numpy as np
import math as m
from matplotlib import pyplot as plt
from modules.Track import Track

class Jet:
    # jets' cone maximum aperture [rad]
    thetaMaxJet = 0.4
    colors = ["blue","orange","green","red","purple","brown","pink","gray","olive","cyan"]

    def __init__(self, PV = [], SV = [], tracksPV = [], tracksSV = [], tracksPileup = []):
        self.PV = np.asarray(PV)
        self.SV = np.asarray(SV)
        self.tracksPV = tracksPV
        self.tracksSV = tracksSV
        self.tracksPileup = tracksPileup

    def print(self):
        print("Jet details:")
        print("- Primary vertex coordinates: [", self.PV[0], ",", self.PV[1],",", self.PV[2],"] mm")
        print("- Secondary vertex coordinates: [", self.SV[0], ",", self.SV[1],",", self.SV[2],"] mm")
        print("- Tracks in primary vertex:", len(self.tracksPV),";")
        for t in self.tracksPV:
            t.print()
        print("- Tracks in secondary vertex:", len(self.tracksSV),";")
        for t in self.tracksSV:
            t.print()
        print("- Pileup Tracks:", len(self.tracksSV),";")
        for t in self.tracksSV:
            t.print()
        print("------------------------------------------")

    def draw(self, couples = [], fittedSV = [], fittedTracks = [], drawJetCone = False):
        # If jet has been fed to vertex finder
        if couples != []:
            fig = plt.figure()
            ax = fig.add_subplot(1,3,1,projection='3d')
            # Draw PV and SV
            ax.scatter(self.PV[0],self.PV[1],self.PV[2],color='blue',label='PV')
            ax.scatter(self.SV[0],self.SV[1],self.SV[2],color='red',label='SV')
            # Draw PV's tracks
            label=False
            for track in self.tracksPV:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='orange',label="PV tracks")
                    label = True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='orange')
            
            # Draw SV's tracks
            label = False
            for track in self.tracksSV:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='green',label="SV tracks")
                    label=True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='green')
            
            # Draw pileup tracks
            label = False
            for track in self.tracksPileup:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='black',label="pile-up")
                    label=True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='black')
            
            # Draw B-meson flight
            ax.plot(np.row_stack((self.PV,self.SV))[:,0],np.row_stack((self.PV,self.SV))[:,1],np.row_stack((self.PV,self.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
            
            if drawJetCone:
                u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
                x = 10*np.cos(u)*np.sin(v)
                y = 10*np.sin(u)*np.sin(v)
                z = np.sqrt(x**2 + y**2)/m.tan(self.thetaMaxJet)
                ax.plot_surface(x, y, z,alpha = 0.1)
            
            plt.legend(loc='best')
            plt.title("Jet Event Display")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            ax.set_zlabel("z [mm]")

            # plotting coupled tracks
            ax = fig.add_subplot(1,3,2,projection='3d')
            # Draw PV and SV
            ax.scatter(self.PV[0],self.PV[1],self.PV[2],color='blue',label='PV')
            ax.scatter(self.SV[0],self.SV[1],self.SV[2],color='red',label='SV')
            
            for i,c in enumerate(couples):
                label = False
                for track in c:
                    points = []
                    for t in np.linspace(0,20,10):
                        points.append(track.evaluate(t))
                    points = np.asarray(points)
                    if not label:
                        labelStr = "couple: " + str(i)
                        ax.plot(points[:,0],points[:,1],points[:,2],color=self.colors[i],label=labelStr)
                        label = True
                    else:
                        ax.plot(points[:,0],points[:,1],points[:,2],color=self.colors[i])
            
            # Draw B-meson flight
            ax.plot(np.row_stack((self.PV,self.SV))[:,0],np.row_stack((self.PV,self.SV))[:,1],np.row_stack((self.PV,self.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
            
            if drawJetCone:
                u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
                x = 10*np.cos(u)*np.sin(v)
                y = 10*np.sin(u)*np.sin(v)
                z = np.sqrt(x**2 + y**2)/m.tan(self.thetaMaxJet)
                ax.plot_surface(x, y, z,alpha = 0.1)
            
            plt.legend(loc='best')
            plt.title("After Vertex Finding")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            ax.set_zlabel("z [mm]")
            # If jet has been fed to vertex fitter

            if fittedTracks != []:
                ax = fig.add_subplot(1,3,3,projection='3d')
                # Draw PV and SV
                ax.scatter(self.PV[0],self.PV[1],self.PV[2],color='blue',label='PV')
                ax.scatter(self.SV[0],self.SV[1],self.SV[2],color='red',label='SV')
                # Draw fitted Vertex
                ax.scatter(fittedSV[0],fittedSV[1],fittedSV[2],color='orange',label='Fitted SV')
                # Draw B-meson flight
                ax.plot(np.row_stack((self.PV,self.SV))[:,0],np.row_stack((self.PV,self.SV))[:,1],np.row_stack((self.PV,self.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
                
                # Draw fitted tracks
                label = False
                for track in fittedTracks:
                    points = []
                    for t in np.linspace(0,20,10):
                        points.append(track.evaluate(t))
                    points = np.asarray(points)
                    if not label:
                        ax.plot(points[:,0],points[:,1],points[:,2],color='orange',label="Tracks selected\n$\chi^2<\chi^2_{threshold}$")
                        label = True
                    else:
                        ax.plot(points[:,0],points[:,1],points[:,2],color='orange')
               
                if drawJetCone:
                    u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
                    x = 10*np.cos(u)*np.sin(v)
                    y = 10*np.sin(u)*np.sin(v)
                    z = np.sqrt(x**2 + y**2)/m.tan(self.thetaMaxJet)
                    ax.plot_surface(x, y, z,alpha = 0.1)
                
                plt.legend(loc='best')
                plt.title("After Vertex Fitting")
                ax.set_xlabel("x [mm]")
                ax.set_ylabel("y [mm]")
                ax.set_zlabel("z [mm]")
        # If jet has not been fed in SSVF algorithms
        else:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            # Draw PV and SV
            ax.scatter(self.PV[0],self.PV[1],self.PV[2],color='blue',label='PV')
            ax.scatter(self.SV[0],self.SV[1],self.SV[2],color='red',label='SV')
            # Draw PV's tracks
            label=False
            for track in self.tracksPV:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='orange',label="PV tracks")
                    label = True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='orange')
            
            # Draw SV's tracks
            label = False
            for track in self.tracksSV:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='green',label="SV tracks")
                    label=True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='green')
            
            # Draw pileup tracks
            label = False
            for track in self.tracksPileup:
                points = []
                for t in np.linspace(0,20,10):
                    points.append(track.evaluate(t))
                points = np.asarray(points)
                if not label:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='black',label="pile-up")
                    label=True
                else:
                    ax.plot(points[:,0],points[:,1],points[:,2],color='black')
            
            # Draw B-meson flight
            ax.plot(np.row_stack((self.PV,self.SV))[:,0],np.row_stack((self.PV,self.SV))[:,1],np.row_stack((self.PV,self.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
            
            if drawJetCone:
                u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
                x = 10*np.cos(u)*np.sin(v)
                y = 10*np.sin(u)*np.sin(v)
                z = np.sqrt(x**2 + y**2)/m.tan(self.thetaMaxJet)
                ax.plot_surface(x, y, z,alpha = 0.1)
            
            plt.legend(loc='best')
            plt.title("Jet Event Display")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            ax.set_zlabel("z [mm]")
