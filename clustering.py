from modules.JetGenerator import JetGenerator
from modules.SSVF import SSVF
from matplotlib import pyplot as plt
import numpy as np
from sklearn.cluster import KMeans as KM
import math as m

if __name__=="__main__":
    
    generator = JetGenerator(1)
    #generator = JetGenerator(1,nTracksPV=10,nTracksSV=15,nTracksPileup=10)
    jet = generator.generate()
    jet = jet[0]
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
    k = 3
    km = KM(k,max_iter=1000)
    km.fit(jetForClustering)
    couples = km.labels_
    centroids = km.cluster_centers_
    print(couples)

    fig,ax = plt.subplots(1,2)
    plt.suptitle("Clustering in the feature space")

    ax[1].set_title("Clustered features")
    ax[1].scatter(jetForClustering[:,0],jetForClustering[:,1],c=couples,label='tracks')
    ax[1].scatter(centroids[:,0],centroids[:,1],marker='*',color='red',label='centroids')
    ax[1].legend(loc='best')
    ax[1].set_xlabel("IP [mm]")
    ax[1].set_ylabel("$\cos(\eta)$")

    ax[0].set_title("MC truth")
    jetForClustering = []
    for t in jet.tracksPV:
        track = [np.abs(t.IP),t.versor[2]]
        jetForClustering.append(track)
    jetForClustering = np.array(jetForClustering)
    ax[0].scatter(jetForClustering[:,0],jetForClustering[:,1],label="PV")
    jetForClustering = []
    for t in jet.tracksSV:
        track = [np.abs(t.IP),t.versor[2]]
        jetForClustering.append(track)
    jetForClustering = np.array(jetForClustering)
    ax[0].scatter(jetForClustering[:,0],jetForClustering[:,1],label="SV")
    jetForClustering = []
    for t in jet.tracksPileup:
        track = [np.abs(t.IP),t.versor[2]]
        jetForClustering.append(track)
    jetForClustering = np.array(jetForClustering)
    ax[0].scatter(jetForClustering[:,0],jetForClustering[:,1],label="Pileup")
    
    ax[0].set_xlabel("IP [mm]")
    ax[0].set_ylabel("$\cos(\eta)$")
    ax[0].legend(loc='best')

    fig = plt.figure()
    ax = fig.add_subplot(1,2,1,projection='3d')
    # Draw PV and SV
    ax.scatter(jet.PV[0],jet.PV[1],jet.PV[2],color='blue',label='PV')
    ax.scatter(jet.SV[0],jet.SV[1],jet.SV[2],color='red',label='SV')
    # Draw PV's tracks
    label=False
    for track in jet.tracksPV:
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
    for track in jet.tracksSV:
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
    for track in jet.tracksPileup:
        points = []
        for t in np.linspace(-10,20,10):
            points.append(track.evaluate(t))
        points = np.asarray(points)
        if not label:
            ax.plot(points[:,0],points[:,1],points[:,2],color='black',label="pile-up")
            label=True
        else:
            ax.plot(points[:,0],points[:,1],points[:,2],color='black')
    
    # Draw B-meson flight
    ax.plot(np.row_stack((jet.PV,jet.SV))[:,0],np.row_stack((jet.PV,jet.SV))[:,1],np.row_stack((jet.PV,jet.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
    
    # If enabled, draw jet's cone
    if True:
        u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
        x = 10*np.cos(u)*np.sin(v)
        y = 10*np.sin(u)*np.sin(v)
        z = np.sqrt(x**2 + y**2)/m.tan(jet.thetaMaxJet)
        ax.plot_surface(x, y, z,alpha = 0.1)
    # legend and axis labels
    plt.legend(loc='best')
    plt.title("Jet Event Display")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    
    ax = fig.add_subplot(1,2,2,projection='3d')
    # Draw PV and SV
    ax.scatter(jet.PV[0],jet.PV[1],jet.PV[2],color='blue',label='PV')
    ax.scatter(jet.SV[0],jet.SV[1],jet.SV[2],color='red',label='SV')
    colors = ["blue","orange","green","red","purple","brown","pink","gray","olive","cyan"]
    # Draw PV's tracks
    label=True
    for i,track in enumerate(jet.tracksPV):
        points = []
        for t in np.linspace(0,20,10):
            points.append(track.evaluate(t))
        points = np.asarray(points)
        if not label:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i]],label="PV tracks")
            label = True
        else:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i]])
    
    # Draw SV's tracks
    label = True
    delta = len(jet.tracksPV)
    for i,track in enumerate(jet.tracksSV):
        points = []
        for t in np.linspace(0,20,10):
            points.append(track.evaluate(t))
        points = np.asarray(points)
        if not label:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i+delta]],label="SV tracks")
            label=True
        else:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i+delta]])
    
    # Draw pileup tracks
    label = True
    delta += len(jet.tracksSV)
    for i,track in enumerate(jet.tracksPileup):
        points = []
        for t in np.linspace(-10,20,10):
            points.append(track.evaluate(t))
        points = np.asarray(points)
        if not label:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i+delta]],label="pile-up")
            label=True
        else:
            ax.plot(points[:,0],points[:,1],points[:,2],color=colors[couples[i+delta]])
    
    # Draw B-meson flight
    ax.plot(np.row_stack((jet.PV,jet.SV))[:,0],np.row_stack((jet.PV,jet.SV))[:,1],np.row_stack((jet.PV,jet.SV))[:,2],color="gray",linestyle="dashed",label="B meson")
    
    # If enabled, draw jet's cone
    if True:
        u, v = np.mgrid[0:2*m.pi:100j, 0:m.pi:80j]
        x = 10*np.cos(u)*np.sin(v)
        y = 10*np.sin(u)*np.sin(v)
        z = np.sqrt(x**2 + y**2)/m.tan(jet.thetaMaxJet)
        ax.plot_surface(x, y, z,alpha = 0.1)
    # legend and axis labels
    plt.legend(loc='best')
    plt.title("K-Mean clustered Event\n$k=$" + str(k))
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    plt.show()