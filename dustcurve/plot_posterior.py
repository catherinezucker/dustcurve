#the normalization component of this code is taken from the plot_fit function, found in the fitsurfs.py program written by Eddie Schlafly
#This can be originally found at ~schlafly/dustsl/python/fitsurfs.py:plot_fit

import numpy as np
import matplotlib.pyplot as plt

def plot_posterior(post_array,xpts,ypts,best,ratio,unique_co):
    #normalize stellar posterior surfaces
    dispsurf = post_array.copy()
    dispsurf[~np.isfinite(dispsurf)] = 0
    norm = np.sum(np.sum(dispsurf, axis=1), axis=1)+1.e-40
    totsurf = np.sum(dispsurf/norm.reshape(-1,1,1), axis=0)

    #setup and plot normalized stacked stellar posterior
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title("Stacked Stellar Posterior w/ Reddening Profile")
    ax.set_xlabel(r'$\mu$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{E} \left( B - V \right)$', fontsize=16)
    if xpts is not None and ypts is not None:
        dx = np.median(xpts[1:]-xpts[:-1])
        dy = np.median(ypts[1:]-ypts[:-1])
        extent = [xpts[0]-dx/2., xpts[-1]+dx/2.,ypts[0]-dy/2., ypts[-1]+dy/2.]
    ax.imshow(totsurf,origin='lower',aspect='auto', cmap='binary', interpolation='nearest',extent=extent, vmin=0,vmax=0.03)
    ax.autoscale(False)

    #determine the reddening profile
    distances=best[0:int(len(best)/2)]
    reddenings=np.cumsum(np.multiply(best[int(len(best)/2):],ratio))
    
    order=np.argsort(distances)
    unique_co=unique_co[0][order]
    distances=distances[order]
    distances=np.around(distances,2)
    reddenings=np.cumsum(np.multiply(unique_co,ratio))
    
    #convert the actual distance and reddening into something plottable
    repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
    profiley=np.repeat(np.insert(reddenings,0,0),repeats.astype(int))
    profilex=np.linspace(4,19,1500)

    #plot
    ax.plot(profilex,profiley)
    
