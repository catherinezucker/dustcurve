#the normalization component of this code is taken from the plot_fit function, found in the fitsurfs.py program written by Eddie Schlafly
#This can be originally found at ~schlafly/dustsl/python/fitsurfs.py:plot_fit

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def plot_posterior(post_array,xpts,ypts,best,ratio,unique_co,x_range=[4,19],y_range=[0,7],normcol=False, vmin=0.0, vmax=0.03):
    #normalize stellar posterior surfaces
    dispsurf = post_array.copy()
    dispsurf[~np.isfinite(dispsurf)] = 0
    norm = np.sum(np.sum(dispsurf, axis=1), axis=1)+1.e-40
    totsurf = np.sum(dispsurf/norm.reshape(-1,1,1), axis=0)

    if normcol==True:
    #normalize so every distance column has same amount of ink
        colmin = np.min(totsurf, axis=0)
        dispsurf = totsurf-colmin.reshape(1,-1)
        dispsurf = dispsurf/np.sum(dispsurf, axis=0).reshape(1,-1)
    
    else:
        dispsurf=totsurf.copy()

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
    ax.imshow(dispsurf,origin='lower',aspect='auto', cmap='binary', interpolation='nearest',extent=extent, vmin=vmin,vmax=vmax)
    ax.autoscale(False)

    #determine the reddening profile
    distances=best[0:int(len(best)/2)]
    coefficients=best[int(len(best)/2):]
    order=np.argsort(distances)
    
    distances=distances[order]
    distances=np.around(distances,2)
    coefficients=coefficients[order]
    
    for i in range(0,len(unique_co)):
    #convert the actual distance and reddening into something plottable
        conversion=np.multiply(coefficients,ratio)
        reddening=np.cumsum(np.multiply(unique_co[i][order],conversion))
        repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
        profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
        profilex=np.linspace(4,19,1500)

    #plot
        ax.plot(profilex,profiley,alpha=0.1)
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)


def plot_samples(post_array,xpts,ypts,best,samples,ratio,unique_co,x_range=[4,19],y_range=[0,7],normcol=False, vmin=0.0, vmax=0.03):
    
    #normalize stellar posterior surfaces so every individual array sums to one
    dispsurf = post_array.copy()
    dispsurf[~np.isfinite(dispsurf)] = 0
    norm = np.sum(np.sum(dispsurf, axis=1), axis=1)+1.e-40
    totsurf = np.sum(dispsurf/norm.reshape(-1,1,1), axis=0)

    if normcol==True:
    #normalize so every distance column has same amount of ink
        colmin = np.min(totsurf, axis=0)
        dispsurf = totsurf-colmin.reshape(1,-1)
        dispsurf = dispsurf/np.sum(dispsurf, axis=0).reshape(1,-1)
    
    else:
        dispsurf=totsurf.copy()

    #setup and plot normalized stacked stellar posterior
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(1,1,1)
    ax.set_title("Stacked Stellar Posterior w/ Reddening Profile")
    ax.set_xlabel(r'$\mu$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{E} \left( B - V \right)$', fontsize=16)
    if xpts is not None and ypts is not None:
        dx = np.median(xpts[1:]-xpts[:-1])
        dy = np.median(ypts[1:]-ypts[:-1])
        extent = [xpts[0]-dx/2., xpts[-1]+dx/2.,ypts[0]-dy/2., ypts[-1]+dy/2.]
    ax.imshow(dispsurf,origin='lower',aspect='auto', cmap='binary', interpolation='nearest',extent=extent, vmin=vmin,vmax=vmax)
    ax.autoscale(False)
    
    #take the average CO of all pixels at every slice
    avg_co=np.mean(unique_co,axis=0)
    
    nsamples=1000
    #draw random samples from traces
    for i in range(0,nsamples):
        index=np.random.randint(0,traces_cold.shape[1])
        step_samp=samples[:,index]
        #determine the reddening profile
        distances=best[0:int(len(step_samp)/2)]
        coefficients=best[int(len(step_samp)/2):]
        order=np.argsort(distances)
    
        distances=distances[order]
        distances=np.around(distances,2)
        coefficients=coefficients[order]

        #convert the actual distance and reddening samples from given step into something plottable
        conversion=np.multiply(coefficients,ratio)
        reddening=np.cumsum(np.multiply(avg_co[order],conversion))
        repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
        profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
        profilex=np.linspace(4,19,1500)

        #plot
        sample_prof=ax.plot(profilex,profiley,c='g',label="Samples", lw=3)
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        
    #now overlay best fit reddening profile
    distances=best[0:int(len(best)/2)]
    coefficients=best[int(len(best)/2):]
    order=np.argsort(distances)
    
    distances=distances[order]
    distances=np.around(distances,2)
    coefficients=coefficients[order]
    
    #convert the actual distance and reddening into something plottable
    conversion=np.multiply(coefficients,ratio)
    reddening=np.cumsum(np.multiply(avg_co[order],conversion))
    repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
    profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
    profilex=np.linspace(4,19,1500)

    #plot
    best_prof=ax.plot(profilex,profiley, c='b',label="Best fit")    
    
    #create legend for samples and best fit
    light_blue_line = mlines.Line2D([], [], color='green',label='Samples',)
    dark_blue_line = mlines.Line2D([], [], color='blue', label='Best fit')
    
    plt.legend(handles=[light_blue_line, dark_blue_line])
