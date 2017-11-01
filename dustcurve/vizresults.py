import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from dustcurve import globalvars as gv
import pandas as pd
from pandas.tools.plotting import table
from astropy.io import fits
import dynesty
from dynesty import plotting as dyplot



xpts=np.linspace(4,19,240)
ypts=np.linspace(0,7,700)
x_range=[4,19]
y_range=[0,7]


def plot_samples(best,traces_cold,dname,cubefile,normcol=True,co_ratio=0.004,num_samples=250,CO_samp=False,weights=None):

    nslices=int((len(best)-3)/2)

    #normalize stellar posterior surfaces so every individual array sums to one                                                                 
    dispsurf=np.vstack((gv.unique_post[dname]))

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
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r'$\mu$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{E} \left( B - V \right)$', fontsize=16)
    if xpts is not None and ypts is not None:
        dx = np.median(xpts[1:]-xpts[:-1])
        dy = np.median(ypts[1:]-ypts[:-1])
        extent = [xpts[0]-dx/2., xpts[-1]+dx/2.,ypts[0]-dy/2., ypts[-1]+dy/2.]
    ax.imshow(np.sqrt(dispsurf),origin='lower',aspect='auto', cmap='binary', interpolation='nearest',extent=extent)       
    ax.autoscale(False)

    #Now overlay locations of peak probability for each individual stellar posterior on distance+reddening
    dispsurf=np.vstack((gv.unique_post[dname]))
    dispsurf_list=[dispsurf[i,:,:] for i in range(dispsurf.shape[0])]
    maxlike=[np.where(dispsurf_list[i] == dispsurf_list[i].max()) for i in range(len(dispsurf_list))]
    maxred=[maxlike[i][0][0] for i in range(len(maxlike))]
    maxdist=[maxlike[i][1][0] for i in range(len(maxlike))]
    maxredplot=np.array(maxred)*(7.0/700.0)
    maxdistplot=4+np.array(maxdist)*(15.0/240.0)
    
    #take the average CO of all pixels at every slice                                                                                   
    avg_co=np.nanmean(gv.unique_co[dname],axis=0)

    if weights is not None:
        weight_draw=weights
    else:
        weight_draw=np.ones((traces_cold.shape[1]))

    if CO_samp==True:
        pull_samps=np.random.randint(0,len(gv.unique_co[dname]),size=num_samples)
    else:
        pull_samps=np.random.choice(np.arange(0,traces_cold.shape[1]),size=num_samples,p=weight_draw)

    for i in pull_samps:

        #Plot median fit for different CO pixels
        if CO_samp==True:
            fore_dist=best[0]
            fore_red=best[1]
            distances=best[3:nslices+3]
            diff_reddening=best[nslices+3:]*co_ratio*gv.unique_co[dname][i,:]
            
            #sort by distance and then determine cumulative reddening along l.o.s.
            order=np.argsort(distances)
            distances=np.hstack((fore_dist,distances[order]))
            reddening=np.cumsum(np.hstack((fore_red,diff_reddening[order])))

            #define plottable profiles
            repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
            profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
            profilex=np.linspace(4,19,len(profiley))

        #Plot samples given average CO pixel
        else:

            fore_dist=traces_cold[:,i][0]
            fore_red=traces_cold[:,i][1]
            distances=traces_cold[:,i][3:nslices+3]
            diff_reddening=traces_cold[:,i][nslices+3:]*co_ratio*avg_co
            
            #sort by distance and then determine cumulative reddening along l.o.s.
            order=np.argsort(distances)
            distances=np.hstack((fore_dist,distances[order]))
            reddening=np.cumsum(np.hstack((fore_red,diff_reddening[order])))

            #define plottable profiles
            repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
            profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
            profilex=np.linspace(4,19,len(profiley))
            
        #plot                                                                                                                                                       
        sample_prof=ax.plot(profilex,profiley,c='blue',label="Filler", lw=1,alpha=0.1,zorder=1)
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)

    fore_dist=best[0]
    fore_red=best[1]
    distances=best[3:nslices+3]
    diff_reddening=best[nslices+3:]*co_ratio*avg_co
            
    #sort by distance and then determine cumulative reddening along l.o.s.
    order=np.argsort(distances)
    distances=np.hstack((fore_dist,distances[order]))
    reddening=np.cumsum(np.hstack((fore_red,diff_reddening[order])))

    #define plottable profiles
    repeats=np.around((np.ediff1d(distances,to_end=19.0-distances[-1],to_begin=distances[0]-4)/0.01),0)
    profiley=np.repeat(np.insert(reddening,0,0),repeats.astype(int))
    profilex=np.linspace(4,19,len(profiley))
    best_prof=ax.plot(profilex,profiley,c='red',lw=5,zorder=2)

    #As last step, overplot peaks of underlying stellar posteriors
    ax.scatter(maxdistplot,maxredplot,marker='o',c='green',s=5,alpha=0.5,zorder=3,edgecolor="None")

    #create legend for samples and best fit                                                                                                                      
    red_line = mlines.Line2D([], [], color='red',label='Median profile with average CO values',lw=5)

    if CO_samp==True:
        blue_line = mlines.Line2D([], [], color='blue', label='Median fit w/ sampled CO values')
    else:
        blue_line = mlines.Line2D([], [], color='blue', label='Samples  w/ average CO values')

    
    ax.legend(handles=[blue_line,red_line])
    
    ax.set_ylim(0,2)

    ax.set_xticks(np.arange(4,19,1))
    
    #setup labels
    labels=['dfore','rfore','P_b']
    
    for n in range(0,nslices):
        labels.append('d{:d}'.format(n+1))
    
    for n in range(0,nslices):
        labels.append('c{:d}'.format(n+1))
    

    #Format and plot tabel of median samples with errorbounds
    if weights is None:
        theta=pd.DataFrame(traces_cold)
        quantile16,quantile50,quantile84,=theta.quantile([0.16,.5,0.84], axis=1).values

    else:
        quantile16=[]
        quantile50=[]
        quantile84=[]

        for i in range(0,traces_cold.T.shape[1]):
            quantile16.append(dyplot._quantile(traces_cold.T[:,i], [0.16], weights=weights))
            quantile50.append(dyplot._quantile(traces_cold.T[:,i], [0.50], weights=weights))
            quantile84.append(dyplot._quantile(traces_cold.T[:,i], [0.84], weights=weights))
    
        quantile16=np.ravel(np.array(quantile16))
        quantile50=np.ravel(np.array(quantile50))
        quantile84=np.ravel(np.array(quantile84))

    lowerbound=np.round(quantile50-quantile16,2)
    upperbound=np.round(quantile84-quantile50,2)
    quantile50=np.round(quantile50,2)

    cell_text_arr=[]
    for (low,med,high) in zip(lowerbound,quantile50,upperbound):
        cell_text_arr.append(r"${:.2f}_{{-{:.2f}}}^{{+{:.2f}}}$".format(med,low,high))

    cell_text_arr=np.array(cell_text_arr)

    df_tableplot=pd.DataFrame(cell_text_arr.reshape(1,cell_text_arr.shape[0]),columns=labels)

    hdr=fits.getheader("/n/fink2/czucker/CO_Input/"+cubefile)

    if nslices>1:
        spectral_axis=np.round(np.arange(hdr['CRVAL3'],hdr['CRVAL3']+(hdr['NAXIS3']*hdr['CDELT3']),hdr['CDELT3']),1)
        spectral_axis=spectral_axis[0:hdr['NAXIS3']]
    else:
        spectral_axis=np.round(hdr['CRVAL3'],1)

    velarr=np.ones((3+nslices*2))*np.nan
    velarr[3:3+nslices]=spectral_axis
    velarr=velarr.astype(str)
    velarr[velarr=='nan']=''
    velarr=velarr.tolist()
    for i in range(3,nslices+3):
        velarr[i]=r"{} km/s".format(velarr[i])

    df_tableplot.loc[-1]=velarr

    tab=table(ax,df_tableplot,
          rowLabels=['',''],
          colLabels=labels,
          cellLoc = 'center', rowLoc = 'center',
          loc='top',bbox=[0.0, 1.05, 1.0, 0.2])

    tab.auto_set_column_width(labels)
    
    return ax
